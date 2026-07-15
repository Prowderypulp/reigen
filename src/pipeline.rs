//! Conversion pipeline: reader → (filters) → writer.
//!
//! # Dispatch
//!
//! Each format has a module implementing `GenoReader` and/or `GenoWriter`.
//! This file opens the correct pair based on input/output `Format`, boxes
//! them as trait objects, and runs the common streaming loop.
//!
//! # Phase 2 scope
//!
//! All SnpMajor↔SnpMajor pairs: PAM, EIGENSTRAT, (PACKEDPED after next slice).
//! SampleMajor formats (TGENO) and ANCESTRYMAP sparse are still stubbed.
//!
//! # Streaming loop
//!
//! ```text
//! read metadata  → Vec<SnpRow>, Vec<IndRow>
//! build filters  → keep_snps, keep_inds bitmasks
//! open reader    → nind/nsnp verified against metadata
//! open writer    → begin(kept_nind, kept_nsnp, 0, 0)
//! for each input record:
//!     if !keep_snps[i]: skip
//!     if all samples kept: pass record through as-is
//!     else: unpack → project → repack
//!     write
//! finish writer, write output metadata
//! ```

use crate::filter::{
    load_bad_snps, load_pop_keep, load_sample_keep, load_sample_remove, load_snp_keep, ChromFilter,
    IndFilter, SnpFilter,
};
use crate::format::{self, Format};
use crate::geno::{codec, GenoReader, GenoWriter, Layout};
use crate::geno::{eigenstrat, packed_am, packed_ped, pgen, tgeno};
use crate::meta::{self, IndRow, SnpRow};
use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Total installed physical memory in bytes, read from `/proc/meminfo`.
/// Returns `None` when it can't be determined (non-Linux, no `/proc`, or an
/// unparseable file), leaving callers to fall back to a fixed budget.
pub fn detect_total_memory() -> Option<u64> {
    let meminfo = std::fs::read_to_string("/proc/meminfo").ok()?;
    parse_meminfo_total(&meminfo)
}

/// Parse the `MemTotal:` line (reported in kB) from `/proc/meminfo` content.
fn parse_meminfo_total(meminfo: &str) -> Option<u64> {
    for line in meminfo.lines() {
        if let Some(rest) = line.strip_prefix("MemTotal:") {
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb.saturating_mul(1024));
        }
    }
    None
}

/// Default cross-layout transpose memory budget: **half of installed RAM**,
/// so a conversion never tries to claim the whole machine. Falls back to 4 GiB
/// when RAM can't be detected, and is floored at 512 MiB so it is always
/// workable. Overridden by the `--max-mem` flag; also used by internal callers
/// (filter, merge) that don't expose the flag.
pub fn default_max_mem() -> u64 {
    const FALLBACK: u64 = 4 * 1024 * 1024 * 1024;
    const FLOOR: u64 = 512 * 1024 * 1024;
    match detect_total_memory() {
        Some(total) => (total / 2).max(FLOOR),
        None => FALLBACK,
    }
}

#[derive(Debug, Clone)]
pub struct ConvertConfig {
    pub geno_in: PathBuf,
    pub snp_in: PathBuf,
    pub ind_in: PathBuf,
    pub out_fmt: Format,
    pub geno_out: PathBuf,
    pub snp_out: PathBuf,
    pub ind_out: PathBuf,
    pub badsnp: Option<PathBuf>,
    pub snps: Option<PathBuf>,
    pub poplist: Option<PathBuf>,
    pub keep: Option<PathBuf>,
    pub remove: Option<PathBuf>,
    pub chrom: Option<String>,
    pub lopos: Option<u64>,
    pub hipos: Option<u64>,
    pub noxdata: bool,
    pub max_miss_snp: Option<f64>,
    pub max_miss_ind: Option<f64>,
    pub maf: Option<f64>,
    pub max_maf: Option<f64>,
    pub hwe: Option<f64>,
    pub numchrom: u32,
    pub hashcheck: bool,
    pub familynames: bool,
    pub outputgroup: bool,
    /// Peak-memory budget (bytes) for the cross-layout transpose. When the
    /// packed matrix fits, the fast single-pass in-memory transpose is used;
    /// otherwise the transpose is banded to stay within this budget (at the
    /// cost of re-streaming the source once per band). See `stream_cross_layout`.
    pub max_mem: u64,
}

pub fn run_convert(cfg: &ConvertConfig) -> Result<()> {
    let in_fmt = format::infer_input_format(&cfg.geno_in).context("inferring input format")?;

    log::info!("input  format: {in_fmt:?}");
    log::info!("output format: {:?}", cfg.out_fmt);

    let t0 = Instant::now();
    let numchrom = cfg.numchrom;

    // --- 1. Metadata.
    let snp_rows = read_input_snp(&cfg.snp_in, in_fmt, numchrom)?;
    let ind_rows = read_input_ind(&cfg.ind_in, in_fmt, cfg.familynames)?;

    // PLINK 2: write the multiallelic drop report next to the output geno
    // (mirrors merge's `.missnp`; plan §4.4 / acceptance #5).
    if in_fmt == Format::Plink2 {
        let dropped = meta::pvar::dropped_multiallelic(&cfg.snp_in)?;
        if !dropped.is_empty() {
            let report = cfg.geno_out.with_extension("pgen-drop.tsv");
            meta::pvar::write_drop_report(&report, &dropped)?;
            log::info!(
                "wrote {} dropped multiallelic variant(s) to {}",
                dropped.len(),
                report.display()
            );
        }
    }
    log::info!(
        "metadata: {} SNPs, {} samples (read in {:.2?})",
        snp_rows.len(),
        ind_rows.len(),
        t0.elapsed()
    );

    // --- 2. Filters.
    let bad_snps = cfg.badsnp.as_deref().map(load_bad_snps).transpose()?;
    let snp_keep = cfg.snps.as_deref().map(load_snp_keep).transpose()?;
    let pop_keep = cfg.poplist.as_deref().map(load_pop_keep).transpose()?;
    let sample_keep = cfg.keep.as_deref().map(load_sample_keep).transpose()?;
    let sample_remove = cfg.remove.as_deref().map(load_sample_remove).transpose()?;
    let chrom_filter = cfg.chrom.as_deref().map(ChromFilter::parse).transpose()?;
    validate_missingness_threshold("geno", cfg.max_miss_snp)?;
    validate_missingness_threshold("mind", cfg.max_miss_ind)?;
    validate_maf_threshold("maf", cfg.maf)?;
    validate_maf_threshold("max-maf", cfg.max_maf)?;
    validate_hwe_threshold("hwe", cfg.hwe)?;
    let x_chrom = u8::try_from(numchrom + 1).context("numchrom too large")?;
    let y_chrom = u8::try_from(numchrom + 2).context("numchrom too large")?;
    let mt_chrom = u8::try_from(numchrom + 3).context("numchrom too large")?;
    let xy_chrom = u8::try_from(numchrom + 4).context("numchrom too large")?;

    let snp_filter = SnpFilter {
        bad: bad_snps.as_ref(),
        snp_keep: snp_keep.as_ref(),
        chrom: chrom_filter,
        lopos: cfg.lopos,
        hipos: cfg.hipos,
        noxdata: cfg.noxdata,
        x_chrom,
        y_chrom,
        mt_chrom,
        xy_chrom,
    };
    let ind_filter = IndFilter {
        pop_keep: pop_keep.as_ref(),
        sample_keep: sample_keep.as_ref(),
        sample_remove: sample_remove.as_ref(),
    };

    let mut keep_snps: Vec<bool> = snp_rows.iter().map(|s| snp_filter.keep(s)).collect();
    let mut keep_inds: Vec<bool> = ind_rows.iter().map(|i| ind_filter.keep(i)).collect();
    let mut kept_snp_count = keep_snps.iter().filter(|&&k| k).count();
    let mut kept_ind_count = keep_inds.iter().filter(|&&k| k).count();

    if cfg.max_miss_snp.is_some()
        || cfg.max_miss_ind.is_some()
        || cfg.maf.is_some()
        || cfg.max_maf.is_some()
        || cfg.hwe.is_some()
    {
        apply_stat_filters(
            cfg,
            in_fmt,
            &mut keep_snps,
            &mut keep_inds,
            ind_rows.len(),
            snp_rows.len(),
        )?;
        kept_snp_count = keep_snps.iter().filter(|&&k| k).count();
        kept_ind_count = keep_inds.iter().filter(|&&k| k).count();
    }
    log::info!(
        "after filters: {} SNPs, {} samples",
        kept_snp_count,
        kept_ind_count
    );
    if kept_snp_count == 0 {
        bail!(
            "all {} SNPs were removed by the filters — nothing to write\n  \
             hint: relax the SNP filters (--maf/--max-maf/--hwe/--max-miss-snp/--chrom/\
             --from-bp/--to-bp/--snps/--badsnp) or check they match this dataset",
            snp_rows.len()
        );
    }
    if kept_ind_count == 0 {
        bail!(
            "all {} samples were removed by the filters — nothing to write\n  \
             hint: relax the sample filters (--mind/--keep/--remove/--poplist) or check \
             the IDs in your keep/remove list match the .ind/.fam",
            ind_rows.len()
        );
    }

    // --- 3. Open reader + writer (boxed trait objects).
    let mut reader = open_reader(in_fmt, &cfg.geno_in, ind_rows.len(), snp_rows.len())?;

    if cfg.hashcheck {
        if let Some((file_ihash, file_shash)) = reader.header_hashes() {
            let in_ind_ids: Vec<&str> = ind_rows.iter().map(|i| i.id.as_str()).collect();
            let in_snp_ids: Vec<&str> = snp_rows.iter().map(|s| s.id.as_str()).collect();
            let exp_ihash = crate::hash::hasharr(&in_ind_ids);
            let exp_shash = crate::hash::hasharr(&in_snp_ids);
            if file_ihash != exp_ihash || file_shash != exp_shash {
                bail!(
                    "hashcheck FAILED for {}:\n  \
                     header ihash={:08x} shash={:08x}\n  \
                     computed ihash={:08x} shash={:08x}\n  \
                     The .ind / .snp files do not match the .geno / .tgeno that was written. \
                     Either regenerate the geno file with the current metadata, or set \
                     `hashcheck: NO` to bypass.",
                    cfg.geno_in.display(),
                    file_ihash,
                    file_shash,
                    exp_ihash,
                    exp_shash,
                );
            }
            log::info!(
                "hashcheck OK (ihash={:08x} shash={:08x})",
                file_ihash,
                file_shash
            );
        }
    }

    // Compute hashes for the OUTPUT geno header from the kept IDs.
    let out_ind_ids: Vec<&str> = ind_rows
        .iter()
        .zip(keep_inds.iter())
        .filter_map(|(i, &k)| if k { Some(i.id.as_str()) } else { None })
        .collect();
    let out_snp_ids: Vec<&str> = snp_rows
        .iter()
        .zip(keep_snps.iter())
        .filter_map(|(s, &k)| if k { Some(s.id.as_str()) } else { None })
        .collect();
    let out_ihash = crate::hash::hasharr(&out_ind_ids);
    let out_shash = crate::hash::hasharr(&out_snp_ids);

    let mut writer = open_writer(cfg.out_fmt, &cfg.geno_out)?;
    writer.begin(kept_ind_count, kept_snp_count, out_ihash, out_shash)?;

    // --- 4. Layout compatibility.
    if reader.layout() == writer.layout() {
        stream_same_layout(
            reader.as_mut(),
            writer.as_mut(),
            &keep_snps,
            &keep_inds,
            ind_rows.len(),
            kept_ind_count,
        )?;
    } else {
        stream_cross_layout(
            reader.as_mut(),
            writer.as_mut(),
            in_fmt,
            &cfg.geno_in,
            &keep_snps,
            &keep_inds,
            ind_rows.len(),
            snp_rows.len(),
            kept_ind_count,
            kept_snp_count,
            cfg.max_mem,
        )?;
    }
    writer.finish()?;

    // --- 6. Output metadata.
    let kept_snps: Vec<SnpRow> = snp_rows
        .into_iter()
        .zip(keep_snps.iter())
        .filter_map(|(s, &k)| if k { Some(s) } else { None })
        .collect();
    let kept_inds: Vec<IndRow> = ind_rows
        .into_iter()
        .zip(keep_inds.iter())
        .filter_map(|(i, &k)| if k { Some(i) } else { None })
        .collect();
    write_output_snp(&cfg.snp_out, cfg.out_fmt, &kept_snps, numchrom)?;
    write_output_ind(&cfg.ind_out, cfg.out_fmt, &kept_inds, cfg.outputgroup)?;

    log::info!("done in {:.2?}", t0.elapsed());
    Ok(())
}

fn validate_missingness_threshold(name: &str, v: Option<f64>) -> Result<()> {
    if let Some(x) = v {
        if !(0.0..=1.0).contains(&x) {
            bail!("--{name} must be in [0,1], got {x}");
        }
    }
    Ok(())
}

fn validate_maf_threshold(name: &str, v: Option<f64>) -> Result<()> {
    if let Some(x) = v {
        if !(0.0..=0.5).contains(&x) {
            bail!("--{name} must be in [0,0.5], got {x}");
        }
    }
    Ok(())
}

fn validate_hwe_threshold(name: &str, v: Option<f64>) -> Result<()> {
    if let Some(x) = v {
        if x <= 0.0 || x > 1.0 {
            bail!("--{name} must be in (0,1], got {x}");
        }
    }
    Ok(())
}

fn apply_stat_filters(
    cfg: &ConvertConfig,
    in_fmt: Format,
    keep_snps: &mut [bool],
    keep_inds: &mut [bool],
    total_inds: usize,
    total_snps: usize,
) -> Result<()> {
    // PLINK applies --mind (per-sample missingness) BEFORE --geno (per-SNP
    // missingness): mind is computed over the full variant set, then geno over
    // the surviving samples. We match that order so combined --geno/--mind
    // results agree with PLINK. Doing geno first would compute mind only over
    // the already-cleaned (low-missingness) SNPs, silently weakening --mind to
    // a near no-op when paired with --geno (B-012).

    // First phase: per-sample missingness (mind) over the current SNP mask
    // (all SNPs passing the structural snp_filter, before any geno drop).
    if let Some(max_miss_ind) = cfg.max_miss_ind {
        let mut reader = open_reader(in_fmt, &cfg.geno_in, total_inds, total_snps)?;
        let (_, ind_missing) = compute_missing_counts(
            reader.as_mut(),
            keep_snps,
            keep_inds,
            total_inds,
            total_snps,
        )?;
        let denom = keep_snps.iter().filter(|&&k| k).count() as f64;
        if denom > 0.0 {
            for (i, keep) in keep_inds.iter_mut().enumerate() {
                if *keep {
                    let miss_frac = (ind_missing[i] as f64) / denom;
                    if miss_frac > max_miss_ind {
                        *keep = false;
                    }
                }
            }
        }
    }

    // Second phase: per-SNP missingness (geno) over the post-mind sample mask.
    if let Some(max_miss_snp) = cfg.max_miss_snp {
        let mut reader = open_reader(in_fmt, &cfg.geno_in, total_inds, total_snps)?;
        let (snp_missing, _) = compute_missing_counts(
            reader.as_mut(),
            keep_snps,
            keep_inds,
            total_inds,
            total_snps,
        )?;
        let denom = keep_inds.iter().filter(|&&k| k).count() as f64;
        if denom > 0.0 {
            for (j, keep) in keep_snps.iter_mut().enumerate() {
                if *keep {
                    let miss_frac = (snp_missing[j] as f64) / denom;
                    if miss_frac > max_miss_snp {
                        *keep = false;
                    }
                }
            }
        }
    }

    // Third phase: apply MAF and HWE filters
    if cfg.maf.is_some() || cfg.max_maf.is_some() || cfg.hwe.is_some() {
        let mut reader = open_reader(in_fmt, &cfg.geno_in, total_inds, total_snps)?;
        let stats = compute_snp_stats(
            reader.as_mut(),
            keep_snps,
            keep_inds,
            total_inds,
            total_snps,
        )?;

        let min_maf = cfg.maf.unwrap_or(0.0);
        let max_maf = cfg.max_maf.unwrap_or(1.0);
        let hwe_thresh = cfg.hwe.unwrap_or(0.0);

        for (j, keep) in keep_snps.iter_mut().enumerate() {
            if *keep {
                let st = &stats[j];
                let maf = st.maf();

                if maf.is_nan() {
                    // All samples missing for this SNP. Drop if any MAF filter is active.
                    if cfg.maf.is_some() || cfg.max_maf.is_some() {
                        *keep = false;
                        continue;
                    }
                } else if maf < min_maf || maf > max_maf {
                    *keep = false;
                    continue;
                }

                let hwe_p = st.hwe_pvalue();
                if !hwe_p.is_nan() && hwe_p < hwe_thresh {
                    *keep = false;
                }
            }
        }
    }

    Ok(())
}

fn compute_missing_counts(
    reader: &mut dyn GenoReader,
    keep_snps: &[bool],
    keep_inds: &[bool],
    total_inds: usize,
    total_snps: usize,
) -> Result<(Vec<u32>, Vec<u32>)> {
    let mut snp_missing = vec![0u32; total_snps];
    let mut ind_missing = vec![0u32; total_inds];
    let mut in_buf = vec![0u8; reader.record_bytes()];

    match reader.layout() {
        Layout::SnpMajor => {
            let mut unpacked = vec![0u8; total_inds];
            let mut snp_idx = 0usize;
            while reader.read_record(&mut in_buf)? {
                if snp_idx >= total_snps {
                    break;
                }
                if keep_snps[snp_idx] {
                    codec::unpack(&in_buf, total_inds, &mut unpacked);
                    for (i, &g) in unpacked.iter().enumerate() {
                        if keep_inds[i] && g == codec::G_MISSING {
                            snp_missing[snp_idx] += 1;
                            ind_missing[i] += 1;
                        }
                    }
                }
                snp_idx += 1;
            }
        }
        Layout::SampleMajor => {
            let mut unpacked = vec![0u8; total_snps];
            let mut ind_idx = 0usize;
            while reader.read_record(&mut in_buf)? {
                if ind_idx >= total_inds {
                    break;
                }
                if keep_inds[ind_idx] {
                    codec::unpack(&in_buf, total_snps, &mut unpacked);
                    for (j, &g) in unpacked.iter().enumerate() {
                        if keep_snps[j] && g == codec::G_MISSING {
                            snp_missing[j] += 1;
                            ind_missing[ind_idx] += 1;
                        }
                    }
                }
                ind_idx += 1;
            }
        }
    }

    Ok((snp_missing, ind_missing))
}

fn compute_snp_stats(
    reader: &mut dyn GenoReader,
    keep_snps: &[bool],
    keep_inds: &[bool],
    total_inds: usize,
    total_snps: usize,
) -> Result<Vec<crate::stats::per_snp::SnpStats>> {
    let mut stats = vec![crate::stats::per_snp::SnpStats::default(); total_snps];
    let mut in_buf = vec![0u8; reader.record_bytes()];

    match reader.layout() {
        Layout::SnpMajor => {
            let mut unpacked = vec![0u8; total_inds];
            let mut snp_idx = 0usize;
            while reader.read_record(&mut in_buf)? {
                if snp_idx >= total_snps {
                    break;
                }
                if keep_snps[snp_idx] {
                    codec::unpack(&in_buf, total_inds, &mut unpacked);
                    for (i, &g) in unpacked.iter().enumerate() {
                        if keep_inds[i] {
                            stats[snp_idx].observe(g);
                        }
                    }
                }
                snp_idx += 1;
            }
        }
        Layout::SampleMajor => {
            let mut unpacked = vec![0u8; total_snps];
            let mut ind_idx = 0usize;
            while reader.read_record(&mut in_buf)? {
                if ind_idx >= total_inds {
                    break;
                }
                if keep_inds[ind_idx] {
                    codec::unpack(&in_buf, total_snps, &mut unpacked);
                    for (j, &g) in unpacked.iter().enumerate() {
                        if keep_snps[j] {
                            stats[j].observe(g);
                        }
                    }
                }
                ind_idx += 1;
            }
        }
    }

    Ok(stats)
}

pub fn read_input_snp(path: &Path, fmt: Format, numchrom: u32) -> Result<Vec<SnpRow>> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::snp::read(path, numchrom)
        }
        Format::PackedPed => meta::bim::read(path, numchrom),
        Format::Plink2 => meta::pvar::read(path, numchrom),
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
    }
}

pub fn read_input_ind(path: &Path, fmt: Format, familynames: bool) -> Result<Vec<IndRow>> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::ind::read(path)
        }
        Format::PackedPed => meta::fam::read(path, familynames),
        Format::Plink2 => meta::psam::read(path),
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
    }
}

pub fn write_output_snp(path: &Path, fmt: Format, rows: &[SnpRow], numchrom: u32) -> Result<()> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::snp::write(path, rows, numchrom)
        }
        Format::PackedPed => meta::bim::write(path, rows, numchrom),
        Format::Plink2 => meta::pvar::write(path, rows, numchrom),
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
    }
}

pub fn write_output_ind(
    path: &Path,
    fmt: Format,
    rows: &[IndRow],
    outputgroup: bool,
) -> Result<()> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::ind::write(path, rows)
        }
        Format::PackedPed => meta::fam::write(path, rows, outputgroup),
        Format::Plink2 => meta::psam::write(path, rows, outputgroup),
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
    }
}

// ======================================================================
// Geno dispatch
// ======================================================================

/// Public reader dispatch — used by `cmd_export`, `stats`, and other modules
/// that need to open genotype files without going through the full convert
/// pipeline.
pub fn open_reader_pub(
    fmt: Format,
    path: &Path,
    nind: usize,
    nsnp: usize,
) -> Result<Box<dyn GenoReader>> {
    open_reader(fmt, path, nind, nsnp)
}

/// Public writer dispatch.
pub fn open_writer_pub(fmt: Format, path: &Path) -> Result<Box<dyn GenoWriter>> {
    open_writer(fmt, path)
}

fn open_reader(fmt: Format, path: &Path, nind: usize, nsnp: usize) -> Result<Box<dyn GenoReader>> {
    match fmt {
        Format::PackedAncestrymap => Ok(Box::new(
            packed_am::PackedAmReader::open(path, nind, nsnp)
                .with_context(|| format!("open {}", path.display()))?,
        )),
        Format::Eigenstrat => Ok(Box::new(
            eigenstrat::EigenstratReader::open(path, nind, nsnp)
                .with_context(|| format!("open {}", path.display()))?,
        )),
        Format::PackedPed => {
            let bim_path = path.with_extension("bim");
            let flip_02_mask = if bim_path.exists() {
                meta::bim::read_flip_02_mask(&bim_path)
                    .with_context(|| format!("reading {} for the A1=0 flip-mask", bim_path.display()))?
            } else {
                Vec::new()
            };
            Ok(Box::new(
                packed_ped::PackedPedReader::open_with_flip_mask(path, nind, nsnp, flip_02_mask)
                    .with_context(|| format!("open {}", path.display()))?,
            ))
        }
        Format::Plink2 => {
            // PGEN needs the multiallelic keep-mask so the genotype stream
            // aligns with the biallelic-only SnpRows produced by `pvar::read`.
            // Build it from the sibling .pvar (same fileset stem) and thread it
            // in — the pipeline owns the mask, the codec just consumes it.
            let pvar_path = path.with_extension("pvar");
            let keep_ondisk = meta::pvar::keep_ondisk_mask(&pvar_path).with_context(|| {
                format!("reading {} for the multiallelic keep-mask", pvar_path.display())
            })?;
            Ok(Box::new(
                pgen::PgenReader::open(path, keep_ondisk)
                    .with_context(|| format!("open {}", path.display()))?,
            ))
        }
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
        Format::Tgeno => Ok(Box::new(
            tgeno::TgenoReader::open(path, nind, nsnp)
                .with_context(|| format!("open {}", path.display()))?,
        )),
        Format::Ancestrymap => bail!("ANCESTRYMAP sparse reader not implemented"),
    }
}

fn open_writer(fmt: Format, path: &Path) -> Result<Box<dyn GenoWriter>> {
    match fmt {
        Format::PackedAncestrymap => Ok(Box::new(packed_am::PackedAmWriter::create(path)?)),
        Format::Eigenstrat => Ok(Box::new(eigenstrat::EigenstratWriter::create(path)?)),
        Format::PackedPed => Ok(Box::new(packed_ped::PackedPedWriter::create(path)?)),
        Format::Plink2 => Ok(Box::new(pgen::PgenWriter::create(path)?)),
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
        Format::Tgeno => Ok(Box::new(tgeno::TgenoWriter::create(path)?)),
        Format::Ancestrymap => bail!("ANCESTRYMAP sparse writer not implemented"),
    }
}

// ======================================================================
// Streaming loop (same layout both sides)
// ======================================================================

fn stream_same_layout(
    reader: &mut dyn GenoReader,
    writer: &mut dyn GenoWriter,
    keep_snps: &[bool],
    keep_inds: &[bool],
    total_inds: usize,
    kept_inds: usize,
) -> Result<()> {
    match reader.layout() {
        Layout::SnpMajor => {
            stream_snp_major(reader, writer, keep_snps, keep_inds, total_inds, kept_inds)
        }
        Layout::SampleMajor => stream_sample_major(reader, writer, keep_snps, keep_inds),
    }
}

/// SnpMajor → SnpMajor. One record per SNP; filter samples inside the record.
fn stream_snp_major(
    reader: &mut dyn GenoReader,
    writer: &mut dyn GenoWriter,
    keep_snps: &[bool],
    keep_inds: &[bool],
    total_inds: usize,
    kept_inds: usize,
) -> Result<()> {
    let in_rec_bytes = reader.record_bytes();
    let out_rec_bytes = (kept_inds * 2 + 7) / 8;
    let mut in_buf = vec![0u8; in_rec_bytes];
    let mut out_buf = vec![0u8; out_rec_bytes];

    let all_kept = kept_inds == total_inds;
    let mut unpacked = vec![0u8; total_inds];
    let mut projected = vec![0u8; kept_inds];

    let mut snp_idx = 0usize;
    let mut written = 0usize;
    let t = Instant::now();

    while reader.read_record(&mut in_buf)? {
        let keep = keep_snps[snp_idx];
        snp_idx += 1;
        if !keep {
            continue;
        }

        if all_kept {
            writer.write_record(&in_buf)?;
        } else {
            codec::unpack(&in_buf, total_inds, &mut unpacked);
            let mut k = 0;
            for (i, &ki) in keep_inds.iter().enumerate() {
                if ki {
                    projected[k] = unpacked[i];
                    k += 1;
                }
            }
            for b in out_buf.iter_mut() {
                *b = 0;
            }
            codec::pack(&projected, &mut out_buf);
            writer.write_record(&out_buf)?;
        }
        written += 1;
    }

    log::info!("streamed {written} SNP records in {:.2?}", t.elapsed());
    Ok(())
}

/// SampleMajor → SampleMajor (TGENO → TGENO). One record per sample;
/// filter SNPs inside the record. Skip whole record if its sample is
/// filtered out.
fn stream_sample_major(
    reader: &mut dyn GenoReader,
    writer: &mut dyn GenoWriter,
    keep_snps: &[bool],
    keep_inds: &[bool],
) -> Result<()> {
    let in_rec_bytes = reader.record_bytes();
    let total_snps = reader.nsnp();
    let kept_snps_count = keep_snps.iter().filter(|&&k| k).count();
    let out_rec_bytes = (kept_snps_count * 2 + 7) / 8;

    let mut in_buf = vec![0u8; in_rec_bytes];
    let mut out_buf = vec![0u8; out_rec_bytes];

    let all_snps_kept = kept_snps_count == total_snps;
    let mut unpacked = vec![0u8; total_snps];
    let mut projected = vec![0u8; kept_snps_count];

    let mut ind_idx = 0usize;
    let mut written = 0usize;
    let t = Instant::now();

    while reader.read_record(&mut in_buf)? {
        let keep = keep_inds[ind_idx];
        ind_idx += 1;
        if !keep {
            continue;
        }

        if all_snps_kept {
            writer.write_record(&in_buf)?;
        } else {
            codec::unpack(&in_buf, total_snps, &mut unpacked);
            let mut k = 0;
            for (j, &ks) in keep_snps.iter().enumerate() {
                if ks {
                    projected[k] = unpacked[j];
                    k += 1;
                }
            }
            for b in out_buf.iter_mut() {
                *b = 0;
            }
            codec::pack(&projected, &mut out_buf);
            writer.write_record(&out_buf)?;
        }
        written += 1;
    }

    log::info!("streamed {written} sample records in {:.2?}", t.elapsed());
    Ok(())
}

/// Cross-layout: SnpMajor ↔ SampleMajor via full-matrix transpose.
///
/// Cross-layout conversion (TGENO/sample-major <-> a SNP-major format) requires
/// transposing the 2-bit-packed matrix. This dispatcher picks the strategy from
/// the `max_mem` budget: if the whole matrix (source + transposed destination)
/// fits, it takes the fast single-pass in-memory path (`_in_memory`, unchanged);
/// otherwise it uses the memory-bounded banded path (`_banded`), which trades a
/// bounded peak for re-streaming the source once per band.
#[allow(clippy::too_many_arguments)]
fn stream_cross_layout(
    reader: &mut dyn GenoReader,
    writer: &mut dyn GenoWriter,
    in_fmt: Format,
    geno_in: &Path,
    keep_snps: &[bool],
    keep_inds: &[bool],
    total_inds: usize,
    total_snps: usize,
    kept_inds: usize,
    kept_snps_count: usize,
    max_mem: u64,
) -> Result<()> {
    let (src_rows, src_cols) = match reader.layout() {
        Layout::SnpMajor => (kept_snps_count, kept_inds),
        Layout::SampleMajor => (kept_inds, kept_snps_count),
    };
    // Destination rows have `src_rows` cells each (the transpose swaps axes).
    let dst_row_bytes = (src_rows * 2 + 7) / 8;

    // Largest band width `B` (source columns / destination rows per pass) whose
    // two working strips fit the budget:
    //   src_strip = src_rows * ceil(B*2/8)   (<= src_rows*B/4 + src_rows)
    //   dst_strip = B * dst_row_bytes
    // Solve B * (src_rows/4 + dst_row_bytes) <= max_mem - src_rows.
    let per_col = (src_rows as u64) / 4 + dst_row_bytes as u64 + 1;
    let budget = max_mem.saturating_sub(src_rows as u64);
    let band = ((budget / per_col.max(1)).max(1) as usize).min(src_cols.max(1));
    let n_bands = src_cols.div_ceil(band.max(1));

    if n_bands <= 1 {
        return stream_cross_layout_in_memory(
            reader,
            writer,
            keep_snps,
            keep_inds,
            total_inds,
            kept_inds,
            kept_snps_count,
        );
    }
    stream_cross_layout_banded(
        writer,
        in_fmt,
        geno_in,
        reader.layout(),
        keep_snps,
        keep_inds,
        total_inds,
        total_snps,
        kept_inds,
        kept_snps_count,
        band,
        n_bands,
    )
}

/// Materializes the (filtered) source matrix in memory, transposes, writes.
/// Fast single-pass path used when the matrix fits the `--max-mem` budget.
/// For AADR scale (1.23M SNPs × 17.6k samples, ~1.2 GB per side) this uses
/// ~2.5 GB RAM total; larger jobs take `stream_cross_layout_banded` instead.
fn stream_cross_layout_in_memory(
    reader: &mut dyn GenoReader,
    writer: &mut dyn GenoWriter,
    keep_snps: &[bool],
    keep_inds: &[bool],
    total_inds: usize,
    kept_inds: usize,
    kept_snps_count: usize,
) -> Result<()> {
    let t_all = Instant::now();

    // Canonical matrix convention for transpose: rows × cols in cells.
    // From reader's perspective, `rows` = records emitted, `cols` = cells per record.
    let (src_rows, src_cols) = match reader.layout() {
        Layout::SnpMajor => (kept_snps_count, kept_inds), // keep rows=SNPs, cols=inds
        Layout::SampleMajor => (kept_inds, kept_snps_count), // rows=inds, cols=SNPs
    };
    let src_row_bytes = (src_cols * 2 + 7) / 8;

    // --- Read phase: materialize kept source records, projecting the
    //     other axis per record.
    let t_read = Instant::now();
    let mut src_matrix = vec![0u8; src_rows * src_row_bytes];
    let in_rec_bytes = reader.record_bytes();
    let mut in_buf = vec![0u8; in_rec_bytes];

    match reader.layout() {
        Layout::SnpMajor => {
            // Each record = one SNP × total_inds. Keep SNPs per keep_snps;
            // project samples per keep_inds.
            let all_inds = kept_inds == total_inds;
            let mut unpacked = vec![0u8; total_inds];
            let mut projected = vec![0u8; kept_inds];

            let mut snp_idx = 0usize;
            let mut out_row = 0usize;
            while reader.read_record(&mut in_buf)? {
                let keep = keep_snps[snp_idx];
                snp_idx += 1;
                if !keep {
                    continue;
                }

                let row_slice =
                    &mut src_matrix[out_row * src_row_bytes..(out_row + 1) * src_row_bytes];
                if all_inds {
                    row_slice.copy_from_slice(&in_buf[..src_row_bytes]);
                } else {
                    codec::unpack(&in_buf, total_inds, &mut unpacked);
                    let mut k = 0;
                    for (i, &ki) in keep_inds.iter().enumerate() {
                        if ki {
                            projected[k] = unpacked[i];
                            k += 1;
                        }
                    }
                    for b in row_slice.iter_mut() {
                        *b = 0;
                    }
                    codec::pack(&projected, row_slice);
                }
                out_row += 1;
            }
            debug_assert_eq!(out_row, kept_snps_count);
        }
        Layout::SampleMajor => {
            // Each record = one sample × nsnp. Keep samples per keep_inds;
            // project SNPs per keep_snps.
            let total_snps = reader.nsnp();
            let all_snps = kept_snps_count == total_snps;
            let mut unpacked = vec![0u8; total_snps];
            let mut projected = vec![0u8; kept_snps_count];

            let mut ind_idx = 0usize;
            let mut out_row = 0usize;
            while reader.read_record(&mut in_buf)? {
                let keep = keep_inds[ind_idx];
                ind_idx += 1;
                if !keep {
                    continue;
                }

                let row_slice =
                    &mut src_matrix[out_row * src_row_bytes..(out_row + 1) * src_row_bytes];
                if all_snps {
                    row_slice.copy_from_slice(&in_buf[..src_row_bytes]);
                } else {
                    codec::unpack(&in_buf, total_snps, &mut unpacked);
                    let mut k = 0;
                    for (j, &ks) in keep_snps.iter().enumerate() {
                        if ks {
                            projected[k] = unpacked[j];
                            k += 1;
                        }
                    }
                    for b in row_slice.iter_mut() {
                        *b = 0;
                    }
                    codec::pack(&projected, row_slice);
                }
                out_row += 1;
            }
            debug_assert_eq!(out_row, kept_inds);
        }
    }
    log::info!(
        "materialized {}x{} source matrix ({} MB) in {:.2?}",
        src_rows,
        src_cols,
        src_matrix.len() / 1_048_576,
        t_read.elapsed()
    );

    // --- Transpose.
    let dst_rows = src_cols;
    let dst_cols = src_rows;
    let dst_row_bytes = (dst_cols * 2 + 7) / 8;
    let mut dst_matrix = vec![0u8; dst_rows * dst_row_bytes];

    let t_transpose = Instant::now();
    crate::transpose::transpose_packed(&src_matrix, src_rows, src_cols, &mut dst_matrix)?;
    log::info!("transposed in {:.2?}", t_transpose.elapsed());

    // Free source matrix memory before writing.
    drop(src_matrix);

    // --- Write phase.
    let t_write = Instant::now();
    for r in 0..dst_rows {
        let row = &dst_matrix[r * dst_row_bytes..(r + 1) * dst_row_bytes];
        writer.write_record(row)?;
    }
    log::info!("wrote {} records in {:.2?}", dst_rows, t_write.elapsed());

    log::info!("cross-layout total: {:.2?}", t_all.elapsed());
    Ok(())
}

/// Memory-bounded cross-layout transpose. Emits the destination in bands of
/// `band` rows (= `band` source columns). For each band it re-opens the source,
/// streams every record while packing only that band's columns into a small
/// strip, transposes the strip with the shared `transpose_packed` kernel, and
/// writes the band's records in order. Peak RAM ≈ `src_strip + dst_strip`,
/// which the caller sized to the `--max-mem` budget.
///
/// Cost: the source is re-streamed `n_bands` times (read + unpack/project
/// amplification). The transpose math and codec are reused verbatim from the
/// single-pass path, so this stays out of the genotype-corruption blast radius;
/// correctness is pinned by differential tests against the in-memory path and
/// the plink2/PAM oracles.
#[allow(clippy::too_many_arguments)]
fn stream_cross_layout_banded(
    writer: &mut dyn GenoWriter,
    in_fmt: Format,
    geno_in: &Path,
    layout: Layout,
    keep_snps: &[bool],
    keep_inds: &[bool],
    total_inds: usize,
    total_snps: usize,
    kept_inds: usize,
    kept_snps_count: usize,
    band: usize,
    n_bands: usize,
) -> Result<()> {
    let t_all = Instant::now();
    let (src_rows, src_cols) = match layout {
        Layout::SnpMajor => (kept_snps_count, kept_inds),
        Layout::SampleMajor => (kept_inds, kept_snps_count),
    };
    let dst_row_bytes = (src_rows * 2 + 7) / 8;

    log::info!(
        "cross-layout: banded transpose of {}x{} matrix in {} band(s) of <= {} col(s) \
         (re-streams source {} time(s), est. peak {} MB)",
        src_rows,
        src_cols,
        n_bands,
        band,
        n_bands,
        (src_rows * ((band * 2 + 7) / 8) + band * dst_row_bytes) / 1_048_576,
    );

    let mut col_start = 0usize;
    let mut dst_written = 0usize;
    while col_start < src_cols {
        let this_band = band.min(src_cols - col_start);
        let this_row_bytes = (this_band * 2 + 7) / 8;
        // src_strip: src_rows records, each holding this band's `this_band` cells.
        let mut src_strip = vec![0u8; src_rows * this_row_bytes];

        // One fresh streaming pass over the whole source for this band.
        let mut reader = open_reader(in_fmt, geno_in, total_inds, total_snps)?;
        let mut in_buf = vec![0u8; reader.record_bytes()];
        let mut out_row = 0usize;

        match layout {
            Layout::SnpMajor => {
                // Records are SNPs; keep per keep_snps, project samples per keep_inds,
                // then take the band's slice of the kept-sample axis.
                let mut unpacked = vec![0u8; total_inds];
                let mut projected = vec![0u8; kept_inds];
                let mut snp_idx = 0usize;
                while reader.read_record(&mut in_buf)? {
                    let keep = keep_snps[snp_idx];
                    snp_idx += 1;
                    if !keep {
                        continue;
                    }
                    codec::unpack(&in_buf, total_inds, &mut unpacked);
                    let mut k = 0;
                    for (i, &ki) in keep_inds.iter().enumerate() {
                        if ki {
                            projected[k] = unpacked[i];
                            k += 1;
                        }
                    }
                    let row_slice =
                        &mut src_strip[out_row * this_row_bytes..(out_row + 1) * this_row_bytes];
                    codec::pack(&projected[col_start..col_start + this_band], row_slice);
                    out_row += 1;
                }
            }
            Layout::SampleMajor => {
                // Records are samples; keep per keep_inds, project SNPs per keep_snps,
                // then take the band's slice of the kept-SNP axis.
                let mut unpacked = vec![0u8; total_snps];
                let mut projected = vec![0u8; kept_snps_count];
                let mut ind_idx = 0usize;
                while reader.read_record(&mut in_buf)? {
                    let keep = keep_inds[ind_idx];
                    ind_idx += 1;
                    if !keep {
                        continue;
                    }
                    codec::unpack(&in_buf, total_snps, &mut unpacked);
                    let mut k = 0;
                    for (j, &ks) in keep_snps.iter().enumerate() {
                        if ks {
                            projected[k] = unpacked[j];
                            k += 1;
                        }
                    }
                    let row_slice =
                        &mut src_strip[out_row * this_row_bytes..(out_row + 1) * this_row_bytes];
                    codec::pack(&projected[col_start..col_start + this_band], row_slice);
                    out_row += 1;
                }
            }
        }
        debug_assert_eq!(
            out_row, src_rows,
            "banded pass produced {out_row} source rows, expected {src_rows}"
        );

        // Transpose this band: [src_rows x this_band] -> [this_band x src_rows],
        // then emit the band's destination records in order.
        let mut dst_strip = vec![0u8; this_band * dst_row_bytes];
        crate::transpose::transpose_packed(&src_strip, src_rows, this_band, &mut dst_strip)?;
        for r in 0..this_band {
            writer.write_record(&dst_strip[r * dst_row_bytes..(r + 1) * dst_row_bytes])?;
            dst_written += 1;
        }
        col_start += this_band;
    }
    debug_assert_eq!(dst_written, src_cols);
    log::info!("cross-layout (banded) total: {:.2?}", t_all.elapsed());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        compute_missing_counts, default_max_mem, detect_total_memory, parse_meminfo_total,
        run_convert, validate_missingness_threshold, ConvertConfig,
    };
    use crate::format::Format;
    use crate::geno::{codec, GenoReader, Layout};
    use anyhow::Result;

    struct MockReader {
        layout: Layout,
        nind: usize,
        nsnp: usize,
        records: Vec<Vec<u8>>,
        idx: usize,
    }

    impl GenoReader for MockReader {
        fn nind(&self) -> usize {
            self.nind
        }
        fn nsnp(&self) -> usize {
            self.nsnp
        }
        fn layout(&self) -> Layout {
            self.layout
        }
        fn read_record(&mut self, dst: &mut [u8]) -> Result<bool> {
            if self.idx >= self.records.len() {
                return Ok(false);
            }
            dst.copy_from_slice(&self.records[self.idx]);
            self.idx += 1;
            Ok(true)
        }
    }

    #[test]
    fn meminfo_total_parses_kb_to_bytes() {
        let sample = "MemTotal:       16106280 kB\nMemFree:  1234 kB\nMemAvailable: 8000000 kB\n";
        assert_eq!(parse_meminfo_total(sample), Some(16106280u64 * 1024));
        // Missing / malformed -> None (callers fall back).
        assert_eq!(parse_meminfo_total("MemFree: 100 kB\n"), None);
        assert_eq!(parse_meminfo_total("MemTotal: notanumber kB\n"), None);
    }

    #[test]
    fn default_max_mem_is_half_ram_and_floored() {
        // On the test host RAM is detectable; the default must be positive,
        // at least the 512 MiB floor, and no more than total RAM.
        let budget = default_max_mem();
        assert!(budget >= 512 * 1024 * 1024, "budget below floor: {budget}");
        if let Some(total) = detect_total_memory() {
            assert!(budget <= total, "budget {budget} exceeds total RAM {total}");
            // Half of RAM (unless the floor kicked in on a tiny machine).
            assert!(budget == (total / 2).max(512 * 1024 * 1024));
        }
    }

    fn pack(vals: &[u8]) -> Vec<u8> {
        let mut out = vec![0u8; (vals.len() * 2 + 7) / 8];
        codec::pack(vals, &mut out);
        out
    }

    #[test]
    fn validate_threshold_bounds() {
        assert!(validate_missingness_threshold("geno", Some(0.0)).is_ok());
        assert!(validate_missingness_threshold("geno", Some(1.0)).is_ok());
        assert!(validate_missingness_threshold("geno", Some(-0.1)).is_err());
        assert!(validate_missingness_threshold("mind", Some(1.1)).is_err());
    }

    #[test]
    fn missing_counts_snp_major_respects_masks() {
        // 3 SNP x 2 samples, SNP-major records.
        // s0: [0,3], s1:[3,3], s2:[2,0]
        let records = vec![pack(&[0, 3]), pack(&[3, 3]), pack(&[2, 0])];
        let mut r = MockReader {
            layout: Layout::SnpMajor,
            nind: 2,
            nsnp: 3,
            records,
            idx: 0,
        };
        let keep_snps = vec![true, true, false];
        let keep_inds = vec![true, true];
        let (snp_miss, ind_miss) =
            compute_missing_counts(&mut r, &keep_snps, &keep_inds, 2, 3).unwrap();
        assert_eq!(snp_miss, vec![1, 2, 0]);
        assert_eq!(ind_miss, vec![1, 2]);
    }

    #[test]
    fn missing_counts_sample_major_respects_masks() {
        // 2 samples x 3 SNPs, sample-major records.
        // i0: [0,3,2], i1:[3,3,0]
        let records = vec![pack(&[0, 3, 2]), pack(&[3, 3, 0])];
        let mut r = MockReader {
            layout: Layout::SampleMajor,
            nind: 2,
            nsnp: 3,
            records,
            idx: 0,
        };
        let keep_snps = vec![true, true, false];
        let keep_inds = vec![true, true];
        let (snp_miss, ind_miss) =
            compute_missing_counts(&mut r, &keep_snps, &keep_inds, 2, 3).unwrap();
        assert_eq!(snp_miss, vec![1, 2, 0]);
        assert_eq!(ind_miss, vec![1, 2]);
    }

    // B-012: --mind must be applied BEFORE --geno (PLINK order). Construct a
    // case where the two orders disagree: sample S1's missingness is entirely
    // on the high-missingness SNPs that --geno would remove.
    //   geno-first (wrong): drop rs2/rs3, then S1 looks fully called → 2 SNP, 2 ind
    //   mind-first (right): S1 is 2/4 missing → dropped, then all 4 SNPs survive
    //                       over the remaining sample → 4 SNP, 1 ind
    #[test]
    fn mind_applied_before_geno_like_plink() {
        let dir = tempfile::tempdir().unwrap();
        let geno_in = dir.path().join("in.geno");
        let snp_in = dir.path().join("in.snp");
        let ind_in = dir.path().join("in.ind");
        // EIGENSTRAT: SNP-major, one line per SNP, one char per sample (0/1/2/9).
        std::fs::write(&geno_in, "00\n00\n09\n09\n").unwrap();
        std::fs::write(
            &snp_in,
            "rs0 1 0.0 1 A C\nrs1 1 0.0 2 A C\nrs2 1 0.0 3 A C\nrs3 1 0.0 4 A C\n",
        )
        .unwrap();
        std::fs::write(&ind_in, "S0 U Pop\nS1 U Pop\n").unwrap();

        let geno_out = dir.path().join("out.geno");
        let snp_out = dir.path().join("out.snp");
        let ind_out = dir.path().join("out.ind");
        let cfg = ConvertConfig {
            geno_in,
            snp_in,
            ind_in,
            out_fmt: Format::Eigenstrat,
            geno_out,
            snp_out: snp_out.clone(),
            ind_out: ind_out.clone(),
            badsnp: None,
            snps: None,
            poplist: None,
            keep: None,
            remove: None,
            chrom: None,
            lopos: None,
            hipos: None,
            noxdata: false,
            max_miss_snp: Some(0.4),
            max_miss_ind: Some(0.4),
            maf: None,
            max_maf: None,
            hwe: None,
            numchrom: 22,
            hashcheck: false,
            familynames: true,
            outputgroup: false,
            max_mem: default_max_mem(),
        };
        run_convert(&cfg).unwrap();

        let n_snp = std::fs::read_to_string(&snp_out).unwrap().lines().count();
        let n_ind = std::fs::read_to_string(&ind_out).unwrap().lines().count();
        assert_eq!(n_ind, 1, "mind must run first and drop S1 (got {n_ind} samples)");
        assert_eq!(n_snp, 4, "all SNPs survive over the post-mind sample (got {n_snp})");
    }
}

pub fn resolve_paths(
    prefix: Option<String>,
    geno: Option<PathBuf>,
    snp: Option<PathBuf>,
    ind: Option<PathBuf>,
    out_format: Option<Format>,
    is_output: bool,
) -> Result<(PathBuf, PathBuf, PathBuf)> {
    if let Some(p) = prefix {
        let p_path = PathBuf::from(p);
        let g = match geno {
            Some(g) => g,
            None if is_output => {
                let (gext, _, _) = out_format
                    .expect("out_format required for output")
                    .default_output_extensions();
                p_path.with_extension(gext)
            }
            None => {
                // Input prefix: pick whichever genotype file actually exists,
                // preferring PLINK2 → PLINK1 → EIGENSTRAT/PAM. If none exists,
                // say exactly what was looked for instead of failing later with
                // a bare "No such file or directory".
                ["pgen", "bed", "geno"]
                    .iter()
                    .map(|ext| p_path.with_extension(ext))
                    .find(|c| c.exists())
                    .ok_or_else(|| {
                        let b = p_path.display();
                        anyhow::anyhow!(
                            "no genotype file found for input prefix '{b}'\n  \
                             looked for: {b}.pgen, {b}.bed, {b}.geno\n  \
                             hint: check the prefix, or pass the file directly with --in-geno"
                        )
                    })?
            }
        };

        let derive = |ext_pgen: &str, ext_bed: &str, ext_default: &str| {
            let gext = g.extension().and_then(|e| e.to_str());
            match gext {
                Some("pgen") => p_path.with_extension(ext_pgen),
                Some("bed") => p_path.with_extension(ext_bed),
                _ => p_path.with_extension(ext_default),
            }
        };

        let s = match snp {
            Some(s) => s,
            None if is_output => {
                let (_, sext, _) = out_format.unwrap().default_output_extensions();
                p_path.with_extension(sext)
            }
            None => derive("pvar", "bim", "snp"),
        };
        let i = match ind {
            Some(i) => i,
            None if is_output => {
                let (_, _, iext) = out_format.unwrap().default_output_extensions();
                p_path.with_extension(iext)
            }
            None => derive("psam", "fam", "ind"),
        };

        if !is_output {
            check_input_exists(&g, &s, &i)?;
        }
        Ok((g, s, i))
    } else if is_output {
        // No prefix on the output side: the user must name every output file.
        let g = geno.ok_or_else(|| {
            anyhow::anyhow!(
                "no output location given: pass -o/--out-prefix <prefix>, or name every \
                 output file with --out-geno/--out-snp/--out-ind"
            )
        })?;
        let s = snp.ok_or_else(|| anyhow::anyhow!("missing --out-snp (or use -o/--out-prefix)"))?;
        let i = ind.ok_or_else(|| anyhow::anyhow!("missing --out-ind (or use -o/--out-prefix)"))?;
        Ok((g, s, i))
    } else {
        let g = geno.ok_or_else(|| {
            anyhow::anyhow!(
                "no input given: pass -i/--in-prefix <prefix>, or name the files with \
                 --in-geno/--in-snp/--in-ind"
            )
        })?;
        let s = snp
            .ok_or_else(|| anyhow::anyhow!("missing --in-snp (or use -i/--in-prefix)"))?;
        let i = ind
            .ok_or_else(|| anyhow::anyhow!("missing --in-ind (or use -i/--in-prefix)"))?;
        check_input_exists(&g, &s, &i)?;
        Ok((g, s, i))
    }
}

/// Verify each resolved input file exists, naming the missing one specifically
/// so the user gets a clear "this file is missing" rather than a downstream
/// "No such file or directory (os error 2)".
fn check_input_exists(geno: &Path, snp: &Path, ind: &Path) -> Result<()> {
    for (path, flag, what) in [
        (geno, "--in-geno", "genotype"),
        (snp, "--in-snp", "SNP"),
        (ind, "--in-ind", "individual/sample"),
    ] {
        if !path.exists() {
            bail!(
                "{what} input file not found: '{}'\n  \
                 hint: check the path (or the {flag} argument / the -i prefix)",
                path.display()
            );
        }
    }
    Ok(())
}
