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
use crate::geno::{eigenstrat, packed_am, packed_ped, tgeno};
use crate::meta::{self, IndRow, SnpRow};
use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};
use std::time::Instant;

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
        bail!("all SNPs filtered out");
    }
    if kept_ind_count == 0 {
        bail!("all samples filtered out");
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
            &keep_snps,
            &keep_inds,
            ind_rows.len(),
            kept_ind_count,
            kept_snp_count,
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
    // First phase: apply per-SNP missingness (geno) on current sample keep-mask.
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

    // Second phase: apply per-sample missingness (mind) on the SNP mask after geno.
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
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
    }
}

pub fn read_input_ind(path: &Path, fmt: Format, familynames: bool) -> Result<Vec<IndRow>> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::ind::read(path)
        }
        Format::PackedPed => meta::fam::read(path, familynames),
        Format::Ped => bail!("PED text format not supported (use PACKEDPED)"),
    }
}

pub fn write_output_snp(path: &Path, fmt: Format, rows: &[SnpRow], numchrom: u32) -> Result<()> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::snp::write(path, rows, numchrom)
        }
        Format::PackedPed => meta::bim::write(path, rows, numchrom),
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
        Format::PackedPed => Ok(Box::new(
            packed_ped::PackedPedReader::open(path, nind, nsnp)
                .with_context(|| format!("open {}", path.display()))?,
        )),
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
/// Materializes the (filtered) source matrix in memory, transposes, writes.
/// For AADR scale (1.23M SNPs × 17.6k samples, ~1.2 GB per side) this uses
/// ~2.5 GB RAM total. Acceptable on modern systems; if someone complains,
/// we can add streaming tiled transpose later.
fn stream_cross_layout(
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

#[cfg(test)]
mod tests {
    use super::{compute_missing_counts, validate_missingness_threshold};
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
        let g = geno.unwrap_or_else(|| {
            if is_output {
                let (gext, _, _) = out_format
                    .expect("out_format required for output")
                    .default_output_extensions();
                return p_path.with_extension(gext);
            }
            // Check for .bed first, then .geno
            let bed = p_path.with_extension("bed");
            if bed.exists() {
                bed
            } else {
                p_path.with_extension("geno")
            }
        });

        let s = snp.unwrap_or_else(|| {
            if is_output {
                let (_, sext, _) = out_format
                    .expect("out_format required for output")
                    .default_output_extensions();
                return p_path.with_extension(sext);
            }
            if g.extension().and_then(|e| e.to_str()) == Some("bed") {
                p_path.with_extension("bim")
            } else {
                p_path.with_extension("snp")
            }
        });

        let i = ind.unwrap_or_else(|| {
            if is_output {
                let (_, _, iext) = out_format
                    .expect("out_format required for output")
                    .default_output_extensions();
                return p_path.with_extension(iext);
            }
            if g.extension().and_then(|e| e.to_str()) == Some("bed") {
                p_path.with_extension("fam")
            } else {
                p_path.with_extension("ind")
            }
        });

        Ok((g, s, i))
    } else {
        let g = geno.ok_or_else(|| {
            anyhow::anyhow!("missing genotype input (provide --in-geno or --in-prefix)")
        })?;
        let s = snp.ok_or_else(|| {
            anyhow::anyhow!("missing SNP input (provide --in-snp or --in-prefix)")
        })?;
        let i = ind.ok_or_else(|| {
            anyhow::anyhow!("missing individual input (provide --in-ind or --in-prefix)")
        })?;
        Ok((g, s, i))
    }
}
