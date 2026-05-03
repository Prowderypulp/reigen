//! `reigen stats` — compute QC and summary statistics.
//!
//! Single-pass streaming computation of per-SNP and per-sample statistics.
//! Optional pairwise IBS matrix (O(nind²) memory).

pub mod ibs;
pub mod per_sample;
pub mod per_snp;

use crate::format::{self, Format};
use crate::geno::codec;
use crate::meta::{self, IndRow, SnpRow};
use anyhow::{bail, Context, Result};
use clap::Args;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args, Debug)]
pub struct StatsArgs {
    /// Input prefix (derives .geno/.snp/.ind or .bed/.bim/.fam)
    #[arg(short = 'i', long)]
    pub in_prefix: Option<String>,

    /// Input genotype file
    #[arg(long)]
    pub geno: Option<PathBuf>,

    /// Input SNP file
    #[arg(long)]
    pub snp: Option<PathBuf>,

    /// Input individual/family file
    #[arg(long)]
    pub ind: Option<PathBuf>,

    /// Output prefix for stats files
    #[arg(short = 'o', long)]
    pub out_prefix: String,

    /// Disable per-SNP statistics (output: <prefix>.snp_stats.tsv)
    #[arg(long)]
    pub no_per_snp: bool,

    /// Disable per-sample statistics (output: <prefix>.sample_stats.tsv)
    #[arg(long)]
    pub no_per_sample: bool,

    /// Compute pairwise IBS distance matrix (output: <prefix>.ibs.tsv, <prefix>.dst.tsv).
    /// Warning: requires O(nind²) memory.
    #[arg(long)]
    pub ibs: bool,

    /// Number of autosomes (default 22)
    #[arg(long, default_value_t = 22)]
    pub numchrom: u32,

    /// Treat PLINK .fam FID column as pop label
    #[arg(long)]
    pub no_familynames: bool,

    /// Population list to keep
    #[arg(long)]
    pub poplist: Option<PathBuf>,

    /// Restrict to a single chromosome
    #[arg(long)]
    pub chrom: Option<i32>,
}

pub fn run_stats(args: StatsArgs) -> Result<()> {
    let t0 = Instant::now();

    // Resolve input paths
    let (geno_in, snp_in, ind_in) = resolve_input_paths(
        args.in_prefix.as_deref(),
        args.geno.as_deref(),
        args.snp.as_deref(),
        args.ind.as_deref(),
    )?;

    let in_fmt = format::infer_input_format(&geno_in).context("inferring input format")?;
    log::info!("input format: {in_fmt:?}");

    let numchrom = args.numchrom;

    // Load metadata
    let snp_rows = read_snp(&snp_in, in_fmt, numchrom)?;
    let ind_rows = read_ind(&ind_in, in_fmt, !args.no_familynames)?;
    log::info!(
        "metadata: {} SNPs, {} samples",
        snp_rows.len(),
        ind_rows.len()
    );

    // Apply filters
    let pop_keep = args
        .poplist
        .as_deref()
        .map(crate::filter::load_pop_keep)
        .transpose()?;

    let x_chrom = u8::try_from(numchrom + 1).context("numchrom too large")?;
    let y_chrom = u8::try_from(numchrom + 2).context("numchrom too large")?;
    let mt_chrom = u8::try_from(numchrom + 3).context("numchrom too large")?;
    let xy_chrom = u8::try_from(numchrom + 4).context("numchrom too large")?;
    let chrom_filter = if let Some(c) = args.chrom {
        let chrom = u8::try_from(c).context("`--chrom` must be in range 0..=255")?;
        Some(crate::filter::ChromFilter::Single(chrom))
    } else {
        None
    };

    let snp_filter = crate::filter::SnpFilter {
        bad: None,
        snp_keep: None,
        chrom: chrom_filter,
        lopos: None,
        hipos: None,
        noxdata: false,
        x_chrom,
        y_chrom,
        mt_chrom,
        xy_chrom,
    };
    let ind_filter = crate::filter::IndFilter {
        pop_keep: pop_keep.as_ref(),
        sample_keep: None,
        sample_remove: None,
    };

    let keep_snps: Vec<bool> = snp_rows.iter().map(|s| snp_filter.keep(s)).collect();
    let keep_inds: Vec<bool> = ind_rows.iter().map(|i| ind_filter.keep(i)).collect();
    let kept_snp_count = keep_snps.iter().filter(|&&k| k).count();
    let kept_ind_count = keep_inds.iter().filter(|&&k| k).count();

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

    // Filtered metadata for output
    let out_snps: Vec<SnpRow> = snp_rows
        .iter()
        .zip(keep_snps.iter())
        .filter_map(|(s, &k)| if k { Some(s.clone()) } else { None })
        .collect();
    let out_inds: Vec<IndRow> = ind_rows
        .iter()
        .zip(keep_inds.iter())
        .filter_map(|(i, &k)| if k { Some(i.clone()) } else { None })
        .collect();

    // Initialize accumulators
    let mut snp_stats = vec![per_snp::SnpStats::default(); kept_snp_count];
    let mut sample_stats = vec![per_sample::SampleStats::default(); kept_ind_count];
    let mut ibs_acc = if args.ibs {
        log::info!(
            "IBS mode: allocating {} pairs ({:.1} MB)",
            kept_ind_count * (kept_ind_count - 1) / 2,
            (kept_ind_count * (kept_ind_count - 1) / 2 * 8) as f64 / 1_048_576.0
        );
        Some(ibs::IbsAccumulator::new(kept_ind_count))
    } else {
        None
    };

    // Open reader and stream
    let total_inds = ind_rows.len();
    let total_snps = snp_rows.len();
    let mut reader = crate::pipeline::open_reader_pub(in_fmt, &geno_in, total_inds, total_snps)?;

    log::info!("computing statistics (single pass)...");
    let in_rec_bytes = reader.record_bytes();
    let mut in_buf = vec![0u8; in_rec_bytes];

    match reader.layout() {
        crate::geno::Layout::SnpMajor => {
            let mut unpacked = vec![0u8; total_inds];
            let mut projected = vec![0u8; kept_ind_count];
            let mut snp_idx = 0usize;
            let mut out_snp_idx = 0usize;

            while reader.read_record(&mut in_buf)? {
                if snp_idx >= total_snps {
                    break;
                }
                if !keep_snps[snp_idx] {
                    snp_idx += 1;
                    continue;
                }

                codec::unpack(&in_buf, total_inds, &mut unpacked);

                // Project to kept samples
                let mut k = 0usize;
                for (i, &ki) in keep_inds.iter().enumerate() {
                    if ki {
                        projected[k] = unpacked[i];
                        k += 1;
                    }
                }

                // Per-SNP stats
                let ref_freq_for_exp_het;
                {
                    let st = &mut snp_stats[out_snp_idx];
                    for &g in &projected {
                        st.observe(g);
                    }
                    ref_freq_for_exp_het = st.ref_freq();
                }

                // Per-sample stats
                let exp_het = if ref_freq_for_exp_het.is_nan() {
                    0.0
                } else {
                    2.0 * ref_freq_for_exp_het * (1.0 - ref_freq_for_exp_het)
                };
                for (k_idx, &g) in projected.iter().enumerate() {
                    sample_stats[k_idx].observe(g);
                    if g != codec::G_MISSING {
                        sample_stats[k_idx].sum_exp_het += exp_het;
                    }
                }

                // IBS accumulation
                if let Some(ref mut acc) = ibs_acc {
                    acc.observe_snp(&projected);
                }

                out_snp_idx += 1;
                snp_idx += 1;
            }
        }
        crate::geno::Layout::SampleMajor => {
            // Each record is one sample across all SNPs.
            // We collect projected genotypes and then process SNP-by-SNP for
            // F-stat (needs per-SNP ref_freq) and IBS (needs all samples at each SNP).
            let mut unpacked = vec![0u8; total_snps];

            // Materialize the kept-sample × kept-SNP matrix.
            // projected_matrix[kept_sample_idx][kept_snp_idx]
            let mut projected_matrix: Vec<Vec<u8>> = Vec::with_capacity(kept_ind_count);
            let mut ind_idx = 0usize;

            while reader.read_record(&mut in_buf)? {
                if ind_idx >= total_inds {
                    break;
                }
                if keep_inds[ind_idx] {
                    codec::unpack(&in_buf, total_snps, &mut unpacked);
                    // Project to kept SNPs
                    let row: Vec<u8> = unpacked
                        .iter()
                        .enumerate()
                        .filter(|(j, _)| keep_snps[*j])
                        .map(|(_, &g)| g)
                        .collect();
                    projected_matrix.push(row);
                }
                ind_idx += 1;
            }

            // Fill missing samples if file was short
            while projected_matrix.len() < kept_ind_count {
                projected_matrix.push(vec![codec::G_MISSING; kept_snp_count]);
            }

            // Now iterate SNP-by-SNP to accumulate stats correctly.
            let mut snp_col = vec![codec::G_MISSING; kept_ind_count];
            for s in 0..kept_snp_count {
                // Gather this SNP's genotypes across all kept samples
                for (k, sample_row) in projected_matrix.iter().enumerate() {
                    snp_col[k] = sample_row[s];
                }

                // Per-SNP stats
                let ref_freq_for_exp_het;
                {
                    let st = &mut snp_stats[s];
                    for &g in &snp_col {
                        st.observe(g);
                    }
                    ref_freq_for_exp_het = st.ref_freq();
                }

                // Per-sample stats
                let exp_het = if ref_freq_for_exp_het.is_nan() {
                    0.0
                } else {
                    2.0 * ref_freq_for_exp_het * (1.0 - ref_freq_for_exp_het)
                };
                for (k_idx, &g) in snp_col.iter().enumerate() {
                    sample_stats[k_idx].observe(g);
                    if g != codec::G_MISSING {
                        sample_stats[k_idx].sum_exp_het += exp_het;
                    }
                }

                // IBS accumulation
                if let Some(ref mut acc) = ibs_acc {
                    acc.observe_snp(&snp_col);
                }
            }
        }
    }

    log::info!("statistics computed in {:.2?}", t0.elapsed());

    // Write outputs
    let out_prefix = PathBuf::from(&args.out_prefix);

    if !args.no_per_snp {
        let snp_path = out_prefix.with_extension("snp_stats.tsv");
        per_snp::write_snp_stats(&snp_path, &out_snps, &snp_stats)?;
        log::info!("per-SNP stats: {}", snp_path.display());
    }

    if !args.no_per_sample {
        let sample_path = out_prefix.with_extension("sample_stats.tsv");
        per_sample::write_sample_stats(&sample_path, &out_inds, &sample_stats)?;
        log::info!("per-sample stats: {}", sample_path.display());
    }

    if let Some(acc) = ibs_acc {
        let ibs_path = out_prefix.with_extension("ibs.tsv");
        let dist_path = out_prefix.with_extension("dst.tsv");
        let ibs_mat = acc.ibs_matrix();
        let dist_mat = acc.distance_matrix();
        ibs::write_ibs_matrix(&ibs_path, &out_inds, &ibs_mat, false)?;
        ibs::write_ibs_matrix(&dist_path, &out_inds, &dist_mat, true)?;
        log::info!("IBS matrix: {}", ibs_path.display());
        log::info!("distance matrix: {}", dist_path.display());
    }

    log::info!("stats done in {:.2?}", t0.elapsed());
    Ok(())
}

fn resolve_input_paths(
    prefix: Option<&str>,
    geno: Option<&std::path::Path>,
    snp: Option<&std::path::Path>,
    ind: Option<&std::path::Path>,
) -> Result<(PathBuf, PathBuf, PathBuf)> {
    if let Some(p) = prefix {
        let p_path = PathBuf::from(p);
        let bed = p_path.with_extension("bed");
        let is_bed = bed.exists();
        let g = geno.map(PathBuf::from).unwrap_or_else(|| {
            if is_bed {
                bed.clone()
            } else {
                p_path.with_extension("geno")
            }
        });
        let s = snp.map(PathBuf::from).unwrap_or_else(|| {
            if is_bed {
                p_path.with_extension("bim")
            } else {
                p_path.with_extension("snp")
            }
        });
        let i = ind.map(PathBuf::from).unwrap_or_else(|| {
            if is_bed {
                p_path.with_extension("fam")
            } else {
                p_path.with_extension("ind")
            }
        });
        Ok((g, s, i))
    } else {
        let g = geno
            .map(PathBuf::from)
            .ok_or_else(|| anyhow::anyhow!("missing --geno or --in-prefix"))?;
        let s = snp
            .map(PathBuf::from)
            .ok_or_else(|| anyhow::anyhow!("missing --snp or --in-prefix"))?;
        let i = ind
            .map(PathBuf::from)
            .ok_or_else(|| anyhow::anyhow!("missing --ind or --in-prefix"))?;
        Ok((g, s, i))
    }
}

fn read_snp(path: &std::path::Path, fmt: Format, numchrom: u32) -> Result<Vec<SnpRow>> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::snp::read(path, numchrom)
        }
        Format::PackedPed => meta::bim::read(path, numchrom),
        Format::Ped => bail!("PED text format not supported"),
    }
}

fn read_ind(path: &std::path::Path, fmt: Format, familynames: bool) -> Result<Vec<IndRow>> {
    match fmt {
        Format::Eigenstrat | Format::PackedAncestrymap | Format::Ancestrymap | Format::Tgeno => {
            meta::ind::read(path)
        }
        Format::PackedPed => meta::fam::read(path, familynames),
        Format::Ped => bail!("PED text format not supported"),
    }
}
