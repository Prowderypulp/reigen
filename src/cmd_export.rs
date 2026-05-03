//! `reigen export` — export internal formats to VCF.
//!
//! Reads any supported genotype format and writes a VCF 4.3 file with
//! biallelic SNP GT fields only.

use crate::format::{self, Format};
use crate::geno::codec;
use crate::meta::{self, IndRow, SnpRow};
use crate::vcf;
use anyhow::{bail, Context, Result};
use clap::Args;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args, Debug)]
pub struct ExportArgs {
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

    /// Output VCF file path
    #[arg(short = 'o', long)]
    pub out: PathBuf,

    /// Prefix chromosomes with "chr" (e.g., "chr1" instead of "1")
    #[arg(long)]
    pub chr_prefix: bool,

    /// Population list to keep
    #[arg(long)]
    pub poplist: Option<PathBuf>,

    /// Restrict to a single chromosome
    #[arg(long)]
    pub chrom: Option<i32>,

    /// Low BP position (for range filter)
    #[arg(long)]
    pub from_bp: Option<u64>,

    /// High BP position (for range filter)
    #[arg(long)]
    pub to_bp: Option<u64>,

    /// Exclude X/Y/MT/unplaced data
    #[arg(long)]
    pub no_xdata: bool,

    /// Number of autosomes (default 22)
    #[arg(long, default_value_t = 22)]
    pub numchrom: u32,

    /// Skip the .geno header hashcheck
    #[arg(long)]
    pub no_hashcheck: bool,

    /// Treat PLINK .fam FID column as pop label
    #[arg(long)]
    pub no_familynames: bool,
}

pub fn run_export(args: ExportArgs) -> Result<()> {
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
    let bad_snps: Option<ahash::AHashSet<String>> = None;
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
        bad: bad_snps.as_ref(),
        snp_keep: None,
        chrom: chrom_filter,
        lopos: args.from_bp,
        hipos: args.to_bp,
        noxdata: args.no_xdata,
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

    // Open reader
    let mut reader =
        crate::pipeline::open_reader_pub(in_fmt, &geno_in, ind_rows.len(), snp_rows.len())?;

    // Read genotypes — we need the full SNP-major matrix in memory for VCF writing.
    // For SampleMajor (TGENO) input, we must transpose to SNP-major first.
    log::info!("reading genotypes...");
    let total_inds = ind_rows.len();
    let total_snps = snp_rows.len();
    let in_rec_bytes = reader.record_bytes();
    let mut in_buf = vec![0u8; in_rec_bytes];

    let genotypes: Vec<Vec<u8>>;

    match reader.layout() {
        crate::geno::Layout::SnpMajor => {
            let mut unpacked = vec![0u8; total_inds];
            let mut genos = Vec::with_capacity(kept_snp_count);
            let mut snp_idx = 0usize;
            while reader.read_record(&mut in_buf)? {
                if snp_idx >= total_snps {
                    break;
                }
                if keep_snps[snp_idx] {
                    codec::unpack(&in_buf, total_inds, &mut unpacked);
                    let row: Vec<u8> = unpacked
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| keep_inds[*i])
                        .map(|(_, &g)| g)
                        .collect();
                    genos.push(row);
                }
                snp_idx += 1;
            }
            while genos.len() < kept_snp_count {
                genos.push(vec![codec::G_MISSING; kept_ind_count]);
            }
            genotypes = genos;
        }
        crate::geno::Layout::SampleMajor => {
            // Each record is one sample × all SNPs.
            // Build a full matrix then transpose to SNP-major.
            let mut unpacked = vec![0u8; total_snps];
            // Collect per-sample genotype rows (only kept samples)
            let mut sample_rows: Vec<Vec<u8>> = Vec::with_capacity(kept_ind_count);
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
                    sample_rows.push(row);
                }
                ind_idx += 1;
            }

            while sample_rows.len() < kept_ind_count {
                sample_rows.push(vec![codec::G_MISSING; kept_snp_count]);
            }

            // Transpose: sample_rows[sample][snp] → genotypes[snp][sample]
            let mut genos: Vec<Vec<u8>> = Vec::with_capacity(kept_snp_count);
            for s in 0..kept_snp_count {
                let mut snp_row = Vec::with_capacity(kept_ind_count);
                for sample in &sample_rows {
                    snp_row.push(sample[s]);
                }
                genos.push(snp_row);
            }
            genotypes = genos;
        }
    }

    // Collect filtered metadata
    let out_snps: Vec<SnpRow> = snp_rows
        .into_iter()
        .zip(keep_snps.iter())
        .filter_map(|(s, &k)| if k { Some(s) } else { None })
        .collect();
    let out_inds: Vec<IndRow> = ind_rows
        .into_iter()
        .zip(keep_inds.iter())
        .filter_map(|(i, &k)| if k { Some(i) } else { None })
        .collect();

    // Write VCF
    let prefix = if args.chr_prefix { "chr" } else { "" };
    log::info!("writing VCF to {}...", args.out.display());
    vcf::write_vcf(
        &args.out, &out_snps, &out_inds, &genotypes, prefix, numchrom,
    )?;

    log::info!(
        "exported {} SNPs × {} samples to VCF in {:.2?}",
        out_snps.len(),
        out_inds.len(),
        t0.elapsed()
    );
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
