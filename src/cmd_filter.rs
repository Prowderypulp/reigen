use crate::format::Format;
use crate::pipeline::{self, ConvertConfig};
use anyhow::Result;
use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct FilterArgs {
    /// Input prefix (derives .geno/.snp/.ind or .bed/.bim/.fam)
    #[arg(short = 'i', long)]
    pub in_prefix: Option<String>,

    /// Input genotype file
    #[arg(long = "in-geno")]
    pub geno: Option<PathBuf>,

    /// Input SNP file
    #[arg(long = "in-snp")]
    pub snp: Option<PathBuf>,

    /// Input individual/family file
    #[arg(long = "in-ind")]
    pub ind: Option<PathBuf>,

    /// Output format
    #[arg(long)]
    pub out_format: Option<Format>,

    /// Output prefix (derives output paths)
    #[arg(short = 'o', long)]
    pub out_prefix: Option<String>,

    /// Output genotype file
    #[arg(long)]
    pub out_geno: Option<PathBuf>,

    /// Output SNP file
    #[arg(long)]
    pub out_snp: Option<PathBuf>,

    /// Output individual/family file
    #[arg(long)]
    pub out_ind: Option<PathBuf>,

    /// Population list to keep
    #[arg(long)]
    pub poplist: Option<PathBuf>,

    /// Keep only samples in this list (by IID or FID IID)
    #[arg(long)]
    pub keep: Option<PathBuf>,

    /// Remove samples in this list (by IID or FID IID)
    #[arg(long)]
    pub remove: Option<PathBuf>,

    /// List of SNPs to exclude
    #[arg(long, alias = "exclude")]
    pub badsnp: Option<PathBuf>,

    /// List of SNPs to keep
    #[arg(long, alias = "extract", alias = "snplist")]
    pub snps: Option<PathBuf>,

    /// Restrict to a chromosome or range (e.g. 1, 1-5, 1,3-5)
    #[arg(long)]
    pub chrom: Option<String>,

    /// Low BP position (for range filter)
    #[arg(long)]
    pub from_bp: Option<u64>,

    /// High BP position (for range filter)
    #[arg(long)]
    pub to_bp: Option<u64>,

    /// Exclude X/Y/MT/unplaced data
    #[arg(long)]
    pub no_xdata: bool,

    /// Minimum minor allele frequency
    #[arg(long)]
    pub maf: Option<f64>,

    /// Maximum minor allele frequency
    #[arg(long)]
    pub max_maf: Option<f64>,

    /// Hardy-Weinberg equilibrium exact test mid-p threshold (drop SNPs with p < value)
    #[arg(long)]
    pub hwe: Option<f64>,

    /// Maximum per-SNP missingness fraction (PLINK-style --geno semantics).
    #[arg(long = "max-miss-snp", alias = "geno", alias = "maxmissfracsnp")]
    pub max_miss_snp: Option<f64>,

    /// Maximum per-sample missingness fraction (PLINK-style --mind semantics).
    #[arg(long = "mind", alias = "maxmissfracind")]
    pub max_miss_ind: Option<f64>,

    /// Number of autosomes (default 22)
    #[arg(long, default_value_t = 22)]
    pub numchrom: u32,

    /// Skip the .geno header hashcheck (default: enabled)
    #[arg(long)]
    pub no_hashcheck: bool,

    /// Treat PLINK .fam FID column as pop label instead of sample id prefix
    /// (default: use FID as family name, matching upstream `familynames: YES`)
    #[arg(long)]
    pub no_familynames: bool,

    /// Emit population groups in .ind/.fam
    #[arg(long)]
    pub outputgroup: bool,
}

pub fn run_filter(args: FilterArgs) -> Result<()> {
    let (geno_in, snp_in, ind_in) =
        pipeline::resolve_paths(args.in_prefix, args.geno, args.snp, args.ind, None, false)?;

    // Infer format if not provided
    let out_format = match args.out_format {
        Some(f) => f,
        None => crate::format::infer_input_format(&geno_in)?,
    };

    let (geno_out, snp_out, ind_out) = pipeline::resolve_paths(
        args.out_prefix,
        args.out_geno,
        args.out_snp,
        args.out_ind,
        Some(out_format),
        true,
    )?;

    let cfg = ConvertConfig {
        geno_in,
        snp_in,
        ind_in,
        out_fmt: out_format,
        geno_out,
        snp_out,
        ind_out,
        poplist: args.poplist,
        keep: args.keep,
        remove: args.remove,
        badsnp: args.badsnp,
        snps: args.snps,
        chrom: args.chrom,
        lopos: args.from_bp,
        hipos: args.to_bp,
        noxdata: args.no_xdata,
        max_miss_snp: args.max_miss_snp,
        max_miss_ind: args.max_miss_ind,
        maf: args.maf,
        max_maf: args.max_maf,
        hwe: args.hwe,
        numchrom: args.numchrom,
        hashcheck: !args.no_hashcheck,
        familynames: !args.no_familynames,
        outputgroup: args.outputgroup,
    };

    pipeline::run_convert(&cfg)
}
