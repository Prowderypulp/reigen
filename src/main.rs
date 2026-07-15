//! reigen — Population genomics toolkit.
//!
//! Subcommands:
//! - convert     Convert between genotype formats
//! - import      Import DTC raw files
//! - merge       Merge 2+ datasets
//! - export      Export to VCF
//! - vcfimport   Import from VCF
//! - stats       Compute QC and summary statistics

use anyhow::{Context, Result};
use clap::{ArgAction, Parser, Subcommand};
use reigen::{cmd_convert, cmd_export, cmd_filter, cmd_vcfimport, import, merge, stats};

#[derive(Parser, Debug)]
#[command(
    name = "reigen",
    version,
    about = "Population genomics toolkit — convert, filter, merge, import and export genotype data",
    long_about = "reigen converts, filters, merges, and QCs genotype datasets across the \
                  EIGENSTRAT/AdmixTools, PLINK1, PLINK2, TGENO and VCF formats, and imports \
                  direct-to-consumer DNA kits.\n\n\
                  Most subcommands accept `-i/--in-prefix <prefix>` to derive the input triple \
                  automatically. By default reigen is quiet; pass -v for progress, -vv for debug.",
    propagate_version = true,
    after_help = "EXAMPLES:\n  \
        # Convert a PLINK1 dataset (data.bed/.bim/.fam) to PACKEDANCESTRYMAP\n  \
        reigen convert -i data --out-format packedancestrymap -o data_pam\n\n  \
        # Filter to chr1-5, MAF >= 1%, keeping a sample list\n  \
        reigen filter -i data -o data_qc --chrom 1-5 --maf 0.01 --keep keep.txt\n\n  \
        # Import a 23andMe kit and write PLINK1\n  \
        reigen import --in kit.txt --out-format packedped -o kit\n\n  \
        # Merge two panels, then export to VCF\n  \
        reigen merge --in panelA --in panelB --out-format packedped -o merged\n  \
        reigen export -i merged -o merged.vcf\n\n\
        Run 'reigen <command> --help' for the options and examples of each command."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Increase logging: -v = progress, -vv = debug (default: warnings only)
    #[arg(short, long, global = true, action = ArgAction::Count)]
    verbose: u8,

    /// Silence all output except errors
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Convert between genotype formats
    Convert(cmd_convert::ConvertArgs),
    /// Import DTC raw files (23andMe, Ancestry, etc)
    Import(import::ImportArgs),
    /// Merge 2+ datasets with strand reconciliation
    Merge(merge::MergeArgs),
    /// Export to VCF format
    Export(cmd_export::ExportArgs),
    /// Import biallelic SNPs from VCF
    Vcfimport(cmd_vcfimport::VcfImportArgs),
    /// Filter and subset a dataset
    Filter(cmd_filter::FilterArgs),
    /// Compute QC and summary statistics
    Stats(stats::StatsArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Quiet by default (standard CLI behaviour): only warnings and errors reach
    // the terminal. `-v` promotes to progress (info), `-vv`+ to debug, `-q` to
    // errors only. An explicit RUST_LOG env var still overrides all of this.
    let level = if cli.quiet {
        "error"
    } else {
        match cli.verbose {
            0 => "warn",
            1 => "info",
            _ => "debug",
        }
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level))
        .format_timestamp(None)
        .init();

    match cli.command {
        Commands::Convert(args) => cmd_convert::run(args).context("convert failed")?,
        Commands::Import(args) => import::run_import(args).context("import failed")?,
        Commands::Merge(args) => merge::run_merge(args).context("merge failed")?,
        Commands::Export(args) => cmd_export::run_export(args).context("export failed")?,
        Commands::Vcfimport(args) => {
            cmd_vcfimport::run_vcfimport(args).context("vcfimport failed")?
        }
        Commands::Filter(args) => cmd_filter::run_filter(args).context("filter failed")?,
        Commands::Stats(args) => stats::run_stats(args).context("stats failed")?,
    }

    Ok(())
}
