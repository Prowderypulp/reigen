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
use clap::{Parser, Subcommand};
use reigen::{cmd_convert, cmd_export, cmd_filter, cmd_vcfimport, import, merge, stats};

#[derive(Parser, Debug)]
#[command(name = "reigen", version, about = "Population genomics toolkit")]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Verbose logging
    #[arg(short, long, global = true)]
    verbose: bool,
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

    let level = if cli.verbose { "debug" } else { "info" };
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
