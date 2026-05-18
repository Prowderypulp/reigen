//! `reigen vcfimport` — import biallelic SNPs from VCF into internal formats.
//!
//! Reads a VCF file and writes to any supported output format. Multi-allelic
//! records and indels are skipped. Optional site filtering via a reference
//! `.snp` file.

use crate::format::Format;
use crate::geno::{codec, GenoWriter};
use crate::meta::{self, IndRow, Sex, SnpRow};
use crate::vcf;
use anyhow::{bail, Context, Result};
use clap::Args;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Args, Debug)]
pub struct VcfImportArgs {
    /// Input VCF file
    #[arg(long = "in", value_name = "VCF")]
    pub input: PathBuf,

    /// Reference .snp or .bim file for site filtering (optional).
    /// If provided, only SNPs at positions in this file are kept.
    #[arg(long)]
    pub ref_snp: Option<PathBuf>,

    /// Output format
    #[arg(long)]
    pub out_format: Format,

    /// Output prefix (derives .geno/.snp/.ind or .bed/.bim/.fam)
    #[arg(short = 'o', long)]
    pub out_prefix: Option<String>,

    #[arg(long)]
    pub out_geno: Option<PathBuf>,
    #[arg(long)]
    pub out_snp: Option<PathBuf>,
    #[arg(long)]
    pub out_ind: Option<PathBuf>,

    /// Number of autosomes (default 22)
    #[arg(long, default_value_t = 22)]
    pub numchrom: u32,

    /// Default population label for all samples (default: "Unknown")
    #[arg(long, default_value = "Unknown")]
    pub default_pop: String,

    /// SNP ID list (one rsID per line). Filter output to these SNPs only.
    #[arg(long)]
    pub snplist: Option<PathBuf>,
}

pub fn run_vcfimport(args: VcfImportArgs) -> Result<()> {
    let t0 = Instant::now();

    if args.numchrom > 251 {
        bail!("numchrom {} is too large (max 251)", args.numchrom);
    }

    let p_path = args.out_prefix.as_ref().map(PathBuf::from);
    let (geno_ext, snp_ext, ind_ext) = args.out_format.default_output_extensions();
    let out_geno = args
        .out_geno
        .or_else(|| p_path.as_ref().map(|p| p.with_extension(geno_ext)))
        .ok_or_else(|| anyhow::anyhow!("missing --out-geno or --out-prefix"))?;
    let out_snp = args
        .out_snp
        .or_else(|| p_path.as_ref().map(|p| p.with_extension(snp_ext)))
        .ok_or_else(|| anyhow::anyhow!("missing --out-snp or --out-prefix"))?;
    let out_ind = args
        .out_ind
        .or_else(|| p_path.as_ref().map(|p| p.with_extension(ind_ext)))
        .ok_or_else(|| anyhow::anyhow!("missing --out-ind or --out-prefix"))?;

    // Optional snplist filter
    let snplist_ids: Option<HashSet<String>> = if let Some(ref p) = args.snplist {
        log::info!("loading snplist from {}...", p.display());
        let ids = load_snplist(p)?;
        log::info!("loaded {} SNP IDs", ids.len());
        Some(ids)
    } else {
        None
    };

    // Optional reference filter
    let ref_filter = if let Some(ref_path) = &args.ref_snp {
        log::info!(
            "loading reference SNP filter from {}...",
            ref_path.display()
        );
        let is_bim = ref_path.extension().and_then(|e| e.to_str()) == Some("bim");
        let ref_snps = if is_bim {
            meta::bim::read(ref_path, args.numchrom)?
        } else {
            meta::snp::read(ref_path, args.numchrom)?
        };
        let mut map = HashMap::with_capacity(ref_snps.len());
        for (i, s) in ref_snps.iter().enumerate() {
            map.insert((s.chrom, s.physical_pos), i);
        }
        log::info!("loaded {} reference positions for filtering", map.len());
        Some(map)
    } else {
        None
    };

    // Read VCF
    log::info!("reading VCF from {}...", args.input.display());
    let (sample_names, records, stats) = vcf::read_vcf(
        &args.input,
        args.numchrom,
        ref_filter.as_ref(),
        snplist_ids.as_ref(),
    )?;

    log::info!("VCF read stats:");
    log::info!("  Total records:           {}", stats.total_records);
    log::info!("  Kept (biallelic SNP):    {}", stats.kept_biallelic_snp);
    log::info!("  Skipped (multi-allelic): {}", stats.skipped_multiallelic);
    log::info!("  Skipped (indel):         {}", stats.skipped_indel);
    log::info!("  Skipped (no GT):         {}", stats.skipped_no_gt);
    log::info!("  Skipped (no alt):        {}", stats.skipped_no_alt);
    if stats.skipped_ref_filtered > 0 {
        log::info!("  Skipped (ref filter):    {}", stats.skipped_ref_filtered);
    }
    if stats.skipped_snplist > 0 {
        log::info!("  Skipped (snplist):       {}", stats.skipped_snplist);
    }

    if records.is_empty() {
        bail!("no biallelic SNP records found in VCF");
    }
    if sample_names.is_empty() {
        bail!("no samples found in VCF header");
    }

    let nsnp = records.len();
    let nind = sample_names.len();
    log::info!("importing {} SNPs × {} samples", nsnp, nind);

    // Build SNP and sample metadata
    let out_snps: Vec<SnpRow> = records
        .iter()
        .map(|r| SnpRow {
            id: r.id.clone(),
            chrom: r.chrom,
            genetic_pos: 0.0,
            physical_pos: r.pos,
            allele1: r.ref_allele,
            allele2: r.alt_allele,
        })
        .collect();

    let out_inds: Vec<IndRow> = sample_names
        .iter()
        .map(|name| IndRow {
            id: name.clone(),
            sex: Sex::Unknown,
            pop: args.default_pop.clone(),
            ignore: false,
        })
        .collect();

    // Compute hashes
    let ind_ids: Vec<&str> = out_inds.iter().map(|i| i.id.as_str()).collect();
    let snp_ids: Vec<&str> = out_snps.iter().map(|s| s.id.as_str()).collect();
    let ihash = crate::hash::hasharr(&ind_ids);
    let shash = crate::hash::hasharr(&snp_ids);

    // Open writer
    log::info!(
        "writing to {} ({:?})...",
        out_geno.display(),
        args.out_format
    );
    let mut writer: Box<dyn GenoWriter> = match args.out_format {
        Format::PackedAncestrymap => {
            Box::new(crate::geno::packed_am::PackedAmWriter::create(&out_geno)?)
        }
        Format::Eigenstrat => Box::new(crate::geno::eigenstrat::EigenstratWriter::create(
            &out_geno,
        )?),
        Format::PackedPed => Box::new(crate::geno::packed_ped::PackedPedWriter::create(&out_geno)?),
        Format::Tgeno => Box::new(crate::geno::tgeno::TgenoWriter::create(&out_geno)?),
        _ => bail!("output format not supported for VCF import"),
    };

    writer.begin(nind, nsnp, ihash, shash)?;

    if matches!(args.out_format, Format::Tgeno) {
        // TGENO is sample-major: need to transpose the data
        // Build SNP-major matrix, then transpose
        let snp_rec_bytes = (nind * 2 + 7) / 8;
        let mut snp_major = vec![0u8; nsnp * snp_rec_bytes];
        for (snp_idx, record) in records.iter().enumerate() {
            let row = &mut snp_major[snp_idx * snp_rec_bytes..(snp_idx + 1) * snp_rec_bytes];
            codec::pack(&record.genotypes, row);
        }
        // Transpose to sample-major
        let sample_rec_bytes = (nsnp * 2 + 7) / 8;
        let mut sample_major = vec![0u8; nind * sample_rec_bytes];
        crate::transpose::transpose_packed(&snp_major, nsnp, nind, &mut sample_major)?;
        for i in 0..nind {
            let row = &sample_major[i * sample_rec_bytes..(i + 1) * sample_rec_bytes];
            writer.write_record(row)?;
        }
    } else {
        // SNP-major output
        let rec_bytes = (nind * 2 + 7) / 8;
        let mut out_buf = vec![0u8; rec_bytes];
        for record in &records {
            for b in out_buf.iter_mut() {
                *b = 0;
            }
            codec::pack(&record.genotypes, &mut out_buf);
            writer.write_record(&out_buf)?;
        }
    }

    writer.finish()?;

    // Write metadata
    match args.out_format {
        Format::PackedPed => {
            meta::bim::write(&out_snp, &out_snps, args.numchrom)?;
            meta::fam::write(&out_ind, &out_inds, false)?;
        }
        _ => {
            meta::snp::write(&out_snp, &out_snps, args.numchrom)?;
            meta::ind::write(&out_ind, &out_inds)?;
        }
    }

    log::info!(
        "VCF import done: {} SNPs × {} samples in {:.2?}",
        nsnp,
        nind,
        t0.elapsed()
    );
    Ok(())
}

fn load_snplist(path: &Path) -> Result<HashSet<String>> {
    let f = File::open(path).with_context(|| format!("open snplist {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut ids = HashSet::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("read {}", path.display()))?;
        let id = line.trim();
        if !id.is_empty() && !id.starts_with('#') {
            ids.insert(id.to_owned());
        }
    }
    Ok(ids)
}
