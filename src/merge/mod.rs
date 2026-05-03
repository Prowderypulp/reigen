//! `reigen merge` — PLINK-style dataset merge.
//!
//! Output SNP set is **union** by default (the "bigger panel"): samples from a
//! dataset that lacks a given SNP get missing calls at that position.
//! Pass `--intersection` to keep only SNPs present in every input dataset.
//!
//! Strand reconciliation via `strand.rs`; ambiguous A/T and C/G SNPs dropped
//! by default (`--allow-ambiguous` to retain).
//!
//! Inputs may be any SNP-major format; EIGENSTRAT text and TGENO sample-major
//! inputs are auto-converted to temporary PAM datasets before merging.

pub mod key;
pub mod plan;
pub mod stream;

use crate::format::Format;
use crate::geno::{codec, GenoWriter};
use crate::meta::{self, SnpRow};
use crate::pipeline;
use anyhow::Result;
use clap::Args;
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use self::key::FlipDecision;
use self::plan::build_plan;
use self::stream::open_seekable;

#[derive(Args, Debug)]
pub struct MergeArgs {
    /// Input dataset, as geno:snp:ind or as a path prefix. Repeat per dataset.
    #[arg(long = "in", value_name = "PREFIX_OR_TRIPLE")]
    pub inputs: Vec<String>,

    /// Output format
    #[arg(long)]
    pub out_format: Format,

    /// Output prefix (derives .geno/.snp/.ind or .bed/.bim/.fam)
    #[arg(short = 'o', long)]
    pub out_prefix: Option<String>,

    /// Output genotype path
    #[arg(long)]
    pub out_geno: Option<PathBuf>,

    /// Output SNP path
    #[arg(long)]
    pub out_snp: Option<PathBuf>,

    /// Output individual/family path
    #[arg(long)]
    pub out_ind: Option<PathBuf>,

    /// Keep A/T and C/G ambiguous SNPs (default: drop).
    #[arg(long, alias = "no-drop-snps")]
    pub allow_ambiguous: bool,

    /// Keep only SNPs present in every dataset (default: union output).
    #[arg(long)]
    pub intersection: bool,

    /// Accept complement matches when reconciling alleles across datasets.
    #[arg(long)]
    pub flip_strand: bool,

    /// Error on duplicate sample IDs instead of auto-renaming id → id.dN.
    #[arg(long)]
    pub strict_ids: bool,

    /// Number of autosomes (default 22)
    #[arg(long, default_value_t = 22)]
    pub numchrom: u32,

    /// Treat FID column in PLINK .fam as part of the population label.
    #[arg(long)]
    pub no_familynames: bool,

    /// Disable writing .missnp and .missnp.tsv reports for dropped SNPs.
    #[arg(long)]
    pub no_missnp: bool,
}

pub fn run_merge(args: MergeArgs) -> Result<()> {
    if args.inputs.len() < 2 {
        anyhow::bail!("merge requires at least two --in datasets");
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

    let parsed_inputs: Vec<(PathBuf, PathBuf, PathBuf)> = args
        .inputs
        .iter()
        .map(|s| parse_input_spec(s))
        .collect::<Result<_>>()?;

    let mut converted_tempdirs = Vec::new();
    let mut merge_inputs = Vec::with_capacity(parsed_inputs.len());
    for (idx, (geno, snp, ind)) in parsed_inputs.into_iter().enumerate() {
        let in_fmt = crate::format::infer_input_format(&geno)?;
        if matches!(in_fmt, Format::Eigenstrat | Format::Tgeno) {
            let td = tempfile::Builder::new()
                .prefix("reigen-merge-autoconv-")
                .tempdir()?;
            let prefix = td.path().join(format!("ds{idx}"));
            let cfg = pipeline::ConvertConfig {
                geno_in: geno,
                snp_in: snp,
                ind_in: ind,
                out_fmt: Format::PackedAncestrymap,
                geno_out: prefix.with_extension("geno"),
                snp_out: prefix.with_extension("snp"),
                ind_out: prefix.with_extension("ind"),
                badsnp: None,
                snps: None,
                poplist: None,
                keep: None,
                remove: None,
                chrom: None,
                lopos: None,
                hipos: None,
                noxdata: false,
                max_miss_snp: None,
                max_miss_ind: None,
                maf: None,
                max_maf: None,
                hwe: None,
                numchrom: args.numchrom,
                hashcheck: true,
                familynames: !args.no_familynames,
                outputgroup: false,
            };
            log::info!(
                "auto-converting merge input {:?} ({in_fmt:?}) to PAM",
                cfg.geno_in
            );
            pipeline::run_convert(&cfg)?;
            merge_inputs.push((cfg.geno_out, cfg.snp_out, cfg.ind_out));
            converted_tempdirs.push(td);
        } else {
            merge_inputs.push((geno, snp, ind));
        }
    }

    log::info!(
        "building merge plan over {} datasets...",
        merge_inputs.len()
    );
    let plan = build_plan(
        merge_inputs,
        args.allow_ambiguous,
        args.intersection,
        args.flip_strand,
        args.numchrom,
        !args.no_familynames,
        args.strict_ids,
    )?;

    log::info!(
        "output: {} SNPs × {} samples",
        plan.snp_plans.len(),
        plan.output_inds.len()
    );
    if !args.no_missnp {
        let missnp_tsv = out_geno.with_extension("missnp.tsv");
        let missnp = out_geno.with_extension("missnp");
        write_missnp_reports(&plan.dropped_snps, &missnp_tsv, &missnp)?;
        log::info!(
            "dropped SNP report: {} record(s) written to {}",
            plan.dropped_snps.len(),
            missnp_tsv.display()
        );
    }

    // Open one seekable reader per dataset; hoist per-dataset read buffers.
    let mut readers = Vec::with_capacity(plan.datasets.len());
    let mut in_bufs: Vec<Vec<u8>> = Vec::with_capacity(plan.datasets.len());
    let mut unpacked_bufs: Vec<Vec<u8>> = Vec::with_capacity(plan.datasets.len());
    for ds in &plan.datasets {
        let r = open_seekable(ds.format, &ds.geno, ds.inds.len(), ds.snps.len())?;
        in_bufs.push(vec![0u8; r.record_bytes()]);
        unpacked_bufs.push(vec![0u8; ds.inds.len()]);
        readers.push(r);
    }

    // Open writer + compute output hashes.
    let out_ind_ids: Vec<&str> = plan.output_inds.iter().map(|i| i.id.as_str()).collect();
    let out_snp_ids: Vec<&str> = plan.snp_plans.iter().map(|s| s.rep_id.as_str()).collect();
    let ihash = crate::hash::hasharr(&out_ind_ids);
    let shash = crate::hash::hasharr(&out_snp_ids);

    let mut writer: Box<dyn GenoWriter> = match args.out_format {
        Format::PackedAncestrymap => {
            Box::new(crate::geno::packed_am::PackedAmWriter::create(&out_geno)?)
        }
        Format::Eigenstrat => Box::new(crate::geno::eigenstrat::EigenstratWriter::create(
            &out_geno,
        )?),
        Format::PackedPed => Box::new(crate::geno::packed_ped::PackedPedWriter::create(&out_geno)?),
        Format::Tgeno => Box::new(crate::geno::tgeno::TgenoWriter::create(&out_geno)?),
        f => anyhow::bail!("output format {f:?} not supported for merge"),
    };
    writer.begin(plan.output_inds.len(), plan.snp_plans.len(), ihash, shash)?;

    // Build merged SNP-major rows first. For TGENO output we transpose and emit sample-major.
    let nsnp_out = plan.snp_plans.len();
    let nind_out = plan.output_inds.len();
    let out_rec_bytes = (nind_out * 2 + 7) / 8;
    let mut out_buf = vec![0u8; out_rec_bytes];
    let mut merged_genos = vec![0u8; plan.output_inds.len()];
    let mut merged_snp_major = if matches!(args.out_format, Format::Tgeno) {
        vec![0u8; nsnp_out * out_rec_bytes]
    } else {
        Vec::new()
    };

    for (snp_idx, snp_plan) in plan.snp_plans.iter().enumerate() {
        let mut offset = 0usize;
        for (ds_idx, decision) in snp_plan.dataset_decisions.iter().enumerate() {
            let nind = plan.datasets[ds_idx].inds.len();
            match decision {
                Some((local_idx, flip)) => {
                    let reader = &mut readers[ds_idx];
                    reader.seek_record(*local_idx)?;
                    let in_buf = &mut in_bufs[ds_idx];
                    if !reader.read_record(in_buf)? {
                        anyhow::bail!(
                            "dataset {} ran out of records at SNP ({},{})",
                            ds_idx,
                            snp_plan.key.chrom,
                            snp_plan.key.pos
                        );
                    }
                    let slice = &mut merged_genos[offset..offset + nind];
                    codec::unpack(in_buf, nind, slice);
                    if *flip == FlipDecision::Flip {
                        for g in slice.iter_mut() {
                            if *g == 0 {
                                *g = 2;
                            } else if *g == 2 {
                                *g = 0;
                            }
                        }
                    }
                }
                None => {
                    // Dataset lacks this SNP → pad its samples with missing.
                    for g in &mut merged_genos[offset..offset + nind] {
                        *g = 3;
                    }
                }
            }
            offset += nind;
        }

        for b in out_buf.iter_mut() {
            *b = 0;
        }
        codec::pack(&merged_genos, &mut out_buf);
        if matches!(args.out_format, Format::Tgeno) {
            let row = &mut merged_snp_major[snp_idx * out_rec_bytes..(snp_idx + 1) * out_rec_bytes];
            row.copy_from_slice(&out_buf);
        } else {
            writer.write_record(&out_buf)?;
        }
    }

    if matches!(args.out_format, Format::Tgeno) {
        let sample_row_bytes = (nsnp_out * 2 + 7) / 8;
        let mut sample_major = vec![0u8; nind_out * sample_row_bytes];
        crate::transpose::transpose_packed(
            &merged_snp_major,
            nsnp_out,
            nind_out,
            &mut sample_major,
        )?;
        for i in 0..nind_out {
            let row = &sample_major[i * sample_row_bytes..(i + 1) * sample_row_bytes];
            writer.write_record(row)?;
        }
    }
    writer.finish()?;

    // Output metadata.
    let out_snp_rows: Vec<SnpRow> = plan
        .snp_plans
        .iter()
        .map(|s| SnpRow {
            id: s.rep_id.clone(),
            chrom: s.key.chrom,
            genetic_pos: s.rep_gpos,
            physical_pos: s.key.pos,
            allele1: s.allele1,
            allele2: s.allele2,
        })
        .collect();

    match args.out_format {
        Format::PackedPed => {
            meta::bim::write(&out_snp, &out_snp_rows, args.numchrom)?;
            meta::fam::write(&out_ind, &plan.output_inds, false)?;
        }
        _ => {
            meta::snp::write(&out_snp, &out_snp_rows, args.numchrom)?;
            meta::ind::write(&out_ind, &plan.output_inds)?;
        }
    }

    log::info!("merge done.");
    drop(converted_tempdirs);
    Ok(())
}

/// Accept either `geno:snp:ind` (explicit triple) or a path prefix.
fn parse_input_spec(s: &str) -> Result<(PathBuf, PathBuf, PathBuf)> {
    let parts: Vec<&str> = s.split(':').collect();
    match parts.len() {
        3 => Ok((
            PathBuf::from(parts[0]),
            PathBuf::from(parts[1]),
            PathBuf::from(parts[2]),
        )),
        1 => {
            let p = PathBuf::from(parts[0]);
            let bed = p.with_extension("bed");
            let (g, s_ext, i_ext) = if bed.exists() {
                (bed, "bim", "fam")
            } else {
                (p.with_extension("geno"), "snp", "ind")
            };
            Ok((g, p.with_extension(s_ext), p.with_extension(i_ext)))
        }
        _ => anyhow::bail!(
            "invalid --in spec {:?} (expected PREFIX or geno:snp:ind)",
            s
        ),
    }
}

fn write_missnp_reports(
    dropped: &[plan::MissnpRecord],
    tsv_path: &PathBuf,
    rsid_path: &PathBuf,
) -> Result<()> {
    let mut tsv = File::create(tsv_path)?;
    writeln!(
        tsv,
        "rsid\tchrom\tpos\tref_a1\tref_a2\tsrc_a1\tsrc_a2\treason"
    )?;
    for rec in dropped {
        let src_a1 = rec
            .src_a1
            .map(|b| (b as char).to_string())
            .unwrap_or_default();
        let src_a2 = rec
            .src_a2
            .map(|b| (b as char).to_string())
            .unwrap_or_default();
        writeln!(
            tsv,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.rsid,
            rec.chrom,
            rec.pos,
            rec.ref_a1 as char,
            rec.ref_a2 as char,
            src_a1,
            src_a2,
            rec.reason
        )?;
    }

    let mut seen = HashSet::new();
    let mut rsid = File::create(rsid_path)?;
    for rec in dropped {
        if seen.insert(rec.rsid.clone()) {
            writeln!(rsid, "{}", rec.rsid)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_missnp_reports;
    use crate::merge::plan::MissnpRecord;

    #[test]
    fn writes_missnp_outputs() {
        let d = tempfile::tempdir().unwrap();
        let tsv = d.path().join("out.missnp.tsv");
        let missnp = d.path().join("out.missnp");
        let dropped = vec![
            MissnpRecord {
                rsid: "rs1".into(),
                chrom: 1,
                pos: 100,
                ref_a1: b'A',
                ref_a2: b'G',
                src_a1: Some(b'C'),
                src_a2: Some(b'T'),
                reason: "unresolvable_alleles",
            },
            MissnpRecord {
                rsid: "rs1".into(),
                chrom: 1,
                pos: 100,
                ref_a1: b'A',
                ref_a2: b'G',
                src_a1: None,
                src_a2: None,
                reason: "missing_in_dataset",
            },
        ];
        write_missnp_reports(&dropped, &tsv, &missnp).unwrap();
        let tsv_text = std::fs::read_to_string(tsv).unwrap();
        let missnp_text = std::fs::read_to_string(missnp).unwrap();
        assert!(tsv_text.contains("unresolvable_alleles"));
        assert_eq!(missnp_text.lines().count(), 1);
        assert_eq!(missnp_text.trim(), "rs1");
    }
}
