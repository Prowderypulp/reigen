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
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use self::key::{FlipDecision, ReconcileOpts};
use self::plan::{build_plan, MergeInputSpec};
use self::stream::open_seekable;

#[derive(Args, Debug)]
pub struct MergeArgs {
    /// Input dataset, as geno:snp:ind or as a path prefix. Repeat per dataset.
    #[arg(long = "in", value_name = "PREFIX_OR_TRIPLE")]
    pub inputs: Vec<String>,

    /// Text file with one PREFIX or geno:snp:ind spec per line.
    /// Blank lines and lines starting with '#' are ignored.
    #[arg(long = "in-list", value_name = "FILE")]
    pub in_list: Option<PathBuf>,

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

    /// Disable A1/A2 swap-based genotype flipping during reconciliation.
    #[arg(long)]
    pub no_flip_reference: bool,

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

    /// Disable writing .idmap.tsv report for duplicate sample ID renaming.
    #[arg(long)]
    pub no_idmap: bool,
}

pub fn run_merge(args: MergeArgs) -> Result<()> {
    let mut input_specs = args.inputs.clone();
    if let Some(path) = args.in_list.as_ref() {
        input_specs.extend(read_input_list(path)?);
    }
    if input_specs.len() < 2 {
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

    let parsed_inputs: Vec<ParsedInputSpec> = input_specs
        .iter()
        .enumerate()
        .map(|(idx, s)| parse_input_spec(s, idx))
        .collect::<Result<_>>()?;

    let mut converted_tempdirs = Vec::new();
    let mut merge_inputs = Vec::with_capacity(parsed_inputs.len());
    for (idx, parsed) in parsed_inputs.into_iter().enumerate() {
        let geno = parsed.geno;
        let snp = parsed.snp;
        let ind = parsed.ind;
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
            merge_inputs.push(MergeInputSpec {
                label: parsed.label,
                geno: cfg.geno_out,
                snp: cfg.snp_out,
                ind: cfg.ind_out,
            });
            converted_tempdirs.push(td);
        } else {
            merge_inputs.push(MergeInputSpec {
                label: parsed.label,
                geno,
                snp,
                ind,
            });
        }
    }

    log::info!(
        "building merge plan over {} datasets...",
        merge_inputs.len()
    );
    let plan = build_plan(
        merge_inputs,
        ReconcileOpts {
            flip_strand: args.flip_strand,
            allow_ambiguous: args.allow_ambiguous,
            allow_flip_reference: !args.no_flip_reference,
        },
        args.intersection,
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
    if !args.no_idmap && !plan.renamed_samples.is_empty() {
        let idmap = out_geno.with_extension("idmap.tsv");
        write_idmap_report(&plan.renamed_samples, &idmap)?;
        log::info!(
            "sample id rename report: {} record(s) written to {}",
            plan.renamed_samples.len(),
            idmap.display()
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

struct ParsedInputSpec {
    label: String,
    geno: PathBuf,
    snp: PathBuf,
    ind: PathBuf,
}

/// Accept either `geno:snp:ind` (explicit triple) or a path prefix.
fn parse_input_spec(s: &str, idx: usize) -> Result<ParsedInputSpec> {
    let parts: Vec<&str> = s.split(':').collect();
    match parts.len() {
        3 => {
            let geno = PathBuf::from(parts[0]);
            Ok(ParsedInputSpec {
                label: dataset_label_from_path(&geno, idx),
                geno,
                snp: PathBuf::from(parts[1]),
                ind: PathBuf::from(parts[2]),
            })
        }
        1 => {
            let p = PathBuf::from(parts[0]);
            let bed = p.with_extension("bed");
            let (g, s_ext, i_ext) = if bed.exists() {
                (bed, "bim", "fam")
            } else {
                (p.with_extension("geno"), "snp", "ind")
            };
            Ok(ParsedInputSpec {
                label: dataset_label_from_path(&p, idx),
                geno: g,
                snp: p.with_extension(s_ext),
                ind: p.with_extension(i_ext),
            })
        }
        _ => anyhow::bail!(
            "invalid --in spec {:?} (expected PREFIX or geno:snp:ind)",
            s
        ),
    }
}

fn dataset_label_from_path(path: &Path, idx: usize) -> String {
    path.file_stem()
        .or_else(|| path.file_name())
        .map(|s| s.to_string_lossy().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| format!("dataset{}", idx + 1))
}

fn read_input_list(path: &Path) -> Result<Vec<String>> {
    let f = File::open(path)?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for line in r.lines() {
        let raw = line?;
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        out.push(trimmed.to_string());
    }
    Ok(out)
}

fn write_missnp_reports(
    dropped: &[plan::MissnpRecord],
    tsv_path: &PathBuf,
    rsid_path: &PathBuf,
) -> Result<()> {
    let mut tsv = File::create(tsv_path)?;
    writeln!(
        tsv,
        "rsid\tchrom\tpos\tref_a1\tref_a2\tsrc_a1\tsrc_a2\tdataset\treason"
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
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.rsid,
            rec.chrom,
            rec.pos,
            rec.ref_a1 as char,
            rec.ref_a2 as char,
            src_a1,
            src_a2,
            rec.dataset_label,
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

fn write_idmap_report(renamed: &[plan::RenamedSample], path: &Path) -> Result<()> {
    let mut f = File::create(path)?;
    writeln!(f, "dataset\toriginal_id\trenamed_id")?;
    for rec in renamed {
        writeln!(
            f,
            "{}\t{}\t{}",
            rec.dataset_label, rec.original_id, rec.renamed_id
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{parse_input_spec, read_input_list, write_missnp_reports};
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
                dataset_label: "ds1".into(),
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
                dataset_label: "ds2".into(),
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

    #[test]
    fn parse_input_prefix_spec() {
        let parsed = parse_input_spec("cohort/a", 0).unwrap();
        assert!(parsed.geno.ends_with("cohort/a.geno"));
        assert_eq!(parsed.label, "a");
    }

    #[test]
    fn reads_input_list_file() {
        let d = tempfile::tempdir().unwrap();
        let list = d.path().join("merge_inputs.txt");
        std::fs::write(&list, "# c1\ncohort/a\n\ncohort/b\n").unwrap();
        let lines = read_input_list(&list).unwrap();
        assert_eq!(lines, vec!["cohort/a".to_string(), "cohort/b".to_string()]);
    }
}
