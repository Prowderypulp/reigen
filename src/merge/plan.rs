//! Build the merge plan: union (default) or intersection SNP selection across
//! input datasets, with per-dataset flip decisions and a per-dataset
//! position→index map for O(1) lookup.
//!
//! Semantics: default output SNP set is the union of all input SNPs (the
//! "bigger panel"). In union mode, samples from datasets that lack a given SNP
//! get missing calls at that position. In intersection mode, SNPs not present
//! in every dataset are dropped.

use crate::format::{self, Format};
use crate::merge::key::{reconcile, FlipDecision, ReconcileError, ReconcileOpts, SnpKey};
use crate::meta::{IndRow, SnpRow};
use anyhow::Result;
use std::collections::HashMap;
use std::path::PathBuf;

pub struct MergeInputSpec {
    pub label: String,
    pub geno: PathBuf,
    pub snp: PathBuf,
    pub ind: PathBuf,
}

pub struct DatasetMetadata {
    pub label: String,
    pub format: Format,
    pub geno: PathBuf,
    pub snps: Vec<SnpRow>,
    pub inds: Vec<IndRow>,
    /// (chrom, pos) → local SNP index. Built once at plan time.
    pub index: HashMap<SnpKey, usize>,
}

pub struct SnpPlan {
    pub key: SnpKey,
    pub allele1: u8,
    pub allele2: u8,
    /// Representative SNP id + genetic position (taken from the first
    /// dataset that carries this SNP).
    pub rep_id: String,
    pub rep_gpos: f64,
    /// `dataset_decisions[i]`: Some((local_idx, decision)) if dataset `i`
    /// contains this SNP; None if the dataset is missing it (→ pad with
    /// missing at write time).
    pub dataset_decisions: Vec<Option<(usize, FlipDecision)>>,
}

pub struct MergePlan {
    pub datasets: Vec<DatasetMetadata>,
    pub snp_plans: Vec<SnpPlan>,
    pub output_inds: Vec<IndRow>,
    pub dropped_snps: Vec<MissnpRecord>,
    pub renamed_samples: Vec<RenamedSample>,
}

pub struct MissnpRecord {
    pub rsid: String,
    pub chrom: u8,
    pub pos: u64,
    pub ref_a1: u8,
    pub ref_a2: u8,
    pub src_a1: Option<u8>,
    pub src_a2: Option<u8>,
    pub dataset_label: String,
    pub reason: &'static str,
}

pub struct RenamedSample {
    pub dataset_label: String,
    pub original_id: String,
    pub renamed_id: String,
}

pub fn build_plan(
    inputs: Vec<MergeInputSpec>,
    reconcile_opts: ReconcileOpts,
    intersection: bool,
    numchrom: u32,
    familynames: bool,
    strict_ids: bool,
) -> Result<MergePlan> {
    // --- 1. Load metadata for every dataset and build per-dataset index. ---
    // Positions that occur more than once *within* a single dataset. The merge
    // keys SNPs by (chrom, pos) only, so a duplicated position is ambiguous:
    // `index` would resolve genotypes to the last occurrence while the union
    // keeps the first occurrence's alleles/id, silently pairing one variant's
    // alleles with another's (possibly flipped) genotypes. Such positions are
    // dropped from the merge with a `.missnp` record (B-016).
    let mut dup_positions: HashMap<SnpKey, String> = HashMap::new();
    let mut datasets = Vec::with_capacity(inputs.len());
    for input in inputs {
        let geno = input.geno;
        let snp = input.snp;
        let ind = input.ind;
        let format = format::infer_input_format(&geno)?;
        let snps = if snp.extension().and_then(|e| e.to_str()) == Some("bim") {
            crate::meta::bim::read(&snp, numchrom)?
        } else {
            crate::meta::snp::read(&snp, numchrom)?
        };
        let inds = if ind.extension().and_then(|e| e.to_str()) == Some("fam") {
            crate::meta::fam::read(&ind, familynames)?
        } else {
            crate::meta::ind::read(&ind)?
        };
        let mut index = HashMap::with_capacity(snps.len());
        for (i, s) in snps.iter().enumerate() {
            let key = SnpKey {
                chrom: s.chrom,
                pos: s.physical_pos,
            };
            if index.insert(key, i).is_some() {
                dup_positions
                    .entry(key)
                    .or_insert_with(|| input.label.clone());
            }
        }
        datasets.push(DatasetMetadata {
            label: input.label,
            format,
            geno,
            snps,
            inds,
            index,
        });
    }

    // --- 2. Build union of SNP keys. First occurrence wins for the output
    //        allele pair + representative ID. Sort by (chrom, pos) for
    //        sequential output. ---
    let mut seen: HashMap<SnpKey, usize> = HashMap::new();
    let mut key_order: Vec<SnpKey> = Vec::new();
    #[derive(Clone)]
    struct RefInfo {
        a1: u8,
        a2: u8,
        rep_id: String,
        rep_gpos: f64,
    }
    let mut ref_info: Vec<RefInfo> = Vec::new();

    for ds in &datasets {
        for s in &ds.snps {
            let key = SnpKey {
                chrom: s.chrom,
                pos: s.physical_pos,
            };
            if !seen.contains_key(&key) {
                seen.insert(key, ref_info.len());
                key_order.push(key);
                ref_info.push(RefInfo {
                    a1: s.allele1,
                    a2: s.allele2,
                    rep_id: s.id.clone(),
                    rep_gpos: s.genetic_pos,
                });
            }
        }
    }
    key_order.sort();

    // --- 3. For each key, reconcile every dataset that carries it. Drop the
    //        SNP entirely if any dataset's alleles cannot be reconciled. ---
    let mut snp_plans = Vec::with_capacity(key_order.len());
    let mut dropped_snps = Vec::new();
    let mut dropped_ambiguous = 0usize;
    let mut dropped_unresolvable = 0usize;
    let mut dropped_missing_for_intersection = 0usize;
    let mut dropped_duplicate_position = 0usize;
    // SNPs that only reconciled after a strand complement (PLINK's --flip step).
    let mut strand_flipped = 0usize;

    for key in key_order {
        let ri = &ref_info[seen[&key]];

        // Drop positions duplicated within any input dataset (see B-016 above):
        // they cannot be merged unambiguously by (chrom, pos) alone.
        if let Some(label) = dup_positions.get(&key) {
            dropped_snps.push(MissnpRecord {
                rsid: ri.rep_id.clone(),
                chrom: key.chrom,
                pos: key.pos,
                ref_a1: ri.a1,
                ref_a2: ri.a2,
                src_a1: None,
                src_a2: None,
                dataset_label: label.clone(),
                reason: "duplicate_position",
            });
            dropped_duplicate_position += 1;
            continue;
        }
        let mut decisions: Vec<Option<(usize, FlipDecision)>> = Vec::with_capacity(datasets.len());
        let mut unresolvable = false;
        let mut ambiguous = false;
        let mut missing_for_intersection = false;

        for ds in &datasets {
            match ds.index.get(&key) {
                Some(&local_idx) => {
                    let s = &ds.snps[local_idx];
                    match reconcile(s.allele1, s.allele2, ri.a1, ri.a2, reconcile_opts) {
                        Ok(d) => {
                            // If it would NOT have reconciled on the same strand,
                            // the match required a complement → count the flip.
                            if crate::strand::decide_flip(s.allele1, s.allele2, ri.a1, ri.a2, false)
                                .is_none()
                            {
                                strand_flipped += 1;
                            }
                            decisions.push(Some((local_idx, d)));
                        }
                        Err(e) => {
                            if e == ReconcileError::Ambiguous {
                                ambiguous = true;
                                dropped_snps.push(MissnpRecord {
                                    rsid: ri.rep_id.clone(),
                                    chrom: key.chrom,
                                    pos: key.pos,
                                    ref_a1: ri.a1,
                                    ref_a2: ri.a2,
                                    src_a1: Some(s.allele1),
                                    src_a2: Some(s.allele2),
                                    dataset_label: ds.label.clone(),
                                    reason: "ambiguous_at_cg",
                                });
                            } else if e == ReconcileError::InvalidAllele {
                                unresolvable = true;
                                dropped_snps.push(MissnpRecord {
                                    rsid: ri.rep_id.clone(),
                                    chrom: key.chrom,
                                    pos: key.pos,
                                    ref_a1: ri.a1,
                                    ref_a2: ri.a2,
                                    src_a1: Some(s.allele1),
                                    src_a2: Some(s.allele2),
                                    dataset_label: ds.label.clone(),
                                    reason: "invalid_allele_code",
                                });
                            } else {
                                unresolvable = true;
                                // Classify why it didn't reconcile:
                                // - a flip (swap or complement) WOULD have aligned
                                //   the alleles but a policy flag disallowed it
                                //   (e.g. --no-flip-reference) → swap_disallowed.
                                // - with strand flipping on (default) and no flip
                                //   possible → genuine 3+-allele site (PLINK's
                                //   "remove if it still doesn't work after flip").
                                // - with --no-flip-strand, an off-strand match is
                                //   indistinguishable from triallelic → mismatch.
                                let reason = if crate::strand::decide_flip(
                                    s.allele1,
                                    s.allele2,
                                    ri.a1,
                                    ri.a2,
                                    reconcile_opts.flip_strand,
                                )
                                .is_some()
                                {
                                    "allele_swap_disallowed"
                                } else if reconcile_opts.flip_strand {
                                    "multiallelic"
                                } else {
                                    "allele_mismatch_no_flip"
                                };
                                dropped_snps.push(MissnpRecord {
                                    rsid: ri.rep_id.clone(),
                                    chrom: key.chrom,
                                    pos: key.pos,
                                    ref_a1: ri.a1,
                                    ref_a2: ri.a2,
                                    src_a1: Some(s.allele1),
                                    src_a2: Some(s.allele2),
                                    dataset_label: ds.label.clone(),
                                    reason,
                                });
                            }
                            break;
                        }
                    }
                }
                None => {
                    if intersection {
                        missing_for_intersection = true;
                        dropped_snps.push(MissnpRecord {
                            rsid: ri.rep_id.clone(),
                            chrom: key.chrom,
                            pos: key.pos,
                            ref_a1: ri.a1,
                            ref_a2: ri.a2,
                            src_a1: None,
                            src_a2: None,
                            dataset_label: ds.label.clone(),
                            reason: "missing_in_dataset",
                        });
                        break;
                    }
                    decisions.push(None);
                }
            }
        }

        if ambiguous {
            dropped_ambiguous += 1;
            continue;
        }
        if unresolvable {
            dropped_unresolvable += 1;
            continue;
        }
        if missing_for_intersection {
            dropped_missing_for_intersection += 1;
            continue;
        }

        snp_plans.push(SnpPlan {
            key,
            allele1: ri.a1,
            allele2: ri.a2,
            rep_id: ri.rep_id.clone(),
            rep_gpos: ri.rep_gpos,
            dataset_decisions: decisions,
        });
    }

    log::info!(
        "merge plan: {} SNPs retained ({} strand-flipped, {} ambiguous, {} multiallelic/unresolvable, {} duplicate-position, {} missing for {})",
        snp_plans.len(),
        strand_flipped,
        dropped_ambiguous,
        dropped_unresolvable,
        dropped_duplicate_position,
        dropped_missing_for_intersection,
        if intersection {
            "intersection"
        } else {
            "union"
        }
    );

    // --- 4. Concatenate samples, with duplicate-ID handling. ---
    let (output_inds, renamed_samples) = concat_inds(&datasets, strict_ids)?;

    Ok(MergePlan {
        datasets,
        snp_plans,
        output_inds,
        dropped_snps,
        renamed_samples,
    })
}

fn concat_inds(
    datasets: &[DatasetMetadata],
    strict_ids: bool,
) -> Result<(Vec<IndRow>, Vec<RenamedSample>)> {
    let total: usize = datasets.iter().map(|d| d.inds.len()).sum();
    let mut out = Vec::with_capacity(total);
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut collisions = 0usize;
    let mut renames = Vec::new();

    for ds in datasets {
        for ind in &ds.inds {
            if seen.contains_key(&ind.id) {
                if strict_ids {
                    anyhow::bail!(
                        "duplicate sample ID {:?} across datasets (use without --strict-ids to auto-rename)",
                        ind.id
                    );
                }
                collisions += 1;
                let renamed = unique_renamed_id(&ind.id, &ds.label, &seen);
                log::warn!(
                    "duplicate sample id {:?} → renamed to {:?}",
                    ind.id,
                    renamed
                );
                let mut clone = ind.clone();
                clone.id = renamed.clone();
                seen.insert(renamed, out.len());
                renames.push(RenamedSample {
                    dataset_label: ds.label.clone(),
                    original_id: ind.id.clone(),
                    renamed_id: clone.id.clone(),
                });
                out.push(clone);
            } else {
                seen.insert(ind.id.clone(), out.len());
                out.push(ind.clone());
            }
        }
    }
    if collisions > 0 {
        log::warn!("{} sample id collision(s) auto-renamed", collisions);
    }
    Ok((out, renames))
}

fn unique_renamed_id(base_id: &str, dataset_label: &str, seen: &HashMap<String, usize>) -> String {
    let mut candidate = format!("{base_id}.{}", sanitize_label(dataset_label));
    let mut k = 1usize;
    while seen.contains_key(&candidate) {
        candidate = format!("{base_id}.{}.{}", sanitize_label(dataset_label), k);
        k += 1;
    }
    candidate
}

fn sanitize_label(label: &str) -> String {
    let mut out = String::with_capacity(label.len());
    for c in label.chars() {
        if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
            out.push(c);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        "dataset".to_string()
    } else {
        out
    }
}

#[cfg(test)]
mod tests {
    use super::{build_plan, unique_renamed_id, MergeInputSpec};
    use crate::merge::key::ReconcileOpts;
    use std::collections::HashMap;

    #[test]
    fn renaming_avoids_secondary_collisions() {
        let mut seen = HashMap::new();
        seen.insert("id.ds-2".to_string(), 0usize);
        seen.insert("id.ds-2.1".to_string(), 1usize);
        let id = unique_renamed_id("id", "ds-2", &seen);
        assert_eq!(id, "id.ds-2.2");
    }

    /// Write a minimal PLINK fileset (`.bed`/`.bim`/`.fam`). Only the `.bim`/
    /// `.fam` are parsed at plan time, so the `.bed` just needs the magic.
    fn write_plink(dir: &std::path::Path, name: &str, bim: &str, nfam: usize) -> MergeInputSpec {
        let bed = dir.join(format!("{name}.bed"));
        let bimp = dir.join(format!("{name}.bim"));
        let famp = dir.join(format!("{name}.fam"));
        std::fs::write(&bed, [0x6c, 0x1b, 0x01]).unwrap();
        std::fs::write(&bimp, bim).unwrap();
        let fam: String = (0..nfam)
            .map(|i| format!("{name} {name}_s{i} 0 0 0 -9\n"))
            .collect();
        std::fs::write(&famp, fam).unwrap();
        MergeInputSpec {
            label: name.to_string(),
            geno: bed,
            snp: bimp,
            ind: famp,
        }
    }

    // B-016: two SNPs at the same (chrom, pos) in one input must NOT silently
    // corrupt the merge (the union kept the first variant's alleles while the
    // index resolved genotypes to the last). The position is dropped instead,
    // with a `duplicate_position` missnp record; unambiguous SNPs survive.
    #[test]
    fn duplicate_position_within_dataset_is_dropped_not_corrupted() {
        let dir = tempfile::tempdir().unwrap();
        // ds A carries two variants at 1:5000 (A/G and C/T) plus a clean 1:1000.
        let a = write_plink(
            dir.path(),
            "A",
            "1\trsX\t0\t1000\tA\tG\n1\trsDUP_first\t0\t5000\tA\tG\n1\trsDUP_second\t0\t5000\tC\tT\n",
            2,
        );
        let b = write_plink(dir.path(), "B", "1\trsX\t0\t1000\tA\tG\n", 2);

        let opts = ReconcileOpts {
            flip_strand: true,
            allow_ambiguous: true,
            allow_flip_reference: true,
        };
        let plan = build_plan(vec![a, b], opts, false, 22, true, false).unwrap();

        // Only the unambiguous 1:1000 survives.
        assert_eq!(plan.snp_plans.len(), 1);
        assert_eq!(plan.snp_plans[0].key.pos, 1000);
        // The duplicated position is reported, not merged.
        assert!(plan
            .dropped_snps
            .iter()
            .any(|m| m.pos == 5000 && m.reason == "duplicate_position"));
    }
}
