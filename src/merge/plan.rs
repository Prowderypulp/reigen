//! Build the merge plan: union (default) or intersection SNP selection across
//! input datasets, with per-dataset flip decisions and a per-dataset
//! position→index map for O(1) lookup.
//!
//! Semantics: default output SNP set is the union of all input SNPs (the
//! "bigger panel"). In union mode, samples from datasets that lack a given SNP
//! get missing calls at that position. In intersection mode, SNPs not present
//! in every dataset are dropped.

use crate::format::{self, Format};
use crate::merge::key::{reconcile, FlipDecision, SnpKey};
use crate::meta::{IndRow, SnpRow};
use anyhow::Result;
use std::collections::HashMap;
use std::path::PathBuf;

pub struct DatasetMetadata {
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
}

pub struct MissnpRecord {
    pub rsid: String,
    pub chrom: u8,
    pub pos: u64,
    pub ref_a1: u8,
    pub ref_a2: u8,
    pub src_a1: Option<u8>,
    pub src_a2: Option<u8>,
    pub reason: &'static str,
}

pub fn build_plan(
    inputs: Vec<(PathBuf, PathBuf, PathBuf)>,
    allow_ambiguous: bool,
    intersection: bool,
    flipstrand: bool,
    numchrom: u32,
    familynames: bool,
    strict_ids: bool,
) -> Result<MergePlan> {
    // --- 1. Load metadata for every dataset and build per-dataset index. ---
    let mut datasets = Vec::with_capacity(inputs.len());
    for (geno, snp, ind) in inputs {
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
            index.insert(
                SnpKey {
                    chrom: s.chrom,
                    pos: s.physical_pos,
                },
                i,
            );
        }
        datasets.push(DatasetMetadata {
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

    for key in key_order {
        let ri = &ref_info[seen[&key]];
        let mut decisions: Vec<Option<(usize, FlipDecision)>> = Vec::with_capacity(datasets.len());
        let mut unresolvable = false;
        let mut ambiguous = false;
        let mut missing_for_intersection = false;

        for ds in &datasets {
            match ds.index.get(&key) {
                Some(&local_idx) => {
                    let s = &ds.snps[local_idx];
                    match reconcile(
                        s.allele1,
                        s.allele2,
                        ri.a1,
                        ri.a2,
                        flipstrand,
                        allow_ambiguous,
                    ) {
                        Some(d) => decisions.push(Some((local_idx, d))),
                        None => {
                            if crate::strand::is_ambiguous(ri.a1, ri.a2) && !allow_ambiguous {
                                ambiguous = true;
                                dropped_snps.push(MissnpRecord {
                                    rsid: ri.rep_id.clone(),
                                    chrom: key.chrom,
                                    pos: key.pos,
                                    ref_a1: ri.a1,
                                    ref_a2: ri.a2,
                                    src_a1: Some(s.allele1),
                                    src_a2: Some(s.allele2),
                                    reason: "ambiguous_at_cg",
                                });
                            } else {
                                unresolvable = true;
                                dropped_snps.push(MissnpRecord {
                                    rsid: ri.rep_id.clone(),
                                    chrom: key.chrom,
                                    pos: key.pos,
                                    ref_a1: ri.a1,
                                    ref_a2: ri.a2,
                                    src_a1: Some(s.allele1),
                                    src_a2: Some(s.allele2),
                                    reason: "unresolvable_alleles",
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
        "merge plan: {} SNPs retained ({} ambiguous, {} unresolvable, {} missing for {})",
        snp_plans.len(),
        dropped_ambiguous,
        dropped_unresolvable,
        dropped_missing_for_intersection,
        if intersection {
            "intersection"
        } else {
            "union"
        }
    );

    // --- 4. Concatenate samples, with duplicate-ID handling. ---
    let output_inds = concat_inds(&datasets, strict_ids)?;

    Ok(MergePlan {
        datasets,
        snp_plans,
        output_inds,
        dropped_snps,
    })
}

fn concat_inds(datasets: &[DatasetMetadata], strict_ids: bool) -> Result<Vec<IndRow>> {
    let total: usize = datasets.iter().map(|d| d.inds.len()).sum();
    let mut out = Vec::with_capacity(total);
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut collisions = 0usize;

    for (ds_idx, ds) in datasets.iter().enumerate() {
        for ind in &ds.inds {
            if seen.contains_key(&ind.id) {
                if strict_ids {
                    anyhow::bail!(
                        "duplicate sample ID {:?} across datasets (use without --strict-ids to auto-rename)",
                        ind.id
                    );
                }
                collisions += 1;
                let renamed = unique_renamed_id(&ind.id, ds_idx, &seen);
                log::warn!(
                    "duplicate sample id {:?} → renamed to {:?}",
                    ind.id,
                    renamed
                );
                let mut clone = ind.clone();
                clone.id = renamed.clone();
                seen.insert(renamed, out.len());
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
    Ok(out)
}

fn unique_renamed_id(base_id: &str, ds_idx: usize, seen: &HashMap<String, usize>) -> String {
    let mut candidate = format!("{base_id}.d{ds_idx}");
    let mut k = 1usize;
    while seen.contains_key(&candidate) {
        candidate = format!("{base_id}.d{ds_idx}.{k}");
        k += 1;
    }
    candidate
}

#[cfg(test)]
mod tests {
    use super::unique_renamed_id;
    use std::collections::HashMap;

    #[test]
    fn renaming_avoids_secondary_collisions() {
        let mut seen = HashMap::new();
        seen.insert("id.d1".to_string(), 0usize);
        seen.insert("id.d1.1".to_string(), 1usize);
        let id = unique_renamed_id("id", 1, &seen);
        assert_eq!(id, "id.d1.2");
    }
}
