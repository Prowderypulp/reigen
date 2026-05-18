//! SNP and sample filters applied during the conversion pipeline.
//!
//! Each filter produces a `Vec<bool>` keep-mask that the pipeline applies
//! while streaming records. No filter ever materializes the full genotype
//! matrix.
//!
//! # SNP filters (convertf parfile keys)
//!
//! - `badsnpname`: newline-separated SNP IDs to drop.
//! - `chrom`: restrict to single chromosome number.
//! - `lopos` / `hipos`: physical position range (inclusive, 1-based).
//! - `noxdata`: drop X/Y/MT/XY SNPs.
//! - `maxmissfracsnp`: requires a missingness pass — two-pass pipeline.
//!
//! # Sample filters
//!
//! - `poplistname`: newline-separated populations to KEEP.
//! - `.ind` column 3 == "Ignore" (case-insensitive): always drop.
//! - `maxmissfracind`: two-pass.

use crate::meta::{IndRow, SnpRow};
use ahash::AHashSet;
use anyhow::{Context, Result};
use std::path::Path;

/// Load `badsnpname` into a hashset of SNP IDs.
pub fn load_bad_snps(path: &Path) -> Result<AHashSet<String>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("read badsnpname {}", path.display()))?;
    Ok(text.split_ascii_whitespace().map(str::to_owned).collect())
}

/// Load `snps` keep list into a hashset of SNP IDs.
pub fn load_snp_keep(path: &Path) -> Result<AHashSet<String>> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("read snps {}", path.display()))?;
    Ok(text.split_ascii_whitespace().map(str::to_owned).collect())
}

/// Load `poplistname`.
pub fn load_pop_keep(path: &Path) -> Result<AHashSet<String>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("read poplistname {}", path.display()))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(str::to_owned)
        .collect())
}

pub type SampleKey = (Option<String>, String);

fn parse_sample_line(line: &str) -> Option<SampleKey> {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return None;
    }
    let mut parts = trimmed.split_whitespace();
    let p1 = parts.next()?;
    if let Some(p2) = parts.next() {
        Some((Some(p1.to_string()), p2.to_string()))
    } else {
        Some((None, p1.to_string()))
    }
}

fn load_sample_list(path: &Path, name: &str) -> Result<AHashSet<SampleKey>> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("read {name} {}", path.display()))?;
    Ok(text.lines().filter_map(parse_sample_line).collect())
}

/// Load `--keep` list into a hashset of SampleKeys.
pub fn load_sample_keep(path: &Path) -> Result<AHashSet<SampleKey>> {
    load_sample_list(path, "keep")
}

/// Load `--remove` list into a hashset of SampleKeys.
pub fn load_sample_remove(path: &Path) -> Result<AHashSet<SampleKey>> {
    load_sample_list(path, "remove")
}

#[derive(Debug, Clone)]
pub enum ChromFilter {
    Single(u8),
    Set(Vec<u8>),
}

impl ChromFilter {
    pub fn parse(s: &str) -> Result<Self> {
        if s.contains(',') || s.contains('-') {
            let mut set = Vec::new();
            for part in s.split(',') {
                if let Some((start, end)) = part.split_once('-') {
                    let start: u8 = start.parse()?;
                    let end: u8 = end.parse()?;
                    for c in start..=end {
                        set.push(c);
                    }
                } else {
                    set.push(part.parse()?);
                }
            }
            set.sort_unstable();
            set.dedup();
            Ok(ChromFilter::Set(set))
        } else {
            Ok(ChromFilter::Single(s.parse()?))
        }
    }

    pub fn contains(&self, c: u8) -> bool {
        match self {
            ChromFilter::Single(x) => *x == c,
            ChromFilter::Set(set) => set.binary_search(&c).is_ok(),
        }
    }
}

pub struct SnpFilter<'a> {
    pub bad: Option<&'a AHashSet<String>>,
    pub snp_keep: Option<&'a AHashSet<String>>,
    pub chrom: Option<ChromFilter>,
    pub lopos: Option<u64>,
    pub hipos: Option<u64>,
    pub noxdata: bool,
    pub x_chrom: u8,  // numchrom+1
    pub y_chrom: u8,  // numchrom+2
    pub mt_chrom: u8, // numchrom+3
    pub xy_chrom: u8, // numchrom+4
}

impl<'a> SnpFilter<'a> {
    pub fn keep(&self, row: &SnpRow) -> bool {
        if let Some(keep) = self.snp_keep {
            if !keep.contains(&row.id) {
                return false;
            }
        }
        if let Some(bad) = self.bad {
            if bad.contains(&row.id) {
                return false;
            }
        }
        if let Some(c) = &self.chrom {
            if !c.contains(row.chrom) {
                return false;
            }
        }
        if let Some(lo) = self.lopos {
            if row.physical_pos < lo {
                return false;
            }
        }
        if let Some(hi) = self.hipos {
            if row.physical_pos > hi {
                return false;
            }
        }
        if self.noxdata
            && (row.chrom == self.x_chrom
                || row.chrom == self.y_chrom
                || row.chrom == self.mt_chrom
                || row.chrom == self.xy_chrom)
        {
            return false;
        }
        true
    }
}

pub struct IndFilter<'a> {
    pub pop_keep: Option<&'a AHashSet<String>>,
    pub sample_keep: Option<&'a AHashSet<SampleKey>>,
    pub sample_remove: Option<&'a AHashSet<SampleKey>>,
}

impl<'a> IndFilter<'a> {
    pub fn keep(&self, row: &IndRow) -> bool {
        if row.ignore {
            return false;
        }
        if let Some(keep) = self.pop_keep {
            if !keep.contains(&row.pop) {
                return false;
            }
        }
        if let Some(keep) = self.sample_keep {
            let mut matched = keep.contains(&(None, row.id.clone()))
                || keep.contains(&(Some(row.pop.clone()), row.id.clone()));

            if !matched && row.id.contains(':') {
                if let Some((fid, iid)) = row.id.split_once(':') {
                    if fid == row.pop {
                        matched = keep.contains(&(None, iid.to_string()))
                            || keep.contains(&(Some(fid.to_string()), iid.to_string()));
                    }
                }
            }

            if !matched {
                return false;
            }
        }
        if let Some(rem) = self.sample_remove {
            let mut matched = rem.contains(&(None, row.id.clone()))
                || rem.contains(&(Some(row.pop.clone()), row.id.clone()));

            if !matched && row.id.contains(':') {
                if let Some((fid, iid)) = row.id.split_once(':') {
                    if fid == row.pop {
                        matched = rem.contains(&(None, iid.to_string()))
                            || rem.contains(&(Some(fid.to_string()), iid.to_string()));
                    }
                }
            }

            if matched {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meta::Sex;

    fn s(id: &str, chrom: u8, pos: u64) -> SnpRow {
        SnpRow {
            id: id.into(),
            chrom,
            genetic_pos: 0.0,
            physical_pos: pos,
            allele1: b'A',
            allele2: b'C',
        }
    }
    fn i(id: &str, pop: &str, ignore: bool) -> IndRow {
        IndRow {
            id: id.into(),
            sex: Sex::Unknown,
            pop: pop.into(),
            ignore,
        }
    }

    #[test]
    fn snp_chrom_and_pos() {
        let f = SnpFilter {
            bad: None,
            snp_keep: None,
            chrom: Some(ChromFilter::Single(1)),
            lopos: Some(100),
            hipos: Some(200),
            noxdata: false,
            x_chrom: 23,
            y_chrom: 24,
            mt_chrom: 25,
            xy_chrom: 26,
        };
        assert!(f.keep(&s("rs1", 1, 150)));
        assert!(!f.keep(&s("rs2", 2, 150))); // wrong chrom
        assert!(!f.keep(&s("rs3", 1, 50))); // below lopos
        assert!(!f.keep(&s("rs4", 1, 250))); // above hipos
    }

    #[test]
    fn noxdata_drops_x() {
        let f = SnpFilter {
            bad: None,
            snp_keep: None,
            chrom: None,
            lopos: None,
            hipos: None,
            noxdata: true,
            x_chrom: 23,
            y_chrom: 24,
            mt_chrom: 25,
            xy_chrom: 26,
        };
        assert!(!f.keep(&s("rsX", 23, 1)));
        assert!(!f.keep(&s("rsY", 24, 1)));
        assert!(!f.keep(&s("rsMT", 25, 1)));
        assert!(!f.keep(&s("rsXY", 26, 1)));
        assert!(f.keep(&s("rs1", 1, 1)));
    }

    #[test]
    fn ignore_drops() {
        let f = IndFilter {
            pop_keep: None,
            sample_keep: None,
            sample_remove: None,
        };
        assert!(!f.keep(&i("S1", "Pop", true)));
        assert!(f.keep(&i("S1", "Pop", false)));
    }

    #[test]
    fn pop_keep_list() {
        let keep: AHashSet<_> = ["French".to_string(), "Han".to_string()]
            .into_iter()
            .collect();
        let f = IndFilter {
            pop_keep: Some(&keep),
            sample_keep: None,
            sample_remove: None,
        };
        assert!(f.keep(&i("S1", "French", false)));
        assert!(!f.keep(&i("S2", "Dinka", false)));
    }

    #[test]
    fn chrom_filter_parse() {
        let cf = ChromFilter::parse("1-3,5,7-8").unwrap();
        assert!(cf.contains(1));
        assert!(cf.contains(2));
        assert!(cf.contains(3));
        assert!(!cf.contains(4));
        assert!(cf.contains(5));
        assert!(!cf.contains(6));
        assert!(cf.contains(7));
        assert!(cf.contains(8));
        assert!(!cf.contains(9));
    }

    #[test]
    fn sample_keep_remove() {
        let mut keep = AHashSet::new();
        keep.insert((None, "IID1".to_string()));
        keep.insert((Some("FID2".to_string()), "IID2".to_string()));

        let mut rem = AHashSet::new();
        rem.insert((None, "BAD".to_string()));

        let f = IndFilter {
            pop_keep: None,
            sample_keep: Some(&keep),
            sample_remove: Some(&rem),
        };

        // Match single IID
        assert!(f.keep(&i("IID1", "Pop1", false)));
        // Match FID+IID
        assert!(f.keep(&i("IID2", "FID2", false)));
        // Miss due to wrong FID
        assert!(!f.keep(&i("IID2", "FID3", false)));
        // Miss completely
        assert!(!f.keep(&i("IID3", "Pop3", false)));

        // Would pass keep list if not in rem
        let f2 = IndFilter {
            pop_keep: None,
            sample_keep: None,
            sample_remove: Some(&rem),
        };
        assert!(!f2.keep(&i("BAD", "PopX", false)));
        assert!(f2.keep(&i("GOOD", "PopY", false)));
    }
}
