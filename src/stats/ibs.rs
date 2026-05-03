//! Pairwise IBS (Identity By State) distance matrix.
//!
//! For each pair of samples (i, j), the IBS proportion is:
//!   IBS(i,j) = (count of matching alleles) / (2 * count of non-missing sites)
//!
//! The IBS distance is 1 - IBS(i,j).
//!
//! Memory: O(nind²) for the accumulator arrays. For AADR scale (17.6k samples),
//! this is ~2.5 GB. The `--ibs` flag is optional for this reason.

use crate::geno::codec;

/// Pairwise IBS accumulators.
pub struct IbsAccumulator {
    nind: usize,
    /// Number of alleles shared IBS (0, 1, or 2 per site) accumulated.
    /// Indexed as `ibs_sum[i * nind + j]` for i < j (upper triangle).
    ibs_sum: Vec<u32>,
    /// Number of sites where both samples are non-missing.
    /// Indexed as `n_compared[i * nind + j]` for i < j.
    n_compared: Vec<u32>,
}

impl IbsAccumulator {
    pub fn new(nind: usize) -> Self {
        let n_pairs = nind * (nind - 1) / 2;
        IbsAccumulator {
            nind,
            ibs_sum: vec![0u32; n_pairs],
            n_compared: vec![0u32; n_pairs],
        }
    }

    /// Index into the upper-triangle flat arrays.
    #[inline]
    fn pair_idx(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < j);
        // Upper triangle index: sum of (nind-1) + (nind-2) + ... + (nind-i) + (j-i-1)
        i * self.nind - i * (i + 1) / 2 + (j - i - 1)
    }

    /// Accumulate one SNP's genotypes for all sample pairs.
    ///
    /// `genos` must have length `nind`, values in {0, 1, 2, 9/G_MISSING}.
    pub fn observe_snp(&mut self, genos: &[u8]) {
        debug_assert_eq!(genos.len(), self.nind);
        let n = self.nind;

        for i in 0..n {
            let gi = genos[i];
            if gi == codec::G_MISSING {
                continue;
            }
            for j in (i + 1)..n {
                let gj = genos[j];
                if gj == codec::G_MISSING {
                    continue;
                }
                let idx = self.pair_idx(i, j);
                self.n_compared[idx] += 1;
                // IBS alleles shared: |gi - gj| gives 0→2 IBS, 1→1 IBS, 2→0 IBS
                let diff = (gi as i16 - gj as i16).unsigned_abs() as u32;
                self.ibs_sum[idx] += 2 - diff;
            }
        }
    }

    /// Compute the full IBS proportion matrix (symmetric, nind × nind).
    /// Diagonal is 1.0. Missing-pair entries (no shared sites) are NaN.
    pub fn ibs_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.nind;
        let mut mat = vec![vec![1.0f64; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let idx = self.pair_idx(i, j);
                let nc = self.n_compared[idx];
                let ibs = if nc == 0 {
                    f64::NAN
                } else {
                    self.ibs_sum[idx] as f64 / (2.0 * nc as f64)
                };
                mat[i][j] = ibs;
                mat[j][i] = ibs;
            }
        }
        mat
    }

    /// Compute 1 - IBS distance matrix.
    pub fn distance_matrix(&self) -> Vec<Vec<f64>> {
        let mat = self.ibs_matrix();
        mat.into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|v| if v.is_nan() { f64::NAN } else { 1.0 - v })
                    .collect()
            })
            .collect()
    }
}

/// Write IBS distance matrix as a square TSV.
pub fn write_ibs_matrix(
    path: &std::path::Path,
    ind_rows: &[crate::meta::IndRow],
    matrix: &[Vec<f64>],
    is_distance: bool,
) -> anyhow::Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    // Header row
    write!(w, "{}", if is_distance { "DIST" } else { "IBS" })?;
    for ind in ind_rows {
        write!(w, "\t{}", ind.id)?;
    }
    writeln!(w)?;

    for (i, row) in matrix.iter().enumerate() {
        write!(w, "{}", ind_rows[i].id)?;
        for val in row {
            if val.is_nan() {
                write!(w, "\tNA")?;
            } else {
                write!(w, "\t{:.6}", val)?;
            }
        }
        writeln!(w)?;
    }

    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ibs_identical_samples() {
        let mut acc = IbsAccumulator::new(2);
        // Both samples have identical genotypes
        acc.observe_snp(&[0, 0]);
        acc.observe_snp(&[1, 1]);
        acc.observe_snp(&[2, 2]);

        let mat = acc.ibs_matrix();
        assert!((mat[0][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn ibs_opposite_samples() {
        let mut acc = IbsAccumulator::new(2);
        // Maximally different
        acc.observe_snp(&[0, 2]);
        acc.observe_snp(&[0, 2]);
        acc.observe_snp(&[0, 2]);

        let mat = acc.ibs_matrix();
        assert!((mat[0][1] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn ibs_het_pair() {
        let mut acc = IbsAccumulator::new(2);
        // One SNP: 0 vs 1 → share 1 allele out of 2 → IBS = 0.5
        acc.observe_snp(&[0, 1]);

        let mat = acc.ibs_matrix();
        assert!((mat[0][1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn ibs_with_missing() {
        let mut acc = IbsAccumulator::new(2);
        acc.observe_snp(&[0, codec::G_MISSING]); // skipped
        acc.observe_snp(&[1, 1]); // IBS = 1.0

        let mat = acc.ibs_matrix();
        assert!((mat[0][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn ibs_three_samples() {
        let mut acc = IbsAccumulator::new(3);
        acc.observe_snp(&[0, 0, 2]);
        acc.observe_snp(&[1, 1, 1]);

        let mat = acc.ibs_matrix();
        // (0,1): SNP0 = same(IBS=1), SNP1 = same(IBS=1) → 1.0
        assert!((mat[0][1] - 1.0).abs() < 1e-10);
        // (0,2): SNP0 = 0 vs 2 (IBS=0), SNP1 = 1 vs 1 (IBS=1) → 0.5
        assert!((mat[0][2] - 0.5).abs() < 1e-10);

        let dist = acc.distance_matrix();
        assert!((dist[0][1] - 0.0).abs() < 1e-10);
        assert!((dist[0][2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn ibs_diagonal_is_one() {
        let acc = IbsAccumulator::new(3);
        let mat = acc.ibs_matrix();
        for i in 0..3 {
            assert!((mat[i][i] - 1.0).abs() < 1e-10);
        }
    }
}
