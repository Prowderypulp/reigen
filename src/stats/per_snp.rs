//! Per-SNP QC statistics: allele frequency, genotype counts, missingness,
//! HWE exact test, and observed heterozygosity.

/// Accumulated counters for a single SNP across all kept samples.
#[derive(Debug, Clone, Default)]
pub struct SnpStats {
    pub n_hom_ref: u32, // genotype = 2 (two copies of allele1 / reference)
    pub n_het: u32,     // genotype = 1
    pub n_hom_alt: u32, // genotype = 0 (zero copies of allele1 / reference)
    pub n_missing: u32, // genotype = 9 / missing
}

impl SnpStats {
    /// Total genotyped (non-missing) samples.
    #[inline]
    pub fn n_called(&self) -> u32 {
        self.n_hom_ref + self.n_het + self.n_hom_alt
    }

    /// Total samples (called + missing).
    #[inline]
    pub fn n_total(&self) -> u32 {
        self.n_called() + self.n_missing
    }

    /// Per-SNP missingness rate.
    pub fn miss_rate(&self) -> f64 {
        let total = self.n_total();
        if total == 0 {
            return 0.0;
        }
        self.n_missing as f64 / total as f64
    }

    /// Reference allele frequency.
    pub fn ref_freq(&self) -> f64 {
        let n = self.n_called();
        if n == 0 {
            return f64::NAN;
        }
        let total_alleles = n as f64 * 2.0;
        (self.n_hom_ref as f64 * 2.0 + self.n_het as f64) / total_alleles
    }

    /// Alternate allele frequency.
    ///
    /// Counted directly from the alternate-allele copies rather than as
    /// `1 - ref_freq`: subtracting from 1.0 introduces a rounding error (e.g.
    /// `1.0 - 0.9 == 0.09999999999999998`) that flips exact-boundary `--maf`
    /// comparisons against plink (B-014). `2*n_hom_alt + n_het` over `2n` is
    /// exact for the common denominators.
    pub fn alt_freq(&self) -> f64 {
        let n = self.n_called();
        if n == 0 {
            return f64::NAN;
        }
        let total_alleles = n as f64 * 2.0;
        (self.n_hom_alt as f64 * 2.0 + self.n_het as f64) / total_alleles
    }

    /// Minor allele frequency.
    ///
    /// `min(ref_freq, alt_freq)` where both are counted directly from integer
    /// allele copies (see `alt_freq`), so a SNP at exactly the threshold
    /// compares bit-identically to plink instead of being spuriously dropped.
    pub fn maf(&self) -> f64 {
        let rf = self.ref_freq();
        if rf.is_nan() {
            return f64::NAN;
        }
        rf.min(self.alt_freq())
    }

    /// Observed heterozygosity.
    pub fn obs_het(&self) -> f64 {
        let n = self.n_called();
        if n == 0 {
            return f64::NAN;
        }
        self.n_het as f64 / n as f64
    }

    /// HWE exact mid-p test p-value (Wigginton, Cutler & Abecasis 2005).
    ///
    /// Returns the two-sided mid-p value. Small values indicate departure
    /// from HWE.
    pub fn hwe_pvalue(&self) -> f64 {
        let n_ab = self.n_het as usize;
        let n_aa = self.n_hom_ref as usize;
        let n_bb = self.n_hom_alt as usize;
        hwe_exact_midp(n_aa, n_ab, n_bb)
    }

    /// Accumulate a single genotype value.
    #[inline]
    pub fn observe(&mut self, g: u8) {
        match g {
            // genotype = count of allele1 (reference): 2 = hom reference,
            // 0 = hom alternate. (See EIGENSTRAT/AdmixTools convention.)
            2 => self.n_hom_ref += 1,
            1 => self.n_het += 1,
            0 => self.n_hom_alt += 1,
            _ => self.n_missing += 1,
        }
    }
}

/// HWE exact mid-p test.
///
/// Implements the algorithm from Wigginton, Cutler & Abecasis (2005)
/// "A Note on Exact Tests of Hardy-Weinberg Equilibrium."
/// Am. J. Hum. Genet. 76:887-893.
///
/// We compute the probability of each possible heterozygote count given
/// the observed allele counts, using the recursion for the probability
/// distribution. The mid-p value is the sum of probabilities of tables
/// less extreme than observed, plus half the probability of the observed
/// table.
fn hwe_exact_midp(n_aa: usize, n_ab: usize, n_bb: usize) -> f64 {
    let n = n_aa + n_ab + n_bb;
    if n == 0 {
        return 1.0;
    }

    let n_a = 2 * n_aa + n_ab; // total count of allele A
    let n_b = 2 * n_bb + n_ab; // total count of allele B

    // NB: a monomorphic locus (n_a == 0 || n_b == 0) is *not* special-cased to
    // 1.0 here. Its distribution is a single point mass at the observed table,
    // for which the mid-p value is 0.5 * P(obs) = 0.5 — matching plink's
    // `--hardy midp`. The general path below yields exactly that; an early
    // `return 1.0` would report the plain exact p instead of the mid-p and
    // diverge from the oracle on every monomorphic SNP (B-015).

    // Maximum possible heterozygotes given allele counts
    let max_het = n_a.min(n_b);
    // Heterozygote count must have same parity as n_a (and n_b)
    // because n_ab = n_a - 2*n_aa, and n_aa >= 0.

    // Build probability table using recurrence.
    // P(n_het) is proportional to choose(n, n_het) * ... but we use
    // the recursion from the original paper for numerical stability.
    //
    // Start from the minimum possible n_het (0 or 1 depending on parity)
    // and compute relative probabilities.
    let start = n_a % 2; // 0 if n_a even, 1 if odd
    let mut probs: Vec<f64> = Vec::with_capacity((max_het - start) / 2 + 1);

    // We compute log-probabilities relative to the probability at start.
    // P(het=k) / P(het=k-2) = (n_aa_k * n_bb_k * 4) / ((k) * (k-1))
    // where n_aa_k = (n_a - k)/2, n_bb_k = (n_b - k)/2.

    // Direct computation using a ratio recurrence:
    // Let het go from start to max_het in steps of 2.
    // prob[i] = relative probability of het = start + 2*i.

    // Accumulate in LOG space, then exp after subtracting the max. The naive
    // `prev_prob * ratio` recurrence overflows f64 around N≈1000 (the peak
    // unnormalized probability relative to het=start exceeds ~1e308), turning
    // the whole table to Inf/NaN and silently breaking `--hwe` on large
    // cohorts. Summing log-ratios and rescaling by the max is overflow-safe for
    // any N. (B-101)
    let mut log_probs: Vec<f64> = Vec::with_capacity((max_het - start) / 2 + 1);
    log_probs.push(0.0); // ln(1) at het=start

    let mut cum_log = 0.0f64;
    let mut het = start + 2;
    while het <= max_het {
        let prev_het = het - 2;
        let aa_at_prev = (n_a - prev_het) / 2;
        let bb_at_prev = (n_b - prev_het) / 2;
        // ratio = P(het) / P(het-2)
        let ratio = (4.0 * aa_at_prev as f64 * bb_at_prev as f64) / (het as f64 * (het - 1) as f64);
        cum_log += ratio.ln(); // ratio == 0 → -inf → prob 0, which is correct
        log_probs.push(cum_log);
        het += 2;
    }

    let max_log = log_probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    for lp in log_probs.iter() {
        probs.push((lp - max_log).exp());
    }

    // Normalize
    let total: f64 = probs.iter().sum();
    if total <= 0.0 {
        return 1.0;
    }
    for p in probs.iter_mut() {
        *p /= total;
    }

    // Find observed index
    if n_ab < start || (n_ab - start) % 2 != 0 {
        // Observed het count not achievable — shouldn't happen with valid data
        return 1.0;
    }
    let obs_idx = (n_ab - start) / 2;
    if obs_idx >= probs.len() {
        return 1.0;
    }

    let obs_prob = probs[obs_idx];

    // Two-sided mid-p: tables strictly *more* extreme than the observed one (i.e.
    // lower probability) count at full weight; tables *as* extreme as the
    // observed (equal probability — the observed table AND any ties) count at
    // HALF weight. A symmetric HWE distribution routinely has a second table
    // with probability identical to the observed (e.g. het=k and het=k+2 both at
    // the mode); counting those ties at full weight — as the old code did, since
    // it added every `p <= obs_prob` at full weight and only halved the single
    // observed index — inflates the p-value by 0.5*P(tie) and diverges from
    // plink's `--hardy midp` (B-015). Compare with a relative tolerance so
    // genuinely-equal tables are recognised despite float rounding.
    let lo = obs_prob * (1.0 - 1e-7);
    let hi = obs_prob * (1.0 + 1e-7);
    let mut p_value = 0.0;
    for &p in probs.iter() {
        if p < lo {
            p_value += p;
        } else if p <= hi {
            p_value += 0.5 * p;
        }
    }

    // Clamp to [0, 1]
    p_value.clamp(0.0, 1.0)
}

/// Write per-SNP stats as TSV.
pub fn write_snp_stats(
    path: &std::path::Path,
    snp_rows: &[crate::meta::SnpRow],
    stats: &[SnpStats],
) -> anyhow::Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    writeln!(
        w,
        "SNP\tCHROM\tPOS\tA1\tA2\tN_CALLED\tN_MISS\tMISS_RATE\tREF_FREQ\tALT_FREQ\tMAF\tOBS_HET\tHWE_P\tN_HOMREF\tN_HET\tN_HOMALT"
    )?;

    for (snp, st) in snp_rows.iter().zip(stats.iter()) {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.4e}\t{}\t{}\t{}",
            snp.id,
            snp.chrom,
            snp.physical_pos,
            snp.allele1 as char,
            snp.allele2 as char,
            st.n_called(),
            st.n_missing,
            st.miss_rate(),
            st.ref_freq(),
            st.alt_freq(),
            st.maf(),
            st.obs_het(),
            st.hwe_pvalue(),
            st.n_hom_ref,
            st.n_het,
            st.n_hom_alt,
        )?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geno::codec;

    #[test]
    fn snp_stats_basic() {
        let mut s = SnpStats::default();
        s.observe(0);
        s.observe(0);
        s.observe(1);
        s.observe(2);
        s.observe(codec::G_MISSING);

        assert_eq!(s.n_called(), 4);
        assert_eq!(s.n_total(), 5);
        assert!((s.miss_rate() - 0.2).abs() < 1e-10);
        // genotypes 0,0,1,2: ref allele (allele1) copies = one hom-ref (g=2)
        // contributes 2, the het contributes 1 → ref freq = (2*1 + 1)/(4*2)
        // = 3/8 = 0.375. (g=0 are hom-alt, contribute 0 reference alleles.)
        assert!((s.ref_freq() - 0.375).abs() < 1e-10);
        assert!((s.alt_freq() - 0.625).abs() < 1e-10);
        assert!((s.maf() - 0.375).abs() < 1e-10);
        assert!((s.obs_het() - 0.25).abs() < 1e-10);
    }

    #[test]
    fn maf_exact_boundary_is_representable() {
        // B-014: MAF must be counted directly, not via `1.0 - ref_freq`, or a
        // SNP whose true MAF is exactly the threshold gets a value of
        // 0.09999999999999998 and is spuriously dropped by `--maf 0.1`.
        // 12 ALT copies out of 120 (48 hom-ref, 0 het? use 54 hom-ref, 12 het,..)
        // Construct 60 samples: alt copies = 12 => 54 hom-ref, 0 het, 6 hom-alt
        // gives alt = 12, ref = 108. maf = 0.1 exactly.
        let s = SnpStats {
            n_hom_ref: 54, // g=2, ref alleles
            n_het: 0,
            n_hom_alt: 6, // g=0, alt alleles => 12 alt copies
            n_missing: 0,
        };
        assert_eq!(s.alt_freq(), 0.1, "alt_freq must be exactly 0.1");
        assert_eq!(s.maf(), 0.1, "maf must be exactly 0.1, not 0.0999…");
        assert!(!(s.maf() < 0.1), "maf must not fall below the 0.1 threshold");
    }

    #[test]
    fn hwe_perfect_equilibrium() {
        // 25 AA, 50 AB, 25 BB → p=0.5, expected HWE
        // Should give p-value close to 1.0
        let pval = hwe_exact_midp(25, 50, 25);
        assert!(pval > 0.9, "perfect HWE p={pval} should be ~1.0");
    }

    #[test]
    fn hwe_excess_het() {
        // All het: 0 AA, 100 AB, 0 BB → extreme excess of hets
        let pval = hwe_exact_midp(0, 100, 0);
        assert!(pval < 0.001, "extreme het excess p={pval} should be tiny");
    }

    #[test]
    fn hwe_no_het() {
        // No het: 50 AA, 0 AB, 50 BB → extreme deficit of hets
        let pval = hwe_exact_midp(50, 0, 50);
        assert!(pval < 0.001, "extreme het deficit p={pval} should be tiny");
    }

    #[test]
    fn hwe_empty() {
        assert_eq!(hwe_exact_midp(0, 0, 0), 1.0);
    }

    #[test]
    fn hwe_large_n_no_overflow() {
        // B-101: the old raw-f64 recurrence overflowed to NaN for large N.
        // Perfect equilibrium at N=50_000 must still yield a finite, sane p.
        let p = hwe_exact_midp(12_500, 25_000, 12_500);
        assert!(p.is_finite() && !p.is_nan(), "p must be finite, got {p}");
        assert!((0.0..=1.0).contains(&p), "p in range, got {p}");
        assert!(p > 0.5, "perfect equilibrium should not be rejected, p={p}");

        // A strong het deficit at large N stays finite and tiny.
        let p2 = hwe_exact_midp(25_000, 1_000, 24_000);
        assert!(p2.is_finite() && !p2.is_nan(), "p2 finite, got {p2}");
        assert!(p2 < 0.001, "large-N het deficit should be significant, p={p2}");
    }

    #[test]
    fn hwe_monomorphic() {
        // B-015: a monomorphic locus has a single-point distribution; its two-
        // sided *mid*-p is 0.5 * P(obs) = 0.5, matching `plink --hardy midp`.
        // (The old code special-cased this to 1.0 — the plain exact p — which
        // silently diverged from the oracle on every monomorphic SNP.)
        assert!((hwe_exact_midp(100, 0, 0) - 0.5).abs() < 1e-9);
        assert!((hwe_exact_midp(0, 0, 100) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn hwe_midp_halves_ties() {
        // B-015: when a second table is exactly as probable as the observed one
        // (the HWE distribution is symmetric, so het=k and het=k+2 tie at the
        // mode), both must be counted at HALF weight. Ground truth from an
        // independent lgamma implementation, matching `plink --hardy midp`.
        //   15 AA / 14 AB / 4 BB  -> mid-p = 0.71458 (obs het=14 ties het=16)
        let p = hwe_exact_midp(15, 14, 4);
        assert!((p - 0.71458).abs() < 1e-4, "midp with tie, got {p}");
        //   26 AA / 28 AB / 6 BB  -> mid-p = 0.77976
        let p2 = hwe_exact_midp(26, 28, 6);
        assert!((p2 - 0.77976).abs() < 1e-4, "midp with tie, got {p2}");
    }
}
