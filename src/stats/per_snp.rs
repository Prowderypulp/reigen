//! Per-SNP QC statistics: allele frequency, genotype counts, missingness,
//! HWE exact test, and observed heterozygosity.

/// Accumulated counters for a single SNP across all kept samples.
#[derive(Debug, Clone, Default)]
pub struct SnpStats {
    pub n_hom_ref: u32, // genotype = 0
    pub n_het: u32,     // genotype = 1
    pub n_hom_alt: u32, // genotype = 2
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

    /// Alternate allele frequency (= 1 - ref_freq).
    pub fn alt_freq(&self) -> f64 {
        let rf = self.ref_freq();
        if rf.is_nan() {
            return f64::NAN;
        }
        1.0 - rf
    }

    /// Minor allele frequency.
    pub fn maf(&self) -> f64 {
        let af = self.alt_freq();
        if af.is_nan() {
            return f64::NAN;
        }
        af.min(1.0 - af)
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
            0 => self.n_hom_ref += 1,
            1 => self.n_het += 1,
            2 => self.n_hom_alt += 1,
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

    if n_a == 0 || n_b == 0 {
        return 1.0;
    }

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

    // Compute all probabilities relative to het=start.
    probs.push(1.0); // prob at het=start is 1.0 (unnormalized)

    let mut het = start + 2;
    while het <= max_het {
        let prev_het = het - 2;
        let aa_at_prev = (n_a - prev_het) / 2;
        let bb_at_prev = (n_b - prev_het) / 2;
        // ratio = P(het) / P(het-2)
        let ratio = (4.0 * aa_at_prev as f64 * bb_at_prev as f64) / (het as f64 * (het - 1) as f64);
        let prev_prob = *probs.last().unwrap();
        probs.push(prev_prob * ratio);
        het += 2;
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

    // Mid-p: sum of probabilities strictly less extreme + 0.5 * obs_prob
    let mut p_value = 0.0;
    for (i, &p) in probs.iter().enumerate() {
        if i != obs_idx && p <= obs_prob + 1e-15 {
            p_value += p;
        }
    }
    p_value += 0.5 * obs_prob;

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
        // ref freq = (2*2 + 1) / (4*2) = 5/8 = 0.625
        assert!((s.ref_freq() - 0.625).abs() < 1e-10);
        assert!((s.alt_freq() - 0.375).abs() < 1e-10);
        assert!((s.maf() - 0.375).abs() < 1e-10);
        assert!((s.obs_het() - 0.25).abs() < 1e-10);
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
    fn hwe_monomorphic() {
        // All same genotype, no polymorphism
        assert_eq!(hwe_exact_midp(100, 0, 0), 1.0);
    }
}
