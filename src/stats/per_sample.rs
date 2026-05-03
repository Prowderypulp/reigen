//! Per-sample QC statistics: total genotyped SNPs, missingness rate,
//! observed heterozygosity, and inbreeding coefficient (F-statistic).

/// Accumulated counters for a single sample across all kept SNPs.
#[derive(Debug, Clone, Default)]
pub struct SampleStats {
    pub n_hom_ref: u32,
    pub n_het: u32,
    pub n_hom_alt: u32,
    pub n_missing: u32,
    /// Sum of expected heterozygosity across all non-missing SNPs for this
    /// sample, used for computing F-statistic. Accumulated as
    /// 2 * p * (1 - p) where p is the ref allele frequency at each SNP.
    pub sum_exp_het: f64,
}

impl SampleStats {
    #[inline]
    pub fn n_called(&self) -> u32 {
        self.n_hom_ref + self.n_het + self.n_hom_alt
    }

    #[inline]
    pub fn n_total(&self) -> u32 {
        self.n_called() + self.n_missing
    }

    pub fn miss_rate(&self) -> f64 {
        let total = self.n_total();
        if total == 0 {
            return 0.0;
        }
        self.n_missing as f64 / total as f64
    }

    /// Observed heterozygosity rate for this sample.
    pub fn obs_het(&self) -> f64 {
        let n = self.n_called();
        if n == 0 {
            return f64::NAN;
        }
        self.n_het as f64 / n as f64
    }

    /// Inbreeding coefficient (F-statistic) using method of moments:
    ///   F = 1 - (obs_het / exp_het)
    /// where exp_het is the mean expected heterozygosity.
    pub fn f_stat(&self) -> f64 {
        let n = self.n_called();
        if n == 0 || self.sum_exp_het <= 0.0 {
            return f64::NAN;
        }
        let mean_exp_het = self.sum_exp_het / n as f64;
        if mean_exp_het <= 0.0 {
            return f64::NAN;
        }
        1.0 - (self.obs_het() / mean_exp_het)
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

/// Write per-sample stats as TSV.
pub fn write_sample_stats(
    path: &std::path::Path,
    ind_rows: &[crate::meta::IndRow],
    stats: &[SampleStats],
) -> anyhow::Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    writeln!(
        w,
        "SAMPLE\tPOP\tN_CALLED\tN_MISS\tMISS_RATE\tOBS_HET\tF_STAT"
    )?;

    for (ind, st) in ind_rows.iter().zip(stats.iter()) {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}",
            ind.id,
            ind.pop,
            st.n_called(),
            st.n_missing,
            st.miss_rate(),
            st.obs_het(),
            st.f_stat(),
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
    fn sample_stats_basic() {
        let mut s = SampleStats::default();
        for _ in 0..40 {
            s.observe(0);
        }
        for _ in 0..10 {
            s.observe(1);
        }
        for _ in 0..5 {
            s.observe(2);
        }
        for _ in 0..5 {
            s.observe(codec::G_MISSING);
        }

        assert_eq!(s.n_called(), 55);
        assert_eq!(s.n_total(), 60);
        assert!((s.miss_rate() - 5.0 / 60.0).abs() < 1e-10);
        assert!((s.obs_het() - 10.0 / 55.0).abs() < 1e-10);
    }

    #[test]
    fn f_stat_outbred() {
        // If obs_het == exp_het, F = 0
        let mut s = SampleStats::default();
        s.n_hom_ref = 25;
        s.n_het = 50;
        s.n_hom_alt = 25;
        // exp_het = 2*0.5*0.5 = 0.5 for each SNP
        s.sum_exp_het = 100.0 * 0.5; // 100 SNPs, each with exp_het = 0.5
        assert!((s.f_stat() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn f_stat_inbred() {
        // No hets at all → F = 1
        let mut s = SampleStats::default();
        s.n_hom_ref = 50;
        s.n_het = 0;
        s.n_hom_alt = 50;
        s.sum_exp_het = 100.0 * 0.5;
        assert!((s.f_stat() - 1.0).abs() < 1e-10);
    }
}
