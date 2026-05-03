//! Sniff the DTC kit vendor from the header/first bytes of a raw file.
//!
//! Each vendor drops identifying text in the comment header (23andMe prefixes
//! lines with `#`, Ancestry with `#AncestryDNA`, FTDNA and MyHeritage have
//! plain CSV headers). Order matters: the more specific magic is checked
//! before the shared `RSID,CHROMOSOME,POSITION,RESULT` row that appears in
//! both FTDNA and MyHeritage files.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::Read;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Vendor {
    TwentyThreeAndMe,
    Ancestry,
    MyHeritage,
    LivingDna,
    Ftdna,
}

/// Sniff the vendor from the first 64 KiB of the file.
pub fn detect_vendor(path: &Path) -> Result<Vendor> {
    let mut f =
        File::open(path).with_context(|| format!("open {} for sniffing", path.display()))?;
    let mut buf = vec![0u8; 65536];
    let n = f.read(&mut buf)?;
    let s = String::from_utf8_lossy(&buf[..n]);
    detect_vendor_from_text(&s, Some(path))
}

pub fn detect_vendor_from_text(s: &str, path: Option<&Path>) -> Result<Vendor> {
    let lower = s.to_ascii_lowercase();
    let path_disp = path
        .map(|p| p.display().to_string())
        .unwrap_or_else(|| "<input>".to_string());

    // Most specific strings first. Vendor names inside header comment blocks
    // are the authoritative signal; the shared CSV header row is a fallback.
    if lower.contains("23andme") {
        return Ok(Vendor::TwentyThreeAndMe);
    }
    if lower.contains("ancestrydna") {
        return Ok(Vendor::Ancestry);
    }
    if lower.contains("family tree dna") || lower.contains("familytreedna") {
        return Ok(Vendor::Ftdna);
    }
    if lower.contains("myheritage") {
        return Ok(Vendor::MyHeritage);
    }
    if lower.contains("living dna") || lower.contains("livingdna") {
        return Ok(Vendor::LivingDna);
    }

    // Structural fallbacks when the header text was stripped.
    // 23andMe / LivingDNA share this TSV schema; require --vendor.
    if s.contains("rsid\tchromosome\tposition\tgenotype") {
        anyhow::bail!(
            "{}: ambiguous TSV header (23andMe and LivingDNA share this format). \
             Pass --vendor twenty-three-and-me or --vendor living-dna to disambiguate.",
            path_disp
        );
    }
    // AncestryDNA: TSV with split allele columns.
    if s.contains("rsid\tchromosome\tposition\tallele1\tallele2") {
        return Ok(Vendor::Ancestry);
    }
    // FTDNA and MyHeritage share this CSV header. We can't reliably
    // distinguish them without a vendor tag; require `--vendor`.
    if s.contains("RSID,CHROMOSOME,POSITION,RESULT") {
        anyhow::bail!(
            "{}: ambiguous CSV header (FTDNA and MyHeritage share this format). \
             Pass --vendor ftdna or --vendor myheritage to disambiguate.",
            path_disp
        );
    }

    anyhow::bail!("Could not detect vendor for {}", path_disp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn ambiguous_tsv_requires_vendor() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "rsid\tchromosome\tposition\tgenotype").unwrap();
        let err = detect_vendor(f.path()).unwrap_err().to_string();
        assert!(err.contains("ambiguous TSV header"));
    }

    #[test]
    fn detects_ancestry_from_vendor_tag() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "#AncestryDNA").unwrap();
        assert_eq!(detect_vendor(f.path()).unwrap(), Vendor::Ancestry);
    }
}
