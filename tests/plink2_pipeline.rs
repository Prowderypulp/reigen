//! End-to-end pipeline tests for PLINK 2 (PGEN/PVAR/PSAM) support — Track B.
//!
//! These drive the real `reigen` binary across the convert pipeline, exercising
//! format inference, the PVAR/PSAM readers (incl. `##` meta-line skipping and
//! the multiallelic keep-mask), and the PGEN reader/writer wiring.
//!
//! Polarity note: like `convertf_golden.rs`, the read-parity assertion compares
//! PGEN output against the *external* golden PAM fixture (not reigen-vs-reigen),
//! so a genome-wide allele flip would be caught.

use std::path::{Path, PathBuf};
use std::process::Command;

fn golden(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(rel)
}

fn reigen() -> Command {
    Command::new(env!("CARGO_BIN_EXE_reigen"))
}

/// PAM `.geno` record bytes, dropping the `rlen`-byte header (counts + hashes).
fn pam_records(path: &Path, nind: usize) -> Vec<u8> {
    let rlen = std::cmp::max(48, (nind * 2 + 7) / 8);
    std::fs::read(path).expect("read .geno")[rlen..].to_vec()
}

/// Read parity: PGEN → PAM must reproduce the golden PAM genotype records
/// byte-for-byte (same `g = count(allele1)`, ALT = allele1).
#[test]
fn pgen_to_pam_matches_golden() {
    let dir = tempfile::tempdir().unwrap();
    let geno = dir.path().join("out.geno");
    let snp = dir.path().join("out.snp");
    let ind = dir.path().join("out.ind");

    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            golden("plink2/gold.pgen").to_str().unwrap(),
            "--in-snp",
            golden("plink2/gold.pvar").to_str().unwrap(),
            "--in-ind",
            golden("plink2/gold.psam").to_str().unwrap(),
            "--out-format",
            "packedancestrymap",
            "--out-geno",
            geno.to_str().unwrap(),
            "--out-snp",
            snp.to_str().unwrap(),
            "--out-ind",
            ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert PGEN->PAM failed");

    assert_eq!(
        pam_records(&geno, 8),
        pam_records(&golden("pam/gold.geno"), 8),
        "PGEN->PAM genotype records must match the golden PAM (polarity + ALT=allele1)"
    );
}

/// `.pgen` input is inferred from the extension (no `--in-format` needed) and
/// the multiallelic (comma-ALT) variant is dropped, keeping SNP rows aligned
/// with the genotype stream.
#[test]
fn multiallelic_variant_is_dropped() {
    let dir = tempfile::tempdir().unwrap();
    let geno = dir.path().join("m.geno");
    let snp = dir.path().join("m.snp");
    let ind = dir.path().join("m.ind");

    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            golden("plink2/multi.pgen").to_str().unwrap(),
            "--in-snp",
            golden("plink2/multi.pvar").to_str().unwrap(),
            "--in-ind",
            golden("plink2/multi.psam").to_str().unwrap(),
            "--out-format",
            "eigenstrat",
            "--out-geno",
            geno.to_str().unwrap(),
            "--out-snp",
            snp.to_str().unwrap(),
            "--out-ind",
            ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert multi.pgen failed");

    // rs1 (ALT=A,T) dropped; rs0/rs2/rs3 kept, in order.
    let ids: Vec<String> = std::fs::read_to_string(&snp)
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().next().unwrap().to_string())
        .collect();
    assert_eq!(ids, ["rs0", "rs2", "rs3"]);

    // Each .geno line is one SNP across all 6 samples; 3 kept SNPs.
    let geno_txt = std::fs::read_to_string(&geno).unwrap();
    let lines: Vec<&str> = geno_txt.lines().collect();
    assert_eq!(lines.len(), 3, "one genotype row per kept SNP");
    for l in &lines {
        assert_eq!(l.trim().len(), 6, "one genotype char per sample");
    }

    // A `.pgen-drop.tsv` report (next to the output geno) lists the dropped
    // multiallelic variant — plan §4.4 / acceptance #5.
    let report = geno.with_extension("pgen-drop.tsv");
    let report_txt = std::fs::read_to_string(&report).expect("drop report written");
    let data: Vec<&str> = report_txt.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(data.len(), 1, "one dropped variant");
    assert!(data[0].contains("rs1"), "rs1 (ALT=A,T) reported");
    assert!(data[0].contains("multiallelic"), "reason recorded");
}

/// Write path: PAM → PGEN (mode 0x02) → PAM round-trips genotypes exactly, and
/// the emitted `.pgen` is a valid mode-0x02 file (magic `6c 1b 02`).
#[test]
fn pam_to_pgen_to_pam_roundtrip() {
    let dir = tempfile::tempdir().unwrap();
    let pgen = dir.path().join("w.pgen");
    let pvar = dir.path().join("w.pvar");
    let psam = dir.path().join("w.psam");

    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            golden("pam/gold.geno").to_str().unwrap(),
            "--in-snp",
            golden("pam/gold.snp").to_str().unwrap(),
            "--in-ind",
            golden("pam/gold.ind").to_str().unwrap(),
            "--out-format",
            "plink2",
            "--out-geno",
            pgen.to_str().unwrap(),
            "--out-snp",
            pvar.to_str().unwrap(),
            "--out-ind",
            psam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert PAM->PGEN failed");

    let bytes = std::fs::read(&pgen).unwrap();
    assert_eq!(&bytes[0..3], &[0x6c, 0x1b, 0x02], "PGEN mode-0x02 magic");

    // PGEN -> PAM and compare to the original golden PAM.
    let geno = dir.path().join("rt.geno");
    let snp = dir.path().join("rt.snp");
    let ind = dir.path().join("rt.ind");
    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            pgen.to_str().unwrap(),
            "--in-snp",
            pvar.to_str().unwrap(),
            "--in-ind",
            psam.to_str().unwrap(),
            "--out-format",
            "packedancestrymap",
            "--out-geno",
            geno.to_str().unwrap(),
            "--out-snp",
            snp.to_str().unwrap(),
            "--out-ind",
            ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert PGEN->PAM (roundtrip) failed");

    assert_eq!(
        pam_records(&geno, 8),
        pam_records(&golden("pam/gold.geno"), 8),
        "PAM->PGEN->PAM must preserve genotypes exactly"
    );
}

/// Merge accepts PGEN inputs (via temp-PAM auto-convert) and emits PLINK2
/// output. Merging the gold fixture with itself auto-renames duplicate sample
/// IDs, yielding 16 samples; the emitted `.pgen` must read back.
#[test]
fn merge_pgen_inputs_to_plink2() {
    let dir = tempfile::tempdir().unwrap();
    let triple = format!(
        "{}:{}:{}",
        golden("plink2/gold.pgen").to_str().unwrap(),
        golden("plink2/gold.pvar").to_str().unwrap(),
        golden("plink2/gold.psam").to_str().unwrap(),
    );
    let pgen = dir.path().join("m.pgen");
    let pvar = dir.path().join("m.pvar");
    let psam = dir.path().join("m.psam");

    let status = reigen()
        .args([
            "merge",
            "--in",
            &triple,
            "--in",
            &triple,
            "--out-format",
            "plink2",
            "--out-geno",
            pgen.to_str().unwrap(),
            "--out-snp",
            pvar.to_str().unwrap(),
            "--out-ind",
            psam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "merge PGEN->PLINK2 failed");

    let bytes = std::fs::read(&pgen).unwrap();
    assert_eq!(&bytes[0..3], &[0x6c, 0x1b, 0x02], "merge emits mode-0x02 PGEN");

    // 16 samples after duplicate-ID auto-rename.
    let n_samples = std::fs::read_to_string(&psam)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .count();
    assert_eq!(n_samples, 16);

    // Output must read back through the PGEN reader.
    let geno = dir.path().join("v.geno");
    let snp = dir.path().join("v.snp");
    let ind = dir.path().join("v.ind");
    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            pgen.to_str().unwrap(),
            "--in-snp",
            pvar.to_str().unwrap(),
            "--in-ind",
            psam.to_str().unwrap(),
            "--out-format",
            "eigenstrat",
            "--out-geno",
            geno.to_str().unwrap(),
            "--out-snp",
            snp.to_str().unwrap(),
            "--out-ind",
            ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "merged PGEN must read back");
}

/// The write→plink2 byte-equivalence oracle (Track C2).
/// Validates that reigen's emitted mode-0x02 PGEN is read by plink2
/// to produce a .bed/.bim that is byte-identical to the golden PLINK 1 output.
#[test]
fn pgen_write_plink2_equivalence_oracle() {
    // Only run if plink2 is available on PATH.
    if Command::new("plink2").arg("--version").output().is_err() {
        println!("Skipping pgen_write_plink2_equivalence_oracle: plink2 not found on PATH");
        return;
    }

    let dir = tempfile::tempdir().unwrap();
    let pgen = dir.path().join("out.pgen");
    let pvar = dir.path().join("out.pvar");
    let psam = dir.path().join("out.psam");

    // 1. Convert PAM golden to PGEN
    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            golden("pam/gold.geno").to_str().unwrap(),
            "--in-snp",
            golden("pam/gold.snp").to_str().unwrap(),
            "--in-ind",
            golden("pam/gold.ind").to_str().unwrap(),
            "--out-format",
            "plink2",
            "--out-geno",
            pgen.to_str().unwrap(),
            "--out-snp",
            pvar.to_str().unwrap(),
            "--out-ind",
            psam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert PAM->PGEN failed");

    // 2. Use plink2 to read the reigen-emitted PGEN and write a .bed
    let plink_out_prefix = dir.path().join("plink_out");
    let plink_status = Command::new("plink2")
        .args([
            "--pfile",
            dir.path().join("out").to_str().unwrap(),
            "--make-bed",
            "--out",
            plink_out_prefix.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(plink_status.success(), "plink2 failed to read reigen's PGEN output");

    let plink_bed = dir.path().join("plink_out.bed");
    let plink_bim = dir.path().join("plink_out.bim");

    let bed_bytes = std::fs::read(&plink_bed).unwrap();
    let golden_bed = std::fs::read(golden("plink/gold.bed")).unwrap();
    assert_eq!(bed_bytes, golden_bed, "byte-identical .bed equivalent");

    // .bim must be semantically identical (same A1/A2 = ALT/REF = allele1/allele2).
    // plink2 emits tab-separated .bim; the golden was produced by plink1 with
    // space-padding, so byte-identity is not expected. Compare field-by-field.
    let plink_bim_text = std::fs::read_to_string(&plink_bim).unwrap();
    let golden_bim_text = std::fs::read_to_string(golden("plink/gold.bim")).unwrap();
    let plink_fields: Vec<Vec<&str>> = plink_bim_text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split_whitespace().collect())
        .collect();
    let golden_fields: Vec<Vec<&str>> = golden_bim_text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split_whitespace().collect())
        .collect();
    assert_eq!(plink_fields.len(), golden_fields.len(), ".bim row count mismatch");
    for (i, (pf, gf)) in plink_fields.iter().zip(golden_fields.iter()).enumerate() {
        assert_eq!(pf.len(), 6, ".bim row {i}: expected 6 fields");
        assert_eq!(gf.len(), 6, ".bim row {i}: expected 6 fields (golden)");
        // Fields: chrom, id, cm, pos, a1, a2
        assert_eq!(pf[0], gf[0], ".bim row {i}: chrom mismatch");
        assert_eq!(pf[1], gf[1], ".bim row {i}: id mismatch");
        // CM is a float — compare numerically (plink2 may truncate trailing zeros).
        let cm_p: f64 = pf[2].parse().unwrap();
        let cm_g: f64 = gf[2].parse().unwrap();
        assert!((cm_p - cm_g).abs() < 1e-6, ".bim row {i}: CM mismatch ({cm_p} vs {cm_g})");
        assert_eq!(pf[3], gf[3], ".bim row {i}: pos mismatch");
        assert_eq!(pf[4], gf[4], ".bim row {i}: A1 (allele1/ALT) mismatch");
        assert_eq!(pf[5], gf[5], ".bim row {i}: A2 (allele2/REF) mismatch");
    }
}

/// Single-SNP anchor assertion explicitly verifying the identity mapping logic.
/// PVAR ALT=A/REF=G, hom-ALT sample => `g = 2`.
#[test]
fn pgen_single_snp_anchor() {
    let dir = tempfile::tempdir().unwrap();
    let geno = dir.path().join("out.geno");
    let snp = dir.path().join("out.snp");
    let ind = dir.path().join("out.ind");

    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            golden("plink2/gold.pgen").to_str().unwrap(),
            "--in-snp",
            golden("plink2/gold.pvar").to_str().unwrap(),
            "--in-ind",
            golden("plink2/gold.psam").to_str().unwrap(),
            "--out-format",
            "eigenstrat",
            "--out-geno",
            geno.to_str().unwrap(),
            "--out-snp",
            snp.to_str().unwrap(),
            "--out-ind",
            ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert failed");

    let geno_str = std::fs::read_to_string(&geno).unwrap();
    let first_line = geno_str.lines().next().unwrap();
    // In our golden plink2 data, rs0 is ALT=A, REF=G.
    // The 5th and 8th samples (indices 4 and 7) are hom-ALT ('2').
    // The first line should be "99002012".
    assert_eq!(first_line, "99002012", "hom-ALT sample => g = 2");
}

// ======================================================================
// Track C4 — Cross-format CLI smoke tests
// ======================================================================

/// Helper: count non-empty lines in a text file (SNP or sample count proxy).
fn count_lines(path: &Path) -> usize {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count()
}

/// Helper: count data lines (skip # header lines) in a PVAR/PSAM/VCF file.
fn count_data_lines(path: &Path) -> usize {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .count()
}

/// Full format loop: PAM → PGEN → PAM → TGENO → EIGENSTRAT → VCF.
/// At each step, verify the output is non-empty, re-readable, and SNP/sample
/// counts are stable (12 SNPs, 8 samples throughout).
#[test]
fn cross_format_loop_preserves_counts() {
    let dir = tempfile::tempdir().unwrap();

    // Step 1: PAM → PGEN
    let pgen = dir.path().join("s1.pgen");
    let pvar = dir.path().join("s1.pvar");
    let psam = dir.path().join("s1.psam");
    let s1 = reigen()
        .args([
            "convert",
            "--in-geno", golden("pam/gold.geno").to_str().unwrap(),
            "--in-snp",  golden("pam/gold.snp").to_str().unwrap(),
            "--in-ind",  golden("pam/gold.ind").to_str().unwrap(),
            "--out-format", "plink2",
            "--out-geno", pgen.to_str().unwrap(),
            "--out-snp",  pvar.to_str().unwrap(),
            "--out-ind",  psam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(s1.success(), "PAM→PGEN failed");
    assert_eq!(count_data_lines(&pvar), 12, "PGEN step: 12 SNPs in PVAR");
    assert_eq!(count_data_lines(&psam), 8, "PGEN step: 8 samples in PSAM");

    // Step 2: PGEN → PAM
    let geno2 = dir.path().join("s2.geno");
    let snp2 = dir.path().join("s2.snp");
    let ind2 = dir.path().join("s2.ind");
    let s2 = reigen()
        .args([
            "convert",
            "--in-geno", pgen.to_str().unwrap(),
            "--in-snp",  pvar.to_str().unwrap(),
            "--in-ind",  psam.to_str().unwrap(),
            "--out-format", "packedancestrymap",
            "--out-geno", geno2.to_str().unwrap(),
            "--out-snp",  snp2.to_str().unwrap(),
            "--out-ind",  ind2.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(s2.success(), "PGEN→PAM failed");
    assert_eq!(count_lines(&snp2), 12, "PAM step: 12 SNPs");
    assert_eq!(count_lines(&ind2), 8, "PAM step: 8 samples");

    // Step 3: PAM → TGENO
    let geno3 = dir.path().join("s3.geno");
    let snp3 = dir.path().join("s3.snp");
    let ind3 = dir.path().join("s3.ind");
    let s3 = reigen()
        .args([
            "convert",
            "--in-geno", geno2.to_str().unwrap(),
            "--in-snp",  snp2.to_str().unwrap(),
            "--in-ind",  ind2.to_str().unwrap(),
            "--out-format", "tgeno",
            "--out-geno", geno3.to_str().unwrap(),
            "--out-snp",  snp3.to_str().unwrap(),
            "--out-ind",  ind3.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(s3.success(), "PAM→TGENO failed");
    assert_eq!(count_lines(&snp3), 12, "TGENO step: 12 SNPs");
    assert_eq!(count_lines(&ind3), 8, "TGENO step: 8 samples");

    // Step 4: TGENO → EIGENSTRAT
    let geno4 = dir.path().join("s4.geno");
    let snp4 = dir.path().join("s4.snp");
    let ind4 = dir.path().join("s4.ind");
    let s4 = reigen()
        .args([
            "convert",
            "--in-geno", geno3.to_str().unwrap(),
            "--in-snp",  snp3.to_str().unwrap(),
            "--in-ind",  ind3.to_str().unwrap(),
            "--out-format", "eigenstrat",
            "--out-geno", geno4.to_str().unwrap(),
            "--out-snp",  snp4.to_str().unwrap(),
            "--out-ind",  ind4.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(s4.success(), "TGENO→EIGENSTRAT failed");
    assert_eq!(count_lines(&snp4), 12, "EIGENSTRAT step: 12 SNPs");
    assert_eq!(count_lines(&ind4), 8, "EIGENSTRAT step: 8 samples");

    // Step 5: EIGENSTRAT → VCF
    let vcf5 = dir.path().join("s5.vcf");
    let s5 = reigen()
        .args([
            "export",
            "--geno", geno4.to_str().unwrap(),
            "--snp",  snp4.to_str().unwrap(),
            "--ind",  ind4.to_str().unwrap(),
            "-o", vcf5.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(s5.success(), "EIGENSTRAT→VCF failed");
    let vcf_data_lines = count_data_lines(&vcf5);
    assert_eq!(vcf_data_lines, 12, "VCF step: 12 SNP data lines");

    // Final parity: convert the EIGENSTRAT output to PAM and compare with golden.
    let geno_final = dir.path().join("final.geno");
    let snp_final = dir.path().join("final.snp");
    let ind_final = dir.path().join("final.ind");
    let s_final = reigen()
        .args([
            "convert",
            "--in-geno", geno4.to_str().unwrap(),
            "--in-snp",  snp4.to_str().unwrap(),
            "--in-ind",  ind4.to_str().unwrap(),
            "--out-format", "packedancestrymap",
            "--out-geno", geno_final.to_str().unwrap(),
            "--out-snp",  snp_final.to_str().unwrap(),
            "--out-ind",  ind_final.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(s_final.success(), "EIGENSTRAT→PAM final step failed");
    assert_eq!(
        pam_records(&geno_final, 8),
        pam_records(&golden("pam/gold.geno"), 8),
        "Full loop: final PAM must match golden PAM genotypes"
    );
}

/// `convert` subcommand on PGEN input → PLINK1 output.
#[test]
fn convert_pgen_to_plink1() {
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("out.bed");
    let bim = dir.path().join("out.bim");
    let fam = dir.path().join("out.fam");

    let status = reigen()
        .args([
            "convert",
            "--in-geno", golden("plink2/gold.pgen").to_str().unwrap(),
            "--in-snp",  golden("plink2/gold.pvar").to_str().unwrap(),
            "--in-ind",  golden("plink2/gold.psam").to_str().unwrap(),
            "--out-format", "packedped",
            "--out-geno", bed.to_str().unwrap(),
            "--out-snp",  bim.to_str().unwrap(),
            "--out-ind",  fam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "convert PGEN→PLINK1 failed");
    assert_eq!(count_lines(&bim), 12);
    assert_eq!(count_lines(&fam), 8);

    // Re-read the output to verify it's valid.
    let geno = dir.path().join("rt.geno");
    let snp = dir.path().join("rt.snp");
    let ind = dir.path().join("rt.ind");
    let rt = reigen()
        .args([
            "convert",
            "--in-geno", bed.to_str().unwrap(),
            "--in-snp",  bim.to_str().unwrap(),
            "--in-ind",  fam.to_str().unwrap(),
            "--out-format", "packedancestrymap",
            "--out-geno", geno.to_str().unwrap(),
            "--out-snp",  snp.to_str().unwrap(),
            "--out-ind",  ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(rt.success(), "re-read converted PLINK1 failed");
    assert_eq!(
        pam_records(&geno, 8),
        pam_records(&golden("pam/gold.geno"), 8),
        "PGEN→PLINK1→PAM must match golden"
    );
}

/// `filter` subcommand on PGEN input.
#[test]
fn filter_pgen_input() {
    let dir = tempfile::tempdir().unwrap();
    let geno = dir.path().join("out.geno");
    let snp = dir.path().join("out.snp");
    let ind = dir.path().join("out.ind");

    let status = reigen()
        .args([
            "filter",
            "--in-geno", golden("plink2/gold.pgen").to_str().unwrap(),
            "--in-snp",  golden("plink2/gold.pvar").to_str().unwrap(),
            "--in-ind",  golden("plink2/gold.psam").to_str().unwrap(),
            "--out-format", "eigenstrat",
            "--out-geno", geno.to_str().unwrap(),
            "--out-snp",  snp.to_str().unwrap(),
            "--out-ind",  ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "filter PGEN failed");
    assert_eq!(count_lines(&snp), 12, "filter: 12 SNPs");
    assert_eq!(count_lines(&ind), 8, "filter: 8 samples");
}

/// `stats` subcommand on PGEN input — must succeed and produce output.
#[test]
fn stats_pgen_input() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("stats");

    let status = reigen()
        .args([
            "stats",
            "-i", golden("plink2/gold").to_str().unwrap(),
            "-o", prefix.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "stats PGEN failed");

    let snp_stats = dir.path().join("stats.snp_stats.tsv");
    let text = std::fs::read_to_string(&snp_stats).unwrap();
    assert!(!text.trim().is_empty(), "stats output must be non-empty");
    // Should have header + 12 data rows (header doesn't start with #).
    let data_lines: Vec<&str> = text.lines().skip(1).filter(|l| !l.trim().is_empty()).collect();
    assert_eq!(data_lines.len(), 12, "stats: 12 SNP rows");
}

/// `export` subcommand on PGEN input → VCF.
#[test]
fn export_pgen_to_vcf() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("out.vcf");

    let status = reigen()
        .args([
            "export",
            "-i", golden("plink2/gold").to_str().unwrap(),
            "-o", vcf.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "export PGEN→VCF failed");

    let vcf_text = std::fs::read_to_string(&vcf).unwrap();
    let data_lines: Vec<&str> = vcf_text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(data_lines.len(), 12, "VCF: 12 SNP data lines");
}

/// Regression: PGEN→PLINK1 must preserve sex/PAR (XY) and mitochondrial (MT)
/// chromosomes. PVAR and BIM use PLINK numeric codes (25=XY, 26=MT) but reigen's
/// internal codes are 25=MT/26=XY; PVAR once skipped that remap and silently
/// swapped XY↔MT on conversion. A hand-built PGEN fileset with X/Y/XY/MT must
/// round-trip the numeric codes unchanged into the `.bim`.
#[test]
fn pgen_chrom_xy_mt_preserved_to_plink1() {
    let dir = tempfile::tempdir().unwrap();
    let pgen = dir.path().join("c.pgen");
    let pvar = dir.path().join("c.pvar");
    let psam = dir.path().join("c.psam");

    // 4 variants (X, Y, XY, MT), 1 sample, fixed-width mode-0x02 PGEN.
    std::fs::write(
        &pvar,
        "#CHROM\tPOS\tID\tREF\tALT\n\
         23\t1\trX\tG\tA\n24\t2\trY\tG\tA\n25\t3\trXY\tG\tA\n26\t4\trMT\tG\tA\n",
    )
    .unwrap();
    std::fs::write(&psam, "#IID\tSEX\nS0\t0\n").unwrap();
    // header: 6c 1b 02 | M=4 (LE u32) | N=1 (LE u32) | 0x80 | 4 record bytes
    std::fs::write(
        &pgen,
        [
            0x6c, 0x1b, 0x02, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x00,
        ],
    )
    .unwrap();

    let bim = dir.path().join("c.bim");
    let status = reigen()
        .args([
            "convert",
            "--in-geno", pgen.to_str().unwrap(),
            "--in-snp",  pvar.to_str().unwrap(),
            "--in-ind",  psam.to_str().unwrap(),
            "--out-format", "packedped",
            "--out-geno", dir.path().join("c.bed").to_str().unwrap(),
            "--out-snp",  bim.to_str().unwrap(),
            "--out-ind",  dir.path().join("c.fam").to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "PGEN→PLINK1 with X/Y/XY/MT failed");

    let chroms: Vec<String> = std::fs::read_to_string(&bim)
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().next().unwrap().to_string())
        .collect();
    assert_eq!(chroms, ["23", "24", "25", "26"], "XY/MT must not be swapped");
}
