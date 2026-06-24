//! Ground-truth regression test against EIGENSOFT `convertf` output.
//!
//! Every other test in this crate is a reigen→reigen round-trip. Such tests
//! are invariant under a *global* genotype/allele inversion, so they cannot
//! detect the genotype-polarity bug that produced spurious qpAdm results
//! (PACKEDANCESTRYMAP `.geno` counts column 5 / allele1, but reigen had once
//! encoded PLINK output as if it counted allele2, via an erroneous `.bim`
//! A1/A2 swap).
//!
//! These tests convert *convertf-produced* fixtures (`tests/golden/`) with the
//! real reigen binary and assert byte-for-byte agreement with convertf. The
//! fixtures are committed, so no external tools are needed at test time.

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

/// PACKEDANCESTRYMAP `.geno` record bytes, dropping the `rlen`-byte header
/// (which carries counts + hashes that legitimately differ from convertf).
fn pam_records(path: &Path, nind: usize) -> Vec<u8> {
    let rlen = std::cmp::max(48, (nind * 2 + 7) / 8);
    let bytes = std::fs::read(path).expect("read .geno");
    bytes[rlen..].to_vec()
}

/// Columns 5 and 6 (alleles) of each `.snp`/`.bim` line, in file order.
fn alleles(path: &Path) -> Vec<(String, String)> {
    std::fs::read_to_string(path)
        .expect("read snp/bim")
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split_whitespace().collect();
            (c[4].to_string(), c[5].to_string())
        })
        .collect()
}

const NIND: usize = 8; // golden dataset (multiple of 4 → no .bed padding)

#[test]
fn pam_to_plink_matches_convertf() {
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("out.bed");
    let bim = dir.path().join("out.bim");
    let fam = dir.path().join("out.fam");

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
            "packedped",
            "--out-geno",
            bed.to_str().unwrap(),
            "--out-snp",
            bim.to_str().unwrap(),
            "--out-ind",
            fam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert PAM->PLINK failed");

    // nind is a multiple of 4, so the `.bed` has no partial trailing byte and
    // must match convertf exactly (magic + SNP-major 2-bit records).
    let got = std::fs::read(&bed).unwrap();
    let want = std::fs::read(golden("plink/gold.bed")).unwrap();
    assert_eq!(
        got, want,
        "reigen .bed must be byte-identical to convertf (genotype polarity + A1 ordering)"
    );

    // Allele columns: convertf keeps allele1 in A1 (col 5) — no swap.
    assert_eq!(
        alleles(&bim),
        alleles(&golden("plink/gold.bim")),
        "reigen .bim alleles (A1=col5=allele1) must match convertf"
    );
}

#[test]
fn plink_to_pam_matches_convertf() {
    let dir = tempfile::tempdir().unwrap();
    let geno = dir.path().join("out.geno");
    let snp = dir.path().join("out.snp");
    let ind = dir.path().join("out.ind");

    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            golden("plink/gold.bed").to_str().unwrap(),
            "--in-snp",
            golden("plink/gold.bim").to_str().unwrap(),
            "--in-ind",
            golden("plink/gold.fam").to_str().unwrap(),
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
    assert!(status.success(), "reigen convert PLINK->PAM failed");

    // Genotype records (post-header) must be byte-identical to convertf: the
    // `.geno` value counts column 5 (allele1), `2` = homozygous allele1.
    assert_eq!(
        pam_records(&geno, NIND),
        pam_records(&golden("pam/gold.geno"), NIND),
        "reigen .geno records must be byte-identical to convertf"
    );
    assert_eq!(
        alleles(&snp),
        alleles(&golden("pam/gold.snp")),
        "reigen .snp alleles must match convertf"
    );
}

/// Independent semantic anchor: a hand-built single-SNP PAM where we *know*
/// the intended genotype, asserting `geno = 2` is homozygous for `.snp`
/// column 5 (not column 6) after conversion to PLINK.
#[test]
fn geno_two_is_homozygous_column_five() {
    let dir = tempfile::tempdir().unwrap();
    // One SNP, col5=A (allele1), col6=G (allele2). 4 samples: 2,2,0,0.
    let g = dir.path().join("t.geno");
    let s = dir.path().join("t.snp");
    let i = dir.path().join("t.ind");
    std::fs::write(&g, "2200\n").unwrap();
    std::fs::write(&s, "rs1  1  0.0  1000 A G\n").unwrap();
    std::fs::write(&i, "I0 U P\nI1 U P\nI2 U P\nI3 U P\n").unwrap();

    let bed = dir.path().join("o.bed");
    let bim = dir.path().join("o.bim");
    let fam = dir.path().join("o.fam");
    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            g.to_str().unwrap(),
            "--in-snp",
            s.to_str().unwrap(),
            "--in-ind",
            i.to_str().unwrap(),
            "--out-format",
            "packedped",
            "--out-geno",
            bed.to_str().unwrap(),
            "--out-snp",
            bim.to_str().unwrap(),
            "--out-ind",
            fam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    // .bim: A1 (col5) = A (allele1), A2 = G.
    let line = std::fs::read_to_string(&bim).unwrap();
    let cols: Vec<&str> = line.split_whitespace().collect();
    assert_eq!(cols[4], "A", "A1 must be allele1 (col5), no swap");
    assert_eq!(cols[5], "G");

    // .bed: magic (3) + one record byte. Samples 0,1 = geno 2 → hom A1 (`00`);
    // samples 2,3 = geno 0 → hom A2 (`11`). PLINK LSB-first:
    // s0=00, s1=00, s2=11, s3=11 → 0b1111_0000 = 0xF0.
    let bed_bytes = std::fs::read(&bed).unwrap();
    assert_eq!(&bed_bytes[0..3], &[0x6c, 0x1b, 0x01], "bed magic");
    assert_eq!(
        bed_bytes[3], 0xF0,
        "geno=2 must encode as PLINK hom-A1 (00); got {:#010b}",
        bed_bytes[3]
    );
}

#[test]
fn plink_to_plink_preserves_iid_without_reprefixing_family_name() {
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("out.bed");
    let bim = dir.path().join("out.bim");
    let fam = dir.path().join("out.fam");

    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            golden("plink/gold.bed").to_str().unwrap(),
            "--in-snp",
            golden("plink/gold.bim").to_str().unwrap(),
            "--in-ind",
            golden("plink/gold.fam").to_str().unwrap(),
            "--out-format",
            "packedped",
            "--out-geno",
            bed.to_str().unwrap(),
            "--out-snp",
            bim.to_str().unwrap(),
            "--out-ind",
            fam.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert PLINK->PLINK failed");

    let got: Vec<Vec<String>> = std::fs::read_to_string(&fam)
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .collect();
    let want: Vec<Vec<String>> = std::fs::read_to_string(golden("plink/gold.fam"))
        .unwrap()
        .lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .collect();

    assert_eq!(got.len(), want.len());
    for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
        assert_eq!(g[0], w[0], ".fam row {i}: FID changed");
        assert_eq!(g[1], w[1], ".fam row {i}: IID changed");
    }
}
