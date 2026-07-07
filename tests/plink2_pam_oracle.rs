//! Targeted PLINK2 (PGEN/PVAR/PSAM) <-> PAM (PACKEDANCESTRYMAP) conversion tests,
//! validated against **standard `plink2`** as the external oracle.
//!
//! These are deliberately narrower and stronger than `plink2_pipeline.rs`: instead
//! of comparing reigen against a committed golden fixture, they use real `plink2`
//! on *both* ends and check every genotype cell against an independently-authored
//! ground-truth matrix. Both directions are covered:
//!
//!   * READ  (PGEN -> PAM): `plink2` builds the PGEN from a VCF we author, reigen
//!     converts PGEN -> PAM, and we decode the PAM and compare cell-by-cell.
//!   * WRITE (PAM -> PGEN): reigen writes the PGEN (via EIGENSTRAT -> PAM -> PGEN),
//!     `plink2` reads it and exports back to VCF, and we compare GT to truth.
//!
//! Polarity (the B-013 trap): reigen's `g = count(allele1)` and `allele1 = ALT`.
//! Truth here is authored as **ALT dosage** and checked against the ALT allele
//! that each tool actually reports (PAM `.snp` col-5, exported-VCF `ALT` col), so a
//! genome-wide genotype flip or an A1/A2 swap fails the test rather than hiding in
//! a self-consistent round trip. Tests self-skip when `plink2` is not on PATH.

use std::path::{Path, PathBuf};
use std::process::Command;

const MISS: u8 = 9;

fn reigen() -> Command {
    Command::new(env!("CARGO_BIN_EXE_reigen"))
}

fn plink2_available() -> bool {
    Command::new("plink2").arg("--version").output().is_ok()
}

/// One biallelic variant plus its per-sample ALT dosage (0/1/2 or `MISS`).
struct Variant {
    id: &'static str,
    chrom: &'static str,
    pos: u32,
    cm: f64,
    reff: char,
    alt: char,
    /// ALT dosage per sample; `MISS` = missing call.
    dose: Vec<u8>,
}

/// Ground-truth matrix: 8 samples x 12 variants spanning every hardcall class
/// (hom-REF, het, hom-ALT, missing), two chromosomes, and every biallelic
/// nucleotide pairing. Column j is sample j across all variants.
fn truth() -> (Vec<String>, Vec<Variant>) {
    let samples: Vec<String> = (0..8).map(|i| format!("S{i}")).collect();
    // dose rows are length 8 (one per sample).
    let v = |id, chrom, pos, cm, reff, alt, dose: [u8; 8]| Variant {
        id,
        chrom,
        pos,
        cm,
        reff,
        alt,
        dose: dose.to_vec(),
    };
    let variants = vec![
        v("rs0", "1", 1000, 0.001, 'G', 'A', [0, 1, 2, 9, 0, 1, 2, 0]),
        v("rs1", "1", 2000, 0.002, 'C', 'T', [2, 2, 1, 1, 9, 0, 0, 1]),
        v("rs2", "1", 3000, 0.003, 'A', 'C', [1, 0, 1, 2, 2, 1, 0, 9]),
        v("rs3", "1", 4000, 0.004, 'T', 'G', [9, 9, 0, 0, 1, 1, 2, 2]),
        v("rs4", "1", 5000, 0.005, 'A', 'T', [0, 0, 0, 0, 0, 0, 0, 1]),
        v("rs5", "1", 6000, 0.006, 'C', 'G', [2, 2, 2, 2, 2, 2, 2, 1]),
        v("rs6", "1", 7000, 0.007, 'G', 'C', [1, 1, 1, 1, 1, 1, 1, 1]),
        v("rs7", "2", 8000, 0.008, 'C', 'T', [9, 9, 9, 9, 1, 2, 0, 1]),
        v("rs8", "2", 9000, 0.009, 'A', 'G', [0, 1, 2, 0, 1, 2, 0, 1]),
        v("rs9", "2", 10000, 0.010, 'T', 'C', [2, 1, 0, 9, 2, 1, 0, 9]),
        v("rs10", "2", 11000, 0.011, 'A', 'C', [1, 2, 0, 1, 2, 0, 1, 2]),
        v("rs11", "2", 12000, 0.012, 'G', 'T', [0, 0, 1, 1, 2, 2, 9, 9]),
    ];
    (samples, variants)
}

// ---------------------------------------------------------------------------
// Independent reference codecs (do not touch reigen internals)
// ---------------------------------------------------------------------------

/// Decode a PAM `.geno` file into per-variant ALT dosages (`g = count(allele1)`,
/// 0/1/2 or `MISS`). PAM records are 2-bit **MSB-first** (sample 0 in bits 7-6),
/// after an `rlen`-byte header: `rlen = max(48, ceil(nind*2/8))`.
fn decode_pam(path: &Path, nind: usize, nsnp: usize) -> Vec<Vec<u8>> {
    // Header is `rlen` bytes, and every data record is *also* padded to `rlen`.
    let rlen = std::cmp::max(48, (nind * 2 + 7) / 8);
    let raw = std::fs::read(path).expect("read .geno");
    let rec = rlen;
    let body = &raw[rlen..];
    (0..nsnp)
        .map(|s| {
            let chunk = &body[s * rec..(s + 1) * rec];
            (0..nind)
                .map(|i| {
                    let two = (chunk[i / 4] >> (6 - 2 * (i % 4))) & 0b11;
                    match two {
                        0b00 => 0,
                        0b01 => 1,
                        0b10 => 2,
                        _ => MISS,
                    }
                })
                .collect()
        })
        .collect()
}

/// Parse `col5` (allele1) and `col6` (allele2) of a reigen `.snp` file.
fn snp_alleles(path: &Path) -> Vec<(char, char)> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            (f[4].chars().next().unwrap(), f[5].chars().next().unwrap())
        })
        .collect()
}

/// Author a minimal VCF (GT-only) from the truth matrix. Missing = `./.`.
fn write_vcf(path: &Path, samples: &[String], variants: &[Variant]) {
    let mut s = String::from("##fileformat=VCFv4.2\n");
    for c in ["1", "2"] {
        s.push_str(&format!("##contig=<ID={c}>\n"));
    }
    s.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for smp in samples {
        s.push('\t');
        s.push_str(smp);
    }
    s.push('\n');
    for v in variants {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t.\t.\t.\tGT",
            v.chrom, v.pos, v.id, v.reff, v.alt
        ));
        for &d in &v.dose {
            let gt = match d {
                0 => "0/0",
                1 => "0/1",
                2 => "1/1",
                _ => "./.",
            };
            s.push('\t');
            s.push_str(gt);
        }
        s.push('\n');
    }
    std::fs::write(path, s).unwrap();
}

/// Author an EIGENSTRAT `.geno`/`.snp`/`.ind` trio with `allele1 = ALT`, so that
/// the canonical `g = count(allele1)` equals the authored ALT dosage.
fn write_eigenstrat(dir: &Path, samples: &[String], variants: &[Variant]) -> (PathBuf, PathBuf, PathBuf) {
    let geno = dir.join("src.geno");
    let snp = dir.join("src.snp");
    let ind = dir.join("src.ind");

    let mut g = String::new();
    for v in variants {
        for &d in &v.dose {
            g.push(match d {
                0 => '0',
                1 => '1',
                2 => '2',
                _ => '9',
            });
        }
        g.push('\n');
    }
    std::fs::write(&geno, g).unwrap();

    let mut s = String::new();
    for v in variants {
        // id chrom cm pos allele1(=ALT) allele2(=REF)
        s.push_str(&format!(
            "{} {} {:.6} {} {} {}\n",
            v.id, v.chrom, v.cm, v.pos, v.alt, v.reff
        ));
    }
    std::fs::write(&snp, s).unwrap();

    let mut i = String::new();
    for smp in samples {
        i.push_str(&format!("{smp} U Pop0\n"));
    }
    std::fs::write(&ind, i).unwrap();

    (geno, snp, ind)
}

/// Run `reigen convert` and assert success.
fn convert(in_geno: &Path, in_snp: &Path, in_ind: &Path, out_fmt: &str, out_geno: &Path, out_snp: &Path, out_ind: &Path) {
    let status = reigen()
        .args([
            "convert",
            "--in-geno",
            in_geno.to_str().unwrap(),
            "--in-snp",
            in_snp.to_str().unwrap(),
            "--in-ind",
            in_ind.to_str().unwrap(),
            "--out-format",
            out_fmt,
            "--out-geno",
            out_geno.to_str().unwrap(),
            "--out-snp",
            out_snp.to_str().unwrap(),
            "--out-ind",
            out_ind.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "reigen convert -> {out_fmt} failed");
}

// ---------------------------------------------------------------------------
// READ: plink2-built PGEN -> reigen PAM, checked cell-by-cell against truth.
// ---------------------------------------------------------------------------

#[test]
fn read_pgen_to_pam_matches_plink2_oracle() {
    if !plink2_available() {
        eprintln!("SKIP read_pgen_to_pam_matches_plink2_oracle: plink2 not on PATH");
        return;
    }
    let (samples, variants) = truth();
    let dir = tempfile::tempdir().unwrap();

    // 1. Author VCF, let *standard plink2* build the PGEN/PVAR/PSAM.
    let vcf = dir.path().join("src.vcf");
    write_vcf(&vcf, &samples, &variants);
    let oracle = dir.path().join("oracle");
    let ok = Command::new("plink2")
        .args(["--vcf", vcf.to_str().unwrap(), "--make-pgen", "--out", oracle.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(ok.success(), "plink2 --make-pgen failed");

    // 2. reigen: PGEN -> PAM.
    let (geno, snp, ind) = (dir.path().join("out.geno"), dir.path().join("out.snp"), dir.path().join("out.ind"));
    convert(
        &oracle.with_extension("pgen"),
        &oracle.with_extension("pvar"),
        &oracle.with_extension("psam"),
        "packedancestrymap",
        &geno,
        &snp,
        &ind,
    );

    // 3. reigen must declare allele1 = ALT, allele2 = REF (no swap).
    let alleles = snp_alleles(&snp);
    assert_eq!(alleles.len(), variants.len(), "snp row count");
    for (v, (a1, a2)) in variants.iter().zip(&alleles) {
        assert_eq!(*a1, v.alt, "{}: allele1 must be ALT (polarity)", v.id);
        assert_eq!(*a2, v.reff, "{}: allele2 must be REF (polarity)", v.id);
    }

    // 4. Every genotype cell (g = ALT dosage) must match truth.
    let got = decode_pam(&geno, samples.len(), variants.len());
    for (v, row) in variants.iter().zip(&got) {
        assert_eq!(row, &v.dose, "{}: PGEN->PAM genotypes differ from truth", v.id);
    }
}

// ---------------------------------------------------------------------------
// WRITE: reigen-written PGEN -> plink2, exported back to VCF, checked vs truth.
// EIGENSTRAT -> PAM -> PGEN keeps PAM in the write path under test.
// ---------------------------------------------------------------------------

#[test]
fn write_pam_to_pgen_matches_plink2_oracle() {
    if !plink2_available() {
        eprintln!("SKIP write_pam_to_pgen_matches_plink2_oracle: plink2 not on PATH");
        return;
    }
    let (samples, variants) = truth();
    let dir = tempfile::tempdir().unwrap();

    // 1. Author EIGENSTRAT (allele1 = ALT), reigen -> PAM.
    let (eg, es, ei) = write_eigenstrat(dir.path(), &samples, &variants);
    let (pgeno, psnp, pind) = (dir.path().join("pam.geno"), dir.path().join("pam.snp"), dir.path().join("pam.ind"));
    convert(&eg, &es, &ei, "packedancestrymap", &pgeno, &psnp, &pind);

    // Sanity: the PAM itself already carries the right genotypes.
    let pam = decode_pam(&pgeno, samples.len(), variants.len());
    for (v, row) in variants.iter().zip(&pam) {
        assert_eq!(row, &v.dose, "{}: EIGENSTRAT->PAM genotypes differ from truth", v.id);
    }

    // 2. reigen: PAM -> PGEN (the write path under test).
    let out = dir.path().join("reigen");
    convert(
        &pgeno,
        &psnp,
        &pind,
        "plink2",
        &out.with_extension("pgen"),
        &out.with_extension("pvar"),
        &out.with_extension("psam"),
    );
    let magic = std::fs::read(out.with_extension("pgen")).unwrap();
    assert_eq!(&magic[0..3], &[0x6c, 0x1b, 0x02], "reigen PGEN mode-0x02 magic");

    // 3. *Standard plink2* reads reigen's PGEN and exports back to VCF.
    let back = dir.path().join("back");
    let ok = Command::new("plink2")
        .args(["--pfile", out.to_str().unwrap(), "--export", "vcf", "--out", back.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(ok.success(), "plink2 failed to read reigen's PGEN");

    // 4. Parse exported VCF; ALT-dosage from GT must match truth, and REF/ALT
    //    labels must be preserved (a swap would flip the dosage semantics).
    let vcf = std::fs::read_to_string(back.with_extension("vcf")).unwrap();
    let rows: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#') && !l.is_empty()).collect();
    assert_eq!(rows.len(), variants.len(), "exported variant count");
    for (v, line) in variants.iter().zip(&rows) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[3], v.reff.to_string(), "{}: exported REF must equal truth REF", v.id);
        assert_eq!(f[4], v.alt.to_string(), "{}: exported ALT must equal truth ALT", v.id);
        for (j, cell) in f[9..].iter().enumerate() {
            let gt = cell.split(':').next().unwrap();
            let got = alt_dosage(gt);
            assert_eq!(got, v.dose[j], "{} sample {}: exported GT {gt} != truth", v.id, samples[j]);
        }
    }
}

/// ALT dosage from a VCF GT string (`0/0`,`0/1`,`1/1`,`./.`, phased or not).
fn alt_dosage(gt: &str) -> u8 {
    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
    if alleles.iter().any(|a| *a == ".") {
        return MISS;
    }
    alleles.iter().filter(|a| **a == "1").count() as u8
}
