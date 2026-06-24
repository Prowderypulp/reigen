//! Edge-case conversion tests for the binary genotype codecs.
//!
//! Strategy (deliberately different from the convertf Python oracle): we
//! hand-author only the **hash-free** PLINK `.bed`/`.bim`/`.fam` fileset, drive
//! the real `reigen` binary, and assert on decoded genotype *dosages*
//! (number of copies of allele1; `9` = missing). Intermediate PAM/PGEN/TGENO
//! files are produced by reigen itself, so their headers/hashes are always
//! valid — we never forge a PAM header hash (which reigen rightly rejects).
//!
//! Ideas exercised (see design-edge-case-tests.md for motivation):
//!  - all 256 four-sample byte patterns round-trip (LUT injectivity, semantic)
//!  - all-missing SNP (PLINK `01` / PAM `11` special pattern)
//!  - monomorphic all-hom-ref / all-hom-alt (LUT value extremes)
//!  - every `nind % 4` remainder incl. single sample (last-byte padding)
//!  - byte-boundary straddle (13 samples)
//!  - `.bim` placeholder A1=0 flip_02, incl. missing-not-flipped edge
//!  - reader tolerance to non-spec PLINK padding bits
//!  - PGEN and TGENO round-trips (formats the oracle doesn't cover)

use std::path::{Path, PathBuf};
use std::process::Command;

const MISS: u8 = 9;

// ---------------------------------------------------------------------------
// PLINK .bed codec — 2 bits/sample, LSB-first.
// dosage = copies of A1:  00=homA1(2)  01=missing(9)  10=het(1)  11=homA2(0)
// ---------------------------------------------------------------------------
fn dos_to_bed2(d: u8) -> u8 {
    match d {
        2 => 0b00,
        MISS => 0b01,
        1 => 0b10,
        0 => 0b11,
        _ => panic!("bad dosage {d}"),
    }
}
fn bed2_to_dos(two: u8) -> u8 {
    match two & 0b11 {
        0b00 => 2,
        0b01 => MISS,
        0b10 => 1,
        _ => 0,
    }
}

fn encode_bed(rows: &[Vec<u8>], nind: usize) -> Vec<u8> {
    let rec = (nind + 3) / 4;
    let mut out = vec![0x6c, 0x1b, 0x01];
    for r in rows {
        assert_eq!(r.len(), nind);
        let mut bytes = vec![0u8; rec];
        for (i, &d) in r.iter().enumerate() {
            bytes[i / 4] |= dos_to_bed2(d) << (2 * (i % 4));
        }
        out.extend(bytes);
    }
    out
}

fn decode_bed(path: &Path, nind: usize, nsnp: usize) -> Vec<Vec<u8>> {
    let raw = std::fs::read(path).unwrap();
    assert_eq!(&raw[..3], &[0x6c, 0x1b, 0x01], "bad .bed magic");
    let rec = (nind + 3) / 4;
    let body = &raw[3..];
    (0..nsnp)
        .map(|s| {
            let c = &body[s * rec..(s + 1) * rec];
            (0..nind).map(|i| bed2_to_dos(c[i / 4] >> (2 * (i % 4)))).collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// PACKEDANCESTRYMAP .geno decoder — 2 bits/sample MSB-first,
// rlen = max(48, ceil(nind/4)); 00=g0 01=g1 10=g2 11=missing (g = copies of A1)
// ---------------------------------------------------------------------------
fn decode_pam(path: &Path, nind: usize, nsnp: usize) -> Vec<Vec<u8>> {
    let raw = std::fs::read(path).unwrap();
    assert_eq!(&raw[..4], b"GENO", "bad PAM magic");
    let rlen = std::cmp::max(48, (nind + 3) / 4);
    let m = |b: u8| match b & 0b11 {
        0b00 => 0,
        0b01 => 1,
        0b10 => 2,
        _ => MISS,
    };
    (0..nsnp)
        .map(|s| {
            let c = &raw[(s + 1) * rlen..(s + 2) * rlen];
            (0..nind).map(|i| m(c[i / 4] >> (6 - 2 * (i % 4)))).collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Metadata writers + reigen runner
// ---------------------------------------------------------------------------
/// SNP = (id, chrom, bp, a1, a2)
type Snp<'a> = (&'a str, &'a str, u32, char, char);

fn write_bim(path: &Path, snps: &[Snp]) {
    let mut s = String::new();
    for (id, chrom, bp, a1, a2) in snps {
        s.push_str(&format!("{chrom}\t{id}\t0\t{bp}\t{a1}\t{a2}\n"));
    }
    std::fs::write(path, s).unwrap();
}
fn write_fam(path: &Path, nind: usize) {
    let mut s = String::new();
    for i in 0..nind {
        s.push_str(&format!("F{i} I{i} 0 0 1 -9\n"));
    }
    std::fs::write(path, s).unwrap();
}

fn reigen() -> Command {
    Command::new(env!("CARGO_BIN_EXE_reigen"))
}
fn convert(in_geno: &Path, in_snp: &Path, in_ind: &Path, fmt: &str, og: &Path, os: &Path, oi: &Path) {
    let st = reigen()
        .arg("convert")
        .args(["--in-geno"]).arg(in_geno)
        .args(["--in-snp"]).arg(in_snp)
        .args(["--in-ind"]).arg(in_ind)
        .args(["--out-format", fmt])
        .args(["--out-geno"]).arg(og)
        .args(["--out-snp"]).arg(os)
        .args(["--out-ind"]).arg(oi)
        .status()
        .unwrap();
    assert!(st.success(), "reigen convert -> {fmt} failed (in={in_geno:?})");
}

/// Default biallelic, non-placeholder SNPs (no flip): alternating A/G, C/T.
fn default_snps(nsnp: usize) -> Vec<(String, String, u32, char, char)> {
    (0..nsnp)
        .map(|j| {
            let (a1, a2) = if j % 2 == 0 { ('A', 'G') } else { ('C', 'T') };
            (format!("rs{j}"), "1".to_string(), (100 + j) as u32, a1, a2)
        })
        .collect()
}

fn p(dir: &Path, name: &str) -> PathBuf {
    dir.join(name)
}

/// Write a PLINK fileset with the given dosages and SNP metadata.
fn setup_plink(dir: &Path, dos: &[Vec<u8>], snps: &[(String, String, u32, char, char)], nind: usize) {
    std::fs::write(p(dir, "in.bed"), encode_bed(dos, nind)).unwrap();
    let snp_refs: Vec<Snp> = snps
        .iter()
        .map(|(id, c, bp, a1, a2)| (id.as_str(), c.as_str(), *bp, *a1, *a2))
        .collect();
    write_bim(&p(dir, "in.bim"), &snp_refs);
    write_fam(&p(dir, "in.fam"), nind);
}

/// PLINK -> `fmt` (intermediate `mid.<geno_ext>`) -> PLINK, return final dosages.
fn roundtrip(dir: &Path, fmt: &str, geno_ext: &str, snp_ext: &str, ind_ext: &str, nind: usize, nsnp: usize) -> Vec<Vec<u8>> {
    let mg = p(dir, &format!("mid.{geno_ext}"));
    let ms = p(dir, &format!("mid.{snp_ext}"));
    let mi = p(dir, &format!("mid.{ind_ext}"));
    convert(&p(dir, "in.bed"), &p(dir, "in.bim"), &p(dir, "in.fam"), fmt, &mg, &ms, &mi);
    convert(&mg, &ms, &mi, "PACKEDPED", &p(dir, "out.bed"), &p(dir, "out.bim"), &p(dir, "out.fam"));
    decode_bed(&p(dir, "out.bed"), nind, nsnp)
}

/// Run a PLINK->fmt->PLINK round trip for all formats and assert dosages survive.
fn assert_roundtrips_all(dos: &[Vec<u8>], nind: usize) {
    let nsnp = dos.len();
    let snps = default_snps(nsnp);
    for (fmt, ge, se, ie) in [
        ("PACKEDANCESTRYMAP", "geno", "snp", "ind"),
        ("EIGENSTRAT", "geno", "snp", "ind"),
        ("TGENO", "tgeno", "snp", "ind"),
        ("PLINK2", "pgen", "pvar", "psam"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        setup_plink(dir.path(), dos, &snps, nind);
        let got = roundtrip(dir.path(), fmt, ge, se, ie, nind, nsnp);
        assert_eq!(got, dos, "round trip via {fmt} corrupted dosages (nind={nind})");
    }
}

// ===========================================================================
// 1. All 256 four-sample byte patterns
// ===========================================================================
#[test]
fn all_256_four_sample_patterns_roundtrip() {
    let vals = [0u8, 1, 2, MISS];
    let mut dos = Vec::with_capacity(256);
    for a in vals {
        for b in vals {
            for c in vals {
                for d in vals {
                    dos.push(vec![a, b, c, d]);
                }
            }
        }
    }
    assert_eq!(dos.len(), 256);
    assert_roundtrips_all(&dos, 4);
}

// ===========================================================================
// 2. All-missing SNP (PLINK 0x55 / PAM 0xff pattern)
// ===========================================================================
#[test]
fn all_missing_snp_roundtrips() {
    for nind in [5usize, 8] {
        assert_roundtrips_all(&[vec![MISS; nind], vec![2; nind]], nind);
    }
}

// Also verify the PAM bytes directly: all-missing must decode to all-9, not all-0.
#[test]
fn all_missing_snp_pam_bytes_are_missing_not_zero() {
    let nind = 8;
    let dir = tempfile::tempdir().unwrap();
    let snps = default_snps(1);
    setup_plink(dir.path(), &[vec![MISS; nind]], &snps, nind);
    convert(&p(dir.path(), "in.bed"), &p(dir.path(), "in.bim"), &p(dir.path(), "in.fam"),
            "PACKEDANCESTRYMAP", &p(dir.path(), "m.geno"), &p(dir.path(), "m.snp"), &p(dir.path(), "m.ind"));
    let g = decode_pam(&p(dir.path(), "m.geno"), nind, 1);
    assert_eq!(g[0], vec![MISS; nind]);
}

// ===========================================================================
// 3. Monomorphic extremes (all hom-ref g=2, all hom-alt g=0)
// ===========================================================================
#[test]
fn monomorphic_all_hom_ref_and_alt_roundtrip() {
    let nind = 8;
    assert_roundtrips_all(&[vec![2; nind], vec![0; nind]], nind);
}

// ===========================================================================
// 4. Every nind % 4 remainder, including single sample
// ===========================================================================
fn cycling(nind: usize, nsnp: usize) -> Vec<Vec<u8>> {
    let vals = [0u8, 1, 2, MISS];
    (0..nsnp)
        .map(|s| (0..nind).map(|i| vals[(i + s) % 4]).collect())
        .collect()
}

macro_rules! padding_test {
    ($name:ident, $nind:expr) => {
        #[test]
        fn $name() {
            assert_roundtrips_all(&cycling($nind, 3), $nind);
        }
    };
}
padding_test!(padding_nind_1, 1); // single sample: 48-byte PAM header + 1-byte record
padding_test!(padding_nind_2, 2);
padding_test!(padding_nind_3, 3);
padding_test!(padding_nind_5, 5);
padding_test!(padding_nind_6, 6);
padding_test!(padding_nind_7, 7);
padding_test!(padding_nind_9, 9);

// ===========================================================================
// 5. Byte-boundary straddle (13 samples = 3 full bytes + 1 sample)
// ===========================================================================
#[test]
fn straddle_13_samples_roundtrip() {
    assert_roundtrips_all(&cycling(13, 6), 13);
}

// ===========================================================================
// 6. Placeholder A1=0 (flip_02), including missing-not-flipped
// ===========================================================================
#[test]
fn placeholder_a1_zero_flips_genotypes_but_not_missing() {
    // .bim: A1 = '0' placeholder, A2 = 'A' real. PLINK stores hom-A2 / missing.
    // dosage-of-A1 here: 0 = hom-A2 (two copies of the real allele A), 9 = missing.
    let nind = 4;
    let dir = tempfile::tempdir().unwrap();
    let dos = vec![vec![0u8, 0, MISS, 0]]; // 3x hom-A2, 1 missing
    std::fs::write(p(dir.path(), "in.bed"), encode_bed(&dos, nind)).unwrap();
    write_bim(&p(dir.path(), "in.bim"), &[("rsP", "1", 100, '0', 'A')]);
    write_fam(&p(dir.path(), "in.fam"), nind);

    convert(&p(dir.path(), "in.bed"), &p(dir.path(), "in.bim"), &p(dir.path(), "in.fam"),
            "PACKEDANCESTRYMAP", &p(dir.path(), "m.geno"), &p(dir.path(), "m.snp"), &p(dir.path(), "m.ind"));

    // convertf convention: real allele A becomes allele1, genotypes count A.
    // hom-A2 (two copies of A) -> g2; missing stays missing (NOT flipped to g0).
    let g = decode_pam(&p(dir.path(), "m.geno"), nind, 1);
    assert_eq!(g[0], vec![2u8, 2, MISS, 2], "flip_02 wrong (missing must not flip)");

    // .snp must show A as allele1 and X (placeholder) as allele2.
    let snp = std::fs::read_to_string(p(dir.path(), "m.snp")).unwrap();
    let cols: Vec<&str> = snp.split_whitespace().collect();
    assert_eq!((cols[4], cols[5]), ("A", "X"), ".snp alleles not normalized: {snp:?}");
}

// ===========================================================================
// 7. Reader tolerance to non-spec PLINK padding bits
// ===========================================================================
#[test]
fn reader_tolerates_corrupted_padding_bits() {
    // nind=5 -> 2 bytes/record. byte0 = samples 0-3, byte1 bits1-0 = sample4,
    // byte1 bits7-2 = padding for non-existent samples 5/6/7. Set them to 1s.
    let nind = 5;
    let dir = tempfile::tempdir().unwrap();
    // valid samples: [hom-A1(2), het(1), hom-A2(0), miss(9), hom-A1(2)]
    // byte0: s0=00 s1=10 s2=11 s3=01  -> 0b01_11_10_00 = 0x78
    // byte1: s4=00 in bits1-0, padding bits7-2 all 1 -> 0b1111_11_00 = 0xFC
    let bed = vec![0x6c, 0x1b, 0x01, 0x78, 0xFC];
    std::fs::write(p(dir.path(), "in.bed"), bed).unwrap();
    write_bim(&p(dir.path(), "in.bim"), &[("rsC", "1", 100, 'A', 'G')]);
    write_fam(&p(dir.path(), "in.fam"), nind);

    convert(&p(dir.path(), "in.bed"), &p(dir.path(), "in.bim"), &p(dir.path(), "in.fam"),
            "PACKEDANCESTRYMAP", &p(dir.path(), "m.geno"), &p(dir.path(), "m.snp"), &p(dir.path(), "m.ind"));
    let g = decode_pam(&p(dir.path(), "m.geno"), nind, 1);
    assert_eq!(g[0], vec![2u8, 1, 0, MISS, 2], "corrupted padding leaked into valid samples");
}

// ===========================================================================
// 8. PGEN / TGENO get explicit single-format round trips too (small nsnp,
//    which is exactly where the TGENO record-length bug used to bite).
// ===========================================================================
#[test]
fn pgen_roundtrip_mixed_patterns_small() {
    let nind = 6;
    let dos = vec![vec![2u8, 0, 1, MISS, 2, 0], vec![0, 0, 0, 0, 0, 0], vec![MISS; 6]];
    let nsnp = dos.len();
    let snps = default_snps(nsnp);
    let dir = tempfile::tempdir().unwrap();
    setup_plink(dir.path(), &dos, &snps, nind);
    let got = roundtrip(dir.path(), "PLINK2", "pgen", "pvar", "psam", nind, nsnp);
    assert_eq!(got, dos);
}

#[test]
fn tgeno_roundtrip_small_nsnp() {
    // nsnp=3 -> TGENO rlen=1; regression guard for the dropped 48-byte floor.
    let nind = 3;
    let dos = vec![vec![2u8, 1, 0], vec![MISS, 2, 0], vec![1, 1, 1]];
    let nsnp = dos.len();
    let snps = default_snps(nsnp);
    let dir = tempfile::tempdir().unwrap();
    setup_plink(dir.path(), &dos, &snps, nind);
    let got = roundtrip(dir.path(), "TGENO", "tgeno", "snp", "ind", nind, nsnp);
    assert_eq!(got, dos);
}
