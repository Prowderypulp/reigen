//! PLINK `.bim` format.
//!
//! Six whitespace-separated columns — **cols 1 and 2 swapped relative to
//! `.snp`**:
//!
//! ```text
//! chrom  snp_id  gen_pos  phys_pos  a1  a2
//! ```
//!
//! # Chrom codes
//!
//! PLINK uses numeric codes: 1..22 autosomes, 23 = X, 24 = Y, 25 = XY
//! (pseudo-autosomal), 26 = MT. Literal `X`/`Y`/`MT`/`XY` also accepted.
//! Internal mapping used by `.snp` is X=numchrom+1, Y=numchrom+2,
//! MT=numchrom+3, XY=numchrom+4. PLINK numeric codes use XY before MT
//! (25=XY, 26=MT for humans), so numeric 25/26 are remapped at the I/O
//! boundary to keep internal representation consistent.
//!
//! # A1 / A2 semantics
//!
//! PLINK `.bim` col 5 = A1, col 6 = A2. The PLINK `.bed` genotype counts
//! copies of A1 (code `00` = homozygous A1). AdmixTools/EIGENSTRAT stores the
//! same relationship: `.snp` col 5 = allele1 and the `.geno` value counts
//! copies of allele1 (`2` = homozygous allele1). The two formats use the
//! *same* column order, so there is **no allele swap** across this boundary:
//!
//! | column |     .snp     |     .bim     |
//! |--------|--------------|--------------|
//! |    5   |  allele1     | allele1 (A1) |
//! |    6   |  allele2     | allele2 (A2) |
//!
//! Reading `.bim`: `bim[5]` → `SnpRow.allele1`, `bim[6]` → `SnpRow.allele2`.
//! Writing `.bim` is the identity. This matches `convertf` byte-for-byte
//! (verified: convertf EIGENSTRAT→PACKEDPED leaves col 5 in A1, and the
//! `.bed` for `geno=2` is `00`=homozygous A1). The genotype 0↔2 polarity is
//! carried entirely by the `.bed` codec LUT in `geno/packed_ped.rs`; do not
//! re-introduce a swap here.

use super::{split_lines, SnpRow};
use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn read(path: &Path, numchrom: u32) -> Result<Vec<SnpRow>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    if file.metadata()?.len() == 0 {
        return Ok(Vec::new());
    }
    let mmap = unsafe { Mmap::map(&file) }.with_context(|| format!("mmap {}", path.display()))?;

    let mut rows = Vec::new();
    for (lineno, line) in split_lines(&mmap).enumerate() {
        let lineno = lineno + 1;
        if line.iter().all(|&b| b.is_ascii_whitespace()) {
            continue;
        }
        let row = parse_bim_line(line, numchrom)
            .with_context(|| format!("{}:{}", path.display(), lineno))?;
        rows.push(row);
    }
    Ok(rows)
}

pub(crate) fn read_flip_02_mask(path: &Path) -> Result<Vec<bool>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    if file.metadata()?.len() == 0 {
        return Ok(Vec::new());
    }
    let mmap = unsafe { Mmap::map(&file) }.with_context(|| format!("mmap {}", path.display()))?;

    let mut mask = Vec::new();
    for (lineno, line) in split_lines(&mmap).enumerate() {
        let lineno = lineno + 1;
        if line.iter().all(|&b| b.is_ascii_whitespace()) {
            continue;
        }
        let flip = parse_bim_flip_02(line).with_context(|| format!("{}:{}", path.display(), lineno))?;
        mask.push(flip);
    }
    Ok(mask)
}

pub fn write(path: &Path, rows: &[SnpRow], numchrom: u32) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(256 * 1024, file);
    for r in rows {
        // PLINK emits chrom as numeric; we follow that convention. Tools
        // that want literal X/Y/MT can post-process.
        // No swap: bim col 5 (A1) = allele1, bim col 6 (A2) = allele2.
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom_to_plink_numeric(r.chrom, numchrom),
            r.id,
            r.genetic_pos,
            r.physical_pos,
            r.allele1 as char,
            r.allele2 as char,
        )?;
    }
    w.flush()?;
    Ok(())
}

fn parse_bim_line(line: &[u8], numchrom: u32) -> Result<SnpRow> {
    let mut cols = line
        .split(|b: &u8| b.is_ascii_whitespace())
        .filter(|c| !c.is_empty());

    let chrom_raw = cols.next().ok_or_else(|| anyhow!("missing chrom"))?;
    let id = cols.next().ok_or_else(|| anyhow!("missing snp id"))?;
    let gen = cols.next().ok_or_else(|| anyhow!("missing genetic pos"))?;
    let phys = cols.next().ok_or_else(|| anyhow!("missing physical pos"))?;
    let a1 = cols.next().ok_or_else(|| anyhow!("missing a1"))?;
    let a2 = cols.next().ok_or_else(|| anyhow!("missing a2"))?;

    let chrom = parse_chrom(chrom_raw, numchrom)?;
    let id = std::str::from_utf8(id)?.to_owned();
    let genetic_pos: f64 = std::str::from_utf8(gen)?
        .parse()
        .map_err(|e| anyhow!("bad genetic_pos: {e}"))?;
    let physical_pos: u64 = std::str::from_utf8(phys)?
        .parse()
        .map_err(|e| anyhow!("bad physical_pos: {e}"))?;

    if a1.len() != 1 || a2.len() != 1 {
        bail!("multi-char allele in .bim not supported");
    }
    let (allele1, allele2, _flip_02) = normalize_bim_alleles(a1[0], a2[0]);
    Ok(SnpRow {
        id,
        chrom,
        genetic_pos,
        physical_pos,
        allele1,
        allele2,
    })
}

fn parse_bim_flip_02(line: &[u8]) -> Result<bool> {
    let mut cols = line
        .split(|b: &u8| b.is_ascii_whitespace())
        .filter(|c| !c.is_empty());
    let _chrom = cols.next().ok_or_else(|| anyhow!("missing chrom"))?;
    let _id = cols.next().ok_or_else(|| anyhow!("missing snp id"))?;
    let _gen = cols.next().ok_or_else(|| anyhow!("missing genetic pos"))?;
    let _phys = cols.next().ok_or_else(|| anyhow!("missing physical pos"))?;
    let a1 = cols.next().ok_or_else(|| anyhow!("missing a1"))?;
    let a2 = cols.next().ok_or_else(|| anyhow!("missing a2"))?;
    if a1.len() != 1 || a2.len() != 1 {
        bail!("multi-char allele in .bim not supported");
    }
    Ok(matches!((a1[0], a2[0]), (b'0', real) if real != b'0'))
}

fn normalize_bim_alleles(a1: u8, a2: u8) -> (u8, u8, bool) {
    match (a1, a2) {
        (b'0', b'0') => (b'X', b'X', false),
        (b'0', real) => (real, b'X', true),
        (real, b'0') => (real, b'X', false),
        _ => (a1, a2, false),
    }
}

/// Parse a PLINK-convention chromosome token (numeric or `X`/`Y`/`MT`/`XY`)
/// into reigen's internal code. Shared with the PVAR reader so `.bim` and
/// `.pvar` agree on the XY/MT remap (see the note below). `pub(crate)`.
pub(crate) fn parse_chrom(raw: &[u8], numchrom: u32) -> Result<u8> {
    let s = std::str::from_utf8(raw)?;
    let up = s.to_ascii_uppercase();
    let v_plink: u32 = match up.as_str() {
        "X" => numchrom + 1,
        "Y" => numchrom + 2,
        // Assign the PLINK *numeric* slot for each token; the remap below then
        // converts to reigen's internal code. PLINK numbering is XY=25, MT=26.
        "XY" => numchrom + 3,
        "MT" | "M" => numchrom + 4,
        num => num.parse().map_err(|e| anyhow!("bad chrom: {e}"))?,
    };
    // PLINK numeric code order differs from our internal mapping:
    // 23=X, 24=Y, 25=XY, 26=MT (PLINK) vs 23=X, 24=Y, 25=MT, 26=XY (internal).
    let v = if v_plink == numchrom + 3 {
        numchrom + 4
    } else if v_plink == numchrom + 4 {
        numchrom + 3
    } else {
        v_plink
    };
    if v > u8::MAX as u32 {
        bail!("chrom {v} out of u8 range");
    }
    Ok(v as u8)
}

/// Convert reigen's internal chromosome code to a PLINK numeric code (applies
/// the XY/MT remap). Shared with the PVAR writer. `pub(crate)`.
pub(crate) fn chrom_to_plink_numeric(chrom_internal: u8, numchrom: u32) -> u32 {
    let c = chrom_internal as u32;
    if c == numchrom + 3 {
        numchrom + 4
    } else if c == numchrom + 4 {
        numchrom + 3
    } else {
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(s: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(s.as_bytes()).unwrap();
        f
    }

    #[test]
    fn reads_bim_no_swap() {
        // PLINK .bim: chrom id gen phys a1 a2
        let f = write_tmp("1\trs1\t0.001\t752566\tG\tA\n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].id, "rs1");
        assert_eq!(rows[0].chrom, 1);
        // No swap: PLINK A1=G → allele1, PLINK A2=A → allele2.
        assert_eq!(rows[0].allele1, b'G');
        assert_eq!(rows[0].allele2, b'A');
    }

    #[test]
    fn handles_plink_numeric_x() {
        // PLINK often uses numeric 23 for X
        let f = write_tmp("23\trsX\t0.0\t1000\tC\tT\n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows[0].chrom, 23);
    }

    #[test]
    fn handles_literal_x() {
        let f = write_tmp("X\trsX\t0.0\t1000\tA\tC\n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows[0].chrom, 23);
    }

    #[test]
    fn maps_plink_numeric_xy_and_mt_to_internal_order() {
        // PLINK numeric: 25=XY, 26=MT.
        let f_xy = write_tmp("25\trsXY\t0.0\t1000\tA\tC\n");
        let rows_xy = read(f_xy.path(), 22).unwrap();
        assert_eq!(rows_xy[0].chrom, 26); // internal XY

        let f_mt = write_tmp("26\trsMT\t0.0\t1000\tA\tC\n");
        let rows_mt = read(f_mt.path(), 22).unwrap();
        assert_eq!(rows_mt[0].chrom, 25); // internal MT
    }

    #[test]
    fn write_then_read_roundtrip_preserves_snp_allele_order() {
        // Start from a SnpRow in AdmixTools convention. Write as .bim, read
        // back. Allele order must match after double-swap.
        let row = SnpRow {
            id: "rs1".into(),
            chrom: 1,
            genetic_pos: 0.001,
            physical_pos: 752566,
            allele1: b'A',
            allele2: b'G',
        };
        let f = tempfile::NamedTempFile::new().unwrap();
        write(f.path(), &[row.clone()], 22).unwrap();
        let got = read(f.path(), 22).unwrap();
        assert_eq!(got[0].allele1, b'A');
        assert_eq!(got[0].allele2, b'G');
    }

    #[test]
    fn bim_write_content() {
        // Spot-check the on-disk text layout.
        let row = SnpRow {
            id: "rs1".into(),
            chrom: 1,
            genetic_pos: 0.001,
            physical_pos: 752566,
            allele1: b'A',
            allele2: b'G',
        };
        let f = tempfile::NamedTempFile::new().unwrap();
        write(f.path(), &[row], 22).unwrap();
        let text = std::fs::read_to_string(f.path()).unwrap();
        // No swap: bim col 5 (A1) = allele1 ('A'); col 6 (A2) = allele2 ('G').
        assert!(text.contains("\tA\tG\n"), "got: {text:?}");
        assert!(text.starts_with("1\trs1"));
    }

    #[test]
    fn write_emits_plink_numeric_xy_mt_order() {
        let rows = vec![
            SnpRow {
                id: "rsMT".into(),
                chrom: 25, // internal MT
                genetic_pos: 0.0,
                physical_pos: 1,
                allele1: b'A',
                allele2: b'C',
            },
            SnpRow {
                id: "rsXY".into(),
                chrom: 26, // internal XY
                genetic_pos: 0.0,
                physical_pos: 2,
                allele1: b'G',
                allele2: b'T',
            },
        ];
        let f = tempfile::NamedTempFile::new().unwrap();
        write(f.path(), &rows, 22).unwrap();
        let text = std::fs::read_to_string(f.path()).unwrap();
        let mut lines = text.lines();
        assert!(lines.next().unwrap().starts_with("26\trsMT\t")); // MT -> 26 in PLINK
        assert!(lines.next().unwrap().starts_with("25\trsXY\t")); // XY -> 25 in PLINK
    }

    #[test]
    fn parse_chrom_string_tokens_match_plink_numeric_codes() {
        // String X/Y/XY/MT must map to the SAME internal code as the numeric
        // PLINK codes 23/24/25/26 (XY=25, MT=26 in PLINK numbering).
        for (s, num) in [("X", "23"), ("Y", "24"), ("XY", "25"), ("MT", "26"), ("M", "26")] {
            assert_eq!(
                parse_chrom(s.as_bytes(), 22).unwrap(),
                parse_chrom(num.as_bytes(), 22).unwrap(),
                "string {s} must match numeric {num}"
            );
        }
        // And a full round trip back to PLINK numeric is order-preserving.
        assert_eq!(chrom_to_plink_numeric(parse_chrom(b"XY", 22).unwrap(), 22), 25);
        assert_eq!(chrom_to_plink_numeric(parse_chrom(b"MT", 22).unwrap(), 22), 26);
    }

    #[test]
    fn normalizes_zero_placeholder_alleles_like_convertf() {
        let f = write_tmp("1\trs1\t0.0\t100\t0\tA\n2\trs2\t0.0\t200\tA\t0\n2\trs3\t0.0\t300\t0\t0\n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!((rows[0].allele1, rows[0].allele2), (b'A', b'X'));
        assert_eq!((rows[1].allele1, rows[1].allele2), (b'A', b'X'));
        assert_eq!((rows[2].allele1, rows[2].allele2), (b'X', b'X'));

        let mask = read_flip_02_mask(f.path()).unwrap();
        assert_eq!(mask, vec![true, false, false]);
    }
}
