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
//! In PLINK `.bim`:
//! - col 5 = A1 (usually minor allele; PLINK historically put the "1"
//!   allele here — matches the "variant" allele in AdmixTools)
//! - col 6 = A2 (usually major; matches AdmixTools "reference" allele)
//!
//! Our `SnpRow` holds `allele1` = AdmixTools allele1 = EIGENSTRAT "reference"
//! = `.snp` col 5 = **A2** in PLINK. The swap happens here at the I/O
//! boundary:
//!
//! | column |     .snp     |     .bim     |
//! |--------|--------------|--------------|
//! |    5   |  allele1 (0) | allele2 (A1) |
//! |    6   |  allele2 (2) | allele1 (A2) |
//!
//! Translated: reading `.bim`, `bim[5]` → `SnpRow.allele2` and
//! `bim[6]` → `SnpRow.allele1`. Writing `.bim` does the same swap.

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

pub fn write(path: &Path, rows: &[SnpRow], numchrom: u32) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(256 * 1024, file);
    for r in rows {
        // PLINK emits chrom as numeric; we follow that convention. Tools
        // that want literal X/Y/MT can post-process.
        // Swap: bim col 5 = allele2, bim col 6 = allele1 (see module doc).
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom_to_plink_numeric(r.chrom, numchrom),
            r.id,
            r.genetic_pos,
            r.physical_pos,
            r.allele2 as char,
            r.allele1 as char,
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
    // Swap at boundary: PLINK A1 (col 5) → SnpRow.allele2.
    Ok(SnpRow {
        id,
        chrom,
        genetic_pos,
        physical_pos,
        allele1: a2[0], // PLINK A2 → AdmixTools allele1
        allele2: a1[0], // PLINK A1 → AdmixTools allele2
    })
}

fn parse_chrom(raw: &[u8], numchrom: u32) -> Result<u8> {
    let s = std::str::from_utf8(raw)?;
    let up = s.to_ascii_uppercase();
    let v_plink: u32 = match up.as_str() {
        "X" => numchrom + 1,
        "Y" => numchrom + 2,
        "MT" | "M" => numchrom + 3,
        "XY" => numchrom + 4,
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

fn chrom_to_plink_numeric(chrom_internal: u8, numchrom: u32) -> u32 {
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
    fn reads_bim_swaps_alleles() {
        // PLINK .bim: chrom id gen phys a1 a2
        let f = write_tmp("1\trs1\t0.001\t752566\tG\tA\n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].id, "rs1");
        assert_eq!(rows[0].chrom, 1);
        // Swap: PLINK A1=G → allele2, PLINK A2=A → allele1
        assert_eq!(rows[0].allele1, b'A');
        assert_eq!(rows[0].allele2, b'G');
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
        // bim col 5 (A1) should be allele2 ('G'); col 6 (A2) should be allele1 ('A').
        assert!(text.contains("\tG\tA\n"), "got: {text:?}");
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
}
