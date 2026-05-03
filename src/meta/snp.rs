//! AdmixTools `.snp` format.
//!
//! Six whitespace-separated columns (last two optional but required for any
//! allele-aware format, which is all we care about):
//!
//! ```text
//! snp_id  chrom  gen_pos  phys_pos  allele1  allele2
//! ```
//!
//! Chrom: numeric OR literal `X` / `Y` / `MT` / `XY`. We map to u8 using
//! `numchrom`-based offsets (X = numchrom+1, Y = numchrom+2, MT = numchrom+3).
//!
//! # Reader
//!
//! `mmap` + `memchr(b'\n')` for line splits. Per line we scan whitespace
//! with a tiny state machine to avoid `str::split` overhead. Only the six
//! columns are parsed; extras ignored.
//!
//! # Writer
//!
//! Fixed-width columns matching upstream `outsnp` in `mcio.c`.

use super::SnpRow;
use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Parse a `.snp` file. Returns rows in file order.
pub fn read(path: &Path, numchrom: u32) -> Result<Vec<SnpRow>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    if file.metadata()?.len() == 0 {
        return Ok(Vec::new());
    }
    // SAFETY: mmap of a regular file; fine for batch CLI.
    let mmap = unsafe { Mmap::map(&file) }.with_context(|| format!("mmap {}", path.display()))?;

    let bytes: &[u8] = &mmap;
    let mut rows = Vec::new();

    for (lineno, line) in split_lines(bytes).enumerate() {
        let lineno = lineno + 1;
        if line.iter().all(|&b| b.is_ascii_whitespace()) {
            continue;
        }

        let row = parse_snp_line(line, numchrom)
            .with_context(|| format!("{}:{}", path.display(), lineno))?;
        rows.push(row);
    }

    Ok(rows)
}

/// Write a `.snp` file in upstream convertf layout.
pub fn write(path: &Path, rows: &[SnpRow], numchrom: u32) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(256 * 1024, file);

    for row in rows {
        let chrom_str = chrom_to_str(row.chrom, numchrom);
        writeln!(
            w,
            "{:<20} {:>2} {:>12.6} {:>12} {} {}",
            row.id,
            chrom_str,
            row.genetic_pos,
            row.physical_pos,
            row.allele1 as char,
            row.allele2 as char,
        )?;
    }

    w.flush()?;
    Ok(())
}

// ------------------------------------------------------------------
// internals
// ------------------------------------------------------------------

fn split_lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut start = 0usize;
    let len = bytes.len();
    std::iter::from_fn(move || {
        if start >= len {
            return None;
        }
        let rest = &bytes[start..];
        match memchr::memchr(b'\n', rest) {
            Some(off) => {
                let mut end = start + off;
                let line_start = start;
                start = end + 1;
                if end > line_start && bytes[end - 1] == b'\r' {
                    end -= 1;
                }
                Some(&bytes[line_start..end])
            }
            None => {
                let line_start = start;
                start = len;
                Some(&bytes[line_start..len])
            }
        }
    })
}

fn parse_snp_line(line: &[u8], numchrom: u32) -> Result<SnpRow> {
    let mut it = ByteCols::new(line);
    let id = it.next_required("snp id")?;
    let chrom_raw = it.next_required("chrom")?;
    let gen_pos = it.next_required("genetic pos")?;
    let phys_pos = it.next_required("physical pos")?;
    let a1 = it.next_optional();
    let a2 = it.next_optional();

    let chrom = parse_chrom(chrom_raw, numchrom).with_context(|| {
        format!(
            "chrom field {:?}",
            std::str::from_utf8(chrom_raw).unwrap_or("<non-utf8>")
        )
    })?;

    let gen_pos: f64 = std::str::from_utf8(gen_pos)?
        .parse()
        .map_err(|e| anyhow!("bad genetic_pos: {e}"))?;
    let phys_pos: u64 = std::str::from_utf8(phys_pos)?
        .parse()
        .map_err(|e| anyhow!("bad physical_pos: {e}"))?;

    let (allele1, allele2) = match (a1, a2) {
        (Some(a), Some(b)) if a.len() == 1 && b.len() == 1 => (a[0], b[0]),
        (None, None) => (b'X', b'X'),
        (Some(a), Some(b)) => bail!(
            "multi-char allele not supported: {:?} / {:?}",
            std::str::from_utf8(a).unwrap_or("<non-utf8>"),
            std::str::from_utf8(b).unwrap_or("<non-utf8>"),
        ),
        _ => bail!("expected 4 or 6 whitespace columns"),
    };

    Ok(SnpRow {
        id: std::str::from_utf8(id)?.to_owned(),
        chrom,
        genetic_pos: gen_pos,
        physical_pos: phys_pos,
        allele1,
        allele2,
    })
}

struct ByteCols<'a> {
    rest: &'a [u8],
}
impl<'a> ByteCols<'a> {
    fn new(line: &'a [u8]) -> Self {
        Self { rest: line }
    }

    fn advance_past_ws(&mut self) {
        while let Some((&b, r)) = self.rest.split_first() {
            if b.is_ascii_whitespace() {
                self.rest = r;
            } else {
                break;
            }
        }
    }

    fn next_optional(&mut self) -> Option<&'a [u8]> {
        self.advance_past_ws();
        if self.rest.is_empty() {
            return None;
        }
        let end = self
            .rest
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .unwrap_or(self.rest.len());
        let (col, r) = self.rest.split_at(end);
        self.rest = r;
        Some(col)
    }

    fn next_required(&mut self, what: &str) -> Result<&'a [u8]> {
        self.next_optional()
            .ok_or_else(|| anyhow!("missing column: {what}"))
    }
}

fn parse_chrom(raw: &[u8], numchrom: u32) -> Result<u8> {
    let s = std::str::from_utf8(raw)?;
    let up = s.to_ascii_uppercase();
    let v: u32 = match up.as_str() {
        "X" => numchrom + 1,
        "Y" => numchrom + 2,
        "MT" | "M" => numchrom + 3,
        "XY" => numchrom + 4,
        num => num.parse().map_err(|e| anyhow!("bad chrom: {e}"))?,
    };
    if v > u8::MAX as u32 {
        bail!("chrom {v} out of u8 range");
    }
    Ok(v as u8)
}

fn chrom_to_str(c: u8, numchrom: u32) -> String {
    let c32 = c as u32;
    if c32 == numchrom + 1 {
        "X".into()
    } else if c32 == numchrom + 2 {
        "Y".into()
    } else if c32 == numchrom + 3 {
        "MT".into()
    } else if c32 == numchrom + 4 {
        "XY".into()
    } else {
        c.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(contents: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        f
    }

    #[test]
    fn reads_minimal_snp() {
        let f = write_tmp(
            "rs1 1 0.001 752566 A G\n\
             rs2 2 0.0022 832918 C T\n",
        );
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id, "rs1");
        assert_eq!(rows[0].chrom, 1);
        assert_eq!(rows[0].physical_pos, 752566);
        assert_eq!(rows[0].allele1, b'A');
        assert_eq!(rows[1].chrom, 2);
    }

    #[test]
    fn handles_x_chrom() {
        let f = write_tmp("rsX X 0.0 1000 A T\n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows[0].chrom, 23);
    }

    #[test]
    fn handles_mt() {
        let f = write_tmp("rsMT MT 0.0 1 A T\n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows[0].chrom, 25);
    }

    #[test]
    fn tolerates_extra_whitespace() {
        let f = write_tmp("  rs1   1   0.001   752566    A   G  \n");
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows[0].id, "rs1");
        assert_eq!(rows[0].allele2, b'G');
    }

    #[test]
    fn skips_blank_lines() {
        let f = write_tmp(
            "\n\n\
             rs1 1 0.001 752566 A G\n\
             \n\
             rs2 1 0.002 800000 C T\n",
        );
        let rows = read(f.path(), 22).unwrap();
        assert_eq!(rows.len(), 2);
    }

    #[test]
    fn roundtrip_write_read() {
        let f = write_tmp("rs1 1 0.001 752566 A G\nrsX X 0.5 5000000 C T\n");
        let rows = read(f.path(), 22).unwrap();

        let out = tempfile::NamedTempFile::new().unwrap();
        write(out.path(), &rows, 22).unwrap();

        let rows2 = read(out.path(), 22).unwrap();
        assert_eq!(rows.len(), rows2.len());
        for (a, b) in rows.iter().zip(rows2.iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.chrom, b.chrom);
            assert_eq!(a.physical_pos, b.physical_pos);
            assert_eq!(a.allele1, b.allele1);
            assert_eq!(a.allele2, b.allele2);
            assert!((a.genetic_pos - b.genetic_pos).abs() < 1e-9);
        }
    }
}
