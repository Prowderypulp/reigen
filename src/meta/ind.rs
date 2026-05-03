//! AdmixTools `.ind` format — three whitespace-separated columns:
//!
//! ```text
//! sample_id  sex  population_or_status
//! ```
//!
//! - Sample ID: upstream caps at 39 chars; we warn on longer but accept.
//! - Sex: M / F / U (case-insensitive). Anything else → Unknown.
//! - Population: free-form string. Literal `Ignore` (any case) means this
//!   sample is dropped from all convertf outputs — we record `ignore=true`
//!   and let `filter::IndFilter` enforce the drop.

use super::{split_lines, IndRow, Sex};
use anyhow::{anyhow, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const IND_ID_WARN_LEN: usize = 39;

pub fn read(path: &Path) -> Result<Vec<IndRow>> {
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

        let row = parse_ind_line(line).with_context(|| format!("{}:{}", path.display(), lineno))?;
        if row.id.len() > IND_ID_WARN_LEN {
            log::warn!(
                "sample id {:?} exceeds {IND_ID_WARN_LEN} chars — \
                        downstream AdmixTools tools may truncate",
                row.id
            );
        }
        rows.push(row);
    }
    Ok(rows)
}

pub fn write(path: &Path, rows: &[IndRow]) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(64 * 1024, file);
    for r in rows {
        // Upstream: "%20s %c %s\n" — right-aligned id, sex char, pop.
        writeln!(w, "{:>20} {} {}", r.id, r.sex.as_ind_char(), r.pop)?;
    }
    w.flush()?;
    Ok(())
}

fn parse_ind_line(line: &[u8]) -> Result<IndRow> {
    let mut cols = line
        .split(|b: &u8| b.is_ascii_whitespace())
        .filter(|c| !c.is_empty());

    let id = cols.next().ok_or_else(|| anyhow!("missing sample id"))?;
    let sex_raw = cols.next().ok_or_else(|| anyhow!("missing sex"))?;
    let pop = cols.next().ok_or_else(|| anyhow!("missing population"))?;

    let id = std::str::from_utf8(id)?.to_owned();
    let pop = std::str::from_utf8(pop)?.to_owned();
    let sex = Sex::from_char(*sex_raw.first().unwrap_or(&b'U'));
    let ignore = pop.eq_ignore_ascii_case("Ignore");

    Ok(IndRow {
        id,
        sex,
        pop,
        ignore,
    })
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
    fn reads_three_columns() {
        let f = write_tmp(
            "S_French-1.DG M French.DG\n\
             BR_Onge-2.DG  F Onge.DG\n\
             some_ignored  U Ignore\n",
        );
        let rows = read(f.path()).unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].id, "S_French-1.DG");
        assert_eq!(rows[0].sex, Sex::Male);
        assert_eq!(rows[0].pop, "French.DG");
        assert!(!rows[0].ignore);
        assert_eq!(rows[1].sex, Sex::Female);
        assert!(rows[2].ignore);
    }

    #[test]
    fn case_insensitive_sex_and_ignore() {
        let f = write_tmp("a m Pop\nb f Pop\nc u iGnOrE\n");
        let rows = read(f.path()).unwrap();
        assert_eq!(rows[0].sex, Sex::Male);
        assert_eq!(rows[1].sex, Sex::Female);
        assert!(rows[2].ignore);
    }

    #[test]
    fn roundtrip() {
        let f = write_tmp("S1 M Pop1\nS2 F Pop2\n");
        let rows = read(f.path()).unwrap();

        let out = tempfile::NamedTempFile::new().unwrap();
        write(out.path(), &rows).unwrap();

        let rows2 = read(out.path()).unwrap();
        assert_eq!(rows.len(), rows2.len());
        for (a, b) in rows.iter().zip(rows2.iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.sex, b.sex);
            assert_eq!(a.pop, b.pop);
        }
    }
}
