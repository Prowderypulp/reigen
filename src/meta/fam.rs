//! PLINK `.fam` format.
//!
//! Six whitespace-separated columns:
//!
//! ```text
//! FID  IID  PID  MID  sex  phenotype
//! ```
//!
//! - FID: family/population ID → maps to `IndRow.pop`
//! - IID: individual ID → maps to `IndRow.id`
//! - PID, MID: parent IDs — convertf writes `0`, reads ignored
//! - sex: `1`=M, `2`=F, `0` or other = unknown
//! - phenotype: `-9` = missing (convertf default); we ignore on read
//!
//! # `familynames` option
//!
//! When reading PLINK input with `familynames: YES` (upstream default),
//! convertf concatenates `FID:IID` into the sample ID so cross-family
//! duplicate IIDs don't collide. With `familynames: NO`, only IID is used.
//!
//! # `outputgroup` option (write path)
//!
//! Upstream convertf writes phenotype column as:
//! - `outputgroup: NO` (default): `-9` (missing)
//! - `outputgroup: YES`: use `IndRow.pop` as the phenotype string
//!
//! We don't yet track case/control distinction; in both modes we write
//! `-9`. Add if/when a user needs it.

use super::{split_lines, IndRow, Sex};
use anyhow::{anyhow, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Read `.fam`. If `familynames == true`, the returned `IndRow.id` is
/// `"FID:IID"`; otherwise just `IID`.
pub fn read(path: &Path, familynames: bool) -> Result<Vec<IndRow>> {
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
        let row = parse_fam_line(line, familynames)
            .with_context(|| format!("{}:{}", path.display(), lineno))?;
        rows.push(row);
    }
    Ok(rows)
}

pub fn write(path: &Path, rows: &[IndRow], _outputgroup: bool) -> Result<()> {
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::with_capacity(64 * 1024, file);
    for r in rows {
        // FID IID PID MID sex pheno
        // We write PID=MID=0, pheno=-9. outputgroup handling deferred.
        writeln!(
            w,
            "{}\t{}\t0\t0\t{}\t-9",
            r.pop,
            r.id,
            r.sex.as_fam_code() as char
        )?;
    }
    w.flush()?;
    Ok(())
}

fn parse_fam_line(line: &[u8], familynames: bool) -> Result<IndRow> {
    let mut cols = line
        .split(|b: &u8| b.is_ascii_whitespace())
        .filter(|c| !c.is_empty());

    let fid = cols.next().ok_or_else(|| anyhow!("missing FID"))?;
    let iid = cols.next().ok_or_else(|| anyhow!("missing IID"))?;
    let _pid = cols.next().ok_or_else(|| anyhow!("missing PID"))?;
    let _mid = cols.next().ok_or_else(|| anyhow!("missing MID"))?;
    let sex_raw = cols.next().ok_or_else(|| anyhow!("missing sex"))?;
    let _pheno = cols.next(); // tolerate missing pheno column

    let fid = std::str::from_utf8(fid)?;
    let iid = std::str::from_utf8(iid)?;

    let id = if familynames {
        format!("{fid}:{iid}")
    } else {
        iid.to_owned()
    };

    let sex = Sex::from_char(*sex_raw.first().unwrap_or(&b'0'));
    let pop = fid.to_owned();
    // PLINK has no "Ignore" convention — always keep.
    Ok(IndRow {
        id,
        sex,
        pop,
        ignore: false,
    })
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
    fn reads_basic_fam() {
        let f = write_tmp(
            "French\tS_French-1.DG\t0\t0\t1\t-9\n\
             Onge\tBR_Onge-2.DG\t0\t0\t2\t-9\n",
        );
        let rows = read(f.path(), false).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id, "S_French-1.DG");
        assert_eq!(rows[0].pop, "French");
        assert_eq!(rows[0].sex, Sex::Male);
        assert_eq!(rows[1].sex, Sex::Female);
    }

    #[test]
    fn familynames_concatenates_fid_iid() {
        let f = write_tmp("Pop1\tSAMPLE01\t0\t0\t1\t-9\n");
        let rows = read(f.path(), true).unwrap();
        assert_eq!(rows[0].id, "Pop1:SAMPLE01");
    }

    #[test]
    fn sex_zero_is_unknown() {
        let f = write_tmp("P\tS1\t0\t0\t0\t-9\n");
        let rows = read(f.path(), false).unwrap();
        assert_eq!(rows[0].sex, Sex::Unknown);
    }

    #[test]
    fn roundtrip_preserves_id_pop_sex() {
        let rows = vec![
            IndRow {
                id: "S1".into(),
                sex: Sex::Male,
                pop: "Pop1".into(),
                ignore: false,
            },
            IndRow {
                id: "S2".into(),
                sex: Sex::Female,
                pop: "Pop2".into(),
                ignore: false,
            },
        ];
        let f = tempfile::NamedTempFile::new().unwrap();
        write(f.path(), &rows, false).unwrap();
        let got = read(f.path(), false).unwrap();
        assert_eq!(got.len(), 2);
        assert_eq!(got[0].id, "S1");
        assert_eq!(got[0].pop, "Pop1");
        assert_eq!(got[0].sex, Sex::Male);
        assert_eq!(got[1].sex, Sex::Female);
    }

    #[test]
    fn tolerates_missing_pheno_column() {
        // Some tools emit 5-column .fam
        let f = write_tmp("Pop\tS1\t0\t0\t1\n");
        let rows = read(f.path(), false).unwrap();
        assert_eq!(rows[0].sex, Sex::Male);
    }
}
