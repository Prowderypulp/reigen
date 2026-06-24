//! PSAM reader/writer — PLINK 2 sample metadata file.
//!
//! PSAM is a TSV with a `#`-prefixed header line. The canonical columns
//! are `#FID IID` (required) plus optional `SID`, `PAT`, `MAT`, `SEX`, `PHENO1`.

use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::meta::{IndRow, Sex};

pub fn read(path: &Path) -> Result<Vec<IndRow>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut lines = reader.lines();

    // PSAM may carry `##` meta-information lines before the column header.
    // Skip them; the header is the first line starting with a single `#`.
    let header = loop {
        let line = lines
            .next()
            .context("psam: no header line found")?
            .context("psam: error reading header")?;
        if line.starts_with("##") {
            continue;
        }
        if !line.starts_with('#') {
            bail!("psam: column header does not start with '#'");
        }
        break line;
    };

    let columns: Vec<String> = header
        .trim_start_matches('#')
        .split('\t')
        .map(|s| s.trim().to_ascii_uppercase().to_string())
        .collect();

    let find = |name: &str| columns.iter().position(|c| c == name);

    // IID is the only required column; FID is optional (plink2 emits `#IID`-only
    // PSAMs when there is no family ID). When FID is absent it defaults to "0",
    // i.e. the sample ID is just the IID.
    let idx_fid = find("FID");
    let idx_iid = find("IID");
    let idx_sex = find("SEX");

    let i_iid = idx_iid.context("psam: missing required column (need #IID)")?;

    let mut rows = Vec::new();
    for (lineno, line) in lines.enumerate() {
        let line = line.with_context(|| format!("psam: error reading line {}", lineno + 2))?;
        let fields: Vec<&str> = line.split('\t').collect();
        let max = idx_fid.map_or(i_iid, |f| std::cmp::max(f, i_iid)) + 1;
        if fields.len() < max {
            bail!("psam:{} too few columns", lineno + 2);
        }

        let fid = idx_fid.and_then(|f| fields.get(f)).map_or("0", |s| s.trim());
        let iid = fields[i_iid].trim();

        let sex = idx_sex
            .and_then(|is| fields.get(is))
            .map(|s| s.trim().as_bytes().first().copied().unwrap_or(0))
            .map(Sex::from_char)
            .unwrap_or(Sex::Unknown);

        let id = if fid == "0" {
            iid.to_string()
        } else {
            format!("{fid}:{iid}")
        };

        rows.push(IndRow {
            id,
            sex,
            pop: fid.to_string(),
            ignore: false,
        });
    }
    Ok(rows)
}

pub fn write(path: &Path, rows: &[IndRow], outputgroup: bool) -> Result<()> {
    let f = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::new(f);
    if outputgroup {
        writeln!(w, "#FID\tIID\tSEX\tPHENO1")?;
    } else {
        writeln!(w, "#FID\tIID\tSEX")?;
    }
    for ind in rows {
        let fid = &ind.pop;
        let iid = plink_iid(ind);
        let sex_str = match ind.sex {
            Sex::Male => "1",
            Sex::Female => "2",
            Sex::Unknown => "0",
        };
        if outputgroup {
            writeln!(w, "{fid}\t{iid}\t{sex_str}\t{}", ind.pop)?;
        } else {
            writeln!(w, "{fid}\t{iid}\t{sex_str}")?;
        }
    }
    w.flush()?;
    Ok(())
}

fn plink_iid(row: &IndRow) -> &str {
    if let Some((fid, iid)) = row.id.split_once(':') {
        if fid == row.pop {
            return iid;
        }
    }
    row.id.as_str()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_roundtrip() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let rows = vec![
            IndRow {
                id: "ID1".into(),
                sex: Sex::Unknown,
                pop: "0".into(),
                ignore: false,
            },
            IndRow {
                id: "ID2".into(),
                sex: Sex::Male,
                pop: "0".into(),
                ignore: false,
            },
        ];
        write(tmp.path(), &rows, false).unwrap();
        let back = read(tmp.path()).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back[0].id, "ID1");
        assert_eq!(back[1].id, "ID2");
        assert_eq!(back[1].sex, Sex::Male);
    }

    #[test]
    fn fid_iid_concat() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), "#FID\tIID\tSEX\nPOP\tID1\t1\n").unwrap();
        let rows = read(tmp.path()).unwrap();
        assert_eq!(rows[0].id, "POP:ID1");
        assert_eq!(rows[0].sex, Sex::Male);
    }

    #[test]
    fn fid_zero_uses_iid_only() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), "#FID\tIID\tSEX\n0\tID1\t2\n").unwrap();
        let rows = read(tmp.path()).unwrap();
        assert_eq!(rows[0].id, "ID1");
    }

    #[test]
    fn tolerates_sid_pat_mat_columns() {
        // PSAM with SID, PAT, MAT columns (plink2 default output).
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "#FID\tIID\tSID\tPAT\tMAT\tSEX\n\
             FAM1\tID1\tS1\t0\t0\t1\n\
             FAM1\tID2\tS2\t0\t0\t2\n",
        )
        .unwrap();
        let rows = read(tmp.path()).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id, "FAM1:ID1");
        assert_eq!(rows[0].sex, Sex::Male);
        assert_eq!(rows[1].id, "FAM1:ID2");
        assert_eq!(rows[1].sex, Sex::Female);
    }

    #[test]
    fn iid_only_header_no_fid() {
        // plink2 emits #IID-only PSAM when there's no family ID.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "#IID\tSEX\nID1\t1\nID2\t0\n",
        )
        .unwrap();
        let rows = read(tmp.path()).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id, "ID1");
        assert_eq!(rows[0].sex, Sex::Male);
        assert_eq!(rows[1].id, "ID2");
        assert_eq!(rows[1].sex, Sex::Unknown);
    }

    #[test]
    fn tolerates_pheno_columns() {
        // PSAM with PHENO1 column (and no FID).
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "#IID\tSEX\tPHENO1\nID1\t1\tCase\nID2\t2\tControl\n",
        )
        .unwrap();
        let rows = read(tmp.path()).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].sex, Sex::Male);
        assert_eq!(rows[1].sex, Sex::Female);
    }

    #[test]
    fn skips_double_hash_meta_lines() {
        // PSAM with ## meta-lines before the header (like PVAR).
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "##PSAMv1.0\n#FID\tIID\tSEX\n0\tID1\t1\n",
        )
        .unwrap();
        let rows = read(tmp.path()).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].id, "ID1");
    }

    #[test]
    fn write_strips_matching_family_prefix_from_iid() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let rows = vec![IndRow {
            id: "FAM1:ID1".into(),
            sex: Sex::Male,
            pop: "FAM1".into(),
            ignore: false,
        }];
        write(tmp.path(), &rows, false).unwrap();
        let text = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(text, "#FID\tIID\tSEX\nFAM1\tID1\t1\n");
    }

    #[test]
    fn outputgroup_writes_population_to_pheno1() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let rows = vec![IndRow {
            id: "ID1".into(),
            sex: Sex::Female,
            pop: "Pop1".into(),
            ignore: false,
        }];
        write(tmp.path(), &rows, true).unwrap();
        let text = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(text, "#FID\tIID\tSEX\tPHENO1\nPop1\tID1\t2\tPop1\n");
    }
}
