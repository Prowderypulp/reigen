//! PVAR reader/writer — PLINK 2 variant metadata file.
//!
//! PVAR is a TSV with a `#`-prefixed header line. The canonical columns
//! are `#CHROM POS ID REF ALT` (plus optional `CM`, `QUAL`, `FILTER`, `INFO`).
//!
//! # Allele convention (load-bearing — see `PGEN_PLAN.md` §2)
//!
//! PVAR ALT = `.bim` A1 = reigen `allele1`. PVAR REF = reigen `allele2`.
//! **No swap.**

use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::meta::SnpRow;

/// A variant is multiallelic iff its ALT field lists more than one allele
/// (comma-separated). reigen's model is biallelic-only, so these are dropped
/// (see `PGEN_PLAN.md` §4.4). The PGEN reader is given a matching keep-mask so
/// the genotype stream stays aligned with the kept `SnpRow`s.
fn is_multiallelic(alt: &str) -> bool {
    alt.contains(',')
}

/// A variant dropped from the biallelic stream, captured with its original
/// PVAR field values for the `.pgen-drop.tsv` report (plan §4.4, acceptance #5).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DroppedVariant {
    pub chrom: String,
    pub pos: String,
    pub id: String,
    pub ref_allele: String,
    pub alt: String,
    pub reason: &'static str,
}

/// Read a PVAR, returning only the **biallelic** variants as `SnpRow`s.
/// Multiallelic (comma-ALT) variants are dropped; use [`read_with_keep`] when
/// you also need the on-disk keep-mask (e.g. to align the PGEN reader).
pub fn read(path: &Path, numchrom: u32) -> Result<Vec<SnpRow>> {
    Ok(read_inner(path, Some(numchrom))?.0)
}

/// The variants dropped as multiallelic (comma-ALT), with their original PVAR
/// field values — used to write the `.pgen-drop.tsv` report.
pub fn dropped_multiallelic(path: &Path) -> Result<Vec<DroppedVariant>> {
    Ok(read_inner(path, None)?.2)
}

/// Write a `.pgen-drop.tsv` report listing variants dropped from the biallelic
/// stream (mirrors merge's `.missnp`). No-op when `dropped` is empty.
pub fn write_drop_report(path: &Path, dropped: &[DroppedVariant]) -> Result<()> {
    if dropped.is_empty() {
        return Ok(());
    }
    let f = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::new(f);
    writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\treason")?;
    for d in dropped {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            d.chrom, d.pos, d.id, d.ref_allele, d.alt, d.reason
        )?;
    }
    w.flush()?;
    Ok(())
}

/// The per-on-disk-variant keep-mask: `keep[i]` is true iff on-disk variant `i`
/// is biallelic (kept). `keep.len()` equals the PVAR data-row count, which must
/// equal the PGEN header's variant count. Threaded into `PgenReader::open` so
/// the reader yields exactly the kept variants, in order. The mask depends only
/// on the ALT column (comma ⇒ multiallelic), so it needs no `numchrom`.
pub fn keep_ondisk_mask(path: &Path) -> Result<Vec<bool>> {
    Ok(read_inner(path, None)?.1)
}

/// Read a PVAR, returning `(biallelic SnpRows, on-disk keep-mask)`. The mask has
/// one entry per data row (true = biallelic/kept); `SnpRows` contains only the
/// kept rows, so `rows.len() == keep.iter().filter(|&&k| k).count()`.
pub fn read_with_keep(path: &Path, numchrom: u32) -> Result<(Vec<SnpRow>, Vec<bool>)> {
    let (rows, keep, _dropped) = read_inner(path, Some(numchrom))?;
    Ok((rows, keep))
}

/// Core PVAR parse. When `numchrom` is `Some`, full `SnpRow`s are built (chrom
/// parsed via the shared [`crate::meta::bim::parse_chrom`] so `.pvar` and `.bim`
/// agree on the PLINK↔internal XY/MT remap). When `None`, only the keep-mask is
/// computed (no chrom/pos parsing) — used by [`keep_ondisk_mask`], which has no
/// `numchrom` and doesn't need one.
#[allow(clippy::type_complexity)]
fn read_inner(
    path: &Path,
    numchrom: Option<u32>,
) -> Result<(Vec<SnpRow>, Vec<bool>, Vec<DroppedVariant>)> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut lines = reader.lines();

    // PVAR may carry VCF-style `##` meta-information lines before the column
    // header. Skip them; the column header is the first line starting with a
    // single `#` (i.e. `#` but not `##`).
    let header = loop {
        let line = lines
            .next()
            .context("pvar: no #CHROM header line found")?
            .context("pvar: error reading header")?;
        if line.starts_with("##") {
            continue;
        }
        if !line.starts_with('#') {
            bail!("pvar: column header does not start with '#'");
        }
        break line;
    };

    let columns: Vec<String> = header
        .trim_start_matches('#')
        .split('\t')
        .map(|s| s.trim().to_ascii_uppercase().to_string())
        .collect();

    let find = |name: &str| columns.iter().position(|c| c == name);

    let idx_chrom = find("CHROM");
    let idx_pos = find("POS");
    let idx_id = find("ID");
    let idx_alt = find("ALT");
    let idx_ref = find("REF");
    // CM is optional
    let idx_cm = find("CM");

    if idx_chrom.is_none() || idx_pos.is_none() || idx_id.is_none() || idx_alt.is_none() || idx_ref.is_none() {
        bail!("pvar: missing required columns (need #CHROM POS ID REF ALT)");
    }

    let i_chrom = idx_chrom.unwrap();
    let i_pos = idx_pos.unwrap();
    let i_id = idx_id.unwrap();
    let i_alt = idx_alt.unwrap();
    let i_ref = idx_ref.unwrap();

    let mut rows = Vec::new();
    let mut keep = Vec::new();
    let mut dropped: Vec<DroppedVariant> = Vec::new();
    for (lineno, line) in lines.enumerate() {
        let line = line.with_context(|| format!("pvar: error reading line {}", lineno + 2))?;
        let fields: Vec<&str> = line.split('\t').collect();
        let max = std::cmp::max(i_alt, i_ref) + 1;
        if fields.len() < max {
            bail!("pvar:{} too few columns (got {}, need {})", lineno + 2, fields.len(), max);
        }

        let alt = fields[i_alt].trim();
        // Multiallelic variants are out of reigen's biallelic model: record a
        // drop in the keep-mask and capture the original fields for the report.
        if is_multiallelic(alt) {
            keep.push(false);
            let get = |i: usize| fields.get(i).map_or("", |s| s.trim()).to_string();
            dropped.push(DroppedVariant {
                chrom: get(i_chrom),
                pos: get(i_pos),
                id: get(i_id),
                ref_allele: get(i_ref),
                alt: alt.to_string(),
                reason: "multiallelic",
            });
            continue;
        }
        keep.push(true);

        // Mask-only mode: don't parse chrom/pos/alleles (and don't require a
        // numchrom). Keeps `keep_ondisk_mask` independent of chromosome coding.
        let numchrom = match numchrom {
            Some(nc) => nc,
            None => continue,
        };

        // Reuse the .bim chrom parser so PVAR and BIM agree on the PLINK↔internal
        // XY/MT remap (PLINK 25=XY/26=MT vs internal 25=MT/26=XY). Parsing it by
        // hand here previously swapped XY and MT on every PGEN↔BED conversion.
        let chrom = crate::meta::bim::parse_chrom(fields[i_chrom].trim().as_bytes(), numchrom)
            .with_context(|| format!("pvar:{} bad CHROM", lineno + 2))?;

        let pos: u64 = fields[i_pos]
            .trim()
            .parse()
            .with_context(|| format!("pvar:{} bad POS", lineno + 2))?;

        let genetic_pos = idx_cm
            .and_then(|ic| fields.get(ic))
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);

        let ref_str = fields[i_ref].trim();

        let allele1 = alt.as_bytes()[0];
        let allele2 = ref_str.as_bytes()[0];

        rows.push(SnpRow {
            id: fields[i_id].trim().to_string(),
            chrom,
            genetic_pos,
            physical_pos: pos,
            allele1,
            allele2,
        });
    }

    // Warn only on the canonical full read (numchrom present); the mask-only
    // and drop-report passes reuse this and would otherwise double/triple-log.
    if !dropped.is_empty() && numchrom.is_some() {
        log::warn!(
            "pvar: dropped {} multiallelic variant(s) (comma-ALT) from {}; \
             reigen handles biallelic hardcalls only (see .pgen-drop.tsv)",
            dropped.len(),
            path.display()
        );
    }
    Ok((rows, keep, dropped))
}

pub fn write(path: &Path, rows: &[SnpRow], numchrom: u32) -> Result<()> {
    let f = File::create(path).with_context(|| format!("create {}", path.display()))?;
    let mut w = BufWriter::new(f);
    writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tCM")?;
    for s in rows {
        // Emit PLINK numeric chrom codes via the shared .bim helper (handles the
        // internal↔PLINK XY/MT remap), so `.pvar` and `.bim` are identical in
        // the chrom column and PGEN↔BED conversions preserve chromosomes.
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.6}",
            crate::meta::bim::chrom_to_plink_numeric(s.chrom, numchrom),
            s.physical_pos,
            s.id,
            s.allele2 as char,
            s.allele1 as char,
            s.genetic_pos
        )?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_roundtrip() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let rows = vec![
            SnpRow {
                id: "rs0".into(),
                chrom: 1,
                genetic_pos: 0.001,
                physical_pos: 1000,
                allele1: b'A',
                allele2: b'G',
            },
            SnpRow {
                id: "rs1".into(),
                chrom: 23,
                genetic_pos: 0.0,
                physical_pos: 5000,
                allele1: b'T',
                allele2: b'C',
            },
        ];
        write(tmp.path(), &rows, 22).unwrap();
        let back = read(tmp.path(), 22).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back[0].id, "rs0");
        assert_eq!(back[0].allele1, b'A');
        assert_eq!(back[0].allele2, b'G');
        assert_eq!(back[1].id, "rs1");
        assert_eq!(back[1].chrom, 23);
    }

    #[test]
    fn skips_double_hash_meta_lines() {
        // Real plink2 PVARs carry ## meta lines before the #CHROM header.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\n1\t1000\trs0\tG\tA\n",
        )
        .unwrap();
        let (rows, keep) = read_with_keep(tmp.path(), 22).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].id, "rs0");
        assert_eq!(rows[0].allele1, b'A'); // ALT
        assert_eq!(rows[0].allele2, b'G'); // REF
        assert_eq!(keep, vec![true]);
    }

    #[test]
    fn drops_multiallelic_and_builds_mask() {
        // comma-ALT variant is dropped; mask records its on-disk position.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "#CHROM\tPOS\tID\tREF\tALT\n\
             1\t1\trs0\tG\tA\n\
             1\t2\trs1\tC\tA,T\n\
             1\t3\trs2\tT\tC\n",
        )
        .unwrap();
        let (rows, keep) = read_with_keep(tmp.path(), 22).unwrap();
        assert_eq!(keep, vec![true, false, true], "rs1 (comma-ALT) dropped in mask");
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id, "rs0");
        assert_eq!(rows[1].id, "rs2");
    }

    #[test]
    fn reads_real_multi_fixture() {
        let path = format!(
            "{}/tests/golden/plink2/multi.pvar",
            env!("CARGO_MANIFEST_DIR")
        );
        let (rows, keep) = read_with_keep(Path::new(&path), 22).unwrap();
        // multi.pvar: 4 on-disk variants, rs1 triallelic (ALT=A,T).
        assert_eq!(keep, vec![true, false, true, true]);
        assert_eq!(rows.len(), 3);
        assert_eq!(rows.iter().map(|r| r.id.as_str()).collect::<Vec<_>>(), ["rs0", "rs2", "rs3"]);
    }

    #[test]
    fn tolerates_extra_columns_qual_filter_info() {
        // PVAR with QUAL, FILTER, INFO columns (plink2 default output).
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t1000\trs0\tG\tA\t30\tPASS\tAC=2\n\
             1\t2000\trs1\tC\tT\t.\t.\t.\n",
        )
        .unwrap();
        let rows = read(tmp.path(), 22).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id, "rs0");
        assert_eq!(rows[0].allele1, b'A');
        assert_eq!(rows[0].allele2, b'G');
        assert_eq!(rows[1].id, "rs1");
    }

    #[test]
    fn column_order_with_sid_and_no_cm() {
        // PVAR with SID column and no CM — plink2 can emit this.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "#CHROM\tPOS\tID\tREF\tALT\n\
             1\t1000\trs0\tG\tA\n\
             2\t5000\trs1\tC\tT\n",
        )
        .unwrap();
        let rows = read(tmp.path(), 22).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].chrom, 1);
        assert_eq!(rows[0].physical_pos, 1000);
        assert_eq!(rows[0].genetic_pos, 0.0, "CM absent defaults to 0.0");
        assert_eq!(rows[1].chrom, 2);
    }

    #[test]
    fn write_read_preserves_alleles_without_cm() {
        // Round-trip where CM is 0.0 (the default when not in input).
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let rows = vec![
            SnpRow {
                id: "rs0".into(),
                chrom: 1,
                genetic_pos: 0.0,
                physical_pos: 1000,
                allele1: b'A',
                allele2: b'G',
            },
        ];
        write(tmp.path(), &rows, 22).unwrap();
        let back = read(tmp.path(), 22).unwrap();
        assert_eq!(back.len(), 1);
        assert_eq!(back[0].allele1, b'A');
        assert_eq!(back[0].allele2, b'G');
        assert_eq!(back[0].genetic_pos, 0.0);
    }

    /// Regression: PVAR must use the same PLINK↔internal XY/MT remap as `.bim`
    /// (PLINK 25=XY/26=MT vs internal 25=MT/26=XY). Parsing XY/MT by hand here
    /// previously swapped the two on every PGEN↔BED conversion.
    #[test]
    fn chrom_xy_mt_matches_bim_convention() {
        // Internal codes: 23=X, 24=Y, 25=MT, 26=XY (numchrom=22).
        let rows = vec![
            SnpRow { id: "rX".into(), chrom: 23, genetic_pos: 0.0, physical_pos: 1, allele1: b'A', allele2: b'G' },
            SnpRow { id: "rY".into(), chrom: 24, genetic_pos: 0.0, physical_pos: 2, allele1: b'A', allele2: b'G' },
            SnpRow { id: "rMT".into(), chrom: 25, genetic_pos: 0.0, physical_pos: 3, allele1: b'A', allele2: b'G' },
            SnpRow { id: "rXY".into(), chrom: 26, genetic_pos: 0.0, physical_pos: 4, allele1: b'A', allele2: b'G' },
        ];
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write(tmp.path(), &rows, 22).unwrap();

        // Writer emits PLINK numeric: X=23, Y=24, MT=26, XY=25.
        let txt = std::fs::read_to_string(tmp.path()).unwrap();
        let codes: Vec<&str> = txt.lines().skip(1).map(|l| l.split('\t').next().unwrap()).collect();
        assert_eq!(codes, ["23", "24", "26", "25"], "PLINK numeric: MT=26, XY=25");

        // Round-trip back to the same internal codes (no XY/MT swap).
        let back = read(tmp.path(), 22).unwrap();
        assert_eq!(back.iter().map(|r| r.chrom).collect::<Vec<_>>(), vec![23, 24, 25, 26]);

        // A PVAR written with numeric PLINK codes parses to the same internal.
        let tmp2 = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp2.path(),
            "#CHROM\tPOS\tID\tREF\tALT\n23\t1\trX\tG\tA\n24\t2\trY\tG\tA\n26\t3\trMT\tG\tA\n25\t4\trXY\tG\tA\n",
        )
        .unwrap();
        let back2 = read(tmp2.path(), 22).unwrap();
        assert_eq!(back2.iter().map(|r| r.chrom).collect::<Vec<_>>(), vec![23, 24, 25, 26]);
    }

    #[test]
    fn dropped_multiallelic_captures_original_fields() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "#CHROM\tPOS\tID\tREF\tALT\n\
             1\t1\trs0\tG\tA\n\
             1\t2\trs1\tC\tA,T\n\
             2\t3\trs2\tT\tG,C,A\n",
        )
        .unwrap();
        let dropped = dropped_multiallelic(tmp.path()).unwrap();
        assert_eq!(dropped.len(), 2);
        assert_eq!(dropped[0].id, "rs1");
        assert_eq!(dropped[0].alt, "A,T");
        assert_eq!(dropped[0].ref_allele, "C");
        assert_eq!(dropped[0].reason, "multiallelic");
        assert_eq!(dropped[1].id, "rs2");
        assert_eq!(dropped[1].alt, "G,C,A");

        // write_drop_report round-trips the captured fields.
        let out = tempfile::NamedTempFile::new().unwrap();
        write_drop_report(out.path(), &dropped).unwrap();
        let txt = std::fs::read_to_string(out.path()).unwrap();
        assert!(txt.starts_with("#CHROM\tPOS\tID\tREF\tALT\treason\n"));
        assert!(txt.contains("1\t2\trs1\tC\tA,T\tmultiallelic"));

        // No file written for an empty drop list.
        let empty = tempfile::NamedTempFile::new().unwrap();
        let p = empty.path().to_path_buf();
        drop(empty);
        write_drop_report(&p, &[]).unwrap();
        assert!(!p.exists(), "no report for zero drops");
    }
}
