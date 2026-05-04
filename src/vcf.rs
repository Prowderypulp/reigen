//! Minimal VCF 4.3 reader and writer for biallelic SNPs (GT field only).
//!
//! This is intentionally narrow in scope: we only handle biallelic SNPs with
//! a single GT FORMAT field. Multi-allelic records, indels, and additional
//! FORMAT fields are skipped (with warnings) on import. On export, only GT
//! is emitted.
//!
//! No external VCF library is used — the format is simple enough for our
//! biallelic-SNP-only use case.

use crate::geno::codec;
use crate::meta::{IndRow, SnpRow};
use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// ======================================================================
// VCF Writer
// ======================================================================

/// Write a VCF 4.3 file from genotype data.
///
/// `genotypes` is a SNP-major matrix: `genotypes[snp_idx][sample_idx]` = 0/1/2/9.
pub fn write_vcf(
    path: &Path,
    snps: &[SnpRow],
    inds: &[IndRow],
    genotypes: &[Vec<u8>],
    chr_prefix: &str,
    numchrom: u32,
) -> Result<()> {
    let file =
        std::fs::File::create(path).with_context(|| format!("create VCF {}", path.display()))?;
    let mut w = BufWriter::new(file);

    // --- Header ---
    writeln!(w, "##fileformat=VCFv4.3")?;
    writeln!(w, "##source=reigen (population genomics toolkit)")?;

    // Contig lines for autosomes + X/Y/MT
    for c in 1..=numchrom {
        writeln!(w, "##contig=<ID={}{}>", chr_prefix, c)?;
    }
    writeln!(w, "##contig=<ID={}X>", chr_prefix)?;
    writeln!(w, "##contig=<ID={}Y>", chr_prefix)?;
    writeln!(w, "##contig=<ID={}MT>", chr_prefix)?;

    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;

    // Column header
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for ind in inds {
        write!(w, "\t{}", ind.id)?;
    }
    writeln!(w)?;

    // --- Records ---
    for (snp_idx, snp) in snps.iter().enumerate() {
        let chrom_str = chrom_to_vcf_string(snp.chrom, numchrom, chr_prefix);
        let ref_allele = if snp.allele1 == b'0' || snp.allele1 == b'X' {
            '.'
        } else {
            snp.allele1 as char
        };
        let alt_allele = if snp.allele2 == b'0' || snp.allele2 == b'X' {
            '.'
        } else {
            snp.allele2 as char
        };

        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t.\tPASS\t.\tGT",
            chrom_str, snp.physical_pos, snp.id, ref_allele, alt_allele
        )?;

        let genos = &genotypes[snp_idx];
        for sample_idx in 0..inds.len() {
            let g = if sample_idx < genos.len() {
                genos[sample_idx]
            } else {
                codec::G_MISSING
            };
            let gt = match g {
                0 => "0/0",
                1 => "0/1",
                2 => "1/1",
                _ => "./.",
            };
            write!(w, "\t{}", gt)?;
        }
        writeln!(w)?;
    }

    w.flush()?;
    Ok(())
}

fn chrom_to_vcf_string(chrom: u8, numchrom: u32, prefix: &str) -> String {
    let Some(nc) = normalized_numchrom(numchrom) else {
        return format!("{}{}", prefix, chrom);
    };

    if chrom >= 1 && chrom <= nc {
        format!("{}{}", prefix, chrom)
    } else if Some(chrom) == nc.checked_add(1) {
        format!("{}X", prefix)
    } else if Some(chrom) == nc.checked_add(2) {
        format!("{}Y", prefix)
    } else if Some(chrom) == nc.checked_add(3) {
        format!("{}MT", prefix)
    } else if Some(chrom) == nc.checked_add(4) {
        format!("{}XY", prefix)
    } else {
        format!("{}{}", prefix, chrom)
    }
}

// ======================================================================
// VCF Reader
// ======================================================================

/// Parsed VCF record (biallelic SNP only).
#[derive(Debug)]
pub struct VcfRecord {
    pub chrom: u8,
    pub pos: u64,
    pub id: String,
    pub ref_allele: u8,
    pub alt_allele: u8,
    /// Genotypes per sample: 0/1/2/9(missing).
    pub genotypes: Vec<u8>,
}

/// Statistics from VCF reading.
#[derive(Debug, Default)]
pub struct VcfReadStats {
    pub total_records: usize,
    pub kept_biallelic_snp: usize,
    pub skipped_multiallelic: usize,
    pub skipped_indel: usize,
    pub skipped_no_gt: usize,
    pub skipped_no_alt: usize,
    pub skipped_ref_filtered: usize,
    pub skipped_snplist: usize,
}

/// Read a VCF file, returning sample names and biallelic SNP records.
///
/// Multi-allelic records, indels, and records without GT are skipped with
/// counters tracked in `VcfReadStats`.
pub fn read_vcf(
    path: &Path,
    numchrom: u32,
    ref_filter: Option<&std::collections::HashMap<(u8, u64), usize>>,
    snplist: Option<&HashSet<String>>,
) -> Result<(Vec<String>, Vec<VcfRecord>, VcfReadStats)> {
    if normalized_numchrom(numchrom).is_none() {
        bail!("numchrom {} is too large (max 251)", numchrom);
    }

    let file = std::fs::File::open(path).with_context(|| format!("open VCF {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut sample_names: Vec<String> = Vec::new();
    let mut records: Vec<VcfRecord> = Vec::new();
    let mut stats = VcfReadStats::default();
    let mut header_parsed = false;

    for line_result in reader.lines() {
        let line = line_result.with_context(|| format!("read line from {}", path.display()))?;
        let line = line.trim_end();

        if line.starts_with("##") {
            continue;
        }

        if line.starts_with("#CHROM") {
            // Parse sample names from header
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 9 {
                sample_names = cols[9..].iter().map(|s| s.to_string()).collect();
            }
            header_parsed = true;
            continue;
        }

        if !header_parsed {
            continue;
        }

        stats.total_records += 1;

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }

        let chrom_str = cols[0];
        let pos: u64 = match cols[1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let id = cols[2];
        let ref_field = cols[3];
        let alt_field = cols[4];
        let format_field = cols[8];

        // Skip if REF or ALT is not a single base (indel)
        if ref_field.len() != 1 || ref_field == "." {
            stats.skipped_indel += 1;
            continue;
        }

        // Skip multi-allelic (ALT has comma)
        if alt_field.contains(',') {
            stats.skipped_multiallelic += 1;
            continue;
        }

        // Skip monomorphic or missing ALT
        if alt_field == "." {
            stats.skipped_no_alt += 1;
            continue;
        }

        // Skip if ALT is not a single base
        if alt_field.len() != 1 {
            stats.skipped_indel += 1;
            continue;
        }

        let ref_allele = ref_field.as_bytes()[0].to_ascii_uppercase();
        let alt_allele = alt_field.as_bytes()[0].to_ascii_uppercase();

        // Must be ACGT
        if !matches!(ref_allele, b'A' | b'C' | b'G' | b'T')
            || !matches!(alt_allele, b'A' | b'C' | b'G' | b'T')
        {
            stats.skipped_indel += 1;
            continue;
        }

        // Parse chromosome
        let chrom = parse_vcf_chrom(chrom_str, numchrom);
        if chrom == 0 {
            continue;
        }

        // Optional: filter to reference positions
        if let Some(ref_map) = ref_filter {
            if !ref_map.contains_key(&(chrom, pos)) {
                stats.skipped_ref_filtered += 1;
                continue;
            }
        }

        // Optional: filter to snplist by rsID
        if let Some(ids) = snplist {
            if id == "." || !ids.contains(id) {
                stats.skipped_snplist += 1;
                continue;
            }
        }

        // Find GT index in FORMAT field
        let format_parts: Vec<&str> = format_field.split(':').collect();
        let gt_idx = match format_parts.iter().position(|&f| f == "GT") {
            Some(idx) => idx,
            None => {
                stats.skipped_no_gt += 1;
                continue;
            }
        };

        // Parse genotypes
        let n_samples = sample_names.len();
        let mut genotypes = vec![codec::G_MISSING; n_samples];

        for (i, &sample_field) in cols.iter().skip(9).enumerate() {
            if i >= n_samples {
                break;
            }
            let parts: Vec<&str> = sample_field.split(':').collect();
            if gt_idx < parts.len() {
                genotypes[i] = parse_gt(parts[gt_idx]);
            }
        }

        let record_id = if id == "." {
            format!("{}:{}", chrom_str, pos)
        } else {
            id.to_string()
        };

        records.push(VcfRecord {
            chrom,
            pos,
            id: record_id,
            ref_allele,
            alt_allele,
            genotypes,
        });
        stats.kept_biallelic_snp += 1;
    }

    Ok((sample_names, records, stats))
}

/// Parse a VCF GT field (e.g., "0/0", "0|1", "1/1", "./.") into 0/1/2/9.
fn parse_gt(gt: &str) -> u8 {
    let sep = if gt.contains('|') { '|' } else { '/' };
    let parts: Vec<&str> = gt.splitn(2, sep).collect();
    if parts.len() != 2 {
        return codec::G_MISSING;
    }
    let a1 = match parts[0] {
        "0" => 0u8,
        "1" => 1u8,
        _ => return codec::G_MISSING,
    };
    let a2 = match parts[1] {
        "0" => 0u8,
        "1" => 1u8,
        _ => return codec::G_MISSING,
    };
    a1 + a2 // 0+0=0, 0+1=1, 1+0=1, 1+1=2
}

/// Parse VCF chromosome string to internal u8 representation.
fn parse_vcf_chrom(s: &str, numchrom: u32) -> u8 {
    // Strip one leading "chr" prefix (case-insensitive).
    let s = s
        .strip_prefix("chr")
        .or_else(|| s.strip_prefix("Chr"))
        .or_else(|| s.strip_prefix("CHR"))
        .unwrap_or(s);

    if let Ok(c) = s.parse::<u8>() {
        return c;
    }

    let Some(nc) = normalized_numchrom(numchrom) else {
        return 0;
    };

    match s.to_ascii_uppercase().as_str() {
        "X" => nc.checked_add(1).unwrap_or(0),
        "Y" => nc.checked_add(2).unwrap_or(0),
        "M" | "MT" => nc.checked_add(3).unwrap_or(0),
        "XY" => nc.checked_add(4).unwrap_or(0),
        _ => 0, // unknown → skip
    }
}

fn normalized_numchrom(numchrom: u32) -> Option<u8> {
    let nc = u8::try_from(numchrom).ok()?;
    nc.checked_add(4)?;
    Some(nc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_gt_variants() {
        assert_eq!(parse_gt("0/0"), 0);
        assert_eq!(parse_gt("0/1"), 1);
        assert_eq!(parse_gt("1/0"), 1);
        assert_eq!(parse_gt("1/1"), 2);
        assert_eq!(parse_gt("./."), codec::G_MISSING);
        assert_eq!(parse_gt(".|."), codec::G_MISSING);
        // Phased
        assert_eq!(parse_gt("0|0"), 0);
        assert_eq!(parse_gt("0|1"), 1);
        assert_eq!(parse_gt("1|1"), 2);
    }

    #[test]
    fn parse_vcf_chrom_variants() {
        assert_eq!(parse_vcf_chrom("1", 22), 1);
        assert_eq!(parse_vcf_chrom("22", 22), 22);
        assert_eq!(parse_vcf_chrom("chr1", 22), 1);
        assert_eq!(parse_vcf_chrom("chrX", 22), 23);
        assert_eq!(parse_vcf_chrom("chrChr1", 22), 0);
        assert_eq!(parse_vcf_chrom("X", 22), 23);
        assert_eq!(parse_vcf_chrom("Y", 22), 24);
        assert_eq!(parse_vcf_chrom("MT", 22), 25);
        assert_eq!(parse_vcf_chrom("chrM", 22), 25);
        assert_eq!(parse_vcf_chrom("XY", 22), 26);
        assert_eq!(parse_vcf_chrom("UNKNOWN", 22), 0);
    }

    #[test]
    fn parse_vcf_chrom_rejects_large_numchrom() {
        assert_eq!(parse_vcf_chrom("X", 252), 0);
    }

    #[test]
    fn chrom_to_vcf_roundtrip() {
        for c in 1..=22u8 {
            let s = chrom_to_vcf_string(c, 22, "");
            assert_eq!(parse_vcf_chrom(&s, 22), c);
        }
        assert_eq!(parse_vcf_chrom(&chrom_to_vcf_string(23, 22, "chr"), 22), 23);
    }

    #[test]
    fn write_and_read_vcf_roundtrip() {
        let snps = vec![
            SnpRow {
                id: "rs1".into(),
                chrom: 1,
                genetic_pos: 0.0,
                physical_pos: 100,
                allele1: b'A',
                allele2: b'G',
            },
            SnpRow {
                id: "rs2".into(),
                chrom: 2,
                genetic_pos: 0.0,
                physical_pos: 200,
                allele1: b'C',
                allele2: b'T',
            },
        ];
        let inds = vec![
            IndRow {
                id: "S1".into(),
                sex: crate::meta::Sex::Male,
                pop: "Pop1".into(),
                ignore: false,
            },
            IndRow {
                id: "S2".into(),
                sex: crate::meta::Sex::Female,
                pop: "Pop2".into(),
                ignore: false,
            },
        ];
        // SNP-major: genotypes[snp][sample]
        let genotypes = vec![
            vec![0u8, 1],              // rs1: S1=hom_ref, S2=het
            vec![2, codec::G_MISSING], // rs2: S1=hom_alt, S2=missing
        ];

        let dir = tempfile::tempdir().unwrap();
        let vcf_path = dir.path().join("test.vcf");
        write_vcf(&vcf_path, &snps, &inds, &genotypes, "", 22).unwrap();

        // Read it back
        let (samples, records, stats) = read_vcf(&vcf_path, 22, None, None).unwrap();
        assert_eq!(samples, vec!["S1", "S2"]);
        assert_eq!(records.len(), 2);
        assert_eq!(stats.kept_biallelic_snp, 2);

        // Check record 0
        assert_eq!(records[0].chrom, 1);
        assert_eq!(records[0].pos, 100);
        assert_eq!(records[0].id, "rs1");
        assert_eq!(records[0].genotypes, vec![0, 1]);

        // Check record 1
        assert_eq!(records[1].chrom, 2);
        assert_eq!(records[1].pos, 200);
        assert_eq!(records[1].id, "rs2");
        assert_eq!(records[1].genotypes, vec![2, codec::G_MISSING]);
    }
}
