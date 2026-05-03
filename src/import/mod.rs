use crate::format::Format;
use crate::geno::{codec, GenoWriter};
use crate::meta::{self, IndRow, Sex, SnpRow};
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

pub mod ancestry;
pub mod detect;
pub mod ftdna;
pub mod livingdna;
pub mod myheritage;
pub mod reconcile;
pub mod twenty_three;
pub mod vendor;

use self::detect::Vendor;
use self::reconcile::ReconcileResult;
use self::vendor::VendorParser;

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum VendorArg {
    Auto,
    TwentyThreeAndMe,
    Ancestry,
    MyHeritage,
    LivingDna,
    Ftdna,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum SexArg {
    M,
    F,
    U,
}

#[derive(Args, Debug)]
pub struct ImportArgs {
    /// Input DTC kit file
    #[arg(long = "in", value_name = "KIT")]
    pub input: PathBuf,

    /// Reference SNP file (.snp or .bim) — defines the output SNP space and enables
    /// strand reconciliation. Omit to import the kit as-is (plink-style).
    #[arg(long)]
    pub ref_snp: Option<PathBuf>,

    /// SNP ID list (one rsID per line). Filter output to these SNPs only.
    /// Works with or without --ref-snp.
    #[arg(long)]
    pub snplist: Option<PathBuf>,

    /// Sample ID written to .ind / .fam.
    #[arg(long)]
    pub sample_id: String,

    /// Population label (default: same as --sample-id).
    #[arg(long)]
    pub sample_pop: Option<String>,

    /// Sample sex.
    #[arg(long, value_enum, default_value_t = SexArg::U)]
    pub sample_sex: SexArg,

    /// Output format.
    #[arg(long)]
    pub out_format: Format,

    /// Output prefix (derives .geno/.snp/.ind or .bed/.bim/.fam).
    #[arg(short = 'o', long)]
    pub out_prefix: Option<String>,

    #[arg(long)]
    pub out_geno: Option<PathBuf>,
    #[arg(long)]
    pub out_snp: Option<PathBuf>,
    #[arg(long)]
    pub out_ind: Option<PathBuf>,

    /// Force a specific vendor (default: auto-detect from header).
    #[arg(long, value_enum, default_value_t = VendorArg::Auto)]
    pub vendor: VendorArg,

    /// Keep A/T and C/G ambiguous SNPs (default: drop).
    #[arg(long)]
    pub allow_ambiguous: bool,

    /// Disable complement-match strand reconciliation. Default is on for
    /// imports because DTC files don't report strand reliably.
    #[arg(long)]
    pub no_flip_strand: bool,

    /// Include X/Y/MT calls (default: autosomes only).
    #[arg(long)]
    pub include_non_autosomal: bool,

    /// Number of autosomes (default 22).
    #[arg(long, default_value_t = 22)]
    pub numchrom: u32,

    /// Keep the full reference SNP space in output, filling unmatched SNPs as
    /// missing. Default behavior is matched-only output.
    #[arg(long)]
    pub full_reference: bool,
}

pub struct ImportConfig {
    pub input: PathBuf,
    pub ref_snp: Option<PathBuf>,
    pub snplist: Option<PathBuf>,
    pub sample_id: String,
    pub sample_pop: String,
    pub sample_sex: Sex,
    pub out_format: Format,
    pub out_geno: PathBuf,
    pub out_snp: PathBuf,
    pub out_ind: PathBuf,
    pub vendor: Option<Vendor>,
    pub allow_ambiguous: bool,
    pub flip_strand: bool,
    pub include_non_autosomal: bool,
    pub numchrom: u32,
    pub full_reference: bool,
}

#[derive(Default)]
struct ImportStats {
    total_calls: usize,
    matched: usize,
    dropped_mismatch: usize,
    dropped_ambiguous: usize,
    dropped_non_autosomal: usize,
    dropped_nocall: usize,
    dropped_strand: usize,
    dropped_not_in_ref: usize,
    dropped_bad_chrom: usize,
}

pub fn run_import(args: ImportArgs) -> Result<()> {
    let p_path = args.out_prefix.as_ref().map(PathBuf::from);
    let (geno_ext, snp_ext, ind_ext) = args.out_format.default_output_extensions();
    let out_geno = args
        .out_geno
        .or_else(|| p_path.as_ref().map(|p| p.with_extension(geno_ext)))
        .ok_or_else(|| anyhow::anyhow!("missing --out-geno or --out-prefix"))?;
    let out_snp = args
        .out_snp
        .or_else(|| p_path.as_ref().map(|p| p.with_extension(snp_ext)))
        .ok_or_else(|| anyhow::anyhow!("missing --out-snp or --out-prefix"))?;
    let out_ind = args
        .out_ind
        .or_else(|| p_path.as_ref().map(|p| p.with_extension(ind_ext)))
        .ok_or_else(|| anyhow::anyhow!("missing --out-ind or --out-prefix"))?;

    let vendor = match args.vendor {
        VendorArg::Auto => None,
        VendorArg::TwentyThreeAndMe => Some(Vendor::TwentyThreeAndMe),
        VendorArg::Ancestry => Some(Vendor::Ancestry),
        VendorArg::MyHeritage => Some(Vendor::MyHeritage),
        VendorArg::LivingDna => Some(Vendor::LivingDna),
        VendorArg::Ftdna => Some(Vendor::Ftdna),
    };

    let sex = match args.sample_sex {
        SexArg::M => Sex::Male,
        SexArg::F => Sex::Female,
        SexArg::U => Sex::Unknown,
    };

    if args.full_reference && args.ref_snp.is_none() {
        anyhow::bail!("--full-reference requires --ref-snp");
    }

    let cfg = ImportConfig {
        input: args.input,
        ref_snp: args.ref_snp,
        snplist: args.snplist,
        sample_id: args.sample_id.clone(),
        sample_pop: args.sample_pop.unwrap_or(args.sample_id),
        sample_sex: sex,
        out_format: args.out_format,
        out_geno,
        out_snp,
        out_ind,
        vendor,
        allow_ambiguous: args.allow_ambiguous,
        flip_strand: !args.no_flip_strand,
        include_non_autosomal: args.include_non_autosomal,
        numchrom: args.numchrom,
        full_reference: args.full_reference,
    };

    let kit_lines = read_kit_lines(&cfg.input)?;
    let sniff_text = kit_lines.join("\n");

    let vendor = if let Some(v) = cfg.vendor {
        v
    } else {
        let v = detect::detect_vendor_from_text(&sniff_text, Some(&cfg.input))?;
        log::info!("detected vendor: {:?}", v);
        v
    };

    let mut parser: Box<dyn VendorParser> = match vendor {
        Vendor::TwentyThreeAndMe => Box::new(twenty_three::TwentyThreeAndMeParser),
        Vendor::Ancestry => Box::new(ancestry::AncestryParser),
        Vendor::MyHeritage => Box::new(myheritage::MyHeritageParser),
        Vendor::LivingDna => Box::new(livingdna::LivingDnaParser),
        Vendor::Ftdna => Box::new(ftdna::FtdnaParser),
    };

    // Optional snplist filter (rsIDs, one per line)
    let snplist_ids: Option<HashSet<String>> = if let Some(ref p) = cfg.snplist {
        log::info!("loading snplist from {}...", p.display());
        let ids = load_snplist(p)?;
        log::info!("loaded {} SNP IDs", ids.len());
        Some(ids)
    } else {
        None
    };

    log::info!("parsing kit {}...", cfg.input.display());
    let mut stats = ImportStats::default();

    let (out_snps, out_genos) = if let Some(ref ref_path) = cfg.ref_snp {
        // --- ref-guided path: match by position, reconcile strand ---
        log::info!("loading reference SNPs from {}...", ref_path.display());
        let ref_fmt = if ref_path.extension().and_then(|e| e.to_str()) == Some("bim") {
            Format::PackedPed
        } else {
            Format::Eigenstrat
        };
        let ref_snps = match ref_fmt {
            Format::PackedPed => meta::bim::read(ref_path, cfg.numchrom)?,
            _ => meta::snp::read(ref_path, cfg.numchrom)?,
        };
        log::info!("loaded {} reference SNPs", ref_snps.len());

        let mut ref_map: HashMap<(u8, u64), usize> = HashMap::with_capacity(ref_snps.len());
        for (i, s) in ref_snps.iter().enumerate() {
            ref_map.insert((s.chrom, s.physical_pos), i);
        }

        let mut genos = vec![codec::G_MISSING; ref_snps.len()];
        let mut matched_ref_mask = vec![false; ref_snps.len()];

        for line in &kit_lines {
            if let Some(call) = parser.parse_line(line)? {
                stats.total_calls += 1;

                if let Some(ref ids) = snplist_ids {
                    if !ids.contains(&call.rsid) {
                        stats.dropped_not_in_ref += 1;
                        continue;
                    }
                }

                let Some(chrom_u8) = parse_kit_chrom(&call.chrom, cfg.numchrom) else {
                    stats.dropped_bad_chrom += 1;
                    continue;
                };

                if chrom_u8 > cfg.numchrom as u8 && !cfg.include_non_autosomal {
                    stats.dropped_non_autosomal += 1;
                    continue;
                }

                if let Some(&ref_idx) = ref_map.get(&(chrom_u8, call.pos)) {
                    let ref_snp = &ref_snps[ref_idx];
                    match reconcile::reconcile_call(
                        &call,
                        ref_snp,
                        cfg.allow_ambiguous,
                        cfg.flip_strand,
                    ) {
                        ReconcileResult::Encoded(g) => {
                            stats.matched += 1;
                            genos[ref_idx] = g;
                            matched_ref_mask[ref_idx] = true;
                        }
                        ReconcileResult::DropAmbiguous => stats.dropped_ambiguous += 1,
                        ReconcileResult::DropMismatch => stats.dropped_mismatch += 1,
                        ReconcileResult::DropNoCall => stats.dropped_nocall += 1,
                        ReconcileResult::DropStrandUnresolvable => stats.dropped_strand += 1,
                    }
                } else {
                    stats.dropped_not_in_ref += 1;
                }
            }
        }

        log::info!("Import finished:");
        log::info!("  Total calls:              {}", stats.total_calls);
        log::info!("  Matched to ref:           {}", stats.matched);
        log::info!("  Dropped (mismatch):       {}", stats.dropped_mismatch);
        log::info!("  Dropped (strand):         {}", stats.dropped_strand);
        log::info!("  Dropped (ambiguous A/T,C/G): {}", stats.dropped_ambiguous);
        log::info!("  Dropped (no-call):        {}", stats.dropped_nocall);
        log::info!("  Dropped (not in ref):     {}", stats.dropped_not_in_ref);
        log::info!("  Dropped (non-autosomal):  {}", stats.dropped_non_autosomal);
        log::info!("  Dropped (bad chrom):      {}", stats.dropped_bad_chrom);
        log::info!(
            "  Output SNP mode:          {}",
            if cfg.full_reference { "full-reference" } else { "matched-only" }
        );

        if stats.matched == 0 {
            anyhow::bail!("no kit variants matched the reference; output would be empty");
        }

        select_output_panel(&ref_snps, &genos, &matched_ref_mask, cfg.full_reference)
    } else {
        // --- no-ref path: import kit as-is, build SNP metadata from calls ---
        let mut out_snps: Vec<SnpRow> = Vec::new();
        let mut out_genos: Vec<u8> = Vec::new();

        for line in &kit_lines {
            if let Some(call) = parser.parse_line(line)? {
                stats.total_calls += 1;

                if let Some(ref ids) = snplist_ids {
                    if !ids.contains(&call.rsid) {
                        stats.dropped_not_in_ref += 1;
                        continue;
                    }
                }

                let Some(chrom_u8) = parse_kit_chrom(&call.chrom, cfg.numchrom) else {
                    stats.dropped_bad_chrom += 1;
                    continue;
                };

                if chrom_u8 > cfg.numchrom as u8 && !cfg.include_non_autosomal {
                    stats.dropped_non_autosomal += 1;
                    continue;
                }

                let a = call.a1.to_ascii_uppercase() as u8;
                let b = call.a2.to_ascii_uppercase() as u8;
                if !matches!(a, b'A' | b'C' | b'G' | b'T')
                    || !matches!(b, b'A' | b'C' | b'G' | b'T')
                {
                    stats.dropped_nocall += 1;
                    continue;
                }

                // Hom: allele1='0' (unknown ref), allele2=observed, geno=2.
                // Het: sort alphabetically so orientation is deterministic; geno=1.
                let (allele1, allele2, geno) = if a == b {
                    (b'0', a, 2u8)
                } else {
                    let (al1, al2) = if a < b { (a, b) } else { (b, a) };
                    (al1, al2, 1u8)
                };

                out_snps.push(SnpRow {
                    id: call.rsid.clone(),
                    chrom: chrom_u8,
                    genetic_pos: 0.0,
                    physical_pos: call.pos,
                    allele1,
                    allele2,
                });
                out_genos.push(geno);
                stats.matched += 1;
            }
        }

        log::info!("Import finished (no-ref mode):");
        log::info!("  Total calls:              {}", stats.total_calls);
        log::info!("  Kept:                     {}", stats.matched);
        log::info!("  Dropped (no-call):        {}", stats.dropped_nocall);
        if snplist_ids.is_some() {
            log::info!("  Dropped (not in snplist): {}", stats.dropped_not_in_ref);
        }
        log::info!("  Dropped (non-autosomal):  {}", stats.dropped_non_autosomal);
        log::info!("  Dropped (bad chrom):      {}", stats.dropped_bad_chrom);

        if stats.matched == 0 {
            anyhow::bail!("no valid kit calls in output; output would be empty");
        }

        (out_snps, out_genos)
    };

    log::info!("  Output SNPs:              {}", out_snps.len());

    // 3. Write output
    log::info!("writing output to {}...", cfg.out_geno.display());
    if matches!(cfg.out_format, Format::PackedAncestrymap) {
        log::warn!(
            "PACKEDANCESTRYMAP uses a fixed 48-byte minimum SNP record size; \
             for single-sample imports, TGENO is usually much smaller."
        );
    }

    // PACKEDANCESTRYMAP needs a hash header
    let ind_ids = vec![cfg.sample_id.as_str()];
    let snp_ids: Vec<&str> = out_snps.iter().map(|s| s.id.as_str()).collect();
    let ihash = crate::hash::hasharr(&ind_ids);
    let shash = crate::hash::hasharr(&snp_ids);

    let mut writer: Box<dyn GenoWriter> = match cfg.out_format {
        Format::PackedAncestrymap => Box::new(crate::geno::packed_am::PackedAmWriter::create(
            &cfg.out_geno,
        )?),
        Format::Eigenstrat => Box::new(crate::geno::eigenstrat::EigenstratWriter::create(
            &cfg.out_geno,
        )?),
        Format::PackedPed => Box::new(crate::geno::packed_ped::PackedPedWriter::create(
            &cfg.out_geno,
        )?),
        Format::Tgeno => Box::new(crate::geno::tgeno::TgenoWriter::create(&cfg.out_geno)?),
        _ => anyhow::bail!("output format not supported for import"),
    };

    writer.begin(1, out_snps.len(), ihash, shash)?;

    if matches!(cfg.out_format, Format::Tgeno) {
        // TGENO is sample-major: one record per sample, each record spans all SNPs.
        let mut record = vec![0u8; (out_snps.len() * 2 + 7) / 8];
        codec::pack(&out_genos, &mut record);
        writer.write_record(&record)?;
    } else {
        // SNP-major outputs: one record per SNP (here one sample per record).
        let mut record = vec![0u8; 1];
        for &g in &out_genos {
            codec::pack(&[g], &mut record);
            writer.write_record(&record)?;
        }
    }
    writer.finish()?;

    // Metadata
    log::info!("writing metadata...");
    match cfg.out_format {
        Format::PackedPed => {
            meta::bim::write(&cfg.out_snp, &out_snps, cfg.numchrom)?;
            meta::fam::write(
                &cfg.out_ind,
                &[IndRow {
                    id: cfg.sample_id,
                    pop: cfg.sample_pop,
                    sex: cfg.sample_sex,
                    ignore: false,
                }],
                false,
            )?;
        }
        _ => {
            meta::snp::write(&cfg.out_snp, &out_snps, cfg.numchrom)?;
            meta::ind::write(
                &cfg.out_ind,
                &[IndRow {
                    id: cfg.sample_id,
                    pop: cfg.sample_pop,
                    sex: cfg.sample_sex,
                    ignore: false,
                }],
            )?;
        }
    }

    log::info!("done.");
    Ok(())
}

fn load_snplist(path: &Path) -> Result<HashSet<String>> {
    let f = File::open(path).with_context(|| format!("open snplist {}", path.display()))?;
    let reader = BufReader::new(f);
    let mut ids = HashSet::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("read {}", path.display()))?;
        let id = line.trim();
        if !id.is_empty() && !id.starts_with('#') {
            ids.insert(id.to_owned());
        }
    }
    Ok(ids)
}

fn read_kit_lines(path: &Path) -> Result<Vec<String>> {
    let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut magic = [0u8; 4];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    if n >= 2 && magic[..2] == [0x1f, 0x8b] {
        let reader = BufReader::new(flate2::read::MultiGzDecoder::new(file));
        return collect_lines(reader, path);
    }

    if n >= 4 && (magic == *b"PK\x03\x04" || magic == *b"PK\x05\x06" || magic == *b"PK\x07\x08") {
        let mut archive =
            zip::ZipArchive::new(file).with_context(|| format!("open zip {}", path.display()))?;
        let mut regular_indices = Vec::new();
        let mut names = Vec::new();
        for i in 0..archive.len() {
            let f = archive.by_index(i)?;
            if f.is_file() {
                regular_indices.push(i);
                names.push(f.name().to_string());
            }
        }
        if regular_indices.len() != 1 {
            anyhow::bail!(
                "{}: zip import requires exactly one regular file entry; found {} ({})",
                path.display(),
                regular_indices.len(),
                names.join(", ")
            );
        }
        let mut entry = archive.by_index(regular_indices[0])?;
        let mut buf = String::new();
        entry
            .read_to_string(&mut buf)
            .with_context(|| format!("read zip entry {}", entry.name()))?;
        return Ok(buf.lines().map(ToOwned::to_owned).collect());
    }

    collect_lines(BufReader::new(file), path)
}

fn collect_lines<R: BufRead>(reader: R, path: &Path) -> Result<Vec<String>> {
    reader
        .lines()
        .collect::<std::io::Result<Vec<_>>>()
        .with_context(|| format!("read lines from {}", path.display()))
}

fn select_output_panel(
    ref_snps: &[crate::meta::SnpRow],
    genos: &[u8],
    matched_ref_mask: &[bool],
    full_reference: bool,
) -> (Vec<crate::meta::SnpRow>, Vec<u8>) {
    if full_reference {
        return (ref_snps.to_vec(), genos.to_vec());
    }
    let mut snps = Vec::new();
    let mut out_genos = Vec::new();
    for ((s, &geno), &matched) in ref_snps
        .iter()
        .zip(genos.iter())
        .zip(matched_ref_mask.iter())
    {
        if matched {
            snps.push(s.clone());
            out_genos.push(geno);
        }
    }
    (snps, out_genos)
}

fn parse_kit_chrom(chrom: &str, numchrom: u32) -> Option<u8> {
    if let Ok(c) = chrom.parse::<u8>() {
        return Some(c);
    }
    match chrom.trim().to_ascii_uppercase().as_str() {
        "X" | "23" => Some(numchrom as u8 + 1),
        "Y" | "24" => Some(numchrom as u8 + 2),
        "MT" | "M" | "25" => Some(numchrom as u8 + 3),
        "XY" => Some(numchrom as u8 + 4),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::{parse_kit_chrom, read_kit_lines, select_output_panel};
    use crate::import::ancestry::AncestryParser;
    use crate::import::vendor::VendorParser;
    use crate::meta::SnpRow;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn parses_xy_chromosome() {
        assert_eq!(parse_kit_chrom("XY", 22), Some(26));
    }

    #[test]
    fn ancestry_parser_skips_unprefixed_header() {
        let mut p = AncestryParser;
        let row = p
            .parse_line("rsid\tchromosome\tposition\tallele1\tallele2")
            .unwrap();
        assert!(row.is_none());
    }

    #[test]
    fn reads_gzip_input_lines() {
        let f = tempfile::NamedTempFile::new().unwrap();
        let out = File::create(f.path()).unwrap();
        let mut enc = flate2::write::GzEncoder::new(out, flate2::Compression::default());
        writeln!(enc, "#hdr").unwrap();
        writeln!(enc, "rs1\t1\t123\tAG").unwrap();
        enc.finish().unwrap();

        let lines = read_kit_lines(f.path()).unwrap();
        assert_eq!(lines.len(), 2);
        assert!(lines[1].contains("rs1"));
    }

    #[test]
    fn reads_single_file_zip_input_lines() {
        let f = tempfile::NamedTempFile::new().unwrap();
        {
            let out = File::create(f.path()).unwrap();
            let mut zw = zip::ZipWriter::new(out);
            let opts = zip::write::FileOptions::default();
            zw.start_file("kit.txt", opts).unwrap();
            writeln!(zw, "#hdr").unwrap();
            writeln!(zw, "rs1\t1\t123\tAG").unwrap();
            zw.finish().unwrap();
        }

        let lines = read_kit_lines(f.path()).unwrap();
        assert_eq!(lines.len(), 2);
        assert!(lines[1].contains("rs1"));
    }

    #[test]
    fn matched_only_output_excludes_unmatched_reference_snps() {
        let ref_snps = vec![
            SnpRow {
                id: "rs1".into(),
                chrom: 1,
                genetic_pos: 0.0,
                physical_pos: 10,
                allele1: b'A',
                allele2: b'G',
            },
            SnpRow {
                id: "rs2".into(),
                chrom: 1,
                genetic_pos: 0.0,
                physical_pos: 20,
                allele1: b'C',
                allele2: b'T',
            },
        ];
        let genos = vec![2u8, 3u8];
        let matched = vec![true, false];
        let (snps, out_genos) = select_output_panel(&ref_snps, &genos, &matched, false);
        assert_eq!(snps.len(), 1);
        assert_eq!(snps[0].id, "rs1");
        assert_eq!(out_genos, vec![2u8]);
    }

    #[test]
    fn full_reference_output_keeps_unmatched_as_missing() {
        let ref_snps = vec![SnpRow {
            id: "rs1".into(),
            chrom: 1,
            genetic_pos: 0.0,
            physical_pos: 10,
            allele1: b'A',
            allele2: b'G',
        }];
        let genos = vec![3u8];
        let matched = vec![false];
        let (snps, out_genos) = select_output_panel(&ref_snps, &genos, &matched, true);
        assert_eq!(snps.len(), 1);
        assert_eq!(out_genos, vec![3u8]);
    }
}
