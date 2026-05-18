use crate::format::Format;
use crate::geno::codec;
use crate::meta::{IndRow, Sex, SnpRow};
use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

pub mod ancestry;
pub mod detect;
pub mod ftdna;
pub mod livingdna;
pub mod myheritage;
pub mod twenty_three;
pub mod vendor;

use self::detect::Vendor;
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

    /// SNP ID list (one rsID per line). Filter output to these SNPs only.
    #[arg(long, alias = "snplist")]
    pub snps: Option<PathBuf>,

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

    /// Include X/Y/MT calls (default: autosomes only).
    #[arg(long)]
    pub include_non_autosomal: bool,

    /// Number of autosomes (default 22).
    #[arg(long, default_value_t = 22)]
    pub numchrom: u32,
}

pub struct ImportConfig {
    pub input: PathBuf,
    pub snps: Option<PathBuf>,
    pub sample_id: String,
    pub sample_pop: String,
    pub sample_sex: Sex,
    pub out_format: Format,
    pub out_geno: PathBuf,
    pub out_snp: PathBuf,
    pub out_ind: PathBuf,
    pub vendor: Option<Vendor>,
    pub include_non_autosomal: bool,
    pub numchrom: u32,
}

#[derive(Default)]
struct ImportStats {
    total_calls: usize,
    matched: usize,
    dropped_nocall: usize,
    dropped_not_in_snps: usize,
    dropped_non_autosomal: usize,
    dropped_bad_chrom: usize,
}

pub fn run_import(args: ImportArgs) -> Result<()> {
    let out_format = args.out_format;

    let (geno_in, snp_in, ind_in) = crate::pipeline::resolve_paths(
        args.out_prefix.clone(),
        args.out_geno,
        args.out_snp,
        args.out_ind,
        Some(out_format),
        true,
    )?;

    let cfg = ImportConfig {
        input: args.input,
        snps: args.snps,
        sample_id: args.sample_id.clone(),
        sample_pop: args.sample_pop.unwrap_or_else(|| args.sample_id.clone()),
        sample_sex: match args.sample_sex {
            SexArg::M => Sex::Male,
            SexArg::F => Sex::Female,
            SexArg::U => Sex::Unknown,
        },
        out_format,
        out_geno: geno_in,
        out_snp: snp_in,
        out_ind: ind_in,
        vendor: match args.vendor {
            VendorArg::Auto => None,
            VendorArg::TwentyThreeAndMe => Some(Vendor::TwentyThreeAndMe),
            VendorArg::Ancestry => Some(Vendor::Ancestry),
            VendorArg::MyHeritage => Some(Vendor::MyHeritage),
            VendorArg::LivingDna => Some(Vendor::LivingDna),
            VendorArg::Ftdna => Some(Vendor::Ftdna),
        },
        include_non_autosomal: args.include_non_autosomal,
        numchrom: args.numchrom,
    };

    // 1. Detect vendor and load lines
    let kit_lines = read_kit_lines(&cfg.input)?;
    let sniff_text = kit_lines[..kit_lines.len().min(100)].join("\n");

    let vendor = match cfg.vendor {
        Some(v) => v,
        None => detect::detect_vendor_from_text(&sniff_text, Some(&cfg.input))?,
    };
    log::info!("detected vendor: {:?}", vendor);

    let mut parser = get_parser(vendor);

    // 2. Process calls
    let snp_keep = cfg
        .snps
        .as_deref()
        .map(crate::filter::load_snp_keep)
        .transpose()?;

    log::info!("parsing kit {}...", cfg.input.display());
    let mut stats = ImportStats::default();

    let mut out_snps: Vec<SnpRow> = Vec::new();
    let mut out_genos: Vec<u8> = Vec::new();

    for line in &kit_lines {
        if let Some(call) = parser.parse_line(line)? {
            stats.total_calls += 1;

            if let Some(ref ids) = snp_keep {
                if !ids.contains(&call.rsid) {
                    stats.dropped_not_in_snps += 1;
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
            if !matches!(a, b'A' | b'C' | b'G' | b'T') || !matches!(b, b'A' | b'C' | b'G' | b'T') {
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

    log::info!("Import finished:");
    log::info!("  Total calls:              {}", stats.total_calls);
    log::info!("  Kept:                     {}", stats.matched);
    log::info!("  Dropped (no-call):        {}", stats.dropped_nocall);
    if snp_keep.is_some() {
        log::info!(
            "  Dropped (not in SNPs list): {}",
            stats.dropped_not_in_snps
        );
    }
    log::info!(
        "  Dropped (non-autosomal):  {}",
        stats.dropped_non_autosomal
    );
    log::info!("  Dropped (bad chrom):      {}", stats.dropped_bad_chrom);

    if stats.matched == 0 {
        anyhow::bail!("no valid kit calls in output; output would be empty");
    }

    log::info!("  Output SNPs:              {}", out_snps.len());

    // 3. Write output
    log::info!("writing output to {}...", cfg.out_geno.display());
    if matches!(cfg.out_format, Format::PackedAncestrymap) {
        log::warn!(
            "PACKEDANCESTRYMAP uses a fixed 48-byte minimum SNP record size; \
             for single-sample imports, TGENO is usually much smaller."
        );
    }

    let mut writer = crate::pipeline::open_writer_pub(cfg.out_format, &cfg.out_geno)?;
    writer.begin(1, out_snps.len(), 0, 0)?;

    match writer.layout() {
        crate::geno::Layout::SnpMajor => {
            // One record per SNP, each record has 1 sample. 1 byte is enough for 1-4 samples.
            let mut buf = vec![0u8; 1];
            for &g in &out_genos {
                codec::pack(&[g], &mut buf);
                writer.write_record(&buf)?;
            }
        }
        crate::geno::Layout::SampleMajor => {
            // One record for the sample, covering all SNPs.
            let mut buf = vec![0u8; (out_snps.len() * 2 + 7) / 8];
            codec::pack(&out_genos, &mut buf);
            writer.write_record(&buf)?;
        }
    }
    writer.finish()?;

    log::info!("writing metadata...");
    crate::pipeline::write_output_snp(&cfg.out_snp, cfg.out_format, &out_snps, cfg.numchrom)?;

    let ind = IndRow {
        id: cfg.sample_id,
        sex: cfg.sample_sex,
        pop: cfg.sample_pop,
        ignore: false,
    };
    crate::pipeline::write_output_ind(&cfg.out_ind, cfg.out_format, &[ind], false)?;

    log::info!("done.");
    Ok(())
}

pub fn get_parser(vendor: Vendor) -> Box<dyn VendorParser> {
    match vendor {
        Vendor::TwentyThreeAndMe => Box::new(twenty_three::TwentyThreeAndMeParser),
        Vendor::Ancestry => Box::new(ancestry::AncestryParser),
        Vendor::MyHeritage => Box::new(myheritage::MyHeritageParser),
        Vendor::LivingDna => Box::new(livingdna::LivingDnaParser),
        Vendor::Ftdna => Box::new(ftdna::FtdnaParser),
    }
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

fn parse_kit_chrom(s: &str, numchrom: u32) -> Option<u8> {
    if let Ok(c) = s.parse::<u8>() {
        return Some(c);
    }
    match s.trim().to_ascii_uppercase().as_str() {
        "X" | "23" => Some(numchrom as u8 + 1),
        "Y" | "24" => Some(numchrom as u8 + 2),
        "MT" | "M" | "25" => Some(numchrom as u8 + 3),
        "XY" => Some(numchrom as u8 + 4),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::{parse_kit_chrom, read_kit_lines};
    use crate::import::ancestry::AncestryParser;
    use crate::import::vendor::VendorParser;
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
}
