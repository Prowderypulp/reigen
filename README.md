# reigen

[![Rust](https://img.shields.io/badge/language-Rust-orange.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](#license)

`reigen` is a pure-Rust population genomics toolkit for format conversion, DTC kit import, dataset merge/reconciliation, filtering, statistics, and VCF interoperability.

## Key Features

- **CLI toolkit**: `convert`, `import`, `merge`, `filter`, `stats`, `export`, `vcfimport`
- **Format sniffing**: Detects PAM/EIGENSTRAT/TGENO/BED from magic + extensions
- **Cross-format conversion**: AdmixTools and PLINK families with shared filters
- **Strand-aware merge**: Allele reconciliation during multi-dataset merge
- **VCF support**: VCF export plus biallelic SNP VCF import
- **Streaming-first core**: Same-layout operations stream efficiently; some sample-major workflows materialize matrices

---

## Quick Start

### 1. Convert formats
```bash
reigen convert \
    -i input_data \
    --out-format PACKEDANCESTRYMAP \
    -o converted_data
```

### 2. Filter a dataset
```bash
reigen filter \
    -i input_data \
    -o filtered_data \
    --chrom 1-5 \
    --maf 0.01 \
    --keep keep_samples.txt
```

### 3. Compute stats
```bash
reigen stats -i input_data -o qc --ibs
```

### 4. Export to VCF
```bash
reigen export -i input_data -o output.vcf --chr-prefix
```

### 5. Import from VCF
```bash
reigen vcfimport \
    --in input.vcf \
    --out-format PACKEDPED \
    -o vcf_as_bed
```

## Important CLI notes

- `convert`/`filter` input file flags are `--in-geno`, `--in-snp`, `--in-ind` (or `-i/--in-prefix`)
- `filter --out-format` is optional; when omitted it defaults to inferred input format
- missingness aliases:
  - per-SNP: `--max-miss-snp` (alias `--geno`)
  - per-sample: `--mind`
- SNP keep-list aliases on `convert`/`filter`: `--snps`, `--extract`, `--snplist`
- `stats` toggles are negative flags: `--no-per-snp`, `--no-per-sample`
- `vcfimport` only supports biallelic SNP records and enforces `--numchrom <= 251`
- `merge` accepts repeated `--in` and `--in-list <file>` input specs; it writes `.missnp*` drop reports and `.idmap.tsv` when IDs are auto-renamed

## Installation (from source)

```bash
cargo build --release
./target/release/reigen --help
```

## Documentation
For a detailed guide on all subcommands, filters, and format specifications, see [DOCUMENTATION.md](DOCUMENTATION.md).

## License
Licensed under the [GNU General Public License v3.0](LICENSE).
