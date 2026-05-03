# reigen

[![Rust](https://img.shields.io/badge/language-Rust-orange.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](#license)

`reigen` is a high-performance population genomics toolkit for format conversion, DTC kit import, dataset merging, and summary statistics. Written in pure Rust, it provides a fast, memory-efficient alternative to legacy tools like AdmixTools `convertf` and `mergeit`, optimized for modern bioinformatics pipelines and large-scale datasets (e.g., AADR).

## Key Features

- **High-Performance**: Pure Rust implementation with streaming I/O and parallel matrix transpose.
- **AADR-Scale**: Optimized for datasets with millions of SNPs and thousands of samples.
- **Pure Rust**: No external dependencies on legacy C libraries or AdmixTools.
- **Memory Efficient**: Processes data record-by-record with a minimal memory footprint.
- **Strand-Aware**: Integrated allele alignment and strand reconciliation for kit imports and merging.
- **Format Sniffing**: Automatically detects input formats (PAM, EIGENSTRAT, TGENO, BED) from file magic and extensions.

---

## Subcommands

### 1. `convert` — Format Conversion
Lossless conversion between major population genetics formats.

```bash
reigen convert \
    --in-prefix input_data \
    --out-format PACKEDANCESTRYMAP \
    --out-prefix converted_data
```

**Common Flags:**
- `-i, --in-prefix`: Input prefix (derives `.geno`/`.snp`/`.ind` or `.bed`/`.bim`/`.fam`).
- `-o, --out-prefix`: Output prefix.
- `--out-format`: Target format (see [Supported Formats](#supported-formats)).
- `--max-miss-snp <float>`: Filter SNPs by maximum missingness (alias: `--geno`).
- `--max-miss-ind <float>`: Filter samples by maximum missingness (alias: `--mind`).
- `--maf <float>`: Filter SNPs by minimum Minor Allele Frequency.

### 2. `import` — DTC Kit Import
Directly import raw data from consumer genomics vendors (23andMe, AncestryDNA, MyHeritage, LivingDNA, FTDNA) and align them to a reference.

```bash
reigen import \
    --in raw_kit.txt \
    --ref-snp reference.snp \
    --sample-id my_sample \
    --out-prefix my_aligned_kit
```

### 3. `merge` — Dataset Merging
Merge multiple datasets with automatic allele alignment and strand reconciliation.

```bash
reigen merge \
    --in dataset1 \
    --in dataset2 \
    --out-prefix merged_output \
    --intersection # Optional: keep only shared SNPs
```

### 4. `filter` — Subset & Filter
Efficiently subset a dataset by population, sample, chromosome, or genomic range.

```bash
reigen filter \
    -i input_prefix \
    -o filtered_output \
    --poplist populations.txt \
    --chrom 1-5,22
```
*Note: `--out-format` is optional and defaults to the input format.*

### 5. `stats` — Summary Statistics
Compute per-SNP and per-sample QC statistics in a single pass.

```bash
reigen stats -i input_prefix -o stats_report
```
**Outputs:**
- `<prefix>.snp_stats.tsv`: MAF, HWE p-values, missingness, and genotype counts.
- `<prefix>.sample_stats.tsv`: Heterozygosity, missingness, and call rates.
- `<prefix>.ibs.tsv`: Pairwise IBS distance matrix (optional with `--ibs`).

---

## Supported Formats

| Format | Extension | Type | Layout |
|--------|-----------|------|--------|
| **PACKEDANCESTRYMAP** | `.geno` | Binary | SNP-Major |
| **EIGENSTRAT** | `.geno` | Text | SNP-Major |
| **PACKEDPED** | `.bed` | Binary | SNP-Major |
| **TGENO** | `.tgeno` | Binary | Sample-Major |
| **VCF** (Export only) | `.vcf` | Text | SNP-Major |

---

## Performance & Design

`reigen` uses a **streaming architecture** to process data. Same-layout conversions (e.g., PAM → BED) are processed record-by-record, using negligible RAM even for the largest datasets. Cross-layout conversions (e.g., PAM → TGENO) utilize a parallelized blocked-matrix transpose.

- **Matrix Transpose**: Multi-threaded using `rayon`, capable of transposing AADR-scale matrices (1.2M SNPs × 17k samples) in seconds.
- **Safety**: Fully implemented in safe Rust with checked arithmetic for chromosome mappings and type-safe metadata handling.

---

## Workflow Examples

### Convert PLINK to EIGENSTRAT with Filtering
```bash
reigen convert \
    --in-prefix my_plink_data \
    --out-format EIGENSTRAT \
    --out-prefix out_eigen \
    --geno 0.05 \
    --maf 0.01
```

### Compute IBS Matrix for a Population
```bash
reigen stats \
    -i dataset \
    -o pop_stats \
    --poplist my_pop.txt \
    --ibs
```

---

## Installation

```bash
# Clone and build
git clone https://github.com/youruser/reigen.git
cd reigen
cargo build --release

# Run
./target/release/reigen --help
```

## Testing
`reigen` includes an extensive test suite, including bit-level roundtrip verification and synthetic smoke tests.
```bash
cargo test
```

## License
Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at your option.
