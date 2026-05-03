# reigen

[![Rust](https://img.shields.io/badge/language-Rust-orange.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](#license)

`reigen` is a high-performance population genomics toolkit for format conversion, DTC kit import, dataset merging, and summary statistics. Written in pure Rust, it provides a fast, memory-efficient alternative to legacy tools like AdmixTools `convertf` and `mergeit`, optimized for modern bioinformatics pipelines and large-scale datasets (e.g., AADR).

## Key Features

- **High-Performance**: Pure Rust implementation with streaming I/O and parallel matrix transpose.
- **AADR-Scale**: Optimized for datasets with millions of SNPs and thousands of samples.
- **Memory Efficient**: Processes data record-by-record with a minimal memory footprint.
- **Pure Rust**: No external dependencies on legacy C libraries or AdmixTools.
- **Strand-Aware Merging**: Automatic allele alignment and strand reconciliation when merging datasets.
- **Format Sniffing**: Automatically detects input formats (PAM, EIGENSTRAT, TGENO, BED) from file magic and extensions.

---

## Quick Start

### 1. Convert Formats
```bash
reigen convert \
    --in-prefix input_data \
    --out-format PACKEDANCESTRYMAP \
    --out-prefix converted_data
```

### 2. Import DTC Kit (23andMe, Ancestry, etc.)
```bash
reigen import \
    --in raw_kit.txt \
    --sample-id my_sample \
    --out-prefix my_kit
```

### 3. Merge Datasets
```bash
reigen merge \
    --in dataset1 \
    --in dataset2 \
    --out-prefix merged_output
```

### 4. Compute Statistics
```bash
reigen stats -i input_prefix -o stats_report
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

## Documentation
For a detailed guide on all subcommands, filters, and format specifications, see [DOCUMENTATION.md](DOCUMENTATION.md).

## License
Licensed under the [GNU General Public License v3.0](LICENSE).
