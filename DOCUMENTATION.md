# reigen Documentation

`reigen` is a pure-Rust toolkit for population genomics workflows: format conversion, dataset filtering, DTC kit import, multi-dataset merge, stats, and VCF import/export.

---

## 1. Supported formats and layouts

| Format | Layout | Typical files | Notes |
|---|---|---|---|
| `PACKEDANCESTRYMAP` | SNP-major | `.geno .snp .ind` | Binary AdmixTools-style |
| `EIGENSTRAT` | SNP-major | `.geno .snp .ind` | Text `0/1/2/9` |
| `PACKEDPED` | SNP-major | `.bed .bim .fam` | PLINK BED family |
| `TGENO` | Sample-major | `.geno .snp .ind` | Binary sample-major |
| `VCF` | SNP-major | `.vcf` | `export` writes VCF, `vcfimport` reads biallelic SNPs |

`reigen` infers input genotype format from file magic/extension.

---

## 2. Shared path conventions

Most subcommands accept either:

1. `-i/--in-prefix <prefix>` (derive genotype/SNP/IND paths), or
2. explicit files (`--in-geno`, `--in-snp`, `--in-ind`) where supported.

Most writers accept `-o/--out-prefix <prefix>` and derive output extensions from `--out-format`.

---

## 3. Subcommands

### 3.1 `convert`
Convert between supported genotype formats.

```bash
reigen convert -i input --out-format PACKEDPED -o out
```

Common filters:
- sample filters: `--poplist`, `--keep`, `--remove`
- SNP filters: `--snps` (aliases: `--extract`, `--snplist`), `--badsnp` (alias: `--exclude`)
- region filters: `--chrom`, `--from-bp`, `--to-bp`, `--no-xdata`
- QC filters: `--maf`, `--max-maf`, `--hwe`, `--max-miss-snp` (alias: `--geno`), `--mind`

### 3.2 `filter`
Subset/filter a dataset (same filtering surface as `convert`).

```bash
reigen filter -i input -o out --chrom 1-5 --maf 0.01
```

Notes:
- `--out-format` is optional; if omitted, output format defaults to inferred input format.
- Input file flags are `--in-geno`, `--in-snp`, `--in-ind`.

### 3.3 `stats`
Compute per-SNP/per-sample stats and optional IBS matrix.

```bash
reigen stats -i input -o report --ibs
```

Outputs:
- `<prefix>.snp_stats.tsv` (unless `--no-per-snp`)
- `<prefix>.sample_stats.tsv` (unless `--no-per-sample`)
- `<prefix>.ibs.tsv` and `<prefix>.dst.tsv` (with `--ibs`)

### 3.4 `export`
Export dataset to VCF.

```bash
reigen export -i input -o out.vcf --chr-prefix
```

### 3.5 `vcfimport`
Import biallelic SNP VCF into internal formats.

```bash
reigen vcfimport --in in.vcf --out-format PACKEDANCESTRYMAP -o out
```

Notes:
- skips multi-allelic records, indels, and records without `GT`
- optional filters: `--ref-snp`, `--snplist`
- `--numchrom` must be `<= 251`

### 3.6 `import`
Import DTC raw kits (23andMe/Ancestry/MyHeritage/LivingDNA/FTDNA).

```bash
reigen import --in kit.txt --sample-id S1 --out-format PACKEDPED -o s1
```

Highlights:
- vendor auto-detection (or force with `--vendor`)
- supports compressed inputs used by vendor exports
- optional SNP ID filter via `--snps` (alias: `--snplist`)

### 3.7 `merge`
Merge 2+ datasets with allele reconciliation.

```bash
reigen merge --in ds1 --in ds2 --out-format PACKEDPED -o merged
```

Useful switches:
- `--intersection` (otherwise union semantics)
- `--flip-strand`
- `--allow-ambiguous`
- `--strict-ids`

---

## 4. Sample/SNP filtering behavior

### Sample ID matching (`--keep`/`--remove`)
- supports one-column IID lists and two-column `FID IID` lists
- with default family-name handling, PLINK `.fam` IDs can be represented internally as `FID:IID`
- matching handles this `FID:IID` form for two-column lists

### Chromosome coding
Internal mapping is:
- autosomes: `1..numchrom`
- `X = numchrom + 1`
- `Y = numchrom + 2`
- `MT = numchrom + 3`
- `XY = numchrom + 4`

---

## 5. Performance characteristics

- same-layout streaming operations are record-oriented and memory efficient
- cross-layout transforms require transpose/materialization steps
- some workflows on sample-major inputs (notably parts of `export`/`stats`) may hold filtered matrices in memory for correctness

---

## 6. Operational tips

- start with `reigen <subcommand> --help` for current flag contract
- use `--verbose` to inspect filter counts and data-flow decisions
- prefer explicit `--numchrom` for non-human builds; keep `<= 251` for VCF import
