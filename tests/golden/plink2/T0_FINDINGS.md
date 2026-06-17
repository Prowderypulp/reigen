# T0 spike — findings (PGEN format facts)

Purpose: empirically pin the two format facts the PGEN codec/writer depend on,
before any reader/writer code is written. Fixture: `t0_lsb.{pgen,pvar,psam}`,
regenerated only by `t0_genesis.sh` (needs real `plink2`). Verified 2026-06-18.

## Fixture
Hand-built fixed-width **mode 0x02** PGEN, 1 variant × 4 samples, genocodes
`[S0=2, S1=1, S2=0, S3=3]` = (double-ALT, het, hom-REF, missing).
PVAR: `1 1000 rs0 REF=G ALT=A`.

Bytes (`xxd t0_lsb.pgen`): `6c1b 02  01000000  04000000  80  c6`
- `6c 1b 02` — magic + storage mode 0x02 (fixed-width, all type-0 records).
- `01000000` — `uint32` variant count = 1 (little-endian).
- `04000000` — `uint32` sample count = 4.
- `80` — header-format byte 11. Bits 6–7 = `0b10` = 2 = **all REF alleles
  provisional** (the honest choice for data with no tracked REF). Bits 0–5 = 0
  (no per-record type/length arrays for fixed-width; no allele-count bytes).
- `c6` — the single record byte (`ceil(4/4)=1` byte).

## Finding 1 — main-track 2-bit order is LSB-first  ✅
`0xc6 = 1100_0110`.
- LSB-first (sample 0 in low bits): `[2,1,0,3]`.
- MSB-first: `[3,0,1,2]`.

`plink2 --pfile t0_lsb --export A` produced, in sample order, `[0,1,2,NA]`
(under column `rs0_G`). Mapping export-value→genocode gives `[2,1,0,3]`, which
matches **LSB-first**. Confirmed.

## Finding 2 — mode-0x02 header layout is correct  ✅
plink2 loaded the 12-byte header + 1 record byte with no error and decoded all
four samples, confirming the layout above (magic/mode, LE u32 dims, byte-11
provisional-REF encoding). The writer (track A5) can emit exactly this.

## Finding 3 — orientation caveat for golden tests  ⚠️
`plink2 --export A` counted the **REF** allele here (column named `rs0_G`,
hom-REF sample exported as 2), which is the **opposite** of reigen's
`g = count(allele1)` with `allele1 = ALT`. So:
- reigen `g` for `[S0..S3]` = `[2,1,0,9]` (count of ALT=A).
- plink2 `--export A` (REF count) = `[0,1,2,NA]`.

Track C golden tests must NOT compare reigen `g` directly to `--export A`
output. Either force the counted allele with `--export-allele` (count ALT),
flip the comparison, or validate through a `.bed` round-trip. This is the
B-013 polarity trap — see [[reigen-pgen-polarity]] and `architecture.md` §1.3.

## Consequences for implementation
- Codec map (read): PGEN genocode `{0,1,2,3}` → reigen `{0,1,2,9}` is the
  **identity** when `allele1 = ALT` (PVAR ALT → SnpRow.allele1). Confirmed.
- Codec map (write): reigen `{0,1,2,9}` → PGEN `{0,1,2,3}`, packed LSB-first.
- Writer target: mode 0x02, header byte-11 = `0x80`.

---

# C1 — golden fixtures provenance (added 2026-06-18)

These fixtures are the **external plink2 oracle** for the PGEN reader/writer.
They were produced by **real `plink2` v2.0.0-a.6.13**, never by reigen.

> ⚠️ **External oracle — NEVER regenerate with reigen.** Per
> `tests/golden/README.md`: regenerating these from reigen's own output would
> make the test self-referential and blind to a genome-wide polarity flip
> (the B-013 class of bug). Regenerate only with real plink2 via
> `fixtures_genesis.sh`, and only when the fixture definitions intentionally
> change (Gary signs off).

Regenerate with `tests/golden/plink2/fixtures_genesis.sh` (needs real plink2 +
python3; deterministic, seed 1337 for the LD fixture). Committed source inputs:
`ld_src.vcf`, `multi_src.vcf` (both emitted by the genesis script).

## Polarity convention (re-confirmed for these fixtures)
`gold.pgen` round-trips byte-identically back to the committed golden `.bed`:
`plink2 --pfile gold --make-bed` reproduces `tests/golden/plink/gold.bed`
byte-for-byte, with `.bim` A1/A2 preserved. PVAR col-5 ALT = `.bim` A1 =
reigen `allele1`; the read map PGEN `{0,1,2,3}` → reigen `{0,1,2,9}` is the
identity. (Same trap as Finding 3: `--export A` counts REF, the opposite of
reigen `g = count(allele1)`; golden tests must round-trip through `.bed`, not
compare to `--export A`.)

## `gold.{pgen,pvar,psam}` — primary fixture
- **Dimensions:** 8 samples × 12 variants (same dataset as the convertf / `.bed`
  golden fixtures, so allele set + missing calls match).
- **Produced by:** `plink2 --bfile tests/golden/plink/gold --make-pgen`.
- **Storage mode:** `0x10` (variable-width) — confirmed via xxd (byte 2 = 0x10).
  Byte 11 = `0x80` (4-bit record types, 1-byte lengths, provisional-ref = 2 =
  all-provisional, since it came from a PLINK1 fileset).
- **Record-type distribution:** all 12 records **type 0 (uncompressed)** — the
  set is too small/dense to compress. Exercises the **variable-width header +
  type-0 record** path only. (This is expected and the reason `ld.*` exists.)

## `ld.{pgen,pvar,psam}` — compressible / LD fixture (R2) — KEY RESULT
- **Dimensions:** 400 samples × 400 variants.
- **Produced by:** generated `ld_src.vcf` (deterministic, seed 1337; LD blocks of
  near-duplicate variants, 0↔2-inverted dups, sparse-ALT, dense-ALT, and
  mostly-missing variants) → `plink2 --vcf ld_src.vcf --make-pgen`.
- **Storage mode:** `0x10` (variable-width). Byte 11 = `0x40` (4-bit types,
  1-byte lengths, provisional-ref = 1 = none-provisional, imported from VCF).
- **Achieved record-type distribution** (vrtype & 0x07), out of 400 records:

  | type | meaning                | count |
  |------|------------------------|-------|
  | 0    | uncompressed           | 49    |
  | 1    | 1-bit                  | 51    |
  | 2    | LD-compressed          | 198   |
  | 3    | LD-compressed inverted | 25    |
  | 4    | difflist (homREF base) | 31    |
  | 6    | difflist (dblALT base) | 26    |
  | 7    | difflist (missing base)| 20    |

  **All 7 hardcall record types (0,1,2,3,4,6,7) are present.** Type 5 is
  reserved/unused by the spec, so it cannot appear. This fixture therefore
  covers every reader decoder path: the uncompressed reader, the 1-bit reader,
  both LD readers (plain + inverted, requiring a maintained last-non-LD anchor),
  and all three difflist-from-common readers — plus the difflist + base-128
  varint machinery they share. No multiallelic/phase/dosage aux bits are set
  (GT-only VCF, all biallelic).

## `multi.{pgen,pvar,psam}` — multiallelic fixture
- **Dimensions:** 6 samples × 4 variants; **rs1 is triallelic** (PVAR
  `ALT = A,T`).
- **Produced by:** generated `multi_src.vcf` → `plink2 --vcf … --make-pgen`.
- **Storage mode:** `0x10`. Byte 11 = `0x40`. `allele_ct_bytes = 0` — plink2
  does **not** store allele counts in this PGEN header; multiallelic status is
  signalled by the comma in PVAR ALT and by **vrtype bit 3 set** on rs1
  (vrtype `0x08`). This matches plan §4.4: the reader derives the drop-mask from
  PVAR comma-ALT, and may also defensively treat bit-3-set records as
  multiallelic.
- **Record types:** all 4 records category type-0; rs0/rs2/rs3 plain biallelic,
  rs1 has bit 3 set. Exercises the **"drop multiallelic, report, realign"** path.
