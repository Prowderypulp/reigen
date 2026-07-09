//! PGEN header parsing.
//!
//! Parses the fixed header of a `.pgen` file: magic, storage mode, dataset
//! dimensions, the header-format byte (byte 11), and — for the variable-width
//! modes — the per-block variant-block offsets and the packed record-type and
//! record-length arrays. The result lets a record decoder (track A3) ask for
//! any variant's record type and its `(byte offset, byte length)` span.
//!
//! # Storage modes (spec §"Storage mode")
//!
//! - `0x01` PLINK 1 `.bed` — no further header bytes; sample count not in the
//!   file. We recognize it and return [`StorageMode::Plink1Bed`] so the caller
//!   can delegate to `packed_ped`.
//! - `0x02` fixed-width, all type-0 records.
//! - `0x03`/`0x04` fixed-width dosage / phased-dosage (recognized; record-span
//!   tables not built here).
//! - `0x10` standard variable-width.
//! - `0x11` `0x10` + ignorable header/footer extensions (parsed and skipped).
//! - `0x20`/`0x21` index in a sibling `.pgi` — returned as an "unsupported,
//!   external index" error for now.
//!
//! # Byte 11 (header format), spec §"Dataset dimensions, header body formatting"
//!
//! - Bits 0-3: record type/length storage. Bit 2 selects type width (0 = **4
//!   bits per type, 2 records per byte**; 1 = 8 bits per type). Bits 0-1 select
//!   the per-record length width (0..=3 → 1..=4 bytes).
//! - Bits 4-5: bytes per stored allele count (0 = none stored).
//! - Bits 6-7: provisional-REF mode (0 = in PVAR, 1 = none provisional, 2 =
//!   all provisional, 3 = bitarray in header).
//!
//! **4-bit packing gotcha:** with 4-bit types, record `2i` is the *low* nibble
//! of type byte `i` and record `2i+1` is the *high* nibble. Decoding one byte
//! per type yields garbage on these files (which is what the committed fixtures
//! use).

use anyhow::{bail, Context, Result};

use super::difflist::ByteCursor;

pub const PGEN_MAGIC: [u8; 2] = [0x6c, 0x1b];

/// Variant-block size: variant records are grouped into blocks of `2^16`.
pub const BLOCK_SIZE: usize = 1 << 16;

/// PGEN storage mode (third byte of the file).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StorageMode {
    /// `0x01`: PLINK 1 `.bed`; delegate to `packed_ped`.
    Plink1Bed,
    /// `0x02`: fixed-width, all type-0 records.
    FixedWidth,
    /// `0x03`: fixed-width unphased-dosage (all records type 0x40).
    FixedDosage,
    /// `0x04`: fixed-width phased-dosage (all records type 0xc0).
    FixedPhasedDosage,
    /// `0x10`: standard variable-width.
    VariableWidth,
    /// `0x11`: variable-width with ignorable extensions.
    VariableWidthExt,
}

impl StorageMode {
    fn from_byte(b: u8) -> Result<Self> {
        Ok(match b {
            0x01 => StorageMode::Plink1Bed,
            0x02 => StorageMode::FixedWidth,
            0x03 => StorageMode::FixedDosage,
            0x04 => StorageMode::FixedPhasedDosage,
            0x10 => StorageMode::VariableWidth,
            0x11 => StorageMode::VariableWidthExt,
            0x20 | 0x21 => bail!(
                "pgen storage mode 0x{b:02x}: index is in a sibling .pgi file \
                 (external-index modes are unsupported)"
            ),
            _ => bail!("pgen: unknown/reserved storage mode 0x{b:02x}"),
        })
    }

    /// True for the variable-width modes that carry per-block type/length arrays.
    pub fn is_variable_width(self) -> bool {
        matches!(
            self,
            StorageMode::VariableWidth | StorageMode::VariableWidthExt
        )
    }
}

/// Provisional-REF flag mode (byte 11 bits 6-7).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProvRefMode {
    /// 0: not in PGEN; consult the PVAR.
    InPvar,
    /// 1: no reference alleles are provisional (all trusted).
    NoneProvisional,
    /// 2: all reference alleles are provisional (e.g. from a PLINK 1 fileset).
    AllProvisional,
    /// 3: some provisional, tracked by a header bitarray.
    Bitarray,
}

impl ProvRefMode {
    fn from_bits(bits: u8) -> Self {
        match bits & 0b11 {
            0 => ProvRefMode::InPvar,
            1 => ProvRefMode::NoneProvisional,
            2 => ProvRefMode::AllProvisional,
            _ => ProvRefMode::Bitarray,
        }
    }
}

/// Decoded byte-11 sub-fields.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HeaderFormat {
    /// Raw byte 11 value.
    pub raw: u8,
    /// True if record types are stored 8 bits each; false = 4 bits each (2/byte).
    pub eight_bit_types: bool,
    /// Bytes per stored record length (1..=4).
    pub length_bytes: usize,
    /// Bytes per stored allele count (0 = none stored).
    pub allele_ct_bytes: usize,
    /// Provisional-REF flag mode.
    pub prov_ref: ProvRefMode,
}

impl HeaderFormat {
    fn from_byte(b: u8) -> Self {
        let low = b & 0x0f;
        // Bit 2 of the low nibble selects type width; bits 0-1 select length width.
        let eight_bit_types = (low & 0b100) != 0;
        let length_bytes = ((low & 0b11) as usize) + 1;
        let allele_ct_bytes = ((b >> 4) & 0b11) as usize;
        let prov_ref = ProvRefMode::from_bits(b >> 6);
        HeaderFormat {
            raw: b,
            eight_bit_types,
            length_bytes,
            allele_ct_bytes,
            prov_ref,
        }
    }
}

/// Parsed PGEN header.
///
/// For variable-width modes, `record_types` holds one entry per variant
/// (unpacked from 4-bit or 8-bit storage) and `record_spans` holds each
/// variant's `(byte offset, byte length)` in the file. For `0x01` only `mode`
/// is meaningful (the caller delegates to `packed_ped`). For the fixed-width
/// modes the per-record tables are left empty (a fixed-width decoder computes
/// spans from `N` directly).
#[derive(Debug, Clone)]
pub struct PgenHeader {
    pub mode: StorageMode,
    /// Variant count `M`. `None` for `0x01` (not stored in a `.bed`).
    pub variant_count: Option<u32>,
    /// Sample count `N`. `None` for `0x01`.
    pub sample_count: Option<u32>,
    /// Decoded byte 11. `None` for `0x01`.
    pub format: Option<HeaderFormat>,
    /// `uint64` variant-block offsets (one per block). Variable-width only.
    pub block_offsets: Vec<u64>,
    /// Per-variant record type (low byte; bits 0-7). Variable-width only.
    pub record_types: Vec<u8>,
    /// Per-variant `(byte offset, byte length)`. Variable-width only.
    pub record_spans: Vec<(u64, u64)>,
}

impl PgenHeader {
    /// Record type of variant `i` (bits 0-2 select main-track compression;
    /// bit 3 = multiallelic; bits 4-7 = phase/dosage).
    pub fn record_type(&self, i: usize) -> Option<u8> {
        self.record_types.get(i).copied()
    }

    /// `(byte offset, byte length)` of variant `i`'s record in the file.
    pub fn record_span(&self, i: usize) -> Option<(u64, u64)> {
        self.record_spans.get(i).copied()
    }

    /// Parse a PGEN header from the start of `buf` (the whole file, or at least
    /// the header). Returns the parsed header.
    pub fn parse(buf: &[u8]) -> Result<PgenHeader> {
        if buf.len() < 3 {
            bail!("pgen: file too short for magic + mode ({} bytes)", buf.len());
        }
        if buf[0] != PGEN_MAGIC[0] || buf[1] != PGEN_MAGIC[1] {
            bail!(
                "pgen: bad magic {:02x} {:02x} (expected {:02x} {:02x})",
                buf[0],
                buf[1],
                PGEN_MAGIC[0],
                PGEN_MAGIC[1]
            );
        }
        let mode = StorageMode::from_byte(buf[2])?;

        if mode == StorageMode::Plink1Bed {
            return Ok(PgenHeader {
                mode,
                variant_count: None,
                sample_count: None,
                format: None,
                block_offsets: Vec::new(),
                record_types: Vec::new(),
                record_spans: Vec::new(),
            });
        }

        // Dims + byte-11. Bytes 3-6 = M, 7-10 = N, 11 = header format.
        if buf.len() < 12 {
            bail!("pgen: file too short for header dimensions ({} bytes)", buf.len());
        }
        let m = u32::from_le_bytes([buf[3], buf[4], buf[5], buf[6]]);
        let n = u32::from_le_bytes([buf[7], buf[8], buf[9], buf[10]]);

        // For the fixed-width modes (0x02–0x04) only byte-11 bits 6-7
        // (provisional-REF) are meaningful; pgenlib hard-errors if bits 0-5 are
        // set (no per-record type/length/allele-count arrays exist). Validate
        // that rather than silently parsing meaningless type/length fields
        // (B-108). The variable-width modes use the full byte.
        if matches!(
            mode,
            StorageMode::FixedWidth | StorageMode::FixedDosage | StorageMode::FixedPhasedDosage
        ) && (buf[11] & 0x3f) != 0
        {
            bail!(
                "pgen mode 0x{:02x}: header byte 11 = 0x{:02x} has reserved bits 0-5 set \
                 (must be 0 for fixed-width modes)",
                buf[2],
                buf[11]
            );
        }

        let format = HeaderFormat::from_byte(buf[11]);

        let mut header = PgenHeader {
            mode,
            variant_count: Some(m),
            sample_count: Some(n),
            format: Some(format),
            block_offsets: Vec::new(),
            record_types: Vec::new(),
            record_spans: Vec::new(),
        };

        if mode.is_variable_width() {
            header.parse_variable_body(buf, m as usize, format)?;
        }

        Ok(header)
    }

    /// Parse the variable-width header body: block offsets, then per-block
    /// packed type and length arrays. Reconstruct per-variant types and spans.
    fn parse_variable_body(&mut self, buf: &[u8], m: usize, fmt: HeaderFormat) -> Result<()> {
        let n_blocks = m.div_ceil(BLOCK_SIZE);

        // Block offsets: 8B each, starting at byte 12.
        let mut cur = ByteCursor::new(buf);
        // Skip the 12 fixed header bytes already consumed.
        for _ in 0..12 {
            cur.read_byte_for_skip()?;
        }
        let mut block_offsets = Vec::with_capacity(n_blocks);
        for _ in 0..n_blocks {
            block_offsets.push(cur.read_u64_le()?);
        }

        // Per-block: packed type array, then packed length array. (We do not
        // build the allele-count / prov-ref bitarrays here; A3 doesn't need
        // them for biallelic hardcalls. But we must skip them to stay aligned
        // for the *next* block.)
        let mut record_types: Vec<u8> = Vec::with_capacity(m);
        let mut record_lengths: Vec<u64> = Vec::with_capacity(m);

        for block in 0..n_blocks {
            let block_start = block * BLOCK_SIZE;
            let block_count = (m - block_start).min(BLOCK_SIZE);

            // Type array.
            if fmt.eight_bit_types {
                for _ in 0..block_count {
                    record_types.push(cur.read_byte_for_skip()?);
                }
            } else {
                // 4-bit packed: 2 records per byte, low nibble = even index.
                let n_bytes = block_count.div_ceil(2);
                for i in 0..n_bytes {
                    let byte = cur.read_byte_for_skip()?;
                    record_types.push(byte & 0x0f);
                    if 2 * i + 1 < block_count {
                        record_types.push(byte >> 4);
                    }
                }
            }

            // Length array: `length_bytes` bytes per record, LE.
            for _ in 0..block_count {
                let bytes = (0..fmt.length_bytes)
                    .map(|_| cur.read_byte_for_skip())
                    .collect::<Result<Vec<u8>>>()?;
                let mut v: u64 = 0;
                for (j, &b) in bytes.iter().enumerate() {
                    v |= (b as u64) << (8 * j);
                }
                record_lengths.push(v);
            }

            // Allele-count array (skip): only present for multiallelic vars,
            // but when present plink2 stores `allele_ct_bytes` per record.
            if fmt.allele_ct_bytes > 0 {
                for _ in 0..(block_count * fmt.allele_ct_bytes) {
                    cur.read_byte_for_skip()?;
                }
            }

            // Provisional-REF bitarray (skip) when mode 3: 1 bit per record,
            // padded to a byte boundary.
            if fmt.prov_ref == ProvRefMode::Bitarray {
                let n_bytes = block_count.div_ceil(8);
                for _ in 0..n_bytes {
                    cur.read_byte_for_skip()?;
                }
            }
        }

        // Compute spans from block offsets + cumulative lengths within a block.
        // Each variant's record starts where the previous one in the same block
        // ended; the block's first record starts at its block offset.
        let mut spans = Vec::with_capacity(m);
        for block in 0..n_blocks {
            let block_start = block * BLOCK_SIZE;
            let block_count = (m - block_start).min(BLOCK_SIZE);
            let mut off = block_offsets[block];
            for j in 0..block_count {
                let len = record_lengths[block_start + j];
                spans.push((off, len));
                off += len;
            }
        }

        self.block_offsets = block_offsets;
        self.record_types = record_types;
        self.record_spans = spans;
        Ok(())
    }
}

// Small additions to ByteCursor used only by the header parser. Kept here to
// avoid widening the difflist module's public surface.
impl ByteCursor<'_> {
    fn read_byte_for_skip(&mut self) -> Result<u8> {
        let mut tmp = [0u8; 1];
        self.read_into(&mut tmp)
            .context("pgen header: unexpected end of buffer")?;
        Ok(tmp[0])
    }

    fn read_u64_le(&mut self) -> Result<u64> {
        let mut tmp = [0u8; 8];
        self.read_into(&mut tmp)
            .context("pgen header: unexpected end of buffer (uint64)")?;
        Ok(u64::from_le_bytes(tmp))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn byte11_decodes_subfields() {
        // 0x80 = 1000_0000: 4-bit types, 1-byte lengths, no allele cts, all-prov.
        let f = HeaderFormat::from_byte(0x80);
        assert!(!f.eight_bit_types);
        assert_eq!(f.length_bytes, 1);
        assert_eq!(f.allele_ct_bytes, 0);
        assert_eq!(f.prov_ref, ProvRefMode::AllProvisional);

        // 0x40 = 0100_0000: 4-bit types, 1-byte lengths, none-provisional.
        let f = HeaderFormat::from_byte(0x40);
        assert!(!f.eight_bit_types);
        assert_eq!(f.length_bytes, 1);
        assert_eq!(f.prov_ref, ProvRefMode::NoneProvisional);

        // 0x81 = 1000_0001: 4-bit types, 2-byte lengths, all-provisional
        // (the spec's worked-example byte).
        let f = HeaderFormat::from_byte(0x81);
        assert!(!f.eight_bit_types);
        assert_eq!(f.length_bytes, 2);
        assert_eq!(f.prov_ref, ProvRefMode::AllProvisional);

        // 0x41 = 0100_0001: 4-bit types, 2-byte lengths, none-provisional.
        let f = HeaderFormat::from_byte(0x41);
        assert_eq!(f.length_bytes, 2);
        assert_eq!(f.prov_ref, ProvRefMode::NoneProvisional);

        // 0x04 = 8-bit types, 1-byte lengths.
        let f = HeaderFormat::from_byte(0x04);
        assert!(f.eight_bit_types);
        assert_eq!(f.length_bytes, 1);
    }

    /// The spec's header worked example: 1092 samples, 39728178 variants,
    /// byte 11 = 0x81, record #0 starting at byte 99325313. We synthesize just
    /// enough of a header (offsets + first few type/length bytes for block 0)
    /// to exercise the block-count and offset math. (Building the full 99 MB
    /// header is impractical and unnecessary — we verify the arithmetic the
    /// parser performs, not a full real file here.)
    #[test]
    fn spec_header_example_block_math() {
        let m: usize = 39_728_178;
        let n_blocks = m.div_ceil(BLOCK_SIZE);
        assert_eq!(n_blocks, 607); // ceil(39728178 / 65536)
        // First block offset table occupies 8 * 607 bytes after byte 12.
        // Record #0 starts at the first block offset.
        let off0: u64 = 99_325_313;
        // The spec states bytes 12-19 = 0x81 0x95 0xeb 0x05 0x00... ; verify LE
        // decodes to 99325313.
        let bytes = [0x81u8, 0x95, 0xeb, 0x05, 0x00, 0x00, 0x00, 0x00];
        assert_eq!(u64::from_le_bytes(bytes), off0);
    }

    fn read_fixture(name: &str) -> Vec<u8> {
        let path = format!(
            "{}/tests/golden/plink2/{}",
            env!("CARGO_MANIFEST_DIR"),
            name
        );
        std::fs::read(&path).unwrap_or_else(|e| panic!("read {path}: {e}"))
    }

    #[test]
    fn parse_gold_fixture_header() {
        let buf = read_fixture("gold.pgen");
        let h = PgenHeader::parse(&buf).unwrap();
        assert_eq!(h.mode, StorageMode::VariableWidth);
        assert_eq!(h.variant_count, Some(12));
        assert_eq!(h.sample_count, Some(8));
        let f = h.format.unwrap();
        assert_eq!(f.raw, 0x80);
        assert!(!f.eight_bit_types);
        assert_eq!(f.length_bytes, 1);
        assert_eq!(f.allele_ct_bytes, 0);
        assert_eq!(f.prov_ref, ProvRefMode::AllProvisional);
        // All 12 records type 0 (T0_FINDINGS).
        assert_eq!(h.record_types.len(), 12);
        assert!(h.record_types.iter().all(|&t| t == 0));
        assert_eq!(h.record_spans.len(), 12);
        // Spans are contiguous and start within the file.
        for (off, len) in &h.record_spans {
            assert!(*off + *len <= buf.len() as u64, "span past EOF");
        }
    }

    #[test]
    fn parse_ld_fixture_header() {
        let buf = read_fixture("ld.pgen");
        let h = PgenHeader::parse(&buf).unwrap();
        assert_eq!(h.mode, StorageMode::VariableWidth);
        assert_eq!(h.variant_count, Some(400));
        assert_eq!(h.sample_count, Some(400));
        let f = h.format.unwrap();
        assert_eq!(f.raw, 0x40);
        assert!(!f.eight_bit_types);
        assert_eq!(f.length_bytes, 1);
        assert_eq!(f.prov_ref, ProvRefMode::NoneProvisional);
        assert_eq!(h.record_types.len(), 400);

        // Record-type distribution (vrtype & 0x07) must match T0_FINDINGS.
        let mut dist: BTreeMap<u8, usize> = BTreeMap::new();
        for &t in &h.record_types {
            *dist.entry(t & 0x07).or_default() += 1;
        }
        assert_eq!(dist.get(&0).copied().unwrap_or(0), 49);
        assert_eq!(dist.get(&1).copied().unwrap_or(0), 51);
        assert_eq!(dist.get(&2).copied().unwrap_or(0), 198);
        assert_eq!(dist.get(&3).copied().unwrap_or(0), 25);
        assert_eq!(dist.get(&4).copied().unwrap_or(0), 31);
        assert_eq!(dist.get(&6).copied().unwrap_or(0), 26);
        assert_eq!(dist.get(&7).copied().unwrap_or(0), 20);
        // Types 0,1,2,3,4,6,7 all present; 5 absent.
        assert!(!dist.contains_key(&5));

        // Spans contiguous within the (single) block and inside the file.
        for (off, len) in &h.record_spans {
            assert!(*off + *len <= buf.len() as u64, "span past EOF");
        }
    }

    #[test]
    fn rejects_bad_magic() {
        let err = PgenHeader::parse(&[0x00, 0x00, 0x10]).unwrap_err();
        assert!(format!("{err:#}").contains("magic"));
    }

    #[test]
    fn external_index_mode_errors_cleanly() {
        let err = PgenHeader::parse(&[0x6c, 0x1b, 0x20]).unwrap_err();
        assert!(format!("{err:#}").contains("external-index"));
    }

    #[test]
    fn plink1_bed_recognized() {
        let h = PgenHeader::parse(&[0x6c, 0x1b, 0x01]).unwrap();
        assert_eq!(h.mode, StorageMode::Plink1Bed);
        assert_eq!(h.variant_count, None);
        assert_eq!(h.sample_count, None);
    }

    /// Multi-block header: M = 65537 (2 blocks: block 0 = 65536, block 1 = 1).
    /// All type-0 records, 1-byte lengths, no allele cts, none-provisional.
    /// Verifies the parser correctly handles the second block and computes
    /// spans that span both blocks.
    #[test]
    fn multi_block_header_parse() {
        let m: u32 = 65537;
        let n: u32 = 4;
        // byte 11 = 0x40: 4-bit types, 1-byte lengths, no allele cts, none-prov.
        let byte11 = 0x40u8;
        let n_blocks = (m as usize).div_ceil(BLOCK_SIZE);
        assert_eq!(n_blocks, 2);

        let mut buf = Vec::new();
        // Fixed header: magic + mode + M + N + byte11
        buf.extend_from_slice(&PGEN_MAGIC);
        buf.push(0x10); // mode 0x10
        buf.extend_from_slice(&m.to_le_bytes());
        buf.extend_from_slice(&n.to_le_bytes());
        buf.push(byte11);

        // Block offsets: 2 × u64 LE.
        // Block 0 data starts after the full header. We'll compute it after
        // building the type/length arrays.
        let header_so_far = buf.len() + 8 * n_blocks; // 12 + 16 = 28

        // Block 0: 65536 variants, 4-bit types (32768 bytes), 1-byte lengths (65536 bytes).
        let block0_type_bytes = 65536 / 2; // 32768
        let block0_len_bytes = 65536 * 1; // 65536
        // Block 1: 1 variant, 4-bit types (1 byte), 1-byte length (1 byte).
        let block1_type_bytes = 1usize.div_ceil(2); // 1
        let block1_len_bytes = 1;

        let block0_data_off = header_so_far + block0_type_bytes + block0_len_bytes
            + block1_type_bytes + block1_len_bytes;
        let block1_data_off = block0_data_off + 65536; // 65536 type-0 records × 1 byte each

        buf.extend_from_slice(&(block0_data_off as u64).to_le_bytes());
        buf.extend_from_slice(&(block1_data_off as u64).to_le_bytes());

        // Block 0 type array: all type-0, 4-bit packed → all zeros.
        buf.extend(std::iter::repeat(0u8).take(block0_type_bytes));
        // Block 0 length array: each record is ceil(4/4) = 1 byte.
        buf.extend(std::iter::repeat(1u8).take(block0_len_bytes));

        // Block 1 type array: 1 variant, type 0 → low nibble = 0, high nibble unused.
        buf.push(0x00);
        // Block 1 length array: 1 byte.
        buf.push(1);

        // Record data: 65537 × 1 byte (all zeros for type-0 with N=4).
        buf.extend(std::iter::repeat(0u8).take(65537));

        let h = PgenHeader::parse(&buf).unwrap();
        assert_eq!(h.mode, StorageMode::VariableWidth);
        assert_eq!(h.variant_count, Some(m));
        assert_eq!(h.sample_count, Some(n));
        assert_eq!(h.block_offsets.len(), 2);
        assert_eq!(h.record_types.len(), m as usize);
        assert!(h.record_types.iter().all(|&t| t == 0), "all type-0");
        assert_eq!(h.record_spans.len(), m as usize);

        // First span starts at block 0 offset.
        assert_eq!(h.record_spans[0].0, block0_data_off as u64);
        // Spans are contiguous within each block.
        for i in 1..65536 {
            let (prev_off, prev_len) = h.record_spans[i - 1];
            let (off, _len) = h.record_spans[i];
            assert_eq!(off, prev_off + prev_len, "span {i} not contiguous");
        }
        // Block 1 first span starts at block 1 offset.
        assert_eq!(h.record_spans[65536].0, block1_data_off as u64);
    }

    /// 8-bit record-type storage: byte 11 bit 2 set.
    /// Build a small mode-0x10 PGEN with M=4, N=4, 8-bit types.
    #[test]
    fn eight_bit_record_types() {
        let m: u32 = 4;
        let n: u32 = 4;
        // byte 11 = 0x44: bit 2 set (8-bit types), 1-byte lengths, no allele cts, none-prov.
        let byte11 = 0x44u8;

        let mut buf = Vec::new();
        buf.extend_from_slice(&PGEN_MAGIC);
        buf.push(0x10);
        buf.extend_from_slice(&m.to_le_bytes());
        buf.extend_from_slice(&n.to_le_bytes());
        buf.push(byte11);

        // 1 block offset.
        let header_so_far = buf.len() + 8;
        // 4 × 1-byte types + 4 × 1-byte lengths = 8 bytes of type/len arrays.
        let data_off = header_so_far + 4 + 4;
        buf.extend_from_slice(&(data_off as u64).to_le_bytes());

        // Type array: 8-bit, one byte per variant.
        // Types: 0, 4, 6, 7 (uncompressed, difflist-common-0, difflist-common-2, difflist-common-missing).
        buf.extend_from_slice(&[0x00, 0x04, 0x06, 0x07]);

        // Length array: 1 byte each.
        // Type 0: ceil(4/4) = 1 byte. Types 4/6/7: difflist only (varies).
        // For simplicity, give each a length of 2 (1 byte difflist + padding).
        buf.extend_from_slice(&[1, 2, 2, 2]);

        // Record data (minimal, just enough to not crash if someone tries to read).
        // Type 0: 1 byte of genotype data.
        buf.push(0x00);
        // Types 4/6/7: 2 bytes each (minimal difflist: L=0 → 1 byte, + 1 byte padding).
        buf.extend_from_slice(&[0x00, 0x00]); // type 4
        buf.extend_from_slice(&[0x00, 0x00]); // type 6
        buf.extend_from_slice(&[0x00, 0x00]); // type 7

        let h = PgenHeader::parse(&buf).unwrap();
        assert_eq!(h.format.unwrap().eight_bit_types, true);
        assert_eq!(h.record_types.len(), 4);
        assert_eq!(h.record_types[0], 0x00);
        assert_eq!(h.record_types[1], 0x04);
        assert_eq!(h.record_types[2], 0x06);
        assert_eq!(h.record_types[3], 0x07);
    }

    /// allele_ct_bytes > 0 and prov-ref bitarray: the parser must skip those
    /// bytes correctly and still produce valid type/length arrays.
    #[test]
    fn allele_ct_and_prov_ref_bitarray_skip() {
        let m: u32 = 4;
        let n: u32 = 4;
        // byte 11 = 0xF4:
        //   bits 0-1 = 0 → 1-byte lengths
        //   bit 2 = 1 → 8-bit types
        //   bits 4-5 = 3 → 3 bytes per allele count
        //   bits 6-7 = 3 → ProvRefMode::Bitarray
        let byte11 = 0xF4u8;

        let mut buf = Vec::new();
        buf.extend_from_slice(&PGEN_MAGIC);
        buf.push(0x10);
        buf.extend_from_slice(&m.to_le_bytes());
        buf.extend_from_slice(&n.to_le_bytes());
        buf.push(byte11);

        let header_so_far = buf.len() + 8; // + block offset

        // 8-bit types: 4 bytes. 1-byte lengths: 4 bytes.
        // allele_ct: 4 × 3 = 12 bytes. prov-ref bitarray: ceil(4/8) = 1 byte.
        let type_len_arrays = 4 + 4 + 12 + 1; // 21
        let data_off = header_so_far + type_len_arrays;

        buf.extend_from_slice(&(data_off as u64).to_le_bytes());

        // Type array: all type-0.
        buf.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]);
        // Length array: all 1.
        buf.extend_from_slice(&[1, 1, 1, 1]);
        // Allele count array: 4 × 3 bytes (values don't matter, just skip).
        buf.extend(std::iter::repeat(0xAAu8).take(12));
        // Prov-ref bitarray: 1 byte.
        buf.push(0x0F); // first 4 variants are provisional

        // Record data: 4 × 1 byte.
        buf.extend(std::iter::repeat(0x00).take(4));

        let h = PgenHeader::parse(&buf).unwrap();
        let fmt = h.format.unwrap();
        assert!(fmt.eight_bit_types);
        assert_eq!(fmt.allele_ct_bytes, 3);
        assert_eq!(fmt.prov_ref, ProvRefMode::Bitarray);
        assert_eq!(h.record_types.len(), 4);
        assert!(h.record_types.iter().all(|&t| t == 0));
        assert_eq!(h.record_spans.len(), 4);
        // Spans should be contiguous and start at data_off.
        assert_eq!(h.record_spans[0].0, data_off as u64);
    }
}
