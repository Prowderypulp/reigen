//! PLINK `.bed` format (PACKEDPED).
//!
//! # On-disk layout
//!
//! - **Magic** (3 bytes): `0x6c 0x1b 0x01`. Third byte = mode; `0x01` =
//!   SNP-major (one record per SNP, samples packed across). convertf
//!   only emits SNP-major; we accept only SNP-major.
//! - **Records**: `nsnp` × `rec_bytes` where `rec_bytes = ceil(nind / 4)`.
//!   2 bits per genotype, **LSB-first within byte** (sample 0 at bits
//!   1-0, sample 1 at bits 3-2, etc.).
//!
//! # Encoding (PLINK convention — differs from PACKEDANCESTRYMAP)
//!
//! | 2 bits | PLINK meaning | AdmixTools 2-bit |
//! |--------|---------------|------------------|
//! | `00`   | hom A1        | `10` (g=2)       |
//! | `01`   | missing       | `11` (missing)   |
//! | `10`   | het           | `01` (g=1)       |
//! | `11`   | hom A2        | `00` (g=0)       |
//!
//! Recall that in `.bim` we already swap A1/A2 so that AdmixTools allele1
//! corresponds to PLINK A2 (see `meta/bim.rs` doc). So:
//! - PLINK "hom A2" = two copies of `.bim` col 6 = two copies of
//!   AdmixTools allele1 = genotype 0 in AdmixTools encoding.
//! - PLINK "hom A1" = two copies of AdmixTools allele2 = genotype 2.
//!
//! Combined with the bit-order reverse, each 8-bit PLINK byte maps to an
//! 8-bit AdmixTools byte via a 256-entry LUT. Same for the inverse.
//!
//! # Padding
//!
//! Last record byte is padded with `01` bits (missing) when `nind % 4 != 0`
//! per PLINK spec. We follow that on write; we mask on read (the canonical
//! AdmixTools buffer is `ceil(nind*2/8)` bytes — last byte's pad bits
//! could already be zero by construction, no need to re-mask).

use super::{GenoReader, GenoWriter, Layout};
use anyhow::{bail, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

pub const BED_MAGIC: [u8; 3] = [0x6c, 0x1b, 0x01];

// ======================================================================
// LUTs: full-byte PLINK ↔ AdmixTools recode + bit-reverse
// ======================================================================

/// For each of 256 input bytes, the PLINK-encoded byte is converted to
/// an AdmixTools-encoded byte. Encoding change: PLINK {00,01,10,11} →
/// AM {10,11,01,00}. Plus bit-order reverse within byte.
static PLINK_TO_AM: [u8; 256] = build_plink_to_am();

/// Inverse: AM-encoded byte → PLINK-encoded byte.
static AM_TO_PLINK: [u8; 256] = build_am_to_plink();

const fn recode_plink_to_am(two: u8) -> u8 {
    // 00→10, 01→11, 10→01, 11→00
    match two & 0b11 {
        0b00 => 0b10,
        0b01 => 0b11,
        0b10 => 0b01,
        _ => 0b00,
    }
}

const fn recode_am_to_plink(two: u8) -> u8 {
    match two & 0b11 {
        0b00 => 0b11,
        0b01 => 0b10,
        0b10 => 0b00,
        _ => 0b01,
    }
}

const fn build_plink_to_am() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        // PLINK LSB-first: sample 0 at bits 1-0, s1 at 3-2, s2 at 5-4, s3 at 7-6.
        let s0 = b & 0b11;
        let s1 = (b >> 2) & 0b11;
        let s2 = (b >> 4) & 0b11;
        let s3 = (b >> 6) & 0b11;
        // AM MSB-first: s0 at bits 7-6, s1 at 5-4, s2 at 3-2, s3 at 1-0.
        let out = (recode_plink_to_am(s0) << 6)
            | (recode_plink_to_am(s1) << 4)
            | (recode_plink_to_am(s2) << 2)
            | recode_plink_to_am(s3);
        t[i] = out;
        i += 1;
    }
    t
}

const fn build_am_to_plink() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        // AM MSB-first: s0 at bits 7-6, s1 at 5-4, s2 at 3-2, s3 at 1-0.
        let s0 = (b >> 6) & 0b11;
        let s1 = (b >> 4) & 0b11;
        let s2 = (b >> 2) & 0b11;
        let s3 = b & 0b11;
        // PLINK LSB-first: s0 at 1-0, s1 at 3-2, s2 at 5-4, s3 at 7-6.
        let out = recode_am_to_plink(s0)
            | (recode_am_to_plink(s1) << 2)
            | (recode_am_to_plink(s2) << 4)
            | (recode_am_to_plink(s3) << 6);
        t[i] = out;
        i += 1;
    }
    t
}

// ======================================================================
// Reader
// ======================================================================

pub struct PackedPedReader {
    mmap: Mmap,
    nind: usize,
    nsnp: usize,
    /// Bytes per record in PLINK layout: ceil(nind / 4).
    plink_rec_bytes: usize,
    /// Bytes per record in AdmixTools canonical layout: ceil(nind * 2 / 8).
    /// Same number; kept separate for clarity.
    am_rec_bytes: usize,
    /// Scratch for the last (possibly partial) byte during recode.
    last_byte_valid_samples: usize,
    next_idx: usize,
}

impl PackedPedReader {
    /// Jump the streaming cursor to SNP index `idx`. Used by merge.
    pub fn set_next_idx(&mut self, idx: usize) {
        self.next_idx = idx;
    }

    pub fn open(path: &Path, nind: usize, nsnp: usize) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;

        // Verify magic without mmap'ing an empty file.
        let mut magic = [0u8; 3];
        file.read_exact(&mut magic)
            .with_context(|| format!("read magic from {}", path.display()))?;
        if magic[0] != BED_MAGIC[0] || magic[1] != BED_MAGIC[1] {
            bail!(
                "{}: not a PLINK .bed file (magic {:02x} {:02x} {:02x} != {:02x} {:02x} {:02x})",
                path.display(),
                magic[0],
                magic[1],
                magic[2],
                BED_MAGIC[0],
                BED_MAGIC[1],
                BED_MAGIC[2]
            );
        }
        if magic[2] != 0x01 {
            bail!(
                "{}: sample-major .bed not supported (only SNP-major, mode byte {:02x})",
                path.display(),
                magic[2]
            );
        }

        // Now mmap the whole file (including magic — offsets start at 3).
        let file = File::open(path)?;
        let mmap =
            unsafe { Mmap::map(&file) }.with_context(|| format!("mmap {}", path.display()))?;

        let plink_rec_bytes = (nind + 3) / 4;
        let expected_len = 3 + plink_rec_bytes * nsnp;
        if mmap.len() < expected_len {
            bail!(
                "{}: file size {} < expected {} (nind={}, nsnp={}, rec_bytes={})",
                path.display(),
                mmap.len(),
                expected_len,
                nind,
                nsnp,
                plink_rec_bytes
            );
        }
        if mmap.len() > expected_len {
            log::warn!(
                "{}: {} trailing bytes past expected end",
                path.display(),
                mmap.len() - expected_len
            );
        }

        let am_rec_bytes = (nind * 2 + 7) / 8;
        // PLINK rec_bytes = ceil(nind/4) = ceil(nind*2/8) = am_rec_bytes.
        debug_assert_eq!(plink_rec_bytes, am_rec_bytes);

        let last_byte_valid_samples = ((nind - 1) % 4) + 1;

        Ok(Self {
            mmap,
            nind,
            nsnp,
            plink_rec_bytes,
            am_rec_bytes,
            last_byte_valid_samples,
            next_idx: 0,
        })
    }
}

impl GenoReader for PackedPedReader {
    fn nind(&self) -> usize {
        self.nind
    }
    fn nsnp(&self) -> usize {
        self.nsnp
    }
    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }
    fn record_bytes(&self) -> usize {
        self.am_rec_bytes
    }

    fn read_record(&mut self, dst: &mut [u8]) -> Result<bool> {
        if self.next_idx >= self.nsnp {
            return Ok(false);
        }
        if dst.len() != self.am_rec_bytes {
            bail!(
                "dst len {} != record_bytes {}",
                dst.len(),
                self.am_rec_bytes
            );
        }
        let start = 3 + self.plink_rec_bytes * self.next_idx;
        let end = start + self.plink_rec_bytes;
        let src = &self.mmap[start..end];

        // Translate byte-at-a-time via LUT.
        for (i, &b) in src.iter().enumerate() {
            dst[i] = PLINK_TO_AM[b as usize];
        }

        // Mask trailing padding bits in the last byte to zero (AM canonical:
        // unused sample positions = 0 bits). Only when nind % 4 != 0.
        if self.nind % 4 != 0 {
            let valid = self.last_byte_valid_samples;
            // Keep high `valid * 2` bits (MSB-first), zero the rest.
            let keep_bits = valid * 2;
            let mask = 0xFFu8 << (8 - keep_bits);
            let last = self.am_rec_bytes - 1;
            dst[last] &= mask;
        }

        self.next_idx += 1;
        Ok(true)
    }
}

// ======================================================================
// Writer
// ======================================================================

pub struct PackedPedWriter {
    w: BufWriter<File>,
    nind: usize,
    nsnp: usize,
    rec_bytes: usize,
    records_written: usize,
    began: bool,
    /// Scratch buffer for recoded output.
    out_buf: Vec<u8>,
}

impl PackedPedWriter {
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
        let w = BufWriter::with_capacity(256 * 1024, file);
        Ok(Self {
            w,
            nind: 0,
            nsnp: 0,
            rec_bytes: 0,
            records_written: 0,
            began: false,
            out_buf: Vec::new(),
        })
    }
}

impl GenoWriter for PackedPedWriter {
    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }

    fn begin(&mut self, nind: usize, nsnp: usize, _ihash: u32, _shash: u32) -> Result<()> {
        if self.began {
            bail!("PackedPedWriter::begin called twice");
        }
        self.nind = nind;
        self.nsnp = nsnp;
        self.rec_bytes = (nind + 3) / 4;
        self.out_buf = vec![0u8; self.rec_bytes];
        self.w.write_all(&BED_MAGIC)?;
        self.began = true;
        Ok(())
    }

    fn write_record(&mut self, src: &[u8]) -> Result<()> {
        if !self.began {
            bail!("write_record before begin");
        }
        if src.len() != self.rec_bytes {
            bail!("record len {} != expected {}", src.len(), self.rec_bytes);
        }
        if self.records_written >= self.nsnp {
            bail!("too many records: expected {}", self.nsnp);
        }

        // Recode canonical AM bytes to PLINK bytes via LUT.
        for (i, &b) in src.iter().enumerate() {
            self.out_buf[i] = AM_TO_PLINK[b as usize];
        }

        // Pad unused sample slots in the last byte with PLINK "missing" bits
        // (01). In AM canonical those slots are 0 bits (AM code = 0), which
        // the LUT translates to PLINK code `11` (hom A2). We must overwrite
        // the unused lanes to `01` instead.
        if self.nind % 4 != 0 {
            let valid = ((self.nind - 1) % 4) + 1; // 1..=3 valid samples in last byte
            let last = self.rec_bytes - 1;
            let mut last_byte = self.out_buf[last];
            // PLINK is LSB-first; valid samples occupy bits [0 .. 2*valid).
            // Zero the high (pad) bits, then OR in 01 at each pad slot.
            let valid_mask = (1u8 << (valid * 2)) - 1;
            last_byte &= valid_mask;
            for slot in valid..4 {
                last_byte |= 0b01 << (slot * 2);
            }
            self.out_buf[last] = last_byte;
        }

        self.w.write_all(&self.out_buf)?;
        self.records_written += 1;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        if self.records_written != self.nsnp {
            bail!(
                "finish: wrote {} records, expected {}",
                self.records_written,
                self.nsnp
            );
        }
        self.w.flush()?;
        Ok(())
    }
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geno::codec;

    #[test]
    fn lut_inverse() {
        // PLINK → AM → PLINK should be identity.
        for i in 0..=255u8 {
            assert_eq!(
                AM_TO_PLINK[PLINK_TO_AM[i as usize] as usize], i,
                "roundtrip failed for byte {i:08b}"
            );
        }
    }

    #[test]
    fn recode_single_samples() {
        // PLINK hom A1 (00) → AM g=2 (10)
        // PLINK missing (01) → AM missing (11)
        // PLINK het    (10) → AM g=1 (01)
        // PLINK hom A2 (11) → AM g=0 (00)
        assert_eq!(recode_plink_to_am(0b00), 0b10);
        assert_eq!(recode_plink_to_am(0b01), 0b11);
        assert_eq!(recode_plink_to_am(0b10), 0b01);
        assert_eq!(recode_plink_to_am(0b11), 0b00);
    }

    #[test]
    fn byte_recode_with_bit_reverse() {
        // 4 samples packed PLINK LSB-first: s0=hom_A2(11), s1=het(10),
        // s2=missing(01), s3=hom_A1(00)
        // Byte = 0b00_01_10_11
        let plink_byte = 0b00_01_10_11u8;
        let am_byte = PLINK_TO_AM[plink_byte as usize];
        // Expected AM MSB-first: s0=0(00), s1=1(01), s2=miss(11), s3=2(10)
        // Byte = 0b00_01_11_10
        assert_eq!(am_byte, 0b00_01_11_10);
    }

    fn ref_bed(path: &Path, nind: usize, nsnp: usize, genotypes: &[Vec<u8>]) {
        // Build a reference .bed file by hand, encoding genotypes PLINK-style.
        // genotypes[s][i] = 0/1/2/9 (AM encoding) for SNP s, sample i.
        use std::io::Write;
        let mut f = File::create(path).unwrap();
        f.write_all(&BED_MAGIC).unwrap();
        let rec_bytes = (nind + 3) / 4;
        for snp in 0..nsnp {
            let mut out = vec![0u8; rec_bytes];
            for i in 0..nind {
                let am = genotypes[snp][i];
                let plink_code = match am {
                    0 => 0b11, // hom A2
                    1 => 0b10, // het
                    2 => 0b00, // hom A1
                    _ => 0b01, // missing
                };
                let byte = i / 4;
                let shift = 2 * (i % 4); // PLINK LSB-first
                out[byte] |= plink_code << shift;
            }
            // Pad trailing sample slots with missing (01).
            if nind % 4 != 0 {
                let valid = ((nind - 1) % 4) + 1;
                let last = rec_bytes - 1;
                for slot in valid..4 {
                    out[last] |= 0b01 << (slot * 2);
                }
            }
            f.write_all(&out).unwrap();
        }
    }

    #[test]
    fn reader_decodes_reference_bed() {
        let nind = 7; // not multiple of 4 — exercises padding
        let nsnp = 3;
        let gs: Vec<Vec<u8>> = vec![
            vec![0, 1, 2, 9, 0, 1, 2],
            vec![2, 2, 2, 2, 9, 9, 9],
            vec![0, 9, 1, 2, 0, 9, 1],
        ];
        let tmp = tempfile::NamedTempFile::new().unwrap();
        ref_bed(tmp.path(), nind, nsnp, &gs);

        let mut r = PackedPedReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        let mut unpacked = vec![0u8; nind];
        for expected in &gs {
            assert!(r.read_record(&mut buf).unwrap());
            codec::unpack(&buf, nind, &mut unpacked);
            assert_eq!(&unpacked, expected);
        }
        assert!(!r.read_record(&mut buf).unwrap());
    }

    #[test]
    fn writer_reader_roundtrip() {
        // Write via PackedPedWriter, read back, compare.
        let nind = 11; // partial tail
        let nsnp = 20;
        let records: Vec<Vec<u8>> = (0..nsnp)
            .map(|i| (0..nind).map(|j| [0u8, 1, 2, 9][(i + j * 3) % 4]).collect())
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = PackedPedWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for rec in &records {
            let mut p = vec![0u8; (nind * 2 + 7) / 8];
            codec::pack(rec, &mut p);
            w.write_record(&p).unwrap();
        }
        w.finish().unwrap();

        let mut r = PackedPedReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        let mut unpacked = vec![0u8; nind];
        for expected in &records {
            assert!(r.read_record(&mut buf).unwrap());
            codec::unpack(&buf, nind, &mut unpacked);
            assert_eq!(&unpacked, expected);
        }
        assert!(!r.read_record(&mut buf).unwrap());
    }

    #[test]
    fn set_next_idx_seeks_to_correct_record() {
        let nind = 11;
        let nsnp = 8;
        let records: Vec<Vec<u8>> = (0..nsnp)
            .map(|i| (0..nind).map(|j| [0u8, 1, 2, 9][(i * 5 + j) % 4]).collect())
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = PackedPedWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for rec in &records {
            let mut p = vec![0u8; (nind * 2 + 7) / 8];
            codec::pack(rec, &mut p);
            w.write_record(&p).unwrap();
        }
        w.finish().unwrap();

        let mut r = PackedPedReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        let mut unpacked = vec![0u8; nind];
        for &idx in &[4usize, 0, 7, 3, 1, 6, 2, 5] {
            r.set_next_idx(idx);
            assert!(r.read_record(&mut buf).unwrap());
            codec::unpack(&buf, nind, &mut unpacked);
            assert_eq!(&unpacked, &records[idx], "mismatch at seek idx {idx}");
        }
    }

    #[test]
    fn rejects_bad_magic() {
        use std::io::Write;
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut f = File::create(tmp.path()).unwrap();
        f.write_all(&[0x00, 0x00, 0x00]).unwrap();
        match PackedPedReader::open(tmp.path(), 4, 1) {
            Ok(_) => panic!("expected magic error"),
            Err(e) => assert!(format!("{e:#}").contains("magic")),
        }
    }

    #[test]
    fn rejects_sample_major() {
        use std::io::Write;
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut f = File::create(tmp.path()).unwrap();
        f.write_all(&[0x6c, 0x1b, 0x00]).unwrap(); // sample-major
        match PackedPedReader::open(tmp.path(), 4, 1) {
            Ok(_) => panic!("expected sample-major rejection"),
            Err(e) => assert!(format!("{e:#}").to_lowercase().contains("sample-major")),
        }
    }
}
