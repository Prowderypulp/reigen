//! EIGENSTRAT text genotype file.
//!
//! - SNP-major. One line per SNP, one ASCII char per sample, `\n` terminator.
//! - Chars: `'0'`, `'1'`, `'2'`, `'9'` (missing). No header record.
//!
//! # Writer (this phase)
//!
//! For each SNP record (2-bit packed, canonical MSB-first), decode to an
//! ASCII line via a 256-entry byte→4-char LUT: each input byte holds 4
//! genotypes (bit pairs from MSB to LSB), and maps to a fixed 4-byte ASCII
//! string. We emit those 4 bytes per input byte, then trim trailing padding
//! to exactly `nind` chars, append `\n`, and `write_all`.
//!
//! `BufWriter` with a 256 KB buffer keeps syscalls rare even for
//! millions of SNPs.
//!
//! # Reader
//!
//! Deferred to Phase 2 (when we also do EIGENSTRAT→PAM). The reader stub
//! here returns an error to keep the trait impl present for dispatch code.

use super::{GenoReader, GenoWriter, Layout};
use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

// ==================================================================
// LUT: byte (4 packed 2-bit genotypes, MSB-first) → 4 ASCII bytes
// ==================================================================

/// For each of 256 possible input bytes, precompute the 4 output ASCII chars
/// in MSB-first order (sample 0 = bits 7-6, sample 1 = bits 5-4, etc.).
static BYTE_TO_ASCII4: [[u8; 4]; 256] = build_byte_lut();

const fn build_byte_lut() -> [[u8; 4]; 256] {
    let mut table = [[0u8; 4]; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        // MSB-first: bits (7-6), (5-4), (3-2), (1-0).
        table[i][0] = two_bit_to_ascii((b >> 6) & 0b11);
        table[i][1] = two_bit_to_ascii((b >> 4) & 0b11);
        table[i][2] = two_bit_to_ascii((b >> 2) & 0b11);
        table[i][3] = two_bit_to_ascii(b & 0b11);
        i += 1;
    }
    table
}

const fn two_bit_to_ascii(t: u8) -> u8 {
    match t & 0b11 {
        0 => b'0',
        1 => b'1',
        2 => b'2',
        _ => b'9', // 3 = missing
    }
}

// ==================================================================
// Writer
// ==================================================================

pub struct EigenstratWriter {
    w: BufWriter<File>,
    nind: usize,
    nsnp: usize,
    record_bytes: usize,
    /// Reused per-record scratch: `nind + 1` bytes (ASCII + `\n`).
    line_buf: Vec<u8>,
    records_written: usize,
    began: bool,
}

impl EigenstratWriter {
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
        Ok(Self {
            w: BufWriter::with_capacity(256 * 1024, file),
            nind: 0,
            nsnp: 0,
            record_bytes: 0,
            line_buf: Vec::new(),
            records_written: 0,
            began: false,
        })
    }
}

impl GenoWriter for EigenstratWriter {
    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }

    fn begin(&mut self, nind: usize, nsnp: usize, _ihash: u32, _shash: u32) -> Result<()> {
        if self.began {
            bail!("EigenstratWriter::begin called twice");
        }
        self.nind = nind;
        self.nsnp = nsnp;
        self.record_bytes = (nind * 2 + 7) / 8;
        self.line_buf = vec![0u8; nind + 1];
        self.line_buf[nind] = b'\n';
        self.began = true;
        Ok(())
    }

    fn write_record(&mut self, src: &[u8]) -> Result<()> {
        if !self.began {
            bail!("write_record before begin");
        }
        if src.len() != self.record_bytes {
            bail!("record len {} != expected {}", src.len(), self.record_bytes);
        }
        if self.records_written >= self.nsnp {
            bail!("too many records: expected {}", self.nsnp);
        }

        // Decode full bytes first (4 samples each), then the tail.
        let full_bytes = self.nind / 4;
        let tail_samples = self.nind % 4;

        for i in 0..full_bytes {
            let chars = &BYTE_TO_ASCII4[src[i] as usize];
            let dst = &mut self.line_buf[i * 4..(i + 1) * 4];
            dst.copy_from_slice(chars);
        }

        if tail_samples > 0 {
            let chars = &BYTE_TO_ASCII4[src[full_bytes] as usize];
            let dst_start = full_bytes * 4;
            self.line_buf[dst_start..dst_start + tail_samples]
                .copy_from_slice(&chars[..tail_samples]);
        }

        // `line_buf[nind]` is already `\n` from begin().
        self.w.write_all(&self.line_buf)?;
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

// ==================================================================
// Reader
// ==================================================================

/// LUT: ASCII byte → 2-bit genotype. Non-{0,1,2} → 0b11 (missing).
static ASCII_TO_2BIT: [u8; 256] = build_ascii_lut();

const fn build_ascii_lut() -> [u8; 256] {
    let mut t = [0b11u8; 256]; // default: missing
    t[b'0' as usize] = 0b00;
    t[b'1' as usize] = 0b01;
    t[b'2' as usize] = 0b10;
    // b'9' stays missing (0b11) — same encoding.
    t
}

pub struct EigenstratReader {
    mmap: memmap2::Mmap,
    nind: usize,
    nsnp: usize,
    record_bytes: usize,
    /// Byte offset of start of each SNP line. Length `nsnp`. Trailing `\n`
    /// is at offset `line_starts[i] + nind`.
    line_starts: Vec<usize>,
    next_idx: usize,
}

impl EigenstratReader {
    pub fn open(path: &Path, nind: usize, nsnp: usize) -> Result<Self> {
        let file = std::fs::File::open(path).with_context(|| format!("open {}", path.display()))?;
        if file.metadata()?.len() == 0 {
            if nsnp == 0 {
                return Ok(Self {
                    mmap: unsafe { memmap2::Mmap::map(&file)? },
                    nind,
                    nsnp,
                    record_bytes: 0,
                    line_starts: Vec::new(),
                    next_idx: 0,
                });
            }
            bail!(
                "EIGENSTRAT {}: empty file but .snp has {nsnp} rows",
                path.display()
            );
        }
        let mmap = unsafe { memmap2::Mmap::map(&file) }
            .with_context(|| format!("mmap {}", path.display()))?;

        // Index line starts in one pass. Each line MUST be exactly nind ASCII
        // chars + '\n' (or '\r\n'). Upstream convertf doesn't tolerate ragged
        // lines either, so we match.
        let line_starts = index_lines(&mmap, nind, nsnp)
            .with_context(|| format!("index lines in {}", path.display()))?;

        let record_bytes = (nind * 2 + 7) / 8;
        Ok(Self {
            mmap,
            nind,
            nsnp,
            record_bytes,
            line_starts,
            next_idx: 0,
        })
    }
}

/// Scan mmap bytes and return one offset per SNP line. Verifies each line
/// has exactly `nind` sample chars.
fn index_lines(bytes: &[u8], nind: usize, expected_nsnp: usize) -> Result<Vec<usize>> {
    let mut starts = Vec::with_capacity(expected_nsnp);
    let mut cursor = 0usize;
    let len = bytes.len();

    while cursor < len {
        // Skip blank lines defensively — upstream convertf doesn't emit them
        // but files produced by hand editing might have them.
        if bytes[cursor] == b'\n' {
            cursor += 1;
            continue;
        }
        if bytes[cursor] == b'\r' && cursor + 1 < len && bytes[cursor + 1] == b'\n' {
            cursor += 2;
            continue;
        }

        let line_start = cursor;
        let nl = memchr::memchr(b'\n', &bytes[cursor..])
            .map(|off| cursor + off)
            .unwrap_or(len);

        // Allow optional trailing '\r'.
        let content_end = if nl > line_start && bytes[nl.saturating_sub(1)] == b'\r' {
            nl - 1
        } else {
            nl
        };
        let content_len = content_end - line_start;

        if content_len != nind {
            bail!(
                "EIGENSTRAT line {} has {} chars, expected {} (one per sample)",
                starts.len() + 1,
                content_len,
                nind
            );
        }
        if starts.len() >= expected_nsnp {
            bail!("EIGENSTRAT has more lines than .snp has rows ({expected_nsnp})");
        }
        starts.push(line_start);
        cursor = nl + 1; // past '\n', or past EOF
    }

    if starts.len() != expected_nsnp {
        bail!(
            "EIGENSTRAT line count {} != .snp row count {}",
            starts.len(),
            expected_nsnp
        );
    }
    Ok(starts)
}

impl GenoReader for EigenstratReader {
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
        self.record_bytes
    }

    fn read_record(&mut self, dst: &mut [u8]) -> Result<bool> {
        if self.next_idx >= self.nsnp {
            return Ok(false);
        }
        if dst.len() != self.record_bytes {
            bail!(
                "dst len {} != record_bytes {}",
                dst.len(),
                self.record_bytes
            );
        }
        let start = self.line_starts[self.next_idx];
        let end = start + self.nind;
        let line = &self.mmap[start..end];

        // Pack ASCII line into canonical 2-bit MSB-first. Process 4 chars
        // per output byte.
        for b in dst.iter_mut() {
            *b = 0;
        }
        let full = self.nind / 4;
        for i in 0..full {
            let off = i * 4;
            let b = (ASCII_TO_2BIT[line[off] as usize] << 6)
                | (ASCII_TO_2BIT[line[off + 1] as usize] << 4)
                | (ASCII_TO_2BIT[line[off + 2] as usize] << 2)
                | ASCII_TO_2BIT[line[off + 3] as usize];
            dst[i] = b;
        }
        let tail = self.nind % 4;
        if tail > 0 {
            let off = full * 4;
            let mut b = 0u8;
            for k in 0..tail {
                b |= ASCII_TO_2BIT[line[off + k] as usize] << (6 - 2 * k);
            }
            // Remaining bits of `dst[full]` stay 0 (padding).
            dst[full] = b;
        }

        self.next_idx += 1;
        Ok(true)
    }
}

// ==================================================================
// Tests
// ==================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geno::codec;
    use std::io::Read;

    #[test]
    fn lut_correctness_sample() {
        // 0b00_01_10_11 = genotypes 0,1,2,missing → "0129"
        assert_eq!(&BYTE_TO_ASCII4[0b00_01_10_11], b"0129");
        // All zeros → "0000"
        assert_eq!(&BYTE_TO_ASCII4[0], b"0000");
        // All missing → "9999"
        assert_eq!(&BYTE_TO_ASCII4[0xFF], b"9999");
    }

    #[test]
    fn writes_small_record() {
        // nind=5: 5 samples = "01290\n"
        let nind = 5;
        let nsnp = 1;
        let genotypes = vec![0u8, 1, 2, 9, 0];
        let mut packed = vec![0u8; (nind * 2 + 7) / 8];
        codec::pack(&genotypes, &mut packed);

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = EigenstratWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        w.write_record(&packed).unwrap();
        w.finish().unwrap();

        let mut out = String::new();
        std::fs::File::open(tmp.path())
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "01290\n");
    }

    #[test]
    fn writes_exact_multiple_of_four() {
        let nind = 8;
        let nsnp = 2;
        let r0: Vec<u8> = vec![0, 1, 2, 9, 0, 1, 2, 9];
        let r1: Vec<u8> = vec![2, 2, 1, 1, 0, 0, 9, 9];
        let pack = |gs: &[u8]| -> Vec<u8> {
            let mut p = vec![0u8; (gs.len() * 2 + 7) / 8];
            codec::pack(gs, &mut p);
            p
        };

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = EigenstratWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        w.write_record(&pack(&r0)).unwrap();
        w.write_record(&pack(&r1)).unwrap();
        w.finish().unwrap();

        let out = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(out, "01290129\n22110099\n");
    }

    #[test]
    fn writes_with_partial_tail() {
        // nind=7 → 1 full byte (4 samples) + 3 tail.
        let nind = 7;
        let nsnp = 1;
        let gs: Vec<u8> = vec![9, 2, 1, 0, 1, 2, 9];
        let mut p = vec![0u8; (nind * 2 + 7) / 8];
        codec::pack(&gs, &mut p);

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = EigenstratWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        w.write_record(&p).unwrap();
        w.finish().unwrap();

        let out = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(out, "9210129\n");
    }

    // --- Reader tests ---

    #[test]
    fn ascii_lut_correctness() {
        assert_eq!(ASCII_TO_2BIT[b'0' as usize], 0b00);
        assert_eq!(ASCII_TO_2BIT[b'1' as usize], 0b01);
        assert_eq!(ASCII_TO_2BIT[b'2' as usize], 0b10);
        assert_eq!(ASCII_TO_2BIT[b'9' as usize], 0b11);
        assert_eq!(ASCII_TO_2BIT[b'X' as usize], 0b11); // garbage → missing
    }

    fn write_tmp_text(s: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut f, s.as_bytes()).unwrap();
        f
    }

    #[test]
    fn reader_basic() {
        let f = write_tmp_text("01290\n22110\n99221\n");
        let mut r = EigenstratReader::open(f.path(), 5, 3).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        let mut unpacked = vec![0u8; 5];

        assert!(r.read_record(&mut buf).unwrap());
        codec::unpack(&buf, 5, &mut unpacked);
        assert_eq!(unpacked, vec![0, 1, 2, 9, 0]);

        assert!(r.read_record(&mut buf).unwrap());
        codec::unpack(&buf, 5, &mut unpacked);
        assert_eq!(unpacked, vec![2, 2, 1, 1, 0]);

        assert!(r.read_record(&mut buf).unwrap());
        codec::unpack(&buf, 5, &mut unpacked);
        assert_eq!(unpacked, vec![9, 9, 2, 2, 1]);

        assert!(!r.read_record(&mut buf).unwrap());
    }

    #[test]
    fn reader_handles_crlf() {
        let f = write_tmp_text("01290\r\n22110\r\n");
        let mut r = EigenstratReader::open(f.path(), 5, 2).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        assert!(r.read_record(&mut buf).unwrap());
        assert!(r.read_record(&mut buf).unwrap());
        assert!(!r.read_record(&mut buf).unwrap());
    }

    #[test]
    fn reader_rejects_ragged_line() {
        let f = write_tmp_text("01290\n2211\n99221\n"); // line 2 has 4 chars
        match EigenstratReader::open(f.path(), 5, 3) {
            Ok(_) => panic!("expected ragged-line error"),
            Err(e) => assert!(format!("{e:#}").contains("expected 5")),
        }
    }

    #[test]
    fn reader_rejects_line_count_mismatch() {
        let f = write_tmp_text("01290\n");
        match EigenstratReader::open(f.path(), 5, 3) {
            Ok(_) => panic!("expected count mismatch"),
            Err(e) => assert!(format!("{e:#}").contains("line count")),
        }
    }

    #[test]
    fn writer_reader_roundtrip() {
        // Write via EigenstratWriter, read via EigenstratReader, compare.
        let nind = 7; // not multiple of 4 — exercises tail
        let nsnp = 10;
        let records: Vec<Vec<u8>> = (0..nsnp)
            .map(|i| {
                (0..nind)
                    .map(|j| [0u8, 1, 2, 9][(i + j * 3) % 4])
                    .collect::<Vec<_>>()
            })
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = EigenstratWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for rec in &records {
            let mut p = vec![0u8; (nind * 2 + 7) / 8];
            codec::pack(rec, &mut p);
            w.write_record(&p).unwrap();
        }
        w.finish().unwrap();

        let mut r = EigenstratReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        let mut unpacked = vec![0u8; nind];
        for expected in &records {
            assert!(r.read_record(&mut buf).unwrap());
            codec::unpack(&buf, nind, &mut unpacked);
            assert_eq!(&unpacked, expected);
        }
        assert!(!r.read_record(&mut buf).unwrap());
    }
}
