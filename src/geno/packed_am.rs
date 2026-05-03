//! PACKEDANCESTRYMAP format.
//!
//! # On-disk layout (from mcio.c `outpack` / `inpack2`)
//!
//! - SNP-major. One record per SNP.
//! - `rlen = max(48, ceil(nind * 2 / 8))` bytes per record.
//! - First record is a header (padded with zeros to `rlen`):
//!     `sprintf(buff, "GENO %d %d %x %x", nind, nsnp, ihash, shash)`
//! - Subsequent `nsnp` records: 2-bit MSB-first, canonical encoding
//!   (`00=0, 01=1, 10=2, 11=missing`).
//!
//! # Read strategy
//!
//! `mmap` the whole file. Header is bytes `[0..rlen]`; each SNP record is
//! `[rlen + i*rlen .. rlen + (i+1)*rlen]`. The in-memory canonical buffer
//! is `ceil(nind*2/8)` bytes — the record's trailing padding (when nind is
//! small enough that the real record is < 48 bytes) is dropped on read.
//!
//! # Write strategy
//!
//! Build header ASCII into a zero-padded `rlen`-byte buffer. For each SNP,
//! write the `record_bytes` canonical payload + zero padding up to `rlen`.
//! Single `BufWriter<File>` with 256 KB buffer.

use super::{GenoReader, GenoWriter, Layout};
use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub const MIN_RLEN: usize = 48;
const HEADER_MAGIC: &[u8] = b"GENO ";

#[inline]
pub fn rlen_for(nind: usize) -> usize {
    std::cmp::max(MIN_RLEN, (nind * 2 + 7) / 8)
}

// ==================================================================
// Reader
// ==================================================================

pub struct PackedAmReader {
    mmap: Mmap,
    nind: usize,
    nsnp: usize,
    rlen: usize,
    record_bytes: usize,
    next_idx: usize,
    /// Parsed from header; 0/0 when hashes were never written (convertf-rs
    /// placeholder output). Callers can consult these when hashcheck: YES.
    pub ihash: u32,
    pub shash: u32,
}

impl PackedAmReader {
    /// Jump the streaming cursor to SNP index `idx`. The next `read_record`
    /// will produce that SNP. Used by merge for random-access reads.
    pub fn set_next_idx(&mut self, idx: usize) {
        self.next_idx = idx;
    }

    /// Open a PACKEDANCESTRYMAP `.geno` file. `nind` and `nsnp` must come
    /// from the companion `.ind` / `.snp` files — they are used to verify
    /// the header and compute `rlen`.
    pub fn open(path: &Path, nind: usize, nsnp: usize) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mmap =
            unsafe { Mmap::map(&file) }.with_context(|| format!("mmap {}", path.display()))?;

        let rlen = rlen_for(nind);
        let expected_len = rlen * (nsnp + 1);
        if mmap.len() < expected_len {
            bail!(
                "PACKEDANCESTRYMAP {}: file size {} < expected {} (rlen={}, nsnp={}, nind={})",
                path.display(),
                mmap.len(),
                expected_len,
                rlen,
                nsnp,
                nind
            );
        }
        if mmap.len() > expected_len {
            // Upstream tolerates trailing bytes silently; we warn.
            log::warn!(
                "PACKEDANCESTRYMAP {}: file has {} trailing bytes past expected end",
                path.display(),
                mmap.len() - expected_len
            );
        }

        let (hdr_nind, hdr_nsnp, ihash, shash) = parse_header(&mmap[..rlen])
            .with_context(|| format!("parse header of {}", path.display()))?;

        if hdr_nind != nind {
            bail!("PACKEDANCESTRYMAP header nind {hdr_nind} != .ind count {nind}");
        }
        if hdr_nsnp != nsnp {
            bail!("PACKEDANCESTRYMAP header nsnp {hdr_nsnp} != .snp count {nsnp}");
        }

        let record_bytes = (nind * 2 + 7) / 8;

        Ok(Self {
            mmap,
            nind,
            nsnp,
            rlen,
            record_bytes,
            next_idx: 0,
            ihash,
            shash,
        })
    }
}

impl GenoReader for PackedAmReader {
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
    fn header_hashes(&self) -> Option<(u32, u32)> {
        Some((self.ihash, self.shash))
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
        let start = self.rlen * (1 + self.next_idx); // skip header
        let end = start + self.record_bytes;
        dst.copy_from_slice(&self.mmap[start..end]);
        self.next_idx += 1;
        Ok(true)
    }
}

/// Parse `"GENO %d %d %x %x"` from the header record.
fn parse_header(hdr: &[u8]) -> Result<(usize, usize, u32, u32)> {
    if !hdr.starts_with(HEADER_MAGIC) {
        bail!("not a PACKEDANCESTRYMAP file (missing 'GENO ' magic)");
    }
    // Strip trailing NULs from the zero-padded ASCII header.
    let end = hdr.iter().position(|&b| b == 0).unwrap_or(hdr.len());
    let text = std::str::from_utf8(&hdr[..end]).map_err(|e| anyhow!("header not UTF-8: {e}"))?;
    let mut it = text.split_ascii_whitespace();
    let _geno = it.next(); // "GENO"
    let nind: usize = it
        .next()
        .ok_or_else(|| anyhow!("header: missing nind"))?
        .parse()
        .map_err(|e| anyhow!("header nind: {e}"))?;
    let nsnp: usize = it
        .next()
        .ok_or_else(|| anyhow!("header: missing nsnp"))?
        .parse()
        .map_err(|e| anyhow!("header nsnp: {e}"))?;
    let ihash: u32 = u32::from_str_radix(
        it.next().ok_or_else(|| anyhow!("header: missing ihash"))?,
        16,
    )
    .map_err(|e| anyhow!("header ihash: {e}"))?;
    let shash: u32 = u32::from_str_radix(
        it.next().ok_or_else(|| anyhow!("header: missing shash"))?,
        16,
    )
    .map_err(|e| anyhow!("header shash: {e}"))?;
    Ok((nind, nsnp, ihash, shash))
}

// ==================================================================
// Writer
// ==================================================================

pub struct PackedAmWriter {
    w: BufWriter<File>,
    nind: usize,
    nsnp: usize,
    rlen: usize,
    record_bytes: usize,
    records_written: usize,
    header_written: bool,
}

impl PackedAmWriter {
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
        let w = BufWriter::with_capacity(256 * 1024, file);
        Ok(Self {
            w,
            nind: 0,
            nsnp: 0,
            rlen: 0,
            record_bytes: 0,
            records_written: 0,
            header_written: false,
        })
    }
}

impl GenoWriter for PackedAmWriter {
    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }

    fn begin(&mut self, nind: usize, nsnp: usize, ihash: u32, shash: u32) -> Result<()> {
        if self.header_written {
            bail!("PackedAmWriter::begin called twice");
        }
        self.nind = nind;
        self.nsnp = nsnp;
        self.rlen = rlen_for(nind);
        self.record_bytes = (nind * 2 + 7) / 8;

        let mut buf = vec![0u8; self.rlen];
        let header = format!("GENO {} {} {:x} {:x}", nind, nsnp, ihash, shash);
        if header.len() > self.rlen {
            bail!(
                "PACKEDANCESTRYMAP header {header:?} exceeds rlen {}",
                self.rlen
            );
        }
        buf[..header.len()].copy_from_slice(header.as_bytes());
        self.w.write_all(&buf)?;
        self.header_written = true;
        Ok(())
    }

    fn write_record(&mut self, src: &[u8]) -> Result<()> {
        if !self.header_written {
            bail!("write_record before begin");
        }
        if src.len() != self.record_bytes {
            bail!("record len {} != expected {}", src.len(), self.record_bytes);
        }
        if self.records_written >= self.nsnp {
            bail!("too many records written: expected {} SNPs", self.nsnp);
        }
        self.w.write_all(src)?;
        if self.rlen > self.record_bytes {
            // Pad to rlen with zeros (only when nind is so small that
            // real record < 48 bytes).
            let pad = self.rlen - self.record_bytes;
            let zeros = [0u8; 64];
            let mut remaining = pad;
            while remaining > 0 {
                let take = remaining.min(zeros.len());
                self.w.write_all(&zeros[..take])?;
                remaining -= take;
            }
        }
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
// Tests
// ==================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geno::codec;

    #[test]
    fn rlen_min_48() {
        assert_eq!(rlen_for(1), 48);
        assert_eq!(rlen_for(192), 48);
        assert_eq!(rlen_for(193), 49);
        assert_eq!(rlen_for(4000), (4000 * 2 + 7) / 8);
    }

    #[test]
    fn parses_header() {
        let mut buf = vec![0u8; 48];
        let s = b"GENO 5 3 1a2b cafe";
        buf[..s.len()].copy_from_slice(s);
        let (nind, nsnp, ih, sh) = parse_header(&buf).unwrap();
        assert_eq!((nind, nsnp, ih, sh), (5, 3, 0x1a2b, 0xcafe));
    }

    #[test]
    fn rejects_wrong_magic() {
        let buf = vec![b'X'; 48];
        assert!(parse_header(&buf).is_err());
    }

    #[test]
    fn writer_reader_roundtrip_small_nind() {
        // nind=5 → rec_bytes = ceil(10/8) = 2, but rlen = 48.
        let nind = 5;
        let nsnp = 3;
        let records: Vec<Vec<u8>> = (0..nsnp)
            .map(|i| {
                let gs: Vec<u8> = (0..nind).map(|j| [0u8, 1, 2, 9][(i + j) % 4]).collect();
                let mut packed = vec![0u8; (nind * 2 + 7) / 8];
                codec::pack(&gs, &mut packed);
                packed
            })
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = PackedAmWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for r in &records {
            w.write_record(r).unwrap();
        }
        w.finish().unwrap();

        // File size: rlen * (1 + nsnp) = 48 * 4 = 192.
        let meta = std::fs::metadata(tmp.path()).unwrap();
        assert_eq!(meta.len(), 48 * 4);

        // Read back.
        let mut r = PackedAmReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        for expected in &records {
            assert!(r.read_record(&mut buf).unwrap());
            assert_eq!(&buf, expected);
        }
        assert!(!r.read_record(&mut buf).unwrap()); // EOF
    }

    #[test]
    fn writer_reader_roundtrip_large_nind() {
        // nind=300 → rec_bytes = 75, which is > 48 → rlen = 75.
        let nind = 300;
        let nsnp = 10;
        let rec = (nind * 2 + 7) / 8;
        assert!(rec > MIN_RLEN);

        let records: Vec<Vec<u8>> = (0..nsnp)
            .map(|i| {
                let gs: Vec<u8> = (0..nind).map(|j| [0u8, 1, 2, 9][(i + j) % 4]).collect();
                let mut packed = vec![0u8; rec];
                codec::pack(&gs, &mut packed);
                packed
            })
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = PackedAmWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for r in &records {
            w.write_record(r).unwrap();
        }
        w.finish().unwrap();

        let mut r = PackedAmReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        for expected in &records {
            assert!(r.read_record(&mut buf).unwrap());
            assert_eq!(&buf, expected);
        }
        assert!(!r.read_record(&mut buf).unwrap());
    }

    #[test]
    fn set_next_idx_seeks_to_correct_record() {
        // 10 samples (small record, triggers 48-byte padding) × 8 SNPs.
        // Write distinct payloads; seek to arbitrary indices and confirm
        // read_record returns the right one, not the sequentially-next one.
        let nind = 10;
        let nsnp = 8;
        let records: Vec<Vec<u8>> = (0..nsnp)
            .map(|i| {
                let gs: Vec<u8> = (0..nind).map(|j| [0u8, 1, 2, 9][(i * 7 + j) % 4]).collect();
                let mut packed = vec![0u8; (nind * 2 + 7) / 8];
                codec::pack(&gs, &mut packed);
                packed
            })
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = PackedAmWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for r in &records {
            w.write_record(r).unwrap();
        }
        w.finish().unwrap();

        let mut r = PackedAmReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];

        // Out-of-order seeks, including backward.
        for &idx in &[3usize, 0, 7, 5, 1, 6, 2, 4] {
            r.set_next_idx(idx);
            assert!(r.read_record(&mut buf).unwrap());
            assert_eq!(&buf, &records[idx], "mismatch at seek idx {idx}");
        }
    }

    #[test]
    fn reader_detects_bad_header_counts() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        {
            let mut w = PackedAmWriter::create(tmp.path()).unwrap();
            w.begin(10, 5, 0, 0).unwrap();
            for _ in 0..5 {
                w.write_record(&vec![0u8; 3]).unwrap();
            }
            w.finish().unwrap();
        }
        // Open with wrong nind — PackedAmReader lacks Debug, so we match.
        match PackedAmReader::open(tmp.path(), 11, 5) {
            Ok(_) => panic!("expected error for mismatched nind"),
            Err(e) => {
                let msg = format!("{e:#}");
                assert!(msg.contains("nind"), "error did not mention nind: {msg}");
            }
        }
    }
}
