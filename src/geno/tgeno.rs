//! TGENO (TRANSPOSE_GENO) format — sample-major packed.
//!
//! # On-disk layout (v8.0.0+ of AdmixTools)
//!
//! - Sample-major. One record per sample.
//! - `rlen = max(48, ceil(nsnp * 2 / 8))` bytes per record.
//! - First record is a header (zero-padded to `rlen`):
//!     `"TGENO %d %d %x %x"` where fields are nind, nsnp, ihash, shash.
//! - Encoding: canonical 2-bit MSB-first (same as PACKEDANCESTRYMAP).
//!
//! # Why it exists
//!
//! PACKEDANCESTRYMAP has a 48-byte minimum record length. For datasets with
//! few samples (say, nind=10), the real payload is 3 bytes but the file
//! stores 48, wasting 94%. TGENO swaps axes so the per-record payload is
//! nsnp bits (lots) instead of nind bits (few). AADR-scale datasets are
//! large in both directions; TGENO is for small cohorts.
//!
//! # Conversion to/from SnpMajor formats
//!
//! TGENO ↔ PACKEDANCESTRYMAP requires a full matrix transpose of a 2-bit
//! matrix. See `transpose.rs`. The pipeline materializes both matrices
//! fully in memory — acceptable at AADR scale (~1.2 GB packed either way,
//! well under modern RAM).

use super::{GenoReader, GenoWriter, Layout};
use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub const MIN_RLEN: usize = 48;
/// TGENO header is a fixed 48 bytes regardless of rlen — different from
/// PACKEDANCESTRYMAP, where the header is padded to rlen.
pub const HEADER_BYTES: usize = 48;
const HEADER_MAGIC: &[u8] = b"TGENO ";

#[inline]
pub fn rlen_for(nsnp: usize) -> usize {
    std::cmp::max(MIN_RLEN, (nsnp * 2 + 7) / 8)
}

// ======================================================================
// Reader
// ======================================================================

pub struct TgenoReader {
    mmap: Mmap,
    nind: usize,
    nsnp: usize,
    rlen: usize,
    record_bytes: usize,
    next_idx: usize,
    pub ihash: u32,
    pub shash: u32,
}

impl TgenoReader {
    pub fn open(path: &Path, nind: usize, nsnp: usize) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mmap =
            unsafe { Mmap::map(&file) }.with_context(|| format!("mmap {}", path.display()))?;

        let rlen = rlen_for(nsnp);
        let expected = HEADER_BYTES + rlen * nind;
        if mmap.len() < expected {
            bail!(
                "TGENO {}: file size {} < expected {} (header=48, rlen={}, nind={}, nsnp={})",
                path.display(),
                mmap.len(),
                expected,
                rlen,
                nind,
                nsnp
            );
        }
        if mmap.len() > expected {
            log::warn!(
                "TGENO {}: {} trailing bytes past expected end",
                path.display(),
                mmap.len() - expected
            );
        }

        let (h_nind, h_nsnp, ihash, shash) = parse_header(&mmap[..HEADER_BYTES])
            .with_context(|| format!("parse header of {}", path.display()))?;

        if h_nind != nind {
            bail!("TGENO header nind {h_nind} != .ind count {nind}");
        }
        if h_nsnp != nsnp {
            bail!("TGENO header nsnp {h_nsnp} != .snp count {nsnp}");
        }

        let record_bytes = (nsnp * 2 + 7) / 8;

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

impl GenoReader for TgenoReader {
    fn nind(&self) -> usize {
        self.nind
    }
    fn nsnp(&self) -> usize {
        self.nsnp
    }
    fn layout(&self) -> Layout {
        Layout::SampleMajor
    }
    fn record_bytes(&self) -> usize {
        self.record_bytes
    }
    fn header_hashes(&self) -> Option<(u32, u32)> {
        Some((self.ihash, self.shash))
    }

    fn read_record(&mut self, dst: &mut [u8]) -> Result<bool> {
        if self.next_idx >= self.nind {
            return Ok(false);
        }
        if dst.len() != self.record_bytes {
            bail!(
                "dst len {} != record_bytes {}",
                dst.len(),
                self.record_bytes
            );
        }
        let start = HEADER_BYTES + self.rlen * self.next_idx;
        let end = start + self.record_bytes;
        dst.copy_from_slice(&self.mmap[start..end]);
        self.next_idx += 1;
        Ok(true)
    }
}

fn parse_header(hdr: &[u8]) -> Result<(usize, usize, u32, u32)> {
    if !hdr.starts_with(HEADER_MAGIC) {
        bail!("not a TGENO file (missing 'TGENO ' magic)");
    }
    let end = hdr.iter().position(|&b| b == 0).unwrap_or(hdr.len());
    let text = std::str::from_utf8(&hdr[..end]).map_err(|e| anyhow!("header not UTF-8: {e}"))?;
    let mut it = text.split_ascii_whitespace();
    let _tgeno = it.next();
    let nind: usize = it
        .next()
        .ok_or_else(|| anyhow!("header missing nind"))?
        .parse()
        .map_err(|e| anyhow!("header nind: {e}"))?;
    let nsnp: usize = it
        .next()
        .ok_or_else(|| anyhow!("header missing nsnp"))?
        .parse()
        .map_err(|e| anyhow!("header nsnp: {e}"))?;
    let ihash: u32 = u32::from_str_radix(
        it.next().ok_or_else(|| anyhow!("header missing ihash"))?,
        16,
    )
    .map_err(|e| anyhow!("header ihash: {e}"))?;
    let shash: u32 = u32::from_str_radix(
        it.next().ok_or_else(|| anyhow!("header missing shash"))?,
        16,
    )
    .map_err(|e| anyhow!("header shash: {e}"))?;
    Ok((nind, nsnp, ihash, shash))
}

// ======================================================================
// Writer
// ======================================================================

pub struct TgenoWriter {
    w: BufWriter<File>,
    nind: usize,
    nsnp: usize,
    rlen: usize,
    record_bytes: usize,
    records_written: usize,
    header_written: bool,
}

impl TgenoWriter {
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

impl GenoWriter for TgenoWriter {
    fn layout(&self) -> Layout {
        Layout::SampleMajor
    }

    fn begin(&mut self, nind: usize, nsnp: usize, ihash: u32, shash: u32) -> Result<()> {
        if self.header_written {
            bail!("TgenoWriter::begin called twice");
        }
        self.nind = nind;
        self.nsnp = nsnp;
        self.rlen = rlen_for(nsnp);
        self.record_bytes = (nsnp * 2 + 7) / 8;

        let mut buf = vec![0u8; HEADER_BYTES];
        let header = format!("TGENO {} {} {:x} {:x}", nind, nsnp, ihash, shash);
        if header.len() > HEADER_BYTES {
            bail!("TGENO header {header:?} exceeds {HEADER_BYTES} bytes");
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
        if self.records_written >= self.nind {
            bail!("too many records written: expected {} samples", self.nind);
        }
        self.w.write_all(src)?;
        if self.rlen > self.record_bytes {
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
        if self.records_written != self.nind {
            bail!(
                "finish: wrote {} records, expected {}",
                self.records_written,
                self.nind
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
    fn rlen_min_48() {
        assert_eq!(rlen_for(1), 48);
        assert_eq!(rlen_for(192), 48);
        assert_eq!(rlen_for(193), 49);
    }

    #[test]
    fn parses_header() {
        let mut buf = vec![0u8; 48];
        let s = b"TGENO 5 7 cafe face";
        buf[..s.len()].copy_from_slice(s);
        let (nind, nsnp, ih, sh) = parse_header(&buf).unwrap();
        assert_eq!((nind, nsnp, ih, sh), (5, 7, 0xcafe, 0xface));
    }

    #[test]
    fn rejects_wrong_magic() {
        let buf = vec![b'Z'; 48];
        assert!(parse_header(&buf).is_err());
    }

    #[test]
    fn writer_reader_roundtrip() {
        // Small nind (3 samples) → rec_bytes small but rlen=48 pads it.
        let nind = 3;
        let nsnp = 17; // partial last byte
        let records: Vec<Vec<u8>> = (0..nind)
            .map(|i| (0..nsnp).map(|j| [0u8, 1, 2, 9][(i + j) % 4]).collect())
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = TgenoWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for rec in &records {
            let mut p = vec![0u8; (nsnp * 2 + 7) / 8];
            codec::pack(rec, &mut p);
            w.write_record(&p).unwrap();
        }
        w.finish().unwrap();

        let mut r = TgenoReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        let mut unpacked = vec![0u8; nsnp];
        for expected in &records {
            assert!(r.read_record(&mut buf).unwrap());
            codec::unpack(&buf, nsnp, &mut unpacked);
            assert_eq!(&unpacked, expected);
        }
        assert!(!r.read_record(&mut buf).unwrap());
    }

    #[test]
    fn writer_reader_roundtrip_large_nsnp() {
        let nind = 2;
        let nsnp = 400; // rec_bytes = 100, > MIN_RLEN
        let records: Vec<Vec<u8>> = (0..nind)
            .map(|i| (0..nsnp).map(|j| [0u8, 1, 2, 9][(i * 3 + j) % 4]).collect())
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let mut w = TgenoWriter::create(tmp.path()).unwrap();
        w.begin(nind, nsnp, 0, 0).unwrap();
        for rec in &records {
            let mut p = vec![0u8; (nsnp * 2 + 7) / 8];
            codec::pack(rec, &mut p);
            w.write_record(&p).unwrap();
        }
        w.finish().unwrap();

        let mut r = TgenoReader::open(tmp.path(), nind, nsnp).unwrap();
        let mut buf = vec![0u8; r.record_bytes()];
        let mut unpacked = vec![0u8; nsnp];
        for expected in &records {
            assert!(r.read_record(&mut buf).unwrap());
            codec::unpack(&buf, nsnp, &mut unpacked);
            assert_eq!(&unpacked, expected);
        }
    }
}
