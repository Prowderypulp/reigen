//! PGEN writer — implements [`GenoWriter`] for `.pgen`.
//!
//! Emits **storage mode 0x02** (fixed-width, all type-0 uncompressed records),
//! the simplest spec-valid PGEN. Writable in a single sequential pass with no
//! header backpatching.
//!
//! # Encoding
//!
//! Input: canonical AM MSB-first packed bytes ({0,1,2,9} → 00,01,10,11).
//! Output: PGEN LSB-first packed bytes ({0,1,2,3} → 00,01,10,11).
//! The 2-bit values are the same (identity map), only the bit order changes.
//! Uses [`AM_TO_PGEN`] for full-byte conversion.

use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::geno::{GenoWriter, Layout};

use super::record::AM_TO_PGEN;

pub const PGEN_FIXED_MAGIC: [u8; 3] = [0x6c, 0x1b, 0x02];

pub const HEADER_CTRL_ALL_PROVISIONAL: u8 = 0x80;

pub struct PgenWriter {
    w: BufWriter<File>,
    nind: usize,
    nsnp: usize,
    rec_bytes: usize,
    pgen_rec_bytes: usize,
    records_written: usize,
    began: bool,
    out_buf: Vec<u8>,
}

impl PgenWriter {
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
        Ok(Self {
            w: BufWriter::with_capacity(256 * 1024, file),
            nind: 0,
            nsnp: 0,
            rec_bytes: 0,
            pgen_rec_bytes: 0,
            records_written: 0,
            began: false,
            out_buf: Vec::new(),
        })
    }
}

impl GenoWriter for PgenWriter {
    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }

    fn begin(&mut self, nind: usize, nsnp: usize, _ihash: u32, _shash: u32) -> Result<()> {
        if self.began {
            bail!("PgenWriter::begin called twice");
        }
        self.nind = nind;
        self.nsnp = nsnp;
        self.rec_bytes = (nind * 2 + 7) / 8;
        self.pgen_rec_bytes = nind.div_ceil(4);
        self.out_buf = vec![0u8; self.pgen_rec_bytes];

        // Header: magic (2) + mode (1) + M (4) + N (4) + header_ctrl (1) = 12 bytes.
        // M = variant count, N = sample count, LE u32. Checked casts: a silent
        // `as u32` truncation would corrupt the dims for >4G variants/samples
        // (B-103). PGEN's on-disk dims are u32, so anything larger is unwritable.
        let m = u32::try_from(nsnp)
            .map_err(|_| anyhow::anyhow!("PGEN variant count {nsnp} exceeds u32::MAX"))?;
        let n = u32::try_from(nind)
            .map_err(|_| anyhow::anyhow!("PGEN sample count {nind} exceeds u32::MAX"))?;
        self.w.write_all(&PGEN_FIXED_MAGIC)?;
        self.w.write_all(&m.to_le_bytes())?;
        self.w.write_all(&n.to_le_bytes())?;
        self.w.write_all(&[HEADER_CTRL_ALL_PROVISIONAL])?;

        self.began = true;
        Ok(())
    }

    fn write_record(&mut self, src: &[u8]) -> Result<()> {
        if !self.began {
            bail!("write_record before begin");
        }
        if src.len() != self.rec_bytes {
            bail!(
                "record len {} != expected {}",
                src.len(),
                self.rec_bytes
            );
        }
        if self.records_written >= self.nsnp {
            bail!("too many records: expected {}", self.nsnp);
        }

        // AM_TO_PGEN LUT: canonical MSB-first → PGEN LSB-first, per-byte.
        for (i, &b) in src.iter().enumerate().take(self.pgen_rec_bytes) {
            self.out_buf[i] = AM_TO_PGEN[b as usize];
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geno::codec;
    use crate::geno::GenoReader;

    #[test]
    fn writer_reader_roundtrip() {
        // 11 samples (partial tail), 20 variants.
        let nind: usize = 11;
        let nsnp: usize = 20;
        let records: Vec<Vec<u8>> = (0..nsnp)
            .map(|snp| {
                (0..nind)
                    .map(|i| match (snp + i * 3) % 4 {
                        0 => 0u8,
                        1 => 1u8,
                        2 => 2u8,
                        _ => 9u8,
                    })
                    .collect()
            })
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        {
            let mut w = PgenWriter::create(tmp.path()).unwrap();
            w.begin(nind, nsnp, 0, 0).unwrap();
            for rec in &records {
                let mut packed = vec![0u8; (nind * 2 + 7) / 8];
                codec::pack(rec, &mut packed);
                w.write_record(&packed).unwrap();
            }
            w.finish().unwrap();
        }

        // Read back with PgenReader and compare.
        let mut r = crate::geno::pgen::reader::PgenReader::open(tmp.path(), vec![]).unwrap();
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
    fn written_pgen_is_valid_mode02() {
        let nind = 4;
        let nsnp = 2;
        let tmp = tempfile::NamedTempFile::new().unwrap();
        {
            let mut w = PgenWriter::create(tmp.path()).unwrap();
            w.begin(nind, nsnp, 0, 0).unwrap();
            let mut p = vec![0u8; (nind * 2 + 7) / 8];
            codec::pack(&[2, 1, 0, 9], &mut p);
            w.write_record(&p).unwrap();
            codec::pack(&[0, 0, 2, 2], &mut p);
            w.write_record(&p).unwrap();
            w.finish().unwrap();
        }

        let data = std::fs::read(tmp.path()).unwrap();
        // Header: magic + mode + M=2 + N=4 + hc=0x80 = 12 bytes
        assert_eq!(&data[0..3], &PGEN_FIXED_MAGIC);
        assert_eq!(u32::from_le_bytes([data[3], data[4], data[5], data[6]]), 2);
        assert_eq!(u32::from_le_bytes([data[7], data[8], data[9], data[10]]), 4);
        assert_eq!(data[11], 0x80);
        // 2 records × ceil(4/4) = 2 bytes
        assert_eq!(data.len(), 14);
        assert_eq!(data[12] & 0b11, 2, "sample 0 should be PGEN 2");
        assert_eq!((data[12] >> 2) & 0b11, 1, "sample 1 should be PGEN 1");
    }
}
