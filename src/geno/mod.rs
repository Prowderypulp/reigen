//! Genotype I/O.
//!
//! # Canonical in-memory encoding
//!
//! All readers normalize to, and all writers accept, **PACKEDANCESTRYMAP
//! 2-bit convention**:
//!
//! | 2 bits | meaning                |
//! |--------|------------------------|
//! | `00`   | genotype 0 (hom ref)   |
//! | `01`   | genotype 1 (het)       |
//! | `10`   | genotype 2 (hom alt)   |
//! | `11`   | missing                |
//!
//! Within a byte, the **MSB** holds the first sample. That matches
//! PACKEDANCESTRYMAP and TGENO on-disk layout, and differs from PLINK `.bed`
//! (LSB-first, different semantics — handled inside `packed_ped.rs`).
//!
//! # Record definition
//!
//! - SNP-major formats (EIGENSTRAT, PACKEDANCESTRYMAP, PACKEDPED): one record
//!   = one SNP × nind samples.
//! - Sample-major formats (TGENO): one record = one sample × nsnp SNPs.
//!
//! A conversion between layouts goes through `transpose.rs`.

pub mod ancestrymap;
pub mod codec;
pub mod eigenstrat;
pub mod packed_am;
pub mod packed_ped;
pub mod tgeno;

use anyhow::Result;

/// Record orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    /// One record per SNP, covering all samples.
    SnpMajor,
    /// One record per sample, covering all SNPs.
    SampleMajor,
}

/// Streaming reader. Produces records in canonical 2-bit encoding.
pub trait GenoReader {
    fn nind(&self) -> usize;
    fn nsnp(&self) -> usize;
    fn layout(&self) -> Layout;

    /// Number of bytes per record in canonical encoding:
    /// `ceil(n * 2 / 8)` where `n` is nind (SnpMajor) or nsnp (SampleMajor).
    /// No 48-byte minimum applies to the in-memory buffer — only on disk.
    fn record_bytes(&self) -> usize {
        let n = match self.layout() {
            Layout::SnpMajor => self.nind(),
            Layout::SampleMajor => self.nsnp(),
        };
        (n * 2 + 7) / 8
    }

    /// Hashes embedded in the geno file header, when the format has them
    /// (PACKEDANCESTRYMAP, TGENO). `None` for header-less formats
    /// (EIGENSTRAT, PACKEDPED). Used by the pipeline to enforce
    /// `hashcheck: YES`.
    fn header_hashes(&self) -> Option<(u32, u32)> {
        None
    }

    /// Read next record into `dst`. `dst.len()` must equal `record_bytes()`.
    /// Returns Ok(false) at EOF.
    fn read_record(&mut self, dst: &mut [u8]) -> Result<bool>;
}

/// Streaming writer. Accepts records in canonical 2-bit encoding.
pub trait GenoWriter {
    /// Which layout this writer consumes — SnpMajor or SampleMajor.
    /// Must match the reader's layout (or go through transpose).
    fn layout(&self) -> Layout;

    /// Called once before any records. ihash/shash may be 0 when
    /// `hashcheck: NO` is acceptable downstream.
    fn begin(&mut self, nind: usize, nsnp: usize, ihash: u32, shash: u32) -> Result<()>;

    fn write_record(&mut self, src: &[u8]) -> Result<()>;

    fn finish(&mut self) -> Result<()>;
}
