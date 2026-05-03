//! ANCESTRYMAP format (sparse text).
//!
//! One line per *non-missing* genotype:
//!
//! ```text
//! snp_id  sample_id  genotype
//! ```
//!
//! Implied-missing: any (snp, sample) pair not listed is missing.
//!
//! # v1 scope
//!
//! Low priority. Input support matters (legacy datasets ship in this form);
//! output rarely requested. We stub both; implement reader first.

// Stubs for the sparse ANCESTRYMAP text format. Phase-G item; not wired
// into the pipeline dispatch. Kept to reserve the module and API shape.
#![allow(dead_code)]

use super::{GenoReader, GenoWriter, Layout};
use anyhow::Result;
use std::path::Path;

pub struct AncestrymapReader {
    _nind: usize,
    _nsnp: usize,
}

impl AncestrymapReader {
    pub fn open(_path: &Path, _nind: usize, _nsnp: usize) -> Result<Self> {
        todo!("two-pass: build {{(snp,sample) -> g}} hashmap, then emit records")
    }
}

impl GenoReader for AncestrymapReader {
    fn nind(&self) -> usize {
        self._nind
    }
    fn nsnp(&self) -> usize {
        self._nsnp
    }
    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }
    fn read_record(&mut self, _dst: &mut [u8]) -> Result<bool> {
        todo!()
    }
}

pub struct AncestrymapWriter {
    _nind: usize,
    _nsnp: usize,
}

impl AncestrymapWriter {
    pub fn create(_path: &Path) -> Result<Self> {
        todo!()
    }
}

impl GenoWriter for AncestrymapWriter {
    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }
    fn begin(&mut self, _n: usize, _m: usize, _ih: u32, _sh: u32) -> Result<()> {
        todo!()
    }
    fn write_record(&mut self, _src: &[u8]) -> Result<()> {
        todo!()
    }
    fn finish(&mut self) -> Result<()> {
        todo!()
    }
}
