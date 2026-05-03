//! Random-access reader wrappers for merge.
//!
//! PAM/BED readers are mmap-backed; seeking is just setting the `next_idx`
//! cursor. Text (EIGENSTRAT) and sample-major (TGENO) inputs are rejected —
//! merge v1 requires SNP-major fixed-width formats.

use crate::format::Format;
use crate::geno::packed_am::PackedAmReader;
use crate::geno::packed_ped::PackedPedReader;
use crate::geno::GenoReader;
use anyhow::Result;
use std::path::Path;

pub trait SeekableGenoReader: GenoReader {
    fn seek_record(&mut self, idx: usize) -> Result<()>;
}

impl SeekableGenoReader for PackedAmReader {
    fn seek_record(&mut self, idx: usize) -> Result<()> {
        self.set_next_idx(idx);
        Ok(())
    }
}

impl SeekableGenoReader for PackedPedReader {
    fn seek_record(&mut self, idx: usize) -> Result<()> {
        self.set_next_idx(idx);
        Ok(())
    }
}

pub fn open_seekable(
    fmt: Format,
    path: &Path,
    nind: usize,
    nsnp: usize,
) -> Result<Box<dyn SeekableGenoReader>> {
    match fmt {
        Format::PackedAncestrymap => Ok(Box::new(PackedAmReader::open(path, nind, nsnp)?)),
        Format::PackedPed => Ok(Box::new(PackedPedReader::open(path, nind, nsnp)?)),
        Format::Eigenstrat => {
            anyhow::bail!("EIGENSTRAT text not yet supported for merge (no line index)")
        }
        Format::Tgeno => {
            anyhow::bail!("TGENO (sample-major) not supported for merge — convert to PAM first")
        }
        _ => anyhow::bail!("format {:?} not supported for merge", fmt),
    }
}
