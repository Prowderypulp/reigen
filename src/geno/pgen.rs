//! PLINK 2 PGEN genotype codec.
//!
//! PGEN (`.pgen`) is PLINK 2's genotype container. Unlike PLINK 1 `.bed`
//! (handled by [`super::packed_ped`]), PGEN supports compressed variable-width
//! variant records and a richer set of auxiliary tracks (phase, dosage,
//! multiallelic). reigen targets the **biallelic unphased hardcall** subset:
//! the always-present main 2-bit track is decoded and everything else is
//! skipped (see `PGEN_PLAN.md` §1).
//!
//! # Polarity (the load-bearing invariant — `PGEN_PLAN.md` §2)
//!
//! PGEN main-track genocodes `{0,1,2,3}` = `{hom-REF, het, double-ALT,
//! missing}`. With the PLINK-1-compatible PVAR header (`ALT` in col 5 = reigen
//! `allele1`), the map to reigen canonical `{0,1,2,9}` is the **identity**.
//! This is **not** the `.bed` mapping — do not reuse `packed_ped`'s LUT.
//!
//! # Bit order
//!
//! The main 2-bit array is **LSB-first** (sample 0 in the low bits), confirmed
//! empirically in `tests/golden/plink2/T0_FINDINGS.md`. (PACKEDANCESTRYMAP is
//! MSB-first; the codec converts.)
//!
//! # Module layout
//!
//! - [`difflist`]: base-128 varints and the difflist decoder, shared by record
//!   types 1, 2, 3, 4, 6, 7.
//! - [`header`]: header parsing (modes, dims, byte 11, block offsets, per-record
//!   type/length arrays and spans).
//! - [`record`]: record decoders for all 7 hardcall record types (0–7), plus the
//!   PGEN↔AM LUTs (bit-order reverse, no semantic recode).

pub mod difflist;
pub mod header;
pub mod reader;
pub mod record;
pub mod writer;

pub use difflist::{decode_difflist, ByteCursor, DiffList, SampleIdWidth};
pub use header::{HeaderFormat, PgenHeader, ProvRefMode, StorageMode};
pub use reader::PgenReader;
pub use record::{PGEN_TO_AM, AM_TO_PGEN, decode_record};
pub use writer::PgenWriter;
