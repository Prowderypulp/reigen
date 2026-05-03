//! Metadata types shared across all formats.
//!
//! `.snp` / `.bim` describe SNPs; `.ind` / `.fam` describe samples.
//! Column conventions differ between AdmixTools and PLINK — see `snp.rs`,
//! `bim.rs`, `ind.rs`, `fam.rs` for exact layouts.

pub mod bim;
pub mod fam;
pub mod ind;
pub mod snp;

/// Split mmap bytes at `\n` without allocating. Strips trailing `\r`.
/// Shared by all line-oriented metadata parsers.
pub(crate) fn split_lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut start = 0usize;
    let len = bytes.len();
    std::iter::from_fn(move || {
        if start >= len {
            return None;
        }
        let rest = &bytes[start..];
        match memchr::memchr(b'\n', rest) {
            Some(off) => {
                let mut end = start + off;
                let line_start = start;
                start = end + 1;
                if end > line_start && bytes[end - 1] == b'\r' {
                    end -= 1;
                }
                Some(&bytes[line_start..end])
            }
            None => {
                let line_start = start;
                start = len;
                Some(&bytes[line_start..len])
            }
        }
    })
}

/// One SNP row. Normalized across `.snp` and `.bim`.
///
/// Coordinate convention: physical position is 1-based (VCF-style), matching
/// upstream AdmixTools on-disk representation.
#[derive(Debug, Clone)]
pub struct SnpRow {
    pub id: String,
    pub chrom: u8,         // 1..=numchrom, numchrom+1 = X, +2 = Y, +3 = MT
    pub genetic_pos: f64,  // Morgans
    pub physical_pos: u64, // 1-based
    pub allele1: u8,       // ASCII 'A'/'C'/'G'/'T' — reference in EIGENSTRAT convention
    pub allele2: u8,       // ASCII — "variant" allele
}

/// One sample row. Normalized across `.ind` and `.fam`.
#[derive(Debug, Clone)]
pub struct IndRow {
    pub id: String,
    pub sex: Sex,
    /// Population / family label. In `.ind` this is column 3; in `.fam` it's column 1 (FID).
    pub pop: String,
    /// `true` if `.ind` column 3 was literally "Ignore" (case-insensitive).
    pub ignore: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

impl Sex {
    pub fn from_char(c: u8) -> Self {
        match c {
            b'M' | b'm' | b'1' => Sex::Male,
            b'F' | b'f' | b'2' => Sex::Female,
            _ => Sex::Unknown,
        }
    }

    pub fn as_ind_char(self) -> char {
        match self {
            Sex::Male => 'M',
            Sex::Female => 'F',
            Sex::Unknown => 'U',
        }
    }

    pub fn as_fam_code(self) -> u8 {
        match self {
            Sex::Male => b'1',
            Sex::Female => b'2',
            Sex::Unknown => b'0',
        }
    }
}
