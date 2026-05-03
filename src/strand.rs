//! Allele flip and strand-check logic.
//!
//! Used when conversion touches PLINK formats or when the user supplies
//! `allelename` (a reference set of (snp_id, a1, a2)). Upstream convertf
//! behavior matrix:
//!
//! - `flipreference: YES` → for each SNP present in `allelename` with
//!   swapped alleles vs input, flip genotype 0 ↔ 2 and swap a1/a2.
//! - `flipstrand: YES` → also accept complement-swap (A↔T, C↔G).
//! - `strandcheck: YES` (default) → drop A/T and C/G SNPs because strand
//!   cannot be resolved.
//!
//! # Complement table
//!
//! ```text
//! A ↔ T
//! C ↔ G
//! ```

/// Complement of an ASCII allele byte. Non-ACGT returns self.
#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

/// True if {a1,a2} is one of the ambiguous strand pairs (A/T or C/G).
#[inline]
pub fn is_ambiguous(a1: u8, a2: u8) -> bool {
    matches!(
        (a1.to_ascii_uppercase(), a2.to_ascii_uppercase()),
        (b'A', b'T') | (b'T', b'A') | (b'C', b'G') | (b'G', b'C')
    )
}

/// Decide whether genotypes need 0↔2 flip given input and reference alleles.
///
/// Returns:
/// - `Some(false)`: alleles match, no flip.
/// - `Some(true)`: alleles are swapped (possibly after strand complement), flip.
/// - `None`: incompatible — SNP should be dropped.
pub fn decide_flip(in_a1: u8, in_a2: u8, ref_a1: u8, ref_a2: u8, flipstrand: bool) -> Option<bool> {
    let (a, b) = (in_a1.to_ascii_uppercase(), in_a2.to_ascii_uppercase());
    let (x, y) = (ref_a1.to_ascii_uppercase(), ref_a2.to_ascii_uppercase());
    if (a, b) == (x, y) {
        return Some(false);
    }
    if (a, b) == (y, x) {
        return Some(true);
    }
    if flipstrand {
        let (ac, bc) = (complement(a), complement(b));
        if (ac, bc) == (x, y) {
            return Some(false);
        }
        if (ac, bc) == (y, x) {
            return Some(true);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_basic() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn ambiguous_pairs() {
        assert!(is_ambiguous(b'A', b'T'));
        assert!(is_ambiguous(b'C', b'G'));
        assert!(!is_ambiguous(b'A', b'C'));
    }

    #[test]
    fn flip_logic() {
        // Same alleles, no flip.
        assert_eq!(decide_flip(b'A', b'C', b'A', b'C', false), Some(false));
        // Swapped.
        assert_eq!(decide_flip(b'C', b'A', b'A', b'C', false), Some(true));
        // Strand flip, no swap.
        assert_eq!(decide_flip(b'T', b'G', b'A', b'C', true), Some(false));
        // Strand flip + allele swap.
        assert_eq!(decide_flip(b'G', b'T', b'A', b'C', true), Some(true));
        // Strand flip disabled → incompatible.
        assert_eq!(decide_flip(b'T', b'G', b'A', b'C', false), None);
    }
}
