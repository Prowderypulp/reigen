//! Bit-level genotype packing/unpacking primitives.
//!
//! # Conventions
//!
//! Canonical 2-bit encoding: `00=0, 01=1, 10=2, 11=missing`, MSB-first within byte.
//!
//! PLINK `.bed` uses a different convention handled in `packed_ped.rs`.
//!
//! # Scalar today, SIMD later
//!
//! Current implementations are scalar. Hot paths that warrant SIMD:
//! - `ascii_to_packed` / `packed_to_ascii` (EIGENSTRAT text I/O) — AVX2 pshufb LUT.
//! - PLINK bit-order reverse + recode — AVX2 shuffle.
//! - 8×8 bit-matrix transpose (transpose.rs uses this).
//!
//! Keep function signatures stable so SIMD impls can be swapped in behind
//! `#[cfg(target_feature = ...)]` without touching callers.

/// Canonical missing value in unpacked (byte-per-genotype) form.
pub const G_MISSING: u8 = 9;

/// Pack a slice of genotypes (0/1/2/9, one per byte) into canonical 2-bit MSB-first.
///
/// `dst.len()` must be `ceil(src.len() * 2 / 8)`. Trailing bits of the last byte
/// are zero-padded.
pub fn pack(src: &[u8], dst: &mut [u8]) {
    let need = (src.len() * 2 + 7) / 8;
    debug_assert!(dst.len() >= need);
    for b in dst[..need].iter_mut() {
        *b = 0;
    }

    for (i, &g) in src.iter().enumerate() {
        let two = g_to_2bit(g);
        let byte = i / 4;
        let shift = 6 - 2 * (i % 4); // MSB-first: sample 0 at bits 7-6
        dst[byte] |= two << shift;
    }
}

/// Unpack `n` genotypes from canonical 2-bit MSB-first into one byte each.
pub fn unpack(src: &[u8], n: usize, dst: &mut [u8]) {
    debug_assert!(dst.len() >= n);
    debug_assert!(src.len() >= (n * 2 + 7) / 8);
    for i in 0..n {
        let byte = i / 4;
        let shift = 6 - 2 * (i % 4);
        let two = (src[byte] >> shift) & 0b11;
        dst[i] = two_bit_to_g(two);
    }
}

/// 0/1/2/9 → 2-bit (any non-{0,1,2} → 3 = missing).
#[inline(always)]
pub fn g_to_2bit(g: u8) -> u8 {
    match g {
        0 => 0b00,
        1 => 0b01,
        2 => 0b10,
        _ => 0b11,
    }
}

/// 2-bit → 0/1/2/9.
#[inline(always)]
pub fn two_bit_to_g(t: u8) -> u8 {
    match t & 0b11 {
        0 => 0,
        1 => 1,
        2 => 2,
        _ => G_MISSING,
    }
}

/// EIGENSTRAT ASCII byte ('0'/'1'/'2'/'9') → 0/1/2/9.
#[inline(always)]
#[allow(dead_code)] // public codec helper, retained for the EIGENSTRAT text path
pub fn ascii_to_g(c: u8) -> u8 {
    match c {
        b'0' => 0,
        b'1' => 1,
        b'2' => 2,
        _ => G_MISSING,
    }
}

/// 0/1/2/9 → EIGENSTRAT ASCII byte.
#[inline(always)]
#[allow(dead_code)] // public codec helper, retained for the EIGENSTRAT text path
pub fn g_to_ascii(g: u8) -> u8 {
    match g {
        0 => b'0',
        1 => b'1',
        2 => b'2',
        _ => b'9',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_scalar() {
        let src: Vec<u8> = (0..100).map(|i| [0u8, 1, 2, 9][i % 4]).collect();
        let mut packed = vec![0u8; (src.len() * 2 + 7) / 8];
        pack(&src, &mut packed);
        let mut back = vec![0u8; src.len()];
        unpack(&packed, src.len(), &mut back);
        assert_eq!(src, back);
    }

    #[test]
    fn msb_first_layout() {
        // 4 genotypes 0, 1, 2, 9 (missing) → bits 00 01 10 11 → 0x1B, MSB-first.
        let src = [0u8, 1, 2, 9];
        let mut packed = [0u8];
        pack(&src, &mut packed);
        assert_eq!(packed[0], 0b00_01_10_11);
    }

    #[test]
    fn ascii_roundtrip() {
        for g in [0u8, 1, 2, 9] {
            assert_eq!(ascii_to_g(g_to_ascii(g)), g);
        }
        // Garbage char → missing.
        assert_eq!(ascii_to_g(b'X'), G_MISSING);
    }
}
