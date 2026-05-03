//! Patterson's array hash (`hasharr` in `admutils.c`).
//!
//! Used by PACKEDANCESTRYMAP / TGENO to detect stale geno files when sample
//! IDs or SNP IDs change. Stored as hex in the header:
//!
//! ```text
//! GENO nind nsnp ihash shash
//! ```
//!
//! Where `ihash = hasharr(sample_ids)` and `shash = hasharr(snp_ids)`.
//!
//! # C reference (Patterson, admutils.c)
//!
//! ```c
//! int hashit(char *str) {
//!     int j, len, hash = 0;
//!     len = strlen(str);
//!     for (j=0; j<len; j++) {
//!         hash *= 23;
//!         hash += (int) str[j];
//!     }
//!     return hash;
//! }
//!
//! int hasharr(char **xarr, int nxarr) {
//!     int hash = 0, thash, i;
//!     for (i=0; i<nxarr; i++) {
//!         thash = hashit(xarr[i]);
//!         hash *= 17;
//!         hash ^= thash;
//!     }
//!     return hash;
//! }
//! ```
//!
//! # Type / overflow semantics
//!
//! - `int` is `i32` on every platform AdmixTools targets (Linux x86_64).
//!   Convertf has been built and run that way for 15 years; this is the
//!   only correct interpretation.
//! - C signed integer overflow is UB, but in practice GCC/Clang at -O2
//!   emit two's-complement wraparound and Patterson relies on it.
//!   We replicate with `i32::wrapping_mul` / `wrapping_add`.
//! - `(int) str[j]` sign-extends signed-`char` bytes. Linux x86_64 GCC
//!   default: `char` is signed. Sample/SNP IDs are ASCII (`< 128`), so
//!   sign-extension is a no-op for real data. We match the C semantics
//!   anyway in case anyone has high-bit chars.
//! - Final return is `int` → printed via `%x` as the bottom 32 bits.
//!   We return `u32` (the bit pattern) for storage in the geno header.

/// Patterson's per-string hash (`hashit`).
pub fn hashit(s: &str) -> i32 {
    let mut h: i32 = 0;
    for &b in s.as_bytes() {
        // C: `(int) str[j]` where char is signed on Linux x86_64.
        // Reinterpret byte as i8 then sign-extend to i32.
        let c = (b as i8) as i32;
        h = h.wrapping_mul(23).wrapping_add(c);
    }
    h
}

/// Patterson's order-dependent array hash (`hasharr`).
///
/// Returns the C `int` bit pattern reinterpreted as `u32`, ready for the
/// `"GENO ... %x %x"` header.
pub fn hasharr(ids: &[&str]) -> u32 {
    let mut hash: i32 = 0;
    for s in ids {
        let thash = hashit(s);
        hash = hash.wrapping_mul(17);
        hash ^= thash;
    }
    hash as u32
}

/// Convenience overload for owned strings (avoids `.iter().map(String::as_str).collect()`
/// at call sites).
#[allow(dead_code)]
pub fn hasharr_owned(ids: &[String]) -> u32 {
    let mut hash: i32 = 0;
    for s in ids {
        let thash = hashit(s);
        hash = hash.wrapping_mul(17);
        hash ^= thash;
    }
    hash as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hashit_empty() {
        assert_eq!(hashit(""), 0);
    }

    #[test]
    fn hashit_single_char() {
        // 0 * 23 + 'A'(65) = 65
        assert_eq!(hashit("A"), 65);
    }

    #[test]
    fn hashit_two_chars() {
        // 0 → *23 +'A'(65) = 65 → *23 +'B'(66) = 65*23 + 66 = 1561
        assert_eq!(hashit("AB"), 1561);
    }

    #[test]
    fn hashit_known_strings() {
        // Hand-traced reference values; if these break the implementation
        // diverged from C (or these expected values are wrong — recompute
        // from the C function above).
        // "ABC" = ((65*23 + 66) * 23 + 67) = 1561*23 + 67 = 35970
        assert_eq!(hashit("ABC"), 35970);
        // Long string to exercise wrapping at i32 overflow.
        // Just confirm it doesn't panic and is deterministic:
        let h1 = hashit("abcdefghijklmnopqrstuvwxyz0123456789");
        let h2 = hashit("abcdefghijklmnopqrstuvwxyz0123456789");
        assert_eq!(h1, h2);
    }

    #[test]
    fn hasharr_empty_zero() {
        let v: Vec<&str> = vec![];
        assert_eq!(hasharr(&v), 0);
    }

    #[test]
    fn hasharr_single_matches_hashit_cast() {
        let v = vec!["A"];
        // hash = 0; thash=hashit("A")=65; hash = 0*17 = 0; hash ^= 65 = 65
        assert_eq!(hasharr(&v), 65u32);
    }

    #[test]
    fn hasharr_order_sensitive() {
        let a = hasharr(&["A", "B"]);
        let b = hasharr(&["B", "A"]);
        assert_ne!(a, b, "order must affect hash");
    }

    #[test]
    fn hasharr_two_strings_known() {
        // hash=0
        // i=0: thash=hashit("A")=65; hash = 0*17 ^ 65 = 65
        // i=1: thash=hashit("B")=66; hash = 65*17 ^ 66 = 1105 ^ 66 = 1043
        assert_eq!(hasharr(&["A", "B"]) as i32, 1043);
    }

    #[test]
    fn hasharr_owned_matches_str_version() {
        let owned: Vec<String> = vec!["foo".into(), "bar".into(), "baz".into()];
        let refs: Vec<&str> = owned.iter().map(String::as_str).collect();
        assert_eq!(hasharr_owned(&owned), hasharr(&refs));
    }
}
