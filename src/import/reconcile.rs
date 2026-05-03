//! Reconcile a DTC kit call against a reference SNP and produce the
//! canonical genotype code in one step.
//!
//! EIGENSTRAT encoding: 0 = hom ref (allele1), 1 = het, 2 = hom alt (allele2).
//! We count copies of the reference's allele2, after optionally
//! complementing kit alleles to bring them onto the reference strand.

use super::vendor::DtcCall;
use crate::meta::SnpRow;
use crate::strand::{complement, is_ambiguous};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReconcileResult {
    /// Successfully encoded. Value is 0/1/2.
    Encoded(u8),
    /// Reference SNP is A/T or C/G and `--allow-ambiguous` not set.
    DropAmbiguous,
    /// Kit alleles don't match reference (even under complement).
    DropMismatch,
    /// Kit reported a no-call (one or both alleles missing).
    DropNoCall,
    /// Strand-flip complement was required but `--flip-strand` is off.
    DropStrandUnresolvable,
}

pub fn reconcile_call(
    call: &DtcCall,
    ref_snp: &SnpRow,
    allow_ambiguous: bool,
    flip_strand: bool,
) -> ReconcileResult {
    let a = call.a1 as u8;
    let b = call.a2 as u8;
    if !is_valid_base(a) || !is_valid_base(b) {
        return ReconcileResult::DropNoCall;
    }
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();

    let x = ref_snp.allele1.to_ascii_uppercase(); // reference allele
    let y = ref_snp.allele2.to_ascii_uppercase(); // variant allele

    if !allow_ambiguous && is_ambiguous(x, y) {
        return ReconcileResult::DropAmbiguous;
    }

    // Try direct-strand match first.
    match count_alt_copies(a, b, x, y) {
        AlleleMatch::Count(n) => return ReconcileResult::Encoded(n),
        AlleleMatch::Incompatible => {}
    }

    // Then try reverse-complement if allowed.
    if flip_strand {
        let ac = complement(a);
        let bc = complement(b);
        match count_alt_copies(ac, bc, x, y) {
            AlleleMatch::Count(n) => return ReconcileResult::Encoded(n),
            AlleleMatch::Incompatible => return ReconcileResult::DropMismatch,
        }
    }

    // A direct match failed and we weren't allowed to try complement —
    // could be genuine mismatch or a strand issue we can't resolve.
    if complement_would_match(a, b, x, y) {
        ReconcileResult::DropStrandUnresolvable
    } else {
        ReconcileResult::DropMismatch
    }
}

enum AlleleMatch {
    Count(u8),
    Incompatible,
}

/// Returns the number of copies of `y` (variant allele) in kit call `(a, b)`.
/// `Incompatible` if either kit allele is neither `x` nor `y`.
#[inline]
fn count_alt_copies(a: u8, b: u8, x: u8, y: u8) -> AlleleMatch {
    let ca = if a == y {
        1
    } else if a == x {
        0
    } else {
        return AlleleMatch::Incompatible;
    };
    let cb = if b == y {
        1
    } else if b == x {
        0
    } else {
        return AlleleMatch::Incompatible;
    };
    AlleleMatch::Count(ca + cb)
}

#[inline]
fn complement_would_match(a: u8, b: u8, x: u8, y: u8) -> bool {
    matches!(
        count_alt_copies(complement(a), complement(b), x, y),
        AlleleMatch::Count(_)
    )
}

#[inline]
fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

#[cfg(test)]
mod tests {
    use super::*;

    fn snp(a1: u8, a2: u8) -> SnpRow {
        SnpRow {
            id: "rs1".into(),
            chrom: 1,
            genetic_pos: 0.0,
            physical_pos: 100,
            allele1: a1,
            allele2: a2,
        }
    }

    fn call(a1: char, a2: char) -> DtcCall {
        DtcCall {
            rsid: "rs1".into(),
            chrom: "1".into(),
            pos: 100,
            a1,
            a2,
        }
    }

    #[test]
    fn hom_ref_is_zero() {
        assert_eq!(
            reconcile_call(&call('A', 'A'), &snp(b'A', b'G'), false, false),
            ReconcileResult::Encoded(0)
        );
    }

    #[test]
    fn het_is_one() {
        assert_eq!(
            reconcile_call(&call('A', 'G'), &snp(b'A', b'G'), false, false),
            ReconcileResult::Encoded(1)
        );
        // Order-insensitive.
        assert_eq!(
            reconcile_call(&call('G', 'A'), &snp(b'A', b'G'), false, false),
            ReconcileResult::Encoded(1)
        );
    }

    #[test]
    fn hom_alt_is_two() {
        assert_eq!(
            reconcile_call(&call('G', 'G'), &snp(b'A', b'G'), false, false),
            ReconcileResult::Encoded(2)
        );
    }

    #[test]
    fn complement_het() {
        // Ref is A/G, kit reports T/C — same variant on opposite strand, het.
        assert_eq!(
            reconcile_call(&call('T', 'C'), &snp(b'A', b'G'), false, true),
            ReconcileResult::Encoded(1)
        );
    }

    #[test]
    fn complement_without_flag_drops() {
        assert_eq!(
            reconcile_call(&call('T', 'C'), &snp(b'A', b'G'), false, false),
            ReconcileResult::DropStrandUnresolvable
        );
    }

    #[test]
    fn ambiguous_drops() {
        assert_eq!(
            reconcile_call(&call('A', 'T'), &snp(b'A', b'T'), false, false),
            ReconcileResult::DropAmbiguous
        );
        // Allowed if flag set.
        assert_eq!(
            reconcile_call(&call('A', 'T'), &snp(b'A', b'T'), true, false),
            ReconcileResult::Encoded(1)
        );
    }

    #[test]
    fn mismatch_drops() {
        // Ref A/G, kit reports A/C — C is not in the ref pair and its
        // complement G is fine but mixing with A/C doesn't match either strand.
        assert_eq!(
            reconcile_call(&call('A', 'C'), &snp(b'A', b'G'), false, true),
            ReconcileResult::DropMismatch
        );
    }

    #[test]
    fn no_call_drops() {
        assert_eq!(
            reconcile_call(&call('-', '-'), &snp(b'A', b'G'), false, false),
            ReconcileResult::DropNoCall
        );
    }
}
