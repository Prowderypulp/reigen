use crate::strand;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SnpKey {
    pub chrom: u8,
    pub pos: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlipDecision {
    Match,
    Flip,
}

/// Reconcile kit alleles against reference alleles.
pub fn reconcile(
    in_a1: u8,
    in_a2: u8,
    ref_a1: u8,
    ref_a2: u8,
    flipstrand: bool,
    allow_ambiguous: bool,
) -> Option<FlipDecision> {
    if !allow_ambiguous && strand::is_ambiguous(ref_a1, ref_a2) {
        return None;
    }

    match strand::decide_flip(in_a1, in_a2, ref_a1, ref_a2, flipstrand) {
        Some(false) => Some(FlipDecision::Match),
        Some(true) => Some(FlipDecision::Flip),
        None => None,
    }
}
