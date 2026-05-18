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

#[derive(Debug, Clone, Copy)]
pub struct ReconcileOpts {
    pub flip_strand: bool,
    pub allow_ambiguous: bool,
    pub allow_flip_reference: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReconcileError {
    Ambiguous,
    InvalidAllele,
    Unresolvable,
}

/// Reconcile kit alleles against reference alleles.
pub fn reconcile(
    in_a1: u8,
    in_a2: u8,
    ref_a1: u8,
    ref_a2: u8,
    opts: ReconcileOpts,
) -> Result<FlipDecision, ReconcileError> {
    let in_a1 = in_a1.to_ascii_uppercase();
    let in_a2 = in_a2.to_ascii_uppercase();
    let ref_a1 = ref_a1.to_ascii_uppercase();
    let ref_a2 = ref_a2.to_ascii_uppercase();

    if !strand::is_acgt_allele(in_a1)
        || !strand::is_acgt_allele(in_a2)
        || !strand::is_acgt_allele(ref_a1)
        || !strand::is_acgt_allele(ref_a2)
    {
        return Err(ReconcileError::InvalidAllele);
    }

    if !opts.allow_ambiguous && strand::is_ambiguous(ref_a1, ref_a2) {
        return Err(ReconcileError::Ambiguous);
    }

    match strand::decide_flip(in_a1, in_a2, ref_a1, ref_a2, opts.flip_strand) {
        Some(false) => Ok(FlipDecision::Match),
        Some(true) if opts.allow_flip_reference => Ok(FlipDecision::Flip),
        Some(true) => Err(ReconcileError::Unresolvable),
        None => Err(ReconcileError::Unresolvable),
    }
}

#[cfg(test)]
mod tests {
    use super::{reconcile, FlipDecision, ReconcileError, ReconcileOpts};

    #[test]
    fn disallows_reference_flip_when_requested() {
        let opts = ReconcileOpts {
            flip_strand: false,
            allow_ambiguous: true,
            allow_flip_reference: false,
        };
        assert_eq!(
            reconcile(b'C', b'A', b'A', b'C', opts),
            Err(ReconcileError::Unresolvable)
        );
    }

    #[test]
    fn rejects_non_acgt_alleles() {
        let opts = ReconcileOpts {
            flip_strand: false,
            allow_ambiguous: true,
            allow_flip_reference: true,
        };
        assert_eq!(
            reconcile(b'0', b'A', b'A', b'C', opts),
            Err(ReconcileError::InvalidAllele)
        );
    }

    #[test]
    fn allows_flip_with_default_behavior() {
        let opts = ReconcileOpts {
            flip_strand: false,
            allow_ambiguous: true,
            allow_flip_reference: true,
        };
        assert_eq!(
            reconcile(b'C', b'A', b'A', b'C', opts),
            Ok(FlipDecision::Flip)
        );
    }
}
