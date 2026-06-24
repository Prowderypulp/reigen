//! PGEN record decoders + LUTs.
//!
//! Decodes the main 2-bit hardcall track for record types 0–7, skipping
//! auxiliary tracks (multiallelic, phase, dosage). Built on top of
//! [`difflist`] (varints + difflist) and [`header`] (type/length/spans).
//!
//! # LUT: PGEN ↔ canonical AM
//!
//! PGEN packs samples **LSB-first** (sample 0 at bits 1-0). Canonical AM
//! packs **MSB-first** (sample 0 at bits 7-6). The 2-bit _values_ are the
//! same ({0b00=0, 0b01=1, 0b10=2, 0b11=missing}) because PGEN counts ALT
//! alleles and AM counts allele1, and PVAR ALT = AM allele1.
//!
//! So the LUT is a pure bit-order reversal — no semantic recoding needed
//! (unlike PLINK→AM which swaps 00↔10).

use anyhow::{bail, Result};

use super::difflist::{decode_difflist, ByteCursor};

// ======================================================================
// LUTs: full-byte PGEN (LSB-first) ↔ canonical AM (MSB-first)
// ======================================================================

/// PGEN byte → canonical AM byte. Bit-order reverse only.
pub static PGEN_TO_AM: [u8; 256] = build_pgen_to_am();

/// Canonical AM byte → PGEN byte.
pub static AM_TO_PGEN: [u8; 256] = build_am_to_pgen();

const fn build_pgen_to_am() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        let s0 = b & 0b11;
        let s1 = (b >> 2) & 0b11;
        let s2 = (b >> 4) & 0b11;
        let s3 = (b >> 6) & 0b11;
        t[i] = (s0 << 6) | (s1 << 4) | (s2 << 2) | s3;
        i += 1;
    }
    t
}

const fn build_am_to_pgen() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        let s0 = (b >> 6) & 0b11;
        let s1 = (b >> 4) & 0b11;
        let s2 = (b >> 2) & 0b11;
        let s3 = b & 0b11;
        t[i] = s0 | (s1 << 2) | (s2 << 4) | (s3 << 6);
        i += 1;
    }
    t
}

// ======================================================================
// Record decoders
// ======================================================================

const RTYPE_MASK: u8 = 0x07;

/// Decode one variant record into `genovec` (N per-sample PGEN values).
///
/// `anchor` holds the last decoded non-LD variant. Updated on types 0,1,4,6,7.
/// `genovec` receives N u8 values in PGEN encoding {0,1,2,3}.
/// `cursor` must be positioned at the first byte of the variant record.
pub fn decode_record(
    cursor: &mut ByteCursor,
    record_type: u8,
    n_samples: usize,
    anchor: &mut [u8],
    genovec: &mut [u8],
) -> Result<()> {
    let cat = record_type & RTYPE_MASK;
    match cat {
        0 => decode_uncompressed(cursor, n_samples, anchor, genovec),
        1 => decode_onebit(cursor, n_samples, anchor, genovec),
        2 => decode_ld(cursor, n_samples, anchor, genovec),
        3 => decode_ld_inverted(cursor, n_samples, anchor, genovec),
        4 | 6 | 7 => decode_difflist_common(cursor, cat, n_samples, anchor, genovec),
        5 => bail!("pgen record type 5: reserved, should not appear"),
        _ => unreachable!(),
    }
}

/// Type 0: uncompressed 2-bit array, `ceil(N/4)` bytes, LSB-first.
fn decode_uncompressed(
    cursor: &mut ByteCursor,
    n_samples: usize,
    anchor: &mut [u8],
    genovec: &mut [u8],
) -> Result<()> {
    let n_bytes = n_samples.div_ceil(4);
    let raw = cursor.take(n_bytes)?;
    for s in 0..n_samples {
        let byte = raw[s / 4];
        genovec[s] = (byte >> (2 * (s % 4))) & 0b11;
    }
    anchor[..n_samples].copy_from_slice(&genovec[..n_samples]);
    Ok(())
}

/// Decode a type-0 record where we must consume exactly `total_len` bytes
/// (skipping any aux tracks beyond the main 2-bit array). Used when the
/// header provides record spans.
pub fn decode_uncompressed_total(
    cursor: &mut ByteCursor,
    total_len: usize,
    n_samples: usize,
    anchor: &mut [u8],
    genovec: &mut [u8],
) -> Result<()> {
    let n_bytes = n_samples.div_ceil(4);
    if n_bytes > total_len {
        bail!("pgen record: expected at least {n_bytes} data bytes, span is {total_len}");
    }
    let raw = cursor.take(n_bytes)?;
    for s in 0..n_samples {
        let byte = raw[s / 4];
        genovec[s] = (byte >> (2 * (s % 4))) & 0b11;
    }
    let skipped = total_len - n_bytes;
    if skipped > 0 {
        cursor.take(skipped)?;
    }
    anchor[..n_samples].copy_from_slice(&genovec[..n_samples]);
    Ok(())
}

/// Type 1: 1-bit array + difflist.
///
/// Mode byte layout: bits 1-0 = diff = (set - unset) mod 4,
/// bits 3-2 = unset_value. Then set_value = (unset + diff) & 3.
fn decode_onebit(
    cursor: &mut ByteCursor,
    n_samples: usize,
    anchor: &mut [u8],
    genovec: &mut [u8],
) -> Result<()> {
    let mode = cursor.take_one()?;
    let unset = (mode >> 2) & 0b11;
    let diff = mode & 0b11;
    let set = (unset + diff) & 0b11;

    if set == unset {
        bail!("pgen record type 1: mode byte specifies set=unset={set}");
    }

    let n_bit_bytes = n_samples.div_ceil(8);
    let bits = cursor.take(n_bit_bytes)?;
    for s in 0..n_samples {
        let bit = (bits[s / 8] >> (s % 8)) & 1;
        genovec[s] = if bit == 0 { unset } else { set };
    }

    let dl = decode_difflist(cursor, n_samples as u32, true)?;
    for (k, &sid) in dl.sample_ids.iter().enumerate() {
        if (sid as usize) < n_samples {
            if let Some(ref cats) = dl.categories {
                genovec[sid as usize] = cats[k] & 0b11;
            }
        }
    }

    anchor[..n_samples].copy_from_slice(&genovec[..n_samples]);
    Ok(())
}

/// Type 2: LD-compressed. Start from the anchor; apply difflist with categories.
fn decode_ld(
    cursor: &mut ByteCursor,
    n_samples: usize,
    anchor: &[u8],
    genovec: &mut [u8],
) -> Result<()> {
    genovec[..n_samples].copy_from_slice(&anchor[..n_samples]);
    let dl = decode_difflist(cursor, n_samples as u32, true)?;
    for (k, &sid) in dl.sample_ids.iter().enumerate() {
        if (sid as usize) < n_samples {
            if let Some(ref cats) = dl.categories {
                genovec[sid as usize] = cats[k] & 0b11;
            }
        }
    }
    Ok(())
}

/// Type 3: LD-compressed inverted. Same as type 2, then swap 0↔2.
fn decode_ld_inverted(
    cursor: &mut ByteCursor,
    n_samples: usize,
    anchor: &[u8],
    genovec: &mut [u8],
) -> Result<()> {
    decode_ld(cursor, n_samples, anchor, genovec)?;
    for g in genovec.iter_mut().take(n_samples) {
        match *g & 0b11 {
            0 => *g = 2,
            2 => *g = 0,
            _ => {}
        }
    }
    Ok(())
}

/// Types 4, 6, 7: difflist from common value.
/// Common: type 4 → 0 (hom-REF), type 6 → 2 (hom-ALT), type 7 → 3 (missing).
fn decode_difflist_common(
    cursor: &mut ByteCursor,
    category: u8,
    n_samples: usize,
    anchor: &mut [u8],
    genovec: &mut [u8],
) -> Result<()> {
    let common = match category {
        4 => 0u8,
        6 => 2u8,
        7 => 3u8,
        _ => unreachable!(),
    };
    genovec[..n_samples].fill(common);
    let dl = decode_difflist(cursor, n_samples as u32, true)?;
    for (k, &sid) in dl.sample_ids.iter().enumerate() {
        if (sid as usize) < n_samples {
            if let Some(ref cats) = dl.categories {
                genovec[sid as usize] = cats[k] & 0b11;
            }
        }
    }
    anchor[..n_samples].copy_from_slice(&genovec[..n_samples]);
    Ok(())
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lut_inverse() {
        for i in 0..=255u8 {
            assert_eq!(
                AM_TO_PGEN[PGEN_TO_AM[i as usize] as usize],
                i,
                "PGEN→AM→PGEN roundtrip failed for byte {i:08b}"
            );
        }
    }

    #[test]
    fn lut_preserves_values_but_reverses_order() {
        let pgen_byte = 0b11_10_01_00u8;
        let am_byte = PGEN_TO_AM[pgen_byte as usize];
        assert_eq!(am_byte, 0b00_01_10_11);
    }

    #[test]
    fn decode_uncompressed_roundtrip() {
        let vals: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let n: usize = 8;
        let rec_bytes = n.div_ceil(4);

        let mut pgen_bytes = vec![0u8; rec_bytes];
        for (i, &v) in vals.iter().enumerate() {
            let byte = i / 4;
            let shift = 2 * (i % 4);
            pgen_bytes[byte] |= v << shift;
        }

        let mut cursor = ByteCursor::new(&pgen_bytes);
        let mut anchor = vec![0u8; n];
        let mut genovec = vec![0u8; n];

        decode_record(&mut cursor, 0, n, &mut anchor, &mut genovec).unwrap();

        assert_eq!(&genovec, &vals);
        assert_eq!(&anchor, &vals);
    }

    #[test]
    fn decode_onebit() {
        // N=8. Mode: unset=0, set=1, diff=1. 1-bit array: sample 3 = 1.
        // Difflist: sample 7→3 (missing).
        let mode_byte = 0x01u8; // unset=0 (bits3-2), diff=1
        let bit_byte = 1u8 << 3;

        // Difflist: L=1, group first s=7, category=3.
        let mut record = vec![mode_byte, bit_byte];
        record.push(0x01); // L=1
        record.push(0x07); // group first s=7
        record.push(0b11); // category 3
        // no deltas (only entry is group-first)

        let mut cursor = ByteCursor::new(&record);
        let mut anchor = vec![0u8; 8];
        let mut genovec = vec![0u8; 8];

        decode_record(&mut cursor, 1, 8, &mut anchor, &mut genovec).unwrap();
        assert_eq!(genovec, vec![0, 0, 0, 1, 0, 0, 0, 3]);
    }

    #[test]
    fn decode_difflist_common_type4() {
        // N=4, common=0. Difflist: samp 1→2, samp 3→1.
        let mut dl_bytes = vec![0x02u8]; // L=2
        dl_bytes.push(0x01); // group first s=1
        // categories: entry0=2 (10), entry1=1 (01) → LSB-first: (1<<2)|2 = 0x06
        dl_bytes.push(0x06);
        dl_bytes.push(0x02); // delta: s3 - s1 = 2

        let mut cursor = ByteCursor::new(&dl_bytes);
        let mut anchor = vec![0u8; 4];
        let mut genovec = vec![0u8; 4];

        decode_difflist_common(&mut cursor, 4, 4, &mut anchor, &mut genovec).unwrap();
        assert_eq!(genovec, vec![0, 2, 0, 1]);
        assert_eq!(&anchor[..4], &genovec[..4]);
    }

    #[test]
    fn decode_ld_inverted_swaps_0_and_2() {
        let anchor = vec![0u8, 1, 2, 3];
        // LD difflist: samp 0→2, samp 2→0.
        let mut dl_bytes = vec![0x02u8]; // L=2
        dl_bytes.push(0x00); // group first s=0
        // categories: cat0=2, cat1=0 → (0<<2)|2 = 2
        dl_bytes.push(0x02);
        dl_bytes.push(0x02); // delta: s2 - s0 = 2

        let mut cursor = ByteCursor::new(&dl_bytes);
        let mut genovec = vec![0u8; 4];
        let orig_anchor = anchor.clone();

        decode_ld_inverted(&mut cursor, 4, &anchor, &mut genovec).unwrap();

        // LD: anchor[0,1,2,3] with patch [s0→2, s2→0] → [2,1,0,3]
        // Invert 0↔2 → [0,1,2,3] = original anchor.
        assert_eq!(&genovec, &orig_anchor);
    }

    #[test]
    fn decode_type5_errors() {
        let mut cursor = ByteCursor::new(&[]);
        let mut anchor = vec![0u8; 4];
        let mut genovec = vec![0u8; 4];
        let err = decode_record(&mut cursor, 5, 4, &mut anchor, &mut genovec).unwrap_err();
        assert!(format!("{err}").contains("reserved"));
    }
}
