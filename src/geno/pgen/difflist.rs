//! PGEN base-128 varints and difflists.
//!
//! A *difflist* is PGEN's sparse representation of "samples that differ from a
//! common value", used by record types 1 (1-bit), 2/3 (LD), and 4/6/7
//! (difflist-from-common). It encodes either a sorted list of sample IDs, or a
//! list of (sample ID, 2-bit genotype category) pairs.
//!
//! Layout (spec §"Difflists"), in order:
//! 1. A base-128 varint `L` = number of entries. `L == 0` ends the list.
//! 2. The first sample ID of each group of 64 (`G = ceil(L/64)` groups),
//!    stored as uint8/16/24/32 chosen by `N`:
//!    - `N <= 2^8`  -> uint8
//!    - `N <= 2^16` -> uint16
//!    - `N <= 2^24` -> uint24 (uint32 with the high byte dropped)
//!    - otherwise   -> uint32
//! 3. `G - 1` bytes; byte `j` is `(final-component byte size of group j) - 63`.
//!    The last group's size is not stored (it runs to the end).
//! 4. *If categories are carried*, a packed 2-bit array of `ceil(L/4)` bytes,
//!    **LSB-first** (entry 0 in the low 2 bits of byte 0).
//! 5. `L - G` base-128 varints: for every `0 < k < L` that is **not** a
//!    multiple of 64, `sampleID[k] - sampleID[k-1]`.
//!
//! Component 3 (per-group byte sizes) lets a reader skip a group without
//! decoding it. We decode every entry sequentially, so we read those bytes for
//! correct positioning but otherwise ignore their values.

use anyhow::{bail, Result};

/// Sample-ID width in a difflist's group-first array, chosen by sample count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleIdWidth {
    U8,
    U16,
    U24,
    U32,
}

impl SampleIdWidth {
    /// Width implied by the sample count `N` (spec §"Difflists" component 1).
    pub fn for_sample_count(n: u32) -> Self {
        if n <= (1u64 << 8) as u32 {
            SampleIdWidth::U8
        } else if n <= (1u64 << 16) as u32 {
            SampleIdWidth::U16
        } else if n as u64 <= (1u64 << 24) {
            SampleIdWidth::U24
        } else {
            SampleIdWidth::U32
        }
    }

    pub fn bytes(self) -> usize {
        match self {
            SampleIdWidth::U8 => 1,
            SampleIdWidth::U16 => 2,
            SampleIdWidth::U24 => 3,
            SampleIdWidth::U32 => 4,
        }
    }
}

/// A little cursor over a byte slice. Tracks position so callers can learn how
/// many bytes a difflist consumed (needed to advance within a variant record).
pub struct ByteCursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> ByteCursor<'a> {
    pub fn new(buf: &'a [u8]) -> Self {
        ByteCursor { buf, pos: 0 }
    }

    pub fn pos(&self) -> usize {
        self.pos
    }

    pub fn remaining(&self) -> usize {
        self.buf.len() - self.pos
    }

    pub fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.buf.len() {
            bail!(
                "pgen difflist: unexpected end of buffer (need {} bytes at offset {}, have {})",
                n,
                self.pos,
                self.buf.len()
            );
        }
        let s = &self.buf[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    pub fn take_one(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    /// Copy the next `dst.len()` bytes into `dst`, advancing the cursor.
    pub fn read_into(&mut self, dst: &mut [u8]) -> Result<()> {
        let src = self.take(dst.len())?;
        dst.copy_from_slice(src);
        Ok(())
    }

    /// Read a base-128 (LEB128, protobuf-style) varint as a `u64`.
    ///
    /// Each byte contributes 7 low bits; the high bit (0x80) is a continuation
    /// flag. Little-endian group order (first byte = least significant).
    pub fn read_varint(&mut self) -> Result<u64> {
        let mut result: u64 = 0;
        let mut shift: u32 = 0;
        loop {
            let byte = self.take_one()?;
            if shift >= 64 {
                bail!("pgen varint: overflow (more than 10 bytes)");
            }
            result |= ((byte & 0x7f) as u64) << shift;
            if byte & 0x80 == 0 {
                break;
            }
            shift += 7;
        }
        Ok(result)
    }

    /// Read a fixed-width little-endian sample ID of the given width.
    fn read_sample_id(&mut self, width: SampleIdWidth) -> Result<u32> {
        let bytes = self.take(width.bytes())?;
        let mut v: u32 = 0;
        for (i, &b) in bytes.iter().enumerate() {
            v |= (b as u32) << (8 * i);
        }
        Ok(v)
    }
}

/// A decoded difflist.
///
/// `sample_ids` is the sorted list of differing sample IDs. `categories`, when
/// present, is the parallel list of 2-bit genotype categories ({0,1,2,3} in
/// PGEN order — the caller maps to canonical values).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiffList {
    pub sample_ids: Vec<u32>,
    pub categories: Option<Vec<u8>>,
}

impl DiffList {
    pub fn len(&self) -> usize {
        self.sample_ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sample_ids.is_empty()
    }
}

/// Decode a difflist from `cursor`, advancing it past the list.
///
/// - `n_samples` selects the group-first ID width.
/// - `has_categories` selects whether the packed 2-bit category array (spec
///   component 4) is present and decoded.
///
/// After return, `cursor.pos()` is just past the last byte of the difflist.
pub fn decode_difflist(
    cursor: &mut ByteCursor,
    n_samples: u32,
    has_categories: bool,
) -> Result<DiffList> {
    let l = cursor.read_varint()? as usize;
    if l == 0 {
        return Ok(DiffList {
            sample_ids: Vec::new(),
            categories: if has_categories { Some(Vec::new()) } else { None },
        });
    }

    let width = SampleIdWidth::for_sample_count(n_samples);
    let n_groups = l.div_ceil(64);

    // Component 2: group-first sample IDs.
    let mut group_firsts = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        group_firsts.push(cursor.read_sample_id(width)?);
    }

    // Component 3: G-1 per-group byte-size bytes. We don't need the values for
    // a sequential decode, but must consume them for correct positioning.
    let _group_sizes = cursor.take(n_groups - 1)?;

    // Component 4: packed 2-bit categories, LSB-first, ceil(L/4) bytes.
    let categories = if has_categories {
        let n_cat_bytes = l.div_ceil(4);
        let cat_bytes = cursor.take(n_cat_bytes)?;
        let mut cats = Vec::with_capacity(l);
        for k in 0..l {
            let byte = cat_bytes[k / 4];
            let shift = 2 * (k % 4);
            cats.push((byte >> shift) & 0b11);
        }
        Some(cats)
    } else {
        None
    };

    // Component 5: L-G varint deltas. Reconstruct sample IDs. Entry k is a
    // group-first iff k % 64 == 0 (taken from group_firsts); otherwise it is
    // prev + next-delta-varint.
    let mut sample_ids = Vec::with_capacity(l);
    for k in 0..l {
        if k % 64 == 0 {
            sample_ids.push(group_firsts[k / 64]);
        } else {
            let delta = cursor.read_varint()? as u32;
            let prev = sample_ids[k - 1];
            sample_ids.push(prev + delta);
        }
    }

    Ok(DiffList {
        sample_ids,
        categories,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn varint_small() {
        let mut c = ByteCursor::new(&[0x4f]);
        assert_eq!(c.read_varint().unwrap(), 79);
        assert_eq!(c.pos(), 1);
    }

    #[test]
    fn varint_multibyte() {
        // 5000 = 0x1388. Spec encodes it as 0x88 0x27.
        let mut c = ByteCursor::new(&[0x88, 0x27]);
        assert_eq!(c.read_varint().unwrap(), 5000);
        assert_eq!(c.pos(), 2);
    }

    #[test]
    fn varint_large() {
        // 39728178 -> base-128 LE. Verify round structure: 300 = 0xAC 0x02.
        let mut c = ByteCursor::new(&[0xac, 0x02]);
        assert_eq!(c.read_varint().unwrap(), 300);
    }

    #[test]
    fn sample_id_width_thresholds() {
        assert_eq!(SampleIdWidth::for_sample_count(256), SampleIdWidth::U8);
        assert_eq!(SampleIdWidth::for_sample_count(257), SampleIdWidth::U16);
        assert_eq!(SampleIdWidth::for_sample_count(65536), SampleIdWidth::U16);
        assert_eq!(SampleIdWidth::for_sample_count(65537), SampleIdWidth::U24);
        assert_eq!(SampleIdWidth::for_sample_count(488377), SampleIdWidth::U24);
        assert_eq!(
            SampleIdWidth::for_sample_count(16_777_216),
            SampleIdWidth::U24
        );
        assert_eq!(
            SampleIdWidth::for_sample_count(16_777_217),
            SampleIdWidth::U32
        );
    }

    /// The spec's worked example: N=488377, L=79, categories present, sample
    /// IDs 5000, 10000, ..., 395000 (step 5000). Reconstruct the exact bytes
    /// the spec lists and assert the decode.
    #[test]
    fn spec_worked_example() {
        let n: u32 = 488_377;
        let l: usize = 79;
        let expected_ids: Vec<u32> = (1..=l).map(|i| (i as u32) * 5000).collect();
        assert_eq!(expected_ids[0], 5000);
        assert_eq!(expected_ids[63], 320_000); // sample #63
        assert_eq!(expected_ids[64], 325_000); // group #1 first (#64)
        assert_eq!(expected_ids[78], 395_000);

        let mut bytes: Vec<u8> = Vec::new();

        // (1) L as varint = 0x4f.
        bytes.push(0x4f);

        // (2) group firsts as uint24 LE: #0 = 5000, #64 = 325000.
        // 5000 = 0x001388 -> 0x88 0x13 0x00
        // 325000 = 0x04F588 -> 0x88 0xf5 0x04
        bytes.extend_from_slice(&[0x88, 0x13, 0x00, 0x88, 0xf5, 0x04]);

        // (3) G-1 = 1 byte: group #0 final-component size (126) - 63 = 63 = 0x3f.
        bytes.push(0x3f);

        // (4) packed categories, ceil(79/4)=20 bytes, LSB-first. Use a varied
        // pattern so a wrong unpack would be caught. cat[k] = k % 4.
        let cats: Vec<u8> = (0..l).map(|k| (k % 4) as u8).collect();
        let mut cat_bytes = vec![0u8; l.div_ceil(4)];
        for (k, &cat) in cats.iter().enumerate() {
            cat_bytes[k / 4] |= cat << (2 * (k % 4));
        }
        assert_eq!(cat_bytes.len(), 20);
        bytes.extend_from_slice(&cat_bytes);

        // (5) L-G varint deltas. G = ceil(79/64) = 2, so 77 deltas, all 5000
        // = 0x88 0x27 (no entry for k=64, which is a group-first).
        let n_groups = l.div_ceil(64);
        assert_eq!(n_groups, 2);
        for _ in 0..(l - n_groups) {
            bytes.extend_from_slice(&[0x88, 0x27]);
        }

        let total_before = bytes.len();
        let mut cursor = ByteCursor::new(&bytes);
        let dl = decode_difflist(&mut cursor, n, true).unwrap();

        assert_eq!(cursor.pos(), total_before, "consumed entire difflist");
        assert_eq!(dl.sample_ids, expected_ids);
        assert_eq!(dl.categories.as_ref().unwrap(), &cats);
    }

    #[test]
    fn empty_difflist() {
        let mut c = ByteCursor::new(&[0x00]);
        let dl = decode_difflist(&mut c, 4, true).unwrap();
        assert!(dl.is_empty());
        assert_eq!(dl.categories.as_ref().unwrap().len(), 0);
        assert_eq!(c.pos(), 1);
    }

    #[test]
    fn difflist_without_categories() {
        // L=3, N=4 (uint8 IDs), one group, no category array.
        // IDs: 0, 2, 3. Deltas: (2-0)=2, (3-2)=1.
        let bytes = [0x03u8, 0x00, /*delta*/ 0x02, 0x01];
        let mut c = ByteCursor::new(&bytes);
        let dl = decode_difflist(&mut c, 4, false).unwrap();
        assert_eq!(dl.sample_ids, vec![0, 2, 3]);
        assert!(dl.categories.is_none());
        assert_eq!(c.pos(), bytes.len());
    }
}
