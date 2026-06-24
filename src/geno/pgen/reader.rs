//! PGEN reader — implements [`GenoReader`] for `.pgen`.
//!
//! Supports storage modes 0x01 (PLINK 1 .bed, delegates to `packed_ped`),
//! 0x02 (fixed-width), 0x10/0x11 (variable-width). Modes 0x03/0x04 (dosage)
//! are recognized but only the main hardcall track is decoded; dosage and
//! phase tracks are skipped.
//!
//! # Multiallelic handling
//!
//! Variants whose PVAR ALT has a comma are dropped by the pipeline before the
//! reader is opened (the `keep_ondisk` mask reflects this). As a defense-in-depth
//! measure, the reader also checks vrtype bit 3 (multiallelic hardcalls present)
//! and drops such variants with a warning. Dropped variants are skipped but
//! counted; the reader records them in a drop report accessible via
//! [`PgenReader::multiallelic_dropped`].

use anyhow::{bail, Context, Result};
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

use crate::geno::codec;
use crate::geno::{GenoReader, Layout};

use super::difflist::ByteCursor;
use super::header::{PgenHeader, StorageMode};
use super::record::{decode_record, decode_uncompressed_total};

pub struct PgenReader {
    mmap: Mmap,
    n_samples: usize,
    ondisk_variants: usize,
    kept_variants: usize,
    rec_bytes: usize,

    keep_ondisk: Vec<bool>,
    record_types: Vec<u8>,
    record_spans: Vec<(u64, u64)>,
    fixed_rec_width: usize,
    data_start: u64,

    anchor: Vec<u8>,
    genovec: Vec<u8>,

    next_ondisk_idx: usize,
    next_kept_idx: usize,

    multiallelic_dropped: Vec<usize>,
}

impl PgenReader {
    /// Open a `.pgen`, returning a reader that yields only variants where
    /// `keep_ondisk[i]` is true. `keep_ondisk.len()` must equal the on-disk
    /// variant count (the PVAR count); if `keep_ondisk` is empty, all
    /// variants are kept.
    pub fn open(path: &Path, keep_mask: Vec<bool>) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let mmap =
            unsafe { Mmap::map(&file) }.with_context(|| format!("mmap {}", path.display()))?;

        let header = PgenHeader::parse(&mmap).context("parsing PGEN header")?;

        match header.mode {
            StorageMode::Plink1Bed => {
                bail!("pgen mode 0x01: use PackedPedReader for PLINK 1 .bed files");
            }
            StorageMode::FixedDosage | StorageMode::FixedPhasedDosage => {
                bail!(
                    "pgen mode 0x{:02x}: dosage-only PGEN not supported (no hardcall track)",
                    mmap[2]
                );
            }
            _ => {}
        }

        let n_samples = header
            .sample_count
            .context("PGEN header missing sample count")? as usize;
        let ondisk_variants = header
            .variant_count
            .context("PGEN header missing variant count")? as usize;

        let rec_bytes = (n_samples * 2 + 7) / 8;

        let mask_provided = !keep_mask.is_empty();
        let mut keep_ondisk = if !mask_provided {
            vec![true; ondisk_variants]
        } else {
            if keep_mask.len() != ondisk_variants {
                bail!(
                    "keep_ondisk length {} != PGEN on-disk variant count {}",
                    keep_mask.len(),
                    ondisk_variants
                );
            }
            keep_mask
        };

        let record_types = header.record_types.clone();
        let record_spans = header.record_spans.clone();

        let (fixed_rec_width, data_start) = if header.mode == StorageMode::FixedWidth {
            let fw = n_samples.div_ceil(4);
            (fw, 12u64)
        } else {
            (0, 0)
        };

        // Defense-in-depth: scan for vrtype bit-3-set variants that pass
        // keep_ondisk. PVAR already drops comma-ALT variants; if one slips
        // through (no comma but bit 3 set), drop it here.
        let mut multiallelic_dropped = Vec::new();
        if !record_types.is_empty() {
            for i in 0..ondisk_variants {
                if record_types[i] & 0x08 != 0 && keep_ondisk[i] {
                    if mask_provided {
                        // The caller's PVAR-derived mask kept this variant, but
                        // the PGEN record marks it multiallelic. Silently
                        // dropping it would desync the genotype stream from the
                        // SnpRows — surface the disagreement instead.
                        bail!(
                            "pgen variant {i}: record marks multiallelic (vrtype bit 3) but the \
                             PVAR keep-mask kept it; PVAR and PGEN disagree on multiallelic status"
                        );
                    }
                    multiallelic_dropped.push(i);
                    keep_ondisk[i] = false;
                }
            }
        }

        let kept_variants = keep_ondisk.iter().filter(|&&k| k).count();

        Ok(PgenReader {
            mmap,
            n_samples,
            ondisk_variants,
            kept_variants,
            rec_bytes,
            keep_ondisk,
            record_types,
            record_spans,
            fixed_rec_width,
            data_start,
            anchor: vec![0u8; n_samples],
            genovec: vec![0u8; n_samples],
            next_ondisk_idx: 0,
            next_kept_idx: 0,
            multiallelic_dropped,
        })
    }

    /// List of on-disk variant indices dropped because vrtype bit 3 is set
    /// (multiallelic hardcalls). The pipeline's PVAR reader already drops
    /// comma-ALT variants, so this should normally be empty; it's a
    /// defense-in-depth check.
    pub fn multiallelic_dropped(&self) -> &[usize] {
        &self.multiallelic_dropped
    }
}

impl GenoReader for PgenReader {
    fn nind(&self) -> usize {
        self.n_samples
    }

    fn nsnp(&self) -> usize {
        self.kept_variants
    }

    fn layout(&self) -> Layout {
        Layout::SnpMajor
    }

    fn record_bytes(&self) -> usize {
        self.rec_bytes
    }

    fn read_record(&mut self, dst: &mut [u8]) -> Result<bool> {
        if dst.len() != self.rec_bytes {
            bail!(
                "dst len {} != record_bytes {}",
                dst.len(),
                self.rec_bytes
            );
        }
        if self.next_kept_idx >= self.kept_variants {
            return Ok(false);
        }

        // Decode on-disk variants in sequence until the next kept one. EVERY
        // record must be decoded, including dropped ones: an LD-compressed
        // record (type 2/3) diffs against the *previous non-LD variant in
        // on-disk order* (spec §"Main data track"), and a dropped variant
        // (multiallelic / comma-ALT) can be that anchor. Skipping its decode
        // would leave `anchor` stale and silently corrupt the following LD
        // record. So we decode to keep `anchor` live and gate only *emission*
        // on `keep_ondisk`.
        loop {
            if self.next_ondisk_idx >= self.ondisk_variants {
                bail!(
                    "pgen: exhausted on-disk variants at idx {} but only {} of {} kept variants read",
                    self.next_ondisk_idx,
                    self.next_kept_idx,
                    self.kept_variants
                );
            }

            let odx = self.next_ondisk_idx;
            self.next_ondisk_idx += 1;

            // Get the record bytes.
            let (off, total_len) = if self.fixed_rec_width > 0 {
                (self.data_start + odx as u64 * self.fixed_rec_width as u64, self.fixed_rec_width as u64)
            } else {
                let (o, l) = self
                    .record_spans
                    .get(odx)
                    .with_context(|| format!("pgen: variant {odx} missing from record_spans"))?;
                (*o, *l)
            };

            let end = off as usize + total_len as usize;
            if end > self.mmap.len() {
                bail!(
                    "pgen: variant {odx} record at offset {off} len {total_len} past EOF (file size {})",
                    self.mmap.len()
                );
            }
            let record_bytes = &self.mmap[off as usize..end];
            let mut cursor = ByteCursor::new(record_bytes);

            let rt = if !self.record_types.is_empty() {
                self.record_types[odx] & 0x07
            } else {
                0
            };

            if rt == 0 && total_len > 0 {
                decode_uncompressed_total(
                    &mut cursor,
                    total_len as usize,
                    self.n_samples,
                    &mut self.anchor,
                    &mut self.genovec,
                )?;
            } else {
                decode_record(
                    &mut cursor,
                    rt,
                    self.n_samples,
                    &mut self.anchor,
                    &mut self.genovec,
                )?;
            }

            // Dropped variants are decoded only to maintain the LD anchor;
            // emit (break) only when this variant is kept.
            if self.keep_ondisk[odx] {
                break;
            }
        }

        // Map PGEN values {0,1,2,3} → reigen canonical {0,1,2,9}.
        for s in 0..self.n_samples {
            self.genovec[s] = match self.genovec[s] & 0b11 {
                0 => 0,
                1 => 1,
                2 => 2,
                _ => codec::G_MISSING,
            };
        }

        // Pack per-sample values MSB-first into dst.
        dst[..self.rec_bytes].fill(0);
        codec::pack(&self.genovec[..self.n_samples], dst);

        self.next_kept_idx += 1;
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geno::codec;

    fn _unpack_record(r: &mut dyn GenoReader, n: usize) -> Vec<u8> {
        let mut buf = vec![0u8; r.record_bytes()];
        let mut out = vec![0u8; n];
        r.read_record(&mut buf).unwrap();
        codec::unpack(&buf, n, &mut out);
        out
    }

    #[test]
    fn read_gold_fixture_all_type0() {
        let path = format!("{}/tests/golden/plink2/gold.pgen", env!("CARGO_MANIFEST_DIR"));
        let mut r = PgenReader::open(Path::new(&path), vec![]).unwrap();

        assert_eq!(r.nind(), 8);
        assert_eq!(r.nsnp(), 12);
        assert_eq!(r.layout(), Layout::SnpMajor);

        // Read all 12 variants; should all decode without error.
        let mut count = 0;
        while r.read_record(&mut vec![0u8; r.record_bytes()]).unwrap() {
            count += 1;
        }
        assert_eq!(count, 12);
        assert!(r.multiallelic_dropped().is_empty());
    }

    #[test]
    fn read_ld_fixture_decodes_all_types() {
        let path = format!(
            "{}/tests/golden/plink2/ld.pgen",
            env!("CARGO_MANIFEST_DIR")
        );
        let mut r = PgenReader::open(Path::new(&path), vec![]).unwrap();
        assert_eq!(r.nind(), 400);
        assert_eq!(r.nsnp(), 400);

        let mut count = 0;
        while r.read_record(&mut vec![0u8; r.record_bytes()]).unwrap() {
            count += 1;
        }
        assert_eq!(count, 400);
    }

    #[test]
    fn read_t0_fixture_bit_order() {
        // t0_lsb: 1 variant, 4 samples, PGEN values [2,1,0,3]
        // Expected reigen values: [2,1,0,9]
        let path = format!(
            "{}/tests/golden/plink2/t0_lsb.pgen",
            env!("CARGO_MANIFEST_DIR")
        );
        let mut r = PgenReader::open(Path::new(&path), vec![]).unwrap();
        let out = _unpack_record(&mut r, 4);
        assert_eq!(out, vec![2, 1, 0, 9]);
    }

    #[test]
    fn keep_mask_drops_variants() {
        // gold fixture: 12 variants. Keep only 0, 2, 4, 6, 8, 10.
        let path = format!(
            "{}/tests/golden/plink2/gold.pgen",
            env!("CARGO_MANIFEST_DIR")
        );
        let mut keep = vec![false; 12];
        for i in (0..12).step_by(2) {
            keep[i] = true;
        }
        let mut r = PgenReader::open(Path::new(&path), keep).unwrap();
        assert_eq!(r.nsnp(), 6);

        let mut count = 0;
        while r.read_record(&mut vec![0u8; r.record_bytes()]).unwrap() {
            count += 1;
        }
        assert_eq!(count, 6);
    }

    #[test]
    fn multiallelic_fixture_drops_bit3_set() {
        // multi fixture: rs1 has vrtype bit 3 set and PVAR ALT has comma.
        let path = format!(
            "{}/tests/golden/plink2/multi.pgen",
            env!("CARGO_MANIFEST_DIR")
        );
        // keep_ondisk empty → all kept by the mask, but multi reader drops bit-3-set.
        let mut r = PgenReader::open(Path::new(&path), vec![]).unwrap();
        let mut count = 0;
        while r.read_record(&mut vec![0u8; r.record_bytes()]).unwrap() {
            count += 1;
        }
        assert_eq!(count, 3, "rs1 (multiallelic) should be dropped, leaving 3 biallelic");
        assert_eq!(r.multiallelic_dropped().len(), 1);
    }

    #[test]
    fn multiallelic_fixture_with_keep_mask_drops_comma_alt() {
        let path = format!(
            "{}/tests/golden/plink2/multi.pgen",
            env!("CARGO_MANIFEST_DIR")
        );
        // Simulate the PVAR reader: keep_ondisk = [true, false, true, true] (rs1 dropped).
        let keep = vec![true, false, true, true];
        let mut r = PgenReader::open(Path::new(&path), keep).unwrap();
        assert_eq!(r.nsnp(), 3);
        let mut count = 0;
        while r.read_record(&mut vec![0u8; r.record_bytes()]).unwrap() {
            count += 1;
        }
        assert_eq!(count, 3);
    }

    /// Regression: an LD-compressed variant must diff against the previous
    /// non-LD variant in *on-disk* order, even when that anchor variant is
    /// dropped by the keep mask. Hand-built mode-0x10 PGEN, 4 samples, 3
    /// variants: V0 type-0 [0,0,0,0] (kept), V1 type-0 [2,2,2,2] (DROPPED),
    /// V2 type-2 LD diffing V1 with sample 2 → cat 0, i.e. [2,2,0,2] (kept).
    /// If the reader skips decoding V1, V2 wrongly anchors on V0=[0,0,0,0]
    /// and decodes to [0,0,0,0]. With the anchor maintained, V2 → [2,2,0,2].
    #[test]
    fn ld_anchor_survives_dropped_variant() {
        let mut f: Vec<u8> = Vec::new();
        // Header: magic + mode 0x10, M=3, N=4, byte11=0x40 (4-bit types,
        // 1-byte lengths, no allele cts, none-provisional).
        f.extend_from_slice(&[0x6c, 0x1b, 0x10]);
        f.extend_from_slice(&3u32.to_le_bytes());
        f.extend_from_slice(&4u32.to_le_bytes());
        f.push(0x40);
        // Block offset (1 block): first record starts at byte 25
        // (12 header + 8 offset + 2 type bytes + 3 length bytes).
        f.extend_from_slice(&25u64.to_le_bytes());
        // Type array, 4-bit packed: V0=0,V1=0 → 0x00; V2=2 → 0x02.
        f.extend_from_slice(&[0x00, 0x02]);
        // Length array, 1 byte each: V0=1, V1=1, V2=3.
        f.extend_from_slice(&[0x01, 0x01, 0x03]);
        // V0 record: [0,0,0,0] LSB-first → 0x00.
        f.push(0x00);
        // V1 record: [2,2,2,2] → each 0b10 → 0xAA.
        f.push(0xAA);
        // V2 record: type-2 LD difflist {L=1, id=2 (uint8), cat=0}.
        f.extend_from_slice(&[0x01, 0x02, 0x00]);
        assert_eq!(f.len(), 30);

        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), &f).unwrap();

        // Drop V1 (the anchor); keep V0 and V2.
        let keep = vec![true, false, true];
        let mut r = PgenReader::open(tmp.path(), keep).unwrap();
        assert_eq!(r.nsnp(), 2);

        let v0 = _unpack_record(&mut r, 4);
        assert_eq!(v0, vec![0, 0, 0, 0]);
        let v2 = _unpack_record(&mut r, 4);
        assert_eq!(
            v2,
            vec![2, 2, 0, 2],
            "LD variant must anchor on the dropped V1, not V0"
        );
    }

    #[test]
    fn roundtrip_gold_through_pgen_reader_matches_plink1_read() {
        // Read gold fixture via PGEN reader, compare to PLINK1 .bed reader.
        let pgen_path = format!(
            "{}/tests/golden/plink2/gold.pgen",
            env!("CARGO_MANIFEST_DIR")
        );
        let bed_path = format!(
            "{}/tests/golden/plink/gold.bed",
            env!("CARGO_MANIFEST_DIR")
        );

        let mut pgen_r = PgenReader::open(Path::new(&pgen_path), vec![]).unwrap();
        let mut bed_r =
            crate::geno::packed_ped::PackedPedReader::open(Path::new(&bed_path), 8, 12).unwrap();

        let mut pbuf = vec![0u8; pgen_r.record_bytes()];
        let mut bbuf = vec![0u8; bed_r.record_bytes()];
        for _ in 0..12 {
            pgen_r.read_record(&mut pbuf).unwrap();
            bed_r.read_record(&mut bbuf).unwrap();
            assert_eq!(pbuf, bbuf, "PGEN vs BED record mismatch");
        }
    }
}
