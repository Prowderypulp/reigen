//! Transpose a 2-bit-packed genotype matrix between SnpMajor and SampleMajor
//! layouts.
//!
//! # Matrix layout
//!
//! Canonical encoding (see `geno::mod`): 2 bits per cell, MSB-first within
//! byte. Input: `rows` records of `ceil(cols * 2 / 8)` bytes.
//! Output: `cols` records of `ceil(rows * 2 / 8)` bytes.
//!
//! PACKEDANCESTRYMAP (SnpMajor) → TGENO (SampleMajor):
//! - rows = nsnp, cols = nind
//! - output rows = nind, output cols = nsnp
//!
//! # Algorithms
//!
//! Two implementations here:
//!
//! 1. `transpose_packed_scalar`: correctness reference. Cell-by-cell extract
//!    and insert. O(rows × cols) with significant per-cell cost. Fine for
//!    small matrices, unusable at AADR scale (~100 GB ops).
//!
//! 2. `transpose_packed_blocked`: the production path. Outer tiling over
//!    32-cell blocks (8 bytes in each dimension), with inner per-tile
//!    transpose that extracts source bytes once and writes destination
//!    bytes once. Parallelized over row-tiles via rayon.
//!
//! The blocked version is the default from `transpose_packed()`.
//!
//! # Correctness
//!
//! Unit tests compare blocked vs scalar on random matrices and check
//! double-transpose identity: `T(T(M)) == M`.

use anyhow::{bail, Result};
use rayon::prelude::*;

/// Block size in cells per axis. Must be a multiple of 4 so that column tiles
/// are 4-aligned and hit the gather4+LUT fast path (only the final edge tile,
/// when `cols % 4 != 0`, falls back to cell-by-cell). 128 is optimal for tall
/// AADR-shaped matrices.
pub const BLOCK: usize = 128;

/// Gather LUT: maps a source byte (4 packed 2-bit cells) to a u32 with
/// each column's cell pre-shifted to bits 7:6 of its respective byte.
///
/// Entry layout for input byte with cells [c0, c1, c2, c3]:
///   bits 31:30 = c0, bits 23:22 = c1, bits 15:14 = c2, bits 7:6 = c3
///
/// Usage: `GATHER_LUT[src_byte]` then shift right by `2 * row_index`
/// before OR-combining 4 rows into one output byte.
const GATHER_LUT: [u32; 256] = {
    let mut lut = [0u32; 256];
    let mut b = 0u32;
    while b < 256 {
        let c0 = (b >> 6) & 0b11;
        let c1 = (b >> 4) & 0b11;
        let c2 = (b >> 2) & 0b11;
        let c3 = b & 0b11;
        lut[b as usize] = (c0 << 30) | (c1 << 22) | (c2 << 14) | (c3 << 6);
        b += 1;
    }
    lut
};

/// Public entry point. Dispatches to blocked implementation.
pub fn transpose_packed(src: &[u8], rows: usize, cols: usize, dst: &mut [u8]) -> Result<()> {
    let row_bytes = (cols * 2 + 7) / 8;
    let col_bytes = (rows * 2 + 7) / 8;
    if src.len() < row_bytes * rows {
        bail!(
            "transpose: src len {} < rows({}) * row_bytes({})",
            src.len(),
            rows,
            row_bytes
        );
    }
    if dst.len() < col_bytes * cols {
        bail!(
            "transpose: dst len {} < cols({}) * col_bytes({})",
            dst.len(),
            cols,
            col_bytes
        );
    }
    transpose_packed_blocked(src, rows, cols, dst, row_bytes, col_bytes);
    Ok(())
}

// ======================================================================
// Cell-level accessors (used by scalar reference + blocked edge handling)
// ======================================================================

/// Extract the 2-bit cell at column index `c` from a row buffer.
/// MSB-first within byte.
#[inline(always)]
fn get_cell(row: &[u8], c: usize) -> u8 {
    let byte = c / 4;
    let shift = 6 - 2 * (c % 4);
    (row[byte] >> shift) & 0b11
}

/// Set the 2-bit cell at column index `c` in a row buffer.
/// The slot must already be zero (as in a freshly allocated dst buffer).
#[inline(always)]
fn set_cell(row: &mut [u8], c: usize, two_bit: u8) {
    let byte = c / 4;
    let shift = 6 - 2 * (c % 4);
    row[byte] |= (two_bit & 0b11) << shift;
}

// ======================================================================
// Scalar reference — correct, slow, used by tests
// ======================================================================

/// Naive cell-by-cell transpose. Correct for any `rows`/`cols`. Used as
/// test oracle.
#[allow(dead_code)]
pub fn transpose_packed_scalar(src: &[u8], rows: usize, cols: usize, dst: &mut [u8]) {
    let row_bytes = (cols * 2 + 7) / 8;
    let col_bytes = (rows * 2 + 7) / 8;
    for b in dst[..col_bytes * cols].iter_mut() {
        *b = 0;
    }

    for r in 0..rows {
        let src_row = &src[r * row_bytes..(r + 1) * row_bytes];
        for c in 0..cols {
            let v = get_cell(src_row, c);
            if v != 0 {
                let dst_row = &mut dst[c * col_bytes..(c + 1) * col_bytes];
                set_cell(dst_row, r, v);
            }
        }
    }
}

// ======================================================================
// Blocked transpose — production path
// ======================================================================

/// Blocked transpose. Partitions the output into row-tiles of `BLOCK` output
/// rows (= `BLOCK` input columns), processes tiles in parallel via rayon.
/// Each worker writes a disjoint output slice — no locking needed.
fn transpose_packed_blocked(
    src: &[u8],
    rows: usize,
    cols: usize,
    dst: &mut [u8],
    row_bytes: usize,
    col_bytes: usize,
) {
    for b in dst[..col_bytes * cols].iter_mut() {
        *b = 0;
    }

    // Partition output by output-row blocks. Each output row corresponds
    // to one input column, so an output-row block of `BLOCK` rows writes
    // to `cols` output rows of... wait, reconsider:
    //
    // dst layout: dst[c * col_bytes .. (c+1) * col_bytes] is output column c
    // (which == input row c if we think of the original matrix, but since
    // we're transposing, output column c corresponds to *input column c*
    // only when rows == cols. In general:
    //
    // M_in[r, c] == M_out[c, r]
    //
    // So output row c has `rows` cells, packed into `col_bytes` bytes.
    // `dst[c * col_bytes .. (c+1) * col_bytes]` is that output row.
    //
    // Partitioning by contiguous ranges of output rows (= input columns)
    // gives us disjoint output slices. Perfect for rayon.

    // Compute tile ranges in the output-row dimension (= input cols).
    let n_col_tiles = (cols + BLOCK - 1) / BLOCK;

    // Slice dst into chunks of `BLOCK * col_bytes` bytes (one tile's worth
    // of output rows), then par_iter over them.
    dst[..col_bytes * cols]
        .par_chunks_mut(BLOCK * col_bytes)
        .enumerate()
        .for_each(|(tile_idx, dst_tile)| {
            let col_start = tile_idx * BLOCK;
            let col_end = (col_start + BLOCK).min(cols);
            let _ = n_col_tiles;

            // For each row-tile of the input.
            let mut r_start = 0;
            while r_start < rows {
                let r_end = (r_start + BLOCK).min(rows);
                transpose_tile(
                    src, row_bytes, r_start, r_end, col_start, col_end, dst_tile, col_bytes,
                );
                r_start = r_end;
            }
        });
}

/// Transpose one rectangular tile of the matrix.
///
/// Reads `src[r_start..r_end]` rows, `src[.., col_start..col_end]` columns.
/// Writes into `dst_tile` which covers output rows `col_start..col_end` — so
/// within `dst_tile`, output row `c` starts at byte `(c - col_start) * col_bytes`.
///
/// Tile sizes ≤ BLOCK in each direction; may be smaller on matrix edges.
///
/// For aligned tiles (`col_start` and tile width both multiples of 4), uses the
/// gather4+LUT kernel: it reverses the inner loop to gather one byte from each
/// of 4 source rows and emit 4 sequential output bytes, turning the scatter
/// (4 random writes per source byte) into sequential writes through 4 output
/// rows. Falls back to cell-by-cell for edge tiles with non-4-aligned columns.
#[inline]
fn transpose_tile(
    src: &[u8],
    row_bytes: usize,
    r_start: usize,
    r_end: usize,
    col_start: usize,
    col_end: usize,
    dst_tile: &mut [u8],
    col_bytes: usize,
) {
    let tile_cols = col_end - col_start;
    let aligned = (col_start % 4 == 0) && (tile_cols % 4 == 0);

    if aligned {
        transpose_tile_gather4_lut(
            src, row_bytes, r_start, r_end, col_start, col_end, dst_tile, col_bytes,
        );
    } else {
        // Edge tile fallback: cell-by-cell.
        for r in r_start..r_end {
            let src_row = &src[r * row_bytes..(r + 1) * row_bytes];
            for c in col_start..col_end {
                let v = get_cell(src_row, c);
                if v != 0 {
                    // Output row index within the tile:
                    let tile_row = c - col_start;
                    let dst_row = &mut dst_tile[tile_row * col_bytes..(tile_row + 1) * col_bytes];
                    set_cell(dst_row, r, v);
                }
            }
        }
    }
}

/// Transpose one tile using gather4+LUT. Processes 4 input columns and 4 input
/// rows simultaneously, producing 4 complete output bytes per inner iteration.
/// Writes are sequential through the 4 active output rows.
///
/// Requires: `col_start % 4 == 0` and `(col_end - col_start) % 4 == 0`.
fn transpose_tile_gather4_lut(
    src: &[u8],
    row_bytes: usize,
    r_start: usize,
    r_end: usize,
    col_start: usize,
    col_end: usize,
    dst_tile: &mut [u8],
    col_bytes: usize,
) {
    let col_byte_start = col_start / 4;
    let col_byte_end = col_end / 4;

    for cb in col_byte_start..col_byte_end {
        let c_base = cb * 4;
        let tile_row0 = c_base - col_start;
        let dst_base0 = tile_row0 * col_bytes;
        let dst_base1 = (tile_row0 + 1) * col_bytes;
        let dst_base2 = (tile_row0 + 2) * col_bytes;
        let dst_base3 = (tile_row0 + 3) * col_bytes;

        // Main loop: process 4 input rows at a time.
        let mut r = r_start;
        while r + 4 <= r_end {
            let out_byte = r / 4;

            // Gather: read one byte from each of 4 source rows.
            let g0 = GATHER_LUT[src[r * row_bytes + cb] as usize];
            let g1 = GATHER_LUT[src[(r + 1) * row_bytes + cb] as usize];
            let g2 = GATHER_LUT[src[(r + 2) * row_bytes + cb] as usize];
            let g3 = GATHER_LUT[src[(r + 3) * row_bytes + cb] as usize];

            // Combine: each row contributes at a different bit position.
            let combined = g0 | (g1 >> 2) | (g2 >> 4) | (g3 >> 6);

            // Write: extract each output byte.
            let ob0 = (combined >> 24) as u8;
            let ob1 = (combined >> 16) as u8;
            let ob2 = (combined >> 8) as u8;
            let ob3 = combined as u8;

            if ob0 != 0 { dst_tile[dst_base0 + out_byte] |= ob0; }
            if ob1 != 0 { dst_tile[dst_base1 + out_byte] |= ob1; }
            if ob2 != 0 { dst_tile[dst_base2 + out_byte] |= ob2; }
            if ob3 != 0 { dst_tile[dst_base3 + out_byte] |= ob3; }

            r += 4;
        }

        // Remainder: handle leftover rows (< 4) when tile height isn't a
        // multiple of 4 (only the final row-tile when `rows % 4 != 0`).
        while r < r_end {
            let sb = src[r * row_bytes + cb];
            if sb != 0 {
                let dst_byte = r / 4;
                let dst_shift = 6 - 2 * (r % 4);
                dst_tile[dst_base0 + dst_byte] |= ((sb >> 6) & 0b11) << dst_shift;
                dst_tile[dst_base1 + dst_byte] |= ((sb >> 4) & 0b11) << dst_shift;
                dst_tile[dst_base2 + dst_byte] |= ((sb >> 2) & 0b11) << dst_shift;
                dst_tile[dst_base3 + dst_byte] |= (sb & 0b11) << dst_shift;
            }
            r += 1;
        }
    }
}

// ======================================================================
// Tests
// ======================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn build_matrix(rows: usize, cols: usize, seed: u64) -> Vec<u8> {
        let row_bytes = (cols * 2 + 7) / 8;
        let mut m = vec![0u8; rows * row_bytes];
        let mut state = seed.wrapping_mul(0x9E3779B97F4A7C15);
        for r in 0..rows {
            for c in 0..cols {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let v = (state >> 48) as u8 & 0b11;
                let byte = r * row_bytes + c / 4;
                let shift = 6 - 2 * (c % 4);
                m[byte] |= v << shift;
            }
        }
        m
    }

    fn matrices_equal(a: &[u8], b: &[u8], rows: usize, cols: usize, row_bytes: usize) -> bool {
        for r in 0..rows {
            let ra = &a[r * row_bytes..(r + 1) * row_bytes];
            let rb = &b[r * row_bytes..(r + 1) * row_bytes];
            for c in 0..cols {
                if get_cell(ra, c) != get_cell(rb, c) {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn scalar_tiny() {
        // 2 rows × 4 cols: [[0,1,2,3],[3,2,1,0]]
        // row0 byte: 00 01 10 11 = 0x1B
        // row1 byte: 11 10 01 00 = 0xE4
        let src = vec![0x1Bu8, 0xE4];
        let mut dst = vec![0u8; 4]; // 4 output rows × 1 byte each
        transpose_packed_scalar(&src, 2, 4, &mut dst);
        // out row 0: col0 of src = (0, 3) → MSB-first 2 cells: 00 11 _ _ = 0x30
        // out row 1: col1 of src = (1, 2) → 01 10 _ _ = 0x60
        // out row 2: col2 of src = (2, 1) → 10 01 _ _ = 0x90
        // out row 3: col3 of src = (3, 0) → 11 00 _ _ = 0xC0
        assert_eq!(dst, vec![0x30, 0x60, 0x90, 0xC0]);
    }

    #[test]
    fn blocked_matches_scalar_small() {
        let rows = 10;
        let cols = 15;
        let row_bytes = (cols * 2 + 7) / 8;
        let col_bytes = (rows * 2 + 7) / 8;
        let src = build_matrix(rows, cols, 42);
        let mut dst_scalar = vec![0u8; col_bytes * cols];
        let mut dst_blocked = vec![0u8; col_bytes * cols];
        transpose_packed_scalar(&src, rows, cols, &mut dst_scalar);
        transpose_packed(&src, rows, cols, &mut dst_blocked).unwrap();
        assert_eq!(dst_scalar, dst_blocked, "blocked differs from scalar");
        let _ = row_bytes;
    }

    #[test]
    fn blocked_matches_scalar_block_aligned() {
        let rows = BLOCK;
        let cols = BLOCK;
        let col_bytes = (rows * 2 + 7) / 8;
        let src = build_matrix(rows, cols, 1);
        let mut dst_s = vec![0u8; col_bytes * cols];
        let mut dst_b = vec![0u8; col_bytes * cols];
        transpose_packed_scalar(&src, rows, cols, &mut dst_s);
        transpose_packed(&src, rows, cols, &mut dst_b).unwrap();
        assert_eq!(dst_s, dst_b);
    }

    #[test]
    fn blocked_matches_scalar_multi_tile() {
        // 3 tiles in each direction — exercises parallel loop and tile
        // boundaries, some partial tiles.
        let rows = BLOCK * 3 + 5;
        let cols = BLOCK * 2 + 11;
        let col_bytes = (rows * 2 + 7) / 8;
        let src = build_matrix(rows, cols, 12345);
        let mut dst_s = vec![0u8; col_bytes * cols];
        let mut dst_b = vec![0u8; col_bytes * cols];
        transpose_packed_scalar(&src, rows, cols, &mut dst_s);
        transpose_packed(&src, rows, cols, &mut dst_b).unwrap();
        assert_eq!(dst_s, dst_b);
    }

    #[test]
    fn double_transpose_identity() {
        // T(T(M)) == M
        let rows = 73;
        let cols = 41;
        let row_bytes = (cols * 2 + 7) / 8;
        let col_bytes = (rows * 2 + 7) / 8;
        let src = build_matrix(rows, cols, 99);
        let mut once = vec![0u8; col_bytes * cols];
        transpose_packed(&src, rows, cols, &mut once).unwrap();
        // Second transpose: dimensions flip.
        let mut twice = vec![0u8; row_bytes * rows];
        transpose_packed(&once, cols, rows, &mut twice).unwrap();
        assert!(
            matrices_equal(&src, &twice, rows, cols, row_bytes),
            "T(T(M)) != M"
        );
    }

    #[test]
    fn rejects_short_buffers() {
        let src = vec![0u8; 10];
        let mut dst = vec![0u8; 10];
        let err = transpose_packed(&src, 100, 100, &mut dst).unwrap_err();
        let s = format!("{err:#}");
        assert!(s.contains("transpose"));
    }
}
