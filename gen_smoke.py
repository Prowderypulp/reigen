#!/usr/bin/env python3
"""
Generate a small synthetic PACKEDANCESTRYMAP dataset and the EIGENSTRAT
text file it should convert to. Used to smoke-test convertf-rs without
needing the upstream C convertf binary.

Output:
    /tmp/smoke/input.geno   (PACKEDANCESTRYMAP, 2-bit MSB-first)
    /tmp/smoke/input.snp
    /tmp/smoke/input.ind
    /tmp/smoke/expected.geno.txt  (EIGENSTRAT text — the ground truth)

Encoding: 0/1/2 = allele count, 9 = missing; stored 2-bit as 00/01/10/11
MSB-first. rlen = max(48, ceil(nind*2/8)).
"""
import os, random

NIND = 7           # intentionally not a multiple of 4 — tests tail handling
NSNP = 12
SEED = 42

os.makedirs("/tmp/smoke", exist_ok=True)
rng = random.Random(SEED)

def pack_record(gs):
    """gs: list of 0/1/2/9. Returns canonical 2-bit MSB-first bytes."""
    n = len(gs)
    out = bytearray((n * 2 + 7) // 8)
    for i, g in enumerate(gs):
        two = {0:0b00, 1:0b01, 2:0b10}.get(g, 0b11)
        byte = i // 4
        shift = 6 - 2 * (i % 4)
        out[byte] |= two << shift
    return bytes(out)

def gen_record(i):
    # Deterministic pattern mix + some missing.
    return [ ([0,1,2,9][(i + j*3) % 4]) for j in range(NIND) ]

def ascii_line(gs):
    return "".join({0:"0",1:"1",2:"2"}.get(g, "9") for g in gs) + "\n"

# --- .ind ---
with open("/tmp/smoke/input.ind", "w") as f:
    for i in range(NIND):
        pop = "PopA" if i < 3 else "PopB"
        sex = "MFU"[i % 3]
        f.write(f"SAMPLE{i:03d} {sex} {pop}\n")

# --- .snp ---
with open("/tmp/smoke/input.snp", "w") as f:
    for i in range(NSNP):
        chrom = (i % 22) + 1
        gpos = i * 0.001
        ppos = 1000 + i * 100000
        f.write(f"rs{i:04d} {chrom} {gpos:.6f} {ppos} A C\n")

# --- .geno (PAM) + expected EIGENSTRAT ---
rlen = max(48, (NIND * 2 + 7) // 8)
record_bytes = (NIND * 2 + 7) // 8

with open("/tmp/smoke/input.geno", "wb") as f:
    header = f"GENO {NIND} {NSNP} 0 0".encode()
    buf = bytearray(rlen)
    buf[:len(header)] = header
    f.write(bytes(buf))
    for i in range(NSNP):
        gs = gen_record(i)
        rec = pack_record(gs)
        out = bytearray(rlen)
        out[:len(rec)] = rec
        f.write(bytes(out))

with open("/tmp/smoke/expected.geno.txt", "w") as f:
    for i in range(NSNP):
        f.write(ascii_line(gen_record(i)))

print(f"generated nind={NIND} nsnp={NSNP} rlen={rlen} rec_bytes={record_bytes}")
print("files in /tmp/smoke/")
