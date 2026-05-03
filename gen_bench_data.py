import os, random

NIND = 100
NSNP = 10000
OVERLAP = 5000
SEED = 42

os.makedirs("/tmp/reigen_bench", exist_ok=True)
rng = random.Random(SEED)

def pack_record(gs):
    n = len(gs)
    out = bytearray((n * 2 + 7) // 8)
    for i, g in enumerate(gs):
        two = {0:0b00, 1:0b01, 2:0b10}.get(g, 0b11)
        byte = i // 4
        shift = 6 - 2 * (i % 4)
        out[byte] |= two << shift
    return bytes(out)

def gen_dataset(prefix, nind, nsnp, snp_start_idx):
    # .ind
    with open(f"/tmp/reigen_bench/{prefix}.ind", "w") as f:
        for i in range(nind):
            f.write(f"{prefix}_S{i:03d} U Pop\n")
    
    # .snp
    with open(f"/tmp/reigen_bench/{prefix}.snp", "w") as f:
        for i in range(nsnp):
            idx = snp_start_idx + i
            chrom = (idx // 1000) + 1
            ppos = (idx % 1000) * 100 + 1
            f.write(f"rs{idx:06d} {chrom} 0.0 {ppos} A C\n")
            
    # .geno (PAM)
    rlen = max(48, (nind * 2 + 7) // 8)
    with open(f"/tmp/reigen_bench/{prefix}.geno", "wb") as f:
        header = f"GENO {nind} {nsnp} 0 0".encode()
        buf = bytearray(rlen)
        buf[:len(header)] = header
        f.write(bytes(buf))
        for i in range(nsnp):
            gs = [rng.choice([0, 1, 2, 9]) for _ in range(nind)]
            rec = pack_record(gs)
            out = bytearray(rlen)
            out[:len(rec)] = rec
            f.write(bytes(out))

gen_dataset("ds1", NIND, NSNP, 0)
gen_dataset("ds2", NIND, NSNP, NSNP - OVERLAP)

print(f"Generated datasets in /tmp/reigen_bench/")
