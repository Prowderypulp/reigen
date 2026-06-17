#!/bin/bash
# =============================================================================
# fixtures_genesis.sh — regenerate the PLINK2 golden fixtures (task C1)
# =============================================================================
#
#   *** DO NOT RUN THIS CASUALLY. ***
#
# These fixtures are an EXTERNAL ORACLE produced by real `plink2`, NOT by
# reigen. They pin reigen's PGEN reader/writer against the reference tool so a
# future change cannot silently re-invert genotype polarity (the B-013 class of
# bug). Regenerating them with anything other than real plink2 — and ESPECIALLY
# regenerating them from reigen's own output — defeats the entire test and must
# never happen. Only re-run this script if:
#   (a) the fixture *definitions* below intentionally change, AND
#   (b) you have a real `plink2` binary (this was built with v2.0.0-a.6.13), AND
#   (c) Gary has signed off.
# After running, re-read tests/golden/plink2/T0_FINDINGS.md and update the
# provenance / record-type distribution section to match the new output.
#
# Produces (all committed under tests/golden/plink2/):
#   gold.{pgen,pvar,psam}   — primary fixture, from the committed golden .bed
#                             (8 samples x 12 SNPs; same dataset as convertf).
#   ld.{pgen,pvar,psam}     — compressible fixture (400x400) that exercises the
#                             LD / difflist / 1-bit decoders. Built from a
#                             generated VCF (ld_src.vcf, also written here).
#   multi.{pgen,pvar,psam}  — multiallelic fixture (6x4) with one comma-ALT
#                             variant, for the "drop multiallelic" path. Built
#                             from a generated VCF (multi_src.vcf).
#
# Requirements: real `plink2` on PATH, python3. Run from the repo root.
# =============================================================================
set -euo pipefail

HERE="tests/golden/plink2"
BED="tests/golden/plink/gold"   # existing committed golden PLINK1 fileset

command -v plink2 >/dev/null || { echo "ERROR: plink2 not on PATH"; exit 1; }
[ -f "${BED}.bed" ] || { echo "ERROR: ${BED}.bed not found (run from repo root)"; exit 1; }

echo "plink2 version: $(plink2 --version | head -1)"

# -----------------------------------------------------------------------------
# 1. Primary golden fixture — straight from the committed golden .bed.
#    Same 8-sample/12-SNP dataset as the convertf fixtures, so allele set and
#    missing calls match. Yields storage mode 0x10 (variable-width) with all
#    records uncompressed type-0 (dataset too small/dense to compress).
# -----------------------------------------------------------------------------
echo "==> gold.{pgen,pvar,psam}"
plink2 --bfile "${BED}" --make-pgen --out "${HERE}/gold" >/dev/null 2>&1
rm -f "${HERE}/gold.log"

# -----------------------------------------------------------------------------
# 2. Compressible / LD fixture — generated VCF -> --make-pgen.
#    The generator (seed 1337, deterministic) builds:
#      * LD blocks (a base variant + a few near-duplicate variants)   -> type 2
#      * LD-inverted blocks (base + its 0<->2 complement)             -> type 3
#      * very-sparse-ALT variants                                     -> type 4
#      * very-dense-ALT (mostly homALT) variants                      -> type 6
#      * mostly-missing variants                                      -> type 7
#      * plain common/uncompressible variants                  -> type 0 / 1
#    Verified below: all of plink2's hardcall record types (0,1,2,3,4,6,7)
#    actually appear. Type 5 is reserved/unused by the spec.
# -----------------------------------------------------------------------------
echo "==> ld.{pgen,pvar,psam} (generating ld_src.vcf)"
python3 - "$HERE/ld_src.vcf" <<'PYEOF'
import random, sys
random.seed(1337)
N, M = 400, 400
out = sys.argv[1]
samples = [f"S{i}" for i in range(N)]
def draw(maf):
    if random.random() < 0.002: return None
    return (1 if random.random()<maf else 0)+(1 if random.random()<maf else 0)
def invert(d):
    return None if d is None else {0:2,1:1,2:0}[d]
variants=[]; v=0
while v < M:
    bt = random.random()
    if bt < 0.30:                                   # LD block -> type 2
        maf = random.choice([random.uniform(0.05,0.30), random.uniform(0.70,0.92)])
        base=[draw(maf) for _ in range(N)]; variants.append(base[:]); v+=1
        for _ in range(random.randint(2,6)):
            if v>=M: break
            dup=base[:]
            for _ in range(random.randint(1,8)): dup[random.randrange(N)] = draw(maf)
            variants.append(dup); v+=1
    elif bt < 0.45:                                 # LD-inverted -> type 3
        maf = random.uniform(0.10,0.40)
        base=[draw(maf) for _ in range(N)]; variants.append(base[:]); v+=1
        if v<M:
            inv=[invert(d) for d in base]
            for _ in range(random.randint(1,4)): inv[random.randrange(N)] = draw(maf)
            variants.append(inv); v+=1
    elif bt < 0.60:                                 # sparse ALT -> type 4
        variants.append([draw(random.uniform(0.003,0.02)) for _ in range(N)]); v+=1
    elif bt < 0.74:                                 # dense ALT  -> type 6
        variants.append([draw(random.uniform(0.97,0.998)) for _ in range(N)]); v+=1
    elif bt < 0.86:                                 # mostly miss-> type 7
        variants.append([None if random.random()<0.97 else draw(0.5) for _ in range(N)]); v+=1
    else:                                           # common     -> type 0/1
        variants.append([draw(random.uniform(0.1,0.5)) for _ in range(N)]); v+=1
variants=variants[:M]
def gt(d): return "./." if d is None else {0:"0/0",1:"0/1",2:"1/1"}[d]
with open(out,"w") as f:
    f.write("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"x\">\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##contig=<ID=1>\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples)+"\n")
    for i,var in enumerate(variants):
        f.write(f"1\t{(i+1)*100}\trs{i}\tG\tA\t.\tPASS\t.\tGT\t"+"\t".join(gt(d) for d in var)+"\n")
print(f"wrote {out} ({N} samples x {M} variants)")
PYEOF
plink2 --vcf "${HERE}/ld_src.vcf" --make-pgen --out "${HERE}/ld" >/dev/null 2>&1
rm -f "${HERE}/ld.log"

# -----------------------------------------------------------------------------
# 3. Multiallelic fixture — generated VCF with one triallelic (comma-ALT)
#    variant -> --make-pgen. plink2 preserves the multiallelic record (PVAR
#    ALT="A,T", PGEN vrtype bit 3 set), exercising the reader's drop path.
# -----------------------------------------------------------------------------
echo "==> multi.{pgen,pvar,psam} (generating multi_src.vcf)"
python3 - "$HERE/multi_src.vcf" <<'PYEOF'
import sys
out = sys.argv[1]
samples=[f"S{i}" for i in range(6)]
rows=[
 ("1","1000","rs0","G","A",   ["0/0","0/1","1/1","0/0","./.","0/1"]),
 ("1","2000","rs1","C","A,T", ["0/1","1/2","2/2","0/0","0/2","1/1"]),  # MULTIALLELIC
 ("1","3000","rs2","T","C",   ["1/1","0/0","0/1","./.","1/1","0/0"]),
 ("2","4000","rs3","A","G",   ["0/0","0/1","0/0","1/1","0/1","./."]),
]
with open(out,"w") as f:
    f.write("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"x\">\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##contig=<ID=1>\n##contig=<ID=2>\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples)+"\n")
    for c,p,i,r,a,gts in rows:
        f.write(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t.\tGT\t"+"\t".join(gts)+"\n")
print(f"wrote {out} (6 samples x 4 variants, rs1 triallelic)")
PYEOF
plink2 --vcf "${HERE}/multi_src.vcf" --make-pgen --out "${HERE}/multi" >/dev/null 2>&1
rm -f "${HERE}/multi.log"

# -----------------------------------------------------------------------------
# 4. Verify & report achieved record-type distribution (no faking — if a type
#    is absent, the report shows it). Parses the PGEN header directly.
# -----------------------------------------------------------------------------
report_rtypes() {
  python3 - "$1" <<'PYEOF'
import struct, sys, collections
d=open(sys.argv[1],"rb").read()
assert d[0:2]==b"\x6c\x1b"
mode=d[2]
print(f"  storage mode = 0x{mode:02x}")
if mode not in (0x10,0x11):
    print("  fixed-width: all records uncompressed (type 0)"); sys.exit()
M=struct.unpack("<I",d[3:7])[0]; N=struct.unpack("<I",d[7:11])[0]; b11=d[11]
tb = 8 if (b11&0x0f)>=4 else 4
print(f"  M={M} variants, N={N} samples, byte11=0x{b11:02x} ({tb}-bit types)")
B=(M+0xffff)//0x10000; off=12+B*8
if tb==4:
    nb=(M+1)//2; raw=d[off:off+nb]
    ts=[(raw[i//2]>>(0 if i%2==0 else 4))&0x0f for i in range(M)]
else:
    ts=list(d[off:off+M])
cat=collections.Counter(t&7 for t in ts)
nm={0:"0 uncompressed",1:"1 onebit",2:"2 LD",3:"3 LD-inverted",
    4:"4 difflist-homREF",5:"5 reserved",6:"6 difflist-dblALT",7:"7 difflist-missing"}
for c in range(8):
    if cat[c]: print(f"    type {nm[c]:22s}: {cat[c]}")
multi=sum(1 for t in ts if t&0x08)
if multi: print(f"    (bit3 multiallelic set on {multi} record(s))")
PYEOF
}
echo; echo "=== gold.pgen ===";  report_rtypes "${HERE}/gold.pgen"
echo;  echo "=== ld.pgen ===";    report_rtypes "${HERE}/ld.pgen"
echo;  echo "=== multi.pgen ==="; report_rtypes "${HERE}/multi.pgen"
echo
echo "Done. Remember: these are an EXTERNAL ORACLE — never regenerate from reigen."
