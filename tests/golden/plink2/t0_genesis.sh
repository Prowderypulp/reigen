#!/bin/bash
# T0 fixture genesis script — run ONCE with real plink2 to generate
# tests/golden/plink2/t0_lsb.{pgen,pvar,psam}
# Do NOT regenerate unless the fixture definition changes.
set -euo pipefail

python3 -c '
import struct
genos = [2, 1, 0, 3]
lsb_byte = (genos[3] << 6) | (genos[2] << 4) | (genos[1] << 2) | genos[0]
header = b"\x6c\x1b\x02" + struct.pack("<I", 1) + struct.pack("<I", 4) + b"\x80"
with open("tests/golden/plink2/t0_lsb.pgen", "wb") as f:
    f.write(header + struct.pack("<B", lsb_byte))
with open("tests/golden/plink2/t0_lsb.pvar", "w") as f:
    f.write("#CHROM\tPOS\tID\tREF\tALT\n1\t1000\trs0\tG\tA\n")
with open("tests/golden/plink2/t0_lsb.psam", "w") as f:
    f.write("#FID\tIID\tSEX\n")
    for i in range(4):
        f.write(f"1\tS{i}\tNA\n")
print("T0 fixture written.")
'

echo "Verifying with plink2..."
plink2 --pfile tests/golden/plink2/t0_lsb --export A --out /tmp/t0_verify 2>/dev/null
cat /tmp/t0_verify.raw
