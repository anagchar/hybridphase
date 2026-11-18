#!/usr/bin/env python3
"""
Synthetic input for bin/assembly_haplotag.py.

  asm.byname.bam     reads aligned to h1_* / h2_* contigs, grouped by read, with AS tags: a read that crosses
                     k het sites scores 6*k lower on the wrong haplotype (k = 0 -> tie)
  vcf_tagged.bam     the same reads on GRCh38 (chr21, chr22), with HP/PS as from VCF-based haplotagging,
                     CB cell barcodes, and one planted switch error on chr22 at 30 Mb

usage: make_synthetic_asm.py OUTDIR
"""
import os
import sys

import numpy as np
import pysam

rng = np.random.default_rng(11)
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
CHROMS = {"chr21": 46_709_983, "chr22": 50_818_468}
N = 8000

# mixed fragment lengths: single-cell short fragments + long reads
lengths = np.where(rng.random(N) < 0.5, rng.integers(150, 2000, N), rng.integers(3000, 60000, N))
chrom = rng.choice(list(CHROMS), N)
pos = np.array([rng.integers(5_000_000, CHROMS[c] - 100_000) for c in chrom])
true_hp = rng.integers(1, 3, N)
hets = rng.poisson(lengths / 1500)                       # het sites crossed

asm_hdr = {"HD": {"VN": "1.6", "SO": "unsorted"},
           "SQ": [{"SN": f"h{h}_{c}", "LN": l} for h in (1, 2) for c, l in CHROMS.items()]}
ref_hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": c, "LN": l} for c, l in CHROMS.items()]}

with pysam.AlignmentFile(os.path.join(OUT, "asm.byname.bam"), "wb", header=asm_hdr) as fa:
    for i in range(N):
        L = int(lengths[i])
        good = 2 * L - int(rng.integers(0, L // 5 + 1))
        bad = good - 6 * int(hets[i])
        first = true_hp[i] if hets[i] else rng.integers(1, 3)
        for k, h in enumerate([first, 3 - first]):
            r = pysam.AlignedSegment(fa.header)
            r.query_name = f"read{i}"
            r.reference_name = f"h{h}_{chrom[i]}"
            r.reference_start = int(pos[i])
            r.cigartuples = [(0, L)]
            r.mapping_quality = 60 if hets[i] else 0
            r.flag = 0 if k == 0 else 256
            r.set_tag("AS", good if h == true_hp[i] else bad)
            fa.write(r)

order = np.lexsort((pos, [list(CHROMS).index(c) for c in chrom]))
with pysam.AlignmentFile(os.path.join(OUT, "vcf_tagged.bam"), "wb", header=ref_hdr) as fr:
    for i in order:
        r = pysam.AlignedSegment(fr.header)
        r.query_name = f"read{i}"
        r.reference_name = chrom[i]
        r.reference_start = int(pos[i])
        r.cigartuples = [(0, int(lengths[i]))]
        r.mapping_quality = 60
        r.set_tag("CB", f"CELL{rng.integers(0, 40):02d}")
        if hets[i]:
            hp = true_hp[i]
            if chrom[i] == "chr22" and pos[i] > 30_000_000:   # planted switch error in the phased VCF
                hp = 3 - hp
            r.set_tag("HP", int(hp))
            r.set_tag("PS", 1)
        fr.write(r)
pysam.index(os.path.join(OUT, "vcf_tagged.bam"))
print("synthetic assembly data in", OUT)
