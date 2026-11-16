#!/usr/bin/env python3
"""
Masking test for "does SHAPEIT4 connect the read blocks correctly?"

The reads already tell us the true phase INSIDE each read block. We hide that
link in the middle of every big block (give the second half a new PS), let
SHAPEIT4 re-connect the two halves blind, and then check against the reads
whether it joined them the right way round.

  split:  python3 phase_junction_test.py split reads.vcf.gz reads_split.vcf.gz [MIN_SNPS]
  score:  python3 phase_junction_test.py score splits.tsv switch_errors.bed
"""
import sys
import statistics
from collections import defaultdict


def split(inp, out, min_snps=10):
    import pysam

    # Pass 1: which SNPs belong to which read block (PS)
    blocks = defaultdict(list)          # (chrom, PS) -> positions
    phased = []                         # (chrom, pos, PS) for gap statistics
    with pysam.VariantFile(inp) as vf:
        sample = list(vf.header.samples)[0]
        for rec in vf:
            call = rec.samples[sample]
            ps = call.get("PS") if call.phased else None
            if ps is not None:
                blocks[(rec.chrom, ps)].append(rec.pos)
                phased.append((rec.chrom, rec.pos, ps))

    # Split point = first SNP of the second half (only blocks with >= min_snps SNPs)
    splits = {}
    for key, pos in blocks.items():
        if len(pos) >= min_snps:
            pos = sorted(pos)
            half = len(pos) // 2
            splits[key] = (pos[half], pos[half] - pos[half - 1], half, len(pos) - half)

    # Pass 2: relabel the second half of each split block with a new PS
    # (new PS = position of its first SNP, the same convention WhatsHap uses)
    with pysam.VariantFile(inp) as vf, pysam.VariantFile(out, "wz", header=vf.header) as vo:
        sample = list(vf.header.samples)[0]
        for rec in vf:
            call = rec.samples[sample]
            ps = call.get("PS") if call.phased else None
            key = (rec.chrom, ps)
            if ps is not None and key in splits and rec.pos >= splits[key][0]:
                call["PS"] = splits[key][0]
            vo.write(rec)

    with open(out + ".splits.tsv", "w") as fh:
        fh.write("chrom\tPS\tsplit_pos\tgap_bp\tsnps_left\tsnps_right\n")
        for (chrom, ps), (mid, gap, nl, nr) in sorted(splits.items()):
            fh.write(f"{chrom}\t{ps}\t{mid}\t{gap}\t{nl}\t{nr}\n")

    # Gaps at REAL junctions: consecutive phased SNPs that sit in different blocks
    phased.sort()
    real_gaps = [b[1] - a[1] for a, b in zip(phased, phased[1:]) if a[0] == b[0] and a[2] != b[2]]
    test_gaps = [v[1] for v in splits.values()]
    print(f"read blocks: {len(blocks)} | blocks split (>= {min_snps} SNPs): {len(splits)}")
    if test_gaps:
        print(f"median gap at test splits : {statistics.median(test_gaps):,.0f} bp")
    if real_gaps:
        print(f"median gap at real joins  : {statistics.median(real_gaps):,.0f} bp")


def score(splits_tsv, bed):
    splits = []
    with open(splits_tsv) as fh:
        next(fh)
        for line in fh:
            chrom, ps, mid, gap, nl, nr = line.split("\t")
            splits.append((chrom, int(mid), int(gap)))

    # A BED record (chrom, pos_before, pos_after) marks a switch between two
    # consecutive compared SNPs. A split was joined wrongly if a switch spans it.
    switches = defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            f = line.split("\t")
            switches[f[0]].append((int(f[1]), int(f[2])))

    wrong = [(c, m, g) for c, m, g in splits
             if any(s < m <= e for s, e in switches.get(c, []))]
    n, k = len(splits), len(wrong)
    print(f"test junctions            : {n}")
    print(f"joined the WRONG way      : {k}")
    if n:
        print(f"junction error rate       : {100 * k / n:.2f}%")
    if wrong:
        print(f"median gap at wrong joins : {statistics.median(g for _, _, g in wrong):,.0f} bp")


if __name__ == "__main__":
    if len(sys.argv) >= 4 and sys.argv[1] == "split":
        split(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 10)
    elif len(sys.argv) == 4 and sys.argv[1] == "score":
        score(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)