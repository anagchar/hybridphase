#!/usr/bin/env python3
"""
Stand-in for `whatshap compare --tsv-pairwise --switch-error-bed` (truth vs one VCF), for the tests only.
usage: mock_whatshap_compare.py truth.vcf.gz phased.vcf.gz name out.tsv out.bed
"""
import sys
from collections import defaultdict

import pysam

truth_f, vcf_f, name, out_tsv, out_bed = sys.argv[1:]
truth = {(r.chrom, r.pos): r.samples[0]["GT"] for r in pysam.VariantFile(truth_f)}
by_chrom = defaultdict(list)
for r in pysam.VariantFile(vcf_f):
    s = r.samples[0]
    gt = s.get("GT")
    if gt is None or gt[0] == gt[1] or not s.phased or s.get("PS") is None or (r.chrom, r.pos) not in truth:
        continue
    by_chrom[r.chrom].append((r.pos, s["PS"], gt[0] == truth[(r.chrom, r.pos)][0]))

with open(out_tsv, "w") as ft, open(out_bed, "w") as fb:
    ft.write("#sample\tchromosome\tdataset_name0\tdataset_name1\tall_assessed_pairs\tall_switches\t"
             "all_switch_rate\tblockwise_hamming\tblockwise_hamming_rate\tcovered_variants\n")
    for c, v in by_chrom.items():
        pairs = sw = ham = 0
        blocks = defaultdict(list)
        for (p0, ps0, o0), (p1, ps1, o1) in zip(v, v[1:]):
            if ps0 != ps1:
                continue
            pairs += 1
            if o0 != o1:
                sw += 1
                fb.write(f"{c}\t{p0}\t{p1}\n")
        for p, ps, o in v:
            blocks[ps].append(o)
        for b in blocks.values():
            ham += min(sum(b), len(b) - sum(b))
        ft.write(f"SAMPLE\t{c}\ttruth\t{name}\t{pairs}\t{sw}\t{sw / max(pairs, 1)}\t{ham}\t{ham / len(v)}\t{len(v)}\n")
