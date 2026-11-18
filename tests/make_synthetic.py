#!/usr/bin/env python3
"""
Synthetic data to test the helper scripts without a cluster.

Builds a known diploid phase on two chromosomes and derives from it:
  het VCF (unphased), read-phased VCF (short blocks, random orientation, rare internal errors),
  population-phased VCF (panel subset, long-range switch errors), QC tables for the report.

usage: make_synthetic.py OUTDIR
"""
import gzip
import os
import sys

import numpy as np
import pysam

rng = np.random.default_rng(7)
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
CHROMS = {"chr21": 46_709_983, "chr22": 50_818_468}
ALLELES = "ACGT"

header = pysam.VariantHeader()
for c, l in CHROMS.items():
    header.add_line(f"##contig=<ID={c},length={l}>")
header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
header.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">')
header.add_sample("SAMPLE")

sites = {}
for c, l in CHROMS.items():
    n = l // 3000
    pos = np.unique(rng.integers(5_000_000 if c == "chr21" else 16_000_000, l - 100_000, n))
    hap = rng.integers(0, 2, len(pos))                       # allele on truth haplotype 1
    sites[c] = (pos, hap)


def write(path, records):
    tmp = path[:-3]
    with pysam.VariantFile(tmp, "w", header=header) as vo:
        for c, p, gt, phased, ps in records:
            r = vo.new_record(contig=c, start=p - 1, stop=p, alleles=("A", "G"))
            r.samples["SAMPLE"]["GT"] = gt
            r.samples["SAMPLE"].phased = phased
            if ps is not None:
                r.samples["SAMPLE"]["PS"] = int(ps)
            vo.write(r)
    pysam.tabix_index(tmp, preset="vcf", force=True)


het, reads, pop = [], [], []
for c, (pos, hap) in sites.items():
    # read blocks: geometric number of SNPs, random orientation, 0.2% internal flips, 5% singletons
    i = 0
    while i < len(pos):
        k = int(rng.geometric(1 / 25))
        flip = rng.random() < 0.5
        ps = pos[i]
        for j in range(i, min(i + k, len(pos))):
            h = hap[j] ^ flip ^ (rng.random() < 0.002)
            if k == 1:
                reads.append((c, pos[j], (0, 1), False, None))
            else:
                reads.append((c, pos[j], (h, 1 - h), True, ps))
        i += k
    # population phase: 85% of SNPs in the panel, a long-range switch every ~8 Mb
    sw = np.sort(rng.integers(pos[0], pos[-1], max(1, (pos[-1] - pos[0]) // 8_000_000)))
    state = np.searchsorted(sw, pos) % 2
    for j, p in enumerate(pos):
        het.append((c, p, (0, 1), False, None))
        if rng.random() < 0.85:
            h = hap[j] ^ state[j] ^ (rng.random() < 0.003)
            pop.append((c, p, (h, 1 - h), True, None))

write(os.path.join(OUT, "het.vcf.gz"), het)
# truth from the "assembly": the true haplotypes, phased, no PS (like dipcall)
write(os.path.join(OUT, "truth.vcf.gz"), [(c, p, (int(h), 1 - int(h)), True, None)
                                          for c, (pos, hap) in sites.items() for p, h in zip(pos, hap)])
write(os.path.join(OUT, "reads.vcf.gz"), reads)
write(os.path.join(OUT, "pop.vcf.gz"), pop)

with open(os.path.join(OUT, "ref.fa.fai"), "w") as fh:
    for c, l in CHROMS.items():
        fh.write(f"{c}\t{l}\t0\t60\t61\n")

# QC inputs for the report
data = os.path.join(OUT, "data")
os.makedirs(data, exist_ok=True)
n = sum(len(p) for p, _ in sites.values())
with gzip.open(os.path.join(data, "het_calls.tsv.gz"), "wt") as fh, gzip.open(os.path.join(data, "het_kept.tsv.gz"), "wt") as fk:
    fh.write("chrom\tpos\tqual\tdp\tad_ref\tad_alt\n")
    fk.write("chrom\tpos\n")
    for c, (pos, _) in sites.items():
        for p in pos:
            dp = rng.poisson(30)
            alt = rng.binomial(dp, 0.5)
            q = rng.gamma(4, 6)
            fh.write(f"{c}\t{p}\t{q:.2f}\t{dp}\t{dp - alt}\t{alt}\n")
            if q >= 20 and dp >= 15 and dp and 0.25 <= alt / dp <= 0.75:
                fk.write(f"{c}\t{p}\n")
with open(os.path.join(data, "het_funnel.tsv"), "w") as fh:
    fh.write(f"step\tn\nall_records\t{int(n * 2.6)}\npass\t{int(n * 2.2)}\npass_biallelic_snps\t{int(n * 1.9)}\n"
             f"het\t{n}\nhom_alt\t{int(n * 0.7)}\nhet_kept\t{int(n * 0.8)}\n")

# haplotag counts: LOH on part of chr22
with open(os.path.join(data, "haplotag_counts.tsv"), "w") as fh:
    fh.write("chrom\tbin\thp\treads\tbases\n")
    for c, l in CHROMS.items():
        for b in range(l // 1_000_000):
            f1 = 0.85 if (c == "chr22" and 30 <= b < 42) else 0.5
            t = rng.poisson(900)
            h1 = rng.binomial(t, f1)
            un = rng.poisson(120)
            for hp, v in ((1, h1), (2, t - h1), (0, un)):
                fh.write(f"{c}\t{b}\t{hp}\t{v}\t{v * 15000}\n")
print("synthetic data in", OUT)
