#!/usr/bin/env python3
"""
Masking test for "does the population phaser connect the read blocks correctly?"
(no truth set needed)

The reads already tell us the true phase INSIDE each read block. We hide that link in the
middle of every big block (the second half gets a new PS), let the population phaser
re-connect the two halves blind, and check against the reads whether it joined them the
right way round.

  split:  junction_test.py split reads.vcf.gz reads_split.vcf.gz [--min-snps 10]
          -> reads_split.vcf.gz + reads_split.vcf.gz.splits.tsv
  score:  junction_test.py score splits.tsv reads.vcf.gz rejoined.vcf.gz junctions.tsv
          rejoined.vcf.gz = split blocks re-joined with the population phase only
"""
import argparse
import statistics
import sys
from collections import defaultdict

import pysam


def phased_calls(path):
    """(chrom, pos) -> (gt, PS) for phased het calls"""
    out = {}
    with pysam.VariantFile(path) as vf:
        for rec in vf:
            call = rec.samples[0]
            gt = call.get("GT")
            ps = call.get("PS") if "PS" in call else None
            if gt is None or None in gt or len(gt) != 2 or gt[0] == gt[1] or not call.phased or ps is None:
                continue
            out[(rec.chrom, rec.pos)] = (gt, ps)
    return out


def split(inp, out, min_snps):
    blocks = defaultdict(list)
    phased = []
    for (chrom, pos), (gt, ps) in phased_calls(inp).items():
        blocks[(chrom, ps)].append(pos)
        phased.append((chrom, pos, ps))

    splits = {}
    for key, pos in blocks.items():
        if len(pos) >= min_snps:
            pos = sorted(pos)
            half = len(pos) // 2
            splits[key] = (pos[half], pos[half] - pos[half - 1], half, len(pos) - half)

    with pysam.VariantFile(inp) as vf, pysam.VariantFile(out, "wz", header=vf.header) as vo:
        for rec in vf:
            call = rec.samples[0]
            ps = call.get("PS") if (call.phased and "PS" in call) else None
            key = (rec.chrom, ps)
            if ps is not None and key in splits and rec.pos >= splits[key][0]:
                call["PS"] = splits[key][0]
            vo.write(rec)

    with open(out + ".splits.tsv", "w") as fh:
        fh.write("chrom\tps\tsplit_pos\tgap_bp\tsnps_left\tsnps_right\n")
        for (chrom, ps), (mid, gap, nl, nr) in sorted(splits.items()):
            fh.write(f"{chrom}\t{ps}\t{mid}\t{gap}\t{nl}\t{nr}\n")

    phased.sort()
    real_gaps = [b[1] - a[1] for a, b in zip(phased, phased[1:]) if a[0] == b[0] and a[2] != b[2]]
    test_gaps = [v[1] for v in splits.values()]
    print(f"read blocks: {len(blocks)} | blocks split (>= {min_snps} SNPs): {len(splits)}", file=sys.stderr)
    if test_gaps:
        print(f"median gap at test splits : {statistics.median(test_gaps):,.0f} bp", file=sys.stderr)
    if real_gaps:
        print(f"median gap at real joins  : {statistics.median(real_gaps):,.0f} bp", file=sys.stderr)


def score(splits_tsv, reads_vcf, rejoined_vcf, out):
    reads = phased_calls(reads_vcf)
    rej = phased_calls(rejoined_vcf)
    members = defaultdict(list)                     # original block -> positions
    for (chrom, pos), (gt, ps) in reads.items():
        members[(chrom, ps)].append(pos)

    rows = []
    with open(splits_tsv) as fh:
        next(fh)
        for line in fh:
            chrom, ps, mid, gap, nl, nr = line.rstrip("\n").split("\t")
            mid = int(mid)
            sides = {"L": [], "R": []}
            for pos in members.get((chrom, int(ps)), []):
                r = rej.get((chrom, pos))
                if r is None:
                    continue
                same = r[0][0] == reads[(chrom, pos)][0][0]
                sides["L" if pos < mid else "R"].append((same, r[1]))
            if not sides["L"] or not sides["R"]:
                status = "not_joined"
            else:
                pss = {p for side in sides.values() for _, p in side}
                if len(pss) != 1:
                    status = "not_joined"
                else:
                    o = {k: sum(s for s, _ in v) * 2 > len(v) for k, v in sides.items()}
                    status = "correct" if o["L"] == o["R"] else "wrong"
            rows.append((chrom, ps, mid, gap, nl, nr, status))

    with open(out, "w") as fh:
        fh.write("chrom\tps\tsplit_pos\tgap_bp\tsnps_left\tsnps_right\tstatus\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    n = len(rows)
    k = sum(r[-1] == "wrong" for r in rows)
    j = sum(r[-1] != "not_joined" for r in rows)
    print(f"test junctions: {n} | re-joined: {j} | joined the WRONG way: {k}"
          + (f" ({100 * k / j:.2f}%)" if j else ""), file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("reads")
    s.add_argument("out")
    s.add_argument("--min-snps", type=int, default=10)
    c = sub.add_parser("score")
    c.add_argument("splits")
    c.add_argument("reads")
    c.add_argument("rejoined")
    c.add_argument("out")
    a = ap.parse_args()
    if a.cmd == "split":
        split(a.reads, a.out, a.min_snps)
    else:
        score(a.splits, a.reads, a.rejoined, a.out)


if __name__ == "__main__":
    main()
