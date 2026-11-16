#!/usr/bin/env python3
"""
Per-chromosome phasing statistics and block list of a phased VCF (independent of WhatsHap's
column names, so the report always reads the same format).

usage: phase_stats.py phased.vcf.gz stats.tsv blocklist.tsv
"""
import sys
from collections import OrderedDict

import pysam


def n50(spans):
    s = sorted(spans, reverse=True)
    half, acc = sum(s) / 2, 0
    for x in s:
        acc += x
        if acc >= half:
            return x
    return 0


def main(vcf, stats_out, blocks_out):
    chroms = OrderedDict()
    with pysam.VariantFile(vcf) as vf:
        for rec in vf:
            call = rec.samples[0]
            gt = call.get("GT")
            if gt is None or None in gt or len(gt) != 2 or gt[0] == gt[1]:
                continue
            c = chroms.setdefault(rec.chrom, dict(het=0, phased=0, blocks=OrderedDict()))
            c["het"] += 1
            ps = call.get("PS") if "PS" in call else None
            if call.phased and ps is not None:
                c["phased"] += 1
                b = c["blocks"].setdefault(ps, [rec.pos, rec.pos, 0])
                b[1] = rec.pos
                b[2] += 1

    cols = ["chrom", "het_snvs", "phased", "phased_pct", "blocks", "singletons", "block_n50_bp",
            "largest_block_bp", "largest_block_snvs", "phased_span_bp"]
    tot = dict(het=0, phased=0, spans=[], sizes=[], singles=0)
    with open(stats_out, "w") as fs, open(blocks_out, "w") as fb:
        fs.write("\t".join(cols) + "\n")
        fb.write("chrom\tps\tstart\tend\tspan_bp\tn_snvs\n")
        for chrom, c in chroms.items():
            bl = [b for b in c["blocks"].values() if b[2] >= 2]
            spans = [b[1] - b[0] + 1 for b in bl]
            big = max(bl, key=lambda b: b[1] - b[0]) if bl else [0, -1, 0]
            singles = sum(b[2] == 1 for b in c["blocks"].values())
            fs.write("\t".join(map(str, [
                chrom, c["het"], c["phased"], f"{100 * c['phased'] / c['het']:.3f}" if c["het"] else "0",
                len(bl), singles, n50(spans), big[1] - big[0] + 1, big[2], sum(spans)])) + "\n")
            for ps, b in c["blocks"].items():
                fb.write(f"{chrom}\t{ps}\t{b[0]}\t{b[1]}\t{b[1] - b[0] + 1}\t{b[2]}\n")
            tot["het"] += c["het"]
            tot["phased"] += c["phased"]
            tot["spans"] += spans
            tot["sizes"] += [b[2] for b in bl]
            tot["singles"] += singles
        sp = tot["spans"]
        fs.write("\t".join(map(str, [
            "ALL", tot["het"], tot["phased"], f"{100 * tot['phased'] / tot['het']:.3f}" if tot["het"] else "0",
            len(sp), tot["singles"], n50(sp), max(sp) if sp else 0,
            tot["sizes"][sp.index(max(sp))] if sp else 0, sum(sp)])) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
