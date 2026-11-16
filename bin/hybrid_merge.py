#!/usr/bin/env python3
"""
Combine read-based phase blocks with a population-based phase (one chromosome or more).

Every read block (PS) is compared with the population phase at the SNPs both share:
  agree    = the read haplotype 1 carries the same allele as population haplotype 1
  disagree = the opposite
The block is then flipped as a whole to the majority orientation ("joined") if it shares at
least --min-snps SNPs with the population phase and the vote is not a tie. Joined blocks and
population-only SNPs get one phase set per chromosome; blocks that cannot be oriented keep
their own PS. Inside a read block the read phase is never changed.

  --blocks  per-block table (always written): how strongly the population phase supports
            each read block, which is a truth-free QC of both phasers.
  --out     merged VCF (omit to only write the table, e.g. when WhatsHap does the merge).

Without --pop only the block table (without population columns) is written.
"""
import argparse
import sys
from collections import OrderedDict

import pysam

COLS = ["chrom", "ps", "start", "end", "span_bp", "n_snps", "n_shared", "n_agree", "n_disagree",
        "concordance", "flipped", "joined"]


def het_gt(call):
    gt = call.get("GT")
    if gt is None or None in gt or len(gt) != 2 or gt[0] == gt[1]:
        return None
    return gt


def load_pop(path):
    pop = {}
    with pysam.VariantFile(path) as vf:
        for rec in vf:
            call = rec.samples[0]
            gt = het_gt(call)
            if gt is None or not call.phased or not rec.alts:
                continue
            pop[(rec.chrom, rec.pos, rec.ref, rec.alts[0])] = gt
    return pop


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", required=True, help="read-phased VCF with PS tags")
    ap.add_argument("--pop", help="population-phased VCF")
    ap.add_argument("--out", help="merged VCF (.vcf.gz)")
    ap.add_argument("--blocks", required=True, help="per-block orientation table (TSV)")
    ap.add_argument("--min-snps", type=int, default=2, help="min. shared SNPs to orient a block [2]")
    a = ap.parse_args()

    pop = load_pop(a.pop) if a.pop else {}

    # pass 1: per-block votes
    blocks = OrderedDict()
    with pysam.VariantFile(a.reads) as vf:
        for rec in vf:
            call = rec.samples[0]
            gt = het_gt(call)
            ps = call.get("PS") if "PS" in call else None
            if gt is None or not call.phased or ps is None:
                continue
            b = blocks.setdefault((rec.chrom, ps), dict(start=rec.pos, end=rec.pos, n=0, agree=0, disagree=0))
            b["end"] = rec.pos
            b["n"] += 1
            p = pop.get((rec.chrom, rec.pos, rec.ref, rec.alts[0])) if rec.alts else None
            if p is not None:
                if gt[0] == p[0]:
                    b["agree"] += 1
                else:
                    b["disagree"] += 1

    for b in blocks.values():
        shared = b["agree"] + b["disagree"]
        b["joined"] = bool(pop) and shared >= a.min_snps and b["agree"] != b["disagree"]
        b["flipped"] = b["joined"] and b["disagree"] > b["agree"]
        b["conc"] = max(b["agree"], b["disagree"]) / shared if shared else float("nan")

    with open(a.blocks, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for (chrom, ps), b in blocks.items():
            shared = b["agree"] + b["disagree"]
            row = [chrom, ps, b["start"], b["end"], b["end"] - b["start"] + 1, b["n"]]
            if pop:
                row += [shared, b["agree"], b["disagree"], f"{b['conc']:.4f}" if shared else "NA",
                        int(b["flipped"]), int(b["joined"])]
            else:
                row += ["NA"] * 6
            fh.write("\t".join(map(str, row)) + "\n")

    if not a.out:
        return
    if not pop:
        sys.exit("--out needs --pop")

    # pass 2: write merged VCF
    n_flip = n_pop_only = 0
    with pysam.VariantFile(a.reads) as vf:
        if "PS" not in vf.header.formats:
            vf.header.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">')
        vf.header.add_line("##hybridphase_merge=read blocks oriented by population phase "
                           f"(min shared SNPs {a.min_snps}); pop={a.pop}")
        chrom_ps = {}
        with pysam.VariantFile(a.out, "wz", header=vf.header) as vo:
            for rec in vf:
                call = rec.samples[0]
                gt = het_gt(call)
                ps = call.get("PS") if "PS" in call else None
                if gt is not None:
                    b = blocks.get((rec.chrom, ps)) if (call.phased and ps is not None) else None
                    if b is not None:
                        if b["joined"]:
                            cps = chrom_ps.setdefault(rec.chrom, rec.pos)
                            call["GT"] = (gt[1], gt[0]) if b["flipped"] else gt
                            call.phased = True
                            call["PS"] = cps
                            n_flip += b["flipped"]
                    else:
                        p = pop.get((rec.chrom, rec.pos, rec.ref, rec.alts[0])) if rec.alts else None
                        if p is not None:
                            cps = chrom_ps.setdefault(rec.chrom, rec.pos)
                            call["GT"] = p
                            call.phased = True
                            call["PS"] = cps
                            n_pop_only += 1
                vo.write(rec)

    nj = sum(b["joined"] for b in blocks.values())
    print(f"read blocks: {len(blocks)} | joined: {nj} | flipped SNPs: {n_flip} | population-only SNPs: {n_pop_only}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
