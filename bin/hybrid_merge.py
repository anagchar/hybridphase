#!/usr/bin/env python3
"""
Combine read-based phase blocks with a population-based phase (one chromosome or more).

Every read block (PS) is compared with the population phase at the SNPs both share:
  agree    = the read haplotype 1 carries the same allele as population haplotype 1
  disagree = the opposite
The block is then flipped as a whole to the majority orientation ("joined") if it shares at
least --min-snps SNPs with the population phase and the vote is not a tie. Joined blocks and
population-only SNPs get one phase set per chromosome; blocks that cannot be oriented keep
their own PS.

--block-error E (like SHAPEIT4 --use-PS): read blocks are soft constraints. Along each block the
agree/disagree pattern is decoded with a two-state model (orientation "as population" / "flipped"):
a switch inside the block costs E (expected error rate of the read phase sets), each disagreeing
SNP costs the population phaser's per-SNP error rate (estimated from the data). Where the decoded
orientation changes, the block is split and the parts are oriented separately, i.e. a long enough
run of SNPs where the population confidently disagrees overrides the reads. Single disagreeing
SNPs never do. Without --block-error the read phase inside a block is never changed.

  --blocks  per-block table (always written): how strongly the population phase supports
            each read block, which is a truth-free QC of both phasers.
  --out     merged VCF (omit to only write the table, e.g. when WhatsHap does the merge).

Without --pop only the block table (without population columns) is written.
"""
import argparse
import bisect
import math
import sys
from collections import OrderedDict

import pysam

COLS = ["chrom", "ps", "start", "end", "span_bp", "n_snps", "n_shared", "n_agree", "n_disagree",
        "concordance", "flipped", "joined", "parent_ps", "part"]


def viterbi(agree, eps, beta):
    """agree: list of bools along a block -> decoded state per SNP (0 = as population, 1 = flipped)"""
    ls, lw = math.log(1 - eps), math.log(eps)
    hit, miss = math.log(1 - beta), math.log(beta)
    em = lambda st, a: hit if a == (st == 0) else miss
    v = [em(0, agree[0]), em(1, agree[0])]
    back = []
    for a in agree[1:]:
        b0 = 0 if v[0] + ls >= v[1] + lw else 1
        b1 = 1 if v[1] + ls >= v[0] + lw else 0
        v = [max(v[0] + ls, v[1] + lw) + em(0, a), max(v[1] + ls, v[0] + lw) + em(1, a)]
        back.append((b0, b1))
    st = [0 if v[0] >= v[1] else 1]
    for b in reversed(back):
        st.append(b[st[-1]])
    return st[::-1]


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
    ap.add_argument("--block-error", type=float, default=None,
                    help="expected error rate of the read blocks (soft constraints, like SHAPEIT4 --use-PS); "
                         "omit to keep every read block intact")
    a = ap.parse_args()

    pop = load_pop(a.pop) if a.pop else {}

    # pass 1: per-block votes, in order along the block
    blocks = OrderedDict()
    with pysam.VariantFile(a.reads) as vf:
        for rec in vf:
            call = rec.samples[0]
            gt = het_gt(call)
            ps = call.get("PS") if "PS" in call else None
            if gt is None or not call.phased or ps is None:
                continue
            b = blocks.setdefault((rec.chrom, ps), dict(pos=[], votes=[]))
            b["pos"].append(rec.pos)
            p = pop.get((rec.chrom, rec.pos, rec.ref, rec.alts[0])) if rec.alts else None
            if p is not None:
                b["votes"].append((rec.pos, gt[0] == p[0]))

    # population per-SNP error rate: votes against the decoded orientation in well-covered blocks.
    # First estimate counts read switch errors as population errors too, so decode and re-estimate once.
    beta = None
    if pop and a.block_error:
        beta = 0.01
        for _ in range(2):
            mis = tot = 0
            for b in blocks.values():
                if len(b["votes"]) >= 10:
                    ag = [v for _, v in b["votes"]]
                    st = viterbi(ag, a.block_error, beta)
                    mis += sum(v != (x == 0) for v, x in zip(ag, st))
                    tot += len(ag)
            beta = min(0.2, max(0.002, mis / tot if tot else 0.01))

    # split blocks where the decoded orientation changes -> parts with their own votes
    parts = OrderedDict()                     # (chrom, ps) -> list of part dicts; bounds for lookup
    n_split = 0
    for key, b in blocks.items():
        bounds = [b["pos"][0]]
        if beta is not None and len(b["votes"]) >= 2:
            st = viterbi([v for _, v in b["votes"]], a.block_error, beta)
            for (p0, _), s0, (p1, _), s1 in zip(b["votes"], st, b["votes"][1:], st[1:]):
                if s1 != s0:
                    bounds.append(p1)
        n_split += len(bounds) > 1
        plist = []
        for k, lo in enumerate(bounds):
            hi = bounds[k + 1] if k + 1 < len(bounds) else float("inf")
            pos = [x for x in b["pos"] if lo <= x < hi]
            vt = [v for x, v in b["votes"] if lo <= x < hi]
            ag, dis = sum(vt), len(vt) - sum(vt)
            joined = bool(pop) and len(vt) >= a.min_snps and ag != dis
            plist.append(dict(start=pos[0], end=pos[-1], n=len(pos), agree=ag, disagree=dis, joined=joined,
                              flipped=joined and dis > ag, conc=max(ag, dis) / len(vt) if vt else float("nan"),
                              ps=key[1] if k == 0 else pos[0], part=k))
        parts[key] = (bounds, plist)

    with open(a.blocks, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for (chrom, ps), (_, plist) in parts.items():
            for b in plist:
                shared = b["agree"] + b["disagree"]
                row = [chrom, b["ps"], b["start"], b["end"], b["end"] - b["start"] + 1, b["n"]]
                if pop:
                    row += [shared, b["agree"], b["disagree"], f"{b['conc']:.4f}" if shared else "NA",
                            int(b["flipped"]), int(b["joined"])]
                else:
                    row += ["NA"] * 6
                fh.write("\t".join(map(str, row + [ps, b["part"]])) + "\n")

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
                           f"(min shared SNPs {a.min_snps}; block error {a.block_error or 'none'}); pop={a.pop}")
        chrom_ps = {}
        with pysam.VariantFile(a.out, "wz", header=vf.header) as vo:
            for rec in vf:
                call = rec.samples[0]
                gt = het_gt(call)
                ps = call.get("PS") if "PS" in call else None
                if gt is not None:
                    pb = parts.get((rec.chrom, ps)) if (call.phased and ps is not None) else None
                    b = pb[1][bisect.bisect_right(pb[0], rec.pos) - 1] if pb is not None else None
                    if b is not None:
                        if b["joined"]:
                            cps = chrom_ps.setdefault(rec.chrom, rec.pos)
                            call["GT"] = (gt[1], gt[0]) if b["flipped"] else gt
                            call.phased = True
                            call["PS"] = cps
                            n_flip += b["flipped"]
                        elif b["part"] > 0:
                            call["PS"] = b["ps"]          # unjoined part of a split block: its own PS
                    else:
                        p = pop.get((rec.chrom, rec.pos, rec.ref, rec.alts[0])) if rec.alts else None
                        if p is not None:
                            cps = chrom_ps.setdefault(rec.chrom, rec.pos)
                            call["GT"] = p
                            call.phased = True
                            call["PS"] = cps
                            n_pop_only += 1
                vo.write(rec)

    nj = sum(b["joined"] for _, pl in parts.values() for b in pl)
    print(f"read blocks: {len(blocks)} | split: {n_split} (population error {beta}) | parts joined: {nj} | "
          f"flipped SNPs: {n_flip} | population-only SNPs: {n_pop_only}", file=sys.stderr)


if __name__ == "__main__":
    main()
