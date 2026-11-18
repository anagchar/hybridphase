#!/usr/bin/env python3
"""
Haplotag reads by competitive mapping to a haplotype-resolved diploid assembly.

The reads are aligned to hap1 + hap2 together (contigs prefixed h1_ / h2_, minimap2 with secondary
alignments). A read goes to the haplotype with the better alignment score (AS) if the difference is at
least --min-diff; reads that align equally well to both (homozygous regions, short reads without a
het site) stay unassigned.

  assign:  assembly_haplotag.py assign asm.bam assignments.tsv.gz summary_prefix [--min-diff 6]
           asm.bam must be grouped by read (minimap2 output order, NOT coordinate-sorted)
  tag:     assembly_haplotag.py tag assignments.tsv.gz in.bam out.bam summary_prefix [--bin 1000000]
           writes HP (1/2), ha (assembly contig) and hd (AS difference) on every record of an assigned
           read. If in.bam is already haplotagged from the phased VCF (HP/PS), those tags are read for the
           agreement table and then replaced.
"""
import argparse
import gzip
import sys
from collections import defaultdict

import pysam

LEN_BINS = [0, 500, 1000, 2000, 5000, 10000, 20000, 50000, float("inf")]
LEN_LABS = ["<500", "500–1k", "1–2k", "2–5k", "5–10k", "10–20k", "20–50k", ">50k"]


def len_bin(n):
    for k in range(len(LEN_BINS) - 1):
        if n < LEN_BINS[k + 1]:
            return LEN_LABS[k]
    return LEN_LABS[-1]


def hap_of(contig):
    if contig.startswith("h1_"):
        return 1
    if contig.startswith("h2_"):
        return 2
    return None


def assign(a):
    out = gzip.open(a.out, "wt")
    out.write("read\thp\tas_h1\tas_h2\tcontig\tread_len\n")
    by_len = defaultdict(lambda: defaultdict(int))
    diff_hist = defaultdict(int)

    def flush(name, recs):
        best = {1: None, 2: None}
        contig = {1: ".", 2: "."}
        qlen = 0
        for r in recs:
            if not r.is_secondary and not r.is_supplementary:
                qlen = r.infer_read_length() or r.query_length or 0
            if r.is_unmapped or r.is_supplementary:
                continue
            h = hap_of(r.reference_name)
            if h is None or not r.has_tag("AS"):
                continue
            s = r.get_tag("AS")
            if best[h] is None or s > best[h]:
                best[h], contig[h] = s, r.reference_name
        if best[1] is None and best[2] is None:
            hp, d = "unmapped", None
        else:
            s1 = best[1] if best[1] is not None else -10 ** 9
            s2 = best[2] if best[2] is not None else -10 ** 9
            d = s1 - s2
            hp = "H1" if d >= a.min_diff else "H2" if -d >= a.min_diff else "ambiguous"
            if best[1] is not None and best[2] is not None:
                diff_hist[max(-200, min(200, d))] += 1
        c = contig[1] if hp == "H1" else contig[2] if hp == "H2" else (contig[1] if contig[1] != "." else contig[2])
        out.write(f"{name}\t{hp}\t{'' if best[1] is None else best[1]}\t{'' if best[2] is None else best[2]}\t{c}\t{qlen}\n")
        by_len[len_bin(qlen)][hp] += 1

    with pysam.AlignmentFile(a.bam, check_sq=False) as f:
        name, recs = None, []
        for r in f:
            if r.query_name != name:
                if name is not None:
                    flush(name, recs)
                name, recs = r.query_name, []
            recs.append(r)
        if name is not None:
            flush(name, recs)
    out.close()

    with open(a.prefix + ".by_length.tsv", "w") as fh:
        fh.write("length_bin\tclass\treads\n")
        for lb in LEN_LABS:
            for k, v in by_len[lb].items():
                fh.write(f"{lb}\t{k}\t{v}\n")
    with open(a.prefix + ".as_diff.tsv", "w") as fh:
        fh.write("as_diff\treads\n")
        for k in sorted(diff_hist):
            fh.write(f"{k}\t{diff_hist[k]}\n")
    tot = sum(sum(v.values()) for v in by_len.values())
    n = defaultdict(int)
    for v in by_len.values():
        for k, x in v.items():
            n[k] += x
    print(f"reads: {tot:,} | H1 {n['H1']:,} | H2 {n['H2']:,} | ambiguous {n['ambiguous']:,} | "
          f"unmapped {n['unmapped']:,}", file=sys.stderr)


def tag(a):
    hp = {}
    info = {}
    with gzip.open(a.assign, "rt") as fh:
        next(fh)
        for line in fh:
            name, h, s1, s2, contig, _ = line.rstrip("\n").split("\t")
            if h in ("H1", "H2"):
                hp[name] = 1 if h == "H1" else 2
                info[name] = (contig, int(s1 or 0) - int(s2 or 0) if s1 and s2 else 0)

    bins = defaultdict(int)                    # (chrom, bin, hp) -> primary reads
    agree = defaultdict(lambda: [0, 0])        # (chrom, PS, bin) -> [same, opposite] vs VCF-based HP
    cells = defaultdict(lambda: [0, 0, 0])     # (CB, chrom) -> [H1, H2, unassigned]
    with pysam.AlignmentFile(a.bam) as f:
        header = f.header.to_dict()
        header.setdefault("PG", []).append({"ID": "assembly_haplotag", "PN": "assembly_haplotag.py",
                                            "CL": " ".join(sys.argv)})
        with pysam.AlignmentFile(a.out, "wb", header=header) as o:
            for r in f:
                vcf_hp = r.get_tag("HP") if r.has_tag("HP") else None
                vcf_ps = r.get_tag("PS") if r.has_tag("PS") else None
                for t in ("HP", "PS", "ha", "hd"):
                    if r.has_tag(t):
                        r.set_tag(t, None)
                h = hp.get(r.query_name)
                if h is not None:
                    r.set_tag("HP", h, "i")
                    r.set_tag("ha", info[r.query_name][0], "Z")
                    r.set_tag("hd", info[r.query_name][1], "i")
                if not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate):
                    b = r.reference_start // a.bin
                    bins[(r.reference_name, b, h or 0)] += 1
                    if h is not None and vcf_hp in (1, 2):
                        agree[(r.reference_name, vcf_ps, b)][0 if vcf_hp == h else 1] += 1
                    if r.has_tag("CB"):
                        cells[(r.get_tag("CB"), r.reference_name)][(h or 3) - 1] += 1
                o.write(r)

    with open(a.prefix + ".bins.tsv", "w") as fh:
        fh.write("chrom\tbin\thp\treads\n")
        for (c, b, h), n in sorted(bins.items(), key=lambda x: (str(x[0][0]), x[0][1], x[0][2])):
            fh.write(f"{c}\t{b}\t{h}\t{n}\n")
    with open(a.prefix + ".vs_vcf_haplotag.tsv", "w") as fh:
        fh.write("chrom\tps\tbin\tsame\topposite\n")
        for (c, ps, b), (s, op) in sorted(agree.items(), key=lambda x: (str(x[0][0]), str(x[0][1]), x[0][2])):
            fh.write(f"{c}\t{ps}\t{b}\t{s}\t{op}\n")
    if cells:
        with gzip.open(a.prefix + ".per_cell.tsv.gz", "wt") as fh:
            fh.write("cell\tchrom\tH1\tH2\tunassigned\n")
            for (cb, c), (x1, x2, x0) in sorted(cells.items(), key=lambda x: (x[0][0], str(x[0][1]))):
                fh.write(f"{cb}\t{c}\t{x1}\t{x2}\t{x0}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("assign")
    s.add_argument("bam")
    s.add_argument("out")
    s.add_argument("prefix")
    s.add_argument("--min-diff", type=int, default=6, help="min. AS difference to assign a haplotype [6]")
    t = sub.add_parser("tag")
    t.add_argument("assign")
    t.add_argument("bam")
    t.add_argument("out")
    t.add_argument("prefix")
    t.add_argument("--bin", type=int, default=1_000_000)
    a = ap.parse_args()
    assign(a) if a.cmd == "assign" else tag(a)


if __name__ == "__main__":
    main()
