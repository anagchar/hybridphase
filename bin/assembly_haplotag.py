#!/usr/bin/env python3
"""
Haplotag reads with a haplotype-resolved diploid assembly, from HipHap output.  (UNDER DEVELOPMENT)

HipHap (Heinz, Marin, Meyerson & Li, bioRxiv 2026) compares each read's alignments to hap1 and hap2 and
keeps the better-supported one, tagged HP (1/2) and hq (HapQ, Phred-scaled confidence; 0 = tie).

  table:  assembly_haplotag.py table merged.bam assignments.tsv.gz summary_prefix [--min-hapq 20]
          merged.bam = HipHap merged output (read order as from minimap2)
  tag:    assembly_haplotag.py tag assignments.tsv.gz in.bam out.bam summary_prefix [--min-hapq 20] [--bin 1000000]
          writes HP (1/2), hq (HapQ) and ha (assembly contig) on every record of a read with HapQ >= min-hapq.
          If in.bam is already haplotagged from the phased VCF (HP/PS), those tags are read for the
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


def table(a):
    out = gzip.open(a.out, "wt")
    out.write("read\thp\thapq\tcontig\tread_len\n")
    by_len = defaultdict(lambda: defaultdict(int))
    hapq_hist = defaultdict(int)
    with pysam.AlignmentFile(a.bam, check_sq=False) as f:
        for r in f:
            if r.is_secondary or r.is_supplementary:
                continue
            qlen = r.infer_read_length() or r.query_length or 0
            if r.is_unmapped:
                cls, hp, q, contig = "unmapped", 0, 0, "."
            else:
                hp = r.get_tag("HP") if r.has_tag("HP") else 0
                q = r.get_tag("hq") if r.has_tag("hq") else 0
                contig = r.reference_name
                hapq_hist[min(q, 60)] += 1
                cls = ("H1" if hp == 1 else "H2") if q >= a.min_hapq and hp in (1, 2) else \
                    "tie (HapQ 0)" if q == 0 else f"HapQ < {a.min_hapq}"
            out.write(f"{r.query_name}\t{hp}\t{q}\t{contig}\t{qlen}\n")
            by_len[len_bin(qlen)][cls] += 1
    out.close()

    with open(a.prefix + ".by_length.tsv", "w") as fh:
        fh.write("length_bin\tclass\treads\n")
        for lb in LEN_LABS:
            for k, v in by_len[lb].items():
                fh.write(f"{lb}\t{k}\t{v}\n")
    with open(a.prefix + ".hapq.tsv", "w") as fh:
        fh.write("hapq\treads\n")
        for k in sorted(hapq_hist):
            fh.write(f"{k}\t{hapq_hist[k]}\n")
    n = defaultdict(int)
    for v in by_len.values():
        for k, x in v.items():
            n[k] += x
    print("reads: " + " | ".join(f"{k} {v:,}" for k, v in sorted(n.items())), file=sys.stderr)


def tag(a):
    hp = {}
    info = {}
    with gzip.open(a.assign, "rt") as fh:
        next(fh)
        for line in fh:
            name, h, q, contig, _ = line.rstrip("\n").split("\t")
            if h in ("1", "2") and int(q) >= a.min_hapq:
                hp[name] = int(h)
                info[name] = (contig, int(q))

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
                for t in ("HP", "PS", "ha", "hq"):
                    if r.has_tag(t):
                        r.set_tag(t, None)
                h = hp.get(r.query_name)
                if h is not None:
                    r.set_tag("HP", h, "i")
                    r.set_tag("ha", info[r.query_name][0], "Z")
                    r.set_tag("hq", info[r.query_name][1], "i")
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
    s = sub.add_parser("table")
    s.add_argument("bam")
    s.add_argument("out")
    s.add_argument("prefix")
    s.add_argument("--min-hapq", type=int, default=20, help="min. HapQ to count a read as assigned [20]")
    t = sub.add_parser("tag")
    t.add_argument("assign")
    t.add_argument("bam")
    t.add_argument("out")
    t.add_argument("prefix")
    t.add_argument("--min-hapq", type=int, default=20)
    t.add_argument("--bin", type=int, default=1_000_000)
    a = ap.parse_args()
    table(a) if a.cmd == "table" else tag(a)


if __name__ == "__main__":
    main()
