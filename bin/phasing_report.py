#!/usr/bin/env python3
"""
Truth-free phasing report: figures (PNG), tables (TSV) and one self-contained HTML page.

Reads everything the pipeline drops into --data:
  run_info.json                       run parameters
  het_funnel.tsv, het_calls.tsv.gz, het_kept.tsv.gz      variant calling / het SNP filter
  <id>.phase_stats.tsv, <id>.blocklist.tsv               per phased VCF (id = <read>_<pop> or <read>_reads)
  orientation_<combo>.tsv             read block vs population phase (hybrid_merge.py)
  junctions_<combo>.tsv               masking test (junction_test.py)
  <combo>.reads_vs_final.tsv/.bed     whatshap compare reads-only vs final
  method_concordance.tsv              whatshap compare between all final VCFs
  haplotag_counts.tsv                 reads per chromosome / bin / haplotype
  single_cell_allele_qc.tsv           optional, cellsnp-lite per-cell QC

Every section is skipped when its input is missing.
"""
import argparse
import base64
import glob
import html
import json
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# ---- palette (dataviz reference palette, light mode) -------------------------------------------
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#9a9993"
INK, INK2, MUTED, GRID, SURFACE, TRACK = "#0b0b0b", "#52514e", "#8a8984", "#e4e3de", "#ffffff", "#efeeea"
SEQ = {"hi": "#184f95", "mid": "#5598e7", "lo": "#b7d3f6"}       # sequential blue, concordance bins
NOTJOINED = "#c9c8c1"
H1, H2, UNTAGGED = SERIES[0], SERIES[1], "#c9c8c1"

plt.rcParams.update({
    "font.family": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 9,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 10.5,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.size": 0, "ytick.major.size": 0, "legend.frameon": False, "legend.fontsize": 8.5,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})


def ap():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--out", default=".")
    p.add_argument("--min-qual", type=float, default=20)
    p.add_argument("--min-dp", type=float, default=15)
    p.add_argument("--min-vaf", type=float, default=0.25)
    p.add_argument("--max-vaf", type=float, default=0.75)
    return p.parse_args()


# ---- helpers -------------------------------------------------------------------------------------
def read_tsv(path, **kw):
    try:
        return pd.read_csv(path, sep="\t", **kw) if path and os.path.exists(path) and os.path.getsize(path) else None
    except pd.errors.EmptyDataError:
        return None


def by_suffix(data, suffix, prefix=""):
    out = {}
    for f in sorted(glob.glob(os.path.join(data, f"{prefix}*{suffix}"))):
        out[os.path.basename(f)[len(prefix):-len(suffix)]] = f
    return out


def n50(x):
    x = np.sort(np.asarray(x, dtype=float))[::-1]
    if not len(x):
        return 0.0
    return float(x[np.searchsorted(np.cumsum(x), x.sum() / 2)])


def fmt_bp(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    x = float(x)
    if x >= 1e6:
        return f"{x / 1e6:.2f} Mb"
    if x >= 1e3:
        return f"{x / 1e3:.1f} kb"
    return f"{x:.0f} bp"


def pretty(id_):
    """whatshap_shapeit4 -> WhatsHap + SHAPEIT4 ; whatshap_reads -> WhatsHap (reads only)"""
    names = {"whatshap": "WhatsHap", "longphase": "LongPhase", "shapeit4": "SHAPEIT4", "beagle": "Beagle"}
    r, p = id_.split("_", 1)
    if p in ("reads", "none"):
        return f"{names.get(r, r)} (reads only)"
    return f"{names.get(r, r)} + {names.get(p, p)}"


class Report:
    def __init__(self, out):
        self.out = out
        self.sections = []
        os.makedirs(os.path.join(out, "figures"), exist_ok=True)
        os.makedirs(os.path.join(out, "tables"), exist_ok=True)

    def fig(self, fig, name):
        path = os.path.join(self.out, "figures", f"{name}.png")
        fig.savefig(path, dpi=170, bbox_inches="tight")
        plt.close(fig)
        return path

    def table(self, df, name):
        df.to_csv(os.path.join(self.out, "tables", f"{name}.tsv"), sep="\t", index=False)

    def add(self, title, intro, items):
        self.sections.append((title, intro, items))


def legend_below(ax, handles=None, ncol=4, y=-0.14, **kw):
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, y), ncol=ncol, **kw)


def kfmt(v, _=None):
    return f"{v / 1e3:g}k" if abs(v) >= 1e3 else f"{v:g}"


# ---- sections ------------------------------------------------------------------------------------
def section_variants(R, d, a, chroms, lens, kpi):
    funnel = read_tsv(os.path.join(d, "het_funnel.tsv"))
    calls = read_tsv(os.path.join(d, "het_calls.tsv.gz"))
    kept = read_tsv(os.path.join(d, "het_kept.tsv.gz"))
    if funnel is None and calls is None:
        return
    items = []
    if funnel is not None:
        f = funnel.set_index("step")["n"]
        kpi.append(("Germline het SNPs kept", f"{int(f.get('het_kept', 0)):,}",
                    f"of {int(f.get('het', 0)):,} PASS het SNPs"))
        labels = {"all_records": "All Clair3 records", "pass": "PASS", "pass_biallelic_snps": "PASS biallelic SNPs",
                  "het": "Heterozygous", "het_kept": "Het after QUAL/DP/VAF filter"}
        fu = f.reindex([k for k in labels if k in f.index])
        fig, ax = plt.subplots(figsize=(7.2, 2.4))
        y = np.arange(len(fu))[::-1]
        ax.barh(y, fu.values, color=SERIES[0], height=0.55)
        for yy, v in zip(y, fu.values):
            ax.text(v, yy, f"  {int(v):,}", va="center", ha="left", color=INK, fontsize=8.5)
        ax.set_yticks(y, [labels[k] for k in fu.index])
        ax.set_xlim(0, fu.max() * 1.22)
        ax.grid(axis="y", visible=False)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1e6:.1f}M" if v >= 1e6 else f"{v:,.0f}"))
        ax.set_title("From variant calls to the SNPs used for phasing")
        items.append(("fig", R.fig(fig, "01_het_funnel"), "Counts at each step of the germline het SNP filter."))
        R.table(funnel, "01_het_funnel")

    if calls is not None and len(calls):
        calls["vaf"] = calls.ad_alt / (calls.ad_ref + calls.ad_alt).replace(0, np.nan)
        if kept is not None:
            k = set(zip(kept.chrom, kept.pos))
            calls["kept"] = [(c, p) in k for c, p in zip(calls.chrom, calls.pos)]
        else:
            calls["kept"] = ((calls.qual >= a.min_qual) & (calls.dp >= a.min_dp)
                             & calls.vaf.between(a.min_vaf, a.max_vaf))
        fig, axes = plt.subplots(1, 3, figsize=(11, 2.9))
        specs = [("vaf", "VAF (alt / (ref + alt))", np.linspace(0, 1, 51), [a.min_vaf, a.max_vaf]),
                 ("dp", "Depth (DP)", np.arange(0, np.nanpercentile(calls.dp, 99.5) + 2, 1), [a.min_dp]),
                 ("qual", "QUAL", np.linspace(0, np.nanpercentile(calls.qual, 99.5), 60), [a.min_qual])]
        for ax, (col, lab, bins, thr) in zip(axes, specs):
            ax.hist(calls[col].dropna(), bins=bins, color=SERIES[0], alpha=0.25, label="all PASS het SNPs")
            ax.hist(calls.loc[calls.kept, col].dropna(), bins=bins, histtype="step", color=SERIES[0], lw=1.6,
                    label="kept")
            for t in thr:
                ax.axvline(t, color=INK2, lw=1, ls=(0, (3, 2)))
            ax.set_xlabel(lab)
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(kfmt))
        axes[0].set_ylabel("het SNPs")
        axes[0].set_title("Filter thresholds (dashed)")
        legend_below(axes[0], ncol=2)
        fig.tight_layout()
        items.append(("fig", R.fig(fig, "02_het_filter_distributions"),
                      "VAF, depth and QUAL of all PASS heterozygous SNPs; the outline shows the SNPs that pass. "
                      "A germline sample should give one VAF peak at 0.5."))

        kc = calls[calls.kept].groupby("chrom").size().reindex(chroms).fillna(0)
        dens = kc / (lens.reindex(chroms) / 1e6)
        fig, ax = plt.subplots(figsize=(11, 2.6))
        x = np.arange(len(chroms))
        ax.bar(x, dens.values, color=SERIES[0], width=0.62)
        ax.axhline(np.nanmedian(dens.values), color=INK2, lw=1, ls=(0, (3, 2)))
        ax.text(len(chroms) - 0.4, np.nanmedian(dens.values), " median", va="bottom", ha="right", color=INK2, fontsize=8)
        ax.set_xticks(x, [c.replace("chr", "") for c in chroms])
        ax.grid(axis="x", visible=False)
        ax.set_ylabel("het SNPs per Mb")
        ax.set_title("Het SNP density per chromosome")
        items.append(("fig", R.fig(fig, "03_het_density"),
                      "Kept het SNPs per Mb of chromosome. Low values flag LOH or poorly covered chromosomes."))
        R.table(pd.DataFrame({"chrom": chroms, "het_snps_kept": kc.values.astype(int),
                              "per_mb": dens.values.round(2)}), "03_het_density")
    R.add("1 · Variant calling and het SNP selection",
          "Clair3 calls, filtered to biallelic PASS heterozygous SNPs with the QUAL, depth and VAF thresholds of the run.",
          items)


def section_completeness(R, d, chroms, lens, ids, colors, primary, kpi):
    stats = {i: read_tsv(f) for i, f in by_suffix(d, ".phase_stats.tsv").items()}
    stats = {i: s for i, s in stats.items() if s is not None}
    blocks = {i: read_tsv(f) for i, f in by_suffix(d, ".blocklist.tsv").items()}
    if not stats:
        return
    items = []
    order = [i for i in ids if i in stats]
    rows = []
    for i in order:
        s = stats[i].set_index("chrom")
        al = s.loc["ALL"] if "ALL" in s.index else None
        if al is None:
            continue
        rows.append(dict(method=pretty(i), id=i, het_snvs=int(al.het_snvs), phased=int(al.phased),
                         phased_pct=round(float(al.phased_pct), 2), blocks=int(al.blocks),
                         block_n50=fmt_bp(al.block_n50_bp), largest_block=fmt_bp(al.largest_block_bp)))
    summ = pd.DataFrame(rows)
    R.table(summ, "04_phasing_summary")
    items.append(("table", summ.drop(columns=["id"]), "Genome-wide phasing statistics for each phased VCF."))
    for i in order:
        if i.endswith("_reads"):
            al = stats[i].set_index("chrom").loc["ALL"]
            kpi.append((f"Read-block N50 · {pretty(i).split(' ')[0]}", fmt_bp(al.block_n50_bp),
                        f"{float(al.phased_pct):.1f}% of het SNPs phased by reads"))
    if primary in stats:
        al = stats[primary].set_index("chrom").loc["ALL"]
        kpi.append((f"Phased · {pretty(primary)}", f"{float(al.phased_pct):.1f}%",
                    f"block N50 {fmt_bp(al.block_n50_bp)}"))

    # phased % per chromosome
    fig, ax = plt.subplots(figsize=(11, 3.0))
    x = np.arange(len(chroms))
    for i in order:
        s = stats[i].set_index("chrom").reindex(chroms)
        ls = (0, (3, 2)) if i.endswith("_reads") else "-"
        ax.plot(x, s.phased_pct, color=colors[i], ls=ls, marker="o", ms=4, lw=1.6, label=pretty(i))
    ax.set_xticks(x, [c.replace("chr", "") for c in chroms])
    ax.set_ylabel("% het SNPs phased")
    ax.grid(axis="x", visible=False)
    ax.set_title("Phasing completeness per chromosome")
    legend_below(ax, ncol=3)
    items.append(("fig", R.fig(fig, "05_phased_fraction"),
                  "Dashed = read-only blocks, solid = after the population step. "
                  "Read-only SNPs that fall in no block (singletons) are left unphased."))

    # chromosome-scale N50 per chromosome
    fig, ax = plt.subplots(figsize=(11, 3.0))
    for i in order:
        s = stats[i].set_index("chrom").reindex(chroms)
        ls = (0, (3, 2)) if i.endswith("_reads") else "-"
        ax.plot(x, s.block_n50_bp / 1e6, color=colors[i], ls=ls, marker="o", ms=4, lw=1.6, label=pretty(i))
    ax.plot(x, lens.reindex(chroms) / 1e6, color=MUTED, lw=1, ls=":", label="chromosome length")
    ax.set_yscale("log")
    ax.set_xticks(x, [c.replace("chr", "") for c in chroms])
    ax.set_ylabel("block N50 (Mb, log)")
    ax.grid(axis="x", visible=False)
    ax.set_title("Phase block N50 per chromosome")
    legend_below(ax, ncol=3)
    items.append(("fig", R.fig(fig, "06_block_n50"),
                  "Read blocks are limited by read length and het SNP spacing; the population step should lift "
                  "N50 close to the chromosome length (dotted)."))

    # read-block length distribution
    rb = {i: blocks[i] for i in order if i.endswith("_reads") and blocks.get(i) is not None}
    if rb:
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.0))
        allspan = np.concatenate([b.span_bp.values for b in rb.values()])
        bins = np.logspace(np.log10(max(allspan.min(), 10)), np.log10(allspan.max() + 1), 50)
        for i, b in rb.items():
            b = b[b.n_snvs >= 2]
            axes[0].hist(b.span_bp, bins=bins, histtype="step", lw=1.8, color=colors[i], label=pretty(i))
            nn = n50(b.span_bp)
            axes[0].axvline(nn, color=colors[i], lw=1, ls=(0, (3, 2)))
            srt = b.sort_values("span_bp", ascending=False)
            cum = srt.n_snvs.cumsum() / srt.n_snvs.sum() * 100
            axes[1].plot(srt.span_bp, cum, color=colors[i], label=pretty(i))
        axes[0].set_xscale("log")
        axes[0].set_xlabel("read block span (bp, log)")
        axes[0].set_ylabel("blocks")
        axes[0].set_title("Read block lengths (dashed = N50)")
        axes[1].set_xscale("log")
        axes[1].invert_xaxis()
        axes[1].set_xlabel("block span ≥ x (bp, log)")
        axes[1].set_ylabel("% of phased SNPs")
        axes[1].set_title("Phased SNPs in blocks at least this long")
        legend_below(axes[0], ncol=2)
        fig.tight_layout()
        items.append(("fig", R.fig(fig, "07_read_block_lengths"), "Read-only phase blocks (≥ 2 SNPs)."))

        # gaps between consecutive read blocks = what the population phaser has to bridge
        gaps = []
        for i, b in rb.items():
            b = b[b.n_snvs >= 2].sort_values(["chrom", "start"])
            g = (b.start.shift(-1) - b.end)[b.chrom == b.chrom.shift(-1)]
            gaps.append(pd.DataFrame({"id": i, "gap_bp": g.dropna().values}))
        R.real_gaps = pd.concat(gaps)
    R.add("2 · Phasing completeness and block lengths",
          "How much of the het SNPs each method phases and how long the phase blocks are, for the read-only "
          "blocks and for every read × population combination.", items)
    return blocks


def ideogram(R, d, chroms, lens, combos, colors):
    ori = {c: read_tsv(f) for c, f in by_suffix(d, ".tsv", "orientation_").items()}
    ori = {c: o for c, o in ori.items() if o is not None and c in combos and not c.endswith("_none")}
    if not ori:
        return []
    beds = {}
    for c in ori:
        b = read_tsv(os.path.join(d, f"{c}.reads_vs_final.bed"), header=None, usecols=[0, 1, 2], names=["chrom", "s", "e"])
        beds[c] = b
    items = []
    n = len(ori)
    xmax = lens.reindex(chroms).max() / 1e6
    fig, axes = plt.subplots(1, n, figsize=(6.0 * n, 0.24 * len(chroms) + 1.6), sharey=True, squeeze=False)
    for ax, (c, o) in zip(axes[0], ori.items()):
        o = o.copy()
        o["conc"] = pd.to_numeric(o.concordance, errors="coerce")
        o["col"] = np.where(o.joined == 1,
                            np.where(o.conc >= 0.9, SEQ["hi"], np.where(o.conc >= 0.7, SEQ["mid"], SEQ["lo"])),
                            NOTJOINED)
        for y, ch in enumerate(chroms):
            ax.barh(y, lens.get(ch, 0) / 1e6, height=0.62, color=TRACK, lw=0)
            sub = o[o.chrom == ch]
            if len(sub):
                ax.barh(np.full(len(sub), y), (sub.span_bp.clip(lower=2e4)) / 1e6, left=sub.start / 1e6,
                        height=0.62, color=sub.col, lw=0)
            b = beds.get(c)
            if b is not None:
                sw = b[b.chrom == ch]
                if len(sw):
                    ax.vlines((sw.s + sw.e) / 2e6, y - 0.5, y - 0.33, color=INK, lw=0.5)
        ax.set_title(pretty(c))
        ax.set_xlim(0, xmax * 1.02)
        ax.set_ylim(len(chroms) - 0.4, -0.9)
        ax.set_xlabel("position (Mb)")
        ax.grid(axis="y", visible=False)
    axes[0][0].set_yticks(range(len(chroms)), [c.replace("chr", "") for c in chroms])
    handles = [Patch(color=SEQ["hi"], label="oriented, ≥ 90% SNPs agree"),
               Patch(color=SEQ["mid"], label="oriented, 70–90%"),
               Patch(color=SEQ["lo"], label="oriented, < 70%"),
               Patch(color=NOTJOINED, label="not oriented (too few panel SNPs / tie)"),
               Patch(color=TRACK, label="no read block")]
    fig.legend(handles=handles, loc="lower left", ncol=5, bbox_to_anchor=(0.01, -0.04), frameon=False)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    items.append(("fig", R.fig(fig, "08_block_orientation_ideogram"),
                  "Each bar is a read block, coloured by how many of its SNPs agree with the population phase in "
                  "the block's best orientation. Black ticks above a chromosome: switches between the read-only "
                  "and the final phasing (should be rare)."))
    return items


def section_orientation(R, d, chroms, lens, combos, colors, kpi):
    items = ideogram(R, d, chroms, lens, combos, colors)
    ori = {c: read_tsv(f) for c, f in by_suffix(d, ".tsv", "orientation_").items()}
    ori = {c: o for c, o in ori.items() if o is not None and c in combos and not c.endswith("_none")}
    if ori:
        bins = [0, 1e3, 5e3, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, np.inf]
        labs = ["<1k", "1–5k", "5–20k", "20–50k", "50–100k", "100–200k", "200–500k", "0.5–1M", ">1M"]
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.1))
        rows = []
        for c, o in ori.items():
            o = o.copy()
            o["conc"] = pd.to_numeric(o.concordance, errors="coerce")
            o["bin"] = pd.cut(o.span_bp, bins, labels=labs)
            g = o.groupby("bin", observed=False)
            joined = g.joined.mean() * 100
            conc = o[o.joined == 1].groupby("bin", observed=False).conc.median() * 100
            x = np.arange(len(labs))
            axes[0].plot(x, joined.values, color=colors[c], marker="o", ms=4, label=pretty(c))
            axes[1].plot(x, conc.values, color=colors[c], marker="o", ms=4, label=pretty(c))
            snp_w = (o.joined * o.n_snps).sum() / o.n_snps.sum() * 100
            rows.append(dict(method=pretty(c), read_blocks=len(o), oriented=int(o.joined.sum()),
                             snps_in_oriented_blocks_pct=round(snp_w, 2),
                             median_concordance_pct=round(o.loc[o.joined == 1, "conc"].median() * 100, 2),
                             blocks_below_70pct=int(((o.joined == 1) & (o.conc < 0.7)).sum())))
        for ax in axes:
            ax.set_xticks(np.arange(len(labs)), labs, rotation=30, ha="right")
            ax.set_xlabel("read block span (bp)")
            ax.grid(axis="x", visible=False)
        axes[0].set_ylabel("% of blocks oriented")
        axes[0].set_title("Blocks the population phase can orient")
        axes[1].set_ylabel("median % SNPs agreeing")
        axes[1].set_title("Agreement inside oriented blocks")
        axes[1].set_ylim(48, 101)
        legend_below(axes[0], ncol=2, y=-0.42)
        fig.tight_layout()
        items.append(("fig", R.fig(fig, "09_orientation_vs_block_size"),
                      "Small blocks share few SNPs with the panel and are harder to orient. Inside an oriented "
                      "block, agreement below 100% means read phase and population phase disagree at some SNPs "
                      "(population switch errors, or read phasing errors)."))
        t = pd.DataFrame(rows)
        R.table(t, "09_block_orientation")
        items.append(("table", t, "Read blocks vs population phase."))

    # reads vs final
    cmp_ = {c: read_tsv(f) for c, f in by_suffix(d, ".reads_vs_final.tsv").items()}
    cmp_ = {c: t for c, t in cmp_.items() if t is not None}
    if cmp_:
        fig, ax = plt.subplots(figsize=(11, 2.8))
        x = np.arange(len(chroms))
        w = 0.8 / len(cmp_)
        rows = []
        for k, (c, t) in enumerate(cmp_.items()):
            t = t.set_index("chromosome").reindex(chroms)
            rate = 100 * t.all_switches / t.all_assessed_pairs
            ax.bar(x - 0.4 + w * (k + 0.5), rate.values, width=w * 0.9, color=colors.get(c, OTHER), label=pretty(c))
            rows.append(dict(method=pretty(c), assessed_pairs=int(t.all_assessed_pairs.sum()),
                             switches=int(t.all_switches.sum()),
                             switch_rate_pct=round(100 * t.all_switches.sum() / t.all_assessed_pairs.sum(), 4)))
        ax.set_xticks(x, [c.replace("chr", "") for c in chroms])
        ax.set_ylabel("switch rate (%)")
        ax.grid(axis="x", visible=False)
        ax.set_title("Final phasing vs read blocks (switches inside read blocks)")
        legend_below(ax, ncol=3)
        items.append(("fig", R.fig(fig, "10_reads_vs_final"),
                      "Switches between the read-only blocks and the final VCF, counted inside read blocks. With "
                      "--hybrid_merge orient this is 0 by construction; with whatshap it shows where the population "
                      "pseudo-read overruled the reads."))
        t = pd.DataFrame(rows)
        R.table(t, "10_reads_vs_final")
        items.append(("table", t, ""))
    if items:
        R.add("3 · Read blocks vs population phase",
              "Without a truth set, the read blocks are the best local truth: they are compared with the "
              "population phase that is used to connect them.", items)


def section_junctions(R, d, combos, colors, kpi):
    jt = {c: read_tsv(f) for c, f in by_suffix(d, ".tsv", "junctions_").items()}
    jt = {c: t for c, t in jt.items() if t is not None and len(t)}
    if not jt:
        return
    items = []
    bins = [0, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, np.inf]
    labs = ["<2k", "2–5k", "5–10k", "10–20k", "20–50k", "50–100k", ">100k"]
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.2))
    x = np.arange(len(labs))
    w = 0.8 / len(jt)
    for k, (c, t) in enumerate(jt.items()):
        res = t[t.status != "not_joined"].copy()
        res["bin"] = pd.cut(res.gap_bp, bins, labels=labs)
        g = res.groupby("bin", observed=False)
        rate = g.status.apply(lambda s: 100 * (s == "wrong").mean() if len(s) else np.nan)
        axes[0].bar(x - 0.4 + w * (k + 0.5), rate.values, width=w * 0.9, color=colors.get(c, OTHER), label=pretty(c))
        nw = int((t.status == "wrong").sum())
        nj = int((t.status != "not_joined").sum())
        rows.append(dict(method=pretty(c), test_junctions=len(t), rejoined=nj,
                         rejoined_pct=round(100 * nj / len(t), 2), wrong=nw,
                         wrong_join_pct=round(100 * nw / nj, 3) if nj else np.nan,
                         median_test_gap=fmt_bp(t.gap_bp.median())))
    axes[0].set_xticks(x, labs, rotation=30, ha="right")
    axes[0].set_xlabel("gap between the two halves (bp)")
    axes[0].set_ylabel("% joined the wrong way")
    axes[0].grid(axis="x", visible=False)
    axes[0].set_title("Masked junctions re-joined wrongly")
    legend_below(axes[0], ncol=2, y=-0.42)

    # test gaps vs real gaps: the test is only representative where they overlap
    tg = pd.concat(jt.values()).gap_bp
    hb = np.logspace(0, 7, 50)
    axes[1].hist(tg, bins=hb, histtype="step", lw=1.8, color=INK2, density=True, label="test junctions")
    if getattr(R, "real_gaps", None) is not None and len(R.real_gaps):
        for i, g in R.real_gaps.groupby("id"):
            axes[1].hist(g.gap_bp.clip(lower=1), bins=hb, histtype="step", lw=1.8, color=colors.get(i, OTHER),
                         density=True, label=f"real gaps · {pretty(i)}")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("gap (bp, log)")
    axes[1].set_ylabel("density")
    axes[1].set_title("Gaps: test junctions vs real block junctions")
    legend_below(axes[1], ncol=1)
    fig.tight_layout()
    items.append(("fig", R.fig(fig, "11_junction_test"),
                  "Every read block with enough SNPs is split in the middle and the population phaser re-joins the "
                  "halves blind; the reads say which way is right. Real junctions span larger gaps than these test "
                  "junctions (right), so read the error rate per gap bin, not only the overall rate."))
    t = pd.DataFrame(rows)
    R.table(t, "11_junction_test")
    items.append(("table", t, ""))
    best = t.sort_values("wrong_join_pct").iloc[0]
    kpi.append(("Wrong joins (masking test)", f"{best.wrong_join_pct:.2f}%", f"best: {best.method}"))
    R.add("4 · Junction test (truth-free benchmark)",
          "How often does the population phaser connect two read blocks the wrong way round? Measured by "
          "hiding links that the reads already know.", items)


def section_concordance(R, d, combos, colors):
    t = read_tsv(os.path.join(d, "method_concordance.tsv"))
    if t is None or not len(t):
        return
    t = t[t.chromosome != "ALL"] if "chromosome" in t else t
    names = [c for c in combos if c in set(t.dataset_name0) | set(t.dataset_name1)]
    m_sw = pd.DataFrame(np.nan, index=names, columns=names)
    m_ham = m_sw.copy()
    for (a, b), g in t.groupby(["dataset_name0", "dataset_name1"]):
        sw = 100 * g.all_switches.sum() / g.all_assessed_pairs.sum()
        ham = 100 * g.blockwise_hamming.sum() / g.covered_variants.sum() if "blockwise_hamming" in g else np.nan
        m_sw.loc[a, b] = m_sw.loc[b, a] = sw
        m_ham.loc[a, b] = m_ham.loc[b, a] = ham
    fig, axes = plt.subplots(1, 2, figsize=(4.2 * 2 + 1.5, 0.55 * len(names) + 2.2))
    for ax, m, title in [(axes[0], m_sw, "Switch rate between methods (%)"),
                         (axes[1], m_ham, "Block-wise Hamming between methods (%)")]:
        v = m.values.astype(float)
        vmax = np.nanmax(v) if np.isfinite(v).any() else 1
        ax.imshow(np.ma.masked_invalid(v), cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
            "b", ["#f2f7fd", SEQ["mid"], SEQ["hi"]]), vmin=0, vmax=vmax or 1)
        for i in range(len(names)):
            for j in range(len(names)):
                if np.isfinite(v[i, j]):
                    ax.text(j, i, f"{v[i, j]:.2f}", ha="center", va="center", fontsize=8,
                            color="white" if v[i, j] > 0.6 * vmax else INK)
        ax.set_xticks(range(len(names)), [pretty(n) for n in names], rotation=35, ha="right")
        ax.set_yticks(range(len(names)), [pretty(n) for n in names])
        ax.grid(False)
        ax.set_title(title)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.tight_layout()
    items = [("fig", R.fig(fig, "12_method_concordance"),
              "Pairwise comparison of the final VCFs. Methods that agree are not necessarily right, but large "
              "Hamming distances show where at least one of them makes long-range errors.")]
    R.table(m_sw.reset_index().rename(columns={"index": "method"}), "12_method_switch_rate")
    R.table(m_ham.reset_index().rename(columns={"index": "method"}), "12_method_hamming")
    R.add("5 · Agreement between methods", "", items)


def section_haplotag(R, d, chroms, lens, kpi, primary, tool, bin_mb=1):
    t = read_tsv(os.path.join(d, "haplotag_counts.tsv"))
    if t is None or not len(t):
        return
    t = t[t.chrom.isin(chroms)]
    items = []
    per = t.pivot_table(index="chrom", columns="hp", values="reads", aggfunc="sum").reindex(chroms).fillna(0)
    for h in (0, 1, 2):
        if h not in per:
            per[h] = 0
    tot = per.sum(axis=1)
    pct = per.div(tot.replace(0, np.nan), axis=0) * 100
    tagged = 100 * (per[1].sum() + per[2].sum()) / tot.sum()
    kpi.append(("Reads haplotagged", f"{tagged:.1f}%", f"H1 : H2 = {per[1].sum() / max(per[2].sum(), 1):.3f}"))
    fig, ax = plt.subplots(figsize=(11, 2.9))
    x = np.arange(len(chroms))
    bottom = np.zeros(len(chroms))
    for h, col, lab in [(1, H1, "H1"), (2, H2, "H2"), (0, UNTAGGED, "untagged")]:
        ax.bar(x, pct[h].values, bottom=bottom, color=col, width=0.62, label=lab, edgecolor=SURFACE, lw=1)
        bottom += pct[h].fillna(0).values
    ax.set_xticks(x, [c.replace("chr", "") for c in chroms])
    ax.set_ylabel("% primary reads")
    ax.set_ylim(0, 100)
    ax.grid(axis="x", visible=False)
    ax.set_title(f"Haplotagged reads per chromosome · {pretty(primary)} · {tool}")
    legend_below(ax, ncol=3)
    items.append(("fig", R.fig(fig, "13_haplotag_per_chrom"), "Primary mapped reads assigned to H1, H2 or neither."))
    out = per.rename(columns={0: "untagged", 1: "H1", 2: "H2"}).reset_index()
    R.table(out, "13_haplotag_per_chrom")

    # H1 fraction along the genome
    b = t[t.hp.isin([1, 2])].pivot_table(index=["chrom", "bin"], columns="hp", values="reads", aggfunc="sum").fillna(0)
    b = b[(b[1] + b[2]) >= 20].reset_index()
    if len(b):
        b["frac"] = b[1] / (b[1] + b[2])
        offs, o = {}, 0
        for c in chroms:
            offs[c] = o
            o += lens.get(c, 0)
        binsize = float(bin_mb) * 1e6
        b["gpos"] = [offs[c] + (k + 0.5) * binsize for c, k in zip(b.chrom, b.bin)]
        fig, ax = plt.subplots(figsize=(11, 2.8))
        for k, c in enumerate(chroms):
            if k % 2:
                ax.axvspan(offs[c] / 1e6, (offs[c] + lens.get(c, 0)) / 1e6, color=TRACK, lw=0)
        ax.scatter(b.gpos / 1e6, b.frac, s=3, color=SERIES[0], lw=0)
        ax.axhline(0.5, color=INK2, lw=1, ls=(0, (3, 2)))
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlim(0, o / 1e6)
        ax.set_xticks([(offs[c] + lens.get(c, 0) / 2) / 1e6 for c in chroms], [c.replace("chr", "") for c in chroms])
        ax.grid(False)
        ax.set_ylabel("H1 / (H1 + H2)")
        ax.set_title(f"Haplotype 1 read fraction along the genome ({binsize / 1e6:g} Mb bins)")
        items.append(("fig", R.fig(fig, "14_haplotype_fraction_genome"),
                      "A diploid normal sits at 0.5. Copy-number changes and LOH shift whole segments away from 0.5 "
                      "(mirrored if the phase is correct over long distances); jumps between mirrored levels inside a "
                      "segment are long-range switch errors."))
        R.table(b[["chrom", "bin", 1, 2, "frac"]].rename(columns={1: "H1_reads", 2: "H2_reads", "frac": "H1_fraction"}),
                "14_haplotype_fraction_bins")
    R.add("6 · Haplotagged BAM", "Reads assigned to haplotypes with the final phased VCF.", items)


def section_cells(R, d):
    q = read_tsv(os.path.join(d, "single_cell_allele_qc.tsv"))
    if q is None or not len(q):
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 2.9))
    axes[0].hist(q.covered_snps.clip(lower=1), bins=np.logspace(0, np.log10(q.covered_snps.max() + 1), 50),
                 color=SERIES[0])
    axes[0].set_xscale("log")
    axes[0].set_xlabel("covered phased SNPs per cell (log)")
    axes[0].set_ylabel("cells")
    axes[0].set_title(f"Phased SNPs covered per cell (n = {len(q):,} cells)")
    axes[1].scatter(q.total_depth.clip(lower=1), q.covered_snps.clip(lower=1), s=4, color=SERIES[0], lw=0)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("total allele depth")
    axes[1].set_ylabel("covered SNPs")
    axes[1].set_title("Depth vs covered SNPs")
    fig.tight_layout()
    desc = q[["covered_snps", "total_depth", "median_cov_per_covered_snp"]].describe().T.reset_index()
    R.add("7 · Single-cell allele counts", "cellsnp-lite on the final phased SNPs (input for ASCAT.sc).",
          [("fig", R.fig(fig, "15_cell_qc"), ""), ("table", desc.round(1), "")])


# ---- HTML ----------------------------------------------------------------------------------------
CSS = """
:root{--bg:#f7f7f5;--card:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--muted:#8a8984;--line:#e4e3de;--accent:#2a78d6}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#121211;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#8f8e88;--line:#2e2d2a;--accent:#3987e5}}
:root[data-theme="dark"]{--bg:#121211;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#8f8e88;--line:#2e2d2a;--accent:#3987e5}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:1180px;margin:0 auto;padding:28px 16px 60px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:18px;margin:40px 0 6px;padding-top:8px;border-top:1px solid var(--line)}
.sub{color:var(--ink2);margin:0 0 18px}.intro{color:var(--ink2);margin:0 0 14px;max-width:820px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px;margin:18px 0}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.kpi .l{font-size:12px;color:var(--ink2)}.kpi .v{font-size:26px;font-weight:650;margin:2px 0;font-variant-numeric:tabular-nums}.kpi .n{font-size:12px;color:var(--muted)}
figure{margin:14px 0 22px;background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px}
figure img{width:100%;height:auto;display:block}figcaption{color:#52514e;font-size:12.5px;margin-top:8px;max-width:900px}
.tw{overflow-x:auto;margin:10px 0 22px;background:var(--card);border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:12.5px;font-variant-numeric:tabular-nums}
th,td{padding:7px 12px;text-align:right;border-bottom:1px solid var(--line);white-space:nowrap}
th:first-child,td:first-child{text-align:left}th{color:var(--ink2);font-weight:600}tr:last-child td{border-bottom:none}
.cap{color:var(--ink2);font-size:12.5px;padding:8px 12px 0}
nav{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:13px;margin:6px 0 0}nav a{color:var(--accent);text-decoration:none}
details{margin-top:40px;color:var(--ink2)}pre{white-space:pre-wrap;font-size:12px;background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px}
"""


def img64(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def write_html(R, info, kpi, path):
    title = f"Phasing report · {html.escape(str(info.get('sample', '')))}"
    parts = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
             f"<meta name='viewport' content='width=device-width,initial-scale=1'><title>{title}</title>"
             f"<style>{CSS}</style></head><body><main>",
             f"<h1>{title}</h1>",
             f"<p class='sub'>{html.escape(str(info.get('variant_calls', '')))} · "
             f"{len(info.get('chroms', []))} chromosomes · merge: {html.escape(str(info.get('hybrid_merge', '')))} · "
             f"run {html.escape(str(info.get('run_name', '')))} ({html.escape(str(info.get('started', ''))[:16])})</p>"]
    parts.append("<nav>" + "".join(f"<a href='#s{i}'>{html.escape(t)}</a>" for i, (t, _, _) in enumerate(R.sections)) + "</nav>")
    if kpi:
        parts.append("<div class='kpis'>" + "".join(
            f"<div class='kpi'><div class='l'>{html.escape(l)}</div><div class='v'>{html.escape(v)}</div>"
            f"<div class='n'>{html.escape(n)}</div></div>" for l, v, n in kpi) + "</div>")
    for i, (t, intro, items) in enumerate(R.sections):
        parts.append(f"<h2 id='s{i}'>{html.escape(t)}</h2>")
        if intro:
            parts.append(f"<p class='intro'>{html.escape(intro)}</p>")
        for kind, obj, cap in items:
            if kind == "fig":
                parts.append(f"<figure><img alt='{html.escape(cap or os.path.basename(obj))}' "
                             f"src='data:image/png;base64,{img64(obj)}'>"
                             + (f"<figcaption>{html.escape(cap)}</figcaption>" if cap else "") + "</figure>")
            else:
                parts.append("<div class='tw'>" + (f"<div class='cap'>{html.escape(cap)}</div>" if cap else "")
                             + obj.to_html(index=False, border=0, na_rep="–", escape=True) + "</div>")
    parts.append("<details><summary>Run parameters</summary><pre>" + html.escape(json.dumps(info, indent=1)) + "</pre></details>")
    parts.append("</main></body></html>")
    with open(path, "w") as fh:
        fh.write("\n".join(parts))


def main():
    a = ap()
    d = a.data
    info = {}
    if os.path.exists(os.path.join(d, "run_info.json")):
        info = json.load(open(os.path.join(d, "run_info.json")))
    lens = pd.read_csv(a.fai, sep="\t", header=None, usecols=[0, 1], names=["c", "l"]).set_index("c").l
    chroms = info.get("chroms") or [c for c in lens.index if re.fullmatch(r"(chr)?\d+", c)]
    combos = info.get("combos") or sorted(by_suffix(d, ".tsv", "orientation_"))
    reads_ids = sorted({i for i in by_suffix(d, ".phase_stats.tsv") if i.endswith("_reads")})
    ids = combos + reads_ids
    # colour follows the entity: fixed slot per method, never re-assigned
    colors = {i: (SERIES[k] if k < len(SERIES) else OTHER) for k, i in enumerate(ids)}
    primary = info.get("haplotag_combo") or (combos[0] if combos else None)

    R = Report(a.out)
    R.real_gaps = None
    kpi = []
    section_variants(R, d, a, chroms, lens, kpi)
    section_completeness(R, d, chroms, lens, ids, colors, primary, kpi)
    section_orientation(R, d, chroms, lens, combos, colors, kpi)
    section_junctions(R, d, combos, colors, kpi)
    section_concordance(R, d, combos, colors)
    section_haplotag(R, d, chroms, lens, kpi, primary, info.get("haplotag_tool", "whatshap"), info.get("haplotag_bin_mb", 1))
    section_cells(R, d)
    write_html(R, info, kpi, os.path.join(a.out, "phasing_report.html"))
    print("wrote", os.path.join(a.out, "phasing_report.html"))


if __name__ == "__main__":
    main()
