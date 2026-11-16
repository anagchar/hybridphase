#!/usr/bin/env python3
"""
Draws the metro / tube map of the pipeline (assets/splongget_tube_map.svg + .png).

Edit STATIONS / LINES below when the pipeline changes and re-run:
    python3 assets/make_tube_map.py
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch, Polygon
from matplotlib.path import Path

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "splongget_tube_map")

# ---- colours: one per line, tool names in TOOL ---------------------------------------------------
SHARED = "#1B2944"      # path every run takes
WHATSHAP = "#2a78d6"
LONGPHASE = "#1baf7a"
SHAPEIT4 = "#eb6834"
BEAGLE = "#eda100"
QC = "#6250d6"
CELLS = "#e87ba4"
TOOL = "#1c7c54"
INK = "#0b0b0b"

W, H = 1850, 1080
LW = 13                  # track width (pt)


def rounded(points, r=46):
    """polyline -> Path with rounded corners (quadratic Bezier at every inner vertex)"""
    verts, codes = [points[0]], [Path.MOVETO]
    for (x0, y0), (x1, y1), (x2, y2) in zip(points, points[1:], points[2:]):
        d1 = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5
        d2 = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
        k1, k2 = min(r, d1 / 2) / d1, min(r, d2 / 2) / d2
        a = (x1 - (x1 - x0) * k1, y1 - (y1 - y0) * k1)
        b = (x1 + (x2 - x1) * k2, y1 + (y2 - y1) * k2)
        verts += [a, (x1, y1), b]
        codes += [Path.LINETO, Path.CURVE3, Path.CURVE3]
    verts.append(points[-1])
    codes.append(Path.LINETO)
    return Path(verts, codes)


def track(ax, pts, color, dashed=False, z=2):
    ax.add_patch(PathPatch(rounded(pts), fill=False, edgecolor=color, lw=LW, capstyle="butt",
                           joinstyle="round", zorder=z, linestyle=(0, (1.4, 0.9)) if dashed else "-"))


def file_icon(ax, x, y, label, w=74, h=84):
    """document with a folded corner and a black label band, centred on (x, y)"""
    l, t, f = x - w / 2, y - h / 2, 18
    ax.add_patch(Polygon([(l, t), (l + w - f, t), (l + w, t + f), (l + w, t + h), (l, t + h)], closed=True,
                         facecolor="white", edgecolor=INK, lw=3, zorder=12, joinstyle="round"))
    ax.add_patch(Polygon([(l + w - f, t), (l + w - f, t + f), (l + w, t + f)], closed=True,
                         facecolor="white", edgecolor=INK, lw=3, zorder=13, joinstyle="round"))
    ax.add_patch(FancyBboxPatch((l - 8, y - 15), w + 16, 30, boxstyle="round,pad=0,rounding_size=5",
                                facecolor=INK, edgecolor=INK, zorder=14))
    ax.text(x, y, label, ha="center", va="center", color="white", fontsize=13, fontweight="bold", zorder=15)


def label(ax, x, y, text, tools="", va="bottom", ha="center"):
    """station name (black) + tools (green) as one block anchored at (x, y)"""
    lines = text.split("\n") + (tools.split("\n") if tools else [])
    n = len(text.split("\n")) if text else 0
    if not text:
        lines = tools.split("\n")
    lh = 21
    y0 = y - (len(lines) - 1) * lh if va == "bottom" else y - (len(lines) - 1) * lh / 2 if va == "center" else y
    for i, t in enumerate(lines):
        ax.text(x, y0 + i * lh, t, ha=ha, va="center", fontsize=12.5, fontweight="bold", zorder=11,
                color=INK if i < n else TOOL)


def stop(ax, x, y):
    ax.add_patch(Circle((x, y), 15, facecolor="white", edgecolor=INK, lw=4.5, zorder=10))


def end_bar(ax, x, y, vertical):
    if vertical:
        ax.plot([x - 26, x + 26], [y, y], color=INK, lw=6, solid_capstyle="butt", zorder=3)
    else:
        ax.plot([x, x], [y - 26, y + 26], color=INK, lw=6, solid_capstyle="butt", zorder=3)


def main():
    plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"]})
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")

    Y1, Y2 = 380, 820        # main rows
    UP, DN = 280, 480        # the two method tracks around row 1
    XR = 1500                # right-hand turn
    QY = 650                 # QC line

    # ---------------- tracks
    track(ax, [(110, Y1), (600, Y1)], SHARED)
    track(ax, [(600, Y1), (655, UP), (865, UP), (920, Y1)], WHATSHAP)
    track(ax, [(600, Y1), (655, DN), (865, DN), (920, Y1)], LONGPHASE)
    track(ax, [(920, Y1), (975, UP), (1185, UP), (1240, Y1)], SHAPEIT4)
    track(ax, [(920, Y1), (975, DN), (1185, DN), (1240, Y1)], BEAGLE)
    track(ax, [(1240, Y1), (XR, Y1), (XR, Y2), (430, Y2)], SHARED)
    track(ax, [(XR, 600), (XR - 60, QY), (760, QY), (700, Y2)], QC, z=1)
    track(ax, [(XR, 690), (1660, 690), (1660, 1000), (1580, 1000)], CELLS, dashed=True, z=1)
    track(ax, [(345, Y1), (345, 245)], SHARED, z=1)                 # Clair3 VCF out
    track(ax, [(XR, 560), (1580, 560)], SHARED, z=1)                # phased VCF out
    track(ax, [(1330, Y2), (1330, 880)], SHARED, z=1)               # haplotagged BAM out

    # ---------------- row 1
    file_icon(ax, 70, Y1, "BAM")
    label(ax, 70, Y1 + 58, "long-read BAM\n(CB tags)", va="top")
    for x, name, tool in [(210, "Duplicate\nfilter", "samtools"), (345, "Variant\ncalling", "Clair3"),
                          (495, "Germline het\nSNP filter", "bcftools")]:
        stop(ax, x, Y1)
        label(ax, x, Y1 + 32, name, tool, va="top")
    end_bar(ax, 345, 245, True)
    file_icon(ax, 345, 193, "VCF")
    label(ax, 345, 138, "Clair3 variants")

    ax.text(760, UP - 75, "Read-based phasing", ha="center", fontsize=14.5, fontweight="bold", color=INK)
    ax.text(1080, UP - 75, "Population phasing", ha="center", fontsize=14.5, fontweight="bold", color=INK)
    ax.text(1080, UP - 54, "1kGP panel + genetic maps", ha="center", fontsize=11.5, color="#52514e")
    for x, y, tool, va in [(760, UP, "WhatsHap", "bottom"), (760, DN, "LongPhase", "top"),
                           (1080, UP, "SHAPEIT4 --use-PS", "bottom"), (1080, DN, "Beagle 5", "top")]:
        stop(ax, x, y)
        label(ax, x, y - 30 if va == "bottom" else y + 32, "", tool, va=va)

    stop(ax, 1360, Y1)
    label(ax, 1360, Y1 - 30, "Hybrid merge", "block orientation\nor WhatsHap")

    # ---------------- right side: phased VCF + single-cell branch
    stop(ax, XR, 560)
    label(ax, XR - 26, 560, "Chromosome-scale\nphased VCF", "bcftools concat", va="center", ha="right")
    end_bar(ax, 1580, 560, False)
    file_icon(ax, 1632, 560, "VCF")
    for y, name, tool in [(790, "Allele counts\nper cell", "cellsnp-lite"), (900, "ASCAT.sc\ntables", "R")]:
        stop(ax, 1660, y)
        label(ax, 1686, y, name, tool, va="center", ha="left")
    end_bar(ax, 1580, 1000, False)
    file_icon(ax, 1528, 1000, "TSV")
    label(ax, 1528, 938, "per-cell\nallele counts")

    # ---------------- QC line
    for x, name, tool in [(1250, "Junction test", "masked block joins"),
                          (1000, "Phase stats &\nmethod agreement", "WhatsHap compare")]:
        stop(ax, x, QY)
        label(ax, x, QY - 30, name)
        ax.text(x, QY + 32, tool, ha="center", va="center", fontsize=12.5, fontweight="bold", color=TOOL)

    # ---------------- row 2 (right -> left)
    stop(ax, 1330, Y2)
    label(ax, 1330, Y2 - 30, "Haplotag reads", "WhatsHap or LongPhase")
    end_bar(ax, 1330, 880, True)
    file_icon(ax, 1330, 932, "BAM")
    label(ax, 1330, 990, "Haplotagged BAM", va="top")
    stop(ax, 1120, Y2)
    label(ax, 1120, Y2 - 30, "H1 / H2 read\ncounts", "samtools")
    stop(ax, 700, Y2)
    label(ax, 700, Y2 + 32, "Truth-free QC\nreport", "Python / matplotlib", va="top")
    end_bar(ax, 430, Y2, False)
    file_icon(ax, 378, Y2, "HTML")
    label(ax, 378, Y2 + 58, "Phasing report\n+ figures, tables", va="top")

    # ---------------- legend + title
    lx, ly = 70, 520
    ax.text(lx, ly, "Legend", fontsize=20, fontweight="bold", color=SHARED, va="center")
    items = [(SHARED, "Shared path", False), (WHATSHAP, "Read phasing: WhatsHap", False),
             (LONGPHASE, "Read phasing: LongPhase", False), (SHAPEIT4, "Population phasing: SHAPEIT4", False),
             (BEAGLE, "Population phasing: Beagle", False), (QC, "QC and benchmarking", False),
             (CELLS, "Optional: single cell (--barcodes)", True)]
    for i, (c, lab, dsh) in enumerate(items):
        yy = ly + 36 + i * 26
        ax.add_line(Line2D([lx, lx + 60], [yy, yy], color=c, lw=LW * 0.8, solid_capstyle="butt",
                           linestyle=(0, (1.4, 0.9)) if dsh else "-"))
        ax.text(lx + 74, yy, lab, va="center", fontsize=12, fontweight="bold", color=INK)
    ax.text(70, 985, "splongget-phasing", fontsize=34, fontweight="bold", color=SHARED, va="center")
    ax.text(70, 1025, "every read × population combination runs; the report compares them", fontsize=12.5,
            color="#52514e", va="center")

    fig.savefig(OUT + ".svg", facecolor="white")
    fig.savefig(OUT + ".png", facecolor="white", dpi=150)
    print("wrote", OUT + ".svg", OUT + ".png")


if __name__ == "__main__":
    main()
