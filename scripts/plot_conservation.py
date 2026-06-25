#!/usr/bin/env python
"""4-panel conservation figure, modeled on the SARS-CoV-2 reference paper.

Layout (top to bottom):
    A) Genome map: rectangles for each gene/ORF, labeled with the gene name.
    B) Per-nucleotide identity scatter plot (0.90 - 1.00 y-axis, 0 - 30 kb x-axis).
    C) Close-up of a chosen conserved region (auto-pick the largest), with
       primer annotations and a sequence logo / alignment detail.
    D) Multi-sequence alignment detail (subset of sequences and the conserved
       region columns). Dots for identity, red letters for mismatches.

The function is self-contained and can be called with only the
identity vector and a list of conserved regions; if ``msa_records`` is
provided, panel D uses the MSA. Otherwise, panel D uses a simple
pairwise-stacked representation from ``records``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

from Bio import SeqIO
from Bio.Seq import Seq


# ---------------------------------------------------------------------------
# color palette (matched to the reference figure as closely as possible)
# ---------------------------------------------------------------------------
COLORS = {
    "A": "#56B4E9",   # sky blue (ORF1a)
    "B": "#56B4E9",   # sky blue (ORF1b)
    "S": "#F08080",   # light coral / pink
    "3": "#FFD580",   # peach (3a)
    "E": "#90EE90",   # light green
    "M": "#C8A2C8",   # light purple
    "6": "#A9A9A9",   # gray
    "7": "#FFA07A",   # salmon (7a)
    "8": "#FFA07A",   # salmon (8a)
    "9": "#A0522D",   # sienna (N/9b)
    "default": "#9DD3DF",
}
SCATTER_COLOR = "#B5179E"   # magenta (panel B)
PRIMER_F_COLOR = "#1F77B4"  # blue
PRIMER_R_COLOR = "#D62728"  # red
ALIGN_CONSERVED = "#2CA02C"  # green for conserved nucleotide
ALIGN_MISMATCH = "#D62728"   # red for mismatch


def _draw_genome_map(ax, ref_len: int, annotation: dict):
    """Panel A: stacked gene rectangles with labels."""
    ax.set_xlim(0, ref_len)
    ax.set_ylim(-2, 6)
    ax.set_yticks([])
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_position(("data", 0))
    # Two tracks: top for polyprotein cleavage products, bottom for ORFs
    genes = annotation.get("genes", [])
    for gene in genes:
        start = int(gene["start"])
        end = int(gene["end"])
        name = str(gene.get("name", ""))
        track = int(gene.get("track", 0))   # 0 = top, 1 = bottom
        y = 4 if track == 0 else 0
        h = 1.6
        color = gene.get("color") or _color_for_gene(name)
        ax.add_patch(Rectangle((start, y), end - start + 1, h,
                               facecolor=color, edgecolor="black", linewidth=0.4))
        ax.text((start + end) / 2, y + h / 2, name, ha="center", va="center",
                fontsize=8, fontweight="bold", color="black")

    # x-axis ticks (kb)
    ticks = np.linspace(0, ref_len, 7)
    labels = [f"{int(t / 1000):,}" for t in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8)
    ax.tick_params(axis="x", which="both", length=0)
    ax.set_xlabel("Position (nt)", fontsize=9)


def _color_for_gene(name: str) -> str:
    n = name.upper().replace("ORF", "").replace("GENE", "")
    n = n.strip()
    for key in COLORS:
        if key != "default" and key in n:
            return COLORS[key]
    return COLORS["default"]


def _draw_identity_scatter(ax, identity: np.ndarray, ref_len: int, threshold: float):
    """Panel B: per-nucleotide identity scatter."""
    pos = np.arange(1, ref_len + 1)
    ax.scatter(pos, identity, s=2, c=SCATTER_COLOR, alpha=0.6, rasterized=True, linewidths=0)
    ax.set_xlim(0, ref_len)
    ax.set_ylim(0.90, 1.001)
    ax.set_ylabel("Percentage identity per nucleotide", fontsize=9)
    ax.set_yticks(np.arange(0.90, 1.01, 0.01))
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="x", labelbottom=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.axhline(threshold, color="gray", linestyle="--", linewidth=0.6, alpha=0.6)


def _draw_region_closeup(ax, ref_len: int, regions: List[Tuple[int, int, float]],
                        chosen_idx: int, identity: np.ndarray):
    """Panel C: zoom into a conserved region; show per-position identity
    and a connector line down from panel B."""
    if not regions:
        ax.text(0.5, 0.5, "no conserved regions", ha="center", va="center",
                transform=ax.transAxes, fontsize=9)
        ax.set_axis_off()
        return
    start, end, mean_id = regions[chosen_idx]
    pad = max(10, (end - start) // 4)
    lo = max(1, start - pad)
    hi = min(ref_len, end + pad)
    sub_pos = np.arange(lo, hi + 1)
    sub_id = identity[lo - 1:hi]
    ax.plot(sub_pos, sub_id, "-", color=SCATTER_COLOR, linewidth=1.0, alpha=0.7)
    ax.fill_between(sub_pos, 0.90, sub_id, color=SCATTER_COLOR, alpha=0.15)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0.90, 1.001)
    ax.set_title(f"Conserved region {chosen_idx + 1}  ({start}-{end} nt, "
                 f"mean id={mean_id:.3f})", fontsize=9)
    ax.set_xlabel("Position (nt)", fontsize=8)
    ax.set_ylabel("Identity", fontsize=8)
    ax.tick_params(labelsize=7)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _render_alignment_panel(ax, ref_seq_str: str, records: Dict[str, SeqIO.SeqRecord],
                            ref_acc: str, region: Tuple[int, int, float],
                            max_rows: int = 8):
    """Panel D: alignment detail for a conserved region.

    For each sequence, show the subsequence and highlight mismatches vs
    the reference. Dots = match, red letter = mismatch. Reference is the
    top row in bold green.
    """
    from Bio import pairwise2
    start, end, _ = region
    ref_sub = ref_seq_str[start - 1:end]
    L = len(ref_sub)

    # collect sequences (truncate to max_rows including ref)
    other_accs = [a for a in records if a != ref_acc]
    other_accs = other_accs[: max_rows - 1]

    rows = [(ref_acc, ref_sub, True)]
    for acc in other_accs:
        rec_seq = str(records[acc].seq).upper()
        # quick local alignment to map this region's coords
        # for simplicity: if the record is long enough, take the slice
        if len(rec_seq) >= end:
            rows.append((acc, rec_seq[start - 1:end], False))
        else:
            # align full sequence vs reference substring
            aln = pairwise2.align.globalms(rec_seq, ref_sub, 2, -1, -5, -0.5, one_alignment_only=True)
            if not aln:
                continue
            rows.append((acc, aln[0].seqB, False))

    # plot
    n_rows = len(rows)
    row_h = 0.7
    ax.set_xlim(0, L)
    ax.set_ylim(-0.5, n_rows * row_h)
    ax.set_axis_off()
    # column ruler
    for x in range(0, L + 1, 10):
        ax.axvline(x, color="lightgray", linewidth=0.3, zorder=0)

    for i, (acc, sub, is_ref) in enumerate(rows):
        y = (n_rows - i - 1) * row_h + row_h / 2
        # label
        ax.text(-1.5, y, acc[:14], ha="right", va="center", fontsize=6,
                color="black", fontweight="bold" if is_ref else "normal")
        for j, base in enumerate(sub[:L]):
            x = j + 0.5
            if is_ref:
                ax.text(x, y, base, ha="center", va="center", fontsize=6,
                        fontweight="bold", color=ALIGN_CONSERVED)
            else:
                ref_base = ref_sub[j] if j < len(ref_sub) else "-"
                if base == ref_base or base == "-":
                    ax.text(x, y, ".", ha="center", va="center", fontsize=5,
                            color="gray")
                else:
                    ax.text(x, y, base, ha="center", va="center", fontsize=6,
                            color=ALIGN_MISMATCH, fontweight="bold")


def build_figure(
    ref_len: int,
    identity: np.ndarray,
    win_mean: np.ndarray,
    win_pos: np.ndarray,
    regions: List[Tuple[int, int, float]],
    annotation: dict,
    threshold: float,
    window_size: int,
    msa_records: Optional[list],
    records: Dict[str, SeqIO.SeqRecord],
    out_path: Path,
    logger: Optional[logging.Logger] = None,
):
    """Build and save the 4-panel conservation figure."""
    log = logger or logging.getLogger("ncbi_conservation")

    fig = plt.figure(figsize=(11, 14))
    gs = fig.add_gridspec(
        nrows=4, ncols=1,
        height_ratios=[1.0, 3.0, 1.6, 1.6],
        hspace=0.35,
    )
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    axC = fig.add_subplot(gs[2])
    axD = fig.add_subplot(gs[3])

    _draw_genome_map(axA, ref_len, annotation)
    _draw_identity_scatter(axB, identity, ref_len, threshold)

    # Connector: highlight chosen region in panel B
    chosen_idx = 0
    if regions:
        # pick the region with the highest mean identity for the close-up
        chosen_idx = int(np.argmax([r[2] for r in regions]))
        start, end, _ = regions[chosen_idx]
        axB.axvspan(start, end, color="gray", alpha=0.25, zorder=0)
        # add a connector from B down to C
        transB = axB.transData
        transFig = fig.transFigure.inverted()
        x0, y0 = transB.transform((start, 0.91))
        x1, y1 = transB.transform((end, 0.91))
        fx0, fy0 = transFig.transform((x0, y0))
        fx1, fy1 = transFig.transform((x1, y1))
        line = Line2D([fx0, fx1], [fy0, fy0], transform=fig.transFigure,
                      color="gray", linewidth=1.0, linestyle="-")
        fig.add_artist(line)
    _draw_region_closeup(axC, ref_len, regions, chosen_idx, identity)

    # Panel D: alignment detail
    if regions:
        ref_seq_str = str(records[next(iter(records))].seq).upper()
        # locate reference record
        ref_acc_in_records = None
        for acc in records:
            if str(records[acc].seq).upper() == ref_seq_str or len(records[acc].seq) == ref_len:
                ref_acc_in_records = acc
                break
        if ref_acc_in_records is None:
            ref_acc_in_records = next(iter(records))
        _render_alignment_panel(axD, ref_seq_str, records, ref_acc_in_records,
                                regions[chosen_idx])
    else:
        axD.text(0.5, 0.5, "no conserved regions to display",
                 ha="center", va="center", transform=axD.transAxes, fontsize=10)
        axD.set_axis_off()

    # Letter labels in the top-left corner of each panel
    for ax, label in zip([axA, axB, axC, axD], ["A", "B", "C", "D"]):
        ax.text(-0.08, 1.05, label, transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="bottom", ha="right")

    fig.suptitle("Nucleotide conservation analysis", fontsize=13, y=0.995)
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info(f"figure saved to {out_path}")
