# Figure Layout Specification

This document defines the 4-panel layout used by `plot_conservation.py`.
The layout is modeled on the SARS-CoV-2 conservation figure in the reference paper (Figure 4 in the user-supplied screenshot).

## Panel grid

```
+--------------------------------------------------+
|  A: genome map (height ratio 1.0)                |
+--------------------------------------------------+
|  B: per-nucleotide identity scatter (3.0)        |
+--------------------------------------------------+
|  C: conserved region close-up (1.6)               |
+--------------------------------------------------+
|  D: alignment detail (1.6)                       |
+--------------------------------------------------+
```

- Figure size: 11 in × 14 in.
- DPI: 200 (so the saved PNG is ~2200 × 2800 px).
- Vertical spacing (`hspace`) between panels: 0.35.
- Letter labels (A, B, C, D) at the top-left of each panel, fontsize 14, bold.

## Panel A — Genome map

- Two horizontal tracks:
  - Track 0 (top, y=4, height 1.6): polyprotein cleavage products + functional domains (PLpro, 3CLpro, RdRp, Hel, ExoN).
  - Track 1 (bottom, y=0, height 1.6): mature ORFs (ORF1a, ORF1b, S, 3a, E, M, 6, 7a, 7b, 8a, 8b, N, 9b).
- Each rectangle is filled with a color from the COLORS palette and labeled with the gene name in 8 pt bold.
- x-axis: 7 evenly spaced ticks labeled in kb (0, 5, 10, …).
- y-axis: hidden (no ticks, no labels).
- Top/right spines hidden.

## Panel B — Per-nucleotide identity scatter

- Scatter plot, s=2, alpha=0.6, color magenta (#B5179E).
- y-axis: 0.90 - 1.00 with 0.01 ticks; label "Percentage identity per nucleotide".
- x-axis: 0 - reference length; tick labels hidden (panel C shows them).
- Threshold dashed line at the configured identity threshold.
- Highlight bar over the chosen conserved region (gray, alpha 0.25).
- Connector line drawn in figure coordinates from the bar down to panel C.

## Panel C — Conserved region close-up

- Title: "Conserved region N (start-end nt, mean id=X.XXX)".
- Line plot of per-position identity in the padded region (pad = max(10, len/4)).
- Filled area below the line, alpha 0.15.
- y-axis 0.90 - 1.00; x-axis in nt.
- The "chosen" region is the one with the highest mean identity (auto-selected).

## Panel D — Alignment detail

- One row per sequence, max 8 rows (including reference on top).
- Reference row: bold green letters.
- Other rows: dot for match, red bold letter for mismatch, gap as '-'.
- Column ruler every 10 columns.
- Sequence labels on the left, truncated to 14 chars.

## Colors

| Element | Color |
|---|---|
| ORF1a / ORF1b | #56B4E9 (sky blue) |
| S | #F08080 (light coral) |
| 3a | #FFD580 (peach) |
| E | #90EE90 (light green) |
| M | #C8A2C8 (light purple) |
| 6 | #A9A9A9 (gray) |
| 7a / 7b / 8a / 8b | #FFA07A (salmon) |
| N / 9b | #A0522D (sienna) |
| Identity scatter | #B5179E (magenta) |
| Reference in alignment | #2CA02C (bold green) |
| Mismatch in alignment | #D62728 (bold red) |

## How to customize

To add a new virus annotation, create a JSON file with the same structure as `sars-cov-2-annotation.json` and pass it via `--annotation /path/to/your.json`. Each gene entry must have `name`, `start`, `end`, `track` (0 or 1) and optional `color`.
