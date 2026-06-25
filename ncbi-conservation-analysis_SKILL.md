---
name: ncbi-conservation-analysis
description: This skill should be used when the user wants to perform batch nucleotide sequence conservation analysis: download sequences from NCBI by keyword, run multiple sequence alignment, scan the genome for conserved regions using a 50 bp sliding window, generate a multi-panel figure (genome map + per-nucleotide identity scatter plot + conserved region close-up + alignment detail), and export the conserved regions as FASTA. Triggers include "NCBI 核酸下载", "多序列比对", "保守序列", "滑动窗口", "conservation analysis", "MSA + identity plot", "找保守区", "qRT-PCR 引物设计前导分析".
agent_created: true
---

# ncbi-conservation-analysis

End-to-end workflow for batch nucleotide conservation analysis: from keyword-based retrieval on NCBI, through multiple sequence alignment, to sliding-window conservation scoring, conserved-region extraction, and a publication-style multi-panel figure (modeled on SARS-CoV-2 reference papers).

## When to use this skill

Use this skill whenever the user wants to:

- Download a batch of nucleotide sequences from NCBI by keyword (e.g. "SARS-CoV-2 RdRp", "MERS-CoV spike").
- Run a multiple sequence alignment on the downloaded sequences.
- Identify conserved regions genome-wide with a 50 bp sliding window.
- Generate a multi-panel figure with: (A) genome map with annotated ORFs / genes, (B) per-nucleotide identity scatter plot, (C) conserved-region close-up, (D) alignment detail.
- Export the conserved regions as a FASTA file (typically for downstream qRT-PCR primer design).

## Inputs the user must provide

1. **Keyword** for NCBI search (e.g. "SARS-CoV-2 ORF1ab").
2. **Output folder** to save downloads, alignments, and figures.
3. **Reference sequence accession** (the sequence used to draw panels A and B; usually the longest or the type strain). Optional but strongly recommended.
4. **Optional gene annotations** (JSON or simple list) for panel A. If absent, the script will use a built-in SARS-CoV-2 / coronavirus annotation, or fall back to ORF1ab/ORF1a/S/3a/E/M/6/7a/7b/8a/8b/N/9b default for SARS-CoV-2.

## Workflow overview

1. **Download** — Query NCBI Nucleotide via Biopython's `Entrez.esearch` + `Entrez.efetch`; save per-accession FASTA files into `<output>/sequences/`.
2. **Select reference** — Pick the reference sequence (user-provided accession, or longest sequence).
3. **Pairwise alignment vs reference** — Use Biopython `pairwise2` or MAFFT (`--auto`) to align every sequence against the reference. Compute per-position identity; gaps are excluded from identity denominator (matches / (matches + mismatches)). Result: one identity value per reference position, averaged across all pairwise alignments.
4. **Sliding-window conservation** — Slide a 50 bp window (step 1 nt) over the per-position identity vector; report the mean identity in each window. Windows with identity ≥ threshold (default 0.99) are marked as conserved.
5. **Multi-sequence alignment** — Run a multiple alignment (MUSCLE or MAFFT) on the full set, restricted to the reference's coordinate system.
6. **Extract conserved regions** — Merge adjacent conserved windows; extract the corresponding subsequence (from the reference) for each region; write `conserved_regions.fasta`.
7. **Figure** — Build a 4-panel matplotlib figure (see `scripts/plot_conservation.py`).
8. **Cleanup** — Emit a summary CSV (`conservation_summary.csv`) and console report.

## Script entry point

`scripts/run_pipeline.py` orchestrates the full workflow. It accepts a JSON config file or CLI flags. Below is the canonical config:

```json
{
  "keyword": "SARS-CoV-2 RdRp",
  "output_dir": "./results",
  "reference_accession": "NC_045512.2",
  "max_sequences": 20,
  "window_size": 50,
  "window_step": 1,
  "identity_threshold": 0.99,
  "alignment_tool": "mafft",
  "gene_annotations": "sars-cov-2-default"
}
```

For convenience, the figure-generation logic is split into `scripts/plot_conservation.py` so it can be re-run with different thresholds without re-downloading sequences.

## Reusable contents

- `scripts/run_pipeline.py` — full orchestrator (download → align → score → plot → export).
- `scripts/download_ncbi.py` — NCBI retrieval with rate limiting and retries.
- `scripts/sliding_window_identity.py` — per-nucleotide identity + 50 bp sliding window.
- `scripts/plot_conservation.py` — multi-panel figure (panels A/B/C/D).
- `scripts/extract_conserved.py` — merge adjacent windows, export FASTA.
- `references/figure_layout.md` — specification of the 4-panel layout (font sizes, colors, panel ratios).
- `references/sars-cov-2-annotation.json` — default ORF annotations for SARS-CoV-2.
- `references/dependencies.md` — required Python packages and CLI tools.

## Operational notes

- **NCBI rate limit**: enforce ≥0.4 s between requests; set `Entrez.email` from a CLI flag.
- **Tool preference**: if `mafft` is on PATH, use it; otherwise fall back to Biopython's `MuscleCommandline`; if neither is available, fall back to a built-in progressive pairwise aligner (slower but works offline).
- **Long genomes**: per-nucleotide identity array for a 30 kb genome × 20 sequences is ~600 KB float64; trivially fits in memory.
- **Gap handling**: in pairwise identity, treat gaps in the *non-reference* sequence as missing data (do not penalize the reference). The script logs this choice clearly.
- **Threshold sensitivity**: 0.99 is a reasonable default for SARS-CoV-2 within-host diversity; for cross-species analyses (e.g. sarbecovirus spike), drop to 0.85.
- **Figure styling**: matplotlib defaults are overridden to mimic the reference figure — viridis-like gradient for genome map, magenta/purple dots for panel B, bold colored text for the alignment in panel D.

## Quick start

```bash
# Install dependencies (one time)
pip install biopython matplotlib numpy pandas

# Run the pipeline
python scripts/run_pipeline.py \
    --keyword "SARS-CoV-2 complete genome" \
    --output-dir ./results \
    --reference-accession NC_045512.2 \
    --max-sequences 20 \
    --window-size 50 \
    --identity-threshold 0.99
```

Outputs land in `./results/`:
- `sequences/*.fasta` — raw NCBI downloads
- `identity_per_nt.csv` — per-position identity vector
- `conservation_summary.csv` — sliding-window summary
- `conserved_regions.fasta` — extracted conserved subsequences
- `figure.png` — 4-panel publication figure
- `pipeline.log` — full run log

## Limitations

- Pairwise-vs-reference identity (used for panel B) is not identical to a column-wise MSA identity. For very divergent sequences the two can differ; the script logs both when an MSA is available.
- The default gene annotation is hard-coded for SARS-CoV-2; for other viruses, supply a custom `gene_annotations` JSON.
- This skill does not design primers — it only exports the conserved FASTA so a downstream tool (Primer3, Primer-BLAST) can be invoked separately.
