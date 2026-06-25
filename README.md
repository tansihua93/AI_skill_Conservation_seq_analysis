# Nucleic_acid_Conservation_seq_analysis
End-to-end workflow for batch nucleotide conservation analysis: from keyword-based retrieval on NCBI, through multiple sequence alignment, to sliding-window conservation scoring, conserved-region extraction, and a publication-style multi-panel figure.

1. **Download** — Query NCBI Nucleotide via Biopython's `Entrez.esearch` + `Entrez.efetch`; save per-accession FASTA files into `<output>/sequences/`.   
2. **Select reference** — Pick the reference sequence (user-provided accession, or longest sequence).
3. **Pairwise alignment vs reference** — Use Biopython `pairwise2` or MAFFT (`--auto`) to align every sequence against the reference. Compute per-position identity; gaps are excluded from identity denominator (matches / (matches + mismatches)). Result: one identity value per reference position, averaged across all pairwise alignments.
4. **Sliding-window conservation** — Slide a 50 bp window (step 1 nt) over the per-position identity vector; report the mean identity in each window. Windows with identity ≥ threshold (default 0.99) are marked as conserved.
5. **Multi-sequence alignment** — Run a multiple alignment (MUSCLE or MAFFT) on the full set, restricted to the reference's coordinate system.
6. **Extract conserved regions** — Merge adjacent conserved windows; extract the corresponding subsequence (from the reference) for each region; write `conserved_regions.fasta`.
7. **Figure** — Build a 4-panel matplotlib figure (see `scripts/plot_conservation.py`).
8. **Cleanup** — Emit a summary CSV (`conservation_summary.csv`) and console report.
