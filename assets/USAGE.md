# Example: Run on SARS-CoV-2

This directory contains a ready-to-use config for analyzing SARS-CoV-2 conservation. Copy it to your project, edit the email, and run:

```bash
python scripts/run_pipeline.py --config assets/config_template.json
```

The pipeline will:
1. Search NCBI for "SARS-CoV-2 complete genome".
2. Download the first 20 hits (or fewer if matches are limited).
3. Pick NC_045512.2 as the reference (or fall back to the longest sequence).
4. Compute per-nucleotide pairwise identity across all 20 sequences.
5. Slide a 50 bp window and mark regions with mean identity >= 0.99.
6. Run a MAFFT MSA (if installed) on the full set.
7. Extract conserved regions to `conserved_regions.fasta`.
8. Build a 4-panel PNG figure in `./results/figure.png`.

Expected runtime on a modern laptop: 1 - 3 minutes (dominated by NCBI downloads and pairwise alignments).

## Customizing for a different virus

To analyze MERS-CoV, influenza HA, or any other gene:

```bash
python scripts/run_pipeline.py \
    --keyword "MERS-CoV spike glycoprotein complete cds" \
    --output-dir ./mers_results \
    --max-sequences 30 \
    --window-size 50 \
    --identity-threshold 0.85
```

You will also need to provide a custom gene annotation JSON. See `references/figure_layout.md` for the schema.
