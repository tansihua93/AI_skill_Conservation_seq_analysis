# Pipeline usage notes

Detailed recipes for common scenarios:

1. **SARS-CoV-2 complete genome** — `assets/config_template.json` is ready to use.
2. **MERS-CoV / SARS-CoV / influenza HA** — supply a custom `--annotation` JSON; see `references/figure_layout.md` for schema.
3. **Large batch (50+ sequences)** — increase `--max-sequences`, optionally provide an `--api-key` for higher NCBI rate limit.
4. **Cross-species comparison** — lower `--identity-threshold` to 0.80-0.90.

For primer-design follow-up, feed `conserved_regions.fasta` into Primer3 or Primer-BLAST.
