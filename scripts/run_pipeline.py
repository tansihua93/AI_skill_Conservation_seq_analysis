#!/usr/bin/env python
"""End-to-end pipeline for batch nucleotide conservation analysis.

Workflow:
1. Download sequences from NCBI by keyword.
2. Pick / validate reference sequence.
3. Pairwise-align every sequence against reference; compute per-position identity.
4. Sliding-window conservation (50 bp default).
5. Multiple sequence alignment of full set.
6. Extract conserved regions to FASTA.
7. Build 4-panel figure (genome map / identity scatter / close-up / alignment detail).
8. Emit CSV summary and logs.

Usage:
    python run_pipeline.py --keyword "SARS-CoV-2 complete genome" \\
        --output-dir ./results --reference-accession NC_045512.2 \\
        --max-sequences 20 --window-size 50 --identity-threshold 0.99
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

from download_ncbi import download_by_keyword
from sliding_window_identity import (
    compute_pairwise_identity,
    sliding_window_mean,
    merge_conserved_windows,
)
from extract_conserved import write_conserved_fasta
from plot_conservation import build_figure


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Batch nucleotide conservation pipeline (NCBI -> MSA -> 50bp window -> figure)."
    )
    p.add_argument("--keyword", required=True, help="NCBI search keyword.")
    p.add_argument("--output-dir", required=True, help="Output folder.")
    p.add_argument("--reference-accession", default=None,
                   help="Reference accession (FASTA header token). If omitted, longest sequence is used.")
    p.add_argument("--max-sequences", type=int, default=20)
    p.add_argument("--window-size", type=int, default=50)
    p.add_argument("--window-step", type=int, default=1)
    p.add_argument("--identity-threshold", type=float, default=0.99)
    p.add_argument("--alignment-tool", default="auto", choices=["auto", "mafft", "muscle", "builtin"])
    p.add_argument("--email", default="researcher@example.com",
                   help="NCBI Entrez email (required by NCBI policy).")
    p.add_argument("--api-key", default=None, help="NCBI API key for higher rate limit.")
    p.add_argument("--config", default=None, help="Optional JSON config overriding CLI flags.")
    p.add_argument("--annotation", default="sars-cov-2-default",
                   help="Gene annotation preset or path to JSON file.")
    return p.parse_args()


def load_config(args: argparse.Namespace) -> dict:
    cfg = vars(args).copy()
    if args.config and Path(args.config).is_file():
        with open(args.config, "r", encoding="utf-8") as f:
            override = json.load(f)
        cfg.update(override)
    return cfg


def setup_logging(output_dir: Path) -> logging.Logger:
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "pipeline.log"
    logger = logging.getLogger("ncbi_conservation")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    sh = logging.StreamHandler(sys.stdout)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh.setFormatter(fmt); sh.setFormatter(fmt)
    logger.addHandler(fh); logger.addHandler(sh)
    return logger


def main() -> int:
    args = parse_args()
    cfg = load_config(args)
    output_dir = Path(cfg["output_dir"]).resolve()
    logger = setup_logging(output_dir)
    logger.info("=== ncbi-conservation-analysis pipeline start ===")
    for k, v in cfg.items():
        logger.info(f"config: {k} = {v}")

    # ---- Step 1: download ----
    seq_dir = output_dir / "sequences"
    seq_dir.mkdir(parents=True, exist_ok=True)
    accessions, ref_record = download_by_keyword(
        keyword=cfg["keyword"],
        out_dir=seq_dir,
        max_n=cfg["max_sequences"],
        email=cfg["email"],
        api_key=cfg["api_key"],
        reference_accession=cfg.get("reference_accession"),
        logger=logger,
    )
    if not accessions:
        logger.error("No sequences downloaded. Aborting.")
        return 1
    logger.info(f"Downloaded {len(accessions)} sequences: {accessions}")

    # ---- Step 2: load all FASTA ----
    from Bio import SeqIO
    records = {}
    for acc in accessions:
        fasta_path = seq_dir / f"{acc}.fasta"
        if not fasta_path.is_file():
            # find a file containing this accession
            for f in seq_dir.glob("*.fasta"):
                if acc in f.stem:
                    fasta_path = f
                    break
        rec = next(SeqIO.parse(fasta_path, "fasta"))
        records[acc] = rec
    ref_acc = cfg.get("reference_accession") or ref_record or max(
        records, key=lambda a: len(records[a].seq)
    )
    ref_seq = records[ref_acc].seq
    ref_len = len(ref_seq)
    logger.info(f"Reference: {ref_acc}  length={ref_len} bp")

    # ---- Step 3: per-position pairwise identity ----
    identity = compute_pairwise_identity(ref_seq, records, ref_acc, logger=logger)
    import numpy as np
    import pandas as pd
    id_arr = np.array(identity, dtype=float)
    id_path = output_dir / "identity_per_nt.csv"
    pd.DataFrame({"position": np.arange(1, ref_len + 1), "identity": id_arr}).to_csv(
        id_path, index=False
    )
    logger.info(f"Wrote per-position identity -> {id_path}")

    # ---- Step 4: sliding window ----
    win_mean = sliding_window_mean(id_arr, window=cfg["window_size"], step=cfg["window_step"])
    win_pos = np.arange(
        cfg["window_size"] // 2 + 1,
        ref_len - cfg["window_size"] // 2 + 1,
        cfg["window_step"],
    )[: len(win_mean)]
    summary_path = output_dir / "conservation_summary.csv"
    pd.DataFrame({
        "window_center": win_pos,
        "window_start": win_pos - cfg["window_size"] // 2,
        "window_end": win_pos + cfg["window_size"] // 2,
        "mean_identity": win_mean,
    }).to_csv(summary_path, index=False)
    logger.info(f"Wrote sliding-window summary -> {summary_path}")

    # ---- Step 5: multiple sequence alignment (best-effort) ----
    msa_records = None
    try:
        if cfg["alignment_tool"] in ("auto", "mafft"):
            from Bio.Align.Applications import MafftCommandline
            from Bio import AlignIO
            combined = output_dir / "_combined.fasta"
            with open(combined, "w") as f:
                for rec in records.values():
                    f.write(f">{rec.id}\n{str(rec.seq)}\n")
            mafft_cline = MafftCommandline(input=str(combined), auto=True)
            stdout, _ = mafft_cline()
            aln_path = output_dir / "aligned.fasta"
            aln_path.write_text(stdout)
            msa_records = list(AlignIO.read(aln_path, "fasta"))
            logger.info(f"MAFFT MSA -> {aln_path}  ({len(msa_records)} seqs, {msa_records[0].get_alignment_length()} cols)")
    except Exception as e:
        logger.warning(f"MSA step failed ({e}); continuing with pairwise identity only.")

    # ---- Step 6: extract conserved regions ----
    regions = merge_conserved_windows(
        win_mean, win_pos, cfg["window_size"], cfg["identity_threshold"]
    )
    fasta_out = output_dir / "conserved_regions.fasta"
    write_conserved_fasta(ref_seq, ref_acc, regions, fasta_out, logger)
    logger.info(f"Extracted {len(regions)} conserved regions -> {fasta_out}")

    # ---- Step 7: figure ----
    # load annotation
    annot = None
    if cfg["annotation"] and Path(cfg["annotation"]).is_file():
        annot = json.loads(Path(cfg["annotation"]).read_text(encoding="utf-8"))
    else:
        # default SARS-CoV-2
        default_path = Path(__file__).parent.parent / "references" / "sars-cov-2-annotation.json"
        if default_path.is_file():
            annot = json.loads(default_path.read_text(encoding="utf-8"))
    if annot is None:
        annot = {"genes": []}

    fig_path = output_dir / "figure.png"
    build_figure(
        ref_len=ref_len,
        identity=id_arr,
        win_mean=win_mean,
        win_pos=win_pos,
        regions=regions,
        annotation=annot,
        threshold=cfg["identity_threshold"],
        window_size=cfg["window_size"],
        msa_records=msa_records,
        records=records,
        out_path=fig_path,
        logger=logger,
    )
    logger.info(f"Figure -> {fig_path}")

    # ---- Step 8: summary ----
    logger.info("=== summary ===")
    logger.info(f"  reference: {ref_acc}  ({ref_len} bp)")
    logger.info(f"  sequences aligned vs reference: {len(records) - 1}")
    logger.info(f"  conserved regions (>= {cfg['identity_threshold']:.2f}): {len(regions)}")
    total_conserved_bp = sum(end - start + 1 for start, end, _ in regions)
    logger.info(f"  total conserved bases: {total_conserved_bp} bp "
                f"({100.0 * total_conserved_bp / ref_len:.2f}% of reference)")
    logger.info("=== done ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
