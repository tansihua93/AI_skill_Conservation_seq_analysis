#!/usr/bin/env python
"""Extract conserved regions and write a FASTA file.

Regions are defined as merged sliding-window blocks with mean identity
above the threshold. We extract the reference subsequence (1-based
inclusive coordinates) for each region and emit a FASTA record per
region with a descriptive header.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

from Bio.Seq import Seq


def write_conserved_fasta(
    ref_seq: Seq,
    ref_acc: str,
    regions: List[Tuple[int, int, float]],
    out_path: Path,
    logger: logging.Logger = None,
) -> int:
    """Write conserved subsequences to a FASTA file. Returns number written."""
    log = logger or logging.getLogger("ncbi_conservation")
    ref_str = str(ref_seq)
    n_written = 0
    with open(out_path, "w", encoding="utf-8") as f:
        for i, (start, end, mean_id) in enumerate(regions, start=1):
            # 1-based inclusive to 0-based slice
            sub = ref_str[start - 1:end]
            length = end - start + 1
            header = (
                f">{ref_acc}_conserved_region_{i} "
                f"start={start} end={end} length={length} mean_identity={mean_id:.4f}"
            )
            # split sequence into 80-char lines
            f.write(header + "\n")
            for j in range(0, len(sub), 80):
                f.write(sub[j:j + 80] + "\n")
            n_written += 1
    log.info(f"  wrote {n_written} conserved regions to {out_path}")
    return n_written
