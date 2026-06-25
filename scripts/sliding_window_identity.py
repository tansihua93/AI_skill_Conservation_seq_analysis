#!/usr/bin/env python
"""Per-nucleotide pairwise identity and 50 bp sliding-window conservation.

Design:
- For each non-reference sequence, perform a global pairwise alignment vs
  the reference. (Biopython `pairwise2` is used; for long sequences
  consider MAFFT --auto in a streaming mode.)
- For each aligned column, count the number of sequences whose aligned
  base matches the reference base, divided by the number of sequences
  that have a non-gap base at that column. Gaps are excluded from the
  denominator (they do not count as mismatch).
- Sliding-window mean over the resulting identity vector.
- Merge adjacent windows above the threshold into contiguous regions.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Sequence, Tuple

import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq


def _global_align(ref: str, qry: str) -> Tuple[str, str]:
    """Global pairwise alignment. Falls back to identity match if alignment fails."""
    try:
        aln = pairwise2.align.globalms(ref, qry, 2, -1, -5, -0.5, one_alignment_only=True)[0]
        return aln.seqA, aln.seqB
    except Exception:
        # length-mismatch fallback: pad shorter with '-'
        n = max(len(ref), len(qry))
        return ref.ljust(n, "-"), qry.ljust(n, "-")


def compute_pairwise_identity(
    ref_seq: Seq,
    records: Dict[str, SeqIO.SeqRecord],
    ref_acc: str,
    logger: logging.Logger = None,
) -> List[float]:
    """Return per-reference-position identity in [0, 1]. Length == len(ref_seq).

    Gaps in the non-reference sequence are treated as missing data and
    excluded from the denominator. Identity at position i is
        matches_i / (matches_i + mismatches_i)
    where a sequence contributes to the denominator only if its aligned
    base at that column is not '-'.
    """
    log = logger or logging.getLogger("ncbi_conservation")
    ref_str = str(ref_seq).upper()
    n = len(ref_str)
    matches = np.zeros(n, dtype=np.int64)
    covered = np.zeros(n, dtype=np.int64)

    for acc, rec in records.items():
        if acc == ref_acc:
            continue
        qry_str = str(rec.seq).upper()
        log.info(f"  aligning {acc} (len={len(qry_str)}) vs reference ...")
        a, b = _global_align(ref_str, qry_str)
        if len(a) != len(b):
            log.warning(f"  alignment length mismatch for {acc}; skipping")
            continue
        # Walk the alignment, track ref position
        ref_pos = -1
        for ca, cb in zip(a, b):
            if ca != "-":
                ref_pos += 1
            if ref_pos < 0 or ref_pos >= n:
                continue
            if cb == "-":
                continue  # gap in query -> missing data
            covered[ref_pos] += 1
            if ca == cb:
                matches[ref_pos] += 1

    with np.errstate(divide="ignore", invalid="ignore"):
        identity = np.where(covered > 0, matches / np.maximum(covered, 1), np.nan)
    # fill NaN with 1.0 (uncovered positions = no information, default to conserved)
    identity = np.nan_to_num(identity, nan=1.0)
    return identity.tolist()


def sliding_window_mean(values: np.ndarray, window: int = 50, step: int = 1) -> np.ndarray:
    """Compute a moving average over ``values`` using a uniform window.

    Returns an array of length (N - window) // step + 1.
    """
    if window <= 0 or window > len(values):
        raise ValueError(f"window={window} invalid for N={len(values)}")
    if step <= 0:
        raise ValueError("step must be positive")
    n = len(values)
    n_out = (n - window) // step + 1
    if n_out <= 0:
        return np.array([], dtype=float)
    # cumulative sum for fast windowed mean
    csum = np.concatenate([[0.0], np.cumsum(values, dtype=float)])
    out = np.empty(n_out, dtype=float)
    for i in range(n_out):
        s = i * step
        e = s + window
        out[i] = (csum[e] - csum[s]) / window
    return out


def merge_conserved_windows(
    win_mean: np.ndarray,
    win_pos: np.ndarray,
    window_size: int,
    threshold: float,
) -> List[Tuple[int, int, float]]:
    """Merge adjacent conserved windows into contiguous regions.

    Returns a list of (start, end, mean_identity) tuples in reference
    coordinates (1-based inclusive). Adjacent windows are merged if their
    mean identity is >= threshold.
    """
    half = window_size // 2
    regions: List[Tuple[int, int, float]] = []
    cur_start = None
    cur_end = None
    cur_vals: List[float] = []
    for center, mean in zip(win_pos, win_mean):
        if mean >= threshold:
            start = int(center - half)
            end = int(center + half - (window_size % 2 == 0))
            if cur_start is None:
                cur_start, cur_end = start, end
                cur_vals = [mean]
            else:
                # merge if overlapping or adjacent
                if start <= cur_end + 1:
                    cur_end = max(cur_end, end)
                    cur_vals.append(mean)
                else:
                    regions.append((cur_start, cur_end, float(np.mean(cur_vals))))
                    cur_start, cur_end = start, end
                    cur_vals = [mean]
    if cur_start is not None:
        regions.append((cur_start, cur_end, float(np.mean(cur_vals))))
    return regions
