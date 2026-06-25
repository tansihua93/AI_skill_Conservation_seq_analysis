#!/usr/bin/env python
"""NCBI Nucleotide retrieval with rate limiting and per-accession FASTA caching.

Pulls FASTA records via Biopython's Entrez wrapper. Writes one FASTA per
accession to <out_dir>/<accession>.fasta. Honors NCBI rate-limit policy:
3 req/s without API key, 10 req/s with one.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Optional, Tuple, List

from Bio import Entrez, SeqIO


def _throttle(has_api_key: bool):
    time.sleep(0.34 if not has_api_key else 0.1)


def _search(keyword: str, max_n: int, email: str, api_key: Optional[str]) -> List[str]:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.esearch(db="nucleotide", term=keyword, retmax=max_n, usehistory="n")
    record = Entrez.read(handle)
    handle.Close()
    return record.get("IdList", [])


def _fetch_fasta(uid: str, email: str, api_key: Optional[str]) -> Optional[SeqIO.SeqRecord]:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    _throttle(bool(api_key))
    try:
        handle = Entrez.efetch(db="nucleotide", id=uid, rettype="fasta", retmode="text")
        rec = next(SeqIO.parse(handle, "fasta"), None)
        handle.Close()
        return rec
    except Exception as e:
        logging.getLogger("ncbi_conservation").warning(f"efetch failed for {uid}: {e}")
        return None


def _accession_from_record(rec: SeqIO.SeqRecord) -> str:
    # Prefer the versioned accession (.1, .2) when available.
    for feat in rec.features:
        pass
    acc = rec.id or rec.name
    # Biopython splits "NC_045512.2" -> id "NC_045512.2"; if no version, append .1
    if "." not in acc and rec.id:
        acc = f"{rec.id}.1"
    return acc.replace(".", "_").replace("|", "_")


def download_by_keyword(
    keyword: str,
    out_dir: Path,
    max_n: int,
    email: str,
    api_key: Optional[str] = None,
    reference_accession: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[str], Optional[str]]:
    """Download up to ``max_n`` sequences matching ``keyword``.

    Returns (accessions_list, picked_reference_accession_or_None).
    """
    log = logger or logging.getLogger("ncbi_conservation")
    out_dir.mkdir(parents=True, exist_ok=True)
    uids = _search(keyword, max_n, email, api_key)
    if not uids:
        log.warning(f"No UIDs returned for keyword='{keyword}'")
        return [], None

    accessions: List[str] = []
    picked_ref: Optional[str] = None
    for uid in uids:
        rec = _fetch_fasta(uid, email, api_key)
        if rec is None or len(rec.seq) < 200:
            continue
        acc = _accession_from_record(rec)
        out_path = out_dir / f"{acc}.fasta"
        SeqIO.write(rec, out_path, "fasta")
        accessions.append(acc)
        log.info(f"  fetched {acc}  ({len(rec.seq)} bp)  desc='{rec.description[:80]}'")
        if reference_accession and reference_accession.replace(".", "_") == acc:
            picked_ref = acc

    if reference_accession and picked_ref is None:
        log.warning(f"Reference accession {reference_accession} not found in result set; will pick longest.")

    return accessions, picked_ref
