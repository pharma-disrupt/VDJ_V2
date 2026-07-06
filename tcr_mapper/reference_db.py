"""
reference_db.py
===============

Loads the bundled IMGT germline FASTA files from `data/germline/` into
in-memory Python dictionaries on first access (cold start).

Public API:
    get_reference_db() -> ReferenceDB
    ReferenceDB:
        .v_genes: dict[str, dict[str, str]]   # locus -> allele -> seq
        .j_genes: dict[str, dict[str, str]]
        .d_genes: dict[str, dict[str, str]]   # currently empty (D calling is heuristic)
        .c_genes: dict[str, dict[str, str]]
        .constant_exemplars: dict[str, tuple[str, str]]  # locus -> (allele, seq)
        .release_info: dict[str, Any]
"""

from __future__ import annotations

import os
import threading
from pathlib import Path
from typing import Dict, Tuple, Any, List

from Bio import SeqIO


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

# Allow override via env var for non-Vercel deployments.
_DEFAULT_DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "germline"
GERMLINE_DATA_DIR = Path(os.environ.get("TCR_GERMLINE_DIR", str(_DEFAULT_DATA_DIR)))


# Locus -> FASTA file name mapping
_V_FILES = {
    "alpha": "TRAV.fasta",
    "beta":  "TRBV.fasta",
    "delta": "TRDV.fasta",
    "gamma": "TRGV.fasta",
}
_J_FILES = {
    "alpha": "TRAJ.fasta",
    "beta":  "TRBJ.fasta",
    "delta": "TRDJ.fasta",
    "gamma": "TRGJ.fasta",
}
_C_FILES = {
    "alpha":  "TRAC.fasta",
    "beta":   ["TRBC1.fasta", "TRBC2.fasta"],  # two C genes for beta
    "delta":  "TRDC.fasta",
    "gamma":  "TRGC.fasta",
}


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def _parse_header(header: str) -> Tuple[str, str, Dict[str, str]]:
    """
    Parse an IMGT-style FASTA header into (allele, gene, attrs).

    Expected format:  >ALLELE|gene=GENE|category=V|species=Homo+sapiens|...
    Falls back gracefully if the header is just the allele.
    """
    parts = header.split("|")
    allele = parts[0].lstrip(">").strip()
    attrs: Dict[str, str] = {}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            attrs[k.strip()] = v.strip()
    gene = attrs.get("gene", allele.split("*")[0] if "*" in allele else allele)
    return allele, gene, attrs


def _load_fasta(path: Path) -> Dict[str, Tuple[str, str]]:
    """
    Load a FASTA file into {allele: (gene, seq)}.

    If the same allele appears multiple times (e.g. duplicate IMGT records),
    the LAST occurrence wins.
    """
    if not path.exists():
        return {}
    out: Dict[str, Tuple[str, str]] = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        allele, gene, _ = _parse_header(rec.description)
        seq = str(rec.seq).upper().strip()
        if not seq:
            continue
        out[allele] = (gene, seq)
    return out


# ---------------------------------------------------------------------------
# ReferenceDB singleton
# ---------------------------------------------------------------------------

class ReferenceDB:
    """Holds all germline data in memory. Constructed once, reused."""

    def __init__(self, data_dir: Path) -> None:
        self.data_dir: Path = data_dir
        self.v_genes: Dict[str, Dict[str, str]] = {}   # locus -> allele -> seq
        self.j_genes: Dict[str, Dict[str, str]] = {}
        self.d_genes: Dict[str, Dict[str, str]] = {}   # placeholder; D calling is heuristic
        self.c_genes: Dict[str, Dict[str, str]] = {}
        self.v_gene_of_allele: Dict[str, str] = {}     # allele -> gene
        self.j_gene_of_allele: Dict[str, str] = {}
        self.c_gene_of_allele: Dict[str, str] = {}
        self.constant_exemplars: Dict[str, Tuple[str, str]] = {}
        self.release_info: Dict[str, Any] = {}
        self._load()

    # ------------------------------------------------------------------ load

    def _load(self) -> None:
        for locus, fname in _V_FILES.items():
            d = _load_fasta(self.data_dir / fname)
            self.v_genes[locus] = {allele: seq for allele, (gene, seq) in d.items()}
            for allele, (gene, _seq) in d.items():
                self.v_gene_of_allele[allele] = gene

        for locus, fname in _J_FILES.items():
            d = _load_fasta(self.data_dir / fname)
            self.j_genes[locus] = {allele: seq for allele, (gene, seq) in d.items()}
            for allele, (gene, _seq) in d.items():
                self.j_gene_of_allele[allele] = gene

        for locus, fnames in _C_FILES.items():
            if isinstance(fnames, str):
                fnames = [fnames]
            merged: Dict[str, str] = {}
            for fname in fnames:
                d = _load_fasta(self.data_dir / fname)
                for allele, (gene, seq) in d.items():
                    merged[allele] = seq
                    self.c_gene_of_allele[allele] = gene
            self.c_genes[locus] = merged

        # Constant exemplars for chain-type disambiguation
        ex_path = self.data_dir / "CONSTANT_exemplars.fasta"
        if ex_path.exists():
            for rec in SeqIO.parse(str(ex_path), "fasta"):
                allele, gene, attrs = _parse_header(rec.description)
                locus = attrs.get("locus", "").lower()
                if locus:
                    self.constant_exemplars[locus] = (allele, str(rec.seq).upper())

        # If no exemplar file, derive from c_genes (pick first allele per locus)
        for locus in ("alpha", "beta", "delta", "gamma"):
            if locus not in self.constant_exemplars and self.c_genes.get(locus):
                first_allele = next(iter(self.c_genes[locus]))
                self.constant_exemplars[locus] = (first_allele, self.c_genes[locus][first_allele])

        self.release_info = self._compute_release_info()

    # ------------------------------------------------------- release metadata

    def _compute_release_info(self) -> Dict[str, Any]:
        n_v = sum(len(d) for d in self.v_genes.values())
        n_j = sum(len(d) for d in self.j_genes.values())
        n_c = sum(len(d) for d in self.c_genes.values())
        return {
            "source": "IMGT/GENE-DB (embedded representative subset)",
            "species": "Homo sapiens",
            "data_dir": str(self.data_dir),
            "counts": {
                "V_alleles": n_v,
                "J_alleles": n_j,
                "C_alleles": n_c,
                "per_locus": {
                    "V": {k: len(v) for k, v in self.v_genes.items()},
                    "J": {k: len(v) for k, v in self.j_genes.items()},
                    "C": {k: len(v) for k, v in self.c_genes.items()},
                },
            },
            "notes": [
                "Embedded canonical IMGT functional allele translations.",
                "For full coverage, replace data/germline/*.fasta with the",
                "complete IMGT set downloaded from http://www.imgt.org/ligmdb/.",
            ],
        }

    # --------------------------------------------------------------- queries

    def get_v_panel(self, locus: str) -> Dict[str, str]:
        return self.v_genes.get(locus, {})

    def get_j_panel(self, locus: str) -> Dict[str, str]:
        return self.j_genes.get(locus, {})

    def get_c_panel(self, locus: str) -> Dict[str, str]:
        return self.c_genes.get(locus, {})

    def get_d_panel(self, locus: str) -> Dict[str, str]:
        # D-gene calling from protein is unreliable; we do not store D refs.
        # gene_mapper.py performs a heuristic search using the unrearranged
        # D-region mini-sequences defined locally.
        return {}


# ---------------------------------------------------------------------------
# Singleton accessor (thread-safe, lazy)
# ---------------------------------------------------------------------------

_DB: ReferenceDB | None = None
_DB_LOCK = threading.Lock()


def get_reference_db() -> ReferenceDB:
    """Return the singleton ReferenceDB, loading it on first call."""
    global _DB
    if _DB is None:
        with _DB_LOCK:
            if _DB is None:
                _DB = ReferenceDB(GERMLINE_DATA_DIR)
    return _DB


def reload_reference_db(data_dir: Path | None = None) -> ReferenceDB:
    """Force a reload (used by tests)."""
    global _DB
    with _DB_LOCK:
        _DB = ReferenceDB(data_dir or GERMLINE_DATA_DIR)
    return _DB
