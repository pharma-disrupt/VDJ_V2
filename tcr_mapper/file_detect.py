"""
file_detect.py
==============

Best-guess detection of uploaded file format and likely molecule content.

Strategy:
  1. Inspect file extension.
  2. Sniff the first ~2 KB of content to confirm format (PDB ATOM/HEADER,
     mmCIF _atom_site, FASTA '>'-prefixed headers).
  3. Suggest a molecule type based on format + content hints.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Tuple, List

from tcr_mapper.models import FileTypeGuess


# Mapping file extension -> likely format
_EXT_MAP = {
    ".pdb":  "pdb",
    ".ent":  "pdb",
    ".cif":  "mmcif",
    ".mmcif": "mmcif",
    ".fasta": "fasta",
    ".fa":   "fasta",
    ".faa":  "fasta",
    ".fna":  "fasta",
    ".txt":  "fasta",   # ambiguous, sniff to confirm
}


def sniff_content(content_bytes: bytes) -> Tuple[str, List[str]]:
    """
    Inspect the first ~2KB of file content to determine format.

    Returns (format, notes).
    """
    notes: List[str] = []
    if not content_bytes:
        return "unknown", ["empty content"]

    try:
        head = content_bytes[:2048].decode("utf-8", errors="ignore")
    except Exception:  # noqa: BLE001
        return "unknown", ["unreadable content"]

    stripped = head.lstrip()
    if not stripped:
        return "unknown", ["only whitespace"]

    # FASTA: starts with '>' after optional whitespace
    if stripped.startswith(">"):
        return "fasta", ["starts with '>'"]

    # PDB: look for HEADER/ATOM/HETATM/COMPND records
    pdb_markers = ("HEADER", "ATOM  ", "HETATM", "COMPND", "CRYST1", "SEQRES")
    pdb_hits = [m for m in pdb_markers if m in head[:600]]
    if pdb_hits:
        return "pdb", [f"contains PDB record: {','.join(pdb_hits)}"]

    # mmCIF: look for _atom_site. or data_ block
    if "_atom_site." in head or head.lstrip().startswith("data_"):
        return "mmcif", ["contains mmCIF _atom_site or data_ block"]

    # Heuristic: if it's mostly ACGT and lines look like sequences, guess fasta
    alnum_chars = sum(1 for c in head[:512] if c.isalnum())
    if alnum_chars > 0.5 * len(head[:512]) and any(c in head for c in "ACGTU"):
        # Could be a sequence without '>' header — accept as fasta
        return "fasta", ["no '>' but mostly alphanumeric (treating as raw sequence)"]

    return "unknown", ["no recognized markers"]


def suggest_molecule(file_format: str, content_bytes: bytes) -> Tuple[str, List[str]]:
    """
    Suggest molecule content type based on format + content hints.
    Returns (molecule, notes).
    """
    notes: List[str] = []
    if file_format == "fasta":
        return "fasta", ["FASTA always maps to TCR sequence"]

    if file_format in ("pdb", "mmcif"):
        try:
            head = content_bytes[:8192].decode("utf-8", errors="ignore").upper()
        except Exception:  # noqa: BLE001
            head = ""

        # Look for MHC chain identifiers in TITLE/COMPND lines.
        has_mhc = ("MHC" in head) or ("HLA" in head) or ("H-2" in head)
        has_peptide_hint = "PEPTIDE" in head or "B2M" in head

        # Chain count heuristic: count ATOM lines that introduce a new chain id.
        # (Simplified — the pdb_parser does the real chain extraction.)
        if has_mhc or has_peptide_hint:
            return "pmhc_tcr_complex", [
                "TITLE/COMPND mentions MHC/HLA/B2M/peptide"
            ]

        # Default for structure without obvious MHC markers: TCR-only guess
        return "tcr_only", [
            "no MHC/HLA markers in TITLE/COMPND (defaulting to TCR-only)"
        ]

    return "unknown", ["unrecognized format"]


def detect_file(filename: str, content_bytes: bytes) -> FileTypeGuess:
    """Top-level detection combining extension + content sniff."""
    ext = Path(filename).suffix.lower()
    ext_format = _EXT_MAP.get(ext)
    sniff_format, sniff_notes = sniff_content(content_bytes)

    # Prefer content sniff if extension is missing/ambiguous; otherwise trust ext
    # unless sniff strongly disagrees.
    if ext_format and sniff_format != "unknown" and ext_format != sniff_format:
        # Disagreement — trust the sniff (more reliable) but log a note.
        file_format = sniff_format
        notes = sniff_notes + [f"extension suggested {ext_format} but content looks like {sniff_format}"]
    elif ext_format:
        file_format = ext_format
        notes = sniff_notes if sniff_format != "unknown" else ["extension-based; content unrecognized"]
    else:
        file_format = sniff_format
        notes = sniff_notes + [f"unrecognized extension '{ext}'"]

    molecule, mol_notes = suggest_molecule(file_format, content_bytes)
    notes.extend(mol_notes)

    return FileTypeGuess(
        file_format=file_format,
        suggested_molecule=molecule,
        sniffed_notes=notes,
        filename=filename,
        size_bytes=len(content_bytes),
    )
