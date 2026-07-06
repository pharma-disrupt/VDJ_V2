"""
fasta_parser.py
================

Parses a FASTA file (or raw sequence text) into TCRChainInput objects.

Handles:
  - Standard FASTA (>header\nSEQ).
  - Bare sequence text (no '>' header).
  - Nucleotide sequences (auto-translated in the best frame).
  - Multi-record FASTA (one chain per record).
"""

from __future__ import annotations

from typing import List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import io

from tcr_mapper.models import TCRChainInput


# Codon table fallback for non-ATGC chars
_AMBIGUOUS_AA = "X"


def _looks_like_nucleotide(seq: str) -> bool:
    """Heuristic: a sequence is nucleotide if it's >90% ACGTN."""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    if not seq:
        return False
    # If it contains characters outside ACGTN-, it's protein
    non_nt = sum(1 for c in seq if c not in "ACGTN-")
    return non_nt / len(seq) < 0.10


def _translate(seq: str) -> str:
    """
    Translate a nucleotide sequence to protein, trying all three forward
    frames and returning the longest ORF that starts with M or contains
    a long stretch without stops.
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")
    # Strip non-ACGTN
    seq = "".join(c for c in seq if c in "ACGTN")
    if len(seq) < 3:
        return ""

    best_protein = ""
    for frame in range(3):
        sub = seq[frame:]
        # Pad to multiple of 3
        if len(sub) % 3:
            sub = sub[: len(sub) - (len(sub) % 3)]
        try:
            protein = str(Seq(sub).translate(to_stop=False))
        except Exception:  # noqa: BLE001
            continue
        # Replace stop codons with nothing for scoring; pick frame with fewest stops
        # and longest continuous stretch.
        stops = protein.count("*")
        # Score: prefer frames that start with M, have few stops, and are long
        score = len(protein) - 5 * stops
        if protein.startswith("M"):
            score += 10
        if score > len(best_protein):
            best_protein = protein
    return best_protein


def parse_fasta(content: str, source: str = "fasta") -> Tuple[List[TCRChainInput], List[str]]:
    """
    Parse FASTA content into TCRChainInput objects.

    Args:
        content: raw text (FASTA or bare sequence).
        source: the source label to attach (default "fasta").

    Returns:
        (list of TCRChainInput, list of warnings/notes)
    """
    warnings: List[str] = []
    chains: List[TCRChainInput] = []

    content_stripped = content.strip()
    if not content_stripped:
        return chains, ["empty FASTA content"]

    # Wrap in StringIO for SeqIO
    has_header = content_stripped.startswith(">")
    if not has_header:
        # Bare sequence — synthesize a header
        content = ">chain_1\n" + content_stripped
        warnings.append("input had no FASTA header; assigned chain_id='chain_1'")

    handle = io.StringIO(content)
    records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        return chains, ["no FASTA records parsed"]

    for idx, rec in enumerate(records):
        seq = str(rec.seq).upper().strip()
        if not seq:
            warnings.append(f"record {idx} ({rec.id}) has empty sequence; skipped")
            continue

        # Decide nt vs aa
        if _looks_like_nucleotide(seq):
            translated = _translate(seq)
            if not translated:
                warnings.append(f"record {rec.id}: nucleotide but translation failed; skipped")
                continue
            warnings.append(f"record {rec.id}: detected nucleotide, translated ({len(seq)} nt -> {len(translated)} aa)")
            seq = translated

        # Clean: keep only standard AA alphabet
        seq_clean = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in seq)

        chain_id = rec.id or f"chain_{idx+1}"
        # Use first word of id to avoid spaces
        chain_id = chain_id.split()[0] if chain_id else f"chain_{idx+1}"

        chains.append(TCRChainInput(
            chain_id=chain_id,
            sequence=seq_clean,
            source=source,
            inferred_chain_type=None,
        ))

    return chains, warnings
