"""
pipeline.py
===========

Top-level orchestrator. Takes raw file content + process options, dispatches
to the right parser, runs gene mapping, returns a PipelineResult.
"""

from __future__ import annotations

from typing import List, Optional, Dict, Any

from tcr_mapper.models import (
    TCRChainInput,
    GeneCallResult,
    PipelineResult,
    ProcessOptions,
)
from tcr_mapper.file_detect import detect_file
from tcr_mapper.fasta_parser import parse_fasta
from tcr_mapper.pdb_parser import parse_structure
from tcr_mapper.gene_mapper import map_chain
from tcr_mapper.reference_db import get_reference_db


def run_pipeline(
    filename: str,
    content_bytes: bytes,
    options: ProcessOptions,
) -> PipelineResult:
    """
    Run the full TCR germline mapping pipeline.

    Args:
        filename: original filename (used for logging only)
        content_bytes: raw file content as bytes
        options: confirmed user options (file_format, molecule, tcr_chain_ids)
    """
    result = PipelineResult(
        file_format=options.file_format,
        molecule=options.molecule,
    )

    # Decode content to text
    try:
        content_text = content_bytes.decode("utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        try:
            content_text = content_bytes.decode("latin-1", errors="replace")
        except Exception:  # noqa: BLE001
            result.warnings.append("could not decode file content as text")
            return result

    # ----- Parse into TCRChainInput objects -------------------------------
    chains: List[TCRChainInput] = []
    fmt = options.file_format.lower()

    if fmt == "fasta":
        chains, parse_warnings = parse_fasta(content_text, source="fasta")
        result.warnings.extend(parse_warnings)
    elif fmt in ("pdb", "mmcif"):
        chains, parse_warnings = parse_structure(
            content_text,
            file_format=fmt,
            molecule=options.molecule,
            tcr_chain_ids=options.tcr_chain_ids,
        )
        result.warnings.extend(parse_warnings)
    else:
        # Auto-detect from content
        guess = detect_file(filename, content_bytes)
        if guess.file_format == "fasta":
            chains, parse_warnings = parse_fasta(content_text, source="fasta")
        elif guess.file_format in ("pdb", "mmcif"):
            chains, parse_warnings = parse_structure(
                content_text,
                file_format=guess.file_format,
                molecule=options.molecule,
                tcr_chain_ids=options.tcr_chain_ids,
            )
        else:
            result.warnings.append(f"unrecognized file format: {guess.file_format}")
            return result
        result.warnings.extend(parse_warnings)
        result.file_format = guess.file_format

    if not chains:
        result.warnings.append("no TCR chains parsed from input")
        result.reference_info = get_reference_db().release_info
        return result

    # ----- Run gene mapping on each chain ---------------------------------
    for chain in chains:
        try:
            gene_call = map_chain(chain)
            result.chains.append(gene_call)
        except Exception as exc:  # noqa: BLE001
            result.warnings.append(
                f"chain {chain.chain_id}: gene mapping failed: {exc}"
            )
            result.chains.append(GeneCallResult(
                chain_id=chain.chain_id,
                chain_type=chain.inferred_chain_type or "unknown",
                notes=[f"gene mapping failed: {exc}"],
            ))

    # ----- Attach reference info ------------------------------------------
    result.reference_info = get_reference_db().release_info
    return result
