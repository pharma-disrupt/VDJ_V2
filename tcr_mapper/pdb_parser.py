"""
pdb_parser.py
=============

Parses PDB / mmCIF structures into TCRChainInput objects.

Two cases handled:
  Case A — TCR-only PDB: parse all chains, extract AA sequences, classify
           each chain as alpha/beta/gamma/delta via constant-region alignment.
  Case B — pMHC-TCR complex PDB: by default pick chains D & E (IMGT/STCRDab
           convention) as the TCR chains. Validate by length heuristic
           (TCR ~180-215 aa, MHC ~180-275 aa, peptide ~8-15 aa). Fall back
           to scanning all chains if D/E don't fit; surface a warning.
"""

from __future__ import annotations

import io
import warnings as _warnings
from typing import List, Tuple, Dict, Optional

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqUtils import seq1

from tcr_mapper.models import TCRChainInput
from tcr_mapper.gene_mapper import classify_chain_type


# ---------------------------------------------------------------------------
# Length heuristics (residue counts)
# ---------------------------------------------------------------------------

_TCR_CHAIN_LEN_RANGE = (170, 230)   # covers both V+C of one TCR chain
_MHC_CHAIN_LEN_RANGE = (160, 300)
_PEPTIDE_LEN_RANGE = (6, 25)
_TCR_FULL_LEN_RANGE = (180, 230)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _extract_chain_sequences(structure, use_seqres: bool = False) -> Dict[str, str]:
    """
    Return {chain_id: amino_acid_sequence} for all chains in the structure.

    Uses ATOM records (coordinates present). Biopython's PDBParser does not
    expose SEQRES by default — we use ATOM records as the source of truth.
    """
    chain_seqs: Dict[str, str] = {}
    for model in structure:
        for chain in model:
            residues = []
            for residue in chain:
                # Skip waters and hetero atoms
                hetflag = residue.id[0]
                if hetflag.strip() != "":
                    continue
                resname = residue.get_resname()
                if not resname:
                    continue
                aa = seq1(resname, undef_code="X")
                residues.append(aa)
            if residues:
                chain_seqs[chain.id] = "".join(residues)
        break  # only first model
    return chain_seqs


def parse_structure(
    content: str,
    file_format: str,
    molecule: str,
    tcr_chain_ids: Optional[List[str]] = None,
) -> Tuple[List[TCRChainInput], List[str]]:
    """
    Parse a PDB/mmCIF string into TCRChainInput objects.

    Args:
        content: the file content as a string.
        file_format: "pdb" or "mmcif".
        molecule: "tcr_only" | "pmhc_tcr_complex" | "auto".
        tcr_chain_ids: explicit list of chain IDs to extract (overrides D/E default).

    Returns:
        (list of TCRChainInput, list of warnings/notes)
    """
    warnings_list: List[str] = []

    # Silence Biopython's noisy PDBConstructionWarning about missing atoms
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore", PDBConstructionWarning)

        if file_format == "mmcif":
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("tcr", io.StringIO(content))
        else:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("tcr", io.StringIO(content))

    chain_seqs = _extract_chain_sequences(structure, use_seqres=False)
    if not chain_seqs:
        return [], ["no chains with amino-acid residues found in structure"]

    warnings_list.append(
        f"structure has {len(chain_seqs)} chain(s): {sorted(chain_seqs.keys())}"
    )

    # ------------------------------------------------------ explicit override
    if tcr_chain_ids:
        selected = {}
        for cid in tcr_chain_ids:
            if cid in chain_seqs:
                selected[cid] = chain_seqs[cid]
            else:
                warnings_list.append(
                    f"requested chain '{cid}' not found in structure; available: {sorted(chain_seqs.keys())}"
                )
        if selected:
            return _finalize_chains(selected, "pdb_pmhc_complex" if molecule != "tcr_only" else "pdb_tcr_only", warnings_list)
        # else fall through to default logic

    # ----------------------------------------------------------- dispatch
    if molecule == "tcr_only":
        return _case_tcr_only(chain_seqs, warnings_list)
    elif molecule == "pmhc_tcr_complex":
        return _case_pmhc_complex(chain_seqs, warnings_list)
    elif molecule == "auto":
        return _case_auto(chain_seqs, warnings_list)
    else:
        # Default: try TCR-only mode
        warnings_list.append(f"unknown molecule mode '{molecule}'; defaulting to tcr_only")
        return _case_tcr_only(chain_seqs, warnings_list)


# ---------------------------------------------------------------------------
# Case A: TCR-only PDB
# ---------------------------------------------------------------------------

def _case_tcr_only(
    chain_seqs: Dict[str, str],
    warnings_list: List[str],
) -> Tuple[List[TCRChainInput], List[str]]:
    """
    Treat all chains as candidate TCR chains. Classify each via constant-region
    alignment; keep ones that look like TCR (alpha/beta/gamma/delta).
    """
    chains: List[TCRChainInput] = []
    for cid, seq in chain_seqs.items():
        # Heuristic length filter: TCR chains are 170-230 aa after V+J+C
        if not (_TCR_CHAIN_LEN_RANGE[0] <= len(seq) <= _TCR_CHAIN_LEN_RANGE[1] + 30):
            warnings_list.append(
                f"chain {cid} length {len(seq)} outside TCR range ({_TCR_CHAIN_LEN_RANGE}); "
                "still attempting classification"
            )

        ctype, score, notes = classify_chain_type(seq)
        warnings_list.extend([f"chain {cid}: {n}" for n in notes])

        if ctype is None:
            warnings_list.append(
                f"chain {cid} could not be confidently classified as TCR; skipped"
            )
            continue

        chains.append(TCRChainInput(
            chain_id=cid,
            sequence=seq,
            source="pdb_tcr_only",
            inferred_chain_type=ctype,
        ))

    if not chains:
        warnings_list.append("no chains classified as TCR in tcr_only mode")
    return chains, warnings_list


# ---------------------------------------------------------------------------
# Case B: pMHC-TCR complex PDB
# ---------------------------------------------------------------------------

def _case_pmhc_complex(
    chain_seqs: Dict[str, str],
    warnings_list: List[str],
) -> Tuple[List[TCRChainInput], List[str]]:
    """
    Default assumption: chains D & E are TCR (IMGT/STCRDab convention).
    Validate by length. Fall back to scanning all chains if D/E don't fit.
    """
    default_ids = ["D", "E"]
    selected: Dict[str, str] = {}
    for cid in default_ids:
        if cid in chain_seqs:
            seq = chain_seqs[cid]
            if _TCR_CHAIN_LEN_RANGE[0] <= len(seq) <= _TCR_CHAIN_LEN_RANGE[1] + 30:
                selected[cid] = seq
            else:
                warnings_list.append(
                    f"chain {cid} length {len(seq)} outside TCR range "
                    f"({_TCR_CHAIN_LEN_RANGE}); not auto-selected"
                )

    if len(selected) < 2:
        warnings_list.append(
            "default D/E chains did not yield 2 TCR chains; scanning all chains"
        )
        return _scan_all_chains_for_tcr(chain_seqs, warnings_list)

    return _finalize_chains(selected, "pdb_pmhc_complex", warnings_list)


# ---------------------------------------------------------------------------
# Case auto
# ---------------------------------------------------------------------------

def _case_auto(
    chain_seqs: Dict[str, str],
    warnings_list: List[str],
) -> Tuple[List[TCRChainInput], List[str]]:
    """
    Auto-detect: try D/E first (pMHC convention); if that yields 2 TCR chains,
    treat as pMHC-TCR complex. Otherwise treat as TCR-only.
    """
    default_ids = ["D", "E"]
    selected: Dict[str, str] = {}
    for cid in default_ids:
        if cid in chain_seqs:
            seq = chain_seqs[cid]
            if _TCR_CHAIN_LEN_RANGE[0] <= len(seq) <= _TCR_CHAIN_LEN_RANGE[1] + 30:
                selected[cid] = seq

    if len(selected) == 2:
        warnings_list.append("auto-detect: D/E chains look like TCR; treating as pMHC-TCR complex")
        return _finalize_chains(selected, "pdb_pmhc_complex", warnings_list)

    # Otherwise, look at all chains — if exactly 2 fit TCR length, use them
    tcr_like = {
        cid: seq for cid, seq in chain_seqs.items()
        if _TCR_CHAIN_LEN_RANGE[0] <= len(seq) <= _TCR_CHAIN_LEN_RANGE[1] + 30
    }
    if len(tcr_like) >= 2:
        warnings_list.append(
            f"auto-detect: {len(tcr_like)} chains fit TCR length range; "
            "taking the two longest as TCR alpha/beta"
        )
        # Take two longest
        sorted_chains = sorted(tcr_like.items(), key=lambda kv: len(kv[1]), reverse=True)
        picked = dict(sorted_chains[:2])
        return _finalize_chains(picked, "pdb_tcr_only", warnings_list)

    warnings_list.append("auto-detect: falling back to TCR-only mode on all chains")
    return _case_tcr_only(chain_seqs, warnings_list)


# ---------------------------------------------------------------------------
# Fallback: scan all chains for TCR-like sequences
# ---------------------------------------------------------------------------

def _scan_all_chains_for_tcr(
    chain_seqs: Dict[str, str],
    warnings_list: List[str],
) -> Tuple[List[TCRChainInput], List[str]]:
    """Pick the 2 chains whose length best fits TCR (180-230 aa)."""
    candidates = {
        cid: seq for cid, seq in chain_seqs.items()
        if _TCR_CHAIN_LEN_RANGE[0] <= len(seq) <= _TCR_CHAIN_LEN_RANGE[1] + 30
    }
    if not candidates:
        warnings_list.append(
            "no chains fit TCR length range; cannot extract TCR from complex"
        )
        return [], warnings_list

    # Take two longest
    sorted_chains = sorted(candidates.items(), key=lambda kv: len(kv[1]), reverse=True)
    picked = dict(sorted_chains[:2])
    warnings_list.append(
        f"selected chains {list(picked.keys())} as TCR (longest TCR-like)"
    )
    return _finalize_chains(picked, "pdb_pmhc_complex", warnings_list)


# ---------------------------------------------------------------------------
# Finalize: classify each selected chain
# ---------------------------------------------------------------------------

def _finalize_chains(
    selected: Dict[str, str],
    source: str,
    warnings_list: List[str],
) -> Tuple[List[TCRChainInput], List[str]]:
    """Attach chain-type (alpha/beta/gamma/delta) to each selected chain."""
    chains: List[TCRChainInput] = []
    used_types: List[str] = []

    # Try to assign distinct chain types (alpha vs beta) so we don't get two alphas
    # from a D/E pair. Try each chain with its best type, but if both pick the
    # same type, force the second to its second-best.
    ranked: List[Tuple[str, List[Tuple[str, float]]]] = []
    for cid, seq in selected.items():
        ctype, score, notes = classify_chain_type(seq, return_ranked=True)
        warnings_list.extend([f"chain {cid}: {n}" for n in notes])
        if isinstance(ctype, list):
            ranked.append((cid, ctype))
        else:
            ranked.append((cid, [(ctype or "unknown", score)]))

    for cid, ranking in ranked:
        chosen_type = None
        for cand_type, _score in ranking:
            if cand_type and cand_type not in used_types:
                chosen_type = cand_type
                break
        if chosen_type is None and ranking:
            chosen_type = ranking[0][0]
        if chosen_type:
            used_types.append(chosen_type)

        chains.append(TCRChainInput(
            chain_id=cid,
            sequence=selected[cid],
            source=source,
            inferred_chain_type=chosen_type,
        ))

    return chains, warnings_list
