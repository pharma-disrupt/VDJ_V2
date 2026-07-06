"""
gene_mapper.py
==============

The core scientific module: assigns V, D, J, C germline genes to a TCR chain
sequence and extracts the CDR3.

Strategy:
  1. chain-type disambiguation (alpha/beta/gamma/delta) by aligning the
     constant region of the input chain against one exemplar C-region per
     locus using Biopython's pure-Python `Bio.Align.PairwiseAligner` (local).
  2. V gene call: align the chain against the locus's V-gene panel; best hit
     by alignment score wins. Report identity % over the aligned region.
  3. J gene call: same against the J-gene panel.
  4. D gene call (beta/delta only): heuristic search for a short D-region
     motif between V and J alignments. Marked LOW confidence.
  5. C gene call: align the chain's C-terminal region against the locus's
     C-gene panel.
  6. CDR3 extraction: find the conserved Cys at the end of V (the 104-105 IMGT
     position) and the F/W-G-X-G motif at the start of J. Slice out CDR3.
"""

from __future__ import annotations

import re
from typing import Dict, List, Tuple, Optional, Any

from Bio.Align import PairwiseAligner, Alignment

from tcr_mapper.models import TCRChainInput, GeneCallResult
from tcr_mapper.reference_db import get_reference_db


# ---------------------------------------------------------------------------
# PairwiseAligner configuration (constructed once, reused)
# ---------------------------------------------------------------------------

_ALIGNER: Optional[PairwiseAligner] = None


def _get_aligner() -> PairwiseAligner:
    """Build (once) a local-aligner tuned for protein sequence matching."""
    global _ALIGNER
    if _ALIGNER is not None:
        return _ALIGNER

    a = PairwiseAligner()
    a.mode = "local"
    # Standard BLOSUM62-like scoring: match +2, mismatch -1, gap -5/-1
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -5.0
    a.extend_gap_score = -1.0
    _ALIGNER = a
    return a


# ---------------------------------------------------------------------------
# D-region mini-sequences (heuristic — used only for "best guess" D calling)
# ---------------------------------------------------------------------------
# These are the IMGT D-REGION amino-acid translations for the common human
# TRBD and TRDD genes. The D region in vivo is heavily trimmed and N-added,
# so protein-level D calling is unreliable. We use these for a best-effort
# annotation that is always labelled "low confidence".

_TRBD_SEQUENCES: Dict[str, str] = {
    "TRBD1*01": "GTGTGCCAGTGTGCC",
    "TRBD2*01": "GTTGTGAGTGTGCC",
}

_TRDD_SEQUENCES: Dict[str, str] = {
    "TRDD1*01": "GTGTTAAACTGGTACGAC",
    "TRDD2*01": "TTATGACACACAGGTGTGCC",
    "TRDD3*01": "GTGGTTACGCCACTATGCC",
}


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def _align(query: str, target: str) -> Alignment:
    """Run a local alignment and return the alignment object."""
    return _get_aligner().align(query, target)[0]


def _identity(aln: Alignment) -> float:
    """Compute identity % over the aligned region."""
    # alignment has .aligned: tuple of (query_ranges, target_ranges)
    q_range, t_range = aln.aligned
    if not q_range.size:
        return 0.0
    q_seq = str(aln[0]).replace("-", "")
    t_seq = str(aln[1]).replace("-", "")
    if not q_seq or not t_seq:
        return 0.0
    matches = sum(1 for a, b in zip(q_seq, t_seq) if a == b)
    return 100.0 * matches / max(len(q_seq), len(t_seq))


def _aligned_query_seq(aln: Alignment) -> str:
    return str(aln[0]).replace("-", "")


def _aligned_target_seq(aln: Alignment) -> str:
    return str(aln[1]).replace("-", "")


def _query_span(aln: Alignment) -> Tuple[int, int]:
    """Return (start, end) of the aligned region in the query sequence."""
    q_range, _t_range = aln.aligned
    if not q_range.size:
        return (0, 0)
    return (int(q_range[0][0]), int(q_range[0][1]))


# ---------------------------------------------------------------------------
# Chain-type disambiguation
# ---------------------------------------------------------------------------

def classify_chain_type(
    sequence: str,
    return_ranked: bool = False,
) -> Tuple[Optional[str], float, List[str]]:
    """
    Align the chain's C-terminal region against one constant-region exemplar
    per locus (alpha/beta/gamma/delta). Best hit = chain type.

    Args:
        sequence: amino-acid sequence of the chain (V+J+C, or full extracellular).
        return_ranked: if True, return a list of (locus, score) sorted desc
                       instead of a single best (locus, score).

    Returns:
        (chain_type, best_score, notes)
        OR if return_ranked=True:
        ([(locus, score), ...], 0.0, notes)
    """
    db = get_reference_db()
    notes: List[str] = []
    if not db.constant_exemplars:
        notes.append("no constant exemplars available; cannot classify chain type")
        if return_ranked:
            return ([], 0.0, notes)
        return (None, 0.0, notes)

    # Use the last ~120 residues (constant region) for chain-type calling.
    c_region = sequence[-120:] if len(sequence) > 120 else sequence

    scores: List[Tuple[str, float, float]] = []  # (locus, score, identity)
    for locus, (allele, ref_seq) in db.constant_exemplars.items():
        try:
            aln = _align(c_region, ref_seq)
            score = float(aln.score)
            ident = _identity(aln)
            scores.append((locus, score, ident))
        except Exception as exc:  # noqa: BLE001
            notes.append(f"alignment failed for locus {locus}: {exc}")
            continue

    if not scores:
        notes.append("all constant exemplar alignments failed")
        if return_ranked:
            return ([], 0.0, notes)
        return (None, 0.0, notes)

    scores.sort(key=lambda x: x[1], reverse=True)
    best_locus, best_score, best_ident = scores[0]

    if best_ident < 25.0:
        notes.append(
            f"best constant-region identity only {best_ident:.1f}% — "
            "chain may not be a TCR"
        )
        if return_ranked:
            return ([(s[0], s[1]) for s in scores], 0.0, notes)
        return (None, best_score, notes)

    notes.append(
        f"chain-type classified as {best_locus} "
        f"(C-region alignment score={best_score:.1f}, identity={best_ident:.1f}%)"
    )

    if return_ranked:
        return ([(s[0], s[1]) for s in scores], best_score, notes)
    return (best_locus, best_score, notes)


# ---------------------------------------------------------------------------
# V / J / C gene calling
# ---------------------------------------------------------------------------

def _best_hit_in_panel(
    query: str,
    panel: Dict[str, str],
    min_identity: float = 30.0,
) -> Optional[Tuple[str, float, float, Alignment]]:
    """
    Find the best-scoring allele in a panel.

    Returns (allele, score, identity, alignment) or None if no allele
    reaches min_identity.
    """
    best: Optional[Tuple[str, float, float, Alignment]] = None
    for allele, ref in panel.items():
        try:
            aln = _align(query, ref)
        except Exception:  # noqa: BLE001
            continue
        score = float(aln.score)
        ident = _identity(aln)
        if ident < min_identity:
            continue
        if best is None or score > best[1]:
            best = (allele, score, ident, aln)
    return best


def _call_v_gene(sequence: str, locus: str) -> Tuple[Optional[str], Optional[str], float, int, int, List[str]]:
    """
    Returns (allele, gene, identity_pct, q_start, q_end, notes).

    q_start/q_end is the span of the chain that aligned to the V reference
    (used to bound the CDR3 search).
    """
    db = get_reference_db()
    panel = db.get_v_panel(locus)
    notes: List[str] = []
    if not panel:
        notes.append(f"no V-gene panel for locus '{locus}'")
        return (None, None, 0.0, 0, 0, notes)

    # Search the N-terminal portion (V is at the start of the chain)
    query = sequence[: min(len(sequence), 230)]
    hit = _best_hit_in_panel(query, panel, min_identity=40.0)
    if hit is None:
        notes.append("no V gene reached 40% identity threshold")
        return (None, None, 0.0, 0, 0, notes)

    allele, _score, ident, aln = hit
    gene = db.v_gene_of_allele.get(allele, allele.split("*")[0])
    q_start, q_end = _query_span(aln)
    notes.append(f"V gene: {allele} (identity {ident:.1f}%, span [{q_start}:{q_end}])")
    return (allele, gene, ident, q_start, q_end, notes)


def _call_j_gene(sequence: str, locus: str, search_start: int = 0) -> Tuple[Optional[str], Optional[str], float, int, int, List[str]]:
    """
    Returns (allele, gene, identity_pct, q_start, q_end, notes).

    q_start/q_end is the span of the chain aligned to the J reference.
    """
    db = get_reference_db()
    panel = db.get_j_panel(locus)
    notes: List[str] = []
    if not panel:
        notes.append(f"no J-gene panel for locus '{locus}'")
        return (None, None, 0.0, 0, 0, notes)

    # J is located right after V (with CDR3 in between), so search from
    # v_end onward. Cap the search window to ~60 residues past v_end so we
    # don't drag in C-region matches that bias the alignment.
    if search_start:
        s = max(0, search_start - 5)  # small overlap to catch J's 5' anchor
    else:
        s = 0
    # J region + CDR3 is typically <= 40 aa total; give some headroom.
    e = min(len(sequence), s + 70)
    query = sequence[s:e]
    hit = _best_hit_in_panel(query, panel, min_identity=35.0)
    if hit is None:
        # Retry with a wider window (full chain from s onward)
        hit = _best_hit_in_panel(sequence[s:], panel, min_identity=35.0)
        if hit is None:
            # Last resort: search the whole sequence
            hit = _best_hit_in_panel(sequence, panel, min_identity=35.0)
            if hit is None:
                notes.append("no J gene reached 35% identity threshold")
                return (None, None, 0.0, 0, 0, notes)

    allele, _score, ident, aln = hit
    gene = db.j_gene_of_allele.get(allele, allele.split("*")[0])
    q_start, q_end = _query_span(aln)
    # Adjust to absolute coords if we sliced
    if s:
        q_start += s
        q_end += s
    notes.append(f"J gene: {allele} (identity {ident:.1f}%, span [{q_start}:{q_end}])")
    return (allele, gene, ident, q_start, q_end, notes)


def _call_c_gene(sequence: str, locus: str) -> Tuple[Optional[str], Optional[str], float, List[str]]:
    db = get_reference_db()
    panel = db.get_c_panel(locus)
    notes: List[str] = []
    if not panel:
        notes.append(f"no C-gene panel for locus '{locus}'")
        return (None, None, 0.0, notes)

    # Use the last ~140 residues
    query = sequence[-140:] if len(sequence) > 140 else sequence
    hit = _best_hit_in_panel(query, panel, min_identity=30.0)
    if hit is None:
        notes.append("no C gene reached 30% identity threshold")
        return (None, None, 0.0, notes)
    allele, _score, ident, _aln = hit
    gene = db.c_gene_of_allele.get(allele, allele.split("*")[0])
    notes.append(f"C gene: {allele} (identity {ident:.1f}%)")
    return (allele, gene, ident, notes)


# ---------------------------------------------------------------------------
# D gene calling (beta / delta only, low confidence)
# ---------------------------------------------------------------------------

def _call_d_gene(
    sequence: str,
    v_end: int,
    j_start: int,
    locus: str,
) -> Tuple[Optional[str], Optional[str], float, str, List[str]]:
    """
    Heuristic D-gene call. Returns (allele, gene, identity, confidence, notes).

    For beta/delta chains, search the substring between v_end and j_start for
    a short D-region motif. Always labelled "best guess, low confidence".
    """
    notes: List[str] = []
    if locus not in ("beta", "delta"):
        return (None, None, 0.0, "", notes)

    if v_end >= j_start or j_start - v_end < 3:
        notes.append("no room for D region between V and J alignments")
        return (None, None, 0.0, "", notes)

    junction = sequence[v_end:j_start]
    panel = _TRBD_SEQUENCES if locus == "beta" else _TRDD_SEQUENCES

    if not panel:
        return (None, None, 0.0, "", notes)

    best_allele: Optional[str] = None
    best_score = -1.0
    best_gene: Optional[str] = None
    for allele, ref in panel.items():
        try:
            aln = _align(junction, ref)
            score = float(aln.score)
            if score > best_score:
                best_score = score
                best_allele = allele
                best_gene = allele.split("*")[0]
        except Exception:  # noqa: BLE001
            continue

    if best_allele is None:
        return (None, None, 0.0, "", notes)

    # Compute identity over the aligned region
    try:
        aln = _align(junction, panel[best_allele])
        ident = _identity(aln)
    except Exception:  # noqa: BLE001
        ident = 0.0

    notes.append(
        f"D gene (best guess): {best_allele} "
        f"(identity {ident:.1f}% over junction of length {len(junction)})"
    )
    return (best_allele, best_gene, ident, "best guess, low confidence", notes)


# ---------------------------------------------------------------------------
# CDR3 extraction
# ---------------------------------------------------------------------------

# Conserved motifs:
#   - End of V: a Cys (C) immediately before CDR3 (IMGT position ~104)
#   - Start of J: F/W - G - X - G  (the FGXG / WGXG motif)
_J_MOTIF_REGEX = re.compile(r"[FW]G.XG")


def _extract_cdr3(
    sequence: str,
    v_end: int,
    j_start: int,
) -> Tuple[Optional[str], List[str]]:
    """
    Extract the CDR3 sequence.

    Rules:
      - The CDR3 starts at the conserved Cys near the end of V (the last C
        before position v_end). If no Cys is found in the last 30 residues
        of the V alignment, fall back to v_end as the start.
      - The CDR3 ends at the F/W-G-X-G motif at the start of J. If found,
        the CDR3 ends just BEFORE the F/W. If not found, fall back to
        j_start as the end.

    The CDR3 INCLUDES the conserved Cys and EXCLUDES the F/W-G-X-G motif
    (per IMGT definition).
    """
    notes: List[str] = []
    if v_end <= 0 or j_start <= 0 or v_end >= j_start:
        # Fall back to scanning the whole sequence for the Cys + J motif
        return _cdr3_scan_whole(sequence, notes)

    # Find the conserved Cys within the V-aligned region.
    # Search the last 30 residues up to and including v_end.
    search_start = max(0, v_end - 30)
    v_window = sequence[search_start:v_end + 5]  # +5 in case Cys is just past v_end
    cys_pos_in_window = v_window.rfind("C")
    if cys_pos_in_window < 0:
        notes.append("conserved Cys at end of V not found; using v_end as CDR3 start")
        cdr3_start = v_end
    else:
        cdr3_start = search_start + cys_pos_in_window

    # Find the F/W-G-X-G motif at the start of J.
    # Search a window starting a bit before j_start.
    j_window_start = max(0, j_start - 10)
    j_window_end = min(len(sequence), j_start + 30)
    j_window = sequence[j_window_start:j_window_end]
    m = _J_MOTIF_REGEX.search(j_window)
    if m:
        cdr3_end = j_window_start + m.start()  # exclude the motif itself
    else:
        notes.append("F/W-G-X-G motif at start of J not found; using j_start as CDR3 end")
        cdr3_end = j_start

    if cdr3_end <= cdr3_start:
        notes.append(f"CDR3 bounds invalid (start={cdr3_start}, end={cdr3_end}); cannot extract")
        return (None, notes)

    cdr3 = sequence[cdr3_start:cdr3_end]
    notes.append(f"CDR3 extracted: {cdr3} (length {len(cdr3)})")
    return (cdr3, notes)


def _cdr3_scan_whole(sequence: str, notes: List[str]) -> Tuple[Optional[str], List[str]]:
    """Fallback: scan the whole sequence for Cys ... F/W-G-X-G."""
    # Find the last C in the first half
    half = len(sequence) // 2
    cys_pos = sequence.rfind("C", 0, half + 10)
    if cys_pos < 0:
        notes.append("no Cys found in first half of sequence; cannot extract CDR3")
        return (None, notes)

    m = _J_MOTIF_REGEX.search(sequence, cys_pos)
    if m is None:
        notes.append("no F/W-G-X-G motif after the Cys; cannot extract CDR3")
        return (None, notes)

    cdr3 = sequence[cys_pos:m.start()]
    notes.append(f"CDR3 (whole-seq fallback): {cdr3} (length {len(cdr3)})")
    return (cdr3, notes)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def map_chain(chain: TCRChainInput) -> GeneCallResult:
    """
    Map one TCR chain to its germline V/D/J/C genes and CDR3.
    """
    notes: List[str] = []
    sequence = chain.sequence.upper()

    # 1. Determine chain type (use inferred_chain_type if provided by the parser)
    chain_type = chain.inferred_chain_type
    if chain_type is None or chain_type == "unknown":
        ctype, _score, ctype_notes = classify_chain_type(sequence)
        notes.extend(ctype_notes)
        if ctype is None:
            return GeneCallResult(
                chain_id=chain.chain_id,
                chain_type="unknown",
                notes=notes + ["chain type could not be determined; skipping gene mapping"],
            )
        chain_type = ctype

    result = GeneCallResult(chain_id=chain.chain_id, chain_type=chain_type)

    # 2. V gene
    v_allele, v_gene, v_ident, v_start, v_end, v_notes = _call_v_gene(sequence, chain_type)
    notes.extend(v_notes)
    result.v_allele = v_allele
    result.v_gene = v_gene
    result.v_identity_pct = round(v_ident, 1)

    # 3. J gene
    j_allele, j_gene, j_ident, j_start, j_end, j_notes = _call_j_gene(sequence, chain_type, search_start=v_end)
    notes.extend(j_notes)
    result.j_allele = j_allele
    result.j_gene = j_gene
    result.j_identity_pct = round(j_ident, 1)

    # 4. D gene (beta/delta only)
    if chain_type in ("beta", "delta") and v_end and j_start:
        d_allele, d_gene, d_ident, d_conf, d_notes = _call_d_gene(sequence, v_end, j_start, chain_type)
        notes.extend(d_notes)
        result.d_allele = d_allele
        result.d_gene = d_gene
        result.d_identity_pct = round(d_ident, 1)
        result.d_confidence = d_conf

    # 5. C gene
    c_allele, c_gene, c_ident, c_notes = _call_c_gene(sequence, chain_type)
    notes.extend(c_notes)
    result.c_allele = c_allele
    result.c_gene = c_gene
    result.c_identity_pct = round(c_ident, 1)

    # 6. CDR3
    cdr3, cdr3_notes = _extract_cdr3(sequence, v_end, j_start)
    notes.extend(cdr3_notes)
    result.cdr3_sequence = cdr3

    result.notes = notes
    return result
