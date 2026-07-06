"""
Dataclasses that flow through the TCR mapper pipeline.

These are intentionally plain `@dataclass` objects (not pydantic models) so
the processing core has no hard dependency on pydantic. The web layer wraps
these into pydantic models where serialization is needed.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any


@dataclass
class TCRChainInput:
    """Normalized internal representation of one TCR chain ready for mapping."""

    chain_id: str
    sequence: str
    source: str  # "pdb_tcr_only" | "pdb_pmhc_complex" | "fasta"
    inferred_chain_type: Optional[str] = None  # "alpha" | "beta" | "gamma" | "delta" | None


@dataclass
class GeneCallResult:
    """Result of germline gene assignment for one TCR chain."""

    chain_id: str
    chain_type: str  # "alpha" | "beta" | "gamma" | "delta"
    v_gene: Optional[str] = None
    v_allele: Optional[str] = None
    v_identity_pct: float = 0.0
    j_gene: Optional[str] = None
    j_allele: Optional[str] = None
    j_identity_pct: float = 0.0
    d_gene: Optional[str] = None
    d_allele: Optional[str] = None
    d_identity_pct: float = 0.0
    d_confidence: Optional[str] = None  # "best guess, low confidence" | None
    c_gene: Optional[str] = None
    c_allele: Optional[str] = None
    c_identity_pct: float = 0.0
    cdr3_sequence: Optional[str] = None
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class FileTypeGuess:
    """Best-guess at file format/molecule content, returned by /api/upload."""

    file_format: str  # "pdb" | "mmcif" | "fasta" | "unknown"
    suggested_molecule: str  # "tcr_only" | "pmhc_tcr_complex" | "fasta" | "unknown"
    sniffed_notes: List[str] = field(default_factory=list)
    filename: str = ""
    size_bytes: int = 0


@dataclass
class ProcessOptions:
    """Confirmed options from the user, sent to /api/process."""

    file_format: str = "fasta"           # "pdb" | "mmcif" | "fasta"
    molecule: str = "tcr_only"           # "tcr_only" | "pmhc_tcr_complex" | "auto" (PDB only)
    tcr_chain_ids: Optional[List[str]] = None  # explicit chain IDs, e.g. ["D", "E"]


@dataclass
class PipelineResult:
    """Top-level result of the full pipeline."""

    chains: List[GeneCallResult] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    file_format: str = ""
    molecule: str = ""
    reference_info: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "chains": [c.to_dict() for c in self.chains],
            "warnings": self.warnings,
            "file_format": self.file_format,
            "molecule": self.molecule,
            "reference_info": self.reference_info,
        }
