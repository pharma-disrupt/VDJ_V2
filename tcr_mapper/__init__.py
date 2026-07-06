"""
tcr_mapper
==========

Pure-Python TCR germline gene assignment core.

This package has ZERO dependency on FastAPI / Mangum / any web framework.
It can be imported and unit-tested on its own.

Public API:
    from tcr_mapper.pipeline import run_pipeline
    from tcr_mapper.models import TCRChainInput, GeneCallResult, PipelineResult
"""

from tcr_mapper.models import (
    TCRChainInput,
    GeneCallResult,
    PipelineResult,
    FileTypeGuess,
    ProcessOptions,
)
from tcr_mapper.pipeline import run_pipeline

__all__ = [
    "TCRChainInput",
    "GeneCallResult",
    "PipelineResult",
    "FileTypeGuess",
    "ProcessOptions",
    "run_pipeline",
]

__version__ = "0.1.0"
