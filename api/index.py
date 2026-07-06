"""
api/index.py
============

FastAPI entry point for the TCR Germline Gene Mapper.

Wrapped with Mangum so it runs as a Vercel serverless Python function
(also compatible with AWS Lambda). Routes:

  POST /api/upload            - upload file, get back best-guess type for confirmation
  POST /api/process           - run the full pipeline given file + confirmed options
  GET  /api/health            - liveness check
  GET  /api/reference-info    - metadata about the bundled IMGT germline DB
"""

from __future__ import annotations

import base64
import json
import os
import time
from typing import Optional, List

from fastapi import FastAPI, File, Form, UploadFile, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, HTMLResponse
from mangum import Mangum
from pydantic import BaseModel, Field

# Allow importing tcr_mapper whether this file is run from /api/ (Vercel) or
# from the project root (uvicorn).
import sys
from pathlib import Path
_THIS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _THIS_DIR.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from tcr_mapper.file_detect import detect_file
from tcr_mapper.models import ProcessOptions
from tcr_mapper.pipeline import run_pipeline
from tcr_mapper.reference_db import get_reference_db


# ---------------------------------------------------------------------------
# App setup
# ---------------------------------------------------------------------------

app = FastAPI(
    title="TCR Germline Gene Mapper",
    version="0.1.0",
    description="Map T-cell receptor chain sequences to their germline V(D)J/C gene origins.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# Vercel Hobby tier body size cap (slightly below 4.5MB to be safe).
MAX_BODY_BYTES = 4_500_000


# ---------------------------------------------------------------------------
# Pydantic models for request/response
# ---------------------------------------------------------------------------

class UploadResponse(BaseModel):
    filename: str
    size_bytes: int
    file_format: str
    suggested_molecule: str
    sniffed_notes: List[str] = Field(default_factory=list)


class ProcessRequest(BaseModel):
    """
    The process request can be sent either:
      - as multipart/form-data (file + form fields) — preferred for browsers
      - as application/json with a base64-encoded `file_b64` — for clients
        that can't easily send multipart (e.g. some SDKs)
    """
    filename: str
    file_b64: Optional[str] = None
    file_format: str = "fasta"
    molecule: str = "fasta"
    tcr_chain_ids: Optional[List[str]] = None


class ProcessResponse(BaseModel):
    chains: List[dict]
    warnings: List[str]
    file_format: str
    molecule: str
    reference_info: dict
    elapsed_ms: int


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.exception_handler(HTTPException)
async def http_exc_handler(request: Request, exc: HTTPException):
    return JSONResponse(
        status_code=exc.status_code,
        content={"error": exc.detail},
    )


@app.get("/api/health")
async def health():
    """Liveness check."""
    return {
        "status": "ok",
        "service": "tcr-germline-gene-mapper",
        "version": "0.1.0",
        "timestamp": int(time.time()),
    }


@app.get("/api/reference-info")
async def reference_info():
    """Metadata about the bundled IMGT germline DB."""
    db = get_reference_db()
    return db.release_info


@app.post("/api/upload", response_model=UploadResponse)
async def upload(file: UploadFile = File(...)):
    """
    Step 1: upload a file and get back a best-guess at its type so the
    frontend can ask the user to confirm format + molecule content before
    /api/process runs the heavy pipeline.
    """
    content = await file.read()
    if len(content) > MAX_BODY_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"file too large ({len(content)} bytes); max {MAX_BODY_BYTES} bytes on Vercel Hobby tier",
        )

    guess = detect_file(file.filename or "upload.bin", content)
    return UploadResponse(
        filename=guess.filename,
        size_bytes=guess.size_bytes,
        file_format=guess.file_format,
        suggested_molecule=guess.suggested_molecule,
        sniffed_notes=guess.sniffed_notes,
    )


@app.post("/api/process", response_model=ProcessResponse)
async def process(
    file: UploadFile = File(...),
    file_format: str = Form("fasta"),
    molecule: str = Form("fasta"),
    tcr_chain_ids: Optional[str] = Form(None),
):
    """
    Step 2: run the full pipeline. Expects multipart/form-data with the file
    and the user's confirmed options.
    """
    t0 = time.time()
    content = await file.read()
    if len(content) > MAX_BODY_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"file too large ({len(content)} bytes); max {MAX_BODY_BYTES} bytes on Vercel Hobby tier",
        )

    # Parse tcr_chain_ids (comma-separated string -> list)
    chain_ids: Optional[List[str]] = None
    if tcr_chain_ids:
        chain_ids = [c.strip() for c in tcr_chain_ids.split(",") if c.strip()]

    options = ProcessOptions(
        file_format=file_format,
        molecule=molecule,
        tcr_chain_ids=chain_ids,
    )

    result = run_pipeline(file.filename or "upload.bin", content, options)

    elapsed_ms = int((time.time() - t0) * 1000)
    return ProcessResponse(
        chains=[c.to_dict() for c in result.chains],
        warnings=result.warnings,
        file_format=result.file_format,
        molecule=result.molecule,
        reference_info=result.reference_info,
        elapsed_ms=elapsed_ms,
    )


@app.post("/api/process-json", response_model=ProcessResponse)
async def process_json(req: ProcessRequest):
    """
    Alternative JSON-based entry point (file_b64 encoded). Useful for
    programmatic clients.
    """
    t0 = time.time()
    if not req.file_b64:
        raise HTTPException(status_code=400, detail="file_b64 is required")

    try:
        content = base64.b64decode(req.file_b64)
    except Exception:
        raise HTTPException(status_code=400, detail="file_b64 is not valid base64")

    if len(content) > MAX_BODY_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"file too large ({len(content)} bytes); max {MAX_BODY_BYTES} bytes",
        )

    options = ProcessOptions(
        file_format=req.file_format,
        molecule=req.molecule,
        tcr_chain_ids=req.tcr_chain_ids,
    )
    result = run_pipeline(req.filename, content, options)
    elapsed_ms = int((time.time() - t0) * 1000)
    return ProcessResponse(
        chains=[c.to_dict() for c in result.chains],
        warnings=result.warnings,
        file_format=result.file_format,
        molecule=result.molecule,
        reference_info=result.reference_info,
        elapsed_ms=elapsed_ms,
    )


# Fallback: serve a minimal HTML page at /api so hitting the function root
# shows something useful rather than 404.
@app.get("/api", response_class=HTMLResponse)
async def api_root():
    return """
    <html><body>
    <h1>TCR Germline Gene Mapper API</h1>
    <p>Available endpoints:</p>
    <ul>
      <li>POST /api/upload — upload a file, get back a best-guess type</li>
      <li>POST /api/process — run the full pipeline (multipart)</li>
      <li>POST /api/process-json — run the full pipeline (JSON+base64)</li>
      <li>GET /api/health — liveness check</li>
      <li>GET /api/reference-info — bundled IMGT germline DB metadata</li>
    </ul>
    <p><a href="/">Go to the app</a></p>
    </body></html>
    """


# ---------------------------------------------------------------------------
# Mangum handler (Vercel / AWS Lambda entry point)
# ---------------------------------------------------------------------------

handler = Mangum(app, lifespan="off")
