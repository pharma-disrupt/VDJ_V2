# TCR Germline Gene Mapper

A serverless web application that maps T-cell receptor (TCR) chain sequences to their
germline V, D, J, and C gene origins and extracts the inferred CDR3 loop.

Deployable 100% for free on Vercel. No OS-level binaries required — alignment is done
with Biopython's pure-Python `PairwiseAligner` instead of BLAST+.

---

## Features

- **Accepts three input types**
  - TCR-only PDB structure
  - pMHC–TCR complex PDB structure (auto-extracts TCR chains, defaults to IMGT/STCRDab `D`/`E` convention)
  - FASTA (amino acid or nucleotide; nucleotide is translated in-frame)
- **V(D)J germline assignment** against bundled IMGT reference alleles using local alignment
- **CDR3 extraction** via the conserved Cys (end of V) and F/W-G-X-G motif (start of J)
- **Two-step UX** — upload returns a best-guess type, then the user confirms before processing
- **Pure serverless** — FastAPI + Mangum, runs in a single Vercel Python function

---

## Architecture

```
public/ (HTML/JS/CSS, static)  ──►  /api/*  ──►  api/index.py (FastAPI + Mangum)
                                                    │
                                                    ▼
                                         tcr_mapper/ (pure-Python core)
                                            ├── file_detect
                                            ├── pdb_parser / fasta_parser
                                            ├── reference_db (in-memory IMGT dicts)
                                            ├── gene_mapper (PairwiseAligner)
                                            └── pipeline
```

The processing core (`tcr_mapper/`) has **zero dependency on the web layer** — it is
importable and unit-testable on its own.

---

## Local Development

```bash
# 1. Install dependencies (use a venv)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. Run the API
uvicorn api.index:app --reload --port 8000

# 3. Serve the frontend (in another terminal)
cd public && python -m http.server 5173
# Open http://localhost:5173
```

Alternatively, use Vercel's dev server which mirrors production routing:

```bash
npm i -g vercel
vercel dev
```

---

## Deployment to Vercel

1. Push this repo to GitHub.
2. In Vercel dashboard: **New Project → Import** the repo.
3. Framework preset: **Other**. Build command: none. Output directory: `public`.
4. Vercel auto-detects `api/index.py` as a Python serverless function via `@vercel/python`.
5. `vercel.json` is already configured — no further setup needed.
6. Deploy. The app is live at `https://<your-project>.vercel.app`.

---

## API Surface

| Method | Route                  | Purpose                                                       |
|--------|------------------------|---------------------------------------------------------------|
| POST   | `/api/upload`          | Upload a file, get back a best-guess at type for confirmation |
| POST   | `/api/process`         | Run the full pipeline given the file + confirmed options      |
| GET    | `/api/health`          | Liveness check                                                |
| GET    | `/api/reference-info`  | Metadata about the bundled IMGT germline DB                   |

---

## Reference Database

The `data/germline/` directory contains IMGT-style FASTA files of functional
TCR germline alleles (TRAV, TRBV, TRDV, TRGV, TRAJ, TRBJ, TRDJ, TRGJ, TRAC,
TRBC1, TRBC2, TRDC, TRGC).

To refresh from IMGT/GENE-DB:

```bash
python scripts/build_germline_db.py --out data/germline
```

The script attempts to download from IMGT's public FTP/HTTP interface. If
network access is unavailable, the bundled representative alleles are used
as a fallback so the pipeline remains functional.

---

## Known Limitations

- **D-gene calls** are inherently low-confidence from protein sequence alone
  (D regions are short and heavily trimmed/N-added in vivo). They are always
  labelled "best guess, low confidence" or omitted if identity is too low.
- **PDB gaps**: missing residues in crystal structures can truncate extracted
  sequences. `ATOM` records are used by default; cross-check against `SEQRES`
  if available.
- **Serverless payload limit**: Vercel Hobby tier caps request bodies at
  ~4.5 MB. Larger files are rejected with HTTP 413.
- **Cold starts**: the first request after inactivity takes ~2–3 s to load
  the Python runtime and germline DB into memory.
- **This is a research tool, not a clinical or diagnostic device.**

---

## License

MIT — see `LICENSE` (or use the code as you see fit; IMGT data is subject to
IMGT's own academic use terms, see <http://www.imgt.org>).
