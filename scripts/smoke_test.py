"""
Smoke test: build a synthetic TCR-beta chain from known germline alleles,
run the full pipeline on it, and print a structured result so we can verify
V/J/C calls and CDR3 extraction look right.
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tcr_mapper.pipeline import run_pipeline
from tcr_mapper.reference_db import get_reference_db
from tcr_mapper.models import ProcessOptions


def main():
    db = get_reference_db()

    # Pick known alleles to construct a realistic test chain
    v_allele = "TRBV2*01"
    j_allele = "TRBJ2-1*01"
    c_allele = "TRBC1*01"

    v_seq = db.v_genes["beta"][v_allele]
    j_seq = db.j_genes["beta"][j_allele]
    c_seq = db.c_genes["beta"][c_allele]

    # Insert a recognizable CDR3 between V and J
    # Real CDR3 starts with C (conserved end of V) and ends before F/W-G-X-G (start of J)
    cdr3 = "CASSSYEQY"  # canonical beta CDR3 pattern
    fake_chain = v_seq + cdr3 + j_seq + c_seq

    fasta_content = f">test_beta_chain\n{fake_chain}\n"
    options = ProcessOptions(file_format="fasta", molecule="fasta")
    result = run_pipeline("test_beta.fasta", fasta_content.encode(), options)

    print("=" * 70)
    print("TCR Germline Mapper — smoke test")
    print("=" * 70)
    print(f"Input alleles (ground truth):")
    print(f"  V: {v_allele}")
    print(f"  J: {j_allele}")
    print(f"  C: {c_allele}")
    print(f"  CDR3 inserted: {cdr3}")
    print(f"  Full chain length: {len(fake_chain)} aa")
    print()
    print(f"Pipeline result: {len(result.chains)} chain(s)")
    for chain in result.chains:
        print(f"\n--- Chain: {chain.chain_id} ---")
        print(f"  chain_type:    {chain.chain_type}")
        print(f"  V gene/allele: {chain.v_gene} / {chain.v_allele}  (identity {chain.v_identity_pct}%)")
        print(f"  J gene/allele: {chain.j_gene} / {chain.j_allele}  (identity {chain.j_identity_pct}%)")
        print(f"  D gene/allele: {chain.d_gene or '-'} / {chain.d_allele or '-'}  "
              f"(identity {chain.d_identity_pct}%, conf: {chain.d_confidence or '-'})")
        print(f"  C gene/allele: {chain.c_gene} / {chain.c_allele}  (identity {chain.c_identity_pct}%)")
        print(f"  CDR3:          {chain.cdr3_sequence}")
        if chain.notes:
            print("  Notes:")
            for n in chain.notes:
                print(f"    - {n}")

    if result.warnings:
        print("\nWarnings:")
        for w in result.warnings:
            print(f"  - {w}")

    print()
    print("Reference DB info:")
    print(json.dumps(result.reference_info, indent=2))

    # Verify expectations
    print()
    print("=" * 70)
    print("Verification:")
    print("=" * 70)
    chain = result.chains[0]
    checks = [
        ("chain_type == beta", chain.chain_type == "beta"),
        ("V allele detected", chain.v_allele is not None),
        ("V allele matches expected", chain.v_allele == v_allele),
        ("J allele detected", chain.j_allele is not None),
        ("J allele matches expected", chain.j_allele == j_allele),
        ("C allele detected", chain.c_allele is not None),
        ("C allele matches expected", chain.c_allele == c_allele),
        ("CDR3 extracted", chain.cdr3_sequence is not None),
        ("CDR3 starts with C (conserved)", chain.cdr3_sequence and chain.cdr3_sequence.startswith("C")),
        ("CDR3 contains inserted motif",
         chain.cdr3_sequence and "SYEQY" in chain.cdr3_sequence),
        ("D confidence labelled low", chain.d_confidence == "best guess, low confidence"),
    ]
    all_pass = True
    for label, ok in checks:
        marker = "✓" if ok else "✗"
        print(f"  {marker} {label}")
        if not ok:
            all_pass = False

    print()
    print("Overall:", "ALL CHECKS PASSED" if all_pass else "SOME CHECKS FAILED")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
