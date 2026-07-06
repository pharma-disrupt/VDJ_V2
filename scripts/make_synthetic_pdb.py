"""
Generate a tiny synthetic PDB file with two TCR-like chains (D and E) for
testing the pMHC-TCR complex parsing path. Each chain is just a short poly-
alanine + a few key residues — enough for the parser to detect chain length
and assign chain IDs, but NOT enough for real gene calling (that's fine;
this is a parser test, not a biology test).
"""

from pathlib import Path

OUT = Path("/home/z/my-project/tests/fixtures/synthetic_tcr.pdb")

PDB_TEMPLATE = """\
HEADER    T CELL RECEPTOR                  01-JAN-25   XXXX
TITLE     SYNTHETIC TCR FOR PARSER TESTING
COMPND    MOL_ID: 1;
COMPND   2 MOLECULE: TCR BETA CHAIN;
COMPND   3 CHAIN: D;
COMPND   4 MOL_ID: 2;
COMPND   5 MOLECULE: TCR ALPHA CHAIN;
COMPND   6 CHAIN: E
SOURCE    MOL_ID: 1;
SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS
SEQRES   1 D  200  ALA ALA ALA ALA ALA ALA ALA ALA ALA
SEQRES   1 E  200  ALA ALA ALA ALA ALA ALA ALA ALA ALA
ATOM      1  N   ALA D   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA D   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA D   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA D   1       1.251   2.390   0.000  1.00 20.00           O
{atoms_d}
ATOM  {nplus1:5d}  N   ALA E   1      50.000   0.000   0.000  1.00 20.00           N
ATOM  {nplus2:5d}  CA  ALA E   1      51.458   0.000   0.000  1.00 20.00           C
ATOM  {nplus3:5d}  C   ALA E   1      52.009   1.420   0.000  1.00 20.00           C
ATOM  {nplus4:5d}  O   ALA E   1      51.251   2.390   0.000  1.00 20.00           O
{atoms_e}
TER
END
"""

def make_chain_atoms(chain_id: str, start_serial: int, n_residues: int = 200):
    """Generate ATOM records for a chain of n ALA residues."""
    lines = []
    serial = start_serial
    # Simplified backbone: N, CA, C, O per residue
    for i in range(1, n_residues + 1):
        x = (i - 1) * 3.84 + (50 if chain_id == "E" else 0)
        lines.append(
            f"ATOM  {serial:5d}  N   ALA {chain_id}{i:4d}    {x:8.3f}   0.000   0.000  1.00 20.00           N"
        )
        serial += 1
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA {chain_id}{i:4d}    {x+1.458:8.3f}   0.000   0.000  1.00 20.00           C"
        )
        serial += 1
        lines.append(
            f"ATOM  {serial:5d}  C   ALA {chain_id}{i:4d}    {x+2.009:8.3f}   1.420   0.000  1.00 20.00           C"
        )
        serial += 1
        lines.append(
            f"ATOM  {serial:5d}  O   ALA {chain_id}{i:4d}    {x+1.251:8.3f}   2.390   0.000  1.00 20.00           O"
        )
        serial += 1
    return "\n".join(lines), serial


atoms_d, next_serial = make_chain_atoms("D", start_serial=5, n_residues=200)
atoms_e, _ = make_chain_atoms("E", start_serial=next_serial + 4, n_residues=200)

content = PDB_TEMPLATE.format(
    atoms_d=atoms_d,
    atoms_e=atoms_e,
    nplus1=next_serial,
    nplus2=next_serial + 1,
    nplus3=next_serial + 2,
    nplus4=next_serial + 3,
)

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(content)
print(f"Wrote {OUT} ({len(content)} bytes)")
print(f"Chain D: 200 residues (poly-ALA)")
print(f"Chain E: 200 residues (poly-ALA)")
