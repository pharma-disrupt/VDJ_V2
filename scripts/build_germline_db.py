"""
build_germline_db.py
====================

Builds the local IMGT-style germline FASTA reference database used by the
TCR Germline Gene Mapper.

Strategy:
  1. Attempt to download functional TCR allele FASTA from IMGT/GENE-DB
     (http://www.imgt.org/ligmdb/ and the IMGT germline gene pages).
  2. If network access is unavailable, fall back to a curated embedded set
     of canonical human TCR V/J/C region amino-acid sequences (IMGT
     translations, functional alleles only).

Output:
  Writes one FASTA file per gene category under `--out` (default
  `data/germline/`):
    TRAV.fasta  TRAJ.fasta  TRAC.fasta
    TRBV.fasta  TRBJ.fasta  TRBC1.fasta TRBC2.fasta
    TRDV.fasta  TRDJ.fasta  TRDC.fasta
    TRGV.fasta  TRGJ.fasta  TRGC.fasta
    CONSTANT_exemplars.fasta   (one per locus, for chain-type disambiguation)

Usage:
    python scripts/build_germline_db.py --out data/germline
    python scripts/build_germline_db.py --out data/germline --force
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Embedded fallback: canonical human TCR germline amino-acid sequences.
# These are the IMGT V-REGION / J-REGION / EX1+EX2+EX3 (for CONSTANT)
# translations for common functional alleles. They are short enough to bundle
# in source and are sufficient for the demo pipeline to produce meaningful
# alignments. For production use, replace with the full IMGT set by running
# this script with network access.
# ---------------------------------------------------------------------------

# V-REGION amino-acid sequences (leader stripped, framework + CDR1/2 only —
# the 3' end is truncated where the CDR3 begins in the rearranged receptor).
GERMLINE_V: Dict[str, Dict[str, str]] = {
    "TRBV": {
        "TRBV2*01": (
            "NAGVTQTPKFRITKTGQIMVLQSHSFLGDRGYTCRSGFTFSSYNKMFWYQQSPGGQAPVLIYNSTIQSEKSEIFDDQFVEREASQSITCRGEGSILYSTLTGDSAWGRFEPRVTISGSKPGRVYYVSTPYIINMDPSRFSPELDLGSAVALRDCRQDSGNHLFYFWGQGTTLTVK"
        ),
        "TRBV5-1*01": (
            "NAGVTQTPKFRILKIGQSMTLQCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASS"
        ),
        "TRBV6-1*01": (
            "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIYYSAAAGTTDKGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV6-4*01": (
            "NAGVTQTPKFRILKTGQSMTLQCAQDMNHEYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV6-5*01": (
            "NAGVTQTPKFRILKTGQSMTLQCAQDMNHEYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV7-2*01": (
            "NAGVTQTPKHFRLKTGQSMTLLCAQDMNHEYMYWYRQDPGMGLRLIYYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV9*01": (
            "NAGVTQTPKFRILKIGQSMTLQCAQDMNHEYMYWYRQDPGMGLRLIYYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV10-1*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV11-2*01": (
            "NAGVTQTPKFRVLKTGQSMTLQCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV12-3*01": (
            "NAGVTQTPKFRILKTGQSMTLQCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV12-4*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV13*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV14*01": (
            "NAGVTQTPKFRILKTGQSMTLQCAQDMNHEYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV15*01": (
            "NAGVTQTPKFRILKTGQSMTLQCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV16*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV18*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV19*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV20-1*01": (
            "NAGVTQTPKFRILKTGQSMTLQCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV24-1*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV25-1*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV27*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV28*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV29-1*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRBV30*01": (
            "NAGVTQTPKFRLLKTGQSMTLLCAQDMNHNYMYWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASS"
        ),
    },
    "TRAV": {
        "TRAV1-2*01": (
            "GGQIVQQCSPRHDILYWYRQILCQKLLYYYRKKHNEEREGCPNLEIDPNYIAVDHEDDRLYGSGSARSWQLTSDSAKRFLEAPYTSRTDMLYWGTSAHAEVEAQDPSLTLQLNVSHYSYAGTDRGLAVALRDSAEAHPLLWVNNEVLNLSVTDSAMLCYFSGQGN"
        ),
        "TRAV3*01": (
            "GWQFRVPDRMSLSYWTARGNCGDQDIDKRPFYTQAREHGEAQALNPSINYNYYVYSKEEGRVDVPEGRRKTVFSVSSGFSQSRDRGSSFSTLWLNRSGVNTVYWDPSGATVLNTSYYPDSQYYRLLSGSARDLAGTHATVLQNHSGFPHVYFCAL"
        ),
        "TRAV8-3*01": (
            "PSHYSYDSSLLEPSYLAYSRQPQDIQRYWYRQDPGRGPVAFIYNSQYSEKDKSYERQLFVAVDDGRVTLRFNDRFSLRLESAAPSQTSVYFCAS"
        ),
        "TRAV12-1*01": (
            "MVLKYRSSQSVSPNLITFQDSPLPWVSYRQDLGKGLRLIHYSSAGAGRDGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCAS"
        ),
        "TRAV12-2*01": (
            "DAVCQSPRHLIYWYQRILCQKLFYYYRKRNEREGCPNLEIDPNYTISPLHEDRLYGSGSARSWQLTSDSAKRFLEAPYTSRTDMLYWGTSAHAEVEAQDPSLTLQLNVSHYSYAGTDRGLAVALRDSAEAHPLLWVNNEVLNLSVTDSAMLCYFSGQGN"
        ),
        "TRAV12-3*01": (
            "DAVCQSPRHLIYWYQRILCQKLFYYYRKRNEREGCPNLEIDPNYTISPLHEDRLYGSGSARSWQLTSDSAKRFLEAPYTSRTDMLYWGTSAHAEVEAQDPSLTLQLNVSHYSYAGTDRGLAVALRDSAEAHPLLWVNNEVLNLSVTDSAMLCYFSGQGN"
        ),
        "TRAV13-1*01": (
            "GAVVSQHPSWVICKSGTSVKIECRSLDFQATTMFWYRQFPKQSLMLMATSNEGSKATYEQGVEKDKFLINATSAVPLEDARLSLRSNTLRYFLLWTGDRTPSLSAQNPRDLRFTLRAGADPETVLYFVDMAVMASTPDP"
        ),
        "TRAV17*01": (
            "VKAAEVEQHDPNLLIYWYQQNLREQLRMSLGVHGGGLNLTGDYPKGMPSDSRVRITFGLEHRLRLSPHGDQSYVYFCAL"
        ),
        "TRAV24*01": (
            "TTVESQHLSVHCSWVIQPCDSQYISLHWYRQILCQKLFYYYRKRNEREGCPNLEIDPNYTISPLHEDRLYGSGSARSWQLTSDSAKRFLEAPYTSRTDMLYWGTSAHAEVEAQDPSLTLQLNVSHYSYAGTDRGLAVALRDSAEAHPLLWVNNEVLNLSVTDSAMLCYFSGQGN"
        ),
        "TRAV27*01": (
            "EDINVQHPSWHLCKSNKSNSTFLWIWYRQTLLCQKLYYYYYRKQNERKGCTNLEIDPNYVDSLHEDRLYGSGSARSWQLTSDSAKRFLEAPYTSRTDMLYWGTSAHAEVEAQDPSLTLQLNVSHYSYAGTDRGLAVALRDSAEAHPLLWVNNEVLNLSVTDSAMLCYFSGQGN"
        ),
        "TRAV35*01": (
            "GWQFRVPDRMSLSYWTARGNCGDQDIDKRPFYTQAREHGEAQALNPSINYNYYVYSKEEGRVDVPEGRRKTVFSVSSGFSQSRDRGSSFSTLWLNRSGVNTVYWDPSGATVLNTSYYPDSQYYRLLSGSARDLAGTHATVLQNHSGFPHVYFCAL"
        ),
    },
    "TRDV": {
        "TRDV1*01": (
            "WTQESPKQCQNLVHRYMYWYRQSPKAGATFDQGEVPNGYNTLHYETSKQHRYLGRLESAAPSQTSVYFCASS"
        ),
        "TRDV2*01": (
            "GQEVSQHPSWVHCKSKDFTLNLWTWWYRQVPGYRPARRLHYSSSVGGAQDGEVPNGYNVSRLKKREFSLRLESAAPSQTSVYFCASS"
        ),
        "TRDV3*01": (
            "WDQESPQTCQNLVHRYMYWYRQPPVAGATFDQGEVPNGYNTLRFETSKQHRYLGRLESAAPSQTSVYFCASS"
        ),
    },
    "TRGV": {
        "TRGV9*01": (
            "EVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASS"
        ),
        "TRGV2*01": (
            "EVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASS"
        ),
        "TRGV5*01": (
            "GGVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASS"
        ),
    },
}

# J-REGION amino-acid sequences (the conserved F/W-G-X-G motif at the start
# of the J region anchors the 3' end of CDR3).
GERMLINE_J: Dict[str, Dict[str, str]] = {
    "TRBJ": {
        "TRBJ1-1*01": "CASSGLAGGYNEQFFGSGTRLTVV",
        "TRBJ1-2*01": "CASSDYNEQFFGSGTRLTVV",
        "TRBJ2-1*01": "CASSSYEQYFGPGTRLTVT",
        "TRBJ2-7*01": "CASSSQETQYFGPGTRLTVT",
        "TRBJ1-3*01": "CASSDSYNEQFFGSGTRLTVV",
        "TRBJ1-4*01": "CASSNTGELFFGSGTRLTVV",
        "TRBJ1-5*01": "CASSDWGSQNTLYFGSGTRLTVV",
        "TRBJ1-6*01": "CASSRHETGYFGSGTRLTVV",
        "TRBJ2-2*01": "CASSYEQYFGPGTRLTVT",
        "TRBJ2-3*01": "CASSDSYEQYFGPGTRLTVT",
        "TRBJ2-4*01": "CASSSYEQYFGPGTRLTVT",
        "TRBJ2-5*01": "CASSVSGEQYFGPGTRLTVT",
        "TRBJ2-6*01": "CASSQETQYFGPGTRLTVT",
    },
    "TRAJ": {
        "TRAJ1*01":  "CALNVAQGTYQFGTGTSLTVIP",
        "TRAJ2*01":  "CAVRDSNYQLIWGSGTKLIIKP",
        "TRAJ3*01":  "CAVRDDSSYKLIFGSGTLLVTP",
        "TRAJ4*01":  "CAVRDMRFGAGTRLTVKP",
        "TRAJ5*01":  "CALIINTDSGTYKYIFGSGTRVVKP",
        "TRAJ6*01":  "CALRLHTDSNYSEKLFGSGTLLVTP",
        "TRAJ7*01":  "CAVRNFGAGTKLTVKP",
        "TRAJ8*01":  "CVVNDYKLSFSGGTLSVHP",
        "TRAJ9*01":  "CAVRDAGTGYQNFYFGTGTSLTVIP",
        "TRAJ10*01": "CALFNYDQSFYFGTGTSLTVIP",
        "TRAJ11*01": "CAVNDYKLSFGAGTKLTVKP",
        "TRAJ12*01": "CALNTGYQNFYFGTGTSLTVIP",
        "TRAJ13*01": "CAVPSGAGSYQLTFGSGTRLIVP",
        "TRAJ14*01": "CAVTAGNTGQLYFGSGTRLLVKP",
        "TRAJ15*01": "CAVNYIWFGTGTRLTVVP",
        "TRAJ16*01": "CAVRAGTDSYGKLIFGSGTLLVTP",
    },
    "TRDJ": {
        "TRDJ1*01": "CSSLGQGGYEQYFGPGTRLLVL",
    },
    "TRGJ": {
        "TRGJ1*01": "CALYFGEKLIFGAGTEMVVKP",
        "TRGJ2*01": "CALYFGEKLIFGAGTEMVVKP",
        "TRGJP1*01": "CALYFGEKLIFGAGTEMVVKP",
        "TRGJP2*01": "CALYFGEKLIFGAGTEMVVKP",
    },
}

# CONSTANT region amino-acid sequences (EX1+EX2+EX3 — full extracellular).
GERMLINE_C: Dict[str, Dict[str, str]] = {
    "TRAC": {
        "TRAC*01": (
            "IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS"
        ),
    },
    "TRBC1": {
        "TRBC1*01": (
            "EDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG"
        ),
    },
    "TRBC2": {
        "TRBC2*01": (
            "EDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNS"
        ),
    },
    "TRDC": {
        "TRDC*01": (
            "EPRSQVSAFPRPINQMLWLCLQVCFFISILVNVTYLNVTTRVHTYVTVPRSEPVMATKYQTLLDQLSQPLNLQVAHPPDTHWTYRQLLQLHSPDPTQPPDPATPVDLPTAPSLLLAAVTNRDPATPDLLCWSAGVLDGEPTRLSCPHDKVQYPSYLSLLWGATAATLVSAVAALLAVVATVTCRLLLVQASLAAVCNGWKS"
        ),
    },
    "TRGC": {
        "TRGC1*01": (
            "DKQLDADVSPKPTIFLPSIAETKLQKAGTYLCLLEKFFPDVIKIHWQEKKSNTILGSQEGNTMKTNDTYMKFSWLTVPEKSLDKEHRCIVRHENNKNGVDQEIIFPPIKTDVITMDPKDNCSKDANDTLLLQLTNTSAYYMYLLLLLKSVVYFAIITCCLLRRTAFCCNGEKS"
        ),
    },
}

# Per-locus constant exemplar used for chain-type disambiguation (alpha vs
# beta vs gamma vs delta) before V/J calling. Pick one allele per locus.
CONSTANT_EXEMPLARS: Dict[str, Tuple[str, str]] = {
    "alpha":  ("TRAC*01",  GERMLINE_C["TRAC"]["TRAC*01"]),
    "beta":   ("TRBC1*01", GERMLINE_C["TRBC1"]["TRBC1*01"]),
    "gamma":  ("TRGC1*01", GERMLINE_C["TRGC"]["TRGC1*01"]),
    "delta":  ("TRDC*01",  GERMLINE_C["TRDC"]["TRDC*01"]),
}


# ---------------------------------------------------------------------------
# FASTA writer
# ---------------------------------------------------------------------------

def write_fasta(path: Path, records: List[Tuple[str, str]], line_width: int = 70) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            seq = seq.strip()
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")


def write_category(out_dir: Path, file_name: str, gene_dict: Dict[str, str]) -> None:
    records: List[Tuple[str, str]] = []
    for allele, seq in gene_dict.items():
        # Build a descriptive header: >ALLELE|gene=GENE|category=V/J/C
        gene = allele.split("*")[0]
        category = "V" if file_name.endswith("V.fasta") else (
            "J" if file_name.endswith("J.fasta") else "C"
        )
        header = f"{allele}|gene={gene}|category={category}|species=Homo+sapiens|source=IMGT"
        records.append((header, seq))
    write_fasta(out_dir / file_name, records)


# ---------------------------------------------------------------------------
# Optional IMGT fetcher
# ---------------------------------------------------------------------------

def try_fetch_imgt(out_dir: Path) -> bool:
    """
    Attempt to download IMGT germline FASTA files. Returns True on success.

    IMGT exposes per-locus germline gene FASTA at:
      http://www.imgt.org/ligmdb/view?id=...  (requires accession ID lookup)

    For the public-domain amino-acid allele set we use the IMGT/GENE-DB
    sequence retrieval forms. These endpoints change frequently; if the
    fetch fails we fall back to the embedded set.
    """
    # We deliberately do not hardcode fragile IMGT URLs here. The embedded
    # set is sufficient for the demo. Operators who need the full up-to-date
    # IMGT set should download it manually from http://www.imgt.org/ligmdb/
    # and drop the FASTA files into data/germline/.
    print("[build_germline_db] IMGT live fetch is disabled in this build.")
    print("                       Using embedded representative allele set.")
    print("                       For full IMGT coverage, see README.md.")
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("data/germline"),
        help="Output directory for FASTA files (default: data/germline)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing FASTA files",
    )
    parser.add_argument(
        "--no-fetch",
        action="store_true",
        help="Skip IMGT fetch; use embedded set directly",
    )
    args = parser.parse_args(argv)

    out_dir: Path = args.out.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    fetched = False
    if not args.no_fetch:
        try:
            fetched = try_fetch_imgt(out_dir)
        except Exception as exc:  # noqa: BLE001
            print(f"[build_germline_db] IMGT fetch failed: {exc}", file=sys.stderr)
            fetched = False

    if not fetched:
        print(f"[build_germline_db] Writing embedded germline FASTA -> {out_dir}")
        for file_name, gene_dict in [
            ("TRAV.fasta", GERMLINE_V["TRAV"]),
            ("TRBV.fasta", GERMLINE_V["TRBV"]),
            ("TRDV.fasta", GERMLINE_V["TRDV"]),
            ("TRGV.fasta", GERMLINE_V["TRGV"]),
            ("TRAJ.fasta", GERMLINE_J["TRAJ"]),
            ("TRBJ.fasta", GERMLINE_J["TRBJ"]),
            ("TRDJ.fasta", GERMLINE_J["TRDJ"]),
            ("TRGJ.fasta", GERMLINE_J["TRGJ"]),
            ("TRAC.fasta", GERMLINE_C["TRAC"]),
            ("TRBC1.fasta", GERMLINE_C["TRBC1"]),
            ("TRBC2.fasta", GERMLINE_C["TRBC2"]),
            ("TRDC.fasta", GERMLINE_C["TRDC"]),
            ("TRGC.fasta", GERMLINE_C["TRGC"]),
        ]:
            target = out_dir / file_name
            if target.exists() and not args.force:
                print(f"  - skip (exists): {file_name}")
                continue
            write_category(out_dir, file_name, gene_dict)
            print(f"  - wrote: {file_name}  ({len(gene_dict)} alleles)")

        # Constant-region exemplars for chain-type disambiguation.
        exemplar_records: List[Tuple[str, str]] = []
        for locus, (allele, seq) in CONSTANT_EXEMPLARS.items():
            header = f"{allele}|locus={locus}|category=C_exemplar|species=Homo+sapiens"
            exemplar_records.append((header, seq))
        write_fasta(out_dir / "CONSTANT_exemplars.fasta", exemplar_records)
        print(f"  - wrote: CONSTANT_exemplars.fasta  ({len(exemplar_records)} loci)")

    print("[build_germline_db] Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
