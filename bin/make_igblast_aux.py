#!/usr/bin/env python3
# Written by Ayelet Peres and released under the MIT license.
"""Rebuild IgBLAST's germline-specific auxiliary files with receptor_utils.

IgBLAST ships internal_data/<species>/<species>.ndm.imgt (V FWR/CDR boundaries)
and optional_file/<species>_gl.aux (J reading frame, conserved PHE/TRP). Both
encode coordinates of the germline set they were built from, so a custom
reference needs its own or annotation is wrong.

Reads <reference>/<species>/vdj/*_<CHAIN>.fasta, which sourcerer writes
IMGT-gapped -- igblast_base/fasta/ cannot be used instead, as those are
gap-stripped and aggregated per locus class.
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from receptor_utils import aux_formats
from receptor_utils import simple_bio_seq as simple

V_CODE = {"IGHV": "VH", "IGKV": "VK", "IGLV": "VL",
          "TRAV": "VA", "TRBV": "VB", "TRDV": "VD", "TRGV": "VG"}
J_CHAINS = ["IGHJ", "IGKJ", "IGLJ", "TRAJ", "TRBJ", "TRDJ", "TRGJ"]


def allele_name(header):
    # IgBLAST matches ndm/aux records to the database by name, and those are built
    # from bare allele names; aux_formats also derives chain_type from chars 2-3.
    name = header.split()[0]
    return name.split("|")[1] if "|" in name else name


def cleaned(ref, dest):
    simple.write_fasta(str(dest), {allele_name(k): v
                                   for k, v in simple.read_fasta(str(ref)).items()})
    return dest


def parse(path):
    # Both tools write one "#" header; aux output is CRLF. Concatenating verbatim
    # would repeat the header.
    if not path.exists():
        return "", []
    lines = path.read_text().replace("\r", "").splitlines()
    header = next((l for l in lines if l.startswith("#")), "")
    return header, [l for l in lines if l.strip() and not l.startswith("#")]


def build(vdj, work, cdr_coords):
    ndm_header, aux_header, ndm, aux = "", "", [], []
    fasta, out = work / "in.fasta", work / "out"
    for chain, code in V_CODE.items():
        for ref in sorted(vdj.glob(f"*_{chain}.fasta")):
            out.unlink(missing_ok=True)
            aux_formats.ndm_from_fasta(str(cleaned(ref, fasta)), str(out), code, cdr_coords)
            header, rows = parse(out)
            ndm_header, ndm = ndm_header or header, ndm + rows
    for chain in J_CHAINS:
        for ref in sorted(vdj.glob(f"*_{chain}.fasta")):
            out.unlink(missing_ok=True)
            aux_formats.aux_from_fasta(str(cleaned(ref, fasta)), str(out), False)
            header, rows = parse(out)
            aux_header, aux = aux_header or header, aux + rows
    return ndm_header, ndm, aux_header, aux


def annotation_db(out, species):
    # This ships as NCBI's IMGT mirror. An allele missing from it gets no chain
    # type, and IgBLAST then skips the D and J search entirely.
    for locus_class, name in (("ig", f"{species}_V"), ("tr", f"{species}_TR_V")):
        src = out / "fasta" / f"{species}_{locus_class}_v.fasta"
        if not src.is_file() or src.stat().st_size == 0:
            continue
        dest = out / "internal_data" / species / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["makeblastdb", "-parse_seqids", "-dbtype", "nucl",
                        "-in", str(src), "-out", str(dest)],
                       check=True, stdout=subprocess.DEVNULL)
        print(f"{species}: rebuilt {name}")


def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(([header] if header else []) + rows) + "\n")


def self_test():
    import tempfile
    assert allele_name("J00256|IGHJ1*01|Homo_sapiens|F|J-REGION|723..774|52") == "IGHJ1*01"
    assert allele_name("IGHJ1*01") == "IGHJ1*01"
    assert allele_name("IGHV1-18*01 nt|1|") == "IGHV1-18*01"
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "x.ndm"
        src.write_text("#header\r\nIGHV1-2*01\t1\t2\r\n\r\nIGHV1-3*01\t3\t4\r\n")
        header, rows = parse(src)
        assert header == "#header", header
        assert rows == ["IGHV1-2*01\t1\t2", "IGHV1-3*01\t3\t4"], rows
        assert parse(Path(tmp) / "missing") == ("", [])
        out = Path(tmp) / "d" / "o.aux"
        write(out, header, rows)
        assert out.read_text() == "#header\nIGHV1-2*01\t1\t2\nIGHV1-3*01\t3\t4\n"
    print("self-test ok")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--igblast", type=Path, help="igblast_base to copy and update")
    p.add_argument("-r", "--reference", type=Path, help="reference_base with <species>/vdj")
    p.add_argument("-o", "--out", type=Path, default=Path("igblast_base"))
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()

    if a.self_test:
        return self_test()
    if not a.igblast or not a.reference:
        p.error("-i and -r are required")
    cdr_coords = [79, 114, 166, 195, 313]   # IMGT numbering

    shutil.copytree(a.igblast, a.out, symlinks=True, dirs_exist_ok=True)
    work = a.out / ".aux_work"
    work.mkdir(exist_ok=True)

    species = [d for d in sorted(a.reference.iterdir()) if (d / "vdj").is_dir()]
    if not species:
        sys.exit(f"No <species>/vdj directories under {a.reference}")

    for d in species:
        ndm_header, ndm, aux_header, aux = build(d / "vdj", work, cdr_coords)
        if not ndm or not aux:
            sys.exit(f"{d.name}: derived {len(ndm)} V delineations and {len(aux)} J "
                     f"annotations from {d / 'vdj'}. Are the V references IMGT-gapped?")
        write(a.out / "internal_data" / d.name / f"{d.name}.ndm.imgt", ndm_header, ndm)
        write(a.out / "optional_file" / f"{d.name}_gl.aux", aux_header, aux)
        annotation_db(a.out, d.name)
        print(f"{d.name}: {len(ndm)} ndm records, {len(aux)} aux records")

    shutil.rmtree(work)


if __name__ == "__main__":
    main()
