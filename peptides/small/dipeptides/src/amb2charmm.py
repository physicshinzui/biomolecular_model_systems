#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path


ACE_RENAME = {
    "CH3": "CAY",
    "C": "CY",
    "O": "OY",

    # AMBER variants
    "HH31": "HY1",
    "HH32": "HY2",
    "HH33": "HY3",
    "1HH3": "HY1",
    "2HH3": "HY2",
    "3HH3": "HY3",
    "H1": "HY1",
    "H2": "HY2",
    "H3": "HY3",
}

NME_RENAME = {
    "N": "NT",
    "H": "HNT",
    "CH3": "CAT",

    # AMBER variants
    "HH31": "HT1",
    "HH32": "HT2",
    "HH33": "HT3",
    "1HH3": "HT1",
    "2HH3": "HT2",
    "3HH3": "HT3",
    "H1": "HT1",
    "H2": "HT2",
    "H3": "HT3",
}

RESNAME_AMBER_TO_CHARMM = {
    "HID": "HSD",
    "HIE": "HSE",
    "HIP": "HSP",
    "CYS": "CYS",
    "CYX": "CYX",
    "ASP": "ASP",
    "ASH": "ASH",
    "GLU": "GLU",
    "GLH": "GLH",
    "LYS": "LYS",
    "LYN": "LSN",
}

GENERIC_RENAME = {
    "H": "HN",
}

RESIDUE_ATOM_RENAME = {
    "HSD": {
        "HD1": "HD1",
        "HE2": "HE2",
    },
    "HSE": {
        "HD1": "HD1",
        "HE2": "HE2",
    },
    "HSP": {
        "HD1": "HD1",
        "HE2": "HE2",
    },
    "CYS": {
        "HG": "HG1",
    },
    "CYX": {},
    "ASH": {
        "HD2": "HD2",
    },
    "GLH": {
        "HE2": "HE2",
    },
    "LSN": {
        "HZ1": "HZ1",
        "HZ2": "HZ2",
    },
}


def atom_name(line: str) -> str:
    return line[12:16].strip()


def residue_key(line: str) -> tuple[str, str, int, str]:
    return (
        line[17:20].strip(),
        line[21].strip() or "A",
        int(line[22:26]),
        line[26].strip(),
    )


def replace_fields(
    line: str,
    *,
    atom: str,
    resname: str,
    chain: str,
    resseq: int,
    serial: int,
) -> str:
    if len(atom) > 4:
        raise ValueError(f"atom name too long: {atom}")
    if len(resname) > 3:
        raise ValueError(f"residue name too long: {resname}")

    line = line.rstrip("\n")
    line = line[:6] + f"{serial:5d}" + line[11:]
    line = line[:12] + f"{atom:>4s}" + line[16:]
    line = line[:17] + f"{resname:>3s}" + line[20:]
    line = line[:21] + f"{chain:1s}" + line[22:]
    line = line[:22] + f"{resseq:4d}" + line[26:]
    return line + "\n"


def grouped_residues(atom_lines: list[str]) -> list[tuple[tuple[str, str, int, str], list[str]]]:
    out = []
    current_key = None
    current = []

    for line in atom_lines:
        key = residue_key(line)
        if key != current_key:
            if current:
                out.append((current_key, current))
            current_key = key
            current = [line]
        else:
            current.append(line)

    if current:
        out.append((current_key, current))

    return out


def convert_atom_name(resname: str, old: str) -> str:
    if old in RESIDUE_ATOM_RENAME.get(resname, {}):
        return RESIDUE_ATOM_RENAME[resname][old]
    return GENERIC_RENAME.get(old, old)


def convert(in_pdb: Path, out_pdb: Path) -> None:
    atom_lines = [
        line for line in in_pdb.read_text().splitlines(keepends=True)
        if line.startswith(("ATOM  ", "HETATM"))
    ]

    residues = grouped_residues(atom_lines)

    if len(residues) != 3:
        raise ValueError(f"expected exactly ACE-X-NME, found {len(residues)} residues")

    (ace_key, ace), (mid_key, mid), (nme_key, nme) = residues

    if ace_key[0] != "ACE":
        raise ValueError(f"first residue must be ACE, found {ace_key[0]}")
    if nme_key[0] not in {"NME", "NMA"}:
        raise ValueError(f"last residue must be NME/NMA, found {nme_key[0]}")

    amber_resname = mid_key[0]
    charmm_resname = RESNAME_AMBER_TO_CHARMM.get(amber_resname, amber_resname)

    if amber_resname == "CYX":
        raise ValueError(
            "CYX denotes disulphide-bonded cysteine. "
            "A single ACE-CYX-NME molecule has no disulphide partner; "
            "use CYS for thiol cysteine, or provide a real disulphide system."
        )

    chain = mid_key[1]
    resseq = 1
    serial = 1
    out = []

    for line in ace:
        old = atom_name(line)
        if old not in ACE_RENAME:
            raise ValueError(f"unexpected ACE atom: {old}")
        out.append(replace_fields(
            line, atom=ACE_RENAME[old], resname=charmm_resname,
            chain=chain, resseq=resseq, serial=serial
        ))
        serial += 1

    for line in mid:
        old = atom_name(line)
        new = convert_atom_name(charmm_resname, old)
        out.append(replace_fields(
            line, atom=new, resname=charmm_resname,
            chain=chain, resseq=resseq, serial=serial
        ))
        serial += 1

    for line in nme:
        old = atom_name(line)
        if old not in NME_RENAME:
            raise ValueError(f"unexpected NME atom: {old}")
        out.append(replace_fields(
            line, atom=NME_RENAME[old], resname=charmm_resname,
            chain=chain, resseq=resseq, serial=serial
        ))
        serial += 1

    out_pdb.write_text("".join(out) + "TER\nEND\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert AMBER ACE-X-NME PDB files to CHARMM/OpenMM-compatible PDB files."
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input PDB file or directory containing PDB files.",
    )
    parser.add_argument(
        "output",
        type=Path,
        help="Output PDB file or output directory.",
    )
    parser.add_argument(
        "--suffix",
        default="_charmm",
        help="Suffix for output files in directory mode. Default: _charmm",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )

    args = parser.parse_args()

    if args.input.is_file():
        if args.output.exists() and args.output.is_dir():
            out = args.output / f"{args.input.stem}{args.suffix}{args.input.suffix}"
        else:
            out = args.output

        if out.exists() and not args.overwrite:
            raise FileExistsError(f"output exists: {out}")

        convert(args.input, out)
        print(f"converted: {args.input} -> {out}")

    elif args.input.is_dir():
        args.output.mkdir(parents=True, exist_ok=True)

        pdb_files = sorted(args.input.glob("*.pdb"))

        if not pdb_files:
            raise FileNotFoundError(f"no .pdb files found in {args.input}")

        n_ok = 0
        n_fail = 0

        for in_pdb in pdb_files:
            out_pdb = args.output / f"{in_pdb.stem}{args.suffix}.pdb"

            if out_pdb.exists() and not args.overwrite:
                print(f"skip existing: {out_pdb}")
                continue

            try:
                convert(in_pdb, out_pdb)
                print(f"converted: {in_pdb} -> {out_pdb}")
                n_ok += 1
            except Exception as e:
                print(f"failed: {in_pdb}: {e}")
                n_fail += 1

        print(f"done: {n_ok} converted, {n_fail} failed")

    else:
        raise FileNotFoundError(f"input not found: {args.input}")
