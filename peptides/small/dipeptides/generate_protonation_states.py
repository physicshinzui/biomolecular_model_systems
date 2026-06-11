#!/usr/bin/env python3
"""Generate Amber/OpenMM-compatible protonation-state PDB variants."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "protonation_states"


@dataclass(frozen=True)
class Atom:
    record: str
    serial: int
    name: str
    residue: str
    residue_id: int
    x: float
    y: float
    z: float
    element: str
    charge: str = ""

    @property
    def xyz(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z


def parse_pdb(path: Path) -> list[Atom]:
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            Atom(
                record=line[0:6].strip(),
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                residue=line[17:20].strip(),
                residue_id=int(line[22:26]),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                element=line[76:78].strip(),
                charge=line[78:80].strip(),
            )
        )
    return atoms


def atom_name_field(name: str, element: str) -> str:
    if len(name) == 4 or (name and name[0].isdigit()):
        return f"{name:<4}"
    if len(element) == 1:
        return f" {name:<3}"
    return f"{name:<4}"


def format_atom(atom: Atom) -> str:
    return (
        f"{atom.record:<6}{atom.serial:5d} "
        f"{atom_name_field(atom.name, atom.element)}"
        f" {atom.residue:>3}  {atom.residue_id:4d}    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}          "
        f"{atom.element:>2}{atom.charge:>2}"
    )


def format_conect(atoms: list[Atom]) -> list[str]:
    covalent_radii = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}
    bonds = []
    for index, first in enumerate(atoms):
        for second in atoms[index + 1 :]:
            if first.element == second.element == "H":
                continue
            distance = math.sqrt(
                sum(value * value for value in subtract(first.xyz, second.xyz))
            )
            maximum = covalent_radii[first.element] + covalent_radii[second.element] + 0.45
            if distance <= maximum:
                bonds.append((first.serial, second.serial))
    return [f"CONECT{first:5d}{second:5d}" for first, second in bonds]


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(a, factor):
    return tuple(x * factor for x in a)


def normalize(a):
    length = math.sqrt(sum(x * x for x in a))
    return scale(a, 1.0 / length)


def new_hydrogen(
    atoms: list[Atom],
    parent_name: str,
    neighbor_names: tuple[str, ...],
    hydrogen_name: str,
    bond_length: float,
) -> Atom:
    sidechain = {atom.name: atom for atom in atoms if atom.residue_id == 2}
    parent = sidechain[parent_name]
    direction = (0.0, 0.0, 0.0)
    for neighbor_name in neighbor_names:
        direction = add(
            direction, normalize(subtract(parent.xyz, sidechain[neighbor_name].xyz))
        )
    position = add(parent.xyz, scale(normalize(direction), bond_length))
    return Atom(
        record="ATOM",
        serial=0,
        name=hydrogen_name,
        residue=parent.residue,
        residue_id=parent.residue_id,
        x=position[0],
        y=position[1],
        z=position[2],
        element="H",
    )


def make_variant(
    source: str,
    residue_name: str,
    remove: set[str] | None = None,
    add_spec: tuple[str, tuple[str, ...], str, float] | None = None,
    charges: dict[str, str] | None = None,
) -> list[Atom]:
    atoms = parse_pdb(ROOT / source)
    remove = remove or set()
    charges = charges or {}
    result = []
    added = None
    if add_spec:
        added = new_hydrogen(atoms, *add_spec)

    for atom in atoms:
        if atom.residue_id == 2 and atom.name in remove:
            continue
        if atom.residue_id == 2:
            atom = replace(
                atom,
                residue=residue_name,
                charge=charges.get(atom.name, ""),
            )
        result.append(atom)
        if added and atom.residue_id == 2 and atom.name == add_spec[0]:
            result.append(replace(added, residue=residue_name))

    return [replace(atom, serial=index) for index, atom in enumerate(result, start=1)]


VARIANTS = {
    "ASH": dict(
        source="D.pdb",
        residue_name="ASH",
        add_spec=("OD2", ("CG", "OD1"), "HD2", 0.96),
    ),
    "GLH": dict(
        source="E.pdb",
        residue_name="GLH",
        add_spec=("OE2", ("CD", "OE1"), "HE2", 0.96),
    ),
    "HID": dict(
        source="H.pdb",
        residue_name="HID",
        remove={"HE2"},
        add_spec=("ND1", ("CG", "CE1"), "HD1", 1.01),
    ),
    "HIE": dict(source="H.pdb", residue_name="HIE"),
    "HIP": dict(
        source="H.pdb",
        residue_name="HIP",
        add_spec=("ND1", ("CG", "CE1"), "HD1", 1.01),
        charges={"ND1": "1+"},
    ),
    "LYN": dict(source="K.pdb", residue_name="LYN", remove={"1HZ"}),
    "CYM": dict(
        source="C.pdb",
        residue_name="CYM",
        remove={"HG"},
        charges={"SG": "1-"},
    ),
}


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    for name, specification in VARIANTS.items():
        atoms = make_variant(**specification)
        records = [format_atom(atom) for atom in atoms]
        records.extend(format_conect(atoms))
        text = "\n".join(records) + "\nEND\n"
        (OUTPUT_DIR / f"{name}.pdb").write_text(text)
        print(f"wrote {OUTPUT_DIR / f'{name}.pdb'} ({len(atoms)} atoms)")


if __name__ == "__main__":
    main()
