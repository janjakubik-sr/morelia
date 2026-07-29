#!/usr/bin/env python3
"""
measure_residue_centroid_ligand_atom_distances.py

General PDB/mmCIF batch script for measuring distances from selected ligand
atoms to the centroid of any selected receptor residue.

Main use case:
  - choose one receptor residue by residue name, chain and residue number
  - calculate the centroid of that residue, normally from all heavy atoms
  - choose ligand atoms either by exact atom names or by element counts
    such as S:1,O:2, which is robust to CIF atom-name renumbering
  - write combined CSV files for all PDB/mmCIF structures

Outputs:
  <out-prefix>_long.csv
      one row per measured ligand atom
  <out-prefix>_summary.csv
      one row per ligand residue / structure, with one column per requested atom
  <out-prefix>_failures.csv
      written only if some structures could not be processed

Examples:
  python measure_residue_centroid_ligand_atom_distances.py \
      --input-glob '*.cif' \
      --ref-chain A --ref-resname TRP --ref-resseq 223 \
      --ligand-resname LIG1 --ligand-chain B --ligand-resseq 1 \
      --ligand-element-counts S:1,O:2 \
      --out-prefix trp223_ligand_centroid

  python measure_residue_centroid_ligand_atom_distances.py \
      --input-glob 'M1*/*.yob.pdb' \
      --ref-chain A --ref-resname TRP --ref-resseq 378 \
      --ligand-resname UNL --ligand-chain L --ligand-resseq 1 \
      --ligand-atoms S,O,O1 \
      --symmetric-groups O,O1 \
      --out-prefix trp378_ligand_centroid
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import math
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SCRIPT_VERSION = "2026-07-23-v9-wx-gui"

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASX", "GLX", "SEC", "PYL", "MSE",
}
WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
COMMON_IONS = {
    "AL", "BA", "BR", "CA", "CD", "CL", "CO", "CS", "CU", "F", "FE", "HG",
    "IOD", "K", "LI", "MG", "MN", "NA", "NI", "RB", "SR", "ZN",
}
COMMON_SOLVENTS_BUFFERS = {
    "ACE", "ACT", "BME", "DMS", "DMSO", "EDO", "FMT", "GOL", "HEP", "MES",
    "PEG", "PO4", "SO4", "TRS",
}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


@dataclass(frozen=True)
class Atom:
    record: str
    serial: int
    atom_name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    element: str

    @property
    def coord(self) -> Tuple[float, float, float]:
        return self.x, self.y, self.z

    @property
    def residue_key(self) -> Tuple[str, int, str, str]:
        return self.chain, self.resseq, self.icode, self.resname


@dataclass(frozen=True)
class SelectedLigandAtom:
    label: str
    selector: str
    atom: Optional[Atom]
    selection_mode: str


def infer_element(atom_name: str) -> str:
    name = atom_name.strip()
    name = re.sub(r"^[0-9]+", "", name)
    letters = re.sub(r"[^A-Za-z]", "", name)
    if not letters:
        return ""

    two = letters[:2].upper()
    if two in {
        "CL", "BR", "NA", "MG", "CA", "MN", "FE", "CO", "NI", "CU", "ZN",
        "CD", "HG", "LI", "AL", "SI",
    }:
        return two

    return letters[0].upper()


def open_structure_text(path: Path):
    """Open a plain-text or gzipped structure file."""
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def altloc_score(atom: Atom) -> Tuple[int, float]:
    if atom.altloc == "":
        return 3, atom.occupancy
    if atom.altloc == "A":
        return 2, atom.occupancy
    return 1, atom.occupancy


def collapse_altlocs(raw_atoms: List[Atom]) -> List[Atom]:
    """
    Collapse alternate locations per atom. Prefer blank altloc, then A,
    then highest occupancy.
    """
    best: Dict[Tuple[str, str, int, str, str, str], Atom] = {}
    for atom in raw_atoms:
        key = (
            atom.record,
            atom.chain,
            atom.resseq,
            atom.icode,
            atom.resname,
            atom.atom_name,
        )
        old = best.get(key)
        if old is None or altloc_score(atom) > altloc_score(old):
            best[key] = atom
    return list(best.values())


def parse_pdb_atom_line(line: str) -> Optional[Atom]:
    if not line.startswith(("ATOM  ", "HETATM")):
        return None

    try:
        atom_name = line[12:16].strip()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = infer_element(atom_name)

        occ_text = line[54:60].strip()
        occupancy = float(occ_text) if occ_text else 0.0

        serial_text = line[6:11].strip()
        serial = int(serial_text) if serial_text else 0

        return Atom(
            record=line[0:6].strip(),
            serial=serial,
            atom_name=atom_name,
            altloc=line[16:17].strip(),
            resname=line[17:20].strip().upper(),
            chain=line[21:22].strip(),
            resseq=int(line[22:26]),
            icode=line[26:27].strip(),
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
            occupancy=occupancy,
            element=element,
        )
    except Exception:
        return None


def read_pdb_first_model(path: Path) -> List[Atom]:
    raw_atoms: List[Atom] = []
    seen_model = False

    with open_structure_text(path) as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if seen_model:
                    break
                seen_model = True
                continue
            if line.startswith("ENDMDL") and seen_model:
                break

            atom = parse_pdb_atom_line(line)
            if atom is not None:
                raw_atoms.append(atom)

    return collapse_altlocs(raw_atoms)


def cif_missing(value: Optional[str]) -> bool:
    return value is None or value in {".", "?", ""}


def cif_get(row: Dict[str, str], tag: str) -> Optional[str]:
    return row.get(tag) or row.get("_atom_site." + tag)


def clean_cif_value(value: Optional[str], default: str = "") -> str:
    if cif_missing(value):
        return default
    return str(value)


def parse_cif_int(value: Optional[str], default: int = 0) -> int:
    if cif_missing(value):
        return default
    text = str(value)
    try:
        return int(text)
    except ValueError:
        try:
            return int(float(text))
        except ValueError:
            return default


def parse_cif_float(value: Optional[str], default: float = 0.0) -> float:
    if cif_missing(value):
        return default
    try:
        return float(str(value))
    except ValueError:
        return default


def tokenize_cif_line(line: str) -> List[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def iter_cif_atom_site_rows(path: Path) -> Iterable[Dict[str, str]]:
    """Yield dictionaries from the first _atom_site loop in an mmCIF file."""
    with open_structure_text(path) as handle:
        lines = list(handle)

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped != "loop_":
            i += 1
            continue

        i += 1
        tags: List[str] = []
        while i < len(lines):
            tag_line = lines[i].strip()
            if not tag_line.startswith("_"):
                break
            tags.append(tag_line.split()[0])
            i += 1

        if not tags or not all(tag.startswith("_atom_site.") for tag in tags):
            continue

        n_tags = len(tags)
        pending: List[str] = []
        while i < len(lines):
            line = lines[i]
            stripped = line.strip()

            if not stripped:
                i += 1
                continue
            if stripped == "#":
                break
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                i -= 1
                break
            if stripped.startswith(";"):
                i += 1
                while i < len(lines) and not lines[i].startswith(";"):
                    i += 1
                i += 1
                continue

            pending.extend(tokenize_cif_line(line))
            while len(pending) >= n_tags:
                values = pending[:n_tags]
                pending = pending[n_tags:]
                yield dict(zip(tags, values))
            i += 1
        return


def atom_from_cif_row(row: Dict[str, str], cif_id_mode: str = "auth") -> Optional[Atom]:
    group = clean_cif_value(cif_get(row, "group_PDB"), "ATOM").upper()
    if group not in {"ATOM", "HETATM"}:
        return None

    model_num = clean_cif_value(cif_get(row, "pdbx_PDB_model_num"), "1")
    if model_num not in {"1", ""}:
        return None

    if cif_id_mode == "label":
        chain = clean_cif_value(cif_get(row, "label_asym_id"))
        resseq = parse_cif_int(cif_get(row, "label_seq_id"), default=0)
        resname = clean_cif_value(cif_get(row, "label_comp_id")).upper()
        atom_name = clean_cif_value(cif_get(row, "label_atom_id"))
    else:
        chain = clean_cif_value(cif_get(row, "auth_asym_id"), clean_cif_value(cif_get(row, "label_asym_id")))
        resseq = parse_cif_int(cif_get(row, "auth_seq_id"), default=parse_cif_int(cif_get(row, "label_seq_id"), default=0))
        resname = clean_cif_value(cif_get(row, "auth_comp_id"), clean_cif_value(cif_get(row, "label_comp_id"))).upper()
        atom_name = clean_cif_value(cif_get(row, "auth_atom_id"), clean_cif_value(cif_get(row, "label_atom_id")))

    element = clean_cif_value(cif_get(row, "type_symbol")).upper()
    if not element:
        element = infer_element(atom_name)

    icode = clean_cif_value(cif_get(row, "pdbx_PDB_ins_code"))
    altloc = clean_cif_value(cif_get(row, "label_alt_id"))

    return Atom(
        record=group,
        serial=parse_cif_int(cif_get(row, "id"), default=0),
        atom_name=atom_name,
        altloc=altloc,
        resname=resname,
        chain=chain,
        resseq=resseq,
        icode=icode,
        x=parse_cif_float(cif_get(row, "Cartn_x")),
        y=parse_cif_float(cif_get(row, "Cartn_y")),
        z=parse_cif_float(cif_get(row, "Cartn_z")),
        occupancy=parse_cif_float(cif_get(row, "occupancy"), default=0.0),
        element=element,
    )


def read_cif_first_model(path: Path, cif_id_mode: str = "auth") -> List[Atom]:
    raw_atoms: List[Atom] = []
    for row in iter_cif_atom_site_rows(path):
        atom = atom_from_cif_row(row, cif_id_mode=cif_id_mode)
        if atom is not None:
            raw_atoms.append(atom)
    return collapse_altlocs(raw_atoms)


def is_cif_path(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz"))


def read_structure_first_model(path: Path, cif_id_mode: str = "auth") -> List[Atom]:
    if is_cif_path(path):
        return read_cif_first_model(path, cif_id_mode=cif_id_mode)
    return read_pdb_first_model(path)


def dist_xyz(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(
        (a[0] - b[0]) ** 2 +
        (a[1] - b[1]) ** 2 +
        (a[2] - b[2]) ** 2
    )


def centroid(atoms: Iterable[Atom]) -> Tuple[float, float, float]:
    atom_list = list(atoms)
    if not atom_list:
        raise ValueError("Cannot calculate centroid from zero atoms.")
    return (
        sum(a.x for a in atom_list) / len(atom_list),
        sum(a.y for a in atom_list) / len(atom_list),
        sum(a.z for a in atom_list) / len(atom_list),
    )


def group_residues(atoms: Iterable[Atom]) -> Dict[Tuple[str, int, str, str], List[Atom]]:
    grouped: Dict[Tuple[str, int, str, str], List[Atom]] = {}
    for atom in atoms:
        grouped.setdefault(atom.residue_key, []).append(atom)
    return grouped


def find_reference_residue(
    atoms: List[Atom],
    ref_chain: Optional[str],
    ref_resseq: int,
    ref_resname: Optional[str],
    ref_icode: str,
) -> Tuple[Tuple[str, int, str, str], List[Atom]]:
    residues = group_residues(atoms)
    ref_resname_upper = ref_resname.upper() if ref_resname else None

    matches = [
        (key, residue_atoms)
        for key, residue_atoms in residues.items()
        if key[1] == ref_resseq
        and key[2] == ref_icode
        and (ref_chain is None or key[0] == ref_chain)
        and (ref_resname_upper is None or key[3] == ref_resname_upper)
    ]

    if not matches:
        name_text = f" {ref_resname_upper}" if ref_resname_upper else ""
        chain_text = f" chain {ref_chain}" if ref_chain else ""
        raise ValueError(f"Reference residue{name_text} {ref_resseq}{ref_icode}{chain_text} was not found.")

    if len(matches) > 1:
        labels = ", ".join(f"{k[3]}:{k[0] or '<blank>'}:{k[1]}{k[2]}" for k, _ in matches)
        raise ValueError(
            f"Found multiple matching reference residues: {labels}. "
            "Use --ref-chain and/or --ref-resname to choose one."
        )

    return matches[0]


def parse_name_list(text: str) -> List[str]:
    values = [x.strip() for x in text.split(",") if x.strip()]
    if not values:
        raise ValueError("The comma-separated atom list is empty.")
    return values


def parse_optional_name_list(text: Optional[str]) -> Optional[List[str]]:
    if text is None or text.strip() == "":
        return None
    return parse_name_list(text)


def reference_centroid_atoms(residue_atoms: List[Atom], mode: str, ref_atoms: Optional[List[str]]) -> List[Atom]:
    if mode == "heavy":
        selected = [a for a in residue_atoms if a.element != "H"]
    elif mode == "all":
        selected = list(residue_atoms)
    elif mode == "sidechain-heavy":
        selected = [a for a in residue_atoms if a.element != "H" and a.atom_name not in BACKBONE_ATOMS]
    elif mode == "backbone-heavy":
        selected = [a for a in residue_atoms if a.element != "H" and a.atom_name in BACKBONE_ATOMS]
    elif mode == "named":
        if not ref_atoms:
            raise ValueError("--ref-centroid-atoms named requires --ref-atoms, e.g. --ref-atoms CG,CD1,NE1.")
        wanted = set(ref_atoms)
        selected = [a for a in residue_atoms if a.atom_name in wanted]
        found = {a.atom_name for a in selected}
        missing = [name for name in ref_atoms if name not in found]
        if missing:
            raise ValueError("Reference residue is missing requested centroid atom(s): " + ",".join(missing))
    else:
        raise ValueError(f"Unknown --ref-centroid-atoms value: {mode}")

    if not selected:
        raise ValueError(f"No atoms selected for reference residue centroid using mode '{mode}'.")
    return selected


def is_likely_ligand_residue(res_atoms: List[Atom]) -> bool:
    if not res_atoms:
        return False

    resname = res_atoms[0].resname
    if resname in AA3:
        return False
    if resname in WATER_NAMES:
        return False
    if resname in COMMON_IONS:
        return False
    if resname in COMMON_SOLVENTS_BUFFERS:
        return False
    if not all(a.record == "HETATM" for a in res_atoms):
        return False

    heavy = [a for a in res_atoms if a.element != "H"]
    return len(heavy) >= 3


def select_ligands(
    atoms: List[Atom],
    ligand_resname: Optional[str],
    ligand_chain: Optional[str],
    ligand_resseq: Optional[int],
    ligand_icode: Optional[str],
    only_largest_ligand: bool,
) -> List[Tuple[Tuple[str, int, str, str], List[Atom]]]:
    residues = group_residues(atoms)
    matches: List[Tuple[Tuple[str, int, str, str], List[Atom]]] = []
    ligand_resname_upper = ligand_resname.upper() if ligand_resname else None

    for key, res_atoms in residues.items():
        chain, resseq, icode, resname = key

        if ligand_resname_upper is not None and resname != ligand_resname_upper:
            continue
        if ligand_chain is not None and chain != ligand_chain:
            continue
        if ligand_resseq is not None and resseq != ligand_resseq:
            continue
        if ligand_icode is not None and icode != ligand_icode:
            continue

        if ligand_resname_upper is None and not is_likely_ligand_residue(res_atoms):
            continue

        matches.append((key, res_atoms))

    matches.sort(key=lambda item: sum(1 for atom in item[1] if atom.element != "H"), reverse=True)

    if only_largest_ligand and matches:
        return matches[:1]
    return matches


def ligand_atom_map(ligand_atoms: List[Atom]) -> Dict[str, Atom]:
    mapping: Dict[str, Atom] = {}
    for atom in ligand_atoms:
        if atom.atom_name not in mapping or altloc_score(atom) > altloc_score(mapping[atom.atom_name]):
            mapping[atom.atom_name] = atom
    return mapping


def safe_column_name(atom_name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", atom_name.strip())
    return safe or "atom"


def fmt_distances(values: Sequence[float]) -> str:
    return ";".join(f"{value:.3f}" for value in values)


def fmt_names(values: Sequence[str]) -> str:
    return ";".join(values)


def parse_symmetric_groups(values: Sequence[str]) -> List[List[str]]:
    groups: List[List[str]] = []
    for value in values:
        names = parse_name_list(value)
        if len(names) < 2:
            raise ValueError("Each --symmetric-groups entry must contain at least two atom names.")
        groups.append(names)
    return groups


def normalize_element_symbol(text: str) -> str:
    value = text.strip()
    if not value:
        raise ValueError("Empty element symbol in ligand element-count specification.")
    # Preserve conventional capitalization: Cl, Br, Na, S, O, etc.
    if len(value) == 1:
        return value.upper()
    return value[0].upper() + value[1:].lower()


def parse_ligand_element_counts(text: Optional[str]) -> Optional[List[Tuple[str, int]]]:
    """
    Parse ligand atom selection by element. Accepted forms:
      S:1,O:2     explicit counts
      S,O,O       repeated elements, equivalent to S:1,O:2
      S,O         equivalent to S:1,O:1
    The order is retained for stable output columns.
    """
    if text is None or text.strip() == "":
        return None

    parts = [part.strip() for part in text.split(",") if part.strip()]
    if not parts:
        raise ValueError("--ligand-element-counts is empty.")

    ordered: List[Tuple[str, int]] = []
    index_by_element: Dict[str, int] = {}

    for part in parts:
        if ":" in part:
            element_text, count_text = part.split(":", 1)
            element = normalize_element_symbol(element_text)
            try:
                count = int(count_text)
            except ValueError:
                raise ValueError(f"Invalid count in --ligand-element-counts: {part}")
            if count < 1:
                raise ValueError(f"Element count must be at least 1 in --ligand-element-counts: {part}")
        else:
            element = normalize_element_symbol(part)
            count = 1

        if element in index_by_element:
            idx = index_by_element[element]
            old_element, old_count = ordered[idx]
            ordered[idx] = (old_element, old_count + count)
        else:
            index_by_element[element] = len(ordered)
            ordered.append((element, count))

    return ordered


def element_sort_key(atom: Atom, ref_centroid: Sequence[float], mode: str):
    if mode == "distance":
        return (dist_xyz(atom.coord, ref_centroid), atom.serial, atom.atom_name)
    if mode == "serial":
        return (atom.serial, atom.atom_name)
    if mode == "name":
        return (atom.atom_name, atom.serial)
    raise ValueError(f"Unknown --ligand-element-sort value: {mode}")


def select_ligand_atoms_for_measurement(
    ligand_atoms: List[Atom],
    ref_centroid: Sequence[float],
    ligand_atom_names: Optional[List[str]],
    element_counts: Optional[List[Tuple[str, int]]],
    element_sort: str,
) -> Tuple[str, str, List[SelectedLigandAtom], List[str]]:
    """
    Select ligand atoms for measurement. Exact atom-name selection is retained
    for backward compatibility. Element-count selection is preferred when atom
    names vary among PDB/mmCIF files, e.g. --ligand-element-counts S:1,O:2.
    """
    if ligand_atom_names is not None and element_counts is not None:
        raise ValueError("Use either --ligand-atoms or --ligand-element-counts, not both.")

    if element_counts is not None:
        selected: List[SelectedLigandAtom] = []
        warnings: List[str] = []
        heavy_atoms = [a for a in ligand_atoms if a.element != "H"]

        for element, count in element_counts:
            candidates = [a for a in heavy_atoms if normalize_element_symbol(a.element) == element]
            candidates = sorted(candidates, key=lambda a: element_sort_key(a, ref_centroid, element_sort))

            if len(candidates) < count:
                warnings.append(f"requested {element}:{count} but found {len(candidates)}")
            elif len(candidates) > count:
                names = ";".join(a.atom_name for a in candidates[count:])
                warnings.append(f"{element}: selected {count} of {len(candidates)}; unselected={names}")

            for idx in range(1, count + 1):
                atom = candidates[idx - 1] if idx <= len(candidates) else None
                label = f"{element}{idx}"
                selected.append(SelectedLigandAtom(
                    label=label,
                    selector=f"{element}:{idx}",
                    atom=atom,
                    selection_mode="element-counts",
                ))

        requested = ",".join(f"{element}:{count}" for element, count in element_counts)
        return "element-counts", requested, selected, warnings

    if ligand_atom_names is not None:
        by_name = ligand_atom_map(ligand_atoms)
        selected = []
        for idx, atom_name in enumerate(ligand_atom_names, start=1):
            atom = by_name.get(atom_name)
            selected.append(SelectedLigandAtom(
                label=f"atom{idx}",
                selector=atom_name,
                atom=atom,
                selection_mode="atom-names",
            ))
        return "atom-names", fmt_names(ligand_atom_names), selected, []

    selected = []
    heavy_atoms = sorted([a for a in ligand_atoms if a.element != "H"], key=lambda a: (a.serial, a.atom_name))
    for idx, atom in enumerate(heavy_atoms, start=1):
        selected.append(SelectedLigandAtom(
            label=f"atom{idx}",
            selector=atom.atom_name,
            atom=atom,
            selection_mode="all-heavy",
        ))
    return "all-heavy", "all non-hydrogen ligand atoms", selected, []


def analyze_file(path: Path, args: argparse.Namespace) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    atoms = read_structure_first_model(path, cif_id_mode=args.cif_id_mode)
    ref_key, ref_atoms_all = find_reference_residue(
        atoms=atoms,
        ref_chain=args.ref_chain,
        ref_resseq=args.ref_resseq,
        ref_resname=args.ref_resname,
        ref_icode=args.ref_icode,
    )

    ref_atom_names = parse_optional_name_list(args.ref_atoms)
    centroid_atoms = reference_centroid_atoms(ref_atoms_all, args.ref_centroid_atoms, ref_atom_names)
    ref_centroid = centroid(centroid_atoms)

    ligand_atom_names = parse_name_list(args.ligand_atoms) if args.ligand_atoms else None
    ligand_element_counts = parse_ligand_element_counts(args.ligand_element_counts)
    symmetric_groups = parse_symmetric_groups(args.symmetric_groups)

    ligand_residues = select_ligands(
        atoms=atoms,
        ligand_resname=args.ligand_resname,
        ligand_chain=args.ligand_chain,
        ligand_resseq=args.ligand_resseq,
        ligand_icode=args.ligand_icode,
        only_largest_ligand=args.only_largest_ligand,
    )

    if not ligand_residues:
        ligand_desc = args.ligand_resname or "auto-detected ligand"
        raise ValueError(f"No matching ligand residue found for {ligand_desc}.")

    long_rows: List[Dict[str, str]] = []
    summary_rows: List[Dict[str, str]] = []

    ref_chain, ref_resseq, ref_icode, ref_resname = ref_key
    ref_id = f"{ref_resname}:{ref_chain}:{ref_resseq}{ref_icode}"
    ref_centroid_atom_names = fmt_names([a.atom_name for a in centroid_atoms])

    for lig_key, lig_atoms in ligand_residues:
        lig_chain, lig_resseq, lig_icode, lig_resname = lig_key

        selection_mode, requested_selection, selected_atoms, selection_warnings = select_ligand_atoms_for_measurement(
            ligand_atoms=lig_atoms,
            ref_centroid=ref_centroid,
            ligand_atom_names=ligand_atom_names,
            element_counts=ligand_element_counts,
            element_sort=args.ligand_element_sort,
        )

        row: Dict[str, str] = {
            "structure_file": path.name,
            "structure_path": str(path),
            "reference_residue_id": ref_id,
            "reference_resname": ref_resname,
            "reference_chain": ref_chain,
            "reference_resseq": str(ref_resseq),
            "reference_icode": ref_icode,
            "reference_centroid_mode": args.ref_centroid_atoms,
            "reference_centroid_atoms": ref_centroid_atom_names,
            "reference_centroid_x": f"{ref_centroid[0]:.3f}",
            "reference_centroid_y": f"{ref_centroid[1]:.3f}",
            "reference_centroid_z": f"{ref_centroid[2]:.3f}",
            "ligand_resname": lig_resname,
            "ligand_chain": lig_chain,
            "ligand_resseq": str(lig_resseq),
            "ligand_icode": lig_icode,
            "ligand_selection_mode": selection_mode,
            "requested_ligand_atoms": requested_selection,
            "selection_warnings": fmt_names(selection_warnings),
        }

        distances_by_name: Dict[str, float] = {}
        distances_by_label: Dict[str, float] = {}
        missing_labels: List[str] = []
        element_groups: Dict[str, List[Tuple[str, str, float]]] = {}

        for idx, selected in enumerate(selected_atoms, start=1):
            atom = selected.atom
            row[f"atom{idx}_selector"] = selected.selector
            row[f"atom{idx}_label"] = selected.label

            if atom is None:
                missing_labels.append(selected.label)
                row[f"atom{idx}_found"] = "no"
                row[f"atom{idx}_matched_name"] = ""
                row[f"atom{idx}_element"] = ""
                row[f"atom{idx}_distance_to_centroid_A"] = ""
                row[f"{safe_column_name(selected.label)}_atom_name"] = ""
                row[f"{safe_column_name(selected.label)}_distance_to_centroid_A"] = ""
                continue

            value = dist_xyz(atom.coord, ref_centroid)
            distances_by_name[atom.atom_name] = value
            distances_by_label[selected.label] = value

            row[f"atom{idx}_found"] = "yes"
            row[f"atom{idx}_matched_name"] = atom.atom_name
            row[f"atom{idx}_element"] = atom.element
            row[f"atom{idx}_distance_to_centroid_A"] = f"{value:.3f}"
            row[f"{safe_column_name(selected.label)}_atom_name"] = atom.atom_name
            row[f"{safe_column_name(selected.label)}_distance_to_centroid_A"] = f"{value:.3f}"

            normalized_element = normalize_element_symbol(atom.element)
            element_groups.setdefault(normalized_element, []).append((selected.label, atom.atom_name, value))

            long_rows.append({
                "structure_file": path.name,
                "structure_path": str(path),
                "reference_residue_id": ref_id,
                "reference_resname": ref_resname,
                "reference_chain": ref_chain,
                "reference_resseq": str(ref_resseq),
                "reference_icode": ref_icode,
                "reference_centroid_mode": args.ref_centroid_atoms,
                "reference_centroid_atoms": ref_centroid_atom_names,
                "reference_centroid_x": f"{ref_centroid[0]:.3f}",
                "reference_centroid_y": f"{ref_centroid[1]:.3f}",
                "reference_centroid_z": f"{ref_centroid[2]:.3f}",
                "ligand_resname": lig_resname,
                "ligand_chain": lig_chain,
                "ligand_resseq": str(lig_resseq),
                "ligand_icode": lig_icode,
                "ligand_selection_mode": selection_mode,
                "ligand_selector": selected.selector,
                "ligand_output_label": selected.label,
                "ligand_atom_name": atom.atom_name,
                "ligand_element": atom.element,
                "distance_to_centroid_A": f"{value:.3f}",
            })

        row["missing_ligand_atoms"] = fmt_names(missing_labels)
        row["n_missing_ligand_atoms"] = str(len(missing_labels))

        # Automatically provide sorted distance groups for element-count selections.
        # This is the robust replacement for atom-name based symmetric groups in CIFs.
        for element, entries in sorted(element_groups.items()):
            sorted_entries = sorted(entries, key=lambda item: item[2])
            prefix = f"element_{safe_column_name(element)}"
            row[f"{prefix}_labels_by_distance"] = fmt_names([label for label, _, _ in sorted_entries])
            row[f"{prefix}_atom_names_by_distance"] = fmt_names([atom_name for _, atom_name, _ in sorted_entries])
            row[f"{prefix}_distances_sorted_A"] = fmt_distances([value for _, _, value in sorted_entries])

        # Backward-compatible exact-name symmetric groups.
        for group_index, group in enumerate(symmetric_groups, start=1):
            available = [(name, distances_by_name[name]) for name in group if name in distances_by_name]
            missing = [name for name in group if name not in distances_by_name]
            available_sorted = sorted(available, key=lambda item: item[1])
            prefix = f"symgroup{group_index}"
            row[f"{prefix}_atoms"] = fmt_names(group)
            row[f"{prefix}_missing_atoms"] = fmt_names(missing)
            row[f"{prefix}_atom_names_by_distance"] = fmt_names([name for name, _ in available_sorted])
            row[f"{prefix}_distances_sorted_A"] = fmt_distances([value for _, value in available_sorted])

        summary_rows.append(row)

    return long_rows, summary_rows


def read_path_list(path: Path) -> List[str]:
    items: List[str] = []
    with path.open("r", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            items.append(text)
    return items


def looks_like_glob(text: str) -> bool:
    return any(char in text for char in "*?[")


def split_patterns(pattern: str) -> List[str]:
    return [part.strip() for part in pattern.split(",") if part.strip()]


def expand_directory(directory: Path, pattern: str, recursive: bool) -> List[Path]:
    paths: List[Path] = []
    for one_pattern in split_patterns(pattern):
        if recursive:
            paths.extend(directory.rglob(one_pattern))
        else:
            paths.extend(directory.glob(one_pattern))
    return sorted(paths)


def expand_paths(
    specs: Sequence[str],
    input_dirs: Sequence[str],
    input_globs: Sequence[str],
    pattern: str,
    recursive: bool,
    input_list: Optional[str],
) -> List[Path]:
    all_specs: List[str] = []
    all_specs.extend(specs)
    all_specs.extend(input_dirs)
    all_specs.extend(input_globs)

    if input_list is not None:
        all_specs.extend(read_path_list(Path(input_list)))

    paths: List[Path] = []
    for spec in all_specs:
        candidate = Path(spec)
        if candidate.is_dir():
            paths.extend(expand_directory(candidate, pattern, recursive))
            continue
        if looks_like_glob(spec):
            paths.extend(Path(x) for x in glob.glob(spec, recursive=recursive))
            continue
        if candidate.exists():
            paths.append(candidate)

    return sorted(set(path.resolve() for path in paths if path.is_file()))


def get_fieldnames(rows: List[Dict[str, str]], preferred_fields: Sequence[str]) -> List[str]:
    extra_fields = sorted({key for row in rows for key in row.keys()} - set(preferred_fields))
    return list(preferred_fields) + extra_fields


def write_csv(path: Path, rows: List[Dict[str, str]], preferred_fields: Sequence[str]) -> None:
    fieldnames = get_fieldnames(rows, preferred_fields)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Measure selected ligand-atom distances to the centroid of any selected receptor residue in PDB/mmCIF files.",
        epilog=(
            "Example: python measure_residue_centroid_ligand_atom_distances.py "
            "--input-glob '*.cif' --ref-chain A --ref-resname TRP --ref-resseq 223 "
            "--ligand-resname LIG1 --ligand-chain B --ligand-resseq 1 "
            "--ligand-element-counts S:1,O:2 --out-prefix distances"
        ),
    )

    parser.add_argument(
        "structures",
        nargs="*",
        help="Input PDB/mmCIF files, directories, or glob patterns. Directories are searched using --pattern.",
    )
    parser.add_argument(
        "--input-dir",
        dest="input_dirs",
        action="append",
        default=[],
        help="Directory containing PDB/mmCIF files. Can be used multiple times.",
    )
    parser.add_argument(
        "--input-glob",
        action="append",
        default=[],
        help="Glob pattern for input files, e.g. 'M1*/*.yob.pdb' or '*.cif'. Quote the glob.",
    )
    parser.add_argument(
        "--pattern",
        default="*.pdb,*.pdb.gz,*.cif,*.mmcif,*.cif.gz,*.mmcif.gz",
        help="Filename pattern used when an input is a directory. Comma-separated patterns are allowed.",
    )
    parser.add_argument("--recursive", action="store_true", help="Search input directories recursively.")
    parser.add_argument("--input-list", default=None, help="Text file with one input path, directory, or glob per line.")

    parser.add_argument("--version", action="version", version=f"measure_residue_centroid_ligand_atom_distances.py {SCRIPT_VERSION}")
    parser.add_argument("--out-prefix", default="residue_centroid_ligand_distances")

    parser.add_argument("--ref-chain", default=None, help="Chain containing the reference receptor residue, e.g. A.")
    parser.add_argument("--ref-resname", default=None, help="Reference residue name, e.g. TRP, TYR, ASP. Optional but recommended.")
    parser.add_argument(
        "--ref-resseq",
        "--ref-resnum",
        "--ref-residue-number",
        type=int,
        default=None,
        help="Residue number of the reference receptor residue, e.g. 223 or 378. Required.",
    )
    parser.add_argument("--ref-icode", default="", help="Insertion code of the reference residue, if any.")
    parser.add_argument(
        "--ref-centroid-atoms",
        choices=["heavy", "all", "sidechain-heavy", "backbone-heavy", "named"],
        default="heavy",
        help="Atoms used for the reference residue centroid. Default: heavy.",
    )
    parser.add_argument(
        "--ref-atoms",
        default=None,
        help="Comma-separated reference-residue atom names used when --ref-centroid-atoms named, e.g. CG,CD1,NE1.",
    )
    parser.add_argument(
        "--cif-id-mode",
        choices=["auth", "label"],
        default="auth",
        help="For mmCIF, use auth_* or label_* chain/residue/atom identifiers. Default: auth.",
    )

    parser.add_argument("--ligand-resname", default=None, help="Ligand residue name, e.g. UNL or LIG1. If omitted, ligands are auto-detected.")
    parser.add_argument("--ligand-chain", default=None, help="Ligand chain, e.g. L or B.")
    parser.add_argument("--ligand-resseq", type=int, default=None, help="Ligand residue number.")
    parser.add_argument("--ligand-icode", default=None, help="Ligand insertion code.")
    parser.add_argument("--only-largest-ligand", action="store_true", help="Analyze only the largest auto-detected ligand.")
    parser.add_argument(
        "--ligand-atoms",
        default=None,
        help=(
            "Comma-separated exact ligand atom names to measure, e.g. S,O1,O2 or S45,O43,O44. "
            "Use this only when atom names are stable. If omitted and --ligand-element-counts is also omitted, "
            "all non-hydrogen ligand atoms are measured."
        ),
    )
    parser.add_argument(
        "--ligand-element-counts",
        "--ligand-elements",
        default=None,
        help=(
            "Robust ligand atom selection by element counts, e.g. S:1,O:2. "
            "Also accepts repeated elements such as S,O,O. This is recommended for CIF files "
            "whose ligand atom names vary between structures."
        ),
    )
    parser.add_argument(
        "--ligand-element-sort",
        choices=["distance", "serial", "name"],
        default="distance",
        help=(
            "How to order multiple selected atoms of the same element. Default: distance, "
            "which makes symmetric atoms such as two oxygens comparable across structures."
        ),
    )
    parser.add_argument(
        "--symmetric-groups",
        action="append",
        default=[],
        help=(
            "Optional comma-separated group of interchangeable ligand atoms, e.g. --symmetric-groups O1,O2. "
            "Can be used multiple times. The summary CSV gets sorted distances for each group."
        ),
    )

    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.ref_resseq is None:
        print("Error: --ref-resseq is required.", file=sys.stderr)
        return 1

    structure_paths = expand_paths(
        specs=args.structures,
        input_dirs=args.input_dirs,
        input_globs=args.input_glob,
        pattern=args.pattern,
        recursive=args.recursive,
        input_list=args.input_list,
    )

    if not structure_paths:
        print(
            "No input PDB/mmCIF files found. Provide files/globs, a directory, --input-dir, --input-glob, or --input-list.",
            file=sys.stderr,
        )
        return 1

    all_long: List[Dict[str, str]] = []
    all_summary: List[Dict[str, str]] = []
    failures: List[Dict[str, str]] = []

    print(f"Processing {len(structure_paths)} structure file(s).")

    for path in structure_paths:
        try:
            long_rows, summary_rows = analyze_file(path, args)
            all_long.extend(long_rows)
            all_summary.extend(summary_rows)
        except Exception as exc:
            failures.append({"structure_path": str(path), "error": str(exc)})

    long_fields = [
        "structure_file", "structure_path",
        "reference_residue_id", "reference_resname", "reference_chain", "reference_resseq", "reference_icode",
        "reference_centroid_mode", "reference_centroid_atoms", "reference_centroid_x", "reference_centroid_y", "reference_centroid_z",
        "ligand_resname", "ligand_chain", "ligand_resseq", "ligand_icode",
        "ligand_selection_mode", "ligand_selector", "ligand_output_label",
        "ligand_atom_name", "ligand_element", "distance_to_centroid_A",
    ]

    summary_fields = [
        "structure_file", "structure_path",
        "reference_residue_id", "reference_resname", "reference_chain", "reference_resseq", "reference_icode",
        "reference_centroid_mode", "reference_centroid_atoms", "reference_centroid_x", "reference_centroid_y", "reference_centroid_z",
        "ligand_resname", "ligand_chain", "ligand_resseq", "ligand_icode",
        "ligand_selection_mode", "requested_ligand_atoms", "selection_warnings",
        "atom1_selector", "atom1_label", "atom1_found", "atom1_matched_name", "atom1_element", "atom1_distance_to_centroid_A",
        "atom2_selector", "atom2_label", "atom2_found", "atom2_matched_name", "atom2_element", "atom2_distance_to_centroid_A",
        "atom3_selector", "atom3_label", "atom3_found", "atom3_matched_name", "atom3_element", "atom3_distance_to_centroid_A",
        "missing_ligand_atoms", "n_missing_ligand_atoms",
    ]

    out_prefix = Path(args.out_prefix)
    long_csv = out_prefix.with_name(out_prefix.name + "_long.csv")
    summary_csv = out_prefix.with_name(out_prefix.name + "_summary.csv")
    failures_csv = out_prefix.with_name(out_prefix.name + "_failures.csv")

    write_csv(long_csv, all_long, long_fields)
    write_csv(summary_csv, all_summary, summary_fields)

    if failures:
        write_csv(failures_csv, failures, ["structure_path", "error"])

    print(f"Wrote {len(all_long)} atom-distance rows to {long_csv}")
    print(f"Wrote {len(all_summary)} ligand-summary rows to {summary_csv}")
    if failures:
        print(f"Warning: {len(failures)} file(s) failed. See {failures_csv}", file=sys.stderr)

    return 0 if not failures else 2


# ---------------------------------------------------------------------------
# wxPython graphical front end
# ---------------------------------------------------------------------------

def run_analysis_from_namespace(args, emit=None):
    """Run the same batch analysis using an argparse-like namespace."""
    if emit is None:
        emit = print

    if args.ref_resseq is None:
        raise ValueError("Reference residue number is required.")

    structure_paths = expand_paths(
        specs=args.structures,
        input_dirs=args.input_dirs,
        input_globs=args.input_glob,
        pattern=args.pattern,
        recursive=args.recursive,
        input_list=args.input_list,
    )

    if not structure_paths:
        raise ValueError(
            "No input PDB/mmCIF files found. Provide files, directories, glob patterns, or an input-list file."
        )

    all_long: List[Dict[str, str]] = []
    all_summary: List[Dict[str, str]] = []
    failures: List[Dict[str, str]] = []

    emit(f"Processing {len(structure_paths)} structure file(s).")

    for idx, path in enumerate(structure_paths, start=1):
        emit(f"[{idx}/{len(structure_paths)}] {path}")
        try:
            long_rows, summary_rows = analyze_file(path, args)
            all_long.extend(long_rows)
            all_summary.extend(summary_rows)
        except Exception as exc:
            failures.append({"structure_path": str(path), "error": str(exc)})
            emit(f"  FAILED: {exc}")

    long_fields = [
        "structure_file", "structure_path",
        "reference_residue_id", "reference_resname", "reference_chain", "reference_resseq", "reference_icode",
        "reference_centroid_mode", "reference_centroid_atoms", "reference_centroid_x", "reference_centroid_y", "reference_centroid_z",
        "ligand_resname", "ligand_chain", "ligand_resseq", "ligand_icode",
        "ligand_selection_mode", "ligand_selector", "ligand_output_label",
        "ligand_atom_name", "ligand_element", "distance_to_centroid_A",
    ]

    summary_fields = [
        "structure_file", "structure_path",
        "reference_residue_id", "reference_resname", "reference_chain", "reference_resseq", "reference_icode",
        "reference_centroid_mode", "reference_centroid_atoms", "reference_centroid_x", "reference_centroid_y", "reference_centroid_z",
        "ligand_resname", "ligand_chain", "ligand_resseq", "ligand_icode",
        "ligand_selection_mode", "requested_ligand_atoms", "selection_warnings",
        "atom1_selector", "atom1_label", "atom1_found", "atom1_matched_name", "atom1_element", "atom1_distance_to_centroid_A",
        "atom2_selector", "atom2_label", "atom2_found", "atom2_matched_name", "atom2_element", "atom2_distance_to_centroid_A",
        "atom3_selector", "atom3_label", "atom3_found", "atom3_matched_name", "atom3_element", "atom3_distance_to_centroid_A",
        "missing_ligand_atoms", "n_missing_ligand_atoms",
    ]

    out_prefix = Path(args.out_prefix).expanduser()
    if out_prefix.suffix.lower() == ".csv":
        out_prefix = out_prefix.with_suffix("")
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    long_csv = out_prefix.with_name(out_prefix.name + "_long.csv")
    summary_csv = out_prefix.with_name(out_prefix.name + "_summary.csv")
    failures_csv = out_prefix.with_name(out_prefix.name + "_failures.csv")

    write_csv(long_csv, all_long, long_fields)
    write_csv(summary_csv, all_summary, summary_fields)

    if failures:
        write_csv(failures_csv, failures, ["structure_path", "error"])

    emit(f"Wrote {len(all_long)} atom-distance rows to {long_csv}")
    emit(f"Wrote {len(all_summary)} ligand-summary rows to {summary_csv}")
    if failures:
        emit(f"Warning: {len(failures)} file(s) failed. See {failures_csv}")

    return {
        "long_csv": str(long_csv),
        "summary_csv": str(summary_csv),
        "failures_csv": str(failures_csv) if failures else "",
        "n_structures": len(structure_paths),
        "n_long_rows": len(all_long),
        "n_summary_rows": len(all_summary),
        "n_failures": len(failures),
    }


def run_wx_app() -> int:
    try:
        import threading
        import traceback
        from types import SimpleNamespace
        import wx
    except ImportError as exc:
        print("This graphical version requires wxPython.", file=sys.stderr)
        print("On Debian/Ubuntu, try: sudo apt install python3-wxgtk4.0", file=sys.stderr)
        print("Or install with pip in a compatible environment: pip install wxPython", file=sys.stderr)
        print(f"Import error: {exc}", file=sys.stderr)
        return 1

    def blank_to_none(value: str) -> Optional[str]:
        value = value.strip()
        return value if value else None

    def parse_optional_int(value: str, field_name: str) -> Optional[int]:
        value = value.strip()
        if not value:
            return None
        try:
            return int(value)
        except ValueError as exc:
            raise ValueError(f"{field_name} must be an integer.") from exc

    def split_multiline_specs(value: str) -> List[str]:
        return [line.strip() for line in value.splitlines() if line.strip()]

    def split_symmetric_groups(value: str) -> List[str]:
        groups: List[str] = []
        for chunk in re.split(r"[\n;]+", value):
            chunk = chunk.strip()
            if chunk:
                groups.append(chunk)
        return groups

    class DistanceFrame(wx.Frame):
        def __init__(self):
            super().__init__(None, title=f"Ligand atoms to residue centroid distances ({SCRIPT_VERSION})", size=(930, 820))
            self.worker = None
            self._build_ui()
            self.Centre()

        def _make_text(self, parent, value="", style=0, size=(-1, -1)):
            return wx.TextCtrl(parent, value=value, style=style, size=size)

        def _add_row(self, sizer, label, ctrl, help_text=None):
            lab = wx.StaticText(self.panel, label=label)
            if help_text:
                lab.SetToolTip(help_text)
                ctrl.SetToolTip(help_text)
            sizer.Add(lab, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
            sizer.Add(ctrl, 1, wx.EXPAND | wx.ALL, 4)

        def _section(self, parent, title):
            box = wx.StaticBox(parent, label=title)
            sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
            grid = wx.FlexGridSizer(cols=2, hgap=6, vgap=4)
            grid.AddGrowableCol(1, 1)
            sizer.Add(grid, 1, wx.EXPAND | wx.ALL, 6)
            return sizer, grid

        def _build_ui(self):
            root = wx.Panel(self)
            top = wx.BoxSizer(wx.VERTICAL)
            root.SetSizer(top)

            self.scroller = wx.ScrolledWindow(root, style=wx.VSCROLL)
            self.scroller.SetScrollRate(8, 8)
            top.Add(self.scroller, 1, wx.EXPAND)

            self.panel = self.scroller
            main = wx.BoxSizer(wx.VERTICAL)
            self.panel.SetSizer(main)

            # Input files
            input_box = wx.StaticBox(self.panel, label="Input structures")
            input_sizer = wx.StaticBoxSizer(input_box, wx.VERTICAL)
            self.input_specs = self._make_text(
                self.panel,
                value="",
                style=wx.TE_MULTILINE,
                size=(-1, 82),
            )
            self.input_specs.SetToolTip("One file, directory, or glob pattern per line. Examples: *.cif, M1*/*.yob.pdb")
            input_sizer.Add(wx.StaticText(self.panel, label="Files, directories, or glob patterns, one per line:"), 0, wx.LEFT | wx.RIGHT | wx.TOP, 6)
            input_sizer.Add(self.input_specs, 0, wx.EXPAND | wx.ALL, 6)

            row = wx.BoxSizer(wx.HORIZONTAL)
            self.btn_files = wx.Button(self.panel, label="Add files...")
            self.btn_dir = wx.Button(self.panel, label="Choose directory...")
            self.btn_list = wx.Button(self.panel, label="Input-list file...")
            row.Add(self.btn_files, 0, wx.RIGHT, 6)
            row.Add(self.btn_dir, 0, wx.RIGHT, 6)
            row.Add(self.btn_list, 0, wx.RIGHT, 6)
            input_sizer.Add(row, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)

            grid = wx.FlexGridSizer(cols=2, hgap=6, vgap=4)
            grid.AddGrowableCol(1, 1)
            self.input_dir = self._make_text(self.panel, value="")
            self.pattern = self._make_text(self.panel, value="*.pdb,*.pdb.gz,*.cif,*.mmcif,*.cif.gz,*.mmcif.gz")
            self.recursive = wx.CheckBox(self.panel, label="Search input directories recursively")
            self.input_list = self._make_text(self.panel, value="")
            self._add_row(grid, "Optional directory:", self.input_dir)
            self._add_row(grid, "Directory filename pattern:", self.pattern)
            self._add_row(grid, "Input-list file:", self.input_list)
            grid.Add(wx.StaticText(self.panel, label="Recursive:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
            grid.Add(self.recursive, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
            input_sizer.Add(grid, 0, wx.EXPAND | wx.ALL, 6)
            main.Add(input_sizer, 0, wx.EXPAND | wx.ALL, 8)

            # Reference residue
            ref_sizer, ref_grid = self._section(self.panel, "Reference receptor residue")
            self.ref_chain = self._make_text(self.panel, value="A")
            self.ref_resname = self._make_text(self.panel, value="TRP")
            self.ref_resseq = self._make_text(self.panel, value="")
            self.ref_icode = self._make_text(self.panel, value="")
            self.ref_centroid_atoms = wx.Choice(self.panel, choices=["heavy", "all", "sidechain-heavy", "backbone-heavy", "named"])
            self.ref_centroid_atoms.SetStringSelection("heavy")
            self.ref_atoms = self._make_text(self.panel, value="")
            self.cif_id_mode = wx.Choice(self.panel, choices=["auth", "label"])
            self.cif_id_mode.SetStringSelection("auth")
            self._add_row(ref_grid, "Reference chain:", self.ref_chain)
            self._add_row(ref_grid, "Reference residue name:", self.ref_resname)
            self._add_row(ref_grid, "Reference residue number:", self.ref_resseq, "Required. Example: 223, 378.")
            self._add_row(ref_grid, "Reference insertion code:", self.ref_icode)
            self._add_row(ref_grid, "Centroid atoms:", self.ref_centroid_atoms)
            self._add_row(ref_grid, "Named centroid atoms:", self.ref_atoms, "Used only when Centroid atoms = named. Example: CG,CD1,NE1")
            self._add_row(ref_grid, "CIF identifier mode:", self.cif_id_mode, "auth usually matches chain/residue numbers seen in molecular viewers.")
            main.Add(ref_sizer, 0, wx.EXPAND | wx.ALL, 8)

            # Ligand
            lig_sizer, lig_grid = self._section(self.panel, "Ligand residue and atom selection")
            self.ligand_resname = self._make_text(self.panel, value="")
            self.ligand_chain = self._make_text(self.panel, value="")
            self.ligand_resseq = self._make_text(self.panel, value="")
            self.ligand_icode = self._make_text(self.panel, value="")
            self.only_largest_ligand = wx.CheckBox(self.panel, label="Analyze only the largest auto-detected ligand")
            self.selection_mode = wx.RadioBox(
                self.panel,
                label="Ligand atom selection mode",
                choices=["By element counts", "By exact atom names", "All non-hydrogen atoms"],
                majorDimension=1,
                style=wx.RA_SPECIFY_COLS,
            )
            self.selection_mode.SetSelection(0)
            self.ligand_element_counts = self._make_text(self.panel, value="S:1,O:2")
            self.ligand_atoms = self._make_text(self.panel, value="")
            self.ligand_element_sort = wx.Choice(self.panel, choices=["distance", "serial", "name"])
            self.ligand_element_sort.SetStringSelection("distance")
            self.symmetric_groups = self._make_text(self.panel, value="", style=wx.TE_MULTILINE, size=(-1, 56))

            self._add_row(lig_grid, "Ligand residue name:", self.ligand_resname, "Optional. Example: UNL or LIG1. Leave blank for auto-detection.")
            self._add_row(lig_grid, "Ligand chain:", self.ligand_chain)
            self._add_row(lig_grid, "Ligand residue number:", self.ligand_resseq)
            self._add_row(lig_grid, "Ligand insertion code:", self.ligand_icode)
            lig_grid.Add(wx.StaticText(self.panel, label="Only largest ligand:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
            lig_grid.Add(self.only_largest_ligand, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
            lig_grid.Add(wx.StaticText(self.panel, label="Selection mode:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
            lig_grid.Add(self.selection_mode, 1, wx.EXPAND | wx.ALL, 4)
            self._add_row(lig_grid, "Element counts:", self.ligand_element_counts, "Recommended for CIF. Example: S:1,O:2 or S,O,O")
            self._add_row(lig_grid, "Exact atom names:", self.ligand_atoms, "Example: S,O,O1 or S45,O43,O44")
            self._add_row(lig_grid, "Sort same-element atoms by:", self.ligand_element_sort)
            self._add_row(lig_grid, "Symmetric atom-name groups:", self.symmetric_groups, "Optional, one group per line, e.g. O1,O2. Mostly useful only when atom names are stable.")
            main.Add(lig_sizer, 0, wx.EXPAND | wx.ALL, 8)

            # Output
            out_sizer, out_grid = self._section(self.panel, "Output")
            self.out_prefix = self._make_text(self.panel, value="residue_centroid_ligand_distances")
            self.btn_out = wx.Button(self.panel, label="Choose output prefix...")
            out_row = wx.BoxSizer(wx.HORIZONTAL)
            out_row.Add(self.out_prefix, 1, wx.EXPAND | wx.RIGHT, 6)
            out_row.Add(self.btn_out, 0)
            out_grid.Add(wx.StaticText(self.panel, label="Output prefix:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 4)
            out_grid.Add(out_row, 1, wx.EXPAND | wx.ALL, 4)
            main.Add(out_sizer, 0, wx.EXPAND | wx.ALL, 8)

            # Buttons and log
            button_row = wx.BoxSizer(wx.HORIZONTAL)
            self.run_button = wx.Button(self.panel, label="Run analysis")
            self.close_button = wx.Button(self.panel, label="Close")
            button_row.AddStretchSpacer(1)
            button_row.Add(self.run_button, 0, wx.RIGHT, 8)
            button_row.Add(self.close_button, 0)
            main.Add(button_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

            self.log = wx.TextCtrl(self.panel, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.HSCROLL, size=(-1, 160))
            main.Add(wx.StaticText(self.panel, label="Log:"), 0, wx.LEFT | wx.RIGHT, 8)
            main.Add(self.log, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)

            self.btn_files.Bind(wx.EVT_BUTTON, self.on_add_files)
            self.btn_dir.Bind(wx.EVT_BUTTON, self.on_choose_dir)
            self.btn_list.Bind(wx.EVT_BUTTON, self.on_choose_list)
            self.btn_out.Bind(wx.EVT_BUTTON, self.on_choose_output)
            self.run_button.Bind(wx.EVT_BUTTON, self.on_run)
            self.close_button.Bind(wx.EVT_BUTTON, lambda evt: self.Close())

        def append_specs(self, paths: Sequence[str]):
            existing = self.input_specs.GetValue().strip()
            addition = "\n".join(paths)
            self.input_specs.SetValue((existing + "\n" + addition).strip() if existing else addition)

        def on_add_files(self, event):
            wildcard = "Structure files (*.pdb;*.cif;*.mmcif;*.pdb.gz;*.cif.gz)|*.pdb;*.cif;*.mmcif;*.pdb.gz;*.cif.gz|All files (*.*)|*.*"
            with wx.FileDialog(self, "Choose structure files", wildcard=wildcard, style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_FILE_MUST_EXIST) as dlg:
                if dlg.ShowModal() == wx.ID_OK:
                    self.append_specs(dlg.GetPaths())

        def on_choose_dir(self, event):
            with wx.DirDialog(self, "Choose input directory", style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST) as dlg:
                if dlg.ShowModal() == wx.ID_OK:
                    self.input_dir.SetValue(dlg.GetPath())

        def on_choose_list(self, event):
            with wx.FileDialog(self, "Choose input-list text file", wildcard="Text files (*.txt)|*.txt|All files (*.*)|*.*", style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as dlg:
                if dlg.ShowModal() == wx.ID_OK:
                    self.input_list.SetValue(dlg.GetPath())

        def on_choose_output(self, event):
            with wx.FileDialog(self, "Choose output prefix", wildcard="CSV prefix (*.csv)|*.csv|All files (*.*)|*.*", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
                if dlg.ShowModal() == wx.ID_OK:
                    path = Path(dlg.GetPath())
                    if path.suffix.lower() == ".csv":
                        path = path.with_suffix("")
                    self.out_prefix.SetValue(str(path))

        def collect_args(self):
            structures = split_multiline_specs(self.input_specs.GetValue())
            input_dirs = []
            if self.input_dir.GetValue().strip():
                input_dirs.append(self.input_dir.GetValue().strip())

            mode = self.selection_mode.GetSelection()
            ligand_element_counts = None
            ligand_atoms = None
            if mode == 0:
                ligand_element_counts = self.ligand_element_counts.GetValue().strip()
                if not ligand_element_counts:
                    raise ValueError("Element-count selection requires a value such as S:1,O:2.")
            elif mode == 1:
                ligand_atoms = self.ligand_atoms.GetValue().strip()
                if not ligand_atoms:
                    raise ValueError("Exact atom-name selection requires atom names such as S,O1,O2.")

            out_prefix = self.out_prefix.GetValue().strip()
            if not out_prefix:
                raise ValueError("Output prefix is required.")

            ref_resseq = parse_optional_int(self.ref_resseq.GetValue(), "Reference residue number")
            if ref_resseq is None:
                raise ValueError("Reference residue number is required.")

            return SimpleNamespace(
                structures=structures,
                input_dirs=input_dirs,
                input_glob=[],
                pattern=self.pattern.GetValue().strip() or "*.pdb,*.pdb.gz,*.cif,*.mmcif,*.cif.gz,*.mmcif.gz",
                recursive=self.recursive.GetValue(),
                input_list=blank_to_none(self.input_list.GetValue()),
                out_prefix=out_prefix,
                ref_chain=blank_to_none(self.ref_chain.GetValue()),
                ref_resname=blank_to_none(self.ref_resname.GetValue()),
                ref_resseq=ref_resseq,
                ref_icode=self.ref_icode.GetValue().strip(),
                ref_centroid_atoms=self.ref_centroid_atoms.GetStringSelection(),
                ref_atoms=blank_to_none(self.ref_atoms.GetValue()),
                cif_id_mode=self.cif_id_mode.GetStringSelection(),
                ligand_resname=blank_to_none(self.ligand_resname.GetValue()),
                ligand_chain=blank_to_none(self.ligand_chain.GetValue()),
                ligand_resseq=parse_optional_int(self.ligand_resseq.GetValue(), "Ligand residue number"),
                ligand_icode=blank_to_none(self.ligand_icode.GetValue()),
                only_largest_ligand=self.only_largest_ligand.GetValue(),
                ligand_atoms=ligand_atoms,
                ligand_element_counts=ligand_element_counts,
                ligand_element_sort=self.ligand_element_sort.GetStringSelection(),
                symmetric_groups=split_symmetric_groups(self.symmetric_groups.GetValue()),
            )

        def log_line(self, text: str):
            self.log.AppendText(text + "\n")

        def on_run(self, event):
            if self.worker is not None and self.worker.is_alive():
                wx.MessageBox("Analysis is already running.", "Busy", wx.OK | wx.ICON_INFORMATION)
                return
            try:
                args = self.collect_args()
            except Exception as exc:
                wx.MessageBox(str(exc), "Invalid input", wx.OK | wx.ICON_ERROR)
                return

            self.log.SetValue("")
            self.run_button.Disable()
            self.log_line("Starting analysis...")

            def emit(message: str):
                wx.CallAfter(self.log_line, message)

            def worker():
                try:
                    result = run_analysis_from_namespace(args, emit=emit)
                    wx.CallAfter(self.on_done, result, None)
                except Exception:
                    wx.CallAfter(self.on_done, None, traceback.format_exc())

            self.worker = threading.Thread(target=worker, daemon=True)
            self.worker.start()

        def on_done(self, result, error):
            self.run_button.Enable()
            if error:
                self.log_line("\nERROR:\n" + error)
                wx.MessageBox("Analysis failed. See the log for details.", "Error", wx.OK | wx.ICON_ERROR)
                return

            self.log_line("\nDone.")
            msg = (
                f"Processed {result['n_structures']} structure file(s).\n"
                f"Summary rows: {result['n_summary_rows']}\n"
                f"Atom-distance rows: {result['n_long_rows']}\n\n"
                f"Summary CSV:\n{result['summary_csv']}\n\n"
                f"Long CSV:\n{result['long_csv']}"
            )
            if result.get("failures_csv"):
                msg += f"\n\nFailures CSV:\n{result['failures_csv']}"
            wx.MessageBox(msg, "Analysis complete", wx.OK | wx.ICON_INFORMATION)

    app = wx.App(False)
    frame = DistanceFrame()
    frame.Show()
    app.MainLoop()
    return 0


if __name__ == "__main__":
    raise SystemExit(run_wx_app())
