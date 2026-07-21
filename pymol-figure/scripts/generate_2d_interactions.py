"""
Generate a LigPlot-inspired 2D interaction diagram using RDKit and Pillow.

Primary docking workflow:

    python generate_2d_interactions.py \
        --protein receptor.pdb \
        --ligand docked.sdf \
        --pose-index 1 \
        --output output_dir

SDF and MOL2 are preferred because they preserve ligand bond orders and charges.
"""

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import rdDepictor

import auto_detect_interactions as detection


STYLE = {
    "hbond": ("Hydrogen bond", "#169b55"),
    "pi-pi": ("Pi stacking", "#d49a00"),
    "salt-bridge": ("Salt bridge", "#a23db5"),
    "hydrophobic": ("Hydrophobic", "#d7191c"),
    "pocket": ("Pocket context", "#555555"),
}

PRIORITY = {
    "salt-bridge": 0,
    "hbond": 1,
    "pi-pi": 2,
    "hydrophobic": 3,
    "pocket": 4,
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate a LigPlot-inspired 2D docking interaction diagram")
    parser.add_argument("--protein", required=True, help="Receptor PDB file")
    parser.add_argument(
        "--ligand", required=True,
        help="Docked ligand poses in SDF or MOL2 format")
    parser.add_argument(
        "--pose-index", type=int, default=1,
        help="One-based docking pose index (default: 1)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--max-residues", type=int, default=0,
        help="Optional limit for pocket-context-only labels; 0 shows the complete pocket. Direct interactions are never omitted.")
    parser.add_argument(
        "--pocket-cutoff", type=float, default=5.0,
        help="Show standard amino acids within this distance of the ligand "
             "(default: 5.0 A)")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1100)
    parser.add_argument(
        "--opaque-background", action="store_true",
        help="Use a white background instead of the default transparent PNG")
    parser.add_argument(
        "--draw-hydrophobic-lines", action="store_true",
        help="Draw optional red dotted hydrophobic contact lines and distances")
    return parser.parse_args(argv)


def load_protein(path):
    path = Path(path)
    if path.suffix.lower() not in {".pdb", ".ent"}:
        raise ValueError("The receptor must be a PDB file")
    mol = Chem.MolFromPDBFile(
        str(path), removeHs=False, sanitize=True, proximityBonding=True)
    if mol is None:
        raise ValueError(f"Could not parse receptor: {path}")
    return mol


def load_ligand_poses(path):
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".sdf", ".sd"}:
        supplier = Chem.SDMolSupplier(
            str(path), removeHs=False, sanitize=True)
        poses = [mol for mol in supplier if mol is not None]
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(
            str(path), removeHs=False, sanitize=True)
        poses = [mol] if mol is not None else []
    else:
        raise ValueError(
            "Ligand poses must be SDF or MOL2. Convert PDBQT with Meeko so "
            "bond orders and formal charges are preserved.")
    if not poses:
        raise ValueError(f"No valid ligand poses found in: {path}")
    return poses


def assign_ligand_residue_info(mol):
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol()}{atom.GetIdx() + 1}"[:4].rjust(4))
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId("L")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    return mol


def combine_complex(protein, ligand):
    ligand = assign_ligand_residue_info(ligand)
    combined = Chem.CombineMols(protein, ligand)
    Chem.GetSymmSSSR(combined)
    protein_atoms = list(range(protein.GetNumAtoms()))
    ligand_atoms = list(range(protein.GetNumAtoms(), combined.GetNumAtoms()))
    return combined, protein_atoms, ligand_atoms


def _residue(atom):
    info = atom.GetPDBResidueInfo()
    return (
        info.GetResidueName().strip(),
        int(info.GetResidueNumber()),
        info.GetChainId().strip(),
    )


def _record(residue, interaction_type, distance, ligand_indices):
    name, number, chain = residue
    return {
        "name": name,
        "number": number,
        "chain": chain,
        "label": f"{name}{number}" + (f"({chain})" if chain else ""),
        "type": interaction_type,
        "distance": round(float(distance), 3),
        "ligand_atoms": sorted(set(int(i) for i in ligand_indices)),
    }


def detect_hbonds(mol, protein_atoms, ligand_atoms, ligand_offset, cutoff=3.5):
    best = {}
    for protein_index in protein_atoms:
        protein_atom = mol.GetAtomWithIdx(protein_index)
        if not detection._is_polar_heavy(protein_atom):
            continue
        info = protein_atom.GetPDBResidueInfo()
        if not info:
            continue
        residue = _residue(protein_atom)
        atom_name = info.GetName().strip()
        donors, acceptors = detection._protein_donor_acceptor_names(residue[0])
        protein_is_donor = atom_name in donors
        protein_is_acceptor = atom_name in acceptors
        if not (protein_is_donor or protein_is_acceptor):
            continue
        for ligand_index in ligand_atoms:
            ligand_atom = mol.GetAtomWithIdx(ligand_index)
            if not detection._is_polar_heavy(ligand_atom):
                continue
            compatible = (
                protein_is_donor and detection._ligand_is_acceptor(ligand_atom)
            ) or (
                protein_is_acceptor
                and detection._ligand_is_donor(ligand_atom, False)
            )
            if not compatible:
                continue
            dist = detection.distance(
                detection.atom_pos(mol, protein_index),
                detection.atom_pos(mol, ligand_index),
            )
            effective_cutoff = min(cutoff, 3.2) if atom_name in {"N", "O"} else cutoff
            if dist <= effective_cutoff and (
                    residue not in best or dist < best[residue][0]):
                best[residue] = (dist, ligand_index - ligand_offset)
    return [
        _record(residue, "hbond", dist, [ligand_index])
        for residue, (dist, ligand_index) in best.items()
    ]


def detect_hydrophobic(mol, protein_atoms, ligand_atoms, ligand_offset,
                       cutoff=4.0):
    groups = detection.get_res_groups(mol, protein_atoms)
    ligand_contact_atoms = [
        index for index in ligand_atoms
        if mol.GetAtomWithIdx(index).GetAtomicNum() in (6, 16, 17, 35, 53)
    ]
    records = []
    for residue, atom_indices in groups.items():
        atom_map = detection.get_atom_map(mol, atom_indices)
        candidates = detection.HYDROPHOBIC_ATOMS.get(residue[0], [])
        best = None
        for atom_name in candidates:
            protein_index = atom_map.get(atom_name)
            if protein_index is None:
                continue
            for ligand_index in ligand_contact_atoms:
                dist = detection.distance(
                    detection.atom_pos(mol, protein_index),
                    detection.atom_pos(mol, ligand_index),
                )
                if dist <= cutoff and (best is None or dist < best[0]):
                    best = (dist, ligand_index - ligand_offset)
        if best:
            records.append(_record(
                residue, "hydrophobic", best[0], [best[1]]))
    return records


def detect_pocket_context(mol, protein_atoms, ligand_atoms, ligand_offset,
                          cutoff=5.0):
    groups = detection.get_res_groups(mol, protein_atoms)
    ligand_heavy = [
        index for index in ligand_atoms
        if mol.GetAtomWithIdx(index).GetAtomicNum() > 1
    ]
    records = []
    for residue, atom_indices in groups.items():
        if residue[0] not in detection.PROTEIN_LIKE_AA:
            continue
        best = None
        for protein_index in atom_indices:
            if mol.GetAtomWithIdx(protein_index).GetAtomicNum() <= 1:
                continue
            for ligand_index in ligand_heavy:
                dist = detection.distance(
                    detection.atom_pos(mol, protein_index),
                    detection.atom_pos(mol, ligand_index),
                )
                if dist <= cutoff and (best is None or dist < best[0]):
                    best = (dist, ligand_index - ligand_offset)
        if best:
            records.append(_record(
                residue, "pocket", best[0], [best[1]]))
    return records


def detect_pi_stacking(mol, protein_atoms, ligand_atoms, ligand_offset,
                       cutoff=5.5):
    groups = detection.get_res_groups(mol, protein_atoms)
    ligand_set = set(ligand_atoms)
    ligand_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6) or not set(ring).issubset(ligand_set):
            continue
        if all(mol.GetAtomWithIdx(index).GetAtomicNum() in (6, 7)
               for index in ring):
            positions = [detection.atom_pos(mol, index) for index in ring]
            ligand_rings.append((ring, detection.centroid(positions)))

    records = []
    for residue, atom_indices in groups.items():
        template = detection.PROTEIN_ARO_RINGS.get(residue[0])
        if not template:
            continue
        atom_map = detection.get_atom_map(mol, atom_indices)
        ring_indices = [
            atom_map[name] for name in template if name in atom_map]
        if len(ring_indices) < 4:
            continue
        protein_center = detection.centroid([
            detection.atom_pos(mol, index) for index in ring_indices])
        best = None
        for ring, ligand_center in ligand_rings:
            dist = detection.distance(protein_center, ligand_center)
            if dist <= cutoff and (best is None or dist < best[0]):
                best = (dist, ring)
        if best:
            records.append(_record(
                residue,
                "pi-pi",
                best[0],
                [index - ligand_offset for index in best[1]],
            ))
    return records


def detect_salt_bridges(mol, protein_atoms, ligand_atoms, ligand_offset,
                        cutoff=4.0):
    groups = detection.get_res_groups(mol, protein_atoms)
    records = []
    for residue, atom_indices in groups.items():
        atom_map = detection.get_atom_map(mol, atom_indices)
        protein_names = (
            detection.CATIONS.get(residue[0])
            or detection.ANIONS.get(residue[0])
            or []
        )
        protein_positive = residue[0] in detection.CATIONS
        best = None
        for atom_name in protein_names:
            protein_index = atom_map.get(atom_name)
            if protein_index is None:
                continue
            for ligand_index in ligand_atoms:
                ligand_atom = mol.GetAtomWithIdx(ligand_index)
                ligand_charge = ligand_atom.GetFormalCharge()
                compatible = (
                    protein_positive
                    and ligand_atom.GetAtomicNum() == 8
                    and ligand_charge <= 0
                ) or (
                    not protein_positive
                    and ligand_atom.GetAtomicNum() == 7
                    and ligand_charge >= 0
                )
                if not compatible:
                    continue
                dist = detection.distance(
                    detection.atom_pos(mol, protein_index),
                    detection.atom_pos(mol, ligand_index),
                )
                if dist <= cutoff and (best is None or dist < best[0]):
                    best = (dist, ligand_index - ligand_offset)
        if best:
            records.append(_record(
                residue, "salt-bridge", best[0], [best[1]]))
    return records


def detect_interactions(protein, ligand, pocket_cutoff=5.0):
    mol, protein_atoms, ligand_atoms = combine_complex(protein, ligand)
    offset = protein.GetNumAtoms()
    records = []
    records.extend(detect_hbonds(
        mol, protein_atoms, ligand_atoms, offset))
    records.extend(detect_pi_stacking(
        mol, protein_atoms, ligand_atoms, offset))
    records.extend(detect_salt_bridges(
        mol, protein_atoms, ligand_atoms, offset))
    records.extend(detect_hydrophobic(
        mol, protein_atoms, ligand_atoms, offset))
    interacting_residues = {
        (item["name"], item["number"], item["chain"]) for item in records
    }
    records.extend(
        item for item in detect_pocket_context(
            mol, protein_atoms, ligand_atoms, offset, pocket_cutoff)
        if (item["name"], item["number"], item["chain"])
        not in interacting_residues
    )
    records.sort(key=lambda item: (
        PRIORITY[item["type"]], item["distance"],
        item["chain"], item["number"]))
    return records


def _font(size, bold=False):
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else
        "C:/Windows/Fonts/arial.ttf",
        "Arial Bold.ttf" if bold else "Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


def _prepare_ligand(ligand):
    source = Chem.Mol(ligand)
    for atom in source.GetAtoms():
        atom.SetIntProp("_source_index", atom.GetIdx())
    drawing_mol = Chem.RemoveHs(source)
    Chem.AssignStereochemistry(drawing_mol, cleanIt=True, force=True)
    rdDepictor.Compute2DCoords(drawing_mol, clearConfs=True)
    Chem.WedgeMolBonds(drawing_mol, drawing_mol.GetConformer())
    stereo_report = verify_stereochemistry(ligand, drawing_mol)
    if not stereo_report["passed"]:
        raise ValueError(
            "2D stereochemistry verification failed: "
            + "; ".join(stereo_report["errors"]))
    mapping = {
        atom.GetIntProp("_source_index"): atom.GetIdx()
        for atom in drawing_mol.GetAtoms()
        if atom.HasProp("_source_index")
    }
    return drawing_mol, mapping


def _chiral_centers_by_index(mol):
    """Return RDKit's assigned/unassigned tetrahedral centres by atom index."""
    return dict(Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, includeCIP=True))


def verify_stereochemistry(source_mol, drawing_mol):
    """Check that assigned input chirality survives 2D preparation visibly.

    A stereocentre is considered renderable only when its R/S assignment is
    unchanged and an adjacent bond has a solid or hashed wedge direction.  The
    check deliberately does not require wedges for unassigned centres: drawing
    a wedge in that case would invent stereochemical information.
    """
    source_centers = _chiral_centers_by_index(source_mol)
    drawing_centers = _chiral_centers_by_index(drawing_mol)
    source_to_drawing = {
        atom.GetIntProp("_source_index"): atom.GetIdx()
        for atom in drawing_mol.GetAtoms()
        if atom.HasProp("_source_index")
    }
    wedge_directions = {Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH}
    errors = []
    assigned = {}
    unassigned = []

    for source_index, cip in source_centers.items():
        drawing_index = source_to_drawing.get(source_index)
        if drawing_index is None:
            errors.append(f"chiral atom {source_index} was removed during preparation")
            continue
        if cip in {"R", "S"}:
            assigned[source_index] = cip
            if drawing_centers.get(drawing_index) != cip:
                errors.append(
                    f"chiral atom {source_index} changed from {cip} to "
                    f"{drawing_centers.get(drawing_index, 'unassigned')}")
            if not any(
                    bond.GetBondDir() in wedge_directions
                    for bond in drawing_mol.GetAtomWithIdx(drawing_index).GetBonds()):
                errors.append(
                    f"chiral atom {source_index} ({cip}) has no wedge/hash bond")
        else:
            unassigned.append(source_index)
            if any(
                    bond.GetBondDir() in wedge_directions
                    for bond in drawing_mol.GetAtomWithIdx(drawing_index).GetBonds()):
                errors.append(
                    f"unassigned chiral atom {source_index} was given a wedge/hash bond")

    wedge_bonds = sum(
        bond.GetBondDir() in wedge_directions for bond in drawing_mol.GetBonds())
    return {
        "passed": not errors,
        "assigned_centers": assigned,
        "unassigned_centers": unassigned,
        "wedge_bond_count": wedge_bonds,
        "errors": errors,
    }


def _layout_ligand(mol, width, height, ligand_scale=1.0):
    coords = mol.GetConformer().GetPositions()
    min_x = min(point[0] for point in coords)
    max_x = max(point[0] for point in coords)
    min_y = min(point[1] for point in coords)
    max_y = max(point[1] for point in coords)
    source_width = max(max_x - min_x, 1.0)
    source_height = max(max_y - min_y, 1.0)
    scale = ligand_scale * min(
        width * 0.54 / source_width,
        height * 0.58 / source_height,
    )
    center_x, center_y = width / 2, height / 2 - 10
    return {
        atom.GetIdx(): (
            center_x
            + (coords[atom.GetIdx()][0] - (min_x + max_x) / 2) * scale,
            center_y
            - (coords[atom.GetIdx()][1] - (min_y + max_y) / 2) * scale,
        )
        for atom in mol.GetAtoms()
    }


def _draw_bond(draw, start, end, order):
    x1, y1 = start
    x2, y2 = end
    length = max(math.hypot(x2 - x1, y2 - y1), 1.0)
    nx, ny = -(y2 - y1) / length, (x2 - x1) / length
    offsets = [0] if order == 1 else ([-4, 4] if order == 2 else [-6, 0, 6])
    for offset in offsets:
        draw.line(
            (
                x1 + nx * offset, y1 + ny * offset,
                x2 + nx * offset, y2 + ny * offset,
            ),
            fill="#5526a3",
            width=6,
        )


def _draw_wedge_bond(draw, start, end, hashed=False):
    """Draw a stereobond with its narrow end at ``start``."""
    x1, y1 = start
    x2, y2 = end
    length = max(math.hypot(x2 - x1, y2 - y1), 1.0)
    nx, ny = -(y2 - y1) / length, (x2 - x1) / length
    half_width = 12
    color = "#5526a3"

    if not hashed:
        draw.polygon(
            (
                (x1, y1),
                (x2 + nx * half_width, y2 + ny * half_width),
                (x2 - nx * half_width, y2 - ny * half_width),
            ),
            fill=color,
        )
        return

    # A hashed wedge widens away from the stereocentre. The short transverse
    # strokes make the bond direction unambiguous without relying on opacity.
    for step in range(1, 9):
        fraction = step / 9
        x = x1 + (x2 - x1) * fraction
        y = y1 + (y2 - y1) * fraction
        width = max(1.5, half_width * fraction)
        draw.line(
            (
                x - nx * width, y - ny * width,
                x + nx * width, y + ny * width,
            ),
            fill=color,
            width=3,
        )


def _draw_ligand(draw, mol, points):
    for bond in mol.GetBonds():
        start = points[bond.GetBeginAtomIdx()]
        end = points[bond.GetEndAtomIdx()]
        direction = bond.GetBondDir()
        if direction == Chem.BondDir.BEGINWEDGE:
            _draw_wedge_bond(draw, start, end, hashed=False)
            continue
        if direction == Chem.BondDir.BEGINDASH:
            _draw_wedge_bond(draw, start, end, hashed=True)
            continue
        order = 2 if bond.GetIsAromatic() else max(
            1, min(3, int(round(bond.GetBondTypeAsDouble()))))
        _draw_bond(
            draw,
            start,
            end,
            order,
        )
    colors = {
        "C": ("#050505", "#ffffff"),
        "N": ("#1747d1", "#ffffff"),
        "O": ("#e32636", "#ffffff"),
        "S": ("#f4c430", "#111111"),
        "P": ("#e67e22", "#ffffff"),
        "F": ("#24a69a", "#ffffff"),
        "Cl": ("#2f9e44", "#ffffff"),
        "Br": ("#9c4a1a", "#ffffff"),
        "I": ("#6f42c1", "#ffffff"),
    }
    font = _font(16, bold=True)
    for atom in mol.GetAtoms():
        x, y = points[atom.GetIdx()]
        symbol = atom.GetSymbol()
        fill, text = colors.get(symbol, ("#777777", "#ffffff"))
        radius = 12 if symbol != "C" else 10
        draw.ellipse(
            (x - radius, y - radius, x + radius, y + radius),
            fill=fill, outline="#111111", width=1)
        if symbol != "C":
            draw.text((x, y), symbol, fill=text, font=font, anchor="mm")


def _dashed_line(draw, start, end, color, width=3, dash=10, gap=7):
    x1, y1 = start
    x2, y2 = end
    length = math.hypot(x2 - x1, y2 - y1)
    if length == 0:
        return
    ux, uy = (x2 - x1) / length, (y2 - y1) / length
    position = 0
    while position < length:
        segment_end = min(position + dash, length)
        draw.line(
            (
                x1 + ux * position, y1 + uy * position,
                x1 + ux * segment_end, y1 + uy * segment_end,
            ),
            fill=color, width=width)
        position += dash + gap


def _hydrophobic_fan(draw, center, angle_to_ligand):
    cx, cy = center
    radius = 28
    center_angle = math.degrees(angle_to_ligand)
    angles = (center_angle - 70, center_angle + 70)
    draw.arc(
        (cx - radius, cy - radius, cx + radius, cy + radius),
        start=angles[0], end=angles[1], fill=STYLE["hydrophobic"][1], width=4)
    for angle in range(
            int(round(angles[0])), int(round(angles[1])) + 1, 20):
        radians = math.radians(angle)
        inner = (
            cx + radius * math.cos(radians),
            cy + radius * math.sin(radians),
        )
        outer = (
            cx + (radius + 13) * math.cos(radians),
            cy + (radius + 13) * math.sin(radians),
        )
        draw.line((*inner, *outer), fill=STYLE["hydrophobic"][1], width=3)


def _boxes_overlap(first, second, padding=0):
    return not (
        first[2] + padding <= second[0]
        or second[2] + padding <= first[0]
        or first[3] + padding <= second[1]
        or second[3] + padding <= first[1]
    )


def _box_inside_canvas(box, width, height, margin=12):
    return (
        box[0] >= margin and box[1] >= margin
        and box[2] <= width - margin and box[3] <= height - margin
    )


def _verify_annotation_layout(
        residue_boxes, distance_boxes, ligand_box, fan_boxes, width, height):
    """Return a strict, machine-readable pre-export annotation QC report."""
    errors = []
    for index, box in enumerate(residue_boxes):
        if not _box_inside_canvas(box, width, height):
            errors.append(f"residue label {index} is outside the canvas")
        if _boxes_overlap(box, ligand_box, padding=8):
            errors.append(f"residue label {index} overlaps the ligand drawing")
        for other_index, other in enumerate(residue_boxes[:index]):
            if _boxes_overlap(box, other, padding=6):
                errors.append(
                    f"residue labels {other_index} and {index} overlap")
        for fan_index, fan_box in enumerate(fan_boxes):
            if _boxes_overlap(box, fan_box, padding=3):
                errors.append(
                    f"residue label {index} overlaps hydrophobic symbol {fan_index}")

    for index, box in enumerate(distance_boxes):
        if not _box_inside_canvas(box, width, height, margin=8):
            errors.append(f"distance label {index} is outside the canvas")
        if _boxes_overlap(box, ligand_box, padding=3):
            errors.append(f"distance label {index} overlaps the ligand drawing")
        for residue_index, residue_box in enumerate(residue_boxes):
            if _boxes_overlap(box, residue_box, padding=3):
                errors.append(
                    f"distance label {index} overlaps residue label {residue_index}")
        for other_index, other in enumerate(distance_boxes[:index]):
            if _boxes_overlap(box, other, padding=3):
                errors.append(
                    f"distance labels {other_index} and {index} overlap")
        for fan_index, fan_box in enumerate(fan_boxes):
            if _boxes_overlap(box, fan_box, padding=3):
                errors.append(
                    f"distance label {index} overlaps hydrophobic symbol {fan_index}")
    return {"passed": not errors, "errors": errors}


def _radial_label_layout(
        draw, entries, points, center, width, height, key_font, context_font):
    atom_boxes = [
        (x - 28, y - 28, x + 28, y + 28) for x, y in points.values()
    ]
    # Use the same complete ligand footprint as the export quality gate.
    ligand_footprint = (
        min(x for x, _ in points.values()) - 18,
        min(y for _, y in points.values()) - 18,
        max(x for x, _ in points.values()) + 18,
        max(y for _, y in points.values()) + 18,
    )
    occupied = []
    angle_offsets = [0] + [value for step in range(12, 181, 12)
                           for value in (step, -step)]
    canvas_scale = max(1.0, min(width / 1600.0, height / 1100.0))
    radii = [round(radius * canvas_scale)
             for radius in (145, 185, 235, 295, 365, 450, 550, 670, 810)]

    def desired_angle(entry):
        dx = entry["anchor"][0] - center[0]
        dy = entry["anchor"][1] - center[1]
        if math.hypot(dx, dy) < 35:
            key = entry["key"]
            seed = key[2] * 137.508 + sum(ord(char) for char in key[0] + key[1])
            return math.radians(seed % 360)
        return math.atan2(dy, dx)

    for entry in entries:
        entry["desired_angle"] = desired_angle(entry)

    for entry in sorted(
            entries,
            key=lambda item: (
                min(PRIORITY[value["type"]] for value in item["items"]),
                item["desired_angle"],
            )):
        is_key = any(
            item["type"] != "pocket" for item in entry["items"])
        font = key_font if is_key else context_font
        entry["is_key"] = is_key
        entry["font"] = font
        label = entry["items"][0]["label"]
        text_box = draw.textbbox((0, 0), label, font=font, anchor="mm")
        text_width = text_box[2] - text_box[0]
        text_height = text_box[3] - text_box[1]
        decoration_padding = (
            110 if any(item["type"] == "hydrophobic"
                      for item in entry["items"]) else 10
        )
        best = None
        for radius in radii:
            for offset in angle_offsets:
                angle = entry["desired_angle"] + math.radians(offset)
                x = entry["anchor"][0] + radius * math.cos(angle)
                y = entry["anchor"][1] + radius * math.sin(angle)
                box = (
                    x - text_width / 2 - decoration_padding,
                    y - text_height / 2 - decoration_padding,
                    x + text_width / 2 + decoration_padding,
                    y + text_height / 2 + decoration_padding,
                )
                boundary_penalty = (
                    max(35 - box[0], 0)
                    + max(35 - box[1], 0)
                    + max(box[2] - (width - 35), 0)
                    + max(box[3] - (height - 35), 0)
                ) * 1000
                overlap_penalty = sum(
                    100000 for other in occupied
                    if _boxes_overlap(box, other, padding=10)
                )
                ligand_penalty = (
                    250000 if _boxes_overlap(
                        box, ligand_footprint, padding=8) else 0
                ) + sum(
                    25000 for atom_box in atom_boxes
                    if _boxes_overlap(box, atom_box, padding=8)
                )
                invalid = (
                    boundary_penalty > 0
                    or any(_boxes_overlap(box, other, padding=10)
                           for other in occupied)
                    or _boxes_overlap(box, ligand_footprint, padding=8)
                    or any(_boxes_overlap(box, atom_box, padding=8)
                           for atom_box in atom_boxes)
                )
                score = (
                    boundary_penalty + overlap_penalty + ligand_penalty
                    + radius * 2 + abs(offset) * 12
                )
                if invalid:
                    score += 10000000
                if best is None or score < best[0]:
                    best = (score, x, y, box)
        entry["label_x"], entry["label_y"] = best[1], best[2]
        occupied.append(best[3])


def _line_end_at_label(anchor, label, text_width, text_height):
    dx = label[0] - anchor[0]
    dy = label[1] - anchor[1]
    length = max(math.hypot(dx, dy), 1.0)
    ux, uy = dx / length, dy / length
    limits = []
    if abs(ux) > 1e-6:
        limits.append((text_width / 2 + 16) / abs(ux))
    if abs(uy) > 1e-6:
        limits.append((text_height / 2 + 12) / abs(uy))
    offset = min(limits) if limits else 20
    return label[0] - ux * offset, label[1] - uy * offset


def _place_distance_label(
        draw, text, start, end, occupied, avoid_boxes, font, width, height):
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = max(math.hypot(dx, dy), 1.0)
    nx, ny = -dy / length, dx / length
    candidates = []
    for fraction in (0.52, 0.64, 0.40, 0.76, 0.28, 0.86):
        for normal_offset in (24, -24, 36, -36, 48, -48, 64, -64):
            x = start[0] + dx * fraction + nx * normal_offset
            y = start[1] + dy * fraction + ny * normal_offset
            box = draw.textbbox(
                (x, y), text, font=font, anchor="mm")
            boundary_penalty = (
                max(12 - box[0], 0)
                + max(12 - box[1], 0)
                + max(box[2] - (width - 12), 0)
                + max(box[3] - (height - 12), 0)
            ) * 1000
            overlap_penalty = sum(
                100000 for other in occupied
                if _boxes_overlap(box, other, padding=5)
            )
            obstruction_penalty = sum(
                15000 for other in avoid_boxes
                if _boxes_overlap(box, other, padding=3)
            )
            score = (
                boundary_penalty + overlap_penalty + obstruction_penalty
                + abs(normal_offset - 16) * 3
                + abs(fraction - 0.52) * 50
            )
            candidates.append((score, x, y, box))
    _, x, y, box = min(candidates, key=lambda item: item[0])
    occupied.append(box)
    return x, y, box


def select_groups(records, max_residues):
    # Never suppress direct interaction evidence to make the figure fit.
    groups = defaultdict(list)
    for record in records:
        groups[(record["chain"], record["name"], record["number"])].append(record)
    ranked = sorted(
        groups.items(),
        key=lambda pair: (
            min(PRIORITY[item["type"]] for item in pair[1]),
            min(item["distance"] for item in pair[1]),
            pair[0],
        ))
    direct = [(key, items) for key, items in ranked
              if any(item["type"] != "pocket" for item in items)]
    context = [(key, items) for key, items in ranked
               if all(item["type"] == "pocket" for item in items)]
    return direct + (context[:max_residues] if max_residues > 0 else context)


def render_png(ligand, records, output_path, width, height,
               max_residues, opaque_background=False,
               draw_hydrophobic_lines=False, _quality_attempt=0):
    drawing_mol, source_to_drawing = _prepare_ligand(ligand)
    ligand_scale = (1.0, 0.84, 0.70)[min(_quality_attempt, 2)]
    points = _layout_ligand(drawing_mol, width, height, ligand_scale)
    background = (
        (255, 255, 255, 255) if opaque_background
        else (255, 255, 255, 0)
    )
    image = Image.new("RGBA", (width, height), background)
    draw = ImageDraw.Draw(image)
    _draw_ligand(draw, drawing_mol, points)

    center = (width / 2, height / 2)
    placed = []
    for key, items in select_groups(records, max_residues):
        atom_indices = sorted({
            atom_index for item in items
            for atom_index in item["ligand_atoms"]
            if atom_index in source_to_drawing
        })
        anchors = [
            points[source_to_drawing[index]] for index in atom_indices]
        anchor = (
            (
                sum(point[0] for point in anchors) / len(anchors),
                sum(point[1] for point in anchors) / len(anchors),
            )
            if anchors else center
        )
        placed.append({
            "key": key,
            "items": items,
            "anchor": anchor,
        })

    key_residue_font = _font(29, bold=True)
    context_residue_font = _font(21)
    distance_font = _font(17)
    _radial_label_layout(
        draw, placed, points, center, width, height,
        key_residue_font, context_residue_font)

    atom_boxes = [
        (x - 22, y - 22, x + 22, y + 22) for x, y in points.values()
    ]
    ligand_box = (
        min(x for x, _ in points.values()) - 18,
        min(y for _, y in points.values()) - 18,
        max(x for x, _ in points.values()) + 18,
        max(y for _, y in points.values()) + 18,
    )
    residue_label_boxes = []
    for entry in placed:
        label = entry["items"][0]["label"]
        residue_label_boxes.append(draw.textbbox(
            (entry["label_x"], entry["label_y"]),
            label, font=entry["font"], anchor="mm"))
    occupied_distance_boxes = []
    pending_distance_labels = []

    def item_anchor(item):
        anchors = [
            points[source_to_drawing[index]]
            for index in item["ligand_atoms"]
            if index in source_to_drawing
        ]
        if not anchors:
            return center
        return (
            sum(point[0] for point in anchors) / len(anchors),
            sum(point[1] for point in anchors) / len(anchors),
        )

    fan_boxes = []
    for entry in placed:
        label_x = entry["label_x"]
        label_y = entry["label_y"]
        label = entry["items"][0]["label"]
        text_box = draw.textbbox(
            (0, 0), label, font=entry["font"], anchor="mm")
        hydrophobic_item = next(
            (item for item in entry["items"]
             if item["type"] == "hydrophobic"), None)
        entry["fan_center"] = None
        entry["fan_angle"] = None
        if hydrophobic_item:
            hydrophobic_anchor = item_anchor(hydrophobic_item)
            dx = hydrophobic_anchor[0] - label_x
            dy = hydrophobic_anchor[1] - label_y
            direction_length = max(math.hypot(dx, dy), 1.0)
            ux, uy = dx / direction_length, dy / direction_length
            hydrophobic_line_end = _line_end_at_label(
                hydrophobic_anchor,
                (label_x, label_y),
                text_box[2] - text_box[0],
                text_box[3] - text_box[1],
            )
            entry["fan_angle"] = math.atan2(dy, dx)
            entry["fan_center"] = (
                # Keep the fan clear of its own residue label as well as of
                # neighbouring labels; the quality gate verifies this later.
                hydrophobic_line_end[0] + ux * 64,
                hydrophobic_line_end[1] + uy * 64,
            )
            fan_x, fan_y = entry["fan_center"]
            fan_boxes.append((
                fan_x - 48, fan_y - 48, fan_x + 48, fan_y + 48))

    for entry in placed:
        label_x = entry["label_x"]
        label_y = entry["label_y"]
        label = entry["items"][0]["label"]
        residue_font = entry["font"]
        text_box = draw.textbbox(
            (0, 0), label, font=residue_font, anchor="mm")
        line_items = [
            item for item in entry["items"]
            if item["type"] != "pocket"
            and (
                item["type"] != "hydrophobic"
                or draw_hydrophobic_lines
            )
        ]
        fan_center = entry["fan_center"]
        angle_to_ligand = entry["fan_angle"]

        for index, item in enumerate(line_items):
            color = STYLE[item["type"]][1]
            offset = (index - (len(line_items) - 1) / 2) * 7
            raw_start = item_anchor(item)
            if fan_center is not None:
                fan_dx = fan_center[0] - raw_start[0]
                fan_dy = fan_center[1] - raw_start[1]
                fan_distance = max(math.hypot(fan_dx, fan_dy), 1.0)
                raw_end = (
                    fan_center[0] - fan_dx / fan_distance * 34,
                    fan_center[1] - fan_dy / fan_distance * 34,
                )
            else:
                raw_end = _line_end_at_label(
                    raw_start,
                    (label_x, label_y),
                    text_box[2] - text_box[0],
                    text_box[3] - text_box[1],
                )
            dx = raw_end[0] - raw_start[0]
            dy = raw_end[1] - raw_start[1]
            line_length = max(math.hypot(dx, dy), 1.0)
            nx, ny = -dy / line_length, dx / line_length
            start = (
                raw_start[0] + nx * offset,
                raw_start[1] + ny * offset,
            )
            end = (
                raw_end[0] + nx * offset,
                raw_end[1] + ny * offset,
            )
            if item["type"] == "hydrophobic":
                _dashed_line(
                    draw, start, end, color, width=2, dash=5, gap=5)
            else:
                _dashed_line(draw, start, end, color)
            distance_text = f"{item['distance']:.2f}"
            distance_x, distance_y, distance_box = _place_distance_label(
                draw,
                distance_text,
                start,
                end,
                occupied_distance_boxes,
                atom_boxes + residue_label_boxes + fan_boxes + [ligand_box],
                distance_font,
                width,
                height,
            )
            pending_distance_labels.append({
                "text": distance_text,
                "x": distance_x,
                "y": distance_y,
                "box": distance_box,
                "color": color,
            })

        if fan_center is not None:
            _hydrophobic_fan(draw, fan_center, angle_to_ligand)

        label_color = "#d7191c" if entry["is_key"] else "#2468b4"
        draw.text(
            (label_x, label_y),
            label,
            fill=label_color,
            font=residue_font,
            anchor="mm",
        )

    for distance_label in []:  # distance labels hidden for complete interaction maps
        box = distance_label["box"]
        draw.rectangle(
            (box[0] - 3, box[1] - 2, box[2] + 3, box[3] + 2),
            fill=background,
        )
        draw.text(
            (distance_label["x"], distance_label["y"]),
            distance_label["text"],
            fill=distance_label["color"],
            font=distance_font,
            anchor="mm",
        )

    layout_report = _verify_annotation_layout(
        residue_label_boxes,
        [],  # numeric distance labels are intentionally hidden
        ligand_box,
        fan_boxes,
        width,
        height,
    )
    if not layout_report["passed"]:
        if _quality_attempt < 4:
            return render_png(
                ligand, records, output_path,
                math.ceil(width * 1.35), math.ceil(height * 1.35),
                max_residues, opaque_background, draw_hydrophobic_lines,
                _quality_attempt + 1,
            )
        raise ValueError(
            "2D layout quality check failed after automatic relayout: "
            + "; ".join(layout_report["errors"]))

    stereo_report = verify_stereochemistry(ligand, drawing_mol)
    if not stereo_report["passed"]:
        raise ValueError(
            "2D stereochemistry verification failed before export: "
            + "; ".join(stereo_report["errors"]))

    quality_report = {
        "passed": True,
        "layout": layout_report,
        "stereochemistry": stereo_report,
        "canvas": {"width": width, "height": height},
        "layout_attempt": _quality_attempt + 1,
    }
    image.save(output_path, dpi=(300, 300))
    return quality_report


def write_spec(records, path):
    selected = {}
    for record in records:
        if record["type"] in {"hydrophobic", "pocket"}:
            continue
        key = (record["chain"], record["name"], record["number"])
        current = selected.get(key)
        if current is None or PRIORITY[record["type"]] < PRIORITY[current]:
            selected[key] = record["type"]
    parts = []
    for (chain, name, number), interaction_type in sorted(
            selected.items(), key=lambda item: (item[0][0], item[0][2])):
        residue = f"{chain}/{name}/{number}" if chain else f"{name} {number}"
        parts.append(f"{residue} {interaction_type}")
    spec = ", ".join(parts)
    Path(path).write_text(spec + "\n", encoding="utf-8")
    return spec


def main(argv=None):
    args = parse_args(argv)
    try:
        if args.pose_index < 1:
            raise ValueError("--pose-index must be 1 or greater")
        protein = load_protein(args.protein)
        poses = load_ligand_poses(args.ligand)
        if args.pose_index > len(poses):
            raise ValueError(
                f"Pose {args.pose_index} requested, but only "
                f"{len(poses)} pose(s) were found")
        ligand = poses[args.pose_index - 1]
        records = detect_interactions(
            protein, ligand, pocket_cutoff=args.pocket_cutoff)
        if not records:
            raise ValueError("RDKit detected no interactions for this pose")

        selected_groups = select_groups(records, args.max_residues)
        selected_keys = {key for key, _ in selected_groups}
        all_keys = {(item["chain"], item["name"], item["number"]) for item in records}
        direct_keys = {
            key for key, items in select_groups(records, 0)
            if any(item["type"] != "pocket" for item in items)
        }
        omitted_keys = sorted(all_keys - selected_keys)
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        png_path = output_dir / "interaction_2d.png"
        json_path = output_dir / "interaction_2d.json"
        spec_path = output_dir / "interaction_spec.txt"
        quality_report = render_png(
            ligand, records, png_path, args.width, args.height,
            args.max_residues, args.opaque_background,
            args.draw_hydrophobic_lines)
        json_path.write_text(
            json.dumps({
                "protein": str(Path(args.protein).resolve()),
                "ligand": str(Path(args.ligand).resolve()),
                "pose_index": args.pose_index,
                "transparent_background": not args.opaque_background,
                "hydrophobic_lines": args.draw_hydrophobic_lines,
                "interactions": records,
                "coverage": {
                    "detected_residue_count": len(all_keys),
                    "displayed_residue_count": len(selected_keys),
                    "omitted_residue_count": len(omitted_keys),
                    "omitted_residues": [
                        {"chain": key[0], "name": key[1], "number": key[2]}
                        for key in omitted_keys
                    ],
                    "direct_interactions_complete": direct_keys.issubset(selected_keys),
                    "distance_labels_shown": False,
                },
                "quality_check": quality_report,
            }, indent=2),
            encoding="utf-8",
        )
        spec = write_spec(records, spec_path)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"2D PNG: {png_path}")
    print(f"Interaction JSON: {json_path}")
    print(f"PyMOL interaction spec: {spec_path}")
    print(f"3D render spec: {spec}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
