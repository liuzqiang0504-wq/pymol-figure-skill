"""
PyMOL rendering engine for interaction and macro figures.

Run inside PyMOL's own Python (e.g. D:\\PyMOL\\python.exe):
  pymol.exe pymol_render.py --input file.maegz --interactions "PHE 1703 pi-pi, TYR 1719 pi-pi"

Architecture:
  - Load structure 鈫?identify protein + ligand objects
  - Setup interaction scene (浣滅敤鍥? 鈫?render 3 angles
  - Reset 鈫?setup macro scene (瀹忚鍥? 鈫?render 3 angles
  - Save .pse session
"""

import argparse
import json
import math
import os
import sys
import textwrap

# ring_detection.py lives next to this file
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)
import ring_detection

# 鈹€鈹€ PyMOL imports 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
import pymol
from pymol import cmd, stored, util

pymol.finish_launching(['pymol', '-cq'])

# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
# Constants
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲

SAFE_PYMOL_COLORS = {
    'white', 'black', 'red', 'green', 'blue', 'yellow', 'cyan', 'magenta',
    'orange', 'pink', 'teal', 'gray', 'palecyan', 'lightblue', 'salmon',
    'slate', 'violet', 'purple', 'purpleblue', 'warmpink', 'lightmagenta', 'limon',
    'forest', 'firebrick', 'ruby', 'wheat', 'paleyellow', 'marine',
    'olive', 'smudge', 'tv_red', 'tv_green', 'tv_blue', 'tv_yellow',
    'tv_orange', 'deepsalmon', 'density', 'chocolate', 'brown', 'carbon',
}

# PyMOL CBA (color by atom) index for "set 3, #6" 鈥?pink carbon base
LIGAND_CBA_INDEX = 6
RESIDUE_CBA_INDEX = 11

# 鈹€鈹€ Protein ring lookup 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
PROTEIN_RING_ATOMS = ring_detection.AROMATIC_RING_ATOMS

# 鈹€鈹€ H-bond donor/acceptor tables 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
HBOND_DONORS = {
    'ARG': ['NH1', 'NH2'],
    'LYS': ['NZ'],
    'ASN': ['ND2'],
    'GLN': ['NE2'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TYR': ['OH'],
    'HIS': ['ND1', 'NE2'],
    'TRP': ['NE1'],
}

HBOND_ACCEPTORS = {
    'ASP': ['OD1', 'OD2'],
    'GLU': ['OE1', 'OE2'],
    'ASN': ['OD1'],
    'GLN': ['OE1'],
}

HBOND_BACKBONE_DONOR = 'N'
HBOND_BACKBONE_ACCEPTOR = 'O'

# Geometry cutoffs
PI_PI_MAX_DIST = 4.5       # ring centroid distance (A)
PI_PI_MAX_OFFSET = 2.0     # ring plane offset (A)
HBOND_MAX_DIST = 3.5       # donor-acceptor distance (A)
HBOND_MIN_ANGLE = 120.0    # D-H...A angle (degrees)
METAL_CUTOFF = 3.0         # metal-donor distance (A)
CATION_PI_MAX_DIST = 6.0   # cation-ring centroid distance (A)
SALT_BRIDGE_MAX_DIST = 4.0 # charged group distance (A)
CONTACT_MAX_DIST = 4.0     # close heavy-atom contact distance (A)
POCKET_RADIUS = 5.0        # surface around ligand (A)

# 鈹€鈹€ Metal residues 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
METAL_RESIDUES = {
    'ZN', 'MG', 'CA', 'MN', 'FE', 'CO', 'NI', 'CU', 'CD', 'HG',
    'NA', 'K', 'PT', 'AU', 'AG', 'RU', 'OS', 'PD', 'IR', 'RH',
}

# 鈹€鈹€ Salt bridge charged atoms 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
SALT_BRIDGE_ACIDIC = {
    'ASP': ['OD1', 'OD2'],
    'GLU': ['OE1', 'OE2'],
}

SALT_BRIDGE_BASIC = {
    'ARG': ['NH1', 'NH2', 'NE'],
    'LYS': ['NZ'],
    'HIS': ['ND1', 'NE2'],
}

# 鈹€鈹€ Cation-pi cation sources (residue-side) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
CATION_PI_CATIONS = {
    'ARG': ['NH1', 'NH2', 'NE', 'CZ'],
    'LYS': ['NZ'],
    'HIS': ['ND1', 'NE2'],
}

# 鈹€鈹€ 300 DPI ray trace 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
RAY_WIDTH = 2400
RAY_HEIGHT = 1800

DEFAULT_STYLE = {
    "background": "white",
    "render": {
        "ray_shadow": "off",
        "ray_trace_mode": 1,
        "ray_trace_gain": 0.2,
        "ambient": 0.5,
        "antialias": 2,
    },
    "protein_cartoon_color": "palecyan",
    "hide_interacting_residue_cartoon": True,
    "interaction_zoom_buffer": 1.6,
    "ligand_cba_index": LIGAND_CBA_INDEX,
    "residue_cba_index": RESIDUE_CBA_INDEX,
    "dash": {
        "hbond": {"color": "green", "gap": 0.3, "width": 2.5},
        "pi-pi": {"color": "yellow", "gap": 0.3, "width": 2.5},
        "metal": {"color": "purple", "gap": 0.3, "width": 2.5},
        "cation-pi": {"color": "purple", "gap": 0.3, "width": 2.5},
        "salt-bridge": {"color": "purple", "gap": 0.3, "width": 2.5},
        "contact": {"color": "gray", "gap": 0.3, "width": 2.0},
    },
    "labels": {
        "font_paths": ["C:/Windows/Fonts/arial.ttf", "arial.ttf", "Arial.ttf"],
        "font_size": 60,
        "color": [0, 0, 0, 255],
        "outline_width": 0,
        "outline_color": [255, 255, 255, 255],
        "collision_padding": 10,
        "max_shift_px": 220,
        "offset_distance": 5.0,
    },
    "interaction_angles": [
        [-10, 0, "interaction_1.png"],
        [-15, 60, "interaction_2.png"],
        [-10, 120, "interaction_3.png"],
        [-20, 180, "interaction_4.png"],
        [-10, 240, "interaction_5.png"],
        [60, 0, "interaction_6.png"],
    ],
    "macro": {
        "protein_cartoon_color": "palecyan",
        "surface_color": "lightblue",
        "surface_transparency": 0.5,
        "ambient": 0.7,
        "ray_trace_gain": 0.08,
        "two_sided_lighting": 1,
        "paired_to_interaction": True,
        "paired_zoom_buffer": 10.0,
        "auto_angles": True,
        "keep_angle_previews": False,
        "preview_size": [600, 450],
        "candidate_angles": [
            [-20, 0], [-10, 25], [0, 55], [10, 85],
            [-15, 120], [5, 150], [15, 185], [-5, 215],
            [10, 250], [-15, 285], [0, 320], [20, 345],
            [25, 40], [-25, 140], [25, 220], [-25, 315],
        ],
        "angles": [
            [-5, 0, "macro_1.png"],
            [10, 60, "macro_2.png"],
            [-5, 120, "macro_3.png"],
            [10, 180, "macro_4.png"],
            [-5, 240, "macro_5.png"],
            [15, 300, "macro_6.png"],
        ],
    },
}

STYLE = DEFAULT_STYLE
INTERACTION_VIEWS = []


# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
# Argument parsing
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲

def parse_args():
    p = argparse.ArgumentParser(
        description='PyMOL interaction + macro figure generator')
    p.add_argument('--input', required=True, help='Structure file (.maegz, .pdb, etc.)')
    p.add_argument('--output', default=None,
                   help='Output directory (default: <input_dir>/pymol_figures/)')
    p.add_argument('--interactions', required=True,
                   help='Comma-separated interactions: "PHE 1703 pi-pi, ASP 189 hbond"')
    p.add_argument('--include-close-contacts', action='store_true',
                   help='Explicitly allow gray contact dashes (disabled by default)')
    p.add_argument('--style-profile', default=None,
                   help='Optional JSON style profile derived from a reference PyMOL image')
    return p.parse_args()


def _merge_dict(base, override):
    """Recursively merge a user style profile onto the default style."""
    result = dict(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _merge_dict(result[key], value)
        else:
            result[key] = value
    return result


def load_style_profile(path):
    if not path:
        return DEFAULT_STYLE
    with open(path, 'r', encoding='utf-8') as handle:
        user_style = json.load(handle)
    return _merge_dict(DEFAULT_STYLE, user_style)


def apply_close_contact_policy(interactions, include_close_contacts=False):
    """Keep gray close contacts only after an explicit renderer opt-in."""
    contacts = [item for item in interactions if item.get('type') == 'contact']
    if not contacts or include_close_contacts:
        return interactions

    print(
        f"NOTE: skipped {len(contacts)} close-contact interaction(s). "
        "Gray contact dashes are disabled by default; add "
        "--include-close-contacts only when the user explicitly requests them."
    )
    return [item for item in interactions if item.get('type') != 'contact']


def _safe_color(name, fallback='white'):
    if not name:
        return fallback
    if isinstance(name, str) and name.lower() in SAFE_PYMOL_COLORS:
        return name.lower()
    print(f"WARNING: unsupported PyMOL color {name!r}; using {fallback!r}")
    return fallback


def _apply_background():
    bg = STYLE.get('background', 'white')
    if isinstance(bg, str) and bg.lower() in ('transparent', 'none', 'alpha'):
        cmd.bg_color('white')
        cmd.set('ray_opaque_background', 0)
    else:
        cmd.set('ray_opaque_background', 1)
        cmd.bg_color(_safe_color(bg, 'white'))


def _apply_render_settings():
    render = STYLE.get('render', {})
    cmd.set('ray_shadow', render.get('ray_shadow', 'off'))
    cmd.set('ray_trace_mode', int(render.get('ray_trace_mode', 1)))
    cmd.set('ray_trace_gain', float(render.get('ray_trace_gain', 0.2)))
    cmd.set('ambient', float(render.get('ambient', 0.5)))
    cmd.set('antialias', int(render.get('antialias', 2)))


def _apply_macro_render_settings():
    macro_style = STYLE.get('macro', {})
    cmd.set('ambient', float(macro_style.get('ambient', 0.7)))
    cmd.set('ray_trace_gain', float(macro_style.get('ray_trace_gain', 0.08)))
    cmd.set('two_sided_lighting', int(macro_style.get('two_sided_lighting', 1)))


def _set_dash_style(dname, interaction_type):
    style = STYLE.get('dash', {}).get(interaction_type, {})
    cmd.set('dash_color', _safe_color(style.get('color'), 'purple'), dname)
    cmd.set('dash_gap', float(style.get('gap', 0.3)), dname)
    cmd.set('dash_width', float(style.get('width', 2.5)), dname)


def parse_interactions(spec_string):
    """
    Parse interaction spec string.

    Supported formats:
      "PHE 1703 pi-pi"
      "A/PHE/1703 pi-pi" (chain-aware)
      "PHE 1703 pi-pi, ASP 189 hbond" (multiple, comma-separated)

    Returns list of dicts: [{resn, chain, resi, type}, ...]
    """
    interactions = []
    if not spec_string or not spec_string.strip():
        return interactions

    for part in spec_string.split(','):
        part = part.strip()
        if not part:
            continue
        tokens = part.split()

        # Handle two formats:
        #   "PHE 1703 pi-pi"        鈫?3 tokens: resn, resi, type
        #   "A/PHE/1703 pi-pi"      鈫?2 tokens: chain/resn/resi, type
        if len(tokens) == 2 and '/' in tokens[0]:
            # Chain-aware shorthand: first token = "chain/resn/resi"
            res_id = tokens[0]
            itype_raw = tokens[1].lower()
        elif len(tokens) >= 3:
            res_id = tokens[0]   # "PHE" or "A/PHE/1703"
            resi = tokens[1]     # "1703" (ignored if chain-aware)
            itype_raw = tokens[2].lower()
        else:
            print(f"WARNING: skipping malformed interaction spec: {part!r}")
            continue

        # Normalize interaction type
        if itype_raw in ('pi-pi', 'pi-stacking', 'pipi', 'pi_pi'):
            itype = 'pi-pi'
        elif itype_raw in ('hbond', 'h-bond', 'hb', 'hydrogen', 'hydrogenbond'):
            itype = 'hbond'
        elif itype_raw in ('metal', 'metal-coord', 'metalcoord', 'coordination'):
            itype = 'metal'
        elif itype_raw in ('cation-pi', 'cationpi', 'cation_pi', 'cation'):
            itype = 'cation-pi'
        elif itype_raw in ('salt-bridge', 'saltbridge', 'salt_bridge', 'salt'):
            itype = 'salt-bridge'
        elif itype_raw in ('contact', 'close-contact', 'close_contact', 'hydrophobic',
                           'hydrophobic-contact', 'hydrophobic_contact'):
            itype = 'contact'
        else:
            print(f"WARNING: unknown interaction type {itype_raw!r}, treating as contact")
            itype = 'contact'

        # Parse chain-aware format: A/PHE/1703
        if '/' in res_id:
            chain_parts = res_id.split('/')
            if len(chain_parts) == 3:
                chain, resn, resi = chain_parts
            else:
                print(f"WARNING: malformed chain spec: {res_id!r}")
                continue
        else:
            resn = res_id
            chain = ''

        interactions.append({
            'resn': resn.upper(),
            'chain': chain,
            'resi': resi,
            'type': itype,
        })

    return interactions


# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
# Structure loading
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲

def load_structure(input_path):
    """Load structure file. Returns (protein_obj, ligand_obj)."""
    base = os.path.basename(input_path)
    name_no_ext = os.path.splitext(base)[0]
    # Handle double extensions like .maegz
    if name_no_ext.endswith('.mae'):
        name_no_ext = os.path.splitext(name_no_ext)[0]

    cmd.load(input_path)
    objects = cmd.get_object_list()
    print(f"Loaded {len(objects)} object(s): {objects}")

    if len(objects) == 0:
        raise RuntimeError(f"No objects loaded from {input_path}")

    # For .maegz files: protein is the first object (more atoms), ligand is second
    if len(objects) >= 2:
        # Sort by atom count: larger = protein, smaller = ligand
        atom_counts = [(o, cmd.count_atoms(o)) for o in objects]
        atom_counts.sort(key=lambda x: x[1], reverse=True)
        protein_obj = atom_counts[0][0]
        ligand_obj = atom_counts[1][0] if len(atom_counts) > 1 else None
    else:
        protein_obj = objects[0]
        # Try to find ligand by organic selection
        cmd.select('_lig_tmp', f'{protein_obj} and organic')
        if cmd.count_atoms('_lig_tmp') > 0:
            cmd.create('ligand', '_lig_tmp')
            ligand_obj = 'ligand'
            protein_obj = objects[0]
            cmd.delete('_lig_tmp')
        else:
            ligand_obj = None
            cmd.delete('_lig_tmp')

    print(f"Protein: {protein_obj}  |  Ligand: {ligand_obj}")
    return protein_obj, ligand_obj


# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
# Geometry helpers
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲

def get_atom_coords(selection):
    """Return list of [x, y, z] for all atoms matching selection."""
    stored.tmp = []
    cmd.iterate_state(1, selection, 'stored.tmp.append((x, y, z))')
    result = [list(c) for c in stored.tmp]
    stored.tmp = []
    return result


def get_protein_ring_centroid(obj, chain, resi, resn):
    """Compute centroid of the aromatic ring in a protein residue."""
    ring_atoms = PROTEIN_RING_ATOMS.get(resn.upper())
    if not ring_atoms:
        print(f"WARNING: {resn} is not a known aromatic residue")
        return None

    atom_list = '+'.join(ring_atoms)
    if chain:
        sel = f'{obj} and chain {chain} and resi {resi} and name {atom_list}'
    else:
        sel = f'{obj} and resi {resi} and name {atom_list}'

    coords = get_atom_coords(sel)
    if not coords:
        print(f"WARNING: no ring atoms found for {resn} {resi} (sel: {sel})")
        return None

    return ring_detection.centroid(coords)


def find_ligand_rings(ligand_obj):
    """Detect all rings in the ligand and return their centroids + atom indices."""
    stored.tmp = []
    cmd.iterate_state(1, f'{ligand_obj}',
                      'stored.tmp.append((elem, x, y, z, ID, name))')
    atom_data = []
    for elem, x, y, z, aid, aname in stored.tmp:
        atom_data.append({
            'elem': elem,
            'x': x, 'y': y, 'z': z,
            'idx': aid - 1,   # PyMOL IDs are 1-based
            'name': aname,
        })
    stored.tmp = []

    rings = ring_detection.detect_rings(atom_data)
    print(f"Ligand rings detected: {len(rings)}")

    # Build ring 鈫?centroid + atom names mapping
    ring_info = []
    for ring_indices in rings:
        positions = [[atom_data[i]['x'], atom_data[i]['y'], atom_data[i]['z']]
                     for i in ring_indices]
        ring_centroid = ring_detection.centroid(positions)
        ring_atom_ids = [atom_data[i]['idx'] + 1 for i in ring_indices]
        ring_info.append({
            'indices': ring_indices,        # 0-based in atom_data
            'atom_ids': ring_atom_ids,      # 1-based PyMOL IDs
            'centroid': ring_centroid,
            'size': len(ring_indices),
        })
    return ring_info


def distance(a, b):
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def find_closest_ligand_ring(protein_centroid, ligand_rings):
    """Return the ligand ring closest to the protein ring centroid."""
    if not ligand_rings:
        return None
    best = None
    best_dist = float('inf')
    for ring in ligand_rings:
        d = distance(protein_centroid, ring['centroid'])
        if d < best_dist:
            best_dist = d
            best = ring
    return best, best_dist


# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
# Scene builders
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲

def setup_interaction_scene(protein_obj, ligand_obj, interactions, output_dir):
    """Build the 浣滅敤鍥?(interaction close-up) scene."""
    global INTERACTION_VIEWS
    INTERACTION_VIEWS = []

    # IRON RULE 2: hide everything FIRST, then show
    cmd.hide('everything')

    # 鈹€鈹€ Protein cartoon 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.show('cartoon', protein_obj)
    # IRON RULE 1: color_deep for Maestro files
    cmd.color_deep(_safe_color(STYLE.get('protein_cartoon_color'), 'palecyan'), protein_obj, 0)

    # 鈹€鈹€ Global settings 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    _apply_background()
    cmd.set('stick_radius', 0.2)
    cmd.set('sphere_scale', 0.3)
    cmd.set('valence', 1)
    cmd.set('specular', 0.2)
    cmd.set('depth_cue', 0)

    # 鈹€鈹€ Ray trace settings 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    _apply_render_settings()

    # 鈹€鈹€ Label settings 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.set('label_color', 'black')
    cmd.set('label_size', 24)
    cmd.set('label_z_target', 0)    # always render on top, never occluded

    # 鈹€鈹€ Ligand 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.show('sticks', ligand_obj)
    # Color by element (CBA index 5274 鈥?pink carbon base, "set 3 #6")
    util.cba(int(STYLE.get('ligand_cba_index', LIGAND_CBA_INDEX)), ligand_obj, _self=cmd)

    # Hide non-polar hydrogens
    cmd.select('_lig_h_np', f'{ligand_obj} and elem H and (neighbor elem C)')
    cmd.hide('everything', '_lig_h_np')

    # Polar hydrogens in white sticks
    cmd.select('_lig_h_p', f'{ligand_obj} and elem H and not (neighbor elem C)')
    cmd.show('sticks', '_lig_h_p')
    cmd.color('white', '_lig_h_p')

    # 鈹€鈹€ Delete previous distance objects (IRON RULE 5) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    for dname in cmd.get_names('objects'):
        if any(dname.startswith(p) for p in ('pi_', 'hb_', 'metal_', 'catpi_', 'sb_', 'contact_')):
            cmd.delete(dname)

    # 鈹€鈹€ Interacting residue sticks 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    residue_names = []
    for ix, inter in enumerate(interactions):
        resn, chain, resi, itype = inter['resn'], inter['chain'], inter['resi'], inter['type']

        # Collective selection for all interacting residues
        if chain:
            res_sel = f'{protein_obj} and chain {chain} and resi {resi}'
        else:
            res_sel = f'{protein_obj} and resi {resi}'

        if STYLE.get('hide_interacting_residue_cartoon', True):
            cmd.hide('cartoon', res_sel)
        cmd.show('sticks', res_sel)
        util.cba(int(STYLE.get('residue_cba_index', RESIDUE_CBA_INDEX)), res_sel, _self=cmd)
        cmd.hide('everything', f'({res_sel}) and elem H and (neighbor elem C)')
        residue_names.append(f'{resn} {resi}')

        if itype == 'pi-pi':
            _setup_pi_interaction(protein_obj, ligand_obj, resn, chain, resi, ix)
        elif itype == 'hbond':
            _setup_hbond_interaction(protein_obj, ligand_obj, resn, chain, resi, ix)
        elif itype == 'metal':
            _setup_metal_interaction(protein_obj, ligand_obj, resn, chain, resi, ix)
        elif itype == 'cation-pi':
            _setup_cation_pi_interaction(protein_obj, ligand_obj, resn, chain, resi, ix)
        elif itype == 'salt-bridge':
            _setup_salt_bridge_interaction(protein_obj, ligand_obj, resn, chain, resi, ix)
        elif itype == 'contact':
            _setup_contact_interaction(protein_obj, ligand_obj, resn, chain, resi, ix)

    # 鈹€鈹€ Clean up temp selections 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.delete('_lig_h_np')
    cmd.delete('_lig_h_p')

    # 鈹€鈹€ Zoom to interacting region 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    # Collect all interacting residues for zoom
    if interactions:
        res_sels = []
        for inter in interactions:
            chain = inter['chain']
            resi = inter['resi']
            if chain:
                res_sels.append(f'({protein_obj} and chain {chain} and resi {resi})')
            else:
                res_sels.append(f'({protein_obj} and resi {resi})')
        zoom_sel = ' or '.join(res_sels)
    else:
        zoom_sel = ligand_obj
    interaction_zoom_buffer = float(STYLE.get('interaction_zoom_buffer', 1.6))
    cmd.zoom(f'{zoom_sel} or {ligand_obj}', interaction_zoom_buffer)

    # 鈹€鈹€ Collect label info for post-render compositing 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    label_data = []
    for name in cmd.get_names('objects'):
        if name.startswith('lbl_'):
            stored.tmp = []
            cmd.iterate_state(1, name, 'stored.tmp.append((x, y, z))')
            if stored.tmp:
                pos = stored.tmp[0]
            else:
                pos = (0, 0, 0)
            stored.tmp = []
            # Get label text
            stored.tmp = []
            cmd.iterate_state(1, name, 'stored.tmp.append(label)')
            txt = stored.tmp[0] if stored.tmp else ''
            stored.tmp = []
            label_data.append({'pos': pos, 'text': txt})
    # Hide labels before rendering
    for name in cmd.get_names('objects'):
        if name.startswith('lbl_'):
            cmd.hide('everything', name)

    # 鈹€鈹€ Render 3 varied angles 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    # 6 evenly-spaced views: every 60掳 azimuth, alternating tilt
    interaction_angles = STYLE.get('interaction_angles', DEFAULT_STYLE['interaction_angles'])
    for angle_idx, (turn_x, turn_y, fname) in enumerate(interaction_angles, 1):
        cmd.zoom(f'{zoom_sel} or {ligand_obj}', interaction_zoom_buffer)
        cmd.turn('x', turn_x)
        cmd.turn('y', turn_y)
        view = cmd.get_view()
        INTERACTION_VIEWS.append({
            'view': view,
            'interaction_file': fname,
            'macro_file': f'macro_{angle_idx}.png',
        })
        cmd.ray(RAY_WIDTH, RAY_HEIGHT)
        temp_path = os.path.join(output_dir, f'_temp_{angle_idx}.png')
        cmd.png(temp_path, dpi=300)
        final_path = os.path.join(output_dir, fname)

        # Get camera view matrix for 3D鈫?D projection
        view = cmd.get_view()
        _composite_labels_pil(temp_path, final_path, label_data, view,
                              RAY_WIDTH, RAY_HEIGHT)
        os.remove(temp_path)
        print(f"  Saved: {fname}")


def _world_to_screen(pos, view, width, height):
    """Project a 3D world position to 2D screen pixel coordinates using the
    PyMOL view matrix (18-element list from cmd.get_view()).

    PyMOL view matrix layout (0-indexed):
      [0..2]   = rotation row 0 (right)
      [3..5]   = rotation row 1 (up, inverted)
      [6..8]   = rotation row 2 (forward)
      [9..11]  = camera position (origin)
      [12..14] = rotation origin (scene center)
      [15]     = front clip
      [16]     = rear clip
      [17]     = orthoscopic flag
    """
    # Camera position
    cx, cy, cz = view[9], view[10], view[11]
    # Scene center (what the camera looks at)
    ox, oy, oz = view[12], view[13], view[14]

    # Forward vector (from camera toward scene center)
    fx = ox - cx
    fy = oy - cy
    fz = oz - cz
    fmag = math.sqrt(fx*fx + fy*fy + fz*fz)
    if fmag > 1e-6:
        fx, fy, fz = fx/fmag, fy/fmag, fz/fmag

    # Right vector (row 0 of rotation)
    rx, ry, rz = view[0], view[1], view[2]
    # Up vector (row 1 of rotation, PyMOL stores it inverted)
    ux, uy, uz = -view[3], -view[4], -view[5]

    # Vector from scene center to point
    px = pos[0] - ox
    py = pos[1] - oy
    pz = pos[2] - oz

    # Project onto camera basis
    x_cam = px*rx + py*ry + pz*rz   # horizontal displacement
    y_cam = px*ux + py*uy + pz*uz   # vertical displacement

    # Vector from camera to point (for depth check)
    vx = pos[0] - cx
    vy = pos[1] - cy
    vz = pos[2] - cz
    depth = vx*fx + vy*fy + vz*fz  # positive = in front of camera

    # Approximate field of view from zoom level
    slab = abs(view[16] - view[15])
    if slab < 1:
        slab = 20
    scale = height / slab

    screen_x = width / 2 + x_cam * scale
    screen_y = height / 2 + y_cam * scale

    return int(screen_x), int(screen_y), depth


def _label_python_commands():
    """Return Python command candidates for the Pillow label overlay step."""
    import shlex
    candidates = []
    for env_name in ("PYMOL_FIGURE_PYTHON", "PILLOW_PYTHON", "PYTHON"):
        value = os.environ.get(env_name)
        if value:
            candidates.append(shlex.split(value))
    candidates.extend([
        [sys.executable],
        ["py", "-3.12"],
        ["py", "-3.11"],
        ["python"],
        ["python3"],
    ])

    unique = []
    seen = set()
    for cmdline in candidates:
        if not cmdline:
            continue
        key = tuple(cmdline)
        if key in seen:
            continue
        seen.add(key)
        unique.append(cmdline)
    return unique


def _run_label_overlay_script(script):
    import subprocess
    errors = []
    for cmdline in _label_python_commands():
        try:
            subprocess.run(cmdline + ["-c", script], check=True,
                           capture_output=True, text=True, timeout=30)
            return True, " ".join(cmdline)
        except (subprocess.CalledProcessError, FileNotFoundError, OSError) as exc:
            detail = getattr(exc, "stderr", "") or str(exc)
            errors.append(f"{' '.join(cmdline)}: {detail.strip()}")
    return False, "; ".join(errors[-3:])


def _run_python_capture(script, timeout=30):
    import subprocess
    errors = []
    for cmdline in _label_python_commands():
        try:
            result = subprocess.run(
                cmdline + ["-c", script],
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
            return True, result.stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError, OSError) as exc:
            detail = getattr(exc, "stderr", "") or str(exc)
            errors.append(f"{' '.join(cmdline)}: {detail.strip()}")
    return False, "; ".join(errors[-3:])


def _score_macro_preview(image_path):
    """Score a low-resolution macro preview.

    Lower scores are better. The score penalizes large dark regions and rewards
    visible ligand-colored pixels, which tends to avoid back-facing or pocket-
    buried overview angles.
    """
    path_json = json.dumps(os.path.abspath(image_path))
    script = textwrap.dedent(f"""
from PIL import Image
import json

path = json.loads(r'''{path_json}''')
img = Image.open(path).convert('RGBA')
pixels = list(img.getdata())
visible = [p for p in pixels if p[3] > 0]
if not visible:
    print(json.dumps({{"score": 9999, "dark_fraction": 1.0, "ligand_fraction": 0.0}}))
    raise SystemExit
sample_step = max(1, len(visible) // 120000)
visible = visible[::sample_step]
dark = 0
near_black = 0
ligand_like = 0
for r, g, b, a in visible:
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    if lum < 70:
        dark += 1
    if lum < 35:
        near_black += 1
    if r > 170 and g > 140 and b < 90:
        ligand_like += 1
n = len(visible)
dark_fraction = dark / n
near_black_fraction = near_black / n
ligand_fraction = ligand_like / n
score = dark_fraction * 1000 + near_black_fraction * 1400 - ligand_fraction * 1800
print(json.dumps({{
    "score": score,
    "dark_fraction": dark_fraction,
    "near_black_fraction": near_black_fraction,
    "ligand_fraction": ligand_fraction,
}}))
""")
    ok, output = _run_python_capture(script, timeout=30)
    if not ok:
        return None
    try:
        return json.loads(output)
    except json.JSONDecodeError:
        return None


def _composite_labels_pil(temp_path, output_path, label_data, view, width, height):
    """Render labels onto the rendered image using PIL for always-on-top text.

    Projects each label's 3D position to 2D screen coordinates using the
    PyMOL camera matrix, then draws black text at those pixel positions.
    """
    import json
    # Prepare label screen positions
    screen_labels = []
    for lbl in label_data:
        sx, sy, depth = _world_to_screen(lbl['pos'], view, width, height)
        # Skip labels behind the camera
        if depth < 0:
            continue
        screen_labels.append({'x': sx, 'y': sy, 'text': lbl['text']})

    if not screen_labels:
        import shutil
        shutil.copy(temp_path, output_path)
        return

    label_style = STYLE.get('labels', {})
    labels_json = json.dumps(screen_labels)
    font_paths_json = json.dumps(label_style.get('font_paths', DEFAULT_STYLE['labels']['font_paths']))
    font_size = int(label_style.get('font_size', 60))
    font_color_json = json.dumps(label_style.get('color', [0, 0, 0, 255]))
    outline_width = int(label_style.get('outline_width', 0))
    outline_color_json = json.dumps(label_style.get('outline_color', [255, 255, 255, 255]))
    collision_padding = int(label_style.get('collision_padding', 10))
    max_shift_px = int(label_style.get('max_shift_px', 160))
    script = textwrap.dedent(fr"""
from PIL import Image, ImageFont, ImageDraw
import json

labels = json.loads(r'''{labels_json}''')
font_paths = json.loads(r'''{font_paths_json}''')
font_size = {font_size}
font_color = tuple(json.loads(r'''{font_color_json}'''))
outline_width = {outline_width}
outline_color = tuple(json.loads(r'''{outline_color_json}'''))
collision_padding = {collision_padding}
max_shift_px = {max_shift_px}
img = Image.open(r'{temp_path}').convert('RGBA')
draw = ImageDraw.Draw(img)

font = None
for font_path in font_paths:
    try:
        font = ImageFont.truetype(font_path, font_size)
        break
    except OSError:
        pass
if font is None:
    font = ImageFont.load_default()

def padded_bbox(x, y, txt):
    bbox = draw.textbbox((x, y), txt, font=font)
    return (
        bbox[0] - collision_padding,
        bbox[1] - collision_padding,
        bbox[2] + collision_padding,
        bbox[3] + collision_padding,
    )

def overlap_area(a, b):
    left = max(a[0], b[0])
    top = max(a[1], b[1])
    right = min(a[2], b[2])
    bottom = min(a[3], b[3])
    if right <= left or bottom <= top:
        return 0
    return (right - left) * (bottom - top)

def clamp_position(x, y, txt):
    bbox = draw.textbbox((x, y), txt, font=font)
    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    x = max(0, min(x, img.width - width - 1))
    y = max(0, min(y, img.height - height - 1))
    return x, y

def place_labels(labels):
    placed = []
    result = []
    for lbl in sorted(labels, key=lambda item: (item['y'], item['x'])):
        base_x, base_y, txt = int(lbl['x']), int(lbl['y']), lbl['text']
        raw_bbox = draw.textbbox((base_x, base_y), txt, font=font)
        label_w = raw_bbox[2] - raw_bbox[0]
        label_h = raw_bbox[3] - raw_bbox[1]
        step_x = max(45, int(label_w * 0.45))
        step_y = max(36, int(label_h + collision_padding * 2))
        offsets = [(0, 0)]
        max_radius = max(2, int(max_shift_px / max(1, min(step_x, step_y))) + 1)
        for radius in range(1, max_radius + 1):
            for gx in range(-radius, radius + 1):
                for gy in range(-radius, radius + 1):
                    if max(abs(gx), abs(gy)) != radius:
                        continue
                    dx = int(gx * step_x * 0.75)
                    dy = int(gy * step_y * 0.75)
                    travel = (dx * dx + dy * dy) ** 0.5
                    if 0 < travel <= max_shift_px:
                        offsets.append((dx, dy))
        offsets = sorted(set(offsets), key=lambda item: (item[0] * item[0] + item[1] * item[1]) ** 0.5)
        best = None
        for dx, dy in offsets:
            x, y = clamp_position(base_x + dx, base_y + dy, txt)
            bbox = padded_bbox(x, y, txt)
            overlap = sum(overlap_area(bbox, existing) for existing in placed)
            travel = (dx * dx + dy * dy) ** 0.5
            score = overlap * 5000 + travel * 12
            if best is None or score < best[0]:
                best = (score, x, y, bbox)
            if overlap == 0:
                break
        _, x, y, bbox = best
        placed.append(bbox)
        result.append({{'x': x, 'y': y, 'text': txt}})
    return result

for lbl in place_labels(labels):
    x, y, txt = lbl['x'], lbl['y'], lbl['text']
    if outline_width > 0:
        for dx in range(-outline_width, outline_width + 1):
            for dy in range(-outline_width, outline_width + 1):
                if dx == 0 and dy == 0:
                    continue
                draw.text((x + dx, y + dy), txt, fill=outline_color, font=font)
    draw.text((x, y), txt, fill=font_color, font=font)

img.save(r'{output_path}', dpi=(300, 300))
""")
    ok, detail = _run_label_overlay_script(script)
    if not ok:
        import shutil
        shutil.copy(temp_path, output_path)
        print('  (label overlay skipped - Pillow Python not found)')
        if detail:
            print(f'   tried: {detail}')


def _setup_pi_interaction(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Create 蟺-蟺 stacking dashes + residue label for one interaction."""
    # 1. Protein ring centroid
    prot_cent = get_protein_ring_centroid(protein_obj, chain, resi, resn)
    if not prot_cent:
        print(f"SKIP: cannot find ring centroid for {resn} {resi}")
        return

    # 2. Find ligand rings
    lig_rings = find_ligand_rings(ligand_obj)

    # 3. Find closest ligand ring to protein centroid
    result = find_closest_ligand_ring(prot_cent, lig_rings)
    if not result:
        print(f"SKIP: no ligand ring found for {resn} {resi} pi-pi")
        return
    closest_ring, dist = result
    print(f"  {resn} {resi} pi-pi: protein-ligand ring distance = {dist:.2f} A")

    # 4. Create centroid pseudoatoms
    pname = f'cent_{ix}'
    lname = f'cent_lig_{ix}'
    cmd.pseudoatom(pname, pos=prot_cent)
    cmd.pseudoatom(lname, pos=closest_ring['centroid'])
    cmd.set('sphere_scale', 0.32, pname)
    cmd.set('sphere_scale', 0.32, lname)
    cmd.color('yellow', pname)
    cmd.color('yellow', lname)

    # 5. Dashed line between centroids
    dname = f'pi_{ix}'
    cmd.distance(dname, pname, lname)
    _set_dash_style(dname, 'pi-pi')
    cmd.hide('labels', dname)
    # Hide centroid pseudoatoms 鈥?only the dashed line should be visible
    cmd.hide('everything', pname)
    cmd.hide('everything', lname)

    # 6. Residue label (radiates outward from ligand center)
    label_text = f'{resn} {resi}'
    label_pt = _label_offset(ligand_obj, prot_cent)
    lbl_name = f'lbl_{ix}'
    cmd.pseudoatom(lbl_name, pos=label_pt)
    cmd.label(lbl_name, f'"{label_text}"')
    # IRON RULE 3: hide spheres, NOT everything
    cmd.hide('spheres', lbl_name)


def _setup_hbond_interaction(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Create H-bond dashes for one interaction."""
    # Find donor/acceptor atoms in the residue
    donors = HBOND_DONORS.get(resn.upper(), [])
    acceptors = HBOND_ACCEPTORS.get(resn.upper(), [])

    # Also consider backbone
    if resn.upper() not in HBOND_DONORS and resn.upper() not in HBOND_ACCEPTORS:
        donors = [HBOND_BACKBONE_DONOR]
        acceptors = [HBOND_BACKBONE_ACCEPTOR]

    hbond_count = 0

    # Residue as H-bond donor 鈫?find acceptors in ligand
    for donor_atom in donors:
        if chain:
            d_sel = f'{protein_obj} and chain {chain} and resi {resi} and name {donor_atom}'
        else:
            d_sel = f'{protein_obj} and resi {resi} and name {donor_atom}'
        d_coords = get_atom_coords(d_sel)
        if not d_coords:
            continue

        # Find close N/O atoms in ligand
        a_coords_all = get_atom_coords(f'{ligand_obj} and (elem N or elem O)')
        for a_coord in a_coords_all:
            d = distance(d_coords[0], a_coord)
            if d <= HBOND_MAX_DIST:
                a_sel = f'{ligand_obj} and (elem N or elem O) and x > {a_coord[0]-0.1} and x < {a_coord[0]+0.1} and y > {a_coord[1]-0.1} and y < {a_coord[1]+0.1} and z > {a_coord[2]-0.1} and z < {a_coord[2]+0.1}'
                dname = f'hb_{ix}_{hbond_count}'
                cmd.distance(dname, d_sel, a_sel)
                _set_dash_style(dname, 'hbond')
                cmd.hide('labels', dname)
                hbond_count += 1

    # Residue as H-bond acceptor 鈫?find donors in ligand
    for acc_atom in acceptors:
        if chain:
            a_sel = f'{protein_obj} and chain {chain} and resi {resi} and name {acc_atom}'
        else:
            a_sel = f'{protein_obj} and resi {resi} and name {acc_atom}'
        a_coords = get_atom_coords(a_sel)
        if not a_coords:
            continue

        d_coords_all = get_atom_coords(f'{ligand_obj} and (elem N or elem O)')
        for d_coord in d_coords_all:
            d = distance(a_coords[0], d_coord)
            if d <= HBOND_MAX_DIST:
                d_sel = f'{ligand_obj} and (elem N or elem O) and x > {d_coord[0]-0.1} and x < {d_coord[0]+0.1} and y > {d_coord[1]-0.1} and y < {d_coord[1]+0.1} and z > {d_coord[2]-0.1} and z < {d_coord[2]+0.1}'
                dname = f'hb_{ix}_{hbond_count}'
                cmd.distance(dname, d_sel, a_sel)
                _set_dash_style(dname, 'hbond')
                cmd.hide('labels', dname)
                hbond_count += 1

    # Label
    prot_cent = None
    ring_atoms = PROTEIN_RING_ATOMS.get(resn.upper())
    if ring_atoms:
        prot_cent = get_protein_ring_centroid(protein_obj, chain, resi, resn)
    if not prot_cent:
        # Use CA atom as fallback
        if chain:
            ca_sel = f'{protein_obj} and chain {chain} and resi {resi} and name CA'
        else:
            ca_sel = f'{protein_obj} and resi {resi} and name CA'
        ca_coords = get_atom_coords(ca_sel)
        prot_cent = ca_coords[0] if ca_coords else [0, 0, 0]

    label_text = f'{resn} {resi}'
    label_pt = _label_offset(ligand_obj, prot_cent)
    lbl_name = f'lbl_{ix}'
    cmd.pseudoatom(lbl_name, pos=label_pt)
    cmd.label(lbl_name, f'"{label_text}"')
    cmd.hide('spheres', lbl_name)

    print(f"  {resn} {resi} H-bond: {hbond_count} contact(s)")


def _setup_contact_interaction(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Create gray dashes for close heavy-atom contacts."""
    if chain:
        res_sel = f'{protein_obj} and chain {chain} and resi {resi}'
    else:
        res_sel = f'{protein_obj} and resi {resi}'

    prot_sel = f'({res_sel}) and not elem H'
    lig_sel = f'({ligand_obj}) and not elem H'

    stored.tmp = []
    cmd.iterate_state(1, prot_sel, 'stored.tmp.append((ID, x, y, z))')
    prot_atoms = list(stored.tmp)
    stored.tmp = []
    cmd.iterate_state(1, lig_sel, 'stored.tmp.append((ID, x, y, z))')
    lig_atoms = list(stored.tmp)
    stored.tmp = []

    pairs = []
    for pid, px, py, pz in prot_atoms:
        p = [px, py, pz]
        for lid, lx, ly, lz in lig_atoms:
            d = distance(p, [lx, ly, lz])
            if d <= CONTACT_MAX_DIST:
                pairs.append((d, pid, lid))
    pairs.sort(key=lambda x: x[0])

    for j, (d, pid, lid) in enumerate(pairs[:4]):
        dname = f'contact_{ix}_{j}'
        cmd.distance(dname, f'{protein_obj} and id {pid}', f'{ligand_obj} and id {lid}')
        _set_dash_style(dname, 'contact')
        cmd.hide('labels', dname)

    if chain:
        ca_sel = f'{protein_obj} and chain {chain} and resi {resi} and name CA'
    else:
        ca_sel = f'{protein_obj} and resi {resi} and name CA'
    ca_coords = get_atom_coords(ca_sel)
    prot_cent = ca_coords[0] if ca_coords else (get_atom_coords(res_sel) or [[0, 0, 0]])[0]

    label_text = f'{resn} {resi}'
    label_pt = _label_offset(ligand_obj, prot_cent)
    lbl_name = f'lbl_{ix}'
    cmd.pseudoatom(lbl_name, pos=label_pt)
    cmd.label(lbl_name, f'"{label_text}"')
    cmd.hide('spheres', lbl_name)

    print(f"  {resn} {resi} contact: {min(len(pairs), 4)} close contact(s)")


def _setup_metal_interaction(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Create metal coordination dashes + metal sphere display."""
    # Find metal atom in protein
    if chain:
        metal_sel = f'{protein_obj} and chain {chain} and resi {resi}'
    else:
        metal_sel = f'{protein_obj} and resi {resi}'

    # Verify this is actually a metal (check residue name, most reliable)
    is_metal_res = resn.upper() in METAL_RESIDUES
    if not is_metal_res:
        # Also try element-level check
        stored.tmp = []
        cmd.iterate_state(1, metal_sel, 'stored.tmp.append(elem)')
        elements = set(stored.tmp)
        stored.tmp = []
        is_metal_res = any(e.upper() in METAL_RESIDUES for e in elements)

    if not is_metal_res:
        print(f"WARNING: {resn} {resi} does not appear to be a metal ion")

    # Show metal as sphere, color by chain
    cmd.show('spheres', metal_sel)
    cmd.set('sphere_scale', 0.3, metal_sel)
    util.cnc(metal_sel, _self=cmd)

    # Find coordinating atoms (N, O, S) in ligand within cutoff
    coord_sel = f'{ligand_obj} and (elem N or elem O or elem S)'
    cmd.select('_metal_coord', f'{coord_sel} within {METAL_CUTOFF} of {metal_sel}')

    if cmd.count_atoms('_metal_coord') == 0:
        print(f"  {resn} {resi} metal: no coordinating atoms found within {METAL_CUTOFF} A")
        cmd.delete('_metal_coord')
        # Still add label
        _add_metal_label(protein_obj, ligand_obj, resn, chain, resi, ix)
        return

    # Draw purple dashes from metal to each coordinating atom
    stored.tmp = []
    cmd.iterate_state(1, '_metal_coord', 'stored.tmp.append((ID, x, y, z))')
    coord_atoms = [(aid, [x, y, z]) for aid, x, y, z in stored.tmp]
    stored.tmp = []

    # Get metal position
    stored.tmp = []
    cmd.iterate_state(1, metal_sel, 'stored.tmp.append((x, y, z))')
    metal_coords = stored.tmp[0] if stored.tmp else [0, 0, 0]
    stored.tmp = []

    for j, (aid, coord) in enumerate(coord_atoms):
        d = distance(metal_coords, coord)
        dname = f'metal_{ix}_{j}'
        cmd.distance(dname, metal_sel, f'id {aid}')
        _set_dash_style(dname, 'metal')
        cmd.hide('labels', dname)

    cmd.delete('_metal_coord')
    print(f"  {resn} {resi} metal: {len(coord_atoms)} coordinating atom(s) found")

    # Add label
    _add_metal_label(protein_obj, ligand_obj, resn, chain, resi, ix)


def _add_metal_label(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Add a label for a metal ion, offset from its position."""
    if chain:
        metal_sel = f'{protein_obj} and chain {chain} and resi {resi}'
    else:
        metal_sel = f'{protein_obj} and resi {resi}'

    stored.tmp = []
    cmd.iterate_state(1, metal_sel, 'stored.tmp.append((x, y, z))')
    if not stored.tmp:
        return
    pos = stored.tmp[0]
    stored.tmp = []
    label_pt = _label_offset(ligand_obj, pos)
    lbl_name = f'lbl_{ix}'
    cmd.pseudoatom(lbl_name, pos=label_pt)
    cmd.label(lbl_name, f'"{resn} {resi}"')
    cmd.hide('spheres', lbl_name)


def _setup_cation_pi_interaction(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Create cation-pi dashes between a cation and an aromatic ring."""
    # Determine whether specified residue is the cation or the pi system
    is_metal = resn.upper() in METAL_RESIDUES
    is_cation_res = resn.upper() in CATION_PI_CATIONS
    is_aromatic = resn.upper() in PROTEIN_RING_ATOMS

    if is_metal or is_cation_res:
        # Specified residue is the cation 鈥?find nearby aromatic ring
        if chain:
            cation_sel = f'{protein_obj} and chain {chain} and resi {resi}'
        else:
            cation_sel = f'{protein_obj} and resi {resi}'

        # Get cation position
        stored.tmp = []
        cmd.iterate_state(1, cation_sel, 'stored.tmp.append((x, y, z))')
        if not stored.tmp:
            print(f"SKIP: cannot find cation for {resn} {resi}")
            return
        cation_pos = list(stored.tmp[0])
        stored.tmp = []

        # Find closest aromatic ring in protein or ligand
        found = _find_nearest_aromatic_ring(protein_obj, ligand_obj, cation_pos)
        if not found:
            print(f"SKIP: no aromatic ring found near {resn} {resi}")
            return

        ring_centroid, ring_label, ring_source = found
        print(f"  {resn} {resi} cation-pi: cation-aromatic distance = {distance(cation_pos, ring_centroid):.2f} A (with {ring_label})")

        # Pseudoatoms for dashes
        cname = f'cat_{ix}'
        rname = f'catpi_ring_{ix}'
        cmd.pseudoatom(cname, pos=cation_pos)
        cmd.pseudoatom(rname, pos=ring_centroid)
        cmd.set('sphere_scale', 0.32, cname)
        cmd.set('sphere_scale', 0.32, rname)
        cmd.color('purple', cname)
        cmd.color('purple', rname)

    elif is_aromatic:
        # Specified residue is the pi system 鈥?find nearby cation
        if chain:
            ring_sel_base = f'{protein_obj} and chain {chain} and resi {resi}'
        else:
            ring_sel_base = f'{protein_obj} and resi {resi}'

        # Get ring centroid
        ring_centroid = get_protein_ring_centroid(protein_obj, chain, resi, resn)
        if not ring_centroid:
            print(f"SKIP: cannot find ring centroid for {resn} {resi}")
            return

        # Find closest cation
        found = _find_nearest_cation(protein_obj, ligand_obj, ring_centroid)
        if not found:
            print(f"SKIP: no cation found near {resn} {resi}")
            return

        cation_pos, cation_label, cation_source = found
        print(f"  {resn} {resi} cation-pi: ring-cation distance = {distance(ring_centroid, cation_pos):.2f} A (with {cation_label})")

        # Pseudoatoms
        rname = f'catpi_ring_{ix}'
        cname = f'cat_{ix}'
        cmd.pseudoatom(rname, pos=ring_centroid)
        cmd.pseudoatom(cname, pos=cation_pos)
        cmd.set('sphere_scale', 0.32, rname)
        cmd.set('sphere_scale', 0.32, cname)
        cmd.color('purple', rname)
        cmd.color('purple', cname)

    else:
        print(f"SKIP: {resn} {resi} is neither a known cation nor aromatic residue")
        return

    # Draw purple dashes
    dname = f'catpi_{ix}'
    cmd.distance(dname, f'cat_{ix}', f'catpi_ring_{ix}')
    _set_dash_style(dname, 'cation-pi')
    cmd.hide('labels', dname)
    # Hide centroid pseudoatoms 鈥?only the dashed line should be visible
    cmd.hide('everything', f'cat_{ix}')
    cmd.hide('everything', f'catpi_ring_{ix}')

    # Show sphere for metal cation, color by chain
    if is_metal:
        cmd.show('spheres', cation_sel)
        cmd.set('sphere_scale', 0.3, cation_sel)
        util.cnc(cation_sel, _self=cmd)

    # Label
    label_text = f'{resn} {resi}'
    lbl_pt = [ring_centroid[0] + 2.5, ring_centroid[1], ring_centroid[2]]
    lbl_name = f'lbl_{ix}'
    cmd.pseudoatom(lbl_name, pos=lbl_pt)
    cmd.label(lbl_name, f'"{label_text}"')
    cmd.hide('spheres', lbl_name)


def _find_nearest_aromatic_ring(protein_obj, ligand_obj, ref_pos):
    """Find the aromatic ring (protein or ligand) closest to a reference position.
    Returns (centroid, label_string, source) or None."""
    best = None
    best_dist = float('inf')

    # Check protein aromatic residues
    for resn, ring_atoms in PROTEIN_RING_ATOMS.items():
        atom_list = '+'.join(ring_atoms)
        sel = f'{protein_obj} and resn {resn} and name {atom_list}'
        if cmd.count_atoms(sel) == 0:
            continue
        # Get residues matching this pattern
        stored.tmp = []
        cmd.iterate_state(1, sel, 'stored.tmp.append((resi, chain, x, y, z))')
        # Group by residue
        residues = {}
        for rsi, ch, x, y, z in stored.tmp:
            key = (ch, rsi)
            residues.setdefault(key, []).append([x, y, z])
        stored.tmp = []
        for (ch, rsi), positions in residues.items():
            if len(positions) < 3:
                continue
            c = ring_detection.centroid(positions)
            d = distance(ref_pos, c)
            if d < best_dist and d <= CATION_PI_MAX_DIST:
                best_dist = d
                best = (c, f'{resn} {rsi}', 'protein')
    return best


def _find_nearest_cation(protein_obj, ligand_obj, ref_pos):
    """Find the cation closest to a reference position.
    Returns (position, label_string, source) or None."""
    best = None
    best_dist = float('inf')

    # Check metal ions in protein
    metal_sel_parts = ' or '.join(f'elem {m}' for m in METAL_RESIDUES)
    sel = f'{protein_obj} and ({metal_sel_parts})'
    if cmd.count_atoms(sel) > 0:
        stored.tmp = []
        cmd.iterate_state(1, sel, 'stored.tmp.append((elem, resi, chain, x, y, z))')
        for elem, rsi, ch, x, y, z in stored.tmp:
            d = distance(ref_pos, [x, y, z])
            if d < best_dist and d <= CATION_PI_MAX_DIST:
                best_dist = d
                best = ([x, y, z], f'{elem} {rsi}', 'metal')
        stored.tmp = []

    # Check cationic residues (ARG, LYS, HIS)
    for resn, cat_atoms in CATION_PI_CATIONS.items():
        atom_list = '+'.join(cat_atoms)
        sel2 = f'{protein_obj} and resn {resn} and name {atom_list}'
        if cmd.count_atoms(sel2) == 0:
            continue
        stored.tmp = []
        cmd.iterate_state(1, sel2, 'stored.tmp.append((resi, chain, x, y, z))')
        residues = {}
        for rsi, ch, x, y, z in stored.tmp:
            key = (ch, rsi)
            residues.setdefault(key, []).append([x, y, z])
        stored.tmp = []
        for (ch, rsi), positions in residues.items():
            c = ring_detection.centroid(positions)
            d = distance(ref_pos, c)
            if d < best_dist and d <= CATION_PI_MAX_DIST:
                best_dist = d
                best = (c, f'{resn} {rsi}', 'cationic_residue')
    return best


def _setup_salt_bridge_interaction(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Create salt bridge dashes between a charged residue and the ligand."""
    is_acidic = resn.upper() in SALT_BRIDGE_ACIDIC
    is_basic = resn.upper() in SALT_BRIDGE_BASIC

    if not is_acidic and not is_basic:
        print(f"WARNING: {resn} {resi} is not a known charged residue for salt bridges")
        return

    if is_acidic:
        charged_atoms = SALT_BRIDGE_ACIDIC[resn.upper()]
        if chain:
            prot_sel_base = f'{protein_obj} and chain {chain} and resi {resi} and name {"+".join(charged_atoms)}'
        else:
            prot_sel_base = f'{protein_obj} and resi {resi} and name {"+".join(charged_atoms)}'
        # Find basic (N) atoms in ligand
        lig_counter_sel = f'{ligand_obj} and elem N'
        role = 'acidic'
    else:
        charged_atoms = SALT_BRIDGE_BASIC[resn.upper()]
        if chain:
            prot_sel_base = f'{protein_obj} and chain {chain} and resi {resi} and name {"+".join(charged_atoms)}'
        else:
            prot_sel_base = f'{protein_obj} and resi {resi} and name {"+".join(charged_atoms)}'
        # Find acidic (O) atoms in ligand
        lig_counter_sel = f'{ligand_obj} and elem O'
        role = 'basic'

    # Find close contacts
    cmd.select('_sb_contacts', f'{lig_counter_sel} within {SALT_BRIDGE_MAX_DIST} of {prot_sel_base}')
    if cmd.count_atoms('_sb_contacts') == 0:
        print(f"  {resn} {resi} salt-bridge ({role}): no complementary charged atoms within {SALT_BRIDGE_MAX_DIST} A")
        cmd.delete('_sb_contacts')
        _add_salt_bridge_label(protein_obj, ligand_obj, resn, chain, resi, ix)
        return

    # Get charged atom coords in protein
    stored.tmp = []
    cmd.iterate_state(1, prot_sel_base, 'stored.tmp.append((ID, x, y, z))')
    prot_atoms = [(aid, [x, y, z]) for aid, x, y, z in stored.tmp]
    stored.tmp = []

    # Get counter atom coords in ligand
    stored.tmp = []
    cmd.iterate_state(1, '_sb_contacts', 'stored.tmp.append((ID, x, y, z))')
    lig_atoms = [(aid, [x, y, z]) for aid, x, y, z in stored.tmp]
    stored.tmp = []

    bridge_count = 0
    for paid, pc in prot_atoms:
        for laid, lc in lig_atoms:
            d = distance(pc, lc)
            if d <= SALT_BRIDGE_MAX_DIST:
                dname = f'sb_{ix}_{bridge_count}'
                cmd.distance(dname, f'{protein_obj} and id {paid}', f'{ligand_obj} and id {laid}')
                _set_dash_style(dname, 'salt-bridge')
                cmd.hide('labels', dname)
                bridge_count += 1

    cmd.delete('_sb_contacts')
    print(f"  {resn} {resi} salt-bridge ({role}): {bridge_count} contact(s)")
    _add_salt_bridge_label(protein_obj, ligand_obj, resn, chain, resi, ix)


def _add_salt_bridge_label(protein_obj, ligand_obj, resn, chain, resi, ix):
    """Add a label for a salt bridge residue."""
    # Find a reference atom for label placement
    charged_atoms = SALT_BRIDGE_ACIDIC.get(resn.upper()) or SALT_BRIDGE_BASIC.get(resn.upper())
    if not charged_atoms:
        return
    if chain:
        sel = f'{protein_obj} and chain {chain} and resi {resi} and name {charged_atoms[0]}'
    else:
        sel = f'{protein_obj} and resi {resi} and name {charged_atoms[0]}'
    stored.tmp = []
    cmd.iterate_state(1, sel, 'stored.tmp.append((x, y, z))')
    if not stored.tmp:
        return
    pos = stored.tmp[0]
    stored.tmp = []
    label_pt = _label_offset(ligand_obj, pos)
    lbl_name = f'lbl_{ix}'
    cmd.pseudoatom(lbl_name, pos=label_pt)
    cmd.label(lbl_name, f'"{resn} {resi}"')
    cmd.hide('spheres', lbl_name)


def _get_ligand_centroid(ligand_obj):
    """Compute the geometric centroid of the ligand."""
    coords = get_atom_coords(ligand_obj)
    if not coords:
        return [0.0, 0.0, 0.0]
    n = len(coords)
    return [sum(c[i] for c in coords) / n for i in range(3)]


def _label_offset(ligand_obj, prot_cent):
    """Compute label position that radiates outward from ligand centroid.

    Places the label on the far side of the residue relative to the ligand,
    so labels naturally spread apart in different directions.

    Args:
        ligand_obj: PyMOL object name for the ligand
        prot_cent: [x, y, z] centroid of the protein residue

    Returns:
        [x, y, z] absolute label position
    """
    lig_cent = _get_ligand_centroid(ligand_obj)
    offset_distance = float(STYLE.get('labels', {}).get('offset_distance', 5.0))
    dx = prot_cent[0] - lig_cent[0]
    dy = prot_cent[1] - lig_cent[1]
    dz = prot_cent[2] - lig_cent[2]
    d = math.sqrt(dx**2 + dy**2 + dz**2)
    if d > 0.01:
        return [prot_cent[0] + dx / d * offset_distance,
                prot_cent[1] + dy / d * offset_distance,
                prot_cent[2] + dz / d * offset_distance]
    return [prot_cent[0] + offset_distance, prot_cent[1], prot_cent[2]]


def _cleanup_temporary_objects():
    """Remove pseudoatoms, distance objects, and selections from interaction scene."""
    for name in cmd.get_names('objects'):
        if any(name.startswith(p) for p in ('cent_', 'lbl_', 'cat_', 'catpi_')):
            cmd.delete(name)
    for name in cmd.get_names('objects'):
        if any(name.startswith(p) for p in ('pi_', 'hb_', 'metal_', 'catpi_', 'sb_', 'contact_')):
            cmd.delete(name)
    # Delete leftover selections
    for sel_name in ('_lig_h_np', '_lig_h_p', '_metal_coord', '_sb_contacts',
                     '_lig_h_np_m', '_lig_h_p_m', '_pocket', '_lig_tmp'):
        try:
            cmd.delete(sel_name)
        except Exception:
            pass


def _angle_separation(a, b):
    dx = abs(float(a[0]) - float(b[0]))
    dy = abs((float(a[1]) - float(b[1]) + 180) % 360 - 180)
    return dx + dy * 0.5


def _select_macro_angles(protein_obj, ligand_obj, output_dir):
    """Choose macro overview angles by scoring low-resolution previews."""
    macro_style = STYLE.get('macro', {})
    fallback = macro_style.get('angles', DEFAULT_STYLE['macro']['angles'])
    if not macro_style.get('auto_angles', True):
        return fallback

    candidates = macro_style.get('candidate_angles', DEFAULT_STYLE['macro']['candidate_angles'])
    if not candidates:
        return fallback

    preview_size = macro_style.get('preview_size', [600, 450])
    try:
        preview_w = int(preview_size[0])
        preview_h = int(preview_size[1])
    except (TypeError, ValueError, IndexError):
        preview_w, preview_h = 600, 450

    preview_dir = os.path.join(output_dir, '_macro_angle_previews')
    os.makedirs(preview_dir, exist_ok=True)

    base_view = cmd.get_view()
    scored = []
    print("  Selecting macro angles from preview renders...")
    for idx, pair in enumerate(candidates, 1):
        try:
            turn_x, turn_y = float(pair[0]), float(pair[1])
        except (TypeError, ValueError, IndexError):
            continue
        cmd.set_view(base_view)
        cmd.zoom(f'{protein_obj} or {ligand_obj}', buffer=3)
        cmd.turn('x', turn_x)
        cmd.turn('y', turn_y)
        preview_path = os.path.join(preview_dir, f'candidate_{idx:02d}.png')
        cmd.ray(preview_w, preview_h)
        cmd.png(preview_path, dpi=72)
        metrics = _score_macro_preview(preview_path)
        if metrics is None:
            print("  Macro auto-angle scoring skipped - Pillow Python not found")
            cmd.set_view(base_view)
            return fallback
        scored.append({
            'turn_x': turn_x,
            'turn_y': turn_y,
            'score': float(metrics.get('score', 9999)),
            'dark_fraction': float(metrics.get('dark_fraction', 1.0)),
            'ligand_fraction': float(metrics.get('ligand_fraction', 0.0)),
        })

    cmd.set_view(base_view)
    if not macro_style.get('keep_angle_previews', False):
        try:
            import shutil
            shutil.rmtree(preview_dir)
        except OSError:
            pass
    if len(scored) < 6:
        return fallback

    selected = []
    for item in sorted(scored, key=lambda row: row['score']):
        if len(selected) >= 6:
            break
        if all(_angle_separation((item['turn_x'], item['turn_y']),
                                 (other['turn_x'], other['turn_y'])) >= 28
               for other in selected):
            selected.append(item)

    if len(selected) < 6:
        for item in sorted(scored, key=lambda row: row['score']):
            if len(selected) >= 6:
                break
            if item not in selected:
                selected.append(item)

    selected.sort(key=lambda row: row['turn_y'])
    angles = []
    for idx, item in enumerate(selected[:6], 1):
        print(f"  macro_{idx}: x={item['turn_x']:.0f}, y={item['turn_y']:.0f}, dark={item['dark_fraction']:.3f}")
        angles.append([item['turn_x'], item['turn_y'], f'macro_{idx}.png'])
    return angles


def setup_macro_scene(protein_obj, ligand_obj, output_dir):
    """Build the macro overview scene."""
    # IRON RULE 2: hide everything FIRST
    cmd.hide('everything')

    # 鈹€鈹€ Protein cartoon 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.show('cartoon', protein_obj)
    cmd.color_deep(_safe_color(STYLE.get('macro', {}).get('protein_cartoon_color'), 'palecyan'), protein_obj, 0)

    # 鈹€鈹€ Global settings 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    _apply_background()
    cmd.set('stick_radius', 0.2)
    cmd.set('sphere_scale', 0.3)
    cmd.set('valence', 1)
    cmd.set('specular', 0.2)
    cmd.set('depth_cue', 0)

    # 鈹€鈹€ Ray trace settings 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    _apply_render_settings()

    # 鈹€鈹€ Metal ions as spheres (auto-detect) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    _apply_macro_render_settings()
    metal_elements = ' or '.join(f'elem {m}' for m in sorted(METAL_RESIDUES))
    cmd.select('_metals', f'{protein_obj} and ({metal_elements})')
    if cmd.count_atoms('_metals') > 0:
        cmd.show('spheres', '_metals')
        cmd.set('sphere_scale', 0.3, '_metals')
        util.cnc('_metals', _self=cmd)
    cmd.delete('_metals')

    # 鈹€鈹€ Ligand 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.show('sticks', ligand_obj)
    util.cba(int(STYLE.get('ligand_cba_index', LIGAND_CBA_INDEX)), ligand_obj, _self=cmd)
    cmd.select('_lig_h_np_m', f'{ligand_obj} and elem H and (neighbor elem C)')
    cmd.hide('everything', '_lig_h_np_m')
    cmd.select('_lig_h_p_m', f'{ligand_obj} and elem H and not (neighbor elem C)')
    cmd.show('sticks', '_lig_h_p_m')
    cmd.color('white', '_lig_h_p_m')

    # 鈹€鈹€ Binding pocket surface (5 脜 around ligand) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.select('_pocket', f'{protein_obj} within {POCKET_RADIUS} of {ligand_obj}')
    cmd.show('surface', '_pocket')
    macro_style = STYLE.get('macro', {})
    cmd.set('surface_color', _safe_color(macro_style.get('surface_color'), 'lightblue'), '_pocket')
    cmd.set('transparency', float(macro_style.get('surface_transparency', 0.4)), '_pocket')

    # 鈹€鈹€ Clean up temp selections 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    cmd.delete('_lig_h_np_m')
    cmd.delete('_lig_h_p_m')
    cmd.delete('_pocket')

    # 鈹€鈹€ Render 6 selected macro angles 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if macro_style.get('paired_to_interaction', True) and INTERACTION_VIEWS:
        print("  Rendering macro views paired to interaction views...")
        paired_zoom_buffer = float(macro_style.get('paired_zoom_buffer', 10.0))
        for idx, item in enumerate(INTERACTION_VIEWS, 1):
            fname = item.get('macro_file', f'macro_{idx}.png')
            cmd.set_view(item['view'])
            cmd.zoom(f'{protein_obj} or {ligand_obj}', buffer=paired_zoom_buffer)
            cmd.ray(RAY_WIDTH, RAY_HEIGHT)
            cmd.png(os.path.join(output_dir, fname), dpi=300)
            paired_name = item.get('interaction_file', f'interaction_{idx}.png')
            print(f"  Saved: {fname} (paired with {paired_name})")
        return

    macro_angles = _select_macro_angles(protein_obj, ligand_obj, output_dir)
    macro_base_view = cmd.get_view()
    for turn_x, turn_y, fname in macro_angles:
        cmd.set_view(macro_base_view)
        cmd.zoom(f'{protein_obj} or {ligand_obj}', buffer=3)
        cmd.turn('x', turn_x)
        cmd.turn('y', turn_y)
        cmd.ray(RAY_WIDTH, RAY_HEIGHT)
        cmd.png(os.path.join(output_dir, fname), dpi=300)
        print(f"  Saved: {fname}")


# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
# Main
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲

def main():
    global STYLE
    args = parse_args()

    # Validate input
    if not os.path.exists(args.input):
        print(f"ERROR: input file not found: {args.input}")
        sys.exit(1)

    try:
        STYLE = load_style_profile(args.style_profile)
    except Exception as exc:
        print(f"ERROR: could not load style profile: {exc}")
        sys.exit(1)

    # Parse interactions
    interactions = parse_interactions(args.interactions)
    interactions = apply_close_contact_policy(
        interactions, args.include_close_contacts)
    if not interactions:
        print(
            "ERROR: no drawable interactions remain. Close contacts require "
            "the explicit --include-close-contacts option."
        )
        sys.exit(1)
    print(f"Interactions: {interactions}")

    # Determine output directory
    if args.output:
        output_dir = args.output
    else:
        input_dir = os.path.dirname(os.path.abspath(args.input))
        output_dir = os.path.join(input_dir, 'pymol_figures')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output: {output_dir}")

    # 鈹€鈹€ Load structure 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    protein_obj, ligand_obj = load_structure(args.input)
    if ligand_obj is None:
        print("ERROR: could not identify ligand object")
        sys.exit(1)

    # 鈹€鈹€ Interaction figure (浣滅敤鍥? 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    print("\n=== Interaction Figure ===")
    setup_interaction_scene(protein_obj, ligand_obj, interactions, output_dir)

    # 鈹€鈹€ Clean up temporary objects from interaction scene 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    print("\n=== Macro Figure ===")
    _cleanup_temporary_objects()
    setup_macro_scene(protein_obj, ligand_obj, output_dir)

    # 鈹€鈹€ Save session 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    session_path = os.path.join(output_dir, 'session.pse')
    cmd.save(session_path)
    print(f"\nSession saved: {session_path}")
    print("Done. All figures generated.")


if __name__ == '__main__':
    main()

