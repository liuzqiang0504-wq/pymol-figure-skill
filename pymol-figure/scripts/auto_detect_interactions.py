"""
Auto-detect protein-ligand interactions from a PDB file.

RDKit is required for automatic interaction detection. If the current Python
cannot import RDKit, the script tries configured RDKit-capable Python commands
such as PYMOL_FIGURE_RDKIT_PYTHON and py -3.12, then exits with an actionable
error if none works.

Usage:
    py -3.12 auto_detect_interactions.py <input.pdb> [--ligand MGP]
       [--chain B] [--output interactions.txt]

Output: interaction spec string compatible with pymol_render.py
    e.g. "B/TYR/1719 pi-pi, B/GLN/1707 hbond, B/LEU/1772 hbond"
"""
import sys, io, os, argparse, math, subprocess, shlex
from collections import defaultdict

try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError:
    Chem = None
    Point3D = None
    RDKIT_AVAILABLE = False
else:
    RDKIT_AVAILABLE = True


def _python_command_from_env(value):
    if not value:
        return None
    try:
        parts = shlex.split(value)
    except ValueError:
        parts = [value]
    return parts or None


def _read_user_env_var(name):
    """Read HKCU user environment values even when Codex has stale process env."""
    if os.name != "nt":
        return None
    try:
        import winreg
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as key:
            value, _ = winreg.QueryValueEx(key, name)
    except OSError:
        return None
    if isinstance(value, str) and value.strip():
        return value.strip()
    return None


def _looks_like_windowsapps_python(cmdline):
    if not cmdline:
        return False
    exe = os.path.normcase(str(cmdline[0]))
    return "\\microsoft\\windowsapps\\" in exe


def _add_unique_command(preferred, delayed, seen, cmdline):
    if not cmdline:
        return
    key = tuple(cmdline)
    if key in seen:
        return
    seen.add(key)
    if _looks_like_windowsapps_python(cmdline):
        delayed.append(cmdline)
    else:
        preferred.append(cmdline)


def _rdkit_python_candidates():
    preferred = []
    delayed = []
    seen = set()

    for value in (
        _read_user_env_var("PYMOL_FIGURE_RDKIT_PYTHON"),
        os.environ.get("PYMOL_FIGURE_RDKIT_PYTHON"),
    ):
        _add_unique_command(preferred, delayed, seen, _python_command_from_env(value))

    if os.name == "nt":
        local_appdata = os.environ.get("LOCALAPPDATA")
        for path in (
            os.path.join(local_appdata or "", "Programs", "Python", "Python312", "python.exe"),
            r"C:\Python312\python.exe",
            r"C:\Program Files\Python312\python.exe",
        ):
            if path and os.path.exists(path):
                _add_unique_command(preferred, delayed, seen, [path])

    for cmdline in (["py", "-3.12"], ["python3.12"]):
        _add_unique_command(preferred, delayed, seen, cmdline)

    return preferred + delayed

def _reexec_with_rdkit_python(argv):
    """Re-run this detector with a Python that has RDKit, if one is available."""
    if os.environ.get("PYMOL_FIGURE_RDKIT_REEXEC") == "1":
        return False
    check_code = "from rdkit import Chem; print('RDKit OK')"
    for cmdline in _rdkit_python_candidates():
        try:
            probe = subprocess.run(
                cmdline + ["-c", check_code],
                capture_output=True,
                text=True,
                timeout=12,
            )
        except (FileNotFoundError, OSError, subprocess.TimeoutExpired):
            continue
        if probe.returncode != 0:
            continue
        env = os.environ.copy()
        env["PYMOL_FIGURE_RDKIT_REEXEC"] = "1"
        args = list(sys.argv[1:] if argv is None else argv)
        print(f"Current Python has no RDKit; rerunning auto-detection with: {' '.join(cmdline)}")
        result = subprocess.run(cmdline + [os.path.abspath(__file__)] + args, env=env)
        sys.exit(result.returncode)
    return False


def _fail_missing_rdkit():
    print("ERROR: RDKit is required for automatic interaction detection.", file=sys.stderr)
    print("", file=sys.stderr)
    print("Tried current Python plus:", file=sys.stderr)
    for cmdline in _rdkit_python_candidates():
        print(f"  {' '.join(cmdline)}", file=sys.stderr)
    print("", file=sys.stderr)
    print("Fix options:", file=sys.stderr)
    print("  1. Run auto-detection with: py -3.12 auto_detect_interactions.py <input.pdb>", file=sys.stderr)
    print("  2. Or set PYMOL_FIGURE_RDKIT_PYTHON to a Python executable that has RDKit.", file=sys.stderr)
    print("  3. Or install RDKit into the Python used by Codex.", file=sys.stderr)
    sys.exit(2)


# 鈹€鈹€ Aromatic ring definitions by residue name 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
PROTEIN_ARO_RINGS = {
    "PHE": ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"],
    "TYR": ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"],
    "TRP": ["CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2"],
    "HIS": ["CG", "ND1", "CE1", "NE2", "CD2"],
}

# Charged sidechain atoms
CATIONS = {
    "ARG": ["NH1", "NH2", "NE"],
    "LYS": ["NZ"],
    "HIS": ["ND1", "NE2"],
}

ANIONS = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

# Hydrophobic contact atoms
HYDROPHOBIC_ATOMS = {
    "ALA": ["CB"],
    "VAL": ["CG1", "CG2"],
    "LEU": ["CD1", "CD2"],
    "ILE": ["CG2", "CD1"],
    "PRO": ["CB", "CG", "CD"],
    "MET": ["CE"],
}

# Common ligand names (HETATM residues that are NOT water/ions/buffer)
LIGAND_LIKE = {
    # Nucleotides & analogs
    "MGP", "GTP", "GDP", "GMP", "ATP", "ADP", "AMP", "CTP", "UTP", "TTP",
    "M7G", "MGT", "7MG", "GOL", "NMN",
    # Cofactors
    "NAD", "NAI", "NAP", "FAD", "FMN", "COA", "HEM", "HEA", "SAM", "SAH",
    "PLP", "TPP", "THF", "BH4",
    # Sugars
    "NAG", "BMA", "MAN", "GLC", "GAL", "FUC", "XYS", "SIA",
    # Lipids
    "PLM", "LDA", "OLA", "STE", "MYR",
    # Common drug-like fragments
    "BEN", "ACT", "EDO", "PEG", "SO4", "PO4",
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Auto-detect protein-ligand interactions from PDB")
    p.add_argument("input", help="PDB file path")
    p.add_argument("--ligand", "-l", default=None,
                   help="Ligand residue name (e.g. MGP). Auto-detected if omitted.")
    p.add_argument("--chain", "-c", default=None,
                   help="Chain ID for protein (auto-detected if omitted)")
    p.add_argument("--hbond-cutoff", type=float, default=3.5,
                   help="H-bond distance cutoff in A (default: 3.5)")
    p.add_argument("--pipi-cutoff", type=float, default=5.5,
                   help="pi-pi centroid distance cutoff in A (default: 5.5)")
    p.add_argument("--salt-cutoff", type=float, default=4.0,
                   help="Salt bridge distance cutoff in A (default: 4.0)")
    return p.parse_args(argv)


def load_pdb(filepath):
    """Load PDB with RDKit, classify atoms."""
    mol = Chem.MolFromPDBFile(filepath, removeHs=False, proximityBonding=True)
    if mol is None:
        print("ERROR: Could not read PDB", file=sys.stderr)
        sys.exit(1)
    return mol


def classify_atoms(mol):
    """Separate atoms into protein, ligand, water, ion categories."""
    groups = defaultdict(list)
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        rname = info.GetResidueName().strip()
        if rname == "HOH":
            groups["water"].append(atom.GetIdx())
        elif rname in ("NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE"):
            groups["ion"].append(atom.GetIdx())
        else:
            groups[f"res_{rname}"].append(atom.GetIdx())
    return groups


def detect_ligand(groups, specified=None):
    """Identify the ligand residue. Uses user-specified name or auto-detects."""
    if specified:
        key = f"res_{specified}"
        if key in groups:
            return specified, groups[key]
        else:
            print(f"WARNING: specified ligand '{specified}' not found, auto-detecting...",
                  file=sys.stderr)

    # Auto-detect: find HETATM residues that look like ligands
    candidates = []
    for key, atoms in groups.items():
        if key.startswith("res_"):
            rname = key[4:]
            if rname in LIGAND_LIKE:
                candidates.append((rname, atoms))
    if candidates:
        # Pick the largest one
        candidates.sort(key=lambda x: len(x[1]), reverse=True)
        return candidates[0][0], candidates[0][1]

    return None, []


def detect_protein_chain(mol, groups, ligand_atoms):
    """Auto-detect protein chain ID by finding standard amino acids."""
    STANDARD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                   "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                   "TYR", "VAL"}
    chains = defaultdict(int)
    for key, atoms in groups.items():
        if not key.startswith("res_"):
            continue
        rname = key[4:]
        if rname not in STANDARD_AA:
            continue
        for idx in atoms:
            info = mol.GetAtomWithIdx(idx).GetPDBResidueInfo()
            if info:
                chains[info.GetChainId().strip()] += 1
    if chains:
        return max(chains, key=chains.get)
    return ""


def get_res_groups(mol, atom_indices):
    """Group atom indices by residue label."""
    res_map = defaultdict(list)
    for idx in atom_indices:
        info = mol.GetAtomWithIdx(idx).GetPDBResidueInfo()
        key = (info.GetResidueName().strip(), info.GetResidueNumber(),
               info.GetChainId().strip())
        res_map[key].append(idx)
    return res_map


def get_atom_map(mol, atom_list):
    """Create {atom_name: atom_idx} for a set of atom indices."""
    return {mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetName().strip(): i
            for i in atom_list}


def atom_pos(mol, idx):
    p = mol.GetConformer().GetAtomPosition(idx)
    return (p.x, p.y, p.z)


def distance(p1, p2):
    return math.sqrt(sum((a-b)**2 for a, b in zip(p1, p2)))


def centroid(positions):
    n = len(positions)
    return tuple(sum(p[i] for p in positions) / n for i in range(3))


def find_aromatic_rings(res_groups, mol):
    """Find aromatic rings in protein residues using name-based templates."""
    rings = []
    for (resn, resi, chain), atom_list in res_groups.items():
        if resn not in PROTEIN_ARO_RINGS:
            continue
        atom_map = get_atom_map(mol, atom_list)
        template = PROTEIN_ARO_RINGS[resn]
        positions = []
        for aname in template:
            if aname in atom_map:
                positions.append(atom_pos(mol, atom_map[aname]))
        if len(positions) >= 4:
            rings.append((resn, resi, chain, positions))
    return rings


def find_ligand_rings(mol, ligand_atoms, lig_resn, lig_resi, lig_chain):
    """Find aromatic rings in ligand using geometric ring detection.
    Uses 5/6-membered ring detection from PDB CONECT/atom names."""
    rings = []
    ri = mol.GetRingInfo()
    lig_set = set(ligand_atoms)
    for ring in ri.AtomRings():
        if len(ring) not in (5, 6):
            continue
        if not set(ring).issubset(lig_set):
            continue
        is_aro = all(mol.GetAtomWithIdx(i).GetAtomicNum() in (6, 7)
                     for i in ring)
        if is_aro:
            positions = [atom_pos(mol, i) for i in ring]
            rings.append((lig_resn, lig_resi, lig_chain, positions))
    return rings


def detect_hbonds(mol, protein_atoms, ligand_atoms):
    """Detect H-bonds: protein donor/acceptor <-> ligand donor/acceptor.
    Uses explicit H positions from PDB for angle validation."""
    hbonds = []
    conf = mol.GetConformer()

    def _pos(idx):
        p = conf.GetAtomPosition(idx)
        return Point3D(p.x, p.y, p.z)

    for la_idx in ligand_atoms:
        la = mol.GetAtomWithIdx(la_idx)
        linfo = la.GetPDBResidueInfo()
        lres = f"{linfo.GetResidueName().strip()}{linfo.GetResidueNumber()}"

        # Ligand as donor (has attached H)
        if la.GetAtomicNum() in (7, 8):
            for h in la.GetNeighbors():
                if h.GetAtomicNum() != 1:
                    continue
                for pa_idx in protein_atoms:
                    pa = mol.GetAtomWithIdx(pa_idx)
                    if pa.GetAtomicNum() not in (7, 8):
                        continue
                    h_pos = _pos(h.GetIdx())
                    a_pos = _pos(pa_idx)
                    d = h_pos.Distance(a_pos)
                    if d > 2.8:
                        continue
                    d_pos = _pos(la_idx)
                    v1 = Point3D(d_pos.x-h_pos.x, d_pos.y-h_pos.y, d_pos.z-h_pos.z)
                    v2 = Point3D(a_pos.x-h_pos.x, a_pos.y-h_pos.y, a_pos.z-h_pos.z)
                    v1m = math.sqrt(v1.x**2+v1.y**2+v1.z**2)
                    v2m = math.sqrt(v2.x**2+v2.y**2+v2.z**2)
                    if v1m > 1e-6 and v2m > 1e-6:
                        cos_ang = (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/(v1m*v2m)
                        ang = math.degrees(math.acos(max(-1, min(1, cos_ang))))
                        if ang >= 110:
                            pinfo = pa.GetPDBResidueInfo()
                            pres = f"{pinfo.GetResidueName().strip()}{pinfo.GetResidueNumber()}"
                            chain = pinfo.GetChainId().strip()
                            hbonds.append((pres, pinfo.GetResidueName().strip(),
                                           pinfo.GetResidueNumber(), chain, d, ang))

        # Ligand as acceptor
        if la.GetAtomicNum() in (7, 8):
            for pa_idx in protein_atoms:
                pa = mol.GetAtomWithIdx(pa_idx)
                if pa.GetAtomicNum() not in (7, 8):
                    continue
                for h in pa.GetNeighbors():
                    if h.GetAtomicNum() != 1:
                        continue
                    h_pos = _pos(h.GetIdx())
                    l_pos = _pos(la_idx)
                    d = h_pos.Distance(l_pos)
                    if d > 2.8:
                        continue
                    d_pos = _pos(pa_idx)
                    v1 = Point3D(d_pos.x-h_pos.x, d_pos.y-h_pos.y, d_pos.z-h_pos.z)
                    v2 = Point3D(l_pos.x-h_pos.x, l_pos.y-h_pos.y, l_pos.z-h_pos.z)
                    v1m = math.sqrt(v1.x**2+v1.y**2+v1.z**2)
                    v2m = math.sqrt(v2.x**2+v2.y**2+v2.z**2)
                    if v1m > 1e-6 and v2m > 1e-6:
                        cos_ang = (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/(v1m*v2m)
                        ang = math.degrees(math.acos(max(-1, min(1, cos_ang))))
                        if ang >= 110:
                            pinfo = pa.GetPDBResidueInfo()
                            pres = f"{pinfo.GetResidueName().strip()}{pinfo.GetResidueNumber()}"
                            chain = pinfo.GetChainId().strip()
                            hbonds.append((pres, pinfo.GetResidueName().strip(),
                                           pinfo.GetResidueNumber(), chain, d, ang))
    return hbonds


def detect_pipi(protein_rings, ligand_rings, cutoff=5.5):
    """Detect pi-pi stacking by ring centroid distances."""
    pipi = []
    for presn, presi, pchain, ppos in protein_rings:
        pc = centroid(ppos)
        for lresn, lresi, lchain, lpos in ligand_rings:
            lc = centroid(lpos)
            d = distance(pc, lc)
            if d <= cutoff:
                pipi.append((f"{presn}{presi}", presn, presi, pchain, d))
    return pipi


def detect_salt_bridges(res_groups, mol, ligand_atoms, cutoff=4.0):
    """Detect salt bridges: protein charged <-> ligand charged."""
    bridges = []
    for (resn, resi, chain), atom_list in res_groups.items():
        atom_map = get_atom_map(mol, atom_list)
        label = f"{resn}{resi}"

        # Protein cations -> ligand anions? (ligand phosphates, carboxylates)
        if resn in CATIONS:
            for aname in CATIONS[resn]:
                if aname not in atom_map:
                    continue
                cat_pos = atom_pos(mol, atom_map[aname])
                for la_idx in ligand_atoms:
                    la = mol.GetAtomWithIdx(la_idx)
                    lname = la.GetPDBResidueInfo().GetName().strip()
                    # Check for phosphate oxygens or carboxylate oxygens
                    if la.GetAtomicNum() != 8:
                        continue
                    d = distance(cat_pos, atom_pos(mol, la_idx))
                    if d <= cutoff:
                        bridges.append((label, resn, resi, chain, aname,
                                        lname, d, "prot+ -> lig-"))
                        break

        # Protein anions -> ligand cations
        if resn in ANIONS:
            for aname in ANIONS[resn]:
                if aname not in atom_map:
                    continue
                an_pos = atom_pos(mol, atom_map[aname])
                for la_idx in ligand_atoms:
                    la = mol.GetAtomWithIdx(la_idx)
                    # Check for ligand N+ (e.g., 7-methylguanosine N7)
                    if la.GetAtomicNum() != 7:
                        continue
                    lname = la.GetPDBResidueInfo().GetName().strip()
                    d = distance(an_pos, atom_pos(mol, la_idx))
                    if d <= cutoff:
                        bridges.append((label, resn, resi, chain, aname,
                                        lname, d, "prot- -> lig+"))
                        break
    return bridges


def format_interaction_spec(hbonds, pipi_results, salt_bridges, chain):
    """Format detected interactions into the spec string for pymol_render.py."""
    # Track unique residues by interaction type
    hbond_res = set()
    pipi_res = set()
    salt_res = set()

    for pres, resn, resi, pchain, d, ang in hbonds:
        hbond_res.add((resn, resi, pchain))

    for pres, resn, resi, pchain, d in pipi_results:
        pipi_res.add((resn, resi, pchain))

    for label, resn, resi, bchain, aname, lname, d, desc in salt_bridges:
        salt_res.add((resn, resi, bchain))

    # Build spec: pi-pi first (most important), then hbond, then salt bridges
    specs = []

    # pi-pi: use chain-aware format
    for resn, resi, pchain in sorted(pipi_res, key=lambda x: x[1]):
        specs.append(f"{pchain}/{resn}/{resi} pi-pi")

    # hbond: exclude residues already listed as pi-pi
    for resn, resi, pchain in sorted(hbond_res, key=lambda x: x[1]):
        if (resn, resi, pchain) not in pipi_res:
            specs.append(f"{pchain}/{resn}/{resi} hbond")

    # salt bridges
    for resn, resi, bchain in sorted(salt_res, key=lambda x: x[1]):
        if (resn, resi, bchain) not in pipi_res and (resn, resi, bchain) not in hbond_res:
            specs.append(f"{bchain}/{resn}/{resi} salt-bridge")

    return ", ".join(specs)


def main(argv=None):
    if not RDKIT_AVAILABLE:
        _reexec_with_rdkit_python(argv)
        _fail_missing_rdkit()

    args = parse_args(argv)
    mol = load_pdb(args.input)

    # Classify atoms
    groups = classify_atoms(mol)
    print(f"Found {sum(len(v) for v in groups.values())} atoms in "
          f"{len([k for k in groups if k.startswith('res_')])} residue types")

    # Detect ligand
    lig_name, lig_atoms = detect_ligand(groups, args.ligand)
    if not lig_atoms:
        print("ERROR: no ligand found. Specify with --ligand RESNAME", file=sys.stderr)
        sys.exit(1)
    print(f"Ligand: {lig_name} ({len(lig_atoms)} atoms)")

    # Get ligand residue info
    lig_info = mol.GetAtomWithIdx(lig_atoms[0]).GetPDBResidueInfo()
    lig_resi = lig_info.GetResidueNumber()
    lig_chain = lig_info.GetChainId().strip()

    # Detect protein chain
    chain = args.chain or detect_protein_chain(mol, groups, lig_atoms)
    if not chain:
        print("ERROR: could not determine protein chain. Use --chain", file=sys.stderr)
        sys.exit(1)
    print(f"Protein chain: {chain}")

    # Collect protein atoms (exclude ligand, water, ions)
    protein_atoms = []
    exclude_keys = {f"res_{lig_name}", "water", "ion"}
    for key, atoms in groups.items():
        if key not in exclude_keys and key.startswith("res_"):
            protein_atoms.extend(atoms)

    # Also exclude atoms from the same chain as ligand if they're HETATM
    lig_chain_atoms = set()
    for idx in lig_atoms:
        info = mol.GetAtomWithIdx(idx).GetPDBResidueInfo()
        if info:
            lig_chain_atoms.add(info.GetChainId().strip())
    protein_atoms = [i for i in protein_atoms
                     if mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetChainId().strip()
                     not in lig_chain_atoms
                     or mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName().strip()
                     in ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                         "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                         "TYR", "VAL")]

    print(f"Protein atoms: {len(protein_atoms)}")

    # Group protein by residue
    res_groups = get_res_groups(mol, protein_atoms)

    # Detect interactions
    print("\n--- Detecting interactions ---")

    # H-bonds
    hbonds = detect_hbonds(mol, protein_atoms, lig_atoms)
    print(f"\nH-bonds found: {len(hbonds)}")
    seen = set()
    for pres, resn, resi, pchain, d, ang in sorted(hbonds, key=lambda x: x[4]):
        key = (pres, pchain)
        if key not in seen:
            print(f"  {pchain}/{resn}/{resi}  d={d:.1f}A  angle={ang:.0f}deg")
            seen.add(key)

    # pi-pi
    protein_rings = find_aromatic_rings(res_groups, mol)
    print(f"\nProtein aromatic rings: {len(protein_rings)}")
    ligand_rings = find_ligand_rings(mol, lig_atoms, lig_name, lig_resi, lig_chain)
    print(f"Ligand aromatic rings: {len(ligand_rings)}")
    pipi = detect_pipi(protein_rings, ligand_rings, args.pipi_cutoff)
    print(f"pi-pi interactions: {len(pipi)}")
    for pres, resn, resi, pchain, d in sorted(pipi, key=lambda x: x[4]):
        print(f"  {pchain}/{resn}/{resi}  centroid dist={d:.1f}A")

    # Salt bridges
    salt = detect_salt_bridges(res_groups, mol, lig_atoms, args.salt_cutoff)
    print(f"\nSalt bridges: {len(salt)}")
    for s in salt:
        print(f"  {s[4]}:{s[0]} <-> {s[5]}  d={s[6]:.1f}A  ({s[7]})")

    # Format output
    spec = format_interaction_spec(hbonds, pipi, salt, chain)
    print(f"\n=== Interaction spec ===")
    print(spec)
    print()

    return spec


if __name__ == "__main__":
    main()

