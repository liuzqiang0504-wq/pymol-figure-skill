"""
Ring detection for protein residues and ligand molecules.

Pure Python — no PyMOL dependency. Two backends:
  - RDKit (preferred): if rdkit is importable
  - Geometric (fallback): covalent-radius adjacency + DFS + planarity filter

Public API:
  - AROMATIC_RING_ATOMS: dict[str, list[str]] — hardcoded protein ring lookup
  - detect_rings(atom_data) -> list[list[int]] — returns ring atom index sets
  - centroid(positions) -> list[float] — ring centroid from list of [x,y,z] lists
"""

import math
import itertools

# ── Protein aromatic ring atom lookup ──────────────────────────────────────
# Do NOT attempt to detect rings in protein residues. Use this table.

AROMATIC_RING_ATOMS: dict[str, list[str]] = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2'],  # 6-membered ring
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
}

# TRP also has a 5-membered pyrrole ring (CD2, CE2, NE1, CE3, CZ2)
# Default to the 6-membered ring unless specified otherwise.

TRP_PYRROLE_ATOMS = ['CD2', 'CE2', 'NE1', 'CE3', 'CZ2']

# ── Covalent radii (Angstroms) ─────────────────────────────────────────────
# From Alvarez (2008) — used for adjacency detection in geometric fallback.

COVALENT_RADII: dict[str, float] = {
    'H': 0.31,  'B': 0.84,  'C': 0.76,  'N': 0.71,  'O': 0.66,
    'F': 0.57,  'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02,
    'Br': 1.20, 'I': 1.39,
}

BOND_TOLERANCE = 1.25  # multiplier on sum of covalent radii


# ── Distance helpers ───────────────────────────────────────────────────────

def _distance(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def centroid(positions: list[list[float]]) -> list[float]:
    """Centroid of a set of 3D positions."""
    n = len(positions)
    if n == 0:
        return [0.0, 0.0, 0.0]
    return [sum(p[i] for p in positions) / n for i in range(3)]


def _smallest_eigenpair(cov: list[list[float]]) -> tuple[list[float], float]:
    """Return (eigenvector, eigenvalue) for the smallest eigenvalue of a 3x3
    symmetric matrix. Uses power iteration with deflation."""
    # Power iteration for dominant eigenpair
    v = [1.0, 2.0, 3.0]
    for _ in range(30):
        mv = [
            cov[0][0] * v[0] + cov[0][1] * v[1] + cov[0][2] * v[2],
            cov[1][0] * v[0] + cov[1][1] * v[1] + cov[1][2] * v[2],
            cov[2][0] * v[0] + cov[2][1] * v[1] + cov[2][2] * v[2],
        ]
        norm = math.sqrt(mv[0]**2 + mv[1]**2 + mv[2]**2)
        if norm < 1e-12:
            break
        v = [mv[0] / norm, mv[1] / norm, mv[2] / norm]
    # Rayleigh quotient for eigenvalue
    mv = [
        cov[0][0] * v[0] + cov[0][1] * v[1] + cov[0][2] * v[2],
        cov[1][0] * v[0] + cov[1][1] * v[1] + cov[1][2] * v[2],
        cov[2][0] * v[0] + cov[2][1] * v[1] + cov[2][2] * v[2],
    ]
    lambda1 = v[0] * mv[0] + v[1] * mv[1] + v[2] * mv[2]
    # Deflate: C' = C - λ1 * v1 * v1^T
    cov2 = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            cov2[i][j] = cov[i][j] - lambda1 * v[i] * v[j]
    # Power iteration on deflated matrix for second eigenpair
    v2 = [3.0, 2.0, 1.0]
    for _ in range(30):
        mv2 = [
            cov2[0][0] * v2[0] + cov2[0][1] * v2[1] + cov2[0][2] * v2[2],
            cov2[1][0] * v2[0] + cov2[1][1] * v2[1] + cov2[1][2] * v2[2],
            cov2[2][0] * v2[0] + cov2[2][1] * v2[1] + cov2[2][2] * v2[2],
        ]
        norm2 = math.sqrt(mv2[0]**2 + mv2[1]**2 + mv2[2]**2)
        if norm2 < 1e-12:
            break
        v2 = [mv2[0] / norm2, mv2[1] / norm2, mv2[2] / norm2]
    mv2 = [
        cov2[0][0] * v2[0] + cov2[0][1] * v2[1] + cov2[0][2] * v2[2],
        cov2[1][0] * v2[0] + cov2[1][1] * v2[1] + cov2[1][2] * v2[2],
        cov2[2][0] * v2[0] + cov2[2][1] * v2[1] + cov2[2][2] * v2[2],
    ]
    lambda2 = v2[0] * mv2[0] + v2[1] * mv2[1] + v2[2] * mv2[2]
    # Third eigenvector via cross product
    v3 = [
        v[1] * v2[2] - v[2] * v2[1],
        v[2] * v2[0] - v[0] * v2[2],
        v[0] * v2[1] - v[1] * v2[0],
    ]
    norm3 = math.sqrt(v3[0]**2 + v3[1]**2 + v3[2]**2)
    if norm3 < 1e-12:
        return (v3, 0.0)
    v3 = [v3[0] / norm3, v3[1] / norm3, v3[2] / norm3]
    mv3 = [
        cov[0][0] * v3[0] + cov[0][1] * v3[1] + cov[0][2] * v3[2],
        cov[1][0] * v3[0] + cov[1][1] * v3[1] + cov[1][2] * v3[2],
        cov[2][0] * v3[0] + cov[2][1] * v3[1] + cov[2][2] * v3[2],
    ]
    lambda3 = v3[0] * mv3[0] + v3[1] * mv3[1] + v3[2] * mv3[2]
    # Return the eigenpair with smallest eigenvalue
    pairs = [(v, lambda1), (v2, lambda2), (v3, lambda3)]
    pairs.sort(key=lambda x: x[1])
    return pairs[0]


def _plane_rmsd(positions: list[list[float]]) -> float:
    """RMS deviation from the best-fit plane through the given positions."""
    if len(positions) < 3:
        return float('inf')
    c = centroid(positions)
    # Build covariance matrix
    xx = yy = zz = xy = xz = yz = 0.0
    for p in positions:
        dx, dy, dz = p[0] - c[0], p[1] - c[1], p[2] - c[2]
        xx += dx * dx
        yy += dy * dy
        zz += dz * dz
        xy += dx * dy
        xz += dx * dz
        yz += dy * dz
    cov = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
    normal, lambda_min = _smallest_eigenpair(cov)
    # RMSD = sqrt(λ_min / n)
    variance = max(0.0, lambda_min)
    return math.sqrt(variance / len(positions))


# ── RDKit backend ──────────────────────────────────────────────────────────

def _detect_rings_rdkit(atom_data: list[dict]) -> list[list[int]]:
    """Use RDKit to detect all rings (any size, aromatic or not)."""
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("RDKit not available")
    # Sort atoms by index if present
    atoms_sorted = sorted(enumerate(atom_data),
                          key=lambda x: x[1].get('idx', x[0]))
    mol = Chem.EditableMol(Chem.Mol())
    for _, a in atoms_sorted:
        elem = a.get('elem', 'C')
        atom = Chem.Atom(elem)
        mol.AddAtom(atom)
    # Add bonds based on distance
    coords = []
    for _, a in atoms_sorted:
        coords.append([a.get('x', 0.0), a.get('y', 0.0), a.get('z', 0.0)])
    n = len(coords)
    for i in range(n):
        for j in range(i + 1, n):
            r1 = COVALENT_RADII.get(atoms_sorted[i][1].get('elem', 'C'), 0.76)
            r2 = COVALENT_RADII.get(atoms_sorted[j][1].get('elem', 'C'), 0.76)
            cutoff = (r1 + r2) * BOND_TOLERANCE
            if _distance(coords[i], coords[j]) <= cutoff:
                mol.AddBond(i, j, Chem.BondType.SINGLE)
    mol = mol.GetMol()
    Chem.SanitizeMol(mol, catchErrors=True)
    ring_info = mol.GetRingInfo()
    rings = [list(r) for r in ring_info.AtomRings()]
    # Filter: only rings up to size 7, and check planarity
    result = []
    for ring in rings:
        if len(ring) > 7:
            continue
        ring_positions = [coords[idx] for idx in ring]
        if _plane_rmsd(ring_positions) < 0.1:
            result.append(ring)
    return result


# ── Geometric backend ──────────────────────────────────────────────────────

def _build_adjacency(atom_data: list[dict]) -> list[list[int]]:
    """Build adjacency list from covalent-radii distance checks."""
    n = len(atom_data)
    adj = [[] for _ in range(n)]
    elements = [a.get('elem', 'C') for a in atom_data]
    coords = [[a.get('x', 0.0), a.get('y', 0.0), a.get('z', 0.0)]
              for a in atom_data]
    for i in range(n):
        ri = COVALENT_RADII.get(elements[i], 0.76)
        for j in range(i + 1, n):
            rj = COVALENT_RADII.get(elements[j], 0.76)
            cutoff = (ri + rj) * BOND_TOLERANCE
            if _distance(coords[i], coords[j]) <= cutoff:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _dfs_find_cycles(adj: list[list[int]], max_size: int = 7
                     ) -> list[list[int]]:
    """DFS-based cycle enumeration up to max_size."""
    n = len(adj)
    cycles: list[list[int]] = []
    visited = [False] * n

    def _dfs(start: int, current: int, parent: int, depth: int, path: list[int]):
        if depth > max_size:
            return
        visited[current] = True
        for nb in adj[current]:
            if nb == parent:
                continue
            if nb == start and depth >= 2:
                # Found a cycle
                cycles.append(list(path))
                continue
            if visited[nb]:
                continue
            path.append(nb)
            _dfs(start, nb, current, depth + 1, path)
            path.pop()
        visited[current] = False

    for start in range(n):
        visited[start] = True
        _dfs(start, start, -1, 1, [start])
        visited[start] = False

    # Deduplicate: a cycle of length k appears k times (once per starting node)
    # and in both directions. Normalize: rotate so smallest index is first.
    unique: dict[tuple, list[int]] = {}
    for cycle in cycles:
        if len(cycle) < 3:
            continue
        # Rotate so min element is first
        m = min(cycle)
        idx = cycle.index(m)
        rotated = cycle[idx:] + cycle[:idx]
        # Choose direction: compare second element
        if len(rotated) > 1:
            if rotated[-1] < rotated[1]:
                rotated = [rotated[0]] + list(reversed(rotated[1:]))
        key = tuple(rotated)
        unique[key] = rotated
    return list(unique.values())


def _detect_rings_geometric(atom_data: list[dict]) -> list[list[int]]:
    """Geometric fallback: adjacency + DFS + planarity filter."""
    coords = [[a.get('x', 0.0), a.get('y', 0.0), a.get('z', 0.0)]
              for a in atom_data]
    adj = _build_adjacency(atom_data)
    cycles = _dfs_find_cycles(adj, max_size=7)
    # Filter by planarity
    result = []
    for cycle in cycles:
        ring_positions = [coords[idx] for idx in cycle]
        if _plane_rmsd(ring_positions) < 0.1:
            result.append(cycle)
    return result


# ── Public API ─────────────────────────────────────────────────────────────

def detect_rings(atom_data: list[dict]) -> list[list[int]]:
    """
    Detect aromatic rings in a ligand or small molecule.

    Parameters
    ----------
    atom_data : list[dict]
        Each dict must have keys: 'elem', 'x', 'y', 'z'.
        Optional: 'idx' for explicit indexing.

    Returns
    -------
    list[list[int]]
        Each inner list is a set of atom indices forming a ring.
    """
    # Try RDKit first
    try:
        rings = _detect_rings_rdkit(atom_data)
        if rings:
            return rings
    except ImportError:
        pass
    # Fallback to geometric
    return _detect_rings_geometric(atom_data)


# ── Self-test ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    # Benzene ring test: 6 carbons in a plane, 1.4 Å apart
    import math
    r = 1.4  # C-C bond
    benzene = []
    for i in range(6):
        angle = 2 * math.pi * i / 6
        benzene.append({
            'elem': 'C',
            'x': r * math.cos(angle),
            'y': r * math.sin(angle),
            'z': 0.0,
        })
    rings = detect_rings(benzene)
    print(f"Benzene test: found {len(rings)} rings")
    for ring in rings:
        print(f"  Ring size {len(ring)}: {ring}")
    assert len(rings) >= 1, "Should find at least one ring"
    assert 6 in [len(r) for r in rings], "Should find a 6-membered ring"
    print("OK: benzene ring detected")

    # Protein lookup test
    assert 'PHE' in AROMATIC_RING_ATOMS
    assert AROMATIC_RING_ATOMS['PHE'] == ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    print("OK: AROMATIC_RING_ATOMS lookup correct")

    print("All tests passed.")
