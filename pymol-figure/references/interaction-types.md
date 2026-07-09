# Interaction Types

Detection rules and PyMOL selection algebra for non-covalent interactions.

## Protein Aromatic Ring Atom Lookup

Do NOT attempt to detect rings in protein residues. Use this hardcoded table:

| Residue | Ring name | Atom names |
|---------|-----------|------------|
| PHE | phenyl | CG, CD1, CD2, CE1, CE2, CZ |
| TYR | phenol | CG, CD1, CD2, CE1, CE2, CZ |
| TRP | indole (6-ring) | CG, CD1, CD2, CE2, CE3, CZ2 |
| TRP | pyrrole (5-ring) | CD2, CE2, NE1, CE3, CZ2 |
| HIS | imidazole | CG, ND1, CD2, CE1, NE2 |

For TRP pi-stacking interactions, default to the 6-membered ring unless specified otherwise.

## π-π Stacking Detection

### Geometry criteria
- Ring centroid distance: ≤ 4.5 Å
- Ring plane offset: ≤ 2.0 Å
- Inter-planar angle (dihedral): ≤ 30° (parallel) or 60-90° (T-shaped)

### PyMOL centroid calculation
```python
stored.tmp = []
cmd.iterate_state(1, f'{obj} and chain {chain} and resi {resi} and name {atom_names}',
                  'stored.tmp.append((x, y, z))')
n = len(stored.tmp)
centroid = [sum(c[i] for c in stored.tmp) / n for i in range(3)]
```

### Dashed line setup
```python
cmd.distance('pi_interaction', 'centroid_a', 'centroid_b')
cmd.set('dash_color', 'yellow', 'pi_interaction')
cmd.set('dash_gap', 0.3, 'pi_interaction')
cmd.set('dash_width', 2.5, 'pi_interaction')
cmd.hide('labels', 'pi_interaction')
```

## Hydrogen Bond Detection

### Geometry criteria
- Donor-Acceptor distance: ≤ 3.5 Å
- D-H...A angle: ≥ 120°

### Common H-bond donors/acceptors in proteins

| Role | Residues | Atoms |
|------|----------|-------|
| Donor | ARG | NH1, NH2 |
| Donor | LYS | NZ |
| Donor | ASN | ND2 |
| Donor | GLN | NE2 |
| Donor | SER | OG |
| Donor | THR | OG1 |
| Donor | TYR | OH |
| Donor | HIS | ND1, NE2 |
| Donor | TRP | NE1 |
| Donor | Backbone | N |
| Acceptor | ASP | OD1, OD2 |
| Acceptor | GLU | OE1, OE2 |
| Acceptor | ASN | OD1 |
| Acceptor | GLN | OE1 |
| Acceptor | Backbone | O |

### Dashed line setup
```python
cmd.distance('hbond_name', 'donor_sel', 'acceptor_sel')
cmd.set('dash_color', 'green', 'hbond_name')
cmd.set('dash_gap', 0.3, 'hbond_name')
cmd.set('dash_width', 2.5, 'hbond_name')
cmd.hide('labels', 'hbond_name')
```

## Ligand Ring Detection

Run `scripts/ring_detection.py` against the ligand atom data to identify aromatic rings.

The module provides:
- `detect_rings(atom_data: list[dict]) -> list[list[int]]` — returns list of ring atom index sets
- Uses RDKit backend if available, falls back to geometric graph-based detection
- Discovers 5-, 6-, and 7-membered planar rings

## Metal Coordination Detection

### Metal residues in proteins

Common metal ions found in PDB structures:

| Metal | Residue Name | Typical coordination geometry |
|-------|-------------|-------------------------------|
| Zinc | ZN | Tetrahedral (4), distance 2.0-2.5 A |
| Magnesium | MG | Octahedral (6), distance 2.0-2.5 A |
| Calcium | CA | Octahedral/pentagonal bipyramidal, distance 2.3-2.8 A |
| Manganese | MN | Octahedral, distance 2.1-2.6 A |
| Iron | FE | Octahedral, distance 2.0-2.6 A |
| Cobalt | CO | Octahedral, distance 2.0-2.5 A |
| Nickel | NI | Square planar/octahedral, distance 2.0-2.5 A |
| Copper | CU | Square planar/trigonal bipyramidal, distance 2.0-2.5 A |
| Cadmium | CD | Tetrahedral, distance 2.3-2.8 A |
| Sodium | NA | Variable, distance 2.3-3.0 A |
| Potassium | K | Variable, distance 2.7-3.3 A |

### Metal-coordinating atoms (ligand side)

| Element | Typical context |
|---------|----------------|
| N | Amine, imine, heterocycle (imidazole, pyridine) |
| O | Hydroxyl, carboxylate, carbonyl, water |
| S | Thiol, thioether, thiophene |

### Geometry criteria
- Metal-donor distance: <= 3.0 A (general cutoff for most metals)
- Donor atom must be N, O, or S

### Sphere display
```python
cmd.show('spheres', metal_selection)
cmd.set('sphere_scale', 0.3, metal_selection)
util.cnc(metal_selection, _self=cmd)  # color by chain
```

### Dashed line setup
```python
cmd.distance('metal_name', 'metal_sel', 'donor_sel')
cmd.set('dash_color', 'purple', 'metal_name')
cmd.set('dash_gap', 0.3, 'metal_name')
cmd.set('dash_width', 2.5, 'metal_name')
cmd.hide('labels', 'metal_name')
```

## Cation-π Interaction Detection

### Cation sources (protein side)
- Metal ions: ZN, MG, CA, MN, FE, etc.
- Protonated basic residues: ARG (NH1, NH2, NE), LYS (NZ), HIS (ND1, NE2)

### Cation sources (ligand side)
- Protonated amines (NH3+, NH2+, NH+)
- Quaternary ammonium (NR4+)
- Guanidinium groups
- Metal ions in metallodrugs

### π systems (protein side)
Use hardcoded AROMATIC_RING_ATOMS table (PHE, TYR, TRP, HIS).

### Geometry criteria
- Cation-ring centroid distance: <= 6.0 A
- Preferred: cation directly above ring plane (offset <= 2.0 A)

### Dashed line setup
```python
cmd.distance('catpi_name', 'cation_sel', 'centroid_sel')
cmd.set('dash_color', 'purple', 'catpi_name')
cmd.set('dash_gap', 0.3, 'catpi_name')
cmd.set('dash_width', 2.5, 'catpi_name')
cmd.hide('labels', 'catpi_name')
```

## Salt Bridge Detection

### Acidic residues (negative charge)
| Residue | Charged atoms |
|---------|---------------|
| ASP | OD1, OD2 |
| GLU | OE1, OE2 |

### Basic residues (positive charge)
| Residue | Charged atoms |
|---------|---------------|
| ARG | NH1, NH2, NE |
| LYS | NZ |
| HIS | ND1, NE2 |

### Geometry criteria
- Distance between charged heavy atoms: <= 4.0 A

### Dashed line setup
```python
cmd.distance('sb_name', 'acidic_sel', 'basic_sel')
cmd.set('dash_color', 'purple', 'sb_name')
cmd.set('dash_gap', 0.3, 'sb_name')
cmd.set('dash_width', 2.5, 'sb_name')
cmd.hide('labels', 'sb_name')
```

## Optional Close Contacts

Close contacts are pocket-context hints, not specific chemical bonds.

- Detection cutoff: heavy-atom distance <= 4.0 A
- Display: gray dashed lines
- Default policy: report but do not render
- Opt-in policy: render only when the user explicitly requests hydrophobic or
  close contacts and both detector and renderer receive
  `--include-close-contacts`

Do not reconstruct the render interaction list from the detector's diagnostic
contact report. Use the exact final string printed under
`=== Interaction spec ===`.
