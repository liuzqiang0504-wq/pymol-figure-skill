# Optional 2D Interaction Diagrams

Use this workflow only after the user accepts the optional 2D output or
explicitly requests it.

## Inputs

Require:

- receptor `.pdb`
- docked ligand poses in `.sdf` or `.mol2`
- optional one-based pose index, default 1

Prefer SDF or MOL2 because they preserve ligand bond orders and formal charges.
Do not pass PDBQT directly. Convert it with a chemistry-aware tool such as
Meeko first, then use the converted SDF.

## Command

Run with the same Python environment that contains RDKit and Pillow:

```powershell
& $env:PYMOL_FIGURE_RDKIT_PYTHON `
  "SKILL_DIR/scripts/generate_2d_interactions.py" `
  --protein "receptor.pdb" `
  --ligand "docked.sdf" `
  --pose-index 1 `
  --output "<DIR>\2d"
```

If `PYMOL_FIGURE_RDKIT_PYTHON` is unset, use the verified RDKit-capable Python,
such as `py -3.12`.

## Default Visual Rules

- Show all standard amino acids within 5.0 A as pocket context.
- Anchor each residue to its nearest ligand atom and arrange labels radially
  with collision avoidance.
- Treat label positions as a schematic, not a 3D projection or distance scale.
- Render detected interaction residues in larger bold red Arial text.
- Render non-interacting pocket-context residues in smaller blue Arial text.
- Draw hydrogen bonds, pi stacking, and salt bridges with distances.
- Show hydrophobic contacts with red LigPlot-style fan symbols.
- Hide hydrophobic contact lines by default. Add
  `--draw-hydrophobic-lines` only when the user explicitly requests them.
- Preserve assigned tetrahedral stereochemistry from SDF/MOL2 input.
- Draw `BEGINWEDGE` stereobonds as solid wedges and `BEGINDASH` stereobonds as
  hashed wedges, with the narrow end at the stereocenter.
- Do not assign a wedge to an unassigned stereocenter or infer missing
  stereochemistry from the 2D layout.
- **Pre-export quality gate (mandatory):** before the PNG is written, validate
  that residue and distance labels do not overlap one another, the ligand
  drawing, or hydrophobic symbols, and remain inside the canvas. If needed,
  automatically re-layout with a larger canvas (up to two retries); fail rather
  than exporting an unresolved collision.
- **Pre-export stereochemistry gate (mandatory):** confirm each assigned R/S
  centre retains its CIP assignment and is represented by an adjacent solid or
  hashed wedge. Confirm unassigned centres were not given a stereobond.
- Export a transparent PNG by default. Add `--opaque-background` for white.

This output is LigPlot-inspired but does not invoke, redistribute, or require
LigPlot+.

## Outputs And Verification

Require all three files:

- `interaction_2d.png`
- `interaction_2d.json`
- `interaction_spec.txt`

Verify that the PNG has an alpha channel, residue names are readable, key and
context residues remain visually distinct, and distance labels do not cross
lines, atoms, labels, or hydrophobic symbols. When the ligand contains assigned
chiral centers, also verify that RDKit reports the expected CIP assignments and
that the prepared drawing molecule contains the expected solid/hashed
stereobonds.

`interaction_2d.json` includes a `quality_check` record with the final canvas
size, number of layout attempts, collision-check result, and stereochemistry
verification report. Treat a missing or failed record as a failed render.

After changing the 2D renderer, run:

```powershell
py -3.12 "SKILL_DIR/scripts/test_2d_stereochemistry.py"
```
