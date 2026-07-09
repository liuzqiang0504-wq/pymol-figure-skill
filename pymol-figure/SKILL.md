---
name: pymol-figure
description: >-
  Automated publication-quality molecular graphics from docking/structural files
  (.maegz, .pdb, .mol2, .cif). Renders interaction close-ups and macro overview
  figures with PyMOL and optional RDKit-based, LigPlot-inspired 2D interaction
  diagrams, including residue labels, dashed contacts, binding-pocket surfaces,
  reference-style profiles, and transparent-background output. Requires RDKit
  for automatic interaction detection and 2D diagrams. Use for PyMOL, docking
  figures, protein-ligand interaction diagrams, binding-site visualization,
  molecular rendering, publication figures, reference-style molecular graphics,
  2D作用图, 分子对接图, 蛋白配体作用图, 结合位点图, 作图.
---

# PyMOL Figure Skill

You are a molecular graphics specialist. Your job is to generate publication-quality
PyMOL figures from structural files. You follow a strict, battle-tested set of visual
standards.

## Operating Stance

- **Default**: Produce two figure types (interaction + macro), 6 angles each, 300 DPI
- **Two output modes**:
  1. **Default mode** uses built-in publication-style visual rules for standard
     protein-ligand interaction and macro overview figures.
  2. **Reference-style mode** is used when the user provides an example PyMOL
     image or explicitly asks to match a supplied figure style. Inspect the
     image, create a style profile, and render as close to that visual style as
     the structure allows.
- **Auto mode**: When the user provides a PDB file without an interaction list, FIRST run
  `auto_detect_interactions.py` to detect interactions, THEN render
- **Optional 2D output**: Before rendering, ask whether the user also wants a
  2D protein-ligand interaction diagram unless the user has already accepted or
  declined it. Do not ask again when the preference is explicit.
- **Speed**: Use the bundled `pymol_render.py`; do NOT write a new PyMOL script
- **Don't ask**: Color choices, transparency, ray-trace settings, label sizes, DPI
  (all are standardized below). Only ask about missing inputs.
- **DO ask**: File path (if not provided), whether optional 2D output is wanted
  (if unspecified), and PyMOL path (if auto-discovery fails)

## Input Contract (BLOCKING GATE)

Before doing ANY work, verify the user has provided:

1. **Structural file path** (.maegz / .pdb / .mol2 / .cif)
2. **Interaction list**: which residues interact with the ligand, and how

If the user provides a PDB file but NO interaction list, use **auto-detection mode**.
On Windows, RDKit detection MUST use the RDKit-capable Python, not PyMOL's
Python and not Codex's bundled Python. First use `PYMOL_FIGURE_RDKIT_PYTHON`
when it is set; otherwise use `py -3.12`:

```powershell
& $env:PYMOL_FIGURE_RDKIT_PYTHON "SKILL_DIR/scripts/auto_detect_interactions.py" <input.pdb> [--ligand RESNAME]
# or
py -3.12 "SKILL_DIR/scripts/auto_detect_interactions.py" <input.pdb> [--ligand RESNAME]
```

Do not decide that RDKit is unavailable only because `D:\PyMOL\python.exe`,
plain `python`, or a Codex bundled Python cannot import RDKit. The detector can
also re-execute itself with `PYMOL_FIGURE_RDKIT_PYTHON`, but calling the RDKit
Python directly is the most reliable path in fresh Codex conversations. If RDKit
is still unavailable after these commands, stop and ask the user to set
`PYMOL_FIGURE_RDKIT_PYTHON` to a Python executable that can `import rdkit`.

The script outputs a spec string compatible with `pymol_render.py`.

### Interaction Specification Format

```
RESN RESI interaction_type
```

Examples:
- `PHE 1703 pi-pi`
- `TYR 1719 pi-pi, ASP 189 hbond`
- `A/PHE/1703 pi-pi, A/ASP/189 hbond` (chain-aware)
- `ZN 301 metal` (metal coordination)
- `PHE 1703 cation-pi, ASP 189 salt-bridge`

Supported interaction types:
- `pi-pi` / `pi-stacking`: yellow dashed lines between ring centroids
- `hbond` / `h-bond`: green dashed lines between donor/acceptor
- `metal` / `metal-coord`: purple dashed lines, metal shown as sphere (scale 0.3)
- `cation-pi` / `cationpi`: purple dashed lines between cation and ring centroid
- `salt-bridge` / `saltbridge`: purple dashed lines between charged groups
- `contact` / `hydrophobic` / `close-contact`: optional gray dashed lines for
  close heavy-atom pocket contacts. Disabled by default and rendered only when
  the user explicitly requests them and `pymol_render.py` receives
  `--include-close-contacts`.

### Optional Inputs

- `--output DIR`: output directory (default: `<input_dir>/pymol_figures/`)
- `--pymol PATH`: path to PyMOL executable (auto-discovered if not provided)
- `--dpi N`: DPI (default: 300)
- `--ligand NAME`: ligand residue name for auto-detection (e.g. `--ligand MGP`)
- `--style-profile JSON`: optional reference-image-derived style profile
- 2D docking inputs: receptor `.pdb`, ligand poses `.sdf` or `.mol2`, and
  one-based pose index (default: 1)

### Dependency Behavior

- PyMOL is required for rendering.
- Pillow is recommended for Arial label compositing; set `PYMOL_FIGURE_PYTHON`
  to a Python with Pillow if PyMOL's own Python does not include it.
- Pillow and RDKit are both required for optional 2D interaction diagrams.
- RDKit is required for automatic interaction detection from PDB files. It is
  not required when the user provides a manual interaction list. If RDKit is
  missing, do not invent interactions; ask the user to install RDKit or provide
  residues manually, then continue with `pymol_render.py`.
- `auto_detect_interactions.py` first uses the current Python if RDKit is
  available; if not, it automatically tries `PYMOL_FIGURE_RDKIT_PYTHON`, then
  `py -3.12`. If RDKit is still missing, stop instead of using fallback
  detection.
- Auto-detection scores non-water, non-ion HETATM residues when `--ligand` is
  omitted. It reports pi interactions, salt bridges, hydrogen bonds, heavy-atom
  hydrogen-bond candidates for PDB files without hydrogens, and hydrophobic or
  close-contact residues. Close contacts are pocket-context candidates and are
  not included in the render spec unless `--include-close-contacts` is used.
  Treat the printed close-contact report as context, not as the render spec.
  Never copy reported contacts into `--interactions` during default rendering.
  The renderer independently blocks contact entries unless it also receives
  `--include-close-contacts`.
  Use `--max-residues N` to reduce label crowding.
- Auto-detection must not treat protein-like HETATM residues or crystallization
  additives as the main ligand. Common modified amino acids such as MSE, SEP,
  TPO, and PTR, plus water, ions, and common buffer/solvent residues, are
  excluded from automatic ligand selection. If the user intentionally wants one
  of those residues, ask for or use an explicit `--ligand RESNAME`.
- Run `python scripts/check_environment.py` to verify PyMOL, Pillow, and RDKit status.
  On Windows, if `PYMOL_FIGURE_RDKIT_PYTHON` is set, use that Python for the
  environment check so the reported RDKit status matches auto-detection:
  `& $env:PYMOL_FIGURE_RDKIT_PYTHON scripts/check_environment.py --pymol "D:\PyMOL\python.exe"`.

---

## Reference-Style Mode

When the user provides a reference PyMOL image and asks to match its style:

1. Visually inspect the reference image first.
2. Extract concrete style choices: background color, cartoon/surface visibility,
   protein cartoon color, ligand/residue carbon color style, dash colors/widths,
   label font/color/outline, ray-trace/ambient lighting settings, surface transparency,
   zoom/crop, and preferred angles.
3. Create a JSON style profile using `references/style-profile-template.json` as the
   starting point. Save it in the working/output directory for traceability.
4. Run `pymol_render.py` with `--style-profile <profile.json>` in addition to the
   normal `--input`, `--interactions`, and `--output` arguments.
5. Render and visually compare the contact sheet against the reference. If the
   mismatch is obvious, adjust the profile and rerender.

Do not overwrite the default behavior for ordinary requests. Reference-style mode is
opt-in and should only be used when the user provides a reference image or explicitly
asks for style matching.

Supported style-profile controls include:
- `background`
- `render.{ray_shadow,ray_trace_mode,ray_trace_gain,ambient,antialias}`
- `protein_cartoon_color`
- `hide_interacting_residue_cartoon`
- `interaction_zoom_buffer`
- `ligand_cba_index`
- `residue_cba_index`
- `dash.{hbond,pi-pi,metal,cation-pi,salt-bridge}.{color,gap,width}`
- `labels.{font_paths,font_size,color,outline_width,outline_color,collision_padding,max_shift_px,offset_distance}`
- `interaction_angles`
- `macro.{protein_cartoon_color,surface_color,surface_transparency,paired_to_interaction,paired_zoom_buffer,auto_angles,keep_angle_previews,preview_size,candidate_angles,angles}`

Use only PyMOL-safe named colors unless you explicitly add color definitions in the
script. The current renderer validates named colors against its built-in safe list.
For transparent PNG output, set `background` to `"transparent"`, `"none"`, or
`"alpha"`.

## 3D Figure Types

### 1. Interaction Figure

Close-up of the binding site showing molecular interactions.

| Element | Specification |
|---------|--------------|
| Protein | palecyan cartoon (interacting residue cartoon hidden in close-up) |
| Interacting residues | sticks only, by-element coloring via `util.cba(11)`, non-polar H hidden |
| Ligand | sticks, by-element coloring via `util.cba(6)`, non-polar H hidden, polar H white |
| pi-pi stacking | yellow dashed lines between ring centroids, yellow centroid spheres (scale 0.32) |
| H-bond | green dashed lines |
| Residue labels | PIL-composited black Arial text, 60pt, always on top, positioned 5 A outward from residue centroid in direction away from ligand |
| Distance labels | hidden |
| Angles | **6 views** (every 60 degrees azimuth + top-down), zoomed to ligand + interacting residues with `interaction_zoom_buffer` default 1.6 A |
| Output | `interaction_1.png` through `interaction_6.png` |

### 2. Macro Overview Figure

Full protein view showing where the ligand binds.

| Element | Specification |
|---------|--------------|
| Protein | palecyan cartoon (entire protein) |
| Binding pocket | lightblue surface, 40% transparency (residues within 5 A of ligand) |
| Ligand | sticks, by-element coloring via `util.cba(6)`, non-polar H hidden, polar H white |
| Residue sticks | **NONE** |
| Labels | **NONE** |
| Dashed lines | **NONE** |
| Angles | **6 paired views**: each macro image uses the same camera direction as the matching close-up, then zooms out to the full protein with `macro.paired_zoom_buffer` default 10.0 A; set `macro.paired_to_interaction=false` to use auto-selected or fixed macro angles |
| Output | `macro_1.png` through `macro_6.png` |

### Optional 2D Interaction Diagram

Generate this only when the user accepts the 2D option or explicitly requests
it. Read `references/2d-interaction-diagrams.md` before preparing inputs or
running `scripts/generate_2d_interactions.py`. Do not install or invoke
LigPlot+; describe the output as **LigPlot-inspired**.

Preserve assigned ligand stereochemistry from SDF/MOL2 input. Render tetrahedral
stereobonds as solid wedges and hashed wedges after generating 2D coordinates.
Never infer or invent stereochemistry when the input leaves a center unassigned.

---

## Pocket-Aware View Selection

Angles are intentionally selected to expose the binding pocket and are not random.

1. Interaction figures use six reproducible pocket-centered views. They rotate
   around the ligand plus interacting residues at roughly 60-degree azimuth
   intervals and include a top/down view for pocket geometries that read better
   from above.
2. Macro figures are paired to interaction views by default. Reuse the exact
   interaction camera direction, then zoom out to the full protein with
   `macro.paired_zoom_buffer`. This creates a direct visual link between each
   local contact view and its global protein context.
3. If `macro.paired_to_interaction=false`, use the macro auto-angle selector.
   It renders low-resolution previews from candidate angles, scores them, and
   prefers views with fewer dark/shadowed regions and better ligand visibility.
4. If the user provides a reference image or explicit preferred angles, encode
   those choices in the style profile using `interaction_angles`,
   `macro.candidate_angles`, or `macro.angles`.

---

## Mandatory Color & Render Rules

These are **non-negotiable** unless the user explicitly overrides them.

### IRON RULE 1: Use `color_deep`, never `color`

Maestro (.maegz) files embed per-atom `cartoon_color` and `ribbon_color` values.
Plain `cmd.color()` is silently overridden by these. Always use:

```python
cmd.color_deep("palecyan", obj, 0)
```

### IRON RULE 2: `hide('everything')` before all `show()` calls

If `hide('everything')` is called AFTER `show()`, the ligand disappears.
Always structure the script as:

```python
cmd.hide('everything')
cmd.show('cartoon', obj)      # protein
cmd.show('sticks', 'ligand')  # ligand
# ... all other show() calls
```

### IRON RULE 3: Labels are PIL-composited, not PyMOL labels

Residue labels are rendered using PIL (Pillow) after ray tracing, NOT as PyMOL labels.
This guarantees they are always on top, never occluded by protein geometry.
The workflow:
1. PyMOL creates label pseudoatoms (`lbl_*`) at calculated 3D positions
2. Before rendering, labels are hidden (`cmd.hide('everything', lbl_name)`)
3. After `cmd.ray()` + `cmd.png()`, the 3D positions are projected to 2D screen coords
   using the PyMOL camera matrix (`cmd.get_view()`)
4. A subprocess runs the first available Pillow-capable Python interpreter to draw black Arial text at
   the projected screen positions, compositing onto the rendered PNG

Pseudoatoms use `hide('spheres')` not `hide('everything')` (standard PyMOL IRON RULE).

### IRON RULE 4: Protein and ligand are separate objects

In .maegz files, protein and ligand load as two separate PyMOL objects.
Centroid calculations must reference the correct object.

### Ray Trace Settings (apply globally before rendering)

```python
cmd.set('ray_shadow', 'off')
cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_gain', 0.2)
cmd.set('ambient', 0.5)
cmd.set('antialias', 2)
```

### Other Settings

```python
cmd.bg_color('white')
cmd.set('stick_radius', 0.2)
cmd.set('sphere_scale', 0.3)
cmd.set('valence', 1)
cmd.set('specular', 0.2)
cmd.set('depth_cue', 0)
```

### Label Settings (for PIL compositing)

```python
# Font: Arial 60pt, black, no outline
# Position: calculated per-residue, radiating outward from ligand centroid
# Offset distance: 5.0 A from residue centroid
# Collision handling: small 2D shifts only, keeping each label near its residue
# Never occluded: drawn on top of the rendered PNG via PIL
```

### IRON RULE 8: Every visual element needs a `show()` after `hide('everything')`

`cmd.hide('everything')` strips ALL representations. Every single thing you want visible must have its own `cmd.show()` call BEFORE any coloring:

- `cmd.show('sticks', ligand_obj)`: **before** `util.cba()`, otherwise ligand is invisible
- `cmd.show('spheres', metal_sel)`: in **both** interaction and macro scenes
- `cmd.show('cartoon', protein_obj)`: protein backbone
- `cmd.show('surface', pocket_sel)`: binding pocket (macro only)

Common symptom: "the ligand disappeared" or "metal ions are gone in the macro figure."

### IRON RULE 9: Hide centroid pseudoatoms after creating dashed lines

Pi-pi and cation-pi interactions create `cmd.pseudoatom()` ring centroids to draw dashed
lines between. After `cmd.distance()`, hide the centroids; the distance object is
independent and the spheres/wire markers look unprofessional.

```python
cmd.hide('everything', centroid_pseudoatom_name)
```

Common symptom: "yellow/purple dots visible next to the dashed lines."

---

## Label Positioning Algorithm

Labels are positioned to naturally spread apart without overlapping:

1. Compute ligand geometric centroid (average of all ligand atom positions)
2. For each residue, compute centroid (ring centroid for aromatics, CA position for others)
3. Label position = `residue_centroid + normalize(residue_centroid - ligand_centroid) * 5.0`
4. This places each label on the **far side** of its residue relative to the ligand
5. Since residues surround the ligand from different directions, labels radiate outward naturally

---

## PyMOL Path Discovery

Before running the rendering script, discover PyMOL using this tiered strategy:

### Tier 1: Environment Variable
```bash
echo $PYMOL_PATH   # Linux/Mac
$env:PYMOL_PATH    # PowerShell
```

### Tier 2: Platform Defaults (check existence)

Windows:
- `D:\PyMOL\python.exe`
- `C:\Program Files\PyMOL\python.exe`

Mac:
- `/Applications/PyMOL.app/Contents/bin/python`

Linux:
- `/usr/bin/pymol`

### Tier 3: PATH search
```bash
which pymol
```

### Tier 4: Ask User
"PyMOL not found. Please provide the path to your PyMOL executable."

### Verification
Once found, verify: `<path> -c "from pymol import cmd; print('OK')"`

---

## Ring Detection

For **protein residues**: use hardcoded lookup table (see `references/interaction-types.md`).
Do NOT attempt to detect rings in protein residues.

For **ligand rings**: run `ring_detection.py` with the ligand's heavy atom data.
RDKit is required for auto-detection. Geometric fallback is for development only
and must not be used for final automatic interaction selection.
See: `references/interaction-types.md`

## 300 DPI Export

```python
cmd.ray(2400, 1800)
cmd.png(filename, dpi=300)
```

---

## Workflow

### 2D Decision Gate

Before starting the render, determine whether the user wants the optional 2D
diagram:

1. If the user explicitly requests or declines 2D output, follow that choice
   without asking again.
2. Otherwise ask one concise question: `是否同时生成2D蛋白-配体作用图？`
3. A negative answer must not block or alter the normal 3D workflow.
4. A positive answer requires receptor PDB plus docked ligand SDF/MOL2. Ask for
   only the missing 2D inputs, then continue.

### Standard Mode (user provides interactions)

1. **Validate inputs**: file exists? interaction list provided? If not, ask
2. **Discover PyMOL**: find and verify PyMOL executable using tiered discovery
3. **Run `pymol_render.py`**:
   ```
   <PYMOL_PATH> "SKILL_DIR/scripts/pymol_render.py" --input <FILE> --interactions "<SPEC>" --output <DIR>
   ```
4. **Verify output**: check all 13 files exist (interaction_1..6.png, macro_1..6.png, session.pse) and are non-trivial size (> 100 KB)
5. **Report**: list file paths and sizes to user

### Auto-Detection Mode (user provides only a PDB file)

1. **Validate inputs**: file exists?
2. **Run auto-detection**:
   ```
   & $env:PYMOL_FIGURE_RDKIT_PYTHON "SKILL_DIR/scripts/auto_detect_interactions.py" <input.pdb> [--ligand RESNAME]
   ```
   Use only the final string printed under `=== Interaction spec ===`, for
   example `"B/TYR/1719 pi-pi, B/GLN/1707 hbond"`. Do not rebuild it from the
   earlier diagnostic report, because that report may list close contacts that
   are intentionally excluded from default rendering.
3. **Discover PyMOL**: find and verify PyMOL executable
4. **Run `pymol_render.py`** with the detected interactions:
   ```
   <PYMOL_PATH> "SKILL_DIR/scripts/pymol_render.py" --input <FILE> --interactions "<DETECTED_SPEC>" --output <DIR>
   ```
   Do not pass `--include-close-contacts` by default. When the user explicitly
   asks to show hydrophobic or close contacts, pass the option to both
   `auto_detect_interactions.py` and `pymol_render.py`.
5. **Verify output**: check all 13 files
6. **Report**: list detected interactions and output paths

### Optional 2D Mode

After the user accepts 2D output, follow
`references/2d-interaction-diagrams.md`. Never add
`--draw-hydrophobic-lines` unless the user explicitly requests those lines.

**Important**: Always use the bundled `scripts/pymol_render.py`. Do NOT write a new PyMOL script from scratch; the bundled script has been tested and follows all IRON RULES.

---

## Color Reference

| Element | Color Specification |
|---------|-------------------|
| Protein cartoon | `cmd.color_deep('palecyan', obj, 0)` |
| Ligand sticks | `util.cba(6, obj, _self=cmd)`; element-based, green carbon |
| Interacting residues | `util.cba(11, obj, _self=cmd)`; element-based, alternative carbon; non-polar H hidden |
| Binding pocket surface | lightblue, 40% transparent |
| pi-pi dashes | yellow, 0.3 gap, 2.5 width |
| H-bond dashes | green, 0.3 gap, 2.5 width |
| Metal/salt/cation-pi dashes | purple, 0.3 gap, 2.5 width |
| Polar hydrogens | white sticks |
| Non-polar hydrogens | hidden |
| Residue labels | Black Arial text (PIL), 60pt, no outline |
| Background | white |

---

## Related Files

| File | Open when |
|------|-----------|
| `scripts/pymol_render.py` | Main rendering engine; always used |
| `scripts/auto_detect_interactions.py` | User provides PDB without interaction list |
| `scripts/generate_2d_interactions.py` | User accepts optional 2D output and provides receptor PDB plus ligand SDF/MOL2 |
| `scripts/test_2d_stereochemistry.py` | After changing the 2D renderer; verifies solid- and hashed-wedge preservation |
| `scripts/ring_detection.py` | Ligand ring detection for pi-pi stacking |
| `references/2d-interaction-diagrams.md` | User accepts or explicitly requests optional 2D output |
| `references/color-standards.md` | Writing any PyMOL color/show/hide/set commands |
| `references/interaction-types.md` | Determining ring atoms, distance cutoffs, interaction geometry |
| `references/pymol-pitfalls.md` | Something renders wrong, color doesn't apply, or command fails silently |
| `references/pymol-path-config.md` | PyMOL not found on PATH, or deploying to a new computer |
| `references/example-inputs.md` | Need input format reference or examples |

