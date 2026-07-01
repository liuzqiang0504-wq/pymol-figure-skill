# PyMOL Pitfalls

Known bugs, silent failures, and IRON RULE workarounds discovered through painful trial and error.

## IRON RULE 1: Always use `color_deep`, never `color`

**Problem**: Maestro (.maegz) files embed per-atom `cartoon_color` and `ribbon_color`
values. When you call `cmd.color("palecyan", obj)`, these embedded colors silently
override the command. The protein remains colorful (secondary structure colors).

**Fix**:
```python
cmd.color_deep("palecyan", obj, 0)
```

This clears all per-atom color overrides and writes to the underlying atom color index.
The output will show: "Setting: ribbon_color unset for N atoms" and "Colored N atoms
and 1 object."

**Applies to**: Any file from Schrodinger software (Maestro, Glide, Prime, etc.).
Also affects some PDB files with embedded color data.

## IRON RULE 2: `hide('everything')` must come BEFORE all `show()` calls

**Problem**: If you call `hide('everything')` AFTER showing representations, the ligand
disappears and never comes back. The order of operations matters because PyMOL processes
hide/show commands sequentially.

**Correct pattern**:
```python
cmd.hide('everything')
cmd.show('cartoon', obj)       # protein
cmd.show('sticks', 'ligand')   # ligand
# ... all other show calls AFTER hide
```

**Wrong pattern**:
```python
cmd.show('sticks', 'ligand')
cmd.hide('everything')         # KILLS the ligand sticks!
cmd.show('cartoon', obj)       # protein comes back, ligand doesn't
```

## IRON RULE 3: Label pseudoatoms need `hide('spheres')` not `hide('everything')`

**Problem**: When you create a pseudoatom to hold a label and then call
`hide('everything', name)`, the label is destroyed along with the sphere.

**Fix**:
```python
cmd.pseudoatom('lbl_name', pos=[x, y, z])
cmd.label('lbl_name', '"LABEL TEXT"')
cmd.hide('spheres', 'lbl_name')   # correct — keeps the label
```

## IRON RULE 4: Protein and ligand are separate objects in .maegz

**Problem**: Maestro docking output files load as TWO PyMOL objects:
- `filename.maegz` or `filename_6XYA_hbond-opt` → protein structure
- `filename.LK-1-10.mol` → docked ligand

**Fix**: Always use `cmd.get_object_list()` and identify which is which by atom count
or organic selection:
```python
objects = cmd.get_object_list()
protein_obj = objects[0]   # typically the first (more atoms)
ligand_obj = objects[1]    # typically the second
```

## IRON RULE 5: `cmd.distance` creates persistent objects

**Problem**: Each `cmd.distance()` call creates a named distance object. If you re-run
the script or add new distances, old ones remain visible unless explicitly deleted.

**Fix**: Delete all previous distance objects before creating new ones, or use unique names.

## IRON RULE 6: Use `cmd.iterate_state` for coordinate access, not `cmd.iterate`

**Problem**: `cmd.iterate` does not provide `x, y, z` variables. Only
`cmd.iterate_state` or `cmd.alter_state` have access to atomic coordinates.

**Fix**:
```python
cmd.iterate_state(1, selection, 'stored.data.append((x, y, z))')
```

## IRON RULE 7: Disable non-existent colors gracefully

**Problem**: Some color names work in the GUI but not in scripting (e.g., "darkorange",
"lightgreen", "purple"). These cause `Error: Unknown color` and abort the script.

**Fix**: Use only well-known PyMOL color names or RGB tuples:
```python
cmd.set_color('my_color', [r, g, b])   # define custom color first
cmd.color('my_color', selection)
```

Safe color names: white, black, red, green, blue, yellow, cyan, magenta, orange,
pink, teal, gray, palecyan, lightblue, salmon, slate, violet, purpleblue, warmpink.

## IRON RULE 8: After `hide('everything')`, every visual element needs its own `show()` call

**Problem**: `cmd.hide('everything')` removes ALL representations from all objects.
If you call `util.cba()` or `cmd.color()` without a preceding `cmd.show()`, the target
stays invisible. The coloring succeeds silently but nothing appears on screen.

**Fix for ligand**: Always `cmd.show('sticks', ligand_obj)` BEFORE `util.cba()`:
```python
cmd.hide('everything')
cmd.show('cartoon', protein_obj)       # protein cartoon
cmd.show('sticks', ligand_obj)         # ligand — REQUIRED before CBA
util.cba(5274, ligand_obj, _self=cmd)
```

**Fix for metal ions in macro figure**: Metals must be explicitly shown as spheres:
```python
cmd.hide('everything')
cmd.show('cartoon', protein_obj)       # protein
cmd.show('sticks', ligand_obj)         # ligand
# Auto-detect and show metals:
cmd.select('_metals', f'{protein_obj} and (elem ZN or elem MG or elem MN or ...)')
if cmd.count_atoms('_metals') > 0:
    cmd.show('spheres', '_metals')
    cmd.set('sphere_scale', 0.3, '_metals')
    util.cnc('_metals', _self=cmd)
cmd.delete('_metals')
```

**Checklist when writing `setup_interaction_scene` or `setup_macro_scene`:**
- [ ] `show('cartoon', protein_obj)` — protein
- [ ] `show('sticks', ligand_obj)` — ligand (BEFORE util.cba)
- [ ] `show('sticks', res_sel)` — each interacting residue
- [ ] `show('spheres', metal_sel)` — each metal ion
- [ ] `show('surface', pocket_sel)` — binding pocket (macro only)

**Applies to**: Every figure. The most common symptom is "the ligand is invisible"
or "the metal ions disappeared in the macro figure."

## IRON RULE 9: Hide centroid pseudoatoms after creating dashed lines

**Problem**: π-π and cation-pi interactions use `cmd.pseudoatom()` to create ring
centroid markers, then draw dashed lines between them. If the centroid pseudoatoms
are not hidden after creating the distance, they remain visible as spheres or wire
points on the final figure — unprofessional and distracting.

**Fix**: Hide centroid pseudoatoms immediately after the distance is created:

```python
# π-π stacking
cmd.pseudoatom(pname, pos=prot_cent)
cmd.pseudoatom(lname, pos=lig_cent)
cmd.distance(dname, pname, lname)
cmd.hide('labels', dname)
cmd.hide('everything', pname)    # HIDE centroid marker
cmd.hide('everything', lname)    # HIDE centroid marker

# cation-pi
cmd.pseudoatom(f'cat_{ix}', pos=cation_pos)
cmd.pseudoatom(f'catpi_ring_{ix}', pos=ring_centroid)
cmd.distance(dname, f'cat_{ix}', f'catpi_ring_{ix}')
cmd.hide('labels', dname)
cmd.hide('everything', f'cat_{ix}')         # HIDE centroid marker
cmd.hide('everything', f'catpi_ring_{ix}')  # HIDE centroid marker
```

The distance object is independent — hiding the pseudoatoms does not affect the
dashed line visibility.

**Applies to**: π-π stacking and cation-pi interactions. The most common symptom is
"yellow/purple dots visible next to the dashed lines."
