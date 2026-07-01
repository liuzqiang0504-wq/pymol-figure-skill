# Color Standards

Single source of truth for all PyMOL visual parameters.

## Color Commands

### Protein Cartoon

```python
cmd.color_deep("palecyan", obj, 0)
```

**IMPORTANT**: Always use `color_deep`, never `color`. Maestro (.maegz) files embed
per-atom `cartoon_color` and `ribbon_color` values. `color` is silently overridden.

### Interacting Residue Sticks

```python
cmd.color_deep("lightblue", selection, 0)
```

### Ligand (by element)

```python
from pymol import util
util.cba(5274, 'ligand', _self=cmd)

# Hide non-polar H
cmd.select('lig_h_np', f'{lig_obj} and elem H and (neighbor elem C)')
cmd.hide('everything', 'lig_h_np')

# Polar H white
cmd.show('sticks', f'{lig_obj} and elem H and not (neighbor elem C)')
cmd.color('white', f'{lig_obj} and elem H and not (neighbor elem C)')
```

### π-π Stacking Dashes

```python
cmd.distance('pi_name', 'centroid1', 'centroid2')
cmd.set('dash_color', 'yellow', 'pi_name')
cmd.set('dash_gap', 0.3, 'pi_name')
cmd.set('dash_width', 2.5, 'pi_name')
cmd.hide('labels', 'pi_name')
```

### H-bond Dashes

```python
cmd.distance('hbond_name', 'donor_atom', 'acceptor_atom')
cmd.set('dash_color', 'green', 'hbond_name')
cmd.set('dash_gap', 0.3, 'hbond_name')
cmd.set('dash_width', 2.5, 'hbond_name')
cmd.hide('labels', 'hbond_name')
```

### Metal Coordination Dashes (purple)

```python
cmd.distance('metal_name', 'metal_sel', 'donor_sel')
cmd.set('dash_color', 'purple', 'metal_name')
cmd.set('dash_gap', 0.3, 'metal_name')
cmd.set('dash_width', 2.5, 'metal_name')
cmd.hide('labels', 'metal_name')
```

### Metal Ion Sphere

```python
cmd.show('spheres', metal_selection)
cmd.set('sphere_scale', 0.3, metal_selection)
util.cnc(metal_selection, _self=cmd)  # color by chain
```

### Cation-pi Dashes (purple)

```python
cmd.distance('catpi_name', 'cation_sel', 'centroid_sel')
cmd.set('dash_color', 'purple', 'catpi_name')
cmd.set('dash_gap', 0.3, 'catpi_name')
cmd.set('dash_width', 2.5, 'catpi_name')
cmd.hide('labels', 'catpi_name')
```

### Salt Bridge Dashes (purple)

```python
cmd.distance('sb_name', 'acidic_sel', 'basic_sel')
cmd.set('dash_color', 'purple', 'sb_name')
cmd.set('dash_gap', 0.3, 'sb_name')
cmd.set('dash_width', 2.5, 'sb_name')
cmd.hide('labels', 'sb_name')
```

### Cation-pi Centroid Pseudoatoms (HIDDEN after distance creation)

```python
cmd.pseudoatom('cent_name', pos=[x, y, z])
# ... create distance ...
cmd.hide('everything', 'cent_name')  # only the dashed line should be visible
```

### Ring Centroid Pseudoatoms (HIDDEN after distance creation)

```python
cmd.pseudoatom('cent_name', pos=[x, y, z])
# ... create distance ...
cmd.hide('everything', 'cent_name')  # only the dashed line should be visible
```

### Binding Pocket Surface (Macro Figure)

```python
cmd.select('pocket_5A', f'{obj} within 5 of {lig_obj}')
cmd.show('surface', 'pocket_5A')
cmd.set('surface_color', 'lightblue', 'pocket_5A')
cmd.set('transparency', 0.4, 'pocket_5A')
```

### Residue Labels

```python
cmd.set('label_color', 'black')
cmd.set('label_size', 26)
# Place at pseudoatom offset from residue to avoid blocking dashes
cmd.pseudoatom('lbl_res', pos=[x_offset, y_offset, z_offset])
cmd.label('lbl_res', '"PHE 1703"')
cmd.hide('spheres', 'lbl_res')  # hide sphere, keep label
```

## Ray Trace Settings

```python
cmd.set('ray_shadow', 'off')
cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_gain', 0.2)
cmd.set('ambient', 0.5)
cmd.set('antialias', 2)
```

## Global Settings

```python
cmd.bg_color('white')
cmd.set('stick_radius', 0.2)
cmd.set('sphere_scale', 0.3)
cmd.set('valence', 1)
cmd.set('specular', 0.2)
cmd.set('depth_cue', 0)
```

## 300 DPI Export

```python
cmd.ray(2400, 1800)
cmd.png(filename, dpi=300)
```

## 3-Angle Rendering

```python
cmd.zoom(selection, buffer)
# Angle 1
cmd.turn('x', -15); cmd.turn('y', 30)
cmd.ray(2400, 1800); cmd.png('prefix_1.png', dpi=300)
# Angle 2 (+120 degrees)
cmd.turn('y', 120)
cmd.ray(2400, 1800); cmd.png('prefix_2.png', dpi=300)
# Angle 3 (+120 degrees)
cmd.turn('y', 120)
cmd.ray(2400, 1800); cmd.png('prefix_3.png', dpi=300)
```
