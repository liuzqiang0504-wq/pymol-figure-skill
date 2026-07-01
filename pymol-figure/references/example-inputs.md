# Example Inputs

## Canonical Input Format

```
[FILE PATH] [INTERACTION SPEC]
```

Where INTERACTION SPEC is a comma-separated list of:
```
RESN RESI interaction_type
```

## Examples

### Example 1: Single π-π stacking
```
D:\docking\result.maegz PHE 1703 pi-pi
```
→ Renders interaction figure with PHE 1703 + macro figure

### Example 2: Multiple interactions
```
D:\docking\result.maegz PHE 1703 pi-pi, TYR 1719 pi-pi, ASP 189 hbond
```
→ Renders both figure types with all interactions shown

### Example 3: Chain-aware (PDB files with multiple chains)
```
D:\data\complex.pdb A/PHE/1703 pi-pi, A/ASP/189 hbond, B/LYS/42 hbond
```

### Example 4: Custom output directory
```
D:\docking\result.maegz PHE 1703 pi-pi --output D:\figures\paper1
```

### Example 5: Custom PyMOL path
```
D:\docking\result.maegz PHE 1703 pi-pi --pymol D:\PyMOL\python.exe
```

### Example 6: Metal coordination
```
D:\docking\result.maegz ZN 301 metal
```
→ Detects ligand donor atoms (N/O/S) within 3.0 A of ZN 301, shows metal as purple sphere

### Example 7: Cation-pi interaction
```
D:\docking\result.maegz PHE 1703 cation-pi
```
→ Finds nearby cation (metal / protonated amine) and draws purple dash to ring centroid

### Example 8: Salt bridge
```
D:\docking\result.maegz ASP 189 salt-bridge
```
→ Finds complementary charged groups and draws purple dashes

### Example 9: Mixed interaction types
```
D:\docking\result.maegz PHE 1703 pi-pi, ASP 189 hbond, ZN 301 metal, HIS 1722 cation-pi, LYS 42 salt-bridge
```

## Supported Interaction Types

| Keyword | Meaning | Dash color |
|---------|---------|------------|
| `pi-pi`, `pi-stacking` | π-π aromatic stacking | yellow |
| `hbond`, `h-bond` | Hydrogen bond | green |
| `metal`, `metal-coord` | Metal coordination | purple |
| `cation-pi`, `cationpi` | Cation-π interaction | purple |
| `salt-bridge`, `saltbridge` | Salt bridge | purple |

## What the Skill Produces

For any valid input, the skill generates:

```
<output_dir>/
├── interaction_1.png    # 作用图 angle 1
├── interaction_2.png    # 作用图 angle 2
├── interaction_3.png    # 作用图 angle 3
├── macro_1.png          # 宏观图 angle 1
├── macro_2.png          # 宏观图 angle 2
├── macro_3.png          # 宏观图 angle 3
└── session.pse          # PyMOL session file
```

## Default Output Directory

If `--output` is not specified, files go to `<input_file_dir>/pymol_figures/`.
