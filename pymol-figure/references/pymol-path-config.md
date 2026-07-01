# PyMOL Path Configuration

How to find PyMOL on any system.

## Four-Tier Discovery

### Tier 1: Environment Variable

Check if `PYMOL_PATH` is set:
```powershell
$env:PYMOL_PATH        # Windows PowerShell
```
```bash
echo $PYMOL_PATH        # Linux/Mac
```

### Tier 2: Platform-Specific Defaults

Check these paths in order and use the first one that exists:

**Windows**:
- `D:\PyMOL\python.exe`
- `C:\Program Files\PyMOL\python.exe`
- `C:\Program Files\PyMOL2\python.exe`
- `C:\PyMOL\python.exe`

**macOS**:
- `/Applications/PyMOL.app/Contents/bin/python`
- `/opt/homebrew/bin/pymol`

**Linux**:
- `/usr/bin/pymol`
- `/usr/local/bin/pymol`
- `/opt/pymol/bin/python`

### Tier 3: PATH Search

```bash
which pymol        # bash
where pymol        # PowerShell
```

### Tier 4: Ask User

If all tiers fail, ask:
"PyMOL not found. Please provide the full path to your PyMOL executable
(e.g., `D:\PyMOL\python.exe` on Windows or `/usr/bin/pymol` on Linux)."

## Verification

Once found, verify with:
```bash
"<PYMOL_PATH>" -c "from pymol import cmd; print('OK')"
```

If this fails, the path is wrong — escalate to Tier 4.

## Site-packages Location

PyMOL's Python has its own site-packages. Common locations:
- Windows: `<PyMOL dir>\Lib\site-packages`
- macOS: `/Applications/PyMOL.app/Contents/lib/python3.x/site-packages`
- Linux: `/usr/lib/python3.x/site-packages/pymol`

To use PyMOL's Python environment: always invoke via PyMOL's own python executable,
not the system python.

## Usage in Scripts

```python
import sys
sys.path.insert(0, r'<PYMOL_DIR>\Lib\site-packages')
import pymol
from pymol import cmd, stored, util
pymol.finish_launching(['pymol', '-cq'])
```
