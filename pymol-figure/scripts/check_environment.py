"""
Environment check for the PyMOL figure skill.

Run with any Python:
    python scripts/check_environment.py

Optional:
    python scripts/check_environment.py --pymol "D:/PyMOL/python.exe"
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Check PyMOL figure dependencies")
    parser.add_argument("--pymol", default=None, help="Path to PyMOL Python/executable")
    return parser.parse_args()


def _candidate_pymol_paths(explicit: str | None):
    candidates = []
    for value in (explicit, os.environ.get("PYMOL_PATH")):
        if value:
            candidates.append(value)
    candidates.extend([
        r"D:\PyMOL\python.exe",
        r"C:\Program Files\PyMOL\python.exe",
        "/Applications/PyMOL.app/Contents/bin/python",
        "/usr/bin/pymol",
    ])
    for name in ("pymol", "pymol.exe"):
        found = shutil.which(name)
        if found:
            candidates.append(found)

    unique = []
    seen = set()
    for path in candidates:
        key = os.path.normcase(os.path.abspath(path)) if os.path.sep in path else path
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def _check_module(module: str):
    try:
        __import__(module)
        return True, None
    except Exception as exc:
        return False, str(exc)


def _split_command(value: str | None):
    if not value:
        return None
    try:
        parts = shlex.split(value)
    except ValueError:
        parts = [value]
    return parts or None


def _read_user_env_var(name: str):
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


def _candidate_rdkit_pythons():
    preferred = []
    delayed = []
    seen = set()

    for value in (
        _read_user_env_var("PYMOL_FIGURE_RDKIT_PYTHON"),
        os.environ.get("PYMOL_FIGURE_RDKIT_PYTHON"),
    ):
        _add_unique_command(preferred, delayed, seen, _split_command(value))

    if os.name == "nt":
        local_appdata = os.environ.get("LOCALAPPDATA")
        for path in (
            os.path.join(local_appdata or "", "Programs", "Python", "Python312", "python.exe"),
            r"C:\Python312\python.exe",
            r"C:\Program Files\Python312\python.exe",
        ):
            if path and os.path.exists(path):
                _add_unique_command(preferred, delayed, seen, [path])

    for cmdline in ([sys.executable], ["py", "-3.12"], ["python3.12"], ["python"]):
        _add_unique_command(preferred, delayed, seen, cmdline)

    return preferred + delayed

def _check_module_with_command(cmdline, module: str):
    try:
        result = subprocess.run(
            cmdline + ["-c", f"import {module}; print('OK')"],
            capture_output=True,
            text=True,
            timeout=12,
        )
    except Exception as exc:
        return False, str(exc)
    if result.returncode == 0:
        return True, None
    return False, (result.stderr or result.stdout).strip()


def _check_pymol(candidates):
    script = "from pymol import cmd; print('PyMOL import OK')"
    for candidate in candidates:
        if os.path.sep in candidate and not os.path.exists(candidate):
            continue
        try:
            result = subprocess.run(
                [candidate, "-c", script],
                capture_output=True,
                text=True,
                timeout=20,
            )
        except Exception:
            continue
        if result.returncode == 0:
            return True, candidate
    return False, None


def main():
    args = parse_args()
    print("PyMOL figure skill environment check")
    print(f"Current Python: {sys.executable}")

    pymol_ok, pymol_path = _check_pymol(_candidate_pymol_paths(args.pymol))
    if pymol_ok:
        print(f"[OK] PyMOL found: {pymol_path}")
    else:
        print("[MISSING] PyMOL not found. Set PYMOL_PATH or pass --pymol.")

    pillow_ok, pillow_error = _check_module("PIL")
    if pillow_ok:
        print("[OK] Pillow available for Arial label overlay")
    else:
        print(f"[RECOMMENDED] Pillow not available: {pillow_error}")
        print("              Install with: python -m pip install Pillow")

    rdkit_found = False
    rdkit_errors = []
    for cmdline in _candidate_rdkit_pythons():
        ok, error = _check_module_with_command(cmdline, "rdkit")
        if ok:
            print(f"[OK] RDKit available for higher-fidelity auto detection via: {' '.join(cmdline)}")
            rdkit_found = True
            break
        rdkit_errors.append(f"{' '.join(cmdline)}: {error}")
    if not rdkit_found:
        print("[MISSING] RDKit not found in checked Python interpreters.")
        if rdkit_errors:
            print(f"           Last check: {rdkit_errors[-1]}")
        print("           Automatic interaction detection requires RDKit.")
        print("           If RDKit is installed elsewhere, set PYMOL_FIGURE_RDKIT_PYTHON.")

    if pymol_ok and rdkit_found:
        print("\nReady to render. For best label quality, install Pillow in one Python")
        print("and set PYMOL_FIGURE_PYTHON to that interpreter if PyMOL lacks Pillow.")
    else:
        if not pymol_ok:
            print("\nEnvironment is not ready for rendering until PyMOL is available.")
        if not rdkit_found:
            print("\nEnvironment is not ready for auto-detection until RDKit is available.")
        sys.exit(2)


if __name__ == "__main__":
    main()

