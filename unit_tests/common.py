#!/usr/bin/env python3
"""Shared build/run helpers for the cmfgen_jfnk unit tests.

Every test compiles the exact production Fortran source files from the
tree (never copies) together with a small CLI driver, runs the driver
through subprocess, and checks the output.  Plots go to
unit_tests/plots/, binaries to unit_tests/bin/, scratch files to
unit_tests/work/.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

UT = Path(__file__).resolve().parent  # unit_tests/
TREE = UT.parent  # cmfgen_jfnk/
REF = TREE.parent / "cmfgen_reference" / "cur_cmf"
BINDIR = UT / "bin"
PLOTDIR = UT / "plots"
WORKDIR = UT / "work"
F90DIR = UT / "f90"

for d in (BINDIR, PLOTDIR, WORKDIR):
    d.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# Compiler selection.  ifort is preferred (production compiler); gfortran
# is the fallback.  Fixed-form .f files need no special flag for either.
# ----------------------------------------------------------------------
IFORT = shutil.which("ifort")
GFORT = shutil.which("gfortran")

if IFORT:
    FC = IFORT
    FFLAGS = [
        "-O0",
        "-g",
        "-traceback",
        "-extend-source",
        "-diag-disable=10448,10412,10434,5462,5199,7712,7346",
        "-warn",
        "none",
    ]
    OMP = ["-qopenmp"]
    MODFLAG = "-module"
    LAPACK = ["-mkl=sequential"]
elif GFORT:
    FC = GFORT
    FFLAGS = [
        "-O0",
        "-g",
        "-fbacktrace",
        "-ffixed-line-length-132",
        "-w",
        "-fallow-argument-mismatch",
        "-std=legacy",
    ]
    OMP = ["-fopenmp"]
    MODFLAG = "-J"
    LAPACK = ["-llapack", "-lblas"]
else:
    sys.exit("unit_tests: no Fortran compiler found (need ifort or gfortran)")

MODDIR = BINDIR / "mod"
MODDIR.mkdir(exist_ok=True)


def compile_driver(name, sources, extra_flags=None, moddir=None):
    """Compile `sources` (paths, exact tree files + driver) into
    bin/name.  Recompiles when any source is newer than the binary.
    Returns the binary path."""
    binary = BINDIR / name
    sources = [Path(s) for s in sources]
    if binary.exists():
        newest = max(p.stat().st_mtime for p in sources)
        if binary.stat().st_mtime > newest:
            return binary
    md = Path(moddir) if moddir else MODDIR
    md.mkdir(parents=True, exist_ok=True)
    cmd = [FC, *FFLAGS, *OMP, MODFLAG, str(md)]
    if extra_flags:
        cmd += extra_flags
    cmd += [str(p) for p in sources] + LAPACK + ["-o", str(binary)]
    print(f"  compile {name}: {' '.join(str(c) for c in cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-8000:])
        print(r.stderr[-8000:])
        raise SystemExit(f"compile of {name} failed")
    return binary


def compile_objs(name, sources, moddir=None):
    """Compile `sources` to .o files inside bin/obj_<name>/ (no link).
    Used when the same driver must link against alternative objects
    (new vs reference source files with identical routine names)."""
    objdir = BINDIR / f"obj_{name}"
    objdir.mkdir(parents=True, exist_ok=True)
    md = Path(moddir) if moddir else objdir
    md.mkdir(parents=True, exist_ok=True)
    objs = []
    for i, src in enumerate(sources):
        src = Path(src)
        obj = objdir / (src.stem + ".o")
        if not obj.exists() or obj.stat().st_mtime <= src.stat().st_mtime:
            cmd = [FC, *FFLAGS, *OMP, MODFLAG, str(md), "-c", str(src), "-o", str(obj)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print(f"  compile {src} failed:")
                print(r.stdout[-6000:])
                print(r.stderr[-6000:])
                raise SystemExit(f"compile of {src} failed")
        objs.append(obj)
    return objs


def link(name, objects, moddir=None):
    """Link .o files into bin/name."""
    binary = BINDIR / name
    objects = [Path(o) for o in objects]
    if binary.exists():
        newest = max(o.stat().st_mtime for o in objects)
        if binary.stat().st_mtime > newest:
            return binary
    md = Path(moddir) if moddir else MODDIR
    cmd = [FC, *OMP, f"-I{md}", *objects, *LAPACK, "-o", str(binary)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-6000:])
        print(r.stderr[-6000:])
        raise SystemExit(f"link of {name} failed")
    return binary


def run(binary, *args, cwd=None, check=True, timeout=300):
    """Run the driver.  Returns CompletedProcess.  With check=False the
    caller inspects returncode (used for expected-STOP tests)."""
    r = subprocess.run(
        [str(binary)] + [str(a) for a in args],
        capture_output=True,
        text=True,
        cwd=cwd,
        timeout=timeout,
    )
    if check and r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit(f"{binary} {args} failed rc={r.returncode}")
    return r


def max_meaningful_diff(a, b):
    """Largest |a-b| ignoring denormal garbage from uninit memory.

    Some production routines leave an occasional output cell
    uninitialized; it prints as a denormal that jitters run to run.
    Anything below ~1e-290 cannot carry physical information."""
    m = (np.abs(a) > 1e-290) | (np.abs(b) > 1e-290)
    return np.abs(a[m] - b[m]).max() if m.any() else 0.0


def read_vec(path):
    return np.loadtxt(path, dtype=float)


def write_vec(path, arr):
    np.savetxt(path, np.asarray(arr, dtype=float), fmt="%.16e")


class Report:
    """Collect PASS/FAIL lines; exit status = number of failures."""

    def __init__(self, title):
        self.title = title
        self.rows = []
        print(f"\n{'=' * 60}\n{title}\n{'=' * 60}")

    def check(self, name, ok, detail=""):
        tag = "PASS" if ok else "FAIL"
        self.rows.append((name, bool(ok)))
        print(f"  [{tag}] {name}" + (f"  ({detail})" if detail else ""))

    def summary(self):
        nfail = sum(1 for _, ok in self.rows if not ok)
        print(f"{'-' * 60}")
        print(f"{self.title}: {len(self.rows) - nfail} passed, " f"{nfail} failed")
        return nfail
