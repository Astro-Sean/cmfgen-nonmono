#!/usr/bin/env python3
"""Shared build/run helpers for the cmfgen_jfnk unit tests.

Tests link the exact production object files: members of
lib/libsubs.a and lib/libmod.a (the same bytes the production
executables link), or the sibling .o left next to sources in
non-archived dirs (obs/, new_main/).  Small CLI drivers, test
adapters, stubs, and reference-tree sources are still compiled
fresh; reference compiles use the probed production flag set so
new-vs-ref bit-identity comparisons stay apples-to-apples.

Requires a completed `make DEBUG=nodebug all` (build.sh runs it
before the suite).  Production objects and .mod files are
compiler-specific: this mode needs ifort.

Tests run the production Fortran source files from the tree (never
copies) through small CLI drivers, run the driver through
subprocess, and check the output.  Plots go to unit_tests/plots/,
binaries to unit_tests/bin/, scratch files to unit_tests/work/.
"""

import os
import shlex
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
        "-fpp",
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
        "-cpp",
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
    # cwd=md: ifort resolves .mod files from the cwd ahead of -I
    # paths, so compiling inside the module dir keeps stray .mod
    # files in the caller's directory from shadowing production
    # module interfaces.
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=md)
    if r.returncode != 0:
        print(r.stdout[-8000:])
        print(r.stderr[-8000:])
        raise SystemExit(f"compile of {name} failed")
    return binary


def compile_objs(name, sources, moddir=None, extra_flags=None, prod=False):
    """Compile `sources` to .o files inside bin/obj_<name>/ (no link).
    Used when the same driver must link against alternative objects
    (new vs reference source files with identical routine names).
    With prod=True the sources compile with the probed production
    flag set - for the reference side of new-vs-ref comparisons;
    production tree sources should go through prod_obj() instead."""
    objdir = BINDIR / f"obj_{name}"
    objdir.mkdir(parents=True, exist_ok=True)
    md = Path(moddir) if moddir else objdir
    md.mkdir(parents=True, exist_ok=True)
    objs = []
    for i, src in enumerate(sources):
        src = Path(src)
        obj = objdir / (src.stem + ".o")
        # prod=True: prod_comp_flags returns the full probed compile
        # line, compiler token included, so FC must not be prepended.
        prefix = prod_comp_flags(src.suffix) if prod else [FC, *FFLAGS]
        if not obj.exists() or obj.stat().st_mtime <= src.stat().st_mtime:
            cmd = [
                *prefix, *OMP, MODFLAG, str(md),
                *(extra_flags or []), "-c", str(src), "-o", str(obj),
            ]
            # cwd=md: ifort resolves .mod files from the cwd ahead of
            # -I paths; compiling inside the module dir keeps stray
            # .mod files in the caller's directory from shadowing the
            # production module interfaces in mod/.
            r = subprocess.run(
                cmd, capture_output=True, text=True, cwd=md
            )
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


def ref_close(a, b, rtol=1e-13):
    """ULP-level agreement check for new-vs-ref output comparison.

    Returns (ok, max_rel_diff).  Elements below ~1e-290 in magnitude
    are ignored, same policy as max_meaningful_diff.

    The new side links production -O3 archive objects while the
    reference side compiles the frozen sources with the probed
    production flag set.  The two source versions are not identical
    (the new tree carries flag-gated branches), so a few-ULP
    codegen difference is expected and is not a behavioural change.
    rtol=1e-13 sits ~3 decades above machine noise and far below any
    real regression; it replaces the old bit-exact array_equal
    checks, which assumed identical sources on both sides."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return False, np.inf
    m = (np.abs(a) > 1e-290) | (np.abs(b) > 1e-290)
    if not m.any():
        return True, 0.0
    av, bv = a[m], b[m]
    rel = np.abs(av - bv) / np.maximum(np.abs(av), np.abs(bv))
    return bool(rel.max() <= rtol), float(rel.max())


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


# ----------------------------------------------------------------------
# Production objects.
#
# The production build archives every compiled routine into
# lib/libsubs.a, lib/libpgplt.a, lib/libmod.a, and the other lib/*.a
# files, and deletes the sibling .o files in archived dirs; non-
# archived dirs (obs/, new_main/, misc/, ...) keep theirs and link
# the .o directly into their executables.  prod_obj() resolves a
# tree source to the exact object the production link uses: the
# sibling .o when present, else the archive member extracted to
# bin/prod_obj/<archive>/.
#
# Freshness: an object counts as current when its timestamp is not
# older than the source file.  Archive member timestamps are read
# from the ar headers at full precision (`ar tv` only shows minutes,
# and `ar x` does not preserve mtimes).  `make -q` per member is not
# usable here: several Makefiles archive the same member name
# (new_main/mod_subs cross-compiles plane/ sources, plane/ archives
# them again later), so member-vs-member ordering inside the
# dependency graph reports false staleness.  The member-vs-source
# check catches the case that matters - a source edited after the
# last build - while transitive module-dependency staleness remains
# the build system's job (build.sh runs make first).
#
# Duplicate member names across archives resolve in production link
# order (libsubs.a first, as in the new_main/obs link lines); the
# duplicates present today are byte-identical anyway.
#
# ifort is required: production .mod files are compiler-specific.
# ----------------------------------------------------------------------

PROD_LIB = TREE / "lib"
PROD_MOD_DIR = TREE / "mod"
PROD_OBJDIR = BINDIR / "prod_obj"
PROD_MOD_INC = [f"-I{PROD_MOD_DIR}"]

# Archive search order follows the production link lines:
# new_main links -lsubs -lpgplt -lmod, obs links
# -lsubs -lpgplt -ldisp -lmod.  No member is shared between disp and
# mod, so this order matches both exes for every duplicate that
# exists (mod_ray_mom_store.o, opac_mod.o in subs+mod).
_PROD_ARC_ORDER = [
    "libsubs.a",
    "libpgplt.a",
    "libmod.a",
    "libdisp.a",
    "libspec.a",
    "libpgplot_stubs.a",
    "libmy_lapack.a",
    "libmy_blas.a",
]

_arc_members = {}
_prod_flag_cache = {}
_prod_link_cache = None


def _require_ifort():
    if FC != IFORT:
        raise SystemExit(
            "production-object tests require ifort: the production "
            ".o/.mod files are compiler-specific"
        )


def _prod_arcs():
    """Production archives in link-resolution order."""
    arcs = []
    for name in _PROD_ARC_ORDER:
        p = PROD_LIB / name
        if p.exists():
            arcs.append(p)
    for p in sorted(PROD_LIB.glob("lib*.a")):
        if p not in arcs:
            arcs.append(p)
    return arcs


def _ar_header_members(arc):
    """{member_name: epoch_mtime} parsed from ar member headers.

    Handles the '//' long-name string table used for members whose
    names exceed 15 chars (e.g. jfnk_solver_v1.o).  The cache is
    keyed on the archive mtime so a rebuild invalidates it."""
    arc = Path(arc)
    mt = arc.stat().st_mtime
    cached = _arc_members.get(arc)
    if cached is not None and cached[0] == mt:
        return cached[1]
    data = arc.read_bytes()
    if data[:8] != b"!<arch>\n":
        raise SystemExit(f"{arc} is not an ar archive")
    hdrs = []
    table = b""
    off = 8
    while off + 60 <= len(data):
        hdr = data[off : off + 60]
        if hdr[58:60] != b"`\n":
            break
        name = hdr[0:16].decode("ascii", "replace").strip()
        try:
            mtime = int(hdr[16:28].decode("ascii", "replace").strip() or 0)
            size = int(hdr[48:58].decode("ascii", "replace").strip() or 0)
        except ValueError:
            break
        if name == "//":
            table = data[off + 60 : off + 60 + size]
        else:
            hdrs.append((name, mtime))
        off += 60 + size + (size & 1)
    members = {}
    for name, mtime in hdrs:
        if name.startswith("/") and name[1:].isdigit():
            i = int(name[1:])
            end = table.find(b"/", i)
            name = table[i : end if end != -1 else None].decode(
                "ascii", "replace"
            )
        elif name.endswith("/"):
            name = name[:-1]
        if name:
            members[name] = float(mtime)
    _arc_members[arc] = (mt, members)
    return members


def prod_obj(src):
    """Return the production .o for tree source `src` (a Path).

    A sibling .o next to the source wins (obs/, new_main/ and other
    non-archived dirs link their .o files directly).  Otherwise the
    archive member is extracted once into bin/prod_obj/<archive>/ and
    re-extracted when the archive is rebuilt.  Raises SystemExit when
    the object is older than the source - rebuild with
    `make DEBUG=nodebug` or `./build.sh` first.
    """
    _require_ifort()
    src = Path(src)
    src_mt = src.stat().st_mtime
    sib = src.with_suffix(".o")
    if sib.exists() and sib.is_file():
        if sib.stat().st_mtime + 1.0 < src_mt:
            raise SystemExit(
                f"stale object {sib}: older than {src} - run "
                "make DEBUG=nodebug (or ./build.sh) first"
            )
        return sib
    stem = src.stem + ".o"
    arc = None
    member_mt = None
    for a in _prod_arcs():
        members = _ar_header_members(a)
        if stem in members:
            arc, member_mt = a, members[stem]
            break
    if arc is None:
        raise SystemExit(
            f"no production object for {src}: not a sibling .o and "
            f"{stem} is in no lib/*.a archive - is the source "
            "listed in a dir Makefile?"
        )
    if member_mt + 1.0 < src_mt:
        raise SystemExit(
            f"stale archive member {stem} in {arc.name}: older than "
            f"{src} - run make DEBUG=nodebug (or ./build.sh) first"
        )
    outdir = PROD_OBJDIR / arc.stem
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / stem
    if not out.exists() or out.stat().st_mtime < arc.stat().st_mtime:
        r = subprocess.run(["ar", "p", str(arc), stem], capture_output=True)
        if r.returncode != 0 or not r.stdout:
            raise SystemExit(f"ar p {arc} {stem} failed: {r.stderr!r}")
        out.write_bytes(r.stdout)
    return out


def prod_objs(srcs):
    return [prod_obj(s) for s in srcs]


def is_prod_source(path):
    """True when `path` is a production tree source with a .o in the
    production build - i.e. under TREE but outside unit_tests/ and
    tests/.  Drivers, stubs, adapters (unit_tests/f90, tests/f90)
    and reference-tree files return False and must be compiled
    fresh."""
    try:
        rel = Path(path).resolve().relative_to(TREE.resolve())
    except ValueError:
        return False
    return rel.parts[0] not in ("unit_tests", "tests")


def _prod_mod_mtime():
    return max(
        (p.stat().st_mtime for p in PROD_MOD_DIR.glob("*.mod")),
        default=0.0,
    )


def compile_and_link_prod(name, sources, moddir=None, extra_flags=None):
    """Build bin/name from a mixed source list.

    Production tree sources link as their production .o (prod_obj);
    test-only sources (drivers, stubs, adapters under unit_tests/ or
    tests/) and reference-tree files are compiled fresh against the
    production .mod dir.  extra_flags applies to the fresh compiles
    only (e.g. -D defines for a driver) - production objects are
    already built.  Everything links with the production link line
    (link_prod).  Freshly compiled test objects are recompiled when
    mod/ is newer than they are, so a module-interface change cannot
    leave a stale driver object behind.
    """
    prod, test = [], []
    for s in sources:
        s = Path(s)
        (prod if is_prod_source(s) else test).append(s)
    objs = prod_objs(prod) if prod else []
    if test:
        tflags = list(extra_flags or []) + PROD_MOD_INC
        tobjs = compile_objs(
            name + "_src", test, moddir=moddir, extra_flags=tflags
        )
        mod_mt = _prod_mod_mtime()
        stale = [o for o in tobjs if o.stat().st_mtime < mod_mt]
        for o in stale:
            o.unlink()
        if stale:
            tobjs = compile_objs(
                name + "_src", test, moddir=moddir, extra_flags=tflags
            )
        objs += tobjs
    return link_prod(name, objs, moddir=moddir)


_SRC_INDEX = None


def _source_index():
    """{source_stem: [paths]} over the production tree.

    Excludes test trees, the frozen staging copies
    (jfnk_solver_files/, knop_solver_files/ - identical content but
    not compiled by any Makefile), and build outputs.  Used to map
    archive member names back to candidate sources for the
    whole-archive freshness check."""
    global _SRC_INDEX
    if _SRC_INDEX is not None:
        return _SRC_INDEX
    skip_trees = {
        "unit_tests",
        "tests",
        ".git",
        ".devin",
        "cmfgen_obsidian",
        "jfnk_solver_files",
        "knop_solver_files",
    }
    skip_dirs = {"lib", "mod", "exe", "plots", "work", "bin", "plot"}
    idx = {}
    for root, dirs, files in os.walk(TREE):
        dirs[:] = [
            d
            for d in dirs
            if d not in skip_trees and d not in skip_dirs
            and not d.startswith(".")
        ]
        for fn in files:
            if fn.endswith((".f", ".f90", ".F", ".f77", ".for")):
                idx.setdefault(Path(fn).stem, []).append(
                    Path(root) / fn
                )
    _SRC_INDEX = idx
    return idx


def stale_archive_members():
    """[(archive, member, detail)] for archive members older than
    their source.

    Member timestamps come from the ar headers; the comparison source
    is the newest same-stem file in the production tree (conservative
    when several dirs carry the same routine name - the last archive
    writer wins the member slot, and we cannot tell which dir's copy
    produced it).  Members with no matching source (asm/C objects)
    are skipped."""
    idx = _source_index()
    stale = []
    for arc in _prod_arcs():
        for member, mts in _ar_header_members(arc).items():
            stem = member[:-2] if member.endswith(".o") else member
            cands = idx.get(stem)
            if not cands:
                continue
            src_mt = max(p.stat().st_mtime for p in cands)
            if mts + 1.0 < src_mt:
                stale.append((arc.name, member, src_mt - mts))
    return stale


def require_fresh_archives():
    """Fail when any lib/*.a member is older than its source, else
    return silently.  For tests that link whole archives rather than
    individual prod_obj() members."""
    stale = stale_archive_members()
    if stale:
        lines = "\n".join(
            f"  {a}:{m} (source newer by {d:.0f}s)"
            for a, m, d in stale[:20]
        )
        raise SystemExit(
            "stale production archives - run make DEBUG=nodebug "
            f"(or ./build.sh) first:\n{lines}"
        )


def _probe_make_var(var, makefile_dir):
    """Ask a dir Makefile for an expanded variable, e.g. '$(FG)'."""
    probe = WORKDIR / f"probe_{var}.mk"
    probe.write_text(f"probe:\n\t@echo {var}=$({var})\n")
    r = subprocess.run(
        [
            "make", "-s", "-C", str(makefile_dir),
            "-f", "Makefile", "-f", str(probe), "probe",
        ],
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        raise SystemExit(f"make probe for {var} failed: {r.stderr.strip()}")
    for ln in r.stdout.splitlines():
        if ln.startswith(f"{var}="):
            return ln[len(var) + 1 :].strip()
    raise SystemExit(f"make probe for {var} produced no value")


def _strip_mod_path_flags(flagstr):
    """Drop -I/-module tokens (they point at the production mod/ dir,
    which must not leak into reference compiles)."""
    out = []
    skip = False
    for tok in shlex.split(flagstr):
        if skip:
            skip = False
            continue
        if tok == "-module":
            skip = True
            continue
        if tok.startswith("-I"):
            continue
        out.append(tok)
    return out


def prod_comp_flags(suffix):
    """Full production compile command flags for a file of the given
    suffix: the probed $(F90) driver line plus $(FG) for '.f' or
    $(FFREE) for '.f90', with the production -I/-module paths removed
    (callers supply their own module dir).  Probed from the tree
    Makefiles so it tracks the active host block in
    Makefile_definitions - use for reference-tree sources that must
    codegen-match the shipped objects for bit-identity checks."""
    _require_ifort()
    key = "ffree" if suffix.lower() == ".f90" else "fg"
    if key not in _prod_flag_cache:
        var = "FFREE" if key == "ffree" else "FG"
        f90 = shlex.split(_probe_make_var("F90", TREE / "subs"))
        opt = _strip_mod_path_flags(_probe_make_var(var, TREE / "subs"))
        _prod_flag_cache[key] = f90 + opt
    return list(_prod_flag_cache[key])


def prod_link_flags():
    """(compiler_tokens, lapack_tokens) of the production link line,
    probed once from the tree Makefiles."""
    global _prod_link_cache
    _require_ifort()
    if _prod_link_cache is None:
        f90 = shlex.split(_probe_make_var("F90", TREE / "subs"))
        lapack = shlex.split(_probe_make_var("LAPACK", TREE / "subs"))
        blas = shlex.split(_probe_make_var("BLAS", TREE / "subs"))
        _prod_link_cache = (f90, lapack + blas)
    return _prod_link_cache


def link_prod(name, objects, moddir=None):
    """Link objects into bin/name with the production link line
    (ifort -qopenmp -mkl ... plus the MKL/LAPACK/BLAS lists)."""
    _require_ifort()
    binary = BINDIR / name
    objects = [Path(o) for o in objects]
    if binary.exists():
        newest = max(o.stat().st_mtime for o in objects)
        if binary.stat().st_mtime > newest:
            return binary
    f90, libs = prod_link_flags()
    md = Path(moddir) if moddir else MODDIR
    cmd = [*f90, f"-I{md}", *[str(o) for o in objects], *libs, "-o", str(binary)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-6000:])
        print(r.stderr[-6000:])
        raise SystemExit(f"prod link of {name} failed")
    return binary
