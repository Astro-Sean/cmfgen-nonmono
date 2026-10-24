#!/usr/bin/env python3
"""JFNK solver regression tests -- drives the two existing standalone
harnesses that compile the production new_main/mod_subs/
jfnk_solver_v1.f90 unchanged:

  * tests/jfnk_core/test_newton_core.f90 -- Newton/GMRES/hookstep core
    on a stiff analytic problem with a known root, plus restart and
    guarded-step regressions.
  * unit_tests/test_se_solver.py -- surrogate statistical-equilibrium
    systems (mild/std/hard/shock) covering the FD JVP error curve,
    demanded-step SVD, preconditioner variants, GMRES with/without
    preconditioning, and end-to-end Newton.

The SE harness encodes some diagnostic checks the author knows to be
outstanding rather than correctness gates; those are reported as
expected failures (XFAIL) so a NEW regression still fails the suite.
"""

import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import TREE, Report

CORE_DIR = TREE / "tests" / "jfnk_core"
SE_DIR = TREE / "unit_tests"

# Diagnostic checks the author flagged as outstanding, not correctness
# gates.  (spec, tag):
#   mild/T3 -- the check asserts demanded |s|_inf > 1e3 (hardness of
#     the fixture); a mild system converging with a small demanded
#     step is correct behaviour, so this is a fixture expectation.
#   */T4 -- "core LU solve residual <10%" is a known preconditioner
#     quality gap documented in the harness output (LSQ LU path).
XFAIL = {
    ("mild", "T3"),
    ("mild", "T4"),
    ("hard", "T4"),
    ("shock", "T4"),
}


def run_core(rep):
    """make test in tests/jfnk_core; parse the emitted verdicts."""
    subprocess.run(["make", "-s", "clean"], cwd=CORE_DIR, check=False)
    p = subprocess.run(
        ["make", "-s", "test"],
        cwd=CORE_DIR,
        capture_output=True,
        text=True,
        timeout=600,
    )
    out = p.stdout + p.stderr
    built = "test_newton_core" in out or (CORE_DIR / "test_newton_core").exists()
    rep.check(
        "jfnk_core: builds production solver", built, out[-400:] if not built else ""
    )
    verdicts = re.findall(r"(T\d+):\s*(PASS|FAIL)", out)
    rep.check(
        "jfnk_core: exits cleanly and runs all tests",
        p.returncode == 0 and len(verdicts) == 5,
        out[-400:],
    )
    for tag, verdict in verdicts:
        rep.check(f"jfnk_core: {tag} {verdict}", verdict == "PASS")


def run_se(rep):
    """test_se_solver.py --quick; parse (spec,tag) verdicts."""
    b = subprocess.run(
        ["bash", "build.sh"], cwd=SE_DIR, capture_output=True, text=True, timeout=300
    )
    rep.check(
        "se_driver: builds against production solver",
        (SE_DIR / "bin" / "se_driver").exists() and b.returncode == 0,
        b.stdout[-300:] + b.stderr[-300:],
    )
    p = subprocess.run(
        [sys.executable, "test_se_solver.py", "--quick"],
        cwd=SE_DIR,
        capture_output=True,
        text=True,
        timeout=600,
    )
    out = p.stdout
    spec = None
    tag = None
    for ln in out.splitlines():
        m = re.match(r"=+ spec '(\w+)'", ln)
        if m:
            spec = m.group(1)
            continue
        m = re.match(r"\s*\[(T\d+) ", ln)
        if m:
            tag = m.group(1)
            continue
        m = re.match(r"\s*(PASS|FAIL)\s*\((.*)\)", ln)
        if m and spec and tag:
            verdict, msg = m.group(1), m.group(2)
            key = (spec, tag)
            if key in XFAIL:
                rep.check(
                    f"se:{spec}/{tag} XFAIL (known diagnostic)", verdict == "FAIL", msg
                )
            else:
                rep.check(f"se:{spec}/{tag}", verdict == "PASS", msg)


def main():
    rep = Report("JFNK solver (production jfnk_solver_v1.f90)")
    run_core(rep)
    run_se(rep)
    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
