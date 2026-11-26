#!/usr/bin/env python3
"""Matrix-free augmented-operator battery (FIX 4, phase 2 stage 2).

drv_mom_augmv builds the same manufactured captured systems as
drv_mom_augres, packs six test vectors on the interleaved layout
z(k) = (x(1),h(1),...,x(ND)), calls the production MOM_AUG_MATVEC,
and dumps NU_*, every Z, and every AZ. This battery rebuilds the
numpy augmented operator from the same dump (the test_mom_himpl
assemble oracle, corrected for the boundary EPSP sign) and checks
the Fortran apply row by row, then the semantics: skipped columns
return zero rows, and on the manufactured-exact cases A z1 = b so
the committed field annihilates the augmented residual.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import SRC, WORKDIR, Report, compile_and_link_prod, run
from test_mom_augres import HELPERS, ND, NCF
from test_mom_himpl import NW, NZ, assemble, zpack, zxi

F90DIR = Path(__file__).resolve().parents[1] / "f90"
CASES = [
    "exact", "epsmix", "lagb", "him", "clip", "bnd",
    "epscap", "mix", "offg", "skip",
]
NV = 6


def build():
    return compile_and_link_prod(
        "mom_augmv",
        [SRC / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_mom_augmv.f"],
    )


def load_mv(tag):
    """Same section grammar as test_mom_augres.load plus ZVn/AZVn."""
    d = {"nd": ND, "ncf": NCF}
    names2r = ["JNEXT", "TA", "TB", "TC", "CUP", "CUPU", "VB", "VC",
               "LAGSCL", "SRC"]
    names2r_m = ["HNEXT", "HU", "HL", "HS", "EPS", "EPSU", "EPSP",
                 "EPSPU"]
    names2i = ["KUP", "KUH"]
    names1r = ["T1B"]
    names1i = ["OUTBC", "HCHK", "NMJ", "NMH"]
    zsecs = {f"ZV{i}": np.zeros(NZ) for i in range(1, NV + 1)}
    asecs = {f"AZV{i}": np.zeros(NZ) for i in range(1, NV + 1)}
    arrs = {}
    sec = None
    with open(WORKDIR / f"{tag}_mv.txt") as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith("#"):
                p = ln.split()
                if p[1] == "SHAPE":
                    assert (int(p[2]), int(p[3])) == (ND, NCF)
                    sec = None
                    continue
                if p[1] == "FLAGS":
                    sec = None
                    continue
                sec = p[1]
                if sec in names2r:
                    arrs[sec] = np.zeros((ND, NCF))
                elif sec in names2r_m:
                    arrs[sec] = np.zeros((ND - 1, NCF))
                elif sec in names2i:
                    arrs[sec] = np.zeros(
                        (ND if sec == "KUP" else ND - 1, NCF), dtype=int
                    )
                elif sec in names1r:
                    arrs[sec] = np.zeros(NCF)
                elif sec in names1i:
                    n = NCF if sec in ("OUTBC", "HCHK") else ND
                    arrs[sec] = np.zeros(n, dtype=int)
                continue
            p = ln.split()
            if sec in names2r + names2r_m:
                arrs[sec][int(p[0]) - 1, int(p[1]) - 1] = float(p[2])
            elif sec in names2i:
                arrs[sec][int(p[0]) - 1, int(p[1]) - 1] = int(p[2])
            elif sec in names1r + names1i:
                arrs[sec][int(p[0]) - 1] = float(p[1])
            elif sec in zsecs:
                zsecs[sec][int(p[0]) - 1] = float(p[1])
            elif sec in asecs:
                asecs[sec][int(p[0]) - 1] = float(p[1])
    for n in names2r + names2r_m + names2i + names1r + names1i:
        d[n.lower()] = arrs[n]
    for i in range(1, NV + 1):
        d[f"zv{i}"] = zsecs[f"ZV{i}"]
        d[f"azv{i}"] = asecs[f"AZV{i}"]
    return d


def main():
    rep = Report("test_mom_augmv")
    exe = build()
    for case in CASES:
        r = run(exe, case, f"mv_{case}", cwd=WORKDIR)
        rep.check(f"{case}: driver ran", "RESULT_OK" in r.stdout)

    for case in CASES:
        d = load_mv(f"mv_{case}")
        _, _, _, a_aug, b_aug = assemble(d)
        # Skipped (uncaptured) columns: the operator is undefined on
        # them; production convention is zero rows.
        skip_cols = [
            k for k in range(NCF) if int(d["outbc"][k]) < 0
        ]
        for k in skip_cols:
            a_aug[k * NW : (k + 1) * NW, :] = 0.0
        for iv in range(1, NV + 1):
            z = d[f"zv{iv}"]
            az_np = a_aug @ z
            az_f = d[f"azv{iv}"]
            num = np.linalg.norm(az_f - az_np)
            den = max(np.linalg.norm(az_np), 1.0e-30)
            rep.check(
                f"{case}: AZ{iv} matches oracle",
                num / den < 1.0e-10,
                f"rel {num / den:.3e} num {num:.3e}",
            )
        for k in skip_cols:
            blk = d["azv1"][k * NW : (k + 1) * NW]
            rep.check(
                f"{case}: skipped col {k + 1} zero rows",
                np.max(np.abs(blk)) == 0.0,
                f"max {np.max(np.abs(blk)):.3e}",
            )
        # Committed field annihilates A z - b where manufacture is
        # exact and no nonlinear guard fired (cap/clip cases differ
        # by construction).
        if case in ("exact", "epsmix", "lagb", "him", "bnd"):
            zc = zpack(d)
            res = np.linalg.norm(a_aug @ zc - b_aug) / max(
                np.linalg.norm(b_aug), 1.0e-30
            )
            rep.check(
                f"{case}: committed field solves A z = b",
                res < 1.0e-10,
                f"res {res:.3e}",
            )
    rep.summary()


if __name__ == "__main__":
    main()
