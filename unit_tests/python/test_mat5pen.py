#!/usr/bin/env python3
"""Test the production MAT5PEN block-Thomas kernel (P5-006/007).

Files under test:
  subs/mat5pen.f - the Schur-complement / right-hand-side update
  used by the TRI branch of CMF_BLKBAND_V3.

Production semantics exercised through the driver:
  DNOTZERO=T: WRKMAT = F.C - E ; A = A + D.WRKMAT
  always:     A = A - B.C

so A_out = A_in + DNOTZERO * D.(F.C - E) - B.C with ENOTZERO
controlling whether the -E term is subtracted inside WRKMAT.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from common import (
    F90DIR,
    PROD_MOD_INC,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link_prod,
    prod_objs,
    run,
)

N, NR = 5, 3


def write_case(name, mats, enotzero, dnotzero):
    a, b, c, d, e, f = mats
    path = WORKDIR / name
    with open(path, "w") as fh:
        fh.write(
            f"{N} {NR} {'T' if enotzero else 'F'} {'T' if dnotzero else 'F'}\n"
        )
        for m in (a, b, c, d, e, f):
            for row in m:
                fh.write(" ".join(f"{x:.16e}" for x in row) + "\n")
    return path


def expected(mats, enotzero, dnotzero):
    a, b, c, d, e, f = mats
    out = a.copy()
    if dnotzero:
        wrk = f @ c
        if enotzero:
            wrk = wrk - e
        out = out + d @ wrk
    return out - b @ c


def parse_a(out):
    lines = out.splitlines()
    i = lines.index("A_OUT")
    return np.array(
        [[float(x) for x in ln.split()] for ln in lines[i + 1 : i + 1 + N]]
    )


def main():
    rep = Report("MAT5PEN block-Thomas kernel")
    drv_obj = compile_objs(
        "mat5pen_drv",
        [F90DIR / "drv_mat5pen.f90"],
        extra_flags=PROD_MOD_INC,
    )
    binary = link_prod(
        "drv_mat5pen",
        prod_objs(
            [
                TREE / "tools" / "set_kind_module.f",
                TREE / "subs" / "mat5pen.f",
            ]
        )
        + drv_obj,
    )

    rng = np.random.default_rng(11)
    mats = [rng.uniform(-1.0, 1.0, (N, NR)) for _ in (0,)]
    a = mats[0]
    b = rng.uniform(-1.0, 1.0, (N, N))
    c = rng.uniform(-1.0, 1.0, (N, NR))
    d = rng.uniform(-1.0, 1.0, (N, N))
    e = rng.uniform(-1.0, 1.0, (N, NR))
    f = rng.uniform(-1.0, 1.0, (N, N))
    mats = (a, b, c, d, e, f)

    for eno, dno in ((True, True), (False, True), (True, False), (False, False)):
        path = write_case(f"mat5_{eno}_{dno}.txt", mats, eno, dno)
        got = parse_a(run(binary, path).stdout)
        exp = expected(mats, eno, dno)
        err = np.abs(got - exp).max()
        rep.check(
            f"ENOTZERO={eno} DNOTZERO={dno}: A == A + D(F.C-E) - B.C",
            err < 1e-12,
            f"maxerr={err:.2e}",
        )

    # TRI-shaped use: m[k] = C - B.alpha and phi' = L - B.phi - A.phi''
    # on a two-block system; the kernel calls match cmf_blkband_v3.f
    # (:661-674).  With D=E=0 and ENOTZERO irrelevant, A_out=A-B.C.
    n2, nr2 = N, N  # vector case NR=N for clarity
    c2 = rng.uniform(-0.5, 0.5, (N, NR))
    b2 = rng.uniform(-0.5, 0.5, (N, N))
    phi_prev = rng.uniform(-0.5, 0.5, (N, NR))
    m2 = (c2.copy(), b2, phi_prev, np.zeros((N, N)), np.zeros((N, NR)), np.zeros((N, N)))
    path = write_case("mat5_mblock.txt", m2, True, False)
    got = parse_a(run(binary, path).stdout)
    exp = c2 - b2 @ phi_prev
    err = np.abs(got - exp).max()
    rep.check(
        "Schur update m[k]=C[k]-B[k].phi[k-1] matches DGEMM",
        err < 1e-12,
        f"maxerr={err:.2e}",
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
