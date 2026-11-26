#!/usr/bin/env python3
"""Unit tests for the X < -700 overflow bound in exponx.f (SOB_check WP3).

Finding: EXP(-X) overflows float64 for X <~ -709.78.  A weak line
inversion (GAMH ~ -0.5) at sigma_true in (-1.0005, -1) hits the 1e-4
tau-floor arm, giving GAM ~ -5000 -> EXPONX -> -Inf in the EW block.

Bound: for X < -700 the exponential's argument is pinned at 700 while
the outer X is kept, so the return stays finite with the 1/X shape;
results for |X| <= 700 are bit-identical.  Module counters
N_EXPONX_GUARD / N_DEXPONX_GUARD count activations.

Method: drv_exponx links twice - against the production exponx.o and
against the reference source plus f90/stub_exponx_guard.f90 (the
reference has no module), so identical inputs give a direct
new-vs-reference comparison.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from common import (
    F90DIR,
    PROD_MOD_INC,
    REF,
    SRC,
    WORKDIR,
    Report,
    compile_objs,
    link_prod,
    prod_objs,
    ref_close,
    run,
)

HELPERS = [SRC / "tools" / "set_kind_module.f"]
STUB = F90DIR / "stub_exponx_guard.f90"
DRV = F90DIR / "drv_exponx.f90"


def build():
    """Driver + exponx twice: production object vs reference source."""
    md_help = compile_objs(
        "exponx_help", HELPERS, extra_flags=PROD_MOD_INC
    )
    new_exp = prod_objs([SRC / "subs" / "exponx.f"])
    ref_exp = compile_objs(
        "exponx_ref",
        [REF / "subs" / "exponx.f", STUB],
        prod=True,
        moddir=md_help[0].parent,
        extra_flags=PROD_MOD_INC,
    )
    drv_obj = compile_objs(
        "exponx_drv",
        [DRV],
        moddir=md_help[0].parent,
        extra_flags=PROD_MOD_INC,
    )
    bins = {}
    for tag, exp in (("new", new_exp), ("ref", ref_exp)):
        bins[tag] = link_prod(
            f"exponx_{tag}",
            md_help + drv_obj + exp,
            moddir=md_help[0].parent,
        )
    return bins


def write_case(name, xs):
    path = WORKDIR / name
    with open(path, "w") as f:
        f.write(f"{len(xs)}\n")
        for x in xs:
            f.write(f"{x:.16e}\n")
    return path


def parse_rows(out):
    rows, guard = [], None
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "GUARD":
            guard = (int(t[1]), int(t[2]))
        else:
            rows.append([float(x) for x in t[1:]])
    return np.array(rows), guard


def main():
    rep = Report("EXPONX X < -700 overflow bound (WP3)")
    bins = build()

    # -- T1: |X| <= 700 band is bit-identical to upstream, and matches
    # the analytic formula on both branches.
    xs_band = (
        list(np.linspace(-700.0, -1.0, 40))
        + list(np.geomspace(1e-4, 60.0, 40))
        + [0.0, -1e-4, 1e-4, -0.5, 39.9, 40.0, 41.0, 100.0]
    )
    case = write_case("exponx_band.txt", xs_band)
    a, ga = parse_rows(run(bins["new"], case).stdout)
    b, gb = parse_rows(run(bins["ref"], case).stdout)
    x, f_new, d_new = a[:, 0], a[:, 1], a[:, 2]
    f_ref, d_ref = b[:, 1], b[:, 2]

    s = np.abs(x) < 1e-3
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        f_ana = np.where(
            s,
            1.0 - x * (0.5 - x / 6.0 * (1.0 - x / 4.0)),
            np.where(x < 40.0, (1.0 - np.exp(-x)) / x, 1.0 / x),
        )
        y = np.where(x < 40.0, np.exp(-x), 0.0)
        d_ana = (y - (1.0 - y) / x) / x
    d_ana[s] = -0.5 + x[s] * (1.0 - x[s] * (0.375 - 0.1 * x[s])) / 3.0

    ok_f, rel_f = ref_close(f_new, f_ana)
    ok_d, rel_d = ref_close(d_new, d_ana)
    ok = (
        np.array_equal(f_new, f_ref)
        and np.array_equal(d_new, d_ref)
        and ok_f
        and ok_d
        and ga == (0, 0)
    )
    rep.check(
        "T1 |X|<=700: bit-identical to upstream + analytic, no guard hits",
        ok,
        f"rel_f={rel_f:.2e} rel_d={rel_d:.2e} guard={ga}",
    )

    # -- T2: guarded Xs stay finite where the reference overflows.
    xs_big = [-701.0, -705.0, -708.0, -709.78, -710.0, -800.0,
              -5000.0, -1.0e6, -1.0e10, -1.0e300]
    case = write_case("exponx_big.txt", xs_big)
    a_big, ga = parse_rows(run(bins["new"], case).stdout)
    b, _ = parse_rows(run(bins["ref"], case).stdout)
    ok = np.all(np.isfinite(a_big[:, 1])) and np.all(np.isfinite(a_big[:, 2]))
    rep.check(
        "T2 X<-700: all outputs finite",
        ok,
        f"new f range [{a_big[:,1].min():.2e},{a_big[:,1].max():.2e}] "
        f"ref nonfinite={np.sum(~np.isfinite(b[:,1:]))} cells",
    )

    # -- T3: WP3 scenario mapping - sigma=-1.0001, GAMH=-0.5, radial ray.
    # T1 floor pins T1 at 1e-4 -> GAM = -0.5/1e-4 = -5000 -> EXPONX finite.
    i5 = xs_big.index(-5000.0)
    rep.check(
        "T3 sigma=-1.0001/GAMH=-0.5/radial: GAM=-5000 finite",
        np.isfinite(a_big[i5, 1]) and a_big[i5, 1] > 1e290,
        f"EXPONX(-5000)={a_big[i5,1]:.3e} (ref={b[i5,1]:.3e})",
    )

    # -- T4: activation counters = number of guarded evaluations.
    ok = ga == (len(xs_big), len(xs_big))
    rep.check(
        "T4 counters count guarded evals",
        ok,
        f"GUARD={ga} expected=({len(xs_big)},{len(xs_big)})",
    )

    # -- T5: continuity across the -700 boundary (pinned form matches
    # the regular branch just inside the boundary).
    xs_edge = [-700.0, -700.0 * (1.0 + 1e-12)]
    case = write_case("exponx_edge.txt", xs_edge)
    a, _ = parse_rows(run(bins["new"], case).stdout)
    rel_f = abs(a[1, 1] - a[0, 1]) / abs(a[0, 1])
    rel_d = abs(a[1, 2] - a[0, 2]) / abs(a[0, 2])
    rep.check(
        "T5 continuous at X=-700",
        rel_f < 1e-9 and rel_d < 1e-9,
        f"rel_f={rel_f:.2e} rel_d={rel_d:.2e}",
    )

    # -- T6: guarded return stays huge-positive and keeps 1/X ordering.
    ok = (
        a_big[i5, 1] > 0.0
        and a_big[1, 1] > 0.0
        and abs(a_big[1, 1]) < abs(a_big[0, 1])  # |X| larger -> smaller
    )
    rep.check("T6 bounded return keeps 1/X shape", ok, "")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
