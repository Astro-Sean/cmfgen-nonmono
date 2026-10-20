#!/usr/bin/env python3
"""Unit tests for new_main/subs/rv_sn_model_sniin.f (VEL_LAW=12).

Three behaviors under test:

1. RDINR path (RDINR=T): unchanged beta-law setup --
   V(I) = VCORE*(R(I)/R(ND))**BETA1, SIGMA = BETA1-1.  New and
   reference builds must produce identical grids.

2. Hydro path (RDINR=F), monotone-V table: the fixed routine calls
   LIN_INTERP(R,V,ND,R_HYDRO,V_HYDRO,NHYDRO) so V is interpolated
   from the hydro table onto the tau_es-derived R grid (the hydro
   arrays must first be genuinely reversed -- the old I=1,NHYDRO swap
   loop was a no-op).  Checks: run completes on every vfill, V equals
   the piecewise-linear hydro table evaluated at R, and SIGMA matches
   the centred-difference formula applied to the *returned* V.
   The reference build still aborts in LIN_INTERP (swapped call).

3. Hydro path, repeated-V table (velocity plateau): W values need not
   be strictly monotonic in LIN_INTERP -- only the abscissa must.  A
   table with a flat V segment must interpolate cleanly and produce
   SIGMA = -1 across the plateau (dV/dR = 0).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (
    BINDIR,
    F90DIR,
    PLOTDIR,
    REF,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link,
    run,
)

TREE_HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "subs/lin_interp.f",
]

ND = 20
RMAX, RP, VCORE, BETA1 = 5.0e4, 1.0e3, 10.0, 1.5



def build():
    bins = {}
    for tag, root in (("new", TREE), ("ref", REF)):
        objs = compile_objs(
            f"sniin_{tag}",
            [root / h for h in TREE_HELPERS]
            + [
                root / "new_main" / "subs" / "rv_sn_model_sniin.f",
                F90DIR / "drv_sniin.f90",
            ],
        )
        bins[tag] = link(f"sniin_{tag}", objs, moddir=BINDIR / "obj_sniin_help")
    return bins


def write_rdinr(d):
    d.mkdir(parents=True, exist_ok=True)
    r = np.geomspace(RMAX * 0.98, RP, ND)
    with open(d / "RDINR", "w") as f:
        f.write("!Format date: synthetic unit-test grid\n")
        f.write("1.0 1.0 1 %d\n" % ND)
        for i in range(ND):
            f.write(f"{r[i]:16.8e} 1.0 1.0 1.0\n")
            f.write("0.5\n")


R_HYD = np.linspace(1.0e3, 1.0e5, 8)  # file order: increasing radius
TAU_HYD = np.geomspace(1.0e-4, 30.0, 8)
V_MONO = np.linspace(9000.0, 300.0, len(R_HYD))  # monotone decreasing
V_PLAT = np.array([9000.0, 6000.0, 4000.0, 4000.0, 4000.0, 2500.0, 1200.0, 300.0])


def write_hydro(d, v_row=None):
    """Synthetic input_hydro.dat.

    The production row read is
        READ(LU,*) R_HYDRO(I), V_HYDRO, XD, XT, XL, XTAU, TAU_HYDRO(I)
    where V_HYDRO is the *whole array* (an upstream quirk present in
    both trees), so one READ consumes NHYDRO+6 values.  Giving each
    row exactly 14 columns with NHYDRO=8 makes each read take exactly
    one line: R_HYDRO(i)=col1, V_HYDRO(1:8)=cols2-9 (re-read every
    row, so keep them identical), TAU_HYDRO(i)=col14.
    """
    d.mkdir(parents=True, exist_ok=True)
    if v_row is None:
        v_row = np.linspace(9000.0, 300.0, len(R_HYD))
    n = len(R_HYD)
    with open(d / "input_hydro.dat", "w") as f:
        f.write("F\n")  # ADD_ACC_ZONE
        f.write("header1\n")
        f.write(f"{n}\n")
        f.write("header2\n")
        for i in range(n):
            cols = [R_HYD[i], *v_row, 1.0, 1.0, 1.0, 1.0, TAU_HYD[i]]
            assert len(cols) == 14
            f.write(" ".join(f"{x:.8e}" for x in cols) + "\n")


def parse(out):
    rows = {}
    for ln in out.splitlines():
        t = ln.split()
        if t and t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    return np.array([rows[i] for i in sorted(rows)]) if rows else None


def sigma_expected(r, v):
    """The corrected formula + one-sided boundaries."""
    s = np.zeros(len(r))
    for i in range(1, len(r) - 1):
        t1 = (v[i - 1] - v[i]) / (r[i - 1] - r[i])
        t2 = (v[i] - v[i + 1]) / (r[i] - r[i + 1])
        dvdr = (t1 * (r[i] - r[i + 1]) + t2 * (r[i - 1] - r[i])) / (r[i - 1] - r[i + 1])
        s[i] = r[i] * dvdr / v[i] - 1.0
    s[0] = r[0] * (v[0] - v[1]) / (r[0] - r[1]) / v[0] - 1.0
    s[-1] = r[-1] * (v[-2] - v[-1]) / (r[-2] - r[-1]) / v[-1] - 1.0
    return s


def hydro_v_expected(r):
    """Piecewise-linear hydro V evaluated at the returned R grid."""
    return np.interp(r, R_HYD, V_MONO)


def main():
    rep = Report("RV_SN_MODEL_SNIIN")
    bins = build()

    # -- RDINR path: identical in new and ref ---------------------
    d = WORKDIR / "sniin_rdinr"
    write_rdinr(d)
    got = {}
    for tag in ("new", "ref"):
        r = run(bins[tag], "T", RMAX, RP, VCORE, BETA1, "Z", cwd=d, check=False)
        assert "RESULT_OK" in r.stdout, r.stdout[-2000:]
        got[tag] = parse(r.stdout)
    a, b = got["new"], got["ref"]
    rep.check(
        "RDINR: new == ref grid",
        np.array_equal(a, b),
        f"maxdiff={np.abs(a - b).max():.2e}",
    )
    r_g, v_g, s_g = a[:, 0], a[:, 1], a[:, 2]
    vexp = VCORE * (r_g / r_g[-1]) ** BETA1
    rep.check(
        "RDINR: V == VCORE*(R/RP)^BETA1", np.abs(v_g - vexp).max() / v_g.max() < 1e-10
    )
    rep.check(
        "RDINR: SIGMA == BETA1-1 everywhere", np.abs(s_g - (BETA1 - 1.0)).max() < 1e-12
    )
    rep.check(
        "RDINR: R(1)=RMAX, R(ND)=RP",
        abs(r_g[0] - RMAX) < 1e-8 and abs(r_g[-1] - RP) < 1e-8,
    )

    # -- Hydro path, monotone-V table ------------------------------
    d = WORKDIR / "sniin_hydro"
    write_hydro(d)
    results = {}
    for tag in ("new", "ref"):
        r = run(
            bins[tag], "F", RMAX, RP, VCORE, BETA1, "Z",
            cwd=d, check=False, timeout=30,
        )
        results[tag] = parse(r.stdout) if "RESULT_OK" in r.stdout else None
    rep.check(
        "hydro mono-V: new build completes",
        results["new"] is not None,
    )
    # The reference (pre-fix) build either aborts or emits a V grid
    # that does not reproduce the hydro table.
    bad_ref = results["ref"] is None or not np.allclose(
        results["ref"][:, 1], np.interp(results["ref"][:, 0], R_HYD, V_MONO),
        rtol=1e-6, atol=1e-6,
    )
    rep.check(
        "hydro mono-V: ref build is broken (swapped call)",
        bad_ref,
    )
    a = results.get("new")
    if a is not None:
        r_h, v_h, s_h = a[:, 0], a[:, 1], a[:, 2]
        rep.check(
            "hydro mono-V: R strictly decreasing",
            (np.diff(r_h) < 0).all(),
        )
        v_exp = hydro_v_expected(r_h)
        err = np.abs(v_h - v_exp).max() / v_exp.max()
        rep.check(
            "hydro mono-V: V == hydro table interp at R",
            err < 1e-6,
            f"relerr={err:.2e}",
        )
        s_exp = sigma_expected(r_h, v_h)
        err = np.abs(s_h - s_exp).max()
        rep.check(
            "hydro mono-V: SIGMA == centred diff of returned V",
            err < 1e-9,
            f"maxerr={err:.2e}",
        )
        # vfill is inert now: 'I'/'D' fills must give identical grids.
        r2 = run(
            bins["new"], "F", RMAX, RP, VCORE, BETA1, "I",
            cwd=d, check=False, timeout=30,
        )
        b = parse(r2.stdout) if "RESULT_OK" in r2.stdout else None
        rep.check(
            "hydro mono-V: vfill=I output identical (V is output)",
            b is not None and np.array_equal(a, b),
        )

    # -- Hydro path, repeated-V table (velocity plateau) -----------
    d = WORKDIR / "sniin_hydro_plat"
    write_hydro(d, v_row=V_PLAT)
    r = run(
        bins["new"], "F", RMAX, RP, VCORE, BETA1, "Z", cwd=d, check=False, timeout=30
    )
    b = parse(r.stdout) if "RESULT_OK" in r.stdout else None
    rep.check("hydro rep-V: new build completes", b is not None)
    if b is not None:
        r_p, v_p, s_p = b[:, 0], b[:, 1], b[:, 2]
        v_exp = np.interp(r_p, R_HYD, V_PLAT)
        rep.check(
            "hydro rep-V: V == plateau table interp at R",
            np.abs(v_p - v_exp).max() / v_exp.max() < 1e-6,
        )
        rep.check(
            "hydro rep-V: SIGMA == centred diff of returned V",
            np.abs(s_p - sigma_expected(r_p, v_p)).max() < 1e-9,
        )
        # A cell fully inside the plateau has zero dV/dR -> SIGMA = -1.
        rep.check(
            "hydro rep-V: plateau cells carry SIGMA == -1",
            (np.abs(s_p + 1.0) < 1e-9).any(),
            f"min|sig+1|={np.abs(s_p + 1.0).min():.2e}",
        )

    # -- Plot: RDINR beta-law + hydro interp -----------------------
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog(r_g, v_g, "ro-", ms=3, label="driver output")
    ax[0].loglog(r_g, vexp, "k--", lw=0.8, label="beta law")
    ax[0].set_xlabel("R")
    ax[0].set_ylabel("V")
    ax[0].legend()
    ax[0].set_title("RDINR beta-law velocity")
    if a is not None:
        ax[1].plot(r_h, v_h, "ro-", ms=3, label="interp output")
        ax[1].plot(R_HYD, V_MONO, "ks", ms=5, label="hydro table")
        if b is not None:
            ax[1].plot(r_p, v_p, "bx-", ms=3, label="plateau table")
            ax[1].plot(R_HYD, V_PLAT, "b^", ms=4, mfc="none")
        ax[1].set_xscale("log")
        ax[1].set_xlabel("R")
        ax[1].set_ylabel("V")
        ax[1].legend()
        ax[1].set_title("hydro table interpolation")
    fig.tight_layout()
    out = PLOTDIR / "rv_sniin.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
