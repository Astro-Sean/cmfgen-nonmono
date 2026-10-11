#!/usr/bin/env python3
"""Unit tests for new_main/subs/rv_sn_model_sniin.f (VEL_LAW=12).

Three behaviors under test:

1. RDINR path (RDINR=T): unchanged beta-law setup --
   V(I) = VCORE*(R(I)/R(ND))**BETA1, SIGMA = BETA1-1.  New and
   reference builds must produce identical grids.

2. Hydro path (RDINR=F): the routine calls
   LIN_INTERP(V,TA,ND,V_HYDRO,R_HYDRO,NHYDRO), i.e. it uses the
   caller-supplied V array as the *target coordinate* for the velocity
   interpolation.  If V arrives uninitialized (zeros), the new strict
   monotonicity guards abort loudly (pre-existing bug made loud);
   an increasing V mismatches V_HYDRO's direction (the upstream
   reversal loop swaps each pair twice, so V_HYDRO stays decreasing)
   and also aborts.

3. Corrected SIGMA denominator (the actual code change): when V is
   pre-filled strictly increasing inside V_HYDRO's range, the guard
   passes and the routine reaches the SIGMA loop.  The new code
   computes the non-uniform centred difference
       dVdR = [T1*(R_i-R_{i+1}) + T2*(R_{i-1}-R_i)] / (R_{i-1}-R_{i+1})
   while the reference divides by (R_{i-1}-R_i) only.  Both drivers
   run to completion on the same V so the fix is checked cell by cell.
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
VFILL = np.linspace(8900.0, 310.0, ND)  # matches driver 'D'


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


def write_hydro(d, n=8):
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
    r = np.linspace(1.0e3, 1.0e5, n)
    tau = np.geomspace(1.0e-4, 30.0, n)
    v_row = np.linspace(9000.0, 300.0, n)  # decreasing with R index
    with open(d / "input_hydro.dat", "w") as f:
        f.write("F\n")  # ADD_ACC_ZONE
        f.write("header1\n")
        f.write(f"{n}\n")
        f.write("header2\n")
        for i in range(n):
            cols = [r[i], *v_row, 1.0, 1.0, 1.0, 1.0, tau[i]]
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


def sigma_reference(r, v):
    """The pre-fix formula (wrong denominator) for comparison."""
    s = np.zeros(len(r))
    for i in range(1, len(r) - 1):
        t1 = (v[i - 1] - v[i]) / (r[i - 1] - r[i])
        t2 = (v[i] - v[i + 1]) / (r[i] - r[i + 1])
        dvdr = (t1 * (r[i] - r[i + 1]) + t2 * (r[i - 1] - r[i])) / (r[i - 1] - r[i])
        s[i] = r[i] * dvdr / v[i] - 1.0
    s[0] = r[0] * (v[0] - v[1]) / (r[0] - r[1]) / v[0] - 1.0
    s[-1] = r[-1] * (v[-2] - v[-1]) / (r[-2] - r[-1]) / v[-1] - 1.0
    return s


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

    # -- Hydro path, vfill=Z: guard aborts in new build -----------
    d = WORKDIR / "sniin_hydro"
    write_hydro(d)
    r = run(
        bins["new"], "F", RMAX, RP, VCORE, BETA1, "Z", cwd=d, check=False, timeout=30
    )
    aborted = "RESULT_OK" not in r.stdout and (
        "Error in LIN_INTERP" in r.stdout or "monotonic" in r.stdout
    )
    rep.check(
        "hydro vfill=Z: uninit V aborts via guard",
        aborted,
        "guard fired" if aborted else r.stdout[-300:],
    )

    # -- Hydro path, vfill=I: direction mismatch aborts ------------
    r = run(
        bins["new"], "F", RMAX, RP, VCORE, BETA1, "I", cwd=d, check=False, timeout=30
    )
    aborted = "RESULT_OK" not in r.stdout and "Error" in r.stdout
    rep.check("hydro vfill=I: wrong-direction V aborts", aborted)

    # -- Hydro path, vfill=D: reaches SIGMA loop -------------------
    # Both new and ref complete; SIGMA differs by the denominator fix.
    results = {}
    for tag in ("new", "ref"):
        r = run(
            bins[tag], "F", RMAX, RP, VCORE, BETA1, "D", cwd=d, check=False, timeout=30
        )
        if "RESULT_OK" in r.stdout:
            results[tag] = parse(r.stdout)
    rep.check(
        "hydro vfill=D: both builds reach SIGMA loop",
        results.get("new") is not None and results.get("ref") is not None,
        str({k: (v is not None) for k, v in results.items()}),
    )
    if results.get("new") is not None:
        a = results["new"]
        r_h, s_new = a[:, 0], a[:, 2]
        s_exp = sigma_expected(r_h, VFILL)
        err = np.abs(s_new - s_exp).max()
        rep.check(
            "hydro vfill=D: SIGMA == corrected centred diff",
            err < 1e-9,
            f"maxerr={err:.2e}",
        )
        if results.get("ref") is not None:
            s_ref_out = results["ref"][:, 2]
            s_old = sigma_reference(r_h, VFILL)
            err_ref = np.abs(s_ref_out - s_old).max()
            diff = np.abs(s_new - s_ref_out).max()
            rep.check(
                "hydro vfill=D: ref == old-denominator formula",
                err_ref < 1e-9 and diff > 1e-2,
                f"err_ref={err_ref:.2e} new_vs_ref={diff:.2e}",
            )

    # -- Plot: RDINR beta-law + SIGMA formulas ---------------------
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog(r_g, v_g, "ro-", ms=3, label="driver output")
    ax[0].loglog(r_g, vexp, "k--", lw=0.8, label="beta law")
    ax[0].set_xlabel("R")
    ax[0].set_ylabel("V")
    ax[0].legend()
    ax[0].set_title("RDINR beta-law velocity")
    if results.get("new") is not None:
        ax[1].plot(r_h, s_new, "ro-", ms=3, label="new (fixed)")
        ax[1].plot(r_h, s_old, "bx-", ms=3, label="old denominator")
        ax[1].plot(r_h, s_exp, "k--", lw=0.8, label="exact")
        ax[1].set_xscale("log")
        ax[1].set_xlabel("R")
        ax[1].set_ylabel("SIGMA")
        ax[1].legend()
        ax[1].set_title("SIGMA denominator fix")
    fig.tight_layout()
    out = PLOTDIR / "rv_sniin.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
