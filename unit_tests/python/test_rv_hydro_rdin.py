#!/usr/bin/env python3
"""Unit tests for SET_RV_HYDRO_MODEL_V3 in its RDINR branch (VEL_LAW=11).

Under RDINR=T the routine reads the R grid (and a V column) from
RDINR, rescales R to the hydro inner edge, then rebuilds V on the
hydro-table MON_INT_FUNS_V2 cubic and sets SIGMA to the cubic's
analytic derivative.

Item-2 fix under test: when a nonmono flag is set, that SIGMA is
rebuilt by SIGMA_TRUE_CORE -- centred differences on the final
(R,V) pair with the same stagnation guard and |sigma| cap as
SIGMA_TRUE, so the work array diagnostics consume is consistent
with the array the nonmono solvers use.

Checks:
  * stock (flags off): new == reference build bit-for-bit -- the
    flag-gated rebuild leaves the stock RDINR path unchanged.
  * KNOP: SIGMA == centred FD of the returned (R,V) pair,
    one-sided at both boundaries.
  * KNOP: the LUER note announcing the rebuild is printed.
  * KNOP: a non-monotonic hydro V table keeps sigma < -1 cells
    (no clamp is applied to SIGMA under KNOP).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from common import (
    BINDIR,
    F90DIR,
    REF,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link,
    run,
)

SHARED = [
    "tools/set_kind_module.f",
    "new_main/mod_subs/control_variable_mod.f",
    "unix/term_io.f",
    "tools/get_lu.f",
    "subs/mon_interp.f",
    "subs/mon_int_funs_v2.f",
    "new_main/subs/adjust_sn_r_grid.f",
    "new_main/subs/recompute_sigma_true.f",
]

NX = 40
ND = 30

R_HYD = np.geomspace(1.0e5, 1.0e3, NX)  # decreasing
_R = np.geomspace(8.0e4, 1.0e3, ND)  # RDINR grid, R(ND)=R_HYDRO(NX)


def v_table(r):
    """Homologous baseline with a compression dip inside.

    V = 250*(R/1e3) times a localized dip; the dip's outer flank has
    dV/dR < 0 so the true sigma drops below -1 there.
    """
    v = 250.0 * (r / 1.0e3)
    dip = 1.0 - 0.55 * np.exp(-((np.log(r / 2.5e4)) / 0.25) ** 2)
    return v * dip


V_HYD = v_table(R_HYD)


def write_inputs(d):
    """SN_HYDRO_DATA (minimal keyword stream) + RDINR."""
    d.mkdir(parents=True, exist_ok=True)
    rho = 1.0e-13 * (R_HYD / R_HYD[-1]) ** -3.0
    sig = np.zeros(NX)
    with open(d / "SN_HYDRO_DATA", "w") as f:
        f.write("! synthetic SN_HYDRO_DATA for the RDINR unit test\n")
        f.write(f"Number of data points: {NX}\n")
        f.write("Number of mass fractions: 1\n")
        f.write("Time(days) since explosion: 39.0\n")
        f.write("Radius grid\n")
        f.write(" ".join(f"{x:.10e}" for x in R_HYD) + "\n")
        f.write("Velocity\n")
        f.write(" ".join(f"{x:.10e}" for x in V_HYD) + "\n")
        f.write("Sigma\n")
        f.write(" ".join(f"{x:.10e}" for x in sig) + "\n")
        f.write("Density\n")
        f.write(" ".join(f"{x:.10e}" for x in rho) + "\n")
        f.write("Kappa\n")
        f.write(" ".join("1.0e-01" for _ in R_HYD) + "\n")
        f.write("Temperature\n")
        f.write(" ".join("1.0e+00" for _ in R_HYD) + "\n")
        f.write("Mass fractions\n")
        f.write("1.0\n")

    # RDINR: format-date line, header record (two reals + NOLD + ND),
    # then per depth a 6-column row (R + 4 dummies + V) followed by a
    # single-value record (NOLD=1).
    with open(d / "RDINR", "w") as f:
        f.write(" 24-FEB-2004                   !Format date\n")
        f.write(f"{_R[-1]:16.8e} 1.0 1 {ND}\n")
        for i in range(ND):
            f.write(f"{_R[i]:16.8e} 1.0 1.0 1.0 1.0 {v_table(_R[i]):12.5e}\n")
            f.write("1.0\n")


def parse(out):
    rows = {}
    for ln in out.splitlines():
        t = ln.split()
        if t and t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    return np.array([rows[i] for i in sorted(rows)]) if rows else None


def sigma_fd(r, v):
    """SIGMA_TRUE_CORE construction: hybrid log/linear quadratic on
    non-uniform 3-point stencils (30-Oct-2026 scheme); mirrors
    sigma_dvdr in test_sigma_true.py."""
    s = np.zeros(len(r))
    for i in range(1, len(r) - 1):
        if (
            v[i - 1] * v[i + 1] > 0.0
            and v[i] * v[i - 1] > 0.0
            and min(abs(v[i - 1]), abs(v[i]), abs(v[i + 1]))
            > 0.1 * max(abs(v[i - 1]), abs(v[i]), abs(v[i + 1]))
        ):
            hp = np.log(r[i - 1] / r[i])
            hm = np.log(r[i] / r[i + 1])
            s[i] = (
                np.log(abs(v[i - 1])) * hm / (hp * (hp + hm))
                + np.log(abs(v[i])) * (hp - hm) / (hp * hm)
                - np.log(abs(v[i + 1])) * hp / (hm * (hp + hm))
            ) - 1.0
        elif v[i] != 0.0:
            hp = r[i - 1] - r[i]
            hm = r[i] - r[i + 1]
            s[i] = (r[i] / v[i]) * (
                v[i - 1] * hm / (hp * (hp + hm))
                + v[i] * (hp - hm) / (hp * hm)
                - v[i + 1] * hp / (hm * (hp + hm))
            ) - 1.0
    if (
        v[0] * v[1] > 0.0
        and v[1] * v[2] > 0.0
        and min(abs(v[0]), abs(v[1]), abs(v[2]))
        > 0.1 * max(abs(v[0]), abs(v[1]), abs(v[2]))
    ):
        hp = np.log(r[0] / r[1])
        hm = np.log(r[1] / r[2])
        s[0] = (
            np.log(abs(v[0])) * (2.0 * hp + hm) / (hp * (hp + hm))
            - np.log(abs(v[1])) * (hp + hm) / (hp * hm)
            + np.log(abs(v[2])) * hp / (hm * (hp + hm))
        ) - 1.0
    elif v[0] != 0.0:
        hp = r[0] - r[1]
        hm = r[1] - r[2]
        s[0] = (r[0] / v[0]) * (
            v[0] * (2.0 * hp + hm) / (hp * (hp + hm))
            - v[1] * (hp + hm) / (hp * hm)
            + v[2] * hp / (hm * (hp + hm))
        ) - 1.0
    if (
        v[-1] * v[-2] > 0.0
        and v[-2] * v[-3] > 0.0
        and min(abs(v[-1]), abs(v[-2]), abs(v[-3]))
        > 0.1 * max(abs(v[-1]), abs(v[-2]), abs(v[-3]))
    ):
        hp = np.log(r[-2] / r[-1])
        hm = np.log(r[-3] / r[-2])
        s[-1] = -(
            np.log(abs(v[-3])) * hp / (hm * (hp + hm))
            - np.log(abs(v[-2])) * (hp + hm) / (hp * hm)
            + np.log(abs(v[-1])) * (2.0 * hp + hm) / (hp * (hp + hm))
        ) - 1.0
    elif v[-1] != 0.0:
        hp = r[-2] - r[-1]
        hm = r[-3] - r[-2]
        s[-1] = (r[-1] / v[-1]) * (-(
            v[-3] * hp / (hm * (hp + hm))
            - v[-2] * (hp + hm) / (hp * hm)
            + v[-1] * (2.0 * hp + hm) / (hp * (hp + hm))
        )) - 1.0
    return np.clip(s, -1.0e6, 1.0e6)


def build():
    # Shared objects (driver + helpers + our control_variable_mod) are
    # compiled once; only set_rv_hydro_model_v3 differs between trees.
    # The per-tag compile must read module files from the common dir.
    common_moddir = BINDIR / "obj_rvhd_common"
    objs_common = compile_objs(
        "rvhd_common",
        [TREE / h for h in SHARED] + [F90DIR / "drv_rv_hydro_rdin.f90"],
    )
    bins = {}
    for tag, root in (("new", TREE), ("ref", REF)):
        objs = compile_objs(
            f"rvhd_{tag}",
            [root / "new_main" / "subs" / "set_rv_hydro_model_v3.f"],
            moddir=common_moddir,
        )
        bins[tag] = link(
            f"rvhd_{tag}", objs_common + objs, moddir=common_moddir
        )
    return bins


def main():
    rep = Report("SET_RV_HYDRO_MODEL_V3/RDINR")
    bins = build()
    d = WORKDIR / "rvhd_rdinr"
    write_inputs(d)

    out = {}
    for tag in ("new", "ref"):
        r = run(bins[tag], "F", ND, cwd=d, check=False, timeout=60)
        rep.check(
            f"stock {tag}: RDINR run completes",
            "RESULT_OK" in r.stdout,
            r.stdout[-300:] if "RESULT_OK" not in r.stdout else "",
        )
        out[tag] = parse(r.stdout)

    if out.get("new") is not None and out.get("ref") is not None:
        rep.check(
            "stock: new == ref (flag-gated rebuild invisible)",
            np.array_equal(out["new"], out["ref"]),
            f"maxdiff={np.abs(out['new'] - out['ref']).max():.2e}",
        )

    r = run(bins["new"], "T", ND, cwd=d, check=False, timeout=60)
    rep.check("knop: RDINR run completes", "RESULT_OK" in r.stdout)
    rep.check(
        "knop: rebuild note on LUER",
        "SIGMA rebuilt" in r.stdout and "nonmono" in r.stdout,
    )
    a = parse(r.stdout)
    if a is not None:
        r_g, v_g, s_g = a[:, 0], a[:, 1], a[:, 2]
        rep.check(
            "knop: R grid strictly decreasing",
            (np.diff(r_g) < 0).all(),
        )
        err = np.abs(s_g - sigma_fd(r_g, v_g)).max()
        rep.check(
            "knop: SIGMA == hybrid FD of returned V",
            err < 1e-8,
            f"maxerr={err:.2e}",
        )
        rep.check(
            "knop: compression dip keeps sigma < -1",
            (s_g < -1.0).any(),
            f"min={s_g.min():.3f}",
        )
        rep.check(
            "knop: V keeps the non-monotonic dip",
            (np.diff(v_g) > 0).any(),
            f"V range {v_g.min():.1f}..{v_g.max():.1f}",
        )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
