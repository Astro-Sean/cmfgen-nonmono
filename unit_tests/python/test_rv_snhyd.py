#!/usr/bin/env python3
"""Unit tests for new_main/subs/rv_sn_model_snhyd.f (VEL_LAW=13).

This is a new file with no reference counterpart.  It reads
SN_HYDRO_DATA, expands the hydro grid (R -> R + V*dt, with the
matching density scaling), builds an R grid uniform in LOG(tau_es)
(TAU_OPT=0) or in a blended tau/log-R coordinate (TAU_OPT=1), packs
NBND_INS=3 points at each boundary, and interpolates V with the
production monotone cubic -- SIGMA comes from the cubic derivative.

Checks:
  * all newly added input guards abort (ND, missing file, missing
    header counts, non-decreasing input R, non-positive V, non-positive
    expansion factor, non-decreasing expanded R, TAU_OPT=1 with no
    positive Kappa).
  * valid runs: endpoints pinned to the (expanded) hydro profile,
    strict grid monotonicity, V from the cubic matching the hydro
    nodes, SIGMA consistent with dlnV/dlnR-1 (homologous -> 0,
    beta-law -> ~BETA-1), velocity-reversal profiles produce
    SIGMA < -1 (the feature VEL_LAW=13 exists to keep).
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
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link,
    run,
)

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "subs/lin_interp.f",
    "subs/mon_interp.f",
    "subs/mon_int_funs_v2.f",
]

NX = 60
ND = 30
OLD_AGE = 39.0


def build():
    objs = compile_objs(
        "snhyd",
        [TREE / h for h in HELPERS]
        + [
            TREE / "new_main" / "subs" / "rv_sn_model_snhyd.f",
            F90DIR / "drv_snhyd.f90",
        ],
    )
    return link("snhyd", objs, moddir=BINDIR / "obj_snhyd")


def write_snhyd(d, r, v, sig, dens, kap, temp=None, age=OLD_AGE, mass_frac=True):
    """Write a synthetic SN_HYDRO_DATA in the production layout."""
    d.mkdir(parents=True, exist_ok=True)
    n = len(r)
    temp = np.full(n, 2.0) if temp is None else temp

    def put_block(f, name, arr):
        f.write(f"{name}\n")
        for i in range(0, n, 4):
            f.write(" ".join(f"{x:20.12e}" for x in arr[i : i + 4]) + "\n")

    with open(d / "SN_HYDRO_DATA", "w") as f:
        f.write(f"Number of data points:    {n}\n")
        f.write("Number of mass fractions:    1\n")
        f.write("Number of isotopes:    0\n")
        f.write(f"Time(days) since explosion:    {age:.7f}\n")
        f.write("Energy Lum [lsun]:    1.0E+08\n")
        f.write("\n")
        put_block(f, "Radius grid (10^10cm)", r)
        put_block(f, "Velocity (km/s)", v)
        put_block(f, "Sigma (dlnV/dlnr-1)", sig)
        put_block(f, "Temperature (10^4 K)", temp)
        put_block(f, "Density (gm/cm^3)", dens)
        put_block(f, "Kappa (cm^2/gm)", kap)
        if mass_frac:
            put_block(f, "HYD mass fraction", np.ones(n))
    return d


def profiles(kind):
    """R strictly decreasing (outer->inner); V > 0."""
    r = np.geomspace(1.0e5, 1.0e3, NX)
    if kind == "hubble":
        v = 9000.0 * r / r[0]  # V propto R -> sigma = 0
    elif kind == "beta":
        v = 3000.0 + 6000.0 * (r / r[0]) ** 1.5
    elif kind == "reversal":
        v = 3000.0 + 6000.0 * r / r[0]
        v = v + 8000.0 * np.exp(-(((r - 0.18 * r[0]) / (0.05 * r[0])) ** 2))
    else:
        raise ValueError(kind)
    # analytic sigma column: dlnV/dlnR - 1 (for the expansion factor)
    dvdlogr = np.gradient(np.log(v), np.log(r))
    sig = dvdlogr - 1.0
    dens = 1.0e-12 * (r / r[0]) ** -3.0
    kap = np.full(NX, 0.34)
    return r, v, sig, dens, kap


def expand(r, v, sig, dens, dt_days):
    t1 = 0.864 * dt_days
    t2 = 1.0 + t1 * v / r * (sig + 1.0)
    r2 = r + t1 * v
    d2 = dens / t2 / (1.0 + t1 * v / r) ** 2
    return r2, d2


def parse(out):
    lines = out.splitlines()
    rows = {}
    rmax = rcore = None
    for i, ln in enumerate(lines):
        t = ln.split()
        if len(t) == 2 and rmax is None:
            rmax, rcore = float(t[0]), float(t[1])
        elif t and t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    arr = np.array([rows[i] for i in sorted(rows)]) if rows else None
    return rmax, rcore, arr


def run_case(exe, d, age, nd=ND, tau_opt=0):
    return run(exe, age, nd, tau_opt, cwd=d, check=False, timeout=30)


# ----------------------------------------------------------------------
# Python replica of the TAU_OPT=1 grid construction in
# RV_SN_MODEL_SNHYD.  Constants must mirror the PARAMETER block of the
# production routine; the replica covers the blended U coordinate,
# structure-aware insertion with the U_GAP_FAC floor, the linear
# U->logR map, and boundary packing.  _eval_v_sigma additionally
# reproduces the Steffen cubic V/SIGMA evaluation.

NBND_INS = 3
ALPHA_R = 1.5
LOG_DEC_CAP = 4.0
RHO_SPIKE_FAC = 1.5
U_GAP_FAC = 5.0e-2


def _sig_grad(r, v):
    return np.gradient(np.log(v), np.log(r)) - 1.0


def double_shock(nx=200):
    """Two fast shells (sigma < -1 dips) plus a detached dense shell and
    an inner density step: multiple insertion features spread over U."""
    r = np.geomspace(1.0e5, 1.0e3, nx)
    v = 4000.0 + 5000.0 * r / r[0]
    v = v + 7000.0 * np.exp(-(((r - 0.31 * r[0]) / (0.02 * r[0])) ** 2))
    v = v + 5000.0 * np.exp(-(((r - 0.60 * r[0]) / (0.025 * r[0])) ** 2))
    dens = 1.0e-12 * (r / r[0]) ** -3.0
    dens = dens * (1.0 + 30.0 * np.exp(-(((r - 0.50 * r[0]) / (0.01 * r[0])) ** 2)))
    dens = dens * np.where(r < 0.62 * r[0], 8.0, 1.0)
    kap = np.full(nx, 0.34)
    return r, v, _sig_grad(r, v), dens, kap


def sharp_shell(nx=400):
    """Single narrow fast shell whose sigma entry and sigma minimum
    land on adjacent uniform cells bracketing a shared cell boundary:
    the pair collapses the U spacing to ~0.04 of a cell without the
    U_GAP_FAC floor."""
    r = np.geomspace(1.0e5, 1.0e3, nx)
    v = 4000.0 + 5000.0 * r / r[0]
    v = v + 8000.0 * np.exp(-(((r - 0.610 * r[0]) / (0.004 * r[0])) ** 2))
    dens = 1.0e-12 * (r / r[0]) ** -3.0
    kap = np.full(nx, 0.34)
    return r, v, _sig_grad(r, v), dens, kap


def _predict_ugrid(r, v, sig, dens, kap, nd=ND, gap_fac=U_GAP_FAC):
    """Replica of the TAU_OPT=1 coordinate path.  Returns the packed
    output R grid and the uniform-U targets actually used."""
    n = len(r)
    kap_eff = np.maximum(kap, 1.0e-6 * kap.max()) * dens * 1.0e10
    tau = np.zeros(n)
    for i in range(1, n):
        tau[i] = tau[i - 1] + 0.5 * (kap_eff[i - 1] + kap_eff[i]) * (r[i - 1] - r[i])
    tau[0] = 0.01 * tau[1]
    ltau = np.log(tau)
    logr = np.log(r)
    tau_lo = max(ltau[1], ltau[-1] - LOG_DEC_CAP)
    u = np.zeros(n)
    for i in range(1, n):
        u[i] = (
            u[i - 1]
            + max(ltau[i], tau_lo)
            - max(ltau[i - 1], tau_lo)
            + ALPHA_R * (logr[i - 1] - logr[i])
        )
    mnd = nd - 2 * NBND_INS
    dlog_u = u[-1] / (mnd - 1)
    nt = dlog_u * np.arange(mnd)
    used = np.zeros(mnd, bool)
    for j in range(1, n - 1):
        ins = (
            (sig[j] < -1.0 and sig[j - 1] >= -1.0)
            or (sig[j] < -1.0 and sig[j] <= sig[j - 1] and sig[j] <= sig[j + 1])
            or (
                dens[j] > RHO_SPIKE_FAC * dens[j - 1]
                and dens[j] > RHO_SPIKE_FAC * dens[j + 1]
            )
        )
        if ins and 0.0 < u[j] < u[-1]:
            i0 = int(np.floor(u[j] / dlog_u + 0.5)) + 1
            i0 = max(2, min(mnd - 1, i0)) - 1
            if gap_fac == 0.0:
                ok = not used[i0] and nt[i0 - 1] < u[j] < nt[i0 + 1]
            else:
                ok = (
                    not used[i0]
                    and u[j] - nt[i0 - 1] > gap_fac * dlog_u
                    and nt[i0 + 1] - u[j] > gap_fac * dlog_u
                )
            if ok:
                nt[i0] = u[j]
                used[i0] = True
    ta = np.exp(np.interp(nt, u, logr))
    ta[0] = r[0]
    ta[-1] = r[-1]
    rout = np.zeros(nd)
    rout[0] = ta[0]
    rout[-1] = ta[-1]
    for i in range(2, mnd):
        rout[i + 2] = ta[i - 1]
    rout[1] = ta[0] - (ta[0] - ta[1]) / 20.0
    rout[2] = ta[0] - (ta[0] - ta[1]) / 8.0
    rout[3] = ta[0] - (ta[0] - ta[1]) / 3.0
    rout[-2] = rout[-1] + (ta[-2] - ta[-1]) / 20.0
    rout[-3] = rout[-1] + (ta[-2] - ta[-1]) / 8.0
    rout[-4] = rout[-1] + (ta[-2] - ta[-1]) / 3.0
    return rout, nt, dlog_u


def _mon_cubic_coefs(y, x):
    """Steffen (1990) monotone cubic, replicating MON_INT_FUNS_V2."""
    n = len(x)
    h = np.diff(x)
    s = np.diff(y) / h
    d = np.zeros(n)
    d[0] = s[0] + (s[0] - s[1]) * h[0] / (h[0] + h[1])
    d[1:-1] = (s[:-1] * h[1:] + s[1:] * h[:-1]) / (h[:-1] + h[1:])
    d[-1] = s[-1] + (s[-1] - s[-2]) * h[-1] / (h[-2] + h[-1])
    d[0] = (np.sign(s[0]) + np.sign(d[0])) * min(abs(s[0]), 0.5 * abs(d[0]))
    for i in range(1, n - 1):
        d[i] = (np.sign(s[i - 1]) + np.sign(s[i])) * min(
            abs(s[i - 1]), abs(s[i]), 0.5 * abs(d[i])
        )
    d[-1] = (np.sign(s[-1]) + np.sign(d[-1])) * min(abs(s[-1]), 0.5 * abs(d[-1]))
    c = np.zeros((n, 4))
    c[:-1, 0] = (d[:-1] + d[1:] - 2.0 * s) / h**2
    c[:-1, 1] = (3.0 * s - 2.0 * d[:-1] - d[1:]) / h
    c[:-1, 2] = d[:-1]
    c[:-1, 3] = y[:-1]
    c[-1, 2] = d[-1]
    c[-1, 3] = y[-1]
    return c


def _eval_v_sigma(rout, r, v):
    """Cubic V and SIGMA=cubic'(R)*R/V-1 on the output grid, replicating
    the production evaluation loop (segment search marches inward)."""
    n = len(r)
    c = _mon_cubic_coefs(v, r)
    vout = np.zeros(len(rout))
    sout = np.zeros(len(rout))
    j = 0
    for i in range(len(rout)):
        while rout[i] < r[j + 1] and j < n - 2:
            j += 1
        t = rout[i] - r[j]
        vout[i] = ((c[j, 0] * t + c[j, 1]) * t + c[j, 2]) * t + c[j, 3]
        sout[i] = (
            ((3.0 * c[j, 0] * t + 2.0 * c[j, 1]) * t + c[j, 2])
            * rout[i]
            / vout[i]
            - 1.0
        )
    return vout, sout


def main():
    rep = Report("RV_SN_MODEL_SNHYD")
    exe = build()

    # ---------------- guard cases -------------------------------
    d = WORKDIR / "snhyd_ok"
    r, v, sig, dens, kap = profiles("beta")
    write_snhyd(d, r, v, sig, dens, kap)

    p = run_case(exe, d, OLD_AGE, nd=8)
    rep.check(
        "ND < 2*NBND_INS+4 aborts",
        "RESULT_OK" not in p.stdout and "ND too small" in p.stdout,
    )

    d = WORKDIR / "snhyd_nofile"
    d.mkdir(parents=True, exist_ok=True)
    p = run_case(exe, d, OLD_AGE)
    rep.check(
        "missing SN_HYDRO_DATA aborts",
        "RESULT_OK" not in p.stdout and "file not found" in p.stdout.lower(),
    )

    d = WORKDIR / "snhyd_noheader"
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "SN_HYDRO_DATA", "w") as f:
        f.write("Radius grid (10^10cm)\n1.0 2.0 3.0\n")
    p = run_case(exe, d, OLD_AGE)
    rep.check("missing NSP/NX/age header aborts", "RESULT_OK" not in p.stdout)

    d = WORKDIR / "snhyd_badR"
    r2 = r.copy()
    r2[10] = r2[9]  # equal pair, not decreasing
    write_snhyd(d, r2, v, sig, dens, kap)
    p = run_case(exe, d, OLD_AGE)
    rep.check(
        "non-decreasing input R aborts",
        "RESULT_OK" not in p.stdout and "not strictly decreasing" in p.stdout,
    )

    d = WORKDIR / "snhyd_badV"
    v2 = v.copy()
    v2[20] = 0.0
    write_snhyd(d, r, v2, sig, dens, kap)
    p = run_case(exe, d, OLD_AGE)
    rep.check(
        "non-positive V aborts",
        "RESULT_OK" not in p.stdout and "must be positive" in p.stdout,
    )

    # expansion factor <= 0: large dt + sigma << -1
    d = WORKDIR / "snhyd_badexp"
    sig2 = sig.copy()
    sig2[15] = -10.0
    write_snhyd(d, r, v, sig2, dens, kap)
    p = run_case(exe, d, OLD_AGE + 200.0)
    rep.check(
        "non-positive expansion factor aborts",
        "RESULT_OK" not in p.stdout and "non-positive" in p.stdout,
    )

    # expanded R ordering: non-monotonic V + big dt, sigma kept >= -1
    d = WORKDIR / "snhyd_expandorder"
    rv, vv, sv, dv, kv = profiles("reversal")
    write_snhyd(d, rv, vv, np.maximum(sv, -0.9), dv, kv)
    p = run_case(exe, d, OLD_AGE + 200.0)
    rep.check(
        "non-decreasing expanded R aborts",
        "RESULT_OK" not in p.stdout and "Expanded" in p.stdout,
    )

    # TAU_OPT=1 with zero kappa column
    d = WORKDIR / "snhyd_zerokap"
    write_snhyd(d, r, v, sig, dens, np.zeros(NX))
    p = run_case(exe, d, OLD_AGE, tau_opt=1)
    rep.check(
        "TAU_OPT=1 all-zero Kappa aborts",
        "RESULT_OK" not in p.stdout and "positive Kappa" in p.stdout,
    )

    # missing mass-fraction section -> EOF
    d = WORKDIR / "snhyd_nofrac"
    write_snhyd(d, r, v, sig, dens, kap, mass_frac=False)
    p = run_case(exe, d, OLD_AGE)
    rep.check("no mass-fraction section aborts (EOF)", "RESULT_OK" not in p.stdout)

    # zero first density -> nonzero check aborts
    d = WORKDIR / "snhyd_zerodens"
    d2 = dens.copy()
    d2[0] = 0.0
    write_snhyd(d, r, v, sig, d2, kap)
    p = run_case(exe, d, OLD_AGE)
    rep.check("zero leading density aborts", "RESULT_OK" not in p.stdout)

    # ---------------- valid runs ---------------------------------
    # homologous V -> SIGMA clamped to 0
    d = WORKDIR / "snhyd_hubble"
    r, v, sig, dens, kap = profiles("hubble")
    write_snhyd(d, r, v, sig, dens, kap)
    p = run_case(exe, d, OLD_AGE)
    ok = "RESULT_OK" in p.stdout
    rep.check("hubble: completes", ok, p.stdout[-200:] if not ok else "")
    if ok:
        rmax, rcore, a = parse(p.stdout)
        rep.check(
            "hubble: RMAX=hydro edge, RCORE=hydro base",
            abs(rmax - r[0]) / r[0] < 1e-10 and abs(rcore - r[-1]) / r[-1] < 1e-10,
        )
        rep.check("hubble: R strictly decreasing", np.all(np.diff(a[:, 0]) < 0))
        rep.check(
            "hubble: V(1),V(ND) match hydro nodes",
            abs(a[0, 1] - v[0]) / v[0] < 1e-10
            and abs(a[-1, 1] - v[-1]) / v[-1] < 1e-10,
        )
        rep.check("hubble: SIGMA == 0 (Hubble clamp)", np.abs(a[:, 2]).max() == 0.0)

    # beta-law V, no expansion
    d = WORKDIR / "snhyd_beta"
    r, v, sig, dens, kap = profiles("beta")
    write_snhyd(d, r, v, sig, dens, kap)
    p = run_case(exe, d, OLD_AGE)
    ok = "RESULT_OK" in p.stdout
    rep.check("beta TAU_OPT=0: completes", ok, p.stdout[-200:] if not ok else "")
    if ok:
        rmax, rcore, a = parse(p.stdout)
        rep.check("beta: R strictly decreasing", np.all(np.diff(a[:, 0]) < 0))
        rep.check(
            "beta: interior V within hydro range",
            a[:, 1].min() >= v.min() * 0.98 and a[:, 1].max() <= v.max() * 1.02,
        )
        # cubic on smooth beta-law points: sigma ~ dlnV/dlnR-1 of beta
        sig_exp = 1.5 * 6000.0 * (a[:, 0] / r[0]) ** 1.5 / a[:, 1] - 1.0
        err = np.abs(a[:, 2] - sig_exp)
        rep.check(
            "beta: SIGMA ~ analytic dlnV/dlnR-1",
            np.median(err) < 0.05 and err.max() < 0.4,
            f"med={np.median(err):.3f} max={err.max():.3f}",
        )

    # beta with expansion: RMAX/RCORE must follow R+V*dt
    d = WORKDIR / "snhyd_beta_exp"
    write_snhyd(d, r, v, sig, dens, kap)
    dt = 10.0
    p = run_case(exe, d, OLD_AGE + dt)
    ok = "RESULT_OK" in p.stdout
    rep.check("beta expanded: completes", ok, p.stdout[-200:] if not ok else "")
    if ok:
        rex, _ = expand(r, v, sig, dens, dt)
        rmax, rcore, a = parse(p.stdout)
        rep.check(
            "beta expanded: RMAX=R1+V1*dt, RCORE=RN+VN*dt",
            abs(rmax - rex[0]) / rex[0] < 1e-10
            and abs(rcore - rex[-1]) / rex[-1] < 1e-10,
        )
        rep.check(
            "beta expanded: endpoints still on hydro V nodes",
            abs(a[0, 1] - v[0]) / v[0] < 1e-10
            and abs(a[-1, 1] - v[-1]) / v[-1] < 1e-10,
        )

    # velocity reversal: must keep SIGMA < -1 region
    sig_rev = {}
    for opt in (0, 1):
        d = WORKDIR / f"snhyd_rev{opt}"
        rv, vv, sv, dv, kv = profiles("reversal")
        write_snhyd(d, rv, vv, sv, dv, kv)
        p = run_case(exe, d, OLD_AGE, tau_opt=opt)
        ok = "RESULT_OK" in p.stdout
        rep.check(
            f"reversal TAU_OPT={opt}: completes", ok, p.stdout[-200:] if not ok else ""
        )
        if ok:
            _, _, a = parse(p.stdout)
            sig_rev[opt] = a
            rep.check(
                f"reversal TAU_OPT={opt}: SIGMA < -1 present",
                a[:, 2].min() < -1.0,
                f"min={a[:, 2].min():.3f}",
            )
            rep.check(
                f"reversal TAU_OPT={opt}: R strictly decreasing",
                np.all(np.diff(a[:, 0]) < 0),
            )

    # ---------------- TAU_OPT=1 grid construction -------------------
    # Double shock + detached dense shell: multiple insertion features
    # (sigma<-1 entries/minima, density spike, density step).  The
    # output grid must equal the U-coordinate replica node for node.
    d = WORKDIR / "snhyd_dblshock"
    r, v, sig, dens, kap = double_shock()
    write_snhyd(d, r, v, sig, dens, kap)
    p = run_case(exe, d, OLD_AGE, tau_opt=1)
    ok = "RESULT_OK" in p.stdout
    rep.check(
        "double-shock TAU_OPT=1: completes", ok, p.stdout[-200:] if not ok else ""
    )
    if ok:
        _, _, a = parse(p.stdout)
        rexp, nt, _ = _predict_ugrid(r, v, sig, dens, kap)
        rep.check(
            "double-shock: R matches U-coordinate replica",
            np.allclose(a[:, 0], rexp, rtol=1e-9, atol=0.0),
            f"max rel err={np.abs(a[:, 0] - rexp).max() / rexp.min():.2e}",
        )
        rep.check(
            "double-shock: R strictly decreasing", np.all(np.diff(a[:, 0]) < 0)
        )
        rep.check(
            "double-shock: >=2 cells with SIGMA < -1 resolved",
            np.sum(a[:, 2] < -1.0) >= 2,
            f"n={np.sum(a[:, 2] < -1.0)} min={a[:, 2].min():.2f}",
        )
        vexp, sexp = _eval_v_sigma(rexp, r, v)
        rep.check(
            "double-shock: V,SIGMA match Steffen-cubic replica",
            np.allclose(a[:, 1], vexp, rtol=1e-9, atol=0.0)
            and np.allclose(a[:, 2], sexp, rtol=1e-8, atol=1e-8),
            f"maxVerr={np.abs(a[:, 1] - vexp).max():.2e} "
            f"maxSerr={np.abs(a[:, 2] - sexp).max():.2e}",
        )

    # Sharp shell tuned so the sigma entry and sigma minimum map to
    # adjacent uniform cells bracketing the shared boundary.  Without
    # the U_GAP_FAC floor the pair collapses the cell to ~0.04 of a
    # uniform cell; with it the second insertion is rejected and the
    # uniform target survives.
    d = WORKDIR / "snhyd_gapfloor"
    r, v, sig, dens, kap = sharp_shell()
    write_snhyd(d, r, v, sig, dens, kap)
    p = run_case(exe, d, OLD_AGE, tau_opt=1)
    ok = "RESULT_OK" in p.stdout
    rep.check(
        "gap-floor TAU_OPT=1: completes", ok, p.stdout[-200:] if not ok else ""
    )
    if ok:
        _, _, a = parse(p.stdout)
        rfloor, _, _ = _predict_ugrid(r, v, sig, dens, kap, gap_fac=U_GAP_FAC)
        rnofloor, _, _ = _predict_ugrid(r, v, sig, dens, kap, gap_fac=0.0)
        rep.check(
            "gap-floor: R matches floored replica",
            np.allclose(a[:, 0], rfloor, rtol=1e-9, atol=0.0),
            f"max rel err={np.abs(a[:, 0] - rfloor).max() / rfloor.min():.2e}",
        )
        rep.check(
            "gap-floor: rejected insertion leaves uniform cell",
            not np.allclose(a[:, 0], rnofloor, rtol=1e-9, atol=0.0),
            "max rel diff vs no-floor="
            f"{np.abs(a[:, 0] - rnofloor).max() / rnofloor.min():.2e}",
        )

    # ---------------- plot ---------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    rv, vv, sv, dv, kv = profiles("reversal")
    ax[0].plot(rv, vv, "k-", lw=1, label="hydro V")
    for opt, mk in ((0, "ro"), (1, "bx")):
        if opt in sig_rev:
            ax[0].plot(
                sig_rev[opt][:, 0],
                sig_rev[opt][:, 1],
                mk,
                ms=4,
                label=f"grid TAU_OPT={opt}",
            )
    ax[0].set_xscale("log")
    ax[0].set_xlabel("R (1e10 cm)")
    ax[0].set_ylabel("V (km/s)")
    ax[0].legend()
    ax[0].set_title("Velocity reversal preserved")
    if 0 in sig_rev:
        ax[1].plot(sig_rev[0][:, 0], sig_rev[0][:, 2], "ro-", ms=4, label="TAU_OPT=0")
    if 1 in sig_rev:
        ax[1].plot(sig_rev[1][:, 0], sig_rev[1][:, 2], "bx-", ms=4, label="TAU_OPT=1")
    ax[1].axhline(-1, c="k", ls="--", lw=0.8)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("R (1e10 cm)")
    ax[1].set_ylabel("SIGMA")
    ax[1].legend()
    ax[1].set_title("SIGMA = R/V dV/dR - 1")
    fig.tight_layout()
    out = PLOTDIR / "rv_snhyd.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
