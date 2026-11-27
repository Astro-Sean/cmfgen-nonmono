#!/usr/bin/env python3
"""VEL_LAW=14 grid diagnostics on the real m15 SN_HYDRO_DATA.

Runs the production RV_SN_MODEL_SNHYD driver (unit_tests/bin/snhyd)
on the staged m15 2.82d input at production ND=70 / TAU_OPT=1, then
measures the grid-quality quantities demanded by the velaw14.md
checklist: tau-cell statistics, |1+SIGMA| zero-gradient census,
feature budget, extremum fidelity, U->R inversion error, quadrature
error, ltau-shift and homology-rescaling invariance.

Writes docs/audits/vel14_grid_diagnostics.json.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from common import BINDIR, TREE, WORKDIR, run
from test_rv_snhyd import (
    NBND_INS,
    ALPHA_R,
    LOG_DEC_CAP,
    RHO_SPIKE_FAC,
    U_GAP_FAC,
    _mon_cubic_coefs,
    _predict_ugrid,
    _eval_v_sigma,
    write_snhyd,
)

SNHYD_EXE = BINDIR / "snhyd"
M15 = (
    TREE.parents[0]
    / "cmfgen_validation"
    / "m15_2p82d_VAL_sob_diag_ngoff"
    / "cmfgen"
    / "2.82days"
    / "SN_HYDRO_DATA"
)
OUT_JSON = TREE / "docs" / "audits" / "vel14_grid_diagnostics.json"
SN_AGE = 2.8240800  # VADAT [SN_AGE]; equals the file's Time column
ND = 70


def read_sn_hydro(path):
    """Parse the production SN_HYDRO_DATA column blocks."""
    txt = Path(path).read_text().splitlines()
    nx = n = None
    blocks = {}
    i = 0
    while i < len(txt):
        ln = txt[i]
        if "Number of data points:" in ln:
            nx = int(ln.split(":")[1])
        if "Number of" in ln:
            i += 1
            continue
        if "Radius grid" in ln or (
            any(
                k in ln
                for k in (
                    "Velocity",
                    "Sigma",
                    "Temperature",
                    "Density",
                    "Kappa",
                    "mass fraction",
                )
            )
        ):
            name = ln.strip()
            vals = []
            i += 1
            while i < len(txt) and len(vals) < nx:
                vals.extend(float(x) for x in txt[i].split())
                i += 1
            blocks[name] = np.array(vals[:nx])
            continue
        i += 1
    return nx, blocks


def run_driver(d, age=SN_AGE, nd=ND, tau_opt=1):
    p = run(SNHYD_EXE, age, nd, tau_opt, cwd=d, check=False, timeout=60)
    rows = {}
    rmax = rcore = None
    for ln in p.stdout.splitlines():
        t = ln.split()
        if len(t) == 2 and rmax is None:
            rmax, rcore = float(t[0]), float(t[1])
        elif t and t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    a = np.array([rows[i] for i in sorted(rows)]) if rows else None
    return p.stdout, rmax, rcore, a


def u_coordinate(ltau, logr):
    tau_lo = max(ltau[1], ltau[-1] - LOG_DEC_CAP)
    u = np.zeros(len(ltau))
    for i in range(1, len(ltau)):
        u[i] = (
            u[i - 1]
            + max(ltau[i], tau_lo)
            - max(ltau[i - 1], tau_lo)
            + ALPHA_R * (logr[i - 1] - logr[i])
        )
    return u, tau_lo


def tau_trap(krho, r):
    t = np.zeros(len(r))
    for i in range(1, len(r)):
        t[i] = t[i - 1] + 0.5 * (krho[i - 1] + krho[i]) * (r[i - 1] - r[i])
    t[0] = 0.01 * t[1]
    return t


def tau_ref(krho, r, refine=20):
    """High-order tau reference: monotone-cubic interpolate krho onto
    a refine x finer grid, trapezoid there, then subsample back to
    the hydro nodes.  Converges the quadrature without relying on a
    non-uniform Simpson rule."""
    from scipy.interpolate import PchipInterpolator

    f = PchipInterpolator(r[::-1], krho[::-1])
    # rf descends outer -> inner so tau accumulates inward like the
    # production trapezoid
    rf = np.geomspace(r[0], r[-1], (len(r) - 1) * refine + 1)
    kf = f(rf)
    tf = np.zeros(len(rf))
    for i in range(1, len(rf)):
        tf[i] = tf[i - 1] + 0.5 * (kf[i - 1] + kf[i]) * (rf[i - 1] - rf[i])
    out = np.interp(r, rf[::-1], tf[::-1])
    out[0] = 0.01 * out[1]
    return out


def feature_census(sig, dens):
    n = len(sig)
    entries, minima, spikes = [], [], []
    for j in range(1, n - 1):
        if sig[j] < -1.0 and sig[j - 1] >= -1.0:
            entries.append(j)
        if sig[j] < -1.0 and sig[j] <= sig[j - 1] and sig[j] <= sig[j + 1]:
            minima.append(j)
        if dens[j] > RHO_SPIKE_FAC * dens[j - 1] and dens[j] > RHO_SPIKE_FAC * dens[j + 1]:
            spikes.append(j)
    return entries, minima, spikes


def main():
    d = WORKDIR / "vel14_m15"
    d.mkdir(parents=True, exist_ok=True)
    nx, b = read_sn_hydro(M15)
    r = b["Radius grid (10^10cm)"]
    v = b["Velocity (km/s)"]
    sig = b["Sigma (dlnV/dlnr-1)"]
    dens = b["Density (gm/cm^3)"]
    kap = b["Kappa (cm^2/gm)"]

    # production driver run on the staged file (SN_AGE == file age:
    # the expansion is a no-op)
    import shutil

    shutil.copy(M15, d / "SN_HYDRO_DATA")
    out, rmax, rcore, a = run_driver(d)
    rep = {"driver_rc_ok": "RESULT_OK" in out, "rmax": rmax, "rcore": rcore}
    if a is None:
        print(out)
        return 1
    rg, vg, sg = a[:, 0], a[:, 1], a[:, 2]

    # U coordinate on the hydro columns (expansion is identity here)
    kap_eff = np.maximum(kap, 1.0e-6 * kap.max()) * dens * 1.0e10
    tau = tau_trap(kap_eff, r)
    ltau = np.log(tau)
    logr = np.log(r)
    u, tau_lo = u_coordinate(ltau, logr)
    rep["nx_file"] = int(nx)
    rep["nd"] = int(len(rg))
    rep["tau_lo_ln"] = float(tau_lo)
    rep["tau_max"] = float(tau[-1])
    rep["tau_outer_file"] = float(tau[0])
    rep["u_range"] = float(u[-1])

    # C9 grid-quality block
    krho_new = np.exp(
        np.interp(np.log(rg), logr[::-1], np.log(kap_eff)[::-1])
    )
    tau_new = np.zeros(len(rg))
    for i in range(1, len(rg)):
        tau_new[i] = tau_new[i - 1] + 0.5 * (krho_new[i - 1] + krho_new[i]) * (
            rg[i - 1] - rg[i]
        )
    tau_new[0] = 0.01 * tau_new[1]
    # cell 1 inherits the artificial TAU(1)=0.01*TAU(2) seed; the
    # spacing stats start at cell 2
    dlt = np.diff(np.log(tau_new[1:]))
    dlr = np.diff(np.log(rg))  # negative
    rep["grid"] = {
        "min_dlog_tau": float(np.min(dlt)),
        "max_dlog_tau": float(np.max(dlt)),
        "ratio_dlog_tau": float(np.max(dlt) / np.min(dlt)),
        "min_abs_dlogR": float(np.min(np.abs(dlr))),
        "n_abs1pSigma_lt_1e-6": int(np.sum(np.abs(1.0 + sg) < 1e-6)),
        "max_abs_sigma": float(np.abs(sg).max()),
        "n_sigma_lt_m1": int(np.sum(sg < -1.0)),
        "tau_inner": float(tau_new[-1]),
        "tau_outer": float(tau_new[0]),
        "v_outer": float(vg[0]),
        "v_inner": float(vg[-1]),
        "n_below_v100": int(np.sum(vg < 100.0)),
        "n_above_v1000": int(np.sum(vg > 1000.0)),
    }

    # B3 zero-gradient set on the production grid
    zg = np.where(np.abs(1.0 + sg) < 1e-3)[0]
    rep["zero_gradient_nodes_1e-3"] = [
        {
            "i": int(i + 1),
            "r": float(rg[i]),
            "v": float(vg[i]),
            "sigma": float(sg[i]),
            "abs1ps": float(abs(1.0 + sg[i])),
        }
        for i in zg
    ]

    # file-side extremum fidelity (B4): peak V in file vs the
    # Steffen-interpolated V on the new grid
    jpk = int(np.argmax(v))
    rep["extremum"] = {
        "file_peak_v": float(v[jpk]),
        "file_peak_r": float(r[jpk]),
        "grid_peak_v": float(vg.max()),
        "grid_peak_r": float(rg[int(np.argmax(vg))]),
        "v_err_frac": float(abs(vg.max() - v[jpk]) / v[jpk]),
        "r_err_frac": float(abs(rg[int(np.argmax(vg))] - r[jpk]) / r[jpk]),
        "note": "m15 peak sits at the inner file boundary; the "
        "interior-extremum check is the synthetic bump below",
    }

    # B4 interior extremum: log-Gaussian velocity bump on a
    # homologous background; sigma column is the analytic
    # dlnV/dlnR-1
    xr = np.geomspace(1.0e5, 1.0e3, 80)
    lgc = np.log(xr / 1.0e4)
    bump = 2000.0 + 6000.0 * np.exp(-0.5 * (lgc / 0.25) ** 2)
    dvb = (
        -6000.0
        * np.exp(-0.5 * (lgc / 0.25) ** 2)
        * lgc
        / 0.25**2
    )  # dV/dlnR
    sig_b = (xr / bump) * dvb - 1.0
    write_snhyd(
        WORKDIR / "vel14_bump",
        xr,
        bump,
        sig_b,
        1.0e-12 * (xr / 1.0e5) ** -3.0,
        np.full(len(xr), 0.34),
        age=SN_AGE,
    )
    _, _, _, ab = run_driver(WORKDIR / "vel14_bump")
    jbg = int(np.argmax(ab[:, 1]))
    rep["extremum_bump"] = {
        "file_peak_v": float(bump.max()),
        "file_peak_r": float(xr[int(np.argmax(bump))]),
        "grid_peak_v": float(ab[:, 1].max()),
        "grid_peak_r": float(ab[jbg, 0]),
        "v_err_frac": float(abs(ab[:, 1].max() - bump.max()) / bump.max()),
        "r_err_frac": float(abs(ab[jbg, 0] - xr[int(np.argmax(bump))]) / xr[int(np.argmax(bump))]),
        "reversals_preserved": int(np.sum(ab[:, 2] < -1.0)),
    }

    # C4 feature budget on the real file
    ent, mn, spk = feature_census(sig, dens)
    rep["features"] = {
        "reversal_entries": len(ent),
        "sigma_minima": len(mn),
        "density_spikes": len(spk),
        "entry_idx": [int(j + 1) for j in ent],
        "minima_idx": [int(j + 1) for j in mn],
        "spike_idx": [int(j + 1) for j in spk],
    }

    # reversals preserved on output: cells where the cubic sigma < -1
    rep["reversals_preserved"] = int(np.sum(sg < -1.0))

    # C5 U->R inversion error: rebuild U(R) and invert the uniform-U
    # targets by Brent on the piecewise-linear U(R) map
    mnd = ND - 2 * NBND_INS
    dlog_u = u[-1] / (mnd - 1)
    targets = dlog_u * np.arange(mnd)
    # forward map U(R) is monotone piecewise-linear in index; invert
    # by bisection on the tabulated function
    from scipy.optimize import brentq

    def u_of_r(rr):
        return float(np.interp(rr, r[::-1], u[::-1]))

    inv_err = []
    for ut in targets:
        try:
            rb = brentq(lambda x: u_of_r(x) - ut, r[-1], r[0], xtol=1e-14)
        except ValueError:
            rb = r[0] if ut <= 0 else r[-1]
        rl = np.exp(np.interp(ut, u, logr))  # production linear interp
        inv_err.append(abs(rl - rb) / rb)
    rep["u_inversion_max_rel_err"] = float(np.max(inv_err))

    # C3 quadrature error: trapezoid tau vs fine-grid reference
    tsim = tau_ref(kap_eff, r)
    rep["tau_quad"] = {
        "max_rel_trap_vs_ref": float(
            np.max(np.abs(tau[1:] - tsim[1:]) / tsim[1:])
        ),
        "med_rel_trap_vs_ref": float(
            np.median(np.abs(tau[1:] - tsim[1:]) / tsim[1:])
        ),
    }

    # B5: max |dSIGMA/dlnR| jump at inserted nodes on the output grid
    dsl = np.diff(sg) / np.diff(np.log(rg))
    jump = np.abs(np.diff(dsl))
    rep["sigma_c2_max_jump"] = float(jump.max())

    # A4 power-law sigma error at production ND on a power-law file
    write_snhyd(
        WORKDIR / "vel14_plaw",
        np.geomspace(1.0e5, 1.0e3, nx),
        5000.0 * (np.geomspace(1.0e5, 1.0e3, nx) / 1.0e5) ** 1.5,
        np.full(nx, 0.5),  # sigma column: b-1 for b=1.5
        1.0e-12 * (np.geomspace(1.0e5, 1.0e3, nx) / 1.0e5) ** -3.0,
        np.full(nx, 0.34),
        age=SN_AGE,
    )
    _, _, _, ap = run_driver(WORKDIR / "vel14_plaw")
    rep["power_law_b1.5"] = {
        "max_abs_sigma_err": float(np.max(np.abs(ap[:, 2] - 0.5))),
        "med_abs_sigma_err": float(np.median(np.abs(ap[:, 2] - 0.5))),
        "max_rel_v_err": float(
            np.max(
                np.abs(ap[:, 1] - 5000.0 * (ap[:, 0] / 1.0e5) ** 1.5)
                / (5000.0 * (ap[:, 0] / 1.0e5) ** 1.5)
            )
        ),
    }

    # B7 convergence order: monotone analytic profile
    # V = V0*(R/R0)^b*(1+a*sin(k*ln(R/R0))), a*k<b so dV/dR>0.
    # File resolution scales with ND (the Steffen interpolant is
    # fixed by the input knots, so ND alone cannot refine it).
    # sigma_analytic = b - 1 + a*k*cos(k*ln(R/R0))/(1+a*sin(...)).
    b_, a_, k_ = 0.8, 0.1, 3.0
    conv = {}
    prev = None
    for n in (70, 140, 280, 560):
        xr = np.geomspace(1.0e5, 1.0e3, n)
        lr = np.log(xr / 1.0e4)
        xv = 2000.0 * (xr / 1.0e4) ** b_ * (1.0 + a_ * np.sin(k_ * lr))
        xs = (
            b_
            + a_ * k_ * np.cos(k_ * lr) / (1.0 + a_ * np.sin(k_ * lr))
            - 1.0
        )
        write_snhyd(
            WORKDIR / f"vel14_conv{n}",
            xr,
            xv,
            xs,
            1.0e-12 * (xr / 1.0e5) ** -3.0,
            np.full(n, 0.34),
            age=SN_AGE,
        )
        _, _, _, ac = run_driver(WORKDIR / f"vel14_conv{n}", nd=n)
        rr = ac[:, 0]
        v_an = 2000.0 * (rr / 1.0e4) ** b_ * (
            1.0 + a_ * np.sin(k_ * np.log(rr / 1.0e4))
        )
        s_an = (
            b_
            + a_
            * k_
            * np.cos(k_ * np.log(rr / 1.0e4))
            / (1.0 + a_ * np.sin(k_ * np.log(rr / 1.0e4)))
            - 1.0
        )
        ev = float(np.max(np.abs(ac[:, 1] - v_an) / v_an))
        es = float(np.max(np.abs(ac[:, 2] - s_an)))
        ent = {"max_rel_v_err": ev, "max_abs_sigma_err": es}
        if prev:
            ent["order_v"] = float(np.log(prev[1] / ev) / np.log(2.0))
            ent["order_sigma"] = float(np.log(prev[2] / es) / np.log(2.0))
        prev = (n, ev, es)
        conv[str(n)] = ent
    rep["convergence_b7"] = conv

    # B6 production Steffen vs scipy PCHIP on a monotone segment:
    # evaluate both on the conv280 output R nodes.  The limiters
    # differ intentionally (Steffen 1990 vs Fritsch-Carlson); knots
    # are exact on both, interior points deviate by the limiter
    # choice - reported, not asserted equal.
    from scipy.interpolate import PchipInterpolator as _PC

    n = 280
    xr = np.geomspace(1.0e5, 1.0e3, n)
    lr = np.log(xr / 1.0e4)
    xv = 2000.0 * (xr / 1.0e4) ** b_ * (1.0 + a_ * np.sin(k_ * lr))
    pc = _PC(xr[::-1], xv[::-1])
    _, _, _, ac = run_driver(WORKDIR / "vel14_conv280", nd=n)
    v_pc = pc(ac[:, 0])
    s_pc = (ac[:, 0] / ac[:, 1]) * pc.derivative()(ac[:, 0]) - 1.0
    rep["pchip_vs_steffen"] = {
        "max_rel_v_diff": float(
            np.max(np.abs(v_pc - ac[:, 1]) / ac[:, 1])
        ),
        "max_abs_sigma_diff": float(np.max(np.abs(s_pc - ac[:, 2]))),
        "note": "Steffen vs Fritsch-Carlson limiter divergence; both "
        "monotone, neither wrong - see audit doc B6",
    }

    # C7 ltau-shift invariance: kappa rescales
    inv = {}
    for fac in (1.0e-3, 1.0e-1, 10.0, 1.0e3):
        write_snhyd(
            WORKDIR / "vel14_kscale",
            r,
            v,
            sig,
            dens,
            kap * fac,
            age=SN_AGE,
        )
        _, _, _, af = run_driver(WORKDIR / "vel14_kscale")
        inv[str(fac)] = float(
            np.max(np.abs(af[:, 0] - a[:, 0]) / a[:, 0])
        )
    rep["ltau_shift_max_rel_grid_move"] = inv

    # C8 homology rescaling: t -> lam*t with R -> lam*R,
    # rho -> rho/lam^3, V invariant.  Feed the rescaled file at the
    # rescaled age (file Time column tracks the scale so the
    # expansion is identity); require sigma invariant and tau ->
    # tau/lam^2.
    lam = 2.0
    write_snhyd(
        WORKDIR / "vel14_hrescale",
        r * lam,
        v,
        sig,
        dens / lam**3,
        kap,
        age=SN_AGE,  # keep expansion identity: file Time = SN_AGE
    )
    _, _, _, ah = run_driver(WORKDIR / "vel14_hrescale")
    rep["homology_rescale"] = {
        "grid_r_ratio_expected": float(lam),
        "grid_r_ratio": float(ah[0, 0] / a[0, 0]),
        "max_abs_sigma_move": float(np.max(np.abs(ah[:, 2] - a[:, 2]))),
    }
    tau_h = np.zeros(len(ah))
    kap_h = np.exp(np.interp(np.log(ah[:, 0]), np.log(r * lam)[::-1],
                             np.log(kap * dens / lam**3 * 1.0e10)[::-1]))
    for i in range(1, len(ah)):
        tau_h[i] = tau_h[i - 1] + 0.5 * (kap_h[i - 1] + kap_h[i]) * (
            ah[i - 1, 0] - ah[i, 0]
        )
    rep["homology_rescale"]["tau_inner_ratio_vs_expected_0.25"] = float(
        tau_h[-1] / tau_new[-1]
    )

    # D1 continuity: steady-wind Mdot residual on the new grid.
    # File rho is tabulated on file R; interpolate rho on the new
    # grid (log-log) and form R^2 rho V.
    rho_new = np.exp(np.interp(np.log(rg), logr[::-1], np.log(dens)[::-1]))
    mdot = rg**2 * rho_new * vg
    rep["continuity"] = {
        "mdot_rel_spread": float((mdot.max() - mdot.min()) / mdot.mean()),
        "note": "SN homologous field: not a steady-wind profile; "
        "spread measures grid interpolation only",
    }

    OUT_JSON.write_text(json.dumps(rep, indent=2) + "\n")
    print(json.dumps(rep, indent=2))
    print(f"-> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
