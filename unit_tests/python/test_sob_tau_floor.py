#!/usr/bin/env python3
"""Unit tests for the Dessart, Audit & Hillier (2015) tau floor.

Files under test (identical edit in each):
  subs/sobjbar_sim.f, subs/sobew.f,
  subs/sobew_grad_v2.f, subs/sobew_grad_v3.f

For sigma < 0 the factor |1+sigma*mu^2| can vanish along a ray; the
edit replaces it by max(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4).
For sigma >= -1 the expression is bit-identical to upstream.

Method: one driver is linked twice per routine, once against the
jfnk-tree source and once against ../cmfgen_reference/cur_cmf, so the
same inputs give a direct new-vs-reference comparison.  For SOBJBAR_SIM
the angle-averaged escape probability BETA = sum_L w_L*EXF(GAM_L) is
recomputed analytically in Python and compared exactly.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import BINDIR, F90DIR, PLOTDIR, REF, TREE, WORKDIR, \
    Report, compile_objs, link, run

HELPERS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "mod_trap_derivatives.f",
    TREE / "unix" / "term_io.f",
    TREE / "subs" / "zalongp.f",
    TREE / "subs" / "nordtau.f",
    TREE / "subs" / "tcompd.f",
    TREE / "subs" / "xvecd.f",
    TREE / "subs" / "thomas.f",
    TREE / "subs" / "thomas_rh.f",
    TREE / "subs" / "exponx.f",
]

ROUTINES = {"sobjsim": "sobjbar_sim.f", "sobew": "sobew.f",
            "sobewg2": "sobew_grad_v2.f", "sobewg3": "sobew_grad_v3.f"}

ND, NC, NP, NSIM = 25, 4, 4, 1


def exponx(x):
    """Python copy of subs/exponx.f."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    s = np.abs(x) < 1.0e-3
    out[s] = 1.0 - x[s] * (0.5 - x[s] / 6.0 * (1.0 - x[s] / 4.0))
    m = (~s) & (x < 40.0)
    out[m] = (1.0 - np.exp(-x[m])) / x[m]
    out[~s & ~m] = 1.0 / x[~s & ~m]
    return out


def build():
    """Compile helpers + driver once, each sob file twice (new/ref)."""
    help_objs = compile_objs("sob_help", HELPERS)
    drv_obj = compile_objs("sob_drv", [F90DIR / "drv_sob_floor.f90"],
                           moddir=BINDIR / "obj_sob_help")
    bins = {}
    for tag, root in (("new", TREE), ("ref", REF)):
        sob = []
        for key, fname in ROUTINES.items():
            sob += compile_objs(f"sob_{tag}_{key}",
                                [root / "subs" / fname],
                                moddir=BINDIR / "obj_sob_help")
        for key in ROUTINES:
            bins[(key, tag)] = link(
                f"sobfloor_{tag}",
                help_objs + drv_obj + sob,
                moddir=BINDIR / "obj_sob_help")
    return bins


def write_case(name, sigma, p_frac=None, chil0=1.0e6):
    """Write a case file; returns (path, gamh, p, w)."""
    r = np.geomspace(1.0e5, 1.0e3, ND)      # outer -> inner
    v = np.linspace(8.0e3, 1.0e3, ND)
    sig = np.full(ND, sigma)
    if p_frac is None:
        p_frac = [0.50, 0.75, 0.90, 0.99]
    p = np.array(p_frac) * r[-1]
    w = np.full(NP, 1.0 / NP)
    chi = np.full(ND, 1.0e-8)
    esec = np.full(ND, 1.0e-9)
    chil = np.full(ND, chil0)
    etal = chil * 2.0                     # S = ETAL/CHIL = 2
    src = np.full(ND, 2.0)
    path = WORKDIR / name
    with open(path, "w") as f:
        f.write(f"{ND} {NC} {NP} {NSIM}\n")
        f.write("1.0 F T F 0.0 0.0 1.0e40\n")
        f.write("LINEAR ZERO_FLUX\n")
        for i in range(ND):
            f.write(f"{r[i]:.10e} {v[i]:.10e} {sig[i]:.10e} "
                    f"{chi[i]:.4e} {esec[i]:.4e} {chil[i]:.6e} "
                    f"{etal[i]:.6e} {src[i]:.6e}\n")
        for j in range(NP):
            f.write(f"{p[j]:.10e} {w[j]:.6e} {w[j]:.6e}\n")
        for _ in range(NSIM):
            f.write(f"{chil0:.6e} {2.0*chil0:.6e} 1.0\n")
    gamh = chil0 * 3.0e-10 * r / v / 1.0   # FL=1
    return path, r, p, w, gamh


def beta_predict(sig, r, p, w, gamh, floor=True):
    mu2 = 1.0 - (p[None, :] / r[:, None]) ** 2
    t1 = np.abs(1.0 + sig * mu2)
    if floor and sig < 0.0:
        t1 = np.maximum(t1, max(-0.2 * (1.0 + sig), 1.0e-4))
    gam = gamh[:, None] / t1
    return (exponx(gam) * w[None, :]).sum(axis=1)


def max_meaningful_diff(a, b):
    """Largest |a-b| ignoring denormal garbage from uninit memory."""
    m = (np.abs(a) > 1e-290) | (np.abs(b) > 1e-290)
    return np.abs(a[m] - b[m]).max() if m.any() else 0.0


def parse_rows(out):
    head, rows = {}, {}
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "EW_CF":
            head["EW"] = float(t[1])
        elif t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    n = max(rows) + 1
    return np.array([rows[i] for i in range(n)]), head


def main():
    rep = Report("Sobolev tau floor (Dessart+2015)")
    bins = build()

    # -- T1: sigma >= 0 gives bit-identical BETA for all four routines
    # (the edited branch is not reached).
    for key in ROUTINES:
        case, r, p, w, gamh = write_case(f"case_{key}_pos.txt", 2.0)
        a, _ = parse_rows(run(bins[(key, "new")], key, case).stdout)
        b, _ = parse_rows(run(bins[(key, "ref")], key, case).stdout)
        ok = np.array_equal(a, b) or max_meaningful_diff(a, b) == 0.0
        rep.check(f"T1 {key}: sigma=+2 identical to upstream", ok,
                  f"maxdiff={np.abs(a - b).max():.2e}")

    # -- T2: sigma=-0.5 (> -1): floor can never engage, still identical.
    for key in ROUTINES:
        case, r, p, w, gamh = write_case(f"case_{key}_mhalf.txt", -0.5)
        a, _ = parse_rows(run(bins[(key, "new")], key, case).stdout)
        b, _ = parse_rows(run(bins[(key, "ref")], key, case).stdout)
        ok = np.array_equal(a, b) or max_meaningful_diff(a, b) == 0.0
        rep.check(f"T2 {key}: sigma=-0.5 identical to upstream", ok,
                  f"maxdiff={np.abs(a - b).max():.2e}")

    # -- T3: sigma=-2 (< -1): new BETA matches the floored formula
    # exactly, and differs from the unclamped reference.
    # CHIL=1e9 puts GAMH ~ 0.3-4, large enough that EXF(GAM) is
    # sensitive to the floor.
    case, r, p, w, gamh = write_case("case_sim_m2.txt", -2.0,
                                     chil0=1.0e9)
    a, _ = parse_rows(run(bins[("sobjsim", "new")], "sobjsim",
                          case).stdout)
    b, _ = parse_rows(run(bins[("sobjsim", "ref")], "sobjsim",
                          case).stdout)
    beta_pred = beta_predict(-2.0, r, p, w, gamh, floor=True)
    beta_ref_pred = beta_predict(-2.0, r, p, w, gamh, floor=False)
    err_new = np.abs(a[:, 0] - beta_pred).max()
    err_ref = np.abs(b[:, 0] - beta_ref_pred).max()
    ok = err_new < 1e-9 and err_ref < 1e-9 and \
        np.abs(a[:, 0] - b[:, 0]).max() > 1e-3
    rep.check("T3 sobjsim sigma=-2: BETA == floored formula, "
              "ref == unclamped formula", ok,
              f"err_new={err_new:.2e} err_ref={err_ref:.2e} "
              f"diff={np.abs(a[:,0]-b[:,0]).max():.2e}")

    # -- T4: at sigma=-2 the new BETA is bounded and finite at every
    # depth; the reference drives (1-e^-GAM)/GAM -> 0 along rays whose
    # mu^2 hits the pole (BETA undercounted).
    ok = np.all(np.isfinite(a[:, 0])) and np.all(a[:, 0] > 0.0) and \
        np.all(a[:, 0] >= b[:, 0] - 1e-15)
    rep.check("T4 sobjsim sigma=-2: floored BETA bounded, "
              ">= reference", ok,
              f"new min {a[:,0].min():.3f} ref min {b[:,0].min():.3f}")

    # -- T5: same behaviour in the EW/force routines: sigma=-2 gives
    # different-but-finite output vs reference.
    for key in ("sobew", "sobewg2", "sobewg3"):
        case, r, p, w, gamh = write_case(f"case_{key}_m2.txt", -2.0,
                                         chil0=1.0e9)
        a, ha = parse_rows(run(bins[(key, "new")], key, case).stdout)
        b, hb = parse_rows(run(bins[(key, "ref")], key, case).stdout)
        ok = np.all(np.isfinite(a)) and \
            max_meaningful_diff(a, b) > 0.0 and ha["EW"] != hb["EW"]
        rep.check(f"T5 {key}: sigma=-2 outputs differ, stay finite",
                  ok, f"dEW={ha['EW']-hb['EW']:.3e} "
                  f"maxdiff={np.abs(a-b).max():.2e}")

    # -- T6: sigma=-1 edge case: the 1e-4 term bounds the mu->1 rays.
    case, r, p, w, gamh = write_case("case_sim_m1.txt", -1.0)
    a, _ = parse_rows(run(bins[("sobjsim", "new")], "sobjsim",
                          case).stdout)
    beta_pred = beta_predict(-1.0, r, p, w, gamh, floor=True)
    err = np.abs(a[:, 0] - beta_pred).max()
    rep.check("T6 sobjsim sigma=-1: 1e-4 bound matches formula",
              err < 1e-9, f"err={err:.2e}")

    # -- Plot: BETA(I) new vs ref vs analytic for sigma=-2.
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    pred = beta_predict(-2.0, r, p, w, gamh, floor=True)
    ax[0].plot(r, a[:, 0], "ro-", ms=3, label="new (floored)")
    ax[0].plot(r, b[:, 0], "bx-", ms=3, label="reference")
    ax[0].plot(r, pred, "k--", lw=0.8, label="analytic floor")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("R [1e10 cm]")
    ax[0].set_ylabel(r"BETA(I)")
    ax[0].set_title(r"$\sigma=-2$: escape probability per depth")
    ax[0].legend()
    mu2 = 1.0 - (p[None, :] / r[:, None]) ** 2
    t1n = np.maximum(np.abs(1.0 - 2.0 * mu2),
                     max(-0.2 * (1.0 - 2.0), 1e-4))
    ax[1].pcolormesh(r, p, np.log10((gamh[:, None] / t1n).T),
                     shading="nearest", cmap="viridis")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("R [1e10 cm]")
    ax[1].set_ylabel("P [1e10 cm]")
    ax[1].set_title(r"new $\log_{10}$ GAM per (depth, ray)")
    fig.tight_layout()
    out = PLOTDIR / "sob_tau_floor.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
