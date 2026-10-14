#!/usr/bin/env python3
"""Unit tests for new_main/subs/recompute_sigma_true.f.

Covers the dev7 SIGMA_TRUE machinery: the unclamped centred-difference
sigma = dlnV/dlnR - 1 used by DO_NONMONO_SOB and DO_NONMONO_KNOP, its
boundary one-sided differences, the stagnation-cell branch, the |sigma|
<= 1E6 clamp, the ND<3 degenerate-grid guard, and the DO_NONMONO_SOB
re-clamp of the work array SIGMA at >= -0.999.

Compiles the exact production files:
  tools/set_kind_module.f
  new_main/mod_subs/control_variable_mod.f
  unix/term_io.f
  new_main/subs/recompute_sigma_true.f
with the unit_tests/f90/drv_sigma_true.f90 CLI driver.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import F90DIR, PLOTDIR, TREE, WORKDIR, Report, compile_driver, run

TREE_SRC = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "new_main" / "mod_subs" / "control_variable_mod.f",
    TREE / "unix" / "term_io.f",
    TREE / "new_main" / "subs" / "recompute_sigma_true.f",
]


def get_bin():
    return compile_driver("drv_sigma_true", TREE_SRC + [F90DIR / "drv_sigma_true.f90"])


def write_rv(name, r, v, sig=None):
    p = WORKDIR / name
    sig = np.zeros_like(r) if sig is None else sig
    with open(p, "w") as f:
        f.write(f"{len(r)}\n")
        for a, b, c in zip(r, v, sig):
            f.write(f"{a:.16e} {b:.16e} {c:.16e}\n")
    return p


def parse_sigma(out):
    rows = {}
    n_undef = None
    for ln in out.splitlines():
        t = ln.split()
        if t and t[0] == "N_UNDEF":
            n_undef = int(t[1])
        elif t and t[0].isdigit():
            rows[int(t[0]) - 1] = [float(x) for x in t[1:]]
    n = max(rows) + 1
    arr = np.array([rows[i] for i in range(n)])
    return arr, n_undef


def parse_rows_present(out):
    """True when the driver printed at least one numeric result row."""
    for ln in out.splitlines():
        t = ln.split()
        if t and t[0].isdigit():
            return True
    return False


def sigma_analytic(r, v):
    """Centred-difference sigma = (R/V)dV/dR - 1, one-sided at ends,
    matching the production stencil."""
    s = np.zeros(len(r))
    s[0] = (r[0] / v[0]) * (v[0] - v[1]) / (r[0] - r[1]) - 1.0
    s[-1] = (r[-1] / v[-1]) * (v[-2] - v[-1]) / (r[-2] - r[-1]) - 1.0
    for i in range(1, len(r) - 1):
        s[i] = (r[i] / v[i]) * (v[i - 1] - v[i + 1]) / (r[i - 1] - r[i + 1]) - 1.0
    return s


def main():
    rep = Report("SIGMA_TRUE / RECOMPUTE_SIGMA_TRUE")
    bin_ = get_bin()

    # -- T1: homologous expansion v = a*R -> sigma = 0 everywhere.
    nd = 60
    r = np.linspace(5.0e4, 1.0e3, nd)
    v = 2.0e-4 * r
    f = write_rv("rv_hom.txt", r, v)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    ok = nun == 0 and np.max(np.abs(arr[:, 2])) < 1e-10
    rep.check(
        "T1 homologous v=aR gives sigma_true=0",
        ok,
        f"max|sig|={np.max(np.abs(arr[:, 2])):.2e}, " f"N_UNDEF={nun}",
    )

    # -- T2: coasting shell v=const -> sigma=-1 interior.
    v2 = np.full(nd, 3000.0)
    f = write_rv("rv_const.txt", r, v2)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    interior = arr[1:-1, 2]
    ok = nun == 0 and np.allclose(interior, -1.0, atol=1e-12)
    rep.check(
        "T2 constant v gives sigma_true=-1 interior",
        ok,
        f"interior range [{interior.min():.6f},{interior.max():.6f}]",
    )

    # -- T3: beta-law profile vs analytic centred difference.
    r3 = np.geomspace(8.0e4, 5.0e2, nd)
    v3 = 1.0e4 * (1.0 - 5.0e2 / r3) ** 0.8 + 500.0
    f = write_rv("rv_beta.txt", r3, v3)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    ref = sigma_analytic(r3, v3)
    err = np.abs(arr[:, 2] - ref)
    ok = nun == 0 and np.allclose(arr[:, 2], ref, rtol=0, atol=1e-12)
    rep.check(
        "T3 beta-law matches analytic centred difference",
        ok,
        f"max err {err.max():.2e}",
    )

    # -- T4: velocity reversal -> sigma < -1 at the shock, unclamped
    # (infrastructure SIGMA is clamped elsewhere; SIGMA_TRUE must keep
    # the true gradient for the Dessart+2015 tau floor).
    r4 = np.linspace(5.0e4, 1.0e3, nd)
    v4 = np.where(r4 > 2.0e4, 1.2e4 - 0.15 * (r4 - 2.0e4), 8.0e3)
    v4 = np.clip(v4, 500.0, None)
    f = write_rv("rv_shock.txt", r4, v4)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    isig = arr[:, 2]
    ok = nun == 0 and isig.min() < -1.0 and np.all(isig > -1.0e6)
    rep.check(
        "T4 dv/dr<0 yields sigma_true<-1 (unclamped)", ok, f"min sigma {isig.min():.2f}"
    )

    # -- T5: stagnation cell |V| below VTOL -> sign-adopted, bounded.
    v5 = np.linspace(1.0e4, 3.0e2, nd)
    v5[20] = 1.0e-12  # effectively stagnant
    f = write_rv("rv_stag.txt", r, v5)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    s20 = arr[20, 2]
    ok = np.isfinite(s20) and abs(s20) <= 1.0e6
    rep.check("T5 stagnation cell bounded |sigma|<=1e6", ok, f"sigma(21)={s20:.3e}")

    # -- T6: |sigma| clamp at 1E6 for 1/V blow-up on smooth profiles.
    v6 = np.linspace(1.0e4, 3.0e2, nd)
    v6[30] = 1.0e-6 * v6.max() * 2.0  # just above VTOL -> 1/V amplif.
    f = write_rv("rv_clamp.txt", r, v6)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    ok = np.all(np.abs(arr[:, 2]) <= 1.0e6)
    rep.check(
        "T6 |sigma|>1e6 cells clamped", ok, f"max |sigma| {np.abs(arr[:, 2]).max():.3e}"
    )

    # -- T7: degenerate grids ND<3 report all cells undefined.
    f = write_rv("rv_nd2.txt", r[:2], v[:2])
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    ok = nun == 2 and np.all(arr[:, 2] == 0.0)
    rep.check("T7 ND=2 -> N_UNDEF=2, sigma_true=0", ok, f"N_UNDEF={nun}")

    # -- T8: flag gating.  off: SIGMA untouched AND ST sentinel remains.
    sig_in = np.full(nd, -5.0)
    f = write_rv("rv_gate.txt", r, v, sig_in)
    arr, _ = parse_sigma(run(bin_, "recomp", f, "off").stdout)
    ok = np.all(arr[:, 3] == -999.0) and np.all(arr[:, 2] == -5.0)
    rep.check("T8 flags off -> early return, nothing touched", ok)

    # -- T9: DO_NONMONO_SOB re-clamps SIGMA>=-0.999, SIGMA_TRUE unclamped.
    arr, _ = parse_sigma(run(bin_, "recomp", f, "sob").stdout)
    ok = (
        np.all(arr[:, 2] >= -0.999)
        and np.all(arr[:, 3] != -999.0)
        and np.abs(arr[1, 3]) < 1e-9
    )
    rep.check(
        "T9 SOB: SIGMA clamped >=-0.999, SIGMA_TRUE filled",
        ok,
        f"SIGMA[0]={arr[0,2]:.4f} SIGMA_TRUE[1]={arr[1,3]:.4f}",
    )

    # -- T10: DO_NONMONO_KNOP leaves the work SIGMA unclamped.
    arr, _ = parse_sigma(run(bin_, "recomp", f, "knop").stdout)
    ok = np.all(arr[:, 2] == -5.0) and np.all(arr[:, 3] != -999.0)
    rep.check("T10 KNOP: work SIGMA left unclamped", ok, f"SIGMA[0]={arr[0,2]:.4f}")

    # -- T11: negative velocities handled (infall) without crash.
    v11 = -np.linspace(1.0e4, 3.0e2, nd)
    f = write_rv("rv_infall.txt", r, v11)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    ref = sigma_analytic(r, v11)
    ok = np.all(np.isfinite(arr[:, 2])) and np.allclose(
        arr[:, 2], ref, rtol=0, atol=1e-12
    )
    rep.check("T11 negative V (infall) finite, consistent sign", ok)

    # -- T12: N_UNDEF>0 under a non-mono flag aborts by default.
    # Three equal consecutive radii make the centred denominator of
    # the middle-neighbour cell vanish -> N_UNDEF=1.
    r12 = r.copy()
    r12[9] = r12[10] = r12[11] = r[9]
    f = write_rv("rv_undef.txt", r12, v)
    p = run(bin_, "recomp", f, "sob", check=False)
    rep.check(
        "T12 N_UNDEF>0 aborts under DO_NONMONO_SOB",
        "ALLOW_SIGMA_UNDEF" in p.stdout and not parse_rows_present(p.stdout),
    )

    # -- T13: ALLOW_SIGMA_UNDEF=T downgrades the abort to a warning.
    p = run(bin_, "recomp", f, "sob", "allow", check=False)
    arr, _ = parse_sigma(p.stdout)
    rep.check(
        "T13 ALLOW_SIGMA_UNDEF warns and completes",
        "no usable velocity gradient" in p.stdout and arr[10, 3] == 0.0,
    )

    # -- Plot: sigma_true across a non-monotonic profile with a shock
    # plus a stagnation cell -- the regimes the flag was built for.
    rp = np.geomspace(6.0e4, 4.0e2, 160)
    vp = 1.1e4 - 0.22 * (rp - 1.0e4)
    vp = np.where(rp < 1.0e4, 8.8e3 + 0.05 * (1.0e4 - rp), vp)
    vp = np.clip(vp, 600.0, None)
    vp[80] = 1.0e-9
    f = write_rv("rv_plot.txt", rp, vp)
    arr, _ = parse_sigma(run(bin_, "core", f).stdout)
    fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    ax[0].plot(rp, vp, "b-")
    ax[0].set_ylabel("v [km/s]")
    ax[0].set_title("non-monotonic velocity profile")
    ax[1].plot(rp, arr[:, 2], "r-")
    ax[1].axhline(-1.0, color="k", ls="--", lw=0.8, label="infrastructure clamp -1")
    ax[1].axhline(0.0, color="grey", ls=":", lw=0.8)
    ax[1].set_xlabel("R [1e10 cm]")
    ax[1].set_ylabel(r"$\sigma_{true}=d\ln V/d\ln R-1$")
    ax[1].legend()
    fig.tight_layout()
    out = PLOTDIR / "sigma_true_nonmono.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
