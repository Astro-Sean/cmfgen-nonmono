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


def sigma_dvdr(r, v):
    """Independent mirror of the production hybrid stencil (audit
    30-Oct-2026): sigma = d ln|V|/d lnR - 1 via the non-uniform
    quadratic in lnR where the stencil carries one sign and |V| spans
    less than a decade (exact for
    V ~ R^b on any grid); the linear-space quadratic dV/dR at sign
    changes; zero V(i) -> undefined.  Kept in python so the test
    double-checks the formula rather than compiling it twice."""
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
            dvdr = (
                v[i - 1] * hm / (hp * (hp + hm))
                + v[i] * (hp - hm) / (hp * hm)
                - v[i + 1] * hp / (hm * (hp + hm))
            )
            s[i] = (r[i] / v[i]) * dvdr - 1.0
        else:
            s[i] = 0.0
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
        dvdr = (
            v[0] * (2.0 * hp + hm) / (hp * (hp + hm))
            - v[1] * (hp + hm) / (hp * hm)
            + v[2] * hp / (hm * (hp + hm))
        )
        s[0] = (r[0] / v[0]) * dvdr - 1.0
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
        dvdr = (
            v[-3] * hp / (hm * (hp + hm))
            - v[-2] * (hp + hm) / (hp * hm)
            + v[-1] * (2.0 * hp + hm) / (hp * (hp + hm))
        )
        s[-1] = (r[-1] / v[-1]) * (-dvdr) - 1.0
    return np.clip(s, -1.0e6, 1.0e6)


def sigma_true_exact(r, v, dvdr):
    """sigma = (R/V)dV/dR - 1 with the TRUE analytic dV/dR."""
    return (r / v) * dvdr - 1.0


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

    # -- T3: beta-law profile on a log-spaced (non-uniform) grid.
    # The 30-Oct-2026 audit replaced the secant stencil with the
    # non-uniform quadratic; on this grid the secant carries an O(dx)
    # bias vs the true analytic sigma, so the reference is the true
    # derivative, not the old stencil.  Cells where |sigma| saturates
    # the 1e6 clip (near the r=500 singular edge) are compared against
    # the clip, not the unbounded truth.
    r3 = np.geomspace(8.0e4, 5.0e2, nd)
    v3 = 1.0e4 * (1.0 - 5.0e2 / r3) ** 0.8 + 500.0
    f = write_rv("rv_beta.txt", r3, v3)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    dvdr3 = 1.0e4 * 0.8 * (1.0 - 5.0e2 / r3) ** (-0.2) * 5.0e2 / r3**2
    ref3 = sigma_true_exact(r3, v3, dvdr3)
    ref3 = np.clip(ref3, -1.0e6, 1.0e6)
    # Resolved-region check: near the r=500 singular edge the true
    # dV/dR diverges and any 3-point stencil has O(10%) truncation;
    # away from it the log-space quadratic must track truth closely.
    resolved = (r3 > 2.0e3) & (np.abs(ref3) < 9.0e5)
    err3 = np.abs(arr[:, 2] - ref3)[resolved]
    rel3 = err3 / np.maximum(np.abs(ref3)[resolved], 1.0)
    ok = nun == 0 and rel3.max() < 5.0e-3
    rep.check(
        "T3 beta-law on log grid vs true sigma (resolved region)",
        ok,
        f"max rel err {rel3.max():.2e} on {resolved.sum()} cells",
    )
    # the hybrid stencil must also match its independent python mirror
    ref3s = sigma_dvdr(r3, v3)
    err3s = np.abs(arr[:, 2] - ref3s)
    ok = np.allclose(arr[:, 2], ref3s, rtol=0, atol=1e-8)
    rep.check(
        "T3b production stencil == independent hybrid mirror",
        ok,
        f"max err {err3s.max():.2e}",
    )

    # -- T3c: exponential profile v = a*exp(b*R) on a non-uniform
    # grid: the quadratic stencil must track the true sigma, not just
    # the old secant's answer.  (R/V)dV/dR - 1 = b*R - 1 exactly.
    r3c = np.geomspace(8.0e4, 5.0e2, nd)
    b3c = 3.0e-5
    v3c = 2.0e3 * np.exp(b3c * r3c)
    f = write_rv("rv_exp.txt", r3c, v3c)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    ref3c = np.clip(b3c * r3c - 1.0, -1.0e6, 1.0e6)
    err3c = np.abs(arr[:, 2] - ref3c) / np.maximum(np.abs(ref3c), 1.0)
    ok = nun == 0 and err3c.max() < 5.0e-2
    rep.check(
        "T3c exponential v on log grid vs true sigma=bR-1",
        ok,
        f"max rel err {err3c.max():.2e}",
    )

    # -- T3d: smooth zero crossing v = a*(R-R0): sigma=(R/V)dV/dR-1
    # diverges at the crossing with sign sign(V)*sign(dV/dR); the
    # signed-stagnation branch must keep that sign on both sides.
    r3d = np.linspace(5.0e4, 1.0e3, nd)
    r0 = 2.0e4
    v3d = 1.5e-2 * (r3d - r0)
    f = write_rv("rv_xing.txt", r3d, v3d)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    dvdr3d = np.full(nd, 1.5e-2)
    ref3d = np.clip(
        np.where(v3d != 0.0, (r3d / v3d) * dvdr3d - 1.0, 0.0),
        -1.0e6,
        1.0e6,
    )
    i_x = int(np.argmin(np.abs(v3d)))
    # Stencil-equality vs the independent hybrid mirror, plus the
    # classification sign on both sides of the crossing: sigma keeps
    # the algebraic-limit sign sign(V)*sign(dV/dR).  Within ~3 cells
    # of the crossing ln|V| is kinked and the discretization carries
    # O(10%) magnitude error on the log-space cells (documented in
    # the routine header); the sign and bound are what matter there.
    ref3m = sigma_dvdr(r3d, v3d)
    mirr_err = np.abs(arr[:, 2] - ref3m)
    # Cells 3-5 from the crossing use the log-space branch (same
    # sign, <1 decade) where ln|V| still carries the kink; that is
    # the documented trade-off for power-law exactness elsewhere.
    # Assert: mirror-exact everywhere, bounded relative error near
    # the crossing, near-exact far from it, correct sign at +-1.
    far = np.abs(np.arange(nd) - i_x) >= 6
    far_err = np.abs(arr[:, 2] - ref3d)[far]
    near = np.abs(np.arange(nd) - i_x) < 6
    rel_near = (
        np.abs(arr[:, 2] - ref3d)[near] / np.maximum(np.abs(ref3d)[near], 1.0)
    )
    ok = (
        mirr_err.max() < 1.0e-8
        and far_err.max() < 5.0e-2
        and rel_near.max() < 0.15
        and np.sign(arr[i_x - 1, 2]) == np.sign(ref3d[i_x - 1])
        and np.sign(arr[i_x + 1, 2]) == np.sign(ref3d[i_x + 1])
        and np.all(np.abs(arr[:, 2]) <= 1.0001e6)
    )
    rep.check(
        "T3d zero crossing: mirror-exact, signed, bounded",
        ok,
        f"mirror {mirr_err.max():.2e} far {far_err.max():.2e} "
        f"near-rel {rel_near.max():.2e} "
        f"sigma(xing+/-1)={arr[i_x - 1, 2]:.2e}/{arr[i_x + 1, 2]:.2e}",
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
    # Same-sign negative stencils take the log-space branch (|V|
    # spans <1 decade here): sigma = d ln|V|/d lnR - 1, which equals
    # (R/V)dV/dR - 1 algebraically.  Compare vs the hybrid mirror.
    v11 = -np.linspace(1.0e4, 3.0e2, nd)
    f = write_rv("rv_infall.txt", r, v11)
    arr, nun = parse_sigma(run(bin_, "core", f).stdout)
    ref = sigma_dvdr(r, v11)
    dvdr11 = np.gradient(v11, r)
    ref11 = np.clip(sigma_true_exact(r, v11, dvdr11), -1.0e6, 1.0e6)
    err11 = np.abs(arr[:, 2] - ref11) / np.maximum(np.abs(ref11), 1.0)
    ok = (
        np.all(np.isfinite(arr[:, 2]))
        and np.allclose(arr[:, 2], ref, rtol=0, atol=1e-12)
        and err11.max() < 5.0e-2
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
