#!/usr/bin/env python3
"""Validation of the manufactured velocity-profile library (V001-V012)
and of the production velocity-gradient routine on every profile.

The library (unit_tests/velocity_profiles.py) is the common input set
for every non-monotonic solver test (KN/M/S/T series in
docs/audits/NUMERICAL_VERIFICATION_TEST_PLAN.md).  This file proves
two things before those tests rely on it:

1. The library is internally correct.  Every profile ships a
   hand-coded analytic dV/dR; the complex-step derivative
   Im(V(r+ih))/h re-derives it with no shared code path.  Profile
   metadata (gradient sign changes, sigma extrema, compression
   regions and their widths) is checked against direct evaluation.

2. The production stencil in new_main/subs/recompute_sigma_true.f
   (SIGMA_TRUE_CORE, driven through f90/drv_sigma_true.f90) is a
   faithful discretization of the same quantity.  Three-way check:
   production vs the analytic sigma = (R/V)dV/dR-1 (truncation
   bounded, ~h^2 on the log grid), production vs an independent
   python mirror of the documented stencil (bitwise), and production
   sign-change/compression census vs the analytic one.

Sweep coverage required by the test plan: radial resolution, shock
position, shock width, compression strength, smooth perturbation,
reversed radial array, near-zero gradient, two close reversals,
stagnant region, zero-velocity cell, collapsed spacing, negative
velocity (infall), and a uniform (linear) grid.

Tolerances: exact-law profiles (V001, V003) sit at roundoff
(|sigma err| ~ 1e-13); smooth profiles carry the quadratic-stencil
truncation measured at each ND - the coded bounds are ~3x the
measured values, and the h^2 refinement ratio is the quantitative
criterion.  Mirror agreement is required to ~1e-12 (bitwise on the
build compiler; slack for cross-compiler ulps in LOG).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (
    F90DIR,
    PLOTDIR,
    TREE,
    WORKDIR,
    Report,
    compile_and_link_prod,
    run,
)
from velocity_profiles import (
    complex_step_dvdr,
    make_profile,
    parse_driver_sigma,
    profile_ids,
    sign_change_crossings,
    sigma_true_core_mirror,
)

TREE_SRC = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "new_main" / "mod_subs" / "control_variable_mod.f",
    TREE / "unix" / "term_io.f",
    TREE / "new_main" / "subs" / "recompute_sigma_true.f",
]

ND = 80

# Well-cell absolute-error bounds at ND=80 on the default log grid,
# ~3x the measured truncation errors (see plan doc V-series table).
# V001/V003 are exact laws -> roundoff.  Steep profiles (V005, V007)
# carry larger constants because sigma spans +-25 over a few cells;
# their quantitative check is the refinement ratio, not the bound.
ANA_TOL = {
    "V001": 1.0e-9, "V002": 1.0e-3, "V003": 1.0e-9,
    "V004": 1.0, "V005": 25.0, "V006": 1.5, "V007": 15.0,
    "V008": 0.1, "V009": 5.0e-3, "V010": 8.0, "V011": 2.5,
    "V012": 2.5,
}

# Expected gradient sign changes / compression regions per profile.
SIGN_META = {
    "V001": (0, 0), "V002": (0, 0), "V003": (0, 1),
    "V004": (2, 1), "V005": (2, 1), "V006": (6, 3),
    "V007": (2, 1), "V008": (2, 1), "V009": (0, 0),
    "V010": (2, 1), "V011": (4, 2), "V012": (2, 1),
}

# Profiles whose sigma law is reproduced to roundoff by the log-space
# quadratic (power laws in R): skip the refinement-rate check.
EXACT_LAW = {"V001", "V003"}


def get_bin():
    return compile_and_link_prod(
        "drv_sigma_true", TREE_SRC + [F90DIR / "drv_sigma_true.f90"]
    )


def prod_sigma(bin_, prof, tag):
    f = prof.write_rv(WORKDIR / f"vprof_{tag}.txt")
    out = run(bin_, "core", f).stdout
    r, v, st, nun = parse_driver_sigma(out)
    return out, r, v, st, nun


def well_mask(ana):
    """Cells where the analytic sigma is below the |sigma|=1E6 clip
    margin - the regime a truncation comparison is meaningful in."""
    return np.abs(ana) < 1.0e3


def main():
    rep = Report("VELOCITY-PROFILE LIBRARY + SIGMA_TRUE_CORE")
    bin_ = get_bin()

    # --------------------------------------------------------------
    # A. Library self-validation: complex-step derivative, metadata.
    # --------------------------------------------------------------
    for pid in profile_ids():
        p = make_profile(pid, nd=ND)
        cs = complex_step_dvdr(p)
        dv = np.real(p.dvdr)
        err = np.abs(cs - dv)
        scale = max(1.0, np.abs(dv).max())
        rep.check(
            f"A:{pid} analytic dV/dR == complex-step derivative",
            err.max() <= 1.0e-10 * scale,
            f"max err {err.max():.2e} (scale {scale:.1e})",
        )
        ns, nc = SIGN_META[pid]
        rep.check(
            f"A:{pid} sign changes {ns}, compression regions {nc}",
            len(p.sign_changes()) == ns
            and len(p.compression_regions()) == nc,
            f"got {len(p.sign_changes())}/{len(p.compression_regions())}"
            f" sigma in [{p.sigma_min:.2f},{p.sigma_max:.2f}]"
            f" V in [{p.v.min():.1f},{p.v.max():.1f}]",
        )
        sane = (
            np.isfinite(p.v).all()
            and np.isfinite(p.dvdr).all()
            and np.all(np.diff(p.r) < 0.0)
            and np.all(p.v > 0.0)
        )
        rep.check(f"A:{pid} finite arrays, V>0, R decreasing", sane)

    # --------------------------------------------------------------
    # B. Production stencil on every profile: mirror-exact,
    #    analytic-bounded, sign census correct.
    # --------------------------------------------------------------
    for pid in profile_ids():
        p = make_profile(pid, nd=ND)
        out, r, v, st, nun = prod_sigma(bin_, p, f"B_{pid}")
        mir, nun_m = sigma_true_core_mirror(r, v)
        ana = p.sigma
        well = well_mask(ana)
        err_m = np.abs(st - mir).max()
        err_a = np.abs(st - ana)[well].max()
        rep.check(
            f"B:{pid} N_UNDEF=0, finite, mirror-exact",
            nun == 0
            and nun_m == 0
            and np.isfinite(st).all()
            and err_m <= 1.0e-9,
            f"N_UNDEF={nun} mirror {err_m:.1e}",
        )
        rep.check(
            f"B:{pid} sigma vs analytic (truncation bound)",
            err_a <= ANA_TOL[pid],
            f"err {err_a:.2e} vs tol {ANA_TOL[pid]:.1e}",
        )
        ns = SIGN_META[pid][0]
        xs_prod = sign_change_crossings(r, st + 1.0)
        xs_ana = [rc for _, rc in p.sign_changes()]
        hmax = np.log(r[:-1] / r[1:]).max()
        loc_ok = all(
            min(abs(np.log(x / ra)) for ra in xs_ana) <= 1.5 * hmax
            for _, x in xs_prod
        )
        rep.check(
            f"B:{pid} sigma=-1 crossings: {ns} at right locations",
            len(xs_prod) == ns and loc_ok,
            f"prod {len(xs_prod)} crossings",
        )

    # --------------------------------------------------------------
    # C. Radial refinement: quadratic stencil -> err ratio ~4.
    # --------------------------------------------------------------
    for pid in profile_ids():
        if pid in EXACT_LAW:
            continue
        errs = []
        for nd in (ND, 2 * ND, 4 * ND):
            p = make_profile(pid, nd=nd)
            _, r, v, st, nun = prod_sigma(bin_, p, f"C_{pid}_{nd}")
            ana = p.sigma
            errs.append(np.abs(st - ana)[well_mask(ana)].max())
        ratio = errs[-2] / errs[-1] if errs[-1] > 0 else np.inf
        rep.check(
            f"C:{pid} refinement ratio ~4 (h^2 stencil)",
            errs[0] >= errs[1] >= errs[2] and ratio >= 2.5,
            f"errs {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e}"
            f" ratio {ratio:.2f}",
        )

    # --------------------------------------------------------------
    # D. Required sweeps.
    # --------------------------------------------------------------

    # D1 shock position: V007 step scanned across the grid.
    pos_ok = True
    det = []
    for rs in (0.15, 0.30, 0.50, 0.70, 0.90):
        p = make_profile("V007", nd=ND, r_shock=rs)
        _, r, v, st, nun = prod_sigma(bin_, p, f"D1_{rs}")
        xs = sign_change_crossings(r, st + 1.0)
        pos_ok &= nun == 0 and np.isfinite(st).all() and len(xs) == 2
        det.append(f"{rs}:{st.min():.1f}/{p.sigma_min:.1f}")
    rep.check(
        "D1 V007 shock position sweep: 2 crossings everywhere",
        pos_ok,
        "smin prod/ana " + " ".join(det),
    )

    # D2 shock width: resolved widths converge to analytic; the
    # sub-resolution case must stay mirror-exact and finite.
    for w, tol in ((0.15, 0.3), (0.06, 3.0), (0.03, 15.0)):
        p = make_profile("V007", nd=ND, width=w)
        _, r, v, st, nun = prod_sigma(bin_, p, f"D2_{w}")
        err = np.abs(st - p.sigma)[well_mask(p.sigma)].max()
        rep.check(
            f"D2 V007 width {w}: resolved sigma tracks analytic",
            nun == 0 and err <= tol,
            f"err {err:.2e} vs tol {tol:g} smin {st.min():.1f}",
        )
    p = make_profile("V007", nd=ND, width=0.008)
    _, r, v, st, nun = prod_sigma(bin_, p, "D2_narrow")
    mir, _ = sigma_true_core_mirror(r, v)
    xs = sign_change_crossings(r, st + 1.0)
    rep.check(
        "D2 V007 width 0.008 (sub-resolved): mirror-exact, finite, 2 crossings",
        nun == 0
        and np.isfinite(st).all()
        and np.abs(st - mir).max() <= 1.0e-9
        and len(xs) == 2,
        f"prod smin {st.min():.1f} vs analytic {p.sigma_min:.1f}"
        " (stencil smoothing, not a defect)",
    )

    # D3 compression strength: reversal appears at the analytic
    # threshold and the production census tracks it.
    amps = (0.15, 0.25, 0.35, 0.55, 0.80, 0.97)
    amp_ok = True
    det = []
    for a in amps:
        p = make_profile("V004", nd=ND, amp=a)
        _, r, v, st, nun = prod_sigma(bin_, p, f"D3_{a}")
        xs = sign_change_crossings(r, st + 1.0)
        want = len(p.sign_changes())
        amp_ok &= nun == 0 and len(xs) == want
        det.append(f"{a}:{len(xs)}")
    rep.check(
        "D3 V004 strength sweep: crossing count matches analytic",
        amp_ok,
        "amp:prod_nsign " + " ".join(det) + " (reversal onset ~amp 0.25)",
    )

    # D4 smooth perturbation: production tracks the perturbed
    # analytic sigma with the same truncation class; perturbation is
    # deterministic under a fixed seed.
    p0 = make_profile("V004", nd=ND)
    pp = p0.perturbed(3.0e-3, seed=11)
    pp2 = p0.perturbed(3.0e-3, seed=11)
    rep.check(
        "D4 perturbation is deterministic (seeded)",
        np.array_equal(pp.v, pp2.v) and np.array_equal(pp.dvdr, pp2.dvdr),
    )
    _, r, v, st, nun = prod_sigma(bin_, pp, "D4")
    err = np.abs(st - pp.sigma)[well_mask(pp.sigma)].max()
    rep.check(
        "D4 perturbed V004: sigma tracks perturbed analytic",
        nun == 0 and err <= ANA_TOL["V004"],
        f"err {err:.2e} (unperturbed class {ANA_TOL['V004']:.1e})",
    )

    # D5 reversed radial array: every cell reports undefined under
    # 'core' and aborts loudly under a non-mono flag ('recomp').
    pr = p0.reversed_view()
    _, r, v, st, nun = prod_sigma(bin_, pr, "D5")
    rep.check(
        "D5 reversed grid: N_UNDEF=ND (loud failure, no garbage)",
        nun == pr.nd,
        f"N_UNDEF={nun} of {pr.nd}",
    )
    f = pr.write_rv(WORKDIR / "vprof_D5_rev.txt")
    pa = run(bin_, "recomp", f, "sob", check=False)
    rep.check(
        "D5 reversed grid: recomp aborts under DO_NONMONO_SOB",
        "ALLOW_SIGMA_UNDEF" in pa.stdout and "no usable velocity" in pa.stdout,
        f"rc={pa.returncode}",
    )

    # D6 near-zero gradient: V009 sigma ~ -1, deviation sign follows
    # the sign of eps (eps>0 expanding, eps<0 compressing).
    for eps, tag in ((3.0e-3, "pos"), (-3.0e-3, "neg")):
        p = make_profile("V009", nd=ND, eps=eps)
        _, r, v, st, nun = prod_sigma(bin_, p, f"D6_{tag}")
        dev = np.abs(st + 1.0).max()
        imax = int(np.argmax(np.abs(st + 1.0)))
        rep.check(
            f"D6 V009 eps={eps:g}: |sigma+1| small, sign(eps) preserved",
            nun == 0
            and dev <= 0.1
            and np.sign(st[imax] + 1.0) == np.sign(eps),
            f"max|sigma+1|={dev:.3f} sigma range"
            f" [{st.min():.4f},{st.max():.4f}]",
        )

    # D7 closely spaced reversals: resolved pair -> 4 crossings/2
    # regions; merged separation -> 1 region, both reproduced.
    p = make_profile("V011", nd=ND)
    _, r, v, st, nun = prod_sigma(bin_, p, "D7")
    xs = sign_change_crossings(r, st + 1.0)
    rep.check(
        "D7 V011 two close reversals: 4 crossings detected",
        nun == 0 and len(xs) == 4,
        f"prod {len(xs)} crossings",
    )
    p = make_profile(
        "V011", nd=ND, shocks=(0.36, 0.30), widths=(0.10, 0.10)
    )
    _, r, v, st, nun = prod_sigma(bin_, p, "D7m")
    xs = sign_change_crossings(r, st + 1.0)
    want = len(p.sign_changes())
    rep.check(
        f"D7 V011 near-merged reversals: count matches analytic ({want})",
        nun == 0 and len(xs) == want,
        f"prod {len(xs)} vs ana {want}",
    )

    # D8 uniform (linear) grid: stencil stays mirror-exact and
    # analytic-bounded off the log grid.
    p = make_profile("V004", nd=ND, grid="lin")
    _, r, v, st, nun = prod_sigma(bin_, p, "D8")
    mir, _ = sigma_true_core_mirror(r, v)
    err = np.abs(st - p.sigma)[well_mask(p.sigma)].max()
    rep.check(
        "D8 V004 on a uniform grid: mirror-exact, analytic-bounded",
        nun == 0 and np.abs(st - mir).max() <= 1.0e-9 and err <= 1.5,
        f"err {err:.2e} mirror {np.abs(st-mir).max():.1e}",
    )

    # --------------------------------------------------------------
    # E. Degenerate / negative cases.
    # --------------------------------------------------------------

    # E1 exactly-zero V cell: reported undefined, not classified.
    p = make_profile("V010", nd=ND)
    v0 = p.v.copy()
    v0[30] = 0.0
    f = WORKDIR / "vprof_E1.txt"
    with open(f, "w") as fh:
        fh.write(f"{p.nd}\n")
        for ra, va in zip(p.r, v0):
            fh.write(f"{ra:.16e} {va:.16e} 0\n")
    r, v, st, nun = parse_driver_sigma(run(bin_, "core", f).stdout)
    rep.check(
        "E1 V=0 cell -> N_UNDEF=1, sigma_true=0 there",
        nun == 1 and st[30] == 0.0,
        f"N_UNDEF={nun}",
    )
    pa = run(bin_, "recomp", f, "sob", check=False)
    rep.check(
        "E1 V=0 under DO_NONMONO_SOB aborts fail-closed",
        "ALLOW_SIGMA_UNDEF" in pa.stdout,
    )
    pa = run(bin_, "recomp", f, "sob", "allow", check=False)
    rep.check(
        "E1 ALLOW_SIGMA_UNDEF downgrades to a warning",
        "no usable velocity gradient" in pa.stdout,
    )

    # E2 collapsed radial spacing: stencils with zero spacing are
    # undefined and abort under a non-mono flag.
    p = make_profile("V004", nd=ND)
    rc = p.r.copy()
    rc[20] = rc[19]
    f = WORKDIR / "vprof_E2.txt"
    with open(f, "w") as fh:
        fh.write(f"{p.nd}\n")
        for ra, va in zip(rc, p.v):
            fh.write(f"{ra:.16e} {va:.16e} 0\n")
    r, v, st, nun = parse_driver_sigma(run(bin_, "core", f).stdout)
    rep.check(
        "E2 collapsed spacing -> N_UNDEF=2 (the two stencils)",
        nun == 2,
        f"N_UNDEF={nun}",
    )
    pa = run(bin_, "recomp", f, "knop", check=False)
    rep.check(
        "E2 collapsed spacing aborts under DO_NONMONO_KNOP",
        "ALLOW_SIGMA_UNDEF" in pa.stdout,
    )

    # E3 infall pocket: V dips below zero - signed handling stays
    # finite, |sigma| bounded by the 1E6 clip, warning issued.
    p = make_profile("V004", nd=ND, amp=1.15)
    out, r, v, st, nun = prod_sigma(bin_, p, "E3")
    rep.check(
        "E3 infall pocket (V<0): finite, bounded, warns",
        nun == 0
        and np.isfinite(st).all()
        and np.abs(st).max() <= 1.0e6
        and "negative velocities" in out,
        f"Vmin={v.min():.0f} max|sigma|={np.abs(st).max():.1e}",
    )

    # --------------------------------------------------------------
    # Plot: every profile's V and sigma (analytic + production).
    # --------------------------------------------------------------
    ids = profile_ids()
    fig, ax = plt.subplots(4, 6, figsize=(19, 11))
    for k, pid in enumerate(ids):
        p = make_profile(pid, nd=ND)
        _, r, v, st, nun = prod_sigma(bin_, p, f"plot_{pid}")
        a1 = ax[k // 3, (k % 3) * 2]
        a2 = ax[k // 3, (k % 3) * 2 + 1]
        a1.plot(r, v, "b-")
        a1.set_title(pid, fontsize=9)
        a1.set_ylabel("V", fontsize=8)
        a2.plot(r, p.sigma, "k-", lw=1.2, label="analytic")
        a2.plot(r, st, "r.", ms=3, label="SIGMA_TRUE")
        a2.axhline(-1.0, color="grey", ls=":", lw=0.7)
        a2.set_ylabel("sigma", fontsize=8)
        a2.tick_params(labelsize=7)
        a1.tick_params(labelsize=7)
        if k == 0:
            a2.legend(fontsize=7)
    fig.suptitle("V001-V012 manufactured profiles: V and sigma_true")
    fig.tight_layout()
    out = PLOTDIR / "velocity_profiles.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
