#!/usr/bin/env python3
"""Regression test: short-ray (NI=1,2) stale-GAM fix in SOBEW family.

Bug: in subs/sobew.f, sobew_grad.f, sobew_grad_v2.f and
sobew_grad_v3.f the NI=1 and NI=2 branches never called ZALONGP or
computed GAM, so EXPONX(GAM(I)) on the two outermost tangent rays
used the GAM array left by the previous (NI=3) ray.  Upstream
cmfgen_reference has the identical bug in all four routines, so
linking the driver against jfnk-tree vs reference source gives a
direct fixed-vs-buggy comparison per routine.

Routines exercised through drv_sob_floor:
    sobew   -> SOBEW        (EW diagnostic, production-reachable)
    sobewg1 -> SOBEW_GRAD   (force/EW; no current production caller)
    sobewg2 -> SOBEW_GRAD_V2 (force/EW; called by obs + disp tools)
    sobewg3 -> SOBEW_GRAD_V3 (emiss/EW; no current production caller)

Geometry: NC core rays have NI=ND; shell ray LS>NC sees
NI=ND-(LS-NC-1) grid points and is tangent at P(LS)=R(NI), so with
NP=NC+ND the last two rays have NI=2 and NI=1.

Decisive observable: LINE_FLUX(i) accumulates
0.5*w(LS)*EXPONX(GAM(i,LS))*bracket(i,LS) over rays that reach depth i;
with a thin continuum the per-ray brackets and the per-depth weights
are reproduced exactly on the Python side, so each new/ref LINE_FLUX
ratio is predictable under both current-ray and stale-ray semantics.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from common import (
    BINDIR,
    F90DIR,
    PROD_MOD_INC,
    REF,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link_prod,
    prod_objs,
    run,
)

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

ND, NC, NP = 25, 4, 29  # NP = NC + ND -> last rays are NI=2, NI=1

ROUTINES = ["sobew", "sobewg1", "sobewg2", "sobewg3"]


def exponx(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    s = np.abs(x) < 1.0e-3
    out[s] = 1.0 - x[s] * (0.5 - x[s] / 6.0 * (1.0 - x[s] / 4.0))
    m = (~s) & (x < 40.0)
    out[m] = (1.0 - np.exp(-x[m])) / x[m]
    out[~s & ~m] = 1.0 / x[~s & ~m]
    return out


def geom():
    r = np.geomspace(1.0e5, 1.0e3, ND)  # outer -> inner
    v = np.linspace(8.0e3, 1.0e3, ND)
    p = np.empty(NP)
    p[:NC] = np.linspace(0.3, 0.9, NC) * r[-1]  # core rays, P < R(ND)
    for ls in range(NC, NP):  # shell rays: tangent at R(NI)
        ni = ND - (ls - NC)  # NI = ND+NC+1-LS-1+1
        p[ls] = r[ni - 1]
    return r, v, p


def t1_of(mu2, sigma, floored=True):
    t1 = np.abs(1.0 + sigma * mu2)
    if floored and sigma < 0.0:
        t1 = max(t1, max(-0.2 * (1.0 + sigma), 1.0e-4))
    return t1


def gam_at(i, ls, r, p, gamh, sigma, stale, floored=True):
    """GAM(I,LS) under correct (stale=False) or buggy semantics.

    Buggy semantics: NI=1/2 rays (LS=NP-1,NP) reuse the GAM computed
    on the last NI>2 ray (LS=NP-2).  `floored` selects the jfnk
    (DAH15-floored) vs reference (unclamped ABS) denominator.
    """
    ni = ND - (ls - NC) if ls >= NC else ND
    if stale and ni <= 2:
        ls = NC + ND - 3  # last ray with NI>2
    z2 = r[i] ** 2 - p[ls] ** 2
    t1 = t1_of(z2 / r[i] ** 2, sigma, floored)
    return gamh[i] / t1


def bracket(i, ls, r, p, noes):
    """(T3+T4) factor for depth I on ray LS with a thin continuum.

    T4=0 for core rays (LS<=NC).  AV,TB are O(1e-7) here and dropped.
    """
    ni = ND - (ls - NC) if ls >= NC else ND
    z = np.sqrt(np.maximum(r[: ni + 1] ** 2 - p[ls] ** 2, 0.0))
    dtau = np.empty(ni)
    dtau[0] = noes[0] * r[0]
    for k in range(1, ni):
        dtau[k] = dtau[k - 1] + 0.5 * (noes[k - 1] + noes[k]) * (
            z[k - 1] - z[k]
        )
    t3 = np.exp(-dtau[i])
    t4 = 0.0 if ls < NC else np.exp(-(2.0 * dtau[ni - 1] - dtau[i]))
    return t3 + t4


def s_sum(i, r, p, w, gamh, sigma, stale, noes, floored=None):
    if floored is None:
        floored = not stale  # ref (stale) build lacks the DAH15 floor
    tot = 0.0
    for ls in range(NP):
        ni = ND - (ls - NC) if ls >= NC else ND
        if ni < i + 1:  # ray does not reach depth I
            continue
        tot += (
            w[ls]
            * exponx(
                np.array([gam_at(i, ls, r, p, gamh, sigma, stale, floored)])
            )[0]
            * bracket(i, ls, r, p, noes)
        )
    return float(tot)


def write_case(
    name, sigma, chil0, w_last=None, p_last=None, swap_last=False,
    w_map=None,
):
    r, v, p = geom()
    w = np.full(NP, 1.0 / NP)
    if w_map is not None:
        w[:] = 0.0
        for k, v_ in w_map.items():
            w[k] = v_
    if p_last is not None:
        p[-2], p[-1] = p_last
    if w_last is not None:
        w[-2], w[-1] = w_last
    if swap_last:  # same physical rays + weights, reversed order
        p[-2], p[-1] = p[-1], p[-2]
        w[-2], w[-1] = w[-1], w[-2]
    chi = np.full(ND, 1.0e-12)
    esec = np.full(ND, 1.0e-13)
    chil = np.full(ND, chil0)
    etal = np.abs(chil) * 2.0
    src = np.full(ND, 2.0)
    path = WORKDIR / name
    with open(path, "w") as f:
        f.write(f"{ND} {NC} {NP} 1\n")
        f.write("1.0 F F F 0.0 0.0 1.0e40\n")
        f.write("LINEAR ZERO_FLUX\n")
        for i in range(ND):
            f.write(
                f"{r[i]:.10e} {v[i]:.10e} {sigma:.10e} {chi[i]:.4e} "
                f"{esec[i]:.4e} {chil[i]:.6e} {etal[i]:.6e} {src[i]:.6e}\n"
            )
        for j in range(NP):
            f.write(f"{p[j]:.10e} {w[j]:.6e} {w[j]:.6e}\n")
        f.write(f"{chil0:.6e} {2.0 * abs(chil0):.6e} 1.0\n")
    gamh = chil0 * 3.0e-10 * r / v
    noes = chi - esec
    return path, r, p, w, gamh, noes


def parse(out):
    ew, cf, lf = None, None, {}
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "EW_CF":
            ew = float(t[1])
            cf = float(t[-1])
        elif t[0].lstrip("-").isdigit():
            lf[int(t[0]) - 1] = float(t[-1])
    return ew, cf, np.array([lf[i] for i in sorted(lf)])


def lf_pred(i, cf_printed, r, p, w, gamh, sigma, noes, etal0):
    """Absolute LINE_FLUX(i) prediction for the fixed (current-ray,
    floored) build, using the printed CONT_FLUX for normalization.

    The code scales CONT_FLUX by R(1)^2*13.19868 before dividing
    (sobew.f), so norm = 13.19868*(2.997924e-12/FL^2)/cf_printed;
    case files use FL=1.  s_sum omits the 0.5 prefactor of the
    accumulation, so it is restored here.
    """
    t1norm = 13.19868 * 2.997924e-12 / cf_printed
    s = 0.5 * s_sum(i, r, p, w, gamh, sigma, False, noes, floored=True)
    return t1norm * s * etal0 * r[i] ** 3


def run_rt(bins, rt, case):
    """Return (ew_n, cf_n, lf_n, ew_r, cf_r, lf_r)."""
    ew_n, cf_n, lf_n = parse(run(bins["new"], rt, case).stdout)
    ew_r, cf_r, lf_r = parse(run(bins["ref"], rt, case).stdout)
    return ew_n, cf_n, lf_n, ew_r, cf_r, lf_r


def disc_check(i, bins, rt, case, r, p, w, gamh, sigma, noes,
               ref_usable=True):
    """Stale-vs-current discriminator at depth i.

    Primary: measured new/ref LINE_FLUX(i) ratio vs the predicted
    current-ray/stale-ray ratio.  When the reference build returns a
    nonfinite or zero value - or `ref_usable` is False because the
    reference denominator semantics diverge (upstream sobew_grad uses a
    *signed* (1+sigma*mu^2), corrupted for sigma<-1 independently of the
    stale-ray bug) - fall back to an absolute current-ray prediction
    for the fixed build.
    Returns (ok, measured, predicted, mode).
    """
    ew_n, cf_n, lf_n, ew_r, cf_r, lf_r = run_rt(bins, rt, case)
    s_c = s_sum(i, r, p, w, gamh, sigma, False, noes)
    s_s = s_sum(i, r, p, w, gamh, sigma, True, noes)
    if (
        ref_usable
        and np.isfinite(lf_r[i])
        and lf_r[i] != 0.0
        and np.isfinite(lf_n[i])
    ):
        return (
            abs((lf_n[i] / lf_r[i]) / (s_c / s_s) - 1.0) < 1e-3,
            lf_n[i] / lf_r[i],
            s_c / s_s,
            "ratio",
        )
    etal0 = 2.0e9  # all cases use |CHIL|=1e9 -> ETAL=2|CHIL|
    pred = lf_pred(i, cf_n, r, p, w, gamh, sigma, noes, etal0)
    ok = np.isfinite(lf_n[i]) and abs(lf_n[i] / pred - 1.0) < 1e-3
    return ok, lf_n[i], pred, "abs"


def run_suite(rep, bins, rt, cases):
    """All per-routine checks.  `cases` carries pre-generated case
    files plus the geometry/weight data needed for predictions."""

    sigma = -1.05
    # Upstream sobew_grad (v1) uses a signed (1+sigma*mu^2) denominator,
    # so for sigma < -1 its output is corrupted independently of the
    # stale-ray bug and the new-vs-ref ratio is not interpretable.
    ref_ok = not (rt == "sobewg1" and sigma < -1.0)

    # -- T1 (key regression): only the NI=1 ray has weight, so
    # LINE_FLUX(1) isolates a single EXPONX(GAM) factor: the ratio
    # new/ref is a pure single-ray stale-vs-current discriminator.
    case, r, p, w, gamh, noes = cases["m105_ni1"]
    ok, meas, pred, mode = disc_check(
        0, bins, rt, case, r, p, w, gamh, sigma, noes, ref_usable=ref_ok
    )
    g_c = gam_at(0, NP - 1, r, p, gamh, sigma, False, True)
    g_s = gam_at(0, NP - 1, r, p, gamh, sigma, True, False)
    rep.check(
        f"{rt} T1 NI=1 ray uses current-ray GAM",
        ok,
        f"pred={pred:.4f} meas={meas:.4f} [{mode}] "
        f"(GAM new={g_c:.3f} stale={g_s:.3f})",
    )

    # -- T2: only the NI=2 ray has weight: LINE_FLUX(2) sees it alone.
    case, r, p, w, gamh, noes = cases["m105_ni2"]
    ok, meas, pred, mode = disc_check(
        1, bins, rt, case, r, p, w, gamh, sigma, noes, ref_usable=ref_ok
    )
    g_c2 = gam_at(1, NP - 2, r, p, gamh, sigma, False, True)
    g_s2 = gam_at(1, NP - 2, r, p, gamh, sigma, True, False)
    rep.check(
        f"{rt} T2 NI=2 ray uses current-ray GAM",
        ok,
        f"pred={pred:.4f} meas={meas:.4f} [{mode}] "
        f"(GAM new={g_c2:.3f} stale={g_s2:.3f})",
    )

    # -- T3: for sigma=+2 the floor never engages and the only
    # new-vs-ref difference is the stale-ray fix: depths I>=3 are
    # fed only by NI>2 rays -> bit-identical.
    c3 = cases["pos3"][0]
    _, _, lf3_n, _, _, lf3_r = run_rt(bins, rt, c3)
    ok = np.allclose(lf3_n[2:], lf3_r[2:], rtol=1e-12, atol=0.0)
    rep.check(
        f"{rt} T3 sigma=+2 depths I>=3 identical new vs ref",
        ok,
        f"maxdiff={np.abs(lf3_n[2:] - lf3_r[2:]).max():.2e}",
    )

    # -- T4: ray-ordering — swap the last two rays (distinct P,
    # distinct weights, sparse so a single ray decides).  Depth 2 is
    # reached only by the ray at the NI=2 position: the new/ref
    # LINE_FLUX(2) ratio per case is a pure EXPONX(GAM) ratio whose
    # prediction differs between orderings under both semantics.
    res = {}
    for tag, cs in (("A", cases["ordA"]), ("B", cases["ordB"])):
        res[tag] = disc_check(
            1, bins, rt, cs[0], cs[1], cs[2], cs[3], cs[4], sigma,
            cs[5], ref_usable=ref_ok,
        )
    okA, mA, prA, mdA = res["A"]
    okB, mB, prB, mdB = res["B"]
    ok = okA and okB and abs(prA - prB) > 1e-3 * abs(prA)
    rep.check(
        f"{rt} T4 reordering: per-case check matches current-ray pred",
        ok,
        f"A {mA:.4f}/{prA:.4f} [{mdA}] B {mB:.4f}/{prB:.4f} [{mdB}]",
    )

    # -- T5: monotonic sigma=+2 — fix must still use current ray.
    case, r, p, w, gamh, noes = cases["pos"]
    ok, meas, pred, mode = disc_check(
        0, bins, rt, case, r, p, w, gamh, 2.0, noes
    )
    rep.check(
        f"{rt} T5 sigma=+2 current-ray GAM",
        ok,
        f"pred={pred:.4f} meas={meas:.4f} [{mode}]",
    )

    # -- T6: sigma=-1 edge: same NI=1 discriminator as T1 — the
    # 1e-4 floor must be retained on the short ray (floored new vs
    # unfloored ref prediction).
    case, r, p, w, gamh, noes = cases["m1_ni1"]
    ok, meas, pred, mode = disc_check(
        0, bins, rt, case, r, p, w, gamh, -1.0, noes
    )
    g_c = gam_at(0, NP - 1, r, p, gamh, -1.0, stale=False)
    rep.check(
        f"{rt} T6 sigma=-1 NI=1 ray: floored current-ray GAM",
        ok,
        f"pred={pred:.4f} meas={meas:.4f} [{mode}] "
        f"(GAM(1,last)={g_c:.4f})",
    )

    # -- T7: negative GAMH (inversion), sigma=+2 so floor is irrelevant
    # and only the stale-ray fix differs: sign and EXPONX gain retained.
    case, r, p, w, gamh, noes = cases["neg"]
    ew_n, cf_n, lf_n, ew_r, cf_r, lf_r = run_rt(bins, rt, case)
    s_c = s_sum(0, r, p, w, gamh, 2.0, False, noes, floored=False)
    s_s = s_sum(0, r, p, w, gamh, 2.0, True, noes, floored=False)
    rp = s_c / s_s
    rm = lf_n[0] / lf_r[0]
    ok = np.isfinite(lf_n).all() and abs(rm / rp - 1.0) < 1e-3
    rep.check(
        f"{rt} T7 negative GAMH: gain retained, current-ray GAM",
        ok,
        f"pred={rp:.4f} meas={rm:.4f} EW_new={ew_n:.4e}",
    )

    # -- T8: sigma below -1 (strong reversal regime), NI=1 ray.
    case, r, p, w, gamh, noes = cases["m2_ni1"]
    ok, meas, pred, mode = disc_check(
        0, bins, rt, case, r, p, w, gamh, -2.0, noes,
        ref_usable=(rt != "sobewg1"),
    )
    rep.check(
        f"{rt} T8 sigma=-2 NI=1 ray: floored current-ray GAM",
        ok,
        f"pred={pred:.4f} meas={meas:.4f} [{mode}]",
    )


def main():
    rep = Report("SOBEW short-ray GAM (NI=1,2): SOBEW + GRAD family")
    md = BINDIR / "obj_sob_help"
    help_objs = prod_objs(HELPERS)
    drv = compile_objs(
        "sob_drv",
        [F90DIR / "drv_sob_floor.f90"],
        moddir=md,
        extra_flags=PROD_MOD_INC,
    )
    SOB_FILES = [
        "sobew.f",
        "sobjbar_sim.f",
        "sobew_grad.f",
        "sobew_grad_v2.f",
        "sobew_grad_v3.f",
    ]
    new_sob = prod_objs([TREE / "subs" / f for f in SOB_FILES])
    ref_sob = compile_objs(
        "sob_ref_ewsr",
        [REF / "subs" / f for f in SOB_FILES],
        prod=True,
        moddir=md,
        extra_flags=PROD_MOD_INC,
    )
    bins = {
        "new": link_prod("ewsr_new", help_objs + drv + new_sob, moddir=md),
        "ref": link_prod("ewsr_ref", help_objs + drv + ref_sob, moddir=md),
    }

    r0, v0, p0 = geom()
    cases = {
        "m105_ni1": write_case(
            "sr_m105.txt", -1.05, 1.0e9, w_map={0: 1.0e-3, NP - 1: 1.0}
        ),
        "m105_ni2": write_case(
            "sr_m105b.txt", -1.05, 1.0e9, w_map={0: 1.0e-3, NP - 2: 1.0}
        ),
        "m1_ni1": write_case(
            "sr_m1.txt", -1.0, 1.0e9, w_map={0: 1.0e-3, NP - 1: 1.0}
        ),
        "m2_ni1": write_case(
            "sr_m2.txt", -2.0, 1.0e9, w_map={0: 1.0e-3, NP - 1: 1.0}
        ),
        "pos3": write_case("sr_pos3.txt", 2.0, 1.0e9),
        "ordA": write_case(
            "sr_ordA.txt", -1.05, 1.0e9,
            w_map={0: 1.0e-3, NP - 2: 0.02, NP - 1: 0.08},
            p_last=(0.5 * r0[1], 0.7 * r0[0]),
        ),
        "ordB": write_case(
            "sr_ordB.txt", -1.05, 1.0e9,
            w_map={0: 1.0e-3, NP - 2: 0.02, NP - 1: 0.08},
            p_last=(0.5 * r0[1], 0.7 * r0[0]),
            swap_last=True,
        ),
        "pos": write_case("sr_pos.txt", 2.0, 1.0e9),
        "neg": write_case("sr_neg.txt", 2.0, -1.0e9),
    }

    for rt in ROUTINES:
        run_suite(rep, bins, rt, cases)

    rep.summary()


if __name__ == "__main__":
    main()
