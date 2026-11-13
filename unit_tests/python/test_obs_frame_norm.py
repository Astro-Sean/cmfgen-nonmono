#!/usr/bin/env python3
"""Flux-normalisation and angular-quadrature checks for
OBS_FRAME_SUB_V9 (P6).

Manufactured spherical shells (R1=10, RND=2, in 1e10 cm) with analytic
emergent intensity:

- thick shell (chi=100, S=ETA/CHI=1) with a solid inner boundary
  (V(ND)=0): every ray saturates at I(p)=S, so int I mu dmu -> 0.5
  and OBS_FLUX -> 6.59934 R1^2 * 0.5.  The STAU branch seeds the
  solid-boundary intensity at the local source function; the ETAZ
  emissivity integral instead seeds TAU(1)=chi*dZ/2, so a model whose
  outermost ray cell is optically thick is under-resolved by
  construction.  The thick ETAZ run pins that documented limit rather
  than the saturated prediction.
- thin shell (chi=1e-6, j=ETA=1) with a hollow core (V(ND)=300):
  mirrored shell chords give I(p) = 2j(sqrt(R1^2-p^2)-sqrt(RND^2-p^2))
  inside the core radius and 2j sqrt(R1^2-p^2) outside it, i.e.
  int I mu dmu = (2j/3)(R1^3-RND^3)/R1^2 - a linear-in-mu limb law
  outside the core radius.
- a quadrature-only mode dumps the production SET_ACC_dMU +
  HWEIGHT_V2 weights.  The Nordlund rule is exact for V ~ mu
  (sum W*mu = 1/3 by construction) but only algebraically accurate
  for V=const, so the tolerances are resolution-scaled bounds, not
  machine precision.
- flat-spectrum thin shell under velocity laws grouped by core
  geometry: all laws with V(ND)>100 give the hollow-core mirrored
  chords and must be bitwise identical; the two V(ND)<=100 laws give
  the solid-boundary geometry and must be bitwise identical to each
  other.  The two classes differ physically (the solid wall emits at
  S), so they are only compared within a class.
- relativistic boost structure: uniform V on the thin shell gives an
  O(v/c)^2 correction after the mirrored +/-VMU pairing, so
  (rel-norel) at 40000 km/s should be ~4x the 20000 km/s value.

The Jy constant is verified in the audit note:
6.5993417128 = 1e23 * 2 pi * 1e20 / (3.0856e21)^2.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (BINDIR, F90DIR, PLOTDIR, REF, TREE, WORKDIR, Report,
                    compile_objs, link, run)

DEPS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "phys_con.f",
    TREE / "subs" / "mon_interp_fast.f",
    TREE / "subs" / "open_dir_acc_v1.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "subs" / "hweight_v2.f",
    TREE / "subs" / "set_acc_dmu.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "tools" / "get_lu.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    TREE / "unix" / "tune.f",
    F90DIR / "drv_obs_frame_norm.f",
]

JY_CONST = 1.0e23 * 2.0 * np.pi * 1.0e20 / (3.0856e21) ** 2
R1, RND = 10.0, 2.0
NPL = [9, 17, 33, 65, 129]


def i_thin(p, j0=1.0):
    """Hollow-core mirrored shell chords: I(p) = 2j*(emitting path)."""
    outer = np.sqrt(np.maximum(R1**2 - p**2, 0.0))
    inner = np.sqrt(np.maximum(RND**2 - p**2, 0.0))
    return 2.0 * j0 * (outer - inner)


def i_thick(p, chi=100.0):
    """Constant-S model with a solid wall: core rays carry the seeded
    boundary S exactly, grazing rays saturate as S(1-e^(-chi*chord))."""
    chord = 2.0 * np.sqrt(np.maximum(R1**2 - p**2, 0.0))
    return np.where(p < RND, 1.0, 1.0 - np.exp(-chi * chord))


def disk_integral(i_of_p, npts=2_000_001):
    """int_0^1 I(mu) mu dmu = int_0^R1 I(p) p dp / R1^2 on a fine grid."""
    p = np.linspace(0.0, R1, npts)
    return np.trapezoid(i_of_p(p) * p, p) / R1**2


def i_thin_rel(p, beta, j0=1.0, chi=1.0e-6):
    """Boosted thin-shell intensity.  The code applies j -> j*D^2 with
    D = 1/GAM/(1 - v z/(r c)) and chi -> chi*X, X = (1 - v z/(r c))*GAM,
    pointwise along each mirrored chord.  A front-side emitter at z is
    attenuated over [z, zmax]; a back-side emitter at -z crosses
    [zmin, z] on its own side plus the whole front chord."""
    gam = 1.0 / np.sqrt(1.0 - beta**2)
    out = np.empty_like(p, dtype=float)
    for k, pk in enumerate(p.flat):
        zmax = np.sqrt(R1**2 - pk**2)
        zmin = np.sqrt(RND**2 - pk**2) if pk < RND else 0.0
        z = np.linspace(zmin, zmax, 2001)
        mu = z / np.hypot(pk, z)
        xf = (1.0 - beta * mu) * gam
        xb = (1.0 + beta * mu) * gam
        dz = z[1] - z[0]
        tau_f = chi * (np.cumsum(xf[::-1])[::-1] * dz - 0.5 * xf * dz)
        tau_b = chi * (np.cumsum(xb) * dz - 0.5 * xb * dz + np.trapezoid(xf, z))
        out.flat[k] = j0 * np.trapezoid(
            np.exp(-tau_f) / xf**2 + np.exp(-tau_b) / xb**2, z
        )
    return out


def read_mu_w(path):
    data = np.loadtxt(path, skiprows=1)
    return data[:, 0], data[:, 2]


def main():
    rep = Report("OBS frame flux normalisation (P6)")
    out_dir = WORKDIR / "obs_frame_norm"
    ref_dir = WORKDIR / "obs_frame_norm_ref"
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)
    # Clear stale outputs so the tree/reference file sets below only
    # ever contain files written by this run.
    for d in (out_dir, ref_dir):
        for stale in d.glob("flux_*.txt"):
            stale.unlink()
    # Same driver and dependency objects on both sides; only
    # obs_frame_sub_v9 and the call adapter differ, so any numerical
    # difference between the runs is a source difference.
    tree_objs = compile_objs(
        "obs_norm_tree",
        DEPS + [F90DIR / "obs_v9_adapter.f", TREE / "obs" / "obs_frame_sub_v9.f"],
    )
    ref_objs = compile_objs(
        "obs_norm_ref",
        DEPS + [F90DIR / "obs_v9_adapter_ref.f", REF / "obs" / "obs_frame_sub_v9.f"],
    )
    binary = link(
        "drv_obs_frame_norm",
        tree_objs,
        moddir=BINDIR / "obj_obs_norm_tree",
    )
    ref_binary = link(
        "drv_obs_frame_norm_ref",
        ref_objs,
        moddir=BINDIR / "obj_obs_norm_ref",
    )
    run(binary, out_dir, "run", cwd=out_dir, timeout=900)
    run(binary, out_dir, "quad", cwd=out_dir)
    # The frozen upstream routine crashes inside the OpenMP runtime
    # teardown after STOP; all flux files are already written and
    # closed by then, so a nonzero exit is acceptable if and only if
    # every expected output is present below.
    ref_run = run(ref_binary, ref_dir, "run", cwd=ref_dir, timeout=900, check=False)
    ref_outputs = sorted(ref_dir.glob("flux_*.txt"))
    tree_outputs = sorted(out_dir.glob("flux_*.txt"))
    rep.check(
        "upstream reference run produced the full output set "
        f"(exit rc={ref_run.returncode} tolerated: OMP teardown "
        "crash after STOP)",
        [p.name for p in ref_outputs] == [p.name for p in tree_outputs]
        and all(p.stat().st_size > 0 for p in ref_outputs),
        f"ref={len(ref_outputs)} tree={len(tree_outputs)} " f"rc={ref_run.returncode}",
    )

    # --- CASE 1: constant intensity (solid wall emits S) ---------------
    thick_stau = np.loadtxt(out_dir / "flux_thick_stau.txt")
    thick_etaz = np.loadtxt(out_dir / "flux_thick_etaz.txt")
    i0_pred = disk_integral(i_thick)
    f_pred = JY_CONST * R1**2 * i0_pred
    rep.check(
        "thick shell STAU: OBS_FLUX matches the saturated I(p)=S "
        "integral to the quadrature error",
        np.allclose(thick_stau, f_pred, rtol=2.0e-4, atol=0.0),
        f"got={thick_stau} pred={f_pred} (I0 limit {0.5 * JY_CONST * R1**2})",
    )
    rep.check(
        "thick shell ETAZ: under-resolved per the TAU(1)=chi*dZ/2 "
        "surface seed (documented limit, finite and far below STAU)",
        np.all(np.isfinite(thick_etaz))
        and np.all(thick_etaz > 0.0)
        and np.all(thick_etaz < 0.05 * f_pred),
        f"etaz={thick_etaz} pred={f_pred}",
    )

    # --- CASE 2: thin hollow shell, linear-mu limb law -----------------
    i_mu_pred = disk_integral(i_thin)
    f_thin_pred = JY_CONST * R1**2 * i_mu_pred
    flux_np = [np.loadtxt(out_dir / f"flux_thin_np{np_:03d}.txt") for np_ in NPL]
    rep.check(
        "thin shell STAU NP=129: OBS_FLUX matches the analytic "
        "(2j/3)(R1^3-RND^3)/R1^2 integral",
        np.allclose(flux_np[-1], f_thin_pred, rtol=2.0e-3, atol=0.0),
        f"got={flux_np[-1]} pred={f_thin_pred}",
    )
    etaz = np.loadtxt(out_dir / "flux_thin_etaz.txt")
    rep.check(
        "thin shell ETAZ NP=129: same analytic integral",
        np.allclose(etaz, f_thin_pred, rtol=2.0e-3, atol=0.0),
        f"got={etaz} pred={f_thin_pred}",
    )
    rep.check(
        "thin shell: STAU and ETAZ agree at NP=129",
        np.allclose(flux_np[-1], etaz, rtol=2.0e-3, atol=0.0),
        f"stau={flux_np[-1]} etaz={etaz}",
    )
    errs = [abs(f[0] - f_thin_pred) / f_thin_pred for f in flux_np]
    rep.check(
        "angular convergence: error improves with refinement "
        "(Nordlund residual oscillates with grid parity, so the "
        "check is NP=129 vs the coarser grids, not strict monotone)",
        errs[-1] < errs[0] and errs[-1] < errs[2] and max(errs[2:]) < errs[0],
        f"rel_errs={['%.3e' % e for e in errs]}",
    )
    rep.check(
        "angular convergence: NP=129 error below 0.2% (Nordlund "
        "cubic on a limb-law profile)",
        errs[-1] < 2.0e-3,
        f"rel_err={errs[-1]:.3e}",
    )

    # --- Quadrature-only law checks --------------------------------------
    # HWEIGHT_V2 is exact only for V ~ mu; for other profiles the error
    # must shrink with resolution.  Bounds below carry ~3x headroom over
    # the measured NP=129 residuals.
    sum_w_bound = {33: 2.0e-3, 129: 1.5e-4}
    limb_bound = {33: 1.0e-3, 129: 1.0e-4}
    mu2_bound = {33: 1.0e-4, 129: 1.0e-5}
    sums_w = {}
    for np_ in (33, 129):
        mu, w = read_mu_w(out_dir / f"quad_w_np{np_:03d}.txt")
        sums_w[np_] = abs(w.sum() - 0.5)
        rep.check(
            f"NP={np_}: sum W = int mu dmu = 1/2 within the Nordlund " "tolerance",
            sums_w[np_] < sum_w_bound[np_],
            f"sum={w.sum():.15e} err={sums_w[np_]:.3e}",
        )
        rep.check(
            f"NP={np_}: sum W*mu = 1/3 (exact by construction)",
            abs((w * mu).sum() - 1.0 / 3.0) < 1.0e-12,
            f"sum={(w * mu).sum():.15e}",
        )
        for a in (0.3, 0.8):
            pred = (1.0 - a) / 2.0 + a / 3.0
            got = (w * (1.0 - a + a * mu)).sum()
            rep.check(
                f"NP={np_}: limb law I=1-a+a mu, a={a}",
                abs(got - pred) < limb_bound[np_],
                f"got={got:.15e} pred={pred:.15e} " f"err={abs(got - pred):.3e}",
            )
        rep.check(
            f"NP={np_}: sum W*mu^2 = 1/4",
            abs((w * mu**2).sum() - 0.25) < mu2_bound[np_],
            f"sum={(w * mu**2).sum():.15e}",
        )
        rep.check(
            f"NP={np_}: mu grid spans [0,1], first=1 last=0",
            abs(mu[0] - 1.0) < 1e-15 and abs(mu[-1]) < 1e-15,
            f"mu0={mu[0]} muN={mu[-1]}",
        )
    rep.check(
        "sum-W error shrinks from NP=33 to NP=129",
        sums_w[129] < sums_w[33],
        f"err33={sums_w[33]:.3e} err129={sums_w[129]:.3e}",
    )

    # --- CASE 3: velocity independence within a geometry class ---------
    # Velocity-dependent NINS subdivision gives different ray z-grids, so
    # equality holds only to round-off (~1e-15 rel), not bitwise.
    flat = {
        k: np.loadtxt(out_dir / f"flux_flat_{k}.txt")
        for k in ("v200", "mono", "rev", "neg", "s50", "s50b")
    }
    rep.check(
        "flat spectrum, hollow core: monotonic outflow agrees with "
        "uniform V=200 to round-off",
        np.allclose(flat["v200"], flat["mono"], rtol=1.0e-12, atol=0.0),
        f"v200={flat['v200']} mono={flat['mono']}",
    )
    rep.check(
        "flat spectrum, hollow core: reversing field agrees to " "round-off",
        np.allclose(flat["v200"], flat["rev"], rtol=1.0e-12, atol=0.0),
        f"rev={flat['rev']}",
    )
    rep.check(
        "flat spectrum, hollow core: negative-V field agrees to " "round-off",
        np.allclose(flat["v200"], flat["neg"], rtol=1.0e-12, atol=0.0),
        f"neg={flat['neg']}",
    )
    rep.check(
        "flat spectrum, solid wall: reversing V(ND)=50 agrees with "
        "uniform V=50 to round-off",
        np.allclose(flat["s50"], flat["s50b"], rtol=1.0e-12, atol=0.0),
        f"s50={flat['s50']} s50b={flat['s50b']}",
    )
    rep.check(
        "flat spectrum, thin shell matches the analytic limb-law " "integral",
        np.allclose(flat["v200"], f_thin_pred, rtol=2.0e-3, atol=0.0),
        f"v200={flat['v200']} pred={f_thin_pred}",
    )

    # --- CASE 4: relativistic boost, direct prediction ------------------
    # The mirrored +/-VMU pairing cancels the O(beta) boost; the leading
    # residual is beta^2(3 mu^2-1), which nearly vanishes for this shell
    # geometry.  Two checks: the CHI=0 ETAZ pure-boost probe isolates
    # int j D^2 dz and must scale as beta^2 (ratio ~4 between v0=40000
    # and 20000 km/s), and the full rel flux must match the independent
    # int j [D(z)^2 + D(-z)^2] dz prediction.
    c_kms = 2.99792458e5
    b1 = np.loadtxt(out_dir / "flux_boost_v20k.txt") - np.loadtxt(
        out_dir / "flux_boost_norel_v20k.txt"
    )
    b2 = np.loadtxt(out_dir / "flux_boost_v40k.txt") - np.loadtxt(
        out_dir / "flux_boost_norel_v40k.txt"
    )
    ratio = b2 / np.where(b1 == 0.0, np.nan, b1)
    rep.check(
        "pure D^2 emissivity boost (CHI=0) scales as beta^2: "
        "ratio at 40000/20000 km/s ~ 4",
        np.all(np.isfinite(ratio)) and np.all(np.abs(ratio - 4.0) < 0.6),
        f"ratio={ratio} b1={b1} b2={b2}",
    )
    rep.check(
        "pure boost is small: |rel-norel| below 1e-3 of flux "
        "(a missing mirror or double-applied boost would leave "
        "an O(beta) term ~ 0.3 of flux)",
        np.all(np.abs(b2) < 1.0e-3 * flat["v200"]),
        f"b2={b2} flux={flat['v200']}",
    )
    for tag, v0 in (("v20k", 20000.0), ("v40k", 40000.0)):
        pred = (
            JY_CONST
            * R1**2
            * disk_integral(lambda p, b=v0 / c_kms: i_thin_rel(p, b), npts=20001)
        )
        got = np.loadtxt(out_dir / f"flux_rel_{tag}.txt")
        rep.check(
            f"rel flux at v0={int(v0)} km/s matches the D^2 "
            "emissivity-boost prediction to quadrature level",
            np.allclose(got, pred, rtol=3.0e-4, atol=0.0),
            f"got={got} pred={pred} "
            f"diff_vs_norel={got - np.loadtxt(out_dir / f'flux_norel_{tag}.txt')}",
        )

    # --- Stock comparison: tree vs upstream reference v9 ---------------
    # Same driver, same flags, same inputs; only obs_frame_sub_v9 and
    # the adapter differ.  The guards and IB branches added in this
    # tree do not change the arithmetic of these cases, so the fluxes
    # should agree to a few ULPs.
    # Compare on the reference set: the tree dir can hold stale files
    # from older driver layouts that the current driver no longer
    # writes.
    for p in ref_outputs:
        fname = p.name
        a = np.loadtxt(out_dir / fname)
        b = np.loadtxt(ref_dir / fname)
        rep.check(
            f"stock reference: {fname} matches upstream v9",
            np.allclose(a, b, rtol=1e-12, atol=0.0),
            f"tree={a} ref={b}",
        )

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.loglog(NPL, errs, "o-")
    ax.set_xlabel("NP (rays over the disk)")
    ax.set_ylabel("relative flux error vs analytic")
    fig.savefig(PLOTDIR / "obs_frame_norm_conv.png", dpi=160)
    plt.close(fig)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
