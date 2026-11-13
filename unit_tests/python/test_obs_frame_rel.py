#!/usr/bin/env python3
"""test_obs_frame_rel.py - relativistic transformation audit (P7).

Verifies the complete relativistic chain in OBS_FRAME_SUB_V9:
the zero-velocity limit, the Doppler convention through uniform-
velocity flux predictions, the sign symmetry under +/-V, the
non-monotonic path, the aberrated inner-boundary seed, and the
fail-closed velocity validation.  Runs under OMP_NUM_THREADS=1,2,4,8
for the determinism check.
"""

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from common import BINDIR, F90DIR, TREE, WORKDIR, Report, compile_objs, link, run

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
    F90DIR / "drv_obs_frame_rel.f",
]

JY_CONST = 1.0e23 * 2.0 * np.pi * 1.0e20 / (3.0856e21) ** 2
R1 = 10.0
RND = 2.0
C_KMS = 299792.458
NP = 65


def i_thin_rel(p, beta, j0=1.0, chi=1.0e-6):
    """Boosted thin-shell intensity, identical geometry to the P6
    test: front emitter at z attenuated over [z, zmax], back emitter
    at -z crossing [zmin, z] on its own side plus the front chord.
    chi_obs = chi*X with X = (1-beta*mu)*gam, j_obs = j/X**2."""
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


def disk_integral(i_of_p, npts=40001):
    """int_0^1 I(mu) mu dmu = int_0^R1 I(p) p dp / R1^2."""
    p = np.linspace(0.0, R1, npts)
    return np.trapezoid(i_of_p(p) * p, p) / R1**2


def ib_predict(pgrid, wgrid, nc, beta, aberrated=True, chi=1.0e-6):
    """Predicted OBS_FLUX for the diffusion-boundary manufactured
    model (ETA=0, IB_J=1, IB_H=0.2 constant).  The seed at the
    boundary point is D^3*(J + 3 mu H) attenuated by the near-side
    chord; grazing rays contribute zero.  mu is aberrated when
    aberrated=True, matching the audited convention."""
    gam = 1.0 / np.sqrt(1.0 - beta**2)
    total = 0.0
    for ls in range(nc):
        p = pgrid[ls]
        mu_b = np.sqrt(max(0.0, RND**2 - p**2)) / RND
        mu_c = (mu_b - beta) / (1.0 - beta * mu_b) if aberrated else mu_b
        dop = 1.0 / gam / (1.0 - beta * mu_b)
        zb = np.sqrt(RND**2 - p**2)
        z1 = np.sqrt(R1**2 - p**2)
        z = np.linspace(zb, z1, 2001)
        xf = (1.0 - beta * z / np.hypot(p, z)) * gam
        tau = chi * np.trapezoid(xf, z)
        seed = dop**3 * (1.0 + 3.0 * mu_c * 0.2) * np.exp(-tau)
        total += wgrid[ls] * seed
    return total * JY_CONST * R1**2


def main():
    rep = Report("OBS frame relativistic transformations (P7)")
    out_dir = WORKDIR / "obs_frame_rel"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("flux_*.txt"):
        stale.unlink()
    for stale in out_dir.glob("quad_grid.txt"):
        stale.unlink()

    objs = compile_objs(
        "obs_rel_tree",
        DEPS + [F90DIR / "obs_v9_adapter.f", TREE / "obs" / "obs_frame_sub_v9.f"],
    )
    binary = link("drv_obs_frame_rel", objs, moddir=BINDIR / "obj_obs_rel_tree")

    run(binary, out_dir, "quad", cwd=out_dir)
    prop = run(binary, out_dir, "prop", cwd=out_dir)
    run(binary, out_dir, "run", cwd=out_dir, timeout=900)

    # --- Property checks: closed-form Doppler, round-trip, invariants
    for tag in ("doppler_closed_form", "aberration_roundtrip", "invariant"):
        line = next((ln for ln in prop.stdout.splitlines() if tag in ln), "")
        rep.check(
            f"property check {tag}: zero failures",
            line.endswith("= 0") or "=  0" in line or "=   0" in line,
            line or "marker not found",
        )

    # --- Quadrature grid for the boundary prediction ---------------
    q = np.loadtxt(out_dir / "quad_grid.txt", skiprows=1)
    pgrid, wgrid = q[:, 0], q[:, 2]
    with open(out_dir / "quad_grid.txt") as fh:
        nc = int(fh.readline().split()[1])
    rep.check(
        "quad grid has core rays for the boundary test",
        nc >= 5,
        f"NC={nc}",
    )

    # --- TEST 1: zero velocity --------------------------------------
    rel0 = np.loadtxt(out_dir / "flux_rel_v0.txt")
    norel0 = np.loadtxt(out_dir / "flux_norel_v0.txt")
    rep.check(
        "V=0 with DO_REL=T reproduces DO_REL=F bitwise "
        "(GAM=1 and NU_ON_NU0=1 are exact)",
        np.array_equal(rel0, norel0),
        f"rel={rel0} norel={norel0}",
    )

    # --- TEST 2/4/5: uniform-velocity thin shell --------------------
    for v0, tag in ((20000.0, "u20"), (40000.0, "u40")):
        got = np.loadtxt(out_dir / f"flux_rel_{tag}.txt")
        pred = (
            JY_CONST
            * R1**2
            * disk_integral(lambda p, b=v0 / C_KMS: i_thin_rel(p, b), npts=20001)
        )
        rep.check(
            f"uniform V={int(v0)} km/s thin shell matches the "
            "independent D^2/attenuation prediction",
            np.allclose(got, pred, rtol=3.0e-4, atol=0.0),
            f"got={got} pred={pred}",
        )

    # --- TEST 5: moderately thick slab ------------------------------
    got = np.loadtxt(out_dir / "flux_slab_u20.txt")
    pred = (
        JY_CONST
        * R1**2
        * disk_integral(
            lambda p: i_thin_rel(p, 20000.0 / C_KMS, chi=0.5),
            npts=20001,
        )
    )
    rep.check(
        "uniform V=20000 slab (chi=0.5) matches the prediction with "
        "the contracted path L'=X*L (tolerance is the code's own "
        "ray-grid discretisation at tau~1, the prediction is exact)",
        np.allclose(got, pred, rtol=2.0e-3, atol=0.0),
        f"got={got} pred={pred}",
    )

    # --- TEST 6: mirror pair ----------------------------------------
    mpos = np.loadtxt(out_dir / "flux_mirror_pos.txt")
    mneg = np.loadtxt(out_dir / "flux_mirror_neg.txt")
    rep.check(
        "chi=0 mirror pair (+V / -V with identical |V| and HOLLOW "
        "core) agrees to ulp-level accumulation (a missing mirror "
        "segment would shift the flux by O(1))",
        np.allclose(mpos, mneg, rtol=1.0e-8, atol=0.0),
        f"pos={mpos} neg={mneg} rel={abs(mpos - mneg) / mpos}",
    )

    # --- TEST 7: non-monotonic under rel ----------------------------
    nm = np.loadtxt(out_dir / "flux_rel_nonmono.txt")
    nm0 = np.loadtxt(out_dir / "flux_norel_nonmono.txt")
    rep.check(
        "non-monotonic field (reversal + negative cells) completes "
        "with finite positive flux under DO_REL",
        np.all(np.isfinite(nm)) and np.all(nm > 0.0),
        f"rel={nm}",
    )
    rep.check(
        "non-monotonic rel flux stays near the non-rel value "
        "(boost is a small correction at beta~0.05)",
        np.allclose(nm, nm0, rtol=2.0e-2, atol=0.0),
        f"rel={nm} norel={nm0} diff={nm - nm0}",
    )

    # --- TEST 8: inner-boundary aberration --------------------------
    beta99 = 99.9 / C_KMS
    ib_rel_p = np.loadtxt(out_dir / "flux_ib_rel_p.txt")
    ib_norel_p = np.loadtxt(out_dir / "flux_ib_norel_p.txt")
    ib_rel_m = np.loadtxt(out_dir / "flux_ib_rel_m.txt")
    ib_rel_v0 = np.loadtxt(out_dir / "flux_ib_rel_v0.txt")
    ib_norel_v0 = np.loadtxt(out_dir / "flux_ib_norel_v0.txt")
    ib_etaz_p = np.loadtxt(out_dir / "flux_ib_etaz_p.txt")
    ib_etaz_m = np.loadtxt(out_dir / "flux_ib_etaz_m.txt")

    pred_abr_p = ib_predict(pgrid, wgrid, nc, beta99, aberrated=True)
    pred_raw_p = ib_predict(pgrid, wgrid, nc, beta99, aberrated=False)
    pred_abr_m = ib_predict(pgrid, wgrid, nc, -beta99, aberrated=True)

    rep.check(
        "IB +99.9 km/s boundary matches the ABERRATED prediction "
        "(J + 3*mu_cmf*H evaluated in the comoving frame)",
        np.allclose(ib_rel_p, pred_abr_p, rtol=2.0e-5, atol=0.0),
        f"got={ib_rel_p} pred={pred_abr_p}",
    )
    rep.check(
        "the aberration is measurable: flux differs from the "
        "un-aberrated prescription by the expected O(beta) shift",
        not np.allclose(ib_rel_p, pred_raw_p, rtol=1.0e-6, atol=0.0)
        and abs(ib_rel_p - pred_raw_p)[0] > 0.0,
        f"got={ib_rel_p} unaberrated_pred={pred_raw_p} "
        f"shift={ib_rel_p - pred_raw_p}",
    )
    rep.check(
        "IB -99.9 km/s boundary matches the aberrated prediction "
        "with the negative-beta sign",
        np.allclose(ib_rel_m, pred_abr_m, rtol=2.0e-5, atol=0.0),
        f"got={ib_rel_m} pred={pred_abr_m}",
    )
    rep.check(
        "+V seed exceeds -V seed: the D^3 boost dominates the "
        "aberration at beta~3e-4 (aberration acts on the 3*mu*H "
        "term, D^3 on the whole J-seeded intensity)",
        np.all(ib_rel_p > ib_rel_m),
        f"pos={ib_rel_p} neg={ib_rel_m}",
    )
    rep.check(
        "beta=0 boundary recovers the non-rel result exactly",
        np.array_equal(ib_rel_v0, ib_norel_v0),
        f"rel={ib_rel_v0} norel={ib_norel_v0}",
    )
    rep.check(
        "with DO_REL=F the flat boundary flux is beta-independent "
        "(no aberration path exists outside rel corrections)",
        np.allclose(ib_norel_p, ib_norel_v0, rtol=1.0e-10, atol=0.0),
        f"p={ib_norel_p} v0={ib_norel_v0}",
    )
    rep.check(
        "ETAZ and STAU boundary seeds agree under rel corrections",
        np.allclose(ib_etaz_p, ib_rel_p, rtol=1.0e-4, atol=0.0)
        and np.allclose(ib_etaz_m, ib_rel_m, rtol=1.0e-4, atol=0.0),
        f"stau_p={ib_rel_p} etaz_p={ib_etaz_p} "
        f"stau_m={ib_rel_m} etaz_m={ib_etaz_m}",
    )

    # --- TEST 9: superluminal rejection -----------------------------
    # The literal 299792.458 can sit a few ulps below the fp value of
    # SPEED_OF_LIGHT()*1e-5 used by the code, so the superluminal
    # probes use 299793 - unambiguously >= c under any rounding.
    cval = C_KMS
    for v, name in (
        (299793.0, "v=c+"),
        (1.5 * cval, "v=1.5c"),
        (-299793.0, "v=-c-"),
        (-1.5 * cval, "v=-1.5c"),
    ):
        r = run(binary, out_dir, "sl", repr(v), cwd=out_dir, check=False)
        rep.check(
            f"{name} aborts with SUPERLUMINAL_VELOCITY under DO_REL",
            "SUPERLUMINAL_VELOCITY" in (r.stdout + r.stderr),
            f"rc={r.returncode} out={r.stdout[-120:]}",
        )
    for v, name in ((0.99 * cval, "v=0.99c"), (-0.99 * cval, "v=-0.99c")):
        r = run(binary, out_dir, "sl", repr(v), cwd=out_dir, check=False)
        rep.check(
            f"{name} passes validation and completes",
            "SL completed" in r.stdout and "SUPERLUMINAL" not in (r.stdout + r.stderr),
            f"rc={r.returncode} out={r.stdout[-120:]}",
        )
    r = run(binary, out_dir, "nan", cwd=out_dir, check=False)
    rep.check(
        "NaN velocity aborts with NON_FINITE_VELOCITY",
        "NON_FINITE_VELOCITY" in (r.stdout + r.stderr),
        f"rc={r.returncode} out={r.stdout[-120:]}",
    )

    # --- TEST 10: OpenMP determinism --------------------------------
    omp_dirs = {}
    for nt in (1, 2, 4, 8):
        d = out_dir / f"omp{nt}"
        d.mkdir(exist_ok=True)
        for stale in d.glob("flux_*.txt"):
            stale.unlink()
        env = dict(os.environ, OMP_NUM_THREADS=str(nt))
        r = subprocess.run(
            [str(binary), str(d), "run"],
            cwd=d,
            env=env,
            capture_output=True,
            text=True,
            timeout=900,
        )
        if r.returncode != 0:
            rep.check(f"OMP {nt} threads completes", False, r.stderr[-200:])
            omp_dirs[nt] = None
            continue
        omp_dirs[nt] = d
    if all(omp_dirs.values()):
        base = omp_dirs[1]
        diffs = []
        for nt in (2, 4, 8):
            for f in sorted(base.glob("flux_*.txt")):
                a = np.loadtxt(f)
                b = np.loadtxt(omp_dirs[nt] / f.name)
                diffs.append(np.abs(a - b).max())
        rep.check(
            "OMP determinism: flux outputs identical at " "1/2/4/8 threads",
            max(diffs) == 0.0,
            f"max_abs_diff={max(diffs):.3e}",
        )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
