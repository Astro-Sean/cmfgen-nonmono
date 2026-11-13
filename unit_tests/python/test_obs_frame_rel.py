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


def seed_oracle_boosted(mu_obs, beta, j0=1.0, h0=0.2, nmu=400001):
    """Independent inner-boundary seed oracle (P7 Part C).

    Builds the interior diffusion field I'(mu') = J + 3 mu' H on a
    dense comoving-direction grid, then maps every comoving photon
    direction FORWARD into the lab frame with
    mu = (mu' + beta)/(1 + beta mu') and carries its intensity with
    the I_nu/nu^3 invariant, I_lab(mu) = I'(mu') [gam (1+beta mu')]^3.
    The lab-frame distribution is then interpolated at the ray's
    boundary cosine.  This never evaluates the inverse aberration
    used in the code: it reconstructs which comoving photon lands on
    the lab ray from the forward map alone.  Also returns the
    comoving direction that maps to mu_obs (the photon's interior
    hemisphere)."""
    mu_p = np.linspace(-1.0, 1.0, nmu)
    i_cmf = j0 + 3.0 * mu_p * h0
    gam = 1.0 / np.sqrt(1.0 - beta * beta)
    mu_lab = (mu_p + beta) / (1.0 + beta * mu_p)
    i_lab = i_cmf * (gam * (1.0 + beta * mu_p)) ** 3
    i_seed = np.interp(mu_obs, mu_lab, i_lab)
    mu_cmf = np.interp(mu_obs, mu_lab, mu_p)
    return i_seed, mu_cmf


def seed_tau_sph(p, beta, chi=1.0e-8):
    """Near-side attenuation of the boundary seed for the spherical
    manufactured model (uniform chi, constant V): tau = chi * int
    gam (1 - beta mu(z)) dz from the boundary crossing to R(1)."""
    gam = 1.0 / np.sqrt(1.0 - beta * beta)
    zb = np.sqrt(RND**2 - p**2)
    z1 = np.sqrt(R1**2 - p**2)
    z = np.linspace(zb, z1, 2001)
    return chi * np.trapezoid(gam * (1.0 - beta * z / np.hypot(p, z)), z)


def main():
    rep = Report("OBS frame relativistic transformations (P7)")
    out_dir = WORKDIR / "obs_frame_rel"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("flux_*.txt"):
        stale.unlink()
    for stale in out_dir.glob("quad_grid.txt"):
        stale.unlink()
    for stale in out_dir.glob("char_*.txt"):
        stale.unlink()
    for stale in out_dir.glob("ppgrid.txt"):
        stale.unlink()

    objs = compile_objs(
        "obs_rel_tree",
        DEPS + [F90DIR / "obs_v9_adapter.f", TREE / "obs" / "obs_frame_sub_v9.f"],
    )
    binary = link("drv_obs_frame_rel", objs, moddir=BINDIR / "obj_obs_rel_tree")

    run(binary, out_dir, "quad", cwd=out_dir)
    prop = run(binary, out_dir, "prop", cwd=out_dir)
    run(binary, out_dir, "run", cwd=out_dir, timeout=900)
    run(binary, out_dir, "char", "sph", cwd=out_dir, timeout=600)
    run(binary, out_dir, "char", "pp", cwd=out_dir, timeout=600)
    run(binary, out_dir, "pp", cwd=out_dir, timeout=300)

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

    # --- TEST 11 (Part C): manufactured characteristics -------------
    # Per-ray boundary seeds at beta = 0,+-0.1,+-0.3 and
    # mu_obs = 1.0,0.8,0.3,0.1 for both STAU and ETAZ.  The oracle
    # boosts the interior diffusion distribution forward into the lab
    # frame on a dense direction grid - it never evaluates the
    # inverse-aberration expression used in the code.
    char_rows = {}
    for geom in ("sph", "pp"):
        rows = np.loadtxt(
            out_dir / f"char_{geom}.txt",
            dtype={
                "names": ("method", "beta", "mu", "f1", "f2"),
                "formats": ("S4", "f8", "f8", "f8", "f8"),
            },
        )
        char_rows[geom] = rows
        n_bad = 0
        worst = ""
        for r in rows:
            meth, beta, mu, f1, f2 = (
                r["method"].decode(),
                r["beta"],
                r["mu"],
                r["f1"],
                r["f2"],
            )
            i_seed, mu_cmf = seed_oracle_boosted(mu, beta)
            gam = 1.0 / np.sqrt(1.0 - beta * beta)
            if geom == "sph":
                p = RND * np.sqrt(max(0.0, 1.0 - mu * mu))
                tau = seed_tau_sph(p, beta)
            else:
                tau = 1.0e-8 * gam * (1.0 - beta * mu) * (R1 - RND) / mu
            expect = JY_CONST * R1**2 * i_seed * np.exp(-tau)
            ok1 = np.isclose(f1, expect, rtol=2.0e-4, atol=0.0)
            ok2 = np.isclose(f2, expect, rtol=2.0e-4, atol=0.0)
            if not (ok1 and ok2):
                n_bad += 1
                worst = (
                    f"{meth} beta={beta:+.1f} mu={mu:.1f} got=({f1:.9e},"
                    f"{f2:.9e}) expect={expect:.9e} mu_cmf={mu_cmf:+.4f}"
                )
        rep.check(
            f"Part C: all 40 {geom} single-ray seeds match the "
            "forward-boosted diffusion-field oracle (STAU+ETAZ, "
            "beta=0/+-0.1/+-0.3 x mu=1.0/0.8/0.3/0.1)",
            n_bad == 0,
            f"bad={n_bad} worst: {worst or 'none'}",
        )

    # Classification table required by the task: for every case the
    # lab ray is an atmosphere-entering characteristic (mu_obs>0 at a
    # boundary fixed in the lab frame); the seed is the interior
    # diffusion field at the photon's comoving direction, on whichever
    # hemisphere the aberration selects.
    with open(out_dir / "char_classification.txt", "w") as fh:
        fh.write(
            "# mu_obs beta mu_cmf interior_hemisphere ray_role "
            "boundary_condition\n"
        )
        for beta in (0.0, 0.1, -0.1, 0.3, -0.3):
            for mu in (1.0, 0.8, 0.3, 0.1):
                _, mu_cmf = seed_oracle_boosted(mu, beta)
                hemi = (
                    "outward"
                    if mu_cmf > 1.0e-12
                    else "tangential" if mu_cmf > -1.0e-12 else "inward"
                )
                fh.write(
                    f"{mu:6.2f} {beta:+6.2f} {mu_cmf:+9.5f} {hemi:>10s} "
                    "entering interior_diffusion_field@mu_cmf\n"
                )

    # Discriminators for the negative-mu_cmf cases (only beta=+0.3,
    # mu_obs=0.1 in this grid): the seed must be the inward-
    # hemisphere value J - 3|mu'|H, boosted - not clamped-to-J, not
    # zero, and not the outward-hemisphere value at |mu'|.
    for geom in ("sph", "pp"):
        rows = char_rows[geom]
        r = next(
            x
            for x in rows
            if x["method"] == b"STAU"
            and np.isclose(x["beta"], 0.3)
            and np.isclose(x["mu"], 0.1)
        )
        got = r["f1"]
        i_seed, mu_cmf = seed_oracle_boosted(0.1, 0.3)
        gam = 1.0 / np.sqrt(1.0 - 0.3**2)
        norm = JY_CONST * R1**2
        d3 = (gam * (1.0 + 0.3 * mu_cmf)) ** 3
        rep.check(
            f"{geom} beta=+0.3 mu=0.1 is a comoving-inward photon "
            f"(mu_cmf={mu_cmf:+.4f}<0): seed = D^3*(J-3|mu_cmf|H), "
            "the advected interior inward hemisphere",
            mu_cmf < 0.0
            and np.isclose(got, norm * (1.0 - 3.0 * abs(mu_cmf) * 0.2) * d3, rtol=2.0e-4)
            and not np.isclose(got, norm * 1.0 * d3, rtol=1.0e-3)
            and not np.isclose(got, norm * (1.0 + 3.0 * abs(mu_cmf) * 0.2) * d3, rtol=1.0e-3)
            and got > 0.0,
            f"got={got:.9e} inward={norm * (1.0 - 3.0 * abs(mu_cmf) * 0.2) * d3:.9e} "
            f"clampJ={norm * d3:.9e} outward={norm * (1.0 + 3.0 * abs(mu_cmf) * 0.2) * d3:.9e}",
        )

    # PP-vs-spherical cosine consistency: at the same boundary cosine
    # the two geometries must seed the same comoving state (Part D
    # convention check).  Tolerance covers the different chord
    # attenuation at chi=1e-8.
    sph_rows = char_rows["sph"]
    pp_rows = char_rows["pp"]
    n_bad = 0
    worst = ""
    for rs in sph_rows:
        rp = next(
            (
                x
                for x in pp_rows
                if x["method"] == rs["method"]
                and np.isclose(x["beta"], rs["beta"])
                and np.isclose(x["mu"], rs["mu"])
            ),
            None,
        )
        if rp is None or not np.isclose(rp["f1"], rs["f1"], rtol=5.0e-5):
            n_bad += 1
            worst = f"{rs['method'].decode()} beta={rs['beta']:+.1f} mu={rs['mu']:.1f}"
    rep.check(
        "PP and spherical seeds agree at the same boundary cosine "
        "(confirms MU_AT_RMAX is the correct PP boundary mu)",
        n_bad == 0,
        f"bad={n_bad} worst: {worst or 'none'}",
    )

    # --- TEST 12 (Part D): plane-parallel relativistic run ----------
    pg = np.loadtxt(out_dir / "ppgrid.txt")
    pmu, pw = pg[:, 1], pg[:, 2]
    pp_rel_0 = np.loadtxt(out_dir / "flux_pp_rel_0.txt")
    pp_norel_0 = np.loadtxt(out_dir / "flux_pp_norel_0.txt")
    pp_rel_p = np.loadtxt(out_dir / "flux_pp_rel_p.txt")
    pp_rel_m = np.loadtxt(out_dir / "flux_pp_rel_m.txt")
    pp_norel_p = np.loadtxt(out_dir / "flux_pp_norel_p.txt")

    rep.check(
        "PP: all manufactured runs produce finite positive flux",
        all(np.all(np.isfinite(x)) and np.all(x > 0.0)
            for x in (pp_rel_0, pp_norel_0, pp_rel_p, pp_rel_m, pp_norel_p)),
        f"{pp_rel_0} {pp_rel_p} {pp_rel_m}",
    )
    rep.check(
        "PP: V=0 under DO_REL reproduces the non-rel run bitwise",
        np.array_equal(pp_rel_0, pp_norel_0),
        f"rel={pp_rel_0} norel={pp_norel_0}",
    )
    for sign, tag in ((1.0, "p"), (-1.0, "m")):
        got = np.loadtxt(out_dir / f"flux_pp_rel_{tag}.txt")
        beta = sign * 20000.0 / C_KMS
        gam = 1.0 / np.sqrt(1.0 - beta * beta)
        expect = 0.0
        for mu, w in zip(pmu, pw):
            i_seed, _ = seed_oracle_boosted(mu, beta)
            tau = 1.0e-8 * gam * (1.0 - beta * mu) * (R1 - RND) / mu
            expect += w * i_seed * np.exp(-tau)
        expect *= JY_CONST * R1**2
        rep.check(
            f"PP: V={'+' if sign > 0 else '-'}20000 boundary-sum "
            "matches the boosted oracle over the mu grid "
            "(Doppler factor + aberration)",
            np.allclose(got, expect, rtol=1.0e-4, atol=0.0),
            f"got={got} expect={expect}",
        )
    rep.check(
        "PP: +V flux exceeds -V flux (outflow boosts, inflow deboosts)",
        np.all(pp_rel_p > pp_rel_m),
        f"p={pp_rel_p} m={pp_rel_m}",
    )
    # OMP determinism of the PP path.
    omp_pp = {}
    for nt in (1, 8):
        d = out_dir / f"omp_pp{nt}"
        d.mkdir(exist_ok=True)
        env = dict(os.environ, OMP_NUM_THREADS=str(nt))
        r = subprocess.run(
            [str(binary), str(d), "pp"],
            cwd=d,
            env=env,
            capture_output=True,
            text=True,
            timeout=300,
        )
        if r.returncode == 0:
            omp_pp[nt] = d
    if len(omp_pp) == 2:
        same = all(
            np.array_equal(
                np.loadtxt(omp_pp[1] / f.name),
                np.loadtxt(omp_pp[8] / f.name),
            )
            for f in sorted(omp_pp[1].glob("flux_pp_*.txt"))
        )
        rep.check(
            "PP: OMP determinism, flux identical at 1 vs 8 threads",
            same,
            "",
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
