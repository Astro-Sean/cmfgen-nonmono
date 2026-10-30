#!/usr/bin/env python3
"""Runtime verification of MOM_J_CMF_NONMONO_V1's discrete operator.

drv_mom_solve.f runs the production moment solver on manufactured
atmospheres whose sigma_true profiles put controlled compression bands
on the expanded grid, then dumps every internal coefficient array.
This test recomputes each coefficient from first principles (the Knop
et al. 2009 Paper I split, V11 moment structure) and compares against
the production values; it also checks field admissibility, the
tridiagonal residual, the lagged fixed-point drift, and the zero-GAM
resonance behaviour (MOM_GAM_FAC_MIN floor + exact-zero exception).
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    F90DIR,
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
    "unix/gen_asci_open.f",
    "unix/verbose_info.f",
    "subs/set_case.f",
    "subs/phys_con.f",
    "subs/mon_int_funs_v2.f",
    "subs/mon_interp.f",
    "subs/mod_trap_derivatives.f",
    "subs/nordtau.f",
    "subs/qfromf.f",
    "subs/integrate.f",
    "subs/write_vec.f",
    "unix/dir_acc_pars_gen.f",
    "plane/mod_ray_mom_store.f",
    "new_main/mod_subs/knop_status_mod.f90",
    "plane/get_moms_non_rel.f",
    "new_main/mod_subs/mod_cmfgen.f",
    "new_main/mod_subs/control_variable_mod.f",
    "new_main/mod_subs/jfnk_cmfgen_int_mod.f",
    "new_main/mod_subs/jfnk_solver_v1.f90",
    "new_main/mod_subs/check_param_consistency.f",
    "newsubs/rd_store_log.f",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
]

CGAM = 3.33564e-6  # 1/c in the code's km/s/cm units
MOM_GAM_FAC_MIN = 1.0e-2
MOM_DLN_MIN = 1.0e-6


def build():
    objs = compile_objs(
        "mom_solve",
        [TREE / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_mom_solve.f"],
    )
    return link("mom_solve", objs)


def load(tag):
    d = np.loadtxt(WORKDIR / f"{tag}_nm.txt")
    return {
        "r": d[:, 1], "v": d[:, 2], "sigt": d[:, 3], "konj": d[:, 4],
        "nmj": d[:, 5], "nmh": d[:, 6], "chi": d[:, 7], "eta": d[:, 8],
        "esec": d[:, 9], "chihat": d[:, 10], "gams": d[:, 11],
        "q": d[:, 12], "w": d[:, 13], "eps": d[:, 14], "psi": d[:, 15],
        "psiprev": d[:, 16], "hu": d[:, 17], "hl": d[:, 18],
        "hs": d[:, 19], "ta": d[:, 20], "tb": d[:, 21], "tc": d[:, 22],
        "dtau": d[:, 23], "dtauonq": d[:, 24], "qf": d[:, 25],
        "rsqj": d[:, 26], "rsqh": d[:, 27], "xmsave": d[:, 28],
        "rsqjp": d[:, 29], "rsqhp": d[:, 30],
    }


def gam_fac(sigt, k):
    """Production factor: SIGN(MAX(|T1|,MIN),T1); T1==0 -> 0."""
    t1 = 1.0 + sigt * k
    f = np.sign(t1) * np.maximum(np.abs(t1), MOM_GAM_FAC_MIN)
    f[t1 == 0.0] = 0.0
    return f


def check_coefficients(rep, name, tag, dln=0.05):
    d = load(tag)
    n = len(d["r"])
    konj = d["konj"]

    fac = gam_fac(d["sigt"], konj)
    con_gam = CGAM * d["v"] / d["r"]
    gams_ref = con_gam * fac
    rep.check(
        f"{name}: signed GAM = CON_GAM*GAM_FAC",
        np.allclose(d["gams"], gams_ref, rtol=1e-10, atol=1e-14),
        f"max diff {np.abs(d['gams']-gams_ref).max():.2e}",
    )
    chihat_ref = d["chi"] + np.abs(gams_ref) / max(dln, MOM_DLN_MIN)
    rep.check(
        f"{name}: chi_hat = chi + |GAM|/DLN (Knop Eq 7)",
        np.allclose(d["chihat"], chihat_ref, rtol=1e-10, atol=1e-14),
        f"max diff {np.abs(d['chihat']-chihat_ref).max():.2e}",
    )
    gam_ref = con_gam / chihat_ref / dln
    q_ref = gam_ref * fac
    rep.check(
        f"{name}: Q_NM = (CON_GAM/chi_hat/DLN)*GAM_FAC (signed)",
        np.allclose(d["q"], q_ref, rtol=1e-9, atol=1e-14),
        f"max diff {np.abs(d['q']-q_ref).max():.2e}",
    )
    psi_ref = d["dtauonq"] * np.abs(q_ref)
    i2 = slice(1, n - 1)
    rep.check(
        f"{name}: PSI = DTAUONQ*|Q| >= 0 implicit diagonal",
        np.allclose(d["psi"][i2], psi_ref[i2], rtol=1e-9, atol=1e-14),
        f"max diff {np.abs(d['psi'][i2]-psi_ref[i2]).max():.2e}",
    )
    rep.check(
        f"{name}: PSI/PSIPREV non-negative everywhere",
        bool((d["psi"] >= 0).all() and (d["psiprev"] >= 0).all()),
    )
    rep.check(
        f"{name}: MOM_GAM_PIN -> PSIPREV == PSI on interior",
        np.allclose(d["psiprev"][i2], d["psi"][i2], rtol=1e-12, atol=0),
        f"max diff {np.abs(d['psiprev'][i2]-d['psi'][i2]).max():.2e}",
    )

    # W = GAMH*|GAMH_FAC|, GAMH = CON_GAMH/DLN/(CHIHAT_i+CHIHAT_{i+1})
    av_sig = 0.5 * (d["sigt"][:-1] + d["sigt"][1:])
    gfac_h = gam_fac(av_sig, d["nmh"][:-1])
    con_gamh = 2.0 * CGAM * (d["v"][:-1] + d["v"][1:]) / (
        d["r"][:-1] + d["r"][1:]
    )
    gamh_ref = con_gamh / dln / (chihat_ref[:-1] + chihat_ref[1:])
    w_ref = gamh_ref * np.abs(gfac_h)
    rep.check(
        f"{name}: W = GAMH*|1+av_sig*g| (H-coupling)",
        np.allclose(d["w"][:-1], w_ref, rtol=1e-9, atol=1e-14),
        f"max diff {np.abs(d['w'][:-1]-w_ref).max():.2e}",
    )
    eps_ref = (
        gamh_ref * av_sig * d["nmj"][:-1] / (1.0 + w_ref)
    )
    rep.check(
        f"{name}: EPS = GAMH*av_sig*NMID_ON_J/(1+W)",
        np.allclose(d["eps"][:-1], eps_ref, rtol=1e-9, atol=1e-14),
        f"max diff {np.abs(d['eps'][:-1]-eps_ref).max():.2e}",
    )

    # HU/HL/HS midpoint H coefficients
    hu_ref = (
        konj[1:] * d["qf"][1:] / (1.0 + w_ref) / d["dtau"][:-1]
    )
    hl_ref = (
        konj[:-1] * d["qf"][:-1] / (1.0 + w_ref) / d["dtau"][:-1]
    )
    rep.check(
        f"{name}: HU/HL midpoint coefficients",
        np.allclose(d["hu"][:-1], hu_ref, rtol=1e-9, atol=1e-14)
        and np.allclose(d["hl"][:-1], hl_ref, rtol=1e-9, atol=1e-14),
        f"max diffs {np.abs(d['hu'][:-1]-hu_ref).max():.2e}/"
        f"{np.abs(d['hl'][:-1]-hl_ref).max():.2e}",
    )
    gamh_sgn = con_gamh * gfac_h
    hs_ref = np.where(
        gamh_sgn < 0, w_ref / (1.0 + w_ref), w_ref / (1.0 + w_ref)
    )
    rep.check(
        f"{name}: HS = W/(1+W) (pinned; both branches)",
        np.allclose(d["hs"][:-1], hs_ref, rtol=1e-9, atol=1e-14),
        f"max diff {np.abs(d['hs'][:-1]-hs_ref).max():.2e}",
    )

    # Interior matrix rows: TB = DTAUONQ*(1-COH) + PSI + HL + HU_prev
    # with the EPS depth-upwind splits folded in.
    coh = d["esec"] / d["chi"]
    tb_ref = (
        d["dtauonq"] * np.maximum(0.0, 1.0 - coh)
        + d["psi"] + d["hl"] + np.roll(d["hu"], 1)
    )
    ta_ref = -np.roll(d["hl"], 1)
    tc_ref = -d["hu"]
    eps_i = d["eps"]
    for i in range(1, n - 1):
        e_prev, e_cur = eps_i[i - 1], eps_i[i]
        if e_prev >= 0:
            ta_ref[i] -= 2.0 * e_prev
        else:
            tb_ref[i] -= 2.0 * e_prev
        if e_cur >= 0:
            tb_ref[i] += 2.0 * e_cur
        else:
            tc_ref[i] += 2.0 * e_cur
    rep.check(
        f"{name}: TA/TB/TC interior assembly (incl EPS depth-upwind)",
        np.allclose(d["ta"][i2], ta_ref[i2], rtol=1e-9, atol=1e-14)
        and np.allclose(d["tb"][i2], tb_ref[i2], rtol=1e-9, atol=1e-14)
        and np.allclose(d["tc"][i2], tc_ref[i2], rtol=1e-9, atol=1e-14),
        f"max diffs {np.abs(d['ta'][i2]-ta_ref[i2]).max():.2e}/"
        f"{np.abs(d['tb'][i2]-tb_ref[i2]).max():.2e}/"
        f"{np.abs(d['tc'][i2]-tc_ref[i2]).max():.2e}",
    )
    rep.check(
        f"{name}: interior diagonals positive (M-matrix)",
        bool((d["tb"][i2] > 0).all() and (d["ta"][i2] <= 0).all()
             and (d["tc"][i2] <= 0).all()),
    )
    return d


def check_admissible(rep, name, tag):
    d = load(tag)
    fin = np.isfinite(d["rsqj"]).all() and np.isfinite(d["rsqh"]).all()
    rep.check(f"{name}: RSQJ/RSQH finite", bool(fin))
    rep.check(
        f"{name}: J >= 0 (no significant negative moments)",
        bool((d["rsqj"] >= -1e-30).all()),
        f"min {d['rsqj'].min():.3e}",
    )
    jn = np.loadtxt(WORKDIR / f"{tag}_jnext.txt")
    rep.check(
        f"{name}: J_NEXT_FRQ store finite and non-negative",
        bool(np.isfinite(jn[:, 1:]).all() and (jn[:, 1:] >= 0).all()),
    )


def residual_check(rep, name, tag):
    """A.x = b on the dumped last-frequency solve (pre-repair)."""
    d = load(tag)
    n = len(d["r"])
    x, b = d["rsqj"], d["xmsave"]
    res = (
        d["ta"] * np.roll(x, 1) + d["tb"] * x + d["tc"] * np.roll(x, -1)
        - b
    )
    rel = np.abs(res[1 : n - 1]) / np.maximum(
        np.abs(b[1 : n - 1]), 1e-6 * np.abs(b).max()
    )
    rep.check(
        f"{name}: tridiagonal residual |Ax-b| small",
        rel.max() < 1e-6,
        f"max rel {rel.max():.2e}",
    )


def kv(out, key):
    for line in out.splitlines():
        if line.strip().startswith(key + "="):
            return float(line.split("=")[1].strip().split()[0])
    return None


def main():
    rep = Report("MOM_J_CMF_NONMONO_V1 discrete operator")
    exe = build()

    # --- Coefficient-level verification on two field shapes --------
    d1 = run(exe, "comp1", "nonj", "4", "1", "diff", "none", "0",
             "t_comp1", cwd=WORKDIR)
    rep.check("comp1: driver completed", "RESULT_OK" in d1.stdout)
    check_coefficients(rep, "comp1", "t_comp1")
    residual_check(rep, "comp1", "t_comp1")
    check_admissible(rep, "comp1", "t_comp1")

    d2 = run(exe, "homo", "nonj", "4", "1", "diff", "none", "0",
             "t_homo", cwd=WORKDIR)
    rep.check("homo: driver completed", "RESULT_OK" in d2.stdout)
    check_coefficients(rep, "homo", "t_homo")
    residual_check(rep, "homo", "t_homo")
    check_admissible(rep, "homo", "t_homo")

    # --- Upwind direction: compressing cells read the lagged map ---
    # J_NEXT_FRQ column k+1 feeds the compressing cells of solve k.
    # With MOM_GAM_PIN the explicit weight equals PSI; check the store
    # is actually being read (nonzero lagged contribution) where Q<0.
    jn = np.loadtxt(WORKDIR / "t_comp1_jnext.txt")
    dd = load("t_comp1")
    ncmp = int((dd["q"] < 0).sum())
    rep.check(
        "comp1: compressing rows exist on the expanded grid",
        ncmp > 0,
        f"{ncmp} rows",
    )
    rep.check(
        "comp1: lagged store populated at compressing slots",
        bool((jn[:, 1:] > 0).all()),
    )

    # --- Fixed-point drift: max change must contract across sweeps -
    it_lines = [l for l in d1.stdout.splitlines() if l.startswith("IT=")]
    chgs = [float(l.split("MOM_GS_MAXCHG=")[1].split()[0])
            for l in it_lines]
    rep.check(
        "comp1: lagged-map drift contracts over sweeps",
        len(chgs) >= 4 and chgs[-1] < chgs[1],
        f"trajectory {[f'{c:.2e}' for c in chgs]}",
    )

    # --- Zero-GAM resonance: sigma_true = -1/K band ------------------
    dz = run(exe, "zero", "nonj", "6", "1", "diff", "none", "0",
             "t_zero", cwd=WORKDIR)
    rep.check("zero: driver completed", "RESULT_OK" in dz.stdout)
    d = load("t_zero")
    t1 = 1.0 + d["sigt"] * d["konj"]
    exact = t1 == 0.0
    near = (np.abs(t1) > 0.0) & (np.abs(t1) < 1e-8)
    rep.check(
        "zero: exact resonant cells exist (1+sigma*K=0)",
        int(exact.sum()) > 0,
        f"{int(exact.sum())} exact, {int(near.sum())} near cells",
    )
    # Exact-zero exception: GAM_FAC=0 -> Q=0 -> chi_hat=chi at the
    # resonant cells; just off exact zero |fac|=1e-2 floor.
    rep.check(
        "zero: Q=0 at exact resonance (exact-zero exception)",
        bool((d["q"][exact] == 0.0).all())
        and bool((d["chihat"][exact] == d["chi"][exact]).all()),
    )
    if near.any():
        # The floor biases |GAM_FAC| to 1e-2 on near-resonant cells:
        # |Q| = |GAM|/chi_hat/DLN with |GAM|=CON_GAM*1e-2 -- a 1%-of-
        # floor systematic.  Bound it rather than demand zero.
        rep.check(
            "zero: near-resonant |Q| ~ floor*|GAM| level (bounded)",
            bool((np.abs(d["q"][near]) < 0.02).all()),
            f"max|Q| {np.abs(d['q'][near]).max():.2e}",
        )
    # Physical admissibility at and around the band.
    check_admissible(rep, "zero", "t_zero")
    rep.check(
        "zero: no solver failures at resonance",
        kv(dz.stdout, "N_DGTSV") == 0 and kv(dz.stdout, "N_NONF") == 0,
    )

    # Resonance continuity: cell-to-cell J jumps across the band edge
    # stay at the smooth-slope level -- no step discontinuity from the
    # floored/exact-zero factor switch.
    band = np.where(np.abs(t1) < 1e-8)[0]
    if len(band) > 0:
        lo, hi = band[0], band[-1]
        seg = d["rsqj"][max(0, lo - 2) : min(len(d["r"]), hi + 3)]
        step = np.abs(np.diff(seg)) / np.maximum(seg[:-1], 1e-30)
        rep.check(
            "zero: no step discontinuity in J across resonance",
            step.max() < 0.05,
            f"max cell-to-cell step {step.max():.2%}",
        )

    # --- Sign-flip gauntlet and resolved shock ----------------------
    for mode in ("sigx", "jump", "strong", "dec"):
        r = run(exe, mode, "nonj", "6", "1", "diff", "none", "0",
                f"t_{mode}", cwd=WORKDIR)
        rep.check(f"{mode}: driver completed",
                  "RESULT_OK" in r.stdout)
        rep.check(
            f"{mode}: no DGTSV/non-finite failures",
            kv(r.stdout, "N_DGTSV") == 0
            and kv(r.stdout, "N_NONF") == 0,
        )
        rep.check(
            f"{mode}: no significant negative J",
            kv(r.stdout, "NEGJ_LRG") == 0,
            f"NEGJ_LRG={kv(r.stdout, 'NEGJ_LRG')}",
        )
        check_admissible(rep, mode, f"t_{mode}")

    # --- Boundary treatments: OUT_BC_TYPE=2 row + inner methods ----
    for inbm in ("diff", "zfx"):
        r = run(exe, "comp1", "nonj", "4", "2", inbm, "none", "0",
                f"t_bc2_{inbm}", cwd=WORKDIR)
        rep.check(f"outbc=2/{inbm}: completed",
                  "RESULT_OK" in r.stdout,
                  r.stdout[-200:] if "RESULT_OK" not in r.stdout else "")
        check_admissible(rep, f"outbc=2/{inbm}", f"t_bc2_{inbm}")

    # --- G_ONLY path ------------------------------------------------
    rg = run(exe, "comp1", "gonly", "4", "1", "diff", "none", "0",
             "t_gonly", cwd=WORKDIR)
    rep.check("gonly: completed", "RESULT_OK" in rg.stdout)
    check_admissible(rep, "gonly", "t_gonly")

    # --- Checkpoint round-trip (MOMJNXT_GS) -------------------------
    rc = run(exe, "ckpt", "nonj", "1", "1", "diff", "none", "0",
             "t_ckpt", cwd=WORKDIR)
    rep.check(
        "ckpt: round-trip restores J/H_NEXT bit-exact",
        "CKPT_OK= T" in rc.stdout and "CKPT_MAXDIFF=  0.00000000E+00"
        in rc.stdout,
    )
    rep.check(
        "ckpt: grid fingerprint mismatch rejected, store zeroed",
        "CKPT_FP_OK= F" in rc.stdout
        and "CKPT_FP_STORE=  0.00000000E+00" in rc.stdout,
    )
    rep.check(
        "ckpt: stale tag rejected",
        "CKPT_TAG_OK= F" in rc.stdout,
    )

    rep.summary()


if __name__ == "__main__":
    main()
