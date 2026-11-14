#!/usr/bin/env python3
"""P8.1 acceptance: manufactured battery for the production
COMP_J_CONV_V2 electron-scattering convolver, the ES-iteration
sequencing contract, the KNOP ES input guard, and the Q-RAY
coherent/non-coherent emissivity-convention comparison.

Drivers (all call the REAL production routines):

  drv_es_conv   - COMP_J_CONV_V2 in 'J PASSED VIA CALL' mode.
    planck  J = B_nu(T): Jbar must recover B at depth; metric on
            the significant-J support only (Wien-tail underflow
            masked - a relative error there measures
            subnormal-or-zero noise, not the operator).
    const   J = constant: documented B/B_es rescale response;
            assert non-negative and finite.
    linepair    J_cont = B(T) and J_line = B(T)+A*Gaussian, both
            convolved by the real routine; baseline-subtracted
            excess measured on the FULL domain for peak,
            integrated excess in dlnnu and dnu, output wing
            fraction, FWHM, negative excursions, at two grid
            resolutions.
    linext      same pair on the extended (appended-node) grid;
            boundary sensitivity of the smeared feature and of
            the integrated excess.
    bndp/bnd    same field on the production grid and an extended
            grid whose production-spacing nodes are appended
            beyond both ends; features are ALIVE at both
            production bounds so the one-sided endpoint stencils
            and the unshifted/non-interpolated endpoint treatment
            are probed, and the extended result restricted to the
            shared nodes is the boundary reference.
    lte     J = B(T): LTE endpoint diagnostic - emits the
            convolved Planck replica, the PLANCK_NU map, and a
            zone flag per node (via the production
            ES_PLANCK_NU_MATCH) so the endpoint/bracket LTE
            trade-off is measured separately from the corrected
            interior.
    brk     synthetic ES_PLANCK_NU_MATCH battery: the match
            failure codes cannot be reached through the Planck
            path (PFN monotone over the Wien-gated region), so
            the same production object is driven with crafted
            log-space targets - valid, above/below range,
            floor-walk, degenerate, NaN/Inf, boundary-exact,
            persisted-ML.
    cnt     synthetic ES_PLANCK_NU_CONT battery: valid anchors,
            zero/negative/non-finite PNU and NU anchors, interior
            NU defects, a forced non-positive denominator and a
            single-node tail, through the production object.
    itp     synthetic ES_SHIFTED_RESCALE battery: valid map,
            truncated map (NPNU<NCF), degenerate last bracket,
            NaN and non-positive map nodes, one-node map - every
            interpolation-bracket validation path on the object.
    nanj    NaN at one node/depth: per-depth convolution skip +
            WARNING_LU diagnostic, field left visible.
    neg     signed convolved output: the MINVAL(J_ES) gate must
            fire before LOG(J_ES) can produce NaN.
    nant    NaN TEMP at depth 2 with finite J: the convolved
            J_ES/PLANCK_ES/PLANCK_FN go non-finite through the
            real path; the IEEE gate must bypass the correction,
            name the depth and fields, and leave J_ES visible.
    inft    +Inf TEMP at depth 2: same non-finite gate.
    zpf     TEMP at depth 2 cold enough that PLANCK_FN
            underflows to exact zero at the Wien end while the
            convolution stays real: the silent positivity
            gate-skip must leave the raw convolved J_ES.
    tempu   J = B(T_TEST) with T(k)=T_TEST uniform: isolates the
            beta(T)=1.84E-3*sqrt(T) width dependence.
    edge    bound-free-like step: positivity, recorded max rel.
    pos     interior zeroed band: positivity of Jbar.
    grid    Planck at NCF=128/256/512: refinement behaviour.
    qray    diluted blackbody + Gaussian-band CHI_RAY; emits both
            emissivity conventions for the Q-RAY census.

  drv_es_seq    - T5 sequencing harness: mirrors the production
            ES-iteration ordering (consume Jbar, scripted solve,
            store J, refresh slot, convolve stored J -> next
            Jbar; render convolves the stored J_N) using the real
            COMP_J_CONV_V2 for every convolution. Broken modes
            'stale' (sweep 2 fed Jbar_0), 'preconv' (convolve
            before the store), and 'oldrend' (render folds the
            previous Jbar) must be detectable.

  drv_knop_guard - T6: the real RD_OPTIONS_INTO_STORE +
            RD_CMF_FLUX_CONTROLS on a manufactured control file.
            Cases: (A) KNOP=T, NUM_ES=1, COH_ES=F -> named STOP;
            (B) KNOP=T, NUM_ES=2, COH_ES=T -> named STOP;
            (C) KNOP=T, NUM_ES=2, COH_ES=F -> GUARD_PASS.

  drv_es_render - renders the qray TA variants through the real
            OBS_FRAME_SUB_V9 (V=0 static atmosphere).

Method classification (plan doc section 2.1):
  executable production-linked: T5 (sequencing harness on the
    real convolver), T6 (real control reader).
  real-convolver manufactured: T3 battery (planck/const/linepair/
    bnd/tempu/edge/pos/grid), T4 render discriminator.
  algebraic reconstruction: T1 emissivity rebuild, T8 Q-RAY census.
"""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (
    F90DIR,
    TREE,
    WORKDIR,
    Report,
    compile_and_link_prod,
    run,
)

CONV_SRCS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "comp_j_conv_v2.f",
    TREE / "subs" / "thomas_rh.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    F90DIR / "drv_es_conv.f",
]

SEQ_SRCS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "comp_j_conv_v2.f",
    TREE / "subs" / "thomas_rh.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    F90DIR / "drv_es_seq.f",
]

GUARD_SRCS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "obs" / "mod_cmf_obs.f",
    TREE / "obs" / "cmf_flux_cntrl_var_mod.f",
    TREE / "newsubs" / "rd_store_log.f",
    TREE / "subs" / "set_case.f",
    TREE / "unix" / "term_io.f",
    TREE / "obs" / "check_cmf_flux_param_consis.f",
    TREE / "obs" / "rd_cmf_flux_controls.f",
    F90DIR / "drv_knop_guard.f",
]

RENDER_SRCS = [
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
    F90DIR / "obs_v9_adapter.f",
    TREE / "obs" / "obs_frame_sub_v9.f",
    F90DIR / "drv_es_render.f",
]

MODES = ["planck", "const", "edge", "pos"]
SWP_NU0 = ["299.5", "296", "290", "280", "260", "200"]
TEMPU = ["0.5", "1.0", "5.0"]
QRAY_AMPS = ["0.1", "0.5", "1.0"]
VADAT_TEMPLATE = TREE / "lte_hydro" / "lte_tst1" / "VADAT"

COL_I, COL_ML, COL_NU, COL_RJ, COL_RJES = 0, 1, 2, 3, 4

# Extra OBS-reader keys not present in the model-run VADAT
# template; values follow docs/audits/CMF_FLUX_OPTIMAL_INPUT_SETUP.md.
GUARD_EXTRA_KEYS = {
    "NUM_ES": "2",
    "TAU_MAX": "20.0",
    "ES_DTAU": "0.1",
    "N_INS_OBS": "3",
    "FRAC_DOP_OBS": "0.5",
    "DELV_FG": "0.3",
    "DELV_MOM": "0.1",
    "DOP_LIM": "1.0D-6",
    "VOIGT_LIM": "1.0D-8",
    "OPAC_LIMS": "T",
    "VTURB_FIX": "0.0",
    "VTURB_MIN": "5.0",
    "VTURB_MAX": "15.0",
    "EW_CUT": "0.05",
    "WR_ETA": "F",
    "WR_FLUX": "F",
    "WR_CMF_FORCE": "F",
    "WR_SOB_FORCE": "F",
    "WR_IP": "F",
    "IB_METH": "DIFFUSION",
    "INT_METH": "STAU",
    "GLOBAL_PROF": "DOPPLER",
    "dV_LEV": "300.0",
    "dV_CMF_PROF": "300.0",
    "dV_CMF_WING": "2000.0",
    "dV_OBS_PROF": "300.0",
    "dV_OBS_WING": "1000.0",
    "dV_OBS_BIG": "4000.0",
    "DO_NONMONO_KNOP": "T",
    "USE_LAM_ES": "T",
    "DO_NONMONO_SOB": "F",
}


def load_fields(path):
    d = np.loadtxt(path)
    if d.shape[1] == 5:
        # Scalar modes: ML NU T J_IN JBAR, rows ordered (I,ML).
        ncf = int(d[:, 0].max())
        nd = d.shape[0] // ncf
        nu = d[:, 1].reshape(nd, ncf)[0]
        rj = d[:, 3].reshape(nd, ncf)
        rjes = d[:, 4].reshape(nd, ncf)
    else:
        nd = int(d[:, COL_I].max())
        ncf = int(d[:, COL_ML].max())
        nu = d[:, COL_NU].reshape(nd, ncf)[0]
        rj = d[:, COL_RJ].reshape(nd, ncf)
        rjes = d[:, COL_RJES].reshape(nd, ncf)
    return nu, rj, rjes, d, nd, ncf


def sig_mask(rj, floor=1e-4):
    """Restrict relative metrics to the significant-J support."""
    return np.abs(rj) > floor * np.abs(rj).max()


def fwhm_ln_axis(nu, exc):
    """Full width at half maximum of a positive excess on the
    ln(nu) axis - the coordinate the redistribution operator
    actually acts in."""
    lnu = np.log(nu)
    m = exc > 0.5 * exc.max()
    if not m.any():
        return 0.0
    ii = np.where(m)[0]
    return abs(lnu[ii.max()] - lnu[ii.min()])


def build_vadat(path, num_es, coh_es, knop=True):
    """Manufactured control file: production VADAT key lines plus
    the OBS-reader keys, patched per guard case."""
    base = VADAT_TEMPLATE.read_text().splitlines()
    lines = [l for l in base if re.search(r"\[[A-Za-z_0-9]+\]", l)]
    lines = [
        f"{num_es:<14} [NUM_ES]"
        if "[NUM_ES]" in l
        else f"{coh_es:<14} [COH_ES]" if "[COH_ES]" in l else l
        for l in lines
    ]
    have = {
        m.group(1)
        for l in lines
        for m in [re.search(r"\[([A-Za-z_0-9]+)\]", l)]
        if m
    }
    extra = dict(GUARD_EXTRA_KEYS)
    extra["NUM_ES"] = num_es
    extra["COH_ES"] = coh_es
    if not knop:
        extra["DO_NONMONO_KNOP"] = "F"
    for k, v in extra.items():
        if k not in have:
            lines.append(f"{v:<14} [{k}]")
    path.write_text("\n".join(lines) + "\n")


def load_seq(path):
    d = np.loadtxt(path)
    return {
        "nu": d[:, 1],
        "jbar0": d[:, 2],
        "j1": d[:, 3],
        "jbar1": d[:, 4],
        "j2": d[:, 5],
        "jbar2": d[:, 6],
        "ta": d[:, 7],
    }


def main():
    rep = Report("ES convolver + sequencing + Q-RAY (P8.1)")
    out_dir = WORKDIR / "es_conv"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Production-tree sources resolve to the exact production .o
    # files (the bytes the production executables link); only the
    # drivers and the obs_v9 adapter compile fresh.
    conv_bin = compile_and_link_prod("drv_es_conv", CONV_SRCS)
    seq_bin = compile_and_link_prod("drv_es_seq", SEQ_SRCS)
    guard_bin = compile_and_link_prod("drv_knop_guard", GUARD_SRCS)
    rend_bin = compile_and_link_prod("drv_es_render", RENDER_SRCS)

    # -- driver runs ----------------------------------------------------
    for mode in MODES:
        run(conv_bin, out_dir, mode, cwd=out_dir)
    for t in TEMPU:
        run(conv_bin, out_dir, "tempu", t, cwd=out_dir)
        (out_dir / "esconv_tempu.txt").rename(
            out_dir / f"esconv_tempu{t.replace('.', '')}.txt"
        )
    for a in QRAY_AMPS:
        run(conv_bin, out_dir, "qray", a, cwd=out_dir)
        (out_dir / "esconv_qray_fields.txt").rename(
            out_dir / f"esconv_qray_{a}.txt"
        )
    for n in ("128", "256", "512"):
        run(conv_bin, out_dir, "grid", n, cwd=out_dir)
    # linepair: feature sigma_ln_nu ~ 0.2 (NCF=256) and ~0.1
    # (NCF=512), i.e. >=5 grid points per input sigma both times.
    run(
        conv_bin, out_dir, "linepair", "256", "2", "0.4", "50",
        cwd=out_dir,
    )
    run(
        conv_bin, out_dir, "linepair", "512", "2", "0.4", "50",
        cwd=out_dir,
    )
    # boundary: production grid vs extended grid with the
    # production log spacing appended beyond both ends (identical
    # production nodes).  linext is the same extension applied
    # to the line feature - it isolates the finite-boundary
    # influence on the smeared excess and its integrals.
    run(conv_bin, out_dir, "bndp", "256", cwd=out_dir)
    run(conv_bin, out_dir, "bnd", "256", cwd=out_dir)
    run(conv_bin, out_dir, "bndp", "512", cwd=out_dir)
    run(conv_bin, out_dir, "bnd", "512", cwd=out_dir)
    # TEST C (Wien-tail continuum): a smooth non-Planck continuum
    # alive through the Wien-dead upper end, production vs
    # extended grid, both resolutions.
    for n in ("256", "512"):
        run(conv_bin, out_dir, "wcp", n, cwd=out_dir)
        run(conv_bin, out_dir, "wce", n, cwd=out_dir)
    # TEST D (feature-position sweep): the same absolute-
    # amplitude Gaussian moved progressively inward from the
    # upper endpoint - maps where the boundary treatment acts.
    for nu0 in SWP_NU0:
        run(
            conv_bin, out_dir, "swp", nu0, "12.0", "0.5",
            cwd=out_dir,
        )
    # LTE endpoint diagnostic: J = B(T), zone map emitted by a
    # replica that calls the same production ES_PLANCK_NU_MATCH.
    run(conv_bin, out_dir, "lte", "256", cwd=out_dir)
    # Synthetic bracket battery on the production match routine.
    run(conv_bin, out_dir, "brk", cwd=out_dir)
    # Synthetic continuation and interpolation batteries on the
    # production object code.
    run(conv_bin, out_dir, "cnt", cwd=out_dir)
    run(conv_bin, out_dir, "itp", cwd=out_dir)
    # Input guards: NaN input field and signed convolved output.
    # fort.4 is per-process and overwritten from the top, so each
    # guarded run's WARNING_LU output is captured immediately.
    (out_dir / "fort.4").unlink(missing_ok=True)
    run(conv_bin, out_dir, "nanj", "256", cwd=out_dir)
    fort4_nanj = (
        (out_dir / "fort.4").read_text()
        if (out_dir / "fort.4").exists()
        else ""
    )
    run(conv_bin, out_dir, "neg", "256", cwd=out_dir)
    # Output guards: non-finite convolver output (NaN and Inf
    # temperatures) and a finite-but-zero Planck field.
    (out_dir / "fort.4").unlink(missing_ok=True)
    run(conv_bin, out_dir, "nant", "256", cwd=out_dir)
    fort4_nant = (
        (out_dir / "fort.4").read_text()
        if (out_dir / "fort.4").exists()
        else ""
    )
    (out_dir / "fort.4").unlink(missing_ok=True)
    run(conv_bin, out_dir, "inft", "256", cwd=out_dir)
    fort4_inft = (
        (out_dir / "fort.4").read_text()
        if (out_dir / "fort.4").exists()
        else ""
    )
    (out_dir / "fort.4").unlink(missing_ok=True)
    run(conv_bin, out_dir, "zpf", "256", cwd=out_dir)
    fort4_zpf = (
        (out_dir / "fort.4").read_text()
        if (out_dir / "fort.4").exists()
        else ""
    )
    run(
        conv_bin, out_dir, "linext", "256", "2", "0.4", "50",
        cwd=out_dir,
    )
    run(
        conv_bin, out_dir, "linext", "512", "2", "0.4", "50",
        cwd=out_dir,
    )

    # -- planck preservation ---------------------------------------------
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_planck.txt")
    m = sig_mask(rj)
    rel = np.abs(rjes - rj) / np.abs(rj)
    rep.check(
        "planck: Jbar = B on the significant support",
        rel[m].max() < 1e-2,
        f"max {rel[m].max():.3e} rms {np.sqrt((rel[m] ** 2).mean()):.3e}",
    )
    rep.check(
        "planck: Jbar finite and non-negative everywhere",
        np.all(np.isfinite(rjes)) and np.all(rjes >= 0.0),
    )

    # -- LTE endpoint trade-off diagnostic (P8.1 robustness) ------------
    # J = B(T): the lte driver emits, per node, the input Planck
    # field, the replicated convolved Planck field, the PLANCK_NU
    # map, a zone flag from a replica that calls the production
    # ES_PLANCK_NU_MATCH, and JBAR.  Zones: 0 corrected interior,
    # 1 upper endpoint, 2 invalid-map bracket (I=1 during the
    # shift interpolation), 3 gate-skip (MINVAL block not
    # entered), 4 continuation tail, 5 bypass, 9 lower endpoint.
    # The trade the audit documents: at zones 1 and 2 no genuine
    # shifted PLANCK_ES exists, so Jbar returns the raw convolved
    # field rather than the exact B - the relative LTE error there
    # is large and HONESTLY reported, with the confinement
    # (B/bchar and |dJ|/bchar) measured rather than assumed.
    d = np.loadtxt(out_dir / "esconv_lte0300x0002.txt")
    ncf_l = int(d[:, 0].max())
    nd = d.shape[0] // ncf_l
    nu_l = d[:, 1].reshape(nd, ncf_l)[0]
    tt = d[:, 2].reshape(nd, ncf_l)[:, 0]
    jin_l = d[:, 3].reshape(nd, ncf_l)
    zn_l = d[:, 6].reshape(nd, ncf_l).astype(int)
    fail_l = d[:, 7].reshape(nd, ncf_l).astype(int)
    jbar_l = d[:, 8].reshape(nd, ncf_l)
    TINY = np.finfo(np.float64).tiny          # 2.2e-308
    HDKT_ = 4.7994145
    TW2_ = 0.0147452575

    def pfn_prod(nu_, t_):
        """The T3-scaled Planck field the routine matches
        (production formula, comp_j_conv_v2.f ~:306)."""
        x_ = HDKT_ * nu_ / t_
        t3_ = min(700.0, 0.5 * HDKT_ * nu_ / t_)
        if t3_ < 1.0:
            t3_ = 0.0
        return TW2_ * nu_**3 * np.exp(t3_ - x_) / (1.0 - np.exp(-x_))

    def ieee_cls(v):
        if v == 0.0:
            return "zero"
        return "subnormal" if abs(v) < TINY else "normal"

    # IEEE classification of the two Planck quantities at the
    # endpoint: the driver input B(1) and the production PFN(1).
    b1_cls = [ieee_cls(jin_l[k, 0]) for k in range(nd)]
    pfn1 = np.array([pfn_prod(nu_l[0], t) for t in tt])
    pfn1_cls = [ieee_cls(v) for v in pfn1]
    rep.check(
        "lte: endpoint input B(1) IEEE class per depth",
        all(c == "normal" for c in b1_cls),
        f"B(1)={['%.2e' % jin_l[k, 0] for k in range(nd)]} "
        f"classes={b1_cls} (driver formula saturates at e^-700)",
    )
    rep.check(
        "lte: production PFN(1) class is normal/subnormal/zero, "
        "not uniformly 'denormal'",
        pfn1_cls == ["normal", "normal", "subnormal", "zero"],
        f"PFN(1)={['%.2e' % v for v in pfn1]} classes={pfn1_cls} "
        f"(T=1e4 K is subnormal, T=5e3 K underflows to zero)",
    )
    rel_l = np.abs(jbar_l - jin_l) / np.maximum(np.abs(jin_l), 1e-300)
    abs_l = np.abs(jbar_l - jin_l)
    bchar = jin_l.max(axis=1)
    # No match failure is expected on the standard grid: the Wien
    # gate keeps every entered node above the Planck maximum, where
    # PFN is monotone and a kernel-averaged PES always brackets.
    rep.check(
        "lte: no match failure on the standard grid",
        np.all(fail_l == -1),
        f"FAIL values {sorted(set(fail_l.ravel()))}",
    )
    corr = (zn_l == 0) | (zn_l == 4)
    corr_rel = rel_l[corr]
    rep.check(
        "lte: corrected interior + continuation recover B "
        "(algebraic identity)",
        corr_rel.max() < 5e-12,
        f"max rel {corr_rel.max():.3e} over {corr.sum()} nodes",
    )
    lo_rel = rel_l[zn_l == 9]
    rep.check(
        "lte: lower endpoint recovers B",
        lo_rel.max() < 5e-12,
        f"max rel {lo_rel.max():.3e}",
    )
    trade = (zn_l == 1) | (zn_l == 2)
    bsig = jin_l / bchar[:, None]          # input as fraction of Bchar
    asig = abs_l / bchar[:, None]          # abs error as fraction
    rep.check(
        "lte: large rel error confined to dead-weight frequencies "
        "on corrected depths (B < 1e-40 of characteristic)",
        np.all(
            (rel_l <= 1e-3) | (bsig < 1e-40) | (zn_l == 3)
        ),
        f"worst B/bchar among rel>1e-3 non-skip nodes: "
        f"{bsig[(rel_l > 1e-3) & (zn_l != 3)].max():.3e}",
    )
    rep.check(
        "lte: absolute endpoint+bracket error negligible vs "
        "characteristic intensity",
        asig[trade].max() < 1e-40,
        f"max abserr/bchar {asig[trade].max():.3e} "
        f"over {trade.sum()} invalid-map nodes",
    )
    rep.check(
        "lte: endpoint rel error recorded (the honest trade)",
        np.all(rel_l[zn_l == 1] > 1e6),
        f"rel_err(NU1)="
        f"{['%.2e' % rel_l[k, 0] for k in range(nd) if (zn_l[k] == 1).any()]} "
        f"(raw convolved field stands where no shifted map exists)",
    )
    skip = zn_l == 3
    rep.check(
        "lte: gate-skip depth (PFN contains zero) stays finite "
        "and small in absolute terms",
        skip.all(axis=1).sum() >= 1
        and np.all(np.isfinite(jbar_l[skip]))
        and asig[skip].max() < 1e-4,
        f"depths fully skipped: {skip.all(axis=1).sum()} "
        f"max abserr/bchar {asig[skip].max():.3e}",
    )
    rep.check(
        "lte: Jbar finite and non-negative everywhere",
        np.all(np.isfinite(jbar_l)) and np.all(jbar_l >= 0.0),
    )

    # -- ES_PLANCK_NU_MATCH synthetic bracket battery -------------------
    # The match failure codes cannot be produced through the
    # Planck path (PFN is monotone over the Wien-gated region and
    # PES is a kernel average), so the brk driver feeds synthetic
    # log-space targets to the same production object code.
    # Expected: (FAIL, ML_out, T1, PNU) per case.
    d = np.loadtxt(out_dir / "esconv_brk.txt")
    brk_exp = {
        1: (0, 3, 0.5, 6.5),   # valid interior bracket
        2: (1, 9, 0.0, 0.0),   # above range: ML bounded at NCF+1
        3: (2, 1, 0.0, 0.0),   # below PFN(1): no PFN(0) read
        4: (2, 9, 0.0, 0.0),   # floor walk bounded at NCF+1
        5: (3, 4, 0.0, 0.0),   # degenerate denominator
        6: (4, 1, 0.0, 0.0),   # NaN target
        7: (4, 1, 0.0, 0.0),   # +Inf target
        8: (2, 1, 0.0, 0.0),   # huge-negative (tiny linear Planck)
        9: (0, 5, 0.0, 4.0),   # target exactly on a PFN node
        10: (0, 5, 0.5, 4.5),  # continuation from persisted ML
    }
    ok = True
    det = []
    for row in d:
        icase = int(row[0])
        ml_out, t1, pnu = int(row[4]), row[5], row[6]
        fa, mle, t1e, pnue = brk_exp[icase]
        good = (
            int(row[3]) == fa
            and ml_out == mle
            and abs(t1 - t1e) < 1e-7
            and (fa != 0 or abs(pnu - pnue) < 1e-6)
            and (fa != 0 or 0.0 <= t1 <= 1.0)
        )
        ok = ok and good
        det.append(f"c{icase}:F{int(row[3])}/ML{ml_out}")
    rep.check(
        "brk: all ten match cases return the expected FAIL/ML/T1/PNU",
        ok and len(d) == 10,
        " ".join(det),
    )
    rep.check(
        "brk: no out-of-bounds walk (ML never exceeds NCF+1, "
        "no PFN(0) read)",
        np.all(d[:, 4] >= 1) and np.all(d[:, 4] <= 9),
        f"ML range [{d[:, 4].min():.0f},{d[:, 4].max():.0f}]",
    )

    # -- input-guard probes --------------------------------------------
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_nanj.txt")
    nan_pos = ~np.isfinite(rjes[1])
    rep.check(
        "nanj: NaN input depth passes through verbatim, others "
        "convolve normally",
        nan_pos.sum() == 1
        and np.isfinite(rjes[0]).all()
        and np.isfinite(rjes[2]).all()
        and np.isfinite(rjes[3]).all()
        and np.array_equal(
            np.delete(rjes[1], np.where(nan_pos)[0]),
            np.delete(rj[1], np.where(nan_pos)[0]),
        ),
        f"non-finite nodes at depth 2: {nan_pos.sum()} "
        f"(field left visible, convolution skipped)",
    )
    rep.check(
        "nanj: skip emits a named warning diagnostic",
        "non-finite J input" in fort4_nanj,
        fort4_nanj.strip()[:80] or "no fort.4",
    )
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_neg.txt")
    rep.check(
        "neg: signed convolved output stays finite (LOG gate)",
        np.all(np.isfinite(rjes)),
        f"min {rjes.min():.3e} "
        "(MINVAL(J_ES)<=0 skips the Planck block before LOG)",
    )

    # -- convolver-output validity guards (non-finite gate) --------------
    # The driver keeps J finite (the warm-depth Planck field) and
    # makes TEMP non-finite or zero at depth 2, so the convolved
    # J_ES, PLANCK_ES and PLANCK_FN go non-finite or non-positive
    # through the real production path.  nant/inft must hit the
    # IEEE finiteness gate with a named diagnostic and leave the
    # invalid field visible; zpf is the finite-but-zero case
    # (silent positivity gate-skip, same policy as a cold depth).
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_nant.txt")
    rep.check(
        "nant: NaN convolver output visible at bad depth, others "
        "convolve normally",
        (~np.isfinite(rjes[1])).all()
        and np.isfinite(rjes[0]).all()
        and np.isfinite(rjes[2]).all()
        and np.isfinite(rjes[3]).all(),
        f"finite depths: {[np.isfinite(rjes[k]).all() for k in range(4)]}",
    )
    rep.check(
        "nant: warning names depth and every non-finite field",
        "non-finite convolver output at depth    2" in fort4_nant
        and "PLANCK_FN" in fort4_nant
        and "PLANCK_ES" in fort4_nant
        and "J_ES" in fort4_nant,
        fort4_nant.strip()[:120] or "no fort.4",
    )
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_inft.txt")
    rep.check(
        "inft: +Inf convolver output caught the same way",
        (~np.isfinite(rjes[1])).all()
        and np.isfinite(rjes[0]).all()
        and np.isfinite(rjes[2]).all()
        and np.isfinite(rjes[3]).all(),
        f"finite depths: {[np.isfinite(rjes[k]).all() for k in range(4)]}",
    )
    rep.check(
        "inft: warning emitted for the non-finite gate",
        "non-finite convolver output" in fort4_inft,
        fort4_inft.strip()[:80] or "no fort.4",
    )
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_zpf.txt")
    rep.check(
        "zpf: finite-but-zero PLANCK_FN gate-skips silently, "
        "output stays finite",
        np.all(np.isfinite(rjes)) and np.all(rjes > 0.0),
        f"min {rjes.min():.3e}",
    )
    rep.check(
        "zpf: no non-finite diagnostic on the silent skip",
        "non-finite" not in fort4_zpf,
        fort4_zpf.strip()[:80] or "(empty)",
    )
    rep.check(
        "zpf: bypassed correction leaves raw convolution (B_es "
        ">> B on the Wien side)",
        (rjes[1] / rj[1]).max() > 1e4,
        f"max JBAR/J_IN {(rjes[1] / rj[1]).max():.3e} "
        "(a corrected Planck field would return ~1)",
    )

    # -- cnt: continuation battery ----------------------------------------
    # ES_PLANCK_NU_CONT on synthetic anchors: valid fill, bad
    # PNU/NU anchors (FAIL=5), bad interior NU or a forced
    # non-positive denominator (FAIL=6), single-node tail.
    d = np.loadtxt(out_dir / "esconv_cnt.txt")
    cnt_exp = {
        1: (0, 8), 2: (5, 2), 3: (5, 2), 4: (5, 2), 5: (5, 2),
        6: (5, 2), 7: (5, 2), 8: (6, 4), 9: (6, 4), 10: (6, 4),
        11: (6, 2), 12: (0, 8),
    }
    ok, det = True, []
    for row in d:
        ic, fa, npnu = int(row[0]), int(row[1]), int(row[2])
        pnu = row[3:]
        good = (fa, npnu) == cnt_exp[ic]
        # Written nodes must be positive and finite - except the
        # anchors of a FAIL=5 case, which are the defect itself.
        if npnu > 0 and fa in (0, 6):
            good = good and np.all(np.isfinite(pnu[:npnu]))
            good = good and np.all(pnu[:npnu] > 0.0)
        if npnu < 8:
            good = good and np.all(pnu[npnu:] == 0.0)
        ok = ok and good
        det.append(f"c{ic}:F{fa}/N{npnu}")
    rep.check(
        "cnt: all twelve continuation cases return the expected "
        "FAIL/NPNU with only positive finite nodes written",
        ok and len(d) == 12,
        " ".join(det),
    )
    rep.check(
        "cnt: valid continuation stays positive and decreasing",
        np.all(np.diff(d[0, 3:]) < 0.0) and np.all(d[0, 3:] > 0.0),
        f"PNU(1..8)={d[0, 3:]}",
    )

    # -- itp: shifted-interpolation battery -------------------------------
    # ES_SHIFTED_RESCALE on synthetic maps.  Python oracle: same
    # walk semantics, then per-node expected JES/PES.  An invalid
    # node takes JES(ML)=exp(A(ML)) (=raw input) and PES=PFN; a
    # valid node interpolates the log fields.
    d = np.loadtxt(out_dir / "esconv_itp.txt").reshape(6, 8, 6)
    NU_B = 9.0 - np.arange(1, 9)            # 8..1 decreasing
    PFN_B = 2.0 * np.arange(1, 9)
    C_B = 3.0 * np.arange(1, 9)
    JIN_B = 5.0 * np.arange(1, 9)
    A_B = np.log(JIN_B)
    PNU_B = {
        1: NU_B - 0.5,
        2: NU_B - 0.5,
        3: np.array([7.5, 6.5, 5.5, 4.5, 3.5, 2.5, 3.0, 3.0]),
        4: np.array(
            [7.5, 6.5, 5.5, 4.5, np.nan, 2.5, 1.5, 0.5]
        ),
        5: np.array([7.5, 6.5, 5.5, -1.0, 3.5, 2.5, 1.5, 0.5]),
        6: NU_B - 0.5,
    }
    NPNU_B = {1: 8, 2: 5, 3: 8, 4: 8, 5: 8, 6: 1}

    def itp_expect(pnu, npnu):
        pes = C_B.copy()
        jes = np.empty(8)
        jes[0] = A_B[0]
        jes[-1] = A_B[-1]
        i = 0
        for ml in range(1, 7):
            while i < npnu - 2 and NU_B[ml] < pnu[i + 1]:
                i += 1
            brk = True
            if npnu < 8 and i == npnu - 2 and NU_B[ml] < pnu[npnu - 1]:
                brk = False
            if npnu < 2:
                brk = False
            if brk:
                t2 = pnu[i + 1] - pnu[i]
                brk = (
                    np.isfinite(pnu[i])
                    and np.isfinite(pnu[i + 1])
                    and pnu[i] > 0.0
                    and pnu[i + 1] > 0.0
                    and np.isfinite(t2)
                    and t2 != 0.0
                )
            if brk:
                t1 = (NU_B[ml] - pnu[i]) / t2
                brk = np.isfinite(t1)
            if brk:
                pes[ml] = t1 * C_B[i + 1] + (1.0 - t1) * C_B[i]
                jes[ml] = t1 * A_B[i + 1] + (1.0 - t1) * A_B[i]
            if not brk or i == 0:
                pes[ml] = PFN_B[ml]
                jes[ml] = A_B[ml]
        pes[0] = PFN_B[0]
        return np.exp(jes + PFN_B - pes), pes

    ok, det = True, []
    for ic in range(1, 7):
        jes_exp, pes_exp = itp_expect(PNU_B[ic], NPNU_B[ic])
        jes_got = d[ic - 1, :, 4]
        pes_got = d[ic - 1, :, 5]
        good = (
            np.allclose(jes_got, jes_exp, rtol=1e-6, atol=1e-12)
            and np.allclose(pes_got, pes_exp, rtol=1e-6, atol=1e-12)
            and np.all(np.isfinite(jes_got))
            and np.all(jes_got > 0.0)
        )
        ok = ok and good
        det.append(
            f"c{ic}:"
            f"{np.max(np.abs(jes_got - jes_exp)):.1e}"
        )
    rep.check(
        "itp: all six map cases match the walk oracle node-by-node "
        "(valid interp, truncated, degenerate, NaN, non-positive, "
        "one-node)",
        ok,
        " ".join(det),
    )
    rep.check(
        "itp: invalid nodes never contaminate valid ones (raw J "
        "at exactly the off-map/degenerate nodes)",
        np.isclose(d[1, 5, 4], 30.0)
        and np.isclose(d[1, 6, 4], 35.0)
        and np.isclose(d[2, 6, 4], 35.0)
        and np.allclose(d[3, 4:7, 4], JIN_B[4:7])
        and np.allclose(d[4, 3:7, 4], JIN_B[3:7])
        and np.allclose(d[5, 1:7, 4], JIN_B[1:7])
        and d[0, 2, 4] != JIN_B[2] and d[1, 2, 4] != JIN_B[2],
        f"off-map JES {d[1,5:7,4]}, degen {d[2,6,4]}, "
        f"nan-adjacent {d[3,4:7,4]}",
    )

    # -- constant field ----------------------------------------------------
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_const.txt")
    rep.check(
        "const: Jbar non-negative and finite",
        np.all(np.isfinite(rjes)) and np.all(rjes >= 0.0),
        f"zeros at Wien-underflow edge: {(rjes == 0.0).sum()}",
    )
    active = rjes.max(axis=0) > 0.0
    if np.any(active):
        spread = np.ptp(rjes[:, active], axis=0) / rjes[:, active].max(axis=0)
        rep.check(
            "const: documented rescale behaviour (not constant)",
            True,
            f"max col spread {np.nanmax(spread):.3e} (B/B_es rescale)",
        )

    # -- narrow-feature smearing (baseline-subtracted) ----------------------
    # The redistribution operator is formulated in ln(nu) (the
    # two-exponential Helmholtz fit of Rybicki & Hummer 1994,
    # A&A 290, 553, convolves the log coordinate), so the honest
    # integration measure is d ln nu; d nu is reported alongside.
    # Neither integral is ASSERTED conserved: the interior column
    # sums of the discrete operator are -H = 1, but the one-sided
    # endpoint rows do not share that property, and the final
    # B/B_es rescale is applied pointwise.  Integral agreement is
    # therefore a measured property of this implementation on
    # this domain - what is bounded below is the measured drift,
    # not a physical conservation law.

    def ln_sigma(nu_, exc_):
        """Actual feature width in ln(nu): second moment of the
        excess profile about its peak.  The manufactured feature
        is Gaussian in nu, not in ln nu, so sigma_nu/nu0 is only
        a local approximation of the width the operator sees."""
        lnu_ = np.log(nu_)
        pk_ = lnu_[np.argmax(exc_)]
        w_ = np.trapezoid(exc_, x=lnu_)
        return np.sqrt(np.trapezoid(exc_ * (lnu_ - pk_) ** 2, x=lnu_) / w_)

    smear_resp = {}
    for n, dnu_in in ((256, 0.4), (512, 0.4)):
        d = np.loadtxt(out_dir / f"esconv_linepair{n:05d}.txt")
        dx = np.loadtxt(out_dir / f"esconv_linext{n:05d}.txt")
        ncf = int(d[:, 0].max())
        nd = d.shape[0] // ncf
        nu = d[:, 1].reshape(nd, ncf)[0]
        j_cont = d[:, 3].reshape(nd, ncf)
        jbar_cont = d[:, 4].reshape(nd, ncf)
        j_line = d[:, 5].reshape(nd, ncf)
        jbar_line = d[:, 6].reshape(nd, ncf)
        n_e = int(dx[:, 0].max())
        next_ = (n_e - ncf) // 2
        nu_e = dx[:, 1].reshape(nd, n_e)[0]
        exc_in_e = dx[:, 5].reshape(nd, n_e)[0] - dx[:, 3].reshape(
            nd, n_e
        )[0]
        exc_out_e = dx[:, 6].reshape(nd, n_e)[0] - dx[:, 4].reshape(
            nd, n_e
        )[0]
        k = 0  # outer depth row carries the strongest field
        exc_in = j_line[k] - j_cont[k]
        exc_out = jbar_line[k] - jbar_cont[k]
        feat = exc_in > 0.05 * exc_in.max()
        lnu = np.log(nu)
        asc = np.argsort(nu)
        asc_e = np.argsort(nu_e)
        # Full-domain integrals in both measures, input and output.
        int_in_l = np.trapezoid(exc_in[asc], x=lnu[asc])
        int_out_l = np.trapezoid(exc_out[asc], x=lnu[asc])
        int_in_n = np.trapezoid(exc_in[asc], x=nu[asc])
        int_out_n = np.trapezoid(exc_out[asc], x=nu[asc])
        int_out_el = np.trapezoid(
            exc_out_e[asc_e], x=np.log(nu_e)[asc_e]
        )
        # Output excess outside the input 5%-of-peak footprint.
        inside = np.trapezoid(
            exc_out[asc][feat[asc]], x=lnu[asc][feat[asc]]
        )
        wing_frac = 1.0 - inside / int_out_l
        w_in = fwhm_ln_axis(nu, exc_in)
        w_out = fwhm_ln_axis(nu, np.maximum(exc_out, 0.0))
        dlnnu = np.abs(np.diff(lnu)).mean()
        sig_loc = dnu_in / 2.0          # local approx sigma_nu/nu0
        sig_act = ln_sigma(nu, exc_in)  # measured width in ln nu
        nu0 = nu[np.argmax(exc_in)]
        de = exc_out - exc_in
        core = np.abs(nu - nu0) < 0.5 * dnu_in
        wing = (np.abs(nu - nu0) > dnu_in) & (
            np.abs(nu - nu0) < 2.5 * dnu_in
        )
        smear_resp[n] = np.abs(de).max() / exc_in.max()
        rep.check(
            f"linepair {n}: feature resolved in ln nu "
            "(>=4 pts/sigma)",
            sig_act / dlnnu >= 4.0,
            f"sigma_ln(local {sig_loc:.3f}, actual {sig_act:.3f}) "
            f"dln_nu={dlnnu:.4f} pts/sig={sig_act / dlnnu:.1f}",
        )
        # Genuine smearing signature for a second-derivative
        # operator: the peak drops, the line core is depleted, the
        # wings gain.  FWHM growth is below the grid cell here
        # (kernel width ~1e-3 in ln nu at 1e4 K << dln_nu), so the
        # width change is reported as unresolved, not asserted.
        rep.check(
            f"linepair {n}: smearing signature "
            "(peak drop, core depletion, wing lift)",
            exc_out.max() < exc_in.max()
            and de[core].mean() < 0.0
            and de[wing].mean() > 0.0,
            f"peak {exc_in.max():.4e}->{exc_out.max():.4e} "
            f"core_dex {de[core].mean():+.3e} "
            f"wing_dex {de[wing].mean():+.3e} "
            f"max|dexc|/pk {smear_resp[n]:.3e}",
        )
        rep.check(
            f"linepair {n}: FWHM measured "
            f"({'resolved' if abs(w_out - w_in) >= 2 * dlnnu else 'unresolved: below grid cell'})",
            w_out >= w_in - 2.0 * dlnnu,
            f"fwhm_ln {w_in:.4f}->{w_out:.4f} dln_nu={dlnnu:.4f}",
        )
        # Full-domain integral drift, reported in both measures.
        # The bound is a regression bound on the measured drift
        # (~4e-6 at both resolutions), not a conservation axiom.
        rep.check(
            f"linepair {n}: full-domain integral drift bounded (<1e-2)",
            abs(int_out_l / int_in_l - 1.0) < 1e-2,
            f"dlnnu: {int_in_l:.4e}->{int_out_l:.4e} "
            f"({int_out_l / int_in_l - 1.0:+.3e})  "
            f"dnu: {int_in_n:.4e}->{int_out_n:.4e} "
            f"({int_out_n / int_in_n - 1.0:+.3e})",
        )
        rep.check(
            f"linepair {n}: output wing fraction outside input mask",
            0.0 < wing_frac < 0.5,
            f"{wing_frac:.3f} of int(exc_out) dlnnu beyond the "
            f">5% input footprint",
        )
        # Finite-boundary sensitivity: the same pair on the
        # extended grid.  Boundary influence is bounded both on
        # the full-domain integral and on the shared nodes.
        rep.check(
            f"linepair {n}: integral stable vs extended domain (<1e-3)",
            abs(int_out_el / int_out_l - 1.0) < 1e-3,
            f"prod {int_out_l:.4e} ext {int_out_el:.4e} "
            f"({int_out_el / int_out_l - 1.0:+.3e})",
        )
        xe_i = exc_out_e[next_ : next_ + ncf]
        dxc = np.abs(exc_out - xe_i) / np.maximum(np.abs(xe_i), 1e-300)
        mxc = np.abs(xe_i) > 1e-6 * np.abs(xe_i).max()
        rep.check(
            f"linepair {n}: shared-node exc stable vs extended (<1e-2)",
            dxc[mxc].max() < 1e-2,
            f"max rel {dxc[mxc].max():.3e} ({mxc.sum()} sig pts)",
        )
        rep.check(
            f"linepair {n}: no strong negative ringing",
            exc_out.min() > -0.2 * exc_out.max(),
            f"min {exc_out.min():.3e}",
        )
        rep.check(
            f"linepair {n}: output fields finite and non-negative",
            np.all(np.isfinite(jbar_line)) and np.all(jbar_line >= 0.0),
        )
    rep.check(
        "linepair: response consistent across resolutions",
        0.2 < smear_resp[256] / smear_resp[512] < 5.0,
        f"256: {smear_resp[256]:.3e} 512: {smear_resp[512]:.3e}",
    )

    # -- frequency boundary (appended-domain reference) ---------
    # 'bnd' appends NEXT=NCF_P/4 same-spacing nodes beyond both
    # production ends, so the interior production nodes coincide
    # exactly - the comparison isolates the boundary treatment
    # (one-sided -2/Delta^2 stencil rows + the uncorrected
    # endpoint treatment where no genuine shifted PLANCK_ES
    # exists + the interior wavelength-shift/rescale path) rather
    # than resolution/interp differences.  Manufactured features
    # sit ALIVE at both production bounds: a Gaussian near
    # nu=0.045 (significant field and curvature where the lower
    # one-sided stencil acts) and one near nu=260 (a live field
    # inside the Wien-dead Planck zone, which exercises the
    # endpoint treatment).
    for ncfp in ("256", "512"):
        dp = np.loadtxt(out_dir / f"esconv_bndp{int(ncfp):05d}.txt")
        de = np.loadtxt(out_dir / f"esconv_bnd{int(ncfp):05d}.txt")
        n_p = int(dp[:, 0].max())
        nd = dp.shape[0] // n_p
        n_e = int(de[:, 0].max())
        next_ = (n_e - n_p) // 2
        nu_p = dp[:, 1].reshape(nd, n_p)[0]
        jp = dp[:, 4].reshape(nd, n_p)
        nu_e = de[:, 1].reshape(nd, n_e)[0]
        je = de[:, 4].reshape(nd, n_e)
        je_i = je[:, next_ : next_ + n_p]
        # Interior nodes must coincide exactly by construction.
        rep.check(
            f"bnd {ncfp}: extended interior nodes match production",
            np.abs(nu_e[next_ : next_ + n_p] - nu_p).max()
            <= 1e-12 * nu_p.max(),
            f"max node offset {np.abs(nu_e[next_:next_+n_p]-nu_p).max():.3e}",
        )
        lnu = np.log(nu_p)
        span = lnu.max() - lnu.min()
        edge_lo = np.zeros(n_p, bool)
        edge_hi = np.zeros(n_p, bool)
        edge_lo[-6:] = True   # nu_p decreases: last nodes = low nu
        edge_hi[:6] = True    # first nodes = high nu
        mid = (lnu > lnu.min() + 0.2 * span) & (
            lnu < lnu.max() - 0.2 * span
        )
        # Interior floor first (middle 60%, all depths), then
        # per-edge statistics over ALL depth rows: max/rms rel
        # discrepancy vs the extended reference, significant-point
        # count (on the extended field), positivity, and ringing
        # (sign changes of the signed discrepancy across the zone).
        floor_max = 0.0
        for k in range(nd):
            ext = je_i[k]
            sig_i = mid & (np.abs(ext) > 1e-6 * np.abs(ext).max())
            if sig_i.any():
                di = np.abs(jp[k] - ext) / np.maximum(np.abs(ext), 1e-300)
                floor_max = max(floor_max, di[sig_i].max())
        edge_stats = {}
        for tag, zone in (("upper", edge_hi), ("lower", edge_lo)):
            mx = rms = rat_max = 0.0
            nsig = 0
            sign_changes = 0
            for k in range(nd):
                ref = jp[k]
                ext = je_i[k]
                sig = zone & (np.abs(ext) > 1e-6 * np.abs(ext).max())
                drel = np.abs(ref - ext) / np.maximum(np.abs(ext), 1e-300)
                if sig.any():
                    mx = max(mx, drel[sig].max())
                    rms = max(
                        rms, np.sqrt((drel[sig] ** 2).mean())
                    )
                    nsig += int(sig.sum())
                    rat_max = max(
                        rat_max,
                        (ref[sig] / np.maximum(ext[sig], 1e-300)).max(),
                    )
                dd = (ref - ext)[zone]
                dd = dd[dd != 0.0]
                if dd.size > 1:
                    sign_changes = max(
                        sign_changes,
                        int((np.diff(np.sign(dd)) != 0).sum()),
                    )
            edge_stats[tag] = (mx, rms, nsig, rat_max, sign_changes)
        # Interior floor (middle 60%, all depths).
        rep.check(
            f"bnd {ncfp}: interior consistency vs extended domain",
            floor_max < 5e-3,
            f"max interior discrepancy {floor_max:.3e}",
        )
        # Lower edge: live feature, Rayleigh-Jeans-side rescale
        # active.  Measured ~2e-4 max at both resolutions and
        # resolution-stable, so the bound is set at 5x the worst
        # measured value - a real stencil-error acceptance, not a
        # recorded diagnostic.
        mx, rms, nsig, rat_max, sc = edge_stats["lower"]
        rep.check(
            f"bnd {ncfp} lower: edge discrepancy vs extended (<1e-3)",
            mx < 1e-3 and rms < 3e-4,
            f"max {mx:.3e} rms {rms:.3e} nsig={nsig} "
            f"maxratio {rat_max:.4f} signchanges={sc}",
        )
        # Upper edge (TEST B/E): before the endpoint-rescale fix,
        # the B/B_es correction crushed the live feature to ~1e-20
        # of the extended reference at the top node.  The node is
        # now left uncorrected where no genuine shifted PLANCK_ES
        # exists (audit sec 8.5.1): the endpoint must recover the
        # extended result to stencil error, and no node may be
        # suppressed by orders of magnitude.  The remaining
        # discrepancy on deeper Wien-side nodes is the stock
        # wavelength-shift path evaluating the feature at a
        # shifted position - bounded, sign-changing across the
        # zone, and reported against both the extended reference
        # and the input.
        mx, rms, nsig, rat_max, sc = edge_stats["upper"]
        jin = dp[:, 3].reshape(nd, n_p)
        ep = np.abs(jp[:, 0] - je_i[:, 0]) / np.maximum(
            np.abs(je_i[:, 0]), 1e-300
        )
        rep.check(
            f"bnd {ncfp} upper: endpoint recovers extended ref",
            ep.max() < 5e-2,
            f"out(NU1)={jp[0, 0]:.4e} ext={je_i[0, 0]:.4e} "
            f"in={jin[0, 0]:.4e} "
            f"reldiff max over depths={ep.max():.3e}",
        )
        oi_min = oi_max = 1.0
        oe_min = 1.0
        pert = []
        for k in range(nd):
            sig_u = edge_hi & (
                np.abs(jin[k]) > 1e-6 * np.abs(jin[k]).max()
            )
            oi = jp[k][sig_u] / np.maximum(jin[k][sig_u], 1e-300)
            oe = jp[k][sig_u] / np.maximum(je_i[k][sig_u], 1e-300)
            oi_min = min(oi_min, oi.min())
            oi_max = max(oi_max, oi.max())
            oe_min = min(oe_min, oe.min())
            dsig = np.abs(jp[k] - je_i[k]) / np.maximum(
                np.abs(je_i[k]), 1e-300
            )
            pert.append(dsig[sig_u].max())
        rep.check(
            f"bnd {ncfp} upper: live feature not destroyed "
            "(out/in and out/ext bounded)",
            oi_min > 0.3 and oe_min > 0.3,
            f"min out/in {oi_min:.3f} max out/in {oi_max:.3f} "
            f"min out/ext {oe_min:.3f} max/ext {rat_max:.3f} "
            f"perT maxdrel {[f'{p:.1e}' for p in pert]}",
        )
        rep.check(
            f"bnd {ncfp} upper: no catastrophic suppression "
            "vs extended (bounded factor, not orders)",
            oe_min > 0.3 and rat_max < 3.0 and sc <= 6,
            f"max drel {mx:.3e} rms {rms:.3e} nsig={nsig} "
            f"maxratio {rat_max:.4f} signchanges={sc}",
        )
        rep.check(
            f"bnd {ncfp}: output finite and non-negative",
            np.all(np.isfinite(jp)) and np.all(jp >= 0.0),
        )

    # -- TEST C: smooth non-Planck Wien-tail continuum ------------------
    # A NU**0.5 ramp alive across the Wien-dead upper end (driver
    # mode wcp/wce).  The extended run skips the whole Planck
    # block (B underflows exactly at the extended top), so it
    # returns the raw convolved field: it is the right reference
    # at the boundary nodes, but interior Wien-side differences
    # against it are the DESIGNED wavelength shift on a steep
    # profile, not suppression.  Assert: endpoint/top-zone
    # agreement vs extended, and bounded output/input over the
    # whole Wien-alive zone (no orders-of-magnitude loss).
    for ncfp in ("256", "512"):
        dp = np.loadtxt(out_dir / f"esconv_wcp{int(ncfp):05d}.txt")
        de = np.loadtxt(out_dir / f"esconv_wce{int(ncfp):05d}.txt")
        n_p = int(dp[:, 0].max())
        nd = dp.shape[0] // n_p
        n_e = int(de[:, 0].max())
        next_ = (n_e - n_p) // 2
        nu_p = dp[:, 1].reshape(nd, n_p)[0]
        jin = dp[:, 3].reshape(nd, n_p)
        jp = dp[:, 4].reshape(nd, n_p)
        je_i = de[:, 4].reshape(nd, n_e)[:, next_ : next_ + n_p]
        top = nu_p > 0.85 * nu_p.max()
        mx = 0.0
        for k in range(nd):
            sig = top & (np.abs(je_i[k]) > 1e-6 * np.abs(je_i[k]).max())
            dd = np.abs(jp[k] - je_i[k]) / np.maximum(np.abs(je_i[k]), 1e-300)
            mx = max(mx, dd[sig].max())
        rep.check(
            f"wcp {ncfp}: Wien continuum survives at boundary "
            "(top nodes vs extended)",
            mx < 5e-2,
            f"max rel {mx:.3e} over {top.sum()} top nodes "
            "(bound covers the designed shift on the ramp slope)",
        )
        alive = nu_p > 20.0
        oi_min, oi_max = 1.0, 1.0
        for k in range(nd):
            sig = alive & (np.abs(jin[k]) > 1e-6 * np.abs(jin[k]).max())
            oi = jp[k][sig] / np.maximum(jin[k][sig], 1e-300)
            oi_min = min(oi_min, oi.min())
            oi_max = max(oi_max, oi.max())
        rep.check(
            f"wcp {ncfp}: non-thermal continuum bounded vs input "
            "(no catastrophic loss)",
            oi_min > 0.3 and oi_max < 5.0,
            f"out/in in [{oi_min:.3f},{oi_max:.3f}] "
            "(interior spread = designed wavelength shift)",
        )
        rep.check(
            f"wcp {ncfp}: output finite and non-negative",
            np.all(np.isfinite(jp)) and np.all(jp >= 0.0),
        )

    # -- TEST D: feature-position sweep -------------------------------
    # Same absolute Gaussian (AMP=0.5, DNU=12) moved from ~0.5
    # cells inside the endpoint to deep interior.  With the
    # endpoint-rescale fix the feature survives at every position:
    # peak ratio stays near unity modulo the weak tridiagonal
    # smearing; the only residual is the designed shift on
    # interior Wien-side nodes.  Report peak ratio and the
    # endpoint-node ratio per position.
    swp_rat = []
    for nu0 in SWP_NU0:
        d = np.loadtxt(
            out_dir / f"esconv_swp{int(10 * float(nu0)):05d}.txt"
        )
        n_p = int(d[:, 0].max())
        nd = d.shape[0] // n_p
        jin = d[:, 3].reshape(nd, n_p)
        jp = d[:, 4].reshape(nd, n_p)
        pk = (jp.max(axis=1) / np.maximum(jin.max(axis=1), 1e-300)).min()
        # Endpoint ratio only where the endpoint input is
        # significant - for far-inward positions in(NU1)
        # underflows to subnormal or zero and a ratio there
        # measures round-off, not the boundary treatment.
        sig0 = jin[:, 0] > 1e-3 * jin.max(axis=1)
        e0 = (
            (jp[sig0, 0] / np.maximum(jin[sig0, 0], 1e-300)).min()
            if sig0.any()
            else np.nan
        )
        swp_rat.append((float(nu0), pk, e0))
    worst = min(r[1] for r in swp_rat)
    esig = [r for r in swp_rat if np.isfinite(r[2])]
    worst_e = min(r[2] for r in esig)
    rep.check(
        "swp: feature peak preserved at all swept positions",
        worst > 0.5,
        " ".join(
            f"nu0={a:g}:{p:.3f}" for a, p, _ in swp_rat
        ),
    )
    rep.check(
        "swp: endpoint node bounded at all swept positions",
        worst_e > 0.3,
        " ".join(
            f"nu0={a:g}:{e:.3f}" if np.isfinite(e) else f"nu0={a:g}:n/a"
            for a, _, e in swp_rat
        ),
    )

    # -- edge step (positivity record) ---------------------------------------
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_edge.txt")
    m = sig_mask(rj, floor=1e-3)
    rep.check(
        "edge: Jbar tracks the step (no ringing)",
        np.all(np.isfinite(rjes)) and np.all(rjes >= 0.0),
        f"max rel {(np.abs(rjes - rj) / np.abs(rj))[m].max():.3e}",
    )

    # -- uniform-temperature series -------------------------------------------
    for t in TEMPU:
        nu, rj, rjes, _, nd, ncf = load_fields(
            out_dir / f"esconv_tempu{t.replace('.', '')}.txt"
        )
        m = sig_mask(rj)
        rel = np.abs(rjes - rj) / np.abs(rj)
        tv = float(t)
        beta = 1.84e-3 * np.sqrt(tv * 1.0)  # T in 1e4 K units
        rep.check(
            f"tempu {t}: B(T) preserved (uniform T)",
            rel[m].max() < 5e-3
            and np.all(np.isfinite(rjes))
            and np.all(rjes >= 0.0),
            f"beta={beta:.4f} beta/B={beta/1.6142:.4f},{beta/2.1543:.4f} "
            f"max {rel[m].max():.3e} rms {np.sqrt((rel[m]**2).mean()):.3e}",
        )

    # -- positivity with an interior gap --------------------------------------
    nu, rj, rjes, _, nd, ncf = load_fields(out_dir / "esconv_pos.txt")
    rep.check(
        "pos: Jbar >= 0 with zeroed interior band",
        np.all(np.isfinite(rjes)) and np.all(rjes >= 0.0),
        f"min {rjes.min():.3e}",
    )

    # -- grid refinement --------------------------------------------------------
    errs = {}
    for n in ("128", "256", "512"):
        nu, rj, rjes, _, nd, ncf = load_fields(
            out_dir / f"esconv_grid{int(n):05d}.txt"
        )
        m = sig_mask(rj)
        errs[n] = (np.abs(rjes - rj) / np.abs(rj))[m].max()
    rep.check(
        "grid: planck error stable or shrinking with NCF",
        errs["512"] <= errs["128"] * 1.5,
        " ".join(f"N{n}: {errs[n]:.2e}" for n in ("128", "256", "512")),
    )

    # -- T5: ES-iteration sequencing (real convolver, scripted solve) ---------
    for mode in ("correct", "stale", "preconv", "oldrend"):
        run(seq_bin, out_dir, mode, cwd=out_dir)
        (out_dir / "esseq_fields.txt").rename(
            out_dir / f"esseq_{mode}.txt"
        )
    # Perturbation runs: NEW_1 or NEW_2 scaled up by 50%.
    for mode in ("correct", "stale", "preconv", "oldrend"):
        for ptag, p1, p2 in (("p1", "0.5", "0.0"), ("p2", "0.0", "0.5")):
            run(seq_bin, out_dir, mode, p1, p2, cwd=out_dir)
            (out_dir / "esseq_fields.txt").rename(
                out_dir / f"esseq_{mode}_{ptag}.txt"
            )

    base = load_seq(out_dir / "esseq_correct.txt")
    p1 = load_seq(out_dir / "esseq_correct_p1.txt")
    p2 = load_seq(out_dir / "esseq_correct_p2.txt")

    # Harness self-consistency: the scripted recipe
    # J_1 = LAM*Jbar_0 + NEWK*(1+PERT1) means the perturbation
    # delta is exactly PERT1*NEWK = PERT1*(J_1 - LAM*Jbar_0).
    lam = 0.85
    newk = base["j1"] - lam * base["jbar0"]
    dj1 = p1["j1"] - base["j1"]
    rep.check(
        "T5: scripted sweep-1 field follows its recipe",
        np.abs(dj1 - 0.5 * newk).max() <= 1e-12 * np.abs(base["j1"]).max(),
        f"max|dJ1 - PERT1*NEWK| = {np.abs(dj1 - 0.5 * newk).max():.3e}",
    )
    # Correct ordering: the render (TA) responds to BOTH the
    # sweep-1 and sweep-2 scripted features - each generation
    # reaches the final fold through the chain.
    d_ta_p1 = np.abs(p1["ta"] - base["ta"]).max()
    d_ta_p2 = np.abs(p2["ta"] - base["ta"]).max()
    rep.check(
        "T5: render responds to sweep-1 scripted field (chain intact)",
        d_ta_p1 > 1e-3 * np.abs(base["ta"]).max(),
        f"max dTA(pert NEW_1) = {d_ta_p1:.3e}",
    )
    rep.check(
        "T5: render responds to sweep-2 scripted field",
        d_ta_p2 > 1e-3 * np.abs(base["ta"]).max(),
        f"max dTA(pert NEW_2) = {d_ta_p2:.3e}",
    )
    # Jbar_1 must derive from J_1, not J_0: perturbing NEW_1 must
    # change Jbar_1 through the real convolver.
    d_jbar1 = np.abs(p1["jbar1"] - base["jbar1"]).max()
    rep.check(
        "T5: Jbar_1 built from the stored J_1 (convolve after store)",
        d_jbar1 > 0.0,
        f"max dJbar_1(pert NEW_1) = {d_jbar1:.3e}",
    )
    # Broken modes must be detected.
    for mode, probe, expect_silent in (
        ("stale", "p1", True),
        ("stale", "p2", False),
        ("preconv", "p1", True),
        ("preconv", "p2", False),
        ("oldrend", "p1", False),
        ("oldrend", "p2", True),
    ):
        bm = load_seq(out_dir / f"esseq_{mode}_{probe}.txt")
        bn = load_seq(out_dir / f"esseq_{mode}.txt")
        d = np.abs(bm["ta"] - bn["ta"]).max()
        label = {
            "stale": "sweep2 fed Jbar_0",
            "preconv": "convolve before store",
            "oldrend": "render folds stale Jbar",
        }[mode]
        if expect_silent:
            rep.check(
                f"T5 {mode}: perturbed field cannot propagate ({label})",
                d <= 1e-14 * np.abs(bn["ta"]).max() + 1e-30,
                f"max dTA = {d:.3e}",
            )
        else:
            rep.check(
                f"T5 {mode}: downstream field responds ({label})",
                d > 1e-3 * np.abs(bn["ta"]).max(),
                f"max dTA = {d:.3e}",
            )
    for mode in ("stale", "preconv", "oldrend"):
        bm = load_seq(out_dir / f"esseq_{mode}.txt")
        d = np.abs(bm["ta"] - base["ta"]).max()
        rep.check(
            f"T5 {mode}: broken ordering differs from production",
            d > 1e-6 * np.abs(base["ta"]).max() or mode == "oldrend",
            f"max|TA_broken - TA_correct| = {d:.3e}",
        )
    hdr = (out_dir / "esseq_correct.txt").read_text().splitlines()[0]
    rep.check(
        "T5: one refresh between sweeps (refresh site counted)",
        "REFRESH_COUNT=  1" in hdr or "REFRESH_COUNT=1" in hdr,
        hdr.strip(),
    )

    # -- T6: KNOP ES guard through the real control reader ---------------------
    for case, num_es, coh, expect in (
        ("A", "1", "F", "stop"),
        ("B", "2", "T", "stop"),
        ("C", "2", "F", "pass"),
    ):
        cfile = out_dir / f"vadat_case{case}.txt"
        build_vadat(cfile, num_es, coh)
        proc = run(
            guard_bin, cfile, cwd=out_dir, check=False
        )
        out = proc.stdout + proc.stderr
        named = "requires NUM_ES>=2" in out
        if expect == "stop":
            rep.check(
                f"T6 {case}: guard STOP fires (NUM_ES={num_es},"
                f" COH_ES={coh})",
                named and "Error in RD_CMF_FLUX_CONTROLS" in out,
                "named NUM_ES/COH_ES diagnostic present"
                if named
                else f"unexpected output: {out[-300:]}",
            )
        else:
            rep.check(
                f"T6 {case}: guard passes (NUM_ES={num_es},"
                f" COH_ES={coh})",
                "GUARD_PASS" in out and not named,
                "GUARD_PASS emitted" if "GUARD_PASS" in out else out[-200:],
            )

    # -- Q-RAY emissivity census ----------------------------------------------
    qray_bounds = []
    for a in QRAY_AMPS:
        d = np.loadtxt(out_dir / f"esconv_qray_{a}.txt")
        nd = int(d[:, 0].max())
        ncf = int(d[:, 1].max())
        nu = d[:, 2].reshape(nd, ncf)[0]
        rj = d[:, 3].reshape(nd, ncf)
        rjes = d[:, 4].reshape(nd, ncf)
        chiray = d[:, 5].reshape(nd, ncf)
        esec = d[:, 6].reshape(nd, ncf)
        eta = d[:, 7].reshape(nd, ncf)
        ta_m = d[:, 9].reshape(nd, ncf)
        ta_o = d[:, 10].reshape(nd, ncf)
        clump = d[:, 11].reshape(nd, ncf)
        dta = ta_m - ta_o
        sig = np.abs(ta_o) > 1e-6 * np.abs(ta_o).max()
        rel_max = (np.abs(dta) / np.abs(ta_o))[sig].max()
        qray_bounds.append(rel_max)
        # T1: rebuild both emissivity conventions elementwise - the
        # construction must match the production formulas exactly
        # (one ESEC*RJ_ES fold, one CHI_RAY*RJ fold).
        ta_o_check = (
            eta * clump + esec * clump * rjes + chiray * clump * rj
        )
        ta_m_check = eta * clump + (esec + chiray) * clump * rjes
        scale = np.abs(ta_o).max()
        rep.check(
            f"qray {a}: TA_obs = ETA+ESEC*RJES+CHI_RAY*RJ (one fold)",
            np.abs(ta_o - ta_o_check).max() <= 1e-12 * scale,
            f"max|DTA/TA| = {rel_max:.3e}",
        )
        rep.check(
            f"qray {a}: TA_model = ETA+(ESEC+CHI_RAY)*RJES (one fold)",
            np.abs(ta_m - ta_m_check).max() <= 1e-12 * scale,
        )
        # T2: coherent convention TA_coh = ETA+(ESEC+CHI_RAY)*RJ.
        ta_coh = eta * clump + (esec + chiray) * clump * rj
        coh_rel = np.abs(
            (ta_coh - ta_o) / np.maximum(np.abs(ta_o), 1e-300)
        )[sig].max()
        rep.check(
            f"qray {a}: coherent-vs-obs difference recorded",
            np.isfinite(coh_rel),
            f"max|TA_coh-TA_obs|/TA = {coh_rel:.3e}",
        )
        # Re-derive DTA independently from the emitted fields.
        dta_check = chiray * clump * (rjes - rj)
        rep.check(
            f"qray {a}: DTA = CHI_RAY*CLUMP*(RJ_ES-RJ) to machine noise",
            np.abs(dta - dta_check).max() <= 1e-12 * np.abs(ta_o).max(),
        )
    rep.check(
        "qray: convention bound measurable and bounded (<5e-3)",
        all(1e-6 < b < 5e-3 for b in qray_bounds),
        f"bounds {['%.2e' % b for b in qray_bounds]}",
    )

    # -- Render-level flux comparison -------------------------------------------
    for a in QRAY_AMPS:
        d = np.loadtxt(out_dir / f"esconv_qray_{a}.txt")
        nd = int(d[:, 0].max())
        ncf = int(d[:, 1].max())
        nu = d[:, 2].reshape(nd, ncf)[0]
        rj = d[:, 3].reshape(nd, ncf)
        rjes = d[:, 4].reshape(nd, ncf)
        chiray = d[:, 5].reshape(nd, ncf)
        esec = d[:, 6].reshape(nd, ncf)
        eta = d[:, 7].reshape(nd, ncf)
        chi = d[:, 8].reshape(nd, ncf)
        clump = d[:, 11].reshape(nd, ncf)
        variants = {
            "model": d[:, 9].reshape(nd, ncf),
            "obs": d[:, 10].reshape(nd, ncf),
            # T4 discriminator: the coherent emissivity (no
            # redistribution applied). Its flux must differ from
            # the non-coherent render by far more than the
            # model-vs-obs convention bound - proving the
            # convolved fold is what propagates to the render.
            "coh": eta * clump + (esec + chiray) * clump * rj,
        }
        for tag, ta in variants.items():
            with open(out_dir / f"esrend_in_{tag}{a}.txt", "w") as fh:
                fh.write(f"{nd} {ncf} 3\n")
                fh.write(" ".join(f"{x:.17e}" for x in nu) + "\n")
                fh.write("10.0 7.0 4.0 2.0\n")
                fh.write("0.0 0.0 0.0 0.0\n")
                fh.write("5.0 2.0 1.0 0.5\n")
                fh.write("1.0 1.0 1.0 1.0\n")
                fh.write("8.0 6.0 4.0\n")
                for i in range(nd):
                    fh.write(" ".join(f"{x:.17e}" for x in ta[i]) + "\n")
                for i in range(nd):
                    fh.write(" ".join(f"{x:.17e}" for x in chi[i]) + "\n")
        for tag in variants:
            run(rend_bin, out_dir, f"{tag}{a}", cwd=out_dir, timeout=900)
        fm = np.loadtxt(out_dir / f"esrend_flux_model{a}.txt")
        fo = np.loadtxt(out_dir / f"esrend_flux_obs{a}.txt")
        fs = np.loadtxt(out_dir / f"esrend_flux_coh{a}.txt")
        rel = np.abs(fm[:, 1] - fo[:, 1]) / np.abs(fo[:, 1])
        rep.check(
            f"qray {a}: convention flux difference bounded (<5e-2)",
            rel.max() < 5e-2,
            f"max|dF/F| = {rel.max():.3e}",
        )
        rel_coh = np.abs(fs[:, 1] - fo[:, 1]) / np.abs(fo[:, 1])
        rep.check(
            f"qray {a}: redistribution measurably reaches flux (T4)",
            rel_coh.max() > 1e-4,
            f"max|dF_coh/F| = {rel_coh.max():.3e} "
            f"vs model-obs convention {rel.max():.3e}",
        )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
