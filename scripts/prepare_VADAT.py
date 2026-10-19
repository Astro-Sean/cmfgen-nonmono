#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

# Module-level default; overridden by --model-type CLI arg in main()
model_type = "static_snapshot"

"""
HPC-side script: update CMFGEN VADAT and write RDINR from SN_HYDRO_DATA.

Only SN_HYDRO_DATA is copied to the HPC (STELLA-to-CMFGEN conversion runs off the HPC and writes
SN_HYDRO_DATA only). This script does everything needed on the HPC:
  - Parse SN_HYDRO_DATA
  - Update VADAT (scalars, mass, abundances, profile, LSTAR)
  - Diagnose the outer hydro tail (steady wind vs shell edge) and report
    whether a wind extension is hydro-derivable; warn on IB_METH/tau mismatches
  - Write RDINR (R grid + V for RD_IN_R_GRID=T; use --no-write-rdinr to skip)
  - Write input_hydro.dat in SNIIN format (VEL_LAW=12)
  - Write T_IN and GAMMAS_IN from Temperature/Electron density when present (LTE initial)
  - Update MODEL_SPEC (ND, NP) to match SN_HYDRO_DATA

Type IIn / metal-rich CSM workflow (single-final LTE cold-start or final warm restart):
  - Single-final (LTE cold-start, default): Use VADAT_INIT as the template with your
    IIn/CSM settings (abundances, VEL_LAW=12, SN_T_OPT=USE_HYDRO, etc.). Run
    prepare_VADAT --single-final: it builds VADAT from VADAT_INIT + SN_HYDRO_DATA
    and applies the single_final profile (LTE_EST=T for LTE initial level
    populations, T fixed to hydro, alternating LAMBDA/linearization convergence).
  - Final (warm restart, time-series): After a previous day's run completes, run
    prepare_VADAT --final. It renames *OUT to *_IN, then updates VADAT with the
    final profile (LTE_EST=F, LIN_INT=T, reads *_IN populations from previous run).

  Temperature policy (three modes, cmfgen config keys):
  - thin (DEFAULT): FIX_T=F + FIX_INB_T=T + FIX_X_DPTH=<tau=2/3 photosphere
    index, auto-filled when omitted or 0>.  The fixpop_in_ba_v3 inner-boundary
    mask pins T at depths ND+1-FIX_X_DPTH..ND and lets T solve over the outer
    optically thin region where the lines form.
  - pinned (cmfgen.fix_temperature=true): FIX_T=T -- T fixed to the input
    SN_HYDRO_DATA profile at all depths; the pre-default policy.  Use this for
    A/B comparisons that must stay consistent with older pinned runs.
  - free (cmfgen.solve_temperature=true): no T stamps at all; whatever keys the
    user puts in vadat_overrides rule verbatim (irradiated ejecta with a
    PNT_SRCE inner source, or hand-picked domains).
  Do NOT use FIX_T_AUTO for the thin domain: VARFIXT pins T where
  tau_es < TAU_SCL_T (the thin region) and solves it deep, the opposite of
  the default domain.  DO_JFNK must be F for FIX_X_DPTH runs: JFNK mirrors
  only the MOD_FIX_T domain and would leave the inner T rows live.

  Electron density: FIX_NE=F for both profiles (NLTE self-consistent).
  Ne is free to allow self-consistent ionization in the optically thin CSM.

  INC_RAD_DECAYS is T when SN_HYDRO_DATA lists one or more isotope mass-fraction blocks
  (Number of isotopes > 0). CMFGEN must read NUC_DECAY_DATA in that case so RD_SN_DATA can
  match isotope names to the ISO table; with INC_RAD_DECAYS F, NUM_ISOTOPES stays zero and
  RD_SN_DATA stops with "isotope data not recognized". Otherwise INC_RAD_DECAYS is F.

Profile behavior:
- single_final: base is VADAT_INIT if present, else VADAT_INIT.bak; writes to VADAT.
- final: base is existing VADAT; overwrites VADAT with updated values.

If init fails with "J(mom) and J(ray) differ by more than 100%" then CMF_BLKBAND_V3 / "bad
initial population guesses": init profile uses N_OB_INS=12 and VFRAC_FG=0.6 for OB stability (coarser rays are usually fine).
If it still fails, try a bootstrap run with FIX_T=T and FIX_X_DPTH set to your ND (fix T at all
depths) in VADAT for the first run only; then restore init profile and rerun from the new POINT1.

Mass M = ∫ 4π r^2 ρ(r) dr from SN_HYDRO_DATA; RSTAR, RMAX, VINF, SN_AGE, etc. from hydro;
abundances (DUM/X) from depth-averaged mass fractions. LSTAR from "Energy Lum [lsun]"
in SN_HYDRO_DATA (or --lum-xg fallback). LSTAR is set after the profile so it is not overwritten.

When "Energy Lum [lsun]" is present and > 0, RSTAR and RMAX are set from the same τ=2/3
photosphere used for FIX_X_DPTH: RSTAR = R_phot (10^10 cm), RMAX = R_outer/R_phot so the
physical outer radius matches the hydro grid while the reference radius is the photosphere.
LSTAR is taken from the file (not recomputed from a blackbody). Use --no-lum-rst-geometry
to keep the old behaviour (RSTAR = inner grid radius, RMAX = R_outer/R_inner).
RD_IN_R_GRID=T so CMFGEN uses the R grid from RDINR; N_IB_INS, N_OB_INS, and N_INS are set to 0 (AMR handles boundary refinement).

Why velocity_law.f can see R(ID)=R(ID+1): Consecutive radii in RDINR come from
SN_HYDRO_DATA. They can be duplicate or nearly identical due to (1) limited
precision when SN_HYDRO_DATA is written or parsed, (2) STELLA fine zoning, or
(3) an older file written before STELLA-to-CMFGEN enforced min spacing. prepare_VADAT
enforces a minimum spacing when writing RDINR (MIN_DR_RDINR_1E10CM); velocity_law.f
also guards with DR_MIN so CMFGEN does not crash (SIGFPE).
"""

import argparse
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from log_format import Tee, print_header, print_subheader
except ImportError:
    try:
        from scripts.log_format import Tee, print_header, print_subheader
    except ImportError:
        # Fallback: keep the script usable if log_format is not importable.
        def print_header(title: str = "") -> None:
            print(f"\n=== {title} ===")

        def print_subheader(title: str = "") -> None:
            print(f"\n--- {title} ---")

        class Tee:
            def __init__(self, stdout, f):
                self._stdout = stdout
                self._file = f

            def write(self, s):
                self._stdout.write(s)
                self._file.write(s)
                self._file.flush()

            def flush(self):
                self._stdout.flush()
                self._file.flush()

# ---------------- relativistic config helper ----------------
def get_rel_config() -> dict:
    """Return static-only solver config for static_snapshot models.

    All models are STATIC snapshots with non-monotonic velocity.
    VEL_LAW=11 reads R/V/SIGMA from SN_HYDRO_DATA via SET_RV_HYDRO_MODEL_V3.
    No relativistic effects: REL_OBS=F, REL_CMF=F.
    No time evolution: USE_DJDT_RTE=F, DO_DDT=F, DO_CO_MOV_DDT=F.
    No radioactive decay heating: SN_MODEL=.FALSE. so cmfgen_sub.f:2972 is skipped.

    USE_LAM_ES=T: Lambda iteration (J=TC formal solution).  This bypasses
    the moment solver (MOM_J_CMF_V11) which clamps SIGMA to >= -0.999.
    The formal solver FG_J_CMF_NONMONO_V1 uses SIGMA_TRUE (unclamped) directly.
    USE_J_REL=F: No relativistic moment solver (not needed for static model).
    INCL_REL=F: no relativistic terms in transfer equation.
    INCL_ADV_TRANS=F: no advection (static snapshot).
    N_TYPE: always N_ON_J (set by apply_profile). N_ON_J is stable at velocity
    reversals where H=0 (N/J is well-defined since J>0 always). G_ONLY (N/H)
    crashes when H=0. N_ON_J is also required for OBS (cmf_flux.exe).

    Line treatment: GLOBAL_LINE=SOB (Sobolev approximation, set by enforce_sn_iin_physics).
    """
    return {
        "TS_NO": "1",               # First/static model
        "USE_DJDT_RTE": "F",
        "DO_CO_MOV_DDT": "F",
        "DO_DDT": "F",
        "USE_J_REL": "F",           # No relativistic moment solver (static model)
        "INCL_REL": "F",            # No relativistic terms in transfer eq
        "INCL_ADV_TRANS": "F",      # No advection (static snapshot)
        "USE_FRM_REL": "F",
        "USE_LAM_ES": "T",          # Lambda iteration (bypasses sigma clamp in moment solver)
        # N_TYPE is always N_ON_J (set by apply_profile, not here).
        # N_ON_J is stable at velocity reversals (H=0, J>0 always).
        "N_TYPE": "N_ON_J",        # Always N_ON_J: stable when H=0 at non-monotonic velocity.
        "CHK_J": "ABS_VAL",         # Patch negative J at shock fronts
    }

# ---------------- constants ----------------
MSUN_CGS = 1.98847e33  # g
DAY = 86400.0          # s
LSUN_ERG_S = 3.828e33  # L_sun in erg/s (for LSTAR)


def _parse_lstar_lsun(token):
    """Parse a VADAT-format LSTAR token into L_sun; None if not numeric."""
    if token is None:
        return None
    try:
        return float(str(token).replace("D", "E").replace("d", "e"))
    except (TypeError, ValueError):
        return None


def is_luminous_inner_boundary(dif_val, ib_meth_val):
    """True only when the inner boundary transports the LSTAR flux.

    DIF=T + IB_METH=DIFFUSION injects the hydro luminosity at the inner
    edge; IB_METH=PNT_SRCE injects the point-source blackbody instead
    (both require a luminous header LSTAR for RLUMST normalisation).
    ZERO_FLUX, HOLLOW and SCHUSTER carry no net flux, so a luminous
    model converges with no energy input (X_SOB 40d radiated ~9e39 erg/s
    against a header LSTAR of 8.8e41 erg/s; L=0 at both boundaries all
    1692 iterations).
    """
    ib = str(ib_meth_val).strip().upper()
    if "PNT_SRCE" in ib:
        return True
    return str(dif_val).strip().upper().startswith("T") and "DIFFUSION" in ib


def check_inner_boundary_luminosity(
    dif_val,
    ib_meth_val,
    tau_inner,
    lstar_lsun,
    allow_dead=False,
    line_treatment="BLANK",
):
    """Decide whether the inner boundary choice may discard LSTAR.

    Returns (fatal, messages).  fatal is True when a luminous model
    (lstar_lsun > 0) on an optically thick inner edge selects a
    flux-discarding boundary without cmfgen.allow_dead_inner_boundary.
    Callers print messages and abort when fatal is True.
    """
    if is_luminous_inner_boundary(dif_val, ib_meth_val):
        return False, []
    thresh = 25.0 if str(line_treatment).strip().upper() == "SOB" else 50.0
    if tau_inner is None or tau_inner < thresh:
        return False, []
    msgs = [
        f"  {dif_val}/{ib_meth_val} inner boundary with "
        f"tau_inner={tau_inner:.1f} >= {thresh:.0f}: the inner edge is "
        "optically thick but this boundary carries no flux.",
        "  DIFFUSION (DIF=T, IB_METH=DIFFUSION) is the standard choice "
        "for feeding the LSTAR flux through the boundary.",
    ]
    if lstar_lsun is not None and lstar_lsun > 0.0 and not allow_dead:
        msgs.append(
            f"  LSTAR={lstar_lsun:.3e} L_sun is nonzero: a "
            "flux-discarding boundary leaves the model with no energy "
            "source, and the populations converge to a dead solution "
            "radiating only stored heat."
        )
        msgs.append(
            "  Set DIF=T + IB_METH=DIFFUSION to feed the luminosity, or "
            "cmfgen.allow_dead_inner_boundary: true to run the dead "
            "model deliberately."
        )
        return True, msgs
    return False, msgs


# Keys that MUST be present in VADAT when either DO_NOMONO_KNOP=T or
# DO_NOMONO_SOB=T.  Keep in sync with
# run_model_grid.py:enforce_cmfgen_sn_iin_settings().
#
# Two non-monotonic modes are supported (selected by
# cmfgen.non_mono_mode in the YAML config):
#
#   KNOP  — Knop et al. (2009) CMF formal solver (FG_J_CMF_NONMONO_V1).
#           Uses BLANK line treatment (continuum + line transfer together).
#           The non-monotonic velocity field is handled directly in the
#           continuum radiation transfer.  This is the default.
#           Maps to: DO_NOMONO_KNOP=T, DO_NOMONO_SOB=F, GLOBAL_LINE=BLANK
#
#   SOB   - Sobolev escape-probability line treatment (local |dv/dr|).
#           Uses the standard CMF formal solver (FG_J_CMF_V13) with
#           Sobolev line treatment (SOB).  The non-monotonicity is handled
#           at the population-solver level via SIGMA_TRUE and the
#           non-monotonic BA solver (SOLVEBA_NONMONO_V1).
#           Maps to: DO_NOMONO_KNOP=F, DO_NOMONO_SOB=T, GLOBAL_LINE=SOB
#
# Both modes share the non-monotonic population solver infrastructure
# (SOLVE_FOR_POPS_NONMONO_V1 / SOLVE_FOR_POPS_JFNK_V1,
# CMF_BLKBAND_NONMONO_V1, RECOMPUTE_SIGMA_TRUE, EXTEND_SIGMA_TRUE).
#
# The line-treatment keys (GLOBAL_LINE, DO_SOB_LINES, SOB_FREQ_IN_OBS,
# TRANS_*) differ between the two modes; all other forced keys are shared.
# The mode-selection keys (DO_NOMONO_KNOP, DO_NOMONO_SOB) also differ.

# Shared forced keys (identical for KNOP and SOB modes).
_NON_MONO_FORCED_SHARED = {
    'LAM_SET': 'F',
    'USE_J_REL': 'F',
    'USE_FRM_REL': 'F',
    'USE_DJDT_RTE': 'F',
    # USE_LAM_ES is NOT forced here.  For KNOP (BLANK) mode it stays
    # user-overridable: USE_LAM_ES=F is a consistent alternative since
    # COMP_J_BLANK calls MOM_J_CMF_NONMONO_V1 (SIGMA_TRUE, Knop
    # chi_hat discretization) rather than the clamped-sigma
    # MOM_J_CMF_V11.  For SOB mode it defaults to T in
    # _NON_MONO_LINE_SOB but is user-overridable since 12-Oct-2026:
    # with an explicit USE_LAM_ES=F override the Fortran warns instead
    # of STOPping, runs MOM_J_CMF_NONMONO_V1 for the J update, and
    # DO_VAR_CONT uses fixed-J + GET_MOM_ALO_DIAG (the KNOP+LAM_ES=F
    # recipe) rather than the inconsistent V11/V12 pairing.
    # SIGMA_TRUE still drives the Sobolev line path.
    'INC_GRID': 'F',
    'N_TYPE': 'N_ON_J',
    'NEG_OPAC_OPT': 'SRCE_CHK',
    'PURE_HUB': 'F',
    'VEL_LAW': '11',
    'RD_IN_R_GRID': 'T',
    # VEL_LAW=11 + RD_IN_R_GRID=T reads the R grid verbatim from RDINR,
    # so MODEL_SPEC ND must equal the RDINR point count exactly.
    # Non-zero boundary insertion raises MODEL_SPEC ND to
    # NHYDRO+N_OB_INS+N_IB_INS and the reader aborts with
    # 'NDOLD and ND are not equal in RDINR'.
    'N_OB_INS': '0',
    'N_IB_INS': '0',
    'FG_OPT': 'INT/INS',
    'USE_FIXED_J': 'F',
    'PP_MOD': 'F',
    'PP_NOV': 'F',
    'AUTO_SMOOTH': 'F',
    # NON_MONO_SOLVER is kept for backwards compatibility with the Fortran
    # code, which derives it from DO_NOMONO_KNOP .OR. DO_NOMONO_SOB.
    # Setting it here ensures the VADAT has the keyword even if the
    # Fortran reader encounters it before the new keywords.
    'NON_MONO_SOLVER': 'T',
    # 2026-09-13: MAX_CHNG removed from the forced set.  It is a run-abort
    # safety cap, not a solver-routing requirement; with FIX_NE=F on a cold
    # (~1500 K) snapshot the first damped corrections are ~1E+170% and a
    # forced 1D100 aborts every iteration-1 (CMFGEN STOPs before writing
    # RVTJ, the job exits, and the afterok OBS fires on missing input).
    # The profile defaults still supply 1.0D150; a deliberate
    # vadat_overrides MAX_CHNG now wins.
    # 2026-09-14: LAM_SCALE_OPT was forced to LIMIT to prevent unbounded
    # LAMBDA corrections at the optically thin inner boundary.
    # 2026-09-14 (revised): LAM_SCALE_OPT removed from the forced set.
    # The user can now control this via vadat_overrides.  With
    # LAM_SCALE_OPT=NONE, MAX_LAM_COR (MAX_LAM in VADAT) is the sole
    # limiter on LAMBDA iteration jumps.  Set MAX_LAM=1.0D10 in
    # vadat_overrides to allow large jumps (LTE->NLTE transition).
    # Set LAM_SCALE_OPT=LIMIT in vadat_overrides to re-enable the
    # +-0.999 per-level clip if runaway corrections occur.
}

# Per-mode line-treatment keys.
_NON_MONO_LINE_SOB = {
    'DO_NOMONO_KNOP': 'F',
    # dev7 spelling: V7 reads DO_NONMONO_KNOP (V6 reads the
    # DO_NOMONO_KNOP above; each tree ignores the other's key).
    'DO_NONMONO_KNOP': 'F',
    'DO_NOMONO_SOB': 'T',
    # dev7 spelling: the V7 Fortran reads DO_NONMONO_SOB (V6 reads the
    # DO_NOMONO_SOB above; each tree ignores the other's key).
    'DO_NONMONO_SOB': 'T',
    # USE_LAM_ES=T is the validated default (J from the FG ray solve).
    # USE_LAM_ES=F selects stock MOM_J_CMF_V11 on clamped SIGMA plus
    # the VAR_MOM_J_CMF_V12 continuum linearization -- opt-in endgame
    # variant (12-Oct-2026); an explicit vadat_overrides key wins.
    'USE_LAM_ES': 'T',
    # SOL_METH=DIAG (NUM_BNDS=1) is the SOB default.  (The earlier claim
    # that TRI selects USE_SELECTIVE_JVAR/VAR_MOM_J_CMF_NONMONO_V1 is
    # stale: neither exists in this tree, and do_var_cont.f zeroes
    # TX/TVX under DO_NONMONO_SOB so the BA matrix is block-diagonal
    # either way.)  depth_independent_solver=false opts into
    # TRI/NUM_BNDS=3.  NUM_BNDS follows SOL_METH via
    # update_model_spec_nd().
    'SOL_METH': 'DIAG',
    'GLOBAL_LINE': 'SOB',
    'DO_SOB_LINES': 'T',
    'SOB_FREQ_IN_OBS': 'T',
    'TRANS_HI': 'SOB',
    'TRANS_HeI': 'SOB',
    'TRANS_He2': 'SOB',
    'TRANS_OI': 'SOB',
    'TRANS_O2': 'SOB',
    'TRANS_OIII': 'SOB',
    'TRANS_SkI': 'SOB',
    'TRANS_Sk2': 'SOB',
    'TRANS_SkIII': 'SOB',
    'TRANS_SkIV': 'SOB',
    'TRANS_CaI': 'SOB',
    'TRANS_Ca2': 'SOB',
    'TRANS_CaIII': 'SOB',
    'TRANS_FeI': 'SOB',
    'TRANS_Fe2': 'SOB',
    'TRANS_FeIII': 'SOB',
    'TRANS_FeIV': 'SOB',
    'TRANS_FeV': 'SOB',
}

_NON_MONO_LINE_BLANK = {
    'DO_NOMONO_KNOP': 'T',
    # dev7 spelling: V7 reads DO_NONMONO_KNOP (V6 reads the
    # DO_NOMONO_KNOP above; each tree ignores the other's key).
    'DO_NONMONO_KNOP': 'T',
    'DO_NOMONO_SOB': 'F',
    # Explicitly off in KNOP/BLANK mode so a stale DO_NONMONO_SOB=T in a
    # VADAT_INIT template cannot leak into a non-SOB run.
    'DO_NONMONO_SOB': 'F',
    # USE_LAM_ES is deliberately NOT forced in KNOP mode.  The profile
    # default is T (FG_J_CMF_NONMONO_V1 formal-solution J); an explicit
    # USE_LAM_ES=F selects MOM_J_CMF_NONMONO_V1, the non-monotonic
    # variable-Eddington moment solve (valid in V7 since 19-Sep-2026;
    # V6 likewise).  run_model_grid.py's cmfgen.use_lam_es key pairs
    # T with SOL_METH=DIAG and F with SOL_METH=TRI/NUM_BNDS=3 — the
    # verified TypeIIn_10d setups.
    'GLOBAL_LINE': 'BLANK',
    'DO_SOB_LINES': 'F',
    'SOB_FREQ_IN_OBS': 'F',
    'TRANS_HI': 'BLANK',
    'TRANS_HeI': 'BLANK',
    'TRANS_He2': 'BLANK',
    'TRANS_OI': 'BLANK',
    'TRANS_O2': 'BLANK',
    'TRANS_OIII': 'BLANK',
    'TRANS_SkI': 'BLANK',
    'TRANS_Sk2': 'BLANK',
    'TRANS_SkIII': 'BLANK',
    'TRANS_SkIV': 'BLANK',
    'TRANS_CaI': 'BLANK',
    'TRANS_Ca2': 'BLANK',
    'TRANS_CaIII': 'BLANK',
    'TRANS_FeI': 'BLANK',
    'TRANS_Fe2': 'BLANK',
    'TRANS_FeIII': 'BLANK',
    'TRANS_FeIV': 'BLANK',
    'TRANS_FeV': 'BLANK',
}


# SOB-mode recommended settings -- the Dessart, Audit & Hillier (2015,
# arXiv:1503.05463 sec 3.1) recipe as run on the m15 model by L. Dessart
# (converged VADAT extracted at cmfgen_validation/jnk_extracted/VADAT)
# and validated on dev7 job 12001805 (27-Sep-2026): stable Lambda
# contraction from a cold LTE start where the earlier stack diverged
# from it1 (job 12000611).
#
# Use CMFGEN's stock Lambda-to-Newton schedule. Residual acceptance,
# raw-demand probes and damped-full steps remain opt-in controls.
#
# These are defaults, not requirements: an explicit vadat_overrides
# value still wins (explicit_user_keys are popped before stamping,
# the same mechanism the forced table uses).
NON_MONO_SOB_RECOMMENDED = {
    'ALLOW_OL': 'F',          # geometric line overlap off (Luc's setting)
    'DC_METH': 'LTE',         # cold start from LTE pops; the non-monotonic
                            # T check re-injects 'R' via vadat_overrides,
                            # which beats this stamp
    'THK_LINE': 'T',          # thick-line outer boundary condition
    'DO_DIS': 'T',            # level dissolution (occupation probabilities)
    'INC_PEN': 'T',           # Penning ionization
    'INC_RAY': 'T',           # Rayleigh scattering
    'FIX_T': 'T',             # static snapshot: hold hydro T while pops settle
    'FIX_T_AUTO': 'F',
    'FIX_NE': 'F',            # solve Ne from charge conservation
    'FIX_IMP': 'F',           # iterate trace ions; freezing them zeroed BA corrections
    'DO_LAM_IT': 'T',
    'DO_LAM_AUTO': 'T',
    'LAM_VAL': '2.0D2',       # stock Lambda-to-Newton change threshold
    'NUM_LAM': '2',           # Lambda calls between local Newton attempts
    'BA_PROBE_LAM': '0',      # use the stock NUM_LAM cadence
    'FRESH_BA_EXIT': 'T',     # assemble a fresh matrix on return to Newton
    'N_FULL_TREL': '0',
    'T_REL_CHG': '3.0D0',
    'COMP_BA': 'T',
    'FIX_BA': '1.0D1',        # reuse the matrix only below 10% change
    'N_FIX_BA': '2',          # refresh local Newton matrix every two calls
    'MAX_LAM': '1.0D10',      # effectively unclamped Lambda correction factor
    'MAX_LIN': '1.0D1',       # manual-typical 10x cap on ordinary full steps
    'BA_CHK_FAC': '1.0D-05',
    'EPS_TERM': '5.0D-02',
    'BA_MAXCH_FAIL': '1.0D30',
    'DO_BA_RES_ACC': 'F',
    'BA_ENTER_RAW': '1.0D30',
    'BA_DAMP_CAP': '0.0D0',
    # DO_NG off by default (13-Oct-2026, X_SOB 40d): scheduled NG kicks at
    # the T-peak/interaction depth overshot by up to ~30% every ITS/NG
    # iterations, pinning MAXCH at ~0.5% above EPS_TERM once the dead-level
    # metric noise was gated. The map already contracts ~0.65x/iter without
    # acceleration. Explicit vadat_overrides DO_NG=T still enables it.
    'DO_NG': 'F',
    'CHK_NG': 'T',            # T-correction consistency check
    'NG_ABSTAIN': 'T',        # per-component hold in NG_MIT_BAND_OPT_V1
    'BEG_NG': '5.0D0',
    'IBEG_NG': '15',
    'BW_NG': '40',
    'ITS/NG': '8',
    'T_INIT_TAU': '1.5D0',
    'SIG_GAU_KMS': '3.0D3',
    # Local Dessart Sobolev treatment only; RH78 non-local coupling is a
    # separate experiment and did not cure the TypeIIn divergence.
    'DO_NONLOC_SOB': 'F',
}


def _build_non_mono_sob_recommended(out_lines, explicit_user_keys):
    explicit = set(explicit_user_keys)
    recommended = dict(NON_MONO_SOB_RECOMMENDED)
    for key in explicit:
        recommended.pop(key, None)

    lam_val = recommended.get('LAM_VAL')
    if lam_val is None:
        lam_val = _get_key_token(out_lines, 'LAM_VAL') or '2.0D2'
    if 'LAM_ENTER' not in explicit:
        recommended['LAM_ENTER'] = lam_val
    if 'LAM_EXIT' not in explicit:
        recommended['LAM_EXIT'] = lam_val
    if 'LAM_N_STABLE' not in explicit:
        recommended['LAM_N_STABLE'] = '1'
    return recommended


def build_non_mono_forced_vadat(line_treatment='BLANK'):
    """Build the NON_MONO_FORCED_VADAT dict for the given non-monotonic mode.

    Parameters
    ----------
    line_treatment : str
        'BLANK' for Knop et al. (2009) mode (default — the Knop formal
        solver handles non-monotonic velocity fields directly in the
        continuum, with BLANK line treatment) or 'SOB' for Sobolev line
        treatment (standard CMF formal
        solver with Sobolev line treatment, non-monotonicity handled at
        the population-solver level).

    Returns
    -------
    dict
        The complete forced-key table for the selected non-monotonic mode.
        Includes DO_NOMONO_KNOP/DO_NOMONO_SOB mode-selection keys.
    """
    lt = str(line_treatment).strip().upper()
    if lt == 'SOB':
        line_keys = _NON_MONO_LINE_SOB
    elif lt == 'BLANK':
        line_keys = _NON_MONO_LINE_BLANK
    else:
        raise ValueError(
            f"non_mono_line_treatment='{line_treatment}' is not valid; "
            f"must be 'SOB' or 'BLANK'")
    return {**_NON_MONO_FORCED_SHARED, **line_keys}


# Default forced table (SOB mode, for backwards compatibility with code
# that imports NON_MONO_FORCED_VADAT directly).
NON_MONO_FORCED_VADAT = build_non_mono_forced_vadat('BLANK')


def _apply_forced_vadat_keys(out_lines, forced_settings, tag="NON_MONO",
                             reason="required for non-monotonic static solver"):
    """Replace or append VADAT keys. Any TRANS_* key is forced to the
    line-treatment mode (SOB for Sobolev, BLANK for blanketing)."""
    vadat_re = re.compile(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]')
    found = set()
    new_lines = []
    force_trans_blank = any(k.startswith('TRANS_') for k in forced_settings)
    # Determine the TRANS_* catch-all value from the first TRANS_* entry.
    trans_default = 'SOB'
    for k, v in forced_settings.items():
        if k.startswith('TRANS_'):
            trans_default = str(v)
            break
    for line in out_lines:
        m = vadat_re.match(line)
        if not m:
            new_lines.append(line)
            continue
        key = m.group(2)
        if key in forced_settings:
            found.add(key)
            forced_val = forced_settings[key]
            current = m.group(1).strip()
            if current.upper() != str(forced_val).upper():
                print(f"  [{tag}] {key} forced: {current} -> {forced_val} "
                      f"({reason})")
            new_lines.append(f"{forced_val} [{key}]\n")
        elif force_trans_blank and key.startswith('TRANS_'):
            found.add(key)
            current = m.group(1).strip()
            if current.upper() != trans_default:
                print(f"  [{tag}] {key} forced: {current} -> {trans_default} "
                      f"({reason})")
            new_lines.append(f"{trans_default} [{key}]\n")
        else:
            new_lines.append(line)
    for key, val in forced_settings.items():
        if key not in found:
            new_lines.append(f"{val} [{key}]\n")
            print(f"  [{tag}] {key} added: {val} "
                  f"({reason})")
    return new_lines


def require_chg_exch_data(vadat_lines, wdir):
    """Require CHG_EXCH_DATA in the run directory when INC_CHG=T.

    RD_CHG_EXCH_V3 opens CHG_EXCH_DATA with status='OLD' in the run CWD and
    STOPs on failure, so a T setting without the file is a guaranteed crash.
    When the run dir lacks the file, stage it from the atomic tree
    ($ATOMIC/misc/chg_exch_data.*, newest dated copy); refuse only when no
    source file exists.
    """
    tok = _get_key_token(vadat_lines, 'INC_CHG')
    if tok is None or tok.strip().upper() != 'T':
        return
    data = Path(wdir) / 'CHG_EXCH_DATA'
    if data.exists():
        return
    atomic_dir = Path(os.environ.get('ATOMIC', '/u/sbrennan/atomic'))
    candidates = sorted(
        (atomic_dir / 'misc').glob('chg_exch_data.*'),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    for src in candidates:
        if src.is_file():
            os.symlink(src, data)
            print(f"  [INC_CHG] linked {data.name} -> {src}")
            return
    raise FileNotFoundError(
        f"INC_CHG=T requested but CHG_EXCH_DATA is missing from {wdir} "
        f"and no chg_exch_data.* exists under {atomic_dir}/misc. Charge "
        f"exchange needs that reaction file: stage it into the atomic "
        f"tree or set INC_CHG=F in vadat_overrides.")


# MODEL_SPEC / VADAT ion-tag prefixes -> atomic-tree element directories.
# ION_ID tags do not all use chemical symbols: Sk=Si, Tk=Ti, Nk=Ni, and
# stage suffixes mix roman (III), digit (2 -> II) and spelled (SIX -> VI)
# forms. Longest prefix wins (SkIII -> Si III, not S + kIII).
_ION_TAG_ELEM_DIR = {
    "H": "HYD", "He": "HE", "C": "CARB", "N": "NIT", "O": "OXY",
    "F": "FLU", "Ne": "NEON", "Na": "NA", "Mg": "MG", "Al": "AL",
    "Si": "SIL", "Sk": "SIL", "P": "PHOS", "S": "SUL", "Cl": "CHL",
    "Ar": "ARG", "K": "POT", "Ca": "CA", "Sc": "SCAN", "Ti": "TIT",
    "Tk": "TIT", "V": "VAN", "Cr": "CHRO", "Mn": "MAN", "Fe": "FE",
    "Co": "COB", "Ni": "NICK", "Nk": "NICK", "Ba": "BAR",
}

_DIE_SPELLED_STAGE = {
    "FOUR": "IV", "FIVE": "V", "SIX": "VI", "SEV": "VII",
    "EIGHT": "VIII", "NINE": "IX", "TEN": "X",
}
_ROMAN_STAGES = {
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
}
_DIGIT_TO_ROMAN = {
    "2": "II", "3": "III", "4": "IV", "5": "V", "6": "VI",
    "7": "VII", "8": "VIII", "9": "IX",
}


def _ion_tag_to_atomic_stage(tag):
    """Map a CMFGEN ion tag to (atomic element dir, roman stage dir).

    Handles the ION_ID suffix conventions: roman (CIII -> CARB/III),
    digit (Sk2 -> SIL/II) and spelled (OSIX -> OXY/VI, FeSEV -> FE/VII).
    Returns None when the tag cannot be resolved.
    """
    for prefix in sorted(_ION_TAG_ELEM_DIR, key=len, reverse=True):
        if not tag.startswith(prefix):
            continue
        suffix = tag[len(prefix):]
        if not suffix:
            continue
        stage = suffix.upper()
        if stage.isdigit():
            stage = _DIGIT_TO_ROMAN.get(stage)
        else:
            stage = _DIE_SPELLED_STAGE.get(stage, stage)
        if stage in _ROMAN_STAGES:
            return _ION_TAG_ELEM_DIR[prefix], stage
    return None


def require_die_data(vadat_lines, wdir):
    """Require DIE<tag> files for every enabled DIE_<tag> in VADAT.

    When a DIE_<ion> flag carries a T in either position, the ion's
    dielectronic reader opens 'DIE'//ION_ID with status='OLD' in the run
    CWD, so an enabled flag without the file is a guaranteed crash. Stage
    missing files from the atomic tree
    ($ATOMIC/<elem>/<stage>/<vintage>/die_data, newest dated copy);
    refuse only when no source file exists or the tag is unknown.
    """
    requested = []
    for line in vadat_lines:
        m = re.match(r"^\s*([^\s]+)\s+\[DIE_(\w+)\]", line)
        if not m or m.group(2) == "AS_LINE":
            continue
        if "T" in m.group(1).upper():
            requested.append(m.group(2))
    if not requested:
        return
    atomic_dir = Path(os.environ.get("ATOMIC", "/u/sbrennan/atomic"))
    for tag in sorted(set(requested)):
        data = Path(wdir) / f"DIE{tag}"
        if data.exists():
            continue
        resolved = _ion_tag_to_atomic_stage(tag)
        if resolved is None:
            raise ValueError(
                f"DIE_{tag} is enabled in VADAT but the ion tag '{tag}' "
                f"does not map to an atomic <elem>/<stage> directory. "
                f"Disable the flag or stage DIE{tag} into {wdir} manually.")
        elem_dir, stage_dir = resolved
        candidates = sorted(
            atomic_dir.glob(f"{elem_dir}/{stage_dir}/*/die_data"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        for src in candidates:
            if src.is_file():
                os.symlink(src, data)
                print(f"  [DIE_{tag}] linked {data.name} -> {src}")
                break
        else:
            raise FileNotFoundError(
                f"DIE_{tag} is enabled in VADAT but DIE{tag} is missing "
                f"from {wdir} and no {elem_dir}/{stage_dir}/*/die_data "
                f"exists under {atomic_dir}. Dielectronic data for that "
                f"ion are unavailable: set DIE_{tag}=F,F or stage the "
                f"file into the atomic tree.")


# ==============================================================================
# PHYSICS ENFORCER: Bulletproof settings for STELLA / SOB / Type IIn
# ==============================================================================
def enforce_sn_iin_physics(vadat_lines, line_treatment='BLANK',
                           user_overrides=None):
    """
    Overwrites or appends critical VADAT settings for STATIC STELLA snapshots
    of interacting transients.

    Line treatment is selected by ``line_treatment``:
      'BLANK' — Blanketing / line-overlap mode (default — the Sobolev
                approximation assumes monotonic velocity fields and is
                invalid for non-monotonic flow).
      'SOB'   - Sobolev line approximation (local |dv/dr|).

    This is a STATIC-ONLY model: no time-dependent radiative transfer,
    no time-dependent population terms, no old-time radiation-field terms.

    Configuration:
      - VEL_LAW=11 (SN hydrodynamic structure from SN_HYDRO_DATA/RDINR)
      - PURE_HUB=F (preserve non-monotonic velocity from STELLA)
      - GLOBAL_LINE=SOB or BLANK (selected by line_treatment)
      - NON_MONO_SOLVER=T (Knop et al. 2009 formal solver for sigma < -1)
      - USE_LAM_ES=T (Lambda iteration; bypasses moment solver that clamps sigma)
      - FIX_T=T (hold temperature from hydrodynamic snapshot)
      - IB_METH=DIFFUSION (for tau_inner >= 50)

    This runs AFTER all CLI arguments and YAML overrides are processed,
    ensuring no contradictory flags can break the physics.

    Settings that are user-configurable via vadat_overrides are respected:
    DIF, IB_METH, IB_STAB, MAX_LAM, MAX_LIN, BEG_NG, IBEG_NG, BW_NG,
    N_OB_INS, N_IB_INS, INC_GRID, ALL_FREQ, SCL_LN, DO_AV, ITS/NG,
    BA_CHK_FAC, T_INIT_TAU, MIN_CF, LAM_VAL, NUM_LAM, N_FIX_BA, FIX_BA,
    NEW_LINE_BA, IND_BA_METH, ST_INT, END_INT, IT_ON_T,
    T_ON_R, SCL_SL_OPAC, NON_MONO_SOLVER, FIX_IMP, NG_ABSTAIN,
    DO_BA_RES_ACC, BA_RES_GROW, BA_TRUST, BA_REJ_FRAC, BA_NREJ_LIM,
    BA_ENTER_RAW, BA_PROBE_LAM, BA_DAMP_CAP, BA_SAT_RAW, BA_SAT_LIM,
    LAM_ENTER, LAM_EXIT, LAM_N_STABLE, FRESH_BA_EXIT, N_FULL_TREL,
    T_REL_CHG.
    FIX_NE is NOT user-configurable — it is forced to F for all models
    (Ne solved from charge conservation; enforced in main()).
    AUTO_SMOOTH is NOT user-configurable — it is forced to F for all models
    (population smoothing can destabilize non-monotonic velocity models).
    """
    lt = str(line_treatment).strip().upper()
    if lt == 'SOB':
        _line_defaults = _NON_MONO_LINE_SOB
        _line_label = 'SOBOLEV (SOB)'
    elif lt == 'BLANK':
        _line_defaults = _NON_MONO_LINE_BLANK
        _line_label = 'BLANKETING (BLANK)'
    else:
        raise ValueError(
            f"line_treatment='{line_treatment}' is not valid; "
            f"must be 'SOB' or 'BLANK'")

    # Default physics settings for STATIC SN models.
    # Line-treatment keys are mode-dependent; all others are shared.
    # These are applied UNLESS the user has explicitly set them via
    # vadat_overrides (in which case the value in the VADAT lines
    # already reflects the user's choice and should be preserved).
    defaults = {
        'VEL_LAW': '11',              # SN hydrodynamic structure reader
        'PURE_HUB': 'F',              # Do NOT force Hubble flow; preserve STELLA velocity
        'T_ON_R': 'T',
        'DC_METH': 'TR',
        'INC_GRID': 'F',
        'ALL_FREQ': 'F',
        'IB_METH': 'DIFFUSION',       # For tau_inner >= 50
        'DIF': 'T',                   # Diffusion BC (consistent with IB_METH=DIFFUSION)
        'USE_LAM_ES': 'T',            # Lambda iteration (bypasses moment solver sigma clamp)
        'N_TYPE': 'N_ON_J',
        'SOL_METH': 'DIAG',
        'REL_OBS': 'F',
        'REL_CMF': 'F',
        # Static-only: no time-dependent terms
        'DO_DDT': 'F',
        'INC_ADV': 'F',
        'INC_AD': 'F',
        'INCL_DJDT': 'F',
        'USE_DJDT_RTE': 'F',
        'USE_Dr4JDT': 'F',
        # Line treatment: mode-dependent (SOB or BLANK).
        'LAM_SET': 'F',              # Keep GLOBAL_LINE active; suppress overlap warning
        # Numerical stability
        'COH_ES': 'T',                # Coherent ES for convergence
        'CHK_L_POS': 'T',
        'NEG_OPAC_OPT': 'SRCE_CHK',
        'METHOD': 'ZERO',
        'THK_CONT': 'F',              # Thin outer continuum BC — wind models are optically thin at R[ND]
        # SN line cooling / optical depth
        'SCL_SL_OPAC': 'T',
        'SCL_LN': 'F',
        'NEW_LINE_BA': 'T',
        'IND_BA_METH': '45',
        # No clumping
        'DO_CL': 'F',
        # Initial model settings
        'LIN_INT': 'F',               # New hydrodynamic structure
        'FIX_T': 'T',                 # Hold temperature from hydro
        'FIX_T_AUTO': 'F',            # Do not auto-adjust T (static initial)
        'DO_LAM_IT': 'T',
        'DO_LAM_AUTO': 'T',
        'MAX_LAM': '1.0D10',  # max Lambda-iteration correction factor (1e10 = effectively unclamped)
        # DO_NG off by default (13-Oct-2026, X_SOB 40d): scheduled NG kicks
        # overshot at the interaction depth and pinned MAXCH above EPS_TERM
        # once the dead-level demand was gated. Opt back in via
        # vadat_overrides DO_NG=T.
        'DO_NG': 'F',
        'CHK_NG': 'T',
        'NG_ABSTAIN': 'T',  # hold unsafe components rather than vetoing the depth
        'BEG_NG': '5.0D0',  # start Ng acceleration when max % change falls below this
        'IBEG_NG': '15',
        'BW_NG': '40',
        'ITS/NG': '8',
        'MAX_LIN': '1.0D1',  # manual-typical 10x cap on ordinary full steps
        'AUTO_SMOOTH': 'F',          # Forced off — smoothing destabilizes non-monotonic velocity models
        # Required by the solver (GET_KEY_STRING aborts without it); 0
        # matches the running MOM/SOB configs and the recommended default.
        'N_FULL_TREL': '0',
    }
    # Merge mode-dependent line-treatment keys into defaults.
    defaults.update(_line_defaults)

    # Keys that the user can override via vadat_overrides.
    # If these are already set in the VADAT lines (by the override step),
    # we preserve them instead of forcing the default.
    user_overridable = {
        'DIF', 'IB_METH', 'IB_STAB', 'MAX_LAM', 'MAX_LIN', 'MAX_dT',
        'BEG_NG', 'IBEG_NG', 'BW_NG', 'N_OB_INS', 'N_IB_INS',
        'INC_GRID', 'ALL_FREQ', 'SCL_LN', 'DO_AV', 'ITS/NG',
        'BA_CHK_FAC', 'T_INIT_TAU', 'MIN_CF', 'LAM_VAL', 'NUM_LAM',
        'N_FIX_BA', 'FIX_BA', 'NEW_LINE_BA', 'IND_BA_METH', 'COMP_BA',
        'ST_INT', 'END_INT', 'IT_ON_T', 'T_ON_R', 'SCL_SL_OPAC',
        'NON_MONO_SOLVER', 'FIX_NE', 'FIX_IMP', 'DO_NG', 'NG_ABSTAIN',
        'DO_LAM_IT', 'DO_LAM_AUTO', 'USE_LAM_ES',
        'DO_BA_RES_ACC', 'BA_RES_GROW', 'BA_TRUST', 'BA_REJ_FRAC',
        'BA_NREJ_LIM', 'BA_ENTER_RAW', 'BA_PROBE_LAM', 'BA_DAMP_CAP',
        'BA_SAT_RAW', 'BA_SAT_LIM', 'LAM_ENTER', 'LAM_EXIT',
        'LAM_N_STABLE', 'FRESH_BA_EXIT', 'N_FULL_TREL', 'T_REL_CHG',
        'OB_RAT', 'THK_CONT', 'STOP_IF_BP', 'REV_RGRID',
        'LTE_EST', 'FIX_T', 'FIX_T_AUTO',
        'SOL_METH', 'VEL_LAW',
        'VFRAC_MOM', 'VFRAC_FG',
        'LAM_SCALE_OPT',
        # JFNK solver controls (V4, 2026-10)
        'DO_JFNK', 'USE_JFNK', 'JFNK_NITS', 'JFNK_M', 'JFNK_MAXGMRES',
        'JFNK_GTOL',
        'JFNK_TOL', 'JFNK_DEL', 'JFNK_MNDL', 'JFNK_MXDL',
        'JFNK_EPS', 'JFNK_FIXT', 'JFNK_VERB', 'JFNK_WARMUP', 'MAX_JFNK',
        'JFNK_RIDGE', 'JFNK_REFINE', 'JFNK_LEFTPC', 'JFNK_PC_RF',
        'JFNK_LAM_THRESH', 'JFNK_PRECOND', 'JFNK_RSCALE', 'JFNK_MAX_RESTART',
        'JFNK_LS', 'JFNK_LS_MAX', 'JFNK_LS_BETA', 'JFNK_LS_C', 'JFNK_CAUCHY',
        'JFNK_DIAG', 'JFNK_USE_BANDED', 'JFNK_SELECTIVE_JVAR',
        'JFNK_ROW_TINY', 'JFNK_FRZ_GATE', 'JFNK_COR_GATE',
        # NOTE: GLOBAL_LINE, DO_SOB_LINES, PURE_HUB are NOT user-overridable
        # when NON_MONO_SOLVER=T — they are forced by the non-monotonic
        # enforcement block below.
    }

    # Keys the user set explicitly via vadat_overrides keep their value —
    # even for keys outside user_overridable (e.g. DO_NOMONO_KNOP, REL_OBS).
    # Needed to run the standard/reference solver stack on trees whose
    # Fortran reads DO_NOMO_* (dev6), where the unconditional KNOP=T stamp
    # would otherwise activate the non-monotonic solver path.
    explicit_user_keys = set()
    if user_overrides:
        explicit_user_keys = {str(k).strip() for k in user_overrides}

    vadat_re = re.compile(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]')
    out_lines = []
    found_keys = set()

    for line in vadat_lines:
        m = vadat_re.match(line)
        if m:
            key = m.group(2)
            found_keys.add(key)
            if (key in defaults and key not in user_overridable
                    and key not in explicit_user_keys):
                out_lines.append(f"{defaults[key]} [{key}]\n")
                continue
        out_lines.append(line)

    for key, val in defaults.items():
        if (key not in found_keys and key not in user_overridable
                and key not in explicit_user_keys):
            out_lines.append(f"{val} [{key}]\n")

    # ── Non-monotonic solver enforcement ───────────────────────────────
    # When NON_MONO_SOLVER=T, the Knop et al. 2009 CMF formal solver
    # (FG_J_CMF_NONMONO_V1) handles the radiation field using the true
    # (unclamped) SIGMA_TRUE values.  This is required for non-monotonic
    # velocity fields where sigma < -1 (reverse shocks, compression regions,
    # cold dense shells).
    #
    # The static SOB configuration for non-monotonic models is:
    #   GLOBAL_LINE=SOB   (Sobolev line approximation)
    #   DO_SOB_LINES=T    (Sobolev equivalent widths)
    #   SOB_FREQ_IN_OBS=T (Sobolev nu in observed spectrum)
    #   NON_MONO_SOLVER=T (continuum formal solver — Knop et al. 2009)
    #   USE_J_REL=F        (no relativistic moment solver)
    #   USE_LAM_ES=T       (Lambda iteration; bypasses MOM_J_CMF_V11 sigma clamp)
    #   PURE_HUB=F         (preserve non-monotonic velocity from STELLA)
    #   VEL_LAW=11         (SN hydrodynamic structure reader)
    #   FIX_T=T            (hold temperature from hydro snapshot)
    #   NEG_OPAC_OPT=SRCE_CHK (handle negative opacity via source check)
    #
    # With USE_LAM_ES=T, the moment solver (MOM_J_CMF_V11) is bypassed
    # (J=TC formal solution), so its sigma clamping does not affect the
    # solution.  The formal solver FG_J_CMF_NONMONO_V1 uses SIGMA_TRUE
    # (unclamped) directly.
    #
    # These settings are ALWAYS forced when NON_MONO_SOLVER=T, regardless
    # of what the profile defaults or user overrides say, because they are
    # required for the non-monotonic solver to function correctly.
    non_mono_active = False
    for line in out_lines:
        m = vadat_re.match(line)
        if m and m.group(2) == 'NON_MONO_SOLVER':
            if m.group(1).strip().upper() == 'T':
                non_mono_active = True
            break

    if non_mono_active:
        # Settings that MUST be forced for non-monotonic solver consistency.
        # These are not user-overridable when NON_MONO_SOLVER=T because
        # incorrect values cause crashes or physically inconsistent results.
        #
        # Forced because they are required for the Knop et al. (2009) CMF
        # formal solver AND for consistent line formation.  YAML / CLI
        # cannot override these when NON_MONO_SOLVER=T:
        #   GLOBAL_LINE=SOB/BLANK + LAM_SET=F + all TRANS_*=SOB/BLANK
        #     — lines stay in the selected frequency set; LAM_SET=T would
        #       rewrite long/short wavelengths even with SOB/BLANK.
        #   DO_SOB_LINES / SOB_FREQ_IN_OBS — set per line_treatment mode
        #   FG_OPT=INT/INS — DIFF path does not implement chi_hat
        #   VEL_LAW=11 + RD_IN_R_GRID=T — SET_RV_HYDRO_MODEL_V3 reads the
        #     conversion R grid; RD_IN_R_GRID=F would rebuild R and PCHIP-
        #     smooth V, destroying sigma < -1 extrema
        #   INC_GRID=F / N_TYPE=N_ON_J — Fortran STOP otherwise.
        #   USE_LAM_ES=T is the SOB default but user-overridable to F
        #   (12-Oct-2026: Fortran warns, moment solver + VAR_MOM linearization)
        #   USE_FIXED_J / USE_FRM_REL / PP_MOD / PP_NOV = F — those paths
        #     bypass FG_J_CMF_NONMONO_V1
        #
        # DIF and IB_METH are NOT forced — they depend on CMFGEN grey tau.
        # FIX_T, DO_LAM_IT, DO_LAM_AUTO remain iteration choices.
        # SOL_METH remains DIAG vs TRI unless USE_JFNK=T (forced DIAG below).
        _forced = build_non_mono_forced_vadat(lt)
        for _k in explicit_user_keys:
            _forced.pop(_k, None)
        out_lines = _apply_forced_vadat_keys(out_lines, _forced)

        if lt == 'SOB':
            # Keep the tested stock Lambda-to-Newton schedule unless the
            # user explicitly selects an experimental BA control.
            _rec = _build_non_mono_sob_recommended(
                out_lines, explicit_user_keys)
            out_lines = _apply_forced_vadat_keys(
                out_lines, _rec, tag="NON_MONO_SOB",
                reason="Dessart+2015 Sobolev recipe default")

            # DIF=F + JBAR_W_EDD=T is a guaranteed STOP on the first BA
            # iteration: SUB_SOB_LINE_V3 -> DO_VAR_CONT -> PERTJFEAU_IBC,
            # which only has the diffusion inner BC installed (job
            # 11988551).  Route the line dJ linearization through
            # NEWPERTJD unless the user pinned JBAR_W_EDD.
            _dif_tok = _get_key_token(out_lines, 'DIF')
            if _dif_tok is None or not _dif_tok.strip().upper().startswith('T'):
                _jb_tok = _get_key_token(out_lines, 'JBAR_W_EDD')
                if 'JBAR_W_EDD' in explicit_user_keys:
                    if _jb_tok is not None and _jb_tok.strip().upper().startswith('T'):
                        print("  [NON_MONO_SOB] WARNING: JBAR_W_EDD=T with DIF=F --")
                        print("    the line-section dJ linearization reaches")
                        print("    PERTJFEAU_IBC, which STOPs under DIF=F")
                        print("    (only diffusion approximation installed).")
                else:
                    out_lines = _apply_forced_vadat_keys(
                        out_lines, {'JBAR_W_EDD': 'F'}, tag="NON_MONO_SOB",
                        reason="DIF=F routes line dJ through NEWPERTJD; "
                               "PERTJFEAU_IBC requires DIF=T")

            # LSTAR=0 leaves the inner boundary dead (DTDR=0): the m15
            # model diverged from it1 under this setting (job 12000611).
            # The header-LSTAR stamp in main() normally prevents this;
            # warn if an explicit override still leaves it zero.
            _lstar_tok = _get_key_token(out_lines, 'LSTAR')
            _lstar_zero = False
            if _lstar_tok is not None:
                try:
                    _lstar_zero = abs(float(
                        _lstar_tok.replace('D', 'E').replace('d', 'e'))) < 1.0e-30
                except ValueError:
                    pass
            if _lstar_zero:
                print("  [NON_MONO_SOB] WARNING: LSTAR=0 -- dead inner boundary")
                print("    (DTDR=0).  The m15 SOB run diverged from it1 under")
                print("    this setting; use the SN_HYDRO_DATA 'Energy Lum'")
                print("    value with DIF=T/IB_METH=DIFFUSION (Dessart+2015).")

        print(f"  [NON_MONO] Static non-monotonic {lt} configuration:")
        print("    Formal solver: FG_J_CMF_NONMONO_V1 (Knop et al. 2009)")
        if lt == 'SOB':
            print("    Line treatment: SOB (Sobolev approximation)")
        else:
            print("    Line treatment: BLANK (blanketing mode, legacy)")
        print("    FG_OPT: INT/INS (chi_hat short characteristics)")
        print("    Velocity field: arbitrary / non-monotonic (PURE_HUB=F, VEL_LAW=11, RD_IN_R_GRID=T)")
        print("    Sigma handling: allow sigma < -1 (SIGMA_TRUE preserved)")
        # Show actual solver/iteration settings (user-overridable)
        use_lam_es_val = 'T'
        sol_meth_val = 'TRIDIAG'
        n_type_val = 'N_ON_J'
        fix_t_val = 'T'
        do_lam_it_val = 'T'
        for line in out_lines:
            m = vadat_re.match(line)
            if m:
                if m.group(2) == 'USE_LAM_ES':
                    use_lam_es_val = m.group(1).strip()
                elif m.group(2) == 'SOL_METH':
                    sol_meth_val = m.group(1).strip()
                elif m.group(2) == 'N_TYPE':
                    n_type_val = m.group(1).strip()
                elif m.group(2) == 'FIX_T':
                    fix_t_val = m.group(1).strip()
                elif m.group(2) == 'DO_LAM_IT':
                    do_lam_it_val = m.group(1).strip()
        if use_lam_es_val.upper() == 'T':
            print(f"    Moment solver: bypassed (USE_LAM_ES=T)")
        else:
            print(f"    Moment solver: enabled (USE_LAM_ES=F, dual-sigma with clamped SIGMA)")
        print(f"    BA solver: {sol_meth_val} (SOL_METH={sol_meth_val})")
        print(f"    Eddington factor: {n_type_val} (N_TYPE={n_type_val})")
        if fix_t_val.upper() == 'T':
            print("    Temperature: fixed from hydro (FIX_T=T)")
        else:
            print("    Temperature: variable (FIX_T=F, solved in BA matrix)")
        if do_lam_it_val.upper() == 'T':
            print("    Iteration: pure Lambda (DO_LAM_IT=T)")
        else:
            print("    Iteration: auto-switch Lambda->BA (DO_LAM_IT=F, DO_LAM_AUTO=T)")
        print("    Time evolution: disabled (static model)")

    # JFNK GMRES preconditioner inverts the DIAG BA by default.  When
    # JFNK_USE_BANDED=T, the preconditioner uses a banded LU that
    # requires SOL_METH=TRI.  Force DIAG only when JFNK is active AND
    # banded mode is NOT requested.
    jfnk_active = False
    jfnk_banded = False
    for line in out_lines:
        m = vadat_re.match(line)
        if m and m.group(2) in ('DO_JFNK', 'USE_JFNK'):
            if m.group(1).strip().upper() == 'T':
                jfnk_active = True
        if m and m.group(2) == 'JFNK_USE_BANDED':
            if m.group(1).strip().upper() == 'T':
                jfnk_banded = True
    if jfnk_active and not jfnk_banded:
        out_lines = _apply_forced_vadat_keys(
            out_lines, {'SOL_METH': 'DIAG'}, tag="JFNK"
        )
        print("  [JFNK] SOL_METH forced to DIAG (JFNK preconditioner is DIAG-only)")
    elif jfnk_active and jfnk_banded:
        out_lines = _apply_forced_vadat_keys(
            out_lines, {'SOL_METH': 'TRI'}, tag="JFNK"
        )
        print("  [JFNK] SOL_METH forced to TRI (JFNK_USE_BANDED=T — banded LU)")

    # Non-coherent ES companions: when the effective COH_ES=F survives to
    # the written VADAT (explicit user override), the frequency grid must
    # carry the extended e.s. wings and the linearization must not mix in a
    # coherent fraction.
    coh_es_tok = None
    for line in out_lines:
        m = vadat_re.match(line)
        if m and m.group(2) == 'COH_ES':
            coh_es_tok = m.group(1).strip()
            break
    if coh_es_tok is not None and coh_es_tok.upper() == 'F':
        for _c in apply_noncoherent_es_companions(
                out_lines, explicit_user_keys=explicit_user_keys):
            print(f"  [ES_NONCOH] {_c}")

    return out_lines

# Map human-readable element symbols to CMFGEN "DUM" labels.
# Full composition for complex interacting supernova systems: H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe, Ni, Co
# CMFGEN species names from cmfgen.f source code
CMFGEN_LABELS = {
    "H": "HYD", "HE": "HE", "He": "HE",
    "C": "CARB", "N": "NIT", "O": "OXY", "Ne": "NEON",
    "Mg": "MAG", "Al": "ALUM", "Si": "SIL", "S": "SUL", "Ar": "ARG",
    "Ca": "CAL", "Sc": "SCAN", "Ti": "TIT", "Cr": "CHRO", "Mn": "MAN",
    "Fe": "IRON", "Ni": "NICK", "Co": "COB", "Na": "SOD", "P": "PHOS",
    "Cl": "CHL",
}

# (Z, mean atomic weight in amu) for DUM keys in SN_HYDRO_DATA mass-fraction headers.
# Used by write_gammas_in_from_hydro so GETELEC_V2's sum_s gamma_s * n_s matches n_e.
# Full composition for complex interacting supernova systems: H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe, Ni, Co
# CMFGEN species names from cmfgen.f source code
DUM_ATOMIC: Dict[str, Tuple[int, float]] = {
    "HYD": (1, 1.008),
    "HE": (2, 4.002602),
    "CARB": (6, 12.011),
    "NIT": (7, 14.007),
    "OXY": (8, 15.999),
    "NEON": (10, 20.180),
    "SOD": (11, 22.990),
    "MAG": (12, 24.305),
    "ALUM": (13, 26.982),
    "SIL": (14, 28.085),
    "PHOS": (15, 30.974),
    "SUL": (16, 32.06),
    "CHL": (17, 35.45),
    "ARG": (18, 39.948),
    "CAL": (20, 40.078),
    "SCAN": (21, 44.956),
    "TIT": (22, 47.867),
    "CHRO": (24, 51.996),
    "MAN": (25, 54.938),
    "IRON": (26, 55.845),
    "COB": (27, 58.933),
    "NICK": (28, 58.693),
}

AMU_G = 1.66053906660e-24  # g, unified atomic mass unit

_NUM_TOKEN_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$")
_VADAT_LINE_RE = re.compile(r"^(\s*)(.+?)(\s+)\[([^\]]+)\](.*)$")

# Which VADAT keyword to set for mass (edit if your CMFGEN expects a different key)
MASS_KEY = "MASS"   # e.g. "MEJ", "MSTAR", "MASS", etc.

# ---------- .xg file parsing (same format as rho.xg) ----------
_NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"
_XG_TIME_RE = re.compile(
    r"^\s*#?\s*(?:time|t)\s*(?:\(\s*s\s*\))?\s*[:=]\s*(" + _NUM + r")\s*$",
    re.IGNORECASE,
)


def _xg_as_float(x: str) -> float:
    return float(str(x).replace("D", "E").replace("d", "e"))


def _xg_read_pairs(lines: List[str]) -> tuple:
    a_list, b_list = [], []
    for line in lines:
        if not line or line[0] in "#!cC":
            continue
        nums = re.findall(_NUM, line)
        if len(nums) >= 2:
            a, b = _xg_as_float(nums[0]), _xg_as_float(nums[1])
            if np.isfinite(a) and np.isfinite(b):
                a_list.append(a)
                b_list.append(b)
    return np.asarray(a_list, float), np.asarray(b_list, float)


def load_xg_profiles(path: Path) -> Dict[float, tuple]:
    """Parse an .xg file (same format as rho.xg). Returns dict time_sec -> (col_a, col_b)."""
    if not path.exists():
        raise FileNotFoundError(path)
    raw_lines = [ln.replace('"', "").strip() for ln in path.read_text(errors="ignore").splitlines()]
    profiles: Dict[float, tuple] = {}
    cur_t, cur_lines = None, []
    for line in raw_lines:
        if not line:
            continue
        mt = _XG_TIME_RE.match(line)
        if mt:
            if cur_t is not None and cur_lines:
                a, b = _xg_read_pairs(cur_lines)
                if a.size:
                    profiles[cur_t] = (a, b)
            cur_t, cur_lines = _xg_as_float(mt.group(1)), []
            continue
        cur_lines.append(line)
    if cur_t is not None and cur_lines:
        a, b = _xg_read_pairs(cur_lines)
        if a.size:
            profiles[cur_t] = (a, b)
    if not profiles:
        a, b = _xg_read_pairs(raw_lines)
        if a.size:
            return {0.0: (a, b)}
    return dict(sorted(profiles.items()))


def read_lum_at_outer_from_lum_xg(lum_xg_path: Path, time_days: float) -> Optional[float]:
    """
    Read lum.xg (same format as rho.xg), pick snapshot at time_days, return luminosity at outer
    boundary (shell just behind start of grid) in erg/s. Returns None if file missing or parse fails.
    """
    try:
        profiles = load_xg_profiles(lum_xg_path)
    except Exception:
        return None
    times_sec = np.array(list(profiles.keys()), float)
    if not times_sec.size:
        return None
    idx = int(np.argmin(np.abs(times_sec / DAY - float(time_days))))
    t_sec = float(times_sec[idx])
    col1, col2 = profiles[t_sec]
    col1, col2 = np.asarray(col1, float), np.asarray(col2, float)
    if not col1.size or not col2.size:
        return None
    # Assume first column is radius (cm); outer boundary = max radius
    i_outer = int(np.argmax(col1))
    return float(col2[i_outer])


def _strip_inline_comment(line: str) -> str:
    if "!" in line:
        line = line.split("!", 1)[0]
    return line.rstrip()


def _format_like(old_token: str, new_value: float) -> str:
    old = old_token.strip()
    exp_letter = "D" if "D" in old.upper() else "E"

    sci_match = re.search(r"([eEdD])([+-]?)(\d+)$", old)
    mantissa = re.split(r"[eEdD]", old, flags=re.I)[0]
    decimals = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0

    if sci_match:
        exp_digits = len(sci_match.group(3))
        base = f"{new_value:.{decimals}E}"
        mant, exp = base.split("E")
        exp_i = int(exp)
        sign = "+" if exp_i >= 0 else "-"
        exp_abs = abs(exp_i)
        exp_str = f"{exp_abs:0{exp_digits}d}"
        return f"{mant}{exp_letter}{sign}{exp_str}"

    return f"{new_value:.{decimals}f}"


@dataclass
class HydroData:
    params: Dict[str, float]
    mass_fractions: Dict[str, np.ndarray]
    radius_cm: np.ndarray
    rho_cgs: np.ndarray
    vel_cms: np.ndarray  # cm/s, same order as radius_cm (inner->outer)
    temp_104k: Optional[np.ndarray] = None  # Temperature in 10^4 K (when present in SN_HYDRO_DATA)
    kappa_cgs: Optional[np.ndarray] = None  # Total opacity (cm^2/g) when present in SN_HYDRO_DATA
    ed_cgs: Optional[np.ndarray] = None  # Electron density (1/cm^3) when present in SN_HYDRO_DATA
    na_cgs: Optional[np.ndarray] = None  # Atom density (1/cm^3) when present in SN_HYDRO_DATA


# ---------- SN_HYDRO_DATA parsing ----------
def sn_hydro_header_num_isotopes(filepath: Path, max_lines: int = 400) -> int:
    """
    Return the integer after 'Number of isotopes:' in SN_HYDRO_DATA, or 0 if absent.
    Scans an early prefix only. Header field order is not guaranteed across writers,
    so this is used both from parse_sn_hydro and as a belt-and-suspenders check before CMFGEN.
    """
    try:
        lines = filepath.read_text().splitlines()
    except OSError:
        return 0
    for line in lines[:max_lines]:
        m = re.search(r"Number of isotopes:\s*(\d+)", line, flags=re.IGNORECASE)
        if m:
            return int(m.group(1))
    return 0


def parse_sn_hydro(filepath: Path) -> HydroData:
    from grid_mapping.sn_hydro_io import read_sn_hydro_data
    _d = read_sn_hydro_data(filepath)
    if (not np.all(np.isfinite(_d.r_1e10cm)) or np.any(_d.r_1e10cm <= 0)
            or np.any(np.diff(_d.r_1e10cm) >= 0)):
        raise ValueError("SN_HYDRO_DATA radius must be finite, positive, "
                         "and strictly decreasing (outer->inner)")
    for _name, _a in (("Temperature", _d.T_1e4K),
                      ("Density", _d.density_gcm3),
                      ("Electron density", _d.electron_density)):
        if _a is not None and (not np.all(np.isfinite(_a)) or np.any(_a <= 0)):
            raise ValueError(f"SN_HYDRO_DATA {_name} must be finite and positive")

    lines = filepath.read_text().splitlines()

    # Scan a fixed early prefix so we capture NDATA, time, and isotope count regardless
    # of line order (some files put Time(days) before Number of isotopes).
    header_limit = min(200, len(lines))
    ndata: Optional[int] = None
    time_days: Optional[float] = None
    num_isotopes_header: Optional[int] = None
    for line in lines[:header_limit]:
        if ndata is None:
            m = re.search(r"Number of data points:\s+(\d+)", line)
            if m:
                ndata = int(m.group(1))
        if time_days is None:
            m = re.search(r"Time\(days\)\s+since\s+explosion:\s+([0-9EeDd+\-\.]+)", line)
            if m:
                time_days = float(m.group(1).replace("D", "E").replace("d", "e"))
        if num_isotopes_header is None:
            m = re.search(r"Number of isotopes:\s*(\d+)", line, flags=re.IGNORECASE)
            if m:
                num_isotopes_header = int(m.group(1))
    if ndata is None or time_days is None:
        raise ValueError("Failed to parse NDATA or Time(days) from SN_HYDRO_DATA")
    
    # Ensure positive time for CMFGEN/SEDONA compatibility
    # [DEPRECATED] HERACLES may have negative times for pre-explosion phases
    # HERACLES is deprecated; the project now uses the STELLA -> CMFGEN pipeline.
    if time_days < 0:
        print(f"[SN_HYDRO_DATA Parse] Negative time detected: {time_days:.2f} days")
        print(f"  Converting to positive time since explosion bomb introduction")
        time_days = abs(time_days)
    
    # Optional: Energy Lum [lsun] (inner boundary luminosity; written by STELLA-to-CMFGEN)
    energy_lum_lsun: Optional[float] = None
    for line in lines:
        m = re.search(r"Energy\s+Lum\s*\[lsun\]\s*:\s*([0-9EeDd+\-\.]+)", line, re.IGNORECASE)
        if m:
            energy_lum_lsun = float(m.group(1).replace("D", "E").replace("d", "e"))
            break

    # Optional: Inner boundary Rosseland tau (written by convert_stella_swd_to_cmfgen.py)
    # This is computed from constant kappa_es=0.34, which matches what CMFGEN uses
    # for its opacity calculation. The STELLA NLTE kappa column can give unphysically
    # low opacities in the recombination zone (T < 5300 K), producing tau ~2.7
    # when the actual Rosseland tau is ~28.
    tau_inner_header: Optional[float] = None
    for line in lines:
        m = re.search(r"Inner\s+boundary\s+Rosseland\s+tau\s*:\s*([0-9EeDd+\-\.]+)", line, re.IGNORECASE)
        if m:
            tau_inner_header = float(m.group(1).replace("D", "E").replace("d", "e"))
            break

    # Optional: Inner boundary chi-column tau (written by convert_stella_swd_to_cmfgen.py)
    # This is the tau computed by integrating the actual chi (kappa*rho) column in the
    # file, using STELLA NLTE kappa. It will differ from tau_inner_header (constant
    # kappa_es=0.34) when the inner boundary is in cool neutral gas. Both values are
    # written for transparency; the header Rosseland tau is the one used for the
    # boundary condition.
    tau_inner_chi_header: Optional[float] = None
    for line in lines:
        m = re.search(r"Inner\s+boundary\s+chi-column\s+tau\s*:\s*([0-9EeDd+\-\.]+)", line, re.IGNORECASE)
        if m:
            tau_inner_chi_header = float(m.group(1).replace("D", "E").replace("d", "e"))
            break

    def find_start(keywords: List[str]) -> int:
        for i, ln in enumerate(lines):
            low = ln.lower()
            if all(kw in low for kw in keywords):
                return i + 1
        raise ValueError(f"Section not found for keywords {keywords}")

    def find_start_optional(keywords: List[str]) -> Optional[int]:
        for i, ln in enumerate(lines):
            low = ln.lower()
            if all(kw in low for kw in keywords):
                return i + 1
        return None

    def parse_block(start_index: int) -> np.ndarray:
        from grid_mapping.sn_hydro_io import _read_numeric_block
        return _read_numeric_block(lines, start_index, ndata)[0]

    radius_start = find_start(["radius grid"])
    vel_start = find_start(["velocity", "km/s"])
    # Density: accept "Density (gm/cm^3)" or "Density (g/cm^3)" or "Density" with common unit variants
    try:
        rho_start = find_start(["density", "gm/cm"])
    except ValueError:
        try:
            rho_start = find_start(["density", "g/cm"])
        except ValueError:
            rho_start = find_start(["density"])

    radius_cm10 = parse_block(radius_start)
    vel_kms = parse_block(vel_start)
    rho = parse_block(rho_start)

    N = min(len(radius_cm10), len(vel_kms), len(rho), ndata)
    if N <= 1:
        raise ValueError("Parsed too few hydro points from SN_HYDRO_DATA")

    radius_cm = radius_cm10[:N] * 1e10  # cm
    vel_cms = vel_kms[:N] * 1e5         # cm/s
    rho_cgs = rho[:N]                   # g/cm^3

    # File order may be outer->inner (radius decreasing). Sort to inner->outer
    # so index 0 = inner boundary, index -1 = outer boundary (for T_OB_104 etc.).
    sort_idx = np.argsort(radius_cm)
    radius_cm = radius_cm[sort_idx]
    vel_cms = vel_cms[sort_idx]
    rho_cgs = rho_cgs[sort_idx]

    temp_104k: Optional[np.ndarray] = None
    kappa_cgs: Optional[np.ndarray] = None
    # Use ["temperature", "10^4"] to avoid matching header comment "Inner boundary temperature [K]"
    temp_start = find_start_optional(["temperature", "10^4"])
    if temp_start is not None:
        temp_block = parse_block(temp_start)
        if len(temp_block) >= N:
            temp_104k = temp_block[:N][sort_idx].copy()  # same inner->outer order

    # Optional: Opacity (cm^2/g or similar) if present; used to build tau profile consistently
    # with external plotting scripts that use kappa(r)*rho(r). Accept both "Opacity" and "Kappa"
    # headers, and both "cm^2/g" and "cm^2/gm" unit spellings.
    kappa_start = find_start_optional(["opacity", "cm^2/g"])
    if kappa_start is None:
        kappa_start = find_start_optional(["opacity", "cm^2/gm"])
    if kappa_start is None:
        kappa_start = find_start_optional(["kappa", "cm^2/g"])
    if kappa_start is None:
        kappa_start = find_start_optional(["kappa", "cm^2/gm"])
    if kappa_start is not None:
        kappa_block = parse_block(kappa_start)
        if len(kappa_block) >= N:
            kappa_cgs = kappa_block[:N][sort_idx].copy()

    ed_cgs: Optional[np.ndarray] = None
    ed_start = find_start_optional(["electron density"])
    if ed_start is not None:
        ed_block = parse_block(ed_start)
        if len(ed_block) >= N:
            ed_cgs = ed_block[:N][sort_idx].copy()

    na_cgs: Optional[np.ndarray] = None
    na_start = find_start_optional(["atom density"])
    if na_start is not None:
        na_block = parse_block(na_start)
        if len(na_block) >= N:
            na_cgs = na_block[:N][sort_idx].copy()

    # geometry scalars
    R_inner_cm = float(radius_cm[0])
    R_outer_cm = float(radius_cm[-1])

    # compute mass (g) and convert to MSUN
    integrand = 4.0 * np.pi * (radius_cm ** 2) * rho_cgs  # g/cm
    if hasattr(np, "trapezoid"):
        m_g = float(np.trapezoid(integrand, radius_cm))
    else:
        m_g = float(np.trapz(integrand, radius_cm))
    m_msun = m_g / MSUN_CGS

    params: Dict[str, float] = {
        "ND": float(N),
        "time_days": float(time_days),
        "RSTAR": float(R_inner_cm / 1e10),      # in 1e10 cm
        "RMAX": float(R_outer_cm / R_inner_cm), # ratio
        "VINF": float(np.max(vel_kms[:N])),     # km/s
        "RCUB_RHO": float(R_inner_cm / 1e10),
        "RHO_CORE": float(rho_cgs[0]),
        "VEL_CORE": float(vel_cms[0] / 1e5),    # km/s
        "MASS": float(m_msun),
    }
    if num_isotopes_header is not None:
        params["num_isotopes"] = float(num_isotopes_header)
    if energy_lum_lsun is not None and np.isfinite(energy_lum_lsun):
        params["energy_lum_lsun"] = float(energy_lum_lsun)
    if tau_inner_header is not None and np.isfinite(tau_inner_header):
        params["tau_inner_header"] = float(tau_inner_header)
    if tau_inner_chi_header is not None and np.isfinite(tau_inner_chi_header):
        params["tau_inner_chi_header"] = float(tau_inner_chi_header)
    if temp_104k is not None and temp_104k.size > 0:
        params["T_MIN_104"] = float(np.min(temp_104k))
        params["T_MAX_104"] = float(np.max(temp_104k))
        params["T_OB_104"] = float(temp_104k[-1])   # outer boundary (index -1 after inner->outer sort)

    # mass fractions
    mass_fractions: Dict[str, np.ndarray] = {}
    i = 0
    while i < len(lines):
        low = lines[i].lower()
        m = re.search(r"^(\w+)\s+mass fraction", low)
        if not m:
            i += 1
            continue

        element_token = m.group(1)
        upper_tok = element_token.upper()

        dum_label: Optional[str] = None
        if upper_tok in {v.upper() for v in CMFGEN_LABELS.values()}:
            dum_label = upper_tok
        else:
            for sym, dum in CMFGEN_LABELS.items():
                if sym.upper() == upper_tok:
                    dum_label = dum.upper()
                    break

        j = i + 1
        vals: List[str] = []
        while j < len(lines) and re.search(r"[0-9]", lines[j]):
            vals.extend(re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?", lines[j]))
            j += 1

        arr = parse_block(i + 1)[:N]
        if dum_label is not None and len(arr) == N:
            mass_fractions[dum_label] = arr[sort_idx]  # apply same inner->outer sort as radius/rho/vel

        i = j

    return HydroData(
        params=params,
        mass_fractions=mass_fractions,
        radius_cm=radius_cm,
        rho_cgs=rho_cgs,
        vel_cms=vel_cms,
        temp_104k=temp_104k,
        kappa_cgs=kappa_cgs,
        ed_cgs=ed_cgs,
        na_cgs=na_cgs,
    )


def print_hydro_summary(hydro: HydroData, prefix: str = "SN_HYDRO_DATA") -> None:
    """Print a short summary of parsed or written hydro data for verification."""
    nd = int(hydro.params["ND"])
    r_cm = hydro.radius_cm
    v_kms = hydro.vel_cms / 1e5  # preserve V sign (no abs) for infall support
    r_1e10 = r_cm / 1e10
    print(f"  [{prefix}] ND={nd}  SN_AGE={hydro.params['time_days']:.4f} d  "
          f"RSTAR={hydro.params['RSTAR']:.6e} (1e10 cm)  RMAX={hydro.params['RMAX']:.6e}  "
          f"VINF={hydro.params['VINF']:.2f} km/s  MASS={hydro.params['MASS']:.6f} M_sun")
    print(f"  [{prefix}] R(1e10 cm): min={r_1e10.min():.6e} max={r_1e10.max():.6e}  "
          f"R(1)={r_1e10[0]:.6e} R(ND)={r_1e10[-1]:.6e}")
    print(f"  [{prefix}] V(km/s):    min={v_kms.min():.2f} max={v_kms.max():.2f}  "
          f"V(1)={v_kms[0]:.2f} V(ND)={v_kms[-1]:.2f}")
    print(f"  [{prefix}] rho(g/cm^3): min={hydro.rho_cgs.min():.4e} max={hydro.rho_cgs.max():.4e}")
    if hydro.temp_104k is not None and hydro.temp_104k.size > 0:
        t_k = hydro.temp_104k * 1e4
        print(f"  [{prefix}] T(K):       min={t_k.min():.0f} max={t_k.max():.0f}")
    print()


def recommended_t_exc(hydro: HydroData, t_exc_floor: float = 0.7) -> Optional[float]:
    """
    Recommend T_EXC (in 10^4 K) from the STELLA hydro snapshot.

    Per full_descr.pdf: "Setting T_EXC to approximately the photospheric
    temperature (or a little lower) can prevent overflows when first beginning
    a SN model. Populations in the outer optically thin regions are more closely
    related to the radiation field (T_EXC) than the local temperature."

    T_EXC sets a MINIMUM excitation temperature for the LTE initial guess.
    When T_local < T_EXC, populations are computed at T_EXC instead of T_local.
    This is physically motivated: in optically thin regions, the radiation
    field (set by the photosphere) determines populations, not T_local.

    We measure T_photosphere from the STELLA output as the temperature at the
    Rosseland tau=2/3 surface. This is the radiation field temperature that
    characterizes the outer optically thin regions.

    A minimum floor of ``t_exc_floor`` (default 0.7 = 7000 K) is enforced to
    ensure high-n levels have non-negligible LTE populations (e.g. He I n=5 at
    23.97 eV has exp(-24/k*5000) ~ 10^-25, still small but 15 orders better
    than at 3000K).  For cold systems (T_gas < 3000 K) where the proxy
    departure coefficients already account for the T_exc/T_gas gap, a lower
    floor (e.g. 0.3 = 3000 K) prevents the VADAT T_EXC from over-exciting
    relative to the proxy's own excitation model.
    For very hot photospheres (>20000 K), we cap at 2.0 (20000 K) to avoid
    over-ionizing ground states of abundant species.
    """
    if hydro.temp_104k is None or hydro.temp_104k.size == 0:
        print("  [T_EXC] No temperature data in hydro; using default 1.0 (10000 K)")
        return 1.0

    # Compute tau from outer boundary inward
    tau_arr = compute_tau_from_outer(hydro)
    r_cm = hydro.radius_cm
    n = len(r_cm)

    # Find photosphere: index where tau is closest to 2/3
    tau_phot = 2.0 / 3.0
    i_phot = int(np.argmin(np.abs(tau_arr - tau_phot)))

    # Temperature at the photosphere (hydro.temp_104k is inner->outer, same as r_cm)
    t_phot_104k = float(hydro.temp_104k[i_phot])
    t_phot_K = t_phot_104k * 1e4

    # Enforce reasonable bounds:
    #   Floor: t_exc_floor (default 0.7 = 7000 K) — per user request and
    #          VADAT_INIT template default.  T_EXC sets the minimum excitation
    #          temperature for LTE initial populations. 7000K gives robust LTE
    #          populations in outer optically thin regions where
    #          T_local << T_radiation.  For cold systems, lower this via the
    #          YAML ``cmfgen.t_exc_floor`` key.
    #   Cap:   2.0 (20000 K) — avoids over-ionizing abundant ground states
    t_exc = max(float(t_exc_floor), min(t_phot_104k, 2.0))

    print(f"  [T_EXC] T_photosphere (tau=2/3) = {t_phot_K:.0f} K "
          f"(depth {i_phot}/{n}); T_EXC = {t_exc:.4f} (10^4 K) = {t_exc*1e4:.0f} K "
          f"[floor={t_exc_floor*1e4:.0f} K, cap=20000 K]")
    return t_exc


def recommended_t_exc_init(hydro: HydroData, t_exc_floor: float = 0.7) -> Optional[float]:
    """
    Recommend T_EXC (in 10^4 K) for the *init* run.
    Same value as recommended_t_exc — T_EXC ≈ T_photosphere is the manual's
    recommendation for all profiles. No ramping needed: T_EXC only affects
    the initial guess (iteration 0); once iterating, the real rate equations
    at T_local take over. A converged model is identical regardless of T_EXC.
    """
    return recommended_t_exc(hydro, t_exc_floor=t_exc_floor)


# ---------- optical depth for FIX_INB_T / FIX_INB_NE / FIX_X_DPTH ----------
# CMFGEN (solve_for_pops.f, fixpop_in_ba_v3.f): when FIX_INB_T=F and FIX_X_DPTH>0,
# T is fixed only at the inner FIX_X_DPTH depths; outer T varies. Same for Ne when
# FIX_INB_NE=T. "Fix at all depths" is not applied when these are set.
TAU_PHOT = 2.0 / 3.0  # Photosphere: fix T where tau >= this; vary T where tau < this
KAPPA_ES_CGS = 0.34   # Electron scattering opacity (cm^2/g), approximate


def compute_tau_from_outer(hydro: HydroData) -> np.ndarray:
    """
    τ(r_i) = ∫ κρ ds from the outer surface inward to radius r_i.
    Hydro radius_cm is sorted ascending (index 0 = inner, -1 = outer); τ[-1] = 0, τ[0] = max.

    The STELLA NLTE kappa column can give unphysically low opacities in the
    recombination zone (T < 5300 K) — e.g. kappa=0.0005 when electron scattering
    alone gives 0.34 cm^2/g. We floor kappa at KAPPA_ES_CGS to ensure the tau
    profile matches what CMFGEN computes with its own opacity calculation.
    """
    r_cm = hydro.radius_cm
    rho = hydro.rho_cgs
    n = len(r_cm)
    if hydro.kappa_cgs is not None and len(hydro.kappa_cgs) == len(r_cm):
        kappa = np.asarray(hydro.kappa_cgs, float)
        # Floor at electron scattering opacity — STELLA NLTE kappa is unreliable
        # in the recombination zone (T < 5300 K) where H recombines and the
        # NLTE Rosseland mean drops below electron scattering.
        kappa = np.maximum(kappa, KAPPA_ES_CGS)
    else:
        kappa = np.full_like(rho, KAPPA_ES_CGS, dtype=float)
    tau = np.zeros(n, dtype=float)
    if n < 2:
        return tau
    for j in range(n - 2, -1, -1):
        dr_cm = abs(r_cm[j + 1] - r_cm[j])
        if dr_cm <= 0.0:
            tau[j] = tau[j + 1]
            continue
        chimid = 0.5 * (kappa[j + 1] * rho[j + 1] + kappa[j] * rho[j])
        tau[j] = tau[j + 1] + max(chimid * dr_cm, 0.0)
    return tau


def photospheric_radius_cm(
    r_cm: np.ndarray,
    tau: np.ndarray,
    tau_phot: float,
) -> float:
    """
    Interpolate R where τ (from outside in) equals tau_phot.
    r_cm ascending: inner at 0, outer at -1. tau[0] >= tau[1] >= ... >= tau[-1]=0.
    """
    r_cm = np.asarray(r_cm, float)
    tau = np.asarray(tau, float)
    n = int(r_cm.size)
    if n < 1:
        return 1.0
    if n < 2 or float(tau[0]) < float(tau_phot):
        return float(r_cm[0])
    if float(tau_phot) >= float(tau[0]):
        return float(r_cm[0])
    for i in range(n - 1):
        t0, t1 = float(tau[i]), float(tau[i + 1])
        if t0 >= float(tau_phot) > t1:
            r0, r1 = float(r_cm[i]), float(r_cm[i + 1])
            den = t1 - t0
            if abs(den) < 1.0e-99:
                return r0
            f = (float(tau_phot) - t0) / den
            return float(r0 + f * (r1 - r0))
    return float(r_cm[-1])


MSUN_YR_CGS = 6.3019e25  # 1 M_sun/yr in g/s (for outer-tail mdot diagnostic)


def diagnose_outer_boundary(
    hydro: HydroData,
    tau_arr: Optional[np.ndarray] = None,
    tail_frac: float = 0.15,
    min_tail: int = 5,
) -> Dict[str, float]:
    """Diagnose the outer SN_HYDRO_DATA tail for grid-boundary decisions.

    Two questions matter for the CMFGEN grid:

    1. Is the outer tail a steady wind (rho ~ r^-2, Mdot = 4 pi r^2 rho v
       roughly constant)?  Only then can a wind extension beyond the hydro
       grid (a VEL_LAW=12-style RSTAR/RMAX/MDOT construction, e.g. the
       Dessart .fl reference runs) be derived from the hydro alone.  A
       shell edge (rho falling much faster than r^-2) is a truncation of
       resolved structure — any continuation is an external model
       prescription, not a hydro result.

    2. How much optical depth does the outer tail carry?  With
       VEL_LAW=11 + RD_IN_R_GRID=T (the non-monotonic solver path) the
       CMFGEN grid IS the hydro grid, so the question is whether the grid
       edge clips anything optically significant.

    Returns a dict; printing is left to report_outer_boundary().
    """
    r = np.asarray(hydro.radius_cm, float)
    v = np.asarray(hydro.vel_cms, float)
    rho = np.asarray(hydro.rho_cgs, float)
    n = r.size
    diag: Dict[str, float] = {
        'n_tail': 0, 'slope': np.nan, 'mdot_edge_msun_yr': np.nan,
        'mdot_spread': np.nan, 'v_edge_kms': np.nan, 'tau_tail': np.nan,
        'is_wind': False, 'is_shell_edge': False,
    }
    if n < 4:
        return diag
    nt = int(np.clip(int(round(tail_frac * n)), min_tail, n - 1))
    r_t, v_t, rho_t = r[-nt:], v[-nt:], rho[-nt:]
    diag['n_tail'] = nt

    good = np.isfinite(r_t) & np.isfinite(rho_t) & (r_t > 0) & (rho_t > 0)
    if int(good.sum()) >= 3:
        diag['slope'] = float(np.polyfit(np.log(r_t[good]),
                                         np.log(rho_t[good]), 1)[0])
    if np.isfinite(v[-1]):
        diag['v_edge_kms'] = float(abs(v[-1])) / 1e5
    if np.isfinite(r[-1]) and np.isfinite(rho[-1]) and np.isfinite(v[-1]):
        diag['mdot_edge_msun_yr'] = float(
            4.0 * np.pi * r[-1] ** 2 * rho[-1] * abs(v[-1]) / MSUN_YR_CGS)
    mdot_t = 4.0 * np.pi * r_t ** 2 * rho_t * np.abs(v_t) / MSUN_YR_CGS
    mdot_t = mdot_t[np.isfinite(mdot_t) & (mdot_t > 0)]
    if mdot_t.size >= 3:
        diag['mdot_spread'] = float(mdot_t.max() / mdot_t.min())
    if tau_arr is not None and np.asarray(tau_arr).size == n:
        diag['tau_tail'] = float(np.asarray(tau_arr, float)[n - nt])

    slope = diag['slope']
    if np.isfinite(slope):
        spread = diag['mdot_spread']
        # steady wind: rho ~ r^-2 with a roughly constant mass-loss rate
        diag['is_wind'] = bool(
            -2.6 <= slope <= -1.4
            and (not np.isfinite(spread) or spread <= 3.0))
        # shell edge: density falls far faster than a wind continuation
        diag['is_shell_edge'] = bool(slope < -3.0)
    return diag


def report_outer_boundary(diag: Dict[str, float]) -> None:
    """Print the outer-boundary diagnostic and its VADAT implications."""
    if not diag.get('n_tail'):
        print("  [Outer boundary] insufficient grid points for tail analysis")
        return
    print("  [Outer boundary] SN_HYDRO_DATA tail analysis:")
    slope = diag['slope']
    slope_str = f"{slope:.2f}" if np.isfinite(slope) else "n/a"
    spread = diag['mdot_spread']
    spread_str = f", Mdot spread x{spread:.1f}" if np.isfinite(spread) else ""
    print(f"    outer {diag['n_tail']} cells: dln(rho)/dln(r) = {slope_str} "
          f"(steady wind = -2.0), v_edge = {diag['v_edge_kms']:.1f} km/s, "
          f"Mdot_edge = {diag['mdot_edge_msun_yr']:.3e} Msun/yr{spread_str}")
    if np.isfinite(diag['tau_tail']):
        print(f"    optical depth contributed by the outer tail: "
              f"tau ~ {diag['tau_tail']:.2e}")
    if diag['is_wind']:
        print("    -> tail IS a steady wind: a wind extension would be "
              "hydro-derivable (MDOT=Mdot_edge, VCORE=v_edge, RSTAR=R_outer).")
    elif diag['is_shell_edge']:
        print("    -> tail is a SHELL EDGE (rho falls faster than r^-2): no "
              "steady-wind continuation exists; grid-edge RMAX is self-consistent.")
        print("       A dense CSM wind beyond the hydro grid (e.g. a reference-run"
              " style atmosphere) is an external prescription, not a hydro result"
              " -- supply MDOT/RSTAR/RMAX via vadat_overrides only if it is part"
              " of the model definition.")
    else:
        print("    -> tail is neither a clean steady wind nor a sharp shell edge;"
              " grid-edge RMAX kept (no hydro-derived extension).")
    if np.isfinite(diag['tau_tail']) and diag['tau_tail'] > 1.0:
        print(f"    WARNING: outer tail carries tau ~ {diag['tau_tail']:.1f} > 1 --")
        print("             the grid edge may clip part of the line-formation"
              "             region; check the hydro domain is adequate.")
    # The non-monotonic solver path (VEL_LAW=11 + RD_IN_R_GRID=T) uses the
    # hydro grid verbatim, so this diagnostic is informational there; a wind
    # extension is only realizable in the VEL_LAW=12 construction.
    print("    (VEL_LAW=11 + RD_IN_R_GRID=T -- forced by the non-monotonic solver --"
          " keeps the hydro grid verbatim; a wind extension is only available"
          " via the VEL_LAW=12 construction.)")
    print()


def maybe_set_rst_from_energy_lum(
    hydro: HydroData,
    tau_phot: float = TAU_PHOT,
) -> None:
    """
    If SN_HYDRO_DATA header gave Energy Lum [lsun], set RSTAR = R_phot / 1e10 and
    RMAX = R_outer / R_phot so CMFGEN's reference radius matches the τ=tau_phot surface.
    LSTAR itself is still applied later from the same header (not overwritten here).
    """
    lum = hydro.params.get("energy_lum_lsun")
    if lum is None or not np.isfinite(lum) or float(lum) <= 0.0:
        return
    r_cm = hydro.radius_cm
    if r_cm.size < 2:
        return
    tau_arr = compute_tau_from_outer(hydro)
    r_phot_cm = photospheric_radius_cm(r_cm, tau_arr, float(tau_phot))
    r_outer_cm = float(r_cm[-1])
    r_phot_cm = max(min(r_phot_cm, r_outer_cm), float(r_cm[0]))
    hydro.params["RSTAR"] = float(r_phot_cm / 1.0e10)
    hydro.params["RMAX"] = float(r_outer_cm / max(r_phot_cm, 1.0e-200))
    print(
        f"  [SN_HYDRO_DATA] Energy Lum (LSTAR) = {float(lum):.6e} L_sun; "
        f"RSTAR/RMAX from tau={float(tau_phot):g} photosphere "
        f"(R_phot = {r_phot_cm / 1.0e10:.6e} 1e10 cm, R_outer/R_phot = {hydro.params['RMAX']:.6f})."
    )


def estimate_tau_es_and_fix_x_dpth(
    hydro: HydroData,
    tau_phot: float = TAU_PHOT,
) -> Tuple[np.ndarray, int, int]:
    """
    Estimate optical depth from outer to inner on the hydro grid.

    Preference order for opacity:
      1) If SN_HYDRO_DATA provides an opacity column (kappa_cgs in cm^2/g), use tau from kappa(r)*rho(r).
      2) Otherwise, fall back to pure electron-scattering tau with constant KAPPA_ES_CGS.

    Hydro order: index 0 = inner (core), index ND-1 = outer. Integrate tau from outer in.
    Returns (tau at each point, FIX_X_DPTH, hydro index of photosphere).
    FIX_X_DPTH = number of inner depths to fix T and Ne (tau >= tau_phot); T and Ne
    vary in outer (tau < tau_phot). CMFGEN honours this as "fix inner only".
    """
    r_cm = hydro.radius_cm
    n = len(r_cm)
    if n < 2:
        return np.zeros(n), 0, -1

    tau = compute_tau_from_outer(hydro)

    # Photosphere: index where tau is closest to tau_phot (2/3), consistent with external
    # analysis scripts that use argmin(|tau - tau_phot|).
    i_phot = int(np.argmin(np.abs(tau - float(tau_phot))))
    # Fix T at the inner region (indices 0..i_phot inclusive) where tau >= tau_phot.
    fix_x_dpth = int(np.clip(i_phot + 1, 1, n))
    return tau, fix_x_dpth, i_phot


# ---------- VADAT helpers ----------
def read_lines(path: Path) -> List[str]:
    return path.read_text().splitlines(True)

def write_lines(path: Path, lines: List[str]) -> None:
    path.write_text("".join(lines))

def backup_if_exists(path: Path, suffix: str) -> Optional[Path]:
    if not path.exists():
        return None
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    bak = path.with_name(f"{path.name}.bak_{suffix}_{ts}")
    bak.write_text(path.read_text())
    return bak

def update_keywords_in_lines(
    lines: List[str],
    key_to_value: Dict[str, object],
    numeric_format_like_old: bool = True,
) -> List[str]:
    changes: List[str] = []
    for idx, raw in enumerate(lines):
        line_nocom = _strip_inline_comment(raw.rstrip("\n"))
        m = _VADAT_LINE_RE.match(line_nocom)
        if not m:
            continue

        pre, old_token, mid, key, post = m.groups()
        key = key.strip()

        if key not in key_to_value:
            continue

        new_obj = key_to_value[key]
        if isinstance(new_obj, (float, int)):
            if numeric_format_like_old and _NUM_TOKEN_RE.match(old_token.strip()):
                new_token = _format_like(old_token, float(new_obj))
            else:
                new_token = f"{float(new_obj):.8E}".replace("E", "D")
        else:
            new_token = str(new_obj)

        if old_token.strip() == new_token.strip():
            continue

        raw_line = lines[idx].rstrip("\n")
        lb = raw_line.find("[")
        if lb < 0:
            continue
        prefix = raw_line[:lb]
        suffix = raw_line[lb:]
        leading_ws = len(prefix) - len(prefix.lstrip(" \t"))
        lines[idx] = prefix[:leading_ws] + f"{new_token} " + suffix + "\n"
        changes.append(f"{key}: {old_token.strip()} -> {new_token.strip()}")

    return changes

def find_existing_keys(lines: List[str]) -> Dict[str, int]:
    key_idx: Dict[str, int] = {}
    for i, raw in enumerate(lines):
        line = _strip_inline_comment(raw).strip()
        m = _VADAT_LINE_RE.match(line)
        if not m:
            continue
        key = m.group(4).strip()
        if key not in key_idx:
            key_idx[key] = i
    return key_idx

def set_or_append_keyword(lines: List[str], key: str, value_token: str) -> None:
    key_idx = find_existing_keys(lines)
    if key in key_idx:
        i = key_idx[key]
        raw_line = lines[i].rstrip("\n")
        lb = raw_line.find("[")
        if lb >= 0:
            prefix = raw_line[:lb]
            suffix = raw_line[lb:]
            leading_ws = len(prefix) - len(prefix.lstrip(" \t"))
            lines[i] = prefix[:leading_ws] + f"{value_token} " + suffix + "\n"
        else:
            lines[i] = f"{value_token} [{key}]\n"
    else:
        lines.append(f"{value_token} [{key}]\n")

def clean_vadat_lines(lines: List[str]) -> List[str]:
    out: List[str] = []
    seen_keys: set = set()
    for raw in lines:
        s = _strip_inline_comment(raw).strip()
        if not s:
            continue
        m = _VADAT_LINE_RE.match(s)
        if not m:
            continue
        token = m.group(2).strip()
        key = m.group(4).strip()
        if key in seen_keys:
            continue
        seen_keys.add(key)
        out.append(f"{token} [{key}]\n")
    return out


def _get_key_token(lines: List[str], key: str) -> Optional[str]:
    """Return the token (value) for key in keyword-style lines, or None if missing."""
    key_idx = find_existing_keys(lines)
    if key not in key_idx:
        return None
    raw = _strip_inline_comment(lines[key_idx[key]]).strip()
    m = _VADAT_LINE_RE.match(raw)
    return m.group(2).strip() if m else None


# Floors for the electron-scattering wing extents when COH_ES=F.  With
# non-coherent ES, line photons are redistributed far into the red wing;
# ins_line_v5.f extends the frequency grid to
#   1 - (ES_WING_EXT + R_CMF_WING_EXT*VINF) / c  on the red side,
# so ES_WING_EXT must exceed VINF or the down-scattered wing is truncated.
ES_WING_EXT_FLOOR = 12000.0   # km/s; > VINF=10472 km/s for the IIn grid
R_CMF_WING_FLOOR = 4.0        # in Vinf units; valid range is 2-20


def _vadat_float_token(tok: Optional[str]) -> Optional[float]:
    """Parse a Fortran-style numeric token (D exponent) or return None."""
    if tok is None:
        return None
    try:
        return float(str(tok).strip().upper().replace("D", "E"))
    except ValueError:
        return None


def apply_noncoherent_es_companions(
    lines: List[str],
    explicit_user_keys: Optional[set] = None,
) -> List[str]:
    """Set the companion keys a COH_ES=F (non-coherent ES) model requires.

    Call this only when the effective COH_ES value is F.  With COH_ES=F the
    electron-scattering source term uses RJ_ES -- J convolved with the e.s.
    redistribution kernel, stored in ES_J_CONV (cmfgen_sub.f,
    comp_j_blank.f).  The companions:

    - ES_WING_EXT / R_CMF_WING_EXT: floored so the frequency grid covers the
      down-scattered wings (see ES_WING_EXT_FLOOR comment).
    - MIX_COH=F: pure non-coherent ES in the linearization; ES_FAC mixing
      only applies when MIX_COH=T.

    OLD_J is deliberately not managed here: it reads a literal OLD_J_FILE
    (comp_j_conv_v2.f STOPs when the file is missing or record-incompatible)
    and nothing in the pipeline produces that file, so it stays a manual
    user choice.

    Keys listed in explicit_user_keys keep their value.
    Returns a list of change-description strings.
    """
    explicit = {str(k).strip() for k in (explicit_user_keys or set())}
    changes: List[str] = []

    def _floor(key: str, floor: float, fmt: str) -> None:
        if key in explicit:
            return
        cur = _get_key_token(lines, key)
        cur_v = _vadat_float_token(cur)
        if cur_v is None or cur_v < floor:
            new = fmt.format(floor)
            set_or_append_keyword(lines, key, new)
            changes.append(f"{key}: {cur} -> {new} (COH_ES=F)")

    _floor("ES_WING_EXT", ES_WING_EXT_FLOOR, "{:.1f}D0")
    _floor("R_CMF_WING_EXT", R_CMF_WING_FLOOR, "{:.1f}D0")

    if "MIX_COH" not in explicit:
        cur = _get_key_token(lines, "MIX_COH")
        if cur is None or cur.upper() != "F":
            set_or_append_keyword(lines, "MIX_COH", "F")
            changes.append(f"MIX_COH: {cur} -> F (pure non-coherent ES source)")

    return changes


# ---------- IN_ITS helpers ----------
def update_in_its_for_profile(wdir: Path, profile: str,
                              iteration_strategy=None) -> List[str]:
    """
    Keep IN_ITS and VADAT consistent for LAMBDA control.

    Policy (static SOB model with FIX_T=T):
      - All profiles:  DO_LAM_IT=T, DO_LAM_AUTO=T, D2_EQ_D1=F
      - D2_EQ_D1=F: depth-2 populations stay independent of depth 1.
      - DO_LAM_AUTO=T: auto-switch from LAMBDA to full iterations once corrections
        are small (MAXCH < 50%).
      - FIX_T=T + FIX_T_AUTO=F: T is fixed to the hydrodynamic snapshot for the
        entire run (static model).  This is the initial convergence test mode.
      - DO_T_AUTO=F when FIX_T=T AND FIX_T_AUTO=F (T fixed to SN_HYDRO_DATA).
      - DO_GT_AUTO=F when FIX_T=T AND FIX_T_AUTO=F (no grey T iteration needed).

    For bootstrap mode (FIX_T=F or FIX_T_AUTO=T):
      - DO_T_AUTO=T (T varies in BA matrix or auto-released after stabilization)
      - DO_GT_AUTO=T (grey temperature iteration for radiative equilibrium)

    DO_LAM_IT=T: start with LAMBDA iterations (stable cold-start from LTE).
    With DO_LAM_AUTO=T: auto-switches to full BA iterations when MAXCH < 50%.
    FIX_IMP=T (VADAT) → RD_FIX_IMP=.TRUE. → FIX_IMPURITY=RD_FIX_IMP in all
    code paths (source fix) → MOD_FIX_IMPURITY=.TRUE. → FIXPOP_IN_BA_V3 active.
    """
    in_its_path = wdir / "IN_ITS"
    if not in_its_path.exists():
        return []

    # Check if VADAT has FIX_T=F (T allowed to vary via vadat_overrides)
    # or FIX_T_AUTO=T (bootstrap: T fixed initially, auto-released when pops stabilize)
    vadat_path = wdir / "VADAT"
    fix_t_is_f = False
    fix_t_auto = False
    is_diag = False
    if vadat_path.is_file():
        for vline in read_lines(vadat_path):
            if "[FIX_T]" in vline:
                val = vline.split("[")[0].strip()
                if val.upper() == "F":
                    fix_t_is_f = True
            if "[FIX_T_AUTO]" in vline:
                val = vline.split("[")[0].strip()
                if val.upper() == "T":
                    fix_t_auto = True
            if "[SOL_METH]" in vline:
                val = vline.split("[")[0].strip()
                if val.upper() == "DIAG":
                    is_diag = True

    lines = read_lines(in_its_path)
    # Read final VADAT values for LAMBDA control so IN_ITS matches
    # vadat_overrides (e.g. DO_LAM_IT=F for a BA start).
    do_lam_val = "T"
    # Default matches the policy above (all profiles DO_LAM_AUTO=T); the
    # written VADAT normally supplies the key, this only covers a VADAT
    # that lacks it entirely.
    do_auto_val = "T"
    if vadat_path.is_file():
        for vline in read_lines(vadat_path):
            if "[DO_LAM_IT]" in vline:
                do_lam_val = vline.split("[")[0].strip().upper() or "T"
            if "[DO_LAM_AUTO]" in vline:
                do_auto_val = vline.split("[")[0].strip().upper() or "T"

    # DO_LAM_IT: start in LAMBDA (T) or BA (F) mode.
    # DO_LAM_AUTO: auto-switch between LAMBDA and BA on MAXCH.
    want_lam = do_lam_val
    want_auto = do_auto_val

    # If the existing IN_ITS already has DO_LAM_IT=F, the model has already
    # transitioned from LAMBDA to FULL/BA mode (via DO_LAM_AUTO or explicit
    # IT_SPECIFIER).  Preserve this across requeues — resetting to T would
    # force the model back into LAMBDA mode with BA-mode populations, causing
    # catastrophic divergence (MAXCH=1E+07).
    existing_lam = _get_key_token(lines, "DO_LAM_IT")
    if existing_lam is not None and existing_lam.upper() == "F":
        want_lam = "F"

    # If an explicit IT_SPECIFIER iteration strategy is present and it
    # contains any FULL / FULL(FIXT) cycles, force DO_LAM_IT=F in IN_ITS.
    # Otherwise cmfgen_sub.f re-reads IN_ITS after SPECIFY_IT_CYCLE_V3 and
    # clobbers the FULL(FIXT) setting by setting LAMBDA_ITERATION=.TRUE.
    # Only apply this when iteration_strategy is specified in the YAML —
    # when it is absent, DO_LAM_AUTO controls the transition and
    # DO_LAM_IT should stay T (start in LAMBDA mode).
    if iteration_strategy is not None:
        it_spec_path = wdir / "IT_SPECIFIER"
        if it_spec_path.is_file():
            try:
                it_text = it_spec_path.read_text()
                if re.search(r"(?i)FULL(?:\s*\(|\b)", it_text) and want_lam == "T":
                    want_lam = "F"
            except OSError:
                pass
    # DO_GT_AUTO: grey temperature iteration.  When FIX_T=T (static model with
    # temperature fixed to hydro), DO_GT_AUTO must be F — there is no need to
    # iterate the grey temperature since T is held fixed.
    want_gt = "F" if (not fix_t_is_f and not fix_t_auto) else "T"
    # DO_T_AUTO=T when FIX_T=F (T varies in BA matrix) OR FIX_T_AUTO=T (bootstrap:
    # auto-release T once populations have stabilized)
    want_t = "T" if (fix_t_is_f or fix_t_auto) else "F"
    # D2_EQ_D1=F: independent depth-2 populations (pipeline policy; pinning
    # them to depth 1 masked genuine outer-boundary corrections).
    want_d2 = "F"
    # NOTE: REVISE_R_GRID is a VADAT keyword (REV_RGRID), not an IN_ITS keyword

    before_lam = _get_key_token(lines, "DO_LAM_IT")
    before_auto = _get_key_token(lines, "DO_LAM_AUTO")
    before_gt = _get_key_token(lines, "DO_GT_AUTO")
    before_t = _get_key_token(lines, "DO_T_AUTO")
    before_d2 = _get_key_token(lines, "D2_EQ_D1")

    set_or_append_keyword(lines, "DO_LAM_IT", want_lam)
    set_or_append_keyword(lines, "DO_LAM_AUTO", want_auto)
    set_or_append_keyword(lines, "DO_GT_AUTO", want_gt)
    set_or_append_keyword(lines, "DO_T_AUTO", want_t)
    set_or_append_keyword(lines, "D2_EQ_D1", want_d2)

    changes: List[str] = []
    after_lam = _get_key_token(lines, "DO_LAM_IT")
    after_auto = _get_key_token(lines, "DO_LAM_AUTO")
    after_gt = _get_key_token(lines, "DO_GT_AUTO")
    after_t = _get_key_token(lines, "DO_T_AUTO")
    after_d2 = _get_key_token(lines, "D2_EQ_D1")

    if before_lam != after_lam:
        changes.append(f"DO_LAM_IT: {before_lam if before_lam is not None else '(missing)'} -> {after_lam}")
    if before_auto != after_auto:
        changes.append(f"DO_LAM_AUTO: {before_auto if before_auto is not None else '(missing)'} -> {after_auto}")
    if before_gt != after_gt:
        changes.append(f"DO_GT_AUTO: {before_gt if before_gt is not None else '(missing)'} -> {after_gt}")
    if before_t != after_t:
        changes.append(f"DO_T_AUTO: {before_t if before_t is not None else '(missing)'} -> {after_t}")
    if before_d2 != after_d2:
        changes.append(f"D2_EQ_D1: {before_d2 if before_d2 is not None else '(missing)'} -> {after_d2}")

    if changes:
        # Keep IN_ITS in the same "clean" format used by VADAT.
        clean = clean_vadat_lines(lines)
        write_lines(in_its_path, clean)

    return changes


# Minimum spacing between consecutive R points in RDINR (in 10^10 cm). Ensures
# velocity_law.f never sees R(ID)-R(ID+1) ≈ 0 (division-by-zero / SIGFPE).
MIN_DR_RDINR_1E10CM = 1.0e-6


def downsample_hydro_gradient_aware(
    hydro: HydroData, target_nd: int, max_dv_kms: Optional[float] = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Regrid the hydro R/V profile to target_nd points using gradient-based AMR.

    Selects which hydro R points to KEEP (not interpolate) based on local
    structure. This is critical for non-monotonic velocity profiles where
    uniform-V spacing fails (R(V) is not single-valued at velocity plateaus
    and reverse shocks).

    Strategy:
    1. Compute importance weights from |dlog(rho)/dlog(R)|, |dlog(T)/dlog(R)|,
       |dV/dR|, and |dsigma/dR| at each hydro point.
    2. Always keep boundary points and points where |sigma| > 1.
    3. Thin points from low-gradient regions (velocity plateaus) first.
    4. Ensure minimum dr/R spacing.

    Returns (r_1e10_outer_to_inner, v_kms_outer_to_inner, indices_into_hydro)
    where r_1e10 is in units of 10^10 cm and v_kms is in km/s, both ordered
    outer (largest R) to inner (smallest R).
    """
    n_hydro = len(hydro.radius_cm)

    # Hydro data is stored inner->outer (radius_cm increasing). We need outer->inner.
    r_oi = hydro.radius_cm[::-1] / 1e10  # outer->inner, 10^10 cm
    v_oi = hydro.vel_cms[::-1] / 1e5  # outer->inner, km/s, preserve V sign (no abs)
    rho_oi = hydro.rho_cgs[::-1]  # outer->inner
    t_oi = hydro.temp_104k[::-1]  # outer->inner

    if n_hydro <= target_nd:
        # No thinning needed — return as-is
        indices = np.arange(n_hydro)
        return r_oi.copy(), v_oi.copy(), indices

    # Compute gradients (outer->inner, so diff is negative for decreasing values)
    log_r = np.log(r_oi)
    log_rho = np.log(np.maximum(rho_oi, 1e-30))
    log_t = np.log(np.maximum(t_oi, 1e-4))

    # Local gradients (centered difference)
    dlogr = np.zeros(n_hydro)
    dlogrho = np.zeros(n_hydro)
    dlogt = np.zeros(n_hydro)
    dv = np.zeros(n_hydro)

    dlogr[1:-1] = (log_r[:-2] - log_r[2:]) / 2.0
    dlogrho[1:-1] = np.abs((log_rho[:-2] - log_rho[2:]) / 2.0)
    dlogt[1:-1] = np.abs((log_t[:-2] - log_t[2:]) / 2.0)
    dv[1:-1] = np.abs((v_oi[:-2] - v_oi[2:]) / 2.0)

    # Sigma = dlnV/dlnR - 1 (compute from hydro data)
    # For non-monotonic V, use abs(dV) to avoid log of negative
    log_v = np.log(np.maximum(v_oi, 0.1))  # floor at 0.1 km/s
    dlogv = np.zeros(n_hydro)
    dlogv[1:-1] = np.abs((log_v[:-2] - log_v[2:]) / 2.0)
    dsigma = np.abs(dlogv - dlogr)  # |d(dlnV/dlnR)/dR| proxy

    # Importance weight: high where any gradient is large
    # Normalize each gradient to [0, 1] and combine
    def norm(x):
        mx = x.max()
        return x / mx if mx > 0 else x

    weight = norm(dlogrho) + norm(dlogt) + norm(dv) + norm(dsigma)

    # Boost weight where |sigma| is large (reverse shocks)
    sigma = np.zeros(n_hydro)
    sigma[1:-1] = (np.log(v_oi[:-2]) - np.log(v_oi[2:])) / (log_r[:-2] - log_r[2:]) - 1.0
    weight += norm(np.abs(sigma)) * 2.0  # sigma gets extra weight

    # Always keep boundary points
    weight[0] = 1e10
    weight[-1] = 1e10

    # Always keep points where |sigma| > 10 (strong reverse shocks)
    strong_sigma = np.abs(sigma) > 10
    weight[strong_sigma] += 1e6

    # Boost (but don't protect) points where velocity changes rapidly
    rapid_dv = dv > np.median(dv[dv > 0]) * 5 if (dv > 0).any() else np.zeros(n_hydro, dtype=bool)
    weight[rapid_dv] += 1e5

    # Iteratively remove the lowest-weight point until we reach target_nd
    keep = np.ones(n_hydro, dtype=bool)
    n_current = n_hydro
    skip_points = set()  # points that can't be removed due to gap constraint

    while n_current > target_nd:
        # Find candidate points to remove (not boundary, not strong_sigma, not skipped)
        candidates = np.where(keep)[0]
        candidates = candidates[(candidates > 0) & (candidates < n_hydro - 1)]

        # Don't remove points adjacent to strong_sigma points
        protected = strong_sigma.copy()
        protected[1:] |= strong_sigma[:-1]
        protected[:-1] |= strong_sigma[1:]
        candidates = candidates[~protected[candidates]]

        # Remove skipped points from candidates
        candidates = np.array([c for c in candidates if c not in skip_points])

        if len(candidates) == 0:
            break

        # Find the lowest-weight candidate
        weights = weight[candidates]
        worst = candidates[np.argmin(weights)]

        # Check that removing this point doesn't create too large a gap
        kept_indices = np.where(keep)[0]
        idx_in_kept = np.searchsorted(kept_indices, worst)
        prev_idx = kept_indices[idx_in_kept - 1] if idx_in_kept > 0 else worst
        next_idx = kept_indices[idx_in_kept + 1] if idx_in_kept < len(kept_indices) - 1 else worst

        # Calculate gap if we remove this point
        r_prev = r_oi[prev_idx]
        r_next = r_oi[next_idx]
        gap_ratio = abs(r_prev - r_next) / r_oi[worst]

        # Don't remove if it creates a gap > 20% of R — skip this point permanently
        if gap_ratio > 0.20:
            skip_points.add(worst)
            continue

        # Don't remove if it creates a velocity jump CMFGEN cannot integrate.
        # The gap test above is purely geometric, so on its own it allows a
        # cell across which the velocity changes by many Doppler widths --
        # which is what left a 2331 km/s jump (117 Doppler widths) in the
        # SN_B_v2_M3_V5000 40 d grid. MOM_JREL_V9 compares V(I)-V(I+1)
        # against VFRAC_MOM*VDOP_VEC, so respect the same budget here.
        if max_dv_kms is not None:
            if abs(v_oi[prev_idx] - v_oi[next_idx]) > max_dv_kms:
                skip_points.add(worst)
                continue

        keep[worst] = False
        n_current -= 1

    # Extract kept points
    indices = np.where(keep)[0]
    r_out = r_oi[indices]
    v_out = v_oi[indices]

    # Ensure strict monotonicity (R decreasing outer->inner)
    for i in range(len(r_out) - 1):
        dr = r_out[i] - r_out[i + 1]
        if dr < MIN_DR_RDINR_1E10CM:
            r_out[i + 1] = r_out[i] - MIN_DR_RDINR_1E10CM

    print(f"[AMR] Downsampled {n_hydro} -> {len(indices)} points (gradient-aware)")
    print(f"[AMR] Kept {strong_sigma.sum()} strong-sigma points, {rapid_dv.sum()} rapid-dV points")
    print(f"[AMR] R range: {r_out.min():.4e} to {r_out.max():.4e} (1e10 cm)")
    print(f"[AMR] V range: {v_out.min():.4e} to {v_out.max():.4e} (km/s)")

    return r_out, v_out, indices


def downsample_hydro_velocity_aware(
    hydro: HydroData, target_nd: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Regrid the hydro R/V profile to target_nd points with uniform velocity spacing.

    Places points uniformly in velocity space (equal dV between adjacent zones) by
    interpolating R(V) from the hydro data. This prevents MOM_JREL_V9 from inserting
    extra points (which happens when V(I)-V(I+1) > VFRAC_MOM*VDOP_VEC).

    ALWAYS regrids, even when target_nd equals the hydro ND. This is necessary because
    the SN_HYDRO_DATA may already have target_nd points but with non-uniform velocity
    spacing (e.g., from STELLA→CMFGEN conversion that uses tau-based or log-R spacing,
    producing large dV in high-velocity outer regions).

    With target_nd=150 and a velocity span of ~1400 km/s, the average dV per
    interval is ~9.4 km/s, well below the typical VDOP_MIN of ~24 km/s
    (VFRAC_MOM=1.0, VDOP_VEC=23.77 km/s).

    Returns (r_1e10_outer_to_inner, v_kms_outer_to_inner, indices_into_hydro)
    where r_1e10 is in units of 10^10 cm and v_kms is in km/s, both ordered
    outer (largest R, largest V) to inner (smallest R, smallest V).
    """
    n_hydro = len(hydro.radius_cm)

    # Hydro data is stored inner->outer (radius_cm increasing). We need outer->inner.
    r_all = hydro.radius_cm[::-1] / 1e10  # outer->inner, 10^10 cm
    v_all = hydro.vel_cms[::-1] / 1e5  # outer->inner, km/s, preserve V sign (no abs)

    # Uniform spacing in velocity: place target_nd points equally spaced between
    # V_min and V_max. This ensures |V(I)-V(I+1)| = (V_max-V_min)/(ND-1) for all
    # adjacent zones, which is what MOM_JREL_V9 checks against VDOP_MIN.
    #
    # We do NOT use cumulative velocity distance (sum of |dV|) because for
    # non-monotonic velocity profiles, that would produce dV larger than
    # (V_max-V_min)/(ND-1), causing CMFGEN to insert extra points.
    # The non-monotonicity is handled by USE_LAM_ES=F (standard moment solver with SIGMA clamp), not
    # by the grid spacing.
    v_min = v_all.min()
    v_max = v_all.max()
    target_v = np.linspace(v_max, v_min, target_nd)  # outer (large V) -> inner (small V)

    # Interpolate R at the target velocities
    # np.interp requires increasing x, so sort V and R together
    sort_idx = np.argsort(v_all)
    v_sorted = v_all[sort_idx]
    r_sorted = r_all[sort_idx]
    r_downsampled = np.interp(target_v, v_sorted, r_sorted)

    # Ensure strict monotonicity (R decreasing outer->inner) and minimum spacing.
    # For non-monotonic velocity, R(V) may not be monotonic, so we force it.
    for i in range(len(r_downsampled) - 1):
        dr = r_downsampled[i] - r_downsampled[i + 1]
        if dr < MIN_DR_RDINR_1E10CM:
            r_downsampled[i + 1] = r_downsampled[i] - MIN_DR_RDINR_1E10CM

    # Use the target velocities directly as the V values.
    # We placed them uniformly in V space, so dV = (V_max - V_min)/(ND-1) is
    # guaranteed to be uniform and below VDOP_MIN. Re-interpolating V from R
    # would break this for non-monotonic velocity (V(R) is not single-valued).
    v_downsampled = target_v.copy()

    # Find nearest hydro indices for reference
    indices = np.array([np.argmin(np.abs(r_all - r)) for r in r_downsampled])

    return r_downsampled, v_downsampled, indices


def refine_rdinr_velocity_jumps(
    r_oi: np.ndarray, v_oi: np.ndarray, max_dv_kms: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Insert R points so no cell of the written grid carries |dV| > max_dv_kms.

    r_oi/v_oi are ordered outer->inner (r decreasing). Interior points are
    placed at equal-V steps within each violating cell, i.e. the R position is
    found by linear inverse interpolation of V(R) inside the cell. The
    Fortran side re-derives V(I) and SIGMA(I) from the SN_HYDRO_DATA spline at
    each R(I), so only the R placement controls the discretization; the
    written V column is diagnostic.

    Returns (r_new, v_new) outer->inner.
    """
    if max_dv_kms is None or max_dv_kms <= 0:
        return r_oi, v_oi
    r_out = [r_oi[0]]
    v_out = [v_oi[0]]
    n_added = 0
    for i in range(len(r_oi) - 1):
        dv = abs(v_oi[i + 1] - v_oi[i])
        n_sub = int(np.ceil(dv / max_dv_kms)) if max_dv_kms > 0 else 1
        for k in range(1, n_sub):
            f = k / n_sub
            r_out.append(r_oi[i] + (r_oi[i + 1] - r_oi[i]) * f)
            v_out.append(v_oi[i] + (v_oi[i + 1] - v_oi[i]) * f)
            n_added += 1
        r_out.append(r_oi[i + 1])
        v_out.append(v_oi[i + 1])
    if n_added:
        print(f"[RDINR refine] max_dV={max_dv_kms} km/s: inserted {n_added} "
              f"points ({len(r_oi)} -> {len(r_out)})")
    return np.asarray(r_out), np.asarray(v_out)


def refine_rdinr_boundaries(
    r_oi: np.ndarray, v_oi: np.ndarray,
    n_ib: int = 0, n_ob: int = 0, ratio: float = 2.5,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Geometric boundary refinement inside the outermost and innermost cells,
    mirroring SET_SN_R_GRID's N_OB_INS/N_IB_INS scheme (the boundary insertion
    used in Luc's SN models, dead under RD_IN_R_GRID=T).  The outermost cell
    is split into n_ob+1 sub-cells and the innermost into n_ib+1, with log-R
    widths forming a geometric series (finest sub-cell against the boundary).
    For n_ib=3, ratio=2.5 the inserted fractions are 0.039/0.137/0.381 of the
    cell's log-width -- close to SET_SN_R_GRID's hardcoded 0.06/0.16/0.4.

    Only the R placement controls the discretization: CMFGEN re-derives V and
    SIGMA from the SN_HYDRO_DATA spline at each radius, and V here is
    linear-interpolated inside the split cell purely as a diagnostic column.
    r_oi/v_oi are outer->inner (r decreasing).  Returns (r_new, v_new).
    """
    n_ib = int(n_ib or 0)
    n_ob = int(n_ob or 0)
    if len(r_oi) < 3 or (n_ib <= 0 and n_ob <= 0):
        return r_oi, v_oi

    def _split(r_out, r_in, v_out, v_in, n, boundary):
        # dlr is the log-width of the cell being split; w are the n+1
        # sub-cell widths as fractions of dlr, ascending away from the
        # boundary (finest at the boundary).
        dlr = np.log(r_out) - np.log(r_in)
        w = ratio ** np.arange(n + 1, dtype=float)
        w = w / w.sum()
        fracs = np.cumsum(w[:n])  # inserted points, fraction of dlr from boundary
        if boundary == 'inner':
            # Points sit above r_in in R; emit in outer->inner order.
            r_new = np.exp(np.log(r_in) + fracs * dlr)[::-1]
            f = 1.0 - fracs[::-1]  # interp fraction measured outer->inner
        else:
            r_new = np.exp(np.log(r_out) - fracs * dlr)
            f = fracs
        return r_new, v_out + (v_in - v_out) * f

    r_out = [r_oi[0]]
    v_out = [v_oi[0]]
    if n_ob > 0:
        ri, vi = _split(r_oi[0], r_oi[1], v_oi[0], v_oi[1], n_ob, 'outer')
        r_out += list(ri)
        v_out += list(vi)
    r_out += list(r_oi[1:-1])
    v_out += list(v_oi[1:-1])
    if n_ib > 0:
        ri, vi = _split(r_oi[-2], r_oi[-1], v_oi[-2], v_oi[-1], n_ib, 'inner')
        r_out += list(ri)
        v_out += list(vi)
    r_out.append(r_oi[-1])
    v_out.append(v_oi[-1])
    r_new = np.asarray(r_out)
    v_new = np.asarray(v_out)
    print(f"[RDINR boundary] inserted {n_ob} outer + {n_ib} inner boundary "
          f"points ({len(r_oi)} -> {len(r_new)}); outermost cell "
          f"{r_oi[0]:.4e}->{r_new[1]:.4e}, innermost "
          f"{r_new[-2]:.4e}->{r_oi[-1]:.4e} (1e10 cm)")
    return r_new, v_new


def write_rdinr(hydro: HydroData, path: Path, target_nd: Optional[int] = None,
                max_dv_kms: Optional[float] = None,
                n_ib: int = 0, n_ob: int = 0,
                bnd_ratio: float = 2.5) -> tuple[int, np.ndarray, np.ndarray]:
    """
    Write RDINR for RD_IN_R_GRID=T (same format as departure coefficient files).
    SET_RV_HYDRO_MODEL_V3 reads: '!Format date' line, header (R_ND, LUM, NOLD, ND),
    then for each depth: one line of 6 values (R,4 dummies,V), then READ(T1,J=1,NOLD).
    Some compilers advance to next record when NOLD=0, causing misalignment. Use NOLD=1
    and one dummy line per depth so the second read consumes a record. R in 10^10 cm,
    V in km/s; order outer->inner.

    If target_nd is specified, downsamples the hydro grid to target_nd points using
    velocity-aware spacing (see downsample_hydro_velocity_aware). This requires
    VEL_LAW=11 with RD_IN_R_GRID=T so CMFGEN uses the RDINR grid directly.

    Returns (nd_written, r_oi, v_oi) where r_oi/v_oi are the arrays actually written
    (outer->inner order, r in 10^10 cm, v in km/s).
    """
    nd_hydro = int(hydro.params["ND"])
    if len(hydro.radius_cm) < nd_hydro or len(hydro.vel_cms) < nd_hydro:
        raise ValueError(
            f"RDINR: need at least {nd_hydro} points, got radius_cm={len(hydro.radius_cm)}, vel_cms={len(hydro.vel_cms)}"
        )

    if target_nd is not None and int(target_nd) != nd_hydro:
        # Only subset (never interpolate in V) when the hydro file is larger
        # than the requested ND.  Uniform-dV interpolation destroys
        # non-monotonic velocity profiles — R(V) is not single-valued at a
        # reverse shock — so we keep original hydro (R,V) pairs.
        r_oi, v_oi, _ = downsample_hydro_gradient_aware(hydro, int(target_nd))
        nd = len(r_oi)
        print(f"[RDINR] Subset SN_HYDRO_DATA {nd_hydro} -> {nd} points "
              f"(original STELLA/hydro values, no V-interpolation)")
    else:
        r_oi = (hydro.radius_cm[::-1] / 1e10).astype(float, copy=True)
        v_oi = hydro.vel_cms[::-1] / 1e5  # preserve V sign (no abs)
        nd = nd_hydro
        if target_nd is not None:
            print(f"[RDINR] Using SN_HYDRO_DATA grid as-is (ND={nd}; "
                  f"same mapping as the converter)")

    # Refine cells whose per-cell Doppler shift exceeds max_dv_kms.  The base
    # grid carries the shock compression (sigma ~ -120) as a single ~1650 km/s
    # cell, which is where the lagged CMF coupling goes non-contractive.
    r_oi, v_oi = refine_rdinr_velocity_jumps(r_oi, v_oi, max_dv_kms)

    # Geometric boundary refinement: Luc's N_IB_INS/N_OB_INS scheme from
    # SET_SN_R_GRID, applied RDINR-side because RD_IN_R_GRID=T bypasses the
    # Fortran boundary insertion.  Extra inner points resolve the steep tau
    # rise feeding the diffusion inner boundary (where the stiff depth-ND
    # ionization residuals live); extra outer points resolve the
    # free-streaming transition in the thin shell edge.
    r_oi, v_oi = refine_rdinr_boundaries(r_oi, v_oi, n_ib=n_ib, n_ob=n_ob,
                                         ratio=bnd_ratio)
    nd = len(r_oi)

    # Enforce minimum spacing (outer->inner: r_oi[i] > r_oi[i+1]); avoids velocity_law.f SIGFPE
    for i in range(nd - 1):
        dr = r_oi[i] - r_oi[i + 1]
        if dr < MIN_DR_RDINR_1E10CM:
            r_oi[i + 1] = r_oi[i] - MIN_DR_RDINR_1E10CM

    with open(path, "w", newline="\n") as f:
        f.write(" 24-FEB-2004                    !Format date\n")
        # NOLD=1 so code does READ(T1,J=1,1) and consumes one value per depth
        f.write(f" {float(r_oi[-1]):22.15E}    {1.0:11.4E}     {1:4d}    {nd:4d}\n")
        for i in range(nd):
            f.write(
                f" {float(r_oi[i]):22.15E} {0.0:22.15E} {0.0:22.15E} {0.0:22.15E} {0.0:22.15E} {float(v_oi[i]):22.15E}\n"
            )
            # One dummy value so READ(T1,J=1,NOLD) with NOLD=1 consumes a record
            f.write(" 0.0\n")

    return nd, r_oi, v_oi


def write_input_hydro_dat(hydro: HydroData, path: Path, use_hydro_grid: bool = False) -> None:
    """
    Write input_hydro.dat (R, rho, T per depth) for workflows that use it.
    Format: USE_HYDRO_GRID ADD_ACC_ZONE, header line, NHYDRO, bogus line, then NHYDRO lines of
    R(10^10 cm), dummy, rho(g/cm^3), T(K) in order outer->inner.
    Units: T is in Kelvin (not 10^4 K). HydroData stores T in 10^4 K; we convert when writing.
    Note: With VEL_LAW=12 and RD_IN_R_GRID=T, CMFGEN gets T/ED from SN_HYDRO_DATA (10^4 K)
    and does not read this file for T; this file is still written for consistency and other use.
    """
    r_cm = hydro.radius_cm
    rho = hydro.rho_cgs
    n = len(r_cm)
    if n == 0 or len(rho) < n:
        raise ValueError("write_input_hydro_dat: need radius and density arrays")
    # Temperature: HydroData has T in 10^4 K; file format is T in Kelvin (convert * 1e4). Fallback 5000 K.
    if hydro.temp_104k is not None and len(hydro.temp_104k) >= n:
        t_k = hydro.temp_104k[:n] * 1e4
    else:
        t_k = np.full(n, 3000.0)
    r_1e10 = r_cm / 1e10
    with open(path, "w", newline="\n") as f:
        f.write(f" {'T' if use_hydro_grid else 'F'} F\n")  # USE_HYDRO_GRID ADD_ACC_ZONE
        f.write("! R(1e10 cm) dummy rho(g/cm^3) T(K) from SN_HYDRO_DATA\n")
        f.write(f" {n}\n")
        f.write("! NHYDRO lines below (outer -> inner)\n")
        # Outer -> inner so set_abund_clump stores first line at R_HYDRO(NHYDRO), last at R_HYDRO(1)
        for i in range(n - 1, -1, -1):
            f.write(f" {r_1e10[i]:22.15E} {0.0:22.15E} {rho[i]:22.15E} {t_k[i]:22.15E}\n")


# Electron-scattering opacity (cm^2/g) fallback when Kappa section missing in SN_HYDRO_DATA
KAPPA_ES_CGS = 0.34


def smooth_velocity_jumps(
    r_cm: np.ndarray, v_kms: np.ndarray, rho: np.ndarray, t_k: np.ndarray,
    *, max_dv_ratio: float = 2.0, n_insert: int = 8
) -> tuple:
    """
    Insert intermediate zones at velocity discontinuities to smooth the jump.

    The interacting-SN hydro grid can have a 35x velocity jump at the shell-ejecta
    interface. CMFGEN's Lambda iteration cannot converge with such a discontinuity
    because the Sobolev sigma (dV/dR - V/R) becomes enormous. This function inserts
    logarithmically-spaced intermediate zones that transition velocity, density, and
    temperature smoothly across the jump.

    Parameters
    ----------
    r_cm, v_kms, rho, t_k : arrays in inner-to-outer order (CMFGEN convention,
        index 0 = innermost). The file is written outer-to-inner, but internally
        HydroData stores inner-to-outer.
    max_dv_ratio : if V_max/V_min between consecutive zones exceeds this, insert zones.
    n_insert : number of intermediate zones to insert at each detected jump.

    Returns
    -------
    (r_new, v_new, rho_new, t_new, n_jumps) : smoothed arrays and count of jumps found.
    """
    n = len(r_cm)
    if n < 2:
        return r_cm, v_kms, rho, t_k, 0

    # Detect jumps: consecutive zones where velocity ratio exceeds threshold
    jump_indices = []
    for i in range(n - 1):
        v_lo = min(abs(v_kms[i]), abs(v_kms[i + 1]))
        v_hi = max(abs(v_kms[i]), abs(v_kms[i + 1]))
        if v_lo > 0 and v_hi / v_lo > max_dv_ratio:
            jump_indices.append(i)

    if not jump_indices:
        return r_cm, v_kms, rho, t_k, 0

    print(f"[smooth_velocity_jumps] Found {len(jump_indices)} jump(s) with ratio > {max_dv_ratio}")
    for ji in jump_indices:
        v_lo = min(abs(v_kms[ji]), abs(v_kms[ji + 1]))
        v_hi = max(abs(v_kms[ji]), abs(v_kms[ji + 1]))
        print(f"  idx {ji}->{ji+1}: V={v_kms[ji]:.1f} -> {v_kms[ji+1]:.1f} km/s "
              f"(ratio={v_hi/v_lo:.1f}), R={r_cm[ji]:.4e} -> {r_cm[ji+1]:.4e} cm, "
              f"rho={rho[ji]:.4e} -> {rho[ji+1]:.4e}")

    # Build new arrays with inserted zones
    # Work in inner-to-outer order (index 0 = innermost)
    r_parts = []
    v_parts = []
    rho_parts = []
    t_parts = []
    prev = 0
    for ji in jump_indices:
        # Add zones up to (but not including) the jump point
        r_parts.append(r_cm[prev:ji + 1])
        v_parts.append(v_kms[prev:ji + 1])
        rho_parts.append(rho[prev:ji + 1])
        t_parts.append(t_k[prev:ji + 1])

        # Insert n_insert zones between ji and ji+1 using log interpolation
        # Log interpolation is physical: V ~ R^-1 for homologous, density ~ R^-3
        r0, r1 = r_cm[ji], r_cm[ji + 1]
        v0, v1 = v_kms[ji], v_kms[ji + 1]
        rho0, rho1 = rho[ji], rho[ji + 1]
        t0, t1 = t_k[ji], t_k[ji + 1]

        # Log-spaced in radius (geometric progression)
        log_frac = np.linspace(0, 1, n_insert + 2)[1:-1]  # exclude endpoints
        r_mid = r0 * (r1 / r0) ** log_frac
        # Log interpolation for V, rho, T (handle sign for V)
        v_mid = np.sign(v0) * np.abs(v0) * (np.abs(v1) / np.abs(v0)) ** log_frac
        rho_mid = rho0 * (rho1 / rho0) ** log_frac
        t_mid = t0 * (t1 / t0) ** log_frac

        r_parts.append(r_mid)
        v_parts.append(v_mid)
        rho_parts.append(rho_mid)
        t_parts.append(t_mid)

        prev = ji + 1

    # Add remaining zones
    r_parts.append(r_cm[prev:])
    v_parts.append(v_kms[prev:])
    rho_parts.append(rho[prev:])
    t_parts.append(t_k[prev:])

    r_new = np.concatenate(r_parts)
    v_new = np.concatenate(v_parts)
    rho_new = np.concatenate(rho_parts)
    t_new = np.concatenate(t_parts)

    print(f"[smooth_velocity_jumps] Inserted {len(jump_indices) * n_insert} zones: "
          f"{n} -> {len(r_new)} total")

    return r_new, v_new, rho_new, t_new, len(jump_indices)


# Solar composition mass fractions (Asplund et al. 2009, Grevesse et al. 2010)
SOLAR_MASS_FRACTIONS = {
    "HYD": 0.7381,
    "HE": 0.2485,
    "CARB": 0.0030,
    "OXY": 0.0096,
    "SIL": 0.0066,
    "CAL": 0.00064,
    "IRON": 0.0126,
    "NICK": 0.00007,
    "COB": 0.00002,
}


def blend_composition_to_solar(
    hydro: HydroData, sigma: float = 5.0, solar_weight: float = 0.0,
    min_mass_fraction: float = 0.0,
) -> HydroData:
    """
    Blend the chemical composition across discontinuities by Gaussian-smoothing
    the mass fractions in log-space. This removes the 100% chemical jumps
    (e.g. H: 0->0.7, Fe: 0.44->0.0) that prevent Lambda-iteration convergence.

    Parameters:
        sigma: Gaussian smoothing width in zones (applied to log10 of mass fractions)
        solar_weight: If > 0, blend toward solar composition with this weight (0-1).
                      0 = pure Gaussian smoothing of existing profile,
                      1 = replace with solar composition entirely.
        min_mass_fraction: If > 0, set a floor on every species at every depth.
                           Species with mass fractions below this value are raised
                           to this value, then all mass fractions are renormalized
                           to sum to 1. This ensures that every species has enough
                           atoms/levels to affect the J field when populations change.
                           Default: 0 (no floor).

    Returns a new HydroData with blended mass fractions. All other arrays unchanged.
    """
    from scipy.ndimage import gaussian_filter1d

    n = len(hydro.radius_cm)
    new_mf = {}
    all_keys = set(hydro.mass_fractions.keys()) | set(SOLAR_MASS_FRACTIONS.keys())

    for key in all_keys:
        arr = hydro.mass_fractions.get(key, np.full(n, 1e-30))
        solar = SOLAR_MASS_FRACTIONS.get(key, 1e-30)

        if solar_weight > 0:
            # Blend toward solar: new = (1-w)*original + w*solar
            blended = arr * (1 - solar_weight) + solar * solar_weight
        else:
            # Pure Gaussian smoothing of the existing profile
            # Work in log-space to handle the many orders of magnitude
            log_arr = np.log10(np.maximum(arr, 1e-30))
            log_smoothed = gaussian_filter1d(log_arr, sigma=sigma)
            blended = 10.0 ** log_smoothed

        # Ensure non-negative
        blended = np.maximum(blended, 1e-30)
        new_mf[key] = blended

    # Apply minimum mass fraction floor
    if min_mass_fraction > 0:
        # First renormalize, then floor, then renormalize again
        total = np.zeros(n)
        for key in new_mf:
            total += new_mf[key]
        for key in new_mf:
            new_mf[key] = new_mf[key] / total
        # Now apply floor
        for key in new_mf:
            below = new_mf[key] < min_mass_fraction
            n_below = int(np.sum(below))
            if n_below > 0:
                new_mf[key] = np.where(below, min_mass_fraction, new_mf[key])
        print(f"[blend_composition] Applied min mass fraction floor: {min_mass_fraction:.1e}")
        # Renormalize again after flooring
        total = np.zeros(n)
        for key in new_mf:
            total += new_mf[key]
        for key in new_mf:
            new_mf[key] = new_mf[key] / total
        # Final check: ensure floor is maintained (iterative)
        for key in new_mf:
            below = new_mf[key] < min_mass_fraction
            if np.any(below):
                new_mf[key] = np.maximum(new_mf[key], min_mass_fraction)
        total = np.zeros(n)
        for key in new_mf:
            total += new_mf[key]
        for key in new_mf:
            new_mf[key] = new_mf[key] / total

    # Renormalize each zone so mass fractions sum to 1
    total = np.zeros(n)
    for key in new_mf:
        total += new_mf[key]
    for key in new_mf:
        new_mf[key] = new_mf[key] / total

    new_hydro = HydroData(
        params=hydro.params,
        mass_fractions=new_mf,
        radius_cm=hydro.radius_cm.copy(),
        rho_cgs=hydro.rho_cgs.copy(),
        vel_cms=hydro.vel_cms.copy(),
        temp_104k=hydro.temp_104k.copy() if hydro.temp_104k is not None else None,
        kappa_cgs=hydro.kappa_cgs.copy() if hydro.kappa_cgs is not None else None,
        ed_cgs=hydro.ed_cgs.copy() if hydro.ed_cgs is not None else None,
        na_cgs=hydro.na_cgs.copy() if hydro.na_cgs is not None else None,
    )

    # Report the blending
    print(f"[blend_composition] Blended composition: sigma={sigma}, solar_weight={solar_weight}")
    print(f"[blend_composition] Mass fraction ranges after blending:")
    for key in ["HYD", "HE", "OXY", "SIL", "IRON", "CARB", "CAL"]:
        if key in new_mf:
            print(f"  {key}: {new_mf[key].min():.4e} to {new_mf[key].max():.4e}")

    return new_hydro


def monotonize_velocity(hydro: HydroData) -> HydroData:
    """
    Replace non-monotonic velocity regions with a smooth monotonic interpolation.

    In CMFGEN's convention (outer to inner), velocity should decrease monotonically.
    Where the velocity reverses (goes from low back to high), replace it with a
    log-linear interpolation that preserves monotonicity. This removes negative-sigma
    (dv/dr < 0) regions that can cause convergence problems.

    All other arrays (density, temperature, etc.) are preserved — only velocity is modified.
    """
    n = len(hydro.radius_cm)
    v = hydro.vel_cms.copy()  # inner-to-outer order

    # In inner-to-outer order, velocity should INCREASE monotonically for a homologous
    # ejecta. The non-monotonicity is where V decreases going outward.
    # Find the reversal region: where v[i+1] < v[i] (decrease going outward)
    reversals = []
    in_reversal = False
    start = 0
    for i in range(n - 1):
        if v[i + 1] < v[i] - 1.0:  # significant decrease
            if not in_reversal:
                start = i
                in_reversal = True
        else:
            if in_reversal:
                reversals.append((start, i))
                in_reversal = False
    if in_reversal:
        reversals.append((start, n - 1))

    if not reversals:
        print("[monotonize_velocity] No velocity reversals found — already monotonic")
        return hydro

    for (s, e) in reversals:
        # Interpolate from v[s] to v[e+1] (or v[e] if at boundary) in log-space
        end = min(e + 1, n - 1)
        if end <= s:
            continue
        v0, v1 = v[s], v[end]
        if v0 <= 0:
            v0 = 1.0
        if v1 <= 0:
            v1 = 1.0
        frac = np.linspace(0, 1, end - s + 1)
        v[s:end + 1] = v0 * (v1 / v0) ** frac
        print(f"[monotonize_velocity] Replaced reversal at idx {s}-{end}: "
              f"V was {v0/1e5:.1f}->{v1/1e5:.1f} km/s (now smooth monotonic)")

    # Second pass: fix any remaining non-monotonic points by capping
    for i in range(n - 1):
        if v[i + 1] < v[i]:
            v[i + 1] = v[i] * 1.0001  # ensure strictly increasing outward

    new_hydro = HydroData(
        params=hydro.params,
        mass_fractions=hydro.mass_fractions,
        radius_cm=hydro.radius_cm.copy(),
        rho_cgs=hydro.rho_cgs.copy(),
        vel_cms=v,
        temp_104k=hydro.temp_104k.copy() if hydro.temp_104k is not None else None,
        kappa_cgs=hydro.kappa_cgs.copy() if hydro.kappa_cgs is not None else None,
        ed_cgs=hydro.ed_cgs.copy() if hydro.ed_cgs is not None else None,
        na_cgs=hydro.na_cgs.copy() if hydro.na_cgs is not None else None,
    )
    return new_hydro


def insert_zones_at_velocity_reversals(
    hydro: HydroData, *, n_insert: int = 50, min_dv_km_s: float = 100.0
) -> HydroData:
    """
    Insert log-interpolated zones at velocity reversal regions (negative sigma).

    In a non-monotonic velocity profile, the velocity decreases going outward at
    certain points (dv/dr < 0). These regions need extra grid resolution to
    properly resolve the radiative transfer. This function inserts n_insert zones
    at each velocity reversal, log-interpolating all arrays.

    Unlike smooth_hydrodata_jumps, this PRESERVES the non-monotonic velocity
    profile — it just adds more grid points across the reversal.

    Parameters:
        n_insert: Number of zones to insert at each reversal.
        min_dv_km_s: Minimum velocity decrease (km/s) to count as a reversal.
                     Default: 100 km/s (ignores small numerical glitches).
    """
    n_old = len(hydro.radius_cm)
    v = hydro.vel_cms.copy()  # inner-to-outer

    # Find velocity reversals (where v decreases going outward by > min_dv)
    jump_indices = []
    for i in range(n_old - 1):
        if (v[i] - v[i + 1]) > min_dv_km_s * 1e5:  # convert km/s to cm/s
            jump_indices.append(i)

    if not jump_indices:
        print(f"[insert_reversals] No velocity reversals found")
        return hydro

    print(f"[insert_reversals] Found {len(jump_indices)} velocity reversal(s), "
          f"inserting {n_insert} zones each")

    # Build new arrays with interpolated zones at reversals
    def build_new(arr, is_log=True):
        if arr is None:
            return None
        parts = []
        prev = 0
        for ji in jump_indices:
            parts.append(arr[prev:ji + 1])
            a0, a1 = arr[ji], arr[ji + 1]
            frac = np.linspace(0, 1, n_insert + 2)[1:-1]
            if is_log and a0 > 0 and a1 > 0:
                mid = a0 * (a1 / a0) ** frac
            else:
                mid = a0 * (1 - frac) + a1 * frac
            parts.append(mid)
            prev = ji + 1
        parts.append(arr[prev:])
        return np.concatenate(parts)

    new_mf = {}
    for key in hydro.mass_fractions:
        new_mf[key] = build_new(hydro.mass_fractions[key], is_log=False)

    # Renormalize
    n_new = len(new_mf[list(new_mf.keys())[0]])
    total = np.zeros(n_new)
    for key in new_mf:
        total += new_mf[key]
    for key in new_mf:
        new_mf[key] = new_mf[key] / total

    new_r = build_new(hydro.radius_cm, is_log=True)
    new_v = build_new(hydro.vel_cms, is_log=True)
    new_rho = build_new(hydro.rho_cgs, is_log=True)
    new_temp = build_new(hydro.temp_104k, is_log=True) if hydro.temp_104k is not None else None
    new_kappa = build_new(hydro.kappa_cgs, is_log=True) if hydro.kappa_cgs is not None else None
    new_ed = build_new(hydro.ed_cgs, is_log=True) if hydro.ed_cgs is not None else None
    new_na = build_new(hydro.na_cgs, is_log=True) if hydro.na_cgs is not None else None

    n_new = len(new_r)
    new_hydro = HydroData(
        params={**hydro.params, "ND": float(n_new)},
        mass_fractions=new_mf,
        radius_cm=new_r,
        rho_cgs=new_rho,
        vel_cms=new_v,
        temp_104k=new_temp,
        kappa_cgs=new_kappa,
        ed_cgs=new_ed,
        na_cgs=new_na,
    )

    print(f"[insert_reversals] Grid expanded: {n_old} -> {n_new} zones")
    return new_hydro


def insert_zones_at_chemical_interfaces(
    hydro: HydroData, *, n_insert: int = 20, jump_threshold: float = 0.1
) -> HydroData:
    """
    Insert log-interpolated zones at chemical composition discontinuities.

    Detects zones where mass fractions change by more than jump_threshold (summed
    across all species) and inserts n_insert intermediate zones with log-interpolated
    composition, density, temperature, velocity, and radius.

    This ensures adequate grid resolution where opacity changes rapidly due to
    chemical transitions (e.g., H/He shell → O/Si/Fe ejecta interface).
    """
    n_old = len(hydro.radius_cm)
    mf = hydro.mass_fractions

    # Detect chemical jumps (inner-to-outer order)
    jump_indices = []
    for i in range(n_old - 1):
        total_jump = 0.0
        for key in mf:
            total_jump += abs(mf[key][i] - mf[key][i + 1])
        if total_jump > jump_threshold:
            jump_indices.append(i)

    if not jump_indices:
        print(f"[insert_chemical] No chemical jumps > {jump_threshold} found")
        return hydro

    print(f"[insert_chemical] Found {len(jump_indices)} chemical jump(s), "
          f"inserting {n_insert} zones each")

    # Build new arrays with interpolated zones at chemical jumps
    def build_new(arr, is_log=True):
        if arr is None:
            return None
        parts = []
        prev = 0
        for ji in jump_indices:
            parts.append(arr[prev:ji + 1])
            a0, a1 = arr[ji], arr[ji + 1]
            frac = np.linspace(0, 1, n_insert + 2)[1:-1]
            if is_log and a0 > 0 and a1 > 0:
                mid = a0 * (a1 / a0) ** frac
            else:
                mid = a0 * (1 - frac) + a1 * frac
            parts.append(mid)
            prev = ji + 1
        parts.append(arr[prev:])
        return np.concatenate(parts)

    # Build new mass fractions with log-interpolated values
    new_mf = {}
    for key in mf:
        new_mf[key] = build_new(mf[key], is_log=False)  # linear for mass fractions

    # Renormalize mass fractions to sum=1 at each zone
    n_new = len(new_mf[list(new_mf.keys())[0]])
    total = np.zeros(n_new)
    for key in new_mf:
        total += new_mf[key]
    for key in new_mf:
        new_mf[key] = new_mf[key] / total

    new_r = build_new(hydro.radius_cm, is_log=True)
    new_v = build_new(hydro.vel_cms, is_log=True)
    new_rho = build_new(hydro.rho_cgs, is_log=True)
    new_temp = build_new(hydro.temp_104k, is_log=True) if hydro.temp_104k is not None else None
    new_kappa = build_new(hydro.kappa_cgs, is_log=True) if hydro.kappa_cgs is not None else None
    new_ed = build_new(hydro.ed_cgs, is_log=True) if hydro.ed_cgs is not None else None
    new_na = build_new(hydro.na_cgs, is_log=True) if hydro.na_cgs is not None else None

    n_new = len(new_r)
    new_hydro = HydroData(
        params={**hydro.params, "ND": float(n_new)},
        mass_fractions=new_mf,
        radius_cm=new_r,
        rho_cgs=new_rho,
        vel_cms=new_v,
        temp_104k=new_temp,
        kappa_cgs=new_kappa,
        ed_cgs=new_ed,
        na_cgs=new_na,
    )

    print(f"[insert_chemical] Grid expanded: {n_old} -> {n_new} zones")
    return new_hydro


def smooth_hydrodata_jumps(
    hydro: HydroData, *, max_dv_ratio: float = 2.0, n_insert: int = 8,
    smooth_rho_sigma: float = 0.0, smooth_temp_sigma: float = 0.0,
) -> HydroData:
    """
    Smooth velocity discontinuities in a HydroData object by inserting
    log-interpolated intermediate zones. All arrays (radius, velocity, density,
    temperature, kappa, ed, na, mass_fractions) are updated consistently.

    If smooth_rho_sigma > 0, also apply Gaussian smoothing to the density
    profile in log-space across the interface to reduce the 400x density jump.
    If smooth_temp_sigma > 0, apply Gaussian smoothing to temperature.

    This must be called BEFORE any file writing (RDINR, T_IN, GAMMAS_IN,
    FIN_CAL_GRID, input_hydro.dat) so all files share the same smoothed grid.
    """
    n_old = len(hydro.radius_cm)
    v_kms = hydro.vel_cms / 1e5

    # Detect jump indices (same logic as smooth_velocity_jumps)
    jump_indices = []
    for i in range(n_old - 1):
        v_lo = min(abs(v_kms[i]), abs(v_kms[i + 1]))
        v_hi = max(abs(v_kms[i]), abs(v_kms[i + 1]))
        if v_lo > 0 and v_hi / v_lo > max_dv_ratio:
            jump_indices.append(i)

    if not jump_indices:
        return hydro

    print(f"[smooth_hydrodata_jumps] Found {len(jump_indices)} jump(s), inserting {n_insert} zones each")

    # Build new arrays by iterating through old data and inserting zones at jumps
    def build_new(arr, is_log=True):
        """Build new array with interpolated zones at jumps."""
        if arr is None:
            return None
        parts = []
        prev = 0
        for ji in jump_indices:
            parts.append(arr[prev:ji + 1])
            a0, a1 = arr[ji], arr[ji + 1]
            frac = np.linspace(0, 1, n_insert + 2)[1:-1]
            if is_log and a0 > 0 and a1 > 0:
                mid = a0 * (a1 / a0) ** frac
            else:
                mid = a0 * (1 - frac) + a1 * frac
            parts.append(mid)
            prev = ji + 1
        parts.append(arr[prev:])
        return np.concatenate(parts)

    # Build all new arrays
    new_r = build_new(hydro.radius_cm, is_log=True)
    new_v = build_new(hydro.vel_cms, is_log=True)
    new_rho = build_new(hydro.rho_cgs, is_log=True)
    new_temp = build_new(hydro.temp_104k, is_log=True) if hydro.temp_104k is not None else None
    new_kappa = build_new(hydro.kappa_cgs, is_log=True) if hydro.kappa_cgs is not None else None
    new_ed = build_new(hydro.ed_cgs, is_log=True) if hydro.ed_cgs is not None else None
    new_na = build_new(hydro.na_cgs, is_log=True) if hydro.na_cgs is not None else None

    new_mf = {}
    for key, arr in hydro.mass_fractions.items():
        new_mf[key] = build_new(arr, is_log=False)  # linear for mass fractions (can be zero)

    n_new = len(new_r)
    # Optional Gaussian smoothing of density and temperature in log-space
    # This reduces the 400x density jump at the shell interface that prevents
    # Lambda-iteration convergence even after velocity smoothing.
    if smooth_rho_sigma > 0 and new_rho is not None:
        from scipy.ndimage import gaussian_filter1d
        log_rho = np.log10(np.maximum(new_rho, 1e-30))
        log_rho_smoothed = gaussian_filter1d(log_rho, sigma=smooth_rho_sigma)
        new_rho = 10.0 ** log_rho_smoothed
        # Preserve total mass: scale so integral of rho * r^2 dr is conserved
        # (approximate — just rescale to match original total)
        old_mass_integral = np.sum(new_rho * new_r**2 * np.gradient(new_r))
        new_mass_integral = np.sum(new_rho * new_r**2 * np.gradient(new_r))
        if new_mass_integral > 0:
            new_rho *= old_mass_integral / new_mass_integral
        print(f"[smooth_hydrodata_jumps] Gaussian density smoothing: sigma={smooth_rho_sigma}")

    if smooth_temp_sigma > 0 and new_temp is not None:
        from scipy.ndimage import gaussian_filter1d
        log_t = np.log10(np.maximum(new_temp, 1e-3))
        log_t_smoothed = gaussian_filter1d(log_t, sigma=smooth_temp_sigma)
        new_temp = 10.0 ** log_t_smoothed
        print(f"[smooth_hydrodata_jumps] Gaussian temperature smoothing: sigma={smooth_temp_sigma}")

    new_hydro = HydroData(
        params={**hydro.params, "ND": float(n_new)},
        mass_fractions=new_mf,
        radius_cm=new_r,
        rho_cgs=new_rho,
        vel_cms=new_v,
        temp_104k=new_temp,
        kappa_cgs=new_kappa,
        ed_cgs=new_ed,
        na_cgs=new_na,
    )

    print(f"[smooth_hydrodata_jumps] HydroData smoothed: {n_old} -> {n_new} zones")
    return new_hydro


def write_input_hydro_dat_sniin(
    hydro: HydroData, path: Path, use_hydro_grid: bool = False,
    *, smooth_jumps: bool = False, smooth_max_dv_ratio: float = 2.0,
    smooth_n_insert: int = 8
) -> None:
    """
    Write input_hydro.dat in the format read by RV_SN_MODEL_SNIIN (VEL_LAW=12) and
    SET_ABUND_CLUMP (cols 1-4 per line: R(1e10 cm), dummy, RHO(g/cm^3), T(K)).
    Format: USE_HYDRO_GRID ADD_ACC_ZONE, header, NHYDRO, bogus line, then per point (outer->inner):
    R, V, RHO, T(K), XL, XTAU, TAU.

    CRITICAL: The TAU column must use the same constant kappa (0.34 cm^2/g) that the
    STELLA-to-CMFGEN converter uses for inner boundary selection. The converter selects
    zones where tau(kappa=0.34) reaches cmfgen_tau_max, but the actual STELLA NLTE
    Rosseland kappa can be 60x lower in cold neutral gas. If input_hydro.dat uses the
    NLTE kappa, RV_SN_MODEL_SNIIN sees tau=1.1 instead of tau=67, constructs a grid
    far too shallow for the diffusion boundary condition, and CMFGEN halts with
    "Error in MON_INTERP - values outside range" when HYDRO_TERMS_V5 tries to find
    the photosphere at tau=2/3 outside the grid.
    """
    r_cm = hydro.radius_cm
    rho = hydro.rho_cgs
    n = len(r_cm)
    if n == 0 or len(rho) < n:
        raise ValueError("write_input_hydro_dat_sniin: need radius and density")
    if hydro.temp_104k is None or len(hydro.temp_104k) < n:
        raise ValueError("write_input_hydro_dat_sniin: SN_HYDRO_DATA must have Temperature section for VEL_LAW=12")
    r_1e10 = r_cm / 1e10
    v_kms = hydro.vel_cms / 1e5  # cm/s -> km/s, preserve V sign (no abs) for infall support
    t_k = hydro.temp_104k[:n] * 1e4  # 10^4 K -> K
    # Use constant electron-scattering kappa (0.34 cm^2/g) for the tau column, NOT the
    # STELLA NLTE kappa. This must match the kappa used by the converter for grid selection.
    # The NLTE kappa (hydro.kappa_cgs) can be 60x lower (e.g. 0.003 vs 0.34) in cold
    # neutral gas, producing tau=1.1 instead of tau=67. RV_SN_MODEL_SNIIN uses this tau
    # to construct the R grid; a too-low tau makes the grid too shallow for the diffusion
    # boundary condition (requires tau>=50).
    kappa = np.full(n, KAPPA_ES_CGS)
    # Tau from outer to inner: tau[n-1]=0, tau[j] = tau[j+1] + 0.5*(k*rho)[j,j+1] * (r[j+1]-r[j])
    # The outer tau must be > 0 because RV_SN_MODEL_SNIIN takes LOG(TAU_HYDRO(1)) to build a
    # geometric tau grid; a zero outer tau produces NaN/Inf and SIGSEGV in LIN_INTERP.
    # CRITICAL: The outer tau must NOT be too small relative to the second point's tau.
    # RV_SN_MODEL_SNIIN builds a geometric grid from TAU(outer) to TAU(inner). If the outer
    # tau is 1e-30 but the second point is ~1e-4, the geometric grid puts ~80% of its points
    # in the empty gap [1e-30, 1e-4] where LIN_INTERP has no hydro data, producing NaN R
    # values and SIGFPE in the SIGMA computation (division by V=0 from extrapolated R).
    # Fix: set the outer tau to a fraction of the second point's tau so the geometric grid
    # spans a reasonable dynamic range.
    tau_oi = np.zeros(n)
    for j in range(n - 2, -1, -1):
        dr = r_cm[j + 1] - r_cm[j]  # > 0 (outer to inner)
        krho_avg = 0.5 * (kappa[j] * rho[j] + kappa[j + 1] * rho[j + 1])
        tau_oi[j] = tau_oi[j + 1] + krho_avg * dr
    if n >= 2 and tau_oi[-2] > 0:
        # Set outer tau to 1% of the second point's tau (geometric grid needs >0 start).
        # This keeps the tau dynamic range manageable for RV_SN_MODEL_SNIIN's geometric grid.
        tau_oi[-1] = max(tau_oi[-1], tau_oi[-2] * 1.0e-2)
    else:
        tau_oi[-1] = max(tau_oi[-1], 1.0e-30)
    # Print so we can verify the data being written matches what was read from SN_HYDRO_DATA
    print(f"[input_hydro.dat] NHYDRO={n} (from SN_HYDRO_DATA)")
    print(f"[input_hydro.dat] R(1e10 cm): min={r_1e10.min():.6e} max={r_1e10.max():.6e}")
    print(f"[input_hydro.dat] RHO(g/cm^3): min={rho.min():.6e} max={rho.max():.6e}")
    print(f"[input_hydro.dat] T(K): min={t_k.min():.2f} max={t_k.max():.2f}")
    with open(path, "w", newline="\n") as f:
        f.write(f" {'T' if use_hydro_grid else 'F'} F\n")  # USE_HYDRO_GRID ADD_ACC_ZONE
        f.write("R(1e10cm) V(km/s) RHO T(K) XL XTAU TAU\n")
        f.write(f" {n}\n")
        f.write("NHYDRO lines below (outer -> inner)\n")
        for i in range(n - 1, -1, -1):
            xl, xtau = 0.0, 0.0
            tau_val = tau_oi[i]
            f.write(
                f" {r_1e10[i]:22.15E} {v_kms[i]:22.15E} {rho[i]:22.15E} {t_k[i]:22.15E} {xl:22.15E} {xtau:22.15E} {tau_val:22.15E}\n"
            )
    print(f"Wrote {path.name} (VEL_LAW=12 SNIIN format): ND={n}  R(1e10 cm) outer->inner, TAU outer={tau_oi[-1]:.4e} inner={tau_oi[0]:.4e}")


def write_t_in_from_hydro(hydro: HydroData, path: Path, *, nold: int = 0, t_max_104: float = 1.0, t_min_104: float = 0.1) -> None:
    """
    Write T_IN from HydroData (format read by REGRID_T_ED / INIT_TEMP_V2).
    Format matches star_init/NT35000_logg425: leading blank, '!Format date' line,
    header (R_core, LUMOLD, NOLD, ND), then per depth: one line of 8 columns
    (R, dummy, ED, T, ION_FRAC, VEL, CLUMP, depth_index) in 12.7E style, then
    NOLD values (5 per line) or single " 0.0" when NOLD=0. Order outer->inner.

    Args:
        t_max_104: Maximum temperature in units of 10^4 K (default: 1.0 = 10,000 K)
        t_min_104: Minimum temperature in units of 10^4 K (default: 0.1 = 1000 K).
            Must match the converter's t_floor_K so that proxy departure
            coefficients (computed relative to LTE at the STELLA temperature)
            are applied to LTE at the same temperature in CMFGEN.
    """
    n = len(hydro.radius_cm)
    if n == 0:
        raise ValueError("write_t_in_from_hydro: empty hydro")
    if hydro.temp_104k is None or len(hydro.temp_104k) < n:
        raise ValueError("write_t_in_from_hydro: SN_HYDRO_DATA must have Temperature section")
    if hydro.ed_cgs is None or len(hydro.ed_cgs) < n:
        raise ValueError("write_t_in_from_hydro: SN_HYDRO_DATA must have Electron density section")
    r_1e10 = hydro.radius_cm / 1e10
    v_kms = hydro.vel_cms / 1e5  # preserve V sign (no abs)
    t_104 = hydro.temp_104k[:n]
    ed = hydro.ed_cgs[:n]
    r_core = float(r_1e10[0])
    lum_old = float(hydro.params.get("energy_lum_lsun", 1.0))
    if lum_old <= 0.0:
        lum_old = 1.0
    with open(path, "w", newline="\n") as f:
        f.write("\n")
        f.write(" 24-FEB-2004                           !Format date\n")
        f.write("\n")
        f.write(f"   {r_core:22.15E}     {lum_old:22.15E} {nold:8d} {n:8d}\n")
        f.write("\n")
        for i in range(n - 1, -1, -1):
            r = max(r_1e10[i], 1e-30)
            e = max(ed[i], 1e-30)
            t = min(max(t_104[i], t_min_104), t_max_104)  # floor at t_min_104 (matches converter t_floor_K)
            v = v_kms[i]
            depth_idx = n - i
            f.write(
                f"   {r:22.15E} {v:22.15E} {e:22.15E} {t:22.15E} {1.0:22.15E} {v:22.15E} {1.0:22.15E} {depth_idx:6d} \n"
            )
            if nold > 0:
                for k in range(0, nold, 5):
                    chunk = [0.0] * min(5, nold - k)
                    f.write(" " + " ".join(f"{x:12.7E}" for x in chunk) + "\n")
            else:
                f.write(" 0.0\n")
    print(f"Wrote {path.name}: ND={n}  R(10^10 cm) {r_1e10[-1]:.4e} -> {r_core:.4e}  LUMOLD={lum_old:.4e} L_sun (for INIT_TEMP_V2)")


def write_fin_cal_grid_from_hydro(hydro: HydroData, path: Path) -> None:
    """
    Write FIN_CAL_GRID from HydroData (format read by SET_NEW_MODEL_ESTIMATES
    when LTE_EST=T and SN_HYDRO_MODEL=.FALSE.).

    Read format (set_new_model_estimates.f:222):
      READ(LUIN,*)T1,T1,T1,T1,T(I),T2,ED(I)
      ED(I)=ED(I)*T2

    We write: index, R(km), tau_dummy, V(km/s), T(10^4 K), Na(1/cm^3),
    Ne/Na, chi_dummy, chi_ratio_dummy, gamma_dummy

    CMFGEN reads: T1=index, T1=R, T1=tau, T1=V, T=T(10^4 K), T2=Ne/Na,
    ED=Na. Then ED(I) = Na * (Ne/Na) = Ne.

    Order: outer-to-inner (depth 1 = outermost, CMFGEN convention).
    """
    n = len(hydro.radius_cm)
    if n == 0:
        raise ValueError("write_fin_cal_grid_from_hydro: empty hydro")
    if hydro.temp_104k is None or len(hydro.temp_104k) < n:
        raise ValueError("write_fin_cal_grid_from_hydro: SN_HYDRO_DATA must have Temperature section")
    if hydro.ed_cgs is None or len(hydro.ed_cgs) < n:
        raise ValueError("write_fin_cal_grid_from_hydro: SN_HYDRO_DATA must have Electron density section")

    r_km = hydro.radius_cm / 1e5
    v_kms = hydro.vel_cms / 1e5  # preserve V sign (no abs)
    t_104 = hydro.temp_104k[:n]
    ne = hydro.ed_cgs[:n]

    # Use Na (atom density) directly from SN_HYDRO_DATA when available.
    # CMFGEN reads FIN_CAL_GRID as: T1=index, T1=R, T1=tau, T1=V, T=T(10^4 K),
    # T2=Ne/Na, ED=Na. Then ED(I) = Na * (Ne/Na) = Ne.
    # Using the actual Na from STELLA (rather than Na=Ne) gives the correct
    # ionization fraction (Ne/Na ~ 0.0006 at 4335 K) which CMFGEN uses for
    # its initial population estimate and optical depth computation.
    if hydro.na_cgs is not None and len(hydro.na_cgs) >= n:
        na = hydro.na_cgs[:n]
    else:
        # Fallback: compute Na from rho and mean atomic weight
        # Na = rho / (mu * m_u), mu ~ 1.3 for H/He mixture
        na = hydro.rho_cgs[:n] / (1.3 * 1.66e-24)
    ne_over_na = ne / np.maximum(na, 1e-30)

    with open(path, "w", newline="\n") as f:
        f.write("Index            R          TAU            V            T           Na        Ne/Na        Xross       Xr/Xes          Gam\n")
        for i in range(n - 1, -1, -1):  # outer-to-inner
            depth_idx = n - i
            r = max(r_km[i], 1e-30)
            v = v_kms[i]
            t = max(t_104[i], 1e-4)
            na_val = max(na[i], 1e-30)
            ne_ratio = max(ne_over_na[i], 1e-30)
            f.write(
                f"  {depth_idx:3d}  {r:12.4E}  {1e-4:12.4E}  {v:12.4E}  {t:12.4E}  "
                f"{na_val:12.4E}  {ne_ratio:12.4E}  {1e-8:12.4E}  {1.0:12.4E}  {0.87:12.4E}\n"
            )
    print(f"Wrote {path.name}: ND={n}  T(10^4 K) {t_104[0]:.4f} -> {t_104[-1]:.4f}  (for LTE_EST=T initial T/ED)")


def write_gammas_in_from_hydro(hydro: HydroData, path: Path) -> None:
    """
    Write GAMMAS_IN from HydroData (format read by GETELEC_V2).
    Format matches star_init/NT35000_logg425: leading blank, '!Number of depth points',
    '!Electron density', ED (10 per line, .3E), '!Radius (10^10cm)', '!Temperature (10^4 K)',
    then per species '!Atomic Number and descriptor.' and gamma rows (10 per line, .3E).

    GETELEC_V2 forms an initial log(ED) estimate as sum_s OLD_GAM(J,s) * n_s at each depth.
    If only hydrogen is listed but helium/iron carry most electrons, that sum can sit far
    below the template electron density at every depth; the overlap window then has NX=0 and
    LIN_INTERP is called with ND=0 (fatal). We therefore write one block per species in
    hydro.mass_fractions with gamma_s = Z_s * n_e / max(sum_k Z_k n_k, tiny), so
    sum_s gamma_s * n_s = n_e at each depth (fully stripped upper estimate for the denominator).
    """
    n = len(hydro.radius_cm)
    if n == 0:
        raise ValueError("write_gammas_in_from_hydro: empty hydro")
    if hydro.ed_cgs is None or len(hydro.ed_cgs) < n:
        raise ValueError("write_gammas_in_from_hydro: SN_HYDRO_DATA must have Electron density section")
    ed = hydro.ed_cgs[:n]
    r_1e10 = hydro.radius_cm / 1e10
    ed_oi = ed[::-1]
    r_oi = r_1e10[::-1]
    t_oi = (
        hydro.temp_104k[:n][::-1]
        if hydro.temp_104k is not None and len(hydro.temp_104k) >= n
        else np.ones(n)
    )
    rho = hydro.rho_cgs[:n]
    rho_oi = rho[::-1]

    # (Z, dum, n_cm3 outer-first) for each species with mass fractions and known (Z,A).
    species_rows: List[Tuple[int, str, np.ndarray]] = []
    for dum, arr in hydro.mass_fractions.items():
        if dum not in DUM_ATOMIC:
            continue
        z_el, a_amu = DUM_ATOMIC[dum]
        x_oi = np.clip(np.asarray(arr[:n], float)[::-1], 0.0, 1.0)
        n_cm3 = rho_oi * x_oi / (a_amu * AMU_G)
        species_rows.append((z_el, dum, n_cm3))

    species_rows.sort(key=lambda t: t[0])

    charge_sum = np.zeros(n, dtype=float)
    for z_el, _dum, n_cm3 in species_rows:
        charge_sum = charge_sum + z_el * n_cm3
    denom = np.maximum(charge_sum, np.maximum(1e-99, 1e-40 * np.max(charge_sum)))

    _fmt10 = lambda vals: "    " + "   ".join(f"{x:.3E}" for x in vals)

    with open(path, "w", newline="\n") as f:
        f.write("\n")
        f.write(f"  {n}                    !Number of depth points\n")
        f.write("\n")
        f.write(" !Electron density\n")
        for i in range(0, n, 10):
            f.write(_fmt10(ed_oi[i : i + 10]) + "\n")
        f.write("\n")
        f.write(" !Radius (10^10cm)\n")
        for i in range(0, n, 10):
            f.write(_fmt10(r_oi[i : i + 10]) + "\n")
        f.write("\n")
        f.write(" !Temperature (10^4 K)\n")
        for i in range(0, n, 10):
            f.write(_fmt10(t_oi[i : i + 10]) + "\n")
        f.write("\n")

        if not species_rows:
            f.write("   1.     HYD                           !Atomic Number and descriptor.\n")
            gam = np.ones(n)
            for i in range(0, n, 10):
                f.write(_fmt10(gam[i : i + 10]) + "\n")
            f.write("\n")
            print(f"Wrote {path.name}: ND={n}  (HYD-only gammas; no mass fractions in hydro)")
            return

        for z_el, dum, n_cm3 in species_rows:
            f.write(f"   {float(z_el):4.1f}     {dum:<29} !Atomic Number and descriptor.\n")
            gam = z_el * (ed_oi / denom)
            gam = np.clip(gam, 1e-50, 1e50)
            for i in range(0, n, 10):
                f.write(_fmt10(gam[i : i + 10]) + "\n")
            f.write("\n")

    dlist = ",".join(d for _z, d, _n in species_rows)
    print(f"Wrote {path.name}: ND={n}  species={dlist}")


def update_model_spec(
    model_spec_path: Path,
    nd: int,
    nc: Optional[int] = None,
    num_bnds: Optional[int] = None,
) -> List[str]:
    """
    Update MODEL_SPEC so ND matches SN_HYDRO_DATA, NP = ND + NC, and optionally NUM_BNDS.
    NUM_BNDS must match solver (VADAT SOL_METH): 1 = DIAG, 3 = TRIDIAG; must be 1 or 3 per CMFGEN.
    If nc is None, read NC from the file; otherwise use the provided value.
    Returns list of change description strings.
    """
    lines = read_lines(model_spec_path)
    key_idx = find_existing_keys(lines)
    changes: List[str] = []

    if "ND" not in key_idx:
        return changes

    if nc is None:
        tok = _get_key_token(lines, "NC")
        if tok is None:
            return changes
        try:
            nc = int(float(tok))
        except (ValueError, TypeError):
            return changes

    np_val = nd + nc

    old_nd = _get_key_token(lines, "ND")
    old_np = _get_key_token(lines, "NP")
    if old_nd is not None:
        set_or_append_keyword(lines, "ND", str(nd))
        changes.append(f"ND: {old_nd} -> {nd}")
    if "NP" in key_idx and old_np is not None:
        set_or_append_keyword(lines, "NP", str(np_val))
        changes.append(f"NP: {old_np} -> {np_val}  (ND+NC)")

    if num_bnds is not None and num_bnds in (1, 3):
        old_nb = _get_key_token(lines, "NUM_BNDS")
        set_or_append_keyword(lines, "NUM_BNDS", str(num_bnds))
        if old_nb is not None and str(old_nb).strip() != str(num_bnds):
            changes.append(f"NUM_BNDS: {old_nb} -> {num_bnds}  (consistent with N_TYPE/SOL_METH)")
        elif old_nb is None:
            changes.append(f"NUM_BNDS: (added) -> {num_bnds}  (consistent with N_TYPE/SOL_METH)")

    if changes:
        write_lines(model_spec_path, lines)
    return changes


def remove_species_from_model_spec(model_spec_path: Path, species_to_remove: List[str]) -> List[str]:
    """
    Comment out species in MODEL_SPEC by prefixing the ISF line with '!'.
    species_to_remove is a list of ISF key names, e.g. ['FeIV_ISF', 'FeV_ISF'].
    Returns list of change description strings.
    """
    if not species_to_remove:
        return []
    lines = read_lines(model_spec_path)
    changes: List[str] = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('!'):
            continue  # already commented
        for sp in species_to_remove:
            if f'[{sp}]' in stripped:
                lines[i] = '!' + line.lstrip('!')
                changes.append(f"Removed species {sp} (commented out in MODEL_SPEC)")
                break
    if changes:
        write_lines(model_spec_path, lines)
    return changes


def ensure_die_keys_for_model_spec(lines: List[str], model_spec_path: Path) -> List[str]:
    """
    Ensure VADAT contains DIE_* keys for ions enabled in MODEL_SPEC.
    Missing keys can trigger GET_KEY_STRING failures at startup.
    """
    if not model_spec_path.exists():
        return []

    text = model_spec_path.read_text(errors="ignore")
    changes: List[str] = []
    # RD_CONTROL_VARIABLES requires DIE_<ion> for every ion present in
    # MODEL_SPEC: TMP_KEY='DIE_'//ION_ID(ID) is read unconditionally for
    # all IDs where XzV_PRES is true. Derive keys generically from every
    # <prefix>_ISF tag rather than a fixed map so arbitrary ion sets work.
    for ion_tag in sorted(set(re.findall(r"\[(\w+)_ISF\]", text))):
        die_key = f"DIE_{ion_tag}"
        if _get_key_token(lines, die_key) is None:
            set_or_append_keyword(lines, die_key, "F,F")
            changes.append(f"{die_key}: (added) -> F,F")
    return changes


# ---------- abundance updates ----------
def depth_averaged_mass_fractions(mf: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    Compute depth-averaged mass fractions from SN_HYDRO_DATA.

    Very tiny positive fractions (e.g. 1e-50 from interpolation noise) are
    floored before normalisation so VADAT does not get absurd values like
    1e-100 in DUM/X. Such values are physically irrelevant but can cause
    numerical issues in CMFGEN.
    """
    if not mf:
        return {}
    ave: Dict[str, float] = {}
    min_frac = 1.0e-8  # floor for tiny but non-zero species
    for dum, arr in mf.items():
        if arr.size == 0:
            continue
        a = np.asarray(arr, float)
        # Treat non-finite or negative entries as zero, then average
        a = np.where(np.isfinite(a) & (a > 0.0), a, 0.0)
        if not np.any(a > 0.0):
            continue
        mean_val = float(np.mean(a))
        if mean_val <= 0.0:
            continue
        ave[dum] = max(mean_val, min_frac)
    total = sum(ave.values())
    if total <= 0.0:
        return {}
    return {dum: val / total for dum, val in ave.items()}

def update_abundances_mass_fraction(lines: List[str], ave_mass_frac: Dict[str, float], model_spec_path: Path = None) -> List[str]:
    changes: List[str] = []
    key_idx = find_existing_keys(lines)

    # Filter: only write abundances for elements that have ISF entries in MODEL_SPEC.
    # This prevents "ghost species" — elements with abundances but no model atom,
    # which CMFGEN would try to solve with default/minimal levels, causing convergence issues.
    allowed_dums = None
    if model_spec_path is not None and model_spec_path.exists():
        ms_text = model_spec_path.read_text(errors="ignore")
        # Map DUM labels to element names for ISF matching
        dum_to_isf_prefix = {
            "HYD": "HI", "HE": "HeI", "OXY": "OI", "SIL": "Sk", "CAL": "Ca",
            "IRON": "Fe", "NICK": "Nk", "COB": "Co", "CARB": "C", "NIT": "N",
            "NEON": "Ne", "MAG": "Mg", "ALUM": "Al", "SUL": "S",
            "ARG": "Ar", "TIT": "Tk", "CHRO": "Cr", "SCAN": "Sc",
            "PHOS": "P", "CHL": "Cl", "MAN": "Mn", "SOD": "Na",
        }
        allowed_dums = set()
        for dum, prefix in dum_to_isf_prefix.items():
            if f"[{prefix}" in ms_text and "_ISF]" in ms_text:
                allowed_dums.add(dum)
        # Always include NICK and COB for radioactive decay even if no ISF
        allowed_dums.add("NICK")
        allowed_dums.add("COB")

    for dum, frac in ave_mass_frac.items():
        if allowed_dums is not None and dum not in allowed_dums:
            continue
        key = f"{dum}/X"
        new_val = -float(frac)

        if key in key_idx:
            i = key_idx[key]
            raw = _strip_inline_comment(lines[i]).rstrip("\n")
            m = _VADAT_LINE_RE.match(raw.strip())
            if m:
                old_token = m.group(2)
                if _NUM_TOKEN_RE.match(old_token.strip()):
                    new_token = _format_like(old_token, new_val)
                else:
                    new_token = f"{new_val:.8E}".replace("E", "D")
                old_display = old_token.strip()
            else:
                new_token = f"{new_val:.8E}".replace("E", "D")
                old_display = "(unknown)"

            set_or_append_keyword(lines, key, new_token)
            changes.append(f"{key}: {old_display} -> {new_token}")
        else:
            new_token = f"{new_val:.7E}".replace("E", "D")
            set_or_append_keyword(lines, key, new_token)
            changes.append(f"{key}: (added) -> {new_token}")

    return changes


def renormalize_vadat_mass_fractions_to_unity(lines: List[str]) -> List[str]:
    """
    Ensure all DUM/X entries that are mass fractions (negative values) sum to 1 in absolute value.
    SET_ABUND_CLUMP requires 'Sum of mass fractions should be unity'. If the template had extra
    species or SN_HYDRO_DATA species were updated, the total can be != 1; scale all mass-frac
    entries so the sum is 1.
    """
    key_idx = find_existing_keys(lines)
    mass_frac_keys: List[str] = []
    mass_frac_vals: List[float] = []
    for key, i in key_idx.items():
        if not key.endswith("/X"):
            continue
        tok = _get_key_token(lines, key)
        if tok is None:
            continue
        try:
            v = float(tok.replace("D", "E").replace("d", "e"))
        except ValueError:
            continue
        if v < 0:  # mass fraction
            mass_frac_keys.append(key)
            mass_frac_vals.append(abs(v))
    if not mass_frac_keys:
        return []
    total = sum(mass_frac_vals)
    if total <= 0 or abs(total - 1.0) < 0.001:
        return []
    scale = 1.0 / total
    for key, abs_val in zip(mass_frac_keys, mass_frac_vals):
        new_abs = abs_val * scale
        new_token = f"{-new_abs:.8E}".replace("E", "D")
        set_or_append_keyword(lines, key, new_token)
    return [f"Mass fractions renormalized to sum = 1 (was {total:.6f}) for SET_ABUND_CLUMP"]


# ---------- init/final profiles ----------
# TRIDIAG solver (NUM_BNDS=3) re-enabled 2026-08-12: the previous revert was
# because DO_LAM_IT=T (continuous LAMBDA) never did full BA solves, so TRIDIAG
# couldn't couple depths. With alternating LAMBDA/linearization (the standard
# CMFGEN approach), TRIDIAG works correctly:
# the BA solve couples adjacent depth points, stabilising ionisation fronts.
# Reference: ttt.astro.su.se/~maer0651/M15-He-O-Ca-evo-new-4-all/1/ (converges
# in ~40 iterations with TRIDIAG/N_ON_J/NUM_BNDS=3).
def _num_bnds_for_n_type(n_type: str) -> int:
    """Return NUM_BNDS matching the solver for a given N_TYPE.

    N_ON_J -> 3 (TRIDIAG: couples adjacent depths in the BA matrix)
    G_ONLY -> 1 (DIAG: each depth solved independently)

    Note: N_TYPE and SOL_METH are decoupled when depth_independent_solver=True
    with non-monotonic velocity — N_ON_J is used for Eddington factor stability
    (H=0 at velocity reversals) while DIAG is used for solver stability.
    """
    n = n_type.strip().upper()
    if n == "N_ON_J":
        return 3
    return 1


def apply_profile(
    lines: List[str],
    profile: str,
    convergence_mode: str = "default",
    num_isotopes: int = 0,
    depth_independent_solver: bool = True,
    user_overrides: Optional[dict] = None,
    allow_t_release: bool = False,
    thin_t_release: bool = False,
) -> Tuple[List[str], int]:
    """
    Apply final VADAT profile for CMFGEN convergence run.

    Two profiles remain after the INIT/STABLE removal:
      - single_final: LTE cold-start (first day or single snapshot; no prior populations)
      - final:        warm restart from *_IN files (subsequent days in a time series)

    num_isotopes: from SN_HYDRO_DATA header. If > 0, INC_RAD_DECAYS must be T so CMFGEN
    reads NUC_DECAY_DATA and RD_SN_DATA can match isotope blocks to the ISO table.

    depth_independent_solver: if True (default), use DIAG (NUM_BNDS=1) for SOL_METH.
        N_TYPE is always N_ON_J (stable at velocity reversals where H=0).
        DIAG prevents inter-depth instability propagation; TRIDIAG couples adjacent depths.
        If False, use TRIDIAG (NUM_BNDS=3, N_TYPE=N_ON_J) which couples adjacent depths.
    """
    if profile not in {"final", "single_final"}:
        raise ValueError(f"Unknown profile: {profile} (only 'final' and 'single_final' supported)")
    if convergence_mode not in {"default", "stable"}:
        raise ValueError(f"Unknown convergence_mode: {convergence_mode}")

    inc_rad_decays = "T" if num_isotopes > 0 else "F"

    # ---------------------------------------------------------------------
    # Runtime overrides for stability experiments / batch resubmissions.
    #
    # These are intentionally environment-variable driven so Slurm job scripts
    # (or submit wrappers) can adjust behavior without editing model files.
    #
    # - CMFGEN_N_OB_INS:
    #     Overrides VADAT [N_OB_INS] (outer-boundary insertion points). CMFGEN
    #     suggests a finer outer grid as the preferred fix for J-instabilities.
    # - CMFGEN_N_IB_INS:
    #     Overrides VADAT [N_IB_INS] (inner-boundary insertion points). For AMR
    #     stability, use 1 instead of 3 to avoid boundary artefacts with steep gradients.
    # ---------------------------------------------------------------------
    # Environment variable overrides removed - use hardcoded values in profile dictionaries

    if profile == "single_final":
        # Single-final: skip INIT entirely; start from LTE populations.
        # LTE_EST=T: forces DC_INTERP_METHOD='LTE' (rd_control_variables.f:847),
        # so DC_METH is irrelevant. With SN_HYDRO_MODEL=.FALSE. (VEL_LAW=12),
        # T/ED come from FIN_CAL_GRID (created by write_fin_cal_grid_from_hydro).
        # Pure Lambda iteration (no Non-Lambda auto-switch), non-coherent ES.
        # T fixed to hydro, TRIDIAG solver (couples adjacent depths).
        # Ng acceleration enabled: activates after IBEG_NG Lambda iterations once
        # the max population change falls below BEG_NG (guards early blowups).
        desired = {
            "LTE_EST": "T",
            "DC_METH": "TR",
            "SN_T_OPT": "USE_HYDRO",
            "INC_RAD_DECAYS": inc_rad_decays,
            "INC_SHOCK_POWER": "F",  # add shock heating beyond hydro T (F = heating already in T(r))
            "INC_AD": "F",           # Static snapshot: no adiabatic cooling (DO_DDT=F, DO_HYDRO=F)
            "INC_ADV": "F",          # Static snapshot: no wind advection
            "ALL_FREQ": "F",  # apply the grid-accuracy test only where needed, not at all nu
            "LAM_SET": "F",  # keep GLOBAL_LINE=SOB active; suppress overlap warning

            "RD_CF_FILE": "F",
            "MIN_CF": "3.0D-3",  # lowest continuum frequency in 1e15 Hz (far-IR edge)
            "MAX_CF": "2.0D2",  # highest continuum frequency in 1e15 Hz (far-UV edge)
            "FRAC_SP": "1.10D0",  # max ratio of adjacent continuum-grid frequencies
            "AMP_FAC": "1.05D0",  # allowed growth of FRAC_SP toward high nu
            "MAX_BF": "0.10D0",  # max grid spacing near bound-free edges (1e15 Hz)
            "CROSS": "F",
            "V_CROSS": "500.0D0",
            "SIG_GAU_KMS": "500.0D0",

            # Relativistic configuration (all models are static snapshots)
            **get_rel_config(),

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "20.0D0",  # turbulent velocity (km/s); sets line Doppler width
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)
            "JC_W_EDD": "T",         # MUST be T: opens EDDFACTOR file (unit 35) via OPEN_RW_EDDFACTOR.
            "JBAR_W_EDD": "T",  # Eddington factors for line Jbar
                                     # Without this, cmfgen_sub.f:1348 condition (ACCURATE.OR.EDD_CONT.OR.EDD_LINECONT)
                                     # is FALSE (INC_GRID=F, JBAR_W_EDD=F), so unit 35 is never opened.
                                     # COMP_J_BLANK then WRITE(LU_EDD,REC=...) to unopened unit → fort.35 crash
                                     # (forrtl severe 37: inconsistent record length).
            "THK_CONT": "F",  # Thin outer continuum boundary: wind/SN-ejecta models are
                                     # optically thin at the outer grid edge.  A thick outer BC
                                     # (free-streaming suppressed) is the wrong physics; set
                                     # THK_CONT=T via vadat_overrides for optically thick models.
            "THK_LINE": "F",            # No thick line BC — inner boundary is optically thin
            "DIF": "F",               # No diffusion approximation — ZERO_FLUX inner boundary (avoids negative DBB from non-monotonic T)
            "IB_METH": "ZERO_FLUX",   # Zero flux at inner boundary: carries no LSTAR flux (dead boundary; X_SOB 40d converged to L=0)
            "IB_STAB": "0.1D0",       # Inner boundary stabilization
            "STOP_IF_MW": "F",        # Don't halt on "poor flux" warnings
            "COH_ES": "T",  # electron-scattering source term (T = coherent, F = non-coherent)
            "MIX_COH": "F",

            # Static snapshot: T fixed to hydro (FIX_T=T); Ne solved from
            # charge conservation (FIX_NE=F is the project default for all
            # runs, adopted 25-Sep-2026). The hydro Ne is a Saha/LTE guess,
            # not the NLTE ionization balance. The user's vadat_overrides
            # take precedence over these via the user_overrides parameter.
            "FIX_T": "T",            # T fixed to hydro (static snapshot)
            "FIX_T_AUTO": "F",       # No auto-release (static snapshot)
            # Hydro Ne is a rough estimate (Saha/ LTE guess in the input
            # file, not NLTE); holding it fixed propagates that error
            # into every ionization balance.  Let CMFGEN solve Ne from
            # charge conservation instead (user directive 18-Sep-2026).
            "FIX_NE": "F",
            "FIX_IMP": "F",  # impurity populations iterated like major species
            "FIX_INB_T": "F",
            "FIX_X_DPTH": "0",
            "TAU_SCL_T": "0.0D0",    # No tau-scaling for T (FIX_T=T; T from hydro)

            # TRIDIAG solver: NUM_BNDS=3 (couples adjacent depths)
            "SCALE_OPT": "MAJOR",     # MAJOR: major species drive scaling
            "LAM_SCALE_OPT": "NONE",  # NONE: disable LIMIT clipping — allow large decreases for LTE→NLTE transition in non-LTE systems
            "MAX_LIN": "1.0D1",  # max correction factor allowed per linearization step
            "MAX_LAM": "1.0D10",  # max Lambda-iteration correction factor (1e10 = effectively unclamped)
                                     # SN snapshots have wildly wrong LTE initial populations — need large MAX_LAM to traverse the gap.
            "MAX_CHNG": "1.0D150",  # abort if a correction exceeds this % change; MAX_LAM is the real limiter
            "AUTO_SMOOTH": "F",     # Disable population smoothing for zones with STEQ>0.998 or <-1e5
            "COMP_BA": "F",
            "N_FIX_BA": "2",  # recompute the BA matrix every N iterations
            "BA_CHK_FAC": "1.0D-05",  # fractional accuracy demanded of the BA solution
            "FIX_BA": "1.0D1",  # switch to BA solver once max % change falls below this
            "LAM_VAL": "2.0D2",  # keep Lambda-iterating while max % change exceeds this
            "NUM_LAM": "2",  # Lambda iterations between full linearizations

            "GLOBAL_LINE": "SOB",
            "GF_CUT": "2.0D-3",  # discard lines weaker than this normalized gf-cut
            "GF_LEV_CUT": "9",
            "AT_CUT": "10",  # apply GF_CUT only to species with atomic number >= this
            "MIN_TRANS": "9",  # discard levels with fewer bound-bound transitions than this
            "CHK_L_POS": "T",
            "NEG_OPAC_OPT": "SRCE_CHK",
            # SCL_SL_OPAC=T: scale SL line opacities/emissivities so line cooling from
            # a Super-Level is zero when the net rate is zero. Prevents artificial
            # heating/cooling spikes at contact discontinuities (SOB line mode).
            # MUST set SCL_LN=F: source code (rd_control_variables.f:745) aborts if
            # both SCL_SL_OPAC and SCL_LN are T.
            "SCL_SL_OPAC": "T",
            "SCL_LN": "F",
            "SCL_LN_FAC": "0.5D0",           # Manual default (unused when SCL_LN=F; kept for reference)
            "ES_FAC": "0.6D0",  # electron-scattering source weight (only used when MIX_COH=T)

            # T_ON_R=T: interpolate T on R grid (not spherical tau grid). Required for
            # non-monotonic T profiles at shock density inversions (SOB line mode).
            # When DC_METH='TR', the source default is already T, but we set it explicitly.
            "T_ON_R": "T",
            # LTE cold start: LTE_EST=F (reads T_IN/GAMMAS_IN), DC_METH=LTE (LTE DCs)
            "LIN_INT": "F",           # DO NOT interpolate from old CMFGEN model — use STELLA snapshot T/pops
            "IT_ON_T": "F",           # No initial T iteration — hydro T used as-is (user directive 18-Sep-2026)
            "POP_SCALE": "T",
            # T_EXC: set from SN_HYDRO_DATA photospheric temperature (computed below)
            "T_MIN": "2.0D-01",  # temperature floor in 1e4 K units (2000 K)
            "GREY_TAU": "0.5D0",  # tau depth anchoring the grey-T estimate
            "METHOD": "ZERO",  # depth-correction interpolation (ZERO = constant between depths)
            "INC_RAY": "T",  # include Rayleigh scattering opacity
            # Linear interpolation when reading SN_HYDRO_DATA; the manual recommends
            # this for shell/interacting profiles with rapid density/abundance changes.
            "LIN_RD_SN_DATA": "T",

            "REV_RGRID": "F",       # Disabled: no runtime R-grid revision — SN_HYDRO_DATA arrives pre-optimized (user directive 18-Sep-2026)
            "RD_IN_R_GRID": "F",
            "N_IB_INS": "0",          # AMR handles boundary refinement
            "N_OB_INS": "0",          # No boundary insertions — SN_HYDRO_DATA grid arrives pre-optimized (user directive 18-Sep-2026)
            "OB_RAT": "1.5D0",
            "IB_RAT": "2.0D0",
            # N_TYPE is set by get_rel_config() — do not override here

            # INC_GRID=T would create an enhanced spatial grid for solving the RTE.
            # ALL_FREQ=F: apply the enhanced spatial grid only where needed.
            # ST_INT/END_INT: bracket the forward/reverse shock region (depth indices,
            # CMFGEN outside-in convention). Defaults 40-55; override via vadat_overrides.
            # WARNING: with non-monotonic velocity, EXTEND_VTSIGMA interpolates sigma
            # at ALL depths. Monitor OUTGEN for sigma artefacts at velocity reversals.
            "INC_GRID": "F",
            "VFRAC_FG": "20000.0D0", # Max velocity spacing along ray in km/s (formal-solution path; set at line ~2860)
            "OB_METH": "HONJ",

            "ACC_END": "5.0D0",
            "N_INS": "0",            # No grid insertion — prevents extra depths at inner boundary causing BA singularity
            "ST_INT": "40",          # Enhanced grid start depth (bracket forward/reverse shock)
            "END_INT": "55",         # Enhanced grid end depth (bracket forward/reverse shock)
            "ND_QUAD": "0",
            "INTERP_TYPE": "LOG",
            "N_PAR": "2000",  # frequency intervals between BA recomputations

            "MAX_DOP": "6.0D0",  # resonance-zone half-width in Doppler units
            "FRAC_DOP": "1.0D0",  # frequency spacing inside resonance zone (Doppler units)
            "dV_CMF_PROF": "200.0D0",  # comoving spacing in line cores (km/s)
            "dV_CMF_WING": "300.0D0",  # comoving spacing in line wings (km/s)
            "ES_WING_EXT": "12000.0D0",  # Keep 12000: wide e.s. wings for IIn-like models
            "R_CMF_WING_EXT": "3.0D0",  # coherent-ES wing extent in VINF units; pipeline floors 4.0 when COH_ES=F

            "OBS_EXT_RAT": "1.1D0",  # extend obs grid this factor beyond the bluest line edge
            "dV_OBS_PROF": "100.0D0",  # observer spacing in line cores (km/s)
            "dV_OBS_WING": "400.0D0",  # observer spacing in line wings (km/s)
            "dV_OBS_BIG": "2000.0D0",  # observer spacing in line-free regions (km/s)

            "EPS_TERM": "5.0D-02",  # iterate until the max % population change falls below this
            "ACC_F": "1.0D-04",  # accuracy demanded of the Eddington-factor solve

            # Alternating LAMBDA/linearization (standard CMFGEN approach).
            # Pure LAMBDA mode: DO_LAM_IT=T, DO_LAM_AUTO=F.  The generic
            # profile leaves FIX_IMP=F; the SOB recipe also leaves it F so
            # FIXPOP_IN_BA_V3 cannot zero trace-population corrections.
            "USE_LAM_ES": "T",      # J=B(T) Lambda iteration — bypasses moment solver for non-monotonic velocity
            "DO_LAM_IT": "T",       # Start in LAMBDA mode (stable cold-start from LTE)
            "DO_LAM_AUTO": "T",  # auto-switch from Lambda to full iterations once corrections shrink
            # Ng acceleration DISABLED by default (13-Oct-2026, X_SOB 40d):
            # scheduled kicks fired every ITS/NG iterations regardless of
            # residual size and overshot at the interaction depth by up to
            # ~30%, pinning MAXCH near 0.5% above EPS_TERM once dead-level
            # demand was gated. The unaccelerated map contracted ~0.65x/iter.
            # Re-enable per-model via vadat_overrides DO_NG=T; the safety
            # mechanisms remain: per-depth rejection, global cancellation,
            # post-Ng LAMBDA fallback, CHK_NG=T, NG_CLEAN>=4, NG_ABSTAIN=T.
            "DO_NG": "F",             # Disabled: NG kicks pinned the near-convergence
                                        # residual floor on X_SOB; opt in via overrides.
            "CHK_NG": "T",            # Check NG acceleration validity
            "BEG_NG": "5.0D0",  # start Ng acceleration when max % change falls below this
            "IBEG_NG": "15",  # minimum iterations before first Ng acceleration
            "BW_NG": "40",  # depths per Ng acceleration band
            "ITS/NG": "8",  # run Ng every N-th iteration

            # DO_AV is NOT set: with alternating LAMBDA/linearization, population
            # averaging is not needed (the BA solve couples depths and prevents
            # the oscillations that required DO_AV in continuous LAMBDA mode).
        }

    elif profile == "final":
        # Warm restart from *_IN files (time-series: subsequent days seed from previous day).
        # LTE_EST=F: read populations from *_IN files (renamed from previous run's *_OUT).
        # LIN_INT=T: interpolate populations between depth points.
        # Same damping/solver as single_final for consistency.
        desired = {
            "SN_T_OPT": "USE_HYDRO",
            "INC_RAD_DECAYS": inc_rad_decays,
            "INC_SHOCK_POWER": "F",
            "INC_AD": "F",           # Static snapshot: no adiabatic cooling (DO_DDT=F, DO_HYDRO=F)
            "INC_ADV": "F",          # Static snapshot: no wind advection
            "ALL_FREQ": "F",  # apply the grid-accuracy test only where needed, not at all nu
            "LAM_SET": "F",  # keep GLOBAL_LINE=SOB active; suppress overlap warning

            # Continuum frequency grid limits and spacing
            "RD_CF_FILE": "F",      # Build frequency grid internally
            "MIN_CF": "3.0D-3",  # lowest continuum frequency in 1e15 Hz (far-IR edge)
            "MAX_CF": "2.0D2",  # highest continuum frequency in 1e15 Hz (far-UV edge)
            "FRAC_SP": "1.10D0",  # max ratio of adjacent continuum-grid frequencies
            "AMP_FAC": "1.05D0",  # allowed growth of FRAC_SP toward high nu
            "MAX_BF": "0.10D0",  # max grid spacing near bound-free edges (1e15 Hz)
            "CROSS": "F",          # Sparse chi/eta with interpolation
            "V_CROSS": "500.0D0",   # 500 km/s for 500 km/s feature resolution
            "SIG_GAU_KMS": "500.0D0",  # match V_CROSS for consistency

            # Relativistic configuration (all models are static snapshots)
            **get_rel_config(),
            "JC_W_EDD": "T",         # MUST be T: opens EDDFACTOR file (unit 35) — see init profile comment
            "JBAR_W_EDD": "T",  # Eddington factors for line Jbar
            "THK_CONT": "F",  # Thin outer continuum boundary: wind/SN-ejecta models are
                                     # optically thin at the outer grid edge (see single_final).
            "THK_LINE": "F",            # No thick line BC at outer boundary
            "DIF": "F",               # No diffusion approximation — ZERO FLUX inner boundary for self-powered models
            "IB_METH": "ZERO_FLUX",   # Zero flux inner boundary (I+=I-); carries no LSTAR flux - dead boundary, needs LSTAR=0 or DIFFUSION
            "IB_STAB": "0.1D0",       # Inner boundary stabilization (0.1 for ZERO_FLUX; manual default for non-diffusion)
            "STOP_IF_MW": "F",        # Don't halt on "poor flux" warnings (diagnostic with ZERO_FLUX + VEL_LAW=12; SN_MODEL=.FALSE. normalizes by L[ND])
            "COH_ES": "T",  # electron-scattering source term (T = coherent, F = non-coherent)
            "MIX_COH": "F",
            "ES_FAC": "0.6D0",

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "20.0D0",  # turbulent velocity (km/s); sets line Doppler width
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)

            # Static snapshot: T fixed to hydro (FIX_T=T); Ne solved from
            # charge conservation (FIX_NE=F is the project default for all
            # runs, adopted 25-Sep-2026). The user's vadat_overrides take
            # precedence over these via the user_overrides parameter.
            "FIX_T": "T",            # T fixed to hydro (static snapshot)
            "FIX_T_AUTO": "F",       # No auto-release (static snapshot)
            # Hydro Ne is a rough estimate; solve Ne from charge
            # conservation (user directive 18-Sep-2026, same as
            # single_final).
            "FIX_NE": "F",
            "FIX_IMP": "F",  # impurity populations iterated like major species
            "FIX_INB_T": "F",
            "FIX_INB_NE": "F",       # Ne free at all depths in final
            "FIX_X_DPTH": "0",        # T fixed at ALL depths
            "TAU_SCL_T": "0.0D0",    # No tau-scaling for T (FIX_T=T; T from hydro)

            # DIAG solver: NUM_BNDS=1 (each depth independent; TRIDIAG reverted 2026-08-01)
            "SCALE_OPT": "MAJOR",     # MAJOR: major species drive scaling
            "LAM_SCALE_OPT": "NONE",  # NONE: disable LIMIT clipping — allow large decreases for LTE→NLTE transition
            "MAX_LIN": "1.0D1",  # max correction factor allowed per linearization step
            "MAX_LAM": "1.0D10",  # max Lambda-iteration correction factor (1e10 = effectively unclamped)
            "MAX_CHNG": "1.0D150",  # abort if a correction exceeds this % change; MAX_LAM is the real limiter
            "AUTO_SMOOTH": "F",     # Disable population smoothing for zones with STEQ>0.998 or <-1e5
            "COMP_BA": "F",
            "N_FIX_BA": "2",  # recompute the BA matrix every N iterations
            "BA_CHK_FAC": "1.0D-05",  # fractional accuracy demanded of the BA solution
            "FIX_BA": "1.0D1",  # switch to BA solver once max % change falls below this
            "LAM_VAL": "2.0D2",  # keep Lambda-iterating while max % change exceeds this
            "NUM_LAM": "2",  # Lambda iterations between full linearizations

            "GLOBAL_LINE": "SOB",
            "GF_CUT": "2.0D-3",  # discard lines weaker than this normalized gf-cut
            "GF_LEV_CUT": "9",
            "AT_CUT": "10",  # apply GF_CUT only to species with atomic number >= this
            "MIN_TRANS": "9",  # discard levels with fewer bound-bound transitions than this
            "THK_LINE": "F",            # No thick line BC at outer boundary
            "CHK_L_POS": "T",
            "NEG_OPAC_OPT": "SRCE_CHK",
            # SCL_SL_OPAC=T: scale SL line opacities/emissivities for heating/cooling consistency.
            # MUST set SCL_LN=F: source code aborts if both SCL_SL_OPAC and SCL_LN are T.
            "SCL_SL_OPAC": "T",
            "SCL_LN": "F",
            "SCL_LN_FAC": "0.5D0",           # Manual default (unused when SCL_LN=F; kept for reference)
            "MIX_COH": "F",
            "ES_FAC": "0.6D0",  # electron-scattering source weight (only used when MIX_COH=T)

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "20.0D0",  # turbulent velocity (km/s); sets line Doppler width
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)

            # T_ON_R=T: interpolate T on R grid (not spherical tau grid). Required for
            # non-monotonic T profiles at shock density inversions (SOB line mode).
            "T_ON_R": "T",
            "LTE_EST": "F",           # read population levels from *_IN files
            "USE_FIXED_J": "F",     # Let J vary in FINAL — populations are settled, large MAX_LAM allows convergence
            "IT_ON_T": "F",           # No initial T iteration — hydro T used as-is (user directive 18-Sep-2026)
            "POP_SCALE": "T",
            "T_MIN": "2.0D-01",  # temperature floor in 1e4 K units (2000 K)
            # T_EXC: set from SN_HYDRO_DATA photospheric temperature (computed below)
            "GREY_TAU": "0.5D0",  # tau depth anchoring the grey-T estimate
            "METHOD": "ZERO",  # depth-correction interpolation (ZERO = constant between depths)
            "INC_RAY": "T",  # include Rayleigh scattering opacity

            "REV_RGRID": "F",       # Disabled: no runtime R-grid revision — SN_HYDRO_DATA arrives pre-optimized (user directive 18-Sep-2026)
            "RD_IN_R_GRID": "F",
            "N_IB_INS": "0",          # AMR handles boundary refinement
            "N_OB_INS": "0",          # No boundary insertions — SN_HYDRO_DATA grid arrives pre-optimized (user directive 18-Sep-2026)
            "OB_RAT": "1.5D0",
            "IB_RAT": "2.0D0",
            # N_TYPE is set by get_rel_config() — do not override here

            # INC_GRID=T: create enhanced spatial grid for solving the RTE.
            # WARNING: with non-monotonic velocity, EXTEND_VTSIGMA interpolates sigma
            # at ALL depths. Monitor OUTGEN for sigma artefacts at velocity reversals.
            "INC_GRID": "T",
            "VFRAC_FG": "20000.0D0", # Max velocity spacing along ray in km/s (formal-solution path)
            "OB_METH": "HONJ",

            "ALL_FREQ": "F",  # apply the grid-accuracy test only where needed, not at all nu
            "ACC_END": "5.0D0",
            "N_INS": "0",            # No grid insertion
            "ST_INT": "40",          # Enhanced grid start depth (bracket forward/reverse shock)
            "END_INT": "55",         # Enhanced grid end depth (bracket forward/reverse shock)
            "ND_QUAD": "0",
            "INTERP_TYPE": "LOG",
            "N_PAR": "2000",  # frequency intervals between BA recomputations

            "MAX_DOP": "6.0D0",  # resonance-zone half-width in Doppler units
            "FRAC_DOP": "1.0D0",  # frequency spacing inside resonance zone (Doppler units)
            "dV_CMF_PROF": "200.0D0",  # comoving spacing in line cores (km/s)
            "dV_CMF_WING": "300.0D0",  # comoving spacing in line wings (km/s)
            "ES_WING_EXT": "12000.0D0",  # Keep 12000: wide e.s. wings for IIn-like models
            "R_CMF_WING_EXT": "3.0D0",  # coherent-ES wing extent in VINF units; pipeline floors 4.0 when COH_ES=F

            "OBS_EXT_RAT": "1.1D0",  # extend obs grid this factor beyond the bluest line edge
            "dV_OBS_PROF": "100.0D0",  # observer spacing in line cores (km/s)
            "dV_OBS_WING": "400.0D0",  # observer spacing in line wings (km/s)
            "dV_OBS_BIG": "2000.0D0",  # observer spacing in line-free regions (km/s)

            "EPS_TERM": "5.0D-02",  # iterate until the max % population change falls below this
            "ACC_F": "1.0D-04",  # accuracy demanded of the Eddington-factor solve

            "DO_LAM_IT": "T",
            "DO_LAM_AUTO": "T",  # auto-switch from Lambda to full iterations once corrections shrink

            # Ng acceleration DISABLED for FINAL (13-Oct-2026): same default
            # as single_final -- scheduled kicks overshot near convergence on
            # X_SOB and pinned MAXCH above EPS_TERM. Safety machinery
            # (per-depth rejection, global cancellation, post-Ng LAMBDA
            # fallback, CHK_NG=T, NG_CLEAN, NG_ABSTAIN=T) stays armed for
            # models that opt back in via vadat_overrides DO_NG=T.
            "DO_NG": "F",             # Disabled: NG kick overshoot pinned the X_SOB
                                        # residual floor near convergence; opt in via overrides.
            "CHK_NG": "T",            # Check NG acceleration validity
            "BEG_NG": "5.0D0",  # start Ng acceleration when max % change falls below this
            "IBEG_NG": "15",  # minimum iterations before first Ng acceleration
            "BW_NG": "40",  # depths per Ng acceleration band
            "ITS/NG": "8",  # run Ng every N-th iteration
            "DO_AV": "F",  # population averaging between iterations
            "NOSC_AV": "4",          # Manual default: 4 consecutive oscillations before averaging
            "ITS/AV": "8",           # Manual default: 8 iterations between averaging attempts
        }

    if convergence_mode == "stable":
        # Conservative damping for NLTE-heavy IIn models where LIMIT clipping repeats
        # and lambda corrections remain large for many iterations.
        # These overrides apply on top of the profile defaults for extra damping.
        # Stable-mode overrides: aligned with CMFGEN manual (full_descr.pdf) recommendations.
        # These overrides apply on top of the profile defaults for both final and single_final
        # (both use --convergence-mode stable in submit_cmfgen.sh).
        #
        # Previous overrides (MAX_LAM=1.1D0, NUM_LAM=40) caused an infinite Lambda
        # loop: trace species needing 10000% correction could only change 10% per iteration,
        # and without Ng acceleration the model never converged (549 iterations, X_2 10d job).
        # Ng is now enabled with robust safety (BEG_NG=5%, IBEG_NG=15, per-depth 10x rejection).
        #
        # CMFGEN manual (full_descr.pdf) recommendations for slow convergence:
        #   MAX_LAM: "Typically adopt 10^5. For a new species with very low abundance,
        #             this can be set as high as 10^20."
        #   MAX_LIN: "Typically adopt 10."
        #   LAM_VAL: "Typically do a Lambda-iteration if % change is greater than 400."
        #   NUM_LAM: "Typically choose 2."
        #   N_FIX_BA: "recommended value is 2."
        #   BA_CHK_FAC: "Typically we adopt 5% for the largest fractional change."
        #   DO_NG: "Ng acceleration is an important tool for obtaining convergence"
        #          (especially with NUM_BNDS=1 / DIAG solver).
        #   IBEG_NG: "at least 15 (or more) iterations have been completed"
        #   BEG_NG: "Typically choose 5% to 10%"
        #   ITS/NG: "We typically adopt 8."
        #
        # Strategy: Very slow LTE→NLTE transition via MAX_LAM=1.01.
        # MAX_LAM=1.01D0 (1.01x) for all profiles — limits trace species to 1.01x growth
        #   per iteration. The J-population cascade is fundamental: trace species need
        #   to grow ~1e14x from LTE, and their line opacity changes destabilize J at
        #   other depths. With MAX_LAM=3, cascade at iter 22; MAX_LAM=1.5, cascade at
        #   iter 58; USE_FIXED_J=T (holds continuum opacity only), cascade at iter 58.
        #   With MAX_LAM=1.01, each population correction is tiny (1.01x), so the model
        #   can absorb the cascade gradually without exploding. Convergence takes
        #   ~4000 iterations (~20 hours at 19 sec/iter) but the 24h wall time + requeue
        #   mechanism handles this.
        #   LIT_LIM=1-MAX_LAM=-0.01 caps trace species increases at 1.01x.
        #   LIMIT only caps DECREASES (STEQ>1.1→0.999); INCREASES (STEQ<0) are capped
        #   by LIT_LIM. MAJOR scaling excludes trace species from the scale factor.
        # AUTO_SMOOTH=T — enable population smoothing when STEQ<-1e5 or >0.998.
        # smooth_pops_as_we_iterate.f extrapolates from deeper zones for wild corrections.
        # FIX_IMP=F for all profiles — freezing trace species suppresses the
        # BA corrections those populations need during the LTE-to-NLTE move.
        # FIX_NE=F for all profiles — STELLA Saha Ne is unreliable for X_1 (predates Yat*AZ fix).
        # FIX_T=T for all profiles — T fixed to hydro (moment solver + SIGMA clamp handles non-monotonic velocity).
        stable_overrides = {
            "single_final": {
                # Generic stable profile: alternating Lambda/full with NG
                # enabled.  SOB runs stamp the optimized
                # NON_MONO_SOB_RECOMMENDED schedule after this block, so
                # values such as FIX_IMP, NUM_LAM and the BA gates differ
                # there by design.
                #
                # Settings follow CMFGEN manual (full_descr.pdf) defaults:
                #   MAX_LAM=1D5: "Typically adopt 10^5" — allows fast population adjustment.
                #   MAX_LIN=1D1: "Typically adopt 10" — 10x max correction for BA mode.
                #   LAM_VAL=4D2: "Typically do a Lambda-iteration if % change > 400"
                #   NUM_LAM=2:   "Typically choose 2" — frequent BA coupling for stability.
                #   N_FIX_BA=2:  generic Jacobian-refresh interval; SOB uses 10.
                #   AUTO_SMOOTH=F: Manual warns "in some cases it can make convergence worse"
                #   DO_AV=F:     Manual default is FALSE — our issue is divergence, not oscillation.
                #   DO_NG=T:     Ng polynomial extrapolation after BEG_NG/IBEG_NG.
                "MAX_LIN": "1.0D1",  # max correction factor allowed per linearization step
                "MAX_LAM": "1.0D10",   # raised to 1e10 per user directive 15-Sep-2026; large jumps for LTE→NLTE transition
                "LAM_VAL": "2.0D2",  # keep Lambda-iterating while max % change exceeds this
                "NUM_LAM": "2",  # Lambda iterations between full linearizations
                "N_FIX_BA": "2",  # recompute the BA matrix every N iterations
                "BA_CHK_FAC": "1.0D-05",  # fractional accuracy demanded of the BA solution
                "FIX_BA": "1.0D1",  # switch to BA solver once max % change falls below this

                "MAX_CHNG": "1.0D150",  # abort if a correction exceeds this % change; MAX_LAM is the real limiter
                "AUTO_SMOOTH": "F",     # Disable population smoothing for zones with STEQ>0.998 or <-1e5
                "FIX_NE": "F",          # Ne free — STELLA Saha Ne is unreliable
                "FIX_IMP": "F",  # impurity populations iterated like major species
                "FIX_T": "T",            # T pinned to the hydro profile for all runs (14-Oct-2026 policy)
                "FIX_T_AUTO": "F",       # No T auto-release
                "TAU_SCL_T": "0.0D0",    # No tau-scaling for T (FIX_T=T)
                "FIX_INB_T": "F",        # Don't fix T at inner boundary separately
                "DO_NG": "F",           # Disabled by default (13-Oct-2026): NG kicks
                                        # overshot near convergence on X_SOB and held
                                        # MAXCH above EPS_TERM. Opt in via
                                        # vadat_overrides DO_NG=T.
                "CHK_NG": "T",          # Check NG acceleration validity
                "BEG_NG": "5.0D0",  # start Ng acceleration when max % change falls below this
                "IBEG_NG": "15",  # minimum iterations before first Ng acceleration
                "BW_NG": "40",  # depths per Ng acceleration band
                "ITS/NG": "8",  # run Ng every N-th iteration
                "DO_LAM_IT": "T",       # Start in LAMBDA mode (stable cold-start from LTE)
                "DO_LAM_AUTO": "T",    # Auto-switch from LAMBDA to full iterations once corrections are small
                "DO_AV": "F",  # population averaging between iterations
                "METHOD": "ZERO",  # depth-correction interpolation (ZERO = constant between depths)
                "SCALE_OPT": "MAJOR",   # Manual: "Preferred option is MAJOR"
                "LAM_SCALE_OPT": "NONE", # NONE: disable LIMIT clipping — allow large decreases for LTE→NLTE transition
                "REV_RGRID": "F",       # Disabled: no runtime R-grid revision — SN_HYDRO_DATA arrives pre-optimized (user directive 18-Sep-2026)
            },
            "final": {
                "MAX_LIN": "1.0D1",  # max correction factor allowed per linearization step
                "MAX_LAM": "1.0D10",   # raised to 1e10 per user directive 15-Sep-2026; large jumps for LTE→NLTE transition
                "LAM_VAL": "2.0D2",  # keep Lambda-iterating while max % change exceeds this
                "NUM_LAM": "2",  # Lambda iterations between full linearizations
                "N_FIX_BA": "2",  # recompute the BA matrix every N iterations
                "BA_CHK_FAC": "1.0D-05",  # fractional accuracy demanded of the BA solution

                "FIX_BA": "1.0D1",  # switch to BA solver once max % change falls below this
                "MAX_CHNG": "1.0D150",  # abort if a correction exceeds this % change; MAX_LAM is the real limiter
                "AUTO_SMOOTH": "F",     # Disable population smoothing for zones with STEQ>0.998 or <-1e5
                "FIX_NE": "F",          # Ne free — STELLA Saha Ne unreliable
                "FIX_IMP": "F",  # impurity populations iterated like major species
                "FIX_T": "T",            # T pinned to the hydro profile for all runs (14-Oct-2026 policy)
                "FIX_T_AUTO": "F",       # No T auto-release
                "TAU_SCL_T": "0.0D0",    # No tau-scaling for T (FIX_T=T)
                "FIX_INB_T": "F",        # Don't fix T at inner boundary separately
                "DO_NG": "F",           # Disabled by default (13-Oct-2026): NG kick
                                        # overshoot pinned the X_SOB residual floor
                                        # near convergence; opt in via DO_NG=T.
                "CHK_NG": "T",          # Check NG acceleration validity
                "BEG_NG": "5.0D0",  # start Ng acceleration when max % change falls below this
                "IBEG_NG": "15",  # minimum iterations before first Ng acceleration
                "BW_NG": "40",  # depths per Ng acceleration band
                "ITS/NG": "8",  # run Ng every N-th iteration
                "DO_AV": "F",  # population averaging between iterations
                "METHOD": "ZERO",  # depth-correction interpolation (ZERO = constant between depths)
                "SCALE_OPT": "MAJOR",     # MAJOR: major species drive scaling
                "LAM_SCALE_OPT": "NONE", # NONE: disable LIMIT clipping — allow large decreases for LTE→NLTE transition
            },
        }
        if profile in stable_overrides:
            desired.update(stable_overrides[profile])

    # Solver selection: controlled by depth_independent_solver flag.
    #
    # N_TYPE is always N_ON_J for non-monotonic velocity models:
    #   H (flux) goes through zero at velocity reversals, so G_ONLY (N/H)
    #   divides by zero -> garbage Eddington factors -> population divergence.
    #   N_ON_J (N/J) is stable since J > 0 always (from formal solution).
    #   Required for OBS (cmf_flux.exe) which crashes with G_ONLY.
    #   Ref: CMFGEN manual (full_descr.txt): N_ON_J is "preferred option".
    #   N_ON_J remains well-defined where the flux H changes sign.
    #
    # SOL_METH is controlled by depth_independent_solver:
    #   True  -> DIAG (NUM_BNDS=1): each depth solved independently.
    #            Trace species instability stays local, no inter-depth propagation.
    #   False -> TRIDIAG (NUM_BNDS=3): couples adjacent depths in BA matrix.
    #            More physically accurate but unstable with non-monotonic velocity.
    #
    # With the default USE_LAM_ES=T, N_TYPE has no effect on the main solver
    # (RJ=TC in comp_j_blank.f); it only affects the Eddington factor
    # computation in FG_J_CMF (used by Sobolev Jbar).  With the opt-in
    # USE_LAM_ES=F path, MOM_J_CMF_V11 does run and N_TYPE selects its
    # boundary variable (N_ON_J stable where H=0; G_ONLY is diagnostic-only).
    desired["N_TYPE"] = "N_ON_J"
    if depth_independent_solver:
        num_bnds = 1
        sol_meth = "DIAG"
    else:
        num_bnds = 3
        sol_meth = "TRIDIAG"
    desired["SOL_METH"] = sol_meth
    n_type_val = desired["N_TYPE"]

    before = {k: None for k in desired}
    key_idx = find_existing_keys(lines)
    for k in desired:
        if k in key_idx:
            raw = _strip_inline_comment(lines[key_idx[k]]).strip()
            m = _VADAT_LINE_RE.match(raw)
            before[k] = m.group(2).strip() if m else None

    for k, v in desired.items():
        if user_overrides and k in user_overrides:
            continue
        set_or_append_keyword(lines, k, v)

    changes: List[str] = []
    for k, v in desired.items():
        if user_overrides and k in user_overrides:
            changes.append(f"{k}: preserved (user override = {user_overrides[k]})")
            continue
        old = before[k]
        if old is None:
            changes.append(f"{k}: (added) -> {v}")
        elif old != v:
            changes.append(f"{k}: {old} -> {v}")

    # T policy stamps by mode (explicit vadat_overrides always win):
    #   pinned (fix_temperature): all-depths hydro pin.
    #   thin (default): FIX_T=F + FIX_INB_T=T; FIX_X_DPTH is auto-filled
    #     later from the tau=2/3 photosphere, so T solves only where the
    #     ejecta is optically thin.
    #   free (solve_temperature): no stamps -- user keys rule verbatim
    #     (irradiated ejecta with a PNT_SRCE source, or manual domains).
    if thin_t_release:
        t_stamps = (("FIX_T", "F"), ("FIX_T_AUTO", "F"), ("FIX_INB_T", "T"))
    elif allow_t_release:
        t_stamps = ()
    else:
        t_stamps = (("FIX_T", "T"), ("FIX_T_AUTO", "F"))
    for k, v in t_stamps:
        if user_overrides and k in user_overrides:
            changes.append(f"{k}: preserved (user override = {user_overrides[k]})")
            continue
        set_or_append_keyword(lines, k, v)
        changes.append(f"{k}: -> {v} (T policy)")

    # Non-coherent electron scattering companions.  The effective COH_ES is
    # the user override when present, else the profile default.  With
    # COH_ES=F the ES source term uses RJ_ES (J convolved with the e.s.
    # kernel, ES_J_CONV), so the frequency grid must carry the extended
    # scattering wings and the linearization must not mix in a coherent
    # fraction.
    coh_es_eff = desired.get("COH_ES", "T")
    if user_overrides and "COH_ES" in user_overrides:
        coh_es_eff = str(user_overrides["COH_ES"])
    if coh_es_eff.strip().upper() == "F":
        es_changes = apply_noncoherent_es_companions(
            lines,
            explicit_user_keys=set(user_overrides) if user_overrides else set(),
        )
        changes.extend(f"[ES_NONCOH] {c}" for c in es_changes)

    if convergence_mode == "stable":
        changes.append("Convergence mode: stable (extra damping and longer lambda-phase enabled)")
    return changes, num_bnds


# ---------- OUT to IN (same as cur_cmf/com/out2in.sh) ----------
def out2in(wdir: Path) -> List[Tuple[str, str]]:
    """
    Copy *OUT files to *_IN in wdir (population levels from converged run).
    Matches cur_cmf/com/out2in.sh: each name ending with OUT becomes *_IN.
    Preserves previous *_IN files by backing them up to *_IN.bak.
    Uses copy instead of rename so *OUT files are also preserved
    (like cpmod behaviour).
    Returns list of (old_name, new_name) for reporting.
    """
    import shutil
    renames: List[Tuple[str, str]] = []
    for f in sorted(wdir.iterdir()):
        if not f.is_file():
            continue
        name = f.name
        if not name.endswith("OUT"):
            continue
        # Two CMFGEN naming conventions exist:
        #   HIOUT      -> HI_IN
        #   CFDAT_OUT  -> CFDAT_IN
        # Handle the explicit "_OUT" suffix first to avoid producing "CFDAT__IN".
        if name.endswith("_OUT"):
            in_name = name[:-4] + "_IN"
        else:
            in_name = name[:-3] + "_IN"
        dest = wdir / in_name
        try:
            # Backup existing *_IN file if it exists
            if dest.exists():
                backup = wdir / (in_name + ".bak")
                shutil.copy2(str(dest), str(backup))
            # Copy (not move) so *OUT is preserved
            shutil.copy2(str(f), str(dest))
            renames.append((f.name, in_name))
        except OSError:
            pass
    return renames


# ---------- NEW base selection + rename behavior ----------
def select_base_vadat_path(wdir: Path, profile: str, vadat_out: Path, vadat_init: Path) -> Path:
    """
    Select base VADAT file for the given profile:
      single_final  -> base is VADAT_INIT (LTE cold start; no prior VADAT needed)
      final         -> base is existing VADAT (warm restart from previous run)
    """
    if profile == "single_final":
        # prefer VADAT_INIT, else fall back to VADAT_INIT.bak
        if vadat_init.exists():
            return vadat_init
        vadat_init_bak = vadat_init.with_name(vadat_init.name + ".bak")
        if vadat_init_bak.exists():
            return vadat_init_bak
        raise FileNotFoundError(f"Missing {vadat_init.name} and {vadat_init_bak.name} in {wdir}")
    else:
        # final: use existing VADAT only (warm restart)
        if vadat_out.exists():
            return vadat_out
        raise FileNotFoundError(f"Missing {vadat_out.name} in {wdir} ({profile} profile requires it.)")

# def maybe_rename_vadat_inint_to_init_bak(wdir: Path, profile: str, vadat_init: Path) -> Optional[Path]:
#     """
#     After writing VADAT for init profile, rename VADAT_ININT -> VADAT_INIT.bak.

#     Safe behavior:
#     - If VADAT_ININT does not exist: do nothing.
#     - If VADAT_INIT.bak exists: back it up (timestamped) before overwrite.
#     """
#     if profile != "init":
#         return None

#     src = wdir / "VADAT_ININT"
#     if not src.exists():
#         return None

#     # dst = wdir / (vadat_init.name + ".bak")  # VADAT_INIT.bak
#     # if dst.exists():
#     #     backup_if_exists(dst, "pre_rename")

#     dst.write_text(src.read_text())
#     src.unlink()
#     return dst


def compute_optimal_vfrac_mom(hydro_data, safety_factor: float = 1.5) -> tuple[float, bool]:
    """
    Compute optimal VFRAC_MOM for the RT grid refinement in the moment solver.

    VFRAC_MOM (DELV_FRAC_MOM) is a RADIATIVE TRANSFER GRID RESOLUTION parameter,
    NOT a momentum correction or velocity damping parameter. It controls how
    many extra depth points the moment solver (MOM_J_CMF_V11 / MOM_JREL_V9 /
    PP_MOM_CMF_V1) inserts between the user's depth grid points for the
    ray-by-ray moment equation solve:

      T2 = VFRAC_MOM * min(VDOP_VEC)        [km/s]
      N_extra(I) = INT( |V(I) - V(I+1)| / T2 )    (rounded up to even)

    A smaller VFRAC_MOM means more RT grid points (better resolution of
    velocity structure along rays, more accurate J, but slower and more
    memory). A larger VFRAC_MOM means fewer RT grid points (faster but
    less accurate at velocity reversals).

    VFRAC_MOM=0 is INVALID (division by zero in the solver).
    VFRAC_MOM does NOT modify the physical velocity field.

    This function computes a VFRAC_MOM that prevents excessive ND expansion
    (which slows down the computation and can cause memory issues) while
    still providing adequate resolution at velocity reversals. The key
    trade-off is:
      - Too small: CMFGEN inserts many extra points → ND grows → slow, memory
      - Too large: CMFGEN inserts no extra points → J inaccurate at reversals

    Since VDOP_VEC varies during runtime (it depends on populations and T),
    we use a conservative estimate based on the maximum velocity jump.

    **IMPORTANT**: This function analyzes the *hydro* velocity profile. When VEL_LAW=12
    (RV_SN_MODEL_SNIIN) is used with RD_IN_R_GRID=F, CMFGEN creates its own tau_ES-based
    grid that has MUCH larger velocity jumps than the hydro grid (because the tau grid
    has large R spacing in low-density outer regions). The VFRAC_MOM computed here may
    be insufficient for VEL_LAW=12. For precise grid control, use --target-nd with
    VEL_LAW=11 (RD_IN_R_GRID=T), which reads the R grid from RDINR directly.

    Parameters:
    -----------
    hydro_data : HydroData
        Hydro data containing velocity profile (from parse_sn_hydro)
    safety_factor : float
        Safety factor to apply (default: 1.5 for 50% margin, higher due to VDOP variability)

    Returns:
    --------
    tuple[float, bool]: (optimal_vfrac_mom, recommend_no_relativistic)
    """
    # hydro.vel_cms is in cm/s, convert to km/s
    velocities = np.abs(hydro_data.vel_cms) / 1e5  # Convert to km/s, outer->inner order
    
    # Calculate velocity jumps
    dv = np.abs(np.diff(velocities))
    max_dv = dv.max() if len(dv) > 0 else 0.0
    
    # Check for non-monotonic velocity profile (sign of velocity changes)
    vel_sign_changes = np.sum(np.diff(np.sign(np.diff(velocities))) != 0)
    is_non_monotonic = vel_sign_changes > 0
    
    # Conservative approach: ensure VFRAC_MOM is large enough even for small VDOP values
    # CMFGEN typically uses VDOP in range 20-40 km/s, so we use VDOP=20 as worst case
    conservative_vdop = 20.0  # Conservative estimate of minimum VDOP CMFGEN might use
    required_vfrac_mom = (max_dv * safety_factor) / conservative_vdop
    
    # Apply reasonable bounds.
    # VFRAC_MOM is a RT grid resolution parameter (Doppler-width fraction).
    #   VFRAC_MOM=0 is INVALID (division by zero in MOM_J_CMF_V11:279, etc.)
    #   VFRAC_MOM < 0.5: very fine RT grid — accurate but slow, may cause ND explosion
    #   VFRAC_MOM ~ 1.0: standard for well-resolved velocity grids (e.g. VEL_LAW=11)
    #   VFRAC_MOM ~ 5-15: coarse RT grid — fast, adequate when velocity jumps are large
    #   VFRAC_MOM = 20000: effectively no RT grid refinement (lte_hydro example)
    # For non-monotonic velocity with sharp reversals (this project), we need
    # VFRAC_MOM large enough to prevent ND explosion but small enough to
    # resolve the velocity structure at the interaction region.
    min_vfrac_mom = 0.5   # Minimum: VFRAC_MOM=0 is invalid; 0.5 is very fine
    max_vfrac_mom = 5.0   # Cap: with USE_LAM_ES=T, moment solver is secondary;
                          # 5.0 balances speed vs J accuracy at velocity reversals

    optimal_vfrac_mom = max(min_vfrac_mom, min(required_vfrac_mom, max_vfrac_mom))
    
    # For non-monotonic profiles, USE_LAM_ES=F (standard moment solver with SIGMA clamp) bypasses
    # MOM_JREL_V9 entirely, avoiding sigma cap instabilities. No need to disable relativistic terms.
    recommend_no_relativistic = False
    
    print(f"[VFRAC_MOM Optimization] Velocity analysis:")
    print(f"  Max velocity jump: {max_dv:.1f} km/s")
    print(f"  Velocity profile: {'non-monotonic' if is_non_monotonic else 'monotonic'} ({vel_sign_changes} sign changes)")
    print(f"  Conservative VDOP estimate: {conservative_vdop:.1f} km/s (worst case)")
    print(f"  Safety factor: {safety_factor:.1f}")
    print(f"  Required VFRAC_MOM: {required_vfrac_mom:.2f}")
    print(f"  Optimal VFRAC_MOM: {optimal_vfrac_mom:.2f}")
    
    # Check if current VFRAC_MOM would cause issues
    current_vfrac_mom = 0.5  # Default value
    current_threshold = current_vfrac_mom * conservative_vdop
    violations = np.sum(dv > current_threshold)
    
    if violations > 0:
        print(f"  Current VFRAC_MOM={current_vfrac_mom} would cause {violations} ND expansion violations")
        print(f"  → Upgrading to VFRAC_MOM={optimal_vfrac_mom:.2f} to prevent expansion")
        print(f"  Threshold with VFRAC_MOM={optimal_vfrac_mom:.2f}: {optimal_vfrac_mom * conservative_vdop:.1f} km/s")
    else:
        print(f"  Current VFRAC_MOM={current_vfrac_mom} is sufficient (no violations)")
    
    if is_non_monotonic:
        print(f"  → Non-monotonic velocity detected: VFRAC_MOM controls RT grid resolution; "
              f"USE_LAM_ES controls solver choice (independent parameters)")
    
    return optimal_vfrac_mom, recommend_no_relativistic


# ---------- main ----------
def main() -> int:
    p = argparse.ArgumentParser(
        description="Update CMFGEN VADAT from SN_HYDRO_DATA and write VADAT (clean)."
    )
    p.add_argument("directory", nargs="?", default=".", help="Model directory (default: current).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--final", action="store_true", help="Final run (warm restart): base existing VADAT, apply final profile (LTE_EST=F, LIN_INT=T, reads *_IN populations).")
    g.add_argument("--single-final", action="store_true", help="Single-final LTE cold-start: base VADAT_INIT, write VADAT with single_final profile (LTE_EST=T, no prior run needed).")

    p.add_argument("--vadat-init", default="VADAT_INIT", help="Template VADAT init file base name (default: VADAT_INIT).")
    p.add_argument("--vadat-out", default="VADAT", help="Output VADAT filename (default: VADAT).")
    p.add_argument("--hydro", default="SN_HYDRO_DATA", help="Hydro file name (default: SN_HYDRO_DATA).")
    p.add_argument("--lum-xg", default="lum.xg", help="lum.xg path relative to directory for LSTAR (default: lum.xg); set empty to disable.")
    p.add_argument(
        "--no-lum-rst-geometry",
        action="store_true",
        help=(
            "When SN_HYDRO_DATA has Energy Lum [lsun], still use inner-grid RSTAR and "
            "RMAX=R_outer/R_inner (legacy). Default: RSTAR=R_phot(tau=2/3), RMAX=R_outer/R_phot."
        ),
    )
    p.add_argument("--no-write-rdinr", action="store_true", help="Do not write RDINR (default: write RDINR from SN_HYDRO_DATA so CMFGEN uses this grid).")
    p.add_argument("--no-tin-gammas", action="store_true", help="Do not write T_IN or GAMMAS_IN (default: write both when Temperature and Electron density are present in SN_HYDRO_DATA).")
    p.add_argument("--model-spec", default="MODEL_SPEC", help="MODEL_SPEC file to update so ND/NP match SN_HYDRO_DATA (default: MODEL_SPEC).")
    p.add_argument("--no-model-spec", action="store_true", help="Do not update MODEL_SPEC.")
    p.add_argument("--no-outgen", action="store_true", help="Do not append output to OUTGEN (default: append to OUTGEN).")
    p.add_argument("--no-out2in", action="store_true", help="Do not rename *OUT to *_IN for --final (default: do out2in so init populations become IN).")
    p.add_argument(
        "--convergence-mode",
        choices=("default", "stable"),
        default="stable",
        help=(
            "Convergence tuning: stable (default) applies stronger damping and longer Lambda phases "
            "for NLTE-heavy SN models; use default only if you need the lighter preset."
        ),
    )
    p.add_argument(
        "--fix_ne",
        action="store_true",
        help="DEPRECATED, ignored: FIX_NE is always F -- CMFGEN solves Ne from "
             "charge conservation. Accepted for backward compatibility only.",
    )
    p.add_argument(
        "--vadat-overrides",
        default=None,
        help="Path to YAML config file containing vadat_overrides section to apply after auto-configuration.",
    )
    p.add_argument(
        "--remove-species",
        default=None,
        help="Comma-separated ISF keys to comment out in MODEL_SPEC (e.g. 'FeIV_ISF,FeV_ISF'). "
             "Can also be read from 'remove_species' key in the YAML config's cmfgen section.",
    )
    p.add_argument(
        "--vel-law",
        type=int,
        default=12,
        help="CMFGEN VEL_LAW value (default: 12 for IIn tau-based grid; 11 for direct SN_HYDRO_DATA R/V/SIGMA). "
             "VEL_LAW=12 (RV_SN_MODEL_SNIIN) is designed for non-monotonic velocity flows: builds a tau-based "
             "grid from input_hydro.dat (RD_IN_R_GRID=F). VEL_LAW=11 reads R/V from RDINR (RD_IN_R_GRID=T). "
             "When --target-nd is specified, VEL_LAW is forced to 11 unless vel_law=12 is explicitly requested.",
    )
    p.add_argument(
        "--target-nd",
        type=int,
        default=None,
        help="Target number of depth points for the CMFGEN grid (default: use all hydro points). "
             "When set, downsamples the hydro R/V grid to target_nd points using velocity-aware spacing. "
             "With VEL_LAW=11: forces RD_IN_R_GRID=T so CMFGEN uses the downsampled RDINR grid directly. "
             "With VEL_LAW=12: target_nd only affects RDINR (which CMFGEN ignores; input_hydro.dat has full grid). "
             "Example: --vel-law 11 --target-nd 150",
    )
    p.add_argument(
        "--t-exc-floor",
        type=float,
        default=0.7,
        help="Minimum T_EXC (10^4 K) for the VADAT LTE excitation floor (default: 0.7 = 7000 K). "
             "For cold systems (T_gas < 3000 K) where the proxy departure coefficients already "
             "model the T_exc/T_gas gap, lower this (e.g. 0.3 = 3000 K) to prevent the VADAT "
             "T_EXC from over-exciting relative to the proxy's own excitation model.",
    )
    p.add_argument(
        "--rdinr-max-dv-kms",
        type=float,
        default=None,
        help="Refine the RDINR R grid so that no cell carries |dV| larger than this "
             "(km/s). Interior points are placed at equal-V steps inside each "
             "violating cell; CMFGEN recomputes V and SIGMA from the hydro spline "
             "at the new radii, so the physical profile is unchanged. Targets the "
             "comoving-frame discretization instability at strong compressions "
             "(per-cell Doppler shift >> frequency spacing). Only meaningful with "
             "VEL_LAW=11 / RD_IN_R_GRID=T.",
    )
    p.add_argument(
        "--rdinr-ib-ins",
        type=int,
        default=0,
        help="Split the innermost RDINR cell into this many geometric sub-cells "
             "(log-R widths, finest at the boundary). Mirrors the N_IB_INS scheme "
             "in SET_SN_R_GRID that RD_IN_R_GRID=T bypasses; resolves the steep "
             "inner tau rise feeding the diffusion boundary. Only meaningful "
             "with VEL_LAW=11 / RD_IN_R_GRID=T.",
    )
    p.add_argument(
        "--rdinr-ob-ins",
        type=int,
        default=0,
        help="Same geometric splitting applied to the outermost RDINR cell; "
             "resolves the free-streaming transition at the thin outer edge.",
    )
    p.add_argument(
        "--rdinr-bnd-ratio",
        type=float,
        default=2.5,
        help="Geometric width ratio for --rdinr-ib-ins/--rdinr-ob-ins sub-cells "
             "(default 2.5, close to SET_SN_R_GRID's 0.4/0.16/0.06 fractions).",
    )
    p.add_argument(
        "--use-hydro-grid",
        action="store_true",
        default=False,
        help="Write USE_HYDRO_GRID=T in input_hydro.dat so RV_SN_MODEL_SNIIN uses the "
             "hydro R/V grid directly, bypassing the log-uniform tau regridding. "
             "Requires NHYDRO==ND (no downsampling). Gives complete control of the "
             "input grid to the caller. Only effective with VEL_LAW=12.",
    )
    p.add_argument(
        "--smooth-velocity-jumps",
        action="store_true",
        default=False,
        help="Insert intermediate zones at velocity discontinuities (e.g. the 35x "
             "shell-ejecta interface jump) to enable Lambda-iteration convergence. "
             "Log-interpolated zones reduce per-step sigma while preserving shell structure.",
    )
    p.add_argument(
        "--smooth-max-dv-ratio",
        type=float,
        default=2.0,
        help="Velocity ratio threshold for jump detection (default: 2.0). "
             "If V_max/V_min between consecutive zones exceeds this, insert zones.",
    )
    p.add_argument(
        "--smooth-n-insert",
        type=int,
        default=8,
        help="Number of intermediate zones to insert at each detected velocity jump (default: 8).",
    )
    p.add_argument(
        "--smooth-rho-sigma",
        type=float,
        default=0.0,
        help="Gaussian smoothing sigma for density profile in log-space (default: 0 = off). "
             "Reduces the 400x density jump at the shell interface.",
    )
    p.add_argument(
        "--smooth-temp-sigma",
        type=float,
        default=0.0,
        help="Gaussian smoothing sigma for temperature profile in log-space (default: 0 = off).",
    )
    p.add_argument(
        "--blend-composition",
        action="store_true",
        default=False,
        help="Blend chemical composition across discontinuities by Gaussian-smoothing "
             "mass fractions. Removes the 100% chemical jumps (H: 0->0.7, Fe: 0.44->0) "
             "that prevent Lambda-iteration convergence.",
    )
    p.add_argument(
        "--blend-sigma",
        type=float,
        default=5.0,
        help="Gaussian smoothing width for composition blending (default: 5 zones).",
    )
    p.add_argument(
        "--blend-solar-weight",
        type=float,
        default=0.0,
        help="Weight for blending toward solar composition (0=pure smoothing, 1=full solar). "
             "Default: 0 (smooth existing profile only).",
    )
    p.add_argument(
        "--min-mass-fraction",
        type=float,
        default=0.0,
        help="Minimum mass fraction floor for all species at all depths (default: 0 = off). "
             "Ensures every species has enough atoms/levels to affect the J field. "
             "Example: --min-mass-fraction 1e-6",
    )
    p.add_argument(
        "--monotonize-velocity",
        action="store_true",
        default=False,
        help="Replace non-monotonic velocity regions with smooth monotonic interpolation. "
             "Removes negative-sigma (dv/dr < 0) regions that cause convergence problems.",
    )
    p.add_argument(
        "--chemical-n-insert",
        type=int,
        default=0,
        help="Number of zones to insert at each chemical composition discontinuity (default: 0 = off). "
             "Ensures adequate grid resolution where opacity changes rapidly due to chemical transitions.",
    )
    p.add_argument(
        "--chemical-jump-threshold",
        type=float,
        default=0.1,
        help="Total mass-fraction jump threshold for detecting chemical interfaces (default: 0.1).",
    )
    p.add_argument(
        "--reversal-n-insert",
        type=int,
        default=0,
        help="Number of zones to insert at each velocity reversal (negative sigma region). "
             "Default: 0 (no insertion). Set to e.g. 50 to resolve the non-monotonic velocity "
             "profile with extra grid points without modifying the velocity structure.",
    )
    p.add_argument(
        "--no-relativistic",
        action="store_true",
        help="Disable time evolution (USE_DJDT_RTE=F, DO_DDT=F, DO_CO_MOV_DDT=F). "
             "USE_LAM_ES=F (standard moment solver with SIGMA clamp) handles non-monotonic velocity.",
    )
    p.add_argument(
        "--model-type",
        default="static_snapshot",
        help="CMFGEN model type (all models are static snapshots). "
             "Determines relativistic solver configuration via get_rel_config().",
    )
    p.add_argument(
        "--inc-rad-decays",
        default="1",
        help="Include radioactive decay isotopes: '1' (default) or '0' (force F to avoid isotope consistency issues). "
             "Was previously read from CMFGEN_INC_RAD_DECAYS env var.",
    )

    args = p.parse_args()
    # Set model_type from CLI arg (was previously read from CMFGEN_MODEL_TYPE env var)
    global model_type
    model_type = args.model_type
    wdir = Path(args.directory).resolve()
    if args.single_final:
        profile = "single_final"
    else:
        profile = "final"

    # Load VADAT overrides from YAML config if provided
    vadat_overrides = {}
    cmfgen_cfg = {}
    allow_t_release = False
    thin_t_release = False
    depth_independent_solver = True  # default: DIAG (depth-independent)
    iteration_strategy = None  # IT_SPECIFIER content from YAML
    adjust_corrections = None  # ADJUST_CORRECTIONS content from YAML
    non_mono_line_treatment = 'BLANK'  # default: blanketing (Sobolev invalid for non-monotonic flow)
    if args.vadat_overrides:
        try:
            import yaml
            config_path = Path(args.vadat_overrides)
            if not config_path.exists():
                print(f"  ERROR: --vadat-overrides file not found: {config_path}", file=sys.stderr)
                print(f"  vadat_overrides (THK_CONT, OB_RAT, STOP_IF_BP, etc.) are critical for correct CMFGEN setup.", file=sys.stderr)
                print(f"  Refusing to continue with missing config — would produce wrong boundary conditions.", file=sys.stderr)
                sys.exit(1)
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
            if config and 'cmfgen' in config:
                cmfgen_cfg = config['cmfgen']
                # T policy (see header recipe):
                #   default            -> thin_t_release: FIX_T=F +
                #     FIX_INB_T=T + FIX_X_DPTH=<tau=2/3 index>; T solved
                #     only in the optically thin region.
                #   fix_temperature    -> all-depths hydro pin (the pre-
                #     default policy; keep for pinned A/B comparisons).
                #   solve_temperature  -> fully free T: no pins stamped,
                #     user keys in vadat_overrides rule verbatim.
                fix_temperature = bool(cmfgen_cfg.get('fix_temperature', False))
                solve_temperature = bool(cmfgen_cfg.get('solve_temperature', False))
                allow_t_release = not fix_temperature
                thin_t_release = allow_t_release and not solve_temperature
                if fix_temperature and solve_temperature:
                    print(
                        "  WARNING: cmfgen.fix_temperature and "
                        "solve_temperature both set -- fix_temperature wins "
                        "(all-depths pin)"
                    )
                if 'vadat_overrides' in cmfgen_cfg:
                    vadat_overrides = cmfgen_cfg['vadat_overrides']
                    print(f"  Loaded vadat_overrides from {config_path}: {list(vadat_overrides.keys())}")
                if 'depth_independent_solver' in cmfgen_cfg:
                    depth_independent_solver = bool(cmfgen_cfg['depth_independent_solver'])
                    print(f"  depth_independent_solver={depth_independent_solver} from config")
                # Iteration strategy (IT_SPECIFIER) and correction scaling
                # (ADJUST_CORRECTIONS) are optional convergence control files
                # read by CMFGEN every iteration. See docs/full_descr.txt
                # sections "IT SPECIFIER" and "ADJUST CORRECTIONS".
                iteration_strategy = cmfgen_cfg.get('iteration_strategy')
                adjust_corrections = cmfgen_cfg.get('adjust_corrections')
                _nms = str(vadat_overrides.get('NON_MONO_SOLVER', '')).strip().upper()
                if _nms.startswith('T') and iteration_strategy is not None:
                    print("  [NON_MONO] keeping iteration_strategy "
                          "(explicit IT_SPECIFIER controls Lambda→JFNK switch)")
                if 'remove_species' in cmfgen_cfg and not args.remove_species:
                    sp_list = cmfgen_cfg['remove_species']
                    if isinstance(sp_list, list):
                        args.remove_species = ','.join(sp_list)
                    else:
                        args.remove_species = str(sp_list)
                    print(f"  remove_species from config: {args.remove_species}")
                # Non-monotonic solver mode selection.
                # Primary key: non_mono_mode ('KNOP' or 'SOB')
                # Legacy key: non_mono_line_treatment ('BLANK' or 'SOB')
                # BLANK maps to KNOP mode, SOB maps to SOB mode.
                # Default: KNOP (Knop et al. 2009 formal solver + BLANK).
                non_mono_mode = cmfgen_cfg.get('non_mono_mode', None)
                if non_mono_mode is not None:
                    _nm = str(non_mono_mode).strip().upper()
                    if _nm not in ('KNOP', 'SOB'):
                        print(f"  WARNING: non_mono_mode="
                              f"'{non_mono_mode}' is not KNOP or SOB; "
                              f"defaulting to KNOP")
                        _nm = 'KNOP'
                    non_mono_mode = _nm
                    # Derive line_treatment from mode
                    non_mono_line_treatment = 'BLANK' if non_mono_mode == 'KNOP' else 'SOB'
                    print(f"  non_mono_mode={non_mono_mode} "
                          f"(line_treatment={non_mono_line_treatment})")
                else:
                    # Legacy: use non_mono_line_treatment directly
                    non_mono_line_treatment = cmfgen_cfg.get(
                        'non_mono_line_treatment', 'BLANK')
                    _lt = str(non_mono_line_treatment).strip().upper()
                    if _lt not in ('SOB', 'BLANK'):
                        print(f"  WARNING: non_mono_line_treatment="
                              f"'{non_mono_line_treatment}' is not SOB or BLANK; "
                              f"defaulting to BLANK")
                        _lt = 'BLANK'
                    non_mono_line_treatment = _lt
                    non_mono_mode = 'KNOP' if _lt == 'BLANK' else 'SOB'
                    print(f"  non_mono_line_treatment={non_mono_line_treatment} "
                          f"(mode={non_mono_mode})")
            else:
                print(f"  WARNING: No cmfgen section found in {config_path}")
        except SystemExit:
            raise
        except Exception as e:
            print(f"  ERROR: Failed to load config from {args.vadat_overrides}: {e}", file=sys.stderr)
            sys.exit(1)

    # Tee stdout to OUTGEN so prep output appears in the same log as CMFGEN
    outgen_path = wdir / "OUTGEN"
    outgen_file = None
    old_stdout = sys.stdout
    if not args.no_outgen:
        try:
            outgen_file = open(outgen_path, "a")
            sys.stdout = Tee(old_stdout, outgen_file)
        except OSError:
            outgen_file = None

    def _restore_stdout():
        if outgen_file is not None:
            sys.stdout = old_stdout
            outgen_file.close()

    try:
        print_header(f"prepare_VADAT.py [{profile}]")
        print(f"  started:   {datetime.now().isoformat()}")
        print(f"  directory: {wdir.resolve()}")

        hydro_path = wdir / args.hydro
        vadat_init_path = wdir / args.vadat_init
        vadat_out_path = wdir / args.vadat_out

        if not hydro_path.exists():
            print(f"ERROR: missing {hydro_path}")
            return 1

        try:
            base_vadat_path = select_base_vadat_path(wdir, profile, vadat_out_path, vadat_init_path)
        except FileNotFoundError as e:
            print(f"ERROR: {e}")
            return 1

        print_subheader("Input files")
        print(f"  hydro:     {hydro_path.name}")
        print(f"  base VADAT: {base_vadat_path.name} -> {vadat_out_path.name}")

        # For final (warm restart): rename *OUT -> *_IN (previous run populations become IN)
        # single_final starts from LTE so no *OUT files exist yet; skip out2in.
        out2in_renames: List[Tuple[str, str]] = []
        if profile == "final" and not args.no_out2in:
            out2in_renames = out2in(wdir)
            if out2in_renames:
                print(f"  OUT -> IN (population files from previous run):")
                for old, new in out2in_renames:
                    print(f"    {old} -> {new}")

        hydro = parse_sn_hydro(hydro_path)

        if (
            args.smooth_velocity_jumps
            or args.monotonize_velocity
            or args.reversal_n_insert > 0
            or args.chemical_n_insert > 0
            or args.blend_composition
        ):
            print_subheader("Hydro grid processing")

        # Smooth velocity discontinuities by inserting intermediate zones.
        # This must happen BEFORE any file writing (RDINR, T_IN, GAMMAS_IN,
        # FIN_CAL_GRID, input_hydro.dat) so all files share the same smoothed grid.
        # The interacting-SN grid has a 35x velocity jump at the shell-ejecta
        # interface that prevents Lambda-iteration convergence.
        if args.smooth_velocity_jumps:
            hydro = smooth_hydrodata_jumps(
                hydro,
                max_dv_ratio=args.smooth_max_dv_ratio,
                n_insert=args.smooth_n_insert,
                smooth_rho_sigma=args.smooth_rho_sigma,
                smooth_temp_sigma=args.smooth_temp_sigma,
            )

        # Monotonize velocity (remove negative-sigma regions).
        # Done BEFORE zone insertion so the velocity profile is smooth throughout.
        if args.monotonize_velocity:
            hydro = monotonize_velocity(hydro)

        # Insert extra zones at velocity reversals (negative sigma regions).
        # This PRESERVES the non-monotonic velocity profile — it just adds
        # more grid points across the reversal for better resolution.
        # Done BEFORE chemical insertion so reversals are detected on the
        # original velocity profile.
        if args.reversal_n_insert > 0:
            hydro = insert_zones_at_velocity_reversals(
                hydro, n_insert=args.reversal_n_insert,
            )

        # Insert extra zones at chemical composition interfaces.
        # Done BEFORE blending so jumps are sharp and fewer are detected.
        if args.chemical_n_insert > 0:
            hydro = insert_zones_at_chemical_interfaces(
                hydro, n_insert=args.chemical_n_insert,
                jump_threshold=args.chemical_jump_threshold,
            )

        # Blend chemical composition across discontinuities.
        # The interacting-SN model has 100% chemical jumps at the shell/ejecta
        # interface (H: 0->0.7, Fe: 0.44->0) that prevent Lambda convergence.
        # Done AFTER chemical insertion so the blended profile has good resolution.
        if args.blend_composition:
            hydro = blend_composition_to_solar(
                hydro, sigma=args.blend_sigma, solar_weight=args.blend_solar_weight,
                min_mass_fraction=args.min_mass_fraction,
            )

        print_subheader("Boundary optical depth diagnostics")

        # Validate inner boundary optical depth for diffusion approximation
        # CMFGEN manual: "Rosseland optical depth at the inner boundary > 50 (ideally close to 100)"
        #
        # Prefer the header's "Inner boundary Rosseland tau" when available — it's
        # computed from constant kappa_es=0.34 by convert_stella_swd_to_cmfgen.py,
        # matching CMFGEN's opacity. The STELLA NLTE kappa column can give
        # unphysically low opacities in the recombination zone (T < 5300 K).
        _tau_check = compute_tau_from_outer(hydro)
        _tau_at_ib_kappa = float(_tau_check[0])  # from kappa column
        _tau_at_ib_header = hydro.params.get("tau_inner_header")
        _tau_at_ib_chi_header = hydro.params.get("tau_inner_chi_header")
        if _tau_at_ib_header is not None and _tau_at_ib_header > _tau_at_ib_kappa:
            _tau_at_ib = float(_tau_at_ib_header)
            print(f"  [Tau] Using header Rosseland tau={_tau_at_ib:.1f} (kappa column gives { _tau_at_ib_kappa:.1f} — STELLA NLTE kappa unreliable in recombination zone)")
        else:
            _tau_at_ib = _tau_at_ib_kappa
        # Report chi-column tau from header if available (for consistency documentation)
        if _tau_at_ib_chi_header is not None:
            _tau_ratio = _tau_at_ib / max(_tau_at_ib_chi_header, 1e-10)
            print(f"  [Tau] Chi-column tau (from header): {_tau_at_ib_chi_header:.1f}  (ratio: {_tau_ratio:.1f}x)")
        if _tau_at_ib < 50.0:
            print(f"\n  *** NOTE: tau at inner boundary = {_tau_at_ib:.1f} < 50 ***")
            print(f"      DIFFUSION boundary with LSTAR — tau>50 recommended for LTE validity but not required.")
            print(f"      For LTE validity at the inner boundary, consider increasing cmfgen_tau_max.\n")

        # Outer-boundary diagnostic: classify the hydro tail (steady wind vs
        # shell edge) and report the optical depth it carries.  Determines
        # whether a wind extension could be hydro-derived (VEL_LAW=12 path)
        # or whether grid-edge RMAX is the only self-consistent choice.
        report_outer_boundary(diagnose_outer_boundary(hydro, tau_arr=_tau_check))

        print_subheader("RT grid resolution estimate (VFRAC_MOM)")

        # Automatic VFRAC_MOM optimization: compute the RT grid resolution
        # parameter that prevents excessive ND expansion while resolving
        # velocity reversals. VFRAC_MOM does NOT modify the velocity field;
        # it only controls how many extra RT grid points the moment solver
        # inserts between depth points (see compute_optimal_vfrac_mom docs).
        optimal_vfrac_mom, recommend_no_relativistic = compute_optimal_vfrac_mom(hydro, safety_factor=1.5)
        
        lines = read_lines(base_vadat_path)
        model_spec_path = wdir / args.model_spec

        # Photospheric RSTAR/RMAX from Energy Lum path.
        if not args.no_lum_rst_geometry:
            maybe_set_rst_from_energy_lum(hydro, tau_phot=TAU_PHOT)

        print_subheader("SN_HYDRO_DATA summary")
        print("Data read from SN_HYDRO_DATA (verify correct file/epoch):")
        print_hydro_summary(hydro, "read")

        # Update scalar keywords if present
        scalar_map = {
            "RMAX": hydro.params["RMAX"],
            "RSTAR": hydro.params["RSTAR"],
            "VINF": hydro.params["VINF"],
            "SN_AGE": hydro.params["time_days"],
            "RCUB_RHO": hydro.params.get("RCUB_RHO", 0.0),
            "RHO_CORE": hydro.params.get("RHO_CORE", 0.0),
            "VEL_CORE": hydro.params.get("VEL_CORE", 100.0),
            "VCORE": hydro.params.get("VEL_CORE", 100.0),   # km/s; required for VEL_LAW=12
            "BETA1": 1.0,                                    # homologous flow; required for VEL_LAW=12
        }
        # single_final: use T at mid-depth so outer (NLTE) half gets excitation floor, not full LTE.
        t_exc = recommended_t_exc_init(hydro, t_exc_floor=args.t_exc_floor) if profile == "single_final" else recommended_t_exc(hydro, t_exc_floor=args.t_exc_floor)
        if t_exc is not None:
            scalar_map["T_EXC"] = t_exc
        scalar_changes = update_keywords_in_lines(lines, scalar_map, numeric_format_like_old=True)

        # update/append mass keyword in MSUN
        mass_msun = float(hydro.params["MASS"])
        mass_token = f"{mass_msun:.8E}".replace("E", "D")
        set_or_append_keyword(lines, MASS_KEY, mass_token)
        mass_changes = [f"{MASS_KEY}: (set) -> {mass_token}   ! computed from int [4 pi r^2 rho] dr"]

        # abundance update behavior
        ave_mf = depth_averaged_mass_fractions(hydro.mass_fractions)
        abund_changes: List[str] = []
        if ave_mf:
            abund_changes = update_abundances_mass_fraction(lines, ave_mf, model_spec_path=model_spec_path)
        # SET_ABUND_CLUMP requires sum of mass fractions = 1; renormalize in case template + hydro != 1
        renorm_changes = renormalize_vadat_mass_fractions_to_unity(lines)
        abund_changes.extend(renorm_changes)

        n_iso_hydro = sn_hydro_header_num_isotopes(hydro_path)
        n_iso = max(int(hydro.params.get("num_isotopes", 0) or 0), n_iso_hydro)
        profile_changes, num_bnds = apply_profile(
            lines, profile, convergence_mode=args.convergence_mode, num_isotopes=n_iso,
            depth_independent_solver=depth_independent_solver,
            user_overrides=vadat_overrides if vadat_overrides else None,
            allow_t_release=allow_t_release,
            thin_t_release=thin_t_release,
        )
        
        print_subheader("VADAT keywords: RT grid and velocity law")

        # Set VFRAC_MOM after apply_profile (single source — not in profile dicts).
        #
        # VFRAC_MOM (DELV_FRAC_MOM) is a RADIATIVE TRANSFER GRID RESOLUTION parameter
        # for the moment-equation solver (MOM_J_CMF_V11 / MOM_JREL_V9 / PP_MOM_CMF_V1).
        # It controls the maximum velocity spacing (in Doppler-width fractions) between
        # adjacent depth points along rays. If |V(I)-V(I+1)| > VFRAC_MOM * min(VDOP_VEC),
        # the solver inserts extra RT grid points between depths I and I+1.
        #
        # VFRAC_MOM does NOT modify the physical velocity field. It only controls RT
        # grid refinement. VFRAC_MOM=0 is INVALID (division by zero).
        #
        # When --target-nd is specified with VEL_LAW=11, force VFRAC_MOM=1.0 (the velocity-aware
        # RDINR grid has dV << VDOP_MIN, so no ND expansion occurs).
        # When VEL_LAW=12 is used, the tau-based grid is NOT velocity-aware — it has large
        # velocity jumps at the shock front. Forcing VFRAC_MOM=1.0 causes CMFGEN to
        # insert extra points (ND: 150→218), which then produce J values exploding to -1e84
        # at the inserted depths. Use the optimized VFRAC_MOM instead to prevent ND expansion.
        # However, if the user explicitly set VFRAC_MOM in vadat_overrides, respect it.
        _override_vfrac_mom = vadat_overrides.get("VFRAC_MOM") if vadat_overrides else None
        if _override_vfrac_mom:
            vfrac_mom_token = str(_override_vfrac_mom)
            print(f"  [VFRAC_MOM] Using explicit override from YAML: {vfrac_mom_token}")
        elif args.target_nd is not None and args.vel_law != 12:
            vfrac_mom_token = "1.0D0"
        else:
            vfrac_mom_token = f"{optimal_vfrac_mom:.2f}D0"

        # Validation: VFRAC_MOM must be > 0 (VFRAC_MOM=0 causes division by zero
        # in MOM_J_CMF_V11:279, MOM_JREL_V9:368, PP_MOM_CMF_V1:259)
        _vfrac_mom_val = float(vfrac_mom_token.replace("D", "E").replace("d", "e"))
        if _vfrac_mom_val <= 0.0:
            print(f"  [VFRAC_MOM] ERROR: VFRAC_MOM={vfrac_mom_token} is invalid (must be > 0). "
                  f"Forcing to 0.5D0 (minimum safe value).")
            vfrac_mom_token = "0.5D0"
            _vfrac_mom_val = 0.5
        if _vfrac_mom_val < 0.5:
            print(f"  [VFRAC_MOM] WARNING: VFRAC_MOM={vfrac_mom_token} is very fine (may cause "
                  f"ND explosion). Consider a larger value for coarse grids.")
        if _vfrac_mom_val > 50.0:
            print(f"  [VFRAC_MOM] WARNING: VFRAC_MOM={vfrac_mom_token} is very coarse (J may be "
                  f"inaccurate at velocity reversals). Consider a smaller value for accuracy.")

        vfrac_mom_changes = update_keywords_in_lines(lines, {"VFRAC_MOM": vfrac_mom_token}, numeric_format_like_old=True)
        if vfrac_mom_changes:
            profile_changes.extend(vfrac_mom_changes)
            print(f"  [VFRAC_MOM] Set to {vfrac_mom_token} "
                  f"(RT grid resolution: max {optimal_vfrac_mom:.1f} Doppler widths between depths)")
        
        # For non-monotonic velocity profiles, recommend --no-relativistic option
        # Note: We cannot automatically apply this here as it requires command-line flag
        if recommend_no_relativistic:
            print(f"  [MOMENT_SOLVER] WARNING: Non-monotonic velocity detected")
            print(f"    VFRAC_MOM RT grid resolution optimized - USE_LAM_ES handles non-monotonic velocity")
            print(f"    If issues persist, consider: run_model_grid.py configs/X1.yaml --no-relativistic")
        
        # RD_SN_DATA matches isotope MF lines to ISO(:) from NUC_DECAY_DATA; that table is only
        # loaded when INC_RAD_DECAYS=T. Force T from the file header so a missed parse cannot
        # leave INC_RAD_DECAYS=F (which yields NUM_ISOTOPES=0 and "isotope data not recognized").
        # However, if --inc-rad-decays=0, force F to avoid isotope consistency issues.
        inc_rad_decays_val = args.inc_rad_decays
        if n_iso_hydro > 0 and inc_rad_decays_val == "1":
            set_or_append_keyword(lines, "INC_RAD_DECAYS", "T")
            profile_changes.append(
                f"INC_RAD_DECAYS=T (SN_HYDRO_DATA declares {n_iso_hydro} isotope block(s); "
                "required for RD_SN_DATA / NUC_DECAY_DATA ISO table)"
            )
        elif n_iso_hydro > 0 and inc_rad_decays_val == "0":
            set_or_append_keyword(lines, "INC_RAD_DECAYS", "F")
            profile_changes.append(
                f"INC_RAD_DECAYS=F (forced by --inc-rad-decays=0 to avoid isotope consistency issues)"
            )

        # Determine whether the non-monotonic solver is enabled (from YAML
        # vadat_overrides or already set in the VADAT lines). This controls
        # VFRAC_FG: the non-monotonic formal-solution solver
        # (FG_J_CMF_NONMONO_V1) needs fine ray resolution to handle velocity
        # reversals; a coarse VFRAC_FG (e.g. 20000) causes J to blow up.
        # Check for DO_NOMONO_KNOP, DO_NOMONO_SOB, or legacy NON_MONO_SOLVER.
        non_mono_solver = False
        do_nomono_knop = False
        if vadat_overrides:
            for _key in ("DO_NOMONO_KNOP", "DO_NOMONO_SOB", "NON_MONO_SOLVER"):
                _val = str(vadat_overrides.get(_key, "")).upper()
                if _val.startswith("T"):
                    non_mono_solver = True
                    if _key == "DO_NOMONO_KNOP":
                        do_nomono_knop = True
                    break
        if not non_mono_solver:
            for _l in lines:
                for _key in ("DO_NOMONO_KNOP", "DO_NOMONO_SOB", "NON_MONO_SOLVER"):
                    if f"[{_key}]" in _l:
                        _val = _l.split("[")[0].strip().upper()
                        if _val.startswith("T"):
                            non_mono_solver = True
                            if _key == "DO_NOMONO_KNOP":
                                do_nomono_knop = True
                        break
                if non_mono_solver:
                    break
        # Legacy compatibility: NON_MONO_SOLVER=T without DO_NOMONO_KNOP
        # maps to DO_NOMONO_KNOP=T in the Fortran code.
        if non_mono_solver and not do_nomono_knop:
            # Check if legacy NON_MONO_SOLVER is the one that's set
            _legacy = False
            if vadat_overrides:
                _val = str(vadat_overrides.get("NON_MONO_SOLVER", "")).upper()
                if _val.startswith("T") and not str(vadat_overrides.get("DO_NOMONO_KNOP", "")).upper().startswith("T") and not str(vadat_overrides.get("DO_NOMONO_SOB", "")).upper().startswith("T"):
                    _legacy = True
            if not _legacy:
                for _l in lines:
                    if "[NON_MONO_SOLVER]" in _l:
                        _val = _l.split("[")[0].strip().upper()
                        if _val.startswith("T"):
                            _has_knop = any(f"[DO_NOMONO_KNOP]" in _l2 and _l2.split("[")[0].strip().upper().startswith("T") for _l2 in lines)
                            _has_sob = any(f"[DO_NOMONO_SOB]" in _l2 and _l2.split("[")[0].strip().upper().startswith("T") for _l2 in lines)
                            if not _has_knop and not _has_sob:
                                _legacy = True
                        break
            if _legacy:
                do_nomono_knop = True

        # enforce_sn_iin_physics() stamps DO_NOMONO_KNOP=T (BLANK mode) or
        # DO_NOMONO_SOB=T (SOB mode) into every generated VADAT just before
        # it is written, so a non-monotonic solver mode is always active —
        # even when neither the YAML overrides nor the incoming template
        # requests one. The fine VFRAC_FG is therefore unconditional; the
        # detection above still determines do_nomono_knop for the legacy
        # NON_MONO_SOLVER mapping and for diagnostics.
        non_mono_solver = True

        # VFRAC_FG: maximum projected-velocity spacing along each ray in the
        # formal solution (FG_J_CMF_*/FG_J_CMF_NONMONO_V1), in units of the
        # minimum Doppler width. The solvers insert extra ray points wherever
        # |d(V*mu)| > VFRAC_FG*min(VDOP) — on BOTH signs of the velocity step
        # (ABS difference since 2026-09-18), so velocity reversals are refined
        # rather than skipped.
        #
        # Fine value (0.5) is required for ANY non-monotonic solver mode:
        #   - DO_NOMONO_KNOP -> FG_J_CMF_NONMONO_V1 (multi-root resonance solver)
        #   - DO_NOMONO_SOB  -> FG_J_CMF_V13 running on a non-monotonic V field
        # A coarse value (e.g. 20000, the manual's "omit" setting for monotonic
        # winds) leaves the reversal under-resolved and J can blow up at
        # non-monotonic depths — v7_bench VADATs show KNOP=T + VFRAC_FG=20000
        # reached production this way. The manual's 20000 no-refinement
        # setting (full_descr.txt:4114-4129) is only safe for a genuinely
        # standard-solver VADAT, which this pipeline no longer emits.
        # An explicit vadat_overrides VFRAC_FG always wins.
        _vfrac_fg_override = vadat_overrides.get("VFRAC_FG") if vadat_overrides else None
        if _vfrac_fg_override is not None:
            vfrac_fg_token = str(_vfrac_fg_override)
            print(f"  [VFRAC_FG] Using explicit override from YAML: {vfrac_fg_token}")
        else:
            vfrac_fg_token = "0.50D0" if non_mono_solver else "20000.0D0"
        set_or_append_keyword(lines, "VFRAC_FG", vfrac_fg_token)

        # Set VEL_LAW from CLI argument (default: 12 for IIn models with tau-based grid)
        # VEL_LAW=12 (RV_SN_MODEL_SNIIN) is designed for non-monotonic velocity flows
        # (interacting/IIn SNe): it reads input_hydro.dat (SNIIN format) and builds its own
        # tau_ES-based grid (RD_IN_R_GRID=F), ignoring the RDINR R grid.
        # VEL_LAW=11 (SET_RV_HYDRO_MODEL_V3) reads R/V directly from RDINR (RD_IN_R_GRID=T).
        # When --target-nd is specified AND vel_law != 12, force VEL_LAW=11 so the downsampled
        # RDINR grid is used directly. When vel_law=12 is explicitly requested, respect it
        # even with target_nd (target_nd only affects RDINR, which VEL_LAW=12 ignores).
        effective_vel_law = args.vel_law
        if args.target_nd is not None and effective_vel_law not in (12, 13, 14):
            effective_vel_law = 11
            set_or_append_keyword(lines, "VEL_LAW", "11")
            set_or_append_keyword(lines, "RD_IN_R_GRID", "T")
            # vfrac_mom_token already resolves to the explicit
            # vadat_overrides value when the user set one, else "1.0D0"
            # for this branch. Do not stomp the override here.
            set_or_append_keyword(lines, "VFRAC_MOM", vfrac_mom_token)
            profile_changes.append(
                f"VEL_LAW=11, RD_IN_R_GRID=T, VFRAC_MOM={vfrac_mom_token} (--target-nd={args.target_nd}: "
                f"VEL_LAW=11 reads R grid from RDINR; velocity-aware {args.target_nd}-point grid "
                f"has dV << VDOP_MIN, preventing CMFGEN ND adjustment)"
            )
            print(f"  [target-nd={args.target_nd}] VEL_LAW=11, RD_IN_R_GRID=T, VFRAC_MOM={vfrac_mom_token}")
        elif effective_vel_law == 11:
            # VEL_LAW=11 without target_nd: use the full SN_HYDRO_DATA grid as RDINR.
            # SET_RV_HYDRO_MODEL_V3 reads R from RDINR and interpolates V from SN_HYDRO_DATA.
            # Only R must be monotonic (R decreasing outer->inner); V can be non-monotonic.
            # This preserves the conversion grid exactly — CMFGEN does NOT regrid.
            set_or_append_keyword(lines, "VEL_LAW", "11")
            set_or_append_keyword(lines, "RD_IN_R_GRID", "T")
            # VFRAC_MOM: with the full hydro grid, velocity jumps are small (set by conversion).
            # Use the optimized value to prevent ND expansion if any jumps exceed VDOP_MIN.
            set_or_append_keyword(lines, "VFRAC_MOM", vfrac_mom_token)
            profile_changes.append(
                f"VEL_LAW=11, RD_IN_R_GRID=T (SET_RV_HYDRO_MODEL_V3: reads R grid from RDINR, "
                f"interpolates V from SN_HYDRO_DATA; non-monotonic V is OK — only R must be monotonic). "
                f"Full {hydro.params['ND']}-point grid used (no CMFGEN regridding)."
            )
            print(f"  [VEL_LAW=11] RD_IN_R_GRID=T, full grid ND={hydro.params['ND']}, "
                  f"VFRAC_MOM={vfrac_mom_token}, VFRAC_FG={vfrac_fg_token}")
        elif effective_vel_law == 12:
            # VEL_LAW=12: RV_SN_MODEL_SNIIN for non-monotonic flows (IIn/interacting SNe)
            set_or_append_keyword(lines, "VEL_LAW", "12")
            set_or_append_keyword(lines, "RD_IN_R_GRID", "F")
            # Accuracy options for VEL_LAW=12 tau-based grid with non-monotonic velocity:
            #   METHOD: ZERO for all modes (zero-order depth interpolation; manual
            #     p.83 also recommends ZERO for ionization-front convergence).
            #     Do NOT override here — stable_overrides set METHOD=ZERO.
            #   VFRAC_MOM: grid-adaptive value from optimal_vfrac_mom so the moment
            #     solver inserts extra RT grid points only where the velocity jump
            #     actually requires it (VFRAC_MOM * min(VDOP_VEC) threshold).
            #   VFRAC_FG: set above (shared with all VEL_LAW paths): 0.50 for
            #     any non-monotonic solver mode, 20000 otherwise, YAML
            #     vadat_overrides:VFRAC_FG always wins.
            #   TRAP_J=T: trapezoidal weights for J integration (already set in profiles)
            set_or_append_keyword(lines, "VFRAC_MOM", vfrac_mom_token)
            set_or_append_keyword(lines, "TRAP_J", "T")
            profile_changes.append(
                f"VEL_LAW=12, RD_IN_R_GRID=F (RV_SN_MODEL_SNIIN: tau-based grid from input_hydro.dat, "
                f"designed for non-monotonic velocity flows)"
            )
            _method_val = "ZERO"
            profile_changes.append(
                f"  Accuracy: METHOD={_method_val}, VFRAC_MOM={vfrac_mom_token} (grid-adaptive), "
                f"VFRAC_FG={vfrac_fg_token} (formal-solution path), TRAP_J=T"
            )
            if args.target_nd is not None:
                profile_changes.append(
                    f"  (target-nd={args.target_nd} only affects RDINR, which VEL_LAW=12 ignores; "
                    f"input_hydro.dat written with full SN_HYDRO_DATA grid)"
                )
                print(f"  [target-nd={args.target_nd}] VEL_LAW=12: target_nd affects RDINR only (ignored by CMFGEN)")
            print(f"  VEL_LAW=12, RD_IN_R_GRID=F (non-monotonic flow support)")
            _method_print = "ZERO"
            print(f"  Accuracy: METHOD={_method_print}, VFRAC_MOM={vfrac_mom_token}, VFRAC_FG={vfrac_fg_token}, TRAP_J=T")
        elif effective_vel_law in (13, 14):
            # VEL_LAW=13 (RV_SN_MODEL_SNHYD): same tau-uniform grid recipe as
            # VEL_LAW=12, but sourced directly from SN_HYDRO_DATA so no
            # input_hydro.dat is needed. V/SIGMA come from the monotone-cubic
            # interpolation used by SET_RV_HYDRO_MODEL_V3, so non-monotonic
            # sigma structure is preserved.
            # VEL_LAW=14 is identical except the grid coordinate integrates
            # the file's Kappa (Rosseland) column instead of a constant
            # kappa_es -- composition/ionization-aware, for H-poor models.
            tau_desc = "tau-ES" if effective_vel_law == 13 else "tau-Rosseland"
            set_or_append_keyword(lines, "VEL_LAW", str(effective_vel_law))
            set_or_append_keyword(lines, "RD_IN_R_GRID", "F")
            set_or_append_keyword(lines, "VFRAC_MOM", vfrac_mom_token)
            set_or_append_keyword(lines, "TRAP_J", "T")
            profile_changes.append(
                f"VEL_LAW={effective_vel_law}, RD_IN_R_GRID=F (RV_SN_MODEL_SNHYD: "
                f"{tau_desc} grid from SN_HYDRO_DATA; monotone-cubic V/SIGMA; "
                "no input_hydro.dat)"
            )
            print(
                f"  VEL_LAW={effective_vel_law}, RD_IN_R_GRID=F "
                f"({tau_desc} grid from SN_HYDRO_DATA)"
            )
        else:
            set_or_append_keyword(lines, "VEL_LAW", str(effective_vel_law))
            # VEL_LAW=11 reads R/V directly from RDINR (RD_IN_R_GRID=T)
            if effective_vel_law == 11:
                set_or_append_keyword(lines, "RD_IN_R_GRID", "T")
                profile_changes.append(f"VEL_LAW={effective_vel_law} (from --vel-law CLI argument), RD_IN_R_GRID=T")
            else:
                profile_changes.append(f"VEL_LAW={effective_vel_law} (from --vel-law CLI argument)")

        # Enforce INC_GRID=F when USE_LAM_ES=T (non-monotonic solver).
        # INC_GRID=T creates an enhanced spatial grid via EXTEND_VTSIGMA,
        # which interpolates sigma at ALL depths. With non-monotonic velocity
        # and USE_LAM_ES=T (J=TC formal solution), the enhanced grid can
        # introduce sigma artefacts at velocity reversals that destabilize
        # the radiative transfer solution.
        if non_mono_solver:
            _inc_grid_val = _get_key_token(lines, "INC_GRID")
            if _inc_grid_val is None or _inc_grid_val.upper().startswith("T"):
                set_or_append_keyword(lines, "INC_GRID", "F")
                if _inc_grid_val is None:
                    profile_changes.append("INC_GRID: (added) F (USE_LAM_ES=T requires INC_GRID=F)")
                else:
                    profile_changes.append(f"INC_GRID: {_inc_grid_val} -> F (USE_LAM_ES=T requires INC_GRID=F)")
            # Always print warning when USE_LAM_ES=T and INC_GRID=F
            print(f"  WARNING: USE_LAM_ES=T (non-monotonic solver) requires INC_GRID=F.")
            print(f"           INC_GRID=T with non-monotonic velocity creates sigma artefacts")
            print(f"           at velocity reversals that destabilize the RT solution.")
            print(f"           INC_GRID is {'F (already correct)' if _inc_grid_val and _inc_grid_val.upper().startswith('F') else 'F (forced)'}")

        # --no-relativistic: disable time-dependent solvers.
        # get_rel_config sets USE_LAM_ES=F + USE_J_REL=F + N_TYPE=N_ON_J
        # (J=TC formal solution, bypasses moment equations).
        # Only time-evolution keys are disabled; USE_LAM_ES is preserved.
        if args.no_relativistic:
            rel_keys = ["USE_DJDT_RTE", "DO_DDT", "DO_CO_MOV_DDT"]
            for key in rel_keys:
                set_or_append_keyword(lines, key, "F")
            profile_changes.append(
                "--no-relativistic: time-evolution solvers disabled; "
                "USE_LAM_ES=F (standard moment solver with SIGMA clamp) preserved for radiation transport."
            )

        die_key_changes = ensure_die_keys_for_model_spec(lines, model_spec_path)
        profile_changes.extend(die_key_changes)

        # Init/Stable/Final: ensure LAMBDA is enabled (DO_LAM_IT=T) and allow very large MAXCH so we don't abort.
        # DO_LAM_AUTO is conditional on solver type (controlled in IN_ITS):
        #   DIAG (NUM_BNDS=1): F — pure LAMBDA (block-diagonal BA adds no inter-depth info)
        #   TRIDIAG (NUM_BNDS=3): T — auto-switch to BA solve when MAXCH < LAM_VAL
        # MAX_LAM=1.0D10 for all profiles per user directive 15-Sep-2026 (was 1.5).
        # NOTE — superseded cascade history, watch for recurrence on X_1 restarts:
        # MAX_LAM=100 caused trace species explosion (100x/iter → 1e100 in 50 iters).
        # MAX_LAM=3 caused J-population cascade: trace species grew 3^20=3.5e9x in 20
        # iters, their opacity change destabilized J at other depths (iter 22 explosion
        # at depth 139, 6.4e19 → BA matrix singular at iter 119).
        # With MAX_LAM=1.5: LIT_LIM=-0.5 → 1.5x increase max; 1.5^20=3325x in 20 iters
        # — 1e6x slower than MAX_LAM=3, J has time to adjust smoothly.
        # Radiation solver: USE_LAM_ES=F (standard moment solver with SIGMA clamp) bypasses MOM_JREL_V9,
        #   which causes instabilities with negative sigma at velocity reversals.
        # MAX_CHNG: 1e50 to prevent convergence aborts due to large population changes
        # FIX_NE control:
        # - All profiles: FIX_NE=F (STELLA Saha Ne unreliable for X_1; predates Yat*AZ fix)
        # - FIX_IMP=F keeps trace-population corrections active in BA solves
        force_fix_ne = args.fix_ne

        if profile == "final" and out2in_renames:
            profile_changes.append(f"OUT->IN: {len(out2in_renames)} files (*OUT renamed to *_IN)")

        # Temperature control (both --single-final and --final): compute the
        # tau=2/3 photosphere index used by the thin-release default to pin T
        # at the inner depths (FIX_T=F, FIX_INB_T=T, FIX_X_DPTH=fix_x_dpth);
        # outer layers solve T.  See the header recipe for the three modes.
        nd = int(hydro.params["ND"])
        tau_arr, fix_x_dpth, i_phot = estimate_tau_es_and_fix_x_dpth(
            hydro, tau_phot=TAU_PHOT
        )
        tau_inner = float(tau_arr[0]) if tau_arr.size else 0.0
        tau_ob = float(tau_arr[-1]) if tau_arr.size else 0.0
        tau_phot_val = float(tau_arr[i_phot]) if tau_arr.size and 0 <= i_phot < tau_arr.size else 0.0

        # FIX_NE is always F (user directive 04-Oct-2026): CMFGEN solves the
        # electron density from charge conservation in every run.  The hydro
        # Ne is a Saha/LTE guess, not the NLTE ionization balance; holding it
        # fixed propagates that error into every ionization balance.  Any
        # FIX_NE=T / FIX_INB_NE=T request (--fix_ne flag or vadat_overrides)
        # is overridden here with a warning rather than silently ignored.
        _ne_t_requests = []
        if force_fix_ne:
            _ne_t_requests.append("--fix_ne")
        if str(vadat_overrides.get("FIX_NE", "")).strip().upper() == "T":
            _ne_t_requests.append("vadat_overrides FIX_NE=T")
        if str(vadat_overrides.get("FIX_INB_NE", "")).strip().upper() == "T":
            _ne_t_requests.append("vadat_overrides FIX_INB_NE=T")
        set_or_append_keyword(lines, "FIX_NE", "F")
        set_or_append_keyword(lines, "FIX_INB_NE", "F")
        if _ne_t_requests:
            _req = ", ".join(_ne_t_requests)
            profile_changes.append(
                f"FIX_NE=F, FIX_INB_NE=F (T requested via {_req} -- overridden; FIX_NE is always F)")
            print(f"  Ne control [{profile}]: FIX_NE=F "
                  f"(WARNING: {_req} overridden -- FIX_NE is always F)")
        else:
            profile_changes.append("FIX_NE=F, FIX_INB_NE=F (set in VADAT)")
            print(f"  Ne control [{profile}]: FIX_NE=F (set in VADAT)")

        print()

        print_subheader("Inner boundary (LSTAR / DIF / IB_METH)")

        # LSTAR: explicit cmfgen.inner_boundary_luminosity_lsun wins;
        # otherwise take it from SN_HYDRO_DATA (Energy Lum [lsun]).
        # The header value is the blackbody luminosity at the inner boundary
        # (zone behind tau_max=30), computed by convert_stella_swd_to_cmfgen.py.
        lstar_token: Optional[str] = None
        _lstar_cfg = cmfgen_cfg.get("inner_boundary_luminosity_lsun")
        if _lstar_cfg is not None:
            try:
                val = float(_lstar_cfg)
            except (TypeError, ValueError):
                print(
                    f"  ERROR: cmfgen.inner_boundary_luminosity_lsun="
                    f"{_lstar_cfg!r} is not numeric",
                    file=sys.stderr,
                )
                sys.exit(1)
            profile_changes = [c for c in profile_changes if not (c.strip().startswith("LSTAR:") and " -> " in c)]
            lstar_token = f"{val:.6E}".replace("E", "D")
            set_or_append_keyword(lines, "LSTAR", lstar_token)
            profile_changes.append(f"LSTAR: from cmfgen.inner_boundary_luminosity_lsun = {val:.6e} L_sun (explicit)")
        elif "energy_lum_lsun" in hydro.params:
            profile_changes = [c for c in profile_changes if not (c.strip().startswith("LSTAR:") and " -> " in c)]
            val = float(hydro.params["energy_lum_lsun"])
            lstar_token = f"{val:.6E}".replace("E", "D")
            set_or_append_keyword(lines, "LSTAR", lstar_token)
            profile_changes.append(f"LSTAR: from SN_HYDRO_DATA Energy Lum [lsun] = {val:.6e} L_sun (inner boundary above tau_max=30)")
        elif args.lum_xg:
            lum_xg_path = (wdir / args.lum_xg) if not Path(args.lum_xg).is_absolute() else Path(args.lum_xg)
            if lum_xg_path.exists():
                profile_changes = [c for c in profile_changes if not (c.strip().startswith("LSTAR:") and " -> " in c)]
                L_erg_s = read_lum_at_outer_from_lum_xg(lum_xg_path, float(hydro.params["time_days"]))
                if L_erg_s is not None and np.isfinite(L_erg_s) and L_erg_s > 0:
                    L_lsun = L_erg_s / LSUN_ERG_S
                    lstar_token = f"{L_lsun:.6E}".replace("E", "D")
                    set_or_append_keyword(lines, "LSTAR", lstar_token)
                    profile_changes.append(f"LSTAR: from {lum_xg_path.name} (outer) = {L_lsun:.6e} L_sun")

        # Inner boundary condition: respect config overrides for DIF/IB_METH.
        #
        # The diffusion approximation (DIF=T, IB_METH=DIFFUSION) computes
        # DBB = dB/dR ~ dT/dR at the inner boundary. When T is decreasing
        # inward (as in interacting SNe where the inner boundary is in a
        # cold dense shell below the T peak at the shell-ejecta interface),
        # DTDR < 0 makes DBB < 0, which makes I_CORE = (ETA + DBB*...)/CHI
        # negative, producing negative mean intensities and halting CMFGEN.
        #
        # ZERO_FLUX (DIF=F, IB_METH=ZERO_FLUX) sets I_CORE = I_M(NI,LS),
        # the local mean intensity — no DBB, no DTDR dependency. It is
        # only consistent for a genuinely dead inner boundary (LSTAR=0):
        # the boundary carries no flux, so a nonzero LSTAR is discarded
        # and the model converges with no energy input (X_SOB 40d: L=0
        # at both boundaries for all 1692 iterations).
        #
        # However, the converged X model used DIF=T/DIFFUSION with tau_inner=67.
        # If the config explicitly sets DIF/IB_METH, respect those overrides.
        # Otherwise default to ZERO_FLUX for safety with non-monotonic T.
        config_dif = vadat_overrides.get("DIF")
        config_ib_meth = vadat_overrides.get("IB_METH")
        if config_dif is not None or config_ib_meth is not None:
            # Config explicitly overrides — respect it
            ib_override = {}
            if config_dif is not None:
                ib_override["DIF"] = config_dif
            if config_ib_meth is not None:
                ib_override["IB_METH"] = config_ib_meth
            ib_changes = update_keywords_in_lines(lines, ib_override)
            if ib_changes:
                profile_changes.extend(ib_changes)
            dif_val = config_dif or "F"
            ib_val = config_ib_meth or "ZERO_FLUX"
            print(f"  [IB] Config override: DIF={dif_val}, IB_METH={ib_val} "
                  f"(tau_inner={_tau_at_ib:.1f}, LSTAR={lstar_token})")
            # Warn if DIF=T with low tau — diffusion approximation requires
            # Rosseland optical depth >> 1 (typically > 50). The tau_inner
            # from SN_HYDRO_DATA header is the electron scattering optical
            # depth, which can be much larger than the Rosseland optical
            # depth. With low Rosseland tau, DIF=T causes the grey
            # temperature solver to flatten dT/dR (DTDR → 0), which
            # eventually destabilizes the formal solution (J blows up).
            # SOB mode relaxes the cutoff: Dessart+2015 (arXiv:1503.05463)
            # truncates at Rosseland tau~30 with a diffusion boundary, so
            # only warn below ~20 there.
            _dif_warn_tau = 20.0 if non_mono_line_treatment == 'SOB' else 50.0
            if dif_val.upper().startswith("T") and _tau_at_ib < _dif_warn_tau:
                print(f"  WARNING: DIF=T with tau_inner={_tau_at_ib:.1f} < {_dif_warn_tau:.0f}.")
                print(f"           The diffusion approximation requires Rosseland optical depth >> 1.")
                print(f"           The SN_HYDRO_DATA tau is electron scattering; the actual Rosseland")
                print(f"           optical depth may be much lower (CMFGEN reports it in OUTGEN).")
                print(f"           If DTDR decreases toward zero and J blows up, switch to:")
                print(f"             DIF: F")
                print(f"             IB_METH: ZERO_FLUX")
                print(f"             IB_STAB: 0.1D0")
            # The opposite mismatch is now guarded, not just warned: an
            # explicitly flux-discarding boundary on an optically thick
            # inner edge leaves a luminous model with no energy source
            # (the m15 SOB run diverged on ZERO_FLUX at tau_es~37, job
            # 12001805 converged with DIFFUSION + header LSTAR; the X_SOB
            # 40d ZERO_FLUX run converged with L=0 and rendered ~100x
            # under-luminous).  A nonzero effective LSTAR refuses unless
            # cmfgen.allow_dead_inner_boundary: true.
            _lstar_eff = _parse_lstar_lsun(
                vadat_overrides.get("LSTAR") or lstar_token
            )
            _ib_fatal, _ib_msgs = check_inner_boundary_luminosity(
                dif_val,
                ib_val,
                _tau_at_ib,
                _lstar_eff,
                allow_dead=bool(cmfgen_cfg.get("allow_dead_inner_boundary", False)),
                line_treatment=non_mono_line_treatment,
            )
            if _ib_msgs:
                _sev = "ERROR" if _ib_fatal else "WARNING"
                print(f"  {_sev}: {_ib_msgs[0].strip()}")
                for _m in _ib_msgs[1:]:
                    print(_m)
            if _ib_fatal:
                sys.exit(1)
        elif lstar_token is not None:
            # Auto-decide DIF/IB_METH based on the inner boundary optical depth.
            # DIF=T/IB_METH=DIFFUSION is the physically correct choice when
            # tau_inner >= 50 (deep in the diffusion regime).  The previous
            # default of always using DIF=F/ZERO_FLUX was a workaround for
            # non-monotonic T causing negative DBB, but with the e-scattering
            # tau now correctly reported in the header (>= 50 for optically
            # thick boundaries), DIF=T is both physically correct and
            # numerically stable (the converged X model used DIF=T with
            # tau_inner=67).  DIF=F/ZERO_FLUX is reserved for genuinely
            # optically thin inner boundaries.
            #
            # SOB mode uses a lower cutoff (tau >= 25): Dessart+2015
            # (arXiv:1503.05463 sec 3.1) truncates at Rosseland tau~30 with
            # a diffusion boundary.  The tau>=50 rule left the m15 model
            # (tau_es~37) on ZERO_FLUX which diverged from it1; with
            # DIFFUSION + header LSTAR it converged (job 12001805).
            _dif_tau_thresh = 25.0 if non_mono_line_treatment == 'SOB' else 50.0
            if _tau_at_ib >= _dif_tau_thresh:
                ib_changes = update_keywords_in_lines(lines, {
                    "DIF": "T",
                    "IB_METH": "DIFFUSION",
                })
                if ib_changes:
                    profile_changes.extend(ib_changes)
                print(f"  [IB] DIFFUSION boundary (tau_inner={_tau_at_ib:.1f} >= {_dif_tau_thresh:.0f}, LSTAR={lstar_token})")
                print(f"       Inner boundary is optically thick — diffusion approximation is valid.")
            else:
                ib_changes = update_keywords_in_lines(lines, {
                    "DIF": "F",
                    "IB_METH": "ZERO_FLUX",
                })
                if ib_changes:
                    profile_changes.extend(ib_changes)
                print(f"  [IB] ZERO_FLUX boundary (tau_inner={_tau_at_ib:.1f} < {_dif_tau_thresh:.0f}, LSTAR={lstar_token})")
                print(f"       Inner boundary is optically thin — zero-flux BC (no diffusion approximation).")

        print_subheader("CMFGEN grid files (RDINR / input_hydro.dat)")

        # Build all CMFGEN inputs from SN_HYDRO_DATA (only input file needed)
        vel_law_tok = _get_key_token(lines, "VEL_LAW")
        vel_law_int: Optional[int] = None
        if vel_law_tok is not None:
            try:
                vel_law_int = int(float(vel_law_tok))
            except (ValueError, TypeError):
                vel_law_int = None
        if vel_law_int is None:
            print("WARNING: could not parse VEL_LAW in VADAT; skipping strict hydro-input consistency checks.")
        elif vel_law_int not in (11, 12, 13, 14):
            print(f"WARNING: VEL_LAW={vel_law_int} is uncommon for SN_HYDRO_DATA-driven runs (expected 11, 12, 13 or 14).")

        nd_written = int(hydro.params["ND"])  # default; overwritten by write_rdinr if called
        if not args.no_write_rdinr:
            rdinr_path = wdir / "RDINR"
            # For VEL_LAW=12: skip velocity-aware downsampling (RDINR is ignored,
            # RD_IN_R_GRID=F; input_hydro.dat is the actual grid source).
            rdinr_target_nd = args.target_nd if effective_vel_law != 12 else None
            # CLI wins; fall back to cmfgen.rdinr_*_ins in the YAML so a
            # pending/requeued job picks up config edits on regeneration.
            rdinr_ib = args.rdinr_ib_ins or int(cmfgen_cfg.get('rdinr_ib_ins', 0) or 0)
            rdinr_ob = args.rdinr_ob_ins or int(cmfgen_cfg.get('rdinr_ob_ins', 0) or 0)
            rdinr_ratio = float(cmfgen_cfg.get('rdinr_bnd_ratio', args.rdinr_bnd_ratio))
            nd_written, r_oi, v_oi = write_rdinr(
                hydro, rdinr_path, target_nd=rdinr_target_nd,
                max_dv_kms=args.rdinr_max_dv_kms,
                n_ib=rdinr_ib, n_ob=rdinr_ob,
                bnd_ratio=rdinr_ratio)
            profile_changes.append(f"Wrote {rdinr_path.name} (R grid from SN_HYDRO_DATA, ND={nd_written})")
            print(f"RDINR written: ND={nd_written}  R(1)={r_oi[0]:.6e} R(ND)={r_oi[-1]:.6e} (1e10 cm)  "
                  f"V(1)={v_oi[0]:.2f} V(ND)={v_oi[-1]:.2f} (km/s)")
            # Verify velocity jumps on the written grid
            dv = np.abs(np.diff(v_oi))
            max_dv = dv.max() if len(dv) > 0 else 0.0
            # VDOP_VEC = 12.85 * sqrt(TDOP/AMASS_DOP + (VTURB/12.85)^2)
            # With TDOP=4, AMASS_DOP=1e40, VTURB=20: VDOP_VEC = 20.0 km/s (VTURB-dominated)
            vdop_vec_est = 20.0  # VTURB-dominated (AMASS_DOP=1e40 per manual recommendation)
            vfrac_mom_est = 1.0 if args.target_nd is not None else optimal_vfrac_mom
            vdop_min_est = vfrac_mom_est * vdop_vec_est
            violations = int(np.sum(dv > vdop_min_est))
            print(f"  Grid velocity check: max_dV={max_dv:.2f} km/s, VDOP_MIN~{vdop_min_est:.2f} km/s, "
                  f"violations={violations}")
            if violations > 0:
                print(f"  WARNING: {violations} velocity jumps exceed VDOP_MIN — CMFGEN may insert extra points!")
            print()

        # input_hydro.dat in SNIIN format is required for VEL_LAW=12.
        # Written unconditionally (not gated on --no-write-rdinr) because
        # RV_SN_MODEL_SNIIN and SET_ABUND_CLUMP read it for R/V/RHO/T.
        if vel_law_int == 12:
            input_hydro_path = wdir / "input_hydro.dat"
            use_hydro_grid = args.use_hydro_grid
            # ── Non-monotonic solver: force USE_HYDRO_GRID=T ───────────
            # The tau-based regridding path in rv_sn_model_sniin.f uses
            # LIN_INTERP (linear interpolation) to regrid V onto a new
            # tau-uniform grid.  This smooths velocity extrema and destroys
            # the extreme sigma values (sigma < -1) that the non-monotonic
            # formal solver (FG_J_CMF_NONMONO_V1) is designed to handle.
            # By using the hydro grid directly (USE_HYDRO_GRID=T), SIGMA_TRUE
            # is computed from the original V(R) without interpolation,
            # preserving the full sigma dynamic range from SN_HYDRO_DATA.
            # Requires NHYDRO == ND (no downsampling).
            if non_mono_solver and not use_hydro_grid:
                if args.target_nd is None:
                    use_hydro_grid = True
                    print("  [NON_MONO] USE_HYDRO_GRID forced: F -> T "
                          "(preserve true sigma, avoid interpolation smoothing)")
                else:
                    print("  [NON_MONO] WARNING: cannot force USE_HYDRO_GRID=T "
                          "with --target-nd (requires NHYDRO==ND). "
                          "Sigma values will be smoothed by tau regridding.")
            if use_hydro_grid and args.target_nd is not None:
                print("  WARNING: --use-hydro-grid with --target-nd: USE_HYDRO_GRID requires NHYDRO==ND, "
                      "but target_nd downsamples. Disabling USE_HYDRO_GRID.")
                use_hydro_grid = False
            if use_hydro_grid:
                print("  USE_HYDRO_GRID=T: RV_SN_MODEL_SNIIN will use hydro R/V grid directly (no tau regridding)")
                # Force N_OB_INS=0 and N_IB_INS=0 in VADAT: USE_HYDRO_GRID
                # bypasses boundary insertion (R(1:ND)=R_HYDRO(1:ND) directly).
                # Non-zero N_OB_INS would cause ND != NHYDRO and the
                # USE_HYDRO_GRID path would abort.
                set_or_append_keyword(lines, "N_OB_INS", "0")
                set_or_append_keyword(lines, "N_IB_INS", "0")
                profile_changes.append("USE_HYDRO_GRID=T: N_OB_INS=0, N_IB_INS=0 (no boundary insertion)")

            # For VEL_LAW=12, rv_sn_model_sniin.f computes MND=ND-2*NBND_INS.
            # If input_hydro.dat has exactly MND points, the Fortran code uses the
            # Python AMR grid directly (preserving resolution at sigma<-1 regions).
            # Otherwise it creates a log-uniform tau grid that may extend beyond the
            # hydro data R range, causing LIN_INTERP failures in set_abund_clump.f.
            # So downsample to MND = target_nd - 2*NBND_INS when target_nd is provided.
            # Read NBND_INS (=N_OB_INS) from the VADAT that was already written by apply_profile.
            n_ob_ins_tok = _get_key_token(lines, "N_OB_INS")
            nbnd_ins_sniin = 3  # default in rv_sn_model_sniin.f
            if n_ob_ins_tok is not None:
                try:
                    nbnd_ins_sniin = int(float(n_ob_ins_tok))
                    # Honour an explicit 0: USE_HYDRO_GRID / stella_subset
                    # grids must not be silently shortened by 2*NBND_INS.
                except (ValueError, TypeError):
                    pass
            if args.target_nd is not None:
                n_hydro = len(hydro.radius_cm)
                if nbnd_ins_sniin < 1:
                    mnd = int(args.target_nd)
                else:
                    mnd = int(args.target_nd) - 2 * nbnd_ins_sniin
                if mnd > 10 and mnd < n_hydro:
                    # Use gradient-aware downsampling for non-monotonic velocity.
                    # Uniform-V spacing fails when R(V) is not single-valued
                    # (velocity plateaus, reverse shocks).
                    # Budget the velocity jump so downsampling cannot create a
                    # cell CMFGEN's comoving-frame discretisation can't
                    # integrate.  Computed locally rather than reusing
                    # vdop_min_est, which is only defined in the
                    # "not args.no_write_rdinr" branch above.
                    # VDOP_VEC is VTURB-dominated (AMASS_DOP=1e40), and this
                    # branch always has args.target_nd set, so VFRAC_MOM=1.
                    hydro_max_dv_kms = 1.0 * 20.0
                    r_oi, v_oi, indices = downsample_hydro_gradient_aware(
                        hydro, mnd, max_dv_kms=hydro_max_dv_kms)
                    # r_oi is in 1e10 cm units (outer->inner), v_oi is in km/s (outer->inner)
                    r_new = r_oi[::-1]  # inner->outer, 1e10 cm
                    v_new = v_oi[::-1]  # inner->outer, km/s
                    # Use actual hydro values at kept indices (not interpolated)
                    # indices are into the outer->inner array; hydro is inner->outer
                    hydro_idx_oi_to_inner = len(hydro.radius_cm) - 1 - indices
                    # rho/t/kappa from hydro are inner->outer; indices are outer->inner
                    # so hydro_idx_oi_to_inner gives values in outer->inner order.
                    # r_new/v_new are reversed to inner->outer, so reverse rho/t/kappa too.
                    rho_new = hydro.rho_cgs[hydro_idx_oi_to_inner][::-1]  # inner->outer
                    t_new = hydro.temp_104k[hydro_idx_oi_to_inner][::-1]  # inner->outer
                    if hydro.kappa_cgs is not None:
                        kappa_new = hydro.kappa_cgs[hydro_idx_oi_to_inner][::-1]  # inner->outer
                    else:
                        kappa_new = None

                    # Orientation guard: verify r_new is inner->outer (increasing)
                    # and rho/t match the hydro data at the same radii.
                    assert np.all(np.diff(r_new) > 0), \
                        f"r_new must be increasing (inner->outer), got dr[0]={np.diff(r_new)[0]}"
                    # Spot-check: first point (innermost) should match hydro[0]
                    assert np.isclose(r_new[0], hydro.radius_cm[0] / 1e10, rtol=1e-10), \
                        f"r_new[0]={r_new[0]} != hydro R[0]={hydro.radius_cm[0]/1e10}"
                    assert np.isclose(rho_new[0], hydro.rho_cgs[0], rtol=1e-10), \
                        f"rho_new[0]={rho_new[0]} != hydro rho[0]={hydro.rho_cgs[0]}"
                    assert np.isclose(t_new[0], hydro.temp_104k[0], rtol=1e-10), \
                        f"t_new[0]={t_new[0]} != hydro T[0]={hydro.temp_104k[0]}"
                    n = mnd
                    r_1e10 = r_new  # already in 1e10 cm
                    v_kms = np.abs(v_new)  # already in km/s
                    t_k = t_new * 1e4
                    if kappa_new is not None:
                        kappa = kappa_new
                    else:
                        kappa = np.full(n, KAPPA_ES_CGS)
                    # Use constant kappa for tau (same as write_input_hydro_dat_sniin)
                    kappa = np.full(n, KAPPA_ES_CGS)
                    tau_oi = np.zeros(n)
                    for j in range(n - 2, -1, -1):
                        dr = (r_new[j + 1] - r_new[j]) * 1e10  # 1e10 cm -> cm
                        krho_avg = 0.5 * (kappa[j] * rho_new[j] + kappa[j + 1] * rho_new[j + 1])
                        tau_oi[j] = tau_oi[j + 1] + krho_avg * dr
                    if n >= 2 and tau_oi[-2] > 0:
                        tau_oi[-1] = max(tau_oi[-1], tau_oi[-2] * 1.0e-2)
                    else:
                        tau_oi[-1] = max(tau_oi[-1], 1.0e-30)
                    print(f"[input_hydro.dat] NHYDRO={n} (downsampled to MND={mnd}, target_nd={args.target_nd})")
                    print(f"[input_hydro.dat] R(1e10 cm): min={r_1e10.min():.6e} max={r_1e10.max():.6e}")
                    with open(input_hydro_path, "w", newline="\n") as f:
                        f.write(f" {'T' if use_hydro_grid else 'F'} F\n")  # USE_HYDRO_GRID ADD_ACC_ZONE
                        f.write("R(1e10cm) V(km/s) RHO T(K) XL XTAU TAU\n")
                        f.write(f"{n}\n")
                        f.write("! bogus line\n")
                        for j in range(n - 1, -1, -1):
                            xl = 0.0
                            f.write(f"  {r_1e10[j]:.12E}  {v_kms[j]:.12E}  {rho_new[j]:.12E}  {t_k[j]:.12E}  {xl:.12E}  {xl:.12E}  {tau_oi[j]:.12E}\n")
                    profile_changes.append(f"Wrote {input_hydro_path.name} (VEL_LAW=12, MND={mnd} for NHYDRO==MND path)")
                else:
                    write_input_hydro_dat_sniin(
                        hydro, input_hydro_path, use_hydro_grid=use_hydro_grid,
                        smooth_jumps=args.smooth_velocity_jumps,
                        smooth_max_dv_ratio=args.smooth_max_dv_ratio,
                        smooth_n_insert=args.smooth_n_insert,
                    )
                    profile_changes.append(f"Wrote {input_hydro_path.name} (VEL_LAW=12 RV_SN_MODEL_SNIIN)")
            else:
                write_input_hydro_dat_sniin(
                    hydro, input_hydro_path, use_hydro_grid=use_hydro_grid,
                    smooth_jumps=args.smooth_velocity_jumps,
                    smooth_max_dv_ratio=args.smooth_max_dv_ratio,
                    smooth_n_insert=args.smooth_n_insert,
                )
                profile_changes.append(f"Wrote {input_hydro_path.name} (VEL_LAW=12 RV_SN_MODEL_SNIIN)")
        else:
            profile_changes.append(
                f"Skipped input_hydro.dat write (VEL_LAW={vel_law_int}; SNIIN format is for VEL_LAW=12)"
            )

        print_subheader("Initial-condition files (T_IN / GAMMAS_IN / FIN_CAL_GRID)")

        # T_IN, GAMMAS_IN, and FIN_CAL_GRID for LTE initial conditions
        # (VEL_LAW=12 + LTE_EST=T + SN_HYDRO_MODEL=.FALSE. requires FIN_CAL_GRID for T/ED)
        if args.no_tin_gammas:
            print("  (T_IN/GAMMAS_IN/FIN_CAL_GRID skipped: --no-tin-gammas requested)")
        else:
            if hydro.temp_104k is not None and hydro.ed_cgs is not None:
                # Use actual max T from hydro (no artificial 10000 K cap);
                # the cap was causing bad LTE populations at inner depths where T > 10000 K.
                hydro_t_max_104 = float(np.max(hydro.temp_104k))
                hydro_t_min_104 = float(np.min(hydro.temp_104k))
                write_t_in_from_hydro(hydro, wdir / "T_IN", t_max_104=hydro_t_max_104, t_min_104=hydro_t_min_104)
                write_gammas_in_from_hydro(hydro, wdir / "GAMMAS_IN")
                write_fin_cal_grid_from_hydro(hydro, wdir / "FIN_CAL_GRID")
                profile_changes.append("Wrote T_IN, GAMMAS_IN, and FIN_CAL_GRID (LTE initial T/ED from SN_HYDRO_DATA)")
            else:
                print("  (T_IN/GAMMAS_IN skipped: SN_HYDRO_DATA needs Temperature and Electron density sections)")
        print()

        print_subheader("Consistency enforcement and vadat_overrides")

        # Enforce DC_METH=R when temperature is non-monotonic.
        # REGRID_LOG_DC_V1 with option TR interpolates in T space, but has a
        # known issue with non-monotonic T (regrid_log_dc_v1.f:284:
        # "This code may have issues if T is very nonmonotonic"). With multiple
        # T matches per depth, TR picks the wrong interval, producing pathological
        # DCs. R (radius interpolation) is safe because R is always monotonic.
        # hydro.temp_104k is in inner->outer order; REGRID reads _IN files in
        # outer->inner (depth 1 = outer), so we check non-monotonicity in both
        # orderings to catch the cannon/hot-shell T reversal.
        if hydro is not None and hydro.temp_104k is not None and hydro.temp_104k.size > 1:
            _t = hydro.temp_104k
            # Check in inner->outer order (T increasing outward)
            _n_non_mono_io = int(np.sum(np.diff(_t) > 0))
            # Check in outer->inner order (T increasing inward)
            _n_non_mono_oi = int(np.sum(np.diff(_t[::-1]) > 0))
            _n_non_mono = max(_n_non_mono_io, _n_non_mono_oi)
            if _n_non_mono > len(_t) * 0.1:
                print(f"  [DC_METH] T is non-monotonic (inner->outer: {_n_non_mono_io}, "
                      f"outer->inner: {_n_non_mono_oi}/{len(_t)-1} zones), "
                      f"enforcing DC_METH=R to avoid REGRID_LOG_DC_V1 TR interpolation failure")
                # Force DC_METH=R in the VADAT lines directly
                for _i, _ln in enumerate(lines):
                    if "[DC_METH]" in _ln:
                        lines[_i] = "R [DC_METH]"
                        break
                # Also set it in vadat_overrides so it survives any later override
                if vadat_overrides is not None:
                    vadat_overrides["DC_METH"] = "R"
                else:
                    vadat_overrides = {"DC_METH": "R"}

        # Apply vadat_overrides as final step (after all profile-specific settings
        # including --no-relativistic) so config overrides take precedence.
        if vadat_overrides:
            # SOL_METH is set by apply_profile (TRI for TRIDIAG); reject overrides
            # that would conflict with the NUM_BNDS/SOL_METH pairing.
            if "SOL_METH" in vadat_overrides:
                requested_sol_meth = str(vadat_overrides["SOL_METH"]).strip().upper()
                if requested_sol_meth not in ("DIAG", "TRI", "TRIDIAG"):
                    print(
                        f"  WARNING: vadat_overrides SOL_METH={requested_sol_meth} rejected. "
                        f"Only DIAG and TRI are supported."
                    )
                    vadat_overrides = {k: v for k, v in vadat_overrides.items() if k != "SOL_METH"}
            # Radiation solver keys are set by get_rel_config() — reject overrides
            # that would conflict with the USE_LAM_ES=F / N_ON_J / USE_J_REL=F configuration.
            # USE_LAM_ES is NOT protected — the converged X model used USE_LAM_ES=T,
            # and config overrides should be respected.
            _protected_solver_keys = {
                "INCL_REL",
                "USE_FRM_REL", "INCL_ADV_TRANS",
            }
            _conflicting = _protected_solver_keys & set(vadat_overrides.keys())
            if _conflicting:
                print(
                    f"  WARNING: vadat_overrides contains radiation solver keys "
                    f"({', '.join(sorted(_conflicting))}) — rejected. "
                    f"These are set by get_rel_config() (USE_LAM_ES=F, N_ON_J, USE_J_REL=F)."
                )
                vadat_overrides = {k: v for k, v in vadat_overrides.items()
                                   if k not in _conflicting}
            # FIX_NE / FIX_INB_NE are always F (Ne solved from charge
            # conservation): the Ne-control block above already stamped F and
            # warned on a T request.  Strip the keys here too so this late
            # re-application cannot undo the stamp.
            _ne_strip = {"FIX_NE", "FIX_INB_NE"} & set(vadat_overrides.keys())
            if _ne_strip:
                _ne_bad = sorted(k for k in _ne_strip
                                 if str(vadat_overrides.get(k, "")).strip().upper() != "F")
                if _ne_bad:
                    print(
                        f"  WARNING: vadat_overrides {', '.join(_ne_bad)} "
                        f"rejected -- FIX_NE is always F."
                    )
                vadat_overrides = {k: v for k, v in vadat_overrides.items()
                                   if k not in _ne_strip}
            # T is always fixed to the SN_HYDRO_DATA profile: these input
            # models rely on the hydrodynamic temperature, and a T release
            # (FIX_T=F, FIX_T_AUTO/DO_T_AUTO/DO_GT_AUTO=T) would replace it
            # with a radiative-equilibrium T.  Strip any request for a
            # non-fixed value so this late re-application cannot undo the
            # apply_profile stamp.  Skipped when cmfgen.solve_temperature
            # opts the model out of the pin (irradiated-ejecta or
            # FIX_INB_T/FIX_X_DPTH partial-depth release runs).
            if not allow_t_release:
                _t_pinned = {"FIX_T": "T", "FIX_T_AUTO": "F",
                             "DO_T_AUTO": "F", "DO_GT_AUTO": "F"}
                _t_strip = sorted(k for k in _t_pinned
                                  if k in vadat_overrides
                                  and str(vadat_overrides[k]).strip().upper()
                                      != _t_pinned[k])
                if _t_strip:
                    print(
                        f"  WARNING: vadat_overrides {', '.join(_t_strip)} "
                        f"rejected -- T is pinned to the hydro profile "
                        f"(FIX_T=T, FIX_T_AUTO/DO_T_AUTO/DO_GT_AUTO=F)."
                    )
                    vadat_overrides = {k: v for k, v in vadat_overrides.items()
                                       if k not in _t_strip}
            # Thin-region T release: when FIX_INB_T=T is in effect (stamped
            # by the thin default or requested via override) and the user did
            # not give FIX_X_DPTH (or left it at 0), fill it with the tau=2/3
            # photosphere index computed above -- pins T at the inner depths
            # only, exactly the "solve T where optically thin" domain.
            _fxinb_eff = str(
                vadat_overrides.get(
                    "FIX_INB_T", "T" if thin_t_release else "F"
                )
            ).strip().upper()
            if allow_t_release and _fxinb_eff == "T":
                _xdp = str(vadat_overrides.get("FIX_X_DPTH", "")).strip()
                if _xdp in ("", "0", "0.0", "0.0D0"):
                    vadat_overrides["FIX_X_DPTH"] = str(fix_x_dpth)
                    print(
                        f"  FIX_X_DPTH auto = {fix_x_dpth} "
                        f"(tau=2/3 photosphere at hydro index {i_phot}; "
                        f"override FIX_X_DPTH in vadat_overrides to pin a "
                        f"different number of inner depths)")
                _fixt_ov = vadat_overrides.get("FIX_T")
                _fixt_cur = str(_fixt_ov if _fixt_ov is not None
                                else _get_key_token(lines, "FIX_T") or "")
                if _fixt_cur.strip().upper() == "T":
                    print(
                        "  WARNING: FIX_INB_T=T has no effect while "
                        "FIX_T=T (the all-depths pin wins).  Set FIX_T=F in "
                        "vadat_overrides for a partial-depth T release.")
                _jfnk_eff = str(vadat_overrides.get(
                    "DO_JFNK", _get_key_token(lines, "DO_JFNK") or "F"
                )).strip().upper()
                if _jfnk_eff == "T":
                    print(
                        "  WARNING: DO_JFNK=T with a FIX_X_DPTH inner T pin: "
                        "JFNK mirrors only the MOD_FIX_T domain and leaves "
                        "the inner T rows live.  Set DO_JFNK=F.")
            final_overrides = update_keywords_in_lines(lines, vadat_overrides, numeric_format_like_old=True)
            # Also append any override keys that don't exist in the VADAT template
            existing = find_existing_keys(lines)
            for k, v in vadat_overrides.items():
                if k not in existing:
                    set_or_append_keyword(lines, k, str(v))
                    final_overrides.append(f"{k}: (added) {v}")

            # Consistency enforcement: if USE_J_REL was overridden to F, also disable
            # INCL_REL and USE_FRM_REL (CMFGEN check_param_consistency.f requires
            # USE_J_REL=T if INCL_REL=T, otherwise CMFGEN exits with
            # "Can only include relativistic terms if USE_J_REL is set to TRUE").
            use_j_rel_val = _get_key_token(lines, "USE_J_REL")
            if use_j_rel_val is not None and use_j_rel_val.strip().upper() in ("F", "FALSE"):
                auto_fixes = {}
                incl_rel_val = _get_key_token(lines, "INCL_REL")
                if incl_rel_val is not None and incl_rel_val.strip().upper() not in ("F", "FALSE"):
                    auto_fixes["INCL_REL"] = "F"
                frm_rel_val = _get_key_token(lines, "USE_FRM_REL")
                if frm_rel_val is not None and frm_rel_val.strip().upper() not in ("F", "FALSE"):
                    auto_fixes["USE_FRM_REL"] = "F"
                if auto_fixes:
                    fix_changes = update_keywords_in_lines(lines, auto_fixes, numeric_format_like_old=True)
                    if fix_changes:
                        print(f"  [Consistency] USE_J_REL=F: auto-disabled {', '.join(auto_fixes.keys())} "
                              f"(CMFGEN requires USE_J_REL=T for relativistic terms)")
                        final_overrides.extend(fix_changes)

            if final_overrides:
                print(f"  Applied vadat_overrides (final, post-profile):")
                for change in final_overrides:
                    print(f"    {change}")
                profile_changes.extend(final_overrides)

            # Align JFNK_FIXT with FIX_T so T is a Newton unknown iff FIX_T=F.
            # CMFGEN also enforces this at VADAT read time; keep the written
            # VADAT consistent for provenance.
            fix_t_tok = _get_key_token(lines, "FIX_T")
            if fix_t_tok is not None:
                fix_t_val = fix_t_tok.strip().upper()
                if fix_t_val in ("T", "TRUE", "F", "FALSE"):
                    jfnk_fixt_want = "T" if fix_t_val in ("T", "TRUE") else "F"
                    jfnk_fixt_now = _get_key_token(lines, "JFNK_FIXT")
                    if jfnk_fixt_now is None or jfnk_fixt_now.strip().upper() != jfnk_fixt_want:
                        set_or_append_keyword(lines, "JFNK_FIXT", jfnk_fixt_want)
                        msg = f"JFNK_FIXT={jfnk_fixt_want} (aligned to FIX_T={jfnk_fixt_want})"
                        print(f"  [Consistency] {msg}")
                        profile_changes.append(msg)

        # NUM_BNDS/SOL_METH set by apply_profile (N_ON_J always; DIAG or
        # TRIDIAG per depth_independent_solver).  TRIDIAG is the CMFGEN-
        # recommended setting (manual: "best compromise between convergence
        # and memory requirements").  DIAG is slower but simpler.
        # KNOWN ISSUE: VAR_MOM_J_CMF_V12 uses clamped SIGMA with non-monotonic
        # Eddington factors — this can produce zero BA corrections in strongly
        # reversing flows regardless of DIAG/TRIDIAG (see BA_AUDIT_REPORT.txt).
        # The proper fix is to use consistent Eddington factors for the BA
        # variation, not to switch from TRIDIAG to DIAG.
        if True:
            n_type_final = _get_key_token(lines, "N_TYPE") or "N_ON_J"
            sol_meth_final = _get_key_token(lines, "SOL_METH") or ("TRI" if num_bnds == 3 else "DIAG")
            override_msg = (
                f"NUM_BNDS={num_bnds}, SOL_METH={sol_meth_final} "
                f"(N_TYPE={n_type_final} for Eddington factor)"
            )
            print(f"  {override_msg}")
            profile_changes.append(override_msg)

        # # Always back up existing output VADAT before overwriting
        # bak = backup_if_exists(vadat_out_path, f"pre_{profile}")
        # if bak:
        #     print(f"Backed up existing {vadat_out_path.name} -> {bak.name}")

        print_subheader("Physics enforcement and VADAT write")

        # ======================================================================
        # RE-APPLY NON-MONOTONIC FORCED SETTINGS (post-override)
        # The non-monotonic forced settings (MAX_CHNG, LAM_SCALE_OPT, etc.)
        # are applied earlier, but user vadat_overrides can overwrite them
        # or add NON_MONO_SOLVER=T after the initial check. Re-check and
        # re-apply the critical non-monotonic settings here so they always
        # take precedence, preventing accidental misconfiguration.
        # ======================================================================
        _non_mono_post = False
        for line in lines:
            m = re.match(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]', line)
            if m and m.group(2) == 'NON_MONO_SOLVER':
                if m.group(1).strip().upper() == 'T':
                    _non_mono_post = True
                break
        if _non_mono_post:
            _reapply_keys = build_non_mono_forced_vadat(non_mono_line_treatment)
            # Keys the user set explicitly in vadat_overrides keep their
            # value — same explicit_user_keys semantics as
            # enforce_sn_iin_physics(), which pops them from the forced
            # table before stamping.
            for _key in (vadat_overrides or {}):
                _reapply_keys.pop(str(_key).strip(), None)
            for _key, _forced_val in _reapply_keys.items():
                _current = _get_key_token(lines, _key)
                if _current is not None and _current.strip().upper() != _forced_val.upper():
                    update_keywords_in_lines(lines, {_key: _forced_val}, numeric_format_like_old=True)
                    print(f"  [NON_MONO] {_key} re-forced post-override: {_current} -> {_forced_val}")
                elif _current is None:
                    set_or_append_keyword(lines, _key, _forced_val)
                    print(f"  [NON_MONO] {_key} added post-override: {_forced_val}")

        _jfnk_post = False
        _jfnk_banded_post = False
        for line in lines:
            m = re.match(r'^\s*([^\s]+)\s+\[([A-Za-z0-9_/]+)\]', line)
            if m and m.group(2) in ('DO_JFNK', 'USE_JFNK'):
                if m.group(1).strip().upper() == 'T':
                    _jfnk_post = True
            if m and m.group(2) == 'JFNK_USE_BANDED':
                if m.group(1).strip().upper() == 'T':
                    _jfnk_banded_post = True
        if _jfnk_post:
            _current = _get_key_token(lines, 'SOL_METH')
            _target = 'TRI' if _jfnk_banded_post else 'DIAG'
            if _current is None:
                set_or_append_keyword(lines, 'SOL_METH', _target)
                print(f"  [JFNK] SOL_METH added post-override: {_target}")
            elif _current.strip().upper() != _target:
                update_keywords_in_lines(lines, {'SOL_METH': _target}, numeric_format_like_old=True)
                print(f"  [JFNK] SOL_METH re-forced post-override: {_current} -> {_target}")

        clean = clean_vadat_lines(lines)
        # ======================================================================
        # ENFORCE PHYSICS BEFORE WRITING
        # This neutralizes any contradictory CLI flags or YAML overrides.
        # ======================================================================
        clean = enforce_sn_iin_physics(clean, line_treatment=non_mono_line_treatment,
                                       user_overrides=vadat_overrides)
        require_chg_exch_data(clean, wdir)
        require_die_data(clean, wdir)
        write_lines(vadat_out_path, clean)

        # ── NUM_BNDS / SOL_METH reconciliation ──────────────────────────
        # enforce_sn_iin_physics() no longer forces SOL_METH (it's user-
        # overridable).  But vadat_overrides may have changed it after
        # apply_profile() set num_bnds.  Reconcile num_bnds with the final
        # VADAT SOL_METH so MODEL_SPEC stays consistent with VADAT.
        final_sol_meth = _get_key_token(clean, "SOL_METH")
        if final_sol_meth is not None:
            final_sol_meth_upper = final_sol_meth.strip().upper()
            expected_num_bnds = 1 if final_sol_meth_upper in ("DIAG",) else 3
            if expected_num_bnds != num_bnds:
                print(f"  [RECONCILE] SOL_METH={final_sol_meth_upper} in final VADAT "
                      f"but num_bnds={num_bnds} from apply_profile(). "
                      f"Adjusting num_bnds -> {expected_num_bnds} for MODEL_SPEC consistency.")
                num_bnds = expected_num_bnds

        print_subheader("Iteration control (IN_ITS / IT_SPECIFIER / ADJUST_*)")

        # Keep IN_ITS consistent with VADAT for lambda control.
        in_its_changes = update_in_its_for_profile(wdir, profile,
                                                     iteration_strategy=iteration_strategy)
        if in_its_changes:
            print(f"IN_ITS updates (IN_ITS):")
            for c in in_its_changes:
                print(f"  - {c}")
            print()

        # Create ADJUST_R_DEFAULTS. Runtime R-grid refinement enabled.
        # adjust_r_grid_v4.f has been patched with a sigma guard
        # (warn + clamp to -0.999 instead of STOP) so non-monotonic velocity
        # models can safely use runtime regridding at ionization fronts.
        adjust_r_path = wdir / "ADJUST_R_DEFAULTS"
        adjust_r_content = """! ADJUST_R_DEFAULTS - Runtime R-grid refinement enabled
! Tau-based regridding at ionization fronts (sigma guard patched).
! Format: value BEFORE [KEY] (same as VADAT). No leading spaces for character values.
!
! Conservative parameters for non-monotonic velocity:
! - N_ITS=5: only 5 grid revisions total, then grid is frozen
! - STRT_ITS=15: let non-monotonic solver handle initial transient first
! - FREQ_ITS=10: regrid every 10 iterations (not too frequent)
!
         5 [N_ITS]
! Number of grid revisions to perform (conservative for non-monotonic velocity)
        15 [STRT_ITS]
! First iteration to start grid revision (after non-monotonic solver stabilizes)
        10 [FREQ_ITS]
! Frequency of grid revisions (every 10 iterations)
TAU_SPACE [GRID_METH]
! Regridding method: TAU_SPACE (tau-based, does not directly modify velocity)
UNIFORM [GRID_TYPE]
! Grid type: UNIFORM, MODUN, FIX_NX, REFINE, VTAU
! UNIFORM: equal spacing in log(tau) across full range — handles multiple ionization fronts
  -2.0D0 [FG_MIN]
! Min log(Tau) for fine grid: targets tau~0.01 (deep ejecta)
   0.5D0 [FG_MAX]
! Max log(Tau) for fine grid: up to tau~3 (shock front region)
! Parameters for outer/inner boundary handling
DEFAULT [OB_OPT]
! Outer boundary option: DEFAULT or SPECIFY
DEFAULT [IB_OPT]
! Inner boundary option: DEFAULT or SPECIFY
"""
        adjust_r_path.write_text(adjust_r_content)
        print(f"Wrote ADJUST_R_DEFAULTS (N_ITS=5; runtime R-grid refinement enabled)")
        print()

        # ------------------------------------------------------------------
        # IT_SPECIFIER: iteration cycle control file
        # ------------------------------------------------------------------
        # When present, forces a particular iteration sequence (LAM, FULL,
        # FULL(FIXT), FULL(FIXBA), DO_NG, NORM, etc.). Read every iteration
        # by SPECIFY_IT_CYCLE_V3. Can be modified while the model is running.
        #
        # YAML config: cmfgen.iteration_strategy
        #   - If a list of strings, each string is one line in IT_SPECIFIER
        #     (e.g. ["FULL(FIXT) 5", "DO_NG", "FULL(FIXT) 5", "DO_NG", "NORM"])
        #   - If a string, treated as a single-line strategy
        #   - If None (absent), write "NORM" (let CMFGEN decide via DO_LAM_AUTO)
        #
        it_spec_path = wdir / "IT_SPECIFIER"
        if iteration_strategy is not None:
            if isinstance(iteration_strategy, str):
                it_lines = [iteration_strategy]
            elif isinstance(iteration_strategy, list):
                it_lines = [str(line) for line in iteration_strategy]
            else:
                print(f"  WARNING: iteration_strategy has unexpected type "
                      f"({type(iteration_strategy).__name__}), ignoring")
                it_lines = ["NORM"]
            it_content = "! IT_SPECIFIER - Generated by prepare_VADAT.py from YAML config\n"
            it_content += "! Read every iteration by SPECIFY_IT_CYCLE_V3\n"
            it_content += "! Can be modified while the model is running.\n!\n"
            for line in it_lines:
                stripped = line.strip()
                if stripped and not stripped.startswith("!"):
                    # A bare FULL cycle sets FIXED_T=.FALSE. in
                    # SPECIFY_IT_CYCLE_V3, releasing T at runtime.  T is
                    # pinned to the hydro profile for all runs, so rewrite
                    # bare FULL to FULL(FIXT); FULL(...FIXT...) and
                    # FULL(FIXBA) already keep T fixed and pass through.
                    opt = stripped.split()[0].upper()
                    if opt == "FULL":
                        stripped = stripped.replace(
                            "FULL", "FULL(FIXT)", 1)
                        print(
                            "  WARNING: IT_SPECIFIER 'FULL' rewritten to "
                            "'FULL(FIXT)' -- T is pinned to the hydro "
                            "profile."
                        )
                    elif opt.startswith("FULL(") and "FIXT" not in opt:
                        stripped = re.sub(
                            r"FULL\(", "FULL(FIXT,", stripped,
                            count=1, flags=re.IGNORECASE)
                        print(
                            "  WARNING: IT_SPECIFIER 'FULL(...)' without "
                            "FIXT rewritten to FULL(FIXT,...) -- T is "
                            "pinned to the hydro profile."
                        )
                    it_content += f"{stripped}\n"
                else:
                    it_content += f"{line}\n"
            it_spec_path.write_text(it_content)
            print(f"Wrote IT_SPECIFIER ({len(it_lines)} lines): {it_lines}")
            print()
        else:
            # iteration_strategy absent: let CMFGEN's DO_LAM_AUTO handle
            # the Lambda→BA transition automatically. Remove any existing
            # IT_SPECIFIER so that SPECIFY_IT_CYCLE_V3 returns immediately
            # (line 61: IF IOS .NE. 0 RETURN) and the DO_LAM_AUTO logic
            # in cmfgen_sub.f and solve_for_pops.f controls the iteration
            # mode. This is the standard CMFGEN approach: Lambda iterations
            # until MAXCH < 50%, then auto-switch to Full BA with LAM_VAL/
            # NUM_LAM fallback.
            if it_spec_path.exists():
                it_spec_path.unlink()
                print(f"Removed existing IT_SPECIFIER (DO_LAM_AUTO will control iteration mode)")
            else:
                print(f"No IT_SPECIFIER written (DO_LAM_AUTO will control iteration mode)")
            print()

        # ------------------------------------------------------------------
        # ADJUST_CORRECTIONS: runtime correction scaling file
        # ------------------------------------------------------------------
        # When present, provides relaxation parameters (RELAX, T_LIM,
        # MAX_CHNG, CONSIS_CNT) to scale population corrections at specific
        # depth ranges. Read every iteration by fiddle_pop_corrections_v2.f.
        # Only affects BA (full) iterations, NOT Lambda iterations.
        #
        # YAML config: cmfgen.adjust_corrections
        #   - If a dict, keys map to ADJUST_CORRECTIONS keywords:
        #     l_st, l_end, relax, t_lim, max_chng, consis_cnt
        #   - If None (absent), write default values (no scaling effect)
        #
        # Use MODEL_SPEC ND (hydro ND + 2*N_OB_INS) for L_END cap, not hydro ND alone.
        # CMFGEN's fiddle_pop_corrections_v2.f resets L_END=ND if L_END>ND, so
        # setting L_END to the MODEL_SPEC ND is safe — it covers all depths including
        # boundary insertions, and CMFGEN will clamp it to the actual runtime ND.
        n_points_hydro = int(getattr(hydro, 'params', {}).get('ND', 0))
        n_ob_ins_val = 0
        for _l in lines:
            if "[N_OB_INS]" in _l:
                try:
                    n_ob_ins_val = int(float(_l.split("[")[0].strip()))
                except (ValueError, TypeError):
                    pass
                break
        n_points = n_points_hydro + 2 * max(n_ob_ins_val, 0) if n_points_hydro > 0 else 0
        ac_path = wdir / "ADJUST_CORRECTIONS"
        if adjust_corrections is not None:
            if isinstance(adjust_corrections, dict):
                ac_lines = [
                    "! ADJUST_CORRECTIONS - Generated by prepare_VADAT.py from YAML config",
                    "! Read every iteration by fiddle_pop_corrections_v2.f",
                    "! Only affects BA (full) iterations, NOT Lambda iterations.",
                    "!",
                ]
                # Map YAML keys to ADJUST_CORRECTIONS format
                ac_key_map = {
                    "l_st":       ("L_ST",       "1",     "Beginning depth"),
                    "l_end":      ("L_END",      "ND",    "End depth"),
                    "relax":      ("RELAX",      "1.0D0", "Relaxation parameter (0<Z<2)"),
                    "t_lim":      ("T_LIM",      "0.2D0", "Max fractional T change (0<Z<0.2)"),
                    "max_chng":   ("MAX_CHNG",   "1.0D3", "1+max fractional pop increase (>1)"),
                    "consis_cnt": ("CONSIS_CNT", "0",     "Consistency check frequency"),
                }
                for ykey, (ckey, default, desc) in ac_key_map.items():
                    val = adjust_corrections.get(ykey, default)
                    # Validate numeric parameters against CMFGEN's own validity checks
                    # (see fiddle_pop_corrections_v2.f lines 84-87, 112-140).
                    # CMFGEN silently resets invalid values to defaults, but we warn
                    # loudly here so the user fixes the YAML config instead of
                    # getting silent behaviour changes mid-run.
                    if ckey in ("RELAX", "T_LIM", "MAX_CHNG"):
                        try:
                            fv = float(val)
                        except (ValueError, TypeError):
                            fv = None
                        if fv is not None:
                            if ckey == "RELAX" and not (0.0 < fv < 2.0):
                                print(f"  WARNING: adjust_corrections.relax={val} is "
                                      f"out of range (0 < Z < 2). Clamping to 1.0.")
                                val = 1.0
                            elif ckey == "T_LIM" and not (0.0 <= fv <= 0.2):
                                print(f"  WARNING: adjust_corrections.t_lim={val} is "
                                      f"out of range (0 <= Z <= 0.2). Clamping to 0.2.")
                                val = 0.2
                            elif ckey == "MAX_CHNG" and fv <= 1.0:
                                # CMFGEN check: POP_LIM_VARIABLE .LE. 1.0 => error.
                                # max_chng=1.0 means "0x population increase allowed"
                                # which is physically nonsensical. The likely intent
                                # was "limit changes" — use 3.0 (allow 3x increase)
                                # as a safe conservative default.
                                print(f"  WARNING: adjust_corrections.max_chng={val} is "
                                      f"invalid (must be > 1.0; CMFGEN rejects <= 1.0). "
                                      f"Clamping to 3.0.")
                                val = 3.0
                    if ckey in ("L_ST", "L_END") and n_points > 0:
                        try:
                            if isinstance(val, str) and val.upper() in ("ND", "NDMAX", ""):
                                v = n_points
                            else:
                                v = int(float(val))
                        except (ValueError, TypeError):
                            v = n_points if ckey == "L_END" else 1
                        v = max(1, min(v, n_points))
                        val = v
                    # Format numeric values with D notation if they look like floats
                    val_str = str(val)
                    if ckey in ("RELAX", "T_LIM", "MAX_CHNG") and "D" not in val_str.upper():
                        try:
                            float(val_str)
                            # Format as Fortran D notation: e.g. 0.5 -> 5.0D-01
                            val_str = f"{float(val_str):.4E}".replace("E", "D")
                        except ValueError:
                            pass
                    ac_lines.append(f"{val_str} [{ckey}]")
                    ac_lines.append(f"! {desc}")
                ac_content = "\n".join(ac_lines) + "\n"
                ac_path.write_text(ac_content)
                print(f"Wrote ADJUST_CORRECTIONS: {adjust_corrections}")
                print()
            else:
                print(f"  WARNING: adjust_corrections has unexpected type "
                      f"({type(adjust_corrections).__name__}), ignoring")
        else:
            # No adjust_corrections in YAML: do not write the file at all.
            # fiddle_pop_corrections_v2.f treats a missing file as "no
            # fiddling" — the old default block (T_LIM=0.2) was NOT inert.
            # Per user directive 18-Sep-2026, ADJUST_CORRECTIONS is excluded
            # by default; set cmfgen.adjust_corrections in YAML to enable.
            # Remove a stale copy if one was carried over from base files.
            if ac_path.exists():
                ac_path.unlink()
                print("Removed ADJUST_CORRECTIONS (not configured; corrections unclamped)")
            else:
                print("No ADJUST_CORRECTIONS (default: corrections unclamped)")
            print()

        # init-profile requested rename
        # renamed = maybe_rename_vadat_inint_to_init_bak(wdir, profile, vadat_init_path)
        # if renamed:
        #     print(f"Renamed VADAT_ININT -> {renamed.name}")

        keys = []
        for ln in clean:
            m = _VADAT_LINE_RE.match(ln.strip())
            if m:
                keys.append(m.group(4).strip())
        dup = sorted({k for k in keys if keys.count(k) > 1})
        if dup:
            print("WARNING: duplicate keys in clean VADAT:", ", ".join(dup))

        print_subheader("Final VADAT readback")

        # ── Final configuration, read back from the file just written ──
        # Intermediate prints describe the branch that was active at that
        # point in main(); later overrides (VEL_LAW selection, enforce_
        # sn_iin_physics, vadat_overrides) can legitimately change the
        # values again.  Only the file is fact, so the closing report is
        # generated from the written VADAT itself and cannot disagree.
        try:
            _final_lines = vadat_out_path.read_text().splitlines()
        except OSError:
            _final_lines = clean
        print("[Final configuration] values read back from written VADAT:")
        _final_keys = [
            "VEL_LAW", "RD_IN_R_GRID", "USE_LAM_ES", "DO_NONMONO_SOB",
            "DO_NONMONO_KNOP", "SOL_METH", "DC_METH", "DIF", "IB_METH",
            "VFRAC_MOM", "VFRAC_FG", "FIX_T", "FIX_T_AUTO", "LTE_EST",
            "RSTAR", "RMAX", "VINF", "SN_AGE",
        ]
        for _k in _final_keys:
            _tok = _get_key_token(_final_lines, _k)
            print(f"    {_k:<16} = {_tok if _tok is not None else '(absent)'}")
        # vadat_overrides win over hydro-derived scalars; say so where the
        # file differs from the hydro computation rather than repeating the
        # earlier message that implied the hydro value survived.
        if vadat_overrides:
            for _k in ("RSTAR", "RMAX", "VINF"):
                if _k in vadat_overrides:
                    _tok = _get_key_token(_final_lines, _k)
                    _hv = hydro.params.get(_k)
                    print(f"    NOTE: {_k}={_tok} comes from vadat_overrides; "
                          f"the hydro-derived value was {_hv}")
        # An RSTAR/RMAX override invalidates the tau=2/3 photosphere anchor
        # computed from the hydro luminosity, so flag the conflict rather
        # than letting the earlier photosphere print stand as the record.
        _rstar_tok = _get_key_token(_final_lines, "RSTAR")
        _rstar_fin = _vadat_float_token(_rstar_tok)
        _rstar_hyd = hydro.params.get("RSTAR")
        if (_rstar_fin is not None and _rstar_hyd is not None
                and abs(_rstar_fin - _rstar_hyd) > 1e-6 * max(abs(_rstar_hyd), 1.0)):
            print(f"    WARNING: final RSTAR={_rstar_fin:.6e} differs from the "
                  f"hydro/photosphere value {_rstar_hyd:.6e}; the tau=2/3 "
                  "photosphere anchoring printed earlier no longer applies.")
        _rmax_tok = _get_key_token(_final_lines, "RMAX")
        _rmax_fin = _vadat_float_token(_rmax_tok)
        _rmax_hyd = hydro.params.get("RMAX")
        if (_rmax_fin is not None and _rmax_hyd is not None
                and abs(_rmax_fin - _rmax_hyd) > 1e-6 * max(abs(_rmax_hyd), 1.0)):
            print(f"    WARNING: final RMAX={_rmax_fin:.6e} differs from the "
                  f"hydro/photosphere value {_rmax_hyd:.6e}; the tau=2/3 "
                  "photosphere anchoring printed earlier no longer applies.")
        print()

        print()
        print(f"OK: base={base_vadat_path.name} -> wrote {vadat_out_path.name} (no comments).")
        if profile == "single_final":
            print("  (Single-final: LTE cold-start; T fixed to hydro; Ne free for NLTE.)")
        else:
            print("  (Final: LTE_EST=F, LIN_INT=T; warm restart from *_IN populations.)")
            if out2in_renames:
                print(f"  Renamed {len(out2in_renames)} *OUT -> *_IN (previous run populations as starting point).")
        # One-line confirmation of key scalars written to VADAT (from hydro)
        print(f"  VADAT scalars from hydro: RSTAR={hydro.params['RSTAR']:.6e} RMAX={hydro.params['RMAX']:.6e} "
              f"VINF={hydro.params['VINF']:.2f} SN_AGE={hydro.params['time_days']:.4f} MASS={hydro.params['MASS']:.6f}")
        print()
        if scalar_changes:
            print("Scalar updates:")
            for c in scalar_changes:
                print("  -", c)
            print()
        if mass_changes:
            print("Mass update:")
            for c in mass_changes:
                print("  -", c)
            print()
        if abund_changes:
            print("Abundance updates:")
            for c in abund_changes:
                print("  -", c)
            print()
        if profile_changes:
            print(f"Profile '{profile}' updates:")
            for c in profile_changes:
                print("  -", c)
            print()

        print_subheader("MODEL_SPEC")

        # Update MODEL_SPEC so ND (and NP = ND+NC) match the grid size.
        # When --target-nd is specified, ND = target_nd (the downsampled grid size).
        if model_spec_path.exists() and not args.no_model_spec:
            nd_hydro = int(hydro.params["ND"])
            # For VEL_LAW=12: always use nd_hydro (input_hydro.dat grid size),
            # never nd_written (which may be velocity-downsampled RDINR size).
            # NBND_INS=0 in rv_sn_model_sniin.f: use AMR grid directly without
            # boundary insertions. ND = NHYDRO so MND = ND = NHYDRO, triggering
            # the direct AMR grid path. Boundary resolution is handled by the
            # Python AMR converter's tau-based gridding.
            if effective_vel_law == 12:
                # When USE_HYDRO_GRID=T, rv_sn_model_sniin.f uses the hydro
                # grid directly (R(1:ND)=R_HYDRO(1:ND)) with no boundary
                # insertion.  ND must equal NHYDRO exactly.  N_OB_INS is
                # irrelevant (only used in the tau-based path) but set to 0
                # for consistency.
                if use_hydro_grid:
                    nd_for_model_spec = nd_hydro
                    # Force N_OB_INS=0 and N_IB_INS=0 for consistency
                    set_or_append_keyword(lines, "N_OB_INS", "0")
                    set_or_append_keyword(lines, "N_IB_INS", "0")
                    print(f"  [USE_HYDRO_GRID] MODEL_SPEC ND={nd_hydro} (=NHYDRO, no boundary insertion)")
                else:
                    # rv_sn_model_sniin.f uses NBND_INS=N_OB_INS from VADAT (passed via
                    # cmfgen_sub.f). When target_nd is provided, input_hydro.dat is
                    # downsampled to MND=target_nd-2*NBND_INS points, and we set
                    # ND=target_nd so that NHYDRO==MND triggers the direct AMR grid
                    # path (no log-uniform tau regridding that could extend beyond
                    # the hydro R range).
                    #
                    # Re-read N_OB_INS from VADAT lines AFTER overrides have been
                    # applied (line ~3133), so that config-level N_OB_INS overrides
                    # are respected. The earlier read (line ~3006) used the pre-override
                    # profile default and is only used for input_hydro.dat downsampling.
                    n_ob_ins_after_override = _get_key_token(lines, "N_OB_INS")
                    nbnd_ins_sniin_final = 3  # default in rv_sn_model_sniin.f
                    if n_ob_ins_after_override is not None:
                        try:
                            nbnd_ins_sniin_final = int(float(n_ob_ins_after_override))
                            if nbnd_ins_sniin_final < 0:
                                nbnd_ins_sniin_final = 3  # fallback to default
                        except (ValueError, TypeError):
                            pass
                    NBND_INS_SNIIN = nbnd_ins_sniin_final
                    expected_nd = nd_hydro + 2 * NBND_INS_SNIIN
                    if args.target_nd is not None:
                        if int(args.target_nd) != expected_nd:
                            print(f"[WARNING] target_nd={args.target_nd} does not match actual "
                                  f"SN_HYDRO_DATA ND={nd_hydro} + 2*NBND_INS={2*NBND_INS_SNIIN} "
                                  f"= {expected_nd}. Using actual value to prevent LIN_INTERP crash.")
                        nd_for_model_spec = expected_nd
                    else:
                        nd_for_model_spec = expected_nd
            else:
                # nd_written already reflects target_nd downsampling and any
                # --rdinr-max-dv-kms refinement; it equals nd_hydro when neither
                # applied, so it is the correct VEL_LAW=11 count in all cases.
                nd_for_model_spec = nd_written
            model_spec_changes = update_model_spec(model_spec_path, nd_for_model_spec, nc=None, num_bnds=num_bnds)
            if model_spec_changes:
                print()
                print(f"MODEL_SPEC updates ({model_spec_path.name}):")
                for c in model_spec_changes:
                    print("  -", c)
                if num_bnds == 3:
                    print("  (NUM_BNDS=3 uses TRIDIAG solver; couples adjacent depths for ionisation-front stability.)")
                elif num_bnds == 1:
                    print("  (NUM_BNDS=1 uses DIAG solver; each depth independent.)")
                print()
            # Remove species per --remove-species or config remove_species
            if args.remove_species:
                sp_to_remove = [s.strip() for s in args.remove_species.split(',') if s.strip()]
                if sp_to_remove:
                    sp_changes = remove_species_from_model_spec(model_spec_path, sp_to_remove)
                    if sp_changes:
                        print(f"MODEL_SPEC species removal ({model_spec_path.name}):")
                        for c in sp_changes:
                            print("  -", c)
                        print()
            else:
                print(f"MODEL_SPEC: no updates (file must contain ND and NC keywords).")
        elif args.no_model_spec:
            print("MODEL_SPEC: skipped (--no-model-spec).")
        elif not model_spec_path.exists():
            print(f"MODEL_SPEC: {model_spec_path.name} not found; skipping.")

        print_subheader(f"prepare_VADAT.py [{profile}] complete")
        return 0
    finally:
        _restore_stdout()


if __name__ == "__main__":
    sys.exit(main())
