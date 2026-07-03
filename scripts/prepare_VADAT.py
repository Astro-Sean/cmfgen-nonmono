#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

"""
HPC-side script: update CMFGEN VADAT and write RDINR from SN_HYDRO_DATA.

Only SN_HYDRO_DATA is copied to the HPC (SNEC-to-CMFGEN runs off the HPC and writes
SN_HYDRO_DATA only). This script does everything needed on the HPC:
  - Parse SN_HYDRO_DATA
  - Update VADAT (scalars, mass, abundances, profile, LSTAR)
  - Write RDINR (R grid + V for RD_IN_R_GRID=T; use --no-write-rdinr to skip)
  - Write input_hydro.dat in SNIIN format (VEL_LAW=12)
  - Write T_IN and GAMMAS_IN from Temperature/Electron density when present (LTE initial)
  - Update MODEL_SPEC (ND, NP) to match SN_HYDRO_DATA

Type IIn / metal-rich CSM workflow (LTE initial model then final run):
  - Initial model (LTE): Use VADAT_INIT as the template with your IIn/CSM settings
    (abundances, VEL_LAW=12, SN_T_OPT=USE_HYDRO, etc.). Run prepare_VADAT --init:
    it builds VADAT from VADAT_INIT + SN_HYDRO_DATA and applies the init profile
    (LTE_EST=T for LTE initial level populations, T/Ne fixed to hydro, loose convergence).
    Only initial conditions are set here; the INIT job produces the first model.
  - Final run: After INIT completes, run prepare_VADAT --final. It renames *OUT
    to *_IN (same as cur_cmf/com/out2in.sh), then updates VADAT with the final
    profile (LTE_EST=F, LIN_INT=T, tighter convergence). With LTE_EST=F, CMFGEN
    reads the *_IN files (init population levels); then run the FINAL job.

  Temperature: FIX_T=T always. Temperature is fixed to the input SN_HYDRO_DATA at all depths
  to prevent CMFGEN from solving for temperature and use the hydro temperature instead.

  Electron density: for --init we also fix Ne (FIX_NE=T, FIX_INB_NE=T) to prevent
  huge first-iteration population corrections ("bad initial population guesses")
  in stiff Type IIn/CSM setups. For --final we keep FIX_NE=F so Ne can converge
  self-consistently in NLTE.

  INC_RAD_DECAYS is T when SN_HYDRO_DATA lists one or more isotope mass-fraction blocks
  (Number of isotopes > 0). CMFGEN must read NUC_DECAY_DATA in that case so RD_SN_DATA can
  match isotope names to the ISO table; with INC_RAD_DECAYS F, NUM_ISOTOPES stays zero and
  RD_SN_DATA stops with "isotope data not recognized". Otherwise INC_RAD_DECAYS is F.
  Init uses a higher MAX_CHNG cap than stable for the first population solve.

Profile behavior:
- INIT: base is VADAT_INIT if present, else VADAT_INIT.bak; writes to VADAT.
- FINAL: base is existing VADAT; overwrites VADAT with updated values.

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
RD_IN_R_GRID=T so CMFGEN uses the R grid from RDINR; N_IB_INS and N_OB_INS are set by the init/final profile (e.g. 2 and 5).

Why velocity_law.f can see R(ID)=R(ID+1): Consecutive radii in RDINR come from
SN_HYDRO_DATA. They can be duplicate or nearly identical due to (1) limited
precision when SN_HYDRO_DATA is written or parsed, (2) SNEC fine zoning, or
(3) an older file written before SNEC-to-CMFGEN enforced min spacing. prepare_VADAT
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

# ---------------- constants ----------------
MSUN_CGS = 1.98847e33  # g
DAY = 86400.0          # s
LSUN_ERG_S = 3.828e33  # L_sun in erg/s (for LSTAR)

# Map human-readable element symbols to CMFGEN "DUM" labels.
# Full composition for complex interacting supernova systems: H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe, Ni, Co
# CMFGEN species names from cmfgen.f source code
CMFGEN_LABELS = {
    "H": "HYD", "HE": "HE", "He": "HE",
    "C": "CARB", "N": "NIT", "O": "OXY", "Ne": "NEON",
    "Mg": "MAG", "Si": "SIL", "S": "SUL", "Ar": "ARG", "Ca": "CAL",
    "Ti": "TIT", "Cr": "CHRO", "Fe": "IRON", "Ni": "NICK", "Co": "COB",
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
    "MAG": (12, 24.305),
    "SIL": (14, 28.085),
    "SUL": (16, 32.06),
    "ARG": (18, 39.948),
    "CAL": (20, 40.078),
    "TIT": (22, 47.867),
    "CHRO": (24, 51.996),
    "IRON": (26, 55.845),
    "COB": (27, 58.933),
    "NICK": (28, 58.693),
}

AMU_G = 1.66053906660e-24  # g, unified atomic mass unit

_NUM_TOKEN_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$")
_VADAT_LINE_RE = re.compile(r"^(\s*)(\S+)(\s+)\[([^\]]+)\](.*)$")

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
    """Parse an SNEC .xg file (same format as rho.xg). Returns dict time_sec -> (col_a, col_b)."""
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
    # HERACLES may have negative times for pre-explosion phases
    if time_days < 0:
        print(f"[SN_HYDRO_DATA Parse] Negative time detected: {time_days:.2f} days")
        print(f"  Converting to positive time since explosion bomb introduction")
        time_days = abs(time_days)
    
    # Optional: Energy Lum [lsun] (inner boundary luminosity; written by SNEC-to-CMFGEN)
    energy_lum_lsun: Optional[float] = None
    for line in lines:
        m = re.search(r"Energy\s+Lum\s*\[lsun\]\s*:\s*([0-9EeDd+\-\.]+)", line, re.IGNORECASE)
        if m:
            energy_lum_lsun = float(m.group(1).replace("D", "E").replace("d", "e"))
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
        vals: List[str] = []
        i = start_index
        while i < len(lines) and re.search(r"[0-9]", lines[i]):
            vals.extend(re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?", lines[i]))
            i += 1
        return np.array([float(v.replace("D", "E").replace("d", "e")) for v in vals], dtype=float)

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

        arr = np.array([float(v.replace("D", "E").replace("d", "e")) for v in vals], dtype=float)[:N]
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
    )


def print_hydro_summary(hydro: HydroData, prefix: str = "SN_HYDRO_DATA") -> None:
    """Print a short summary of parsed or written hydro data for verification."""
    nd = int(hydro.params["ND"])
    r_cm = hydro.radius_cm
    v_kms = np.abs(hydro.vel_cms) / 1e5
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


def recommended_t_exc(hydro: HydroData) -> Optional[float]:
    """
    Recommend T_EXC (in 10^4 K) from SN_HYDRO_DATA when present.
    CMFGEN uses T_EXC in the *outer region* as minimum excitation T for LTE
    initial estimates, so use the outer-boundary temperature (not the global
    minimum). Slightly below outer T so the real structure is not clamped,
    with a floor for numerical stability.
    """
    t_ob = hydro.params.get("T_OB_104")
    if t_ob is None:
        return None
    # Slightly below outer-boundary T; floor at 0.15 (1500 K) to avoid underflows
    return max(0.30, min(t_ob * 0.95, 10.0))  # 3000 K floor for NLTE coverage


def recommended_t_exc_init(hydro: HydroData) -> Optional[float]:
    """
    Recommend T_EXC (in 10^4 K) for the *init* run when ~half the model is NLTE.
    Uses T at mid-depth so that the inner (LTE) half keeps full LTE at local T,
    and the outer (NLTE) half gets the excitation floor in SET_DC_LTE_V2 instead
    of full LTE—giving a softer initial guess and often faster convergence.
    """
    if hydro.temp_104k is None or hydro.temp_104k.size == 0:
        return recommended_t_exc(hydro)
    nd = int(hydro.params["ND"])
    i_mid = nd // 2
    t_mid = float(hydro.temp_104k[i_mid])
    return max(2.0, min(t_mid, 10.0))


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
    """
    r_cm = hydro.radius_cm
    rho = hydro.rho_cgs
    n = len(r_cm)
    if hydro.kappa_cgs is not None and len(hydro.kappa_cgs) == len(r_cm):
        kappa = np.asarray(hydro.kappa_cgs, float)
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
    for raw in lines:
        s = _strip_inline_comment(raw).strip()
        if not s:
            continue
        m = _VADAT_LINE_RE.match(s)
        if not m:
            continue
        token = m.group(2).strip()
        key = m.group(4).strip()
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


# ---------- IN_ITS helpers ----------
def update_in_its_for_profile(wdir: Path, profile: str) -> List[str]:
    """
    Keep IN_ITS and VADAT consistent for LAMBDA control.

    We have seen runs become unstable when VADAT has DO_LAM_IT=T but IN_ITS has DO_LAM_IT=F
    (or vice-versa), effectively changing the iteration scheme mid-run.

    Policy:
      - init:   DO_LAM_IT=T, DO_LAM_AUTO=F, DO_GT_AUTO=T, DO_T_AUTO=F, D2_EQ_D1=T
      - stable: DO_LAM_IT=T, DO_LAM_AUTO=F, DO_GT_AUTO=T, DO_T_AUTO=F, D2_EQ_D1=T
      - final:  DO_LAM_IT=T, DO_LAM_AUTO=T,  DO_GT_AUTO=T, DO_T_AUTO=F, D2_EQ_D1=T

    Temperature auto-updates (DO_T_AUTO=F) are disabled since FIX_T=T fixes T to SN_HYDRO_DATA.
    DO_GT_AUTO=T and D2_EQ_D1=T are enabled for better convergence in stiff Type IIn/CSM setups.
    """
    in_its_path = wdir / "IN_ITS"
    if not in_its_path.exists():
        return []

    lines = read_lines(in_its_path)
    want_lam = "T"
    want_auto = "T" if profile == "final" else "F"
    want_gt = "T"
    want_t = "F"
    want_d2 = "T"  # Always enable D2_EQ_D1 for stability (tie depth-2 pops to depth-1)
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
MIN_DR_RDINR_1E10CM = 1.0e-10


def write_rdinr(hydro: HydroData, path: Path) -> None:
    """
    Write RDINR for RD_IN_R_GRID=T (same format as departure coefficient files).
    SET_RV_HYDRO_MODEL_V3 reads: '!Format date' line, header (R_ND, LUM, NOLD, ND),
    then for each depth: one line of 6 values (R,4 dummies,V), then READ(T1,J=1,NOLD).
    Some compilers advance to next record when NOLD=0, causing misalignment. Use NOLD=1
    and one dummy line per depth so the second read consumes a record. R in 10^10 cm,
    V in km/s; order outer->inner.
    
    Note: Grid refinement now happens in convert_heracles_to_cmfgen.py (Dessart method),
    so this function just writes the already-refined grid from SN_HYDRO_DATA.
    """
    nd = int(hydro.params["ND"])
    if len(hydro.radius_cm) < nd or len(hydro.vel_cms) < nd:
        raise ValueError(
            f"RDINR: need at least {nd} points, got radius_cm={len(hydro.radius_cm)}, vel_cms={len(hydro.vel_cms)}"
        )
    r_oi = (hydro.radius_cm[::-1] / 1e10).astype(float, copy=True)
    v_oi = np.abs(hydro.vel_cms[::-1]) / 1e5
    
    # Enforce minimum spacing (outer->inner: r_oi[i] > r_oi[i+1]); avoids velocity_law.f SIGFPE
    for i in range(nd - 1):
        dr = r_oi[i] - r_oi[i + 1]
        if dr < MIN_DR_RDINR_1E10CM:
            r_oi[i + 1] = r_oi[i] - MIN_DR_RDINR_1E10CM
    
    with open(path, "w", newline="\n") as f:
        f.write(" 24-FEB-2004                    !Format date\n")
        # NOLD=1 so code does READ(T1,J=1,1) and consumes one value per depth
        f.write(f" {float(r_oi[-1]):15.7E}    {1.0:11.4E}     {1:4d}    {nd:4d}\n")
        for i in range(nd):
            f.write(
                f" {float(r_oi[i]):15.7E} {0.0:15.7E} {0.0:15.7E} {0.0:15.7E} {0.0:15.7E} {float(v_oi[i]):15.7E}\n"
            )
            # One dummy value so READ(T1,J=1,NOLD) with NOLD=1 consumes a record
            f.write(" 0.0\n")


def write_input_hydro_dat(hydro: HydroData, path: Path) -> None:
    """
    Write input_hydro.dat (R, rho, T per depth) for workflows that use it.
    Format: ADD_ACC_ZONE (F), header line, NHYDRO, bogus line, then NHYDRO lines of
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
        f.write(" F\n")  # ADD_ACC_ZONE
        f.write("! R(1e10 cm) dummy rho(g/cm^3) T(K) from SN_HYDRO_DATA\n")
        f.write(f" {n}\n")
        f.write("! NHYDRO lines below (outer -> inner)\n")
        # Outer -> inner so set_abund_clump stores first line at R_HYDRO(NHYDRO), last at R_HYDRO(1)
        for i in range(n - 1, -1, -1):
            f.write(f" {r_1e10[i]:15.7E} {0.0:15.7E} {rho[i]:15.7E} {t_k[i]:15.7E}\n")


# Electron-scattering opacity (cm^2/g) fallback when Kappa section missing in SN_HYDRO_DATA
KAPPA_ES_CGS = 0.34


def write_input_hydro_dat_sniin(hydro: HydroData, path: Path) -> None:
    """
    Write input_hydro.dat in the format read by RV_SN_MODEL_SNIIN (VEL_LAW=12) and
    SET_ABUND_CLUMP (cols 1-4 per line: R(1e10 cm), dummy, RHO(g/cm^3), T(K)).
    Format: ADD_ACC_ZONE (F), header, NHYDRO, bogus line, then per point (outer->inner):
    R, V, RHO, T(K), XL, XTAU, TAU.
    """
    r_cm = hydro.radius_cm
    rho = hydro.rho_cgs
    n = len(r_cm)
    if n == 0 or len(rho) < n:
        raise ValueError("write_input_hydro_dat_sniin: need radius and density")
    if hydro.temp_104k is None or len(hydro.temp_104k) < n:
        raise ValueError("write_input_hydro_dat_sniin: SN_HYDRO_DATA must have Temperature section for VEL_LAW=12")
    r_1e10 = r_cm / 1e10
    v_kms = np.abs(hydro.vel_cms) / 1e5  # cm/s -> km/s
    t_k = hydro.temp_104k[:n] * 1e4  # 10^4 K -> K
    if hydro.kappa_cgs is not None and len(hydro.kappa_cgs) >= n:
        kappa = hydro.kappa_cgs[:n]
    else:
        kappa = np.full(n, KAPPA_ES_CGS)
    # Tau from outer to inner: tau[n-1]=0, tau[j] = tau[j+1] + 0.5*(k*rho)[j,j+1] * (r[j+1]-r[j])
    # The outer tau must be > 0 because RV_SN_MODEL_SNIIN takes LOG(TAU_HYDRO(1)) to build a
    # geometric tau grid; a zero outer tau produces NaN/Inf and SIGSEGV in LIN_INTERP.
    tau_oi = np.zeros(n)
    for j in range(n - 2, -1, -1):
        dr = r_cm[j + 1] - r_cm[j]  # > 0 (outer to inner)
        krho_avg = 0.5 * (kappa[j] * rho[j] + kappa[j + 1] * rho[j + 1])
        tau_oi[j] = tau_oi[j + 1] + krho_avg * dr
    tau_oi[-1] = max(tau_oi[-1], 1.0e-30)
    # Print so we can verify the data being written matches what was read from SN_HYDRO_DATA
    print(f"[input_hydro.dat] NHYDRO={n} (from SN_HYDRO_DATA)")
    print(f"[input_hydro.dat] R(1e10 cm): min={r_1e10.min():.6e} max={r_1e10.max():.6e}")
    print(f"[input_hydro.dat] RHO(g/cm^3): min={rho.min():.6e} max={rho.max():.6e}")
    print(f"[input_hydro.dat] T(K): min={t_k.min():.2f} max={t_k.max():.2f}")
    with open(path, "w", newline="\n") as f:
        f.write(" F\n")  # ADD_ACC_ZONE
        f.write("! R(1e10 cm) V(km/s) RHO(g/cm^3) T(K) XL XTAU TAU (cols 3,4 for SET_ABUND_CLUMP)\n")
        f.write(f" {n}\n")
        f.write("! NHYDRO lines below (outer -> inner)\n")
        for i in range(n - 1, -1, -1):
            xl, xtau = 0.0, 0.0
            tau_val = tau_oi[i]
            f.write(
                f" {r_1e10[i]:15.7E} {v_kms[i]:15.7E} {rho[i]:15.7E} {t_k[i]:15.7E} {xl:15.7E} {xtau:15.7E} {tau_val:15.7E}\n"
            )
    print(f"Wrote {path.name} (VEL_LAW=12 SNIIN format): ND={n}  R(1e10 cm) outer->inner, TAU outer={tau_oi[-1]:.4e} inner={tau_oi[0]:.4e}")


def write_t_in_from_hydro(hydro: HydroData, path: Path, *, nold: int = 0, t_max_104: float = 1.0) -> None:
    """
    Write T_IN from HydroData (format read by REGRID_T_ED / INIT_TEMP_V2).
    Format matches star_init/NT35000_logg425: leading blank, '!Format date' line,
    header (R_core, LUMOLD, NOLD, ND), then per depth: one line of 8 columns
    (R, dummy, ED, T, ION_FRAC, VEL, CLUMP, depth_index) in 12.7E style, then
    NOLD values (5 per line) or single " 0.0" when NOLD=0. Order outer->inner.

    Args:
        t_max_104: Maximum temperature in units of 10^4 K (default: 1.0 = 10,000 K)
    """
    n = len(hydro.radius_cm)
    if n == 0:
        raise ValueError("write_t_in_from_hydro: empty hydro")
    if hydro.temp_104k is None or len(hydro.temp_104k) < n:
        raise ValueError("write_t_in_from_hydro: SN_HYDRO_DATA must have Temperature section")
    if hydro.ed_cgs is None or len(hydro.ed_cgs) < n:
        raise ValueError("write_t_in_from_hydro: SN_HYDRO_DATA must have Electron density section")
    r_1e10 = hydro.radius_cm / 1e10
    v_kms = np.abs(hydro.vel_cms) / 1e5
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
        f.write(f"   {r_core:12.7E}     {lum_old:12.7E} {nold:8d} {n:8d}\n")
        f.write("\n")
        for i in range(n - 1, -1, -1):
            r = max(r_1e10[i], 1e-30)
            e = max(ed[i], 1e-30)
            t = min(max(t_104[i], 0.3), t_max_104)  # floor at 3000 K (0.3 x 10^4) for NLTE consistency with SEDONA
            v = v_kms[i]
            depth_idx = n - i
            f.write(
                f"   {r:12.7E} {v:12.7E} {e:12.7E} {t:12.7E} {1.0:12.7E} {v:12.7E} {1.0:12.7E} {depth_idx:6d} \n"
            )
            if nold > 0:
                for k in range(0, nold, 5):
                    chunk = [0.0] * min(5, nold - k)
                    f.write(" " + " ".join(f"{x:12.7E}" for x in chunk) + "\n")
            else:
                f.write(" 0.0\n")
            f.write("\n")
    print(f"Wrote {path.name}: ND={n}  R(10^10 cm) {r_1e10[-1]:.4e} -> {r_core:.4e}  LUMOLD={lum_old:.4e} L_sun (for INIT_TEMP_V2)")


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


def ensure_die_keys_for_model_spec(lines: List[str], model_spec_path: Path) -> List[str]:
    """
    Ensure VADAT contains DIE_* keys for ions enabled in MODEL_SPEC.
    Missing keys can trigger GET_KEY_STRING failures at startup.
    """
    if not model_spec_path.exists():
        return []

    text = model_spec_path.read_text(errors="ignore")
    changes: List[str] = []
    ion_to_die = {
        "SkI_ISF": "DIE_SkI",
        "Sk2_ISF": "DIE_Sk2",
        "FeI_ISF": "DIE_FeI",
        "Fe2_ISF": "DIE_Fe2",
    }
    for ion_tag, die_key in ion_to_die.items():
        if f"[{ion_tag}]" in text and _get_key_token(lines, die_key) is None:
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

def update_abundances_mass_fraction(lines: List[str], ave_mass_frac: Dict[str, float]) -> List[str]:
    changes: List[str] = []
    key_idx = find_existing_keys(lines)

    for dum, frac in ave_mass_frac.items():
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
# N_TYPE determines NUM_BNDS: G_ONLY -> 1 (diagonal), N_ON_J/MIXED -> 3 (band). SOL_METH must match.
def _num_bnds_for_n_type(n_type: str) -> int:
    return 1 if (n_type or "").strip().upper() == "G_ONLY" else 3


def apply_profile(
    lines: List[str],
    profile: str,
    convergence_mode: str = "default",
    num_isotopes: int = 0,
) -> Tuple[List[str], int]:
    """
    Apply init, stable, or final VADAT profile. Recommended procedure for SNEC input:
      - Init:   fix T/Ne for stability, LTE initial populations (3 its)
      - Stable: pure Lambda iterations to settle populations before Non-Lambda (20 its)
      - Final:  fix T below photosphere; allow Ne to vary; full Non-Lambda convergence

    num_isotopes: from SN_HYDRO_DATA header. If > 0, INC_RAD_DECAYS must be T so CMFGEN
    reads NUC_DECAY_DATA and RD_SN_DATA can match isotope blocks to the ISO table.
    """
    if profile not in {"init", "stable", "final", "single_final"}:
        raise ValueError(f"Unknown profile: {profile}")
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
        # Single-final: skip INIT entirely; start from LTE populations (LTE_EST=T).
        # Uses Lambda iterations initially; CMFGEN auto-switches to Non-Lambda once stable.
        # Dessart+2015: T fixed to hydro, coherent ES, G_ONLY/DIAG solver throughout.
        # DO_NG=F during LTE bootstrap to avoid divergent Ng on unstable initial populations.
        desired = {
            "LTE_EST": "T",
            "SN_T_OPT": "USE_T_IN",
            "LSTAR": "0.0D0",
            "INC_RAD_DECAYS": inc_rad_decays,
            "INC_SHOCK_POWER": "F",   # Dessart+2015: shock heating already in hydro T(r)
            "ALL_FREQ": "T",

            "RD_CF_FILE": "F",
            "MIN_CF": "1.0D-1",        # Raised to 0.1 to avoid MOM_J_REL_V9 NaN errors at extreme far-IR
            "MAX_CF": "13.5D0",
            "FRAC_SP": "1.08D0",
            "AMP_FAC": "1.04D0",
            "MAX_BF": "0.05D0",
            "CROSS": "F",
            "V_CROSS": "500.0D0",
            "SIG_GAU_KMS": "500.0D0",

            # Dessart+2015: relativistic moment solver required for SN models
            "USE_JREL": "T",
            "USE_J_REL": "T",
            "INCL_REL": "T",
            "INCL_ADV_TRANS": "T",
            "USE_FRM_REL": "T",
            "NMONO_J": "F",           # monotonic J for LTE bootstrap stability

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "75.0D0",      # [OPT] Elevated turbulence for shock regions (not in STELLA T)
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)
            "JC_W_EDD": "T",
            "THK_CONT": "F",
            "USE_LAM_ES": "T",        # [OPT] was F; enable Lambda iteration for J (stable for SN models)
            "COH_ES": "T",            # Dessart+2015: coherent ES throughout
            "MIX_COH": "F",

            # Dessart+2015: T fixed to hydro; Ne free for self-consistency
            "FIX_T": "T",
            "FIX_NE": "F",            # Ne free: self-consistent with radiation field (was T)
            "FIX_IMP": "T",
            "FIX_INB_T": "F",
            "FIX_INB_NE": "F",
            "FIX_X_DPTH": "0",

            # N_ON_J/TRIDIAG: consistent with NUM_BNDS=3 (off-diagonal corrections for blanketed IIn)
            # NOTE: TRIDIAG couples nearest-neighbour depths, propagating corrections across shock regions
            # more effectively than DIAG, roughly halving iteration count for complex velocity fields.
            "SCALE_OPT": "MAJOR",
            "MAX_LIN": "2.0D0",
            "MAX_CHNG": "1.0D4",
            "COMP_BA": "T",
            "N_FIX_BA": "1",          # [OPT] was 50; recompute BA every iteration for N_ON_J stability
            "BA_CHK_FAC": "1.0D-02", # [OPT] was 1.0D-05; standard threshold for N_ON_J
            "FIX_BA": "1.0D-02",     # [OPT] was 5.0D-04; standard 1% threshold for N_ON_J
            "LAM_VAL": "2.0D2",
            "NUM_LAM": "40",

            "GLOBAL_LINE": "SOB",
            "GF_CUT": "5.0D-3",
            "GF_LEV_CUT": "9",
            "AT_CUT": "22",
            "MIN_TRANS": "15",
            "THK_LINE": "F",
            "CHK_L_POS": "T",
            "NEG_OPAC_OPT": "SRCE_CHK",
            "SCL_LN": "T",
            "SCL_LN_FAC": "0.1D0",
            "ES_FAC": "0.8D0",        # [OPT] was 0.6; standard damping for robustness

            # LTE cold start: LTE_EST=T, no *_IN files required
            "LIN_INT": "F",
            "IT_ON_T": "F",
            "POP_SCALE": "T",
            "T_MIN": "3.0D-01",
            "T_EXC": "3.00D-1",
            "GREY_TAU": "10.0D0",
            "METHOD": "ZERO",
            "INC_RAY": "F",

            "REV_RGRID": "F",
            "RD_IN_R_GRID": "F",
            "N_IB_INS": "10",         # CMFGEN hard limit is 16; 10 is safe
            "N_OB_INS": "10",
            "N_TYPE": "G_ONLY",      # Changed from N_ON_J due to convergence issues (excessive iteration count with NaN errors)

            "INC_GRID": "F",        # Disabled to avoid SIGMA > -1 limitation in EXTEND_VTSIGMA for interacting SNe
            "VFRAC_FG": "0.50D0",    # [OPT] was 0.4; coarser ray spacing for SN (prevents ND adjustment)
            "VFRAC_MOM": "0.50D0",   # [OPT] was 0.4; coarser moment spacing for SN (prevents ND adjustment)
            "OB_METH": "HONJ",

            "ACC_END": "5.0D0",
            "N_INS": "2",
            "ST_INT": "1",
            "END_INT": "30",
            "ND_QUAD": "30",
            "INTERP_TYPE": "LOG",
            "N_PAR": "150",

            "MAX_DOP": "8.0D0",
            "FRAC_DOP": "1.5D0",
            "dV_CMF_PROF": "80.0D0",
            "dV_CMF_WING": "200.0D0",
            "ES_WING_EXT": "12000.0D0",
            "R_CMF_WING_EXT": "4.0D0",

            "OBS_EXT_RAT": "1.3D0",
            "dV_OBS_PROF": "40.0D0",
            "dV_OBS_WING": "150.0D0",
            "dV_OBS_BIG": "800.0D0",

            "EPS_TERM": "1.0D-01",    # relaxed convergence for DIAG stability (was 0.05; 0.1 gives more room)
            "ACC_F": "5.0D-02",       # relaxed Eddington convergence for DIAG stability (was 0.01)

            # Lambda auto-switch; no Ng on unstable LTE start
            "DO_LAM_IT": "T",
            "DO_NG": "F",             # disabled for LTE bootstrap; re-enable once populations settle
            "IBEG_NG": "50",          # start Ng after 50 iterations of Lambda settling
            "BW_NG": "30",
            "ITS/NG": "4",
        }

    elif profile == "init":
        # Very loose parameters for initial run to maximise chance of convergence (final profile tightens).
        # Dessart+2015: T fixed to hydro, coherent ES, G_ONLY/DIAG solver.
        # FIX_NE=T during init: prevents huge first-iteration corrections from unstable LTE start.
        # NMONO_J=F: monotonic J during LTE_EST phase; non-monotonic only in final run.
        desired = {
            "SN_T_OPT": "USE_T_IN",
            "LSTAR": "0.0D0",         # SN: 0 (no central L; heating from hydro T/ED)
            "INC_RAD_DECAYS": inc_rad_decays,
            "INC_SHOCK_POWER": "F",   # Dessart+2015: shock heating already in hydro T(r)
            "ALL_FREQ": "T",

            # [OPT] Continuum frequency grid for IIn models
            "RD_CF_FILE": "F",      # Build frequency grid internally
            # Avoid an observed edge-case: "HNU zero for frequency: 3.00000000E-03" at the low-nu edge.
            # Raise MIN_CF so the solver doesn't hit the extreme far-IR tail during bootstrap.
            "MIN_CF": "1.0D-1",        # Raised to 0.1 to avoid MOM_J_REL_V9 NaN errors at extreme far-IR
            # Avoid an edge-case observed in OUTGEN: "HNU zero for frequency: ~1.5E+01".
            # Empirically this can trigger at the high-frequency end of the continuum grid.
            # Back off further from the upper edge to keep the init bootstrap stable.
            "MAX_CF": "13.5D0",
            "FRAC_SP": "1.08D0",    # [OPT] was 1.10; smoother low-nu spacing
            "AMP_FAC": "1.04D0",    # [OPT] was 1.05; smoother high-nu spacing
            "MAX_BF": "0.05D0",     # [OPT] was 0.10; tighter BF edge resolution
            "CROSS": "F",          # Sparse chi/eta with interpolation
            "V_CROSS": "500.0D0",   # 500 km/s for 500 km/s feature resolution
            "SIG_GAU_KMS": "500.0D0",  # match V_CROSS for consistency

            # SN models require USE_J_REL (moment solver) or USE_DJDT_RTE. For USE_J_REL, CMFGEN
            # also requires INCL_REL and INCL_ADV_TRANS to be enabled (see CHECK_PARAM_CONSISTENCY).
            # Set both USE_JREL and USE_J_REL spellings for compatibility with different builds.
            "USE_JREL": "T",
            "USE_J_REL": "T",
            "INCL_REL": "T",
            "INCL_ADV_TRANS": "T",
            # When using the relativistic moment solver, also enable formal-relativistic flux evaluation.
            "USE_FRM_REL": "T",
            "NMONO_J": "F",          # F = use monotonic J during LTE_EST; non-monotonic J enabled in final profile
            "JC_W_EDD": "T",
            "THK_CONT": "F",         # Thick outer BC during init

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "75.0D0",      # [OPT] Elevated turbulence for shock regions (not in STELLA T)
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)
            # INIT bootstrap: disable USE_LAM_ES. We have observed "HNU zero for frequency" failures
            # early in great iteration 1 along the formal-solution path when USE_LAM_ES is enabled.
            # FINAL retains USE_LAM_ES=T for the intended Type IIn treatment once restartable state exists.
            "USE_LAM_ES": "F",
            "COH_ES": "T",
            "MIX_COH": "F",          # F = no coherent/incoherent mix in linearisation (more stable)

            # Dessart+2015: T fixed to hydro at all depths; Ne free for self-consistency
            "FIX_T": "T",
            "FIX_NE": "F",           # Ne free: self-consistent with radiation field (was T)
            "FIX_IMP": "T",          # T = auto-fix impurity species during init
            "FIX_INB_T": "F",        # inner boundary T from hydro (FIX_T=T covers all depths)
            "FIX_INB_NE": "F",
            "FIX_X_DPTH": "0",       # T fixed at ALL depths (consistent with FIX_T=T)

            # G_ONLY/DIAG: consistent with NUM_BNDS=1 (Dessart+2015 style)
            # NOTE: With DIAG (no depth coupling), BA must be VERY accurate before freezing
            # or the uncorrected errors cause runaway oscillations (10^8% divergences).
            "SCALE_OPT": "MAJOR",
            "MAX_LIN": "2.0D0",
            "MAX_CHNG": "1.0D4",     # finite hard cap prevents runaway corrections
            "COMP_BA": "T",
            "N_FIX_BA": "50",        # many BA recomputes (was 5; need >25 for DIAG stability)
            "BA_CHK_FAC": "1.0D-05", # check BA every iteration when frozen (tighter than 5e-4)
            "FIX_BA": "5.0D-04",     # freeze only at < 0.05% error (was 5%; 100x tighter for DIAG)
            "LAM_VAL": "600.0D0",    # raise Lambda fallback threshold for IIn
            "NUM_LAM": "2",

            # [OPT] Bound-bound & blanketing for IIn models
            "GLOBAL_LINE": "SOB",   # Sobolev rates
            "GF_CUT": "5.0D-3",      # [OPT] was 1e-3; remove weaker lines
            "GF_LEV_CUT": "9",       # Apply GF_CUT above this level
            "AT_CUT": "22",          # Min atomic number for GF cut
            "MIN_TRANS": "15",       # [OPT] was 10; fewer marginal levels
            "THK_LINE": "F",         # Thick-line outer BC off
            "CHK_L_POS": "T",        # Prevent negative total opacity
            "NEG_OPAC_OPT": "SRCE_CHK", # Limit source when chi negative
            "SCL_LN": "T",           # Ignore lines in continuum opacity
            "SCL_LN_FAC": "0.1D0",   # Line scale-out limit (Doppler frac)
            "MIX_COH": "F",          # [OPT] was T; disable for convergence phase
            "ES_FAC": "0.6D0",       # Increased damping for ES oscillation

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "75.0D0",      # [OPT] Elevated turbulence for shock regions (not in STELLA T)
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)

            # LTE initial populations; when ~half the model is NLTE, T_EXC from mid-depth softens outer-half guess (see recommended_t_exc_init).
            "LTE_EST": "T",          # T = set departure coefficients to LTE (outer half uses T_EXC floor in SET_DC_LTE_V2)
            "IT_ON_T": "F",          # F = do not iterate on initial T (skip grey temperature calculation)
            "POP_SCALE": "T",        # T = scale populations so conservation eq. satisfied (consistent with abundances)
            # Raise the temperature floor. Very low T floors (e.g. 1e-2 -> 100 K) can trigger
            # underflow/zeroing in high-frequency Planck factors, which shows up as
            # "HNU zero for frequency" errors in OUTGEN.
            "T_MIN": "3.0D-01",      # 3000 K (unified floor for stability)
            "T_EXC": "3.00D-1",      # Excitation temperature (3000 K)
            "GREY_TAU": "10.0D0",     # Grey optical depth (SNA stable setting)
            "METHOD": "ZERO",        # Continuum tau method (SNA stable setting)
            "INC_RAY": "F",         # Rayleigh scattering off (SNA stable setting)

            "REV_RGRID": "T",       # Enable automatic grid revision to handle problematic regions (was F)
            "RD_IN_R_GRID": "F",
            "N_IB_INS": "15",        # More boundary points for DIAG solver stability (was 10; CMFGEN limit is 16)
            "N_OB_INS": "15",        # More boundary points for DIAG solver stability (was 10; CMFGEN limit is 16)
            "OB_RAT": "1.5D0",       # Outer boundary step ratio for grid stability
            "IB_RAT": "2.0D0",       # Inner boundary step ratio for grid stability
            # NOTE: ADJUST_R_DEFAULTS file is created separately, not a VADAT keyword
            # More robust bootstrap for fixed-T SN models.
            # NOTE: CMFGEN may print "HNU zero for frequency" when H_nu ~ 0 at some depths; this is
            # often transient and not necessarily fatal.
            "N_TYPE": "G_ONLY",      # Changed from N_ON_J due to convergence issues (excessive iteration count with NaN errors)

            # J stability: enable grid revision for init to ensure monotonic grid
            "INC_GRID": "F",        # Disabled to avoid SIGMA > -1 limitation in EXTEND_VTSIGMA for interacting SNe
            "VFRAC_FG": "0.50D0",    # [OPT] was 0.4; coarser ray spacing for SN (prevents ND adjustment)
            "VFRAC_MOM": "0.50D0",   # [OPT] was 0.4; coarser moment spacing for SN (prevents ND adjustment)
            "OB_METH": "HONJ",   # Use HONJ for outer boundary (required by MOM_JREL_V9)

            # [OPT] Spatial grid settings for IIn models
            "ALL_FREQ": "T",       # Apply accuracy at all nu
            "ACC_END": "5.0D0",     # Extend enhanced accuracy to more frequencies
            "N_INS": "2",            # Extra depths (if INC_GRID=T)
            "ST_INT": "1",          # First depth for insertion
            "END_INT": "30",         # Last depth for insertion
            "ND_QUAD": "30",         # [OPT] was 50; fewer inner quadrature depths
            "INTERP_TYPE": "LOG",  # Log R interpolation
            "N_PAR": "150",          # [OPT] was 300; fewer BA update intervals

            # [OPT] Line broadening for high-VINF IIn models
            "MAX_DOP": "8.0D0",         # [OPT] was 6.0; wider resonance zone for strong line overlap
            "FRAC_DOP": "1.5D0",        # [OPT] was 1.0; coarser in-resonance sampling
            "dV_CMF_PROF": "80.0D0",    # [OPT] was 200; resolve 100 km/s features
            "dV_CMF_WING": "200.0D0",   # [OPT] was 300; resolve 100 km/s features
            "ES_WING_EXT": "12000.0D0", # [OPT] was 2500; MUST exceed VINF=10472 km/s
            "R_CMF_WING_EXT": "4.0D0",  # [OPT] was 3.0; wider coherent ES wing

            # [OPT] Observed-frame spectrum for high-VINF IIn
            "OBS_EXT_RAT": "1.3D0",   # [OPT] was 1.1; wider observed-frame coverage
            "dV_OBS_PROF": "40.0D0",  # [OPT] was 100; resolve 100 km/s features
            "dV_OBS_WING": "150.0D0", # [OPT] was 400; resolve 100 km/s features
            "dV_OBS_BIG": "800.0D0",  # [OPT] was 2000; resolve 100 km/s features

            # Convergence settings for init
            "EPS_TERM": "5.0D-02",   # relaxed but not too loose; consistent with stable phase
            "ACC_F": "1.0D-02",      # tighter Eddington convergence reduces late-iteration oscillations
            "ACC_END": "5.0D0",

            "DO_LAM_IT": "T",        # use Lambda iterations during init

            # INIT: disable Ng acceleration (can diverge on bad LTE initial populations)
            "DO_NG": "F",
            "IBEG_NG": "50",         # if DO_NG later enabled, start after 50 iterations
            "BW_NG": "30",
            "ITS/NG": "4",
        }

    elif profile == "stable":
        # Pure Lambda-iteration phase between init and final.
        # Uses populations from init, runs pure Lambda to settle NLTE populations
        # before allowing Non-Lambda iterations in final.
        desired = {
            "SN_T_OPT": "USE_T_IN",
            "LSTAR": "0.0D0",
            "INC_RAD_DECAYS": inc_rad_decays,
            "INC_SHOCK_POWER": "F",
            "ALL_FREQ": "T",

            # [OPT] Continuum frequency grid for IIn models
            "RD_CF_FILE": "F",      # Build frequency grid internally
            "MIN_CF": "1.0D-1",        # Raised to 0.1 to avoid MOM_J_REL_V9 NaN errors at extreme far-IR
            "MAX_CF": "13.5D0",
            "FRAC_SP": "1.08D0",    # [OPT] was 1.10; smoother low-nu spacing
            "AMP_FAC": "1.04D0",    # [OPT] was 1.05; smoother high-nu spacing
            "MAX_BF": "0.05D0",     # [OPT] was 0.10; tighter BF edge resolution
            "CROSS": "F",          # Sparse chi/eta with interpolation
            "V_CROSS": "500.0D0",   # 500 km/s for 500 km/s feature resolution
            "SIG_GAU_KMS": "500.0D0",  # match V_CROSS for consistency

            "USE_JREL": "T",
            "USE_J_REL": "T",
            "INCL_REL": "T",
            "INCL_ADV_TRANS": "T",
            # Keep monotonic J for stability during Lambda settling
            "USE_FRM_REL": "T",
            "NMONO_J": "F",
            "JC_W_EDD": "T",
            # Electron scattering: lambda iteration for J (Type IIn)
            "USE_LAM_ES": "F",
            "COH_ES": "T",
            "MIX_COH": "F",
            "ES_FAC": "0.6D0",       # Increased damping for ES oscillation

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "75.0D0",      # [OPT] Elevated turbulence for shock regions (not in STELLA T)
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)

            # Dessart+2015: T fixed to hydro; Ne free to converge self-consistently
            "FIX_T": "T",
            "FIX_NE": "F",            # Ne free in stable phase (populations from init are now reasonable)
            "FIX_IMP": "T",
            "FIX_INB_T": "F",
            "FIX_X_DPTH": "0",        # T fixed at ALL depths

            # G_ONLY/DIAG: consistent with NUM_BNDS=1 (Dessart+2015 style)
            # NOTE: With DIAG (no depth coupling), BA must be VERY accurate before freezing
            # or the uncorrected errors cause runaway oscillations (10^8% divergences).
            "SCALE_OPT": "MAJOR",
            "MAX_LIN": "2.0D0",       # tighter than 3.0; reduces LIMIT-clipping
            "MAX_CHNG": "1.0D4",      # finite hard cap
            "COMP_BA": "T",
            "N_FIX_BA": "50",        # many BA recomputes (was 5; need >25 for DIAG stability)
            "BA_CHK_FAC": "1.0D-05", # check BA every iteration when frozen (tighter than 5e-4)
            "FIX_BA": "5.0D-04",     # freeze only at < 0.05% error (was 5%; 100x tighter for DIAG)
            "LAM_VAL": "600.0D0",
            "NUM_LAM": "2",

            # [OPT] Bound-bound & blanketing for IIn models
            "GLOBAL_LINE": "SOB",   # Sobolev rates
            "GF_CUT": "5.0D-3",      # [OPT] was 1e-3; remove weaker lines
            "GF_LEV_CUT": "9",       # Apply GF_CUT above this level
            "AT_CUT": "22",          # Min atomic number for GF cut
            "MIN_TRANS": "15",       # [OPT] was 10; fewer marginal levels
            "THK_LINE": "F",         # Thick-line outer BC off
            "CHK_L_POS": "T",        # Prevent negative total opacity
            "NEG_OPAC_OPT": "SRCE_CHK", # Limit source when chi negative
            "SCL_LN": "T",           # Ignore lines in continuum opacity
            "SCL_LN_FAC": "0.1D0",   # Line scale-out limit (Doppler frac)
            "MIX_COH": "F",          # [OPT] was T; disable for convergence phase

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "75.0D0",      # [OPT] Elevated turbulence for shock regions (not in STELLA T)
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)

            # Use interpolated populations from init (not LTE)
            "LTE_EST": "F",
            "LIN_INT": "T",
            "THK_CONT": "F",

            "REV_RGRID": "T",       # Enable automatic grid revision for problematic regions
            "RD_IN_R_GRID": "F",
            "N_IB_INS": "15",        # More boundary points for DIAG solver stability (was 10; CMFGEN limit is 16)
            "N_OB_INS": "15",        # More boundary points for DIAG solver stability (was 10; CMFGEN limit is 16)
            "OB_RAT": "1.5D0",       # Outer boundary step ratio for grid stability
            "IB_RAT": "2.0D0",       # Inner boundary step ratio for grid stability
            # NOTE: ADJUST_R_DEFAULTS file is created separately, not a VADAT keyword
            # Use G_ONLY (NUM_BNDS=1, SOL_METH=DIAG)
            "N_TYPE": "G_ONLY",      # Changed from N_ON_J due to convergence issues (excessive iteration count with NaN errors)
            "POP_SCALE": "T",        # SNA stable setting
            "GREY_TAU": "10.0D0",     # SNA stable setting
            "METHOD": "ZERO",        # SNA stable setting
            "INC_RAY": "F",         # SNA stable setting
            "T_EXC": "3.00D-1",      # 3000 K unified floor
            "T_MIN": "3.0D-01",      # 3000 K unified floor

            "INC_GRID": "F",        # Disable automatic grid revision
            "VFRAC_FG": "0.50D0",    # [OPT] was 0.4; coarser ray spacing for SN    # Finer ray spacing for J stability
            "VFRAC_MOM": "0.50D0",   # [OPT] was 0.4; coarser moment spacing for SN   # Finer moment spacing for J stability
            "OB_METH": "HONJ",   # Use HONJ for outer boundary (required by MOM_JREL_V9)

            # [OPT] Spatial grid settings for IIn models
            "ALL_FREQ": "T",       # Apply accuracy at all nu
            "ACC_END": "5.0D0",     # Extend enhanced accuracy to more frequencies
            "N_INS": "2",            # Extra depths (if INC_GRID=T)
            "ST_INT": "1",          # First depth for insertion
            "END_INT": "30",         # Last depth for insertion
            "ND_QUAD": "30",         # [OPT] was 50; fewer inner quadrature depths
            "INTERP_TYPE": "LOG",  # Log R interpolation
            "N_PAR": "150",          # [OPT] was 300; fewer BA update intervals

            # [OPT] Line broadening for high-VINF IIn models
            "MAX_DOP": "8.0D0",         # [OPT] was 6.0; wider resonance zone for strong line overlap
            "FRAC_DOP": "1.5D0",        # [OPT] was 1.0; coarser in-resonance sampling
            "dV_CMF_PROF": "80.0D0",    # [OPT] was 200; resolve 100 km/s features
            "dV_CMF_WING": "200.0D0",   # [OPT] was 300; resolve 100 km/s features
            "ES_WING_EXT": "12000.0D0", # [OPT] was 2500; MUST exceed VINF=10472 km/s
            "R_CMF_WING_EXT": "4.0D0",  # [OPT] was 3.0; wider coherent ES wing

            # [OPT] Observed-frame spectrum for high-VINF IIn
            "OBS_EXT_RAT": "1.3D0",   # [OPT] was 1.1; wider observed-frame coverage
            "dV_OBS_PROF": "40.0D0",  # [OPT] was 100; resolve 100 km/s features
            "dV_OBS_WING": "150.0D0", # [OPT] was 400; resolve 100 km/s features
            "dV_OBS_BIG": "800.0D0",  # [OPT] was 2000; resolve 100 km/s features

            # Convergence settings for stable phase
            "EPS_TERM": "5.0D-02",
            "ACC_F": "1.0D-02",       # tighter Eddington convergence reduces late oscillations
            "ACC_END": "5.0D0",

            # Stay in pure Lambda; high LAM_VAL prevents Non-Lambda switch
            "DO_LAM_IT": "T",

            # Ng disabled during Lambda settling phase to avoid amplifying oscillations
            "DO_NG": "F",
            "IBEG_NG": "50",
            "BW_NG": "30",
            "ITS/NG": "4",
        }

    elif profile == "final":
        # Use all INIT parameters except: LTE_EST=F, MAX_CHNG=1.0D4, DO_NG=T, DO_LAM_AUTO=T
        # Reads population levels from *_IN files (from previous run)
        desired = {
            "SN_T_OPT": "USE_T_IN",
            "LSTAR": "0.0D0",         # SN: 0 (no central L; heating from hydro T/ED, not CMFGEN rad decays)
            "INC_RAD_DECAYS": inc_rad_decays,
            "INC_SHOCK_POWER": "F",
            "ALL_FREQ": "T",

            # [OPT] Continuum frequency grid for IIn models
            "RD_CF_FILE": "F",      # Build frequency grid internally
            "MIN_CF": "1.0D-1",        # Raised to 0.1 to avoid MOM_J_REL_V9 NaN errors at extreme far-IR
            "MAX_CF": "13.5D0",
            "FRAC_SP": "1.08D0",    # [OPT] was 1.10; smoother low-nu spacing
            "AMP_FAC": "1.04D0",    # [OPT] was 1.05; smoother high-nu spacing
            "MAX_BF": "0.05D0",     # [OPT] was 0.10; tighter BF edge resolution
            "CROSS": "F",          # Sparse chi/eta with interpolation
            "V_CROSS": "500.0D0",   # 500 km/s for 500 km/s feature resolution
            "SIG_GAU_KMS": "500.0D0",  # match V_CROSS for consistency

            "USE_JREL": "T",
            "USE_J_REL": "T",
            "INCL_REL": "T",
            "INCL_ADV_TRANS": "T",
            "USE_FRM_REL": "T",
            "NMONO_J": "F",
            "JC_W_EDD": "T",
            "THK_CONT": "F",
            "USE_LAM_ES": "F",
            "COH_ES": "T",
            "MIX_COH": "F",
            "ES_FAC": "0.6D0",

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "75.0D0",      # [OPT] Elevated turbulence for shock regions (not in STELLA T)
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)

            # Dessart+2015: T fixed to hydro; Ne free in final
            "FIX_T": "T",
            "FIX_NE": "F",            # Ne free in final (populations now well-established)
            "FIX_IMP": "T",
            "FIX_INB_T": "F",
            "FIX_X_DPTH": "0",        # T fixed at ALL depths

            # G_ONLY/DIAG: consistent with NUM_BNDS=1
            # NOTE: With DIAG (no depth coupling), BA must be VERY accurate before freezing
            # or the uncorrected errors cause runaway oscillations (10^8% divergences).
            "SCALE_OPT": "MAJOR",
            "MAX_LIN": "2.0D0",
            "MAX_CHNG": "1.0D4",
            "COMP_BA": "T",
            "N_FIX_BA": "50",        # many BA recomputes (was 5; need >25 for DIAG stability)
            "BA_CHK_FAC": "1.0D-05", # check BA every iteration when frozen (tighter than 5e-4)
            "FIX_BA": "5.0D-04",     # freeze only at < 0.05% error (was 5%; 100x tighter for DIAG)
            "LAM_VAL": "600.0D0",
            "NUM_LAM": "2",

            "GLOBAL_LINE": "SOB",
            "GF_CUT": "5.0D-3",
            "GF_LEV_CUT": "9",
            "AT_CUT": "22",
            "MIN_TRANS": "15",
            "THK_LINE": "F",
            "CHK_L_POS": "T",
            "NEG_OPAC_OPT": "SRCE_CHK",
            "SCL_LN": "T",
            "SCL_LN_FAC": "0.1D0",
            "MIX_COH": "F",
            "ES_FAC": "0.6D0",        # Dessart+2015: 0.6 throughout for IIn

            # [OPT] IIn CSM interaction settings (critical even with FIX_T from STELLA)
            "INC_CHG": "F",          # [OPT] was T; charge exchange disabled (CHG_EXCH_DATA not available)
            "VTURB": "75.0D0",      # [OPT] Elevated turbulence for shock regions (not in STELLA T)
            "ALLOW_OL": "T",         # [OPT] Line overlap for dense CSM (geometric, not thermal)

            "LTE_EST": "F",           # read population levels from *_IN files
            "IT_ON_T": "F",
            "POP_SCALE": "T",
            "T_MIN": "3.0D-01",
            "T_EXC": "3.00D-1",
            "GREY_TAU": "10.0D0",
            "METHOD": "ZERO",
            "INC_RAY": "F",

            "REV_RGRID": "T",
            "RD_IN_R_GRID": "F",
            "N_IB_INS": "15",         # More boundary points for DIAG solver stability (was 10; CMFGEN limit is 16)
            "N_OB_INS": "15",
            "OB_RAT": "1.5D0",
            "IB_RAT": "2.0D0",
            "N_TYPE": "G_ONLY",      # Changed from N_ON_J due to convergence issues (excessive iteration count with NaN errors)

            "INC_GRID": "F",        # Disabled to avoid SIGMA > -1 limitation in EXTEND_VTSIGMA for interacting SNe
            "VFRAC_FG": "0.50D0",    # [OPT] was 0.4; coarser ray spacing for SN
            "VFRAC_MOM": "0.50D0",   # [OPT] was 0.4; coarser moment spacing for SN
            "OB_METH": "HONJ",

            "ALL_FREQ": "T",
            "ACC_END": "5.0D0",
            "N_INS": "2",
            "ST_INT": "1",
            "END_INT": "30",
            "ND_QUAD": "30",
            "INTERP_TYPE": "LOG",
            "N_PAR": "150",

            "MAX_DOP": "8.0D0",
            "FRAC_DOP": "1.5D0",
            "dV_CMF_PROF": "80.0D0",
            "dV_CMF_WING": "200.0D0",
            "ES_WING_EXT": "12000.0D0",
            "R_CMF_WING_EXT": "4.0D0",

            "OBS_EXT_RAT": "1.3D0",
            "dV_OBS_PROF": "40.0D0",
            "dV_OBS_WING": "150.0D0",
            "dV_OBS_BIG": "800.0D0",

            "EPS_TERM": "1.0D-02",    # tightest criterion for final convergence
            "ACC_F": "1.0D-02",       # tighter Eddington convergence

            "DO_LAM_IT": "T",

            # Enable Ng acceleration in final phase once populations are settled
            "DO_NG": "T",
            "IBEG_NG": "50",          # start Ng after 50 Lambda iterations
            "BW_NG": "30",
            "ITS/NG": "4",
        }

    if convergence_mode == "stable":
        # Conservative damping for NLTE-heavy IIn models where LIMIT clipping repeats
        # and lambda corrections remain large for many iterations (Dessart+2015 style).
        # These overrides apply on top of the profile defaults for extra damping.
        # Stable-mode overrides: tighter step sizes and more Lambda iterations.
        # CRITICAL: N_FIX_BA must be large (>=30) and FIX_BA must be small (<=1e-2)
        # to keep CMFGEN in lambda-iteration until populations are well-converged.
        # Switching to BA too early (N_FIX_BA=5, FIX_BA=15%) causes the
        # "maximum change > 1" divergence in SOLVEBA_V13 because the Newton step
        # is taken while corrections are still O(10-100%), far outside the linear regime.
        stable_overrides = {
            "single_final": {
                "MAX_LIN": "1.0D0",     # tighter step size for LTE bootstrap stability
                "MAX_LAM": "1.0D5",     # 1e5 cap; increased to allow large LTE→NLTE corrections during cold-start
                "LAM_VAL": "2.0D2",
                "NUM_LAM": "40",
                "N_FIX_BA": "30",       # stay in lambda until corrections are small
                "FIX_BA": "1.0D-02",   # switch to BA only when max change < 1%
                "MAX_CHNG": "1.0D4",
                "IBEG_NG": "80",        # delay Ng start to iteration 80
                "SCALE_OPT": "MAJOR",
            },
            "init": {
                "MAX_LIN": "1.0D0",     # tighter than 2.0 to reduce LIMIT-clipped corrections
                "MAX_LAM": "1.0D5",     # 1e5 cap; increased to allow large LTE→NLTE corrections during cold-start
                "LAM_VAL": "2.0D2",     # stay in Lambda longer on init
                "NUM_LAM": "40",        # more Lambda iterations before switch
                "N_FIX_BA": "30",       # stay in lambda until corrections are small
                "FIX_BA": "1.0D-02",   # switch to BA only when max change < 1%
                "MAX_CHNG": "1.0D4",
                "IBEG_NG": "80",
                "SCALE_OPT": "MAJOR",
            },
            "stable": {
                "MAX_LIN": "1.5D0",     # tighter than 2.0; cuts LIMIT-clipping frequency
                "MAX_LAM": "1.0D3",     # 1000x; populations partially converged from init checkpoint
                "LAM_VAL": "1.0D3",     # stay in Lambda even longer for shock-rich models
                "NUM_LAM": "50",        # more Lambda iterations
                "N_FIX_BA": "30",       # stay in lambda until corrections are small
                "FIX_BA": "1.0D-02",   # switch to BA only when max change < 1%
                "IBEG_NG": "100",       # delay Ng further; populations still settling
                "SCALE_OPT": "MAJOR",
            },
            "final": {
                "MAX_LIN": "1.5D0",     # moderately tight for final convergence
                "MAX_LAM": "1.0D2",     # 100x cap; reduced to prevent population instabilities during final convergence
                "LAM_VAL": "2.0D2",     # keep Lambda fallback threshold
                "NUM_LAM": "10",        # moderate Lambda iterations
                "N_FIX_BA": "30",       # stay in lambda until corrections are small
                "FIX_BA": "1.0D-02",   # switch to BA only when max change < 1%
                "IBEG_NG": "50",        # start Ng after 50 Lambda iterations
                "SCALE_OPT": "MAJOR",
            },
        }
        if profile in stable_overrides:
            desired.update(stable_overrides[profile])

    # NUM_BNDS and SOL_METH must match N_TYPE (MODEL_SPEC and VADAT consistent)
    n_type = desired.get("N_TYPE", "G_ONLY")
    num_bnds = _num_bnds_for_n_type(n_type)
    sol_meth = "DIAG" if num_bnds == 1 else "TRIDIAG"
    desired["SOL_METH"] = sol_meth

    before = {k: None for k in desired}
    key_idx = find_existing_keys(lines)
    for k in desired:
        if k in key_idx:
            raw = _strip_inline_comment(lines[key_idx[k]]).strip()
            m = _VADAT_LINE_RE.match(raw)
            before[k] = m.group(2).strip() if m else None

    for k, v in desired.items():
        set_or_append_keyword(lines, k, v)

    changes: List[str] = []
    for k, v in desired.items():
        old = before[k]
        if old is None:
            changes.append(f"{k}: (added) -> {v}")
        elif old != v:
            changes.append(f"{k}: {old} -> {v}")
    if convergence_mode == "stable":
        changes.append("Convergence mode: stable (extra damping and longer lambda-phase enabled)")
    return changes, num_bnds


# ---------- OUT to IN (same as cur_cmf/com/out2in.sh) ----------
def out2in(wdir: Path) -> List[Tuple[str, str]]:
    """
    Rename *OUT files to *_IN in wdir (population levels from converged run).
    Matches cur_cmf/com/out2in.sh: each name ending with OUT becomes *_IN.
    Returns list of (old_name, new_name) for reporting.
    """
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
            f.rename(dest)
            renames.append((f.name, in_name))
        except OSError:
            pass
    return renames


# ---------- NEW base selection + rename behavior ----------
def select_base_vadat_path(wdir: Path, profile: str, vadat_out: Path, vadat_init: Path) -> Path:
    """
    Requested behavior:
      init          -> base is VADAT_INIT if exists else VADAT_INIT.bak if exists (else error)
      single_final  -> same as init (LTE cold start; no prior VADAT needed)
      stable        -> base is existing VADAT (else error)
      final         -> base is existing VADAT (else error)
    """
    if profile in ("init", "single_final"):
        # prefer VADAT_INIT, else fall back to VADAT_INIT.bak
        if vadat_init.exists():
            return vadat_init
        vadat_init_bak = vadat_init.with_name(vadat_init.name + ".bak")
        if vadat_init_bak.exists():
            return vadat_init_bak
        raise FileNotFoundError(f"Missing {vadat_init.name} and {vadat_init_bak.name} in {wdir}")
    else:
        # stable/final: use existing VADAT only
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


# ---------- main ----------
def main() -> int:
    p = argparse.ArgumentParser(
        description="Update CMFGEN VADAT from SN_HYDRO_DATA and write VADAT (clean)."
    )
    p.add_argument("directory", nargs="?", default=".", help="Model directory (default: current).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--init", action="store_true", help="LTE initial model: base VADAT_INIT, write VADAT with init profile (LTE_EST=T, T fixed, loose convergence).")
    g.add_argument("--stable", action="store_true", help="Stable Lambda-only run: base existing VADAT, pure Lambda iterations to settle populations before final (LTE_EST=F, LIN_INT=T).")
    g.add_argument("--final", action="store_true", help="Final run: base existing VADAT (after init/stable), apply final profile (LTE_EST=F, LIN_INT=T, tight convergence).")
    g.add_argument("--single-final", action="store_true", help="Single-final LTE cold-start: base VADAT_INIT, write VADAT with single_final profile (LTE_EST=T, DO_LAM_AUTO=T, no prior INIT needed).")

    p.add_argument("--vadat-init", default="VADAT_INIT", help="Template VADAT init file base name (default: VADAT_INIT).")
    p.add_argument("--vadat-out", default="VADAT", help="Output VADAT filename (default: VADAT).")
    p.add_argument("--hydro", default="SN_HYDRO_DATA", help="Hydro file name (default: SN_HYDRO_DATA).")
    p.add_argument("--lum-xg", default="lum.xg", help="lum.xg path relative to directory for LSTAR (default: lum.xg); set empty to disable.")
    p.add_argument(
        "--keep-lstar-zero",
        action="store_true",
        help="Keep LSTAR=0 when base/template has LSTAR=0, even if SN_HYDRO_DATA contains Energy Lum.",
    )
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
        "--thin-outer-bc",
        action="store_true",
        help="Treat outer boundary as optically thin: set THK_CONT=F in VADAT for this profile.",
    )
    p.add_argument(
        "--fix_ne",
        action="store_true",
        help="Fix electron density (Ne) in init/stable profiles. Default: allow Ne to vary.",
    )
    p.add_argument(
        "--vadat-overrides",
        default=None,
        help="Path to YAML config file containing vadat_overrides section to apply after auto-configuration.",
    )
    p.add_argument(
        "--vel-law",
        type=int,
        default=12,
        help="CMFGEN VEL_LAW value (default: 12 for IIn tau-based grid; 11 for direct SN_HYDRO_DATA R/V/SIGMA).",
    )
    p.add_argument(
        "--no-relativistic",
        action="store_true",
        help="Disable relativistic moment solver and corrections for faster debug runs. Sets USE_JREL=F, USE_J_REL=F, INCL_REL=F, INCL_ADV_TRANS=F, USE_FRM_REL=F.",
    )

    args = p.parse_args()
    wdir = Path(args.directory).resolve()
    if args.init:
        profile = "init"
    elif args.stable:
        profile = "stable"
    elif args.single_final:
        profile = "single_final"
    else:
        profile = "final"

    # Load VADAT overrides from YAML config if provided
    vadat_overrides = {}
    if args.vadat_overrides:
        try:
            import yaml
            config_path = Path(args.vadat_overrides)
            if config_path.exists():
                with open(config_path, 'r') as f:
                    config = yaml.safe_load(f)
                if config and 'cmfgen' in config and 'vadat_overrides' in config['cmfgen']:
                    vadat_overrides = config['cmfgen']['vadat_overrides']
                    print(f"  Loaded vadat_overrides from {config_path}: {list(vadat_overrides.keys())}")
        except Exception as e:
            print(f"  WARNING: Failed to load vadat_overrides from {args.vadat_overrides}: {e}")

    # Tee stdout to OUTGEN so prep output appears in the same log as CMFGEN
    outgen_path = wdir / "OUTGEN"
    outgen_file = None
    old_stdout = sys.stdout
    if not args.no_outgen:
        try:
            outgen_file = open(outgen_path, "a")
            outgen_file.write(
                f"\n=== prepare_VADAT.py ({profile}) {datetime.now().isoformat()} ===\n\n"
            )
            outgen_file.flush()

            class _Tee:
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

            sys.stdout = _Tee(old_stdout, outgen_file)
        except OSError:
            pass

    def _restore_stdout():
        if outgen_file is not None:
            sys.stdout = old_stdout
            outgen_file.close()

    try:
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

        # For stable/final: rename *OUT -> *_IN (previous run population levels become IN for next run)
        # single_final starts from LTE so no *OUT files exist yet; skip out2in.
        out2in_renames: List[Tuple[str, str]] = []
        if profile in ("stable", "final") and not args.no_out2in:
            out2in_renames = out2in(wdir)
            if out2in_renames:
                src = "init" if profile == "stable" else "stable/init"
                print(f"OUT -> IN (population files from {src} run):")
                for old, new in out2in_renames:
                    print(f"  {old} -> {new}")
                print()

        hydro = parse_sn_hydro(hydro_path)
        
        lines = read_lines(base_vadat_path)

        # Type IIn: if base VADAT has LSTAR=0, keep it (do not overwrite with Energy Lum)
        base_lstar_tok = _get_key_token(lines, "LSTAR")
        base_lstar_zero = False
        if base_lstar_tok is not None:
            try:
                base_lstar_zero = abs(float(base_lstar_tok.replace("D", "E").replace("d", "e"))) < 1.0e-30
            except (ValueError, TypeError):
                pass

        # Photospheric RSTAR/RMAX from Energy Lum path: skip if user forces LSTAR=0 from template.
        if not args.no_lum_rst_geometry and not (args.keep_lstar_zero and base_lstar_zero):
            maybe_set_rst_from_energy_lum(hydro, tau_phot=TAU_PHOT)

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
        # Init: use T at mid-depth so outer (NLTE) half gets excitation floor, not full LTE.
        t_exc = recommended_t_exc_init(hydro) if profile in ("init", "single_final") else recommended_t_exc(hydro)
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
            abund_changes = update_abundances_mass_fraction(lines, ave_mf)
        # SET_ABUND_CLUMP requires sum of mass fractions = 1; renormalize in case template + hydro != 1
        renorm_changes = renormalize_vadat_mass_fractions_to_unity(lines)
        abund_changes.extend(renorm_changes)

        n_iso_hydro = sn_hydro_header_num_isotopes(hydro_path)
        n_iso = max(int(hydro.params.get("num_isotopes", 0) or 0), n_iso_hydro)
        profile_changes, num_bnds = apply_profile(
            lines, profile, convergence_mode=args.convergence_mode, num_isotopes=n_iso
        )
        
        # Apply VADAT overrides from config (if provided)
        if vadat_overrides:
            # Validate SOL_METH/N_TYPE consistency before applying
            if "SOL_METH" in vadat_overrides:
                n_type_in_lines = _get_key_token(lines, "N_TYPE") or "G_ONLY"
                required_sol_meth = "DIAG" if _num_bnds_for_n_type(n_type_in_lines) == 1 else "TRIDIAG"
                requested_sol_meth = str(vadat_overrides["SOL_METH"]).strip().upper()
                if requested_sol_meth != required_sol_meth:
                    print(
                        f"  WARNING: vadat_overrides SOL_METH={requested_sol_meth} is inconsistent with "
                        f"N_TYPE={n_type_in_lines} (requires {required_sol_meth}). "
                        f"Override ignored to prevent CMFGEN crash."
                    )
                    vadat_overrides = {k: v for k, v in vadat_overrides.items() if k != "SOL_METH"}
            override_changes = update_keywords_in_lines(lines, vadat_overrides, numeric_format_like_old=True)
            if override_changes:
                print(f"  Applied vadat_overrides:")
                for change in override_changes:
                    print(f"    {change}")
                profile_changes.extend(override_changes)
        # RD_SN_DATA matches isotope MF lines to ISO(:) from NUC_DECAY_DATA; that table is only
        # loaded when INC_RAD_DECAYS=T. Force T from the file header so a missed parse cannot
        # leave INC_RAD_DECAYS=F (which yields NUM_ISOTOPES=0 and "isotope data not recognized").
        # However, if CMFGEN_INC_RAD_DECAYS is set to 0, force F to avoid isotope consistency issues.
        inc_rad_decays_env = os.environ.get("CMFGEN_INC_RAD_DECAYS", "1")
        if n_iso_hydro > 0 and inc_rad_decays_env == "1":
            set_or_append_keyword(lines, "INC_RAD_DECAYS", "T")
            profile_changes.append(
                f"INC_RAD_DECAYS=T (SN_HYDRO_DATA declares {n_iso_hydro} isotope block(s); "
                "required for RD_SN_DATA / NUC_DECAY_DATA ISO table)"
            )
        elif n_iso_hydro > 0 and inc_rad_decays_env == "0":
            set_or_append_keyword(lines, "INC_RAD_DECAYS", "F")
            profile_changes.append(
                f"INC_RAD_DECAYS=F (forced by CMFGEN_INC_RAD_DECAYS=0 to avoid isotope consistency issues)"
            )

        # Set VEL_LAW from CLI argument (default: 12 for IIn models with tau-based grid)
        set_or_append_keyword(lines, "VEL_LAW", str(args.vel_law))
        profile_changes.append(f"VEL_LAW={args.vel_law} (from --vel-law CLI argument)")

        # Disable relativistic solvers for faster debug runs (--no-relativistic)
        if args.no_relativistic:
            rel_keys = ["USE_JREL", "USE_J_REL", "INCL_REL", "INCL_ADV_TRANS", "USE_FRM_REL"]
            for key in rel_keys:
                set_or_append_keyword(lines, key, "F")
            profile_changes.append(
                "Relativistic solvers disabled (--no-relativistic): "
                f"{', '.join(rel_keys)}=F"
            )

        model_spec_path = wdir / args.model_spec
        die_key_changes = ensure_die_keys_for_model_spec(lines, model_spec_path)
        profile_changes.extend(die_key_changes)

        # Init/Stable: ensure LAMBDA is enabled (DO_LAM_IT=T) and allow very large MAXCH so we don't abort.
        # DO_LAM_AUTO is controlled via IN_ITS (init/stable: F, final: T).
        # MAX_LAM caps the per-step Lambda population change: too high causes population explosions
        # at the cold-start outer boundary (STEQ ~ 1e9), which inflates populations by MAX_LAM-fold
        # and crashes STEQ in the next iteration with a SIGSEGV. Values per profile:
        #   init / single_final: 1e5 (cold-start, farthest from NLTE solution - increased to allow large LTE→NLTE corrections)
        #   stable:              1e3 (partially converged from init checkpoint)
        #   final:               1e2 (well-converged from stable/init checkpoint - reduced for stability)
        if profile == "init":
            set_or_append_keyword(lines, "DO_LAM_IT", "T")
            set_or_append_keyword(lines, "SCALE_OPT", "MAJOR")
            set_or_append_keyword(lines, "MAX_LAM", "1.0D5")
            set_or_append_keyword(lines, "MAX_CHNG", "1.0D30")
            set_or_append_keyword(lines, "LIN_INT", "F")
            # Keep RD_IN_R_GRID=T to use RDINR, but LIN_INT=F to avoid HI_IN read error
            profile_changes.append(
                "Init profile: DO_LAM_IT=T, SCALE_OPT=MAJOR, MAX_LAM=1.0D5, MAX_CHNG=1.0D30, LIN_INT=F  "
                "(LAMBDA iters; MAX_LAM=1e5 allows large LTE→NLTE corrections during cold-start)"
            )
        elif profile == "single_final":
            set_or_append_keyword(lines, "DO_LAM_IT", "T")
            set_or_append_keyword(lines, "SCALE_OPT", "MAJOR")
            set_or_append_keyword(lines, "MAX_LAM", "1.0D5")
            set_or_append_keyword(lines, "MAX_CHNG", "1.0D30")
            profile_changes.append(
                "single_final profile: DO_LAM_IT=T, SCALE_OPT=MAJOR, MAX_LAM=1.0D5, MAX_CHNG=1.0D30  "
                "(LAMBDA iters; MAX_LAM=1e5 allows large LTE→NLTE corrections during cold-start)"
            )
        elif profile == "stable":
            set_or_append_keyword(lines, "DO_LAM_IT", "T")
            set_or_append_keyword(lines, "SCALE_OPT", "MAJOR")
            set_or_append_keyword(lines, "MAX_LAM", "1.0D3")
            set_or_append_keyword(lines, "MAX_CHNG", "1.0D30")
            profile_changes.append(
                "stable profile: DO_LAM_IT=T, SCALE_OPT=MAJOR, MAX_LAM=1.0D3, MAX_CHNG=1.0D30  "
                "(LAMBDA iters; MAX_LAM=1000 for partially-converged populations from init)"
            )
        elif profile == "final":
            set_or_append_keyword(lines, "SCALE_OPT", "MAJOR")
            set_or_append_keyword(lines, "MAX_LAM", "1.0D2")
            set_or_append_keyword(lines, "MAX_CHNG", "1.0D30")
            profile_changes.append(
                "Final profile: SCALE_OPT=MAJOR, MAX_LAM=1.0D2, MAX_CHNG=1.0D30  "
                "(MAX_LAM=100 for stable final convergence with well-converged populations)"
            )
        if profile == "final" and out2in_renames:
            profile_changes.append(f"OUT->IN: {len(out2in_renames)} files (*OUT renamed to *_IN)")

        # Temperature control (both --init and --final):
        # Fix T only in the inner region (depths with tau_es >= tau_phot), i.e.
        # FIX_T=T, FIX_INB_T=F, FIX_X_DPTH=fix_x_dpth. Outer layers are allowed
        # to adjust (important for optically thin CSM).
        nd = int(hydro.params["ND"])
        tau_arr, fix_x_dpth, i_phot = estimate_tau_es_and_fix_x_dpth(
            hydro, tau_phot=TAU_PHOT
        )
        tau_inner = float(tau_arr[0]) if tau_arr.size else 0.0
        tau_ob = float(tau_arr[-1]) if tau_arr.size else 0.0
        tau_phot_val = float(tau_arr[i_phot]) if tau_arr.size and 0 <= i_phot < tau_arr.size else 0.0

        # Avoid confusing diffs from apply_profile: we always override FIX_T/FIX_INB_T/FIX_INB_NE/FIX_X_DPTH here.
        profile_changes = [
            c
            for c in profile_changes
            if not (
                c.strip().startswith("FIX_T:")
                or c.strip().startswith("FIX_INB_T:")
                or c.strip().startswith("FIX_INB_NE:")
                or c.strip().startswith("FIX_X_DPTH:")
                or c.strip().startswith("TAU_SCL_T:")
            )
        ]

        # Force: T fixed to input SN_HYDRO_DATA at all depths. FIX_X_DPTH=0.
        set_or_append_keyword(lines, "FIX_T", "T")
        set_or_append_keyword(lines, "FIX_INB_T", "F")
        set_or_append_keyword(lines, "FIX_X_DPTH", "0")
        # Electron density control:
        # - Default: allow Ne to vary. With T fixed, hard-fixing Ne can make early iterations
        #   extremely stiff and prone to luminosity blow-ups / NaNs.
        # - Optional: --fix_ne forces Ne fixed for init/stable when desired.
        if profile in ("init", "stable", "single_final") and args.fix_ne:
            set_or_append_keyword(lines, "FIX_NE", "T")
            set_or_append_keyword(lines, "FIX_INB_NE", "T")
        else:
            set_or_append_keyword(lines, "FIX_NE", "F")
            set_or_append_keyword(lines, "FIX_INB_NE", "F")
        set_or_append_keyword(lines, "FIX_X_DPTH", "0")
        set_or_append_keyword(lines, "TAU_SCL_T", "0.0D0")
        if profile in ("init", "stable", "single_final"):
            label = "init" if profile == "init" else ("single_final" if profile == "single_final" else "stable")
            ne_mode = "T" if args.fix_ne else "F"
            profile_changes.append(
                f"FIX_T=T, FIX_NE={ne_mode}, FIX_INB_T=F, FIX_INB_NE={ne_mode}, FIX_X_DPTH=0  "
                f"(T fixed to SN_HYDRO_DATA; Ne {'fixed' if args.fix_ne else 'varies'} for {label})"
            )
            print(
                f"  Forced T/Ne fix [{profile}]: FIX_T=T, FIX_NE={ne_mode}, FIX_INB_T=F, FIX_INB_NE={ne_mode}, FIX_X_DPTH=0 "
                f"(T fixed to SN_HYDRO_DATA; Ne {'fixed' if args.fix_ne else 'varies'} for {label})"
            )
        else:  # final
            profile_changes.append(
                f"FIX_T=T, FIX_NE=F, FIX_INB_T=F, FIX_INB_NE=F, FIX_X_DPTH=0  "
                f"(T fixed to SN_HYDRO_DATA; Ne varies)"
            )
            print(
                f"  Forced T fix [{profile}]: FIX_T=T, FIX_X_DPTH=0 (T fixed to SN_HYDRO_DATA; Ne varies everywhere)"
            )
        print()

        # If requested, force a "thin outer boundary" treatment by setting THK_CONT=F
        # regardless of the default profile. This is useful when the hydro outer
        # layers are clearly optically thin and the thick-continuum assumption
        # causes BA / DTDR issues.
        if args.thin_outer_bc:
            prev_thk = _get_key_token(lines, "THK_CONT")
            set_or_append_keyword(lines, "THK_CONT", "F")
            if prev_thk is None:
                profile_changes.append("THK_CONT: (added) -> F  (thin outer boundary override)")
            elif prev_thk.strip().upper() != "F":
                profile_changes.append(f"THK_CONT: {prev_thk.strip()} -> F  (thin outer boundary override)")

        # LSTAR after profile so it is not overwritten by profile's 0.0D0.
        # By default, set LSTAR from SN_HYDRO_DATA (Energy Lum) when present.
        # Use --keep-lstar-zero to preserve template/base LSTAR=0.
        lstar_token: Optional[str] = None
        if args.keep_lstar_zero and base_lstar_zero:
            profile_changes = [c for c in profile_changes if not (c.strip().startswith("LSTAR:") and " -> " in c)]
            set_or_append_keyword(lines, "LSTAR", "0.0D0")
            profile_changes.append("LSTAR: kept 0 (--keep-lstar-zero and base/template had LSTAR=0)")
        elif "energy_lum_lsun" in hydro.params:
            profile_changes = [c for c in profile_changes if not (c.strip().startswith("LSTAR:") and " -> " in c)]
            val = float(hydro.params["energy_lum_lsun"])
            lstar_token = f"{val:.6E}".replace("E", "D")
            set_or_append_keyword(lines, "LSTAR", lstar_token)
            profile_changes.append(f"LSTAR: from SN_HYDRO_DATA Energy Lum [lsun] = {val:.6e} L_sun (inner boundary)")
        elif base_lstar_zero:
            # No Energy Lum in file; keep LSTAR=0 from base (Type IIn; power from hydro/decays only).
            profile_changes = [c for c in profile_changes if not (c.strip().startswith("LSTAR:") and " -> " in c)]
            profile_changes.append("LSTAR: kept 0 (base had LSTAR=0, no Energy Lum in SN_HYDRO_DATA)")
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
        elif vel_law_int not in (11, 12):
            print(f"WARNING: VEL_LAW={vel_law_int} is uncommon for SN_HYDRO_DATA-driven runs (expected 11 or 12).")

        if not args.no_write_rdinr:
            rdinr_path = wdir / "RDINR"
            write_rdinr(hydro, rdinr_path)
            profile_changes.append(f"Wrote {rdinr_path.name} (R grid from SN_HYDRO_DATA)")
            nd = int(hydro.params["ND"])
            r_oi = (hydro.radius_cm[::-1] / 1e10).astype(float, copy=True)
            v_oi = np.abs(hydro.vel_cms[::-1]) / 1e5
            print(f"RDINR written: ND={nd}  R(1)={r_oi[0]:.6e} R(ND)={r_oi[-1]:.6e} (1e10 cm)  "
                  f"V(1)={v_oi[0]:.2f} V(ND)={v_oi[-1]:.2f} (km/s)")
            # input_hydro.dat in SNIIN format is only relevant for VEL_LAW=12.
            if vel_law_int == 12:
                input_hydro_path = wdir / "input_hydro.dat"
                write_input_hydro_dat_sniin(hydro, input_hydro_path)
                profile_changes.append(f"Wrote {input_hydro_path.name} (VEL_LAW=12 RV_SN_MODEL_SNIIN)")
            else:
                profile_changes.append(
                    f"Skipped input_hydro.dat write (VEL_LAW={vel_law_int}; SNIIN format is for VEL_LAW=12)"
                )
            # T_IN and GAMMAS_IN for LTE initial conditions (VEL_LAW=12 + LTE_EST when FIN_CAL_GRID absent)
            if args.no_tin_gammas:
                print("  (T_IN/GAMMAS_IN skipped: --no-tin-gammas requested)")
            else:
                if hydro.temp_104k is not None and hydro.ed_cgs is not None:
                    write_t_in_from_hydro(hydro, wdir / "T_IN")
                    write_gammas_in_from_hydro(hydro, wdir / "GAMMAS_IN")
                    profile_changes.append("Wrote T_IN and GAMMAS_IN (LTE initial T/ED from SN_HYDRO_DATA)")
                else:
                    print("  (T_IN/GAMMAS_IN skipped: SN_HYDRO_DATA needs Temperature and Electron density sections)")
            print()

        # Outer boundary: if N_OB_INS/N_IB_INS are set > CMFGEN limit, cap them (profile uses 0 = no insertion)
        nd = int(hydro.params["ND"])
        # Allow large boundary insertion counts for SN models when requested.
        # Default cap was too low for ND=100 and N_*_INS=30.
        # Increased to 60 for J-stability at outer boundary in Type IIn models.
        N_OB_INS_MAX = 60
        max_ob_ins = N_OB_INS_MAX
        max_ib_ins = N_OB_INS_MAX
        for key, max_val in (("N_OB_INS", max_ob_ins), ("N_IB_INS", max_ib_ins)):
            tok = _get_key_token(lines, key)
            if tok is not None:
                try:
                    val = int(float(tok))
                    if val > max_val:
                        set_or_append_keyword(lines, key, str(max_val))
                        profile_changes.append(f"{key}: capped {val} -> {max_val}")
                except (ValueError, TypeError):
                    pass

        # # Always back up existing output VADAT before overwriting
        # bak = backup_if_exists(vadat_out_path, f"pre_{profile}")
        # if bak:
        #     print(f"Backed up existing {vadat_out_path.name} -> {bak.name}")

        clean = clean_vadat_lines(lines)
        write_lines(vadat_out_path, clean)

        # Keep IN_ITS consistent with VADAT for lambda control.
        in_its_changes = update_in_its_for_profile(wdir, profile)
        if in_its_changes:
            print(f"IN_ITS updates (IN_ITS):")
            for c in in_its_changes:
                print(f"  - {c}")
            print()

        # Create ADJUST_R_DEFAULTS for runtime grid refinement (handles problematic regions)
        adjust_r_path = wdir / "ADJUST_R_DEFAULTS"
        adjust_r_content = """! ADJUST_R_DEFAULTS - Runtime R-grid refinement for CMFGEN stability
! Targets the deep ejecta/shock interaction region where convergence
! problems concentrate (e.g. depth ~285-295 in a 311-point grid).
! Dessart+2015: grid refinement near shock front critical for IIn models.
!
[N_ITS] 5
! Number of grid revisions to perform
[STRT_ITS] 1
! First iteration to start grid revision
[FREQ_ITS] 10
! Frequency of grid revisions (every N iterations)
[GRID_METH] TAU_SPACE
! Regridding method: TAU_SPACE, VEL_SPACE, FULL_R, OB_ONLY
[GRID_TYPE] REFINE
! Grid type: UNIFORM, MODUN, FIX_NX, REFINE, VTAU
[FG_MIN] -2.0D0
! Min log(Tau) for fine grid: targets tau~0.01 (deep ejecta)
[FG_MAX] 0.5D0
! Max log(Tau) for fine grid: up to tau~3 (shock front region)
! Parameters for outer/inner boundary handling
[OB_OPT] DEFAULT
! Outer boundary option: DEFAULT or SPECIFY
[IB_OPT] DEFAULT
! Inner boundary option: DEFAULT or SPECIFY
"""
        adjust_r_path.write_text(adjust_r_content)
        print(f"Wrote ADJUST_R_DEFAULTS (runtime grid refinement for problematic regions)")
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

        print()
        print(f"OK: base={base_vadat_path.name} -> wrote {vadat_out_path.name} (no comments).")
        if profile == "init":
            print("  (Init: LTE initial populations; T allowed to vary at all depths; Ne from hydro.)")
        elif profile == "stable":
            print("  (Stable: pure Lambda settling; T allowed to vary at all depths; Ne varies.)")
        else:
            print("  (Final: LTE_EST=F, LIN_INT=T; T allowed to vary at all depths; Ne varies.)")
            if out2in_renames:
                print(f"  Renamed {len(out2in_renames)} *OUT -> *_IN (init populations as starting point).")
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

        # Update MODEL_SPEC so ND (and NP = ND+NC) match SN_HYDRO_DATA
        if model_spec_path.exists() and not args.no_model_spec:
            nd_hydro = int(hydro.params["ND"])
            model_spec_changes = update_model_spec(model_spec_path, nd_hydro, nc=None, num_bnds=num_bnds)
            if model_spec_changes:
                print()
                print(f"MODEL_SPEC updates ({model_spec_path.name}):")
                for c in model_spec_changes:
                    print("  -", c)
                if num_bnds == 1:
                    print("  (NUM_BNDS=1 uses DIAG solver; avoids DGETRF/band-solver failures with SN hydro.)")
                print()
            else:
                print(f"MODEL_SPEC: no updates (file must contain ND and NC keywords).")
        elif args.no_model_spec:
            print("MODEL_SPEC: skipped (--no-model-spec).")
        elif not model_spec_path.exists():
            print(f"MODEL_SPEC: {model_spec_path.name} not found; skipping.")

        return 0
    finally:
        _restore_stdout()


if __name__ == "__main__":
    sys.exit(main())
