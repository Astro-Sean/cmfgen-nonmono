#!/usr/bin/env python3
"""
Prepare the SN_B initial model, v2: power-law CSM (rho ~ r^-n) setup.

This is a variant of prepare_SN_B_model.py in which the 1/r^2 constant-
velocity "wind" between the inner shell and the old SN ejecta is replaced
by an analytical power-law circumstellar medium:

  * Density:  rho(r) = rho15 * (r / 1e15 cm)^-n,  set via --csm-powerlaw-n
    (n = 1, 2, 3 ...).  With density matching (default with ejecta) the
    law is anchored at the ejecta inner edge, so rho is continuous at
    the interface and the CSM is a smooth inward continuation of the
    SN_A envelope.  An optional sharp outer cutoff
    1 / (1 + (r/Rk)^cutoff_n) is available -- mainly useful in
    --no-ejecta mode, since with an ejecta component the old SN_A ejecta
    itself plays the role of the outer envelope.

  * Velocity: homologous expansion  u(r) = r / t_cs
    where t_cs is the kinematic age of the envelope.  An optional
    radiative acceleration term u_a = u0 * (Rs/r)^2 is available.

  * Density normalisation (three mutually exclusive options):
      - density-matched to the ejecta inner edge (default with ejecta):
        rho_csm(r_interface) = rho_ejecta(r_interface).  This makes the
        CSM a smooth inward continuation of the SN_A envelope.
      - --csm-rho15: rho at 1e15 cm in g/cc.
      - --csm-w:     w = 4 pi r^2 rho at 1e15 cm in g/cm.

  * Kinematic age t_cs: if --csm-age-days is 0 (default), t_cs is set to
    the SN_A ejecta age (source_age + t_start) so that the homologous CSM
    velocity joins the ejecta velocity continuously at the interface.
    Set it explicitly for a slower (younger-age) CSM.

Structure (inner -> outer):
  1. H/He shell (Dessart 2010jl flat-core + power-law envelope, same as v1)
  2. Shell-CSM buffer (log-linear density ramp, same as v1)
  3. Power-law CSM (rho ~ r^-n, homologous u)
  4. Thin contact-discontinuity ramp (v, T, composition; density continuous
     when density-matched)
  5. Old SN_A ejecta (the "envelope")

Velocity-driven convention identical to v1: EBurst=0 and EKO=0 are enforced
in the .dat file so STELLA preserves the constructed velocity field.

The CLI is drop-in compatible with prepare_SN_B_model.py: all v1 arguments
are accepted.  Wind-specific arguments (--wind-mdot, --wind-velocity,
--wind-duration-days, --wind-mdot-bg, --we-transition-*) are parsed but
ignored in favour of the --csm-* family; --wind-T doubles as the default
CSM temperature and --wind-r-outer-cm as the default CSM outer radius in
--no-ejecta mode.  --density-matched-wind maps onto CSM density matching.

Modes:
  default      : shell + power-law CSM + SN_A ejecta
  --no-ejecta  : shell + power-law CSM only (needs
                 --csm-r-outer-cm or --wind-r-outer-cm and an explicit
                 density normalisation)
  --no-wind    : shell in immediate contact with ejecta (identical to v1;
                 the CSM is skipped entirely)

Usage:
    python prepare_SN_B_model_v2.py <source_dir> <output_dir> [options]

    # SN_A ejecta as the envelope, density-matched CSM (default):
    python prepare_SN_B_model_v2.py SN_A SN_B_v2 --shell-composition solar

    # Example: ~1.5 yr kinematic age, explicit density:
    python prepare_SN_B_model_v2.py SN_A SN_B_v2 \\
        --csm-age-days 550 --csm-rho15 6e-15 --no-csm-density-match
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# Reuse the v1 machinery: constants, shell builder, ejecta reader,
# buffer/grid helpers, boundary cleaners, validators and file writers.
sys.path.insert(0, str(Path(__file__).parent))
from prepare_SN_B_model import (  # noqa: E402
    A_RAD,
    K_B,
    KAPPA_ES,
    M_H,
    MSUN,
    SOLAR_COMPOSITION,
    YR_TO_S,
    _build_optimized_shell_grid,
    _cap_buffer_zones,
    _shell_mass_factor,
    clean_combined_outer_boundary,
    combine_shell_ejecta,
    create_hhe_shell,
    plot_ejecta_smoothing,
    read_source_ejecta,
    validate_stella_grid,
    write_x_abn,
    write_x_dat,
    write_x_hyd,
    write_x_xni,
)

R_REF_1E15 = 1.0e15  # reference radius (cm)
M_P = 1.6726e-24  # proton mass (g), for n ~ rho/(mu*m_p) prints


# ----------------------------------------------------------------------
# Power-law CSM profile
# ----------------------------------------------------------------------


def _csm_cutoff_factor(csm, r):
    """Sharp outer cutoff factor 1/(1+(r/Rk)^n); 1 when disabled."""
    rk = csm.get("cutoff_r_cm", 0.0)
    if rk and rk > 0:
        return 1.0 / (1.0 + (np.asarray(r) / rk) ** csm["cutoff_n"])
    return np.ones_like(np.asarray(r, dtype=float))


def _csm_rho(csm, r):
    """CSM density law times optional cutoff.

    powerlaw: rho15 * (r/1e15)^slope.
    beta wind: rho = Mdot / (4 pi r^2 v(r)) with the beta velocity law.
    """
    r = np.asarray(r, dtype=float)
    if csm.get("law") == "beta":
        v = _csm_v(csm, r)
        rho = csm["mdot_gs"] / (4.0 * np.pi * r**2 * np.maximum(v, 1e-10))
        rho = rho * csm.get("clumping_D", 1.0)
    else:
        rho = csm["rho15"] * (r / R_REF_1E15) ** csm["slope"]
    return rho * _csm_cutoff_factor(csm, r)


def _csm_v(csm, r):
    """CSM velocity law.

    powerlaw: homologous u = r/t_cs plus optional radiative term
        u_a = u0 (Rs/r)^2, with Rs at the shell/CSM interface radius.
    beta wind: standard beta velocity law
        v = v_inf (1 - R*/r)^beta (Castor, Abbott & Klein 1975 form).
    """
    r = np.asarray(r, dtype=float)
    if csm.get("law") == "beta":
        fac = np.clip(1.0 - csm["rstar_cm"] / np.maximum(r, 1e-30), 1e-3, 1.0)
        return csm["vinf_cms"] * fac ** csm["beta"]
    v = r / csm["t_cs_s"]
    if csm.get("ua0_cms", 0.0) > 0:
        v = v + csm["ua0_cms"] * (csm["r_acc_cm"] / r) ** 2
    return v


def _csm_log_slope(csm, r):
    """Local d log rho / d log r of the CSM density law at r."""
    r = np.asarray(r, dtype=float)
    if csm.get("law") == "beta":
        # rho = Mdot/(4 pi r^2 v); dv/dr>0 makes the profile steeper
        # than r^-2 near R*.  dlog v/dlog r = beta (R*/r)/(1 - R*/r).
        fac = np.clip(1.0 - csm["rstar_cm"] / np.maximum(r, 1e-30), 1e-3, 1.0)
        return -2.0 - csm["beta"] * (csm["rstar_cm"] / r) / fac
    return np.full(np.asarray(r).shape, csm["slope"], dtype=float)


def _csm_T(csm, r):
    """CSM temperature: T15 * (rho/rho15)^T_gamma, floored at T_min.

    An optional temperature law uses T ~ rho^alpha
    (alpha ~ 0.32) throughout the star and the attached envelope.  The
    same law applied to the CSM keeps the dense part near the
    photoionized value T15 while the tenuous cutoff tail cools toward
    T_min instead of staying hot in near-vacuum.
    """
    gamma = csm.get("T_gamma", 0.0)
    if gamma == 0.0:
        T = np.full(np.asarray(r).shape, csm["T"], dtype=float)
    else:
        T = csm["T"] * (_csm_rho(csm, r) / csm["rho15"]) ** gamma
    return np.maximum(T, csm.get("T_min", 0.0))


def _csm_set_match(csm, r_target, rho_target):
    """Rescale the density normalisation so rho_csm(r_target)=rho_target.

    For the beta wind this solves for Mdot instead of rho15.
    """
    r_target = float(r_target)
    if csm.get("law") == "beta":
        v = float(_csm_v(csm, r_target))
        csm["mdot_gs"] = (
            rho_target * 4.0 * np.pi * r_target**2 * max(v, 1e-10)
            / csm.get("clumping_D", 1.0)
        )
        csm["norm_mode"] = "match"
        return
    law = (r_target / R_REF_1E15) ** csm["slope"] * float(
        _csm_cutoff_factor(csm, r_target)
    )
    csm["rho15"] = rho_target / max(law, 1e-300)


def _hermite_log_ramp(r, r_lo, r_hi, rho_lo, rho_hi, s_lo, s_hi):
    """C1-smooth density ramp: matches endpoint densities AND slopes.

    Hermite cubic in (log r, log rho) with endpoint log-log slopes
    s_lo, s_hi, so the buffer joins the shell envelope and the CSM
    power law without a visible corner.  Falls back to a log-linear
    ramp if the cubic would go non-monotonic.
    """
    x = np.log10(np.asarray(r, dtype=float))
    x0, x1 = np.log10(r_lo), np.log10(r_hi)
    dx = max(x1 - x0, 1e-30)
    t = np.clip((x - x0) / dx, 0.0, 1.0)
    y0 = np.log10(max(rho_lo, 1e-300))
    y1 = np.log10(max(rho_hi, 1e-300))
    m0, m1 = s_lo * dx, s_hi * dx
    t2, t3 = t * t, t * t * t
    y = (
        (2.0 * t3 - 3.0 * t2 + 1.0) * y0
        + (t3 - 2.0 * t2 + t) * m0
        + (-2.0 * t3 + 3.0 * t2) * y1
        + (t3 - t2) * m1
    )
    if len(y) > 1 and np.any(np.diff(y) * np.sign(y1 - y0) < -1e-10):
        y = np.linspace(y0, y1, len(y))
    return 10.0**y


def _shell_edge_slope(shell, nfit=5):
    """d log rho / d log r at the shell outer edge (envelope tail)."""
    nfit = min(nfit, len(shell["r"]))
    return float(
        np.polyfit(
            np.log10(shell["r"][-nfit:]),
            np.log10(np.maximum(shell["rho"][-nfit:], 1e-300)),
            1,
        )[0]
    )


def create_powerlaw_csm(
    r_inner_cm,
    r_outer_cm,
    slope=-0.4,
    age_days=None,
    rho15=None,
    w_gcm=None,
    T_K=5000.0,
    composition="solar",
    xh=0.7,
    yhe=0.3,
    ua0_kms=0.0,
    r_acc_cm=None,
    cutoff_r_cm=0.0,
    cutoff_n=12.0,
    T_gamma=0.0,
    T_min=0.0,
    match_shell=False,
    law="powerlaw",
    mdot_msun_yr=0.0,
    vinf_kms=0.0,
    rstar_cm=0.0,
    beta=1.0,
    clumping_D=1.0,
):
    """Create the power-law CSM profile.

    rho(r) = rho15 * (r/1e15)^slope, u(r) = r/t_cs, T = T_K isothermal.

    Normalisation: rho15 is the density at 1e15 cm BEFORE the optional
    cutoff factor.  Give either rho15 directly or w_gcm = 4 pi r^2 rho at
    1e15 cm; density matching to the ejecta is applied later by the
    combiner via _csm_set_match.

    Parameters:
        r_inner_cm, r_outer_cm: radial extent of the CSM region.
        slope:        density power-law index s (rho ~ r^s, s < 0).
        age_days:     kinematic age t_cs in days (homologous u = r/t_cs).
        rho15:        density at 1e15 cm in g/cc (None if matched later).
        w_gcm:        w = 4 pi r^2 rho at 1e15 cm in g/cm (alternative norm).
        T_K:          isothermal CSM temperature.
        composition:  'solar' or 'hhe'.
        ua0_kms:      radiative acceleration amplitude u0 (0 = off).
        r_acc_cm:     reference radius Rs for u_a = u0(Rs/r)^2
                      (default: r_inner_cm, i.e. the CDS radius).
        cutoff_r_cm:  sharp outer cutoff radius Rk (0 = no cutoff).
        cutoff_n:     cutoff steepness.
    """
    if age_days is None or age_days <= 0:
        raise ValueError(
            "create_powerlaw_csm requires age_days > 0 " "(kinematic age t_cs)."
        )

    csm = {
        "r_inner": r_inner_cm,
        "r_outer": r_outer_cm,
        "slope": slope,
        "t_cs_s": age_days * 86400.0,
        "rho15": 0.0,
        "T": T_K,
        "ua0_cms": ua0_kms * 1e5,
        "r_acc_cm": r_acc_cm if r_acc_cm else r_inner_cm,
        "cutoff_r_cm": cutoff_r_cm,
        "cutoff_n": cutoff_n,
        "composition": composition,
        "xh": xh,
        "yhe": yhe,
        "T_gamma": T_gamma,
        "T_min": T_min,
        "match_shell": match_shell,
        "law": law,
        "mdot_gs": mdot_msun_yr * MSUN / YR_TO_S,
        "vinf_cms": vinf_kms * 1e5,
        "rstar_cm": rstar_cm if rstar_cm > 0 else r_inner_cm,
        "beta": beta,
        "clumping_D": clumping_D,
    }

    # Normalisation
    if law == "beta":
        csm["norm_mode"] = "mdot" if csm["mdot_gs"] > 0 else "match"
    elif rho15 is not None and rho15 > 0:
        csm["rho15"] = rho15
        csm["norm_mode"] = "rho15"
    elif w_gcm is not None and w_gcm > 0:
        csm["rho15"] = w_gcm / (4.0 * np.pi * R_REF_1E15**2)
        csm["norm_mode"] = "w"
    else:
        csm["norm_mode"] = "match"  # set later via _csm_set_match

    return csm


def sample_powerlaw_csm(csm, r):
    """Sample the CSM on a radius grid; returns dict of arrays."""
    r = np.asarray(r, dtype=float)
    n = len(r)
    out = {
        "r": r,
        "rho": _csm_rho(csm, r),
        "v": _csm_v(csm, r),
        "T": _csm_T(csm, r),
    }
    if csm["composition"] == "solar":
        sc = SOLAR_COMPOSITION
        for k in ("fH", "fHe", "fO", "fSi", "fCa", "fFe"):
            out[k] = np.full(n, sc[k])
    else:
        out["fH"] = np.full(n, csm["xh"])
        out["fHe"] = np.full(n, csm["yhe"])
        for k in ("fO", "fSi", "fCa", "fFe"):
            out[k] = np.zeros(n)
    return out


def print_csm_diagnostics(csm, r_lo_cm, r_hi_cm):
    """Print CSM diagnostics."""
    rr = np.logspace(np.log10(r_lo_cm), np.log10(r_hi_cm), 2000)
    rho = _csm_rho(csm, rr)
    v = _csm_v(csm, rr)
    r_lo = np.concatenate([[r_lo_cm], rr[:-1]])
    dm = (4.0 / 3.0) * np.pi * (rr**3 - r_lo**3) * rho
    dm = np.maximum(dm, 0.0)
    m_csm = np.sum(dm)
    e_csm = 0.5 * np.sum(dm * v**2)
    tau_es = (
        KAPPA_ES * np.trapezoid(rho, rr)
        if hasattr(np, "trapezoid")
        else KAPPA_ES * np.trapz(rho, rr)
    )

    rho_15 = _csm_rho(csm, R_REF_1E15)
    n_est = rho_15 / (1.4 * M_P)  # n ~ rho/(mu m_p), mu~1.4
    w_15 = 4.0 * np.pi * R_REF_1E15**2 * rho_15
    t_yr = csm["t_cs_s"] / YR_TO_S
    mdot_equiv = (
        (
            4.0
            * np.pi
            * r_lo_cm**2
            * float(_csm_rho(csm, r_lo_cm))
            * float(_csm_v(csm, r_lo_cm))
        )
        * YR_TO_S
        / MSUN
    )
    mrate = m_csm / MSUN / max(t_yr, 1e-30)  # M_csm / t_cs, Msun/yr

    if csm.get("law") == "beta":
        print(
            f"  CSM diagnostics (beta wind: v_inf={csm['vinf_cms']/1e5:.0f} km/s, "
            f"beta={csm['beta']:g}, R*={csm['rstar_cm']:.2e} cm, "
            f"Mdot={csm['mdot_gs']*YR_TO_S/MSUN:.3e} Msun/yr):"
        )
    else:
        print(f"  CSM diagnostics (rho ~ r^{csm['slope']:.2f}):")
    print(
        f"    rho(1e15) = {rho_15:.3e} g/cc  ->  n ~ {n_est:.2e} cm^-3"
    )
    print(f"    w(1e15) = 4pi r^2 rho = {w_15:.3e} g/cm")
    print(
        f"    M_csm = {m_csm/MSUN:.3f} Msun, "
        f"E_kin = {e_csm:.2e} erg"
    )
    print(
        f"    t_cs = {csm['t_cs_s']/86400.0:.1f} d = {t_yr:.2f} yr; "
        f"u = {v[0]/1e5:.0f} .. {v[-1]/1e5:.0f} km/s"
    )
    print(
        f"    tau_es(CSM) ~ {tau_es:.2f}; "
        f"M_csm/t_cs = {mrate:.3f} Msun/yr"
    )
    t_prof = _csm_T(csm, rr)
    print(
        f"    T_csm = {t_prof[0]:.0f} .. {t_prof[-1]:.0f} K "
        f"(T15={csm['T']:.0f} K, T~rho^{csm['T_gamma']}, "
        f"floor={csm['T_min']:.0f} K)"
    )
    print(f"    equivalent steady-wind Mdot at r_inner = {mdot_equiv:.3e} " f"Msun/yr")


# ----------------------------------------------------------------------
# Combiners
# ----------------------------------------------------------------------


def combine_shell_csm_ejecta(
    shell,
    csm,
    ejecta,
    nzon_total=500,
    sw_transition_dex=0.5,
    sc_v_sharp_zones=0,
    cd_transition_zones=5,
    ce_v_ramp_zones=0,
    v_blend_zones=0,
    ejecta_T_floor=1000.0,
    ejecta_T_const=0.0,
    min_dlogr=0.005,
    max_dlogr=None,
    density_matched=True,
):
    """Combine shell + power-law CSM + ejecta.

    Layout (inner -> outer):
      shell | shell-CSM buffer | CSM | CD ramp | ejecta

    Same zone-allocation and grid-spacing logic as v1's
    combine_shell_wind_ejecta_density_matched, but the middle component
    follows the CSM laws (rho ~ r^s, homologous u = r/t_cs).

    density_matched=True: rho15 was rescaled so rho_csm(r_ej_inner) =
    rho_ejecta(r_ej_inner); the CSM law then continues through the thin
    CD ramp and the density is continuous by construction.
    density_matched=False (explicit rho15/w): the density jump at the
    interface is spread log-linearly over enough zones to keep
    |dlog rho/dlog r| < 16 (the CD ramp grows into the CSM tail if the
    5-zone default is too thin).

    ce_v_ramp_zones > 0: the last N zones before the ejecta form an
    acceleration zone - the velocity ramps (smoothstep) from the local
    CSM/wind value up to the ejecta inner-edge velocity, so the
    interface stays velocity-continuous when v_inf < v_ejecta.
    """
    from scipy.interpolate import PchipInterpolator

    buffer_max_dlogr = (
        min(max_dlogr, 2 * min_dlogr) if max_dlogr is not None else 2 * min_dlogr
    )

    ej_order = np.argsort(ejecta["r"])
    ej_r_sorted = ejecta["r"][ej_order]
    ej_rho_sorted = ejecta["rho"][ej_order]
    ej_v_sorted = ejecta["v"][ej_order]
    ej_T_sorted = ejecta["T"][ej_order]

    r_min = shell["r_inner"]
    r_shell_outer = shell["r_outer"]
    r_ej_inner = ej_r_sorted[0] * 0.99
    r_max = ej_r_sorted[-1]

    log_r_min = np.log10(r_min)
    log_r_shell = np.log10(r_shell_outer)
    log_r_ej = np.log10(r_ej_inner)
    log_r_max = np.log10(r_max)

    # --- Density normalisation at the interface ---
    log_r_ej_src = np.log10(ej_r_sorted)
    log_rho_ej_src = np.log10(np.maximum(ej_rho_sorted, 1e-300))
    rho_ej_at_interface = 10.0 ** np.interp(
        np.log10(r_ej_inner), log_r_ej_src, log_rho_ej_src
    )

    if csm.get("match_shell"):
        # Anchor the power law at the shell edge: the SC buffer ramp is
        # then collinear with the CSM, so the whole gap is one clean
        # r^-n line.  The ejecta interface keeps a density step that
        # the CD ramp spreads over its trailing zones.
        rho_shell_edge = float(shell["rho"][-1])
        _csm_set_match(csm, r_shell_outer, rho_shell_edge)
        density_matched = False
        print(
            f"  CSM shell match: rho15 -> {csm['rho15']:.3e} g/cc "
            f"(rho_shell(r_outer)={rho_shell_edge:.3e})"
        )
    if density_matched:
        if rho_ej_at_interface < 1e-30:
            raise ValueError(
                f"Ejecta density at the CSM-ejecta interface is "
                f"{rho_ej_at_interface:.3e} g/cc; a matched CSM would be "
                f"a vacuum. Use explicit --csm-rho15/--csm-w instead."
            )
        _csm_set_match(csm, r_ej_inner, rho_ej_at_interface)
        if csm.get("law") == "beta":
            print(
                f"  CSM density match: Mdot -> "
                f"{csm['mdot_gs']*YR_TO_S/MSUN:.3e} Msun/yr "
                f"(rho_ejecta(r_interface)={rho_ej_at_interface:.3e})"
            )
        else:
            print(
                f"  CSM density match: rho15 -> {csm['rho15']:.3e} g/cc "
                f"(rho_ejecta(r_interface)={rho_ej_at_interface:.3e})"
            )

    rho_csm_at_interface = _csm_rho(csm, r_ej_inner)
    cd_contrast = abs(
        np.log10(max(rho_ej_at_interface, 1e-300) / max(rho_csm_at_interface, 1e-300))
    )
    print(
        f"  CSM-ejecta interface: rho_csm={rho_csm_at_interface:.3e} vs "
        f"rho_ej={rho_ej_at_interface:.3e} g/cc "
        f"(contrast={cd_contrast:.2f} dex)"
    )

    # --- Shell-CSM buffer (same construction as v1 SW buffer) ---
    total_gap_dex = log_r_ej - log_r_shell
    max_buffer_dex = total_gap_dex / 3.0
    sw_transition_dex = min(sw_transition_dex, max_buffer_dex)

    r_sc_buffer_end = r_shell_outer * (10**sw_transition_dex)
    rho_sc_start = shell["rho"][-1]
    rho_sc_end = _csm_rho(csm, r_sc_buffer_end)
    sc_contrast = np.log10(max(rho_sc_start, 1e-300) / max(rho_sc_end, 1e-300))

    nzon_sc_buffer = max(20, int(sc_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sc_buffer = min(nzon_sc_buffer, 60)
    nzon_sc_buffer = _cap_buffer_zones(
        nzon_sc_buffer,
        sw_transition_dex,
        min_dlogr=min_dlogr,
        max_dlogr=buffer_max_dlogr,
        label="SC buffer",
    )
    print(
        f"  Shell-CSM buffer: {nzon_sc_buffer} zones, "
        f"r=[{r_shell_outer:.3e}, {r_sc_buffer_end:.3e}] cm"
    )
    print(
        f"    rho: {rho_sc_start:.3e} -> {rho_sc_end:.3e} g/cc "
        f"(contrast={sc_contrast:.1f} dex)"
    )

    # --- Zone allocation (same scheme as v1 density-matched combiner) ---
    nzon_shell = max(80, int(nzon_total * 0.20))
    nzon_ejecta = max(100, int(nzon_total * 0.55))
    nzon_csm = (
        nzon_total - nzon_shell - nzon_ejecta - nzon_sc_buffer - cd_transition_zones
    )
    nzon_csm = max(40, nzon_csm)

    csm_span_dex = log_r_ej - np.log10(r_sc_buffer_end)
    ej_span_dex = log_r_max - log_r_ej
    nzon_csm_capped = _cap_buffer_zones(
        nzon_csm + cd_transition_zones,
        csm_span_dex,
        min_dlogr=min_dlogr,
        max_dlogr=buffer_max_dlogr,
        label="CSM+CD",
    )
    nzon_csm = max(nzon_csm_capped - cd_transition_zones, 10)
    nzon_ejecta = _cap_buffer_zones(
        nzon_ejecta,
        ej_span_dex,
        min_dlogr=min_dlogr,
        max_dlogr=buffer_max_dlogr,
        label="Ejecta",
    )

    n_freed = (
        nzon_total
        - nzon_shell
        - nzon_sc_buffer
        - nzon_csm
        - cd_transition_zones
        - nzon_ejecta
    )
    if n_freed < 0:
        n_excess = -n_freed
        n_from_ej = min(n_excess, max(0, nzon_ejecta - 100))
        nzon_ejecta -= n_from_ej
        n_excess -= n_from_ej
        n_from_csm = min(n_excess, max(0, nzon_csm - 10))
        nzon_csm -= n_from_csm
        n_excess -= n_from_csm
        n_from_shell = min(n_excess, max(0, nzon_shell - 80))
        nzon_shell -= n_from_shell
        n_excess -= n_from_shell
        if n_from_ej or n_from_csm or n_from_shell:
            print(
                f"  Buffer max_dlogr enforcement grew buffers by "
                f"{-n_freed} zones; reduced ejecta by {n_from_ej}, "
                f"CSM by {n_from_csm}, shell by {n_from_shell}"
            )
        if n_excess > 0:
            print(
                f"  WARNING: {n_excess} zones could not be removed - "
                f"total will be {nzon_total + n_excess}"
            )
        n_freed = 0
    if n_freed > 0:
        shell_span_dex = log_r_shell - log_r_min
        max_shell = max(int(shell_span_dex / min_dlogr), nzon_shell)
        n_to_shell = min(n_freed, max_shell - nzon_shell)
        if n_to_shell > 0:
            nzon_shell += n_to_shell
            n_freed -= n_to_shell
            print(f"  Redistributed {n_to_shell} freed zones to shell")
        if n_freed > 0:
            print(
                f"  WARNING: {n_freed} zones could not be placed - "
                f"total will be {nzon_total - n_freed}"
            )

    n_total = nzon_shell + nzon_sc_buffer + nzon_csm + cd_transition_zones + nzon_ejecta
    print(
        f"  Zone allocation: {nzon_shell} shell + {nzon_sc_buffer} SC buf "
        f"+ {nzon_csm} CSM + {cd_transition_zones} CD ramp "
        f"+ {nzon_ejecta} ejecta = {n_total}"
    )

    # --- Shell grid (split core/envelope as in v1) ---
    r0_shell = shell.get("r0", None)
    if r0_shell is not None and r0_shell > shell["r_inner"]:
        shell_r = _build_optimized_shell_grid(
            r_min,
            r0_shell,
            r_shell_outer,
            nzon_shell,
            core_fraction=0.25,
            min_dlogr=min_dlogr,
        )
    else:
        shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    log_r_src = np.log10(shell["r"])
    log_rho_src = np.log10(np.maximum(shell["rho"], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell["v"])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell["T"])

    # Shell mass conservation (same as v1)
    if "rho_0" in shell and "r0" in shell:
        n_rho_val = shell.get("n_rho", 8.0)
        frac_v0 = shell["r0"] / shell["r_outer"]
        factor = _shell_mass_factor(n_rho_val, frac_v0)
        M_shell_target = (4 / 3) * np.pi * shell["r0"] ** 3 * shell["rho_0"] * factor
    else:
        r_lo_src = np.concatenate([[shell["r_inner"]], shell["r"][:-1]])
        M_shell_target = (
            (4 / 3) * np.pi * (shell["r"] ** 3 - r_lo_src**3) * shell["rho"]
        ).sum()
    r_lo_shell = np.concatenate([[shell_r[0] * 0.5], shell_r[:-1]])
    M_shell_interp = ((4 / 3) * np.pi * (shell_r**3 - r_lo_shell**3) * shell_rho).sum()
    if M_shell_interp > 0:
        shell_rho *= M_shell_target / M_shell_interp
        print(
            f"  Shell mass conservation: target={M_shell_target/MSUN:.5f} "
            f"Msun, renorm={M_shell_target/M_shell_interp:.4f}"
        )

    # --- Shell-CSM buffer: slope-matched density ramp ---
    sc_buf_r = np.logspace(log_r_shell, np.log10(r_sc_buffer_end), nzon_sc_buffer + 1)[
        1:
    ]
    sc_buf_rho = _hermite_log_ramp(
        sc_buf_r,
        r_shell_outer,
        r_sc_buffer_end,
        rho_sc_start,
        rho_sc_end,
        _shell_edge_slope(shell),
        float(_csm_log_slope(csm, r_sc_buffer_end)),
    )
    if sc_v_sharp_zones > 0:
        # Sharp velocity drop: the buffer gas IS CSM, so it follows the
        # homologous u=r/t_cs law everywhere except the first few zones,
        # which resolve the shell->CSM jump like a shock front.  Only the
        # velocity is sharpened; the density still needs the full ramp
        # to keep |dlog rho/dlog r| under control.
        k = min(sc_v_sharp_zones, nzon_sc_buffer - 2)
        sc_buf_v = _csm_v(csm, sc_buf_r)
        sc_buf_v[:k] = np.linspace(shell["v"][-1], sc_buf_v[k], k)
    else:
        sc_buf_v = np.linspace(
            shell["v"][-1], float(_csm_v(csm, r_sc_buffer_end)), nzon_sc_buffer
        )
    # Temperature blends use smoothstep (zero slope at both ends), like
    # the composition blending, to avoid kinks in T at the interfaces.
    sc_tw = np.linspace(0.0, 1.0, nzon_sc_buffer)
    sc_tw = sc_tw * sc_tw * (3.0 - 2.0 * sc_tw)
    sc_buf_T = shell["T"][-1] * (1.0 - sc_tw) + float(
        _csm_T(csm, r_sc_buffer_end)
    ) * sc_tw

    # --- CSM grid: pure CSM + CD ramp share one log-spaced span ---
    csm_r = np.logspace(
        np.log10(r_sc_buffer_end), log_r_ej, nzon_csm + cd_transition_zones + 1
    )[1:]
    csm_rho = _csm_rho(csm, csm_r)
    csm_v = _csm_v(csm, csm_r)
    csm_T = _csm_T(csm, csm_r)

    # CD temperature ramp: smoothstep from the CSM temperature at the CD
    # start to the ejecta inner-edge T across the thin CD segment
    cd_tw = np.linspace(0.0, 1.0, cd_transition_zones)
    cd_tw = cd_tw * cd_tw * (3.0 - 2.0 * cd_tw)
    ej_T_end = (
        ejecta_T_const
        if ejecta_T_const > 0
        else max(ej_T_sorted[0], ejecta_T_floor)
    )
    csm_T[nzon_csm:] = csm_T[nzon_csm] * (1.0 - cd_tw) + ej_T_end * cd_tw

    # Acceleration zone: ramp the last ce_v_ramp_zones of the CSM+CD
    # velocity up to the ejecta inner-edge velocity (smoothstep, zero
    # slope at both ends). Density is left on the CSM law - the ramp is
    # a resolved transition region, like the shell-side buffer.
    if ce_v_ramp_zones > 0:
        ramp_n = min(ce_v_ramp_zones, len(csm_v) - 1)
        i0 = len(csm_v) - ramp_n
        tw = np.linspace(0.0, 1.0, ramp_n)
        tw = tw * tw * (3.0 - 2.0 * tw)
        v_start = csm_v[i0]
        csm_v[i0:] = v_start * (1.0 - tw) + ej_v_sorted[0] * tw
        print(
            f"  CE velocity ramp: {ramp_n} zones, "
            f"{v_start/1e5:.1f} -> {ej_v_sorted[0]/1e5:.1f} km/s"
        )

    # --- CD density handling ---
    if not density_matched:
        # Spread the density jump log-linearly over enough trailing zones
        # of the CSM+CD grid to keep |dlog rho/dlog r| < 16.
        grid_dlogr = np.diff(np.log10(csm_r)).mean()
        need = int(np.ceil(cd_contrast / 16.0 / max(grid_dlogr, 1e-6)))
        ramp_n = max(cd_transition_zones, min(need, len(csm_r) - 1))
        i0 = len(csm_r) - ramp_n
        rho_start = csm_rho[i0]
        csm_rho[i0:] = 10.0 ** np.linspace(
            np.log10(max(rho_start, 1e-300)),
            np.log10(max(rho_ej_at_interface, 1e-300)),
            ramp_n,
        )
        print(
            f"  CD density ramp (unmatched): {ramp_n} zones, "
            f"{rho_start:.3e} -> {rho_ej_at_interface:.3e} g/cc"
        )
    # matched case: csm_rho is already continuous through the CD zones.

    csm_r_pure = csm_r[:nzon_csm]
    cd_r = csm_r[nzon_csm:]
    csm_rho_pure = csm_rho[:nzon_csm]
    cd_rho = csm_rho[nzon_csm:]
    csm_v_pure = csm_v[:nzon_csm]
    cd_v = csm_v[nzon_csm:]
    csm_T_pure = csm_T[:nzon_csm]
    cd_T = csm_T[nzon_csm:]

    # --- Ejecta grid (PCHIP, same as v1) ---
    ej_r_new = np.logspace(log_r_ej, log_r_max, nzon_ejecta + 1)[1:]

    def log_interp(y, r_target):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = PchipInterpolator(log_r_ej_src, log_y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return 10.0 ** ip(lr)

    def lin_interp(y, r_target):
        ip = PchipInterpolator(log_r_ej_src, y, extrapolate=False)
        lr = np.clip(np.log10(r_target), log_r_ej_src[0], log_r_ej_src[-1])
        return ip(lr)

    ej_rho = log_interp(ej_rho_sorted, ej_r_new)
    ej_v = lin_interp(ej_v_sorted, ej_r_new)
    if ejecta_T_const > 0:
        # Flat ejecta temperature requested: a smooth, isothermal initial
        # condition - the SN_A source T profile is ignored.
        ej_T = np.full(nzon_ejecta, ejecta_T_const)
    else:
        ej_T = np.maximum(log_interp(ej_T_sorted, ej_r_new), ejecta_T_floor)

    # Optional velocity blend at the interface (off by default; CMFGEN
    # handles non-monotonic velocity and STELLA only needs smooth rho).
    # Use the actual last CD-zone velocity so the ce_v_ramp_zones
    # acceleration is reflected in the diagnostic.
    u_cd = float(csm_v[-1])
    if v_blend_zones > 1:
        v_blend_zones = min(v_blend_zones, nzon_ejecta // 3)
        blend_w = np.linspace(0.0, 1.0, v_blend_zones)
        full_w = np.concatenate([blend_w, np.ones(nzon_ejecta - v_blend_zones)])
        ej_v = u_cd * (1.0 - full_w) + ej_v * full_w
        print(
            f"  Velocity blend: {v_blend_zones} zones, "
            f"u_csm={u_cd/1e5:.1f} km/s -> v_ej(r)"
        )
    else:
        dv_kms = (ej_v_sorted[0] - u_cd) / 1e5
        print(
            f"  Velocity blend: SKIPPED "
            f"(u_csm={u_cd/1e5:.1f} vs v_ej={ej_v_sorted[0]/1e5:.1f} km/s, "
            f"jump={dv_kms:+.1f} km/s at interface)"
        )

    # Light ejecta T and rho smoothing (same as v1)
    if len(ej_T) > 5:
        kernel = np.ones(5) / 5.0
        ej_T_smooth = np.convolve(ej_T, kernel, mode="same")
        ej_T_smooth[:2] = ej_T[:2]
        ej_T_smooth[-2:] = ej_T[-2:]
        ej_T = np.maximum(ej_T_smooth, ejecta_T_floor)

    for _ in range(3):
        log_ej_rho = np.log10(np.maximum(ej_rho, 1e-300))
        log_ej_r = np.log10(ej_r_new)
        grads = np.abs(np.diff(log_ej_rho) / np.diff(log_ej_r))
        smoothed = ej_rho.copy()
        for i in range(1, len(ej_rho) - 1):
            if grads[i - 1] > 16 or grads[i] > 16:
                smoothed[i] = np.sqrt(max(ej_rho[i - 1] * ej_rho[i + 1], 1e-300))
        ej_rho = smoothed

    # Pin the first ejecta zone to the last CD-zone density for exact
    # continuity (a 1-2% jump still triggers AQ oscillations in STELLA).
    if len(ej_rho) > 0 and len(cd_rho) > 0:
        ej_rho[0] = cd_rho[-1]
        for i in range(1, min(5, len(ej_rho))):
            ej_rho[i] = np.sqrt(max(ej_rho[i - 1] * ej_rho[i], 1e-300))

    # --- Combine ---
    r_combined = np.concatenate([shell_r, sc_buf_r, csm_r_pure, cd_r, ej_r_new])
    rho_combined = np.concatenate([shell_rho, sc_buf_rho, csm_rho_pure, cd_rho, ej_rho])
    v_combined = np.concatenate([shell_v, sc_buf_v, csm_v_pure, cd_v, ej_v])
    T_combined = np.concatenate([shell_T, sc_buf_T, csm_T_pure, cd_T, ej_T])

    # --- Composition ---
    all_elem_keys = ["fH", "fHe", "fO", "fSi", "fCa", "fFe"]
    ej_comp = {}
    for key in all_elem_keys:
        if key in ejecta:
            ej_comp[key] = lin_interp(ejecta[key][ej_order], ej_r_new)
        else:
            ej_comp[key] = np.zeros(nzon_ejecta)

    ej_comp_inner = {}
    for key in all_elem_keys:
        ej_comp_inner[key] = float(ejecta[key][ej_order][0]) if key in ejecta else 0.0

    shell_comp_vals = {k: float(shell.get(k, [0.0])[0]) for k in all_elem_keys}
    csm_sample = sample_powerlaw_csm(csm, np.array([r_sc_buffer_end]))
    csm_comp_vals = {k: float(csm_sample[k][0]) for k in all_elem_keys}

    sc_u = np.linspace(0.0, 1.0, nzon_sc_buffer)
    sc_w = sc_u * sc_u * (3.0 - 2.0 * sc_u)
    cd_w = np.linspace(0.0, 1.0, cd_transition_zones)

    comp_combined = {}
    for key in all_elem_keys:
        sv = shell_comp_vals[key]
        cv = csm_comp_vals[key]
        ev = ej_comp_inner[key]
        sc_buf_comp = sv * (1.0 - sc_w) + cv * sc_w
        cd_comp = cv * (1.0 - cd_w) + ev * cd_w
        comp_combined[key] = np.concatenate(
            [
                np.full(nzon_shell, sv),
                sc_buf_comp,
                np.full(nzon_csm, cv),
                cd_comp,
                ej_comp[key],
            ]
        )

    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    # --- Masses and gradient checks ---
    r_lo = np.concatenate([[r_combined[0] * 0.5], r_combined[:-1]])
    dm = (4 / 3) * np.pi * (r_combined**3 - r_lo**3) * rho_combined
    dm = np.maximum(dm, 0.0)
    M_total = np.sum(dm)

    idx_shell = nzon_shell
    idx_sc = idx_shell + nzon_sc_buffer
    idx_csm = idx_sc + nzon_csm
    idx_cd = idx_csm + cd_transition_zones

    print(
        f"  Mass: shell={np.sum(dm[:idx_shell])/MSUN:.3f} Msun, "
        f"SC buf={np.sum(dm[idx_shell:idx_sc])/MSUN:.4f} Msun, "
        f"CSM={np.sum(dm[idx_sc:idx_csm])/MSUN:.4f} Msun, "
        f"CD={np.sum(dm[idx_csm:idx_cd])/MSUN:.6f} Msun, "
        f"ejecta={np.sum(dm[idx_cd:])/MSUN:.3f} Msun, "
        f"total={M_total/MSUN:.3f} Msun"
    )

    log_r = np.log10(r_combined)
    log_rho_arr = np.log10(np.maximum(rho_combined, 1e-300))
    if len(log_r) > 1:
        gradients = np.abs(np.diff(log_rho_arr) / np.diff(log_r))
        grad_interior = gradients[:-5]
        max_grad_interior = np.max(grad_interior) if len(grad_interior) else 0
        print(
            f"  Max |dlog rho/dlog r|: {np.max(gradients):.1f} "
            f"(interior: {max_grad_interior:.1f})"
        )
        if max_grad_interior > 16:
            print(
                f"  WARNING: Interior density gradient "
                f"{max_grad_interior:.1f} > 16 - may cause solver stall"
            )

    return {
        "r": r_combined,
        "rho": rho_combined,
        "v": v_combined,
        "T": T_combined,
        "idx_ej": idx_cd,
        "idx_csm": idx_sc,
        **comp_combined,
    }


def combine_shell_csm(
    shell,
    csm,
    nzon_total=500,
    sw_transition_dex=0.5,
    sc_v_sharp_zones=0,
    min_dlogr=0.005,
    max_dlogr=None,
):
    """Combine shell + power-law CSM only (--no-ejecta mode).

    Layout: shell | shell-CSM buffer | CSM out to csm['r_outer'].
    Mirrors v1's combine_shell_wind with the CSM replacing the wind.
    """
    buffer_max_dlogr = (
        min(max_dlogr, 2 * min_dlogr) if max_dlogr is not None else 2 * min_dlogr
    )

    r_min = shell["r_inner"]
    r_shell_outer = shell["r_outer"]
    r_max = csm["r_outer"]

    if csm.get("match_shell"):
        rho_shell_edge = float(shell["rho"][-1])
        _csm_set_match(csm, r_shell_outer, rho_shell_edge)
        print(
            f"  CSM shell match: rho15 -> {csm['rho15']:.3e} g/cc "
            f"(rho_shell(r_outer)={rho_shell_edge:.3e})"
        )

    log_r_min = np.log10(r_min)
    log_r_shell = np.log10(r_shell_outer)
    log_r_max = np.log10(r_max)

    total_gap_dex = log_r_max - log_r_shell
    sw_transition_dex = min(sw_transition_dex, total_gap_dex / 3.0)

    r_sc_buffer_end = r_shell_outer * (10**sw_transition_dex)
    rho_sc_start = shell["rho"][-1]
    rho_sc_end = _csm_rho(csm, r_sc_buffer_end)
    sc_contrast = np.log10(max(rho_sc_start, 1e-300) / max(rho_sc_end, 1e-300))

    nzon_sc_buffer = max(20, int(sc_contrast / 16 * sw_transition_dex / 0.02))
    nzon_sc_buffer = min(nzon_sc_buffer, 60)
    nzon_sc_buffer = _cap_buffer_zones(
        nzon_sc_buffer,
        sw_transition_dex,
        min_dlogr=min_dlogr,
        max_dlogr=buffer_max_dlogr,
        label="SC buffer",
    )
    print(
        f"  Shell-CSM buffer: {nzon_sc_buffer} zones, "
        f"r=[{r_shell_outer:.3e}, {r_sc_buffer_end:.3e}] cm"
    )
    print(
        f"    rho: {rho_sc_start:.3e} -> {rho_sc_end:.3e} g/cc "
        f"(contrast={sc_contrast:.1f} dex)"
    )

    nzon_shell = max(100, int(nzon_total * 0.40))
    nzon_csm = nzon_total - nzon_shell - nzon_sc_buffer
    nzon_csm = max(60, nzon_csm)

    csm_span_dex = log_r_max - np.log10(r_sc_buffer_end)
    nzon_csm = _cap_buffer_zones(
        nzon_csm,
        csm_span_dex,
        min_dlogr=min_dlogr,
        max_dlogr=buffer_max_dlogr,
        label="CSM",
    )

    n_freed = nzon_total - nzon_shell - nzon_sc_buffer - nzon_csm
    if n_freed < 0:
        n_excess = -n_freed
        n_from_shell = min(n_excess, max(0, nzon_shell - 100))
        nzon_shell -= n_from_shell
        n_excess -= n_from_shell
        n_from_csm = min(n_excess, max(0, nzon_csm - 60))
        nzon_csm -= n_from_csm
        if n_excess > 0:
            print(
                f"  WARNING: {n_excess} zones could not be removed - "
                f"total will be {nzon_total + n_excess}"
            )
        n_freed = 0
    if n_freed > 0:
        shell_span_dex = log_r_shell - log_r_min
        max_shell = max(int(shell_span_dex / min_dlogr), nzon_shell)
        n_to_shell = min(n_freed, max_shell - nzon_shell)
        if n_to_shell > 0:
            nzon_shell += n_to_shell
            n_freed -= n_to_shell
            print(f"  Redistributed {n_to_shell} freed zones to shell")
        if n_freed > 0:
            print(
                f"  WARNING: {n_freed} zones could not be placed - "
                f"total will be {nzon_total - n_freed}"
            )

    n_total = nzon_shell + nzon_sc_buffer + nzon_csm
    print(
        f"  Zone allocation: {nzon_shell} shell + {nzon_sc_buffer} SC buf "
        f"+ {nzon_csm} CSM = {n_total}"
    )

    r0_shell = shell.get("r0", None)
    if r0_shell is not None and r0_shell > shell["r_inner"]:
        shell_r = _build_optimized_shell_grid(
            r_min,
            r0_shell,
            r_shell_outer,
            nzon_shell,
            core_fraction=0.25,
            min_dlogr=min_dlogr,
        )
    else:
        shell_r = np.logspace(log_r_min, log_r_shell, nzon_shell)
    log_r_src = np.log10(shell["r"])
    log_rho_src = np.log10(np.maximum(shell["rho"], 1e-300))
    shell_rho = 10.0 ** np.interp(np.log10(shell_r), log_r_src, log_rho_src)
    shell_v = np.interp(np.log10(shell_r), log_r_src, shell["v"])
    shell_T = np.interp(np.log10(shell_r), log_r_src, shell["T"])

    if "rho_0" in shell and "r0" in shell:
        n_rho_val = shell.get("n_rho", 8.0)
        frac_v0 = shell["r0"] / shell["r_outer"]
        factor = _shell_mass_factor(n_rho_val, frac_v0)
        M_shell_target = (4 / 3) * np.pi * shell["r0"] ** 3 * shell["rho_0"] * factor
    else:
        r_lo_src = np.concatenate([[shell["r_inner"]], shell["r"][:-1]])
        M_shell_target = (
            (4 / 3) * np.pi * (shell["r"] ** 3 - r_lo_src**3) * shell["rho"]
        ).sum()
    r_lo_shell = np.concatenate([[shell_r[0] * 0.5], shell_r[:-1]])
    M_shell_interp = ((4 / 3) * np.pi * (shell_r**3 - r_lo_shell**3) * shell_rho).sum()
    if M_shell_interp > 0:
        shell_rho *= M_shell_target / M_shell_interp
        print(
            f"  Shell mass conservation: target={M_shell_target/MSUN:.5f} "
            f"Msun, renorm={M_shell_target/M_shell_interp:.4f}"
        )

    sc_buf_r = np.logspace(log_r_shell, np.log10(r_sc_buffer_end), nzon_sc_buffer + 1)[
        1:
    ]
    sc_buf_rho = _hermite_log_ramp(
        sc_buf_r,
        r_shell_outer,
        r_sc_buffer_end,
        rho_sc_start,
        rho_sc_end,
        _shell_edge_slope(shell),
        float(_csm_log_slope(csm, r_sc_buffer_end)),
    )
    if sc_v_sharp_zones > 0:
        k = min(sc_v_sharp_zones, nzon_sc_buffer - 2)
        sc_buf_v = _csm_v(csm, sc_buf_r)
        sc_buf_v[:k] = np.linspace(shell["v"][-1], sc_buf_v[k], k)
    else:
        sc_buf_v = np.linspace(
            shell["v"][-1], float(_csm_v(csm, r_sc_buffer_end)), nzon_sc_buffer
        )
    sc_tw = np.linspace(0.0, 1.0, nzon_sc_buffer)
    sc_tw = sc_tw * sc_tw * (3.0 - 2.0 * sc_tw)
    sc_buf_T = shell["T"][-1] * (1.0 - sc_tw) + float(
        _csm_T(csm, r_sc_buffer_end)
    ) * sc_tw

    csm_r = np.logspace(np.log10(r_sc_buffer_end), log_r_max, nzon_csm + 1)[1:]
    csm_rho = _csm_rho(csm, csm_r)
    csm_v = _csm_v(csm, csm_r)
    csm_T = _csm_T(csm, csm_r)

    r_combined = np.concatenate([shell_r, sc_buf_r, csm_r])
    rho_combined = np.concatenate([shell_rho, sc_buf_rho, csm_rho])
    v_combined = np.concatenate([shell_v, sc_buf_v, csm_v])
    T_combined = np.concatenate([shell_T, sc_buf_T, csm_T])

    all_elem_keys = ["fH", "fHe", "fO", "fSi", "fCa", "fFe"]
    shell_comp_vals = {k: float(shell.get(k, [0.0])[0]) for k in all_elem_keys}
    csm_sample = sample_powerlaw_csm(csm, np.array([r_sc_buffer_end]))
    csm_comp_vals = {k: float(csm_sample[k][0]) for k in all_elem_keys}

    sc_u = np.linspace(0.0, 1.0, nzon_sc_buffer)
    sc_w = sc_u * sc_u * (3.0 - 2.0 * sc_u)
    comp_combined = {}
    for key in all_elem_keys:
        sv = shell_comp_vals[key]
        cv = csm_comp_vals[key]
        comp_combined[key] = np.concatenate(
            [
                np.full(nzon_shell, sv),
                sv * (1.0 - sc_w) + cv * sc_w,
                np.full(nzon_csm, cv),
            ]
        )

    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    r_lo = np.concatenate([[r_combined[0] * 0.5], r_combined[:-1]])
    dm = (4 / 3) * np.pi * (r_combined**3 - r_lo**3) * rho_combined
    dm = np.maximum(dm, 0.0)
    print(
        f"  Mass: shell={np.sum(dm[:nzon_shell])/MSUN:.3f} Msun, "
        f"SC buf={np.sum(dm[nzon_shell:nzon_shell+nzon_sc_buffer])/MSUN:.4f} "
        f"Msun, CSM={np.sum(dm[nzon_shell+nzon_sc_buffer:])/MSUN:.4f} Msun, "
        f"total={np.sum(dm)/MSUN:.3f} Msun"
    )

    return {
        "r": r_combined,
        "rho": rho_combined,
        "v": v_combined,
        "T": T_combined,
        "idx_csm": nzon_shell,
        **comp_combined,
    }


def write_v2_input_profile(fname, g):
    """Write 12-column input_profile."""
    r = g["r"]
    rho = g["rho"]
    T = g["T"]

    fH = g["fH"]
    fHe = g["fHe"]
    fO = g["fO"]
    fSi = g["fSi"]
    fFe = g["fFe"]
    fCa = g.get("fCa", np.zeros_like(fH))

    mu = 1.0 / (fH / 1.0 + fHe / 4.0 + fO / 16.0 + fSi / 28.0 + fCa / 40.0 + fFe / 56.0)
    mu = np.maximum(mu, 0.5)
    p_gas = rho * K_B * T / (mu * M_H)
    p_rad = A_RAD * T**4 / 3.0
    p = p_gas + p_rad

    header = (
        "# input_profile generated by prepare_SN_B_model_v2.py\n"
        "# columns: r_cm rho_cgs p_cgs u_cms T_K fx_H fx_He fx_O "
        "fx_Si fx_Ca fx_Fe\n"
        "# Inner: H/He shell; Middle: power-law CSM "
        "(flat rho, homologous u); Outer: SN_A ejecta\n"
    )
    with open(fname, "w") as f:
        f.write(header)
        for i in range(len(r)):
            f.write(
                f"{r[i]:.8e} {rho[i]:.8e} {p[i]:.8e} {g['v'][i]:.8e} "
                f"{T[i]:.8e} "
                f"{fH[i]:.8e} {fHe[i]:.8e} {fO[i]:.8e} "
                f"{fSi[i]:.8e} {fCa[i]:.8e} {fFe[i]:.8e} {0.0:.8e}\n"
            )
    print(f"Written {fname}  ({len(r)} zones)")


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "source_dir",
        nargs="?",
        default=None,
        help="Source SN model dir (SN_A); unused with --no-ejecta",
    )
    parser.add_argument("output_dir", help="SN_B model output directory")
    parser.add_argument("--model-name", default=None)

    # Shell (same as v1)
    parser.add_argument("--shell-mass", type=float, default=1.0)
    parser.add_argument("--shell-velocity", type=float, default=5000)
    parser.add_argument("--shell-xh", type=float, default=0.7)
    parser.add_argument("--shell-yhe", type=float, default=0.3)
    parser.add_argument(
        "--shell-composition", type=str, default="hhe", choices=["hhe", "solar"]
    )
    parser.add_argument("--n-rho", type=float, default=8.0)
    parser.add_argument("--v0", type=float, default=3000.0)
    parser.add_argument("--n-T", type=float, default=0.4)
    parser.add_argument("--T-interface", type=float, default=5000.0)
    parser.add_argument("--T-max", type=float, default=65000.0)
    parser.add_argument("--t-start", type=float, default=10.0)
    parser.add_argument("--tcura-days", type=float, default=None)
    parser.add_argument("--ejecta-T-floor", type=float, default=1000.0)
    parser.add_argument(
        "--ejecta-T-const",
        type=float,
        default=0.0,
        help="Force a constant ejecta temperature in K, ignoring the "
        "source T profile (0 = off). For smooth initial conditions.",
    )
    parser.add_argument("--ejecta-rho-min", type=float, default=0.0)
    parser.add_argument("--ejecta-r-max-cm", type=float, default=0.0)
    parser.add_argument("--ejecta-v-max-kms", type=float, default=0.0)
    parser.add_argument(
        "--ejecta-r-inner-cm",
        type=float,
        default=0.0,
        help="Drop ejecta zones inside this radius (cm), "
        "moving the CSM-ejecta interface OUTWARD. "
        "Discards tenuous inner ejecta material and "
        "enlarges the CSM region",
    )

    # Wind args: accepted for CLI compatibility with run_model_grid.py;
    # only --wind-T and --wind-r-outer-cm are reused (CSM defaults).
    parser.add_argument(
        "--wind-mdot",
        type=float,
        default=0.0,
        help="Wind mass-loss rate in Msun/yr; used only with "
        "--beta-wind (0 = density-matched to the ejecta)",
    )
    parser.add_argument(
        "--wind-velocity",
        type=float,
        default=0.0,
        help="Wind terminal velocity v_inf in km/s; used only with "
        "--beta-wind (0 = ejecta inner-edge velocity)",
    )
    parser.add_argument(
        "--wind-T", type=float, default=5000, help="Fallback for --csm-T"
    )
    parser.add_argument(
        "--wind-r-outer-cm",
        type=float,
        default=0.0,
        help="Fallback for --csm-r-outer-cm (--no-ejecta)",
    )
    parser.add_argument(
        "--wind-duration-days", type=float, default=0.0, help="IGNORED in v2"
    )
    parser.add_argument(
        "--wind-mdot-bg", type=float, default=1e-7, help="IGNORED in v2"
    )
    parser.add_argument(
        "--we-transition-dex",
        type=float,
        default=0.5,
        help="IGNORED in v2 (thin CD ramp instead)",
    )
    parser.add_argument(
        "--we-transition-k", type=float, default=3.0, help="IGNORED in v2"
    )
    parser.add_argument(
        "--sw-transition-k", type=float, default=3.0, help="IGNORED in v2"
    )
    parser.add_argument(
        "--density-matched-wind",
        action="store_true",
        help="v2: same as default CSM density matching",
    )
    parser.add_argument(
        "--se-transition-dex",
        type=float,
        default=0.5,
        help="Shell-ejecta buffer width (--no-wind mode)",
    )

    # Run control (same as v1)
    parser.add_argument("--target-days", type=float, default=45.0)
    parser.add_argument("--snapshot-day", type=float, default=30.0)
    parser.add_argument("--nzon", type=int, default=500)
    parser.add_argument("--source-age", type=float, default=400.0)
    parser.add_argument("--explosion-energy", type=float, default=0.0)
    parser.add_argument("--nstmax", type=int, default=2000000)
    parser.add_argument("--aq", type=float, default=10.0)
    parser.add_argument("--fit-tau", type=float, default=5.0)
    parser.add_argument("--meth", type=int, default=3)
    parser.add_argument("--temp-floor", type=float, default=0.03)
    parser.add_argument("--floor-rad", type=float, default=0.01)
    parser.add_argument("--hmin", type=float, default=1e-14)
    parser.add_argument("--hmax", type=float, default=3600.0)
    parser.add_argument("--wacc-r", type=float, default=1.0)
    parser.add_argument("--wacc-v", type=float, default=1.0)
    parser.add_argument("--wacc-t", type=float, default=0.5)
    parser.add_argument("--wacc-rad", type=float, default=0.5)
    parser.add_argument("--tau-tol", type=float, default=1.3)
    parser.add_argument("--iout", type=int, default=10)
    parser.add_argument("--disable-scattering", action="store_true")
    parser.add_argument("--min-dlogr", type=float, default=0.005)
    parser.add_argument("--max-dlogr", type=float, default=0.15)
    parser.add_argument("--strict-grid", action="store_true")
    parser.add_argument("--sw-transition-dex", type=float, default=0.5)
    parser.add_argument(
        "--sc-v-sharp-zones",
        type=int,
        default=0,
        help="Resolve the shell->CSM velocity drop over this many "
        "buffer zones (0 = gradual ramp across the whole buffer). "
        "The rest of the buffer follows the homologous u=r/t_cs law",
    )
    parser.add_argument(
        "--no-wind",
        action="store_true",
        help="Skip CSM entirely: immediate shell-ejecta " "contact (v1 behaviour)",
    )
    parser.add_argument(
        "--immediate-contact", action="store_true", help="Alias for --no-wind"
    )
    parser.add_argument(
        "--no-ejecta",
        action="store_true",
        help="Shell + power-law CSM only",
    )
    parser.add_argument("--cd-transition-zones", type=int, default=5)
    parser.add_argument(
        "--ce-v-ramp-zones",
        type=int,
        default=0,
        help="Zones for the CSM-ejecta acceleration ramp: the last N "
        "zones before the ejecta ramp the velocity smoothly up to the "
        "ejecta inner-edge velocity (0 = off)",
    )
    parser.add_argument("--v-blend-zones", type=int, default=0)
    parser.add_argument("--force-homologous-velocity", action="store_true")
    parser.add_argument("--no-smooth-ejecta", action="store_true")
    parser.add_argument("--smooth-window", type=int, default=15)
    parser.add_argument("--smooth-polyorder", type=int, default=3)
    parser.add_argument("--comp-smooth-window", type=int, default=7)

    # Power-law CSM (new in v2)
    parser.add_argument(
        "--csm-powerlaw-n",
        type=float,
        default=2.0,
        help="CSM density power-law index n in rho ~ r^-n "
        "(n=1,2,3 controls how steeply the density falls "
        "between the shell and the ejecta)",
    )
    parser.add_argument(
        "--beta-wind",
        action="store_true",
        help="Beta-law steady wind instead of the power-law CSM: "
        "v(r) = v_inf (1 - R*/r)^beta, rho = Mdot/(4 pi r^2 v). "
        "v_inf comes from --wind-velocity (0 = ejecta inner-edge "
        "velocity), Mdot from --wind-mdot (0 = density-matched "
        "to the ejecta), R* from --wind-rstar-cm (0 = shell "
        "outer edge), beta from --wind-beta",
    )
    parser.add_argument(
        "--wind-beta",
        type=float,
        default=1.0,
        help="Beta velocity-law index (0.5-1.0 for hot massive "
        "stars); only used with --beta-wind",
    )
    parser.add_argument(
        "--wind-clumping",
        type=float,
        default=1.0,
        help="Wind clumping factor D: rho -> D*rho (rho^2-sensitive "
        "diagnostics); only used with --beta-wind",
    )
    parser.add_argument(
        "--wind-rstar-cm",
        type=float,
        default=0.0,
        help="Wind launch radius R* for the beta law (cm); "
        "0 = shell outer edge. Must be inside the CSM inner edge",
    )
    parser.add_argument(
        "--csm-match-shell",
        action="store_true",
        help="Anchor the CSM power law at the shell edge density "
        "instead of the ejecta: the whole gap follows one clean "
        "rho ~ r^-n line; the ejecta interface keeps a density "
        "step spread over the CD ramp. Overrides density match "
        "and --csm-rho15/--csm-w",
    )
    parser.add_argument(
        "--csm-age-days",
        type=float,
        default=0.0,
        help="CSM kinematic age t_cs in days; 0 = "
        "source_age + t_start (velocity-continuous "
        "with the ejecta)",
    )
    parser.add_argument(
        "--csm-rho15",
        type=float,
        default=0.0,
        help="CSM density at 1e15 cm in g/cc "
        "g/cc; overrides density match",
    )
    parser.add_argument(
        "--csm-w",
        type=float,
        default=0.0,
        help="w = 4 pi r^2 rho at 1e15 cm in g/cm "
        "g/cm; alternative normalisation",
    )
    parser.add_argument(
        "--csm-density-match",
        dest="csm_density_match",
        action="store_true",
        default=True,
    )
    parser.add_argument(
        "--no-csm-density-match",
        dest="csm_density_match",
        action="store_false",
        help="Do not match CSM density to the ejecta "
        "inner edge; use --csm-rho15/--csm-w",
    )
    parser.add_argument(
        "--csm-T",
        type=float,
        default=None,
        help="CSM temperature at 1e15 cm in K (default: --wind-T); "
        "with --csm-T-gamma>0 this is the power-law normalization "
        "T ~ rho^gamma, isothermal if gamma=0",
    )
    parser.add_argument(
        "--csm-T-gamma",
        type=float,
        default=0.3,
        help="Exponent for the CSM temperature law T ~ rho^gamma "
        "(0 = isothermal). Cools "
        "the tenuous cutoff tail toward --csm-T-min",
    )
    parser.add_argument(
        "--csm-T-min",
        type=float,
        default=None,
        help="CSM temperature floor in K for the T ~ rho^gamma law "
        "(default: --ejecta-T-floor, so the tail joins the ejecta)",
    )
    parser.add_argument(
        "--csm-composition",
        type=str,
        default=None,
        choices=["hhe", "solar"],
        help="CSM composition (default: --shell-composition; "
        "solar gives solar abundances)",
    )
    parser.add_argument(
        "--csm-ua0",
        type=float,
        default=0.0,
        help="Radiative acceleration u0 in km/s " "(0 = off)",
    )
    parser.add_argument(
        "--csm-r-outer-cm",
        type=float,
        default=0.0,
        help="CSM outer radius (--no-ejecta mode only; "
        "cm). With ejecta, the CSM always "
        "reaches the ejecta inner edge - use "
        "--csm-cutoff-r-cm for a bounded envelope",
    )
    parser.add_argument(
        "--csm-cutoff-r-cm",
        type=float,
        default=0.0,
        help="Sharp outer cutoff radius Rk in cm "
        "(0 = no cutoff). With ejecta, "
        "the CSM density drops steeply at Rk and the "
        "CD ramp bridges the tail to the ejecta",
    )
    parser.add_argument(
        "--csm-cutoff-n",
        type=float,
        default=12.0,
        help="Cutoff steepness n",
    )

    args = parser.parse_args()

    if args.immediate_contact:
        args.no_wind = True
    if args.no_ejecta and args.no_wind:
        parser.error(
            "--no-ejecta with --no-wind leaves a shell-only "
            "model, which is not supported"
        )

    source_dir = Path(args.source_dir) if args.source_dir else None
    output_dir = Path(args.output_dir)
    model_name = args.model_name if args.model_name else output_dir.name

    stella_dir = output_dir / "stella"
    modmake_dir = stella_dir / "modmake"
    plots_dir = output_dir / "plots"
    modmake_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    csm_T = args.csm_T if args.csm_T is not None else args.wind_T
    csm_T_min = (
        args.csm_T_min if args.csm_T_min is not None else args.ejecta_T_floor
    )
    csm_comp = args.csm_composition if args.csm_composition else args.shell_composition

    print(f"\n=== Preparing SN_B v2 initial model (rho ~ r^-n CSM) ===")
    print(f"  Source directory: {source_dir}")
    print(f"  Output directory: {output_dir}")
    print(
        f"  Shell: M={args.shell_mass} Msun, v={args.shell_velocity} km/s, "
        f"X={args.shell_xh}, Y={args.shell_yhe}, "
        f"composition={args.shell_composition}"
    )
    if args.no_wind:
        print(
            f"  Mode: immediate shell-ejecta contact (--no-wind); "
            f"CSM skipped"
        )
    elif args.no_ejecta:
        print(f"  Mode: shell + power-law CSM only (--no-ejecta)")
    else:
        print(f"  Mode: shell + power-law CSM + SN_A ejecta envelope")
    if not args.no_wind:
        age_note = (
            f"{args.csm_age_days} d"
            if args.csm_age_days > 0
            else f"auto (= SN_A age {args.source_age + args.t_start} d)"
        )
        if args.no_ejecta:
            norm_note = (
                f"rho15={args.csm_rho15:.3e}"
                if args.csm_rho15 > 0
                else (
                    f"w={args.csm_w:.3e}" if args.csm_w > 0 else "default rho15=6e-15"
                )
            )
        else:
            norm_note = (
                "match to ejecta"
                if args.csm_density_match
                else (
                    f"rho15={args.csm_rho15:.3e}"
                    if args.csm_rho15 > 0
                    else (f"w={args.csm_w:.3e}" if args.csm_w > 0 else "MATCH REQUIRED")
                )
            )
        if args.beta_wind:
            vinf_note = (
                f"{args.wind_velocity:.0f} km/s"
                if args.wind_velocity > 0
                else "auto (ejecta edge)"
            )
            mdot_note = (
                f"{args.wind_mdot:.3e} Msun/yr"
                if args.wind_mdot > 0
                else "ejecta-matched"
            )
            print(
                f"  CSM: beta wind v_inf={vinf_note}, "
                f"beta={args.wind_beta:g}, Mdot={mdot_note}, "
                f"T={csm_T:.0f} K, comp={csm_comp}"
            )
        else:
            print(
                f"  CSM: rho ~ r^-{args.csm_powerlaw_n:g}, t_cs={age_note}, "
                f"norm={norm_note}, T={csm_T:.0f} K, comp={csm_comp}"
            )
        if not args.beta_wind and (args.wind_mdot > 0 or args.wind_velocity > 0):
            print(
                f"  NOTE: --wind-mdot/--wind-velocity only apply to "
                f"--beta-wind mode and are ignored"
            )
    print(f"  Source age: {args.source_age} days")
    print(
        f"  STELLA target: {args.target_days} days, "
        f"CMFGEN snapshot: {args.snapshot_day} days"
    )

    # --- Step 1: source ejecta ---
    ejecta = None
    if not args.no_ejecta:
        effective_age = args.source_age + args.t_start
        print(
            f"\n--- Step 1: Reading source ejecta at {effective_age} days "
            f"({args.source_age}d at SN_B + {args.t_start}d t_start) ---"
        )
        ejecta = read_source_ejecta(
            source_dir,
            source_age=args.source_age,
            t_start_days=args.t_start,
            force_homologous_velocity=args.force_homologous_velocity,
            smooth_ejecta=not args.no_smooth_ejecta,
            smooth_window=args.smooth_window,
            smooth_polyorder=args.smooth_polyorder,
            comp_window=args.comp_smooth_window,
        )
        print(
            f"  Ejecta: {len(ejecta['r'])} zones, "
            f"r=[{ejecta['r'].min():.3e}, {ejecta['r'].max():.3e}] cm"
        )
        print(
            f"  Ejecta: v=[{ejecta['v'].min():.1f}, " f"{ejecta['v'].max():.1f}] km/s"
        )
        print(
            f"  Ejecta: rho=[{ejecta['rho'].min():.3e}, "
            f"{ejecta['rho'].max():.3e}] g/cc"
        )

        ejecta["v"] = ejecta["v"] * 1e5  # km/s -> cm/s

        # Ejecta truncation options (same as v1)
        if args.ejecta_rho_min > 0:
            mask = ejecta["rho"] >= args.ejecta_rho_min
            n_before = len(ejecta["r"])
            if n_before - mask.sum() > 0:
                r_max_before = ejecta["r"].max()
                for key in ejecta:
                    if (
                        isinstance(ejecta[key], np.ndarray)
                        and len(ejecta[key]) == n_before
                    ):
                        ejecta[key] = ejecta[key][mask]
                print(
                    f"  Ejecta truncation (rho_min={args.ejecta_rho_min:.1e}): "
                    f"{n_before} -> {mask.sum()} zones, "
                    f"r_max {r_max_before:.3e} -> {ejecta['r'].max():.3e} cm"
                )
        if args.ejecta_r_max_cm > 0:
            mask = ejecta["r"] <= args.ejecta_r_max_cm
            n_before = len(ejecta["r"])
            if n_before - mask.sum() > 0:
                for key in ejecta:
                    if (
                        isinstance(ejecta[key], np.ndarray)
                        and len(ejecta[key]) == n_before
                    ):
                        ejecta[key] = ejecta[key][mask]
                print(
                    f"  Ejecta radius truncation: r_max -> "
                    f"{ejecta['r'].max():.3e} cm ({mask.sum()} zones)"
                )
        if args.ejecta_v_max_kms > 0:
            v_max_cms = args.ejecta_v_max_kms * 1e5
            # Radius cut at the first outward v = v_max crossing.  A
            # plain v <= v_max mask also keeps the decelerated pile-up
            # tail at the outer edge, leaving a hole in r and a
            # non-monotonic velocity for the combiner to interpolate.
            order = np.argsort(ejecta["r"])
            v_asc = ejecta["v"][order]
            over = np.where(v_asc > v_max_cms)[0]
            if len(over):
                keep = np.zeros(len(v_asc), dtype=bool)
                keep[order[: int(over[0])]] = True
                n_before = len(ejecta["r"])
                for key in ejecta:
                    if (
                        isinstance(ejecta[key], np.ndarray)
                        and len(ejecta[key]) == n_before
                    ):
                        ejecta[key] = ejecta[key][keep]
                print(
                    f"  Ejecta velocity truncation: v_max -> "
                    f"{ejecta['v'].max()/1e5:.0f} km/s, r_max -> "
                    f"{ejecta['r'].max():.3e} cm ({keep.sum()} zones)"
                )
        if args.ejecta_r_inner_cm > 0:
            mask = ejecta["r"] >= args.ejecta_r_inner_cm
            n_before = len(ejecta["r"])
            if n_before - mask.sum() > 0:
                r_min_before = ejecta["r"].min()
                for key in ejecta:
                    if (
                        isinstance(ejecta[key], np.ndarray)
                        and len(ejecta[key]) == n_before
                    ):
                        ejecta[key] = ejecta[key][mask]
                print(
                    f"  Ejecta inner-edge truncation: r_inner "
                    f"{r_min_before:.3e} -> {ejecta['r'].min():.3e} cm "
                    f"({mask.sum()} zones kept)"
                )

    # --- Step 2: H/He shell ---
    r_ej_inner = (ejecta["r"].min() * 0.99) if ejecta is not None else None
    if args.no_wind:
        print(f"\n--- Step 2: Creating H/He shell (immediate contact) ---")
        shell = create_hhe_shell(
            args.shell_mass,
            args.shell_velocity,
            args.shell_xh,
            args.shell_yhe,
            n_rho=args.n_rho,
            v0_kms=args.v0,
            n_T=args.n_T,
            T_interface_K=args.T_interface,
            T_max_K=args.T_max,
            t_start_days=args.t_start,
            r_outer_override=r_ej_inner,
            composition=args.shell_composition,
        )
    else:
        print(f"\n--- Step 2: Creating H/He shell (Dessart structure) ---")
        shell = create_hhe_shell(
            args.shell_mass,
            args.shell_velocity,
            args.shell_xh,
            args.shell_yhe,
            n_rho=args.n_rho,
            v0_kms=args.v0,
            n_T=args.n_T,
            T_interface_K=args.T_interface,
            T_max_K=args.T_max,
            t_start_days=args.t_start,
            composition=args.shell_composition,
        )

    # Grid inversion check (same as v1)
    if not args.no_wind and not args.no_ejecta:
        r_shell_outer = shell["r_outer"]
        r_ej_inner_check = ejecta["r"].min()
        if r_shell_outer >= r_ej_inner_check:
            raise ValueError(
                f"Grid inversion: shell outer radius has overtaken the "
                f"ejecta inner radius.\n"
                f"  r_shell_outer = {r_shell_outer:.3e} cm, "
                f"r_ejecta_inner = {r_ej_inner_check:.3e} cm\n"
                f"  Reduce --t-start below "
                f"{r_ej_inner_check / (args.shell_velocity * 1e5) / 86400.0:.2f} d"
            )
        gap_dex = np.log10(r_ej_inner_check) - np.log10(r_shell_outer)
        print(f"  Grid inversion check: PASSED (gap = {gap_dex:.2f} dex)")

    # --- Step 3: Power-law CSM + combine ---
    csm = None
    if args.no_wind:
        print(f"\n--- Step 3: SKIPPING CSM (--no-wind, immediate contact) ---")
        combined = combine_shell_ejecta(
            shell,
            ejecta,
            nzon_total=args.nzon,
            se_transition_dex=args.se_transition_dex,
            immediate_contact=True,
            ejecta_T_floor=args.ejecta_T_floor,
            min_dlogr=args.min_dlogr,
            max_dlogr=args.max_dlogr,
        )
    elif args.no_ejecta:
        print(f"\n--- Step 3: Creating power-law CSM (no ejecta) ---")
        r_csm_outer = (
            args.csm_r_outer_cm if args.csm_r_outer_cm > 0 else args.wind_r_outer_cm
        )
        if r_csm_outer <= 0:
            raise ValueError(
                "--no-ejecta mode needs a CSM outer radius: set "
                "--csm-r-outer-cm."
            )
        csm_r_inner = shell["r_outer"] * 1.01
        if csm_r_inner >= r_csm_outer:
            raise ValueError(
                f"CSM grid inversion: r_inner={csm_r_inner:.3e} >= "
                f"r_outer={r_csm_outer:.3e} cm."
            )
        t_cs = (
            args.csm_age_days
            if args.csm_age_days > 0
            else args.source_age + args.t_start
        )
        rho15 = args.csm_rho15 if args.csm_rho15 > 0 else None
        w_gcm = args.csm_w if args.csm_w > 0 else None
        if rho15 is None and w_gcm is None:
            rho15 = 6.0e-15
            print(
                f"  No density normalisation given; defaulting to the "
                f"explicit rho(1e15) = {rho15:.1e} g/cc"
            )
        csm = create_powerlaw_csm(
            csm_r_inner,
            r_csm_outer,
            slope=-args.csm_powerlaw_n,
            age_days=t_cs,
            rho15=rho15,
            w_gcm=w_gcm,
            T_K=csm_T,
            composition=csm_comp,
            xh=args.shell_xh,
            yhe=args.shell_yhe,
            ua0_kms=args.csm_ua0,
            r_acc_cm=shell["r_outer"],
            cutoff_r_cm=args.csm_cutoff_r_cm,
            cutoff_n=args.csm_cutoff_n,
            T_gamma=args.csm_T_gamma,
            T_min=csm_T_min,
            match_shell=args.csm_match_shell,
        )
        print_csm_diagnostics(csm, csm_r_inner, r_csm_outer)
        combined = combine_shell_csm(
            shell,
            csm,
            nzon_total=args.nzon,
            sw_transition_dex=args.sw_transition_dex,
            sc_v_sharp_zones=args.sc_v_sharp_zones,
            min_dlogr=args.min_dlogr,
            max_dlogr=args.max_dlogr,
        )
    else:
        print(f"\n--- Step 3: Creating power-law CSM + combining ---")
        csm_r_inner = shell["r_outer"] * 1.01
        csm_r_outer = ejecta["r"].min() * 0.99
        if csm_r_inner >= csm_r_outer:
            raise ValueError(
                f"CSM grid inversion: r_inner={csm_r_inner:.3e} >= "
                f"r_outer={csm_r_outer:.3e} cm. Reduce --t-start."
            )
        t_cs = (
            args.csm_age_days
            if args.csm_age_days > 0
            else args.source_age + args.t_start
        )
        if args.csm_age_days <= 0:
            print(
                f"  t_cs = {t_cs:.1f} d (auto = SN_A ejecta age; "
                f"velocity-continuous at the interface)"
            )
        rho15 = args.csm_rho15 if args.csm_rho15 > 0 else None
        w_gcm = args.csm_w if args.csm_w > 0 else None
        explicit_norm = (rho15 is not None) or (w_gcm is not None)
        if args.csm_match_shell and explicit_norm:
            parser.error(
                "--csm-match-shell cannot be combined with "
                "--csm-rho15/--csm-w (it sets the normalisation itself)"
            )
        density_matched = (
            args.csm_density_match
            and not explicit_norm
            and not args.csm_match_shell
        )
        if args.beta_wind:
            # v_inf = ejecta inner-edge velocity keeps the interface
            # velocity-continuous, like t_cs=auto does for the power law
            ej_order = np.argsort(ejecta["r"])
            vinf_kms = (
                args.wind_velocity
                if args.wind_velocity > 0
                else float(ejecta["v"][ej_order][0]) / 1e5
            )
            law, mdot = "beta", args.wind_mdot
        else:
            vinf_kms, law, mdot = 0.0, "powerlaw", 0.0
        csm = create_powerlaw_csm(
            csm_r_inner,
            csm_r_outer,
            slope=-args.csm_powerlaw_n,
            age_days=t_cs,
            rho15=rho15,
            w_gcm=w_gcm,
            T_K=csm_T,
            composition=csm_comp,
            xh=args.shell_xh,
            yhe=args.shell_yhe,
            ua0_kms=args.csm_ua0,
            r_acc_cm=shell["r_outer"],
            cutoff_r_cm=args.csm_cutoff_r_cm,
            cutoff_n=args.csm_cutoff_n,
            T_gamma=args.csm_T_gamma,
            T_min=csm_T_min,
            match_shell=args.csm_match_shell,
            law=law,
            mdot_msun_yr=mdot,
            vinf_kms=vinf_kms,
            rstar_cm=args.wind_rstar_cm,
            beta=args.wind_beta,
            clumping_D=args.wind_clumping,
        )
        combined = combine_shell_csm_ejecta(
            shell,
            csm,
            ejecta,
            nzon_total=args.nzon,
            sw_transition_dex=args.sw_transition_dex,
            sc_v_sharp_zones=args.sc_v_sharp_zones,
            cd_transition_zones=args.cd_transition_zones,
            ce_v_ramp_zones=args.ce_v_ramp_zones,
            v_blend_zones=args.v_blend_zones,
            ejecta_T_floor=args.ejecta_T_floor,
            ejecta_T_const=args.ejecta_T_const,
            min_dlogr=args.min_dlogr,
            max_dlogr=args.max_dlogr,
            density_matched=density_matched,
        )
        print_csm_diagnostics(csm, csm_r_inner, csm_r_outer)

    print(f"  Combined: {len(combined['r'])} zones")
    print(f"  r: [{combined['r'].min():.3e}, {combined['r'].max():.3e}] cm")
    print(
        f"  v: [{combined['v'].min()/1e5:.1f}, " f"{combined['v'].max()/1e5:.1f}] km/s"
    )
    print(f"  rho: [{combined['rho'].min():.3e}, " f"{combined['rho'].max():.3e}] g/cc")

    # --- Outer boundary cleanup (ejecta mode only, same as v1) ---
    if not args.no_ejecta:
        combined = clean_combined_outer_boundary(
            combined, ej_start_idx=combined.get("idx_ej")
        )
        # Re-apply the outer velocity clip: the cleaner extrapolates the
        # edge velocity with a polynomial fit and enforces a monotone
        # increase, which re-inflates a clipped ejecta edge past
        # ejecta_v_max.  Cap those zones at v_max rather than dropping
        # them - dropping would shrink both the zone count and r_max.
        if args.ejecta_v_max_kms > 0:
            v_cap = args.ejecta_v_max_kms * 1e5
            ej0 = combined.get("idx_ej") or 0
            over = np.where(combined["v"][ej0:] > v_cap)[0]
            if len(over):
                combined["v"][ej0 + over] = v_cap
                print(
                    f"  Re-applied ejecta v_max clip: capped {len(over)} "
                    f"re-inflated edge zones at {args.ejecta_v_max_kms:.0f} km/s"
                )

    # --- Grid validation ---
    validate_stella_grid(
        combined["r"],
        min_dlogr=args.min_dlogr,
        max_dlogr=args.max_dlogr,
        strict=args.strict_grid,
    )

    # --- Step 4: write STELLA inputs ---
    print(f"\n--- Step 4: Writing STELLA input files ---")
    hyd_path = modmake_dir / f"{model_name}.hyd"
    abn_path = modmake_dir / f"{model_name}.abn"
    xni_path = modmake_dir / f"{model_name}.xni"
    dat_path = modmake_dir / f"{model_name}.dat"
    ip_path = modmake_dir / "input_profile"

    write_x_hyd(str(hyd_path), combined, t_start_days=args.t_start)
    write_x_abn(str(abn_path), combined)
    write_x_xni(str(xni_path), len(combined["r"]))
    if args.explosion_energy > 0:
        print(
            f"  WARNING: --explosion-energy={args.explosion_energy} adds "
            f"a thermal bomb; SN_B is velocity-driven, keep it 0."
        )
    write_x_dat(
        str(dat_path),
        model_name,
        args.target_days,
        args.snapshot_day,
        explosion_energy=args.explosion_energy,
        nstmax=args.nstmax,
        aq=args.aq,
        fit_tau=args.fit_tau,
        meth=args.meth,
        temp_floor=args.temp_floor,
        floor_rad=args.floor_rad,
        hmin=args.hmin,
        hmax=args.hmax,
        t_start_days=args.t_start,
        tcura_days=args.tcura_days,
        r_max_cm=combined["r"].max(),
        v_max_cm_s=combined["v"].max(),
        wacc_r=args.wacc_r,
        wacc_v=args.wacc_v,
        wacc_t=args.wacc_t,
        wacc_rad=args.wacc_rad,
        tau_tol=args.tau_tol,
        iout=args.iout,
        disable_scattering=args.disable_scattering,
    )
    write_v2_input_profile(str(ip_path), combined)

    import shutil

    shutil.copy2(str(ip_path), str(stella_dir / "input_profile"))

    # --- Step 5: plots ---
    print(f"\n--- Step 5: Generating plots ---")
    if ejecta is not None and "raw_ejecta" in ejecta and "smoothed_ejecta" in ejecta:
        try:
            plot_ejecta_smoothing(
                ejecta["raw_ejecta"],
                ejecta["smoothed_ejecta"],
                str(plots_dir / "ejecta_smoothing.png"),
                source_age=args.source_age,
            )
        except Exception as e:
            print(f"  Ejecta smoothing plot failed: {e}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        ax_rho_dens, ax_r = axes[0]
        ax_v, ax_rho = axes[1]

        r_s = combined["r"]
        v_k = combined["v"] / 1e5

        r_shell_lo = shell["r_inner"]
        r_shell_hi = shell["r_outer"]
        r_max_cm = combined["r"][-1]

        # Component masses for the legend labels (4 pi r^2 rho dr per zone)
        dm = 4.0 * np.pi * r_s**2 * combined["rho"] * np.gradient(r_s) / MSUN
        r_split = (
            r_ej_inner
            if (r_ej_inner is not None and not args.no_ejecta)
            else r_max_cm
        )
        m_shell = dm[r_s <= r_shell_hi].sum()
        m_csm = dm[(r_s > r_shell_hi) & (r_s <= r_split)].sum()
        m_ej = dm[r_s > r_split].sum()
        # Effective density exponent for the legend: the beta wind has
        # no fixed n (slope -> -2 only at r >> R*), so use the mean
        # local log-slope across the CSM span.
        csm_mask = (r_s > r_shell_hi) & (r_s <= r_split)
        if csm is None:
            n_eff = 0.0
            vel_law = ""
        elif csm.get("law") == "beta":
            n_eff = -float(
                np.mean(_csm_log_slope(csm, r_s[csm_mask]))
            )
            vel_law = "beta law"
        else:
            n_eff = abs(csm["slope"])
            vel_law = "u=r/t_cs"
        csm_label = f"rho_wind~r^-{n_eff:.1f} ({m_csm:.2f} Msun)"

        ax_rho_dens.plot(r_s, combined["rho"], "k-", linewidth=1.5)
        ax_rho_dens.axvspan(
            r_shell_lo,
            r_shell_hi,
            color="red",
            alpha=0.15,
            label=f"shell ({m_shell:.2f} Msun)",
        )
        if args.no_wind:
            ax_rho_dens.axvspan(
                r_shell_hi,
                r_ej_inner,
                color="orange",
                alpha=0.10,
                label=f"buffer ({m_csm:.2f} Msun)",
            )
            ax_rho_dens.axvspan(
                r_ej_inner,
                r_max_cm,
                color="green",
                alpha=0.15,
                label=f"ejecta ({m_ej:.2f} Msun)",
            )
        elif args.no_ejecta:
            ax_rho_dens.axvspan(
                r_shell_hi, r_max_cm, color="cyan", alpha=0.15, label=csm_label
            )
        else:
            ax_rho_dens.axvspan(
                r_shell_hi, r_ej_inner, color="cyan", alpha=0.15, label=csm_label
            )
            ax_rho_dens.axvspan(
                r_ej_inner,
                r_max_cm,
                color="green",
                alpha=0.15,
                label=f"SN_A ejecta ({m_ej:.2f} Msun)",
            )
        ax_rho_dens.set_xlabel("Radius (cm)")
        ax_rho_dens.set_ylabel("Density (g/cm^3)")
        ax_rho_dens.set_title("Density Profile")
        ax_rho_dens.set_xscale("log")
        ax_rho_dens.set_yscale("log")
        ax_rho_dens.grid(True, alpha=0.3)
        ax_rho_dens.legend(fontsize=9, loc="best")

        ax_r.plot(r_s, v_k, "b-", linewidth=1.5)
        ax_r.set_xlabel("Radius (cm)")
        ax_r.set_ylabel("Velocity (km/s)")
        ax_r.set_title("Velocity Profile")
        ax_r.set_xscale("log")
        ax_r.grid(True, alpha=0.3)
        ax_r.axvspan(
            r_shell_lo,
            r_shell_hi,
            color="red",
            alpha=0.10,
            label="shell",
        )
        ax_r.axvspan(
            r_shell_hi,
            r_split,
            color="cyan",
            alpha=0.10,
            label=f"v_wind {vel_law}".rstrip(),
        )
        if r_ej_inner is not None and not args.no_ejecta:
            ax_r.axvspan(
                r_ej_inner,
                r_max_cm,
                color="green",
                alpha=0.10,
                label="SN_A ejecta",
            )
        ax_r.axvline(r_shell_hi, color="r", linestyle="--", alpha=0.5)
        if r_ej_inner is not None and not args.no_ejecta:
            ax_r.axvline(r_ej_inner, color="g", linestyle="--", alpha=0.5)
        ax_r.legend(fontsize=9)

        ax_v.plot(r_s, combined["T"], "r-", linewidth=1.5)
        ax_v.set_xlabel("Radius (cm)")
        ax_v.set_ylabel("Temperature (K)")
        ax_v.set_title("Temperature Profile")
        ax_v.set_xscale("log")
        ax_v.set_yscale("log")
        ax_v.grid(True, alpha=0.3)

        for key, color, label in [
            ("fH", "blue", "H"),
            ("fHe", "orange", "He"),
            ("fO", "green", "O"),
            ("fSi", "red", "Si"),
            ("fFe", "purple", "Fe"),
            ("fCa", "brown", "Ca"),
        ]:
            if np.any(combined[key] > 1e-10):
                ax_rho.plot(r_s, combined[key], color=color, label=label, linewidth=1.5)
        ax_rho.set_xlabel("Radius (cm)")
        ax_rho.set_ylabel("Mass Fraction")
        ax_rho.set_title("Composition Profile")
        ax_rho.set_xscale("log")
        ax_rho.set_ylim(-0.05, 1.05)
        ax_rho.legend(fontsize=8)
        ax_rho.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(
            str(plots_dir / "sn_b_v2_initial.png"), dpi=150, bbox_inches="tight"
        )
        plt.close()
        print(f"  Saved: {plots_dir}/sn_b_v2_initial.png")
    except Exception as e:
        print(f"  Plot generation failed: {e}")

    print(f"\n=== SN_B v2 preparation complete ===")
    print(f"  STELLA input files in: {modmake_dir}")
    print(f"  Run STELLA with: python scripts/run_model_grid.py {output_dir}/config.yaml")


if __name__ == "__main__":
    main()
