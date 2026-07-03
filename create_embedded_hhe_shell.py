#!/usr/bin/env python3
"""
Create initial HERACLES profile for embedded H/He-rich supernova model
Based on arXiv:1503.05463 - Numerical simulations of super-luminous supernovae of type IIn
This script sets up a multi-shell structure:
1. Inner region (R < Rt): Ejecta A - H/He-rich progenitor envelope (M_ejecta_A, fx_H_A)
2. Middle region (Rt < R < 1e16 cm): H/He-rich CSM wind
3. Outer region (R > 1e16 cm): Ejecta B - metal-rich stripped envelope SN (M_ejecta_B, fx_H_B)
Usage:
    python create_embedded_hhe_shell.py --days 400 --output profile_embedded_hhe.dat --plots
"""
import numpy as np
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
def create_embedded_hhe_profile(days_ejecta_A=11.6, days_ejecta_B=400, n_cells=10000,
                                 fx_Fe_ejecta=0.03, fx_Si_ejecta=0.05, fx_He_ejecta=0.40, fx_H_ejecta=0.01,
                                 M_ejecta_A=9.8, Mdot_wind=1e-3, M_ejecta_B=0.1,
                                 E_kin_ejecta_A=1e51, E_kin_ejecta_B=0,
                                 fx_H_A=0.70, fx_He_A=0.28,
                                 fx_O_A=0.009, fx_Si_A=0.0007, fx_Ca_A=0.0002, fx_Fe_A=0.0013,
                                 fx_H_wind=0.70, fx_He_wind=0.28,
                                 fx_H_B=0.01, fx_He_B=0.40,
                                 fx_O_B=0.15, fx_Si_B=0.05, fx_Ca_B=0.02, fx_Fe_B=0.03, fx_Ni_B=0.02,
                                 temp_profile_type='adiabatic', N_T_ejecta=0.6, T_inner_max=100000,
                                 transition_width_fraction=0.2, transition_power_law_exponent=0.5,
                                 sigmoid_k=20,
                                 stationary_ejecta_A=False,
                                 R_star_cm=None, n_shell=2.0, wind_velocity_type='beta_law',
                                 V_inf_kms=100.0, beta_wind=1.0, wind_time_days=400.0,
                                 R_CSM_MAX=None, V_ejecta_B_max_kms=10000.0,
                                 N_rho_ejecta_B=7):
    """
    Create a multi-shell HERACLES profile with embedded H/He system.
    Structure (NEW LAYOUT):
    1. SN ejecta A (R < R_inner_shell): H/He-rich, stationary shell (RSG mode) or homologous expansion
    2. SN ejecta B (R_inner_shell < R < R_crossover): Metal-rich, stripped envelope
    3. Wind (R_crossover < R < R_CSM_MAX): H/He-rich wind (beta-law or homologous expansion)
    Where R_crossover is determined by velocity cutoff (V_ejecta_B_max_kms) instead of density comparison.
    Metal-rich ejecta B extends to V_ejecta_B_max_kms, then uniform wind is appended.
    Evolution timeline:
    1. Ejecta B evolves for days_ejecta_B days (initial stripped envelope)
    2. Ejecta A explodes and both evolve for days_ejecta_A days
    Parameters:
    -----------
    days_ejecta_A : float
        H/He-rich ejecta A age in days (default: 11.6 for model X)
    days_ejecta_B : float
        Metal-rich ejecta B initial evolution in days (default: 400)
    n_cells : int
        Number of radial cells
    E_kin_ejecta_A : float
        Kinetic energy of ejecta A in erg (default: 1e51)
    E_kin_ejecta_B : float
        Kinetic energy of ejecta B in erg (default: 0, meaning no separate kinetic energy constraint)
    fx_Fe_ejecta : float
        Iron mass fraction in outer ejecta (default: 0.03, SN2011dh-like)
    fx_Si_ejecta : float
        Silicon mass fraction in outer ejecta (default: 0.05, SN2011dh-like)
    fx_He_ejecta : float
        Helium mass fraction in outer ejecta (default: 0.40, SN2011dh-like Type IIb)
    fx_H_ejecta : float
        Hydrogen mass fraction in outer ejecta (default: 0.01, SN2011dh-like Type IIb)
    M_ejecta_A : float
        Mass of ejecta A in Msun (default: 9.8)
    Mdot_wind : float
        Mass loss rate for wind in Msun/yr (default: 1e-3)
    M_ejecta_B : float
        Mass of ejecta B in Msun (default: 0.1)
    fx_H_A : float
        Hydrogen mass fraction in ejecta A (default: 0.70)
    fx_He_A : float
        Helium mass fraction in ejecta A (default: 0.28)
    fx_H_wind : float
        Hydrogen mass fraction in wind (default: 0.70)
    fx_He_wind : float
        Helium mass fraction in wind (default: 0.28)
    fx_H_B : float
        Hydrogen mass fraction in ejecta B (default: 0.01)
    fx_He_B : float
        Helium mass fraction in ejecta B (default: 0.40)
    temp_profile_type : str
        Temperature profile type for wind: 'isothermal' or 'adiabatic' (default: 'adiabatic')
    N_T_ejecta : float
        Temperature power law exponent for ejecta (default: 0.6, more realistic than 0.4)
    T_inner_max : float
        Maximum temperature in inner ejecta in K (default: 100000 K)
    transition_width_fraction : float
        Fraction of wind region to use for velocity transition layer (default: 0.2)
    transition_power_law_exponent : float
        Power law exponent for velocity increase in transition layer (v ∝ r^α, default: 0.5)
    R_CSM_MAX : float, optional
        Outer boundary of wind region in cm (default: None, computed from ejecta B expansion)
    V_ejecta_B_max_kms : float
        Maximum velocity for metal-rich ejecta B in km/s (default: 10000.0)
        The metal-rich ejecta extends to this velocity, then wind is appended.
    wind_velocity_type : str
        Wind velocity profile type: 'beta_law' or 'homologous' (default: 'beta_law')
        - beta_law: v(r) = V_inf * (1 - R_star/r)^beta
        - homologous: v(r) = r / t_wind (where t_wind = wind_time_days)
    V_inf_kms : float
        Terminal wind velocity in km/s (for beta_law mode, default: 100.0)
    beta_wind : float
        Beta-law exponent (for beta_law mode, default: 1.0)
    wind_time_days : float
        Characteristic time for homologous expansion in days (for homologous mode, default: 400.0)
    Returns:
    --------
    dict : Dictionary containing r, rho, p, u, and mass fractions
    """
    # Physical constants
    Msun = 1.989e33  # g
    Mpc = 3.086e24  # cm
    km_s = 1e5  # cm/s
    day = 86400  # s
    def sigmoid(x, x0, k=10):
        """Sigmoid function for smooth transitions.
        Parameters:
        -----------
        x : array
            Input values (typically radius)
        x0 : float
            Center of transition
        k : float
            Steepness of transition (higher = sharper)
        Returns:
        --------
        array
            Sigmoid values between 0 and 1
        """
        return 1 / (1 + np.exp(-k * (np.log10(x) - np.log10(x0))))
    # Determine mode: RSG shell + beta-law wind, or legacy homologous ejecta A
    use_rsg_mode = R_star_cm is not None
    # Physical constants for density/velocity setup
    N_rho = N_rho_ejecta_B  # density power law exponent for ejecta B (default 7 for stripped SN envelope)
    N_T = N_T_ejecta
    T_CSM = 5000  # K constant wind temperature
    T_max_inner = T_inner_max
    # Convert Mdot from Msun/yr to g/s (only if wind is present)
    if Mdot_wind > 0:
        Mdot_CSM_in = Mdot_wind * Msun / (365.25 * day)  # g/s
    else:
        Mdot_CSM_in = 0  # No wind
    # Inner shell boundary
    M_ejecta_A_g = M_ejecta_A * Msun
    R_ejecta_min = 1e7  # cm (innermost grid point)
    if use_rsg_mode:
        # ---- RSG MODE ----
        # Inner shell: r < R_star  (stationary power-law density)
        # Wind:        R_star < r < R_cross  (beta-law or homologous velocity, rho = Mdot/4pi r^2 v)
        R_inner_shell = R_star_cm          # shell surface = RSG stellar radius
        V_inf = V_inf_kms * km_s           # terminal wind velocity in cm/s
        v_wind_floor = 1.0 * km_s          # floor to avoid rho divergence at r=R_star
        t_wind = wind_time_days * day      # characteristic time for homologous expansion
        print(f"  RSG mode: R_star = {R_star_cm:.2e} cm, wind_velocity_type = {wind_velocity_type}")
        if wind_velocity_type == 'beta_law':
            print(f"    Beta-law: V_inf = {V_inf_kms} km/s, beta = {beta_wind}")
        else:
            print(f"    Homologous: t_wind = {wind_time_days} days")
        print(f"    n_shell = {n_shell}")
    else:
        # ---- LEGACY MODE ----
        R_t = 1e15  # cm (ejecta/CSM interface radius)
        t_ejecta_A = days_ejecta_A * day
        t_ejecta_A_geom = t_ejecta_A if t_ejecta_A > 0 else 7 * day
        V_max = R_t / t_ejecta_A_geom
        V0 = V_max / 3
        v_CSM = 100 * km_s
        R_inner_shell = R_t  # shell surface = R_t for legacy mode
    # Compute t_ejecta_B_total for velocity profile (needed regardless of R_CSM_MAX)
    t_ejecta_B_total = (days_ejecta_B + (0 if use_rsg_mode else days_ejecta_A)) * day
    # R_CSM_max: outer edge from ejecta B homologous expansion or explicit R_CSM_MAX
    if R_CSM_MAX is not None:
        R_CSM_max = float(R_CSM_MAX)
        print(f"  Using explicit R_CSM_MAX = {R_CSM_max:.3e} cm")
    else:
        V_max_ejecta_B = 10000 * km_s
        R_CSM_max = V_max_ejecta_B * t_ejecta_B_total
        print(f"  Computed R_CSM_max from ejecta B expansion = {R_CSM_max:.3e} cm")
    # Wind velocity function for crossover search (used in both modes for rho_wind)
    if use_rsg_mode:
        if wind_velocity_type == 'beta_law':
            def wind_velocity(r_arr):
                x = np.maximum(1.0 - R_star_cm / r_arr, 0.0)
                return np.maximum(V_inf * x**beta_wind, v_wind_floor)
        elif wind_velocity_type == 'homologous':
            def wind_velocity(r_arr):
                # Homologous expansion: v(r) = r / t_wind
                # This gives v = 0 at r=0 and increases linearly with radius
                return r_arr / t_wind
        else:
            raise ValueError(f"Unknown wind_velocity_type: {wind_velocity_type}. Use 'beta_law' or 'homologous'")
    else:
        def wind_velocity(r_arr):
            return np.full_like(r_arr, v_CSM)
    # Find R_crossover: where ejecta B velocity reaches V_ejecta_B_max_kms
    # Metal-rich ejecta extends to this velocity, then wind is appended
    V_ejecta_B_max = V_ejecta_B_max_kms * km_s
    R_crossover = V_ejecta_B_max * t_ejecta_B_total
    R_crossover = max(R_crossover, R_inner_shell * 1.1)
    R_crossover = min(R_crossover, R_CSM_max * 0.9)
    print(f"  R_crossover = {R_crossover:.3e} cm (where ejecta B velocity = {V_ejecta_B_max_kms} km/s)")
    if not use_rsg_mode:
        R_transition_start = R_inner_shell + (1 - transition_width_fraction) * (R_crossover - R_inner_shell)
    # ---- Build radial grid (NEW LAYOUT: A → B → Wind) ----
    n_shell_cells = int(0.40 * n_cells)  # 40% inner shell (Ejecta A)
    n_ejecta_B    = int(0.40 * n_cells)  # 40% ejecta B (increased from 30%)
    n_wind        = n_cells - n_shell_cells - n_ejecta_B  # 20% wind (decreased from 30%)
    r_shell   = np.logspace(np.log10(R_ejecta_min), np.log10(R_inner_shell), n_shell_cells)
    r_ejecta_B = np.logspace(np.log10(R_inner_shell), np.log10(R_crossover), n_ejecta_B)
    r_wind     = np.logspace(np.log10(R_crossover), np.log10(R_CSM_max), n_wind)
    # Drop first point of each subsequent segment to avoid duplicate boundary points
    # (logspace includes both endpoints, so R_inner_shell and R_crossover would repeat)
    r = np.concatenate([r_shell, r_ejecta_B[1:], r_wind[1:]])
    # Initialize arrays
    rho = np.zeros_like(r)
    p   = np.zeros_like(r)
    u   = np.zeros_like(r)
    T   = np.zeros_like(r)
    fx_H  = np.zeros_like(r)
    fx_He = np.zeros_like(r)
    fx_O  = np.zeros_like(r)
    fx_Si = np.zeros_like(r)
    fx_Ca = np.zeros_like(r)
    fx_Fe = np.zeros_like(r)
    fx_Ni = np.zeros_like(r)
    # ================================================================
    # INNER SHELL  (r < R_inner_shell)
    # ================================================================
    if use_rsg_mode:
        # RSG envelope: Dessart et al. (2015) homologous ejecta profile
        # Flat inner core (v < v_t) + steep outer power law (v > v_t, rho ∝ v^{-n})
        # v_t and rho_0 are determined simultaneously from M_ejecta_A and E_kin_ejecta_A.
        #
        # In homologous flow r = v * t_A, so at epoch t_A (= days_ejecta_A * 86400):
        #   rho(r) = rho_0              for r < R_t  (= v_t * t_A)
        #   rho(r) = rho_0*(r/R_t)^-n  for r < r < R_star
        #
        # Normalisation integrals (at fixed t_A):
        #   M  = 4π t_A^3 [rho_0 * v_t^3/3  +  rho_0 * v_t^n * integral(v_t, v_max, v^(2-n) dv)]
        #   E  = (1/2) * 4π t_A^3 * t_A^2 [rho_0 * v_t^5/5  +  rho_0 * v_t^n * integral(v_t, v_max, v^(4-n) dv)]
        #
        # Solving: E/M = v_t^2 * f(n) where f(n) depends on n and v_max/v_t ratio.
        # For n > 5 and v_max >> v_t (outer power law dominates energy):
        #   f(n) = [1/5 + 1/(n-5)] / [1/3 + 1/(n-3)]  (approximate for v_max -> inf)
        # We solve numerically for v_t from E_kin/M_ej ratio, then set rho_0 from mass.
        n_s = n_shell
        t_A = days_ejecta_A * 86400.0
        v_max_A = R_star_cm / t_A   # maximum velocity = R_star / t_A
        if E_kin_ejecta_A > 0 and not stationary_ejecta_A:
            # Solve for v_t from E_kin/M ratio (only relevant for moving ejecta)
            # f(v_t) = E_kin/M_ej  where both integrals depend on v_t
            def _mass_energy(v_t):
                if abs(n_s - 3.0) < 1e-6:
                    I_m_env = v_t**3 * np.log(v_max_A / v_t)
                else:
                    I_m_env = v_t**n_s * (v_max_A**(3-n_s) - v_t**(3-n_s)) / (3 - n_s)
                I_m = v_t**3 / 3.0 + I_m_env
                if abs(n_s - 5.0) < 1e-6:
                    I_e_env = v_t**5 * np.log(v_max_A / v_t)
                else:
                    I_e_env = v_t**n_s * (v_max_A**(5-n_s) - v_t**(5-n_s)) / (5 - n_s)
                I_e = v_t**5 / 5.0 + I_e_env
                M_model = 4 * np.pi * t_A**3 * I_m   # in units of rho_0
                E_model = 0.5 * 4 * np.pi * t_A**5 * I_e
                return M_model, E_model
            # Binary search for v_t: E/M = E_kin_ejecta_A / M_ejecta_A_g
            target = E_kin_ejecta_A / M_ejecta_A_g
            v_lo, v_hi = v_max_A * 1e-6, v_max_A * 0.999
            for _ in range(80):
                v_mid = 0.5 * (v_lo + v_hi)
                M_, E_ = _mass_energy(v_mid)
                ratio = E_ / M_
                if ratio > target:
                    v_hi = v_mid
                else:
                    v_lo = v_mid
            v_t = 0.5 * (v_lo + v_hi)
            M_norm, _ = _mass_energy(v_t)
            rho_0_shell = M_ejecta_A_g / (4 * np.pi * t_A**3 * (v_t**3 / 3.0 +
                          (v_t**n_s * (v_max_A**(3-n_s) - v_t**(3-n_s)) / (3 - n_s)
                           if abs(n_s-3)>1e-6 else v_t**3 * np.log(v_max_A/v_t))))
        else:
            # No kinetic energy constraint: place v_t at 10% of v_max (Dessart fiducial)
            v_t = 0.1 * v_max_A
            if abs(n_s - 3.0) < 1e-6:
                I_m = v_t**3 / 3.0 + v_t**3 * np.log(v_max_A / v_t)
            else:
                I_m = v_t**3 / 3.0 + v_t**n_s * (v_max_A**(3-n_s) - v_t**(3-n_s)) / (3 - n_s)
            rho_0_shell = M_ejecta_A_g / (4 * np.pi * t_A**3 * I_m)
        R_t = v_t * t_A   # transition radius at epoch t_A
        # Density profile: flat core inside R_t, power law outside
        rho_ejecta = np.where(r <= R_t, rho_0_shell,
                              rho_0_shell * (r / R_t)**(-n_s))
        # Stationary in RSG mode — STELLA injects explosion energy via thermal bomb
        # BUT: give shell a small non-zero velocity to prevent infinite du/u
        # at the A/B boundary, which crashes STELLA's integrator with scattering ON.
        # u rises gradually from ~0 at center to match ejecta B velocity at R_star.
        # The thermal bomb still dominates — this is purely for numerical stability.
        u_boundary = R_inner_shell / t_ejecta_B_total  # ejecta B velocity at boundary
        u_ejecta = u_boundary * (r / R_inner_shell)**0.1
        u_floor = 1e3  # cm/s  (~0.01 km/s) absolute floor
        u_ejecta = np.maximum(u_ejecta, u_floor)
        print(f"  RSG shell (Dessart+2015): rho_0={rho_0_shell:.2e} g/cm3, "
              f"v_t={v_t/km_s:.0f} km/s, n={n_s}, M_shell={M_ejecta_A} Msun")
        print(f"    Shell velocity: floor={u_floor:.0f} cm/s, "
              f"boundary={u_boundary:.2e} cm/s (gradual power-law rise)")
    else:
        # Legacy: constant core + power-law envelope
        R0 = R_inner_shell / 3
        n_rho = N_rho
        if n_rho != 3:
            integral_inner = (R0**3 - R_ejecta_min**3) / 3
            integral_outer = R0**n_rho * (R_inner_shell**(3-n_rho) - R0**(3-n_rho)) / (3 - n_rho)
            rho_0 = M_ejecta_A_g / (4 * np.pi * (integral_inner + integral_outer))
        else:
            integral_inner = (R0**3 - R_ejecta_min**3) / 3
            integral_outer = R0**3 * np.log(R_inner_shell / R0)
            rho_0 = M_ejecta_A_g / (4 * np.pi * (integral_inner + integral_outer))
        rho_ejecta = np.where(r < R0, rho_0, rho_0 * (r / R0)**(-N_rho))
        u_ejecta = r / R_inner_shell * V_max
        if stationary_ejecta_A:
            u_ejecta = np.zeros_like(r)
            print(f"  stationary_ejecta_A=True: zeroing inner H/He ejecta velocities (r < R_t = {R_inner_shell:.2e} cm)")
    # Shell composition (solar-like H/He envelope - 5-element composition)
    fx_H_ejecta_arr  = fx_H_A  * np.ones_like(r)
    fx_He_ejecta_arr = fx_He_A * np.ones_like(r)
    fx_O_ejecta_arr  = fx_O_A * np.ones_like(r)
    fx_Si_ejecta_arr = fx_Si_A * np.ones_like(r)
    fx_Ca_ejecta_arr = fx_Ca_A * np.ones_like(r)
    fx_Fe_ejecta_arr = fx_Fe_A * np.ones_like(r)
    fx_Ni_ejecta_arr = 0.0 * np.ones_like(r)  # Ni is negligible in solar composition
    # ================================================================
    # EJECTA B  (R_inner_shell < r < R_crossover) - NEW LAYOUT
    # ================================================================
    M_ejecta_B_g = M_ejecta_B * Msun
    # t_ejecta_B_total already computed above (accounts for RSG vs legacy mode)
    u_ejecta_B = r / t_ejecta_B_total
    # Density profile for ejecta B: homologous expansion
    # In homologous flow r = v * t, so rho(r) = (1/t^3) * dM/dV where dV = 4*pi*v^2*dv
    # For a power-law mass distribution dM/dv ∝ v^(2 - N_rho):
    #   rho(r) ∝ r^(-N_rho) / t^3
    # Normalised by total mass between V_inner = R_inner_shell/t and V_max = R_crossover/t:
    #   M_B = integral_{R_inner}^{R_cross} rho(r) 4*pi*r^2 dr
    # => rho(r) = rho_0 * (r / R_inner_shell)^(-N_rho)
    # with rho_0 chosen so total mass = M_ejecta_B.
    # N_rho_ejecta_B controls the steepness (default 7 for stripped-envelope SN outer ejecta).
    # Using N_rho < 3 gives a profile that rises outward (unphysical), so floor at 3.
    N_rho_eff = max(N_rho, 3 + 1e-6)
    if abs(N_rho_eff - 3.0) > 1e-6:
        integral_B = R_inner_shell**N_rho_eff * \
                     (R_crossover**(3 - N_rho_eff) - R_inner_shell**(3 - N_rho_eff)) / (3 - N_rho_eff)
    else:
        integral_B = R_inner_shell**3 * np.log(R_crossover / R_inner_shell)
    rho_0_B_final = M_ejecta_B_g / (4 * np.pi * integral_B)
    rho_ejecta_B = rho_0_B_final * (r / R_inner_shell)**(-N_rho_eff)
    print(f"  Ejecta B: homologous rho ∝ r^(-{N_rho_eff:.0f}), "
          f"rho at R_inner={rho_0_B_final:.2e}, rho at R_cross="
          f"{rho_0_B_final*(R_crossover/R_inner_shell)**(-N_rho_eff):.2e} g/cm3")
    # Composition: metal-rich stripped envelope (ejecta B) - SN2011dh-like (5 elements + Ni)
    fx_H_ejecta_B  = fx_H_B * np.ones_like(r)
    fx_He_ejecta_B = fx_He_B * np.ones_like(r)
    fx_O_ejecta_B  = fx_O_B * np.ones_like(r)
    fx_Si_ejecta_B = fx_Si_B * np.ones_like(r)
    fx_Ca_ejecta_B = fx_Ca_B * np.ones_like(r)
    fx_Fe_ejecta_B = fx_Fe_B * np.ones_like(r)
    fx_Ni_ejecta_B = fx_Ni_B * np.ones_like(r)
    # ================================================================
    # WIND REGION  (R_star_cm < r < R_CSM_MAX) - Wind starts from RSG surface
    # ================================================================
    if use_rsg_mode:
        # Beta-law: v(r) = V_inf * (1 - R_star/r)^beta, floored at v_wind_floor
        # Wind starts from R_star_cm (RSG surface)
        x_wind = np.maximum(1.0 - R_star_cm / r, 0.0)
        u_wind = np.maximum(V_inf * x_wind**beta_wind, v_wind_floor)
    else:
        # Legacy: constant velocity with power-law ramp in outer wind
        u_wind = np.where(r < R_transition_start,
                          v_CSM,
                          v_CSM * (r / R_transition_start)**transition_power_law_exponent)
    # Wind density from mass continuity: rho = Mdot / (4 pi r^2 v)
    # Wind is calculated starting from R_star_cm (not from crossover)
    # Only calculate if wind is present
    if Mdot_wind > 0:
        rho_wind = Mdot_CSM_in / (4 * np.pi * r**2 * u_wind)
    else:
        rho_wind = np.zeros_like(r)  # No wind density
    # Wind composition (solar-like H/He - 5-element composition)
    fx_H_wind_arr  = fx_H_wind  * np.ones_like(r)
    fx_He_wind_arr = fx_He_wind * np.ones_like(r)
    fx_O_wind_arr  = fx_O_A * np.ones_like(r)
    fx_Si_wind_arr = fx_Si_A * np.ones_like(r)
    fx_Ca_wind_arr = fx_Ca_A * np.ones_like(r)
    fx_Fe_wind_arr = fx_Fe_A * np.ones_like(r)
    fx_Ni_wind_arr = 0.0 * np.ones_like(r)  # Ni is negligible in solar composition
    # ================================================================
    # SMOOTH TRANSITION LAYERS (NEW LAYOUT: A → B → Wind)
    # Wind is calculated from R_star_cm (RSG surface) outward
    # Physical structure:
    # - Shell (R_ejecta_min to R_inner_shell): stationary RSG shell
    # - Ejecta B (R_inner_shell to R_crossover): metal-rich ejecta dominates wind
    # - Wind (R_crossover to R_CSM_MAX): wind dominates ejecta B
    # Note: Wind density is calculated from R_star_cm, but only used where it dominates
    # ================================================================
    # Create masks for each region
    mask_shell = r < R_inner_shell
    # Wind region only exists if Mdot_wind > 0
    if Mdot_wind > 0:
        # Use density comparison to determine crossover:
        # wind dominates where rho_wind > rho_ejecta_B
        # This is physically correct and avoids vacuum zones
        wind_dominates = rho_wind > rho_ejecta_B
        mask_wind = (r >= R_inner_shell) & wind_dominates
        mask_ejecta_B = (r >= R_inner_shell) & ~wind_dominates
    else:
        # No wind - ejecta B extends to outer edge
        mask_wind = np.zeros_like(r, dtype=bool)
        mask_ejecta_B = r >= R_inner_shell
    # Blending at A/B boundary (R_inner_shell) using sigmoid for smooth transition
    delta_R_AB = 0.02 * R_inner_shell
    mask_blend_AB = (r >= R_inner_shell - delta_R_AB) & (r < R_inner_shell + delta_R_AB)
    # Use sigmoid for smooth transition (sigmoid_k controls steepness)
    blend_alpha_AB = sigmoid(r, R_inner_shell, k=sigmoid_k)
    log_rho_shell = np.log(rho_ejecta)
    log_rho_ejecta_B = np.log(rho_ejecta_B)
    rho_blend_AB = np.exp(log_rho_shell * (1.0 - blend_alpha_AB) + log_rho_ejecta_B * blend_alpha_AB)
    u_blend_AB = u_ejecta * (1.0 - blend_alpha_AB) + u_ejecta_B * blend_alpha_AB
    fx_H_blend_AB  = fx_H_ejecta_arr  * (1.0 - blend_alpha_AB) + fx_H_ejecta_B  * blend_alpha_AB
    fx_He_blend_AB = fx_He_ejecta_arr * (1.0 - blend_alpha_AB) + fx_He_ejecta_B * blend_alpha_AB
    fx_O_blend_AB  = fx_O_ejecta_arr  * (1.0 - blend_alpha_AB) + fx_O_ejecta_B  * blend_alpha_AB
    fx_Si_blend_AB = fx_Si_ejecta_arr * (1.0 - blend_alpha_AB) + fx_Si_ejecta_B * blend_alpha_AB
    fx_Ca_blend_AB = fx_Ca_ejecta_arr * (1.0 - blend_alpha_AB) + fx_Ca_ejecta_B * blend_alpha_AB
    fx_Fe_blend_AB = fx_Fe_ejecta_arr * (1.0 - blend_alpha_AB) + fx_Fe_ejecta_B * blend_alpha_AB
    fx_Ni_blend_AB = fx_Ni_ejecta_arr * (1.0 - blend_alpha_AB) + fx_Ni_ejecta_B * blend_alpha_AB
    # Blending at B/Wind boundary - only if wind exists
    # Use actual density-based crossover radius, not the fixed R_crossover
    if Mdot_wind > 0:
        # Find actual crossover radius where wind density = ejecta B density
        wind_dom_indices = np.where((r >= R_inner_shell) & (rho_wind > rho_ejecta_B))[0]
        R_crossover_actual = r[wind_dom_indices[0]] if len(wind_dom_indices) > 0 else R_crossover
        delta_R_BW = 0.02 * R_crossover_actual
        mask_blend_BW = (r >= R_crossover_actual - delta_R_BW) & (r < R_crossover_actual + delta_R_BW)
        # Use sigmoid for smooth transition (sigmoid_k controls steepness)
        blend_alpha_BW = sigmoid(r, R_crossover_actual, k=sigmoid_k)
        log_rho_ejecta_B_2 = np.log(rho_ejecta_B)
        log_rho_wind_2 = np.log(rho_wind)
        rho_blend_BW = np.exp(log_rho_ejecta_B_2 * (1.0 - blend_alpha_BW) + log_rho_wind_2 * blend_alpha_BW)
        u_blend_BW = u_ejecta_B * (1.0 - blend_alpha_BW) + u_wind * blend_alpha_BW
    else:
        # No wind - no blending needed
        mask_blend_BW = np.zeros_like(r, dtype=bool)
        rho_blend_BW = rho_ejecta_B
        u_blend_BW = u_ejecta_B
    # Composition blending at B/Wind boundary - only if wind exists
    if Mdot_wind > 0:
        fx_H_blend_BW  = fx_H_ejecta_B  * (1.0 - blend_alpha_BW) + fx_H_wind_arr  * blend_alpha_BW
        fx_He_blend_BW = fx_He_ejecta_B * (1.0 - blend_alpha_BW) + fx_He_wind_arr * blend_alpha_BW
        fx_O_blend_BW  = fx_O_ejecta_B  * (1.0 - blend_alpha_BW) + fx_O_wind_arr  * blend_alpha_BW
        fx_Si_blend_BW = fx_Si_ejecta_B * (1.0 - blend_alpha_BW) + fx_Si_wind_arr * blend_alpha_BW
        fx_Ca_blend_BW = fx_Ca_ejecta_B * (1.0 - blend_alpha_BW) + fx_Ca_wind_arr * blend_alpha_BW
        fx_Fe_blend_BW = fx_Fe_ejecta_B * (1.0 - blend_alpha_BW) + fx_Fe_wind_arr * blend_alpha_BW
        fx_Ni_blend_BW = fx_Ni_ejecta_B * (1.0 - blend_alpha_BW) + fx_Ni_wind_arr * blend_alpha_BW
    else:
        # No wind - use ejecta B composition
        fx_H_blend_BW  = fx_H_ejecta_B
        fx_He_blend_BW = fx_He_ejecta_B
        fx_O_blend_BW  = fx_O_ejecta_B
        fx_Si_blend_BW = fx_Si_ejecta_B
        fx_Ca_blend_BW = fx_Ca_ejecta_B
        fx_Fe_blend_BW = fx_Fe_ejecta_B
        fx_Ni_blend_BW = fx_Ni_ejecta_B
    # Apply blending with priority: BW > AB > base regions
    rho = np.where(mask_blend_BW, rho_blend_BW,
                   np.where(mask_blend_AB, rho_blend_AB,
                            np.where(mask_shell, rho_ejecta,
                                     np.where(mask_ejecta_B, rho_ejecta_B, rho_wind))))
    # Velocity: homologous (u = r/t) is enforced everywhere in r >= R_inner_shell
    # up to R_crossover (the velocity-based boundary). Beyond R_crossover, use wind velocity
    # only if wind is present (Mdot_wind > 0). With no wind, continue homologous expansion.
    mask_ejecta_B_vel = (r >= R_inner_shell) & (r <= R_crossover)
    if Mdot_wind > 0:
        u_outer = u_wind        # wind velocity beyond R_crossover
    else:
        u_outer = u_ejecta_B    # continue homologous beyond R_crossover (no wind material)
    u = np.where(mask_blend_AB, u_blend_AB,
                 np.where(mask_shell, u_ejecta,
                          np.where(mask_ejecta_B_vel, u_ejecta_B, u_outer)))
    fx_H  = np.where(mask_blend_BW, fx_H_blend_BW,
                     np.where(mask_blend_AB, fx_H_blend_AB,
                              np.where(mask_shell, fx_H_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_H_ejecta_B, fx_H_wind_arr))))
    fx_He = np.where(mask_blend_BW, fx_He_blend_BW,
                     np.where(mask_blend_AB, fx_He_blend_AB,
                              np.where(mask_shell, fx_He_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_He_ejecta_B, fx_He_wind_arr))))
    fx_O  = np.where(mask_blend_BW, fx_O_blend_BW,
                     np.where(mask_blend_AB, fx_O_blend_AB,
                              np.where(mask_shell, fx_O_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_O_ejecta_B, fx_O_wind_arr))))
    fx_O  = np.where(mask_blend_BW, fx_O_blend_BW,
                     np.where(mask_blend_AB, fx_O_blend_AB,
                              np.where(mask_shell, fx_O_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_O_ejecta_B, fx_O_wind_arr))))
    fx_Si = np.where(mask_blend_BW, fx_Si_blend_BW,
                     np.where(mask_blend_AB, fx_Si_blend_AB,
                              np.where(mask_shell, fx_Si_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_Si_ejecta_B, fx_Si_wind_arr))))
    fx_Ca = np.where(mask_blend_BW, fx_Ca_blend_BW,
                     np.where(mask_blend_AB, fx_Ca_blend_AB,
                              np.where(mask_shell, fx_Ca_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_Ca_ejecta_B, fx_Ca_wind_arr))))
    fx_Fe = np.where(mask_blend_BW, fx_Fe_blend_BW,
                     np.where(mask_blend_AB, fx_Fe_blend_AB,
                              np.where(mask_shell, fx_Fe_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_Fe_ejecta_B, fx_Fe_wind_arr))))
    fx_Ni = np.where(mask_blend_BW, fx_Ni_blend_BW,
                     np.where(mask_blend_AB, fx_Ni_blend_AB,
                              np.where(mask_shell, fx_Ni_ejecta_arr,
                                       np.where(mask_ejecta_B, fx_Ni_ejecta_B, fx_Ni_wind_arr))))

    # ================================================================
    # POST-PROCESSING: Gaussian smoothing of velocity at transitions
    # ================================================================
    # The sigmoid blending reduces discontinuities but can still leave
    # sharp gradients that blow up sigma when grid refinement inserts
    # points. Apply a gentle boxcar convolution over ~5-11 points at
    # the A/B and B/Wind transition boundaries.
    def boxcar_smooth_1d(x, n_window=7):
        """Apply centered moving-average (boxcar) smoothing."""
        if len(x) < n_window:
            return x.copy()
        kernel = np.ones(n_window) / n_window
        pad = n_window // 2
        x_padded = np.pad(x, (pad, pad), mode='edge')
        return np.convolve(x_padded, kernel, mode='valid')

    # Identify A/B and B/Wind transition zones by rapid change in dlnv/dlnr
    # Homologous: slope ≈ 1, Stationary: slope ≈ 0. Transitions have intermediate
    # slopes that change rapidly between adjacent points.
    u_log = np.log(np.maximum(u, 1e-30))
    r_log = np.log(np.maximum(r, 1e-30))
    # First derivative (slope) on a coarser stencil to avoid point-to-point noise
    slope = np.zeros_like(u)
    for i in range(2, len(u)-2):
        dlnu = u_log[i+2] - u_log[i-2]
        dlnr = r_log[i+2] - r_log[i-2]
        if abs(dlnr) > 1e-30:
            slope[i] = dlnu / dlnr
    # Second derivative (curvature) of slope
    curvature = np.zeros_like(slope)
    curvature[1:-1] = slope[2:] - 2*slope[1:-1] + slope[:-2]
    # Transition: curvature is large AND slope is between 0.1 and 0.9
    # (avoids catching pure homologous or pure stationary regions)
    transition_mask = (np.abs(curvature) > 0.1) & (slope > 0.1) & (slope < 0.9)
    # Expand mask by neighbours
    n_expand = 7
    for _ in range(n_expand):
        transition_mask[1:] |= transition_mask[:-1]
        transition_mask[:-1] |= transition_mask[1:]
    n_trans = np.sum(transition_mask)
    if n_trans > 0:
        print(f"  [Smoothing] Applying velocity smoothing at {n_trans} transition points")
        u_smoothed = boxcar_smooth_1d(u, n_window=7)
        u_orig = u.copy()
        u = np.where(transition_mask, u_smoothed, u)
        dv_max = np.max(np.abs(u - u_orig))
        print(f"    Max velocity change from smoothing: {dv_max:.2e} cm/s ({dv_max/1e5:.2f} km/s)")

    # Calculate sigma (velocity dispersion) = dlnV/dlnr - 1 on final blended u
    u_safe = np.maximum(u, 1e-10)
    sigma = np.zeros_like(r)
    for i in range(1, len(r)-1):
        dlnv = np.log(u_safe[i+1]) - np.log(u_safe[i-1])
        dlnr = np.log(r[i+1]) - np.log(r[i-1])
        sigma[i] = dlnv / dlnr - 1
    sigma[0] = sigma[1]
    sigma[-1] = sigma[-2]
    # ================================================================
    # TEMPERATURE
    # ================================================================
    T_min = 1000.0  # K minimum floor
    T_interface = 2000  # K at shell/wind boundary
    if temp_profile_type == 'adiabatic':
        GAMMA = 5.0 / 3.0
        # Use the last point inside the shell for rho_at_inner to avoid picking up
        # the tenuous wind density, which would give unphysically high shell temps
        idx_shell = np.where(r < R_inner_shell)[0]
        if len(idx_shell) > 0:
            idx_inner = idx_shell[-1]
        else:
            idx_inner = np.argmin(np.abs(r - R_inner_shell))
        rho_at_inner = rho[idx_inner]
        T_shell = T_interface * (rho / rho_at_inner)**(GAMMA - 1.0)
        # Wind and ejecta B: constant T_CSM (5000 K for metal-rich ejecta)
        T_csm = np.full_like(r, T_CSM)
        # Blend temperature across A/B boundary with same sigmoid as density/velocity
        # to eliminate the contact discontinuity that crashes STELLA's integrator
        T = T_shell * (1.0 - blend_alpha_AB) + T_csm * blend_alpha_AB
    else:
        T = np.full_like(r, T_CSM)
    T = np.maximum(T, T_min)
    # Apply density floor to prevent numerical singularities
    rho_floor = 1e-30  # g/cm^3
    rho = np.maximum(rho, rho_floor)
    # Calculate pressure from temperature and density using ideal gas law
    # Calculate mean molecular weight from actual composition
    # For neutral gas: mu = 1 / (sum(fx_i / A_i)) where A_i is atomic mass
    # H: 1, He: 4, O: 16, Si: 28, Ca: 40, Fe: 56
    mu = 1.0 / (fx_H/1.0 + fx_He/4.0 + fx_O/16.0 + fx_Si/28.0 + fx_Ca/40.0 + fx_Fe/56.0)
    k_B = 1.381e-16  # erg/K
    m_H = 1.673e-24  # g
    p = rho * k_B * T / (mu * m_H)
    # No normalization needed - we set composition explicitly
    return {
        'r': r,
        'rho': rho,
        'p': p,
        'u': u,
        'T': T,
        'sigma': sigma,
        'fx_H': fx_H,
        'fx_He': fx_He,
        'fx_O': fx_O,
        'fx_Si': fx_Si,
        'fx_Ca': fx_Ca,
        'fx_Fe': fx_Fe,
        'fx_Ni': fx_Ni,
        'parameters': {
            'M_ejecta_A_Msun': M_ejecta_A_g / Msun,
            'M_ejecta_B_Msun': M_ejecta_B_g / Msun,
            'Mdot_wind_Msun_yr': Mdot_wind,
            'R_inner_shell_cm': R_inner_shell,
            'R_crossover_cm': R_crossover,
            'R_CSM_max_cm': R_CSM_max,
            'Mdot_CSM_in_Msun_yr': Mdot_CSM_in * 365.25 * day / Msun,
            'days_ejecta_B': days_ejecta_B,
            'n_cells': n_cells,
            'use_rsg_mode': use_rsg_mode,
            'R_star_cm': R_star_cm if use_rsg_mode else None,
            'V_inf_kms': V_inf_kms if use_rsg_mode else None,
            'beta_wind': beta_wind if use_rsg_mode else None,
            'n_shell': n_shell if use_rsg_mode else None,
        }
    }
def write_heracles_profile(profile, output_file):
    """
    Write profile to HERACLES format.
    Parameters:
    -----------
    profile : dict
        Dictionary containing r, rho, p, u, and mass fractions
    output_file : str or Path
        Output file path
    """
    output_file = Path(output_file)
    with open(output_file, 'w') as f:
        # Write header
        f.write("# HERACLES profileinit file for embedded H/He-rich supernova model\n")
        f.write("# Based on arXiv:1503.05463 approach\n")
        f.write("# columns: r_cm rho_cgs p_cgs u_cgs fx_H fx_He fx_O fx_Si fx_Ca fx_Fe fx_Ni\n")
        f.write(f"# Model parameters:\n")
        params = profile['parameters']
        f.write(f"#   M_ejecta_A = {params['M_ejecta_A_Msun']:.2f} Msun\n")
        f.write(f"#   M_ejecta_B = {params['M_ejecta_B_Msun']:.2f} Msun\n")
        f.write(f"#   Mdot_wind = {params['Mdot_wind_Msun_yr']:.3e} Msun/yr\n")
        f.write(f"#   R_inner_shell = {params['R_inner_shell_cm']:.2e} cm\n")
        f.write(f"#   R_crossover = {params['R_crossover_cm']:.2e} cm (where wind > ejecta B)\n")
        f.write(f"#   R_CSM_max = {params['R_CSM_max_cm']:.2e} cm\n")
        f.write(f"#   Mdot_CSM_in = {params['Mdot_CSM_in_Msun_yr']:.3e} Msun/yr\n")
        if params.get('use_rsg_mode'):
            f.write(f"#   RSG mode: R_star={params['R_star_cm']:.2e} cm, V_inf={params['V_inf_kms']} km/s, beta={params['beta_wind']}, n_shell={params['n_shell']}\n")
        f.write(f"#   Ejecta B age = {params['days_ejecta_B']:.1f} days\n")
        f.write(f"#   n_cells = {params['n_cells']}\n")
        f.write("# Row order: strictly increasing r_cm (required by HERACLES init.f90 interpolation).\n")
        f.write("#\n")
        # Write data (HERACLES format: r rho p u fx_H fx_He fx_O fx_Si fx_Ca fx_Fe fx_Ni)
        n = len(profile['r'])
        for i in range(n):
            f.write(f"{profile['r'][i]:.8e} {profile['rho'][i]:.8e} "
                   f"{profile['p'][i]:.8e} {profile['u'][i]:.8e} "
                   f"{profile['fx_H'][i]:.8e} {profile['fx_He'][i]:.8e} "
                   f"{profile['fx_O'][i]:.8e} "
                   f"{profile['fx_Si'][i]:.8e} {profile['fx_Ca'][i]:.8e} "
                   f"{profile['fx_Fe'][i]:.8e} {profile['fx_Ni'][i]:.8e}\n")
    print(f"Profile written to {output_file}")
    print(f"Total cells: {n}")
    print(f"Radius range: {profile['r'][0]:.2e} to {profile['r'][-1]:.2e} cm")
    print(f"Density range: {profile['rho'].min():.2e} to {profile['rho'].max():.2e} g/cm^3")
    print(f"Velocity range: {profile['u'].min():.2e} to {profile['u'].max():.2e} cm/s")
    print(f"Temperature range: {profile['T'].min():.2e} to {profile['T'].max():.2e} K")
    # Verify kinetic energy
    dr = np.diff(profile['r'], prepend=profile['r'][0])
    dm = 4 * np.pi * profile['r']**2 * profile['rho'] * dr
    E_kin_actual = 0.5 * np.sum(dm * profile['u']**2)
    print(f"Kinetic energy: {E_kin_actual:.2e} erg (target: {1e51:.2e} erg)")
def plot_profile(profile, output_dir='plots'):
    """
    Create plots of the profile similar to arXiv:1503.05463 paper style.
    Parameters:
    -----------
    profile : dict
        Dictionary containing r, rho, p, u, and mass fractions
    output_dir : str or Path
        Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    r = profile['r']
    rho = profile['rho']
    p = profile['p']
    u = profile['u']
    T = profile['T']
    sigma = profile['sigma']
    fx_H = profile['fx_H']
    fx_He = profile['fx_He']
    fx_O = profile['fx_O']
    fx_Si = profile['fx_Si']
    fx_Ca = profile['fx_Ca']
    fx_Fe = profile['fx_Fe']
    fx_Ni = profile['fx_Ni']
    params = profile['parameters']
    # Convert radius to units of 10^15 cm for most plots, and Rsun for structure overview
    R15 = 1e15  # cm per unit
    Rsun = 6.957e10  # cm per solar radius
    r_plt = r / R15  # radius in units of 10^15 cm
    r_plt_rsun = r / Rsun  # radius in units of solar radii
    use_rsg_mode = params.get('use_rsg_mode', False)
    R_inner_shell_cm = params.get('R_inner_shell_cm', params.get('R_t_cm', 1e15))
    r_t_plt = R_inner_shell_cm / R15
    r_t_plt_rsun = R_inner_shell_cm / Rsun
    # NEW LAYOUT: use R_crossover instead of R_CSM_transition
    r_crossover_plt = params['R_crossover_cm'] / R15
    r_crossover_plt_rsun = params['R_crossover_cm'] / Rsun
    r_CSM_max_plt = params['R_CSM_max_cm'] / R15
    r_CSM_max_plt_rsun = params['R_CSM_max_cm'] / Rsun
    xlabel_r = r'Radius ($10^{15}$ cm)'
    xlabel_r_rsun = r'Radius ($R_{\odot}$)'
    # Mode-aware labels (NEW LAYOUT: A → B → Wind)
    inner_boundary_label = r'RSG surface ($R_\star$)' if use_rsg_mode else 'Ejecta A/B boundary'
    inner_region_label   = 'RSG Shell' if use_rsg_mode else 'Ejecta A'
    middle_region_label  = 'Ejecta B'
    outer_region_label   = 'Wind'
    if use_rsg_mode:
        title_suffix = f'RSG: R_star={params["R_star_cm"]:.1e} cm, V_inf={params["V_inf_kms"]} km/s'
    else:
        title_suffix = f'A: {params.get("days_ejecta_A", 0):.0f}d, B: {params["days_ejecta_B"]:.0f}d'
    # Dynamic x-axis limits based on actual data range with padding
    r_min_plt = r_plt.min()
    r_max_plt = r_plt.max()
    xlim_min = r_min_plt * 0.5  # 50% padding on lower end
    xlim_max = r_max_plt * 2.0  # 100% padding on upper end
    # Log x-axis limits for Rsun units — use actual inner radius
    xlim_min_rsun = r_plt_rsun.min() * 0.5
    xlim_max_rsun = r_plt_rsun.max() * 2.0
    # Set up plot style similar to paper
    plt.style.use('seaborn-v0_8-darkgrid')
    # Unified region label helper — always uses axes-fraction x positions
    # computed from log(r) midpoints of each region, so labels are centred
    # correctly regardless of the x-axis scale (log or linear) or unit choice.
    def _label_positions_logr(r_lo, r_b1, r_b2, r_hi):
        """Return three axes-fraction x positions at log-midpoints of three regions.
        r_lo/r_b1/r_b2/r_hi are boundary radii in the SAME units as the plot x-axis.
        Uses log midpoints; clamps with minimum separation of 0.07 axes units."""
        eps = 1e-30
        x1 = (np.log10(max(r_lo, eps)) + np.log10(max(r_b1, eps))) / 2
        x2 = (np.log10(max(r_b1, eps)) + np.log10(max(r_b2, eps))) / 2
        x3 = (np.log10(max(r_b2, eps)) + np.log10(max(r_hi, eps))) / 2
        lo = np.log10(max(r_lo, eps))
        hi = np.log10(max(r_hi, eps))
        span = hi - lo if hi > lo else 1.0
        ax1f = np.clip((x1 - lo) / span, 0.02, 0.93)
        ax2f = np.clip((x2 - lo) / span, ax1f + 0.07, 0.93)
        ax3f = np.clip((x3 - lo) / span, ax2f + 0.07, 0.98)
        return ax1f, ax2f, ax3f
    def add_region_labels(ax):
        """Labels for log-x plots in 10^15 cm units (NEW LAYOUT: A → B → Wind).
        Uses the actual plotted xlim boundaries so positions are always correct."""
        a1, a2, a3 = _label_positions_logr(xlim_min, r_t_plt,
                                            r_crossover_plt, xlim_max)
        for xf, txt, col in [(a1, inner_region_label, 'steelblue'),
                              (a2, middle_region_label, 'green'),
                              (a3, outer_region_label,  'firebrick')]:
            ax.text(xf, 0.97, txt, fontsize=9, fontweight='bold', color=col,
                    ha='center', va='top', transform=ax.transAxes)
    def add_region_labels_rsun(ax, x_min, x_max):
        """Labels for log-x plots in R_sun units (NEW LAYOUT: A → B → Wind)."""
        a1, a2, a3 = _label_positions_logr(x_min, r_t_plt_rsun,
                                            r_crossover_plt_rsun, x_max)
        for xf, txt, col in [(a1, inner_region_label, 'steelblue'),
                              (a2, middle_region_label, 'green'),
                              (a3, outer_region_label,  'firebrick')]:
            ax.text(xf, 0.95, txt, fontsize=9, fontweight='bold', color=col,
                    ha='center', va='top', transform=ax.transAxes)
    # Figure 1: Density and velocity profiles
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    # Density plot
    ax1.loglog(r_plt, rho, 'b-', linewidth=2, label='Density')
    ax1.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.5,
                label=inner_boundary_label)
    ax1.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.5,
                label='Wind onset')
    ax1.set_xlabel(xlabel_r, fontsize=12)
    ax1.set_ylabel('Density (g cm$^{-3}$)', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(xlim_min, xlim_max)
    ax1.set_ylim(rho.min() * 0.1, rho.max() * 10)
    add_region_labels(ax1)
    # Velocity plot (km/s)
    ax2.loglog(r_plt, u / 1e5, 'r-', linewidth=2, label='Velocity')
    ax2.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.5,
                label=inner_boundary_label)
    ax2.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.5,
                label='Wind onset')
    ax2.set_xlabel(xlabel_r, fontsize=12)
    ax2.set_ylabel('Velocity (km s$^{-1}$)', fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(xlim_min, xlim_max)
    # Handle zero velocity for log scale: use minimum positive value
    u_km_s = u / 1e5
    u_min_positive = u_km_s[u_km_s > 0].min() if np.any(u_km_s > 0) else 1e-3
    ax2.set_ylim(u_min_positive * 0.1, u_km_s.max() * 10)
    add_region_labels(ax2)
    plt.tight_layout()
    plt.savefig(output_dir / 'density_velocity_profiles.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved density/velocity plot to {output_dir / 'density_velocity_profiles.png'}")
    # Figure 2: Composition profiles
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.semilogx(r_plt, fx_H, label='H',  linewidth=2)
    ax.semilogx(r_plt, fx_He, label='He', linewidth=2)
    ax.semilogx(r_plt, fx_O, label='O', linewidth=2)
    ax.semilogx(r_plt, fx_Si, label='Si', linewidth=2)
    ax.semilogx(r_plt, fx_Ca, label='Ca', linewidth=2)
    ax.semilogx(r_plt, fx_Fe, label='Fe', linewidth=2)
    ax.semilogx(r_plt, fx_Ni, label='Ni', linewidth=2)
    ax.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.5,
               label=inner_boundary_label)
    ax.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.5,
               label='Wind onset')
    ax.set_xlabel(xlabel_r, fontsize=12)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.legend(fontsize=10, loc='center left')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.1)
    ax.set_xlim(xlim_min, xlim_max)
    add_region_labels(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'composition_profile.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved composition plot to {output_dir / 'composition_profile.png'}")
    # Compute enclosed mass for reuse in multiple plots
    dr = np.diff(r, prepend=r[0])
    dm = 4 * np.pi * r**2 * rho * dr
    enclosed_mass_msun = np.cumsum(dm) / 1.989e33  # Msun
    # Figure 3: Combined structure overview (2x3 grid) - Rsun units, log x-axis
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    ax1, ax2, ax3 = axes[0]
    ax4, ax5, ax6 = axes[1]
    def _vlines_rsun(ax):
        ax.axvline(r_t_plt_rsun, color='steelblue', linestyle='--', alpha=0.5)
        ax.axvline(r_crossover_plt_rsun, color='green', linestyle='--', alpha=0.5)
    # Density (log-log)
    ax1.loglog(r_plt_rsun, rho, 'b-', linewidth=2)
    _vlines_rsun(ax1)
    ax1.set_xlabel(xlabel_r_rsun)
    ax1.set_ylabel('Density (g cm$^{-3}$)')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(xlim_min_rsun, xlim_max_rsun)
    ax1.set_ylim(rho.min() * 0.1, rho.max() * 10)
    add_region_labels_rsun(ax1, xlim_min_rsun, xlim_max_rsun)
    # Velocity (log x, log y — handles zero-vel inner shell)
    u_km = u / 1e5
    u_pos = np.where(u_km > 0, u_km, np.nan)
    ax2.loglog(r_plt_rsun, u_pos, 'r-', linewidth=2)
    _vlines_rsun(ax2)
    ax2.set_xlabel(xlabel_r_rsun)
    ax2.set_ylabel('Velocity (km s$^{-1}$)')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(xlim_min_rsun, xlim_max_rsun)
    add_region_labels_rsun(ax2, xlim_min_rsun, xlim_max_rsun)
    # Pressure (log-log)
    ax3.loglog(r_plt_rsun, p, 'g-', linewidth=2)
    _vlines_rsun(ax3)
    ax3.set_xlabel(xlabel_r_rsun)
    ax3.set_ylabel('Pressure (erg cm$^{-3}$)')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(xlim_min_rsun, xlim_max_rsun)
    ax3.set_ylim(p.min() * 0.1, p.max() * 10)
    add_region_labels_rsun(ax3, xlim_min_rsun, xlim_max_rsun)
    # Composition (log x)
    ax4.semilogx(r_plt_rsun, fx_H, label='H',  linewidth=2)
    ax4.semilogx(r_plt_rsun, fx_He, label='He',  linewidth=2)
    ax4.semilogx(r_plt_rsun, fx_O, label='O',  linewidth=2)
    ax4.semilogx(r_plt_rsun, fx_Si, label='Si',  linewidth=2)
    ax4.semilogx(r_plt_rsun, fx_Ca, label='Ca',  linewidth=2)
    ax4.semilogx(r_plt_rsun, fx_Fe, label='Fe',  linewidth=2)
    ax4.semilogx(r_plt_rsun, fx_Ni, label='Ni',  linewidth=2)
    _vlines_rsun(ax4)
    ax4.set_xlabel(xlabel_r_rsun)
    ax4.set_ylabel('Mass Fraction')
    ax4.legend(fontsize=8, ncol=2, loc='center right')
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, 1.1)
    ax4.set_xlim(xlim_min_rsun, xlim_max_rsun)
    add_region_labels_rsun(ax4, xlim_min_rsun, xlim_max_rsun)
    # Enclosed mass vs radius (log x)
    ax5.semilogx(r_plt_rsun, enclosed_mass_msun, 'k-', linewidth=2)
    _vlines_rsun(ax5)
    ax5.axhline(params['M_ejecta_A_Msun'], color='k', linestyle=':', alpha=0.5,
                label=f'M_ejecta_A = {params["M_ejecta_A_Msun"]:.1f} M$_\\odot$')
    ax5.set_xlabel(xlabel_r_rsun)
    ax5.set_ylabel('Enclosed Mass (Msun)')
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim(xlim_min_rsun, xlim_max_rsun)
    ax5.set_ylim(enclosed_mass_msun.min() * 0.9, enclosed_mass_msun.max() * 1.1)
    add_region_labels_rsun(ax5, xlim_min_rsun, xlim_max_rsun)
    # Temperature (log-log)
    ax6.loglog(r_plt_rsun, T, 'g-', linewidth=2)
    _vlines_rsun(ax6)
    ax6.set_xlabel(xlabel_r_rsun)
    ax6.set_ylabel('Temperature (K)')
    ax6.grid(True, alpha=0.3)
    ax6.set_xlim(xlim_min_rsun, xlim_max_rsun)
    ax6.set_ylim(T.min() * 0.1, T.max() * 10)
    add_region_labels_rsun(ax6, xlim_min_rsun, xlim_max_rsun)
    plt.suptitle(f'Progenitor Structure — {title_suffix}',
                 fontsize=15, fontweight='bold')
    # Single shared legend for the boundary vlines
    legend_handles = [
        Line2D([0], [0], color='steelblue', linestyle='--', alpha=0.7, label=inner_boundary_label),
        Line2D([0], [0], color='firebrick', linestyle='--', alpha=0.7, label='Ejecta B onset'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=2, fontsize=10,
               bbox_to_anchor=(0.5, 0.0), framealpha=0.8)
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(output_dir / 'structure_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved structure overview to {output_dir / 'structure_overview.png'}")
    # Figure 4: Enclosed mass vs radius (dedicated, publication-quality)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogx(r_plt, enclosed_mass_msun, 'k-', linewidth=2, label='Enclosed mass')
    ax.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.8,
               label=inner_boundary_label)
    ax.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.8,
               label='Wind onset')
    idx_t = np.argmin(np.abs(r - R_inner_shell_cm))
    idx_crossover = np.argmin(np.abs(r - params['R_crossover_cm']))
    M_shell_plot  = enclosed_mass_msun[idx_t]
    M_ejecta_B_plot = enclosed_mass_msun[idx_crossover] - enclosed_mass_msun[idx_t]
    M_wind_plot   = enclosed_mass_msun[-1] - enclosed_mass_msun[idx_crossover]
    ax.axhline(M_shell_plot, color='steelblue', linestyle=':', alpha=0.6,
               label=f'M_{inner_region_label.replace(" ","_")} = {M_shell_plot:.3f} M$_\\odot$')
    ax.axhline(M_shell_plot + M_ejecta_B_plot, color='green', linestyle=':', alpha=0.6,
               label=f'M_ejecta_B = {M_ejecta_B_plot:.3f} M$_\\odot$')
    ax.axhline(enclosed_mass_msun[-1], color='firebrick', linestyle=':', alpha=0.6,
               label=f'M_wind = {M_wind_plot:.3f} M$_\\odot$')
    ax.set_xlabel(xlabel_r, fontsize=12)
    ax.set_ylabel('Enclosed Mass (Msun)', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(enclosed_mass_msun.min() * 0.9, enclosed_mass_msun.max() * 1.1)
    # Region labels (line 1) + mass labels (line 2) at top of axes (NEW LAYOUT: A → B → Wind)
    a1, a2, a3 = _label_positions_logr(xlim_min, r_t_plt, r_crossover_plt, xlim_max)
    for xf, lbl, mass, col in [
            (a1, inner_region_label, M_shell_plot,   'steelblue'),
            (a2, middle_region_label, M_ejecta_B_plot, 'green'),
            (a3, outer_region_label,  M_wind_plot,     'firebrick')]:
        ax.text(xf, 0.97, lbl,
                fontsize=9, fontweight='bold', color=col,
                ha='center', va='top', transform=ax.transAxes)
        ax.text(xf, 0.91, f'{mass:.2f} M$_\\odot$',
                fontsize=8, color=col,
                ha='center', va='top', transform=ax.transAxes)
    plt.tight_layout()
    plt.savefig(output_dir / 'enclosed_mass.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved enclosed mass plot to {output_dir / 'enclosed_mass.png'}")
    # Figure 5: Optical depth vs radius (outer to inner boundary)
    # Calculate optical depth from outer boundary inward
    kappa = 0.3  # cm²/g (typical opacity for SN ejecta)
    dr = np.diff(r, prepend=r[0])
    # Optical depth element
    dtau = kappa * rho * dr
    # Cumulative sum from outer boundary (reverse order)
    tau = np.cumsum(dtau[::-1])[::-1]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.loglog(r_plt, tau, 'k-', linewidth=2, label='Optical depth')
    ax.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.8,
               label=inner_boundary_label)
    ax.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.8,
               label='Wind onset')
    ax.axhline(1.0, color='g', linestyle=':', alpha=0.6,
               label='τ = 1 (photosphere)')
    ax.axhline(2/3, color='orange', linestyle=':', alpha=0.6,
               label='τ = 2/3')
    ax.set_xlabel(xlabel_r, fontsize=12)
    ax.set_ylabel('Optical Depth τ', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(tau.min() * 0.1, tau.max() * 10)
    add_region_labels(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'optical_depth.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved optical depth plot to {output_dir / 'optical_depth.png'}")
    # Figure 6: Sigma (velocity dispersion) vs radius (CMFGEN microturbulence)
    fig, ax = plt.subplots(figsize=(10, 6))
    # Plot sigma values
    ax.semilogx(r_plt, sigma, 'm-', linewidth=2, label='Sigma (dlnV/dlnr-1)')
    ax.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.8,
               label=inner_boundary_label)
    ax.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.8,
               label='Wind onset')
    ax.axhline(y=-0.99, color='k', linestyle=':', alpha=0.5, label='σ=-0.99 (clipping limit)')
    ax.axhline(y=-1.0, color='gray', linestyle='--', alpha=0.5, label='σ=-1.0 (CMFGEN limit)')
    ax.set_xlabel(xlabel_r, fontsize=12)
    ax.set_ylabel('Sigma (dimensionless)', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(xlim_min, xlim_max)
    # Use symlog y-axis to show strong negative values (same as convert_heracles_to_cmfgen.py)
    ax.set_yscale('symlog', linthresh=1.0)
    # Set reasonable y-axis limits
    sigma_valid = sigma[np.isfinite(sigma)]
    if len(sigma_valid) > 0:
        sigma_min = sigma_valid.min()
        sigma_max = sigma_valid.max()
        ax.set_ylim(sigma_min * 1.1, sigma_max * 1.1)
    else:
        ax.set_ylim(-10, 10)
    add_region_labels(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'sigma.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved sigma plot to {output_dir / 'sigma.png'}")
    # Figure 7: dr and dr/r vs radius — zone resolution diagnostic
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    dr_arr = np.diff(r)
    r_mid  = 0.5 * (r[:-1] + r[1:])
    r_mid_plt = r_mid / R15
    ax1.loglog(r_mid_plt, dr_arr, 'b-', linewidth=1.5, label='dr (cm)')
    ax1.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.5, label=inner_boundary_label)
    ax1.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.5, label='Wind onset')
    ax1.set_ylabel('Zone width dr (cm)', fontsize=12)
    ax1.legend(fontsize=9, loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(xlim_min, xlim_max)
    add_region_labels(ax1)
    dr_over_r = dr_arr / r_mid
    ax2.loglog(r_mid_plt, dr_over_r, 'g-', linewidth=1.5, label='dr/r')
    ax2.axvline(r_t_plt, color='steelblue', linestyle='--', alpha=0.5)
    ax2.axvline(r_crossover_plt, color='green', linestyle='--', alpha=0.5)
    ax2.axhline(dr_over_r.mean(), color='k', linestyle=':', alpha=0.6,
                label=f'mean dr/r = {dr_over_r.mean():.3f}')
    ax2.set_xlabel(xlabel_r, fontsize=12)
    ax2.set_ylabel('dr / r  (fractional zone width)', fontsize=12)
    ax2.legend(fontsize=9, loc='upper left')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(xlim_min, xlim_max)
    add_region_labels(ax2)
    plt.suptitle(f'Zone Resolution ({title_suffix})', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / 'zone_resolution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved zone resolution plot to {output_dir / 'zone_resolution.png'}")
    # Figure 8: Stacked mass fractions vs radius (log Rsun)
    fig, ax = plt.subplots(figsize=(10, 8))
    # Stack the mass fractions from bottom to top
    ax.stackplot(r_plt_rsun, fx_H, fx_He, fx_O, fx_Si, fx_Ca, fx_Fe, fx_Ni,
                 labels=['H', 'He', 'O', 'Si', 'Ca', 'Fe', 'Ni'],
                 colors=['blue', 'orange', 'red', 'purple', 'brown', 'gray', 'pink'],
                 alpha=0.7)
    ax.axvline(r_t_plt_rsun, color='steelblue', linestyle='--', alpha=0.8,
               label=inner_boundary_label)
    ax.axvline(r_crossover_plt_rsun, color='green', linestyle='--', alpha=0.8,
               label='Wind onset')
    ax.set_xlabel(xlabel_r_rsun, fontsize=12)
    ax.set_ylabel('Mass Fraction', fontsize=12)
    ax.legend(fontsize=10, loc='lower center', bbox_to_anchor=(0.5, 1.05), ncol=4)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(xlim_min_rsun, xlim_max_rsun)
    ax.set_ylim(0, 1.1)
    ax.set_xscale('log')
    plt.tight_layout()
    plt.savefig(output_dir / 'stacked_composition.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved stacked composition plot to {output_dir / 'stacked_composition.png'}")
def main():
    parser = argparse.ArgumentParser(
        description='Create HERACLES profile for embedded H/He system with wind and stripped envelope'
    )
    parser.add_argument('--days_ejecta_A', type=float, default=11.6,
                        help='H/He-rich ejecta A age in days (default: 11.6 for model X)')
    parser.add_argument('--days_ejecta_B', type=float, default=400,
                        help='Metal-rich ejecta B initial evolution in days (default: 400)')
    parser.add_argument('--E_kin_ejecta_A', type=float, default=1e51,
                        help='Kinetic energy of ejecta A in erg (default: 1e51)')
    parser.add_argument('--E_kin_ejecta_B', type=float, default=0,
                        help='Kinetic energy of ejecta B in erg (default: 0, homologous expansion only)')
    parser.add_argument('--n_cells', type=int, default=10000,
                        help='Number of radial cells (default: 10000)')
    parser.add_argument('--output', type=str, 
                        default='profile_embedded_hhe.dat',
                        help='Output file name (default: profile_embedded_hhe.dat)')
    parser.add_argument('--plots', action='store_true',
                        help='Generate diagnostic plots in plots/ directory')
    parser.add_argument('--plot_dir', type=str, default='plots',
                        help='Directory for plots (default: plots)')
    parser.add_argument('--fx_Fe', type=float, default=0.03,
                        help='Iron mass fraction in outer ejecta (default: 0.03, SN2011dh-like)')
    parser.add_argument('--fx_Si', type=float, default=0.05,
                        help='Silicon mass fraction in outer ejecta (default: 0.05, SN2011dh-like)')
    parser.add_argument('--fx_O', type=float, default=0.15,
                        help='Oxygen mass fraction in outer ejecta (default: 0.15, SN2011dh-like)')
    parser.add_argument('--fx_Ca', type=float, default=0.02,
                        help='Calcium mass fraction in outer ejecta (default: 0.02, SN2011dh-like)')
    parser.add_argument('--fx_He', type=float, default=0.40,
                        help='Helium mass fraction in outer ejecta (default: 0.40, SN2011dh-like Type IIb)')
    parser.add_argument('--fx_H', type=float, default=0.01,
                        help='Hydrogen mass fraction in outer ejecta (default: 0.01, SN2011dh-like Type IIb)')
    parser.add_argument('--M_ejecta_A', type=float, default=9.8,
                        help='Mass of ejecta A in Msun (default: 9.8)')
    parser.add_argument('--Mdot_wind', type=float, default=1e-3,
                        help='Mass loss rate for wind in Msun/yr (default: 1e-3)')
    parser.add_argument('--M_ejecta_B', type=float, default=0.1,
                        help='Mass of ejecta B in Msun (default: 0.1)')
    parser.add_argument('--fx_H_A', type=float, default=0.70,
                        help='Hydrogen mass fraction in ejecta A (default: 0.70)')
    parser.add_argument('--fx_He_A', type=float, default=0.28,
                        help='Helium mass fraction in ejecta A (default: 0.28)')
    parser.add_argument('--fx_H_wind', type=float, default=0.70,
                        help='Hydrogen mass fraction in wind (default: 0.70)')
    parser.add_argument('--fx_He_wind', type=float, default=0.28,
                        help='Helium mass fraction in wind (default: 0.28)')
    parser.add_argument('--fx_H_B', type=float, default=0.01,
                        help='Hydrogen mass fraction in ejecta B (default: 0.01)')
    parser.add_argument('--fx_He_B', type=float, default=0.40,
                        help='Helium mass fraction in ejecta B (default: 0.40)')
    parser.add_argument('--fx_O_B', type=float, default=0.15,
                        help='Oxygen mass fraction in ejecta B (default: 0.15)')
    parser.add_argument('--fx_Ca_B', type=float, default=0.02,
                        help='Calcium mass fraction in ejecta B (default: 0.02)')
    args = parser.parse_args()
    # Create profile
    print(f"Creating embedded H/He profile...")
    print(f"  Ejecta A age: {args.days_ejecta_A} days")
    print(f"  Ejecta B initial evolution: {args.days_ejecta_B} days")
    profile = create_embedded_hhe_profile(days_ejecta_A=args.days_ejecta_A,
                                          days_ejecta_B=args.days_ejecta_B,
                                          n_cells=args.n_cells,
                                          fx_Fe_ejecta=args.fx_Fe,
                                          fx_Si_ejecta=args.fx_Si,
                                          fx_He_ejecta=args.fx_He,
                                          fx_H_ejecta=args.fx_H,
                                          M_ejecta_A=args.M_ejecta_A,
                                          Mdot_wind=args.Mdot_wind,
                                          M_ejecta_B=args.M_ejecta_B,
                                          E_kin_ejecta_A=args.E_kin_ejecta_A,
                                          E_kin_ejecta_B=args.E_kin_ejecta_B,
                                          fx_H_A=args.fx_H_A,
                                          fx_He_A=args.fx_He_A,
                                          fx_O_A=0.009,
                                          fx_Si_A=0.0007,
                                          fx_Ca_A=0.0002,
                                          fx_Fe_A=0.0013,
                                          fx_H_wind=args.fx_H_wind,
                                          fx_He_wind=args.fx_He_wind,
                                          fx_H_B=args.fx_H_B,
                                          fx_He_B=args.fx_He_B,
                                          fx_O_B=args.fx_O_B,
                                          fx_Si_B=args.fx_Si,
                                          fx_Ca_B=args.fx_Ca_B,
                                          fx_Fe_B=args.fx_Fe,
                                          fx_Ni_B=0.02)
    # Write to file
    write_heracles_profile(profile, args.output)
    print("\nModel structure (paper model X):")
    print(f"  SN ejecta A (inner): {profile['parameters']['M_ejecta_A_Msun']:.2f} Msun")
    print(f"  SN ejecta B (outer): {profile['parameters']['M_ejecta_B_Msun']:.2f} Msun")
    print(f"  CSM (inner): Mdot = {profile['parameters']['Mdot_CSM_in_Msun_yr']:.3f} Msun/yr")
    print(f"  V_max = {profile['parameters']['V_max_km_s']:.1f} km/s")
    print(f"  Ejecta A age: {profile['parameters'].get('days_ejecta_A', 0):.1f} days")
    print(f"  Ejecta B initial evolution: {profile['parameters']['days_ejecta_B']:.1f} days")
    print(f"  v_CSM = {profile['parameters']['v_CSM_km_s']:.1f} km/s")
    # Calculate masses of each region
    r = profile['r']
    rho = profile['rho']
    dr = np.diff(r, prepend=r[0])
    dm = 4 * np.pi * r**2 * rho * dr
    enclosed_mass = np.cumsum(dm) / 1.989e33  # Msun
    R_t = profile['parameters'].get('R_t_cm', profile['parameters'].get('R_inner_shell_cm', 1e15))
    R_CSM_transition = profile['parameters']['R_CSM_transition_cm']
    R_CSM_max = profile['parameters']['R_CSM_max_cm']
    # Find indices for each region
    idx_t = np.argmin(np.abs(r - R_t))
    idx_transition = np.argmin(np.abs(r - R_CSM_transition))
    # Mass of inner ejecta (r < R_t)
    M_ejecta_inner = enclosed_mass[idx_t]
    # Mass of CSM inner (R_t < r < R_CSM_transition)
    M_CSM_inner = enclosed_mass[idx_transition] - enclosed_mass[idx_t]
    # Mass of CSM outer / ejecta B (r > R_CSM_transition)
    M_CSM_outer = enclosed_mass[-1] - enclosed_mass[idx_transition]
    print(f"\nRegion masses:")
    print(f"  Ejecta A (inner H/He, r < R_t): {M_ejecta_inner:.3f} Msun")
    print(f"  Wind (CSM inner, R_t < r < R_transition): {M_CSM_inner:.3f} Msun")
    print(f"  Ejecta B (outer metal-rich, r > R_transition): {M_CSM_outer:.3f} Msun")
    print(f"  Total mass: {enclosed_mass[-1]:.3f} Msun")
    print(f"Ejecta composition (SN2011dh-like):")
    print(f"  fx_H = {args.fx_H}")
    print(f"  fx_He = {args.fx_He}")
    print(f"  fx_C = {args.fx_C}")
    print(f"  Remaining = {1.0 - args.fx_Fe - args.fx_Si - args.fx_C - args.fx_Ca - args.fx_He - args.fx_H:.4f}")
    print(f"  fx_Si = {args.fx_Si}")
    print(f"  fx_Fe = {args.fx_Fe}")
    # Generate plots if requested
    if args.plots:
        print(f"\nGenerating plots in {args.plot_dir}/...")
        plot_profile(profile, output_dir=args.plot_dir)
        print("Plot generation complete.")
if __name__ == '__main__':
    main()
