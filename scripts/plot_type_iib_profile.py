#!/usr/bin/env python3
"""
Type IIb supernova ejecta profile in velocity space with homologous evolution.

Models the density and velocity structure of a Type IIb supernova at
50, 100, 200, and 400 days after explosion, assuming homologous expansion.

Core assumptions:
- Homologous expansion after the early post-explosion phase.
- Each mass shell moves at constant velocity: r(v,t) = v * t.
- Density scales with time: rho(v,t) = rho(v,t0) * (t0/t)^3.
- Profile is defined in velocity space, not fixed radius space.
- Temperature: radiation-dominated adiabatic, T proportional to t^-1.

Profile form (broken power law in velocity):
- Inner ejecta (v < v_t): shallow slope, rho proportional to v^(-m)
- Outer ejecta (v > v_t): steep slope, rho proportional to v^(-n)
- H-rich material confined to outer low-mass layers (v > v_H_cut)
- He and metal-rich material dominate the inner ejecta (v < v_H_cut)

Usage:
    python scripts/plot_type_iib_profile.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Physical constants
MSUN = 1.989e33      # g
DAY = 86400.0        # s
KM_S = 1e5           # cm/s
K_B = 1.381e-16      # erg/K
M_H = 1.673e-24      # g


def broken_powerlaw_density(v, rho_0, v_t, m, n):
    """
    Broken power-law density in velocity space.

    rho(v) = rho_0 * (v/v_t)^(-m)   for v <= v_t  (inner, shallow)
    rho(v) = rho_0 * (v/v_t)^(-n)   for v >  v_t  (outer, steep)

    Parameters
    ----------
    v : array
        Velocity in cm/s.
    rho_0 : float
        Density at v_t in g/cm^3 (at reference time t0).
    v_t : float
        Transition velocity in cm/s.
    m : float
        Inner slope exponent (shallow, e.g. 1-3).
    n : float
        Outer slope exponent (steep, e.g. 7-10).

    Returns
    -------
    rho : array
        Density at reference time in g/cm^3.
    """
    f = np.where(v <= v_t, (v / v_t)**(-m), (v / v_t)**(-n))
    return rho_0 * f


def mass_energy_integrals(v_grid, v_t, m, n):
    """
    Compute mass and energy integrals (in units of rho_0) for broken power law.

    M = 4*pi*t0^3 * rho_0 * I_m
    E = 0.5 * 4*pi*t0^5 * rho_0 * I_e

    where I_m = integral(f(v) * v^2 dv) and I_e = integral(f(v) * v^4 dv).
    """
    f = np.where(v_grid <= v_t, (v_grid / v_t)**(-m), (v_grid / v_t)**(-n))
    I_m = np.trapz(f * v_grid**2, v_grid)
    I_e = np.trapz(f * v_grid**4, v_grid)
    return I_m, I_e


def solve_v_t_for_energy(v_grid, m, n, E_kin, M_ej, t0, v_max):
    """
    Solve for v_t such that E_kin/M_ej = 0.5 * I_e / I_m.

    In homologous expansion:
        KE = 0.5 * 4*pi*t0^3 * rho_0 * I_e
        M  = 4*pi*t0^3 * rho_0 * I_m
    So E/M = 0.5 * I_e / I_m  (t0^3 and rho_0 cancel).

    Binary search over v_t. Returns (v_t, rho_0).
    """
    target = E_kin / M_ej  # specific kinetic energy (erg/g = cm^2/s^2)
    v_lo, v_hi = v_grid[0] * 1.1, v_max * 0.999
    for _ in range(100):
        v_mid = 0.5 * (v_lo + v_hi)
        I_m, I_e = mass_energy_integrals(v_grid, v_mid, m, n)
        ratio = 0.5 * I_e / I_m  # E/M (cm^2/s^2)
        if ratio > target:
            v_hi = v_mid
        else:
            v_lo = v_mid
    v_t = 0.5 * (v_lo + v_hi)
    I_m, _ = mass_energy_integrals(v_grid, v_t, m, n)
    rho_0 = M_ej / (4.0 * np.pi * t0**3 * I_m)
    return v_t, rho_0


def type_iib_profile(t_days, M_ej=2.0, E_kin=1.0e51,
                     v_min_kms=100.0, v_t_kms=5000.0, v_max_kms=30000.0,
                     m=2.0, n=7.0,
                     t0_days=1.0, T0_inner=2.0e6, gamma=4.0/3.0,
                     v_H_cut_kms=8000.0,
                     fx_H_inner=0.02, fx_He_inner=0.38,
                     fx_O_inner=0.18, fx_Si_inner=0.07,
                     fx_Ca_inner=0.02, fx_Fe_inner=0.07,
                     fx_H_outer=0.60, fx_He_outer=0.38,
                     fx_O_outer=0.01, fx_Si_outer=0.003,
                     fx_Ca_outer=0.001, fx_Fe_outer=0.001):
    """
    Generate Type IIb ejecta profile at time t_days.

    Composition is velocity-shell based:
    - v > v_H_cut: H-rich outer layer (thin H envelope of Type IIb)
    - v <= v_H_cut: He/metal-rich inner ejecta (He core + deeper layers)

    Returns dict with r, v, rho, T, composition arrays, and diagnostics.
    """
    t = t_days * DAY
    t0 = t0_days * DAY

    v_min = v_min_kms * KM_S
    v_t = v_t_kms * KM_S
    v_max = v_max_kms * KM_S
    v_H_cut = v_H_cut_kms * KM_S

    # Velocity grid: log-spaced for accurate integration
    v = np.logspace(np.log10(v_min), np.log10(v_max), 5000)

    # Solve for v_t from E_kin/M ratio, then compute rho_0 from mass constraint
    if E_kin > 0:
        v_t, rho_0 = solve_v_t_for_energy(v, m, n, E_kin, M_ej * MSUN, t0, v_max)
    else:
        # No KE constraint: use specified v_t
        v_t = v_t_kms * KM_S
        I_m, _ = mass_energy_integrals(v, v_t, m, n)
        rho_0 = M_ej * MSUN / (4.0 * np.pi * t0**3 * I_m)
    rho_t0 = broken_powerlaw_density(v, rho_0, v_t, m, n)

    # Homologous evolution to time t
    r = v * t
    rho = rho_t0 * (t0 / t)**3

    # Temperature: radiation-dominated adiabatic, T proportional to t^-1
    # At t0, adiabatic profile: T(v) = T0_inner * (rho(v)/rho_max)^(gamma-1)
    T_t0 = T0_inner * (rho_t0 / rho_t0.max())**(gamma - 1.0)
    T = T_t0 * (t0 / t)

    # Composition: velocity-shell based (Type IIb structure)
    mask_H = v >= v_H_cut
    fx_H = np.where(mask_H, fx_H_outer, fx_H_inner)
    fx_He = np.where(mask_H, fx_He_outer, fx_He_inner)
    fx_O = np.where(mask_H, fx_O_outer, fx_O_inner)
    fx_Si = np.where(mask_H, fx_Si_outer, fx_Si_inner)
    fx_Ca = np.where(mask_H, fx_Ca_outer, fx_Ca_inner)
    fx_Fe = np.where(mask_H, fx_Fe_outer, fx_Fe_inner)

    # Normalize compositions to sum to 1
    comp_sum = fx_H + fx_He + fx_O + fx_Si + fx_Ca + fx_Fe
    fx_H /= comp_sum
    fx_He /= comp_sum
    fx_O /= comp_sum
    fx_Si /= comp_sum
    fx_Ca /= comp_sum
    fx_Fe /= comp_sum

    # Diagnostics: mass and KE
    dr = np.diff(r)
    r_mid = 0.5 * (r[:-1] + r[1:])
    rho_mid = 0.5 * (rho[:-1] + rho[1:])
    v_mid = 0.5 * (v[:-1] + v[1:])
    dm = 4.0 * np.pi * r_mid**2 * rho_mid * dr
    mass_total = np.sum(dm)
    ke = 0.5 * np.sum(dm * v_mid**2)

    # H mass fraction
    fx_H_mid = 0.5 * (fx_H[:-1] + fx_H[1:])
    M_H = np.sum(dm * fx_H_mid)

    return {
        'v': v,
        'r': r,
        'rho': rho,
        'T': T,
        'fx_H': fx_H,
        'fx_He': fx_He,
        'fx_O': fx_O,
        'fx_Si': fx_Si,
        'fx_Ca': fx_Ca,
        'fx_Fe': fx_Fe,
        'rho_0': rho_0,
        'v_t': v_t,
        'mass_total': mass_total,
        'ke': ke,
        'M_H': M_H,
    }


def plot_profiles(times_days, output_path):
    """Plot density vs radius and density vs velocity at multiple times."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = plt.cm.viridis(np.linspace(0, 1, len(times_days)))

    for t, c in zip(times_days, colors):
        prof = type_iib_profile(t)
        r = prof['r']
        v = prof['v']
        rho = prof['rho']
        T = prof['T']

        axes[0, 0].loglog(r, rho, color=c, label=f'{t} d', lw=1.5)
        axes[0, 1].loglog(v / KM_S, rho, color=c, label=f'{t} d', lw=1.5)
        axes[1, 0].loglog(r, T, color=c, label=f'{t} d', lw=1.5)
        axes[1, 1].semilogx(v / KM_S, prof['fx_H'], color=c, label=f'{t} d H', lw=1.5)
        axes[1, 1].semilogx(v / KM_S, prof['fx_He'], color=c, ls='--', label=f'{t} d He', lw=1.5)

    axes[0, 0].set_xlabel('Radius r [cm]')
    axes[0, 0].set_ylabel(r'Density $\rho$ [g cm$^{-3}$]')
    axes[0, 0].set_title('Density vs Radius (homologous expansion)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, which='both', ls='--', alpha=0.5)

    axes[0, 1].set_xlabel('Velocity v [km s$^{-1}$]')
    axes[0, 1].set_ylabel(r'Density $\rho$ [g cm$^{-3}$]')
    axes[0, 1].set_title('Density vs Velocity (profile shape preserved)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, which='both', ls='--', alpha=0.5)

    axes[1, 0].set_xlabel('Radius r [cm]')
    axes[1, 0].set_ylabel('Temperature T [K]')
    axes[1, 0].set_title('Temperature vs Radius (adiabatic, T proportional to t^-1)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, which='both', ls='--', alpha=0.5)

    axes[1, 1].set_xlabel('Velocity v [km s$^{-1}$]')
    axes[1, 1].set_ylabel('Mass fraction')
    axes[1, 1].set_title('Composition vs Velocity (H-rich outer, He/metal inner)')
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(True, which='both', ls='--', alpha=0.5)

    fig.suptitle('Type IIb SN Ejecta: Homologous Evolution (broken power law in velocity space)',
                 fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"Saved plot to {output_path}")


def write_table(times_days, output_path):
    """Write a table with t, v, r, rho, T for each time."""
    with open(output_path, 'w') as f:
        f.write("# t_days  v_kms  r_cm  rho_g_cm3  T_K  fx_H  fx_He  fx_O\n")
        for t in times_days:
            prof = type_iib_profile(t)
            idx = np.linspace(0, len(prof['v']) - 1, 200, dtype=int)
            for i in idx:
                f.write(f"{t:.1f}  {prof['v'][i]/KM_S:.3e}  "
                        f"{prof['r'][i]:.5e}  {prof['rho'][i]:.5e}  "
                        f"{prof['T'][i]:.5e}  {prof['fx_H'][i]:.4f}  "
                        f"{prof['fx_He'][i]:.4f}  {prof['fx_O'][i]:.4f}\n")
    print(f"Saved table to {output_path}")


if __name__ == '__main__':
    times = [50, 100, 200, 400]
    outdir = Path('/viper/ptmp1/sbrennan/SN2024pba/plots')
    outdir.mkdir(parents=True, exist_ok=True)

    plot_path = outdir / 'type_iib_homologous_profiles.png'
    table_path = outdir / 'type_iib_homologous_profiles.dat'

    plot_profiles(times, plot_path)
    write_table(times, table_path)

    print("\n=== Type IIb Ejecta Profile Summary ===")
    print("Broken power law: rho(v) = rho_0 * (v/v_t)^(-m) for v<v_t, (v/v_t)^(-n) for v>v_t")
    print("  m=2 (shallow inner), n=7 (steep outer), v_t=5000 km/s")
    print("  M_ej=2.0 Msun, E_kin=1e51 erg, v_max=30000 km/s")
    print("  Composition: H-rich for v>8000 km/s, He/metal-rich for v<8000 km/s")
    print()
    for t in times:
        prof = type_iib_profile(t)
        print(f"t = {t} days:")
        print(f"  M = {prof['mass_total']/MSUN:.3f} Msun (target: 2.0)")
        print(f"  KE = {prof['ke']/1e51:.3f} x 10^51 erg (target: 1.0)")
        print(f"  M_H = {prof['M_H']/MSUN:.4f} Msun")
        print(f"  R range = {prof['r'].min():.3e} - {prof['r'].max():.3e} cm")
        print(f"  v range = {prof['v'].min()/KM_S:.0f} - {prof['v'].max()/KM_S:.0f} km/s")
        print(f"  rho range = {prof['rho'].min():.3e} - {prof['rho'].max():.3e} g/cm3")
        print(f"  T range = {prof['T'].min():.0f} - {prof['T'].max():.0f} K")
        print()

    print("Note: The ejecta become more extended in radius and more dilute in")
    print("density with time. The profile shape in velocity space is preserved")
    print("(homologous expansion). Temperature decreases as t^-1 (radiation-dominated).")
