#!/usr/bin/env python3
"""
Generate synthetic SN_HYDRO_DATA for a Type II-P supernova.
Creates a realistic hydrodynamic profile with 100 grid points.
"""

import numpy as np
import sys

def generate_type2p_profile(n_points=100):
    """
    Generate a synthetic Type II-P SN hydrodynamic profile.
    
    Parameters:
    -----------
    n_points : int
        Number of radial grid points (default: 100)
    
    Returns:
    --------
    dict containing all profile data
    """
    # Physical constants (based on example SN_HYDRO_DATA)
    R_MIN = 2.3e14  # cm (inner boundary)
    R_MAX = 4.7e15  # cm (outer boundary)
    V_INF = 10000.0  # km/s (outer velocity)
    T_INNER = 1.5  # 10^4 K (inner temperature)
    T_OUTER = 0.5  # 10^4 K (outer temperature)
    RHO_INNER = 6e-13  # g/cm^3 (inner density)
    RHO_OUTER = 8e-17  # g/cm^3 (outer density)
    
    # Create logarithmic radial grid (CMFGEN expects outer to inner)
    r = np.logspace(np.log10(R_MAX), np.log10(R_MIN), n_points) / 1e10  # in 10^10 cm
    
    # Velocity profile (Hubble-like expansion)
    # Velocity decreases from outer to inner
    v = V_INF * (r / r[0])**0.9  # km/s
    
    # Sigma (velocity gradient)
    # Example shows sigma = 0.0 everywhere
    sigma = 0.0 * np.ones(n_points)
    
    # Temperature profile (hot inner, cool outer)
    # Since r goes from outer to inner, we need to invert the profile
    # Inner (small r) should be hot, outer (large r) should be cool
    t_profile = T_INNER * (r / r[-1])**(-0.3)
    t_profile = np.maximum(t_profile, T_OUTER)
    
    # Density profile (steep gradient)
    # Inner (small r) should be dense, outer (large r) should be tenuous
    rho = RHO_INNER * (r / r[-1])**(-2.5)
    rho = np.maximum(rho, RHO_OUTER)
    
    # Electron density (assuming partial ionization)
    # For a Type II-P, ionization fraction ~0.1-0.5
    # Inner (small r) should be more ionized due to higher temperature
    ionization_frac = 0.3 * (r / r[-1])**(-0.2)
    ionization_frac = np.clip(ionization_frac, 0.1, 0.8)
    
    # Convert density to number density (assuming mean molecular weight ~1.3)
    mean_molecular_weight = 1.3
    m_p = 1.6726e-24  # proton mass in g
    atom_density = rho / (mean_molecular_weight * m_p)
    electron_density = ionization_frac * atom_density
    
    # Kappa (opacity) - more realistic model
    # Electron scattering + bound-free with proper temperature dependence
    kappa_es = 0.34  # cm^2/g (electron scattering)
    # Bound-free opacity: scales as T^-3.5, with floor for cool regions
    kappa_bf = 1.0e6 * (t_profile * 1e4)**(-3.5)
    kappa_bf = np.clip(kappa_bf, 0.01, 100.0)
    # Total opacity with smooth transition
    kappa = kappa_es + kappa_bf
    # Ensure reasonable range
    kappa = np.clip(kappa, 0.3, 50.0)
    
    # Mass fractions (based on example SN_HYDRO_DATA)
    # Only HYD, HE, OXY, CAL have non-zero values
    # HYD and HE vary with radius (H decreases inward, He increases inward)
    # OXY and CAL are small and constant
    hyd_mf = 0.70 * (r / r[0])**0.1  # H decreases slightly inward
    hyd_mf = np.clip(hyd_mf, 0.4, 0.75)
    he_mf = 0.28 * (r / r[-1])**0.1  # He increases slightly inward
    he_mf = np.clip(he_mf, 0.25, 0.60)
    oxy_mf = 0.02 * np.ones(n_points)  # O constant
    cal_mf = 5e-5 * np.ones(n_points)  # Ca constant
    
    mass_fractions = {
        'HYD': hyd_mf,
        'HE': he_mf,
        'OXY': oxy_mf,
        'CAL': cal_mf
    }
    
    # Rosseland mean opacity (10^-10 cm^-1)
    # Based on example: scales with density and temperature
    rosseland = kappa * rho * 1e10  # Convert from cm^2/g to 10^-10 cm^-1
    
    return {
        'n_points': n_points,
        'n_mass_fractions': 26,  # Total species count (most zero)
        'n_isotopes': 0,
        'time_days': 53.8,  # Based on example
        'lum_lsun': 1e8,
        'r': r,
        'v': v,
        'sigma': sigma,
        't': t_profile,
        'rho': rho,
        'electron_density': electron_density,
        'atom_density': atom_density,
        'kappa': kappa,
        'rosseland': rosseland,
        'mass_fractions': mass_fractions
    }

def write_sn_hydro_data(profile, filename):
    """
    Write profile to SN_HYDRO_DATA format matching Dessart's reference format.

    Header titles and data format match the reference SN_HYDRO_DATA from
    https://ttt.astro.su.se/~maer0651/M15-He-O-Ca-evo-new-4-all/34/SN_HYDRO_DATA

    Parameters:
    -----------
    profile : dict
        Profile data from generate_type2p_profile
    filename : str
        Output filename
    """
    def write_block(f, label, values, per_line=8):
        f.write(f"{label}\n")
        for i in range(0, len(values), per_line):
            chunk = values[i:i+per_line]
            f.write(''.join(f"    {v:14.7E}" for v in chunk) + '\n')
        f.write('\n')

    with open(filename, 'w') as f:
        # Header (SN_HYDRO_DATA format)
        f.write(f"Number of data points:{profile['n_points']:>10d}\n")
        f.write(f"Number of mass fractions:{profile['n_mass_fractions']:>10d}\n")
        f.write(f"Number of isotopes:{profile['n_isotopes']:>10d}\n")
        f.write(f"Time(days) since explosion:{profile['time_days']:>15.7f}\n")
        # Extra header line used by our pipeline; read by prepare_VADAT.py,
        # run_model_grid.py, fix_sedona_param.py.
        f.write(f"Energy Lum [lsun]:{profile['lum_lsun']:>19.6E}\n")
        f.write("\n")

        # Data blocks (SN_HYDRO_DATA order and titles)
        write_block(f, "Radius grid (10^10cm)", profile['r'])
        write_block(f, "Velocity (km/s)",       profile['v'])
        write_block(f, "Sigma (dlnV/dlnr-1)",   profile['sigma'])
        write_block(f, "Temperature (10^4 K)",  profile['t'])
        write_block(f, "Density (gm/cm^3)",     profile['rho'])
        write_block(f, "Atom density (/cm^3)",  profile['atom_density'])
        write_block(f, "Electron density (/cm^3)", profile['electron_density'])
        write_block(f, "Rosseland mean opacty (10^{-10} cm^{-1})", profile['rosseland'])
        write_block(f, "Kappa (cm^2/gm)",       profile['kappa'])

        # Mass fractions (26 species total, most zero)
        species_order = ['HYD', 'HE', 'CARB', 'NIT', 'OXY', 'FLU', 'NEON', 'SOD',
                         'MAG', 'ALUM', 'SIL', 'PHOS', 'SUL', 'CHL', 'ARG', 'POT',
                         'CAL', 'SCAN', 'TIT', 'VAN', 'CHRO', 'MAN', 'IRON', 'COB', 'NICK', 'BAR']
        for species in species_order:
            if species in profile['mass_fractions']:
                mf = profile['mass_fractions'][species]
            else:
                mf = np.zeros(profile['n_points'])
            write_block(f, f"{species} mass fraction", mf)

def main():
    """Main function."""
    if len(sys.argv) > 1:
        output_file = sys.argv[1]
    else:
        output_file = "SN_HYDRO_DATA"
    
    print(f"Generating synthetic Type II-P SN_HYDRO_DATA with 100 grid points...")
    profile = generate_type2p_profile(n_points=100)
    write_sn_hydro_data(profile, output_file)
    print(f"Written to {output_file}")
    
    # Print some statistics
    print(f"\nProfile statistics:")
    print(f"  Radius range: {profile['r'][0]:.2e} to {profile['r'][-1]:.2e} (10^10 cm)")
    print(f"  Velocity range: {profile['v'][0]:.2e} to {profile['v'][-1]:.2e} km/s")
    print(f"  Temperature range: {profile['t'][0]:.2e} to {profile['t'][-1]:.2e} (10^4 K)")
    print(f"  Density range: {profile['rho'][0]:.2e} to {profile['rho'][-1]:.2e} g/cm^3")
    print(f"  Electron density range: {profile['electron_density'][0]:.2e} to {profile['electron_density'][-1]:.2e} cm^-3")

if __name__ == "__main__":
    main()
