#!/usr/bin/env python3
"""
Smooth Sigma values in existing SN_HYDRO_DATA file to prevent CMFGEN segmentation faults.
Uses Savitzky-Golay filtering to spread sharp velocity gradients over more grid points.
"""

import sys
import numpy as np
from scipy.signal import savgol_filter


def read_sn_hydro_data(filename):
    """
    Read SN_HYDRO_DATA file and return all data arrays.
    
    Returns:
    --------
    Dictionary containing:
        - n_points: Number of data points
        - n_mass_fractions: Number of mass fractions
        - n_isotopes: Number of isotopes
        - time_days: Time since explosion (days)
        - energy_lum: Energy luminosity (L_sun)
        - radius: Radius array (10^10 cm)
        - velocity: Velocity array (km/s)
        - sigma: Sigma array (dlnV/dlnr - 1)
        - temperature: Temperature array (K)
        - density: Density array (g/cm³)
        - mass_fractions: Dictionary of mass fraction arrays
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header lines up to the first data section. Ejecta-only files
    # (e.g. Stockholm grids) omit the Energy Lum header entirely.
    n_points = n_mass_fractions = n_isotopes = None
    time_days = energy_lum = None
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith('Radius grid'):
            break
        if s.startswith('Number of data points'):
            n_points = int(s.split(':')[1].strip())
        elif s.startswith('Number of mass fractions'):
            n_mass_fractions = int(s.split(':')[1].strip())
        elif s.startswith('Number of isotopes'):
            n_isotopes = int(s.split(':')[1].strip())
        elif s.startswith('Time(days) since explosion'):
            time_days = float(s.split(':')[1].strip())
        elif s.startswith('Energy Lum'):
            energy_lum = float(s.split(':')[1].strip())
        i += 1
    assert n_points is not None, "Missing 'Number of data points' header"
    assert n_mass_fractions is not None, "Missing 'Number of mass fractions' header"
    assert n_isotopes is not None, "Missing 'Number of isotopes' header"
    assert time_days is not None, "Missing 'Time(days) since explosion' header"

    # Read radius
    assert 'Radius grid' in lines[i], f"Expected 'Radius grid' at line {i}"
    i += 1
    radius_data = []
    while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Velocity'):
        radius_data.extend([float(x) for x in lines[i].split()])
        i += 1
    radius = np.array(radius_data)
    while i < len(lines) and not lines[i].strip():
        i += 1  # skip blank separators
    
    # Read velocity
    assert 'Velocity' in lines[i], f"Expected 'Velocity' at line {i}"
    i += 1
    velocity_data = []
    while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Sigma'):
        velocity_data.extend([float(x) for x in lines[i].split()])
        i += 1
    velocity = np.array(velocity_data)
    while i < len(lines) and not lines[i].strip():
        i += 1  # skip blank separators
    
    # Read sigma
    assert 'Sigma' in lines[i], f"Expected 'Sigma' at line {i}"
    i += 1
    sigma_data = []
    while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Temperature'):
        sigma_data.extend([float(x) for x in lines[i].split()])
        i += 1
    sigma = np.array(sigma_data)
    while i < len(lines) and not lines[i].strip():
        i += 1  # skip blank separators
    
    # Read temperature
    assert 'Temperature' in lines[i], f"Expected 'Temperature' at line {i}"
    i += 1
    temperature_data = []
    while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Density'):
        temperature_data.extend([float(x) for x in lines[i].split()])
        i += 1
    temperature = np.array(temperature_data)
    while i < len(lines) and not lines[i].strip():
        i += 1  # skip blank separators
    
    # Read density
    assert 'Density' in lines[i], f"Expected 'Density' at line {i}"
    i += 1
    density_data = []
    while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Atom') \
            and not lines[i].strip().startswith('Electron') and not lines[i].strip().startswith('Rosseland') \
            and not lines[i].strip().startswith('Kappa') and not ('mass fraction' in lines[i].lower()):
        density_data.extend([float(x) for x in lines[i].split()])
        i += 1
    density = np.array(density_data)
    while i < len(lines) and not lines[i].strip():
        i += 1  # skip blank separators

    # Skip atom density section (comes before electron density in Dessart format)
    if i < len(lines) and 'Atom' in lines[i]:
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Electron') \
                and not lines[i].strip().startswith('Rosseland') and not lines[i].strip().startswith('Kappa') \
                and not ('mass fraction' in lines[i].lower()):
            i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1  # skip blank separators

    # Skip electron density section
    if i < len(lines) and 'Electron' in lines[i]:
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Rosseland') \
                and not lines[i].strip().startswith('Kappa') and not ('mass fraction' in lines[i].lower()):
            i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1  # skip blank separators

    # Skip Rosseland mean opacity section
    if i < len(lines) and 'Rosseland' in lines[i]:
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('Kappa') \
                and not ('mass fraction' in lines[i].lower()):
            i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1  # skip blank separators

    # Skip kappa section
    if i < len(lines) and 'Kappa' in lines[i]:
        i += 1
        while i < len(lines) and lines[i].strip() and not ('mass fraction' in lines[i].lower()):
            i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1  # skip blank separators
    
    # Read mass fractions (individual sections like "HYD mass fraction")
    mass_fractions = {}
    for _ in range(n_mass_fractions):
        # Find next mass fraction section
        while i < len(lines) and not ('mass fraction' in lines[i].lower()):
            i += 1
        
        if i >= len(lines):
            break
        
        # Extract element name from line like "HYD mass fraction"
        elem_line = lines[i].strip()
        elem_name = elem_line.split(' mass fraction')[0].strip()
        i += 1
        
        # Read data until empty line or next section
        mass_frac_data = []
        while i < len(lines) and lines[i].strip() and not ('mass fraction' in lines[i].lower()):
            mass_frac_data.extend([float(x) for x in lines[i].split()])
            i += 1
            # Stop if we've collected enough data points
            if len(mass_frac_data) >= n_points:
                break
        
        mass_fractions[elem_name] = np.array(mass_frac_data[:n_points])
        
        # Skip empty line if present
        if i < len(lines) and not lines[i].strip():
            i += 1
    
    return {
        'n_points': n_points,
        'n_mass_fractions': n_mass_fractions,
        'n_isotopes': n_isotopes,
        'time_days': time_days,
        'energy_lum': energy_lum,
        'radius': radius,
        'velocity': velocity,
        'sigma': sigma,
        'temperature': temperature,
        'density': density,
        'mass_fractions': mass_fractions
    }


def write_sn_hydro_data(filename, data):
    """
    Write SN_HYDRO_DATA file with smoothed Sigma values.
    Uses the SN_HYDRO_DATA header titles and data format.
    """
    def write_block(f, label, values, per_line=8):
        f.write(f"{label}\n")
        for i in range(0, len(values), per_line):
            chunk = values[i:i+per_line]
            f.write(''.join(f"    {v:14.7E}" for v in chunk) + '\n')
        f.write('\n')

    with open(filename, 'w') as f:
        # Header (SN_HYDRO_DATA format)
        f.write(f"Number of data points:{data['n_points']:>10d}\n")
        f.write(f"Number of mass fractions:{data['n_mass_fractions']:>10d}\n")
        f.write(f"Number of isotopes:{data['n_isotopes']:>10d}\n")
        f.write(f"Time(days) since explosion:{data['time_days']:>15.7f}\n")
        if data.get('energy_lum') is not None:
            f.write(f"Energy Lum [lsun]:{data['energy_lum']:>19.6E}\n")
        f.write("\n")

        # Data blocks (SN_HYDRO_DATA order and titles)
        write_block(f, "Radius grid (10^10cm)", data['radius'])
        write_block(f, "Velocity (km/s)",       data['velocity'])
        write_block(f, "Sigma (dlnV/dlnr-1)",   data['sigma'])
        write_block(f, "Temperature (10^4 K)",  data['temperature'])
        write_block(f, "Density (gm/cm^3)",     data['density'])

        # Mass fractions
        for elem in sorted(data['mass_fractions'].keys()):
            write_block(f, f"{elem} mass fraction", data['mass_fractions'][elem])


def smooth_sigma_savgol(r, velocity, window_length=9, polyorder=3, threshold=2.0, target_max=5.0):
    """
    Apply Savitzky-Golay smoothing to velocity profile and recalculate Sigma.
    
    Parameters:
    -----------
    r : array
        Radius values (10^10 cm)
    velocity : array
        Velocity values (km/s)
    window_length : int
        Length of the filter window (must be odd, default: 9)
    polyorder : int
        Order of the polynomial used to fit the samples (default: 3)
    threshold : float
        Sigma threshold above which velocity smoothing is applied (default: 2.0)
    target_max : float
        Target maximum Sigma value after smoothing (default: 5.0)
    
    Returns:
    --------
    sigma_recalculated : array
        Recalculated Sigma values from smoothed velocity
    velocity_smoothed : array
        Smoothed velocity values
    """
    # Convert velocity to cm/s for Sigma calculation
    v_cgs = velocity * 1e5  # km/s to cm/s
    r_cgs = r * 1e10  # 10^10 cm to cm
    
    # Calculate original Sigma
    sigma_original = np.zeros_like(r)
    for i in range(1, len(r)-1):
        dlnv = np.log(v_cgs[i+1]) - np.log(v_cgs[i-1])
        dlnr = np.log(r_cgs[i+1]) - np.log(r_cgs[i-1])
        if abs(dlnr) > 1e-10:
            sigma_original[i] = dlnv / dlnr - 1
    sigma_original[0] = sigma_original[1]
    sigma_original[-1] = sigma_original[-2]
    
    velocity_smoothed = v_cgs.copy()
    
    # Ensure window_length is odd
    if window_length % 2 == 0:
        window_length += 1
    
    # Ensure window_length is not larger than array
    if window_length >= len(velocity_smoothed):
        window_length = len(velocity_smoothed) - 1
        if window_length % 2 == 0:
            window_length -= 1
    
    # Identify shock regions where Sigma > threshold
    shock_mask = np.abs(sigma_original) > threshold
    
    if not np.any(shock_mask):
        print("[Sigma Smoothing] No shock regions detected (Sigma < threshold)")
        return sigma_original, velocity
    
    # Find contiguous shock regions
    shock_regions = []
    in_shock = False
    shock_start = 0
    
    for i in range(len(sigma_original)):
        if shock_mask[i] and not in_shock:
            in_shock = True
            shock_start = i
        elif not shock_mask[i] and in_shock:
            in_shock = False
            shock_regions.append((shock_start, i))
    
    if in_shock:
        shock_regions.append((shock_start, len(sigma_original)))
    
    # Apply Savitzky-Golay smoothing to velocity in each shock region with padding
    for start, end in shock_regions:
        # Add padding around shock region for smooth transitions
        pad = window_length // 2
        padded_start = max(0, start - pad)
        padded_end = min(len(velocity_smoothed), end + pad)
        
        # Apply Savitzky-Golay filter to padded region
        region_velocity = velocity_smoothed[padded_start:padded_end]
        if len(region_velocity) >= window_length:
            smoothed_region = savgol_filter(region_velocity, window_length, polyorder)
            velocity_smoothed[padded_start:padded_end] = smoothed_region
    
    # Print statistics before peak reduction
    print(f"[Sigma Smoothing] Applied Savitzky-Golay smoothing to velocity")
    print(f"[Sigma Smoothing] Window length: {window_length}, Polynomial order: {polyorder}")
    print(f"[Sigma Smoothing] Threshold: {threshold}")
    print(f"[Sigma Smoothing] Found {len(shock_regions)} shock regions")
    print(f"[Sigma Smoothing] Original Sigma range: {sigma_original.min():.3f} to {sigma_original.max():.3f}")
    print(f"[Sigma Smoothing] Velocity range before: {v_cgs.min()/1e5:.3f} to {v_cgs.max()/1e5:.3f} km/s")
    print(f"[Sigma Smoothing] Velocity range after:  {velocity_smoothed.min()/1e5:.3f} to {velocity_smoothed.max()/1e5:.3f} km/s")
    
    # Recalculate Sigma from smoothed velocity
    sigma_recalculated = np.zeros_like(r)
    for i in range(1, len(r)-1):
        dlnv = np.log(velocity_smoothed[i+1]) - np.log(velocity_smoothed[i-1])
        dlnr = np.log(r_cgs[i+1]) - np.log(r_cgs[i-1])
        if abs(dlnr) > 1e-10:
            sigma_recalculated[i] = dlnv / dlnr - 1
    sigma_recalculated[0] = sigma_recalculated[1]
    sigma_recalculated[-1] = sigma_recalculated[-2]
    
    print(f"[Sigma Smoothing] Recalculated Sigma range:  {sigma_recalculated.min():.3f} to {sigma_recalculated.max():.3f}")
    
    # Apply additional Savitzky-Golay smoothing to the recalculated Sigma values
    # This further reduces sharp gradients that remain after velocity smoothing
    shock_mask = np.abs(sigma_recalculated) > threshold
    if np.any(shock_mask):
        for start, end in shock_regions:
            pad = window_length // 2
            padded_start = max(0, start - pad)
            padded_end = min(len(sigma_recalculated), end + pad)
            
            region_sigma = sigma_recalculated[padded_start:padded_end]
            if len(region_sigma) >= window_length:
                smoothed_region = savgol_filter(region_sigma, window_length, polyorder)
                sigma_recalculated[padded_start:padded_end] = smoothed_region
        
        print(f"[Sigma Smoothing] Applied additional Sigma smoothing")
        print(f"[Sigma Smoothing] Sigma range after additional smoothing:  {sigma_recalculated.min():.3f} to {sigma_recalculated.max():.3f}")
    
    # Reduce peak values to target_max by further smoothing velocity
    if np.max(np.abs(sigma_recalculated)) > target_max:
        # Scale down velocity in regions with high Sigma
        scale_factor = target_max / np.max(np.abs(sigma_recalculated))
        mask = np.abs(sigma_recalculated) > target_max
        # Scale velocity toward the local mean to reduce Sigma
        for i in np.where(mask)[0]:
            # Use local average of velocity
            start_idx = max(0, i - window_length//2)
            end_idx = min(len(velocity_smoothed), i + window_length//2 + 1)
            local_mean = np.mean(velocity_smoothed[start_idx:end_idx])
            velocity_smoothed[i] = velocity_smoothed[i] * (1 - 0.5 * (1 - scale_factor)) + local_mean * 0.5 * (1 - scale_factor)
        
        # Recalculate Sigma after velocity adjustment
        sigma_recalculated = np.zeros_like(r)
        for i in range(1, len(r)-1):
            dlnv = np.log(velocity_smoothed[i+1]) - np.log(velocity_smoothed[i-1])
            dlnr = np.log(r_cgs[i+1]) - np.log(r_cgs[i-1])
            if abs(dlnr) > 1e-10:
                sigma_recalculated[i] = dlnv / dlnr - 1
        sigma_recalculated[0] = sigma_recalculated[1]
        sigma_recalculated[-1] = sigma_recalculated[-2]
        
        print(f"[Sigma Smoothing] Further velocity adjustment applied")
        print(f"[Sigma Smoothing] Final Sigma range:  {sigma_recalculated.min():.3f} to {sigma_recalculated.max():.3f}")
        print(f"[Sigma Smoothing] Final velocity range:  {velocity_smoothed.min()/1e5:.3f} to {velocity_smoothed.max()/1e5:.3f} km/s")
    
    return sigma_recalculated, velocity_smoothed / 1e5  # Return velocity in km/s


def main():
    if len(sys.argv) < 2:
        print("Usage: python smooth_sigma_hydro_data.py <SN_HYDRO_DATA_file> [window_length] [polyorder] [threshold] [target_max]")
        print("  window_length: Savitzky-Golay window length (default: 9)")
        print("  polyorder: Polynomial order (default: 3)")
        print("  threshold: Sigma threshold for smoothing (default: 2.0)")
        print("  target_max: Target maximum Sigma value after smoothing (default: 5.0)")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Parse optional parameters
    window_length = 9
    polyorder = 3
    threshold = 2.0
    target_max = 5.0
    
    if len(sys.argv) > 2:
        window_length = int(sys.argv[2])
    if len(sys.argv) > 3:
        polyorder = int(sys.argv[3])
    if len(sys.argv) > 4:
        threshold = float(sys.argv[4])
    if len(sys.argv) > 5:
        target_max = float(sys.argv[5])
    
    print(f"Reading SN_HYDRO_DATA from {input_file}...")
    data = read_sn_hydro_data(input_file)
    
    print(f"  Number of points: {data['n_points']}")
    print(f"  Time since explosion: {data['time_days']:.2f} days")
    print(f"  Sigma range before smoothing: {data['sigma'].min():.3f} to {data['sigma'].max():.3f}")
    
    # Apply smoothing
    print(f"\nApplying Sigma smoothing...")
    sigma_smoothed, velocity_smoothed = smooth_sigma_savgol(data['radius'], data['velocity'], 
                                                              window_length=window_length, 
                                                              polyorder=polyorder, 
                                                              threshold=threshold, 
                                                              target_max=target_max)
    
    # Update data
    data['sigma'] = sigma_smoothed
    data['velocity'] = velocity_smoothed
    
    # Write output
    output_file = input_file.replace('.dat', '_smoothed.dat').replace('SN_HYDRO_DATA', 'SN_HYDRO_DATA_smoothed')
    print(f"\nWriting smoothed data to {output_file}...")
    write_sn_hydro_data(output_file, data)
    
    print(f"\nDone! Smoothed Sigma range: {sigma_smoothed.min():.3f} to {sigma_smoothed.max():.3f}")
    print(f"Smoothed velocity range: {velocity_smoothed.min():.3f} to {velocity_smoothed.max():.3f} km/s")
    print(f"Original file preserved at {input_file}")


if __name__ == '__main__':
    main()
