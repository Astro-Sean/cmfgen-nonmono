#!/usr/bin/env python3
"""
Check the status of a running STELLA simulation.

Run this script in a model directory (e.g., /ptmp/sbrennan/SN2024pba/models/X_STELA)
to see the current state of the STELLA job.
"""

import os
import re
import subprocess
import sys
from pathlib import Path


def find_stella_output():
    """Find the most recent STELLA output file in the stella/ subdirectory."""
    stella_dir = Path("stella")
    if not stella_dir.exists():
        print(f"ERROR: stella/ directory not found in {Path.cwd()}")
        sys.exit(1)
    
    # Look for stella_*.out files
    out_files = list(stella_dir.glob("stella_*.out"))
    if not out_files:
        print(f"ERROR: No stella_*.out files found in {stella_dir}")
        sys.exit(1)
    
    # Return the most recently modified file
    return max(out_files, key=lambda f: f.stat().st_mtime)


def find_stella_dat_file():
    """Find the STELLA .dat file (in strad/ or modmake/ subdirectory)."""
    stella_dir = Path("stella")
    
    # Try strad/ first, then modmake/
    for subdir in ["strad", "modmake"]:
        dat_dir = stella_dir / subdir
        if dat_dir.exists():
            dat_files = list(dat_dir.glob("*.dat"))
            if dat_files:
                return dat_files[0]  # Return the first .dat file found
    
    return None


def parse_terminal_time(dat_file):
    """Parse the terminal time from the STELLA .dat file."""
    if not dat_file or not dat_file.exists():
        return None
    
    try:
        with open(dat_file, 'r') as f:
            lines = f.readlines()
        
        # Look for NTO section (output times)
        # Format: "    NTO then    DATA TO" followed by count and list of times
        nto_section_started = False
        nto_count = 0
        times = []
        
        for line in lines:
            line = line.strip()
            if 'NTO then' in line or 'NTO' in line:
                nto_section_started = True
                continue
            if nto_section_started:
                # Try to parse as a number (time in days)
                try:
                    time_val = float(line)
                    times.append(time_val)
                except ValueError:
                    # If we hit a non-numeric line, we're done with the NTO section
                    if times:
                        break
        
        # Terminal time is the last value in the NTO list
        if times:
            return times[-1]
        
        # Fallback: try to parse TcurB from line 5-6
        # Format: "   -9000     -1                      0.0000   65.0000"
        for line in lines:
            if 'TcurA' in line or 'TcurB' in line:
                # Next line should have the values
                idx = lines.index(line)
                if idx + 1 < len(lines):
                    next_line = lines[idx + 1].strip()
                    parts = next_line.split()
                    if len(parts) >= 4:
                        try:
                            return float(parts[-1])  # TcurB is the last value
                        except ValueError:
                            pass
        
    except Exception as e:
        print(f"Warning: Could not parse terminal time from .dat file: {e}")
    
    return None


def parse_stella_output(out_file):
    """Parse the STELLA output file to extract current status."""
    # Use strings to handle binary output files
    result = subprocess.run(
        ["strings", str(out_file)],
        capture_output=True,
        text=True
    )
    
    lines = result.stdout.split('\n')
    
    # Extract latest timestep info (lines with "Hu=" and "td=")
    # These lines look like: "     775    5759  1 Hu= 1.305E-04 Ro= 1.130E-05 td=  18.4"
    timestep_lines = []
    for line in lines:
        if 'Hu=' in line and 'td=' in line and 'Tau=' in line:
            timestep_lines.append(line.strip())
    
    # Use td from timestep line as simulation time (in days)
    time_days = None
    time_seconds = None
    latest_timestep = None
    step = None
    iterations = None
    hu = None
    ro = None
    td = None
    ncond = None
    tau = None
    
    if timestep_lines:
        latest_timestep = timestep_lines[-1]
        # Parse: step, iterations, Hu, Ro, td, Ncnd, Tau
        # Format: "     4976     789  2 Hu= 4.644E-03 Ro= 2.382E-05 td=  27.4     Ncnd=330 Tau=    13."
        try:
            # Extract step and iterations (first two numbers)
            parts = latest_timestep.split()
            step = int(parts[0])
            iterations = int(parts[1])
            # Extract labeled values using regex
            hu_match = re.search(r'Hu=\s*([0-9.eE+-]+)', latest_timestep)
            ro_match = re.search(r'Ro=\s*([0-9.eE+-]+)', latest_timestep)
            td_match = re.search(r'td=\s*([0-9.eE+-]+)', latest_timestep)
            ncond_match = re.search(r'Ncnd=\s*([0-9]+)', latest_timestep)
            tau_match = re.search(r'Tau=\s*([0-9.eE+-]+)', latest_timestep)
            if hu_match:
                hu = float(hu_match.group(1))
            if ro_match:
                ro = float(ro_match.group(1))
            if td_match:
                td = float(td_match.group(1))
                time_days = td
                time_seconds = td * 86400.0
            if ncond_match:
                ncond = int(ncond_match.group(1))
            if tau_match:
                tau = float(tau_match.group(1))
        except (IndexError, ValueError):
            pass

    # Extract timestep size (lines with "H=" and "Hused=")
    # These lines look like: "NSTEP,MAXER,NQUUSED,' H=',H,' Hused=',HUSED,"
    timestep_size_lines = []
    for line in lines:
        if 'H=' in line and 'Hused=' in line:
            timestep_size_lines.append(line.strip())
    
    h_value = None
    hused_value = None
    if timestep_size_lines:
        latest_h_line = timestep_size_lines[-1]
        # Parse: H=... Hused=...
        try:
            h_match = re.search(r'H=\s*([0-9.eE+-]+)', latest_h_line)
            hused_match = re.search(r'Hused=\s*([0-9.eE+-]+)', latest_h_line)
            if h_match:
                h_value = float(h_match.group(1))
            if hused_match:
                hused_value = float(hused_match.group(1))
        except (IndexError, ValueError):
            pass
    
    # Extract photospheric info (lines with "primag Tph Rph Iph=")
    photosphere_lines = []
    for line in lines:
        if 'primag Tph Rph Iph=' in line:
            photosphere_lines.append(line.strip())
    
    latest_photosphere = None
    tph = None
    rph_cm = None
    if photosphere_lines:
        latest_photosphere = photosphere_lines[-1]
        # Parse: primag Tph Rph Iph
        # Format: "primag Tph Rph Iph=         472  0.157922286969477"
        # Note: Tph in STELLA output is in units of 10^4 K (UTp = 1e4)
        # Note: Rph in STELLA output is in units of solar radii (R_sun)
        try:
            # Split on '=' and take the second part
            values_str = latest_photosphere.split('=')[1].strip()
            parts = values_str.split()
            if len(parts) >= 2:
                tph_internal = float(parts[0])  # In units of 10^4 K
                tph = tph_internal * 1e4  # Convert to Kelvin
                rph_solar = float(parts[1])  # In units of solar radii
                rph_cm = rph_solar * 6.957e10  # Convert to cm (solar radius)
                iph = float(parts[2]) if len(parts) >= 3 else None
        except (IndexError, ValueError):
            pass
    
    # Check for job status via squeue
    job_id = None
    job_match = re.search(r'stella_\w+_(\d+)\.out', out_file.name)
    if job_match:
        job_id = job_match.group(1)
    
    job_status = "UNKNOWN"
    job_time = ""
    if job_id:
        try:
            result = subprocess.run(
                ["squeue", "-j", job_id, "-o", "%T %M"],
                capture_output=True,
                text=True
            )
            if result.returncode == 0 and len(result.stdout.strip().split('\n')) > 1:
                status_line = result.stdout.strip().split('\n')[1]
                parts = status_line.split()
                if len(parts) >= 2:
                    job_status = parts[0]
                    job_time = parts[1]
        except Exception:
            pass
    
    # Calculate photospheric velocity if we have Rph and time
    vph_cm_s = None
    vph_km_s = None
    if rph_cm is not None and time_seconds > 0:
        vph_cm_s = rph_cm / time_seconds  # cm/s
        vph_km_s = vph_cm_s / 1e4  # km/s (STELLA uses 10^4 cm/s as velocity unit)
    
    # Parse terminal time from .dat file
    dat_file = find_stella_dat_file()
    terminal_time = parse_terminal_time(dat_file)
    
    return {
        'time_seconds': time_seconds,
        'time_days': time_days,
        'timestep': latest_timestep,
        'step': step,
        'iterations': iterations,
        'hu': hu,
        'ro': ro,
        'td': td,
        'ncond': ncond,
        'tau': tau,
        'h_value': h_value,
        'hused_value': hused_value,
        'photosphere': latest_photosphere,
        'tph': tph,
        'rph_cm': rph_cm,
        'vph_cm_s': vph_cm_s,
        'vph_km_s': vph_km_s,
        'job_id': job_id,
        'job_status': job_status,
        'job_time': job_time,
        'out_file': out_file,
        'terminal_time': terminal_time,
        'dat_file': dat_file
    }


def print_status(status):
    """Print a formatted status summary."""
    print("=" * 70)
    print(f"STELLA Model Status")
    print("=" * 70)
    print(f"Output file: {status['out_file'].name}")
    print(f"Job ID: {status['job_id']}")
    print(f"Job status: {status['job_status']} ({status['job_time']})")
    print()
    
    # Print progress message with terminal time
    current_days = status['time_days']
    terminal_days = status.get('terminal_time')
    if current_days is not None:
        if terminal_days is not None:
            print(f"STELLA is currently at {current_days:.3f} days and will terminate on {terminal_days:.1f} days")
        else:
            print(f"STELLA is currently at {current_days:.3f} days")
    print()
    
    if status['timestep']:
        print("Latest timestep:")
        # Parse and display individual components
        if status.get('step') is not None:
            print(f"  Step number: {status['step']}")
        if status.get('iterations') is not None:
            print(f"  Iterations: {status['iterations']}")
        if status.get('hu') is not None:
            print(f"  Hu (timestep size): {status['hu']:.3e} s")
        if status.get('ro') is not None:
            print(f"  Ro (inner density): {status['ro']:.3e} g/cm^3")
        if status.get('ncond') is not None:
            print(f"  Ncnd (condensed zones): {status['ncond']}")
        if status.get('tau') is not None:
            print(f"  Tau (optical depth): {status['tau']:.1f}")
        print()
    
    if status['photosphere']:
        print("Latest photosphere:")
        if status.get('tph'):
            print(f"  Tph (temperature): {status['tph']:.0f} K")
        if status.get('rph_cm'):
            print(f"  Rph (radius): {status['rph_cm']:.2e} cm")
        if status.get('vph_km_s'):
            print(f"  Vph (expansion velocity): {status['vph_km_s']:.1f} km/s")
        print()
    
    print("=" * 70)


def main():
    if len(sys.argv) > 1:
        # Allow specifying a different directory
        target_dir = Path(sys.argv[1])
        os.chdir(target_dir)
    
    out_file = find_stella_output()
    print(f"Checking STELLA output: {out_file}")
    
    status = parse_stella_output(out_file)
    if status:
        print_status(status)
    else:
        print("ERROR: Could not parse STELLA output")
        sys.exit(1)


if __name__ == '__main__':
    main()
