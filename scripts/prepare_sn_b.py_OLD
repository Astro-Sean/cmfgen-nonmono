#!/usr/bin/env python3
"""
Prepare the SN_B initial model: embed a H/He shell inside the 349-day SN_A ejecta.

This script creates the initial conditions for SN_B, which represents a second
explosion (H/He shell launch) inside the 349-day-old SN_A ejecta.

The SN_B model consists of:
  1. Inner zones: A new H/He shell (1.0 Msun, X=0.7, Y=0.3) launched at 5000 km/s
  2. Outer zones: The 349-day homologously evolved SN_A ejecta (the "CSM")

The H/He shell is placed at the center with an initial velocity of 5000 km/s.
The SN_A ejecta retains its homologous velocity structure.

STELLA starts at t=0 (the shell launch time). The SN_A ejecta is already
349 days old at this point — its structure represents the pre-existing CSM.

Outputs:
  - SN_B.hyd: STELLA hydrodynamics input
  - SN_B.abn: STELLA composition input
  - SN_B.xni: Ni56 distribution (zero — no radioactive decay in SN_B)
  - SN_B.dat: STELLA control file
  - input_profile: 12-column profile for convert_profile_to_stella.py

Usage:
    python prepare_sn_b.py <sn_a_dir> <sn_b_dir> [options]

    --shell-mass 1.0        H/He shell mass in Msun
    --shell-velocity 5000   Shell launch velocity in km/s
    --shell-xh 0.7          H mass fraction in shell
    --shell-yhe 0.3         He mass fraction in shell
    --target-days 50        STELLA simulation duration (days)
    --snapshot-day 40       CMFGEN snapshot day after shell launch
    --nzon 500              Number of STELLA zones
    --sn-a-age 349          Age of SN_A ejecta at shell launch (days)
"""

import argparse
import sys
import struct
import numpy as np
from pathlib import Path

# Import shared utilities
sys.path.insert(0, str(Path(__file__).parent))
from convert_stella_swd_to_cmfgen import read_swd
from homologous_evolution import find_closest_epoch, homologous_evolve

MSUN = 1.989e33
K_B = 1.3807e-16
M_H = 1.6735e-24
A_RAD = 7.5657e-15
KAPPA_ES = 0.2 * (1 + 0.7) / (1 + 0.7 + 0.25 * 0.3)  # electron scattering opacity


def read_sn_a_ejecta(sn_a_dir, sn_a_age=349, tolerance=20):
    """Read the SN_A ejecta profile at sn_a_age days.

    Reads the STELLA .swd file from SN_A, finds the epoch closest to sn_a_age,
    and homologously evolves it if necessary.

    Returns dict with r, rho, v, T, composition (outer->inner order).
    """
    sn_a_dir = Path(sn_a_dir)
    model_name = sn_a_dir.name

    # Find .swd file
    swd_file = sn_a_dir / 'stella' / 'strad' / f'{model_name}.swd'
    if not swd_file.exists():
        swd_file = sn_a_dir / 'stella' / f'{model_name}.swd'
    if not swd_file.exists():
        raise FileNotFoundError(f"SN_A .swd not found: {swd_file}")

    print(f"  Reading SN_A .swd: {swd_file}")
    epochs = read_swd(str(swd_file))
    available = sorted(epochs.keys())
    print(f"  Available epochs: {available[0]:.2f} ... {available[-1]:.2f} days")

    # Find closest epoch
    t_source, epoch_data = find_closest_epoch(epochs, sn_a_age, tolerance)

    # If the closest epoch is not at sn_a_age, homologously evolve
    if abs(t_source - sn_a_age) > 1.0:
        print(f"  Homologously evolving from {t_source:.2f} to {sn_a_age} days...")
        epoch_data = homologous_evolve(epoch_data, t_source, sn_a_age)
    else:
        print(f"  Using epoch at {t_source:.2f} days (matches requested age)")

    return epoch_data


def create_hhe_shell(shell_mass_msun, shell_velocity_kms, xh, yhe, nzon_shell=100,
                     r_inner_cm=1e8, T_shell_K=1e5):
    """Create the H/He shell profile.

    The shell is a uniform-density sphere with:
    - Mass: shell_mass_msun
    - Velocity: shell_velocity_kms (homologous: v = v0 * r / r_outer)
    - Composition: X=xh, Y=yhe
    - Temperature: T_shell_K (isothermal)

    The shell is divided into nzon_shell zones, log-spaced in radius.
    """
    M_shell = shell_mass_msun * MSUN
    v_shell = shell_velocity_kms * 1e5  # cm/s

    # Calculate outer radius from mass and density
    # For a uniform sphere: M = (4/3) * pi * r^3 * rho
    # Choose rho so that the shell is reasonably compact
    # At 5000 km/s, the shell will expand to ~r_outer * (t / t_dyn) in homologous expansion
    # For a 1 Msun shell at 5000 km/s, the dynamical timescale is:
    # t_dyn ~ r / v ~ 1e9 / 5e8 ~ 2 seconds
    # We want the shell to start compact (r_inner ~ 1e8 cm, r_outer ~ 1e10 cm)
    # so that it's much smaller than the SN_A ejecta inner radius

    # Use a density that gives r_outer ~ 1e10 cm for 1 Msun
    # M = (4/3) * pi * r^3 * rho => rho = M / ((4/3) * pi * r^3)
    r_outer = (3 * M_shell / (4 * np.pi * 1e7)) ** (1/3)  # rho ~ 1e7 g/cc
    r_outer = max(r_outer, 1e9)  # At least 1e9 cm
    r_inner = max(r_inner_cm, 1e7)

    print(f"  H/He shell: M={shell_mass_msun} Msun, v={shell_velocity_kms} km/s")
    print(f"    r_inner={r_inner:.3e} cm, r_outer={r_outer:.3e} cm")
    print(f"    X_H={xh}, Y_He={yhe}")

    # Log-spaced zones
    r = np.logspace(np.log10(r_inner), np.log10(r_outer), nzon_shell)
    r_lo = np.concatenate([[r_inner * 0.5], r[:-1]])
    dm = (4/3) * np.pi * (r**3 - r_lo**3)
    # Uniform density
    rho = M_shell / ((4/3) * np.pi * (r_outer**3 - r_inner**3))
    rho_arr = np.full(nzon_shell, rho)

    # Homologous velocity: v = v_shell * r / r_outer
    v_arr = v_shell * r / r_outer

    # Isothermal temperature
    T_arr = np.full(nzon_shell, T_shell_K)

    # Composition
    fH = np.full(nzon_shell, xh)
    fHe = np.full(nzon_shell, yhe)
    fC = np.zeros(nzon_shell)
    fO = np.zeros(nzon_shell)
    fSi = np.zeros(nzon_shell)
    fCa = np.zeros(nzon_shell)
    fFe = np.zeros(nzon_shell)
    fNi = np.zeros(nzon_shell)

    return {
        'r': r, 'rho': rho_arr, 'v': v_arr, 'T': T_arr,
        'fH': fH, 'fHe': fHe, 'fC': fC, 'fO': fO,
        'fSi': fSi, 'fCa': fCa, 'fFe': fFe, 'fNi': fNi,
        'r_inner': r_inner, 'r_outer': r_outer,
    }


def combine_shell_and_csm(shell, csm, nzon_total=500):
    """Combine the H/He shell (inner) with the SN_A ejecta (outer).

    The shell occupies the inner zones, the CSM occupies the outer zones.
    Zone allocation is proportional to mass but with a minimum for each component.
    """
    M_shell = (4/3) * np.pi * (shell['r_outer']**3 - shell['r_inner']**3) * shell['rho'][0]
    # CSM mass: sort by radius (inner->outer) for correct mass calculation
    csm_order = np.argsort(csm['r'])
    csm_r_sorted = csm['r'][csm_order]
    csm_rho_sorted = csm['rho'][csm_order]
    csm_r_lo = np.concatenate([[csm_r_sorted[0]*0.9], csm_r_sorted[:-1]])
    csm_dm = (4/3) * np.pi * (csm_r_sorted**3 - csm_r_lo**3) * csm_rho_sorted
    M_csm = np.sum(np.maximum(csm_dm, 0))
    M_total = M_shell + M_csm

    # Allocate zones: at least 30% to each component, rest proportional to mass
    nzon_shell = max(50, int(nzon_total * max(0.3, M_shell / M_total)))
    nzon_csm = nzon_total - nzon_shell
    print(f"  Zone allocation: {nzon_shell} shell + {nzon_csm} CSM = {nzon_total} total")
    print(f"  Mass: shell={M_shell/MSUN:.3f} Msun, CSM={M_csm/MSUN:.3f} Msun, total={M_total/MSUN:.3f} Msun")

    # Rezone the shell
    shell_r = np.logspace(np.log10(shell['r_inner']), np.log10(shell['r_outer']), nzon_shell)
    shell_r_lo = np.concatenate([[shell['r_inner'] * 0.5], shell_r[:-1]])
    shell_rho = np.full(nzon_shell, M_shell / ((4/3) * np.pi * (shell['r_outer']**3 - shell['r_inner']**3)))
    shell_v = shell['v'][0] * shell_r / shell['r_outer']  # homologous
    shell_T = np.full(nzon_shell, shell['T'][0])

    # Rezone the CSM (log-spaced, preserving structure)
    csm_r_orig = csm['r']  # outer->inner (decreasing)
    # Reverse to inner->outer for the CSM
    order = np.argsort(csm_r_orig)
    csm_r_sorted = csm_r_orig[order]
    csm_rho_sorted = csm['rho'][order]
    csm_v_sorted = csm['v'][order]
    csm_T_sorted = csm['T'][order]

    # CSM starts just outside the shell
    csm_r_inner = shell['r_outer'] * 1.01  # small gap
    csm_r_outer = csm_r_sorted[-1]

    # Log-spaced CSM grid
    csm_r_new = np.logspace(np.log10(csm_r_inner), np.log10(csm_r_outer), nzon_csm)

    # Interpolate CSM profiles
    from scipy.interpolate import interp1d
    log_r_orig = np.log10(csm_r_sorted)

    def log_interp(y, name):
        log_y = np.log10(np.maximum(y, 1e-300))
        ip = interp1d(log_r_orig, log_y, kind='linear',
                      bounds_error=False, fill_value=(log_y[0], log_y[-1]))
        return 10.0 ** ip(np.log10(csm_r_new))

    def lin_interp(y, name):
        ip = interp1d(log_r_orig, y, kind='linear',
                      bounds_error=False, fill_value=(y[0], y[-1]))
        return ip(np.log10(csm_r_new))

    csm_rho_new = log_interp(csm_rho_sorted, 'rho')
    csm_v_new = lin_interp(csm_v_sorted, 'v')
    csm_T_new = np.maximum(log_interp(csm_T_sorted, 'T'), 1000.0)

    # Interpolate CSM composition
    csm_comp = {}
    for key in ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe', 'fNi']:
        if key in csm:
            csm_comp[key] = lin_interp(csm[key][order], key)
        else:
            csm_comp[key] = np.zeros(nzon_csm)

    # Combine: shell (inner) + CSM (outer)
    r_combined = np.concatenate([shell_r, csm_r_new])
    rho_combined = np.concatenate([shell_rho, csm_rho_new])
    v_combined = np.concatenate([shell_v, csm_v_new])
    T_combined = np.concatenate([shell_T, csm_T_new])

    comp_combined = {}
    for key in ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe', 'fNi']:
        shell_val = shell.get(key, np.zeros(nzon_shell))
        # For shell, use uniform composition
        if key == 'fH':
            shell_val = np.full(nzon_shell, shell.get('fH', [0.7])[0])
        elif key == 'fHe':
            shell_val = np.full(nzon_shell, shell.get('fHe', [0.3])[0])
        else:
            shell_val = np.zeros(nzon_shell)
        comp_combined[key] = np.concatenate([shell_val, csm_comp[key]])

    # Renormalize composition
    comp_sum = sum(comp_combined.values())
    comp_sum = np.maximum(comp_sum, 1e-10)
    for key in comp_combined:
        comp_combined[key] /= comp_sum

    return {
        'r': r_combined,
        'rho': rho_combined,
        'v': v_combined,
        'T': T_combined,
        **comp_combined,
    }


def write_sn_b_hyd(fname, g, mass_cut_msun=0.0):
    """Write STELLA .hyd file for SN_B."""
    r = g['r']; rho = g['rho']; T = g['T']; v = g['v']
    nzon = len(r)

    r_lo = np.concatenate([[r[0] * 0.5], r[:-1]])
    dm = (4/3) * np.pi * (r**3 - r_lo**3) * rho
    dm = np.maximum(dm, 0.0)
    M_r = np.cumsum(dm)

    rcen = r[0] * 0.5
    rho_cen = rho[0]

    with open(fname, 'w') as f:
        f.write(f"{0.0:12.3E}{nzon:6d}{mass_cut_msun:13.5E}{rcen:13.5E}{rho_cen:13.5E}\n")
        for km in range(nzon):
            Mr = M_r[km] / MSUN + mass_cut_msun
            f.write(f"{km+1:6d}{r[km]:13.5E}{rho[km]:13.5E}{T[km]:13.5E}"
                    f"{v[km]:13.5E}{Mr:13.5E}\n")
    print(f"Written {fname}  ({nzon} zones, M_total={M_r[-1]/MSUN:.4f} Msun)")


def write_sn_b_abn(fname, g):
    """Write STELLA .abn file for SN_B."""
    n = len(g['r'])
    with open(fname, 'w') as f:
        for km in range(n):
            f.write(f"{km+1:4d}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fH'][km]:10.3E}{g['fHe'][km]:10.3E}{g['fC'][km]:10.3E}"
                    f"{0.0:10.3E}{g['fO'][km]:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}{0.0:10.3E}"
                    f"{g['fSi'][km]:10.3E}"
                    f"{0.0:10.3E}{0.0:10.3E}{g['fCa'][km]:10.3E}{g['fFe'][km]:10.3E}"
                    f"{0.0:10.3E}{g['fNi'][km]:10.3E}\n")
    print(f"Written {fname}  ({n} zones)")


def write_sn_b_xni(fname, nzon):
    """Write a zero .xni file (no radioactive decay in SN_B)."""
    n_floats = 2 + 2 * nzon
    floats = [0.0] * n_floats
    floats[0] = 5.605194e-42
    floats[1] = -1.417409e-34
    data = struct.pack('f' * n_floats, *floats)
    with open(fname, 'wb') as f:
        f.write(struct.pack('i', len(data)))
        f.write(data)
        f.write(struct.pack('i', len(data)))
    print(f"Written {fname}  ({len(data)} bytes, zero Ni56)")


def write_sn_b_dat(fname, model_name, target_days, snapshot_day,
                   explosion_energy=0.0, explosion_mass=0.0, explosion_time=1.0,
                   nstmax=2000000, eps=0.003, hmax=3600.0, meth=3,
                   aq=4.0, bq=1.0, iout=1, snapshot_interval=5.0):
    """Write STELLA .dat control file for SN_B.

    SN_B uses EKO=0 (preserve velocities) since the shell already has
    its launch velocity encoded in the .hyd file.
    """
    t_start = 0.0
    t_end = target_days

    # Output epochs
    KOMAX = 80
    # Dense sampling around the snapshot day, plus early epochs
    early = [0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]
    if snapshot_day not in early:
        early.append(snapshot_day)
    late = np.linspace(max(35.0, snapshot_day), target_days,
                       max(5, int((target_days - 35) / snapshot_interval))).tolist()
    epochs = sorted(set(early + late))[:KOMAX]

    with open(fname, 'w') as f:
        f.write(f"   EPS   (STIFF ACCURACY)   HMIN (MIN.STEP(SEC)) HMAX (MAX.STEP(SEC))\n")
        f.write(f"   {eps:.4g}                    1e-14               {hmax:.4g}\n")
        f.write(f"   METH  ADAMS(1)-GEAR(2)-GEAR+BGH(3)   JSTART  (0-START_ELSE_ORDER) MAXORD KNadap\n")
        f.write(f"    {meth}                        0                            4        -4\n")
        f.write(f"   abs(NSTA)>abs(NSTB)=>reads last   TcurA     TcurB\n")
        f.write(f"   -9000     -1                      {t_start:.4f}   {t_end:.4f}\n")
        f.write(f"   NSTMAX  NDebug  NOUT    IOUT (0 -full <0 print zones)   MBATCH\n")
        # NOUT must be < 200 to avoid per-step PRIBAL calls (stradsep5tt.f:640)
        nout_val = len(epochs)
        f.write(f"   {nstmax}   100000  {nout_val}      {iout}                             900\n")
        f.write(f"   AMNI (CONTAM. MASS)     XMNI (56NI MASS)\n")
        f.write(f"   0.0d0                  0.0d0\n")
        f.write(f"   AMHT(Solar)  EBurst(1e50) tBurst(s)  tstart -- Heated Core, Energy & time\n")
        f.write(f"   {explosion_mass:.4g}d0      {explosion_energy:.4g}d0      {explosion_time:.4g}d0      0.0d0\n")
        f.write(f"   EKO (kin.energy,1e50) mass fract.tri. u profile  us (out +1)\n")
        # EKO=0: preserve .hyd velocities (shell has v=5000 km/s, CSM has homologous v)
        f.write(f"   0.0d0                      0.d0  9.d-1  1.d0          +1.\n")
        f.write(f"   THRMAT   CRAP     CONV     EDTM     CHNCND  Givdtl\n")
        f.write(f"   1.D-30   1.D0     F        T        T       T\n")
        f.write(f"   FLOOR(R) FLOOR(V) FLOOR(T) FLOOR(RADIAT)    FLOOR(UCONV)\n")
        f.write(f"   1.D-12   1.D-03   1.D-03D+00   1.D-02           1.D-04\n")
        f.write(f"   Wacc(R)  Wacc(V)  Wacc(T)  Wacc(RADIAT)\n")
        f.write(f"   1d-0    1d-0    1.0D-00   1D-00\n")
        f.write(f"   FitTau   TAUTOL   Rvis   AQ    BQ       DRT     NRT     SCAT\n")
        f.write(f"   15.0D0   1.1D0     1.D0   {aq:.1f}D0  {bq:.1f}D+00   1.      1       T\n")
        f.write(f"    NTO then    DATA TO\n")
        f.write(f"    {len(epochs)}\n")
        for ep in epochs:
            f.write(f"     {ep:.4f}\n")
    print(f"Written {fname}  ({len(epochs)} output epochs, target={target_days} days)")


def write_sn_b_input_profile(fname, g):
    """Write 12-column input_profile for convert_profile_to_stella.py."""
    r = g['r']; rho = g['rho']; T = g['T']

    fH = g['fH']; fHe = g['fHe']; fC = g['fC']; fO = g['fO']
    fSi = g['fSi']; fFe = g['fFe']; fCa = g.get('fCa', np.zeros_like(fH))
    fNi = g.get('fNi', np.zeros_like(fH))

    mu = 1.0 / (fH/1.0 + fHe/4.0 + fO/16.0 + fSi/28.0 + fCa/40.0 + fFe/56.0 + fNi/58.0)
    mu = np.maximum(mu, 0.5)
    p_gas = rho * K_B * T / (mu * M_H)
    p_rad = A_RAD * T**4 / 3.0
    p = p_gas + p_rad

    fO_profile = fO + fC  # Merge C into O

    header = ("# input_profile generated by prepare_sn_b.py\n"
              "# columns: r_cm rho_cgs p_cgs u_cms T_K fx_H fx_He fx_O fx_Si fx_Ca fx_Fe fx_Ni\n"
              "# Inner zones: H/He shell; Outer zones: SN_A ejecta (349-day CSM)\n")
    with open(fname, 'w') as f:
        f.write(header)
        for i in range(len(r)):
            f.write(f"{r[i]:.8e} {rho[i]:.8e} {p[i]:.8e} {g['v'][i]:.8e} "
                    f"{T[i]:.8e} "
                    f"{fH[i]:.8e} {fHe[i]:.8e} {fO_profile[i]:.8e} "
                    f"{fSi[i]:.8e} {fCa[i]:.8e} {fFe[i]:.8e} {fNi[i]:.8e}\n")
    print(f"Written {fname}  ({len(r)} zones)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sn_a_dir', help='SN_A model directory (e.g. SN_A)')
    parser.add_argument('sn_b_dir', help='SN_B output directory (e.g. SN_B)')
    parser.add_argument('--shell-mass', type=float, default=1.0,
                        help='H/He shell mass in Msun (default: 1.0)')
    parser.add_argument('--shell-velocity', type=float, default=5000,
                        help='Shell launch velocity in km/s (default: 5000)')
    parser.add_argument('--shell-xh', type=float, default=0.7,
                        help='H mass fraction in shell (default: 0.7)')
    parser.add_argument('--shell-yhe', type=float, default=0.3,
                        help='He mass fraction in shell (default: 0.3)')
    parser.add_argument('--target-days', type=float, default=50.0,
                        help='STELLA simulation duration in days (default: 50)')
    parser.add_argument('--snapshot-day', type=float, default=40.0,
                        help='CMFGEN snapshot day after shell launch (default: 40)')
    parser.add_argument('--nzon', type=int, default=500,
                        help='Total number of STELLA zones (default: 500)')
    parser.add_argument('--sn-a-age', type=float, default=349.0,
                        help='Age of SN_A ejecta at shell launch in days (default: 349)')
    parser.add_argument('--explosion-energy', type=float, default=0.0,
                        help='EBurst in 1e50 erg (0=no bomb, shell velocity from .hyd)')
    parser.add_argument('--nstmax', type=int, default=2000000)
    parser.add_argument('--aq', type=float, default=4.0,
                        help='Artificial viscosity (default: 4.0 for shell-CSM interaction)')
    args = parser.parse_args()

    sn_a_dir = Path(args.sn_a_dir)
    sn_b_dir = Path(args.sn_b_dir)
    model_name = sn_b_dir.name

    # Create output directories
    stella_dir = sn_b_dir / 'stella'
    modmake_dir = stella_dir / 'modmake'
    plots_dir = sn_b_dir / 'plots'
    modmake_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Preparing SN_B initial model ===")
    print(f"  SN_A directory: {sn_a_dir}")
    print(f"  SN_B directory: {sn_b_dir}")
    print(f"  Shell: M={args.shell_mass} Msun, v={args.shell_velocity} km/s, X={args.shell_xh}, Y={args.shell_yhe}")
    print(f"  SN_A age at launch: {args.sn_a_age} days")
    print(f"  STELLA target: {args.target_days} days, CMFGEN snapshot: {args.snapshot_day} days")

    # Step 1: Read SN_A ejecta
    print(f"\n--- Step 1: Reading SN_A ejecta at {args.sn_a_age} days ---")
    csm = read_sn_a_ejecta(sn_a_dir, sn_a_age=args.sn_a_age)

    # Extract composition from CSM
    if csm.get('composition') is not None:
        comp = csm['composition']
        csm['fH'] = comp[:, 0]
        csm['fHe'] = comp[:, 1]
        csm['fC'] = comp[:, 2]
        csm['fO'] = comp[:, 4]
        csm['fSi'] = comp[:, 10]
        csm['fCa'] = comp[:, 12]
        csm['fFe'] = comp[:, 13]
        csm['fNi'] = comp[:, 14]
    else:
        # No composition data — set defaults
        for key in ['fH', 'fHe', 'fC', 'fO', 'fSi', 'fCa', 'fFe', 'fNi']:
            csm[key] = np.zeros(len(csm['r']))
        csm['fHe'] = np.ones(len(csm['r']))  # Default to He

    print(f"  CSM: {len(csm['r'])} zones, r=[{csm['r'].min():.3e}, {csm['r'].max():.3e}] cm")
    print(f"  CSM: v=[{csm['v'].min():.1f}, {csm['v'].max():.1f}] km/s")
    print(f"  CSM: rho=[{csm['rho'].min():.3e}, {csm['rho'].max():.3e}] g/cc")

    # Convert CSM velocities from km/s to cm/s for .hyd file
    csm['v'] = csm['v'] * 1e5  # km/s -> cm/s

    # Step 2: Create H/He shell
    print(f"\n--- Step 2: Creating H/He shell ---")
    shell = create_hhe_shell(args.shell_mass, args.shell_velocity,
                             args.shell_xh, args.shell_yhe)

    # Step 3: Combine shell + CSM
    print(f"\n--- Step 3: Combining shell + CSM ---")
    combined = combine_shell_and_csm(shell, csm, nzon_total=args.nzon)

    print(f"  Combined: {len(combined['r'])} zones")
    print(f"  r: [{combined['r'].min():.3e}, {combined['r'].max():.3e}] cm")
    print(f"  v: [{combined['v'].min():.1f}, {combined['v'].max():.1f}] cm/s "
          f"([{combined['v'].min()/1e5:.1f}, {combined['v'].max()/1e5:.1f}] km/s)")
    print(f"  rho: [{combined['rho'].min():.3e}, {combined['rho'].max():.3e}] g/cc")

    # Step 4: Write STELLA input files
    print(f"\n--- Step 4: Writing STELLA input files ---")
    hyd_path = modmake_dir / f'{model_name}.hyd'
    abn_path = modmake_dir / f'{model_name}.abn'
    xni_path = modmake_dir / f'{model_name}.xni'
    dat_path = modmake_dir / f'{model_name}.dat'
    ip_path = modmake_dir / 'input_profile'

    write_sn_b_hyd(str(hyd_path), combined)
    write_sn_b_abn(str(abn_path), combined)
    write_sn_b_xni(str(xni_path), args.nzon)
    write_sn_b_dat(str(dat_path), model_name, args.target_days, args.snapshot_day,
                   explosion_energy=args.explosion_energy,
                   nstmax=args.nstmax, aq=args.aq)
    write_sn_b_input_profile(str(ip_path), combined)

    # Also copy input_profile to stella dir
    import shutil
    shutil.copy2(str(ip_path), str(stella_dir / 'input_profile'))

    # Step 5: Generate plots
    print(f"\n--- Step 5: Generating plots ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Structure plot
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'{model_name} - Initial Model (Shell + SN_A CSM)', fontsize=14, fontweight='bold')
        r_rsun = combined['r'] / 6.96e10

        axes[0, 0].loglog(r_rsun, combined['rho'], 'b-', linewidth=0.8)
        axes[0, 0].set_xlabel('Radius (R_sun)'); axes[0, 0].set_ylabel('Density (g/cm³)')
        axes[0, 0].set_title('Density Profile'); axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].semilogx(r_rsun, combined['v'] / 1e5, 'r-', linewidth=0.8)
        axes[0, 1].set_xlabel('Radius (R_sun)'); axes[0, 1].set_ylabel('Velocity (km/s)')
        axes[0, 1].set_title('Velocity Profile'); axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].axhline(args.shell_velocity, color='r', linestyle='--', alpha=0.5, label=f'Shell v={args.shell_velocity} km/s')
        axes[0, 1].legend()

        axes[1, 0].loglog(r_rsun, combined['T'], 'g-', linewidth=0.8)
        axes[1, 0].set_xlabel('Radius (R_sun)'); axes[1, 0].set_ylabel('Temperature (K)')
        axes[1, 0].set_title('Temperature Profile'); axes[1, 0].grid(True, alpha=0.3)

        # Composition
        colors = {'H': 'blue', 'He': 'green', 'C': 'orange', 'O': 'red',
                  'Si': 'purple', 'Fe': 'brown', 'Ni': 'gray', 'Ca': 'pink'}
        for elem in ['H', 'He', 'C', 'O', 'Si', 'Fe']:
            farr = combined[f'f{elem}']
            if np.any(farr > 1e-10):
                axes[1, 1].semilogx(r_rsun, farr, label=elem, color=colors[elem], linewidth=1.0)
        axes[1, 1].set_xlabel('Radius (R_sun)'); axes[1, 1].set_ylabel('Mass Fraction')
        axes[1, 1].set_title('Composition'); axes[1, 1].legend(fontsize=8)
        axes[1, 1].set_ylim(1e-4, 1.5); axes[1, 1].set_yscale('log')
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(str(plots_dir / 'sn_b_input_structure.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Structure plot: {plots_dir / 'sn_b_input_structure.png'}")

        # Mass fraction distribution
        r_lo = np.concatenate([[combined['r'][0]*0.5], combined['r'][:-1]])
        dm = (4/3) * np.pi * (combined['r']**3 - r_lo**3) * combined['rho']
        m_enc = np.cumsum(dm) / MSUN

        fig, ax = plt.subplots(figsize=(12, 6))
        fig.suptitle(f'{model_name} - Mass Fraction Distribution', fontsize=14, fontweight='bold')
        for elem in ['H', 'He', 'C', 'O', 'Si', 'Fe']:
            farr = combined[f'f{elem}']
            if np.any(farr > 1e-10):
                ax.plot(m_enc, farr, label=elem, color=colors[elem], linewidth=1.0)
        ax.set_xlabel('Enclosed Mass (Msun)'); ax.set_ylabel('Mass Fraction')
        ax.set_ylim(1e-4, 1.5); ax.set_yscale('log')
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(str(plots_dir / 'sn_b_input_mass_fractions.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Mass fraction plot: {plots_dir / 'sn_b_input_mass_fractions.png'}")

    except Exception as ex:
        print(f"  WARNING: Plot generation failed: {ex}")

    print(f"\n=== Done. SN_B files written to {modmake_dir}/ ===")
    print(f"  {model_name}.hyd, {model_name}.abn, {model_name}.xni, {model_name}.dat, input_profile")
    print(f"  Plots in {plots_dir}/")


if __name__ == '__main__':
    main()
