#!/usr/bin/env python3
"""Plot STELLA structure overview from .swd file (similar to [DEPRECATED] HERACLES structure_overview.png).

HERACLES is deprecated; the project now uses the STELLA -> CMFGEN pipeline.

STELLA .swd file format (ASCII):
  t_day   km   log10(M_ext_Msun)   log10(R_cm)   V_kms
  log10(T_gas_K)   log10(T_rad_K)   log10(rho_gcc)
  log10(P_cgs)   log10(Qv)   log10(eng)   Flum_1e40   WRKX

Usage:
    python plot_stella_structure.py SWD_FILE [OPTIONS]
      --abn-file PATH    STELLA .abn file for mass fractions
      --epoch T          Epoch to plot in days (default: last)
      --output PATH      Output plot path (default: structure_overview.png)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_swd(swd_path: str) -> dict:
    """Read STELLA .swd file and return dict of epochs."""
    data = []
    current_epoch = None
    epoch_data = []
    
    with open(swd_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            
            cols = stripped.split()
            if len(cols) < 13:
                continue
            
            # First column is time for header lines (km=1 or km=Nzon)
            t_val = float(cols[0])
            km = int(cols[1])
            
            if km == 1 and t_val > 0:
                # New epoch header
                if epoch_data:
                    data.append((current_epoch, np.array(epoch_data)))
                current_epoch = t_val
                epoch_data = []
            
            epoch_data.append([float(x) for x in cols])
    
    # Don't forget last epoch
    if epoch_data:
        data.append((current_epoch, np.array(epoch_data)))
    
    return {t: arr for t, arr in data}


def read_abn(abn_path: str) -> dict:
    """Read STELLA .abn file for mass fractions.
    
    Format: km  dum  dum  dum  H  He  C  N  O  Ne  Na  Mg  Al  Si  S  Ar  Ca  FePeak  Ni58  Ni56
    """
    data = []
    with open(abn_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            cols = [float(x) for x in stripped.split()]
            if len(cols) >= 20:
                data.append(cols)
    
    arr = np.array(data)
    return {
        'km': arr[:, 0].astype(int),
        'H': arr[:, 4],
        'He': arr[:, 5],
        'C': arr[:, 6],
        'N': arr[:, 7],
        'O': arr[:, 8],
        'Ne': arr[:, 9],
        'Na': arr[:, 10],
        'Mg': arr[:, 11],
        'Al': arr[:, 12],
        'Si': arr[:, 13],
        'S': arr[:, 14],
        'Ar': arr[:, 15],
        'Ca': arr[:, 16],
        'Fe': arr[:, 17],
        'Ni58': arr[:, 18],
        'Ni56': arr[:, 19],
    }


def plot_structure_overview(swd_data: np.ndarray, abn_data: dict, epoch_days: float,
                            output_path: str, model_name: str):
    """Create 6-panel structure overview plot."""
    # Extract data from swd (columns 0-12)
    # t_day, km, log10(M), log10(R), V(1e8 cm/s), log10(Tg), log10(Tr), log10(rho/1e-6), log10(P), log10(Qv), log10(eng), Flum(1e40), WRKX
    km = swd_data[:, 1].astype(int)
    R_cm = 10**swd_data[:, 3]
    V_kms = swd_data[:, 4] * 1000.0  # Convert from 1e8 cm/s to km/s
    T_gas = 10**swd_data[:, 5]
    T_rad = 10**swd_data[:, 6]
    rho = (10**swd_data[:, 7]) * 1e-6  # Convert from log10(rho/1e-6) to g/cm^3
    P = 10**swd_data[:, 8]
    
    # Extract mass fractions from abn (match by km index)
    fx_H = np.array([abn_data['H'][i-1] for i in km])
    fx_He = np.array([abn_data['He'][i-1] for i in km])
    fx_C = np.array([abn_data['C'][i-1] for i in km])
    fx_O = np.array([abn_data['O'][i-1] for i in km])
    fx_Si = np.array([abn_data['Si'][i-1] for i in km])
    fx_Fe = np.array([abn_data['Fe'][i-1] for i in km])
    
    # Calculate enclosed mass
    dr = np.diff(R_cm, prepend=R_cm[0])
    r_mid = 0.5 * (R_cm[1:] + R_cm[:-1])
    r_mid = np.concatenate([[R_cm[0]], r_mid])
    dm = 4 * np.pi * r_mid**2 * rho * dr
    M_enc = np.cumsum(dm) / 1.989e33  # Convert to Msun
    
    # Radius in 10^15 cm for plotting
    r_plt = R_cm / 1e15
    
    # Create 6-panel figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    # 1. Density (log-log)
    ax1 = axes[0]
    ax1.loglog(r_plt, rho, 'b-', linewidth=2)
    ax1.set_xlabel('Radius (10$^{15}$ cm)')
    ax1.set_ylabel('Density (g cm$^{-3}$)')
    ax1.set_title('Density')
    ax1.grid(True, alpha=0.3)
    
    # 2. Velocity (log-log, 10^4 km/s)
    ax2 = axes[1]
    ax2.loglog(r_plt, V_kms/1e4, 'r-', linewidth=2)
    ax2.set_xlabel('Radius (10$^{15}$ cm)')
    ax2.set_ylabel('Velocity ($10^4$ km s$^{-1}$)')
    ax2.set_title('Velocity')
    ax2.grid(True, alpha=0.3)
    
    # 3. Pressure (log-log)
    ax3 = axes[2]
    ax3.loglog(r_plt, P, 'g-', linewidth=2)
    ax3.set_xlabel('Radius (10$^{15}$ cm)')
    ax3.set_ylabel('Pressure (erg cm$^{-3}$)')
    ax3.set_title('Pressure')
    ax3.grid(True, alpha=0.3)
    
    # 4. Composition (log-x)
    ax4 = axes[3]
    ax4.semilogx(r_plt, fx_H, label='H', linewidth=2)
    ax4.semilogx(r_plt, fx_He, label='He', linewidth=2)
    ax4.semilogx(r_plt, fx_O, label='O', linewidth=2)
    ax4.semilogx(r_plt, fx_Si, label='Si', linewidth=2)
    ax4.semilogx(r_plt, fx_Fe, label='Fe', linewidth=2)
    ax4.set_xlabel('Radius (10$^{15}$ cm)')
    ax4.set_ylabel('Mass Fraction')
    ax4.set_title('Composition')
    ax4.legend(fontsize=8, ncol=2)
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, 1.1)
    
    # 5. Enclosed mass (linear y-axis)
    ax5 = axes[4]
    ax5.plot(r_plt, M_enc, 'k-', linewidth=2)
    ax5.set_xlabel('Radius (10$^{15}$ cm)')
    ax5.set_xscale('log')
    ax5.set_yscale('linear')
    ax5.set_ylabel(r'Enclosed Mass (M$_\odot$)')
    ax5.set_title(r'Enclosed Mass (Total: ' + f'{M_enc[-1]:.2f}' + r' M$_\odot$)')
    ax5.grid(True, alpha=0.3)
    
    # 6. Temperature (log x, linear y)
    ax6 = axes[5]
    ax6.semilogx(r_plt, T_gas, 'g-', linewidth=2, label='T_gas')
    ax6.semilogx(r_plt, T_rad, 'm--', linewidth=1.5, label='T_rad')
    ax6.set_xlabel('Radius (10$^{15}$ cm)')
    ax6.set_ylabel('Temperature (K)')
    ax6.set_title('Temperature')
    ax6.set_yscale('linear')
    ax6.legend(fontsize=8)
    ax6.grid(True, alpha=0.3)
    
    plt.suptitle(f'STELLA Structure Overview — {model_name} @ {epoch_days:.2f} days',
                 fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Structure overview saved to {output_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('swd_file', help='STELLA .swd file')
    parser.add_argument('--abn-file', required=True, help='STELLA .abn file for mass fractions')
    parser.add_argument('--epoch', type=float, default=None,
                        help='Epoch to plot in days (default: last)')
    parser.add_argument('--output', default=None,
                        help='Output plot path (default: structure_overview.png)')
    args = parser.parse_args()
    
    swd_path = Path(args.swd_file)
    if not swd_path.exists():
        print(f"ERROR: swd file not found: {swd_path}", file=sys.stderr)
        sys.exit(1)
    
    abn_path = Path(args.abn_file)
    if not abn_path.exists():
        print(f"ERROR: abn file not found: {abn_path}", file=sys.stderr)
        sys.exit(1)
    
    model_name = swd_path.stem
    output_path = args.output or str(swd_path.parent / 'structure_overview.png')
    
    print(f"Reading {swd_path}...")
    swd_epochs = read_swd(str(swd_path))
    print(f"  Found {len(swd_epochs)} epochs: {list(swd_epochs.keys())}")
    
    valid_epochs = [t for t in swd_epochs.keys() if t is not None]
    if args.epoch is None:
        epoch_days = max(valid_epochs)
    else:
        epoch_days = min(valid_epochs, key=lambda t: abs(t - args.epoch))
    
    print(f"  Plotting epoch: {epoch_days:.3f} days")
    swd_data = swd_epochs[epoch_days]
    
    print(f"Reading {abn_path}...")
    abn_data = read_abn(str(abn_path))
    print(f"  {len(abn_data['km'])} zones")
    
    plot_structure_overview(swd_data, abn_data, epoch_days, output_path, model_name)


if __name__ == '__main__':
    main()
