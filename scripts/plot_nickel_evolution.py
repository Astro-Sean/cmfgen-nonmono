#!/usr/bin/env python3
"""
Plot Ni56 -> Co56 -> Fe56 evolution for a STELLA run.

STELLA does not advect/decay composition fields: XNifor (initial Ni56 mass
fraction per zone) is static and enters only the radioactive heating rate.
The Ni56/Co56/Fe56 evolution is therefore reconstructed analytically on the
Lagrangian zones:

    Ni56(t) = X0 * exp(-t/tau_Ni)
    Co56(t) = X0 * (tau_Co/(tau_Co-tau_Ni)) * (exp(-t/tau_Co) - exp(-t/tau_Ni))
    Fe56(t) = X0 - Ni56 - Co56            (Bateman chain, tau_Ni=8.8d, tau_Co=111.3d)

Spatial evolution is taken from the .swd: each zone's radius at every epoch.
The decay age is PROPER_T ~= t_ob - t_ob(first epoch), i.e. the run's internal
clock (t = 0 corresponds to the .xni initial distribution).

Usage:
    python plot_nickel_evolution.py MODEL_DIR [OPTIONS]

Options:
    --xni-file PATH          Path to .xni file (default: stella/eve/MODEL.xni)
    --abn-dir PATH           Directory containing input .abn (initial Ni56 col check)
    --swd-file PATH          Path to .swd file for radius data
    --n-epochs N             Number of epochs to plot (default: 10)
    --output PATH            Output plot path (default: nickel_evolution.png)
"""

import struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from pathlib import Path
import sys
import argparse


# Decay e-folding times (STELLA: volenpumnoint.f; also used by
# convert_stella_swd_to_cmfgen.py for the SN_HYDRO_DATA isotope blocks)
TAU_NI = 8.8    # days, Ni56 -> Co56
TAU_CO = 111.3  # days, Co56 -> Fe56


def read_xni(xni_path: str) -> np.ndarray:
    """Read a STELLA .xni file -> Ni56 mass fraction per zone.

    Two on-disk formats exist:
      * Fortran unformatted single record of nzon REAL*8 values
        (written by generate_xni_file.py and eve2): 4-byte record markers
        bracket nzon doubles.
      * Legacy float32 layout: 2 header floats then nzon (value, pad) pairs.

    Returns:
        ni56 mass fraction array, inner -> outer zone order.
    """
    data = Path(xni_path).read_bytes()

    # Fortran unformatted REAL*8 record: [rec_len][rec_len bytes][rec_len]
    if len(data) >= 8:
        rec_len = struct.unpack('i', data[:4])[0]
        tail = struct.unpack('i', data[-4:])[0]
        if rec_len > 0 and rec_len % 8 == 0 and 4 + rec_len + 4 == len(data) \
                and tail == rec_len:
            n = rec_len // 8
            return np.array(struct.unpack(f'{n}d', data[4:4 + rec_len]))

    # Legacy float32: 2 header floats + nzon pairs
    n_floats = len(data) // 4
    floats = struct.unpack('f' * n_floats, data)
    nzon = (n_floats - 2) // 2
    return np.array([floats[2 + 2 * z] for z in nzon])


def read_swd_epochs(swd_path: str):
    """Parse a STELLA .swd into per-epoch zone data.

    Format: an epoch starts on a line whose first column is the epoch time
    (observer time t_ob = PROPER_T - t_retard); its second column is 1.
    Following rows have 0 in column 0 and the zone index in column 1.
    Relevant columns: [0]=t or 0, [1]=zone, [2]=log10(M_ext/Msun),
    [3]=log10(R/cm), [4]=V/1e8 cm/s.

    Returns list of (t_ob_days, radius_cm[nzones], velocity_cms[nzones]).
    """
    epochs = []
    t_cur = None
    rad, vel = [], []
    with open(swd_path) as f:
        for line in f:
            cols = line.split()
            if len(cols) < 13:
                continue
            try:
                c0 = float(cols[0])
                zone = int(cols[1])
            except ValueError:
                continue
            if c0 != 0.0:
                if t_cur is not None and rad:
                    epochs.append((t_cur, np.array(rad), np.array(vel)))
                t_cur = c0
                rad, vel = [10.0 ** float(cols[3])], [float(cols[4]) * 1e8]
            else:
                rad.append(10.0 ** float(cols[3]))
                vel.append(float(cols[4]) * 1e8)
    if t_cur is not None and rad:
        epochs.append((t_cur, np.array(rad), np.array(vel)))
    return epochs


def decay_chain(X0, age_days):
    """Bateman solution for Ni56 -> Co56 -> Fe56 given initial fraction X0."""
    t = np.maximum(age_days, 0.0)
    ni = X0 * np.exp(-t / TAU_NI)
    co = X0 * (TAU_CO / (TAU_CO - TAU_NI)) * (np.exp(-t / TAU_CO) - np.exp(-t / TAU_NI))
    fe = X0 - ni - co
    return ni, co, fe


def plot_nickel_evolution(xni_path, abn_dir, swd_path, n_epochs, output_path,
                          model_name):
    """Ni56 and Co56 mass fraction vs radius at several epochs.

    The .xni initial distribution is placed on the Lagrangian zones and decayed
    analytically; each epoch's zone radii come from the .swd so the curves show
    the distribution moving outward while decaying.
    """
    print(f"Reading initial Ni56 from {xni_path}")
    ni56_0 = read_xni(xni_path)
    nzon = len(ni56_0)
    print(f"  {nzon} zones, max X(Ni56) = {ni56_0.max():.4e}, "
          f"{(ni56_0 > 1e-30).sum()} non-zero")

    epochs = read_swd_epochs(swd_path)
    print(f"Found {len(epochs)} epochs in {swd_path}")
    if not epochs:
        print("ERROR: no epochs parsed from .swd")
        return

    # Decay clock: PROPER_T ~= t_ob - t_ob(first).  Zone count may differ
    # between epochs (Nzmod); interpolate onto the .xni zone count when needed.
    t_first = epochs[0][0]

    if len(epochs) > n_epochs:
        pick = np.linspace(0, len(epochs) - 1, n_epochs, dtype=int)
        epochs = [epochs[i] for i in pick]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=True)
    species = [
        ('Ni56', r'$^{56}$Ni'),
        ('Co56', r'$^{56}$Co'),
        ('Fe56', r'$^{56}$Fe'),
    ]

    norm = Normalize(vmin=0.0, vmax=max(t - t_first for t, _, _ in epochs))
    cmap = cm.viridis

    # Initial Ni56 distribution on the first epoch's radius grid
    r0 = epochs[0][1]
    if len(r0) != nzon:
        r0 = np.interp(np.linspace(0, 1, nzon), np.linspace(0, 1, len(r0)), r0)
    for ax in axes:
        ax.plot(r0, ni56_0, 'k--', lw=2, alpha=0.8,
                label=r'initial $^{56}$Ni (.xni)')

    for t_ob, radius, _vel in epochs:
        age = t_ob - t_first
        ni, co, fe = decay_chain(ni56_0, age)
        r = radius
        if len(r) != nzon:
            # resample zone values onto the epoch's zone-fraction axis
            xi = np.linspace(0, 1, nzon)
            xq = np.linspace(0, 1, len(r))
            ni = np.interp(xq, xi, ni)
            co = np.interp(xq, xi, co)
            fe = np.interp(xq, xi, fe)
        color = cmap(norm(age))
        for ax, (key, _label), prof in zip(axes, species, (ni, co, fe)):
            ax.plot(r, np.maximum(prof, 1e-12), color=color, lw=1.8, alpha=0.9)
        axes[0].plot([], [], color=color, lw=1.8,
                     label=f't = {age:.1f} d')

    for ax, (_key, label) in zip(axes, species):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-8, 5e-1)
        ax.set_xlabel('Radius (cm)')
        ax.set_title(label)
        ax.grid(True, alpha=0.3, which='both')
    axes[0].set_ylabel('Mass fraction')
    axes[0].legend(fontsize=8, loc='upper right')

    fig.suptitle(
        f'{model_name}: radioactive decay chain on Lagrangian zones '
        f'(decay age = PROPER_T; $\\tau_{{Ni}}$={TAU_NI} d, '
        f'$\\tau_{{Co}}$={TAU_CO} d)')
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.subplots_adjust(left=0.06, right=0.90, top=0.82, wspace=0.12)
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.67])
    fig.colorbar(sm, cax=cax, label='time since explosion (d)')

    plt.savefig(output_path, dpi=150)
    print(f"Saved plot to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot Ni56/Co56 evolution from STELLA')
    parser.add_argument('model_dir', help='Model directory path')
    parser.add_argument('--xni-file', help='Path to .xni file')
    parser.add_argument('--abn-dir',
                        help='Directory containing .abn files (unused; kept '
                             'for backward-compatible invocation)')
    parser.add_argument('--swd-file', help='Path to .swd file for radius data')
    parser.add_argument('--n-epochs', type=int, default=10,
                        help='Number of epochs to plot (default: 10)')
    parser.add_argument('--output', default='nickel_evolution.png',
                        help='Output plot path (default: nickel_evolution.png)')

    args = parser.parse_args()

    model_dir = Path(args.model_dir)
    model_name = model_dir.name

    xni_path = args.xni_file or str(model_dir / 'stella' / 'eve' / f'{model_name}.xni')
    swd_path = args.swd_file or str(model_dir / 'stella' / f'{model_name}.swd')

    if not Path(xni_path).exists():
        print(f"ERROR: .xni file not found: {xni_path}")
        sys.exit(1)
    if not Path(swd_path).exists():
        print(f"ERROR: .swd file not found: {swd_path}")
        sys.exit(1)

    plot_nickel_evolution(xni_path, args.abn_dir, swd_path, args.n_epochs,
                          args.output, model_name)


if __name__ == '__main__':
    main()
