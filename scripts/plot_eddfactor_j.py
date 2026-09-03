#!/usr/bin/env python3
"""
Plot J field from an EDDFACTOR file and optionally compare with W*B(T,nu).

Reads an EDDFACTOR binary file (produced by generate_eddfactor.py or CMFGEN)
and plots:
  1. J(nu, r) as a 2D heatmap (frequency vs depth)
  2. J at selected frequencies vs depth
  3. J at selected depths vs frequency
  4. If --grey: comparison of CMFGEN-computed J vs W*B(T,nu) at selected frequencies

Usage:
  python scripts/plot_eddfactor_j.py <EDDFACTOR_path> [--grey]
  python scripts/plot_eddfactor_j.py X_v2/cmfgen/19.59days/EDDFACTOR --grey
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add scripts dir to path for generate_eddfactor imports
sys.path.insert(0, str(Path(__file__).parent))
from generate_eddfactor import (
    read_eddfactor, read_sn_hydro_data,
    compute_dilution_factor, planck_function, compute_grey_j,
)


def plot_j_field(edd_path, output_dir=None, show_grey=False):
    """Generate J field plots from an EDDFACTOR file."""
    edd_path = Path(edd_path)
    info_path = edd_path.parent / (edd_path.name + '_INFO')
    hydro_path = edd_path.parent / 'SN_HYDRO_DATA'

    if output_dir is None:
        output_dir = edd_path.parent / 'plots'
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading EDDFACTOR: {edd_path}")
    edd = read_eddfactor(edd_path, info_path)
    nd = edd['nd']
    ncf = edd['ncf']
    freqs = edd['cont_frequencies']
    j_cont = edd['cont_j']
    r = edd['r']
    v = edd['v']

    print(f"  ND={nd}, NCF={ncf}")
    print(f"  Frequency range: {freqs[0]:.4e} to {freqs[-1]:.4e} (10^15 Hz)")
    print(f"  R range: {r[0]:.4e} to {r[-1]:.4e} cm")

    # Read SN_HYDRO_DATA for grey comparison
    hydro = None
    grey_j = None
    if show_grey and hydro_path.exists():
        print(f"Reading SN_HYDRO_DATA: {hydro_path}")
        hydro = read_sn_hydro_data(hydro_path)
        t_hydro = hydro['t'][:nd]
        r_hydro = hydro['r'][:nd]
        print(f"  T range: {t_hydro[0]:.1f} to {t_hydro[-1]:.1f} K")
        print(f"  Computing grey J = W * B(T, nu)...")
        grey_j = compute_grey_j(freqs, r_hydro, t_hydro)

    # Depth axis: use radius (outer to inner, CMFGEN convention)
    # Convert to 10^15 cm for readable axis
    r_1e15 = r / 1e15

    # ── Plot 1: J(nu, r) 2D heatmap ──────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    # Use log10 of J (J can span many orders of magnitude)
    j_log = np.log10(np.maximum(j_cont, 1e-300))
    im = ax.imshow(j_log, aspect='auto', origin='upper',
                   extent=[freqs[0], freqs[-1], nd-0.5, 0.5],
                   cmap='viridis')
    ax.set_xlabel('Frequency (10^15 Hz)')
    ax.set_ylabel('Depth (1=outer, ND=inner)')
    ax.set_title(f'log10 J(nu, r) from EDDFACTOR\n{edd_path.parent.name}')
    plt.colorbar(im, ax=ax, label='log10 J [erg/cm^2/s/Hz/sr]')
    fig.tight_layout()
    out1 = output_dir / 'j_field_heatmap.png'
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out1}")

    # ── Plot 2: J at selected frequencies vs depth ───────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    # Select frequencies: near-IR, optical, UV, far-UV
    target_freqs = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    colors = plt.cm.tab10(np.linspace(0, 1, len(target_freqs)))
    for tf, color in zip(target_freqs, colors):
        idx = np.argmin(np.abs(freqs - tf))
        if abs(freqs[idx] - tf) / tf > 0.2:
            continue  # skip if no close match
        label = f'nu={freqs[idx]:.3f}e15 Hz'
        ax.semilogy(range(1, nd+1), np.maximum(j_cont[idx, :], 1e-300),
                    label=label, color=color)
        if grey_j is not None:
            ax.semilogy(range(1, nd+1), np.maximum(grey_j[idx, :], 1e-300),
                        '--', color=color, alpha=0.5)
    ax.set_xlabel('Depth (1=outer, ND=inner)')
    ax.set_ylabel('J [erg/cm^2/s/Hz/sr]')
    title = f'J vs depth at selected frequencies\n{edd_path.parent.name}'
    if grey_j is not None:
        title += '\n(solid=EDDFACTOR, dashed=W*B(T,nu))'
    ax.set_title(title)
    ax.legend(fontsize=8, loc='best')
    ax.set_ylim(bottom=1e-100)
    fig.tight_layout()
    out2 = output_dir / 'j_vs_depth.png'
    fig.savefig(out2, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out2}")

    # ── Plot 3: J at selected depths vs frequency ────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    target_depths = [0, nd//4, nd//2, 3*nd//4, nd-1]
    for d_idx in target_depths:
        label = f'depth {d_idx+1} (R={r[d_idx]:.2e} cm)'
        ax.loglog(freqs, np.maximum(j_cont[:, d_idx], 1e-300),
                  label=label)
        if grey_j is not None:
            ax.loglog(freqs, np.maximum(grey_j[:, d_idx], 1e-300),
                      '--', alpha=0.5)
    ax.set_xlabel('Frequency (10^15 Hz)')
    ax.set_ylabel('J [erg/cm^2/s/Hz/sr]')
    title = f'J vs frequency at selected depths\n{edd_path.parent.name}'
    if grey_j is not None:
        title += '\n(solid=EDDFACTOR, dashed=W*B(T,nu))'
    ax.set_title(title)
    ax.legend(fontsize=8, loc='best')
    ax.set_ylim(bottom=1e-100)
    fig.tight_layout()
    out3 = output_dir / 'j_vs_frequency.png'
    fig.savefig(out3, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out3}")

    # ── Plot 4: Ratio EDDFACTOR / grey at selected frequencies ───────
    if grey_j is not None:
        fig, ax = plt.subplots(figsize=(10, 6))
        for tf, color in zip(target_freqs, colors):
            idx = np.argmin(np.abs(freqs - tf))
            if abs(freqs[idx] - tf) / tf > 0.2:
                continue
            ratio = j_cont[idx, :] / np.maximum(grey_j[idx, :], 1e-300)
            label = f'nu={freqs[idx]:.3f}e15 Hz'
            ax.semilogy(range(1, nd+1), np.maximum(ratio, 1e-10),
                        label=label, color=color)
        ax.axhline(1.0, color='black', linestyle=':', alpha=0.5)
        ax.set_xlabel('Depth (1=outer, ND=inner)')
        ax.set_ylabel('J_EDDFACTOR / J_grey')
        ax.set_title(f'J ratio: EDDFACTOR / W*B(T,nu)\n{edd_path.parent.name}')
        ax.legend(fontsize=8, loc='best')
        ax.set_ylim(1e-5, 1e5)
        fig.tight_layout()
        out4 = output_dir / 'j_ratio.png'
        fig.savefig(out4, dpi=150)
        plt.close(fig)
        print(f"  Saved: {out4}")

    print(f"\nAll plots saved to: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot J field from an EDDFACTOR file."
    )
    parser.add_argument('edd_path', help='Path to EDDFACTOR file')
    parser.add_argument('--grey', action='store_true',
                        help='Also compute and plot W*B(T,nu) from SN_HYDRO_DATA for comparison')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory for plots (default: <edd_dir>/plots)')
    args = parser.parse_args()

    if not Path(args.edd_path).exists():
        print(f"ERROR: {args.edd_path} not found", file=sys.stderr)
        sys.exit(1)

    plot_j_field(args.edd_path, output_dir=args.output_dir, show_grey=args.grey)


if __name__ == "__main__":
    main()
