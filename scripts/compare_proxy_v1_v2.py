#!/usr/bin/env python3
"""Compare proxy_populations.py (V1) and proxy_populations_v2.py (V2) output.

Runs both proxies on the same model directory (or uses existing output dirs),
then generates diagnostic plots and a machine-readable comparison table.

Plots produced (in <output_dir>/plots/):
  1. DI comparison per species (V1 vs V2 vs existing converged) vs radius
  2. DC comparison per species (V1 vs V2 vs existing) vs radius
  3. ED comparison (V1 vs V2 vs GAMMAS_IN vs existing) vs radius
  4. T comparison (V1 vs V2 vs existing) vs radius
  5. DI ratio (V2/V1) per species vs radius
  6. Summary table (text + JSON)

Usage:
  python scripts/compare_proxy_v1_v2.py <model_dir> [--v1-dir DIR] [--v2-dir DIR]
      [--existing-as-prev] [--output-dir DIR]

If --v1-dir / --v2-dir are omitted, the script runs V1 and V2 on the model
and writes to <output_dir>/v1_out and <output_dir>/v2_out respectively.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Make proxy_populations_v2 importable for read_in_file
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proxy_populations_v2 as pp  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ELEMENT_NAMES = {
    'HYD': 'H', 'HE': 'He', 'OXY': 'O', 'SIL': 'Si',
    'CAL': 'Ca', 'IRON': 'Fe', 'NICK': 'Ni', 'COB': 'Co',
}
ROMAN = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI', 7: 'VII'}


def species_display_name(name: str) -> str:
    """Convert 'HI' -> 'H I', 'FeIII' -> 'Fe III', etc."""
    elem, stage = pp.parse_species_name(name)
    if elem is None:
        return name
    return f"{elem} {ROMAN.get(stage, str(stage))}"


def load_in_file(path: str) -> Optional[pp.InFileData]:
    """Wrapper that returns None on missing/unparseable files."""
    return pp.read_in_file(path)


def read_gammas_in(path: str, nd: int) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Read ED and R from GAMMAS_IN. Returns (ed, radius) or (None, None)."""
    if not os.path.exists(path):
        return None, None
    with open(path) as f:
        lines = f.readlines()
    ed = []
    radius = []
    section = None
    for line in lines:
        s = line.strip()
        if s.startswith('!'):
            if 'Electron density' in s:
                section = 'ed'
            elif 'Radius' in s:
                section = 'radius'
            else:
                section = None
            continue
        if not s or section is None:
            continue
        for tok in s.split():
            try:
                val = float(tok)
            except ValueError:
                break
            if section == 'ed':
                ed.append(val)
            elif section == 'radius':
                radius.append(val)
    ed = np.array(ed[:nd]) if len(ed) >= nd else None
    radius = np.array(radius[:nd]) if len(radius) >= nd else None
    return ed, radius


def run_proxy(script: str, model_dir: str, output_dir: str,
              extra_args: List[str] = None) -> int:
    """Run a proxy script and return its exit code."""
    cmd = [sys.executable, script, model_dir, "--output-dir", output_dir]
    if extra_args:
        cmd.extend(extra_args)
    hydro_data = os.path.join(model_dir, "SN_HYDRO_DATA")
    if os.path.exists(hydro_data) and "proxy_populations.py" in script:
        # V1 needs --hydro-data explicitly when not two levels up
        cmd.extend(["--hydro-data", hydro_data])
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"  ERROR (exit {result.returncode}): {result.stderr[-500:]}")
    return result.returncode


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _setup_matplotlib():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def plot_di_comparison(species_data: Dict, hydro_radius: np.ndarray,
                       output_dir: str, plt):
    """Plot DI vs radius for each species: V1, V2, existing."""
    n = len(species_data)
    n_cols = 3
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4.5 * n_rows),
                             squeeze=False)
    fig.suptitle('DI (next-ion density) comparison: V1 vs V2 vs existing',
                 fontsize=14, y=1.01)
    r_cm = hydro_radius * 1e10  # 10^10 cm -> cm
    for idx, (name, data) in enumerate(sorted(species_data.items())):
        ax = axes[idx // n_cols][idx % n_cols]
        label = species_display_name(name)
        for source, color, ls in [('existing', 'gray', '--'),
                                   ('v1', 'blue', '-'),
                                   ('v2', 'red', '-')]:
            d = data.get(source)
            if d is not None:
                ax.semilogy(r_cm, np.maximum(d.di, 1e-300), color=color,
                            linestyle=ls, linewidth=1.2, label=source)
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.set_xlabel('Radius (cm)', fontsize=9)
        ax.set_ylabel('DI (cm$^{-3}$)', fontsize=9)
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3, which='both')
        ax.set_xscale('log')
    # hide unused
    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)
    plt.tight_layout()
    path = os.path.join(output_dir, 'compare_DI_vs_radius.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")


def plot_dc_comparison(species_data: Dict, hydro_radius: np.ndarray,
                       output_dir: str, plt):
    """Plot ground-state DC vs radius for each species."""
    n = len(species_data)
    n_cols = 3
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4.5 * n_rows),
                             squeeze=False)
    fig.suptitle('Ground-state DC comparison: V1 vs V2 vs existing',
                 fontsize=14, y=1.01)
    r_cm = hydro_radius * 1e10
    for idx, (name, data) in enumerate(sorted(species_data.items())):
        ax = axes[idx // n_cols][idx % n_cols]
        label = species_display_name(name)
        for source, color, ls in [('existing', 'gray', '--'),
                                   ('v1', 'blue', '-'),
                                   ('v2', 'red', '-')]:
            d = data.get(source)
            if d is not None and d.nold > 0:
                ax.semilogy(r_cm, np.maximum(np.abs(d.dc[0]), 1e-300),
                            color=color, linestyle=ls, linewidth=1.2,
                            label=source)
        ax.axhline(1.0, color='green', linestyle=':', alpha=0.5, label='LTE (DC=1)')
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.set_xlabel('Radius (cm)', fontsize=9)
        ax.set_ylabel('DC (ground state)', fontsize=9)
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3, which='both')
        ax.set_xscale('log')
    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)
    plt.tight_layout()
    path = os.path.join(output_dir, 'compare_DC_ground_vs_radius.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")


def plot_ed_comparison(v1_data: Optional[pp.InFileData],
                       v2_data: Optional[pp.InFileData],
                       existing_data: Optional[pp.InFileData],
                       gammas_ed: Optional[np.ndarray],
                       hydro_ed: Optional[np.ndarray],
                       hydro_radius: np.ndarray,
                       output_dir: str, plt):
    """Plot ED comparison from all sources."""
    r_cm = hydro_radius * 1e10
    fig, ax = plt.subplots(figsize=(10, 6))
    if existing_data is not None:
        ax.semilogy(r_cm, np.maximum(existing_data.ed, 1e-300),
                    'gray', '--', linewidth=1.5, label='Existing _IN (converged)')
    if v1_data is not None:
        ax.semilogy(r_cm, np.maximum(v1_data.ed, 1e-300),
                    'blue', '-', linewidth=1.5, label='V1 _IN')
    if v2_data is not None:
        ax.semilogy(r_cm, np.maximum(v2_data.ed, 1e-300),
                    'red', '-', linewidth=1.5, label='V2 _IN (GAMMAS_IN)')
    if gammas_ed is not None:
        ax.semilogy(r_cm, np.maximum(gammas_ed, 1e-300),
                    'green', '-.', linewidth=1.5, label='GAMMAS_IN (source)')
    if hydro_ed is not None:
        ax.semilogy(r_cm, np.maximum(hydro_ed, 1e-300),
                    'orange', ':', linewidth=1.5, label='SN_HYDRO_DATA')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('Electron density (cm$^{-3}$)', fontsize=11)
    ax.set_title('ED source comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xscale('log')
    plt.tight_layout()
    path = os.path.join(output_dir, 'compare_ED_sources.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")


def plot_temperature_comparison(v1_data: Optional[pp.InFileData],
                                v2_data: Optional[pp.InFileData],
                                existing_data: Optional[pp.InFileData],
                                hydro_t: np.ndarray,
                                hydro_radius: np.ndarray,
                                output_dir: str, plt):
    """Plot T comparison from all sources."""
    r_cm = hydro_radius * 1e10
    fig, ax = plt.subplots(figsize=(10, 6))
    if existing_data is not None:
        ax.plot(r_cm, existing_data.t * 1e4, 'gray', '--', linewidth=1.5,
                label='Existing _IN')
    if v1_data is not None:
        ax.plot(r_cm, v1_data.t * 1e4, 'blue', '-', linewidth=1.5, label='V1 _IN')
    if v2_data is not None:
        ax.plot(r_cm, v2_data.t * 1e4, 'red', '-', linewidth=1.5, label='V2 _IN')
    ax.plot(r_cm, hydro_t * 1e4, 'orange', ':', linewidth=1.5,
            label='SN_HYDRO_DATA (T_local)')
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('Temperature (K)', fontsize=11)
    ax.set_title('Temperature column comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xscale('log')
    plt.tight_layout()
    path = os.path.join(output_dir, 'compare_T_sources.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")


def plot_di_ratio(species_data: Dict, hydro_radius: np.ndarray,
                  output_dir: str, plt):
    """Plot DI ratio V2/V1 per species."""
    n = 0
    ratios = {}
    for name, data in sorted(species_data.items()):
        v1 = data.get('v1')
        v2 = data.get('v2')
        if v1 is not None and v2 is not None:
            r = v2.di / np.maximum(v1.di, 1e-300)
            ratios[name] = r
            n += 1
    if n == 0:
        return
    n_cols = 3
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4.5 * n_rows),
                             squeeze=False)
    fig.suptitle('DI ratio (V2 / V1) vs radius', fontsize=14, y=1.01)
    r_cm = hydro_radius * 1e10
    for idx, (name, r) in enumerate(ratios.items()):
        ax = axes[idx // n_cols][idx % n_cols]
        label = species_display_name(name)
        ax.semilogy(r_cm, np.maximum(np.abs(r), 1e-10), 'k-', linewidth=1.2)
        ax.axhline(1.0, color='green', linestyle=':', alpha=0.5, label='V2=V1')
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.set_xlabel('Radius (cm)', fontsize=9)
        ax.set_ylabel('DI$_{V2}$ / DI$_{V1}$', fontsize=9)
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3, which='both')
        ax.set_xscale('log')
    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)
    plt.tight_layout()
    path = os.path.join(output_dir, 'compare_DI_ratio_v2_v1.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")


def plot_summary_table(species_data: Dict, output_dir: str, plt):
    """Write a text summary table and a JSON comparison file."""
    rows = []
    for name, data in sorted(species_data.items()):
        row = {'species': name}
        for source in ['existing', 'v1', 'v2']:
            d = data.get(source)
            if d is not None:
                row[f'{source}_di_min'] = float(np.min(d.di))
                row[f'{source}_di_max'] = float(np.max(d.di))
                row[f'{source}_dc_gs_min'] = float(np.min(np.abs(d.dc[0]))) if d.nold > 0 else None
                row[f'{source}_dc_gs_max'] = float(np.max(np.abs(d.dc[0]))) if d.nold > 0 else None
                row[f'{source}_nlev'] = d.nold
            else:
                for k in ['di_min', 'di_max', 'dc_gs_min', 'dc_gs_max', 'nlev']:
                    row[f'{source}_{k}'] = None
        # ratios
        v1 = data.get('v1')
        v2 = data.get('v2')
        if v1 is not None and v2 is not None:
            safe = np.maximum(v1.di, 1e-300)
            row['di_ratio_v2_v1_min'] = float(np.min(v2.di / safe))
            row['di_ratio_v2_v1_max'] = float(np.max(v2.di / safe))
            row['di_log10_diff_max'] = float(np.max(np.abs(
                np.log10(np.maximum(v2.di, 1e-300) / safe))))
        rows.append(row)

    # JSON
    json_path = os.path.join(output_dir, 'comparison_table.json')
    with open(json_path, 'w') as f:
        json.dump(rows, f, indent=2)
    print(f"  JSON table: {json_path}")

    # Text table
    txt_path = os.path.join(output_dir, 'comparison_table.txt')
    with open(txt_path, 'w') as f:
        f.write("V1 vs V2 proxy comparison table\n")
        f.write("=" * 120 + "\n\n")
        f.write(f"{'Species':8s} ")
        for src in ['exist', 'v1', 'v2']:
            f.write(f"| {src} DI min      {src} DI max      {src} DCgs min     {src} DCgs max     {src} NLEV  ")
        f.write("| V2/V1 DI min   V2/V1 DI max   max log10 diff\n")
        f.write("-" * 200 + "\n")
        for row in rows:
            f.write(f"{row['species']:8s} ")
            for src in ['existing', 'v1', 'v2']:
                for key in ['di_min', 'di_max', 'dc_gs_min', 'dc_gs_max', 'nlev']:
                    val = row.get(f'{src}_{key}')
                    if val is None:
                        f.write(f" {'n/a':>13s}")
                    elif key == 'nlev':
                        f.write(f" {int(val):>13d}")
                    else:
                        f.write(f" {val:>13.3e}")
            if 'di_ratio_v2_v1_min' in row:
                f.write(f" {row['di_ratio_v2_v1_min']:>13.3e}"
                        f" {row['di_ratio_v2_v1_max']:>13.3e}"
                        f" {row['di_log10_diff_max']:>13.3f}")
            f.write("\n")
    print(f"  Text table: {txt_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compare V1 and V2 proxy population outputs')
    parser.add_argument('model_dir', help='CMFGEN model directory')
    parser.add_argument('--v1-dir', help='Existing V1 output dir (skip running V1)')
    parser.add_argument('--v2-dir', help='Existing V2 output dir (skip running V2)')
    parser.add_argument('--output-dir', default=None,
                        help='Output dir for comparison plots (default: model_dir/proxy_comparison)')
    parser.add_argument('--v1-extra', nargs='*', default=None,
                        help='Extra args for V1 (e.g. --lte-floor --nlte-solver)')
    parser.add_argument('--v2-extra', nargs='*', default=None,
                        help='Extra args for V2 (e.g. --strategy GREY)')
    args = parser.parse_args()

    model_dir = os.path.abspath(args.model_dir)
    output_dir = args.output_dir or os.path.join(model_dir, 'proxy_comparison')
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    # --- Run proxies if needed ---
    v1_dir = args.v1_dir
    v2_dir = args.v2_dir
    if v1_dir is None:
        v1_dir = os.path.join(output_dir, 'v1_out')
        print("Running V1 (proxy_populations.py)...")
        rc = run_proxy(os.path.join(os.path.dirname(__file__), 'proxy_populations.py'),
                       model_dir, v1_dir, extra_args=args.v1_extra or [])
        if rc != 0:
            print(f"V1 failed (exit {rc}), using whatever output was produced")
    if v2_dir is None:
        v2_dir = os.path.join(output_dir, 'v2_out')
        print("Running V2 (proxy_populations_v2.py)...")
        rc = run_proxy(os.path.join(os.path.dirname(__file__), 'proxy_populations_v2.py'),
                       model_dir, v2_dir, extra_args=args.v2_extra or [])
        if rc != 0:
            print(f"V2 failed (exit {rc}), using whatever output was produced")

    # --- Load hydro data for radius axis ---
    hydro_path = os.path.join(model_dir, 'SN_HYDRO_DATA')
    if not os.path.exists(hydro_path):
        # V1 default: two levels up
        hydro_path = os.path.join(os.path.dirname(os.path.dirname(model_dir)),
                                  'SN_HYDRO_DATA')
    if not os.path.exists(hydro_path):
        print(f"ERROR: SN_HYDRO_DATA not found for {model_dir}")
        sys.exit(1)
    hydro = pp.parse_sn_hydro_data(hydro_path)
    hydro_radius = hydro.radius  # 10^10 cm, outer->inner
    hydro_ed = hydro.ed
    hydro_t = hydro.temperature

    # --- Load GAMMAS_IN ED ---
    gammas_ed, gammas_r = read_gammas_in(os.path.join(model_dir, 'GAMMAS_IN'),
                                         hydro.nd)

    # --- Determine species list from V2 output (or MODEL) ---
    def _in_stem(path: Path) -> str:
        """Extract species name from e.g. 'HI_IN' -> 'HI'."""
        name = path.name
        if name.endswith('_IN'):
            return name[:-3]
        return path.stem

    v2_ins = sorted(Path(v2_dir).glob("*_IN"))
    v2_ins = [f for f in v2_ins if f.name != 'GAMMAS_IN']
    species_names = [_in_stem(f) for f in v2_ins]
    if not species_names:
        # fall back to V1
        v1_ins = sorted(Path(v1_dir).glob("*_IN"))
        v1_ins = [f for f in v1_ins if f.name != 'GAMMAS_IN']
        species_names = [_in_stem(f) for f in v1_ins]
    if not species_names:
        print("ERROR: no _IN files found in V1 or V2 output")
        sys.exit(1)

    print(f"\nComparing {len(species_names)} species: {species_names}")

    # --- Load all _IN files ---
    species_data: Dict[str, Dict] = {}
    for name in species_names:
        data = {}
        # V2
        v2_path = os.path.join(v2_dir, f"{name}_IN")
        d = load_in_file(v2_path)
        if d is not None:
            data['v2'] = d
        # V1
        v1_path = os.path.join(v1_dir, f"{name}_IN")
        d = load_in_file(v1_path)
        if d is not None:
            data['v1'] = d
        # Existing (in model_dir, the converged or previous-epoch file)
        exist_path = os.path.join(model_dir, f"{name}_IN")
        if os.path.exists(exist_path) and os.path.abspath(exist_path) != os.path.abspath(v2_path):
            d = load_in_file(exist_path)
            if d is not None:
                data['existing'] = d
        if data:
            species_data[name] = data

    # --- Generate plots ---
    print(f"\nGenerating plots in {plots_dir}...")
    plt = _setup_matplotlib()

    plot_di_comparison(species_data, hydro_radius, plots_dir, plt)
    plot_dc_comparison(species_data, hydro_radius, plots_dir, plt)
    plot_di_ratio(species_data, hydro_radius, plots_dir, plt)

    # ED and T comparison (use HI as representative)
    hi_data = species_data.get('HI', {})
    v1_hi = hi_data.get('v1')
    v2_hi = hi_data.get('v2')
    exist_hi = hi_data.get('existing')
    plot_ed_comparison(v1_hi, v2_hi, exist_hi, gammas_ed, hydro_ed,
                       hydro_radius, plots_dir, plt)
    plot_temperature_comparison(v1_hi, v2_hi, exist_hi, hydro_t,
                                hydro_radius, plots_dir, plt)

    # Summary table
    plot_summary_table(species_data, output_dir, plt)

    print(f"\nDone. Comparison output in {output_dir}")


if __name__ == '__main__':
    main()
