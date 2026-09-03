#!/usr/bin/env python3
"""Plot level populations under different proxy_populations assumptions.

Computes the ground-state level population for representative species under
each of the four proxy_populations approaches that were tried, using the SAME SN_HYDRO_DATA
so the only difference is the assumption (DI source, DC formula, temperature).

The population formula (CMFGEN's LTEPOP_WLD_V2 with FIX_T=T, so T_local is
always used regardless of what T is written in the _IN file):

    n_gs(j) ∝ DC_gs(j) * DI(j) * ED(j) * T_local(j)^{-1.5}
              * exp(-edge_gs * HDKT / T_local(j))

where DI is the next-higher-stage ion density (cm^-3), DC is the departure
coefficient, ED is electron density, T_local is the gas temperature (10^4 K),
and edge_gs is the ground-state edge frequency (10^15 Hz).

The four assumptions compared:

1. **LTE (CMFGEN LTE_EST)**: DI from full Saha at T_local, DC from
   SET_DC_LTE_V2 (DC=1 for T > TMIN_LOC, depleted below).
2. **Nebular approximation**: DI from dilute Saha at T_ion with W, DC=1.
   (The current proxy_populations.py approach.)
3. **Dilute LTE + analytic DC**: DI from dilute Saha at T_ion with W,
   DC from compute_dc_dilute_lte (exp(edge*(1/T_local - 1/T_ion))).
   (The "exact analytic DC" that gave MAXCH=9.53e149.)
4. **Old geometric-W proxy**: DI from full Saha at T_local, DC from
   compute_dc with T_exc solved from ground-state DC=W.
   (The original proxy approach that gave MAXCH~1e23.)
"""

import os
import sys
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Import proxy_populations functions
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proxy_populations as pig


# Representative species to plot (name, element_label, atomic_mass, ion_stage)
# Chosen to span ionization energies and abundance regimes
PLOT_SPECIES = [
    ('HI',    'HYD',  1.0,  1, 'H I',    'H'),
    ('HeI',   'HE',   4.0,  1, 'He I',   'He'),
    ('Ca2',   'CAL', 40.1,  2, 'Ca II',  'Ca'),
    ('Fe2',   'IRON',55.8,  2, 'Fe II',  'Fe'),
    ('FeIII', 'IRON',55.8,  3, 'Fe III', 'Fe'),
]

# Method labels and colors
METHODS = [
    ('LTE (CMFGEN LTE_EST)',       '#1f77b4', '-'),   # blue
    ('Nebular (DC=1, dilute Saha)', '#2ca02c', '-'),   # green
    ('Dilute LTE + analytic DC',   '#d62728', '-'),   # red
    ('Old geometric-W proxy',      '#ff7f0e', '--'),  # orange dashed
]


def compute_populations(hydro, w, t_ion, tau, species_info, model_dir):
    """Compute ground-state populations under all four assumptions.

    Returns dict: method_name -> np.array of ground-state population (cm^-3).
    """
    name, elem_label, atomic_mass, ion_stage = species_info[:4]

    # Parse atomic data
    oscdat_path = os.path.join(model_dir, f'{name}_F_OSCDAT')
    if not os.path.exists(oscdat_path):
        return None
    atom_data = pig.parse_f_oscdat(oscdat_path)

    # Parse MODEL to get N_F
    model_path = os.path.join(model_dir, 'MODEL')
    if os.path.exists(model_path):
        species_list = pig.parse_model_species(model_path)
    else:
        species_list = pig.parse_model_spec_species(os.path.join(model_dir, 'MODEL_SPEC'))

    nlev = None
    for sp in species_list:
        if sp['name'] == name:
            nlev = sp['n_f']
            break
    if nlev is None:
        nlev = atom_data['n_levels']

    fedge = atom_data['fedge'][:nlev]
    edge_gs = fedge[0]  # ground state edge frequency in 10^15 Hz

    nd = hydro['nd']
    t_local = hydro['temperature']  # 10^4 K
    ed = np.array(hydro['ed'])
    t_local_arr = np.array(t_local)

    # --- DI: next higher stage ion density ---
    # 1. Full Saha at T_local (for LTE and old geometric-W)
    di_full_tlocal, _ = pig.compute_saha_ion_frac(
        elem_label, ion_stage + 1, hydro, atomic_mass,
        t_override=None, min_ion_frac=1e-10)
    di_full_tlocal = np.array(di_full_tlocal)

    # 2. Dilute Saha at T_ion with W (for nebular and dilute LTE + analytic DC)
    di_dilute_tion, _ = pig.compute_saha_ion_frac(
        elem_label, ion_stage + 1, hydro, atomic_mass,
        t_override=t_ion, min_ion_frac=1e-10, w_dilute=w)
    di_dilute_tion = np.array(di_dilute_tion)

    # --- DC: ground-state departure coefficient ---
    # 1. LTE DC (SET_DC_LTE_V2)
    dc_lte = pig.compute_dc_lte_cmfgen(nlev, nd, t_local, fedge, t_exc_min=0.3)
    dc_gs_lte = np.array([dc_lte[0][i] for i in range(nd)])

    # 2. Nebular: DC=1
    dc_gs_nebular = np.ones(nd)

    # 3. Dilute LTE + analytic DC (compute_dc_dilute_lte)
    # Need di_ratio = DI_full/DI_dilute per depth
    # di_ratio accounts for W dilution in ionization balance
    dc_dilute = pig.compute_dc_dilute_lte(
        nlev, nd, t_local, t_ion, fedge, w, t_exc_min=0.3, di_ratio=None)
    dc_gs_dilute = np.array([dc_dilute[0][i] for i in range(nd)])

    # 4. Old geometric-W proxy: DC from compute_dc with T_exc from ground-state DC=W
    t_exc = [pig.get_excite_temp(t_local[i], w[i], edge_gs) for i in range(nd)]
    dc_old = pig.compute_dc(nlev, nd, t_local, t_exc, fedge, w, t_exc_min=0.3)
    dc_gs_old = np.array([dc_old[0][i] for i in range(nd)])

    # --- Compute populations ---
    # n_gs ∝ DC_gs * DI * ED * T^{-1.5} * exp(-edge_gs * HDKT / T)
    # (constant prefactor is the same across all methods, so omitted for comparison)
    boltz = np.array([math.exp(-edge_gs * pig.HDKT / max(t, 0.01)) for t in t_local])
    t_factor = t_local_arr ** (-1.5)

    pops = {}
    pops['LTE (CMFGEN LTE_EST)'] = dc_gs_lte * di_full_tlocal * ed * t_factor * boltz
    pops['Nebular (DC=1, dilute Saha)'] = dc_gs_nebular * di_dilute_tion * ed * t_factor * boltz
    pops['Dilute LTE + analytic DC'] = dc_gs_dilute * di_dilute_tion * ed * t_factor * boltz
    pops['Old geometric-W proxy'] = dc_gs_old * di_full_tlocal * ed * t_factor * boltz

    return pops


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Plot proxy_populations population comparison')
    parser.add_argument('--model-dir', default='models/X_1/cmfgen/5.0days',
                        help='CMFGEN model directory with F_OSCDAT and MODEL files')
    parser.add_argument('--hydro-data', default=None,
                        help='SN_HYDRO_DATA path (default: <model-dir>/SN_HYDRO_DATA)')
    parser.add_argument('--output', default='plots/proxy_population_comparison.png',
                        help='Output plot path')
    parser.add_argument('--tau-trans', type=float, default=1.0)
    args = parser.parse_args()

    model_dir = args.model_dir
    hydro_path = args.hydro_data or os.path.join(model_dir, 'SN_HYDRO_DATA')
    output_path = args.output

    print(f"Model dir:  {model_dir}")
    print(f"Hydro data: {hydro_path}")
    print(f"Output:     {output_path}")

    # Parse hydro data
    hydro = pig.parse_sn_hydro_data(hydro_path)
    nd = hydro['nd']
    print(f"  ND={nd}, L={hydro['lum_lsun']:.3e} Lsun, "
          f"t={hydro['time_days']:.2f} days")

    # Compute radiation field quantities (same as generate_proxy_populations)
    tau = pig.compute_tau(hydro)
    r_ph_10, ph_idx = pig.find_photosphere(hydro, tau, tau_target=2.0/3.0)
    w = pig.compute_w_eff(hydro['radius'], r_ph_10, tau, tau_trans=args.tau_trans)
    t_rad = pig.compute_t_rad(hydro, tau=tau, r_ph_10cm=r_ph_10, tau_trans=args.tau_trans)
    t_ion = pig.compute_t_ion(hydro['temperature'], t_rad, tau, tau_trans=args.tau_trans)

    # Convert to inner->outer for plotting (x-axis increases left-to-right)
    r_plt = np.array(hydro['radius'][::-1]) * 1e10  # cm
    t_plt = np.array(hydro['temperature'][::-1]) * 1e4  # K
    tau_plt = np.array(tau[::-1])
    w_plt = np.array(w[::-1])
    t_ion_plt = np.array(t_ion[::-1]) * 1e4
    t_rad_plt = np.array(t_rad[::-1]) * 1e4

    # Compute populations for each species
    all_pops = {}
    for sp_info in PLOT_SPECIES:
        name = sp_info[0]
        print(f"  Computing populations for {name}...")
        pops = compute_populations(hydro, w, t_ion, tau, sp_info, model_dir)
        if pops is not None:
            all_pops[name] = pops

    if not all_pops:
        print("ERROR: No populations computed. Check model_dir and F_OSCDAT files.")
        sys.exit(1)

    # --- Plot: 2-column layout ---
    # Left column: 2 context panels + 3 species
    # Right column: 3 species + 1 legend/info panel
    n_species = len(all_pops)
    n_rows = max(4, (n_species + 2 + 1) // 2)  # enough rows for all panels
    fig, axes = plt.subplots(n_rows, 2, figsize=(20, 4.5 * n_rows),
                             squeeze=False)
    fig.suptitle('Proxy-In Level Populations: Different Assumptions\n'
                 f'(t={hydro["time_days"]:.1f} d, ND={nd}, '
                 f'R_ph={r_ph_10*1e10:.2e} cm)',
                 fontsize=14, y=0.998)

    xlim = (r_plt.min()*0.95, r_plt.max()*1.05)

    # Panel (0,0): Temperature profiles
    ax1 = axes[0][0]
    ax1.semilogy(r_plt, t_plt, 'k-', linewidth=2, label='T_local')
    ax1.semilogy(r_plt, t_rad_plt, 'b--', linewidth=1.5, label='T_rad')
    ax1.semilogy(r_plt, t_ion_plt, 'r-', linewidth=1.5, alpha=0.7, label='T_ion')
    ax1.axvline(r_ph_10 * 1e10, color='gray', linestyle=':', alpha=0.5, label='R_ph (τ=2/3)')
    ax1.set_ylabel('Temperature (K)', fontsize=11)
    ax1.set_xlim(xlim)
    ax1.legend(loc='lower right', fontsize=9)
    ax1.set_title('Temperature Structure', fontsize=12)
    ax1.grid(True, alpha=0.3)

    # Panel (0,1): W and tau profiles
    ax2 = axes[0][1]
    ax2b = ax2.twinx()
    ax2.semilogy(r_plt, w_plt, 'g-', linewidth=2, label='W_eff')
    ax2b.semilogy(r_plt, tau_plt, 'k--', linewidth=1.5, label='τ')
    ax2.set_ylabel('W_eff (dilution)', fontsize=11, color='g')
    ax2b.set_ylabel('Optical depth τ', fontsize=11)
    ax2.set_xlim(xlim)
    ax2.axvline(r_ph_10 * 1e10, color='gray', linestyle=':', alpha=0.5)
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc='center right', fontsize=9)
    ax2.set_title('Radiation Dilution & Optical Depth', fontsize=12)
    ax2.grid(True, alpha=0.3)

    # Species population panels — fill remaining slots
    panel_idx = 1  # next available row in left column
    col = 0
    legend_drawn = False
    for name, pops in all_pops.items():
        if panel_idx >= n_rows:
            panel_idx = 0
            col = 1
        ax = axes[panel_idx][col]
        sp_info = next(s for s in PLOT_SPECIES if s[0] == name)
        display_name = sp_info[4]

        for (method, color, ls) in METHODS:
            if method in pops:
                pop_arr = np.array(pops[method][::-1])  # inner->outer
                pop_arr = np.maximum(pop_arr, 1e-300)
                ax.semilogy(r_plt, pop_arr, color=color, linestyle=ls,
                           linewidth=1.5, label=method)

        ax.axvline(r_ph_10 * 1e10, color='gray', linestyle=':', alpha=0.5)
        ax.set_ylabel('n_gs (cm⁻³)', fontsize=10)
        ax.set_xlim(xlim)
        ion_energy = pig.IONIZATION_ENERGIES_CM.get(
            (sp_info[1], sp_info[3]), 0)
        chi_ev = ion_energy * 1.2398e-4  # cm^-1 to eV
        ax.set_title(f'{display_name} (χ={chi_ev:.1f} eV)', fontsize=12)
        ax.grid(True, alpha=0.3)
        if not legend_drawn:
            ax.legend(loc='best', fontsize=8, ncol=2)
            legend_drawn = True
        panel_idx += 1

    # Label x-axis on bottom panels
    for col in range(2):
        axes[n_rows - 1][col].set_xlabel('Radius (cm)', fontsize=11)

    # Hide unused panels
    for row in range(n_rows):
        for col in range(2):
            ax = axes[row][col]
            if not ax.has_data() and ax.get_title() == '':
                ax.set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.992])
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {output_path}")

    # Also save a per-species comparison table
    print("\n--- Population ratio summary (max/min across methods at each depth) ---")
    for name, pops in all_pops.items():
        sp_info = next(s for s in PLOT_SPECIES if s[0] == name)
        display_name = sp_info[4]
        # Compute max ratio between methods at any depth
        arrays = [np.array(pops[m]) for m, _, _ in METHODS if m in pops]
        if len(arrays) < 2:
            continue
        # Ratio of max to min at each depth
        stacked = np.stack(arrays)
        max_vals = np.max(stacked, axis=0)
        min_vals = np.min(stacked[stacked > 0], axis=0)
        ratios = max_vals / np.maximum(min_vals, 1e-300)
        print(f"  {display_name:8s}: max ratio = {ratios.max():.2e} "
              f"(at depth {np.argmax(ratios)})")


if __name__ == '__main__':
    main()
