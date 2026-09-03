#!/usr/bin/env python3
"""Generate SN_B grid configs varying shell mass and velocity.

Grid: 4x4x2
  Shell mass:    3.0, 5.0, 7.0, 10.0 Msun
  Shell velocity: 3000, 5000, 7000, 10000 km/s
  Wind:          with wind and no-wind (_NW) variants

Source: SN_A_g_T20k_Ni06 (best photometric match to SN 2011dh)
Shell composition: solar (Asplund+2009, X=0.738 Y=0.249 Z=0.013)
Codes: STELLA + SuperLite (no CMFGEN/SEDONA/SIROCCO)
STELLA evolution: 60 days
"""
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_B_grid"
SOURCE_DIR = str(BASE_DIR / "SN_A" / "SN_A_g_T20k_Ni06")

MASSES = [3.0, 5.0, 7.0, 10.0]
VELOCITIES = [3000, 5000, 7000, 10000]
TARGET_DAYS = 60.0


def make_config(model_name, shell_mass, shell_velocity, no_wind=False):
    wind_desc = "no wind (shell in immediate contact with SN_A ejecta)" if no_wind else "wind"
    return {
        # Header comments preserved via description
        "model_name": model_name,
        "description": (
            f"SN_B grid: H/He shell (solar comp) + {wind_desc} + 400d-old SN_A_g_T20k_Ni06 ejecta. "
            f"M_shell={shell_mass} Msun, v_shell={shell_velocity} km/s. STELLA evolved to {TARGET_DAYS:.0f} days."
        ),

        # Code switches: STELLA + SuperLite only
        "run_cmfgen": False,
        "run_sirocco": False,
        "run_sedona": False,
        "run_superlite": True,

        # Workflow control: clean rerun of STELLA
        "skip_stella": False,
        "reuse_stella": False,
        "cancel_jobs": True,
        "force_cancel": True,
        "plot_lightcurve": True,
        "debug": False,
        "debug_days": 19.59,

        "restart_stella": False,

        # X model: shell + wind + old SN_A ejecta
        "x_model_from_sn": {
            "source_dir": SOURCE_DIR,
            "source_age_days": 400.0,
            "shell_mass_msun": shell_mass,
            "shell_velocity_kms": shell_velocity,
            "shell_composition": "solar",
            "shell_xh": 0.7381,        # overridden by composition=solar
            "shell_yhe": 0.2485,        # overridden by composition=solar
            # Dessart+2015 (SN 2010jl) power-law shell structure
            "n_rho": 8.0,
            "v0_kms": 3000.0,
            "n_T": 0.4,
            "T_interface_K": 2000.0,
            "T_max_K": 65000.0,
            "t_start_days": 5.0,
            "ejecta_T_floor_K": 2000.0,
            "ejecta_rho_min": 1.0e-20,
            "no_wind": no_wind,
            "se_transition_dex": 0.5,
            "wind_mdot_msun_yr": 1.0e-2,
            "wind_velocity_kms": 100,
            "wind_T_K": 2000,
            "explosion_energy": 0.0,
            "aq": 10.0,
            "fit_tau": 5.0,
            "meth": 3,
            "temp_floor": 0.03,
            "floor_rad": 0.01,
            "hmin": 1.0e-14,
            "hmax": 3600.0,
        },

        # Ejecta parameters (unused for X model mode, required by schema)
        "ejecta": {
            "M_ejecta_A": 0.0,
            "R_star_cm": 0.0,
            "n_shell": 0.0,
            "V_ejecta_A_initial_kms": 0.0,
        },

        "hydro_backend": "stella",
        "stella": {
            "convection": False,
            "target_days": TARGET_DAYS,
            "nzon": 500,
            "kappa": 0.34,
            "snapshot_interval": 1.0,
            "iout": 10,
            "r_min": 0.0,
            "rho_min": 0.0,
            "stella_root": "/ptmp/sbrennan/stella/OpenStella",
            "partition": "general",
            "memory": "16G",
            "time": "23:59:00",
            "outer_fraction": 0.3,
            "transition_radius": 1e13,
            "explosion_energy": 0.0,
            "explosion_mass": 0.0,
            "explosion_time": 1.0,
            "explosion_delay": 0.0,
            "initial_kinetic_energy": 0.0,  # EKO=0: preserve .hyd velocities
            "nstmax": 2000000,
            "eps": 0.003,
            "hmax": 3600.0,
            "hmin": 1.0e-14,
            "meth": 3,
            "aq": 10.0,
            "bq": 1.0,
            "fit_tau": 5.0,
            "tau_tol": 1.3,
            "disable_scattering": False,
            "stella_tau_min": 1e-4,
            "wacc_r": 1.0,
            "wacc_v": 1.0,
            "wacc_t": 0.5,
            "wacc_rad": 0.5,
            "temp_floor": 0.03,
            "floor_rad": 0.01,
            "use_ronfict": True,
            "nickel56": {
                "enabled": False,
            },
        },

        "temperature": {
            "temp_profile_type": "adiabatic",
            "N_T_ejecta": 0.6,
            "T_inner_max": 100000,
            "T_csm": 5000,
        },

        # CMFGEN conversion (not used since run_cmfgen=False, but kept for consistency)
        "cmfgen_conversion": {
            "heracles_dir": "heracles",
            "snapshot_days": [19.59],
            "epoch_tolerance": 2.0,
            "stella_extra_days": 0,
            "output_file": "cmfgen/SN_HYDRO_DATA",
            "cmfgen_tau_max": 30,
            "cmfgen_tau_min": 0.01,
            "n_grid": 100,
            "adaptive_grid": True,
            "gradient_threshold": 0.5,
            "velocity_max_kms": 10000,
            "v_min_km_s": 2500.0,
            "opacity_model": "stella_nlte",
            "kappa_constant": 0.34,
            "enable_grid_thinning": False,
            "grid_thinning_threshold": 0.005,
            "enable_adaptive_grid": True,
            "use_tau_amr": True,
            "base_grid_method": "log_tau",
            "tau_spacing": "auto",
            "composition_aware": True,
            "composition_grid_insertion": False,
            "composition_gradient_threshold": 0.5,
            "base_grid_target": 100,
            "max_points": 150,
            "max_thin_zones": 0,
            "max_thick_tau": 1e10,
            "min_points": 80,
            "min_nd": 80,
            "sigma_conformance_headroom": 30,
            "max_radius_jump": 0.05,
            "max_temp_jump": 0.10,
            "max_density_jump": 0.15,
            "smooth_velocity": False,
            "no_velocity_smoothing": True,
            "smooth_wiggles": False,
            "wiggle_window": 5,
            "velocity_smooth_window_kms": 50.0,
            "enable_velocity_spacing_optimization": False,
            "include_ne_amr": False,
            "max_ne_jump": 0.10,
            "ionization_front_aware": True,
            "ionization_front_threshold": 0.5,
            "ionization_front_min_points": 5,
            "ionization_front_width": 3,
            "ionization_front_cap": 0.15,
            "resolve_recombination_fronts": True,
            "recombination_t_jump": 0.3,
            "recombination_front_width": 5,
            "max_logtau_jump": 0.25,
            "max_velocity_jump": 0.30,
            "max_sigma_jump": 2.0,
            "max_cell_ratio": 2.0,
            "thin_grid_tau": 1.0,
            "thin_grid_priority": 1.0,
            "outer_boundary_extra_points": 5,
            "outer_boundary_ratio": 1.5,
            "min_dr_r": 1.0e-4,
            "max_dr_r": 0.15,
            "sigma_priority_threshold": 1.0,
            "sigma_priority": 2.0,
            "dsigma_priority_threshold": 0.3,
            "dsigma_priority": 2.0,
            "smooth_sigma": False,
            "smooth_sigma_scatter": False,
            "max_dsigma": 1.0,
            "sigma_smooth_passes": 10,
            "sigma_smooth_window": 0.3,
            "smooth_sigma_global": False,
            "sigma_global_window": 7,
            "sigma_global_max_rel_change": 0.05,
            "sigma_global_preserve_threshold": 2.0,
            "trace_threshold": 1.0e-6,
        },

        # CMFGEN (disabled, kept for future use)
        "cmfgen": {
            "base_dir": "./cmfgen_base",
            "model_type": "static_snapshot",
            "init_num_its": 1000,
            "single_final": True,
            "proxy_populations": True,
            "proxy_populations_v2": False,
            "proxy_populations_min_ion_frac": 1.0e-15,
            "proxy_populations_nlte_solver": True,
            "inc_rad_decays": True,
            "vturb": 20.0,
            "incl_rel": False,
            "vel_law": 12,
            "depth_independent_solver": True,
            "vadat_overrides": {
                "OB_RAT": "1.0D0",
                "STOP_IF_BP": "F",
                "THK_CONT": "F",
                "IB_STAB": "0.1D0",
                "DIF": "T",
                "IB_METH": "DIFFUSION",
                "FIX_T": "T",
                "FIX_T_AUTO": "F",
                "FIX_IMP": "T",
                "DO_LAM_AUTO": "F",
                "DO_LAM_IT": "T",
                "N_IB_INS": "0",
                "N_OB_INS": "5",
                "BA_CHK_FAC": "5.0D-02",
                "REV_RGRID": "F",
                "MIN_CF": "1.0D-2",
                "T_INIT_TAU": "20.0D0",
                "MAX_LAM": "1.2D0",
                "MAX_LIN": "1.2D0",
                "LAM_VAL": "1.0D2",
                "NUM_LAM": "6",
                "N_FIX_BA": "5",
                "ALL_FREQ": "T",
                "AUTO_SMOOTH": "F",
                "DO_AV": "T",
                "DO_NG": "F",
                "USE_LAM_ES": "T",
            },
            "partition": "general",
            "time": "23:59:00",
            "mem": "16G",
            "init_mem": "16G",
            "init_partition": "general",
            "init_time": "23:59:00",
            "n_ob_ins": 30,
            "final_requeue_when": "finish",
            "final_requeue_max": 500,
            "auto_set_rstar": False,
            "auto_set_lstar": False,
            "phot_tau": 0.6666667,
            "bb_teff_k": 9100,
            "lstar_radius_mode": "RMAX",
        },

        # SuperLite (MC-aware, direct .swd path)
        "superlite": {
            "partition": "general",
            "ntasks": 1,
            "cpus_per_task": 8,
            "memory": "16G",
            "time": "06:00:00",
            "install_dir": "/ptmp/sbrennan/superlite",
            "exe": "/ptmp/sbrennan/superlite/superlite",
            "truncate_tau": 50,
            "input_source": "swd",
            "max_cells": 250,
            "t_floor_K": 2000.0,
        },

        "output": {
            "model_dir": f"./SN_B_grid/{model_name}",
            "stella_output": "stella",
            "cmfgen_output": "cmfgen",
            "sedona_output": "sedona",
            "superlite_output": "superlite",
            "plots": "plots",
            "outputs": "output",
        },
    }


def main():
    GRID_DIR.mkdir(parents=True, exist_ok=True)
    configs = []
    for mass in MASSES:
        for vel in VELOCITIES:
            m_str = f"{mass:.1f}"
            v_str = f"{vel}"
            for no_wind in (False, True):
                suffix = "_NW" if no_wind else ""
                model_name = f"M{m_str}_V{v_str}{suffix}"
                model_dir = GRID_DIR / model_name
                model_dir.mkdir(parents=True, exist_ok=True)

                config = make_config(model_name, mass, vel, no_wind=no_wind)
                config_path = model_dir / "config.yaml"

                wind_label = "NO WIND" if no_wind else "with wind"
                header = (
                    f"# SN_B grid model: {model_name} ({wind_label})\n"
                    f"# Shell: M={mass} Msun, v={vel} km/s, solar composition (Asplund+2009)\n"
                    f"# Source: SN_A_g_T20k_Ni06 (best photometric match to SN 2011dh)\n"
                    f"# Structure: H/He shell + {wind_label} + 400-day-old SN_A ejecta\n"
                    f"# STELLA evolution: {TARGET_DAYS:.0f} days\n"
                    f"# Codes: STELLA + SuperLite\n"
                    f"\n"
                )
                with open(config_path, 'w') as f:
                    f.write(header)
                    yaml.dump(config, f, default_flow_style=False, sort_keys=False)
                configs.append(str(config_path))
                print(f"  Written: {config_path}")

    print(f"\nGenerated {len(configs)} configs in {GRID_DIR}/")
    print("\nTo run the full grid (2 at a time via dependency chaining):")
    print("  bash scripts/submit_sn_b_grid.sh")


if __name__ == '__main__':
    main()
