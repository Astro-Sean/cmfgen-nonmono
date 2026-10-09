#!/usr/bin/env python3
"""Generate SN_A grid configs for matching the SN 2024pba first-peak lightcurve.

Based on SN_A_v4 (best SN 2011dh match: KEPLER s12.0 stripped, 1.6 Msun He
core + 0.10 Msun H envelope at 200 Rsun, T_base=50kK, thin CSM wind).

SN 2024pba SN_A peaks at about the same luminosity as SN_A_v4 (M_r ~ -17.3)
but declines slower (~0.02 vs ~0.04 mag/day) and is redder.  The grid varies
the two knobs that most directly control peak width and tail slope:

  Grid: 3x3
    E_burst (stella.explosion_energy): 5, 9, 13   (~0.4 / 0.7 / 1.0 foe KE)
    Ni56 mass (stella.nickel56.total_mass_msun): 0.04, 0.075, 0.15 Msun

Models are named E<ee>_Ni<nnn>:  E05_Ni04, E09_Ni075, E13_Ni15, ...
E09_Ni075 is the SN_A_v4 baseline (STELLA-only).

All models use the parallel (PARDISO) OpenStella build; run_model_grid.py
also defaults to it.  Codes: STELLA only — the lightcurve comparison uses
the STELLA .tt UBVRI output plus the post-run pysynphot ugriz extraction.

Observational target: observational_data/SN_A_lightcurve.dat (SDSS g/r;
no i-band data exist for this event).
"""
import os
import copy
import yaml
from pathlib import Path

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_A_pba_grid"
TEMPLATE = BASE_DIR / "SN_A_v4" / "config.yaml"
PRESN_SOURCE = BASE_DIR / "SN_A" / "s12.0-presn-stripped"
PARALLEL_STELLA = "/ptmp/sbrennan/stella/OpenStella_parallel"

E_BURSTS = [5.0, 9.0, 13.0]
NI56_MASSES = [0.04, 0.075, 0.15]


def tag_e(e):
    return f"E{int(round(e)):02d}"


def tag_ni(ni):
    return f"Ni{ni:g}".replace(".", "p")


def make_config(base, model_name, e_burst, ni56):
    cfg = copy.deepcopy(base)

    cfg["model_name"] = model_name
    cfg["description"] = (
        f"SN_A_pba grid: SN_A_v4 physics (1.6 Msun He core + 0.10 Msun H env "
        f"at 200 Rsun) with E_burst={e_burst} (~{e_burst/13:.1f} foe) and "
        f"Ni56={ni56} Msun. Target: SN 2024pba first-peak g/r lightcurve. "
        f"Parallel STELLA, 65 d."
    )

    # Lightcurve grid: STELLA only (UBVRI .tt + pysynphot ugriz post-processing)
    cfg["run_cmfgen"] = False
    cfg["run_sirocco"] = False
    cfg["run_sedona"] = False
    cfg["run_superlite"] = False

    cfg["skip_stella"] = False
    cfg["reuse_stella"] = False
    cfg["restart_stella"] = False
    cfg["cancel_jobs"] = True
    cfg["force_cancel"] = True
    cfg["plot_lightcurve"] = True
    cfg["debug"] = False

    # Grid axes
    cfg["stella"]["explosion_energy"] = e_burst
    cfg["stella"]["nickel56"]["total_mass_msun"] = ni56

    # Parallel (PARDISO) OpenStella build
    cfg["stella"]["stella_root"] = PARALLEL_STELLA
    cfg["stella"]["cpus_per_task"] = 3

    # Per-model paths
    cfg["kepler_presn"]["presn_file"] = (
        f"/ptmp/sbrennan/SN2024pba/SN_A_pba_grid/{model_name}/s12.0-presn-stripped")
    cfg["output"]["model_dir"] = f"./SN_A_pba_grid/{model_name}"

    return cfg


def main():
    with open(TEMPLATE) as f:
        base = yaml.safe_load(f)

    GRID_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Generating {len(E_BURSTS) * len(NI56_MASSES)} SN_A_pba grid configs "
          f"in {GRID_DIR}/\n")
    print(f"{'Model':<16} {'E_burst':>8} {'Ni56':>6}")
    print("-" * 34)

    configs = []
    for e in E_BURSTS:
        for ni in NI56_MASSES:
            model_name = f"{tag_e(e)}_{tag_ni(ni)}"
            model_dir = GRID_DIR / model_name
            model_dir.mkdir(parents=True, exist_ok=True)

            # Per-model symlink to the KEPLER pre-SN (read-only input;
            # converter writes via --outdir into the model dir)
            kep_link = model_dir / "s12.0-presn-stripped"
            if not kep_link.exists():
                os.symlink(PRESN_SOURCE, kep_link)

            cfg = make_config(base, model_name, e, ni)
            config_path = model_dir / "config.yaml"

            header = (
                f"# SN_A_pba grid model: {model_name}\n"
                f"# Base: SN_A_v4 (KEPLER s12.0 stripped + 0.10 Msun H env at 200 Rsun)\n"
                f"# Grid axis: E_burst={e}, Ni56={ni} Msun\n"
                f"# Target: SN 2024pba first-peak lightcurve (obs SN_A_lightcurve.dat, g/r)\n"
                f"# Codes: STELLA only (parallel OpenStella build)\n\n"
            )
            with open(config_path, "w") as f:
                f.write(header)
                yaml.dump(cfg, f, default_flow_style=False, sort_keys=False,
                          width=120)

            configs.append(str(config_path))
            print(f"{model_name:<16} {e:8.1f} {ni:6.3f}")

    print(f"\nGenerated {len(configs)} configs.")
    print("\nTo submit the grid (STELLA jobs auto-chained, max 2 concurrent):")
    print("  python scripts/run_model_grid.py SN_A_pba_grid/*/config.yaml")


if __name__ == "__main__":
    main()
