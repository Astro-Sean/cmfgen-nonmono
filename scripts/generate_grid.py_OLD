#!/usr/bin/env python3
"""Generate a grid of SN_12C model configs for SN_A matching.

Grid parameters (3x3x3 = 27 models):
  - H surface mass: [0.03, 0.05, 0.07] Msun
  - Explosion energy (EBurst): [5.0, 6.0, 7.0] (target ~0.4, 0.5, 0.6 foe KE)
  - Ni56 mass: [0.04, 0.06, 0.08] Msun

Based on SN_12C_CSM5 config. Changes:
  - STELLA target_days: 50 (from 100)
  - CMFGEN snapshot at 35 days post-SBO (from 30)
  - mdot fixed at 1e-10 (very low, stability only)
"""
import yaml
import copy
import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path("/viper/ptmp1/sbrennan/SN2024pba")
TEMPLATE_CONFIG = PROJECT_ROOT / "configs" / "SN_12C_CSM5.yaml"
OUTPUT_DIR = PROJECT_ROOT / "configs"

# Grid definition
H_MASSES = [0.03, 0.05, 0.07]
EBURSTS = [5.0, 6.0, 7.0]
NI56_MASSES = [0.04, 0.06, 0.08]

# Target KE estimates (EBurst * ~0.083 foe per unit, approximate)
# EBurst=5.0 -> ~0.4 foe, EBurst=6.0 -> ~0.5 foe, EBurst=7.0 -> ~0.6 foe


def fmt_h(h):
    return f"H{int(h*100):02d}"


def fmt_e(e):
    return f"E{int(e*10):02d}"


def fmt_n(n):
    return f"N{int(n*100):02d}"


def generate_grid():
    with open(TEMPLATE_CONFIG) as f:
        base = yaml.safe_load(f)

    configs = []
    for h_mass in H_MASSES:
        for eburst in EBURSTS:
            for ni56 in NI56_MASSES:
                cfg = copy.deepcopy(base)
                name = f"SN_12C_GRID_{fmt_h(h_mass)}_{fmt_e(eburst)}_{fmt_n(ni56)}"
                cfg["model_name"] = name
                cfg["description"] = (
                    f"Grid model: H={h_mass} Msun, EBurst={eburst} "
                    f"(~{eburst*0.083:.2f} foe KE), Ni56={ni56} Msun. "
                    f"STELLA 50d, CMFGEN at 35d post-SBO."
                )

                # H surface mass
                cfg["kepler_presn"]["h_surface_mass_msun"] = h_mass

                # Explosion energy
                cfg["stella"]["explosion_energy"] = eburst

                # Ni56 mass
                cfg["stella"]["nickel56"]["total_mass_msun"] = ni56

                # STELLA: 50 days
                cfg["stella"]["target_days"] = 50.0

                # CMFGEN: 35 days post-SBO
                cfg["cmfgen_conversion"]["snapshot_days"] = [35]

                # mdot stays very low (1e-10)
                cfg["kepler_presn"]["csm_wind"]["mdot_msun_yr"] = 1.0e-10

                # Output directory
                cfg["output"]["model_dir"] = f"./{name}"

                # Only run STELLA + CMFGEN
                cfg["run_cmfgen"] = True
                cfg["run_sedona"] = False
                cfg["run_sirocco"] = False
                cfg["run_superlite"] = False
                cfg["reuse_stella"] = False
                cfg["restart_stella"] = False

                configs.append((name, cfg))

    return configs


def write_configs(configs):
    written = []
    for name, cfg in configs:
        path = OUTPUT_DIR / f"{name}.yaml"
        with open(path, "w") as f:
            f.write(f"# {name} - grid model for SN_A matching\n")
            f.write(f"# Based on SN_12C_CSM5 with varied H mass, explosion energy, Ni56 mass\n")
            f.write(f"# STELLA 50d, CMFGEN at 35d post-SBO, mdot=1e-10\n\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False, width=120)
        written.append((name, path))
        print(f"  Written: {path}")
    return written


def submit_jobs(written, skip_stella=False, run_cmfgen=True):
    """Submit all grid models through run_model_grid.py."""
    for name, path in written:
        cmd = [
            "python", "scripts/run_model_grid.py",
            str(path),
        ]
        if skip_stella:
            cmd.append("--skip-stella")
        if run_cmfgen:
            cmd.append("--run-cmfgen")
        cmd.extend(["--no-sedona", "--no-sirocco", "--no-superlite"])

        print(f"\nSubmitting {name}...")
        print(f"  Command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd, cwd=str(PROJECT_ROOT),
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            # Extract job IDs from output
            for line in result.stdout.splitlines():
                if "Job ID" in line or "submitted" in line.lower():
                    print(f"  {line.strip()}")
        else:
            print(f"  ERROR: {result.stderr[:500]}")
            print(f"  stdout: {result.stdout[:500]}")


def main():
    print("=== SN_12C Grid Generation for SN_A Matching ===")
    print(f"Grid: {len(H_MASSES)} H masses x {len(EBURSTS)} energies x {len(NI56_MASSES)} Ni56 = "
          f"{len(H_MASSES)*len(EBURSTS)*len(NI56_MASSES)} models")
    print(f"H masses: {H_MASSES} Msun")
    print(f"EBurst: {EBURSTS} (target KE ~{[e*0.083 for e in EBURSTS]} foe)")
    print(f"Ni56: {NI56_MASSES} Msun")
    print()

    configs = generate_grid()
    print(f"Generated {len(configs)} configs")

    written = write_configs(configs)
    print(f"\nWrote {len(written)} config files to {OUTPUT_DIR}")

    if "--submit" in sys.argv:
        skip_stella = "--skip-stella" in sys.argv
        print(f"\nSubmitting jobs (skip_stella={skip_stella})...")
        submit_jobs(written, skip_stella=skip_stella)
    else:
        print("\nDry run — use --submit to submit all jobs")
        print("Example: python scripts/generate_grid.py --submit")


if __name__ == "__main__":
    main()
