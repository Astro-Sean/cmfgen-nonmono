#!/usr/bin/env python3
"""Generate the BX batch: diagnose the RSG_M5_V07k_W1m4G_windctl flash.

RSG_M5_V07k_W1m4G_windctl (compact n_rho=11 shell, mdot=1e-4 wind, no
ejecta) is the only control that produces a bright r-band event:
a ~-12.5 floor for ~20 d, a ~4 mag jump to -17.6 at t~24 d, then a
slow decay.  The flash coincides with the shock reaching ~2.2e15 cm --
exactly the graft seam where the windctl extension was spliced onto the
parent's wind->ejecta transition ramp (a ~40% rho bump plus a velocity
kick 850->960 km/s).

So the flash is either:
  a) physical - shock emergence once the thin-wind column ahead can no
     longer trap the accumulated internal energy (K&K breakout class), or
  b) an artifact - the shock crossing the graft seam.

These no-ejecta models decide between them:

  W1m4L : seamless mdot=1e-4 wind to 4.3e16 (same outer radius the
          grafted windctl had, but a single clean create_wind -- no
          seam).  Flash persists -> physical; gone -> seam artifact.
  W1m4S : same wind ending at 2.4e15 (edge at the seam radius).
          Shock exiting into vacuum there -> clean edge-breakout
          reference for comparison with the other two.
  W2m4L : mdot=2e-4 seamless long wind - mdot sensitivity (does a
          slightly thicker column smother the emergence?).
  W05m4L: mdot=5e-5 seamless long wind - thinner still: earlier or
          brighter flash?

All are diagnostics of the flash mechanism only -- they cannot be
SN_B solutions (no ejecta -> no metal lines).
"""
import os
import sys

import yaml

BASE = "/viper/ptmp1/sbrennan/SN2024pba"
CONFIGS = os.path.join(BASE, "configs")
PARENT = os.path.join(CONFIGS, "RSG_M5_V07k_W1m4G.yaml")

R_LONG = 4.3e16  # same outer edge the paired windctl grafted to

# (name, mdot, wind_r_outer_cm, rationale)
VARIANTS = [
    ("BX_M5_V07k_W1m4L", 1e-4, R_LONG,
     "seamless long thin wind, no graft seam: flash persists -> "
     "physical emergence, gone -> seam artifact"),
    ("BX_M5_V07k_W1m4S", 1e-4, 2.4e15,
     "thin wind ending at the seam radius: clean edge-breakout "
     "reference"),
    ("BX_M5_V07k_W2m4L", 2e-4, R_LONG,
     "2x denser wind: does the extra column smother the emergence?"),
    ("BX_M5_V07k_W05m4L", 5e-5, R_LONG,
     "2x thinner wind: earlier/brighter emergence flash?"),
]


def main():
    with open(PARENT) as f:
        parent = yaml.safe_load(f)

    for name, mdot, r_out, why in VARIANTS:
        if len(name) > 28:
            sys.exit(f"name too long: {name} ({len(name)} chars)")
        cfg = dict(parent)
        cfg["model_name"] = name
        cfg["description"] = (
            f"BX flash diagnostic (no ejecta): compact n_rho=11 shell "
            f"(M=5 Msun, v=7000 km/s), mdot={mdot:g} wind to "
            f"r_outer={r_out:.2e} cm.  {why}."
        )
        x = dict(parent["x_model_from_sn"])
        x.update(
            {
                "no_ejecta": True,
                "wind_mdot_msun_yr": mdot,
                "wind_r_outer_cm": r_out,
                "paired_wind_control": False,
            }
        )
        cfg["x_model_from_sn"] = x

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_cocoon/{name}"
        cfg["output"] = out

        path = os.path.join(CONFIGS, f"{name}.yaml")
        with open(path, "w") as f:
            f.write("# BX batch: diagnose the windctl flash -- physical "
                    "thin-wind emergence or graft-seam artifact?\n")
            f.write(f"# {why}\n")
            f.write(f"# Parent config: {PARENT}\n")
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        print("wrote", path)


if __name__ == "__main__":
    main()
