# SN_A / SN_B Double-Explosion Simulation Plan

## Overview

This project models a double-explosion scenario involving a stripped-envelope
star (SN_A) followed by a secondary H/He shell ejection (SN_B) into the
SN_A ejecta. The goal is to produce CMFGEN spectra at two epochs:

1. **SN_A spectrum**: At ~30 days post-shock breakout (SBO) of the first explosion
2. **SN_B spectrum**: At ~40 days after the H/He shell is launched into the
   349-day-old SN_A ejecta

---

## SN_A: Type IIb Explosion of a Stripped Star

### Initial Model
- **Pre-SN model**: KEPLER `he8.00-presn` (5.6 Msun He star, located at
  `SN_A/he8.00-presn`)
- **Mass cut**: 1.4 Msun (compact remnant / neutron star)
- **Ejecta mass**: ~4.2 Msun (above the mass cut)
- **H surface layer**: Thin H envelope added to the He star surface
  - H mass fraction: 0.20 (20% H by mass in the surface layer)
  - Layer mass: 0.10 Msun (covers the outermost zones)
  - Other mass fractions rescaled by (1 - 0.20) = 0.80

### Explosion Parameters
- **Explosion energy**: 1.0e51 erg (EBurst=10 in units of 1e50 erg)
- **Explosion mechanism**: Thermal bomb (AMHT=0.5 Msun heated core, tBurst=1 s)
- **EKO**: 0 (preserve .hyd velocities; thermal bomb sets the explosion)
- **Ni56**: 0.06 Msun (SN 2011dh-like), distributed over inner 0.5 Msun
- **Artificial viscosity**: AQ=2.0, BQ=1.0

### STELLA Evolution
- **Duration**: 60 days post-explosion
- **Zones**: 500 (log-spaced in radius)
- **Integrator**: GEAR+BGH (meth=3), EPS=0.003
- **Timestep limits**: HMAX=3600 s (1 hour), HMIN=1e-14 s
- **Max timesteps**: 2,000,000

### CMFGEN Snapshot
- **Snapshot epoch**: 30 days post-SBO (from STELLA .swd)
- **Tolerance**: ±5 days (closest STELLA epoch used)

### Homologous Evolution
- **Target age**: 349 days post-explosion
- **Method**: Scale from the last STELLA epoch using:
  - r(t) = r(t0) * (t / t0)
  - v(t) = v(t0) (constant velocity)
  - rho(t) = rho(t0) * (t0 / t)^3
  - T(t) = T(t0) * (t0 / t) (adiabatic)
- **Purpose**: Provides the CSM structure for SN_B

### Output Location
All SN_A files are in `SN_A/`:
- `SN_A/config.yaml`: Configuration
- `SN_A/he8.00-presn`: KEPLER pre-SN model
- `SN_A/stella/`: STELLA input and output
- `SN_A/cmfgen/`: CMFGEN input and output
- `SN_A/plots/`: Diagnostic plots

---

## SN_B: H/He Shell Explosion Inside SN_A Ejecta

### Concept
At 349 days post-explosion, SN_A is in homologous expansion. A new H/He
shell is placed at the center and launched at 5000 km/s into the SN_A
ejecta (which now acts as CSM).

### Shell Parameters
- **Shell mass**: 1.0 Msun
- **Shell velocity**: 5000 km/s (homologous: v = v0 * r / r_outer)
- **Shell composition**: X_H = 0.7, Y_He = 0.3 (H-dominated)
- **Shell radius**: ~1e9 cm (compact, much smaller than CSM inner radius)
- **No thermal bomb**: Shell velocity is set in the .hyd file (EKO=0)
- **No Ni56**: No radioactive decay in SN_B (zero .xni file)

### CSM (SN_A Ejecta at 349 Days)
- **Source**: SN_A STELLA .swd, homologously evolved to 349 days
- **Structure**: Retains the homologous velocity profile from SN_A
- **Composition**: He-rich ejecta with thin H surface layer

### Initial Model Construction (`prepare_sn_b.py`)
1. Read SN_A STELLA `.swd` output
2. Find epoch closest to 349 days (homologously evolve if needed)
3. Create H/He shell (1.0 Msun, 5000 km/s, X=0.7, Y=0.3)
4. Combine: shell (inner ~150 zones) + CSM (outer ~350 zones) = 500 zones
5. Write `.hyd`, `.abn`, `.xni`, `.dat`, `input_profile`

### STELLA Evolution
- **Duration**: 50 days after shell launch
- **Zones**: 500
- **EKO**: 0 (preserve both shell and CSM velocities)
- **Artificial viscosity**: AQ=4.0 (higher, for shell-CSM shock interaction)
- **Scattering**: Enabled (H-rich shell)

### CMFGEN Snapshot
- **Snapshot epoch**: 40 days after shell launch
- **Tolerance**: ±5 days

### Output Location
All SN_B files are in `SN_B/`:
- `SN_B/config.yaml`: Configuration
- `SN_B/stella/`: STELLA input and output
- `SN_B/cmfgen/`: CMFGEN input and output
- `SN_B/plots/`: Diagnostic plots

---

## Diagnostic Plots

Plots are generated at every major step using `scripts/plot_model.py`:

### Plot Types
1. **Structure plot** (4 panels): Density, velocity, temperature vs radius
2. **Composition plot**: Mass fractions (H, He, C, O, Si, Fe) vs radius
3. **Lightcurve plot**: Luminosity vs time (from STELLA .swd)
4. **Mass fraction distribution**: Mass fractions vs enclosed mass

### SN_A Plot Steps
| Step | Description | When |
|------|-------------|------|
| `stella_input` | KEPLER pre-SN converted to STELLA .hyd | Before STELLA run |
| `stella_output` | STELLA .swd at final epoch | After STELLA completes |
| `homologous` | 349-day homologously evolved snapshot | After homologous evolution |
| `cmfgen_input` | CMFGEN SN_HYDRO_DATA | After .swd -> SN_HYDRO_DATA conversion |

### SN_B Plot Steps
| Step | Description | When |
|------|-------------|------|
| `sn_b_input` | Shell + CSM initial model | Before STELLA run |
| `stella_output` | STELLA .swd at final epoch | After STELLA completes |
| `cmfgen_input` | CMFGEN SN_HYDRO_DATA | After .swd -> SN_HYDRO_DATA conversion |

---

## Execution Sequence

### Step 1: Run SN_A
```bash
cd /viper/ptmp1/sbrennan/SN2024pba
python scripts/run_model_grid.py SN_A/config.yaml
```
This:
1. Converts KEPLER pre-SN to STELLA format (with H surface layer)
2. Submits STELLA job (runs to 60 days)
3. After STELLA: converts 30-day snapshot to CMFGEN SN_HYDRO_DATA
4. After STELLA: homologously evolves to 349 days, converts to SN_HYDRO_DATA
5. Submits CMFGEN jobs for both epochs

### Step 2: Generate SN_A Plots
```bash
python scripts/plot_model.py SN_A --step all
```

### Step 3: Run SN_B (after SN_A STELLA completes)
```bash
python scripts/run_model_grid.py SN_B/config.yaml
```
This:
1. Reads SN_A STELLA .swd output
2. Homologously evolves to 349 days
3. Creates H/He shell (1.0 Msun, 5000 km/s)
4. Combines shell + CSM into STELLA .hyd format
5. Submits STELLA job (runs 50 days)
6. After STELLA: converts 40-day snapshot to CMFGEN SN_HYDRO_DATA
7. Submits CMFGEN job

### Step 4: Generate SN_B Plots
```bash
python scripts/plot_model.py SN_B --step all
```

---

## Key Scripts

| Script | Purpose |
|--------|---------|
| `scripts/run_model_grid.py` | Main workflow orchestrator |
| `scripts/convert_kepler_presn_to_stella.py` | KEPLER pre-SN -> STELLA conversion |
| `scripts/prepare_sn_b.py` | SN_B initial model (shell + CSM) |
| `scripts/homologous_evolution.py` | Late-time homologous evolution |
| `scripts/convert_stella_swd_to_cmfgen.py` | STELLA .swd -> CMFGEN SN_HYDRO_DATA |
| `scripts/plot_model.py` | Diagnostic plotting |
| `scripts/stella_job_template.slurm` | STELLA SLURM template |

## Key Configuration Files

| File | Purpose |
|------|---------|
| `SN_A/config.yaml` | SN_A model configuration |
| `SN_B/config.yaml` | SN_B model configuration |
| `SN_A/he8.00-presn` | KEPLER pre-SN model (input) |

---

## Physics Notes

### Why EKO=0 for Both Models
EKO (initial_kinetic_energy) in STELLA's .dat file controls velocity handling:
- **EKO > 0**: STELLA overwrites all .hyd velocities with a triangular profile
- **EKO = 0**: STELLA preserves the .hyd velocity profile as-is

For SN_A, the thermal bomb (EBurst=10) provides the explosion energy while
EKO=0 preserves the zero initial velocities (STELLA sets velocities during
the bomb). For SN_B, EKO=0 preserves both the shell velocity (5000 km/s)
and the CSM homologous velocity profile.

### Why AQ=4.0 for SN_B
The H/He shell launched at 5000 km/s creates a strong shock when it hits
the slower CSM (~500-6000 km/s). The higher artificial viscosity (AQ=4.0
vs the default 2.0) stabilizes the GEAR solver during this shock interaction.

### Homologous Evolution Validity
At 349 days post-explosion, the ejecta is in free expansion (v = r/t),
making homologous evolution accurate. The key assumptions:
- Each fluid element moves at constant velocity (no forces)
- Density drops as t^-3 (spherical expansion)
- Temperature drops as t^-1 (adiabatic, radiation-dominated)

These are valid when:
- The expansion is supersonic (v >> sound speed)
- No ongoing energy injection (Ni56 decay is negligible at 349 days)
- No significant CSM interaction (the ejecta is freely expanding)
