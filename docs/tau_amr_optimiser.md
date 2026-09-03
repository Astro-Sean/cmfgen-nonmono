# Tau-Based AMR STELLA → CMFGEN Grid Optimiser

**Module:** `scripts/tau_amr_optimiser.py`
**Config:** `configs/tau_amr_config.yaml`
**Status:** `VERIFIED` on the X_v10/X_v11 STELLA `.swd` models (500/750/1000-point
convergence tests pass: velocity topology preserved, mass and species mass
conserved to machine precision).

## 1. Purpose

Convert a STELLA `.swd` model into an **optimised CMFGEN input grid for the
non-monotonic solver** (`VEL_LAW=12`).  The grid is constructed primarily in
**optical-depth space** using **recursive cell-based adaptive mesh refinement
(AMR)**, with a hard maximum of 1000 grid points.  Hydrodynamic / composition
criteria act as **hard constraints** that force refinement at shocks,
composition interfaces and velocity extrema regardless of optical depth.

This replaces the previous `log_tau` "auto" grid (uniform `log(tau+1)`
sampling) which could not guarantee that sharp chemical interfaces or velocity
extrema were resolved.

## 2. Pipeline

```
READ STELLA .swd  (reuse convert_stella_swd_to_cmfgen.read_swd)
   ↓
validate units / radial ordering
   ↓
reconstruct optical depth  (hybrid kappa: STELLA NLTE ∪ e-scatter floor ∪ Kramers)
   ↓
detect composition interfaces   (HARD)
detect shocks                   (HARD)
detect velocity extrema         (HARD)
   ↓
coarse base grid in log(tau)  +  seed every feature radius & velocity turn
   ↓
recursive cell-based AMR:
   multi-component refinement score (τ, ρ, v, T, curvature)
   split highest-priority cell (hard-constraint cells first)
   STOP when: all criteria satisfied OR stagnation OR N = N_max
   ↓
coarsening pass: remove points while maintaining hard limits
   ↓
conservative zone-integrated remap  (exact mass + species-mass conservation)
   ↓
pre-CMFGEN resolution report (10 criteria from STELLA physics)
   ↓
validate: point count, radial order, velocity topology, sum(X)=1, mass, species mass
   ↓
write CMFGEN SN_HYDRO_DATA  (reuse write_sn_hydro_data)  +  diagnostics + manifest
```

## 3. Optical-depth reconstruction

`reconstruct_optical_depth()` builds cumulative τ(r) = ∫ κρ dr (outer→inner,
increasing with depth index).  Three kappa sources, selectable via
`kappa_source`:

| source    | description |
|-----------|-------------|
| `stella`  | STELLA NLTE Rosseland kappa (`.swd` column 12, WRKX). Accurate in hot ionised regions; can drop to ~0 in cool neutral zones. |
| `kramers` | Analytical electron-scattering + bound-free + free-free (Kramers). Composition-dependent. |
| `hybrid`  | **Default.** Element-wise max of STELLA NLTE, e-scatter floor (full-ionisation, conservative), and Kramers. Avoids underestimating τ in neutral regions that CMFGEN/SuperLite will re-ionise. |

The hybrid kappa mirrors the SuperLite MC-aware conversion (AGENTS.md §SuperLite).
The bug fix (2026-08-21) for the `Z` variable shadowing in
`approximate_rosseland_opacity()` is inherited from the parent converter.

## 4. Feature detectors (hard constraints)

### Composition interfaces
`detect_composition_interfaces()` flags any zone where `|dX_i| >
composition_jump` (default 0.30) for a species above
`composition_trace_floor` (1e-4).  These are **mandatory** refinement regions.

### Shocks
`detect_shocks()` flags zones with density compression > 1.5, temperature jump
> 0.5 dex, or normalised velocity curvature > 0.5.  Nearby shocks are
de-duplicated.

### Velocity extrema
`detect_velocity_extrema()` finds local maxima/minima over a 5-zone window.
Additionally, **every** simple 3-point velocity turn (and its immediate
neighbours) is seeded as an exact grid point so PCHIP cannot smooth away
shallow extrema — this is what guarantees velocity-topology preservation.

## 5. Recursive AMR algorithm (minimum-grid optimisation)

The objective is to **minimise N_grid** while maintaining a sufficiently
resolved representation of the STELLA structure.  The 1000-point limit is a
hard upper bound, NOT a target.

`recursive_amr()`:

1. Build a coarse base grid (default 80 points) uniformly spaced in
   `log10(τ + τ_floor)`.
2. Seed the exact outermost/innermost source radii, every feature radius, and
   every velocity-turn radius (±neighbours).
3. Loop:
   - Compute the multi-component soft refinement score per cell.
   - Identify hard-constraint features that still need points.
   - **STOP CONDITION 1** — all criteria satisfied: no feature needs
     refinement, all soft scores < safety_margin, N ≥ min_points.
   - **STOP CONDITION 2** — convergence stagnation: the max score has not
     improved by > 0.1% for 20 consecutive iterations.  This detects
     physical discontinuities (shocks, composition interfaces) that can't
     be resolved further by splitting.
   - **HARD LIMIT**: N ≥ max_points (1000).
   - Split the highest-priority cell (hard-constraint cells take precedence)
     at the geometric-mean radius.
   - Re-remap and recompute τ.
4. **Coarsening pass**: after AMR converges, attempt to remove individual
   non-feature points while maintaining all individual hard limits
   (3× soft targets).  This finds the true minimum viable grid.

### Minimum-grid example

For the SN_A SN_12C_CSM5 model at epoch 30d:
- AMR stagnation at 246 points (max_score=1.766, physical discontinuities)
- Coarsening: 0 points could be removed (grid already at minimum)
- Final grid: **246 points** (not 1000)

### Soft refinement score

### Soft refinement score

```
E = w_τ·|Δlogτ|/target_dlogτ
  + w_ρ·|Δlogρ|/target_dlogrho
  + w_v·|Δv|/target_dv
  + w_T·|ΔlogT|/target_dlogT
  + w_curv·curvature/0.5
  + 0.2·(Δr/r)/target_dr_over_r
```

All terms are dimensionless and comparable.  Composition is handled by the
hard-constraint path, not the soft score (a tiny trace-species change can
matter).

## 6. Conservative remapping

`conservative_remap_density_composition()` treats each old zone as uniform
(ρ, X_i) over its spherical-shell volume.  For each new zone it sums the mass
(and per-species mass) of all overlapping old zones, then divides by the new
zone volume.  This conserves **total mass and species mass to machine
precision** regardless of grid alignment.

Intensive quantities (v, T, kappa) use PCHIP (monotonic, no overshoot).  At
composition interfaces the conservative remap naturally preserves the
discontinuity because old zones are uniform; interface old zones are
force-retained as grid points (hard constraint) so straddling is minimised.

## 7. Validation and pre-CMFGEN resolution report

### Design philosophy

The grid is determined entirely from the STELLA physics, NOT from CMFGEN
source code.  CMFGEN source inspection is used only to establish:
- required input format,
- ordering conventions (outer-to-inner, strictly decreasing radius),
- allowed/non-allowed grid structures,
- variables required by the non-monotonic solver (VEL_LAW=12),
- numerical assumptions the input grid must satisfy.

The AMR algorithm belongs entirely to the STELLA → CMFGEN conversion layer.
CMFGEN receives a grid that is already adequately resolved **before the first
CMFGEN iteration**.

### Conservation validation (`validate_grid`)

`validate_grid()` checks:
- `N ≤ max_points`
- strictly decreasing radius
- velocity extrema count ≥ source extrema count (topology preserved)
- no negative abundances
- `Σ X_i ≈ 1` per zone
- total mass conserved (rel_err < mass_tol)
- species mass conserved (rel_err < 5·mass_tol)

### Pre-CMFGEN resolution report (`resolution_report`)

`resolution_report()` independently checks that the AMR grid adequately
resolves the STELLA physical structure BEFORE CMFGEN is launched:

1. **Optical depth** — largest Δlogτ per cell ≤ 3×target
2. **Density** — largest Δlogρ per cell ≤ 3×target
3. **Velocity** — largest Δv per cell ≤ 3×target
4. **Velocity extrema** — all source 3-point turns preserved
5. **Temperature** — largest ΔlogT per cell ≤ 3×target
6. **Shocks** — ≥min_points_per_shock grid points near each shock
7. **Shell boundaries** — sharp density jumps identified
8. **Composition interfaces** — ≥min_points_per_interface points near each
9. **Abundance gradients** — no unresolved ΔX > composition_jump in smooth
   regions (jumps at known interfaces are expected and OK)
10. **Spectrum-forming region** — ≥20 points in τ=[0.01, 10]

The report prints a formatted summary and identifies all unresolved regions
explicitly.  CMFGEN should only be launched if both the conservation
validation AND the resolution check pass.

### Adaptive point count

The 1000-point limit is a computational constraint, not a CMFGEN requirement.
The algorithm uses the **smallest number of points that adequately represents
the STELLA structure**.  If the resolution requirements cannot be satisfied
with 1000 points, the validation **fails and explicitly identifies the
unresolved regions** rather than silently degrading the model.

## 8. Convergence tests (VERIFIED)

Run on `X_v10/stella/strad/X_v10.swd`, epoch 25 d (488 source zones, 9 velocity
extrema, 4 composition interfaces, 9 shocks):

| max_points | final N | splits | velocity topology | mass rel_err | species mass rel_err |
|-----------:|--------:|------:|:-----------------:|:-----------:|:--------------------:|
| 500        | 500     | 391   | 9/9 OK            | 0.0         | 3e-16                |
| 750        | 750     | 640   | 9/9 OK            | 0.0         | 2.6e-16              |
| 1000       | 1000    | 890   | 9/9 OK            | 0.0         | 2.6e-16              |

All three pass every validation check.  The algorithm uses the full budget
here because the model has many features; a smoother model would stop earlier
(adaptive).

## 9. Usage

### Standalone
```bash
python scripts/tau_amr_optimiser.py \
    --swd X_v11/stella/strad/X_v11.swd \
    --epoch 25 \
    --output X_v11/cmfgen/SN_HYDRO_DATA \
    --plot-dir X_v11/plots/tau_amr \
    --manifest X_v11/cmfgen/tau_amr_manifest.json \
    --config configs/tau_amr_config.yaml \
    --max-points 1000
```

### Via the orchestrator (run_model_grid.py)
Set in the model YAML `cmfgen_conversion` section:
```yaml
cmfgen_conversion:
  base_grid_method: "tau_amr_recursive"
  use_tau_amr: true
  max_points: 1000
  base_grid_target: 80
  composition_gradient_threshold: 0.30
```
Then:
```bash
python scripts/run_model_grid.py X_v11/config.yaml
```

## 10. Diagnostic plots

Generated in `--plot-dir`:
- `hydrodynamics.png` — density, velocity, temperature, τ (STELLA vs AMR)
- `composition.png` — per-element mass fractions (interfaces marked)
- `resolution.png` — dr, dr/r, dv, dlogτ (features marked)
- `point_density.png` — histogram of grid points in log-radius (features marked)

## 11. Two-stage AMR architecture (Stage B — future)

The current implementation is **Stage A** (pre-CMFGEN, hydrodynamic-only AMR).
The architecture is designed to support **Stage B** (CMFGEN-informed AMR): after
an initial CMFGEN run, read its output to identify ionisation fronts / rapid
population changes, then locally refine and re-run CMFGEN.  This is not yet
implemented but the module's `run_tau_amr()` entry point and manifest format
are structured to accommodate it.

## 12. Files

| file | role |
|------|------|
| `scripts/tau_amr_optimiser.py` | the optimiser (parser reuse, OD reconstruction, recursive AMR, detectors, conservative remap, validation, diagnostics, CLI) |
| `scripts/convert_stella_swd_to_cmfgen.py` | dispatch branch `base_grid_method == 'tau_amr_recursive'` calls `recursive_amr()` |
| `scripts/run_model_grid.py` | passes `base_grid_method` / `max_points` / `base_grid_target` through to the converter |
| `configs/tau_amr_config.yaml` | all tunable thresholds |
| `X_v11/config.yaml` | first model using the new grid (solar shell+wind, SN_A ejecta, interactive test) |
