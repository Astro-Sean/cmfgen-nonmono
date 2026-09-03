# STELLA → SuperLite Pipeline Audit

## 1. Current Implementation

### Two conversion paths

1. **Legacy** (`convert_sn_hydro_to_superlite.py`, 384 lines): CMFGEN `SN_HYDRO_DATA` → `input.str`. 9 grouped elements. No adaptive regridding. Geometric-mean cell boundaries from zone centres. Simple tau truncation using kappa from SN_HYDRO_DATA or grey 0.34 cm²/g.

2. **Preferred MC-aware** (`convert_stella_swd_to_superlite.py`, 1021 lines): STELLA `.swd` → `input.str` directly. 15 individual elements. Hybrid opacity (STELLA native + ES floor + Kramers bf/ff). Tau-based truncation. Tau-weighted downsampling.

### SuperLite source verification (VERIFIED)

SuperLite installed at `/ptmp/sbrennan/superlite/`, git commit `176d0a9` (r1.0.0 + bug fixes). Key physics findings from source inspection:

| Parameter | Source default | Project value | Meaning |
|-----------|---------------|---------------|---------|
| `in_trn_tauddmc` | 5.0 | 3.0 | MFP per cell for IMC→DDMC transition (total: abs+scat) |
| `in_taulump` | 10.0 | 3.0 | MFP for DDMC group lumping (absorption only) |
| `in_grp_ng` | -1 (file) | 6000 | Number of wavelength groups (logarithmic grid) |
| `in_src_n2s` | -1 | 25 (33M) | 2^n source particles |
| `in_nlte` | .false. | .true. | NLTE for H bound-bound only |
| `in_nlte_nelem` | 1 | 1 | Only H supported (hardcoded check) |
| `in_niter` | 1 | 3 | Radiation temperature feedback iterations |
| `in_nebular` | .false. | .true. | Mazzali & Lucy 1993 nebular approximation |
| `in_R_phot` | -1 | computed | Photospheric radius for dilution factor |

**Critical physics facts (VERIFIED from source):**
- SuperLite does **NOT** assume homologous expansion. Uses input `vx_left`/`vx_right` directly for Doppler shifts and velocity divergence (`GRID/grid_setup.f:20,48-51`). No `v=r/t` computation anywhere.
- Source emission is distributed as **T_gas^4** spatially, then **renormalised** to `in_L_bol` (`SOURCE/sourceenergy.f90:32,56`). Per-cell `lum` column is read but **disabled** (commented out).
- `rad_temp` defaults to `temp` if not provided (`GAS/gas_setup.f:35-36`).
- DDMC transition uses **total** (absorption+scattering) optical depth per cell (`TRANSPORT1/transport11.f90:333`).
- NLTE is restricted to **hydrogen only** (`inputparmod.f:411`: `if(in_nlte.and.in_nlte_nelem.ne.1) stop`).

## 2. Current Numerical Settings

### Hardcoded in `input_parameters_reference.par` (NOT in YAML)

| Setting | Value | Affects |
|---------|-------|---------|
| `in_grp_ng` | 6000 | Accuracy, runtime |
| `in_grp_wlmin` | 1d-8 (1 Å) | Accuracy |
| `in_grp_wlmax` | 30d-5 (30000 Å) | Accuracy |
| `in_flx_ndim` | 3500 | Accuracy, output |
| `in_flx_wlmin` | 3d-5 (3000 Å) | Output |
| `in_flx_wlmax` | 2d-4 (20000 Å) | Output |
| `in_src_n2s` | 25 (33M) | Accuracy, runtime |
| `in_trn_tauddmc` | 3.0 | Stability, accuracy |
| `in_taulump` | 3.0 | Stability, accuracy |
| `in_nlte` | true | Accuracy, runtime |
| `in_nlte_nelem` | 1 | Accuracy |
| `in_niter` | 3 | Accuracy, runtime |
| `in_nebular` | true | Accuracy |
| `in_nobbopac` | false | Accuracy |
| `in_nobfopac` | false | Accuracy |
| `in_noffopac` | false | Accuracy |
| `in_nothmson` | false | Accuracy |

### Exposed in YAML (SN_B_grid configs)

| Key | Default | SN_B_grid value |
|-----|---------|-----------------|
| `truncate_tau` | None→20 | **10000** (effectively no truncation) |
| `max_cells` | 150 | 250 |
| `input_source` | sn_hydro_data | swd |
| `t_floor_K` | 2000 | 2000 |
| `cpus_per_task` | 96 | 8 |

### NOT exposed in YAML at all (X1.yaml, SN_A have no superlite section)

All of the above hardcoded parameters plus: `min_cells`, `no_downsample`, `rho_floor`.

## 3. Identified Problems

### P1: Grid method is NOT adaptive (CRITICAL)

`tau_weighted_downsample` uses a **fixed** split:
- 40% of cells to "spectrum region" (tau 0–1), uniform in **zone index**
- 60% to "diffusion region" (tau 1–tau_max), uniform in **tau**

This is **not** adaptive mesh refinement. It does not respond to:
- density gradients
- temperature gradients
- velocity gradients / reversals / shocks
- composition gradients / interfaces
- the metal-rich cannon
- optical-depth gradients within each region

### P2: Severe grid non-uniformity (VERIFIED)

Production grid `M1.0_V10000/30.0days` (250 cells):
- `dr/r` ranges from **7.6e-5** to **0.85** — a factor of **11,000×**
- Giant gap at cell 1: `dr/r=0.85` (inner region, r=2.56e13 cm)
- Giant gap at cell 110: `dr/r=0.32` (r=6.11e15 cm, near cannon inner edge)
- The inner region (cells 0–50) has very few cells despite containing the dense shell

### P3: Metal-rich cannon not detected or preserved (CRITICAL)

The cannon (cells 112–153 in production grid):
- r = 9e15 – 2e16 cm
- v = 2150 – 4800 km/s
- mass = 0.65 Msun
- H-poor (H=0), Si/Fe/O-rich (Si up to 0.56, Fe up to 0.42, O up to 0.92)
- **42 cells** represent this 0.65 Msun component

The converter has **no logic** to detect or preserve this structure. The uniform-in-index downsampling can merge cannon cells with adjacent H-rich ejecta, smearing the composition interface.

### P4: Velocity reversals not resolved (CRITICAL)

Strong velocity reversal at r≈2.45e15 cm: v drops from 3306 → 127 km/s over ~10 cells. This is a physical shock/contact discontinuity. The downsampler has no velocity-gradient awareness and can smear it.

### P5: No conservation tests (CRITICAL)

Mass conservation is **printed** but **not enforced**. No element-mass conservation check. No composition normalisation check per cell. No velocity-structure preservation check. Failed conservation does not stop execution.

### P6: No convergence framework (CRITICAL)

No mechanism to test:
- grid convergence (compare spectra at N=50,100,150,200,300)
- wavelength convergence (compare at 500,1000,2000,3000,6000 groups)
- MC particle convergence (compare at 1e4,1e5,1e6,3e6)
- tau truncation convergence (compare at tau=10,20,50,100)

### P7: Parameters hardcoded (VIOLATES AGENTS.md §15)

All accuracy-affecting parameters (wavelength groups, particle count, DDMC thresholds, NLTE, nebular, flux range) are hardcoded in the template, not in YAML. This violates the requirement that "every value that materially affects spectral accuracy must be explicit."

### P8: truncate_tau=10000 with max_cells=250 (SUPPORTED)

SN_B_grid configs set `truncate_tau: 10000` (no truncation) but `max_cells: 250`. This means the full STELLA grid (~500 zones) is downsampled to 250 cells. The downsampler's uniform-in-index and uniform-in-tau split destroys physical structure. The combination is contradictory: keeping all optically thick material (tau=10000) but then aggressively downsampling.

### P9: No regression tests

No tests verify that the converter:
- preserves ordering
- creates no negative masses
- maintains composition normalisation
- preserves velocity extrema
- preserves composition discontinuities
- does not assume homologous expansion

## 4. Recommended Changes

### R1: New adaptive regridding module (`superlite_adaptive_grid.py`)

Replace `tau_weighted_downsample` with a physics-aware AMR algorithm that refines based on:
- density gradient: `d(log ρ)/d(log r)`
- temperature gradient: `d(log T)/d(log r)`
- velocity gradient: `dv/dr`, especially negative gradients (shocks)
- composition gradients: `dX_i/d(log r)` for H, He, O, Si, Ca, Fe, etc.
- optical depth: cumulative tau structure
- cannon detection: forced refinement at metal-rich/H-poor interfaces

Conservative remapping: mass-weighted averages for intensive quantities, exact mass conservation.

### R2: Cannon detection and preservation

Detect regions where Fe/Si/Ca/O increase sharply while H decreases. Force minimum cells at each interface. Report cannon diagnostics.

### R3: Conservation tests

Enforce (stop on failure in production):
- total mass: `|M_out - M_in|/M_in < tol`
- element masses: `|M_X,out - M_X,in|/M_X,in < tol` for each element
- composition: `|ΣX_i - 1| < tol` per cell
- velocity extrema preserved

### R4: YAML configuration layer

Expose ALL numerical parameters in `superlite:` YAML section, separated into physical/numerical/accuracy/debug groups.

### R5: Convergence framework (`superlite_convergence.py`)

Automated convergence testing with quantitative spectral comparison metrics.

### R6: Regression tests

Test suite covering: smooth ejecta, non-homologous velocity, shocks, dense shells, composition interfaces, cannon, optically thick/thin regions.
