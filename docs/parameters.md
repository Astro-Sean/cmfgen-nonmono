# CMFGEN VADAT Parameter Reference

This document explains every parameter in the VADAT control file used by
CMFGEN for the X_2 supernova model at 5.0 days
(`models/X_2/cmfgen/5.0days/VADAT`).

VADAT is the main "steering wheel" for a CMFGEN run. It tells the code:

1. **What star/supernova to model** (radius, mass, velocity law, composition)
2. **How to discretize the radiation field** (frequency grid, depth grid,
   line profiles)
3. **How to solve the statistical equilibrium equations** (linearization vs
   Lambda iteration, Ng acceleration, what to hold fixed)
4. **Which physical processes to include** (clumping, X-rays, two-photon,
   charge exchange, dielectronic transitions)

The file is read by `RD_CONTROL_VARIABLES` in
`cmfgen_source/new_main/mod_subs/rd_control_variables.f`. Each parameter is
written as `VALUE [NAME]`. The `[NAME]` in brackets is the keyword CMFGEN
searches for, so the order does not matter — only the keywords do.

This document is written for someone learning radiative transfer, so each
section starts with a short "what this is doing physically" explanation
before listing the parameters.

> **Depth convention reminder:** CMFGEN numbers depth from the **outside in**.
> Depth 1 is the outermost shell (lowest density, highest velocity); depth ND
> is the innermost shell (highest density). See `AGENTS.md` for details.

---

## 1. Stellar/SN Structure (lines 1–9)

These describe the **macroscopic structure** of the object being modelled:
its radius, mass, and how the velocity varies with radius.

| Line | Parameter  | Value           | Meaning |
|------|------------|-----------------|---------|
| 1    | `RSTAR`    | 39626.608       | Stellar/ejecta radius at the inner boundary, in units of 10^10 cm (~2.64 AU). For a SN this is the radius of the innermost shell of the ejecta at the chosen epoch. |
| 2    | `RMAX`     | 4.80958855      | Maximum radius in units of RSTAR. So the outer edge sits at 4.81 × RSTAR. The full grid runs from RSTAR to RMAX. |
| 3    | `VEL_LAW`  | 12              | Velocity-law type. Type 12 = "SN hydro model with a beta-law fallback" — the velocity structure is read from `SN_HYDRO_DATA`, but `VCORE` and `BETA1` are still used to fill in any missing pieces. |
| 4    | `VINF`     | 536.83          | Terminal velocity (km/s) of the beta-law component. For a SN this is the asymptotic velocity of the outer ejecta. |
| 5    | `VCORE`    | 196.5           | Core velocity (km/s). The velocity at the inner edge of the beta-law region. |
| 6    | `BETA1`    | 1.0             | Exponent of the beta-law: v(r) = VCORE + (VINF - VCOR) × (1 - RSTAR/r)^BETA1. BETA=1 is the standard linear-in-(1-1/r) law. |
| 7    | `LSTAR`    | 0.0D0           | Stellar luminosity in solar units. **Set to 0 for a SN** because the luminosity is supplied by the ejecta itself (radioactive decay + shock heating), not by a central engine. |
| 8    | `MDOT`     | 1.0D-20         | Mass-loss rate in Msun/yr. Effectively zero — SNe do not have a steady wind. |
| 9    | `MASS`     | 2.08025019D-01  | Ejecta mass in solar units (~0.208 Msun). Used to convert the density profile to a mass-consistent structure. |

### What "velocity law" means physically

In a SN, the ejecta is in **homologous expansion**: every shell moves
outward at a velocity proportional to its radius, `v ∝ r`. The beta-law is
a smooth analytic approximation to this. CMFGEN can either use the analytic
beta-law directly, or read a full hydrodynamical snapshot (which is what
`VEL_LAW=12` does here — the actual velocities come from `SN_HYDRO_DATA`,
with `VCORE`/`BETA1` as fallbacks).

---

## 2. Clumping (lines 10–12)

| Line | Parameter   | Value | Meaning |
|------|-------------|-------|---------|
| 10   | `DO_CL`     | F     | Do not include clumping (density inhomogeneities). |
| 11   | `CL_LAW`    | NONE  | Clumping law name — irrelevant since `DO_CL=F`. |
| 12   | `N_CL_PAR`  | 0     | Number of clumping parameters — irrelevant. |

Clumping multiplies the density locally by a factor 1/f (where f is the
filling factor) to mimic wind inhomogeneities. For an ejecta-only SN model
this is normally off.

---

## 3. Supernova-Specific Controls (lines 13–36)

These only matter when `VEL_LAW` is a SN type (10, 11, or 12). They control
how the SN is aged, whether radioactive decay is included, and how the
ejecta is heated.

| Line | Parameter                  | Value      | Meaning |
|------|----------------------------|------------|---------|
| 13   | `TS_NO`                    | 0          | Time-sequence number — used when chaining multiple SN snapshots. 0 means a standalone run. |
| 14   | `SN_AGE`                   | 4.902      | Age of the supernova in days at this snapshot. ~4.9 days post-explosion. |
| 15   | `SN_T_OPT`                 | USE_T_IN   | Tells CMFGEN to use the temperature already present in `SN_HYDRO_DATA` rather than recomputing it from a grey opacity model. |
| 16   | `SN_DEBUG`                 | T          | Print extra SN-specific diagnostic messages. |
| 17   | `PURE_HUB`                 | F          | Do NOT force pure Hubble flow (v ∝ r everywhere). Allows non-homologous regions (e.g. shell collisions) to be modelled. |
| 18   | `SN_DEBUG_LUM`             | T          | Print extra diagnostics about the SN luminosity budget. |
| 19   | `GAMRAY_TRANS`             | LOCAL      | How gamma-rays (from ^56Ni → ^56Co → ^56Fe decay) are transported. LOCAL = deposit the gamma-ray energy at the depth where it is emitted, no global transport. Cheaper and fine when the ejecta is optically thick to gamma-rays. |
| 20   | `TRT_NON_TE`               | F          | Do not include non-thermal electrons from gamma-ray deposition as a separate population. (When T, the code tracks a non-thermal electron distribution that can ionize H, He, etc.) |
| 21   | `INC_SHOCK_POWER`          | F          | Do not add an extra shock-heating luminosity source at the inner boundary. |
| 22   | `SCL_PWR_BY_FCL`           | F          | Irrelevant since `INC_SHOCK_POWER=F`. |
| 23   | `ADD_SHOCK_POWER_SLOWLY`   | F          | Irrelevant since `INC_SHOCK_POWER=F`. |
| 24   | `DO_HYDRO`                 | F          | Do NOT adjust the hydrostatic structure. The structure is fixed from the input snapshot. |
| 25   | `INC_AD`                   | T          | Include the **advection term** in the radiative transfer equation. This is the term that accounts for the fact that the gas is moving — radiation is advected along with the flow. Important for SNe. |
| 26   | `INC_ADV`                  | F          | Do not include advection in the moment equations separately (already handled by `INC_AD`). |
| 27   | `DO_DDT`                   | F          | Do not include explicit time-derivative (∂J/∂t) terms in the transfer equation. The model is treated as a sequence of steady-state snapshots. |
| 28   | `DO_CO_MOV_DDT`            | F          | Do not use the co-moving time-derivative formulation. |
| 29   | `INC_RAD_DECAYS`           | T          | Include radioactive decay heating from ^56Ni, ^56Co, etc. This is the dominant luminosity source for many SNe at this epoch. Auto-set by `prepare_VADAT.py`: T when `SN_HYDRO_DATA` has isotope mass-fraction blocks, F otherwise. The VADAT_INIT template default is F; `prepare_VADAT.py` overrides to T when isotopes are present (line 2578). |
| 30   | `INS_DEP`                  | F          | Do not treat the gamma-ray energy deposition as instantaneous. Use the steady-state deposition approximation instead. |
| 31   | `RD_CF_FILE`               | F          | Do not read the continuum frequency grid from a file — let CMFGEN compute it. |
| 32   | `MIN_CF`                   | 1.0D-2     | Minimum continuum frequency (in some internal unit). Sets the reddest wavelength CMFGEN will treat as "continuum". |
| 33   | `MAX_CF`                   | 25.0D0     | Maximum continuum frequency. Sets the bluest continuum frequency. |
| 34   | `FRAC_SP`                  | 1.08D0     | Fractional spacing between adjacent frequency points where the spacing is small (near line centers). 1.08 means each step is 8% larger than the previous. |
| 35   | `AMP_FAC`                  | 1.04D0     | Amplification factor for the frequency grid in regions where the spacing is large (line wings / continuum). 1.04 = 4% growth per step. |
| 36   | `MAX_BF`                   | 0.05D0     | Maximum frequency spacing (as a fraction of the edge frequency) close to a bound-free edge. Keeps the grid fine enough to resolve the edge. |

### Why frequency grid parameters matter

CMFGEN solves the transfer equation on a **frequency grid**. Lines are very
narrow (a few km/s wide) but the spectrum spans a huge range, so the grid
must be **adaptive**: fine near line centers and bound-free edges, coarse in
the continuum. `FRAC_SP`, `AMP_FAC`, and `MAX_BF` control that adaptivity.

---

## 4. Level Dissolution & Photoionization (lines 37–45)

| Line | Parameter      | Value    | Meaning |
|------|----------------|----------|---------|
| 37   | `DO_DIS`       | F        | Do not include **level dissolution** (the merging of high-lying bound levels into the continuum when densities are high). |
| 38   | `dV_LEV`       | 200.0D0  | Velocity spacing (km/s) on the low-frequency side of a bound-free edge where dissolution would be evaluated. |
| 39   | `AMP_DIS`      | 1.4D0    | Amplification factor for the dissolution frequency grid. |
| 40   | `MIN_DIS`      | 0.1D0    | Minimum frequency at which dissolution is computed. |
| 41   | `CROSS`        | F        | Do not recompute all photoionization cross-sections at every depth — use interpolated values. |
| 42   | `V_CROSS`      | 500.0D0  | Maximum velocity separation (km/s) between full cross-section evaluations. |
| 43   | `SIG_GAU_KMS`  | 500.0D0  | Sigma (km/s) of the Gaussian used to smooth noisy photoionization cross-section data. |
| 44   | `EXT_LINE_VAR` | 0.5D0    | How far (in units of VINF) beyond the resonance zone a line is allowed to vary. 0.5 = half of VINF. |
| 45   | `ZNET_VAR_LIM` | 0.01D0   | Include a line in the full linearization only when |ZNET - 1| > 0.01, where ZNET is the net rate ratio. Small ZNET means the line is near equilibrium and can be treated more simply. |

### What "level dissolution" is

In dense gas, the electric field of neighbouring ions shifts and broadens
the energy levels of an atom. Levels close to the ionization limit get
pushed into the continuum and effectively disappear — their population
merges into the free-electron pool. This is **level dissolution** (also
called "continuum lowering"). It matters at high densities (early SN
epochs); it is turned off here.

---

## 5. Line Treatment (lines 46–47)

| Line | Parameter | Value   | Meaning |
|------|-----------|---------|---------|
| 46   | `WNET`    | T       | Use Lambda iteration (a simple iterative scheme) for weak lines instead of full linearization. Weak lines have a small effect on the populations, so a cheap treatment is fine. |
| 47   | `WK_LIM`  | 0.10D0  | A line counts as "weak" if its opacity is less than 10% of the continuum opacity at its wavelength. |

### Strong vs weak lines

A **strong line** has an optical depth comparable to or larger than the
continuum. It strongly affects the radiation field and the level
populations, so it must be solved **simultaneously** with the populations
(via linearization). A **weak line** has a small optical depth and only
weakly feeds back on the populations, so it can be solved with a simpler
**Lambda iteration** (fix the populations, solve the line, update the
radiation, repeat). This is much cheaper.

---

## 6. Relativistic & Time-Dependent Transfer (lines 48–57)

| Line | Parameter        | Value | Meaning |
|------|------------------|-------|---------|
| 48   | `USE_J_REL`      | F     | Do not use the relativistic moment solver (`MOM_JREL_V9`). Disabled because MOM_JREL_V9 caps sigma (dlnv/dlnr) at -1.0, which produces NaN luminosity and corrupts the R grid when the velocity profile is non-monotonic (sigma < -1 at velocity reversals, e.g. shell collisions). See AGENTS.md for full history. |
| 49   | `USE_JREL`       | F     | Duplicate/legacy spelling of the above — also off. |
| 50   | `USE_LAM_ES`     | T     | **J=S formal solution**: set J directly from the source function (TC), bypassing the moment equations entirely (`comp_j_blank.f:727`: `RJ(1:ND)=TC(1:ND)`). This is the preferred approach for non-monotonic velocity — it avoids MOM_JREL_V9's sigma cap at -1.0 that causes NaN. With USE_LAM_ES=T, J follows the source function smoothly across velocity reversals, avoiding the sharp formal-solution features that cause population oscillations. |
| 51   | `USE_FRM_REL`    | F     | Do not use the relativistic formal solver (`CMF_FORMAL_REL`) for the flux. Disabled because USE_J_REL=F (CMFGEN's `check_param_consistency.f` requires USE_J_REL=T for relativistic terms). |
| 52   | `INCL_REL`       | F     | Do not include relativistic terms (aberrations, Doppler shifts from the flow) in the transfer equation. Disabled because USE_J_REL=F (CMFGEN requires USE_J_REL=T if INCL_REL=T). |
| 53   | `INCL_ADV_TRANS` | F     | Do not include advection terms in the **transfer equation** (as opposed to the moment equations — see `INC_AD` above). |
| 54   | `INCL_DJDT`      | F     | Do not include ∂J/∂t terms in the transfer equation. |
| 55   | `USE_DJDT_RTE`   | F     | Do not use the time-dependent transfer solver. |
| 56   | `JC_W_EDD`       | T     | Compute the continuum mean intensity J using **Eddington factors** (a closure relation that relates J to the flux H). Cheaper and more stable than a full moment solve. |
| 57   | `JBAR_W_EDD`     | T     | Same as above but for the line (Sobolev) mean intensity Jbar. |

### Eddington factors — a key RT trick

The radiation field has three angle moments:
- J = mean intensity (zeroth moment)
- H = flux (first moment)
- K = radiation pressure (second moment)

The transfer equation gives you a hierarchy: each moment depends on the
next one up. To close it, you need a relation between K and J. The
**Eddington factor** f = K/J is that closure. CMFGEN computes f from a
formal solution along rays, then uses it in the cheaper moment equations.
This is much faster than solving the full angle-dependent transfer
equation at every iteration.

---

## 7. Boundary Conditions & Method (lines 58–67)

| Line | Parameter  | Value     | Meaning |
|------|------------|-----------|---------|
| 58   | `DIF`      | T         | Use the **diffusion approximation** at the inner boundary. Physically: deep inside the ejecta, radiation is trapped and diffuses outward — its intensity is nearly isotropic and proportional to the local Planck function. |
| 59   | `IB_METH`  | DIFFUSION | Inner boundary method. Must match `DIF=T`. |
| 60   | `COH_ES`   | T         | Treat electron scattering as **coherent** (no frequency shift in the electron's rest frame). Valid when kT_e << m_e c^2. |
| 61   | `N_TYPE`   | G_ONLY    | Eddington factor type. G_ONLY = use the continuum (grey) Eddington factor (N/H ratio). This is safe with USE_LAM_ES=T because the moment equations are bypassed (J comes from the formal solution, not from moments). N_ON_J (N/J ratio) is only for USE_J_REL=T (MOM_JREL_V9), which is disabled for this model. Note: G_ONLY uses N/H which is undefined when H=0 at velocity reversal, but this is not a problem with USE_LAM_ES=T since the moment solver is not called. |
| 62   | `FG_OPT`   | INT/INS   | Solution option for the formal solver `FG_J_CMF`. INT/INS = use integration with insertion of extra depth points where needed. |
| 63   | `CHK_J`    | ABS_VAL   | How to check the consistency of J. ABS_VAL = check the absolute value of the residual. |
| 64   | `CHK_XM`   | SET_POS   | Check the moment XM (mass column) and force it to be positive. |
| 65   | `CHK_H`    | MAX_VAL   | Check the maximum value of H (flux) for consistency. |
| 66   | `CHK_NG`   | F         | Do not perform a validity check on Ng-accelerated corrections. |
| 67   | `METHOD`   | LOGMON    | Method for computing the continuum optical depth: LOGMON = log-monotonic interpolation. |
| 68   | `LUM_METH` | ZERO      | Method for computing the luminosity from the emissivity eta. ZERO = use a simple (zeroth-order) approach. |

### Boundary conditions in RT

At the **inner boundary** (deep in the ejecta), the gas is so optically
thick that radiation is in thermal equilibrium with the matter — this is
the **diffusion approximation**, J ≈ B(T). At the **outer boundary** (the
photosphere / surface), no radiation comes in from outside, so the
incoming intensity is zero. CMFGEN has several variants of each; this run
uses the simplest physically reasonable choices.

---

## 8. Velocity Grid for Transfer (lines 69–70)

| Line | Parameter   | Value      | Meaning |
|------|-------------|------------|---------|
| 69   | `VFRAC_FG`  | 20000.0D0  | Maximum velocity spacing (in Doppler units) in the formal solver `FG_J_CMF`. Larger = coarser but faster. |
| 70   | `VFRAC_MOM` | 20000.0D0  | Same but for the moment solver `MOM_J_CMF`. |

These are very large, meaning the depth grid is used as-is for the
radiation solve without inserting extra points. This is fine for a
coarse SN model.

---

## 9. Continuum & Eddington Factor Options (lines 71–77)

| Line | Parameter  | Value    | Meaning |
|------|------------|----------|---------|
| 71   | `THK_CONT` | F        | Do not use the "thick continuum" boundary condition (which would treat the continuum as optically thick everywhere). |
| 72   | `TRAP_J`   | T        | Use trapezoidal weights when integrating over angles to compute J. More accurate than simple midpoint. |
| 73   | `COMP_F`   | T        | Compute new Eddington factors f at each iteration (rather than reusing old ones). |
| 74   | `ACC_F`    | 1.0D-02  | Target accuracy for the Eddington factor: 1% . |
| 75   | `N_IB_INS` | 15       | Number of extra grid points to insert near the inner boundary for the formal solver. |
| 76   | `N_OB_INS` | 30       | Number of extra grid points to insert near the outer boundary. |
| 77   | `REXT_FAC` | 0.0D0    | Do not extend the radius grid beyond RMAX for the thick-continuum treatment. |

---

## 10. Line Profile & Doppler Parameters (lines 78–91)

These control the **line profile shape** — how the opacity of a spectral
line is distributed in frequency. In a moving atmosphere, the profile is
set by thermal broadening (Doppler) and any turbulent velocity.

| Line | Parameter        | Value    | Meaning |
|------|------------------|----------|---------|
| 78   | `FIX_DOP`        | T        | Use the same Doppler width for all species (rather than computing it per ion). |
| 79   | `TDOP`           | 4.0D0    | Temperature (in 10^4 K, so 40000 K) used to compute the thermal Doppler width. |
| 80   | `AMASS_DOP`      | 1.0D+40 | Atomic mass (in amu) for Doppler width. Manual: "set AMASS_DOP=1.0D+40, use VTURB to control width." This makes the thermal Doppler contribution negligible, so VTURB=20 km/s solely controls line broadening. |
| 81   | `VTURB`          | 20.0D0   | Turbulent velocity (km/s) added in quadrature to the thermal velocity to broaden lines. |
| 82   | `MAX_DOP`        | 8.0D0    | Maximum Doppler width (in some internal unit) — caps the broadening. |
| 83   | `FRAC_DOP`       | 1.5D0    | Fractional spacing across the Doppler profile. 1.5 = each frequency step is 50% larger than the previous. |
| 84   | `dV_CMF_PROF`    | 80.0D0   | Velocity spacing (km/s) of the frequency grid in the **core** of a line, in the co-moving frame. |
| 85   | `dV_CMF_WING`    | 200.0D0  | Velocity spacing in the **wings** of a line, co-moving frame. Wings are broader so coarser spacing is OK. |
| 86   | `ES_WING_EXT`    | 12000.0D0| Extent (km/s) of the electron-scattering wing. e-scattering broadens lines via redistribution. |
| 87   | `R_CMF_WING_EXT` | 4.0D0    | Radius (in RSTAR) out to which the CMF wing extension is applied. |
| 88   | `OBS_EXT_RAT`    | 1.3D0    | Ratio by which to extend the line profile in the observer-frame grid. |
| 89   | `dV_OBS_PROF`    | 40.0D0   | Velocity spacing (km/s) in the line core, observer frame. |
| 90   | `dV_OBS_WING`    | 150.0D0  | Velocity spacing in the line wings, observer frame. |
| 91   | `dV_OBS_BIG`     | 800.0D0  | Velocity spacing for very broad line components, observer frame. |

### Doppler broadening — the basics

Atoms in a gas are moving randomly due to their thermal energy. An atom
moving toward you sees the photon blue-shifted; one moving away sees it
red-shifted. The net effect is that a spectral line — which would be a
single frequency in a stationary atom — gets **broadened** into a profile
whose width is set by the thermal velocity v_th = sqrt(2kT/m). Adding a
turbulent velocity v_turb in quadrature gives a wider "effective"
temperature. The line shape is approximately Gaussian with width
sqrt(v_th^2 + v_turb^2).

---

## 11. Formal Solution & Wavelength Range (lines 92–105)

| Line | Parameter          | Value    | Meaning |
|------|--------------------|----------|---------|
| 92   | `FLUX_CAL_ONLY`    | F        | Do a full population solution, not just a flux calculation. (When T, CMFGEN skips the population iteration and only computes the emergent spectrum from fixed populations.) |
| 93   | `EXT_FRM_SOL`      | F        | Do not extend the formal solution beyond the default range. |
| 94   | `INS_F_FRM_SOL`    | T        | Insert extra frequency points in the formal solution where needed for accuracy. |
| 95   | `FRM_OPT`          | INT      | Formal solution option: INT = integration (as opposed to a ray-by-ray method). |
| 96   | `DO_SOB_LINES`     | T        | Include Sobolev (line) transfer. The Sobolev approximation treats a line as interacting with the radiation field only at the resonance point where the Doppler shift matches the line frequency. |
| 97   | `SOB_FREQ_IN_OBS`  | T        | Evaluate the Sobolev line frequencies in the observer frame. |
| 98   | `LAM_SET`          | F        | Do not set the transition type by wavelength. |
| 99   | `F_LAM_BEG`        | 50.0D0   | Starting wavelength (Å) for the flux calculation. 50 Å = extreme UV. |
| 100  | `F_LAM_END`        | 5.0D4    | Ending wavelength (Å) for the flux calculation. 50000 Å = far IR. |
| 101  | `GLOBAL_LINE`      | SOB      | Treat all lines globally as Sobolev lines (as opposed to treating some with a full profile). |
| 102  | `GF_CUT`           | 5.0D-3   | Cutoff on the oscillator strength × abundance (gf value): lines weaker than this are dropped. |
| 103  | `GF_LEV_CUT`       | 9        | Cutoff on the level index: only consider lines originating from the first 9 levels of each ion. |
| 104  | `AT_CUT`           | 22       | Cutoff on the atomic number: only include lines from elements with Z <= 22 (Ti). Heavier elements are ignored. |
| 105  | `MIN_TRANS`        | 15       | Minimum number of transitions to include per ion, regardless of the gf cut. Ensures weak ions still get some lines. |

### The Sobolev approximation

In a rapidly expanding atmosphere (like a SN), each point in the flow only
"sees" a line at one specific radius — the **resonance point** — where the
Doppler shift from the velocity field makes the local line frequency match
the photon frequency. Everywhere else, the line is transparent. This lets
you replace an integral over the line profile with a single point
evaluation: this is the **Sobolev approximation**. It is hugely cheaper
than a full line profile calculation and is very accurate when the
expansion velocity is much larger than the thermal velocity (v_exp >> v_th),
which is true in SNe.

---

## 12. Line Overlap & Opacity Checks (lines 106–114)

| Line | Parameter       | Value     | Meaning |
|------|-----------------|-----------|---------|
| 106  | `THK_LINE`      | F         | Do not treat lines as optically thick. |
| 107  | `CHK_L_POS`     | T         | Check that the line opacity is positive; fix it if not. |
| 108  | `NEG_OPAC_OPT`  | SRCE_CHK  | If negative opacity is found, check the source function and reset it. |
| 109  | `He2_RES=0`     | F         | Do not force the He II resonance lines to zero (i.e. include them normally). |
| 110  | `ALLOW_OL`      | T         | Allow **line overlap** — multiple lines can overlap in frequency and be treated together. |
| 111  | `OL_DIF`        | 50.0D0    | Maximum frequency difference (km/s) for two lines to count as overlapping. |
| 112  | `SCL_LN`        | T         | Scale the line cooling rates to prevent them from diverging when populations are far from convergence. |
| 113  | `SCL_LN_FAC`    | 0.1D0     | Factor by which to scale line cooling rates (10% of the unscaled value). |
| 114  | `SCL_SL_OPAC`   | F         | Do not scale the Sobolev line opacity separately. |

### Line overlap

In a dense spectrum, many lines can overlap in frequency — a photon
emitted by one line can be absorbed by another. CMFGEN handles this by
grouping overlapping lines into "clusters" and solving them together.
`ALLOW_OL=T` turns this on; `OL_DIF` sets the tolerance for what counts
as overlapping.

---

## 13. Electron Scattering & Additional Processes (lines 115–122)

| Line | Parameter   | Value | Meaning |
|------|-------------|-------|---------|
| 115  | `OLD_J`     | F     | Do not use a previous J file to initialize the electron-scattering source function. |
| 116  | `MIX_COH`   | F     | Do not mix coherent and non-coherent electron scattering in the linearization. |
| 117  | `ES_FAC`    | 0.6D0 | Fractional proximity between the radiation field J and the electron-scattering J_ES for the coherent variation. 0.6 is a typical stability choice. |
| 118  | `INC_CHG`   | F     | Do not include charge-exchange reactions (e.g. H + He+ ↔ H+ + He). |
| 119  | `INC_TWO`   | T     | Include **two-photon transitions** (e.g. the 2s → 1s decay in H I, which emits a continuous spectrum). |
| 120  | `INC_PEN`   | F     | Do not include Penning ionization (an excited atom collides with another and ionizes it). |
| 121  | `INC_RAY`   | F     | Do not include Rayleigh scattering. |
| 122  | `INC_XRAYS` | F     | Do not include a separate X-ray emission component. |

### Two-photon emission

Some transitions are forbidden as single-photon processes (e.g. 2s → 1s
in hydrogen-like atoms, because both levels have the same parity). They
instead proceed by **simultaneous emission of two photons** whose
energies sum to the level spacing. This produces a broad continuum bump
in the UV, important for the spectra of H-rich SNe near maximum light.

---

## 14. Radius Grid & Initial Estimates (lines 123–131)

| Line | Parameter    | Value    | Meaning |
|------|--------------|----------|---------|
| 123  | `RD_IN_R_GRID` | F      | Do not read a predetermined radius grid from a file — use the grid from `SN_HYDRO_DATA`. |
| 124  | `LIN_INT`    | F        | Do not use direct linear interpolation for a new model. |
| 125  | `POP_SCALE`  | T        | Scale the populations so that the conservation equations (sum of ion fractions = 1, charge neutrality) are exactly satisfied after each iteration. |
| 126  | `IT_ON_T`    | T        | Improve the initial temperature estimate by iterating on it before the main population solve. |
| 127  | `T_INIT_TAU` | 20.0D0   | Optical depth above which the initial T is set exactly to the spherical (grey) T. 20 = deep in the ejecta. |
| 128  | `GREY_TAU`   | 10.0D0   | Optical depth above which T is set to the grey (Eddington) temperature during the iterative process. |
| 129  | `LTE_EST`    | T        | Use **LTE** (local thermodynamic equilibrium) populations as the initial estimate. Good for a cold start. **Note (2026-08-03):** `LTE_EST=T` is now used in all configs. The alternative `proxy_populations` approach (LTE_EST=F + custom _IN files) was empirically shown to be worse for cold snapshots — see AGENTS.md "proxy_populations vs CMFGEN built-in LTE_EST". |
| 130  | `T_EXC`      | 3.0D-01  | Minimum temperature (in 10^4 K, so 3000 K) used when computing LTE populations. Prevents absurdly low-T divergences. |
| 131  | `T_MIN`      | 3.0D-01  | Minimum allowed electron temperature anywhere in the model (3000 K). |

### LTE vs non-LTE

In **LTE** (Local Thermodynamic Equilibrium), every level population is
given by the Saha-Boltzmann formula at the local T. This is wrong in a SN
atmosphere (the radiation field is not a blackbody, so the populations are
not at their LTE values), but it is a useful **starting guess** for the
iterative solver. CMFGEN then iterates to the correct non-LTE
populations. `LTE_EST=T` just controls the initial guess, not the final
answer.

---

## 15. Composition (lines 132–139)

These are the **mass fractions** of each element, expressed as logarithms
(base 10) of the number fraction relative to hydrogen. So `HYD/X = 0`
would mean X_H = 1 (pure H); a negative value means less than that.

| Line | Parameter  | Value            | Element         |
|------|------------|------------------|-----------------|
| 132  | `HYD/X`    | -2.703159D-02    | Hydrogen (log of H/H = 0 by definition; the small negative value is a numerical offset). |
| 133  | `HE/X`     | -5.14D-01        | Helium. log10(He/H) ≈ -0.514 → He/H ≈ 0.31. |
| 134  | `OXY/X`    | -2.43D-01        | Oxygen. log10(O/H) ≈ -0.243 → O/H ≈ 0.57. |
| 135  | `SIL/X`    | -9.4595584D-02   | Silicon. log10(Si/H) ≈ -0.095 → Si/H ≈ 0.80. |
| 136  | `CAL/X`    | -2.70D-02        | Calcium. |
| 137  | `IRON/X`   | -9.45965845D-02  | Iron. |
| 138  | `NICK/X`   | -9.99999984D-09  | Nickel. Effectively zero (the Ni is created by the decay chain, not present initially). |
| 139  | `COB/X`    | -9.9999998D-09   | Cobalt. Effectively zero. |

These abundances are overridden by the `SN_HYDRO_DATA` file (which carries
the full composition profile from the hydrodynamics code), so the values
here are mostly fallbacks.

---

## 16. Dielectronic Transitions (lines 140–173)

**Dielectronic recombination** is a process where a free electron
recombines with an ion, but instead of emitting a single photon, the
energy goes into exciting a bound electron of the recombined ion. The
ion then stabilizes by emitting one or more photons. The reverse
(dielectronic capture) is an autoionization pathway.

| Line | Parameter      | Value    | Meaning |
|------|----------------|----------|---------|
| 140  | `DIE_AS_LINE`  | F        | Do not treat each dielectronic transition as an individual line. (When T, each is solved as a Sobolev line; when F, they are folded into the recombination rates.) |
| 141  | `VSM_DIE`      | 1000.0D0 | Velocity (km/s) for smoothing dielectronic transitions when treated as lines. |
| 142–173 | `DIE_<ion>` | F,F      | For each ion (HI, HeI, He2, Ni, NII, NIII, NV, OI, ..., FeSEV), two flags: (1) automatically compute dielectronic data, (2) compute via the WI (Weatherall & Ireland) method. All set to F,F = do not include dielectronic transitions for any ion. |

The list of ions covered is:
HI, HeI, He2, Ni, N2, NIII, NV, OI, O2, OIII, OIV, OV, OSIX,
Sk2, SkIII, SkIV, CaI, Ca2, PIII, PIV, PV, PPIX,
SIII, SIV, SV, SSIX, Fe2, FeIII, FeIV, FeV, FeSIX, FeSEV.

---

## 17. Per-Ion Level Fixing (lines 174–211)

For each ionization stage of each element, you can choose to **fix** a
certain number of levels — i.e. hold their populations at the initial
estimate and not let the solver change them. This is useful when:

- You are adding a new ion and want to converge the rest of the model first
- A particular ion is causing instability
- You want a cheap exploratory run

All values here are `0`, meaning **no levels are fixed** — every level is
free to converge.

| Lines | Parameters | Value | Meaning |
|-------|-----------|-------|---------|
| 174–210 | `FIX_HI`, `FIX_HYD`, `FIX_HeI`, `FIX_He2`, `FIX_HE`, `FIX_NI`, `FIX_N2`, ..., `FIX_IRON` | 0 | Number of levels to fix for each ion / species. 0 = none. |
| 211 | `FIX_NE` | F | Do not fix the electron density — let it vary with the ionization balance. |
| 212 | `FIX_IMP` | F | Do not **automatically fix impurity species**. (See note below.) |
| 213 | `FIX_T` | T | **Fix the temperature** at the initial estimate — do not iterate on T. This is common for a single-snapshot run where the T structure is taken from the hydro code. |
| 214 | `FIX_T_AUTO` | F | Do not let CMFGEN automatically decide when to fix/unfix T. |
| 215 | `FIX_INB_T` | F | Do not fix T only at the inner boundary. |
| 216 | `FIX_X_DPTH` | 0 | Do not fix T at the last N depths. |
| 217 | `TAU_SCL_T` | 0.0D0 | Electron-scattering optical depth above which T is held fixed. 0 = use the `FIX_T` flag globally. |

### About `FIX_IMP`

`FIX_IMP=T` tells CMFGEN to automatically "fix" any ion whose abundance is
so low that it cannot affect the electron density or the radiation field —
a so-called **impurity species**. This is normally a safe optimization.
However, in this project it was set to `F` because the automatic
impurity-fixing logic was found to spuriously zero out corrections during
the Lambda-iteration phase of a cold start, causing the solver to declare
false convergence and then crash with `DTDR=NaN` on the next iteration.
See `updates2cmfgen.md` and `AGENTS.md` for the full history.

---

## 18. Solver Method & Convergence (lines 218–247)

These control the **iterative scheme** that solves the statistical
equilibrium equations.

| Line | Parameter      | Value    | Meaning |
|------|----------------|----------|---------|
| 218  | `SOL_METH`     | DIAG     | Matrix solution method. DIAG = diagonal (no inter-depth coupling in the linearization). Cheapest; appropriate when NUM_BNDS=1. Other options: TRI (tridiagonal), PEN (pentadiagonal), GSIT (Gauss-Seidel), MIN (minimum residual). |
| 219  | `SCALE_OPT`    | MAJOR    | Population scaling option. MAJOR = scale each depth independently based on major species only (pop > 1e-10 * Ne), ignoring trace/impurity species that could dominate the scale factor. This is depth-independent (no inter-depth coupling). Other options: LOCAL (all species), NONE (no scaling), GLOBAL (single scale for all depths — introduces inter-depth coupling, not recommended). |
| 220  | `EPS_TERM`     | 5.0D-02  | Convergence criterion: if the maximum fractional population change drops below 5%, declare convergence and stop. Manual: "Typically choose 0.01 to 0.1." |
| 221  | `MAX_LIN`      | 2.0D0    | Maximum fractional change allowed per **linearization** iteration. 2.0 = 200% change. Manual: "Typically adopt 10"; 2.0 is conservative for SN. |
| 222  | `MAX_LAM`      | 1.0D2    | Maximum fractional change allowed per **Lambda iteration**. 100x max change. Manual: "Typically adopt 10^5; up to 10^20 for new species." 1.0D2 is safe for LTE cold start (1e5 caused BA singularity from zero pops). Previous 1.1D0 (10% cap) caused infinite Lambda loop — trace species needing 10000% correction could only change 10% per iteration. |
| 223  | `MAX_CHNG`     | 1.0D100  | If the maximum fractional change ever exceeds this, abort. 10^100 = effectively no limit. Manual: "Typically adopt 1.0×10^40." |
| 224  | `COMP_BA`      | T        | Compute the BA (rate-equation) matrix. The BA matrix is the linearized system that gets solved each iteration. |
| 225  | `STORE_BA_INV` | T        | Store the inverse of the BA matrix between iterations (saves recomputation). |
| 226  | `WR_PRT_INV`   | T        | Store only part of the inverse (for tridiagonal storage savings). |
| 227  | `STORE_BA`     | T        | Store the BA matrix itself between iterations. |
| 228  | `N_FIX_BA`     | 5        | Number of iterations to hold the BA matrix fixed once it is converged enough. Manual: "recommended value is 2." Previous 30 kept BA frozen too long. |
| 229  | `BA_CHK_FAC`   | 5.0D-02  | If the change in J is less than 5%, ignore the correction to the BA matrix. Manual: "Typically we adopt 5%." Previous 1.0D-05 was too tight — BA was never recomputed. |
| 230  | `FIX_BA`       | 1.0D-02  | If the maximum fractional population change drops below 1%, stop recomputing the BA matrix (hold it fixed). |
| 231  | `LAM_VAL`      | 4.0D2    | If the maximum fractional population change exceeds 400, switch to Lambda iteration. Manual: "Typically do a Lambda-iteration if % change is greater than 400." Previous 200 triggered Lambda too early. |
| 232  | `NUM_LAM`      | 2        | Number of Lambda iterations to do when the above triggers. Manual: "Typically choose 2." Previous 10 kept CMFGEN in Lambda-only mode too long between FULL(FIXT) blocks. |
| 233  | `RD_SOL`       | F        | Do not read in a previous solution vector to update populations. |
| 234  | `INC_GRID`     | T        | Increase the depth-grid size to improve accuracy where needed. |
| 235  | `ALL_FREQ`     | T        | Increase accuracy for all frequencies (not just a few key ones). |
| 236  | `ACC_END`      | 5.0D0    | Increase accuracy for all frequencies below this value (in some internal unit). |
| 237  | `N_INS`        | 2        | Number of points to insert in the higher-accuracy grid (1, 2, or 3). |
| 238  | `ST_INT`       | 1        | Interpolate from depth index 1 (outermost). |
| 239  | `END_INT`      | 30       | Interpolate to depth index 30. |
| 240  | `ND_QUAD`      | 30       | Use quadratic interpolation from depth ND-30 to ND (innermost 30 points). |
| 241  | `INTERP_TYPE`  | LOG      | Perform interpolations in log space (better for quantities that span many orders of magnitude). |
| 242  | `N_PAR`        | 150      | Rate of BA incrementation by BA_PAR in the continuum loop (number of frequencies). |
| 243  | `DO_NG`        | T        | **Perform Ng acceleration.** Manual: "Ng acceleration is an important tool for obtaining convergence" especially with NUM_BNDS=1/DIAG. Previous F caused infinite Lambda loop (549 iterations, no convergence). |
| 244  | `BEG_NG`       | 1.0D1    | Percentage accuracy at which to begin Ng acceleration. Manual: "Typically choose 5% to 10%." Previous 1.0D15 effectively disabled Ng (never reached 10^15 % accuracy). |
| 245  | `IBEG_NG`      | 15       | Iteration at which to begin Ng acceleration. Manual: "at least 15 (or more) iterations have been completed." |
| 246  | `BW_NG`        | 30       | Depth band width for Ng acceleration. |
| 247  | `ITS/NG`       | 8        | Number of iterations between successive Ng accelerations. Manual: "We typically adopt 8." |

### Linearization vs Lambda iteration

There are two main ways to solve the non-LTE equations:

1. **Linearization** (Newton's method): linearize the equations around the
   current populations, solve the linear system for the correction, apply
   it, repeat. Fast convergence (quadratic) but each step is expensive and
   can be unstable if the initial guess is bad.

2. **Lambda iteration**: fix the radiation field, solve for the
   populations; then fix the populations, solve for the radiation field;
   repeat. Very stable but very slow convergence (linear).

CMFGEN uses linearization as the main engine, but falls back to Lambda
iteration when the changes are too large (`LAM_VAL`). `MAX_LIN` and
`MAX_LAM` cap the size of each correction to prevent the solver from
overshooting.

### Ng acceleration

**Ng acceleration** is an extrapolation technique: once you have a
sequence of iterates that are converging smoothly, you can extrapolate
ahead to skip many iterations. It is very powerful but can be unstable —
if the iterates are not yet in the smooth-convergence regime, the
extrapolation can blow up. The CMFGEN manual states: "Ng acceleration is
an important tool for obtaining convergence" especially with NUM_BNDS=1
(DIAG solver), where "convergence may not be as clearcut as with
NUM_BNDS=3."

Previously `DO_NG=F` was chosen due to NaN propagation in the cold-start
phase. However, this caused an infinite Lambda loop (549 iterations with
10000% corrections never converging) because without Ng, pure Lambda
converges linearly and takes thousands of iterations for IIn models with
extreme gradients. The current setup enables Ng (`DO_NG=T`) with
`IBEG_NG=15` (manual: "at least 15 iterations") and `BEG_NG=10%`
(manual: "5% to 10%"), which allows Ng to activate only after the initial
Lambda settling phase. See `updates2cmfgen.md` for the full history.

---

## 19. Reference Frame & Final Options (lines 248–259)

| Line | Parameter      | Value | Meaning |
|------|----------------|-------|---------|
| 248  | `REL_OBS`      | T     | Use the full relativistic observer-frame transformation when computing the emergent spectrum. |
| 249  | `REL_CMF`      | T     | Use the full relativistic co-moving-frame transformation for the internal transfer. |
| 250  | `AUTO_ADD`     | F     | Do not automatically add missing ionization-stage files. |
| 251  | `DO_LAM_IT`    | T     | Do a Lambda iteration (i.e. a radiation-field-only update) at the end of each main iteration. Helps smooth out oscillations. |
| 252  | `DIE_SkI`      | F,F   | Dielectronic flags for Sk I (same two-flag format as section 16). |
| 253  | `DIE_FeI`      | F,F   | Dielectronic flags for Fe I. |
| 254  | `FIX_INB_NE`   | F     | Do not fix the electron density at the inner boundary. |
| 255  | `REV_RGRID`    | F     | Do not automatically revise the radius grid. |
| 256  | `OB_METH`      | HONJ  | Outer boundary method. HONJ = the "HonJ" (H-of-J) closure, which uses the Eddington factor H/J. |
| 257  | `DO_LAM_AUTO`  | T     | **Automatically switch between Lambda and Full BA iterations.** When MAXCH < 50%, CMFGEN auto-switches `DO_LAM_IT` from T to F, beginning Full BA iterations. `LAM_VAL=400` / `NUM_LAM=2` provide fallback: if MAXCH > 400, stay Lambda; every 2 Lambda iterations, try one Full BA. This is the standard CMFGEN approach — no IT_SPECIFIER is needed. |
| 258  | `OB_RAT`       | 1.0D0 | Ratio used by the outer boundary condition (1.0 = no extra scaling). |
| 259  | `STOP_IF_BP`   | F     | Do not abort the run if a bad/inconsistent parameter is detected — print a warning and continue. |

---

## Quick Reference: What a Typical SN Run Needs

If you are setting up a new SN snapshot, the parameters you will most
often need to touch are:

- **`SN_AGE`** — the epoch of the snapshot, in days.
- **`RSTAR`, `RMAX`, `VINF`, `VCORE`, `BETA1`** — only if you are not
  using `SN_HYDRO_DATA` for the structure.
- **`MASS`** — the ejecta mass.
- **Composition (`HYD/X`, `HE/X`, ...)** — usually overridden by
  `SN_HYDRO_DATA`, but check.
- **`FIX_T`** — usually T for a single snapshot (take T from the hydro).
- **`LTE_EST`** — T for a cold start (all configs, 2026-08-03). The
  `proxy_populations` alternative (LTE_EST=F + custom _IN files from
  `proxy_populations.py`) is deprecated — empirically worse than
  `LTE_EST=T` for cold snapshots because `FIX_T=T` makes CMFGEN ignore
  the _IN file's T column. See AGENTS.md "proxy_populations vs CMFGEN built-in
  LTE_EST".
- **`DO_NG`, `IBEG_NG`, `ITS/NG`** — Ng acceleration controls. Manual:
  "important tool for obtaining convergence" with NUM_BNDS=1/DIAG.
  Enabled (`DO_NG=T`) with `IBEG_NG=15` (start after 15 iterations),
  `BEG_NG=10%`, `ITS/NG=8`. Previous `DO_NG=F` caused infinite Lambda
  loop (549 iterations, no convergence).
- **`MAX_LIN`, `MAX_LAM`** — convergence aggressiveness. Manual:
  `MAX_LAM` "typically 10^5"; we use 1.0D2 (100x) for LTE cold start
  safety. Previous `MAX_LAM=1.1D0` (10% cap) caused infinite Lambda
  loop — trace species needing 10000% correction could only change
  10% per iteration.
- **`EPS_TERM`** — when to stop iterating.
- **`F_LAM_BEG`, `F_LAM_END`** — the wavelength range of the output
  spectrum.

### Non-monotonic velocity (interacting SNe)

For models with non-monotonic velocity (shell collisions, velocity
reversals, sigma < -1), the following MUST be set consistently in
`vadat_overrides`:

```yaml
vadat_overrides:
  USE_J_REL: "F"       # Disable MOM_JREL_V9 (caps sigma at -1.0 → NaN)
  USE_LAM_ES: "T"      # J=TC formal solution (bypasses moment equations)
  INCL_REL: "F"        # Required: CMFGEN checks USE_J_REL=T if INCL_REL=T
  USE_FRM_REL: "F"     # Consistent with non-relativistic mode
  N_TYPE: "G_ONLY"     # N_ON_J is only for MOM_JREL_V9 (needs USE_J_REL=T)
```

`prepare_VADAT.py` has automatic consistency enforcement: when USE_J_REL
is overridden to F, it auto-disables INCL_REL and USE_FRM_REL even if not
explicitly listed. See AGENTS.md for the full debugging history.

Everything else can usually be left at the defaults from
`prepare_VADAT.py`.

### BA iteration control (DO_LAM_AUTO-only, revised 2026-09-02)

For non-monotonic models, the iteration mode is controlled by CMFGEN's
built-in `DO_LAM_AUTO` mechanism — **no IT_SPECIFIER is used**.

The previous approach used an IT_SPECIFIER schedule with hardcoded
LAMBDA(FIXT) and FULL(FIXT) blocks, but this forced BA iterations too
early (at iteration 12 when MAXCH was still >500%), producing physically
meaningless corrections.

The current approach lets CMFGEN decide when to switch from Lambda to
Full BA:

1. **Lambda phase**: `DO_LAM_IT=T` → all iterations are Lambda until MAXCH < 50%
2. **Auto-switch**: `DO_LAM_AUTO=T` switches `DO_LAM_IT=F` in IN_ITS
3. **Full BA phase**: Full BA with `LAM_VAL=400` / `NUM_LAM=2` fallback
4. **Automatic fallback**: if BA overshoots, `solve_for_pops.f` returns to Lambda

To use this approach, simply omit `iteration_strategy` from the YAML
config. `prepare_VADAT.py` will delete any existing IT_SPECIFIER so
that `SPECIFY_IT_CYCLE_V3` returns immediately.

Verification from OUTGEN:
```
DIAG: iter=   1 IT_SPEC_PRESENT= F
DIAG: LAM_BEFORE= T RD_LAMBDA= T
DIAG: LAM_AFTER= T COMP_BA= T FIX_T= T
DIAG: OVERRIDE branch (RD_LAMBDA=T)
 LAMBDA iteration used
```

During Full BA iterations, the zone-dependent BA is active:
```
DO_VAR_CONT: zone-dependent BA, 23 compression depths of 160 held fixed-J
```

---

## Further Reading

- `AGENTS.md` — project-specific conventions and the history of
  convergence fixes for the X_2 model.
- `updates2cmfgen.md` — detailed log of every source-code change made to
  CMFGEN for this project, with rationale.
- `scripts/prepare_VADAT.py` — the script that generates this VADAT file
  from a YAML config; reading it shows which parameters are tuned per
  model and which are fixed defaults.
- Hillier, D. J. (CMFGEN documentation) — the underlying physics and
  numerical methods.
- Mihalas, *Stellar Atmospheres* — the standard textbook for non-LTE
  radiative transfer.
- Hubeny & Mihalas, *Theory of Stellar Atmospheres* — modern treatment
  of the same material.
