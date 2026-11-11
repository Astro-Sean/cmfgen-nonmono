# KNOP Outer Population-Opacity Drift Diagnosis (P1-T1c, PARTS 5-6)

Status: ANALYSIS + EXPERIMENT DESIGN. The frozen-opacity and
population-feedback experiments (PARTS 5-6) require production
`run_model_grid.py`/Slurm executions that were **not run in this
session**; they are specified precisely below so they can be executed
and the open questions settled. This report does not claim results it
did not measure.

Companion report: `KNOP_PRODUCTION_DEPENDENCY_GRAPH.md` (PARTS 1-2).

Evidence classes used: **[MATH]** code-structure result, **[MANUF]**
manufactured-model measurement, **[PROD]** production-model measurement,
**[HYP]** unverified hypothesis.

---

## 1. Why the drift points to the outer iteration

The radiation-store sub-map is now well characterized:

- The compressing lagged read is strictly slot-ascending
  (`read = write + NM_RAY`, `NM_RAY >= 1` always) - **[MATH]**.
- On the manufactured fixed-opacity model the store dependency graph is
  a verified DAG (Tarjan SCC, all thresholds) and the sweep Jacobian is
  nilpotent; the full sweep+refresh map has `J = 0` - **[MANUF]**.
- The refresh solves each ray's field to `KN_RAY_TOL` and commits
  unconditionally on clean convergence, so the store is slaved to the
  (fixed) opacities - **[MATH]**.

If the fixed-opacity radiation map is contractive/seed-independent, a
persistent `GS_MAXCHG` / field-maturity stall cannot be generated inside
the radiation solve. The remaining degree of freedom that changes every
great iteration is the **opacity/emissivity pair `(chi, eta)`**, which is
rebuilt from the populations each great iteration. The radiation field
is therefore slaved to `(chi,eta)`, and `(chi,eta)` is slaved to the
population vector `n`. Persistent field drift must enter through
`n -> (chi,eta) -> J -> R_ij -> n`. **[MATH - inference]**

This is consistent with, but does not by itself prove, the standing
verdict that the stall lives in the outer population-opacity iteration.
The two candidate mechanisms already on record are:

1. **Population floor-pinned approximate null-space** (Pillar 5,
   `KNOP_SECOND_PASS_AUDIT.md`): `C_MAT(I,J) *= POPS(J,L)` in
   `generate_full_matrix_v3.f` makes floor-pinned levels
   (`POPS ~ 1E-30`) produce approximate-null Jacobian columns;
   `DGEEQU`/`DGETRF` unscale them into `1E21-1E53` fractional-correction
   demands. Confirmed **[MATH]** on a fast model and corroborated on
   `VAL_knop_formal` (luminosity-divergent runaway, it1489).
2. **Field-maturity stall** (`FRESH_knop`, `GS_REL_NOPROG_LIM`): a
   distinct GS freeze/release limit cycle with moderate (~38%) demand -
   **[MATH]** confirmed as a *different* mode from the runaway.

Neither mechanism has been shown to be the *proximate* cause of the
persistent radiation-field drift specifically; that is what PARTS 5-6
measure.

---

## 2. PART 5 - frozen-opacity and coupled-iteration experiments

All three arms start from the **same checkpoint** at a selected
problematic production iteration (e.g. a `FRESH_knop`-class stall or a
`VAL_knop_formal`-class pre-runaway iteration).

### Checkpoint contents (to save)

- `POPS` (all depths x levels), `n`-relevant population residuals.
- `CHI_EXT`, `ETA_EXT` (and `KN_CHI_EXT`, `KN_ETA_EXT`, `KN_FREQ`,
  `KN_DLN_P/N`, `KN_DBB`, `KN_IC` - the sweep-captured opacity set).
- `J` (mean intensity) and the `I_*_DEPTH_NEXT` store.
- Population-solver residuals: `maxraw`, `exm`, `maxapp`, `rcond`,
  `colmin_u`, `nflr` where available (`POP_AUDIT=T`).
- `GS_MAXCHG`, `GS_MAX_ALL`, per-iteration `R_ij` for the dominant
  transitions.

### Arm A - frozen radiation test

Hold `POPS`, `CHI_EXT`, `ETA_EXT` (and boundary conditions) fixed.
Repeat the complete production radiation solve (sweep + refresh) for
several great iterations **without** a population update.

Measure per iteration: store residual `max|dI_NEXT|`, radiation-field
residual `max|dJ|`, formal-solve residual, iteration count, and whether
the solution becomes stationary.

**Prediction if the radiation map is benign:** the field converges to a
fixed point and `max|dJ| -> 0`; any residual drift is bounded by the
refresh tolerance `KN_RAY_TOL`. A *persistent* nonzero `dJ` under frozen
opacities would refute the "radiation map is contractive" conclusion and
re-open the field-side hypothesis. **[HYP - to be measured]**

### Arm B - one-step population perturbation

From the same frozen baseline:

1. Perform **one** population update (single `SOLVEBA`/JFNK step).
2. Record `delta_n_i` per level.
3. Recompute `chi`, `eta`; record `delta_chi_nu`, `delta_eta_nu`.
4. Re-solve radiation; record `delta_J_nu` and `delta_R_ij`.

This isolates the single-step gain `G = dJ/dn` and identifies which
species/levels produce the largest opacity changes.

### Arm C - normal fully coupled iteration

Run the production iteration from the same checkpoint, recording the
same quantities at every great iteration. Compare drift rates and
residual trajectories against Arms A and B.

**Discriminating power:**
- If Arm A is stationary but Arm C drifts -> the drift is injected by
  the population-opacity update (the `n -> (chi,eta)` map), not the
  radiation solve.
- If Arm A itself drifts -> the radiation map is not contractive on this
  geometry (back to PART 3).
- If Arm B shows a single population step produces `delta_J` comparable
  to the per-iteration drift in Arm C -> the population->opacity->field
  gain is O(1) and the outer iteration is the amplifier.

---

## 3. PART 6 - population-opacity feedback diagnostics

For representative problematic species/levels (the `POP_AUDIT` census
candidates: floor-pinned levels like `FeI:3d6_4s2_a5De[2]` @d63 with
`POPS~1E-31`, plus the levels that dominate `BA_DEM_EXR`):

Per great iteration record `delta_n_i`, `delta_chi_nu`, `delta_eta_nu`,
`delta_J_nu`, `delta_R_ij` (absolute and relative). For each species/level
also record its **fractional contribution** to total `chi_nu` and
`eta_nu` so that a large *relative* change in a floor-pinned trace
population is not mistaken for a significant opacity change.

Specific questions to answer:

- Which species dominate `delta_chi`/`delta_eta`?
- Which transitions dominate `delta_R_ij`?
- Does the population residual (`maxraw`/`exm`) decrease while `J` still
  drifts? (decoupled convergence -> field slaved to a moving opacity)
- Does `delta_chi`/`delta_eta` *precede* or *follow* `delta_J`?
  (opacity-driven vs field-driven feedback direction)
- Do the same levels repeatedly drive population-release or demand spikes
  (`BA_DEM_EXR`, `RUNAWAY` channels)?

### Conditioning check (no matrix change)

Recompute per-depth `rcond`, `colmin_u`, `nflr`, and the
`C_MAT *= POPS` column-norm distribution at the checkpoint depths. Test
whether the approximate-null columns identified in Pillar 5 are present
at the *same* depths that dominate `delta_chi`. A correlation between
near-null population columns and opacity-drift hotspots would support
the conditioning mechanism as the drift driver; absence would point to a
genuine slow outer-mode instead. **[HYP]**

---

## 4. Honest status and what would settle it

| Question | Status |
|---|---|
| Radiation store map nilpotent for production geometry | OPEN - needs PART 3 spot-check |
| Refresh `J=0` genuine vs tolerance artefact | OPEN - needs PART 4 multi-seed test |
| Frozen-opacity radiation converges independently | UNTESTED - Arm A |
| Population step -> `delta_chi/eta/J` gain magnitude | UNTESTED - Arm B |
| Drift lives in outer `n->(chi,eta)` map | HYP - strongly suggested, unproven |
| Population null-space causes the field drift | HYP - mechanism confirmed, causal link unproven |

The single most informative next run is **Arm A** (frozen opacity): it
decisively splits "radiation-side" from "population-opacity-side" drift
with one production execution. Arms B and C then quantify the gain and
localize the driver. All three need a `FRESH_knop`-class or
`VAL_knop_formal`-class checkpoint with `POP_AUDIT=T` and OUTGEN
preserved (the `RUNAWAY`/`POP_AUDIT` channels do not survive the
restart-prep consumption noted in `KNOP_SECOND_PASS_AUDIT.md`).

Until those run, the correct statement is: **the radiation-store map is
contractive on the tested geometry and cannot by itself sustain drift;
the persistent field drift is consistent with an outer
population-opacity feedback, most plausibly the floor-pinned
approximate-null-space mechanism, but the causal link and the exact
driver levels remain unverified.**
