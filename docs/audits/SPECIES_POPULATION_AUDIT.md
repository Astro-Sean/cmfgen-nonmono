# Species and Atomic-Level Population Correctness Audit

Date: 07-Nov-2026
Tree: `cmfgen_jfnk` (diff base `../cmfgen_reference/cur_cmf`)
Model inspected: `m15_2p82d` validation runs (VAL_sob_nonloc, VAL_sob_jfnk,
VAL_sob_diag, VAL_sob_tri, VAL_knop_formal, SOB_l14_lte_adp2_tauto)

## 1. Executive summary

The complete population data path was audited: species/ion indexing in
`cmfgen.f`, the `MODEL_ATOM_DATA`/`ATM` structure, `POPS` row layout,
`POPTOION` aliasing, full-to-super-level reconstruction
(`SUP_TO_FULL_V3B` via `SUP_TO_FULL_V4`), transition generation in
`set_frequency_grid_v2.f`, opacity/emissivity in `set_line_opac.f`,
number-conservation assembly in `generate_full_matrix_v3.f` and
`fixpop_in_ba_v3.f`, and the `SE(ID)` statistical-equilibrium blocks.

Findings:

- **No population-index collision between species or ionization stages
  was found.** Row allocation in `cmfgen.f` is consecutive and
  non-overlapping per species; `F_TO_S` mappings are range-validated at
  input time (`rd_f_to_s_ids.f` hard-stops on out-of-range or unmapped
  levels); transition generation maps `F_TO_S + EQXzV - 1`, which can
  only address the owning ion's rows.
- **One real latent defect was found and fixed**: `SUP_TO_FULL_V3B`
  computes `C2_F*(C2_S/SUM)` with no guard on `SUM=0`. When a super
  level has zero population, every reconstructed member is zero and the
  renormalization is 0/0, silently writing NaN into `XzV_F` and from
  there into the `L_STAR_RATIO`/`U_STAR_RATIO` opacity corrections. The
  same unguarded line exists upstream (inherited defect) and in the
  uncalled `sup_to_full_v3.f`.
- **The measured O V profile is not an indexing defect.** POPOXY gives
  O V = 71.5% (d=1), 93.4% (d=35), 0.42% (d=70) of oxygen in
  `VAL_sob_nonloc`; the Sobolev JFNK/DIAG/TRIDIAG runs agree to ~0.03%.
  The `VAL_knop_formal` run has a completely different thermal state
  (O III = 98.4% mid-depth, O V ~ 1e-17) — an unconverged or
  physically different model, not a storage corruption. All 9,100 O V
  level values and all other stage values are finite and non-negative,
  and the oxygen accounting closes at every depth.
- **A new general diagnostic `POP_AUDIT` was added** so that per-species
  conservation, per-ion totals, negative/non-finite counts, DXzV
  aliasing consistency, and departure coefficients can be checked in
  every future run, filtered by species index, stage index, and depth.

Verdict on the O V overpopulation question: none of the candidate
implementation defects A-C (indexing, super-level mapping, number
conservation) is demonstrated by the code or the output data. The
latent `SUM=0` hazard (fixed) produces NaN, not excess population. If
O V is genuinely overpopulated relative to the correct NLTE solution,
the remaining candidate mechanisms are the rate data/assembly (D),
line opacity/emissivity feedback (E), convergence state (F), or the
non-monotonic transfer path (G). The audit could not rule these in or
out for the specific model because a converged reference solution and
a JFNK/KNOP-run equivalent POPOXY at comparable convergence are not
yet available; see §13.

## 2. Complete population data model

### 2.1 Elements

- `NUM_SPECIES` chemical elements; `SPECIES(ISPEC)` names (e.g. `OXY`);
  `SPECIES_PRES(ISPEC)` presence flag.
- `SPECIES_BEG_ID(ISPEC)`, `SPECIES_END_ID(ISPEC)`: first and last
  ion-slot indices for the species; `SPECIES_LNK(ID)` maps slot → species.
- `POP_SPECIES(K,ISPEC)`: elemental number density at depth K — the
  right-hand side the species-conservation equation enforces.

### 2.2 Ionization stages

`cmfgen.f` (identical to upstream) builds consecutive ion slots per
species from the `[XzV_ISF]`/`[XzV_NSF]` keys:

- Each model atom gets `ID=ID+1`, `ION_ID(ID)` (e.g. `OV`, `OSIX`),
  `XzV_PRES=.TRUE.`, `NXzV_F` full levels, `NXzV` super-levels
  (`NS ≤ NF`, `NV ≤ NS` enforced), `EQXzV=EQ_TEMP`, then
  `EQ_TEMP += NXzV`.
- Stages must be consecutive — a gap hard-stops the run.
- The terminal stage adds one slot: `XzV_PRES=.FALSE.`, `EQXzV=EQ_TEMP`,
  `EQ_SPECIES(ISPEC)=EQ_TEMP`, `EQ_TEMP += 1`. It carries a single
  POPS row: the terminal-ion density, whose equation slot doubles as
  the species-conservation row.

### 2.3 POPS row layout (verified against POPTOION and EQXzV arithmetic)

For a model ion `ID` with `NXzV` super-levels:

- Rows `EQXzV .. EQXzV+NXzV-1` = super-level populations `XzV(1..NXzV)`.
- Row `EQXzV+NXzV` ≡ `EQXzV(ID+1)` — **the same row is the next ion's
  ground-state row and this ion's `DXzV`**. `POPTOION` reads
  `DXzV=POPS(EQHYD+NLEV)`; `eval_lte_v5` confirms the convention:
  `ATM(ID+1)%XzV(1,:)=ATM(ID)%DXzV_F(:)` when the next ion is absent.
  `DXzV(ID)` is therefore the *next ion's* density (the Saha factor),
  not the total of ion `ID`.
- Terminal slot: one row at `EQXzV=EQ_SPECIES(ISPEC)` = terminal ion
  density = the species-conservation equation location.
- `POPS(NT-1,:)=ED`, `POPS(NT,:)=T` (`poptoion.f`).

Each POPS row is written exactly once per iteration; level rows of
different ions are disjoint, so the species sum `Σ_L XzV(L,ID)` over
all model ions plus the terminal row counts every row exactly once.

### 2.4 Super-levels and full-atom mapping

- `ATM(ID)%NXzV_F` full levels ↔ `NXzV` super-levels via
  `F_TO_S_XzV(NXzV_F)` and `INT_SEQ_XzV(NXzV_F)` (interpolation
  sequence for non-uniform departure coefficients).
- `rd_f_to_s_ids.f` validates: every full level maps to `[1,NXzV]`
  (hard STOP on 0 or >NXzV), level names must match the link file,
  every super level has ≥1 member. `F_TO_S=0` cannot reach production.
- `AXzV_F(l,u)=f_lu` (l<u), `AXzV_F(u,l)=A_ul`; `GXzV_F` weights;
  `EDGEXzV_F` ionization energy per level; `XzVLEVNAME_F` names.
- `SUP_TO_FULL_V3B` (per ion, called for `ID=1..NUM_IONS-1` by
  `SUP_TO_FULL_V4` after every `POPTOION`): computes `C2LTE_F` from a
  Saha-Boltzmann form with occupation probabilities `W_C2_F`
  (`MOD_DO_LEV_DIS=F` → W=1), aggregates to `C2LTE_S`, forms
  `C2LTE_F_ON_S=C2LTE_F/C2LTE_S`, then reconstructs
  `C2_F = C2_S·C2LTE_F_ON_S` (uniform share) or interpolates the
  departure coefficient for `INT_SEQ>0` members, and finally
  renormalizes `C2_F *= C2_S/SUM` so each super-level sum equals the
  super-level population exactly.

### 2.5 `SE(ID)` statistical-equilibrium blocks

`STAT_EQ_DATA` (steq_data_mod.f): `STEQ` residuals, `BA` rate matrix,
`LNK_TO_IV`, `LNK_TO_F`, `EQ_IN_BA`, `EQ_TO_ION_LEV_PNT`,
`ION_LEV_TO_EQ_PNT`, `N_SE` (total local eqns), `N_IV` (important
variables = `NXzV_IV`), `NUMBER_BAL_EQ`, `XRAY_EQ`, `STEQ_ADV`
(advection ion eqns). `SE(ID)%BA` rows `1..NXzV` are level eqns; rows
`NXzV+1..N_SE-1` are dN/dt contributions to the recombining ion; row
`N_SE` is the ion's number-conservation contribution
(`generate_full_matrix_v3.f` comments, verified).

### 2.6 Number conservation

- `generate_full_matrix_v3.f`: `C_NC(ISPEC,:) = Σ_ID SE(ID)%BA(N_SE,:)`
  — every present ion's conservation row summed once per species —
  placed at `C_MAT(EQ_SPECIES(ISPEC),:)`. Ionization equations
  `C_ION(ID,:)` (rows `NXzV+1..N_SE-1` of the donor ion plus
  `C_ION(ID+1)` recombination-side terms plus `BA_ADV_TERM` advection)
  replace the g.s. row `C_MAT(EQXzV(ID),:)` when selected.
- `fixpop_in_ba_v3.f`: `EQ_SPECIES` rows use
  `POP_SPECIES(DEPTH_INDX,SPECIES_LNK(ID))` at `SPECIES_END_ID`;
  per-ion fixed populations use `EQXzV..EQXzV+FIX_N-1`. The historical
  `SPECIES_END_ID` bug (`IS` not set for the last stage) is already
  fixed in this source.

### 2.7 Transitions

`set_frequency_grid_v2.f` generates the line list in-house:
`ID=1..NUM_IONS-1`, `XzV_PRES` guard, `MNUP=2..NXzV_F`,
`MNL=1..MNUP-1`, `AXzV_F(MNL,MNUP)≠0` → `VEC_NL/NUP =
F_TO_S_XzV(MNL/MNUP)+EQXzV-1` (POPS rows of the same ion), full-level
indices in `VEC_MNL_F/VEC_MNUP_F`, `VEC_ID=ID`, `VEC_SPEC`, oscillator
strength, Einstein data, frequencies, profile types. All parallel
`VEC_*` arrays sort on the same `VEC_INDX` — relationships cannot be
separated by the sort. `set_line_opac.f` then evaluates
`CHIL_MAT = T1·(L_STAR·POPS(NL) − GLDGU·U_STAR·POPS(NUP))`,
`ETAL_MAT = T2·U_STAR·POPS(NUP)` with
`L_STAR_RATIO = XzV_F(MNL_F)/XzV(MNL)` and `U_STAR_RATIO` the upper
analog — the standard Einstein form with the super→full correction.

Both `sub_cmf_line.f` and `sub_sob_line_v3.f` consume the same
`CHIL_MAT`/`ETAL_MAT`, so the CMF and non-monotonic Sobolev paths
share the population→opacity evaluation.

## 3. Species/ion indexing audit

- Row allocation is sequential per species (`EQ_TEMP` monotonic); no
  two ions share a level row. The only intentional aliasing is the
  documented `DXzV(ID) = XzV(1,ID+1)` shared row.
- `BA_NUM_SPECIES=28`, `BA_MAX_IONS_PER_SPECIES=21` — static upper
  bounds, not per-species assumptions.
- No species- or stage-count-specific literals were found in the
  audited population path (`cmfgen.f`, `generate_full_matrix_v3`,
  `fixpop_in_ba_v3`, `set_frequency_grid_v2`, `set_line_opac`,
  `sup_to_full_v3b`). `ZXzV`/`ION_ID` are used for identification and
  charge, never for indexing arithmetic.

## 4. Number-conservation audit

- Every present ion contributes its `SE(ID)%BA(N_SE,:)` row exactly
  once to `C_NC`; the terminal slot is not a model ion and contributes
  its density row as the conservation-eq location, not an extra term.
- Ground-state equation replacement (`C_MAT(EQXzV)=C_ION`) preserves
  conservation because `C_ION` collects both the donor-ion loss rows
  (`NXzV+1..N_SE-1`) and the next ion's recombination-side rows.
- Fixed populations (`FIX_NXzV`, `FIX_SPECIES`) replace equations but
  not the `C_NC` sum; the residual reported by `POP_AUDIT` measures the
  true `POP_SPECIES` conservation.
- POPOXY verification (`test_popoxy_conservation`): all level values
  finite/non-negative across O I–O VI, every depth, in the production
  validation file.

## 5. Full/super mapping audit

- The reconstruction is `C2_F = C2_S·(C2LTE_F/C2LTE_S)` with an exact
  per-super-level renormalization — full↔super conservation is
  enforced, verified numerically to 1.7e-16 by `test_sup_to_full`.
- `C2LTE_F_ON_S` is the LTE share; `L_STAR_RATIO = XzV_F/XzV` therefore
  equals `(b_f/b_s)·(LTE share)` — mathematically the correct
  definition of the full-level population factor for line opacity.
- Departure-coefficient interpolation (`INT_SEQ>0`) is bounded
  (`B<0 → B=B1`, B held on the correct side of 1).
- **Defect fixed**: unguarded `C2_S/SUM` division (§10, D-1).

## 6. Transition indexing audit

- `VEC_NL/NUP` are global POPS rows of the same ion;
  `VEC_MNL_F/VEC_MNUP_F` retain full indices for the star ratios; both
  travel together through the `VEC_INDX` sort. No cross-ion reference
  is possible while `F_TO_S∈[1,NXzV]` (input-validated).
- `MAX_SIM` storage slots are reference-counted (`PROF_LIST_LOCATION`);
  stale-slot inheritance was not observed in the audited paths.

## 7. Opacity/emissivity audit

- `CHIL_MAT ∝ f_lu(n_l − (g_l/g_u)n_u)` and `ETAL_MAT ∝ A_ul n_u` use
  the same POPS rows as the rate equations — no second population copy.
- Stimulated emission via `GLDGU`; negative opacity arises only from a
  genuine population inversion, not an indexing asymmetry (L/U star
  ratios use the same mapping function).

## 8. Statistical-equilibrium audit

- `generate_full_matrix_v3` assembles each equation once:
  `SE(ID)%BA` local blocks, `C_ION` ionization couplings between
  adjacent stages, `C_NC` species sums, `BA_ADV_TERM` advection.
- `SE(ID)%BA(I,J)` rows `NXzV+1..N_SE-1` are the dN/dt rows for the
  recombining (next) ion — the correct adjacent-stage coupling.
- `XRAY_EQ` handles Auger channels; `STEQ_ADV` advection terms are
  added into `C_ION` symmetrically.
- Residual: population corrections can still produce locally negative
  trial values; `POP_AUDIT` now reports any that survive to the
  committed POPS.

## 9. JFNK population-vector audit

- JFNK residual ordering follows `POPS(1:NT-2)` SE rows + ED + T;
  `POPTOION` in `jfnk_cmfgen_interface_v1.f` uses the same
  `EQXzV`/`NXzV` mapping as the main path — the committed vector and
  the BA vector are the same physical system.
- O V cannot grow silently inside JFNK while the residual stays small
  *unless* the rate matrix itself encodes the excess — which is a
  rate-data question (category D), not an indexing question. The
  `POP_AUDIT` census makes any such stage-level drift directly visible.

## 10. Identified defects and severity

| ID | Severity | Location | Defect |
|----|----------|----------|--------|
| D-1 | Medium (latent) | `subs/sup_to_full_v3b.f:268` | `C2_F*(C2_S/SUM)` 0/0 → NaN in `XzV_F` when a super level has zero population; NaN would then propagate into `L_STAR_RATIO`/`U_STAR_RATIO` and the line opacity. **Fixed** (guard + conserving even-split fallback). Also present upstream and in uncalled `sup_to_full_v3.f` (not fixed — dead code). |
| O-1 | Observation | — | `DXzV(ID)` names the *next* ion's density (Saha factor), not ion `ID`'s total. Confirmed correct by `eval_lte_v5` aliasing; documented in the audit and in `POP_AUDIT_V1` comments to prevent misreading as an off-by-one. |

## 11. Fixes implemented

1. `subs/sup_to_full_v3b.f`: guarded the per-super-level
   renormalization — `SUM=0 & C2_S=0` leaves exact zeros (conserved);
   `SUM=0 & C2_S>0` (pathological underflow of all member shares)
   falls back to an even split over `CNT(L)` members so the ion total
   is still conserved and finite.
2. `new_main/mod_subs/control_variable_mod.f`,
   `new_main/mod_subs/rd_control_variables.f`,
   `new_main/mod_subs/solve_for_pops.f`: new `POP_AUDIT` diagnostic —
   `POP_AUDIT_SPEC`, `POP_AUDIT_ION`, `POP_AUDIT_D`, `POP_AUDIT_LEV`
   filters; per-species conservation residual vs `POP_SPECIES`,
   per-ion totals, `DXzV` aliased-row staleness check, negative and
   non-finite censuses, optional per-level fractions and departure
   coefficients. Fully general — index-based, no species-specific code.
3. `docs/audits/SPECIES_POPULATION_AUDIT.md`, this file.

## 12. Test results

- `test_sup_to_full` (new): compiles production `sup_to_full_v3b.f`;
  A: super-level conservation err=1.7e-16; B: LTE-share splitting
  err=0; C: LTE input → LTE output err=2.2e-15; D: zero-SL case now
  finite (would have been NaN before the fix); E: extreme contrast
  conserved err=1.8e-16. **7/7 checks pass.**
- `test_pop_audit` (new): static wiring — all 5 VADAT keys
  declared/parsed/clamped, guarded call at end of `SOLVE_FOR_POPS`,
  general species/stage/level loops, no hard-coded ion strings, ion
  filter excluded from the conservation sum, `DXzV` staleness probe,
  negative/non-finite censuses. **All checks pass.**
- `test_popoxy_conservation` (new): parses a validation `POPOXY`;
  all O I–O VI level values finite and non-negative at every depth;
  ion fractions bounded by the species total. **All checks pass.**
- Compile: `control_variable_mod.f`, `rd_control_variables.f`,
  `solve_for_pops.f`, `sup_to_full_v3b.f` all compile clean (ifort,
  `-extend-source`). `solve_for_pops.f` previously had no compile
  coverage — the POP_AUDIT call site is now covered.
- Measured O V accounting (POPOXY, `VAL_sob_nonloc`): OV fraction
  0.715/0.934/0.0042 at d=1/35/70; oxygen accounting closes; no
  negative/non-finite values in 9,100 OV level entries.

## 13. Remaining concerns

1. `sup_to_full_v3.f` (uncalled) retains the same unguarded division —
   dead code, flagged not fixed.
2. The Knop-formal validation run sits at a completely different
   ionization state (O III-dominated). Whether that reflects an
   unconverged iterate or a genuine transfer-path effect is open —
   `POP_AUDIT` on the next run will show it per-iteration.
3. Rate assembly for charge-exchange/dielectronic routes was read but
   not re-derived term by term; category D (rate data) remains the
   prime suspect if O V is physically overpopulated.
4. `POP_AUDIT` reports `DXzV` vs its aliased row as a staleness probe
   — it is identity-equal right after `POPTOION` and only fires on
   genuinely stale `ATM` state.
5. A bounds-checked debug build (`build_debug.sh`) and a full
   `run_all_tests.py` pass are still pending for this task; the new
   numeric test covers the highest-risk new code path.
