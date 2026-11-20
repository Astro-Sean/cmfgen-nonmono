# CMFGEN KNOP bulk-tolerance live-model validation (Phase 2.15)

Status: provenance, store-commit lifecycle, and production-run
evidence complete. The live-model A/B sweep (population and
spectral differences) could not be run in this session: it needs
full CMFGEN model runs with different VADAT values (slurm jobs,
hours-days). What exists on disk is analysed here; what requires
new runs is listed explicitly in section 7. This is an honest
partial report, not a claim of live validation.

Provenance: HEAD 2f37d53a; model run
`cmfgen_validation/m15_2p82d_FRESH_knop/cmfgen/2.82days`
(VADAT `T [DO_NONMONO_KNOP]`); SN_HYDRO_DATA 70 depths.

## 1. Live-model provenance and sign coverage

SN_HYDRO_DATA: V 48.9 -> 7073.6 km/s, monotonically increasing;
SIGMA (dlnV/dlnr-1) -1.037 to -156.9 -- **all 70 depths
compressing, zero expanding cells**. Same for the two other
available models checked (handoff_smoke/knop_mom, tmp/v13_audit -
identical hydro).

Therefore:

1. A genuine live velocity-inversion model with expanding cells
   is NOT available in the current run set.
2. The Phase 2.14 mixed-sign question (both signs in the seeded
   state) remains untested on a real model - the only mixed data
   is the synthetic `p1t1s_mixed.txt` deck.
3. The m15 deck itself IS the production reference: its seeded
   state is all-compressing, matching the nilpotent-omega=1
   structure proven in Phase 2.11.

## 2. Store-commit lifecycle audit (WRITE_STORE / WRITE_STORE_CAP)

Two distinct commit paths in `KNOP_RAY_SOLVE_V1` (~10320-10415):

- Converged path (`WRITE_STORE`): `STORE = T2 + INEXT_RELAX*(T1 -
  T2)` - the new field blended under the refresh relaxation
  (INEXT_RELAX=0.5 production).
- Provisional path (`WRITE_STORE_CAP`, FG_COMMIT_CONV fail-close):
  `V1 = V2 + SIGN(min(|V1-V2|, FG_COMMIT_NC*max(V2, 1e-4*VMAX_F)),
  V1-V2)` - a per-cell capped ratchet, no blending; floor/NaN
  guard; NSTORE counts large changes.

Downstream lifecycle: the store is read by SEED_FLD at the next
refresh pass and by the population iteration between refreshes
(via I_*_DEPTH_NEXT). A provisional store therefore DOES seed the
next refresh and DOES feed populations - it is not inert
scratch. The distinction between provisional and converged is
carried by the exit class (IXIT/NXIT), not by a store-level
flag: once written, the store itself no longer knows it was
capped.

Answer to the key question: FG_COMMIT_NC is bounded (<=5%
relative per cell per refresh by default) but NOT inert - it
writes real radiation field into the lagged store which
populations then see. Safe as a progress mechanism; NOT safe to
treat as a converged solution.

## 3. Production-run evidence (m15_2p82d_FRESH_knop)

From OUTGEN (run with old `<=0` semantics and VADAT default
6e-2):

- Per refresh: `bulk= 7-18` rays exited via the bad-fraction
  test (~10-20% of 84 rays), `bulk rays: maxres= 1.8e-2 - 3.1e-2`
  (9-15x KN_RAY_TOL=2e-3).
- Unconverged cell fraction per refresh: 1.1e-5 - 2.2e-4.
- The run ended `SOLVE_FOR_POPS: 10 releases ... stopping as
  non-convergent` - the model never converged under the
  production configuration.

These bulk exits committed fields ~10-30x above tolerance into
the store every refresh. Under the corrected Phase 2.14
semantics, `BTOL=0` disables this exit entirely (the documented
zero-semantics); `BTOL=6e-2` retains it as a *named* loose-floor
acceptance.

## 4. What the evidence does and does not show

Established:
- The certification defect (Phase 2.13/2.14) is real and
  frequent in production: ~15% of rays/refresh committed at
  ~10-30x tolerance.
- The provisional store writes real field data that populations
  consume - capped, not converged.
- With `BTOL=0` the same solve reaches ixit=1 at RES=1.9e-4,
  qres=0 in identical pass count on the extracted deck - there
  is no cost argument for the loose floor on this operator.

Not established (requires new model runs):
- Population deltas n_bulk vs n_reference at BTOL=0/1e-3/1e-2/
  3e-2/6e-2 - needs 5 full model runs.
- Spectral deltas on identical observer settings - needs the
  same runs plus obs synthesis.
- Reproducibility/warm-start equivalence on the live model.
- Whether BTOL=6e-2 changes the *final* converged population
  fixed point or only the path length to it.

## 5. Recommendation (bounded by available evidence)

- Keep `KN_RAY_BULK_TOL=0` available and now correctly
  disabling - it is the scientifically clean default and costs
  nothing on the tested operator.
- The production VADAT default of 6e-2 under DO_NONMONO_KNOP
  (set in `rd_control_variables.f:1468`) should be justified by
  a live A/B run before continued use: it accepts fields up to
  30x tolerance as bulk-converged. If retained, treat those
  commits as provisional (FG_COMMIT_NC ratchet) rather than
  converged.
- Do NOT allow a provisional (capped) store to mark a ray
  converged for final population convergence; allow it only as
  a progress mechanism (which is exactly what FG_COMMIT_NC is).

## 6. Answers

1. Genuine live velocity-inversion model: NO - all available
   models are all-compressing m15 hydro.
2. Both signs in seeded state: NO on any available model.
3. FG_COMMIT_NC safe through lifecycle: bounded but not inert -
   writes real capped field that seeds next refresh and feeds
   populations; the store carries no converged/provisional flag
   after the write.
4. Does BTOL=6e-2 produce population differences: not measurable
   here - needs the A/B runs (production evidence shows 15% of
   rays/refresh committed at ~10-30x tolerance).
5. Most-sensitive species/levels: not measurable - needs runs.
6. Spectrum changed: not measurable - needs runs.
7. Reproducible across starts: not measurable on live model.
8. Unconverged cells concentrated: yes - the bulk exits are
   driven by a tiny bad fraction (1e-5 - 2e-4 of cells); the
   41/6227 bad-cell count on the diagnostic deck localises to
   the compressing chain.
9. Defensible production value: BTOL=0 pending live A/B; the
   only argument for >0 is the parked-store tradeoff, which
   FG_COMMIT_NC already handles separately.
10. Provisional stores in final convergence: NO - progress
    mechanism only; the corrected gate already enforces this.

## 7. Required follow-up runs (exact procedure)

For a real validation: run `m15_2p82d_FRESH_knop` (or the
configured KNOP recipe) with VADAT `KN_RAY_BULK_TOL` set to
{0, 1e-3, 1e-2, 3e-2, 6e-2} from identical cold populations,
then obs synthesis on identical settings, and diff
RVTJ/BASIC_MODEL/observed spectra pairwise vs the BTOL=0 case.
The Phase 2.14 `p1t1x` cert mode gives the per-ray exit/status
record for the same data in the diagnostic harness.

## 8. Safety / regression

- No source changes this phase (analysis + audit only);
  knray_p214 remains the frozen Phase 2.14 executable.
- Production defaults, controller, equations, population and
  observer solvers untouched.
