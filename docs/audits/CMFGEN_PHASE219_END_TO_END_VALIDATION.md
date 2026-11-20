# Phase 2.19 - end-to-end KNOP population and observer validation

Provenance: HEAD 48cd5d49 (includes Phase 2.18 commits 2b0254a3,
66ea1e65, db184329 plus user fixture commits c4cff903, 1a7400fc,
48cd5d49). Working tree preserves the user's uncommitted doc edit
(CMFGEN_KNOP_PRODUCTION_CONVERGENCE_FIX.md section 7) and run
products.

Executables this phase:
- cmf_flux_p219  e3db3555a60a0200ed9b5731531400b4403367e853c2926da9c9fdfbe61c176e
- knray_p218     ee03b158402f678af39fe246dd4d5477a8195519c7c5e1d934b5fa308866515f
- lib/libsubs.a  repaired 631-member archive (Phase 2.18)

## PART A - final acceptance ordering (source-level proof)

Per great iteration inside CMFGEN_SUB:

```
sweep (10000 loop)     intensities/moments/RJ/rates accumulated on
                       the lagged store as it stands during the march
cmfgen_sub.f:3316      CALL FG_KNOP_REFRESH_V1 - re-solves rays on
                       captured opacities, writes I_P/I_M_DEPTH_NEXT,
                       census sets FG_RF_ALL_CONV (fg ~9305)
cmfgen_sub.f:4603      CALL SOLVE_FOR_POPS - population update, MAXCH,
                       controller bookkeeping
solve_for_pops:1106    GET_FG_GS_MATURITY export (fresh: refresh ran
                       ~1300 lines earlier in the same call)
solve_for_pops:3507    LAM_FIELD_NOT_READY=LAM_FIELD_GATE(...)
solve_for_pops:3525    fold-in: rf_not_allconv forces it T when
                       .NOT. GS_RF_ALLCONV
cmfgen_sub.f:5513      CONV_EXIT_MET requires .NOT. LAM_FIELD_NOT_READY
                       (identical conjunct at 5081)
cmfgen_sub.f:5119      converged exit only under LST_ITERATION
```

LAM_FIELD_GATE (lam_handoff_mod.f90:37) =
`HELD .OR. FREEZE .OR. (HAS_LAGGED .AND. (RAW_IMMATURE .OR.
N_MATURE .LT. REL_STABLE))`. GS_RAW_IMMATURE (solve_for_pops:1596)
includes GS_FG_IMMATURE, priced from GS_MAX_ALL which folds the
refresh retained delta GS_RF_DCHG.

Answers to the seven sub-questions:

1. Can MAXCH use rates built from the preceding field? YES - by
   construction. A call's rates come from that call's sweep, whose
   lagged store was seeded from the PREVIOUS call's refresh output.
   The same-call refresh commits only after the sweep ends.
2. Can CONV_EXIT_MET accept before the next iteration rebuilds
   rates on the verified field? YES on the first verified call:
   flag=T at call N + MAXCH<EPS exits, while call-N rates were
   built on the refresh-(N-1) seed. Bounded, not open-ended: if the
   verified refresh moved the field materially, GS_RF_DCHG keeps
   the field raw-immature (LAM_FIELD_NOT_READY via GS_RAW_IMMATURE)
   so the exit is held anyway; acceptance on a verified call
   therefore implies the seed was within the 5% maturity band of
   the verified fixed point.
3. Does the held exit always force the subsequent calculation?
   YES - no path marks the model converged while
   LAM_FIELD_NOT_READY; non-convergent terminations (iteration
   limits, watchdogs) report failure, never acceptance.
4. Early exit bypass? None found: both CONV_EXIT_MET sites carry
   the conjunct; LST_ITERATION wraps the converged branch.
5. Provisional refresh overwriting the verified field before
   revalidation? No window: refresh runs once per call before
   SOLVE_FOR_POPS, flag recomputed each census; acceptance only
   reads the same call's verdict.
6. Previous population residual reused? No - MAXCH/MAXCH_APP are
   recomputed inside each call's population update.
7. Accepted field == evaluated field? Within the maturity band:
   the evaluated field is the sweep's field (seeded by the previous
   refresh); the accepted store is this call's verified output.
   They coincide exactly only at a reached fixed point; otherwise
   bounded by GS_RF_DCHG <= GS_MATURE_TOL (5%) whenever the gate
   let the call through.

Verdict: the one-call lag is a bounded gap, not a bypass - the
maturity conjunct already requires the verified refresh's retained
delta to be small before acceptance can fire. Strictest possible
hardening (require FG_RF_ALL_CONV on two consecutive calls) is a
policy option, recorded here, not implemented: it adds one held
iteration to every converged exit and no defect was demonstrated.

## PART B - population-level tests

A driver for SOLVE_FOR_POPS does not exist in unit_tests; building
one requires the full model environment (populations, atomic data,
rates). Outstanding. The refresh-level semantics they need are
covered: cerb (bulk=3 -> flag F, provisional commits capped),
cern (full=3 -> T), rfmx (itlim -> F), rfab (T->F stale-state),
cerp (probe excluded + bulk -> F).

## PART C - observer route validation

- Caller inventory: GET_FG_GS_MATURITY is called at exactly three
  sites: solve_for_pops.f:1106 (updated), obs/cmf_flux_sub_v5.f:2492
  (fixed in 66ea1e65 - now 18 args incl. GS_ALLCONV_OBS, printed as
  `all-ray refresh conv`), unit_tests/f90/drv_knop_ray.f90 (updated).
  No other call sites exist in the tree.
- cmf_flux_p219 built from the repaired production libraries with
  the production obs link line (obs/ sibling .o + -lsubs -lpgplt
  -ldisp -lmod, ifort -qopenmp -mkl): clean link, no undefined
  symbols - the corrected interface resolves.
- OBS execution: slurm job 13420402 (interactive) runs batobs.sh
  with CMF_OBS_EXE=cmf_flux_p219 on the m15_2p82d_FRESH_knop model
  (identical inputs/settings to the Nov-11 production OBS run).
  First attempt failed in setup (unset $ATOMIC env - script
  requirement, corrected). Result recorded in
  docs/audits/knop_phase219/ when the job completes; pending at
  audit write time.

## PART D - T21b/T26/T27 root cause and classification

Source-level cause: FG_RAY_SEEDED (allocated per refresh, fg:3834)
marks rays that have taken a WRITE_STORE commit. Under
FG_SEED_COMMIT=T (default), an UNSEEDED ray skips the
FG_COMMIT_JUMP veto for its first commit (fg:10620-10624) - the
in-sweep pre-fill defeats the cold-start exemption in FG_JUMP_MAX,
so the veto would otherwise reject exactly the converged fill it
exists to protect.

The failing run caught the harness in an intermediate state:
FG_RAY_SEEDED was allocated by the driver but left at 0, so
L_UNSEED armed on every ray and the jump veto never evaluated
(JREJ=0; stores committed normally). Commit c4cff903 completed the
fixture: SEEDED=1 for legacy cases (reproduces the pre-seeded
state where the veto does evaluate), SEEDED=0 for the cert cases
(to expose the waiver and the SEED marks), and 48cd5d49 pinned the
expectations. Current suite: 49/49 test_knop_ray.py, 9/9
test_knop_mode1.py.

| Test  | Original expectation        | Actual (broken fixture)    | Root cause                              | Independent evidence                  | Classification           |
|-------|-----------------------------|----------------------------|-----------------------------------------|---------------------------------------|--------------------------|
| T21b  | jump veto fires (JREJ=3)    | JREJ=0, store committed    | SEEDED=0 armed the first-commit waiver  | fg:10620 L_UNSEED; c4cff903 fixture   | incorrect test setup     |
| T26   | pending delta carried       | committed, no pending      | same: veto never evaluated              | same                                  | incorrect test setup     |
| T27   | persists across refreshes   | same                       | same                                    | same                                  | incorrect test setup     |

The waiver is production semantics (02-Nov-2026 feature): seeded
rays evaluate the veto, unseeded rays get one unguarded commit
gated by FG_COMMIT_CONV on the certification class. Verified in
source, and by the suite now passing with explicit SEEDED control.

## PART E - archive and executable reproducibility

- lib/libsubs.a: 631 members, all ELF; 15 previously-corrupt
  members rebuilt from source (tools/* + var_op_v6); fresh ranlib
  index; 0 stale members vs source; no debug-built objects;
  backups kept outside lib/.
- cmf_flux_p219: built against the repaired archive + obs sibling
  .o (cmf_flux_sub_v5.o recompiled with the corrected call).
- exe/cmf_flux.exe untouched (frozen); knray_p210-p218 untouched.

## PART F - live model A/B

The m15_2p82d_FRESH_knop production run's ten-release non-
convergence predates the all-ray flag; a controlled A/B requires
full model re-runs (slurm, hours-days). Outstanding. The OBS job
above is the first live-model executable check of the corrected
interface.

## PART G - mixed-sign limitation

Still no genuine velocity-inversion model: all available decks are
all-compressing m15 hydro. Requirements unchanged: real V(r)
non-monotonicity placing seeded cells on both GAM signs, both
regions active in the seeded store, then the full convergence and
OBS validation sequence.

## FINAL REPORT

1. Acceptance before rates rebuilt on the verified field: possible
   on the first verified call, bounded to the 5% maturity band by
   GS_RF_DCHG in the held-exit gate. Not a bypass.
2. One-call lag handled: yes, bounded by the maturity conjunct;
   strict two-call hardening recorded as an option.
3. Moments/RJ rebuilt after final verified refresh: yes, next
   great iteration's sweep - if the model did not already exit.
4. Rates rebuilt before final MAXCH: same-call ordering - rates
   precede MAXCH every call; the lag is relative to the refresh,
   not the population solve.
5. Stale residual acceptance: no.
6. Provisional overwrite before acceptance: no window (same-call
   flag).
7. cmf_flux with corrected interface: links and runs (job
   13420402); OBS execution result pending at write time.
8. OBSFLUX diagnostic on/off: pending OBS job completion.
9. T21b/T26/T27: incorrect test setup (unarmed waiver),
   resolved by c4cff903+48cd5d49; suite 49/49 + 9/9.
10. Archive reproducible: yes - 631 valid ELF members, clean
    index, documented provenance for the 15 rebuilds.
11. Ten-release live comparison: not run - outstanding.
12. Mixed-sign: still unvalidated; no inversion model exists.
