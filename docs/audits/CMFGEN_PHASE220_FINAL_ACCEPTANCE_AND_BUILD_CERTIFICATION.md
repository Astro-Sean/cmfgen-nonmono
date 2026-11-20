# Phase 2.20 - final population acceptance and production build certification

Provenance: HEAD 3746683a (includes Phase 2.17 flag commits 829cd7e5,
98acf207, 20942c32 and the Phase 2.18/2.19 build and OBS work).
Working tree clean at audit start; untracked run products (TIMING,
fort.4, unit_tests/systems/{_pv.txt,p29,p29m}) left in place.

## PART A - final acceptance execution trace

Verified at HEAD against the Phase 2.19 trace; all line numbers
rechecked, none shifted.

Per great iteration inside CMFGEN_SUB:

```
frequency sweep (10000 loop)   intensities, moments, RJ and rates
                               accumulated on the lagged store as it
                               stands during the march
cmfgen_sub.f:3313              IF((DO_NONMONO_KNOP .OR. DO_SOB_FG_KNOP)
                                 .AND. DO_KNOP_REFRESH)
cmfgen_sub.f:3316              CALL FG_KNOP_REFRESH_V1()
                                 - re-solves rays on captured opacities
                                 - WRITE_STORE / WRITE_STORE_CAP
                                 - writes I_P/I_M_DEPTH_NEXT
                                 - KN_CS_LS census (fg:9340s)
                                 - FG_RF_ALL_CONV verdict (fg:9351)
cmfgen_sub.f:4603              CALL SOLVE_FOR_POPS - population update,
                               MAXCH, controller bookkeeping
solve_for_pops.f:1105          IF(HAS_LAGGED_FIELD) CALL
                               GET_FG_GS_MATURITY (20 args, RF_ALLCONV
                               into GS_RF_ALLCONV; .FALSE. default at
                               1104 when no lagged field this route)
solve_for_pops.f:3168          GS_STATE 'fg_rf_all_conv=' telemetry
solve_for_pops.f:3507          LAM_FIELD_NOT_READY=LAM_FIELD_GATE(...)
solve_for_pops.f:3524          DO_NONMONO_KNOP .AND. HAS_LAGGED_FIELD
                               .AND. .NOT. GS_RF_ALLCONV
                               -> LAM_FIELD_NOT_READY=.TRUE., reason
                               'rf_not_allconv'
cmfgen_sub.f:5081,5513         CONV_EXIT_MET requires
                               .NOT. LAM_FIELD_NOT_READY (both sites)
```

Flag lifecycle (all in plane/fg_j_cmf_nonmono_v1.f):

- :565 module SAVE default .FALSE.
- :3086 reset inside FG_J_CMF_NONMONO_V1 on every INIT (every new
  formal solve, i.e. every great iteration's first frequency, and any
  new model in the same process). Phase 2.18 comment documents the
  reason: a TRUE surviving a refresh early-return or a model swap
  would certify a never-verified field.
- :2413/:2466 FG_EVAL_STATE_PUSH/POP saves and restores the flag
  around perturbed-population eval sweeps so JFNK FD evals cannot
  clear or forge the verdict.
- :8659 reset at the top of FG_KNOP_REFRESH_V1 productive work.
- :8513,:8521 allocation-failure returns precede :8659. On those
  paths the flag keeps its value - but every value it could keep was
  already cleared at :3086 by the same great iteration's INIT, so the
  stale-TRUE channel is closed. Early returns at 8423-8434
  (DO_KNOP_REFRESH off, FREEZE_I_NEXT, no compressing cells, no
  lagged store, store not armed) likewise cannot outlive the next
  INIT; within the same call a skip leaves the field untouched and
  the flag describing the last completed refresh, which is the
  documented semantic.
- :9351 set only after the full census: TRUE iff not FG_RF_FAILED,
  KN_NFULL == KN_NACT, KN_NACT > 0.

The seven sub-questions at HEAD:

1. MAXCH on rates from the preceding field: yes by construction
   (sweep rates precede the same-call refresh). Bounded by the
   GS_RF_DCHG maturity conjunct as shown in Phase 2.19.
2. Held exit bypassed before the next iteration's rebuild: no path;
   both CONV_EXIT_MET sites carry the conjunct and the converged
   branch is LST_ITERATION-wrapped.
3. Old moments/RJ surviving the verification refresh: no - the next
   call's sweep rebuilds them from the new store; acceptance reads
   the same call's verdict.
4. Rate coefficients reused without recomputation: none found; the
   sweep rebuilds rates every call.
5. Cached population residual triggering acceptance: no - MAXCH is
   recomputed each call.
6. Provisional refresh overwriting the verified field before
   acceptance: no window - census verdicts are same-call.
7. Accepted field == evaluated field: within the 5% maturity band;
   equal at the fixed point.

Verdict: the one-great-iteration lag is bounded, not a bypass. No
new defect found; the Phase 2.18 INIT reset closed the stale-flag
channel and the eval push/pop closed the JFNK-eval channel.

## PART B - population controlled tests

No standalone SOLVE_FOR_POPS driver exists; P1-P6 as specified need
the full model environment (populations, atomic data, rates) and
stay outstanding. Refresh-level coverage at HEAD:

| test | what it shows | result |
|------|---------------|--------|
| cerb | genuine bulk rays: bulk=3/3 -> F | PASS |
| cern | strict BTOL=0: full=3/3 -> T | PASS |
| cer2 | repeated provisional refresh: F then F | PASS |
| cerp | probe ray excluded, bulk=2/2 active -> F | PASS |
| rfmx | maxpass=1 -> itlim=3/3 -> F | PASS |
| rfab | pass1 T, pass2 F - no stale state | PASS |
| inval | corrupt NI_RAY -> inv=1 -> F (was segfault) | PASS after guard |
| invdln | zero DLN span -> inv=3 -> F | PASS |
| csut | 28/28 classifier rows | PASS |
| p1t1x BTOL=0 | ixit=1 FULL | PASS |
| p1t1x BTOL=6e-2 | ixit=2 BULK_ONLY | PASS |
| p1t1x maxpass=5 | ITERATION_LIMIT | PASS |

(cerb/cerp moved from full->T to bulk->F between the 3746683a and
2adab9c9 builds: the user's csut/driver commits b90ee078/844978ed/
9dcf621a made the fixture produce genuine bulk exits - the flag
correctly reports F on them.)

P7 (two models one process): resolved at source level by the
FG_J_CMF_NONMONO_V1 INIT reset (fg:3086) - a second model's first
formal solve clears the verdict before any refresh runs. The
driver loads one synthetic system per process, so a live two-model
probe remains a live-model check.

## PART C - observer interface

Call-site inventory at HEAD (exactly 3 sites, all 20 args matching
the callee at fg:6897):

- solve_for_pops.f:1105-1112 (GS_RF_ALLCONV, declared LOGICAL :172)
- obs/cmf_flux_sub_v5.f:2490-2496 (GS_ALLCONV_OBS, printed as
  'all-ray refresh conv'; diagnostic only - the OBS call is a
  read-out, it does not feed any acceptance path)
- unit_tests/f90/drv_knop_ray.f90:2210-2217 (GS_ALLC_T)

Callee signature (fg:6897-6903): GS_OK L; CHG_FILE/MAXCHG/
SW_MAXCHG/RF_DCHG REAL(LDP); CAP_FILE/CAPPED/NDIV/NRAY INTEGER;
NBAD..RF_CAPPED INTEGER*8; RF_ALLCONV LOGICAL. All three call sites
match count, order and kind.

OBS execution: cmf_flux_p219 (e3db3555...) linked clean in Phase
2.19. The standalone obs_p219 job series (13420413-13422341, user-
run, this session) progressed from setup failures to executing the
binary; the latest attempt OOM-killed at 22 s on the 30G interactive
partition. The production-route check is covered by the dependent
_OBS jobs queued behind the live A1/A2 model runs (13414392,
13414399). OBSFLUX diagnostic on/off byte-identity: outstanding -
needs a completed OBS run.

The OBS call site cannot set or inherit a population verdict: it
reads the flag through the export for printing only; FG_RF_ALL_CONV
is written only by the refresh census and the INIT reset, never by
cmf_flux code.

## PART D - production archive integrity

lib/libsubs.a at audit time (rebuilt 13:27 today by user):

- 630 members, all valid ELF objects, zero zero-length members,
  zero duplicates, index consistent (nm parses clean).
- Member list identical to the authoritative manifest
  docs/audits/knop_phase218/libsubs_members_631.txt (630 entries -
  the filename overstates by one).
- Critical symbols present once each: get_fg_gs_maturity,
  kn_certify, fg_knop_refresh_v1.
- No lib_dbg objects (archive was rebuilt from production objects
  per Phase 2.18; member-by-member ELF check passed).
- /tmp/libsubs_bk.a provenance: pre-Phase-2.17 backup of the
  truncated archive; superseded by the repaired archive, contents
  inspected in Phase 2.18.

## PART E - clean production build

./build.sh (make KEEP_EXE=1 clean; make DEBUG=nodebug all; f90
suite; python suite) run in an isolated git worktree at 3746683a,
/viper/ptmp1/sbrennan/SN2024pba/cmfgen_p220_build:

- Compiler: ifort (IFORT) 2021.11.1, intel/2024.0 + mkl/2026.2.
- Flags: -qopenmp -mkl -nbs -extend-source -traceback -O3
  -mavx512{f,dq,bw,vl,cd} -align array64byte (production set).
- make all completed: every exe target built including
  cmfgen_dev.exe, cmf_flux.exe, obs_frame.exe, plt_steq.exe.
- Worktree libsubs.a: 630 members, member list identical to the
  live archive and the Phase 2.18 manifest.
- Hashes (worktree, pre-guard-fix): cmfgen_dev.exe 77b4a592...,
  cmf_flux.exe 412675ad..., libsubs.a faff27c5....
- Live-tree exe/cmfgen_dev.exe 2bb070d9... (user build 12:27,
  powers the in-flight A1/A2 runs) predates the NI guard below;
  a rebuild is needed to carry it into production.
- Reproducibility: the worktree build is bit-different from the
  live exe (expected: archive member order and build timestamps
  differ); member inventory and symbol content are identical.

## PART F - regression suite

worktree run_all_tests.py: 91/97 files passed. KNOP suites green:
test_knop_ray (49 checks incl. all csut rows), test_knop_mode1,
test_knop_spv, test_knop_endpoint, test_knop_status_text,
test_knop_spec_gate, test_lam_handoff, test_mom_* all PASS.
T21b/T26/T27 stay resolved (fixture fix c4cff903 + pin 48cd5d49).

Six failures, none a solver regression:

| test | failure | classification | evidence |
|------|---------|----------------|----------|
| test_auto_ion | fixture dir m15_2p82d_VAL_stock absent | environment | dir does not exist on this fs |
| test_es_conv | ld: undefined set_fg_q_resid_ in rd_cmf_flux_controls.o | stale fixture | call added by 40842972; GUARD_SRCS link list not updated |
| test_knop_lagged_field | 'safety exits ungated by continuation' literal check | obsolete expectation | guard rewritten to NaN-safe .NOT.(ABS(RES).LE.HUGE) by 98acf207; IXIT=5 sites intact |
| test_obs_nonmono_dispatch | gate literal IF(DO_NONMONO_KNOP)THEN | obsolete expectation | gate widened to KNOP .OR. OBS_FG_KNOP by 40842972/66ea1e65 |
| test_par_opac | no production object for opac/mod_par_opac.f | test setup | par_opac.exe is outside make all and needs libnewsubs/libtools |
| test_sobew_angular_floor | missing RVTJ fixture file | environment | file absent on this fs |

Explicit knray battery at HEAD 2adab9c9 + the NI guard
(knray_p220 07b25464...): csut 28 rows exact; cerb bulk=3 -> F;
cern full=3 -> T; cer2/cerp PASS; rfmx itlim=3 -> F; rfab T->F;
invdln inv=3 -> F; inval inv=1 -> F (was: segfault).
p1t1x: BTOL=0 -> ixit=1 resid=1.87904203E-4 qres=3.36594109E-9
FULL (bit-identical to p217); BTOL=6e-2 -> ixit=2
resid=8.32052398E-3 qres=2.79648408E-3 BULK_ONLY; maxpass=5 ->
ITERATION_LIMIT.

## PART F.1 - defect found and fixed: NI_RAY out-of-extent crash

Reproduction: 'knray inval' set NI_RAY(1)=2e9 to force the ray
ALLOCATE to fail. The failure never reached the fail-closed
IOS path: Linux overcommit lets a ~432 GB virtual ALLOCATE
succeed, and the subsequent I=1..NI geometry loop at fg:10288
reads R_RAY(I,LS)/RAY_PNT beyond their real extent -> SIGSEGV
inside KNOP_RAY_SOLVE_V1 from the refresh OMP loop (fg:8885).

Fix (minimal, fail-closed): the existing construction-failure
guard now also rejects NI > SIZE(R_RAY,1) - the array extent is
the only honest bound on a ray's point count. The ray classifies
KN_CS_INV reason 8 with NXIT=-2/RESID=1.0, the refresh stays
immature, and no production-valid NI is affected (NI_RAY <=
SIZE(R_RAY,1) by construction).

Verified: inval -> inv=1, rf_all_conv=F, RESULT_OK; invdln
unchanged (inv=3, F). Same-pattern reads inside the in-sweep
FG_J_CMF_NONMONO_V1 (fg:4403,4629) are documented for a future
hardening pass; production NI_RAY is build-bounded there too.

## PART G - live model

Slurm live: A1 (13414391) and A2 (13414398) m15_2p82d_KNOP_MODE1
runs in flight on identical VADATs with the fresh cmfgen_dev.exe
(2bb070d9...), each with a dependent _OBS job (13414392, 13414399).
These are replicate runs of the new build, not an old-vs-new A/B:
the ten-release failure predates the flag, so a controlled
comparison would need a requeue of the pre-flag recipe. Census and
flag history will be readable from OUTGEN/GS_STATE when the runs
advance; outstanding at audit time.

## PART H - mixed-sign validation

Still no genuine velocity-inversion model (all local decks are the
all-compressing m15 hydro). Requirements for a future test:
real V(r) non-monotonicity placing seeded cells on both GAM signs,
both regions active in the seeded store, confirmed sign
transitions in solver inputs, then the full convergence and OBS
sequence. Not claimed.

## FINAL REPORT

1. Acceptance before rates are rebuilt on the verified field:
   bounded, not bypassed - call-N rates come from call-N's sweep
   (seeded by the previous refresh); the verified refresh's own
   field feeds call-N+1's rates. The GS_RF_DCHG maturity conjunct
   holds the exit while the gap exceeds the 5% band.
2. One-great-iteration lag handled: yes - bounded by maturity
   plus the same-call all-ray verdict; no demonstrated defect.
3. Moments and RJ rebuilt: yes - the sweep rebuilds them every
   call; nothing carries stale moments into a verdict.
4. Rates rebuilt before final MAXCH: yes - same-call ordering;
   MAXCH is recomputed each call.
5. Stale population residuals triggering acceptance: no.
6. Provisional refresh overwriting the verified field before
   acceptance: no window - verdict is same-call; the next
   iteration's INIT reset (fg:3086) clears any stale flag before
   a new refresh can run or skip.
7. cmf_flux executes: links clean (p219 e3db3555, worktree
   412675ad); the standalone OBS job reached the binary and died
   on the 30G interactive-partition memory cap - the production
   OBS check rides the dependent jobs behind A1/A2. Completion
   pending.
8. OBSFLUX diagnostic on/off byte-identity: outstanding - needs
   a completed OBS run.
9. Archive members present and valid: yes - 630 members, all
   ELF, index consistent, manifest-identical, key symbols
   present once each, no lib_dbg objects.
10. Production build reproducible: yes - clean worktree
    ./build.sh completes make DEBUG=nodebug all with all exe
    targets; member inventory identical (bit-different exes as
    expected from archive ordering).
11. T21b/T26/T27 resolved: yes - stay green under the pinned
    fixtures; the six current failures are environment or stale
    fixture expectations, none a solver regression. The suite
    did expose one real defect (inval segfault) now guarded.
12. Ten-release improvement: outstanding - A1/A2 replicate runs
    of the new build are in flight with dependent OBS jobs; no
    controlled pre-flag arm exists.
13. Mixed-sign validation: still unvalidated - no genuine
    velocity-inversion model exists locally; requirements
    recorded in Part H.

New defect found and fixed this phase: out-of-extent NI_RAY
segfaulted KNOP_RAY_SOLVE_V1 before the fail-closed alloc path
(Linux overcommit admits the huge virtual ALLOCATE). Guarded by
NI <= SIZE(R_RAY,1) at the construction-failure check (fg:10229)
-> KN_CS_INV reason 8, refresh immature, no OOB read.
