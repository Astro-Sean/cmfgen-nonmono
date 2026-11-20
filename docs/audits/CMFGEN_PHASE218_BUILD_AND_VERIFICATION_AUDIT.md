# Phase 2.18 - build-integrity and final-verification audit

Provenance: HEAD b953b83e (end of Phase 2.17) plus this phase's
fail-closed source fixes (committed separately). Working tree
carried concurrent user edits to unit_tests/f90/drv_knop_ray.f90
and unit_tests/python/test_knop_ray.py throughout the audit.

Compiler: ifort 2021.11.1 (intel/2024.0.2). Production flags
(viper_intel_nodebug): `ifort -qopenmp -mkl -nbs
-diag-disable=10448,10412,10434 -extend-source -traceback -O3
-mavx512f -mavx512dq -mavx512bw -mavx512vl -mavx512cd
-align array64byte -I mod/ -module mod/`.

## 1. Archive findings (Parts A, B)

lib/libsubs.a at Phase-2.17 state (sha 5dfc4d92, knray_p217 link
source) was damaged two ways, neither caused by this work:

- **Truncation.** Only 520 of 631 members were recoverable. The
  fuller member set (631, e.g. ave_level_energy.o, colcool_sl_v4.o,
  fg_j_cmf_v12/v13.o, mom_jrel_*.o) reappeared in the archive image
  left behind by a later crashed `ar` atomic write
  (lib/libsubs.a.tmp). knray_p217 did not need the missing 111
  (they are routines outside the driver's dependency set), but a
  production `cmfgen` link would.
- **15 corrupt members** (tools/*.f + newsubs/var_op_v6.f):
  garbage payloads or zero length; `count_occur.o` was 0 bytes.
  The same corruption exists byte-identically in the frozen
  archive exe_frozen/508b12ca..., so it predates Phase 2.17 -
  most likely an older `ar`/`ranlib` fault that wrote index data
  into member slots. lib_dbg copies of the same members are valid
  ELF with debug info.

Impact on knray_p217: none of the 15 is in the driver's 26-source
dependency set except set_kind_module, whose module-init symbol
the link resolved from a sibling object (tools/*.o are deleted
after archiving by the dir Makefile, which is why the resolver
then fell back to the archive member and hit the corruption).
All other archive members are valid ELF; `nm --defined-only`
audit shows only 2 duplicate-definition cases
(cmf_form_sol_v2/v3 share module+entry names by design;
fdg_f_to_s_ns_v1/v2 share `do_sl_adjustement_`) - pre-existing
upstream versioning, resolved in production link order, neither
in the driver set.

Repair performed (Part C):
- All 15 corrupt members rebuilt from source with production
  flags (ELF verified).
- Archive rebuilt from all 631 members + updated
  fg_j_cmf_nonmono_v1.o via `ar crs` + `ranlib` (safer than
  incremental `ar r` on this filesystem, which repeatedly
  bus-errored mid-atomic-write leaving st* temp files; those
  debris were removed).
- Freshness: `stale_archive_members()` reports 0 stale members;
  the two Phase-2.17 members are newer than source.
- lib/libmod.a, libmy_blas.a, libmy_lapack.a, libpgplot_stubs.a,
  libpgplt.a: all members ELF-valid; libdisp.a, libspec.a also
  present and used by other exes.
- Backup of the pre-repair archive preserved outside lib/ at
  /tmp/p218_full/libsubs_pre218.a (a stray libsubs_pre218.a inside
  lib/ trips the whole-archive freshness sweep - documented).

## 2. Clean-build comparison (Part C, F)

knray_p218 (sha ee03b158...) links 26 sources: 23 production
objects resolved from the repaired archives + 3 test objects,
plus the same MKL/link line as production. All extracted members
are ELF-valid with recorded provenance.

Numerical comparison vs frozen knray_p217 is confounded by a
concurrent driver change (uncommitted drv_knop_ray.f90): the
user's in-flight fixture redesign moved the cer* cases to
KN_RAY_TOL=-1.0E-30 + KN_RAY_BULK_TOL=1e30, which now genuinely
produces `bulk` census classes (previously the deck converged to
exact 0 residual and certified FULL). Where the cases overlap the
old fixture, outputs match; where the fixture changed, the
differences trace to driver settings, not object provenance.
All driver-relevant archive members are byte-identical between
the pre-repair and repaired archives except fg_j_cmf_nonmono_v1.o
(this phase's edits) and set_kind_module.o (rebuilt).

## 3. Final-verification source trace (Part D)

Per great iteration inside CMFGEN_SUB:

  10000 sweep ends        -> rates/RJ already accumulated on the
                             sweep's evolving lagged store
  cmfgen_sub.f:3316       -> FG_KNOP_REFRESH_V1() writes
                             I_P/I_M_DEPTH_NEXT, runs the cert
                             census, sets FG_RF_ALL_CONV
  cmfgen_sub.f:4603       -> SOLVE_FOR_POPS: GET_FG_GS_MATURITY at
                             solve_for_pops.f:1106 exports the flag
                             (same call, post-refresh - fresh)
  solve_for_pops.f:3525   -> fold-in: rf_not_allconv forces
                             LAM_FIELD_NOT_READY=T
  cmfgen_sub.f:5081/5513  -> CONV_EXIT_MET requires
                             .NOT. LAM_FIELD_NOT_READY

Established from control flow:

- No stale flag: the flag is written only at census (fg:9305) and
  read only after the same call's refresh. A subsequent BULK
  refresh flips it F in the same call before the exit test - no
  window where a provisional overwrite can precede acceptance.
- No stale residual: MAXCH/MAXCH_APP are computed inside the same
  call's population update on rates built during the same call's
  sweep.
- Rate-evaluation lag (honest caveat): a call's rates are built on
  the field as it stood through that call's sweep - i.e., the
  store seeded from the PREVIOUS refresh, before this call's
  refresh committed. Acceptance can therefore fire on the first
  verified refresh even though that call's rates predate the
  verified field. If the verified refresh materially changed the
  field, the NEXT call's rates see it; a model that exits at the
  first verified refresh has implicitly accepted the fixed point
  only when prior seeds were already close (the nilpotent
  compressing deck reaches it exactly). Full strictness would need
  the flag true on two consecutive calls; documented as a
  limitation, not implemented (would change termination timing).
- Eval sweeps: FG_EVAL_STATE_PUSH_V1/POP_V1 now carry
  FG_RF_ALL_CONV with the rest of the census (this phase's fix) so
  an eval span cannot leave a verified field reading unverified.
- Early returns in FG_KNOP_REFRESH_V1 fall into two classes:
  pre-reset guards (DO_KNOP_REFRESH/FREEZE_I_NEXT/HAS_COMPRESSING/
  allocation/IO) return with the field untouched, where the
  previous verdict still describes the unchanged store; post-reset
  failures set FG_RF_FAILED=T and the flag stays F (fail closed).
  The INIT-time reset added this phase additionally clears the
  flag at every new formal solve, so a skipped refresh can never
  inherit TRUE across calls or across two models in one process.

## 4. Source changes made this phase (all fail-closed hardening)

1. plane/fg_j_cmf_nonmono_v1.f
   a. `FG_RF_ALL_CONV=.FALSE.` added to the INIT block of
      FG_J_CMF_NONMONO_V1 (~3038): clears the verdict whenever a
      new formal solve begins; closes the stale-TRUE window across
      skipped refreshes and multi-model processes.
   b. FG_EVAL_RFAC added to the PUSH/POP telemetry set: the flag
      now survives JFNK eval spans like every other census field.
2. obs/cmf_flux_sub_v5.f - **interface mismatch fix**: the OBS
   maturity diagnostic still called GET_FG_GS_MATURITY with the
   pre-Phase-2.17 17-arg signature (the callee now takes 18). On a
   routine with no explicit interface this is stack corruption on
   the OBS route - the same class of bug as the Phase 2.17 driver
   segfault. Added GS_ALLCONV_OBS as the trailing arg and printed
   `all-ray refresh conv` in the OBS census block.
3. lib/libsubs.a - repaired archive (documented above).

No transfer physics, population equations, tolerances, or the
relaxation controller were touched.

## 5. Test results (Part E, F)

knray driver battery on the clean build (current driver fixture):

  csut  20/20 classifier rows exact (all six classes + reasons)
  cerb  bulk=3  -> rf_all_conv=F   (TEST 2 - now a real BULK exit)
  cern  full=3  -> rf_all_conv=T   (TEST 1)
  rfmx  itlim=3 -> rf_all_conv=F   (TEST 3)
  rfab  T pass1 -> F pass2         (TEST 6 stale-state)
  cerp  bulk=2/2 active, probe out -> F   (TEST 5 inverted: probe
        excluded AND a bulk active ray correctly fails the flag)

test_knop_ray.py suite: 46/49 PASS. Failures T21b/T26/T27 are in
the user's in-flight fixture redesign (FG_RAY_SEEDED=1 now arms
the seed-fill waiver for legacy cases, so the jump veto no longer
fires where the fixture intends; the tests' expectations predate
the waiver being observable). Production solver paths unchanged;
the failures are driver-fixture semantics, not solver regressions.

Part E cases 1/2/4/5/6 (population-level state machine) and TEST 7
(two models, one process) require a live model run - outstanding,
as is TEST 8 (early-return) beyond the source-level proof above.

## 6. Live-model A/B and OBSFLUX (Part G)

Not run - requires the Slurm production model (hours-days).
The flag's release/freeze path is unchanged, so the ten-release
hypothesis is testable but unverified. OBSFLUX reproducibility:
the OBS caller signature fix makes the previous binary state
(strictly) incorrect; identity on OBSFLUX must be re-verified on
the next obs build/run. exe/cmf_flux.exe not relinked here (not
overwriting run products); obs/cmf_flux_sub_v5.o rebuilt in place.

## 7. Answers

1. libsubs.a was NOT consistent at p217: truncated (520/631
   recoverable) and carrying 15 pre-corrupted members. Repaired
   archive: 631 members, all ELF, fresh index, 0 stale.
2. Debug objects in production archive: none found (prod members
   are non-debug ELF; debug-compiled objects live only in
   lib_dbg/).
3. Critical symbols from unexpected objects: no - the only
   ambiguity was set_kind_module, resolved from a sibling .o.
4. knray_p217 reproducible: YES at the functional level - the
   clean build produces identical test outputs on the unchanged
   fixture subset; bit-different binary expected (member order,
   repaired members).
5. Verified refresh forces moment recomputation: indirectly -
   the next sweep rebuilds moments on the updated store; within
   the same call the rates predate the refresh (lag documented).
6. RJ/rate recomputation: yes, next great iteration (same lag).
7. Stale population residual permitting acceptance: no - MAXCH is
   evaluated inside the same call on same-call rates.
8. Later provisional refresh overwriting before acceptance: no
   window - flag and exit test are same-call, recomputed per
   refresh.
9. State-machine tests: all driver-realizable cases pass; the
   population-level cases and live A/B are outstanding.
10. Ten-release cycle in live comparison: not run.
11. OBSFLUX reproducible: untested locally (needs a model); the
    signature fix removes a real stack-corruption path.
12. Mixed-sign validation: still none - no genuine velocity
    inversion model exists.

## 8. Remaining limitations

- Rate/verification one-call lag (section 3) - strict two-call
  verified-field consistency is a policy decision, not a bug.
- The 631-member archive's fuller member set vs the 520-member
  truncated image: the repaired archive is the authoritative
  production library; per-member source-version audit beyond the
  freshness check was not repeated (all members predate this work
  except the two rebuilt Phase-2.17 objects and the 15 rebuilt
  tools members).
- OBSFLUX identity pending a rebuilt obs exe + model run.
- T21b/T26/T27 test expectations await the user's fixture
  conclusion (SEED-waiver observability).
