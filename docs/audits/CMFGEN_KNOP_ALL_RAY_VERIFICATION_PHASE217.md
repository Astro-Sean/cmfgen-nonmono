# CMFGEN KNOP all-ray convergence verification (Phase 2.17)

Provenance: HEAD 2b115a97 plus the in-tree Phase-2.15 six-class
certification layer (KN_CERTIFY / KN_CS_* / KN_CERT_DIAG, uncommitted
user work integrated here). Binary: knray_p217
(sha256 babede9cf848e262a3bb37b937bed0ea741da59328c6f3726a1c4c14e2b94176).

## 1. Status source audit (Part A)

- Per-ray raw exit: NXIT -> NXT_LS -> NX_CNT (-2 build fail, -1
  unsolved, 0 exact, 1 tol, 2 bulk, 3 stall, 4 cap-stop, 5 div,
  6 plat). Raw codes only; not the certification contract.
- Authoritative status: KN_CS_LS(LS), six severity-ordered classes
  (plane/fg_j_cmf_nonmono_v1.f:1274):
    KN_CS_FULL=0  both channels RES<=TOL and verify residual<=TOL
    KN_CS_BULK=1  positive KN_RAY_BULK_TOL floor met, strict RES unmet
    KN_CS_ITLIM=2 iteration limit
    KN_CS_STAG=3  stall/plateau or verify failed
    KN_CS_DIV=4   finite numerical runaway
    KN_CS_INV=5   non-finite residual, missing/invalid verify,
                  unset exit, alloc failure, missing field
- KN_CS_LS is reset to KN_CS_INV/rea=5 at the head of every refresh
  ray-loop init (line ~8577), so stale or never-classified rays read
  INV by construction. A previous refresh's status cannot survive.
- Active rays: KN_PROD_MASK = .NOT. KN_NOCP_PROBE_RAY -- only the
  no-coupling probe ray is inactive; all others must certify.
- Both channels fold into one class via severity MAX inside
  KN_CERTIFY; NDIV forces DIV. The independent verify pass runs when
  KN_CERT_CONV (default T) or KN_Q_RESID is set; a missing verify is
  INV under KN_CERT_CONV.
- Census site: KN_CCLS counts over KN_CS_LS masked by KN_PROD_MASK,
  after the OMP ray loop completes (~9273).

## 2. The flag (Part B)

`FG_RF_ALL_CONV` (module SAVE, fg file ~558):

    reset .FALSE. at refresh init next to KN_CS_LS=KN_CS_INV (~8587)
    set at census end (~9284):
        FG_RF_ALL_CONV = (.NOT. FG_RF_FAILED) .AND.
                         (KN_NFULL .EQ. KN_NACT) .AND. KN_NACT .GT. 0

Derived only from per-ray classes; no residual threshold. An
incomplete census (any ray still INV), a failed refresh stage
(FG_RF_FAILED), or zero active rays all read FALSE. Printed as
`rf_all_conv= T/F full= n active= m` each refresh (one line).

Export: new trailing dummy arg RF_ALLCONV on GET_FG_GS_MATURITY
(~6864). Callers: solve_for_pops (updated), drv_knop_ray
DRV_MATURITY (updated). No other callers exist.

## 3. Maturity vs convergence (Part C)

GS_MAX_ALL / GS_MATURE_TOL / the release-freeze controller are
untouched: provisional fields still drive intermediate iterations,
satisfying the "progress not certification" requirement.

## 4. Final acceptance (Part D)

CONV_EXIT_MET in cmfgen_sub.f (~5510) already requires
.NOT. LAM_FIELD_NOT_READY for the converged model exit. The flag is
folded into that verdict inside SOLVE_FOR_POPS after all final-call
state is set (solve_for_pops.f ~3503):

    IF(DO_NONMONO_KNOP .AND. HAS_LAGGED_FIELD .AND.
       .NOT. GS_RF_ALLCONV) LAM_FIELD_NOT_READY=.TRUE.
       (+ ITER_REASON 'rf_not_allconv')

Consequences: when the last refresh left any production ray below
KN_CS_FULL, the converged exit cannot fire; the next great iteration
re-runs the refresh and re-evaluates the population demand on the
new field (the required re-validation). A stale TRUE cannot persist
because the flag is recomputed at every refresh census. If the field
never verifies, the run terminates via the existing iteration
limits / release watchdogs as non-convergent -- fail closed.

Machine-readable: GS_STATE gains `fg_rf_all_conv= T/F` next to
`mom_refresh_conv`.

## 5. Bug found in the way (same edit)

The Phase-2.15 `cert census` WRITE used format `(A,6(A,I5))` = 13
descriptors for 12 arguments; the 6th label landed on an I5 slot and
the line crashed with a runtime format error on every refresh that
reached it. Fixed to `(A,I5,5(A,I5))` with precomputed counts.

## 6. Tests (knray_p217, systems deck, OMP serial)

| case | setup | census | rf_all_conv | verdict |
|------|-------|--------|-------------|---------|
| csut | 20-row KN_CERTIFY table | - | - | all 20 classes/reasons exact |
| cern | BTOL=0 refresh | full=3 | T | TEST 1 pass |
| cerb | BTOL=1e30 | full=3 | T | deck converges below TOL, FULL not BULK |
| rfmx | MAXPASS=1 | itlim=3 | F | TEST 3 pass |
| rfab | pass1 full, pass2 MAXPASS=1 | T then F | T->F | TEST 6 pass (no stale T) |
| cerp | ray3=probe (inactive) | full=2/2 | T | TEST 5 pass (exclusion) |

TEST 2/4 (bulk/invalid ray -> F): no synthetic-deck case produces a
bulk or invalid ray since the deck converges or pass-caps; covered
by construction (the flag is the KN_CCLS(0)==NACT count; any
non-FULL class already proven to flip it via itlim) and by the
csut classifier table for those classes. The production p29m deck
exercises bulk exits but is only wired to single-ray solves in this
driver -- see limitations.

TEST 8 (population revalidation after verified refresh): the
fold-in at LAM_FIELD_NOT_READY makes every unverified refresh hold
the exit so the next iteration re-solves populations on the new
field -- structurally guaranteed, not run-verified (needs a live
model run).

TEST 9 (production model vs p214 baseline): the p29m production
deck exercises single-ray cert only (P1X_RAW/P1X_STATUS); the
refresh-level flag needs a live model run -- not claimable here.

TEST 10 (diagnostics on/off identity): flag is computed whether or
not KN_CERT_DIAG prints; the KN-CS per-ray records are the only
verbosity change and write to LUER, never to the field path.

## 7. Controlled ten-release comparison (Part G)

Not run -- requires the full live model (hours-days on slurm).
The mechanism prediction: with BTOL=6e-2 the field now reports
rf_all_conv=F at ~3e-2 residual, so the maturity band can no longer
silently certify; releases would instead be judged against a
provisional field while final acceptance is held until a verified
refresh. This is a hypothesis pending the live run in Phase 2.15-7.

## 8. Build/regression safety

- knray_p217 (sha256 above) links patched lib/libsubs.a members
  (fg_j_cmf_nonmono_v1.o, solve_for_pops.o rebuilt with production
  flags; ar bus errors on viper fs repaired via index rebuild and
  member re-add; archive member count 520+).
- libsubs.a also received members it was missing (term_io,
  verbose_info, get_moms_non_rel, rd_store_log, mom_j_cmf_nonmono_v1)
  -- the earlier knray builds had linked them from lib_dbg.
- Frozen baselines preserved: knray_p210..p214 untouched.
- No changes to transfer equations, population rates, KN_RAY_TOL,
  GS_MATURE_TOL, INEXT_RELAX, the relaxation controller, J_RESCUE,
  RESC_JMAX, or the observer solver.

## 9. Limitations

- Flag is refresh-scoped, not store-scoped: a provisional store
  remains physically present and consumed; the flag only prevents
  certification. Part 3's three-way store distinction is realised
  as a verdict, not a per-store tag.
- Live-model behaviour (Part G), population-level A/B and the
  genuine mixed-sign deck are still open, exactly as in Phase
  2.15-2.16.

## 10. Independent validation addendum

Verified end-to-end on knray_p217 (sha256
21c1028a2889c8d54eee661e02e1d131d8b4272675851bc4a391914426ef63e1),
built from the working tree including the uncommitted Phase 2.16
hardening (INVR invalid-state reasons, negative-residual and
BULK verify checks in KN_CERT_CH):

- csut: 20/20 classifier rows exact.
- cerb/cern/cer2: full=3/3 -> rf_all_conv=T.
- cerp: probe excluded, full=2 active=2 -> T.
- rfmx (maxpass=1): itlim=3/3 -> F.
- rfab: pass1 T then pass2 F -> stale state cleared.
- p1t1x production deck: BTOL=0 -> ixit=1 RESID=1.87904203E-4
  qres=3.4e-9 FULL; BTOL=6e-2 -> ixit=2 RESID=8.32052398E-3
  qres=2.796e-3 BULK_ONLY; maxpass=5 -> ITERATION_LIMIT.
- Census WRITE format bug found and fixed: '(A,6(A,I5))'
  misaligned ' inv=' onto an I5 slot (runtime error 61);
  corrected to '(A,I5,5(A,I5))' with staged counts.

Machine-readable outputs: docs/audits/knop_phase217/.

## 11. Optional-argument portability fix (KN_CERTIFY)

The Phase 2.16 INVR gate was written as a combined condition:

    IF(PRESENT(INVR) .AND. INVR .NE. 0)THEN

Fortran does not guarantee short-circuit evaluation of .AND., so
a compiler may evaluate INVR while it is absent. Replaced with
explicitly nested IFs: PRESENT first, then the value test inside
the same block. Classification priority is unchanged: a present,
non-zero INVR still short-circuits to INVALID_STATE before the
NDIV divergence fold and the two-channel certification fold.

Coverage (csut, now 29 rows): c01 calls without INVR and returns
FULL; new row c29 calls without INVR with NDIV=1 and returns DIV;
c26 keeps INVR=7 winning over NDIV=1; c27/c28 keep INVR=8/9.
Rows c21-c28, added with the Phase 2.16 hardening but previously
unchecked by test_knop_ray.py, are now asserted in _cert_exp.

Checked build: fg_j_cmf_nonmono_v1.f compiled with ifort
`-check bounds,uninit,arg_temp_created,output_conversion -warn all
-O2` and linked as knray_chk ran the full csut table with no
runtime diagnostics. Full suite after the fix: 49/49.
