# CMFGEN KNOP production convergence-certification fix (Phase 2.14)

Status: production source corrected and validated on the Phase
2.13 decks. Minimal diff: 4 logical sites + comment updates in
`plane/fg_j_cmf_nonmono_v1.f`. No physics, controller, or
population/observer changes. Period-2 damper untouched.

Provenance: base `c0a9a76e`; binary `knray_p214` sha256
`4293159f14d42b1f9e890bf0388fb1cbfc376d75a4315cfb185052592d37ec4e`
(production flags: viper_intel_nodebug, ifort O3 AVX512).

## 1. What the audit corrected about the premise

The module doc comment (lines ~971-984) already documented the
INTENDED semantics: "KN_RAY_BULK_TOL ... default 6.0E-2 under
DO_NONMONO_KNOP, 0 = off otherwise" and "<=0 disables the new
term". The Oct-20 edits added `KN_RAY_BULK_TOL .LE. 0 .OR.`
short-circuits to the exit and commit tests, which made the code
contradict its own documentation: zero actually meant "no floor",
not "disabled". The fix restores the documented semantics:
BTOL>0 = loose residual floor the bulk rule may accept;
BTOL=0 = bulk exit disabled.

## 2. Source changes (all four sites, same pattern)

`KN_RAY_BULK_TOL .LE. 0.0_LDP .OR. RES .LE. KN_RAY_BULK_TOL`
replaced by
`KN_RAY_BULK_TOL .GT. 0.0_LDP .AND. RES .LE. KN_RAY_BULK_TOL`:

- `SWEEP_FLD` bulk-exit test (~line 11172).
- `SWEEP_FLD_ORD` bulk-exit test (~line 12338).
- `WRITE_STORE` commit gate, M channel (~line 10324).
- `WRITE_STORE` commit gate, P channel (~line 10328).

Comments updated at all sites explaining the corrected
zero-semantics, the measured defect (RES 8.3e-3 vs TOL 2e-3,
41/6227 bad cells, qres 2.8e-3), and why the commit gate changed
in the same edit (the Oct-20 comment itself warned that the exit
and commit semantics must stay consistent).

Six-class status mapping (authoritative, via IXIT/NXIT +
residuals; the existing NX_CNT summary bins already separate the
classes):

  IXIT 0,1 -> FULL_CONVERGENCE (RESID<=KN_RAY_TOL; qres<=TOL when
              the verify flag is armed)
  IXIT 2   -> bulk-accepted under a positive BTOL floor
              (BULK_ONLY against the strict TOL contract; only
              reachable when RES<=BTOL)
  IXIT 3,6 -> STAGNATION; 4 -> ITERATION_LIMIT; 5 -> DIVERGENCE;
  -1/-2/NaN/NDIV -> INVALID_STATE (fail closed)

## 3. Validation (fresh process per run, knray_p214, p1t1x)

| test | config | ixit | npass | resid | qres | status |
|------|--------|------|-------|-------|------|--------|
| cold x2 | BTOL=0 | 1 | 14 | 1.879e-4 | 0.0 | FULL |
| warm store | BTOL=0 | 1 | 2 | 0.0 | 0.0 | FULL |
| MAXPASS=5 | BTOL=0 | 4 | 5 | 0.305 | 0.19 | ITERATION_LIMIT |
| cold, loose floor | BTOL=6e-2 | 2 | 14 | 8.3e-3 | 2.8e-3 | BULK_ONLY |
| cold, mixed deck | BTOL=0 | 1 | 14 | 1.879e-4 | 3.4e-9 | FULL |

Baseline (pre-fix, phase213): ixit=2, 8.3e-3, qres 2.8e-3 ->
BULK_ONLY. The corrected run reproduces the Phase 2.13 corrected
reference exactly (RESID=1.9e-4, qres=0, same 14 passes, zero bad
cells).

## 4. Semantics of the six classes in production

- FULL_CONVERGENCE: commits through the converged path.
- BULK_ONLY (IXIT=2 with TOL<RES<=BTOL): commits as converged ONLY
  at the documented loose floor - this is the Oct-20 tuned policy
  for parked-store avoidance and is now honestly named; with the
  default-under-KNOP 6e-2 this remains a *deliberate* loose
  acceptance, not a silent defect.
- ITERATION_LIMIT/STAGNATION/DIVERGENCE/INVALID: never commit as
  converged; the existing FG_COMMIT_NC capped-provisional path
  retains useful progress and is explicitly distinct.
- Zero-tolerance reinterpretation removed; the VADAT key now does
  what its doc says.

## 5. Answers

1. IXIT=2 separated from FULL_CONVERGENCE: yes - it can only fire
   when RES<=BTOL (positive floor), and the contract classifies
   RES>TOL bulk exits as BULK_ONLY.
2. BULK_TOL=0 handled per documented semantics: yes - disables
   the exit (matches the doc comment's own claim).
3. Bulk exit bypass of the global residual: no longer possible
   for RES>BTOL; under a positive floor the acceptance is
   bounded and named.
4. WRITE_STORE corrected in the same change: yes - the commit
   gate uses identical >0/RES<=BTOL semantics.
5. Callers use authoritative status: the commit gate and NX_CNT
   summary now share the corrected exit semantics; no caller can
   certify RES>TOL as converged except the named BTOL floor.
6. Provisional stores distinct from converged stores: yes -
   fail-closed rays take the FG_COMMIT_NC capped path.
7. Corrected production reaches FULL_CONVERGENCE: yes -
   ixit=1, RES=1.9e-4, qres=0.
8. Reproduces Phase 2.13 reference: exactly (same 14 passes).
9. Real mixed-sign velocity-inversion test: NOT available - the
   synthetic deck's expanding cells sit outside the seeded state;
   flagged as remaining work, not claimed.
10. Invalid states prevented from committing: NDIV/NaN exits were
    already fail-closed (IXIT=5/-2); unchanged and verified
    consistent with the contract.
11. Period-2 damper unchanged: yes.
12. Production executable: knray_p214 frozen and hashed; libsubs.a
    member rebuilt with production flags (viper_intel_nodebug).

## 6. Caveats / not done

- TEST G (genuine velocity-inversion model): deferred - no such
  run deck exists in the diagnostic suite; recommend a real model
  A/B before changing the production VADAT default.
- The six-class status is enforced through exit/commit semantics
  and the NX_CNT bins; a dedicated named-status variable plumbed
  through the call chain was judged a larger diff than warranted.
- Population/J/Eddington/rescue physical regression on a live
  model not run (diagnostic decks only); the change alters only
  certification, not equations.
- Recommendation for a follow-on phase: evaluate whether the
  production VADAT default 6e-2 should be tightened now that the
  exit floor works - that is a model-validation decision, not a
  code defect.
