# CMFGEN KNOP bulk-exit convergence-certification audit (Phase 2.13)

Status: defect characterized end to end; diagnostic convergence
contract implemented driver-side; A/B experiments show a residual
floor on the bulk exit restores true convergence with no physics
or controller changes.

Provenance: HEAD aa6bea81 + uncommitted driver modes; binary
`knray_p213`; decks `p1t1s_in.txt` (all-compressing, NX=1680) and
`p1t1s_mixed.txt` (18 expanding depths, NX=1248).

## 1. IXIT exit-path audit

`SWEEP_FLD`/`SWEEP_FLD_ORD` assign IXIT per direction; the codes:

| IXIT | class | condition |
|------|-------|-----------|
| -1   | unset | initial value |
| 0/1  | exact / residual-tol | RES < KN_RAY_TOL |
| 2    | **bulk** | NPO>=6 AND NB_O<=MAX(50,NT_O/1000) AND (BULK_TOL<=0 OR RES<=BULK_TOL) |
| 3    | stall | RATIO_OM>KN_RAY_STALL twice (NPO>=8), noncontracting |
| 4    | pass cap | NPO>=KN_RAY_MAXPASS |
| 5    | diverged | NDIV>0 or RES>1e3*BEST_RES |
| 6    | plateau | W6 rule, three consecutive ratios in [0.95,1.05] |
| 9    | early | auxiliary early path |

Callers: `KNOP_RAY_SOLVE_V1` maps IXIT_M/IXIT_P -> NXIT
(divergence wins, else the residual-setting channel's class);
`FG_KNOP_REFRESH_V1` bins exit classes into the refresh summary;
the WRITE_STORE commit gate (lines ~10320-10360) admits a commit
when `IXIT in {0,1,2}` AND `RES<=KN_RAY_TOL` OR `IXIT==2` with
`BULK_TOL<=0 OR RES<=BULK_TOL`. The in-source comments (17-20 Oct
edits) document that the gate was deliberately opened for
BULK_TOL<=0 bulk exits to avoid parked stores: **an IXIT=2 exit
with the default BULK_TOL=0 is treated as full convergence for
commit purposes**, and the summary counts it under the
"converged" bulk class.

## 2. Bulk-exit mathematics

- `NTOT_CUR` = every updated cell this pass; `NBAD_CUR` = cells
  whose normalized update T1/DCELL exceeds `KN_RAY_TOL`
  (per-cell 2e-3 here). Counts at ~lines 11723-13265.
- `NB_O<=MAX(50,NT_O/1000)` is a **bad-fraction heuristic**: at
  most 50 cells (or 0.1%) still outside the per-cell tolerance.
  It measures a *subset* convergence, not the global residual RES
  (the ray-level normalized max/L2 metric) nor the independent
  one-pass verify residual KN_Q_RESID.
- `KN_RAY_BULK_TOL=0` semantics (traced, not inferred): the exit
  condition `(BULK_TOL<=0 OR RES<=BULK_TOL)` short-circuits to
  TRUE, so **zero = no residual floor = unconditional acceptance
  once the bad count is small** - not "disabled" and not "exact
  zero". Same short-circuit appears in the commit gate.
- Therefore YES: IXIT=2 can certify RES=8.3e-3 with qres=2.8e-3
  under KN_RAY_TOL=2e-3 - demonstrated live (below).

## 3. Convergence contract (diagnostic)

Statuses computed in `cert` mode from returned quantities only:

- FULL_CONVERGENCE: NXIT in {0,1} AND RESID<=TOL AND qres<=TOL
  (qres is the independent one-pass verify residual).
- BULK_ONLY: NXIT==2 but RESID>TOL or qres>TOL.
- ITERATION_LIMIT: NXIT==4. STAGNATION: NXIT in {3,6}.
- DIVERGENCE: NXIT==5. INVALID_STATE: NDIV>0 / NaN / unset.

Contract rule: a bad-fraction exit is a *subset* certificate;
FULL requires both the global ray residual AND the independent
verify residual to pass the requested tolerance.

## 4. A/B experiments (fresh process per case, cold start)

| case | config | ixit | npass | resid | qres | status |
|------|--------|------|-------|-------|------|--------|
| A prod deck | BULK_TOL=0 (default) | 2 | 14 | 8.3e-3 | 2.8e-3 | **BULK_ONLY** |
| D prod deck | BULK_TOL=1e-3 | 1 | 14 | 1.9e-4 | 0.0 | **FULL_CONVERGENCE** |
| A mixed | BULK_TOL=0 | 2 | 14 | 8.3e-3 | 2.8e-3 | **BULK_ONLY** |
| D mixed | BULK_TOL=1e-3 | 1 | 14 | 1.9e-4 | 3.4e-9 | **FULL_CONVERGENCE** |
| A warm seed | BULK_TOL=0 | 1 | 2 | 2.9e-12 | 0.0 | FULL_CONVERGENCE |

NBAD=41 of NTOT=6227 cells at the bulk-exit pass - inside the
50-cell threshold, so the fraction rule fired while the global
residual was 4x tolerance.

The B/D/E/F split (exit policy vs controller): the bulk-exit
defect is fully isolated - with BULK_TOL=1e-3 the SAME adaptive
controller converges to qres~0 in the same pass count, so no
controller change is needed to fix the certification defect.
(Structure-aware omega remains a separate optimization from
Phase 2.12: it converges at ~13 passes under fixed omega=1; the
exit-policy fix is orthogonal.)

## 5. Period-2 damper audit (for Phase 2.14, not implemented)

The NOSC>=2 branch damps omega*0.7 on consecutive sign flips
(NPO>=4, omega>0.45). On a nilpotent trajectory the residual can
alternate sign during the compositional transient, so the damper
can fire before convergence is even assessable. Recommendations:
(a) hold the damper off during a proven-DAG transient window
(chain-length passes); (b) require sustained sign flips on
comparable magnitudes, not merely NOSC>=2; (c) never re-damp
after a valid FULL_CONVERGENCE certificate; (d) period-2 damping
is appropriate for genuinely cyclic operators (mixed-sign decks)
where alternating modes are real eigenmodes.

## 6. Regression/safety

- No production source edits this phase; all classification is
  driver-side. Production defaults, transfer equation, boundary
  conditions, KN_COUPLING_OK, J_RESCUE/RESC_JMAX, population and
  observer solvers untouched.
- Mixed deck result identical to the all-compressing deck
  (expanding cells do not enter the seeded state on this
  geometry) - mixed-sign validation is structurally present but
  weak on this deck; flagged for a genuine velocity-inversion
  model test.

## 7. Answers

1. IXIT=2 = "bad-fraction exit": <=50 cells (0.1%) failed the
   per-cell update tolerance. A subset criterion, not a global
   residual certificate.
2. Treated as success by: the WRITE_STORE commit gate (admits
   with BULK_TOL<=0) and the refresh summary's converged bin.
   The maturity gate still reads RES, so downstream it is not
   fully hidden - but the commit itself writes the unconverged
   field.
3. NB_O = per-pass count of cells whose normalized update
   exceeds KN_RAY_TOL.
4. BULK_TOL=0 = "no residual floor" (short-circuit TRUE), not
   disabled, not exact-zero.
5. Can certify unconverged: YES - demonstrated: 8.3e-3 vs 2e-3.
6. Valid convergence requires global RESID<=TOL AND independent
   qres<=TOL (the verify-pass residual), not a subset count.
7. Corrected diagnostic exit prevents false convergence: YES -
   classifier marks it BULK_ONLY; with BULK_TOL=1e-3 the real
   solve continues to FULL_CONVERGENCE (ixit=1, qres~0).
8. Preserves valid convergence: YES - warm/nilpotent cases still
   report FULL_CONVERGENCE (qres=0).
9. Mixed-sign safe: yes on this deck (same classification), with
   the caveat that expanding cells sit outside the seeded state.
10. Phase 2.14 recommendation: (a) set a residual floor for the
    bulk exit - either KN_RAY_BULK_TOL > 0 at a documented value
    (1e-3 = TOL/2 worked here; justify against production
    tolerance) or require RES<=KN_RAY_TOL at IXIT=2; (b) surface
    the six-class status in the refresh summary/log so BULK_ONLY
    is visible; (c) separately evaluate the period-2-damper
    integration per Part 8. Do NOT change the exit class meaning
    silently - the commit gate at ~10320 must be updated in the
    same edit so commit semantics and exit semantics stay
    consistent (the 20-Oct comment warns these two were
    inconsistent before).
