# CMFGEN KNOP production convergence-certification fix (Phases 2.14-2.15)

Status: Phase 2.14 corrected the `KN_RAY_BULK_TOL=0` semantics;
Phase 2.15 implemented the authoritative six-class certification
(`KN_CERTIFY`, `KN_CS_LS`/`KN_CS_REA`, `KN_CERT_CONV`/
`KN_CERT_DIAG`) with the independent Q verify evaluated before
the commit decision. See section 7 for the Phase 2.15
implementation, validation evidence and acceptance answers.

Phase 2.14 status: production source corrected and validated on
the Phase 2.13 decks. Minimal diff: 4 logical sites + comment
updates in `plane/fg_j_cmf_nonmono_v1.f`. No physics,
controller, or population/observer changes. Period-2 damper
untouched.

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
  (Superseded 20-Nov-2026: Phase 2.15 implemented exactly that -
  `KN_CS_LS`/`KN_CS_REA` now carry the authoritative per-ray
  status through commit, census and callers; see section 7.)
- Population/J/Eddington/rescue physical regression on a live
  model not run (diagnostic decks only); the change alters only
  certification, not equations.
- Recommendation for a follow-on phase: evaluate whether the
  production VADAT default 6e-2 should be tightened now that the
  exit floor works - that is a model-validation decision, not a
  code defect.

## 7. Phase 2.15 - production certification (implementation)

Provenance: commits `55174211` (six-class layer) through
`48cd5d49` (test-contract updates), plus Phase 2.17 (`96648c78`,
all-ray verdict `FG_RF_ALL_CONV`) and Phase 2.18 (`db184329`,
build integrity). Validated binary: `unit_tests/bin/knray_p215`
sha256
`3504af4a26a543fc8dfa59e6232684ac2a3faa69f695f65096e96e9460a20313`;
archive `lib/libsubs.a` sha256
`e66618945801bfc6d7fc28956cede808f2755b0ebc0e82923c019b1a414e47f1`.
Machine-readable results: `docs/audits/knop_phase215/cert_results.txt`.

### 7.1 Authoritative six-class status

`KN_CERTIFY` / `KN_CERT_CH` (module scope, pure) fold each
channel's exit code, production residual and independent verify
residual into a severity-ordered class; the ray status is the
worse of the two channels. `KN_CS_LS`/`KN_CS_REA` are SAVE
arrays that survive the refresh epilogue (allocated per refresh
like `FG_RAY_VC`), so callers read the same verdict the commit
gate used. Classes: `KN_CS_FULL`, `KN_CS_BULK`, `KN_CS_ITLIM`,
`KN_CS_STAG`, `KN_CS_DIV`, `KN_CS_INV`; reason codes: 0 ok,
1 verify failed, 2 verify missing/invalid, 3 non-finite
residual, 4 allocation failure, 5 unset/inconsistent exit,
6 diverged.

### 7.2 Ordering: solve, verify, classify, commit

Inside `KNOP_RAY_SOLVE_V1` the independent Q verify now runs
before the commit decision: the final production field is copied
to `FLD_Q`, exactly one scratch pass is applied per channel
(same transfer operator, coefficients and frequency ordering),
and `KN_Q_LS(LS)=MAX(Q_M,Q_P)` is recorded with per-channel
argmax detail in `KN_QA_*`. All host-shared sweep outputs the
verify pass would overwrite (divergence and clean-exit flags,
worst-cell records, clamp/negative counts, pass counts, rho
exports, period-2 state) are saved before the verify and
restored after, so production telemetry is untouched. The
commit gate then reads `ICLS` from `KN_CERTIFY`, not `IXIT`:
`FG_COMMIT_CONV .AND. ICLS .NE. KN_CS_FULL` routes every
non-certified exit to `WRITE_STORE_CAP` (bounded by
`FG_COMMIT_NC`) or store retention (NC=0); only `KN_CS_FULL`
reaches the ordinary `WRITE_STORE` + `FG_RAY_SEEDED=1` path.

### 7.3 Controls and telemetry

- `KN_CERT_CONV` (VADAT, default TRUE): requires a valid
  independent verify inside the FULL contract. With it off the
  status is still recorded honestly - a `BULK_ONLY` ray is never
  labelled FULL, matching the legacy commit behaviour on the
  store side only.
- `KN_CERT_DIAG` (VADAT, default FALSE): one `KN-CS` record per
  ray per refresh on LUER carrying ray, cls, rea, ixm, ixp, nx,
  resm, resp, res, qm, qp, tol, btol, np, nbd, nt and the store
  decision (`dec`, FG_RAY_VC semantics). Production output stays
  compact by default.
- Refresh epilogue prints a `cert census:` line alongside the
  raw `ray exits:` census and sets `FG_RF_ALL_CONV` (Phase 2.17)
  only when every production ray is `KN_CS_FULL`.

### 7.4 Status propagation

- Store gate: `KN_CS_FULL` only (7.2).
- Census/telemetry: `KN_CCLS(0:5)` + `cert census:` line,
  masked to production rays via `KN_PROD_MASK`.
- Callers: `FG_RF_ALL_CONV` is exported through
  `GET_FG_GS_MATURITY`; `solve_for_pops` folds it into
  `LAM_FIELD_NOT_READY`, so an unverified refresh holds the
  converged exit (Phase 2.17 audit, section 4). The flag is
  cleared on INIT and carried through eval push/pop
  (Phase 2.18 fixes).
- Probe rays (`KN_NOCP_PROBE`) classify normally but are masked
  out of the production census; their verdict stays 7.

### 7.5 KN_RAY_BULK_TOL defaults (Part 6 reconciliation)

- Module DATA: `KN_RAY_BULK_TOL=0` (off).
- VADAT parse (`rd_control_variables.f`): `KN_RAY_BULK_TOL_V=0`;
  `IF(DO_NONMONO_KNOP)KN_RAY_BULK_TOL_V=6.0E-2_LDP` - the
  documented production default under KNOP.
- Runtime after parsing: 0 outside KNOP, 6.0E-2 under KNOP,
  explicit user value wins. `BTOL=0` disables the bulk exit in
  both sweeps (Phase 2.14 semantics preserved); a positive floor
  can still fire `IXIT=2`, which now classifies `BULK_ONLY`
  unless the strict residual contract is also met.

### 7.6 Test evidence (knray_p215, all on the committed source)

- `csut`: 20-row classifier table, all six classes and reason
  codes exact (49/49 `test_knop_ray.py`, 9/9
  `test_knop_mode1.py`, `test_knop_lagged_field.py` static
  contract checks pass).
- `cerb` (BTOL=1e30 floor, TOL=-1e-30 disables residual exit):
  `bulk=3`, `cert census bulk=3`, all rays `CLS=1`, `VC=3`
  capped provisional, `SEED=0`, `RF_NCCAP=3` - no ordinary
  WRITE_STORE, no seeded mark.
- `cer0` (NC=0): `CLS=1`, `VC=2`, `RF_NCREJ=3`, store retained.
- `cer2` (two refreshes): `RF_NCCAP=6` cumulative - provisional
  commits stay capped and never promote.
- `cern` (strict): `CLS=0`, `VC=0`, `SEED=1` - genuine fixed
  point certifies FULL.
- `cerp` (probe on ray 3): probe `VC=7`, `KN-NOP delta=0`,
  production rays `CLS=1`/`VC=3`, census `bulk=2` (probe
  excluded).
- `rfmx`/`rfab` (Phase 2.17): all-ray flag reads F under
  iteration limits and resets between refreshes.
- Physical regression: no transfer-equation, relaxation or
  period-2 changes in the diff; Phase 2.18 verified archive
  integrity and rebuild reproducibility.

### 7.7 Answers to the Phase 2.15 acceptance questions

1. IXIT=2 as FULL_CONVERGENCE: no. The classifier maps
   `IXIT=2` to `KN_CS_BULK` unless `RES<=KN_RAY_TOL` and the
   verify also passes; a bulk-floored exit above TOL can never
   reach the `WRITE_STORE` branch.
2. BULK_ONLY reaching ordinary WRITE_STORE: no - the gate is
   `ICLS .NE. KN_CS_FULL -> capped provisional or retain`.
3. BULK_ONLY represented as provisional: yes - `CLS=1`,
   `VC=3` under NC>0 (`cerb`), `VC=2` retained under NC=0
   (`cer0`).
4. Q residual before final commit: yes - the verify block
   precedes `KN_CERTIFY` and the gate in `KNOP_RAY_SOLVE_V1`.
5. Missing/invalid Q producing FULL: no - `IREA=2` maps to
   `KN_CS_INV` under `KN_CERT_CONV` (csut c10/c11).
6. One authoritative six-class status: yes - `KN_CS_LS`/
   `KN_CS_REA` written once per ray, read by the gate, the
   census and the callers.
7. All callers use the status: yes - `solve_for_pops` consumes
   `FG_RF_ALL_CONV` (all-rays-FULL only); no caller reads
   `IXIT` for certification.
8. `KN_RAY_BULK_TOL=0`: preserved - bulk exit disabled in both
   sweeps.
9. Actual production default documented: module 0, VADAT
   6.0E-2 under DO_NONMONO_KNOP (7.5).
10. Provisional commits bounded by FG_COMMIT_NC: yes -
    `WRITE_STORE_CAP` only; `cer2` shows six capped commits
    over two refreshes, never a full commit.
11. Invalid states kept from full commitment: yes - NaN/Inf,
    unset exit, alloc failure and divergence classify
    `KN_CS_INV`/`KN_CS_DIV` (csut c08-c11, c17, c20; T12/T13).
12. Phase 2.13 deck FULL_CONVERGENCE: yes on the diagnostic
    deck (`cern`); the live-model residual evidence is in
    `CMFGEN_KNOP_BULK_TOL_LIVE_MODEL_VALIDATION.md` (ixit=1,
    RES=1.9e-4, qres=0).
13. Positive-BTOL produces BULK_ONLY: yes - `cerb` on the
    unit deck; live-model bulk exits sit at ~3e-2 under the
    6e-2 floor (phase215 outgen_ray_exits.txt).
14. Real mixed-sign inversion test: BLOCKED - no live deck
    with both signs in the seeded state exists (Phase 2.15
    audit section 6, item 1-2). The GAM-crossing unit test
    (T15) covers the operator, not a full model.
15. Period-2 damper unchanged: yes - the verify saves and
    restores the p2 state; no diff in the damper code.
16. Physical equations unchanged: yes - only certification
    ordering, classification and commit routing changed.
17. Frozen and hashed executable: yes - `unit_tests/bin/
    knray_p215` + `obj_knray_p215_src`, hashes in
    `docs/audits/knop_phase215/cert_results.txt`.

### 7.8 Still open (carried to later phases)

- Live-model A/B (population fields, spectra) requires slurm
  runs; Phase 2.15-2.16 audits document this as outstanding.
- Genuine mixed-sign velocity-inversion model: BLOCKED, not
  claimed.
- `KN_CERT_CONV` runtime cost: one scratch pass per channel
  per ray per refresh (~2/NPASS extra sweeps); measured zero
  on the unit deck, unmeasured on a live model.
