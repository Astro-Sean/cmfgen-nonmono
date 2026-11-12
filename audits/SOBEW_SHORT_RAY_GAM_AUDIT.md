# SOBEW-family short-ray stale-GAM audit

Base commit for the fix work: `cd0f4374f1be07645a6a79f205a48f4285ea864b`
("knop iso: decisive arm split - M/C converging, R/R2 stuck").
Current HEAD at time of writing: `247c9ad68c146b418d35dc7be890a9cba314f399`
- three knop commits (`control_variable_mod.f`,
`rd_control_variables.f`, `fg_j_cmf_nonmono_v1.f` only) that do not
touch the sobew/test files. The fix diff is unchanged by them.

## 0. Repository-state reconciliation

The requested baseline `a029fd37` is one commit behind `cd0f4374`;
`cd0f4374` adds only `docs/audits/KNOP_FAILURE_ISOLATION_RESULTS.md`, so
`subs/sobew.f` is identical on both commits.

**The earlier "implemented" fix was never committed.** It exists as
uncommitted working-tree changes on top of `cd0f4374`:

```
 M subs/sobew.f                        (the NI=1/2 GAM fix)
 M unit_tests/f90/drv_sob_floor.f90    (exact-shape AQW/HAQW/etc.)
 ?? unit_tests/python/test_sobew_short_ray.py
 ?? audits/                            (audit documents)
```

That is why a fetch of the GitHub commit `cd0f4374` shows the old code:
the remote tree at that SHA is genuinely unfixed; the fix is local
only. This task continues on the same working tree. No commits, resets,
rebases, or clean operations were performed by this task. The
pre-existing knop work (`control_variable_mod.f`,
`rd_control_variables.f`, `fg_j_cmf_nonmono_v1.f`) that was uncommitted
when the fix was written has since been committed externally as
`1c875cd5`, `86fca689`, `247c9ad6`; it does not overlap the sobew/test
files.

## 1. The defect (verified, reproduced)

In each of `subs/sobew.f`, `subs/sobew_grad.f`, `subs/sobew_grad_v2.f`,
`subs/sobew_grad_v3.f` (pre-fix):

- `NI = ND-(LS-NC-1)`, `NI=ND` for `LS<=NC`.
- `GAMH(I)` is initialized for **all** depths before the ray loop
  (`GAMH(I)=CHIL(I)*3.0E-10*R(I)/V(I)/FL`).
- `NI=1` branch: no `ZALONGP`, no `GAM`.
- `NI=2` branch: `CALL ZALONGP` for the continuum (`NORDTAU`) only;
  **no `GAM`**.
- `NI>2` branch: `CALL ZALONGP` then `GAM(I)=GAMH(I)/T1`.
- Common post-branch loop: `T2=EXPONX(GAM(I))` for every `I=1,NI`,
  including the `NI=1` and `NI=2` rays.

Because `GAM` persists across `LS` iterations, the `NI=1`/`NI=2` rays
consumed the `GAM(1:NI)` left by the preceding `NI=3` ray. Since `Z`
(hence `mu=Z/R` and `T1`) depends on `P(LS)`, values from a different
ray are not interchangeable: **stale-ray dependence**, not merely
uninitialized memory. No upstream invariant protects it.

## 2. Reachability per routine

`NI=1`/`NI=2` rays exist iff `NP=ND+NC` (`cmfgen.f:462` allows
`NP=ND+NC` or `ND+NC-2`; the latter stops at `NI=3`). Model
`m15_2p82d_FRESH_sob_dieoff`: `NC=15`, `NP=85=ND(70)+NC(15)` — the
short rays exist.

| routine | call sites | reachable? |
|---|---|---|
| `SOBEW` | `cmfgen_sub.f:2269,2273` (dielectronic/EW), same file last-iteration block; `sub_sob_line.f:145`, `sub_sob_line_v2.f:147`, `sub_sob_line_v3.f:194,198` | YES - model EW diagnostics (`LU_EW/LU_DR/LU_NET/LU_HT`); does not feed rates or the solver |
| `SOBEW_GRAD` (v1) | none in tree (verified by grep) | NO - dead code; fixed for consistency |
| `SOBEW_GRAD_V2` | `obs/cmf_flux_sub_v5.f:2980,2986` (`DO_CMF_EW`/`DO_NONMONO_SOB` ABS_EW + FORCE_MULT), `disp/maingen.f:3956,6486` | YES - obs/disp analysis tools, same NP |
| `SOBEW_GRAD_V3` | none in tree | NO - dead code; fixed for consistency |

## 3. Physical formula preserved

    tau_S(mu) = GAMH / T1,  T1 = |1 + sigma*mu^2|,  mu = Z/R
    sigma < 0:  T1 = max(|1+sigma*mu^2|, -0.2*(1+sigma), 1.0e-4)
    sigma >= 0: T1 = |1 + sigma*mu^2|

Standard LOS-gradient form (`dv_los/ds = (V/R)(1+sigma*mu^2)`;
Lamers, Cerruti-Sola & Perinotto 1987, ApJ 314, 726), identical to the
`SOBJBAR_SIM` prescription. No new floor; `GAMH` (incl. sign)
unchanged.

**Upstream inconsistency found:** the *reference* `sobew_grad.f` uses a
**signed** denominator `GAMH/(1+sigma*mu^2)` (no `ABS`), while ref
`sobew.f`, `sobew_grad_v2.f`, `sobew_grad_v3.f` use `ABS(...)`. The
fork's `sobew_grad.f` already carried the ABS+floor form before this
fix (added earlier with the comment "should not carry the unguarded
form"); this task preserved it. For `sigma < -1` the reference v1
output is corrupted by sign-flipped denominators (NaN / +Inf bursts) —
independent of the stale-ray bug; see the tests document.

## 4. Minimal correction applied (all four routines)

Identical relocation in each file: `CALL ZALONGP(R,Z,P(LS),NI)` and the
complete `T1`/`GAM` block moved to just after `AV(1:NI)=0.0`, before
the `THICK`/`IBOUND` setup and the NI-specific branches (matching
`SOBJBAR_SIM` ordering). The duplicate `ZALONGP` in `NI=2` and the
`ZALONGP`+`GAM` block in `NI>2` were removed. No transfer-solution,
boundary-condition, force, EW, or emissivity changes.

Post-fix anchors:

| file | ZALONGP | GAM loop end | NI=1 | NI=2 | EXPONX(GAM) |
|---|---|---|---|---|---|
| `subs/sobew.f` | 78 | 102 | 121 | 124 | ~189 |
| `subs/sobew_grad.f` | 73 | 94 | 113 | 116 | 181 |
| `subs/sobew_grad_v2.f` | 88 | 108 | 127 | 130 | 233 |
| `subs/sobew_grad_v3.f` | 91 | 111 | 130 | 133 | 236 |

## 5. Reproduction evidence

Instrumented debug copies (test-only, production untouched) dumped
`LS,I,NI,GAM,AV,TB,DTAU` for a `sigma=-1.05` case:

- Reference: `LS=28 (NI=2)` GAM(1:2)=`8.574,4.828` — identical to ray
  27 (`NI=3`), i.e. stale; `LS=29 (NI=1)` GAM(1)=`8.574` still.
- Fixed: `LS=28` GAM=`5.636,3.212`; `LS=29` GAM(1)=`3.750` = `GAMH(1)`
  (tangent at `R(1)`, `mu=0`, `T1=1`). Exact current-ray values.
- Core rays identical except where the pre-existing DAH15 floor
  engages; `AV`, `TB`, `DTAU` identical everywhere — the continuum
  solution is untouched.

## 6. EW / force impact (focused synthetic cases, uniform weights)

New vs reference build for `SOBEW` (ref has neither fix nor floor):

| sigma | chil | EW_new | EW_ref | dEW |
|-------|------|--------|--------|-----|
| +2.0  | +1e9 | 6.834e8 | 6.885e8 | -0.74% (pure stale-GAM fix) |
| -1.05 | +1e9 | 5.359e7 | 4.982e7 | +7.6% (fix + floor) |
| -1.0  | +1e9 | 4.876e7 | 4.531e7 | +7.6% (fix + floor) |
| -1.05 | -1e9 | 2.06e167 | NaN | unsaturated EXPONX gain; both diverge |

`FORCE_MULT` (`GLINE(I)+=TB(I)*EXPONX(GAM(I))*HAQW`) receives the same
corrected `EXPONX(GAM)` on short rays in the grad routines; EW
diagnostics change identically.

## 7. Test-harness correction

`unit_tests/f90/drv_sob_floor.f90` previously declared
`aqw(NDMAX=400,NPMAX=40)`, `haqw(NDMAX,NPMAX)`, `chilm(ndmax,nsmax)`,
`etalm(ndmax,nsmax)`, `bbcor(ndmax,nsmax)`, `znet(ndmax,nsmax)`,
`vb(ndmax,nsmax)`, `vc(ndmax,nsmax)`, `betac(ndmax,nsmax)` and passed
them to routines whose dummies are `(ND=25,*)` / `(ND,NSIM)` shaped —
column-major stride mismatch scrambled every ray weight (`AQW(I,LS)`
read `aqw` at flat offset `I+(LS-1)*25` inside a 400-strided array).
Invisible under uniform weights; broke all weight-sensitive tests.
Fixed: all 2D arrays are now `allocatable`, allocated at exact
`(nd,np)` / `(nd,nsim)` after the header read. A `sobewg1` branch was
added to exercise `SOBEW_GRAD`. `test_sob_tau_floor.py` gained a
`LINK_ONLY` object (`sobew_grad.f`) so the shared driver resolves.

## 8. Results

- `python/test_sobew_short_ray.py`: **32 passed, 0 failed** (8 checks
  per routine: NI=1, NI=2, NI>2-identity, ordering, monotonic sigma,
  sigma=-1, negative GAMH, sigma=-2). See SOBEW_SHORT_RAY_GAM_TESTS.md.
- `python/test_sob_tau_floor.py`: **17 passed, 0 failed** (regression).

## 9. Limitations

- `SOBEW`/`SOBEW_GRAD_V2` are diagnostic paths (EW / force output),
  not the population solver; the stale GAM did not affect rates.
- The production model's reported EWs were computed with stale
  short-ray GAM; real-model magnitude not estimated here.
- `sigma<0` new-vs-ref comparisons conflate the fix with the fork's
  pre-existing floor; the `sigma=+2` checks isolate the fix.
- For `sobewg1` at `sigma<-1` the reference is unusable (signed
  denominator); those checks use an absolute current-ray prediction
  (agreement ~2e-9 relative).
