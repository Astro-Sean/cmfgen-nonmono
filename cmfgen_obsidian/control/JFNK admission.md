# JFNK admission: the three-phase strategy

JFNK replaces the *middle* map's linear solve, not the iteration
structure. It is deliberately an **endgame accelerator**, gated behind
sustained criteria - because each coupled eval costs ~250-400 s and
Newton fired outside its basin is pure waste (the it11 failure, see
[[JFNK solver]]).

```text
Phase 0  Lambda/Knop iterations until the field is mature AND the
         demanded correction is small.
Phase 1  Frozen-J JFNK  - cheap evals, BA ~ exact Jacobian; a strict
         upgrade over one linearized SOLVEBA step per iteration.
Phase 2  Coupled JFNK   - only when |F| is small AND the field is a
         demonstrated fixed point; the true endgame step.
```

![[fig_three_phases.png]]

> [!quote] Knoll & Keyes 2004, JCP 193, 357 ([DOI](https://doi.org/10.1016/j.jcp.2003.08.010))
> JFNK is "most effective near the solution"; Krylov methods "are not
> effective solvers in their own right" without preconditioning.
> That is the whole reason this file exists: the gates are there to
> keep Newton inside its basin and keep a preconditioner on it.

Analogy: three gears on the same engine. Lambda iterations are low
gear - big torque, slow crawl, can't stall. The frozen JFNK tier is
middle gear - fast but only grips once the engine (the field) is
warmed up. Coupled JFNK is top gear - the whole machine in one
calculation - only safe at full speed, near the root. You would never
start a cold engine in top gear; you would never fire Newton on a
field still settling.

This is the literature-standard operator-split -> Newton handoff made
quantitative (Hauschildt & Baron practice; Arramy et al. 2024,
[arXiv:2406.17234](https://arxiv.org/abs/2406.17234), diagnose their
own stalls the same way). Arramy's measured JFNK advantage is "up to
a factor 2" in formal-solution evaluations - which is exactly why
Newton sits behind two sustained gates instead of replacing the loop.

## Admission chain (in order)

```text
JFNK_THIS = USE_JFNK_SOLVER          (VADAT DO_JFNK; USE_JFNK alias)
        AND ITS_SINCE_START > JFNK_WARMUP        (process warmup)
        AND (.NOT. LAMBDA_ITERATION .OR. ANYITER) (live SE call)
        AND (.NOT. DO_NONMONO_KNOP .OR. .NOT. GS_FREEZE)
        AND USE_LAM_ES
```

`DO_JFNK` is admitted under `DO_NONMONO_KNOP` or `DO_NONMONO_SOB`
(read-time STOP otherwise). Under SOB the GS_FREEZE clause is vacuous
(there is no lagged Knop field) and the coupled tier is unreachable -
`JFNK_COUPLED=T` is rejected at VADAT read, so only the frozen-J
residual can fire. `USE_LAM_ES=T` is already mandatory for SOB.

- `GS_FREEZE` - **not** the raw `GS_FG_IMMATURE` - is the field test:
  it subsumes the instantaneous flag, the `N_MATURE/GS_REL_STABLE`
  stability hold, and the `MOM_TRIAL` veto. Newton fires only when the
  lagged field has been a demonstrated fixed point for
  `GS_REL_STABLE` consecutive sweeps and the SE update is live.
- No JFNK during `MOM_TRIAL` or other bounded rollback-protected trial
  iterations.
- Admission also requires `USE_LAM_ES=T`. With `USE_LAM_ES=F`,
  `MOM_J_CMF_NONMONO_V1` updates lagged `J_NEXT_FRQ`/`H_NEXT_FRQ`
  stores in place, with no eval freeze. JFNK therefore stays deferred
  rather than allowing trial sweeps to mutate the stored field.

## Demand tiers

```text
frozen tier : BA_LAST_RAW < JFNK_FRZ_GATE  for JFNK_GATE_NREQ
              consecutive eligible calls ->  JFNK_THIS, coupled=F
coupled tier: BA_LAST_RAW < JFNK_COR_GATE  for NREQ consecutive
              eligible calls             ->  JFNK_COUPLED_THIS=T
```

- `BA_LAST_RAW` is the last **non-LAMBDA** call's demanded
  correction - the residual scale of the SE problem measured on the
  operator Newton would actually see (04-Oct-2026: was `MAXCH`, the
  previous call's demand; in the alternating regime that call is a
  LAMBDA solve whose small demand reflects the LAMBDA attractor.
  The opt run admitted JFNK at it73 on a 95% LAMBDA demand while the
  last full demand was 4.8e4% and |F|_w=2.2e12 vs target 2.2e8).
  Before the first non-LAMBDA call (`BA_LAST_RAW=0`) it falls back
  to `MAXCH`.
- Hit counters (`JFNK_FRZ_HITS`/`JFNK_CPL_HITS`) advance **only on
  otherwise-eligible calls**: during warmup, Lambda, freeze, or
  `USE_LAM_ES=F` the demand is measured against a field that is not an
  admissible JFNK operator (and under freeze the deferred correction
  keeps growing), so those calls must not pre-arm the gate. Counters
  reset on ineligible calls - one transiently low MAXCH cannot trigger
  Newton.
- Demand **contraction** is considered before admission: demand must
  be falling, not just below the gate.
- A BA demand whose residual test is pending, rejected, masked, or
  soft-rejected cannot add admission hits. JFNK waits for a later
  residual-accepted BA measurement; this prevents a small failed probe
  from looking like a contracting Newton basin.

## Hybrid config values (`dessart_m15_2p82d_JFNK_hyb.yaml`)

| Key | Value | Role |
|---|---|---|
| `JFNK_FRZ_GATE` | 1e3 | frozen-tier demand gate (%) |
| `JFNK_COR_GATE` | 1e2 | coupled-tier demand gate (%) |
| `JFNK_GATE_NREQ` | 2 | consecutive sub-gate calls |
| `JFNK_WARMUP` | 6 | calls before JFNK may fire |
| `JFNK_MAX_EVALS` | 40 | eval budget per activation (SOB zlaw14: 200 - 40 starves one GMRES cycle) |
| `JFNK_EPS` | 1e-3 | noise-aware FD floor (coupled) |
| `JFNK_M` | 30 | GMRES restart dim |
| `GS_REL_STABLE` | 2 | stability hold on the field |

## Intended run signature

```text
gs_warmup -> gs_freeze_forced / Lambda steps -> live BA steps
         -> JFNK gate defers (hits growing)
         -> jfnk_frz admissions
         -> jfnk_cpl near the root
```

A run showing `jfnk_cpl` while `frozen=T` or
`nmature < GS_REL_STABLE` in the same call is a bug.

## Cost model (why the tiers exist)

```text
Lambda iteration:      ~1 sweep + SE solves ~ 400-1000 s
Frozen-J JFNK eval:    ~1-5 s   (SE assembly only, no sweep)
Coupled JFNK eval:     ~250-400 s (full Knop sweep re-run)

Frozen activation  ~= (1 + nits*(1+m) + LS) * ~3 s   ~ 2-5 min
Coupled activation ~= (1 + nits*(1+m) + LS + retries) * 300 s
                      ~ 30-50 evals ~ 3-5 h; worst case ~100+ ~ 9 h
```

JFNK pays for itself only if one activation replaces >= ~30 Lambda
iterations - viable only inside the basin with a small Newton count.
Hence: cheap frozen tier first (validates the machinery, cleans the
outer map at ~zero cost), coupled only at the end.

## Design invariants (for future edits)

1. **Newton only on a live SE update** - keep `GS_FREEZE` in the
   admission test; weakening it to `GS_FG_IMMATURE` alone re-opens the
   trial-iteration hole.
2. **Ordering before accuracy** - a converged field spent on a safe
   step beats an accurate Jacobian spent on a transient one.
3. **Hysteresis is load-bearing** - every counter exists because a
   single lucky call is not evidence; new gates need consecutive-call
   semantics and reset-on-ineligible.
4. **Every unchecked step needs a net** - resacc for BA, best-iterate
   + commit checks for JFNK, caps for Lambda.
5. **The two residuals must agree at baseline** - the stored-term
   split must preserve `F(0) = outer STEQ` or the F(0) test (correctly)
   rejects.
6. **Log-space consistency** - every consumer of x uses the same
   `exp(x)-1` map and +/-40 clip.
7. **Coupled evals are the budget** - anything spending formal sweeps
   belongs behind a flag or counter (`JFNK_DEEPDIAG`, `JFNK_MAX_EVALS`).

## Status

- Frozen-J path: functional, exercised as SOLVEBA replacement;
  first production validation on the hybrid run (resubmitted
  03-Oct as job 12073582 on the order-of-ops fixed binary,
  `7221de9`).
- Coupled JFNK: correctly deferred by the demand gate; has never yet
  produced an accepted Newton step on this model - the gate defers
  every full-designated call while raw demand stays saturated.
- Lambda fallback after JFNK failure: validated (seq run kept
  contracting after `JFNK_DONE`); the `JFNK_DISABLED` self-disable
  itself was fixed 03-Oct (`e8abf97`) - it previously printed the
  fallback message but ran Newton anyway.

## Related

- [[JFNK solver]]
- [[Field maturity and GS freeze]]
- [[Failure modes and guards]]
- [[JFNK_SOLVER_AUDIT]] - cost asymmetry + JVP noise floor
