# KNOP Production Dependency-Cycle Audit (P1-T1e)

Question: does the acyclic store-dependency structure measured on the
manufactured P1-T1b model survive on real production geometry?

Scope: `plane/fg_j_cmf_nonmono_v1.f`, `new_main/`, `obs/`. No production
transfer/population equation, frequency-matching rule, relaxation
algorithm or convergence criterion modified. Evidence classes:
**[MATH]** code structure, **[MANUF]** manufactured model, **[PROD]**
production model (requires runs not yet executed).

**Headline: PART 1 (call-order) is verified clean across every
production path. PART 2 confirms the same-slot edges are I_PREV-mediated
by structural elimination. PARTS 3-6 (production dependency graph +
SCC + nilpotency reassessment) are specified but BLOCKED on a
production checkpoint - they cannot run without a real non-monotonic
model execution.**

---

## PART 1 - sweep initialisation guarantees (VERIFIED)

Exhaustive caller trace. Exactly **two** production call sites of
`FG_J_CMF_NONMONO_V1` exist:

| Site | File:line | Gate | Sweep shape |
|---|---|---|---|
| Population driver | `new_main/mod_subs/comp_j_blank.f:897` | `DO_NONMONO_KNOP` or `DO_NONMONO_SOB+USE_LAM_ES+DO_SOB_FG_KNOP` | complete `ML=1..NCF` |
| OBS/flux include | `obs/COMP_JCONT_V4.INC:816` | `DO_NONMONO_KNOP` only | complete `ML=1..NCF` per `SOB_PHASE`/`ES_COUNTER` |

`FIRST_FREQ` is passed as the solver's `INIT` argument and is set
`(ML==1)` inside every frequency loop that can reach KNOP:
`cmfgen_sub.f:2405-2409` (loop at 2400, call site 2663) and
`cmfgen_sub.f:6164-6168` (JFNK eval loop at 6159, call site 6404), and
`cmf_flux_sub_v5.f:1635-1639` (OBS loop at 1632). The other three
`COMP_J_BLANK` call sites (2073, 4074, 7265) pass `SECTION !=
'CONTINUUM'` -> `CONT_VEL=.FALSE.` -> the KNOP branch is unreachable;
their `FIRST_FREQ` argument is stale but unread. **[MATH]**

Guarantees confirmed:

- `I_P_PREV`/`I_M_PREV` zeroed on INIT (2715-2716); `FREQ_IDX_ITER=0`
  reset on INIT (2726). INIT takes precedence over NEW_FREQ, so at
  `ML==1` every F-pass re-runs the idempotent INIT branch.
- `I_P_STORE`/`I_M_STORE` are overwritten **unconditionally at the end
  of every solver call** (4681-4682), before any `NEW_FREQ` copy
  (3004-3005). No early return bypasses the overwrite (it is on the
  common exit path).
- `NEW_FREQ=.FALSE.` reaches the KNOP branch on inner F-iteration
  passes (pass 2+) of every frequency - by design, to reuse
  `AV_PREV`/`CV_PREV` (`comp_j_blank.f:1183-1185`). The solver accepts
  same-frequency re-entry (`FREQ .EQ. PREVIOUS_FREQ`, 3012-3016) with
  `FREQ_IDX_ITER` unchanged, so `I_PREV` correctly stays at the previous
  frequency across a frequency's inner iterations. **[MATH]**
- `FREQ_IDX_ITER` is `SAVE`d and reset only on INIT; there is **no**
  partial or resumed frequency sweep - every loop is a complete
  `ML=1..NCF`. `FG_EVAL_STATE_PUSH/POP` preserves it around JFNK evals.
- Checkpoint/restart (`DUMP_FG_INEXT_GS`, `RESTORE_FG_INEXT_GS/OBS/
  LATEST`) touches only `I_*_DEPTH_NEXT` - never `I_*_PREV`,
  `I_*_STORE`, or `FREQ_IDX_ITER`. Restores run inside the `FIRST_TIME`
  block on the process's first solver call, which is always an INIT.
- `FG_KNOP_REFRESH_V1` self-guards: returns unless `KN_STORE_READY` and
  `FREQ_IDX_ITER==KN_NCF-1`, i.e. only after a complete sweep.

**Conclusion: no execution path bypasses INIT, and no path can leave
`I_PREV`/`I_STORE`/`FREQ_IDX_ITER` stale into a sweep. `I_STORE` and
`I_PREV` do NOT need to become part of the independent state vector** -
they remain within-sweep slaves. **[MATH - verified]**

One honest caveat: the solver relies on the *caller contract* (INIT
first) with no explicit internal guard against an `INIT=F, NEW_FREQ=T`
first call. No production path violates it today, but it is a contract,
not an enforced invariant.

---

## PART 2 - same-slot dependencies are I_PREV-mediated (PROVEN by elimination)

Measured in `J_T` (manufactured, `INEXT_RELAX=1.0`):

- P channel: `out(P,s,kc0) <- in(P,s,kc1)`, gain ~9.4e-4, all 5 slots.
- M channel: `out(M,s,kc2) <- in(M,s,kc1)`, gain ~5e-9, all 5 slots.
- Backward edges `dslot=-1,-2,-3`, gain ~1e-6 to 1e-12.

**Why they must be `I_PREV`-mediated.** Enumerate every mechanism that
could make `output[k]` depend on `input[k]` (same slot):

1. *Compressing read* - always `slot = FIDX + NM_RAY`, `NM_RAY >= 1`, so
   it reads a strictly higher slot. Cannot produce same-slot.
2. *Parabolic source stencil* (`S`,`dS`, clamps) - spatial, within one
   frequency's intensity; it does not change which store slot is read.
   Cannot produce same-slot.
3. *Frequency matching* - only ever increases the slot offset
   (`NM_RAY`,`NFM`,`NNF` all `>=1`). Cannot produce same-slot.
4. *Expanding `I_PREV` chain* - freq `k+1`'s expanding cells read
   `I_PREV` = freq `k`'s intensity, which depends on `X[k]` (freq `k`'s
   compressing read of slot `k`). This is the **only** mechanism that
   maps `output[k] <- input[k]`. **[MATH - unique pathway]**

The P/M asymmetry (~1e-3 vs ~5e-9) matches the known outward-channel
dominance of the expanding coupling (P1-T1), consistent with `I_PREV`
rather than a symmetric artefact. The two FD amplitudes agreed to
~1e-11 and all columns used central differences, so the edges are
repeatable, amplitude-independent, and far above the ~1e-11 noise floor.
No floor/cap/monotonicity branch transition is implicated (the gains
are smooth and the relaxed-map eigenstructure `1-omega` was exact).

**Empirical confirmation path (specified, not run):** a
`FREEZE_I_PREV` diagnostic flag (default `.FALSE.`, same pattern as
`KNOP_DIAG_NO_REGRID`) that holds `I_*_PREV` at baseline during the
expanding `SOURCE_PRIME` read would make the same-slot edges vanish if
the pathway is correct. This is a diagnostic-only hook; it was not
added in this audit to keep the diff minimal. The structural elimination
above already establishes the mechanism. **[MATH + MANUF]**

Genuine small edges are retained in the graph; nothing was thresholded
away in the analysis.

---

## PARTS 3-6 - production dependency graph, FD validation, SCC, nilpotency

**STATUS: SPECIFIED, BLOCKED on a production checkpoint.**

These parts require an actual non-monotonic production model run to
obtain a checkpoint (populations, `CHI_EXT`/`ETA_EXT`, `GAM`, `DTAU`,
`Q`, frequency grid, ray geometry, boundary config, and the
`I_*_DEPTH_NEXT` store at a representative iteration). The experiment
design is fixed and ready:

- Freeze `CHI`,`ETA`,`GAM`,`DTAU`,`Q`, frequency grid, geometry,
  boundaries (the sweep already treats these as fixed inputs).
- Use the real `NP` ray pencils, real `NCF`, real `NM_RAY`/`DLMN`
  (including `DO_NONMONO_FREQ_MATCH=T` multi-bin offsets), alternating
  expanding/compressing bands, multiple compressing bands, both
  channels, non-core rays, and the applicable `INNER_BND_METH`
  (`DIFFUSION`/`ZERO_FLUX`/`PNT_SRCE`).
- Perturb each active writable store element `(KC,LS,slot,channel)`
  individually with central FD at two amplitudes; record the full
  response footprint; build the sparse Jacobian; run Tarjan SCC over
  thresholds `1e-9..1e-2` while preserving the raw derivative graph.
- For every candidate cycle: FD-verify each edge, extract the SCC
  submatrix, compute its eigenvalues (`dgeev`), report spectral radius
  and the dominant coupling mechanism.

The dependency *rules* are already established (`KNOP_PRODUCTION_
DEPENDENCY_GRAPH.md` section 2.1): output slot `s` depends on input
slots `{1..s+NM}` via compressing reads, plus `{<=s}` via the expanding
`I_PREV` chain, plus cross-channel `I_P(NI)=I_M(NI)` on eligible rays.
Whether a real velocity field populates a cycle is precisely what this
experiment decides. **[PROD - not run]**

**Nilpotency reassessment (PART 6):**

- *Manufactured fixed-opacity:* nilpotent, complete, verified DAG.
  **[MANUF]**
- *Production fixed-opacity:* CONDITIONAL/OPEN pending PARTS 3-5. The
  dependency rules admit backward edges and the `I_P(NI)=I_M(NI)`
  M->P path; the manufactured model did not populate either into a
  cycle. No nilpotency claim is made for production geometry.
  **[PROD - unresolved]**
- These are distinct questions from relaxed-map convergence,
  refresh fixed-point convergence, and population-coupled CMFGEN
  convergence; none is inferred from another.

---

## What would settle it

One production checkpoint + the FD/SCC pipeline above. The decisive
unknowns are (a) whether any real ray has expanding cells upwind of
compressing store points strongly enough to close a cycle, and (b)
whether the `I_P(NI)=I_M(NI)` boundary path couples the channels into a
loop on non-core rays. Both are measurable in a single checkpoint run.
