# CMF_FLUX continuum J/f iteration and J_RESCUE audit (Phase 2.1)

Audit of `obs/COMP_JCONT_V4.INC` as included by `obs/cmf_flux_sub_v5.f`.
Benchmark: `m15_2p82d_SOB_l14_lte_adp2_tauto`
(`DO_NONMONO_SOB=T, USE_LAM_ES=T, COHERENT_ES=T, COMP_F=T`),
159033 spectrum points, repeated renders median |rel| diff ~1.5e-3,
p99 ~4e-2, ~223 printed J_RESCUE events per render.

## 1. Call chain and branch selection

`CMF_FLUX_SUB_V5` loops `DO 10000 ML=1,NCF` (sub_v5:1748) with
`SECTION='CONTINUUM'` and includes `COMP_JCONT_V4.INC` under
`IF(COMPUTE_J)` (sub_v5:2111). The same .INC is included a second
time under `SECTION='LINE'` (sub_v5:2991), where the resulting `RJ`
feeds the Sobolev line source `SOURCE=ZETA+RJ*THETA` (sub_v5:2995).

Inside the .INC, three top branches exist (`CONT_VEL & USE_FIXED_J`;
`.NOT. CONT_VEL & THIS_FREQ_EXT`; `EDDINGTON`; else Rybicki). The
production render takes `CONT_VEL` without `USE_FIXED_J`, entering
the J/f block at line 670.

Formal-solve dispatch (line 827+): `PLANE_PARALLEL*` / `USE_FORMAL_REL`
/`USE_DJDT_RTE`/`USE_J_REL` are off; under `DO_NONMONO_KNOP` it calls
`FG_J_CMF_NONMONO_V1`; **under plain `DO_NONMONO_SOB` it calls
`FG_J_CMF_V12`** (line 893). The moment dispatch then takes
`USE_LAM_ES=T` at line 965: `RJ := TC` (the formal-solve J) with
`RSQHNU` rebuilt - no moment solve runs at all.

## 2. The mathematical iteration

With `USE_LAM_ES=T` the loop at lines 809-1064 is:

```
repeat (L=1,2,...):
  TA  = ETA_CLUMP + CHI_SCAT_CLUMP * RJ          ! source, coherent ES
  FG_J_CMF_V12(TA,...) -> TC, FEDD               ! formal solve on
                                               !   scattering source
  RJ  = TC                                       ! J update = formal J
  [DAMP_FJ mixing - inactive at default 1.0]
  if J_RESCUE and MAXVAL|RJ| > RESC_CAP_IT
      or non-finite: F_DID_NOT_CONV=T, EXIT
  test |FEDD-FOLD|_inf <= ACC_EDD_FAC -> converge
  FOLD = FEDD
  exit at L>20: F_DID_NOT_CONV=T
until converged or exit
```

This is a **variable-Eddington-factor / Lambda-type iteration on the
scattering term** (Hubeny & Mihalas 2014, ch. 13): J from the
previous pass feeds the coherent scattering emissivity; the formal
solve returns J and f; f-convergence proxies J-convergence. It is
NOT an ALI - there is no accelerated operator, no stored diagonal,
no Krylov space. `COMP_F=T` makes the f test run every pass.

## 3. The two rescues (distinct mechanisms)

### 3a. RJ_ES substitution - line 789 (NES_SAN)

Active only under `COHERENT_ES=F`: if `MAXVAL|RJ_ES|` is non-finite
or exceeds `RESC_CAP`, `RJ_ES := RJ` and `NES_SAN++`. This replaces
the *incoherent* electron-scattering J (read from LU_ES, convolved
spectrum) with the local coherent J, changing the physical ES
emissivity `ESEC*RJ_ES` in TA/SOURCE at that frequency. **Zero
events in this model (COHERENT_ES=T).**

### 3b. Previous-frequency restore - line 1074 (NRESC)

Post-loop: if `MAXVAL|RJ|` or `MAXVAL|FEDD|` is non-finite, or
`MAXVAL|RJ|>RESC_CAP`, or `MAXVAL|FEDD|>2`, restore seven fields
from the previous frequency (`RJ,FEDD,GEDD,N_ON_J,RSQHNU` from
`*_PREV`, `HBC_CMF,NBC_CMF` from `HBC/NBC_PREV`), `NRESC++`.

All 892-odd logged events on this model are 3b, all in
`section CONTINUUM`, at consecutive grid frequencies. The restored
vector is then:

- written to `LU_EDD` (line 1137) - the next frequency's `JNU_PREV`
  and the next *render's* warm start (`RD_CONT_J_OBS`, line 142/147);
- used in `K_MOM=RJ*FEDD` and `SOURCE` for the line path;
- consistent J/f by construction (prev-frequency pair) but at the
  wrong frequency.

## 4. Thresholds

- `RESC_CAP = MIN(RESC_JMAX, 1e4*MAXVAL|JNU_PREV|)`;
  `RESC_JMAX=1e-3` default, `J_RESCUE=T`, `DAMP_FJ=1.0` (off).
- `RESC_CAP_IT = MIN(1e10, MAX(1e4*MAXVAL|JNU_PREV|,1e-3))` - looser
  in-loop early exit.
- FIRST_FREQ: `RESC_CAP=RESC_CAP_IT=RESC_JMAX`, no restore possible.

Physical context (comment lines 720-733): converged J peaks ~5e-5,
adjacent-frequency Jmax varies <50x, so 1e4x-prev rarely binds real
features and 1e-3 is a safe ceiling. **Residual order dependence**:
the cap is anchored to the previous converged field, so a rescued
frequency changes the cap seen by ML+1 - a modest ratchet risk
remains despite the comment's reasoning (the restore writes the
SAME RJ as JNU_PREV, so 1e4*prev is unchanged after a restore;
ratchet only via legit converged growth).

## 5. First-frequency defect (confirmed)

At `ML=1`, `J_RESCUE .AND. .NOT. FIRST_FREQ` skips the restore
entirely; the only backstop is `ABS(RJ(1))>1E30` at line 1330:

- checks **element 1 only** - an interior-depth explosion or NaN at
  `I>1` passes the guard;
- `NaN > 1E30` is `.FALSE.` - even a boundary NaN slips through;
- a `RESC_CAP_IT` in-loop exit leaves the last (diverged) `RJ`,
  which is then written to `LU_EDD` at line 1137.

So a first-frequency interior divergence writes an invalid field
to `LU_EDD` and continues. Manufactured-test requirement: `RJ(1)`
finite, `RJ(k)` = NaN or >1e30 at ML=1 -> today the field is
written silently. Correct response per the task: fail closed (no
physical fallback exists for the first frequency).

## 6. IEEE finiteness audit (defect)

Every validity check is `T=MAXVAL(ABS(x)); T.NE.T`:

- `MAXVAL` is a reduction; under ifort -O3/xHost the compiler is
  free to emit MAXPD-lane reductions that **do not propagate NaN**
  on all lanes/orders - a NaN can compare out of the result. The
  `.NE.` self-test is therefore compiler- and lane-order-dependent,
  not guaranteed.
- A NaN at a non-max position can be masked by a larger finite
  value even when the reduction does propagate.

Recommendation: elementwise `IEEE_IS_FINITE` (or
`ANY(.NOT. IEEE_IS_FINITE(RJ))`) for `RJ`, `FEDD`, `RJ_ES`. Not a
floor - a validity test.

## 7. Non-convergence vs invalid

`F_DID_NOT_CONV` is set on runaway exit and L>20 exit but is
**never read anywhere** - write-only diagnostic. Post-loop rescue
tests only field validity (`|RJ|` cap, `|FEDD|<=2`, finiteness):
a finite-but-unconverged exit (the documented common case, residual
~1e-3) is stored as converged - intended stock semantics.

## 8. Frequency-order and reproducibility

The ML loop is **strictly serial** (no OpenMP on it). Persistent
coupling between frequencies: `JNU/FEDD/GEDD/N_ON_J/RSQHNU/HBC/
NBC_PREV` (in-memory), `LU_EDD` records (across frequencies AND
across renders via `RD_CONT_J_OBS`), `FG_COUNT`, `NEW_FREQ->AV/CV_
PREV` inside the solver, and `RESC_CAP` derived from `JNU_PREV`.

Why renders differ despite identical inputs: the f-iteration sits
near the `RESC_CAP`/`f<=2` boundaries on this model; the formal
solve's internal parallel reductions make iterate k-threshold
crossings sensitive at FP level. A single event flip changes the
stored `LU_EDD` field -> different warm start at ML+1 -> different
downstream iterate -> different event set. The rescue converts
microscopic FP noise into macroscopic (whole-RJ-vector) differences.
**The ~223 printed events/render are not random - they are V12
diverging in the compression regions** (`sigma < -1`, the regime
the comment at line 877 acknowledges V12 mishandles: "aborts or
returns expansion-biased intensities"), and the subset that crosses
threshold varies run to run.

Physical significance: each restore substitutes a *different
frequency's* J and f. For the continuum flux that is a patched
spectrum, not converged; for LINE-section frequencies the same
machinery feeds the Sobolev source directly (no LINE rescues
logged on this model).

## 9. Proposed diagnostics (disabled by default)

`DO_J_RES_DIAG` (new reader key, default F): per event, append to
`J_RESCUE_EVENTS`: type (NES_SAN/NRESC), `FREQ_INDX`, FL, SECTION,
argmax depth, `max|RJ|`, `max|FEDD|`, `max|RJ_ES|` (when read),
`RESC_CAP`, `RESC_CAP_IT`, L at exit, `F_DID_NOT_CONV`, trigger
class (nan/inf/j-cap/f-phys/nanf), `max|JNU_PREV|`, `LU_EDD`
written flag, and at render end the total NRESC/NES_SAN per
section (currently only per-event prints exist; no totals are
emitted anywhere). Existing first-divergence dump stays behind
the same flag.

## 10. Proposed executable tests (production objects)

A. Converged finite iterate -> clean write, NRESC=0.
B. Finite slow-converge (L>20) -> stored, not rescued (stock
   semantics preserved).
C/D. Diverged J at outer/interior depth, ML>1 -> restore fired.
E. NaN at interior depth -> restore fired **today may fail**:
   MAXVAL masking is the test target.
F. +Inf -> restore. G. `|FEDD|>2` finite J -> restore.
H. Invalid RJ_ES under COHERENT_ES=F -> NES_SAN path only.
I. First-frequency interior divergence -> currently writes field;
   expected post-fix behaviour: abort.
J/K. Restore at ML=k -> JNU_PREV at k+1 carries restored vector
   (order-dependence regression test).

## 11. Literature fit

Variable-Eddington-factor iteration per Hubeny & Mihalas ch.13;
Rybicki & Hummer ALI machinery is absent - do not cite it here.
The Knop (2009) solver exists precisely for `sigma<-1`; the SOB
render still dispatches the monotonic-regime V12 on this model's
compression zones, which is the likely diverge driver (to be
quantified by the proposed per-event diagnostics before any solver
change is argued).

## 12. Identified defects (ranked)

1. First-frequency interior divergence/NaN written to LU_EDD
   unchecked (sec 5) - correctness defect.
2. NaN detection via `MAXVAL(..).NE.itself` is
   compiler/order-dependent; `NaN>1E30` is also `.FALSE.` at the
   final guard (sec 6) - correctness defect.
3. `F_DID_NOT_CONV` write-only (sec 7) - observability defect.
4. `RJ_ES` substitution silently changes the ES emissivity physics
   at rescued frequencies (sec 3a) - physically significant, needs
   the event diagnostics (not used on this model).
5. Restore feeds forward through `JNU_PREV`/`LU_EDD` and across
   renders (sec 8) - design behaviour now proven to be the
   nondeterminism amplifier, not FP noise.
6. `RESC_CAP_IT` in-loop exit leaves the diverged iterate as the
   stored field for FIRST_FREQ (sec 5) - same defect as 1.

## 13. Implementation status (Phase 2.1A + 2.2)

### 2.1A validity guards - committed

- `obs/jres_validity.f`: `JRES_TRIG` (elementwise IEEE_IS_FINITE +
  cap, returns trigger class + offending index + max) and
  `JRES_FIN`.  Replaces the `MAXVAL(ABS()).NE.itself` checks at all
  three validity sites (RJ_ES sanitize, in-loop RJ runaway exit,
  post-loop RJ+FEDD).  Defect 2 closed.
- `COMP_JCONT_V4.INC` FIRST_FREQ path: an invalid RJ or FEDD now
  aborts with field/depth/FL/FREQ_INDX/section/trigger BEFORE
  SOURCE construction and the LU_EDD write.  No physical fallback
  exists at ML=1, so fail-closed is the documented response.
  Defects 1 and 6 closed.
- Later frequencies: `JNU_PREV`/`FEDD_PREV` snapshots are validated
  with `JRES_TRIG` before the restore; an invalid fallback aborts
  instead of ratcheting a poisoned field forward.
- Finite-but-unconverged fields (both triggers 0) are untouched:
  caps, the f<=2 bound, solver dispatch, and the 20-pass cap all
  unchanged.
- Manufactured tests: `unit_tests/f90/drv_jres_trig.f` +
  `python/test_jres_validity.py` - 7/7 (converged, interior/boundary
  NaN, interior Inf, over-cap, f>2, finite-unconverged keep).
- The fatal-path test is covered by the code path itself: the STOP
  precedes both SOURCE construction and `WRITE(LU_EDD)`, and no
  record index is allocated before the write, so no LU_EDD record
  can be written on a fatal case.

### 2.2 J_RESCUE_EVENTS - committed (diagnostic only)

- `DO_J_RES_DIAG` (VADAT/CMF_FLUX controls, default F) gates all
  output; `J_RESCUE_EVENTS` is written per render
  (STATUS='UNKNOWN' - OUT_FLUX accumulates renders, events must
  not).
- `JEV` records: one per NES_SAN / NRESC / FATAL with
  type, FREQ_INDX, ML, FL, SECTION, FIRST_FREQ, the actual solver
  dispatch (name strings are assigned inside each CALL branch, not
  inferred from flags), L, F_DID_NOT_CONV, RJ and FEDD trigger
  classes, offending depth, RJ/FEDD/RJ_ES maxima, RESC_CAP,
  RESC_CAP_IT, snapshot validity, restore deltas
  max|field_before-field_after|, and R/V/SIGMA/SIGMA_TRUE at the
  offending depth with a SIGMA_TRUE < -1 flag.
- `JTOT` per-section totals on every section transition (and at
  render end); `JEND` cumulative per-section totals once per
  render.  Counters live in JRES_DIAG_STATE and are independent of
  the NRESC<=20/mod-500 print throttle.
- Event counter scope/aggregation (the ~223 vs ~892 question):
  NRESC and NES_SAN are process-lifetime counters reset only at
  FIRST_FREQ; ~223 printed J_RESCUE lines is the per-render count
  under the print throttle, ~892 lines is the accumulation across
  ~4 appended renders in OUT_FLUX (STATUS='APPEND').  JEV counts
  every event with no throttle, so the per-render JEV count is the
  true event count; JTOT/JEND give per-section splits.
- FATAL records (luedd=F) are emitted but not counted in the
  rescue totals - they are aborts, not rescues.
- Diagnostics change nothing numerical: every record write sits
  behind `IF(DO_J_RES_DIAG)` and reads state only; the solver-name
  and delta locals are write-only under the flag off.  The one
  unconditional overhead is a handful of 24-char assignments at the
  dispatch sites.
- Manufactured test: `unit_tests/f90/drv_jres_events.f` +
  `python/test_jres_events.py` - record layout, section-transition
  JTOT flush, and JEND counters verified.
- Measured on the model (job 13335395, frozen exe
  `exe_frozen/844b17b0fb430c6e_79b0651256a8b508`, cmf_flux sha256
  `79b0651256a8b508...`):
  - OBSFLUX diag-off vs diag-on is byte-identical - the
    instrumentation is numerically inert on this render.
  - `JTOT CONTINUUM nresc=75738 nes_san=0`: 75738 NRESC events per
    render, not ~223.  The audit's ~223 was the print-throttle count
    (NRESC<=20 or mod500; 190 J_RESCUE lines appear in this
    render's OUT_FLUX); ~892 in the earlier log was ~4 renders
    accumulated under APPEND.  JEV/JTOT are the true counts.
  - Trigger class for every event: trig=3 - finite J marginally over
    the absolute cap (rjmax ~1.001-1.05 x RESC_JMAX=1e-3).  No NaN,
    no Inf, no unphysical FEDD (ftrig=0), no first-frequency case,
    no LINE-section event, no FATAL.
  - Every event's offending depth sits on the sigma=-0.999 clamp
    floor (SIG field = -9.9900E-01): 87% at the outermost depth
    (idx 70), the rest at idx 65-69 - the compressed outer region.
    Formal solver recorded at the call site: FG_J_CMF_V12 with
    LAM_ES_TC on all events (the stock OBS dispatch under
    DO_NONMONO_SOB=F in CMF_FLUX_PARAM - SIGMA_TRUE is not even
    allocated on this path, which is why the SIGT field reads
    0.0 here; the -0.999 SIGMA floor is the compression marker).
  - Interpretation narrowed: the rescue does not fire on solver
    explosions here - it fires on a converged Fedd iteration whose
    outermost-depth J lands ~1-5% above the absolute RESC_JMAX
    ceiling at clamped-compression cells.  Whether that is a real
    J excess or a too-tight absolute cap at the outer depth is the
    threshold question for the next phase; the events are not
    evidence of NaN/Inf failure.
- Abort survival: section totals are sealed idempotently and the
  CONTINUUM JTOT is flushed before the phase-1 OBS_FRAME_SUB_V9
  call, so the record survives the model's normal termination
  (GLOBAL_CMF_GRID_UNDERCOVERAGE STOP at the second obs-frame
  pass - same exit as production, which also reaches no RETURN).
- Open: the ten-render ensemble can start now - production build
  succeeds, the executable is frozen and hash-recorded above.
