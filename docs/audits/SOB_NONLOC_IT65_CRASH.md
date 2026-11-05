# SOB non-local model crash audit — m15_2p82d_VAL_sob_nonloc

Companion to `OUTGEN_ITER57_FAILURE_AUDIT.md` (KNOP it57). This audit
covers the SIGSEGV in the m15 2.82-day `DO_NONMONO_SOB + DO_NONLOC_SOB`
validation model (production job 12689657, exit 139, ~1 h in, last flush
inside it65 setup after `T 15 F`).

## 1. Evidence chain

1. Production `OUTGEN` ends mid it65 setup; stderr shows a bare
   `Segmentation fault` from `batch.sh`, no traceback (release build).
2. Debug rebuild (`exe_dbg/cmfgen_dev.exe`, bounds checking on) rerun on
   the interactive partition (job 12723926, `batch_dbg.slurm`) failed in
   ~50 s with:
   - `GFF_VEC bad T: I= 11 T= NaN RNU= 2.000000E+02 ZION= ...`
   - `forrtl: severe (408): Subscript #2 of the array A has value
     -2147483648 which is less than the lower bound of 1`
3. Fault site `subs/gff.f:95-109`: `T(I)=NaN` -> `LOGT=NaN` -> `U`,`GB`
   NaN -> `IU`,`IGB = INT(NaN) = INT_MIN` -> `A(IU,IGB)` wild index into
   the 12x11 free-free Gaunt table. Release build reads out of bounds
   silently -> SIGSEGV.
4. The NaN T arrived via `SCRTEMP`: the debug rerun read the it64
   checkpoint and hit the same fault inside the it65 STEQ segment, so
   `POPS(NT,11)` was already poisoned when it64 wrote its state.
5. it64 upstream state: population solve exploding in the outer
   atmosphere (max increase `5.9E14%` at depth 6, exempted demand
   `2.1E20%` at depth 8) while `SOB_RH` reported `max dJ/J` to `6.8E7`
   and `max dRATE` to `2.3E8` over ~500 K corrections per call. The NaN
   T at depth 11 sits inside that diverging region.

## 2. Classification

- Crash mechanism: **implementation defect** — an unguarded table index
  on a non-finite T. `GFF_VEC` has a `bad T` diagnostic print but the
  index computation proceeds regardless.
- Upstream origin (NaN T production): **convergence/conditioning
  failure** — the non-local Sobolev population corrections were already
  diverging at depths 6-11 during it64; the bad T is a symptom, not the
  defect. Root-cause work on `SOB_RH` magnitudes remains open.
- Same crash class previously seen on the KNOP path; the DBG-T entry
  guard and `KN_REPAIR_T` were added for it but gated on
  `DO_NONMONO_KNOP`, so this SOB model received neither diagnostic nor
  repair.

## 3. Changes made

- `new_main/cmfgen_sub.f` — DBG-T entry guard gate widened from
  `DO_NONMONO_KNOP` to `DO_NONMONO_KNOP .OR. DO_NONMONO_SOB`. Diagnostic
  is unchanged; `KN_REPAIR_T` opt-in now also reaches SOB runs.
- `new_main/mod_subs/solve_for_pops.f` — the `SOL(NT,:)/POPS(NT,:)`
  bad-T localisation diagnostic gate widened the same way (still behind
  `CONV_DIAG_LVL .GE. 2`).
- No repair or clamp added to `GFF_VEC`: the guard converts the silent
  segfault into a diagnosed stop; the correct fix remains upstream
  (whatever writes NaN into `POPS(NT,I)`).

## 4. Verification

- Interactive job 12724056 after the gate change: 100/100 operator,
  8/8 sweep table, 77/77 VADAT, 13/13 gate, 10/10 ray, FP-trap
  `RESULT_OK`.
- Debug rerun with the extended guards (pending): expected to print
  `DBG-T` at it65 entry identifying depth 11 and, if the run is repeated
  from an earlier `SCRTEMP`, to localise the it64 routine that first
  writes the NaN.

## 5. First invalid state — it64's solve returned an all-NaN depth column

`SCRATCH/STEQ_VALS` preserves the per-iteration `STEQ SOLUTION ARRAY`
dumps (WR2D_V2 format `I(J)#` = `STEQ(I, J..J+9)`; `(11)` page covers
depths 11-20). Findings:

- The last dump (it64's final applied LAMBDA solve) contains exactly
  1839 NaN values, all in the first column of the `(11)` page:
  `STEQ(I,11) = NaN` for every one of the 1839 unknowns. The whole
  depth-11 column came back NaN; depths 12-20 are finite.
- All 63 earlier dumps are NaN-free — the corruption first appears in
  it64's final solve.
- Correction application (`locsolut.f`, `solveba_v*.f`) computes
  `POPS(J,I)*(1-STEQ(J,I)*SCALE)`. Every scale-limit test uses `.GT.`/
  `.LT.` comparisons which a NaN operand fails, so a NaN `STEQ` element
  is applied unclamped: `POPS(:,11)` (all levels plus `POPS(NT,11)`=T)
  went NaN, then `SCRTEMP` carried the poisoned depth to it65.
- The same file's `&` "Radiative Equlibrium Equation" dump shows the
  depth-11 RE residual already at `-9.7E5`, growing to `-5.4E7` by
  depth 20 — the depth-11 matrix block was being fed a wildly
  out-of-balance state (consistent with `SOB_RH` `dJ/J` 6.8E7 /
  `dRATE` 2.3E8 corrections).

## 6. Producer analysis -- CORRECTED

**The DRC-injection hypothesis is disproved.** STEQ_VALS holds the
it64 residual dumps written just before the final solve: the full
`STEQ ARRAY` and the radiative-equilibrium row contain no NaN or
Inf anywhere (depth 11 is a finite -2.57E5). `DRC` enters only the
residual via `T2 = ETAL*ZNET - DRC` (sub_sob_line_v3.f:324-328) and
is deliberately kept out of the BA matrix, so the residual was
finite when the solve ran. The all-NaN `SOL(:,11)` column was
produced inside the BA matrix assembly or the CMF_BLKBAND_V3
block solve -- most consistent with a singular or NaN-poisoned
depth-11 block (a depth-local defect gives exactly an all-NaN
column at that depth and nothing elsewhere). The earlier claims
below are retained for the record but the mechanism link is
unproven:

- `BB_COR = T1/T2` (cmfgen_sub.f:4078-4080) cannot go denormal:
  `FL` is the mean of `FL_SIM` over the line group (cmfgen_sub.f
  ~3901) and lines in a group sit close in frequency, so the ratio
  stays O(1). The denormal-BB_COR path is not available here.
- `dRATE` is an absolute rate in emissivity units; comparing its
  magnitude to the dimensionless `dJ/J` is meaningless, so the
  "blow-up is in /BB_COR" inference was unsound.
- `DRC` still has no finiteness check, and `DJ .LE. 0` cannot catch
  a NaN DJ -- those holes remain real, but there is no evidence
  they fired in this crash.

What stands: the non-finite `DRATE` stop and the dominance census
are diagnostic tripwires; the post-solve `STEQ` finite check in
SOLVEBA_V13 stops at the first invalid state with depth/level
location; the restart from the it63 SCRTEMP record reproduces the
it64 solve cleanly for localisation.

Driver sequencing (unchanged): `dRATE` printed 1E7-2.6E8 across the
WHOLE run -- steady-state magnitude, not degradation.

Coupling topology: 4,830 pairs over 70 depths (~69 partners per
depth -- near all-to-all; RH_NDROP=0, no V<=0 depths). The mu^2<=1
window admits essentially every pair on this velocity field. Each
pair's contribution is bounded (the trapezoid weight cancels the
UPS edge), but the acknowledged missing third-resonance extinction
makes every pair a bounded overestimate -- compounding ~69 per depth
into systematically inflated remote pumping. Outer depths can
legitimately be remote-pumping dominated, but a systematic 10-20%
overestimate x 69 partners is enough to push DRC past the local
rate, escalate populations in the SE solve, and set up the observed
outer-atmosphere explosion (depths 6-11) that preceded the NaN
tip-over.

To quantify the systematic-vs-edge question a dominance census was
added: RH_NDOM counts (I,L) cells where |DRC|>|ETAL*ZNET| and
reports alongside the existing max dJ/J / dRATE line. A large census
means the lagged remote term is driving the residual rather than
perturbing it -- the artefact signature.

## 7. Open questions

- Confirm which hole fired first (denormal `BB_COR` overflow vs
  `Inf*0` in the `DJ` sum) — the `DRATE` finiteness stop now prints
  I, L, DJ, CHIL, BB_COR, NPART on reproduction (resubmit 12724112).
- The dominance census (RH_NDOM) will show whether the remote
  correction is systematically exceeding the local rate (artefact)
  or only at edge cells (plausible physics).
- If systematic dominance is confirmed: the fix is the missing
  third-resonance exp(-t'') extinction (Feldmeier & Nikutta), or a
  physically-motivated bound on DRC relative to the attainable
  remote flux — not an ad-hoc cap.
- The other three SOB validation models (`sob_diag`, `sob_jfnk`,
  `sob_tri`) completed cleanly — non-local coupling is the differentiator.

## 8. it66 resubmit -- kernel exonerated, lagged-map diagnosis (2026-11-05)

The 4-CPU resubmit (job 12729722) ran past the it65 crash point to
it66+ without a SIGSEGV, but the upstream divergence persists
unchanged: the OIV:2s2_4p_2Po[1/2] exempted demand at depths 9-11
oscillates between 1E20 and 1E39 per iteration, and `maxraw`
flip-flops ~1E15/1E22 on non-exempted rows in the same region.

Two new pieces of evidence change the picture:

- `SOB_RH_TEST_S1=T` shows `max dJ(S=1)` flat at ~0.2 while
  `max dJ/J` reaches 5E7.  The two-point kernel is a contraction in
  source space (remote dJ bounded well below the brightest emitter),
  so the blow-up is not geometric overcounting.
- An offline replica of the pair selection on the run's R,V grid
  reproduces the 4,830 pairs exactly and finds only 382 (8%) with a
  third resonance between the endpoints; the geometric weights at the
  problem depths are at the median and come mostly from adjacent
  cells.  Missing third-resonance extinction is real but cannot be
  the dominant inflation channel.

The remaining mechanism is the lagged feedback itself: DRC is
computed from the current populations and re-enters through S_L on
the next call; where remote absorption dominates the local rates the
population response amplifies the next computed DRC -- a
non-contractive fixed-point map.

Fix applied: `SOB_RH_RELAX` (VADAT, default 1.0) under-relaxes the
applied term across calls, applied = applied + w*(computed-applied).
The converged fixed point is unchanged -- the flag damps the
iteration path, not the remote flux.  The self-report now also
prints the (I,L) location of the largest DRATE and JRAT, plus the
computed-applied lag when relaxation is active.  Unit coverage:
`test_sob_rh_couple.py` T8 verifies applied_k = (1-(1-w)^k)*computed
across three phases; suite 14/14.

The resubmit runs with `SOB_RH_RELAX=0.2`.  If the lagged map still
does not contract, the next correction is path extinction along the
chord (third-resonance and continuum), not a cap on the term.
