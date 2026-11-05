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

## 6. Open questions

- Whether the NaN entered the depth-11 matrix block from a `SOB_RH`
  rate correction (overflow/0/0) or from a depth-11 population that was
  already degenerate when the matrix was assembled — needs a rerun with
  the extended diagnostics on an earlier checkpoint or a finiteness
  check inside the STEQ assembly.
- Whether a finiteness guard on the STEQ column before correction
  application (diagnosed abort rather than silent repair) is the right
  minimal defence — it converts checkpoint poisoning into a stop at the
  actual defect site.
- The other three SOB validation models (`sob_diag`, `sob_jfnk`,
  `sob_tri`) completed cleanly — non-local coupling is the differentiator.
