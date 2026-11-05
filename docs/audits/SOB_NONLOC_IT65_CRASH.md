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

## 6. Producer analysis — how NaN enters the depth-11 RHS

`SOB_RH_COUPLE` (subs/sob_rh_couple.f) returns `DRC(ND,NSIM)` which
`sub_sob_line_v3.f:324-328` subtracts directly into the SE residual:
`T2 = ETAL_MAT*ZNET_SIM - DRC_SIM`, then `STEQ(MNUP,K) -= T2/T1`,
`STEQ(MNL,K) += T2/T1`, `STEQ_T(K) -= T2*T3`. A single non-finite
`DRC(11,L)` poisons all of depth 11's residual rows including the
temperature row -- which is exactly how `T(11)` went NaN in the solve.

Three holes let a non-finite DRC through:

1. `DRATE = DJ*CHIL_MAT/BB_COR` (sob_rh_couple.f:250) is guarded only
   by `BB_COR .EQ. 0`. `BB_COR = T1/T2` (cmfgen_sub.f:4078-4080) is a
   Planck rescaling ratio; at cold depths `T1 = FL^3/(exp(HDKT*FL/T)-1)`
   collapses exponentially for high-FL lines, so `BB_COR` can land in
   the denormal-but-nonzero range -> `DRATE` overflows to Inf.
   Diagnostics agree: `max dJ/J` stayed ~1-20 while `max dRATE` hit
   2.3E8 -- the blow-up is in the /BB_COR conditioning, not the kernel
   sum itself.
2. `DJ .LE. 0` (line 248) never catches NaN (NaN comparisons fail), so
   a NaN DJ from an Inf*0 product in the partner sum also passes.
3. `DRC` accumulates with no finiteness check; `RH_MAXDZ` logs only the
   per-correction max, never the accumulated value, so an Inf build-up
   is invisible in the reported diagnostics.

Same guard weakness class at sub_sob_line_v3.f:362 (`SRAT` divides by
`BB_COR` guarded only `.GT. 0`).

Driver sequencing: `dRATE` printed 1E7-1E8 across it55-64 -- the
conditioning degraded steadily until a DRATE crossed to Inf/NaN at
it64's final solve. Outer-atmosphere population explosion (exempted
demand 2.1E20% at depth 8) was running in parallel and is the likely
reason BB_COR/S_L reached the bad range at depth 11.

## 7. Open questions

- Confirm which hole fired first (denormal `BB_COR` overflow vs
  `Inf*0` in the `DJ` sum) — a finiteness check on `DRATE`/`DJ`
  reporting line/depth would name it on the next reproduction.
- Whether a finiteness guard on `DRC` before `STEQ` injection
  (diagnosed abort rather than silent repair) is the right minimal
  defence — it converts checkpoint poisoning into a stop at the
  actual defect site.
- The deeper physical question: whether `dRATE` ~1E8 corrections are
  legitimate thick-coupling pumping or a discrete-kernel artefact of
  the RH78 two-point sum (the `mu^2<=1` cut bounds `D2` only relative
  to velocities; `UPS=1/SQRT(D1*D2)` can still be huge for close
  pairs).
- The other three SOB validation models (`sob_diag`, `sob_jfnk`,
  `sob_tri`) completed cleanly — non-local coupling is the differentiator.
