# KNOP solver stall diagnosis -- m15 2.82d runs, 31-Oct-2026

Evidence from the scheduled runs (OUTGEN in
`cmfgen_validation/m15_2p82d_KNOP_{fresh,cold,warmS}/cmfgen/2.82days/`).

## Outcomes

- `warmS` (resumed run): it866, pinned at 75% for hundreds of great
  iterations at ~412 s/it, then crashed on a forrtl-61 format mismatch
  in the DBG-T diagnostic (`cmfgen_sub.f:1494`, fixed in 53eb0e88).
- `fresh`: it23 at 75%, hit the 2 h wall twice (~620 s/it).
- `cold`: real numerical failure -- NaN CHI_RAY (chi_hat) and ETA_RAY
  at freq index 1, ray point 2, `GAM=-4.96e-11`, `SIGMA_RAY=-1.04`.
  PHASE 2 first-failure diagnostic + JNU_CLAMP_DIAG written.

## The stall mode (warmS, and why KNOP struggles here)

Three coupled observations from the live gates:

1. **The lagged FG field never contracts.** Every sweep prints
   `fg d ~ 6e-2..1.5e-1, bad = 250-630K cells >50%, churn = 1-20M cells
   >5%`, alternating up and down -- a bounded oscillation of the
   lagged `I_NEXT` store, not divergence and not contraction.
   `capped=0` (no non-finite writes), so the sweep is internally
   healthy; the fixed-point map itself is oscillating.

2. **The population channel is off the charts.** `pops demand/applied
   = 1.000E+07 / 1.000E+07` -- the correction demand sits AT the
   1e7 trust cap every iteration. `worst margin ~119-129 dex (exm)`:
   the ex_maxch channel is ~120 dex above its gate, meaning some
   level wants a correction of order 10^120-fold. This is a
   runaway-level mode: a trace (or near-empty) population is being
   driven to a huge value by the BA solve, capped, then driven again
   next iteration. The cap prevents physical damage but guarantees
   the population exit gate can never pass.

3. **The residual is localized in the deep interior.** `res depth`
   lists the top-5 residuals: depths 64-70 (the innermost cells of the
   70-point grid) hold ~60-70% of the column sum (depth 65: 24%,
   66: 21%, 64: 15%, 70: 10%, 69: 9.6%). se residual ~4e13 total,
   ~1.5e13 excluding worst columns.

Diagnosis: the m15 non-monotonic model's deep interior is where the
reverse-shock compression puts the strongest negative-sigma cells,
so the deepest cells are simultaneously (a) optically thick and
line-blanketed -- extreme opacity/emissivity sensitivity to T and
pops -- and (b) in the compressing branch of the Knop frequency
split, where J depends on the lagged next-frequency store. A
population oscillation at those depths changes the source function
enough per great iteration that the lagged I_NEXT map's update is
comparable to the field itself -> sustained >50% per-sweep changes in
300-600K ray cells. The field oscillation feeds back into the
populations (J-dependent rates), closing the loop. The system is in
a stable period-2/limit-cycle of the coupled map, not approaching a
fixed point. Neither the sweep drift gate nor the pops gate can
close, so the run sits at 63-75% indefinitely.

Why SOB-esali converged and KNOP does not: the SOB+LAM_ES path uses
the same FG solve but with the population gates already satisfied on
a different (weaker-blanketing) effective problem; the esali recipe
also converged only at 4.6% under a 5% threshold, i.e., the same
deep-interior mode exists but the exit criterion swallowed it. The
KNOP runs are not wrong -- they are reporting a genuine failure of
the coupled fixed point to contract at the current damping/relaxation
strategy on this model.

## The cold-run NaN (separate, real defect)

`GAM=-4.96e-11` at `SIGMA=-1.04` -- a near-resonance ray cell at the
outer boundary. chi_hat and eta arriving as NaN means the NaN entered
upstream of the ray solve (interpolated CHI/ETA on the ray grid, or
a division in the coefficient assembly). Worth a targeted debug run:
the zero-GAM band at the model edge is where the audit's floor
analysis predicted the only unregularized singular behavior. This is
a robustness defect, distinct from the stall.

## Recommendations (ordered)

1. Identify the runaway level: print the argmax population of the
   ex_maxch channel (ion, level, depth) on the demand>=1e7 iterations.
   If it is a trace level at depths 64-70, the fix is a trace-demand
   floor policy decision, not a solver change.
2. Check whether the deep cells are on the compressing branch for
   most frequencies -- if so, increase JNEXT_RELAX damping there, or
   tighten MOM_DAMP_NM locally, as a controlled experiment.
3. Localize the cold-run NaN (CHI/ETA interpolation onto ray grid vs
   the near-zero-GAM coefficient path).
4. The forrtl-61 DBG-T format bug at cmfgen_sub.f:1494 already has a
   committed fix (53eb0e88); confirm the resubmitted job uses a build
   containing it.
