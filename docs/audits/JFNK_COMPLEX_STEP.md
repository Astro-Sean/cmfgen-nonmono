# Complex-Step JVP for JFNK - Design Note (W22)

Status: stub only. `[JFNK_FDCPLX]` (default F) is registered in
`rd_control_variables.f` and `control_variable_mod.f`; enabling it
fails closed in `JFNK_JVP_V1`. This note is the port plan.

## Why

The forward-difference JVP

    Jv ~= (F(x + eps*v) - F(x)) / eps

loses ~8 digits to cancellation once |eps*v| ~ 1e-8 relative to x.
CMFGEN's coupled residual runs a full radiation sweep per eval, so
eps cannot shrink much below 1e-8 without the difference drowning
in sweep noise - the observed FD-vs-stored mismatch (~17% on the
dominant column) is largely this truncation/cancellation floor.

The complex-step trick (Martins & Ning 2021; Arramy et al. 2024,
Sec 2.3, eqs. 43-45) perturbs in the imaginary direction:

    F(x + i*eps*v) = F(x) + i*eps*Jv + O(eps^2)
    Jv = Im[F(x + i*eps*v)] / eps          eps ~ 1e-200

No subtraction, no cancellation: the JVP is second-order accurate
at machine precision even for ill-scaled rows.

## What must change

1. Type switch: the eval copy of the coupled residual
   (`JFNK_RAD_SWEEP_V1` callback, cmfgen_sub.f eval region) must
   run in COMPLEX(LDP) for x, POPS, SE, and every rate/J nu term
   on the perturbation path. The real residual F(x) is already in
   `jfnk_fx`; only the imaginary part is new work.

2. Linearization boundary: F is linear in the perturbation only to
   O(eps). Any max/min/abs/sign in the eval path must use the
   "complex-safe" forms (max of complex args takes the real-part
   ordering; the standard `zmax = cmplx(max(re,im))` trick is NOT
   needed since eps is imaginary-only and Im max picks the same
   branch as Re max at eps->0 provided no tie).

3. JVP assembly: `JFNK_JVP_V1` currently does one real eval and one
   subtraction. Complex-step does one complex eval and a divide of
   the imaginary part - same count of sweeps per JVP, ~2x cost per
   eval (complex arithmetic), vastly better accuracy.

4. Fixed-unknown masking, sentinels (`JFNK_EVAL_BROKEN`), and the
   eval-budget counter all stay identical - they wrap the eval
   boundary, not the arithmetic.

5. The production CMFGEN routines do not change: only the eval
   copy in cmfgen_sub.f (the JFNK_IN_EVAL region) is duplicated to
   a complex twin, exactly as the real eval copy was made. This is
   the "shared eval" point of Arramy Sec 2.3 - the operator used
   for the JVP and the residual must be bit-compatible apart from
   the perturbation.

## Estimated scope

- One complex twin of the eval-region sweep (~1500 lines, mostly
  mechanical LDP->LDC edits plus min/max/abs review).
- A complex overload for `SE` population store reads in the eval.
- `JFNK_JVP_V1` switch on `JFNK_FD_COMPLEX` choosing
  `JVP = Im(F_eval(x + i*eps*v))/eps`, eps=1e-200.
- Verification: JVP must reproduce the analytic `J_BA` T-column
  within ~1e-10 relative (vs ~1e-4 achievable by FD), checked by
  the existing `JFNK_ROW_ANATOMY_V1` diag.

## Risk

Low physics risk (the flag is off by default and fails closed),
moderate maintenance risk: the complex twin must be re-synced when
the real eval copy changes - the same hazard the real copy already
carries, doubled. Recommend a shared INCLUDE file if ported.
