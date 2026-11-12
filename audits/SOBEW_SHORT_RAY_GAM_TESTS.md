# SOBEW-family short-ray GAM: test description and results

Base commit (pre-fix): `cd0f4374f1be07645a6a79f205a48f4285ea864b`.
Source correction committed on main: `099435c7` (`sobew`),
`e5e0aae9` (`sobew_grad`), `b93f3377` (`sobew_grad_v2`),
`686a5718` (`sobew_grad_v3`). HEAD at review close:
`e6ed21a74dcc8760cdf65f2c1c91978851b24519`.

## Test file

`unit_tests/python/test_sobew_short_ray.py`. Links
`unit_tests/f90/drv_sob_floor.f90` against:

- new: production objects `subs/sobew.f`, `sobew_grad.f`,
  `sobew_grad_v2.f`, `sobew_grad_v3.f`, `sobjbar_sim.f`;
- ref: the same five routines compiled from `cmfgen_reference`
  (pristine upstream: identical stale-GAM bug in all four EW routines,
  no denominator floor, and a *signed* denominator in `sobew_grad.f`).

Throughout, "floored" means this implementation's own
`max(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4)` denominator regularisation
- a numerical safeguard of this codebase, not part of the published
Dessart, Audit & Hillier 2015 expression.

Helpers (production objects per AGENTS.md policy):
`zalongp.f`, `nordtau.f`, `tcompd.f`, `xvecd.f`, `thomas.f`,
`thomas_rh.f`, `exponx.f`, `set_kind_module.f`,
`mod_trap_derivatives.f`, `term_io.f`.

Driver routine keys: `sobew`, `sobewg1` (SOBEW_GRAD), `sobewg2`,
`sobewg3`. `sobewg1` needed a new driver branch (signature uses `DIF`
logical instead of `INNER_BND_METH`).

Geometry: `ND=25, NC=4, NP=29=ND+NC`; `LS<=NC` -> `NI=ND`, shell ray
`LS>NC` -> `NI=ND-(LS-NC-1)`, so `LS=27,28,29` have `NI=3,2,1`. The
Python side reproduces `EXPONX`, the `DTAU` thermal-depth recurrence,
and the `(T3,T4)` brackets exactly, so each `LINE_FLUX(i)` is
predictable under both current-ray (floored) and stale-ray
(unfloored) semantics.

## Checks: 8 LINE_FLUX checks x 4 routines + 3 FORCE_MULT checks
x 3 gradient routines = 41, all pass

| # | content |
|---|---------|
| T1 | `sigma=-1.05`, NI=1 ray isolated by weighting: new/ref LINE_FLUX(1) ratio vs current-ray vs stale-ray predictions (new GAM=3.750 vs stale 8.574) |
| T2 | same at NI=2, LINE_FLUX(2) (new 3.212 vs stale 4.828) |
| T3 | `sigma=+2`: depths `I>=3` (fed only by NI>2 rays) bit-identical new vs ref |
| T4 | swap last two physical rays+weights: per-case check vs prediction for both orderings; asserts orderings discriminate |
| T5 | `sigma=+2` monotonic: current-ray discriminator (no floor) |
| T6 | `sigma=-1`: NI=1 discriminator, 1e-4 floor retained |
| T7 | negative `GAMH` at `sigma=+2`: EXPONX gain + sign retained |
| T8 | `sigma=-2` (< -1): NI=1 floored current-ray GAM |
| F1 | `sigma=-1.05` (-0.5 for `sobewg1`), NI=2 ray isolated: `fm(1)` new/ref ratio vs the EXPONX(GAM) ratio; fm(i>=2) must be exactly 0 |
| F2 | `sigma=+2`, NI=2 isolated: `fm(1)` discriminator |
| F3 | `sigma=+2`, `chil=-1e9` (negative GAMH): `fm(1)` finite, negative, discriminated |

### FORCE_MULT discriminator (F1-F3, sobewg* only)

`FORCE_MULT(i) = K*R(i)^2*CHIL(i)/ESEC(i) * GLINE(i)` with
`GLINE(i) = sum_ls TB(i,ls)*EXPONX(GAM(i,ls))*HAQW(i,ls)` and
`K = 16e20*pi^2/RLUM/3.845e33` shared by all three gradient routines
(`sobew_grad_v2.f:241,300-302`; identical in v1 and v3). With a single
weighted ray, `TB(i)` (the transfer solution, independent of GAM) and
the K-factor cancel exactly in the new/ref ratio, leaving the pure
`EXPONX(GAM_new)/EXPONX(GAM_stale)` prediction. A buggy build returns
meas=1.0 and fails (verified: meas=1.0000 vs pred=1.5161).

Depth choice: the discrimination runs at `i=1`, not the tangent depth.
A tangent ray has `TB=0` at its deepest point under the hollow inner
boundary, so `GLINE(2)=0` on the isolated NI=2 ray and `fm(i>=2)` is
exactly zero - asserted inside F1. The NI=1 ray (tangent at R(1))
gives no force signal at all, so the force checks use NI=2; NI=1 is
already covered by the LINE_FLUX checks.

The driver initialises `fm = 0.0_ldp` before each `SOBEW_GRAD*` call
(`FORCE_MULT` accumulates without internal zeroing).

For `sobewg1` at `sigma<-1` the reference produces NaN `fm`
(signed denominator) and no absolute force prediction is possible
because `TB` is not independently computable; F1 for `sobewg1`
therefore runs at `sigma=-0.5`, where `1+sigma*mu^2 > 0` everywhere so
signed==ABS and the reference is interpretable.

Two discriminator modes:

- **ratio**: measured `LF_new(i)/LF_ref(i)` vs predicted
  `S_current/S_stale`. Used whenever the reference output is finite and
  the reference denominator semantics match the stale model.
- **abs**: absolute `LF_new(i)` vs `T1_norm * 0.5 * S_current *
  ETAL(i) * R(i)^3`, with `T1_norm = 13.19868*2.997924e-12/CF_printed`
  (the code scales `CONT_FLUX` by `R1^2*13.19868` before dividing).
  Agreement ~2e-9 relative. Used for `sobewg1` at `sigma<-1`.

## Why sobewg1 uses the abs fallback for sigma < -1

Upstream `sobew_grad.f` computes `GAMH/(1+sigma*mu^2)` with a **signed**
denominator (no `ABS`), unlike upstream `sobew.f`/`sobew_grad_v2.f`/
`sobew_grad_v3.f` which use `ABS(...)`. For `sigma<-1` the reference
v1 output is corrupted by sign-flipped GAM (NaN at `-1.05`, +1e33
overflow bursts at `-2`) - independent of the stale-ray bug. The
new-vs-ref ratio is not interpretable there, so `disc_check` falls
back to the absolute current-ray prediction, which still detects
stale-ray dependence (a buggy build deviates by ~4x).

At `sigma=+2` and `sigma=-1` signed==abs for `mu^2<=1`, so the ratio
mode remains valid and is retained.

## Why the test detects stale-ray dependence

- Predictions are computed under *both* semantics; measured values must
  match the current-ray one. A buggy build would produce
  `meas/pred = 1/2.233` (T1), `1/1.454` (T2), `1/0.982` (T5),
  `1/13.5` (T8) - all far outside the 1e-3 tolerance, so the suite
  provably fails on the pre-fix implementation (the reference binary
  itself embodies that failure mode).
- `isfinite()` alone would pass on the buggy build.
- Case files are routine-independent; identical inputs run through all
  four routines.

## Harness bug fixed during development

The driver passed `aqw(NDMAX=400,NPMAX=40)` (and `haqw`, `chilm`,
`etalm`, `bbcor`, `znet`, `vb`, `vc`, `betac`) to `(ND,*)` / `(ND,NSIM)`
dummies - column-major stride mismatch scrambled every weight
(invisible under uniform weights, which masked it in the original
tau-floor suite). Now allocated at exact `(nd,np)` / `(nd,nsim)`. This
was a test-harness defect only; production code is untouched.

Also corrected while developing the absolute check: `s_sum` omits the
0.5 accumulation prefactor (harmless in ratios) and `lf_pred` must use
`R(i)^3`, `ETAL(i)`, and `cf_printed` directly.

## Commands and results

```
cd unit_tests
python3 python/test_sobew_short_ray.py   # 41 passed, 0 failed
python3 python/test_sob_tau_floor.py     # 17 passed, 0 failed
```

`test_sob_tau_floor.py` gained a `LINK_ONLY` list (`sobew_grad.f`) so
the shared driver resolves the `sobewg1` branch without adding a new
check row.

## Limitations

- `sigma<0` new-vs-ref ratios conflate the fix with the missing ref
  floor; `sigma=+2` checks isolate the fix alone.
- The `NI=1` ray has `TB(1)=0`, so a small core-ray weight is retained
  to keep `CONT_FLUX` nonzero.
- `NI` is positional; the ordering test uses non-tangent `P` for the
  last two rays deliberately.
- The force checks are ratio-mode only: an absolute `fm` prediction
  needs `TB(i)` per ray, which would mean re-implementing the
  two-point and `TCOMPD` transfer solves in the harness. The ratio
  cancels `TB` exactly, so the check is still a genuine numerical
  discrimination rather than a finite-output test.
- The NI=1 ray cannot be exercised through `fm` (its only depth is
  the tangent point, where `TB=0`); NI=1 coverage stays on the
  LINE_FLUX checks.
- `sobewg1`/`sobewg3` have no production callers; they are exercised
  for completeness and consistency, not because a live path failed.
