# KNOP I_PREV Ablation Audit (P1-T1e, PART 5)

Question: are the weak same-slot cross-depth Jacobian edges measured in
the P1-T1b/P1-T1d manufactured-model experiments genuinely mediated by
the expanding `I_PREV` frequency-feedback channel?

Method: a diagnostic-only `FREEZE_I_PREV` flag ablates that channel,
and the identical `p1t1b` sweep Jacobian is recomputed frozen vs
baseline. Evidence class: **[MANUF]** manufactured fixed-opacity model.
This tests the *origin of a particular Jacobian dependency*; it does
**not** establish production nilpotency or convergence.

## 1. The intervention

`plane/fg_j_cmf_nonmono_v1.f` (commit `1d5aa693`): new module flag
`FREEZE_I_PREV`, default `.FALSE.`. When set, the `NEW_FREQ` branch
skips the `I_*_STORE -> I_*_PREV` copy:

```fortran
ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
  IF(.NOT. FREEZE_I_PREV)THEN
    I_P_PREV(:,:)=I_P_STORE(:,:)
    I_M_PREV(:,:)=I_M_STORE(:,:)
  END IF
END IF
```

`I_*_PREV` is written in exactly two places: INIT zeroing and this
`NEW_FREQ` copy (verified by source grep - no other writes exist). With
`FREEZE_I_PREV=T`, `I_*_PREV` keeps its INIT-zeroed value for the whole
sweep, so expanding cells evaluate
`SOURCE_PRIME = SOURCE_RAY + Q*I_*_PREV` with `I_*_PREV = 0`, i.e.
`SOURCE_PRIME = SOURCE_RAY`. The expanding frequency-feedback term is
ablated while every other part of the operator is untouched:

- `I_*_STORE` still updates every call (untouched).
- `FREQ_IDX_ITER` still increments on `NEW_FREQ` (untouched).
- Compressing `I_*_DEPTH_NEXT` donor reads unchanged (untouched).
- Frequency matching `NM_RAY`/`DLMN`, ray geometry, transfer
  coefficients, refresh, convergence criteria: all unchanged.
- The perturbation is applied to the same store element via the same
  `MAP_F`; the only difference between runs is the `I_PREV` ablation.

`unit_tests/f90/drv_knop_ray.f90`: new case `p1t1bf` reuses the
`p1t1b` Jacobian machinery (NCASE=58) with `FREEZE_I_PREV=T`, restricted
to the sweep-only `swp_w05` variant (IV=2, `INEXT_RELAX=0.5`).

## 2. Runs

Build: `knray` via `unit_tests/common.py::compile_and_link_prod`
against production `lib/libsubs.a` (`fg_j_cmf_nonmono_v1.o` rebuilt by
`make DEBUG=nodebug`).

- Baseline: `knray p1t1b`  -> `/tmp/kb/p1t1b_Jswp_w05_a{1,2}.txt`
- Frozen:   `knray p1t1bf` -> `/tmp/kbf/p1t1b_Jswp_w05_a{1,2}.txt`

State vector `X` = 30 lagged-store elements = 3 compressing depths x
5 frequency slots x 2 channels (P then M), slot-major packing per
`p1t1b_idx.txt`. Jacobian `J(i,j) = dF_i/dX_j` by central FD.

## 3. Measured result

Same-slot (`dslot = 0`) off-diagonal entries, max |J| per channel:

| channel | baseline `p1t1b` | frozen `p1t1bf` |
|---|---|---|
| P | **4.795e-04** | **0.0** |
| M | **3.264e-09** | **0.0** |

Per-edge (P `out(KC1,S) <- in(KC2,S)`, all 5 slots):
4.795/4.769/4.744/4.720/4.695 e-04 baseline -> exactly 0 frozen.
M `out(KC3,S) <- in(KC2,S)`: 3.26/2.84/2.49/2.20/1.95 e-09 -> 0.

Both finite-difference amplitudes (`a1`, `a2`) give identical results;
amplitude repeatability is ~1e-11 (the FD noise floor), so the edges are
real and not numerical noise.

Off-diagonal nonzeros: 30 (baseline) -> 16 (frozen). The 14 removed
edges are precisely the `I_PREV`-mediated ones (10 same-slot + a few
weak backward edges); the 16 persisting are the direct compressing
`I_*_DEPTH_NEXT` donor reads, which are correctly unaffected.

Eigenvalues (`dgeev`, `p1t1b_eig_swp_w05.txt`): all ~0.5
(= `1 - INEXT_RELAX`) in both runs - the relaxed-map diagonal
dominates; ablating the weak `I_PREV` edges does not move the spectral
radius, as expected for ~5e-4 perturbations.

## 4. Interpretation

The same-slot cross-depth dependencies are **causally `I_PREV`-mediated**:
they exist only because freq `k`'s intensity is copied to `I_PREV` and
read by freq `k+1`'s expanding cells. Removing that single copy
drops them to exactly zero. This confirms the P1-T1d structural
elimination argument (expanding `I_PREV` is the only same-slot
pathway) experimentally.

**Scope limits (explicit):**

- This is the *manufactured* fixed-opacity model, sweep-only
  (`swp_w05`), `INEXT_RELAX=0.5`. It says nothing directly about
  production geometry, the full map (sweep+refresh), or
  population-coupled convergence.
- `FREEZE_I_PREV` sets `I_PREV = 0` (a *cold* expanding field), not a
  warm unperturbed baseline. It therefore ablates the entire expanding
  feedback channel rather than freezing it at baseline values. The
  edges going to *exactly* zero shows the mediation is total, but a
  warm-baseline freeze (per-frequency `I_PREV` snapshots) would be a
  more surgical variant if a nonzero frozen field is ever needed.
- `I_PREV=0` changes the base solution (`max|dJ| ~ 1.6e-2` vs
  baseline), so this is a channel-ablation causal test, not a
  perturbation about the same fixed point. The binary question -
  "do the edges require `I_PREV`?" - is nonetheless answered cleanly.

## 5. What this does and does not establish

- **Established [MANUF]:** the same-slot Jacobian edges originate
  exclusively in the expanding `I_PREV` channel. They are not FD noise,
  not parabolic-stencil leakage, not a compressing-read artefact.
- **Not established:** production-geometry nilpotency (still needs the
  PARTS 2-4 checkpoint + structural graph + SCC experiment), refresh
  convergence, or population-coupled stability. None is inferred here.

## 6. Remaining work

- Production structural dependency graph + Tarjan SCC + targeted FD
  (PARTS 2-4): still blocked on a real non-monotonic checkpoint; the
  `FGINEXT_GS` file lacks the opacity/source arrays needed to rebuild a
  sweep offline.
- Refresh `FLD` internal-residual instrumentation (PART 6): specified,
  not yet implemented.
