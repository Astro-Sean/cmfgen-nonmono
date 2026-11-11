# KNOP P1-T1b - global lagged radiation-map eigenanalysis - 11-Nov-2026

Follow-up to `KNOP_P1_T1_FIELD_SENSITIVITY.md`.  P1-T1 measured a *local*
single-cell store gain `g ~ 0.97` at `|Q| -> 1`, giving an effective
relaxed gain `(1-omega) + omega*g ~ 0.985` and a hypothesis of a
near-non-contractive lagged map.  P1-T1b builds the *global* map over
the full lagged `I_NEXT` store and measures its dense Jacobian spectrum
directly, to test whether the network gain approaches or exceeds 1 once
cross-frequency store sharing is included.

Result: **the lagged radiation-store map is nilpotent.**  The local
`~0.985` gain is a per-hop transmission along an acyclic frequency
ladder, not a loop gain.  No eigenvalue of the transport Jacobian is
nonzero; the only nonzero spectrum is the `(1-omega)` relaxation memory
term.  At the production `INEXT_RELAX = 0.5` the spectral radius is
`0.5`, not `~0.985`.  The radiation store cannot sustain a growing or
marginally-contractive mode on its own.

Runtime evidence: `unit_tests/bin/knray p1t1b` (driver `DRV_P1T1B` in
`unit_tests/f90/drv_knop_ray.f90`), products `p1t1b_*.txt` in the run
directory.  Manufactured model `NDA=8, NPA=2, NCFA=6`, one ray
(`NRAYS=1`), `SOLUTION_METHOD='INTEGRAL'`, `KN_RAY_TOL=2e-3`,
`KN_RAY_MAXPASS=60`, `NONMONO_INEXT_CAP=1e20`, `INEXT_RELAX` swept over
`{0.5, 1.0}`.

---

## 1. What was measured

The radiation map `F(X)` is defined on the lagged store state
`X = I_*_DEPTH_NEXT` restricted to the elements the production sweep
actually reads and writes:

- Compressing sparse depths only: `COMP_DEPTH_IDX(K) > 0` for `K` whose
  `J_PNT` lands on a valid ray point.  In this model three compressing
  depths (`KC = 1,2,3`).
- Writable store slots `s = 1..NCF-1` (`s = 1..5`).  Slot `NCF` is a
  read-only phantom boundary (never written); frequency 1 (`FIDX=0`)
  is not a writable lagged slot.  `FREQ_IDX_ITER` runs `0..NCF-1`, so
  frequency `L` writes slot `L-1`.
- Both channels `P` (outward) and `M` (inward), ray index `LS=1`
  (`NP=1`).

State dimension `NX = 3 depths x 5 slots x 2 channels = 30`, ordered
channel-major then slot-major then depth (`ELEM` dump in
`p1t1b_idx.txt`).

`MAP_F(X, XO, IRF)`:

1. Restores the slaved bookkeeping arrays (`I_*_STORE`, `I_*_PREV`)
   to a fixed post-cold-sweep snapshot so every evaluation starts from
   identical conditions.
2. Installs `X` into the active `I_*_DEPTH_NEXT` elements.
3. Runs the ordered frequency sweep `L = 1..NCFA` calling
   `FG_J_CMF_NONMONO_V1` with `INIT = (L==1)`, `NEW_FREQ = .TRUE.`,
   keeping the normal in-sweep relaxed store updates.
4. If `IRF=1`, calls `FG_KNOP_REFRESH_V1()` (the production post-sweep
   refresh); `IRF=0` is sweep-only.
5. Packs the active store elements back into `XO`.

Three map variants were differentiated:

| variant | map | `INEXT_RELAX` |
|---|---|---|
| `fullw05` | sweep + `FG_KNOP_REFRESH_V1` | 0.5 |
| `swp_w05` | sweep only | 0.5 |
| `swp_w10` | sweep only | 1.0 |

For each variant a dense `30x30` Jacobian was built by central finite
differences at two amplitudes (`1e-4` and `1e-5` times
`max(|X_j|, 1e-3)`), falling back to forward differences only where the
minus perturbation would go negative (none did - all columns central).
`DGEEV` (LAPACK/MKL) gave the spectrum; direct nonlinear iteration from
three perturbation directions (dominant eigenvector, all-ones, largest
outward element) recorded the transient.

A diagnostic-only module flag `KNOP_DIAG_NO_REGRID` (default `.FALSE.`)
suppresses a spurious ray-grid rebuild that the compression-set
hysteresis check fires on every `INIT` for this manufactured sigma
profile.  It changes nothing in production; without it the store is
zeroed each call and the map is trivially constant.

## 2. Headline result - the transport Jacobian is nilpotent

Store update (production `fg_j_cmf_nonmono_v1.f` ~4840):

```
I_*_DEPTH_NEXT(KC,LS,FIDX) = T2 + omega*(T1 - T2)
```

with `T2 = SNP` the pre-update store snapshot (= installed `X[s]`, since
each slot is written once per sweep) and `T1` the just-computed
intensity.  Hence

```
output[s] = (1-omega)*X[s] + omega*I_computed[s]
```

The lagged compressing read (`~4114-4152`) pulls
`I_*_DEPTH_NEXT(KC,LS,FIDX+NM_RAY)` with `NM_RAY = 1` at every
compressing ray point (verified `NM=1` at rays 4-7).  Frequency `L`
(`FIDX = L-1`) therefore reads slot `L` but writes slot `L-1`:

```
output slot s  depends on  input slot s+1      (read = write + NM, NM>=1)
```

Every compressing read is strictly slot-ascending.  There is no path by
which an output slot depends on an input at the same or a lower slot, so
the transport Jacobian `J_T` is strictly upper-triangular in the slot
ordering and therefore **nilpotent**: all 30 eigenvalues are exactly
zero, and `J_T^5 = 0` on the five-slot chain.

Measured (`p1t1b_Jswp_w10`, `omega=1` so `J = J_T` directly):

- `max|J_T| = 0.968`, all 30 `dgeev` eigenvalues `= 0` (real parts and
  imaginary parts both `0.0E+00`).
- Dominant entries (`|J|>0.01`, P channel only): each output slot `s`
  reads input slot `s+1` at `KC2` with gain `0.963-0.968`:

  ```
  out(P,s,KC1), out(P,s,KC2) <- in(P,s+1,KC2),  s = 1..4,  gain ~0.96
  ```

  Output slot 5 reads the phantom slot 6 (=0) and is constant; input
  slot 1 is never read.  The chain is a pure frequency-ladder shift.
- M-channel entries are `~1e-9` (inward lagged coupling negligible in
  this model); they cannot create a feedback loop in any case because
  the M read/write obeys the same `read = write + NM` slot rule.

## 3. The three spectra

`eig(J_F) = (1-omega) + omega*eig(J_T) = 1-omega` for every mode, since
`eig(J_T)=0`:

| variant | `J_F` structure | spectrum |
|---|---|---|
| `swp_w10` (`omega=1`) | `J_T` nilpotent | all 30 eigenvalues `0` |
| `swp_w05` (`omega=0.5`) | `0.5*I + 0.5*J_T` | all 30 eigenvalues `0.5` |
| `fullw05` (sweep+refresh) | refresh-determined | all 30 eigenvalues `0`, `J` exactly `0` |

The `0.5` eigenvalues at `omega=0.5` are **pure relaxation memory** -
the `(1-omega)*X[s]` diagonal - not a physical slow mode.  Spectral
radius `rho = 1-omega = 0.5`: the lagged store converges linearly at
`~0.5` per sweep (`~3` sweeps per digit), and at `omega=1` it is
finite-step (nilpotent).  This is the opposite of the P1-T1 `~0.985`
estimate.

Direct nonlinear iteration (`p1t1b_nl.txt`, 25 steps, perturbation
`delta = 1e-3*||X0||`) confirms the spectrum dynamically:

- Sweep-only `omega=0.5`: the step-to-step contraction ratio settles at
  `0.500015` for the dominant-eigenvector and largest-element
  directions (`RD` halves each iteration, `0.075 -> 2.1e-6` over 25
  steps).  The all-ones direction contracts more slowly early
  (`~0.6-0.9`, settling `~0.60`) as it projects onto the nilpotent
  flush, but still converges to `RD = 1.7e-5`.  No near-unity stall.
- Full map (`sweep+refresh`, `omega=0.5`): `F(X0+delta) = X0` on the
  first application for all three directions (`RD = 0` thereafter) -
  one-step convergence, consistent with `J = 0`.

### Reconciliation with the P1-T1 local gain

P1-T1's `g ~ 0.97` is `dI_solved/dI_NEXT` for one compressing cell - a
real, near-unity *per-hop transmission*.  But that transmission feeds
only the *next* store slot; the perturbation then propagates
`s -> s-1 -> ... -> boundary` and exits after at most `NCF-1` hops.  A
near-unity gain on an acyclic chain still has spectral radius zero.
P1-T1's `(1-omega)+omega*g ~ 0.985` treated the local gain as a scalar
loop gain; the correct global object is a nilpotent shift.  **The local
gain is real; the inferred near-non-contractive mode does not exist.**

### The full production map is seed-independent

`fullw05` returned `J = 0` *exactly* (0 entries `>1e-15`), and direct
nonlinear iteration converged back to `X0` in a single map application
from all three perturbation directions (`RN -> 0` after iter 1).  With
opacities held fixed, `FG_KNOP_REFRESH_V1` iterates the store to a
unique fixed point that depends only on the (fixed) opacities, not on
the lagged seed.  The radiation field is therefore **slaved to the
populations**: the lagged `I_NEXT` store carries no independent state
once the refresh runs.

## 4. Interpretation for the field-drift problem

The observed `GS_MAXCHG` / field-maturity stall cannot originate in the
lagged radiation-store iteration:

1. The transport Jacobian is nilpotent (`rho=0`); the relaxed map has
   `rho = 1-omega = 0.5`.  Perturbations flush out in `<= NCF-1` sweeps.
   There is no near-unity or amplifying radiation-store mode.
2. The full map (with refresh) is exactly state-independent
   (`J=0`): given fixed opacities the store is already at its fixed
   point.

Any persistent radiation-field drift must therefore be driven from
**outside** this sub-map - i.e. by the outer population-opacity
iteration changing the transport operator `T` (and hence the refresh
fixed point) each great iteration, or by an incompletely-converged
refresh/population solve.  This is consistent with the second-pass
audit's population-conditioning findings (Pillar 5) and reinforces the
standing verdict that KNOP's blocker is the coupled population system,
not the formal lagged solve.  It also means raising `INEXT_RELAX`
toward 1 is *not* dangerous for this map (it would only remove the
`0.5` damping and expose the already-nilpotent `J_T`); the relaxation
is a nonlinear-robustness choice, not a contraction requirement.

## 5. Caveats and scope

- **Fixed opacities.**  `F(X)` holds `ETA/CHI/ESC` (and hence the
  populations) fixed; population feedback is the separate outer map and
  is not included.  The nilpotent result is a statement about the
  radiation-store sub-map only.
- **Manufactured model.**  `NDA=8`, `NCFA=6`, one ray, `INTEGRAL`
  method.  The nilpotent structure is *structural* (`read_slot =
  write_slot + NM_RAY`, `NM_RAY >= 1` for compressing upwind), so it
  should persist for larger `NCF`/ray counts, but the eigenvalues were
  measured only here.
- **`NM_RAY >= 1` assumption.**  Nilpotence requires the compressing
  read to always target a strictly redder (higher) slot.  `NM_RAY` is
  the matched frequency offset to the characteristic foot; it is `1`
  here and `>=1` by construction for compressing flow.  An `NM_RAY <= 0`
  would introduce same/lower-slot feedback and is not observed.
- **Relaxation bypass.**  The store write falls back to `= T1` (no
  relaxation) when the old value `T2 <= 0`.  At the warm `X0` all active
  elements are positive, so the smooth `(1-omega)X + omega*T` branch
  governs the measured Jacobian; a state driven non-positive would see
  a different (gain-1) local map.
- **M-channel weakness.**  Inward lagged coupling measured `~1e-9`;
  whether that is generic or specific to this model's geometry was not
  isolated.  It does not affect the nilpotent conclusion.

## 6. Method validation

- Repeatability: `REP_MAXDIFF = 0`, `REP_REL = 0` (map is deterministic).
- Warm-up: 30 full-map iterations; residuals `||F(X)-X|| = 0` from iter
  2 on - `X0` is an exact fixed point of every variant.
- FD consistency: `max|J_a1 - J_a2| ~ 1e-11` on every column of every
  variant; zero columns required forward differences.
- `dgeev` `INFO=0`, `NEL=30` for all variants.
- Unit tests link the production `libsubs.a` / `.mod` objects directly;
  `KNOP_DIAG_NO_REGRID` defaults `.FALSE.` so production is unchanged.

## 7. Bottom line for P1-T1b

The question was whether the *network* lagged-map gain approaches or
exceeds 1.  It does not: the global transport Jacobian is nilpotent
(spectral radius 0), and the relaxed production map has spectral radius
`1-INEXT_RELAX = 0.5`.  The P1-T1 `~0.985` figure was a local
transmission coefficient on an acyclic chain, not a loop gain.  The
radiation lagged solve is unconditionally stable and fast; the
field-maturity stall must be sought in the population-opacity outer
iteration (P5) or refresh convergence, not in `I_NEXT` store dynamics.
