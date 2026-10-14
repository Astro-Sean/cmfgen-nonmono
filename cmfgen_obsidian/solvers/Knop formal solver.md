# Knop formal solver (`DO_NONMONO_KNOP`)

The Knop, Hauschildt & Baron (2009) comoving-frame solver for
arbitrary velocity fields - the heart of this tree. Audited
term-by-term against Papers I+II in [[KNOP_MATH_AUDIT]];
convergence strategy vs the literature in
[[KNOP_CONVERGENCE_STRATEGY]].

Flag: `DO_NONMONO_KNOP` selects the Knop population solver in VADAT.
The main VADAT population modes remain mutually exclusive with
`DO_NONMONO_SOB`. As of 13-Oct-2026, the CMF_FLUX OBS Knop formal-J path
has been removed while that code was restored to the stock CMFGEN
reference. `prepare_CMF_PARAM.py` rejects Knop population snapshots for
OBS rather than emitting a spectrum with a different formal J. This does
not remove the Knop population solver from the main CMFGEN run.

## Prerequisites (hard-enforced in both readers, 03-Oct-2026)

`rd_control_variables.f` (VADAT) now STOPs, not warns, on:

- `GLOBAL_LINE` other than `'BLANK'`, or `GLOBAL_LINE='NONE'` with
  any ion present not at `TRANS_*='BLANK'` (NONE needs a TRANS_* key
  per species slot - a missing key aborts in GET_KEY_STRING first).
  Anything else routes lines outside `COMP_J_BLANK`, where the Knop
  solver lives.
- `LAM_SET=T` (rewrites edge lines to SOB and resets
  GLOBAL_LINE_SWITCH to NONE at run time).
- `USE_FIXED_J`, `USE_FORMAL_REL`, `USE_J_REL`, `USE_DJDT_RTE`,
  `PLANE_PARALLEL`, `PLANE_PARALLEL_NO_V` - those paths bypass
  `FG_J_CMF_NONMONO_V1`.
- `INC_GRID=T` / `ACCURATE=T` (extended grid uses the monotonic moment
  solver) - pre-existing check.
- `GAMRAY_TRANS='RAD_TRANS'` (05-Oct-2026): rejected under both
  non-monotonic modes. `GAMRAY_SUB_V3` / `CMF_FORMAL_REL_V4_GAM` is a
  signed-sigma solver and receives the unclamped gradient under KNOP -
  a shock at the outer edge
  (sigma(1)<-1) also makes its beta-law extension diverge. Both modes
  are rejected because cached and NONLOCAL deposition tables are also
  interpolated on a multi-valued V coordinate.
- `GRID_METH='VEL_SPACE'` in ADJUST_R_DEFAULTS is rejected inside
  `ADJUST_R_GRID_V4` under either nonmono flag: `DO_VEL_REGRID_V2`
  interpolates R(V), which assumes a monotonic velocity field.
  TAU_SPACE, FULL_R and OB_ONLY remain allowed.
- `DO_BA_RES_ACC=F` warns only: the verified config uses
  residual-accepted BA steps.

`rd_cmf_flux_controls.f` no longer reads or dispatches
`DO_NONMONO_KNOP` in CMF_FLUX. The OBS stage is stock for monotonic
models and supports the separate `DO_NONMONO_SOB` line-render branch.
Knop population models need a separately validated OBS solver before
they can be reported as complete spectral models.

## The chi_hat discretization (Paper I, nu-form)

For the equation $dI/ds = \eta - \chi I + \Gamma\, dI/d\ln\nu$ with
$\Gamma = (V/R)(1+\sigma\mu^2)/c$:

$$ \hat\chi = \chi + \frac{|\Gamma|}{d\ln\nu} \quad(\text{always} > 0),
\qquad
Q = \frac{\Gamma}{d\ln\nu\,\hat\chi}, \qquad
S' = \frac{\eta}{\hat\chi} + |Q|\, I_{\rm up} $$

![[fig_chihat.png]]

The physical picture: photons advected into a cell from the adjacent
frequency act like *extra absorption* locally. Folding their arrival
rate into the opacity keeps it positive; the photons themselves stay
explicit in the source term `|Q|*I_up`, so nothing is lost - only the
bookkeeping changes.

- **Sign of $\Gamma$ selects the upwind frequency**: expanding cells
  read the live `I(l-1)` computed this sweep (Gauss-Seidel);
  compressing cells read `I_NEXT(l+1)` from the **lagged store**
  written by the previous great iteration (Jacobi, one-iteration lag).
- `d ln nu` uses the **upwind-side spacing**: `dLOG_NU` expanding;
  `DLN_N_SAV = ln(nu_k/nu_{k+1})` compressing (`SET_FG_NU_NEXT` per
  call) - necessary at ~10x line-edge spacing jumps.
- Then the march is stock V13: DTAU, E-functions, A0-A4, monotonic
  slopes, same boundary options (DIFFUSION / PNT_SRCE / ZERO_FLUX).
- `xi=1` (fully implicit wavelength coupling) throughout - what
  Paper II requires for strong coupling:

> [!quote] Hauschildt & Baron 2004, A&A 417, 317 ([DOI](https://doi.org/10.1051/0004-6361:20034473), [astro-ph/0401164](https://arxiv.org/abs/astro-ph/0401164))
> "the interpretation of the entire wavelength discretization as a
> source term is not valid for all optical depths"
- Red-edge boundary: `I_NEXT(:,:,NCF)` stays zero => photons advected
  below the grid leave it (zero-incident wavelength BC). The only
  place chi_hat's absorption lacks its source partner - a documented
  BC choice.
- Thick outer boundary (`THK_CONT=T`): the beta-law extension floors
  its exponent source at `sigma >= -0.999` (05-Oct-2026). Under KNOP
  `SIGMA` is not re-clamped, so a shock at the outer edge would have
  given `BETA<0` and a diverging `V_EXT`; the floor yields a
  near-constant-V extension and prints a warning - the thick BC then
  rests on an approximation and the model should be checked.
- Aberration `3*a_nu*I` is absent - inherited from stock CMFGEN's
  nu-form equation, deliberate.

### Deviations vs the paper (documented, deliberate)

- **Diagonal ALO only**: `ALO_DJDETA`/`ALO_DJDCHI` accumulate the true
  diagonal sensitivity of the computed J (`dAV/dS' = 0.5*(A2(K-1)+A2(K))`
  with the same `JQW` quadrature weights); `do_var_cont.f` zeroes the
  `TX/TVX` off-diagonals under KNOP and fills only the diagonal. The
  paper's Lambda* is tridiagonal in wavelength - a convergence-rate
  choice, same fixed point.
- Domain decomposition -> OpenMP over rays (characteristics are
  independent) + the sparse `COMP_DEPTH_IDX` compressing-depth store.

## The lagged store and the refresh

`I_NEXT` (and `J_NEXT_FRQ`/`H_NEXT_FRQ` on the MOM path) hold the
"not-yet-computed upwind" intensities. One lagged iteration per great
iteration is too slow - red-to-blue information would take ~NCF
iterations to cross a compression. So after each frequency sweep:

> 05-Oct-2026 fix: `FREEZE_I_NEXT` actually holds during JFNK evals
> now. `IF(INIT) FREEZE_I_NEXT=.FALSE.` in `fg_j_cmf_nonmono_v1.f`
> fired on every first-frequency call - including inside eval sweeps -
> so the store was being re-advanced at perturbed populations on every
> residual/JVP eval, making F(x) eval-order-dependent and injecting an
> h-independent drift into the FD JVP (the observed |Jv|~1e18 vs
> |BA.v|~1e13 mismatch). The INIT-clear is now gated on
> `.NOT. JFNK_IN_EVAL`. The MOM-side store (`J_NEXT_FRQ`) has no
> freeze at all - dormant while `USE_LAM_ES=T`.

**`FG_KNOP_REFRESH_V1`** re-solves each ray's (space, frequency)
intensity field on the sweep-captured opacities until the lagged
store reaches its fixed point:

- Per-frequency capture during the sweep: `KN_CHI_EXT`/`KN_ETA_EXT`
  (single precision on the extended depth grid), geometry, boundary
  data - all pass-invariant (coefficient pack, commit `41083d9`).
- `KNOP_RAY_SOLVE_V1` re-runs the same short-characteristics
  discretization per pass; fixed point identical by construction.
- Guards: `DO_KNOP_REFRESH`, `HAS_COMPRESSING`, complete sweep,
  INTEGRAL method only (TRIDIAG formal path not refreshed).
- Convergence machinery: adaptive per-ray omega (~0.5 damp ->
  ~1.45 over-relax, keyed on same-direction residual history),
  divergence abort (`1E3*MIN_RES`), `1E40` retention scan.
- Per-sweep exits (13-Oct-2026 census): `0` exact (no compressing
  cells), `1` residual < `KN_RAY_TOL`, `2` bulk bad-fraction
  (<=50 cells or <=0.1% of significant cells, after pass 6), `3`
  same-direction stall (two consecutive same-direction ratios >0.9
  after pass 8), `4` pass cap `KN_RAY_MAXPASS`, `5` diverged. The
  summary prints a `ray exits` count per class plus the stall/cap
  classes' max residual and bad-cell share, so the refresh log shows
  which stop rule leaves above-tolerance cells on the release gate.
- Store write-back uses `INEXT_RELAX=0.5` under-relaxation and the
  `NONMONO_INEXT_CAP=1e20` overflow guard.

### Guarded Aitken: `KN_ORD_AITKEN` (default F)

Aitken delta-squared extrapolation of the ordered pass iterate,
applied per cell only when the cell still moved by more than
`KN_RAY_TOL`, its two consecutive pass deltas contract
(`|d_n| < |d_{n-1}|`), and the extrapolated value is finite,
positive and below `NONMONO_INEXT_CAP`. Passes 1-2 build the delta
history. Aimed at the marginally contractive `|Q| -> 1` compression
cells that drive the same-direction stall exit; preserves the fixed
point if it converges. A/B test vs the `KN_ORD_AITKEN=F` control
(m15 2.82d, 13-Oct-2026): the per-cell depth-1 secant amplified the
seed transient - 25-26 diverged rays per iteration vs 0 in the
control, residual max 10.5->853 vs 1.0->0.17. Rejected as a
production path; superseded by `KN_ORD_AA`.
#status/failed

### Safeguarded Anderson: `KN_ORD_AA` (default F)

Type-II Anderson acceleration of the ordered pass iterate, mirroring
the runtime-tested MOM refresh AA (`mom_j_cmf_nonmono_v1.f`) plus a
residual-acceptance rollback the MOM version lacks. History pairs are
`(x_i, r_i=g(x_i)-x_i)` per relaxed march; `KN_ORD_AA_M` (default 2)
secants are kept (M+1 pairs). Guards: residual-weighted Gram system,
non-affine cell exclusion (weighted residual/secant >1 keeps the raw
map output), Tikhonov ridge on the Gram diagonal, a garbage veto on
non-finite or `||corr||>1E3||f||` corrections, `10|r|` per-cell step
cap, positivity/`NONMONO_INEXT_CAP` fallback. If the pass after an
applied mix does not improve, the
pre-mix iterate is restored, history cleared and AA disabled for the
rest of that field solve. History collection starts at pass 4 so
fast rays never allocate the (M+1)*2*NI*NCF pair store; stored
residuals are rescaled when the adaptive omega changes the relaxed
map (exact for the affine map). Per-refresh `AA: applied= vetoed=
rejected= solve-fail=` counters print under the `ray exits` census.

Veto bound (13-Oct-2026): the inherited MOM `||corr||>||f||` veto is
mis-scaled for the marginal-mode regime this accelerator targets -
a cell contracting at rate lambda sits |r|/(1-lambda) from its fixed
point, so the correct AA step exceeds ||f|| by ~50-100x for the
measured lambda~0.98-0.99 stall rays. The strict veto rejected ~3x
more mixes than it applied on arm 12186580; the bound is now a
garbage filter at 1E3*||f||, with the per-cell cap and the rollback
remaining the real safeguards.
#status/unverified (built; runtime arm pending)

### Per-cell period-2 damping: `KN_ORD_P2DAMP` (default F)

Collapse of bounded lambda~-1 oscillations inside `SWEEP_FLD_ORD`.
The whole-ray `KN_SGN_SUM` detector added in 0a207b1 keys on the
signed residual sum over all significant cells, so a minority
oscillating set is drowned by the contracting bulk - on the m15_2p82d
kndump arm ~20 near-floor cells flip-flopped ~2x/~0.5x per pass for
60+ iterations at a stable 0.53 max residual, and `OMEGA_RAY` never
moved off 1.0.

With `KN_ORD_P2DAMP=T` each significant cell (above the 1e-4 column
floor) keeps the sign of its last update delta `(VNEW-VOLD)` across
passes. After `KN_ORD_P2N` (default 2) consecutive sign flips the
commit becomes the bracket midpoint - `sqrt(VOLD*VNEW)` when both
endpoints are strictly positive (the observed flip is
multiplicative, ratios ~0.475/2.1, so the geometric midpoint is the
exact fixed point of the pure mode), `0.5*(VOLD+VNEW)` otherwise -
until `KN_ORD_P2REL` (default 3) consecutive flip-free passes
release the cell back to the relaxed update. The fixed point is
unchanged: a midpoint inside the flip bracket is a bounded step of
the same fixed-point problem, not a clamp to a target. State is
three int*1 arrays (~3 B/cell) allocated per ray solve, so the
~80 B/cell `KN_RF_TOT_MB` thread budget is effectively unchanged.
The `nfm` slot does not need freezing: in the ordered march the
compressing upwind offset is computed inline in `KN_ORD_SPV` from
`DRF_L`/`KN_DLN_N` (static per cell), so the map is deterministic
pass to pass.

Diagnostics: per-ray `p2damped` write count accumulates into the
per-iteration `p2damped writes=` line (flag-gated, so the F build
prints nothing new), plus the worst-ray block. Validation arm:
`dessart_m15_2p82d_KNOP_FORMAL_p2damp` warm-restarts the it136
kndump checkpoint.

#status/unverified (built; runtime arm pending)

### Stall threshold: `KN_RAY_STALL` (default 0.98)

Same-direction residual ratio above which a pass counts toward the
`NSTALL>=2` early exit in `SWEEP_FLD`/`SWEEP_FLD_ORD`. The control
arm's exit census (m15 2.82d it1-3) showed ~12 stall rays hold
essentially all unconverged cells; raising toward 1 lets still-
contracting rays use more of `KN_RAY_MAXPASS` instead of quitting
at <10%/pass. The `KN_RAY_STALL=0.98` arm (13-Oct-2026) cut the
unconverged-cell fraction ~40x vs control at matched iterations
(0.21-0.29% vs 11-13%) and cleared its two `KN_RAY_MAXPASS` cap
hits on the it2 restart, so 0.98 is now the default. Set 0.90 in
VADAT to restore the pre-change behaviour.

**Ray-grid rebuild triggers** (`NEW_R_GRID`, evaluated on each INIT):
relative R drift >1E-8, `VDOP_FRAC`/`SOLUTION_OPTIONS` change,
compression-set membership change (with the -0.98/+1km/s exit
hysteresis), and - since 05-Oct-2026 - per-cell V and `SIGMA_TRUE`
drift >1E-6 relative to the `V_EXT`/`SIGMA_EXT` snapshot. The V/sigma
fingerprints close a stale-field hole: a velocity update at fixed R
(e.g. `DO_V_PARAMS`) changed neither R nor necessarily the
compression set, so `GAM`, `chi_hat` and the ray insertion spacing
were silently computed on the old field. Same check `V_SM_SV` does
in `MOM_J_CMF_NONMONO_V1`.

Rebuild must free **every** persistent ray array before the re-ALLOCATE
or the run stops on STAT=151 (already allocated). `NM_RAY`/`DLN_RAY`
(the GS `I_NEXT` slot-match store from the unconditional allocation
block) were missing from the dealloc list, so the first rebuild after
they were allocated crashed. Surfaced 13-Oct-2026 on X_SOB OBS: a
compression-set membership change at depth 51 triggered a rebuild and
`cmf_flux` died inside the first sweep. Fixed by deallocating both in
the `NEW_R_GRID` block; their contents index the old ray grid and are
invalid on rebuild anyway. #status/verified (build; X_SOB OBS resubmit
is the runtime test)

### Ordered march: `KN_RAY_ORDERED` (default T)

The legacy refresh did alternating **global** frequency passes -
half of all strong links were stale every pass, pinning the residual
on marginal `|Q|/(chi+|Q|) -> 1` links.

`SWEEP_FLD_ORD` instead marches **rows along the characteristic** and
frequency inside each row in the row's upwind direction (sign of
`GAM`): expanding rows march l ascending so the same-row
`FLD(j,l-1)` read is fresh; compressing rows march l descending. The
only remaining lag is the weak spatial corner - the Knop
**C-corner** coefficient C+/- `I_{i-1,l+1}` handled through the
four-deep `SPW` ring previewed by `KN_ORD_SPW_FILL`. Theory-derived,
matches Paper II Sec 3.2's "physically guided index order".

> [!quote] Knop et al. 2009, A&A 501, 813 (Paper II), Sec 3.2 - on the linear-source ordered march
> "the solution is then direct and its speed is optimal"

Analogy for the ordering: the alternating global passes are like
painting a wall top-to-bottom and bottom-to-top in turns - each pass
re-wets half the strokes. The ordered march follows the direction the
paint actually runs: each stroke only touches still-wet edges
(fresh neighbours) except at one weak corner, so the wall finishes in
a few passes instead of stalling on the last drips.

### Frequency matching: `DO_NONMONO_FREQ_MATCH` (default F)

Remaps the frequency grid along the ray so the upwind sample lands on
the characteristic foot - the Ceulemans et al. (2025) fix for the
oscillatory eigenmodes that survive any iteration ordering when the
per-cell Doppler drift exceeds the bin spacing. With |dV| up to ~330
bins/cell on m15, recommended ON for compressing models
(`MATCH_FREQ`, `NNF2 = ceil(DRF_L/KN_DLN_N)` capped on-grid).

## The MOM path (`USE_LAM_ES=F`)

`MOM_J_CMF_NONMONO_V1` - **not** Knop's characteristic solver; it is
the V11 variable-Eddington moment solve generalized to sign-changing
frequency coupling:

- Per-equation coupling: `a_J = CON_GAM*(1+sigma_true*K_ON_J)` at
  depths, `a_H = CON_GAMH*(1+sigma_true*NMID_ON_HMID)` at midpoints.
- `chi_hat(mid) = chi_avg + |a_H|/DLN`; `DTAU` uses chi_hat - no
  negative optical depths.
- `PSI = DTAUONQ*|Q|` (implicit diagonal) vs `PSIPREV` (explicit
  upwind). Key fix: compressing branches must use the **same local
  |Q| weight** on both sides (PSI_NM), making the RHS a convex
  combination - the lagged map is then non-expansive. The earlier
  signed-PSIPREV form could inject negative upwind sources after a
  `1+sigma*K_on_J` sign flip.
- Negative-J repair (`XM_FB`): same-frequency lagged `J_NEXT` where
  positive, else local `SOURCE*R^2` - **never** the bluer
  `RSQJNU_PREV` (that cross-frequency fallback caused the it1
  L=2.7e20 ratchet).
- `MOM_KNOP_REFRESH_V1` re-drives the solver itself in alternating
  frequency order; `JNEXT_RF_READ` pass-start snapshot keeps the
  compressing reads Jacobi (matching in-sweep semantics - a live
  descending read created a GS cascade that exploded to 1E13).
- Containment: `1E40` XM cap, `JNEXT_BEST` best-pass restore, abort at
  `PASS_CHG > 1E3*BEST_CHG`, quarantine-age counters
  (`MOM_GS_NPQUAR` at age 3+).

> [!warning] MOM path is a moment analogue, not the paper operator
> Its lagged map is genuinely non-contractive at sigma ~ -37..-55
> compressions (~3% of cells sign-flip every pass, in-sweep too).
> OBS requires `USE_LAM_ES=T` anyway; production runs use the FG/ray
> path. See [[Open issues]].

## Checkpoints

`FGINEXT_GS` (FG store) and `MOMJNXT_GS` (MOM store) written by
`DUMP_*_GS` inside `SCR_RITE_V2` next to POINT1/2; fingerprinted on
restore (bitwise R/V/SIGMA + `COMP_DEPTH_IDX` match). See
[[Checkpoints and state files]].

## Performance notes

- Refresh ~640-750 s of a ~930-1000 s great iteration; memory/
  latency-bound, threads saturated - `cpus_per_task=16` adopted as
  production default (R1 in [[future_fixes]]).
- Remaining candidates: MOM-path coefficient caching (R2), strided
  store access (R3), single-precision chi/eta coefficient pack (R4),
  sweep OMP tail (R5).

## Related

- [[The non-monotonic problem]]
- [[Field maturity and GS freeze]] - what the refresh feeds
- [[JFNK solver]] - the sweep is re-run per eval in coupled mode
- [[Literature]] - Knop I/II, Ceulemans 2025
