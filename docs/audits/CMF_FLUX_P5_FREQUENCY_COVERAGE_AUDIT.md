# CMF_FLUX P5 audit note - observer-frame frequency coverage

Scope: `obs/obs_frame_sub_v9.f` window construction, Doppler extrema,
frequency hunts, edge clamps and census. Written before the P5
implementation; each item states whether it is a confirmed defect.

## 1. The frequency mapping (verified, not a defect)

In-loop, per ray point I and observer frequency ML:

  `T1 = OBS_FREQ(ML) * (1 - VMU_RAY(I)) * GAM_RAY(I)`   (= nu_cmf)

with `VMU_RAY = V_RAY * Z_RAY / R_RAY / C_KMS` (v dot n / c, positive
when the flow moves toward the observer) and
`NU_ON_NU0_RAY = 1/GAM/(1-VMU)` (nu_obs/nu_cmf). Consistency:
`nu_cmf = nu_obs*(1-VMU)*GAM` rearranges to `nu_obs = nu_cmf*NU_ON_NU0`
exactly as stored. For an outward-moving front surface (VMU>0) the CMF
frequency is red-shifted; on the mirrored far side `VMU -> -VMU` gives
the blue-shifted partner. This is the standard invariant form
`eta/nu^2`, `chi*nu`, `I/nu^3` (Mihalas 1978, Sec. on observer/CMF;
Rybicki & Lightman 1979 Sec. 4.9) and matches the powers already
verified in the main audit (emissivity *NU_ON_NU0^2, opacity
/NU_ON_NU0, intensity *NU_ON_NU0^3). Non-rel limit: GAM=1,
NU_ON_NU0=1 - identical formula, no branch asymmetry. Sign convention
matches the upstream CMFGEN convention (outflow toward observer =
blueshift); COMP_OBS_V3 uses the same sign sense on the boundary
intensities.

Mirrored ray points are written into VMU_RAY/GAM_RAY at :620-634
BEFORE the extrema at :705-729, so both ray sides enter every MAX/MIN.
NRAY is NR (no mirror) for core-intersecting non-hollow rays - there is
no uninitialized tail read (:554-556).

## 2. The coverage theorem (verified, not a defect)

For band OUT_ML covering observer indices [m1,m2], the required set is

  R = { OBS_FREQ(ML) * X(I) : ML in [m1,m2], I in [1,NRAY] },
  X(I) = (1 - VMU_RAY(I)) * GAM_RAY(I) > 0.

Since OBS_FREQ is strictly decreasing (enforced upstream in
cmf_flux_sub_v5.f:1262) and X>0:

  min R = OBS_FREQ(m2)*X_min,   max R = OBS_FREQ(m1)*X_max.

The construction computes the bound B = OBS_FREQ(m1)*T2_hi with
T2_hi = 1+1.2*MAX_VMU (non-rel, 20% slack over X_max = 1+max(-VMU)) or
T2_hi = X_max exactly (rel), then hunts J down the decreasing
FREQ_CMF until B >= FREQ_CMF(J). The previous step held
B < FREQ_CMF(J-1), so

  SM_FREQ_CMF(1) = FREQ_CMF(max(1,J-10)) >= FREQ_CMF(J-1) > B >= max R.

The window therefore STRICTLY contains every required T1 whenever the
global grid contains the bound B - even when STRT is pinned at 1,
because then SM_FREQ_CMF(1)=FREQ_CMF(1) >= FREQ_CMF(J-1) > B still
holds. The red end is symmetric.

Consequence: local-window undercoverage is unreachable under the
current construction. Every edge clamp that fires corresponds to the
bound B or R lying outside FREQ_CMF - i.e. GLOBAL undercoverage only.
The +/-10-cell padding only needs to be >=2 for bracket safety; 10 is
generous. NUM_BANDS partitions only the observer index range; each
band still carries the same ray-global extrema, so band count cannot
change coverage, only the size of each SM_FREQ_CMF copy.

What the hunts do NOT protect: the bound itself exceeding FREQ_CMF.
The bounded hunts saturate at the grid edge (:720-737), the window is
then pinned at the global edge, and per-point T1 values beyond it hit
the edge clamps. That is the only reachable undercoverage class.

## 3. Confirmed defects

### D5.1 - 2% edge tolerance is a physical extrapolation, not a tolerance
Sites: :846-868 (loc-1 blue), :883-899 (loc-1 red), :951-967 (loc-2
blue), :997-1013 (loc-2b red).
Each allows T1 up to 2% outside the local window and substitutes the
window-edge opacity/emissivity. Measured on a real non-monotonic model
the overshoot was ~0.1% (dev7 comment, :874-877) - a genuine global
coverage gap, substituted silently. The tolerance was never derived;
the correct round-off scale is set by association-order differences in
`OBS_FREQ(ML)*(1-VMU)*GAM` vs `OBS_FREQ(m1)*MAXVAL((1-VMU)*GAM)`, i.e.
O(eps_mach) ~ 1e-15 relative. Correction: tolerance -> 1e-12 relative
(documented: ~10^3 ulps of headroom, orders below the smallest grid
spacing or physical Doppler gap); anything larger aborts.

### D5.2 - local/global undercoverage conflated in the census
N_EDGE_CLAMP counts all four sites together; CLAMP_GRID_MAX tracks
global overshoot only at the single worst point. Given section 2 the
distinction is moot today (all clamps are global) but the census
should still classify each clamp against FREQ_CMF, so that a future
construction regression reports LOCAL_WINDOW_UNDERCOVERAGE rather than
looking like a grid gap. Correction: split counters N_CLAMP_LOCAL /
N_CLAMP_GLOBAL, keep per-class worst-case (ML,LS,I,T1,edges).

### D5.3 - no bracket-invariant validation
At :908 and :1023 the interpolation fraction
  `T2 = (T1-SM_FREQ_CMF(K))/(SM_FREQ_CMF(K-1)-SM_FREQ_CMF(K))`
is used with no check that K in [2,SM_NCF], the denominator is
positive, T1 is finite, or T2 in [0,1]. A NaN T1 (NaN in V) passes
both hunts silently (every comparison is false), produces NaN ETA/CHI,
and is only caught much later by CHK_OBS_FLUX. Correction: a compact
guard after each bracket hunt - T1 finite, K in [2,SM_NCF],
denominator > 0, T2 in [-eps,1+eps] with eps=1e-12 - aborts with full
context. ETA/CHI finiteness checked right after interpolation.

### D5.4 - no serial verification of the construction
The proof in section 2 relies on the hunts having placed the bounds
inside the grid. There is no runtime assertion of this. Correction:
after SM_FREQ_CMF is loaded for each band, verify
  `OBS_FREQ(m1)*X_max <= SM_FREQ_CMF(1)` and
  `OBS_FREQ(m2)*X_min >= SM_FREQ_CMF(SM_NCF)` within eps;
violations are LOCAL_WINDOW_UNDERCOVERAGE (construction bug) or
GLOBAL_CMF_GRID_UNDERCOVERAGE (bound outside FREQ_CMF) - report and
stop. O(1) per band per ray.

### D5.5 - global undercoverage currently continues with warning only
Per spec and per scientific honesty: when the required frequency is
outside FREQ_CMF the opacity simply does not exist; substituting the
edge value produces a plausible but wrong spectrum. Correction: abort
with GLOBAL_CMF_GRID_UNDERCOVERAGE naming ML, LS, I, T1, window and
grid edges, and the fix (widen MIN_CF/MAX_CF). Behavioural change on
the dev7 clamp path only; upstream never had these clamps (it read
out of bounds). The abort threshold is the derived 1e-12 - real gaps
(~1e-3 measured) abort, association-order ulps (~1e-15) still clamp
and count.

### Not defects (verified)
- Hunt directions: both walks re-bracket from any stale INDX start;
  velocity reversals are handled correctly, at worst extra iterations.
- Mirror symmetry: extrema include the far side; non-hollow core rays
  correctly truncate at the inner boundary.
- Rel branch exact extrema: MAXVAL((1-VMU)*GAM) over the discrete ray
  set is exactly the maximum required factor - coverage needs no slack.
- Band partitioning: each band is independently fully covered; NUM_BANDS
  changes memory use, not the spectrum.
- INDX reuse: stale start index only costs hunt iterations; the two
  while loops re-establish the bracket invariant unconditionally.

## 4. Proposed minimal change set

1. Scalars: `EDGE_TOL=1.0E-12_LDP`, extrema `X_MIN_RAY`,`X_MAX_RAY`
   hoisted out of the OUT_ML loop (loop-invariant today, recomputed per
   band); counters `N_CLAMP_LOCAL`,`N_CLAMP_GLOBAL` + per-class worst
   records.
2. Per-band serial coverage assertion after SM_FREQ_CMF load.
3. Four clamp sites rewritten to the derived tolerance + classification
   + named aborts (LOCAL_WINDOW_UNDERCOVERAGE /
   GLOBAL_CMF_GRID_UNDERCOVERAGE), same comment-block style.
4. Bracket guard after each hunt (two sites).
5. Census report extended with the per-class counts.
6. Tests: drv_obs_frame_cov.f - monotonic/reversing/negative V, global
   undercoverage abort, NUM_BANDS partitioning (NOS>=2000), rel mapping
   and non-rel limit, grid-expansion convergence (covered -> identical
   spectrum), OMP matrix. Static pins extended.
