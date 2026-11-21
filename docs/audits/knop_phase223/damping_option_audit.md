# Phase 2.23 source audit: KN_ORD_P2DAMP and FG_ORD_AA
All line refs: plane/fg_j_cmf_nonmono_v1.f (current working source).

## KN_ORD_P2DAMP (VADAT [KN_ORD_P2DAMP], default F)
- Implementation: inside SWEEP_FLD_ORD commit block, ~12254-12295
  (P-channel MARCH cells) and ~12350-12384 (inner boundary row);
  state decl ~11973; activation ~12052-12061.
- State tested: SIGN of (VNEW - VOLD), the per-cell per-PASS update
  delta inside ONE ray/channel solve. Not refresh-to-refresh.
- Granularity: per cell (J,L), per channel solve. int*1 arrays
  P2_SGN/P2_FLIP/P2_DMP allocated per solve (3 B/cell, per thread).
- Detection: SGN_P != P2_SGN on FG_ORD_P2N=2 CONSECUTIVE flips ->
  P2_DMP=1. Requires the cell to pass the significance gate
  (MAX(VOLD,|VNEW|) > 1E-4*VMAXF) at ~12224 - floor cells below the
  gate are NOT tracked.
- Action while damped: FLD=SQRT(VOLD*VNEW) when both endpoints >0
  (exact fixed point of a pure multiplicative lambda=-1 mode);
  else FLD=0.5*(VOLD+VNEW). Bounded inside the bracket. Replaces the
  relaxed commit VOLD+OMEGA*(VNEW-VOLD). Release after FG_ORD_P2REL=3
  consecutive flip-free passes.
- What it modifies: the ITERATE PATH inside one solve. It does not
  change the fixed-point map T(x); the committed midpoint is inside
  the measured bracket. RES still measures the raw |VNEW-VOLD| defect
  before the commit -> no false-convergence pathway through RES.
- False-positive behaviour: a converging cell keeps a constant delta
  sign -> never damps. A noisy cell at the significance boundary can
  damp but the commit stays inside the bracket -> bounded harm, and
  the residual metric stays honest.
- KN_CERTIFY interaction: none direct - certification reads RES_LS/
  IXIT which are computed from the same march metrics; if the damped
  state sits at the fixed point the ray converges legitimately.

## FG_ORD_AA (VADAT [KN_ORD_AA], default F; KN_ORD_AA_M default 2)
- Implementation: pass-end block ~12907-13030 inside SWEEP_FLD_ORD;
  rollback ~12603-12640; history rescale on omega change ~12753.
- Operation: type-II Anderson on the pass iterate. Secant pairs
  (AA_XH iterates, AA_RH relaxed residuals) depth M+1=3; weighted
  Gram system (weight 1/MAX(|x|,1E-4*colmax), same normalization as
  the march metric), Tikhonov ridge 1E-8, DGESV solve; applied
  correction per cell FLD = x - sum gamma*(x - x_j - r_j).
- Safeguards: non-affine cell exclusion (weighted res/secant >
  AA_GCAP=1.0 excluded and keeps raw map output); veto non-finite or
  ||corr||>1E3*||f|| (KN_AA_VETO); per-cell step cap 10*|r|
  (NONMONO_INEXT_CAP bound too); positivity fallback (XV<=0 keeps T4);
  singular Gram -> KN_AA_SNG, skip; residual-acceptance rollback
  restores the pre-mix iterate and disables AA for that solve
  (KN_AA_REJ) when the next pass does not improve on AA_REF.
- Supersedes Aitken: DO_AIT = FG_ORD_AITKEN .AND. .NOT. FG_ORD_AA.
- Scope: all rays/all cells within a solve; modifies the iterate path,
  not the map. Same honesty argument as P2DAMP for RES.
- KN_CERTIFY: none direct.

## Applicability to the measured stall (evidence-based)
- The stall argmax pool on rays 81P/82P/83P: ~6K-80K floor-dominated
  cells above KN_RAY_TOL at j=2 (P channel), dabs <= ~4E-9 absolute.
- P2DAMP applies only to cells inside the significance gate; argmax
  cells with |VNEW| > 1E-4*VMAXF are tracked. Whether the stall cells
  flip per-pass is UNPROVEN - needs the KN_SGN_SUM/per-cell sign
  history (proposed driver print, not yet implemented).
- AA accelerates marginal (rho~0.99) modes; the measured rho_max
  0.9895 is in its target regime, but floor-dominated cells may fail
  the AA_GCAP affine test and keep the raw map output.
- Neither option changes the floored-denominator RES metric. If the
  working diagnosis holds (metric-level failure on dark cells),
  the more direct test is an absolute-denominator residual variant
  (user Task-4; a solver change - not implemented here).
