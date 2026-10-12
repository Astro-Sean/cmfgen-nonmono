# JFNK review — m15 2.82d SOB probe (jobs 12155376 / 12155427)

Each section follows: hypothesis, evidence, change, new VADAT keys and
defaults, open questions.  Log line references are to the probe OUTGEN
`.../m15_2p82d_SOB_jfnk_probe/cmfgen/2.82days/OUTGEN` unless noted.

## Issue 1 — nits=0 with no reason

- Hypothesis: an unlabeled early return between "JFNK SETUP" and the
  nits print.
- Evidence: the Newton loop in `newton_hook_v1` (jfnk_solver_v1.f90)
  could exit at the entry convergence test, the max-iteration cap,
  GMRES stagnation, or the reach gate, and none tagged a reason.
  COLPROBE/REACH-FD restore state (repeat-vs-pristine |diff|=0
  verified on probe OUTGEN).
- Change: `jfnk_exit_reason` (module var, `jfnk_solver_v1_mod`) is set
  at every solver return path and printed in the JFNK RESULT block and
  at the bail-outs in `jfnk_solve_corr_v1.f`.  Commit 5d68bc0.
- Probe evidence (12155427, it444): exit path now labeled — GMRES
  stagnation -> "finite stagnation; retaining caller mode".

## Issue 2 — residual scaling / convergence metric

- Hypothesis: Euclidean |F| is dominated by converged HI rows; the
  metric should measure relative corrections.
- Evidence: top-10 residual rows are HI n=22/23 at K=65-70, ~100%
  collisional; F_i/BA_ii ~ 1e-6 for those rows (they are done).  The
  anomalous fix: divide each row by its BA diagonal -> residual in
  log-population-correction units.
- Change: new RSCALE mode `DIAG` (`JFNK_RSCALE: DIAG`, default still
  EUCLIDEAN) — pop rows scaled 1/max(|BA_ii|,1e-30*blockmax); the Ne
  (charge) row scaled by 1/ED(K); T row by its BA diagonal.  New
  `jfnk_maxch_est` = max_i |F_i|/|BA_ii| over population rows,
  computed on bare F inside JFNK_RESIDUAL_V1 before the PTC term is
  added; printed on the `nits=` line and in the RESULT block.  New key
  `JFNK_MXCH_TOL` (default 0 = off): adds a second convergence test so
  a stale residual norm cannot stall a finished solve.  Commit 7a4ef65.
- GROSS mode: not implemented — no stored gross-rate vector exists;
  |BA_ii| for a rate-equation row is the total outflow scale, so DIAG
  is a close surrogate.

## Issue 3 — clip before reach test / globalization

- Evidence: REACH/REACH-FD priced the raw M^-1 F direction whose
  |z|_inf was ~800 dex; the commit path then clipped each component to
  +/-ln(MAX_JFNK_COR).  Diagnostics described a step that could never
  be taken.
- Change: in the reach block, `Z_CLIP` = Z_CHK with the same
  per-component bounds (pop/Ne: X_MIN..X_MAX; T rows: X_T bounds,
  honoring MOD_FIX_T_D range).  Stored-BA reach and the FD pricing now
  run on Z_CLIP; the print `REACH clip: max|z| raw -> clipped, nclip=`
  reports the reduction.  The unclipped variant is kept as the
  (|z|>1 zeroed) diagnostic.  Trust-region cap inside the solver
  (jfnk_step_inf_cap) already limits |s|_inf uniformly.  Commit b4a6fa2.
- Probe (pending): rerun with the new keys will show reach priced on
  the clipped step.

## Issue 4 — huge-|z| unknowns

- Evidence: top-|z| components are high-J rows (J=660 at K=43, J=307
  K=34, J=298/290 K=25, J=461/462 K=26, J=1832 K=31, J=304 K=49,
  J=1255 K=21, J=743 K=32) — high-lying levels of mid/high ions, not
  the residual-carrying HI levels.  Species census print added in
  c1a78dd labels ion/level/n/nLTE/BA-diag for the top 50.
- Interpretation: near-empty, radiation-dominated levels (FeSEV,
  ArSIX, CaSIX show rad/F~1 in the anatomy) with near-null BA rows;
  the LU solve returns multi-dex corrections the commit clip flattens.
- Change (issue 4a): `JFNK_POP_FL` (default 0 = off).  At factor time,
  rows with POPS_BASE < JFNK_POP_FL * (ion total at that depth) are
  marked ROW_REG: identity row in the factored block, z_i projected to
  zero at apply.  Newton cannot move them; the residual row itself is
  unchanged.  Commit cd8ac5e.
- Acceptance check: rerun with JFNK_POP_FL ~1D-6 should show max|z|
  drop from ~800 to O(clip bound) and PC residual improve (the
  near-null rows are also the conditioning outliers).

## Issue 5 — per-block ridge logic

- Hypothesis: ridged=70 is a trigger artifact.  REFUTED: the probe
  YAML sets JFNK_RIDGE=1D-4, which becomes jfnk_ridge_floor — every
  block starts ITRY=0 with ridge 1e-4 by design.  The ridge is
  per-row proportional (1e-4*|diag_i|) so it cannot swamp weak rows
  and is applied to charge/conservation rows identically (relative).
- Change: per-depth diagnostics — min/max pivot of the equilibrated
  factorization, DGECON 1-norm rcond, and the ITRY index that produced
  it, printed beside each ridged depth.  New key `JFNK_RIDGE_PIV`
  (default 0 = off): a successful-but-ill-conditioned first
  factorization (pivot ratio below the key) is rejected so the ITRY
  loop escalates the ridge.  Commits 9d0580e, bf2b646 (alloc fix).

## Issue 6 — BA vs FD Jacobian mismatch (topF column, HI J=23 K=70)

- Evidence (12155427 it444): |FDcol|=6.48e18 |BAcol|=7.41e18,
  |diff|/|BA|=14.6%.  Diagonal J=23 agrees to 13% (FD -6.04e18 vs BA
  -6.97e18); J=22 agrees to 6%.  The downward couplings flip sign:
  J=21 FD -4.76e17 vs BA +5.44e15 (~87x, sign flip), J=25 -2.98e16 vs
  +1.33e15, J=20 -2.06e16 vs +1.77e15.
- Interpretation: the FD residual differentiates through the
  dissolution/occupation-probability factor w(n); the stored BA column
  was assembled with w held at the activation baseline
  ("POPION/dissolution factors held at activation baseline").  The
  sign-flipped downward entries are consistent with dw/dn terms BA
  does not carry.  topZ column (Fe/Ne charge) agrees to 1.3e-3 — the
  mismatch is concentrated in the high-n hydrogen shell structure.
- Change: h-scan added on the topF self element (COLPROBE hscan topF),
  matching the existing topZ scan (commit b4a6fa2).  Whether the
  mismatch is real physics (w-derivatives) or cancellation noise will
  show as h-stability in the next probe.  cos(BAz,A_FDz) for the
  clipped step now prints under Issue 3.

## Issue 7 — cos(BAz,A_FDz)=1.85e11

- Cause: JJZ_DP (|BAz|^2) was overwritten inside the colprobe block by
  the repeat-vs-pristine norm before the second REACH-FD print used it
  as a denominator.
- Fix: dedicated scratch variable; JJZ_DP retains |BAz|^2.  Regression
  checker `cmfgen_jfnk_unit_test/verify_cosine_bounds.py` parses OUTGEN
  for |cos|>1.  Commit fd3d155.

## Issue 8 — admission gate

- Machinery already exists and is correct: the gate uses BA_LAST_RAW
  (last non-LAMBDA demand), a contraction test
  (JFNK_DEM_IN <= JFNK_DEM_RISE * last), and JFNK_GATE_NREQ hysteresis
  on consecutive eligible calls.  DEFER/ADMIT blocks print the
  deciding values.  The "gate = 1e30" in the log is the YAML value
  JFNK_FRZ_GATE=1.0D30 (always admit), not a code default.
- Recommendation (config, not code): JFNK_FRZ_GATE ~3.0D1 (admit once
  demanded correction stays below ~30%), JFNK_GATE_NREQ=2..3,
  JFNK_COR_GATE ~1.0D1 for the coupled tier when it exists.
- Luminosity diagnostic: "Luminosity of star = 0" under FIXED_T +
  IB_METH=ZERO_FLUX is meaningless — RLUMST is not integrated in this
  configuration.  Left unchanged (upstream print; minimal diff) — treat
  the line as dead output in this mode.

## Issue 9 — prepare_VADAT.py messaging

- Written VADAT (fact): VEL_LAW=14, RD_IN_R_GRID=F, VFRAC_MOM=1.0,
  USE_LAM_ES=T, DC_METH=R, RSTAR=1.454688D4, RMAX=8.614.
- Sources of contradictions:
  - VEL_LAW/RD_IN_R_GRID: the profile-expectation print at ~line 770
    says "VEL_LAW=11, RD_IN_R_GRID=T" before the vel-law selection
    branch runs; effective_vel_law=14 then writes 14/F correctly with
    its own print at ~line 5083.
  - VFRAC_MOM: the optimizer prints 5.00 first; a later branch
    (target-nd / vel-law path) writes 1.0.
  - USE_LAM_ES: printed as preserved/F in an early message but the
    non-mono profile sets T and the file shows T.
  - DC_METH: early print claims "LTE preserved"; the non-monotonic-T
    enforcement (line ~5560) rewrites the line to DC_METH=R — the file
    is correct, the message was not.
  - RSTAR/RMAX: vadat_overrides replaced the tau=2/3 photosphere
    values (41342->14546, 3.03->8.61) with no warning.
- Change: a "[Final configuration]" block now parses the written VADAT
  back and prints effective values for the physics keys, flags
  override-vs-hydro scalar mismatches, and warns when final
  RSTAR/RMAX differ from the photosphere anchor.  Outer repo commit
  964c2f92.

## Issue 10 — physics/data sanity (report only)

- Velocity field: RDINR is monotone decreasing outward
  (48.9 km/s at the surface rising inward to 7074 km/s) — a
  decelerating profile; sigma = dlnV/dlnR - 1 spans [-120, -1.05].
  18 jumps > 20 km/s confirmed, concentrated at the inner edge
  (R ~ 1.46e14 cm) where dV exceeds 500 km/s per cell.
- Sobolev validity: at the steep inner jumps sigma << -1 gives a
  short Sobolev length (fine).  In the outer near-flat cells
  (sigma ~ -1, V ~ 50-135 km/s, VTURB=50) dv/dr -> 0 makes the
  Sobolev length long; line transfer there is only approximate — the
  DAH tau floor bounds |1+sigma*mu^2| but does not restore physical
  accuracy.  This is a known limitation of the non-monotonic Sobolev
  treatment, not a bug; flag when interpreting outer-region line
  formation.
- SIV_F_OSCDAT line 1690: 3s_3p(3Po)4s_2Po[1/2]-3s2_5s_2Se[1/2]
  (levels 45-46) has f = 1.100E+01 — unphysically large for a single
  transition (f should be O(1) or less); GENOSC_V9's "bad oscillator
  strengths" warning is correct.  Recommend replacing with a NIST/
  TOPbase-consistent value upstream; data file left untouched.
- HeI_F_OSCDAT levels 23/24 (1s5d3De 193917.2568 vs 1s5d1De
  193917.2546 cm^-1) differ by 0.0022 cm^-1 — effectively degenerate;
  the out-of-order warning is benign.  Documented, not edited
  (read-only atomic data).

## GMRES restart behavior (this branch)

- Lucky-breakdown now ends the Arnoldi cycle (flag `lb_hit`), the LSQ
  step is formed, and the true residual re-seeds v(:,1).  Restart cap
  prevents runaway cycling.
- Live evidence (12155427 it444): cycle 1 claimed hess_res 5.67e13 at
  it=21 (a "37% reduction") but RESTART DIAG showed true_res=3.43e14
  (ratio 6.06 — the Hessenberg residual underreports because the FD
  JVP is not a fixed operator across evals).  Cycle 2 died at it=22
  with hess/true in agreement (ratio 1.03); restart cap -> stagnation
  -> clean Lambda fallback.
- The 16% preconditioner residual (||Mz-v||/||v|| on v=F) is genuine
  LU backward error on blocks whose equilibrated conditioning is
  ~1e15+; it is not an audit artifact (raw-vs-factored multiply was
  ruled out — the decomposed audit reports 4.7e12 on pop rows while
  T/Ne rows carry |v| only 4.7e-2).  Reducing row-space conditioning
  via JFNK_POP_FL / higher JFNK_ROW_TINY is the lever, not more ridge.

## Before/after log excerpts

Before (failed job 12153985, OUTGEN ~line 1740, pre-branch binary):

```
 weighted residual norm at entry   :  4.1172E+12
 JFNK REACH: frac=  9.9996E-01 |Jz|/|F|=  1.0001E+00 cos(F,Jz)=  9.9998E-01
 JFNK REACH-FD: |A_FD z|/|F|=  4.6573E+00 cos(F,A_FDz)=  2.4551E-01 cos(BAz,A_FDz)=  2.4551E-01
 JFNK REACH-FD (|z|<1): |A_FD z|/|F|=  8.7615E-01 cos(F,A_FDz)=  7.4728E-01 cos(BAz,A_FDz)=  2.7485E+10   <- impossible cosine
 JFNK: equilibrated block LU   nfact=   1  ok=  70  ridged=  70  Jacobi=   0    <- no pivot/rcond info
 JFNK: nits=  0 res=  4.1172E+12                                                             <- no exit reason
 GMRES: it=  1 res=  3.9912E+12 ... it= 4 (then silence — no restart lines)
 JFNK: GMRES stagnated — aborting Newton (Lambda fallback)
```

After (probe 12155427, OUTGEN ~line 1030, restart cap + PTC fold +
frozen-row reduction + exit reasons; predates Issues 2-6 code):

```
 weighted residual norm at entry   :  1.0334E+14
 JFNK REACH: frac=  9.9759E-01 |Jz|/|F|=  1.0053E+00 cos(F,Jz)=  9.9882E-01
 JFNK REACH-FD: |A_FD z|/|F|=  5.5564E+00 cos(F,A_FDz)=  2.5018E-01 cos(BAz,A_FDz)=  2.5835E-01
 JFNK REACH-FD (|z|<1): ... cos(BAz,A_FDz)=  1.8515E+11   <- JJZ_DP bug still present in this binary
 GMRES RESTART DIAG: it= 21 hess_res=  5.6687E+13 true_res=  3.4338E+14 ratio=  6.0576E+00
 GMRES RESTART DIAG: it= 22 hess_res=  3.4024E+14 true_res=  3.5140E+14 ratio=  1.0328E+00
 GMRES: restart cap (  2) -- treating as stagnation
 JFNK: finite stagnation; retaining caller mode, count   1 of   3    <- labeled exit
 second admission: JFNK: GMRES stagnated; using inexact direction (pred_res/ref=  8.7030E-01)  <- no longer a blind abort
```

## Test status

- `cmfgen_jfnk_unit_test` surrogate suite: all tests pass except the
  four pre-existing failures documented in AGENTS_HISTORY (T3 'mild'
  demand threshold; T4 core-LU residual on mild/hard/shock).  New
  `verify_cosine_bounds.py` correctly flags the three out-of-bounds
  cosines in the pre-fix probe OUTGEN.
- `unit_tests` production-file suite: 10/11 files pass;
  test_interp_guards fails on uninitialized-memory jitter in output
  column 3 (routines untouched since baseline a15dba1, flaky pass set
  varies run to run; the interpolated column is bit-identical).
- Build: `make DEBUG=nodebug all` clean, 0 errors;
  exe/cmfgen_dev.exe, cmf_flux.exe, obs_frame.exe all produced.

## Recommended next probe VADAT keys (non-default, opt-in)

- JFNK_RSCALE: DIAG        # relative-correction metric
- JFNK_POP_FL: 1.0D-6      # freeze near-empty levels in Newton step
- JFNK_MXCH_TOL: 5.0D-2    # converge when max|F/D| < 5%
- JFNK_RIDGE_PIV: 1.0D-10  # escalate ridge on measured ill-conditioning
- JFNK_FRZ_GATE: 3.0D1, JFNK_GATE_NREQ: '3'  # sane admission defaults
