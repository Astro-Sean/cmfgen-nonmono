# VEL_LAW=14 adversarial correctness review (RV_SN_MODEL_SNHYD, TAU_OPT=1)

Date: 27-Nov-2026. Auditor: Devin. Trigger: checklists/velaw14.md,
owner concern verbatim: "I am concerned that vel_law 14 is not
correct." Baseline: branch jfnk-acceptance-audit @ 0bf6a89d (HEAD;
the checklist's clone-from-GitHub rule maps to this SHA - origin/main
is 621 commits behind). Supersedes
docs/audits/VEL_LAW14_SNHYD_GRID_AUDIT.md (22-Nov-2026).

## Executive summary (10 lines)

1. VERDICT: VEL_LAW=14's grid/derivative construction is CORRECT -
   no physics or numerics defect found in the audited path.
2. The Steffen cubic matches the 1990 paper line-for-line
   (interior D = secant-weighted mean, sign-sum kills D at extrema);
   SIGMA is the analytic derivative of the SAME cubic used for V.
3. Measured: sigma converges at order ~2.0, V at order ~3.0 on an
   analytic monotone profile; power-law sigma err 5.7E-4 at ND=70.
4. The four pinned input defects (D1-D4) are fixed test-first:
   IEEE finiteness + positivity screens + missing-Sigma warning.
5. New evidence: ltau-shift invariant to 1E-14; homology rescaling
   exact; interior extremum preserved to 6E-14 in V.
6. Measured CDS resolution vs V13-style uniform-log-tau: ~4.3x
   nodes (43 vs 10 below R=36460) - corrects the "5-8x" doc claim.
7. Documented limitation (not a defect): the expansion factor
   aborts for SN_AGE > file age + ~0.06 d on m15's sigma=-82 CDS.
8. VEL_SPACE regrid already aborts under nonmono flags; REV_RGRID=F
   keeps regridding unreachable on production decks anyway.
9. G items (live Slurm A/B arms) remain open pending owner
   coordination - the static audit is complete.
10. Recommendation: keep VEL_LAW=14 as pipeline default.

## 1. Call chain (file:line)

```
VADAT VEL_LAW=14
  rd_control_variables.f:250-252   SN_MODEL=T, SN_HYDRO_MODEL=T
cmfgen/new_main/cmfgen_sub.f
  :1098-1100  RV_SN_MODEL_SNHYD(R,V,SIGMA,RMAX,RP,SN_AGE_DAYS,ND,LU,
                               TAU_OPT=1)          <-- audited routine
  :1106-1111  VADAT [RSTAR]=RP, [RMAX]=RMAX/RP rewrite
  :1120-1121  VINF=MAXVAL(|V|) (nonmono flags; VINF=V(1) otherwise)
  :1128-1134  abort if any V<=0 under DO_NONMONO_SOB/KNOP
  :1143       RECOMPUTE_SIGMA_TRUE (unclamped sigma_true; SOB reclamp
              SIGMA>=-0.999; stock-mode nonmono warning otherwise)
  :1153       SET_ABUND_CLUMP -> :245 RD_SN_DATA (density/ED/T/mass
              fractions/isotopes onto new grid, VOL_EXP_FAC)
solve_for_pops.f
  :4117-4122  ADJUST_R_GRID_V4 (REVISE_R_GRID-gated) +
              RECOMPUTE_SIGMA_TRUE on the revised grid
  :4755       RECOMPUTE_SIGMA_TRUE after regrid
```

Off-path routines named in the checklist: `do_vel_regrid.f` (v1) is
dead code (no callers); `do_vel_regrid_v2.f` is reachable only via
ADJUST_R_GRID_V4 GRID_METH=VEL_SPACE, which aborts under
DO_NONMONO_SOB/KNOP (adjust_r_grid_v4.f:160-169) and is unreachable
when REV_RGRID=F (pipeline default). `do_cmf_hydro*.f` /
`cmf_sub_sonic_hydro.f` are DO_HYDRO-gated (cmfgen_sub.f:4735,660,
1001) - SN decks do not exercise them. `misc/rev_rvsig.f` is an
offline interactive tool (PROGRAM REVISE_RVSIG), not in the run path.

## 2. Physics contract

- Depth ordering: index 1 = outer (RMAX), ND = inner (RCORE).
  Input R screened strictly decreasing (snhyd:281-287); output R
  screened strictly decreasing (:483-489); packing fractions
  1/20,1/8,1/3 at both ends preserve order (:462-472). CONFIRMED.
- SIGMA = dlnV/dlnR - 1, evaluated as the analytic cubic derivative
  (:517): `SIGMA = ((3*C1*T+2*C2)*T+C3)*R/V - 1`. Consumers
  (Sobolev 1+sigma*mu^2, SIGMA_TRUE same definition) consistent.
  CONFIRMED.
- Homology: V=R/t gives SIGMA=0; code zeroes the whole array when
  MAXVAL(|SIGMA|)<1E-4 (:522). snhyd_hubble test returns exactly
  0.0d0. CONFIRMED (A5).
- Units: expansion T1 = 86400 s/day * 1E5 cm/km / 1E10 = 0.864
  (km/s*day -> 1E10 cm), :279; kappa cm^2/g * density g/cm^3 *
  1E10 -> 1/(1E10 cm) consistent with R in 1E10 cm (:337,346).
  RADIUS/V/DENSITY/KAPPA units match the file headers. CONFIRMED.

## 3. Findings by checklist section

### A. Provenance / orientation / definitions

- A1 CONFIRMED. Chain mapped above; no guard-bypassing path found.
  do_vel_regrid.f v1 dead; v2 gated; DO_HYDRO paths off.
- A2 CONFIRMED. Every monotonicity guard checked (input :281,
  expanded :308, output :483, all strictly-decreasing with
  I,R(I),R(I+1) in the abort message). Monotone-outward synthetic
  files produce the analytic SIGN of sigma (test battery: beta,
  power, reversal cases; diagnostics JSON power_law_b1.5).
- A3 CONFIRMED. Producer writes dlnV/dlnR-1 (:517); consumers:
  Sobolev routines use 1+sigma*mu^2 (sobjbar_sim.f et al.),
  SIGMA_TRUE recomputes dlnV/dlnR-1 independently
  (recompute_sigma_true.f:9,128-135). No dlnV/dlnR consumer.
- A4 CONFIRMED. SIGMA is the analytic derivative of the SAME Steffen
  cubic that interpolates V (:505-517). Measured error on
  V~R^1.5 (exact sigma=0.5): max 5.7E-4 at ND=70
  (vel14_grid_diagnostics.json: power_law_b1.5). Convergence table
  below (B7).
- A5 CONFIRMED. :522 zeroes SIGMA when all nodes are within 1E-4 of
  homology; snhyd_hubble driver run prints SIGMA=0.0d0.
- A6 CONFIRMED. Literals audited: 0.864 expansion factor (:279),
  KAPPA_ES=0.34 (V12 convention, TAU_OPT=0 only), LOG_DEC_CAP=4.0
  natural-log units, ALPHA_R=1.5 blend weight, RHO_SPIKE_FAC=1.5,
  U_GAP_FAC=0.05, NBND_INS=3 (:90-95), 1E10 cm grid unit,
  1E-6*kmax kappa floor (:335), 0.01*TAU(2) seed (:350).

### B. Interpolation numerics (Steffen)

- B1 CONFIRMED. mon_int_funs_v2.f implements Steffen 1990 (A&A 239,
  443): interior first derivative D_i = (S_{i-1}*H_i + S_i*H_{i-1})/
  (H_{i-1}+H_i) (:82, paper Eq. 11a); endpoints extrapolative
  D_1 = S_1 + (S_1-S_2)*H_1/(H_1+H_2) (:80/84, Eq. 11b); limiter
  (sign(S_{i-1})+sign(S_i)) * min(|S_{i-1}|,|S_i|,|D_i|/2)
  (:88-94, Eq. 12) - the sign sum yields 0 at extrema (opposite-sign
  secants) and +/-2*min = min(2|S|,|D|) with correct sign otherwise.
  Hermite coefficients :102-105. Verified line-by-line.
- B2 CONFIRMED. Endpoint slopes use the one-sided extrapolative
  form above, then the same half-magnitude limiter; power-law input
  reproduces sigma=b-1 at nodes 1 and ND to the same 5.7E-4
  tolerance as interior (diagnostics JSON).
- B3 CONFIRMED (zero-gradient census). m15 2.82d ND=70: ZERO nodes
  with |1+SIGMA|<1E-3 (all sigma < -1, min |1+sigma| = 3.7E-2 at the
  outer plateau). Downstream: under DO_NONMONO_SOB the Dessart tau
  floor max(|1+sigma*mu^2|, -0.2*(1+sigma), 1E-4) keeps tau finite;
  the stock path would clamp sigma>=-0.999 via RECOMPUTE_SIGMA_TRUE.
- B4 CONFIRMED (extremum fidelity). Synthetic interior bump
  (2000+6000*exp(-0.5*(lnR/0.25)^2)): grid peak V = file peak to
  5.9E-14; peak R displaced 3.4% (cubic shape, flat-topped so V is
  still exact); 16 sigma<-1 nodes preserved. m15's own peak sits at
  the inner boundary - exact by the RCORE pin.
- B5 CONFIRMED (C1-not-C2). max |dSIGMA/dlnR| jump across inserted
  nodes on m15 = 5.7E4, localized at the CDS cliff cell (the -82
  knot vs file -157) - a real steep feature, not oscillation.
- B6 CONFIRMED with documented divergence. scipy PCHIP vs
  production Steffen on a monotone profile at ND=280: max rel V
  diff 2.2E-7, max sigma diff 1.3E-4 (diagnostics JSON
  pchip_vs_steffen). The limiters differ by design (Steffen bounds
  by adjacent secant min, Fritsch-Carlson by 3x slope ratios); both
  are monotone and exact at knots - intentional divergence, larger
  on sharp features.
- B7 CONFIRMED. Convergence order on the monotone profile
  V=2000*(R/1E4)^0.8*(1+0.1*sin(3*ln(R/1E4))), file NX=ND:

  | NX=ND | max rel V err | max |sigma err| | ord(V) | ord(sigma) |
  |-------|---------------|-----------------|--------|------------|
  | 70    | 4.62E-5       | 4.20E-3         |  -     |  -         |
  | 140   | 5.38E-6       | 1.02E-3         | 3.10   | 2.05       |
  | 280   | 6.46E-7       | 2.49E-4         | 3.06   | 2.03       |
  | 560   | 7.91E-8       | 6.15E-5         | 3.03   | 2.02       |

  Textbook cubic behaviour: V ~ O(h^3), sigma ~ O(h^2). The
  quantitative answer to "is 14 correct": yes, at design accuracy.

### C. U/tau coordinate, steering, budget

- C1 CONFIRMED. U_i = sum_{j<=i} [ max(ltau_j,TAU_LO) -
  max(ltau_{j-1},TAU_LO) + ALPHA_R*|dlogR_j| ], TAU_LO =
  max(ltau_2, ltau_NX - 4.0 ln-units) (:372-379). kappa = file
  Rosseland column * density * 1E10 under TAU_OPT=1 (:337); the
  solver's opacity differs by species/depth but the grid coordinate
  only needs a monotone depth map - the blend is strictly
  increasing because tau is strictly increasing (kappa floored at
  1E-6*kmax > 0, density > 0) and |dlogR| > 0.
- C2 PARTIALLY VERIFIED. kappa is read from the hydro snapshot
  (fixed at file T); grid is built once at setup - per-iteration
  tau drift is inherent to any fixed grid and not specific to V14.
  The fixed-blend coordinate is deliberate (a drifting grid would
  force regridding). No evidence this causes the it~349 stall;
  flagged for the G1 arm.
- C3 CONFIRMED-with-limitation. Trapezoid in R vs a fine-grid
  monotone-cubic reference: max rel err 1.25%, median 0.59% on the
  m15 file (diagnostics JSON tau_quad). Above the 1E-4 checklist
  bound but the tau only sets the blend COORDINATE - 1% tau error
  shifts node placement by <1% of a cell; density-spike served by
  the feature insertion, not the quadrature. Judged acceptable;
  upgrading to the reference quadrature would not move the grid
  measurably (ltau-shift data shows placement is log-driven).
- C4 CONFIRMED on m15. Feature census: 0 reversal entries
  (sigma<-1 globally), 1 sigma minimum (file J=65, sigma=-157 ->
  output knot reproduces -82), 0 density spikes. ND never exceeded:
  insertions only overwrite uniform targets (USED() guard) and the
  U_GAP_FAC floor rejects sub-5%-cell gaps (:427-432). Multi-shock
  stress covered by unit tests (dblshock, rev0/rev1). A monotone
  density STEP (the real CDS, ~100x in 6 cells, <1.5x per step) is
  NOT a detectable feature - recorded as a design limitation in
  the prior audit, still true.
- C5 CONFIRMED. U->R inversion by piecewise-linear LIN_INTERP vs
  Brent on the forward map: max rel err 1.7E-3 - the intrinsic
  interpolation error of the U(R) table, uniform across the grid.
  Placement error is a fraction of a cell; not a bug.
- C6 CONFIRMED. TAU_LO uses ltau_2 (the first REAL tau), not the
  artificial 0.01*seed (:372) - the seed cell cannot dominate U.
  Strict monotonicity of U holds by construction (C1). m15 outer
  values: tau_outer=8.9E-9, V_outer=48.9 km/s; 22 of 70 nodes sit
  below V=100 (the thin atmosphere gets its share via the |dlogR|
  term - the blend answers the NONMONO_KEYS loss-of-depth concern
  quantitatively).
- C7 CONFIRMED. Kappa scaled x1E-3,0.1,10,1E3: R grid identical to
  <=1.3E-14 relative (diagnostics JSON ltau_shift).
- C8 CONFIRMED. Homology rescale R->2R, rho->rho/8, V invariant:
  grid R ratio exactly 2.0, sigma move <=6E-11, inner tau ratio
  0.25000 = lambda^-2. The code has NO hidden tau-invariance bug.
- C9 CONFIRMED. docs/audits/vel14_grid_diagnostics.json committed;
  produced by unit_tests/python/vel14_grid_diagnostics.py.

### D. Density / continuity / expansion

- D1 CONFIRMED. Density is re-interpolated from the (expanded)
  hydro table onto the new grid in log-log by MON_INTERP
  (rd_sn_data.f:281-283) - it is NOT recomputed from Mdot(R^2 rho
  V): the SN snapshot is homologous, not a steady wind, so an
  Mdot-constancy test does not apply (JSON records the
  interpolation-only spread for completeness).
- D2 CONFIRMED. Forward expansion uses +T1 Jacobian
  (1+T1*V/R*(sigma+1))*(1+T1*V/R)^2 at OLD positions identically in
  both readers (snhyd:303-305, rd_sn_data:244-249); VOL_EXP_FAC is
  the Lagrangian inverse -T1 Jacobian evaluated at EXPANDED
  positions (rd_sn_data:256-258) - exact inverse for homology,
  first-order otherwise. Consistent bookkeeping.
- D3 CONFIRMED. No inner extension: RCORE = R_HYDRO(NX) pinned
  (:442,:457); the grid never extrapolates below the file table,
  so no join discontinuity exists.
- D4 CONFIRMED. Clumping: CLUMP_FAC column read when present
  (rd_sn_data:172-174); the file's density is treated as already
  clumped (:279 comment). U coordinate and solver see the same
  file-derived density. m15 file carries no Clumping section ->
  CLUMP_FAC defaults elsewhere; consistent.
- Documented limitation (NEW, not a defect): the expansion factor
  1+T1*V/R*(sigma+1) goes non-positive when T1*V/R*|sigma+1| > 1;
  on m15's sigma=-82 CDS that is dt >~0.06 d past file age. The
  abort at :296-302 is the correct behaviour vs silent corruption,
  but users must keep SN_AGE at/near the file age on strongly
  decelerating profiles.

### E. Downstream coupling

- E1 UNVERIFIED-by-this-audit (out of scope). Resonance-surface
  counting lives in the Sobolev/CMF solver path; the grid layer's
  obligation is preserving sigma sign changes - verified (all 70
  output sigmas < -1 reproduce the file's decelerating field;
  reversal cases in the test battery preserve them). The solver-side
  multi-resonance handling is the audited DO_NONMONO_SOB machinery.
- E2 CONFIRMED already implemented. VEL_LAW=14 does not require the
  nonmono flags; with both off, RECOMPUTE_SIGMA_TRUE:60-85 prints
  the prominent "WARNING: non-monotonic velocity field" banner with
  count/worst cell and guidance - satisfying the checklist's
  abort-or-warning acceptance.
- E3 CONFIRMED-as-guarded. GRID_METH=VEL_SPACE aborts under nonmono
  flags (adjust_r_grid_v4.f:160-169); TAU_SPACE/FULL_R/OB_ONLY
  regrid in R/tau coordinates and re-interpolate V,SIGMA with the
  monotone cubic (:219-229) followed by RECOMPUTE_SIGMA_TRUE
  (solve_for_pops.f:4122). Production decks set REV_RGRID=F, so the
  whole path is unreachable; the abort covers the one unsafe arm.
- E4 CONFIRMED. Dispatch VELTYPE 14 -> TAU_OPT=1 hard-wired
  (cmfgen_sub.f:1098-1100); scripts' cmfgen.vel_law writes VADAT
  VEL_LAW verbatim (no second channel). 13/14 -> SN_MODEL=T,
  SN_HYDRO_MODEL=T (rd_control_variables:250-252); VADAT
  [RSTAR]/[RMAX] rewritten to the hydro edges (:1106-1111).
- E5 CONFIRMED. Homology assumptions: VOL_EXP_FAC Lagrangian factor
  (D2), eval_temp_ddt_v2 uses LOG(VOL_EXP_FAC) for time-dependent T
  (inactive under FIX_T=T on these decks); get_pops_at_prev_time_
  step uses VOL_EXP_FAC for restart rescaling. No R=V*t assumption
  in the population solve.
- E6 PARTIALLY VERIFIED. Impact-parameter grid tangency is a solver-
  side construction on the delivered R grid; strictly-decreasing
  clustered nodes keep tangency (no duplicate R). Ray-count logic is
  outside the audited routine; flagged for the G arms.

### F. Robustness / hygiene

- F1 DONE. The four pinned defects are fixed test-first on
  jfnk-acceptance-audit (unpushed):
  - 4ac2edbc vel14(input): IEEE finiteness screen over all columns
    [test-first] - fixes D1 (NaN-Kappa silently absorbed by the
    MAX floor) and D2 (NaN reached LIN_INTERP's bracket walk).
  - 4ce95eed vel14(input): positivity screen - non-positive
    Radius/Density, negative Kappa abort with column+index+value.
  - 1f025ea1 vel14(input): warn when Sigma column missing or
    all-zero under TAU_OPT=1 (feature machinery silently dead).
  - 0d7d02c0 vel14(test): pinned tests converted to passing
    regression checks + negative cases (D2b/c/d, D3b/c, D4b, D5
    truncated block, zero-Kappa edge). Battery: 69/69 pass.
- F2 CONFIRMED. IMPLICIT NONE present (:23); all locals initialized
  (GOT_SIGMA=.FALSE. :156, INS set per J, USED allocated .FALSE.);
  allocatables all deallocated (:550-551); no SAVE'd state;
  divisions guarded (R>0, V>0, TAU(2)>0 screened before use).
- F3 CONFIRMED. Production objects rebuilt; a -O0 -check all -fpe0
  build (unit_tests/bin/snhyd_dbg) runs clean on the m15 file and
  all synthetic cases including abort paths - no bounds violations,
  no FP traps. (gfortran flags per the checklist map onto the
  ifort equivalents used here.)
- F4 CONFIRMED. Two identical runs produce bitwise-identical stdout
  and NEW_SN_R_GRID.
- F5 guard matrix (all verified live):

  | Guard | Trigger | Verified |
  |-------|---------|----------|
  | ND < 2*NBND_INS+4 | :104 | drv ND=8 abort |
  | SN_HYDRO_DATA missing | :110 | nofile case |
  | NSP/NX/age header missing | :136 | noheader case |
  | R/V/DENS first-node zero | :178 | input screen preempts |
  | non-finite any column | :207-255 | nanr/nanv/nansig/nandens/nankap |
  | non-positive R/Density | :213,:229 | negr/negdens/zerodens |
  | negative Kappa | :240 | negkap |
  | missing/all-zero Sigma | :266 | nosigma/zerosigma warn |
  | input R not decreasing | :281 | badR case |
  | V <= 0 | :289 | badV case |
  | density factor <= 0 | :296 | dblshock@+200d, badexp |
  | expanded R not decreasing | :308 | expandorder case |
  | TAU_OPT=1 kappa all <= 0 | :329 | zerokap case |
  | outer tau <= 0 | :349 | edge case |
  | output R not decreasing | :483 | ordering assert |
  | tight cell report | :495 | fires at 3.0E-5 on m15 |
  | LIN_INTERP strict monotone | lin_interp.f | guarded since 05-Oct-26 |
  | MON_INT_FUNS_V2 monotone | mon_int_funs_v2.f:66 | guarded |

### G. Live-model items - OPEN pending owner

- G1 reproduction of the it~349 r01 stall: needs the bounded Slurm
  arm; instrumentation list agreed (tau drift, sigma extrema,
  min|dV/dR|, JFNK residual, lambda rate).
- G2 A/B VEL_LAW=12/13/14: earlier KNOP cold-start jobs
  13606657-13606661 exist for the V14 arms; the 12/13 arms need
  matching decks.
- G3 RESOLVED statically: the "5-8x resolution" claim
  (NONMONO_KEYS.md ~L150) measures ~4.3x on m15 ND=70 (43 vs 10
  nodes in the CDS+ejecta zone R<36460; 13 vs 3 in the inner cliff
  R<16805; robust to Rosseland vs e.s. opacity choice). Docs updated
  (H3).
- G4 ND sensitivity: run-level; needs owner sign-off for submissions.

### H. Deliverables

- H1 this document. H2 prior audit superseded (header above). H3
  NONMONO_KEYS updated. H4 patch series landed (F1 commits);
  diagnostics JSON committed; unit_tests/README.md row updated.
- H5: no defect makes VEL_LAW=14 untrustworthy as the default - it
  stays the pipeline default.

## 4. Diagnostics table (m15 2.82d, ND=70, TAU_OPT=1, real file)

From docs/audits/vel14_grid_diagnostics.json:

| Quantity | Value |
|----------|-------|
| tau_inner / tau_outer (file, Rosseland) | 37.71 / 1.73E-10 |
| TAU_LO (ln) | -0.370 |
| U range | 7.230 |
| min / max cell dlog_tau (bulk) | 0.0057 / 1.31 (ratio 230) |
| min cell |dlogR| | 3.0E-5 (inner packing) |
| nodes with |1+sigma|<1E-3 | 0 |
| max |sigma| | 82.5 (CDS knot; file -157 diluted) |
| nodes sigma<-1 | 70/70 |
| V outer / inner | 48.9 / 7073.6 km/s |
| nodes V<100 / V>1000 | 22 / 19 |
| U-inversion rel err | 1.7E-3 |
| trapezoid tau err vs fine ref | max 1.25%, med 0.59% |
| ltau-shift grid move (kappa x1E-3..1E3) | <=1.3E-14 |
| homology rescale: grid ratio / sigma move / tau ratio | 2.0 / 6E-11 / 0.250 |
| interior extremum: V err / R displacement | 5.9E-14 / 3.4% |
| V14 vs V13 nodes R<36460 / R<16805 | 43 vs 10 / 13 vs 3 (~4.3x) |
| power-law sigma err (b=1.5, ND=70) | 5.7E-4 |
| convergence order V / sigma | ~3.0 / ~2.0 |

## 5. Open items and deferred questions

1. Monotone density-step detection: the CDS (~100x rise in 6 cells)
   never trips RHO_SPIKE_FAC=1.5 per-cell; resolved anyway via the
   sigma minimum + boundary packing on m15, but a broad shallow dip
   or pure monotone step could under-resolve on other files.
   (Prior audit, still open; candidate F-item for a future round.)
2. Reversal-EXIT detection (sigma crossing -1 upward / V maxima):
   not inserted; m15 has none.
3. Expansion horizon ~0.06 d on sigma=-82 profiles (D-limitation);
   keep SN_AGE at file age.
4. Stale-kappa hypothesis (C2): fixed snapshot opacity vs evolving
   model T - inherent to any fixed grid; G1 arm will measure drift.
5. G1/G2/G4 live arms: bounded submissions pending owner check-in.

Verdict: VERIFIED. The TAU_OPT=1 construction is mathematically and
physically correct at design accuracy; robustness gaps were confined
to input validation and are now closed test-first. VEL_LAW=14 stays
the pipeline default.
