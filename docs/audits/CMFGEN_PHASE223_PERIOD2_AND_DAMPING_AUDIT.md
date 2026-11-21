# Phase 2.23 - KNOP outer-boundary oscillation and damping audit

Date: 20-Nov-2026. Model: m15_2p82d_KNOP (all-compressing ejecta;
v(R) monotone 48.9 -> 7073.6 km/s; NO mixed-sign inversion).
Evidence: A1/A2 OUTGENs (jobs 13414391/13414398, frozen 93364cac,
KN_MAT_FREEZE=T), diagcap_i 13430220 (interactive, OUT_OF_MEMORY at
66 min), diagdump_i 13433532 (interactive, TIMEOUT at 2h,
KN_PSWEEP_DIAG.dat 12 MB + KNRDUMP decks rays 81/82).
Artifacts: docs/audits/knop_phase223/.

---------------------------------------------------------------
0. CORRECTIONS TO THE PHASE 2.22 OFFENDER PICTURE
---------------------------------------------------------------

Two earlier claims need correcting on the new per-pass (PSD)
evidence from jobs 13430220/13433532:

1. "Rays 22/45/50 are persistent non-committing production
   offenders" - WRONG. PSD final records show rays 22/45/50 exit
   IXIT=1 (tol-converged) on both channels every refresh, res
   4E-4..1E-5, in 1-5 passes. Their KN-CEN flagged cells are the
   largest RELATIVE committed movers among converged rays
   (floor-level columns ~1E-28 peak; VR~0.499 vs jump-veto 50) -
   store-scan noise, disjoint from the true offender set.

2. "The three production rays with q>1E-2 are 22/45/50" -
   unsupported inference. The q census counts, it does not name.
   The stall/capped pair is rays 81P and 82P; the third q-flagged
   production ray is most plausibly ray 12 (M-channel res range
   4E-4..1.9E-2 over solves) or ray 80P on marginal refreshes.

Corrected population map (mature state, ~43 solves captured):

  - committed clean (81 rays): all others incl. 22/45/50
  - capped commits, verdict 3 (2 rays): 81P, 82P (IXIT=3 stall)
  - unsolved (1): ray 85 (NI=1 degenerate grazing)
  - probe-suppressed (1): ray 83 (IXIT=3 stall, verdict 7)
  - KN-CEN flagged (store-jump noise): {22,45,50 ch1},
    {83 ch2 probe}

---------------------------------------------------------------
1. PER-PASS RESIDUAL ATTRIBUTION (now resolved via PSD data)
---------------------------------------------------------------

The per-pass transfer residual on the stall rays sits at jmx=2 -
the boundary-adjacent second ray-point on the outward (P) channel -
every pass, while the frequency argmax wanders (lmx 1302..18733
across EUV/soft-X-ray slots). Sign/GAM at the argmax: +1
(expanding branch, ascending frequency order); NFM_L=1.

Final-pass records (last solve of the capture):

  ray/ch res     dmax    dabs     jmx lmx   nbad   ntot   nflr  g nfm
  81P   3.16E-2  1.4E-2  5.1E-9    2  14406    830  183711  9344 + 1
  82P   7.26E-2  3.3E-2  6.2E-10   2  10278   6007  216939  7521 + 1
  83P*  6.89E-1  3.1E-1  3.8E-8    2  14403  80029  162930  5700 + 1
  *probe

  Controls (converged grazing rays, same cell class):
  76P   3.39E-3  1.7E-3  3.9E-9    1  14406     16  545307 18463 + 1
  77P   4.35E-3  2.1E-3  7.1E-26   1   1437     23  490728 16665 + 1
  78P   5.46E-3  2.7E-3  7.3E-26   1   1437     41  436175 14841 + 1
  79P   5.95E-3  2.9E-3  7.3E-26   1   1437     63  381633 13006 + 1
  80P   7.16E-3  3.5E-3  4.0E-27   2   1302    116  327098 11164 + 1

  KN-CEN rays (converge cleanly; argmax at DEEP real cells):
  22M   4.47E-4  jmx=94 lmx=17881  dabs=2.2E-5   nbad=0
  22P   3.76E-4  jmx=99 lmx=17881  dabs=2.0E-5   nbad=0
  45M   9.04E-6  jmx=39 lmx=21328  dabs=5.6E-10  nbad=0
  50M   1.02E-5  jmx=34 lmx=21328  dabs=1.6E-10  nbad=0

Facts:

- The stall set {81P,82P,83P-probe} is exactly the three shortest
  non-degenerate grazing rays (NI=5,4,3; tangent depths 5,4,3).
  M channels of the same rays converge in 2 passes (1E-5..1E-9).
- The argmax is a POPULATION, not a cell: nbad ~830-80K cells above
  tol on the stall rays, all j=2-adjacent floor columns; lmx wanders
  between passes. Absolute defects |dV| <= ~4E-9 (argmax), and the
  q absmax bound over all cells is 1.5E-4 - the map is converged in
  intensity units; the residual is floored-denominator relative
  noise on physically dark cells.
- Converged grazing controls 76-80 carry the same j=1-2 floor class
  at 2.6-7.2E-3 - they squeeze under the 2E-3 bulk criterion; the
  stall rays are the shortest, where the dark-column fraction and
  the boundary-lag coupling are largest.
- Mechanism (from section 4 of the user's stall audit, confirmed):
  the ordered P march sweeps J=NI-1..1; DS(J,L) reads SPW(J-1,L)
  built from row J-1's PRE-march source - a genuine one-pass Jacobi
  lag that lands exactly at j=2 because row 1 is the outermost
  (vacuum-boundary) row. On dark columns this lag cannot contract
  below the floored-denominator tolerance.

---------------------------------------------------------------
2. PERIOD-2 VERDICTS (state-level tests)
---------------------------------------------------------------

Method: D1=x(n)-x(n-1) sign flips and |D2|/|D1| on the available
signed state series. A clean multiplicative period-2 mode shows
~100% flips and |D2| << |D1|.

a) Offender-ray store cells (rays 22/45/50 ch1 flagged cells,
   KN-CEN v2 store values, 263 refreshes): 90% D1 sign flips,
   |D2|/|D1|=0.29, fit x=a+bn+c(-1)^n gives alternating amplitude
   0.031% of mean vs 0.13% total drift. A REAL square-wave
   alternation of the retained store - but at ~1E-40 absolute
   amplitude on floor-level cells of CONVERGED rays. Confirms the
   oscillation machinery exists; not the residual driver.

b) Probe ray 83 retained-store seed (KNSEED_P dumps, KN_DK_XREF,
   264 refreshes): at j=2 (the argmax depth class, real ~1E-8
   cells) monotone drift, |D2|/|D1|~2.0, 0-1 flips - NOT period-2.
   At j=1 floor cells ~45% flips, |D2|/|D1|~1.0-1.3 - jitter
   around pinned values, not a limit cycle.

c) Ray 82P solve-residual sequence (43 solves): 63% D1 flips,
   |D2|/|D1|=0.91 - oscillatory, NOT a clean period-2 orbit.
   Ray 81P: 63% flips, |D2|/|D1|=1.88 - trend+jitter.
   Probe 83P: 49% flips. Converged control 76P: 61% flips around
   its accepted 3.4-4.4E-3 plateau - same jitter character.

d) Ray 82 signed per-cell/store history: UNAVAILABLE in any
   produced telemetry. Requires a v2-format deck replay (section
   5) plus a diagnostic print of KN_SGN_SUM or per-pass iterate
   sign - proposed, not implemented.

Verdict: the refresh-level store oscillation is CONFIRMED
period-2-like at micro-amplitude on floor cells. The ray-82
solve residual is oscillatory (the omega controller ratchets to
the 0.45 floor; NOSC/growth damping fires) but a strict period-2
orbit of its argmax population is NOT demonstrated - consistent
with a large pool of near-equal floor cells rather than a single
locked 2-cycle.

---------------------------------------------------------------
3. RESIDUAL FAMILY SEPARATION (source definitions)
---------------------------------------------------------------

Full detail: knop_phase223/residual_and_metric_definitions.md.

- A. Per-pass transfer residual RES (MARCH ~12214): max normalized
  overwrite defect |VNEW-VOLD|/MAX(VOLD,1E-4*VMAXF) over the pass;
  VMAXF floored at 1E-25. Measured pre-relaxation, uncapped by the
  significance gate.
- B. KN_Q_RESID (~10596): one independent scratch verify pass per
  channel, ||T(x)-x|| same normalization; no store writes.
- C. FG_GS_MAXCHG/GS_RF_DCHG (~6040, ~9853): lagged-store vs fresh
  sweep comparison and retained committed-store delta.
- D. FG_COMMIT_NC cap (~14420): |V1-V2|<=0.05*MAX(V2,1E-4*VMAX_F)
  provisional bound; saturating commit returns exactly 0.05.

They measure different stages of the refresh map. A small committed
store change CAN coexist with a large transfer residual - observed:
rays 81/82 stall at RES 3-7E-2 while their capped stores ratchet
bounded by 0.05, and committed rays sit at <=4.9E-3 retained delta.
GS_RF_DCHG pinned at 0.05 = the saturated cap, not a convergence
certificate; FG_RF_ALL_CONV remains gated on the full bundle.

---------------------------------------------------------------
4. THE THREE DISJOINT POPULATIONS (reconciliation)
---------------------------------------------------------------

  store-drift argmax (FG_GS_MAXCHG): ray 82, depths 1-2
     -> SAME ray as the transfer-residual offender. One cause:
     its stall exits drive capped commits, which produce the
     largest retained delta by construction.
  per-pass residual argmax: j=2 cells on 81P/82P (+probe 83P)
  KN-CEN flagged set: 22/45/50 ch1 depths 2-5 - converged rays,
     floor-noise store movers, unrelated to the residual.

The "common outer-boundary problem" hypothesis survives in a
narrower form: all flagged/argmax cells sit in the outer-edge
shell (dark EUV columns, thin BC, weak |sigma|), and the stall
rays are the shortest grazing rays where the boundary-lag term
and the floor-cell fraction dominate. But rays 22/45/50 do NOT
share the residual failure - the earlier lumping was wrong.

---------------------------------------------------------------
5. DIAGNOSTIC INFRASTRUCTURE STATUS
---------------------------------------------------------------

- KN_PSWEEP (all-ray per-pass capture): present in ALL pinned
  binaries AND in the live d50dcacd build. diagcap_i/diagdump_i
  produced the data used above. diagcap 13429272 (small partition)
  remains PENDING on QOSGrpCpuLimit - preserved, its data would add
  refreshes beyond the ~43 already captured.
- KN_RDUMP replay decks: v1 decks for rays 81/82 exist (14 GB each,
  ~34 refreshes appended). Not replayable: v1 omits
  KN_CHI_EXT/KN_ETA_EXT that BUILD_FREQ's on-node branch reads;
  replay segfaults inside KNOP_RAY_SOLVE_V1 (job 13444839). v2
  writer (magic 20261123) exists only in the uncommitted Phase 2.22
  source - needs a rebuilt exe + a fresh short armed capture.
  Details: knop_phase223/deck_provenance.md.
- knray_replay: built (Nov 20 18:09) but the link env is currently
  incomplete (lib/libsubs.a 117 members, missing unix objects);
  prod_obj cache has all needed members for a manual relink.
- Proposed diagnostic-only addition (NOT implemented): print
  KN_SGN_SUM per pass in the replay driver (or arm KNDUMP in
  replay) to get the signed per-pass/per-cell residual that the
  period-2 proof needs. Smallest change: one WRITE in
  drv_knop_replay.f90 loop mode.

---------------------------------------------------------------
6. DAMPING OPTION AUDITS (source-level; nothing enabled)
---------------------------------------------------------------

Full detail: knop_phase223/damping_option_audit.md.

KN_ORD_P2DAMP (~12052, ~12254-12295, ~12350-12384):
per-cell per-PASS sign-flip detector on significant cells;
FG_ORD_P2N=2 consecutive flips -> damped commit = bracket midpoint
(sqrt(VOLD*VNEW) positive endpoints, else arithmetic). Releases
after 3 flip-free passes. Iteration-path change only; RES still
measures the raw map defect -> no false convergence. Applies only
to cells above the significance gate. VADAT key is wired in the
pinned 93364cac binary (KN_ORD_P2DAMP, KN_ORD_P2N, KN_ORD_P2REL).

FG_ORD_AA (~12753, ~12907-13030, rollback ~12603): safeguarded
type-II Anderson on the pass iterate, M=2 (3 pairs), weighted Gram
+ ridge, veto ||corr||>1E3||f||, 10|r| cap, positivity fallback,
residual-acceptance rollback disabling AA for the solve. Keys
KN_ORD_AA/KN_ORD_AA_M present in the pinned binary.

Neither touches the floored-denominator metric that generates the
reported residual. If the diagnosis "metric artifact on dark
cells" holds, the decisive test is an absolute-denominator
residual variant (a solver change - out of scope here).

---------------------------------------------------------------
7. EXPERIMENT-ARM DESIGN (knop_phase223/experiment_arm_design.md)
---------------------------------------------------------------

B: baseline = A1/A2 + existing PSD captures (no new job).
C: KN_ORD_P2DAMP=T only, on pinned 93364cac (pure config change).
D: FG_ORD_AA=T only, same. E: combined only if C or D helps.
All arms: same checkpoint/store as A2; record ixit census,
verdicts, NP2/AA counters, GS_RF_DCHG, q census, refresh count to
full certification. NOT submitted: QOSGrpCpuLimit holds every
pending job; interactive diagnostic captures are a separate QOS.

Free-matter arm (recap): clone B2 deck, pinned 58a699d0,
KN_MAT_FREEZE=F, identical donor hashes. Not submitted.

---------------------------------------------------------------
8. FINAL REPORT - ANSWERS
---------------------------------------------------------------

1. Is ray 82 demonstrably period-2? Oscillatory: YES. Period-2 at
   the solve level: NOT demonstrated - 63% flips, |D2|/|D1|=0.91
   across 43 solves, argmax wandering over a large floor-cell pool.
   Clean period-2 alternation IS demonstrated in the retained
   store at floor cells (90% flips, micro-amplitude) on converged
   rays. Signed per-cell proof for ray 82: pending v2 deck +
   KN_SGN_SUM print (proposed).

2. Does KN_ORD_P2DAMP act on the state responsible? Partially.
   It damps per-pass sign-flipping cells above the significance
   gate - the stall argmax cells qualify when |VNEW|>1E-4*VMAXF.
   Whether those cells flip per-pass is unproven; the documented
   trigger case ("~20 near-floor cells 2x/0.5x flip-flop") is the
   same cell class, so it is a plausible instrument, not a proven
   one. It cannot act on cells below the gate.

3. Are rays 22/45/50 failing for the same per-pass reason? NO -
   they do not fail at all. They converge IXIT=1 every refresh;
   their flagged cells are store-scan noise on converged floor
   columns. The earlier offender-set assignment is corrected
   (section 0).

4. Mechanism of the real offenders (81P/82P): the downwind
   DS-corner one-pass lag at the outer boundary (row 1's pre-march
   source feeds row 2's slope terms) on dark EUV columns, reported
   by a floored-denominator relative metric. In absolute units the
   map is converged (<=4E-9 argmax, <=1.5E-4 global).

5. Does the per-pass solve converge while the refresh map fails?
   Inverted: the per-pass solve on 81/82 does NOT reach KN_RAY_TOL
   (it stalls on floor-normalized defects); the committed field is
   converged in absolute units. The refresh store drifts because
   capped commits bound the motion, not because the map diverges.

6. Would damping reduce KN_Q_RESID? Unknown; the q tail lives on
   the same floor-cell class. If P2DAMP's bracket collapse reaches
   those cells it could help; AA targets the rho~0.99 marginal
   mode. Neither is established to hit the argmax pool - the sign
   test must run first.

7. Would it improve the certified-ray count? The 2 capped rays are
   the entire non-certified production set; a genuine fix would
   move 81/82 to verdict 0 and let FG_RF_ALL_CONV certify (with the
   full bundle). A damping that only quiets the metric without
   converging the map would be caught by q_resid + capped commits.

8. Reproducibility: the diagcap_i arm replayed A2's refresh-1
   census bit-identically, so an experiment arm on the same
   checkpoint + pinned exe reproduces the baseline exactly.

9. Free-matter: designed, not submitted; frozen-matter arms can
   not test release/freeze behaviour.

10. Missing before mixed-sign trust: a genuine inversion model
    with active seeded cells on both sides of GAM=0, plus the
    per-cell sign-level oscillation evidence above.

---------------------------------------------------------------
9. OUTSTANDING QUESTIONS / NEXT INSTRUMENT
---------------------------------------------------------------

- Signed per-pass delta history on rays 81/82 (v2 deck + driver
  print) - the decisive period-2 datum.
- Absolute-denominator residual variant test (user Task-4).
- Whether the j=2 argmax cells flip sign per-pass (P2DAMP trigger
  condition) - measurable from the same replay.
- B1/B2, diagcap 13429272, observer jobs: still pending; preserve.
- No solver, tolerance, cap, donor, or boundary changes made. No
  commits. Working tree (incl. the user's uncommitted Phase 2.22
  rdump machinery) untouched.
