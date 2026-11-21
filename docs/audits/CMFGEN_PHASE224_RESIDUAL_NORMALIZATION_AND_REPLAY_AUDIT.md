# Phase 2.24 - grazing-ray residual normalization and replay audit

ERRATUM (Phase 2.25): the replay driver passed the deck's KN_NC through
while renumbering the ray as LS=1. Every replayed ray was therefore a
core ray, and the P-march tangent row was set to 1E-30 (NaN from the
DIFFUSION I_CORE branch). The per-cell numbers in sections 2-4 and
6.1-6.2 (dark-cell population, max|d| 1.9E-8, verify 3.37E+3) are
superseded by docs/audits/CMFGEN_PHASE225_PHYSICAL_RESIDUAL_AND_ACCEPTANCE_POLICY.md.
With correct boundaries, most bad cells are bright (I>=1E-9) and
fixed-point errors reach 1-24%. "0 floored denominators" still holds.

Date: 20-Nov-2026. Model: m15_2p82d_KNOP (all-compressing; no
mixed-sign inversion). Evidence: replay of the diagdump_i v1
KNRDUMP decks (rays 81/82) with the repaired knray_replay driver;
KN_PSWEEP_DIAG per-pass records (Phase 2.23); KNSEED store dumps.
Artifacts: docs/audits/knop_phase224/.

---------------------------------------------------------------
1. RESIDUAL NORMALIZATION - EXACT SOURCE EXPRESSION
---------------------------------------------------------------

March cells (P channel ~12211-12229; M symmetric ~12387):
  VNEW = FLD(J+1,L)*PK_A0(J,L)
       + SPW(J+1)*PK_A1(J,L) + SPJ*PK_A2(J,L)
       - DSU*PK_A3(J,L) - DSD*PK_A4(J,L)        [P direction]
  VNEW sanitized: <=0 or NaN -> 1E-30 (NNEG); >INEXT_CAP -> cap.
  T1   = ABS(VNEW-VOLD) / MAX(VOLD, 1.0E-4*VMAXF(L))
  RES       = MAX over marched cells this pass
  RES_FREQ(L) = per-frequency max
  VMAXF(L)  = MAX(MAXVAL(|FLD(:,L)|), 1E-25)  live column peak
  significance gate (NTOT/NBAD only):
    MAX(VOLD,|VNEW|) > 1E-4*VMAXF(L)
  RES includes EVERY marched cell - the gate only masks the
  bad-cell census, not the max.

Boundary row (M-channel j=1 block ~12164): T1 = VOLD/MAX(VOLD,
1E-4*VMAXF(L)) - the residual of a zero-flux reset is the old
field magnitude itself; FLD(1,L)=0 follows.

Acceptance: RES < KN_RAY_TOL=2E-3 (tol-exit ixit=1); bulk exit on
the NBAD/NTOT ratio; stall/plateau on rho and counter history.

Confirmed: R(i)=|new-old|/MAX(old, 1E-4*colpeak) as hypothesized,
plus a 1E-25 absolute floor on the column peak.

---------------------------------------------------------------
2. FLOOR-DOMINANCE - DIRECT MEASUREMENT (replay, r0000 deck)
---------------------------------------------------------------

Replay driver now reads v1 decks (section 5). maxp=n runs
reproduce production passes 1..n (fresh adaptive-omega schedule
each call); committed end states for passes 13-16 captured.

Committed-delta census at pass 15->16 (prod normalization):

  cells > tol (sig-masked):     5592  (j1=2030 j2=2764 j3=798)
  denom floored (VOLD<=1e-4vf):    0
  |d| > 1E-9:                      4
  |d| > 1E-12:                  5016

Per-cell components: knop_phase224/ray82P_percell_components.tsv.

THE RESIDUAL IS NOT FLOOR-DOMINATED - it is SELF-normalized:
denominator = the cell's own value ~1E-20..1E-26. A 3.5% wobble
on a dark cell reports res~0.035; the same wobble on a 1E-8 cell
is identical in this metric. The 1E-4*VMAXF floor never binds on
the bad population.

---------------------------------------------------------------
3. NORMALIZATION SENSITIVITY (diagnostic only)
---------------------------------------------------------------

knop_phase224/normalization_sensitivity.md. All relative variants
(prod, |VNEW|, column-peak) give ~5.3-5.6K cells over 2E-3 - the
bad population is insensitive to which relative denominator is
used. The absolute metric: max 1.87E-8, median 2E-15 intensity
units - physically negligible on every cell.

Interpretation of denominators: prod's denominator tracks the
cell's own magnitude, so it measures RELATIVE iteration error
per cell - appropriate for cells carrying flux, meaningless for
dark cells. Column-peak normalization (|d|/vmaxf) is the "local
radiation-field scale" variant - nearly identical result here
because dark columns have dark peaks. Only the absolute metric
exposes the real accuracy: the map is converged to <2E-8 in
intensity units everywhere.

---------------------------------------------------------------
4. DS-CORNER JACOBI LAG - AUDIT + SIGNED DATA
---------------------------------------------------------------

Source trace: knop_phase224/ds_corner_dependency.md. Confirmed in
code: the P march previews row J-1's source from its PRE-march
state (KN_ORD_SPW_FILL on FLD before the row is marched); DS(J)
for the interior branch reads SEC(J-1)=(SPW(J)-SPW(J-1))/DTAU,
so DSD at row J carries row J-1's previous-iterate source - a
deliberate one-pass Jacobi lag per Knop Paper II Eqs 31/32, not
a bug. It vanishes at a fixed point; it can only slow or perturb
convergence, not create a bias.

Signed evidence (production trajectory):
- bad population confined to j=1-2 (+row-3 tail), zero at j=4:
  the boundary pair where the lag lives.
- j1/j2 deltas move in lockstep same-sign (j2 attenuated) - a
  coupled boundary-pair drift, not an alternating mode.
- d1 sign pattern on the worst cells: +,+,- over passes 13-16;
  |D2|~|D1| - a slow drift with a collective correction event,
  not a period-2 orbit. 1130/5593 and 2143/5593 of the bad cells
  flipped sign at the two transitions - a minority, mixed-sign
  population.

---------------------------------------------------------------
5. REPLAY COMPATIBILITY - FIXED
---------------------------------------------------------------

drv_knop_replay.f90 changes (diagnostic driver only; solver
untouched, git-visible diff):

- Explicit version table: accepts 20261122 (v1) and 20261123
  (v2); any other magic stops with
  RESULT deck_version_unsupported=<magic> reader_supports=...
- v1 decks: KN_CHI_EXT/KN_ETA_EXT reconstructed as
  EXP(KN_CHI_COEF(K,4,L)) - the exact node values of the
  MON_INT_FUNS_V2 monotone cubic built on LOG(EXT) (records
  deck_ext_mode=coef_node). JMP_SIB zero-filled (diagnostic
  prints only; deck_sib_mode=zero). v1 FBITS masked to bits 0-8.
- Documented limitation: v1 never wrote seed/seeded/cert flags;
  replay treats the ray as unseeded.

Manual relink (lib/libsubs.a currently holds 117 members and no
unix objects after the 18:15 rebuild; link uses the prod_obj
cache + the production link line: ifort -qopenmp -mkl -nbs +
MKL lists). knray_replay rebuilt at ~23:00.

Replay results (deck KNRDUMP_L082_r0000.bin, refresh r0000):
  base    res=6.63E-2 ixit=3 npo=25 omg=0.450
          (prod:  res=6.47E-2 ixit=3 npo=15) - stall reproduced
  verify  res=3.37E+3 - one-pass unrelaxed defect at the
          stalled state: a single pathological dark cell
          (max over 225K cells), consistent with the
          denominator issue.
  hst 30  (new mode; KN_REPLAY_HST.bin per-pass committed
          fields): the omega=1 pass-at-a-time trajectory
          CONTRACTS - ~10-62 cells>tol by pass 15-30 vs ~5593 on
          the production trajectory. The stall depends on the
          in-solve controller history (omega ratchet 0.45),
          not only the map.

Signed-history mode hst added to the driver (writes x(n) per
one-pass solve + x(0)=seed); plus KN_REPLAY_END.bin dump in
base/maxp modes for production-trajectory states.

---------------------------------------------------------------
6. FINAL ANSWERS
---------------------------------------------------------------

1. Floor-caused residual? NO - zero floored denominators among
   the 5592 bad cells. Denominator = the cell's own dark value.
   Correcting the earlier "floored-denominator" phrasing: it is
   a self-relative residual on dark cells, not a floor artifact.

2. Absolute defects negligible? YES - max 1.87E-8, median 2E-15
   intensity units over the committed step; argmax cells ~1E-25.

3. DS-corner lag explains the history? PARTIALLY - mechanism
   confirmed in source (intentional C-corner lag); the bad
   population occupies exactly the boundary pair j=1-2 and
   moves in lockstep. Causal share not isolated - needs
   per-cell source-term instrumentation.

4. Oscillation or bias? SLOW COUPLED DRIFT with occasional
   collective corrections - minority sign flips, |D2|~|D1|,
   lockstep j1/j2. NOT a period-2 mode. NOT a fixed-point bias
   (the lag vanishes at convergence).

5. Does the acceptance test distinguish dark-but-converged from
   genuinely-wrong? NO - a self-relative per-cell tolerance
   cannot. The transfer map is converged in intensity units on
   ray 82P; the reported residual is a metric artifact on cells
   that carry no flux. This is a tolerance-design issue, not a
   transfer-accuracy failure, on this ray.

6. Safe explicit version handling? YES - implemented (version
   table, clear diagnostic, v1 reconstruction without
   fabrication). See knop_phase224/deck_format_spec.md.

7. Evidence P2DAMP helps? NONE - no demonstrated per-cell
   period-2 mode; the detector would not trigger on the
   measured sign structure. Also, loop-mode convergence at
   omega=1 says the map contracts - damping a contracting
   drift gains nothing and could mask the metric issue.

8. Evidence AA helps independent residual? INDIRECT - AA
   accelerates contracting marginal modes (rho~0.99 measured);
   the loop data show the map does contract, so AA could
   plausibly shorten the stall - but the "residual" is metric
   noise on dark cells, so it fixes a reporting problem.
   Untested; not enabled.

9. Outstanding before damping tests: decide whether the
   acceptance metric should be changed (absolute-vs-relative
   residual policy - a solver change outside this phase's
   scope); v2 deck capture for ray 81 + a converged control;
   per-cell source instrumentation if the lag's causal share
   must be quantified.

10. Mixed-sign status: unchanged - unvalidated; no genuine
    inversion model exists.

---------------------------------------------------------------
7. PROVENANCE AND LIMITATIONS
---------------------------------------------------------------

- Decks: KNRDUMP_L081/L082_r0000.bin (diagdump_i job 13433532;
  v1 writer in live exe d50dcacd). Replayed refresh r0000 only;
  decks carry ~34 appended refreshes not yet replayed (reader
  uses the first header - later refreshes reachable by seeking
  past per-solve trailers; not implemented).
- Replay binary: knray_replay @ 22:59; drv_knop_replay.f90 with
  the v1/hst/end-dump edits (uncommitted, diagnostic-only);
  production objects from unit_tests/bin/prod_obj (fg_j .o of
  17:30 vintage = v1-writer source, matching deck provenance);
  driver compiled -O0 vs production -O3 - FP-level drift
  expected; replay_vs_prod_maxdiff=1.7E-2 recorded.
- lib/ is mid-rebuild (libsubs.a 117 members, no libpgplt.a);
  the manual link sidesteps it. run_knop_replay.py stays broken
  until the archives are restored - documented, not worked
  around in the build system.
- Converged-control decks do not exist (only rays 81/82 armed);
  control comparison used the PSD per-pass records instead.
- No production solver, tolerance, cap, march, or config change.
  Jobs: A1/A2 RUNNING; B1/B2/diagcap/obs PENDING (QOSGrpCpuLimit).
- M-channel history not captured (hst records FLDP only; ray 82's
  M channel converges anyway).
