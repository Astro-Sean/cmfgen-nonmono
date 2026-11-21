# Phase 2.25 - physically weighted KNOP residual and acceptance-policy audit

Date: 20-Nov-2026. Model m15_2p82d_KNOP (all-compressing; no
mixed-sign inversion). Offline analysis of the diagdump_i v1 decks
(rays 81, 82; refresh r0000) with the replay driver, plus
KN_PSWEEP_DIAG and OUTGEN records. Artifacts: knop_phase225/.
No production solver, threshold, damping or configuration change.

Status: PARTIAL. Two replay-fidelity problems were found. One is
fixed (section 0); the other is a known limit (section 1). Several
required items stay open (section 10).

---------------------------------------------------------------
0. REPLAY DRIVER BUG - PHASE 2.24 NUMBERS ARE SUPERSEDED
---------------------------------------------------------------

drv_knop_replay.f90 renumbers the dumped ray as LS=1 but passed the
deck's KN_NC=15 through. Every replayed ray was therefore a core
ray. For ray 82 (P=1.025E5 > R_ND) the P-march inner row took the
DIFFUSION I_CORE branch (fg_j_cmf_nonmono_v1.f ~12307), where
SQRT(R_ND**2-P**2) is NaN. That was sanitized to 1E-30 on the whole
tangent row instead of IM_NI(L).

Fix (driver only): KN_NC=MERGE(1,0,LSX.LE.KNC). NC_PNT_SRCE_L is a
solver PARAMETER and cannot be remapped, so a PNT_SRCE deck now
prints a warning; these decks use DIFFUSION. After the fix the
replay-vs-production max abs difference falls from 1.7E-2 to
2.4E-4 (ray 82) and from 1.9E-2 to 2.7E-4 (ray 81).

Consequences for Phase 2.24 (erratum added there):
- The "dark-cell, self-relative" picture came from the corrupted
  tangent row. With correct boundaries, 4394 of the 5146 bad cells
  on ray 82P have I >= 1E-9.
- The "max |d| = 1.9E-8, median 2E-15" values are superseded
  (section 2).
- The "0 floored denominators" result still holds after the fix.
- The verify-mode value 3.37E+3 was an artifact; it is now 0.999.

---------------------------------------------------------------
1. REPLAY FIDELITY (corrected binary)
---------------------------------------------------------------

base mode, against the production trailer and FLDP_END in the deck:
  ray 82P: replay res=5.99E-2 ixit=3 npo=27 | prod 6.47E-2 ixit=3 npo=15
  ray 81P: replay res=1.86E-2 ixit=3 npo=23 | prod 2.22E-2 ixit=3 npo=27
  field, cells with I>1E-9: replay end vs prod end, median relative
  difference 3.7E-2 (82) and 3.3E-2 (81), p99 about 5E-2, roughly
  uniform across rows including the M-fed tangent row.
  Production itself moved only 4.3E-4 (median) from seed to end.

The stall class reproduces. The field does not reproduce better
than about 3-4%, and the offset is already present on the
M-channel row. That points to an input difference, not an
iteration difference. Leading hypothesis (unverified): the v1 deck
has no KN_CHI_EXT/KN_ETA_EXT, so on-node values are rebuilt from
the coefficient cubic. That is exact only if the cubic and the
snapshot come from the same iteration; if production used a
different EXT snapshot, the replay solves a slightly different map.
Settling this needs a v2 deck.

Rule for this report: replay numbers are evidence about the map's
behaviour (convergence, sign structure, row localisation), not
production field values at better than about 4%.

---------------------------------------------------------------
2. MAXIMUM ABSOLUTE DEFECT (corrected replay; committed step of
   production-trajectory pass 15->16, maxp runs)
---------------------------------------------------------------

Ray 82P maximum |d| cell:
  j=2, r=1.2530E5, l=18731, lam=1033.30 A (O VI 1032/1038 region)
  vold=4.8088E-6 vnew=4.7720E-6 d=-3.684E-8 |d|=3.68E-8
  column peak=4.817E-6; S=eta/chi=9.06E-5; chi=2.24E-9;
  GAM=+8.55E-10; production denominator=4.809E-6 (not floored);
  res=7.66E-3 > KN_RAY_TOL; in the bad set: YES.
  Local relative defect is 0.77% per committed pass on the cell
  that is the brightest in its column.
Ray 81P maximum: j=1, l=14402, 568.8 A, vold=1.175E-5,
  d=-4.26E-8, res=3.62E-3, in the bad set: YES.
SIGMA is not stored in the v1 deck; see
knop_phase222/ray_geometry.tsv (sigma about -1.04..-1.18 at
depths 1-5).

The 100 largest |d| (knop_phase225/ray8{1,2}P_top100_absdefect.tsv)
span 324-1641 A; 85/100 (82P) and 28/100 (81P) are in the bad set.
So the maximum-absolute-defect cells and the bad population
overlap. The bad set is not a separate population of dark cells.

The bad set by intensity (82P, 5146 cells):
  I>=1E-9: 4394 | 1E-12..1E-9: 524 | 1E-16..1E-12: 148 |
  below 1E-16: 80. Floored denominators: 0.
Classes (reference Imax=1.70E-2): A |d|>1E-6*Imax: 4;
  C mixed: 4578; B |d|<=1E-10*Imax: 564; D non-finite: 0.
Ray 81P: 97 bad cells, 51 of them below 1E-20 (dark), 33 at I>=1E-9.

Per-pass defect vs error. d is omega*(T(x)-x) with omega=0.45.
The distance to the fixed point is about |T(x)-x|/(1-rho). The
measured rho is 0.93-0.99 (late hst-mode |d| ratio median 0.927;
production rho_max 0.9895), which gives a factor of about 30-220
on d. Measured directly: the omega=1 pass-at-a-time replay reaches
max|d| ~2E-11/pass by pass 30, and stagnates there with about 25
cells above tolerance. The production-trajectory state (maxp 16)
then differs from that reference by median 2.9E-4, p99 3.3E-2 and
max 0.24 relative (I>1E-9; rows 1-3; row 4 agrees to 4E-8).
On those cells the pass-to-pass defect understates the error by
one to two orders of magnitude.

Answer: the defect is NOT negligible against the local intensity
scale on the bright EUV/FUV cells of the boundary pair.

---------------------------------------------------------------
3. PHYSICALLY WEIGHTED MEASURES (committed step, d-based)
---------------------------------------------------------------

Weights: actual nonuniform KN_FREQ grid, |gradient(nu)| bin widths.
No floor is used: S_INT>0 on every row. Artifacts:
ray8xP_band_weighted.tsv.

R_INT = sum|d|dnu / sum|I|dnu
  82P: j1 2.6E-6, j2 3.2E-6, j3 4.9E-7, j4 0 | whole ray 1.16E-6
  81P: j1 1.9E-6, j2 1.7E-6, j3 8.8E-8, j4 1.5E-8 | ray 5.1E-7
Bad set: 63% (82P) and 28% (81P) of E_INT, but only 1.4E-4 and
4.2E-5 of S_INT.
Bands: 228-912 A carries 70% (82P) and 79% (81P) of E_INT, with
R=8.5E-6..1.0E-5 (82P). Below 228 A, R<=2.4E-6.

Rate proxies (no cross sections; per-row relative weighted defect):
  photon number 1/nu: 2.0E-6 / 2.5E-6 (82P j1/j2)
  H I-like sigma~(nu/nu0)^-3 above 912 A edge: 2.0E-5 / 2.6E-5
  He II-like above 228 A edge: 9.8E-6 / 1.1E-5
The H I-weighted defect is about 10x the intensity-weighted one.
High-weight bands amplify the defect, as the brief anticipated.
These are proxies only. Real rates use the model's cross sections
and line profiles, which this phase did not evaluate.

All values are per committed pass. With the 30-220x factor above,
the corresponding fixed-point errors are of order 1E-4..5E-3
(H I-weighted).

Ray share of J (estimate). Trapezoid-in-mu weights from
p=R(86-LS) (ray_jqw_estimate.tsv): ray 82 carries about 0.15 of
the [0,1] mu range at depths 1-2. Ray 81 carries about 0.09.
Under an isotropic proxy, a cell error e on ray 82 maps to about
0.15*e in J at that depth and frequency. The production JQW
values were not read; this is an estimate only.

Share of the total intensity-weighted residual from 81P/82P: NOT
MEASURED. Only rays 81 and 82 have decks. The PSD records give an
argmax |d| per ray, not integrals.

---------------------------------------------------------------
4. POPULATION AND SPECTRUM PATHWAYS (source trace)
---------------------------------------------------------------

Rays -> J: JNU_STORE(I) += JQW(I,LS)*AV(K,LS), all rays and depths
(fg_j_cmf_nonmono_v1.f ~6310). The 81/82 grazing rays enter J at
depths 1-5 only.
J -> rates/populations: J moments feed the radiative rates in the
CMFGEN_SUB/SOLVE_FOR_POPS chain. On frozen-matter runs
(KN_MAT_FREEZE=T) populations are held, so no feedback was
observed.
Spectrum: CMF_FLUX_SUB_V5 recomputes its own formal solution
(COMP_JCONT_V4 KNOP dispatch) from the populations.
The population-solver ray field reaches the spectrum only through
populations, so it cannot change it directly. No observer A/B was
run: the obs jobs are pending, and identical-population inputs
would show no difference by construction.
Population sensitivity: NOT MEASURED (needs a free-matter run;
design in knop_phase222/freematter_design.md).

---------------------------------------------------------------
5. DS-CORNER LAG (quantitative status)
---------------------------------------------------------------

Lag cells: every P-march row J uses SPW(J-1) built from row J-1's
pre-march state inside DSD=DS(J). Row 1, the vacuum-boundary row,
is the last row updated.
Measured: committed d on rows 1 and 2 has the same sign in 95.6%
(82P) and 98.3% (81P) of columns. sum|d| row1/row2 = 0.78 (82P)
and 1.11 (81P). Rows 3-4 carry <=10% of the bad cells (row 4:
0 on 82P).
The fraction of the residual caused by the lag is NOT
QUANTIFIED. Separating it needs the pre- and post-march S'(J-1)
and the PK_A3/PK_A4 coefficients per cell, which the driver cannot
see. Proposed diagnostic (not implemented): an existing
shadow-evaluation path (KN_ORD_SPV_SHADOW, the frozen DS
replicate near ~14148) could export DSD(pre)-DSD(post) per cell
under a replay-only flag. That is a solver-file diagnostic, which
needs review before any implementation.
No control-ray decks exist, so the control comparison could not
be made.

---------------------------------------------------------------
6. ACCEPTANCE POLICIES (knop_phase225/policy_comparison.tsv)
---------------------------------------------------------------

Census context (OUTGEN, latest refresh): full=59 bulk=16 stag=8
inv=1 of 84 active rays; rf_all_conv=F. The non-FULL set is 25
rays. The 16 BULK rays (24M-26M, 68P-80P) also block
FG_RF_ALL_CONV. Rays 81/82 are not the only blockers.

A strict (reference): 59 accepted. No false-convergence pathway.
B dual local/global: with a 6E-2 local cap, 75 rays pass. The
  global R_INT ~1E-6 would pass rays whose bright O VI cells carry
  1-24% fixed-point error. That is a false-convergence pathway
  unless the local cap is near KN_RAY_TOL, and then nothing
  changes.
C significance mask by integral share: at share>1E-12 or >1E-9,
  only dark cells are masked (114/764 on 82P). Bright bad cells
  remain, so no verdict changes. At >1E-6, 69 (82P) and 18 (81P)
  bad cells remain, still rejected. Masking does not convert 81/82.
D separate formal and physical verdicts: architecturally the
  safest (section 7). It cannot yet label 81/82 NEGLIGIBLE, because
  rate and population sensitivity is unmeasured and the fixed-point
  error is up to 24% on contributing cells.

---------------------------------------------------------------
7. CERTIFICATION ARCHITECTURE
---------------------------------------------------------------

knop_phase225/certification_architecture.md. The current chain is
ray IXIT -> KN_CS_LS -> FG_RF_ALL_CONV (all FULL) -> GS_RF_ALLCONV
-> LAM_FIELD_NOT_READY -> CONV_EXIT_MET -> knop_spec_gate. A
physical verdict can be inserted as a separate flag consumed only
at the LAM_FIELD_NOT_READY fold-in (solve_for_pops.f ~3511-3527).
KN_CS_* and FG_RF_ALL_CONV would keep their meaning, and both flags
would be logged. This keeps the five levels distinct (ray solve,
refresh, field adequacy, populations, spectrum). Not implemented.

---------------------------------------------------------------
8. ANSWERS
---------------------------------------------------------------

1. Max |d| cell: ray 82P, j=2, l=18731, 1033.3 A, I=4.81E-6,
   d=-3.68E-8 (corrected replay). The Phase 2.24 value of 1.9E-8
   came from the corrupted replay.
2. Large: 0.77% per pass on a column-peak cell. The fixed-point
   error on such cells reaches 1-24%.
3. Not measured. Decks exist only for 81/82.
4. Integrated: small per pass (R_INT ~1E-6), but not negligible
   after the convergence factor, and the bad set holds 63% of the
   82P E_INT.
5. Rate terms: the H I-weighted defect is about 10x larger than
   the intensity-weighted one. Real rates are not computed. Not
   shown negligible.
6. Yes in principle (JQW to J to rates, depths 1-5). Not measured;
   needs free matter.
7. Only through populations. The observer re-solves the transfer.
   No A/B was run.
8. Not quantified. The rows 1-2 lockstep (96-98% same sign) fits
   the mechanism, but that is not an attribution.
9. Keep Policy A for formal certification. Adopt the Policy D
   architecture only after rate and population sensitivity is
   measured. The data rule out B and C as currently framed.
10. Yes, by a separate flag at the LAM_FIELD_NOT_READY fold-in,
    without changing KN_CS_FULL or FG_RF_ALL_CONV.
11. P2DAMP: no. There is still no period-2 mode. AA: the omega=1
    replay contracts (rho about 0.93) and then stagnates near a
    1E-11/pass floor. AA could shorten the transient, so it may
    still be worth measuring as a separate experiment, but no
    physical benefit has been shown.
12. Mixed-sign validation needs a genuine inversion model with
    seeded cells on both sides of GAM=0, plus a replay with full
    fidelity (v2 decks).

---------------------------------------------------------------
9. CHANGES MADE (diagnostic driver only)
---------------------------------------------------------------

unit_tests/f90/drv_knop_replay.f90:
- KN_NC core-index remap (bug fix); PNT_SRCE warning.
- The end-state dump also writes KR_DBG_FLDM.
- The version rejection still works (magic 20261199 rejected).
The knray_replay binary was relinked manually (provenance.txt).
Decks were never modified (mtimes unchanged).

---------------------------------------------------------------
10. LIMITATIONS AND OUTSTANDING ITEMS
---------------------------------------------------------------

- Replay field fidelity is about 3-4% (v1 deck; EXT
  reconstruction suspected). This needs a v2 deck, which needs a
  rebuilt exe and a short armed capture.
- No decks for probe 83 or the control rays: the per-cell
  comparison against converged rays is outstanding.
- Only refresh r0000 has been replayed (about 34 refreshes are
  appended).
- Real rate sensitivity (cross sections, line profiles),
  population sensitivity and observer A/B are not done.
- The DS-lag share is not quantified (the proposed diagnostic
  needs review).
- The reference state (hst 60) itself stagnates, so "fixed-point
  error" is measured against a converged-to-1E-11/pass iterate,
  not an exact fixed point.
