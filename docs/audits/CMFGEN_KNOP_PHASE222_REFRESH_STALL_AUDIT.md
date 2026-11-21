# Phase 2.22 - Production KNOP refresh stall audit

Model: `m15_2p82d_KNOP` (all-compressing ejecta; v(R) monotone
48.9 -> 7073.6 km/s; no mixed-sign inversion).

Primary evidence run: `m15_2p82d_KNOP_A2_MODE1` (job 13414398,
frozen pre-flag exe `93364cac`, `KN_MAT_FREEZE=T` - populations
snapshotted and rolled back every call by construction, so the
field refresh is isolated from population feedback).

Diagnostic arms:

- `m15_2p82d_KNOP_diagcap` (job 13429272, small partition,
  `KN_PSWEEP_DIAG=T`, `KN_PSWEEP_LS=-1`): queued behind
  `QOSGrpCpuLimit`; has not produced `KN_PSWEEP_DIAG.dat`.
- `m15_2p82d_KNOP_diagcap_i` (job 13430220, interactive partition,
  same config via `run_cmfgen_SINGLE_FINAL.slurm`): RUNNING and
  producing `KN_PSWEEP_DIAG.dat`. Cold-start LTE arm on the same
  donor FGINEXT_GS, so its first refresh starts from the mature
  store state: its refresh-1 census (`residual max=6.4698E-02`,
  `capped=3`, `stall=3`, `plat=1`) is bit-identical to A2's first
  recorded refresh, and refreshes 2-3 settle to `capped=2`,
  matching A2's steady state. The arm is therefore a
  deterministic replay of the A2 trajectory and its per-ray
  records answer the capped-ray question directly.

Status: audit in progress. Source locations below are line numbers
in `plane/fg_j_cmf_nonmono_v1.f` at baseline commit `e28833fe`
unless marked otherwise; the running executables were built at
frozen commit `217049ee` - the audited lines are identical between
the two (verified on the sweeps/march/store blocks).

---------------------------------------------------------------
1. THE REFRESH-LEVEL OBSERVED STATE (A2, ~135 refreshes)
---------------------------------------------------------------

From OUTGEN (iteration 159 refresh block):

    residual max= 6.4698E-02  (tol 2.0E-03)
    ray exits: tol=44 bulk=37 stall=3 cap=0 div=0 plat=0
      unsolved=1 (legitimate NI=1 slot)
    bulk rays: maxres=2.618E-02 floor=6.0E-02
    stall/cap rays: maxres=6.963E-01 bad cells=6521
    commit verdicts: committed=81 capped=2 jump-rej=0
      nc-rej=0 div=0 unsolved=1 worst=(ixit 3,v 7)
    retained-field delta: max=5.0735E-01
    store cells>50%: 1
    unconverged cells=19841 of 628950981 (frac 3.15e-5)
    hot frequencies: 43259 in 982 contiguous runs
    matched/compress=17174408/258714053 ; max offset=158

Over 122+ refreshes the residual is stationary
(median ~6.5E-2, range 3.4E-2 to 7.5E-2; one sample below 5%).
GS_RF_DCHG (committed-generation drift) sits at 4.0-4.5E-2:
the store is a near-fixed point while the per-ray residual
plateaus. Both are consistent with the residual being set by
the same small set of cells every refresh.

---------------------------------------------------------------
1a. TASK 1 RESOLVED: THE TWO CAPPED PRODUCTION RAYS ARE 81P/82P
---------------------------------------------------------------

PSD record format (`plane/fg_j_cmf_nonmono_v1.f` ~12408-12424,
write site ~12903): one wrapped line per pass,
`PSD ray ch ref NPO_final IXIT` followed by 18 columns
(npo, omg0, omega, res, dmax, drms, dabs, jmx, lmx, nbad, ntot,
nv, nflr, nnf, ncp, nng, sign(GAM@jmx), NFM_L@jmx). ref is
LEDG_N_RF (great-iteration refresh counter; the three KN-CEN
refreshes captured so far all landed inside ref 0).

Per-ray exit census over the first 3 refreshes (2292+ PSD
records, all 84 solved rays):

  refresh 1 (census capped=3, matching A2 refresh 1):
    ray 81 P  ixit=3  npo=14  res=2.56E-2
    ray 82 P  ixit=3  npo=15  res=6.47E-2   <- = A2 residual max
    ray 83 P  ixit=3  npo=26  res=6.96E-1   <- probe, suppressed
    ray 12 M  ixit=6  npo=10  res=1.90E-2   <- plateau exit
  refresh 2+ (census capped=2, matching A2 steady state):
    ray 81 P  ixit=3  res=2.36E-2 .. 2.22E-2 (npo 17,18,27)
    ray 82 P  ixit=3  res=6.63E-2 .. 7.22E-2 (npo 21,22)
    ray 83 P  ixit=3  res=6.96E-1 .. 6.29E-1 (npo 26,30)
    ray 12 M  ixit=2  (bulk-exited; transient only)
  every other (ray,ch): ixit 1-2. M channels of 81/82/83
  converge in 2-3 passes at res 1e-4..1e-7.

The two capped production rays in the mature state are
therefore ray 81 and ray 82, outgoing (P) channel only.
The KN-CEN flagged set {22,45,50 ch=1} is a DISJOINT
population: those rays exit ixit=1 every refresh; their cells
are merely the largest RELATIVE committed movers (floor-noise
cells, section 2d). Offender set != capped set. A2's own
residual max line is ray 82's stall residual
(6.4698E-02 identical on refresh 1; range 6.3-7.2E-2 over
refreshes = A2's observed 3.4-7.5E-2 plateau band).

Stall-ray cell signature (last-pass records):

  ray   res     dmax    dabs    jmx lmx   nbad  ntot    nflr gsgn nfm
  81P   2.2E-2  1.0E-2  1.0E-30  2  1373    730  263621  9344  +1   1
  82P   6.5E-2  2.9E-2  2.9E-26  2  1370   5932  214003  7521  +1   1
  83P*  6.3E-1  2.8E-1  6.9E-28  2  1635  87145  162930  5700  +1   1
  12M   1.9E-2  1.9E-2  1.6E-17 109 1655    302   84178 1.99M   -1   1

  *probe. ncp=nng=0 on all final passes (no in-sweep cap or
  negative-floor events at exit). Argmax cells wander within
  the floor population between refreshes (lmx 1368..14403)
  but always sit at jmx=2 - the boundary-adjacent inner cell
  on the outward channel, GAM>0 (redder-donor branch), NFM=1.
  Absolute defects at the argmax are 1e-9..1e-30: the ~2-7%
  residuals are floored-denominator relative noise on dark
  cells, not a physical field error. Ray 12M's plateau is the
  same phenomenon one level deeper: ntot collapsed 5.8M->84k
  as cells fell below the significance floor, and the stuck
  residual 1.87e-2 alternates 1.831/1.865/1.896e-2 - a
  period-2/3 limit cycle in the metric, killed by the
  plateau detector after 10 passes.

Omega trajectories on stall rays ratchet to the 0.45 floor
(1.0 -> 0.7 -> 0.805 -> 0.9257 -> 0.648 -> 0.4536 -> 0.45),
with residuals oscillating non-monotonically - the
alternating-sweep period-2 signature predicted in section 4.

---------------------------------------------------------------
2. FLOOR-DOMINATED CELLS: EXACT SEMANTICS (Task 2)
---------------------------------------------------------------

Definitions audited in MARCH (13508+) and SWEEP_FLD_ORD (11738+):

a. Update equation (inward branch, MARCH ~13540):

    VNEW = FLD(I,L)*A0L(I)
         + SP_L(I)*A1L(I) + SP_L(I+1)*A2L(I)
         + DS_L(I)*A3L(I) + DS_L(I+1)*A4L(I)

   then, in order:

    i.  VNEW <= 0 or NaN  ->  VNEW = 1.0E-30     (the "floor";
        NNEG/NNG_P count the events)
    ii. VNEW > NONMONO_INEXT_CAP -> capped to CAP (NCAP/NCP_P)
    iii.DCELL = |VNEW-VOLD| / MAX(VOLD, 1.0E-4*VMAX_L)
        where VMAX_L = MAX(MAXVAL(FLD(1:NI,L)), 1.0E-25)
        (per-frequency column peak of the LIVE solve field,
         floored at 1.0E-25)
    iv. RES  = MAX(RES, DCELL)       <-- OUTSIDE the gate
        RESL = MAX(RESL, DCELL)      <-- OUTSIDE the gate
    v.  IF(MAX(VOLD,|VNEW|) > 1.0E-4*VMAX_L):
            NTOT_CUR++, and
            IF(DCELL > KN_RAY_TOL) NBAD_CUR++
    vi. FLD(I+1,L) = VOLD + OMEGA_RAY*(VNEW-VOLD)   (SOR write)

b. floor-dom flag (KN-Q diagnostic, ~12362):
       floor-dom = (VOLD < 1.0E-4*MAX(VMAXF,1.0E-25))

Answers to the task's ten questions:

1. Variable floored: VNEW (the freshly computed march value),
   not VOLD.
2. Floor value: 1.0E-30.
3. Applied to VNEW only (before the relaxation write).
4. BEFORE relaxation: the floor replaces VNEW, then
   FLD = VOLD + omega*(VNEW-VOLD) relaxes toward the floored
   value.
5. Floor cells enter NBAD only when they pass the significance
   gate (MAX(VOLD,|VNEW|) > 1e-4*VMAX_L); a pure floor cell
   (both below the gate) never enters NBAD.
6. Same for NTOT: gated out.
7. RES: YES - DCELL contributes to RES unconditionally
   (RES=MAX(RES,DCELL) precedes the gate).
8. R_MAP (KN_PSWEEP per-pass map defect): same normalized
   form, SCALE=MAX(|X|,1E-4*VMAXF(L)) (comment ~1447):
   floor cells contribute through the floored denominator.
9. Stored I_NEXT: WRITE_STORE (14288+) re-sanitizes V1:
   non-positive/NaN -> 1.0E-30, >CAP -> CAP, then writes
   unconditionally. WRITE_STORE_CAP (~14420) does the same
   sanitize, then additionally bounds the committed step:
      DVB = FG_COMMIT_NC * MAX(V2, 1.0E-4*VMAX_F(LS2))
      IF(|V1-V2| > DVB) V1 = V2 + SIGN(DVB, V1-V2)
   with FG_COMMIT_NC=0.05, FG_COMMIT_JUMP=50.0.
10. Donor source: SP_L is built from FLD AFTER previous cells'
    writes, so floored values DO feed downstream sources
    within the same pass (in-place Gauss-Seidel in depth;
    frequency donors come from whichever pass state exists at
    SRC_PRIME time - see section 4).

c. Direct cell-level evidence (q argmax record, refresh ~135):

    q argmax j,l= 2 1282 vold,vnew,vmaxraw,vmaxf,denom=
      8.74361E-30 1.02835E-28 9.24789E-30 1.00000E-25 1.00000E-29
      floor-dom=T

   The cell's column peak vmaxraw=9.25E-30 is BELOW the 1.0E-25
   metric floor, so the denominator is 1.0E-4*1.0E-25=1.0E-29
   and the measured q=9.41 = |1.03E-28-8.7E-30|/1.0E-29.
   In absolute terms the defect is ~9.4E-29 - physically zero
   (physical J peaks are ~1e-5). The q-metric's whole tail is
   generated by the floored denominator on dark columns:

    q sig ch=1: sigmax=1.96E-2 absmax=4.11E-4 flrmax=6.65E-4
    q sig ch=2: sigmax=9.41    absmax=1.45E-4 flrmax=0.382

   i.e. the largest ABSOLUTE map defect over all 1.26e9 cells
   is ~1.5E-4; the relative tail (up to 9.4) lives entirely in
   cells whose columns are at/below the metric floor.

d. The flagged-cell population (KN-CEN, stable from refresh 2):

    ray 22 ch=1 (jp,l2) = (2..5, 490)
    ray 45 ch=1 (jp,l2) = (2..5, 2068)
    ray 50 ch=1 (jp,l2) = (2..5, 1040)
    ray 83 ch=2 (jp,l2) = (1..2, 1248..1249)   [probe; keep out
                                              of prod stats]

   All production flagged cells: v1=1.0E-30 (floored fresh),
   v2~1e-34..1e-38 (retained, below the write floor), vms=
   column peak ~1.0e-28, vr~0.4994, rat~1e4..1e7.
   Geometry check: ch=1 flagged cells sit at jp=2..5, i.e. the
   points immediately inside the outer (vacuum) boundary of the
   incoming channel - where the incoming intensity physically
   vanishes.  Compact-depth map qmap2 confirms k=2 ->
   r=1.25303E5 (R(2), outer region). Column peaks ~1e-28 mean
   the entire frequency column is dark, not just the flagged
   cells.

e. FGINEXT_GS parse (rec 9/10 = I_P/I_M_DEPTH_NEXT, dims
   N_COMP=70, NP=85, NCF=56377, tag=295):

    ray 22 slot 489 : cells jp2-5 = 7.87..7.92E-33 ; colmax 1.06E-28
    ray 45 slot 2067: cells jp2-5 = 1.500..1.517E-32; colmax 1.10E-28
    ray 50 slot 1039: cells jp2-5 = 7.25..7.39E-30 ; colmax 1.20E-28
    ray 83 I_P 1247 : cells ~3.57E-28 ; colmax 3.57E-28
    ray 83 I_P 1248 : cells ~5.02E-28 ; colmax 5.02E-28

   The committed store agrees with the KN-CEN records; the
   cells sit at 1e-33..1e-30 inside columns peaking at ~1e-28.

---------------------------------------------------------------
3. WHICH CELLS HOLD THE RESIDUAL PLATEAU (Task 2 cont.)
---------------------------------------------------------------

The 6.5e-2 plateau and q=9.4 are both floor-normalized
artifacts: |DeltaV| ~ 1e-28 to 1e-31 on cells whose columns
peak at ~1e-28 (or below the 1e-25 metric floor), divided by
denominators floored at 1e-4*VMAX. The production flagged
cells contribute vr~0.5 against the 50x jump-veto scale
(FG_COMMIT_JUMP=50), not against KN_RAY_TOL.

The per-pass PSD capture (job 13430220, section 1a) now
confirms it directly: on stall rays 81P/82P/83P the RES argmax
sits at jmx=2 with |Delta V| ~ 1e-9..1e-30 - the same
boundary-adjacent floor-cell class as the KN-CEN flagged
cells, on rays disjoint from the flagged set. The
absolute-defect bound is established twice over: the map is
converged to <=1.5E-4 in intensity units over every cell.

---------------------------------------------------------------
4. THE ACTUAL REFRESH OPERATOR (Task 3)
---------------------------------------------------------------

IMPORTANT CORRECTION (26-Nov-2026): the earlier version of this
section described the unordered SWEEP_FLD path (ASC/DESC
alternating per pass). All 84 production rays run ORD_SWEEP=T -
the ordered march SWEEP_FLD_ORD - which does NOT alternate. The
unordered path remains relevant only as the context of the
ray-12M period-2 observation; it is not the production operator.

Structure of one ordered pass (SWEEP_FLD_ORD, ~11738+;
KNOP_RAY_SOLVE_V1 10039+ calls it twice per refresh - once
inward for the I_M field, once outward for I_P):

    pack build (once per solve, pass-invariant):
      DO L=1,NCF: BUILD_FREQ(L) -> PK_DTAU,PK_A0..PK_A4,
      PK_SR,PK_Q (NI x NCF), PK_ETANI/CHINI per L

    per pass NPO:
      rows marched in characteristic direction
        (P channel: J = NI-1 .. 1, outward from the IM_NI
         boundary row; M channel: J = 2 .. NI inward)
      frequency direction chosen PER ROW by GAM sign:
        GAM(J)>=0 (expanding):   L = 1..NCF   (ascending)
        GAM(J)<0  (compressing): L = NCF..1   (descending)
      helpers per row: KN_ORD_SPW_FILL (ring buffer of this
        row's source values), KN_ORD_SEC (secants between ring
        entries), KN_ORD_DS (monotonic-parabola slope),
        KN_ORD_SPV_IMPL (donor interpolation)
      march update per cell:
        VNEW = FLD(I,L)*A0 + SP_L*A1 + SP_L(I+1)*A2
             + DS_L*A3 + DS_L(I+1)*A4   (P-channel form)
        FLD(I,L) = VOLD + OMEGA_RAY*(VNEW - VOLD)

Remaining lagged dependencies (the operator is NOT a single
triangular solve):

  - DS corner: DS(J,L) reads SEC(J-1,L), built from SPW(J-1),
    which KN_ORD_SPW_FILL filled with the PRE-march (one-pass-
    lagged) source of the row beyond the marching front. Loop:
    (J,L)_n <- DS(J,L) <- SPV(J-1,L) <- FLD(J-1,L-1)_{n-1} - a
    genuine Jacobi-style one-row lag, the main candidate for the
    slow/stall mode on rays 81/82.
  - compressing-cell interpolant donors on the not-yet-marched
    side of the front (JPR1/JPR2/WINX at L+NFM).
  - MATCH_FREQ: NFM_L(I) (default 1, up to ~158 bins) scales the
    compressing upwind slot by DRF_L/KN_DLN_N; with the row-wise
    DESC order the donor is usually current-within-pass, but the
    corner through the next row is lagged as above.

Nonsmooth branches: VNEW floor (1e-30), VNEW cap
(NONMONO_INEXT_CAP), SLOPES monotonicity clamps, MATCH_FREQ
NFM cap - the operator is piecewise affine; the floor branch
keeps dark cells in the iteration at negligible absolute level.

Controller (both sweeps): RATIO_OM growth damping (x0.7, floor
0.45), recovery x1.15 when consecutive-pass ratio < 0.5, NOSC
sign-flip damping, NSTALL>=2 / plateau exits (IXIT=3/6), bulk
exit IXIT=2 (subset criterion + KN_RAY_BULK_TOL floor + optional
FG_BULK_CONTR contraction gate), pass cap IXIT=4, hot-mode
continuation (TRY_HOT) before stall/cap exits fire.

Unordered path (SWEEP_FLD, not production here): alternates
ASC/DESC whole-ray frequency direction per pass; the composed
2-pass map supports period-2 orbits (measured on ray 12M:
residuals alternating 0.0187/0.0183 across passes 4-10, at
dabs=1.6e-17 absolute amplitude). NOSC exists because such
lambda ~ -1 modes were observed on that path.

Measured on A2/cold replay: rho_p90~0.95, rho_max~0.99-1.18
(OUTGEN 'omega est') - the composed refresh map is not uniformly
contracting on the worst rays; whether the ordered march's
one-row lag alone reproduces the stall is the Task-4 question
(replay driver below).

---------------------------------------------------------------
5. COMMIT LOGIC (Task 6)
---------------------------------------------------------------

Verdict tree (FG_COMMIT_CONV block, ~10688-10790):

    JUMP > FG_COMMIT_JUMP(50)      -> verdict 1, reject
    not KN_CS_FULL (any non-converged exit class):
        FG_COMMIT_NC>0 -> WRITE_STORE_CAP, verdict 3
        else           -> verdict 2, reject
    KN_CS_FULL (tol-converged)     -> WRITE_STORE, verdict 0
    probe ray (mode 1)             -> verdict 7, suppressed

Consequences:

- A capped (verdict 3) ray keeps its stall residual in
  RES_LS -> the maturity gate cannot certify it converged.
  Fail-closed by construction.
- WRITE_STORE commits the sanitized fresh field directly;
  WRITE_STORE_CAP bounds each cell's committed step by
  DVB=0.05*max(V2,1e-4*VMAX_F). On floor-dominated cells
  (V2~1e-34, VMAX_F~1e-28) DVB~5e-34, so the committed value
  drifts on a ~1e-34/refresh ratchet - the observed v2 ladder
  1e-38..1e-34 across jp=2..5 is the ratchet mid-flight.
- GS_RF_DCHG is measured on the committed generation
  (retained-field delta line), distinct from the per-ray RES.
  The 4.0-4.5E-2 drift is real committed motion, dominated by
  the same floor-cell population.

---------------------------------------------------------------
6. WORKING DIAGNOSIS (pending per-ray PSD confirmation)
---------------------------------------------------------------

The plateau has two separable pieces:

(a) The metric: RES and the independent q defect are both
    floored-denominator relative measures; their tails sit on
    cells whose fields are ~1e-28 and below - dark/vacuum
    inflow cells. In absolute units the map is converged
    (absmax <= 1.5E-4). These cells are excluded from NBAD/
    NTOT but unconditionally included in RES.

(b) The dynamics: rays 81P and 82P stall-exit (IXIT=3) and
    take capped commits every refresh (resolved, section 1a -
    the capped pair is NOT inside the KN-CEN offender set
    {22,45,50}; the two populations are disjoint). The stall
    argmax sits at j=2, GAM>0, NFM_L=1 cells - the inner
    boundary-adjacent floor population on the outward channel.
    On the evidence of the alternating ASC/DESC sweep +
    SOR(omega->0.45) + floor/cap branches, and of the
    oscillating residual sequences and the ray-12M period-2/3
    metric limit cycle, the working hypothesis is a
    period-2 / slow-mode orbit confined to floor-dominated
    cells: the operator is converged in absolute units
    (dabs <= ~1e-9 on the worst cell) while the floored
    relative metric cannot certify it. That would make this a
    METRIC/controller failure riding on a physically dark
    cell set, not an operator failure - the pending test is
    whether an absolute-denominator variant of the residual
    certifies the same rays (controlled experiment, Task 4).

---------------------------------------------------------------
7. OPEN ITEMS / PENDING DATA
---------------------------------------------------------------

- Job 13430220 (interactive, diagcap_i): PRODUCING DATA - 3
  KN-CEN refreshes + per-ray PSD captured; the capped pair is
  identified (81P, 82P; section 1a). Continues running to
  accumulate later-refresh passes for the stall persistence
  check and Task-4 seed extraction.
- Job 13429272 (diagcap on small): still QOSGrpCpuLimit.
- A1/A2 (13414391/98), B1/B2 (13427111/15 pending), obs arms
  (13427285/6): monitored; B1/B2 never started (same quota).
- Task 4 experiment design (implemented 26-Nov-2026):
  KN_RDUMP_LS/KN_RDUMP_LS2 (VADAT, default -1 = off) arm a
  replay-deck writer inside KNOP_RAY_SOLVE_V1. Phase 1
  (pre-P-sweep) serializes every module input the ordered
  solve consumes: KN_ND/KN_ND_EXT/NCF/KN_NC/N_COMP/ND_ADD,
  KN_RAY_* controller constants, R_EXT/LOG_R_EXT/KN_R/KN_P,
  Z/R_RAY/RAY_PNT/J_PNT/COMP_DEPTH_IDX, GAM/GAM_MAG/dGAMdR,
  KN_CHI_COEF/KN_ETA_COEF, KN_FREQ/KN_DLN_P/KN_DLN_N/KN_IC/
  KN_DBB, KN_INNER_METH, the post-M field, both
  I_*_DEPTH_NEXT store columns and the seeded P field; flag
  bits record FREQ_MATCH/P2DAMP/ORD_SWEEP/AITKEN/AA/probe/
  BULK_CONTR/ADAPT_OMEGA_V2/PLATEAU_DET plus (deck v2) the
  KN_CHI_EXT/KN_ETA_EXT presence bit, FG_SEED_COMMIT,
  FG_RAY_SEEDED(LS), KN_Q_RESID, KN_CERT_CONV and the
  JMP_SIB_M/P presence bit. Deck v2 (magic 20261123) also
  serializes KN_CHI_EXT/KN_ETA_EXT (extension-grid tables the
  on-node cells read in BUILD_FREQ) and JMP_SIB_M/P (the
  per-frequency sibling-column references the commit jump
  veto compares against); v1 decks are rejected by the
  driver on the magic field. Phase 2 appends the
  production outcome (RES_P/NP_P/IXIT_P/OMEGA_RAY + final
  iterate) as the verification record. The offline driver
  unit_tests/f90/drv_knop_replay.f90 (+ runner
  unit_tests/python/run_knop_replay.py) pokes the deck into
  the same module globals and calls the production
  KNOP_RAY_SOLVE_V1 unchanged - the replay runs the real
  ordered march and controller, not a reimplementation.
  Modes: base (must reproduce the production record),
  loop N (one-pass calls with store write-back = caller-side
  pass loop), om/btol/stall/maxp overrides, seed cold|scale,
  verify (unrelaxed pass = independent map defect).
  Armed run: clone m15_2p82d_KNOP_diagdump_i (interactive,
  repo dev exe), config_resolved.yaml carries
  KN_RDUMP_LS=81, KN_RDUMP_LS2=82.
