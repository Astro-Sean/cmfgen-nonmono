!
! Data module for FG_J_CMF. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE FG_J_CMF_NONMONO_MOD_V1
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! The *_STORE routines are used to store the radiation field, as computed
! on the previous call to FG_J_CMF.
! The *_PREV routines are used to store the radiation field, as computed
! for the previous frequency.
!
! The *_PREV routines are updated from the *_STORE routines when NEW_FREQ
! is .TRUE.
!
!***************************************************************************
!***************************************************************************
!
! Data variables required by both routines.
!
! Dimensionded NRAY_MAX,NP
!
	REAL(KIND=LDP), ALLOCATABLE :: R_RAY(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: Z(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: DTAU(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: AV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CV(:,:)
!
! Dimensiond ND_EXT,4
!
	REAL(KIND=LDP), ALLOCATABLE :: V_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_COEF(:,:)
!
	INTEGER, ALLOCATABLE :: J_PNT(:,:)
	INTEGER, ALLOCATABLE :: H_PNT(:,:)
!
! Dimensioned ND+ND_ADD=ND_EXT
!
	REAL(KIND=LDP), ALLOCATABLE :: R_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_R_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: Z_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_ETA_EXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_CHI_EXT(:)
!
! Dimensioned NRAY_MAX,NP
!
	INTEGER, ALLOCATABLE :: RAY_PNT(:,:)
!
	INTEGER, ALLOCATABLE :: NI_RAY(:)
	INTEGER, ALLOCATABLE :: MAX_LS(:)
!
! Dimensioned NRAY_MAX
!
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR(:)
	REAL(KIND=LDP), ALLOCATABLE :: Q(:)
	REAL(KIND=LDP), ALLOCATABLE :: QH(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_RAY(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE_RAY(:)
!
! Non-monotonic velocity-field support.
!
! SIGN(1,GAM(I,LS)) identifies the local frequency-flow direction:
!   +1  => expanding flow  (T2 >= 0, standard Sobolev/CMF resonance)
!   -1  => compressing flow (T2 < 0, reverse shock / cold dense shell)
! GAM_MAG(I,LS) = |V/R| * |1 + SIGMA*MU^2| = |GAM| (always >= 0).
! Use ABS(T1)*ABS(T2), not T1*ABS(T2): T1 carries the sign of V, so
! inflow (V < 0) would otherwise make GAM_MAG negative and corrupt chi_hat.
! GAM_MAG is used to construct the generalized opacity chi_hat (see below)
! and for resonance counting diagnostics.
! The signed GAM is retained as GAM(I,LS) and is used in the explicit
! source term Q = (GAM/dLOG_NU)/chi_hat per Knop et al. 2009
! (Paper I: A&A 496, 295; Paper II: A&A 501, 813).
! chi_hat = chi + |GAM|/dLOG_NU (the generalized opacity) is kept strictly
! positive; the velocity-coupling term is moved into the source
! S' = SOURCE + Q*I_upwind.  CMFGEN solves the frequency (log-nu) form of
! the transfer equation.  Knop's 4*a*I is the I_lambda proportional term;
! its I_nu analogue (+3*a_nu*I, or +4*a_nu*I under a d(nu*I)/dnu split)
! does NOT vanish -- only the invariant intensity I_nu/nu^3 is free of it.
! It is omitted here by the same convention as stock V13/V11: the ray
! equation carries only the GAM*dI/d(ln nu) drift term, and the dilution
! physics enters at the moment level via the advection terms (MOM_J_DDT)
! or relativistic coefficients (MOM_JREL).  No S_tilde source is needed.
! Audit note: the 3a*I_nu term cancels identically between dI/ds and
! a*dI/dln(nu) in the I_nu form, so nothing is dropped; the piece of
! Knop's lambda-form 4*a*I that survives the split is the implicit
! (same-bin) part of his upwind frequency difference, which here is
! the |GAM|/DLN absorbed into chi_hat.  Sign reversal under
! compression is handled by the GAM<0 upwind switch, not by the term.
! edited by Sean Brennan
!
	REAL(KIND=LDP), ALLOCATABLE :: GAM_MAG(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: dGAMdR(:,:)
!
! N_NEG_SIGMA : count of depths with Sigma < -1 (diagnostic only).
! N_NEG_INT   : count of intensity points clamped to >= 0 (diagnostic).
! N_NAN_JNU   : count of NaN/negative JNU_STORE points clamped (diagnostic).
! edited by Sean Brennan
!
	INTEGER :: N_NEG_SIGMA=0
	INTEGER :: N_NEG_INT=0
	INTEGER :: N_NAN_JNU=0
!
! N_NEG_INT and N_NAN_JNU are reset on every call (i.e. every frequency),
! so on their own they are invisible: nothing reported them and the sweep
! finished with the counters holding only the last frequency's value.
! The *_SWEEP totals below accumulate them across a whole frequency sweep
! and are reported at the start of the next sweep.  A non-zero
! N_NEG_INT_SWEEP means the short-characteristics solution produced
! negative intensities that were clamped to 1.0E-30, which is the direct
! observable for an inconsistent source term in the compression zones.
! edited by Sean Brennan
!
	INTEGER :: N_NEG_INT_SWEEP=0
	INTEGER :: N_NAN_JNU_SWEEP=0
	INTEGER :: N_SWEEP_CNT=0
!
! V4 diagnostic (2026-09-09): Depth histogram of JNU clamping.
! JNU_CLAMP_DEPTH(I) counts how many times depth I had JNU_STORE clamped
! across the current sweep.  Written to JNU_CLAMP_DIAG.txt at the end
! of the first sweep of the first great iteration, then cleared.
! edited by Sean Brennan
	INTEGER, ALLOCATABLE, SAVE :: JNU_CLAMP_DEPTH(:)
	INTEGER, SAVE :: JNU_CLAMP_DIAG_WRITTEN=0
	INTEGER, SAVE :: JNU_CLAMP_FREQ_CNT=0
	INTEGER, PARAMETER :: JNU_CLAMP_FREQ_MAX=10000
	REAL(KIND=LDP), SAVE :: JNU_CLAMP_FREQ_VALS(10000)
	REAL(KIND=LDP), SAVE :: JNU_CLAMP_FREQ_LAMS(10000)
!
! Phase 2 first-occurrence diagnostic: prints the FIRST NaN/negative
! value encountered in the formal solver, with full context.
! FIRST_FAIL_FOUND is .TRUE. after the first print, suppressing further
! output.  Protected by OMP CRITICAL in the parallel ray loop.
! edited by Sean Brennan (Phase 2 audit)
!
	LOGICAL, SAVE :: FIRST_FAIL_FOUND=.FALSE.
	INTEGER, SAVE :: FIRST_FAIL_TYPE=0
!
! DIAG_PRINTED : ensures the per-ray resonance diagnostics are printed
! only once per CMFGEN run (on the first frequency), not on every INIT.
!
! edited by Sean Brennan
	LOGICAL, SAVE :: DIAG_PRINTED=.FALSE.
!
! INEXT_MAP_PRINTED : ensures the I_NEXT depth-mapping diagnostic table
! (ray point -> RAY_PNT -> K2 -> COMP_DEPTH_IDX) is printed once per
! grid definition, on the first frequency of a sweep.  Re-armed whenever
! NEW_R_GRID rebuilds the ray geometry (audit diagnostic, 15-Sep-2026).
!
	LOGICAL, SAVE :: INEXT_MAP_PRINTED=.FALSE.
!
! SIGMA_ZERO_WARNED : ensures the all-zero SIGMA_TRUE_IN notice is
! printed once per process, not once per frequency.
!
	LOGICAL, SAVE :: SIGMA_ZERO_WARNED=.FALSE.
!
! Multi-root resonance summary counters (printed once after the LS loop).
!
! edited by Sean Brennan
	INTEGER, SAVE :: N_RAYS_MULTI=0
	INTEGER, SAVE :: NRES_MAX=0
	INTEGER, SAVE :: NRES_MIN=0
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to the DIFFERENCE equation approach for  solving the
! transfer equation.
!	
	REAL(KIND=LDP), ALLOCATABLE :: AV_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: AV_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CV_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CV_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: PAR_AV(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: GAMH(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TA(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TB(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TC(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: GB(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: H(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: XM(:)
	REAL(KIND=LDP), ALLOCATABLE :: U(:)
	REAL(KIND=LDP), ALLOCATABLE :: VB(:)
	REAL(KIND=LDP), ALLOCATABLE :: VC(:)
	REAL(KIND=LDP), ALLOCATABLE :: DIV(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: OLDCHI(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLDCHI_STORE(:)
!
!***************************************************************************
!***************************************************************************
!
! Variables specific to short characteristic approach of solving the
! transfer equation.
!
! Dimensioned NRAY_MAX,NP
!
	REAL(KIND=LDP), ALLOCATABLE :: I_P_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_P_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M_PREV(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M_STORE(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: A0(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A1(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A2(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A3(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: A4(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: I_P(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M(:,:)
!
! Full Knop et al. 2009 sign-dependent upwind frequency coupling.
!
! In non-monotonic flows, compressing points (GAM < 0) couple to the
! NEXT frequency (lower frequency, higher FREQ_INDX), not the previous
! one.  Since CMFGEN sweeps frequencies in decreasing order, the next
! frequency has not yet been computed in the current great iteration.
! We use a Gauss-Seidel approach: I_NEXT comes from the PREVIOUS great
! iteration, stored at depth grid points for all frequencies.
!
! OPTIMIZED (A+B): Sparse depth storage + in-place write to NEXT.
!
! Only depths where SIGMA_TRUE < -1 can ever have GAM < 0 (for any ray
! angle).  We store only those N_COMP depths (Optimization A), using
! COMP_DEPTH_IDX(K) to map full-depth index K -> sparse index (1..N_COMP)
! or 0 if the depth is not compressing.
!
! The frequency sweep increments FREQ_IDX_ITER = 0, 1, ..., NCF-1 (0
! during the first frequency; storage slot s = FREQ_IDX_ITER+1).  At step i
! we READ NEXT[i+1] and WRITE NEXT[i].  Since NEXT[i] was already read
! at step i-1, we can write directly to NEXT, eliminating the separate
! CURR arrays and the INIT copy/zero (Optimization B).
!
! Memory: 2 * N_COMP * NP * NCF_IN * 8 bytes (was 4 * ND * NP * NCF_IN * 8).
!
! edited by Sean Brennan
	REAL(KIND=LDP), ALLOCATABLE :: I_P_DEPTH_NEXT(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: I_M_DEPTH_NEXT(:,:,:)
!
! COMP_DEPTH_IDX(K): sparse index for depth K (1..N_COMP) or 0 if not compressing.
! N_COMP: number of depths with SIGMA_TRUE < -1.
!
! edited by Sean Brennan
	INTEGER, ALLOCATABLE :: COMP_DEPTH_IDX(:)
	INTEGER, SAVE :: N_COMP=0
!
! Compression-set hysteresis (edited by Devin, 29-Sep-2026): the
! NEW_R_GRID test below rebuilds the entire ray geometry -- and with
! it the I_P/I_M_DEPTH_NEXT coupling store, which is reallocated
! zeroed -- whenever the set of compressing depths changes.  A depth
! whose SIGMA_TRUE hovers near -1 (or whose V hovers near 0) can then
! flap the set on each NG-acceleration/grid revision and cold-start
! the Gauss-Seidel field repeatedly, which the SOLVE_FOR_POPS maturity
! gate reads as a permanently immature field (forced LAMBDA forever).
! The mask is therefore sticky on exit only: a depth enters the set on
! the same SIGMA_TRUE < -1 / V < 0 test as before, but once flagged it
! leaves the set only when SIGMA_TRUE rises decisively above -1 (by
! COMP_SIG_HYST) or V rises decisively above zero (COMP_V_HYST, cm/s).
! The exit side is the only physically negligible direction: at
! SIGMA_TRUE > -0.98 the compressing term |1+sigma*mu^2| needs mu >
! ~0.99 to be negative at all, so which side of the flag a marginal
! depth sits on changes no ray source term meaningfully.  The stored
! mask is still built with the plain thresholds below, so the
! FGINEXT_GS fingerprint compare (bitwise COMP_DEPTH_IDX match) is
! unaffected.
!
	REAL(KIND=LDP), PARAMETER :: COMP_SIG_HYST=2.0E-2_LDP
	REAL(KIND=LDP), PARAMETER :: COMP_V_HYST=1.0E+5_LDP
!
! FREQ_IDX_ITER counts frequencies within the current great iteration.
! Reset to 0 at each INIT, incremented at each NEW_FREQ.
!
! edited by Sean Brennan
	INTEGER, SAVE :: FREQ_IDX_ITER=0
!
! HAS_COMPRESSING: .TRUE. if any input SIGMA_TRUE_IN < -1 (compression zones).
! When .FALSE., the large (ND,NP,NCF) Gauss-Seidel intensity arrays are NOT
! allocated (saving up to ~40 GB for ND=300), and all points use the
! standard I_PREV upwind path (identical to monotonic behaviour).
!
! edited by Sean Brennan
	LOGICAL, SAVE :: HAS_COMPRESSING=.FALSE.
!
! FREEZE_I_NEXT: when .TRUE., the in-place writes to I_P_DEPTH_NEXT /
! I_M_DEPTH_NEXT at the end of the formal solve are skipped, so an
! outer solver's finite-difference residual evaluations see the same
! outer Gauss-Seidel I_NEXT state.  The reads of I_NEXT for
! compressing points are unaffected, so perturbed residuals still
! feel the frozen coupling.  Set by JFNK_SOLVE_CORR_V1 around its
! coupled residual/JVP eval sweeps ([JFNK_COUPLED]=T); the default
! path leaves it .FALSE.
!
! edited by Sean Brennan (V4, 2026-10)
	LOGICAL, SAVE :: FREEZE_I_NEXT=.FALSE.
!
! Gauss-Seidel restart/checkpoint state (added 16-Sep-2026).
!
! I_P/I_M_DEPTH_NEXT are process-local (SAVE'd module memory) and were
! previously lost on every SLURM requeue, forcing the N_LAM_RWARM
! LAMBDA warm-up in SOLVE_FOR_POPS while I_NEXT re-converged
! from zero.  They are now checkpointed to the unformatted file
! FGINEXT_GS by DUMP_FG_INEXT_GS, which is called from SCR_RITE_V2 at
! every population checkpoint (always a complete GS generation, since
! the frequency sweep finishes before the population solve).
!
! I_NEXT_GS_OK   : .TRUE. only when RESTORE_FG_INEXT_GS loaded a
!                  validated complete generation in THIS process.
! FG_GS_TAG_EXP  : the NITSF value SCR_READ_V2 restored at startup
!                  (set by SET_FG_GS_TAG from CMFGEN_SUB); the file tag
!                  must equal it, tying the GS state to the same
!                  checkpoint event as SCRTEMP.
! FG_GS_NP/NCF   : NP and NCF_IN at allocation time (saved because the
!                  dump routine, called from SCR_RITE_V2, cannot see
!                  the routine arguments).
!
! FG_GS_MAXCHG  : in-sweep store drift until the per-ray refresh runs,
!                  then the max unrelaxed refresh residual normalized by
!                  max(old intensity, 1e-4 of the frequency-column peak).
!                  Reset at INIT; when no refresh runs, the in-sweep
!                  store-change metric remains as the fallback.
! FG_INEXT_CAPPED: count of I_NEXT writes this sweep that hit the
!                  NONMONO_INEXT_CAP clamp or the non-finite floor.
! FG_GS_CHG_FILE / FG_GS_CAP_FILE: the refresh max residual and capped
!                  write count recorded in the restored checkpoint file
!                  (-1 = no file).  SOLVE_FOR_POPS uses them only to
!                  decide whether to skip the process warm-up; the fresh
!                  bulk census still controls every live population step.
!                  The saved max must be below GS_MATURE_TOL and the
!                  saved cap count must be zero.
!
! edited by Sean Brennan (16-Sep-2026)
	LOGICAL, SAVE :: I_NEXT_GS_OK=.FALSE.
	INTEGER, SAVE :: FG_GS_TAG_EXP=-1
	INTEGER, SAVE :: FG_GS_NP=-1
	INTEGER, SAVE :: FG_GS_NCF=-1
	REAL(KIND=LDP), SAVE :: FG_GS_MAXCHG=0.0_LDP
	INTEGER*8, SAVE :: FG_GS_NCONV=0
	INTEGER*8, SAVE :: FG_INEXT_CAPPED=0
!
! Same-census cap accounting (maturity-gate fix): the sweep-store
! loop and the refresh passes visit different cell sets, so a single
! FG_INEXT_CAPPED numerator against the refresh's FG_GS_NTOT
! denominator mis-rates the cap share.  FG_SW_CAPPED counts the
! sweep cap/NaN events only, FG_SW_CAPTOT counts the sweep-store
! writes evaluated (the sweep census), and FG_RF_CAPPED counts the
! refresh-pass events only.  FG_INEXT_CAPPED keeps accumulating the
! total so the FGINEXT_GS checkpoint record is unchanged.
!
	INTEGER*8, SAVE :: FG_SW_CAPPED=0
	INTEGER*8, SAVE :: FG_SW_CAPTOT=0
	INTEGER*8, SAVE :: FG_RF_CAPPED=0
!
! Compressing-cell sink census (edited by Devin, 18-Oct-2026): a ray
! point with GAM < 0 whose bracketing depths are both non-compressing
! (or which runs while HAS_COMPRESSING is false) keeps
! chi_hat = chi + |GAM|/DLN but drops the |Q|*I_NEXT source term --
! the cell absorbs without the compensating upwind feed.  Counted per
! ray in FG_NOSLOT_RAY (LS-indexed so the OMP ray loop is race-free),
! summed into FG_NOSLOT_CNT once per frequency, and read once per
! great iteration by GET_FG_NOSLOT for the KN_NOSLOT census line.
!
	INTEGER*8, SAVE :: FG_NOSLOT_CNT=0
	INTEGER, ALLOCATABLE :: FG_NOSLOT_RAY(:)
!
! FG_GS_NBAD / FG_GS_NTOT: bulk convergence accounting from the last
! FG_KNOP_REFRESH_V1 pass -- number of significant cells (intensity
! above 1e-4 of the column peak) whose unrelaxed fixed-point residual
! still exceeds KN_RAY_TOL, and the total significant-cell count.
! NTOT=0 means the refresh did not run this iteration, so gates must
! fall back to the max-norm alone.
!
	INTEGER*8, SAVE :: FG_GS_NBAD=-1
	INTEGER*8, SAVE :: FG_GS_NTOT=0
	REAL(KIND=LDP), SAVE :: FG_GS_CHG_FILE=-1.0_LDP
	INTEGER, SAVE :: FG_GS_CAP_FILE=-1
!
! W1-W9 (17-Oct-2026, KNOP maturity-chain repair):
!
! FG_SW_MAXCHG / FG_SW_NCONV: the in-sweep store drift, snapshotted
! by FG_KNOP_REFRESH_V1 before it overwrites FG_GS_MAXCHG/NCONV
! with the refresh fixed-point metrics.  The two semantics measure
! different things (sweep drift is what the populations experience;
! refresh residual measures lagged-map self-consistency), so the
! gate needs both -- MOM keeps the same split in MOM_SW_*.
!
! FG_GS_NDIV: rays whose refresh solve diverged (IXIT=5) on the
! last refresh.  Gated on ray share, not the RESID=1.0 pin the
! diverged ray used to leave in MAXVAL(RES_LS).
!
! FG_GS_NBAD_BLK / FG_GS_NTOT_BLK: best-pass bad-cell census.  For
! stall (IXIT=3) and pass-cap (IXIT=4) exits the last pass still
! carries above-tolerance cells that the gate reads as immature
! even though the ray's best pass was inside the bulk band, so the
! gate reads MIN(FG_GS_NBAD, FG_GS_NBAD_BLK) (W3).
!
! FG_DIV_HEAL: store slots zeroed by the diverged-ray heal (W2):
! a diverged ray skips its store write, so a poisoned slot value
! (non-finite or > FG_DIV_HEAL_CAP * the ray's stored column peak)
! would re-seed and re-diverge the next refresh (MOM heals the
! same class of slot each call).  Zeroed slots take the seeded
! floor path on the next refresh.
!
! FG_COMMIT_REJ / FG_COMMIT_JUMP: whole-ray commit rollback (W9).
! When the max |V1-V2|/max(V2,1e-4*peak) over a ray's refreshed
! stores exceeds FG_COMMIT_JUMP the writes are skipped entirely --
! bounds the committed-near-blowup class (audit: retained delta
! 4.19E-1 -> 5.03E3) while allowing genuine factor-of-few rewrites.
!
	REAL(KIND=LDP), SAVE :: FG_SW_MAXCHG=0.0_LDP
	INTEGER*8, SAVE :: FG_SW_NCONV=0
	INTEGER, SAVE :: FG_GS_NDIV=0
	INTEGER*8, SAVE :: FG_GS_NBAD_BLK=-1
	INTEGER*8, SAVE :: FG_GS_NTOT_BLK=0
!
! FG_GS_NSTAL / FG_GS_NPLAT split the bulk bad-cell census by exit
! class: stall/cap rays (IXIT=3,4) are still slow-converging, plateau
! rays (IXIT=6) stopped improving entirely.  The two feed different
! unsafe-share tiers in SOLVE_FOR_POPS.
!
	INTEGER*8, SAVE :: FG_GS_NSTAL=0
	INTEGER*8, SAVE :: FG_GS_NPLAT=0
	INTEGER, SAVE :: FG_DIV_HEAL=0
	INTEGER, SAVE :: FG_COMMIT_REJ=0
	REAL(KIND=LDP), SAVE :: FG_COMMIT_JUMP=50.0_LDP
	REAL(KIND=LDP), SAVE :: FG_DIV_HEAL_CAP=1.0E+6_LDP
!
! FG_COMMIT_CONV (VADAT [FG_COMMIT_CONV], default F): fail-closed
! refresh commits.  When T, a ray field solve that ended on a stall,
! pass-cap or plateau exit (IXIT 3/4/6 -- residual still above
! KN_RAY_TOL) does not reach WRITE_STORE: the ray's lagged store
! keeps its pre-refresh generation instead of taking the last
! iterate of a non-converged fixed-point solve.  Only the
! FG_COMMIT_JUMP bound guarded those exits, so a capped ray could
! commit a field tens of factors off the sane seed (w19 class:
! residual 0.465 at the cap with retained delta 21).  FG_NC_REJ
! counts these rejections separately from jump rejects.  IXIT 0
! (exact), 1 (tol) and 2 (bulk-maturity) still commit.
!
	LOGICAL, SAVE :: FG_COMMIT_CONV=.FALSE.
	INTEGER, SAVE :: FG_NC_REJ=0
!
! FG_COMMIT_NC (VADAT [FG_COMMIT_NC], default 0): per-cell bound on
! the committed change for a NON-converged ray under FG_COMMIT_CONV.
! FG_COMMIT_CONV alone skips the write on bad exits, but a bulk exit
! (IXIT=2) is classed converged even when the residual is far above
! KN_RAY_TOL, so the non-converged test keys on
! MAX(RES_M,RES_P)>KN_RAY_TOL as well as the exit code.  NC>0 writes
! a capped provisional field (each cell bounded to NC times its
! retained value) instead of skipping outright: the store still
! tracks a contracting field but cannot ratchet by more than NC per
! refresh.  NC=0 keeps the pure reject.  w19 committed bulk-exit
! deltas of 8-27x at residual ~0.47; the geometric ratchet was the
! luminosity amplifier.
! FG_NC_CAP counts capped commits; FG_RAY_VC records the per-ray
! verdict for the summary census and worst-ray line.
!
	REAL(KIND=LDP), SAVE :: FG_COMMIT_NC=0.0_LDP
	INTEGER, SAVE :: FG_NC_CAP=0
! Verdict codes: -1 unsolved, 0 committed, 1 jump reject,
! 2 non-converged reject, 3 capped commit, 4 diverged.
	INTEGER, ALLOCATABLE, SAVE :: FG_RAY_VC(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: FG_RAY_DCH(:)
!
! Per-ray worst-cell record: depth/freq index of the last-pass
! residual max (of the worse direction) and how many consecutive
! passes held the same cell.  A persistent cell means one unstable
! mode dominates; a moving max means different cells reset the
! residual each pass.  Filled for solved rays only (0 otherwise).
!
	INTEGER, ALLOCATABLE, SAVE :: FG_RAY_WI(:),FG_RAY_WL(:)
	INTEGER, ALLOCATABLE, SAVE :: FG_RAY_WS(:)
!
! FG_BULK_CONTR (VADAT [FG_BULK_CONTR], default F): require the
! same-direction residual to be contracting (RATIO_OM<=KN_RAY_STALL)
! before the bulk bad-fraction exit can fire.  Without it a ray
! whose residual is flat at 0.05-0.5 with a sub-0.1% bad fraction
! exits "bulk converged" at pass 6 -- before the stall/plateau
! classifiers (which need pass>=8) ever run -- and commits under
! the converged label.  w19: 49/84 rays took this path at residual
! 0.25-0.49 during the luminosity runaway.  Ordering fix; F keeps
! the old exit.
!
	LOGICAL, SAVE :: FG_BULK_CONTR=.FALSE.
!
! W7 (17-Oct-2026): extended per-ray pass budget for rays with
! contiguous hot-frequency runs (RES_FREQ>KN_RAY_TOL over a run).
! After the standard exits fire, a ray with a hot run continues
! restricted to the run's window +/- 2 frequencies for up to
! KN_RAY_HOTPASSES extra passes; KN_HOT_TOTAL_PASSES caps the
! refresh-wide spend so pathological farms cannot dominate.
! KN_HOT_SPENT is the running refresh-global counter.
!
	INTEGER, SAVE :: KN_RAY_HOTPASSES=120
	INTEGER, SAVE :: KN_HOT_TOTAL_PASSES=480
	INTEGER*8, SAVE :: KN_HOT_SPENT=0
!
! W5 (17-Oct-2026): FG_ADAPT_OMEGA_V2=T prices rho_est from the
! per-pass residual sequence of the current refresh (same-map
! contraction) instead of the cross-iteration KN_RF_RESID ratio,
! which conflates source drift with map contraction.
! FG_PLATEAU_DET=T arms the plateau (|ratio-1|<0.05 x3) and the
! combined-direction stall triggers in both sweep loops.
!
	LOGICAL, SAVE :: FG_ADAPT_OMEGA_V2=.TRUE.
	LOGICAL, SAVE :: FG_PLATEAU_DET=.TRUE.
!
! Boundary-closure audit counters.  HBC = HNU_AT_OB/JNU_STORE(1)
! blows up when the boundary J hits its 1D-30 floor, and the moment
! solver's outer boundary row consumes the ratio directly.  Track
! the worst |HBC|, the smallest boundary J, and the count of
! non-finite or negative closures per sweep; reported by
! DUMP_FG_INEXT_GS.
!
	REAL(KIND=LDP), SAVE :: FG_MAX_HBC=0.0_LDP
	REAL(KIND=LDP), SAVE :: FG_MIN_JBND=1.0E+30_LDP
	INTEGER*8, SAVE :: FG_N_BADHBC=0
!
! Boundary-intensity and identity counters.  FG_N_NEG_BND_I counts
! floored boundary intensities (I+ or I- < 0 or non-finite at the
! outer boundary); FG_N_BND_IDENT counts frequencies where the J/K
! boundary identities (JNU = 0.5*(JPLUS+JMIN)) fail, and
! FG_MAX_BND_RES keeps the worst identity residual.
!
	INTEGER*8, SAVE :: FG_N_NEG_BND_I=0
	INTEGER*8, SAVE :: FG_N_BND_IDENT=0
	REAL(KIND=LDP), SAVE :: FG_MAX_BND_RES=0.0_LDP
!
! OBS-grid restore state (added Oct-2026).  The observation-frame
! driver CMF_FLUX builds its own frequency grid (different NCF and
! spacing from the model grid that wrote FGINEXT_GS), so the checkpoint
! cannot be restored by position.  SET_FG_GS_OBS_GRID (called by
! CMF_FLUX_SUB_V5 before its frequency loop) stores a copy of the OBS
! grid here; the first-call init then dispatches to
! RESTORE_FG_INEXT_OBS, which remaps the stored I_NEXT field in
! frequency instead of cold-starting it.  The (KC,LS) indexing needs
! no remap: OBS reads NP from RVTJ and builds rays via the same IMPAR
! on the same R, so compressing depths and rays coincide with the
! model's.
!
	LOGICAL, SAVE :: GS_OBS_MODE=.FALSE.
	INTEGER, SAVE :: GS_NCF_OBS=0
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: GS_NU_OBS(:)
!
! Diagonal approximate Lambda operator (Olson-Auer-Buchler) extracted
! from the actual non-monotonic formal solution.
!
! For each ray LS and depth grid point I (ray node K=J_PNT(I,LS)):
!   dI_M(K)/dS'(K) = A2(K-1,LS)      (K >= 2; 0 at K=1)
!   dI_P(K)/dS'(K) = A2(K,LS)        (K <= NI-1)
!                  = dI_M/dS' at K=NI when I_P(NI)=I_M(NI) (non-core
!                    ray or ZERO_FLUX/PNT_SRCE inner boundary)
!                  = 0 when I_P(NI)=I_CORE
! and AV(K,LS)=0.5*(I_P+I_M), so dAV/dS' = 0.5*(dI_P/dS' + dI_M/dS').
! With S' = ETA/chi_hat + Q*I_upwind (Q proportional to 1/chi_hat):
!   dS'/dETA = 1/chi_hat   and   dS'/dCHI = -S'/chi_hat
! where chi_hat = CHI_RAY (= chi + |GAM|/dLOG_NU, the opacity actually
! used to build DTAU on this ray).
!
! ALO_DJDETA(K) = sum_rays JQW * dAV/dS' / chi_hat      (dJ/dETA diag)
! ALO_DJDCHI(K) = sum_rays -JQW * dAV/dS' * S' / chi_hat (dJ/dCHI diag)
! LAMSTAR_DIAG(K) = sum_rays JQW * dAV/dS'  (dimensionless Lambda* diag)
!
! Off-diagonal terms (A0 propagation) and the non-differentiable dS
! (Hermite-derivative) terms are neglected.  Arrays are zeroed on every
! call (each frequency); ALO_VALID is set .TRUE. only after a completed
! INTEGRAL-method accumulation.
!
! edited by Sean Brennan
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: ALO_DJDETA(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: ALO_DJDCHI(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: LAMSTAR_DIAG(:)
	LOGICAL, SAVE :: ALO_VALID=.FALSE.
	INTEGER, SAVE :: ALO_ND=0
!
! Per-ray captures of OMP-PRIVATE quantities needed by the serial ALO
! accumulation: chi_hat, and the direction-dependent modified source
! S'_M = S + Q*I_M^upwind (inward) / S'_P = S + Q*I_P^upwind (outward).
! Written inside the parallel ray loop (each thread owns its LS column),
! read in the serial J_LOOP afterwards.
!
	REAL(KIND=LDP), ALLOCATABLE :: ALO_CHIHAT_RAY(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ALO_SPR_M_RAY(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ALO_SPR_P_RAY(:,:)
!
! edited by Sean Brennan
! Dimensioned NRAY_MAX
!
	REAL(KIND=LDP), ALLOCATABLE :: EE(:)
	REAL(KIND=LDP), ALLOCATABLE :: E0(:)
	REAL(KIND=LDP), ALLOCATABLE :: E1(:)
	REAL(KIND=LDP), ALLOCATABLE :: E2(:)
	REAL(KIND=LDP), ALLOCATABLE :: E3(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE_PRIME(:)
	REAL(KIND=LDP), ALLOCATABLE :: S(:)
	REAL(KIND=LDP), ALLOCATABLE :: dS(:)
!
! Knop et al. 2009 (Paper I: A&A 496, 295; Paper II: A&A 501, 813)
! S_tilde REMOVED (07-Sep-2026): Knop's S_tilde carries the I_lambda
! proportional term 4*a_lambda*I.  The removed code transplanted that
! wavelength-form term into the frequency-form scheme where its
! coefficient is wrong.  See the comment above the formal solution for
! the full justification.
!
	REAL(KIND=LDP) PREVIOUS_FREQ
	REAL(KIND=LDP) VDOP_FRAC_SAV
!
	INTEGER ND_EXT
	INTEGER ND_ADD
	INTEGER NP_MAX
	INTEGER NRAY_MAX
	INTEGER IDMIN,IDMAX
!
	CHARACTER(LEN=10) OLD_SOLUTION_OPTIONS
	CHARACTER(LEN=10) SOLUTION_METHOD
	LOGICAL INSERT
	LOGICAL FIRST_TIME
	LOGICAL NEW_R_GRID
	LOGICAL WRITE_IP
	LOGICAL DIF_OR_ZF
!
	DATA VDOP_FRAC_SAV/-1001.0_LDP/ 	!Absurd value.
	DATA FIRST_TIME/.TRUE./
	DATA PREVIOUS_FREQ/0.0_LDP/
!
! Ceiling applied when I_P/I_M_DEPTH_NEXT are stored for the next
! great iteration's compressing-zone source feedback (SOURCE_PRIME
! picks up ABS(Q)*I_NEXT).  Bounds the Gauss-Seidel frequency
! coupling against runaway intensities produced by badly
! unconverged population states.  Read from VADAT
! (NONMONO_INEXT_CAP); default 1.0D20 leaves physical
! intensities untouched while stopping Inf/huge feedback.
!
	REAL(KIND=LDP) NONMONO_INEXT_CAP
	DATA NONMONO_INEXT_CAP/1.0E+20_LDP/
!
! Under-relaxation factor for the lagged I_NEXT Gauss-Seidel write.
! The compressing-branch upwind field is fed back into the source on
! the next great iteration; on models with strong compression the
! undamped update (straight replacement) can drive a bounded limit
! cycle or amplify the SE<->RT loop.  Writing
!   I_NEXT = OLD + INEXT_RELAX*(NEW - OLD)
! damps the fixed-point map while leaving the fully-converged fixed
! point unchanged (NEW=OLD is still stationary).  Cells with no
! prior value (OLD<=0) are filled directly.
!
	REAL(KIND=LDP) INEXT_RELAX
	DATA INEXT_RELAX/0.5_LDP/
!
! KN_ADAPT_RELAX (VADAT [KN_ADAPT_RELAX], default T, set by
! SET_FG_KN_ADAPT): replace the fixed write-back weight with a
! spectral-radius estimate of the outer (sweep+refresh) fixed-point
! map.  rho_est is the ratio of successive refresh residual maxima
! (the dominant-mode contraction factor of the lagged map); the
! weight is omega = 1/(1+rho_est), so a non-contracting or
! sign-flipping dominant mode (rho_est >= 1) is damped to <= 1/2
! while a contracting map recovers toward the full update.
! KN_RF_RESID_PRV carries the previous call's residual maximum.
!
	LOGICAL KN_ADAPT_RELAX
	DATA KN_ADAPT_RELAX/.TRUE./
	REAL(KIND=LDP) KN_RF_RESID_PRV
	DATA KN_RF_RESID_PRV/-1.0_LDP/
!
! ----------------------------------------------------------------------
! Per-ray Knop refresh store (added Oct-2026).
!
! The original lagged design advanced the compressing-branch coupling
! one frequency bin per great iteration: at frequency k the solver
! read I_NEXT(k+1) written during the PREVIOUS sweep.  With ~4e4
! frequency bins and |GAM|/dLOG_NU >> CHI inside the interaction
! shell, the red->blue information chain requires ~NCF great
! iterations to mature -- this is the structural reason the KNOP
! models never reached the 5% GS maturity gate.
!
! Knop, Hauschildt & Baron (2009b, Paper II, Sect. 3.2) instead
! converge the wavelength boundary-value problem inside the formal
! solution for a fixed source function, ordering the solution to
! follow the information flow (blue->red where GAM>=0, red->blue
! where GAM<0).  This block implements that scheme as a post-sweep
! per-ray re-solve: during the main frequency sweep the solver stores
! (per frequency) the extended-grid opacity/emissivity, the
! scattering term, the frequencies and the sign-dependent dLOG_NU
! values.  After the sweep, FG_KNOP_REFRESH_V1 re-solves each ray's
! full (space,frequency) intensity field with alternating-direction
! sweeps until the coupled field is converged, and writes the result
! back into I_P_DEPTH_NEXT/I_M_DEPTH_NEXT (the same in-place store),
! so the next great iteration sees a converged upwind field.
!
! FG_NU_NEXT: nu(k+1) for the current sweep call, set per frequency by
! the caller (SET_FG_NU_NEXT) so compressing cells can use the correct
! asymmetric spacing ln(nu_k/nu_{k+1}).
!
	REAL(KIND=LDP) FG_NU_NEXT
!
! FG_FREQ_MATCH: VADAT DO_NONMONO_FREQ_MATCH (default F, set by
! SET_FG_FREQ_MATCH).  At compressing ray points the lagged upwind
! bin is rematched from the adjacent slot to the slot containing
! the characteristic foot, n = |GAM|*dZ/DLN bins along the ray, and
! the frequency derivative spans the matched log range.  Once the
! per-cell comoving drift exceeds one bin the +1 stencil
! extrapolates the derivative -- the oscillatory-eigenvalue regime
! of Ceulemans et al. (2025, JQSRT 343).  NM_RAY/DLN_RAY are
! OMP-private per-ray work arrays, filled each frequency call;
! defaults (offset 1, DLN_N_SAV) reproduce the pre-match stencil.
!
	LOGICAL FG_FREQ_MATCH
	DATA FG_FREQ_MATCH/.FALSE./
	INTEGER, ALLOCATABLE :: NM_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: DLN_RAY(:)
!
! DO_KNOP_REFRESH: VADAT-controlled switch (default TRUE when
! DO_NONMONO_KNOP).  Set by SET_FG_KNOP_REFRESH from
! RD_CONTROL_VARIABLES.  FALSE reproduces the original lagged
! behaviour for A/B testing.
!
	LOGICAL DO_KNOP_REFRESH
!
! PRINT_KNOP_DEBUG: VADAT-controlled switch ([PRINT_KNOP_DEBUG],
! default FALSE) set by SET_FG_KNOP_DEBUG.  Gates the verbose
! diagnostic tables (e.g. the I_NEXT depth-mapping dump, O(NP*NI)
! lines) that flood OUTGEN on fine ray grids.  The checks behind
! those tables still run; anomalies print an error regardless.
!
	LOGICAL PRINT_KNOP_DEBUG
!
! Captured per-frequency sweep inputs (allocated when the compressing
! store is built; sized by the actual ND_EXT/NCF_IN of the model).
! KN_CHI_EXT/KN_ETA_EXT are the pre-chi_hat opacity/emissivity on the
! extended grid.  Single precision (default REAL): the full-depth,
! full-frequency input store needs ~2*ND_EXT*NCF_IN*4 B which is
! ~300 MB for a typical model -- 4x cheaper than the equivalent
! double store and consistent with the requested memory target.
!
	REAL, ALLOCATABLE :: KN_CHI_EXT(:,:)
	REAL, ALLOCATABLE :: KN_ETA_EXT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_FREQ(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_DLN_P(:)  !ln(nu_{l-1}/nu_l) expanding (backward span)
	REAL(KIND=LDP), ALLOCATABLE :: KN_DLN_N(:)  !ln(nu_l/nu_{l+1}) compressing (forward span)
	REAL(KIND=LDP), ALLOCATABLE :: KN_DBB(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_IC(:)
!
! Monotone-cubic coefficients (log space, extended grid) for the
! captured opacities -- the same kernels the in-sweep solve applies
! via CHI_COEF/ETA_COEF.  Filled by FG_KNOP_REFRESH_V1 each call.
!
	REAL(KIND=LDP), ALLOCATABLE :: KN_CHI_COEF(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_ETA_COEF(:,:,:)
!
! Scalar state captured at grid-build time.
!
	INTEGER KN_ND,KN_NP,KN_NC,KN_NCF,KN_ND_EXT
	CHARACTER(LEN=12) KN_INNER_METH
	CHARACTER(LEN=6) KN_METHOD_SAV
	REAL(KIND=LDP), ALLOCATABLE :: KN_R(:),KN_P(:)
!
	LOGICAL KN_STORE_READY
	LOGICAL KN_RF_WARNED
!
! Convergence controls / diagnostics for the per-ray fixed-point
! iteration.  KN_RAY_TOL: maximum relative overwrite change per pass
! (values below 1e-25 carry no flux and are floored in the metric).
! KN_RAY_MAXPASS: alternating-direction pass cap per direction field.
! KN_MAX_THREADS: OMP thread cap for the refresh loop -- the per-ray
! work arrays scale as NI*NCF and each thread allocates its own copy,
! so this bounds refresh memory at ~KN_MAX_THREADS*NI*NCF*8 B.
!
	REAL(KIND=LDP) KN_RAY_TOL
	INTEGER KN_RAY_MAXPASS,KN_MAX_THREADS
	REAL(KIND=LDP) KN_RF_RESID
!
! FG_RF_DCHG: retained-field delta of the last refresh -- the
! max relative change between the in-sweep lagged store (the
! field the populations saw) and the field committed by
! WRITE_STORE, using the same per-frequency peak floor as the
! residual metrics.  FG-side analogue of MOM_RF_DCHG; the
! population gate reads it when deciding whether the refresh
! actually moved the field it will gate on.
!
	REAL(KIND=LDP) FG_RF_DCHG
!
! KN_RF_TOT_MB: total refresh memory budget (MB) across all OMP
! threads.  Each thread in KNOP_RAY_SOLVE_V1 can hold ~80 B per
! (depth,frequency) cell (8-array coefficient pack + FLD + work), so
! the thread count is reduced when threads*budget would exceed this.
!
	REAL(KIND=LDP) KN_RF_TOT_MB
!
! KN_RAY_MAX_INS: per-interval cap on the Doppler-driven ray-point
! insertion.  Unbounded INT(T1) can overflow the counter on a huge
! velocity jump and produce an arbitrarily long ray.
!
	INTEGER KN_RAY_MAX_INS
!
! KN_RAY_PACK_MB: per-thread memory cap (MB) on the pass-invariant
! coefficient pack used by KNOP_RAY_SOLVE_V1.  The pack is
! 8*NI*NCF*8 B per thread; above this cap (or on allocation failure)
! the per-pass BUILD_FREQ rebuild is used instead, preserving the
! previous behaviour exactly.
!
	REAL(KIND=LDP) KN_RAY_PACK_MB
!
! KN_RAY_OMEGA: under-relaxation applied to each overwrite inside the
! refresh marches.  The unrelaxed fixed-point map develops a ~-1
! eigenvalue limit cycle at isolated cells (observed: residual
! stalls at ~0.14, period-4 oscillation); damping the write leaves
! the fixed point unchanged but contracts the cycle.
!
	REAL(KIND=LDP) KN_RAY_OMEGA
!
! FG_ORD_SWEEP: VADAT [KN_RAY_ORDERED] (default T, set by
! SET_FG_ORD_SWEEP).  Selects the ordered refresh march of
! SWEEP_FLD_ORD: rows are marched along the characteristic in space
! and the per-row wavelength direction follows sign(GAM) -- the
! physically ordered Gauss-Seidel scheme of Knop et al. (2009,
! Paper II Sec. 2).  Only the spatially downwind slope corner
! (the lagged C terms of their Eqs. 31/32) and the deeper half of
! the compressing JPR interpolant stay lagged, instead of half of
! every strong coupling each alternating pass.  Requires the
! coefficient pack and NI>=3; otherwise SWEEP_FLD keeps running.
! F restores the alternating scheme for A/B testing.
!
	LOGICAL FG_ORD_SWEEP
!
! FG_ORD_AITKEN: VADAT [KN_ORD_AITKEN] (default F, set by
! SET_FG_AITKEN).  When T, SWEEP_FLD_ORD applies a guarded Aitken
! delta-squared extrapolation to the field iterate at each pass
! boundary, targeting the marginally contractive (|Q|->1) cells at
! strong compressions that otherwise stall the march against the
! bulk-maturity exit.  Restarted-Aitken convention: the extrapolated
! iterate replaces the pass result, so the next pass marches from
! the accelerated state.  Guards: significant cells only, |delta
! ratio|<1, positive and capped output; allocation failure disables
! the acceleration for that solve silently.
!
	LOGICAL FG_ORD_AITKEN
!
! FG_ORD_AA: VADAT [KN_ORD_AA] (default F, set by SET_FG_ORD_AA).
! Safeguarded type-II Anderson acceleration of the ordered pass
! iterate inside SWEEP_FLD_ORD, mirroring the runtime-tested MOM
! refresh AA (mom_j_cmf_nonmono_v1.f): residual-weighted Gram
! system, non-affine cell exclusion, Tikhonov ridge, ||corr||>||f||
! veto, 10|f| per-cell cap, positivity fallback.  FG_ORD_AA_M is
! the secant depth [KN_ORD_AA_M, default 2]; pairs kept = M+1.
! On top of the MOM guards a residual-acceptance rollback restores
! the pre-mix iterate when the next pass does not improve.
!
	LOGICAL FG_ORD_AA
	INTEGER FG_ORD_AA_M
!
! FG_ORD_P2DAMP: VADAT [KN_ORD_P2DAMP] (default F, set by
! SET_FG_ORD_P2DAMP).  Per-cell period-2 collapse inside
! SWEEP_FLD_ORD: after FG_ORD_P2N [KN_ORD_P2N, default 2]
! consecutive sign flips of (VNEW-VOLD) at a significant cell, the
! commit becomes the bracket midpoint -- sqrt(VOLD*VNEW) when both
! endpoints are positive, 0.5*(VOLD+VNEW) otherwise -- which is the
! exact fixed point of a pure multiplicative lambda=-1 mode and a
! bounded step inside the bracket in general.  Releases after
! FG_ORD_P2REL [KN_ORD_P2REL, default 3] flip-free passes.  The
! state costs 3 B per cell as int*1 arrays, allocated per ray
! solve, so the ~80 B/cell thread budget at KN_RF_TOT_MB is
! effectively unchanged.
!
	LOGICAL FG_ORD_P2DAMP
	INTEGER FG_ORD_P2N,FG_ORD_P2REL
!
! KN_RAY_STALL: VADAT [KN_RAY_STALL] (default 0.98, set by
! SET_FG_STALL_LIM).  Same-direction residual ratio above which a
! pass counts toward the NSTALL early exit.  The 0.90 setting
! quit rays still contracting at up to 10% per pass; those rays
! held nearly all the unconverged cells in the m15 2.82d cold
! validation (exit census), so the default was raised toward the
! flat-pass limit after the KN_RAY_STALL=0.98 arm cut the
! unconverged fraction ~40x at matched iterations.
!
	REAL(KIND=LDP) KN_RAY_STALL
!
! Per-refresh AA counters, accumulated under OMP CRITICAL in
! SWEEP_FLD_ORD and summarised by FG_KNOP_REFRESH_V1.
!
	INTEGER*8 KN_AA_APP,KN_AA_VETO,KN_AA_REJ,KN_AA_SNG
! Per-pass convergence history of the refresh (max over rays of the
! overwrite residual at each alternating-direction pass).  Diagnostic
! only; updated under OMP CRITICAL inside SWEEP_FLD.
!
	REAL(KIND=LDP), ALLOCATABLE :: KN_RES_PASS(:)
!
! Per-frequency count of (ray,field) pairs whose last completed pass
! still had an above-tolerance overwrite residual at that frequency.
! Diagnostic only; accumulated under OMP CRITICAL in
! KNOP_RAY_SOLVE_V1 and summarised by FG_KNOP_REFRESH_V1.  Used to
! locate the strongly-coupled frequency runs for the block-GS
! redesign (cmfgen_obsidian/audits/BLOCK_GS_DESIGN.md).
!
	INTEGER, ALLOCATABLE :: KN_HOT_CNT(:)
!
	SAVE
!
	DATA FG_NU_NEXT /-1.0_LDP/
	DATA DO_KNOP_REFRESH /.TRUE./
	DATA PRINT_KNOP_DEBUG /.FALSE./
	DATA KN_STORE_READY /.FALSE./
	DATA KN_RF_WARNED /.FALSE./
	DATA KN_RAY_TOL /2.0E-03_LDP/
	DATA KN_RAY_OMEGA /1.0_LDP/
	DATA FG_ORD_SWEEP /.TRUE./
	DATA FG_ORD_AITKEN /.FALSE./
	DATA FG_ORD_AA /.FALSE./
	DATA FG_ORD_AA_M /2/
	DATA FG_ORD_P2DAMP /.FALSE./
	DATA FG_ORD_P2N /2/,FG_ORD_P2REL /3/
	DATA KN_RAY_STALL /0.98_LDP/
	DATA KN_AA_APP /0_8/,KN_AA_VETO /0_8/,KN_AA_REJ /0_8/,
	1     KN_AA_SNG /0_8/
	DATA KN_RAY_MAXPASS /60/
	DATA KN_MAX_THREADS /16/
	DATA KN_RAY_PACK_MB /2048.0_LDP/
	DATA KN_RF_TOT_MB /32768.0_LDP/
	DATA KN_RAY_MAX_INS /256/
	DATA KN_ND /0/,KN_NP /0/,KN_NC /0/,KN_NCF /0/,KN_ND_EXT /0/
	DATA KN_RF_RESID /0.0_LDP/
	DATA FG_RF_DCHG /0.0_LDP/
!
	END MODULE FG_J_CMF_NONMONO_MOD_V1
!
! 
!
! Routine to compute the Eddington F, G and RSQN_ON_RSQJ Eddington factors
! for a single frequency. The transfer is done in the comoving-frame with
! first order frequency differencing. A DIFFERENCE or INTEGRAL equation
! approach can be used.
!
!The DIFFERENCE equation approach is at least 30% faster. However,
! the INTEGRAL equation approach is more stable. Because monotonic cubic
! interpolation is used for the Source function in the INTEGRAL equation
! approach, the intensities are guaranteed positive (or zero).
!
! The F, G, and RSQN_ON_RSQJ Eddington factors must be supplied.
!
! NB:
!	F = K / J
!	G=  N / H
!	RSQN_ON_RSQJ(I) = RSQ_N(I)/( RSQ_J(I)+ RQS_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' G is defined at all depths, and
!       RSQN_ON_RSQJ=0 at all depths.
!     IF N_TYPE='N_ON_J' RSQN_ON_RSQJ is defined at all
!       depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' one of G or RSQN_ON_RSQJ is
!       non-zero, and is the value to be used in MOM_J_CMF
!
! Routine also returns I+, so that observers flux can be computed. Note because
! the outer boundary can be thick, I+ is actually I+ - I- and hence is a flux
! like variable (This I+ = 2v at outer boundary).
!
! The radiation field at the previous (bluer) frequency is stored locally.
! Coupling to this frequency is controlled by the variable INIT. If INIT is
! true, the atmosphere is treated with the assumption that V=0.
!
! INIT must be TRUE on the very first call.
!
!
	SUBROUTINE FG_J_CMF_NONMONO_V1(ETA,CHI,ESEC,V,SIGMA,R,P,
	1                  JNU,FEDD,JQW,HQW,KQW,NQW,HMIDQW,NMIDQW,
	1                  RETURNED_IN_HBC,RETURNED_OUT_HBC,IPLUS_P,
	1                  FREQ,dLOG_NU,INNER_BND_METH,DBB,IC,
	1                  VDOP_VEC,VDOP_FRAC,REXT_FAC,
	1                  METHOD,SOLUTION_OPTIONS,
	1                  THK,INIT,NEW_FREQ,NC,NP,ND,SIGMA_TRUE_IN,NCF_IN)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1
	USE MOD_RAY_MOM_STORE
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_IN_EVAL,CONV_DIAG_LVL
	IMPLICIT NONE
!
! Altered 05-Oct-2026: V and SIGMA_TRUE fingerprints added to the NEW_R_GRID test so a
!                         velocity update at fixed R cannot leave the ray geometry on a
!                         stale field (edited by Devin).
! Altered 29-Sep-2026: One-sided hysteresis on the compressing-depth set used by the
!                         NEW_R_GRID test; in-loop bound checks on the J_PNT/H_PNT
!                         ray-bracket searches (edited by Devin).
! Altered 19-Aug-2015: Added check that SIGMA> -1.0D0; Some code improvements (cur_hmi, 28-Jun-2015).
! Altered 30-Dec-2013: Fixed expressions for HPLUS_OB and HPLUS_IB: Change OMP schedule to dynamic.
! Altered 08-Jan-2012: Changed to V12. Added REXT_FAC to call.
! Altered 07-Jul-2011: Changed extrapolation region I to max-value of 10 (rathee than ND/6).
!                         This was done for shell model. Mayu need further changes.
! Altered 21-Apr-2011: Bug fix for HNU_AT_OB when using DIFF option (not generally sed anyway).
! Altered 14-Jan-2010: Added JPLUS_IB, JMIN_IB etc to improve computation of
!                        boundary conditions.
! Altered 14-May-2009: Altered value of IDMAX (for ETA and CHI interpolaton).
!                        IDMAX & IDMIN now stored in FG_J_CMF_MOD_V11.
! Altered 19-Jan-2009: Changed IP_DATA to IP_FG_DATA which must exist for I(p)
!                        to be output.
! Created 17-Dec-2008: Based on FG_J_CMF_V10
!
	INTEGER NC,NP,ND
	INTEGER NCF_IN			!Total number of frequencies (for Knop storage)
	REAL(KIND=LDP) ETA(ND),CHI(ND),ESEC(ND)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND),P(NP)
	REAL(KIND=LDP) SIGMA_TRUE_IN(ND)	!Unclamped SIGMA for non-monotonic physics
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!     QW denotes quadrature weight.
!
	REAL(KIND=LDP) JQW(ND,NP)
	REAL(KIND=LDP) HQW(ND,NP)
	REAL(KIND=LDP) KQW(ND,NP)
	REAL(KIND=LDP) NQW(ND,NP)
!
	REAL(KIND=LDP) HMIDQW(ND,NP)
	REAL(KIND=LDP) NMIDQW(ND,NP)
!
	REAL(KIND=LDP) JNU(ND)
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) IPLUS_P(NP)
!
! VDOP_VEC(I) is the minimum DOPPLER width for all species at depth I. It will include
! both a turbulent, and and thermal contribution for the ionization species with the
! highest mass. VDOP_FRAC is used to set the minimum velocity step size along a ray.
!
	REAL(KIND=LDP) VDOP_VEC(ND)
	REAL(KIND=LDP) VDOP_FRAC
	REAL(KIND=LDP) REXT_FAC		!Factor to scale RMAX by if thick atmosphere.
	REAL(KIND=LDP) RETURNED_IN_HBC
	REAL(KIND=LDP) RETURNED_OUT_HBC
!
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) IC
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) dLOG_NU
!
	CHARACTER*(*) SOLUTION_OPTIONS
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE
	CHARACTER(LEN=*) INNER_BND_METH
!
! Use "Thick" boundary condition. at outer boundary. Only noted when INIT
! is true. All subsequent frequencies will use the same boundary condition
! independent of the passed value (Until INIT is set to TRUE again).
!
	LOGICAL THK
!
! First frequency -- no frequency coupling.

	LOGICAL INIT
!
! Upon leaving this routine the radiation field along each ray is stored. This
! will provide the blue wing information necessary for the next frequency.
! This routine may, however, be used in an iterative loop. In this case the
! "blue wing" information should remain unaltered between calls.
! NEW_FREQ indicates that a new_frequency is being passed, and hence the "blue
! wing" information should be updated.
!
	LOGICAL NEW_FREQ	
!
	INTEGER, PARAMETER :: ND_ADD_MAX=24
	INTEGER, PARAMETER :: NC_PNT_SRCE=2
!
	INTEGER, SAVE :: FREQ_CNT
	INTEGER ACCESS_F
        INTEGER LU_IP
	INTEGER IOS,REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC
	INTEGER N_MAP_BAD
!
! The following arrays do not need to be stored, and hence can be created
! each time.
!
	REAL(KIND=LDP) CV_BOUND(NP)		!Outer boundary V
	REAL(KIND=LDP) I_M_IN_BND(NP)		!Inner boundary
	REAL(KIND=LDP) IBOUND(NP)		!Incident intensity on outer boundary.
!
	INTEGER N_ERR_MAX,FG_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL(KIND=LDP) FG_ERR_ON_FREQ
	INTEGER FG_ERR_TYPE
	COMMON /FG_J_CMF_ERR/FG_ERR_ON_FREQ(N_ERR_MAX),
	1                    FG_ERR_TYPE(N_ERR_MAX),FG_ERR_CNT
	LOGICAL NEG_AV_VALUE
	LOGICAL COMP_REQUIRED
	LOGICAL L_IN_SET
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER KN_DUMP_NS
!
! Local variables.
!
	REAL(KIND=LDP), PARAMETER :: ONE=1
	INTEGER, PARAMETER :: NINS=4
!
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
!
	INTEGER NI_SMALL
	INTEGER I,J,K,K2,KC,LS
	INTEGER KC1,KC2
	INTEGER LK
	REAL(KIND=LDP) W_INEXT,INEXT_V
	REAL(KIND=LDP) DLN_C,DLN_N_SAV
	REAL(KIND=LDP) DENOM_IN,DENOM_SCL
	INTEGER NI
	INTEGER NIR_ALO
	INTEGER N_OUT_LAMSTAR
	INTEGER NP_TMP
	INTEGER NRES
	INTEGER NSIG_NEG,NSIG_NEG1
!
	REAL(KIND=LDP) DBC
	REAL(KIND=LDP) I_CORE
	REAL(KIND=LDP) T1,T2,T3,T4
	REAL(KIND=LDP) VMAX_IP(NP),VMAX_IM(NP)
	REAL(KIND=LDP) T1_DTAU
	REAL(KIND=LDP) TA2,TB2
	REAL(KIND=LDP) DELTA_Z
	REAL(KIND=LDP) ALPHA
	REAL(KIND=LDP) ESEC_POW
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) VINF
	REAL(KIND=LDP) RMAX,DEL_R_FAC
	REAL(KIND=LDP) MU,dZ,PSQ
	REAL(KIND=LDP) DEL_R
	REAL(KIND=LDP) LAM_MIN,LAM_MAX
!
! Diagonal-ALO temporaries (see module comments near ALO_DJDETA).
!
	REAL(KIND=LDP) ALO_DIM,ALO_DIP,ALO_DAV
	REAL(KIND=LDP) LAMSTAR_MIN,LAMSTAR_MAX
!
! Change the following statement to TRUE if running on a VECTOR machine.
!
	LOGICAL, PARAMETER :: VECTOR_MACHINE=.FALSE.
	LOGICAL, PARAMETER :: DEFINE_AT_MID_POINTS=.FALSE.          !TRUE.
!
!
!
! Allocate data for moments which will be used to construct the Eddington
! factors.
!
	LUER=ERROR_LU()
!
! V4 audit fix (2026-09-07): Reset per-call diagnostic counters so the
! user can identify which iteration produced negative/NaN intensities.
! N_NEG_SIGMA is cumulative (printed once), so it is not reset here.
!
! On the first frequency of a sweep, report what the previous sweep
! accumulated and then start fresh.  Reporting here (rather than per
! frequency) keeps OUTGEN small while still surfacing the counts.
! edited by Sean Brennan
!
	IF(INIT)THEN
	  IF(N_NEG_INT_SWEEP .GT. 0 .OR. N_NAN_JNU_SWEEP .GT. 0)THEN
	    WRITE(LUER,'(A)')'FG_J_CMF_NONMONO_V1: previous sweep clamps'
	    IF(N_NEG_INT_SWEEP .GT. 0)THEN
	      WRITE(LUER,'(A,I10,A,I6)')
	1         '    negative intensities=',N_NEG_INT_SWEEP,
	1         '  over N_freq=',N_SWEEP_CNT
	      WRITE(LUER,'(A)')'    (non-zero: source and generalized'//
	1         ' opacity disagree in compression zones)'
	    END IF
	    IF(N_NAN_JNU_SWEEP .GT. 0)
	1     WRITE(LUER,'(A,I10)')
	1         '    NaN/negative JNU points=',N_NAN_JNU_SWEEP
	  END IF
	  N_NEG_INT_SWEEP=0
	  N_NAN_JNU_SWEEP=0
	  N_SWEEP_CNT=0
	END IF
	N_NEG_INT=0
	N_NAN_JNU=0
	IF(INIT)THEN
	  FG_MAX_HBC=0.0_LDP
	  FG_MIN_JBND=1.0E+30_LDP
	  FG_N_BADHBC=0
	  FG_N_NEG_BND_I=0
	  FG_N_BND_IDENT=0
	  FG_MAX_BND_RES=0.0_LDP
	END IF
!
! Allocation/dimension guard for SIGMA_TRUE_IN.  The caller
! (COMP_J_BLANK) passes SIGMA_TRUE.  If the array is unallocated or
! smaller than ND, the unclamped sigma profile used by the
! non-monotonic chi_hat formalism would be invalid, causing
! out-of-bounds reads or garbage in compression-zone detection.
! edited by Sean Brennan (04-Sep-2026)
!
	IF(SIZE(SIGMA_TRUE_IN) .LT. ND)THEN
	  WRITE(LUER,*)
	1    'Error in FG_J_CMF_NONMONO_V1: SIGMA_TRUE_IN has size',
	1    SIZE(SIGMA_TRUE_IN),' but ND=',ND
	  WRITE(LUER,*)
	1    'The unclamped velocity-gradient profile is required for the',
	1    ' non-monotonic formal solver.'
	  WRITE(LUER,*)
	1    'Check that SIGMA_TRUE (or SIGMA_TRUE_EXT) is allocated and',
	1    ' refreshed by RECOMPUTE_SIGMA_TRUE / EXTEND_SIGMA_TRUE.'
	  STOP
	END IF
! Print once per process only: on a homologous model sigma=0 is the
! correct answer, and repeating the warning every frequency floods
! OUTGEN.  The single print still flags a genuinely unpopulated
! SIGMA_TRUE on a non-monotonic model.
!
	IF(MAXVAL(ABS(SIGMA_TRUE_IN(1:ND))) .LT. 1.0E-10_LDP .AND.
	1   .NOT. SIGMA_ZERO_WARNED)THEN
	  SIGMA_ZERO_WARNED=.TRUE.
	  WRITE(LUER,*)
	1    'Warning: SIGMA_TRUE_IN is all zeros in FG_J_CMF_NONMONO_V1.'
	  WRITE(LUER,*)
	1    'The velocity-gradient profile was not populated. Re-derive',
	1    ' from R and V via RECOMPUTE_SIGMA_TRUE before formal transfer.',
	1    ' (printed once; expected for a homologous/monotonic model)'
	END IF
!
! V4 diagnostic (2026-09-09): allocate depth histogram for JNU clamping
	IF(ALLOCATED(JNU_CLAMP_DEPTH))THEN
	  IF(SIZE(JNU_CLAMP_DEPTH) .NE. ND)DEALLOCATE(JNU_CLAMP_DEPTH)
	END IF
	IF(.NOT. ALLOCATED(JNU_CLAMP_DEPTH))THEN
	  ALLOCATE(JNU_CLAMP_DEPTH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)JNU_CLAMP_DEPTH=0
	END IF
!
! The MOD_RAY_MOM_STORE block was allocated once at the first call's
! ND and never resized, while ND_STORE recorded that size for
! GET_MOMS_NON_REL.  A later call with a different ND (grid revision)
! wrote/zeroed out of bounds and left stale radii in R_STORE.  Apply
! the same realloc-on-change pattern as the ALO block below.
!
	IF(ALLOCATED(JNU_STORE))THEN
	  IF(SIZE(JNU_STORE) .NE. ND)THEN
	    DEALLOCATE(JNU_STORE,HNU_STORE,KNU_STORE,NNU_STORE,
	1      R_STORE,GAM_REL_STORE,RMID_STORE,EXT_RMID_STORE,
	2      STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error deallocating JNU_STORE block in'//
	1        ' FG_J_CMF_NONMONO_V1: Status=',IOS
	      STOP
	    END IF
	  END IF
	END IF
	IF(.NOT. ALLOCATED(JNU_STORE))THEN
	  ND_STORE=ND
	  ALLOCATE (JNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (KNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NNU_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (R_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_REL_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RMID_STORE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EXT_RMID_STORE(ND+1),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating JNU_STORE block in FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP	
          END IF
	END IF
!
! Diagonal-ALO storage (Olson-Auer-Buchler operator from the formal
! solution).  Allocated at the current call's ND (ND or ND_EXT for the
! extended-grid call).  Zeroed on EVERY call (each frequency) so a stale
! operator can never leak between frequencies or iterations.  ALO_VALID
! is set .TRUE. only when the INTEGRAL ray sweep has completed.
! edited by Sean Brennan
!
	IF(.NOT. ALLOCATED(ALO_DJDETA))THEN
	  ALLOCATE(ALO_DJDETA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(ALO_DJDCHI(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(LAMSTAR_DIAG(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating ALO arrays in FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP
	  END IF
	  ALO_ND=ND
	ELSE IF(ALO_ND .NE. ND)THEN
	  DEALLOCATE(ALO_DJDETA,ALO_DJDCHI,LAMSTAR_DIAG,STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error deallocating ALO arrays in'//
	1      ' FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP
	  END IF
	  ALLOCATE(ALO_DJDETA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(ALO_DJDCHI(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(LAMSTAR_DIAG(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error re-allocating ALO arrays in FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP
	  END IF
	  ALO_ND=ND
	END IF
	ALO_DJDETA(1:ND)=0.0_LDP
	ALO_DJDCHI(1:ND)=0.0_LDP
	LAMSTAR_DIAG(1:ND)=0.0_LDP
	ALO_VALID=.FALSE.
	IF(ALLOCATED(ALO_CHIHAT_RAY))THEN
	  ALO_CHIHAT_RAY=0.0_LDP
	  ALO_SPR_M_RAY=0.0_LDP
	  ALO_SPR_P_RAY=0.0_LDP
	END IF
!
! V4 audit fix (2026-09-07): reset FREEZE_I_NEXT on INIT so a leftover
! .TRUE. cannot suppress I_NEXT writes in the next formal solve.
! edited by Devin, 05-Oct-2026: must NOT clear during a JFNK eval sweep.
! The driver arms FREEZE_I_NEXT once per activation (jfnk_solve_corr_v1.f)
! precisely so the perturbed-pop evals cannot advance the lagged store;
! clearing here on every INIT made the freeze a no-op -- each eval wrote
! I_NEXT at perturbed pops, so F(x) depended on eval order and the FD JVP
! inherited a drift term (~1E18 vs BA ~1E13).  JFNK_IN_EVAL is set for the
! whole eval callback (cmfgen_sub.f JFNK_RAD_SWEEP_V1); outside evals the
! original stale-freeze cleanup still applies.
!
	IF(INIT)THEN
	  IF(.NOT. JFNK_IN_EVAL)FREEZE_I_NEXT=.FALSE.
	  R_STORE(1:ND)=R(1:ND)
	  GAM_REL_STORE(1:ND)=1.0_LDP/SQRT(1.0_LDP-(V(1:ND)/2.99792458E+10_LDP)**2)
	  DO I=1,ND-1
	    RMID_STORE(I)=0.5_LDP*(R(I)+R(I+1))
	    EXT_RMID_STORE(I+1)=0.5_LDP*(R(I)+R(I+1))
	  END DO
	  EXT_RMID_STORE(1)=R(1); EXT_RMID_STORE(ND+1)=R(ND)
	END IF
!
! Check to see whether we have a new R grid, or the solution options have
! changd. This can only happen when INIT is TRUE.
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. .NOT. FIRST_TIME)THEN
	  DO I=1,ND
!
! Relative-tolerance grid test: an exact .NE. compare fires on
! harmless roundoff (e.g. an R grid reloaded through an ASCII
! restart), which rebuilds the whole ray geometry and re-maps the
! lagged radiation state for nothing.  1E-8 relative is far tighter
! than any physical grid change but immune to representation noise.
!
	    IF(ABS(R(I)-R_EXT(ND_ADD+I)) .GT. 1.0E-8_LDP*
	1      MAX(ABS(R(I)),ABS(R_EXT(ND_ADD+I)),1.0_LDP))THEN
	      NEW_R_GRID=.TRUE.
	      J=ERROR_LU()
	      WRITE(J,*)'Warning: Updating RGRID in FG_J_CMF_NONMONO_V1'
	      EXIT
	    END IF
	  END DO
	  IF(ABS(VDOP_FRAC-VDOP_FRAC_SAV) .GT.
	1       1.0E-12_LDP*MAX(1.0_LDP,ABS(VDOP_FRAC))
	1       .OR.  OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)NEW_R_GRID=.TRUE.
!
! V / SIGMA_TRUE fingerprints (05-Oct-2026).  The build block copies
! V and SIGMA_TRUE_IN into the main-grid segment of V_EXT/SIGMA_EXT,
! so those arrays are already the stored record -- no extra save
! arrays are needed.  A V update at fixed R (the DO_V_PARAMS wind
! revision is the known path) leaves V_COEF/SIGMA_COEF, GAM, chi_hat
! and the ray insertion spacing on the old field while neither the
! R test above nor the compression-set test below necessarily fires.
! MOM_J_CMF_NONMONO_V1 fingerprints V the same way (V_SM_SV).
! Per-cell relative test with a floor of 1.0 in the natural units
! (km/s and sigma): bit-level and ASCII-restart noise (~1E-7
! relative) does not fire, and a near-stagnation cell cannot trigger
! on roundoff.
!
	  IF(.NOT. NEW_R_GRID)THEN
	    DO I=1,ND
	      IF(ABS(V(I)-V_EXT(ND_ADD+I)) .GT. 1.0E-6_LDP*
	1        MAX(ABS(V(I)),ABS(V_EXT(ND_ADD+I)),1.0_LDP))THEN
	        NEW_R_GRID=.TRUE.
	        J=ERROR_LU()
	        WRITE(J,'(A,I5,A,2ES16.8)')
	1          'Warning: V drift at depth',I,
	2          ' -- rebuilding ray grid (V vs stored):',
	3          V(I),V_EXT(ND_ADD+I)
	        EXIT
	      END IF
	      IF(ABS(SIGMA_TRUE_IN(I)-SIGMA_EXT(ND_ADD+I)) .GT.
	1        1.0E-6_LDP*MAX(ABS(SIGMA_TRUE_IN(I)),
	2        ABS(SIGMA_EXT(ND_ADD+I)),1.0_LDP))THEN
	        NEW_R_GRID=.TRUE.
	        J=ERROR_LU()
	        WRITE(J,'(A,I5,A,2ES16.8)')
	1          'Warning: SIGMA_TRUE drift at depth',I,
	2          ' -- rebuilding ray grid (sigma vs stored):',
	3          SIGMA_TRUE_IN(I),SIGMA_EXT(ND_ADD+I)
	        EXIT
	      END IF
	    END DO
	  END IF
	  IF(.NOT. NEW_R_GRID)THEN
	    DO I=1,ND
	      COMP_REQUIRED=.FALSE.
	      DO J=MAX(1,I-1),MIN(ND,I+1)
!
! edited by Devin, 29-Sep-2026: one-sided deadband on the compression
! set.  Depths enter the set on the unchanged SIGMA_TRUE < -1 / V < 0
! test, but a depth already inside the set (COMP_DEPTH_IDX(J) > 0)
! leaves only on a decisive exit -- SIGMA_TRUE above -1+COMP_SIG_HYST
! or V above +COMP_V_HYST.  Without it, a depth sitting near sigma=-1
! or V=0 flaps the set every time NG acceleration or a grid revision
! perturbs V; each flap rebuilt the ray geometry and re-zeroed
! I_P/I_M_DEPTH_NEXT, so the GS field could never mature.
!
	        L_IN_SET=.FALSE.
	        IF(ALLOCATED(COMP_DEPTH_IDX))THEN
	          IF(COMP_DEPTH_IDX(J) .GT. 0)L_IN_SET=.TRUE.
	        END IF
	        IF(L_IN_SET)THEN
	          IF(SIGMA_TRUE_IN(J) .LT. -1.0_LDP+COMP_SIG_HYST .OR.
	1            V(J) .LT. COMP_V_HYST)COMP_REQUIRED=.TRUE.
	        ELSE
	          IF(SIGMA_TRUE_IN(J) .LT. -1.0_LDP .OR.
	1            V(J) .LT. 0.0_LDP)COMP_REQUIRED=.TRUE.
	        END IF
	      END DO
	      IF(ALLOCATED(COMP_DEPTH_IDX))THEN
	        IF(COMP_REQUIRED .NEQV. (COMP_DEPTH_IDX(I) .GT. 0))THEN
	          NEW_R_GRID=.TRUE.
	          J=ERROR_LU()
	          WRITE(J,'(A,I6,A)')'FG_J_CMF_NONMONO_V1: compression'//
	1            ' set changed at depth',I,' -- rebuilding ray grid'
	          EXIT
	        END IF
	      ELSE IF(COMP_REQUIRED)THEN
	        NEW_R_GRID=.TRUE.
	        J=ERROR_LU()
	        WRITE(J,'(A)')'FG_J_CMF_NONMONO_V1: compression set'//
	1          ' appeared -- rebuilding ray grid'
	        EXIT
	      END IF
	    END DO
	  END IF
	END IF
	
! Deallocate all allocated rays if we are using a diferent solution technique.
! This option will only be used when testing, since in CMFGEN we will always use
! the same atmospheric structure.
!
	IF( ALLOCATED(R_EXT) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R_EXT )
	  DEALLOCATE ( LOG_R_EXT )
	  DEALLOCATE ( V_EXT )
	  DEALLOCATE ( Z_EXT )
	  DEALLOCATE ( SIGMA_EXT )
	  DEALLOCATE ( ETA_EXT )
	  DEALLOCATE ( CHI_EXT )
	  IF(ALLOCATED(LOG_ETA_EXT))DEALLOCATE ( LOG_ETA_EXT )
	  IF(ALLOCATED(LOG_CHI_EXT))DEALLOCATE ( LOG_CHI_EXT )
!
	  IF(ALLOCATED(Z))DEALLOCATE ( Z )
	  IF(ALLOCATED(RAY_PNT))DEALLOCATE ( RAY_PNT )
	  IF(ALLOCATED(R_RAY))DEALLOCATE ( R_RAY )
	  IF(ALLOCATED(NI_RAY))DEALLOCATE ( NI_RAY )
	  IF(ALLOCATED(J_PNT))DEALLOCATE ( J_PNT )
	  IF(ALLOCATED(H_PNT))DEALLOCATE ( H_PNT )
	  IF(ALLOCATED(V_COEF))DEALLOCATE ( V_COEF )
	  IF(ALLOCATED(SIGMA_COEF))DEALLOCATE ( SIGMA_COEF )
	  IF(ALLOCATED(ETA_COEF))DEALLOCATE ( ETA_COEF )
	  IF(ALLOCATED(CHI_COEF))DEALLOCATE ( CHI_COEF )
	  IF(ALLOCATED(V_RAY))DEALLOCATE ( V_RAY )
	  IF(ALLOCATED(SIGMA_RAY))DEALLOCATE ( SIGMA_RAY )
	  IF(ALLOCATED(ETA_RAY))DEALLOCATE ( ETA_RAY )
	  IF(ALLOCATED(CHI_RAY))DEALLOCATE ( CHI_RAY )
	  IF(ALLOCATED(dCHIdR_RAY))DEALLOCATE ( dCHIdR_RAY )
	  IF(ALLOCATED(SOURCE_RAY))DEALLOCATE ( SOURCE_RAY )
	  IF(ALLOCATED(GAM_MAG))DEALLOCATE ( GAM_MAG )
	  IF(ALLOCATED(dGAMdR))DEALLOCATE ( dGAMdR )
	  IF(ALLOCATED(dCHIdR))DEALLOCATE ( dCHIdR )
	  IF(ALLOCATED(Q))DEALLOCATE ( Q )
!
! NM_RAY/DLN_RAY (GS I_NEXT slot-match store) sit in the unconditional
! allocation block but were missing here, so any grid rebuild died on
! re-ALLOCATE with STAT=151 (already allocated).  Their contents index
! the old ray grid, so they are invalid on rebuild regardless.
!
	  IF(ALLOCATED(NM_RAY))DEALLOCATE ( NM_RAY )
	  IF(ALLOCATED(DLN_RAY))DEALLOCATE ( DLN_RAY )
	  IF(ALLOCATED(MAX_LS))DEALLOCATE ( MAX_LS )
	  IF(ALLOCATED(AV))DEALLOCATE ( AV )
	  IF(ALLOCATED(CV))DEALLOCATE ( CV )
	  IF(ALLOCATED(DTAU))DEALLOCATE ( DTAU )
	  IF(ALLOCATED(GAM))DEALLOCATE ( GAM )
!
! These arrays are only used when using the DIFFERENCE apporach.
!
	  IF( ALLOCATED(AV_PREV) )THEN
	    DEALLOCATE ( AV_PREV )
	    DEALLOCATE ( AV_STORE )
	    DEALLOCATE ( PAR_AV )
	    DEALLOCATE ( CV_PREV )
	    DEALLOCATE ( CV_STORE )
	    DEALLOCATE ( GAMH )
!
	    DEALLOCATE ( TA )
	    DEALLOCATE ( TB )
	    DEALLOCATE ( TC )
	    DEALLOCATE ( GB )
	    DEALLOCATE ( H )
	    DEALLOCATE ( OLDCHI )
	    DEALLOCATE ( OLDCHI_STORE )
!
	    DEALLOCATE ( XM )
	    DEALLOCATE ( U )
	    DEALLOCATE ( VB )
	    DEALLOCATE ( VC )
	    DEALLOCATE ( QH )
	    DEALLOCATE ( DIV )
	  END IF
!
! These arrays are only used when using the INTEGERAL apporach.
!
	  IF( ALLOCATED(I_P) )THEN
	    DEALLOCATE ( I_P )
	    DEALLOCATE ( I_M )
!
	    DEALLOCATE ( I_P_PREV )
	    DEALLOCATE ( I_P_STORE )
	    DEALLOCATE ( I_M_PREV )
	    DEALLOCATE ( I_M_STORE )
!
	    DEALLOCATE ( A0 )
	    DEALLOCATE ( A1 )
	    DEALLOCATE ( A2 )
	    DEALLOCATE ( A3 )
	    DEALLOCATE ( A4 )
!
	    DEALLOCATE (EE)
	    DEALLOCATE (E0)
	    DEALLOCATE (E1)
	    DEALLOCATE (E2)
	    DEALLOCATE (E3)
!
	    DEALLOCATE (SOURCE_PRIME)
	    DEALLOCATE (S)
	    DEALLOCATE (dS)
	    IF(ALLOCATED(ALO_CHIHAT_RAY))DEALLOCATE ( ALO_CHIHAT_RAY )
	    IF(ALLOCATED(ALO_SPR_M_RAY))DEALLOCATE ( ALO_SPR_M_RAY )
	    IF(ALLOCATED(ALO_SPR_P_RAY))DEALLOCATE ( ALO_SPR_P_RAY )
!
! edited by Sean Brennan
	    IF(ALLOCATED(I_P_DEPTH_NEXT))DEALLOCATE ( I_P_DEPTH_NEXT )
	    IF(ALLOCATED(I_M_DEPTH_NEXT))DEALLOCATE ( I_M_DEPTH_NEXT )
	    IF(ALLOCATED(COMP_DEPTH_IDX))DEALLOCATE ( COMP_DEPTH_IDX )
!
! Per-ray Knop refresh store: invalid on any grid rebuild.
!
	    IF(ALLOCATED(KN_CHI_EXT))DEALLOCATE ( KN_CHI_EXT )
	    IF(ALLOCATED(KN_ETA_EXT))DEALLOCATE ( KN_ETA_EXT )
	    IF(ALLOCATED(KN_FREQ))DEALLOCATE ( KN_FREQ )
	    IF(ALLOCATED(KN_DLN_P))DEALLOCATE ( KN_DLN_P )
	    IF(ALLOCATED(KN_DLN_N))DEALLOCATE ( KN_DLN_N )
	    IF(ALLOCATED(KN_DBB))DEALLOCATE ( KN_DBB )
	    IF(ALLOCATED(KN_IC))DEALLOCATE ( KN_IC )
	    IF(ALLOCATED(KN_CHI_COEF))DEALLOCATE ( KN_CHI_COEF )
	    IF(ALLOCATED(KN_ETA_COEF))DEALLOCATE ( KN_ETA_COEF )
	    IF(ALLOCATED(KN_R))DEALLOCATE ( KN_R )
	    IF(ALLOCATED(KN_P))DEALLOCATE ( KN_P )
	    KN_STORE_READY=.FALSE.
! Re-arm the I_NEXT mapping diagnostic: the ray geometry and
! COMP_DEPTH_IDX are about to be rebuilt for the new grid.
	    INEXT_MAP_PRINTED=.FALSE.
		  END IF
	END IF
	VDOP_FRAC_SAV=VDOP_FRAC
!
! 
!
! Set up the revised grid to improve computational accuracy. Unles we are
! carrying out tests, these need only be constructed once.
!
	IF(FIRST_TIME .OR. .NOT. ALLOCATED(R_EXT) )THEN
!
	  OLD_SOLUTION_OPTIONS=SOLUTION_OPTIONS
	  IF(SOLUTION_OPTIONS .EQ. 'DIFF/INS')THEN
	    SOLUTION_METHOD='DIFFERENCE'
	    INSERT=.TRUE.
	  ELSE IF(SOLUTION_OPTIONS(1:4) .EQ. 'DIFF')THEN
	    SOLUTION_METHOD='DIFFERENCE'
	    INSERT=.FALSE.
	  ELSE IF(SOLUTION_OPTIONS .EQ. 'INT/INS')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.TRUE.
	  ELSE IF(SOLUTION_OPTIONS(1:3) .EQ. 'INT')THEN
	    SOLUTION_METHOD='INTEGRAL'
	    INSERT=.FALSE.
	  ELSE
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_NONMONO_V1: Invalid solution option'
	    WRITE(J,*)SOLUTION_OPTIONS
	    STOP
	  END IF
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_NONMONO_V1:'
	    WRITE(J,*)'The DIFFERENCE method lacks sign-dependent upwinding.'
	    WRITE(J,*)'Use INT or INT/INS for non-monotonic transfer.'
	    STOP
	  END IF
!
	  IF(INNER_BND_METH .NE. 'DIFFUSION' .AND.
	1     INNER_BND_METH .NE. 'ZERO_FLUX' .AND.
	1     INNER_BND_METH .NE. 'PNT_SRCE' .AND.
	1     INNER_BND_METH .NE. 'SCHUSTER')THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_NONMONO_V1: Invalid inner boundary condition'
	    WRITE(J,*)INNER_BND_METH
	    STOP
	  END IF
!
          ND_ADD=0
          IF(THK)ND_ADD=ND_ADD_MAX
          ND_EXT=ND+ND_ADD
!
	  ALLOCATE ( R_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_R_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( V_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( Z_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_ETA_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_CHI_EXT(ND_EXT),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating R_EXT block in FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP	
          END IF
!
! Compute the extended R grid, excluding inserted points.
!
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THK)THEN
	    IF(REXT_FAC .GT. 1.0_LDP .AND. REXT_FAC .LT. 10.0_LDP)THEN
	       RMAX=REXT_FAC*R(1)
	    ELSE IF(V(ND) .LT. 10.0_LDP .AND. R(1)/R(ND) .GE. 9.99_LDP)THEN
	      RMAX=10.0_LDP*R(1)		!Stellar wind
	    ELSE IF(V(ND) .GT. 10.0_LDP .OR. V(1) .GT. 2.0E+04_LDP)THEN
	      RMAX=1.5_LDP*R(1)		!SN model
	    ELSE
	      RMAX=MIN(10.0_LDP,SQRT(R(1)/R(ND)))*R(1)
	    END IF
	    ALPHA=R(1)+(R(1)-R(2))
	    DEL_R_FAC=EXP( LOG(RMAX/ALPHA)/(ND_ADD-3) )
	    R_EXT(1)=RMAX
	    R_EXT(4)=RMAX/DEL_R_FAC
	    R_EXT(2)=R_EXT(1)-0.1_LDP*(R_EXT(1)-R_EXT(4))
	    R_EXT(3)=R_EXT(1)-0.4_LDP*(R_EXT(1)-R_EXT(4))
	    DO I=5,ND_ADD-1
	      R_EXT(I)=R_EXT(I-1)/DEL_R_FAC
	    END DO
	    R_EXT(ND_ADD)=ALPHA
	  END IF
!
! Compute VEXT and R_EXT. We assume a BETA velocity law at large R.
!
	  V_EXT(ND_ADD+1:ND_EXT)=V(1:ND)
!
! Use the TRUE (unclamped) SIGMA for the main-grid portion of the extended
! grid.  This is the single point where the true negative SIGMA values
! (shocks, velocity reversals) enter the formal solver's physics.  All
! downstream quantities (SIGMA_COEF, SIGMA_RAY, GAM, T2, resonance scan)
! are derived from SIGMA_EXT and therefore automatically use the
! true values.  The SIGMA dummy is used only for the outer beta-law
! extension (THK block).
!
! edited by Sean Brennan
	  SIGMA_EXT(ND_ADD+1:ND_EXT)=SIGMA_TRUE_IN(1:ND)
	  IF(THK)THEN
!
! 05-Oct-2026 (Devin): SIGMA is NOT re-clamped under
! DO_NONMONO_KNOP, so a shock/inversion reaching the outer edge
! (sigma(1) < -1) would give BETA < 0 and V_EXT would diverge as
! R_EXT -> RMAX.  The beta law only makes sense for a coasting
! outer flow; floor the exponent source at -0.999, which limits the
! extension to nearly constant V -- the least-wrong extrapolation
! when the outer edge is in compression.  Warn so the user knows the
! thick boundary rests on an approximation.
!
	    T1=MAX(SIGMA(1),-0.999_LDP)
	    IF(SIGMA(1) .LT. -0.999_LDP .OR. V(1) .LE. 0.0_LDP)THEN
	      J=ERROR_LU()
	      WRITE(J,'(A)')'Warning: non-monotonic outer edge in'
	1      //' FG_J_CMF_NONMONO_V1 THK extension'
	      WRITE(J,'(A,2ES16.8)')'  SIGMA(1),V(1)=',SIGMA(1),V(1)
	      WRITE(J,'(A)')'  beta-law extension floored at sigma=-0.999'
	1      //' (near-constant V beyond R(1))'
	    END IF
	    BETA=(T1+1.0_LDP)*(R(1)/R(ND)-1.0_LDP)
            VINF=V(1)/(1-R(ND)/R(1))**BETA
	    DO I=1,ND_ADD
	      V_EXT(I)=VINF*(1.0_LDP-R_EXT(ND_EXT)/R_EXT(I))**BETA
	      SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0_LDP)-1.0_LDP
	    END DO
	    J=ERROR_LU()
	    WRITE(J,'(A)')' '
	    WRITE(J,*)'Using thick boundary condition in FG_J_CMF_NONMONO_V1'
	    WRITE(J,'(2(A,ES16.8,3X))')' R(1)=',R(1),'RMAX=',RMAX
	    WRITE(J,'(2(A,ES16.8,3X))')' V(1)=',V(1),'VMAX=',V_EXT(1)
	  END IF
!
! Define zone used to extrapolate opacities and emissivities.
!
	  IF(ND_ADD .NE. 0)THEN
	    IDMIN=1
	    T1=R_EXT(1)/R(1)
	    IF(T1 .LT. R(1)/R(ND/3))THEN
	      T1=MIN(3.0_LDP,T1)
	      DO I=1,ND
	        IF(R(1)/R(I) .GT. T1)THEN
	          IDMAX=MAX(4,I)
	          EXIT
	        END IF
	      END DO
	      IDMAX=MIN(IDMAX,ND/6)
	    ELSE
!	      IDMAX=ND/6
	      IDMAX=MIN(10,ND/6)
	    END IF
	    WRITE(6,'(A,I4,A)')' Using depth 1 and',IDMAX,' to extrapolate opacities'
	  END IF
!
! Work out the maximim number of points per ray. This will allow us to allocate
! the require memory.
!
	  NRAY_MAX=0
	  DO LS=1,NP
	    NI_SMALL=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND_EXT
!
	    DO I=1,NI_SMALL
	      IF(R_EXT(I) .EQ. P(LS))THEN
	        Z_EXT(I)=0.0_LDP
	      ELSE
	        Z_EXT(I)=SQRT(R_EXT(I)*R_EXT(I)-P(LS)*P(LS))
	      END IF
	   END DO
!
	    K=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD))
	    IF(.NOT.(T2 .GT. 0.0_LDP))THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'FG_J_CMF_NONMONO_V1: invalid Doppler grid'
	1     ,' scale T2=',T2,' (VDOP_FRAC=',VDOP_FRAC,')'
	      STOP
	    END IF
	    DO I=1,NI_SMALL-1
	      T1=ABS(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01_LDP,T1)
	      IF(T1 .GT. KN_RAY_MAX_INS)T1=REAL(KN_RAY_MAX_INS,KIND=LDP)
	      IF(T1 .GT. 1.0_LDP)K=K+INT(T1)
	      K=K+1
	    END DO
            NRAY_MAX=MAX(NRAY_MAX,K)
	  END DO
	  J=ERROR_LU()
	  WRITE(J,'(A,I6)')' Maximum number of points along ray in FG_J_CMF_NONMONO_V1 is',NRAY_MAX
!
	  ALLOCATE ( Z(NRAY_MAX,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( RAY_PNT(NRAY_MAX,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( R_RAY(NRAY_MAX,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( NI_RAY(NP),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating Z block in FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP	
          END IF
!
! 
!
! Define the grid along each ray. NB: RAY_PNT is used to indicate the
! interpolation interval on the original (EXT) grid.
!
	  DO LS=1,NP
	    NI_SMALL=ND_EXT-(LS-NC-1)
	    IF(LS .LE. NC+1)NI_SMALL=ND_EXT
!
	    DO I=1,NI_SMALL
	      IF(R_EXT(I) .EQ. P(LS))THEN
	        Z_EXT(I)=0.0_LDP
	      ELSE
	        Z_EXT(I)=SQRT(R_EXT(I)*R_EXT(I)-P(LS)*P(LS))
	      END IF
	    END DO
!
	    K=1
	    Z(1,LS)=Z_EXT(1)
	    RAY_PNT(1,LS)=1
	    T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:NI_SMALL-ND_ADD))
	    DO I=1,NI_SMALL-1
	      T1=ABS(Z_EXT(I)*V_EXT(I)/R_EXT(I)-Z_EXT(I+1)*V_EXT(I+1)/R_EXT(I+1))/T2
	      IF(.NOT. THK .AND. NI_SMALL .EQ. 2)T1=MAX(2.01_LDP,T1)
	      IF(T1 .GT. KN_RAY_MAX_INS)T1=REAL(KN_RAY_MAX_INS,KIND=LDP)
	      IF(T1 .GT. 1.0_LDP)THEN
	        DELTA_Z=(Z_EXT(I+1)-Z_EXT(I))/(INT(T1)+1)
	        DO J=1,INT(T1)
	          K=K+1
	          Z(K,LS)=Z(K-1,LS)+DELTA_Z
	          RAY_PNT(K,LS)=I
	        END DO
	      END IF
	      K=K+1
	      Z(K,LS)=Z_EXT(I+1)
	      RAY_PNT(K,LS)=I
	    END DO
	    NI_RAY(LS)=K
!
	    PSQ=P(LS)*P(LS)
	    R_RAY(1,LS)=R_EXT(1)
	    DO I=2,NI_RAY(LS)
	      R_RAY(I,LS)=SQRT(Z(I,LS)*Z(I,LS)+PSQ)
	    END DO
	    R_RAY(NI_RAY(LS),LS)=R_EXT(NI_SMALL)
!
	  END DO
!
! 
!
! J_PNT and H_PNT are used to indicate the postion of J(I) amd H(I) along the ray
! so that J and H can be computed on our regular radius grid.
!
	  ALLOCATE ( J_PNT(ND,NP),STAT=IOS); J_PNT(:,:)=0
	  IF(IOS .EQ. 0)ALLOCATE ( H_PNT(ND,NP),STAT=IOS); H_PNT(:,:)=0
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating J_PNT block in FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP	
          END IF
!
! LS=NP with only 1 depth point must be treated separately.
!
	  NP_TMP=NP; IF(NI_RAY(NP) .EQ. 1)NP_TMP=NP-1
	  DO LS=1,NP_TMP
	    NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	    K=1
	    DO I=1,NI_SMALL
!
! edited by Devin, 29-Sep-2026: bound the march in-loop.  The ray
! grid here is hand-built with inserted variable-density points, so
! a depth that brackets no R_RAY interval (endpoint roundoff) would
! otherwise walk K past NI_RAY(LS) and read out of bounds -- an
! undiagnosed hang rather than a crash in an unchecked build.
!
	      DO WHILE(J_PNT(I,LS) .EQ. 0)
	        IF(R(I) .LE. R_RAY(K,LS) .AND. R(I) .GE. R_RAY(K+1,LS))THEN
	          IF( (R_RAY(K,LS)-R(I)) .LT. (R(I)-R_RAY(K+1,LS)) )THEN
	            J_PNT(I,LS)=K
	          ELSE
	            J_PNT(I,LS)=K+1
	          END IF
	        ELSE
	          K=K+1
	          IF(K+1 .GT. NI_RAY(LS))THEN
	            WRITE(LUER,*)'Error setting J_PNT in'//
	1            ' FG_J_CMF_NONMONO_V1 -- no ray bracket for depth'
	            WRITE(LUER,*)'Depth=',I,'Ray=',LS,'R(I)=',R(I)
	            WRITE(LUER,*)'R_RAY(1,LS)=',R_RAY(1,LS),
	1            'R_RAY(NI,LS)=',R_RAY(NI_RAY(LS),LS)
	            STOP
	          END IF
	        END IF
	      END DO
	    END DO
	  END DO
	  IF(NI_RAY(NP) .EQ. 1)J_PNT(1,NP)=1
!
	  DO I=1,ND
	    DO LS=1,NP-I+1
	      K=J_PNT(I,LS)
	      IF(J_PNT(I,LS) .LT. 1 .OR. J_PNT(I,LS) .GT. NI_RAY(LS))THEN
	        WRITE(LUER,*)'Error setting J_PNT in FG_J_CMF_NONMONO_V1 -- invalid values'
	        WRITE(LUER,*)'Depth=',I,'Ray=',LS,'J_PNT value=',J_PNT(I,LS)
	        WRITE(LUER,*)'NP=',NP,'R(1)=',R(1),'R_RAY(1,LS)=',R_RAY(1,LS)
	        WRITE(LUER,*)'NI_SMALL=',NI_SMALL
	        STOP
	      ELSE IF( ABS(R_RAY(K,LS)-R(I))/R(I) .GT. 1.0E-12_LDP)THEN
	        WRITE(LUER,*)'Error setting J_PNT in FG_JCMF_V12 -- invalid values'
	        WRITE(LUER,*)'Fractional difference is',ABS(R_RAY(K,LS)-R(I))/R(I)
	        WRITE(LUER,*)'Depth=',I,'Ray=',LS,'J_PNT value=',J_PNT(I,LS)
	        STOP
	      END IF
	    END DO
	  END DO
!
	  IF( SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
! CV is defined on the nodes. Need to interpolate to the midpoint of R.
!
	    DO LS=1,NP
	      NI_SMALL=ND-(LS-NC-1);   IF(LS .LE. NC+1)NI_SMALL=ND
	      K=1
	      DO I=1,NI_SMALL-1
	        DO WHILE(H_PNT(I,LS) .EQ. 0)
	          T1=0.5_LDP*(R(I)+R(I+1))
	          IF(T1 .LE. R_RAY(K,LS) .AND. T1 .GE. R_RAY(K+1,LS))THEN
	            H_PNT(I,LS)=K
	          ELSE
	            K=K+1
	            IF(K+1 .GT. NI_RAY(LS))THEN
	              WRITE(LUER,*)'Error setting H_PNT in'//
	1              ' FG_J_CMF_NONMONO_V1 -- no ray bracket for midpoint'
	              WRITE(LUER,*)'Depth=',I,'Ray=',LS,'T1=',T1
	              WRITE(LUER,*)'R_RAY(1,LS)=',R_RAY(1,LS),
	1              'R_RAY(NI,LS)=',R_RAY(NI_RAY(LS),LS)
	              STOP
	            END IF
	          END IF
	        END DO
	      END DO
	    END DO
	  ELSE
!
! CV is defined at the mid-points, but these don't, in general, correspond to
! the midpoints of the R grid.
!
	    DO LS=1,NP
	      K=1
	      DO I=1,ND-1
	        DO WHILE(H_PNT(I,LS) .EQ. 0)
	          T1=0.5_LDP*(R(I)+R(I+1))
	          IF(T1 .LE. 0.5_LDP*(R_RAY(K,LS)+R_RAY(K+1,LS)) .AND. T1 .GE. 0.5_LDP*(R_RAY(K+1,LS)+R_RAY(K+2,LS)) )THEN
	            H_PNT(I,LS)=K
	          ELSE
	            K=K+1
	            IF(K+2 .GT. NI_RAY(LS))THEN
	              WRITE(LUER,*)'Error setting H_PNT in'//
	1              ' FG_J_CMF_NONMONO_V1 -- no ray bracket for midpoint'
	              WRITE(LUER,*)'Depth=',I,'Ray=',LS,'T1=',T1
	              WRITE(LUER,*)'R_RAY(1,LS)=',R_RAY(1,LS),
	1              'R_RAY(NI,LS)=',R_RAY(NI_RAY(LS),LS)
	              STOP
	            END IF
	          END IF
	        END DO
	      END DO
	    END DO
	  END IF
!
! 
!
	  ALLOCATE ( V_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF(ND_EXT,4),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF(ND_EXT,4),STAT=IOS )
!
!***************************************************************************
!***************************************************************************
!
	  IF(IOS .EQ. 0)ALLOCATE ( V_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( dCHIdR_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SOURCE_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( GAM_MAG(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( dGAMdR(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( dCHIdR(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( Q(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( NM_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( DLN_RAY(NRAY_MAX),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( MAX_LS(NRAY_MAX),STAT=IOS )
!
	  IF(IOS .EQ. 0)ALLOCATE ( AV(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CV(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( DTAU(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( GAM(NRAY_MAX,NP),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating V_RAY etc in FG_J_CMF_NONMONO_V1: Status=',IOS
	    STOP	
          END IF
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by both DIFFERENCE solution approach.
!
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    ALLOCATE ( AV_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( AV_STORE(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( PAR_AV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( CV_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( CV_STORE(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( GAMH(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( TA(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( TB(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( TC(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( GB(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( H(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( OLDCHI(NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( OLDCHI_STORE(NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( XM(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( U(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( VB(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( VC(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( QH(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( DIV(NRAY_MAX),STAT=IOS )
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error allocating DIFFERENCE block in FG_J_CMF_NONMONO_V1: Status=',IOS
	      STOP	
            END IF
	  END IF
!
!***************************************************************************
!***************************************************************************
!
! Allocation of variables required by INTEGRAL solution.
!
	  IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
!   Check whether the input model has any compressing points (SIGMA < -1).
!   If not, skip the large (ND,NP,NCF) Gauss-Seidel arrays entirely.
!
! edited by Sean Brennan
! GAM=(V/R)(1+sigma*mu^2): V<0 makes GAM<0 at every mu regardless of
! sigma, so backflow cells need the same GS store as sigma<-1 cells.
!
	    HAS_COMPRESSING=ANY(SIGMA_TRUE_IN(1:ND) .LT. -1.0_LDP) .OR.
	1                    ANY(V(1:ND) .LT. 0.0_LDP)
	    N_COMP=0
	    IF(HAS_COMPRESSING)THEN
	      ALLOCATE(COMP_DEPTH_IDX(ND),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error allocating COMP_DEPTH_IDX: Status=',IOS
	        STOP
	      END IF
	      COMP_DEPTH_IDX=0
	      DO I=1,ND
	        IF(SIGMA_TRUE_IN(I) .LT. -1.0_LDP .OR.
	1          V(I) .LT. 0.0_LDP)THEN
	          COMP_DEPTH_IDX(MAX(1,I-1):MIN(ND,I+1))=1
	        END IF
	      END DO
	      DO I=1,ND
	        IF(COMP_DEPTH_IDX(I) .GT. 0)THEN
	          N_COMP=N_COMP+1
	          COMP_DEPTH_IDX(I)=N_COMP
	        END IF
	      END DO
	      WRITE(LUER,'(A)')'FG_J_CMF_NONMONO_V1: Gauss-Seidel storage'
	      WRITE(LUER,'(A,I6,A,I6)')
	1       '    compressing points=',N_COMP,'  of ND=',ND
	      WRITE(LUER,'(A)')
	1       '    sparse in-place write (Opt A+B, no CURR arrays)'
	    ELSE
	      WRITE(LUER,'(A)')
	1       'FG_J_CMF_NONMONO_V1: no compressing points — GS storage skipped'
	    END IF
!
	    ALLOCATE (I_P(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (I_M(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( I_P_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( I_P_STORE(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( I_M_PREV(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( I_M_STORE(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE ( A0(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A1(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A2(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A3(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE ( A4(NRAY_MAX,NP),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE (EE(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E0(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E1(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E2(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (E3(NRAY_MAX),STAT=IOS )
!
	    IF(IOS .EQ. 0)ALLOCATE (SOURCE_PRIME(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (S(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (dS(NRAY_MAX),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (ALO_CHIHAT_RAY(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (ALO_SPR_M_RAY(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (ALO_SPR_P_RAY(NRAY_MAX,NP),STAT=IOS )
	    IF(IOS .EQ. 0)ALLOCATE (FG_NOSLOT_RAY(NP),STAT=IOS )
	    IF(IOS .EQ. 0)FG_NOSLOT_RAY=0
!
! Optimized Knop et al. 2009 depth-point intensity storage.
! Sparse depth storage (Opt A): only N_COMP compressing depths.
! In-place write (Opt B): write directly to NEXT, no CURR arrays.
! Memory: 2 * N_COMP * NP * NCF_IN * 8 bytes.
!
! edited by Sean Brennan
	    IF(IOS .EQ. 0 .AND. HAS_COMPRESSING)THEN
	      ALLOCATE (I_P_DEPTH_NEXT(N_COMP,NP,NCF_IN),STAT=IOS )
	      IF(IOS .EQ. 0)I_P_DEPTH_NEXT=0.0_LDP
	      IF(IOS .EQ. 0)ALLOCATE (I_M_DEPTH_NEXT(N_COMP,NP,NCF_IN),STAT=IOS )
	      IF(IOS .EQ. 0)I_M_DEPTH_NEXT=0.0_LDP
	    END IF
!
! Per-ray Knop refresh store (Oct-2026): captured per-frequency
! sweep inputs used by FG_KNOP_REFRESH_V1 to converge the lagged
! compressing-branch coupling after each frequency sweep.  Same
! lifetime as the compressing store it feeds.
!
	    IF(IOS .EQ. 0 .AND. HAS_COMPRESSING)THEN
	      ALLOCATE (KN_CHI_EXT(ND_EXT,NCF_IN),
	1       KN_ETA_EXT(ND_EXT,NCF_IN),
	2       KN_FREQ(NCF_IN),
	3       KN_DLN_P(NCF_IN),KN_DLN_N(NCF_IN),
	4       KN_DBB(NCF_IN),KN_IC(NCF_IN),
	5       KN_CHI_COEF(ND_EXT,4,NCF_IN),
	6       KN_ETA_COEF(ND_EXT,4,NCF_IN),
	7       KN_R(ND),KN_P(NP),STAT=IOS)
	      IF(IOS .EQ. 0)THEN
	        KN_CHI_EXT=0.0; KN_ETA_EXT=0.0
	        KN_FREQ=0.0_LDP
	        KN_DLN_P=1.0_LDP; KN_DLN_N=1.0_LDP
	        KN_DBB=0.0_LDP; KN_IC=0.0_LDP
	        KN_CHI_COEF=0.0_LDP; KN_ETA_COEF=0.0_LDP
	        KN_R(1:ND)=R(1:ND); KN_P(1:NP)=P(1:NP)
	        KN_ND=ND; KN_NP=NP; KN_NC=NC; KN_NCF=NCF_IN
	        KN_ND_EXT=ND_EXT
	        KN_INNER_METH=INNER_BND_METH
	        KN_METHOD_SAV=METHOD
	      END IF
	    END IF
		    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error allocating INTEGRAL block in FG_J_CMF_NONMONO_V1: Status=',IOS
	      STOP	
            END IF
!
! GS restart (16-Sep-2026): save dims for the checkpoint dump
! (SCR_RITE_V2 cannot see NP/NCF_IN), then attempt to restore a
! validated complete I_NEXT generation from FGINEXT_GS.  On any
! validation failure the arrays stay zeroed (cold start) and
! I_NEXT_GS_OK stays .FALSE., so the N_LAM_RWARM LAMBDA warm-up in
! SOLVE_FOR_POPS still applies as the fallback.
!
	    FG_GS_NP=NP
	    FG_GS_NCF=NCF_IN
	    IF(GS_OBS_MODE)THEN
	      CALL RESTORE_FG_INEXT_OBS(R,V,SIGMA,ND,NP,GS_NU_OBS,NCF_IN)
	    ELSE
	      CALL RESTORE_FG_INEXT_GS(R,V,SIGMA,ND,NP,NCF_IN)
	    END IF
!
	  END IF
!
	  FIRST_TIME=.FALSE.
	ELSE
	  IF(OLD_SOLUTION_OPTIONS .NE. SOLUTION_OPTIONS)THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in FG_J_CMF_NONMONO_V1'
	    WRITE(J,*)'Can''t switch SOLUTION_OPTIONS while runing code'
	    WRITE(J,*)'New setting:',SOLUTION_OPTIONS
	    WRITE(J,*)'Old setting:',OLD_SOLUTION_OPTIONS
	    STOP
	  END IF
	END IF
!
!
!
	NEG_AV_VALUE=.FALSE.
!
! Perform initializations.
!
	IF(INIT)THEN
!
! Insert extra points into radius grid. Not all points will be used along a ray.
! We only insert additional points in the interval between Z(NI)=0 and Z(NI-1),
! and between Z(NI-1) and  Z(NI-2).
!
	  FG_ERR_ON_FREQ(:)=0.0_LDP
	  FG_ERR_TYPE(:)=0
	  FG_ERR_CNT=0
	  IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    AV_PREV(:,:)=0.0_LDP
	    CV_PREV(:,:)=0.0_LDP
          ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=0.0_LDP
	    I_M_PREV(:,:)=0.0_LDP
!
! Optimized Knop et al. 2009: in-place write to NEXT (Optimization B).
! No CURR arrays and no INIT copy/zero needed.  The NEXT arrays already
! contain the intensities from the previous great iteration (written
! in-place during the sweep).  On the first iteration, NEXT is zero
! (allocated + zeroed), so compressing points get SOURCE_PRIME = SOURCE
! (cold start), identical to the original behaviour.
!
! edited by Sean Brennan
	    FREQ_IDX_ITER=0
	    FG_GS_MAXCHG=0.0_LDP
	    FG_GS_NCONV=0
	    FG_INEXT_CAPPED=0
	    FG_SW_CAPPED=0
	    FG_SW_CAPTOT=0
	    FG_RF_CAPPED=0
	    FG_GS_NBAD=-1
	    FG_GS_NTOT=0
!
! W1-W9 resets (17-Oct-2026): every counter the maturity gate reads
! must clear at sweep entry so a single bad iteration cannot veto
! the release path forever.
!
	    FG_SW_MAXCHG=0.0_LDP
	    FG_SW_NCONV=0
	    FG_GS_NDIV=0
	    FG_GS_NBAD_BLK=-1
	    FG_GS_NTOT_BLK=0
	    FG_GS_NSTAL=0
	    FG_GS_NPLAT=0
	    FG_DIV_HEAL=0
	    FG_COMMIT_REJ=0
	    FG_NC_REJ=0
	    FG_NC_CAP=0
	    KN_HOT_SPENT=0
	    KN_STORE_READY=.FALSE.
	    FG_NOSLOT_CNT=0
	  END IF
!
	  LOG_R_EXT(1:ND_EXT)=LOG(R_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(V_COEF,V_EXT,LOG_R_EXT,ND_EXT)
	  CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_EXT,LOG_R_EXT,ND_EXT)
!
! Non-monotonic velocity laws (shocks, cold dense shells, reverse shocks)
! routinely give SIGMA = dlnV/dlnR - 1 < -1 in compression zones.  This is
! EXPECTED when DO_NONMONO_KNOP is active and is not an error.  We count the
! offending depths for diagnostics but do not stop or warn per-point.
!
	  N_NEG_SIGMA=0
	  DO I=1,ND
	    IF(SIGMA_TRUE_IN(I) .LT. -1.0_LDP)N_NEG_SIGMA=N_NEG_SIGMA+1
	  END DO
	  IF(N_NEG_SIGMA .GT. 0 .AND. INIT .AND. .NOT. DIAG_PRINTED)THEN
	    WRITE(LUER,'(A,I6,A)')'FG_J_CMF_NONMONO_V1: ',N_NEG_SIGMA,
	1       ' depths have SIGMA < -1 (non-monotonic flow handled).'
	  END IF
	
	  DO LS=1,NP
!
! The check on SIGMA is to allow for the possibility of a non-monotonic velcoity law in the
! region where V is small, and hece where SIGAMA is unimportant.
!
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      T1=LOG(R_RAY(I,LS)/R_EXT(K))
	      V_RAY(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
	      SIGMA_RAY(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
!
!        No SIGMA_RAY clamp: the Knop et al. 2009 formalism handles
!        arbitrary negative SIGMA by keeping chi_hat = chi + |GAM|/dLOG_NU
!        (always positive opacity) and moving the velocity-coupling
!        term into the modified source S' = SOURCE + Q*I_upwind.
!        The signed GAM (1+SIGMA*MU^2) enters the source term via Q,
!        not the opacity.  GAM_MAG = |GAM| is used to construct chi_hat
!        (line: Q=GAM_MAG/dLOG_NU; CHI_RAY=CHI_RAY+Q) and for resonance
!        counting.  SIGMA_RAY is derived from SIGMA_EXT, which is
!        populated from the unclamped SIGMA_TRUE_IN, so non-monotonic
!        velocity fields (sigma < -1) are correctly handled.
! edited by Sean Brennan
	    END DO
!
!
! Compute GAMMA for a NON-MONOTONIC velocity field.
!
! T2 = 1 + SIGMA*MU^2  is the comoving-frame frequency gradient along the ray
!     (= d(mu*V)/dz up to the V/R prefactor).  In a monotonic wind T2 > 0.
!     In compression zones (SIGMA < -1, |MU| large) T2 < 0: the flow is
!     converging, the Sobolev resonance is approached from the opposite
!     direction, and the naive Sobolev optical depth CHIL*R/V/T2 goes
!     negative.
!
! Non-monotonic handling (replaces the Dessart+2015 softened-|T2| guard and
! the SIGMA >= -0.999 clamps used in FG_J_CMF_V13):
!
!   * GAM(I,LS)  = T1 * T2          (TRUE SIGNED value; retained so the
!                                    source-function integration can account
!                                    for ray direction relative to the flow)
!   * GAM_MAG(I,LS) = |T1| * |T2|      (magnitude; always >= 0; used to
!                                    construct chi_hat and for DTAU guards.
!                                    ABS(T1) is required so V < 0 is safe.)
!   * SIGN(1,GAM(I,LS)) gives the local expanding/compressing direction.
!
! No |T2| floor is applied here.  At a true frequency turning point (T2=0),
! GAM=0 and GAM_MAG=0, so the velocity contribution to chi_hat correctly
! vanishes.  Division by zero is prevented downstream by the chi_hat floor
! (CHI_RAY >= 1.0E-30) and the DTAU zero-guard, not by flooring T2.
! edited by Sean Brennan
!
	NI=NI_RAY(LS)
	DO I=1,NI_RAY(LS)
	  MU=Z(I,LS)/R_RAY(I,LS)
	  T1=3.33564E-06_LDP*V_RAY(I)/R_RAY(I,LS)
	  T2=1.0_LDP+SIGMA_RAY(I)*(MU**2)
	  GAM(I,LS)=T1*T2
	  GAM_MAG(I,LS)=ABS(T1)*ABS(T2)
	END DO
!
! Compute dGAMdR for the INTEGRAL solution method.
!
! This matches FG_J_CMF_V13 lines 975-989, but uses the TRUE (unclamped)
! SIGMA_RAY so the derivative reflects the real velocity gradient.
!
! dGAMdR enters the Euler-Maclaurin correction in DTAU via:
!   dCHIdR_RAY += SIGN(1,GAM) * dGAMdR / dLOG_NU
!
! SIGN(1,GAM) gives the correct sign for d|GAM|/dR:
!   GAM >= 0 (expanding):  d|GAM|/dR = +dGAM/dR
!   GAM < 0 (compressing): d|GAM|/dR = -dGAM/dR
!
! This makes non-monotonic regions consistent with V13 in monotonic
! regions (where SIGN(1,GAM) = +1 and d|GAM|/dR = dGAM/dR, identical
! to V13).  The DTAU clamping handles any cusps at velocity turning
! points where d|GAM|/dR is discontinuous.
!
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL' .AND. NI_RAY(LS) .NE. 1)THEN
	  DO I=1,NI_RAY(LS)
	    MU=Z(I,LS)/R_RAY(I,LS)
	    T1=3.33564E-06_LDP*V_RAY(I)/R_RAY(I,LS)
	    dGAMdR(I,LS)=GAM(I,LS)*SIGMA_RAY(I)/R_RAY(I,LS)
	    J=MAX(I-1,1); K=MIN(NI,I+1)
!
! A degenerate ray (duplicate radii after insertion) can make the
! spacing identically zero; skipping the central difference leaves
! the sigma-scaled and curvature terms, which are still defined.
!
	    T4=R_RAY(J,LS)-R_RAY(K,LS)
	    IF(ABS(T4) .GT. 1.0E-30_LDP*R_RAY(I,LS))THEN
	      dGAMdR(I,LS)=dGAMdR(I,LS)+T1*(MU**2)*
	1       (SIGMA_RAY(J)-SIGMA_RAY(K))/T4
	    END IF
!
! Ray-curvature term: MU=Z/R varies along the ray even for fixed impact
! parameter.  dMU/dR = (1-MU^2)/(MU*R), so
!   T1 * d(SIGMA*MU^2)/dR gets an extra T1 * SIGMA * 2*MU * dMU/dR
!                             = T1 * 2*SIGMA*(1-MU^2)/R.
! This term was missing (audit 2026-09-07, finding #15).  It is required
! so that d(chi_hat)/dR used in the Euler-Maclaurin DTAU correction is the
! true derivative of |GAM| along the ray.
! edited by Devin (10-Sep-2026)
!
	    dGAMdR(I,LS)=dGAMdR(I,LS)+
	1       T1*2.0_LDP*SIGMA_RAY(I)*(1.0_LDP-MU**2)/R_RAY(I,LS)
	  END DO
	END IF
!
! Multi-root resonance scan (diagnostic + directionality bookkeeping).
!
! For a non-monotonic v(r) the comoving frequency nu_cm(r) = nu_obs*(1-V/c)
! can cross a given rest frequency multiple times along a ray.  We bracket
! every zone [I,I+1] where T2 changes sign (i.e. where d(nu_cm)/ds vanishes)
! -- these are the resonance turning points.  No interpolation is performed
! here; the piecewise short-characteristics integration below naturally
! sums the source-function contribution over each monotonic sub-interval
! because DTAU is built from chi_hat = chi + |GAM|/dLOG_NU (always > 0)
! and the signed GAM enters only the source-function advection term
! Q*I_upwind (I_PREV for GAM>=0, I_NEXT for GAM<0). The count is
! recorded for diagnostics so the
! user can verify that multiple resonances were present and handled.
! edited by Sean Brennan
!
	IF(NI_RAY(LS) .GT. 1 .AND. INIT .AND. .NOT. DIAG_PRINTED)THEN
	  NRES=0
	  DO I=1,NI_RAY(LS)-1
	    TA2=1.0_LDP+SIGMA_RAY(I)*(Z(I,LS)/R_RAY(I,LS))**2
	    TB2=1.0_LDP+SIGMA_RAY(I+1)*(Z(I+1,LS)/R_RAY(I+1,LS))**2
	    IF(TA2*TB2 .LE. 0.0_LDP)NRES=NRES+1
	  END DO
	  IF(NRES .GT. 1)N_RAYS_MULTI=N_RAYS_MULTI+1
	  IF(NRES .GT. NRES_MAX)NRES_MAX=NRES
	  IF(NRES_MIN .EQ. 0 .OR. NRES .LT. NRES_MIN)NRES_MIN=NRES
	END IF
!
! Half-point GAMH for the DIFFERENCE solution method.
!
! We use the TRUE SIGMA_RAY (no -0.999 clamp) so the derivative reflects the
! real compression gradient.  GAMH is computed directly from SIGMA_RAY at
! the half-point and then clamped to non-negative (IF GAMH<0: GAMH=-GAMH)
! so the resulting DTAU increments stay non-negative.  GAM_MAG is not used
! here; GAMH has its own absolute-value clamp.
!
	IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	  DO I=1,NI_RAY(LS)-1
	    MU=(Z(I,LS)+Z(I+1,LS))/(R_RAY(I,LS)+R_RAY(I+1,LS))
	    T1=SIGMA_RAY(I)
	    T2=SIGMA_RAY(I+1)
	    GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564E-06_LDP*
	1                   (  1.0_LDP+0.5_LDP*(MU**2)*
	1                   (T1+T2)  )/
	1                   (R_RAY(I,LS)+R_RAY(I+1,LS))
!
!	    Enforce non-negative |GAMH| magnitude for the optical-depth
!	    increment while preserving the sign for directionality.
!
	    IF(GAMH(I,LS) .LT. 0.0_LDP)GAMH(I,LS)=-GAMH(I,LS)
	  END DO
	  GAMH(NI,LS)=0.0_LDP
	END IF
!
	END DO		!LS Loop
!
! Print multi-root resonance summary (once per run) and mark diagnostics done.
!
! edited by Sean Brennan
	IF(INIT .AND. .NOT. DIAG_PRINTED)THEN
!
!	    Knop et al. 2009 diagnostics: sigma and chi_hat statistics.
!	    chi_hat = chi + |GAM|/dLOG_NU is always positive by construction.
!	    sigma values are the unclamped SIGMA_TRUE_IN from the hydro model.
!
! edited by Sean Brennan
	  T1=MINVAL(SIGMA_TRUE_IN(1:ND))
	  T2=MAXVAL(SIGMA_TRUE_IN(1:ND))
	  NSIG_NEG=0
	  NSIG_NEG1=0
	  DO I=1,ND
	    IF(SIGMA_TRUE_IN(I) .LT. 0.0_LDP)NSIG_NEG=NSIG_NEG+1
	    IF(SIGMA_TRUE_IN(I) .LT. -1.0_LDP)NSIG_NEG1=NSIG_NEG1+1
	  END DO
	  WRITE(LUER,'(A)')'FG_J_CMF_NONMONO_V1: non-monotonic formal solver'
	  IF(N_RAYS_MULTI .GT. 0)THEN
	    WRITE(LUER,'(A,I4,A,I3,A,I3,A)')
	1       '    multi-root rays=',N_RAYS_MULTI,
	1       '  (min=',NRES_MIN,'  max=',NRES_MAX,' turning points)'
	  ELSE
	    WRITE(LUER,'(A)')'    no multi-root resonances detected'
	  END IF
	  WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1     '    sigma_min=',T1,'  sigma_max=',T2
	  WRITE(LUER,'(A,I6,A,I6)')
	1     '    N(sigma<0)=',NSIG_NEG,'  N(sigma<-1)=',NSIG_NEG1
	  WRITE(LUER,'(A)')
	1     '    chi_hat = chi + |GAM|/dLOG_NU (Knop Paper I Eq 7)'
	  WRITE(LUER,'(A)')
	1     '    signed Q = (GAM/dLOG_NU)/chi_hat in S_hat (Eq 9)'
	  WRITE(LUER,'(A)')
	1     '    no S_tilde (nu-form 3a_nu*I dropped by CMFGEN convention)'
	  WRITE(LUER,'(A)')
	1     '    upwind: I_PREV (GAM>=0) / I_NEXT (GAM<0, lagged GS)'
	  DIAG_PRINTED=.TRUE.
	END IF
!
! edited by Sean Brennan
! Determine maximum ray index having I points. Used in the
! vector section of the DIFFERENCE approach.
!
	  DO I=1,NRAY_MAX
	    DO LS=1,NP
	      IF(NI_RAY(LS) .GE. I)MAX_LS(I)=LS
	    END DO
	  END DO
!
	ELSE IF(NEW_FREQ)THEN
	  FREQ_IDX_ITER=FREQ_IDX_ITER+1
          IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
	    OLDCHI(1:NP)=OLDCHI_STORE(1:NP)
	    AV_PREV(:,:)=AV_STORE(:,:)
	    CV_PREV(:,:)=CV_STORE(:,:)
	  ELSE IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
	    I_P_PREV(:,:)=I_P_STORE(:,:)
	    I_M_PREV(:,:)=I_M_STORE(:,:)
	  END IF
	  IF(FREQ .GE. PREVIOUS_FREQ)THEN
	    WRITE(ERROR_LU(),*)'Error in FG_J_CMF_NONMONO_V1'
	    WRITE(ERROR_LU(),*)'Frequencies must be monotonically decreasng'
	    STOP
	  END IF
	ELSE IF(FREQ .NE. PREVIOUS_FREQ)THEN
	   WRITE(ERROR_LU(),*)'Error in FG_J_CMF_NONMONO_V1'
	   WRITE(ERROR_LU(),*)
	1      'Frequencies must not change if NEW_FREQ=.FALSE.'
	   STOP
	END IF
	PREVIOUS_FREQ=FREQ
!
! 
!
	CALL TUNE(1,'FG_BIG_SEC')
!
! Compute CHI_EXT, and ETA_EXT. CHI_EXT could be saved, as it doesn't change
! during the iteration procedure. ETA does however change (since it depends
! of J) and thus ETA_EXT must be re-computed on each entry.
!
! The first checks whether we may have negative line opacities due
! to stimulated emission. In such a case we simply assume an 1/r^2
! extrapolation.
!
! We also interpolate in ESEC, since ESEC (in the absence of negative
! absorption) provides a lower bound to the opacity. NB: When CHI is much
! larger then ESEC its variation with r dominates, and it is possible to
! extrapolate CHI below ESEC.
!
	CALL TUNE(1,'FG_CHI_BEG')
	IF(ND_ADD .NE. 0)THEN
	  IF(CHI(IDMIN) .LE. ESEC(IDMIN) .OR. CHI(IDMAX) .LE. ESEC(IDMAX))THEN
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0_LDP)ESEC_POW=2.0_LDP
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      dCHIdR(I)= -ESEC_POW*CHI_EXT(I)/R_EXT(I)
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(IDMAX)-ESEC(IDMAX)) / (CHI(IDMIN)-ESEC(IDMIN)) )
	1          /LOG(R(IDMIN)/R(IDMAX))
	    IF(ALPHA .LT. 2.0_LDP)ALPHA=2.0_LDP
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0_LDP)ESEC_POW=2.0_LDP
   	    DO I=1,ND_ADD
	      T1=(CHI(IDMIN)-ESEC(IDMIN))*(R(IDMIN)/R_EXT(I))**ALPHA
	      T2=ESEC(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      CHI_EXT(I)=T1+T2
	      dCHIdR(I)=(-ALPHA*T1-ESEC_POW*T2)/R_EXT(I)
	    END DO
	  END IF
	  DO I=ND_ADD+1,ND_EXT
	    CHI_EXT(I)=CHI(I-ND_ADD)
	  END DO
!
! We limit alpha to 3.5 to avoid excess envelope emission. If alpha were
! 3 we would get a logarithmic flux divergence as we increase the volume.
!
	  ALPHA=LOG(ETA(IDMAX)/ETA(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	  IF(ALPHA .LT. 3.5_LDP)ALPHA=3.5_LDP
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(IDMIN)*(R(IDMIN)/R_EXT(I))**ALPHA
	    IF(ETA_EXT(I) .LE. 1.0E-280_LDP)ETA_EXT(I)=1.0E-280_LDP
	  END DO
	  DO I=ND_ADD+1,ND_EXT
	    ETA_EXT(I)=ETA(I-ND_ADD)
	  END DO
	ELSE
	  CHI_EXT(1:ND)=CHI(1:ND)		!NB: In this can ND=ND_EXT
	  ETA_EXT(1:ND)=ETA(1:ND)
	END IF
!
	LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	CALL TUNE(2,'FG_CHI_BEG')
!
! Per-ray Knop refresh: capture this frequency's sweep inputs.
! FREQ_IDX_ITER is 0 during the first frequency of the sweep (the
! INIT branch resets it and skips the NEW_FREQ increment), so during
! frequency k it equals k-1 and the lagged I_NEXT store holds
! frequency s+1 in slot s.  The KN_* capture arrays instead use
! slot s = frequency s (LK=FREQ_IDX_ITER+1), the convention the
! refresh worker iterates over; the shift back to the lagged-store
! layout happens only inside SEED_FLD/WRITE_STORE.  The
! sign-dependent dLOG_NU columns give expanding cells the backward
! spacing ln(nu_{k-1}/nu_k) = dLOG_NU and compressing cells the
! forward spacing ln(nu_k/nu_{k+1}); FG_NU_NEXT is nu(k+1) set by
! SET_FG_NU_NEXT before each call.  KN_STORE_READY is armed once
! the last frequency of the sweep is captured.
!
	  IF(ALLOCATED(KN_CHI_EXT))THEN
	    LK=FREQ_IDX_ITER+1
	    IF(LK .GE. 1 .AND. LK .LE. KN_NCF)THEN
	      KN_CHI_EXT(1:ND_EXT,LK)=CHI_EXT(1:ND_EXT)
	      KN_ETA_EXT(1:ND_EXT,LK)=ETA_EXT(1:ND_EXT)
	      KN_FREQ(LK)=FREQ
	      KN_DLN_P(LK)=dLOG_NU
	      IF(FG_NU_NEXT .GT. 0.0_LDP .AND.
	1         FG_NU_NEXT .LT. FREQ)THEN
	        KN_DLN_N(LK)=LOG(FREQ/FG_NU_NEXT)
	      ELSE
	        KN_DLN_N(LK)=dLOG_NU
	      END IF
	      KN_DBB(LK)=DBB
	      KN_IC(LK)=IC
	      IF(LK .EQ. KN_NCF)KN_STORE_READY=.TRUE.
	    END IF
	  END IF
!
	JPLUS_IB=0.0_LDP; HPLUS_IB=0.0_LDP; KPLUS_IB=0.0_LDP; NPLUS_IB=0.0_LDP
	JMIN_IB=0.0_LDP;  HMIN_IB=0.0_LDP; KMIN_IB=0.0_LDP;   NMIN_IB=0.0_LDP
	JPLUS_OB=0.0_LDP;  HPLUS_OB=0.0_LDP; KPLUS_OB=0.0_LDP; NPLUS_OB=0.0_LDP
	JMIN_OB=0.0_LDP;   HMIN_OB=0.0_LDP;  KMIN_OB=0.0_LDP;  NMIN_OB=0.0_LDP
!
! 
!***************************************************************************
!***************************************************************************
!
	IF(SOLUTION_METHOD .EQ. 'DIFFERENCE')THEN
!
! Zero AV and CV matrices.
!
	  AV(:,:)=0.0_LDP
	  CV(:,:)=0.0_LDP
!
!	  CALL TUNE(1,'LS_LOOP')
!
	  DO LS=1,NP
	    NI=NI_RAY(LS)
!
! NB: dCHIdR = LOG(CHI)/LOG(R) * CHI/R
!
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      IF(R_RAY(I,LS) .EQ. R_EXT(K))THEN
	        CHI_RAY(I)=CHI_EXT(K)
	        ETA_RAY(I)=ETA_EXT(K)
	        dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/R_RAY(I,LS)
	      ELSE
	        T1=LOG(R_RAY(I,LS)/R_EXT(K))
	        T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	        CHI_RAY(I)=EXP(T2)
	        T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	        ETA_RAY(I)=EXP(T2)
	        T2=(3.0_LDP*CHI_COEF(K,1)*T1+2.0_LDP*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	        dCHIdR_RAY(I)=T2*CHI_RAY(I)/R_RAY(I,LS)
	       END IF
	     END DO
!
	    SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
!
! Even in the THK case we assume that the intensity incident on the
! outer boundary is zero. In THK case we effectively get IBOUND not
! equal zero at the model atmosphere boundary as a consequence of the
! extension.
!
	    IBOUND(LS)=0.0_LDP
!
! 
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum
! calculation for the first frequency.
!
	    IF(INIT)THEN
	      DO I=1,NI
	        Q(I)=0.0_LDP
	        QH(I)=0.0_LDP
	      END DO
	      OLDCHI(LS)=CHI_RAY(NI)
	    ELSE IF(NEW_FREQ)THEN
	      DO I=1,NI-1
	        IF(CHI_RAY(I) .LT. 1.0E-30_LDP)CHI_RAY(I)=1.0E-30_LDP
	        QH(I)=GAMH(I,LS)*2.0_LDP/((CHI_RAY(I)+CHI_RAY(I+1))*dLOG_NU)
	        Q(I)=GAM(I,LS)/(CHI_RAY(I)*dLOG_NU)
	      END DO
	      QH(NI)=0.0_LDP
	      IF(CHI_RAY(NI) .LT. 1.0E-30_LDP)CHI_RAY(NI)=1.0E-30_LDP
	      Q(NI)=GAM(NI,LS)/(CHI_RAY(NI)*dLOG_NU)
	    END IF
!
	    IF(LS .LE. NC)THEN
	      OLDCHI_STORE(LS)=CHI_RAY(NI)
	      DIF_OR_ZF=.FALSE.
	      IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	        T1=0.0_LDP
	        IF(.NOT. INIT)T1=GAM(NI,LS)/(CHI_RAY(NI)*dLOG_NU)	     !Q(NI)
	        DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI_RAY(NI)
	1            *(1.0_LDP+T1*(1.0_LDP-CHI_RAY(NI)/OLDCHI(LS)))
	        DIF_OR_ZF=.TRUE.
	      ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	        DIF_OR_ZF=.TRUE.
	        DBC=0.0_LDP
	        IF(NC .LE. NC_PNT_SRCE)DBC=IC
	      ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	        DIF_OR_ZF=.TRUE.
	        DBC=0.0_LDP
	      END IF
	    END IF
!
	    IF(NEW_FREQ)THEN
!
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is not necessary when
! we re using monotonic interpolation.
!
	      IF(METHOD .EQ. 'ZERO')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
	        END DO
	      ELSE IF(INSERT .OR. METHOD(4:6) .EQ. 'MON')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0_LDP )
	        END DO
	      ELSE
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0_LDP )
	          IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	            DTAU(I,LS)=MAX(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MIN(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          ELSE
	            DTAU(I,LS)=MIN(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MAX(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          END IF
	       END DO
	      END IF
!
	      CALL TUVGHD_RH(TA(1,LS),TB(1,LS),TC(1,LS),
	1            U,VB,VC,GB(1,LS),H(1,LS),XM,
	1            Q,QH,DTAU(1,LS),SOURCE_RAY,DIF_OR_ZF,DBC,IC,LS,NC,NI)
	      XM(1)=-IBOUND(LS)
!
! Update PAR_AV matrix. U, VC, VC, AV_PREV, and CV_PREV do not change on
! successive calls with the same frequency. Thus PAR_AV can be used on
! successive calls with the same frequency.
!
	      PAR_AV(1,LS)=U(1)*AV_PREV(1,LS)
	      DO I=2,NI-1
	         PAR_AV(I,LS)=U(I)*AV_PREV(I,LS)-
	1             (VB(I)*CV_PREV(I-1,LS)+VC(I)*CV_PREV(I,LS))
	      END DO
	      PAR_AV(NI,LS)=U(NI)*AV_PREV(NI,LS)-VB(NI)*CV_PREV(NI-1,LS)
!
	      TC(NI,LS)=0.0_LDP				!As used.
	      DIV(1)=1.0_LDP/(TC(1,LS)+TB(1,LS))
	      TC(1,LS)=TC(1,LS)*DIV(1)
	      TB(1,LS)=TB(1,LS)*DIV(1)
	      DO I=2,NI
	        DIV(I)=1.0_LDP/(TA(I,LS)*TB(I-1,LS)+TB(I,LS)+TC(I,LS))
	        TB(I,LS)=(TA(I,LS)*TB(I-1,LS)+TB(I,LS))*DIV(I)
	        TC(I,LS)=TC(I,LS)*DIV(I)
	      END DO
	      TB(1:NI,LS)=DIV(1:NI)
	    ELSE
!
	      XM(1)=0.0_LDP
	      DO I=2,NI-1
	        XM(I)=-0.5_LDP*SOURCE_RAY(I)*(DTAU(I-1,LS)+DTAU(I,LS))
	      END DO
!
	      IF(LS .GT. NC)THEN
	        XM(NI)=0.5_LDP*DTAU(NI-1,LS)*SOURCE_RAY(NI)
	      ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	        XM(NI)=DBC
	      ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	        XM(NI)=DBC
	      ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	        XM(NI)=0.0_LDP
	      ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	        XM(NI)=IC
	     END IF
!
	    END IF
!
! Update AV matrix.
!
	    DO I=1,NI
	      AV(I,LS)=XM(I)+PAR_AV(I,LS)
	    END DO
!
	  END DO				!LS
!
! 
!
! Solve for the radiation field along ray for this frequency.
!
	  IF(VECTOR_MACHINE)THEN
!
! This section is for a VECTOR machine. Loop over inner index is outer
! loop to remove a dependency. This is inefficient on scaler machines
! as array is not accessed sequentially.
!
!
! Forward substitution.
!
	    AV(1,1:NP)=AV(1,1:NP)*TB(1,1:NP)
	    DO I=2,NRAY_MAX
	      DO LS=1,MAX_LS(I)
	        IF(I .LE. NI_RAY(LS))
	1           AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY(LS),LS)=-AV(NI_RAY(LS),LS)
	    END DO
	    DO I=NRAY_MAX,1,-1
	      DO LS=1,MAX_LS(I)
	        IF(I .LT. NI_RAY(LS))
	1          AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
	      END DO
	    END DO
	  ELSE
C
C This section of code is for a scaler machine where the dependence
C on a previous computation does not matter. More efficient than previous
C code as array is accessed in correct manner.
C
C Forward substitution.
C
	    DO LS=1,NP
	      AV(1,LS)=AV(1,LS)*TB(1,LS)
	      DO I=2,NI_RAY(LS)
	        AV(I,LS)=(AV(I,LS)+TA(I,LS)*AV(I-1,LS))*TB(I,LS)
	      END DO
	    END DO
C
C Backward substitution.
C
	    DO LS=1,NP
	      AV(NI_RAY(LS),LS)=-AV(NI_RAY(LS),LS)
	      DO I=NI_RAY(LS)-1,1,-1
	        AV(I,LS)=TC(I,LS)*AV(I+1,LS)-AV(I,LS)
	      END DO
	    END DO
	  END IF
!	  CALL TUNE(2,'LS_LOOP')
!
!
!
! Verify validity of AV values. If these go negative, we set them +ve
! but a factor of 10 smaller. We also ensure that the AV are not
! extremely close to zero by comparing with neighboring AV values.
! The CV (fluxes) are computed using the revised values.
!
	  DO LS=1,NP
	   NI=NI_RAY(LS)
!
	    K=2
	    DO I=1,NI
	      T1=1.0E-08_LDP*ABS(AV(K,LS))
	      IF(AV(I,LS) .LE. T1)THEN
	        AV(I,LS)=MAX(0.01_LDP*ABS(AV(I,LS)),T1)
	        NEG_AV_VALUE=.TRUE.
	      ELSE
	        K=K+1
	      END IF
	      IF(I .EQ. 1)K=1
	    END DO
!
! Update C vector (i.e. flux variable).
!
	    DO I=1,NI_RAY(LS)-1
	      CV(I,LS)=GB(I,LS)*(AV(I,LS)-AV(I+1,LS))+H(I,LS)*CV_PREV(I,LS)
	    END DO
!
! Note that V=AV(1)-IBOUND.
!
	    IF(ND_ADD .EQ. 0)THEN
	      CV_BOUND(LS)=AV(1,LS)-IBOUND(LS)
	    ELSE IF(LS .EQ. NP)THEN
	      CV_BOUND(LS)=0.0_LDP
	    ELSE
!
! J_PNT(1,LS)=1 is the common case -- R(1) sits on the outermost ray
! point -- so the unbounded average below read CV(0,LS), an
! out-of-bounds element that silently skewed every boundary moment
! under THK_CONT.  CV is only filled on intervals 1:NI_RAY-1, so the
! interpolation falls back to the nearest valid interval at both ends.
!
	      K=J_PNT(1,LS)
	      IF(K .LT. 1 .OR. K .GT. NI_RAY(LS))THEN
	        WRITE(LUER,*)'Invalid outer J_PNT:',LS,K,NI_RAY(LS)
	        STOP
	      END IF
	      IF(K .LE. 1)THEN
	        CV_BOUND(LS)=CV(1,LS)
	      ELSE IF(K .GE. NI_RAY(LS)-1)THEN
	        CV_BOUND(LS)=CV(NI_RAY(LS)-1,LS)
	      ELSE
	        CV_BOUND(LS)=0.5_LDP*(CV(K-1,LS)+CV(K+1,LS))
	      END IF
	    END IF
	    I_M_IN_BND(LS)=AV(NI,LS)-(AV(NI,LS)-AV(NI-1,LS))/DTAU(NI-1,LS)
	  END DO
!
	  AV_STORE(:,:)=AV(:,:)
	  CV_STORE(:,:)=CV(:,:)
!
	END IF			!Solution method
!
!************************************************************************
!************************************************************************
!
	CALL TUNE(1,'FG_SOL_INT')
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
!
! Forward logarithmic spacing ln(nu_k/nu_{k+1}) for the compressing
! upwind coupling: a per-frequency scalar, hoisted out of the ray
! loops where it was re-evaluated per point (3 LOG calls per
! compressing ray point).  FG_NU_NEXT is set per call by
! SET_FG_NU_NEXT; falls back to dLOG_NU when unavailable.
!
	DLN_N_SAV=dLOG_NU
	IF(FG_NU_NEXT .GT. 0.0_LDP .AND. FG_NU_NEXT .LT. FREQ)
	1   DLN_N_SAV=LOG(FREQ/FG_NU_NEXT)
!
! Enter loop to perform integration along each ray.
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP1 PRIVATE(SOURCE_PRIME,SOURCE_RAY,CHI_RAY,ETA_RAY,dCHIdR_RAY,Q,
!$OMP& EE,E0,E1,E2,E3,S,dS,T1,T2,T1_DTAU,I_CORE,dZ,NI,I,K,K2,KC,
!$OMP& KC1,KC2,W_INEXT,INEXT_V,DLN_C,NM_RAY,DLN_RAY,J)
!
	  DO LS=1,NP
!
	    FG_NOSLOT_RAY(LS)=0
	    IF(NI_RAY(LS) .EQ. 1)THEN
	      I_P(1,LS)=0.0_LDP
	      I_M(1,LS)=0.0_LDP
	      NI=1
	      GOTO 1000
	    END IF
!
! NB: dCHIdR = LOG(CHI)/LOG(R) * CHI/R
!
!	    CALL TUNE(1,'FG_CHI_Q')
	    NI=NI_RAY(LS)
	    DO I=1,NI_RAY(LS)
	      K=RAY_PNT(I,LS)
	      IF(R_RAY(I,LS) .EQ. R_EXT(K))THEN
	        CHI_RAY(I)=CHI_EXT(K)
	        ETA_RAY(I)=ETA_EXT(K)
	        dCHIdR_RAY(I)=CHI_COEF(K,3)*CHI_RAY(I)/R_RAY(I,LS)
	      ELSE
	        T1=LOG(R_RAY(I,LS)/R_EXT(K))
	        T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	        CHI_RAY(I)=EXP(T2)
	        T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	        ETA_RAY(I)=EXP(T2)
	        T2=(3.0_LDP*CHI_COEF(K,1)*T1+2.0_LDP*CHI_COEF(K,2))*T1+CHI_COEF(K,3)
	        dCHIdR_RAY(I)=T2*CHI_RAY(I)/R_RAY(I,LS)
	      END IF
	    END DO
!
! By setting PF(1)=0 when evaluating SOURCE we ensure a pure continuum
! calculation for the first frequency.
!
!        Knop, Hauschildt & Baron 2009 formalism:
!          Paper I:  A&A 496, 295-298 (arXiv:0902.1876)
!          Paper II: A&A 501, 813-820 (arXiv:0909.4264)
!        The wavelength boundary-value problem and the sign-dependent
!        local upwind discretisation originate in
!          Baron & Hauschildt 2004, A&A 427, 987-1000
!          (arXiv:astro-ph/0408212)
!        which Paper II itself calls "Paper I"; throughout this file
!        Paper I/II always means the two Knop et al. 2009 papers above.
!
!        The scheme is the fully implicit xi=1 member of the
!        Hauschildt & Baron (2004) xi-mixed discretisation: Paper I
!        found xi=1 required for alternating velocity fields.
!
!        CMFGEN uses the FREQUENCY (log-nu) comoving-frame formulation.
!        The transfer equation along a ray is:
!
!          dI/ds = eta - chi*I + GAM * dI/d(ln nu)
!
!        where GAM = (V/R)*(1 + SIGMA*MU^2) is the signed comoving-frame
!        frequency gradient (V in km/s, R in cm, 3.33564E-06 converts
!        to cgs).  The frequency derivative is discretised as:
!
!          GAM * (I_upwind - I_local) / dLOG_NU
!
!        Splitting into implicit (opacity) and explicit (source) parts:
!
!          chi_hat = chi + |GAM|/dLOG_NU    (generalized opacity, always > 0)
!          Q      = (GAM/dLOG_NU) / chi_hat (dimensionless, SIGNED)
!          S'     = eta/chi_hat + Q * I_upwind
!
!        so the transfer equation becomes  dI/dtau_hat = I - S',
!        with tau_hat = chi_hat * dz.  chi_hat absorbs the IMPLICIT
!        part of the coupling (the I_local term), while Q*I_upwind
!        handles the EXPLICIT part (coupling to the upwind frequency).
!
!        IMPORTANT: Knop's Paper I is written in the WAVELENGTH
!        formulation where the transfer equation has an additional
!        explicit 4*a*I term (a = V/(cR)) from spherical divergence:
!          dI/ds = eta - chi*I - 4*a*I - a * d(lambda*I)/d(lambda)
!        In that formulation, removing 4*a from the opacity requires
!        adding it back as an explicit source S_tilde = 4*a*I/chi_hat
!        (Paper I Eq 10).
!
!        Status of the 4*a*I term in the frequency formulation:
!        Knop's 4*a_lambda*I is the intensity-proportional part of the
!        comoving-frequency coupling exposed by the d(lambda*I)/d(lambda)
!        split.  It does NOT vanish for I_nu -- only the invariant
!        intensity I_nu/nu^3 (Chen et al. 2007, Eq 10) is free of a
!        proportional term.  The exact O(v/c) ray equation for I_nu is
!          dI/ds = eta - chi*I - a_nu*( dI/d(ln nu) - 3*I )
!        with a_nu = d(ln nu)/ds = -a_lambda.  Written Knop-style the
!        same coupling is +4*a_nu*I - a_nu*d(nu*I)/dnu; for the
!        dI/d(ln nu) split used here the nu-form '4aI' coefficient is 3.
!
!        CMFGEN's convention drops it on the ray: the equation solved
!        here and in stock FG_J_CMF_V13 is
!          dI/ds = eta - chi*I + GAM * dI/d(ln nu)
!        with GAM = (V/cR)(1 + SIGMA*mu^2) the signed O(v/c)
!        frequency-drift coefficient.  This is the inherited CMFGEN
!        truncation, not an omission introduced by this routine; for
!        GAM > 0 this file reduces to V13 exactly.  At the moment level
!        the dilution physics enters through the advection terms of
!        MOM_J_DDT ((1/c r^3) D(r^3 J)/Dt contains 3*v*J/(c r), the
!        exact sigma=0 counterpart) and through the relativistic
!        coefficient structure of MOM_JREL (the gamma^2*(sigma+1)-1
!        coefficients, the H-moment aberration term, and the
!        INCL_ADVEC_TERMS pieces folded into CHI_J).  The
!        (1+SIGMA*K_ON_J) factors are the angular weight of the drift
!        coefficient (1+SIGMA*mu^2); they belong to the derivative
!        term, not to the proportional term.
!
!        The S_tilde removed on 07-Sep-2026 transplanted the lambda-form
!        4*a*I into this nu-form scheme: the literal '4' is specific to
!        the d(lambda*I)/d(lambda) split and the term double-counted
!        the retained drift coupling, so the removal was correct.  If
!        exact Chen/Knop consistency is ever required, the term to add
!        is +3*a_nu*I with the same sign-aware split as the drift term
!        (absorption-side in expansion, source-side in compression);
!        it would change results vs stock V13 and is deliberately not
!        added here.
!        (Sean Brennan, 07-Sep-2026; justification expanded
!         09-Oct-2026 after external review)
!
!        In non-monotonic flows GAM can be negative.  The Knop fix
!        (Paper I Eq 7) is to use |GAM| in the opacity (keeping it
!        positive) while retaining the SIGNED GAM in the source via Q.
!
!        For expanding flow (GAM > 0):  chi_hat = chi + GAM/dLOG_NU,
!        Q > 0 -- identical to V13.  Upwind frequency = previous (I_PREV).
!        For compressing flow (GAM < 0): chi_hat = chi + |GAM|/dLOG_NU
!        (larger than chi), Q < 0 (signed).  Upwind frequency = NEXT
!        (I_NEXT), because the comoving frequency reverses direction.
!
!        The sign-dependent upwind is the KEY difference from V13:
!        compressing points couple to I_NEXT (the next lower frequency),
!        not I_PREV.  Since CMFGEN sweeps frequencies in decreasing
!        order, I_NEXT is not yet computed in the current sweep.  I_NEXT
!        therefore comes from the PREVIOUS great iteration, stored at
!        depth grid points (I_P_DEPTH_NEXT, I_M_DEPTH_NEXT) -- a
!        Gauss-Seidel lag across great iterations.  This lagged coupling
!        is a CMFGEN adaptation; Paper II converges the same wavelength
!        boundary-value problem inside the formal solution with its
!        physically guided GS stepping (Sect. 3.2, Eqs 25-28), which
!        FG_KNOP_REFRESH_V1 implements as the post-sweep re-solve.
!
!        On the first great iteration, I_NEXT = 0, so compressing
!        points get S' = S (no coupling).  This is the cold-start
!        approximation; subsequent iterations refine I_NEXT.
!
!        No Q-splitting or Q limiter is needed -- the implicit
!        absorption into chi_hat provides natural stability, and the
!        signed Q with correct upwind direction provides the correct
!        physics.  An intensity floor (I >= 1.0E-30) is applied after
!        the short-characteristics integration to catch unphysical
!        negative/NaN overshoot from the explicit source term in
!        extreme compression zones.
!
!        GAM_MAG(I,LS) = |GAM(I,LS)| is already computed above.
! edited by Sean Brennan (07-Sep-2026)
!
!        Compressing cells (GAM<0) couple to the NEXT (redder)
!        frequency, so the frequency derivative spanning (k,k+1) must
!        use the forward spacing ln(nu_k/nu_{k+1}), not the backward
!        dLOG_NU = ln(nu_{k-1}/nu_k).  At line edges the grid spacing
!        changes by ~10x and the wrong spacing mis-weights the Knop
!        coupling there.  FG_NU_NEXT (nu_{k+1}) is set per call by
!        SET_FG_NU_NEXT; when unavailable the old dLOG_NU is used.
!        (Oct-2026 asymmetric-spacing fix)
!
!        Matched frequency offsets (VADAT DO_NONMONO_FREQ_MATCH):
!        the compressing upwind is rematched from the adjacent slot
!        to the slot containing the characteristic foot, n =
!        |GAM|*dZ/DLN bins along the ray, and the frequency
!        derivative spans the matched log range ln(nu_k/nu_{k+n}).
!        Defaults reproduce the pre-match stencil: slot +1 and the
!        forward spacing DLN_N_SAV.  Set before the INIT split so the
!        bluest frequency can use the same stencil -- at INIT only
!        compressing cells consume it (expanding cells have no bluer
!        upwind and keep Q=0).
!
	      DO I=1,NI
	        NM_RAY(I)=1
	        DLN_RAY(I)=DLN_N_SAV
	      END DO
!        FREQ_IDX_ITER+1 < NCF_IN: the compressing upwind is the
!        redder slot, so the last real frequency has nothing to
!        rematch to (slot NCF_IN is the phantom).
!
	      IF(FG_FREQ_MATCH .AND. FREQ_IDX_ITER+1 .LT. NCF_IN .AND.
	1       DLN_N_SAV .GT. 0.0_LDP)THEN
	        DO I=1,NI
	          IF(GAM(I,LS) .LT. 0.0_LDP)THEN
	            IF(I .EQ. 1)THEN
	              T1=ABS(Z(1,LS)-Z(2,LS))
	            ELSE IF(I .EQ. NI)THEN
	              T1=ABS(Z(NI-1,LS)-Z(NI,LS))
	            ELSE
	              T1=0.5_LDP*(ABS(Z(I-1,LS)-Z(I,LS))+
	1            ABS(Z(I,LS)-Z(I+1,LS)))
	            END IF
	            T2=GAM_MAG(I,LS)*T1/DLN_N_SAV
!
! Cap K2 at the last REAL lagged-store slot.  During frequency k the
! store index FREQ_IDX_ITER+K2 holds nu_{k+K2} (slot s = frequency
! s+1), so the largest usable offset is NCF_IN-FREQ_IDX_ITER-1: slot
! NCF_IN is a force-zeroed phantom beyond the last frequency.
! (Oct-2026 matched-span fix)
!
	            IF(T2 .GE. DBLE(NCF_IN-FREQ_IDX_ITER-1))THEN
	              K2=NCF_IN-FREQ_IDX_ITER-1
	            ELSE
!
! CEILING keeps the characteristic foot inside the matched
! stencil; NINT could round down and leave the residual
! extrapolating up to half a bin past the foot.
!
	              K2=MAX(1,CEILING(T2))
	            END IF
	            NM_RAY(I)=K2
	            DLN_RAY(I)=K2*DLN_N_SAV
!
! KN_FREQ is indexed by frequency (slot s = frequency s), while the
! lagged I_NEXT store is shifted by one (slot s = frequency s+1).
! The upwind read uses store slot FREQ_IDX_ITER+K2 = frequency
! k+K2, so the matched log span must close on KN_FREQ index
! FREQ_IDX_ITER+1+K2 -- the same bin MATCH_FREQ uses in the refresh
! worker (its NFM offset L+NNF in KN_ indexing).  The previous
! FREQ_IDX_ITER+K2 index measured the span to nu_{k+K2-1}, so the
! in-sweep frequency derivative and the refresh operator disagreed
! by one bin.  (Oct-2026 matched-span fix)
!
	            IF(ALLOCATED(KN_FREQ) .AND.
	1           FREQ_IDX_ITER+1+K2 .LE. NCF_IN)THEN
	              IF(KN_FREQ(FREQ_IDX_ITER+1+K2) .GT. 0.0_LDP .AND.
	1             KN_FREQ(FREQ_IDX_ITER+1+K2) .LT. FREQ)
	1           DLN_RAY(I)=
	1             LOG(FREQ/KN_FREQ(FREQ_IDX_ITER+1+K2))
	            END IF
	          END IF
	        END DO
	      END IF
!
	    IF(INIT)THEN
	      Q(1:NI)=0.0_LDP
!        INIT path: guard chi and eta against non-finite values from
!        interpolation.  The .NOT.(|x| <= HUGE) test catches +Inf,
!        which the old x /= x and x < 1e-30 tests both missed.
! edited by Sean Brennan (04-Sep-2026)
	      DO I=1,NI
	        IF(.NOT.(ABS(CHI_RAY(I)) .LE. HUGE(CHI_RAY(I))) .OR.
	1          CHI_RAY(I) .LT. 1.0E-30_LDP)CHI_RAY(I)=1.0E-30_LDP
	        IF(.NOT.(ABS(ETA_RAY(I)) .LE. HUGE(ETA_RAY(I))))
	1          ETA_RAY(I)=0.0_LDP
	        IF(.NOT.(ABS(GAM(I,LS)) .LE. HUGE(GAM(I,LS))))
	1          GAM(I,LS)=0.0_LDP
	        IF(.NOT.(ABS(GAM_MAG(I,LS)) .LE. HUGE(GAM_MAG(I,LS))))
	1          GAM_MAG(I,LS)=0.0_LDP
	        IF(.NOT.(ABS(dGAMdR(I,LS)) .LE. HUGE(dGAMdR(I,LS))))
	1          dGAMdR(I,LS)=0.0_LDP
	      END DO
!
!        Compressing cells at the bluest frequency: the upwind bin is
!        the REDDER lagged slot (store slot 1 = frequency 2), which
!        exists -- only expanding cells have no upwind here (blue-edge
!        inflow is zero).  Forcing Q=0 on compressing cells drops the
!        advective coupling at the band edge and leaves photons unable
!        to leave the bluest bin; apply the same chi_hat/Q treatment
!        as interior frequencies.  On the first great iteration the
!        lagged store is zero, so the coupling acts as pure outflow --
!        the same cold-start approximation the interior bins use.
!        NM_RAY/DLN_RAY are set above the INIT split, so the update is
!        valid without the HAS_COMPRESSING store census.
!
	      DO I=1,NI
	        IF(GAM(I,LS) .LT. 0.0_LDP)THEN
	          DLN_C=DLN_RAY(I)
	          IF(.NOT.(DLN_C .GT. 0.0_LDP))DLN_C=1.0E-30_LDP
	          T1=GAM_MAG(I,LS)/DLN_C
	          CHI_RAY(I)=CHI_RAY(I)+T1
	          dCHIdR_RAY(I)=dCHIdR_RAY(I)-dGAMdR(I,LS)/DLN_C
	        END IF
	      END DO
!        Post-update guards on the compressing increments: a cubic-
!        overshoot GAM<0 can drive chi_hat to +Inf (uncaught by the
!        entry test above) and leave DTAU infinite downstream.
	      DO I=1,NI
	        IF(.NOT.(ABS(CHI_RAY(I)) .LE. HUGE(CHI_RAY(I))) .OR.
	1          CHI_RAY(I) .LT. 1.0E-30_LDP)CHI_RAY(I)=1.0E-30_LDP
	        IF(.NOT.(ABS(dCHIdR_RAY(I)) .LE. HUGE(dCHIdR_RAY(I))))
	1          dCHIdR_RAY(I)=0.0_LDP
	      END DO
	      DO I=1,NI
	        IF(GAM(I,LS) .LT. 0.0_LDP)THEN
	          DLN_C=DLN_RAY(I)
	          IF(.NOT.(DLN_C .GT. 0.0_LDP))DLN_C=1.0E-30_LDP
	          Q(I)=GAM(I,LS)/(DLN_C*CHI_RAY(I))
	        END IF
	      END DO
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
!        Final belt-and-suspenders guard, matching the ELSE branch.
	      DO I=1,NI
	        IF(.NOT.(ABS(Q(I)) .LE. HUGE(Q(I))))Q(I)=0.0_LDP
	        IF(.NOT.(ABS(SOURCE_RAY(I)) .LE. HUGE(SOURCE_RAY(I))))
	1          SOURCE_RAY(I)=0.0_LDP
	      END DO
	    ELSE
!        Generalized opacity: chi_hat = chi + |GAM|/dLOG_NU (always > 0)
!        NM_RAY/DLN_RAY are set above the INIT split.
! edited by Sean Brennan
!
	      DO I=1,NI
	        DLN_C=dLOG_NU
	        IF(GAM(I,LS) .LT. 0.0_LDP)DLN_C=DLN_RAY(I)
!
! DLN_C is a log-frequency spacing; a degenerate grid (duplicate
! frequencies) or an unseeded DLN_RAY slot can make it zero or
! negative, and 1/DLN_C would then poison chi_hat.  Flooring it
! keeps chi_hat well above CHI_RAY -- locally opaque, the safe
! failure direction.
!
	        IF(.NOT.(DLN_C .GT. 0.0_LDP))DLN_C=1.0E-30_LDP
	        Q(I)=GAM_MAG(I,LS)/DLN_C
	        CHI_RAY(I)=CHI_RAY(I)+Q(I)
	      END DO
!        Guard: ensure chi_hat is strictly positive to avoid 0/0 NaN.
!        At low-opacity boundary points both chi and |GAM| can vanish.
!
!        CRITICAL: The .LT. test alone does NOT catch NaN.  In IEEE-754,
!        NaN .LT. x is always FALSE, so a NaN chi_hat would pass through
!        this guard, making Q, SOURCE_RAY, and DTAU all NaN.
!        The NaN then propagates into the SE rates and STEQ.
!        The .NOT.(|x| <= HUGE) test catches NaN and Inf.
!
! edited by Sean Brennan (04-Sep-2026)
	      DO I=1,NI
	        IF(.NOT. (ABS(CHI_RAY(I)) .LE. HUGE(CHI_RAY(I))) .OR.
	1          CHI_RAY(I) .LT. 1.0E-30_LDP)THEN
!
! Phase 2 first-occurrence diagnostic for chi_hat failure.
! Edited 30-Sep-2026: the test now runs BEFORE the clamp so the
! reported CHI_RAY is the raw failure -- the unconditional clamp
! previously sat above this branch and made it unreachable.  The
! .NOT.(|x| <= HUGE) test also catches Inf, which x /= x misses.
!$OMP CRITICAL(FIRST_FAIL_CRIT)
                IF(.NOT. FIRST_FAIL_FOUND)THEN
                  FIRST_FAIL_FOUND=.TRUE.
                  FIRST_FAIL_TYPE=4
                  J=ERROR_LU()
                  WRITE(J,'(/,A)')' ================== PHASE 2: FIRST NUMERICAL FAILURE (chi_hat) ================='
                  WRITE(J,*)"Type: chi_hat NaN or < 1e-30 (before clamping)"
                  WRITE(J,"(A,I8)")"  Freq index: ",FREQ_IDX_ITER
                  WRITE(J,"(A,ES16.8)")"  Frequency: ",FREQ
                  WRITE(J,"(A,I6)")"  Ray (LS): ",LS
                  WRITE(J,"(A,I6)")"  Ray point (I): ",I
                  WRITE(J,"(A,ES16.8)")"  CHI_RAY (chi_hat, pre-clamp): ",CHI_RAY(I)
                  WRITE(J,"(A,ES16.8)")"  GAM: ",GAM(I,LS)
                  WRITE(J,"(A,ES16.8)")"  GAM_MAG: ",GAM_MAG(I,LS)
                  WRITE(J,"(A,ES16.8)")"  dLOG_NU: ",dLOG_NU
                  WRITE(J,"(A,ES16.8)")"  R_RAY: ",R_RAY(I,LS)
                  WRITE(J,"(A,ES16.8)")"  V_RAY: ",V_RAY(I)
                  WRITE(J,"(A,ES16.8)")"  SIGMA_RAY: ",SIGMA_RAY(I)
                  WRITE(J,"(A,ES16.8)")"  ETA_RAY: ",ETA_RAY(I)
                  WRITE(J,'(/,A)')' ========================== END PHASE 2 FIRST FAILURE =========================='
                END IF
!$OMP END CRITICAL(FIRST_FAIL_CRIT)
	          CHI_RAY(I)=1.0E-30_LDP
	        END IF
	      END DO
!        Update dCHIdR_RAY to include the derivative of chi_hat.
!        d(chi_hat)/dR = dchi/dR + d|GAM|/dR / dLOG_NU
!        d|GAM|/dR = SIGN(1,GAM)*dGAM/dR.
!        This matches V13 where GAM is non-negative in monotonic flow.
!        The DTAU clamping below handles any cusps at velocity turning
!        points where d|GAM|/dR is discontinuous.
! edited by Sean Brennan
	      DO I=1,NI
	        DLN_C=dLOG_NU
	        IF(GAM(I,LS) .LT. 0.0_LDP)DLN_C=DLN_RAY(I)
	        dCHIdR_RAY(I)=dCHIdR_RAY(I)+
	1           SIGN(1.0_LDP,GAM(I,LS))*dGAMdR(I,LS)/DLN_C
	      END DO
!        Dimensionless signed coupling: Q = (GAM/dLOG_NU) / chi_hat
!
!        NaN guards on ETA_RAY and GAM: if the cubic interpolation
!        produced NaN (e.g. from extreme population values or
!        coefficient overflow), the source and Q would be NaN.
!        Replace NaN with zero - this gives a zero source contribution,
!        which is physically neutral (no emission, no coupling).
!
! edited by Sean Brennan (04-Sep-2026)
	      DO I=1,NI
	        IF(ETA_RAY(I) .NE. ETA_RAY(I))ETA_RAY(I)=0.0_LDP
	        IF(GAM(I,LS) .NE. GAM(I,LS))GAM(I,LS)=0.0_LDP
	        IF(GAM_MAG(I,LS) .NE. GAM_MAG(I,LS))GAM_MAG(I,LS)=0.0_LDP
	        IF(dGAMdR(I,LS) .NE. dGAMdR(I,LS))dGAMdR(I,LS)=0.0_LDP
	        IF(dCHIdR_RAY(I) .NE. dCHIdR_RAY(I))dCHIdR_RAY(I)=0.0_LDP
	      END DO
	      DO I=1,NI
	        DLN_C=dLOG_NU
	        IF(GAM(I,LS) .LT. 0.0_LDP)DLN_C=DLN_RAY(I)
	        Q(I)=GAM(I,LS)/(DLN_C*CHI_RAY(I))
	      END DO
	      SOURCE_RAY(1:NI)=ETA_RAY(1:NI)/CHI_RAY(1:NI)
!        Final guard: if Q or SOURCE_RAY is still NaN (shouldn't happen
!        after the above guards, but belt-and-suspenders), zero them.
	      DO I=1,NI
	        IF(Q(I) .NE. Q(I))Q(I)=0.0_LDP
	        IF(SOURCE_RAY(I) .NE. SOURCE_RAY(I))SOURCE_RAY(I)=0.0_LDP
	      END DO
!
! edited by Sean Brennan
	    END IF
!	    CALL TUNE(2,'FG_CHI_Q')
!
!
!
!	    CALL TUNE(1,'FG_NEW_F')
	    IF(NEW_FREQ)THEN
!
! Compute the optical depth increments. This code is from TAU, and NORDTAU. We
! check that the Euler-Mauclarin correction is not too large. This is mainly
! done to prevent negative optical depths. The check is still necessary when
! we re using monotonic interpolation, since the monotonic interpolation only
! applies to CHI. CHI_RAY contains an additional term due to the frequency
! derivative. As non relativistic, DTAU is the same for both directions.
!
	      IF(METHOD .EQ. 'ZERO')THEN
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
	        END DO
	      ELSE
	        DO I=1,NI-1
	          dZ=Z(I,LS)-Z(I+1,LS)
	          DTAU(I,LS)=0.5_LDP*dZ*( CHI_RAY(I)+CHI_RAY(I+1) +
	1            dZ*( dCHIdR_RAY(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1            dCHIdR_RAY(I)*Z(I,LS)/R_RAY(I,LS) )/6.0_LDP )
	          IF( CHI_RAY(I) .LT. CHI_RAY(I+1) )THEN
	            DTAU(I,LS)=MAX(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MIN(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          ELSE
	            DTAU(I,LS)=MIN(CHI_RAY(I)*dZ,DTAU(I,LS))
	            DTAU(I,LS)=MAX(CHI_RAY(I+1)*dZ,DTAU(I,LS))
	          END IF
	        END DO
	      END IF
!
!        DTAU NaN guard: ensure DTAU is non-negative and finite.
!        With the chi_hat and dCHIdR_RAY guards above, DTAU should
!        always be finite, but the Euler-Maclaurin correction can
!        still produce small negative values at discontinuities.
!
!        When the corrected DTAU is invalid (NaN or negative), fall back
!        to the guaranteed-positive trapezoidal estimate
!          0.5*(CHI_RAY(I)+CHI_RAY(I+1))*dZ
!        (the same formula used in the METHOD.EQ.'ZERO' branch) rather
!        than an arbitrary 1.0E-10 constant.  CHI_RAY is floored at
!        1.0E-30 and dZ=Z(I)-Z(I+1)>0, so this fallback is always
!        positive and preserves the actual opacity magnitude of the
!        segment.  At a sharp velocity-gradient reversal (large |GAM|,
!        extreme negative SIGMA) the true optical-depth increment can
!        be large via the |GAM|/dLOG_NU term; substituting a near-zero
!        constant there would falsely mark the segment as optically
!        thin and corrupt the short-characteristics intensity for the
!        rest of that ray.
! edited by Sean Brennan (04-Sep-2026)
      DO I=1,NI-1
        IF(DTAU(I,LS) .NE. DTAU(I,LS) .OR.
     1     DTAU(I,LS) .LT. 0.0_LDP)THEN
!
! Phase 2 first-occurrence diagnostic for DTAU failure.
!$OMP CRITICAL(FIRST_FAIL_CRIT)
          IF(.NOT. FIRST_FAIL_FOUND)THEN
            FIRST_FAIL_FOUND=.TRUE.
            FIRST_FAIL_TYPE=6
            J=ERROR_LU()
            WRITE(J,'(/,A)')' =================== PHASE 2: FIRST NUMERICAL FAILURE (DTAU) ==================='
            WRITE(J,*)"Type: DTAU NaN or negative (before fallback)"
            WRITE(J,"(A,I8)")"  Freq index: ",FREQ_IDX_ITER
            WRITE(J,"(A,ES16.8)")"  Frequency: ",FREQ
            WRITE(J,"(A,I6)")"  Ray (LS): ",LS
            WRITE(J,"(A,I6)")"  Ray point (I): ",I
            WRITE(J,"(A,ES16.8)")"  DTAU (pre-fallback): ",DTAU(I,LS)
            WRITE(J,"(A,ES16.8)")"  CHI_RAY(I): ",CHI_RAY(I)
            WRITE(J,"(A,ES16.8)")"  CHI_RAY(I+1): ",CHI_RAY(I+1)
            WRITE(J,"(A,ES16.8)")"  dCHIdR_RAY(I): ",dCHIdR_RAY(I)
            WRITE(J,"(A,ES16.8)")"  dCHIdR_RAY(I+1): ",dCHIdR_RAY(I+1)
            WRITE(J,"(A,ES16.8)")"  dZ: ",Z(I,LS)-Z(I+1,LS)
            WRITE(J,"(A,ES16.8)")"  GAM(I): ",GAM(I,LS)
            WRITE(J,"(A,ES16.8)")"  GAM(I+1): ",GAM(I+1,LS)
            WRITE(J,"(A,ES16.8)")"  dGAMdR(I): ",dGAMdR(I,LS)
            WRITE(J,"(A,ES16.8)")"  dLOG_NU: ",dLOG_NU
            WRITE(J,'(/,A)')' ========================== END PHASE 2 FIRST FAILURE =========================='
          END IF
!$OMP END CRITICAL(FIRST_FAIL_CRIT)
          DTAU(I,LS)=0.5_LDP*(CHI_RAY(I)+CHI_RAY(I+1))*
     1               (Z(I,LS)-Z(I+1,LS))
        END IF
      END DO
!
!        DTAU is now computed from chi_hat = chi + |GAM|/dLOG_NU, which
!        is always positive.  No scaling is needed - the
!        generalized opacity already contains the implicit velocity
!        coupling.  This is identical to the V13 monotonic solver but
!        with |GAM| instead of GAM in the opacity.
!
! edited by Sean Brennan
! Compute the functions used to evaluate the weighted integral over the
! polynomial fit to the source function.
!
! NB:  En(I)= EXP(-DTAU) {Integral[0 to DTAU] t^n EXP(t) dt }/ DTAU^n
!
	      DO I=1,NI-1
	        T1=DTAU(I,LS)
	        EE(I)=0.0_LDP
	        IF(T1 .LT. 700.0_LDP)EE(I)=EXP(-T1)
	        IF(T1 .GE. 40.0_LDP)THEN
	        ELSE IF(T1 .GT. 0.5_LDP)THEN
	          E0(I)=1.0_LDP-EE(I)
	          E1(I)=1.0_LDP-E0(I)/T1
	          E2(I)=1.0_LDP-2.0_LDP*E1(I)/T1
	          E3(I)=1.0_LDP-3.0_LDP*E2(I)/T1
	        ELSE IF(T1 .GT. 0.1_LDP)THEN
	          E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1               (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1               (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP*
	1               (1.0_LDP-T1/10.0_LDP*(1.0_LDP-T1/11.0_LDP*
	1               (1.0_LDP-T1/12.0_LDP*(1.0_LDP-T1/13.0_LDP)))))))) )
	          E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	          E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	          E0(I)=T1*( 1.0_LDP-E1(I) )
	        ELSE
	          E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1               (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1               (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP) ))))
	          E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	          E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	          E0(I)=T1*( 1.0_LDP-E1(I) )
	        END IF
	      END DO
!
	      DO I=1,NI-1
	        T1=DTAU(I,LS)
	        IF(T1 .GE. 40.0_LDP)THEN
	          A0(I,LS)=EE(I)
	          A1(I,LS)=(6.0_LDP-12.0_LDP/T1)/T1/T1
	          A2(I,LS)=1.0_LDP-A1(I,LS)
	          A3(I,LS)=(2.0_LDP-6.0_LDP/T1)/T1
	          A4(I,LS)=(4.0_LDP-6.0_LDP/T1)/T1-1.0_LDP
	        ELSE
	          A0(I,LS)=EE(I)
	          A1(I,LS)=E0(I)-3.0_LDP*E2(I)+2.0_LDP*E3(I)
	          A2(I,LS)=3.0_LDP*E2(I)-2.0_LDP*E3(I)
	          A3(I,LS)=DTAU(I,LS)*(E1(I)-2.0_LDP*E2(I)+E3(I))
	          A4(I,LS)=DTAU(I,LS)*(E3(I)-E2(I))
	        END IF
	      END DO
!
! S_tilde correction REMOVED (Sean Brennan, 07-Sep-2026).
!
! The previous implementation included an explicit source term
!   S_tilde = 4 * Q * dLOG_NU * I_upwind
! subtracted from the formal solution, following Knop et al. 2009
! Paper I Eq 10.  This is CORRECT in Knop's WAVELENGTH formulation
! where the transfer equation (Eq 1) contains an explicit 4*a*I term:
!   dI/ds = eta - chi*I - 4*a*I - a * d(lambda*I)/d(lambda)
! The 4*a*I term arises from expanding d(lambda*I)/d(lambda) and from
! relativistic aberration corrections (Chen et al. 2007).  Knop's
! S_tilde compensates for removing this 4*a*I from the opacity.
!
! CMFGEN solves the FREQUENCY (log-nu) form and, by convention, drops
! the I_nu proportional term (+3*a_nu*I): the equation differenced
! here and in stock FG_J_CMF_V13 is
!   dI/ds = eta - chi*I - GAM * dI/d(ln nu)
! The I_lambda 4*a*I cannot be carried over as 4*Q*dLOG_NU*I -- wrong
! coefficient for I_nu and a double-count of the drift term.  The
! proportional physics enters at the moment level via the advection
! terms (MOM_J_DDT) and the relativistic coefficient structure
! (MOM_JREL).  See the full justification in the comment block above
! the outward rays formal solution.
!
! The coupling is discretised as
!   GAM * (I_upwind - I_l) / dLOG_NU
! and enters the opacity as chi_hat = chi + |GAM|/dLOG_NU and the
! source as S = eta/chi_hat + Q*I_upwind with signed Q.
!
! Adding S_tilde = 4*Q*dLOG_NU*I in the frequency formulation
! DOUBLE-COUNTS the velocity coupling.  For monotonic flow (GAM > 0)
! the non-monotonic solver would NOT reduce to V13, producing
! inconsistent results.  The S_tilde, ALPHA_L, and BETA_L arrays
! and their computation have been removed.
!
! edited by Sean Brennan (07-Sep-2026)
	    END IF
!
!	    CALL TUNE(2,'FG_NEW_F')
!
! ******************* INWARD DIRECTED RAYS *********************************
!
! Compute the Source function for inward directed rays, and find the
! monotonic interpolating polynomial.
!
!	    CALL TUNE(1,'FG_INT_INT')
!
!        Source function for inward directed rays.
!        Full Knop et al. 2009: sign-dependent upwind frequency.
!        Expanding (GAM >= 0): upwind = previous frequency = I_M_PREV
!          (already computed in the current decreasing-frequency sweep)
!        Compressing (GAM < 0): upwind = next frequency = I_M_DEPTH_NEXT
!          (from the previous great iteration, Gauss-Seidel approach)
!        On the first great iteration, I_M_DEPTH_NEXT is zero, so
!        compressing points get SOURCE_PRIME = SOURCE_RAY (no coupling).
!        This converges to the correct I_NEXT over CMFGEN great iterations.
!
! edited by Sean Brennan
	    DO I=1,NI
	      IF(GAM(I,LS) .GE. 0.0_LDP)THEN
	        SOURCE_PRIME(I)=SOURCE_RAY(I)+Q(I)*I_M_PREV(I,LS)
	      ELSE IF(HAS_COMPRESSING)THEN
	        K=RAY_PNT(I,LS)-ND_ADD
!
!	        Linear interpolation of the lagged I_NEXT between the two
!	        radial depths bracketing the ray point (audit 17-Sep-2026).
!	        The nearest-depth lookup used previously introduced a
!	        spatial discontinuity exactly where the velocity gradient is
!	        strongest.  When only one bracketing depth is a compressing
!	        (stored) depth, that value is used -- the intensity at the
!	        next frequency is only stored at compressing depths, so a
!	        non-compressing neighbour has no I_NEXT to interpolate with.
!
	        KC1=0
	        KC2=0
	        IF(K .GE. 1 .AND. K .LE. ND)KC1=COMP_DEPTH_IDX(K)
	        IF(K+1 .GE. 1 .AND. K+1 .LE. ND)KC2=COMP_DEPTH_IDX(K+1)
	        IF((KC1 .GE. 1 .OR. KC2 .GE. 1) .AND.
	1       FREQ_IDX_ITER+NM_RAY(I) .GE. 1 .AND.
	2       FREQ_IDX_ITER+NM_RAY(I) .LE. NCF_IN)THEN
!
!	          GAM < 0 here, so the upwind point in frequency has already
!	          been switched to the NEXT (redder) frequency.  The upwind
!	          coefficient must therefore be the MAGNITUDE |GAM|/dLOG_NU
!	          divided by chi_hat, matching the +|GAM|/dLOG_NU that
!	          CHI_RAY already carries.  Using the signed Q here would
!	          flip the sign a second time and subtract the upwind
!	          intensity from the source, driving SOURCE_PRIME (and hence
!	          I) negative throughout every compression zone.
!	          edited by Sean Brennan
!
!	          NM_RAY(I) is the matched slot offset (1 = adjacent bin;
!	          DO_NONMONO_FREQ_MATCH may rematch it to the bin holding
!	          the characteristic foot).
!
	          IF(KC1 .GE. 1 .AND. KC2 .GE. 1)THEN
	            IF(ABS(R(K+1)-R(K)) .GT. 1.0E-14_LDP*
	1          MAX(ABS(R(K)),ABS(R(K+1)),R(1)))THEN
	              W_INEXT=(R_RAY(I,LS)-R(K))/(R(K+1)-R(K))
	              W_INEXT=MAX(0.0_LDP,MIN(1.0_LDP,W_INEXT))
	              INEXT_V=(1.0_LDP-W_INEXT)*
	1              I_M_DEPTH_NEXT(KC1,LS,FREQ_IDX_ITER+NM_RAY(I))+
	2              W_INEXT*I_M_DEPTH_NEXT(KC2,LS,FREQ_IDX_ITER+NM_RAY(I))
	            ELSE
	              INEXT_V=
	1              I_M_DEPTH_NEXT(KC1,LS,FREQ_IDX_ITER+NM_RAY(I))
	            END IF
	          ELSE IF(KC1 .GE. 1)THEN
	            INEXT_V=I_M_DEPTH_NEXT(KC1,LS,FREQ_IDX_ITER+NM_RAY(I))
	          ELSE
	            INEXT_V=I_M_DEPTH_NEXT(KC2,LS,FREQ_IDX_ITER+NM_RAY(I))
	          END IF
	          SOURCE_PRIME(I)=SOURCE_RAY(I)+ABS(Q(I))*INEXT_V
	        ELSE
!
!	          GAM < 0 but no I_NEXT slot is available.  Drop the
!	          frequency coupling entirely rather than applying it with a
!	          sign inconsistent with CHI_RAY.
!
	          FG_NOSLOT_RAY(LS)=FG_NOSLOT_RAY(LS)+1
	          SOURCE_PRIME(I)=SOURCE_RAY(I)
	        END IF
	      ELSE
!	        GAM < 0 but HAS_COMPRESSING is false: the depth-grid SIGMA_TRUE
!	        never dropped below -1, so I_NEXT was not allocated, yet cubic
!	        interpolation of SIGMA_RAY can still overshoot below -1.
!	        Signed Q*I_PREV would subtract the upwind intensity (Q<0).
!	        Drop frequency coupling, matching the no-slot branch above.
	        FG_NOSLOT_RAY(LS)=FG_NOSLOT_RAY(LS)+1
	        SOURCE_PRIME(I)=SOURCE_RAY(I)
	      END IF
	    END DO
!
! Capture chi_hat and the inward modified source S'_M for the serial
! diagonal-ALO accumulation (these are OMP-PRIVATE inside this loop).
!
	    ALO_CHIHAT_RAY(1:NI,LS)=CHI_RAY(1:NI)
	    ALO_SPR_M_RAY(1:NI,LS)=SOURCE_PRIME(1:NI)
!
! S_tilde computation REMOVED -- see comment above (frequency formulation
! has no 4*a*I term, so no explicit S_tilde source is needed).
! edited by Sean Brennan (07-Sep-2026)
	    DO I=1,NI-1
	      IF(DTAU(I,LS) .GT. 1.0E-30_LDP)THEN
	        S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/DTAU(I,LS)
	      ELSE
	        S(I)=0.0_LDP
	      END IF
	    END DO
!
! Now compute the derivatives node I.
!
! Guard DTAU sums against zero to prevent 0/0 NaN in dS.
! In compression zones, consecutive optical-depth increments can
! both be zero at low-opacity boundary points.
!
! edited by Sean Brennan (05-Sep-2026)
	    T1_DTAU=DTAU(1,LS)+DTAU(2,LS)
	    IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	      dS(1)=0.0_LDP
	    ELSE
	      dS(1)=S(1) +(S(1)-S(2))*DTAU(1,LS)/T1_DTAU
	    END IF
	    DO I=2,NI-1
	      T1_DTAU=DTAU(I-1,LS)+DTAU(I,LS)
	      IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	        dS(I)=0.0_LDP
	      ELSE
	        dS(I)=(S(I-1)*DTAU(I,LS)+S(I)*DTAU(I-1,LS))/T1_DTAU
	      END IF
	    END DO
!
! NI=2 (allowed on the thinnest insert-free rays) has no NI-2 node;
! constant-slope extrapolation is the only defined derivative.
!
	    IF(NI .GE. 3)THEN
	      T1_DTAU=DTAU(NI-2,LS)+DTAU(NI-1,LS)
	      IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	        dS(NI)=0.0_LDP
	      ELSE
	        dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*DTAU(NI-1,LS)/T1_DTAU
	      END IF
	    ELSE
	      dS(NI)=S(NI-1)
	    END IF
!
! Adjust first derivatives so that function is monotonic  in each interval.
!
	    dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*
	1                      MIN(ABS(S(1)),0.5_LDP*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1               MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1               MIN(ABS(S(NI-1)),0.5_LDP*ABS(dS(NI)))
!
            I_M(1,LS)=0.0_LDP
	    DO I=1,NI-1
	      I_M(I+1,LS)=I_M(I,LS)*A0(I,LS)+ (
	1             SOURCE_PRIME(I)*A1(I,LS)
	1        +    SOURCE_PRIME(I+1)*A2(I,LS)
	1        +    dS(I)*A3(I,LS)
	1        +    dS(I+1)*A4(I,LS) )
	    END DO
!
	    DO I=1,NI-1
	      IF(I_M(I+1,LS) .LT. 0 .OR. I_M(I+1,LS) .NE. I_M(I+1,LS))THEN
!
!	        Non-monotonic flows can drive the short-characteristics
!	        intensity slightly negative in compression zones where the
!	        source-function advection term changes sign.  Enforce the
!	        physical bound I >= 0 by clamping to a small positive floor
!	        (1.0E-30) and record the event for diagnostics.  This is the
!	        INT/INS safety net required by the non-monotonic solver; it
!	        does NOT alter the opacity (DTAU is non-negative since
!	        chi_hat = chi + |GAM|/dLOG_NU > 0) and only caps
!	        unphysical overshoot from the explicit source term.
!
!	        The .NE. self-comparison catches NaN (NaN /= NaN is TRUE in
!	        IEEE-754), which .LT. 0 alone misses.  Without this, NaN
!	        intensities propagate into I_M_DEPTH_NEXT and contaminate
!	        the source function of subsequent great iterations, ultimately
!	        corrupting the BA matrix during FULL/BA iterations.
!
!	        Note: per-occurrence WRITE removed - it caused 20+ GB OUTGEN
!	        files.  The N_NEG_INT counter is now incremented with
!	        !$OMP ATOMIC to prevent the OpenMP race condition that
!	        occurred in the parallel ray loop.
!	        The count is still tracked for end-of-run summary diagnostics.
!
!$OMP ATOMIC
! edited by Sean Brennan
	        N_NEG_INT=N_NEG_INT+1
	        I_M(I+1,LS)=1.0E-30_LDP
!
! Phase 2 first-occurrence diagnostic for I_M failure.
! Prints full context of the FIRST negative/NaN I_M value.
!$OMP CRITICAL(FIRST_FAIL_CRIT)
              IF(.NOT. FIRST_FAIL_FOUND)THEN
                FIRST_FAIL_FOUND=.TRUE.
                FIRST_FAIL_TYPE=1
                J=ERROR_LU()
                WRITE(J,'(/,A)')' ==================== PHASE 2: FIRST NUMERICAL FAILURE (I_M) ==================='
                WRITE(J,*)"Type: I_M negative or NaN (inward ray)"
                WRITE(J,"(A,I8)")"  Freq index: ",FREQ_IDX_ITER
                WRITE(J,"(A,ES16.8)")"  Frequency: ",FREQ
                WRITE(J,"(A,I6)")"  Ray (LS): ",LS
                WRITE(J,"(A,I6)")"  Ray point (I+1): ",I+1
                WRITE(J,"(A,ES16.8)")"  R_RAY: ",R_RAY(I+1,LS)
                WRITE(J,"(A,ES16.8)")"  V_RAY: ",V_RAY(I+1)
                WRITE(J,"(A,ES16.8)")"  MU: ",Z(I+1,LS)/R_RAY(I+1,LS)
                WRITE(J,"(A,ES16.8)")"  SIGMA_RAY: ",SIGMA_RAY(I+1)
                WRITE(J,"(A,ES16.8)")"  GAM: ",GAM(I+1,LS)
                WRITE(J,"(A,ES16.8)")"  GAM_MAG: ",GAM_MAG(I+1,LS)
                WRITE(J,"(A,ES16.8)")"  CHI_RAY (chi_hat): ",CHI_RAY(I+1)
                WRITE(J,"(A,ES16.8)")"  dLOG_NU: ",dLOG_NU
                WRITE(J,"(A,ES16.8)")"  Q: ",Q(I+1)
                WRITE(J,"(A,ES16.8)")"  DTAU: ",DTAU(I,LS)
                WRITE(J,"(A,ES16.8)")"  SOURCE_RAY: ",SOURCE_RAY(I+1)
                WRITE(J,"(A,ES16.8)")"  SOURCE_PRIME: ",SOURCE_PRIME(I+1)
                WRITE(J,"(A,ES16.8)")"  I_M (failed): ",I_M(I+1,LS)
                WRITE(J,"(A,ES16.8)")"  I_M(prev pt): ",I_M(I,LS)
                WRITE(J,"(A,ES16.8)")"  A0: ",A0(I,LS)
                WRITE(J,"(A,ES16.8)")"  A1: ",A1(I,LS)
                WRITE(J,"(A,ES16.8)")"  A2: ",A2(I,LS)
                WRITE(J,"(A,ES16.8)")"  A3: ",A3(I,LS)
                WRITE(J,"(A,ES16.8)")"  A4: ",A4(I,LS)
                WRITE(J,"(A,ES16.8)")"  dS(I): ",dS(I)
                WRITE(J,"(A,ES16.8)")"  dS(I+1): ",dS(I+1)
                WRITE(J,"(A,ES16.8)")"  S(I): ",S(I)
                IF(GAM(I+1,LS) .GE. 0.0_LDP)THEN
                  WRITE(J,"(A)")"  Region: EXPANDING (GAM >= 0)"
                  WRITE(J,"(A,ES16.8)")"  I_M_PREV(upwind): ",I_M_PREV(I+1,LS)
                ELSE
                  WRITE(J,"(A)")"  Region: COMPRESSING (GAM < 0)"
                  K=RAY_PNT(I+1,LS)-ND_ADD
                  IF(K .GE. 1 .AND. K .LE. ND)THEN
                    KC=COMP_DEPTH_IDX(K)
                    IF(KC .GE. 1 .AND. FREQ_IDX_ITER+1 .LE. NCF_IN)THEN
                      WRITE(J,"(A,ES16.8)")"  I_M_DEPTH_NEXT: ",I_M_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER+1)
                    ELSE
                      WRITE(J,"(A)")"  I_M_DEPTH_NEXT: NOT AVAILABLE (no slot)"
                    END IF
                  ELSE
                    WRITE(J,"(A)")"  I_M_DEPTH_NEXT: NOT AVAILABLE (K out of range)"
                  END IF
                END IF
                WRITE(J,'(/,A)')' ========================== END PHASE 2 FIRST FAILURE =========================='
              END IF
!$OMP END CRITICAL(FIRST_FAIL_CRIT)
	      END IF
	    END DO
!
! Upper cap on I_M: prevent the explicit frequency coupling from
! amplifying the intensity beyond NONMONO_INEXT_CAP.  In models with
! extreme velocity gradients (|SIGMA| >> 1), the cumulative growth
! factor sum(DTAU*Q) over all depths can be >> 1, making the explicit
! scheme unstable.  The cap is far above any physical intensity and
! only prevents Inf/overflow.  Applied after the negative/NaN clamp
! so that both bounds are enforced.  The I_M_STORE cap below then
! ensures I_M_PREV (expanding-zone feedback) is also bounded.
! edited by Sean Brennan (14-Sep-2026)
	    DO I=1,NI
	      IF(I_M(I,LS) .GT. NONMONO_INEXT_CAP)I_M(I,LS)=NONMONO_INEXT_CAP
	    END DO
!
! ******************* OUTWARD DIRECTED RAYS *********************************
!
! Compute the Source function for outward directed rays, and find the
! monotonic interpolating polynomial.
!
!        Source function for outward directed rays.
!        Full Knop et al. 2009: sign-dependent upwind frequency.
!        Expanding (GAM >= 0): upwind = previous frequency = I_P_PREV
!        Compressing (GAM < 0): upwind = next frequency = I_P_DEPTH_NEXT
!          (from previous great iteration, Gauss-Seidel approach)
!
! edited by Sean Brennan
	    DO I=1,NI
	      IF(GAM(I,LS) .GE. 0.0_LDP)THEN
	        SOURCE_PRIME(I)=SOURCE_RAY(I)+Q(I)*I_P_PREV(I,LS)
	      ELSE IF(HAS_COMPRESSING)THEN
	        K=RAY_PNT(I,LS)-ND_ADD
!
!	        Linear interpolation of the lagged I_NEXT between the two
!	        radial depths bracketing the ray point (audit 17-Sep-2026).
!	        The nearest-depth lookup used previously introduced a
!	        spatial discontinuity exactly where the velocity gradient is
!	        strongest.  When only one bracketing depth is a compressing
!	        (stored) depth, that value is used -- the intensity at the
!	        next frequency is only stored at compressing depths, so a
!	        non-compressing neighbour has no I_NEXT to interpolate with.
!
	        KC1=0
	        KC2=0
	        IF(K .GE. 1 .AND. K .LE. ND)KC1=COMP_DEPTH_IDX(K)
	        IF(K+1 .GE. 1 .AND. K+1 .LE. ND)KC2=COMP_DEPTH_IDX(K+1)
	        IF((KC1 .GE. 1 .OR. KC2 .GE. 1) .AND.
	1       FREQ_IDX_ITER+NM_RAY(I) .GE. 1 .AND.
	2       FREQ_IDX_ITER+NM_RAY(I) .LE. NCF_IN)THEN
!
!	          See the matching comment in the inward-ray block: with the
!	          upwind frequency point switched, the coefficient must be
!	          |GAM|/dLOG_NU/chi_hat, not the signed Q.
!	          edited by Sean Brennan
!
	          IF(KC1 .GE. 1 .AND. KC2 .GE. 1)THEN
	            IF(ABS(R(K+1)-R(K)) .GT. 1.0E-14_LDP*
	1          MAX(ABS(R(K)),ABS(R(K+1)),R(1)))THEN
	              W_INEXT=(R_RAY(I,LS)-R(K))/(R(K+1)-R(K))
	              W_INEXT=MAX(0.0_LDP,MIN(1.0_LDP,W_INEXT))
	              INEXT_V=(1.0_LDP-W_INEXT)*
	1              I_P_DEPTH_NEXT(KC1,LS,FREQ_IDX_ITER+NM_RAY(I))+
	2              W_INEXT*I_P_DEPTH_NEXT(KC2,LS,FREQ_IDX_ITER+NM_RAY(I))
	            ELSE
	              INEXT_V=
	1              I_P_DEPTH_NEXT(KC1,LS,FREQ_IDX_ITER+NM_RAY(I))
	            END IF
	          ELSE IF(KC1 .GE. 1)THEN
	            INEXT_V=I_P_DEPTH_NEXT(KC1,LS,FREQ_IDX_ITER+NM_RAY(I))
	          ELSE
	            INEXT_V=I_P_DEPTH_NEXT(KC2,LS,FREQ_IDX_ITER+NM_RAY(I))
	          END IF
	          SOURCE_PRIME(I)=SOURCE_RAY(I)+ABS(Q(I))*INEXT_V
	        ELSE
!
!	          GAM < 0 but no I_NEXT slot is available.  Drop the
!	          frequency coupling rather than applying an inconsistent sign.
!
	          FG_NOSLOT_RAY(LS)=FG_NOSLOT_RAY(LS)+1
	          SOURCE_PRIME(I)=SOURCE_RAY(I)
	        END IF
	      ELSE
!	        GAM < 0 without I_NEXT storage (cubic SIGMA overshoot).
!	        Drop frequency coupling rather than applying signed Q*I_PREV.
	        FG_NOSLOT_RAY(LS)=FG_NOSLOT_RAY(LS)+1
	        SOURCE_PRIME(I)=SOURCE_RAY(I)
	      END IF
	    END DO
!
! Capture the outward modified source S'_P for the diagonal ALO.
!
	    ALO_SPR_P_RAY(1:NI,LS)=SOURCE_PRIME(1:NI)
!
! S_tilde computation REMOVED -- see comment above (frequency formulation
! has no 4*a*I term, so no explicit S_tilde source is needed).
! edited by Sean Brennan (07-Sep-2026)
	    DO I=1,NI-1
	      IF(DTAU(I,LS) .GT. 1.0E-30_LDP)THEN
	        S(I)=(SOURCE_PRIME(I+1)-SOURCE_PRIME(I))/DTAU(I,LS)
	      ELSE
	        S(I)=0.0_LDP
	      END IF
	    END DO
!
! Now compute the derivatives at node I.
!
! Guard DTAU sums against zero (same as inward sweep).
! edited by Sean Brennan (05-Sep-2026)
	    T1_DTAU=DTAU(1,LS)+DTAU(2,LS)
	    IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	      dS(1)=0.0_LDP
	    ELSE
	      dS(1)=S(1) +(S(1)-S(2))*DTAU(1,LS)/T1_DTAU
	    END IF
	    DO I=2,NI-1
	      T1_DTAU=DTAU(I-1,LS)+DTAU(I,LS)
	      IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	        dS(I)=0.0_LDP
	      ELSE
	        dS(I)=(S(I-1)*DTAU(I,LS)+S(I)*DTAU(I-1,LS))/T1_DTAU
	      END IF
	    END DO
!
! NI=2 (allowed on the thinnest insert-free rays) has no NI-2 node;
! constant-slope extrapolation is the only defined derivative.
!
	    IF(NI .GE. 3)THEN
	      T1_DTAU=DTAU(NI-2,LS)+DTAU(NI-1,LS)
	      IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	        dS(NI)=0.0_LDP
	      ELSE
	        dS(NI)=S(NI-1)+(S(NI-1)-S(NI-2))*DTAU(NI-1,LS)/T1_DTAU
	      END IF
	    ELSE
	      dS(NI)=S(NI-1)
	    END IF
!
! Adjust the first derivatives so that function is monotonic in each interval.
!
	    dS(1)=( SIGN(ONE,S(1))+SIGN(ONE,dS(1)) )*
	1                      MIN(ABS(S(1)),0.5_LDP*ABS(dS(1)))
	    DO I=2,NI-1
	      dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	    END DO
	    dS(NI)=( SIGN(ONE,S(NI-1))+SIGN(ONE,dS(NI)) )*
	1            MIN(ABS(S(NI-1)),0.5_LDP*ABS(dS(NI)))
!
	    IF(INNER_BND_METH .EQ. 'DIFFUSION' .AND. LS .LE. NC)THEN
	      I_CORE=( ETA_RAY(NI)+ DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND) )/CHI_RAY(NI)
	    ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	      I_CORE=I_M(NI,LS)
	    ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	      I_CORE=I_M(NI,LS)
	      IF(LS .LE. NC_PNT_SRCE)I_CORE=IC+I_M(NI,LS)
	    ELSE
	      I_CORE=IC
	    END IF
	    IF(LS .LE. NC)THEN
	      I_P(NI,LS)=I_CORE
	    ELSE
	      I_P(NI,LS)=I_M(NI,LS)
	    END IF
	    DO I=NI-1,1,-1
	      I_P(I,LS)=I_P(I+1,LS)*A0(I,LS)+ (
	1             SOURCE_PRIME(I+1)*A1(I,LS)
	1        +    SOURCE_PRIME(I)*A2(I,LS)
	1        -    dS(I+1)*A3(I,LS)
	1        -    dS(I)*A4(I,LS) )
	    END DO
!
	    DO I=1,NI
	      IF(I_P(I,LS) .LT. 0 .OR. I_P(I,LS) .NE. I_P(I,LS))THEN
!
!	        Enforce I >= 0 for outward-directed rays (see I_M note above).
!	        The .NE. self-comparison catches NaN (NaN /= NaN is TRUE).
!$OMP ATOMIC
! edited by Sean Brennan
	        N_NEG_INT=N_NEG_INT+1
	        I_P(I,LS)=1.0E-30_LDP
!
! Phase 2 first-occurrence diagnostic for I_P failure.
!$OMP CRITICAL(FIRST_FAIL_CRIT)
              IF(.NOT. FIRST_FAIL_FOUND)THEN
                FIRST_FAIL_FOUND=.TRUE.
                FIRST_FAIL_TYPE=2
                J=ERROR_LU()
                WRITE(J,'(/,A)')' ==================== PHASE 2: FIRST NUMERICAL FAILURE (I_P) ==================='
                WRITE(J,*)"Type: I_P negative or NaN (outward ray)"
                WRITE(J,"(A,I8)")"  Freq index: ",FREQ_IDX_ITER
                WRITE(J,"(A,ES16.8)")"  Frequency: ",FREQ
                WRITE(J,"(A,I6)")"  Ray (LS): ",LS
                WRITE(J,"(A,I6)")"  Ray point (I): ",I
                WRITE(J,"(A,ES16.8)")"  R_RAY: ",R_RAY(I,LS)
                WRITE(J,"(A,ES16.8)")"  V_RAY: ",V_RAY(I)
                WRITE(J,"(A,ES16.8)")"  MU: ",Z(I,LS)/R_RAY(I,LS)
                WRITE(J,"(A,ES16.8)")"  SIGMA_RAY: ",SIGMA_RAY(I)
                WRITE(J,"(A,ES16.8)")"  GAM: ",GAM(I,LS)
                WRITE(J,"(A,ES16.8)")"  GAM_MAG: ",GAM_MAG(I,LS)
                WRITE(J,"(A,ES16.8)")"  CHI_RAY (chi_hat): ",CHI_RAY(I)
                WRITE(J,"(A,ES16.8)")"  dLOG_NU: ",dLOG_NU
                WRITE(J,"(A,ES16.8)")"  Q: ",Q(I)
                WRITE(J,"(A,ES16.8)")"  DTAU: ",DTAU(I,LS)
                WRITE(J,"(A,ES16.8)")"  SOURCE_RAY: ",SOURCE_RAY(I)
                WRITE(J,"(A,ES16.8)")"  SOURCE_PRIME: ",SOURCE_PRIME(I)
                WRITE(J,"(A,ES16.8)")"  I_P (failed): ",I_P(I,LS)
                WRITE(J,"(A,ES16.8)")"  I_P(next pt): ",I_P(I+1,LS)
                WRITE(J,"(A,ES16.8)")"  A0: ",A0(I,LS)
                WRITE(J,"(A,ES16.8)")"  A1: ",A1(I,LS)
                WRITE(J,"(A,ES16.8)")"  A2: ",A2(I,LS)
                WRITE(J,"(A,ES16.8)")"  A3: ",A3(I,LS)
                WRITE(J,"(A,ES16.8)")"  A4: ",A4(I,LS)
                WRITE(J,"(A,ES16.8)")"  dS(I): ",dS(I)
                WRITE(J,"(A,ES16.8)")"  dS(I+1): ",dS(I+1)
                WRITE(J,"(A,ES16.8)")"  S(I): ",S(I)
                IF(GAM(I,LS) .GE. 0.0_LDP)THEN
                  WRITE(J,"(A)")"  Region: EXPANDING (GAM >= 0)"
                  WRITE(J,"(A,ES16.8)")"  I_P_PREV(upwind): ",I_P_PREV(I,LS)
                ELSE
                  WRITE(J,"(A)")"  Region: COMPRESSING (GAM < 0)"
                  K=RAY_PNT(I,LS)-ND_ADD
                  IF(K .GE. 1 .AND. K .LE. ND)THEN
                    KC=COMP_DEPTH_IDX(K)
                    IF(KC .GE. 1 .AND. FREQ_IDX_ITER+1 .LE. NCF_IN)THEN
                      WRITE(J,"(A,ES16.8)")"  I_P_DEPTH_NEXT: ",I_P_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER+1)
                    ELSE
                      WRITE(J,"(A)")"  I_P_DEPTH_NEXT: NOT AVAILABLE (no slot)"
                    END IF
                  ELSE
                    WRITE(J,"(A)")"  I_P_DEPTH_NEXT: NOT AVAILABLE (K out of range)"
                  END IF
                END IF
                WRITE(J,'(/,A)')' ========================== END PHASE 2 FIRST FAILURE =========================='
              END IF
!$OMP END CRITICAL(FIRST_FAIL_CRIT)
	      END IF
	    END DO
!
! Upper cap on I_P: same rationale as the I_M cap above.  Prevents
! the explicit frequency coupling from amplifying the outward-directed
! intensity beyond NONMONO_INEXT_CAP in models with extreme velocity
! gradients.  Applied after the negative/NaN clamp.
! edited by Sean Brennan (14-Sep-2026)
	    DO I=1,NI
	      IF(I_P(I,LS) .GT. NONMONO_INEXT_CAP)I_P(I,LS)=NONMONO_INEXT_CAP
	    END DO
C Note that V=AV(1)-IBOUND.
C
1000	    CONTINUE
	    K=J_PNT(1,LS)
	    CV_BOUND(LS)=0.5_LDP*(I_P(K,LS)-I_M(K,LS))
	    I_M_IN_BND(LS)=I_M(NI,LS)
!
!	    CALL TUNE(2,'FG_INT_INT')
	  END DO			!LS
!$OMP END PARALLEL DO
!
! Compute the mean intensity like variable U at each grid point. Used to
! compute J.
!
!$OMP PARALLEL DO PRIVATE(I)
!
	  DO LS=1,NP
	    DO I=1,NI_RAY(LS)
	      AV(I,LS)=0.5_LDP*( I_P(I,LS)+I_M(I,LS) )
	    END DO
	  END DO
!
!$OMP END PARALLEL DO
!
! We define CV on the gregular grid. We perform the interpolation onto
! the original R grid as we compute H.
!
!$OMP PARALLEL DO PRIVATE(I)
	  DO LS=1,NP
	    DO I=1,NI_RAY(LS)
	      CV(I,LS)= 0.5_LDP*( I_P(I,LS)-I_M(I,LS) )
	    END DO
	  END DO
!$OMP END PARALLEL DO
!
	  I_P_STORE(:,:)=MIN(I_P(:,:),NONMONO_INEXT_CAP)
	  I_M_STORE(:,:)=MIN(I_M(:,:),NONMONO_INEXT_CAP)
!
! Fold this frequency's compressing-cell sink count into the
! great-iteration census (serial region; the ray loop is done).
!
	  FG_NOSLOT_CNT=FG_NOSLOT_CNT+
	1          SUM(INT(FG_NOSLOT_RAY(1:NP),8))
!
! Cap the expanding-zone frequency feedback (I_M_PREV / I_P_PREV) at
! NONMONO_INEXT_CAP, matching the cap already applied to the compressing-
! zone feedback (I_M_DEPTH_NEXT / I_P_DEPTH_NEXT) below.  Without this,
! the explicit frequency coupling SOURCE_PRIME = SOURCE_RAY + Q*I_M_PREV
! (expanding, GAM >= 0) amplifies the intensity exponentially across
! frequencies when the cumulative growth factor sum(DTAU*Q) >> 1, which
! occurs in models with extreme velocity gradients (|SIGMA| >> 1).
! The compressing-zone feedback was already capped; the expanding-zone
! feedback was not, creating an asymmetry that allowed runaway growth.
! The cap is far above any physical intensity (B_nu ~ 1e-6 to 1e3) and
! only prevents Inf/overflow from destabilising the SE/BA solvers.
! edited by Sean Brennan (14-Sep-2026)
!
! Store intensities at depth grid points for Knop GS iteration.
! Optimized (A+B): only store for compressing depths (COMP_DEPTH_IDX > 0),
! and write directly to NEXT (in-place, no CURR copy needed).
! J_PNT(K,LS) maps physical depth K to the nearest ray point index I.
! Only physical depths K=1..ND are stored (not the outer thick-boundary
! extension points, which are monotonic and use I_PREV upwinding).
!
! edited by Sean Brennan
! Skip the in-place I_NEXT writes when FREEZE_I_NEXT is set.  The flag
! exists so a JFNK residual evaluation (JFNK_SOLVE_CORR_V1 sets it
! around the coupled eval sweep) sees the same Gauss-Seidel I_NEXT
! state as the baseline; the default path leaves it .FALSE.
!
	  IF(HAS_COMPRESSING .AND. FREQ_IDX_ITER .GE. 1 .AND.
	1   FREQ_IDX_ITER .LE. NCF_IN .AND.
	1   .NOT. FREEZE_I_NEXT)THEN
! FG_GS_MAXCHG/FG_INEXT_CAPPED accumulate the last-sweep convergence
! metrics written into the GS checkpoint; a restored checkpoint that is
! complete but still changing by >50% (or that had cap/floor events)
! does not skip the restart warm-up (see SOLVE_FOR_POPS).
! This block is serial (after the ray PARALLEL DOs), so the scalar
! accumulators cannot race.
!
! Column peak per ray for the >50% churn count below.  Without this
! floor a cell whose stored value sits at the 1E-30 write floor counts
! as churned whenever the fresh intensity differs at all (|T1-T2|/T2
! exceeds 0.5 for any T1 > 1.5E-30), so FG_GS_NCONV never decays and
! the KNOP maturity gate in SOLVE_FOR_POPS stays shut even when every
! significant cell is converged.  The 1E-4*peak test mirrors the one
! the refresh residual applies in MARCH (VMAX_L).
! edited by Devin (28-Sep-2026)
!
	    DO LS=1,NP
	      VMAX_IP(LS)=MAX(1.0E-25_LDP,
	1        MAXVAL(I_P(1:NI_RAY(LS),LS),
	2        MASK=I_P(1:NI_RAY(LS),LS) .EQ. I_P(1:NI_RAY(LS),LS)))
	      VMAX_IM(LS)=MAX(1.0E-25_LDP,
	1        MAXVAL(I_M(1:NI_RAY(LS),LS),
	2        MASK=I_M(1:NI_RAY(LS),LS) .EQ. I_M(1:NI_RAY(LS),LS)))
	    END DO
	    KN_DUMP_NS=0
	    DO K=1,ND
	      KC=COMP_DEPTH_IDX(K)
	      IF(KC .GE. 1)THEN
	        DO LS=1,NP
	          I=J_PNT(K,LS)
	          IF(I .GE. 1 .AND. I .LE. NI_RAY(LS))THEN
	            FG_SW_CAPTOT=FG_SW_CAPTOT+1
	            IF(I_P(I,LS) .EQ. I_P(I,LS) .AND.
	1              I_P(I,LS) .GT. 0.0_LDP)THEN
	              T1=MIN(I_P(I,LS),NONMONO_INEXT_CAP)
	              IF(I_P(I,LS) .GT. NONMONO_INEXT_CAP)THEN
	                FG_INEXT_CAPPED=FG_INEXT_CAPPED+1
	                FG_SW_CAPPED=FG_SW_CAPPED+1
	              END IF
	            ELSE
	              T1=1.0E-30_LDP
	              IF(.NOT. (I_P(I,LS) .EQ. I_P(I,LS)))THEN
	                FG_INEXT_CAPPED=FG_INEXT_CAPPED+1
	                FG_SW_CAPPED=FG_SW_CAPPED+1
	              END IF
	            END IF
	            T2=I_P_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER)
	            IF(T2 .GT. 0.0_LDP)THEN
!
! MAXCHG denominators are floored at the ray peak (same contract as
! MOM_GS_MAXCHG): a stored cell ~1E-300 vs a fresh 1E-38 reports
! relchg ~1E262 on a cell the NCONV gate already deems
! insignificant, and GS_MAX_ALL=MAX(...) then pins the maturity
! gate at GS_UNSAFE_TOL forever.  Job 12255111 ran with
! sweep drift=1.5E263 from exactly this.
!
	              IF(ABS(T1-T2)/T2 .GT. 0.5_LDP .AND.
	1                MAX(T1,T2) .GT. 1.0E-4_LDP*VMAX_IP(LS))
	2                FG_GS_NCONV=FG_GS_NCONV+1
	              FG_GS_MAXCHG=MAX(FG_GS_MAXCHG,ABS(T1-T2)/
	1                MAX(T2,1.0E-4_LDP*VMAX_IP(LS)))
!
! KNSTOR: fixed-point diagnostic for the lagged I_NEXT store.  A
! cell whose freshly solved intensity still differs from the stored
! value by >25% after relaxation cycles is part of the residual
! floor that holds the field above GS_UNSAFE_TOL.
!
	              IF(PRINT_KNOP_DEBUG .AND. ABS(T1-T2)/T2 .GT.
	1               0.25_LDP .AND.
	2               MAX(T1,T2) .GT. 1.0E-4_LDP*VMAX_IP(LS) .AND.
	3               KN_DUMP_NS .LT. 2000)THEN
	                KN_DUMP_NS=KN_DUMP_NS+1
	                WRITE(LUER,'(A,I6,A,I4,A,I4,A,I5,4(A,ES11.4))')
	1               ' KNSTOR P f=',FREQ_IDX_ITER,' r=',LS,' k=',K,
	2               ' i=',I,' st=',T2,' sv=',T1,' rel=',
	3               ABS(T1-T2)/T2,' sig=',SIGMA_TRUE_IN(K)
	              END IF
	            ELSE IF(T1 .GT. 1.0E-4_LDP*VMAX_IP(LS))THEN
	              FG_GS_NCONV=FG_GS_NCONV+1
	              FG_GS_MAXCHG=MAX(FG_GS_MAXCHG,1.0_LDP)
	            END IF
	            IF(T2 .GT. 0.0_LDP)THEN
	              I_P_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER)=
	1              T2+INEXT_RELAX*(T1-T2)
	            ELSE
	              I_P_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER)=T1
	            END IF
	            FG_SW_CAPTOT=FG_SW_CAPTOT+1
	            IF(I_M(I,LS) .EQ. I_M(I,LS) .AND.
	1              I_M(I,LS) .GT. 0.0_LDP)THEN
	              T1=MIN(I_M(I,LS),NONMONO_INEXT_CAP)
	              IF(I_M(I,LS) .GT. NONMONO_INEXT_CAP)THEN
	                FG_INEXT_CAPPED=FG_INEXT_CAPPED+1
	                FG_SW_CAPPED=FG_SW_CAPPED+1
	              END IF
	            ELSE
	              T1=1.0E-30_LDP
	              IF(.NOT. (I_M(I,LS) .EQ. I_M(I,LS)))THEN
	                FG_INEXT_CAPPED=FG_INEXT_CAPPED+1
	                FG_SW_CAPPED=FG_SW_CAPPED+1
	              END IF
	            END IF
	            T2=I_M_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER)
	            IF(T2 .GT. 0.0_LDP)THEN
	              IF(ABS(T1-T2)/T2 .GT. 0.5_LDP .AND.
	1                MAX(T1,T2) .GT. 1.0E-4_LDP*VMAX_IM(LS))
	2                FG_GS_NCONV=FG_GS_NCONV+1
	              FG_GS_MAXCHG=MAX(FG_GS_MAXCHG,ABS(T1-T2)/
	1                MAX(T2,1.0E-4_LDP*VMAX_IM(LS)))
	              IF(PRINT_KNOP_DEBUG .AND. ABS(T1-T2)/T2 .GT.
	1               0.25_LDP .AND.
	2               MAX(T1,T2) .GT. 1.0E-4_LDP*VMAX_IM(LS) .AND.
	3               KN_DUMP_NS .LT. 2000)THEN
	                KN_DUMP_NS=KN_DUMP_NS+1
	                WRITE(LUER,'(A,I6,A,I4,A,I4,A,I5,4(A,ES11.4))')
	1               ' KNSTOR M f=',FREQ_IDX_ITER,' r=',LS,' k=',K,
	2               ' i=',I,' st=',T2,' sv=',T1,' rel=',
	3               ABS(T1-T2)/T2,' sig=',SIGMA_TRUE_IN(K)
	              END IF
	            ELSE IF(T1 .GT. 1.0E-4_LDP*VMAX_IM(LS))THEN
	              FG_GS_NCONV=FG_GS_NCONV+1
	              FG_GS_MAXCHG=MAX(FG_GS_MAXCHG,1.0_LDP)
	            END IF
	            IF(T2 .GT. 0.0_LDP)THEN
	              I_M_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER)=
	1              T2+INEXT_RELAX*(T1-T2)
	            ELSE
	              I_M_DEPTH_NEXT(KC,LS,FREQ_IDX_ITER)=T1
	            END IF
	          END IF
	        END DO
	      END IF
	    END DO
	  END IF
!
! Code-audit diagnostic (15-Sep-2026): check the I_NEXT depth mapping
! once per grid definition, on the first frequency of a sweep.  The
! mapping itself is verified unconditionally -- a non-finite ray radius
! or coupling coefficient at a compressing-slot point reports an ERROR.
! No radius-vs-depth equality test: most ray points sit between depth
! nodes and the read side interpolates between the bracketing depths
! (KC1/KC2), so R_RAY!=R(K2) is normal.  The full table (O(NP*NI)
! lines) only prints when PRINT_KNOP_DEBUG=T; on fine ray grids it
! floods OUTGEN otherwise.
!
	  IF(HAS_COMPRESSING .AND. FREQ_IDX_ITER .EQ. 1 .AND.
	1   .NOT. INEXT_MAP_PRINTED)THEN
	    INEXT_MAP_PRINTED=.TRUE.
	    J=ERROR_LU()
	    N_MAP_BAD=0
	    IF(PRINT_KNOP_DEBUG)THEN
	      WRITE(J,'(/,A)')
	1     ' ============ I_NEXT depth-mapping diagnostic (first freq of sweep) ============'
	      WRITE(J,'(A)')
	1     'Read-side map (ray point -> sparse I_NEXT slot):'
	    END IF
	    DO LS=1,NP
	      IF(PRINT_KNOP_DEBUG)THEN
	        WRITE(J,'(A,I6,A,I6)')'  Ray LS=',LS,'  NI_RAY=',NI_RAY(LS)
	        WRITE(J,'(A)')
	1       '      I  RAY_PNT      K     K2     KC'//
	1       '       R_RAY(I)       R(K2)          GAM'
	      END IF
	      DO I=1,NI_RAY(LS)
	        K=RAY_PNT(I,LS)-ND_ADD
	        K2=K
	        IF(K .GE. 1 .AND. K .LT. ND)THEN
	          IF(ABS(R_RAY(I,LS)-R(K+1)) .LT.
	1             ABS(R_RAY(I,LS)-R(K)))K2=K+1
	        END IF
	        KC=0
	        IF(K2 .GE. 1 .AND. K2 .LE. ND)KC=COMP_DEPTH_IDX(K2)
	        IF(KC .GE. 1)THEN
	          IF(.NOT.(ABS(R_RAY(I,LS)) .LE. HUGE(1.0_LDP)) .OR.
	1           .NOT.(ABS(GAM(I,LS)) .LE. HUGE(1.0_LDP)))
	2             N_MAP_BAD=N_MAP_BAD+1
	          IF(PRINT_KNOP_DEBUG)WRITE(J,'(5I7,3ES15.6)')
	1              I,RAY_PNT(I,LS),K,K2,KC,R_RAY(I,LS),R(K2),GAM(I,LS)
	        END IF
	      END DO
	    END DO
	    IF(PRINT_KNOP_DEBUG)THEN
	      WRITE(J,'(A)')
	1     'Write-side map (depth K -> J_PNT ray point), LS=1 and LS=NP:'
	      WRITE(J,'(A)')
	1     '      K     KC     LS    J_PNT'//
	1     '       R(K)           R_RAY(J_PNT)'
	      DO K=1,ND
	        KC=COMP_DEPTH_IDX(K)
	        IF(KC .GE. 1)THEN
	          DO LS=1,NP,MAX(1,NP-1)
	            I=J_PNT(K,LS)
	            IF(I .GE. 1 .AND. I .LE. NI_RAY(LS))THEN
	              WRITE(J,'(4I7,2ES15.6)')K,KC,LS,I,R(K),R_RAY(I,LS)
	            END IF
	          END DO
	        END IF
	      END DO
	      WRITE(J,'(/,A)')' ===================== END I_NEXT depth-mapping diagnostic ====================='
	    END IF
	    IF(N_MAP_BAD .GT. 0)THEN
	      WRITE(J,'(A,I8,A)')
	1     'ERROR: I_NEXT depth-mapping found',N_MAP_BAD,
	2     ' anomalous ray point(s): NaN/Inf in R_RAY or GAM'
	      FLUSH(J)
	    END IF
	  END IF
!
! edited by Sean Brennan
!	  IF(NEW_FREQ)THEN
!	     WRITE(117,*)FREQ
!	     DO LS=1,NP
!	       WRITE(117,'(5ES15.6)')P(LS),I_M_PREV(1,LS),I_P_PREV(1,LS), I_M_PREV(2,LS),I_P_PREV(2,LS)
!	     END DO
!	  END IF
!
	END IF				!Solution method
	CALL TUNE(2,'FG_SOL_INT')
!
!
!***************************************************************************
!***************************************************************************
!
!
! Zero boundary conditions.
!
	HBC=0.0_LDP			!H/J at model outer boundary.
	NBC=0.0_LDP			!N/J at model outer boundary.
	IN_HBC=0.0_LDP
!
! Initialize intensity matrices.
!
	JNU_STORE(:)=0.0_LDP			!1:ND
	HNU_STORE(:)=0.0_LDP
	KNU_STORE(:)=0.0_LDP
	NNU_STORE(:)=0.0_LDP
!
! J and K are always evaluated on the nodes.
!
	CALL TUNE(1,'J_LOOP')
	DO LS=1,NP
	  DO I=1,MIN(ND,ND-(LS-NC-1))
	    K=J_PNT(I,LS)
	    JNU_STORE(I)=JNU_STORE(I)+JQW(I,LS)*AV(K,LS)
	    KNU_STORE(I)=KNU_STORE(I)+KQW(I,LS)*AV(K,LS)
	  END DO
	END DO
!
! Diagonal Olson-Auer-Buchler approximate Lambda operator from the
! actual formal solution.  Loop bounds, weights JQW(I,LS) and the
! J_PNT depth->ray mapping are IDENTICAL to the JNU_STORE accumulation
! above, so ALO_DJD* are exactly dJ/dETA and dJ/dCHI diagonal
! sensitivities of the computed JNU_STORE.
!
! For ray node K on ray LS:
!   dI_M(K)/dS'(K) = A2(K-1,LS)  (K>=2; 0 at the outer boundary K=1)
!   dI_P(K)/dS'(K) = A2(K,LS)    (K<=NI-1)
!   at K=NI: non-core rays (LS>NC) and ZERO_FLUX/PNT_SRCE boundaries
!     have I_P(NI)=I_M(NI), so dI_P/dS' = dI_M/dS' = A2(NI-1,LS);
!     core rays with DIFFUSION/RAY inner boundary have I_P(NI)=I_CORE
!     (independent of S'), so dI_P/dS' = 0.
! AV(K,LS)=0.5*(I_P+I_M) => dAV/dS' = 0.5*(dI_P/dS' + dI_M/dS').
! With S'_dir = ETA/chi_hat + Q_dir*I_upwind and both terms proportional
! to 1/chi_hat: dS'_dir/dETA = 1/chi_hat and dS'_dir/dCHI = -S'_dir/chi_hat
! (dir = M,P; the inward and outward modified sources differ because
! their upwind intensities I_M_PREV / I_P_PREV differ).
!
! chi_hat, S'_M and S'_P are OMP-PRIVATE inside the ray loop; the
! per-ray captures ALO_CHIHAT_RAY / ALO_SPR_M_RAY / ALO_SPR_P_RAY are
! written there (each thread owns its LS column) and read here.
!
! This block runs in the serial J_LOOP, so no OpenMP race exists.
! Only accumulated when the INTEGRAL sweep actually ran (AV valid);
! SOLUTION_METHOD guards are inherited via ALO_VALID below.
!
! edited by Sean Brennan
!
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL' .AND.
	1   ALLOCATED(ALO_CHIHAT_RAY))THEN
	  DO LS=1,NP
	    NIR_ALO=NI_RAY(LS)
	    DO I=1,MIN(ND,ND-(LS-NC-1))
	      K=J_PNT(I,LS)
	      IF(K .GE. 1 .AND. K .LE. NIR_ALO)THEN
	        ALO_DIM=0.0_LDP
	        IF(K .GE. 2)ALO_DIM=A2(K-1,LS)
	        ALO_DIP=0.0_LDP
	        IF(K .LE. NIR_ALO-1)THEN
	          ALO_DIP=A2(K,LS)
	        ELSE IF(K .EQ. NIR_ALO .AND.
	1           (LS .GT. NC .OR.
	1            INNER_BND_METH .EQ. 'ZERO_FLUX' .OR.
	1            INNER_BND_METH .EQ. 'PNT_SRCE'))THEN
	          ALO_DIP=ALO_DIM
	        END IF
	        ALO_DAV=0.5_LDP*(ALO_DIP+ALO_DIM)
	        LAMSTAR_DIAG(I)=LAMSTAR_DIAG(I)+JQW(I,LS)*ALO_DAV
	        IF(ALO_CHIHAT_RAY(K,LS) .GT. 0.0_LDP)THEN
	          ALO_DJDETA(I)=ALO_DJDETA(I)+
	1             JQW(I,LS)*ALO_DAV/ALO_CHIHAT_RAY(K,LS)
	          ALO_DJDCHI(I)=ALO_DJDCHI(I)-
	1             JQW(I,LS)*0.5_LDP*
	1             (ALO_DIP*ALO_SPR_P_RAY(K,LS)+
	1              ALO_DIM*ALO_SPR_M_RAY(K,LS))/ALO_CHIHAT_RAY(K,LS)
	        END IF
	      END IF
	    END DO
	  END DO
	  ALO_VALID=.TRUE.
!
	  IF(INIT)THEN
	    LAMSTAR_MIN=LAMSTAR_DIAG(1)
	    LAMSTAR_MAX=LAMSTAR_DIAG(1)
	    N_OUT_LAMSTAR=0
	    DO I=1,ND
	      LAMSTAR_MIN=MIN(LAMSTAR_MIN,LAMSTAR_DIAG(I))
	      LAMSTAR_MAX=MAX(LAMSTAR_MAX,LAMSTAR_DIAG(I))
	      IF(LAMSTAR_DIAG(I) .LE. 0.0_LDP .OR.
	1         LAMSTAR_DIAG(I) .GT. 1.0_LDP)N_OUT_LAMSTAR=N_OUT_LAMSTAR+1
	    END DO
	    IF(N_OUT_LAMSTAR .GT. 0)THEN
	      WRITE(LUER,'(A,ES12.4,A,ES12.4,A,I8)')
	1     '  LAMSTAR_DIAG: min=',LAMSTAR_MIN,' max=',LAMSTAR_MAX,
	1     ' outside (0,1]:',N_OUT_LAMSTAR
	    END IF
	  END IF
	END IF
!
! This procedure works for the INTEGRAL or DIFFERENCE approach.
!
	DO LS=1,NP
	  K=J_PNT(1,LS)
!
! T1/T2 are the outward/inward boundary intensities I+ = AV + CV_BOUND
! and I- = AV - CV_BOUND.  A negative value is unphysical -- the
! formal solution produced a moment decomposition that cannot be a
! real intensity pair -- so the floor is now counted in
! FG_N_NEG_BND_I rather than silently absorbed into JMIN_OB.
! .NOT.(.GE.0) also traps NaN/Inf from a corrupted AV or CV_BOUND.
!
	  T1=AV(K,LS)+CV_BOUND(LS)
	  IF(.NOT.(T1 .GE. 0.0_LDP))THEN
	    FG_N_NEG_BND_I=FG_N_NEG_BND_I+1
	    T1=0.0_LDP
	  END IF
	  T2=AV(K,LS)-CV_BOUND(LS)
	  IF(.NOT.(T2 .GE. 0.0_LDP))THEN
	    FG_N_NEG_BND_I=FG_N_NEG_BND_I+1
	    T2=0.0_LDP
	  END IF
	  JPLUS_OB=JPLUS_OB+JQW(1,LS)*T1
	  HPLUS_OB=HPLUS_OB+HQW(1,LS)*T1
	  KPLUS_OB=KPLUS_OB+KQW(1,LS)*T1
	  NPLUS_OB=NPLUS_OB+NQW(1,LS)*T1
	  JMIN_OB=JMIN_OB+JQW(1,LS)*T2
	  HMIN_OB=HMIN_OB+HQW(1,LS)*T2
	  KMIN_OB=KMIN_OB+KQW(1,LS)*T2
	  NMIN_OB=NMIN_OB+NQW(1,LS)*T2
	END DO
!
! Boundary identity audit: since AV=(I_P+I_M)/2 is accumulated with
! the same JQW weights, JNU_STORE(1) must equal 0.5*(JPLUS+JMIN) and
! KNU_STORE(1) likewise.  Any residual means a floored intensity or
! a corrupted CV_BOUND entered the sums -- counted as a boundary
! failure in FG_N_BND_IDENT.  H and N are not checked this way: on
! the midpoint convention HNU_AT_OB uses H_PNT, not J_PNT, so the
! residual there measures definition, not corruption.
!
	T1=ABS(JNU_STORE(1)-0.5_LDP*(JPLUS_OB+JMIN_OB))/
	1   MAX(JNU_STORE(1),0.5_LDP*(JPLUS_OB+JMIN_OB),1.0E-300_LDP)
	T2=ABS(KNU_STORE(1)-0.5_LDP*(KPLUS_OB+KMIN_OB))/
	1   MAX(ABS(KNU_STORE(1)),0.5_LDP*(KPLUS_OB+KMIN_OB),
	2        1.0E-300_LDP)
	FG_MAX_BND_RES=MAX(FG_MAX_BND_RES,T1,T2)
	IF(T1 .GT. 1.0E-03_LDP .OR. T2 .GT. 1.0E-03_LDP)
	1   FG_N_BND_IDENT=FG_N_BND_IDENT+1
!
	DO LS=1,NC+1
	  K=J_PNT(ND,LS)
	  T1=AV(K,LS)-0.5_LDP*I_M_IN_BND(LS)
	  T2=I_M_IN_BND(LS)
	  JPLUS_IB = JPLUS_IB+JQW(ND,LS)*T1
	  HPLUS_IB = HPLUS_IB+HQW(ND,LS)*T1
	  KPLUS_IB = KPLUS_IB+KQW(ND,LS)*T1
	  NPLUS_IB = NPLUS_IB+NQW(ND,LS)*T1
	   JMIN_IB = JMIN_IB +JQW(ND,LS)*T2
	   HMIN_IB = HMIN_IB +HQW(ND,LS)*T2
	   KMIN_IB = KMIN_IB +KQW(ND,LS)*T2
	   NMIN_IB = NMIN_IB +NQW(ND,LS)*T2
	END DO
!
! Evaluate the H & N moments. These may be evaluated at the nodes,
! or at the node midpoints, depending on the solution technique.
!
! With the integral technique, H & N are evaluated at the nodes.
!
	IF(SOLUTION_METHOD .EQ. 'INTEGRAL')THEN
          IF(DEFINE_AT_MID_POINTS)THEN
	    DO LS=1,NP
	      DO I=1,MIN(ND,ND-(LS-NC-1))-1
                K=H_PNT(I,LS)
                T2=0.5_LDP*(R(I)+R(I+1))
                T1=(T2-R_RAY(K,LS))/(R_RAY(K+1,LS)-R_RAY(K,LS))
                HNU_STORE(I)=HNU_STORE(I)+HMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
                NNU_STORE(I)=NNU_STORE(I)+NMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
              END DO
            END DO
	    HNU_AT_OB=0.0_LDP; NNU_AT_OB=0.0_LDP
	    DO LS=1,NP
	      K=H_PNT(1,LS)
	      HNU_AT_OB=HNU_AT_OB+HQW(1,LS)*CV(K,LS)
	      NNU_AT_OB=NNU_AT_OB+NQW(1,LS)*CV(K,LS)
	    END DO
	    HN_DEF_ON_NODES=.FALSE.
	  ELSE
	    HN_DEF_ON_NODES=.TRUE.
	    DO LS=1,NP
	      DO I=1,MIN(ND,ND-(LS-NC-1))
	        K=J_PNT(I,LS)
	        HNU_STORE(I)=HNU_STORE(I)+HQW(I,LS)*CV(K,LS)
	        NNU_STORE(I)=NNU_STORE(I)+NQW(I,LS)*CV(K,LS)
	      END DO
	      IPLUS_P(LS)=2.0_LDP*CV_BOUND(LS)
	    END DO
	    HNU_AT_OB=HNU_STORE(1); NNU_AT_OB=NNU_STORE(1)
	    HNU_AT_IB=HNU_STORE(ND); NNU_AT_IB=NNU_STORE(ND)
	  END IF
	ELSE
!
! With the DIFFERENCE technique, H & N are evaluated at the midpoints
! of the nodes. We only need to interpolate when we have included
! additional points along each ray.
!
	  HN_DEF_ON_NODES=.FALSE.
	  DO LS=1,NP
	    DO I=1,MIN(ND,ND-(LS-NC-1))-1
	      K=H_PNT(I,LS)
	      T2=0.5_LDP*(R(I)+R(I+1))
	      T1=( 2.0_LDP*T2-(R_RAY(K,LS)+R_RAY(K+1,LS)) )/(R_RAY(K+2,LS)-R_RAY(K,LS))
	      HNU_STORE(I)=HNU_STORE(I)+HMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
	      NNU_STORE(I)=NNU_STORE(I)+NMIDQW(I,LS)*((1.0_LDP-T1)*CV(K,LS)+T1*CV(K+1,LS))
	    END DO
	  END DO			!End do LS
!
	  HNU_AT_OB=0.0_LDP; NNU_AT_OB=0.0_LDP
	  DO LS=1,NP
	    HNU_AT_OB=HNU_AT_OB+CV_BOUND(LS)*HQW(1,LS)
	    NNU_AT_OB=NNU_AT_OB+CV_BOUND(LS)*NQW(1,LS)
	    IPLUS_P(LS)=2.0_LDP*CV_BOUND(LS)
	  END DO
	  IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	    HNU_AT_IB=DBB/R(ND)/CHI_RAY(NI)/3.0_LDP
	    NNU_AT_IB=DBB/R(ND)/CHI_RAY(NI)/5.0_LDP
	  ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	    HNU_AT_IB=0.0_LDP; NNU_AT_IB=0.0_LDP
	    DO LS=1,NC_PNT_SRCE
	      HNU_AT_IB=HNU_AT_IB+0.5_LDP*IC*HQW(ND,LS)
	      NNU_AT_IB=NNU_AT_IB+0.5_LDP*IC*NQW(ND,LS)
	    END DO
	  ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	    HNU_AT_IB=0.0_LDP
	    NNU_AT_IB=0.0_LDP
	  ELSE
	    HNU_AT_IB=0.0_LDP; NNU_AT_IB=0.0_LDP
	    DO LS=1,NC+1
	      HNU_AT_IB=HNU_AT_IB+0.5_LDP*(IC-I_M_IN_BND(LS))*HQW(ND,LS)
	      NNU_AT_IB=NNU_AT_IB+0.5_LDP*(IC-I_M_IN_BND(LS))*NQW(ND,LS)
	    END DO
	  END IF
	END IF
!
! Get boundary conditons for moment calculations.
!
	IN_HBC=0.0_LDP
	DO LS=1,NC+1
	   IN_HBC=IN_HBC + HQW(ND,LS)*I_M_IN_BND(LS)
	END DO
!
	CALL TUNE(2,'J_LOOP')
!
!
!
! Compute the boundary Eddington factors.
!
! Clamp only NaN and NEGATIVE JNU_STORE at all depths BEFORE computing
! FEDD or boundary conditions.  This prevents NaN from propagating into
! FEDD (KNU/JNU) and thence into the BA matrix via VAR_MOM_J_CMF_V12.
! The .NE. self-comparison catches NaN (NaN /= NaN is TRUE in IEEE-754).
!
! ROOT-CAUSE FIX (2026-09-10): The previous condition JNU_STORE < 1.0E-30
! caught legitimately small EUV intensities (e.g. at 120 Angstrom the
! Wien tail is far below 1.0E-30 for T ~ 5000-10000 K).  This produced
! ~600k false "NaN/negative JNU" clamps per sweep, all at the first
! (highest) CMF frequency, and the count INCREASED across great iterations
! because the clamped value (1.0E-30) was orders of magnitude ABOVE the
! physical intensity, feeding back an unphysical EUV boost.  V13 only
! checks JNU_STORE < 0 (negative), not small positive.  This change
! matches V13: only NaN and negative values are clamped; small positive
! values pass through unchanged.  A separate FEDD guard below prevents
! division by a genuinely tiny JNU_STORE.
!
! edited by Sean Brennan
	DO I=1,ND
	  IF(JNU_STORE(I) .LT. 0.0_LDP .OR.
	1    JNU_STORE(I) .NE. JNU_STORE(I))THEN
!$OMP ATOMIC
! edited by Sean Brennan
	    N_NAN_JNU=N_NAN_JNU+1
	    JNU_STORE(I)=1.0E-30_LDP
!
! Phase 2 first-occurrence diagnostic for JNU failure.
! This catches NaN/negative J that propagated through the angular integration.
!$OMP CRITICAL(FIRST_FAIL_CRIT)
      IF(.NOT. FIRST_FAIL_FOUND)THEN
        FIRST_FAIL_FOUND=.TRUE.
        FIRST_FAIL_TYPE=3
        J=ERROR_LU()
        WRITE(J,'(/,A)')' ==================== PHASE 2: FIRST NUMERICAL FAILURE (JNU) ==================='
        WRITE(J,*)"Type: JNU_STORE negative or NaN (after angular integration)"
        WRITE(J,"(A,I8)")"  Freq index: ",FREQ_IDX_ITER
        WRITE(J,"(A,ES16.8)")"  Frequency: ",FREQ
        WRITE(J,"(A,I6)")"  Depth index (I): ",I
        WRITE(J,"(A,ES16.8)")"  JNU_STORE (failed): ",JNU_STORE(I)
        WRITE(J,"(A,ES16.8)")"  R: ",R(I)
        WRITE(J,"(A,ES16.8)")"  V: ",V(I)
        WRITE(J,"(A,ES16.8)")"  SIGMA_TRUE: ",SIGMA_TRUE_IN(I)
        WRITE(J,"(A,ES16.8)")"  IC (inner bnd): ",IC
        WRITE(J,"(A,ES16.8)")"  DBB: ",DBB
        WRITE(J,"(A,I6)")"  N_NEG_INT (this freq): ",N_NEG_INT
        WRITE(J,"(A,I6)")"  N_NAN_JNU (this freq): ",N_NAN_JNU
        WRITE(J,*)"  NOTE: JNU failure after I_M/I_P did NOT trigger first."
        WRITE(J,*)"  This means the angular integration produced the NaN."
        WRITE(J,'(/,A)')' ========================== END PHASE 2 FIRST FAILURE =========================='
      END IF
!$OMP END CRITICAL(FIRST_FAIL_CRIT)
! V4 diagnostic (2026-09-09): accumulate depth histogram
	    IF(ALLOCATED(JNU_CLAMP_DEPTH))JNU_CLAMP_DEPTH(I)=JNU_CLAMP_DEPTH(I)+1
! Sample frequency values (first 10000 events only)
	    IF(JNU_CLAMP_FREQ_CNT .LT. JNU_CLAMP_FREQ_MAX)THEN
	      JNU_CLAMP_FREQ_CNT=JNU_CLAMP_FREQ_CNT+1
	      JNU_CLAMP_FREQ_VALS(JNU_CLAMP_FREQ_CNT)=FREQ
	      JNU_CLAMP_FREQ_LAMS(JNU_CLAMP_FREQ_CNT)=2.99792458E+18_LDP/FREQ/1.0E+15_LDP
	    END IF
	  END IF
	END DO
!
! edited by Sean Brennan
! Boundary-closure audit.  Two distinct quantities share the name
! HBC: (i) the measured closure HNU_AT_OB/JNU_STORE(1), computed on
! the floored boundary J below and passed through to HBC_CMF so
! HFLUX_AT_OB keeps reporting the measured flux; (ii) the validated
! closure the moment solver's boundary row uses (HBC_ROW), where a
! negative or NaN value is replaced by the physical H(1)=0 closure
! and counted in MOM_BND_HBC_FL.  Flooring J first means the
! reported HBC can be finite-but-huge or negative when the formal
! solution underflows; FG_N_BADHBC counts those cases.  The
! asymmetry is deliberate: the raw value is an observable
! diagnostic, the floored value is the operator input.
!
	FG_MIN_JBND=MIN(FG_MIN_JBND,JNU_STORE(1))
	IF(JNU_STORE(1) .LT. 1.0E-30_LDP)JNU_STORE(1)=1.0E-30_LDP
!
! Inner-boundary Eddington IN_HBC = <mu*I_in>/(2*J(ND)-IC).  The old
! guard floored only |denom| < 1e-30, which is scale-blind: a small
! but nonzero signed denominator (J(ND) ~ IC/2 cancellation) produced
! an arbitrarily large or negative IN_HBC that then entered
! TB_NM(ND)=F/DTAU+IN_HBC and XM_NM(ND)=R^2*IC*(0.25+0.5*IN_HBC) in
! the moment solver.  Floor the denominator relative to the terms
! forming it, and validate the ratio into the physical range [0,1]
! (|mu| <= 1); an unusable measurement falls back to 0.5, the
! hemispheric mean of |mu| for incident radiation -- the same
! degenerate-Eddington-replaced-by-physical-closure policy as the
! MOM HBC_ROW floor.  The .NOT.(.GE.) tests also catch NaN.
! (Oct-2026 IN_HBC guard fix)
!
	DENOM_IN=2.0_LDP*JNU_STORE(ND)-IC
	DENOM_SCL=MAX(ABS(IC),2.0_LDP*JNU_STORE(ND),1.0E-30_LDP)
	IF(ABS(DENOM_IN) .LT. 1.0E-6_LDP*DENOM_SCL)THEN
	  IN_HBC=0.5_LDP
	ELSE
	  IN_HBC=IN_HBC/DENOM_IN
	END IF
	IF(.NOT.(IN_HBC .GE. 0.0_LDP))IN_HBC=0.5_LDP
	IF(.NOT.(IN_HBC .LE. 1.0_LDP))IN_HBC=0.5_LDP
	HBC=HNU_AT_OB/JNU_STORE(1)
	NBC=NNU_AT_OB/JNU_STORE(1)
	IF(.NOT. (ABS(HBC) .LE. HUGE(HBC)))THEN
	  FG_N_BADHBC=FG_N_BADHBC+1
	ELSE
	  FG_MAX_HBC=MAX(FG_MAX_HBC,ABS(HBC))
	  IF(HBC .LT. 0.0_LDP)FG_N_BADHBC=FG_N_BADHBC+1
	END IF
	RETURNED_OUT_HBC=HBC
	RETURNED_IN_HBC=IN_HBC
!
! Compute the Eddington factor, F. This required in CMFGEN to
! compute the K moment (from J computed by MOM_J_CMF).
!
! FEDD guard (2026-09-10): With the per-depth clamp removed (only NaN
! and negative are clamped now), JNU_STORE can be arbitrarily small at
! EUV frequencies.  KNU_STORE and JNU_STORE should be the same order
! of magnitude (both ~0 at EUV), so FEDD = KNU/JNU should be ~1/3.
! But if JNU_STORE is exactly zero or sub-tiny while KNU is not, FEDD
! could blow up.  Guard: if JNU_STORE < 1.0E-30, set FEDD = 1/3
! (isotropic limit) and JNU = 1.0E-30 for downstream consumers.
! This matches V13's boundary floor philosophy without the per-depth
! clamp that caused the false EUV cascade.
!
! edited by Sean Brennan
!
	DO I=1,ND
	  IF(JNU_STORE(I) .LT. 1.0E-30_LDP)THEN
	    JNU(I)=1.0E-30_LDP
	    FEDD(I)=1.0_LDP/3.0_LDP
	  ELSE
	    JNU(I)=JNU_STORE(I)
	    FEDD(I)=KNU_STORE(I)/JNU_STORE(I)
	  END IF
	END DO
!
! Store frequencies at which errors occurred, and give an indication of the
! error.
!
	DO J=1,N_ERR_MAX
	  IF(.NOT. NEG_AV_VALUE)EXIT
	  IF(FG_ERR_ON_FREQ(J) .EQ. FREQ)THEN
	    FG_ERR_TYPE(J)=FG_ERR_TYPE(J)+1
	    NEG_AV_VALUE=.FALSE.
	  ELSE IF(J .EQ. FG_ERR_CNT+1)THEN
	    FG_ERR_CNT=J
	    FG_ERR_TYPE(J)=1
	    FG_ERR_ON_FREQ(J)=FREQ
	    NEG_AV_VALUE=.FALSE.
	  END IF
	END DO
!
! This will only output IP_FG_DATA if a file IP_FG_DATA exists.
! This is only for debugging purposes.
!
        ACCESS_F=5
	LU_IP=56
        IF(INIT)THEN
	  INQUIRE(FILE='IP_FG_DATA',EXIST=WRITE_IP)
	  IF(WRITE_IP)THEN
	    CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
            I=WORD_SIZE*(NP+1)/UNIT_SIZE
            CALL WRITE_DIRECT_INFO_V3(NP,I,'20-Aug-2000','IP_FG_DATA',LU_IP)
            OPEN(UNIT=LU_IP,FILE='IP_FG_DATA',FORM='UNFORMATTED',
	1         ACCESS='DIRECT',STATUS='UNKNOWN',RECL=I,IOSTAT=IOS)
!
! A failed direct-access open left WRITE_IP armed: every later
! REC write then hits an unopened unit on the production path.
! Fail closed instead -- IP_FG_DATA is an opt-in diagnostic file,
! so disabling it costs nothing.
!
	    IF(IOS .NE. 0)THEN
	      WRITE_IP=.FALSE.
	      WRITE(LUER,'(A,I8)')'IP_FG_DATA open failed, IOSTAT=',IOS
	    END IF
	    IF(WRITE_IP)THEN
	      FREQ_CNT=0
              WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
              WRITE(LU_IP,REC=ACCESS_F)(P(I),I=1,NP)
	    END IF
	  END IF
	END IF
        IF(WRITE_IP)THEN
	  T1=0.0_LDP
          IF(FREQ_CNT .NE. 0)READ(LU_IP,REC=ACCESS_F+FREQ_CNT)(IBOUND(LS),LS=1,NP),T1
	  IF(T1 .NE. FREQ)FREQ_CNT=FREQ_CNT+1
          WRITE(LU_IP,REC=3)ACCESS_F,FREQ_CNT,NP
          WRITE(LU_IP,REC=ACCESS_F+FREQ_CNT)(CV_BOUND(LS),LS=1,NP),FREQ
        END IF
	CALL TUNE(2,'FG_BIG_SEC')
!
! Roll this frequency's clamp counts into the sweep totals so they are
! reported on the next INIT (see the block near the top of the routine).
! edited by Sean Brennan
!
	N_NEG_INT_SWEEP=N_NEG_INT_SWEEP+N_NEG_INT
	N_NAN_JNU_SWEEP=N_NAN_JNU_SWEEP+N_NAN_JNU
	N_SWEEP_CNT=N_SWEEP_CNT+1
!
! V4 diagnostic (2026-09-09): Write depth histogram and frequency samples
! to JNU_CLAMP_DIAG.txt after the first sweep of the first great iteration.
! This is a one-shot diagnostic - it writes once and never again.
! edited by Sean Brennan
! W18: gated on CONV_DIAG_LVL>=1 -- this write ran unconditionally on
! the production path on every cold start.  IOSTAT checked so a
! failing open never wedges the sweep.
	IF(CONV_DIAG_LVL .GE. 1 .AND.
	1   JNU_CLAMP_DIAG_WRITTEN .EQ. 0 .AND. N_SWEEP_CNT .GE. NCF_IN)THEN
	  JNU_CLAMP_DIAG_WRITTEN=1
	  OPEN(UNIT=999,FILE='JNU_CLAMP_DIAG.txt',STATUS='REPLACE',
	1   ACTION='WRITE',IOSTAT=IOS)
	  IF(IOS .NE. 0)GOTO 9998
	  WRITE(999,'(A)')' ========== JNU Clamp Diagnostic (first sweep, first great iteration) =========='
	  WRITE(999,'(A,I8,A,I6)')'Total clamped: ',N_NAN_JNU_SWEEP,' over ',N_SWEEP_CNT,' frequencies'
	  WRITE(999,*)
	  WRITE(999,'(A)')' ---------- Depth histogram (depth, count, fraction of frequencies) ------------'
	  IF(ALLOCATED(JNU_CLAMP_DEPTH))THEN
	    DO I=1,ND
	      IF(JNU_CLAMP_DEPTH(I) .GT. 0)THEN
	        WRITE(999,'(I6,I8,F10.4)')I,JNU_CLAMP_DEPTH(I),
	1           DBLE(JNU_CLAMP_DEPTH(I))/DBLE(N_SWEEP_CNT)
	      END IF
	    END DO
	  END IF
	  WRITE(999,*)
	  WRITE(999,'(A,I6,A)')' ---------- Frequency samples (first ',JNU_CLAMP_FREQ_CNT,' events) ----------------------------'
	  WRITE(999,'(A)')'  Freq(10^15 Hz)   Lambda(Ang)   Depth_index'
	  DO I=1,MIN(JNU_CLAMP_FREQ_CNT,100)
	    WRITE(999,'(ES14.6,ES14.2)')JNU_CLAMP_FREQ_VALS(I),JNU_CLAMP_FREQ_LAMS(I)
	  END DO
	  WRITE(999,*)
	  WRITE(999,'(A)')' ------------------- Frequency histogram (binned by lambda) --------------------'
	  WRITE(999,'(A)')'  Lambda range (Ang)    Count'
	  DO I=1,20
	    LAM_MIN=120.0_LDP*10.0_LDP**(I*0.3_LDP)
	    LAM_MAX=120.0_LDP*10.0_LDP**((I+1)*0.3_LDP)
	    KC=0
	    DO J=1,JNU_CLAMP_FREQ_CNT
	      IF(JNU_CLAMP_FREQ_LAMS(J) .GE. LAM_MIN .AND.
	1         JNU_CLAMP_FREQ_LAMS(J) .LT. LAM_MAX)KC=KC+1
	    END DO
	    IF(KC .GT. 0)WRITE(999,'(F10.1,A,F10.1,I8)')LAM_MIN,' -',LAM_MAX,KC
	  END DO
	  CLOSE(999)
	  WRITE(LUER,'(A)')'JNU_CLAMP_DIAG.txt written (one-shot diagnostic)'
	END IF
	GOTO 9997
9998	WRITE(LUER,'(A)')'JNU_CLAMP_DIAG.txt open failed -- diagnostic skipped'
!
9997	CONTINUE
	RETURN
	END
!
! ==============================================================================
!  Gauss-Seidel I_NEXT restart/checkpoint routines (added 16-Sep-2026).
!
!  The in-place GS arrays I_P/I_M_DEPTH_NEXT hold a COMPLETE generation only
!  between frequency sweeps.  They are checkpointed at every SCR_RITE_V2 call
!  (which always runs after the population solve, i.e. long after the sweep
!  completed), tagged with the NITSF value written to the POINT files.  On
!  restart, SCR_READ_V2 returns that same NITSF and CMFGEN_SUB passes it to
!  SET_FG_GS_TAG, so a file is accepted only if it was written by the same
!  checkpoint event the population restart used.
!
!  File FGINEXT_GS (unformatted sequential):
!    rec 1: MAGIC, VERSION, TAG
!    rec 2: N_COMP, NP, NCF_IN, ND
!    rec 3-5: R(1:ND), V(1:ND), SIGMA(1:ND)
!    rec 6: COMP_DEPTH_IDX(1:ND)
!    rec 7: FG_GS_MAXCHG, FG_INEXT_CAPPED,
!           [appended W1-W9]: FG_GS_NDIV, FG_COMMIT_REJ, FG_DIV_HEAL,
!           FG_SW_MAXCHG, FG_SW_NCONV, FG_RF_DCHG
!           (prefix-compatible: readers pulling the original two
!           fields still succeed on the longer record)
!    rec 8: KN_FREQ(1:NCF_IN)      (VERSION >= 3 only; the writer's
!           frequency grid -- needed so the OBS driver, which builds a
!           different grid, can remap I_NEXT in frequency)
!    rec 9: I_P_DEPTH_NEXT (N_COMP*NP*NCF_IN values; rec 8-10 skipped
!           if N_COMP=0)
!    rec 10: I_M_DEPTH_NEXT
!    rec 11: -MAGIC
!  Any read/validation failure leaves the arrays zeroed and I_NEXT_GS_OK
!  .FALSE. -> existing N_LAM_RWARM LAMBDA warm-up remains the fallback.
! ==============================================================================
!
	SUBROUTINE SET_FG_GS_TAG(N)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FG_GS_TAG_EXP
	IMPLICIT NONE
	INTEGER N
	FG_GS_TAG_EXP=N
	RETURN
	END
!
	SUBROUTINE SET_FG_GS_OBS_GRID(NU,NCF)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: GS_OBS_MODE,GS_NCF_OBS,GS_NU_OBS
	IMPLICIT NONE
	INTEGER NCF
	REAL(KIND=LDP) NU(NCF)
!
! Called by CMF_FLUX_SUB_V5 once its frequency grid is final and before
! the formal-solve loop starts.  The OBS grid differs from the model
! grid that wrote FGINEXT_GS, so the solver's first-call init must use
! RESTORE_FG_INEXT_OBS (frequency remap) instead of the positional
! RESTORE_FG_INEXT_GS.  NU is copied; the caller's array stays live.
!
	IF(ALLOCATED(GS_NU_OBS))DEALLOCATE(GS_NU_OBS)
	ALLOCATE(GS_NU_OBS(NCF))
	GS_NU_OBS(1:NCF)=NU(1:NCF)
	GS_NCF_OBS=NCF
	GS_OBS_MODE=.TRUE.
	RETURN
	END
!
! Return the restored-checkpoint maturity state for the N_LAM_RWARM
! warm-up decision in SOLVE_FOR_POPS.  Kept as a subroutine (not a
! module USE) because new_main/mod_subs compiles before plane/ in the
! library build order, so module visibility is not available there.
!
!   GS_OK    : a validated FGINEXT_GS generation was restored
!   CHG_FILE : last-sweep max relative change recorded in the file
!              (-1 when no checkpoint was read)
!   CAP_FILE : capped/non-finite write count recorded in the file
!
	SUBROUTINE GET_FG_GS_MATURITY(GS_OK,CHG_FILE,CAP_FILE,
	1                             MAXCHG,CAPPED,NBAD,NTOT,NCONV,
	2                             NDIV,NRAY,NBAD_BLK,
	3                             SW_MAXCHG,SW_NCONV,RF_DCHG,
	4                             NBAD_STAL,NBAD_PLAT,
	5                             SW_CAPPED,SW_CAPTOT,RF_CAPPED)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY:
	1     I_NEXT_GS_OK,FG_GS_CHG_FILE,FG_GS_CAP_FILE,
	2     FG_GS_MAXCHG,FG_INEXT_CAPPED,FG_GS_NBAD,FG_GS_NTOT,
	3     FG_GS_NCONV,FG_GS_NP,FG_GS_NDIV,FG_GS_NBAD_BLK,
	4     FG_SW_MAXCHG,FG_SW_NCONV,FG_RF_DCHG,
	5     FG_GS_NSTAL,FG_GS_NPLAT,
	6     FG_SW_CAPPED,FG_SW_CAPTOT,FG_RF_CAPPED
	IMPLICIT NONE
	LOGICAL GS_OK
	REAL(KIND=LDP) CHG_FILE,MAXCHG,SW_MAXCHG,RF_DCHG
	INTEGER CAP_FILE,CAPPED,NDIV,NRAY
	INTEGER*8 NBAD,NTOT,NCONV,NBAD_BLK,SW_NCONV
	INTEGER*8 NBAD_STAL,NBAD_PLAT
	INTEGER*8 SW_CAPPED,SW_CAPTOT,RF_CAPPED
	GS_OK=I_NEXT_GS_OK
	CHG_FILE=FG_GS_CHG_FILE
	CAP_FILE=FG_GS_CAP_FILE
	MAXCHG=FG_GS_MAXCHG
	CAPPED=FG_INEXT_CAPPED
	NBAD=FG_GS_NBAD
	NTOT=FG_GS_NTOT
	NCONV=FG_GS_NCONV
!
! W1-W9 (17-Oct-2026): appended transfer fields -- the diverged-ray
! count (share-gated in SOLVE_FOR_POPS), the solved-ray count for
! that share, the best-pass bad census (W3 bulk floor), and the
! sweep-semantics drift the refresh used to clobber (W4).
!
	NDIV=FG_GS_NDIV
	NRAY=FG_GS_NP
	NBAD_BLK=FG_GS_NBAD_BLK
	SW_MAXCHG=FG_SW_MAXCHG
	SW_NCONV=FG_SW_NCONV
	RF_DCHG=FG_RF_DCHG
	NBAD_STAL=FG_GS_NSTAL
	NBAD_PLAT=FG_GS_NPLAT
!
! Same-census cap fields: sweep-only caps against sweep-store
! writes, refresh-only caps against the refresh significant-cell
! total.  SOLVE_FOR_POPS picks the pair matching whichever side
! ran this iteration.
!
	SW_CAPPED=FG_SW_CAPPED
	SW_CAPTOT=FG_SW_CAPTOT
	RF_CAPPED=FG_RF_CAPPED
	RETURN
	END
!
	SUBROUTINE GET_FG_NOSLOT(CNT)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FG_NOSLOT_CNT
	IMPLICIT NONE
	INTEGER*8 CNT
	CNT=FG_NOSLOT_CNT
	RETURN
	END
!
	SUBROUTINE GET_FG_ALO_DIAG(DJDCHI,DJDETA,N,VALID)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY:
	1     ALO_DJDCHI,ALO_DJDETA,ALO_VALID,ALO_ND
	IMPLICIT NONE
	INTEGER N
	REAL(KIND=LDP) DJDCHI(N),DJDETA(N)
	LOGICAL VALID
!
! Accessor for the diagonal approximate-Lambda derivatives (dJ/dCHI,
! dJ/dETA at same depth) stored by FG_J_CMF_NONMONO_V1 for the current
! frequency.  Needed by DO_VAR_CONT's fixed-J branch in new_main, which
! compiles before plane and so cannot USE the module directly.
!
	VALID=.FALSE.
	IF(ALLOCATED(ALO_DJDCHI) .AND. ALLOCATED(ALO_DJDETA))THEN
	  IF(ALO_VALID .AND. ALO_ND .EQ. N)THEN
	    DJDCHI(1:N)=ALO_DJDCHI(1:N)
	    DJDETA(1:N)=ALO_DJDETA(1:N)
	    VALID=.TRUE.
	  END IF
	END IF
	RETURN
	END
!
	SUBROUTINE DUMP_FG_INEXT_GS(NITSF,R,V,SIGMA,ND,LU)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1
	USE IFPORT
	IMPLICIT NONE
	INTEGER NITSF,ND,LU
	REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND)
	INTEGER IOS,I,LUER,ERROR_LU
	INTEGER, PARAMETER :: MAGIC=20260916, VER=3
	LOGICAL REN_OK
	EXTERNAL ERROR_LU
!
! Skip entirely when the arrays are not allocated.  This covers the
! startup SCR_RITE (before the first formal solve): writing a header
! there would clobber a still-valid checkpoint file before the restore
! in the first FG call could read it.  With no compressing depths the
! arrays are never allocated and no file is ever written; the restore
! then reports "no checkpoint" and the warm-up behaves as before.
!
	IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT) .OR.
	1   .NOT. ALLOCATED(I_M_DEPTH_NEXT))RETURN
!
! Report the per-sweep Gauss-Seidel convergence of the I_NEXT coupling
! arrays (Knop et al. 2009 Paper II pathway): the largest relative
! change written this sweep plus the count of capped/bad writes
! (sweep write-back plus any refresh clamp/divergence events).  This
! report must not contain the literal string "NaN": submit_cmfgen.sh's
! check_outgen_for_fatal treats it as a crash signature and would
! veto every requeue.
!
	LUER=ERROR_LU()
	WRITE(LUER,'(/,A)')' --------------------------------- KNOP STORE ----------------------------------'
	WRITE(LUER,'(A)')'  I_NEXT sweep update:'
	WRITE(LUER,'(A,ES12.4)')
	1  '    max rel. change=',FG_GS_MAXCHG
	WRITE(LUER,'(A,I10)')'    cap/bad writes=',FG_INEXT_CAPPED
	WRITE(LUER,'(A,I10)')'    cells>50%=',FG_GS_NCONV
	WRITE(LUER,'(A,I8,2(A,I10))')
	1   '    div rays=',FG_GS_NDIV,
	2   '  healed slots=',FG_DIV_HEAL,
	3   '  commit rejs=',FG_COMMIT_REJ+FG_NC_REJ
	WRITE(LUER,'(A,I6)')'    NCF=',FG_GS_NCF
	WRITE(LUER,'(A,ES12.4)')'    max |HBC|=',FG_MAX_HBC
	WRITE(LUER,'(A,ES12.4)')'    min boundary J=',FG_MIN_JBND
	WRITE(LUER,'(A,I10)')'    bad HBC closures=',FG_N_BADHBC
	WRITE(LUER,'(A,I10)')'    floored bnd I+/-=',FG_N_NEG_BND_I
	WRITE(LUER,'(A,I10,A,ES10.2)')'    bnd J/K identity fails=',
	1   FG_N_BND_IDENT,'  max res=',FG_MAX_BND_RES
	WRITE(LUER,'(A)')' -------------------------------------------------------------------------------'
!
! Format version 2 (16-Sep-2026): the fingerprint is the FULL R, V,
! SIGMA and COMP_DEPTH_IDX vectors, not just endpoints/sums -- two
! different velocity structures can share endpoints and checksums while
! the Sigma < -1 regions that define the GS coupling move, so endpoint
! matching was not sufficient to tie the stored I_NEXT to this grid.
!
! Format version 3 (Oct-2026): adds the writer's frequency grid KN_FREQ
! between the metrics record and the intensity records.  The
! observation-frame driver (cmf_flux) builds a different frequency grid
! than the model run, so RESTORE_FG_INEXT_OBS needs the model NU array
! to remap the stored I_NEXT field in frequency.  KN_FREQ precedes the
! intensity records so the OBS restore can read the grid first and
! then stream each intensity array: read, remap, deallocate, keeping
! the transient footprint to one array.
!
! The file is written to FGINEXT_GS.tmp and atomically renamed onto
! FGINEXT_GS so a kill mid-write cannot destroy the previous valid
! checkpoint (STATUS='REPLACE' truncated the old file first).
!
	OPEN(UNIT=LU,FILE='FGINEXT_GS.tmp',FORM='UNFORMATTED',
	1     STATUS='REPLACE',IOSTAT=IOS)
	IF(IOS .NE. 0)RETURN
	WRITE(LU,IOSTAT=IOS)MAGIC,VER,NITSF
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)N_COMP,FG_GS_NP,FG_GS_NCF,ND
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)(R(I),I=1,ND)
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)(V(I),I=1,ND)
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)(SIGMA(I),I=1,ND)
	IF(IOS .EQ. 0)THEN
	  IF(ALLOCATED(COMP_DEPTH_IDX))THEN
	    WRITE(LU,IOSTAT=IOS)(COMP_DEPTH_IDX(I),I=1,ND)
	  ELSE
	    WRITE(LU,IOSTAT=IOS)(0,I=1,ND)
	  END IF
	END IF
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)FG_GS_MAXCHG,
	1     INT(MIN(FG_INEXT_CAPPED,2147483647_8),4),
	2     FG_GS_NDIV,FG_COMMIT_REJ,FG_DIV_HEAL,
	3     FG_SW_MAXCHG,FG_SW_NCONV,FG_RF_DCHG
	IF(IOS .EQ. 0 .AND. N_COMP .GT. 0)THEN
	  WRITE(LU,IOSTAT=IOS)KN_FREQ
	  IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)I_P_DEPTH_NEXT
	  IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)I_M_DEPTH_NEXT
	END IF
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)-MAGIC
	CLOSE(LU)
!
! Only rename onto the live name if every record was written; a torn
! .tmp file is simply left behind and the old checkpoint stays valid.
!
	IF(IOS .EQ. 0)THEN
	  REN_OK=RENAMEFILEQQ('FGINEXT_GS.tmp','FGINEXT_GS')
	END IF
	RETURN
	END
!
	SUBROUTINE RESTORE_FG_INEXT_GS(R,V,SIGMA,ND,NP,NCF_IN)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1
	IMPLICIT NONE
	INTEGER ND,NP,NCF_IN
	REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND)
	INTEGER IOS,LUER,ERROR_LU
	INTEGER MAG,VER,TAG
	INTEGER N_COMP_F,NP_F,NCF_F,ND_F,ENDM,I
	INTEGER, PARAMETER :: MAGIC=20260916, VER_MIN=2, VER_MAX=3
	REAL(KIND=LDP), ALLOCATABLE :: R_F(:),V_F(:),SG_F(:)
	INTEGER, ALLOCATABLE :: CMP_F(:)
	EXTERNAL ERROR_LU
	LOGICAL FILE_THERE
!
	LUER=ERROR_LU()
	I_NEXT_GS_OK=.FALSE.
	FG_GS_CHG_FILE=-1.0_LDP
	FG_GS_CAP_FILE=-1
!
	WRITE(LUER,'(A)')'GS RESTART STATE (FG_J_CMF_NONMONO_V1):'
!
! Hard invariant: NEWMOD (or any state where SCR_READ_V2 did not return
! a positive NITSF) never restores GS state.  This cannot rely on the
! tag check alone: a stale file from a different model could carry any
! tag, so the cold-start case exits before the file is even opened.
!
	IF(FG_GS_TAG_EXP .LE. 0)THEN
	  WRITE(LUER,'(A)')'  cold start (NITSF<=0): no GS restore'
	  WRITE(LUER,'(A)')'  FALLING BACK TO LAMBDA REBUILD (I_NEXT=0)'
	  RETURN
	END IF
!
	INQUIRE(FILE='FGINEXT_GS',EXIST=FILE_THERE)
	IF(.NOT. FILE_THERE)THEN
	  WRITE(LUER,'(A)')'  I_NEXT checkpoint present: NO'
	  WRITE(LUER,'(A)')'  FALLING BACK TO LAMBDA REBUILD (I_NEXT=0)'
	  RETURN
	END IF
	WRITE(LUER,'(A)')'  I_NEXT checkpoint present: YES'
!
	OPEN(UNIT=291,FILE='FGINEXT_GS',FORM='UNFORMATTED',
	1     STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)MAG,VER,TAG
	IF(IOS .NE. 0 .OR. MAG .NE. MAGIC .OR.
	1   VER .LT. VER_MIN .OR. VER .GT. VER_MAX)GOTO 900
	READ(291,IOSTAT=IOS)N_COMP_F,NP_F,NCF_F,ND_F
	IF(IOS .NE. 0)GOTO 900
!
! Validate: generation tag must equal the NITSF SCR_READ_V2 restored
! (ties this file to the same checkpoint event as SCRTEMP), and the
! dimensions must match before the grid fingerprint is read.
!
	WRITE(LUER,'(A,I8)')'  GS tag (checkpoint iter): ',TAG
	WRITE(LUER,'(A,I8)')'  expected tag (NITSF):     ',FG_GS_TAG_EXP
	IF(TAG .NE. FG_GS_TAG_EXP)THEN
	  WRITE(LUER,'(A)')'  I_NEXT valid: NO (stale/foreign checkpoint)'
	  GOTO 800
	END IF
	IF(N_COMP_F .NE. N_COMP .OR. NP_F .NE. NP .OR.
	1   NCF_F .NE. NCF_IN .OR. ND_F .NE. ND)THEN
	  WRITE(LUER,'(A)')'  I_NEXT valid: NO (dimension mismatch)'
	  WRITE(LUER,'(A,4I8)')'   file: N_COMP,NP,NCF,ND =',
	1                     N_COMP_F,NP_F,NCF_F,ND_F
	  WRITE(LUER,'(A,4I8)')'   run :                  =',
	1                     N_COMP,NP,NCF_IN,ND
	  GOTO 800
	END IF
!
! Full-grid fingerprint: every element of R, V and SIGMA must be
! identical to the file's, then the entire COMP_DEPTH_IDX map.  This
! guarantees I_NEXT(KC,LS,K) refers to the same physical depth, ray and
! frequency after the requeue -- endpoint/checksum matching cannot.
! Each unformatted READ consumes a whole record, so the file vectors
! are read into matching-shape locals, not elementwise.
!
	ALLOCATE(R_F(ND),V_F(ND),SG_F(ND),CMP_F(ND),STAT=IOS)
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)R_F
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)V_F
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)SG_F
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)CMP_F
	IF(IOS .NE. 0)GOTO 900
	IF(ANY(R_F .NE. R) .OR. ANY(V_F .NE. V) .OR.
	1   ANY(SG_F .NE. SIGMA))GOTO 810
	IF(ALLOCATED(COMP_DEPTH_IDX))THEN
	  IF(ANY(CMP_F .NE. COMP_DEPTH_IDX))GOTO 820
	ELSE
	  IF(ANY(CMP_F .NE. 0))GOTO 820
	END IF
	DEALLOCATE(R_F,V_F,SG_F,CMP_F)
!
! Last-sweep convergence metrics recorded at dump time.
!
	READ(291,IOSTAT=IOS)FG_GS_CHG_FILE,FG_GS_CAP_FILE
	IF(IOS .NE. 0)GOTO 900
!
! VERSION >= 3 files carry the writer's frequency grid ahead of the
! intensity records.  On this path the grid must equal the run grid
! (same NCF enforced above), so it is read only to keep the stream
! aligned; a non-positive or non-monotone grid marks the file corrupt.
!
	IF(VER .GE. 3 .AND. N_COMP .GT. 0)THEN
	  ALLOCATE(R_F(NCF_F),STAT=IOS)
	  IF(IOS .NE. 0)GOTO 900
	  READ(291,IOSTAT=IOS)R_F
	  IF(IOS .NE. 0)GOTO 900
	  IF(ANY(R_F .LE. 0.0_LDP) .OR.
	1   ANY(R_F(2:NCF_F) .GE. R_F(1:NCF_F-1)))GOTO 830
	  DEALLOCATE(R_F)
	END IF
!
! Validated header.  Read the two arrays (they share one file/write so
! they are always the same GS generation), then the end marker.
!
	IF(N_COMP .GT. 0)THEN
	  IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT))GOTO 900
	  READ(291,IOSTAT=IOS)I_P_DEPTH_NEXT
	  IF(IOS .NE. 0)GOTO 900
	  READ(291,IOSTAT=IOS)I_M_DEPTH_NEXT
	  IF(IOS .NE. 0)GOTO 900
	END IF
	READ(291,IOSTAT=IOS)ENDM
	IF(IOS .NE. 0 .OR. ENDM .NE. -MAGIC)GOTO 900
!
! Value validation (audit 17-Sep-2026): a structurally intact file can
! still carry a corrupted field -- the convergence metric must be
! finite and non-negative, the cap/floor count non-negative, and the
! restored intensities finite and non-negative everywhere.  Any
! violation rejects the whole generation: the next compressing-
! frequency solve would otherwise read the corrupted values as
! explicit upwind source terms.
!
	IF(.NOT.(ABS(FG_GS_CHG_FILE) .LE. HUGE(1.0_LDP)) .OR.
	1   FG_GS_CHG_FILE .LT. 0.0_LDP .OR.
	1   FG_GS_CAP_FILE .LT. 0)GOTO 830
	IF(N_COMP .GT. 0)THEN
	  IF(ANY(.NOT.(ABS(I_P_DEPTH_NEXT) .LE. HUGE(1.0_LDP))) .OR.
	1     ANY(.NOT.(ABS(I_M_DEPTH_NEXT) .LE. HUGE(1.0_LDP))) .OR.
	1     ANY(I_P_DEPTH_NEXT .LT. 0.0_LDP) .OR.
	1     ANY(I_M_DEPTH_NEXT .LT. 0.0_LDP))GOTO 830
	END IF
	CLOSE(291)
!
	I_NEXT_GS_OK=.TRUE.
	WRITE(LUER,'(A)')'  I_NEXT valid: YES'
	WRITE(LUER,'(A,I8,A)')
	1 '  RESTORED COMPLETE I_NEXT GENERATION (tag=',TAG,')'
	WRITE(LUER,'(A,I6,A,3I8)')'  N_COMP=',N_COMP,
	1 '  NP,NCF,ND=',NP,NCF_IN,ND
	WRITE(LUER,'(A,ES10.3,A,I6)')
	1 '  last-sweep max change=',FG_GS_CHG_FILE,
	1 '  cap/floor events=',FG_GS_CAP_FILE
	RETURN
!
! File unreadable/torn mid-record: treated as invalid.
!
900	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (read/format failure)'
	GOTO 800
810	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (grid fingerprint mismatch)'
	GOTO 800
820	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (COMP_DEPTH map mismatch)'
	GOTO 800
830	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (non-finite/negative data)'
800	CONTINUE
	INQUIRE(UNIT=291,OPENED=FILE_THERE)
	IF(FILE_THERE)CLOSE(291)
	I_NEXT_GS_OK=.FALSE.
	FG_GS_CHG_FILE=-1.0_LDP
	FG_GS_CAP_FILE=-1
	IF(ALLOCATED(I_P_DEPTH_NEXT))I_P_DEPTH_NEXT=0.0_LDP
	IF(ALLOCATED(I_M_DEPTH_NEXT))I_M_DEPTH_NEXT=0.0_LDP
	WRITE(LUER,'(A)')'  FALLING BACK TO LAMBDA REBUILD (I_NEXT=0)'
	RETURN
	END
!
! ------------------------------------------------------------------------------
! In-run rollback restore (30-Sep-2026).
!
! RESTORE_FG_INEXT_GS is gated on FG_GS_TAG_EXP, the NITSF restored at
! job start, so it can only run at INIT.  The resacc rejection in
! SOLVEBA needs the same restore mid-run: when a BA step is rolled
! back, the lagged I_NEXT field in memory was solved on the post-step
! (rejected) populations, and re-solving from that seed can leave the
! GS basin entirely -- LAM_dbg it49-56: pops rollback alone left
! se_res2 growing 439 -> 2e7 at held pops.  At the reject point the
! FGINEXT_GS on disk still holds the generation checkpointed at the
! end of the step's own iteration, i.e. the field solved on the
! pre-step populations the rollback restores -- the consistent pair.
!
! The NITSF tag check is meaningless mid-run (the file was written by
! this process), so this variant skips it and relies on the same
! R/V/SIGMA/COMP_DEPTH_IDX fingerprint the INIT restore uses.  A
! foreign or init-rejected FGINEXT_GS with matching dimensions can
! survive in the run directory (the INIT path rejects without
! deleting), and a velocity/grid revision between the last dump and
! the rollback is exactly the case the fingerprint exists to catch:
! accepting the file would pair the restored populations with a
! field solved on a different velocity structure.
!
! On failure the in-memory field is left untouched: a stale field
! after a pops rollback is no worse than the pops-only rollback the
! caller falls back to, while zeroing it (as the INIT path does)
! would guarantee a cold rebuild.
!
	SUBROUTINE RESTORE_FG_INEXT_LATEST(OK)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1
	USE MOD_CMFGEN, ONLY: R,V,SIGMA
	IMPLICIT NONE
	LOGICAL OK
	INTEGER IOS,LUER,ERROR_LU
	LOGICAL LFRQ_OK
	INTEGER MAG,VER,TAG
	INTEGER N_COMP_F,NP_F,NCF_F,ND_F,ENDM
	INTEGER, PARAMETER :: MAGIC=20260916, VER_MIN=2, VER_MAX=3
	REAL(KIND=LDP), ALLOCATABLE :: WRK(:)
	INTEGER, ALLOCATABLE :: IWRK(:)
	EXTERNAL ERROR_LU
	LOGICAL FILE_THERE
!
	OK=.FALSE.
	LUER=ERROR_LU()
!
	IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT) .OR.
	1   .NOT. ALLOCATED(I_M_DEPTH_NEXT))RETURN
	INQUIRE(FILE='FGINEXT_GS',EXIST=FILE_THERE)
	IF(.NOT. FILE_THERE)RETURN
!
	OPEN(UNIT=291,FILE='FGINEXT_GS',FORM='UNFORMATTED',
	1     STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)MAG,VER,TAG
	IF(IOS .NE. 0 .OR. MAG .NE. MAGIC .OR.
	1   VER .LT. VER_MIN .OR. VER .GT. VER_MAX)GOTO 900
	READ(291,IOSTAT=IOS)N_COMP_F,NP_F,NCF_F,ND_F
	IF(IOS .NE. 0)GOTO 900
	IF(N_COMP_F .NE. N_COMP .OR. NP_F .NE. FG_GS_NP .OR.
	1   NCF_F .NE. FG_GS_NCF .OR. ND_F .NE. KN_ND)GOTO 900
!
	ALLOCATE(WRK(ND_F),IWRK(ND_F),STAT=IOS)
	IF(IOS .NE. 0)GOTO 900
!
! Fingerprint records: each must equal the live grid exactly (the
! file holds a verbatim copy of the model arrays, so bitwise equality
! is the right test).  MOD_CMFGEN R/V/SIGMA are the arrays CMFGEN_SUB
! passes to SCR_RITE_V2, which is what DUMP_FG_INEXT_GS wrote.
!
	IF(.NOT. ALLOCATED(R) .OR. .NOT. ALLOCATED(V) .OR.
	1   .NOT. ALLOCATED(SIGMA))GOTO 920
	IF(SIZE(R) .NE. ND_F .OR. SIZE(V) .NE. ND_F .OR.
	1   SIZE(SIGMA) .NE. ND_F)GOTO 920
	READ(291,IOSTAT=IOS)WRK
	IF(IOS .NE. 0)GOTO 910
	IF(ANY(WRK .NE. R))GOTO 920
	READ(291,IOSTAT=IOS)WRK
	IF(IOS .NE. 0)GOTO 910
	IF(ANY(WRK .NE. V))GOTO 920
	READ(291,IOSTAT=IOS)WRK
	IF(IOS .NE. 0)GOTO 910
	IF(ANY(WRK .NE. SIGMA))GOTO 920
	READ(291,IOSTAT=IOS)IWRK
	IF(IOS .NE. 0)GOTO 910
	IF(ALLOCATED(COMP_DEPTH_IDX))THEN
	  IF(ANY(IWRK .NE. COMP_DEPTH_IDX))GOTO 930
	ELSE
	  IF(ANY(IWRK .NE. 0))GOTO 930
	END IF
	DEALLOCATE(WRK,IWRK)
!
	READ(291,IOSTAT=IOS)FG_GS_CHG_FILE,FG_GS_CAP_FILE
	IF(IOS .NE. 0)GOTO 900
!
	IF(VER .GE. 3 .AND. N_COMP .GT. 0)THEN
	  ALLOCATE(WRK(NCF_F),STAT=IOS)
	  IF(IOS .NE. 0)GOTO 900
	  READ(291,IOSTAT=IOS)WRK
	  IF(IOS .NE. 0)THEN
	    DEALLOCATE(WRK)
	    GOTO 900
	  END IF
	  LFRQ_OK=ALLOCATED(KN_FREQ)
	  IF(LFRQ_OK)LFRQ_OK=(SIZE(KN_FREQ) .EQ. NCF_F)
	  IF(LFRQ_OK)LFRQ_OK=.NOT. ANY(WRK .NE. KN_FREQ)
	  IF(.NOT. LFRQ_OK)THEN
	    DEALLOCATE(WRK)
	    GOTO 920
	  END IF
	  DEALLOCATE(WRK)
	END IF
!
	IF(N_COMP .GT. 0)THEN
	  READ(291,IOSTAT=IOS)I_P_DEPTH_NEXT
	  IF(IOS .NE. 0)GOTO 900
	  READ(291,IOSTAT=IOS)I_M_DEPTH_NEXT
	  IF(IOS .NE. 0)GOTO 900
	END IF
	READ(291,IOSTAT=IOS)ENDM
	IF(IOS .NE. 0 .OR. ENDM .NE. -MAGIC)GOTO 900
!
	CLOSE(291)
	I_NEXT_GS_OK=.TRUE.
	OK=.TRUE.
	WRITE(LUER,'(A,I8)')
	1  '  I_NEXT rollback restore accepted (tag=',TAG
	RETURN
!
900	CONTINUE
	INQUIRE(UNIT=291,OPENED=FILE_THERE)
	IF(FILE_THERE)CLOSE(291)
	WRITE(LUER,'(A)')
	1  '  I_NEXT rollback restore failed -- keeping in-memory field'
	RETURN
910	CONTINUE
	DEALLOCATE(WRK,IWRK)
	GOTO 900
920	CONTINUE
	WRITE(LUER,'(A)')
	1  '  I_NEXT rollback restore rejected: grid fingerprint mismatch'
	IF(ALLOCATED(WRK))DEALLOCATE(WRK)
	IF(ALLOCATED(IWRK))DEALLOCATE(IWRK)
	GOTO 900
930	CONTINUE
	WRITE(LUER,'(A)')
	1  '  I_NEXT rollback restore rejected: COMP_DEPTH map mismatch'
	IF(ALLOCATED(WRK))DEALLOCATE(WRK)
	IF(ALLOCATED(IWRK))DEALLOCATE(IWRK)
	GOTO 900
	END
!
! ------------------------------------------------------------------------------
! Observation-frame restore (Oct-2026).
!
! cmf_flux builds its own frequency grid, so the FGINEXT_GS checkpoint
! (written on the model grid) cannot be restored by position.  This
! routine remaps the stored I_NEXT field in frequency onto the OBS grid
! NU_OBS.  Store convention: slot s holds the field of frequency s+1
! on the caller's grid, so OBS slot s is filled with the file field at
! NU_OBS(s+1); slot NCF_IN is never read by the sweep and is zeroed.
! The (KC,LS) indices carry over unchanged: cmf_flux takes NP from RVTJ
! and builds rays via the same IMPAR on the same R, so rays and
! compressing depths coincide with the model's; the R/V/SIGMA/
! COMP_DEPTH_IDX fingerprint below enforces it.
!
! Interpolation is linear in log(I) vs log(nu) between bracketing file
! frequencies, clamped to the file edge values outside the file grid --
! the same log-space convention BUILD_FREQ uses when interpolating
! lagged intensities onto a ray.
!
! Requires a VERSION >= 3 checkpoint (frequency-grid record); older
! files are rejected loudly rather than silently cold-starting.  A
! failed restore leaves I_NEXT zeroed, i.e. the same cold start the
! OBS path had before this routine existed.
!
	SUBROUTINE RESTORE_FG_INEXT_OBS(R,V,SIGMA,ND,NP,NU_OBS,NCF_IN)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1
	IMPLICIT NONE
	INTEGER ND,NP,NCF_IN
	REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND),NU_OBS(NCF_IN)
	INTEGER IOS,LUER,ERROR_LU
	INTEGER MAG,VER,TAG
	INTEGER N_COMP_F,NP_F,NCF_F,ND_F,ENDM
	INTEGER KC,LS,IS,JF
	INTEGER, PARAMETER :: MAGIC=20260916, VER_MIN=3
	REAL(KIND=LDP), ALLOCATABLE :: R_F(:),V_F(:),SG_F(:),NU_F(:)
	REAL(KIND=LDP), ALLOCATABLE :: I_F(:,:,:)
	INTEGER, ALLOCATABLE :: CMP_F(:)
	REAL(KIND=LDP) TGT,W
	REAL(KIND=LDP), PARAMETER :: I_FLOOR=1.0E-40_LDP
	REAL(KIND=LDP), PARAMETER :: FP_TOL=1.0E-6_LDP
	EXTERNAL ERROR_LU
	LOGICAL FILE_THERE
	CHARACTER(LEN=20) FNAME
!
	LUER=ERROR_LU()
	I_NEXT_GS_OK=.FALSE.
	FG_GS_CHG_FILE=-1.0_LDP
	FG_GS_CAP_FILE=-1
!
	WRITE(LUER,'(A)')'GS RESTART STATE (FG_J_CMF_NONMONO_V1, OBS grid):'
!
	IF(SIZE(NU_OBS) .NE. GS_NCF_OBS .OR. NCF_IN .LT. 2)GOTO 840
	IF(NU_OBS(1) .LE. NU_OBS(NCF_IN) .OR.
	1   ANY(NU_OBS .LE. 0.0_LDP))GOTO 840
!
! batobs.sh runs in the model's obs/ subdirectory and reads the parent
! run directory's files (../RVTJ), so try both locations.
!
	FNAME='FGINEXT_GS'
	INQUIRE(FILE=FNAME,EXIST=FILE_THERE)
	IF(.NOT. FILE_THERE)THEN
	  FNAME='../FGINEXT_GS'
	  INQUIRE(FILE=FNAME,EXIST=FILE_THERE)
	END IF
	IF(.NOT. FILE_THERE)THEN
	  WRITE(LUER,'(A)')'  I_NEXT checkpoint present: NO'
	  GOTO 800
	END IF
	WRITE(LUER,'(A,A)')'  I_NEXT checkpoint present: YES  ',FNAME
!
	OPEN(UNIT=291,FILE=FNAME,FORM='UNFORMATTED',
	1     STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)MAG,VER,TAG
	IF(IOS .NE. 0 .OR. MAG .NE. MAGIC)GOTO 900
	IF(VER .LT. VER_MIN)THEN
	  WRITE(LUER,'(A,I2,A)')'  I_NEXT valid: NO (VERSION=',VER,
	1   ' has no frequency grid; rerun FINAL once to rewrite it)'
	  GOTO 800
	END IF
	READ(291,IOSTAT=IOS)N_COMP_F,NP_F,NCF_F,ND_F
	IF(IOS .NE. 0)GOTO 900
!
! No NITSF tag exists in cmf_flux (SCR_READ_V2 is not called); the
! R/V/SIGMA/COMP_DEPTH_IDX fingerprint ties the file to this model
! instead.  NCF_F is expected to differ -- it sizes the file arrays.
!
	WRITE(LUER,'(A,I8)')'  GS tag (checkpoint iter): ',TAG
	IF(N_COMP_F .NE. N_COMP .OR. NP_F .NE. NP .OR.
	1   ND_F .NE. ND)THEN
	  WRITE(LUER,'(A)')'  I_NEXT valid: NO (dimension mismatch)'
	  WRITE(LUER,'(A,4I8)')'   file: N_COMP,NP,NCF,ND =',
	1                     N_COMP_F,NP_F,NCF_F,ND_F
	  WRITE(LUER,'(A,4I8)')'   run :                  =',
	1                     N_COMP,NP,NCF_IN,ND
	  GOTO 800
	END IF
!
	ALLOCATE(R_F(ND),V_F(ND),SG_F(ND),CMP_F(ND),STAT=IOS)
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)R_F
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)V_F
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)SG_F
	IF(IOS .NE. 0)GOTO 900
	READ(291,IOSTAT=IOS)CMP_F
	IF(IOS .NE. 0)GOTO 900
!
! This path re-reads R/V/SIGMA from RVTJ, which is ASCII (~11
! significant digits): the values differ from the binary checkpoint's
! by round-trip error, so a bitwise fingerprint can never match here.
! Compare within a relative tolerance scaled by each vector's own
! maximum instead -- far above round-trip (~1e-11) yet far below any
! real grid change.  The integer COMP_DEPTH_IDX map stays exact.
!
	IF(ANY(ABS(R_F-R) .GT. FP_TOL*MAXVAL(ABS(R))) .OR.
	1   ANY(ABS(V_F-V) .GT. FP_TOL*MAXVAL(ABS(V))) .OR.
	2   ANY(ABS(SG_F-SIGMA) .GT. FP_TOL*MAXVAL(ABS(SIGMA))))THEN
	  WRITE(LUER,'(A,3ES10.2)')'  max scaled grid diff (R,V,SIG):',
	1    MAXVAL(ABS(R_F-R))/MAXVAL(ABS(R)),
	2    MAXVAL(ABS(V_F-V))/MAXVAL(ABS(V)),
	3    MAXVAL(ABS(SG_F-SIGMA))/MAXVAL(ABS(SIGMA))
	  GOTO 810
	END IF
	IF(ALLOCATED(COMP_DEPTH_IDX))THEN
	  IF(ANY(CMP_F .NE. COMP_DEPTH_IDX))GOTO 820
	ELSE
	  IF(ANY(CMP_F .NE. 0))GOTO 820
	END IF
	DEALLOCATE(R_F,V_F,SG_F,CMP_F)
!
	READ(291,IOSTAT=IOS)FG_GS_CHG_FILE,FG_GS_CAP_FILE
	IF(IOS .NE. 0)GOTO 900
	IF(.NOT.(ABS(FG_GS_CHG_FILE) .LE. HUGE(1.0_LDP)) .OR.
	1   FG_GS_CHG_FILE .LT. 0.0_LDP .OR.
	2   FG_GS_CAP_FILE .LT. 0)GOTO 830
!
	IF(N_COMP .GT. 0)THEN
	  ALLOCATE(NU_F(NCF_F),STAT=IOS)
	  IF(IOS .NE. 0)GOTO 900
	  READ(291,IOSTAT=IOS)NU_F
	  IF(IOS .NE. 0)GOTO 900
	  IF(.NOT. ANY(ABS(NU_F) .LE. HUGE(1.0_LDP)) .OR.
	1     ANY(NU_F .LE. 0.0_LDP) .OR.
	2     ANY(NU_F(2:NCF_F) .GE. NU_F(1:NCF_F-1)))GOTO 830
	  IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT) .OR.
	1     .NOT. ALLOCATED(I_M_DEPTH_NEXT))GOTO 900
	  ALLOCATE(I_F(N_COMP,NP,NCF_F),STAT=IOS)
	  IF(IOS .NE. 0)GOTO 900
!
! Stream each intensity array: read, validate, remap, free.  Reading
! both at once would double the transient footprint.
!
	  READ(291,IOSTAT=IOS)I_F
	  IF(IOS .NE. 0)GOTO 900
	  IF(ANY(.NOT.(ABS(I_F) .LE. HUGE(1.0_LDP))) .OR.
	1     ANY(I_F .LT. 0.0_LDP))GOTO 830
!$OMP PARALLEL DO PRIVATE(LS,KC,IS,JF,TGT,W) SCHEDULE(STATIC)
	  DO LS=1,NP
	    DO KC=1,N_COMP
	      JF=1
	      DO IS=1,NCF_IN-1
	        TGT=NU_OBS(IS+1)
	        IF(TGT .GE. NU_F(1))THEN
	          I_P_DEPTH_NEXT(KC,LS,IS)=I_F(KC,LS,1)
	        ELSE
	          DO WHILE(JF .LT. NCF_F)
	            IF(NU_F(JF+1) .LE. TGT)EXIT
	            JF=JF+1
	          END DO
	          IF(JF .GE. NCF_F)THEN
	            I_P_DEPTH_NEXT(KC,LS,IS)=I_F(KC,LS,NCF_F)
	          ELSE
	            W=LOG(NU_F(JF)/TGT)/LOG(NU_F(JF)/NU_F(JF+1))
	            I_P_DEPTH_NEXT(KC,LS,IS)=EXP(
	1             (1.0_LDP-W)*LOG(MAX(I_F(KC,LS,JF),I_FLOOR))+
	2             W*LOG(MAX(I_F(KC,LS,JF+1),I_FLOOR)))
	          END IF
	        END IF
	      END DO
	      I_P_DEPTH_NEXT(KC,LS,NCF_IN)=0.0_LDP
	    END DO
	  END DO
!$OMP END PARALLEL DO
!
	  READ(291,IOSTAT=IOS)I_F
	  IF(IOS .NE. 0)GOTO 900
	  IF(ANY(.NOT.(ABS(I_F) .LE. HUGE(1.0_LDP))) .OR.
	1     ANY(I_F .LT. 0.0_LDP))GOTO 830
!$OMP PARALLEL DO PRIVATE(LS,KC,IS,JF,TGT,W) SCHEDULE(STATIC)
	  DO LS=1,NP
	    DO KC=1,N_COMP
	      JF=1
	      DO IS=1,NCF_IN-1
	        TGT=NU_OBS(IS+1)
	        IF(TGT .GE. NU_F(1))THEN
	          I_M_DEPTH_NEXT(KC,LS,IS)=I_F(KC,LS,1)
	        ELSE
	          DO WHILE(JF .LT. NCF_F)
	            IF(NU_F(JF+1) .LE. TGT)EXIT
	            JF=JF+1
	          END DO
	          IF(JF .GE. NCF_F)THEN
	            I_M_DEPTH_NEXT(KC,LS,IS)=I_F(KC,LS,NCF_F)
	          ELSE
	            W=LOG(NU_F(JF)/TGT)/LOG(NU_F(JF)/NU_F(JF+1))
	            I_M_DEPTH_NEXT(KC,LS,IS)=EXP(
	1             (1.0_LDP-W)*LOG(MAX(I_F(KC,LS,JF),I_FLOOR))+
	2             W*LOG(MAX(I_F(KC,LS,JF+1),I_FLOOR)))
	          END IF
	        END IF
	      END DO
	      I_M_DEPTH_NEXT(KC,LS,NCF_IN)=0.0_LDP
	    END DO
	  END DO
!$OMP END PARALLEL DO
	  DEALLOCATE(I_F,NU_F)
	END IF
	READ(291,IOSTAT=IOS)ENDM
	IF(IOS .NE. 0 .OR. ENDM .NE. -MAGIC)GOTO 900
	CLOSE(291)
!
	I_NEXT_GS_OK=.TRUE.
	WRITE(LUER,'(A)')'  I_NEXT valid: YES'
	WRITE(LUER,'(A,I8,A)')
	1 '  RESTORED COMPLETE I_NEXT GENERATION (tag=',TAG,')'
	WRITE(LUER,'(A,I6,A,I8,A,I8)')
	1 '  N_COMP=',N_COMP,'  OBS remap: file NCF=',NCF_F,
	2 '  -> OBS NCF=',NCF_IN
	WRITE(LUER,'(A,ES10.3,A,I6)')
	1 '  last-sweep max change=',FG_GS_CHG_FILE,
	1 '  cap/floor events=',FG_GS_CAP_FILE
	RETURN
!
! File unreadable/torn mid-record: treated as invalid.
!
900	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (read/format failure)'
	GOTO 800
810	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (grid fingerprint mismatch)'
	GOTO 800
820	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (COMP_DEPTH map mismatch)'
	GOTO 800
830	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (non-finite/negative data)'
	GOTO 800
840	CONTINUE
	WRITE(LUER,'(A)')'  I_NEXT valid: NO (bad OBS frequency grid)'
800	CONTINUE
	INQUIRE(UNIT=291,OPENED=FILE_THERE)
	IF(FILE_THERE)CLOSE(291)
	I_NEXT_GS_OK=.FALSE.
	FG_GS_CHG_FILE=-1.0_LDP
	FG_GS_CAP_FILE=-1
	IF(ALLOCATED(I_P_DEPTH_NEXT))I_P_DEPTH_NEXT=0.0_LDP
	IF(ALLOCATED(I_M_DEPTH_NEXT))I_M_DEPTH_NEXT=0.0_LDP
	IF(ALLOCATED(I_F))DEALLOCATE(I_F)
	IF(ALLOCATED(NU_F))DEALLOCATE(NU_F)
	WRITE(LUER,'(A)')'  FALLING BACK TO COLD START (I_NEXT=0)'
	RETURN
	END
!
! ==============================================================================
!
! Per-ray Knop refresh routines (Oct-2026).
!
! FG_KNOP_REFRESH_V1 is called once per great iteration from CMFGEN_SUB
! after the frequency sweep completes (and only when DO_NONMONO_KNOP
! is set).  For every ray LS it re-solves the full (space,frequency)
! intensity field on the opacities captured during the sweep
! (KN_CHI_EXT/KN_ETA_EXT), iterating alternating-direction frequency
! sweeps: ascending frequency index gives expanding cells (GAM>=0) a
! fresh upwind at l-1, descending gives compressing cells (GAM<0) a
! fresh upwind at l+1 -- the Paper-II ordering that follows the
! information flow along each characteristic.  The converged field is
! written back to I_P_DEPTH_NEXT/I_M_DEPTH_NEXT through the same
! under-relaxed in-place write the sweep uses, so the next great
! iteration's compressing cells see a converged upwind field rather
! than one advanced a single frequency bin.
!
! SET_FG_NU_NEXT passes nu(k+1) to the solver each call so the
! compressing branch uses the correct asymmetric dLOG_NU.
!
! ==============================================================================
!
	SUBROUTINE SET_FG_NU_NEXT(NU_NEXT)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FG_NU_NEXT
	IMPLICIT NONE
	REAL(KIND=LDP) NU_NEXT
	FG_NU_NEXT=NU_NEXT
	RETURN
	END
!
	SUBROUTINE SET_FG_KNOP_REFRESH(DO_RF)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: DO_KNOP_REFRESH
	IMPLICIT NONE
	LOGICAL DO_RF
	DO_KNOP_REFRESH=DO_RF
	RETURN
	END
!
	SUBROUTINE SET_FG_KNOP_DEBUG(DO_DBG)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: PRINT_KNOP_DEBUG
	IMPLICIT NONE
	LOGICAL DO_DBG
	PRINT_KNOP_DEBUG=DO_DBG
	RETURN
	END
!
	SUBROUTINE SET_FG_FREQ_MATCH(DO_MATCH)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FG_FREQ_MATCH
	IMPLICIT NONE
	LOGICAL DO_MATCH
	FG_FREQ_MATCH=DO_MATCH
	RETURN
	END
!
	SUBROUTINE SET_FG_ORD_SWEEP(DO_ORD)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FG_ORD_SWEEP
	IMPLICIT NONE
	LOGICAL DO_ORD
	FG_ORD_SWEEP=DO_ORD
	RETURN
	END
!
	SUBROUTINE SET_FG_AITKEN(DO_AIT)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FG_ORD_AITKEN
	IMPLICIT NONE
	LOGICAL DO_AIT
	FG_ORD_AITKEN=DO_AIT
	RETURN
	END
!
	SUBROUTINE SET_FG_ORD_AA(DO_AA,AA_M)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: FG_ORD_AA,FG_ORD_AA_M
	IMPLICIT NONE
	LOGICAL DO_AA
	INTEGER AA_M
	FG_ORD_AA=DO_AA
	FG_ORD_AA_M=MAX(1,MIN(AA_M,8))
	RETURN
	END
!
	SUBROUTINE SET_FG_ORD_P2DAMP(DO_P2,N_P2,REL_P2)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY:
	1     FG_ORD_P2DAMP,FG_ORD_P2N,FG_ORD_P2REL
	IMPLICIT NONE
	LOGICAL DO_P2
	INTEGER N_P2,REL_P2
	FG_ORD_P2DAMP=DO_P2
	FG_ORD_P2N=MAX(1,MIN(N_P2,100))
	FG_ORD_P2REL=MAX(1,MIN(REL_P2,100))
	RETURN
	END
!
	SUBROUTINE SET_FG_STALL_LIM(T)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: KN_RAY_STALL
	IMPLICIT NONE
	REAL(KIND=LDP) T
	KN_RAY_STALL=MIN(MAX(T,0.5_LDP),0.999_LDP)
	RETURN
	END
!
! SET_FG_KN_OPTS: runtime wiring for the previously compile-time
! FG refresh knobs KN_RAY_TOL, KN_RAY_MAXPASS and the lagged-store
! relaxation factor INEXT_RELAX.
!
	SUBROUTINE SET_FG_KN_OPTS(TOL,MAXP,RLX)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY:
	1     KN_RAY_TOL,KN_RAY_MAXPASS,INEXT_RELAX
	IMPLICIT NONE
	REAL(KIND=LDP) TOL,RLX
	INTEGER MAXP
	KN_RAY_TOL=MAX(TOL,1.0E-12_LDP)
	KN_RAY_MAXPASS=MAX(MAXP,1)
	INEXT_RELAX=MIN(MAX(RLX,1.0E-3_LDP),1.0_LDP)
	RETURN
	END
!
	SUBROUTINE SET_FG_KN_ADAPT(DO_AD)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY: KN_ADAPT_RELAX
	IMPLICIT NONE
	LOGICAL DO_AD
	KN_ADAPT_RELAX=DO_AD
	RETURN
	END
!
! W1-W9 setters (17-Oct-2026): VADAT wiring for the
! maturity-chain-repair knobs.  Same clamp-in-setter pattern as
! SET_FG_KN_OPTS.
!
	SUBROUTINE SET_FG_OMEGA_V2(DO_V2,DO_PLAT)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY:
	1     FG_ADAPT_OMEGA_V2,FG_PLATEAU_DET
	IMPLICIT NONE
	LOGICAL DO_V2,DO_PLAT
	FG_ADAPT_OMEGA_V2=DO_V2
	FG_PLATEAU_DET=DO_PLAT
	RETURN
	END
!
	SUBROUTINE SET_FG_HOTPASS(NRAY_P,NTOT_P)
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY:
	1     KN_RAY_HOTPASSES,KN_HOT_TOTAL_PASSES
	IMPLICIT NONE
	INTEGER NRAY_P,NTOT_P
	KN_RAY_HOTPASSES=MAX(0,NRAY_P)
	KN_HOT_TOTAL_PASSES=MAX(0,NTOT_P)
	RETURN
	END
!
	SUBROUTINE SET_FG_COMMIT(JCAP,HCAP,DO_CC,NCAPC,BCON)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1, ONLY:
	1     FG_COMMIT_JUMP,FG_DIV_HEAL_CAP,FG_COMMIT_CONV,
	2     FG_COMMIT_NC,FG_BULK_CONTR
	IMPLICIT NONE
	REAL(KIND=LDP) JCAP,HCAP,NCAPC
	LOGICAL DO_CC,BCON
	FG_COMMIT_JUMP=MAX(JCAP,1.0_LDP)
	FG_DIV_HEAL_CAP=MAX(HCAP,1.0_LDP)
	FG_COMMIT_CONV=DO_CC
	FG_COMMIT_NC=MAX(NCAPC,0.0_LDP)
	FG_BULK_CONTR=BCON
	RETURN
	END
!
	SUBROUTINE FG_KNOP_REFRESH_V1()
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1
	USE CONTROL_VARIABLE_MOD, ONLY: CONV_DIAG_LVL
	IMPLICIT NONE
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER OMP_GET_MAX_THREADS
	REAL(KIND=LDP) OMP_GET_WTIME
	EXTERNAL OMP_GET_MAX_THREADS,OMP_GET_WTIME
!
	INTEGER LS,L,I,J,K,IOS,N_TH
	REAL(KIND=LDP) T0,T1_W,T2
	REAL(KIND=LDP), ALLOCATABLE :: RES_LS(:)
	INTEGER, ALLOCATABLE :: NST_LS(:),NPS_LS(:),NCP_LS(:)
	INTEGER, ALLOCATABLE :: NNG_LS(:),NDV_LS(:),NMX_LS(:)
	INTEGER, ALLOCATABLE :: NXT_LS(:)
	INTEGER*8, ALLOCATABLE :: NBD_LS(:),NTT_LS(:),NMT_LS(:),NCM_LS(:)
	INTEGER*8, ALLOCATABLE :: NP2_LS(:)
	INTEGER*8, ALLOCATABLE :: NBB_LS(:),NTB_LS(:)
	REAL(KIND=LDP), ALLOCATABLE :: LG_WRK(:)
	REAL(KIND=LDP) RHO_H(4)
!
! Per-exit-class census of the ray solves (-1 = unsolved early
! return, 0-5 as coded in KNOP_RAY_SOLVE_V1).  The residual max
! alone cannot say WHY rays stop short of KN_RAY_TOL: the bulk and
! stall exits both stop above tolerance while the population-release
! gate still counts their leftover bad cells, so the census shows
! which stop rule actually limits refresh maturity.
!
	INTEGER NX_CNT(-1:6)
	INTEGER*8 NX_BAD(-1:6)
	REAL(KIND=LDP) NX_RMX(-1:6)
!
	LUER=ERROR_LU()
	IF(.NOT. DO_KNOP_REFRESH)RETURN
	IF(FREEZE_I_NEXT)RETURN
	IF(.NOT. HAS_COMPRESSING)RETURN
	IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT))RETURN
	IF(.NOT. KN_STORE_READY)THEN
	  IF(.NOT. KN_RF_WARNED)THEN
	    WRITE(LUER,'(A,I8,A,I8)')'FG_KNOP_REFRESH: capture store'
	1      //' not armed (FREQ_IDX_ITER=',FREQ_IDX_ITER,
	2      ', NCF=',KN_NCF,') -- refresh skipped'
	    KN_RF_WARNED=.TRUE.
	  END IF
	  RETURN
	END IF
! FREQ_IDX_ITER is 0-based during the sweep (INIT skips the NEW_FREQ
! increment), so a complete sweep ends at KN_NCF-1.
	IF(FREQ_IDX_ITER .NE. KN_NCF-1)THEN
	  IF(.NOT. KN_RF_WARNED)THEN
	    WRITE(LUER,'(A,I8,A,I8)')'FG_KNOP_REFRESH: incomplete sweep'
	1      //' (FREQ_IDX_ITER=',FREQ_IDX_ITER,', NCF=',KN_NCF,
	2      ') -- refresh skipped'
	    KN_RF_WARNED=.TRUE.
	  END IF
	  RETURN
	END IF
	IF(SOLUTION_METHOD .NE. 'INTEGRAL')THEN
	  IF(.NOT. KN_RF_WARNED)THEN
	    WRITE(LUER,'(A)')'FG_KNOP_REFRESH: SOLUTION_METHOD is not'
	1      //' INTEGRAL -- refresh disabled'
	    KN_RF_WARNED=.TRUE.
	  END IF
	  RETURN
	END IF
!
! Monotone-cubic coefficients (log space, extended grid) for the
! captured per-frequency opacity/emissivity.  Same kernels as the
! in-sweep CHI_COEF/ETA_COEF build.
!
	ALLOCATE(LG_WRK(KN_ND_EXT),STAT=IOS)
	IF(IOS .NE. 0)RETURN
	DO L=1,KN_NCF
	  DO I=1,KN_ND_EXT
	    LG_WRK(I)=LOG(MAX(DBLE(KN_CHI_EXT(I,L)),1.0E-300_LDP))
	  END DO
	  CALL MON_INT_FUNS_V2(KN_CHI_COEF(1,1,L),LG_WRK,LOG_R_EXT,
	1    KN_ND_EXT)
	  DO I=1,KN_ND_EXT
	    LG_WRK(I)=LOG(MAX(DBLE(KN_ETA_EXT(I,L)),1.0E-300_LDP))
	  END DO
	  CALL MON_INT_FUNS_V2(KN_ETA_COEF(1,1,L),LG_WRK,LOG_R_EXT,
	1    KN_ND_EXT)
	END DO
	DEALLOCATE(LG_WRK)
!
	ALLOCATE(RES_LS(KN_NP),NST_LS(KN_NP),NPS_LS(KN_NP),
	1   NCP_LS(KN_NP),NNG_LS(KN_NP),NBD_LS(KN_NP),NTT_LS(KN_NP),
	2   NDV_LS(KN_NP),NMX_LS(KN_NP),NMT_LS(KN_NP),NCM_LS(KN_NP),
	3   NXT_LS(KN_NP),NP2_LS(KN_NP),NBB_LS(KN_NP),NTB_LS(KN_NP),
	4   STAT=IOS)
	IF(IOS .NE. 0)RETURN
	RES_LS=0.0_LDP
	NST_LS=0; NPS_LS=0; NCP_LS=0; NNG_LS=0; NBD_LS=0; NTT_LS=0
	NDV_LS=0; NMX_LS=1; NMT_LS=0; NCM_LS=0; NXT_LS=-1; NP2_LS=0_8
	NBB_LS=0_8; NTB_LS=0_8
	KN_AA_APP=0_8; KN_AA_VETO=0_8; KN_AA_REJ=0_8; KN_AA_SNG=0_8
	FG_RF_DCHG=0.0_LDP
	KN_HOT_SPENT=0
!
! Per-ray commit verdict/delta: each thread writes only its own LS
! slot inside KNOP_RAY_SOLVE_V1.  -1 marks rays that never reached
! the commit decision (unsolved/alloc-fail).
!
	IF(ALLOCATED(FG_RAY_VC))DEALLOCATE(FG_RAY_VC)
	IF(ALLOCATED(FG_RAY_DCH))DEALLOCATE(FG_RAY_DCH)
	ALLOCATE(FG_RAY_VC(KN_NP),FG_RAY_DCH(KN_NP),STAT=IOS)
	IF(IOS .EQ. 0)THEN
	  FG_RAY_VC=-1
	  FG_RAY_DCH=0.0_LDP
	END IF
	IF(ALLOCATED(FG_RAY_WI))DEALLOCATE(FG_RAY_WI)
	IF(ALLOCATED(FG_RAY_WL))DEALLOCATE(FG_RAY_WL)
	IF(ALLOCATED(FG_RAY_WS))DEALLOCATE(FG_RAY_WS)
	ALLOCATE(FG_RAY_WI(KN_NP),FG_RAY_WL(KN_NP),FG_RAY_WS(KN_NP),
	1   STAT=IOS)
	IF(IOS .EQ. 0)THEN
	  FG_RAY_WI=0
	  FG_RAY_WL=0
	  FG_RAY_WS=0
	END IF
!
	IF(ALLOCATED(KN_RES_PASS))DEALLOCATE(KN_RES_PASS)
	ALLOCATE(KN_RES_PASS(2*(KN_RAY_MAXPASS+KN_RAY_HOTPASSES)),
	1   STAT=IOS)
	IF(IOS .NE. 0)RETURN
	KN_RES_PASS=0.0_LDP
!
	IF(ALLOCATED(KN_HOT_CNT))DEALLOCATE(KN_HOT_CNT)
	ALLOCATE(KN_HOT_CNT(KN_NCF),STAT=IOS)
	IF(IOS .EQ. 0)KN_HOT_CNT=0
!
! One ray per thread: each thread owns column LS of the shared module
! arrays (GAM, RAY_PNT, I_*_DEPTH_NEXT), so no atomics are needed.
! KN_MAX_THREADS bounds the per-thread NI*NCF work-array memory.
!
	N_TH=OMP_GET_MAX_THREADS()
	IF(N_TH .GT. KN_MAX_THREADS)N_TH=KN_MAX_THREADS
!
! Global memory budget: KN_MAX_THREADS alone bounds memory only when
! NI*NCF is fixed; scale threads down so N_TH*(~80 B per cell) stays
! under KN_RF_TOT_MB.  MAXVAL(NI_RAY) is the largest per-ray NI.
!
	T2=80.0_LDP*MAXVAL(NI_RAY(1:KN_NP))*KN_NCF/1.0E+6_LDP
	IF(T2 .GT. 0.0_LDP)THEN
	  N_TH=MIN(N_TH,MAX(1,INT(KN_RF_TOT_MB/T2)))
	END IF
	T0=OMP_GET_WTIME()
!$OMP PARALLEL DO NUM_THREADS(N_TH) SCHEDULE(DYNAMIC)
!$OMP& DEFAULT(SHARED) PRIVATE(LS)
	DO LS=1,KN_NP
	  CALL KNOP_RAY_SOLVE_V1(LS,RES_LS(LS),NST_LS(LS),
	1    NPS_LS(LS),NCP_LS(LS),NNG_LS(LS),NBD_LS(LS),NTT_LS(LS),
	2    NDV_LS(LS),NMT_LS(LS),NCM_LS(LS),NMX_LS(LS),NXT_LS(LS),
	3    NP2_LS(LS),NBB_LS(LS),NTB_LS(LS))
	END DO
!$OMP END PARALLEL DO
	T1_W=OMP_GET_WTIME()-T0
!
! The per-ray residual replaces the one-bin-lag GS metric: a converged
! refresh means the compressing-branch field the next sweep reads is
! a fixed point of the wavelength coupling for the current source
! function, so the SE update may proceed immediately.
!
! W2 (17-Oct-2026): diverged rays report RESID=1.0 so a lone
! runaway cannot pin the max-norm backstop -- one diverged ray out
! of hundreds used to hold GS_MAXCHG=1.0 > GS_UNSAFE_TOL on every
! iteration.  They are counted in FG_GS_NDIV and share-gated in
! SOLVE_FOR_POPS instead; the backstop still sees every
! non-diverged residual.
!
	IF(COUNT(NDV_LS .EQ. 0) .GT. 0)THEN
	  KN_RF_RESID=MAXVAL(RES_LS,MASK=NDV_LS .EQ. 0)
	ELSE
	  KN_RF_RESID=MAXVAL(RES_LS)
	END IF
!
! Adaptive under-relaxation of the lagged store (KN_ADAPT_RELAX,
! default T): the residual maxima of successive refresh calls are the
! dominant-mode contraction factors rho_est of the outer
! sweep+refresh fixed-point map.  omega = 1/(1+rho_est) damps an
! expansive or sign-flipping mode (rho_est >= 1) to <= 1/2 and lets a
! contracting map recover toward the undamped write.  INEXT_RELAX is
! smoothed over calls so one noisy residual cannot swing it; the
! [0.05,1] clamp keeps the update causal (bounded damping) and
! prevents a collapsed residual from disabling relaxation forever.
! VADAT INEXT_RELAX then acts as the seed value.  Oct-2026.
!
! W5 (17-Oct-2026, [FG_ADAPT_OMEGA_V2], default T): the old
! cross-iteration ratio miskeys omega -- the refresh residual drifts
! with the populations' source update, so a converging field can
! read rho_est~1 whether the map oscillates (wants omega~0.5) or
! contracts marginally (wants omega~1).  KN_RES_PASS holds the
! per-pass max-over-rays residual sequence of THIS refresh, so the
! median of its last ~4 pass ratios measures the same-map
! contraction factor directly: oscillation (rho>1) still drives
! omega down, smooth contraction (rho~0.9) lets it rise toward the
! unrelaxed commit W4 made safe.  Clamp [0.05,20] as in the old
! rule; F restores the cross-iteration estimate.
!
	IF(KN_ADAPT_RELAX .AND. KN_RF_RESID .GT. 0.0_LDP)THEN
	  IF(FG_ADAPT_OMEGA_V2 .AND. ALLOCATED(KN_RES_PASS))THEN
	    J=SIZE(KN_RES_PASS)
	    DO WHILE(J .GT. 1)
	      IF(KN_RES_PASS(J) .GT. 0.0_LDP)EXIT
	      J=J-1
	    END DO
	    IF(J .GE. 3 .AND. KN_RES_PASS(J-1) .GT. 0.0_LDP)THEN
	      L=0
	      DO I=MAX(2,J-3),J
	        IF(KN_RES_PASS(I-1) .GT. 0.0_LDP)THEN
	          L=L+1
	          RHO_H(L)=KN_RES_PASS(I)/KN_RES_PASS(I-1)
	        END IF
	      END DO
	      IF(L .GT. 0)THEN
	        DO I=1,L-1
	          DO K=I+1,L
	            IF(RHO_H(K) .LT. RHO_H(I))THEN
	              T2=RHO_H(I)
	              RHO_H(I)=RHO_H(K)
	              RHO_H(K)=T2
	            END IF
	          END DO
	        END DO
	        T2=RHO_H((L+1)/2)
	        T2=MIN(MAX(T2,0.05_LDP),20.0_LDP)
	        IF(CONV_DIAG_LVL .GE. 2)WRITE(LUER,'(A,ES10.3,A,I4,
	1        A,ES10.3)')
	2        '    omega est: rho=',T2,' window passes=',L,
	3        ' applied omega=',0.5_LDP*INEXT_RELAX+
	4        0.5_LDP*MIN(1.0_LDP,MAX(0.05_LDP,1.0_LDP/(1.0_LDP+T2)))
	        INEXT_RELAX=0.5_LDP*INEXT_RELAX+
	1          0.5_LDP*MIN(1.0_LDP,MAX(0.05_LDP,1.0_LDP/(1.0_LDP+T2)))
	      END IF
	    END IF
	  ELSE IF(KN_RF_RESID_PRV .GT. 0.0_LDP)THEN
	    T2=KN_RF_RESID/KN_RF_RESID_PRV
	    T2=MIN(MAX(T2,0.0_LDP),50.0_LDP)
	    INEXT_RELAX=0.5_LDP*INEXT_RELAX+
	1      0.5_LDP*MIN(1.0_LDP,MAX(0.05_LDP,1.0_LDP/(1.0_LDP+T2)))
	  END IF
	  KN_RF_RESID_PRV=KN_RF_RESID
	END IF
!
! W4 (17-Oct-2026): keep both metric semantics.  The sweep churn
! (fresh solve vs stored field -- what the populations saw) is
! snapshotted into FG_SW_* before FG_GS_* takes the refresh
! fixed-point values; SOLVE_FOR_POPS now reads both.  MOM keeps
! the identical split in MOM_SW_*.
!
	FG_SW_MAXCHG=FG_GS_MAXCHG
	FG_SW_NCONV=FG_GS_NCONV
	FG_GS_MAXCHG=KN_RF_RESID
	FG_GS_NCONV=SUM(NST_LS)
	FG_RF_CAPPED=SUM(NCP_LS)+SUM(NNG_LS)+SUM(NDV_LS)
	FG_INEXT_CAPPED=FG_INEXT_CAPPED+FG_RF_CAPPED
	FG_GS_NDIV=SUM(NDV_LS)
	FG_GS_NBAD=SUM(NBD_LS)
	FG_GS_NTOT=SUM(NTT_LS)
	FG_GS_NBAD_BLK=SUM(NBB_LS)
	FG_GS_NTOT_BLK=SUM(NTB_LS)
	FG_GS_NSTAL=NX_BAD(3)+NX_BAD(4)
	FG_GS_NPLAT=NX_BAD(6)
	WRITE(LUER,'(/,A)')' --------------------------------- KNOP FIELD ----------------------------------'
	WRITE(LUER,'(A)')'  KNOP I_NEXT refresh (lagged coupling):'
	WRITE(LUER,'(A,ES12.4,A,ES9.1,A)')
	1  '    residual max=',KN_RF_RESID,'  (tol',KN_RAY_TOL,')'
	WRITE(LUER,'(A,I6)')'    passes<=',MAXVAL(NPS_LS)
	WRITE(LUER,'(A,I8)')'    diverged rays=',SUM(NDV_LS)
	WRITE(LUER,'(A,I8)')'    clamp events=',SUM(NNG_LS)
	WRITE(LUER,'(A,I14,A,I14,A,ES9.2)')
	1  '    unconverged cells=',FG_GS_NBAD,' of',FG_GS_NTOT,
	2  '   frac=',DBLE(FG_GS_NBAD)/MAX(1_8,FG_GS_NTOT)
	WRITE(LUER,'(A,I14)')'    store cells>50%=',FG_GS_NCONV
	WRITE(LUER,'(A,I14,A,I14)')
	1  '    best-pass bad cells=',FG_GS_NBAD_BLK,' of',
	2  FG_GS_NTOT_BLK
	WRITE(LUER,'(A,I6,A,I6,A,I6,A,I6)')
	1  '    healed slots=',FG_DIV_HEAL,
	2  '   commit rejects=',FG_COMMIT_REJ,
	3  '   nc rejects=',FG_NC_REJ,
	4  '   nc capped=',FG_NC_CAP
	IF(KN_HOT_SPENT .GT. 0)WRITE(LUER,'(A,I8,A,I6)')
	1  '    hot-pass budget spent=',INT(MIN(KN_HOT_SPENT,2**30_8)),
	2  '  cap=',KN_HOT_TOTAL_PASSES
	WRITE(LUER,'(A,F9.1,A)')'    t=',T1_W,' s'
	WRITE(LUER,'(A,ES12.4)')'    retained-field delta: max=',
	1  FG_RF_DCHG
	IF(KN_ADAPT_RELAX)WRITE(LUER,'(A,ES10.3)')
	1  '    adaptive INEXT_RELAX: omega=',INEXT_RELAX
!
! Ray exit census: which stop rule ended each field solve (codes in
! KNOP_RAY_SOLVE_V1).  Stall and cap exits leave above-tolerance
! cells that the release gate still reads as immature, so their
! counts and bad-cell share decide whether the stall criterion or
! the pass cap is the refresh bottleneck.
!
	NX_CNT=0; NX_BAD=0; NX_RMX=0.0_LDP
	DO LS=1,KN_NP
	  NX_CNT(NXT_LS(LS))=NX_CNT(NXT_LS(LS))+1
	  NX_BAD(NXT_LS(LS))=NX_BAD(NXT_LS(LS))+NBD_LS(LS)
	  NX_RMX(NXT_LS(LS))=MAX(NX_RMX(NXT_LS(LS)),RES_LS(LS))
	END DO
	WRITE(LUER,'(A,I5,6(A,I5))')
	1  '    ray exits: exact=',NX_CNT(0),' tol=',NX_CNT(1),
	2  ' bulk=',NX_CNT(2),' stall=',NX_CNT(3),' cap=',NX_CNT(4),
	3  ' div=',NX_CNT(5),' plat=',NX_CNT(6)
	IF(NX_CNT(3)+NX_CNT(4)+NX_CNT(6) .GT. 0)THEN
	  WRITE(LUER,'(A,ES10.3,A,I14)')
	1    '    stall/cap rays: maxres=',
	2    MAX(NX_RMX(3),MAX(NX_RMX(4),NX_RMX(6))),
	3    ' bad cells=',NX_BAD(3)+NX_BAD(4)
	END IF
	IF(NX_CNT(-1) .GT. 0)THEN
	  WRITE(LUER,'(A,I5)')'    ray exits: unsolved=',NX_CNT(-1)
	END IF
	IF(FG_ORD_AA)THEN
	  WRITE(LUER,'(4(A,I8))')
	1    '    AA: applied=',KN_AA_APP,' vetoed=',KN_AA_VETO,
	2    ' rejected=',KN_AA_REJ,' solve-fail=',KN_AA_SNG
	END IF
	IF(FG_ORD_P2DAMP)THEN
	  WRITE(LUER,'(A,I14)')
	1    '    p2damped writes=',SUM(NP2_LS)
	END IF
!
! Temporary DO_NONMONO_KNOP diagnostic table.  The formal refresh can
! appear converged in aggregate while one ray remains tied to a long
! compressing-frequency stencil.  Report that ray and how often the
! frequency matcher selected more than the adjacent upwind bin.
!
	LS=MAXLOC(RES_LS,DIM=1)
	WRITE(LUER,'(A,I6,A,I6)')
	1  ' FG-KNOP worst ray=',LS,' ni=',NI_RAY(LS)
	WRITE(LUER,'(A,ES11.3,A,I4)')
	1  '    residual=',RES_LS(LS),' passes=',NPS_LS(LS)
	IF(FG_ORD_P2DAMP)WRITE(LUER,'(A,I14)')
	1  '    p2damped=',NP2_LS(LS)
	WRITE(LUER,'(A,I14,A,I14)')
	1  '    matched/compress=',SUM(NMT_LS),'/',SUM(NCM_LS)
	WRITE(LUER,'(A,I6)')'    max offset=',MAXVAL(NMX_LS)
!
! Commit census + worst-ray verdict (17-Oct-2026): the exit census
! alone cannot show whether an above-tolerance ray's field reached
! the store.  Verdict codes: 0 committed, 1 jump reject, 2
! non-converged reject, 3 capped commit (FG_COMMIT_NC), 4 diverged,
! -1 unsolved.  The worst ray's own exit class, verdict and
! retained delta identify which acceptance channel carries the
! largest residual.
!
	IF(ALLOCATED(FG_RAY_VC))THEN
	  WRITE(LUER,'(6(A,I5),A,I5,A,I5,A)')
	1    '    commit verdicts: committed=',COUNT(FG_RAY_VC.EQ.0),
	2    '  jump-rej=',COUNT(FG_RAY_VC.EQ.1),
	3    '  nc-rej=',COUNT(FG_RAY_VC.EQ.2),
	4    '  capped=',COUNT(FG_RAY_VC.EQ.3),
	5    '  div=',COUNT(FG_RAY_VC.EQ.4),
	6    '  unsolved=',COUNT(FG_RAY_VC.EQ.-1),
	7    '  worst=(ixit',NXT_LS(LS),',v',FG_RAY_VC(LS),')'
	  WRITE(LUER,'(A,ES11.3,A,ES11.3,A,ES11.3)')
	1    '    worst-ray delta: commit=',FG_RAY_DCH(LS),
	2    '  max=',MAXVAL(FG_RAY_DCH),
	3    '  committed-max=',
	4    MAX(0.0_LDP,MAXVAL(FG_RAY_DCH,MASK=FG_RAY_VC.EQ.0))
	END IF
	IF(ALLOCATED(FG_RAY_WI))THEN
	  WRITE(LUER,'(A,I5,A,I5,A,I4,A,I4)')
	1    '    worst cell: i=',FG_RAY_WI(LS),' l=',FG_RAY_WL(LS),
	2    ' same-cell passes=',FG_RAY_WS(LS),
	3    ' max same=',MAXVAL(FG_RAY_WS)
	END IF
!
! Per-pass convergence history (max over rays).  I_M and I_P sweeps
! share the array; pass counts can differ between the two fields, so
! the trailing zeros are trimmed on print.  The pass list and the
! individual hot-frequency runs are debug-only; the one-line count
! stays in OUTGEN.
!
	IF(ALLOCATED(KN_RES_PASS))THEN
	  J=SIZE(KN_RES_PASS)
	  DO WHILE(J .GT. 0)
	    IF(KN_RES_PASS(J) .GT. 0.0_LDP)EXIT
	    J=J-1
	  END DO
	  IF(J .GT. 0 .AND. PRINT_KNOP_DEBUG)THEN
	    WRITE(LUER,'(/,A)')' ------------------------------- RAY CONVERGENCE -------------------------------'
	    WRITE(LUER,'(A)')
	1   '    ray convergence (per-pass max residual):'
	    DO I=1,J,8
	      WRITE(LUER,'(8(1X,ES10.3))')
	1        KN_RES_PASS(I:MIN(J,I+7))
	    END DO
	  END IF
	  DEALLOCATE(KN_RES_PASS)
	END IF
!
! Hot-frequency runs: contiguous stretches of KN_HOT_CNT>0 are the
! block-GS candidates (cmfgen_obsidian/audits/BLOCK_GS_DESIGN.md).  A contiguous run
! means the marginal |Q|/(chi+|Q|) links form a solvable block;
! scattered hits mean the residual floor is distributed and the
! block structure buys little.
!
	IF(ALLOCATED(KN_HOT_CNT))THEN
	  L=COUNT(KN_HOT_CNT(1:KN_NCF) .GT. 0)
	  J=0
	  I=1
	  DO WHILE(I .LE. KN_NCF)
	    IF(KN_HOT_CNT(I) .GT. 0)THEN
	      J=J+1
	      LS=I
	      DO WHILE(I .LT. KN_NCF)
	        IF(KN_HOT_CNT(I+1) .LE. 0)EXIT
	        I=I+1
	      END DO
	      IF(J .LE. 12 .AND. PRINT_KNOP_DEBUG)THEN
	        IF(J .EQ. 1)WRITE(LUER,'(/,A)')
	1        ' ------------------------------- HOT FREQUENCIES -------------------------------'
	        WRITE(LUER,'(A,I3,A,I8,A,I8,A,ES10.3,A,ES10.3,A,I5)')
	1        '    hot run',J,': L=',LS,'-',I,
	2        '  nu=',KN_FREQ(LS),'-',KN_FREQ(I),
	3        '  max rays=',MAXVAL(KN_HOT_CNT(LS:I))
	      END IF
	    END IF
	    I=I+1
	  END DO
	  WRITE(LUER,'(A,I8,A,I5,A)')'    hot frequencies:',L,
	1    ' in',J,' contiguous runs (block-GS candidates)'
	  DEALLOCATE(KN_HOT_CNT)
	END IF
	DEALLOCATE(RES_LS,NST_LS,NPS_LS,NCP_LS,NNG_LS,NBD_LS,NTT_LS,
	1   NDV_LS,NMX_LS,NMT_LS,NCM_LS,NXT_LS,NP2_LS,
	2   NBB_LS,NTB_LS)
	RETURN
	END
!
! ==============================================================================
!
! KNOP_RAY_SOLVE_V1: fixed-point solve of the coupled wavelength
! transfer problem for a single ray, using the sweep-captured
! opacities and the same short-characteristics discretisation as the
! in-sweep INTEGRAL path.  The field arrays are local to this routine
! (per-thread under the OMP loop in FG_KNOP_REFRESH_V1).
!
! ==============================================================================
!
	SUBROUTINE KNOP_RAY_SOLVE_V1(LS,RESID,NSTORE,NPASS,NCAP,NNEG,
	1   NBAD,NTOT,NDIV,NMATCH,NCOMP,NMAX,NXIT,NP2,NBBLK,NTBLK)
	USE SET_KIND_MODULE
	USE FG_J_CMF_NONMONO_MOD_V1
	IMPLICIT NONE
!
	INTEGER LS,NSTORE,NPASS,NCAP,NNEG,NDIV,NMAX,NXIT
	INTEGER*8 NBAD,NTOT,NMATCH,NCOMP,NP2,NBBLK,NTBLK
	REAL(KIND=LDP) RESID
!
! NBBLK/NTBLK (W3, 17-Oct-2026): best-pass bad-cell census summed
! over the two field solves.  For stall/cap exits the last-pass
! NBAD still counts every above-tolerance cell, so the gate reads
! MIN(GS_NBAD, GS_NBAD_BLK) -- the best pass is the closest the
! fixed-point map actually reached.
!
! W_LO/W_HI/NHOT (W7, 17-Oct-2026): hot-pass extension state.  When
! a stall or cap exit fires with residual left only inside a
! contiguous hot-frequency run, the sweep continues restricted to
! that window (+/-2) for up to KN_RAY_HOTPASSES extra passes;
! KN_HOT_SPENT bounds the refresh-wide spend.
!
	INTEGER W_LO,W_HI,NHOT
!
! SGN_ACC/SGN_PRV/NOSC_FL (W6, 17-Oct-2026): the KN_SGN_SUM sign-flip
! detector ported to the alternating SWEEP_FLD path -- MARCH adds the
! signed per-cell update into the host-scope SGN_ACC.
!
	REAL(KIND=LDP) SGN_ACC,SGN_PRV
	INTEGER NOSC_FL
!
	INTEGER, PARAMETER :: NC_PNT_SRCE_L=2
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER KN_DUMP_N
	INTEGER NI,NCF
	INTEGER I,K,KC,L,JPT,IOS
	INTEGER, ALLOCATABLE :: IKC1(:),IKC2(:),JPR1(:),JPR2(:),IGRD(:)
	INTEGER NGRD
	REAL(KIND=LDP), ALLOCATABLE :: WINX(:),IM_NI(:),T1V(:)
	REAL(KIND=LDP), ALLOCATABLE :: FLD(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: FLD_P(:,:)
!
! RES_FREQ(L): per-frequency max overwrite residual of the last
! completed pass, accumulated inside MARCH.  Stage-1 instrumentation
! for the block-GS redesign (cmfgen_obsidian/audits/BLOCK_GS_DESIGN.md): the global
! max residual cannot show whether the marginal cells sit on a
! contiguous frequency run (a block) or are scattered -- the
! near-unity |Q|/(chi+|Q|) links live at frequencies whose line
! cores overlap the compression ramp.  Folded into the module-level
! KN_HOT_CNT histogram once per completed field solve.
!
	REAL(KIND=LDP), ALLOCATABLE :: RES_FREQ(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_L(:),ETA_L(:),DCHR_L(:),
	1   Q_L(:),SR_L(:),SP_L(:),S_L(:),DS_L(:),DTAU_L(:),
	2   A0L(:),A1L(:),A2L(:),A3L(:),A4L(:),
	3   EEL(:),E0L(:),E1L(:),E2L(:),E3L(:)
!
! Matched-frequency work arrays (DO_NONMONO_FREQ_MATCH): per-ray
! drift in ln(nu) units DRF_L = |GAM|*ds, refilled per frequency
! into NFM_L (slot offset) and DLNM_L (log span).  MATCH_FREQ
! fills both; defaults keep the +1 slot / KN_DLN_N(L) span.
!
	INTEGER, ALLOCATABLE :: NFM_L(:)
	REAL(KIND=LDP), ALLOCATABLE :: DLNM_L(:),DRF_L(:)
!
! Pass-invariant per-frequency coefficient pack (see the build block
! below): the BUILD_FREQ outputs depend only on the captured
! opacities and static ray geometry, never on the iterated field, so
! they are built once per ray and replayed each pass for both
! direction fields.  8*NI*NCF*8 B plus the two NCF scalars.
!
	REAL(KIND=LDP), ALLOCATABLE :: PK_DTAU(:,:),PK_A0(:,:),
	1   PK_A1(:,:),PK_A2(:,:),PK_A3(:,:),PK_A4(:,:),
	2   PK_SR(:,:),PK_Q(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: PK_ETANI(:),PK_CHINI(:)
!
! SPW: 4-deep ring of the modified source S' used by the ordered
! march (SWEEP_FLD_ORD).  A march step touches four consecutive rows
! only, so MOD(row,4) slots never collide.
!
	REAL(KIND=LDP), ALLOCATABLE :: SPW(:,:)
!
	REAL(KIND=LDP) RES_M,RES_P,DLN,T1,T2,DZ,T1_DTAU,I_CORE
	REAL(KIND=LDP) DCH_M,DCH_P
	REAL(KIND=LDP) VNEW,VOLD,INEXT_V
	REAL(KIND=LDP) OMEGA_RAY,RES_PREV(2),RATIO_OM
!
! Worst-cell tracking (17-Oct-2026): MARCH records the (depth,freq)
! cell holding the running max residual; SWEEP_FLD counts how many
! consecutive passes keep the same worst cell.  A persistent cell
! identifies a single unstable mode; a moving max identifies
! different cells resetting the residual each pass.
!
	INTEGER WRC_I,WRC_L,WRC_PI,WRC_PL,WRC_SAME
	INTEGER WRC_IM,WRC_LM,WRC_SM
	REAL(KIND=LDP) DCELL
	INTEGER*8 NBAD_CUR,NTOT_CUR,NB_M,NB_P,NT_M,NT_P
	INTEGER*8 NP2_CNT
	INTEGER IDIR
	INTEGER NP_M,NP_P
! Exit census of the two field solves: 0 exact (no compressing
! cells), 1 residual tol, 2 bulk bad-fraction, 3 same-direction
! stall, 4 pass cap, 5 diverged.  NXIT reports the field that sets
! the ray residual, so the summary shows which stop rule limits the
! refresh -- the stall and cap exits leave cells above KN_RAY_TOL
! that the population-release gate still counts as immature.
!
	INTEGER IXIT_M,IXIT_P
	INTEGER*8 NB_BM,NT_BM,NB_BP,NT_BP
	LOGICAL HAS_CP,CLEAN_EXIT,DIVERGED,HAVE_PACK,ORD_SWEEP
	LOGICAL L_HOT_OVER
	REAL(KIND=LDP) JUMP
!
	RESID=0.0_LDP
	NSTORE=0; NPASS=0; NCAP=0; NNEG=0; NDIV=0; NMAX=1
	NBAD=0; NTOT=0; NB_M=0; NB_P=0; NT_M=0; NT_P=0
	NMATCH=0; NCOMP=0
	NXIT=-1; IXIT_M=-1; IXIT_P=-1
	RES_M=0.0_LDP; RES_P=0.0_LDP; NP_M=0; NP_P=0
	NBBLK=0_8; NTBLK=0_8; NB_BM=0_8; NT_BM=0_8
	NB_BP=0_8; NT_BP=0_8
	NP2=0_8; NP2_CNT=0_8
	DIVERGED=.FALSE.
	DCH_M=0.0_LDP; DCH_P=0.0_LDP
	W_LO=0; W_HI=0; NHOT=0
	SGN_ACC=0.0_LDP; SGN_PRV=0.0_LDP; NOSC_FL=0
	LUER=ERROR_LU()
	KN_DUMP_N=0
	NI=NI_RAY(LS)
	NCF=KN_NCF
	IF(NI .LT. 2 .OR. NCF .LT. 2)RETURN
!
	ALLOCATE(CHI_L(NI),ETA_L(NI),DCHR_L(NI),Q_L(NI),SR_L(NI),
	1   SP_L(NI),S_L(NI),DS_L(NI),DTAU_L(NI),
	2   A0L(NI),A1L(NI),A2L(NI),A3L(NI),A4L(NI),
	3   EEL(NI),E0L(NI),E1L(NI),E2L(NI),E3L(NI),
	4   IKC1(NI),IKC2(NI),JPR1(NI),JPR2(NI),WINX(NI),T1V(NI),
	5   IGRD(NI),IM_NI(NCF),RES_FREQ(NCF),
	6   NFM_L(NI),DLNM_L(NI),DRF_L(NI),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,'(A,I5,A,I6,A,I8)')
	1   '  FG-KNOP refresh alloc failed: ray LS=',LS,
	2   ' NI=',NI,' bytes~',27_8*NI+3_8*NCF
	  RESID=1.0_LDP
	  RETURN
	END IF
	RES_FREQ=0.0_LDP
!
! Per-ray-point compressing lookup, identical logic to the in-sweep
! I_NEXT read: ray point I sits between depth-grid points K and K+1
! (K=RAY_PNT-ND_ADD); only depths flagged compressing carry an I_NEXT
! slot.  JPR* are the ray points J_PNT maps those depths to, so the
! upwind interpolant reads the ray's own field at exactly the points
! the in-sweep code reads the sparse store.
!
	HAS_CP=.FALSE.
	NGRD=0
	DO I=1,NI
	  IKC1(I)=0; IKC2(I)=0; JPR1(I)=0; JPR2(I)=0
	  WINX(I)=0.0_LDP
	  K=RAY_PNT(I,LS)-ND_ADD
	  IF(K .GE. 1 .AND. K .LE. KN_ND)THEN
	    IF(COMP_DEPTH_IDX(K) .GE. 1)THEN
	      JPT=J_PNT(K,LS)
	      IF(JPT .GE. 1 .AND. JPT .LE. NI)THEN
	        IKC1(I)=1; JPR1(I)=JPT
	      END IF
	    END IF
	  END IF
	  IF(K+1 .GE. 1 .AND. K+1 .LE. KN_ND)THEN
	    IF(COMP_DEPTH_IDX(K+1) .GE. 1)THEN
	      JPT=J_PNT(K+1,LS)
	      IF(JPT .GE. 1 .AND. JPT .LE. NI)THEN
	        IKC2(I)=1; JPR2(I)=JPT
	      END IF
	    END IF
	  END IF
	  IF(IKC1(I) .EQ. 1 .AND. IKC2(I) .EQ. 1)THEN
	    IF(KN_R(K+1) .NE. KN_R(K))THEN
	      WINX(I)=(R_RAY(I,LS)-KN_R(K))/(KN_R(K+1)-KN_R(K))
	      WINX(I)=MAX(0.0_LDP,MIN(1.0_LDP,WINX(I)))
	    END IF
	  END IF
	  IF(R_RAY(I,LS) .EQ. R_EXT(RAY_PNT(I,LS)))THEN
	    NGRD=NGRD+1; IGRD(NGRD)=I
	  END IF
	  IF(GAM(I,LS) .LT. 0.0_LDP .AND.
	1     (IKC1(I) .EQ. 1 .OR. IKC2(I) .EQ. 1))HAS_CP=.TRUE.
!
! Per-point comoving frequency drift across the local cell, used by
! MATCH_FREQ to rematch the compressing upwind slot (offset 1 =
! adjacent bin).  ds is the symmetric Z cell, same convention as
! the in-sweep NM_RAY fill.
!
	  IF(I .EQ. 1)THEN
	    DRF_L(I)=ABS(Z(1,LS)-Z(2,LS))
	  ELSE IF(I .EQ. NI)THEN
	    DRF_L(I)=ABS(Z(NI-1,LS)-Z(NI,LS))
	  ELSE
	    DRF_L(I)=0.5_LDP*(ABS(Z(I-1,LS)-Z(I,LS))+
	1       ABS(Z(I,LS)-Z(I+1,LS)))
	  END IF
	  DRF_L(I)=ABS(GAM(I,LS))*DRF_L(I)
	END DO
!
! Build the pass-invariant coefficient pack once per ray.  Every
! quantity BUILD_FREQ produces (cubic opacity evals, chi_hat, Q, SR,
! DCHR, DTAU, the E/A short-characteristics coefficients) depends only
! on the captured opacities and static ray geometry -- never on the
! iterated field -- so replaying it each pass is bit-identical to the
! previous per-pass rebuild while removing ~70% of the per-pass cost.
! The pack is shared by the I_M and I_P solves, whose per-frequency
! coefficients are identical.
!
! Fallback: when the pack does not fit within KN_RAY_PACK_MB per
! thread or the allocation fails, SWEEP_FLD keeps calling BUILD_FREQ
! inside the pass loop (previous behaviour, unchanged).
!
	HAVE_PACK=.FALSE.
	IF(64.0_LDP*NI*NCF .LE. 1.0E6_LDP*KN_RAY_PACK_MB)THEN
	  ALLOCATE(PK_DTAU(NI,NCF),PK_A0(NI,NCF),PK_A1(NI,NCF),
	1     PK_A2(NI,NCF),PK_A3(NI,NCF),PK_A4(NI,NCF),
	2     PK_SR(NI,NCF),PK_Q(NI,NCF),
	3     PK_ETANI(NCF),PK_CHINI(NCF),STAT=IOS)
	  HAVE_PACK=(IOS .EQ. 0)
	  IF(.NOT. HAVE_PACK)THEN
	    IF(ALLOCATED(PK_DTAU))DEALLOCATE(PK_DTAU)
	    IF(ALLOCATED(PK_A0))DEALLOCATE(PK_A0)
	    IF(ALLOCATED(PK_A1))DEALLOCATE(PK_A1)
	    IF(ALLOCATED(PK_A2))DEALLOCATE(PK_A2)
	    IF(ALLOCATED(PK_A3))DEALLOCATE(PK_A3)
	    IF(ALLOCATED(PK_A4))DEALLOCATE(PK_A4)
	    IF(ALLOCATED(PK_SR))DEALLOCATE(PK_SR)
	    IF(ALLOCATED(PK_Q))DEALLOCATE(PK_Q)
	    IF(ALLOCATED(PK_ETANI))DEALLOCATE(PK_ETANI)
	    IF(ALLOCATED(PK_CHINI))DEALLOCATE(PK_CHINI)
	  END IF
	END IF
	IF(HAVE_PACK)THEN
	  DO L=1,NCF
	    CALL BUILD_FREQ(L)
	    NCOMP=NCOMP+INT(COUNT(GAM(1:NI,LS) .LT. 0.0_LDP),8)
	    NMATCH=NMATCH+INT(COUNT(NFM_L(1:NI) .GT. 1),8)
	    NMAX=MAX(NMAX,MAXVAL(NFM_L(1:NI)))
	    PK_DTAU(1:NI,L)=DTAU_L(1:NI)
	    PK_A0(1:NI,L)=A0L(1:NI)
	    PK_A1(1:NI,L)=A1L(1:NI)
	    PK_A2(1:NI,L)=A2L(1:NI)
	    PK_A3(1:NI,L)=A3L(1:NI)
	    PK_A4(1:NI,L)=A4L(1:NI)
	    PK_SR(1:NI,L)=SR_L(1:NI)
	    PK_Q(1:NI,L)=Q_L(1:NI)
	    PK_ETANI(L)=ETA_L(NI)
	    PK_CHINI(L)=CHI_L(NI)
	  END DO
	END IF
!
! The ordered march (KN_RAY_ORDERED) indexes the pack per (row,
! frequency), so it needs the pack and at least three rows for the
! slope endpoints; otherwise the alternating passes run unchanged.
!
	ORD_SWEEP=FG_ORD_SWEEP .AND. HAVE_PACK .AND. NI .GE. 3
	IF(ORD_SWEEP)THEN
	  ALLOCATE(SPW(NCF,4),STAT=IOS)
	  IF(IOS .NE. 0)ORD_SWEEP=.FALSE.
	END IF
!
! ---------------- I_M (inward) field ----------------
!
	ALLOCATE(FLD(NI,NCF),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,'(A,I5)')'  FG-KNOP refresh alloc failed: FLD, LS=',LS
	  RESID=1.0_LDP
	  GOTO 900
	END IF
	FLD(1:NI,1:NCF)=0.0_LDP
	CALL SEED_FLD(FLD,I_M_DEPTH_NEXT)
	W_LO=0; W_HI=0; NHOT=0
	IF(ORD_SWEEP)THEN
	  CALL SWEEP_FLD_ORD(FLD,.TRUE.,RES_M,NP_M,NB_M,NT_M,IXIT_M,
	1     NB_BM,NT_BM)
	ELSE
	  CALL SWEEP_FLD(FLD,.TRUE.,RES_M,NP_M,NB_M,NT_M,IXIT_M,
	1     NB_BM,NT_BM)
	END IF
!
! Worst-cell record of the I_M sweep; the I_P sweep below overwrites
! the shared WRC_* host variables, so keep the direction that ends
! with the larger residual.
!
	WRC_IM=WRC_I; WRC_LM=WRC_L; WRC_SM=WRC_SAME
	IF(DIVERGED)THEN
!
!	  Diverged I_M: the runaway field must not reach the store, and
!	  I_P's inner boundary reads IM_NI from this field, so the I_P
!	  solve is meaningless too.  Both stores keep their pre-refresh
!	  (in-sweep) state for this ray.
!
	  NDIV=1
	  DEALLOCATE(FLD)
	  GOTO 800
	END IF
	IM_NI(1:NCF)=FLD(NI,1:NCF)
!
! Fold this field's last-pass residuals into the hot-frequency
! histogram before the I_P solve overwrites RES_FREQ.
!
	IF(ALLOCATED(KN_HOT_CNT))THEN
!$OMP CRITICAL (KN_HOT_UPD)
	  DO L=1,NCF
	    IF(RES_FREQ(L) .GT. KN_RAY_TOL)KN_HOT_CNT(L)=KN_HOT_CNT(L)+1
	  END DO
!$OMP END CRITICAL (KN_HOT_UPD)
	END IF
!
! ---------------- I_P (outward) field ----------------
!
! The I_M store is committed only after the I_P solve also succeeds:
! writing I_M first left the stores at different refresh generations
! on an I_P divergence, so the next great iteration would read a
! mixed-generation Gauss-Seidel state.  FLD keeps the I_M result
! buffered; FLD_P carries the I_P solve.
!
	ALLOCATE(FLD_P(NI,NCF),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,'(A,I5)')'  FG-KNOP refresh alloc failed: FLD_P, LS=',LS
	  RESID=1.0_LDP
	  DEALLOCATE(FLD)
	  GOTO 900
	END IF
	FLD_P(1:NI,1:NCF)=0.0_LDP
	CALL SEED_FLD(FLD_P,I_P_DEPTH_NEXT)
	W_LO=0; W_HI=0; NHOT=0
	IF(ORD_SWEEP)THEN
	  CALL SWEEP_FLD_ORD(FLD_P,.FALSE.,RES_P,NP_P,NB_P,NT_P,IXIT_P,
	1     NB_BP,NT_BP)
	ELSE
	  CALL SWEEP_FLD(FLD_P,.FALSE.,RES_P,NP_P,NB_P,NT_P,IXIT_P,
	1     NB_BP,NT_BP)
	END IF
	IF(DIVERGED)THEN
	  NDIV=1
	ELSE
	  IF(ALLOCATED(KN_HOT_CNT))THEN
!$OMP CRITICAL (KN_HOT_UPD)
	    DO L=1,NCF
	      IF(RES_FREQ(L) .GT. KN_RAY_TOL)KN_HOT_CNT(L)=
	1         KN_HOT_CNT(L)+1
	    END DO
!$OMP END CRITICAL (KN_HOT_UPD)
	  END IF
!
! W9 (17-Oct-2026, [FG_COMMIT_JUMP], default 50.0): whole-ray
! commit guard.  A runaway pass that stayed below the divergence
! threshold would once commit a near-blowup into the store (audit:
! retained delta 4.19E-1 -> 5.03E3).  Reject the commit when the
! max |V1-V2|/max(V2,1e-4*peak) over either store exceeds
! FG_COMMIT_JUMP; the pre-refresh slot values stay in place, so
! the next refresh re-seeds from the committed (sane) field.
! Genuine factor-of-few rewrites still pass.  Direct slot fills
! (cells the store left at the floor) are excluded -- those are
! seeds, not overwrites.
!
	  JUMP=MAX(FG_JUMP_MAX(FLD,I_M_DEPTH_NEXT),
	1         FG_JUMP_MAX(FLD_P,I_P_DEPTH_NEXT))
	  IF(JUMP .GT. FG_COMMIT_JUMP)THEN
!$OMP CRITICAL (FG_REJ_UPD)
	    FG_COMMIT_REJ=FG_COMMIT_REJ+1
!$OMP END CRITICAL (FG_REJ_UPD)
	    IF(ALLOCATED(FG_RAY_VC))FG_RAY_VC(LS)=1
	    IF(ALLOCATED(FG_RAY_DCH))FG_RAY_DCH(LS)=JUMP
	  ELSE IF(FG_COMMIT_CONV .AND.
	1    (.NOT. (IXIT_M .GE. 0 .AND. IXIT_M .LE. 2) .OR.
	2     .NOT. (IXIT_P .GE. 0 .AND. IXIT_P .LE. 2) .OR.
	3     MAX(RES_M,RES_P) .GT. KN_RAY_TOL))THEN
!
! edited by Devin, 17-Oct-2026: [FG_COMMIT_CONV] fail-closed
! commits.  A stall, pass-cap or plateau exit (IXIT 3/4/6, or an
! unset code) means the fixed-point solve ended above KN_RAY_TOL;
! committing its last iterate wrote that unconverged field into
! the lagged store bounded only by the FG_COMMIT_JUMP magnitude --
! the w19 runaway class (residual 0.465 vs tol 2E-3 committed with
! retained delta 21).  The residual test also catches the bulk
! exit: IXIT=2 fires on the bad-cell fraction alone, so a ray with
! residual ~0.5 and <0.1% bad cells was still classed converged.
! With FG_COMMIT_NC>0 the ray takes a capped provisional commit
! (per-cell change bounded by NC; the store keeps tracking a
! contracting field but cannot ratchet faster than NC/refresh).
! NC=0 skips the write: the store keeps the pre-refresh
! generation, a sane consistent seed.  The residual still reports
! RES_M/RES_P, so the maturity gate reads the ray as unconverged
! either way.
!
	    IF(FG_COMMIT_NC .GT. 0.0_LDP)THEN
	      CALL WRITE_STORE_CAP(FLD,I_M_DEPTH_NEXT,DCH_M)
	      CALL WRITE_STORE_CAP(FLD_P,I_P_DEPTH_NEXT,DCH_P)
!$OMP CRITICAL (FG_CAP_UPD)
	      FG_NC_CAP=FG_NC_CAP+1
!$OMP END CRITICAL (FG_CAP_UPD)
!$OMP CRITICAL (FG_DCHG_UPD)
	      FG_RF_DCHG=MAX(FG_RF_DCHG,DCH_M,DCH_P)
!$OMP END CRITICAL (FG_DCHG_UPD)
	      IF(ALLOCATED(FG_RAY_VC))FG_RAY_VC(LS)=3
	      IF(ALLOCATED(FG_RAY_DCH))FG_RAY_DCH(LS)=
	1        MAX(DCH_M,DCH_P)
	    ELSE
!$OMP CRITICAL (FG_NC_UPD)
	      FG_NC_REJ=FG_NC_REJ+1
!$OMP END CRITICAL (FG_NC_UPD)
	      IF(ALLOCATED(FG_RAY_VC))FG_RAY_VC(LS)=2
	      IF(ALLOCATED(FG_RAY_DCH))FG_RAY_DCH(LS)=JUMP
	    END IF
	  ELSE
	    CALL WRITE_STORE(FLD,I_M_DEPTH_NEXT,DCH_M)
	    CALL WRITE_STORE(FLD_P,I_P_DEPTH_NEXT,DCH_P)
!$OMP CRITICAL (FG_DCHG_UPD)
	    FG_RF_DCHG=MAX(FG_RF_DCHG,DCH_M,DCH_P)
!$OMP END CRITICAL (FG_DCHG_UPD)
	    IF(ALLOCATED(FG_RAY_VC))FG_RAY_VC(LS)=0
	    IF(ALLOCATED(FG_RAY_DCH))FG_RAY_DCH(LS)=MAX(DCH_M,DCH_P)
	  END IF
	END IF
	DEALLOCATE(FLD,FLD_P)
!
800	CONTINUE
!
! A diverged ray reports a 100% residual so the maturity gate reads
! it as unconverged; NBAD/NTOT keep the last completed pass counts.
!
! W2 (17-Oct-2026): heal the ray's retained stores before the next
! refresh re-seeds from them.  A diverged ray skips the store
! write, so a poisoned slot value (non-finite, or >FG_DIV_HEAL_CAP
! times the column median) would re-seed and re-diverge the next
! refresh -- the mechanism that once held GS_MAXCHG=1.0 forever.
! Zeroed slots take the seeded floor path, mirroring the MOM
! heal (mom_j_cmf_nonmono_v1.f).
!
	IF(NDIV .EQ. 1)THEN
	  CALL HEAL_STORE(I_M_DEPTH_NEXT)
	  CALL HEAL_STORE(I_P_DEPTH_NEXT)
	  RESID=1.0_LDP
	  NXIT=5
	  IF(ALLOCATED(FG_RAY_VC))FG_RAY_VC(LS)=4
	  IF(ALLOCATED(FG_RAY_DCH))FG_RAY_DCH(LS)=0.0_LDP
	ELSE
	  RESID=MAX(RES_M,RES_P)
	  IF(RES_M .GE. RES_P)THEN
	    NXIT=IXIT_M
	  ELSE
	    NXIT=IXIT_P
	  END IF
	END IF
!
! Export the worse direction's worst-cell record for the refresh
! summary; a diverged I_M keeps its own record since I_P never ran.
!
	IF(ALLOCATED(FG_RAY_WI))THEN
	  IF(RES_M .GE. RES_P)THEN
	    FG_RAY_WI(LS)=WRC_IM
	    FG_RAY_WL(LS)=WRC_LM
	    FG_RAY_WS(LS)=WRC_SM
	  ELSE
	    FG_RAY_WI(LS)=WRC_I
	    FG_RAY_WL(LS)=WRC_L
	    FG_RAY_WS(LS)=WRC_SAME
	  END IF
	END IF
	NPASS=MAX(NP_M,NP_P)
	NBAD=NB_M+NB_P
	NTOT=NT_M+NT_P
!
! W3 (17-Oct-2026): the best-pass census.  For stall/cap exits the
! last pass still counts every above-tolerance cell, so the gate
! reads MIN(GS_NBAD, GS_NBAD_BLK) -- the best pass is what the
! fixed-point map actually reached.
!
	NBBLK=NB_BM+NB_BP
	NTBLK=NT_BM+NT_BP
!
900	CONTINUE
	NP2=NP2_CNT
	IF(ALLOCATED(CHI_L))THEN
	  DEALLOCATE(CHI_L,ETA_L,DCHR_L,Q_L,SR_L,SP_L,S_L,DS_L,
	1    DTAU_L,A0L,A1L,A2L,A3L,A4L,EEL,E0L,E1L,E2L,E3L,
	2    IKC1,IKC2,JPR1,JPR2,WINX,T1V,IGRD,IM_NI,RES_FREQ)
	END IF
	IF(ALLOCATED(PK_DTAU))DEALLOCATE(PK_DTAU)
	IF(ALLOCATED(PK_A0))DEALLOCATE(PK_A0)
	IF(ALLOCATED(PK_A1))DEALLOCATE(PK_A1)
	IF(ALLOCATED(PK_A2))DEALLOCATE(PK_A2)
	IF(ALLOCATED(PK_A3))DEALLOCATE(PK_A3)
	IF(ALLOCATED(PK_A4))DEALLOCATE(PK_A4)
	IF(ALLOCATED(PK_SR))DEALLOCATE(PK_SR)
	IF(ALLOCATED(PK_Q))DEALLOCATE(PK_Q)
	IF(ALLOCATED(PK_ETANI))DEALLOCATE(PK_ETANI)
	IF(ALLOCATED(PK_CHINI))DEALLOCATE(PK_CHINI)
	IF(ALLOCATED(SPW))DEALLOCATE(SPW)
	IF(ALLOCATED(FLD))DEALLOCATE(FLD)
	RETURN
!
	CONTAINS
!
! ----------------------------------------------------------------
! Warm-start the per-ray field from the lagged depth store.  The
! in-sweep store layout is shifted: slot s holds frequency s+1
! (FREQ_IDX_ITER is 0 during the first frequency's call), while FLD
! is indexed by true frequency l.  FLD(:,1) has no stored counter-
! part and stays zero.
! ----------------------------------------------------------------
	SUBROUTINE SEED_FLD(FLD,STORE)
	REAL(KIND=LDP) FLD(NI,NCF)
	REAL(KIND=LDP) STORE(:,:,:)
	INTEGER KS,LS2,JP2
	DO KS=1,KN_ND
	  IF(COMP_DEPTH_IDX(KS) .GE. 1)THEN
	    JP2=J_PNT(KS,LS)
	    IF(JP2 .GE. 1 .AND. JP2 .LE. NI)THEN
	      DO LS2=2,NCF
	        FLD(JP2,LS2)=STORE(COMP_DEPTH_IDX(KS),LS,LS2-1)
	      END DO
	    END IF
	  END IF
	END DO
	RETURN
	END
!
! ----------------------------------------------------------------
! Alternating-direction fixed-point iteration on the coupled field.
! Odd passes ascend in frequency (fresh upwind for expanding cells),
! even passes descend (fresh upwind for compressing cells).  A ray
! with no usable compressing slots is exact after the first pass.
! ----------------------------------------------------------------
	SUBROUTINE SWEEP_FLD(FLD,INWARD,RES,NPO,NB_O,NT_O,IXIT,
	1   NB_B,NT_B)
	REAL(KIND=LDP) FLD(NI,NCF),RES
	LOGICAL INWARD
	INTEGER NPO,LB,LE,LDI,IXIT
	INTEGER*8 NB_O,NT_O,NB_B,NT_B
	INTEGER NSTALL,NPLAT(2)
	REAL(KIND=LDP) BEST_RES
	INTEGER*8 BEST_NB,BEST_NT
	REAL(KIND=LDP) RATIO_PV(2)
	LOGICAL ASC,HOTMODE,L_PLAT
	NPO=0
	NSTALL=0
	NPLAT(1)=0; NPLAT(2)=0
	L_PLAT=.FALSE.
	IXIT=-1
	CLEAN_EXIT=.FALSE.
	DIVERGED=.FALSE.
	HOTMODE=.FALSE.
	RES=1.0E+30_LDP
	RES_PREV(1)=0.0_LDP; RES_PREV(2)=0.0_LDP
	RATIO_PV(1)=0.0_LDP; RATIO_PV(2)=0.0_LDP
	RATIO_OM=0.0_LDP
	WRC_I=0; WRC_L=0; WRC_PI=0; WRC_PL=0; WRC_SAME=0
!
! SGN_PRV/NOSC_FL are host variables: without this reset the first
! pass of the I_P sweep compares its sign against the I_M sweep's
! last accumulator, scoring a spurious flip and damping OMEGA_RAY
! on stale state (SWEEP_FLD_ORD already resets KN_SGN_PREV/NOSC at
! its own head).
!
	SGN_ACC=0.0_LDP; SGN_PRV=0.0_LDP; NOSC_FL=0
	OMEGA_RAY=KN_RAY_OMEGA
	BEST_RES=HUGE(RES)
	BEST_NB=0_8; BEST_NT=0_8
	DO
	  NPO=NPO+1
!
! W7 (17-Oct-2026): hot-pass extension.  When a stall or cap exit
! fires with residual left only inside a contiguous hot-frequency
! run, the sweep continues restricted to that window for up to
! KN_RAY_HOTPASSES extra passes; KN_HOT_SPENT bounds the
! refresh-wide spend.
!
	  IF(HOTMODE)THEN
	    NHOT=NHOT+1
!$OMP CRITICAL (KN_HOT_SPD)
	    KN_HOT_SPENT=KN_HOT_SPENT+1
	    L_HOT_OVER=KN_HOT_SPENT .GT. KN_HOT_TOTAL_PASSES
!$OMP END CRITICAL (KN_HOT_SPD)
	    IF(NHOT .GT. KN_RAY_HOTPASSES .OR. L_HOT_OVER)THEN
	      IXIT=4
	      EXIT
	    END IF
	  END IF
	  ASC=(MOD(NPO,2) .EQ. 1)
	  RES=0.0_LDP
	  RES_FREQ(1:NCF)=0.0_LDP
	  NBAD_CUR=0; NTOT_CUR=0
	  SGN_ACC=0.0_LDP
	  IF(ASC)THEN
	    LB=1; LE=NCF; LDI=1
	  ELSE
	    LB=NCF; LE=1; LDI=-1
	  END IF
	  DO L=LB,LE,LDI
	    IF(HOTMODE .AND. (L .LT. W_LO .OR. L .GT. W_HI))CYCLE
	    IF(HAVE_PACK)THEN
!
! Replay the pass-invariant coefficients built once per ray (same
! values BUILD_FREQ would compute; see the pack build above).
! ETA_L/CHI_L are needed at NI only, for the DIFFUSION inner
! boundary in MARCH.
!
	      DTAU_L(1:NI)=PK_DTAU(1:NI,L)
	      A0L(1:NI)=PK_A0(1:NI,L)
	      A1L(1:NI)=PK_A1(1:NI,L)
	      A2L(1:NI)=PK_A2(1:NI,L)
	      A3L(1:NI)=PK_A3(1:NI,L)
	      A4L(1:NI)=PK_A4(1:NI,L)
	      SR_L(1:NI)=PK_SR(1:NI,L)
	      Q_L(1:NI)=PK_Q(1:NI,L)
	      ETA_L(NI)=PK_ETANI(L)
	      CHI_L(NI)=PK_CHINI(L)
	    ELSE
	      CALL BUILD_FREQ(L)
	    END IF
	    CALL MATCH_FREQ(L)
	    CALL SRC_PRIME(FLD,L)
	    CALL SLOPES()
	    CALL MARCH(FLD,L,INWARD,RES)
	  END DO
	  NB_O=NBAD_CUR; NT_O=NTOT_CUR
	  IF(ALLOCATED(KN_RES_PASS) .AND. .NOT. HOTMODE)THEN
	    IF(NPO .LE. SIZE(KN_RES_PASS))THEN
!$OMP CRITICAL (KN_RES_UPD)
	      KN_RES_PASS(NPO)=MAX(KN_RES_PASS(NPO),RES)
!$OMP END CRITICAL (KN_RES_UPD)
	    END IF
	  END IF
!
! Divergence containment: a non-finite residual, or a same-direction
! blow-up (RES more than 1e3 above the best pass so far, checked once
! two passes have completed), means the iteration is diverging.  Exit
! at once: the caller skips the write-back, so the runaway field never
! reaches the store.  No FLD copy is kept (NI*NCF*8 B per thread).
!
	  IF(RES .NE. RES .OR. RES .GT. HUGE(RES) .OR.
	1     (NPO .GE. 3 .AND. RES .GT. 1.0E+3_LDP*BEST_RES))THEN
	    DIVERGED=.TRUE.
	    IXIT=5
	    EXIT
	  END IF
!
! W3 (17-Oct-2026): best-pass census.  The stall/cap exits below
! carry NB_O/NT_O from the last pass; the gate reads
! MIN(GS_NBAD, GS_NBAD_BLK), so track the census of the pass the
! map actually reached.
!
	  IF(RES .LT. BEST_RES)THEN
	    BEST_RES=RES
	    BEST_NB=NBAD_CUR
	    BEST_NT=NTOT_CUR
	  END IF
!
! Adaptive relaxation keyed on the SAME-DIRECTION residual history
! (pass N vs N-2): alternating-direction sweeps interleave a
! compressing-upwind pass with an expanding-upwind pass, so the raw
! pass-to-pass sequence can sawtooth even while each direction
! converges geometrically.  OMEGA_RAY must never exceed 1 here:
! MARCH applies it inside the spatial short-characteristics
! recurrence, where each cell relaxes from the just-relaxed upwind
! value, so an omega>1 error amplification compounds per cell along
! the ray.  Knop Paper II uses plain GS/Jacobi -- no over-relaxation.
! Same-direction growth signals an oscillatory (-1-eigenvalue) cycle
! -> damp; fast contraction lets omega recover toward 1 only.
!
	  IDIR=1
	  IF(.NOT. ASC)IDIR=2
	  IF(RES_PREV(IDIR) .GT. 0.0_LDP)THEN
	    RATIO_OM=RES/RES_PREV(IDIR)
!
! edited by Devin, 14-Oct-2026: suppress damping while the
! same-direction baseline (pass NPO-2) is still inside the
! warm-start seed transient (passes 1-2).  A transient baseline
! reads a spurious growth ratio and can drive OMEGA_RAY to the
! 0.45 floor on seeding artifacts alone, over-damping the ray for
! the whole refresh.  Recovery toward 1 is left unfloored: it only
! undoes earlier damping.
!
	    IF(RATIO_OM .GT. 1.05_LDP .AND. NPO .GE. 5)THEN
	      OMEGA_RAY=MAX(0.45_LDP,OMEGA_RAY*0.7_LDP)
	    ELSE IF(RATIO_OM .LT. 0.5_LDP)THEN
	      OMEGA_RAY=MIN(1.0_LDP,OMEGA_RAY*1.15_LDP)
	    END IF
!
! W6 (17-Oct-2026): sign-flip damping ported from SWEEP_FLD_ORD --
! on a period-2 limit cycle (lambda ~ -1) the signed per-cell
! update sum flips every pass; two consecutive flips damp
! OMEGA_RAY.  Marginal contraction (lambda ~ +0.99) keeps its
! sign and is left to the stall exits.
!
	    IF(SGN_PRV .NE. 0.0_LDP .AND.
	1     SGN_ACC*SGN_PRV .LT. 0.0_LDP)THEN
	      NOSC_FL=NOSC_FL+1
	    ELSE
	      NOSC_FL=0
	    END IF
	    SGN_PRV=SGN_ACC
	    IF(NOSC_FL .GE. 2 .AND. NPO .GE. 5 .AND.
	1     OMEGA_RAY .GT. 0.45_LDP)THEN
	      OMEGA_RAY=MAX(0.45_LDP,OMEGA_RAY*0.7_LDP)
	    END IF
!
! Stall detection: count consecutive same-direction comparisons with
! residual ratio above KN_RAY_STALL (VADAT, default 0.98 = <2%
! improvement).  The lagged field is rebuilt from fresh opacities
! every great iteration, so polishing a slowly-contracting fixed
! point buys little -- NSTALL>=2 (four straight passes of
! stagnation) exits below once past the initial transient (NPO>=8).
! The 0.98 default keeps rays contracting at >=2%/pass alive; at
! 0.90 those rays held nearly all unconverged cells on this model.
!
	    IF(NPO .GE. 8)THEN
	      IF(RATIO_OM .GT. KN_RAY_STALL)THEN
	        NSTALL=NSTALL+1
	      ELSE
	        NSTALL=0
	      END IF
!
! W6 (17-Oct-2026, [FG_PLATEAU_DET], default T): two additions to
! the same-direction-only stall.  (a) Plateau: three consecutive
! same-direction ratios in [0.95,1.05] -- flat, not contracting.
! (b) Combined-direction thrash: MAX ratio >0.98 AND the previous
! same-direction ratio >0.95 -- both directions stagnate, so the
! ray cannot polish out; exiting at pass ~10 lets W3 count the
! best pass instead of riding KN_RAY_MAXPASS.
!
	      IF(FG_PLATEAU_DET)THEN
	        IF(ABS(RATIO_OM-1.0_LDP) .LT. 0.05_LDP)THEN
	          NPLAT(IDIR)=NPLAT(IDIR)+1
	        ELSE
	          NPLAT(IDIR)=0
	        END IF
	        IF(NPLAT(IDIR) .GE. 3)THEN
	          NSTALL=MAX(NSTALL,2)
	          L_PLAT=.TRUE.
	        END IF
	        IF(RATIO_OM .GT. 0.98_LDP .AND.
	1         RATIO_PV(IDIR) .GT. 0.95_LDP)NSTALL=MAX(NSTALL,2)
	      END IF
	      RATIO_PV(IDIR)=RATIO_OM
	    END IF
	  END IF
	  RES_PREV(IDIR)=RES
!
! Worst-cell persistence: if the same (depth,freq) cell held the
! pass max last pass too, one unstable mode dominates; a moving max
! means the residual is being reset by different cells each pass.
!
	  IF(WRC_I .EQ. WRC_PI .AND. WRC_L .EQ. WRC_PL .AND.
	1    WRC_I .GT. 0)THEN
	    WRC_SAME=WRC_SAME+1
	  ELSE
	    WRC_SAME=1
	  END IF
	  WRC_PI=WRC_I; WRC_PL=WRC_L
	  IF(.NOT. HAS_CP)THEN
!
!	    No compressing cells: the single pass is the exact solve, so
!	    the pass-1 change is a seed transient, not a fixed-point
!	    residual.  Report 0 so this ray cannot pin the gate; the
!	    per-frequency residuals are also dropped so this ray cannot
!	    inflate KN_HOT_CNT.
!
	    CLEAN_EXIT=.TRUE.
	    IXIT=0
	    RES=0.0_LDP
	    RES_FREQ(1:NCF)=0.0_LDP
	    EXIT
	  END IF
	  IF(RES .LT. KN_RAY_TOL)THEN
	    CLEAN_EXIT=.TRUE.
	    IXIT=1
	    EXIT
	  END IF
!
! Bulk exit: once the significant-cell bad fraction is at or below
! the maturity-gate level, further passes only polish cells the gate
! ignores.  Require a minimum of 6 passes so the first couple of
! sweep directions cannot exit on a partially-propagated field.
!
! [FG_BULK_CONTR]: the bulk test ignores the residual magnitude, so
! a flat noncontracting ray (RATIO_OM~1) with a sub-0.1% bad
! fraction exited "converged" at pass 6 -- before the stall and
! plateau classifiers (pass>=8) could run.  w19 committed 8-27x
! retained deltas through this door at residual 0.25-0.49.  With the
! flag set the exit also requires same-direction contraction;
! noncontracting rays fall through to the stall/plateau/cap classes.
!
	  IF(NPO .GE. 6 .AND. NT_O .GT. 0 .AND.
	1     NB_O .LE. MAX(50_8,NT_O/1000) .AND.
	2     (.NOT. FG_BULK_CONTR .OR. RATIO_OM .LE. KN_RAY_STALL))THEN
	    IXIT=2
	    EXIT
	  END IF
!
! W7 (17-Oct-2026): hot-pass extension.  A stall or pass-cap exit
! with residual left only inside a contiguous hot-frequency run
! continues restricted to that window instead of exiting -- the
! marginal |Q|/(chi+|Q|) modes live on contiguous compression
! frequencies, so extra passes there cost little per pass.
!
	  IF(NSTALL .GE. 2)THEN
	    IF(.NOT. HOTMODE .AND. TRY_HOT())THEN
	      HOTMODE=.TRUE.
	    ELSE
	      IXIT=3
	      IF(L_PLAT)IXIT=6
	      EXIT
	    END IF
	  END IF
	  IF(NPO .GE. KN_RAY_MAXPASS)THEN
	    IF(.NOT. HOTMODE .AND. TRY_HOT())THEN
	      HOTMODE=.TRUE.
	    ELSE
	      IXIT=4
	      EXIT
	    END IF
	  END IF
	END DO
!
! A ray that exited cleanly (RES<TOL, or the exact single pass of a
! ray with no compressing cells) has no unconverged cells by
! definition; counting its last-pass seed transient would inflate
! the bulk bad-cell fraction the maturity gate reads.
!
	IF(CLEAN_EXIT)NB_O=0
!
! W3 (17-Oct-2026): export the best-pass census.  Stall/cap exits
! carry the last-pass NB_O/NT_O; the gate takes
! MIN(GS_NBAD, GS_NBAD_BLK) so a ray whose best pass was inside the
! bulk band cannot keep the maturity gate pinned.
!
	NB_B=BEST_NB
	NT_B=BEST_NT
	IF(CLEAN_EXIT)NB_B=0_8
	RETURN
	END
!
! ----------------------------------------------------------------
! Modified source S'(r,l) for a single cell of the ordered march.
! Identical stencil to SRC_PRIME: expanding cells take the same-row
! upwind FLD(r,l-1), compressing cells take the JPR-interpolated
! upwind at l+NFM (matched slot under DO_NONMONO_FREQ_MATCH).  The
! caller guarantees r is inside 1..NI and uses the ring slot of r.
! ----------------------------------------------------------------
	SUBROUTINE KN_ORD_SPV(FLD,R,L,SPV)
	REAL(KIND=LDP) FLD(NI,NCF),SPV
	INTEGER R,L,NNF2
	IF(GAM(R,LS) .GE. 0.0_LDP)THEN
	  SPV=PK_SR(R,L)
	  IF(L .GT. 1)SPV=PK_SR(R,L)+PK_Q(R,L)*FLD(R,L-1)
	ELSE
	  NNF2=1
	  IF(FG_FREQ_MATCH .AND. L .LT. NCF .AND.
	1     KN_DLN_N(L) .GT. 0.0_LDP)THEN
	    T1=DRF_L(R)/KN_DLN_N(L)
	    IF(T1 .GE. DBLE(NCF-L))THEN
	      NNF2=NCF-L
	    ELSE
	      NNF2=MAX(1,CEILING(T1))
	    END IF
	  END IF
	  IF(L+NNF2 .LE. NCF)THEN
	    IF(IKC1(R) .EQ. 1 .AND. IKC2(R) .EQ. 1)THEN
	      INEXT_V=(1.0_LDP-WINX(R))*FLD(JPR1(R),L+NNF2)+
	1             WINX(R)*FLD(JPR2(R),L+NNF2)
	    ELSE IF(IKC1(R) .EQ. 1)THEN
	      INEXT_V=FLD(JPR1(R),L+NNF2)
	    ELSE IF(IKC2(R) .EQ. 1)THEN
	      INEXT_V=FLD(JPR2(R),L+NNF2)
	    ELSE
	      INEXT_V=0.0_LDP
	    END IF
	    SPV=PK_SR(R,L)+ABS(PK_Q(R,L))*INEXT_V
	  ELSE
	    SPV=PK_SR(R,L)
	  END IF
	END IF
	RETURN
	END
!
! ----------------------------------------------------------------
! Fill ring slot of row R with S'(R,:) evaluated on the current FLD
! state.  Used for the boundary rows and for the "preview" of the
! not-yet-marched neighbour row -- the lagged C-corner source of
! Knop et al. (2009, Paper II Eqs. 31/32).
! ----------------------------------------------------------------
	SUBROUTINE KN_ORD_SPW_FILL(FLD,R)
	REAL(KIND=LDP) FLD(NI,NCF)
	INTEGER R,L2
	DO L2=1,NCF
	  CALL KN_ORD_SPV(FLD,R,L2,T1)
	  SPW(L2,MOD(R-1,4)+1)=T1
	END DO
	RETURN
	END
!
! ----------------------------------------------------------------
! Secant slope S(r) = (S'(r+1)-S'(r))/DTAU(r) on the SPW ring,
! matching SLOPES including the thin-cell guard.
! ----------------------------------------------------------------
	REAL(KIND=LDP) FUNCTION KN_ORD_SEC(R,L)
	INTEGER R,L
	IF(PK_DTAU(R,L) .GT. 1.0E-30_LDP)THEN
	  KN_ORD_SEC=(SPW(L,MOD(R,4)+1)-SPW(L,MOD(R-1,4)+1))/
	1    PK_DTAU(R,L)
	ELSE
	  KN_ORD_SEC=0.0_LDP
	END IF
	RETURN
	END
!
! ----------------------------------------------------------------
! Monotonic-parabola slope DS(r) on the SPW ring -- the SLOPES
! construction including the endpoint extrapolations and the
! sign-and-min clamps.  Reads S' rows r-2..r+1 (via the secants),
! so the caller must hold those rows in the ring.
! ----------------------------------------------------------------
	REAL(KIND=LDP) FUNCTION KN_ORD_DS(R,L)
	INTEGER R,L
	REAL(KIND=LDP) SA,SB,DD,TDT
	IF(R .EQ. 1)THEN
	  SA=KN_ORD_SEC(1,L)
	  SB=KN_ORD_SEC(2,L)
	  TDT=PK_DTAU(1,L)+PK_DTAU(2,L)
	  IF(ABS(TDT) .LT. 1.0E-30_LDP)THEN
	    DD=0.0_LDP
	  ELSE
	    DD=SA+(SA-SB)*PK_DTAU(1,L)/TDT
	  END IF
	  KN_ORD_DS=(SIGN(1.0_LDP,SA)+SIGN(1.0_LDP,DD))*
	1    MIN(ABS(SA),0.5_LDP*ABS(DD))
	ELSE IF(R .EQ. NI)THEN
	  SA=KN_ORD_SEC(NI-1,L)
	  SB=KN_ORD_SEC(NI-2,L)
	  TDT=PK_DTAU(NI-2,L)+PK_DTAU(NI-1,L)
	  IF(ABS(TDT) .LT. 1.0E-30_LDP)THEN
	    DD=0.0_LDP
	  ELSE
	    DD=SA+(SA-SB)*PK_DTAU(NI-1,L)/TDT
	  END IF
	  KN_ORD_DS=(SIGN(1.0_LDP,SA)+SIGN(1.0_LDP,DD))*
	1    MIN(ABS(SA),0.5_LDP*ABS(DD))
	ELSE
	  SA=KN_ORD_SEC(R-1,L)
	  SB=KN_ORD_SEC(R,L)
	  TDT=PK_DTAU(R-1,L)+PK_DTAU(R,L)
	  IF(ABS(TDT) .LT. 1.0E-30_LDP)THEN
	    DD=0.0_LDP
	  ELSE
	    DD=(SA*PK_DTAU(R,L)+SB*PK_DTAU(R-1,L))/TDT
	  END IF
	  KN_ORD_DS=(SIGN(1.0_LDP,SA)+SIGN(1.0_LDP,SB))*
	1    MIN(ABS(SA),ABS(SB),0.5_LDP*ABS(DD))
	END IF
	RETURN
	END
!
! ----------------------------------------------------------------
! Physically ordered fixed-point iteration (VADAT KN_RAY_ORDERED,
! default T): the (i,l) march of Knop et al. (2009, Paper II Sec.2).
! Rows are marched in the characteristic direction -- j ascending
! for I_M, descending for I_P -- and the frequency order inside a
! row follows sign(GAM(j)): ascending on expanding rows (fresh
! same-row upwind at l-1), descending on compressing rows (fresh
! self-mapped upwind at l+NFM when JPR falls on the row itself).
!
! The only lagged quantities are the source of the one row beyond
! the marching front -- SP(j+1) for I_M, SP(j-1) for I_P, entering
! through the parabolic slope term DS -- plus compressing
! interpolant donors JPR on the unmarched side of the front.  Both
! are the corner terms Knop also leaves lagged (their C+/- coupling
! of Eqs. 31/32); the main |Q| wavelength link never waits a pass.
! The DS corner vanishes where the monotonicity clamp zeroes the
! slopes, i.e. exactly in the steep compression gradients that
! stall the alternating scheme.
!
! Same fixed point, same residual metric, same guards and exits as
! SWEEP_FLD; the pass-order alternation is dropped because the
! ordering is fixed by the characteristic, not by pass parity.
! ----------------------------------------------------------------
	SUBROUTINE SWEEP_FLD_ORD(FLD,INWARD,RES,NPO,NB_O,NT_O,IXIT,
	1   NB_B,NT_B)
	REAL(KIND=LDP) FLD(NI,NCF),RES
	LOGICAL INWARD
	INTEGER NPO,IXIT
	INTEGER*8 NB_O,NT_O,NB_B,NT_B
	INTEGER NSTALL,J,L,LB,LE,LDI
	INTEGER NOSC,NPLAT
	LOGICAL L_PLAT
	REAL(KIND=LDP) KN_SGN_SUM,KN_SGN_PREV
	REAL(KIND=LDP) BEST_RES
	INTEGER*8 BEST_NB,BEST_NT
	REAL(KIND=LDP) RATIO_PV
	LOGICAL HOTMODE
	REAL(KIND=LDP), ALLOCATABLE :: VMAXF(:)
	REAL(KIND=LDP) SPJ,DSD,DSU
!
! FG_ORD_AITKEN (VADAT [KN_ORD_AITKEN], default F): guarded Aitken
! delta-squared extrapolation of the pass iterate, aimed at the
! marginally contractive |Q|->1 compression cells that stall the
! march against the bulk-maturity exit.  FLD_PRV holds the last
! committed iterate, DLT_PRV the previous pass delta; the pair of
! consecutive deltas gives the per-cell contraction ratio.  Restart
! convention: the extrapolated value replaces the pass result, so
! the next pass marches from the accelerated state.
!
	REAL(KIND=LDP), ALLOCATABLE :: FLD_PRV(:,:),DLT_PRV(:,:)
	REAL(KIND=LDP) DEL_CUR,DEL_OLD,XACC
	LOGICAL DO_AIT
!
! FG_ORD_AA (VADAT [KN_ORD_AA], default F; secant depth
! [KN_ORD_AA_M], default 2): type-II Anderson acceleration of the
! pass iterate, mirroring the runtime-tested MOM refresh AA
! (mom_j_cmf_nonmono_v1.f) plus a residual-acceptance rollback.
! History pairs are (x_i, r_i=g(x_i)-x_i) where g is one relaxed
! march; slot 1 is the newest pair.  Guards: residual-weighted
! Gram system, exclusion of cells outside the affine regime
! (weighted residual/secant > AA_GCAP), Tikhonov ridge,
! ||corr||>1E3||f|| veto (garbage bound only -- a marginal mode's
! correct step is ~||f||/(1-lambda)), 10|r| per-cell step cap,
! positivity/cap fallback.  If the pass after an applied mix does
! not improve,
! the pre-mix iterate is restored, the history cleared, and AA
! disabled for the rest of the solve.  Collection only starts at
! NPO>=4 (past the seed transient), so fast rays never pay the
! (M+1)-pair memory.
!
	REAL(KIND=LDP), ALLOCATABLE :: AA_XH(:,:,:),AA_RH(:,:,:)
	REAL(KIND=LDP) GRM_AA(8,8),WGR_AA(8,8),BRH_AA(8),WBR_AA(8)
	REAL(KIND=LDP) GAM_AA(8)
	REAL(KIND=LDP) AA_REF,FNRM,DNRM,DTR,TCOL,STEP,T3,T4
	REAL(KIND=LDP) WV,SC,XV
	INTEGER IPIV_AA(8),NH_AA,MP_AA,INFO_AA,JJ,KK
	INTEGER*8 NAA_SKIP
	REAL(KIND=LDP), PARAMETER :: AA_GCAP=1.0_LDP
	LOGICAL DO_AA,AA_PEND,AA_OFF
	EXTERNAL DGESV
!
! FG_ORD_P2DAMP (VADAT [KN_ORD_P2DAMP], default F): per-cell
! period-2 collapse for the bounded lambda~-1 modes the whole-ray
! KN_SGN_SUM detector cannot see (the m15_2p82d kndump run pinned
! the gate at a stable 0.53 max residual on ~20 near-floor cells
! flip-flopping ~2x/~0.5x per pass).  P2_SGN holds the sign of the
! last update delta at each significant cell, P2_FLIP the
! consecutive flip count, P2_DMP the damped state (0 = off,
! 1..FG_ORD_P2REL = damped; counts flip-free passes until release).
! int*1 keeps the per-ray state at 3 B/cell; allocated per solve
! so each OMP thread owns its own copy.
!
	INTEGER*1, ALLOCATABLE :: P2_SGN(:,:),P2_FLIP(:,:),P2_DMP(:,:)
	INTEGER*1 SGN_P
	LOGICAL DO_P2
!
	NPO=0
	NSTALL=0
	NOSC=0
	NPLAT=0
	L_PLAT=.FALSE.
	RATIO_PV=0.0_LDP
	HOTMODE=.FALSE.
	IXIT=-1
	CLEAN_EXIT=.FALSE.
	DIVERGED=.FALSE.
	RES=1.0E+30_LDP
	RES_PREV(1)=0.0_LDP
	RES_PREV(2)=0.0_LDP
	RATIO_OM=0.0_LDP
	WRC_I=0; WRC_L=0; WRC_PI=0; WRC_PL=0; WRC_SAME=0
	KN_SGN_SUM=0.0_LDP
	KN_SGN_PREV=0.0_LDP
	OMEGA_RAY=KN_RAY_OMEGA
	BEST_RES=HUGE(RES)
	BEST_NB=0_8; BEST_NT=0_8
	ALLOCATE(VMAXF(NCF),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,'(A)')'  FG-KNOP sweep alloc failed: VMAXF'
	  DIVERGED=.TRUE.
	  IXIT=5
	  RETURN
	END IF
! AA supersedes the per-cell secant step; running both on the
! same iterate lets a bad Aitken jump corrupt the AA history.
!
	DO_AIT=FG_ORD_AITKEN .AND. .NOT. FG_ORD_AA
	IF(DO_AIT)THEN
	  ALLOCATE(FLD_PRV(NI,NCF),DLT_PRV(NI,NCF),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    DO_AIT=.FALSE.
	  ELSE
	    FLD_PRV(1:NI,1:NCF)=FLD(1:NI,1:NCF)
	    DLT_PRV(1:NI,1:NCF)=0.0_LDP
	  END IF
	END IF
	DO_AA=FG_ORD_AA
	AA_PEND=.FALSE.
	AA_OFF=.FALSE.
	NH_AA=0
	AA_REF=0.0_LDP
!
	DO_P2=FG_ORD_P2DAMP
	IF(DO_P2)THEN
	  ALLOCATE(P2_SGN(NI,NCF),P2_FLIP(NI,NCF),P2_DMP(NI,NCF),
	1    STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    DO_P2=.FALSE.
	  ELSE
	    P2_SGN=0_1; P2_FLIP=0_1; P2_DMP=0_1
	  END IF
	END IF
!
	DO
	  NPO=NPO+1
!
! W7 (17-Oct-2026): hot-pass extension spend accounting, same rule
! as SWEEP_FLD -- up to KN_RAY_HOTPASSES extra passes per ray,
! KN_HOT_TOTAL_PASSES per refresh.
!
	  IF(HOTMODE)THEN
	    NHOT=NHOT+1
!$OMP CRITICAL (KN_HOT_SPD)
	    KN_HOT_SPENT=KN_HOT_SPENT+1
	    L_HOT_OVER=KN_HOT_SPENT .GT. KN_HOT_TOTAL_PASSES
!$OMP END CRITICAL (KN_HOT_SPD)
	    IF(NHOT .GT. KN_RAY_HOTPASSES .OR. L_HOT_OVER)THEN
	      IXIT=4
	      EXIT
	    END IF
	  END IF
!
! AA history: save the pre-march iterate as the newest pair input.
! Lazy allocation at NPO>=4 keeps fast rays free of the
! (M+1)*2*NI*NCF history memory.
!
	  IF(DO_AA .AND. .NOT. AA_OFF .AND. NPO .GE. 4)THEN
	    IF(.NOT. ALLOCATED(AA_XH))THEN
	      ALLOCATE(AA_XH(NI,NCF,FG_ORD_AA_M+1),
	1        AA_RH(NI,NCF,FG_ORD_AA_M+1),STAT=IOS)
	      IF(IOS .NE. 0)DO_AA=.FALSE.
	    END IF
	    IF(ALLOCATED(AA_XH))THEN
	      DO JJ=MIN(NH_AA+1,FG_ORD_AA_M+1),2,-1
	        AA_XH(1:NI,1:NCF,JJ)=AA_XH(1:NI,1:NCF,JJ-1)
	        AA_RH(1:NI,1:NCF,JJ)=AA_RH(1:NI,1:NCF,JJ-1)
	      END DO
	      AA_XH(1:NI,1:NCF,1)=FLD(1:NI,1:NCF)
	      NH_AA=MIN(NH_AA+1,FG_ORD_AA_M+1)
	    END IF
	  END IF
	  RES=0.0_LDP
	  RES_FREQ(1:NCF)=0.0_LDP
	  NBAD_CUR=0
	  NTOT_CUR=0
	  KN_SGN_SUM=0.0_LDP
!
! Column peaks at pass start for the significance floor, same role
! as VMAX_L in MARCH.
!
	  DO L=1,NCF
	    VMAXF(L)=MAX(MAXVAL(FLD(1:NI,L)),1.0E-25_LDP)
	  END DO
!
	  IF(INWARD)THEN
!
! Boundary row I=1: zero incident intensity, same residual
! accounting as MARCH.
!
	    DO L=1,NCF
	      VOLD=FLD(1,L)
	      IF(VOLD .GT. 0.0_LDP)THEN
	        T1=VOLD/MAX(VOLD,1.0E-4_LDP*VMAXF(L))
	        IF(T1 .GT. RES)THEN
	          WRC_I=1
	          WRC_L=L
	        END IF
	        RES=MAX(RES,T1)
	        RES_FREQ(L)=MAX(RES_FREQ(L),T1)
	      END IF
	      FLD(1,L)=0.0_LDP
	    END DO
	    CALL KN_ORD_SPW_FILL(FLD,1)
	    DO J=2,NI
!
! Preview of row j+1 built once from its pre-march state; it is the
! only not-fresh source this step reads (inside DS(j)).
!
	      IF(J+1 .LE. NI)CALL KN_ORD_SPW_FILL(FLD,J+1)
	      IF(GAM(J,LS) .GE. 0.0_LDP)THEN
	        LB=1; LE=NCF; LDI=1
	      ELSE
	        LB=NCF; LE=1; LDI=-1
	      END IF
	      DO L=LB,LE,LDI
	        IF(HOTMODE .AND. (L .LT. W_LO .OR. L .GT. W_HI))CYCLE
	        CALL KN_ORD_SPV(FLD,J,L,SPJ)
	        SPW(L,MOD(J-1,4)+1)=SPJ
	        DSD=KN_ORD_DS(J-1,L)
	        DSU=KN_ORD_DS(J,L)
	        VNEW=FLD(J-1,L)*PK_A0(J-1,L)+
	1           (SPW(L,MOD(J-2,4)+1)*PK_A1(J-1,L)+
	2           SPJ*PK_A2(J-1,L)+
	3           DSD*PK_A3(J-1,L)+DSU*PK_A4(J-1,L))
	        IF(VNEW .LT. 0.0_LDP .OR. VNEW .NE. VNEW)THEN
	          VNEW=1.0E-30_LDP
	          NNEG=NNEG+1
	        END IF
	        IF(VNEW .GT. NONMONO_INEXT_CAP)THEN
	          VNEW=NONMONO_INEXT_CAP
	          NNEG=NNEG+1
	        END IF
	        VOLD=FLD(J,L)
	        T1=ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAXF(L))
	        IF(T1 .GT. RES)THEN
	          WRC_I=J
	          WRC_L=L
	        END IF
	        RES=MAX(RES,T1)
	        RES_FREQ(L)=MAX(RES_FREQ(L),T1)
	        IF(MAX(VOLD,ABS(VNEW)) .GT. 1.0E-4_LDP*VMAXF(L))THEN
	          NTOT_CUR=NTOT_CUR+1
	          KN_SGN_SUM=KN_SGN_SUM+
	1           (VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAXF(L))
	          IF(T1 .GT. KN_RAY_TOL)NBAD_CUR=NBAD_CUR+1
!
! KN_DUMP: per-cell record of late-pass residual outliers.  Cells
! still changing >20% past the transient passes form the
! non-contractive set that pins the maturity gate; the record
! carries the coupling context so failing modes can be classified
! offline.  Capped per ray so a pathological solve cannot flood
! the log.
!
	          IF(PRINT_KNOP_DEBUG .AND. T1 .GT. 0.2_LDP .AND.
	1           NPO .GE. 8 .AND. KN_DUMP_N .LT. 2000)THEN
	            KN_DUMP_N=KN_DUMP_N+1
	            WRITE(LUER,'(A,I4,A,I4,A,I5,A,I5,4(A,ES11.4),3(A,I4))')
	1           ' KNDUMP p=',NPO,' r=',LS,' j=',J,' l=',L,
	2           ' res=',T1,' vo=',VOLD,' vn=',VNEW,
	3           ' gam=',GAM(J,LS),' jp1=',JPR1(J),' jp2=',JPR2(J),
	4           ' nfm=',NFM_L(J)
	          END IF
!
! P2 tracking (KN_ORD_P2DAMP): the sign of this cell's update delta
! is kept across passes; a significant cell that flips
! FG_ORD_P2N times in a row rides a lambda<0 mode the bulk
! KN_SGN_SUM cannot see.  P2_DMP counts flip-free passes while
! damped, so the collapse releases once the cell settles.
!
	          IF(DO_P2)THEN
	            SGN_P=0_1
	            IF(VNEW .GT. VOLD)SGN_P=1_1
	            IF(VNEW .LT. VOLD)SGN_P=-1_1
	            IF(P2_SGN(J,L) .NE. 0_1 .AND. SGN_P .NE. 0_1 .AND.
	1             SGN_P .NE. P2_SGN(J,L))THEN
	              IF(P2_FLIP(J,L) .LT. 100_1)
	1               P2_FLIP(J,L)=P2_FLIP(J,L)+1_1
	              IF(INT(P2_FLIP(J,L)) .GE. FG_ORD_P2N)
	1               P2_DMP(J,L)=1_1
	            ELSE
	              P2_FLIP(J,L)=0_1
	              IF(P2_DMP(J,L) .GT. 0_1)THEN
	                P2_DMP(J,L)=P2_DMP(J,L)+1_1
	                IF(INT(P2_DMP(J,L)) .GT. FG_ORD_P2REL)
	1                 P2_DMP(J,L)=0_1
	              END IF
	            END IF
	            P2_SGN(J,L)=SGN_P
	          END IF
	        END IF
!
! Damped commit: for a pure multiplicative period-2 mode the
! geometric midpoint of the flip bracket is the fixed point itself
! (F(x)/x* = x*/x), and in general a bounded step inside the
! bracket -- the relaxation path changes, the fixed point does not.
! Arithmetic fallback covers a zero endpoint.
!
	        IF(DO_P2)THEN
	          IF(P2_DMP(J,L) .GT. 0_1)THEN
	            IF(VOLD .GT. 0.0_LDP .AND. VNEW .GT. 0.0_LDP)THEN
	              FLD(J,L)=SQRT(VOLD*VNEW)
	            ELSE
	              FLD(J,L)=0.5_LDP*(VOLD+VNEW)
	            END IF
	            NP2_CNT=NP2_CNT+1
	          ELSE
	            FLD(J,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	          END IF
	        ELSE
	          FLD(J,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	        END IF
	      END DO
	    END DO
	  ELSE
!
! Inner boundary row I=NI: same options and accounting as the
! outward branch of MARCH.
!
	    DO L=1,NCF
	      IF(KN_INNER_METH .EQ. 'DIFFUSION' .AND.
	1         LS .LE. KN_NC)THEN
	        I_CORE=(PK_ETANI(L)+KN_DBB(L)*
	1         SQRT(KN_R(KN_ND)*KN_R(KN_ND)-KN_P(LS)*KN_P(LS))
	2         /KN_R(KN_ND))/PK_CHINI(L)
	      ELSE IF(KN_INNER_METH .EQ. 'ZERO_FLUX')THEN
	        I_CORE=IM_NI(L)
	      ELSE IF(KN_INNER_METH .EQ. 'PNT_SRCE')THEN
	        I_CORE=IM_NI(L)
	        IF(LS .LE. NC_PNT_SRCE_L)I_CORE=KN_IC(L)+IM_NI(L)
	      ELSE
	        I_CORE=KN_IC(L)
	      END IF
	      IF(LS .LE. KN_NC)THEN
	        VNEW=I_CORE
	      ELSE
	        VNEW=IM_NI(L)
	      END IF
	      IF(VNEW .LT. 0.0_LDP .OR. VNEW .NE. VNEW)THEN
	        VNEW=1.0E-30_LDP
	        NNEG=NNEG+1
	      END IF
	      IF(VNEW .GT. NONMONO_INEXT_CAP)THEN
	        VNEW=NONMONO_INEXT_CAP
	        NNEG=NNEG+1
	      END IF
	      VOLD=FLD(NI,L)
	      T1=ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAXF(L))
	      IF(T1 .GT. RES)THEN
	        WRC_I=NI
	        WRC_L=L
	      END IF
	      RES=MAX(RES,T1)
	      RES_FREQ(L)=MAX(RES_FREQ(L),T1)
	      IF(MAX(VOLD,ABS(VNEW)) .GT. 1.0E-4_LDP*VMAXF(L))THEN
	        NTOT_CUR=NTOT_CUR+1
	        IF(T1 .GT. KN_RAY_TOL)NBAD_CUR=NBAD_CUR+1
	        IF(DO_P2)THEN
	          SGN_P=0_1
	          IF(VNEW .GT. VOLD)SGN_P=1_1
	          IF(VNEW .LT. VOLD)SGN_P=-1_1
	          IF(P2_SGN(NI,L) .NE. 0_1 .AND. SGN_P .NE. 0_1 .AND.
	1           SGN_P .NE. P2_SGN(NI,L))THEN
	            IF(P2_FLIP(NI,L) .LT. 100_1)
	1             P2_FLIP(NI,L)=P2_FLIP(NI,L)+1_1
	            IF(INT(P2_FLIP(NI,L)) .GE. FG_ORD_P2N)
	1             P2_DMP(NI,L)=1_1
	          ELSE
	            P2_FLIP(NI,L)=0_1
	            IF(P2_DMP(NI,L) .GT. 0_1)THEN
	              P2_DMP(NI,L)=P2_DMP(NI,L)+1_1
	              IF(INT(P2_DMP(NI,L)) .GT. FG_ORD_P2REL)
	1               P2_DMP(NI,L)=0_1
	            END IF
	          END IF
	          P2_SGN(NI,L)=SGN_P
	        END IF
	      END IF
	      IF(DO_P2)THEN
	        IF(P2_DMP(NI,L) .GT. 0_1)THEN
	          IF(VOLD .GT. 0.0_LDP .AND. VNEW .GT. 0.0_LDP)THEN
	            FLD(NI,L)=SQRT(VOLD*VNEW)
	          ELSE
	            FLD(NI,L)=0.5_LDP*(VOLD+VNEW)
	          END IF
	          NP2_CNT=NP2_CNT+1
	        ELSE
	          FLD(NI,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	        END IF
	      ELSE
	        FLD(NI,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	      END IF
	    END DO
	    CALL KN_ORD_SPW_FILL(FLD,NI)
	    DO J=NI-1,1,-1
	      IF(J-1 .GE. 1)CALL KN_ORD_SPW_FILL(FLD,J-1)
	      IF(GAM(J,LS) .GE. 0.0_LDP)THEN
	        LB=1; LE=NCF; LDI=1
	      ELSE
	        LB=NCF; LE=1; LDI=-1
	      END IF
	      DO L=LB,LE,LDI
	        IF(HOTMODE .AND. (L .LT. W_LO .OR. L .GT. W_HI))CYCLE
	        CALL KN_ORD_SPV(FLD,J,L,SPJ)
	        SPW(L,MOD(J-1,4)+1)=SPJ
	        DSU=KN_ORD_DS(J+1,L)
	        DSD=KN_ORD_DS(J,L)
	        VNEW=FLD(J+1,L)*PK_A0(J,L)+
	1           (SPW(L,MOD(J,4)+1)*PK_A1(J,L)+
	2           SPJ*PK_A2(J,L)-
	3           DSU*PK_A3(J,L)-DSD*PK_A4(J,L))
	        IF(VNEW .LT. 0.0_LDP .OR. VNEW .NE. VNEW)THEN
	          VNEW=1.0E-30_LDP
	          NNEG=NNEG+1
	        END IF
	        IF(VNEW .GT. NONMONO_INEXT_CAP)THEN
	          VNEW=NONMONO_INEXT_CAP
	          NNEG=NNEG+1
	        END IF
	        VOLD=FLD(J,L)
	        T1=ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAXF(L))
	        IF(T1 .GT. RES)THEN
	          WRC_I=J
	          WRC_L=L
	        END IF
	        RES=MAX(RES,T1)
	        RES_FREQ(L)=MAX(RES_FREQ(L),T1)
	        IF(MAX(VOLD,ABS(VNEW)) .GT. 1.0E-4_LDP*VMAXF(L))THEN
	          NTOT_CUR=NTOT_CUR+1
	          KN_SGN_SUM=KN_SGN_SUM+
	1           (VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAXF(L))
	          IF(T1 .GT. KN_RAY_TOL)NBAD_CUR=NBAD_CUR+1
	          IF(PRINT_KNOP_DEBUG .AND. T1 .GT. 0.2_LDP .AND.
	1           NPO .GE. 8 .AND. KN_DUMP_N .LT. 2000)THEN
	            KN_DUMP_N=KN_DUMP_N+1
	            WRITE(LUER,'(A,I4,A,I4,A,I5,A,I5,4(A,ES11.4),3(A,I4))')
	1           ' KNDUMP p=',NPO,' r=',LS,' j=',J,' l=',L,
	2           ' res=',T1,' vo=',VOLD,' vn=',VNEW,
	3           ' gam=',GAM(J,LS),' jp1=',JPR1(J),' jp2=',JPR2(J),
	4           ' nfm=',NFM_L(J)
	          END IF
!
! P2 tracking (KN_ORD_P2DAMP): same per-cell sign-flip test as the
! inward march above.
!
	          IF(DO_P2)THEN
	            SGN_P=0_1
	            IF(VNEW .GT. VOLD)SGN_P=1_1
	            IF(VNEW .LT. VOLD)SGN_P=-1_1
	            IF(P2_SGN(J,L) .NE. 0_1 .AND. SGN_P .NE. 0_1 .AND.
	1             SGN_P .NE. P2_SGN(J,L))THEN
	              IF(P2_FLIP(J,L) .LT. 100_1)
	1               P2_FLIP(J,L)=P2_FLIP(J,L)+1_1
	              IF(INT(P2_FLIP(J,L)) .GE. FG_ORD_P2N)
	1               P2_DMP(J,L)=1_1
	            ELSE
	              P2_FLIP(J,L)=0_1
	              IF(P2_DMP(J,L) .GT. 0_1)THEN
	                P2_DMP(J,L)=P2_DMP(J,L)+1_1
	                IF(INT(P2_DMP(J,L)) .GT. FG_ORD_P2REL)
	1                 P2_DMP(J,L)=0_1
	              END IF
	            END IF
	            P2_SGN(J,L)=SGN_P
	          END IF
	        END IF
	        IF(DO_P2)THEN
	          IF(P2_DMP(J,L) .GT. 0_1)THEN
	            IF(VOLD .GT. 0.0_LDP .AND. VNEW .GT. 0.0_LDP)THEN
	              FLD(J,L)=SQRT(VOLD*VNEW)
	            ELSE
	              FLD(J,L)=0.5_LDP*(VOLD+VNEW)
	            END IF
	            NP2_CNT=NP2_CNT+1
	          ELSE
	            FLD(J,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	          END IF
	        ELSE
	          FLD(J,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	        END IF
	      END DO
	    END DO
	  END IF
!
	  NB_O=NBAD_CUR
	  NT_O=NTOT_CUR
!
! Newest history residual r_k = g(x_k)-x_k (the pass delta).  On a
! rejected-mix pass the pair is discarded with the history below.
!
	  IF(DO_AA .AND. .NOT. AA_OFF .AND. ALLOCATED(AA_XH))
	1    AA_RH(1:NI,1:NCF,1)=FLD(1:NI,1:NCF)-AA_XH(1:NI,1:NCF,1)
!
! Residual acceptance for the pending AA mix: if this pass did not
! improve on the pre-mix residual, the accelerated iterate was
! harmful.  Restore the pre-mix field (the slot-2 map output),
! clear the history and stop mixing this solve.  Checked before
! the divergence exit so a poisoned candidate is rejected instead
! of marking the ray diverged.
!
	  IF(AA_PEND)THEN
	    AA_PEND=.FALSE.
	    IF(RES .NE. RES .OR. RES .GT. AA_REF)THEN
	      FLD(1:NI,1:NCF)=AA_XH(1:NI,1:NCF,2)+AA_RH(1:NI,1:NCF,2)
	      RES=AA_REF
	      NH_AA=0
	      AA_OFF=.TRUE.
!$OMP CRITICAL (KN_AA_UPD)
	      KN_AA_REJ=KN_AA_REJ+1
!$OMP END CRITICAL (KN_AA_UPD)
!
! edited by Devin, 17-Oct-2026: NB_O/NT_O were latched from the
! rejected-mix pass above; the bulk-exit census must describe the
! restored iterate.  The slot-2 history pair is exactly the accepted
! step (x_pre = AA_XH(2), g(x_pre)-x_pre = AA_RH(2)), so rebuild the
! census with the same significance floor and tolerance the march
! used.  Only reachable under [KN_ORD_AA]=T (default F).
!
	      NBAD_CUR=0
	      NTOT_CUR=0
	      DO L=1,NCF
	        DO J=1,NI
	          T1=ABS(AA_RH(J,L,2))/
	1             MAX(ABS(AA_XH(J,L,2)),1.0E-4_LDP*VMAXF(L))
	          IF(MAX(ABS(AA_XH(J,L,2)),
	1             ABS(AA_XH(J,L,2)+AA_RH(J,L,2))) .GT.
	2             1.0E-4_LDP*VMAXF(L))THEN
	            NTOT_CUR=NTOT_CUR+1
	            IF(T1 .GT. KN_RAY_TOL)NBAD_CUR=NBAD_CUR+1
	          END IF
	        END DO
	      END DO
	      NB_O=NBAD_CUR
	      NT_O=NTOT_CUR
	    END IF
	  END IF
	  IF(ALLOCATED(KN_RES_PASS) .AND. .NOT. HOTMODE)THEN
	    IF(NPO .LE. SIZE(KN_RES_PASS))THEN
!$OMP CRITICAL (KN_RES_UPD)
	      KN_RES_PASS(NPO)=MAX(KN_RES_PASS(NPO),RES)
!$OMP END CRITICAL (KN_RES_UPD)
	    END IF
	  END IF
!
! Same divergence containment, relaxation adaptation and exits as
! SWEEP_FLD.  The residual history is now a single sequence (no
! direction parity), so same-direction bookkeeping reduces to the
! previous pass.
!
	  IF(RES .NE. RES .OR. RES .GT. HUGE(RES) .OR.
	1     (NPO .GE. 3 .AND. RES .GT. 1.0E+3_LDP*BEST_RES))THEN
	    DIVERGED=.TRUE.
	    IXIT=5
	    EXIT
	  END IF
!
! W3 (17-Oct-2026): best-pass census, same role as SWEEP_FLD --
! stall/cap exits carry the last-pass NB_O; the gate reads
! MIN(GS_NBAD, GS_NBAD_BLK).
!
	  IF(RES .LT. BEST_RES)THEN
	    BEST_RES=RES
	    BEST_NB=NBAD_CUR
	    BEST_NT=NTOT_CUR
	  END IF
	  IF(RES_PREV(1) .GT. 0.0_LDP)THEN
	    RATIO_OM=RES/RES_PREV(1)
!
! Oscillation detector: on a period-2 limit cycle (lambda ~ -1,
! observed on the compressing hot blocks at a fluctuating ~0.4
! residual plateau) the signed per-cell residual sum flips sign
! every pass.  A monotone marginal mode (lambda ~ +0.99) keeps its
! sign and is left to the stall/patience exits -- damping it would
! only slow a contracting tail.
!
	    IF(KN_SGN_PREV .NE. 0.0_LDP .AND.
	1     KN_SGN_SUM*KN_SGN_PREV .LT. 0.0_LDP)THEN
	      NOSC=NOSC+1
	    ELSE
	      NOSC=0
	    END IF
	    KN_SGN_PREV=KN_SGN_SUM
!
! edited by Devin, 14-Oct-2026: same seed-transient floor as
! SWEEP_FLD -- here the residual history is a single sequence, so
! the baseline is NPO-1 and the floor is NPO>=4 (a baseline from
! passes 1-2 is still the seed transient).  Applies to the
! growth-ratio and the sign-flip damping branches alike.
!
	    IF(RATIO_OM .GT. 1.05_LDP .AND. NPO .GE. 4)THEN
	      T4=OMEGA_RAY
	      OMEGA_RAY=MAX(0.45_LDP,OMEGA_RAY*0.7_LDP)
	    ELSE IF(RATIO_OM .LT. 0.5_LDP)THEN
	      T4=OMEGA_RAY
	      OMEGA_RAY=MIN(1.0_LDP,OMEGA_RAY*1.15_LDP)
	    ELSE IF(NOSC .GE. 2 .AND. NPO .GE. 4 .AND.
	1       OMEGA_RAY .GT. 0.45_LDP)THEN
!
! Consecutive sign flips in the flat band: a -1 mode under an
! under-relaxed update has |1-2*OMEGA|<1, so decaying OMEGA
! converges it to the same fixed point -- no clamping involved.
!
	      T4=OMEGA_RAY
	      OMEGA_RAY=MAX(0.45_LDP,OMEGA_RAY*0.7_LDP)
	    ELSE
	      T4=0.0_LDP
	    END IF
!
! An OMEGA_RAY change alters the relaxed map.  Stored pass deltas
! are omega*(F(x)-x), so for the affine refresh map a rescale to
! the new omega keeps the history pairs consistent instead of
! restarting it.  The same applies to the Aitken delta history
! (DLT_PRV, KN_ORD_AITKEN): an unscaled DEL_OLD/DEL_CUR pair across
! an omega change mixes two different maps, so the geometric-tail
! assumption behind the extrapolation only holds piecewise.
!
	    IF(T4 .GT. 0.0_LDP)THEN
	      IF(NH_AA .GT. 0 .AND. ALLOCATED(AA_XH))
	1        AA_RH(1:NI,1:NCF,1:NH_AA)=
	2        AA_RH(1:NI,1:NCF,1:NH_AA)*OMEGA_RAY/T4
	      IF(ALLOCATED(DLT_PRV))
	1        DLT_PRV(1:NI,1:NCF)=DLT_PRV(1:NI,1:NCF)*OMEGA_RAY/T4
	    END IF
	    IF(NPO .GE. 8)THEN
	      IF(RATIO_OM .GT. KN_RAY_STALL)THEN
	        NSTALL=NSTALL+1
	      ELSE
	        NSTALL=0
	      END IF
!
! W6 (17-Oct-2026, [FG_PLATEAU_DET], default T): ORD's residual is
! a single sequence, so the plateau test runs on consecutive
! passes: three ratios in [0.95,1.05] is flat, and two straight
! >0.95 ratios ending above 0.98 is the stagnation the same-
! direction-only counter can ride past.
!
	      IF(FG_PLATEAU_DET)THEN
	        IF(ABS(RATIO_OM-1.0_LDP) .LT. 0.05_LDP)THEN
	          NPLAT=NPLAT+1
	        ELSE
	          NPLAT=0
	        END IF
	        IF(NPLAT .GE. 3)THEN
	          NSTALL=MAX(NSTALL,2)
	          L_PLAT=.TRUE.
	        END IF
	        IF(RATIO_OM .GT. 0.98_LDP .AND.
	1         RATIO_PV .GT. 0.95_LDP)NSTALL=MAX(NSTALL,2)
	      END IF
	      RATIO_PV=RATIO_OM
	    END IF
	  END IF
	  RES_PREV(1)=RES
	  IF(WRC_I .EQ. WRC_PI .AND. WRC_L .EQ. WRC_PL .AND.
	1    WRC_I .GT. 0)THEN
	    WRC_SAME=WRC_SAME+1
	  ELSE
	    WRC_SAME=1
	  END IF
	  WRC_PI=WRC_I; WRC_PL=WRC_L
	  IF(.NOT. HAS_CP)THEN
!
!	    No compressing cells: pass 1 still carries the stale slope
!	    corner of the seed, so confirm on a second pass before a
!	    clean exit (see SWEEP_FLD for the legacy one-pass exit).
!
	    IF(NPO .GE. 2)THEN
	      CLEAN_EXIT=.TRUE.
	      IXIT=0
	      RES=0.0_LDP
	      RES_FREQ(1:NCF)=0.0_LDP
	      EXIT
	    END IF
	  END IF
	  IF(RES .LT. KN_RAY_TOL)THEN
	    CLEAN_EXIT=.TRUE.
	    IXIT=1
	    EXIT
	  END IF
! [FG_BULK_CONTR]: same ordering fix as SWEEP_FLD -- the bulk exit
! must not fire on a noncontracting pass sequence; RATIO_OM here is
! the consecutive-pass ratio.
!
	  IF(NPO .GE. 6 .AND. NT_O .GT. 0 .AND.
	1     NB_O .LE. MAX(50_8,NT_O/1000) .AND.
	2     (.NOT. FG_BULK_CONTR .OR. RATIO_OM .LE. KN_RAY_STALL))THEN
	    IXIT=2
	    EXIT
	  END IF
!
! W7 (17-Oct-2026): hot-pass extension, same rule as SWEEP_FLD --
! a stall/cap exit with residual left only inside a contiguous
! hot-frequency run continues restricted to that window.
!
	  IF(NSTALL .GE. 2)THEN
	    IF(.NOT. HOTMODE .AND. TRY_HOT())THEN
	      HOTMODE=.TRUE.
	    ELSE
	      IXIT=3
	      IF(L_PLAT)IXIT=6
	      EXIT
	    END IF
	  END IF
	  IF(NPO .GE. KN_RAY_MAXPASS)THEN
	    IF(.NOT. HOTMODE .AND. TRY_HOT())THEN
	      HOTMODE=.TRUE.
	    ELSE
	      IXIT=4
	      EXIT
	    END IF
	  END IF
!
! Aitken delta-squared on the pass iterate (KN_ORD_AITKEN).  Applied
! only to cells whose last pass still moved by more than the march
! tolerance (the above-tolerance overwrite criterion) and whose two
! consecutive deltas show local contraction |d_n|<|d_{n-1}| --
! expanding or sign-incoherent cells are skipped, and the
! extrapolated value is kept only if finite, positive and below the
! store cap.  Passes 1-2 just build the delta history (the seed ->
! pass-1 transient is not a usable contraction estimate); the
! committed iterate is the extrapolated one, so the SPW ring's
! lagged corner terms read a consistent lagged state on the next
! pass either way.
!
	  IF(DO_AIT)THEN
	    DO L=1,NCF
	      DO J=1,NI
	        DEL_CUR=FLD(J,L)-FLD_PRV(J,L)
	        DEL_OLD=DLT_PRV(J,L)
	        IF(NPO .GE. 3 .AND.
	1         ABS(FLD(J,L)) .GT. 1.0E-4_LDP*VMAXF(L) .AND.
	2         ABS(DEL_CUR) .GT. KN_RAY_TOL*ABS(FLD(J,L)) .AND.
	3         ABS(DEL_CUR) .LT. ABS(DEL_OLD) .AND.
	4         ABS(DEL_CUR-DEL_OLD) .GT.
	5         1.0E-30_LDP*VMAXF(L))THEN
	          XACC=FLD(J,L)-DEL_CUR*DEL_CUR/(DEL_CUR-DEL_OLD)
	          IF(XACC .GE. 0.0_LDP .AND. XACC .EQ. XACC .AND.
	1             XACC .LE. NONMONO_INEXT_CAP)FLD(J,L)=XACC
	        END IF
	        DLT_PRV(J,L)=DEL_CUR
	        FLD_PRV(J,L)=FLD(J,L)
	      END DO
	    END DO
	  END IF
!
! Safeguarded type-II Anderson mix (KN_ORD_AA).  Design follows
! the MOM refresh AA (mom_j_cmf_nonmono_v1.f): the Gram system is
! built on residuals weighted by the same normalisation as the
! march metric (1/MAX(|x|,1E-4*column peak)); cells outside the
! affine regime (weighted residual/secant > AA_GCAP) are excluded
! from the fit and keep the raw map output.  A small ridge keeps
! gamma bounded on near-dependent secants; the mixed correction is
! vetoed when non-finite or >1E3||f|| (a marginal mode needs a step
! of order ||f||/(1-lambda), so the MOM-style ||corr||>||f|| bound
! would reject every useful mix); the per-cell step is capped at
! 10|r| and non-positive outputs fall back to the map output.
! An applied mix sets AA_PEND; the next pass either confirms it
! (RES<=AA_REF) or rolls back above.
!
	  IF(DO_AA .AND. .NOT. AA_OFF .AND. ALLOCATED(AA_XH) .AND.
	1     NH_AA .GE. 2)THEN
	    MP_AA=NH_AA-1
	    GRM_AA=0.0_LDP
	    BRH_AA=0.0_LDP
	    NAA_SKIP=0_8
	    DO L=1,NCF
	      TCOL=MAX(MAXVAL(ABS(FLD(1:NI,L))),1.0E-28_LDP)
	      DO J=1,NI
	        WV=1.0_LDP/MAX(ABS(AA_XH(J,L,1)),1.0E-4_LDP*TCOL)
	        SC=AA_RH(J,L,1)*WV
	        T1=ABS(SC)
	        DO JJ=1,MP_AA
	          T1=MAX(T1,ABS(AA_RH(J,L,1)-AA_RH(J,L,JJ+1))*WV)
	        END DO
	        IF(T1 .GT. AA_GCAP)THEN
	          NAA_SKIP=NAA_SKIP+1
	          CYCLE
	        END IF
	        DO JJ=1,MP_AA
	          XV=(AA_RH(J,L,1)-AA_RH(J,L,JJ+1))*WV
	          BRH_AA(JJ)=BRH_AA(JJ)+XV*SC
	          DO KK=1,JJ
	            GRM_AA(KK,JJ)=GRM_AA(KK,JJ)+
	1             XV*(AA_RH(J,L,1)-AA_RH(J,L,KK+1))*WV
	          END DO
	        END DO
	      END DO
	    END DO
	    DO JJ=1,MP_AA
	      DO KK=JJ+1,MP_AA
	        GRM_AA(JJ,KK)=GRM_AA(KK,JJ)
	      END DO
	    END DO
	    DTR=0.0_LDP
	    DO JJ=1,MP_AA
	      DTR=MAX(DTR,GRM_AA(JJ,JJ))
	    END DO
	    DTR=1.0E-8_LDP*MAX(DTR,1.0E-30_LDP)
	    DO JJ=1,MP_AA
	      GRM_AA(JJ,JJ)=GRM_AA(JJ,JJ)+DTR
	    END DO
	    WGR_AA(1:MP_AA,1:MP_AA)=GRM_AA(1:MP_AA,1:MP_AA)
	    WBR_AA(1:MP_AA)=BRH_AA(1:MP_AA)
	    IF(INT(NI,8)*INT(NCF,8)-NAA_SKIP .LT.
	1       MAX(100_8,INT(NI,8)*INT(NCF,8)/1000))THEN
	      INFO_AA=-999
	    ELSE
	      CALL DGESV(MP_AA,1,WGR_AA,8,IPIV_AA,WBR_AA,8,INFO_AA)
	    END IF
	    IF(INFO_AA .NE. 0)THEN
!$OMP CRITICAL (KN_AA_UPD)
	      KN_AA_SNG=KN_AA_SNG+1
!$OMP END CRITICAL (KN_AA_UPD)
	    ELSE
	      GAM_AA(1:MP_AA)=WBR_AA(1:MP_AA)
	      FNRM=0.0_LDP
	      DNRM=0.0_LDP
	      DO L=1,NCF
	        TCOL=MAX(MAXVAL(ABS(FLD(1:NI,L))),1.0E-28_LDP)
	        DO J=1,NI
	          WV=1.0_LDP/MAX(ABS(AA_XH(J,L,1)),1.0E-4_LDP*TCOL)
	          IF(ABS(AA_RH(J,L,1))*WV .GT. AA_GCAP)CYCLE
	          FNRM=FNRM+(AA_RH(J,L,1)*WV)**2
	          SC=0.0_LDP
	          DO JJ=1,MP_AA
	            SC=SC+GAM_AA(JJ)*(FLD(J,L)-AA_XH(J,L,JJ+1)-
	1             AA_RH(J,L,JJ+1))*WV
	          END DO
	          DNRM=DNRM+SC*SC
	        END DO
	      END DO
	      FNRM=SQRT(FNRM)
	      DNRM=SQRT(DNRM)
!
! The MOM refresh vetoes ||corr||>||f||, but that bound is wrong for
! the regime this accelerator targets here: on a cell contracting at
! rate lambda the distance to the fixed point is |r|/(1-lambda), so
! the correct AA step exceeds ||f|| by ~1/(1-lambda).  The marginal
! KNOP rays measure lambda~0.98-0.99 (50-100x ||f||); the strict veto
! rejected nearly every mix on arm 12186580 (vetoed ~3x applied) and
! AA could never act on the modes it was built for.  The bound is
! still needed against ill-conditioned gamma garbage, so veto only
! non-finite or >1E3*||f|| corrections (lambda>0.999 -- beyond what
! KN_RAY_MAXPASS could converge anyway).  The 10|r| per-cell cap and
! the residual rollback remain the real safeguards.
!
	      IF(DNRM .NE. DNRM .OR. DNRM .GT. 1.0E+3_LDP*FNRM)THEN
!$OMP CRITICAL (KN_AA_UPD)
	        KN_AA_VETO=KN_AA_VETO+1
!$OMP END CRITICAL (KN_AA_UPD)
	      ELSE
	        DO L=1,NCF
	          TCOL=MAX(MAXVAL(ABS(FLD(1:NI,L))),1.0E-28_LDP)
	          DO J=1,NI
	            WV=1.0_LDP/MAX(ABS(AA_XH(J,L,1)),1.0E-4_LDP*TCOL)
	            IF(ABS(AA_RH(J,L,1))*WV .GT. AA_GCAP)CYCLE
	            XV=0.0_LDP
	            DO JJ=1,MP_AA
	              XV=XV+GAM_AA(JJ)*(FLD(J,L)-AA_XH(J,L,JJ+1)-
	1               AA_RH(J,L,JJ+1))
	            END DO
	            T4=FLD(J,L)
	            XV=T4-XV
	            T3=1.0E+1_LDP*MAX(ABS(AA_RH(J,L,1)),
	1             1.0E-4_LDP*TCOL)
	            STEP=T4-XV
	            IF(ABS(STEP) .GT. T3)XV=T4-SIGN(T3,STEP)
	            IF(XV .LE. 0.0_LDP .AND. T4 .GT. 0.0_LDP)XV=T4
	            IF(XV .GT. NONMONO_INEXT_CAP)XV=NONMONO_INEXT_CAP
	            FLD(J,L)=XV
	          END DO
	        END DO
	        AA_PEND=.TRUE.
	        AA_REF=RES
!$OMP CRITICAL (KN_AA_UPD)
	        KN_AA_APP=KN_AA_APP+1
!$OMP END CRITICAL (KN_AA_UPD)
	      END IF
	    END IF
	  END IF
	END DO
!
	DEALLOCATE(VMAXF)
	IF(DO_AIT)DEALLOCATE(FLD_PRV,DLT_PRV)
	IF(DO_P2)DEALLOCATE(P2_SGN,P2_FLIP,P2_DMP)
	IF(ALLOCATED(AA_XH))DEALLOCATE(AA_XH)
	IF(ALLOCATED(AA_RH))DEALLOCATE(AA_RH)
	IF(CLEAN_EXIT)NB_O=0
!
! W3 (17-Oct-2026): export the best-pass census -- the gate takes
! MIN(GS_NBAD, GS_NBAD_BLK) for stall/cap exits.
!
	NB_B=BEST_NB
	NT_B=BEST_NT
	IF(CLEAN_EXIT)NB_B=0_8
	RETURN
	END
!
! ----------------------------------------------------------------
! Rebuild the per-frequency discretisation on this ray: ray-point
! opacity/emissivity (same monotone-cubic kernels as the sweep),
! chi_hat/Q/source and the DTAU/E/A short-characteristics
! coefficients.  Mirrors the in-sweep INTEGRAL branch exactly.
! ----------------------------------------------------------------
	SUBROUTINE BUILD_FREQ(L)
	INTEGER L
!
! The branch and the reused scalar temporaries in the original single
! loop created an assumed output dependence that blocked vectorisation,
! so every EXP ran scalar (verified with -qopt-report); the identical
! in-sweep loop does vectorise.  Split into branch-free loops over
! T1V and the cubic evaluations.  Ray points coinciding with grid
! nodes (IGRD, built once per ray) then take the tabulated values
! exactly, as before.
!
	DO I=1,NI
	  T1V(I)=LOG(R_RAY(I,LS)/R_EXT(RAY_PNT(I,LS)))
	END DO
	DO I=1,NI
	  K=RAY_PNT(I,LS)
	  CHI_L(I)=EXP(((KN_CHI_COEF(K,1,L)*T1V(I)+
	1       KN_CHI_COEF(K,2,L))*T1V(I)+KN_CHI_COEF(K,3,L))*T1V(I)+
	1       KN_CHI_COEF(K,4,L))
	  ETA_L(I)=EXP(((KN_ETA_COEF(K,1,L)*T1V(I)+
	1       KN_ETA_COEF(K,2,L))*T1V(I)+KN_ETA_COEF(K,3,L))*T1V(I)+
	1       KN_ETA_COEF(K,4,L))
	  DCHR_L(I)=(3.0_LDP*KN_CHI_COEF(K,1,L)*T1V(I)+
	1       2.0_LDP*KN_CHI_COEF(K,2,L))*T1V(I)+KN_CHI_COEF(K,3,L)
	  DCHR_L(I)=DCHR_L(I)*CHI_L(I)/R_RAY(I,LS)
	END DO
	DO JPT=1,NGRD
	  I=IGRD(JPT)
	  K=RAY_PNT(I,LS)
	  CHI_L(I)=DBLE(KN_CHI_EXT(K,L))
	  ETA_L(I)=DBLE(KN_ETA_EXT(K,L))
	  DCHR_L(I)=KN_CHI_COEF(K,3,L)*CHI_L(I)/R_RAY(I,LS)
	END DO
!
! First frequency: expanding cells keep zero blue-edge inflow
! (Q_L=0, matching the in-sweep INIT branch).  Compressing cells
! instead keep their chi_hat/Q coupling to the lagged redder bin --
! the in-sweep INIT fix: forcing Q=0 there dropped the advective
! outflow from the bluest bin and left the refresh and sweep
! disagreeing at the band edge.  DLNM_L comes from MATCH_FREQ(L=1)
! so the matched span applies here too.
!
	IF(L .EQ. 1)THEN
	  CALL MATCH_FREQ(L)
	  DO I=1,NI
	    IF(.NOT.(ABS(CHI_L(I)) .LE. HUGE(CHI_L(I))) .OR.
	1      CHI_L(I) .LT. 1.0E-30_LDP)CHI_L(I)=1.0E-30_LDP
	    IF(.NOT.(ABS(ETA_L(I)) .LE. HUGE(ETA_L(I))))
	1      ETA_L(I)=0.0_LDP
	    IF(.NOT.(ABS(GAM(I,LS)) .LE. HUGE(GAM(I,LS))))
	1      GAM(I,LS)=0.0_LDP
	    IF(.NOT.(ABS(GAM_MAG(I,LS)) .LE. HUGE(GAM_MAG(I,LS))))
	1      GAM_MAG(I,LS)=0.0_LDP
	    IF(.NOT.(ABS(dGAMdR(I,LS)) .LE. HUGE(dGAMdR(I,LS))))
	1      dGAMdR(I,LS)=0.0_LDP
	  END DO
	  DO I=1,NI
	    Q_L(I)=0.0_LDP
	    IF(GAM(I,LS) .LT. 0.0_LDP)THEN
	      DLN=DLNM_L(I)
	      IF(.NOT.(DLN .GT. 0.0_LDP))DLN=1.0E-30_LDP
	      CHI_L(I)=CHI_L(I)+GAM_MAG(I,LS)/DLN
	      DCHR_L(I)=DCHR_L(I)-dGAMdR(I,LS)/DLN
	    END IF
	  END DO
	  DO I=1,NI
	    IF(.NOT.(ABS(CHI_L(I)) .LE. HUGE(CHI_L(I))) .OR.
	1      CHI_L(I) .LT. 1.0E-30_LDP)CHI_L(I)=1.0E-30_LDP
	    IF(.NOT.(ABS(DCHR_L(I)) .LE. HUGE(DCHR_L(I))))
	1      DCHR_L(I)=0.0_LDP
	  END DO
	  DO I=1,NI
	    IF(GAM(I,LS) .LT. 0.0_LDP)THEN
	      DLN=DLNM_L(I)
	      IF(.NOT.(DLN .GT. 0.0_LDP))DLN=1.0E-30_LDP
	      Q_L(I)=GAM(I,LS)/(DLN*CHI_L(I))
	    END IF
	    SR_L(I)=ETA_L(I)/CHI_L(I)
	  END DO
	  DO I=1,NI
	    IF(.NOT.(ABS(Q_L(I)) .LE. HUGE(Q_L(I))))Q_L(I)=0.0_LDP
	    IF(.NOT.(ABS(SR_L(I)) .LE. HUGE(SR_L(I))))SR_L(I)=0.0_LDP
	  END DO
	ELSE
!
! The matched compressing slot/span (offset 1 and KN_DLN_N(L) when
! DO_NONMONO_FREQ_MATCH is off) is filled per frequency before the
! coefficient loops read it.
!
	  CALL MATCH_FREQ(L)
	  DO I=1,NI
	    DLN=DLNM_L(I)
	    Q_L(I)=GAM_MAG(I,LS)/DLN
	    CHI_L(I)=CHI_L(I)+Q_L(I)
	  END DO
	  DO I=1,NI
	    IF(CHI_L(I) .NE. CHI_L(I) .OR. CHI_L(I) .LT. 1.0E-30_LDP)
	1      CHI_L(I)=1.0E-30_LDP
	  END DO
	  DO I=1,NI
	    DLN=DLNM_L(I)
	    DCHR_L(I)=DCHR_L(I)+
	1      SIGN(1.0_LDP,GAM(I,LS))*dGAMdR(I,LS)/DLN
	  END DO
	  DO I=1,NI
	    IF(ETA_L(I) .NE. ETA_L(I))ETA_L(I)=0.0_LDP
	    IF(GAM(I,LS) .NE. GAM(I,LS))GAM(I,LS)=0.0_LDP
	    IF(GAM_MAG(I,LS) .NE. GAM_MAG(I,LS))GAM_MAG(I,LS)=0.0_LDP
	    IF(dGAMdR(I,LS) .NE. dGAMdR(I,LS))dGAMdR(I,LS)=0.0_LDP
	    IF(DCHR_L(I) .NE. DCHR_L(I))DCHR_L(I)=0.0_LDP
	  END DO
	  DO I=1,NI
	    DLN=DLNM_L(I)
	    Q_L(I)=GAM(I,LS)/(DLN*CHI_L(I))
	    SR_L(I)=ETA_L(I)/CHI_L(I)
	  END DO
	  DO I=1,NI
	    IF(Q_L(I) .NE. Q_L(I))Q_L(I)=0.0_LDP
	    IF(SR_L(I) .NE. SR_L(I))SR_L(I)=0.0_LDP
	  END DO
	END IF
!
	IF(KN_METHOD_SAV .EQ. 'ZERO')THEN
	  DO I=1,NI-1
	    DZ=Z(I,LS)-Z(I+1,LS)
	    DTAU_L(I)=0.5_LDP*(CHI_L(I)+CHI_L(I+1))*DZ
	  END DO
	ELSE
	  DO I=1,NI-1
	    DZ=Z(I,LS)-Z(I+1,LS)
	    DTAU_L(I)=0.5_LDP*DZ*( CHI_L(I)+CHI_L(I+1) +
	1      DZ*( DCHR_L(I+1)*Z(I+1,LS)/R_RAY(I+1,LS) -
	1      DCHR_L(I)*Z(I,LS)/R_RAY(I,LS) )/6.0_LDP )
	    IF(CHI_L(I) .LT. CHI_L(I+1))THEN
	      DTAU_L(I)=MAX(CHI_L(I)*DZ,DTAU_L(I))
	      DTAU_L(I)=MIN(CHI_L(I+1)*DZ,DTAU_L(I))
	    ELSE
	      DTAU_L(I)=MIN(CHI_L(I)*DZ,DTAU_L(I))
	      DTAU_L(I)=MAX(CHI_L(I+1)*DZ,DTAU_L(I))
	    END IF
	  END DO
	END IF
	DO I=1,NI-1
	  IF(DTAU_L(I) .NE. DTAU_L(I) .OR. DTAU_L(I) .LT. 0.0_LDP)
	1    DTAU_L(I)=0.5_LDP*(CHI_L(I)+CHI_L(I+1))*
	1       (Z(I,LS)-Z(I+1,LS))
	END DO
!
	DO I=1,NI-1
	  T1=DTAU_L(I)
	  EEL(I)=0.0_LDP
	  IF(T1 .LT. 700.0_LDP)EEL(I)=EXP(-T1)
	  IF(T1 .GE. 40.0_LDP)THEN
	  ELSE IF(T1 .GT. 0.5_LDP)THEN
	    E0L(I)=1.0_LDP-EEL(I)
	    E1L(I)=1.0_LDP-E0L(I)/T1
	    E2L(I)=1.0_LDP-2.0_LDP*E1L(I)/T1
	    E3L(I)=1.0_LDP-3.0_LDP*E2L(I)/T1
	  ELSE IF(T1 .GT. 0.1_LDP)THEN
	    E3L(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1         (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	2         (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP*
	3         (1.0_LDP-T1/10.0_LDP*(1.0_LDP-T1/11.0_LDP*
	4         (1.0_LDP-T1/12.0_LDP*(1.0_LDP-T1/13.0_LDP)))))))) )
	    E2L(I)=T1*( 1.0_LDP-E3L(I) )/3.0_LDP
	    E1L(I)=T1*( 1.0_LDP-E2L(I) )/2.0_LDP
	    E0L(I)=T1*( 1.0_LDP-E1L(I) )
	  ELSE
	    E3L(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1         (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	2         (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP) ))))
	    E2L(I)=T1*( 1.0_LDP-E3L(I) )/3.0_LDP
	    E1L(I)=T1*( 1.0_LDP-E2L(I) )/2.0_LDP
	    E0L(I)=T1*( 1.0_LDP-E1L(I) )
	  END IF
	END DO
!
	DO I=1,NI-1
	  T1=DTAU_L(I)
	  IF(T1 .GE. 40.0_LDP)THEN
	    A0L(I)=EEL(I)
	    A1L(I)=(6.0_LDP-12.0_LDP/T1)/T1/T1
	    A2L(I)=1.0_LDP-A1L(I)
	    A3L(I)=(2.0_LDP-6.0_LDP/T1)/T1
	    A4L(I)=(4.0_LDP-6.0_LDP/T1)/T1-1.0_LDP
	  ELSE
	    A0L(I)=EEL(I)
	    A1L(I)=E0L(I)-3.0_LDP*E2L(I)+2.0_LDP*E3L(I)
	    A2L(I)=3.0_LDP*E2L(I)-2.0_LDP*E3L(I)
	    A3L(I)=DTAU_L(I)*(E1L(I)-2.0_LDP*E2L(I)+E3L(I))
	    A4L(I)=DTAU_L(I)*(E3L(I)-E2L(I))
	  END IF
	END DO
	RETURN
	END
!
! ----------------------------------------------------------------
! Per-point frequency-derivative spans for frequency L.  Expanding
! cells get the backward span KN_DLN_P(L) = ln(nu_{L-1}/nu_L), the
! same dLOG_NU the in-sweep solve uses for GAM>=0, so the refresh
! coefficients reproduce the sweep exactly; compressing cells get
! the forward span (matched, or +1 bin) below.
! ----------------------------------------------------------------
	SUBROUTINE MATCH_FREQ(L)
!
! DO_NONMONO_FREQ_MATCH: remap the compressing upwind slot from the
! adjacent bin to the bin containing the characteristic foot and
! widen the frequency-derivative span accordingly (Ceulemans et
! al. 2025).  NFM_L = slot offset (capped on-grid, >=1);
! DLNM_L = span ln(nu_L/nu_{L+NFM}) via the captured KN_FREQ grid,
! falling back to NFM*KN_DLN_N(L).  With the flag off compressing
! points keep the original +1 / KN_DLN_N(L) treatment.
!
	INTEGER L,NNF
	DO I=1,NI
	  NFM_L(I)=1
	  DLNM_L(I)=KN_DLN_P(L)
	END DO
	IF(L .GE. NCF)RETURN
	IF(FG_FREQ_MATCH .AND. KN_DLN_N(L) .GT. 0.0_LDP)THEN
	  DO I=1,NI
	    IF(GAM(I,LS) .LT. 0.0_LDP)THEN
	      T1=DRF_L(I)/KN_DLN_N(L)
	      IF(T1 .GE. DBLE(NCF-L))THEN
	        NNF=NCF-L
	      ELSE
	        NNF=MAX(1,CEILING(T1))
	      END IF
	      NFM_L(I)=NNF
	      DLNM_L(I)=NNF*KN_DLN_N(L)
	      IF(ALLOCATED(KN_FREQ))THEN
	        IF(KN_FREQ(L+NNF) .GT. 0.0_LDP .AND.
	1           KN_FREQ(L+NNF) .LT. KN_FREQ(L))
	2          DLNM_L(I)=LOG(KN_FREQ(L)/KN_FREQ(L+NNF))
	      END IF
	    END IF
	  END DO
	ELSE
	  DO I=1,NI
	    IF(GAM(I,LS) .LT. 0.0_LDP)DLNM_L(I)=KN_DLN_N(L)
	  END DO
	END IF
	RETURN
	END
!
! ----------------------------------------------------------------
! Modified source function S' for frequency L.  Expanding cells
! (GAM>=0) take the upwind at L-1 (fresh on ascending passes),
! compressing cells take the spatially-interpolated upwind at
! L+NFM_L(I) (matched slot; fresh on descending passes), matching
! the in-sweep I_NEXT read.
! ----------------------------------------------------------------
	SUBROUTINE SRC_PRIME(FLD,L)
	REAL(KIND=LDP) FLD(NI,NCF)
	INTEGER L
	DO I=1,NI
	  IF(GAM(I,LS) .GE. 0.0_LDP)THEN
	    INEXT_V=0.0_LDP
	    IF(L .GT. 1)INEXT_V=FLD(I,L-1)
	    SP_L(I)=SR_L(I)+Q_L(I)*INEXT_V
	  ELSE IF(L+NFM_L(I) .LE. NCF)THEN
	    IF(IKC1(I) .EQ. 1 .AND. IKC2(I) .EQ. 1)THEN
	      INEXT_V=(1.0_LDP-WINX(I))*FLD(JPR1(I),L+NFM_L(I))+
	1              WINX(I)*FLD(JPR2(I),L+NFM_L(I))
	    ELSE IF(IKC1(I) .EQ. 1)THEN
	      INEXT_V=FLD(JPR1(I),L+NFM_L(I))
	    ELSE IF(IKC2(I) .EQ. 1)THEN
	      INEXT_V=FLD(JPR2(I),L+NFM_L(I))
	    ELSE
	      INEXT_V=0.0_LDP
	    END IF
	    SP_L(I)=SR_L(I)+ABS(Q_L(I))*INEXT_V
	  ELSE
	    SP_L(I)=SR_L(I)
	  END IF
	END DO
	RETURN
	END
!
! ----------------------------------------------------------------
! Monotonic-parabola slopes of S' in optical depth (identical to the
! in-sweep construction, including the monotonicity clamps).
! ----------------------------------------------------------------
	SUBROUTINE SLOPES()
	DO I=1,NI-1
	  IF(DTAU_L(I) .GT. 1.0E-30_LDP)THEN
	    S_L(I)=(SP_L(I+1)-SP_L(I))/DTAU_L(I)
	  ELSE
	    S_L(I)=0.0_LDP
	  END IF
	END DO
!
	S_L(NI)=0.0_LDP
	T1_DTAU=DTAU_L(1)+DTAU_L(2)
	IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	  DS_L(1)=0.0_LDP
	ELSE
	  DS_L(1)=S_L(1)+(S_L(1)-S_L(2))*DTAU_L(1)/T1_DTAU
	END IF
	DO I=2,NI-1
	  T1_DTAU=DTAU_L(I-1)+DTAU_L(I)
	  IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	    DS_L(I)=0.0_LDP
	  ELSE
	    DS_L(I)=(S_L(I-1)*DTAU_L(I)+S_L(I)*DTAU_L(I-1))/T1_DTAU
	  END IF
	END DO
	IF(NI .GT. 2)THEN
	  T1_DTAU=DTAU_L(NI-2)+DTAU_L(NI-1)
	  IF(ABS(T1_DTAU) .LT. 1.0E-30_LDP)THEN
	    DS_L(NI)=0.0_LDP
	  ELSE
	    DS_L(NI)=S_L(NI-1)+
	1        (S_L(NI-1)-S_L(NI-2))*DTAU_L(NI-1)/T1_DTAU
	  END IF
	ELSE
	  DS_L(NI)=0.0_LDP
	END IF
!
	DS_L(1)=( SIGN(1.0_LDP,S_L(1))+SIGN(1.0_LDP,DS_L(1)) )*
	1          MIN(ABS(S_L(1)),0.5_LDP*ABS(DS_L(1)))
	DO I=2,NI-1
	  DS_L(I)=( SIGN(1.0_LDP,S_L(I-1))+SIGN(1.0_LDP,S_L(I)) )*
	1          MIN(ABS(S_L(I-1)),ABS(S_L(I)),0.5_LDP*ABS(DS_L(I)))
	END DO
	DS_L(NI)=( SIGN(1.0_LDP,S_L(NI-1))+SIGN(1.0_LDP,DS_L(NI)) )*
	1          MIN(ABS(S_L(NI-1)),0.5_LDP*ABS(DS_L(NI)))
	RETURN
	END
!
! ----------------------------------------------------------------
! Short-characteristics march for frequency L, tracking the maximum
! relative overwrite change into RES.  The 1e-25 metric floor keeps
! cells with negligible intensity (< ~1e-25) from dominating the
! convergence test; the field values themselves are unchanged.
! ----------------------------------------------------------------
	SUBROUTINE MARCH(FLD,L,INWARD,RES)
	REAL(KIND=LDP) FLD(NI,NCF),RES
	INTEGER L
	LOGICAL INWARD
	REAL(KIND=LDP) VMAX_L,RESL
! Residual is measured relative to the larger of the cell's previous
! value and 1e-4 of the frequency-column peak: cells below ~0.01% of
! the column maximum carry no measurable flux and would otherwise pin
! the max-norm (and hence the maturity gate) at irrelevant floor
! oscillations.  The bad-cell census uses this same normalized
! residual, so its fraction matches the max-norm tolerance.  RESL
! mirrors the max restricted to this frequency for block-GS diagnostics.
	VMAX_L=MAX(MAXVAL(FLD(1:NI,L)),1.0E-25_LDP)
	RESL=0.0_LDP
	IF(INWARD)THEN
	  VNEW=0.0_LDP
	  VOLD=FLD(1,L)
	  IF(ABS(VNEW-VOLD) .GT. 1.0E-25_LDP*VOLD)THEN
	    DCELL=ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAX_L)
	    IF(DCELL .GT. RES)THEN
	      WRC_I=1
	      WRC_L=L
	    END IF
	    RES=MAX(RES,DCELL)
	    RESL=MAX(RESL,DCELL)
	  END IF
	  FLD(1,L)=0.0_LDP
	  DO I=1,NI-1
	    VNEW=FLD(I,L)*A0L(I)+(SP_L(I)*A1L(I)+SP_L(I+1)*A2L(I)
	1         +DS_L(I)*A3L(I)+DS_L(I+1)*A4L(I))
	    IF(VNEW .LT. 0.0_LDP .OR. VNEW .NE. VNEW)THEN
	      VNEW=1.0E-30_LDP
	      NNEG=NNEG+1
	    END IF
	    IF(VNEW .GT. NONMONO_INEXT_CAP)THEN
	      VNEW=NONMONO_INEXT_CAP
	      NNEG=NNEG+1
	    END IF
	    VOLD=FLD(I+1,L)
	    DCELL=ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAX_L)
	    IF(DCELL .GT. RES)THEN
	      WRC_I=I+1
	      WRC_L=L
	    END IF
	    RES=MAX(RES,DCELL)
	    RESL=MAX(RESL,DCELL)
	    IF(MAX(VOLD,ABS(VNEW)) .GT. 1.0E-4_LDP*VMAX_L)THEN
	      NTOT_CUR=NTOT_CUR+1
!
! Signed per-cell update for the W6 sign-flip detector in
! SWEEP_FLD (same metric the ORD loop accumulates inline).
!
	      SGN_ACC=SGN_ACC+DCELL
	      IF(ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAX_L)
	1       .GT. KN_RAY_TOL)NBAD_CUR=NBAD_CUR+1
	    END IF
	    FLD(I+1,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	  END DO
	ELSE
	  IF(KN_INNER_METH .EQ. 'DIFFUSION' .AND. LS .LE. KN_NC)THEN
	    I_CORE=(ETA_L(NI)+KN_DBB(L)*
	1       SQRT(KN_R(KN_ND)*KN_R(KN_ND)-KN_P(LS)*KN_P(LS))
	1       /KN_R(KN_ND))/CHI_L(NI)
	  ELSE IF(KN_INNER_METH .EQ. 'ZERO_FLUX')THEN
	    I_CORE=IM_NI(L)
	  ELSE IF(KN_INNER_METH .EQ. 'PNT_SRCE')THEN
	    I_CORE=IM_NI(L)
	    IF(LS .LE. NC_PNT_SRCE_L)I_CORE=KN_IC(L)+IM_NI(L)
	  ELSE
	    I_CORE=KN_IC(L)
	  END IF
	  IF(LS .LE. KN_NC)THEN
	    VNEW=I_CORE
	  ELSE
	    VNEW=IM_NI(L)
	  END IF
	  IF(VNEW .LT. 0.0_LDP .OR. VNEW .NE. VNEW)THEN
	    VNEW=1.0E-30_LDP
	    NNEG=NNEG+1
	  END IF
	  IF(VNEW .GT. NONMONO_INEXT_CAP)THEN
	    VNEW=NONMONO_INEXT_CAP
	    NNEG=NNEG+1
	  END IF
	  VOLD=FLD(NI,L)
	  DCELL=ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAX_L)
	  IF(DCELL .GT. RES)THEN
	    WRC_I=NI
	    WRC_L=L
	  END IF
	  RES=MAX(RES,DCELL)
	  RESL=MAX(RESL,DCELL)
	  IF(MAX(VOLD,ABS(VNEW)) .GT. 1.0E-4_LDP*VMAX_L)THEN
	    NTOT_CUR=NTOT_CUR+1
	    IF(DCELL .GT. KN_RAY_TOL)NBAD_CUR=NBAD_CUR+1
	  END IF
	  FLD(NI,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	  DO I=NI-1,1,-1
	    VNEW=FLD(I+1,L)*A0L(I)+(SP_L(I+1)*A1L(I)+SP_L(I)*A2L(I)
	1         -DS_L(I+1)*A3L(I)-DS_L(I)*A4L(I))
	    IF(VNEW .LT. 0.0_LDP .OR. VNEW .NE. VNEW)THEN
	      VNEW=1.0E-30_LDP
	      NNEG=NNEG+1
	    END IF
	    IF(VNEW .GT. NONMONO_INEXT_CAP)THEN
	      VNEW=NONMONO_INEXT_CAP
	      NNEG=NNEG+1
	    END IF
	    VOLD=FLD(I,L)
	    DCELL=ABS(VNEW-VOLD)/MAX(VOLD,1.0E-4_LDP*VMAX_L)
	    IF(DCELL .GT. RES)THEN
	      WRC_I=I
	      WRC_L=L
	    END IF
	    RES=MAX(RES,DCELL)
	    RESL=MAX(RESL,DCELL)
	    IF(MAX(VOLD,ABS(VNEW)) .GT. 1.0E-4_LDP*VMAX_L)THEN
	      NTOT_CUR=NTOT_CUR+1
	      SGN_ACC=SGN_ACC+DCELL
	      IF(DCELL .GT. KN_RAY_TOL)NBAD_CUR=NBAD_CUR+1
	    END IF
	    FLD(I,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)
	  END DO
	END IF
	RES_FREQ(L)=RESL
	RETURN
	END
!
! ----------------------------------------------------------------
! Write the converged ray field back to the lagged depth store with
! the same guards the in-sweep write uses, but an UNRELAXED commit
! (W4, 17-Oct-2026): the inner march's OMEGA_RAY is the only
! under-relaxation on the refresh path.  NSTORE counts cells whose
! merged change exceeded 50% (same semantics as the in-sweep
! FG_GS_NCONV); NCAP counts the cap/NaN floor events (same
! semantics as FG_INEXT_CAPPED).
! ----------------------------------------------------------------
	SUBROUTINE WRITE_STORE(FLD,STORE,DMAX)
	REAL(KIND=LDP) FLD(NI,NCF)
	REAL(KIND=LDP) STORE(:,:,:)
	REAL(KIND=LDP) DMAX
	INTEGER KS,LS2,JP2
	REAL(KIND=LDP) V1,V2,DLOC
	REAL(KIND=LDP) VMAX_F(NCF)
!
! Column peak per frequency for the >50% churn count.  Without this
! floor a store cell sitting at the 1E-30 write floor counts as
! churned whenever the refreshed field differs at all (|V1-V2|/V2 is
! then >0.5 for any V1>1.5E-30), so the count never decays to zero
! and the maturity gate in SOLVE_FOR_POPS stays closed even when
! every significant cell is converged.  The 1E-4*peak significance
! test is the same one MARCH applies to the residual metrics
! (VMAX_L), so NSTORE and NTOT now cover the same cell population.
! edited by Devin, 28-Sep-2026
!
	DO LS2=1,NCF
	  VMAX_F(LS2)=MAX(1.0E-25_LDP,
	1    MAXVAL(FLD(1:NI,LS2),
	2    MASK=FLD(1:NI,LS2) .EQ. FLD(1:NI,LS2)))
	END DO
	DMAX=0.0_LDP
	DO KS=1,KN_ND
	  KC=COMP_DEPTH_IDX(KS)
	  IF(KC .GE. 1)THEN
	    JP2=J_PNT(KS,LS)
	    IF(JP2 .GE. 1 .AND. JP2 .LE. NI)THEN
! The in-sweep store holds frequency l in slot l-1; the refreshed
! frequency-1 field has no store slot and is not written back.
	      DO LS2=2,NCF
	        V1=FLD(JP2,LS2)
	        IF(V1 .EQ. V1 .AND. V1 .GT. 0.0_LDP)THEN
	          IF(V1 .GT. NONMONO_INEXT_CAP)THEN
	            V1=NONMONO_INEXT_CAP
	            NCAP=NCAP+1
	          END IF
	        ELSE
	          V1=1.0E-30_LDP
	          IF(.NOT. (FLD(JP2,LS2) .EQ. FLD(JP2,LS2)))
	1            NCAP=NCAP+1
	        END IF
	        V2=STORE(KC,LS,LS2-1)
	        IF(V2 .GT. 0.0_LDP)THEN
	          IF(ABS(V1-V2)/V2 .GT. 0.5_LDP .AND.
	1          MAX(V1,V2) .GT. 1.0E-4_LDP*VMAX_F(LS2))
	2          NSTORE=NSTORE+1
	        ELSE IF(V1 .GT. 1.0E-4_LDP*VMAX_F(LS2))THEN
	          NSTORE=NSTORE+1
	        END IF
! Retained-field delta uses the committed value, not the raw
! solve output -- the gate reads the store, so the metric must
! describe what the store actually holds.  A cell the seed left
! at the floor (or non-finite, healed above) counts 1.0 when the
! refresh writes a significant value.
!
! W4 (17-Oct-2026): unrelaxed commit.  The inner march already
! under-relaxes via OMEGA_RAY; applying INEXT_RELAX here double-
! damped every slot (effective motion ~0.22 of the map with the
! audited omega=0.5) and locked the marginal |Q|~1 modes at
! per-pass contraction ~0.995.  The generation-level INEXT_RELAX
! stays on the sweep write sites only.
!
	        IF(V2 .GT. 0.0_LDP)THEN
	          DLOC=ABS(V1-V2)/MAX(V2,1.0E-4_LDP*VMAX_F(LS2))
	          STORE(KC,LS,LS2-1)=V1
	        ELSE
	          DLOC=0.0_LDP
	          IF(V1 .GT. 1.0E-4_LDP*VMAX_F(LS2))DLOC=1.0_LDP
	          STORE(KC,LS,LS2-1)=V1
	        END IF
	        DMAX=MAX(DMAX,DLOC)
	      END DO
	    END IF
	  END IF
	END DO
	RETURN
	END
!
! ----------------------------------------------------------------
! Capped provisional commit (17-Oct-2026, [FG_COMMIT_NC] under
! [FG_COMMIT_CONV]): same write-back as WRITE_STORE, but each cell's
! committed change is bounded by FG_COMMIT_NC times the retained
! value (significance floored at 1E-4 of the frequency column peak,
! matching DLOC's own normalization).  A seed fill (V2 at/below the
! write floor) may inject at most NC times that significance floor
! -- an unconverged ray cannot drop a large value into an empty
! slot.  NSTORE keeps its >50%-moved semantics; DMAX reports the
! capped delta actually committed.
! ----------------------------------------------------------------
	SUBROUTINE WRITE_STORE_CAP(FLD,STORE,DMAX)
	REAL(KIND=LDP) FLD(NI,NCF)
	REAL(KIND=LDP) STORE(:,:,:)
	REAL(KIND=LDP) DMAX
	INTEGER KS,LS2,JP2
	REAL(KIND=LDP) V1,V2,DLOC,DVB
	REAL(KIND=LDP) VMAX_F(NCF)
!
	DO LS2=1,NCF
	  VMAX_F(LS2)=MAX(1.0E-25_LDP,
	1    MAXVAL(FLD(1:NI,LS2),
	2    MASK=FLD(1:NI,LS2) .EQ. FLD(1:NI,LS2)))
	END DO
	DMAX=0.0_LDP
	DO KS=1,KN_ND
	  KC=COMP_DEPTH_IDX(KS)
	  IF(KC .GE. 1)THEN
	    JP2=J_PNT(KS,LS)
	    IF(JP2 .GE. 1 .AND. JP2 .LE. NI)THEN
	      DO LS2=2,NCF
	        V1=FLD(JP2,LS2)
	        IF(V1 .EQ. V1 .AND. V1 .GT. 0.0_LDP)THEN
	          IF(V1 .GT. NONMONO_INEXT_CAP)THEN
	            V1=NONMONO_INEXT_CAP
	            NCAP=NCAP+1
	          END IF
	        ELSE
	          V1=1.0E-30_LDP
	          IF(.NOT. (FLD(JP2,LS2) .EQ. FLD(JP2,LS2)))
	1            NCAP=NCAP+1
	        END IF
	        V2=STORE(KC,LS,LS2-1)
!
! Cap the committed change at FG_COMMIT_NC times the retained value
! (or the significance floor when the slot held only a seed).  The
! provisional write moves the store toward the refreshed field
! without letting an unconverged iterate ratchet it by an unbounded
! factor per refresh.
!
	        DVB=FG_COMMIT_NC*MAX(V2,1.0E-4_LDP*VMAX_F(LS2))
	        IF(ABS(V1-V2) .GT. DVB)V1=V2+SIGN(DVB,V1-V2)
!
	        IF(V2 .GT. 0.0_LDP)THEN
	          IF(ABS(V1-V2)/V2 .GT. 0.5_LDP .AND.
	1          MAX(V1,V2) .GT. 1.0E-4_LDP*VMAX_F(LS2))
	2          NSTORE=NSTORE+1
	        ELSE IF(V1 .GT. 1.0E-4_LDP*VMAX_F(LS2))THEN
	          NSTORE=NSTORE+1
	        END IF
	        IF(V2 .GT. 0.0_LDP)THEN
	          DLOC=ABS(V1-V2)/MAX(V2,1.0E-4_LDP*VMAX_F(LS2))
	          STORE(KC,LS,LS2-1)=V1
	        ELSE
	          DLOC=0.0_LDP
	          IF(V1 .GT. 1.0E-4_LDP*VMAX_F(LS2))DLOC=1.0_LDP
	          STORE(KC,LS,LS2-1)=V1
	        END IF
	        DMAX=MAX(DMAX,DLOC)
	      END DO
	    END IF
	  END IF
	END DO
	RETURN
	END
!
! ----------------------------------------------------------------
! W7 (17-Oct-2026): hot-pass extension test.  Returns .TRUE. and
! sets the host-scope W_LO/W_HI window when the last-pass residual
! field is confined to a contiguous run of above-tolerance
! frequencies, the per-ray KN_RAY_HOTPASSES budget remains, and the
! refresh-wide KN_HOT_TOTAL_PASSES cap has not been reached.  The
! run window is extended +/-2 frequencies to cover the upwind reads
! of the marginal cells.
! ----------------------------------------------------------------
	LOGICAL FUNCTION TRY_HOT()
	INTEGER L_H
	TRY_HOT=.FALSE.
	IF(KN_RAY_HOTPASSES .LE. 0)RETURN
	IF(KN_HOT_TOTAL_PASSES .LE. 0)RETURN
	IF(KN_HOT_SPENT .GE. KN_HOT_TOTAL_PASSES)RETURN
	W_LO=NCF
	W_HI=1
	DO L_H=1,NCF
	  IF(RES_FREQ(L_H) .GT. KN_RAY_TOL)THEN
	    W_LO=MIN(W_LO,L_H)
	    W_HI=MAX(W_HI,L_H)
	  END IF
	END DO
	IF(W_HI .LT. W_LO)RETURN
	W_LO=MAX(1,W_LO-2)
	W_HI=MIN(NCF,W_HI+2)
	TRY_HOT=.TRUE.
	RETURN
	END
!
! ----------------------------------------------------------------
! W9 (17-Oct-2026): max committed-write jump for the FG_COMMIT_JUMP
! guard.  Measures max |V1-V2|/MAX(V2,column peak/FG_COMMIT_JUMP)
! over the cells the store already holds (V2>0 and finite).  Fills
! into never-committed frequency columns are bounded by the
! sibling-ray peak at the same field frequency (VMAX_X); only a
! column that is at the floor for every ray (true cold start) is
! exempt, so a fresh seed can still land.  The guard bounds the
! committed-near-blowup class (audit: retained delta 4.19E-1 ->
! 5.03E3) while allowing genuine factor-of-few rewrites.
! ----------------------------------------------------------------
	REAL(KIND=LDP) FUNCTION FG_JUMP_MAX(FLD,STORE)
	REAL(KIND=LDP) FLD(NI,NCF)
	REAL(KIND=LDP) STORE(:,:,:)
	INTEGER KS,KS2,LS2,JP2,RS,KC2
	REAL(KIND=LDP) V1,V2,VX
	REAL(KIND=LDP) VMAX_S(NCF)
	REAL(KIND=LDP) VMAX_X(NCF)
	FG_JUMP_MAX=0.0_LDP
!
! Store-column peak per frequency.  The jump test is keyed on the
! retained store's own scale, not the new field's: a write is a
! suspicious jump only when it is both >FG_COMMIT_JUMP times the
! stored value AND larger than the whole stored column peak.  The
! old 1E-4*VMAX_F floor rejected every seed-replacement write
! (floor cells hold +1E-30, not 0, so the V2<=0 skip never saw
! them) and left the store frozen at its seed for the whole run.
!
	DO LS2=1,NCF
	  VMAX_S(LS2)=0.0_LDP
	  VMAX_X(LS2)=-1.0_LDP
	END DO
	DO KS=1,KN_ND
	  KC=COMP_DEPTH_IDX(KS)
	  IF(KC .GE. 1)THEN
	    DO LS2=2,NCF
	      V2=STORE(KC,LS,LS2-1)
	      IF(V2 .EQ. V2)VMAX_S(LS2)=MAX(VMAX_S(LS2),V2)
	    END DO
	  END IF
	END DO
	DO KS=1,KN_ND
	  KC=COMP_DEPTH_IDX(KS)
	  IF(KC .GE. 1)THEN
	    JP2=J_PNT(KS,LS)
	    IF(JP2 .GE. 1 .AND. JP2 .LE. NI)THEN
	      DO LS2=2,NCF
	        V1=FLD(JP2,LS2)
	        IF(V1 .NE. V1)THEN
	          FG_JUMP_MAX=HUGE(FG_JUMP_MAX)
	          RETURN
	        END IF
! WRITE_STORE floors non-positive output to 1E-30; floor-bound
! cells are ordinary events, not whole-ray veto grounds.  Judge
! the value actually written.
	        V1=MAX(V1,1.0E-30_LDP)
	        V1=MIN(V1,NONMONO_INEXT_CAP)
! A frequency column whose store peak is at the write floor has no
! committed field of its own to compare against.  Fills there were
! once exempt entirely, which let a converged-classified ray seed a
! column with an extreme value -- the w19 it90 failure channel
! (committed delta ~3E7 while JUMP stayed under the bound).  Fall
! back to the sibling-ray peak at the same field frequency: rays
! sample the same atmosphere, so a genuine seed lands within
! FG_COMMIT_JUMP of what neighbouring rays hold.  Only when every
! ray's column is at the floor (true cold start) is the fill
! exempt.
	        IF(VMAX_S(LS2) .LE. 1.0E-28_LDP)THEN
	          IF(VMAX_X(LS2) .LT. 0.0_LDP)THEN
	            VX=0.0_LDP
	            DO RS=1,KN_NP
	              DO KS2=1,KN_ND
	                KC2=COMP_DEPTH_IDX(KS2)
	                IF(KC2 .GE. 1)THEN
	                  V2=STORE(KC2,RS,LS2-1)
	                  IF(V2 .EQ. V2)VX=MAX(VX,V2)
	                END IF
	              END DO
	            END DO
	            VMAX_X(LS2)=VX
	          END IF
	          IF(VMAX_X(LS2) .GT. 1.0E-28_LDP .AND.
	1          V1 .GT. FG_COMMIT_JUMP*VMAX_X(LS2))THEN
	            FG_JUMP_MAX=HUGE(FG_JUMP_MAX)
	            RETURN
	          END IF
	          CYCLE
	        END IF
	        V2=STORE(KC,LS,LS2-1)
	        IF(V2 .LE. 0.0_LDP .OR. V2 .NE. V2)CYCLE
	        FG_JUMP_MAX=MAX(FG_JUMP_MAX,ABS(V1-V2)/
	1          MAX(V2,VMAX_S(LS2)/FG_COMMIT_JUMP))
	      END DO
	    END IF
	  END IF
	END DO
	RETURN
	END
!
! ----------------------------------------------------------------
! W2 (17-Oct-2026): diverged-ray store heal.  A diverged ray skips
! the store write, so a poisoned slot would re-seed the next
! refresh and re-diverge -- the mechanism that held GS_MAXCHG=1.0
! for the whole run.  Slots that are non-finite, non-positive, or
! above FG_DIV_HEAL_CAP times the per-depth median of the sane
! values are zeroed; the next refresh then takes the seeded floor
! path for those cells (mirrors the MOM heal).
! ----------------------------------------------------------------
	SUBROUTINE HEAL_STORE(STORE)
	REAL(KIND=LDP) STORE(:,:,:)
	INTEGER KS,LS2,KC,NV,NHL,II,JJ
	REAL(KIND=LDP) VMED,VT,VSORT(NCF)
	NHL=0
	DO KS=1,KN_ND
	  KC=COMP_DEPTH_IDX(KS)
	  IF(KC .GE. 1)THEN
	    NV=0
	    DO LS2=1,NCF-1
	      IF(STORE(KC,LS,LS2) .GT. 0.0_LDP .AND.
	1         STORE(KC,LS,LS2) .LT. HUGE(STORE))THEN
	        NV=NV+1
	        VSORT(NV)=STORE(KC,LS,LS2)
	      END IF
	    END DO
	    IF(NV .EQ. 0)CYCLE
!
! Median of the sane values, insertion sort -- per-depth scale
! reference for the poisoned-slot test.  NCF is O(100), so the
! O(N^2) cost only lands on diverged rays.
!
	    DO II=2,NV
	      VT=VSORT(II)
	      JJ=II-1
	      DO WHILE(JJ .GE. 1)
	        IF(VSORT(JJ) .LE. VT)EXIT
	        VSORT(JJ+1)=VSORT(JJ)
	        JJ=JJ-1
	      END DO
	      VSORT(JJ+1)=VT
	    END DO
	    VMED=VSORT((NV+1)/2)
	    DO LS2=1,NCF-1
	      IF(.NOT. (STORE(KC,LS,LS2) .GT. 0.0_LDP) .OR.
	1         .NOT. (STORE(KC,LS,LS2) .LT. HUGE(STORE)) .OR.
	2         STORE(KC,LS,LS2) .GT. FG_DIV_HEAL_CAP*VMED)THEN
	        STORE(KC,LS,LS2)=0.0_LDP
	        NHL=NHL+1
	      END IF
	    END DO
	  END IF
	END DO
!$OMP CRITICAL (FG_HEAL_UPD)
	FG_DIV_HEAL=FG_DIV_HEAL+NHL
!$OMP END CRITICAL (FG_HEAL_UPD)
	RETURN
	END
!
	END
