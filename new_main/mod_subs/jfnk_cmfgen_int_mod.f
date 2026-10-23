!
! JFNK_CMFGEN_INT_MOD: Module for JFNK-CMFGEN interface state.
!
! This module holds the shared state used by the JFNK callback routines
! (residual, JVP, preconditioner, dot product).  It must be compiled
! before jfnk_cmfgen_interface_v1.f.
!
! POPS_WORK, ED_WORK, T_WORK are pointers to the actual POPS, ED, T
! arrays in the caller (SOLVE_FOR_POPS_JFNK_V1 / cmfgen_sub.f).
! The callbacks access them through this module.
!
! editted by Sean Brennan (V4, 2026-09)
!
	MODULE JFNK_CMFGEN_INT_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
	SAVE
!
! Problem dimensions
!
	INTEGER, SAVE :: JFNK_NT=0
	INTEGER, SAVE :: JFNK_ND=0
	INTEGER, SAVE :: JFNK_NION=0
	INTEGER, SAVE :: JFNK_NUM_BNDS=0
	INTEGER, SAVE :: JFNK_DIAG_INDX=0
	INTEGER, SAVE :: JFNK_LUER=6
!
! POPS_BASE: the original populations from the great iteration
! POPS_SAVE: saved POPS during perturbation
! POPS_WORK: the working POPS array (copied from cmfgen_sub.f's local POPS)
! ED_SAVE, T_SAVE: saved ED and T
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: POPS_BASE(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: POPS_SAVE(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: POPS_WORK(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: ED_SAVE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: T_SAVE(:)
!
! DIERECOM/DIECOOL work arrays (local, not in any module)
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_DIERECOM(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_DIECOOL(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_ADDRECOM(:,:)
!
! Z_POP and DO_LEV_DISSOLUTION (local to cmfgen_sub.f, passed to SUP_TO_FULL_V4)
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_Z_POP(:)
	LOGICAL, SAVE :: JFNK_DO_LEV_DISSOL=.FALSE.
!
! Baseline copies used only when JFNK_LAG_FACTORS is enabled.
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_POPION_LAG(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_B_LEV_DIS_LAG(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_A_LEV_DIS_LAG(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_X_LEV_DIS_LAG(:)
	LOGICAL, SAVE :: JFNK_LAG_DO_LEV_DIS=.FALSE.
	LOGICAL, SAVE :: JFNK_LAG_FACTORS_READY=.FALSE.
	LOGICAL, SAVE :: JFNK_LAG_FACTORS_ACTIVE=.FALSE.
!
! STEQ save/restore arrays (added 2026-09-07):
! JFNK_RESIDUAL_V1 modifies SE(ID)%STEQ, STEQ_T, STEQ_ED, and STEQ_T_EHB
! via JFNK_COMPUTE_STEQ_V1.  These arrays are module-level state in
! STEQ_DATA_MOD and may be read by other code paths between JFNK
! evaluations (e.g., the outer great iteration, Ng acceleration, or
! diagnostic output).  Without save/restore, the perturbed STEQ values
! from the last JVP evaluation would be visible to those code paths,
! producing inconsistent results.
!
! JFNK_STEQ_T_SAVE / JFNK_STEQ_ED_SAVE / JFNK_STEQ_T_EHB_SAVE:
!   1D arrays of length ND, saved before and restored after each
!   residual evaluation.
!
! JFNK_SE_STEQ_SAVE:
!   Flat array holding the assembled STEQ vector (NT*ND), saved/restored
!   around each residual evaluation.  The per-ion SE(ID)%STEQ arrays have
!   variable dimensions (NX varies by ion), so we save the assembled flat
!   vector instead.  The assembly is done by JFNK_ASSEMBLE_STEQ_V1 which
!   maps per-ion STEQ into the NT*ND flat vector.
!
! JFNK_SE_STEQ_ION_SAVE / JFNK_SE_STEQ_OFFSET (Finding #4, 2026-09-08):
!   Flat buffer and offset array for saving/restoring the per-ion
!   SE(ID)%STEQ arrays around each residual evaluation.  The per-ion
!   arrays have variable first dimension (SE(ID)%N_SE varies by ion),
!   so they are packed into a single flat buffer with offsets.
!   JFNK_SE_STEQ_OFFSET(ID) is the start index for ion ID in the flat
!   buffer; JFNK_SE_STEQ_OFFSET(NION+1) is the total size.
!   This makes JFNK_RESIDUAL_V1 fully side-effect-free with respect to
!   SE(ID)%STEQ, preventing perturbed STEQ from leaking to other code
!   paths (Ng acceleration, diagnostics) between JFNK evaluations.
!
! JFNK_SE_ADV_SAVE / JFNK_SE_TEHB_SAVE (order-of-operations audit,
! 03-Oct-2026): the save/restore above missed two module arrays that
! JFNK_COMPUTE_STEQ_V1 also rewrites at perturbed populations --
! SE(ID)%STEQ_ADV (advective SE terms from STEQ_ADVEC_V4 /
! STEQ_CO_MOV_DERIV_V3, folded into the assembled residual by
! GENERATE_FULL_MATRIX_V3 and read by SOLVEBA_V13's resacc metric) and
! SE(ID)%T_EHB (electron heat balance, filled by EVALSE_QWVJ during
! the eval sweep).  Both are dimensioned (ND) per ion, so a flat
! NION*ND buffer suffices; ion ID occupies [(ID-1)*ND+1, ID*ND].
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_STEQ_T_SAVE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_STEQ_ED_SAVE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_STEQ_T_EHB_SAVE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_SE_STEQ_ION_SAVE(:)
	INTEGER, ALLOCATABLE, SAVE :: JFNK_SE_STEQ_OFFSET(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_SE_ADV_SAVE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_SE_TEHB_SAVE(:)
!
! Finite-difference step size for JVP (local, not from control_variable_mod)
!
	REAL(KIND=LDP), SAVE :: JFNK_EPS_JVP_LOC=
	1                               SQRT(EPSILON(1.0_LDP))
!
! Flag: whether to recompute STEQ (T) for the perturbed state
!
	LOGICAL, SAVE :: JFNK_FIXED_T_LOC=.TRUE.
	LOGICAL, SAVE :: JFNK_FIXED_NE_LOC=.TRUE.
!
! Per-activation residual mode (03-Oct-2026, two-stage admission):
! JFNK_COUPLED in CONTROL_VARIABLE_MOD is the VADAT capability flag;
! JFNK_COUPLED_LOC is the mode selected for THIS activation by the
! FRZ_GATE/COR_GATE tier logic in SOLVE_FOR_POPS (frozen-J Newton is
! admitted at the looser gate, coupled evals only near the root).
! Every JFNK_COUPLED branch inside residual/JVP/rebuild reads the LOC
! flag so a single executable can mix frozen and coupled activations.
!
	LOGICAL, SAVE :: JFNK_COUPLED_LOC=.FALSE.
!
! Coupled-SOB reach-fail fallback state ([JFNK_COUPLED_SOB], default
! F; 23-Oct-2026).  Under DO_NONMONO_SOB the coupled tier is never
! demand-admitted -- it fires only inside JFNK_SOLVE_CORR_V1 after
! the frozen tier certifies the residual unreachable, and only if the
! eval residual at x=0 first reproduces the production-assembled STEQ
! to 1D-10 relative (otherwise Newton would optimize a different
! equation than the outer iteration).  A failed check latches
! JFNK_CPL_SOB_DEAD: the sweep callback is detached and every later
! reach skip falls back to caller-mode SOLVEBA exactly as before.
! JFNK_CPL_SOB_CHK_DONE records a passed check so it runs once per
! run unless [JFNK_CPL_SOB_EACH]=T re-arms it on every coupled
! admission.  Run scope: JFNK_RETREAT_RESET_V1 clears both on NEWMOD.
!
	LOGICAL, SAVE :: JFNK_CPL_SOB_DEAD=.FALSE.
	LOGICAL, SAVE :: JFNK_CPL_SOB_CHK_DONE=.FALSE.
! Runaway |z| threshold after M^{-1}; default 10 matches MAX_JFNK.
! editted by Sean Brennan (08-Sep-2026)
	REAL(KIND=LDP), SAVE :: JFNK_PRECOND_ZMAX=1.0E1_LDP
!
! Diagonal of the full fixed-J BA at baseline (NT x ND), captured when
! CMAT_RAD is built / when the LU cache is factored.  Used as a
! damped-Jacobi fallback when the LU inverse returns |z| > ZMAX.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BA_DIAG(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_PRECOND_DIAG(:,:)
! Effective diagonal after LU ridge/identity-row repair; PTC rescales
! against the operator actually applied, not the unregularized BA.
	LOGICAL, SAVE :: JFNK_BA_DIAG_READY=.FALSE.
!
! JFNK_PTC_SGN (NT*ND): per-row sign of the raw BA diagonal captured
! at factorization when [JFNK_PTC_SIGNED] is on (default F).  The
! residual PTC term, the factorization fold and the apply rescale all
! multiply the shift by this sign so |d_i + shift| = |d_i| + sigma_i
! on loss-dominated rows instead of cancelling toward zero.  Frozen
! between factorizations so the residual stays a fixed function
! within a Newton solve; all +1 (READY=F) means unsigned historical
! behaviour, which is also the documented fallback for the first
! activation before any factorization has run.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_PTC_SGN(:)
	LOGICAL, SAVE :: JFNK_PTC_SGN_READY=.FALSE.
!
! Cached DIAG BA LU factorization for the JFNK preconditioner
! (Issue 1 — Knoll & Keyes 2004).  Built once by JFNK_FACTOR_PRECOND_V1
! and reused by JFNK_PRECOND_V1 via DGETRS only.  Refresh when GMRES
! lags (Issue 3).
!
! JFNK_LU_FAC(:,:,K)  — equilibrated LU factors from DGETRF at depth K
! JFNK_BA_A(:,:,K)    — original (ridge-augmented) BA for optional DGERFS
! JFNK_IPIV(:,K)      — pivots
! JFNK_ROW_SF/COL_SF  — DGEEQU scale factors
! JFNK_LU_OK(K)       — .TRUE. if depth K factored successfully
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_LU_FAC(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BA_A(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BA_RIDGE_ADD(:,:)
	INTEGER, ALLOCATABLE, SAVE :: JFNK_IPIV(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_ROW_SF(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_COL_SF(:,:)
	LOGICAL, ALLOCATABLE, SAVE :: JFNK_LU_OK(:)
! JFNK_ROW_REG(J,K): rows replaced by the identity row (or repaired as
! null rows/cols in the banded path) at factor time.  At apply time the
! corresponding output components are projected to zero -- a linear
! operation -- instead of solving z_i=v_i and tripping the |z|>1D8
! nonlinear clip.  Null rows carry no BA information, so giving them
! no Newton update is the correct projection.
	LOGICAL, ALLOCATABLE, SAVE :: JFNK_ROW_REG(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_DEPTH_RIDGE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_LU_PIVRAT(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_LU_RCOND(:)
	INTEGER, ALLOCATABLE, SAVE :: JFNK_LU_ITRY(:)
	LOGICAL, SAVE :: JFNK_LU_READY=.FALSE.
!
! JFNK_PRECOND=SCHUR state: two-level (n,T) Schur-complement PC.
! Population sub-block rows 1..JFNK_NT-2 are factored per depth with
! the LSQ path (LU_FAC holds the Cholesky factor of the normal
! equations); JNT/JTN/JTT hold the T column/row/diagonal couplings
! from the raw BA block, and JFNK_SCHUR_S holds the factored dense
! ND x ND Schur complement S = JTT - JTn Jnn^-1 JnT.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_SCHUR_S(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_SCHUR_JNT(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_SCHUR_JTN(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_SCHUR_JTT(:)
	INTEGER, ALLOCATABLE, SAVE :: JFNK_SCHUR_IPIV(:)
	LOGICAL, SAVE :: JFNK_SCHUR_OK=.FALSE.
!
! JFNK_FACTOR_DONE (audit fix F3, 2026-09-11): set .TRUE. once
! JFNK_FACTOR_PRECOND_V1 has run for the current baseline.  The
! factorization-validity decision in JFNK_PRECOND_V1 must use this
! flag, NOT JFNK_LU_READY — the latter is only set on the LU path,
! so with JFNK_PRECOND_TYPE='JACOBI' (or when all depth LUs failed)
! every psolve call re-ran the full (possibly banded) factorization.
!
	LOGICAL, SAVE :: JFNK_FACTOR_DONE=.FALSE.
!
! JFNK_BA_BUILT (17-Oct-2026, W5B): the admission calls
! JFNK_FACTOR_PRECOND_V1(.FALSE.) so the reach/direction gates run
! before the per-depth LU; the real factor call skips the
! GENERATE_FULL_MATRIX_V3 rebuild when this is set.  Cleared by
! JFNK_INIT_V1 and JFNK_REBUILD_BASELINE_V1 (new baseline -> new
! operator).
!
	LOGICAL, SAVE :: JFNK_BA_BUILT=.FALSE.
!
! JFNK_PC_DIAG_ONLY: while set, JFNK_PRECOND_V1 applies the folded
! diagonal (Jacobi) instead of the LU solve and never triggers the
! lazy factorization -- used by the admission precheck so a gate
! skip pays no factorization cost.
!
	LOGICAL, SAVE :: JFNK_PC_DIAG_ONLY=.FALSE.
!
! JFNK_EVAL_LIM: the activation-local eval budget the residual
! evaluator enforces; the driver sets it from JFNK_MAX_EVALS with
! the budget-rung trim.  Initialized to 0 (=off) for safety when a
! residual eval somehow runs outside the driver.
!
	INTEGER, SAVE :: JFNK_EVAL_LIM=0
!
! JFNK_EVAL_BYPASS: lets a single deliberate residual eval through
! the budget/wall-clock sentinels.  The commit-time applied-state
! verification (JFNK_SOLVE_CORR_V1) runs AFTER the solve has spent
! the budget; without a bypass it receives the Y=1e20 sentinel and
! rolls back steps the line search already verified.  The budget
! exists to bound GMRES burn, not to blind the commit check.
!
	LOGICAL, SAVE :: JFNK_EVAL_BYPASS=.FALSE.
	INTEGER, SAVE :: JFNK_N_LU_FACTORS=0
	INTEGER, SAVE :: JFNK_N_PRECOND_REFRESH=0
	INTEGER, SAVE :: JFNK_N_RIDGE_DEPTHS=0
	INTEGER, SAVE :: JFNK_N_LU_FAIL_DEPTH=0
!
! Consecutive admissions that ended in the reach-gate skip
! (JFNK_INFO=-2 via the JFNK_REACH_GATE test).  SOLVE_FOR_POPS reads
! it to stand down the SOB stable-LAMBDA relief once the frozen tier
! has been repeatedly proved unable to move the residual
! ([JFNK_REL_REACH_MAX]); a committed Newton or any non-reach outcome
! resets it.  Written only inside JFNK_SOLVE_CORR_V1.
!
	INTEGER, SAVE :: JFNK_REACH_SKIP_STREAK=0
!
! V4 audit fix (2026-09-07): Per-depth T-fix range for VARFIXT support
! (audit finding #11).  When JFNK_FIXED_T_LOC=.TRUE., T is fixed at all
! depths.  When .FALSE., T is free at all depths UNLESS a partial range
! is specified via JFNK_FIX_T_D_ST / JFNK_FIX_T_D_END (1-indexed, inclusive).
! JFNK_FIX_T_D_ST=0 means no partial range (use JFNK_FIXED_T_LOC only).
!
	INTEGER, SAVE :: JFNK_FIX_T_D_ST=0
	INTEGER, SAVE :: JFNK_FIX_T_D_END=0
!
! FIX_T_AUTO adaptive damped-T state.  When JFNK_TAUTO_ACTIVE is set
! for an activation, the depth range above is NOT a hard mask: T rows
! of F/JVP stay live and every accepted step has its T components
! scaled by JFNK_TAU_ALPHA and capped at JFNK_TAU_DMAX.  The ladder
! (JFNK_TAU_RUNG, JFNK_TAU_NCOMMIT, JFNK_TAU_RE0) persists across
! activations so T is released gradually as the energy residual
! drops; a committed log-T step above JFNK_TAU_EPS reverts a rung.
!
	LOGICAL, SAVE :: JFNK_TAUTO_ACTIVE=.FALSE.
	REAL(KIND=LDP), SAVE :: JFNK_TAU_ALPHA=0.05_LDP
	REAL(KIND=LDP), SAVE :: JFNK_TAU_RE0=0.0_LDP
	REAL(KIND=LDP), SAVE :: JFNK_TAU_RE_LAST=0.0_LDP
	INTEGER, SAVE :: JFNK_TAU_RUNG=1
	INTEGER, SAVE :: JFNK_TAU_NCOMMIT=0
!
! JFNK_BUDGET_RUNG (17-Oct-2026, replaces pause/re-arm): counted
! failures raise it, committed outcomes reset it; the driver trims
! per-activation budgets toward the JFNK_RUNG_MIN_* floors.  Module
! scope so JFNK_RETREAT_RESET_V1 clears it on NEWMOD.
!
	INTEGER, SAVE :: JFNK_BUDGET_RUNG=0
!
! Effective damped domain for this activation: the VARFIXT range when
! set, else all depths.  Kept separate from JFNK_FIX_T_D_* so the
! mask sites can keep their range semantics under both modes.
!
	INTEGER, SAVE :: JFNK_TAU_DOM_ST=0
	INTEGER, SAVE :: JFNK_TAU_DOM_END=0
!
! Preconditioner state
!
	LOGICAL, SAVE :: JFNK_PRECOND_READY=.FALSE.
	CHARACTER(LEN=4), SAVE :: JFNK_METH_SOL='DIAG'
!
! Per-element residual scaling vector.
! Production default is Euclidean (RSCALE=1).  SQRT_STEQ and DIAG_RDG
! remain diagnostic options via VADAT JFNK_RSCALE.
! editted by Sean Brennan (10-Sep-2026)
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_RSCALE(:)
	LOGICAL, SAVE :: JFNK_RSCALE_READY=.FALSE.
	REAL(KIND=LDP), SAVE :: JFNK_RSCALE_FLOOR=1.0E-12_LDP
!
! Radiative STEQ (frozen-J operator splitting):
!
! JFNK_COMPUTE_STEQ_V1 does not call EVALSE_QWVJ_V8, so it only
! returns collisional / autoionization / charge-exchange / advection
! terms.  The outer great iteration has already accumulated the full
! STEQ (radiative + collisional) at the baseline populations n_0.
!
! With J frozen, radiative rates R_ij(J) are constant, but the
! radiative residual is still linear in the populations:
!   STEQ_rad(n) = STEQ_rad(n_0) + BA_rad * x
! where x is the fractional correction (n = n_0 * (1+x)) and BA_rad
! is the radiative part of the CMFGEN BA matrix (already scaled by
! n so that BA*x is a STEQ increment).
!
! JFNK_STEQ_RAD_FLAT = STEQ_outer(n_0) - STEQ_coll(n_0)
! JFNK_CMAT_RAD      = CMAT_full(n_0)  - CMAT_coll(n_0)
!                    (NT x NT diagonal blocks at each depth)
!
! The residual is then
!   F(x) = STEQ_coll(n_0*(1+x)) + STEQ_rad(n_0) + CMAT_RAD * x_frac
! which is the correct frozen-J statistical-equilibrium residual.
! The finite-difference JVP of this F is J_coll + J_rad, matching
! the full-BA preconditioner.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_STEQ_RAD_FLAT(:)
	LOGICAL, SAVE :: JFNK_RAD_READY=.FALSE.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_CMAT_RAD(:,:,:)
	LOGICAL, SAVE :: JFNK_BA_RAD_READY=.FALSE.
!
! Activations since CMAT_RAD was last rebuilt.  The frozen-J
! operator F = STEQ_coll(n) + STEQ_rad(n0) + BA_rad*x shares the
! root of the full system for any BA_rad because the linearization
! term vanishes at x=0, so a stale radiative Jacobian only alters
! the Newton path, never the answer.  Reuse is capped by
! JFNK_CMAT_MAX_AGE because the convergence rate does degrade as
! the baseline drifts; rebuilding costs two full BA assemblies
! (BA_full then BA_coll) per activation.
!
	INTEGER, SAVE :: JFNK_CMAT_AGE=HUGE(1)
	INTEGER, PARAMETER :: JFNK_CMAT_MAX_AGE=4
!
! Banded BA_rad storage (all bands, not just diagonal).
! JFNK_CMAT_RAD_BANDS(:,:,BAND,K) = NT×NT block for band BAND_INDX at
! depth K.  When JFNK_USE_BANDED=F, only the diagonal band
! (BAND_INDX=DIAG_INDX) is populated and the solver behaves identically
! to the original DIAG-only path.  When JFNK_USE_BANDED=T, all
! NUM_BNDS bands are populated (requires SOL_METH=TRI; the off-diagonal
! BA bands are populated by BA_UPDATE_V7 from the VJ cross-depth
! J-variation and by STEQ_ADVEC_V4 for advection).
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_CMAT_RAD_BANDS(:,:,:,:)
	LOGICAL, SAVE :: JFNK_BA_RAD_BANDS_READY=.FALSE.
!
! Banded preconditioner storage (LAPACK banded format).
! JFNK_BANDED_LUAB(2*KL+KU+1, N) stores the LU factors from DGBTRF.
! JFNK_BANDED_IPIV(N) stores the pivot indices.
! KL = KU = (NUM_BNDS/2) * (2*NT - 1) ≈ 2*NT for NUM_BNDS=3.
! When JFNK_USE_BANDED=F, these are not allocated and the solver
! uses the per-depth DIAG LU / Jacobi fallback (legacy behaviour).
!
	INTEGER, SAVE :: JFNK_BANDED_KL=0
	INTEGER, SAVE :: JFNK_BANDED_KU=0
	INTEGER, SAVE :: JFNK_BANDED_N=0
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BANDED_LUAB(:,:)
	INTEGER, ALLOCATABLE, SAVE :: JFNK_BANDED_IPIV(:)
	LOGICAL, SAVE :: JFNK_BANDED_LU_READY=.FALSE.
!
! Banded equilibration scale factors (DGBEQU).
! JFNK_BANDED_ROW_SF(N): row scale factors
! JFNK_BANDED_COL_SF(N): column scale factors
! Applied as: A_eq = diag(ROW_SF) * A * diag(COL_SF)
! Solve: scale b by ROW_SF, DGBTRS, scale x by COL_SF.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BANDED_ROW_SF(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BANDED_COL_SF(:)
	LOGICAL, SAVE :: JFNK_BANDED_EQUIL_OK=.FALSE.
!
! JFNK_BANDED_ARAW: unscaled copy of the assembled banded matrix,
! saved before equilibration/DGBTRF overwrites LUAB.  Used once per
! factorization for a solve-residual check ||A x - b|| / ||b|| on a
! synthetic RHS, so that a mathematically wrong banded LU is caught
! directly rather than inferred from GMRES non-descent behaviour.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BANDED_ARAW(:,:)
	LOGICAL, SAVE :: JFNK_BANDED_CHK_DONE=.FALSE.
!
! JFNK_BTD_*: block-Thomas factorization of the banded matrix
! [JFNK_BANDED_TDMA, default F].  DGBTRF on the N=NT*ND banded system
! costs O(N*KL*KU) ~ 3E15 flops at NT~2E3/ND~70 (infeasible per
! activation); the block recurrence is ND*(4/3)*NT^3 (seconds once per
! activation) with O(ND*NT^2) applies.  JFNK_BTD_L holds the scaled
! lower blocks, JFNK_BTD_LU the LU factors of the modified diagonals
! D~_k = D_k - L_k C_{k-1}, JFNK_BTD_C the transformed uppers
! C_k = D~_k^{-1} U_k.  All (NT,NT,ND); ~5.7 GB at NT=1839/ND=70.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BTD_L(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BTD_LU(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_BTD_C(:,:,:)
	INTEGER, ALLOCATABLE, SAVE :: JFNK_BTD_IPIV(:,:)
	LOGICAL, SAVE :: JFNK_BTD_READY=.FALSE.
!
! JFNK_USE_BANDED_LOC: runtime flag (set in JFNK_INIT_V1 from the
!   JFNK_USE_BANDED VADAT control variable).  When .TRUE., build banded
!   BA_rad (all bands) and use a banded LU preconditioner (DGBTRF/DGBTRS).
!   Requires SOL_METH=TRI (NUM_BNDS=3).  The off-diagonal bands capture
!   cross-depth coupling from VJ (dJ/dopacity) and advection.
!   When .FALSE., the solver uses the original DIAG-only path.
!
	LOGICAL, SAVE :: JFNK_USE_BANDED_LOC=.FALSE.
!
! Persistent work matrices for hot-path BA assembly (Task 1 optimization):
! JFNK_C_MAT / JFNK_D_MAT replace automatic C_MAT(NT,NT) / D_MAT(NT,NT)
! in JFNK_ASSEMBLE_STEQ_V1, JFNK_BUILD_BA_RAD_V1, and
! JFNK_FACTOR_PRECOND_V1.  Allocated once in JFNK_INIT_V1 to eliminate
! stack allocations from the GMRES Krylov iteration hot path.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_C_MAT(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_D_MAT(:,:)
!
! Persistent workspace for hot-path routines (audit §3, 09-Sep-2026):
! These replace large automatic arrays that were on the stack in
! JFNK_RESIDUAL_V1, JFNK_JVP_V1, JFNK_PRECOND_V1, and
! JFNK_COMPUTE_STEQ_V1.  For NT=200, ND=200, the worst-case call chain
! Newton→GMRES→JVP→RESIDUAL→COMPUTE_STEQ had ~7MB of stack-allocated
! automatic arrays, exceeding the default ifort stack size.
!
! JFNK_W_FRAC(NT,ND)    — FRAC in RESIDUAL
! JFNK_W_STEQ_MAT(NT,ND) — STEQ_MAT in RESIDUAL and PRECOND (sequential, not nested)
! JFNK_W_YRAD(N)        — YRAD in RESIDUAL
! JFNK_W_X_PERT(N)      — X_PERT in JVP
! JFNK_W_F_PERT(N)      — F_PERT in JVP
! JFNK_W_V_USE(N)       — V_USE in JVP
! JFNK_W_V_SAVE(NT,ND)  — V_SAVE in PRECOND
! JFNK_W_CNM(200,200)   — CNM in COMPUTE_STEQ
! JFNK_W_DCNM(200,200)  — DCNM in COMPUTE_STEQ
!
!=======================================================================
! REENTRANCY / SHARED-BUFFER INVARIANT (Finding 4, review 2026-09-10)
!=======================================================================
! These JFNK_W_* and JFNK_C_MAT/JFNK_D_MAT buffers are persistent
! module-level scratch space, NOT thread-safe re-entrant workspace.
!
! The correctness of this design depends on the following invariant:
!
!   No subroutine that uses buffer X may be called (directly or
!   transitively) while another live use of buffer X is still on
!   the call stack.
!
! This currently holds by inspection of the call graph:
!   - JFNK_RESIDUAL_V1 uses JFNK_W_FRAC, JFNK_W_STEQ_MAT, JFNK_W_YRAD,
!     JFNK_W_X_SAVE; it calls JFNK_COMPUTE_STEQ_V1 (uses JFNK_W_CNM,
!     JFNK_W_DCNM) and JFNK_BA_RAD_MATVEC_V1 (uses JFNK_C_MAT).
!   - JFNK_JVP_V1 uses JFNK_W_X_PERT, JFNK_W_F_PERT, JFNK_W_V_USE;
!     it calls JFNK_RESIDUAL_V1 (different buffers).
!   - JFNK_PRECOND_V1 uses JFNK_W_V_SAVE, JFNK_W_STEQ_MAT (sequential
!     with RESIDUAL, not nested); it calls JFNK_FACTOR_PRECOND_V1
!     (uses JFNK_C_MAT).
!   - JFNK_ASSEMBLE_STEQ_V1 uses JFNK_D_MAT.
!
! None of these routines call each other while sharing the same buffer.
!
! WARNING: Do NOT call any of these routines from inside an OMP
! parallel region unless the buffer is made OMP-private/threadprivate
! first.  Several of these buffers (JFNK_W_STEQ_MAT, JFNK_W_V_SAVE)
! are written inside !$OMP PARALLEL DO regions within their owning
! subroutine, which is safe only because the outer call is sequential.
!
! If adding a new call path that uses any of these buffers, verify
! that no other live invocation on the call stack uses the same buffer.
!=======================================================================
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_FRAC(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_STEQ_MAT(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_YRAD(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_X_PERT(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_F_PERT(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_V_USE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_X_SAVE(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_V_SAVE(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_CNM(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JFNK_W_DCNM(:,:)
	INTEGER, PARAMETER :: JFNK_N_S_MAX = 200
!
	TYPE JFNK_BA_HOLD_TYPE
	  REAL(KIND=LDP), ALLOCATABLE :: BA(:,:,:,:)
	  REAL(KIND=LDP), ALLOCATABLE :: BA_PAR(:,:,:)
	END TYPE JFNK_BA_HOLD_TYPE
!
! Diagnostic counters (local to this module)
!
! Coupled-residual machinery (V7 port, 02-Oct-2026):
!
! JFNK_SWEEP_CB: procedure pointer to the internal subroutine
! JFNK_RAD_SWEEP_V1 inside CMFGEN_SUB.  It re-runs the continuum
! frequency loop (opacity -> formal solve -> radiative rates, with the
! Knop I_NEXT store frozen) at a perturbed population state, so the
! finite-difference JVP probes dJ/dn -- the fully coupled SE+RTE
! residual used by Arramy et al. (2024).  Under DO_NONMONO_SOB +
! [JFNK_COUPLED_SOB] it additionally replays the SECTION='LINE' pass
! (SUB_SOB_LINE_V3) so Sobolev line rates enter the eval residual the
! same way they enter the production assembly.  Registered once from
! CMFGEN_SUB when USE_JFNK_SOLVER .AND. JFNK_COUPLED; JFNK_RESIDUAL_V1
! checks ASSOCIATED before calling it, and a failed coupled-SOB
! consistency check NULLIFYs it for the rest of the run
! (JFNK_CPL_SOB_DEAD).  The pointer target is an internal procedure --
! valid only while CMFGEN_SUB is active, which always holds (the
! solver runs inside SOLVE_FOR_POPS).
!
! JFNK_SKIP_SE_ZERO: set .TRUE. around the post-sweep
! JFNK_COMPUTE_STEQ_V1 call so the collisional terms are accumulated on
! top of the radiative terms the sweep just wrote into SE%STEQ instead
! of zeroing them away.
!
	ABSTRACT INTERFACE
	  SUBROUTINE JFNK_SWEEP_PROTO(POPS_EVAL)
	    USE SET_KIND_MODULE
	    REAL(KIND=LDP), INTENT(IN) :: POPS_EVAL(:,:)
	  END SUBROUTINE JFNK_SWEEP_PROTO
	END INTERFACE
	PROCEDURE(JFNK_SWEEP_PROTO), POINTER, SAVE :: JFNK_SWEEP_CB => NULL()
	LOGICAL, SAVE :: JFNK_SKIP_SE_ZERO=.FALSE.
!
! JFNK_EDD_START_REC: the EDDFACTOR/ES_J_CONV direct-access record
! index at which the continuum frequency loop starts.  CMFGEN_SUB
! captures it in the normal path right where ACCESS_F is loaded from
! EDD_CONT_REC (or written under COMPUTE_EDDFAC), and the eval sweep
! restores ACCESS_F to it -- COMP_J_BLANK walks the cursor to EOF
! during a sweep, so without the reset the first eval read lands past
! the last record.
!
	INTEGER, SAVE :: JFNK_EDD_START_REC=-1
!
! JFNK_T_START: wall-clock stamp (OMP_GET_WTIME) of the current
! activation, set in JFNK_INIT_V1; JFNK_RESIDUAL_V1 reads it for the
! [JFNK_MAX_SEC] elapsed-time abort (edited by Devin, 15-Oct-2026).
!
	DOUBLE PRECISION, SAVE :: JFNK_T_START=0.0D0
	LOGICAL, SAVE :: JFNK_SEC_HIT=.FALSE.
	INTEGER, SAVE :: JFNK_N_RESIDUAL_EVALS=0
	INTEGER, SAVE :: JFNK_N_RESIDUAL_NAN=0    ! residual evals with NaN/Inf
	INTEGER, SAVE :: JFNK_N_JVP_EVALS=0
	INTEGER, SAVE :: JFNK_N_JVP_NAN=0         ! NaN/Inf entries in last JVP
	INTEGER, SAVE :: JFNK_N_JVP_FAILURES=0    ! JVP calls with excessive NaN
	INTEGER, SAVE :: JFNK_N_NEW_CLIP=0        ! constraints newly crossed in
!                                           last residual eval (audit
!                                           17-Sep-2026, issue 4): a clip
!                                           that was INACTIVE at the
!                                           baseline point jfnk_x but active
!                                           at the perturbed point.  The FD
!                                           JVP retries with smaller h when
!                                           this is nonzero.
	INTEGER, SAVE :: JFNK_N_JVP_HRETRY=0      ! JVP calls needing h-reduction
!
! Diagnostic-only FD step multiplier for the COLPROBE h-scan.  The
! nominal JVP step (~1.5D-8 in ln n) is below the cancellation floor
! of rows whose partial terms reach 1D15-1D19, so a single-h probe
! cannot separate a small true derivative from roundoff.  JVP_V1
! multiplies H_STEP by this factor; the colprobe sweeps it over a
! few decades and resets it to 1 afterward.
!
	REAL(KIND=LDP), SAVE :: JFNK_PROBE_HSCALE=1.0_LDP
!
! JFNK_PROBE_H_LAST: H_STEP actually used by the most recent
! JFNK_JVP_V1 call, after the safety cap, the HSCALE multiplier and
! any clip-crossing retries.  0 after a failed call.  Lets the
! deferral-path probe report the true perturbation ||h*v|| instead
! of recomputing a value the cap may have changed.
!
	REAL(KIND=LDP), SAVE :: JFNK_PROBE_H_LAST=0.0_LDP
!
! JFNK_HADAPT_MULT: post-cap FD-step multiplier selected by the
! [JFNK_HADAPT] plateau finder in JFNK_JVP_V1.  1.0 until the one-shot
! h-sweep runs; the selected value makes every later JVP repeat the
! scanned stable step on its own direction (multipliers in the sweep
! are relative to the post-cap H_STEP, so this multiplies post-cap too).
!
	REAL(KIND=LDP), SAVE :: JFNK_HADAPT_MULT=1.0_LDP
!
! JFNK_BA_SRC_IT: MAIN_COUNTER at the last COMPUTE_BA population
! call, stamped in SOLVE_FOR_POPS where 'BA matrix computed' prints.
! JFNK_BA_A is rebuilt per admission from SE(ID)%BA, so the stored
! operator's age is the outer-BA age, not the admission count.
!
	INTEGER, SAVE :: JFNK_BA_SRC_IT=-1
!
! JFNK_DEFER_PROBE_N: deferral-path probes already run this run
! ([JFNK_DEFER_PROBE] cap, default 0 = off).
!
	INTEGER, SAVE :: JFNK_DEFER_PROBE_N=0
	INTEGER, SAVE :: JFNK_N_PRECOND_APPLIES=0
	INTEGER, SAVE :: JFNK_N_PRECOND_FAILURES=0
!
! Frozen equation-structure flags for JFNK_ASSEMBLE_STEQ_V1.
! GENERATE_FULL_MATRIX_V3 recomputes REPLACE (ground-state row swaps
! to the ionization-conservation equation) and ZERO_STEQ (FIXPOP
! zeroing) from the CURRENT state on every call.  A finite-difference
! perturbation near the FAC threshold can therefore flip a row's
! algebraic structure inside a JVP step, so the FD measures an O(1/h)
! structural artifact rather than a derivative.  The admission and
! baseline-rebuild drivers reset JFNK_EQ_STRUCT_READY; the first
! baseline eval then harvests both flag sets, and all perturbed evals
! replay them through USE_PASSED_REP=T so the JVP differentiates a
! fixed equation set identical to the stored BA's.
!
	LOGICAL, ALLOCATABLE, SAVE :: JFNK_REPLACE_SAV(:,:)
	LOGICAL, ALLOCATABLE, SAVE :: JFNK_ZSTEQ_SAV(:,:)
	LOGICAL, SAVE :: JFNK_EQ_STRUCT_READY=.FALSE.
!
! W11/W12 (17-Oct-2026): census + preconditioner-enrichment state.
!
! JFNK_CENSUS_REG(i): row-regularization mask computed by
! JFNK_RESIDUAL_CENSUS_V1 on the built (unfactored) JFNK_BA_A -- the
! same ROW_TINY near-null and JFNK_POP_FL population-floor tests the
! LU path applies at factor time, reproduced early so the frozen-tier
! census gate and the coupled-tier anchor selection can read them
! before any factorization runs.  Allocated per admission.
!
! JFNK_ENRICH_*: low-rank field-response enrichment of the
! coupled-tier preconditioner ([JFNK_ENRICH], default F).
! JFNK_ENRICH_U(:,m) = M0^-1 W(:,m) where W(:,m) = (J - BA_A) e_c for
! anchor unknown c = JFNK_ENRICH_IDX(m); JFNK_ENRICH_S holds the
! factored K x K capacitance matrix I + E^T U.  The apply in
! JFNK_PRECOND_V1 is then z -> z - U S^-1 E^T z (Sherman-Morrison-
! Woodbury on the left-preconditioned operator).  KLOC is the live
! anchor count; ACTIVE gates the apply.
!
	LOGICAL, ALLOCATABLE, SAVE :: JFNK_CENSUS_REG(:)
	LOGICAL, SAVE :: JFNK_ENRICH_ACTIVE=.FALSE.
	INTEGER, SAVE :: JFNK_ENRICH_KLOC=0
	INTEGER, ALLOCATABLE, SAVE :: JFNK_ENRICH_IDX(:)
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: JFNK_ENRICH_U(:,:)
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: JFNK_ENRICH_S(:,:)
	INTEGER, ALLOCATABLE, SAVE :: JFNK_ENRICH_IPIV(:)
!
! JFNK_PTC_FOLDED: set by JFNK_FACTOR_PRECOND_V1 when the pseudo-time
! diagonal s_i = max(|F_i|, 1e-8*|F|max) has been folded into the
! factored BA blocks.  The FD JVP differentiates G = F + (x-x0)*s/dtau
! so the operator GMRES sees is J + s*I/dtau; bare BA is nearly
! singular on weak rows and the LU emits z ~ 1e2..1e5 along directions
! where the true operator holds the s floor.  With s inside the
! factored matrix the preconditioned direction matches the FD
! operator, and the post-solve D/(D+sig) rescale in JFNK_PRECOND_V1
! must be skipped to avoid applying the shift twice.
!
	LOGICAL, SAVE :: JFNK_PTC_FOLDED=.FALSE.
!
! JFNK_LIN_STAGE ([JFNK_LINTEST], default F): once-per-activation arm
! flag for the JVP operator-consistency probes in JFNK_JVP_V1.  0 =
! armed; the first post-factorization direction with a valid
! production image runs the repeat/homogeneity/additivity evals at the
! same jfnk_x and flips STAGE to 1.  JFNK_INIT_V1 re-arms it each
! activation.  No direction storage is needed: all probes run inside
! one JVP call so a Newton update of jfnk_x between calls cannot
! contaminate the defect.
!
	INTEGER, SAVE :: JFNK_LIN_STAGE=0
!
	END MODULE JFNK_CMFGEN_INT_MOD
