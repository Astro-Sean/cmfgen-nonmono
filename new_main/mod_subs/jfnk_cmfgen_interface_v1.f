!
! JFNK_CMFGEN_INTERFACE_V1: CMFGEN-specific interface routines for the
! JFNK solver (jfnk_solver_v1.f90).
!
! This file provides the callback routines needed by newton_hook_v1:
!   - JFNK_RESIDUAL_V1:  F(x) = STEQ_coll(n0*(1+x)) + STEQ_rad(n0)
!                        + CMAT_RAD * x_frac   (frozen-J, linearized rad)
!   - JFNK_JVP_V1:       J*v = [F(n+eps*v) - F(n)] / eps
!   - JFNK_PRECOND_V1:   z = M^{-1} * v  (fixed-J BA preconditioner)
!   - JFNK_DOT_V1:       d = a . b
!
! POPS, ED, T are accessed via MOD_CMFGEN (ED, T) and via the
! JFNK_CMFGEN_INT_MOD module (POPS_BASE, POPS_SAVE).
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_INIT_V1(NT, ND, NION, NUM_BNDS, DIAG_INDX, LUER,
	1                        POPS, ED_IN, T_IN, Z_POP_IN,
	1                        DO_LEV_DISSOL_IN,
	1                        FIXED_T_MODE, METH_SOL,
	1                        BA_COMPUTED, EPS_JVP)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_ridge_floor, jfnk_ridge_floor_set,
	1                            jfnk_n_reg_proj, jfnk_ptc_scale_ready,
	2                            jfnk_ptc_anchor_ready
	USE CONTROL_VARIABLE_MOD, ONLY: FIXED_NE, JFNK_USE_BANDED,
	1                                JFNK_SELECTIVE_JVAR
	USE STEQ_DATA_MOD, ONLY: SE
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: NT, ND, NION, NUM_BNDS, DIAG_INDX, LUER
	REAL(KIND=LDP), INTENT(IN) :: POPS(NT, ND)
	REAL(KIND=LDP), INTENT(IN) :: ED_IN(ND), T_IN(ND)
	REAL(KIND=LDP), INTENT(IN) :: Z_POP_IN(NT)
	LOGICAL, INTENT(IN) :: DO_LEV_DISSOL_IN
	LOGICAL, INTENT(IN) :: FIXED_T_MODE
	CHARACTER(LEN=*), INTENT(IN) :: METH_SOL
	LOGICAL, INTENT(IN) :: BA_COMPUTED
	REAL(KIND=LDP), INTENT(IN) :: EPS_JVP
	INTEGER :: IOS
	INTEGER :: ID_ION, SE_TOTAL_SIZE
!
	JFNK_NT = NT
	JFNK_ND = ND
	JFNK_NION = NION
	JFNK_NUM_BNDS = NUM_BNDS
	JFNK_DIAG_INDX = DIAG_INDX
	JFNK_LUER = LUER
	JFNK_FIXED_T_LOC = FIXED_T_MODE
	JFNK_FIXED_NE_LOC = FIXED_NE
	JFNK_METH_SOL = METH_SOL
	JFNK_BA_COMPUTED = BA_COMPUTED
	JFNK_EPS_JVP_LOC = EPS_JVP
	JFNK_BA_RAD_READY = .FALSE.
	JFNK_RAD_READY = .FALSE.
	JFNK_BA_RAD_BANDS_READY = .FALSE.
	JFNK_BANDED_LU_READY = .FALSE.
!
! Wire JFNK_USE_BANDED_LOC from the VADAT control variable.  When T,
! the solver builds banded BA_rad and uses a banded LU preconditioner.
! Requires NUM_BNDS > 1 (SOL_METH=TRI); if NUM_BNDS=1 (DIAG), banded
! mode is silently disabled because there are no off-diagonal bands.
!
! The banded LU captures cross-depth J-coupling from the VJ array
! (dJ/dopacity via BA_UPDATE_V7) which populates the off-diagonal
! BA bands in the full Jacobian.  This provides a better
! preconditioner than block-diagonal LU when the radiative coupling
! between adjacent depths is significant.  No selective J-variation
! is required — the off-diagonal bands are populated by the standard
! BA_UPDATE_V7 path whenever NUM_BNDS>1.
!
	JFNK_USE_BANDED_LOC = JFNK_USE_BANDED .AND. (NUM_BNDS .GT. 1)
	IF(JFNK_USE_BANDED .AND. .NOT. JFNK_USE_BANDED_LOC
	1   .AND. LUER .GT. 0)THEN
	  WRITE(LUER,'(A)')'  JFNK: banded LU requested but NUM_BNDS=1 &
	1— using block-diagonal LU'
	END IF
!
! Explicit cache invalidation (ChatGPT audit fix #6, 2026-09-09):
! the BA_RAD/RAD flags were reset above, but the preconditioner and
! RSCALE cache flags could survive between initialisations and let
! JFNK_PRECOND_V1 use a stale BA diagonal / RSCALE from a previous
! model.  Invalidate all *_READY flags up front; they are re-set as
! the corresponding structures are rebuilt.
!
	JFNK_RSCALE_READY = .FALSE.
	JFNK_LU_READY = .FALSE.
	JFNK_FACTOR_DONE = .FALSE.
	JFNK_BA_DIAG_READY = .FALSE.
!
! Disarm the pseudo-transient scale/anchor (numerical audit, 2026-10):
! both are per-activation state frozen from this activation's F(0) by
! JFNK_SET_PTC_SCALE_V1.  Left armed from the previous activation, the
! F(0) test eval below runs before the re-arm, so the stale anchor
! x0 adds (0 - x0)*s/dtau to F(0): frozen mode then fails the
! identity test on every second and later activation, and coupled
! mode adopts the contaminated baseline so every FD JVP inherits the
! constant divided by h.  A stale jfnk_ptc_scale sized for a previous
! N (R-grid revision) would also be indexed past its end.  Clearing
! the readiness flags makes the residual eval skip the PTC term until
! the driver measures the new baseline.
!
	jfnk_ptc_scale_ready = .FALSE.
	jfnk_ptc_anchor_ready = .FALSE.
!
! Reallocate if the problem size changed (R-grid revision).
!
	IF(ALLOCATED(POPS_BASE))THEN
	  IF(SIZE(POPS_BASE,1).NE.NT .OR. SIZE(POPS_BASE,2).NE.ND)THEN
	    DEALLOCATE(POPS_BASE,POPS_SAVE,POPS_WORK)
	    DEALLOCATE(ED_SAVE,T_SAVE)
	    IF(ALLOCATED(JFNK_DIERECOM))DEALLOCATE(JFNK_DIERECOM)
	    IF(ALLOCATED(JFNK_DIECOOL))DEALLOCATE(JFNK_DIECOOL)
	    IF(ALLOCATED(JFNK_ADDRECOM))DEALLOCATE(JFNK_ADDRECOM)
	    IF(ALLOCATED(JFNK_Z_POP))DEALLOCATE(JFNK_Z_POP)
	    IF(ALLOCATED(JFNK_RSCALE))DEALLOCATE(JFNK_RSCALE)
	    IF(ALLOCATED(JFNK_STEQ_RAD_FLAT))DEALLOCATE(JFNK_STEQ_RAD_FLAT)
	    IF(ALLOCATED(JFNK_STEQ_T_SAVE))DEALLOCATE(JFNK_STEQ_T_SAVE)
	    IF(ALLOCATED(JFNK_STEQ_ED_SAVE))DEALLOCATE(JFNK_STEQ_ED_SAVE)
	    IF(ALLOCATED(JFNK_STEQ_T_EHB_SAVE))
	1       DEALLOCATE(JFNK_STEQ_T_EHB_SAVE)
	    IF(ALLOCATED(JFNK_SE_STEQ_ION_SAVE))
	1       DEALLOCATE(JFNK_SE_STEQ_ION_SAVE)
	    IF(ALLOCATED(JFNK_SE_ADV_SAVE))DEALLOCATE(JFNK_SE_ADV_SAVE)
	    IF(ALLOCATED(JFNK_SE_TEHB_SAVE))DEALLOCATE(JFNK_SE_TEHB_SAVE)
	    IF(ALLOCATED(JFNK_SE_STEQ_OFFSET))
	1       DEALLOCATE(JFNK_SE_STEQ_OFFSET)
	    IF(ALLOCATED(JFNK_CMAT_RAD))DEALLOCATE(JFNK_CMAT_RAD)
	    IF(ALLOCATED(JFNK_BA_DIAG))DEALLOCATE(JFNK_BA_DIAG)
	    IF(ALLOCATED(JFNK_PRECOND_DIAG))
	1       DEALLOCATE(JFNK_PRECOND_DIAG)
	    IF(ALLOCATED(JFNK_LU_FAC))DEALLOCATE(JFNK_LU_FAC)
	    IF(ALLOCATED(JFNK_BA_A))DEALLOCATE(JFNK_BA_A)
	    IF(ALLOCATED(JFNK_IPIV))DEALLOCATE(JFNK_IPIV)
	    IF(ALLOCATED(JFNK_ROW_SF))DEALLOCATE(JFNK_ROW_SF)
	    IF(ALLOCATED(JFNK_COL_SF))DEALLOCATE(JFNK_COL_SF)
	    IF(ALLOCATED(JFNK_LU_OK))DEALLOCATE(JFNK_LU_OK)
	    IF(ALLOCATED(JFNK_ROW_REG))DEALLOCATE(JFNK_ROW_REG)
	    IF(ALLOCATED(JFNK_DEPTH_RIDGE))DEALLOCATE(JFNK_DEPTH_RIDGE)
	    IF(ALLOCATED(JFNK_LU_PIVRAT))DEALLOCATE(JFNK_LU_PIVRAT,
	1       JFNK_LU_RCOND, JFNK_LU_ITRY)
	    IF(ALLOCATED(JFNK_C_MAT))DEALLOCATE(JFNK_C_MAT)
	    IF(ALLOCATED(JFNK_D_MAT))DEALLOCATE(JFNK_D_MAT)
	    IF(ALLOCATED(JFNK_W_FRAC))DEALLOCATE(JFNK_W_FRAC)
	    IF(ALLOCATED(JFNK_W_STEQ_MAT))DEALLOCATE(JFNK_W_STEQ_MAT)
	    IF(ALLOCATED(JFNK_W_YRAD))DEALLOCATE(JFNK_W_YRAD)
	    IF(ALLOCATED(JFNK_W_X_PERT))DEALLOCATE(JFNK_W_X_PERT)
	    IF(ALLOCATED(JFNK_W_F_PERT))DEALLOCATE(JFNK_W_F_PERT)
	    IF(ALLOCATED(JFNK_W_V_USE))DEALLOCATE(JFNK_W_V_USE)
	    IF(ALLOCATED(JFNK_W_X_SAVE))DEALLOCATE(JFNK_W_X_SAVE)
	    IF(ALLOCATED(JFNK_W_V_SAVE))DEALLOCATE(JFNK_W_V_SAVE)
	    IF(ALLOCATED(JFNK_W_CNM))DEALLOCATE(JFNK_W_CNM)
	    IF(ALLOCATED(JFNK_W_DCNM))DEALLOCATE(JFNK_W_DCNM)
	    JFNK_LU_READY=.FALSE.
	    JFNK_BA_DIAG_READY=.FALSE.
	  END IF
	END IF
!
! Reallocate NION-shaped arrays if the ion set changed while the
! grid stayed the same (ChatGPT audit fix #5, 2026-09-09).  In a
! model-grid workflow the atomic model (NION) can change between
! epochs even when NT/ND are identical; without this the
! DIERECOM/DIECOOL/ADDRECOM and SE_STEQ_OFFSET arrays keep stale
! dimensions/state from the previous model.  JFNK_NION has already
! been set to the new NION above, so compare the allocated second
! dimension against it.
!
	IF(ALLOCATED(JFNK_DIERECOM))THEN
	  IF(SIZE(JFNK_DIERECOM,2).NE.NION)THEN
	    IF(ALLOCATED(JFNK_DIERECOM))DEALLOCATE(JFNK_DIERECOM)
	    IF(ALLOCATED(JFNK_DIECOOL))DEALLOCATE(JFNK_DIECOOL)
	    IF(ALLOCATED(JFNK_ADDRECOM))DEALLOCATE(JFNK_ADDRECOM)
	    IF(ALLOCATED(JFNK_SE_STEQ_OFFSET))DEALLOCATE(JFNK_SE_STEQ_OFFSET)
	    IF(ALLOCATED(JFNK_SE_STEQ_ION_SAVE))DEALLOCATE(JFNK_SE_STEQ_ION_SAVE)
	    IF(ALLOCATED(JFNK_SE_ADV_SAVE))DEALLOCATE(JFNK_SE_ADV_SAVE)
	    IF(ALLOCATED(JFNK_SE_TEHB_SAVE))DEALLOCATE(JFNK_SE_TEHB_SAVE)
	    JFNK_LU_READY=.FALSE.
	    JFNK_BA_DIAG_READY=.FALSE.
	  END IF
	END IF
!
! Age-limited CMAT_RAD reuse across activations.  Stale BA_rad keeps
! the same fixed point (see jfnk_cmfgen_int_mod), so skipping the
! rebuild saves a collisional-only BA assembly; STEQ_RAD_FLAT is
! still rebuilt every activation because it depends on the current
! baseline.  The preconditioner is unaffected: FACTOR_PRECOND
! re-derives JFNK_BA_DIAG and the ridge from the live BA each call.
!
	IF(JFNK_CMAT_AGE .LT. JFNK_CMAT_MAX_AGE .AND.
	1   ALLOCATED(JFNK_CMAT_RAD))THEN
	  JFNK_BA_RAD_READY=.TRUE.
	  JFNK_CMAT_AGE=JFNK_CMAT_AGE+1
	  IF(JFNK_USE_BANDED_LOC .AND.
	1     ALLOCATED(JFNK_CMAT_RAD_BANDS))
	1   JFNK_BA_RAD_BANDS_READY=.TRUE.
	  WRITE(LUER,'(A,I2,A)')
	1   '  JFNK: reusing CMAT_RAD (age ',JFNK_CMAT_AGE,
	2   ' activations)'
	END IF
!
	IOS=0
	IF(.NOT. ALLOCATED(POPS_BASE)) ALLOCATE(POPS_BASE(NT, ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc POPS_BASE failed'
	IOS=0
	IF(.NOT. ALLOCATED(POPS_SAVE)) ALLOCATE(POPS_SAVE(NT, ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc POPS_SAVE failed'
	IOS=0
	IF(.NOT. ALLOCATED(POPS_WORK)) ALLOCATE(POPS_WORK(NT, ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc POPS_WORK failed'
	IOS=0
	IF(.NOT. ALLOCATED(ED_SAVE)) ALLOCATE(ED_SAVE(ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc ED_SAVE failed'
	IOS=0
	IF(.NOT. ALLOCATED(T_SAVE)) ALLOCATE(T_SAVE(ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc T_SAVE failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_DIERECOM)) ALLOCATE(JFNK_DIERECOM(ND,NION),
	1   STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc DIERECOM failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_DIECOOL)) ALLOCATE(JFNK_DIECOOL(ND,NION),
	1   STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc DIECOOL failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_ADDRECOM)) ALLOCATE(JFNK_ADDRECOM(ND,NION),
	1   STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc ADDRECOM failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_Z_POP)) ALLOCATE(JFNK_Z_POP(NT),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc Z_POP failed'
!
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_RSCALE))
	1   ALLOCATE(JFNK_RSCALE(NT*ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc RSCALE failed'
	JFNK_RSCALE = 1.0_LDP
	JFNK_RSCALE_READY = .FALSE.
!
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_STEQ_RAD_FLAT))
	1   ALLOCATE(JFNK_STEQ_RAD_FLAT(NT*ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc STEQ_RAD_FLAT failed'
	JFNK_STEQ_RAD_FLAT = 0.0_LDP
!
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_STEQ_T_SAVE))
	1   ALLOCATE(JFNK_STEQ_T_SAVE(ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc STEQ_T_SAVE failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_STEQ_ED_SAVE))
	1   ALLOCATE(JFNK_STEQ_ED_SAVE(ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc STEQ_ED_SAVE failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_STEQ_T_EHB_SAVE))
	1   ALLOCATE(JFNK_STEQ_T_EHB_SAVE(ND),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc STEQ_T_EHB_SAVE failed'
!
! V4 audit fix (Finding #4, 2026-09-08): Allocate per-ion STEQ save
! buffer and offset array.  The per-ion SE(ID)%STEQ arrays have variable
! first dimension (SE(ID)%N_SE varies by ion), so we pack them into a
! single flat buffer with an offset array tracking each ion's start.
! JFNK_SE_STEQ_OFFSET(ID) = start index for ion ID in the flat buffer.
! JFNK_SE_STEQ_OFFSET(NION+1) = total size.
!
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_SE_STEQ_OFFSET))
	1   ALLOCATE(JFNK_SE_STEQ_OFFSET(NION+1),STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc SE_STEQ_OFFSET failed'
	JFNK_SE_STEQ_OFFSET = 0
	SE_TOTAL_SIZE = 0
	DO ID_ION = 1, NION
	  JFNK_SE_STEQ_OFFSET(ID_ION) = SE_TOTAL_SIZE
	  IF(ASSOCIATED(SE(ID_ION)%STEQ))THEN
	    SE_TOTAL_SIZE = SE_TOTAL_SIZE + SIZE(SE(ID_ION)%STEQ)
	  END IF
	END DO
	JFNK_SE_STEQ_OFFSET(NION+1) = SE_TOTAL_SIZE
!
	IOS=0
	IF(ALLOCATED(JFNK_SE_STEQ_ION_SAVE))THEN
	  IF(SIZE(JFNK_SE_STEQ_ION_SAVE) .NE. SE_TOTAL_SIZE)
	1     DEALLOCATE(JFNK_SE_STEQ_ION_SAVE)
	END IF
	IF(.NOT. ALLOCATED(JFNK_SE_STEQ_ION_SAVE))THEN
	  IF(SE_TOTAL_SIZE .GT. 0)THEN
	    ALLOCATE(JFNK_SE_STEQ_ION_SAVE(SE_TOTAL_SIZE),STAT=IOS)
	  ELSE
	    ALLOCATE(JFNK_SE_STEQ_ION_SAVE(1),STAT=IOS)
	  END IF
	END IF
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc SE_STEQ_ION_SAVE failed'
	JFNK_SE_STEQ_ION_SAVE = 0.0_LDP
!
! Per-ion advective/EHB STEQ save buffers (order-of-operations audit,
! 03-Oct-2026): SE%STEQ_ADV and SE%T_EHB are each dimensioned (ND)
! per ion, so NION*ND flat buffers suffice.  Same size-check realloc
! pattern as JFNK_SE_STEQ_ION_SAVE above.
!
	IF(ALLOCATED(JFNK_SE_ADV_SAVE))THEN
	  IF(SIZE(JFNK_SE_ADV_SAVE) .NE. NION*ND)
	1     DEALLOCATE(JFNK_SE_ADV_SAVE)
	END IF
	IF(.NOT. ALLOCATED(JFNK_SE_ADV_SAVE))THEN
	  ALLOCATE(JFNK_SE_ADV_SAVE(MAX(1,NION*ND)),STAT=IOS)
	END IF
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc SE_ADV_SAVE failed'
	JFNK_SE_ADV_SAVE = 0.0_LDP
!
	IF(ALLOCATED(JFNK_SE_TEHB_SAVE))THEN
	  IF(SIZE(JFNK_SE_TEHB_SAVE) .NE. NION*ND)
	1     DEALLOCATE(JFNK_SE_TEHB_SAVE)
	END IF
	IF(.NOT. ALLOCATED(JFNK_SE_TEHB_SAVE))THEN
	  ALLOCATE(JFNK_SE_TEHB_SAVE(MAX(1,NION*ND)),STAT=IOS)
	END IF
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc SE_TEHB_SAVE failed'
	JFNK_SE_TEHB_SAVE = 0.0_LDP
!
! Task 1: Allocate persistent work matrices for hot-path BA assembly.
! These replace automatic C_MAT(NT,NT) / D_MAT(NT,NT) in
! JFNK_ASSEMBLE_STEQ_V1, JFNK_BUILD_BA_RAD_V1, and
! JFNK_FACTOR_PRECOND_V1, eliminating stack allocations from the
! GMRES Krylov iteration and residual evaluation hot paths.
!
	IOS=0
	IF(ALLOCATED(JFNK_C_MAT))THEN
	  IF(SIZE(JFNK_C_MAT,1).NE.NT .OR. SIZE(JFNK_C_MAT,2).NE.NT)
	1   DEALLOCATE(JFNK_C_MAT)
	END IF
	IF(.NOT. ALLOCATED(JFNK_C_MAT))THEN
	  ALLOCATE(JFNK_C_MAT(NT, NT), STAT=IOS)
	  IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc C_MAT failed'
	END IF
	JFNK_C_MAT = 0.0_LDP
!
	IOS=0
	IF(ALLOCATED(JFNK_D_MAT))THEN
	  IF(SIZE(JFNK_D_MAT,1).NE.NT .OR. SIZE(JFNK_D_MAT,2).NE.NT)
	1   DEALLOCATE(JFNK_D_MAT)
	END IF
	IF(.NOT. ALLOCATED(JFNK_D_MAT))THEN
	  ALLOCATE(JFNK_D_MAT(NT, NT), STAT=IOS)
	  IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc D_MAT failed'
	END IF
	JFNK_D_MAT = 0.0_LDP
!
! Allocate persistent workspace arrays (audit §3, 09-Sep-2026):
! These replace large automatic arrays on the stack.
!
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_FRAC))
	1   ALLOCATE(JFNK_W_FRAC(NT, ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_FRAC failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_STEQ_MAT))
	1   ALLOCATE(JFNK_W_STEQ_MAT(NT, ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_STEQ_MAT failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_YRAD))
	1   ALLOCATE(JFNK_W_YRAD(NT*ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_YRAD failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_X_PERT))
	1   ALLOCATE(JFNK_W_X_PERT(NT*ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_X_PERT failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_F_PERT))
	1   ALLOCATE(JFNK_W_F_PERT(NT*ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_F_PERT failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_V_USE))
	1   ALLOCATE(JFNK_W_V_USE(NT*ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_V_USE failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_X_SAVE))
	1   ALLOCATE(JFNK_W_X_SAVE(NT*ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_X_SAVE failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_V_SAVE))
	1   ALLOCATE(JFNK_W_V_SAVE(NT, ND), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_V_SAVE failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_CNM))
	1   ALLOCATE(JFNK_W_CNM(JFNK_N_S_MAX, JFNK_N_S_MAX),
	1   STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_CNM failed'
	IOS=0
	IF(.NOT. ALLOCATED(JFNK_W_DCNM))
	1   ALLOCATE(JFNK_W_DCNM(JFNK_N_S_MAX, JFNK_N_S_MAX),
	1   STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_INIT_V1: alloc W_DCNM failed'
!
	POPS_BASE = POPS
	POPS_WORK = POPS
	ED_SAVE = ED_IN
	T_SAVE = T_IN
	JFNK_Z_POP = Z_POP_IN
	JFNK_DO_LEV_DISSOL = DO_LEV_DISSOL_IN
	JFNK_DIERECOM = 0.0_LDP
	JFNK_DIECOOL = 0.0_LDP
	JFNK_ADDRECOM = 0.0_LDP
!
	JFNK_N_RESIDUAL_EVALS = 0
	JFNK_N_RESIDUAL_NAN = 0
	JFNK_N_JVP_EVALS = 0
	JFNK_N_JVP_NAN = 0
	JFNK_N_JVP_FAILURES = 0
	JFNK_N_PRECOND_APPLIES = 0
	JFNK_N_PRECOND_FAILURES = 0
	JFNK_N_LU_FACTORS = 0
	JFNK_N_PRECOND_REFRESH = 0
	jfnk_n_reg_proj = 0
	JFNK_LU_READY = .FALSE.
!
! Fix #6 (review 2026-09-11): Reset the ridge-floor SAVE state so
! that the first model does not determine the floor for subsequent
! models.  Without this, jfnk_ridge_floor_set=.TRUE. from a previous
! model causes JFNK_FACTOR_PRECOND_V1 to skip the floor capture,
! inheriting a stale ridge floor.
!
	jfnk_ridge_floor_set = .FALSE.
	jfnk_ridge_floor = 0.0_LDP
!
	END SUBROUTINE JFNK_INIT_V1


!
! JFNK_COMPUTE_RSCALE_V1: Compute per-element residual scaling from
! the baseline STEQ vector.
!
! The scaling is:
!   RSCALE(i) = 1 / max(|STEQ_baseline(i)|, FLOOR_USE)
!   FLOOR_USE = max(JFNK_RSCALE_FLOOR, 1e-8*|F|_max, 1e-12)
!
! This normalizes each residual element so that species with small STEQ
! (trace elements like Fe) contribute equally to the norm as species
! with large STEQ (dominant elements like H).  Without this scaling,
! the JFNK solver may converge H while ignoring Fe.
!
! The scaling is computed once from the initial STEQ and reused for
! the entire JFNK solve.  It is applied in JFNK_DOT_V1 as a weighted
! dot product:  dot = sum (RSCALE(i) * A(i)) * (RSCALE(i) * B(i))
!
! added 2026-09-07
!
	SUBROUTINE JFNK_COMPUTE_RSCALE_V1(N, F_BASE)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE MOD_CMFGEN, ONLY: ED
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_RSCALE_MODE, JFNK_BA_RIDGE_FAC
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(IN) :: F_BASE(N)
!
	INTEGER I, K, J
	REAL(KIND=LDP) ABS_FI, FMAX, FLOOR_USE, DABS, DMAX_K, RIDGE_K, RIDGE_FAC
	INTEGER LUER, ERROR_LU
	EXTERNAL ERROR_LU
!
! Production default is Euclidean (RSCALE=1).  The 1/sqrt(|STEQ|)
! weighting (SQRT_STEQ) and the diagnostic sqrt(|D|/(|D|+ridge))
! option (DIAG_RDG) remain available via VADAT JFNK_RSCALE.
! Full-system audit (2026-09-10): Euclidean reduced the unweighted
! STEQ residual ~54%; 1/sqrt(|STEQ|) reduced it ~0%.
! editted by Sean Brennan (10-Sep-2026)
!
	LUER = ERROR_LU()
	JFNK_RSCALE(1:N) = 1.0_LDP
!
	IF(JFNK_RSCALE_MODE(1:8) .EQ. 'SQRT_STE')THEN
	  FMAX = 0.0_LDP
	  DO I = 1, N
	    FMAX = MAX(FMAX, ABS(REAL(F_BASE(I), KIND=LDP)))
	  END DO
	  FLOOR_USE = MAX(JFNK_RSCALE_FLOOR, 1.0E-8_LDP * FMAX, 1.0E-12_LDP)
	  DO I = 1, N
	    ABS_FI = ABS(REAL(F_BASE(I), KIND=LDP))
	    IF(ABS_FI .LT. FLOOR_USE) ABS_FI = FLOOR_USE
	    JFNK_RSCALE(I) = 1.0_LDP / SQRT(ABS_FI)
	  END DO
	ELSE IF(JFNK_RSCALE_MODE(1:4) .EQ. 'DIAG')THEN
	  IF(JFNK_BA_DIAG_READY)THEN
	    DO K = 1, JFNK_ND
	      DMAX_K = 0.0_LDP
	      DO J = 1, JFNK_NT
	        DMAX_K = MAX(DMAX_K, ABS(JFNK_BA_DIAG(J, K)))
	      END DO
	      DO J = 1, JFNK_NT
	        I = (K-1)*JFNK_NT + J
	        IF(I .GT. N) CYCLE
	        IF(J .EQ. JFNK_NT-1)THEN
	          JFNK_RSCALE(I) = 1.0_LDP/MAX(ED(K), 1.0E-30_LDP)
	        ELSE
	          JFNK_RSCALE(I) = 1.0_LDP/
	1             MAX(ABS(JFNK_BA_DIAG(J, K)), 1.0E-30_LDP*DMAX_K)
	        END IF
	      END DO
	    END DO
	  ELSE
	    IF(LUER .GT. 0)WRITE(LUER,'(A)')
	1     '  JFNK: RSCALE=DIAG requested but BA diagonal not yet'
	2     //' available -- keeping EUCLIDEAN for this call'
	  END IF
	ELSE IF(JFNK_RSCALE_MODE(1:8) .EQ. 'DIAG_RDG')THEN
	  RIDGE_FAC = JFNK_BA_RIDGE_FAC
	  IF(RIDGE_FAC .LE. 0.0_LDP) RIDGE_FAC = 1.0E-3_LDP
	  IF(JFNK_BA_DIAG_READY)THEN
	    DO K = 1, JFNK_ND
	      DMAX_K = 0.0_LDP
	      DO J = 1, JFNK_NT
	        DMAX_K = MAX(DMAX_K, ABS(JFNK_BA_DIAG(J, K)))
	      END DO
	      RIDGE_K = RIDGE_FAC * MAX(DMAX_K, 1.0_LDP)
	      DO J = 1, JFNK_NT
	        I = (K-1)*JFNK_NT + J
	        IF(I .GT. N) CYCLE
	        DABS = ABS(JFNK_BA_DIAG(J, K))
	        JFNK_RSCALE(I) = SQRT(DABS / (DABS + RIDGE_K + 1.0E-30_LDP))
	      END DO
	    END DO
	  END IF
	END IF
!
	JFNK_RSCALE_READY = .TRUE.
	IF(LUER .GT. 0)THEN
	  WRITE(LUER,'(A,A)')'  JFNK: RSCALE mode=', TRIM(JFNK_RSCALE_MODE)
	END IF
!
	END SUBROUTINE JFNK_COMPUTE_RSCALE_V1


!
! JFNK_ASSEMBLE_STEQ_V1: Assemble SE(ID)%STEQ into a flat vector.
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_ASSEMBLE_STEQ_V1(STEQ_OUT, NT, ND, NION, NUM_BNDS,
	1                                 DIAG_INDX, POPS)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: LTE_MODEL, USE_ELEC_HEAT_BAL
	USE JFNK_CMFGEN_INT_MOD, ONLY: JFNK_C_MAT, JFNK_D_MAT,
	1     JFNK_REPLACE_SAV, JFNK_ZSTEQ_SAV, JFNK_EQ_STRUCT_READY
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: NT, ND, NION, NUM_BNDS, DIAG_INDX
	REAL(KIND=LDP), INTENT(IN) :: POPS(NT, ND)
	REAL(KIND=LDP), INTENT(OUT) :: STEQ_OUT(NT, ND)
!
	LOGICAL ZERO_STEQ(NT)
	LOGICAL REPLACE(NION)
	LOGICAL FIRST_MATRIX, LAST_MATRIX, USE_PASSED_REP
	INTEGER K, DEPTH_INDX, BAND_INDX
!
	STEQ_OUT = 0.0_LDP
	JFNK_C_MAT = 0.0_LDP
	JFNK_D_MAT = 0.0_LDP
	ZERO_STEQ = .FALSE.
	REPLACE = .FALSE.
	FIRST_MATRIX = .TRUE.
	LAST_MATRIX = .FALSE.
!
! When JFNK_EQ_STRUCT_READY=F (baseline eval after an admission or
! rebuild reset) harvest the REPLACE/ZERO_STEQ decisions GENERATE
! makes per depth; perturbed evals replay the frozen set so the FD
! JVP differentiates a fixed equation structure.  Freezing at the
! baseline is also the equation structure the stored BA was built
! with, keeping M and the FD operator structurally consistent.
!
	IF(.NOT. ALLOCATED(JFNK_REPLACE_SAV))THEN
	  ALLOCATE(JFNK_REPLACE_SAV(NION,ND),JFNK_ZSTEQ_SAV(NT,ND))
	ELSE IF(SIZE(JFNK_REPLACE_SAV,1).NE.NION .OR.
	1        SIZE(JFNK_REPLACE_SAV,2).NE.ND .OR.
	2        SIZE(JFNK_ZSTEQ_SAV,1).NE.NT .OR.
	3        SIZE(JFNK_ZSTEQ_SAV,2).NE.ND)THEN
	  DEALLOCATE(JFNK_REPLACE_SAV,JFNK_ZSTEQ_SAV)
	  ALLOCATE(JFNK_REPLACE_SAV(NION,ND),JFNK_ZSTEQ_SAV(NT,ND))
	END IF
	USE_PASSED_REP = JFNK_EQ_STRUCT_READY
!
	DEPTH_INDX = 1
	BAND_INDX = DIAG_INDX
!
	DO K = 1, ND
	  DEPTH_INDX = K
	  IF(K .EQ. ND) LAST_MATRIX = .TRUE.
	  IF(USE_PASSED_REP)THEN
	    REPLACE(1:NION) = JFNK_REPLACE_SAV(1:NION,K)
	    ZERO_STEQ(1:NT) = JFNK_ZSTEQ_SAV(1:NT,K)
	  END IF
	  CALL GENERATE_FULL_MATRIX_V3(
	1       JFNK_D_MAT, STEQ_OUT(1,K), POPS, REPLACE, ZERO_STEQ,
	1       NT, ND, NION, NUM_BNDS, BAND_INDX, DIAG_INDX, DEPTH_INDX,
	1       FIRST_MATRIX, LAST_MATRIX, USE_PASSED_REP)
	  IF(.NOT. USE_PASSED_REP)THEN
	    JFNK_REPLACE_SAV(1:NION,K) = REPLACE(1:NION)
	    JFNK_ZSTEQ_SAV(1:NT,K) = ZERO_STEQ(1:NT)
	  END IF
	  FIRST_MATRIX = .FALSE.
	END DO
	IF(.NOT. USE_PASSED_REP)JFNK_EQ_STRUCT_READY = .TRUE.
!
! GENERATE_FULL_MATRIX_V3 already writes the charge and energy equations:
!   STEQ_VEC(NT-1) = STEQ_ED(K)
!   STEQ_VEC(NT)   = STEQ_T(K)  (or STEQ_T_EHB when USE_ELEC_HEAT_BAL)
! Do not add STEQ_ED / STEQ_T again: that doubled F_ED(0) and F_T(0)
! relative to the outer-iteration STEQ (and mixed STEQ_T into the EHB
! equation).  Species equations are assembled only from SE(ID)%STEQ.
!
	END SUBROUTINE JFNK_ASSEMBLE_STEQ_V1


!
! JFNK_SYNC_ATM_FROM_POPS: Unpack POPS into ATM(ID)%XzV / DXzV and
! recompute LTE / full-atom populations.  CMFGEN stores super-level
! populations both in the packed POPS array and in ATM(ID)%XzV;
! STEQ_MULTI_V10 reads ATM, so POPS_WORK must be copied before any
! collisional residual evaluation.
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_SYNC_ATM_FROM_POPS(POPS_IN)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE MOD_CMFGEN, ONLY: ATM, NUM_IONS, ED, T, POPION
	USE MOD_LEV_DIS_BLK, ONLY: B_LEV_DIS, A_LEV_DIS, X_LEV_DIS,
	1   MOD_DO_LEV_DIS, MOD_HOLD_LEV_DIS
	IMPLICIT NONE
	REAL(KIND=LDP), INTENT(IN) :: POPS_IN(JFNK_NT, JFNK_ND)
	INTEGER ID, I, J
!
	DO ID = 1, NUM_IONS - 1
	  CALL POPTOION(POPS_IN, ATM(ID)%XzV, ATM(ID)%DXzV, ED, T,
	1     ATM(ID)%EQXzV, ATM(ID)%NXzV, JFNK_NT, JFNK_ND,
	1     ATM(ID)%XzV_PRES)
	END DO
!
! Occupation-probability constants depend on ED, T and the ion density.
! Recompute them before expanding super-levels, matching cmfgen_sub.f.
!
	IF(JFNK_LAG_FACTORS_ACTIVE .AND. JFNK_LAG_FACTORS_READY)THEN
	  POPION = JFNK_POPION_LAG
	  B_LEV_DIS = JFNK_B_LEV_DIS_LAG
	  A_LEV_DIS = JFNK_A_LEV_DIS_LAG
	  X_LEV_DIS = JFNK_X_LEV_DIS_LAG
	  MOD_DO_LEV_DIS = JFNK_LAG_DO_LEV_DIS
	  MOD_HOLD_LEV_DIS = .TRUE.
	ELSE
	  MOD_HOLD_LEV_DIS = .FALSE.
	  DO J = 1, JFNK_ND
	    POPION(J) = 0.0_LDP
	    DO I = 1, JFNK_NT
	      IF(JFNK_Z_POP(I) .GT. 0.01_LDP)
	1         POPION(J) = POPION(J) + POPS_IN(I, J)
	    END DO
	  END DO
	  CALL COMP_LEV_DIS_BLK(ED, POPION, T, JFNK_DO_LEV_DISSOL,
	1     JFNK_ND)
	END IF
	CALL SUP_TO_FULL_V4(POPS_IN, JFNK_Z_POP, JFNK_DO_LEV_DISSOL,
	1   JFNK_ND, JFNK_NT)
	MOD_HOLD_LEV_DIS = .FALSE.
!
	END SUBROUTINE JFNK_SYNC_ATM_FROM_POPS
!
!
!
! Copy the occupation-probability inputs after the baseline sync.  The
! BA includes direct ED and ion-density columns, but not derivatives of
! these aggregate dissolution factors.
	SUBROUTINE JFNK_CAPTURE_LAG_FACTORS_V1(OK)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE MOD_CMFGEN, ONLY: POPION
	USE MOD_LEV_DIS_BLK, ONLY: B_LEV_DIS, A_LEV_DIS, X_LEV_DIS,
	1   MOD_DO_LEV_DIS
	IMPLICIT NONE
	LOGICAL, INTENT(OUT) :: OK
	INTEGER IOS
!
	OK=.FALSE.
	JFNK_LAG_FACTORS_READY=.FALSE.
	IF(.NOT. ALLOCATED(B_LEV_DIS) .OR.
	1   .NOT. ALLOCATED(A_LEV_DIS) .OR.
	2   .NOT. ALLOCATED(X_LEV_DIS) .OR.
	3   .NOT. ALLOCATED(POPION))RETURN
	IF(SIZE(B_LEV_DIS) .NE. JFNK_ND .OR.
	1   SIZE(A_LEV_DIS) .NE. JFNK_ND .OR.
	2   SIZE(X_LEV_DIS) .NE. JFNK_ND .OR.
	3   SIZE(POPION) .NE. JFNK_ND)RETURN
	IF(ALLOCATED(JFNK_POPION_LAG))DEALLOCATE(JFNK_POPION_LAG)
	IF(ALLOCATED(JFNK_B_LEV_DIS_LAG))DEALLOCATE(JFNK_B_LEV_DIS_LAG)
	IF(ALLOCATED(JFNK_A_LEV_DIS_LAG))DEALLOCATE(JFNK_A_LEV_DIS_LAG)
	IF(ALLOCATED(JFNK_X_LEV_DIS_LAG))DEALLOCATE(JFNK_X_LEV_DIS_LAG)
	ALLOCATE(JFNK_POPION_LAG(JFNK_ND),
	1        JFNK_B_LEV_DIS_LAG(JFNK_ND),
	2        JFNK_A_LEV_DIS_LAG(JFNK_ND),
	3        JFNK_X_LEV_DIS_LAG(JFNK_ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  IF(ALLOCATED(JFNK_POPION_LAG))DEALLOCATE(JFNK_POPION_LAG)
	  IF(ALLOCATED(JFNK_B_LEV_DIS_LAG))DEALLOCATE(JFNK_B_LEV_DIS_LAG)
	  IF(ALLOCATED(JFNK_A_LEV_DIS_LAG))DEALLOCATE(JFNK_A_LEV_DIS_LAG)
	  IF(ALLOCATED(JFNK_X_LEV_DIS_LAG))DEALLOCATE(JFNK_X_LEV_DIS_LAG)
	  RETURN
	END IF
	JFNK_POPION_LAG=POPION
	JFNK_B_LEV_DIS_LAG=B_LEV_DIS
	JFNK_A_LEV_DIS_LAG=A_LEV_DIS
	JFNK_X_LEV_DIS_LAG=X_LEV_DIS
	JFNK_LAG_DO_LEV_DIS=MOD_DO_LEV_DIS
	JFNK_LAG_FACTORS_READY=.TRUE.
	OK=.TRUE.
	END SUBROUTINE JFNK_CAPTURE_LAG_FACTORS_V1
!
!
!
! JFNK_SAVE_MOD_STEQ_V1 / JFNK_RESTORE_MOD_STEQ_V1: bracket for any
! code path that calls JFNK_COMPUTE_STEQ_V1 or a JFNK eval sweep, both
! of which rewrite the module-level STEQ state at trial populations.
! Covers the full set: STEQ_T, STEQ_ED, STEQ_T_EHB, per-ion SE%STEQ,
! SE%STEQ_ADV and SE%T_EHB.  Order-of-operations audit (03-Oct-2026):
! the advective and EHB per-ion arrays were not in the original
! residual save/restore, so they leaked perturbed values to
! GENERATE_FULL_MATRIX_V3 / SOLVEBA_V13's resacc metric between evals.
!
	SUBROUTINE JFNK_SAVE_MOD_STEQ_V1()
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE STEQ_DATA_MOD
	IMPLICIT NONE
	INTEGER ID_ION, OFF, OFF2, NSZ
!
	JFNK_STEQ_T_SAVE = STEQ_T
	JFNK_STEQ_ED_SAVE = STEQ_ED
	JFNK_STEQ_T_EHB_SAVE = STEQ_T_EHB
	DO ID_ION = 1, JFNK_NION
	  OFF = JFNK_SE_STEQ_OFFSET(ID_ION)
	  IF(ASSOCIATED(SE(ID_ION)%STEQ))THEN
	    JFNK_SE_STEQ_ION_SAVE(OFF+1:OFF+SIZE(SE(ID_ION)%STEQ)) =
	1       RESHAPE(SE(ID_ION)%STEQ, (/SIZE(SE(ID_ION)%STEQ)/))
	  END IF
	  OFF2 = (ID_ION-1)*JFNK_ND
	  IF(ASSOCIATED(SE(ID_ION)%STEQ_ADV))THEN
	    NSZ = MIN(SIZE(SE(ID_ION)%STEQ_ADV), JFNK_ND)
	    JFNK_SE_ADV_SAVE(OFF2+1:OFF2+NSZ) =
	1       SE(ID_ION)%STEQ_ADV(1:NSZ)
	  END IF
	  IF(ASSOCIATED(SE(ID_ION)%T_EHB))THEN
	    NSZ = MIN(SIZE(SE(ID_ION)%T_EHB), JFNK_ND)
	    JFNK_SE_TEHB_SAVE(OFF2+1:OFF2+NSZ) =
	1       SE(ID_ION)%T_EHB(1:NSZ)
	  END IF
	END DO
!
	END SUBROUTINE JFNK_SAVE_MOD_STEQ_V1
!
!
!
	SUBROUTINE JFNK_RESTORE_MOD_STEQ_V1()
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE STEQ_DATA_MOD
	IMPLICIT NONE
	INTEGER ID_ION, OFF, OFF2, NSZ
!
	STEQ_T = JFNK_STEQ_T_SAVE
	STEQ_ED = JFNK_STEQ_ED_SAVE
	STEQ_T_EHB = JFNK_STEQ_T_EHB_SAVE
	DO ID_ION = 1, JFNK_NION
	  OFF = JFNK_SE_STEQ_OFFSET(ID_ION)
	  IF(ASSOCIATED(SE(ID_ION)%STEQ))THEN
	    SE(ID_ION)%STEQ = RESHAPE(
	1       JFNK_SE_STEQ_ION_SAVE(OFF+1:OFF+SIZE(SE(ID_ION)%STEQ)),
	1       SHAPE(SE(ID_ION)%STEQ))
	  END IF
	  OFF2 = (ID_ION-1)*JFNK_ND
	  IF(ASSOCIATED(SE(ID_ION)%STEQ_ADV))THEN
	    NSZ = MIN(SIZE(SE(ID_ION)%STEQ_ADV), JFNK_ND)
	    SE(ID_ION)%STEQ_ADV(1:NSZ) =
	1       JFNK_SE_ADV_SAVE(OFF2+1:OFF2+NSZ)
	  END IF
	  IF(ASSOCIATED(SE(ID_ION)%T_EHB))THEN
	    NSZ = MIN(SIZE(SE(ID_ION)%T_EHB), JFNK_ND)
	    SE(ID_ION)%T_EHB(1:NSZ) =
	1       JFNK_SE_TEHB_SAVE(OFF2+1:OFF2+NSZ)
	  END IF
	END DO
!
	END SUBROUTINE JFNK_RESTORE_MOD_STEQ_V1


!
! JFNK_BUILD_BA_RAD_V1: Build the frozen-J radiative Jacobian blocks.
!
! CMAT_RAD(i,j,k) = CMAT_full(i,j,k) - CMAT_coll(i,j,k)
! where CMAT_* are the GENERATE_FULL_MATRIX_V3 diagonal BA blocks
! (already scaled by POPS so that CMAT*x is a STEQ increment for
! the fractional correction x).
!
! The full BA from the outer iteration is saved, a collisional-only
! BA is computed, both are expanded to NT x NT, and the difference
! is stored.  The outer-iteration BA is then restored so the GMRES
! preconditioner still sees the full Jacobian.
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_BUILD_BA_RAD_V1()
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
	TYPE(JFNK_BA_HOLD_TYPE), ALLOCATABLE :: BA_HOLD(:)
	REAL(KIND=LDP), ALLOCATABLE :: BA_ED_HOLD(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: BA_T_HOLD(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: BA_T_EHB_HOLD(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: BA_T_PAR_HOLD(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: BA_T_PAR_EHB_HOLD(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: BA_ADV_HOLD(:,:)
	REAL(KIND=LDP) DUMMY_STEQ(JFNK_NT)
	LOGICAL ZERO_STEQ(JFNK_NT)
	LOGICAL REPLACE(JFNK_NION)
	LOGICAL FIRST_MATRIX, LAST_MATRIX, USE_PASSED_REP
	INTEGER ID, J, K, N1, N2, N3, N4, IOS, IB
	INTEGER LUER, ERROR_LU
	EXTERNAL ERROR_LU
!
! editted by Devin, 08-Oct-2026: coll/rad norm split diagnostic.
! The preconditioner factors BA_full = BA_coll + BA_rad (the outer-
! iteration BA is restored at the end of this routine), so the
! interesting audit question is which part dominates and on which
! depths -- not whether the PC is missing radiative physics.
!
	REAL(KIND=LDP) NRM_COLL_F, NRM_RAD_F, NRM_COLL_MX, NRM_RAD_MX
	INTEGER K_RAD_DOM
!
	LUER = ERROR_LU()
	NRM_COLL_F = 0.0_LDP
	NRM_RAD_F = 0.0_LDP
	NRM_COLL_MX = 0.0_LDP
	NRM_RAD_MX = 0.0_LDP
	K_RAD_DOM = 0
	JFNK_BA_RAD_READY = .FALSE.
	JFNK_BA_RAD_BANDS_READY = .FALSE.
	JFNK_BA_DIAG_READY = .FALSE.
!
! Allocate banded BA_rad storage when JFNK_USE_BANDED=T.
! JFNK_CMAT_RAD_BANDS(NT, NT, NUM_BNDS, ND) stores the radiative BA
! block for each (band, depth) pair.  Only allocated when banded mode
! is requested to avoid memory overhead in the default DIAG path.
!
	IF(JFNK_USE_BANDED_LOC)THEN
	  IOS=0
	  IF(ALLOCATED(JFNK_CMAT_RAD_BANDS))THEN
	    IF(SIZE(JFNK_CMAT_RAD_BANDS,1).NE.JFNK_NT .OR.
	1       SIZE(JFNK_CMAT_RAD_BANDS,2).NE.JFNK_NT .OR.
	2       SIZE(JFNK_CMAT_RAD_BANDS,3).NE.JFNK_NUM_BNDS .OR.
	3       SIZE(JFNK_CMAT_RAD_BANDS,4).NE.JFNK_ND)THEN
	      DEALLOCATE(JFNK_CMAT_RAD_BANDS)
	    END IF
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_CMAT_RAD_BANDS))THEN
	    ALLOCATE(JFNK_CMAT_RAD_BANDS(JFNK_NT, JFNK_NT,
	1       JFNK_NUM_BNDS, JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0)THEN
	      WRITE(LUER,*)'JFNK_BUILD_BA_RAD_V1: cannot alloc CMAT_RAD_BANDS'
	      STOP 'JFNK_BUILD_BA_RAD_V1: alloc CMAT_RAD_BANDS failed'
	    END IF
	  END IF
	  JFNK_CMAT_RAD_BANDS = 0.0_LDP
	END IF
!
	IOS=0
	IF(ALLOCATED(JFNK_CMAT_RAD))THEN
	  IF(SIZE(JFNK_CMAT_RAD,1).NE.JFNK_NT .OR.
	1    SIZE(JFNK_CMAT_RAD,2).NE.JFNK_NT .OR.
	2    SIZE(JFNK_CMAT_RAD,3).NE.JFNK_ND)THEN
	    DEALLOCATE(JFNK_CMAT_RAD)
	  END IF
	END IF
	IF(.NOT. ALLOCATED(JFNK_CMAT_RAD))THEN
	  ALLOCATE(JFNK_CMAT_RAD(JFNK_NT, JFNK_NT, JFNK_ND), STAT=IOS)
	  IF(IOS.NE.0)THEN
	    WRITE(LUER,*)'JFNK_BUILD_BA_RAD_V1: cannot allocate CMAT_RAD'
	    STOP 'JFNK_BUILD_BA_RAD_V1: alloc CMAT_RAD failed'
	  END IF
	END IF
	JFNK_CMAT_RAD = 0.0_LDP
!
	IOS=0
	IF(ALLOCATED(JFNK_BA_DIAG))THEN
	  IF(SIZE(JFNK_BA_DIAG,1).NE.JFNK_NT .OR.
	1    SIZE(JFNK_BA_DIAG,2).NE.JFNK_ND)THEN
	    DEALLOCATE(JFNK_BA_DIAG)
	  END IF
	END IF
	IF(.NOT. ALLOCATED(JFNK_BA_DIAG))THEN
	  ALLOCATE(JFNK_BA_DIAG(JFNK_NT, JFNK_ND), STAT=IOS)
	  IF(IOS.NE.0)THEN
	    WRITE(LUER,*)'JFNK_BUILD_BA_RAD_V1: cannot allocate BA_DIAG'
	    STOP 'JFNK_BUILD_BA_RAD_V1: alloc BA_DIAG failed'
	  END IF
	END IF
	JFNK_BA_DIAG = 0.0_LDP
!
! NaN diagnostic: check SE(ID)%BA for NaN BEFORE GENERATE_FULL_MATRIX_V3
	  BLOCK
	    INTEGER :: ID_CHK, K_CHK, N_NAN_SEBA, N_TOTAL_SEBA
	    N_NAN_SEBA = 0; N_TOTAL_SEBA = 0
	    DO ID_CHK = 1, JFNK_NION
	      IF(ASSOCIATED(SE(ID_CHK)%BA))THEN
	        DO K_CHK = 1, SIZE(SE(ID_CHK)%BA,4)
	          IF(SE(ID_CHK)%BA(1,1,1,K_CHK) .NE.
	1               SE(ID_CHK)%BA(1,1,1,K_CHK))THEN
	            N_NAN_SEBA = N_NAN_SEBA + 1
	          END IF
	          N_TOTAL_SEBA = N_TOTAL_SEBA + 1
	        END DO
	      END IF
	    END DO
	    WRITE(LUER,'(A,I6,A,I6)')
	1      '  SE%BA NaN check (before GEN_FULL): NaN=', N_NAN_SEBA,
	1      ' depths out of ', N_TOTAL_SEBA
	    FLUSH(LUER)
	  END BLOCK
!
! --- Expand the full (radiative+collisional) BA into NT x NT blocks ---
!
! When JFNK_USE_BANDED=T, loop over all NUM_BNDS bands and store each
! NT×NT block in JFNK_CMAT_RAD_BANDS.  The diagonal band is always
! passed first (required by GENERATE_FULL_MATRIX_V3).  Off-diagonal
! bands carry cross-depth J-coupling from VJ (dJ/dopacity at
! neighbouring depths, written into SE(ID)%BA off-diagonal bands by
! the moment/ALO update in ba_update_v7) plus advection coupling
! (BA_ADV_TERM).
!
! When JFNK_USE_BANDED=F, only the diagonal band is extracted (legacy
! DIAG-only behaviour).
!
	ZERO_STEQ = .FALSE.
	REPLACE = .FALSE.
	USE_PASSED_REP = .FALSE.
	FIRST_MATRIX = .TRUE.
	LAST_MATRIX = .FALSE.
	DO K = 1, JFNK_ND
	  IF(K .EQ. JFNK_ND) LAST_MATRIX = .TRUE.
	  CALL GENERATE_FULL_MATRIX_V3(
	1       JFNK_C_MAT, DUMMY_STEQ, POPS_BASE, REPLACE, ZERO_STEQ,
	1       JFNK_NT, JFNK_ND, JFNK_NION, JFNK_NUM_BNDS,
	1       JFNK_DIAG_INDX, JFNK_DIAG_INDX, K,
	1       FIRST_MATRIX, LAST_MATRIX, USE_PASSED_REP)
	  JFNK_CMAT_RAD(:,:,K) = JFNK_C_MAT
	  IF(JFNK_USE_BANDED_LOC)THEN
	    JFNK_CMAT_RAD_BANDS(:,:,JFNK_DIAG_INDX,K) = JFNK_C_MAT
	  END IF
	  DO J = 1, JFNK_NT
	    JFNK_BA_DIAG(J, K) = JFNK_C_MAT(J, J)
	  END DO
!
! Zero NaN/Inf in BA_full before storing in CMAT_RAD.  The outer
! iteration's SE(ID)%BA has NaN at compression depths because
! VAR_MOM_J_CMF_V12 writes NaN to KI/BA when SIGMA_TRUE < -1.
! The mask zeros NaN in TX/TVX but not in SE(ID)%BA or KI.
! Zeroing NaN here is consistent with the mask's approach: at
! compression depths, the moment-equation J-variation is
! inconsistent with the non-monotonic formal solver, so the
! radiative J-coupling from the moment equation is discarded.
! The ALO (in TX/TVX) provides the physical depth coupling.
	  BLOCK
	    INTEGER :: N_NAN_ZN
	    N_NAN_ZN = COUNT(JFNK_C_MAT .NE. JFNK_C_MAT)
	    WHERE(JFNK_C_MAT .NE. JFNK_C_MAT) JFNK_C_MAT = 0.0_LDP
	    N_NAN_ZN = N_NAN_ZN +
	1              COUNT(ABS(JFNK_C_MAT) .GT. 1.0E30_LDP)
	    WHERE(ABS(JFNK_C_MAT) .GT. 1.0E30_LDP)
	      JFNK_C_MAT = 0.0_LDP
	    END WHERE
	    IF(N_NAN_ZN .GT. 0)THEN
	      JFNK_CMAT_RAD(:,:,K) = JFNK_C_MAT
	      IF(JFNK_USE_BANDED_LOC)THEN
	        JFNK_CMAT_RAD_BANDS(:,:,JFNK_DIAG_INDX,K) = JFNK_C_MAT
	      END IF
	      DO J = 1, JFNK_NT
	        JFNK_BA_DIAG(J, K) = JFNK_C_MAT(J, J)
	      END DO
	      IF(K .LE. 3 .OR. K .GE. JFNK_ND-2)THEN
	        WRITE(LUER,'(A,I3,A,I6)')
	1        '  BA_FULL zeroed NaN: K=', K,
	1        ' count=', N_NAN_ZN
	        FLUSH(LUER)
	      END IF
	    END IF
	  END BLOCK
	  FIRST_MATRIX = .FALSE.
!
! Off-diagonal bands for this depth (banded mode only).
! Interleaved with diagonal: GENERATE_FULL_MATRIX_V3 requires
! the diagonal band first at each depth (LST_DEPTH_INDX check).
! USE_PASSED_REP=.TRUE. skips equation-replacement (already done
! on the diagonal pass).
!
	  IF(JFNK_USE_BANDED_LOC .AND. JFNK_NUM_BNDS .GT. 1)THEN
	    DO IB = 1, JFNK_NUM_BNDS
	      IF(IB .EQ. JFNK_DIAG_INDX) CYCLE
	      CALL GENERATE_FULL_MATRIX_V3(
	1         JFNK_C_MAT, DUMMY_STEQ, POPS_BASE, REPLACE, ZERO_STEQ,
	1         JFNK_NT, JFNK_ND, JFNK_NION, JFNK_NUM_BNDS,
	1         IB, JFNK_DIAG_INDX, K,
	1         .FALSE., (K .EQ. JFNK_ND), .TRUE.)
! Zero NaN/Inf in off-diagonal BA_full bands (same reason as diagonal)
	      BLOCK
	        WHERE(JFNK_C_MAT .NE. JFNK_C_MAT)
	          JFNK_C_MAT = 0.0_LDP
	        END WHERE
	        WHERE(ABS(JFNK_C_MAT) .GT. 1.0E30_LDP)
	          JFNK_C_MAT = 0.0_LDP
	        END WHERE
	      END BLOCK
	      JFNK_CMAT_RAD_BANDS(:,:,IB,K) = JFNK_C_MAT
	    END DO
	  END IF
	END DO
	JFNK_BA_DIAG_READY = .TRUE.
!
! --- Save the full BA so we can restore it after the collisional rebuild ---
!
	ALLOCATE(BA_HOLD(JFNK_NION), STAT=IOS)
	IF(IOS.NE.0) STOP 'JFNK_BUILD_BA_RAD_V1: alloc BA_HOLD failed'
	DO ID = 1, JFNK_NION
	  IF(ASSOCIATED(SE(ID)%BA))THEN
	    N1 = SIZE(SE(ID)%BA,1)
	    N2 = SIZE(SE(ID)%BA,2)
	    N3 = SIZE(SE(ID)%BA,3)
	    N4 = SIZE(SE(ID)%BA,4)
	    ALLOCATE(BA_HOLD(ID)%BA(N1,N2,N3,N4), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_BUILD_BA_RAD_V1: alloc BA copy failed'
	    BA_HOLD(ID)%BA = SE(ID)%BA
	  END IF
	  IF(ASSOCIATED(SE(ID)%BA_PAR))THEN
	    N1 = SIZE(SE(ID)%BA_PAR,1)
	    N2 = SIZE(SE(ID)%BA_PAR,2)
	    N3 = SIZE(SE(ID)%BA_PAR,3)
	    ALLOCATE(BA_HOLD(ID)%BA_PAR(N1,N2,N3), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_BUILD_BA_RAD_V1: alloc BA_PAR copy failed'
	    BA_HOLD(ID)%BA_PAR = SE(ID)%BA_PAR
	  END IF
	END DO
	IF(ALLOCATED(BA_ED))THEN
	  ALLOCATE(BA_ED_HOLD(SIZE(BA_ED,1),SIZE(BA_ED,2),SIZE(BA_ED,3)))
	  BA_ED_HOLD = BA_ED
	END IF
	IF(ALLOCATED(BA_T))THEN
	  ALLOCATE(BA_T_HOLD(SIZE(BA_T,1),SIZE(BA_T,2),SIZE(BA_T,3)))
	  BA_T_HOLD = BA_T
	END IF
	IF(ALLOCATED(BA_T_EHB))THEN
	  ALLOCATE(BA_T_EHB_HOLD(SIZE(BA_T_EHB,1),SIZE(BA_T_EHB,2),
	1   SIZE(BA_T_EHB,3)))
	  BA_T_EHB_HOLD = BA_T_EHB
	END IF
	IF(ALLOCATED(BA_T_PAR))THEN
	  ALLOCATE(BA_T_PAR_HOLD(SIZE(BA_T_PAR,1),SIZE(BA_T_PAR,2)))
	  BA_T_PAR_HOLD = BA_T_PAR
	END IF
	IF(ALLOCATED(BA_T_PAR_EHB))THEN
	  ALLOCATE(BA_T_PAR_EHB_HOLD(SIZE(BA_T_PAR_EHB,1),
	1   SIZE(BA_T_PAR_EHB,2)))
	  BA_T_PAR_EHB_HOLD = BA_T_PAR_EHB
	END IF
	IF(ALLOCATED(BA_ADV_TERM))THEN
	  ALLOCATE(BA_ADV_HOLD(SIZE(BA_ADV_TERM,1),SIZE(BA_ADV_TERM,2)))
	  BA_ADV_HOLD = BA_ADV_TERM
	END IF
!
! --- Rebuild a collisional-only BA at baseline populations ---
!
! Fix #3 (09-Sep-2026): JFNK_COMPUTE_STEQ_V1 unconditionally zeroes and
! recomputes SE(ID)%STEQ, STEQ_T, STEQ_ED, STEQ_T_EHB (and STEQ_ADV /
! T_EHB).  Without save/restore, those module-level globals are left
! holding the collisional-only residual at baseline populations,
! corrupting any code that inspects them before the next full outer
! STEQ computation.  Shared helper, same bracket as JFNK_RESIDUAL_V1.
!
	CALL JFNK_SAVE_MOD_STEQ_V1()
!
	CALL JFNK_SYNC_ATM_FROM_POPS(POPS_BASE)
	CALL JFNK_COMPUTE_STEQ_V1(JFNK_NT, JFNK_ND, JFNK_NION,
	1     JFNK_NUM_BNDS, JFNK_DIAG_INDX, POPS_BASE,
	1     .TRUE., JFNK_FIXED_T_LOC, .FALSE.)
!
	ZERO_STEQ = .FALSE.
	REPLACE = .FALSE.
	USE_PASSED_REP = .FALSE.
	FIRST_MATRIX = .TRUE.
	LAST_MATRIX = .FALSE.
	DO K = 1, JFNK_ND
	  IF(K .EQ. JFNK_ND) LAST_MATRIX = .TRUE.
	  CALL GENERATE_FULL_MATRIX_V3(
	1       JFNK_C_MAT, DUMMY_STEQ, POPS_BASE, REPLACE, ZERO_STEQ,
	1       JFNK_NT, JFNK_ND, JFNK_NION, JFNK_NUM_BNDS,
	1       JFNK_DIAG_INDX, JFNK_DIAG_INDX, K,
	1       FIRST_MATRIX, LAST_MATRIX, USE_PASSED_REP)
	  JFNK_CMAT_RAD(:,:,K) = JFNK_CMAT_RAD(:,:,K) - JFNK_C_MAT
	  IF(JFNK_USE_BANDED_LOC)THEN
	    JFNK_CMAT_RAD_BANDS(:,:,JFNK_DIAG_INDX,K) =
	1       JFNK_CMAT_RAD_BANDS(:,:,JFNK_DIAG_INDX,K) - JFNK_C_MAT
	  END IF
!
! Accumulate coll/rad norms for the audit summary.  JFNK_C_MAT is
! the collisional block; CMAT_RAD holds BA_full - BA_coll after the
! subtraction above.
!
	  NRM_COLL_F = NRM_COLL_F + SUM(JFNK_C_MAT**2)
	  NRM_RAD_F = NRM_RAD_F + SUM(JFNK_CMAT_RAD(:,:,K)**2)
	  NRM_COLL_MX = MAX(NRM_COLL_MX, MAXVAL(ABS(JFNK_C_MAT)))
	  NRM_RAD_MX = MAX(NRM_RAD_MX, MAXVAL(ABS(JFNK_CMAT_RAD(:,:,K))))
	  IF(MAXVAL(ABS(JFNK_CMAT_RAD(:,:,K))) .GT.
	1     MAXVAL(ABS(JFNK_C_MAT))) K_RAD_DOM = K_RAD_DOM + 1
!
! F(0) diagnostic: check CMAT_RAD (BA_full - BA_coll) for NaN
	  IF(K .LE. 3 .OR. K .GE. JFNK_ND-2)THEN
	    BLOCK
	      INTEGER :: N_NAN_CR
	      N_NAN_CR = COUNT(JFNK_CMAT_RAD(:,:,K) .NE.
	1                   JFNK_CMAT_RAD(:,:,K))
      IF(N_NAN_CR .GT. 0)THEN
        WRITE(LUER,'(A,I3,A,I6,A,I6)')
	1      '  CMAT_RAD NaN check: K=', K,
	1      ' NaN=', N_NAN_CR, ' out of ', JFNK_NT*JFNK_NT
        FLUSH(LUER)
      END IF
	    END BLOCK
	  END IF
	  FIRST_MATRIX = .FALSE.
!
! Off-diagonal bands for collisional-only BA (banded mode).
! Interleaved with diagonal (GENERATE_FULL_MATRIX_V3 requires
! diagonal band first at each depth).
!
	  IF(JFNK_USE_BANDED_LOC .AND. JFNK_NUM_BNDS .GT. 1)THEN
	    DO IB = 1, JFNK_NUM_BNDS
	      IF(IB .EQ. JFNK_DIAG_INDX) CYCLE
	      CALL GENERATE_FULL_MATRIX_V3(
	1         JFNK_C_MAT, DUMMY_STEQ, POPS_BASE, REPLACE, ZERO_STEQ,
	1         JFNK_NT, JFNK_ND, JFNK_NION, JFNK_NUM_BNDS,
	1         IB, JFNK_DIAG_INDX, K,
	1         .FALSE., (K .EQ. JFNK_ND), .TRUE.)
	      JFNK_CMAT_RAD_BANDS(:,:,IB,K) =
	1       JFNK_CMAT_RAD_BANDS(:,:,IB,K) - JFNK_C_MAT
	    END DO
	  END IF
	END DO
	JFNK_BA_RAD_BANDS_READY = JFNK_USE_BANDED_LOC
!
! --- Restore the full BA for the GMRES preconditioner ---
!
	DO ID = 1, JFNK_NION
	  IF(ALLOCATED(BA_HOLD(ID)%BA)) SE(ID)%BA = BA_HOLD(ID)%BA
	  IF(ALLOCATED(BA_HOLD(ID)%BA_PAR))THEN
	    SE(ID)%BA_PAR = BA_HOLD(ID)%BA_PAR
	  END IF
	END DO
	IF(ALLOCATED(BA_ED_HOLD)) BA_ED = BA_ED_HOLD
	IF(ALLOCATED(BA_T_HOLD)) BA_T = BA_T_HOLD
	IF(ALLOCATED(BA_T_EHB_HOLD)) BA_T_EHB = BA_T_EHB_HOLD
	IF(ALLOCATED(BA_T_PAR_HOLD)) BA_T_PAR = BA_T_PAR_HOLD
	IF(ALLOCATED(BA_T_PAR_EHB_HOLD)) BA_T_PAR_EHB = BA_T_PAR_EHB_HOLD
	IF(ALLOCATED(BA_ADV_HOLD)) BA_ADV_TERM = BA_ADV_HOLD
!
! Restore STEQ module-level arrays (Fix #3, 09-Sep-2026; shared
! helper covers STEQ_ADV / T_EHB as of 03-Oct-2026).
!
	CALL JFNK_RESTORE_MOD_STEQ_V1()
!
	DEALLOCATE(BA_HOLD)
	IF(ALLOCATED(BA_ED_HOLD)) DEALLOCATE(BA_ED_HOLD)
	IF(ALLOCATED(BA_T_HOLD)) DEALLOCATE(BA_T_HOLD)
	IF(ALLOCATED(BA_T_EHB_HOLD)) DEALLOCATE(BA_T_EHB_HOLD)
	IF(ALLOCATED(BA_T_PAR_HOLD)) DEALLOCATE(BA_T_PAR_HOLD)
	IF(ALLOCATED(BA_T_PAR_EHB_HOLD)) DEALLOCATE(BA_T_PAR_EHB_HOLD)
	IF(ALLOCATED(BA_ADV_HOLD)) DEALLOCATE(BA_ADV_HOLD)
!
	JFNK_BA_RAD_READY = .TRUE.
	JFNK_CMAT_AGE = 0
	WRITE(LUER,'(A)')
	1  '  JFNK: radiative Jacobian CMAT_RAD = BA_full - BA_coll'
	WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1  '  JFNK: ||BA_coll||_F/max =', SQRT(NRM_COLL_F),
	1  ' /', NRM_COLL_MX
	WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1  '  JFNK: ||BA_rad ||_F/max =', SQRT(NRM_RAD_F),
	1  ' /', NRM_RAD_MX
	WRITE(LUER,'(A,I4,A,I4,A)')
	1  '  JFNK: rad-dominated depth blocks:', K_RAD_DOM,
	1  ' of', JFNK_ND, ' (max-element basis)'
	FLUSH(LUER)
!
	END SUBROUTINE JFNK_BUILD_BA_RAD_V1


!
! JFNK_BA_RAD_MATVEC_V1: y = CMAT_RAD * x_frac
!
! x is the JFNK state (log pops / ED / T under [JFNK_LOGX]=T, else
! fractional pops / ED with log T).  Every column is converted to the
! actual fractional change exp(x)-1 (or x in the linear map) so the BA
! scaling CMAT*(dn/n, dT/T) matches the residual's n(x) map.
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_BA_RAD_MATVEC_V1(X, Y)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_log_unknown
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN)  :: X(JFNK_NT*JFNK_ND)
	DOUBLE PRECISION, INTENT(OUT) :: Y(JFNK_NT*JFNK_ND)
	REAL(KIND=LDP) XFRAC(JFNK_NT)
	REAL(KIND=LDP) YLOC(JFNK_NT)
	INTEGER K, J, IB, L
!
	Y = 0.0D0
	IF(.NOT. JFNK_BA_RAD_READY) RETURN
!
! Banded matvec: y = sum over bands of CMAT_RAD_BANDS(:,:,B,K) * x(L)
! where L = BND_TO_FULL(B,K) is the variable depth for band B at
! observation depth K.  Cross-depth J-coupling enters through the
! off-diagonal bands (B .NE. DIAG_INDX).
!
! BND_TO_FULL(B,K) = K + B - DIAG_INDX
! For NUM_BNDS=3, DIAG_INDX=2: B=1→L=K-1, B=2→L=K, B=3→L=K+1.
! Fix #4 (review 2026-09-11): the old formula used integer division
! (NUM_BNDS/ND)*(DIAG_INDX-K) which is always 0 for NUM_BNDS << ND.
! It worked by accident for the 3-band case but is not a valid
! general band mapping.  The correct offset is simply K + B -
! DIAG_INDX, assuming DIAG_INDX = (NUM_BNDS+1)/2 = BAND_HALF + 1.
!
	IF(JFNK_BA_RAD_BANDS_READY)THEN
	  DO K = 1, JFNK_ND
	    YLOC = 0.0_LDP
	    DO IB = 1, JFNK_NUM_BNDS
	      L = K + IB - JFNK_DIAG_INDX
	      IF(L .LT. 1 .OR. L .GT. JFNK_ND) CYCLE
! CMAT_RAD is d(STEQ)/d(frac), so XFRAC must be the actual fractional
! change.  Under [JFNK_LOGX]=T the residual map is n=n0*exp(x), giving
! frac=exp(x)-1 -- using x directly (as the T column once did for all)
! is only first-order correct and errs by (exp(x)-1)/x at finite x.
! The +/-40 clip matches the residual's population bound.
	      DO J = 1, JFNK_NT - 1
		IF(jfnk_log_unknown)THEN
		  XFRAC(J) = REAL(EXP(MIN(MAX(X((L-1)*JFNK_NT + J),
	1           -40.0D0), 40.0D0)) - 1.0D0, KIND=LDP)
		ELSE
		  XFRAC(J) = REAL(X((L-1)*JFNK_NT + J), KIND=LDP)
		END IF
	      END DO
	      IF(JFNK_FIXED_NE_LOC) XFRAC(JFNK_NT-1) = 0.0_LDP
	      IF(JFNK_FIXED_T_LOC)THEN
		XFRAC(JFNK_NT) = 0.0_LDP
	      ELSE IF(JFNK_FIX_T_D_ST .GT. 0 .AND.
	1         L .GE. JFNK_FIX_T_D_ST .AND.
	2         L .LE. JFNK_FIX_T_D_END .AND.
	3         .NOT. JFNK_TAUTO_ACTIVE)THEN
		XFRAC(JFNK_NT) = 0.0_LDP
	      ELSE
		XFRAC(JFNK_NT) =
	1       REAL(EXP(MIN(MAX(X((L-1)*JFNK_NT + JFNK_NT),
	2           -50.0D0), 50.0D0)) - 1.0D0, KIND=LDP)
	      END IF
	      YLOC = YLOC + MATMUL(JFNK_CMAT_RAD_BANDS(:,:,IB,K), XFRAC)
	    END DO
	    Y((K-1)*JFNK_NT+1 : K*JFNK_NT) = DBLE(YLOC)
	  END DO
	  RETURN
	END IF
!
! DIAG-only matvec (legacy): y(K) = CMAT_RAD(:,:,K) * x(K)
! Same frac conversion as the banded path above.
!
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT - 1
	    IF(jfnk_log_unknown)THEN
	      XFRAC(J) = REAL(EXP(MIN(MAX(X((K-1)*JFNK_NT + J),
	1           -40.0D0), 40.0D0)) - 1.0D0, KIND=LDP)
	    ELSE
	      XFRAC(J) = REAL(X((K-1)*JFNK_NT + J), KIND=LDP)
	    END IF
	  END DO
	  IF(JFNK_FIXED_NE_LOC) XFRAC(JFNK_NT-1) = 0.0_LDP
	  IF(JFNK_FIXED_T_LOC)THEN
	    XFRAC(JFNK_NT) = 0.0_LDP
	  ELSE IF(JFNK_FIX_T_D_ST .GT. 0 .AND.
	1         K .GE. JFNK_FIX_T_D_ST .AND.
	2         K .LE. JFNK_FIX_T_D_END .AND.
	3         .NOT. JFNK_TAUTO_ACTIVE)THEN
	    XFRAC(JFNK_NT) = 0.0_LDP
	  ELSE
	    XFRAC(JFNK_NT) =
	1       REAL(EXP(MIN(MAX(X((K-1)*JFNK_NT + JFNK_NT),
	2           -50.0D0), 50.0D0)) - 1.0D0, KIND=LDP)
	  END IF
	  Y((K-1)*JFNK_NT+1 : K*JFNK_NT) =
	1     DBLE(MATMUL(JFNK_CMAT_RAD(:,:,K), XFRAC))
	END DO
!
	END SUBROUTINE JFNK_BA_RAD_MATVEC_V1


!
! JFNK_COMPUTE_STEQ_V1: Recompute SE(ID)%STEQ from the current POPS
! with J frozen.  Uses module-level ED, T from MOD_CMFGEN.
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_COMPUTE_STEQ_V1(NT, ND, NION, NUM_BNDS, DIAG_INDX,
	1                               POPS_IN, COMPUTE_BA_IN,
	1                               FIXED_T_IN, LST_ITERATION_IN)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN, ONLY: ED, T, ATM, NUM_IONS, ION_ID,
	1   POP_SPECIES, SPECIES_LNK, NUM_SPECIES
	USE STEQ_DATA_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: DO_CO_MOV_DDT, SN_MODEL,
	1   ADVEC_RELAX_PARAM, LINEAR_ADV, INCL_ADVECTION, TIME_SEQ_NO,
	1   LAMBDA_ITERATION
	USE LINE_MOD, ONLY: AVE_ENERGY
	USE JFNK_CMFGEN_INT_MOD, ONLY: JFNK_DIERECOM, JFNK_DIECOOL,
	1   JFNK_W_CNM, JFNK_W_DCNM, JFNK_N_S_MAX, JFNK_SKIP_SE_ZERO
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: NT, ND, NION, NUM_BNDS, DIAG_INDX
	REAL(KIND=LDP), INTENT(IN) :: POPS_IN(NT, ND)
	LOGICAL, INTENT(IN) :: COMPUTE_BA_IN, FIXED_T_IN, LST_ITERATION_IN
!
	INTEGER ID, K, DST, DEND, I, J
	CHARACTER(LEN=80) TMP_STRING
	INTEGER LOC_ID
	REAL(KIND=LDP) T1, T2
	INTEGER LUER, ERROR_LU
	EXTERNAL ERROR_LU
	EXTERNAL OMEGA_GEN_V3
!
	LUER = ERROR_LU()
!
! Zero STEQ arrays (always recomputed from scratch).
! BA arrays are ONLY zeroed when COMPUTE_BA_IN=.TRUE., because when
! COMPUTE_BA_IN=.FALSE. (used during JFNK residual evaluations), the BA
! matrix must be preserved for use by the GMRES preconditioner.
! CORRECTION (2026-09-07): Original code unconditionally zeroed BA arrays,
! destroying the preconditioner after each residual evaluation.
!
! JFNK_SKIP_SE_ZERO (coupled mode, 02-Oct-2026): when TRUE the caller
! has already zeroed the arrays and run the radiation sweep -- the
! radiative terms it added to SE%STEQ/STEQ_T must be preserved and the
! collisional terms below accumulated on top.
!
	IF(.NOT. JFNK_SKIP_SE_ZERO)THEN
!$OMP PARALLEL DO
	  DO ID = 1, NION
	    SE(ID)%STEQ   = 0.0_LDP
	    SE(ID)%T_EHB  = 0.0_LDP
	  END DO
!$OMP END PARALLEL DO
	  STEQ_ED = 0.0_LDP
	  STEQ_T = 0.0_LDP
	  STEQ_T_EHB = 0.0_LDP
	END IF
! DIERECOM/DIECOOL are scratch accumulators for STEQ_AUTO_V2 below --
! always zeroed, in both modes.
	JFNK_DIERECOM = 0.0_LDP
	JFNK_DIECOOL = 0.0_LDP
!
	IF(COMPUTE_BA_IN)THEN
!$OMP PARALLEL DO
	  DO ID = 1, NION
	    SE(ID)%BA     = 0.0_LDP
	    SE(ID)%BA_PAR = 0.0_LDP
	  END DO
!$OMP END PARALLEL DO
	  BA_ED = 0.0_LDP
	  BA_T = 0.0_LDP
	  BA_T_PAR = 0.0_LDP
	  BA_T_EHB = 0.0_LDP
	  BA_T_PAR_EHB = 0.0_LDP
	END IF
!
	DST = 1
	DEND = ND
!
	DO ID = 1, NUM_IONS - 1
	  LOC_ID = ID
	  IF(ATM(ID)%XzV_PRES)THEN
	    TMP_STRING = TRIM(ION_ID(ID))//'_COL_DATA'
	    CALL STEQ_MULTI_V10(JFNK_W_CNM, JFNK_W_DCNM, ED, T,
	1       ATM(ID)%XzV,            ATM(ID)%XzVLTE,         ATM(ID)%dlnXzVLTE_dlnT,
	1       AVE_ENERGY(ATM(ID)%EQXzV),
	1       ATM(ID)%NXzV,           ATM(ID)%DXzV,           ATM(ID)%XzV_F,
	1       ATM(ID)%XzVLTE_F_ON_S,  ATM(ID)%W_XzV_F,        ATM(ID)%AXzV_F,
	1       ATM(ID)%EDGEXzV_F,      ATM(ID)%GXzV_F,         ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%NXzV_F,         ATM(ID)%F_TO_S_XzV,
	1       POP_SPECIES(1,SPECIES_LNK(ID)), ATM(ID+1)%XzV_PRES, ATM(ID)%ZXzV,
	1       LOC_ID, TMP_STRING, OMEGA_GEN_V3,
	1       ATM(ID)%EQXzV, NUM_BNDS, ND, NION, NT,
	1       COMPUTE_BA_IN, FIXED_T_IN, LST_ITERATION_IN, DST, DEND)
!
	    TMP_STRING = TRIM(ION_ID(ID))//'_AUTO_DATA'
	    CALL STEQ_AUTO_V2(ED, T,
	1       ATM(ID)%XzV,        ATM(ID)%NXzV,         ATM(ID)%DXzV,
	1       ATM(ID)%XzV_F,      ATM(ID)%XzVLTE_F,     ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%GXzV_F,     ATM(ID)%XzVLEVNAME_F, ATM(ID)%NXzV_F,
	1       ATM(ID)%F_TO_S_XzV, LOC_ID,
	1       JFNK_DIERECOM(1,ATM(ID)%INDX_XzV),
	1       JFNK_DIECOOL(1,ATM(ID)%INDX_XzV),
	1       TMP_STRING,         NUM_BNDS, ND, COMPUTE_BA_IN, DST, DEND)
	  END IF
	END DO
!
	CALL STEQNE_V4(ED, NT, DIAG_INDX, ND, COMPUTE_BA_IN, DST, DEND)
!
! CORRECTION (2026-09-07): STEQ_BA_TWO_PHOT_RATE_V3 has 7 arguments:
!   (POPS, NT, ND, DIAG_INDX, UPDATE_BA, LU_OUT, WRITE_RATES)
! The original call passed only 6 arguments, with LST_ITERATION_IN
! (a LOGICAL) mapped to LU_OUT (an INTEGER), and WRITE_RATES missing.
! This is a type mismatch and missing-argument bug that could cause
! runtime crashes or corrupt I/O.  Fixed to pass LUER and .FALSE.
!
	CALL STEQ_BA_TWO_PHOT_RATE_V3(POPS_IN, NT, ND,
	1       DIAG_INDX, COMPUTE_BA_IN, LUER, .FALSE.)
!
	CALL STEQ_BA_CHG_EXCH_V3(POPS_IN, T, NT, ND, DIAG_INDX,
	1   COMPUTE_BA_IN)
!
	IF(DO_CO_MOV_DDT .AND. SN_MODEL)THEN
	  CALL STEQ_CO_MOV_DERIV_V3(POPS_IN, ADVEC_RELAX_PARAM, LINEAR_ADV,
	1       DO_CO_MOV_DDT, LAMBDA_ITERATION, COMPUTE_BA_IN,
	1       TIME_SEQ_NO, NUM_BNDS, ND, NT)
	ELSE
	  CALL STEQ_ADVEC_V4(ADVEC_RELAX_PARAM, LINEAR_ADV, NUM_BNDS, ND,
	1       INCL_ADVECTION, LAMBDA_ITERATION, COMPUTE_BA_IN)
	END IF
!
	END SUBROUTINE JFNK_COMPUTE_STEQ_V1


!
! JFNK_SET_PTC_SCALE_V1: freeze the pseudo-time weight s_i = max(|F_i|, floor)
! from the baseline residual.  s is held fixed for the rest of the
! Newton loop so the finite-difference derivative of (x - x0)*s/dtau is
! exactly s/dtau.  Called after F(0) is adopted, and again when a
! coupled rebuild moves the baseline.
!
! edited by Devin, 06-Oct-2026 (SER anchor fix): also (re)set the
! pseudo-time anchor jfnk_ptc_x0 to 0 -- the baseline x = 0 is the
! first anchor, and newton_hook_v1 re-anchors at each accepted iterate
! (Kelley & Keyes 1998; Knoll & Keyes 2004 Sec 4.2).
!
	SUBROUTINE JFNK_SET_PTC_SCALE_V1(N, F0)
	USE SET_KIND_MODULE
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_ptc_scale, jfnk_ptc_scale_ready,
	1     jfnk_ptc_x0, jfnk_ptc_anchor_ready
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(IN) :: F0(N)
	INTEGER I, IOS
	DOUBLE PRECISION FMAX, FLOOR_USE
!
	IF(ALLOCATED(jfnk_ptc_scale))THEN
	  IF(SIZE(jfnk_ptc_scale) .NE. N) DEALLOCATE(jfnk_ptc_scale)
	END IF
	IF(.NOT. ALLOCATED(jfnk_ptc_scale))THEN
	  ALLOCATE(jfnk_ptc_scale(N), STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    jfnk_ptc_scale_ready = .FALSE.
	    jfnk_ptc_anchor_ready = .FALSE.
	    RETURN
	  END IF
	END IF
	IF(ALLOCATED(jfnk_ptc_x0))THEN
	  IF(SIZE(jfnk_ptc_x0) .NE. N) DEALLOCATE(jfnk_ptc_x0)
	END IF
	IF(.NOT. ALLOCATED(jfnk_ptc_x0))THEN
	  ALLOCATE(jfnk_ptc_x0(N), STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    jfnk_ptc_scale_ready = .FALSE.
	    jfnk_ptc_anchor_ready = .FALSE.
	    RETURN
	  END IF
	END IF
	FMAX = 0.0D0
	DO I = 1, N
	  FMAX = MAX(FMAX, ABS(F0(I)))
	END DO
	FLOOR_USE = MAX(1.0D-30, 1.0D-8 * FMAX)
	DO I = 1, N
	  jfnk_ptc_scale(I) = REAL(MAX(ABS(F0(I)), FLOOR_USE), KIND=LDP)
	END DO
	jfnk_ptc_scale_ready = .TRUE.
	jfnk_ptc_x0 = 0.0D0
	jfnk_ptc_anchor_ready = .TRUE.
	END SUBROUTINE JFNK_SET_PTC_SCALE_V1
!
!
! JFNK_RESIDUAL_V1: frozen-J statistical-equilibrium residual.
!
!   F(x) = STEQ_coll(n0*(1+x)) + STEQ_rad(n0) + CMAT_RAD * x_frac
!
! STEQ_coll is recomputed from the perturbed populations with J frozen.
! STEQ_rad(n0) is the constant leftover from the outer iteration.
! CMAT_RAD * x_frac is the linearized radiative increment BA_rad * x
! (CMAT already scaled by POPS so x is a fractional correction;
! temperature uses x_frac(NT) = exp(x_T)-1).
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_RESIDUAL_V1(N, X, Y)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_x, jfnk_log_unknown,
	1     jfnk_ptc_active, jfnk_ptc_dtau, jfnk_ptc_scale,
	2     jfnk_ptc_scale_ready, jfnk_ptc_x0, jfnk_ptc_anchor_ready,
	3     jfnk_resid_nbad, jfnk_maxch_est
	USE MOD_CMFGEN, ONLY: ED, T
	USE STEQ_DATA_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_MAX_EVALS
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(IN) :: X(N)
	DOUBLE PRECISION, INTENT(OUT) :: Y(N)
!
	INTEGER K, J, N_BAD, N_NAN_Y, N_INF_Y, IDX, ID_ION
	INTEGER LUER
	INTEGER, EXTERNAL :: ERROR_LU
	REAL(KIND=LDP) MCE_DABS, MCE_DMAX
!
! Audit 17-Sep-2026 (issue 4): count constraint crossings that are NEW
! relative to the baseline Newton point jfnk_x -- i.e. a clip that is
! inactive at jfnk_x but active at the perturbed point X.  The residual
! is nonsmooth at the population floor and the T-exp cap, so a FD JVP
! that crosses such a boundary no longer approximates the smooth
! Jacobian.  JFNK_N_NEW_CLIP is read by JFNK_JVP_V1, which retries the
! finite difference with a smaller step.
!
	INTEGER N_NEW_CLIP_LOC
	REAL(KIND=LDP) T_BASE_UNCL, X_TRY, X_BASE
!
! FRAC, STEQ_MAT, and YRAD are now persistent module workspace
! (JFNK_W_FRAC, JFNK_W_STEQ_MAT, JFNK_W_YRAD) to avoid stack
! allocation.  audit §3, 09-Sep-2026.
!
	LUER = ERROR_LU()
!
! Eval budget (VADAT [JFNK_MAX_EVALS], default 0=off): a coupled
! eval is a full formal sweep, so an activation that exhausts the
! budget is spending more than a failed solve is worth.  Return a
! sentinel residual (no sweep): GMRES aborts on the JVP NaN-count
! path and trial steps see |F| huge, so the Newton loop rejects,
! the trust region collapses, and the driver takes the
! best-iterate/Lambda fallback instead of burning sweeps.  The
! counter resets in JFNK_INIT_V1, so the budget is per-activation.
!
	IF(JFNK_MAX_EVALS .GT. 0 .AND.
	1     JFNK_N_RESIDUAL_EVALS .GE. JFNK_MAX_EVALS)THEN
	  IF(JFNK_N_RESIDUAL_EVALS .EQ. JFNK_MAX_EVALS)THEN
	    WRITE(LUER,'(A,I6,A)')
	1     '  JFNK: eval budget (',JFNK_MAX_EVALS,
	2     ') exhausted -- forcing rejection path'
	    FLUSH(LUER)
	  END IF
	  Y=1.0D20
	  jfnk_resid_nbad=N
	  JFNK_N_RESIDUAL_EVALS = JFNK_N_RESIDUAL_EVALS + 1
	  RETURN
	END IF
!
! Coupled tier without a registered sweep callback cannot evaluate
! the residual at all: the frozen-J pieces (JFNK_STEQ_RAD_FLAT,
! CMAT_RAD) are only built in frozen mode, so the previous
! collisional-only fallback evaluated a different operator and the
! coupled-mode F(0) adoption masked the mismatch (order-of-operations
! audit, 03-Oct-2026).  Treat as a broken eval -- same sentinel the
! exhausted-budget path returns, before any state is touched, so the
! driver aborts on F(0) and takes the LAMBDA fallback.
!
	IF(JFNK_COUPLED_LOC .AND. .NOT. ASSOCIATED(JFNK_SWEEP_CB))THEN
	  WRITE(LUER,'(A)')
	1   '  JFNK: coupled eval admitted but JFNK_SWEEP_CB is not'//
	2   ' registered -- eval impossible, forcing fallback'
	  FLUSH(LUER)
	  jfnk_resid_nbad=N
	  JFNK_N_RESIDUAL_EVALS = JFNK_N_RESIDUAL_EVALS + 1
	  Y=1.0D20
	  RETURN
	END IF
!
! Task 3: Parallelize depth-wise loops.  Each K writes a distinct
! column of FRAC / POPS_WORK / STEQ_MAT and a distinct block of Y,
! so there are no inter-iteration data dependencies.
!
!$OMP PARALLEL DO PRIVATE(J)
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT
	    JFNK_W_FRAC(J, K) = REAL(X((K-1)*JFNK_NT + J), KIND=LDP)
	  END DO
	END DO
!$OMP END PARALLEL DO
!
! Alias-safe copy of X: JFNK_BA_RAD_MATVEC_V1 (below) reads X after Y
! has been written, so if X and Y alias (as in the F(0) identity test),
! the BA_RAD term would be corrupted.  Save X into module workspace
! immediately, before any Y writes occur.  Fix #5 (09-Sep-2026).
!
	DO IDX = 1, N
	  JFNK_W_X_SAVE(IDX) = REAL(X(IDX), KIND=LDP)
	END DO
!
! Save current state (ED/T from MOD_CMFGEN, POPS_WORK from module)
!
	POPS_SAVE = POPS_WORK
	ED_SAVE = ED
	T_SAVE = T
	N_NEW_CLIP_LOC = 0
!
! Apply fractional correction to POPS_WORK
!
! For populations (indices 1..NT-2) and electron density (NT-1):
!   POPS_new = POPS_base * (1 + x)
! with a floor at 1.0E-100 to prevent negative/zero populations.
!
! For temperature (index NT), when JFNK_FIXED_T_LOC=.FALSE.:
!   T_new = T_base * exp(x)    (log-transform, added 2026-09-07)
!
! The log-transform is used for temperature because:
!   1. It guarantees T > 0 for any finite x (no negative temperatures)
!   2. It treats relative changes symmetrically: a factor of 2 increase
!      (x = ln 2 ≈ 0.69) is the same magnitude as a factor of 2 decrease
!      (x = -ln 2 ≈ -0.69)
!   3. At shock fronts, temperature can change by orders of magnitude
!      between adjacent depth points.  The linear formulation T*(1+x)
!      requires |x| < 1 to keep T positive, which limits the correction
!      to a factor of 2.  The log-transform exp(x) allows arbitrarily
!      large corrections while maintaining positivity.
!   4. The Newton correction x is naturally O(1) in log space because
!      the STEQ for T is d(ln T)/d(ln T) = 1, so the Jacobian is
!      well-scaled.
!
! The POPS(NT,K) entry stores T in the same array as populations, so
! we apply the log-transform to POPS(NT,K) as well.
!
!$OMP PARALLEL DO PRIVATE(J, T_BASE_UNCL, X_TRY, X_BASE) REDUCTION(+:N_NEW_CLIP_LOC)
	DO K = 1, JFNK_ND
! Populations (J=1..NT-2) and Ne (J=NT-1).
! [JFNK_LOGX] (default T): x = ln(n/n0), so n = n0*exp(x).
! The map stays positive without the 1e-100 floor, which made F
! non-differentiable wherever a trial step crossed it.  The +/-40
! bound only stops exp from overflowing; a crossing of that bound
! is counted so the finite-difference step is retried.
! [JFNK_LOGX]=F keeps n = n0*(1+x) with the floor.
!
	  DO J = 1, JFNK_NT - 2
	    IF(jfnk_log_unknown)THEN
	      X_TRY = JFNK_W_FRAC(J, K)
	      X_BASE = REAL(jfnk_x((K-1)*JFNK_NT + J), KIND=LDP)
	      IF(ABS(X_TRY) .GT. 40.0_LDP .AND. ABS(X_BASE) .LE. 40.0_LDP)
	1       N_NEW_CLIP_LOC = N_NEW_CLIP_LOC + 1
	      X_TRY = MIN(MAX(X_TRY, -40.0_LDP), 40.0_LDP)
	      POPS_WORK(J, K) = POPS_BASE(J, K) * EXP(X_TRY)
	    ELSE
	      POPS_WORK(J, K) = POPS_BASE(J, K) * (1.0_LDP + JFNK_W_FRAC(J, K))
	      IF(POPS_WORK(J, K) .LT. 1.0E-100_LDP)THEN
	        T_BASE_UNCL = POPS_BASE(J, K) *
	1         (1.0_LDP + REAL(jfnk_x((K-1)*JFNK_NT + J),KIND=LDP))
	        IF(T_BASE_UNCL .GE. 1.0E-100_LDP)
	1         N_NEW_CLIP_LOC = N_NEW_CLIP_LOC + 1
	        POPS_WORK(J, K) = 1.0E-100_LDP
	      END IF
	    END IF
	  END DO
! ED correction (J=NT-1): skip when FIX_NE (singular BA row otherwise).
	  IF(JFNK_FIXED_NE_LOC)THEN
	    POPS_WORK(JFNK_NT-1, K) = POPS_BASE(JFNK_NT-1, K)
	    ED(K) = ED_SAVE(K)
	  ELSE IF(jfnk_log_unknown)THEN
	    X_TRY = JFNK_W_FRAC(JFNK_NT-1, K)
	    X_BASE = REAL(jfnk_x((K-1)*JFNK_NT + JFNK_NT-1), KIND=LDP)
	    IF(ABS(X_TRY) .GT. 40.0_LDP .AND. ABS(X_BASE) .LE. 40.0_LDP)
	1     N_NEW_CLIP_LOC = N_NEW_CLIP_LOC + 1
	    X_TRY = MIN(MAX(X_TRY, -40.0_LDP), 40.0_LDP)
	    POPS_WORK(JFNK_NT-1, K) = POPS_BASE(JFNK_NT-1, K) * EXP(X_TRY)
	    ED(K) = POPS_WORK(JFNK_NT-1, K)
	  ELSE
	    POPS_WORK(JFNK_NT-1, K) = POPS_BASE(JFNK_NT-1, K) *
	1       (1.0_LDP + JFNK_W_FRAC(JFNK_NT-1, K))
	    IF(POPS_WORK(JFNK_NT-1, K) .LT. 1.0E-100_LDP)THEN
	      T_BASE_UNCL = POPS_BASE(JFNK_NT-1, K) *
	1       (1.0_LDP +
	2       REAL(jfnk_x((K-1)*JFNK_NT + JFNK_NT-1),KIND=LDP))
	      IF(T_BASE_UNCL .GE. 1.0E-100_LDP)
	1       N_NEW_CLIP_LOC = N_NEW_CLIP_LOC + 1
	      POPS_WORK(JFNK_NT-1, K) = 1.0E-100_LDP
	    END IF
	    ED(K) = POPS_WORK(JFNK_NT-1, K)
	  END IF
!
! V4 audit fix (2026-09-07): Support partial FIXED_T depth range
! (audit finding #11).  T is fixed when:
!   - JFNK_FIXED_T_LOC=.TRUE. (global), OR
!   - JFNK_FIX_T_D_ST > 0 AND K is in [ST, END] (partial VARFIXT range)
	  IF(JFNK_FIXED_T_LOC)THEN
	    POPS_WORK(JFNK_NT, K) = POPS_BASE(JFNK_NT, K)
	    T(K) = T_SAVE(K)
	  ELSE IF(JFNK_FIX_T_D_ST .GT. 0 .AND.
	1         K .GE. JFNK_FIX_T_D_ST .AND.
	2         K .LE. JFNK_FIX_T_D_END .AND.
	3         .NOT. JFNK_TAUTO_ACTIVE)THEN
	    POPS_WORK(JFNK_NT, K) = POPS_BASE(JFNK_NT, K)
	    T(K) = T_SAVE(K)
	  ELSE
!   Log-transform for temperature: T_new = T_base * exp(x)
!   This replaces the old T_new = T_base * (1 + x) which could produce
!   negative temperatures when |x| > 1 (common at shock fronts).
!   The exp() guarantees T > 0 for any finite x.
	    IF(ABS(JFNK_W_FRAC(JFNK_NT, K)) .GT. 50.0_LDP .AND.
	1       ABS(jfnk_x((K-1)*JFNK_NT + JFNK_NT)) .LE. 50.0D0)
	2       N_NEW_CLIP_LOC = N_NEW_CLIP_LOC + 1
	    POPS_WORK(JFNK_NT, K) = POPS_BASE(JFNK_NT, K) *
	1       EXP(MIN(MAX(JFNK_W_FRAC(JFNK_NT, K), -50.0_LDP), 50.0_LDP))
	    IF(POPS_WORK(JFNK_NT, K) .LT. 1.0_LDP)
	1       POPS_WORK(JFNK_NT, K) = 1.0_LDP
	    T(K) = POPS_WORK(JFNK_NT, K)
	  END IF
	END DO
!$OMP END PARALLEL DO
!
! Expose the new-constraint-crossing count to JFNK_JVP_V1 (issue 4).
!
	JFNK_N_NEW_CLIP = N_NEW_CLIP_LOC
!
	CALL JFNK_SYNC_ATM_FROM_POPS(POPS_WORK)
!
! Save STEQ arrays before recomputation (added 2026-09-07):
! JFNK_COMPUTE_STEQ_V1 modifies SE(ID)%STEQ, STEQ_T, STEQ_ED, and
! STEQ_T_EHB.  These are module-level state in STEQ_DATA_MOD that may
! be read by other code paths (Ng acceleration, diagnostics, the outer
! great iteration) between JFNK evaluations.  We save the module-level
! 1D arrays and the flat assembled STEQ vector here, and restore them
! after the residual is assembled, so the perturbed STEQ values do not
! leak to those code paths.
!
! Note: The per-ion SE(ID)%STEQ arrays have variable dimensions (NX
! varies by ion), so they are packed into a flat buffer inside
! JFNK_SAVE_MOD_STEQ_V1 along with STEQ_ADV / T_EHB.  This makes the
! residual evaluation fully side-effect-free (Finding #4, 2026-09-08;
! extended 03-Oct-2026 for the advective/EHB arrays).
!
	CALL JFNK_SAVE_MOD_STEQ_V1()
!
! Coupled residual (V7 port, 02-Oct-2026; Arramy et al. 2024 style):
! re-run the formal radiation sweep at the perturbed populations so
! F(x) contains the full pops -> opacity/emissivity -> J -> rates -> SE
! chain, and the finite-difference JVP probes dJ/dn.  The lagged Knop
! I_NEXT store is frozen (FREEZE_I_NEXT=T, set by the driver) so the
! sweep is a deterministic function of the trial populations -- the
! store is outer-iteration state, refreshed by the next real sweep
! after the correction is committed.
!
! Ordering: the sweep accumulates radiative terms into SE%STEQ (which
! the zeroing below clears first) and fills the TWO_PHOT_MOD rates that
! STEQ_BA_TWO_PHOT_RATE_V3 consumes inside JFNK_COMPUTE_STEQ_V1, so the
! sweep MUST run before the compute_steq call.  JFNK_SKIP_SE_ZERO then
! keeps compute_steq's collisional/additive terms on top of the
! radiative ones.  The frozen-J additions (JFNK_STEQ_RAD_FLAT and the
! CMAT_RAD matvec) are skipped automatically: JFNK_RAD_READY and
! JFNK_BA_RAD_READY are only set in frozen mode.
!
! JFNK_SWEEP_CB is guaranteed associated here: the guard at the top
! of this routine already forced a sentinel return when the coupled
! tier was admitted without a registered callback.
!
	IF(JFNK_COUPLED_LOC)THEN
!$OMP PARALLEL DO
	  DO ID_ION = 1, JFNK_NION
	    IF(ASSOCIATED(SE(ID_ION)%STEQ))SE(ID_ION)%STEQ = 0.0_LDP
	    IF(ASSOCIATED(SE(ID_ION)%T_EHB))SE(ID_ION)%T_EHB = 0.0_LDP
	  END DO
!$OMP END PARALLEL DO
	  STEQ_ED = 0.0_LDP
	  STEQ_T = 0.0_LDP
	  STEQ_T_EHB = 0.0_LDP
	  CALL JFNK_SWEEP_CB(POPS_WORK)
	  JFNK_SKIP_SE_ZERO = .TRUE.
	  CALL JFNK_COMPUTE_STEQ_V1(JFNK_NT, JFNK_ND, JFNK_NION,
	1       JFNK_NUM_BNDS, JFNK_DIAG_INDX, POPS_WORK,
	1       .FALSE., JFNK_FIXED_T_LOC, .FALSE.)
	  JFNK_SKIP_SE_ZERO = .FALSE.
	ELSE
	  CALL JFNK_COMPUTE_STEQ_V1(JFNK_NT, JFNK_ND, JFNK_NION,
	1       JFNK_NUM_BNDS, JFNK_DIAG_INDX, POPS_WORK,
	1       .FALSE., JFNK_FIXED_T_LOC, .FALSE.)
	END IF
!
	CALL JFNK_ASSEMBLE_STEQ_V1(JFNK_W_STEQ_MAT, JFNK_NT, JFNK_ND,
	1     JFNK_NION, JFNK_NUM_BNDS, JFNK_DIAG_INDX, POPS_WORK)
!
! STEQ_MAT range diagnostic.  Printed once on the first residual eval
! (baseline residual scale) and again only when a NaN/Inf is present —
! not on every JVP eval, which flooded OUTGEN (audit, 2026-09-11).
	BLOCK
	  INTEGER :: N_NAN_SM, N_INF_SM, JJ_SM, KK_SM
	  REAL(KIND=LDP) :: MAX_SM, MIN_SM, VAL_SM
	  N_NAN_SM = 0; N_INF_SM = 0
	  MAX_SM = -HUGE(1.0_LDP); MIN_SM = HUGE(1.0_LDP)
	  DO KK_SM = 1, JFNK_ND
	    DO JJ_SM = 1, JFNK_NT
	      VAL_SM = JFNK_W_STEQ_MAT(JJ_SM, KK_SM)
	      IF(VAL_SM .NE. VAL_SM)THEN
	        N_NAN_SM = N_NAN_SM + 1
	      ELSE IF(ABS(VAL_SM) .GT. 1.0E30_LDP)THEN
	        N_INF_SM = N_INF_SM + 1
	      ELSE
	        MAX_SM = MAX(MAX_SM, VAL_SM)
	        MIN_SM = MIN(MIN_SM, VAL_SM)
	      END IF
	    END DO
	  END DO
	  IF(JFNK_N_RESIDUAL_EVALS .EQ. 0 .OR. N_NAN_SM+N_INF_SM .GT. 0)THEN
	    WRITE(LUER,'(A,I6,A,I6,A,ES12.4,A,ES12.4)')
	1      '  RESID STEQ_MAT: NaN=', N_NAN_SM,
	1      ' Inf=', N_INF_SM,
	1      ' max=', MAX_SM, ' min=', MIN_SM
	    FLUSH(LUER)
	  END IF
	END BLOCK
!
! Task 3: Fused parallel loop — assemble STEQ into Y and add the
! radiative STEQ leftover in a single pass over depth points.
! Each K writes a distinct block of Y, so there are no race conditions.
!
!$OMP PARALLEL DO PRIVATE(J)
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT
	    Y((K-1)*JFNK_NT + J) = DBLE(JFNK_W_STEQ_MAT(J, K))
	    IF(JFNK_RAD_READY .AND. .NOT. JFNK_COUPLED_LOC)THEN
	      Y((K-1)*JFNK_NT + J) = Y((K-1)*JFNK_NT + J) +
	1         JFNK_STEQ_RAD_FLAT((K-1)*JFNK_NT + J)
	    END IF
	  END DO
	END DO
!$OMP END PARALLEL DO
	IF(JFNK_BA_RAD_READY .AND. .NOT. JFNK_COUPLED_LOC)THEN
	  CALL JFNK_BA_RAD_MATVEC_V1(JFNK_W_X_SAVE, JFNK_W_YRAD)
!
! YRAD range diagnostic: same gating as STEQ_MAT — baseline once, then
! only on NaN/Inf.
	  BLOCK
	    INTEGER :: N_NAN_YR, N_INF_YR, JJ_YR
	    REAL(KIND=LDP) :: MAX_YR, MIN_YR, VAL_YR
	    N_NAN_YR = 0; N_INF_YR = 0
	    MAX_YR = -HUGE(1.0_LDP); MIN_YR = HUGE(1.0_LDP)
	    DO JJ_YR = 1, N
	      VAL_YR = JFNK_W_YRAD(JJ_YR)
	      IF(VAL_YR .NE. VAL_YR)THEN
	        N_NAN_YR = N_NAN_YR + 1
	      ELSE IF(ABS(VAL_YR) .GT. 1.0E30_LDP)THEN
	        N_INF_YR = N_INF_YR + 1
	      ELSE
	        MAX_YR = MAX(MAX_YR, VAL_YR)
	        MIN_YR = MIN(MIN_YR, VAL_YR)
	      END IF
	    END DO
	    IF(JFNK_N_RESIDUAL_EVALS .EQ. 0 .OR. N_NAN_YR+N_INF_YR .GT. 0)THEN
	      WRITE(LUER,'(A,I6,A,I6,A,ES12.4,A,ES12.4)')
	1      '  RESID YRAD: NaN=', N_NAN_YR,
	1      ' Inf=', N_INF_YR,
	1      ' max=', MAX_YR, ' min=', MIN_YR
	      FLUSH(LUER)
	    END IF
	  END BLOCK
	  DO IDX = 1, N
	    Y(IDX) = Y(IDX) + JFNK_W_YRAD(IDX)
	  END DO
	END IF
!
! CMFGEN MAXCH analogue: the largest relative correction still
! demanded by any population row.  |F_i|/|BA_ii| prices row i in
! log-population units, so converged rows (the F-dominated HI levels
! that carry the Euclidean norm) drop out of this metric.  Done here,
! before the pseudo-time term is added, so the metric sees F rather
! than G.  Non-finite entries are skipped; the NaN guard below caps
! them for the solver.
!
	jfnk_maxch_est = -1.0_LDP
	IF(JFNK_BA_DIAG_READY)THEN
	  jfnk_maxch_est = 0.0_LDP
	  DO K = 1, JFNK_ND
	    MCE_DMAX = 0.0_LDP
	    DO J = 1, JFNK_NT
	      MCE_DMAX = MAX(MCE_DMAX, ABS(JFNK_BA_DIAG(J,K)))
	    END DO
	    DO J = 1, JFNK_NT-2
	      IDX = (K-1)*JFNK_NT + J
	      IF(IDX .GT. N) CYCLE
	      IF(Y(IDX) .NE. Y(IDX)) CYCLE
	      MCE_DABS = MAX(ABS(JFNK_BA_DIAG(J,K)),
	1        1.0E-30_LDP*MCE_DMAX)
	      jfnk_maxch_est = MAX(jfnk_maxch_est,
	1        ABS(Y(IDX))/MCE_DABS)
	    END DO
	  END DO
	END IF
!
! Pseudo-time term (Knoll & Keyes 2004 Sec 4.2; [JFNK_PTC], default T).
! G_i = F_i + (x_i - x0_i) * s_i / dtau, with s_i frozen from |F| at the
! baseline and x0 = jfnk_ptc_x0 the anchor -- the current Newton iterate
! x_k (SER form, Kelley & Keyes 1998; edited by Devin, 06-Oct-2026).
! The term vanishes at the anchor, so G(x_k) = F(x_k): the residual
! newton_hook_v1 stores is F after jfnk_ptc_strip_v1, and all Newton
! norms act on |F|.  Previously the term was x_i*s_i/dtau -- anchored
! at the solve baseline x = 0 -- which pulled later iterates back and
! could report |G| convergence with F /= 0.  A finite-difference JVP of
! G is still J plus the constant diagonal s/dtau, so the JFNK_PRECOND_V1
! diagonal correction is unchanged.  Large dtau drops the term and G is
! F.  Added before the fixed-unknown mask so those rows stay zero.
! Skipped until JFNK_SET_PTC_SCALE_V1 has run (the F(0) identity test
! is one such call); a missing/mismatched anchor falls back to x0 = 0.
!
	IF(jfnk_ptc_active .AND. jfnk_ptc_scale_ready)THEN
	  IF(jfnk_ptc_anchor_ready .AND. ALLOCATED(jfnk_ptc_x0) .AND.
	1     SIZE(jfnk_ptc_x0) .EQ. N)THEN
	    DO IDX = 1, N
	      Y(IDX) = Y(IDX) + (X(IDX) - DBLE(jfnk_ptc_x0(IDX))) *
	1         DBLE(jfnk_ptc_scale(IDX)) /
	2         MAX(jfnk_ptc_dtau, 1.0D-8)
	    END DO
	  ELSE
	    DO IDX = 1, N
	      Y(IDX) = Y(IDX) + X(IDX) *
	1         DBLE(jfnk_ptc_scale(IDX)) /
	2         MAX(jfnk_ptc_dtau, 1.0D-8)
	    END DO
	  END IF
	END IF
!
! Fixed T / Ne are not Newton unknowns.  Zero those residual components
! so GMRES does not try to satisfy a singular BA row.
! editted by Sean Brennan (08-Sep-2026)
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, Y)
!
! NaN/Inf policy: do not zero individual components (GMRES would treat
! those equations as already converged).  Cap each bad component at
! +/-1D20 (sign kept where defined) so |F| stays huge and the Hookstep
! trust region shrinks, while the finite components keep their values.
! The previous whole-vector overwrite (Y=1D20 everywhere) erased the
! finite entries, made the F(0) test report rel=1 for every failure,
! and hid which equation actually blew up.  jfnk_resid_nbad records
! the count for the driver so the F(0) test can abort on a broken
! eval rather than on a large lagged-field offset.
!
	N_BAD = 0
	N_NAN_Y = 0
	N_INF_Y = 0
!$OMP SIMD REDUCTION(+:N_BAD,N_NAN_Y,N_INF_Y)
	DO IDX = 1, N
	  IF(Y(IDX) .NE. Y(IDX))THEN
	    N_BAD = N_BAD + 1
	    N_NAN_Y = N_NAN_Y + 1
	  ELSE IF(ABS(Y(IDX)) .GT. 1.0D30)THEN
	    N_BAD = N_BAD + 1
	    N_INF_Y = N_INF_Y + 1
	  END IF
	END DO
!$OMP END SIMD
	jfnk_resid_nbad = N_BAD
	IF(N_BAD .GT. 0)THEN
	  JFNK_N_RESIDUAL_NAN = JFNK_N_RESIDUAL_NAN + 1
	  WRITE(LUER,'(A,I6,A,I6,A,I6)')
	1    '  RESID NaN guard: bad=', N_BAD,
	1    ' NaN=', N_NAN_Y, ' Inf=', N_INF_Y
	  FLUSH(LUER)
	  DO IDX = 1, N
	    IF(Y(IDX) .NE. Y(IDX))THEN
	      Y(IDX) = 1.0D20
	    ELSE IF(ABS(Y(IDX)) .GT. 1.0D30)THEN
	      Y(IDX) = SIGN(1.0D20, Y(IDX))
	    END IF
	  END DO
	END IF
!
! Restore state: POPS_WORK, ED, T, ATM, and STEQ module-level arrays.
! This ensures that after every residual evaluation, the module-level
! state is identical to what it was before, making F(x) evaluation
! side-effect-free (pure) with respect to the global CMFGEN state.
!
	POPS_WORK = POPS_SAVE
	ED = ED_SAVE
	T = T_SAVE
!
! Restore STEQ module-level 1D arrays to their pre-evaluation values.
! Also restore the per-ion SE(ID)%STEQ / STEQ_ADV / T_EHB arrays from
! the flat buffers (Finding #4, 2026-09-08; extended 03-Oct-2026).
!
	CALL JFNK_RESTORE_MOD_STEQ_V1()
!
	CALL JFNK_SYNC_ATM_FROM_POPS(POPS_WORK)
!
	JFNK_N_RESIDUAL_EVALS = JFNK_N_RESIDUAL_EVALS + 1
!
	END SUBROUTINE JFNK_RESIDUAL_V1


!
! JFNK_REBUILD_BASELINE_V1: Apply a Newton correction to POPS_BASE and
! recompute the collisional STEQ, linearized radiative STEQ, BA_RAD,
! preconditioner factorization, and RSCALE so the next JFNK attempt
! starts from a genuinely new physical state.
!
! Fix #2 (review 2026-09-11): the previous retry loop rebuilt the
! preconditioner from the SAME POPS_BASE, producing the same BA.  This
! subroutine applies the best correction found so far to POPS_BASE
! (with the same clipping as the final population update), then
! recomputes everything, giving the next attempt new information.
!
! The radiative STEQ is linearized: STEQ_RAD_new = STEQ_RAD_old +
! CMAT_RAD_old * x_corr.  After the update, CMAT_RAD is rebuilt from
! the new POPS_BASE via JFNK_BUILD_BA_RAD_V1.
!
	SUBROUTINE JFNK_REBUILD_BASELINE_V1(X_CORR, LUER)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_x, jfnk_fx, jfnk_log_unknown,
	1     jfnk_ptc_active, jfnk_ptc_dtau, jfnk_fx_fresh,
	2     jfnk_ptc_scale_ready
	USE CONTROL_VARIABLE_MOD, ONLY: MAX_JFNK_COR, MAX_dT_COR,
	1                                JFNK_FIXED_T, FIXED_NE
	USE MOD_CMFGEN, ONLY: ED, T
	USE STEQ_DATA_MOD, ONLY: SE
	IMPLICIT NONE
!
	DOUBLE PRECISION, INTENT(IN) :: X_CORR(*)
	INTEGER, INTENT(IN) :: LUER
	REAL(KIND=LDP) X_MAX, X_MIN, X_T_MAX, X_T_MIN
	REAL(KIND=LDP) T1_CORR
	DOUBLE PRECISION X_CORR_FULL(JFNK_NT*JFNK_ND)
	INTEGER K, J, IDX
	INTEGER ID_ION, OFF
	DOUBLE PRECISION, EXTERNAL :: JFNK_DOT_V1
!
! The baseline populations are about to move: the next baseline eval
! re-harvests the equation structure for the new state (see
! JFNK_ASSEMBLE_STEQ_V1).
!
	JFNK_EQ_STRUCT_READY=.FALSE.
!
! Clipping parameters (must match the final update in
! solve_for_pops_jfnk_v1.f)
!
	X_MAX = MAX_JFNK_COR - 1.0_LDP
	X_MIN = 1.0_LDP/MAX_JFNK_COR - 1.0_LDP
	IF(jfnk_log_unknown)THEN
	  X_MAX = LOG(MAX_JFNK_COR)
	  X_MIN = -X_MAX
	END IF
	X_T_MAX = LOG(MAX_JFNK_COR)
	X_T_MIN = -X_T_MAX
	IF(MAX_dT_COR .GT. 0.0_LDP .AND. MAX_dT_COR .LT. 1.0_LDP)THEN
	  X_T_MAX = MIN(X_T_MAX, LOG(1.0_LDP + MAX_dT_COR))
	  X_T_MIN = MAX(X_T_MIN, LOG(1.0_LDP - MAX_dT_COR))
	END IF
!
! Step 1: Apply clipped correction to POPS_BASE
!
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT-2
	    T1_CORR = REAL(X_CORR((K-1)*JFNK_NT + J), KIND=LDP)
	    IF(T1_CORR .GT. X_MAX) T1_CORR = X_MAX
	    IF(T1_CORR .LT. X_MIN) T1_CORR = X_MIN
	    IF(jfnk_log_unknown)THEN
	      POPS_BASE(J, K) = POPS_BASE(J, K) * EXP(T1_CORR)
	    ELSE
	      POPS_BASE(J, K) = POPS_BASE(J, K) * (1.0_LDP + T1_CORR)
	      IF(POPS_BASE(J, K) .LT. 1.0E-100_LDP)
	1       POPS_BASE(J, K) = 1.0E-100_LDP
	    END IF
	  END DO
! ED (NT-1)
	  IF(.NOT. FIXED_NE)THEN
	    T1_CORR = REAL(X_CORR((K-1)*JFNK_NT + JFNK_NT-1), KIND=LDP)
	    IF(T1_CORR .GT. X_MAX) T1_CORR = X_MAX
	    IF(T1_CORR .LT. X_MIN) T1_CORR = X_MIN
	    IF(jfnk_log_unknown)THEN
	      POPS_BASE(JFNK_NT-1, K) = POPS_BASE(JFNK_NT-1, K) * EXP(T1_CORR)
	    ELSE
	      POPS_BASE(JFNK_NT-1, K) = POPS_BASE(JFNK_NT-1, K) *
	1         (1.0_LDP + T1_CORR)
	      IF(POPS_BASE(JFNK_NT-1, K) .LT. 1.0E-100_LDP)
	1       POPS_BASE(JFNK_NT-1, K) = 1.0E-100_LDP
	    END IF
	  END IF
	  ED(K) = POPS_BASE(JFNK_NT-1, K)
! Temperature (NT)
	  IF(.NOT. JFNK_FIXED_T_LOC)THEN
	    IF(JFNK_FIX_T_D_ST .GT. 0 .AND.
	1       K .GE. JFNK_FIX_T_D_ST .AND.
	2       K .LE. JFNK_FIX_T_D_END .AND.
	3       .NOT. JFNK_TAUTO_ACTIVE)THEN
! T held fixed in this depth range
	    ELSE
	      T1_CORR = REAL(X_CORR((K-1)*JFNK_NT + JFNK_NT), KIND=LDP)
	      IF(T1_CORR .GT. X_T_MAX) T1_CORR = X_T_MAX
	      IF(T1_CORR .LT. X_T_MIN) T1_CORR = X_T_MIN
	      POPS_BASE(JFNK_NT, K) = POPS_BASE(JFNK_NT, K) *
	1       EXP(T1_CORR)
	      IF(POPS_BASE(JFNK_NT, K) .LT. 1.0_LDP)
	1     POPS_BASE(JFNK_NT, K) = 1.0_LDP
	    END IF
	  END IF
	  T(K) = POPS_BASE(JFNK_NT, K)
	END DO
!
! Step 2: Sync atmosphere from updated POPS_BASE
!
! POPS_WORK must follow POPS_BASE here: JFNK_FACTOR_PRECOND_V1 (step 7)
! builds the BA blocks with POPS_WORK as the scaling argument, and the
! residual restores ATM from POPS_WORK between evals.  Left at the old
! baseline it silently factors the preconditioner at pre-rebuild
! populations.  The coupled branch already does this copy.
!
	CALL JFNK_SYNC_ATM_FROM_POPS(POPS_BASE)
	POPS_WORK = POPS_BASE
!
! Coupled mode ([JFNK_COUPLED]=T): F is the radiation sweep at
! POPS_BASE.  The frozen-J steps below replace F with collisions
! only, then set JFNK_BA_RAD_READY.  The next coupled eval would add
! that matvec on top of a sweep that already contains the radiative
! rates, so Newton differentiates the wrong operator and the step
! that gets committed does not reduce the real residual.  Re-evaluate
! the coupled residual instead.  POPS_WORK must match POPS_BASE first:
! the residual restores ATM from POPS_WORK on the way out.
!
	IF(JFNK_COUPLED_LOC)THEN
	  POPS_WORK=POPS_BASE
	  X_CORR_FULL=0.0D0
!
! The baseline eval must measure pure F at the new POPS_BASE.  The
! PTC anchor still holds the previous Newton attempt's last accepted
! iterate, so an armed PTC term would add (0 - x0_old)*s/dtau to
! jfnk_fx here and the JVP numerators would inherit the offset
! divided by h.  JFNK_SET_PTC_SCALE_V1 below re-arms with the new
! scale and a zero anchor.
!
	  jfnk_ptc_scale_ready=.FALSE.
	  CALL JFNK_RESIDUAL_V1(JFNK_NT*JFNK_ND,X_CORR_FULL,jfnk_fx)
	  CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_NT*JFNK_ND,jfnk_fx)
	  CALL JFNK_COMPUTE_RSCALE_V1(JFNK_NT*JFNK_ND,jfnk_fx)
	  jfnk_x=0.0D0
	  jfnk_fx_fresh=.TRUE.
	  IF(jfnk_ptc_active)THEN
	    CALL JFNK_SET_PTC_SCALE_V1(JFNK_NT*JFNK_ND,jfnk_fx)
	    jfnk_ptc_dtau=1.0D0
	  END IF
	  IF(LUER .GT. 0)THEN
	    WRITE(LUER,'(A)')
	1     '  JFNK_REBUILD: coupled baseline, F from eval sweep'
	    WRITE(LUER,'(A,ES12.4)')'    New |F(0)|_w = ',
	1       SQRT(JFNK_DOT_V1(JFNK_NT*JFNK_ND,jfnk_fx,jfnk_fx))
	    FLUSH(LUER)
	  END IF
	  RETURN
	END IF
!
! Step 3: Linearized radiative STEQ update:
! STEQ_RAD_new = STEQ_RAD_old + CMAT_RAD_old * X_CORR
! (must use OLD CMAT_RAD before rebuild)
!
	IF(JFNK_BA_RAD_READY)THEN
	  DO IDX = 1, JFNK_NT * JFNK_ND
	    X_CORR_FULL(IDX) = X_CORR(IDX)
	  END DO
	  CALL JFNK_BA_RAD_MATVEC_V1(X_CORR_FULL, JFNK_W_YRAD)
	  DO IDX = 1, JFNK_NT * JFNK_ND
	    JFNK_STEQ_RAD_FLAT(IDX) = JFNK_STEQ_RAD_FLAT(IDX) +
	1       REAL(JFNK_W_YRAD(IDX), KIND=LDP)
	  END DO
	END IF
!
! Step 4: Recompute collisional STEQ at new POPS_BASE.  Bracketed:
! compute_steq rewrites the module STEQ arrays; leaving the
! collisional-only state in place would give a fallback SOLVEBA a
! right-hand side missing every radiative term (order-of-operations
! audit, 03-Oct-2026).  The module arrays therefore stay at the outer
! iteration's full STEQ; only jfnk_fx tracks the new baseline.
!
	CALL JFNK_SAVE_MOD_STEQ_V1()
	CALL JFNK_COMPUTE_STEQ_V1(JFNK_NT, JFNK_ND, JFNK_NION,
	1     JFNK_NUM_BNDS, JFNK_DIAG_INDX, POPS_BASE,
	1     .FALSE., JFNK_FIXED_T_LOC, .FALSE.)
!
	CALL JFNK_ASSEMBLE_STEQ_V1(JFNK_W_STEQ_MAT, JFNK_NT, JFNK_ND,
	1     JFNK_NION, JFNK_NUM_BNDS, JFNK_DIAG_INDX, POPS_BASE)
	CALL JFNK_RESTORE_MOD_STEQ_V1()
!
! Step 5: Assemble new jfnk_fx = STEQ_coll + STEQ_RAD
!
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT
	    jfnk_fx((K-1)*JFNK_NT + J) = DBLE(JFNK_W_STEQ_MAT(J, K))
	    IF(JFNK_RAD_READY)THEN
	      jfnk_fx((K-1)*JFNK_NT + J) = jfnk_fx((K-1)*JFNK_NT + J) +
	1       DBLE(JFNK_STEQ_RAD_FLAT((K-1)*JFNK_NT + J))
	    END IF
	  END DO
	END DO
!
! Step 6: Rebuild BA_RAD from new POPS_BASE (updates CMAT_RAD)
!
	CALL JFNK_BUILD_BA_RAD_V1()
!
! Step 7: Factor preconditioner from new BA
!
	CALL JFNK_FACTOR_PRECOND_V1()
!
! Step 8: Recompute RSCALE and mask fixed unknowns
!
	CALL JFNK_MASK_FIXED_UNKNOWNS(JFNK_NT * JFNK_ND, jfnk_fx)
	CALL JFNK_COMPUTE_RSCALE_V1(JFNK_NT * JFNK_ND, jfnk_fx)
!
! Step 9: Reset Newton correction to zero (new baseline)
!
	jfnk_x = 0.0D0
	jfnk_fx_fresh = .TRUE.
!
! editted by Devin, 05-Oct-2026: refresh the pseudo-time scale on the
! new baseline, mirroring the coupled branch above.  jfnk_ptc_scale
! is frozen from |F(0)| of the baseline it was measured on; after
! POPS_BASE moves, keeping the old scale would leave sigma_i =
! s_i/dtau tied to the pre-rebuild residual.  dtau restarts at 1
! for the new Newton sequence on this baseline.
!
	IF(jfnk_ptc_active)THEN
	  CALL JFNK_SET_PTC_SCALE_V1(JFNK_NT*JFNK_ND,jfnk_fx)
	  jfnk_ptc_dtau=1.0D0
	END IF
!
	IF(LUER .GT. 0)THEN
	  WRITE(LUER,'(A)')'  JFNK_REBUILD: baseline updated with best &
	1correction'
	  WRITE(LUER,'(A,ES12.4)')'    New |F(0)|_w = ',
	1     SQRT(JFNK_DOT_V1(JFNK_NT*JFNK_ND, jfnk_fx, jfnk_fx))
	  FLUSH(LUER)
	END IF
!
	END SUBROUTINE JFNK_REBUILD_BASELINE_V1


!
! JFNK_TEMP_COL_AUDIT_V1: Validate the exp(T) transform and its BA
! derivative against finite differences.
!
! Fix #10 (review 2026-09-11): The residual uses T_new = T_base * exp(x)
! but the BA matrix may have been built for the linear formulation
! T_new = T_base * (1+x).  At x=0 both give the same Jacobian
! (dT/dx = T_base), but the BA might be missing a factor of T if it
! was built as dSTEQ/dT rather than dSTEQ/dT * T.
!
! This test:
! 1. Computes F(0) and F(x_T * e_T) for a small x_T
! 2. FD Jacobian T-column: J_FD = (F(x_T) - F(0)) / x_T
! 3. BA Jacobian T-column: J_BA = JFNK_BA_A(:, NT, K)
! 4. Reports max |J_FD - J_BA| / max(|J_FD|, |J_BA|) per depth
!
! A relative error >> 1 indicates the BA T-column is inconsistent with
! the exp(T) transform (e.g. missing a factor of T_base).
!
	SUBROUTINE JFNK_TEMP_COL_AUDIT_V1(LUER)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_fx
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: LUER
!
	DOUBLE PRECISION, ALLOCATABLE :: X_ZERO(:), X_PERT(:)
	DOUBLE PRECISION, ALLOCATABLE :: F_ZERO(:), F_PERT(:)
	DOUBLE PRECISION, ALLOCATABLE :: J_FD(:), J_BA(:)
	REAL(KIND=LDP) :: X_T, REL_ERR, MAX_FD, MAX_BA, MAX_ERR
	REAL(KIND=LDP) :: GLOBAL_MAX_ERR, GLOBAL_MAX_FD
	INTEGER K, J, N, N_BAD, N_CHECKED
	INTEGER, SAVE :: AUDIT_DONE = 0
!
	IF(AUDIT_DONE .GT. 0) RETURN
	IF(.NOT. ALLOCATED(JFNK_BA_A)) RETURN
!
	AUDIT_DONE = 1
	N = JFNK_NT * JFNK_ND
	X_T = 1.0D-5
!
	ALLOCATE(X_ZERO(N), X_PERT(N), F_ZERO(N), F_PERT(N))
	ALLOCATE(J_FD(N), J_BA(N))
!
! Step 1: F(0)
!
	X_ZERO = 0.0D0
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, X_ZERO)
	CALL JFNK_RESIDUAL_V1(N, X_ZERO, F_ZERO)
!
! Step 2: F(x_T * e_T) — perturb only the temperature column (index NT)
!
	X_PERT = 0.0D0
	DO K = 1, JFNK_ND
	  X_PERT((K-1)*JFNK_NT + JFNK_NT) = X_T
	END DO
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, X_PERT)
	CALL JFNK_RESIDUAL_V1(N, X_PERT, F_PERT)
!
! Step 3: FD Jacobian T-column and BA T-column
!
	GLOBAL_MAX_ERR = 0.0_LDP
	GLOBAL_MAX_FD = 0.0_LDP
	N_BAD = 0
	N_CHECKED = 0
!
	DO K = 1, JFNK_ND
	  MAX_FD = 0.0_LDP
	  MAX_BA = 0.0_LDP
	  MAX_ERR = 0.0_LDP
	  DO J = 1, JFNK_NT
	    J_FD((K-1)*JFNK_NT + J) =
	1     (F_PERT((K-1)*JFNK_NT + J) - F_ZERO((K-1)*JFNK_NT + J)) / X_T
	    J_BA((K-1)*JFNK_NT + J) = DBLE(JFNK_BA_A(J, JFNK_NT, K))
	    MAX_FD = MAX(MAX_FD, ABS(REAL(J_FD((K-1)*JFNK_NT + J),KIND=LDP)))
	    MAX_BA = MAX(MAX_BA, ABS(REAL(J_BA((K-1)*JFNK_NT + J),KIND=LDP)))
	    MAX_ERR = MAX(MAX_ERR,
	1     ABS(REAL(J_FD((K-1)*JFNK_NT + J) -
	1     J_BA((K-1)*JFNK_NT + J), KIND=LDP)))
	  END DO
	  IF(MAX_FD .GT. 1.0E-20_LDP .OR. MAX_BA .GT. 1.0E-20_LDP)THEN
	    N_CHECKED = N_CHECKED + 1
	    REL_ERR = MAX_ERR / MAX(MAX_FD, MAX_BA, 1.0E-30_LDP)
	    IF(REL_ERR .GT. 0.1_LDP) N_BAD = N_BAD + 1
	    GLOBAL_MAX_ERR = MAX(GLOBAL_MAX_ERR, MAX_ERR)
	    GLOBAL_MAX_FD = MAX(GLOBAL_MAX_FD, MAX_FD)
	  END IF
	END DO
!
	WRITE(LUER,'(/,A)')' =========================== T-COLUMN AUDIT (fix #10) =========================='
	WRITE(LUER,'(A,ES10.2,A,I5,A,I5)')
	1   '  x_T=', X_T, '  depths checked=', N_CHECKED, '/', JFNK_ND
	WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1   '  max|J_FD|=', GLOBAL_MAX_FD,
	1   '  max|J_FD-J_BA|=', GLOBAL_MAX_ERR
	IF(N_BAD .GT. 0)THEN
	  WRITE(LUER,'(A,I5,A)')
	1     '  WARNING: ', N_BAD, ' depths have >10% T-column mismatch'
	  WRITE(LUER,'(A)')
	1     '  BA T-column may be inconsistent with exp(T) transform'
	ELSE
	  WRITE(LUER,'(A)')'  T-column BA derivative matches FD (OK)'
	END IF
	FLUSH(LUER)
!
	DEALLOCATE(X_ZERO, X_PERT, F_ZERO, F_PERT, J_FD, J_BA)
!
	END SUBROUTINE JFNK_TEMP_COL_AUDIT_V1


!
! JFNK_JVP_V1: Jacobian-vector product by forward finite differences.
!
! Computes the TRUE J*v using a single global step size:
!   h = EPS0 * max(1, ||x||_inf) / ||v||_inf
!   x_pert = x + h * v
!   J*v = (F(x_pert) - F(x)) / h
!
! This is the correct finite-difference approximation of the directional
! derivative J*v.  It is consistent with:
!   - the GMRES RHS (jfnk_fx = raw physical F(x))
!   - the preconditioner (JFNK_PRECOND_V1 = plain M^{-1} in physical space)
!   - the inner product (JFNK_DOT_V1 = RSCALE-weighted, applied uniformly)
!
! Previous bug (fixed 2026-09-09): the JVP used a per-element scaling
! SI = max(|x_i|, 1), producing diag(1/SI)*J*diag(SI)*v — a similarity-
! transformed operator that is NOT J*v.  This was inconsistent with the
! untransformed RHS and preconditioner, so GMRES solved the wrong linear
! system.  The system was correct only on Newton step 0 (jfnk_x=0, SI=1),
! then degraded as different species accumulated different correction
! sizes — matching the reported "works on iteration 0/1, then degrades"
! symptom.
!
! The trace-species perturbation-noise concern that motivated the SI
! transform is real but secondary: a global step may produce larger
! finite-difference error for species with small |x_i| and large |v_i|.
! This is handled by:
!   1. The RSCALE weighting in JFNK_DOT_V1, which equalizes species
!      contributions to the GMRES norm.
!   2. The diagonal BA preconditioner, which equalizes species scales
!      in the Krylov solve.
!   3. The finite-difference error is O(h * ||J''|| * ||v||^2), bounded
!      by O(EPS0 * max(||x||,1) * ||J''||) — small for well-scaled
!      problems and always a perturbation of the CORRECT operator,
!      not a completely wrong operator.
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_JVP_V1(N, V, JV)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_x, jfnk_fx, jfnk_jvp_failed
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_MAX_EVALS
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(IN) :: V(N)
	DOUBLE PRECISION, INTENT(OUT) :: JV(N)
!
	DOUBLE PRECISION EPS0, VNORM, XNORM, H_STEP
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: JV_ANAL(:)
	DOUBLE PRECISION FD_NORM, AN_NORM, DIFF_NORM
	LOGICAL, SAVE :: JVP_DIAG_DONE = .FALSE.
	INTEGER I, KK, JJ, II, NRETRY
!
! Reset per-call NaN counter and failure flag (audit 17-Sep-2026).
!
	JFNK_N_JVP_NAN = 0
	jfnk_jvp_failed = .FALSE.
!
	JFNK_W_V_USE = V
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, JFNK_W_V_USE)
!
	VNORM = 0.0D0
	XNORM = 0.0D0
	DO I = 1, N
	  VNORM = MAX(VNORM, ABS(JFNK_W_V_USE(I)))
	  XNORM = MAX(XNORM, ABS(jfnk_x(I)))
	END DO
	IF(VNORM .LT. 1.0D-30)THEN
	  JV = 0.0D0
	  RETURN
	END IF
!
! Finite-difference step size (global, not per-element):
!
! Optimal forward-difference step (2026-09-11): h ~ sqrt(eps_mach)
! balances truncation error O(h) against cancellation error O(eps/h).
! 03-Oct-2026: JFNK_EPS_JVP_LOC (VADAT [JFNK_EPS]) is honored again as
! a floor in COUPLED mode only.  The coupled residual is produced by
! an iterative formal solve and carries an evaluation-noise floor
! ~delta_F, so the FD step must scale as h ~ sqrt(delta_F/|F|), not
! sqrt(eps_mach) (Kelley; Pernice & Walker 2008 on approximate
! function values).  JFNK_EPS~1e-3..1e-2 is the intended setting for
! coupled evals on the iterative Knop sweep.  The frozen residual is
! deterministic arithmetic on stored arrays, so it keeps the
! machine-eps rule; absent keys also keep machine-eps everywhere.
!
	EPS0 = SQRT(EPSILON(1.0D0))
	IF(JFNK_COUPLED_LOC)EPS0 = MAX(EPS0, JFNK_EPS_JVP_LOC)
!
! Global forward-difference JVP:
!   h = EPS0 * max(1, ||x||_inf) / ||v||_inf
!   J*v ≈ (F(x + h*v) - F(x)) / h
!
! This is the true directional derivative of F at x in direction v.
! No per-element SI scaling — the matvec, preconditioner, and RHS all
! operate in the same physical space.  The RSCALE weighting in
! JFNK_DOT_V1 provides the species-scale equalization for the GMRES
! inner product, and the diagonal BA preconditioner provides it for the
! Krylov solve.  Baking a second, inconsistent scaling into the matvec
! (as the previous SI-transformed code did) breaks the algebraic
! consistency of the GMRES system.
!
! editted by Sean Brennan (09-Sep-2026)
!
	H_STEP = EPS0 * MAX(XNORM, 1.0D0) / VNORM
!
! Safeguard (revised 17-Sep-2026, audit issue 5): keep ONLY the upper
! cap.  The previous MAX(H_STEP,1D-12) floor was not harmless: for a
! huge ||v||_inf it inflated ||h*v|| to a percent-scale perturbation,
! far above the FD scale and deep into the nonlinear clipping regime
! of the residual.  With no floor, ||h*v||_inf = EPS0*max(1,||x||) by
! construction (the cap only binds for small ||v||, making ||h*v||
! even smaller), so the perturbation is always at the intended scale.
! Instead of a floor, require the resulting perturbation to be usable:
! positive, finite, and no larger than the state norm itself.  An
! underflowed or oversized step means the direction is unusable --
! flag JVP failure rather than silently taking a bad derivative.
!
! Cap the actual perturbation norm, not scalar h: h is inversely
! proportional to ||v||. The coupled EPS0 floor is honored for every
! direction norm; frozen mode keeps its machine-eps step.
	H_STEP = MIN(H_STEP, MAX(1.0D-4, EPS0) *
	1       MAX(XNORM,1.0D0) / VNORM)
!
! COLPROBE h-scan hook (see JFNK_CMFGEN_INT_MOD): a diagnostic-only
! multiplier applied after the safety cap so probes can step across
! the FD cancellation floor.  The 0/Inf and ||h*v|| checks below
! still apply to the scaled step.
!
	IF(JFNK_PROBE_HSCALE .NE. 1.0_LDP)
	1   H_STEP = H_STEP * JFNK_PROBE_HSCALE
	IF(.NOT. (H_STEP .GT. 0.0D0) .OR.
	1     H_STEP * VNORM .GT. MAX(XNORM, 1.0D0))THEN
	  JV = 0.0D0
	  jfnk_jvp_failed = .TRUE.
	  JFNK_N_JVP_FAILURES = JFNK_N_JVP_FAILURES + 1
	  IF(JFNK_LUER .GT. 0) WRITE(JFNK_LUER,'(A,2ES12.4)')
	1    'JFNK_JVP_V1: unusable FD step h,||hv||=', H_STEP,
	2    H_STEP * VNORM
	  RETURN
	END IF
!
! Audit 17-Sep-2026 (issue 4): the residual is nonsmooth at the
! population floor and the T-exp cap.  If the perturbation h*v crossed
! a constraint that was inactive at jfnk_x (counted by the residual in
! JFNK_N_NEW_CLIP), the FD difference contains a clipped, artificial
! derivative.  Retry with h reduced by 10x -- up to two retries; if a
! crossing persists the derivative is still usable at the smaller
! step, so the last attempt is kept (GMRES tolerates inexactness).
!
	DO NRETRY = 0, 2
	  DO I = 1, N
	    JFNK_W_X_PERT(I) = jfnk_x(I) + H_STEP * JFNK_W_V_USE(I)
	  END DO
	  CALL JFNK_RESIDUAL_V1(N, JFNK_W_X_PERT, JFNK_W_F_PERT)
	  IF(JFNK_MAX_EVALS .GT. 0 .AND.
	1       JFNK_N_RESIDUAL_EVALS .GT. JFNK_MAX_EVALS)THEN
! Budget sentinel eval returned -- flag JVP failure so GMRES aborts
! immediately rather than building a Krylov vector on the sentinel.
	    JV = 0.0D0
	    jfnk_jvp_failed = .TRUE.
	    JFNK_N_JVP_FAILURES = JFNK_N_JVP_FAILURES + 1
	    RETURN
	  END IF
	  IF(JFNK_N_NEW_CLIP .EQ. 0)EXIT
	  JFNK_N_JVP_HRETRY = JFNK_N_JVP_HRETRY + 1
!
! edited by Devin, 05-Oct-2026: do not shrink the perturbation
! norm h*||v|| below EPS0*max(1,||x||) in coupled mode.  The
! noise-aware JFNK_EPS_JVP_LOC sets that scale from the residual
! evaluation noise (Pernice & Walker; see the EPS0 block above).
! Comparing h alone is wrong because h is inversely proportional to
! ||v||.  Frozen residuals keep deterministic arithmetic and may
! continue shrinking on a clip retry.
!
	  IF(JFNK_COUPLED_LOC .AND.
	1       H_STEP*0.1D0*VNORM .LT.
	2       EPS0*MAX(XNORM,1.0D0))THEN
	    IF(JFNK_LUER .GT. 0) WRITE(JFNK_LUER,'(A,I8,A)')
	1      'JFNK_JVP_V1: ', JFNK_N_NEW_CLIP,
	2      ' crossings; perturbation at EPS0 floor, keeping JV'
	    EXIT
	  END IF
	  IF(JFNK_LUER .GT. 0) WRITE(JFNK_LUER,'(A,I8,A,ES10.3)')
	1    'JFNK_JVP_V1: ', JFNK_N_NEW_CLIP,
	2    ' new constraint crossings; reducing h to', H_STEP*0.1D0
	  H_STEP = H_STEP * 0.1D0
	END DO
!
	DO I = 1, N
	  JV(I) = (JFNK_W_F_PERT(I) - jfnk_fx(I)) / H_STEP
	END DO
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, JV)
!
! One-time diagnostic: Compare FD JVP with the UNTRANSFORMED analytical
! BA_full * v.  The previous self-check applied the same SI transform
! to the analytical side, so it "passed" even when the matvec was the
! wrong operator.  Now the self-check compares FD J*v against the true
! BA_full * v, so a large discrepancy indicates genuine finite-difference
! error (or a bug in the residual evaluation), not a hidden transform.
!
	IF(.NOT. JVP_DIAG_DONE .AND. JFNK_LU_READY)THEN
	  JVP_DIAG_DONE = .TRUE.
	  IF(.NOT. ALLOCATED(JV_ANAL)) ALLOCATE(JV_ANAL(N))
	  JV_ANAL(1:N) = 0.0D0
! Analytical: JV_ANAL = BA_full * v (physical Jacobian, no SI transform,
! no ridge).  The ridge (JFNK_BA_RIDGE_ADD) is an artificial
! regularization added only to the preconditioner diagonal, not to the
! physical Jacobian.  Subtracting it from the analytical comparison
! would make a correct FD J*v appear to disagree with the analytical
! side (review §7).
	  DO KK = 1, JFNK_ND
	    DO II = 1, JFNK_NT
	      DO JJ = 1, JFNK_NT
	        JV_ANAL((KK-1)*JFNK_NT + II) =
     1        JV_ANAL((KK-1)*JFNK_NT + II) +
     2        DBLE(JFNK_BA_A(II, JJ, KK)) *
     3        JFNK_W_V_USE((KK-1)*JFNK_NT + JJ)
	      END DO
	    END DO
	  END DO
	  CALL JFNK_MASK_FIXED_UNKNOWNS(N, JV_ANAL)
!
! JFNK_BA_A is the physical BA (no ridge).  Local ridge is applied
! only to the LU copy used as a preconditioner.
!
	  FD_NORM = 0.0D0
	  AN_NORM = 0.0D0
	  DIFF_NORM = 0.0D0
	  DO I = 1, N
	    FD_NORM = MAX(FD_NORM, ABS(JV(I)))
	    AN_NORM = MAX(AN_NORM, ABS(JV_ANAL(I)))
	    DIFF_NORM = MAX(DIFF_NORM, ABS(JV(I) - JV_ANAL(I)))
	  END DO
	  IF(JFNK_LUER .GT. 0) WRITE(JFNK_LUER,'(A)')
     1    '  JVP DIAGNOSTIC (FD vs untransformed BA_full*v):'
	  IF(JFNK_LUER .GT. 0) WRITE(JFNK_LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
     1    '    |FD|=', FD_NORM, ' |analytical|=', AN_NORM,
     1    ' |diff|/|anal|=', DIFF_NORM/MAX(AN_NORM,1.0D-30)
	  IF(JFNK_LUER .GT. 0 .AND.
     1    DIFF_NORM/MAX(AN_NORM,1.0D-30) .GT. 0.1D0)THEN
	    DO I = 1, MIN(10, N)
	      IF(JFNK_LUER .GT. 0) WRITE(JFNK_LUER,'(A,I6,A,ES14.6,A,ES14.6)')
     1        '    i=', I, ' FD=', JV(I), ' anal=', JV_ANAL(I)
	    END DO
	  END IF
	END IF
!
! NaN/Inf guard for the JVP result (revised 2026-09-07):
!
! The original code silently replaced each NaN/Inf element with 0.0,
! which masked real numerical failures.  A zero JVP element tells GMRES
! that the Jacobian has zero sensitivity in that direction, causing
! premature convergence or stagnation.
!
! Revised approach:
!   1. Count the number of NaN/Inf entries in the JVP.
!   2. If the fraction of bad entries exceeds JFNK_JVP_NAN_FRAC (default
!      10%), the entire JVP is set to zero and a failure counter is
!      incremented.  GMRES will see a zero JVP and the Newton step will
!      be rejected by the trust-region logic (pred <= 0), causing the
!      trust region to shrink.
!   3. If only a few entries are bad (below the threshold), they are
!      individually zeroed — this handles the common case where one or
!      two species have numerical issues at extreme depths.
!   4. A diagnostic message is printed when the failure threshold is
!      exceeded, so the user can identify the root cause.
!
	DO I = 1, N
	  IF(JV(I) .NE. JV(I)) JFNK_N_JVP_NAN = JFNK_N_JVP_NAN + 1
	  IF(ABS(JV(I)) .GT. 1.0D30) JFNK_N_JVP_NAN = JFNK_N_JVP_NAN + 1
	END DO
!
	IF(JFNK_N_JVP_NAN .GT. 0)THEN
	  IF(JFNK_N_JVP_NAN .GT. (N / 10))THEN
! Too many NaN entries — the JVP is unreliable.  Audit 17-Sep-2026:
! do NOT silently hand GMRES a zero J*v (GMRES interprets that as a
! genuine zero-action direction, producing artificial happy
! breakdown/stagnation).  Raise the explicit failure flag; GMRES
! aborts the Krylov solve (info=3) and Newton rejects the direction,
! shrinks the trust region, and requests a preconditioner refresh.
! JV is still zeroed as a last-resort value for any caller that
! ignores the flag.
	    JV = 0.0D0
	    jfnk_jvp_failed = .TRUE.
	    JFNK_N_JVP_FAILURES = JFNK_N_JVP_FAILURES + 1
	    IF(JFNK_LUER .GT. 0) WRITE(JFNK_LUER,'(A,I6,A,I6,A)')
	1      'JFNK_JVP_V1: ', JFNK_N_JVP_NAN, ' NaN/Inf in JVP (N=',
	1      N, ') — JVP FAILED, GMRES will abort this solve'
	  ELSE
! Few bad entries — zero them individually.  This is safe because the
! remaining entries still provide useful directional information to GMRES.
	    DO I = 1, N
	      IF(JV(I) .NE. JV(I)) JV(I) = 0.0D0
	      IF(ABS(JV(I)) .GT. 1.0D30) JV(I) = 0.0D0
	    END DO
	  END IF
	END IF
!
	JFNK_N_JVP_EVALS = JFNK_N_JVP_EVALS + 1
!
	END SUBROUTINE JFNK_JVP_V1



!
! JFNK_FACTOR_PRECOND_V1: Extract the BA diagonal for the JFNK
! diagonal (Jacobi) preconditioner.
!
! The full BA LU solve gives a Newton step s = -BA^{-1}*F that is
! numerically tiny (~1e-18) because BA spans ~13 orders of magnitude
! (H~1e10 vs Fe~1e-3).  The diagonal Jacobi preconditioner gives a
! useful Lambda-like step s = -F/diag(BA) (O(0.1-1)) that GMRES then
! refines using the full Jacobian (via FD JVP).
!
! This routine only needs the BA diagonal and the full BA matrix
! (for the JVP self-check diagnostic).  The previous DGEEQU/DGETRF
! LU machinery was dead weight — JFNK_PRECOND_V1 never called DGETRS.
!
! editted by Sean Brennan (09-Sep-2026)
!
	SUBROUTINE JFNK_FACTOR_PRECOND_V1()
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_BA_RIDGE_FAC, JFNK_PRECOND_TYPE,
	1                              JFNK_ROW_TINY, JFNK_LSQ_LAM,
	2                              JFNK_RIDGE_PIV, JFNK_POP_FLOOR
	USE MOD_CMFGEN, ONLY: ATM, NUM_IONS
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_ridge_floor, jfnk_ridge_floor_set,
	1     jfnk_ptc_active, jfnk_fx
	IMPLICIT NONE
!
	REAL(KIND=LDP) DUMMY_STEQ(JFNK_NT)
	REAL(KIND=LDP) MAX_VAL, RIDGE_ADD
	REAL(KIND=LDP) ROW_CND, COL_CND, AMAX
	REAL(KIND=LDP) EPS, SCALE_K
	REAL(KIND=LDP) ROW_SF_MIN, ROW_SF_MAX, COL_SF_MIN, COL_SF_MAX
	LOGICAL ZERO_STEQ(JFNK_NT)
	LOGICAL REPLACE(JFNK_NION)
	LOGICAL FIRST_MATRIX, LAST_MATRIX, USE_PASSED_REP
	LOGICAL ROW_REG(JFNK_NT)
	INTEGER K, J, I, IOS, LUER, ERROR_LU
	INTEGER JFNK_LU_INFO
	INTEGER ITRY, N_OK, N_FAIL, N_RIDGE, LU_INFO_K
	INTEGER IB, L, II, JJ, BANDED_INFO, BAND_OFFSET
	INTEGER IROW, JCOL
	INTEGER N_REG_ROWS
	INTEGER II_REG, N_REG_TOP, K_REG_TOP
	INTEGER NPB, ID_ION, N_POP_FRZ
	REAL(KIND=LDP) ION_TOT
	INTEGER, ALLOCATABLE :: NREG_DEPTH(:)
	REAL(KIND=LDP) ROWMAX_I, BLKMAX
	REAL(KIND=LDP) BANDED_RIDGE, BANDED_MAX_DIAG
	REAL(KIND=LDP) LSQ_L2
	REAL(KIND=LDP), ALLOCATABLE :: GN_WORK(:,:)
	LOGICAL, ALLOCATABLE :: REPLACE_SAVE(:,:)
	INTEGER N_BANDED_RETRY
	DOUBLE PRECISION FMAX_PTC, S_PTC
	INTEGER IDX_PTC
	DOUBLE PRECISION ANORM_K, RCOND_K, PIV_MIN, PIV_MAX, S_TMP
	DOUBLE PRECISION, ALLOCATABLE :: GECON_WORK(:)
	INTEGER, ALLOCATABLE :: GECON_IWORK(:)
	INTEGER LU_INFO_G
	EXTERNAL ERROR_LU
!
	LUER = ERROR_LU()
	JFNK_LU_READY = .FALSE.
	JFNK_FACTOR_DONE = .FALSE.
!
! Capture the user-supplied ridge as the relaxation floor on the first
! factorization, so the relaxation path cannot push the ridge below the
! value the user chose.  Added 2026-09-10.
!
	IF(.NOT. jfnk_ridge_floor_set)THEN
	  jfnk_ridge_floor = MAX(JFNK_BA_RIDGE_FAC, 0.0_LDP)
	  jfnk_ridge_floor_set = .TRUE.
	END IF
!
	IOS=0
	IF(ALLOCATED(JFNK_BA_A))THEN
	  IF(SIZE(JFNK_BA_A,1).NE.JFNK_NT .OR.
	1    SIZE(JFNK_BA_A,2).NE.JFNK_NT .OR.
	2    SIZE(JFNK_BA_A,3).NE.JFNK_ND)THEN
	    DEALLOCATE(JFNK_BA_A)
	  END IF
	END IF
	IF(.NOT. ALLOCATED(JFNK_BA_A))THEN
	  ALLOCATE(JFNK_BA_A(JFNK_NT,JFNK_NT,JFNK_ND), STAT=IOS)
	  IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc BA_A failed'
	END IF
!
	IF(ALLOCATED(JFNK_BA_DIAG))THEN
	  IF(SIZE(JFNK_BA_DIAG,1).NE.JFNK_NT .OR.
	1    SIZE(JFNK_BA_DIAG,2).NE.JFNK_ND) DEALLOCATE(JFNK_BA_DIAG)
	END IF
	IF(.NOT. ALLOCATED(JFNK_BA_DIAG))THEN
	  ALLOCATE(JFNK_BA_DIAG(JFNK_NT,JFNK_ND), STAT=IOS)
	  IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc BA_DIAG failed'
	END IF
	IF(ALLOCATED(JFNK_PRECOND_DIAG))THEN
	  IF(SIZE(JFNK_PRECOND_DIAG,1).NE.JFNK_NT .OR.
	1     SIZE(JFNK_PRECOND_DIAG,2).NE.JFNK_ND)
	2    DEALLOCATE(JFNK_PRECOND_DIAG)
	END IF
	IF(.NOT. ALLOCATED(JFNK_PRECOND_DIAG))THEN
	  ALLOCATE(JFNK_PRECOND_DIAG(JFNK_NT,JFNK_ND), STAT=IOS)
	  IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc PC_DIAG failed'
	END IF
	IF(ALLOCATED(JFNK_BA_RIDGE_ADD))THEN
	  IF(SIZE(JFNK_BA_RIDGE_ADD,1).NE.JFNK_NT .OR.
	1     SIZE(JFNK_BA_RIDGE_ADD,2).NE.JFNK_ND)
	2    DEALLOCATE(JFNK_BA_RIDGE_ADD)
	END IF
	IF(.NOT. ALLOCATED(JFNK_BA_RIDGE_ADD))THEN
	  ALLOCATE(JFNK_BA_RIDGE_ADD(JFNK_NT, JFNK_ND), STAT=IOS)
	  IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc RIDGE_ADD failed'
	END IF
	JFNK_BA_RIDGE_ADD = 0.0_LDP
!
! Allocate LU factorization workspace (used when JFNK_PRECOND_TYPE='LU').
! JFNK_LU_FAC stores the equilibrated, LU-factored BA matrix per depth.
! JFNK_IPIV stores the pivot indices.  JFNK_ROW_SF / JFNK_COL_SF store
! the DGEEQU row/column scaling factors.  'LSQ' shares this workspace:
! LU_FAC then holds the Cholesky factor of J_eq^T J_eq + lam^2 I and
! IPIV is unused.
!
	IF(JFNK_PRECOND_TYPE(1:2) .EQ. 'LU' .OR.
	1   JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ' .OR.
	2   JFNK_PRECOND_TYPE(1:5) .EQ. 'SCHUR')THEN
	  IF(ALLOCATED(JFNK_LU_FAC))THEN
	    IF(SIZE(JFNK_LU_FAC,1).NE.JFNK_NT .OR.
	1       SIZE(JFNK_LU_FAC,2).NE.JFNK_NT .OR.
	2       SIZE(JFNK_LU_FAC,3).NE.JFNK_ND)THEN
	      DEALLOCATE(JFNK_LU_FAC)
	    END IF
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_LU_FAC))THEN
	    ALLOCATE(JFNK_LU_FAC(JFNK_NT,JFNK_NT,JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc LU_FAC failed'
	  END IF
	  IF(ALLOCATED(JFNK_IPIV))THEN
	    IF(SIZE(JFNK_IPIV,1).NE.JFNK_NT .OR.
	1       SIZE(JFNK_IPIV,2).NE.JFNK_ND) DEALLOCATE(JFNK_IPIV)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_IPIV))THEN
	    ALLOCATE(JFNK_IPIV(JFNK_NT,JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc IPIV failed'
	  END IF
	  IF(ALLOCATED(JFNK_ROW_SF))THEN
	    IF(SIZE(JFNK_ROW_SF,1).NE.JFNK_NT .OR.
	1       SIZE(JFNK_ROW_SF,2).NE.JFNK_ND) DEALLOCATE(JFNK_ROW_SF)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_ROW_SF))THEN
	    ALLOCATE(JFNK_ROW_SF(JFNK_NT,JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc ROW_SF failed'
	  END IF
	  IF(ALLOCATED(JFNK_COL_SF))THEN
	    IF(SIZE(JFNK_COL_SF,1).NE.JFNK_NT .OR.
	1       SIZE(JFNK_COL_SF,2).NE.JFNK_ND) DEALLOCATE(JFNK_COL_SF)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_COL_SF))THEN
	    ALLOCATE(JFNK_COL_SF(JFNK_NT,JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc COL_SF failed'
	  END IF
	  IF(ALLOCATED(JFNK_LU_OK))THEN
	    IF(SIZE(JFNK_LU_OK) .NE. JFNK_ND) DEALLOCATE(JFNK_LU_OK)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_LU_OK))THEN
	    ALLOCATE(JFNK_LU_OK(JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc LU_OK failed'
	  END IF
	  IF(ALLOCATED(JFNK_ROW_REG))THEN
	    IF(SIZE(JFNK_ROW_REG,1).NE.JFNK_NT .OR.
	1       SIZE(JFNK_ROW_REG,2).NE.JFNK_ND) DEALLOCATE(JFNK_ROW_REG)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_ROW_REG))THEN
	    ALLOCATE(JFNK_ROW_REG(JFNK_NT,JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc ROW_REG failed'
	  END IF
	  JFNK_ROW_REG = .FALSE.
	  IF(ALLOCATED(JFNK_DEPTH_RIDGE))THEN
	    IF(SIZE(JFNK_DEPTH_RIDGE) .NE. JFNK_ND)
	1       DEALLOCATE(JFNK_DEPTH_RIDGE)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_DEPTH_RIDGE))THEN
	    ALLOCATE(JFNK_DEPTH_RIDGE(JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc DEPTH_RIDGE failed'
	  END IF
!
! Per-depth LU diagnostics (pivot ratio, DGECON rcond, accepted retry
! index).  Written on every successful DGETRF -- missing the alloc here
! segfaulted the first JFNK-CPL admission on the m15 KNOP probe.
!
	  IF(ALLOCATED(JFNK_LU_PIVRAT))THEN
	    IF(SIZE(JFNK_LU_PIVRAT) .NE. JFNK_ND) DEALLOCATE(JFNK_LU_PIVRAT,
	1       JFNK_LU_RCOND, JFNK_LU_ITRY)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_LU_PIVRAT))THEN
	    ALLOCATE(JFNK_LU_PIVRAT(JFNK_ND), JFNK_LU_RCOND(JFNK_ND),
	1       JFNK_LU_ITRY(JFNK_ND), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc LU diag failed'
	  END IF
	  IF(JFNK_PRECOND_TYPE(1:5) .EQ. 'SCHUR')THEN
	    IF(ALLOCATED(JFNK_SCHUR_S))THEN
	      IF(SIZE(JFNK_SCHUR_S,1).NE.JFNK_ND)THEN
	        DEALLOCATE(JFNK_SCHUR_S, JFNK_SCHUR_IPIV,
	1            JFNK_SCHUR_JNT, JFNK_SCHUR_JTN,
	2            JFNK_SCHUR_JTT)
	      END IF
	    END IF
	    IF(.NOT. ALLOCATED(JFNK_SCHUR_S))THEN
	      ALLOCATE(JFNK_SCHUR_S(JFNK_ND,JFNK_ND),
	1         JFNK_SCHUR_IPIV(JFNK_ND),
	2         JFNK_SCHUR_JNT(JFNK_NT,JFNK_ND),
	3         JFNK_SCHUR_JTN(JFNK_NT,JFNK_ND),
	4         JFNK_SCHUR_JTT(JFNK_ND), STAT=IOS)
	      IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc SCHUR failed'
	    END IF
	  END IF
	END IF
!
	JFNK_BA_DIAG = 0.0_LDP
	JFNK_BA_DIAG_READY = .FALSE.
!
! Allocate per-depth REPLACE save array for banded assembly.
! GENERATE_FULL_MATRIX_V3 determines REPLACE (which ground-state
! equations are replaced by ionization equations) on the diagonal
! pass.  The banded loop needs the same REPLACE for each depth.
!
	IF(ALLOCATED(REPLACE_SAVE))THEN
	  IF(SIZE(REPLACE_SAVE,1).NE.JFNK_NION .OR.
	1    SIZE(REPLACE_SAVE,2).NE.JFNK_ND) DEALLOCATE(REPLACE_SAVE)
	END IF
	IF(.NOT. ALLOCATED(REPLACE_SAVE))THEN
	  ALLOCATE(REPLACE_SAVE(JFNK_NION, JFNK_ND), STAT=IOS)
	  IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc REPLACE_SAVE'
	END IF
!
	ZERO_STEQ = .FALSE.
	REPLACE = .FALSE.
	USE_PASSED_REP = .FALSE.
	FIRST_MATRIX = .TRUE.
	LAST_MATRIX = .FALSE.
!
! PTC fold scale: s_i = max(|F_i|, 1e-8*|F|max), matching
! JFNK_SET_PTC_SCALE_V1.  Computed here because SET_PTC_SCALE arms
! only after this routine returns on the admission path.  dtau is
! 1 whenever this routine runs (reset on admission/rebuild before
! any Newton step), so folding s directly matches the operator the
! FD JVP differentiates.
!
	JFNK_PTC_FOLDED = .FALSE.
	FMAX_PTC = 0.0D0
	IF(jfnk_ptc_active .AND. ALLOCATED(jfnk_fx))THEN
	  IF(SIZE(jfnk_fx) .EQ. JFNK_NT*JFNK_ND)THEN
	    DO IDX_PTC = 1, JFNK_NT*JFNK_ND
	      FMAX_PTC = MAX(FMAX_PTC, ABS(jfnk_fx(IDX_PTC)))
	    END DO
	    IF(FMAX_PTC .GT. 0.0D0)JFNK_PTC_FOLDED = .TRUE.
	  END IF
	END IF
!
	DO K = 1, JFNK_ND
	  IF(K .EQ. JFNK_ND) LAST_MATRIX = .TRUE.
	  CALL GENERATE_FULL_MATRIX_V3(
	1       JFNK_C_MAT, DUMMY_STEQ, POPS_WORK, REPLACE, ZERO_STEQ,
	1       JFNK_NT, JFNK_ND, JFNK_NION, JFNK_NUM_BNDS,
	1       JFNK_DIAG_INDX, JFNK_DIAG_INDX, K,
	1       FIRST_MATRIX, LAST_MATRIX, USE_PASSED_REP)
	  FIRST_MATRIX = .FALSE.
	  REPLACE_SAVE(:, K) = REPLACE
!
	  DO J = 1, JFNK_NT
	    JFNK_BA_DIAG(J, K) = JFNK_C_MAT(J, J)
	    JFNK_PRECOND_DIAG(J, K) = JFNK_C_MAT(J, J)
	  END DO
!
! Store the physical BA block with no global ridge.  Regularization
! is applied per depth only if equilibrated LU fails.
! editted by Sean Brennan (10-Sep-2026)
!
	  DO I = 1, JFNK_NT
	    JFNK_BA_RIDGE_ADD(I, K) = 0.0_LDP
	  END DO
	  JFNK_BA_A(:,:,K) = JFNK_C_MAT
!
! Fold s_i into the factored diagonal.  Bare BA is nearly singular
! on thousands of weak rows; the LU then emits z ~ 1e2..1e5 along
! directions where the FD operator holds the s floor, and the
! post-solve D/(D+sig) rescale cannot repair an exploded LU
! direction.  With s inside the matrix every block is diagonally
! bounded, the ROW_TINY regularization mostly empties, and the
! preconditioned step matches the FD operator.  If dtau later
! grows, the shift retained here overestimates the live one --
! over-damped on weak rows but bounded, which is the safe side.
!
	  IF(JFNK_PTC_FOLDED)THEN
	    DO J = 1, JFNK_NT
	      IDX_PTC = (K-1)*JFNK_NT + J
	      S_PTC = MAX(ABS(jfnk_fx(IDX_PTC)), 1.0D-8*FMAX_PTC)
	      JFNK_BA_A(J,J,K) = JFNK_BA_A(J,J,K) + S_PTC
	      JFNK_BA_DIAG(J,K) = JFNK_BA_A(J,J,K)
	      JFNK_PRECOND_DIAG(J,K) = JFNK_BA_A(J,J,K)
	    END DO
	  END IF
!
! Frozen-row reduction.  The apply solves the full NT block then
! zeroes z on frozen components (T, or Ne when FIXED_NE), which
! discards the column terms the LU solve already accounted for and
! leaves ~16% residual on population rows.  Solve the reduced system
! instead: zero the frozen row/column off-diagonals and set the
! diagonal to 1, so z_p = A_pp^{-1} v_p (s_T=0 enforced) and the
! post-solve projection is a no-op on the answer.
!
	  IF(JFNK_FIXED_T_LOC .OR.
	1    (JFNK_FIX_T_D_ST .GT. 0 .AND.
	2     K .GE. JFNK_FIX_T_D_ST .AND.
	3     K .LE. JFNK_FIX_T_D_END .AND.
	4     .NOT. JFNK_TAUTO_ACTIVE))THEN
	    DO J = 1, JFNK_NT-1
	      JFNK_BA_A(J, JFNK_NT, K) = 0.0_LDP
	      JFNK_BA_A(JFNK_NT, J, K) = 0.0_LDP
	    END DO
	    JFNK_BA_A(JFNK_NT, JFNK_NT, K) = 1.0_LDP
	    JFNK_BA_DIAG(JFNK_NT,K) = 1.0_LDP
	    JFNK_PRECOND_DIAG(JFNK_NT,K) = 1.0_LDP
	  END IF
	  IF(JFNK_FIXED_NE_LOC)THEN
	    DO J = 1, JFNK_NT
	      IF(J .NE. JFNK_NT-1)THEN
	        JFNK_BA_A(J, JFNK_NT-1, K) = 0.0_LDP
	        JFNK_BA_A(JFNK_NT-1, J, K) = 0.0_LDP
	      END IF
	    END DO
	    JFNK_BA_A(JFNK_NT-1, JFNK_NT-1, K) = 1.0_LDP
	    JFNK_BA_DIAG(JFNK_NT-1,K) = 1.0_LDP
	    JFNK_PRECOND_DIAG(JFNK_NT-1,K) = 1.0_LDP
	  END IF
	END DO
!
! ── Equilibrated depth-block LU ────────────────────────────────────
! M_K = BA(:,:,K).  Attempt DGEEQU + DGETRF with zero ridge first.
! If a block is singular/unsafe, add a local ridge
!   λ_K = ε × max(|diag_K|, 1)
! and increase ε until factorization succeeds.  Jacobi fallback is
! per-depth, not a global abort.
!
	N_OK = 0
	N_FAIL = 0
	N_RIDGE = 0
	N_REG_ROWS = 0
	N_POP_FRZ = 0
	ROW_SF_MIN = 1.0E99_LDP
	ROW_SF_MAX = 0.0_LDP
	COL_SF_MIN = 1.0E99_LDP
	COL_SF_MAX = 0.0_LDP
	JFNK_N_RIDGE_DEPTHS = 0
	JFNK_N_LU_FAIL_DEPTH = 0
!
	IF(JFNK_PRECOND_TYPE(1:2) .EQ. 'LU' .OR.
	1   JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ')THEN
	  JFNK_LU_INFO = 0
	  IF(JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ')THEN
	    ALLOCATE(GN_WORK(JFNK_NT,JFNK_NT), STAT=IOS)
	    IF(IOS .NE. 0) STOP 'JFNK_FACTOR_PRECOND: alloc GN_WORK'
	  END IF
	  DO K = 1, JFNK_ND
	    JFNK_LU_OK(K) = .FALSE.
	    JFNK_DEPTH_RIDGE(K) = 0.0_LDP
	    MAX_VAL = 0.0_LDP
	    DO I = 1, JFNK_NT
	      MAX_VAL = MAX(MAX_VAL, ABS(JFNK_BA_A(I, I, K)))
	    END DO
	    SCALE_K = MAX(MAX_VAL, 1.0_LDP)
!
! Row regularization (JFNK_ROW_TINY): compute blockmax once per depth
! from the raw block.  Rows whose maximum |A_ij| is below
! JFNK_ROW_TINY*blockmax carry no physical information and wreck the
! equilibration (row scale factors spanning 1e-12..1e42 have been seen
! in production).  They are replaced by the identity row below, before
! DGEEQU/DGETRF.  Only the preconditioner matrix is modified — never
! the residual or the JVP.
! edited by Sean Brennan (2026-10)
!
	    BLKMAX = 0.0_LDP
	    DO I = 1, JFNK_NT
	      DO J = 1, JFNK_NT
	        BLKMAX = MAX(BLKMAX, ABS(JFNK_BA_A(I, J, K)))
	      END DO
	    END DO
!
	    DO ITRY = 0, 8
	      ROW_REG=.FALSE.
	      IF(ITRY .EQ. 0)THEN
! editted by Devin, 08-Oct-2026: JFNK_BA_RIDGE_FAC (via
! jfnk_ridge_floor) is the user-set LM base ridge -- the
! floor was captured at first factorization but never
! consumed, so every depth always started at EPS=0.
! Ridge stays 0 unless the user asks for it.
	        EPS = jfnk_ridge_floor
	      ELSE IF(ITRY .EQ. 1)THEN
	        EPS = 1.0E-16_LDP
	      ELSE IF(ITRY .EQ. 2)THEN
	        EPS = 1.0E-14_LDP
	      ELSE IF(ITRY .EQ. 3)THEN
	        EPS = 1.0E-12_LDP
	      ELSE IF(ITRY .EQ. 4)THEN
	        EPS = 1.0E-10_LDP
	      ELSE IF(ITRY .EQ. 5)THEN
	        EPS = 1.0E-8_LDP
	      ELSE IF(ITRY .EQ. 6)THEN
	        EPS = 1.0E-6_LDP
	      ELSE IF(ITRY .EQ. 7)THEN
	        EPS = 1.0E-4_LDP
	      ELSE
	        EPS = 1.0E-2_LDP
	      END IF
! editted by Devin, 08-Oct-2026: ridge is per-row proportional
! (LM form J + alpha*diag(J)), not uniform lambda*I -- a uniform
! raw-unit ridge dominates the small trace-ion diagonals after
! equilibration and erases exactly the dynamics that need Newton
! updates.  RIDGE_ADD is now the dimensionless factor.
!
	      RIDGE_ADD = MAX(EPS, jfnk_ridge_floor)
	      JFNK_LU_FAC(:,:,K) = JFNK_BA_A(:,:,K)
	      DO I = 1, JFNK_NT
	        JFNK_LU_FAC(I, I, K) = JFNK_LU_FAC(I, I, K) +
	1           RIDGE_ADD *
	2           MAX(ABS(JFNK_BA_A(I, I, K)), 1.0E-30_LDP*SCALE_K)
	      END DO
!
! Replace near-null rows by the identity row before equilibration.
! Recomputed each retry (identical result every time; counted once).
! The rowmax test reads the un-augmented BA_A row, not LU_FAC: the
! ridge added above must not let a physically-null row escape
! regularization on retry ITRY>0 and invert to an oversized z.
!
	      IF(BLKMAX .GT. 0.0_LDP)THEN
	        DO I = 1, JFNK_NT
	          ROWMAX_I = 0.0_LDP
	          DO J = 1, JFNK_NT
	            ROWMAX_I = MAX(ROWMAX_I, ABS(JFNK_BA_A(I, J, K)))
	          END DO
	          IF(ROWMAX_I .LT. JFNK_ROW_TINY*BLKMAX)THEN
	            ROW_REG(I)=.TRUE.
	            DO J = 1, JFNK_NT
	              JFNK_LU_FAC(I, J, K) = 0.0_LDP
	            END DO
	            JFNK_LU_FAC(I, I, K) = 1.0_LDP
	            IF(ITRY .EQ. 0)N_REG_ROWS = N_REG_ROWS + 1
	          END IF
	        END DO
	      END IF
!
! [JFNK_POP_FL] (default 0 = off): freeze the Newton correction on
! near-empty levels.  Rows with n_i < JFNK_POP_FL * (ion total at this
! depth) dominated the top-|z| list on the m15 probe (high-lying Fe/Ar/
! Ca/Ne/O levels, rad/F ~ 1): their BA rows are near-null, the LU solve
! returns |z|~hundreds of dex, and the commit clip flattens the answer
! anyway.  Marking them regularized gives an identity row in the factor
! and z_i = 0 at apply, so Newton leaves them on the outer iteration.
! The residual row itself is unchanged -- this only bounds the step.
! Recomputed each retry because ROW_REG resets with the matrix.
!
	      IF(JFNK_POP_FLOOR .GT. 0.0_LDP .AND.
	1       ALLOCATED(POPS_BASE))THEN
	        DO ID_ION = 1, NUM_IONS
	          IF(.NOT. ATM(ID_ION)%XzV_PRES)CYCLE
	          ION_TOT = 0.0_LDP
	          DO J = ATM(ID_ION)%EQXzV,
	1             MIN(ATM(ID_ION)%EQXzV+ATM(ID_ION)%NXzV-1, JFNK_NT-2)
	            IF(J .GE. 1)ION_TOT = ION_TOT + POPS_BASE(J,K)
	          END DO
	          IF(ION_TOT .LE. 0.0_LDP)CYCLE
	          DO J = ATM(ID_ION)%EQXzV,
	1             MIN(ATM(ID_ION)%EQXzV+ATM(ID_ION)%NXzV-1, JFNK_NT-2)
	            IF(J .GE. 1 .AND. .NOT. ROW_REG(J) .AND.
	1             POPS_BASE(J,K) .LT. JFNK_POP_FLOOR*ION_TOT)THEN
	              ROW_REG(J)=.TRUE.
	              DO JJ = 1, JFNK_NT
	                JFNK_LU_FAC(J, JJ, K) = 0.0_LDP
	              END DO
	              JFNK_LU_FAC(J, J, K) = 1.0_LDP
	              IF(ITRY .EQ. 0)THEN
	                N_REG_ROWS = N_REG_ROWS + 1
	                N_POP_FRZ = N_POP_FRZ + 1
	              END IF
	            END IF
	          END DO
	        END DO
	      END IF
!
	      CALL DGEEQU(JFNK_NT, JFNK_NT, JFNK_LU_FAC(1,1,K),
	1         JFNK_NT, JFNK_ROW_SF(1,K), JFNK_COL_SF(1,K),
	2         ROW_CND, COL_CND, AMAX, LU_INFO_K)
	      IF(LU_INFO_K .NE. 0) CYCLE
	      DO J = 1, JFNK_NT
	        DO I = 1, JFNK_NT
	          JFNK_LU_FAC(I,J,K) = JFNK_LU_FAC(I,J,K)
	1           * JFNK_ROW_SF(I,K) * JFNK_COL_SF(J,K)
	        END DO
	      END DO
	      IF(JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ')THEN
!
! LSQ factorization: solve min|J_eq y - v'|^2 + lam^2|y|^2 per block.
! J_eq is the same ridged, reg-substituted, equilibrated block the LU
! path factors.  The LU solve amplifies the near-null singular
! directions of the population Jacobian into corrections of 1e24+ in
! log space that GMRES cannot use (cmfgen_jfnk_unit_test T4/T5); the
! Tikhonov normal equations bound |y| <= |v'|/(2 lam) instead, so the
! preconditioned operator stays finite and Krylov-consistent.
! lam = JFNK_LSQ_LAM * max|J_eq|; max|J_eq| ~ 1 after DGEEQU.
!
        LSQ_L2 = (JFNK_LSQ_LAM * MAX(AMAX, 1.0_LDP))**2
        CALL DSYRK('L', 'T', JFNK_NT, JFNK_NT, 1.0_LDP,
	1         JFNK_LU_FAC(1,1,K), JFNK_NT, 0.0_LDP,
	2         GN_WORK, JFNK_NT)
        DO I = 1, JFNK_NT
          GN_WORK(I,I) = GN_WORK(I,I) + LSQ_L2
        END DO
        JFNK_LU_FAC(:,:,K) = GN_WORK
        CALL DPOTRF('L', JFNK_NT, JFNK_LU_FAC(1,1,K),
	1         JFNK_NT, LU_INFO_K)
      ELSE
        ANORM_K = 0.0_LDP
        DO J = 1, JFNK_NT
          S_TMP = 0.0_LDP
          DO I = 1, JFNK_NT
            S_TMP = S_TMP + ABS(JFNK_LU_FAC(I,J,K))
          END DO
          ANORM_K = MAX(ANORM_K, S_TMP)
        END DO
        CALL DGETRF(JFNK_NT, JFNK_NT, JFNK_LU_FAC(1,1,K),
	1         JFNK_NT, JFNK_IPIV(1,K), LU_INFO_K)
      END IF
	      IF(LU_INFO_K .EQ. 0)THEN
	        JFNK_LU_OK(K) = .TRUE.
	        JFNK_DEPTH_RIDGE(K) = RIDGE_ADD
	        JFNK_LU_ITRY(K) = ITRY
	        PIV_MIN = HUGE(1.0_LDP)
	        PIV_MAX = 0.0_LDP
	        DO I = 1, JFNK_NT
	          PIV_MIN = MIN(PIV_MIN, ABS(JFNK_LU_FAC(I,I,K)))
	          PIV_MAX = MAX(PIV_MAX, ABS(JFNK_LU_FAC(I,I,K)))
	        END DO
	        JFNK_LU_PIVRAT(K) = PIV_MIN/MAX(PIV_MAX,1.0E-300_LDP)
	        JFNK_LU_RCOND(K) = -1.0_LDP
	        IF(JFNK_PRECOND_TYPE(1:2) .EQ. 'LU')THEN
	          IF(.NOT. ALLOCATED(GECON_WORK))THEN
	            ALLOCATE(GECON_WORK(4*JFNK_NT),
	1               GECON_IWORK(JFNK_NT), STAT=IOS)
	            IF(IOS .NE. 0) STOP 'alloc GECON_WORK failed'
	          END IF
	          CALL DGECON('I', JFNK_NT, JFNK_LU_FAC(1,1,K), JFNK_NT,
	1           ANORM_K, RCOND_K, GECON_WORK, GECON_IWORK, LU_INFO_G)
	          IF(LU_INFO_G .EQ. 0)JFNK_LU_RCOND(K) = RCOND_K
	        END IF
	        IF(JFNK_PRECOND_TYPE(1:2) .EQ. 'LU' .AND.
	1          JFNK_RIDGE_PIV .GT. 0.0_LDP .AND.
	2          JFNK_LU_PIVRAT(K) .LT. JFNK_RIDGE_PIV)THEN
! DGETRF can return INFO=0 on a near-singular block whose smallest
! pivot is still large enough to avoid an exact zero; the solve then
! inherits the full amplification.  When JFNK_RIDGE_PIV>0 the first
! successful-but-ill-conditioned factorization is rejected so the
! ITRY loop escalates the ridge instead.  Default 0 keeps the
! previous accept-on-INFO=0 behaviour.
	          JFNK_LU_OK(K) = .FALSE.
	          JFNK_DEPTH_RIDGE(K) = 0.0_LDP
	          CYCLE
	        END IF
	        DO I = 1, JFNK_NT
	          JFNK_BA_RIDGE_ADD(I, K) = RIDGE_ADD *
	1           MAX(ABS(JFNK_BA_A(I, I, K)), 1.0E-30_LDP*SCALE_K)
	          JFNK_PRECOND_DIAG(I, K) =
	1             JFNK_BA_DIAG(I, K) + JFNK_BA_RIDGE_ADD(I, K)
	          IF(ROW_REG(I))JFNK_PRECOND_DIAG(I, K)=1.0_LDP
	        END DO
	        N_OK = N_OK + 1
	        IF(RIDGE_ADD .GT. 0.0_LDP) N_RIDGE = N_RIDGE + 1
	        DO I = 1, JFNK_NT
	          ROW_SF_MIN = MIN(ROW_SF_MIN, JFNK_ROW_SF(I,K))
	          ROW_SF_MAX = MAX(ROW_SF_MAX, JFNK_ROW_SF(I,K))
	          COL_SF_MIN = MIN(COL_SF_MIN, JFNK_COL_SF(I,K))
	          COL_SF_MAX = MAX(COL_SF_MAX, JFNK_COL_SF(I,K))
	        END DO
	        EXIT
	      END IF
	    END DO
!
! Record the regularized-row mask for this depth: at apply time these
! components are projected to zero rather than solved (which returned
! z_i~v_i through the identity row and tripped the |z|>1D8 clip).
!
	    DO I = 1, JFNK_NT
	      JFNK_ROW_REG(I, K) = ROW_REG(I)
	    END DO
!
	    IF(.NOT. JFNK_LU_OK(K))THEN
	      N_FAIL = N_FAIL + 1
	      IF(LUER .GT. 0) WRITE(LUER,'(A,I4,A)')
	1       '  JFNK: LU failed at depth', K, ' — Jacobi fallback'
	    END IF
	  END DO
!
	  IF(ALLOCATED(GN_WORK)) DEALLOCATE(GN_WORK)
!
	  JFNK_N_RIDGE_DEPTHS = N_RIDGE
	  JFNK_N_LU_FAIL_DEPTH = N_FAIL
	  JFNK_LU_READY = (N_OK .GT. 0)
	  JFNK_LU_INFO = 0
	  IF(N_OK .EQ. 0) JFNK_LU_INFO = 1
	  IF(LUER .GT. 0)THEN
	    WRITE(LUER,'(A,A,A,I4,A,I4,A,I4,A,I4)')
	1      '  JFNK: equilibrated block ', JFNK_PRECOND_TYPE,
	2      '  nfact=',
	1      JFNK_N_LU_FACTORS+1, '  ok=', N_OK,
	1      '  ridged=', N_RIDGE, '  Jacobi=', N_FAIL
	    IF(N_OK .GT. 0)THEN
	      WRITE(LUER,'(A,2ES10.2,A,2ES10.2)')
	1        '  JFNK: row scale min/max=', ROW_SF_MIN, ROW_SF_MAX,
	1        '  col scale min/max=', COL_SF_MIN, COL_SF_MAX
	    END IF
	    IF(N_REG_ROWS .GT. 0)THEN
	      WRITE(LUER,'(A,I5,A,ES10.2,A)')
	1        '  JFNK: regularized', N_REG_ROWS,
	1        ' near-null rows (JFNK_ROW_TINY=', JFNK_ROW_TINY, ')'
	      IF(N_POP_FRZ .GT. 0)
	1       WRITE(LUER,'(A,I5,A,ES10.2,A)')
	2        '  JFNK: of which', N_POP_FRZ,
	3        ' frozen on population floor (JFNK_POP_FL=',
	4        JFNK_POP_FLOOR, ')'
!
! editted by Devin, 08-Oct-2026: report the 5 depths carrying the
! most regularized rows.  A concentration at compression depths
! would implicate the BA_full NaN-zeroing path as the null-row
! source; a uniform spread points to structural rank deficiency
! (redundant SE / conservation equations) instead.
!
	      IF(ALLOCATED(JFNK_ROW_REG))THEN
	        ALLOCATE(NREG_DEPTH(JFNK_ND))
	        DO K = 1, JFNK_ND
	          NREG_DEPTH(K) = COUNT(JFNK_ROW_REG(:,K))
	        END DO
	        DO II_REG = 1, MIN(5, JFNK_ND)
	          K_REG_TOP = MAXLOC(NREG_DEPTH, 1)
	          N_REG_TOP = NREG_DEPTH(K_REG_TOP)
	          IF(N_REG_TOP .LE. 0)EXIT
	          WRITE(LUER,'(A,I4,A,I5,A,I4)')
	1            '    depth=', K_REG_TOP, '  reg rows=', N_REG_TOP,
	1            ' of', JFNK_NT
	          NREG_DEPTH(K_REG_TOP) = -1
	        END DO
	        DEALLOCATE(NREG_DEPTH)
	      END IF
	    END IF
	    IF(N_RIDGE .GT. 0)THEN
	      WRITE(LUER,'(A)')
	1        '  JFNK: depths requiring local ridge:'
	      DO K = 1, JFNK_ND
	        IF(JFNK_DEPTH_RIDGE(K) .GT. 0.0_LDP)THEN
	          WRITE(LUER,'(A,I4,A,ES10.2,A,ES10.2,A,ES10.2,A,I3)')
	1            '    depth=', K, '  ridge=', JFNK_DEPTH_RIDGE(K),
	2            '  min/max piv=', JFNK_LU_PIVRAT(K),
	3            '  rcond=', JFNK_LU_RCOND(K),
	4            '  itry=', JFNK_LU_ITRY(K)
	        END IF
	      END DO
	    END IF
	  END IF
	END IF
!
! ── Schur (n,T) preconditioner (JFNK_PRECOND=SCHUR) ─────────────
! Per-depth population sub-block (rows 1..NT-2) factored with the
! same ridge/reg/normal-equations chain as 'LSQ' but WITHOUT row
! equilibration: rsf weighting changes the least-squares fit and
! corrupts Jnn^{-1} JnT (unit-test t_schur).  The T
! column/row of the raw BA block is eliminated into a dense ND x ND
! Schur complement S = JTT - JTn Jnn^{-1} JnT, factored once by
! DGETRF.  S is diagonal in depth for v1 -- cross-depth T coupling is
! carried by the banded path, not duplicated here.
!
	JFNK_SCHUR_OK = .FALSE.
	IF(JFNK_PRECOND_TYPE(1:5) .EQ. 'SCHUR')THEN
	  NPB = JFNK_NT - 2
	  N_OK = 0
	  N_FAIL = 0
	  N_RIDGE = 0
	  N_REG_ROWS = 0
	  IF(ALLOCATED(GN_WORK)) DEALLOCATE(GN_WORK)
	  ALLOCATE(GN_WORK(NPB,NPB), STAT=IOS)
	  IF(IOS .NE. 0) STOP 'JFNK_FACTOR_PRECOND: alloc GN_WORK'
	  DO K = 1, JFNK_ND
	    JFNK_LU_OK(K) = .FALSE.
	    JFNK_DEPTH_RIDGE(K) = 0.0_LDP
	    JFNK_ROW_REG(JFNK_NT-1, K) = .FALSE.
	    JFNK_ROW_REG(JFNK_NT,   K) = .FALSE.
	    BLKMAX = 0.0_LDP
	    DO I = 1, NPB
	      DO J = 1, NPB
	        BLKMAX = MAX(BLKMAX, ABS(JFNK_BA_A(I, J, K)))
	      END DO
	    END DO
	    DO ITRY = 0, 8
	      ROW_REG = .FALSE.
	      IF(ITRY .EQ. 0)THEN
	        EPS = jfnk_ridge_floor
	      ELSE IF(ITRY .EQ. 1)THEN
	        EPS = 1.0E-16_LDP
	      ELSE IF(ITRY .EQ. 2)THEN
	        EPS = 1.0E-14_LDP
	      ELSE IF(ITRY .EQ. 3)THEN
	        EPS = 1.0E-12_LDP
	      ELSE IF(ITRY .EQ. 4)THEN
	        EPS = 1.0E-10_LDP
	      ELSE IF(ITRY .EQ. 5)THEN
	        EPS = 1.0E-8_LDP
	      ELSE IF(ITRY .EQ. 6)THEN
	        EPS = 1.0E-6_LDP
	      ELSE IF(ITRY .EQ. 7)THEN
	        EPS = 1.0E-4_LDP
	      ELSE
	        EPS = 1.0E-2_LDP
	      END IF
	      RIDGE_ADD = MAX(EPS, jfnk_ridge_floor)
	      JFNK_LU_FAC(:,:,K) = 0.0_LDP
	      DO J = 1, NPB
	        DO I = 1, NPB
	          JFNK_LU_FAC(I,J,K) = JFNK_BA_A(I,J,K)
	        END DO
	      END DO
	      DO I = 1, NPB
	        JFNK_LU_FAC(I, I, K) = JFNK_LU_FAC(I, I, K) +
	1           RIDGE_ADD *
	2           MAX(ABS(JFNK_BA_A(I, I, K)), 1.0E-30_LDP)
	      END DO
	      IF(BLKMAX .GT. 0.0_LDP)THEN
	        DO I = 1, NPB
	          ROWMAX_I = 0.0_LDP
	          DO J = 1, NPB
	            ROWMAX_I = MAX(ROWMAX_I, ABS(JFNK_BA_A(I, J, K)))
	          END DO
	          IF(ROWMAX_I .LT. JFNK_ROW_TINY*BLKMAX)THEN
	            ROW_REG(I)=.TRUE.
	            DO J = 1, NPB
	              JFNK_LU_FAC(I, J, K) = 0.0_LDP
	            END DO
	            JFNK_LU_FAC(I, I, K) = 1.0_LDP
	            IF(ITRY .EQ. 0)N_REG_ROWS = N_REG_ROWS + 1
	          END IF
	        END DO
	      END IF
	      DO I = 1, NPB
	        JFNK_ROW_SF(I,K) = 1.0_LDP
	        JFNK_COL_SF(I,K) = 1.0_LDP
	      END DO
	      LSQ_L2 = (JFNK_LSQ_LAM * MAX(BLKMAX, 1.0_LDP))**2
	      CALL DSYRK('L', 'T', NPB, NPB, 1.0_LDP,
	1         JFNK_LU_FAC(1,1,K), JFNK_NT, 0.0_LDP,
	2         GN_WORK, NPB)
	      DO I = 1, NPB
	        GN_WORK(I,I) = GN_WORK(I,I) + LSQ_L2
	      END DO
	      JFNK_LU_FAC(:,:,K) = 0.0_LDP
	      DO J = 1, NPB
	        DO I = 1, NPB
	          JFNK_LU_FAC(I,J,K) = GN_WORK(I,J)
	        END DO
	      END DO
	      CALL DPOTRF('L', NPB, JFNK_LU_FAC(1,1,K),
	1         JFNK_NT, LU_INFO_K)
	      IF(LU_INFO_K .EQ. 0)THEN
	        JFNK_LU_OK(K) = .TRUE.
	        JFNK_DEPTH_RIDGE(K) = RIDGE_ADD
	        DO I = 1, NPB
	          JFNK_BA_RIDGE_ADD(I, K) = RIDGE_ADD *
	1           MAX(ABS(JFNK_BA_A(I, I, K)), 1.0E-30_LDP)
	          JFNK_PRECOND_DIAG(I, K) =
	1             JFNK_BA_DIAG(I, K) + JFNK_BA_RIDGE_ADD(I, K)
	          IF(ROW_REG(I))JFNK_PRECOND_DIAG(I, K)=1.0_LDP
	        END DO
	        DO I = 1, NPB
	          JFNK_ROW_REG(I, K) = ROW_REG(I)
	        END DO
	        N_OK = N_OK + 1
	        IF(RIDGE_ADD .GT. 0.0_LDP) N_RIDGE = N_RIDGE + 1
	        EXIT
	      END IF
	    END DO
	    IF(.NOT. JFNK_LU_OK(K)) N_FAIL = N_FAIL + 1
	  END DO
!
! Extract the raw T coupling blocks and form the dense Schur
! complement.  JTN reads the raw BA row (unaugmented, unscaled);
! the same applies to JNT and JTT.
!
	  DO K = 1, JFNK_ND
	    DO J = 1, JFNK_NT
	      JFNK_SCHUR_JNT(J, K) = JFNK_BA_A(J, JFNK_NT, K)
	      JFNK_SCHUR_JTN(J, K) = JFNK_BA_A(JFNK_NT, J, K)
	    END DO
	    JFNK_SCHUR_JTT(K) = JFNK_BA_A(JFNK_NT, JFNK_NT, K)
	    JFNK_PRECOND_DIAG(JFNK_NT, K) =
	1       JFNK_SCHUR_JTT(K) + JFNK_DEPTH_RIDGE(K) *
	2       MAX(ABS(JFNK_SCHUR_JTT(K)), 1.0E-30_LDP)
	    JFNK_PRECOND_DIAG(JFNK_NT-1, K) = 1.0_LDP
	    JFNK_BA_RIDGE_ADD(JFNK_NT,   K) = JFNK_DEPTH_RIDGE(K) *
	1       MAX(ABS(JFNK_SCHUR_JTT(K)), 1.0E-30_LDP)
	    JFNK_BA_RIDGE_ADD(JFNK_NT-1, K) = 0.0_LDP
	  END DO
	  CALL JFNK_SCHUR_FORM_V1(JFNK_ND, JFNK_NT, NPB, JFNK_NT,
	1       JFNK_LU_FAC, JFNK_BA_A, JFNK_BA_RIDGE_ADD,
	2       JFNK_ROW_REG, JFNK_ROW_SF, JFNK_COL_SF,
	3       JFNK_SCHUR_S, JFNK_SCHUR_IPIV, LU_INFO_K)
	  JFNK_SCHUR_OK = (LU_INFO_K .EQ. 0) .AND. (N_OK .GT. 0)
	  JFNK_LU_READY = JFNK_SCHUR_OK
	  JFNK_N_RIDGE_DEPTHS = N_RIDGE
	  JFNK_N_LU_FAIL_DEPTH = N_FAIL
	  IF(ALLOCATED(GN_WORK)) DEALLOCATE(GN_WORK)
	  IF(LUER .GT. 0)THEN
	    WRITE(LUER,'(A,A,I4,A,I4,A,I4,A,L2)')
	1      '  JFNK: SCHUR pop-block LSQ nfact=',
	2      JFNK_N_LU_FACTORS+1, ' ok=', N_OK,
	3      ' ridged=', N_RIDGE, ' S_ok=', JFNK_SCHUR_OK
	    IF(N_REG_ROWS .GT. 0) WRITE(LUER,'(A,I5)')
	1      '  JFNK: SCHUR regularized pop rows=', N_REG_ROWS
	  END IF
	END IF
!
! ── Banded LU preconditioner (JFNK_USE_BANDED=T) ────────────────────
! Assemble the full NT*ND × NT*ND banded BA matrix in LAPACK banded
! storage format and factor with DGBTRF.  This captures cross-depth
! J-coupling from the off-diagonal BA bands (radiative VJ/dJ-dopacity
! coupling and advection), providing a better preconditioner than
! block-diagonal LU when cross-depth coupling is significant.
!
! Banded storage: AB(KL+KU+1+i-j, j) = A(i, j)
! KL = KU = (NUM_BNDS/2) * NT  (sub/super-diagonal bandwidth)
! LDAB = 2*KL + KU + 1
!
! Block (K, L) maps to global rows (K-1)*NT+1..K*NT, cols (L-1)*NT+1..L*NT
! where L = BND_TO_FULL(B, K) for band B at depth K.
!
! If DGBTRF fails, fall back to the per-depth LU/Jacobi above.
!
	JFNK_BANDED_LU_READY = .FALSE.
	IF(JFNK_USE_BANDED_LOC .AND. JFNK_NUM_BNDS .GT. 1)THEN
! Block-tridiagonal bandwidth: each off-diagonal block is NT×NT, so
! the maximum |i-j| for the lower band is (NUM_BNDS/2)*NT + (NT-1).
! The previous code used (NUM_BNDS/2)*NT, which is too small by NT-1
! and caused out-of-bounds writes in the banded assembly below.
! Fixed 2026-09-11 (overnight JFNK/DTAU work).
	  JFNK_BANDED_KL = (JFNK_NUM_BNDS/2) * JFNK_NT + JFNK_NT - 1
	  JFNK_BANDED_KU = JFNK_BANDED_KL
	  JFNK_BANDED_N = JFNK_NT * JFNK_ND
!
	  IOS=0
	  IF(ALLOCATED(JFNK_BANDED_LUAB))THEN
	    IF(SIZE(JFNK_BANDED_LUAB,1) .NE.
	1       2*JFNK_BANDED_KL+JFNK_BANDED_KU+1 .OR.
	2       SIZE(JFNK_BANDED_LUAB,2) .NE. JFNK_BANDED_N)THEN
	      DEALLOCATE(JFNK_BANDED_LUAB)
	    END IF
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_BANDED_LUAB))THEN
	    ALLOCATE(JFNK_BANDED_LUAB(
	1       2*JFNK_BANDED_KL+JFNK_BANDED_KU+1, JFNK_BANDED_N),
	2       STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc BANDED_LUAB failed'
	  END IF
	  IF(ALLOCATED(JFNK_BANDED_IPIV))THEN
	    IF(SIZE(JFNK_BANDED_IPIV) .NE. JFNK_BANDED_N)
	1       DEALLOCATE(JFNK_BANDED_IPIV)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_BANDED_IPIV))THEN
	    ALLOCATE(JFNK_BANDED_IPIV(JFNK_BANDED_N), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc BANDED_IPIV failed'
	  END IF
	  JFNK_BANDED_LUAB = 0.0_LDP
!
! Assemble banded BA from GENERATE_FULL_MATRIX_V3 (all bands).
! For each depth K, first call the diagonal band with
! USE_PASSED_REP=.FALSE. to determine the correct REPLACE for
! that depth (which ground-state equations are replaced by
! ionization equations).  Then call all bands with
! USE_PASSED_REP=.TRUE. so the same REPLACE is used for all
! bands at this depth.  The off-diagonal bands capture the
! cross-depth J-coupling (dJ/dopacity via VJ) and advection
! coupling (BA_ADV_TERM) from the BA matrix.
!
	  DO K = 1, JFNK_ND
	    REPLACE = REPLACE_SAVE(:, K)
	    DO IB = 1, JFNK_NUM_BNDS
	      BAND_OFFSET = K + IB - JFNK_DIAG_INDX
	      L = BAND_OFFSET
	      IF(L .LT. 1 .OR. L .GT. JFNK_ND) CYCLE
	      CALL GENERATE_FULL_MATRIX_V3(
	1         JFNK_C_MAT, DUMMY_STEQ, POPS_WORK, REPLACE, ZERO_STEQ,
	1         JFNK_NT, JFNK_ND, JFNK_NION, JFNK_NUM_BNDS,
	1         IB, JFNK_DIAG_INDX, K,
	1         .FALSE., (K .EQ. JFNK_ND .AND. IB .EQ. JFNK_NUM_BNDS),
	2         .TRUE.)
!
! Same PTC fold as the block-LU path so the banded factorization
! matches the FD operator's s-shifted diagonal.
!
	      IF(IB .EQ. JFNK_DIAG_INDX .AND. JFNK_PTC_FOLDED)THEN
	        DO JJ = 1, JFNK_NT
	          IDX_PTC = (K-1)*JFNK_NT + JJ
	          JFNK_C_MAT(JJ,JJ) = JFNK_C_MAT(JJ,JJ) +
	1           MAX(ABS(jfnk_fx(IDX_PTC)), 1.0D-8*FMAX_PTC)
	        END DO
	      END IF
!
! Copy NT×NT block into LAPACK banded storage.
! Global row i = (K-1)*NT + II, global col j = (L-1)*NT + JJ.
! AB(KL+KU+1+i-j, j) = A(i, j).
!
	      DO JJ = 1, JFNK_NT
	        DO II = 1, JFNK_NT
		  IROW = JFNK_BANDED_KL+JFNK_BANDED_KU+1+
	1           (K-L)*JFNK_NT + (II-JJ)
		  JCOL = (L-1)*JFNK_NT + JJ
		  IF(IROW .GE. 1 .AND. IROW .LE.
	1         SIZE(JFNK_BANDED_LUAB,1) .AND.
	2         JCOL .GE. 1 .AND. JCOL .LE.
	3         SIZE(JFNK_BANDED_LUAB,2))THEN
		    JFNK_BANDED_LUAB(IROW, JCOL) = JFNK_C_MAT(II, JJ)
		  END IF
	        END DO
	      END DO
	    END DO
	  END DO
!
! Equilibrate the banded matrix before factorization.
! The full BA has extreme dynamic range (row scale 1E-12 to 1E+42),
! so equilibration is essential for the banded LU to produce a
! useful preconditioner.  Without it, the banded matrix is nearly
! singular and the Newton step is in the null space.
!
! DGBEQU computes ROW_SF and COL_SF such that
!   A_eq(i,j) = ROW_SF(i) * A(i,j) * COL_SF(j)
! has row and column norms close to 1.
!
! DGBEQU fails (info>0) if any row or column is entirely zero.
! The BA has zero rows for fixed variables (T, NE) — set their
! diagonal to 1.0 before equilibration so DGBEQU succeeds.  These
! rows produce x_i = b_i which is then masked to zero by
! JFNK_MASK_FIXED_UNKNOWNS.
!
	  BLOCK
	    INTEGER :: N_ZROW, N_ZCOL, J_Z, I_Z, JROW, JCOL
	    INTEGER :: N_REG_B, K_B, J_B
	    REAL(KIND=LDP) :: RMAX_B, BLKMAX_B
	    LOGICAL :: ROW_NZ
	    N_ZROW = 0; N_ZCOL = 0; N_REG_B = 0
!
! Block maximum over the assembled banded matrix (headroom rows are
! zeroed at allocation).  Rows whose max|A_ij| is below
! JFNK_ROW_TINY*blockmax are replaced by the identity row, matching
! the per-depth LU path.  Preconditioner only.  Repaired rows/columns
! are recorded in JFNK_ROW_REG so the apply side projects them to
! zero (linear) instead of returning z_i ~ v_i.
! edited by Sean Brennan (2026-10)
!
	    IF(ALLOCATED(JFNK_ROW_REG))THEN
	      IF(SIZE(JFNK_ROW_REG,1).NE.JFNK_NT .OR.
	1         SIZE(JFNK_ROW_REG,2).NE.JFNK_ND) DEALLOCATE(JFNK_ROW_REG)
	    END IF
	    IF(.NOT. ALLOCATED(JFNK_ROW_REG))THEN
	      ALLOCATE(JFNK_ROW_REG(JFNK_NT, JFNK_ND))
	    END IF
	    JFNK_ROW_REG = .FALSE.
	    BLKMAX_B = MAXVAL(ABS(JFNK_BANDED_LUAB))
! Count/fix zero ROWS: row I_Z has all-zero elements across columns
	    DO I_Z = 1, JFNK_BANDED_N
	      RMAX_B = 0.0_LDP
	      DO J_Z = MAX(1,I_Z-JFNK_BANDED_KL),
	1         MIN(JFNK_BANDED_N,I_Z+JFNK_BANDED_KU)
	        JROW = JFNK_BANDED_KL + JFNK_BANDED_KU + 1 + I_Z - J_Z
	        RMAX_B = MAX(RMAX_B, ABS(JFNK_BANDED_LUAB(JROW, J_Z)))
	      END DO
	      K_B = (I_Z-1)/JFNK_NT + 1
	      J_B = I_Z - (K_B-1)*JFNK_NT
	      IF(RMAX_B .EQ. 0.0_LDP)THEN
	        N_ZROW = N_ZROW + 1
	        JFNK_BANDED_LUAB(JFNK_BANDED_KL+JFNK_BANDED_KU+1, I_Z)
	1           = 1.0_LDP
	        JFNK_ROW_REG(J_B, K_B) = .TRUE.
	      ELSE IF(RMAX_B .LT. JFNK_ROW_TINY*BLKMAX_B)THEN
	        DO J_Z = MAX(1,I_Z-JFNK_BANDED_KL),
	1           MIN(JFNK_BANDED_N,I_Z+JFNK_BANDED_KU)
	          JROW = JFNK_BANDED_KL + JFNK_BANDED_KU + 1 + I_Z - J_Z
	          JFNK_BANDED_LUAB(JROW, J_Z) = 0.0_LDP
	        END DO
	        JFNK_BANDED_LUAB(JFNK_BANDED_KL+JFNK_BANDED_KU+1, I_Z)
	1           = 1.0_LDP
	        JFNK_ROW_REG(J_B, K_B) = .TRUE.
	        N_REG_B = N_REG_B + 1
	      END IF
	    END DO
! Count/fix zero COLUMNS: column J_Z has all-zero elements across rows
	    DO J_Z = 1, JFNK_BANDED_N
	      ROW_NZ = .FALSE.
	      DO I_Z = MAX(1,J_Z-JFNK_BANDED_KU),
	1         MIN(JFNK_BANDED_N,J_Z+JFNK_BANDED_KL)
	        JROW = JFNK_BANDED_KL + JFNK_BANDED_KU + 1 + I_Z - J_Z
	        IF(JFNK_BANDED_LUAB(JROW, J_Z) .NE. 0.0_LDP)
	1           ROW_NZ = .TRUE.
	      END DO
	      IF(.NOT. ROW_NZ)THEN
	        N_ZCOL = N_ZCOL + 1
	        JFNK_BANDED_LUAB(JFNK_BANDED_KL+JFNK_BANDED_KU+1, J_Z)
	1           = 1.0_LDP
	        K_B = (J_Z-1)/JFNK_NT + 1
	        J_B = J_Z - (K_B-1)*JFNK_NT
	        JFNK_ROW_REG(J_B, K_B) = .TRUE.
	      END IF
	    END DO
	    IF(LUER .GT. 0 .AND. (N_ZROW .GT. 0 .OR. N_ZCOL .GT. 0))THEN
	      WRITE(LUER,'(A,I6,A,I6)')
	1       '  JFNK: banded zero rows/cols fixed: rows=', N_ZROW,
	2       '  cols=', N_ZCOL
	    END IF
	    IF(LUER .GT. 0 .AND. N_REG_B .GT. 0)THEN
	      WRITE(LUER,'(A,I6,A,ES10.2,A)')
	1       '  JFNK: regularized', N_REG_B,
	2       ' near-null banded rows (JFNK_ROW_TINY=',
	3       JFNK_ROW_TINY, ')'
	    END IF
	  END BLOCK
!
	  IF(ALLOCATED(JFNK_BANDED_ROW_SF))THEN
	    IF(SIZE(JFNK_BANDED_ROW_SF,1).NE.JFNK_BANDED_N)
	1     DEALLOCATE(JFNK_BANDED_ROW_SF)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_BANDED_ROW_SF))THEN
	    ALLOCATE(JFNK_BANDED_ROW_SF(JFNK_BANDED_N), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc BANDED_ROW_SF'
	  END IF
	  IF(ALLOCATED(JFNK_BANDED_COL_SF))THEN
	    IF(SIZE(JFNK_BANDED_COL_SF,1).NE.JFNK_BANDED_N)
	1     DEALLOCATE(JFNK_BANDED_COL_SF)
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_BANDED_COL_SF))THEN
	    ALLOCATE(JFNK_BANDED_COL_SF(JFNK_BANDED_N), STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc BANDED_COL_SF'
	  END IF
!
! Equilibrate the banded matrix.  IMPORTANT (fix 2026-09-11):
! DGBEQU uses the COMPACT banded layout AB(KU+1+i-j,j) with
! LDAB=KL+KU+1, while DGBTRF/DGBTRS use the expanded layout
! AB(KL+KU+1+i-j,j) with LDAB=2*KL+KU+1 (KL rows of pivoting
! headroom).  Calling DGBEQU on a DGBTRF-layout array reads the
! matrix shifted up by KL rows — it then reports "row i is exactly
! zero" (info=i) for rows that are zero only in the headroom, and
! equilibration is silently skipped, leaving DGBTRF to factor a
! matrix with ~1E54 dynamic range.  The resulting LU is numerically
! meaningless and GMRES steps land in the near-null space.
!
! The correct options are (a) pass the shifted base pointer
! LUAB(KL+1,1) so DGBEQU's KU+1+i-j lands on KL+KU+1+i-j, or
! (b) compute the row/column maxima directly in DGBTRF layout.
! Option (b) is used here: it is transparent, O(N*bandwidth), and
! avoids sequence-association subtleties.
!
! Algorithm (same as DGBEQU): R(i)=1/max_j|A(i,j)| on the raw band,
! then C(j)=1/max_i|R(i)*A(i,j)| on the row-scaled band.  The solve
! then uses b'=R*b, A_eq=D_r*A*D_c, x=C*y (already coded in
! JFNK_PRECOND_V1).
!
	  BLOCK
	    INTEGER :: I_E, J_E, EROW, N_ZROW
	    REAL(KIND=LDP) :: RMAX, CMAX
	    REAL(KIND=LDP) :: RMIN, RMAXF, CMIN, CMAXF
	    RMIN = 1.0E300_LDP; RMAXF = 0.0_LDP
	    CMIN = 1.0E300_LDP; CMAXF = 0.0_LDP
	    N_ZROW = 0
	    DO I_E = 1, JFNK_BANDED_N
	      RMAX = 0.0_LDP
	      DO J_E = MAX(1,I_E-JFNK_BANDED_KL),
	1         MIN(JFNK_BANDED_N,I_E+JFNK_BANDED_KU)
	        EROW = JFNK_BANDED_KL + JFNK_BANDED_KU + 1 + I_E - J_E
	        RMAX = MAX(RMAX, ABS(JFNK_BANDED_LUAB(EROW, J_E)))
	      END DO
	      IF(RMAX .GT. 0.0_LDP)THEN
	        JFNK_BANDED_ROW_SF(I_E) = 1.0_LDP / RMAX
	      ELSE
	        JFNK_BANDED_ROW_SF(I_E) = 1.0_LDP
	        JFNK_BANDED_LUAB(JFNK_BANDED_KL+JFNK_BANDED_KU+1, I_E)
	1           = 1.0_LDP
	        N_ZROW = N_ZROW + 1
	      END IF
	      RMIN = MIN(RMIN, RMAX); RMAXF = MAX(RMAXF, RMAX)
	    END DO
	    DO J_E = 1, JFNK_BANDED_N
	      CMAX = 0.0_LDP
	      DO I_E = MAX(1,J_E-JFNK_BANDED_KU),
	1         MIN(JFNK_BANDED_N,J_E+JFNK_BANDED_KL)
	        EROW = JFNK_BANDED_KL + JFNK_BANDED_KU + 1 + I_E - J_E
	        CMAX = MAX(CMAX,
	1         JFNK_BANDED_ROW_SF(I_E)*ABS(JFNK_BANDED_LUAB(EROW, J_E)))
	      END DO
	      IF(CMAX .GT. 0.0_LDP)THEN
	        JFNK_BANDED_COL_SF(J_E) = 1.0_LDP / CMAX
	      ELSE
	        JFNK_BANDED_COL_SF(J_E) = 1.0_LDP
	      END IF
	      CMIN = MIN(CMIN, CMAX); CMAXF = MAX(CMAXF, CMAX)
	    END DO
	    JFNK_BANDED_EQUIL_OK = .TRUE.
	    IF(LUER .GT. 0)THEN
	      WRITE(LUER,'(A,I6,A,ES10.2,A,ES10.2)')
	1       '  JFNK: banded equil (manual) zero_rows=', N_ZROW,
	2       '  rowmax=', RMIN, '-', RMAXF
	      WRITE(LUER,'(A,ES10.2,A,ES10.2)')
	1       '  JFNK: banded equil scaled colmax=', CMIN, '-', CMAXF
	    END IF
	  END BLOCK
!
! Save an unscaled copy of the assembled matrix for the
! solve-residual check below (post zero-row fix, pre-scaling).
!
	  IF(ALLOCATED(JFNK_BANDED_ARAW))THEN
	    IF(SIZE(JFNK_BANDED_ARAW,1).NE.
	1       2*JFNK_BANDED_KL+JFNK_BANDED_KU+1 .OR.
	2       SIZE(JFNK_BANDED_ARAW,2).NE.JFNK_BANDED_N)THEN
	      DEALLOCATE(JFNK_BANDED_ARAW)
	    END IF
	  END IF
	  IF(.NOT. ALLOCATED(JFNK_BANDED_ARAW))THEN
	    ALLOCATE(JFNK_BANDED_ARAW(
	1       2*JFNK_BANDED_KL+JFNK_BANDED_KU+1, JFNK_BANDED_N),
	2       STAT=IOS)
	    IF(IOS.NE.0) STOP 'JFNK_FACTOR_PRECOND: alloc ARAW failed'
	  END IF
	  JFNK_BANDED_ARAW = JFNK_BANDED_LUAB
!
	  IF(JFNK_BANDED_EQUIL_OK)THEN
! Apply equilibration: A_eq(i,j) = ROW_SF(i) * A(i,j) * COL_SF(j)
	    DO J = 1, JFNK_BANDED_N
	      DO I = MAX(1,J-JFNK_BANDED_KU),
	1         MIN(JFNK_BANDED_N,J+JFNK_BANDED_KL)
	        IROW = JFNK_BANDED_KL + JFNK_BANDED_KU + 1 + I - J
	        JFNK_BANDED_LUAB(IROW, J) =
	1           JFNK_BANDED_LUAB(IROW, J) *
	2           JFNK_BANDED_ROW_SF(I) * JFNK_BANDED_COL_SF(J)
	      END DO
	    END DO
	  END IF
!
! Factor the banded matrix.  Try zero ridge first; if DGBTRF fails,
! add a global ridge proportional to the maximum diagonal element
! and retry once.  This avoids unnecessary regularization when the
! matrix is already well-conditioned.
!
	  N_BANDED_RETRY = 0
	  DO
	    CALL DGBTRF(JFNK_BANDED_N, JFNK_BANDED_N,
	1       JFNK_BANDED_KL, JFNK_BANDED_KU,
	2       JFNK_BANDED_LUAB,
	3       2*JFNK_BANDED_KL+JFNK_BANDED_KU+1,
	4       JFNK_BANDED_IPIV, BANDED_INFO)
	    IF(BANDED_INFO .EQ. 0)THEN
	      JFNK_BANDED_LU_READY = .TRUE.
	      IF(LUER .GT. 0)THEN
	        WRITE(LUER,'(A,I6,A,I4,A,I4,A,I2)')
	1       '  JFNK: banded LU (DGBTRF)  N=', JFNK_BANDED_N,
	2       '  KL=', JFNK_BANDED_KL, '  KU=', JFNK_BANDED_KU,
	3       '  retry=', N_BANDED_RETRY
	      END IF
!
! Solve-residual check on EVERY factorization (audit fix F4,
! 2026-09-11) — a failed check now REJECTS the banded LU: solve A_eq y = D_r b on a
! synthetic b via the fresh LU, recover x = D_c y, then compute
! r = A_raw*x - b with a banded matvec on the unscaled matrix.
! Reported as the componentwise backward error
!   max_i |r_i| / (|A||x| + |b|)_i
! which is ~1E-16..1E-12 for a correct LU and O(0.01..1) for a
! broken assembly or equilibration — the raw ||r||/||b|| form is
! not usable on this matrix (row scale 1E-44..1E+11 => ||x|| huge,
! absolute residual astronomically large even for a perfect LU).
! This is the direct algebraic test that the previous run lacked
! (it only saw the downstream non-descent symptom).
!
	      BLOCK
	          DOUBLE PRECISION, ALLOCATABLE :: BT(:), XT(:)
	          DOUBLE PRECISION :: RN, BN, REL_RES, RI
	          INTEGER :: I_C, J_C, EROW_C, TRS_INFO
	          ALLOCATE(BT(JFNK_BANDED_N), XT(JFNK_BANDED_N))
	          BT = 1.0D0
	          XT = BT
	          IF(JFNK_BANDED_EQUIL_OK)THEN
	            DO I_C = 1, JFNK_BANDED_N
	              XT(I_C) = XT(I_C) * DBLE(JFNK_BANDED_ROW_SF(I_C))
	            END DO
	          END IF
	          CALL DGBTRS('N', JFNK_BANDED_N, JFNK_BANDED_KL,
	1             JFNK_BANDED_KU, 1, JFNK_BANDED_LUAB,
	2             2*JFNK_BANDED_KL+JFNK_BANDED_KU+1,
	3             JFNK_BANDED_IPIV, XT, JFNK_BANDED_N, TRS_INFO)
	          IF(TRS_INFO .EQ. 0 .AND. JFNK_BANDED_EQUIL_OK)THEN
	            DO I_C = 1, JFNK_BANDED_N
	              XT(I_C) = XT(I_C) * DBLE(JFNK_BANDED_COL_SF(I_C))
	            END DO
	          END IF
! Componentwise backward error (Oettli-Prager): r_i / (|A||x| + |b|)_i.
! The raw ||Ax-b||/||b|| is meaningless here: with rows spanning
! 1E-44..1E+11, ||x||~1E+44 for b=1, so even a perfect LU shows a
! huge absolute residual.  The componentwise form is ~machine
! epsilon for a correct factorization regardless of scaling.
	          RN = 0.0D0; BN = 0.0D0
	          DO I_C = 1, JFNK_BANDED_N
	            RI = -BT(I_C)
	            REL_RES = 0.0D0
	            DO J_C = MAX(1,I_C-JFNK_BANDED_KL),
	1               MIN(JFNK_BANDED_N,I_C+JFNK_BANDED_KU)
	              EROW_C = JFNK_BANDED_KL+JFNK_BANDED_KU+1+I_C-J_C
	              RI = RI + DBLE(JFNK_BANDED_ARAW(EROW_C,J_C))*XT(J_C)
	              REL_RES = REL_RES +
	1               ABS(DBLE(JFNK_BANDED_ARAW(EROW_C,J_C))*XT(J_C))
	            END DO
	            REL_RES = REL_RES + ABS(BT(I_C))
	            RN = MAX(RN, ABS(RI)/MAX(REL_RES,1.0D-300))
	            BN = MAX(BN, ABS(BT(I_C)))
	          END DO
	          REL_RES = RN
	          IF(LUER .GT. 0)THEN
	            WRITE(LUER,'(A,ES10.2,A,I3)')
	1             '  JFNK: banded LU residual check  '//
	2             'max_i |r_i|/(|A||x|+|b|)_i =',
	3             REL_RES,'  DGBTRS info=',TRS_INFO
	            FLUSH(LUER)
	          END IF
	          IF(TRS_INFO .NE. 0 .OR. REL_RES .GT. 1.0D-4)THEN
	            JFNK_BANDED_LU_READY = .FALSE.
	            IF(LUER .GT. 0)THEN
	              WRITE(LUER,'(A)')
	1         '  JFNK: banded LU residual check FAILED —'
	1         //' retrying with ridge'
	              FLUSH(LUER)
	            END IF
	          END IF
	          DEALLOCATE(BT, XT)
	        END BLOCK
! Only accept the factorization if the residual check passed.
! A rejected (or DGBTRS-failed) solve falls through to the ridge-retry
! path below so the degenerate ~zero rows are regularised — exactly as
! the per-depth path does — rather than immediately giving up.
! (Audit fix F6, 2026-09-11)
	      IF(JFNK_BANDED_LU_READY)THEN
! PTC post-scaling must use the diagonal of the accepted banded
! matrix, including row repairs and ridge, not fallback-LU state.
	        JFNK_BA_RIDGE_ADD=0.0_LDP
	        IF(N_BANDED_RETRY .GT. 0)
	1          JFNK_BA_RIDGE_ADD=BANDED_RIDGE
	        DO K=1,JFNK_ND
	          DO J=1,JFNK_NT
	            IROW=JFNK_BANDED_KL+JFNK_BANDED_KU+1
	            JCOL=(K-1)*JFNK_NT+J
	            JFNK_PRECOND_DIAG(J,K)=
	1              JFNK_BANDED_ARAW(IROW,JCOL)
	          END DO
	        END DO
	        EXIT
	      END IF
	    END IF
	    N_BANDED_RETRY = N_BANDED_RETRY + 1
	    IF(N_BANDED_RETRY .GT. 2)THEN
	      IF(LUER .GT. 0)THEN
	        WRITE(LUER,'(A,I6,A)')
	1       '  JFNK: banded LU FAILED (info=', BANDED_INFO,
	2       ') — falling back to per-depth LU/Jacobi'
	      END IF
	      EXIT
	    END IF
!
! DGBTRF failed — restore the pristine assembled matrix and retry
! with a global ridge.  DGBTRF overwrites LUAB in place with partial
! LU factors, so adding a ridge to the factorized array (the old
! code) re-factored garbage on the retry — exactly the ill-conditioned
! case the retry exists for.  Rebuild sequence:
!   LUAB = ARAW (pristine, unscaled, post zero-row fix)
!   compute max diagonal from the restored raw matrix
!   add ridge to the raw diagonal
!   ARAW = LUAB so the residual check tests the matrix factored
!   re-apply ROW_SF/COL_SF equilibration (destroyed by DGBTRF)
!   loop back to DGBTRF.
! (Audit fix F2, 2026-09-11)
!
	    JFNK_BANDED_LUAB = JFNK_BANDED_ARAW
	    BANDED_MAX_DIAG = 0.0_LDP
	    DO K = 1, JFNK_ND
	      DO I = 1, JFNK_NT
	        BANDED_MAX_DIAG = MAX(BANDED_MAX_DIAG,
	1         ABS(JFNK_BANDED_LUAB(
	2           JFNK_BANDED_KL+JFNK_BANDED_KU+1,
	3           (K-1)*JFNK_NT + I)))
	      END DO
	    END DO
	    BANDED_RIDGE = MAX(JFNK_BA_RIDGE_FAC,
	1       1.0E-10_LDP) * MAX(BANDED_MAX_DIAG, 1.0_LDP)
	    IF(LUER .GT. 0)THEN
	      WRITE(LUER,'(A,I2,A,ES10.2)')
	1       '  JFNK: banded LU retry ', N_BANDED_RETRY,
	2       '  ridge=', BANDED_RIDGE
	    END IF
	    DO K = 1, JFNK_ND
	      DO I = 1, JFNK_NT
	        JFNK_BANDED_LUAB(
	1           JFNK_BANDED_KL+JFNK_BANDED_KU+1,
	2           (K-1)*JFNK_NT + I) =
	3           JFNK_BANDED_LUAB(
	4             JFNK_BANDED_KL+JFNK_BANDED_KU+1,
	5             (K-1)*JFNK_NT + I) + BANDED_RIDGE
	      END DO
	    END DO
	    JFNK_BANDED_ARAW = JFNK_BANDED_LUAB
	    IF(JFNK_BANDED_EQUIL_OK)THEN
	      DO J = 1, JFNK_BANDED_N
	        DO I = MAX(1,J-JFNK_BANDED_KU),
	1           MIN(JFNK_BANDED_N,J+JFNK_BANDED_KL)
	          IROW = JFNK_BANDED_KL + JFNK_BANDED_KU + 1 + I - J
	          JFNK_BANDED_LUAB(IROW, J) =
	1           JFNK_BANDED_LUAB(IROW, J) *
	2           JFNK_BANDED_ROW_SF(I) * JFNK_BANDED_COL_SF(J)
	        END DO
	      END DO
	    END IF
	  END DO
	END IF
!
	JFNK_BA_DIAG_READY = .TRUE.
	JFNK_N_LU_FACTORS = JFNK_N_LU_FACTORS + 1
	JFNK_FACTOR_DONE = .TRUE.
!
	END SUBROUTINE JFNK_FACTOR_PRECOND_V1


!
! JFNK_PRECOND_V1: Apply diagonal (Jacobi) BA preconditioner.
! z = V / diag(BA + ridge).  Zero entries where |z| > ZMAX.
!
! editted by Sean Brennan (08-Sep-2026)
!
	SUBROUTINE JFNK_PRECOND_V1(N, V)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_n_zmax_hit, jfnk_n_reg_proj,
	1     jfnk_pc_clip_event, jfnk_x, jfnk_log_unknown, jfnk_ptc_active,
	2     jfnk_ptc_dtau, jfnk_ptc_scale, jfnk_ptc_scale_ready
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_BA_RIDGE_FAC, JFNK_PRECOND_TYPE
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(INOUT) :: V(N)
!
	REAL(KIND=LDP) DIAG_USE, Jacobi
	INTEGER K, J, NJAC, NZERO, LU_INFO, NZMAX_LOC
	INTEGER NZCLIP, IDX_PC, K_PC, J_PC, NPROJ
	INTEGER I, IOS
	DOUBLE PRECISION Z_BA, X_PC, D_PC, SIG_PC, RHS_PC
	LOGICAL BANDED_APPLIED
	CHARACTER(LEN=1) NO_TRANS
	REAL(KIND=LDP) WI_LSQ
	REAL(KIND=LDP), ALLOCATABLE :: W_LSQ(:,:)
!
! STEQ_MAT and V_SAVE are now persistent module workspace
! (JFNK_W_STEQ_MAT, JFNK_W_V_SAVE) to avoid stack allocation.
! audit §3, 09-Sep-2026.
!
	NJAC = 0
	NZERO = 0
	NZMAX_LOC = 0
	NPROJ = 0
	BANDED_APPLIED = .FALSE.
!
! NZCLIP counts NONLINEAR preconditioner events: nonfinite->0 writes.
! Zeroing a nonfinite component is the only safe action, but it does
! make M^{-1} depend on the input vector, so the GMRES Krylov
! construction (which assumes a linear operator) is invalid -- the
! event flags jfnk_pc_clip_event and GMRES aborts the solve rather
! than silently solving a different problem.  Finite components of
! any magnitude pass through unclipped (revised 08-Oct-2026): with
! factor-time null rows projected by JFNK_ROW_REG, a huge finite z
! is faithful output of a badly conditioned block, and clipping it
! would itself be the nonlinearity.  Deterministic small-diagonal
! zeroing is linear (independent of v) and is not counted.
!
	NZCLIP = 0
!
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, V)
!
! The preconditioner applies plain M^{-1} (physical-space BA solve).
! RSCALE must NOT be applied here.  In weighted-norm GMRES the weight
! (RSCALE) only enters the inner product (JFNK_DOT_V1); the
! preconditioner operates on physical-space vectors, matching the
! physical-space matvec (JFNK_JVP_V1) and residual (JFNK_RESIDUAL_V1).
!
! The earlier "R*M^{-1}*R^{-1}" transform (audit Finding #5 suggestion)
! was incorrect: it is only valid if ALL GMRES vectors are scaled by
! RSCALE (fully scaled space with Euclidean inner product).  Since the
! code keeps vectors physical and only weights the inner product, the
! correct preconditioner is plain M^{-1}.  The R*M^{-1}*R^{-1} form
! corrupted the Krylov space and caused JFNK to fail on the first
! Newton step (runaway |s|, non-descent, or trust-region collapse).
!
! editted by Sean Brennan (09-Sep-2026)
!
! Task 3: Parallelize the depth-wise unpack/compute/pack loops.
! Each K writes a distinct column of STEQ_MAT / V_SAVE and a distinct
! block of V, so there are no inter-iteration data dependencies.
!
!$OMP PARALLEL DO PRIVATE(J) SHARED(JFNK_W_STEQ_MAT, JFNK_W_V_SAVE, V)
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT
	    JFNK_W_STEQ_MAT(J, K) = REAL(V((K-1)*JFNK_NT + J), KIND=LDP)
	    JFNK_W_V_SAVE(J, K) = JFNK_W_STEQ_MAT(J, K)
	  END DO
	END DO
!$OMP END PARALLEL DO
!
	IF(.NOT. JFNK_FACTOR_DONE)THEN
	  CALL JFNK_FACTOR_PRECOND_V1()
	END IF
!
! ── Banded LU solve (JFNK_USE_BANDED=T) ────────────────────────────
! When the banded LU factorization succeeded, solve the full
! NT*ND system with DGBTRS.  This captures cross-depth J-coupling
! (dJ/dopacity and advection) in the preconditioner, improving
! GMRES convergence when the off-diagonal BA bands are significant.
!
	IF(JFNK_BANDED_LU_READY)THEN
!
! Apply equilibration scaling: b' = ROW_SF * b, then solve
! A_eq * x' = b' with DGBTRS, then x = COL_SF * x'.
!
	  IF(JFNK_BANDED_EQUIL_OK)THEN
	    DO J = 1, JFNK_BANDED_N
	      V(J) = V(J) * DBLE(JFNK_BANDED_ROW_SF(J))
	    END DO
	  END IF
	  CALL DGBTRS('N', JFNK_BANDED_N,
	1     JFNK_BANDED_KL, JFNK_BANDED_KU, 1,
	2     JFNK_BANDED_LUAB,
	3     2*JFNK_BANDED_KL+JFNK_BANDED_KU+1,
	4     JFNK_BANDED_IPIV, V, JFNK_BANDED_N, LU_INFO)
	  IF(LU_INFO .EQ. 0)THEN
	    IF(JFNK_BANDED_EQUIL_OK)THEN
	      DO J = 1, JFNK_BANDED_N
	        V(J) = V(J) * DBLE(JFNK_BANDED_COL_SF(J))
	      END DO
	    END IF
! Project regularized null rows/cols to zero -- a linear map, so
! M^{-1} stays consistent.  (The banded path has no |z| clip at all,
! so without this a repaired row returned z_i ~ v_i, which can be
! ~1e12 when the residual is large.)
	    IF(ALLOCATED(JFNK_ROW_REG))THEN
	      DO J = 1, JFNK_BANDED_N
	        K_PC = (J-1)/JFNK_NT + 1
	        J_PC = J - (K_PC-1)*JFNK_NT
	        IF(JFNK_ROW_REG(J_PC, K_PC))THEN
	          V(J) = 0.0_LDP
	          NPROJ = NPROJ + 1
	        END IF
	      END DO
	      jfnk_n_reg_proj = jfnk_n_reg_proj + NPROJ
	    END IF
	    CALL JFNK_MASK_FIXED_UNKNOWNS(N, V)
	    BANDED_APPLIED = .TRUE.
	  ELSE
! Banded solve failed -- fall through to per-depth LU/Jacobi.
	    IF(JFNK_LUER .GT. 0)THEN
	      WRITE(JFNK_LUER,'(A,I6)')
	1        '  JFNK: DGBTRS failed (info=', LU_INFO,
	2        ') -- Jacobi fallback'
	    END IF
	  END IF
	END IF
!
! editted by Devin, 05-Oct-2026: the banded branch used to RETURN
! above, skipping the shared post-processing at the end of the
! routine (log-unknown exp(-x) rescale, PTC diagonal correction,
! final fixed-unknown mask).  A successful DGBTRS now flags
! BANDED_APPLIED and falls through so a banded apply sees identical
! post-processing to the LU path; only the per-depth solve/pack
! block is bypassed because the banded solve already wrote V.
!
	IF(.NOT. BANDED_APPLIED)THEN
!
	NO_TRANS = 'N'
!
	IF((JFNK_PRECOND_TYPE(1:2) .EQ. 'LU' .OR.
	1    JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ')
	2    .AND. JFNK_LU_READY)THEN
!
! ── Full BA block LU solve with equilibration ─────────────────────
! For each depth K, solves (BA_full + ridge) * step = v exactly per
! depth block (up to LU roundoff).  The ROW_SF/COL_SF equilibration
! only improves the conditioning of the LU factorization (DGETRF)
! and does NOT change the linear system being solved.  The sequence
! is:
!   1. Scale RHS:  v' = ROW_SF * v
!   2. DGETRS:     LU_equilibrated * y = v'  →  y = (BA_equil)^{-1} v'
!   3. Scale sol:  step = COL_SF * y
! which recovers step = (BA_full + ridge)^{-1} * v exactly.
! (Finding 5, review 2026-09-10: corrected misleading comment that
!  described the preconditioner as an approximate scaled matrix; it
!  is the exact inverse of the ridge-augmented BA per depth block.)
!
! For 'LSQ' the factored operator is (J_eq^T J_eq + lam^2 I), so
! the solve is instead
!   1. w   = ROW_SF * v'                (= ROW_SF^2 v)
!   2. rhs = COL_SF * B^T w             (= J_eq^T v')
!   3. DPOTRS: y = (J_eq^T J_eq + lam^2 I)^{-1} rhs
!   4. step = COL_SF * y
! where B is the ridged, reg-substituted block before equilibration.
! B is never rebuilt: B^T w is accumulated from the stored raw block
! JFNK_BA_A plus the ridge diagonal and the reg-row corrections
! (verified against the direct solve in cmfgen_jfnk_unit_test lsqchk).
!
	IF(JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ')THEN
	  ALLOCATE(W_LSQ(JFNK_NT, JFNK_ND), STAT=IOS)
	  IF(IOS .NE. 0) STOP 'JFNK_PRECOND: alloc W_LSQ failed'
	END IF
!
!$OMP PARALLEL DO PRIVATE(J, I, LU_INFO, DIAG_USE, WI_LSQ) REDUCTION(+:NJAC, NZERO, NZMAX_LOC, NZCLIP, NPROJ)
	  DO K = 1, JFNK_ND
	    IF(.NOT. ALLOCATED(JFNK_LU_OK) .OR. .NOT. JFNK_LU_OK(K))THEN
	      DO J = 1, JFNK_NT
	        DIAG_USE = JFNK_BA_DIAG(J, K) + JFNK_BA_RIDGE_ADD(J, K)
	        IF(ALLOCATED(JFNK_ROW_REG) .AND.
	1           JFNK_ROW_REG(J, K))THEN
	          JFNK_W_V_SAVE(J, K) = 0.0_LDP
	          NPROJ = NPROJ + 1
	        ELSE IF(ABS(DIAG_USE) .GT. 1.0E-30_LDP)THEN
	          JFNK_W_V_SAVE(J, K) = JFNK_W_STEQ_MAT(J, K) / DIAG_USE
	        ELSE
	          JFNK_W_V_SAVE(J, K) = 0.0_LDP
	          NZERO = NZERO + 1
	        END IF
	        NJAC = NJAC + 1
	      END DO
	    ELSE
	      DO J = 1, JFNK_NT
	        JFNK_W_V_SAVE(J, K) = JFNK_W_V_SAVE(J, K)
	1         * JFNK_ROW_SF(J, K)
	      END DO
	      IF(JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ')THEN
! w = ROW_SF v' and rhs = B^T w.  B = BA_A + ridge + reg-row
! substitution is expanded as
!   B^T w = BA_A^T w + D_r w + sum_{i reg} w_i (e_i - BA_A(i,:)^T)
! - w_i D_r(i) e_i
! so no second NT x NT block is needed to form J_eq^T v'.
	        DO J = 1, JFNK_NT
	          W_LSQ(J, K) = JFNK_W_V_SAVE(J, K) *
	1           JFNK_ROW_SF(J, K)
	        END DO
	        CALL DGEMV('T', JFNK_NT, JFNK_NT, 1.0_LDP,
	1           JFNK_BA_A(1,1,K), JFNK_NT, W_LSQ(1,K), 1,
	2           0.0_LDP, JFNK_W_V_SAVE(1,K), 1)
	        DO J = 1, JFNK_NT
	          JFNK_W_V_SAVE(J, K) = JFNK_W_V_SAVE(J, K) +
	1           JFNK_BA_RIDGE_ADD(J, K) * W_LSQ(J, K)
	        END DO
	        DO I = 1, JFNK_NT
	          IF(JFNK_ROW_REG(I, K))THEN
	            WI_LSQ = W_LSQ(I, K)
	            DO J = 1, JFNK_NT
	              JFNK_W_V_SAVE(J, K) = JFNK_W_V_SAVE(J, K) -
	1               WI_LSQ * JFNK_BA_A(I, J, K)
	            END DO
	            JFNK_W_V_SAVE(I, K) = JFNK_W_V_SAVE(I, K) +
	1             WI_LSQ - WI_LSQ * JFNK_BA_RIDGE_ADD(I, K)
	          END IF
	        END DO
	        DO J = 1, JFNK_NT
	          JFNK_W_V_SAVE(J, K) = JFNK_W_V_SAVE(J, K) *
	1           JFNK_COL_SF(J, K)
	        END DO
	        CALL DPOTRS('L', JFNK_NT, 1, JFNK_LU_FAC(1,1,K),
	1           JFNK_NT, JFNK_W_V_SAVE(1,K), JFNK_NT, LU_INFO)
	      ELSE
	        CALL DGETRS(NO_TRANS, JFNK_NT, 1, JFNK_LU_FAC(1,1,K),
	1           JFNK_NT, JFNK_IPIV(1,K), JFNK_W_V_SAVE(1,K),
	2           JFNK_NT, LU_INFO)
	      END IF
	      IF(LU_INFO .NE. 0)THEN
	        DO J = 1, JFNK_NT
	          DIAG_USE = JFNK_BA_DIAG(J, K) + JFNK_BA_RIDGE_ADD(J, K)
	          IF(ALLOCATED(JFNK_ROW_REG) .AND.
	1             JFNK_ROW_REG(J, K))THEN
	            JFNK_W_V_SAVE(J, K) = 0.0_LDP
	            NPROJ = NPROJ + 1
	          ELSE IF(ABS(DIAG_USE) .GT. 1.0E-30_LDP)THEN
	            JFNK_W_V_SAVE(J, K) = JFNK_W_STEQ_MAT(J, K) / DIAG_USE
	          ELSE
	            JFNK_W_V_SAVE(J, K) = 0.0_LDP
	            NZERO = NZERO + 1
	          END IF
	        END DO
	      ELSE
	        DO J = 1, JFNK_NT
	          JFNK_W_V_SAVE(J, K) = JFNK_W_V_SAVE(J, K)
	1           * JFNK_COL_SF(J, K)
!
! editted by Devin, 04-Oct-2026: zero the fixed-unknown components
! (T under FIX_T/FIX_T range, Ne under FIX_NE) before the clip test.
! DGETRS mixes all NT rows of the depth block, so a masked row is a
! genuine solve output of arbitrary magnitude; it is discarded by
! JFNK_MASK_FIXED_UNKNOWNS at pack time and must not count as a
! nonlinear-preconditioner clip event.
!
! Regularized null rows: project the solved component to zero instead
! of returning z_i~v_i.  Projection is a linear map, so M^{-1} stays
! a consistent linear operator and the Krylov space is unaffected;
! those directions simply get no update.  Finite z of any magnitude
! now passes through unclipped: with null rows projected, a huge
! finite component is faithful M^{-1} output on a badly conditioned
! block, and clipping it is itself the nonlinearity GMRES cannot
! tolerate.  Only nonfinite components are zeroed and flagged.
!
	          IF(ALLOCATED(JFNK_ROW_REG) .AND.
	1             JFNK_ROW_REG(J, K))THEN
	            JFNK_W_V_SAVE(J, K) = 0.0_LDP
	            NPROJ = NPROJ + 1
	          ELSE IF((JFNK_FIXED_NE_LOC .AND. J .EQ. JFNK_NT-1) .OR.
	1            (J .EQ. JFNK_NT .AND. (JFNK_FIXED_T_LOC .OR.
	2             (JFNK_FIX_T_D_ST .GT. 0 .AND.
	3              K .GE. JFNK_FIX_T_D_ST .AND.
	4              K .LE. JFNK_FIX_T_D_END .AND.
	5              .NOT. JFNK_TAUTO_ACTIVE))))THEN
	            JFNK_W_V_SAVE(J, K) = 0.0_LDP
	          ELSE IF(.NOT. (ABS(JFNK_W_V_SAVE(J, K)) .LE.
	1                 HUGE(1.0_LDP)))THEN
	            JFNK_W_V_SAVE(J, K) = 0.0_LDP
	            NZERO = NZERO + 1
	            NZCLIP = NZCLIP + 1
	          ELSE IF(ABS(JFNK_W_V_SAVE(J, K)) .GT.
	1                 JFNK_PRECOND_ZMAX)THEN
	            NZMAX_LOC = NZMAX_LOC + 1
	          END IF
	          NJAC = NJAC + 1
	        END DO
	      END IF
	    END IF
	  END DO
!$OMP END PARALLEL DO
	jfnk_n_zmax_hit = jfnk_n_zmax_hit + NZMAX_LOC
	jfnk_n_reg_proj = jfnk_n_reg_proj + NPROJ
!
	ELSE IF(JFNK_PRECOND_TYPE(1:5) .EQ. 'SCHUR' .AND.
	1        JFNK_SCHUR_OK)THEN
!
! ── Two-level (n,T) Schur apply ───────────────────────────────────
! z_n = Jnn^-1 v_n (per-depth LSQ pop block); eliminate the T column
! into the dense Schur complement, solve S z_T = v_T - JTn z_n, then
! back-substitute z_n -= Jnn^-1 JnT z_T.  All arithmetic lives in
! JFNK_SCHUR_APPLY_V1 (jfnk_solver_v1.f90) so the unit-test driver
! runs the identical production code.
!
	CALL JFNK_SCHUR_APPLY_V1(JFNK_ND, JFNK_NT, JFNK_NT-2,
	1    JFNK_NT, JFNK_NT-1,
	2    JFNK_LU_FAC, JFNK_BA_A, JFNK_BA_RIDGE_ADD, JFNK_ROW_REG,
	3    JFNK_ROW_SF, JFNK_COL_SF,
	4    JFNK_SCHUR_S, JFNK_SCHUR_IPIV, V, LU_INFO)
	IF(LU_INFO .NE. 0 .AND. JFNK_LUER .GT. 0)THEN
	  WRITE(JFNK_LUER,'(A,I6)')
	1     '  JFNK: SCHUR apply failed (info=', LU_INFO
	END IF
!
! Copy the solved vector into V_SAVE with the same reg-row / fixed /
! nonfinite accounting as the LU/LSQ path, so the shared post-chain
! (exp(-x) rescale, PTC fix, mask, clip flag) sees identical input.
!
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT
	    IF(ALLOCATED(JFNK_ROW_REG) .AND.
	1       JFNK_ROW_REG(J, K))THEN
	      JFNK_W_V_SAVE(J, K) = 0.0_LDP
	      NPROJ = NPROJ + 1
	    ELSE IF(.NOT. (ABS(V((K-1)*JFNK_NT + J)) .LE.
	1            HUGE(1.0_LDP)))THEN
	      JFNK_W_V_SAVE(J, K) = 0.0_LDP
	      NZERO = NZERO + 1
	      NZCLIP = NZCLIP + 1
	    ELSE
	      JFNK_W_V_SAVE(J, K) =
	1         REAL(V((K-1)*JFNK_NT + J), KIND=LDP)
	      NJAC = NJAC + 1
	    END IF
	  END DO
	END DO
	jfnk_n_reg_proj = jfnk_n_reg_proj + NPROJ
!
	ELSE
!
! ── Jacobi (diagonal) preconditioner ──────────────────────────────
!
!$OMP PARALLEL DO PRIVATE(J, Jacobi, DIAG_USE) REDUCTION(+:NJAC, NZERO, NZMAX_LOC, NZCLIP, NPROJ)
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT
	    Jacobi=0.0_LDP
	    IF(JFNK_BA_DIAG_READY)THEN
! Preconditioner ridge (ChatGPT audit fix #4, 2026-09-09):
! JFNK_FACTOR_PRECOND_V1 builds the BA diagonal with a ridge
!   RIDGE_ADD = MAX(MAX_VAL,1) * JFNK_BA_RIDGE_FAC
! and stores it in JFNK_BA_RIDGE_ADD(K).  The previous code here
! independently used MAX(JFNK_BA_RIDGE_FAC,1E-2)*MAX(ABS(DIAG),1),
! so the actual preconditioner was NOT the matrix that was built:
! it imposed a 1% ridge even when the configured ridge was zero,
! and used a per-row |DIAG| scale instead of the global MAX_VAL
! scale.  Use the stored ridge so M^{-1} matches the factored M.
	      DIAG_USE = JFNK_BA_DIAG(J, K) + JFNK_BA_RIDGE_ADD(J, K)
	      IF(ABS(DIAG_USE) .GT. 1.0E-30_LDP)THEN
	        Jacobi = JFNK_W_V_SAVE(J, K) / DIAG_USE
	      END IF
	    END IF
	    IF(ALLOCATED(JFNK_ROW_REG) .AND. JFNK_ROW_REG(J, K))THEN
	      JFNK_W_STEQ_MAT(J, K) = 0.0_LDP
	      NJAC = NJAC + 1
	      NPROJ = NPROJ + 1
	    ELSE IF(.NOT. (ABS(Jacobi) .LE. HUGE(1.0_LDP)))THEN
	      JFNK_W_STEQ_MAT(J, K) = 0.0_LDP
	      NZERO = NZERO + 1
	      NZCLIP = NZCLIP + 1
	    ELSE
	      JFNK_W_STEQ_MAT(J, K) = Jacobi
	      NJAC = NJAC + 1
	      IF(ABS(Jacobi) .GT. JFNK_PRECOND_ZMAX) NZMAX_LOC = NZMAX_LOC + 1
	    END IF
	  END DO
	END DO
!$OMP END PARALLEL DO
	jfnk_n_zmax_hit = jfnk_n_zmax_hit + NZMAX_LOC
	jfnk_n_reg_proj = jfnk_n_reg_proj + NPROJ
!
	END IF
! Audit 17-Sep-2026 fix: a nonlinear preconditioner event (NaN->0 or
! |z|>1D8 clip) makes M^{-1} depend on the input, invalidating the
! Krylov linearity assumption.  Flag it; GMRES aborts (info=3).
	IF(NZCLIP .GT. 0)jfnk_pc_clip_event = .TRUE.
!
! Pack the preconditioned vector back into V.
!
	IF(ALLOCATED(W_LSQ)) DEALLOCATE(W_LSQ)
!
!$OMP PARALLEL DO PRIVATE(J)
	DO K = 1, JFNK_ND
	  DO J = 1, JFNK_NT
	    IF(((JFNK_PRECOND_TYPE(1:2) .EQ. 'LU' .OR.
	1        JFNK_PRECOND_TYPE(1:3) .EQ. 'LSQ')
	2        .AND. JFNK_LU_READY) .OR.
	3        (JFNK_PRECOND_TYPE(1:5) .EQ. 'SCHUR' .AND.
	4         JFNK_SCHUR_OK))THEN
	      V((K-1)*JFNK_NT + J) = DBLE(JFNK_W_V_SAVE(J, K))
	    ELSE
	      V((K-1)*JFNK_NT + J) = DBLE(JFNK_W_STEQ_MAT(J, K))
	    END IF
	  END DO
	END DO
!$OMP END PARALLEL DO
!
! No RSCALE multiplication here — see comment at top of subroutine.
! The preconditioner output is plain M^{-1} * V in physical space.
!
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, V)
!
	END IF   ! .NOT. BANDED_APPLIED
!
! Log unknown: BA was built for a fractional correction.  At x = 0
! that derivative equals dF/d(ln n).  Later, divide the BA solution
! by exp(x) so the chain rule is in the preconditioner.
!
! editted by Devin, 05-Oct-2026: pseudo-time continuation adds the
! fixed diagonal sigma_i = s_i/jfnk_ptc_dtau to the Jacobian (the
! G = F + (x - x_k)*s/dtau term in JFNK_RESIDUAL_V1), so the consistent
! preconditioner is (J + sigma)^{-1}.  The previous blend
!   z = (1-alp)*Z_BA + alp*b/sigma,  alp = sigma/(sigma+|d'|),
! is wrong in the scalar limit: with Z_BA = b/d' (d' is the signed
! preconditioner diagonal, including ridge or banded row repair) it
! returns 2*b/(d'+sigma) for d' > 0 and 0 for d' < 0.  BA rate rows are predominantly loss-dominated (d' < 0),
! so the blend annihilated the preconditioner output on much of the
! vector.  Rescale the BA solve instead:
!   z_i = Z_BA_i * d'_i / (d'_i + sigma_i),
! which is scalar-exact for either sign of d'.  jfnk_ptc_dtau ->
! infinity sends sigma -> 0 and the factor to 1 (pure BA solve);
! dtau -> 0 gives z ~ Z_BA*d'/sigma ~ b/sigma, the diagonal PTC
! solve.  s_i and jfnk_ptc_dtau are frozen for a whole GMRES solve,
! so the weight depends on the fixed diagonal and sigma only, not
! on v: M^{-1} stays linear (Krylov validity). When |d'+sigma| is
! near zero, keep Z_BA rather than dividing by a singular value. If
! d' itself is zero, use the original right-hand side over sigma; a
! null-row identity solve must not be annihilated by a zero scale.
!
	IF(jfnk_log_unknown .OR.
	1   (jfnk_ptc_active .AND. jfnk_ptc_scale_ready))THEN
	  DO IDX_PC = 1, N
	    Z_BA = V(IDX_PC)
	    K_PC = (IDX_PC-1)/JFNK_NT + 1
	    J_PC = IDX_PC - (K_PC-1)*JFNK_NT
!
! A regularized null row must stay zero: the |D_PC|<tiny branch below
! would otherwise return RHS/sigma with RHS still holding the input v
! (banded path never rewrites JFNK_W_V_SAVE), resurrecting z_i ~ v_i.
!
	    IF(ALLOCATED(JFNK_ROW_REG))THEN
	      IF(JFNK_ROW_REG(J_PC, K_PC))THEN
	        V(IDX_PC) = 0.0D0
	        CYCLE
	      END IF
	    END IF
	    IF(jfnk_log_unknown .AND. ALLOCATED(jfnk_x))THEN
	      IF(SIZE(jfnk_x) .EQ. N)THEN
	        X_PC = MIN(MAX(jfnk_x(IDX_PC), -40.0D0), 40.0D0)
	        Z_BA = Z_BA * EXP(-X_PC)
	      END IF
	    END IF
	    IF(jfnk_ptc_active .AND. jfnk_ptc_scale_ready .AND.
	1       JFNK_BA_DIAG_READY .AND. .NOT. JFNK_PTC_FOLDED)THEN
	      D_PC = DBLE(JFNK_PRECOND_DIAG(J_PC, K_PC))
	      SIG_PC = DBLE(jfnk_ptc_scale(IDX_PC)) /
	1       MAX(jfnk_ptc_dtau, 1.0D-8)
	      RHS_PC = DBLE(JFNK_W_V_SAVE(J_PC, K_PC))
	      IF(ABS(D_PC) .GT. 1.0D-30)THEN
	        IF(ABS(D_PC + SIG_PC) .GT.
	1           1.0D-30*MAX(1.0D0, ABS(D_PC)))THEN
	          V(IDX_PC) = Z_BA * D_PC / (D_PC + SIG_PC)
	        ELSE
	          V(IDX_PC) = Z_BA
	        END IF
	      ELSE IF(ABS(SIG_PC) .GT. 1.0D-30)THEN
! A null BA diagonal cannot supply b/d; use the PTC diagonal directly
! rather than multiplying an identity-row solve by a zero factor.
	        V(IDX_PC) = RHS_PC / SIG_PC
	      ELSE
	        V(IDX_PC) = Z_BA
	      END IF
	    ELSE
	      V(IDX_PC) = Z_BA
	    END IF
	  END DO
	  CALL JFNK_MASK_FIXED_UNKNOWNS(N, V)
	END IF
!
	IF(NZERO .GT. 0 .AND. JFNK_LUER .GT. 0 .AND.
	1   JFNK_N_PRECOND_APPLIES .LT. 8)THEN
	  WRITE(JFNK_LUER,'(A,A,A,I8,A,I8,A,ES10.2)')
	1    '  JFNK: precond ', JFNK_PRECOND_TYPE,
	1    ' applied=', NJAC,
	1    '  zeroed=', NZERO,
	1    '  ridge=', JFNK_BA_RIDGE_FAC
	END IF
!
! FIX_T_AUTO: under left preconditioning every Krylov basis vector
! passes through here, so scaling the T components of the output by
! the current ladder alpha down-weights the energy rows inside the
! preconditioned operator and biases GMRES toward population-space
! descent while T stays live (not zeroed).  Step-level damping in
! JFNK_TAUTO_DAMP_V1 is still the hard bound.
!
	CALL JFNK_TAUTO_SCALE_V1(N, V)
!
	JFNK_N_PRECOND_APPLIES = JFNK_N_PRECOND_APPLIES + 1
!
	END SUBROUTINE JFNK_PRECOND_V1




!
! Zero T and/or Ne Newton components when those variables are fixed.
! Called on residuals, JVPs, and preconditioner vectors so GMRES never
! sees a singular FIX_T / FIX_NE BA row as an unknown.
! editted by Sean Brennan (08-Sep-2026)
!
	SUBROUTINE JFNK_MASK_FIXED_UNKNOWNS(N, VEC)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(INOUT) :: VEC(N)
	INTEGER K, IDX
!
	DO K = 1, JFNK_ND
	  IF(JFNK_FIXED_NE_LOC)THEN
	    IDX = (K-1)*JFNK_NT + (JFNK_NT-1)
	    IF(IDX .GE. 1 .AND. IDX .LE. N) VEC(IDX) = 0.0D0
	  END IF
	  IF(JFNK_FIXED_T_LOC)THEN
	    IDX = (K-1)*JFNK_NT + JFNK_NT
	    IF(IDX .GE. 1 .AND. IDX .LE. N) VEC(IDX) = 0.0D0
	  ELSE IF(JFNK_FIX_T_D_ST .GT. 0 .AND.
	1         K .GE. JFNK_FIX_T_D_ST .AND.
	2         K .LE. JFNK_FIX_T_D_END .AND.
	3         .NOT. JFNK_TAUTO_ACTIVE)THEN
	    IDX = (K-1)*JFNK_NT + JFNK_NT
	    IF(IDX .GE. 1 .AND. IDX .LE. N) VEC(IDX) = 0.0D0
	  END IF
	END DO
!
	END SUBROUTINE JFNK_MASK_FIXED_UNKNOWNS




!
! Scale the T components of a preconditioned vector by the current
! FIX_T_AUTO ladder alpha.  Unlike the mask, the components stay
! nonzero: the energy equation keeps its weight in the preconditioned
! system, just reduced.  Only applied when JFNK_TAUTO_ACTIVE and
! JFNK_TAU_PCW; harmless no-op otherwise.
!
	SUBROUTINE JFNK_TAUTO_SCALE_V1(N, VEC)
	USE JFNK_CMFGEN_INT_MOD
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
	INTEGER :: N
	REAL(KIND=8) :: VEC(N)
	INTEGER :: K, IDX
	IF(.NOT. JFNK_TAUTO_ACTIVE)RETURN
	IF(.NOT. JFNK_TAU_PCW)RETURN
	IF(JFNK_TAU_DOM_ST .LE. 0)RETURN
	DO K = JFNK_TAU_DOM_ST, JFNK_TAU_DOM_END
	  IDX = (K-1)*JFNK_NT + JFNK_NT
	  IF(IDX .GE. 1 .AND. IDX .LE. N)
	1   VEC(IDX) = VEC(IDX) * DBLE(JFNK_TAU_ALPHA)
	END DO
	RETURN
	END SUBROUTINE JFNK_TAUTO_SCALE_V1


!
! JFNK_DOT_V1: Weighted dot product for the JFNK solver.
!
! Computes  dot = sum_i (RSCALE(i) * A(i)) * (RSCALE(i) * B(i))
!
! The per-element residual scaling RSCALE normalizes each residual
! element to O(1) so that all species contribute equally to the norm.
! This prevents the solver from being dominated by species with the
! largest STEQ values (e.g., H) while ignoring trace species (e.g., Fe).
!
! If JFNK_RSCALE_READY is .FALSE. (scaling not yet computed), the
! function falls back to the unweighted dot product.
!
! editted by Sean Brennan (V4, 2026-09)
! CORRECTION (2026-09-07): Added per-element residual scaling.
!
	DOUBLE PRECISION FUNCTION JFNK_DOT_V1(N, A, B)
	USE JFNK_CMFGEN_INT_MOD, ONLY: JFNK_RSCALE, JFNK_RSCALE_READY
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(IN) :: A(N), B(N)
	INTEGER I
	DOUBLE PRECISION S
!
	S = 0.0D0
	IF(JFNK_RSCALE_READY)THEN
! Task 4: SIMD reduction for the weighted inner product.  The SA/SB
! temporaries are inlined so the reduction expression is explicit and
! the compiler can vectorize without private-temporary analysis.
! KIND=LDP precision is preserved: DBLE(JFNK_RSCALE(I)) casts the
! LDP scalar to DOUBLE PRECISION for the dot product accumulation,
! matching the original SA*SB product.
!
!$OMP SIMD REDUCTION(+:S)
	  DO I = 1, N
	    S = S + DBLE(JFNK_RSCALE(I)) * A(I) *
	1        DBLE(JFNK_RSCALE(I)) * B(I)
	  END DO
!$OMP END SIMD
	ELSE
!$OMP SIMD REDUCTION(+:S)
	  DO I = 1, N
	    S = S + A(I) * B(I)
	  END DO
!$OMP END SIMD
	END IF
	JFNK_DOT_V1 = S
!
	END FUNCTION JFNK_DOT_V1


!
! JFNK_GMRES_AUDIT_V1: Controlled single-depth toy GMRES experiment.
!
! Extracts the BA block at a specified depth K, constructs a toy
! problem A=BA_k, M=BA_k+ridge*I, and runs a self-contained GMRES
! with four different RSCALE choices to isolate the effect of
! residual weighting on GMRES convergence.
!
! This is a DIAGNOSTIC ONLY — it does not modify any production
! solver state.  Called once after the preconditioner is factored
! and RSCALE is computed.
!
! editted by Sean Brennan (10-Oct-2026)
!
	SUBROUTINE JFNK_GMRES_AUDIT_V1(TEST_K)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_BA_RIDGE_FAC
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: TEST_K
!
	INTEGER, PARAMETER :: M_TOY = 40
	INTEGER, PARAMETER :: MAXITS = 100
!
	REAL(KIND=LDP), ALLOCATABLE :: A_MAT(:,:), M_LU(:,:), M_DIAG(:)
	REAL(KIND=LDP), ALLOCATABLE :: B_VEC(:), X_VEC(:), V_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: RS(:), H(:,:), HH(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: V_BASIS(:,:), Y_LS(:), P_LS(:)
	REAL(KIND=LDP), ALLOCATABLE :: WORK(:)
	INTEGER, ALLOCATABLE :: IPIV_TOY(:), PIV_LS(:)
!
	REAL(KIND=LDP) BETA, RES, RES0, T1, T2
	REAL(KIND=LDP) RIDGE_TOY, MAX_DIAG
	INTEGER I, J, K, JJ, IOS, RANK, INFO_LU, LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=1) NO_TRANS
	LOGICAL DONE
!
! RSCALE choices:
!   1 = Euclidean (R=1)
!   2 = Current (R=1/sqrt(|STEQ|))
!   3 = Proposed (R=sqrt(|diag|/(|diag|+ridge)))
!   4 = Inverse diagnostic (R=sqrt((|diag|+ridge)/|diag|))
!
	REAL(KIND=LDP), ALLOCATABLE :: R_EUCL(:), R_CURR(:), R_PROP(:), R_INV(:)
	REAL(KIND=LDP), ALLOCATABLE :: R_USE(:)
	CHARACTER(LEN=20) :: R_NAME
	INTEGER R_CASE
	LOGICAL, SAVE :: AUDIT_DONE = .FALSE.
	REAL(KIND=LDP) EIG_VAL, EIG_MAX, EIG_MIN
!
	LUER = ERROR_LU()
	NO_TRANS = 'N'
!
! Only run once
!
	IF(AUDIT_DONE) RETURN
	AUDIT_DONE = .TRUE.
!
	IF(.NOT. JFNK_LU_READY) RETURN
	IF(TEST_K .LT. 1 .OR. TEST_K .GT. JFNK_ND) RETURN
!
	WRITE(LUER,'(/,A)')' ==================== GMRES AUDIT (single-depth toy problem) ==================='
	WRITE(LUER,'(A,I4,A,I6)') '  Test depth K=', TEST_K, '  NT=', JFNK_NT
	FLUSH(LUER)
!
! Extract BA block at TEST_K (without ridge — the true Jacobian)
!
	ALLOCATE(A_MAT(JFNK_NT, JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(M_LU(JFNK_NT, JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(M_DIAG(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(B_VEC(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(X_VEC(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(V_VEC(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(R_EUCL(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(R_CURR(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(R_PROP(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(R_INV(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(R_USE(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(IPIV_TOY(JFNK_NT), STAT=IOS)
	IF(IOS .NE. 0) RETURN
!
! A = BA_k (true Jacobian, no ridge)
! M = BA_k + ridge*I (preconditioner)
!
! JFNK_BA_A includes the ridge on its diagonal.  Subtract it to
! get the true Jacobian A.
!
	DO J = 1, JFNK_NT
	  DO I = 1, JFNK_NT
	    A_MAT(I, J) = JFNK_BA_A(I, J, TEST_K)
	  END DO
	  A_MAT(J, J) = A_MAT(J, J) - JFNK_BA_RIDGE_ADD(J, TEST_K)
	END DO
!
! Build M = A + ridge*I (same ridge as production)
!
	RIDGE_TOY = JFNK_BA_RIDGE_ADD(1, TEST_K)  ! global ridge
	MAX_DIAG = 0.0_LDP
	DO I = 1, JFNK_NT
	  M_DIAG(I) = ABS(A_MAT(I, I))
	  MAX_DIAG = MAX(MAX_DIAG, M_DIAG(I))
	END DO
	WRITE(LUER,'(A,ES12.4)') '  ridge=', RIDGE_TOY
	WRITE(LUER,'(A,ES12.4)') '  max|diag|=', MAX_DIAG
	WRITE(LUER,'(A,ES12.4)') '  ridge/max|diag|=', RIDGE_TOY/MAX_DIAG
!
	M_LU = A_MAT
	DO I = 1, JFNK_NT
	  M_LU(I, I) = M_LU(I, I) + RIDGE_TOY
	END DO
!
! Factor M with LU (no equilibration — toy problem)
!
	CALL DGETRF(JFNK_NT, JFNK_NT, M_LU, JFNK_NT, IPIV_TOY, INFO_LU)
	IF(INFO_LU .NE. 0)THEN
	  WRITE(LUER,'(A,I4)') '  LU factorization failed, info=', INFO_LU
	  GO TO 100
	END IF
!
! Build RSCALE choices
!
	R_EUCL = 1.0_LDP
	DO I = 1, JFNK_NT
	  IF(JFNK_RSCALE_READY)THEN
	    R_CURR(I) = JFNK_RSCALE((TEST_K-1)*JFNK_NT + I)
	  ELSE
	    R_CURR(I) = 1.0_LDP
	  END IF
	  IF(M_DIAG(I) .GT. 1.0E-30_LDP .AND. RIDGE_TOY .GT. 0.0_LDP)THEN
	    R_PROP(I) = SQRT(M_DIAG(I) / (M_DIAG(I) + RIDGE_TOY))
	    R_INV(I)  = SQRT((M_DIAG(I) + RIDGE_TOY) / M_DIAG(I))
	  ELSE
	    R_PROP(I) = 1.0E-15_LDP
	    R_INV(I)  = 1.0E15_LDP
	  END IF
	END DO
!
! Print RSCALE summaries
!
	T1 = 0.0_LDP; T2 = 0.0_LDP
	DO I = 1, JFNK_NT
	  T1 = MIN(T1, R_CURR(I)); T2 = MAX(T2, R_CURR(I))
	END DO
	WRITE(LUER,'(A,ES10.3,A,ES10.3)') '  R_CURR: min=', T1, ' max=', T2
	T1 = 0.0_LDP; T2 = 0.0_LDP
	DO I = 1, JFNK_NT
	  T1 = MIN(T1, R_PROP(I)); T2 = MAX(T2, R_PROP(I))
	END DO
	WRITE(LUER,'(A,ES10.3,A,ES10.3)') '  R_PROP: min=', T1, ' max=', T2
	T1 = 0.0_LDP; T2 = 0.0_LDP
	DO I = 1, JFNK_NT
	  T1 = MIN(T1, R_INV(I)); T2 = MAX(T2, R_INV(I))
	END DO
	WRITE(LUER,'(A,ES10.3,A,ES10.3)') '  R_INV:  min=', T1, ' max=', T2
	FLUSH(LUER)
!
! Direct test: M^{-1} * A * v - v for a random v
! If M ≈ A, this should be ≈ 0
!
	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(V_VEC)
	V_VEC = V_VEC * 2.0_LDP - 1.0_LDP  ! [-1, 1]
!
! Compute A*v
!
	B_VEC = MATMUL(A_MAT, V_VEC)
!
! Compute M^{-1} * (A*v)
!
	CALL DGETRS(NO_TRANS, JFNK_NT, 1, M_LU, JFNK_NT, IPIV_TOY, B_VEC,
	1   JFNK_NT, INFO_LU)
!
! B_VEC now holds M^{-1}*A*v.  Compare to v.
!
	T1 = 0.0_LDP; T2 = 0.0_LDP; RES = 0.0_LDP
	DO I = 1, JFNK_NT
	  T1 = MAX(T1, ABS(B_VEC(I) - V_VEC(I)))
	  T2 = MAX(T2, ABS(V_VEC(I)))
	  RES = RES + (B_VEC(I) - V_VEC(I))**2
	END DO
	WRITE(LUER,'(A)') '  Direct M^{-1}*A*v - v test:'
	WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)') '    ||M^{-1}Av-v||_inf=',
	1   T1, ' ||v||_inf=', T2, ' rel=', T1/MAX(T2,1.0E-30_LDP)
	WRITE(LUER,'(A,ES12.4)') '    ||M^{-1}Av-v||_2  =', SQRT(RES)
!
! Eigenvalue spread estimate (diagonal approximation)
!
	EIG_MAX = 1.0_LDP; EIG_MIN = 1.0E30_LDP
	DO I = 1, JFNK_NT
	  IF(M_DIAG(I) .GT. 1.0E-30_LDP)THEN
	    EIG_VAL = M_DIAG(I) / (M_DIAG(I) + RIDGE_TOY)
	    EIG_MAX = MAX(EIG_MAX, EIG_VAL)
	    EIG_MIN = MIN(EIG_MIN, EIG_VAL)
	  END IF
	END DO
	WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1   '  Eigenvalue spread: max=', EIG_MAX, ' min=', EIG_MIN,
	1   ' ratio=', EIG_MAX/MAX(EIG_MIN,1.0E-30_LDP)
	FLUSH(LUER)
!
! Run toy GMRES for each RSCALE choice
!
	ALLOCATE(H(M_TOY+1, M_TOY), STAT=IOS)
	ALLOCATE(V_BASIS(JFNK_NT, M_TOY+1), STAT=IOS)
	ALLOCATE(Y_LS(M_TOY+1), STAT=IOS)
	ALLOCATE(P_LS(M_TOY+1), STAT=IOS)
	ALLOCATE(HH(M_TOY+1, M_TOY), STAT=IOS)
	ALLOCATE(PIV_LS(M_TOY+1), STAT=IOS)
	ALLOCATE(WORK(6*M_TOY+16), STAT=IOS)
!
	DO R_CASE = 1, 4
	  IF(R_CASE .EQ. 1)THEN
	    R_USE = R_EUCL; R_NAME = 'Euclidean (R=1)'
	  ELSE IF(R_CASE .EQ. 2)THEN
	    R_USE = R_CURR; R_NAME = 'Current (1/sqrt|STEQ|)'
	  ELSE IF(R_CASE .EQ. 3)THEN
	    R_USE = R_PROP; R_NAME = 'Proposed (sqrt(d/(d+r)))'
	  ELSE
	    R_USE = R_INV; R_NAME = 'Inverse (sqrt((d+r)/d))'
	  END IF
!
! Toy GMRES: solve A*x = b with M preconditioner and R_USE weighting
! b = random, x_0 = 0
!
	  CALL RANDOM_NUMBER(B_VEC)
	  B_VEC = B_VEC * 2.0_LDP - 1.0_LDP
	  X_VEC = 0.0_LDP
!
! Initial residual: r = b - A*x = b
! Weighted norm: ||r||_w = sqrt(sum(R_USE^2 * r^2))
!
	  BETA = 0.0_LDP
	  DO I = 1, JFNK_NT
	    BETA = BETA + (R_USE(I) * B_VEC(I))**2
	  END DO
	  BETA = SQRT(BETA)
	  RES0 = BETA
!
	  V_BASIS = 0.0_LDP
	  DO I = 1, JFNK_NT
	    V_BASIS(I, 1) = B_VEC(I) / BETA
	  END DO
!
	  H = 0.0_LDP
	  RES = BETA
	  DONE = .FALSE.
!
	  WRITE(LUER,'(A,A)') '  --- RSCALE: ', TRIM(R_NAME)
	  WRITE(LUER,'(A,ES12.4)') '    Initial ||r||_w = ', RES0
!
	  DO J = 1, MIN(MAXITS, M_TOY)
!     w = A * M^{-1} * v_j  (right preconditioning)
	    DO I = 1, JFNK_NT
	      V_VEC(I) = V_BASIS(I, J)
	    END DO
!     M^{-1} * v_j
	    CALL DGETRS(NO_TRANS, JFNK_NT, 1, M_LU, JFNK_NT, IPIV_TOY,
	1       V_VEC, JFNK_NT, INFO_LU)
!     A * (M^{-1} * v_j)
	    B_VEC = MATMUL(A_MAT, V_VEC)
!
!     Modified Gram-Schmidt with weighted inner product
	    DO I = 1, J
	      H(I, J) = 0.0_LDP
	      DO K = 1, JFNK_NT
	        H(I, J) = H(I, J) + R_USE(K)**2 * B_VEC(K) * V_BASIS(K, I)
	      END DO
	      DO K = 1, JFNK_NT
	        B_VEC(K) = B_VEC(K) - H(I, J) * V_BASIS(K, I)
	      END DO
	    END DO
!     Re-orthogonalization
	    DO I = 1, J
	      T1 = 0.0_LDP
	      DO K = 1, JFNK_NT
	        T1 = T1 + R_USE(K)**2 * B_VEC(K) * V_BASIS(K, I)
	      END DO
	      H(I, J) = H(I, J) + T1
	      DO K = 1, JFNK_NT
	        B_VEC(K) = B_VEC(K) - T1 * V_BASIS(K, I)
	      END DO
	    END DO
!     h(j+1,j) = ||w||_w
	    H(J+1, J) = 0.0_LDP
	    DO K = 1, JFNK_NT
	      H(J+1, J) = H(J+1, J) + (R_USE(K) * B_VEC(K))**2
	    END DO
	    H(J+1, J) = SQRT(H(J+1, J))
!
	    IF(H(J+1, J) .LT. 1.0E-14_LDP * BETA)THEN
	      V_BASIS(:, J+1) = 0.0_LDP
	    ELSE
	      DO K = 1, JFNK_NT
	        V_BASIS(K, J+1) = B_VEC(K) / H(J+1, J)
	      END DO
	    END IF
!
!     Least-squares: minimize ||beta*e1 - H*y||_2
	    P_LS(1) = BETA
	    P_LS(2:J+1) = 0.0_LDP
	    HH(1:J+1, 1:J) = H(1:J+1, 1:J)
	    PIV_LS = 0
	    CALL DGELSY(J+1, J, 1, HH, M_TOY+1, P_LS, M_TOY+1, PIV_LS,
	1       1.0E-14_LDP, RANK, WORK, 6*M_TOY+16, INFO_LU)
	    Y_LS = P_LS
!
!     Residual: ||beta*e1 - H*y||_2
	    P_LS(1:J+1) = -MATMUL(H(1:J+1, 1:J), Y_LS(1:J))
	    P_LS(1) = P_LS(1) + BETA
	    RES = SQRT(DOT_PRODUCT(P_LS(1:J+1), P_LS(1:J+1)))
!
	    WRITE(LUER,'(A,I3,A,ES12.4,A,ES8.2)')
	1       '    it=', J, ' res=', RES, ' rel=', RES/RES0
	    FLUSH(LUER)
!
	    IF(RES .LT. 1.0E-10_LDP * RES0 .OR. J .EQ. M_TOY) EXIT
	  END DO
!
	  WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1       '  Final: res=', RES, ' reduction=', RES/RES0
	  FLUSH(LUER)
	END DO  ! R_CASE
!
	DEALLOCATE(H, V_BASIS, Y_LS, P_LS, HH, PIV_LS, WORK)
!
100	CONTINUE
	DEALLOCATE(A_MAT, M_LU, M_DIAG, B_VEC, X_VEC, V_VEC)
	DEALLOCATE(R_EUCL, R_CURR, R_PROP, R_INV, R_USE, IPIV_TOY)
!
	WRITE(LUER,'(/,A)')' =============================== END GMRES AUDIT ==============================='
	FLUSH(LUER)
!
	END SUBROUTINE JFNK_GMRES_AUDIT_V1


!
! JFNK_RSCALE_AUDIT_V1: Full multi-depth RSCALE comparison.
!
! Runs the production GMRES three times at the same saved CMFGEN
! state with different RSCALE choices:
!   1 = Current (1/sqrt(|STEQ|))
!   2 = Euclidean (R=1)
!   3 = Proposed (sqrt(|diag|/(|diag|+ridge)))
!
! For each choice, records:
!   - Weighted GMRES residual (production convergence metric)
!   - Unweighted ||F(-s)||_inf (physical nonlinear residual)
!   - ||s||_w, ||s||_inf, max population correction
!   - ZMAX/clipping events
!
! Then reports the 20 smallest-RSCALE components for the proposed
! case and investigates depth 78 vs neighbors.
!
! DIAGNOSTIC ONLY — does not modify production solver state.
! The production JFNK solve proceeds normally after this audit.
!
! editted by Sean Brennan (10-Oct-2026)
!
	SUBROUTINE JFNK_RSCALE_AUDIT_V1()
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_M_GMRES, JFNK_MAX_GMRES,
	1   JFNK_LEFT_PRECOND
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_x, jfnk_fx,
	1   jfnk_use_left_pc, jfnk_n_zmax_hit, jfnk_total_jvp_evals
! Fix #7: JFNK_N_PRECOND_APPLIES is in JFNK_CMFGEN_INT_MOD (USEd above)
	IMPLICIT NONE
!
	LOGICAL, SAVE :: AUDIT_DONE = .FALSE.
	INTEGER LUER, ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER I, J, K, R_CASE, N, M_TOY, ITS_TOY, INFO_TOY
	INTEGER II, JJ, IOS
	REAL(KIND=LDP) RES_W, PRED_RES, DEL_TOY
	REAL(KIND=LDP) F0_W, F0_INF, F0_L2
	REAL(KIND=LDP) S_W, S_INF, S_MAXPOP
	REAL(KIND=LDP) FS_W, FS_INF, FS_L2
	REAL(KIND=LDP) T1, T2
	LOGICAL SAVE_LEFT_PC
	REAL(KIND=LDP), ALLOCATABLE :: S_VEC(:), F_AFTER(:)
	REAL(KIND=LDP), ALLOCATABLE :: H_WORK(:), V_WORK(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RSCALE_SAVE(:)
	REAL(KIND=LDP), ALLOCATABLE :: R_PROP(:)
!
	EXTERNAL JFNK_JVP_V1, JFNK_PRECOND_V1
	EXTERNAL gmres_m_v1
	DOUBLE PRECISION JFNK_DOT_V1
	EXTERNAL JFNK_DOT_V1
	EXTERNAL JFNK_RESIDUAL_V1, JFNK_MASK_FIXED_UNKNOWNS
!
	CHARACTER(LEN=30) :: R_NAME
!
	IF(AUDIT_DONE) RETURN
	AUDIT_DONE = .TRUE.
!
	IF(.NOT. JFNK_LU_READY) RETURN
	IF(.NOT. JFNK_RSCALE_READY) RETURN
!
	LUER = ERROR_LU()
	N = JFNK_NT * JFNK_ND
	M_TOY = JFNK_M_GMRES
	IF(M_TOY .LT. 1) M_TOY = 10
!
	ALLOCATE(S_VEC(N), F_AFTER(N), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(H_WORK((M_TOY+1)*M_TOY), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(V_WORK(N, M_TOY+1), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(RSCALE_SAVE(N), STAT=IOS)
	IF(IOS .NE. 0) RETURN
	ALLOCATE(R_PROP(N), STAT=IOS)
	IF(IOS .NE. 0) RETURN
!
! Save current RSCALE and solver state
!
	RSCALE_SAVE(1:N) = JFNK_RSCALE(1:N)
!
! Compute proposed RSCALE: R_i = sqrt(|D_i| / (|D_i| + ridge))
!
	DO K = 1, JFNK_ND
	  DO I = 1, JFNK_NT
	    T1 = ABS(JFNK_BA_A(I, I, K) - JFNK_BA_RIDGE_ADD(I, K))
	    T2 = ABS(JFNK_BA_RIDGE_ADD(I, K))
	    IF(T1 .GT. 1.0E-30_LDP .AND. T2 .GT. 0.0_LDP)THEN
	      R_PROP((K-1)*JFNK_NT + I) = SQRT(T1 / (T1 + T2))
	    ELSE
	      R_PROP((K-1)*JFNK_NT + I) = 1.0E-15_LDP
	    END IF
	  END DO
	END DO
!
! Compute initial residual norms (unweighted)
!
	F0_INF = 0.0_LDP
	F0_L2 = 0.0_LDP
	DO I = 1, N
	  F0_INF = MAX(F0_INF, ABS(jfnk_fx(I)))
	  F0_L2 = F0_L2 + jfnk_fx(I)**2
	END DO
	F0_L2 = SQRT(F0_L2)
	F0_W = SQRT(JFNK_DOT_V1(N, jfnk_fx, jfnk_fx))
!
	WRITE(LUER,'(/,A)')
	1   ' =========================== MULTI-DEPTH RSCALE AUDIT =========================='
	WRITE(LUER,'(A,I6,A,I4,A,I8)')
	1   '  NT=', JFNK_NT, ' ND=', JFNK_ND, ' N=', N
	WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1   '  Initial ||F(0)||_w  =', F0_W,
	1   ' ||F(0)||_inf=', F0_INF, ' ||F(0)||_2  =', F0_L2
	FLUSH(LUER)
!
! Save left_pc setting and force right preconditioning for the RSCALE
! comparison below.  T1/T2 are reused as RSCALE scratch inside the
! loop, so the flag must live in its own logical -- the old real
! save was clobbered by the min() pass and the restore read the
! RSCALE minimum instead of the saved state.
!
	SAVE_LEFT_PC = jfnk_use_left_pc
	jfnk_use_left_pc = .FALSE.
!
	DO R_CASE = 1, 3
	  IF(R_CASE .EQ. 1)THEN
	    R_NAME = 'Current (1/sqrt|STEQ|)'
	    JFNK_RSCALE(1:N) = RSCALE_SAVE(1:N)
	  ELSE IF(R_CASE .EQ. 2)THEN
	    R_NAME = 'Euclidean (R=1)'
	    JFNK_RSCALE(1:N) = 1.0_LDP
	  ELSE
	    R_NAME = 'Proposed (sqrt(d/(d+r)))'
	    JFNK_RSCALE(1:N) = R_PROP(1:N)
	  END IF
	  JFNK_RSCALE_READY = .TRUE.
!
! Recompute F0_W with this RSCALE
!
	  F0_W = SQRT(JFNK_DOT_V1(N, jfnk_fx, jfnk_fx))
!
	  WRITE(LUER,'(A)')
	1   ' '
	  WRITE(LUER,'(A,A)') '--- RSCALE: ', TRIM(R_NAME)
!
! RSCALE summary
!
	  T1 = JFNK_RSCALE(1); T2 = JFNK_RSCALE(1)
	  DO I = 1, N
	    T1 = MIN(T1, JFNK_RSCALE(I))
	    T2 = MAX(T2, JFNK_RSCALE(I))
	  END DO
	  WRITE(LUER,'(A,ES10.3,A,ES10.3)')
	1     '  RSCALE: min=', T1, ' max=', T2
	  WRITE(LUER,'(A,ES12.4)')
	1     '  Initial ||F(0)||_w =', F0_W
	  FLUSH(LUER)
!
! Run GMRES: solve J*s = F(0), x=0
! tol=0 forces full iteration history; del=0 disables trust region
!
	  S_VEC = 0.0_LDP
	  ITS_TOY = JFNK_MAX_GMRES
	  IF(ITS_TOY .LT. M_TOY) ITS_TOY = M_TOY
	  INFO_TOY = 1  ! verbose
	  DEL_TOY = 0.0_LDP
	  jfnk_n_zmax_hit = 0
!
	  CALL gmres_m_v1(M_TOY, N, S_VEC, jfnk_fx,
	1     JFNK_JVP_V1, JFNK_PRECOND_V1, JFNK_DOT_V1,
	1     H_WORK, V_WORK, RES_W, DEL_TOY, PRED_RES,
	1     ITS_TOY, INFO_TOY, LUER)
	  FLUSH(LUER)
!
! Record GMRES results
!
	  WRITE(LUER,'(A,ES12.4,A,I5,A,I2)')
	1     '  GMRES final: res_w=', RES_W,
	1     ' its=', ITS_TOY, ' info=', INFO_TOY
	  WRITE(LUER,'(A,I6)')
	1     '  ZMAX hits:', jfnk_n_zmax_hit
!
! Newton correction is -s (GMRES solves J*s=F(0), Newton uses delta=-s)
! Compute correction norms
!
	  S_W = SQRT(JFNK_DOT_V1(N, S_VEC, S_VEC))
	  S_INF = 0.0_LDP
	  S_MAXPOP = 0.0_LDP
	  DO I = 1, N
	    S_INF = MAX(S_INF, ABS(S_VEC(I)))
	    S_MAXPOP = MAX(S_MAXPOP, ABS(S_VEC(I)))
	  END DO
	  WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1     '  Correction: ||s||_w=', S_W,
	1     ' ||s||_inf=', S_INF,
	1     ' max|s|=', S_MAXPOP
!
! Compute nonlinear residual F(-s) by calling JFNK_RESIDUAL_V1
! The correction applied is -s (Newton step: x_new = x_old - s)
!
	  S_VEC = -S_VEC
	  CALL JFNK_RESIDUAL_V1(N, S_VEC, F_AFTER)
	  CALL JFNK_MASK_FIXED_UNKNOWNS(N, F_AFTER)
	  S_VEC = -S_VEC  ! restore
!
	  FS_INF = 0.0_LDP
	  FS_L2 = 0.0_LDP
	  DO I = 1, N
	    FS_INF = MAX(FS_INF, ABS(F_AFTER(I)))
	    FS_L2 = FS_L2 + F_AFTER(I)**2
	  END DO
	  FS_L2 = SQRT(FS_L2)
	  FS_W = SQRT(JFNK_DOT_V1(N, F_AFTER, F_AFTER))
!
	  WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1     '  F(-s): ||F||_w  =', FS_W,
	1     ' ||F||_inf=', FS_INF, ' ||F||_2  =', FS_L2
	  WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1     '  Reduction: weighted  =', FS_W/MAX(F0_W,1.0E-30_LDP),
	1     ' unweighted=', FS_INF/MAX(F0_INF,1.0E-30_LDP)
	  FLUSH(LUER)
	END DO  ! R_CASE
!
! Restore RSCALE and solver state
!
	JFNK_RSCALE(1:N) = RSCALE_SAVE(1:N)
	JFNK_RSCALE_READY = .TRUE.
	jfnk_use_left_pc = SAVE_LEFT_PC
!
! Report 20 smallest-RSCALE components for proposed case
!
	WRITE(LUER,'(/,A)')
	1   ' ------------------- 20 Smallest-RSCALE Components (Proposed) ------------------'
	WRITE(LUER,'(A)')
	1   '  idx  depth  sp  |D|          ridge        RSCALE       |F(0)|        |s|'
!
	DO II = 1, 20
	  T1 = 1.0E30_LDP
	  K = 0
	  DO JJ = 1, N
	    IF(R_PROP(JJ) .GT. 0.0_LDP .AND. R_PROP(JJ) .LT. T1)THEN
	      T1 = R_PROP(JJ)
	      K = JJ
	    END IF
	  END DO
	  IF(K .GT. 0)THEN
	    R_PROP(K) = 1.0E30_LDP  ! mark as visited
	    J = (K-1)/JFNK_NT + 1  ! depth
	    I = MOD(K-1, JFNK_NT) + 1  ! species
	    T1 = ABS(JFNK_BA_A(I, I, J) - JFNK_BA_RIDGE_ADD(I, J))
	    T2 = ABS(JFNK_BA_RIDGE_ADD(I, J))
	    WRITE(LUER,'(I6,I5,I4,5ES13.4)')
	1       K, J, I, T1, T2, SQRT(T1/(T1+T2)),
	1       ABS(jfnk_fx(K)), ABS(S_VEC(K))
	  END IF
	END DO
	FLUSH(LUER)
!
! Investigate depth 78 vs neighbors
!
	WRITE(LUER,'(/,A)')
	1   ' --------------------- Depth 78 vs Neighbors: BA Structure ---------------------'
	WRITE(LUER,'(A)')
	1   '  K   max|diag|    min|diag|    diag_range   n_nearzero   max_ridge/diag'
!
	DO K = 76, 80
	  IF(K .LT. 1 .OR. K .GT. JFNK_ND) CYCLE
	  T1 = 0.0_LDP  ! max|diag|
	  T2 = 1.0E30_LDP  ! min|diag|
	  I = 0  ! n_nearzero
	  DEL_TOY = 0.0_LDP  ! max ridge/diag
	  DO J = 1, JFNK_NT
	    F0_INF = ABS(JFNK_BA_A(J, J, K) - JFNK_BA_RIDGE_ADD(J, K))
	    IF(F0_INF .GT. 1.0E-30_LDP)THEN
	      T1 = MAX(T1, F0_INF)
	      T2 = MIN(T2, F0_INF)
	    ELSE
	      I = I + 1
	    END IF
	    F0_L2 = ABS(JFNK_BA_RIDGE_ADD(J, K))
	    IF(F0_INF .GT. 1.0E-30_LDP)THEN
	      DEL_TOY = MAX(DEL_TOY, F0_L2 / F0_INF)
	    END IF
	  END DO
	  WRITE(LUER,'(I4,5ES13.4,I6,ES13.4)')
	1     K, T1, T2, T1/MAX(T2,1.0E-30_LDP), I, DEL_TOY
	END DO
	FLUSH(LUER)
!
	DEALLOCATE(S_VEC, F_AFTER, H_WORK, V_WORK, RSCALE_SAVE, R_PROP)
!
	WRITE(LUER,'(/,A)')
	1   ' ========================= END MULTI-DEPTH RSCALE AUDIT ========================'
	FLUSH(LUER)
!
	END SUBROUTINE JFNK_RSCALE_AUDIT_V1


!
! JFNK_DIAGNOSTIC_V1: Called at the end of each Newton iteration.
!
! editted by Sean Brennan (V4, 2026-09)
!
	SUBROUTINE JFNK_DIAGNOSTIC_V1()
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_PC_REFRESH_FRAC,
	1                              JFNK_DEEPDIAG
	IMPLICIT NONE
!
! V4 audit (2026-10): Run preconditioner/RSCALE/BA audit once
! at the first Newton iteration (nits=0).
!
! JFNK_DEEPDIAG (03-Oct-2026): COMP_DIAGNOSTIC and TEMP_COL_AUDIT
! spend ~12 residual evaluations -- a full formal sweep each in
! coupled mode -- so they are gated on [JFNK_DEEPDIAG] (default F).
! The PRECOND_AUDIT only applies the cached LU and is kept
! unconditional.
!
	IF(jfnk_nits .EQ. 0) THEN
	  CALL JFNK_PRECOND_AUDIT_V1(JFNK_NT * JFNK_ND,
	1     DBLE(jfnk_fx), JFNK_LUER)
	  IF(JFNK_DEEPDIAG)THEN
	    CALL JFNK_COMP_DIAGNOSTIC_V1(JFNK_NT * JFNK_ND,
	1       DBLE(jfnk_fx), JFNK_LUER)
	    CALL JFNK_TEMP_COL_AUDIT_V1(JFNK_LUER)
	  END IF
	END IF
!
	WRITE(JFNK_LUER,'(A)')'  JFNK diagnostics:'
	WRITE(JFNK_LUER,'(A,I6)')'    Newton-level F evaluations:',
	1     jfnk_total_f_evals
	WRITE(JFNK_LUER,'(A,I6)')'    Residual evals incl JVPs: ',
	1     JFNK_N_RESIDUAL_EVALS
	WRITE(JFNK_LUER,'(A,I6)')'    Total JVP evaluations:',
	1     jfnk_total_jvp_evals
	WRITE(JFNK_LUER,'(A,I6)')'    Preconditioner calls: ',
	1     JFNK_N_PRECOND_APPLIES
	WRITE(JFNK_LUER,'(A,I6)')'    Preconditioner fails: ',
	1     jfnk_precond_fails
	WRITE(JFNK_LUER,'(A,I6)')'    NaN rejected steps:   ',
	1     jfnk_nan_rejects
	WRITE(JFNK_LUER,'(A,I6)')'    Residual NaN evals:   ',
	1     JFNK_N_RESIDUAL_NAN
	WRITE(JFNK_LUER,'(A,I6)')'    LU factorizations:    ',
	1     JFNK_N_LU_FACTORS
	WRITE(JFNK_LUER,'(A,I6)')'    Precond refreshes:    ',
	1     JFNK_N_PRECOND_REFRESH
	WRITE(JFNK_LUER,'(A)')'  Stabilization counters (review §6):'
	WRITE(JFNK_LUER,'(A,I6)')'    Step scaling (|s|>10):  ',
	1     jfnk_n_step_scaled
	WRITE(JFNK_LUER,'(A,I6)')'    Precond ZMAX hits:       ',
	1     jfnk_n_zmax_hit
	WRITE(JFNK_LUER,'(A,I6)')'    Reg-row projections:     ',
	1     jfnk_n_reg_proj
	WRITE(JFNK_LUER,'(A,I6)')'    Cauchy fallback used:    ',
	1     jfnk_n_cauchy_used
	WRITE(JFNK_LUER,'(A,I6)')'    Line search used:        ',
	1     jfnk_n_ls_used
	WRITE(JFNK_LUER,'(A,I6)')'    Trust-region rejections: ',
	1     jfnk_n_tr_rejected
	WRITE(JFNK_LUER,'(A,I6)')'    Descent retries:         ',
	1     jfnk_n_descent_retry
!
! Issue 3: if last GMRES used most of its Krylov budget, bump the
! ridge and re-factor the cached DIAG LU (preconditioner lagging).
! Relaxation (audit §13, 09-Sep-2026): when GMRES used few iterations
! (well below the refresh threshold), relax the ridge back down so
! the preconditioner does not drift permanently away from the physical
! BA.  The ridge only increases when GMRES is struggling, and decreases
! when GMRES is converging easily.
!
	IF(jfnk_need_precond_refresh .NE. 0)THEN
	  WRITE(JFNK_LUER,'(A,I4,A,F5.2)')
	1    '  JFNK: REBUILDING PRECONDITIONER (GMRES its=',
	1    jfnk_last_gmres_its, ' frac=', JFNK_PC_REFRESH_FRAC, ')'
	  CALL JFNK_FACTOR_PRECOND_V1()
	  JFNK_N_PRECOND_REFRESH = JFNK_N_PRECOND_REFRESH + 1
	  jfnk_need_precond_refresh = 0
	END IF
!
	END SUBROUTINE JFNK_DIAGNOSTIC_V1


!
! JFNK_PRECOND_AUDIT_V1: One-shot preconditioner/RSCALE/BA audit.
!
! Reports:
!   A. RSCALE range (min, max, mean, dynamic range)
!   B. BA block scales per species (|diag|, ridge, ridge/|diag|, |offdiag|)
!   C. Preconditioner accuracy: ||M*(M^-1 v) - v|| / ||v|| for test vectors
!   D. Preconditioned operator: ||J*(M^-1 v)|| / ||v|| (how much J*M^-1 amplifies)
!
! Called once at the first JFNK Newton iteration.
!
! V4 audit, 2026-10.
!
	SUBROUTINE JFNK_PRECOND_AUDIT_V1(N, F_BASE, LUER)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_BA_RIDGE_FAC
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: N, LUER
	DOUBLE PRECISION, INTENT(IN) :: F_BASE(N)
!
	REAL(KIND=LDP) R_MIN, R_MAX, R_MEAN, R_DRANGE
	REAL(KIND=LDP) D_MIN, D_MAX, D_MEAN, DG, OFFDG, OFFMAX
	REAL(KIND=LDP) RIDGE_I, RATIO, RATIO_MAX, RATIO_MIN
	INTEGER I, J, K, KK, II, N_SP_REPORT
	LOGICAL, SAVE :: AUDIT_DONE = .FALSE.
!
! Workspace for preconditioner accuracy test
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: V_TEST(:), Z_TEST(:), R_TEST(:)
	INTEGER IOS_A
	REAL(KIND=LDP) V_NORM, R_NORM, REL_ERR, Z_NORM
	REAL(KIND=LDP) JZ_NORM, JM_RATIO
	REAL(KIND=LDP) R_NORM_P, V_NORM_T
	INTEGER JR
!
	IF(AUDIT_DONE) RETURN
	AUDIT_DONE = .TRUE.
!
	IF(.NOT. JFNK_LU_READY) THEN
	  WRITE(LUER,'(A)') '  PRECOND AUDIT: LU not ready, skipping'
	  RETURN
	END IF
!
	IF(.NOT. ALLOCATED(V_TEST)) THEN
	  ALLOCATE(V_TEST(N), Z_TEST(N), R_TEST(N), STAT=IOS_A)
	  IF(IOS_A .NE. 0) THEN
	    WRITE(LUER,'(A)') '  PRECOND AUDIT: alloc failed, skipping'
	    RETURN
	  END IF
	END IF
!
	WRITE(LUER,'(/,A)')' ====================== PRECONDITIONER / RSCALE / BA AUDIT ====================='
!
! ── A. RSCALE range ──────────────────────────────────────────────
!
	IF(JFNK_RSCALE_READY) THEN
	  R_MIN = JFNK_RSCALE(1)
	  R_MAX = JFNK_RSCALE(1)
	  R_MEAN = 0.0_LDP
	  DO I = 1, N
	    R_MIN = MIN(R_MIN, JFNK_RSCALE(I))
	    R_MAX = MAX(R_MAX, JFNK_RSCALE(I))
	    R_MEAN = R_MEAN + JFNK_RSCALE(I)
	  END DO
	  R_MEAN = R_MEAN / DBLE(N)
	  R_DRANGE = R_MAX / MAX(R_MIN, 1.0E-30_LDP)
	  WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)')
	1    '  RSCALE: min=', R_MIN, ' max=', R_MAX,
	1    ' mean=', R_MEAN, ' dynamic_range=', R_DRANGE
	ELSE
	  WRITE(LUER,'(A)') '  RSCALE: not ready (uniform weighting)'
	END IF
!
! ── B. BA block scales per species (first depth only as representative) ──
!
	K = 1
	D_MIN = ABS(JFNK_BA_DIAG(1, K))
	D_MAX = ABS(JFNK_BA_DIAG(1, K))
	D_MEAN = 0.0_LDP
	RATIO_MAX = 0.0_LDP
	RATIO_MIN = 1.0E+30_LDP
	DO J = 1, JFNK_NT
	  DG = ABS(JFNK_BA_DIAG(J, K))
	  D_MIN = MIN(D_MIN, DG)
	  D_MAX = MAX(D_MAX, DG)
	  D_MEAN = D_MEAN + DG
	  RIDGE_I = ABS(JFNK_BA_RIDGE_ADD(J, K))
	  IF(DG .GT. 1.0E-30_LDP) THEN
	    RATIO = RIDGE_I / DG
	    RATIO_MAX = MAX(RATIO_MAX, RATIO)
	    RATIO_MIN = MIN(RATIO_MIN, RATIO)
	  END IF
	END DO
	D_MEAN = D_MEAN / DBLE(JFNK_NT)
	WRITE(LUER,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)')
	1    '  BA diag (depth ', K, '): min=', D_MIN,
	1    ' max=', D_MAX, ' mean=', D_MEAN
	WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1    '  ridge/diag ratio: min=', RATIO_MIN,
	1    ' max=', RATIO_MAX,
	1    ' (ridge_fac=', JFNK_BA_RIDGE_FAC
!
! Report per-species diag and ridge for first 10 and last 10 species
!
	N_SP_REPORT = MIN(JFNK_NT, 10)
	WRITE(LUER,'(A)') '  Per-species |diag| and ridge/diag (first 10):'
	DO J = 1, N_SP_REPORT
	  DG = ABS(JFNK_BA_DIAG(J, K))
	  RIDGE_I = ABS(JFNK_BA_RIDGE_ADD(J, K))
	  IF(DG .GT. 1.0E-30_LDP) THEN
	    RATIO = RIDGE_I / DG
	  ELSE
	    RATIO = -1.0_LDP
	  END IF
	  WRITE(LUER,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)')
	1    '    sp', J, ': |diag|=', DG,
	1    ' ridge=', RIDGE_I, ' ridge/diag=', RATIO
	END DO
	IF(JFNK_NT .GT. 10) THEN
	  WRITE(LUER,'(A)') '  Per-species |diag| and ridge/diag (last 10):'
	  DO J = MAX(N_SP_REPORT+1, JFNK_NT-9), JFNK_NT
	    DG = ABS(JFNK_BA_DIAG(J, K))
	    RIDGE_I = ABS(JFNK_BA_RIDGE_ADD(J, K))
	    IF(DG .GT. 1.0E-30_LDP) THEN
	      RATIO = RIDGE_I / DG
	    ELSE
	      RATIO = -1.0_LDP
	    END IF
	    WRITE(LUER,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)')
	1      '    sp', J, ': |diag|=', DG,
	1      ' ridge=', RIDGE_I, ' ridge/diag=', RATIO
	  END DO
	END IF
!
! ── C. Preconditioner accuracy: ||M*(M^-1 v) - v|| / ||v|| ──────
!
! Test vector 1: the residual F itself
!
	V_TEST = F_BASE
	V_NORM = 0.0_LDP
	DO I = 1, N
	  V_NORM = MAX(V_NORM, ABS(V_TEST(I)))
	END DO
	IF(V_NORM .LT. 1.0E-30_LDP) THEN
	  WRITE(LUER,'(A)') '  PC accuracy: F is zero, skipping'
	  GO TO 100
	END IF
!
! Apply M^{-1} to get z = M^{-1} * v
!
	Z_TEST = V_TEST
	CALL JFNK_PRECOND_V1(N, Z_TEST)
!
! Compute M*z using JFNK_BA_A (the ridge-augmented BA)
!
	R_TEST = 0.0_LDP
	DO KK = 1, JFNK_ND
	  DO II = 1, JFNK_NT
	    DO J = 1, JFNK_NT
	      R_TEST((KK-1)*JFNK_NT + II) =
	1       R_TEST((KK-1)*JFNK_NT + II) +
	2       JFNK_BA_A(II, J, KK) *
	3       Z_TEST((KK-1)*JFNK_NT + J)
	    END DO
	  END DO
	END DO
!
! r = M*z - v
!
	R_NORM = 0.0_LDP
	Z_NORM = 0.0_LDP
	DO I = 1, N
	  R_TEST(I) = R_TEST(I) - V_TEST(I)
	  R_NORM = MAX(R_NORM, ABS(R_TEST(I)))
	  Z_NORM = MAX(Z_NORM, ABS(Z_TEST(I)))
	END DO
	REL_ERR = R_NORM / MAX(V_NORM, 1.0E-30_LDP)
	WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)')
	1    '  PC accuracy (v=F): ||Mz-v||/||v||=', REL_ERR,
	1    ' ||v||_inf=', V_NORM,
	1    ' ||z||_inf=', Z_NORM,
	1    ' ||Mz-v||_inf=', R_NORM
!
! Row-resolved decomposition.  The apply projects frozen rows
! (J=NT T, J=NT-1 Ne) out of z, so the full-space metric mixes the
! intended projection with real LU error.  Report the residual on
! population rows alone and the share of v carried by the frozen
! rows to separate the two.
!
	R_NORM_P = 0.0_LDP
	V_NORM_T = 0.0_LDP
	DO I = 1, N
	  JR = MOD(I-1, JFNK_NT) + 1
	  IF(JR .LT. JFNK_NT-1)THEN
	    R_NORM_P = MAX(R_NORM_P, ABS(R_TEST(I)))
	  ELSE
	    V_NORM_T = MAX(V_NORM_T, ABS(V_TEST(I)))
	  END IF
	END DO
	WRITE(LUER,'(A,ES12.4,A,ES12.4)')
	1    '  PC accuracy (v=F, pop rows): ||Mz-v||_inf=', R_NORM_P,
	1    ' |v|_inf on T/Ne rows=', V_NORM_T
!
! Test vector 2: unit vector in the first species at first depth
!
	V_TEST = 0.0_LDP
	V_TEST(1) = 1.0_LDP
	V_NORM = 1.0_LDP
	Z_TEST = V_TEST
	CALL JFNK_PRECOND_V1(N, Z_TEST)
	R_TEST = 0.0_LDP
	DO KK = 1, JFNK_ND
	  DO II = 1, JFNK_NT
	    DO J = 1, JFNK_NT
	      R_TEST((KK-1)*JFNK_NT + II) =
	1       R_TEST((KK-1)*JFNK_NT + II) +
	2       JFNK_BA_A(II, J, KK) *
	3       Z_TEST((KK-1)*JFNK_NT + J)
	    END DO
	  END DO
	END DO
	R_NORM = 0.0_LDP
	DO I = 1, N
	  R_TEST(I) = R_TEST(I) - V_TEST(I)
	  R_NORM = MAX(R_NORM, ABS(R_TEST(I)))
	END DO
	REL_ERR = R_NORM / V_NORM
	WRITE(LUER,'(A,ES12.4)')
	1    '  PC accuracy (v=e1): ||Mz-v||/||v||=', REL_ERR
!
! Test vector 3: unit vector in the last species at first depth
!
	V_TEST = 0.0_LDP
	V_TEST(JFNK_NT) = 1.0_LDP
	V_NORM = 1.0_LDP
	Z_TEST = V_TEST
	CALL JFNK_PRECOND_V1(N, Z_TEST)
	R_TEST = 0.0_LDP
	DO KK = 1, JFNK_ND
	  DO II = 1, JFNK_NT
	    DO J = 1, JFNK_NT
	      R_TEST((KK-1)*JFNK_NT + II) =
	1       R_TEST((KK-1)*JFNK_NT + II) +
	2       JFNK_BA_A(II, J, KK) *
	3       Z_TEST((KK-1)*JFNK_NT + J)
	    END DO
	  END DO
	END DO
	R_NORM = 0.0_LDP
	DO I = 1, N
	  R_TEST(I) = R_TEST(I) - V_TEST(I)
	  R_NORM = MAX(R_NORM, ABS(R_TEST(I)))
	END DO
	REL_ERR = R_NORM / V_NORM
	WRITE(LUER,'(A,ES12.4)')
	1    '  PC accuracy (v=e_NT): ||Mz-v||/||v||=', REL_ERR
!
100	CONTINUE
	WRITE(LUER,'(/,A)')' ============================== END PRECOND AUDIT =============================='
!
	END SUBROUTINE JFNK_PRECOND_AUDIT_V1
!
!=======================================================================
! JFNK_COMP_DIAGNOSTIC_V1: Comprehensive compression-zone diagnostic.
!
! Called once at the first JFNK Newton iteration.  Performs:
!   1. Residual linearity test (is F(x) linear in x?)
!   2. FD JVP convergence study (EPS = 1e-4 .. 1e-8)
!   3. FD-vs-BA decomposition (per-depth worst discrepancies)
!   4. LU preconditioner residual per depth (||M*z - v||)
!   5. Compression-zone profile (|F|, BA diag, SIGMA, V, dV/dr per depth)
!
! All state (jfnk_x, jfnk_fx, ED, T, POPS_WORK, STEQ) is saved and
! restored so the diagnostic is side-effect-free.
!
! review §1-§5, 2026-09-10.
!=======================================================================
	SUBROUTINE JFNK_COMP_DIAGNOSTIC_V1(N, F_BASE, LUER)
	USE SET_KIND_MODULE
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_x, jfnk_fx
	USE MOD_CMFGEN, ONLY: R, V, SIGMA, SIGMA_TRUE, ED, T
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_BA_RIDGE_FAC,
	1   JFNK_FIXED_T, FIXED_NE
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: N, LUER
	DOUBLE PRECISION, INTENT(IN) :: F_BASE(N)
!
	LOGICAL, SAVE :: DIAG_DONE = .FALSE.
	INTEGER I, J, K, KK, II, JJ
	INTEGER IOS_D
	REAL(KIND=LDP) VNORM, XNORM
!
! Workspace
!
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: V_TEST(:), F0(:), F1(:)
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: F2(:), F3(:)
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: JV_FD(:), JV_BA(:)
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: Z_TEST(:), R_TEST(:)
	DOUBLE PRECISION, ALLOCATABLE, SAVE :: X_SAVE(:), FX_SAVE(:)
!
! Linearity test variables
!
	DOUBLE PRECISION LIN_ERR1, LIN_ERR2, F0N, F1N, F2N, F3N
	DOUBLE PRECISION LIN_REL1, LIN_REL2
!
! FD convergence study variables
!
	INTEGER I_EPS
	DOUBLE PRECISION EPS_LIST(5), EPS_VAL, H_STEP
	DOUBLE PRECISION JV_NORM, JV_PREV, JV_RATIO
!
! Per-depth preconditioner variables
!
	DOUBLE PRECISION V_NORM_K, R_NORM_K, Z_NORM_K, REL_ERR_K
	DOUBLE PRECISION WORST_REL_ERR, DG_MAX, DG_MIN
	REAL(KIND=LDP) RIDGE_RATIO_4B
	INTEGER WORST_K, I_WORST
!
! Compression-zone variables
!
	REAL(KIND=LDP) F_ABS_MAX, F_ABS_K, DG_K, SIGMA_K, V_K, DVDR_K
	REAL(KIND=LDP) RAD_K, T_K, ED_K
	INTEGER K_FMAX
!
	IF(DIAG_DONE) RETURN
	DIAG_DONE = .TRUE.
!
	IF(.NOT. JFNK_LU_READY) THEN
	  WRITE(LUER,'(A)') '  COMP DIAGNOSTIC: LU not ready, skipping'
	  RETURN
	END IF
!
	IOS_D = 0
	IF(.NOT. ALLOCATED(V_TEST)) THEN
	  ALLOCATE(V_TEST(N), F0(N), F1(N), F2(N), F3(N),
	1          JV_FD(N), JV_BA(N), Z_TEST(N), R_TEST(N),
	2          X_SAVE(N), FX_SAVE(N), STAT=IOS_D)
	END IF
	IF(IOS_D .NE. 0) THEN
	  WRITE(LUER,'(A)') '  COMP DIAGNOSTIC: alloc failed, skipping'
	  RETURN
	END IF
!
! Save JFNK state
!
	X_SAVE = jfnk_x
	FX_SAVE = jfnk_fx
!
	WRITE(LUER,'(/,A)')' ==============================================================================='
	WRITE(LUER,'(A)')'                   COMPRESSION-ZONE COMPREHENSIVE DIAGNOSTIC                    '
	WRITE(LUER,'(A)')' ==============================================================================='
!
!=======================================================================
! 1. RESIDUAL LINEARITY TEST
!=======================================================================
! With FIX_NE=T, FIX_T=T, no advection, no charge exchange, and frozen J,
! F(x) should be linear: F(x) = F(0) + J*x where J = BA_full.
! Test: F(2v) - F(0) should equal 2*(F(v) - F(0)).
!=======================================================================
!
	WRITE(LUER,'(/,A)')' -------------------------- 1. Residual Linearity Test -------------------------'
	WRITE(LUER,'(A,L1,A,L1)') '  FIX_T=', JFNK_FIXED_T,
	1   '  FIX_NE=', FIXED_NE
!
! Use F_BASE as the test direction (representative of the actual RHS)
!
	V_TEST = F_BASE
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, V_TEST)
	VNORM = 0.0_LDP
	DO I = 1, N
	  VNORM = MAX(VNORM, ABS(V_TEST(I)))
	END DO
	IF(VNORM .LT. 1.0D-30) THEN
	  WRITE(LUER,'(A)') '  F_BASE is zero, skipping linearity test'
	  GO TO 200
	END IF
!
! Normalize V_TEST to O(1) to avoid large perturbations
!
	V_TEST = V_TEST / VNORM
!
! F(0)
!
	jfnk_x = 0.0D0
	CALL JFNK_RESIDUAL_V1(N, jfnk_x, F0)
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, F0)
!
! F(v)
!
	jfnk_x = V_TEST
	CALL JFNK_RESIDUAL_V1(N, jfnk_x, F1)
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, F1)
!
! F(2v)
!
	jfnk_x = 2.0D0 * V_TEST
	CALL JFNK_RESIDUAL_V1(N, jfnk_x, F2)
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, F2)
!
! F(3v)
!
	jfnk_x = 3.0D0 * V_TEST
	CALL JFNK_RESIDUAL_V1(N, jfnk_x, F3)
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, F3)
!
! Linearity check: F(2v) - F(0) should equal 2*(F(v) - F(0))
!                 F(3v) - F(0) should equal 3*(F(v) - F(0))
!
	LIN_ERR1 = 0.0D0
	LIN_ERR2 = 0.0D0
	F0N = 0.0D0
	F1N = 0.0D0
	F2N = 0.0D0
	F3N = 0.0D0
	DO I = 1, N
	  LIN_ERR1 = MAX(LIN_ERR1, ABS(F2(I) - F0(I) - 2.0D0*(F1(I) - F0(I))))
	  LIN_ERR2 = MAX(LIN_ERR2, ABS(F3(I) - F0(I) - 3.0D0*(F1(I) - F0(I))))
	  F0N = MAX(F0N, ABS(F0(I)))
	  F1N = MAX(F1N, ABS(F1(I) - F0(I)))
	  F2N = MAX(F2N, ABS(F2(I) - F0(I)))
	  F3N = MAX(F3N, ABS(F3(I) - F0(I)))
	END DO
!
	LIN_REL1 = LIN_ERR1 / MAX(F2N, 1.0D-30)
	LIN_REL2 = LIN_ERR2 / MAX(F3N, 1.0D-30)
!
	WRITE(LUER,'(A,ES14.4)') '  ||F(0)||_inf = ', F0N
	WRITE(LUER,'(A,ES14.4)') '  ||F(v)-F(0)|| = ', F1N
	WRITE(LUER,'(A,ES14.4)') '  ||F(2v)-F(0)|| = ', F2N
	WRITE(LUER,'(A,ES14.4)') '  ||F(3v)-F(0)|| = ', F3N
	WRITE(LUER,'(A,ES14.4,A,ES10.2)')
	1   '  Linearity err (2v): ', LIN_ERR1,
	1   '  rel=', LIN_REL1
	WRITE(LUER,'(A,ES14.4,A,ES10.2)')
	1   '  Linearity err (3v): ', LIN_ERR2,
	1   '  rel=', LIN_REL2
!
	IF(LIN_REL1 .LT. 1.0D-10 .AND. LIN_REL2 .LT. 1.0D-10) THEN
	  WRITE(LUER,'(A)') '  => F(x) is LINEAR (to machine precision)'
	ELSE IF(LIN_REL1 .LT. 1.0D-4 .AND. LIN_REL2 .LT. 1.0D-4) THEN
	  WRITE(LUER,'(A)') '  => F(x) is WEAKLY nonlinear'
	ELSE
	  WRITE(LUER,'(A)') '  => F(x) is STRONGLY nonlinear'
	  WRITE(LUER,'(A)')
	1   '  => Nonlinear terms present even with FIX_T/FIX_NE'
	END IF
!
!=======================================================================
! 2. FD JVP CONVERGENCE STUDY
!=======================================================================
! Compute J*v ≈ [F(eps*v) - F(0)] / eps for eps = 1e-4 .. 1e-8
! Report how the FD JVP changes with eps, and identify the optimal eps.
!=======================================================================
!
200	CONTINUE
	WRITE(LUER,'(/,A)')' ------------------------- 2. FD JVP Convergence Study -------------------------'
!
! Use the same normalized V_TEST direction
!
	EPS_LIST = (/1.0D-4, 1.0D-5, 1.0D-6, 1.0D-7, 1.0D-8/)
!
! F(0) is already computed above; recompute if we skipped linearity test
!
	IF(VNORM .LT. 1.0D-30) THEN
	  jfnk_x = 0.0D0
	  CALL JFNK_RESIDUAL_V1(N, jfnk_x, F0)
	  CALL JFNK_MASK_FIXED_UNKNOWNS(N, F0)
	  V_TEST = F_BASE
	  CALL JFNK_MASK_FIXED_UNKNOWNS(N, V_TEST)
	  VNORM = 0.0_LDP
	  DO I = 1, N
	    VNORM = MAX(VNORM, ABS(V_TEST(I)))
	  END DO
	  IF(VNORM .LT. 1.0D-30) GO TO 300
	  V_TEST = V_TEST / VNORM
	END IF
!
	JV_PREV = 0.0D0
	DO I_EPS = 1, 5
	  EPS_VAL = EPS_LIST(I_EPS)
	  H_STEP = EPS_VAL
!
	  jfnk_x = H_STEP * V_TEST
	  CALL JFNK_RESIDUAL_V1(N, jfnk_x, F1)
	  CALL JFNK_MASK_FIXED_UNKNOWNS(N, F1)
!
	  JV_NORM = 0.0D0
	  DO I = 1, N
	    JV_FD(I) = (F1(I) - F0(I)) / H_STEP
	    JV_NORM = MAX(JV_NORM, ABS(JV_FD(I)))
	  END DO
	  CALL JFNK_MASK_FIXED_UNKNOWNS(N, JV_FD)
!
	  IF(I_EPS .GT. 1) THEN
	    JV_RATIO = 0.0D0
	    DO I = 1, N
	      JV_RATIO = MAX(JV_RATIO, ABS(JV_FD(I) - JV_BA(I)))
	    END DO
	    WRITE(LUER,'(A,ES8.1,A,ES14.4,A,ES12.4)')
	1     '  eps=', EPS_VAL, '  ||Jv||_inf=', JV_NORM,
	1     '  change_from_prev=', JV_RATIO
	  ELSE
	    WRITE(LUER,'(A,ES8.1,A,ES14.4)')
	1     '  eps=', EPS_VAL, '  ||Jv||_inf=', JV_NORM
	  END IF
!
! Save for comparison with next eps and with BA
!
	  JV_BA = JV_FD
	  JV_PREV = JV_NORM
	END DO
!
!=======================================================================
! 3. FD-vs-BA COMPARISON (per-depth worst)
!=======================================================================
! Compare the best FD JVP (eps=1e-6) with BA_full * v.
! Report the worst depths and the overall relative discrepancy.
!=======================================================================
!
300	CONTINUE
	WRITE(LUER,'(/,A)')' ---------------------- 3. FD-vs-BA Comparison (per depth) ---------------------'
!
! Recompute FD JVP at eps=1e-6
!
	EPS_VAL = 1.0D-6
	H_STEP = EPS_VAL
	jfnk_x = H_STEP * V_TEST
	CALL JFNK_RESIDUAL_V1(N, jfnk_x, F1)
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, F1)
	DO I = 1, N
	  JV_FD(I) = (F1(I) - F0(I)) / H_STEP
	END DO
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, JV_FD)
!
! Analytical BA_full * v (physical Jacobian; ridge is preconditioner-only)
!
	JV_BA = 0.0D0
	DO KK = 1, JFNK_ND
	  DO II = 1, JFNK_NT
	    DO JJ = 1, JFNK_NT
	      JV_BA((KK-1)*JFNK_NT + II) =
	1       JV_BA((KK-1)*JFNK_NT + II) +
	2       DBLE(JFNK_BA_A(II, JJ, KK)) *
	3       V_TEST((KK-1)*JFNK_NT + JJ)
	    END DO
	  END DO
	END DO
	CALL JFNK_MASK_FIXED_UNKNOWNS(N, JV_BA)
!
! Overall discrepancy
!
	JV_NORM = 0.0D0
	WORST_REL_ERR = 0.0D0
	WORST_K = 1
	DO I = 1, N
	  JV_NORM = MAX(JV_NORM, ABS(JV_FD(I)))
	END DO
!
	WRITE(LUER,'(A)') '  Per-depth ||FD-BA||/||BA|| (worst 10):'
!
! Find and report the 10 worst depths
!
	DO KK = 1, JFNK_ND
	  V_NORM_K = 0.0D0
	  R_NORM_K = 0.0D0
	  DO II = 1, JFNK_NT
	    I = (KK-1)*JFNK_NT + II
	    V_NORM_K = MAX(V_NORM_K, ABS(JV_BA(I)))
	    R_NORM_K = MAX(R_NORM_K, ABS(JV_FD(I) - JV_BA(I)))
	  END DO
	  REL_ERR_K = R_NORM_K / MAX(V_NORM_K, 1.0D-30)
	  IF(REL_ERR_K .GT. WORST_REL_ERR) THEN
	    WORST_REL_ERR = REL_ERR_K
	    WORST_K = KK
	  END IF
	END DO
!
! Report worst depths (just the single worst for brevity)
!
	KK = WORST_K
	V_NORM_K = 0.0D0
	R_NORM_K = 0.0D0
	DO II = 1, JFNK_NT
	  I = (KK-1)*JFNK_NT + II
	  V_NORM_K = MAX(V_NORM_K, ABS(JV_BA(I)))
	  R_NORM_K = MAX(R_NORM_K, ABS(JV_FD(I) - JV_BA(I)))
	END DO
	WRITE(LUER,'(A,I4,A,ES12.4,A,ES12.4,A,ES12.4)')
	1   '  Worst depth K=', KK,
	1   '  ||BA*v||_inf=', V_NORM_K,
	1   '  ||FD-BA||_inf=', R_NORM_K,
	1   '  rel=', R_NORM_K/MAX(V_NORM_K,1.0D-30)
!
! Report overall
!
	R_NORM_K = 0.0D0
	V_NORM_K = 0.0D0
	DO I = 1, N
	  R_NORM_K = MAX(R_NORM_K, ABS(JV_FD(I) - JV_BA(I)))
	  V_NORM_K = MAX(V_NORM_K, ABS(JV_BA(I)))
	END DO
	WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1   '  Overall: ||FD-BA||/||BA||=',
	1   R_NORM_K/MAX(V_NORM_K,1.0D-30),
	1   '  ||FD||_inf=', JV_NORM,
	1   '  ||BA||_inf=', V_NORM_K
!
!=======================================================================
! 4. LU PRECONDITIONER RESIDUAL PER DEPTH
!=======================================================================
! For v = F_BASE, compute z = M^{-1}*v, then ||M*z - v|| per depth.
! Report the worst depths.
!=======================================================================
!
	WRITE(LUER,'(/,A)')' ------------------ 4. LU Preconditioner Residual (per depth) ------------------'
!
	Z_TEST = F_BASE
	CALL JFNK_PRECOND_V1(N, Z_TEST)
!
! Compute M*z using JFNK_BA_A (ridge-augmented BA)
!
	R_TEST = 0.0D0
	DO KK = 1, JFNK_ND
	  DO II = 1, JFNK_NT
	    DO JJ = 1, JFNK_NT
	      R_TEST((KK-1)*JFNK_NT + II) =
	1       R_TEST((KK-1)*JFNK_NT + II) +
	2       DBLE(JFNK_BA_A(II, JJ, KK)) *
	3       Z_TEST((KK-1)*JFNK_NT + JJ)
	    END DO
	  END DO
	END DO
!
! r = M*z - v, per depth
!
	WORST_REL_ERR = 0.0D0
	WORST_K = 1
	DO KK = 1, JFNK_ND
	  V_NORM_K = 0.0D0
	  R_NORM_K = 0.0D0
	  Z_NORM_K = 0.0D0
	  DO II = 1, JFNK_NT
	    I = (KK-1)*JFNK_NT + II
	    R_TEST(I) = R_TEST(I) - F_BASE(I)
	    V_NORM_K = MAX(V_NORM_K, ABS(F_BASE(I)))
	    R_NORM_K = MAX(R_NORM_K, ABS(R_TEST(I)))
	    Z_NORM_K = MAX(Z_NORM_K, ABS(Z_TEST(I)))
	  END DO
	  REL_ERR_K = R_NORM_K / MAX(V_NORM_K, 1.0D-30)
	  IF(REL_ERR_K .GT. WORST_REL_ERR) THEN
	    WORST_REL_ERR = REL_ERR_K
	    WORST_K = KK
	  END IF
	END DO
!
	WRITE(LUER,'(A,I4,A,ES12.4)')
	1   '  Worst depth K=', WORST_K,
	1   '  ||Mz-v||/||v||=', WORST_REL_ERR
!
! Report worst 5 depths
!
	WRITE(LUER,'(A)') '  Per-depth residuals (worst 5):'
	DO I_WORST = 1, 5
	  IF(I_WORST .EQ. 1) KK = WORST_K
	  IF(I_WORST .EQ. 2) KK = MAX(1, WORST_K-1)
	  IF(I_WORST .EQ. 3) KK = MIN(JFNK_ND, WORST_K+1)
	  IF(I_WORST .EQ. 4) KK = MAX(1, WORST_K-2)
	  IF(I_WORST .EQ. 5) KK = MIN(JFNK_ND, WORST_K+2)
	  V_NORM_K = 0.0D0
	  R_NORM_K = 0.0D0
	  Z_NORM_K = 0.0D0
	  DG_MAX = 0.0_LDP
	  DG_MIN = 1.0E+30_LDP
	  DO II = 1, JFNK_NT
	    I = (KK-1)*JFNK_NT + II
	    V_NORM_K = MAX(V_NORM_K, ABS(F_BASE(I)))
	    R_NORM_K = MAX(R_NORM_K, ABS(R_TEST(I)))
	    Z_NORM_K = MAX(Z_NORM_K, ABS(Z_TEST(I)))
	    DG_K = ABS(JFNK_BA_DIAG(II, KK))
	    DG_MAX = MAX(DG_MAX, DG_K)
	    DG_MIN = MIN(DG_MIN, DG_K)
	  END DO
	  REL_ERR_K = R_NORM_K / MAX(V_NORM_K, 1.0D-30)
	  WRITE(LUER,'(A,I4,A,ES10.2,A,ES12.4,A,ES10.2,A,ES10.2)')
	1   '    K=', KK,
	1   '  rel_err=', REL_ERR_K,
	1   '  ||v||=', V_NORM_K,
	1   '  ||z||=', Z_NORM_K,
	1   '  diag_range=', DG_MAX/MAX(DG_MIN,1.0E-30_LDP)
	END DO
!
!=======================================================================
! 4b. WORST-LU-DEPTH BA STRUCTURE DIAGNOSTIC
!=======================================================================
! Print per-species BA diagonal, LU solution, residual, and
! equilibration factors at the worst-LU depth to identify which
! species/level causes the ZMAX hit and ill-conditioning.
!=======================================================================
!
	WRITE(LUER,'(/,A,I4,A)')
	1   ' ---------- 4b. BA Structure at Worst-LU Depth K=',WORST_K,
	2   ' --------------------------'
	WRITE(LUER,'(A)')
	1   '  sp  |diag|        ridge/diag   |z|          |Mz-v|       '//
	1   'ROW_SF       COL_SF'
	DO II = 1, JFNK_NT
	  I = (WORST_K-1)*JFNK_NT + II
	  DG_K = ABS(JFNK_BA_DIAG(II, WORST_K))
	  RIDGE_RATIO_4B = JFNK_BA_RIDGE_ADD(II, WORST_K) /
	1               MAX(DG_K, 1.0E-30_LDP)
	  WRITE(LUER,'(I4,1P,6(1X,ES10.3))')
	1   II, DG_K, RIDGE_RATIO_4B, ABS(Z_TEST(I)),
	1   ABS(R_TEST(I)), JFNK_ROW_SF(II, WORST_K),
	1   JFNK_COL_SF(II, WORST_K)
	END DO
!
! Identify species with |z| > 1.0 (near or at ZMAX cap)
!
	WRITE(LUER,'(A)') '  Species near ZMAX cap (|z|>1.0):'
	DO II = 1, JFNK_NT
	  I = (WORST_K-1)*JFNK_NT + II
	  IF(ABS(Z_TEST(I)) .GT. 1.0D0)THEN
	    DG_K = ABS(JFNK_BA_DIAG(II, WORST_K))
	    WRITE(LUER,'(A,I4,A,1P,4(ES10.3))')
	1     '    sp', II, ': |z|=', ABS(Z_TEST(I)),
	1     DG_K,
	1     JFNK_BA_RIDGE_ADD(II, WORST_K)/MAX(DG_K,1.0E-30_LDP),
	1     ABS(R_TEST(I))
	  END IF
	END DO
	WRITE(LUER,'(/,A)')' ------------------------------------ End 4b -----------------------------------'
!
!=======================================================================
! 5. COMPRESSION-ZONE PROFILE
!=======================================================================
! For each depth, report |F|, BA diagonal, SIGMA, velocity, dV/dr.
! Focus on depths 35-65 (the compression zone).
!=======================================================================
!
	WRITE(LUER,'(/,A)')' ------------------------- 5. Compression-Zone Profile -------------------------'
	WRITE(LUER,'(A)')
	1   '  K    R(10^10cm)  V(km/s)  SIGMA   dV/dr   T(K)    '//
	1   'ED(cm-3)  |F|_max  |diag|_max  |diag|_min'
!
	K_FMAX = 1
	F_ABS_MAX = 0.0_LDP
!
	DO K = 1, JFNK_ND
	  F_ABS_K = 0.0_LDP
	  DG_MAX = 0.0_LDP
	  DG_MIN = 1.0E+30_LDP
	  DO II = 1, JFNK_NT
	    I = (K-1)*JFNK_NT + II
	    F_ABS_K = MAX(F_ABS_K, ABS(F_BASE(I)))
	    DG_K = ABS(JFNK_BA_DIAG(II, K))
	    DG_MAX = MAX(DG_MAX, DG_K)
	    DG_MIN = MIN(DG_MIN, DG_K)
	  END DO
!
	  IF(F_ABS_K .GT. F_ABS_MAX) THEN
	    F_ABS_MAX = F_ABS_K
	    K_FMAX = K
	  END IF
!
! Report depths 35-65 and the worst depth
!
	  IF((K .GE. 35 .AND. K .LE. 65) .OR.
	1    K .EQ. 1 .OR. K .EQ. JFNK_ND) THEN
	    SIGMA_K = SIGMA_TRUE(K)
	    V_K = V(K)
	    RAD_K = R(K)
	    T_K = T(K)
	    ED_K = ED(K)
	    IF(K .GT. 1 .AND. K .LT. JFNK_ND) THEN
	      DVDR_K = (V(K+1) - V(K-1)) / (R(K+1) - R(K-1))
	    ELSE
	      DVDR_K = 0.0_LDP
	    END IF
	    WRITE(LUER,'(I4,F12.4,F10.2,F8.2,ES10.2,F10.0,ES12.3,'//
	1       'ES12.4,ES12.4,ES12.4)')
	1     K, RAD_K, V_K, SIGMA_K, DVDR_K, T_K, ED_K,
	1     F_ABS_K, DG_MAX, DG_MIN
	  END IF
	END DO
!
	WRITE(LUER,'(A,I4,A,ES12.4)')
	1   '  Largest |F| at depth K=', K_FMAX,
	1   '  |F|_max=', F_ABS_MAX
!
! Detailed report for the worst-F depth
!
	K = K_FMAX
	WRITE(LUER,'(A,I4)') '  Per-species |F| at worst-F depth K=', K
	DO II = 1, MIN(JFNK_NT, 20)
	  I = (K-1)*JFNK_NT + II
	  DG_K = ABS(JFNK_BA_DIAG(II, K))
	  WRITE(LUER,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)')
	1   '    sp', II, ': |F|=', ABS(F_BASE(I)),
	1   '  |diag|=', DG_K,
	1   '  |F/diag|=', ABS(F_BASE(I))/MAX(DG_K,1.0E-30_LDP)
	END DO
!
!=======================================================================
! Restore JFNK state
!=======================================================================
!
	jfnk_x = X_SAVE
	jfnk_fx = FX_SAVE
!
	WRITE(LUER,'(/,A)')' ==============================================================================='
	WRITE(LUER,'(A)')'                        END COMPRESSION-ZONE DIAGNOSTIC                         '
	WRITE(LUER,'(A)')' ==============================================================================='
	WRITE(LUER,'(A)') ' '
	FLUSH(LUER)
!
	END SUBROUTINE JFNK_COMP_DIAGNOSTIC_V1
!
!
!
! JFNK_ROW_ANATOMY_V1: residual anatomy at admission.  Prints the
! collisional / stored-radiation split of F, the per-ion rollup of
! |F|^2, and the top JFNK_ROW_DIAG rows labeled by ion, level and
! depth.  Jz_i (the row content of the linear preconditioned step) is
! printed when the reach precheck supplies it: a large |F_i| with
! Jz_i ~ 0 marks a row no population step can reach, which is what
! separates radiation-lag rows from population-driven ones and
! decides whether a population-space residual can ever converge.
!
! jfnk_fx = coll + rad at admission (x=0), so coll_i = F_i - rad_i is
! exact without touching JFNK_W_STEQ_MAT -- the preconditioner reuses
! that buffer as scratch, so reading it here would print solve
! internals, not the collisional STEQ.  Unweighted (physical) norms
! are used: RSCALE rescaling would hide which equations actually
! carry the residual.
!
	SUBROUTINE JFNK_ROW_ANATOMY_V1(LUER,JZ_VEC,N_JZ)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD, ONLY: JFNK_ROW_DIAG
	USE JFNK_CMFGEN_INT_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_fx
	IMPLICIT NONE
!
	INTEGER LUER,N_JZ
	DOUBLE PRECISION JZ_VEC(*)
!
	INTEGER I,J,K,ID,II,N,NTOP,NRK,NGRP,GRP
	INTEGER, ALLOCATABLE :: J2ID(:),ORD(:),GORD(:)
	DOUBLE PRECISION, ALLOCATABLE :: GRP_F2(:),GRP_C2(:),
	1   GRP_R2(:),GRP_JZ2(:),GRP_BA2(:)
	DOUBLE PRECISION F2T,C2T,R2T,JZ2T
	REAL(KIND=LDP) C_I,R_I,BRN
	LOGICAL HAVE_BA,HAVE_JZ
	CHARACTER*14 LBL
!
	IF(LUER .LE. 0)RETURN
	IF(JFNK_ROW_DIAG .LE. 0)RETURN
	IF(JFNK_NT .LE. 0 .OR. JFNK_ND .LE. 0)RETURN
	IF(.NOT. ALLOCATED(jfnk_fx))RETURN
	N=JFNK_NT*JFNK_ND
	IF(SIZE(jfnk_fx) .LT. N)RETURN
	HAVE_JZ=(N_JZ .GE. N)
	HAVE_BA=ALLOCATED(JFNK_BA_A) .AND. JFNK_LU_READY
!
! Variable index J follows the POPS ordering: levels
! EQXzV..EQXzV+NXzV-1 belong to ion ID; the highest stage of each
! species contributes its single conservation row at EQXzV.  J=NT-1
! is the charge (Ne) row and J=NT the energy (T) row.
!
	ALLOCATE(J2ID(JFNK_NT))
	J2ID=0
	DO ID=1,JFNK_NION
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO J=ATM(ID)%EQXzV,
	1        MIN(ATM(ID)%EQXzV+ATM(ID)%NXzV-1,JFNK_NT-2)
	      IF(J .GE. 1)J2ID(J)=ID
	    END DO
	  ELSE
	    IF(ATM(ID)%EQXzV .GE. 1 .AND.
	1       ATM(ID)%EQXzV .LE. JFNK_NT-2)J2ID(ATM(ID)%EQXzV)=ID
	  END IF
	END DO
!
! Groups 1..NION are the ions, NION+1 the charge rows, NION+2 the
! energy rows, NION+3 any row the J->ion map did not claim.
!
	NGRP=JFNK_NION+3
	ALLOCATE(GRP_F2(NGRP),GRP_C2(NGRP),GRP_R2(NGRP),
	1   GRP_JZ2(NGRP),GRP_BA2(NGRP))
	GRP_F2=0.0D0
	GRP_C2=0.0D0
	GRP_R2=0.0D0
	GRP_JZ2=0.0D0
	GRP_BA2=0.0D0
	F2T=0.0D0
	C2T=0.0D0
	R2T=0.0D0
	JZ2T=0.0D0
!
	DO I=1,N
	  K=(I-1)/JFNK_NT+1
	  J=MOD(I-1,JFNK_NT)+1
	  R_I=0.0_LDP
	  IF(ALLOCATED(JFNK_STEQ_RAD_FLAT))R_I=JFNK_STEQ_RAD_FLAT(I)
	  C_I=REAL(jfnk_fx(I),KIND=LDP)-R_I
	  BRN=0.0_LDP
	  IF(HAVE_BA)THEN
	    DO II=1,JFNK_NT
	      BRN=BRN+JFNK_BA_A(J,II,K)**2
	    END DO
	    BRN=SQRT(BRN)
	  END IF
	  IF(J .EQ. JFNK_NT)THEN
	    GRP=JFNK_NION+2
	  ELSE IF(J .EQ. JFNK_NT-1)THEN
	    GRP=JFNK_NION+1
	  ELSE IF(J2ID(J) .GT. 0)THEN
	    GRP=J2ID(J)
	  ELSE
	    GRP=JFNK_NION+3
	  END IF
	  GRP_F2(GRP)=GRP_F2(GRP)+jfnk_fx(I)**2
	  GRP_C2(GRP)=GRP_C2(GRP)+DBLE(C_I)**2
	  GRP_R2(GRP)=GRP_R2(GRP)+DBLE(R_I)**2
	  GRP_BA2(GRP)=GRP_BA2(GRP)+DBLE(BRN)**2
	  IF(HAVE_JZ)GRP_JZ2(GRP)=GRP_JZ2(GRP)+JZ_VEC(I)**2
	  F2T=F2T+jfnk_fx(I)**2
	  C2T=C2T+DBLE(C_I)**2
	  R2T=R2T+DBLE(R_I)**2
	  IF(HAVE_JZ)JZ2T=JZ2T+JZ_VEC(I)**2
	END DO
!
	WRITE(LUER,'(A,ES12.4,A,ES12.4,A,ES12.4)')
	1   '  JFNK ANATOMY: |F|=',REAL(SQRT(F2T)),
	2   ' |coll|=',REAL(SQRT(C2T)),
	3   ' |rad|=',REAL(SQRT(R2T))
	IF(HAVE_JZ)WRITE(LUER,'(A,ES12.4)')
	1   '    |Jz|=',REAL(SQRT(JZ2T))
!
! Per-ion rollup, largest |F| share first.
!
	ALLOCATE(GORD(NGRP))
	DO I=1,NGRP
	  GORD(I)=I
	END DO
	NRK=MIN(8,NGRP)
	DO I=1,NRK
	  DO K=I+1,NGRP
	    IF(GRP_F2(GORD(K)) .GT. GRP_F2(GORD(I)))THEN
	      II=GORD(I)
	      GORD(I)=GORD(K)
	      GORD(K)=II
	    END IF
	  END DO
	END DO
	WRITE(LUER,'(A)')
	1   '    ion           Fshare    coll/F    rad/F    Jz/F'
	2   //'     BArow/F'
	DO I=1,NRK
	  GRP=GORD(I)
	  IF(GRP_F2(GRP) .LE. 0.0D0)GOTO 100
	  IF(GRP .LE. JFNK_NION)THEN
	    LBL=ION_ID(GRP)
	  ELSE IF(GRP .EQ. JFNK_NION+1)THEN
	    LBL='Ne (charge)'
	  ELSE IF(GRP .EQ. JFNK_NION+2)THEN
	    LBL='T (energy)'
	  ELSE
	    LBL='unlabeled'
	  END IF
	  WRITE(LUER,'(4X,A12,5ES10.2)')LBL,
	1     REAL(SQRT(GRP_F2(GRP)/MAX(F2T,1.0D-300))),
	2     REAL(SQRT(GRP_C2(GRP)/MAX(GRP_F2(GRP),1.0D-300))),
	3     REAL(SQRT(GRP_R2(GRP)/MAX(GRP_F2(GRP),1.0D-300))),
	4     REAL(SQRT(GRP_JZ2(GRP)/MAX(GRP_F2(GRP),1.0D-300))),
	5     REAL(SQRT(GRP_BA2(GRP)/MAX(GRP_F2(GRP),1.0D-300)))
	END DO
  100	CONTINUE
!
! Largest residual rows by |F_i|.
!
	NTOP=MIN(JFNK_ROW_DIAG,N)
	ALLOCATE(ORD(N))
	DO I=1,N
	  ORD(I)=I
	END DO
	DO I=1,NTOP
	  DO K=I+1,N
	    IF(jfnk_fx(ORD(K))**2 .GT. jfnk_fx(ORD(I))**2)THEN
	      II=ORD(I)
	      ORD(I)=ORD(K)
	      ORD(K)=II
	    END IF
	  END DO
	END DO
	WRITE(LUER,'(A)')
	1   '    rnk    idx     K     J  ion          lev'
	2   //'     F_i        coll_i     rad_i      Jz_i       BArown'
	DO I=1,NTOP
	  II=ORD(I)
	  K=(II-1)/JFNK_NT+1
	  J=MOD(II-1,JFNK_NT)+1
	  R_I=0.0_LDP
	  IF(ALLOCATED(JFNK_STEQ_RAD_FLAT))R_I=JFNK_STEQ_RAD_FLAT(II)
	  C_I=REAL(jfnk_fx(II),KIND=LDP)-R_I
	  BRN=0.0_LDP
	  IF(HAVE_BA)THEN
	    DO ID=1,JFNK_NT
	      BRN=BRN+JFNK_BA_A(J,ID,K)**2
	    END DO
	    BRN=SQRT(BRN)
	  END IF
	  IF(J .EQ. JFNK_NT)THEN
	    LBL='T (energy)'
	    ID=0
	  ELSE IF(J .EQ. JFNK_NT-1)THEN
	    LBL='Ne (charge)'
	    ID=0
	  ELSE IF(J2ID(J) .GT. 0)THEN
	    LBL=ION_ID(J2ID(J))
	    ID=J-ATM(J2ID(J))%EQXzV+1
	  ELSE
	    LBL='unlabeled'
	    ID=0
	  END IF
	  IF(HAVE_JZ)THEN
	    WRITE(LUER,'(I5,I8,I6,I6,2X,A12,I5,5ES11.3)')I,II,K,J,
	1       LBL,ID,jfnk_fx(II),C_I,R_I,JZ_VEC(II),BRN
	  ELSE
	    WRITE(LUER,'(I5,I8,I6,I6,2X,A12,I5,4ES11.3)')I,II,K,J,
	1       LBL,ID,jfnk_fx(II),C_I,R_I,BRN
	  END IF
	END DO
	DEALLOCATE(J2ID,ORD,GORD,GRP_F2,GRP_C2,GRP_R2,GRP_JZ2,GRP_BA2)
	FLUSH(LUER)
	RETURN
	END SUBROUTINE JFNK_ROW_ANATOMY_V1
