!
! File to check the consistency of CMFGEN parameters
!
! Altered: 23-May-2010 - Fixed pontential compile problems with .NE. and logical operators.
! Created: 07-Jan-2008
!
	SUBROUTINE CHECK_PARAM_CONSISTENCY( )
	USE SET_KIND_MODULE
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	CHARACTER(LEN=10) ER_LAB
!
	ER_LAB='Warning:'
	IF(STOP_IF_BAD_PARAM)ER_LAB='ERROR:'
	LUER=ERROR_LU()
!	
	IF(SN_MODEL)THEN
	  IF(.NOT. DO_FULL_REL_OBS)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' DO_FULL_RELL_OBS should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(.NOT. DO_FULL_REL_CMF)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' DO_FULL_REL_CMF should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
!
! edited 22-Oct-2026: the DO_NONMONO_* flags legitimately
! run SN models with all three off -- under USE_LAM_ES=F the J update
! comes from MOM_J_CMF_NONMONO_V1 (Knop sign-split moment solve on
! SIGMA_TRUE), and rd_control_variables announces that pairing with
! the 'explicit USE_LAM_ES=F honoured' banner.  Exempt the nonmono
! modes so an honoured F is not rejected here; other SN checks and
! the stock path are unchanged.
!
	  IF(.NOT. USE_J_REL .AND. .NOT. USE_DJDT_RTE .AND.
	1     .NOT. USE_LAM_ES .AND.
	2     .NOT. (DO_NONMONO_KNOP .OR. DO_NONMONO_SOB))THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' USE_JREL or USE_DJDT_RTE should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	   END IF
	  IF(USE_J_REL .AND. .NOT. INCL_REL_TERMS)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' INCL_REL should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(USE_J_REL .AND. .NOT. INCL_ADVEC_TERMS_IN_TRANS_EQ)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' INCL_ADV_TRANS should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
!
	  IF(USE_DJDT_RTE .AND. .NOT. INCL_ADIABATIC)THEN
	    WRITE(LUER,*)'Error: inconsistency in control parameters in VADAT'
	    WRITE(LUER,*)'SN model has USE_DJDT_RTE set to TRUE.'
	    WRITE(LUER,*)'INCL_ADIABATIC must also be set to TRUE.'
	    STOP
	  END IF
	  IF(USE_DJDT_RTE .AND. .NOT. DO_CO_MOV_DDT)THEN
	    WRITE(LUER,*)'Error: inconsistency in control parameters in VADAT'
	    WRITE(LUER,*)'SN model has USE_DJDT_RTE set to TRUE.'
	    WRITE(LUER,*)'DO_CO_MOV_DDT should also be set to TRUE.'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
!
! AS TIME_SEQ_NO is not actually used, this check is not necessary.
!
!	  IF(USE_DJDT_RTE .AND. (TIME_SEQ_NO .EQ. 1) )THEN
!	    WRITE(LUER,*)'Error: inconsistency in control parameters in VADAT'
!	    WRITE(LUER,*)'SN model has TS_NO=1 implying initial model'
!	    WRITE(LUER,*)TRIM(ER_LAB),'USE_JREL=TRUE and USE_DJDT_RTE=FALSE for an initial model'
!	    STOP
!	  END IF
!
	  IF(USE_J_REL .AND. USE_DJDT_RTE)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),'Error: USE_JREL or USE_DJDT_RTE cannot be true at the same time.'
	    WRITE(LUER,*)TRIM(ER_LAB),' For initial time seqnece model, only USE_JRL should be true.'
	    STOP
	  END IF
	  IF(PLANE_PARALLEL  .OR. PLANE_PARALLEL_NO_V)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'PP_NOV and PP_MOD cannot be TRUE for a SN model'
	    STOP
	  END IF
	END IF
!
	IF(PLANE_PARALLEL  .AND. PLANE_PARALLEL_NO_V)THEN
	  WRITE(LUER,*)'Error in control parameters in VADAT'
	  WRITE(LUER,*)'PP_NOV and PP_MOD cannot both be TRUE'
	  WRITE(LUER,*)'Use PP_NOV for a plane-parallel model with no velocity field'
	  WRITE(LUER,*)'Use PP_MOD for a plane-parallel model with a velocity field'
	  STOP
	END IF
!
! The [BA_PROBE_LAM] escape probe arms only when residual acceptance
! judges the trial (DO_BA_RES_ACC) and a factor bound holds it
! (BA_DAMP_CAP); without both, the key reads >0 yet no probe can
! ever fire. Warn rather than stop -- the run is legal, the dead
! key is just misleading in OUTGEN.
!
	IF(BA_PROBE_LAM .GT. 0 .AND.
	1   (.NOT. DO_BA_RES_ACC .OR. BA_DAMP_CAP .LE. 0.0_LDP))THEN
	  WRITE(LUER,*)TRIM(ER_LAB),' BA_PROBE_LAM > 0 but the probe'//
	1     ' cannot fire: it requires DO_BA_RES_ACC=T and'//
	2     ' BA_DAMP_CAP > 0'
	END IF
!
! DO_JFNK + USE_LAM_ES=F is a permanently dead combination (audit
! T-03/F-02, closed 04-Nov-2026): admission in SOLVE_FOR_POPS requires
! USE_LAM_ES because the eval-time lagged-field freeze covers only
! the Knop I_NEXT store (FREEZE_I_NEXT in FG_J_CMF_NONMONO_V1). The
! MOM_J_CMF_NONMONO_V1 J_NEXT_FRQ/H_NEXT_FRQ store -- the J path on
! every nonmono LAM_ES=F route -- has no equivalent, so each
! residual/JVP eval would overwrite the committed lagged coupling.
! Opening it needs a MOM eval freeze plus coupled-only admission
! (the frozen tier is the same fixed-J system SOLVEBA solves) and a
! fresh determinism audit of the MOM module state -- substantial
! work for unproven gain, so the combination is rejected here
! instead of reporting armed-but-always-deferred. DO_JFNK already
! requires a nonmono flag, and DO_SOB_FG_KNOP forces USE_LAM_ES=T,
! so this predicate is exactly the MOM-store routes. Stock configs
! are untouched (USE_JFNK_SOLVER is F there or already stopped).
!
	IF(USE_JFNK_SOLVER .AND. .NOT. USE_LAM_ES)THEN
	  WRITE(LUER,*)'Error in control parameters in VADAT'
	  WRITE(LUER,*)'DO_JFNK=T requires USE_LAM_ES=T: the JFNK'//
	1     ' residual evals have no freeze for the'//
	2     ' MOM_J_CMF_NONMONO_V1 J_NEXT/H_NEXT lagged store.'
	  WRITE(LUER,*)'Set USE_LAM_ES=T (Knop formal or V13 path),'//
	1     ' or set DO_JFNK=F.'
	  STOP
	END IF
!
	RETURN
	END	
