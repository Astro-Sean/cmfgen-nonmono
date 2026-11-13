!!
! Write the Doppler-profile opacity and emissivity of each
! TRANS_TYPE='SOB' line of the current blend group into the stored
! CMF-frame arrays ETA_ST/CHI_ST.  The CMF blanketing solve skips SOB
! lines (comoving-frame line transfer is invalid on a non-monotonic
! velocity field), so they reach the observed spectrum only through
! the observer's-frame formal integral OBS_FRAME_SUB_V9, which samples
! these arrays at each ray point's local comoving frequency.  Called
! per blend group when DO_SOB_FLUX_OBS=T; with the default F the
! stored arrays are left untouched.
! This is the only copy of the routine: cmf_flux_sub_v5.f INCLUDEs
! this file, and the unit tests include the same text, so the test
! and production paths cannot diverge.
!
! The profile is a normalized Doppler core of width
! FL_SIM*VDOP_MIN/C_KMS evaluated on the (descending) CMF grid NU.
! Only grid cells inside +/-3 Doppler widths carry weight, and the
! weights are normalized by the local cell width so the deposited
! opacity integrated over frequency equals CHIL_MAT independent of
! how coarsely the CMF grid resolves the line.  The stored arrays
! are unclumped; the caller folds CLUMP_FAC into the render copies
! once, so the deposit must not carry it.  The
! observed-frame profile shape comes from the ray geometry (the
! isovelocity surfaces), not from this profile width, so a coarse CMF
! grid does not broaden the rendered line.
!
! CHIL_MAT can be negative at depths where stimulated emission
! dominates (masering lines).  A formal integral cannot represent
! maser amplification (chi<0 grows the intensity exponentially and
! chi->0+ leaves a pathologically large two-level source function),
! so a depth where this line masers contributes neither opacity nor
! emissivity: the cell stays transparent to the continuum and other
! lines.  Depths flagged by NEG_OPACITY (the blend-group summed
! Sobolev optical depth was strongly negative, so CHECK_LINE_OPAC
! replaced each member's CHIL_MAT with 1/NUM_SIM_LINES) are dropped
! on the same footing: their surviving ETAL_MAT is real spontaneous
! emissivity, but the artificial opacity would inject it as pure
! thin emission with no chance of self-absorption.  N_DROP
! accumulates the dropped (line,depth) count and DROP_SUM the
! dropped emissivity for the caller's diagnostic print.
! ETAL_MAT is non-negative by construction.
!
	SUBROUTINE ADD_SOB_LINE_PROF(ETA_ST,CHI_ST,NU,NCF,
     1	      CHIL_MAT,ETAL_MAT,CLUMP_FAC,FL_SIM,SIM_POINTER,
     2	      VDOP_MIN,ARAD_VEC,C4_VEC,ED_IN,
     3	      TRANS_TYPE,NEG_OPACITY,
     4	      NUM_SIM,N_LINES,MAX_SIM,ND,C_KMS,
     5	      N_DROP,DROP_SUM,INJ_SUM)
	USE SET_KIND_MODULE
	USE CMF_FLUX_CNTRL_VAR_MOD
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: NCF,ND,NUM_SIM,N_LINES,MAX_SIM
	INTEGER, INTENT(INOUT) :: N_DROP
	REAL(KIND=LDP), INTENT(INOUT) :: DROP_SUM,INJ_SUM
	REAL(KIND=LDP), INTENT(INOUT) :: ETA_ST(ND,NCF)
	REAL(KIND=LDP), INTENT(INOUT) :: CHI_ST(ND,NCF)
	REAL(KIND=LDP), INTENT(IN) :: NU(NCF)
	REAL(KIND=LDP), INTENT(IN) :: CHIL_MAT(ND,MAX_SIM)
	REAL(KIND=LDP), INTENT(IN) :: ETAL_MAT(ND,MAX_SIM)
	REAL(KIND=LDP), INTENT(IN) :: CLUMP_FAC(ND)
	REAL(KIND=LDP), INTENT(IN) :: FL_SIM(MAX_SIM)
	REAL(KIND=LDP), INTENT(IN) :: VDOP_MIN(N_LINES)
	REAL(KIND=LDP), INTENT(IN) :: ARAD_VEC(N_LINES)
	REAL(KIND=LDP), INTENT(IN) :: C4_VEC(N_LINES)
	REAL(KIND=LDP), INTENT(IN) :: ED_IN(ND)
	REAL(KIND=LDP), INTENT(IN) :: C_KMS
	LOGICAL, INTENT(IN) :: NEG_OPACITY(ND)
	INTEGER, INTENT(IN) :: SIM_POINTER(MAX_SIM)
	CHARACTER(LEN=*), INTENT(IN) :: TRANS_TYPE(N_LINES)
!
	INTEGER I,K,S,LP,K0,KA,KB,LO,HI,MID
	INTEGER N_DROP_L
	REAL(KIND=LDP) F0,DNU_D,W,XNORM,WT,DN_K,X,T1
	REAL(KIND=LDP) A_VOIGT,GAM_COL_LINE,TMP_ED,PI
	REAL(KIND=LDP) VOIGT
	EXTERNAL VOIGT
	PARAMETER(PI=3.1415926535897932_LDP)
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! The binary search below assumes NU decreases strictly with index
! (the CMF grid is built that way).  A non-monotonic NU would place
! every line at the wrong frequencies, so check once per call and
! abort rather than deposit silently wrong profiles.
!
	LUER=ERROR_LU()
	DO K=1,NCF-1
	  IF(NU(K+1) .GE. NU(K))THEN
	    WRITE(LUER,*)'ADD_SOB_LINE_PROF: NU not strictly'
	1   //' decreasing; first bad index =',K+1
	    WRITE(LUER,*)'   NU(1),NU(K+1),NU(NCF) =',
	1   NU(1),NU(K+1),NU(NCF)
	    STOP
	  END IF
	END DO
!
	DO S=1,NUM_SIM
	  LP=SIM_POINTER(S)
	  IF(TRANS_TYPE(LP)(1:3) .NE. 'SOB')CYCLE
	  F0=FL_SIM(S)
	  DNU_D=F0*VDOP_MIN(LP)/C_KMS
	  IF(DNU_D .LE. 0.0_LDP)CYCLE
	  W=3.0_LDP*DNU_D
!
! DO_SOB_DEP_WINGS (default F): the full CMF solve carries Voigt/
! damping wings on each line, so the pure-Gaussian deposit renders
! line wings systematically too narrow.  Extend the deposit window
! to SOB_DEP_WING_EXT km/s and continue the profile as 1/x^2 tails
! beyond 3 Doppler widths (a standard Voigt-wing surrogate; the
! wing amplitude inherits the line's opacity so only strong lines
! grow visible wings).
!
	  IF(DO_SOB_DEP_WINGS)W=MAX(W,F0*SOB_DEP_WING_EXT/C_KMS)
!
! NU decreases with index.  K0 is the last index with NU(K0) >= F0,
! so F0 lies in [NU(K0+1),NU(K0)).  Lines outside the grid cannot
! deposit anything; name them in LUER so a missing observed line is
! traceable to the CMF window rather than the populations.
!
	  IF(NU(1) .LT. F0 .OR. NU(NCF) .GE. F0)THEN
	    WRITE(LUER,'(A,I6,A,3ES16.7,A)')
	1     'ADD_SOB_LINE_PROF: SOB line',LP,
	2     ' F0,NU(NCF),NU(1) =',F0,NU(NCF),NU(1),
	3     '; outside CMF grid, skipped'
	    CYCLE
	  END IF
	  LO=1
	  HI=NCF
	  DO WHILE(HI-LO .GT. 1)
	    MID=(LO+HI)/2
	    IF(NU(MID) .GE. F0)THEN
	      LO=MID
	    ELSE
	      HI=MID
	    END IF
	  END DO
	  K0=LO
!
! Expand to all cell centers inside F0 +/- W.
!
	  KA=K0
	  DO WHILE(KA .GT. 1 .AND. NU(KA-1) .LE. F0+W)
	    KA=KA-1
	  END DO
	  KB=K0+1
	  DO WHILE(KB .LT. NCF .AND. NU(KB+1) .GE. F0-W)
	    KB=KB+1
	  END DO
!
! Normalize the Gaussian by the local cell widths so the stored
! profile integrates to unity over frequency.  NU is in units of
! 10^15 Hz while CHIL_MAT/ETAL_MAT are per-Hz quantities (cf. the
! 1.0D-15 factor in the LINE_PROF_SIM normalization in SET_PROF_V6),
! so the weight carries an explicit 1.0D-15.
!
! The Gaussian-only path (DO_SOB_DEP_WINGS=F) shares one profile for
! all depths, so XNORM is computed once here.  The Voigt-wing path
! carries a depth-dependent damping parameter, so its normalization
! is done per depth inside the loop below.
!
	  XNORM=0.0_LDP
	  IF(.NOT. DO_SOB_DEP_WINGS)THEN
	    DO K=KA,KB
	      X=(NU(K)-F0)/DNU_D
	      WT=EXP(-X*X)
	      IF(K .EQ. 1)THEN
	        DN_K=NU(1)-NU(2)
	      ELSE IF(K .EQ. NCF)THEN
	        DN_K=NU(NCF-1)-NU(NCF)
	      ELSE
	        DN_K=0.5_LDP*(NU(K-1)-NU(K+1))
	      END IF
	      XNORM=XNORM+WT*DN_K
	    END DO
	    IF(XNORM .LE. 0.0_LDP)CYCLE
	  END IF
!
! Per-line collisional damping constant, same expression as
! SET_PROF_V6: LOC_GAM_COL = 1.55E4*|C4|^0.667, scaled by local ED.
!
	  IF(DO_SOB_DEP_WINGS)
	1   GAM_COL_LINE=1.55E+04_LDP*(ABS(C4_VEC(LP))**0.667_LDP)
!
	  N_DROP_L=0
	  DO I=1,ND
!
! NEG_OPACITY cells had their CHIL_MAT crushed to 1/NUM_SIM_LINES by
! the blend-group maser check above; the surviving ETAL_MAT would
! otherwise be injected as thin emission with an artificial opacity.
!
	    IF(CHIL_MAT(I,S) .LE. 0.0_LDP .OR. NEG_OPACITY(I))THEN
	      N_DROP=N_DROP+1
	      N_DROP_L=N_DROP_L+1
	      IF(.NOT. DO_SOB_DEP_FLOOR)THEN
	        DROP_SUM=DROP_SUM+ETAL_MAT(I,S)
	        CYCLE
	      END IF
	    END IF
	    INJ_SUM=INJ_SUM+ETAL_MAT(I,S)
!
! DO_SOB_DEP_WINGS: true Voigt profile per depth with the same
! A_VOIGT formula the stock CMF profiles use (SET_PROF_V6), so the
! deposited wings match the damping width the blanketing solve would
! have carried.  Normalize per depth since A_VOIGT depends on ED(I).
!
	    IF(DO_SOB_DEP_WINGS)THEN
	      TMP_ED=MIN(ED_IN(I),MAX_PROF_ED)
	      A_VOIGT=0.25E-15_LDP*(ARAD_VEC(LP)+
	1              GAM_COL_LINE*TMP_ED)/PI/DNU_D
	      XNORM=0.0_LDP
	      DO K=KA,KB
	        X=(NU(K)-F0)/DNU_D
	        WT=VOIGT(A_VOIGT,X)
	        IF(K .EQ. 1)THEN
	          DN_K=NU(1)-NU(2)
	        ELSE IF(K .EQ. NCF)THEN
	          DN_K=NU(NCF-1)-NU(NCF)
	        ELSE
	          DN_K=0.5_LDP*(NU(K-1)-NU(K+1))
	        END IF
	        XNORM=XNORM+WT*DN_K
	      END DO
	      IF(XNORM .LE. 0.0_LDP)CYCLE
	    END IF
	    DO K=KA,KB
	      X=(NU(K)-F0)/DNU_D
	      IF(DO_SOB_DEP_WINGS)THEN
	        WT=VOIGT(A_VOIGT,X)
	      ELSE
	        WT=EXP(-X*X)
	      END IF
	      WT=1.0E-15_LDP*WT/XNORM
!
! DO_SOB_DEP_FLOOR: masing cells (CHIL_MAT<=0 or the blend-group
! NEG_OPACITY flag) carry real emissivity but the formal integral
! cannot amplify; deposit them as optically thin emitters by
! substituting SOB_DEP_CHI_FLOOR for the negative opacity.
!
	      IF(DO_SOB_DEP_FLOOR .AND.
	1        (CHIL_MAT(I,S) .LE. 0.0_LDP .OR. NEG_OPACITY(I)))THEN
	        T1=SOB_DEP_CHI_FLOOR
	      ELSE
	        T1=CHIL_MAT(I,S)
	      END IF
	      CHI_ST(I,K)=CHI_ST(I,K)+
     1	          T1*WT
	      ETA_ST(I,K)=ETA_ST(I,K)+
     1	          ETAL_MAT(I,S)*WT
	    END DO
	  END DO
!
! A line masing at more than 10% of the depths renders as a lower
! bound only; name it so the run log shows which transition lost
! flux, not just the group totals.
!
	  IF(N_DROP_L*10 .GT. ND)THEN
	    WRITE(LUER,'(A,I6,A,I6,A,I6,A,ES16.7)')
	1     'ADD_SOB_LINE_PROF: line',LP,' dropped',N_DROP_L,
	2     ' of',ND,' depths (masing/neg-opacity), F0=',F0
	  END IF
	END DO
!
	RETURN
	END
