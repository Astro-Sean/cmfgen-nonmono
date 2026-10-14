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
!
! The profile is a normalized Doppler core of width
! FL_SIM*VDOP_MIN/C_KMS evaluated on the (descending) CMF grid NU.
! Only grid cells inside +/-3 Doppler widths carry weight, and the
! weights are normalized by the local cell width so the deposited
! opacity integrated over frequency equals CHIL_MAT*CLUMP_FAC
! independent of how coarsely the CMF grid resolves the line.  The
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
! dropped clumped emissivity for the caller's diagnostic print.
! ETAL_MAT is non-negative by construction.
!
	SUBROUTINE ADD_SOB_LINE_PROF(ETA_ST,CHI_ST,NU,NCF,
     1	      CHIL_MAT,ETAL_MAT,CLUMP_FAC,FL_SIM,SIM_POINTER,
     2	      VDOP_MIN,TRANS_TYPE,NEG_OPACITY,
     3	      NUM_SIM,N_LINES,MAX_SIM,ND,C_KMS,
     4	      N_DROP,DROP_SUM,INJ_SUM)
	USE SET_KIND_MODULE
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
	REAL(KIND=LDP), INTENT(IN) :: C_KMS
	LOGICAL, INTENT(IN) :: NEG_OPACITY(ND)
	INTEGER, INTENT(IN) :: SIM_POINTER(MAX_SIM)
	CHARACTER(LEN=*), INTENT(IN) :: TRANS_TYPE(N_LINES)
!
	INTEGER I,K,S,LP,K0,KA,KB,LO,HI,MID
	REAL(KIND=LDP) F0,DNU_D,W,XNORM,WT,DN_K
!
	DO S=1,NUM_SIM
	  LP=SIM_POINTER(S)
	  IF(TRANS_TYPE(LP)(1:3) .NE. 'SOB')CYCLE
	  F0=FL_SIM(S)
	  DNU_D=F0*VDOP_MIN(LP)/C_KMS
	  IF(DNU_D .LE. 0.0_LDP)CYCLE
	  W=3.0_LDP*DNU_D
!
! NU decreases with index.  K0 is the last index with NU(K0) >= F0,
! so F0 lies in [NU(K0+1),NU(K0)).  Lines outside the grid are skipped.
!
	  IF(NU(1) .LT. F0)CYCLE
	  IF(NU(NCF) .GE. F0)CYCLE
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
	  XNORM=0.0_LDP
	  DO K=KA,KB
	    WT=EXP(-((NU(K)-F0)/DNU_D)**2)
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
!
	  DO I=1,ND
!
! NEG_OPACITY cells had their CHIL_MAT crushed to 1/NUM_SIM_LINES by
! the blend-group maser check above; the surviving ETAL_MAT would
! otherwise be injected as thin emission with an artificial opacity.
!
	    IF(CHIL_MAT(I,S) .LE. 0.0_LDP .OR. NEG_OPACITY(I))THEN
	      N_DROP=N_DROP+1
	      DROP_SUM=DROP_SUM+ETAL_MAT(I,S)*CLUMP_FAC(I)
	      CYCLE
	    END IF
	    INJ_SUM=INJ_SUM+ETAL_MAT(I,S)*CLUMP_FAC(I)
	    DO K=KA,KB
	      WT=1.0E-15_LDP*EXP(-((NU(K)-F0)/DNU_D)**2)/XNORM
	      CHI_ST(I,K)=CHI_ST(I,K)+
     1	          CHIL_MAT(I,S)*CLUMP_FAC(I)*WT
	      ETA_ST(I,K)=ETA_ST(I,K)+
     1	          ETAL_MAT(I,S)*CLUMP_FAC(I)*WT
	    END DO
	  END DO
	END DO
!
	RETURN
	END
