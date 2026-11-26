!
! Subroutine to reset population after the velocity has been adjusted.
! It is assumed the POP_SPEC and POP_ATOM have been reset.
!
! This routine was written assuming the radisu grid (but not the
! Boundary R values) have changed.
!
	SUBROUTINE WIND_SCALE_POPS_V1(POPS,R_OLD,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created: 29-Jan-2021 : PArtially based on SET_NEW_MODEL_ESTIMATES
!
	INTEGER ND
	INTEGER NT
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) R_OLD(ND)
	REAL(KIND=LDP) Z_POP(NT)
	LOGICAL DO_LEV_DISSOLUTION
!
! Local variables and vectors
!
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP) LOG_R_OLD(ND)
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) GAM(ND)
	REAL(KIND=LDP) OLD_POPATOM(ND)
	REAL(KIND=LDP) T1
!
	INTEGER ID
	INTEGER ISPEC
	INTEGER I,J
	INTEGER, PARAMETER :: IONE=1
!
	LOGICAL FIRST
!
	LOG_R=LOG(R)
	LOG_R_OLD=LOG(R_OLD)
!
! Compute revised electron density. We assume the fractioanl ioization
! is constant.
!
! We only sum to NT-2 as NT-1 is Ne and NT is T.
!
	DO I=1,ND
	  T1=0.0_LDP
	  DO J=1,NT-2
	    T1=T1+POPS(J,I)
	  END DO
	  OLD_POPATOM(I)=T1
	  GAM(I)=ED(I)/T1
	END DO
!
	CALL MON_INTERP(ED,ND,IONE,LOG_R,ND,GAM,ND,LOG_R_OLD,ND)
	ED=ED*POP_ATOM
!
	TA=T
	CALL MON_INTERP(T,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
!
! Compute vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common
! block, and are required by SUP_TO_FULL and LTE_POP_WLD.
!
! As a first estimate of POPION, we simply assume the fractional ioization is
! constant.
!
	DO J=1,ND
	  TA(J)=0.0_LDP
	  DO I=1,NT-2
	    IF(Z_POP(I) .GT. 0.01_LDP)TA(J)=TA(J)+POPS(I,J)
	  END DO
	END DO
        DO I=1,ND
          POPION(I)=POP_ATOM(I)*TA(I)/OLD_POPATOM(I)
        END DO
        CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! Compute LOG(DC)
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
            IF(ATM(ID)%XzV_PRES)ATM(ID)%XZV_F=LOG(ATM(ID)%XzV_F)-ATM(ID)%LOG_XzVLTE_F
	  END DO
	END DO
!
! Interpolate onto the new R grid.
!
	DO ISPEC=1,NUM_SPECIES
!
	  IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	    ID=SPECIES_END_ID(ISPEC)-1
	    TA=LOG(ATM(ID)%DxZV_F)
	    CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	    ATM(ID)%DXzV_F=EXP(TB)
	  END IF
!
	  DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO I=1,ATM(ID)%NXzV_F
	         TA(:)=ATM(ID)%XzV_F(I,:)
	         CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
!	        ATM(ID)%XzV_F(I,:)=EXP(TB(:))
	      END DO
	    END IF
	  END DO
!
	  FIRST=.TRUE.
	  DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL LTEPOP_WLD_V2(ATM(ID)%XzVLTE_F, ATM(ID)%LOG_XzVLTE_F,ATM(ID)%W_XzV_F,
	1              ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,
	1              ATM(ID)%ZXzV,      ATM(ID)%GIONXzV_F,
	1              ATM(ID)%NXzV_F,    ATM(ID)%DXzV_F,     ED,T,ND)
	      CALL CNVT_FR_DC_V2(ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
	1              ATM(ID)%DXzV_F,    ATM(ID)%NXzV_F,
	1              TB,                TA,ND,
	1              FIRST,             ATM(ID+1)%XzV_PRES)
              IF(ID .NE.  SPECIES_BEG_ID(ISPEC))ATM(ID-1)%DXzV_F(1:ND)=TB(1:ND)
   	    END IF
	  END DO
!
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    IF(ATM(ID)%XzV_PRES)THEN
	       CALL SCALE_POPS(ATM(ID)%XzV_F,ATM(ID)%DXzV_F,
	1               POP_SPECIES(1,ISPEC),TA,ATM(ID)%NXzV_F,ND)
	    END IF
	  END DO
	END DO
!
	DO ID=NUM_IONS-1,1,-1
          CALL FULL_TO_SUP(
	1      ATM(ID)%XzV,   ATM(ID)%NXzV,       ATM(ID)%DXzV, ATM(ID)%XzV_PRES,
	1      ATM(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,
	1      ATM(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, ND)
        END DO
!
! Store all quantities in POPS array. This is done here as it enables POPION
! to be readily computed. It also ensures that POPS is correct if we  don't
! iterate on T.
!
	DO ID=1,NUM_IONS-1
	  CALL IONTOPOP(POPS,  ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1          ATM(ID)%EQXzV, ATM(ID)%NXzV, NT,ND, ATM(ID)%XzV_PRES)
	END DO
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	DO J=1,ND
	  POPION(J)=0.0_LDP
	  DO I=1,NT
	    IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	  END DO
	END DO
!
! Evaluates LTE populations for both the FULL atom, and super levels.
!
	CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
!
	RETURN
	END
