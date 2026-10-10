	MODULE BA_EHB_MOD_V1
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Module which contains dJdN integrated over a small frequency band over
! which the continuum cross-sections are assumed not to change.
!
! VJ_R = Int ( dJ*EXP(_hv/kT)/v) dv
! VJ_P = Int ( dJ/v)
! VJ_C = Int ( dJ)
! RJ_SUM = Int (J)
!
	REAL(KIND=LDP), ALLOCATABLE :: VJ_R(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: VJ_P(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: VJ_R_EHB(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: VJ_P_EHB(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RJ_SUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: RJ_SUM_EHB(:)
	SAVE
!
	END MODULE BA_EHB_MOD_V1
!
! Subroutine to increment the variation matrix BA_T_EHB and BA_T_PAR_EHBR
! due to the variation of J in the bound-free heating terms..
!
! BA_T_PAR_EHB is used for the diagonal variation only. It is updated on each call
! rather than BA_T_EHB to improve numerical stability. BA_PAR should contain terms
! of similar size. BA_T_PAR_EHB  will need to be added to BA_T_EHB after every approximately
! every N frequencies. In this way the matrices should suffer less cancelation
! effects due to the addition of large positive and negative terms.
!
!
	SUBROUTINE BA_EHB_BF_UPDATE_V1(VJ,ETA_CONT,CHI_CONT,POPS,RJ,
	1              NU,FQW,NEW_CONT,FINAL_FREQ,DO_SRCE_VAR_ONLY,
	1              NION,NT,NUM_BNDS,ND,DST,DEND)
	USE SET_KIND_MODULE
	USE BA_EHB_MOD_V1
	USE STEQ_DATA_MOD
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 14-Jul-2019
!
	INTEGER NION
  	INTEGER NT,NUM_BNDS,ND,DST,DEND
	REAL(KIND=LDP) VJ(NT,NUM_BNDS,ND)
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) VCHI(NT,ND)
	REAL(KIND=LDP) VETA(NT,ND)
	REAL(KIND=LDP) RJ(ND)
!
	REAL(KIND=LDP) ETA_CONT(ND)
	REAL(KIND=LDP) CHI_CONT(ND)
	REAL(KIND=LDP) ESEC(ND)
!
	REAL(KIND=LDP) NU
	REAL(KIND=LDP) FQW
	REAL(KIND=LDP) STORE(ND)
!
! NEW_CONT indicates that this is the first frequency of a new continuum band
! in which the continuum cross-sections are constant. FINAL_FREQ indicates
! that it is the last frequency of a continuum band.
!
	LOGICAL FINAL_FREQ
	LOGICAL NEW_CONT
	LOGICAL DO_SRCE_VAR_ONLY
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL(KIND=LDP) PLANCKS_CONSTANT,PC
	EXTERNAL PLANCKS_CONSTANT
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) RJ_RAD(ND)
	REAL(KIND=LDP) LOC_QFV_P,LOC_QFV_R
	REAL(KIND=LDP) LOC_QFV_P_EHB,LOC_QFV_R_EHB
	INTEGER I,J,K,L,LS,IOS,ID
	INTEGER DIAG_INDX,BNDST,BNDEND
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
	CALL TUNE(IONE,'BA_EHB_T_UP')
	DIAG_INDX=(NUM_BNDS+1)/2
!
	IF(.NOT. ALLOCATED(VJ_R))THEN
	  ALLOCATE (VJ_R(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VJ_P(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VJ_R_EHB(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VJ_P_EHB(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RJ_SUM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RJ_SUM_EHB(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in BE_EHB_BF_UPDATE_V1'
	    WRITE(I,*)'Unable to allocate required dynamic memory'
	    STOP
	  END IF
	END IF
	PC=1.0E+15_LDP*PLANCKS_CONSTANT()
!
	IF(DO_SRCE_VAR_ONLY)THEN
	  DO L=DST,DEND
	    RJ_RAD(L)=ETA_CONT(L)/(CHI_CONT(L)-ESEC(L))
	  END DO
	ELSE
	  DO L=DST,DEND
	    RJ_RAD(L)=RJ(L)
	  END DO
	END IF
!
! Perform the frequency integral of dJ over. Procedure depends on whether
! this is a new frequency of part of a band. In oder to minimize computation,
! a band which is a single frequency (i.e. NEW_CONT and FINAL_FREQ both TRUE)
! is treated as a special case.
!
	IF(NEW_CONT)THEN
	  T2=FQW/NU
!$OMP PARALLEL DO PRIVATE(T1)
	  DO L=DST,DEND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    VJ_R(:,:,L)=T1*VJ(:,:,L)
	    VJ_P(:,:,L)=T2*VJ(:,:,L)
	    VJ_R_EHB(:,:,L)=T1*NU*VJ(:,:,L)
	    VJ_P_EHB(:,:,L)=T2*NU*VJ(:,:,L)
	    RJ_SUM(L)=T2*RJ(L)
	    RJ_SUM_EHB(L)=T2*NU*RJ(L)
	  END DO
!$OMP END PARALLEL DO
	ELSE
	  T2=FQW/NU
!$OMP PARALLEL DO PRIVATE(T1)
	  DO L=DST,DEND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    VJ_R(:,:,L)=VJ_R(:,:,L)+T1*VJ(:,:,L)
	    VJ_P(:,:,L)=VJ_P(:,:,L)+T2*VJ(:,:,L)
	    VJ_R_EHB(:,:,L)=VJ_R_EHB(:,:,L)+T1*NU*VJ(:,:,L)
	    VJ_P_EHB(:,:,L)=VJ_P_EHB(:,:,L)+T2*NU*VJ(:,:,L)
	    RJ_SUM(L)=RJ_SUM(L)+T2*RJ(L)
	    RJ_SUM_EHB(L)=RJ_SUM_EHB(L)+T2*NU*RJ(L)
	  END DO
!$OMP END PARALLEL DO
	END IF
	CALL TUNE(2,'BA_EHB_T_UP')
!
	IF(FINAL_FREQ)THEN
	  CALL TUNE(1,'BA_EHB_T_FF')
	  STORE=BA_T_PAR_EHB(NT,1:ND)
	  DO ID=1,NION
            IF(SE(ID)%XzV_PRES)THEN
!
!$OMP PARALLEL DO PRIVATE(T1,T2,L,K,J,BNDST,BNDEND,LOC_QFV_P,LOC_QFV_R,LOC_QFV_P_EHB,LOC_QFV_R_EHB)
	      DO L=DST,DEND
	        BNDST=MAX( 1+DIAG_INDX-L, 1 )
	        BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	        LOC_QFV_P=PC*SUM(SE(ID)%QFV_P(:,L))
	        LOC_QFV_R=PC*SUM(SE(ID)%QFV_R(:,L))
	        LOC_QFV_P_EHB=PC*SUM(SE(ID)%QFV_P_EHB(:,L))
		LOC_QFV_R_EHB=PC*SUM(SE(ID)%QFV_R_EHB(:,L))
!
	        DO K=BNDST,BNDEND
                  DO J=1,NT
	            T1=LOC_QFV_P*VJ_P_EHB(J,K,L)-LOC_QFV_R*VJ_R_EHB(J,K,L)
	            T2=LOC_QFV_P_EHB*VJ_P(J,K,L)-LOC_QFV_R_EHB*VJ_R(J,K,L)
	            IF(K .EQ. DIAG_INDX)THEN
	              BA_T_PAR_EHB(J,L)=BA_T_PAR_EHB(J,L)+(T1+T2)
	            ELSE
	              BA_T_EHB(J,K,L)=BA_T_EHB(J,K,L)+(T1+T2)
	            END IF
	          END DO
	        END DO
	      END DO
!$OMP END PARALLEL DO
	    END IF
	  END DO
	   CALL TUNE(2,'BA_EHB_T_FF')
	END IF
!
	RETURN
	END
