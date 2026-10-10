!
! Subroutine to increment the EHB variation matrix BA due to the variation
! of J associated with free-free radiation and level populations.
!
! The variation equation is not updated for LAMBDA iterations.
!
! BA_T_PAR_EHB is used for the diagonal variation only. It is updated on each call
! rather than BA to improve numerical stability. BA_PAR should contain terms
! of similar size. BA_T_PAR_EHB will need to be added to BA_T_PAR after every
! every N frequencies. In this way the matrices should suffer less cancelation
! effects due to the addition of large positive and negative ! terms.
!
! Utilizesg the fact that consecutive frequency terms should be correlated, and
! hence similar in size. While some minor cancellation, should be much less
! then adding to full BA matrix in which terms have arbitrary size.
!
	SUBROUTINE BA_EHB_FF_UPDATE_V1(VJ,VCHI_FF,VETA_FF,
	1              ETA_FF,CHI_FF,T,POPS,RJ,
	1              FQW,NEW_CONT,FINAL_FREQ,DO_SRCE_VAR_ONLY,
	1              NION,NT,NUM_BNDS,ND,DST,DEND)
	USE SET_KIND_MODULE
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Created 14-Jul-2019
!
	REAL(KIND=LDP), SAVE, ALLOCATABLE ::  VJ_T(:,:,:)
	INTEGER, SAVE :: CNT=0
!
	INTEGER NION
  	INTEGER NT,NUM_BNDS,ND,DST,DEND
	REAL(KIND=LDP) VJ(NT,NUM_BNDS,ND)
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) VCHI_FF(NT,ND)
	REAL(KIND=LDP) VETA_FF(NT,ND)
	REAL(KIND=LDP) RJ(ND)
!
	REAL(KIND=LDP) ETA_FF(ND)
	REAL(KIND=LDP) CHI_FF(ND)
	REAL(KIND=LDP) T(ND)
!
	REAL(KIND=LDP) NU
	REAL(KIND=LDP) FQW
!
! NEW_CONT indicates that this is the first frequency of a new continuum band
! in which the continuum cross-sections are constant. FINAL_FREQ indicates
! that it is the last frequency of a continuum band.
!
	LOGICAL FINAL_FREQ
	LOGICAL NEW_CONT
	LOGICAL DO_SRCE_VAR_ONLY
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL(KIND=LDP) FOUR_PI
	REAL(KIND=LDP) T1,T2,QFV_T
	REAL(KIND=LDP) RJ_RAD(ND)
	INTEGER I,J,K,L,LS,IOS,JJ,ID
	INTEGER DIAG_INDX,BNDST,BNDEND
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
	CALL TUNE(IONE,'BA_EHB_UP')
	DIAG_INDX=(NUM_BNDS+1)/2
!
	FOUR_PI=16.0_LDP*ATAN(1.0_LDP)
	IF(.NOT. ALLOCATED(VJ_T))THEN
	  ALLOCATE (VJ_T(NT,NUM_BNDS,ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in BA_EHB_UPDATE_V1'
	    WRITE(I,*)'Unable to allocate required dynamic memory'
	    STOP
	  END IF
	END IF
!
	IF(DO_SRCE_VAR_ONLY)THEN
	  DO L=DST,DEND
	    RJ_RAD(L)=ETA_FF(L)/CHI_FF(L)
	  END DO
	ELSE
	  DO L=DST,DEND
	    RJ_RAD(L)=RJ(L)
	  END DO
	END IF
!
! Perform the frequency integral of dJ over. Procedure depends on whether
! this is a new frequency of part of a band.
!
	T1=1.0E-10_LDP*FOUR_PI*FQW
	IF(NEW_CONT)THEN
!$OMP PARALLEL DO PRIVATE(L,K)
	  DO L=DST,DEND
	    DO K=1,NUM_BNDS
	      IF(K .EQ. DIAG_INDX)THEN
	        VJ_T(:,K,L)=T1*( RJ_RAD(L)*VCHI_FF(:,L) - VETA_FF(:,L) + CHI_FF(L)*VJ(:,K,L) )
	      ELSE
	        VJ_T(:,K,L)=T1*CHI_FF(L)*VJ(:,K,L)
	      END IF
	    END DO
	  END DO
!$OMP END PARALLEL DO
	ELSE
!$OMP PARALLEL DO PRIVATE(I,L,K)
	  DO L=1,ND                        !DST,DEND
	    DO K=1,NUM_BNDS
	      IF(K .EQ. DIAG_INDX)THEN
	        DO I=1,NT
	           VJ_T(I,K,L)= VJ_T(I,K,L) +
	1              T1*( RJ_RAD(L)*VCHI_FF(I,L) - VETA_FF(I,L)+ CHI_FF(L)*VJ(I,K,L) )
	        END DO
	      ELSE
	        DO I=1,NT
	           VJ_T(I,K,L)= VJ_T(I,K,L) + T1*CHI_FF(L)*VJ(I,K,L)
	        END DO
	      END IF
	    END DO
	  END DO
!$OMP END PARALLEL DO
	END IF
!
!	IF(CNT .EQ. 0)WRITE(296,*)NION,NT,ND,DST,DEND
!	CNT=CNT+1
!	WRITE(296,'(I8,11ES14.4)')CNT,ETA_FF(5),CHI_FF(5),RJ_RAD(5),VCHI_FF(NT-1,5),
!	1         VETA_FF(NT-1,5),RJ_RAD(5)*VCHI_FF(NT-1,5),CHI_FF(5)*VJ(NT-1,DIAG_INDX,5),
!	1              T1*( RJ_RAD(5)*VCHI_FF(NT-1,5) - VETA_FF(NT-1,5)+ CHI_FF(5)*VJ(NT-1,DIAG_INDX,5) ),
!	1               VJ_T(NT-1,DIAG_INDX,5),T1*(CHI_FF(5)*RJ_RAD(5)-ETA_FF(5)),T1
!
	IF(FINAL_FREQ)THEN
!
!$OMP PARALLEL DO PRIVATE(L,J,K,BNDST,BNDEND)
!
	  DO L=DST,DEND					!S.E. equation depth
	    BNDST=MAX( 1+DIAG_INDX-L, 1 )
	    BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	    DO K=BNDST,BNDEND	  			!Variable depth.
	      IF(K .EQ. DIAG_INDX)THEN
   	        DO  J=1,NT	  	  	  	!Variable
	          BA_T_PAR_EHB(J,L)=BA_T_PAR_EHB(J,L) + VJ_T(J,K,L)
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	  	!Variable
	          BA_T_EHB(J,K,L)=BA_T_EHB(J,K,L) + VJ_T(J,K,L)
	        END DO
	      END IF
	    END DO
	  END DO
!$OMP END PARALLEL DO
!
	END IF
	CALL TUNE(ITWO,'BA_EHB_UP')
!
	RETURN
	END
