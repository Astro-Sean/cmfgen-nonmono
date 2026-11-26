! This will be a subroutine that will be used to compute the normalization constants left out of the scattering redistribution.
! This code will find the normalization by conserving photons instead of energy, since energy is lost from compton scattering.
!
! Version 2 is using a new method for calculating the deposition
!
! Created Jan. 23, 2015
!
! v4 now has the indicies for ETA and INTENSITY switched from (mu,nu,depth)->(mu,depth,nu)
!------------------------------------------------------------------------------------------------------------------------------
!
	SUBROUTINE GAMMA_ENERGY_DEP_V7(INTENSITY,ETA,CHI_SCAT,CHI_ABS,&
			E_DEP,DECAY_KIN_E,NU_GRID,NF,V,R,NA,ND,GAMRAY_EMISS,SN_AGE_DAYS)
	USE SET_KIND_MODULE
	USE GAM_MU_MOD
	IMPLICIT NONE
!
! Altered 21-Nov-2021 : Modied subroutine call to be consistent with file name.
! Altered 19-Nov-2021 : Added SN_AGE_DAYS to call (version not changed).
!
	INTEGER :: I,J,K,L
	INTEGER :: NF
	INTEGER :: ND
	INTEGER :: NA
!
	REAL(KIND=LDP) :: E_DEP(ND)
	REAL(KIND=LDP) :: DECAY_KIN_E(ND)
	REAL(KIND=LDP) :: V(ND)
	REAL(KIND=LDP) :: R(ND)
	REAL(KIND=LDP) :: NU_GRID(NF)
	REAL(KIND=LDP) :: CHI_SCAT(NF,ND)
	REAL(KIND=LDP) :: CHI_ABS(NF,ND)
	REAL(KIND=LDP) :: CHI_TOT(NF,ND)
	REAL(KIND=LDP) :: GAMRAY_EMISS(ND)
!
	REAL(KIND=LDP) :: INTENSITY(NA,ND,NF)
	REAL(KIND=LDP) :: ETA(NA,ND,NF)
!
	REAL(KIND=LDP) :: TEMP_INT1(NF,ND)
	REAL(KIND=LDP) :: TEMP_INT2(ND)
	REAL(KIND=LDP) :: TEMP_INT3(ND)
	REAL(KIND=LDP) :: TEMP_ETA1(NF,ND)
	REAL(KIND=LDP) :: TEMP_ETA2(ND)
	REAL(KIND=LDP) :: TEMP_FLUX(NF)
!
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
	REAL(KIND=LDP), ALLOCATABLE :: TB(:)
	REAL(KIND=LDP), ALLOCATABLE :: TC(:)
	REAL(KIND=LDP), ALLOCATABLE :: TD(:)
!
	REAL(KIND=LDP) :: FLUX
	REAL(KIND=LDP) :: T1,T2,T3,T4
	REAL(KIND=LDP) :: PI
	REAL(KIND=LDP) :: FOURPI
	REAL(KIND=LDP) :: SN_AGE_DAYS
!
	WRITE(6,*)'Using gamray_energy_deposit_v7.f90'
	PI=ACOS(-1.0_LDP)
	FOURPI=4.0_LDP*PI
	TEMP_INT1=0.0_LDP
	TEMP_INT2=0.0_LDP
	TEMP_INT3=0.0_LDP
	TEMP_ETA1=0.0_LDP
	TEMP_ETA2=0.0_LDP
	OPEN(UNIT=7,FILE='./data/DIAGN_EDEP',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
	WRITE(7,'(5(A17,3X))')'Depth:','SCAT. INTEGRAL:','ETA INTEGRAL:','ABS. INTEGRAL:','DECAY_KIN_E:'
	DO I=1,ND
	  L=R_MU(I)%MU_PTS
	  ALLOCATE(TA(L),TB(L))
	  DO J=1,NF
	    DO K=1,L
	      TA(K)=INTENSITY(K,I,J)
	      TB(K)=ETA(K,I,J)
	    END DO
	    CALL LUM_FROM_ETA(TA,R_MU(I)%MU_VECTOR,L)
	    CALL LUM_FROM_ETA(TB,R_MU(I)%MU_VECTOR,L)
	    TEMP_INT1(J,I)=SUM(TA)
	    TEMP_ETA1(J,I)=SUM(TB)
	    TEMP_ETA1(J,I)=TEMP_ETA1(J,I)*2.0_LDP*PI !+4.0D0*PI*ETA_ISO(I,J) !addding it outside angle loop since isotropic
	    TEMP_INT1(J,I)=TEMP_INT1(J,I)*2.0_LDP*PI
	  END DO
!
! I am combining the opacities in case there is numerical inaccuracy by
! treating it as two terms in the integration. I would effectively
! double the error for the same quadrature.
!
	  ALLOCATE(TC(NF),TD(NF))
	  CHI_TOT(1:NF,I)=CHI_SCAT(1:NF,I)+CHI_ABS(1:NF,I)
	  DO J=1,NF
	    TC(J)=CHI_TOT(J,I)*TEMP_INT1(J,I)
	    TD(J)=TEMP_ETA1(J,I)
	  END DO
	  CALL LUM_FROM_ETA(TC,NU_GRID,NF)
	  CALL LUM_FROM_ETA(TD,NU_GRID,NF)
	  TEMP_INT2(I)=SUM(TC)
	  TEMP_ETA2(I)=SUM(TD)
	  WRITE(7,'(14X,I3,3X,ES16.6,4X,ES16.6,4X,ES16.6)')I,TEMP_INT2(I),TEMP_ETA2(I),&
			DECAY_KIN_E(I)
	  E_DEP(I)=TEMP_INT2(I)-TEMP_ETA2(I)+DECAY_KIN_E(I)
	  DEALLOCATE(TA,TB,TC,TD)
	END DO
	CLOSE(UNIT=7)
!
	ALLOCATE(TA(ND),TB(ND),TC(ND))
	TA=0.0_LDP
	TB=0.0_LDP
	DO I=1,ND
	  TA(I)=E_DEP(I)*R(I)*R(I)
	  TB(I)=GAMRAY_EMISS(I)*R(I)*R(I)
	  TC(I)=DECAY_KIN_E(I)*R(I)*R(I)
	END DO
	CALL LUM_FROM_ETA(TA,R,ND);   T1=SUM(TA)*FOURPI*1.0E+30_LDP
	CALL LUM_FROM_ETA(TB,R,ND);   T2=SUM(TB)*FOURPI*1.0E+30_LDP
	CALL LUM_FROM_ETA(TC,R,ND);   T3=SUM(TC)*FOURPI*1.0E+30_LDP
!
	OPEN(UNIT=7,FILE='GAMRAY_ENERGY_DEP',STATUS='UNKNOWN',ACTION='WRITE')
	WRITE(7,'(A)')'!'
	WRITE(7,'(A)')'! Energy deposited from scattering and photoelectric abs. ergs/s/cm^3'
	WRITE(7,'(A)')'!'
	WRITE(7,'(A,1X,ES16.7,4X,ES12.4,A)')'!                 Total Energy Emitted:',T2,T2/3.8286D+33,' Lsun'
	WRITE(7,'(A,1X,ES16.7,4X,ES12.4,A)')'!                Total Energy Absorbed:',T1,T1/3.828D+33,' Lsun'
	WRITE(7,'(A,1X,ES16.7)')'!          Fraction of energy absorbed:',T1/T2
	WRITE(7,'(A,1X,ES16.7)')'! Fraction absorbed that is kinetic is:',T3/T1
	WRITE(7,'(A)')'!'
	WRITE(7,'(I5,T30,A)')ND,'!Number of depth points'
	WRITE(7,'(ES14.8,T30,A)')SN_AGE_DAYS,'!Current time after explosion '
	WRITE(7,'(A)')'!'
	WRITE(7,'(A1,1X,5(A20,2X))')'!','Radius (10^10 cm)','V (km/s)','Energy dep.',&
		'Local emission','Decay K.E.'
	DO I=1,ND
	  WRITE(7,'(T3,5(ES20.10,2X))')R(I),V(I),E_DEP(I),GAMRAY_EMISS(I),DECAY_KIN_E(I)
	END DO
	CLOSE(UNIT=7)
!
	RETURN
	END SUBROUTINE
