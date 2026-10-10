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
	SUBROUTINE GAMRAY_ENERGY_DEPOSIT_V5(INTENSITY,ETA,CHI_SCAT,CHI_ABS,&
	USE SET_KIND_MODULE
			E_DEP,DECAY_KIN_E,NU_GRID,NF,V,R,NA,ND,GAMRAY_EMISS)
	USE GAM_MU_MOD
!
	IMPLICIT NONE
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
	REAL(KIND=LDP) :: FLUX
	REAL(KIND=LDP) :: T1,T2,T3,T4
	REAL(KIND=LDP) :: PI
	REAL(KIND=LDP) :: FOURPI
!
	WRITE(6,*)'Using GAMRAY_ENERGY_DEPOSIT_V5'
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
	  DO J=1,NF
	    DO K=1,L-1
	      T1=R_MU(I)%MU_VECTOR(K) 	!Remember, the grid goes from 1 to -1
	      T2=R_MU(I)%MU_VECTOR(K+1)
	      TEMP_INT1(J,I)=TEMP_INT1(J,I)+0.5_LDP*(T1-T2)*(INTENSITY(K,I,J)+INTENSITY(K+1,I,J))
	      TEMP_ETA1(J,I)=TEMP_ETA1(J,I)+0.5_LDP*(T1-T2)*(ETA(K,I,J)+ETA(K+1,I,J))
	    END DO
	    TEMP_ETA1(J,I)=TEMP_ETA1(J,I)*2.0_LDP*PI !+4.0D0*PI*ETA_ISO(I,J) !addding it outside angle loop since isotropic
	    TEMP_INT1(J,I)=TEMP_INT1(J,I)*2.0_LDP*PI
	  END DO
	  DO J=1,NF-1
	    T1=NU_GRID(J)	!Remember, NU_GRID is decreasing from hi energy to low energy
	    T2=NU_GRID(J+1)
	    TEMP_INT2(I)=TEMP_INT2(I)+0.5_LDP*(T1-T2)*(CHI_SCAT(J,I)*TEMP_INT1(J,I)&
				+CHI_SCAT(J+1,I)*TEMP_INT1(J+1,I))
	    TEMP_INT3(I)=TEMP_INT3(I)+0.5_LDP*(T1-T2)*(CHI_ABS(J,I)*TEMP_INT1(J,I)&
				+CHI_ABS(J+1,I)*TEMP_INT1(J+1,I))
	    TEMP_ETA2(I)=TEMP_ETA2(I)+0.5_LDP*(T1-T2)*(TEMP_ETA1(J,I)+TEMP_ETA1(J+1,I))
	  END DO
	  WRITE(7,'(14X,I3,3X,ES16.6,4X,ES16.6,4X,ES16.6,4X,ES16.6)')I,TEMP_INT2(I),TEMP_ETA2(I),&
			TEMP_INT3(I),DECAY_KIN_E(I)
	  E_DEP(I)=(TEMP_INT2(I)-TEMP_ETA2(I))+TEMP_INT3(I)+DECAY_KIN_E(I)
	END DO
	CLOSE(UNIT=7)
!
	TEMP_FLUX=0.0_LDP
	DO I=1,NF
	  DO J=1,R_MU(1)%MU_PTS-1
	    T1=R_MU(1)%MU_VECTOR(J)
	    T2=R_MU(1)%MU_VECTOR(J+1)
	    TEMP_FLUX(I)=TEMP_FLUX(I)+0.5_LDP*(T1-T2)*(T1*INTENSITY(J,1,I)+T2*INTENSITY(J,1,I))
	  END DO
	END DO
	FLUX=0.0_LDP
	DO I=1,NF-1
	  T1=NU_GRID(I)-NU_GRID(I+1)
	  FLUX=FLUX+0.5_LDP*T1*(TEMP_FLUX(I)+TEMP_FLUX(I+1))
	END DO
	FLUX=FLUX/2.0_LDP
!
	T1=0.0_LDP
	DO I=1,ND-1
	  T2=R(I)-R(I+1)
	  T3=R(I)*R(I)
	  T4=R(I+1)*R(I+1)
	  T1=T1+0.5_LDP*T2*(E_DEP(I)*T3+E_DEP(I+1)*T4)
	END DO
	T1=T1*1.0E+30_LDP*FOURPI
!
	OPEN(UNIT=7,FILE='./data/GAMRAY_E_DEP',STATUS='UNKNOWN',ACTION='WRITE')
	WRITE(7,'(A,1X,I3)')'! ND:',ND
	WRITE(7,'(A)')'! Energy deposited from scattering and photoelectric abs. ergs/s/cm^3'
	WRITE(7,'(A,1X,ES16.6)')'! Total Energy Absorbed:',T1
	WRITE(7,'(A,1X,ES16.6)')'! Total Eddington Flux at OB:',FLUX
	WRITE(7,'(A1,1X,A17,3X,4(A18,4X))')'!','Radius (10^10 cm)','V (km/s)','Energy dep.',&
		'local emission','decay K.E.'
	DO I=1,ND
	  WRITE(7,'(5(ES18.8,4X))')R(I),V(I),E_DEP(I),GAMRAY_EMISS(I),DECAY_KIN_E(I)
	END DO
	CLOSE(UNIT=7)
!
	RETURN
	END SUBROUTINE
