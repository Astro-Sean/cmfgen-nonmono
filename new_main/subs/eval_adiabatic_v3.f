!
! Routine to adjust the Radiative equilibrium equation for adiabatic cooling.
! This routine:
!
!  (1) Increments STEQ if adiabatic cooling is to be allowed for
!  (2) Increments the variation matrix [BA] if it is being computed, and
!       if adiabatic cooling is being included.
!  (3) Evaluates the adiabatic cooling terms, splitting them into 2 terms ---
!       the dTdR term, and the velocity term. This is done for diagnostic
!       purposes.
!
 	SUBROUTINE EVAL_ADIABATIC_V3(AD_CR_V,AD_CR_DT,
	1                     POPS,AVE_ENERGY,HDKT,COMPUTE_BA,INCL_ADIABATIC,
	1                     DIAG_INDX,NUM_BNDS,NT,ND)
!
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
 	USE STEQ_DATA_MOD
	USE CONTROL_VARIABLE_MOD, ONLY : USE_ELEC_HEAT_BAL, COMP_STEQ_T_EHB
 	IMPLICIT NONE
!
! Altered  13-Jun-2-25 : Changed dimension of SUM_EN to NT from ND.
! Altered  29-Dec-2023 : Fixed evaluation of STEQ_T_EHB and its variation.
! Altered     Dec-2023 : Added evaluation of STEQ_T_EHB. 
! Altered  21-Jun-2004 : Changed to version V3.
!                          Changed to use simple linear differencing.
!                          Changed to incorporate advection terms from rate equations
!                            in electron cooling equation.
!                          Changed to allow for changes in internal (excitation) energy.
! Altered  16-Dec-2002 : DTDR computed using forward differencing.
! Altered  12-Apr-2001 : Changed to use STEQ_DATA_MOD
!                        Changed to V2 (call changed).
! Altered  15-Mar-2001 : Minor bug fix in BA matrix at d=1.
! Created 256-Jul-1994
!
	INTEGER NT
	INTEGER ND
	INTEGER DIAG_INDX
	INTEGER NUM_BNDS
!
! Output:
!
	REAL(KIND=LDP) AD_CR_V(ND)
	REAL(KIND=LDP) AD_CR_DT(ND)
!
! Input:
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) AVE_ENERGY(NT)
	REAL(KIND=LDP) HDKT
!
! Local vectors.
!
	REAL(KIND=LDP) A(ND)
	REAL(KIND=LDP) B(ND)
	REAL(KIND=LDP) C(ND)
	REAL(KIND=LDP) D(ND)
	REAL(KIND=LDP) GAMMA(ND)
	REAL(KIND=LDP) INT_EN(ND)
	REAL(KIND=LDP) COL_EN(ND)
	REAL(KIND=LDP) WORK(ND)
!
	REAL(KIND=LDP) ION_EN(NT)
	REAL(KIND=LDP) TOT_ENERGY(NT)
!
! Local variables.
!
	LOGICAL COMPUTE_BA,INCL_ADIABATIC
!
	INTEGER ERROR_LU
	REAL(KIND=LDP) BOLTZMANN_CONSTANT,FUN_PI
	EXTERNAL BOLTZMANN_CONSTANT,FUN_PI,ERROR_LU
!
	REAL(KIND=LDP) SCALE,T1,T2,PI
	INTEGER I,J,L
	INTEGER LUER
	INTEGER ISPEC
	INTEGER ID
	LOGICAL WRITE_CHK
!
! A full linearization is now obsolete, but check to make sure.
!
	LUER=ERROR_LU()
	IF(NUM_BNDS .EQ. ND)THEN
	  WRITE(LUER,*)'Error --- EVAL_ADIABATIC can''t handle a full linearization'
	  STOP
	END IF
!
! Compute the total excitation energy of each level.
!
	TOT_ENERGY(1:NT)=0.0_LDP
	ION_EN(1:NT)=0.0_LDP
	DO ISPEC=1,NUM_SPECIES
	  T1=0.0_LDP
	  T2=0.0_LDP
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    T2=T2+AVE_ENERGY(ATM(ID)%EQXzV)
	    DO I=1,ATM(ID)%NXzV
	      J=ATM(ID)%EQXzV+I-1
	      TOT_ENERGY(J)=(AVE_ENERGY(ATM(ID)%EQXzV)-AVE_ENERGY(J))+T1
	      ION_EN(J)=T2
	    END DO
	    J=ATM(ID)%EQXzV
	    T1=T1+AVE_ENERGY(J)			!Adding on ionization energy
	  END DO
	  ID=SPECIES_END_ID(ISPEC)-1
	  IF(ID .GT. 0)THEN
	    J=ATM(ID)%EQXzV
	    TOT_ENERGY(J+ATM(ID)%NXzV)=T1
	    ION_EN(J+ATM(ID)%NXzV)=T2
	  END IF
	END DO
!
! Compute the mean energy per atom. At first it is units of 10^15Hz.
!
	INT_EN(:)=0.0_LDP
	COL_EN(:)=0.0_LDP
	DO I=1,ND
	  DO J=1,NT-2
	     INT_EN(I)=INT_EN(I)+POPS(J,I)*TOT_ENERGY(J)
	     COL_EN(I)=COL_EN(I)+POPS(J,I)*ION_EN(J)
	  END DO
	END DO
	INT_EN=HDKT*INT_EN/POP_ATOM
	COL_EN=HDKT*COL_EN/POP_ATOM
!
	IF(USE_ELEC_HEAT_BAL .OR. COMP_STEQ_T_EHB)THEN
	  SCALE=0.1_LDP*BOLTZMANN_CONSTANT()
	  CALL UPDATE_BA_ST(BA_T_EHB,STEQ_T_EHB,COL_EN,ION_EN)
	END IF
!
! For historical resions STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the adiabatic cooling rate by that
! term. Not also that the R units are unimportant, since R.Chi is dimensionless.
! The 10^9 arises since T is in units of 10^4K, and V in units of 10^5  km/s.
!
	PI=FUN_PI()
	SCALE=1.0E+09_LDP*BOLTZMANN_CONSTANT()/4.0_LDP/PI
	CALL UPDATE_BA_ST(BA_T,STEQ_T,INT_EN,TOT_ENERGY)

!
! Now compute the adiabatic cooling rate (in ergs/cm^3/sec) for diagnostic
! purposes. The rate is output to the COOLGEN file.
!
! The factor of 4.0D-10*PI arises from the fact that A, B etc were computed
! for the radiative equilibrium equation. That equation has units a factor of 10^10/4Pi
! larger (10^10 because of opacity definition, 4Pi as we don't scale J).
!
! We split the adiabatic terms into 2 parts: The velocity term, and the
! dTdR (internal energy) term. This split was useful for diagnostic purposes,
! but has now been kept for simplicity so that GENCOOL does not need to be changed.
!
	DO I=1,ND
	  WRITE(277,'(10ES14.4)')A(I),B(I),C(I),D(I),T(I),GAMMA(I),COL_EN(I)
	  FLUSH(UNIT=277)
	END DO
!
	T1=4.0E-10_LDP*PI
	DO I=1,ND-1
	  AD_CR_V(I)= B(I)*T(I)
	  AD_CR_DT(I)=A(I)*(T(I)-T(I+1))+
	1               C(I)*T(I)*(GAMMA(I)-GAMMA(I+1)) +
	1               D(I)*(COL_EN(I)-COL_EN(I+1))
	END DO
	AD_CR_V(ND)=B(ND)*T(ND)
	AD_CR_DT(ND)=A(ND)*(T(ND-1)-T(ND))+
	1               C(ND)*T(ND)*(GAMMA(ND-1)-GAMMA(ND)) +
	1               D(ND)*(COL_EN(ND-1)-COL_EN(ND))
!
	AD_CR_V=AD_CR_V*T1
	AD_CR_DT=AD_CR_DT*T1
!
	WRITE_CHK=.TRUE.
	IF(WRITE_CHK)THEN
	  OPEN(UNIT=7,FILE='ADIABAT_CHK',STATUS='UNKNOWN')
	    WRITE(7,'(A)')' '
	    WRITE(7,'(A)')'  Scaling is for STEQ_VALS(NT,:). This is cgs units scaled'
	    WRITE(7,'(A)')'  by a factor of 10^10 on 4PI. The 10^9 comes from V . T'
	    WRITE(7,'(A)')' '
	    WRITE(7,'(A,ES12.4)')' SCALE=1.0D+09*BOLTZMANN_CONSTANT()/4.0D0/PI=',SCALE
	    WRITE(7,'(A)')' T1=R(I)-R(I+1)'
	    WRITE(7,'(A)')' GAMMA=ED/POP_ATOM'
	    WRITE(7,'(A)')' A=1.5D0*SCALE*(POP_ATOM+ED)*V/T1 * (T(I)-T(I+1)'
	    WRITE(7,'(A)')' B=SCALE*(POP_ATOM+ED)*V*(3.0D0+SIGMA)/R * T(I)'
	    WRITE(7,'(A)')' C=1.5D0*SCALE*POP_ATOM*V/T1 * T(I)*(GAMMA(I)-GAMMA(I+1))'
	    WRITE(7,'(A)')' D=SCALE*POP_ATOM*V/T1* (INT_EN(I)-INT_EN(I+1))'
	    WRITE(7,'(A)')' '
	    WRITE(7,'(8X,A,6X,8(7X,A))')'R','     V',' SIGMA','     T',
	1                 'NU_ION','     A','     B','     C','     D'
	    DO I=1,ND-1
	      WRITE(7,'(ES15.7,8(1X,ES12.4))')
	1              R(I),V(I),SIGMA(I),T(I),INT_EN(I)/HDKT,
	1              A(I)*(T(I)-T(I+1)),B(I)*T(I),
	1              C(I)*T(I)*(GAMMA(I)-GAMMA(I+1)),D(I)*(INT_EN(I)-INT_EN(I+1))
	    END DO
	    WRITE(7,'(ES15.7,8(1X,ES12.4))')
	1              R(ND),V(ND),SIGMA(ND),T(ND),INT_EN(ND)/HDKT,
	1              A(ND)*(T(ND-1)-T(ND)),B(ND)*T(ND),
	1              C(ND)*T(ND)*(GAMMA(ND-1)-GAMMA(ND)),D(ND)*(INT_EN(ND-1)-INT_EN(ND))
	  CLOSE(UNIT=7)
	END IF
!
	RETURN
!
	CONTAINS

	SUBROUTINE UPDATE_BA_ST(NEW_BA_T,NEW_STEQ_T,MEAN_EN,SUM_EN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	REAL(KIND=LDP) NEW_BA_T(NT,NUM_BNDS,ND)
	REAL(KIND=LDP) NEW_STEQ_T(ND)
	REAL(KIND=LDP) MEAN_EN(ND)
	REAL(KIND=LDP) SUM_EN(NT)
!
! We now compute constants for each of the 4 terms. These make
! it simpler and cleaner for the evaluation of the linearization.
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the adiabatic cooling rate by that
! term. Not also that the R units are unimportant, since R.Chi is dimensionless.
! The 10^9 arises since T is in units of 10^4K, and V in units of 10^5 km/s.
!
	DO I=1,ND
	  IF(I .EQ. ND)THEN
	    T1=R(ND-1)-R(ND)
	  ELSE
	    T1=R(I)-R(I+1)
	  END IF
	  A(I)=1.5_LDP*SCALE*(POP_ATOM(I)+ED(I))*V(I)/T1
	  B(I)=SCALE*(POP_ATOM(I)+ED(I))*V(I)*(3.0_LDP+SIGMA(I))/R(I)
	  C(I)=1.5_LDP*SCALE*POP_ATOM(I)*V(I)/T1
	  D(I)=SCALE*POP_ATOM(I)*V(I)/T1
	  GAMMA(I)=ED(I)/POP_ATOM(I)
	END DO
!
	IF(INCL_ADIABATIC)THEN
	  DO I=1,ND-1
 	    WORK(I)=A(I)*(T(I)-T(I+1)) + B(I)*T(I) +
	1              C(I)*T(I)*(GAMMA(I)-GAMMA(I+1)) +
	1              D(I)*(INT_EN(I)-INT_EN(I+1))
	  END DO
 	  WORK(ND)=A(ND)*(T(ND-1)-T(ND)) + B(ND)*T(ND) +
	1              C(ND)*T(ND)*(GAMMA(ND-1)-GAMMA(ND)) +
	1              D(ND)*(MEAN_EN(ND-1)-MEAN_EN(ND))
!
	  DO I=1,ND
	    NEW_STEQ_T(I)=NEW_STEQ_T(I)-WORK(I)
	  END DO
	END IF
!
	IF(INCL_ADIABATIC .AND. COMPUTE_BA)THEN
	  DO I=1,ND-1
!
! Diagonal terms.
!
	    L=DIAG_INDX
	    NEW_BA_T(NT,L,I)=NEW_BA_T(NT,L,I)-A(I)-B(I)-C(I)*(GAMMA(I)-GAMMA(I+1))
	    NEW_BA_T(NT-1,L,I)=NEW_BA_T(NT-1,L,I)-(A(I)+B(I))/(POP_ATOM(I)+ED(I))-
	1                                  C(I)*T(I)/POP_ATOM(I)
	    DO J=1,NT-2
	      NEW_BA_T(J,L,I)=NEW_BA_T(J,L,I)-HDKT*D(I)*SUM_EN(J)/POP_ATOM(I)
	    END DO
!
! Upper diagonal terms.
!
	    IF(NUM_BNDS .GE. 3)THEN
	      L=DIAG_INDX+1
	      NEW_BA_T(NT,L,I)=NEW_BA_T(NT,L,I)+A(I)
	      NEW_BA_T(NT-1,L,I)=NEW_BA_T(NT-1,L,I)+C(I)*T(I)/POP_ATOM(I+1)
	      DO J=1,NT-1
	        NEW_BA_T(J,L,I)=NEW_BA_T(J,L,I)+HDKT*D(I)*SUM_EN(J)/POP_ATOM(I+1)
	      END DO
	    END IF
!
	  END DO	!Loop of depth.
!
! Need to do special case of I=ND
!
	  L=DIAG_INDX
	  NEW_BA_T(NT,L,ND)=NEW_BA_T(NT,L,ND)+A(ND)-B(ND)-C(ND)*(GAMMA(ND-1)-GAMMA(ND))
	  NEW_BA_T(NT-1,L,ND)=NEW_BA_T(NT-1,L,ND)-(A(ND)+B(ND))/(POP_ATOM(ND)+ED(ND)) +
	1                                  C(ND)*T(ND)/POP_ATOM(ND)
	  DO J=1,NT-2
	    NEW_BA_T(J,L,ND)=NEW_BA_T(J,L,ND)+HDKT*D(ND)*SUM_EN(J)/POP_ATOM(ND)
	  END DO
!
	  IF(NUM_BNDS .GE. 3)THEN
	    L=DIAG_INDX-1
	    NEW_BA_T(NT,L,ND)=NEW_BA_T(NT,L,ND)-A(ND)
	    NEW_BA_T(NT-1,L,ND)=NEW_BA_T(NT-1,L,ND)-C(ND)*T(ND)/POP_ATOM(ND-1)
	    DO J=1,NT-2
	      NEW_BA_T(J,L,ND)=NEW_BA_T(J,L,ND)-HDKT*D(ND)*SUM_EN(J)/POP_ATOM(ND-1)
	    END DO
	  END IF
!
	END IF            !End COMPUTE_BA
!
	RETURN
	END SUBROUTINE UPDATE_BA_ST
 	END SUBROUTINE EVAL_ADIABATIC_V3
