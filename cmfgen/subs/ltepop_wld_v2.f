!
! General routine to:
!
!       Compute the LTE populations of the levels in the FULL atom given
!       ED (electron density), T (electron temperature) and the ion density
!       DI C2. Level dissolution is taken into account.
!
	SUBROUTINE LTEPOP_WLD_V2(C2LTE,LOG_C2LTE,W_C2,EDGEC2,GC2,
	1             ZC2,GION_C2,NC2,DIC2,ED,T,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered: 8-Nov-2023: Use fortran RANGE function, and changed MAX_LOG_LTE_POP to MAX_LN_LTE_POP.
! Altered: 5-Apr-2011: MAX_LOG_LTE_POP parameter introduced.
!                      Based on LTEPOP_SL_V1 (original coding early 2011)
!                      Call changed as LOG_C2LTE variable introduced.
!
! Altered 18-Feb-2010 : Take log of density to extend range of LTE populations.
!                         We use NEW_METHOD to allow easy change to previous
!                         version. Only necessary if something untoward crops up
! Altered 24-May-1996 : ND_MAX removed (was unused).
! Created 30-May-1995 : Based on LTEPOP.
!
	INTEGER ND
	REAL(KIND=LDP) ED(ND)			!Electron density
	REAL(KIND=LDP) T(ND)			!Temperature 10^4K
	REAL(KIND=LDP) DIC2(ND)			!Ion density (Full model atom)
!
	INTEGER NC2
	REAL(KIND=LDP) C2LTE(NC2,ND)
	REAL(KIND=LDP) LOG_C2LTE(NC2,ND)
	REAL(KIND=LDP) W_C2(NC2,ND)
	REAL(KIND=LDP) EDGEC2(NC2)
	REAL(KIND=LDP) GC2(NC2)
	REAL(KIND=LDP) GION_C2			!Statistical weight of ion groun state.
	REAL(KIND=LDP) ZC2			!Ion charge
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER I,K
	REAL(KIND=LDP) X,Y,RGU
	REAL(KIND=LDP) MAX_LN_LTE_POP
!
! RANGE retruns the maximum exponent range for the FLOATING POINT variable. We
!   multiply by Log 10 as we compare with natural logarithms. We subtract
!   10 as an extra precaution.
!
	X=10.0_LDP
	MAX_LN_LTE_POP=0.9_LDP*RANGE(X)*LOG(X)
!
! Compute the occupation probabilities.
!
	CALL OCCUPATION_PROB(W_C2,EDGEC2,ZC2,NC2,ND)
!
! Compute the LTE populations of the levels in the full atom, taking level
! dissolution into account.
!
	C2LTE=0.0_LDP
	DO K=1,ND
	 X=HDKT/T(K)
	 RGU=2.07078E-22_LDP*ED(K)*DIC2(K)*( T(K)**(-1.5_LDP) )/GION_C2
	 IF(RGU .LE. 0)THEN
	   WRITE(6,*)'Error in SETPOP_WLD_V2 -- zero RGU'
	   WRITE(6,*)NC2,K,ED(K),DIC2(K),T(K),GION_C2
	   STOP
	 END IF
	 RGU=LOG(RGU)
	 DO I=1,NC2
	   LOG_C2LTE(I,K)=LOG(W_C2(I,K)*GC2(I))+EDGEC2(I)*X+RGU
	   IF(LOG_C2LTE(I,K) .LE. MAX_LN_LTE_POP)C2LTE(I,K)=EXP(LOG_C2LTE(I,K))
	 END DO
	END DO
!
	RETURN
	END
