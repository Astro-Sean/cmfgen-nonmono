C
C Subroutine to compute the occuption probabilities for individual levels
C ignoring the detailed atomic structure. It is assumed that an approximate
C treatement will give better treatment than ignoring level dissolution
C completely.
C
C At present the code set the occupation probability to 1.0D0 exactly
C when neff =< 2Z.
C
C
C NB: The vectors X_LEV_DIS, A_LEV_DIS and B_LEV_DIS must have been previously
C        computed.
C
	SUBROUTINE OCCUPATION_PROB(W_C2,EDGE_C2,ZC2,NC2,ND)
	USE SET_KIND_MODULE
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
C
C Altered 12-Oct-2023 - Occupation prob set to 1 if >  0.99999999D0.
C Altered 08-Sep-2015 - Limit NEFF to be less than 31. All levels we include
C                         have n< 30. If neff> 30, mut have different core.
C                         Temporary fix -- we need to provide core information.
C Altered 03-Jul-2001 - Only computes NEFF for levels below the ionization
C                          limit. We need to include a check on the core.
C Altered 15-Dec-1997 - MOD_LEV_DIS_BLK replaces include file. Level
C                         dissolution can be switched off completely.
C Altered 27-Aug-1995 - K was declared integer instead of REAL
C Created 23-May-1995
C
	INTEGER NC2,ND
	REAL(KIND=LDP) W_C2(NC2,ND)		!Occupation probability
	REAL(KIND=LDP) EDGE_C2(NC2)		!Ionization frequency (10^15 Hz)
	REAL(KIND=LDP) ZC2			!Charge in ion
C
C Local variables.
C
	INTEGER I,LEV
	REAL(KIND=LDP) NEFF,F,Y,BETA,REAL_K,T1
C
	IF(MOD_DO_LEV_DIS)THEN
	  DO LEV=1,NC2
	    T1=3.289395_LDP*ZC2*ZC2/EDGE_C2(LEV)
	    NEFF=0.0_LDP
	    IF(T1 .GT. 0.0_LDP .AND. T1 .LT. 961.0_LDP)NEFF=SQRT(T1)
	    IF(NEFF .LE. 2.01_LDP*ZC2)THEN			!0.01 allows for atomic mass
	      DO I=1,ND
	        W_C2(LEV,I)=1.0_LDP
	      END DO
	    ELSE
	      REAL_K=16.0_LDP*NEFF/(1+NEFF)/(1+NEFF)/3 .0_LDP
	      IF(NEFF .LE. 3.0_LDP)REAL_K=1.0_LDP
	      REAL_K=( REAL_K/ZC2*(ZC2/NEFF)**4 )**1.5_LDP
	      DO I=1,ND
	        Y=1.091_LDP*(X_LEV_DIS(I)+4.0_LDP*(ZC2-1.0_LDP)*A_LEV_DIS(I))
	        BETA=REAL_K*B_LEV_DIS(I)
	        F=Y*BETA*BETA/(7.782_LDP+X_LEV_DIS(I)*BETA)
	        W_C2(LEV,I)=F/(1.0_LDP+F)
	        IF(W_C2(LEV,I) .GT. 0.99999999_LDP)W_C2(LEV,I)=1.0_LDP
	      END DO
	    END IF
	  END DO
	ELSE
	  W_C2(:,:)=1.0_LDP
	END IF
C
	RETURN
	END
