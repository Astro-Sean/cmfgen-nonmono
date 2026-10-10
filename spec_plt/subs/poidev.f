!
! Returns as a floating point number an INTEGER VALUE that is a random deviate
! drawn from a Posson distribution of mean XM using RAN2_DP as a soure of
! random deviates.
!
! Routine is from Numerical Recipes.
!
      FUNCTION POIDEV(XM,IDUM)
	USE SET_KIND_MODULE
      IMPLICIT NONE
!
! Altered: 16-Apr-2022 : Error - ALXM was not being correctly saved.
!                           Not apparent incontinuous calls to POIDEV if
!                           values distinct.
! Altered: 10-Apr-2015 : RAN2_DP accessed (there was a precision issue)
!                        SAVE variables fixed.
!
      REAL(KIND=LDP) XM
      REAL(KIND=LDP) POIDEV
      INTEGER IDUM
!
      REAL(KIND=LDP) T
      REAL(KIND=LDP) Y
      REAL(KIND=LDP) GAMMLN
      REAL(KIND=LDP) EM
!
      REAL(KIND=LDP) RAN2_DP
      REAL(KIND=LDP), SAVE :: G,SQ,ALXM
      REAL(KIND=LDP), SAVE :: OLDM=-1
      REAL(KIND=LDP), PARAMETER :: PI=3.141592654_LDP
!
      IF (XM .LT. 12.0_LDP)THEN
        IF (XM.NE.OLDM) THEN
          OLDM=XM
          G=EXP(-XM)
        ENDIF
        EM=-1
        T=1.0_LDP
2       EM=EM+1.0_LDP
        T=T*RAN2_DP(IDUM)
        IF (T .GT. G) GO TO 2
      ELSE
        IF (XM.NE.OLDM) THEN
          OLDM=XM
          SQ=SQRT(2*XM)
          ALXM=LOG(XM)
          G=XM*ALXM-GAMMLN(XM+1.0_LDP)
        ENDIF
1       Y=TAN(PI*RAN2_DP(IDUM))
        EM=SQ*Y+XM
        IF (EM .LT. 0.0_LDP) GO TO 1
        EM=INT(EM)
        T=0.9_LDP*(1.0_LDP+Y**2)*EXP(EM*ALXM-GAMMLN(EM+1.0_LDP)-G)
        IF (RAN2_DP(IDUM) .GT. T) GO TO 1
      ENDIF
      POIDEV=EM
!
      RETURN
      END
!
      FUNCTION GAMMLN(XX)
	USE SET_KIND_MODULE
      IMPLICIT NONE
      REAL(KIND=LDP) XX
      REAL(KIND=LDP) COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      REAL(KIND=LDP) GAMMLN
      INTEGER J
!
      DATA COF,STP/76.18009173_LDP,-86.50532033_LDP,24.01409822_LDP,
     *    -1.231739516_LDP,.120858003E-2_LDP,-.536382E-5_LDP,2.50662827465_LDP/
      DATA HALF,ONE,FPF/0.5_LDP,1.0_LDP,5.5_LDP/
 !
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
