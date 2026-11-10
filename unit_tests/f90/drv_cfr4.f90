! Driver for CMF_FORMAL_REL_V4 (plane/cmf_formal_rel_v4.f).
!
! Alterations under test (bug fixes, not flag-gated):
!   1. Coastal-ray SQRT guards: MU at the outer/inner boundary uses
!      SQRT(MAX(0,(R-P)*(R+P))); a P numerically above R(1) or R(ND)
!      made upstream evaluate SQRT of a negative number - NaN into
!      MU_VALUE_CHK (outer) and FREQ_CONV_FAC (inner, HOLLOW only).
!   2. N_STORE floor under INNER_BND_METH='HOLLOW': the hollow-core
!      intensity store sizes itself as NINT(2*V(ND)); a coasting inner
!      velocity (V(ND)<0.5 km/s) truncates to 0, and MOD(CUR_LOC+1,0)
!      then faulted inside the frequency wrap.  Prod floors N_STORE at
!      1; upstream crashes.
!
! Cases:
!   ctl   - HOLLOW boundary, normal expanding V.  Identical path;
!           prod must match ref bit-for-bit.
!   coast - V(ND)=0.1 so NINT(2*V(ND))=0.  Prod completes through the
!           floored store; upstream hits MOD(x,0).
!   nudge - outermost tangent ray's P pushed 1e-7 above R(1).  Prod's
!           characteristics clamp + SQRT guard keep the run finite;
!           upstream feeds NaN Z/mu through the ray.
!
      PROGRAM DRV_CFR4
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: ND=12, NC=6, NP=ND+NC
      REAL(KIND=LDP) :: ETA(ND),CHI(ND),ESEC(ND),V(ND),SIGMA(ND),R(ND)
      REAL(KIND=LDP) :: P(NP)
      REAL(KIND=LDP) :: JNU(ND),FEDD(ND),IPLUS(NP)
      REAL(KIND=LDP) :: VDOP_VEC(ND)
      REAL(KIND=LDP) :: RP,DELMU,DBB,IC,B_PLANCK,FREQ,dLOG_NU
      REAL(KIND=LDP) :: REXT_FAC,VDOP_FRAC
      REAL(KIND=LDP) :: RET_HNU_AT_IB,RET_HNU_AT_OB
      CHARACTER(LEN=32) :: CNAME
      INTEGER :: I
!
      CALL GETARG(1,CNAME)
!
      DO I=1,ND
        R(I)=20.0_LDP-18.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        V(I)=2000.0_LDP-1800.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        SIGMA(I)=-0.5_LDP
        ETA(I)=1.0E-6_LDP
        CHI(I)=1.0E-8_LDP
        ESEC(I)=1.0E-10_LDP
        VDOP_VEC(I)=30.0_LDP
        JNU(I)=1.0E-8_LDP
        FEDD(I)=1.0_LDP/3.0_LDP
      END DO
      IPLUS(:)=0.0_LDP
!
      RP=R(ND)
      P(1)=0.0_LDP
      DELMU=1.0_LDP/NC
      DO I=2,NC
        P(I)=RP*SQRT(1.0_LDP-(DELMU*(NC-I+1))**2)
      END DO
      DO I=NC+1,NP
        P(I)=R(ND-I+NC+1)
      END DO
!
      IC=1.0_LDP
      FREQ=1.0_LDP
      dLOG_NU=1.0E-6_LDP
      REXT_FAC=1.0_LDP
      VDOP_FRAC=0.5_LDP
      DBB=0.0_LDP
      B_PLANCK=1.0E-4_LDP
!
      IF(CNAME .EQ. 'coast')THEN
        V(ND)=0.1_LDP
      ELSE IF(CNAME .EQ. 'nudge')THEN
        P(NP)=R(1)*1.0000001_LDP
      END IF
!
      CALL CMF_FORMAL_REL_V4(ETA,CHI,ESEC,V,SIGMA,R,P, &
         JNU,FEDD,RET_HNU_AT_IB,RET_HNU_AT_OB,IPLUS, &
         FREQ,dLOG_NU,B_PLANCK,DBB, &
         'HOLLOW',.FALSE., &
         VDOP_VEC,VDOP_FRAC,REXT_FAC, &
         'ZERO',.TRUE.,.TRUE.,NC,NP,ND)
!
      WRITE(6,'(A,*(ES14.5))')'JNU=',JNU
      WRITE(6,'(A,*(ES14.5))')'FEDD=',FEDD
      WRITE(6,'(A,2ES14.5)')'HNU=',RET_HNU_AT_IB,RET_HNU_AT_OB
      WRITE(6,*)'DONE_CFR4'
!
      END
