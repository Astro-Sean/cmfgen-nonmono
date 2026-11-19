! drv_knop_endpoint.f - independent verification of the KNOP
! P-direction outer-endpoint (j=1) source-slope reconstruction.
!
! The production scheme (KN_ORD_DS, mirrored in KN_SH_EVAL_H) is:
!   SA = (S(2)-S(1))/dt(1)
!   SB = (S(3)-S(2))/dt(2)
!   DD = SA + (SA-SB)*dt(1)/(dt(1)+dt(2))
!   DS = (sign(SA)+sign(DD)) * min(|SA|, 0.5*|DD|)
! a secant extrapolation to the node with a two-factor monotonicity
! limiter. The PROD column below re-evaluates exactly that formula.
!
! The REF/REF3 columns are the independent reference: the derivative
! at tau(1) of the Lagrange interpolating polynomial through the
! first 3 (resp. 4) grid points, evaluated by direct term-by-term
! differentiation of the cardinal basis -- a different algebraic
! path with no limiter and no secant extrapolation.
!
! Manufactured cases A-E per the phase-5.5 audit task. The tau
! grid is geometric: tau(j) = T0 + S0*(R^(j-1)-1), so dt(1)=S0*(R-1)
! halves per refinement level while the stretch ratio stays fixed.
! DUPD reports |DS_err|*|A4REP| with A4REP=1: a unit-
! representative formal-update discrepancy, NOT the actual
! PK_A4(1,L)-weighted production transfer error (the real
! coefficient varies per cell and is not measured here).
!
      PROGRAM DRV_KNOP_ENDPOINT
      IMPLICIT NONE
      INTEGER, PARAMETER :: LDP=SELECTED_REAL_KIND(15,307)
      INTEGER, PARAMETER :: NI=64,NLV=4,NC=5
      REAL(KIND=LDP), PARAMETER :: T0=1.0_LDP,R=1.05_LDP
      REAL(KIND=LDP), PARAMETER :: A4REP=1.0_LDP
      REAL(KIND=LDP), ALLOCATABLE :: TAU(:),S(:)
      REAL(KIND=LDP) S0,DT1,DT2,TDT,SA,SB,DD,DSP,D2,D3
      REAL(KIND=LDP) X1,X2,X3,X4,D12,D13,D23
      REAL(KIND=LDP) ANA,EP,ER,E3,EPP(NLV),ERR(NLV),HH(NLV)
      REAL(KIND=LDP) SGN_P,SGN_D
      INTEGER IC,LV,J
      CHARACTER*1 CN(5)
      DATA CN/'A','B','C','D','E'/
      WRITE(6,'(A)')'CASE LV      H             ANA          '//
     1 ' PROD          ERR_P     REL_P       REF           '//
     2 'ERR_R     REL_R       ERR_3       DUPD'
      DO IC=1,NC
        DO LV=1,NLV
          S0=0.05_LDP/2.0_LDP**(LV-1)
          ALLOCATE(TAU(NI),S(NI))
          DO J=1,NI
            TAU(J)=T0+S0*(R**(J-1)-1.0_LDP)
            S(J)=SFUN(IC,TAU(J),T0)
          END DO
          DT1=TAU(2)-TAU(1)
          DT2=TAU(3)-TAU(2)
          TDT=DT1+DT2
          ANA=DFUN(IC,TAU(1),T0)
!
! Production scheme under test.
!
          SA=(S(2)-S(1))/DT1
          SB=(S(3)-S(2))/DT2
          IF(ABS(TDT) .LT. 1.0E-30_LDP)THEN
            DD=0.0_LDP
          ELSE
            DD=SA+(SA-SB)*DT1/TDT
          END IF
          SGN_P=SIGN(1.0_LDP,SA)
          SGN_D=SIGN(1.0_LDP,DD)
          DSP=(SGN_P+SGN_D)*MIN(ABS(SA),0.5_LDP*ABS(DD))
!
! Independent reference: d/dtau at tau(1) of the Lagrange basis
! through points 1..3 (D2) and 1..4 (D3). Basis derivatives at
! x=x1: l1' = ((x1-x2)+(x1-x3))/((x1-x2)(x1-x3)),
! l2' = (x1-x3)/((x2-x1)(x2-x3)),
! l3' = (x1-x2)/((x3-x1)(x3-x2)).
!
          D12=TAU(1)-TAU(2)
          D13=TAU(1)-TAU(3)
          D23=TAU(2)-TAU(3)
          D2=S(1)*(D12+D13)/(D12*D13)
     1       +S(2)*D13/((TAU(2)-TAU(1))*D23)
     2       +S(3)*D12/((TAU(3)-TAU(1))*(TAU(3)-TAU(2)))
          X1=TAU(1); X2=TAU(2); X3=TAU(3); X4=TAU(4)
          D3=S(1)*(((X1-X2)*(X1-X3))+
     1       ((X1-X2)+(X1-X3))*(X1-X4))/ 
     2       ((X1-X2)*(X1-X3)*(X1-X4))
     3       +S(2)*(X1-X3)*(X1-X4)/
     4       ((X2-X1)*(X2-X3)*(X2-X4))
     5       +S(3)*(X1-X2)*(X1-X4)/
     6       ((X3-X1)*(X3-X2)*(X3-X4))
     7       +S(4)*(X1-X2)*(X1-X3)/
     8       ((X4-X1)*(X4-X2)*(X4-X3))
          EP=ABS(DSP-ANA)
          ER=ABS(D2-ANA)
          E3=ABS(D3-ANA)
          EPP(LV)=EP
          ERR(LV)=ER
          HH(LV)=DT1
          WRITE(6,'(A3,I2,1X,ES12.5,1X,ES13.6,1X,ES13.6,1X,ES11.3,
     1     1X,ES11.3,1X,ES13.6,1X,ES11.3,1X,ES11.3,1X,ES11.3,
     2     1X,ES11.3)')
     3     ' '//CN(IC)//' ',LV,DT1,ANA,DSP,EP,
     4     EP/MAX(ABS(ANA),1.0E-30_LDP),D2,ER,
     5     ER/MAX(ABS(ANA),1.0E-30_LDP),E3,
     6     EP*A4REP
          DEALLOCATE(TAU,S)
        END DO
        DO LV=2,NLV
          IF(ERR(LV-1) .GT. 0.0_LDP .AND. EPP(LV-1) .GT.
     1       0.0_LDP)THEN
            WRITE(6,'(A,I2,A,2(F6.2,1X))')'  '//CN(IC)//
     1       ' order lv',LV,':',
     2       LOG(EPP(LV-1)/EPP(LV))/LOG(HH(LV-1)/HH(LV)),
     3       LOG(ERR(LV-1)/ERR(LV))/LOG(HH(LV-1)/HH(LV))
          END IF
        END DO
      END DO
!
      CONTAINS
!
      REAL(KIND=LDP) FUNCTION SFUN(IC,T,T0)
      INTEGER IC
      REAL(KIND=LDP) T,T0,U
      U=T-T0
      IF(IC .EQ. 1)THEN
        SFUN=2.0_LDP
      ELSE IF(IC .EQ. 2)THEN
        SFUN=0.5_LDP*T
      ELSE IF(IC .EQ. 3)THEN
        SFUN=0.3_LDP*T*T-0.1_LDP*T
      ELSE IF(IC .EQ. 4)THEN
        SFUN=EXP(0.7_LDP*T)
      ELSE
        SFUN=1.0_LDP-0.4_LDP*COS(4.0_LDP*3.141592653589793_LDP*U/0.5_LDP)
      END IF
      RETURN
      END FUNCTION SFUN
!
      REAL(KIND=LDP) FUNCTION DFUN(IC,T,T0)
      INTEGER IC
      REAL(KIND=LDP) T,T0,U
      U=T-T0
      IF(IC .EQ. 1)THEN
        DFUN=0.0_LDP
      ELSE IF(IC .EQ. 2)THEN
        DFUN=0.5_LDP
      ELSE IF(IC .EQ. 3)THEN
        DFUN=0.6_LDP*T-0.1_LDP
      ELSE IF(IC .EQ. 4)THEN
        DFUN=0.7_LDP*EXP(0.7_LDP*T)
      ELSE
        DFUN=0.4_LDP*(4.0_LDP*3.141592653589793_LDP/0.5_LDP)*
     1    SIN(4.0_LDP*3.141592653589793_LDP*U/0.5_LDP)
      END IF
      RETURN
      END FUNCTION DFUN
!
      END
