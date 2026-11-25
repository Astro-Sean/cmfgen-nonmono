! Standalone driver for MOM_AUG_RESID_V1 in MOD_MOM_J_NONMONO_V1.
! Builds a manufactured captured column set (NU_* arrays, J_NEXT_FRQ,
! H_NEXT_FRQ, MOM_NM_J/H, stabiliser flags) on a small grid, calls the
! augmented-residual evaluator, and dumps every input array so the
! python battery can recompute the residuals independently.
!
! For the "exact" family of cases the driver also manufactures the
! committed field: NU_SRC is set so every assembled (implemented)
! moment row is satisfied to round-off, and H_NEXT_FRQ is filled
! column-by-column in upwind order so every H row is exact.
!
! Usage: drv_mom_augres <case> <tag>
!
      PROGRAM DRV_MOM_AUGRES
      USE SET_KIND_MODULE
      USE MOD_MOM_J_NONMONO_V1
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE
!
      INTEGER ND,NCF
      PARAMETER (ND=5,NCF=10)
!
      CHARACTER(LEN=32) CASE,TAG
      INTEGER I,K,U,U1,U2,D,IOS,ITER
      INTEGER OUTBC
      REAL(KIND=LDP) T1B,THI,THB,TE1,TE2,TV,XUP,SCL,T5,T4,XMP
      REAL(KIND=LDP) F0,CKSUM
      LOGICAL LCLIP
!
      CALL GETARG(1,CASE)
      CALL GETARG(2,TAG)
      IF(LEN_TRIM(CASE) .EQ. 0 .OR. LEN_TRIM(TAG) .EQ. 0)
     1   STOP 'DRV_MOM_AUGRES: need <case> <tag>'
!
! Allocate every array the evaluator reads.
!
      ALLOCATE(NU_TA(ND,NCF),NU_TB(ND,NCF),NU_TC(ND,NCF),
     1   NU_X0(ND,NCF),NU_CUP(ND,NCF),NU_CUPU(ND,NCF),
     2   NU_KUP(ND,NCF),
     3   NU_HU(ND-1,NCF),NU_HL(ND-1,NCF),NU_HS(ND-1,NCF),
     4   NU_VB(ND,NCF),NU_VC(ND,NCF),
     5   NU_EPS(ND-1,NCF),NU_EPSU(ND-1,NCF),
     6   NU_EPSP(ND-1,NCF),NU_EPSPU(ND-1,NCF),
     7   NU_LAGSCL(ND,NCF),NU_SRC(ND,NCF),NU_T1B(NCF),
     8   NU_KUH(ND-1,NCF),NU_OUTBC(NCF),NU_HCHK(NCF),
     9   J_NEXT_FRQ(ND,NCF),H_NEXT_FRQ(ND,NCF),
     1   MOM_NM_J(ND),MOM_NM_H(ND))
!
      MOM_NM_J(1:ND)=1
      MOM_NM_H(1:ND)=1
      MOM_NM_J(4)=2; MOM_NM_J(5)=2
      MOM_NM_H(4)=2
!
! Stabiliser flags: all on by default; individual cases toggle.
!
      MOM_XM_POS=.TRUE.
      MOM_EPS_LIM=.TRUE.
      MOM_H_IMPL=.FALSE.
      IF(CASE .EQ. 'him')MOM_H_IMPL=.TRUE.
!
! Base coefficient set. TA/TB/TC hold the implemented depth-upwind
! EPS stencil (built below from the post-cap NU_EPS values), plus
! the MOM_H_IMPL diagonal mass when that case runs it.
!
      DO K=1,NCF
        DO I=1,ND
          NU_TA(I,K)=-0.40_LDP-0.05_LDP*I-0.01_LDP*K
          NU_TB(I,K)=3.0_LDP+0.10_LDP*I+0.02_LDP*K
          NU_TC(I,K)=-0.30_LDP-0.04_LDP*I-0.01_LDP*K
          NU_CUP(I,K)=0.0_LDP
          NU_CUPU(I,K)=0.0_LDP
          NU_KUP(I,K)=0
          NU_VB(I,K)=0.0_LDP
          NU_VC(I,K)=0.0_LDP
          NU_LAGSCL(I,K)=1.0_LDP
          NU_SRC(I,K)=0.0_LDP
          NU_X0(I,K)=0.0_LDP
          J_NEXT_FRQ(I,K)=1.0_LDP+0.30_LDP*I+0.05_LDP*K+
     1      0.10_LDP*SIN(0.7_LDP*I*K)
          H_NEXT_FRQ(I,K)=0.0_LDP
        END DO
        DO I=1,ND-1
          NU_HU(I,K)=0.50_LDP+0.03_LDP*I+0.01_LDP*K
          NU_HL(I,K)=0.40_LDP+0.02_LDP*I+0.01_LDP*K
          NU_HS(I,K)=0.20_LDP+0.01_LDP*I
          NU_EPS(I,K)=0.0_LDP
          NU_EPSU(I,K)=0.0_LDP
          NU_EPSP(I,K)=0.0_LDP
          NU_EPSPU(I,K)=0.0_LDP
          NU_KUH(I,K)=0
        END DO
        NU_OUTBC(K)=1
        NU_T1B(K)=0.0_LDP
        NU_HCHK(K)=0
      END DO
      DO I=2,ND-1
        DO K=1,NCF
          NU_VB(I,K)=-0.10_LDP-0.01_LDP*I
          NU_VC(I,K)=0.12_LDP+0.01_LDP*I
        END DO
      END DO
!
! Frequency directions. Default: cols 1-4 expanding (upwind k-1),
! cols 5-8 compressing (upwind k+NM), cols 9-10 expanding. The
! compressing stretch stays far enough from the red edge that every
! needed index is in grid.
!
      DO K=1,NCF
        IF(K .GE. 5 .AND. K .LE. 8)THEN
          DO I=2,ND-1
            U=K+MOM_NM_J(I)
            IF(U .LE. NCF)THEN
              NU_KUP(I,K)=U
              NU_CUPU(I,K)=0.45_LDP+0.02_LDP*I
              NU_CUP(I,K)=NU_CUPU(I,K)
            END IF
          END DO
          U=K+MOM_NM_J(1)
          IF(U .LE. NCF)THEN
            NU_KUP(1,K)=U
            NU_CUPU(1,K)=-0.50_LDP
            NU_CUP(1,K)=-0.50_LDP
          END IF
          DO I=1,ND-1
            U=K+MOM_NM_H(I)
            IF(U .LE. NCF)THEN
              NU_KUH(I,K)=U
            ELSE
              NU_HS(I,K)=0.0_LDP
            END IF
          END DO
          NU_OUTBC(K)=2
          NU_T1B(K)=0.8_LDP+0.05_LDP*K
        ELSE
          DO I=2,ND-1
            IF(K .GT. 1)THEN
              NU_KUP(I,K)=K-1
              NU_CUPU(I,K)=0.35_LDP+0.02_LDP*I
              NU_CUP(I,K)=NU_CUPU(I,K)
            END IF
          END DO
          IF(K .GT. 1)THEN
            NU_KUP(1,K)=K-1
            NU_CUPU(1,K)=0.55_LDP+0.01_LDP*K
            NU_CUP(1,K)=NU_CUPU(1,K)
          END IF
          DO I=1,ND-1
            IF(K .GT. 1)THEN
              NU_KUH(I,K)=K-1
            ELSE
              NU_HS(I,K)=0.0_LDP
            END IF
          END DO
        END IF
      END DO
!
! EPS coefficients: nonzero except for the pure "exact" baseline.
! Midpoint signs alternate with depth; the unbounded (uncapped)
! partners differ on some midpoints so the target residual sees the
! cap effect.
!
      IF(CASE .NE. 'exact')THEN
        DO K=1,NCF
          DO I=1,ND-1
            F0=0.05_LDP+0.01_LDP*I+0.004_LDP*K
            IF(MOD(I,2) .EQ. 0)THEN
              NU_EPSP(I,K)=F0
              NU_EPS(I,K)=0.8_LDP*F0
            ELSE
              NU_EPSP(I,K)=-F0
              NU_EPS(I,K)=-0.8_LDP*F0
            END IF
            NU_EPSPU(I,K)=NU_EPSP(I,K)
            NU_EPSU(I,K)=NU_EPS(I,K)
          END DO
        END DO
!
! A subset of midpoints is "capped": the post-limiter coefficient is
! a fraction of the uncapped one, so IMPLEMENTED (post-cap) and
! TARGET (uncapped) differ there.
!
        IF(CASE .EQ. 'epsmix' .OR. CASE .EQ. 'lagb' .OR.
     1     CASE .EQ. 'him' .OR. CASE .EQ. 'pert' .OR.
     2     CASE .EQ. 'perth' .OR. CASE .EQ. 'mix' .OR.
     3     CASE .EQ. 'nan' .OR. CASE .EQ. 'offg' .OR.
     4     CASE .EQ. 'bnd' .OR. CASE .EQ. 'epscap' .OR.
     5     CASE .EQ. 'skip' .OR. CASE .EQ. 'state' .OR.
     6     CASE .EQ. 'clip')THEN
          DO K=1,NCF
            DO I=1,ND-1
              IF(MOD(I+K,3) .EQ. 0)THEN
                NU_EPSP(I,K)=0.4_LDP*NU_EPSPU(I,K)
                NU_EPS(I,K)=0.4_LDP*NU_EPSU(I,K)
              END IF
            END DO
          END DO
        END IF
      END IF
!
! Fold the implemented depth-upwind EPS terms and the optional
! MOM_H_IMPL diagonal mass into the captured tridiagonal, exactly
! as the sweep's assembly does.
!
      DO K=1,NCF
        DO I=2,ND-1
          IF(NU_EPS(I-1,K) .GE. 0.0_LDP)THEN
            NU_TA(I,K)=NU_TA(I,K)-2.0_LDP*NU_EPS(I-1,K)
          ELSE
            NU_TB(I,K)=NU_TB(I,K)-2.0_LDP*NU_EPS(I-1,K)
          END IF
          IF(NU_EPS(I,K) .GE. 0.0_LDP)THEN
            NU_TB(I,K)=NU_TB(I,K)+2.0_LDP*NU_EPS(I,K)
          ELSE
            NU_TC(I,K)=NU_TC(I,K)+2.0_LDP*NU_EPS(I,K)
          END IF
          IF(MOM_H_IMPL)NU_TB(I,K)=NU_TB(I,K)+
     1      ABS(NU_VB(I,K))+ABS(NU_VC(I,K))
        END DO
        IF(NU_OUTBC(K) .GT. 1)THEN
          IF(NU_EPS(1,K) .GE. 0.0_LDP)THEN
            NU_TB(1,K)=NU_TB(1,K)-2.0_LDP*NU_EPS(1,K)/NU_T1B(K)
          ELSE
            NU_TC(1,K)=NU_TC(1,K)-2.0_LDP*NU_EPS(1,K)/NU_T1B(K)
          END IF
        END IF
      END DO
!
! Case mutations.
!
      IF(CASE .EQ. 'lagb')THEN
        DO K=1,NCF
          DO I=2,ND-1
            IF(MOD(I+K,4) .EQ. 0)NU_LAGSCL(I,K)=0.30_LDP
          END DO
        END DO
      ELSE IF(CASE .EQ. 'clip')THEN
        MOM_XM_POS=.TRUE.
      ELSE IF(CASE .EQ. 'mix' .OR. CASE .EQ. 'offg' .OR.
     1        CASE .EQ. 'nan' .OR. CASE .EQ. 'skip')THEN
!
! Arbitrary graph: 2-cycles, chains, off-grid reads, stale cols.
!
        NU_KUP(2,2)=3; NU_CUP(2,2)=0.6_LDP; NU_CUPU(2,2)=0.6_LDP
        NU_KUP(2,3)=2; NU_CUP(2,3)=0.5_LDP; NU_CUPU(2,3)=0.5_LDP
        NU_KUP(3,6)=4; NU_CUP(3,6)=0.4_LDP; NU_CUPU(3,6)=0.4_LDP
        NU_KUP(3,7)=6; NU_CUP(3,7)=0.4_LDP; NU_CUPU(3,7)=0.4_LDP
        NU_KUH(1,9)=0
        NU_KUH(2,5)=0
        IF(CASE .EQ. 'skip')THEN
          NU_OUTBC(9)=-1
          NU_OUTBC(10)=-1
        END IF
        IF(CASE .EQ. 'nan')THEN
          J_NEXT_FRQ(3,4)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
        END IF
      ELSE IF(CASE .EQ. 'epscap')THEN
! Boundary row where the EPS cap must fire: EPSP(1)>0 and the
! accumulated positive balance XMP small.
        K=6
        NU_OUTBC(K)=2
        NU_T1B(K)=0.9_LDP
        NU_EPSP(1,K)=0.30_LDP
        NU_EPSPU(1,K)=0.30_LDP
        NU_HS(1,K)=0.0_LDP
        NU_KUP(1,K)=0
        NU_CUP(1,K)=0.0_LDP
        NU_CUPU(1,K)=0.0_LDP
      END IF
!
! Manufacture the committed field. For every case the driver
! computes the implemented RHS exactly (NU_SRC per row, H per
! column in upwind order) -- the "exact" cases therefore return
! residuals at round-off, and the perturbation/stabiliser cases
! isolate exactly the intended defect. Cases with deliberately
! non-exact structure (mix/offg/nan/skip) skip the manufacture and
! keep the generic field.
!
      IF(CASE .EQ. 'exact' .OR. CASE .EQ. 'epsmix' .OR.
     1   CASE .EQ. 'lagb' .OR. CASE .EQ. 'him' .OR.
     2   CASE .EQ. 'pert' .OR. CASE .EQ. 'perth' .OR.
     3   CASE .EQ. 'clip' .OR. CASE .EQ. 'bnd' .OR.
     4   CASE .EQ. 'state')THEN
!
! H field first: columns in the order each upwind index resolves
! (expanding columns ascending, compressing columns descending).
! Mixed-direction upwind graphs can ring (here h(.,8) reads h(.,10)
! while h(.,9..10) read h(.,8)), so both passes are iterated until
! the map contracts to the fixed point.
!
        DO ITER=1,60
        DO K=1,NCF
          IF(K .GE. 5 .AND. K .LE. 8)CYCLE
          DO I=1,ND-1
            THB=0.0_LDP
            U=NU_KUH(I,K)
            IF(U .GE. 1 .AND. U .LE. NCF)
     1        THB=NU_HS(I,K)*H_NEXT_FRQ(I,U)
            T5=0.0_LDP
            IF(NU_EPSP(I,K) .NE. 0.0_LDP)THEN
              IF(NU_EPSP(I,K) .GE. 0.0_LDP)THEN
                D=I; U=K-1
              ELSE
                D=I+1; U=0
                IF(K+MOM_NM_J(D) .LE. NCF)U=K+MOM_NM_J(D)
              END IF
              IF(U .GE. 1 .AND. U .LE. NCF)
     1          T5=2.0_LDP*NU_EPSP(I,K)*J_NEXT_FRQ(D,U)
            END IF
            T4=0.0_LDP
            IF(NU_EPS(I,K) .NE. 0.0_LDP)THEN
              IF(NU_EPS(I,K) .GE. 0.0_LDP)THEN
                D=I
              ELSE
                D=I+1
              END IF
              T4=2.0_LDP*NU_EPS(I,K)*J_NEXT_FRQ(D,K)
            END IF
            H_NEXT_FRQ(I,K)=NU_HU(I,K)*J_NEXT_FRQ(I+1,K)-
     1        NU_HL(I,K)*J_NEXT_FRQ(I,K)+THB+T5-T4
          END DO
        END DO
        DO K=NCF,1,-1
          IF(K .LT. 5 .OR. K .GT. 8)CYCLE
          DO I=1,ND-1
            THB=0.0_LDP
            U=NU_KUH(I,K)
            IF(U .GE. 1 .AND. U .LE. NCF)
     1        THB=NU_HS(I,K)*H_NEXT_FRQ(I,U)
            T5=0.0_LDP
            IF(NU_EPSP(I,K) .NE. 0.0_LDP)THEN
              IF(NU_EPSP(I,K) .GE. 0.0_LDP)THEN
                D=I; U=K-1
              ELSE
                D=I+1; U=0
                IF(K+MOM_NM_J(D) .LE. NCF)U=K+MOM_NM_J(D)
              END IF
              IF(U .GE. 1 .AND. U .LE. NCF)
     1          T5=2.0_LDP*NU_EPSP(I,K)*J_NEXT_FRQ(D,U)
            END IF
            T4=0.0_LDP
            IF(NU_EPS(I,K) .NE. 0.0_LDP)THEN
              IF(NU_EPS(I,K) .GE. 0.0_LDP)THEN
                D=I
              ELSE
                D=I+1
              END IF
              T4=2.0_LDP*NU_EPS(I,K)*J_NEXT_FRQ(D,K)
            END IF
            H_NEXT_FRQ(I,K)=NU_HU(I,K)*J_NEXT_FRQ(I+1,K)-
     1        NU_HL(I,K)*J_NEXT_FRQ(I,K)+THB+T5-T4
          END DO
        END DO
        END DO
!
! NU_SRC so each implemented row is exact.
!
        DO K=1,NCF
          DO I=1,ND
            XUP=0.0_LDP
            IF(NU_KUP(I,K) .GE. 1 .AND. NU_KUP(I,K) .LE. NCF)
     1        XUP=J_NEXT_FRQ(I,NU_KUP(I,K))
            SCL=NU_LAGSCL(I,K)
            IF(I .EQ. 1)THEN
              IF(NU_OUTBC(K) .LE. 1)THEN
                IF(XUP .NE. 0.0_LDP)THEN
                  NU_CUP(1,K)=(NU_TB(1,K)*J_NEXT_FRQ(1,K)+
     1              NU_TC(1,K)*J_NEXT_FRQ(2,K))/XUP
                ELSE
                  NU_CUP(1,K)=0.0_LDP
                  NU_TB(1,K)=-NU_TC(1,K)*J_NEXT_FRQ(2,K)/
     1              MAX(J_NEXT_FRQ(1,K),1.0E-30_LDP)
                END IF
! No limiter touched this row in the manufacture, so the uncapped
! partner used by the target row matches the used coefficient.
                NU_CUPU(1,K)=NU_CUP(1,K)
              ELSE
                THB=0.0_LDP
                U=NU_KUH(1,K)
                IF(NU_HS(1,K) .NE. 0.0_LDP .AND.
     1             U .GE. 1 .AND. U .LE. NCF)
     2            THB=-NU_HS(1,K)*H_NEXT_FRQ(1,U)/NU_T1B(K)
                T4=0.0_LDP
                IF(NU_EPSP(1,K) .NE. 0.0_LDP)THEN
                  IF(NU_EPSP(1,K) .GE. 0.0_LDP)THEN
                    D=1; U=K-1
                  ELSE
                    D=2; U=0
                    IF(K+MOM_NM_J(D) .LE. NCF)U=K+MOM_NM_J(D)
                  END IF
                  IF(U .GE. 1 .AND. U .LE. NCF)
     1              T4=2.0_LDP*NU_EPSP(1,K)*J_NEXT_FRQ(D,U)/
     2                 NU_T1B(K)
                END IF
!
! The evaluator caps the boundary T5 against the balance it sees,
! XMP = LHS + scl*T5_final. Solve that fixed point so the
! manufactured row is exact whether or not the cap fires.
!
                T5=T4
                IF(MOM_EPS_LIM .AND. T4 .GT. 0.0_LDP)THEN
                  XMP=NU_TB(1,K)*J_NEXT_FRQ(1,K)+
     1              NU_TC(1,K)*J_NEXT_FRQ(2,K)
                  IF(T4 .LE. XMP+SCL*T4)THEN
                    T5=T4
                  ELSE IF(XMP .LE. 0.0_LDP)THEN
                    T5=0.0_LDP
                  ELSE IF(SCL .LT. 1.0_LDP)THEN
                    T5=MIN(T4,XMP/(1.0_LDP-SCL))
                  END IF
                END IF
                NU_SRC(1,K)=NU_TB(1,K)*J_NEXT_FRQ(1,K)+
     1            NU_TC(1,K)*J_NEXT_FRQ(2,K)-NU_CUP(1,K)*XUP-
     2            SCL*THB+SCL*T5
              END IF
            ELSE IF(I .EQ. ND)THEN
              NU_SRC(ND,K)=NU_TA(ND,K)*J_NEXT_FRQ(ND-1,K)+
     1          NU_TB(ND,K)*J_NEXT_FRQ(ND,K)
            ELSE
              THI=0.0_LDP
              U=NU_KUH(I-1,K)
              IF(U .GE. 1 .AND. U .LE. NCF)
     1          THI=THI+NU_VB(I,K)*H_NEXT_FRQ(I-1,U)
              U=NU_KUH(I,K)
              IF(U .GE. 1 .AND. U .LE. NCF)
     1          THI=THI+NU_VC(I,K)*H_NEXT_FRQ(I,U)
              TE1=0.0_LDP; TE2=0.0_LDP
              IF(NU_EPSP(I,K) .NE. 0.0_LDP)THEN
                IF(NU_EPSP(I,K) .GE. 0.0_LDP)THEN
                  D=I; U=K-1
                ELSE
                  D=I+1; U=0
                  IF(K+MOM_NM_J(D) .LE. NCF)U=K+MOM_NM_J(D)
                END IF
                IF(U .GE. 1 .AND. U .LE. NCF)
     1            TE1=2.0_LDP*NU_EPSP(I,K)*J_NEXT_FRQ(D,U)
              END IF
              IF(NU_EPSP(I-1,K) .NE. 0.0_LDP)THEN
                IF(NU_EPSP(I-1,K) .GE. 0.0_LDP)THEN
                  D=I-1; U=K-1
                ELSE
                  D=I; U=0
                  IF(K+MOM_NM_J(D) .LE. NCF)U=K+MOM_NM_J(D)
                END IF
                IF(U .GE. 1 .AND. U .LE. NCF)
     1            TE2=2.0_LDP*NU_EPSP(I-1,K)*J_NEXT_FRQ(D,U)
              END IF
              NU_SRC(I,K)=NU_TA(I,K)*J_NEXT_FRQ(I-1,K)+
     1          NU_TB(I,K)*J_NEXT_FRQ(I,K)+
     2          NU_TC(I,K)*J_NEXT_FRQ(I+1,K)-NU_CUP(I,K)*XUP-
     3          SCL*(THI+TE1-TE2)
            END IF
          END DO
        END DO
      END IF
!
! clip case: after exact manufacture, rebuild row (3,4) so the raw
! xm is negative while the LHS is zero -- the implemented clip
! lands the row exactly on xm=0, the target sees the debit.
!
      IF(CASE .EQ. 'clip')THEN
        I=3; K=4
        XUP=0.0_LDP
        IF(NU_KUP(I,K) .GE. 1 .AND. NU_KUP(I,K) .LE. NCF)
     1    XUP=J_NEXT_FRQ(I,NU_KUP(I,K))
        SCL=NU_LAGSCL(I,K)
        THI=0.0_LDP
        U=NU_KUH(I-1,K)
        IF(U .GE. 1 .AND. U .LE. NCF)
     1    THI=THI+NU_VB(I,K)*H_NEXT_FRQ(I-1,U)
        U=NU_KUH(I,K)
        IF(U .GE. 1 .AND. U .LE. NCF)
     1    THI=THI+NU_VC(I,K)*H_NEXT_FRQ(I,U)
        TE1=0.0_LDP; TE2=0.0_LDP
        IF(NU_EPSP(I,K) .NE. 0.0_LDP)THEN
          IF(NU_EPSP(I,K) .GE. 0.0_LDP)THEN
            D=I; U=K-1
          ELSE
            D=I+1; U=0
            IF(K+MOM_NM_J(D) .LE. NCF)U=K+MOM_NM_J(D)
          END IF
          IF(U .GE. 1 .AND. U .LE. NCF)
     1      TE1=2.0_LDP*NU_EPSP(I,K)*J_NEXT_FRQ(D,U)
        END IF
        IF(NU_EPSP(I-1,K) .NE. 0.0_LDP)THEN
          IF(NU_EPSP(I-1,K) .GE. 0.0_LDP)THEN
            D=I-1; U=K-1
          ELSE
            D=I; U=0
            IF(K+MOM_NM_J(D) .LE. NCF)U=K+MOM_NM_J(D)
          END IF
          IF(U .GE. 1 .AND. U .LE. NCF)
     1      TE2=2.0_LDP*NU_EPSP(I-1,K)*J_NEXT_FRQ(D,U)
        END IF
        NU_TB(I,K)=-(NU_TA(I,K)*J_NEXT_FRQ(I-1,K)+
     1    NU_TC(I,K)*J_NEXT_FRQ(I+1,K))/J_NEXT_FRQ(I,K)
        NU_SRC(I,K)=-2.0_LDP-NU_CUP(I,K)*XUP-
     1    SCL*(THI+TE1-TE2)
      END IF
!
! pert/perth: manufacture the exact system, then perturb one cell.
!
      IF(CASE .EQ. 'pert')J_NEXT_FRQ(3,5)=J_NEXT_FRQ(3,5)*1.4_LDP+
     1   0.30_LDP
      IF(CASE .EQ. 'perth')H_NEXT_FRQ(2,6)=H_NEXT_FRQ(2,6)+
     1   0.25_LDP
!
      NU_AUG_RDY=.TRUE.
      CKSUM=0.0_LDP
      DO K=1,NCF
        DO I=1,ND
          CKSUM=CKSUM+J_NEXT_FRQ(I,K)+H_NEXT_FRQ(I,K)
        END DO
      END DO
      WRITE(*,'(A,ES25.16)')'PRECK=',CKSUM
!
      CALL MOM_AUG_RESID_V1(6)
!
      CKSUM=0.0_LDP
      DO K=1,NCF
        DO I=1,ND
          CKSUM=CKSUM+J_NEXT_FRQ(I,K)+H_NEXT_FRQ(I,K)
        END DO
      END DO
      WRITE(*,'(A,ES25.16)')'POSTCK=',CKSUM
!
      WRITE(*,'(A,2(ES25.16))')'AR_R2M=',MOM_AR_R2M(1),MOM_AR_R2M(2)
      WRITE(*,'(A,2(ES25.16))')'AR_R2H=',MOM_AR_R2H(1),MOM_AR_R2H(2)
      WRITE(*,'(A,2(ES25.16))')'AR_S2M=',MOM_AR_S2M(1),MOM_AR_S2M(2)
      WRITE(*,'(A,2(ES25.16))')'AR_S2H=',MOM_AR_S2H(1),MOM_AR_S2H(2)
      WRITE(*,'(A,2(ES25.16))')'AR_MXM=',MOM_AR_MXM(1),MOM_AR_MXM(2)
      WRITE(*,'(A,2(ES25.16))')'AR_MXH=',MOM_AR_MXH(1),MOM_AR_MXH(2)
      WRITE(*,'(A,4(I8))')'AR_CELLM=',MOM_AR_IM(1),MOM_AR_KM(1),
     1   MOM_AR_IM(2),MOM_AR_KM(2)
      WRITE(*,'(A,4(I8))')'AR_CELLH=',MOM_AR_IH(1),MOM_AR_KH(1),
     1   MOM_AR_IH(2),MOM_AR_KH(2)
      WRITE(*,'(A,8(I10))')'AR_CNT=',MOM_AR_NCLIP,MOM_AR_NLAGB,
     1   MOM_AR_NEPSL,MOM_AR_NOFF,MOM_AR_NNAN,MOM_AR_NHCLP,
     2   MOM_AR_NSKIP,MOM_AR_NBNDR
      WRITE(*,'(A,ES25.16)')'AR_RBN2=',MOM_AR_RBN2
!
      OPEN(UNIT=310,FILE=TRIM(TAG)//'_ar.txt',STATUS='REPLACE',
     1  IOSTAT=IOS)
      IF(IOS .NE. 0)STOP 'DRV_MOM_AUGRES: dump open failed'
      WRITE(310,'(A,2I5)')'# SHAPE',ND,NCF
      WRITE(310,'(A,3L3)')'# FLAGS',MOM_XM_POS,MOM_EPS_LIM,
     1  MOM_H_IMPL
      WRITE(310,'(A)')'# JNEXT'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,J_NEXT_FRQ(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# HNEXT'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,H_NEXT_FRQ(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# TA'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_TA(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# TB'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_TB(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# TC'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_TC(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# CUP'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_CUP(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# CUPU'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_CUPU(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# KUP'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,I8)')I,K,NU_KUP(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# HU'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,NU_HU(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# HL'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,NU_HL(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# HS'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,NU_HS(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# VB'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_VB(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# VC'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_VC(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# EPS'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,NU_EPS(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# EPSU'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,NU_EPSU(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# EPSP'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,NU_EPSP(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# EPSPU'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,ES25.16)')I,K,NU_EPSPU(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# LAGSCL'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_LAGSCL(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# SRC'
      DO K=1,NCF
        DO I=1,ND
          WRITE(310,'(2I5,ES25.16)')I,K,NU_SRC(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# T1B'
      DO K=1,NCF
        WRITE(310,'(I5,ES25.16)')K,NU_T1B(K)
      END DO
      WRITE(310,'(A)')'# KUH'
      DO K=1,NCF
        DO I=1,ND-1
          WRITE(310,'(2I5,I8)')I,K,NU_KUH(I,K)
        END DO
      END DO
      WRITE(310,'(A)')'# OUTBC'
      DO K=1,NCF
        WRITE(310,'(I5,I8)')K,NU_OUTBC(K)
      END DO
      WRITE(310,'(A)')'# HCHK'
      DO K=1,NCF
        WRITE(310,'(I5,I8)')K,NU_HCHK(K)
      END DO
      WRITE(310,'(A)')'# NMJ'
      DO I=1,ND
        WRITE(310,'(I5,I8)')I,MOM_NM_J(I)
      END DO
      WRITE(310,'(A)')'# NMH'
      DO I=1,ND
        WRITE(310,'(I5,I8)')I,MOM_NM_H(I)
      END DO
      CLOSE(310)
!
      WRITE(*,'(A)')'RESULT_OK'
      END
