! drv_knop_spv.f90 - runtime ABI/purity test for KN_ORD_SPV_IMPL.
! Calls the module-scope worker with L_CNT=.TRUE. (production
! entry) and L_CNT=.FALSE. (shadow entry) across the full branch
! matrix and prints KEY=val lines for the python harness.
! Verifies: bitwise-identical SPV for the two entries, the
! production donor counter increments exactly on expanding
! L>1 reads, the shadow leaves it untouched, compressing and
! GAM=0 branches never touch it, and the unallocated-counter
! path is safe under bounds checking.
!
      PROGRAM DRV_KNOP_SPV
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      IMPLICIT NONE
      INTEGER, PARAMETER :: NIX=4,NCFX=6
      REAL(KIND=LDP) SR_A(NIX,NCFX),Q_A(NIX,NCFX)
      REAL(KIND=LDP) DRF_A(NIX),WINX_A(NIX)
      REAL(KIND=LDP) FLD_A(NIX,NCFX),SVT,SVF
      INTEGER JPR1_A(NIX),JPR2_A(NIX),IKC1_A(NIX),IKC2_A(NIX)
      INTEGER I,L,R,NC,ND0,NDP,NDS
      INTEGER*8 BT,BF
      LOGICAL BEQ
      CHARACTER*8 CN
      ALLOCATE(GAM(NIX,NCFX),KN_DLN_N(NCFX),KN_DBG_NDON(NIX,NCFX))
      GAM=0.0_LDP
      KN_DLN_N=0.2_LDP
      KN_DBG_NDON=0
      KN_DBG_LS=1
      FG_FREQ_MATCH=.FALSE.
      DO I=1,NIX
        DO L=1,NCFX
          SR_A(I,L)=1.5_LDP+0.01_LDP*I+0.001_LDP*L
          Q_A(I,L)=0.3_LDP+0.01_LDP*I
          FLD_A(I,L)=0.1_LDP*I+0.01_LDP*L
        END DO
        DRF_A(I)=0.4_LDP+0.1_LDP*I
        WINX_A(I)=0.25_LDP
        JPR1_A(I)=MAX(1,I-1)
        JPR2_A(I)=I
        IKC1_A(I)=1
        IKC2_A(I)=1
      END DO
!
! One case per loop pass: shadow first (never counts), then
! production (expected count printed per case).
!
      DO NC=1,10
        CN='        '
        R=2; L=3
!       Per-case state reset: GAM is indexed (depth, ray) and the
!       worker reads GAM(R,LS) with LS=1, so all coupling control
!       lives in column 1; IKC/FM flags also persist otherwise.
        GAM=0.0_LDP
        IKC1_A=1
        IKC2_A=1
        FG_FREQ_MATCH=.FALSE.
        IF(NC .EQ. 1)THEN
          CN='gam0'
        ELSE IF(NC .EQ. 2)THEN
          CN='gamp'
          GAM(2,1)=0.5_LDP
        ELSE IF(NC .EQ. 3)THEN
          CN='gamp1'
          R=2; L=1
          GAM(2,1)=0.5_LDP
        ELSE IF(NC .EQ. 4)THEN
          CN='gamn2'
          R=2; L=2
          GAM(2,1)=-0.4_LDP
        ELSE IF(NC .EQ. 5)THEN
          CN='gamn1'
          R=3; L=2
          GAM(3,1)=-0.4_LDP
          IKC2_A(3)=0
        ELSE IF(NC .EQ. 6)THEN
          CN='gamnn'
          R=4; L=2
          GAM(4,1)=-0.4_LDP
          IKC1_A(4)=0
          IKC2_A(4)=0
        ELSE IF(NC .EQ. 7)THEN
          CN='gamnb'
          R=2; L=6
          GAM(2,1)=-0.4_LDP
        ELSE IF(NC .EQ. 8)THEN
          CN='fmtch'
          R=2; L=4
          GAM(2,1)=-0.4_LDP
          FG_FREQ_MATCH=.TRUE.
        ELSE IF(NC .EQ. 9)THEN
          CN='fmtcb'
          R=2; L=5
          GAM(2,1)=-0.4_LDP
          FG_FREQ_MATCH=.TRUE.
        ELSE
          CN='unall'
          R=2; L=3
          GAM(2,1)=0.5_LDP
          DEALLOCATE(KN_DBG_NDON)
        END IF
        IF(ALLOCATED(KN_DBG_NDON))THEN
          ND0=KN_DBG_NDON(R,L)
        ELSE
          ND0=-1
        END IF
        CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
          WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVF,.FALSE.)
        IF(ALLOCATED(KN_DBG_NDON))THEN
          NDS=KN_DBG_NDON(R,L)-ND0
        ELSE
          NDS=-1
        END IF
        CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
          WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVT,.TRUE.)
        IF(ALLOCATED(KN_DBG_NDON))THEN
          NDP=KN_DBG_NDON(R,L)-ND0-NDS
        ELSE
          NDP=-1
        END IF
        BT=TRANSFER(SVT,BT)
        BF=TRANSFER(SVF,BF)
        BEQ=BT .EQ. BF
        WRITE(*,'(A,A8,1X,A,I2,1X,A,I2,1X,A,ES22.14,1X,A,ES22.14, &
          1X,A,L1,1X,A,I3,1X,A,I3)') &
          'SPCASE ',CN,'R',R,'L',L,'SVT',SVT,'SVF',SVF, &
          'BEQ',BEQ,'NDP',NDP,'NDS',NDS
      END DO
!
! Operand-integrity checksum: the routine must not write its
! array arguments.
!
      WRITE(*,'(A,ES22.14)')'SUMSR ',SUM(SR_A)
      WRITE(*,'(A,ES22.14)')'SUMQ  ',SUM(Q_A)
      WRITE(*,'(A,ES22.14)')'SUMFLD',SUM(FLD_A)
      WRITE(*,'(A,I6)')       'SUMJPR',SUM(JPR1_A)+SUM(JPR2_A)
      WRITE(*,'(A,I6)')       'SUMIKC',SUM(IKC1_A)+SUM(IKC2_A)
      END
