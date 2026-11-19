! drv_knop_spv.f90 - runtime ABI/purity test for KN_ORD_SPV_IMPL.
! Calls the module-scope worker with L_CNT=.TRUE. (production
! entry) and L_CNT=.FALSE. (shadow entry) across the full branch
! matrix and prints KEY=val lines for the python harness.
!
! Operand integrity: checksums (shape + SUM) of every input array
! are printed before and after the matrix; KN_DBG_NDON is compared
! elementwise against a saved copy so the single expected change
! (the production expanding L>1 increment) is verified explicitly
! and every other counter entry is proven unchanged.
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
      INTEGER I,L,R,NC,ND0,NDP,NDS,NDIF
      INTEGER NDON0(NIX,NCFX)
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
! Pre-call operand state: shape plus SUM for every array the
! worker reads, including the module-scope operands.
!
      NDON0=KN_DBG_NDON
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'PRE ','SR',SIZE(SR_A,1),SIZE(SR_A,2),'SUM',SUM(SR_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'PRE ','Q',SIZE(Q_A,1),SIZE(Q_A,2),'SUM',SUM(Q_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'PRE ','FLD',SIZE(FLD_A,1),SIZE(FLD_A,2),'SUM',SUM(FLD_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'PRE ','DRF',SIZE(DRF_A),0,'SUM',SUM(DRF_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'PRE ','WINX',SIZE(WINX_A),0,'SUM',SUM(WINX_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'PRE ','JPR1',SIZE(JPR1_A),0,'SUM',SUM(JPR1_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'PRE ','JPR2',SIZE(JPR2_A),0,'SUM',SUM(JPR2_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'PRE ','IKC1',SIZE(IKC1_A),0,'SUM',SUM(IKC1_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'PRE ','IKC2',SIZE(IKC2_A),0,'SUM',SUM(IKC2_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'PRE ','GAM',SIZE(GAM,1),SIZE(GAM,2),'SUM',SUM(GAM)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'PRE ','DLN',SIZE(KN_DLN_N),0,'SUM',SUM(KN_DLN_N)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'PRE ','NDON',SIZE(KN_DBG_NDON,1),SIZE(KN_DBG_NDON,2), &
        'SUM',SUM(KN_DBG_NDON)
!
! Branch matrix: shadow first (never counts), then production.
! Per-case state reset: GAM is indexed (depth, ray) and the
! worker reads GAM(R,LS) with LS=1, so all coupling control
! lives in column 1; IKC/FM flags also persist otherwise.
!
      DO NC=1,10
        CN='        '
        R=2; L=3
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
          CN='gamnj2'
          R=3; L=2
          GAM(3,1)=-0.4_LDP
          IKC1_A(3)=0
        ELSE IF(NC .EQ. 7)THEN
          CN='gamnn'
          R=4; L=2
          GAM(4,1)=-0.4_LDP
          IKC1_A(4)=0
          IKC2_A(4)=0
        ELSE IF(NC .EQ. 8)THEN
          CN='gamnb'
          R=2; L=6
          GAM(2,1)=-0.4_LDP
        ELSE IF(NC .EQ. 9)THEN
          CN='fmtch'
          R=2; L=4
          GAM(2,1)=-0.4_LDP
          FG_FREQ_MATCH=.TRUE.
        ELSE
          CN='fmtcb'
          R=2; L=5
          GAM(2,1)=-0.4_LDP
          FG_FREQ_MATCH=.TRUE.
        END IF
        ND0=KN_DBG_NDON(R,L)
        CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
          WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVF,.FALSE.)
        NDS=KN_DBG_NDON(R,L)-ND0
        CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
          WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVT,.TRUE.)
        NDP=KN_DBG_NDON(R,L)-ND0-NDS
        BT=TRANSFER(SVT,BT)
        BF=TRANSFER(SVF,BF)
        BEQ=BT .EQ. BF
        WRITE(*,'(A,A8,1X,A,I2,1X,A,I2,1X,A,ES22.14,1X,A,ES22.14, &
          1X,A,L1,1X,A,I3,1X,A,I3)') &
          'SPCASE ',CN,'R',R,'L',L,'SVT',SVT,'SVF',SVF, &
          'BEQ',BEQ,'NDP',NDP,'NDS',NDS
      END DO
!
! Post-call operand state. GAM is driver-controlled state (the
! worker only reads it), so it is reset to the pre-call fill
! first - the POST checksum then proves the worker wrote nothing.
!
      GAM=0.0_LDP
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'POST','SR',SIZE(SR_A,1),SIZE(SR_A,2),'SUM',SUM(SR_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'POST','Q',SIZE(Q_A,1),SIZE(Q_A,2),'SUM',SUM(Q_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'POST','FLD',SIZE(FLD_A,1),SIZE(FLD_A,2),'SUM',SUM(FLD_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'POST','DRF',SIZE(DRF_A),0,'SUM',SUM(DRF_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'POST','WINX',SIZE(WINX_A),0,'SUM',SUM(WINX_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'POST','JPR1',SIZE(JPR1_A),0,'SUM',SUM(JPR1_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'POST','JPR2',SIZE(JPR2_A),0,'SUM',SUM(JPR2_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'POST','IKC1',SIZE(IKC1_A),0,'SUM',SUM(IKC1_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,I10)') &
        'POST','IKC2',SIZE(IKC2_A),0,'SUM',SUM(IKC2_A)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'POST','GAM',SIZE(GAM,1),SIZE(GAM,2),'SUM',SUM(GAM)
      WRITE(*,'(A,A6,1X,I4,1X,I4,1X,A,ES22.14)') &
        'POST','DLN',SIZE(KN_DLN_N),0,'SUM',SUM(KN_DLN_N)
!
! KN_DBG_NDON elementwise diff against the pre-call copy. The
! only permitted change is +1 at (2,3) from the gamp production
! call.
!
      NDIF=0
      DO I=1,NIX
        DO L=1,NCFX
          IF(KN_DBG_NDON(I,L) .NE. NDON0(I,L))THEN
            NDIF=NDIF+1
            WRITE(*,'(A,I2,1X,I2,1X,I4)')'NDELT ',I,L, &
              KN_DBG_NDON(I,L)-NDON0(I,L)
          END IF
        END DO
      END DO
      WRITE(*,'(A,I4)')'NDCHG ',NDIF
!
! Unallocated-counter path, run last: nested guard must keep
! both entries safe with KN_DBG_NDON deallocated.
!
      DEALLOCATE(KN_DBG_NDON)
      R=2; L=3
      GAM=0.0_LDP
      GAM(2,1)=0.5_LDP
      CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
        WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVF,.FALSE.)
      CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
        WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVT,.TRUE.)
      BT=TRANSFER(SVT,BT)
      BF=TRANSFER(SVF,BF)
      BEQ=BT .EQ. BF
      WRITE(*,'(A,A8,1X,A,I2,1X,A,I2,1X,A,ES22.14,1X,A,ES22.14, &
        1X,A,L1,1X,A,I3,1X,A,I3)') &
        'SPCASE ','unall','R',R,'L',L,'SVT',SVT,'SVF',SVF, &
        'BEQ',BEQ,'NDP',-1,'NDS',-1
      END
