! drv_knop_spv.f90 - runtime ABI/purity test for KN_ORD_SPV_IMPL.
! Calls the module-scope worker with L_CNT=.TRUE. (production
! entry) and L_CNT=.FALSE. (shadow entry) across the full branch
! matrix and prints KEY=val lines for the python harness.
!
! Operand integrity: every input array is snapshotted immediately
! before EACH worker invocation and compared elementwise right
! after - NWS/NWP count changed elements per shadow/production
! call. GAM is included in its branch-configured state (the
! driver sets it per case; the worker must never write it).
! KN_DBG_NDON is the single permitted mutable input-side array:
! it is diffed elementwise against a saved copy so the one
! expected increment (+1 at (2,3), production expanding L>1) is
! verified explicitly and every other entry proven unchanged.
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
      INTEGER I,L,R,NC,ND0,NDP,NDS,NDIF,NWS,NWP
      INTEGER NDON0(NIX,NCFX)
      REAL(KIND=LDP) SR0(NIX,NCFX),Q0(NIX,NCFX),FLD0(NIX,NCFX)
      REAL(KIND=LDP) DRF0(NIX),WINX0(NIX)
      REAL(KIND=LDP) GAM0(NIX,NCFX),DLN0(NCFX)
      INTEGER JPR10(NIX),JPR20(NIX),IKC10(NIX),IKC20(NIX)
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
! Matrix-level pre/post checksums are retained as a secondary
! guard; the per-invocation elementwise compare is the primary
! operand-integrity evidence.
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
        CALL SNAP_INP()
        CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
          WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVF,.FALSE.)
        NWS=NWRT_INP()
        NDS=KN_DBG_NDON(R,L)-ND0
        CALL SNAP_INP()
        CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
          WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVT,.TRUE.)
        NWP=NWRT_INP()
        NDP=KN_DBG_NDON(R,L)-ND0-NDS
        BT=TRANSFER(SVT,BT)
        BF=TRANSFER(SVF,BF)
        BEQ=BT .EQ. BF
        WRITE(*,'(A,A8,1X,A,I2,1X,A,I2,1X,A,ES22.14,1X,A,ES22.14, &
          1X,A,L1,1X,A,I3,1X,A,I3,1X,A,I3,1X,A,I3)') &
          'SPCASE ',CN,'R',R,'L',L,'SVT',SVT,'SVF',SVF, &
          'BEQ',BEQ,'NDP',NDP,'NDS',NDS,'NWS',NWS,'NWP',NWP
      END DO
!
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
      CALL SNAP_INP()
      CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
        WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVF,.FALSE.)
      NWS=NWRT_INP()
      CALL SNAP_INP()
      CALL KN_ORD_SPV_IMPL(SR_A,Q_A,DRF_A,JPR1_A,JPR2_A, &
        WINX_A,IKC1_A,IKC2_A,1,NCFX,FLD_A,R,L,SVT,.TRUE.)
      NWP=NWRT_INP()
      BT=TRANSFER(SVT,BT)
      BF=TRANSFER(SVF,BF)
      BEQ=BT .EQ. BF
      WRITE(*,'(A,A8,1X,A,I2,1X,A,I2,1X,A,ES22.14,1X,A,ES22.14, &
        1X,A,L1,1X,A,I3,1X,A,I3,1X,A,I3,1X,A,I3)') &
        'SPCASE ','unall','R',R,'L',L,'SVT',SVT,'SVF',SVF, &
        'BEQ',BEQ,'NDP',-1,'NDS',-1,'NWS',NWS,'NWP',NWP
!
      CONTAINS
!
      SUBROUTINE SNAP_INP()
      SR0=SR_A; Q0=Q_A; FLD0=FLD_A
      DRF0=DRF_A; WINX0=WINX_A
      JPR10=JPR1_A; JPR20=JPR2_A; IKC10=IKC1_A; IKC20=IKC2_A
      GAM0=GAM; DLN0=KN_DLN_N
      END SUBROUTINE SNAP_INP
!
      FUNCTION NWRT_INP()
      INTEGER NWRT_INP
      NWRT_INP=COUNT(SR_A .NE. SR0)+COUNT(Q_A .NE. Q0)+ &
        COUNT(FLD_A .NE. FLD0)+COUNT(DRF_A .NE. DRF0)+ &
        COUNT(WINX_A .NE. WINX0)+COUNT(JPR1_A .NE. JPR10)+ &
        COUNT(JPR2_A .NE. JPR20)+COUNT(IKC1_A .NE. IKC10)+ &
        COUNT(IKC2_A .NE. IKC20)+COUNT(GAM .NE. GAM0)+ &
        COUNT(KN_DLN_N .NE. DLN0)
      END FUNCTION NWRT_INP
!
      END
