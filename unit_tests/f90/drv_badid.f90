! Driver for the round-2 phase-2 FG diagnostics in
! plane/fg_j_cmf_nonmono_v1.f:
!
!   * FG_DSYM_VAL -- symmetric store-drift formula (analytic values).
!   * FG_SYM_ACC  -- the single symmetric-max accumulator: value,
!     argmax cell and channel bookkeeping.
!   * GET_FG_SW_SYM -- the snapshot exporter.
!   * KN_I8SORT   -- the INTEGER*8 heapsort used by the census merge
!     (INDEXX sorts LDP and would lose bits above 2^53).
!   * KN_BADID_MERGE_V1 -- the serial census merge/diff: identical,
!     disjoint and partially overlapping censuses, the per-channel
!     split, and the partial-census branch a failed ray channel
!     (N=-1 sentinel) produces.
!
! Prints KEY=value lines for test_fg_sw_sym.py. The merge runs on
! hand-filled KN_BADRAY lists so no field solve is needed; the
! production append path is exercised end-to-end by drv_knop_ray
! cases badid/bido.
!
      PROGRAM DRV_BADID
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      IMPLICIT NONE
      INTEGER*8 NBDQ(2)
      LOGICAL MSKQ(2)
      INTEGER*8 BNC,BNP,BINT,BNEW,BGON,BNCC(2),BNIC(2),BNUN
      INTEGER*8 IDS(6)
      REAL(KIND=LDP) BJAC,SYM
      INTEGER SK,SR,SL,SC,I,IOS
      LOGICAL SRTD
!
! ---- 1. symmetric drift values ----------------------------------
!
      WRITE(*,'(A,ES16.8)')'DSYM_A=', &
        FG_DSYM_VAL(2.0_LDP,1.0_LDP,1.0E-4_LDP)
      WRITE(*,'(A,ES16.8)')'DSYM_B=', &
        FG_DSYM_VAL(1.0E10_LDP,1.0_LDP,1.0E-4_LDP)
      WRITE(*,'(A,ES16.8)')'DSYM_C=', &
        FG_DSYM_VAL(1.0E-6_LDP,1.0E-30_LDP,1.0E-4_LDP)
      WRITE(*,'(A,ES16.8)')'DSYM_D=', &
        FG_DSYM_VAL(-5.0_LDP,1.0_LDP,1.0E-4_LDP)
      WRITE(*,'(A,ES16.8)')'DSYM_E=', &
        FG_DSYM_VAL(0.0_LDP,0.0_LDP,1.0E-4_LDP)
      WRITE(*,'(A,ES16.8)')'DSYM_F=', &
        FG_DSYM_VAL(1.0E-8_LDP,1.0E-9_LDP,1.0E-4_LDP)
!
! ---- 2. symmetric max accumulator -------------------------------
!
      FG_GS_MAXSYM=0.0_LDP
      FG_GS_SYMK=0; FG_GS_SYMR=0; FG_GS_SYML=0; FG_GS_SYMC=0
      CALL FG_SYM_ACC(0.3_LDP,5,2,100,1)
      CALL FG_SYM_ACC(0.7_LDP,9,3,200,2)
      CALL FG_SYM_ACC(0.4_LDP,1,1,50,1)
      WRITE(*,'(A,ES16.8)')'ACC_SYM=',FG_GS_MAXSYM
      WRITE(*,'(A,4I6)')'ACC_ARG=',FG_GS_SYMK,FG_GS_SYMR, &
        FG_GS_SYML,FG_GS_SYMC
!
! Snapshot exporter echoes module state.
!
      FG_SW_MAXSYM=0.42_LDP
      FG_SW_SYMK=4; FG_SW_SYMR=7; FG_SW_SYML=88; FG_SW_SYMC=1
      CALL GET_FG_SW_SYM(SYM,SK,SR,SL,SC)
      WRITE(*,'(A,ES16.8)')'GET_SYM=',SYM
      WRITE(*,'(A,4I6)')'GET_ARG=',SK,SR,SL,SC
!
! ---- 3. int8 heapsort -------------------------------------------
!
      IDS=(/5_8,1_8,4_8,1_8,3_8,2_8/)
      CALL KN_I8SORT(IDS)
      SRTD=.TRUE.
      DO I=2,6
        IF(IDS(I) .LT. IDS(I-1))SRTD=.FALSE.
      END DO
      WRITE(*,'(A,L2)')'SORTED=',SRTD
!
! ---- 4. census merge --------------------------------------------
!
      CALL SET_FG_BADID(.FALSE.)
      CALL SET_FG_BADID(.TRUE.)
      ALLOCATE(KN_BADRAY(2,2),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT_FAIL badray_alloc'
        STOP
      END IF
      KN_BADRAY%N=0_8
!
! Ray 1: M channel 3 cells (l=5 i=2, l=5 i=3, l=9 i=1),
! P channel 1 cell (l=7 i=4). Ray 2: M channel 1 cell (l=5 i=2).
!
      ALLOCATE(KN_BADRAY(1,1)%V(3))
      KN_BADRAY(1,1)%V(1)=PK(1,5,2,1)
      KN_BADRAY(1,1)%V(2)=PK(1,5,3,1)
      KN_BADRAY(1,1)%V(3)=PK(1,9,1,1)
      KN_BADRAY(1,1)%N=3_8
      ALLOCATE(KN_BADRAY(1,2)%V(1))
      KN_BADRAY(1,2)%V(1)=PK(1,7,4,2)
      KN_BADRAY(1,2)%N=1_8
      ALLOCATE(KN_BADRAY(2,1)%V(1))
      KN_BADRAY(2,1)%V(1)=PK(2,5,2,1)
      KN_BADRAY(2,1)%N=1_8
      ALLOCATE(KN_BADRAY(2,2)%V(1))
      KN_BADRAY(2,2)%N=0_8
      NBDQ=(/4_8,1_8/)
      MSKQ=(/.TRUE.,.TRUE./)
!
! First census: no previous -> JAC=-1.
!
      CALL KN_BADID_MERGE_V1(NBDQ,MSKQ,2,6)
      CALL GET_FG_BADID(BNC,BNP,BINT,BNEW,BGON,BJAC,BNCC,BNIC,BNUN)
      WRITE(*,'(A,I14)')'M1_NC=',BNC
      WRITE(*,'(A,I14)')'M1_NP=',BNP
      WRITE(*,'(A,I14)')'M1_INT=',BINT
      WRITE(*,'(A,I14)')'M1_NEW=',BNEW
      WRITE(*,'(A,I14)')'M1_GONE=',BGON
      WRITE(*,'(A,ES16.8)')'M1_JAC=',BJAC
      WRITE(*,'(A,I14)')'M1_MCUR=',BNCC(1)
      WRITE(*,'(A,I14)')'M1_PCUR=',BNCC(2)
      WRITE(*,'(A,I14)')'M1_UNTR=',BNUN
!
! Second census, identical lists -> JAC=1, all persistent.
!
      CALL KN_BADID_MERGE_V1(NBDQ,MSKQ,2,6)
      CALL GET_FG_BADID(BNC,BNP,BINT,BNEW,BGON,BJAC,BNCC,BNIC,BNUN)
      WRITE(*,'(A,I14)')'M2_NC=',BNC
      WRITE(*,'(A,I14)')'M2_INT=',BINT
      WRITE(*,'(A,I14)')'M2_NEW=',BNEW
      WRITE(*,'(A,I14)')'M2_GONE=',BGON
      WRITE(*,'(A,ES16.8)')'M2_JAC=',BJAC
      WRITE(*,'(A,I14)')'M2_MINT=',BNIC(1)
      WRITE(*,'(A,I14)')'M2_PINT=',BNIC(2)
!
! Third census, ray 2's cell moved to a different frequency ->
! inter=4, new=1, gone=1, JAC=4/6.
!
      KN_BADRAY(2,1)%V(1)=PK(2,11,2,1)
      CALL KN_BADID_MERGE_V1(NBDQ,MSKQ,2,6)
      CALL GET_FG_BADID(BNC,BNP,BINT,BNEW,BGON,BJAC,BNCC,BNIC,BNUN)
      WRITE(*,'(A,I14)')'M3_NC=',BNC
      WRITE(*,'(A,I14)')'M3_INT=',BINT
      WRITE(*,'(A,I14)')'M3_NEW=',BNEW
      WRITE(*,'(A,I14)')'M3_GONE=',BGON
      WRITE(*,'(A,ES16.8)')'M3_JAC=',BJAC
!
! Fourth census, ray 2 M channel marked untracked (the sentinel an
! allocation failure leaves) -> partial census: ray 2 contributes
! nothing, NUNTR=1.
!
      KN_BADRAY(2,1)%N=-1_8
      CALL KN_BADID_MERGE_V1(NBDQ,MSKQ,2,6)
      CALL GET_FG_BADID(BNC,BNP,BINT,BNEW,BGON,BJAC,BNCC,BNIC,BNUN)
      WRITE(*,'(A,I14)')'M4_NC=',BNC
      WRITE(*,'(A,I14)')'M4_UNTR=',BNUN
      WRITE(*,'(A,ES16.8)')'M4_JAC=',BJAC
      WRITE(*,'(A)')'RESULT_OK'
!
      CONTAINS
!
      INTEGER*8 FUNCTION PK(LS,L,I,CH)
      INTEGER LS,L,I,CH
      PK=(LS-1_8)*17179869184_8+(L-1_8)*8192_8+(I-1_8)*2_8+ &
        (CH-1_8)
      END FUNCTION PK
!
      END PROGRAM DRV_BADID
