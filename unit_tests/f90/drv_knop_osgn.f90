! Driver for the WP-5 signed-Aitken store relaxation
! (KN_OMEGA_SGN; audit KNOP_REMEDIATION_WP5_2026-11-26.md).
!
! The production machinery is exercised at two levels:
!
!   scalar: KN_OSGN_STEP fed by a one-dimensional lagged map
!   x_{k+1} = x_k + omega_k*(a + lam*(x_k-a) - x_k). For a linear
!   map the estimator must return the fixed-point omega
!   1/(1-lam) from a single (dx,dr) pair, and the relaxed
!   iteration must reach the root.
!
!   vector: the module arrays themselves -- FG_OSGN_ALLOC +
!   column-buffer pokes replicating what the sweep write sites
!   capture, FG_OSGN_COL_FIN at each column close, FG_OSGN_UPDATE
!   at the refresh tail. Two independent slots sharing one lam
!   recover the same 1/(1-lam); a coefficient-fingerprint move
!   must clear the slot history so the next update holds omega;
!   a second UPDATE inside one sweep generation must fall back.
!
! Cases:
!   sc   contracting map lam=0.3 -- converge, omega -> 10/7
!   so   oscillating map lam=-0.8 -- converge, omega -> 5/9
!   sd   expansion under omega0=0.5 (lam=1.5) -- halve once,
!        then ICODE=3 (retry budget spent)
!   sc0  collapsed denominator -- ICODE=0 hold
!   scl  clamp: lam=0.55 -> omega* = 20/9 hits KN_OSGN_MAX;
!        lam=-1.5 against OMIN=0.5 hits the floor
!   vec  two-slot vector map, sweep-level machinery
!   vfp  fingerprint move clears the history (ICODE=0 hold)
!   vgen UPDATE twice in one generation -- second falls back
!
      PROGRAM DRV_KNOP_OSGN
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      IMPLICIT NONE
!
      INTEGER NCASE,I,J,ICODE,NRET,NSTEP
      INTEGER IS,IC
      REAL(KIND=LDP) LAM,RT,OM,X,XP,R,RN,DX,DR
      REAL(KIND=LDP) OMIN2,OMAX2,EPSG2,ONEW
      REAL(KIND=LDP) XS(2),RS(2),XSN(2)
      REAL(KIND=LDP) T1,T2
      LOGICAL LFALL
      CHARACTER(LEN=32) CNAME
!
      CALL GETARG(1,CNAME)
      NCASE=0
      IF(CNAME .EQ. 'sc')  NCASE=1
      IF(CNAME .EQ. 'so')  NCASE=2
      IF(CNAME .EQ. 'sd')  NCASE=3
      IF(CNAME .EQ. 'sc0') NCASE=4
      IF(CNAME .EQ. 'scl') NCASE=5
      IF(CNAME .EQ. 'vec') NCASE=6
      IF(CNAME .EQ. 'vfp') NCASE=7
      IF(CNAME .EQ. 'vgen')NCASE=8
      IF(NCASE .EQ. 0)THEN
        WRITE(*,'(A)')'RESULT_FAIL unknown_case'
        STOP
      END IF
!
      OMIN2=KN_OSGN_MIN
      OMAX2=KN_OSGN_MAX
      EPSG2=KN_OSGN_EPSG
      RT=1.0_LDP
!
      IF(NCASE .EQ. 4)THEN
! Denominator collapse: dr numerically zero must hold omega.
        NRET=0;  OM=0.5_LDP;  ONEW=-9.0_LDP
        CALL KN_OSGN_STEP(0.0_LDP,0.0_LDP,1.0_LDP,2.0_LDP,OM, &
      &OMIN2,OMAX2,EPSG2,NRET,ONEW,ICODE)
        WRITE(*,'(A,I3)')'OSGN_ICODE=',ICODE
        WRITE(*,'(A,ES16.9)')'OSGN_OMEGA=',ONEW
        WRITE(*,'(A)')'RESULT_OK'
        STOP
      END IF
!
      IF(NCASE .EQ. 5)THEN
! Clamp behaviour, one kernel call each way.
        NRET=0;  OM=0.5_LDP
        R=1.0_LDP
        LAM=0.55_LDP
        DX=OM*R
        DR=OM*(LAM-1.0_LDP)*R
        RN=R+DR
        CALL KN_OSGN_STEP(DX*DR,DR*DR,RN*RN,R*R,OM, &
      &OMIN2,OMAX2,EPSG2,NRET,ONEW,ICODE)
        WRITE(*,'(A,I3)')'OSGN_ICODE_HI=',ICODE
        WRITE(*,'(A,ES16.9)')'OSGN_OMEGA_HI=',ONEW
        NRET=0;  OM=0.5_LDP
        R=1.0_LDP
        LAM=-1.5_LDP
        DX=OM*R
        DR=OM*(LAM-1.0_LDP)*R
        RN=R+DR
        CALL KN_OSGN_STEP(DX*DR,DR*DR,RN*RN,R*R,OM, &
      &0.5_LDP,OMAX2,EPSG2,NRET,ONEW,ICODE)
        WRITE(*,'(A,I3)')'OSGN_ICODE_LO=',ICODE
        WRITE(*,'(A,ES16.9)')'OSGN_OMEGA_LO=',ONEW
        WRITE(*,'(A)')'RESULT_OK'
        STOP
      END IF
!
      IF(NCASE .LE. 3)THEN
! Scalar lagged-map run: apply the returned omega each step.
        IF(NCASE .EQ. 1)LAM=0.3_LDP
        IF(NCASE .EQ. 2)LAM=-0.8_LDP
        IF(NCASE .EQ. 3)LAM=1.5_LDP
        X=2.0_LDP
        OM=0.5_LDP
        NRET=0
        NSTEP=0
        R=(LAM-1.0_LDP)*(X-RT)
        DO I=1,60
          X=X+OM*R
          RN=(LAM-1.0_LDP)*(X-RT)
          DX=OM*R
          DR=RN-R
          CALL KN_OSGN_STEP(DX*DR,DR*DR,RN*RN,R*R,OM, &
      &OMIN2,OMAX2,EPSG2,NRET,ONEW,ICODE)
          WRITE(*,'(A,I3,A,I3,A,ES16.9,A,ES16.9)') &
      &'OSGN_STEP ',I,' ICODE=',ICODE,' OMEGA=',ONEW, &
      &' R=',RN
          NSTEP=I
          IF(ICODE .EQ. 3)EXIT
          IF(ICODE .EQ. 1 .OR. ICODE .EQ. 2)OM=ONEW
          R=RN
          IF(ABS(X-RT) .LT. 1.0E-12_LDP)EXIT
        END DO
        WRITE(*,'(A,I4)')'OSGN_NSTEP=',NSTEP
        WRITE(*,'(A,ES16.9)')'OSGN_RESID=',ABS(X-RT)
        WRITE(*,'(A,ES16.9)')'OSGN_OMEGA=',OM
        WRITE(*,'(A,I3)')'OSGN_ICODE=',ICODE
        WRITE(*,'(A)')'RESULT_OK'
        STOP
      END IF
!
! Vector cases: two slots, one column of three frequencies is
! enough -- emulate the sweep loop the write hooks implement.
!
      CALL FG_OSGN_ALLOC(2,2,3)
      IF(.NOT. ALLOCATED(FG_OSGN_XP))THEN
        WRITE(*,'(A)')'RESULT_FAIL alloc'
        STOP
      END IF
      KN_ND_EXT=2
      KN_NCF=3
      IF(ALLOCATED(KN_CHI_EXT))DEALLOCATE(KN_CHI_EXT)
      IF(ALLOCATED(KN_ETA_EXT))DEALLOCATE(KN_ETA_EXT)
      ALLOCATE(KN_CHI_EXT(2,3),KN_ETA_EXT(2,3))
      KN_CHI_EXT=1.0
      KN_ETA_EXT=0.5
!
      LAM=0.3_LDP
      XS=(/2.0_LDP,3.0_LDP/)
      RS=(/0.0_LDP,0.0_LDP/)
      OM=0.5_LDP
!
      DO IS=1,3
! Sweep head: the column wrap in the sweep write bumps the
! generation and clears the dots.
        FG_OSGN_GEN=FG_OSGN_GEN+1
        FG_OSGN_ACC=0.0_LDP
        DO IC=1,3
          IF(FG_OSGN_CF .GE. 1)CALL FG_OSGN_COL_FIN(FG_OSGN_CF)
          FG_OSGN_CF=IC
! Slot 1 writes its P sample this column.
          T2=XS(1)
          T1=RT+LAM*(T2-RT)
          FG_OSGN_C1P(1,1)=T1
          FG_OSGN_C2P(1,1)=T2
          FG_OSGN_DP(1,1)=T2
          FG_OSGN_TP(1,1)=.TRUE.
          XS(1)=T2+OM*(T1-T2)
          RS(1)=T1-T2
! Slot 2 writes its M sample this column.
          T2=XS(2)
          T1=2.0_LDP*RT+LAM*(T2-2.0_LDP*RT)
          FG_OSGN_C1M(2,2)=T1
          FG_OSGN_C2M(2,2)=T2
          FG_OSGN_DM(2,2)=T2
          FG_OSGN_TM(2,2)=.TRUE.
          XS(2)=T2+OM*(T1-T2)
        END DO
        IF(NCASE .EQ. 7 .AND. IS .EQ. 3)KN_CHI_EXT(1,1)=2.0
        CALL FG_OSGN_UPDATE(OM,ONEW,LFALL,ICODE)
        WRITE(*,'(A,I3,A,I3,A,ES16.9,A,L2)') &
      &'OSGN_SWEEP ',IS,' ICODE=',ICODE,' OMEGA=',ONEW, &
      &' LFALL=',LFALL
        IF(ICODE .EQ. 1 .OR. ICODE .EQ. 2)OM=ONEW
      END DO
!
      WRITE(*,'(A,ES16.9)')'OSGN_OMEGA=',OM
      WRITE(*,'(A,ES16.9)')'OSGN_RESID=', &
      &MAX(ABS(XS(1)-RT),ABS(XS(2)-2.0_LDP*RT))
      WRITE(*,'(A,I3)')'OSGN_ICODE=',ICODE
      WRITE(*,'(A,I3)')'OSGN_NRET=',FG_OSGN_NRET
      WRITE(*,'(A,ES16.9)')'OSGN_RC2=',FG_OSGN_W(2)
      WRITE(*,'(A,ES16.9)')'OSGN_RP2=',FG_OSGN_W(3)
!
      IF(NCASE .EQ. 8)THEN
! Second update inside the same generation must fall back.
        CALL FG_OSGN_UPDATE(OM,ONEW,LFALL,ICODE)
        WRITE(*,'(A,I3)')'OSGN_ICODE2=',ICODE
        WRITE(*,'(A,L2)')'OSGN_LFALL2=',LFALL
        WRITE(*,'(A,ES16.9)')'OSGN_OMEGA2=',ONEW
      END IF
!
      DEALLOCATE(KN_CHI_EXT,KN_ETA_EXT)
      WRITE(*,'(A)')'RESULT_OK'
      STOP
!
      END
