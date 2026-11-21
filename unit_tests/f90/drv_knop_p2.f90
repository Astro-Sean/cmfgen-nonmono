! Driver for the KN_P2_* state primitives (phase 2.29A).
! Exercises the shared tracker/clear/midpoint routines that the
! three ordered-march call sites use, printing KEY=value lines for
! the python harness. P2N=2, P2REL=3 throughout.
!
      PROGRAM DRV_KNOP_P2
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1, ONLY: KN_P2_TRACK,KN_P2_CLEAR, &
          KN_P2_MID
      IMPLICIT NONE
      INTEGER*1 SG,FL,DM
      INTEGER P2N,P2RL,K
      REAL(KIND=LDP) VO,VN,X,XLIM
      LOGICAL FIN
!
      P2N=2; P2RL=3
!
! R1/R8: genuine alternation activates damping on the P2N-th
! consecutive sign reversal.
!
      SG=0_1; FL=0_1; DM=0_1
      VO=1.0_LDP; VN=1.1_LDP
      CALL KN_P2_TRACK(SG,FL,DM,VO,VN,P2N,P2RL)   ! + : fresh sign
      PRINT '(A,3(I3),L2)','R1a ',INT(SG),INT(FL),INT(DM), &
        DM .EQ. 0_1 .AND. SG .EQ. 1_1 .AND. FL .EQ. 0_1
      VO=1.1_LDP; VN=1.0_LDP
      CALL KN_P2_TRACK(SG,FL,DM,VO,VN,P2N,P2RL)   ! - : flip 1
      PRINT '(A,3(I3),L2)','R1b ',INT(SG),INT(FL),INT(DM), &
        FL .EQ. 1_1 .AND. DM .EQ. 0_1
      VO=1.0_LDP; VN=1.05_LDP
      CALL KN_P2_TRACK(SG,FL,DM,VO,VN,P2N,P2RL)   ! + : flip 2 -> arm
      PRINT '(A,3(I3),L2)','R1c ',INT(SG),INT(FL),INT(DM), &
        FL .EQ. 2_1 .AND. DM .EQ. 1_1
!
! R9: a reversal while damped re-arms (DM=1), it does not release.
!
      VO=1.05_LDP; VN=1.0_LDP
      CALL KN_P2_TRACK(SG,FL,DM,VO,VN,P2N,P2RL)   ! - : flip, re-arm
      PRINT '(A,3(I3),L2)','R9a ',INT(SG),INT(FL),INT(DM), &
        DM .EQ. 1_1 .AND. FL .EQ. 3_1
!
! R2/R9: P2REL flip-free same-sign passes release the damping.
!
      DO K=1,3
        VO=VN; VN=VN-0.01_LDP
        CALL KN_P2_TRACK(SG,FL,DM,VO,VN,P2N,P2RL)
      END DO
      PRINT '(A,3(I3),L2)','R2  ',INT(SG),INT(FL),INT(DM), &
        DM .EQ. 0_1 .AND. FL .EQ. 0_1
!
! R3/R5: KN_P2_CLEAR on the insignificant branch resets all three
! state slots; the next significant update starts a clean history.
!
      SG=1_1; FL=2_1; DM=1_1
      CALL KN_P2_CLEAR(SG,FL,DM)
      PRINT '(A,3(I3),L2)','R3  ',INT(SG),INT(FL),INT(DM), &
        SG .EQ. 0_1 .AND. FL .EQ. 0_1 .AND. DM .EQ. 0_1
      VO=1.0_LDP; VN=0.9_LDP
      CALL KN_P2_TRACK(SG,FL,DM,VO,VN,P2N,P2RL)   ! - : fresh sign
      PRINT '(A,3(I3),L2)','R5  ',INT(SG),INT(FL),INT(DM), &
        SG .EQ. -1_1 .AND. FL .EQ. 0_1 .AND. DM .EQ. 0_1
!
! Zero-delta update is a non-reversal: FLIP resets and release
! counting advances while damped.
!
      SG=1_1; FL=1_1; DM=1_1
      VO=1.0_LDP; VN=1.0_LDP
      CALL KN_P2_TRACK(SG,FL,DM,VO,VN,P2N,P2RL)
      PRINT '(A,3(I3),L2)','RZ  ',INT(SG),INT(FL),INT(DM), &
        SG .EQ. 0_1 .AND. FL .EQ. 0_1 .AND. DM .EQ. 2_1
!
! R14: midpoint finiteness across the supported positive range.
!
      FIN=.TRUE.
      VO=1.0E-283_LDP; VN=1.1E-283_LDP
      X=KN_P2_MID(VO,VN)
      XLIM=SQRT(VO)*SQRT(VN)
      PRINT '(A,2(ES12.5),L2)','M1  ',X,XLIM, &
        X .GT. 0.0_LDP .AND. ABS(X/XLIM-1.0_LDP) .LT. 1.0E-12_LDP
      VO=1.0E+200_LDP; VN=1.2E+200_LDP
      X=KN_P2_MID(VO,VN)
      FIN=ABS(X/(SQRT(VO)*SQRT(VN))-1.0_LDP) .LT. 1.0E-12_LDP
      PRINT '(A,ES12.5,L2)','M2  ',X,FIN
!     product form underflows to zero here; sqrt form must not
      VO=1.0E-200_LDP; VN=4.0E-200_LDP
      X=KN_P2_MID(VO,VN)
      PRINT '(A,ES12.5,L2)','M3  ',X,ABS(X-2.0E-200_LDP) .LT. &
        1.0E-212_LDP .OR. ABS(X/2.0E-200_LDP-1.0_LDP) .LT. 1.0E-6_LDP
!     zero endpoint -> arithmetic fallback
      VO=0.0_LDP; VN=1.0_LDP
      X=KN_P2_MID(VO,VN)
      PRINT '(A,ES12.5,L2)','M4  ',X,ABS(X-0.5_LDP) .LT. 1.0E-12_LDP
!     non-positive guard: negative endpoint takes the same branch
      VO=-1.0_LDP; VN=1.0_LDP
      X=KN_P2_MID(VO,VN)
      PRINT '(A,ES12.5,L2)','M5  ',X,ABS(X) .LT. 1.0E-12_LDP
      STOP
      END PROGRAM DRV_KNOP_P2
