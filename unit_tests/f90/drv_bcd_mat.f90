! Driver for READ_BCD_MAT / WRITE_BCD_MAT (new_main/subs/read_bcd_mat.f).
!
! Alterations under test:
!   - READ_BCD_MAT now STOPs on an unrecognized TYPE instead of
!     printing the error and returning with output arrays untouched.
!   - WRITE_BCD_MAT reports its own name in the unrecognized-TYPE
!     error (upstream printed 'Error in READ_BCD_MAT') and also
!     STOPs rather than falling through to the IOSTAT check.
!
! Usage: drv CASE TYP
!   rt   : write TYPE=TYP scratch file then read it back in the same
!          process; every output array must return bitwise identical.
!   rd   : READ_BCD_MAT only (missing file or bad TYPE paths).
!   wr   : WRITE_BCD_MAT only (bad TYPE path).
!
! The routines derive the scratch filename as
!   TRIM(SCRATCH_DIR)//TYP//'SCRATCH'//I4.4
! so SCRATCH_DIR (env) must point at the case directory with a
! trailing '/'.
!
      PROGRAM DRV_BCD_MAT
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: N=4, NION=3, DEPTH_INDX=7
      REAL(KIND=LDP) B(N,N),C(N,N),D(N,N)
      REAL(KIND=LDP) ROW_SF(N),COL_SF(N),POPS(N)
      REAL(KIND=LDP) B0(N,N),C0(N,N),D0(N,N)
      REAL(KIND=LDP) ROW_SF0(N),COL_SF0(N),POPS0(N)
      INTEGER IPIVOT(N),IPIVOT0(N)
      LOGICAL REPLACE_EQ(NION),ZERO_STEQ(N)
      LOGICAL REPLACE_EQ0(NION),ZERO_STEQ0(N)
      CHARACTER(LEN=32) CASE
      CHARACTER(LEN=8) TYP
      INTEGER I,J
      LOGICAL OK
      CALL GETARG(1,CASE)
      CALL GETARG(2,TYP)
      IF(CASE .EQ. 'wr')THEN
        CALL FILL(B,C,D,ROW_SF,COL_SF,POPS,IPIVOT,REPLACE_EQ,ZERO_STEQ,&
     &            N,NION)
        CALL WRITE_BCD_MAT(B,C,D,ROW_SF,COL_SF,IPIVOT,&
     &            POPS,REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,TYP)
        WRITE(*,'(A)')'WROTE'
      ELSE IF(CASE .EQ. 'rd')THEN
        B=0.0_LDP; C=0.0_LDP; D=0.0_LDP
        ROW_SF=0.0_LDP; COL_SF=0.0_LDP; POPS=0.0_LDP
        IPIVOT=0; REPLACE_EQ=.FALSE.; ZERO_STEQ=.FALSE.
        CALL READ_BCD_MAT(B,C,D,ROW_SF,COL_SF,IPIVOT,&
     &            POPS,REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,TYP)
        WRITE(*,'(A)')'READ'
      ELSE IF(CASE .EQ. 'rt')THEN
        CALL FILL(B,C,D,ROW_SF,COL_SF,POPS,IPIVOT,REPLACE_EQ,ZERO_STEQ,&
     &            N,NION)
        B0=B; C0=C; D0=D
        ROW_SF0=ROW_SF; COL_SF0=COL_SF; POPS0=POPS
        IPIVOT0=IPIVOT; REPLACE_EQ0=REPLACE_EQ; ZERO_STEQ0=ZERO_STEQ
        CALL WRITE_BCD_MAT(B,C,D,ROW_SF,COL_SF,IPIVOT,&
     &            POPS,REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,TYP)
        B=0.0_LDP; C=0.0_LDP; D=0.0_LDP
        ROW_SF=0.0_LDP; COL_SF=0.0_LDP; POPS=0.0_LDP
        IPIVOT=0; REPLACE_EQ=.FALSE.; ZERO_STEQ=.FALSE.
        CALL READ_BCD_MAT(B,C,D,ROW_SF,COL_SF,IPIVOT,&
     &            POPS,REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,TYP)
        OK=.TRUE.
        IF(TYP .EQ. 'BC')OK=OK .AND. ALL(B .EQ. B0)
        IF(TYP .NE. 'D')THEN
          OK=OK .AND. ALL(C .EQ. C0)
          OK=OK .AND. ALL(ROW_SF .EQ. ROW_SF0)
          OK=OK .AND. ALL(COL_SF .EQ. COL_SF0)
          OK=OK .AND. ALL(POPS .EQ. POPS0)
          OK=OK .AND. ALL(IPIVOT .EQ. IPIVOT0)
          OK=OK .AND. ALL(REPLACE_EQ .EQV. REPLACE_EQ0)
          OK=OK .AND. ALL(ZERO_STEQ .EQV. ZERO_STEQ0)
        END IF
        IF(TYP .EQ. 'D')OK=OK .AND. ALL(D .EQ. D0)
        IF(OK)THEN
          WRITE(*,'(A)')'RT_OK'
        ELSE
          WRITE(*,'(A)')'RT_MISMATCH'
        END IF
      ELSE
        WRITE(*,'(A)')'BAD_CASE'
      END IF
      END
!
      SUBROUTINE FILL(B,C,D,ROW_SF,COL_SF,POPS,IPIVOT,REPLACE_EQ,&
     &                ZERO_STEQ,N,NION)
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER N,NION
      REAL(KIND=LDP) B(N,N),C(N,N),D(N,N)
      REAL(KIND=LDP) ROW_SF(N),COL_SF(N),POPS(N)
      INTEGER IPIVOT(N)
      LOGICAL REPLACE_EQ(NION),ZERO_STEQ(N)
      INTEGER I,J
      DO J=1,N
        DO I=1,N
          B(I,J)=I*10.0_LDP+J*0.25_LDP-3.0_LDP
          C(I,J)=-I*1.5_LDP+J*7.0_LDP
          D(I,J)=(I-J)*0.125_LDP+1.0_LDP
        END DO
        ROW_SF(I)=1.0_LDP+0.5_LDP*I
        COL_SF(I)=2.0_LDP-0.1_LDP*I
        POPS(I)=I*1.0E+10_LDP-2.5_LDP
        IPIVOT(I)=N+1-I
        ZERO_STEQ(I)=(MOD(I,2) .EQ. 0)
      END DO
      DO I=1,NION
        REPLACE_EQ(I)=(MOD(I,2) .EQ. 1)
      END DO
      END
