!
! Driver for the single-precision monotone interpolators
! (SP_MON_INTERP in pgplt/do_vec_op.f, MON_INTERP_SP in
! pgplt/mon_interp_sp.f; both altered 05-Oct-2026 with the same
! grid-validation guards).  argv(1) selects the case; argv(2)
! selects the routine (1=SP_MON_INTERP, 2=MON_INTERP_SP, default 1).
! Bad grids must STOP inside the routine; the valid case prints
! 'SPM <i> <value>' lines for Python checks.
!
!	1  valid monotone R, ordered QZR        -> completes
!	2  repeated R endpoints                 -> STOP
!	3  reversed R                           -> STOP
!	4  unordered QZR                        -> STOP
!	5  QZR outside R range                  -> STOP
!	6  valid on a decreasing R grid         -> completes
!
	PROGRAM DRV_SP_MON
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: ND=10, NX=7, LIN=3
	REAL*4 :: QZ(NX,LIN), QZR(NX), VARRAY(ND,LIN), R(ND)
	CHARACTER*32 :: ARG
	INTEGER :: IC, WH, I, J
!
	IC=0
	WH=1
	CALL GET_COMMAND_ARGUMENT(1,ARG)
	READ(ARG,*) IC
	CALL GET_COMMAND_ARGUMENT(2,ARG)
	IF(ARG .NE. ' ')READ(ARG,*) WH
!
	DO I=1,ND
	  R(I)=1.0E0+0.5E0*(I-1)          !1.0 .. 5.5 increasing
	END DO
	DO J=1,LIN
	  DO I=1,ND
	    VARRAY(I,J)=R(I)**J          !J=1 linear (exact), J=2,3 smooth
	  END DO
	END DO
!
	SELECT CASE(IC)
	CASE(1)
	  DO I=1,NX
	    QZR(I)=1.2E0+0.6E0*(I-1)      !1.2 .. 4.8 inside [1,5.5]
	  END DO
	CASE(2)
	  R(ND)=R(ND-1)                   !repeated endpoint
	  DO I=1,NX
	    QZR(I)=1.2E0+0.6E0*(I-1)
	  END DO
	CASE(3)
	  R=R(ND:1:-1)                    !reversed grid
	  DO I=1,NX
	    QZR(I)=1.2E0+0.6E0*(I-1)
	  END DO
	CASE(4)
	  QZR(1)=1.2E0
	  QZR(2)=2.0E0
	  QZR(3)=1.9E0                     !unordered output
	  QZR(4)=2.5E0
	  QZR(5)=3.0E0
	  QZR(6)=3.5E0
	  QZR(7)=4.0E0
	CASE(5)
	  QZR(7)=6.0E0                    !outside R range
	  DO I=1,NX-1
	    QZR(I)=1.2E0+0.6E0*(I-1)
	  END DO
	CASE(6)
	  R=R(ND:1:-1)                    !decreasing grid is valid
	  VARRAY=VARRAY(ND:1:-1,:)        !values follow the same order
	  DO I=1,NX
	    QZR(I)=4.8E0-0.6E0*(I-1)      !ordered decreasing inside range
	  END DO
	CASE DEFAULT
	  STOP 'bad case'
	END SELECT
!
	IF(WH .EQ. 2)THEN
	  CALL MON_INTERP_SP(QZ,NX,LIN,QZR,NX,VARRAY,ND,R,ND)
	ELSE
	  CALL SP_MON_INTERP(QZ,NX,LIN,QZR,NX,VARRAY,ND,R,ND)
	END IF
!
	DO J=1,LIN
	  DO I=1,NX
	    WRITE(*,'(A,I3,I3,ES14.7)')'SPM ',I,J,QZ(I,J)
	  END DO
	END DO
!
	END PROGRAM DRV_SP_MON
