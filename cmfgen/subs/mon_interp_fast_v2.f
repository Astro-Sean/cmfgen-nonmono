!
! Subroutine to interpolate an array onto a new grid. The grid vector must be
! either a monotonically decreasing or increasing function. A modified cubic
! polynomial is used to do the interpolation. Instead of using
! the exact cubic estimates for the first derivative at the two nodes,
! we use revised estimates which insure that the interpolating function
! is monotonic in the interpolating interval.
!
! Disadvantages: The interpolating weights can only be defined when the
!                function is known. In principal could use these modified
!                first derivatives to compute an accurate integration
!                formulae. However, the integration weights cannot be defined
!                independently of the function values, as desired in many
!                situations.
!
! Ref: Steffen. M, 1990, A/&A, 239, 443-450
!
	SUBROUTINE MON_INTERP_FAST_V2(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND,IDENT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 06-Oct-2014 : Changed to V2. Added IDENT keyword.
!   Moved 03-Jan-2013 : Moved from obs to subs directory.
! Created 09-Dec-1998 : Based on MON_INTERP.
!                       Design to be fast when from creating a large
!                           array from a much smaller array.
!
	INTEGER NQ,LIN_END,NX,NV,ND
	REAL(KIND=LDP) QZ(NQ,LIN_END),QZR(NX)
	REAL(KIND=LDP) VARRAY(NV,LIN_END),R(ND)
!
	REAL(KIND=LDP) S(ND)		!Slopes
	REAL(KIND=LDP) H(ND)
	INTEGER LST_INTERVAL
	INTEGER ND_SM
	INTEGER IVEC(NX)
	CHARACTER(LEN=*) IDENT
!
	REAL(KIND=LDP) ONE
	PARAMETER (ONE=1.0_LDP)
	INTEGER I,J,ML
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) A(ND)
	REAL(KIND=LDP) B(ND)
	REAL(KIND=LDP) C(ND)
	REAL(KIND=LDP) D(ND)		!Used for derivative at I.
	REAL(KIND=LDP) E(ND)
	REAL(KIND=LDP) SGN
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! The array R may be either monotonically increasing, or decreasing.
! Reject invalid coordinates before the interval walk; a reversed query can
! otherwise advance I beyond the source grid, and repeated R values divide by zero.
! This unconditional guard leaves valid monotone inputs unchanged and stops bad grids.
!
	LUER=ERROR_LU()
	IF(ND .LT. 3 .OR. NX .LT. 1)THEN
	  WRITE(LUER,*)'Error in MON_INTERP_FAST_V2 - invalid grid lengths'
	  WRITE(LUER,*)'ND,NX=',ND,NX
	  WRITE(LUER,*)'IDENT='//TRIM(IDENT)
	  STOP
	END IF
	SGN=SIGN(ONE,R(ND)-R(1))
	IF(SGN .EQ. 0.0_LDP)THEN
	  WRITE(LUER,*)'Error in MON_INTERP_FAST_V2 - repeated input endpoints'
	  WRITE(LUER,*)'IDENT='//TRIM(IDENT)
	  STOP
	END IF
	IF( (SGN*QZR(1) .LT. SGN*R(1)) .OR.
	1   (SGN*QZR(NX) .GT. SGN*R(ND)) )THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in MON_INTERP_FAST_V2 - values outside range'
	  WRITE(LUER,*)'IDENT='//TRIM(IDENT)
	  WRITE(LUER,*)'Ranges shown below'
	  WRITE(LUER,*)'R:',R(1),R(ND)
	  WRITE(LUER,*)'QZR:',QZR(1),QZR(NX)
	  STOP
	END IF
!
! Determine intervals and slopes to minimize computational effort.
!
	DO I=1,ND-1
	  H(I)=R(I+1)-R(I)
	END DO
!
! Check that R and QZR vectors are monotonic.
!
	DO I=1,ND-1
	  IF(SGN*H(I) .LE. 0.0_LDP)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in MON_INTERP_FAST_V2: IDENT='//TRIM(IDENT)
	    WRITE(LUER,*)'R values must be strictly monotonic'
	    WRITE(LUER,*)'I,R(I),R(I+1)=',I,R(I),R(I+1)
	    STOP
	  END IF
	END DO
	DO I=1,NX-1
	  IF(SGN*(QZR(I+1)-QZR(I)) .LE. 0.0_LDP)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in MON_INTERP_FAST_V2: IDENT='//TRIM(IDENT)
	    WRITE(LUER,*)'QZR values must be strictly ordered with R'
	    WRITE(LUER,*)'I,QZR(I),QZR(I+1)=',I,QZR(I),QZR(I+1)
	    STOP
	  END IF
	END DO
!
! Determine the interval (R(I) to R(I+1)) containing QZR(J).
!
	I=1
	DO J=1,NX
500	  IF( SGN*QZR(J) .LE. SGN*R(I+1))THEN
	    IVEC(J)=I
	  ELSE
	    I=I+1
	    GOTO 500
	  END IF
	END DO
	LST_INTERVAL=MAXVAL(IVEC)
	ND_SM=LST_INTERVAL+1
!
! Loop over frequency space.
!
!$OMP PARALLEL PRIVATE(D,S,A,B,C,E,T1,I,J,ML)
!$OMP DO
	DO ML=1,LIN_END
!
! Compute the slopes.
!
	  DO I=1,MIN(ND-1,ND_SM)
	    S(I)=(VARRAY(I+1,ML)-VARRAY(I,ML))/H(I)
	  END DO
!
! Compute the first derivatives at node I.
!
          D(1)=S(1) +(S(1)-S(2))*H(1)/(H(1)+H(2))
	  DO I=2,MIN(ND-1,ND_SM)
            D(I)=(S(I-1)*H(I)+S(I)*H(I-1))/(H(I-1)+H(I))
	  END DO
!
! Adjust first derivatives so that function is monotonic  in each interval.
!
	  D(1)=( SIGN(ONE,S(1))+SIGN(ONE,D(1)) )*MIN(ABS(S(1)),0.5_LDP*ABS(D(1)))
	  DO I=2,MIN(ND-1,ND_SM)
	    D(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(D(I)))
	  END DO
!
! Special treatment for computing slope at end point of R.
!
	  IF(ND .EQ. ND_SM)THEN
            D(ND)=S(ND-1)+(S(ND-1)-S(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
	    D(ND)=( SIGN(ONE,S(ND-1))+SIGN(ONE,D(ND)) )*
	1      MIN(ABS(S(ND-1)),0.5_LDP*ABS(D(ND)))
	  END IF
!
! Determine the ciefficients of the monotonic cubic polynomial.
!
! If T1=X-R(I) then
!             Y=A(I)*T1^3 + B(I)*T1^3 + C(I)*T1 + E(I)
!
	  DO I=1,ND_SM-1
            A(I)=(D(I)+D(I+1)-2.0_LDP*S(I))/H(I)/H(I)
	    B(I)=(3.0_LDP*S(I)-2.0_LDP*D(I)-D(I+1))/H(I)
	    C(I)=D(I)
	    E(I)=VARRAY(I,ML)
	  END DO
!
! Perform the interpolations.
!
	  DO J=1,NX
	    I=IVEC(J)
	    T1=(QZR(J)-R(I))
            QZ(J,ML)=((A(I)*T1+B(I))*T1+C(I))*T1+E(I)
	  END DO
!
	END DO
!$OMP END DO
!$OMP END PARALLEL
!
	RETURN
	END
