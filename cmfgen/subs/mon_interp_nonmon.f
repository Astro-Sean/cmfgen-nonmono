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
	SUBROUTINE MON_INTERP_NONMON(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND,IDENT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created: 21-Jun-2021 : The array QZR does not need to be monotonic (R must
!                           still be momonotonic). Based on mon_interp_fast.f
! Moved:   03-Jan-2013 : Moved from obs to subs directory.
! Created: 09-Dec-1998 : Based on MON_INTERP.
!                        Design to be fast when from creating a large
!                           array from a much smaller array.
!
	INTEGER NQ,LIN_END,NX,NV,ND
	REAL(KIND=LDP) QZ(NQ,LIN_END),QZR(NX)
	REAL(KIND=LDP) VARRAY(NV,LIN_END),R(ND)
	CHARACTER(LEN=*) IDENT
!
	REAL(KIND=LDP) S(ND)		!Slopes
	REAL(KIND=LDP) H(ND)
	INTEGER LST_INTERVAL
	INTEGER ND_SM
	INTEGER IVEC(NX)
!
	REAL(KIND=LDP) ONE
	PARAMETER (ONE=1.0_LDP)
	INTEGER I,J,ML
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) RMIN,RMAX
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
!
	RMAX=MAXVAL(QZR)
	RMIN=MINVAL(QZR)
	SGN=SIGN(ONE,R(ND)-R(1))
	IF( (RMIN .LT. SGN*R(1)) .OR.
	1   (RMAX .GT. SGN*R(ND)) )THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in MON_INTERP_NOMON- values outside range'
	  WRITE(LUER,'(2(A,ES18.8))')'  New R(1)=',RMIN,'   Old R(1)=',R(1)
	  WRITE(LUER,'(2(A,ES18.8))')' New R(NX)=',RMAX,'  Old R(ND)=',R(ND)
	  WRITE(LUER,'(A)')' IDENT=',TRIM(IDENT)
	  STOP
	END IF
!
! Determine intervals and slopes to minimize computational effort.
!
	DO I=1,ND-1
	  H(I)=R(I+1)-R(I)
	END DO
!
! Check that R vector is monotonic.
!
	DO I=1,ND-2
	  IF(H(I)*H(I+1) .LE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in MON_INTERP_NOMON'
	    WRITE(LUER,*)'R values must be monotonic'
	    WRITE(LUER,*)'I value is',I
	    WRITE(LUER,'(A)')' IDENT=',TRIM(IDENT)
	    H(ND)=0.0_LDP
	    WRITE(LUER,*)'J, R, dR follows'
	    DO J=1,ND
	      WRITE(LUER,*)J,R(J),H(J)
	    END DO
	    STOP
	  END IF
	END DO
!
! Determine the interval (R(I) to R(I+1)) containing QZR(J).
!
	I=1
	DO J=1,NX
500	  IF( SGN*QZR(J) .LT. SGN*R(I))THEN
	    I=I-1
	    GOTO 500
	  ELSE IF( SGN*QZR(J) .GT. SGN*R(I+1))THEN
	    I=I+1
	    GOTO 500
	  ELSE
	    IVEC(J)=I
	  END IF
	END DO
	IF(MINVAL(IVEC) .LT. 1 .OR. MAXVAL(IVEC) .GT. ND)THEN
	  WRITE(6,*)'Error setting pinter in MON_INTERP_NONMON'
	  WRITE(6,*)MINVAL(IVEC),MAXVAL(IVEC),ND
	END IF
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
