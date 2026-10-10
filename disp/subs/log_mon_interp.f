!
! Subroutine to interpolate an array onto a new grid. The grid vector must be
! either a monotonically decreasing or increasing function. A modified cubic
! polynomial is used to do the interpolation. Instead of using
! the excact cubic estiamtes for the first derivative at the two nodes,
! we use revised estimates which insure that the interpolating function
! is mononotonic in the interpolating interval.
!
! The techniques is somewhat similar to that suggested by Nordulund.
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
	SUBROUTINE LOG_MON_INTERP(QZ,NQ,LIN_END,QZRIN,NX,VIN,NV,RIN,ND,LOGX,LOGY)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 05-Oct-2026 : Validate positive/monotonic coordinates before
!                         taking LOG(X) and interpolating (Devin).
! Altered 24-May-1996 : ERROR_LU installed
! Created 01-Apr-1992 : Code may need recoding for optimal speed, and for
!                         vectorization.
!
	INTEGER NQ,LIN_END,NX,NV,ND
!
	REAL(KIND=LDP) QZRIN(NX)
	REAL(KIND=LDP) VIN(NV,LIN_END),RIN(ND)
!
	REAL(KIND=LDP) QZ(NQ,LIN_END),QZR(NX)
	REAL(KIND=LDP) VARRAY(NV,LIN_END),R(ND)
!
	LOGICAL LOGX,LOGY
!
	REAL(KIND=LDP) ONE
	PARAMETER (ONE=1.0_LDP)
	INTEGER I,J,M
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) HI,HIM1,HIP1
	REAL(KIND=LDP) SI,SIM1,SIP1
	REAL(KIND=LDP) A,B,C,D,DYI,DYIP1,SGN
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	IF(LOGX)THEN
	  DO I=1,ND
	    IF(RIN(I) .LE. 0.0_LDP)THEN
	      WRITE(6,*)'Error in LOG_MON_INTERP - input coordinate not positive'
	      WRITE(6,*)'I,RIN(I)=',I,RIN(I)
	      STOP
	    END IF
	  END DO
	  DO I=1,NX
	    IF(QZRIN(I) .LE. 0.0_LDP)THEN
	      WRITE(6,*)'Error in LOG_MON_INTERP - output coordinate not positive'
	      WRITE(6,*)'I,QZRIN(I)=',I,QZRIN(I)
	      STOP
	    END IF
	  END DO
	  R(1:ND)=LOG(RIN(1:ND))
	  QZR(1:NX)=LOG(QZRIN(1:NX))
	ELSE
	  R(1:ND)=RIN(1:ND)
	  QZR(1:NX)=QZRIN(1:NX)
	END IF
	IF(LOGY)THEN
	  DO J=1,LIN_END
	    DO I=1,ND
	      IF(VIN(I,J) .GT. 0.0_LDP)THEN
	        VARRAY(I,J)=LOG(VIN(I,J))
	      ELSE IF(VIN(I,J) .EQ. 0.0_LDP)THEN
	        VARRAY(I,J)=1.0E-60_LDP
	      ELSE
	        VARRAY(I,J)=LOG(ABS(VIN(I,J)/10.0_LDP))
	      END IF
	    END DO
	  END DO
	ELSE
	  VARRAY=VIN
	END IF
!
! The array R may be either monotonically increasing, or decreasing.
! A repeated or reversed input coordinate has no unique inverse branch;
! QZR may repeat but must preserve the same ordering.
!
	IF(ND .LT. 3 .OR. NX .LT. 1)THEN
	  WRITE(6,*)'Error in LOG_MON_INTERP - invalid grid lengths'
	  WRITE(6,*)'ND,NX=',ND,NX
	  STOP
	END IF
	SGN=SIGN(ONE,R(ND)-R(1))
	IF(SGN .EQ. 0.0_LDP)THEN
	  WRITE(6,*)'Error in LOG_MON_INTERP - repeated input endpoints'
	  STOP
	END IF
	DO I=1,ND-1
	  IF(SGN*(R(I+1)-R(I)) .LE. 0.0_LDP)THEN
	    WRITE(6,*)'Error in LOG_MON_INTERP - input grid not monotonic'
	    WRITE(6,*)'I,R(I),R(I+1)=',I,RIN(I),RIN(I+1)
	    STOP
	  END IF
	END DO
	DO I=1,NX-1
	  IF(SGN*(QZR(I+1)-QZR(I)) .LT. 0.0_LDP)THEN
	    WRITE(6,*)'Error in LOG_MON_INTERP - output grid not ordered'
	    WRITE(6,*)'I,QZR(I),QZR(I+1)=',I,QZRIN(I),QZRIN(I+1)
	    STOP
	  END IF
	END DO
	IF( (SGN*QZR(1) .LT. SGN*R(1)) .OR.
	1   (SGN*QZR(NX) .GT. SGN*R(ND)) )THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in MON_INTERP - values outside range'
	  WRITE(LUER,'(1X,2(A,ES20.12),A,I4)')'  R(1)=',R(1),'    R(ND)= ',R(ND),'ND=',ND
	  WRITE(LUER,'(1X,2(A,ES20.12),A,I4)')'QZR(1)=',QZR(1),'  QZR(NX)=',QZR(NX),'NX=',NX
	  STOP
	END IF
!
	DO I=1,ND
	  WRITE(171,*)I,RIN(I),R(I)
	END DO
	DO I=1,NX
	  WRITE(171,*)I,QZRIN(I),QZR(I)
	END DO
!
! M is the Index in new interpolated array
!
	I=1
	DO M=1,NX
500	  IF( SGN*QZR(M) .LE. SGN*R(I+1))THEN
	    IF(I .EQ. 1)THEN
	      HI=R(2)-R(1)
              HIP1=R(3)-R(2)
              DO J=1,LIN_END
                SI=(VARRAY(2,J)-VARRAY(1,J))/HI
                SIP1=(VARRAY(3,J)-VARRAY(2,J))/HIP1
                DYI=SI +(SI-SIP1)*HI/(HI+HIP1)
                DYIP1=(SI*HIP1+SIP1*HI)/(HI+HIP1)
	        DYI=( SIGN(ONE,SI)+SIGN(ONE,DYI) )*
	1            MIN(ABS(SI),0.5_LDP*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,SIP1) )*
	1            MIN(ABS(SI),ABS(SIP1),0.5_LDP*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0_LDP*SI)/HI/HI
	        B=(3.0_LDP*SI-2.0_LDP*DYI-DYIP1)/HI
	        C=DYI
	        D=VARRAY(I,J)
                QZ(M,J)=((A*T1+B)*T1+C)*T1+D
	      END DO
	    ELSE IF(I .EQ. ND-1)THEN
	      HI=R(ND)-R(ND-1)
              HIM1=R(ND-1)-R(ND-2)
              DO J=1,LIN_END
                SIM1=(VARRAY(ND-1,J)-VARRAY(ND-2,J))/HIM1
                SI=(VARRAY(ND,J)-VARRAY(ND-1,J))/HI
                DYI=(SIM1*HI+SI*HIM1)/(HIM1+HI)
                DYIP1=SI+(SI-SIM1)*HI/(HIM1+HI)
	        DYI=( SIGN(ONE,SIM1)+SIGN(ONE,SI) )*
	1            MIN(ABS(SIM1),ABS(SI),0.5_LDP*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,DYIP1) )*
	1            MIN(ABS(SI),0.5_LDP*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0_LDP*SI)/HI/HI
	        B=(3.0_LDP*SI-2.0_LDP*DYI-DYIP1)/HI
	        C=DYI
	        D=VARRAY(I,J)
                QZ(M,J)=((A*T1+B)*T1+C)*T1+D
	      END DO
	    ELSE
	      IF(I .GT. ND-1)THEN
	        WRITE(6,*)R(ND-2:ND)
	        WRITE(6,*)QZR(NX-2:NX)
	      END IF
	      HI=R(I+1)-R(I)
              HIM1=R(I)-R(I-1)
              HIP1=R(I+2)-R(I+1)
              DO J=1,LIN_END
                SIM1=(VARRAY(I,J)-VARRAY(I-1,J))/HIM1
                SI=(VARRAY(I+1,J)-VARRAY(I,J))/HI
                SIP1=(VARRAY(I+2,J)-VARRAY(I+1,J))/HIP1
                WRITE(6,'(4ES16.6)')SIM1,HI,SI,HIM1
	        DYI=(SIM1*HI+SI*HIM1)/(HIM1+HI)
                DYIP1=(SI*HIP1+SIP1*HI)/(HI+HIP1)
	        DYI=( SIGN(ONE,SIM1)+SIGN(ONE,SI) )*
	1            MIN(ABS(SIM1),ABS(SI),0.5_LDP*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,SIP1) )*
	1            MIN(ABS(SI),ABS(SIP1),0.5_LDP*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0_LDP*SI)/HI/HI
	        B=(3.0_LDP*SI-2.0_LDP*DYI-DYIP1)/HI
	        C=DYI
	        D=VARRAY(I,J)
                QZ(M,J)=((A*T1+B)*T1+C)*T1+D
	      END DO
	    END IF
	  ELSE
	    I=I+1
	    GOTO 500
	  END IF
!
	END DO
	IF(LOGY)QZ=EXP(QZ)
!
	RETURN
	END
