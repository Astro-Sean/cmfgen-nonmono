C
C Subroutine to perform liner interpolations for RDPOP.
C The edges of the new mesh must lie on or inside the edges of
C the old mesh. R may be monotonically increasing, or decreasing.
C
C Improved error handling (5-OCT-1983)
C
C Altered 05-Oct-2026 - Validate the interpolation coordinates. The
C                       bracket walk requires strictly monotonic R and U
C                       in the same direction (Devin).
C Altered 11-Apr-1989 - Implicit NONE INSTALLED. Main call changed to
C                       LIN_INTERP, and LINPOP and LININT entry points
C                       installed. All three routines are identical
C                       and were all included so save editing a lot of
C                       routines that call LINPOP and LINIT. LIN_INTERP
C                       is the preferred call.
C
	SUBROUTINE LIN_INTERP(R,V,ND,U,W,NIN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER NIN,ND
	REAL(KIND=LDP) R(ND),V(ND),U(NIN),W(NIN)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local varaibles.
C
	REAL(KIND=LDP) T1
	INTEGER J,L,LUER
	LOGICAL INCREASING
C
	ENTRY LININT(R,V,ND,U,W,NIN)
	ENTRY LINPOP(R,V,ND,U,W,NIN)
C
	IF(ND .LT. 2 .OR. NIN .LT. 2)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in LIN_INTERP/LINPOP'
	  WRITE(LUER,*)'At least two grid points are required'
	  STOP
	END IF
	INCREASING=(R(2)-R(1) .GT. 0.0_LDP)
	DO J=1,ND-1
	  IF(INCREASING)THEN
	    IF(R(J+1) .LE. R(J))THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in LIN_INTERP/LINPOP'
	      WRITE(LUER,*)'R values must be strictly increasing'
	      STOP
	    END IF
	  ELSE IF(R(J+1) .GE. R(J))THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in LIN_INTERP/LINPOP'
	    WRITE(LUER,*)'R values must be strictly decreasing'
	    STOP
	  END IF
	END DO
	DO J=1,NIN-1
	  IF(INCREASING)THEN
	    IF(U(J+1) .LE. U(J))THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in LIN_INTERP/LINPOP'
	      WRITE(LUER,*)'U values must be strictly increasing'
	      STOP
	    END IF
	  ELSE IF(U(J+1) .GE. U(J))THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in LIN_INTERP/LINPOP'
	    WRITE(LUER,*)'U values must be strictly decreasing'
	    STOP
	  END IF
	END DO
C
	IF(INCREASING)THEN
	  IF(R(ND) .GT. U(NIN) .OR. R(1) .LT. U(1))THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in LIN_INTERP/LINPOP - values outside range'
	    WRITE(LUER,*)R(ND),'should be LE',U(NIN)
	    WRITE(LUER,*)R(1),'should be GE ',U(1)
	    WRITE(LUER,*)'ND=',ND,'NIN=',NIN
	    STOP
	  END IF
	  L=2
	  DO J=1,ND
50	    IF(R(J) .LE. U(L))THEN
	      T1=(U(L)-R(J))/(U(L)-U(L-1))
	      V(J)=T1*W(L-1)+(1.0_LDP-T1)*W(L)
	    ELSE
	      L=L+1
	      GOTO 50
	    END IF
	  END DO
	ELSE
C
C**********************************************************************
C**********************************************************************
C
	  IF(R(1) .GT. U(1) .OR. R(ND) .LT. U(NIN))THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in LIN_INTERP/LINPOP - values outside range'
	    WRITE(LUER,*)R(ND),'should be GE',U(NIN)
	    WRITE(LUER,*)R(1),'should be LE ',U(1)
	    WRITE(LUER,*)'ND=',ND,'NIN=',NIN
	    STOP
	  END IF
	  L=NIN-1		!-1 to stay in array bounds
	  DO J=ND,1,-1
500	    IF(R(J) .LE. U(L))THEN
	      T1=(U(L+1)-R(J))/(U(L+1)-U(L))
	      V(J)=T1*W(L)+(1.0_LDP-T1)*W(L+1)
	    ELSE
	      L=L-1
	      GOTO 500
	    END IF
	  END DO
C
	END IF
C
	RETURN
	END
