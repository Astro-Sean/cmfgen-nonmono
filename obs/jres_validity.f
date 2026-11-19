! jres_validity.f
!
! Elementwise field-validity helpers for the J_RESCUE machinery in
! COMP_JCONT_V4.INC.  The previous checks used
! T=MAXVAL(ABS(field)) followed by T.NE.T to catch NaNs; a MAXVAL
! reduction is not required to propagate NaN on every compiler
! lane order, so a non-finite element could be masked by a larger
! finite one.  These routines test every element with
! IEEE_IS_FINITE and classify the failure so the caller can report
! the trigger, not just that something failed.
!
! JRES_TRIG returns:
!   0 - field valid (all elements finite, max|a| <= CAP)
!   1 - NaN present (BAD = first non-finite index)
!   2 - Inf present (BAD = first non-finite index)
!   3 - finite but max|a| > CAP (BAD = argmax, TMAX = max|a|)
! A non-finite element takes precedence over a magnitude
! violation: the argmax of a poisoned field is meaningless.
!
! edited by Sean Brennan
!
	INTEGER FUNCTION JRES_TRIG(A,N,CAP,BAD,TMAX)
	USE SET_KIND_MODULE
	USE IEEE_ARITHMETIC
	IMPLICIT NONE
	INTEGER N,BAD,I
	REAL(KIND=LDP) A(N),CAP,TMAX
	JRES_TRIG=0
	BAD=0
	TMAX=0.0_LDP
	DO I=1,N
	  IF(.NOT. IEEE_IS_FINITE(A(I)))THEN
	    BAD=I
	    IF(IEEE_IS_NAN(A(I)))THEN
	      JRES_TRIG=1
	    ELSE
	      JRES_TRIG=2
	    END IF
	    RETURN
	  END IF
	END DO
	DO I=1,N
	  IF(ABS(A(I)) .GT. TMAX)THEN
	    TMAX=ABS(A(I))
	    BAD=I
	  END IF
	END DO
	IF(TMAX .GT. CAP)JRES_TRIG=3
	RETURN
	END
!
! JRES_FIN: finite-only check for fields whose validity has no
! magnitude bound (e.g. GEDD, RSQHNU snapshots).
!
	LOGICAL FUNCTION JRES_FIN(A,N)
	USE SET_KIND_MODULE
	USE IEEE_ARITHMETIC
	IMPLICIT NONE
	INTEGER N,I
	REAL(KIND=LDP) A(N)
	JRES_FIN=.TRUE.
	DO I=1,N
	  IF(.NOT. IEEE_IS_FINITE(A(I)))THEN
	    JRES_FIN=.FALSE.
	    RETURN
	  END IF
	END DO
	RETURN
	END
