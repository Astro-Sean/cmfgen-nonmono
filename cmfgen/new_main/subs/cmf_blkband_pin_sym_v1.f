!
! Symmetric elimination of prescribed unknowns from a dense block.
!
! For each J with PIN(J) true, the system
!              C_MAT . X = RHS
! is reduced by fixing X(J)=XBAR(J): the J-th column's contribution
! is carried into the RHS, row J and column J are zeroed, the
! diagonal is set to 1 and RHS(J) to XBAR(J). The result is the
! standard symmetric row/column elimination for a known value --
! unlike pinning the row alone, the remaining equations no longer
! see the eliminated column, so a near-null direction cannot leak
! into the factorisation and degrade the whole block.
!
! Carries are done first against the unmodified matrix so the
! transformation is order-independent when several unknowns are
! pinned.
!
! Created 27-Nov-2026 (checklist WP-7c): companion to
! CMF_BLKBAND_PIN_NULL_V1; called from CMF_BLKBAND_V3 under
! [BA_SYM_PIN], default F.
!
	SUBROUTINE CMF_BLKBAND_PIN_SYM_V1(C_MAT,RHS,XBAR,PIN,N,N_PIN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N
	INTEGER, INTENT(OUT) :: N_PIN
	REAL(KIND=LDP), INTENT(INOUT) :: C_MAT(N,N)
	REAL(KIND=LDP), INTENT(INOUT) :: RHS(N)
	REAL(KIND=LDP), INTENT(IN) :: XBAR(N)
	LOGICAL, INTENT(IN) :: PIN(N)
	INTEGER I,J
C
	DO J=1,N
	  IF(PIN(J))THEN
	    DO I=1,N
	      RHS(I)=RHS(I)-C_MAT(I,J)*XBAR(J)
	    END DO
	  END IF
	END DO
	DO J=1,N
	  IF(PIN(J))THEN
	    DO I=1,N
	      C_MAT(I,J)=0.0_LDP
	      C_MAT(J,I)=0.0_LDP
	    END DO
	    C_MAT(J,J)=1.0_LDP
	    RHS(J)=XBAR(J)
	  END IF
	END DO
C
	N_PIN=COUNT(PIN(1:N))
	RETURN
	END
