	SUBROUTINE CMF_BLKBAND_PIN_NULL_V1(C_MAT,RHS,N,N_PIN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: N
	INTEGER, INTENT(OUT) :: N_PIN
	REAL(KIND=LDP), INTENT(INOUT) :: C_MAT(N,N)
	REAL(KIND=LDP), INTENT(IN) :: RHS(N)
	INTEGER I
C
C Only a zero row and column with zero RHS form an independent 0=0
C equation. Assigning its free correction to zero leaves all other
C equations unchanged; singular or inconsistent systems still fail
C through CMF_BLKBAND_V3's existing factorization path.
C
	N_PIN=0
	DO I=1,N
	  IF(MAXVAL(ABS(C_MAT(I,1:N))) .EQ. 0.0_LDP .AND.
	1    MAXVAL(ABS(C_MAT(1:N,I))) .EQ. 0.0_LDP .AND.
	2    RHS(I) .EQ. 0.0_LDP)THEN
	    C_MAT(I,I)=1.0_LDP
	    N_PIN=N_PIN+1
	  END IF
	END DO
	RETURN
	END
