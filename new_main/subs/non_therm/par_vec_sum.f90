!
! The simple routine computes:
!
!  Y=Y+X
!
! where Y and X are vectors of length NKT.
!
	SUBROUTINE PAR_VEC_SUM(Y,X,NKT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER NKT
	REAL(KIND=LDP) Y(NKT)
	REAL(KIND=LDP) X(NKT)
!
	INTEGER NOUT,NIN
	INTEGER J,IKT
!
	NOUT=16
        NIN=1+(NKT-1)/NOUT
!
!$OMP PARALLEL DO
        DO J=1,NOUT
          DO IKT=(J-1)*NIN+1,MIN(J*NIN,NKT)
            Y(IKT)=Y(IKT)+X(IKT)
          END DO
        END DO
!
	RETURN
	END
