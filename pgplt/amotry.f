      FUNCTION AMOTRY(P,Y,PSUM,MP,NP,NDIM,FUNK,IHI,FAC)
	USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER IHI,MP,NDIM,NP
      REAL(KIND=LDP) AMOTRY,FAC,P(MP,NP),PSUM(NP),Y(MP),FUNK
      EXTERNAL FUNK
!
! Altered 01-Feb-2022: PTRY is now automatically dimensioned.
!
! Uses function FUNK
!
      INTEGER J
      REAL(KIND=LDP) FAC1,FAC2,YTRY
      REAL(KIND=LDP) PTRY(NDIM)
!
      FAC1=(1.D0-FAC)/NDIM
      FAC2=FAC1-FAC
      DO J=1,NDIM
        PTRY(J)=PSUM(J)*FAC1-P(IHI,J)*FAC2
      END DO
      YTRY=FUNK(PTRY)
!
      IF (YTRY .LT. Y(IHI)) THEN
        Y(IHI)=YTRY
        DO J=1,NDIM
          PSUM(J)=PSUM(J)-P(IHI,J)+PTRY(J)
          P(IHI,J)=PTRY(J)
	END DO
      END IF
      AMOTRY=YTRY
      RETURN
!
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 5
