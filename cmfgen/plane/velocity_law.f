!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE VELOCITY_LAW(RVAL,ID,R,V,ND,VEL,BETA,DBETADR,GAMMA)
	USE SET_KIND_MODULE
      IMPLICIT NONE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Subroutine to calculate the velocity, beta, dbetadr and gamma
! for the given radial point for a given velocity law.  The forms are:
!
!--------------------------------------------------------------------
!

      INTEGER ND
      INTEGER ID
!
      REAL(KIND=LDP) RVAL
      REAL(KIND=LDP) R(ND)
      REAL(KIND=LDP) V(ND)
!
      REAL(KIND=LDP) VEL
      REAL(KIND=LDP) BETA
      REAL(KIND=LDP) DBETADR
      REAL(KIND=LDP) GAMMA
!
      REAL(KIND=LDP), PARAMETER :: C_KMS=2.99792458E+05_LDP
      REAL(KIND=LDP) T1
!
!--------------------------------------------------------------------
!
! We use simple linear interpolaiton. Thus this procedure is valid for
! any velocity law --- not just a hubble expansion.
!
      T1=(RVAL-R(ID+1))/(R(ID)-R(ID+1))
      VEL=T1*V(ID)+(1.0_LDP-T1)*V(ID+1)
      BETA=VEL/C_KMS
      dBETAdR=(V(ID)-V(ID+1))/(R(ID)-R(ID+1))/C_KMS
      GAMMA=1.0_LDP/SQRT(1.0_LDP-BETA*BETA)
!
      RETURN
      END
