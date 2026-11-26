!
! This routine solves for the mean intensity as a function of depth for
! a plane-parallel atmosphere using the Feautrier Technique. A Schuster
! or diffusion approaximation is used for the lower boundary condition.
! This routine may need to be in a loop so that the f values are iterated
! to convergence.
!
	SUBROUTINE FCOMP_PP_V2(R,NEWRJ,NEWRK,SOURCE,CHI,IPLUS,
	1             HBC_J,HBC_S,INBCNEW,DBB,IC,THK,DIFF,ND,NP,METHOD)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 30-Mar-2003
!
	INTEGER ND
	INTEGER NP
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) NEWRJ(ND)
	REAL(KIND=LDP) NEWRK(ND)
	REAL(KIND=LDP) SOURCE(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) IPLUS(NP)
!
! Locally defined arrays.
!
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND),XM(ND)
	REAL(KIND=LDP) DTAU(ND),DTAU_RAD(ND)
	REAL(KIND=LDP) dCHIdR(ND)
!
	REAL(KIND=LDP), SAVE, ALLOCATABLE :: MU(:)
	REAL(KIND=LDP), SAVE, ALLOCATABLE :: JQW(:)
	REAL(KIND=LDP), SAVE, ALLOCATABLE :: KQW(:)
	INTEGER, SAVE :: NP_SAV=0
	LOGICAL, SAVE :: FIRST_TIME
!
	REAL(KIND=LDP), PARAMETER :: RZERO=0.0_LDP
	REAL(KIND=LDP), PARAMETER :: RONE=1.0_LDP
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER I,LS,NP_LOC
	REAL(KIND=LDP) DBB,DBC,IBOUND,TOR,HBC_J,HBC_S,INBCNEW,IC,E1,E2,E3
	LOGICAL DIFF,THK
	CHARACTER*6 METHOD
!
	IF(FIRST_TIME .OR. NP .NE. NP_SAV)THEN
          IF(ALLOCATED(MU))DEALLOCATE(MU,JQW,KQW)
          ALLOCATE (MU(NP),JQW(NP),KQW(NP))
	  CALL GAULEG(RZERO,RONE,MU,JQW,NP)
          KQW=JQW*MU*MU
	  NP_SAV=NP
	  FIRST_TIME=.FALSE.
	END IF
!
! Zero all arrays
!
	NEWRJ(1:ND)=0.0_LDP
	NEWRK(1:ND)=0.0_LDP
	HBC_J=0.0_LDP
	HBC_S=0.0_LDP
	INBCNEW=0.0_LDP
!
! Compute radial optical depth scale.
!
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
	DO I=1,ND-1
	  DTAU_RAD(I)=0.5_LDP*(R(I)-R(I+1))*( CHI(I)+CHI(I+1)+(R(I)-R(I+1))*
	1     (dCHIdR(I+1)-dCHIdR(I))/6.0_LDP)
	END DO
!
! Enter loop for each angle.
!
	DO 2000 LS=1,NP
	  IF(MU(LS) .NE. 0)THEN
	    DBC=DBB*MU(LS)/CHI(ND)
!
! As plane parallel, we assume atmosphere is exponential in R.
!
	    IF(THK)THEN
	      TOR=CHI(1)*(R(1)-R(3))/LOG(CHI(3)/CHI(1))
	      IF(TOR .LT. 0)TOR=0.0_LDP
	      IBOUND=SOURCE(1)*(1.0_LDP-EXP(-TOR/MU(LS)))
	    ELSE
	      IBOUND=0.0_LDP
	    END IF
!
! Compute optical depth along ray.
!
	    DO I=1,ND-1
	      DTAU(I)=DTAU_RAD(I)/MU(LS)
	    END DO
!
! Compute XM. Compute T ( a tridiagonal matrix) and store it as three vectors
! TA,TB and TC . This code is a combined version of XVECD and TCOMPD.
!
	    XM(1)=-IBOUND
	    TA(1)=0.0_LDP
	    TC(1)=1.0_LDP/DTAU(1)
	    TB(1)=-1.0_LDP-TC(1)
	    DO I=2,ND-1
	      TA(I)=TC(I-1)
	      TC(I)=1.0_LDP/DTAU(I)
	      TB(I)=-0.5_LDP*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5_LDP
	    END DO
!
	    IF(DIFF)THEN
	      TA(ND)=-TC(ND-1)
	      TB(ND)=TC(ND-1)
	      XM(ND)=DBC
	    ELSE
	      TA(ND)=-TC(ND-1)
	      TB(ND)=1.0_LDP-TA(ND)
	      XM(ND)=IC
	    END IF
	    TC(ND)=0.0_LDP
!
! Solve the tridiagonal system of equations.
!
	    CALL THOMAS(TA,TB,TC,XM,ND,IONE)
!
! Update the FA and FB matrices (see notes).
!
	    DO I=1,ND
	      NEWRJ(I)=NEWRJ(I)+JQW(LS)*XM(I)
	      NEWRK(I)=NEWRK(I)+KQW(LS)*XM(I)
	    END DO
!
	    HBC_J=HBC_J+JQW(LS)*XM(1)*MU(LS)
	    HBC_S=HBC_S+JQW(LS)*IBOUND*MU(LS)
	    INBCNEW=INBCNEW + JQW(LS)*MU(LS)*
	1      (XM(ND)-(XM(ND)-XM(ND-1))/DTAU(ND-1))
	    IPLUS(LS)=2.0_LDP*(XM(1)-IBOUND)
	  ELSE
	    DO I=1,ND
	      NEWRJ(I)=NEWRJ(I)+JQW(LS)*SOURCE(I)
	      NEWRK(I)=NEWRK(I)+KQW(LS)*SOURCE(I)
	    END DO
	    HBC_J=HBC_J+JQW(LS)*SOURCE(1)*MU(LS)
	    HBC_S=HBC_S+JQW(LS)*IBOUND*MU(LS)
	    INBCNEW=INBCNEW + JQW(LS)*MU(LS)*
	1       (SOURCE(ND)-(SOURCE(ND)-SOURCE(ND-1))/DTAU(ND-1))
	    IPLUS(LS)=0.0_LDP
	  END IF
!
2000	CONTINUE
!
! Compute the factor for the outer boundary condition.
!
	HBC_J=HBC_J/NEWRJ(1)
	HBC_S=HBC_S/SOURCE(1)
	INBCNEW=INBCNEW/(2.0_LDP*NEWRJ(ND)-IC)
!
! Compute the new Feautrier factors.
!
	DO I=1,ND
	  NEWRK(I)=NEWRK(I)/NEWRJ(I)
	END DO
!
	RETURN
	END
