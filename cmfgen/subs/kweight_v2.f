!
! Subroutine to compute the quadrature weights for a modified cubic rule.
! Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
! The routine makes specific assumptions concerning the behaviour of the
! mean intensity at the boundary, and would need to be moified to obtain
! other quantities other than 2nd moment of the intensity. The program assumes
! that the first point corresponds to mu=1.0 .
!
	SUBROUTINE KWEIGHT_V2(X,dX,W,N)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered  14-Jun-2023 dX added to call. Changed to V2
! Altered  24-May-1996 Call to DP_ZERO removed
!                      ERROR_LU etc installed.
! Modified 25-Feb-1987 (Extrapolation to mu=0 included for so that
!                      the angle integration at the outer boundary is
!                      handled correctly. Program assumes that U is of the
!                      for U= a + b*mu*mu .
! Created 25-Nov-1986 (Based on NORDWEIGHT)
!
	INTEGER N
	REAL(KIND=LDP) X(N),dX(N),W(N)
!
! Local data
!
	REAL(KIND=LDP) H,HN,RE,RF,SUM
	INTEGER I,ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	W(:)=0.0_LDP
!
	H=dX(1)
	HN=dX(2)
	RF=H+HN
!
	W(1)=0.5_LDP*H*X(1)*(X(1)-H/3.0_LDP)
	1        -H*H/RF/12.0_LDP*((2.0_LDP+HN/H)*(1.0_LDP-H*H/10.0_LDP)
	1        -X(2)*X(2)+H*H/10.0_LDP)
	W(2)=0.5_LDP*H*X(2)*(X(2)+H/3.0_LDP)
	1        +H*H/RF/12.0_LDP
	1        *(2.0_LDP+H/HN+HN/H)*(1.0_LDP-H*H/10.0_LDP)
	W(3)=H*H/RF/12.0_LDP*(H*H/10.0_LDP-X(2)*X(2)
	1        -H/HN*(1.0_LDP-H*H/10.0_LDP))
!
	DO I=3,N-1
	  H=dX(I-1)
	  RF=dX(I-2)+dX(I-1)
	  RE=dX(I-1)+dX(I)
	  W(I-2)=W(I-2)-(X(I-1)*X(I-1)-H*H/10.0_LDP)*H*H/12.0_LDP/RF
	  W(I-1)=W(I-1)+0.5_LDP*H*X(I-1)*(X(I-1)-H/3.0_LDP)
	1          +(X(I)*X(I)-H*H/10.0_LDP)*H*H/12.0_LDP/RE	
	  W(I)=W(I)+0.5_LDP*H*X(I)*(X(I)+H/3.0_LDP)
	1          +(X(I-1)*X(I-1)-H*H/10.0_LDP)*H*H/12.0_LDP/RF	
	  W(I+1)=W(I+1)-(X(I)*X(I)-H*H/10.0_LDP)*H*H/12.0_LDP/RE
	END DO
!
! Assume that dU/dmu=0 at mu=0, and that U depend on MU*MU. Could
! also have this assumption in previous equation as well.
! Subtracting the X to get dX is fine for high N.
!
	IF(X(N) .EQ. 0)THEN
	  H=X(N-1)-X(N)
	  W(N-1)=W(N-1)+0.5_LDP*H*(X(N-1)*X(N-1)-X(N-1)*H/3.0_LDP)
	1        -H/6.0_LDP*(X(N-1)*X(N-1)-H*H/10.0_LDP)
	  W(N)=W(N)+H/6.0_LDP*(X(N-1)*X(N-1)-H*H/10.0_LDP)
	ELSE
	  H=15.0_LDP*(X(N-1)-X(N))*(X(N-1)+X(N))
	  W(N-1)=W(N-1)+2.0_LDP*(X(N-1)**5)/H
	  W(N)=W(N)+(X(N-1)**3)*(3.0_LDP*X(N-1)*X(N-1)-5.0_LDP*X(N)*X(N))/H
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - Extrapolation to zero required in KWEIGHT'
	END IF
!
! Ensure that the weights have the correct normalization (shouldnt be
! necessary.
!
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	SUM=1.0_LDP/3.0_LDP/SUM
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights required normalization in KWEIGHT'
	END IF
	DO I=1,N
	  W(I)=W(I)*SUM
	END DO
!
	RETURN
	END
