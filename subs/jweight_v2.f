!
! Subroutine to compute the quadrature weights for a modified cubic rule.
! Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
! The routine makes specific assumptions concerning the behaviour of the
! mean intensity at the boundary, and would need to be moified to obtain
! other quantities other than the mean intensity. The program assumes that
! the first point corresponds to mu=1.0 .
!
! D1, and R2 are work vectors.
!
	SUBROUTINE JWEIGHT_V2(X,dX,W,N)
	USE SET_KIND_MODULE
!
! Altered 18-Jun-2023 - Fixed typo.
! Created 14-Jun-2023 - Changed to V2. dX added to call.
! Altered 24-May-1996 - DP_ZERO deleted
!                       ERROR_LU etc installed.
! Altered 18-Feb-1987 - Quadratic extrapolation to zero if X(N) .NE. 0
! Altered 25-Feb-1987 - Normalization implemented.
!
	IMPLICIT NONE
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
	W(1)=0.5_LDP*H-H*H/12.0_LDP/RF-HN*H/RF/12.0_LDP
	W(2)=0.5_LDP*H+H*H/RF/6.0_LDP*(1.0_LDP+0.5_LDP*H/HN)+HN*H/RF/12.0_LDP
	W(3)=H*H/RF/12.0_LDP*(-1.0_LDP-H/HN)
!
	DO I=3,N-2
	  H=dX(I-1)
	  RF=dX(I-2)+dX(I-1)
	  RE=dX(I-1)+dX(I)
	  W(I-2)=W(I-2)-H*H/12.0_LDP/RF	
	  W(I-1)=W(I-1)+0.5_LDP*H+H*H/12.0_LDP/RE
	  W(I)=W(I)+0.5_LDP*H+H*H/12.0_LDP/RF	
	  W(I+1)=W(I+1)-H*H/12.0_LDP/RE
	END DO
!
	IF(X(N) .EQ. 0)THEN
!
! Assumes a quadratic variation for U. The is the integral from X(N-2) to zero.
! Subtracting the X to get dX is fine for high N.
!
	  RF=X(N-3)-X(N-1)
	  H=X(N-2)-X(N-1)
	  HN=X(N-1)-X(N)
	  W(N-3)=W(N-3)-H*H/RF/12.0_LDP
	  W(N-2)=W(N-2)+0.5_LDP*H
	  W(N-1)=W(N-1)+0.5_LDP*(H+2.0_LDP*HN/3.0_LDP)+H*H/HN/6.0_LDP+H*H/12.0_LDP/RF
	  W(N)=W(N)+2.0_LDP*HN/3.0_LDP-H*H/HN/6.0_LDP
	ELSE
!
! For outer boundary. Extrapolate to zero assuming a quadratic
! variation. i.e. U=A + B*mu*mu. Issue a warning indicating
! this. This integration is from X(N-1) to zero. The quadrature
! weight for X(N-1) will be negative. Due to the limit on the
! early DO loop, we also need to compute the integral form X(N-2) to X(N-1).
!
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - Angle extrapolation in JWEIGHT required'
!
! Inegral from X(N-2) to X(N-1)
!
	  H=X(N-2)-X(N-1)
	  RF=X(N-3)-X(N-1)
	  RE=X(N-2)-X(N)
	  W(N-3)=W(N-3)-H*H/12.0_LDP/RF	
	  W(N-2)=W(N-2)+0.5_LDP*H+H*H/12.0_LDP/RE
	  W(N-1)=W(N-1)+0.5_LDP*H+H*H/12.0_LDP/RF	
	  W(N)=W(N)-H*H/12.0_LDP/RE
!
! Integral from X(N-1) to zero.
!
	  RF=(X(N-1)-X(N))*(X(N-1)+X(N))
	  W(N-1)=W(N-1)-X(N-1)*(X(N)*X(N)-X(N-1)*X(N-1)/3.0_LDP)/RF
	  W(N)=W(N)+(X(N-1)**3)*2.0_LDP/RF/3.0_LDP
	END IF
!
! Ensure that the weights have the correct normalization (shouldnt be
! necessary.
!
	SUM=0.0_LDP
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	SUM=1.0_LDP/SUM
	IF(ABS(SUM-1.0_LDP) .GT. 1.0E-12_LDP)THEN
	   LUER=ERROR_LU()
	   WRITE(LUER,*)' Warning - weights were normalized in JWEIGHT'
	END IF
	DO I=1,N
	  W(I)=W(I)*SUM
	END DO
!
	RETURN
	END
