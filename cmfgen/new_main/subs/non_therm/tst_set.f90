	PROGRAM TST_SET
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: N=51
	INTEGER I
	REAL(KIND=LDP) X0,X1
	REAL(KIND=LDP) X(N)
	REAL(KIND=LDP) dX(N)
!
	X0=2.5
	X1=112
	call set_xkt_array(x0,x1,n,x,dx,'lin')
!
	do i=1,n
	  write(100,*)i,x(i),dx(i)
	end do
!
	STOP
	END
