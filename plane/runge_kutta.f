c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine runge_kutta(y,dydx,n,x,h,beta,dbetadr,gamma,deriv)
	USE SET_KIND_MODULE
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Numerical Recipies routine to calculate the 4th order Runge-Kutta
c solution of a set of ordinary differential equations
c
c Altered 2/21/96 DLM Updated to F90 standard
c
c--------------------------------------------------------------------
c
      implicit none
c
      external deriv
c
c Number of equations to solve: passed from calling routine
c
      integer :: n
c
c Local variables
c
      integer :: i
      integer, parameter :: nmax=2
      REAL(KIND=LDP), dimension(nmax) :: dym,dyt,yt
      REAL(KIND=LDP) h6,hh,xh
c
c  Value to be determine
c
      REAL(KIND=LDP),  dimension(n) :: dydx,y
      REAL(KIND=LDP) h,x,beta,dbetadr,gamma
c
      hh=h*0.5_LDP
      h6=h/6.0_LDP
      xh=x+hh
      x=x+h
      do i=1,n
        yt(i)=y(i)+hh*dydx(i)
      enddo
      call deriv(xh,yt,dyt,beta,dbetadr,gamma)
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)
      enddo
      call deriv(xh,yt,dym,beta,dbetadr,gamma)
      do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
      enddo
      call deriv(x,yt,dyt,beta,dbetadr,gamma)
      do i=1,n
        y(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0_LDP*dym(i))
      enddo
c
      return
      end
c
