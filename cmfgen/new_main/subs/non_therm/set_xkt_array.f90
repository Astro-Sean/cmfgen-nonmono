!
! Subroutine that sets the XKT vector (which contains the electron
! energies in eV) for the computation of the non-thermal electron
! spectrum.
!
! Liner or logarithmic spacing can be set -- linear is recommened.
! The routine also retuns the quadrature weights for integratinng
! over the energy grid.
!
    subroutine set_xkt_array(x0,x1,n,x,dx,method)
	USE SET_KIND_MODULE
      implicit none
!
! Created 13-Sep-2011 : Based on set_array (from Luc).
!
      integer n              !Number of points in grid
      real(kind=LDP) x0 	     !Minimum energy in eV
      real(kind=LDP) x1              !Maximm energy in eV
      real(kind=LDP) x(n)            !Energy vector (returned)
      real(kind=LDP) dx(n)           !Quadrature weight
      character*3 method     ! `lin' or 'log'
!
      integer i
      real(kind=LDP) t1,x11
!
      if ((method.ne.'lin').and.(method.ne.'log')) then
         write(6,*) ' method in set array not properly set'
         stop
      endif
!
      if (n.eq.1) then
        write(6,*) 'n must be greater than 1'
        stop
      else if (method.eq.'lin') then
!
! Linear spacing and trapazoidal rule
!
        x(1) = x0
        t1 = (x1-x0) / dble(n-1)
        do i=2,n
          x(i) = x0 +(I-1)*t1
        enddo
        dx(2:n-1) = t1
        dx(1) = t1/2.0_LDP
        dx(n) = t1/2.0_LDP
!
! Logarithmic spacing and trapazoidal rule.
! A smaller step size is used near the maximum energy.
!
      else
        if (x0.eq.0.0_LDP) then
          write(6,*) 'x0 is zero = logarithmic grid not possible'
          stop
        endif
        x11 = 0.99999_LDP*x1
        t1 = log(x11/x0) / dble(n-2)
        x(n) = x1
        dx(n) = x1-x11
        x(n-1) = x11
        do i=n-2,2,-1
          x(i) = x(i+1) * exp(-t1)
        enddo
	x(1)=x0
        dx(1) = 0.5_LDP*(x(2)-x(1))
        do i=2,n-1
          dx(i) = 0.5_LDP*(x(i+1)- x(i-1))
        enddo
      endif
!
      return
    end subroutine set_xkt_array
