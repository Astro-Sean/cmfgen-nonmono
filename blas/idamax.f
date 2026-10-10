      integer function idamax(n,dx,incx)
      use set_kind_module
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      REAL(KIND=LDP) dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = abs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = abs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = abs(dx(1))
      do 30 i = 2,n
         if(abs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = abs(dx(i))
   30 continue
      return
      end
