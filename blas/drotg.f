      subroutine drotg(da,db,c,s)
	USE SET_KIND_MODULE
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      REAL(KIND=LDP) da,db,c,s,roe,scale,r,z
c
      roe = db
      if( abs(da) .gt. abs(db) ) roe = da
      scale = abs(da) + abs(db)
      if( scale .ne. 0.0d0 ) go to 10
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         z = 0.0d0
         go to 20
   10 r = scale*sqrt((da/scale)**2 + (db/scale)**2)
      r = sign(r,roe)     !sign(1.0d0,roe)*r
      c = da/r
      s = db/r
      z = 1.0d0
      if( abs(da) .gt. abs(db) ) z = s
      if( abs(db) .ge. abs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
   20 da = r
      db = z
      return
      end
