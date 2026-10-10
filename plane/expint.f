      FUNCTION expint(n,x)
	USE SET_KIND_MODULE
      INTEGER n,MAXIT
      REAL(KIND=LDP)  expint,x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.E-7_LDP,FPMIN=1.E-200_LDP,EULER=.5772156649_LDP)
      INTEGER i,ii,nm1
      REAL(KIND=LDP) a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1_LDP)))then
        pause 'bad arguments in expint'
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0._LDP)then
        expint=1._LDP/nm1
      else if(x.gt.1._LDP)then
        b=x+n
        c=1._LDP/FPMIN
        d=1._LDP/b
        h=d
        do 11 i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1._LDP/(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1._LDP).lt.EPS)then
            expint=h*exp(-x)
            return
          endif
11      continue
        pause 'continued fraction failed in expint'
      else
        if(nm1.ne.0)then
          expint=1._LDP/nm1
        else
          expint=-log(x)-EULER
        endif
        fact=1.
        do 13 i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-EULER
            do 12 ii=1,nm1
              psi=psi+1._LDP/ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          if(abs(del).lt.abs(expint)*EPS) return
13      continue
        pause 'series failed in expint'
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 5.)61&0X40(9p#ms21.1.
