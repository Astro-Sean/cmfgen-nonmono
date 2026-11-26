      function dcabs1(z)
      USE SET_KIND_MODULE
      REAL(KIND=LDP) dcabs1
      double complex z,zz
      REAL(KIND=LDP) t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = abs(t(1)) + abs(t(2))
      return
      end
