! Stub GET_FG13_ALO_DIAG / GET_FG13_ALO_TRI for the reference-side
! binary in test_fg13.py.  The ALO store is fork-only code - the
! pristine fg_j_cmf_v13.f has no getters - so the ref driver links
! these instead and gets VALID=.FALSE., which the test asserts to
! document that the store is absent upstream.
!
      SUBROUTINE GET_FG13_ALO_DIAG(DJDCHI,DJDETA,N,VALID)
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER N
      REAL(KIND=LDP) DJDCHI(N),DJDETA(N)
      LOGICAL VALID
      VALID=.FALSE.
      RETURN
      END
!
      SUBROUTINE GET_FG13_ALO_TRI(DJDCHI_M,DJDETA_M,DJDCHI_P,DJDETA_P, &
                                  N,VALID)
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER N
      REAL(KIND=LDP) DJDCHI_M(N),DJDETA_M(N)
      REAL(KIND=LDP) DJDCHI_P(N),DJDETA_P(N)
      LOGICAL VALID
      VALID=.FALSE.
      RETURN
      END
