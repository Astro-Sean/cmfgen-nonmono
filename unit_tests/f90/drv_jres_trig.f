! drv_jres_trig.f - manufactured validity tests for JRES_TRIG
! (obs/jres_validity.f, Phase 2.1A).  Each case prints one RESULT
! line: name, trigger code, offending index, max|a|.  The python
! wrapper asserts the classification of converged, NaN (interior
! and boundary), Inf, over-cap, and first-frequency-fatal inputs,
! plus the finite-but-unconverged pass case.
!
      PROGRAM DRV_JRES_TRIG
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=10
      REAL(KIND=LDP) A(ND),CAP,TMAX
      INTEGER JRES_TRIG,BAD,I,K
      EXTERNAL JRES_TRIG
      CHARACTER(LEN=12) NAME
!
! Case matrix: name, cap, then the field pattern.
! bad_case codes: 0 all finite under cap, 1 NaN at I, 2 +Inf at I,
! 3 finite above cap at I, 4 finite-but-'unconverged' (large but
! under cap).
!
      DO K=1,7
        A(1:ND)=1.0E-5_LDP
        CAP=1.0E-3_LDP
        BAD=0
        NAME=' '
        SELECT CASE(K)
        CASE(1)
          NAME='converged'
        CASE(2)
          NAME='nan_interior'
          A(5)=0.0_LDP/0.0_LDP
        CASE(3)
          NAME='nan_first'
          A(1)=0.0_LDP/0.0_LDP
        CASE(4)
          NAME='inf_interior'
          A(7)=1.0_LDP/0.0_LDP
        CASE(5)
          NAME='overcap_int'
          A(6)=1.0E25_LDP
        CASE(6)
          NAME='fedd_over2'
          A(3)=2.5_LDP
          CAP=2.0_LDP
        CASE(7)
          NAME='unconverged'
          A(4)=9.9E-4_LDP
        END SELECT
        I=JRES_TRIG(A,ND,CAP,BAD,TMAX)
        WRITE(6,'(A,1X,I1,1X,I3,1X,ES12.4)')TRIM(NAME),I,BAD,TMAX
      END DO
      END
