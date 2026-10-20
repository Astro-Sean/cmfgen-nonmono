! Driver for SET_RV_HYDRO_MODEL_V3 (VEL_LAW=11) in its RDINR branch.
!
!   drv_rv_hydro_rdin <knop T|F> <nd>
!
! Reads SN_HYDRO_DATA and RDINR from the cwd, calls the production
! routine, and prints R,V,SIGMA.  Arg 1 toggles DO_NONMONO_KNOP so
! the test can compare the stock SIGMA (MON_INT_FUNS cubic
! derivative) against the item-2 rebuild (SIGMA_TRUE_CORE centred
! differences on the returned R,V pair).
!
      PROGRAM DRV_RV_HYDRO_RDIN
      USE SET_KIND_MODULE
      USE CONTROL_VARIABLE_MOD, ONLY: DO_NONMONO_KNOP
      IMPLICIT NONE
!
      INTEGER ND
      REAL(KIND=LDP), ALLOCATABLE :: R(:),V(:),SIGMA(:)
      REAL(KIND=LDP) RMAX,RCORE,RMAX_ON_RCORE,SN_AGE_DAYS
      LOGICAL PURE_HUBBLE_FLOW,RDINR
      INTEGER N_IB_INS,N_OB_INS,LU
      CHARACTER(LEN=8) :: ARG
      CHARACTER(LEN=32) :: SARG
      INTEGER I
!
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      DO_NONMONO_KNOP=(ARG(1:1) .EQ. 'T')
      CALL GET_COMMAND_ARGUMENT(2,SARG)
      READ(SARG,*)ND
!
      ALLOCATE(R(ND),V(ND),SIGMA(ND))
      R=0.0_LDP; V=0.0_LDP; SIGMA=0.0_LDP
!
      RMAX=0.0_LDP
      RCORE=0.0_LDP
      RMAX_ON_RCORE=0.0_LDP
      SN_AGE_DAYS=39.0_LDP
      PURE_HUBBLE_FLOW=.FALSE.
      N_IB_INS=0
      N_OB_INS=0
      RDINR=.TRUE.
      LU=10
!
      CALL SET_RV_HYDRO_MODEL_V3(R,V,SIGMA,RMAX,RCORE,RMAX_ON_RCORE, &
           SN_AGE_DAYS,PURE_HUBBLE_FLOW,N_IB_INS,N_OB_INS,RDINR,ND,LU)
!
      WRITE(*,'(a)')'RESULT_OK'
      DO I=1,ND
        WRITE(*,'(i6,3(1x,es24.16e3))') I,R(I),V(I),SIGMA(I)
      END DO
!
      END
