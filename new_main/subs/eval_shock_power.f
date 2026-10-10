!
! Routine to adjust the radiative equilibrium equation for energy added to the
! gas by the magnetar. Routine taken from EVAL_RAD_DECAY_V1
!
      SUBROUTINE EVAL_SHOCK_POWER(dE_SHOCK_POWER,ND,SHOCK_POWER_FAC)
	USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      USE SHOCK_POWER_MOD
      IMPLICIT NONE
!
! Output: STEQ_T in MOD_CMFGEN is modified.
!
      INTEGER ND
      REAL(KIND=LDP) SCALE
      REAL(KIND=LDP) PI
      REAL(KIND=LDP) dE_SHOCK_POWER(ND)
      REAL(KIND=LDP) SHOCK_POWER_FAC
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the radiactive heating by 4pi.
! Also, since chi, and eta are a factor of 10^10 too large, we need to
! scale by 10^10. The BA matrix does not need to be altered.
!
      PI             = ACOS(-1.0_LDP)
      SCALE          = 1.0E+10_LDP/4.0_LDP/PI
      STEQ_T         = STEQ_T + SCALE*SHOCK_POWER*SHOCK_POWER_FAC
      dE_SHOCK_POWER = SHOCK_POWER*SHOCK_POWER_FAC

      write(6,*) 'Inc magnetar power in STEQ_T -- non-th proc OFF'

      RETURN
      END
