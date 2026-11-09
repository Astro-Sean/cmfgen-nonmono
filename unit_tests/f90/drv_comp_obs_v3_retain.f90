!
! Driver for COMP_OBS_V3 NU_STORE retention regression (09-Nov-2026).
! Same scenario as drv_comp_obs_retain.f90 but calls COMP_OBS_V3, the
! variant cmf_flux_sub_v5 actually uses; the identical MAX_FREQ_CONV_FAC
! retention fix lives in both routines and V3 had no coverage.
!
! Prints 'OBSF <k> <flux>' per observer frequency for Python comparison.
!
	PROGRAM DRV_COMP_OBS_V3_RETAIN
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NP=5, NST=64, N_OBS=40, NCMF=250
	REAL(KIND=LDP) :: NEW_IPLUS(NP), IPLUS_STORE(NST,NP), NU_STORE(NST)
	REAL(KIND=LDP) :: P(NP), MU(NP), WTS(NP)
	REAL(KIND=LDP) :: OBS_FREQ(N_OBS), OBS_FLUX(N_OBS)
	REAL(KIND=LDP) :: VINF, RMAX, FL, I0
	LOGICAL :: FIRST_OBS
	INTEGER :: I, LS
!
	MU=(/0.0_LDP,0.3_LDP,0.65_LDP,0.95_LDP,-0.5_LDP/)
	P =MU
	WTS=(/0.20_LDP,0.25_LDP,0.25_LDP,0.20_LDP,0.10_LDP/)
!
! Descending observer grid 110 -> 60 (linear, 40 bins).
!
	DO I=1,N_OBS
	  OBS_FREQ(I)=110.0_LDP-(I-1)*(50.0_LDP/(N_OBS-1))
	END DO
!
	VINF=3.0E+04_LDP			!beta~0.1 -> conv range ~[0.91,1.055]
	RMAX=1.0E+04_LDP
	I0=2.5_LDP
!
	FIRST_OBS=.TRUE.
	OBS_FLUX=0.0_LDP
	DO I=1,NCMF
	  FL=120.0_LDP*1.004_LDP**(-(I-1))	!Descending CMF grid 120->~44
	  DO LS=1,NP
	    NEW_IPLUS(LS)=I0+0.03_LDP*FL+0.1_LDP
	  END DO
	  CALL COMP_OBS_V3(NEW_IPLUS,FL, &
	       IPLUS_STORE,NU_STORE,NST, &
	       P,MU,WTS,OBS_FREQ,OBS_FLUX,N_OBS, &
	       VINF,RMAX,'IPLUS','LIN_INT',.TRUE., &
	       FIRST_OBS,.FALSE.,NP)
	END DO
!
	DO I=1,N_OBS
	  WRITE(*,'(A,I2,ES24.16)')'OBSF ',I,OBS_FLUX(I)
	END DO
!
	END PROGRAM DRV_COMP_OBS_V3_RETAIN
