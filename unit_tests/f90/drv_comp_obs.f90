!
! Driver for COMP_OBS_V2 regression tests.  One case per argv(1):
!
!   const : constant IPLUS field, VINF=0 (zero-velocity limit)
!   lin   : IPLUS = A*NU + B, relativistic transform on
!   norel : same linear field, DO_FULL_REL=F (frequency shift only)
!   rmax2 : constant field with 2x RMAX (flux scales x4 vs const)
!   nan   : one non-finite IPLUS entry -> expect STOP diagnostic
!   range : OBS_FREQ*FREQ_CONV_FAC(mu=0) > NU_STORE(1) -> STOP
!
! Prints 'OBSF <k> <flux>' for each observer frequency so the Python
! test can compare against independently evaluated expectations.
!
	PROGRAM DRV_COMP_OBS
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NP=4, NST=64, N_OBS=4, NCMF=21
	REAL(KIND=LDP) :: NEW_IPLUS(NP), IPLUS_STORE(NST,NP), NU_STORE(NST)
	REAL(KIND=LDP) :: MU(NP), WTS(NP), OBS_FREQ(N_OBS), OBS_FLUX(N_OBS)
	REAL(KIND=LDP) :: VINF, RMAX, FL, A, B, I0, Z
	LOGICAL :: FIRST_OBS, DO_REL
	CHARACTER(LEN=32) :: MODE
	INTEGER :: I, LS
!
	MU=(/0.0_LDP,0.3_LDP,0.65_LDP,0.95_LDP/)
	WTS=(/0.20_LDP,0.30_LDP,0.30_LDP,0.20_LDP/)
	OBS_FREQ=(/57.0_LDP,55.0_LDP,53.0_LDP,51.0_LDP/)
!
	CALL GET_COMMAND_ARGUMENT(1,MODE)
	VINF=0.0_LDP
	RMAX=1.0E+04_LDP
	DO_REL=.TRUE.
	A=0.0_LDP; B=0.0_LDP; I0=2.5_LDP
!
	SELECT CASE (TRIM(MODE))
	CASE('const')
	CASE('lin')
	  VINF=1.0E+04_LDP
	  A=0.03_LDP; B=0.1_LDP
	CASE('norel')
	  VINF=1.0E+04_LDP
	  A=0.03_LDP; B=0.1_LDP
	  DO_REL=.FALSE.
	CASE('rmax2')
	  RMAX=2.0E+04_LDP
	CASE('nan')
	CASE('range')
	  VINF=3.0E+04_LDP
	  OBS_FREQ(1)=59.97_LDP
	CASE DEFAULT
	  WRITE(*,*)'unknown mode ',TRIM(MODE)
	  STOP
	END SELECT
!
	FIRST_OBS=.TRUE.
	DO I=1,NCMF
	  FL=60.0_LDP-(I-1)		!Descending CMF grid: 60..40
	  DO LS=1,NP
	    NEW_IPLUS(LS)=I0+A*FL+B
	  END DO
	  IF(TRIM(MODE) .EQ. 'nan' .AND. I .EQ. 1)THEN
	    Z=0.0_LDP
	    NEW_IPLUS(2)=1.0_LDP/Z	!+Inf: must hit the finite check
	  END IF
	  CALL COMP_OBS_V2(NEW_IPLUS,FL, &
	       IPLUS_STORE,NU_STORE,NST, &
	       MU,WTS,OBS_FREQ,OBS_FLUX,N_OBS, &
	       VINF,RMAX,'IPLUS','LIN_INT',DO_REL, &
	       FIRST_OBS,NP)
	END DO
!
	DO I=1,N_OBS
	  WRITE(*,'(A,I2,ES24.16)')'OBSF ',I,OBS_FLUX(I)
	END DO
!
	END PROGRAM DRV_COMP_OBS
