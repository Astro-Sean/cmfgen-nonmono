! Subroutine to call the routines necessary to set up the
! gamma-ray routine.
!
!
	SUBROUTINE GAMRAY_SUB_V3(ND,NC,NP,P,R,V,SIGMA,VDOP_VEC,CLUMP_FAC,    &
		MU_AT_RMAX,HQW_AT_RMAX,DELV_FRAC_FG,REXT_FAC,METHOD,INST_DECAY,SN_AGE_DAYS)
	USE SET_KIND_MODULE
	USE MOD_GAMMA_V3
        USE GAMMA_NUC_DECAY_V2
        USE GAM_MU_MOD
	USE MOD_RD_GAMRAY_CNTRL_VARIABLES
	USE NUC_ISO_MOD, ONLY : RADIOACTIVE_DECAY_ENERGY
	IMPLICIT NONE
!
! Altered 20-May-2023 : When now check if the data directory (data) exists.
! Altered 02-Jun-2022 : Fixed bug with format statement (extra )).
! Altered 21-Nov-2021 : Added access to RADIOACTIVE_DECAY_ENERGY
! Altered 19-Nov-2021 : Added SN_AGE_DAYS to call. Added SN_AGE_DAYS to gamma_energy_dep_v7 call.
! Altered 18-Nov-2021 : THK_CONT initialized to FALSE.
! Altered 17-Mar-2019 : Fixed allocation of TC_WRK (DJH)
!
	INTEGER :: I,J,K
	INTEGER :: FI,ID
	INTEGER :: IOS
	INTEGER :: ND,NC,NP
	INTEGER :: NA,NA_MON
	INTEGER, PARAMETER :: IZERO =0
	INTEGER, PARAMETER :: IONE  =1
	INTEGER, PARAMETER :: ITWO  =2
	INTEGER, PARAMETER :: ITHREE=3
	INTEGER, PARAMETER :: IFOUR =4
	INTEGER, PARAMETER :: IFIVE =5
	INTEGER, PARAMETER :: ISIX  =6
	REAL(KIND=LDP) :: T1,T2,T3,T4,T5
	REAL(KIND=LDP) :: P(NP)
	REAL(KIND=LDP) :: R(ND)
	REAL(KIND=LDP) :: V(ND)
	REAL(KIND=LDP) :: SIGMA(ND)
	REAL(KIND=LDP) :: VDOP_VEC(ND)
	REAL(KIND=LDP) :: CLUMP_FAC(ND)
	REAL(KIND=LDP) :: MU_AT_RMAX(NP)
	REAL(KIND=LDP) :: HQW_AT_RMAX(NP)
!
	REAL(KIND=LDP) :: SN_AGE_DAYS
	REAL(KIND=LDP) :: DELV_FRAC_FG
	REAL(KIND=LDP) :: REXT_FAC
	REAL(KIND=LDP) :: CHI_SCAT_CLUMP(ND)
	REAL(KIND=LDP) :: PI
	REAL(KIND=LDP) :: FOURPI
!
	REAL(KIND=LDP), ALLOCATABLE :: TA_WRK(:)
	REAL(KIND=LDP), ALLOCATABLE :: TB_WRK(:)
	REAL(KIND=LDP), ALLOCATABLE :: TC_WRK(:)
!
	LOGICAL :: INST_DECAY
	LOGICAL :: THK_CONT=.FALSE.
	CHARACTER(LEN=6) :: METHOD
	CHARACTER(LEN=40) :: FILENAME
	CHARACTER(LEN=30) :: OPTION
	CHARACTER(LEN=80) :: STRING
!
!	REAL(KIND=LDP), PARAMETER :: PLANCK = 4.135668E-18 	! keV*s
	REAL(KIND=LDP) :: H_PL 	! keV*s
	REAL(KIND=LDP), PARAMETER :: SOL = 299792 		! Units of km/s
	REAL(KIND=LDP), PARAMETER :: ERGS_TO_MEV = 624150.9_LDP	
	REAL(KIND=LDP), PARAMETER :: HoMC2 = 8.09330118E-21_LDP  	! HoMC2 = h/mc^2 in units of seconds
!
! Transfer routine variables
!
	REAL(KIND=LDP) :: dLOG_NU
	REAL(KIND=LDP) :: BNUE
	REAL(KIND=LDP) :: DBB
	REAL(KIND=LDP) :: FL
	REAL(KIND=LDP) :: HFLUX_AT_IB
	REAL(KIND=LDP) :: HFLUX_AT_OB
	REAL(KIND=LDP) :: IPLUS(NP)
	REAL(KIND=LDP) :: FEDD(NP)
!^
	INTEGER :: LUER
	INTEGER :: ERROR_LU
	EXTERNAL :: ERROR_LU
	INTEGER, PARAMETER :: LU_GAM=57
	INTEGER, PARAMETER :: LU_GAMOBS=88
!
	LOGICAL :: FILE_EXISTS
	LOGICAL :: FIRST_OBS_COMP
	LOGICAL :: FIRST_FREQ
	LOGICAL :: NEW_FREQ
	LOGICAL :: TRUE=.TRUE., FALSE=.FALSE.
!
	PI=ACOS(-1.0_LDP)
	FOURPI=4.0_LDP*PI
	LUER=ERROR_LU()
	H_PL=PLANCK*1.0E+3_LDP
!
	INQUIRE(FILE='data',EXIST=FILE_EXISTS)
	IF(.NOT. FILE_EXISTS)THEN
	  WRITE(6,*)' Error in GAMRAY_SUB_V3'
	  WRITE(6,*)' When VERBOSE_GAMMA is TRUE the directore data has to exist'
	  WRITE(6,*)' Most of the gamma-ray diagnistics are output to the data directory'
	  STOP
	END IF
!
! RD_GAMRAY_CNTRL to read in the control parameters (GAMRAY_PARAMS) for the gamma-ray code
!
	FILENAME='GAMRAY_PARAMS'
	CALL RD_GAMRAY_CNTRL_V2(FILENAME)
!
	ALLOCATE (NU_GRID_VEC(NU_GRID_MAX))
!
	CALL GAM_NUC_DECAY_DATA_SUB_V3(VERBOSE_GAMMA)
!
	ALLOCATE (NU_VEC(N_GAMMA))
!
	NU_VEC=0.0_LDP
	CALL GAMMA_NU_GRID_V16(N_GAMMA,NU_VEC,NU_GRID_VEC,NF_GRID_PTS)
	ALLOCATE(NU_VEC_15(NF_GRID_PTS))
	NU_VEC_15(1:NF_GRID_PTS)=NU_GRID_VEC(1:NF_GRID_PTS)/1.0E15_LDP
	ALLOCATE(GAMMA_TAU(NF_GRID_PTS),XRAY_TAU(NF_GRID_PTS))
	ALLOCATE(NU_END(NF_GRID_PTS))
!	ALLOCATE(GAM_NU_MAX(NF_GRID_PTS))
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/gamma_nu_grid.dat',ACTION='WRITE',STATUS='UNKNOWN')
	  CALL WRITV_V2(NU_GRID_VEC,NF_GRID_PTS,6,'Gamma ray freq grid',7)
	  CLOSE(UNIT=7)
	END IF
!
	OPEN(UNIT=LU_GAM,FILE='kevin_testing',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	IF (IOS .NE. 0 ) STOP '*****Cannot open kevin_testing*****'
	CALL SET_LINE_BUFFERING(LU_GAM)
!
	ALLOCATE(GAM_E_ME(NF_GRID_PTS))
	ALLOCATE(GAM_NU_MAX(NF_GRID_PTS))
	ALLOCATE(KLEIN_ARRAY(NF_GRID_PTS,NF_GRID_PTS))
	ALLOCATE(ETA_ISO(ND,NF_GRID_PTS))
	ALLOCATE(GAMMA_FLUX(NF_GRID_PTS))
	ALLOCATE(DECAY_KIN_E(ND))
	ALLOCATE(ED_TOT(ND))
!
! KLEIN_SCAT is a subroutine to set up an array for the calculation of the scattering emissivity
!
	GAM_E_ME=HoMC2*NU_GRID_VEC
	CALL KLEIN_SCAT(KLEIN_ARRAY,NF_GRID_PTS,GAM_E_ME)
!
! This part is necessary for the scattering routine.
! For this next piece, I have to break it up into parts. The limit of
! backscattered photon energy is equal to half of the electron rest mass energy
! because of compton scattering. Therefore, I will take the upper limit for an outgoing
! frequency higher than half of the rest mass to be the first frequency in my vector
! which will represent infinity.
!
! Edit: this is no longer relevant since scattering emissivity is solved
! for all downscattered frequencies. 26-Oct-2017 KDW25
!
!	DO I=1,NF_GRID_PTS
!	  IF (GAM_E_ME(I) .LT. 0.5D0) THEN
!	    GAM_NU_MAX(I)=GAM_E_ME(I)/(1.0D0-2.0D0*GAM_E_ME(I))
!	  ELSE
!	    GAM_NU_MAX(I)=GAM_E_ME(1)
!	  END IF
!	END DO
!
! These next pieces are just diagnosing and supplemental pieces. They
! can be removed if needed.
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/scattering_diff.dat',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	    DO I=1,NF_GRID_PTS
	      WRITE(7,'(ES18.8,2X,ES18.8,2X,ES18.8)')NU_GRID_VEC(I),NU_GRID_VEC(I)*H_PL,&
		(NU_GRID_VEC(I)-NU_GRID_VEC(I)/(1.0_LDP+2.0_LDP*GAM_E_ME(I)))/NU_GRID_VEC(I)*SOL
	    END DO
	  CLOSE(7)
	END IF
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/velocity_step.dat',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	    DO I=1,NF_GRID_PTS-1
	      WRITE(7,'(ES18.8,2X,ES18.8,2X,ES18.8)')NU_GRID_VEC(I),NU_GRID_VEC(I)*H_PL,&
		(NU_GRID_VEC(I)-NU_GRID_VEC(I+1))/NU_GRID_VEC(I)*SOL
	    END DO
	  CLOSE(7)
	END IF
!
! ^
!
!
! Calculating the intrinsic local emission from decays
!
	CALL TUNE(IONE,'ISO EMISSION')
	CALL GAMMA_INT_EMISS_V9(ETA_ISO,NU_GRID_VEC,NF_GRID_PTS,R,ND,V_GAUSS,&
				DECAY_KIN_E,NORM_GAM_LINES,NORM_LINES_ONLY,INST_DECAY)
	CALL TUNE(ITWO,'ISO EMISSION')
!
	ALLOCATE(GAMRAY_EMISS(ND))
	ALLOCATE(TOTAL_DECAY_LUM(ND))
!
	ALLOCATE(TA_WRK(NF_GRID_PTS))
	TA_WRK=0.0_LDP
!
	GAMRAY_EMISS=0.0_LDP
	DO I=1,ND
	  DO K=1,NF_GRID_PTS
	    TA_WRK(K)=ETA_ISO(I,K)
	  END DO
	  CALL LUM_FROM_ETA(TA_WRK,NU_GRID_VEC,NF_GRID_PTS)
	  GAMRAY_EMISS(I)=SUM(TA_WRK)*FOURPI
	  TOTAL_DECAY_LUM(I)=GAMRAY_EMISS(I)+DECAY_KIN_E(I)
	END DO
!
	T3=0.0_LDP
	ALLOCATE(TB_WRK(ND))
	TB_WRK=0.0_LDP
	DO I=1,ND
	  TB_WRK(I)=TOTAL_DECAY_LUM(I)*R(I)*R(I)
	END DO
	CALL LUM_FROM_ETA(TB_WRK,R,ND)
	T3=SUM(TB_WRK)
	T3=T3*FOURPI*1.0E+30_LDP
	GAMRAY_EMISS=GAMRAY_EMISS + DECAY_KIN_E
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/gamma_ray_local_emission.dat',STATUS='UNKNOWN',&
		ACTION='WRITE',IOSTAT=IOS)
	    WRITE(7,'(A,ES20.10)')'Total Nuclear Luminosity in ergs/s:',T3
	    WRITE(7,'(A,ES20.10)')'Total Nuclear Luminosity in Lsun:',T3/3.826D+33
	    WRITE(7,'(A)')'!Gamma ray energies from the different lines and K.E. for each depth in ergs/cm^3/s'
	    WRITE(7,'(A)')'!Local emission is both line energies and positron/electron K.E.'
	    WRITE(7,'(4(A20,1X))')'! Radius','Velocity km/s','Local Emiss','Decay K.E.'
	    DO I=1,ND
	      WRITE(7,'(4(ES20.10,4X))')R(I),V(I),GAMRAY_EMISS(I),DECAY_KIN_E(I)
	    END DO
	  CLOSE(UNIT=7)
	END IF
!
	IF(VERBOSE_GAMMA)THEN
	  FILENAME='./data/ETA_ISO_'
	  OPTION='ETA_ISO'
	  CALL WRITE_ARRAY_ISO(ETA_ISO,ND,NF_GRID_PTS,NU_GRID_VEC,FILENAME)
	END IF
!
	CALL ELECTRON_DENSITY_CALC_V2(ED_TOT,ND)
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/electron_density.dat',ACTION='WRITE',&
		STATUS='UNKNOWN')
	  WRITE(7,'(A6,4X,A18)') 'Depth:','Electron density:'
	  DO I=1,ND
	    WRITE(7,'(2X,I3,5X,ES16.6)')I,ED_TOT(I)
	  END DO
	  CLOSE(UNIT=7)
	END IF
!
	NA=2*NP-1  ! Max number of mu angles
	NA_MON=NA*ANG_MULT ! Max number of linear mu grid angles
!
!
! NU_END is used to find the NU_GRID_VEC index of the back-scattered frequecny of all the frequencies in the grid
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/nu_end.dat',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Trouble opening ./data/nu_end.dat; IOSTAT=',IOS
	    STOP
	  END IF
	  WRITE(7,'(A5,2X,A5,2X,5(A16,2X))')'!  I','IEND','NU','NU_BACK_SCAT','NU(IEND-1)','NU(IEND)','NU(IEND+1)'
	END IF
	DO I=1,NF_GRID_PTS
	  T1=NU_GRID_VEC(I)/(1.0_LDP+2.0_LDP*HoMC2*NU_GRID_VEC(I))
	  J=1
	  DO WHILE(NU_GRID_VEC(J) .GT. T1)
	    J=J+1
	  END DO
	  IF(I .GE. J-1)THEN
	    NU_END(I)=I
	  ELSE IF(I .LT. J-1)THEN
	    NU_END(I)=J-1
	  END IF
	  IF(VERBOSE_GAMMA)THEN
	    IF(NU_END(I) .EQ.1)THEN
	      WRITE(7,'(I5,2X,I5,2X,2(ES16.6,2X),A18,2(ES16.6,2X))')I,NU_END(I),NU_GRID_VEC(I), &
	   	T1,' ',NU_GRID_VEC(NU_END(I)),NU_GRID_VEC(NU_END(I)+1)
	    ELSE IF(NU_END(I) .NE. 1 .AND. NU_END(I) .NE. NF_GRID_PTS)THEN
	      WRITE(7,'(I5,2X,I5,2X,5(ES16.6,2X))')I,NU_END(I),NU_GRID_VEC(I),T1, &
	        NU_GRID_VEC(NU_END(I)-1),NU_GRID_VEC(NU_END(I)),NU_GRID_VEC(NU_END(I)+1)
	    ELSE IF(I .EQ. NF_GRID_PTS .OR. NU_END(I) .EQ. NF_GRID_PTS)THEN
	      WRITE(7,'(I5,2X,I5,2X,4(ES16.6,2X),A18)')I,NU_END(I),NU_GRID_VEC(I),T1, &
	        NU_GRID_VEC(NU_END(I)-1),NU_GRID_VEC(NU_END(I)),' '
	    END IF
	  END IF
	END DO
	CLOSE(UNIT=7)
!
	ALLOCATE(GAM_I(NA,ND,NF_GRID_PTS),GAM_ETA_SCAT(NA,ND,NF_GRID_PTS))
	ALLOCATE(GAM_I_T(ND,NA,NF_GRID_PTS),GAM_ETA_SCAT_T(ND,NA,NF_GRID_PTS))
	ALLOCATE(GAMMA_J(NF_GRID_PTS,ND))
	ALLOCATE(GAM_ETA_MUAVG(NF_GRID_PTS,ND))
!
! ^^^^  "_T" means transpose on arrays to use them in the transfer routine  ^^^^
!
	GAM_I=0.0_LDP !gamma-ray intensity
	GAM_I_T=0.0_LDP !gamma-ray intensity
	GAM_ETA_SCAT=0.0_LDP !gamma-ray scattering array
	GAM_ETA_SCAT_T=0.0_LDP !gamma-ray scattering array
!
	ALLOCATE(GAM_ETA(ND,NA),GAM_INT(ND,NA))
	GAM_ETA=0.0_LDP
	GAM_INT=0.0_LDP
!
	ALLOCATE(XCHI(ND),E_SCAT_CHI(ND),GAM_OPAC(ND),GAM_OPAC_COPY(ND))
	ALLOCATE(GAM_OPAC_CLUMP(ND))
	ALLOCATE(ETA_NORM_CONST(ND))
	ALLOCATE(E_SCAT_ARRAY(NF_GRID_PTS,ND))
	ALLOCATE(CHI_ARRAY_ABS(NF_GRID_PTS,ND))
	ALLOCATE(GAM_ENERGY_DEP(ND))
	ALLOCATE(GAMRAY_LUM_VEC(ND))
	ALLOCATE(GAM_ETA_NUAVG(NA,ND))
!
	GAM_ENERGY_DEP=0.0_LDP
	ETA_NORM_CONST=1.0_LDP
	XCHI=0.0_LDP
	E_SCAT_CHI=0.0_LDP
	GAM_OPAC=0.0_LDP
	GAM_OPAC_COPY=0.0_LDP
	GAM_OPAC_CLUMP=0.0_LDP
!
	WRITE(LUER,'(/,A)')' Set up gamma-ray transfer routine variables...'
!
!************************************************************************************
!************************************************************************************
!
! Now to solve the radiative transfer equation
!
!************************************************************************************
!************************************************************************************
!
	FIRST_FREQ=.FALSE.
	dLOG_NU=0.0_LDP
	DO FI=1,NF_GRID_PTS
	  WRITE(LU_GAM,'(A5,I7)')"FI:",FI
	  BNUE=0.0_LDP
	  DBB=0.0_LDP
	  GAM_INT=0.0_LDP
	  GAM_ETA=0.0_LDP
	  FL=NU_GRID_VEC(FI)/1.0E15_LDP
	  IF(FI .EQ. 1)THEN
	    FIRST_FREQ=.TRUE.
	  ELSE
	    FIRST_FREQ=.FALSE.
	    dLOG_NU=LOG(NU_GRID_VEC(FI-1)/NU_GRID_VEC(FI))
	  END IF
	  NEW_FREQ=.TRUE.
!
! Setting up the opacities each frequencies
!
!         CALL GAMMA_XCROSS_V3(NU_GRID_VEC(FI),ND,XCHI,E_SCAT_CHI,ED_TOT)
	  CALL GAMMA_XCROSS_V4_TEST(NU_GRID_VEC(FI),ND,XCHI,E_SCAT_CHI,ED_TOT)
!
! E_SCAT_CHI is not in units CMFGEN uses so I multiply by the factor of 10^10
!
	  GAM_OPAC(1:ND)=XCHI(1:ND) + 1.0E+10_LDP*E_SCAT_CHI(1:ND)
	  GAM_OPAC_CLUMP(1:ND)=GAM_OPAC(1:ND)*CLUMP_FAC(1:ND)
	  GAM_OPAC_COPY=GAM_OPAC
!
! I store the separate opacities as a function of frequency and depth to
! use later on after Intensity is solved. Note these arrays are in cm^-1
!
	  E_SCAT_ARRAY(FI,1:ND)=E_SCAT_CHI(1:ND)
	  CHI_ARRAY_ABS(FI,1:ND)=XCHI(1:ND)*1D-10
!
!
! Using the array GAM_ETA to be the 2D array to pass to the subroutines
! along with GAM_INT (for gamma intensity). Initialized the total
! emissivity as the scattered plus always the isotropic emission.
!
	  DO J=1,NA ! NA=2*NP-1
	    DO I=1,ND
	      GAM_ETA(I,J)=ETA_ISO(I,FI)+GAM_ETA_SCAT(J,I,FI)
	    END DO
	  END DO
	  GAM_ETA=GAM_ETA*1.0E+10_LDP
!
! Now we solve for the intensity for all rays for the current frequency
! We are solving it from blue to red since Compton scattering downgrades
! photon energy.
!
	  CALL TUNE(IONE,'CMF_FORMAL SOLVER')
	  CALL CMF_FORMAL_REL_V4_GAM( &
		GAM_OPAC_COPY,GAM_OPAC_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P, &
		GAM_OPAC_COPY,FEDD,HFLUX_AT_IB,HFLUX_AT_OB,IPLUS, &
		FL,dLOG_NU,BNUE,DBB, &
		GAM_IB_METH,THK_CONT, &
		VDOP_VEC,DELV_FRAC_FG,REXT_FAC, &
		METHOD,FIRST_FREQ,NEW_FREQ,NC,NP,ND, &
		NA,FI,GAM_ETA,GAM_INT,ANG_MULT)
	  CALL TUNE(ITWO,'CMF_FORMAL SOLVER')
!
! We need to save the intensity from GAM_INT(ND,NA)
!
	  DO I=1,ND
	    DO J=1,NA
	      GAM_I(J,I,FI)=GAM_INT(I,J)
	    END DO
	  END DO
!
! Next we calculate the scattering emissivity from all downscattered
! photons from the current frequency (assuming no emission at the
! current frequency)
!
	  CALL ETA_SCAT_V6(NF_GRID_PTS,NU_GRID_VEC,FI,NU_END(FI),&
		GAM_ETA_SCAT,GAM_I,ND,NA,NA_MON,KLEIN_ARRAY,&
		GAM_E_ME,CHEB_ORDER,ED_TOT,DO_NORM_ETA)
!
	END DO ! FI (over frequency)
!
!************************************************************************************
!************************************************************************************
!
! We have now solved for the specific intensity, so we can now get an
! observers frame flux. Note this code currently uses the same CMF
! frequencies as it does the OBS frequencies
!
	GAMMA_FLUX=0.0_LDP
	ALLOCATE(GAMMA_NU_STORE(NF_GRID_PTS))
	ALLOCATE(GAM_I_STORE(NF_GRID_PTS,NP))
	GAMMA_NU_STORE=0.0_LDP
	GAM_I_STORE=0.0_LDP
	FIRST_OBS_COMP=.TRUE.
!	OBS_REL_FULL=.TRUE.
	DO I=1,NF_GRID_PTS
	  FL=NU_VEC_15(I)
	  CALL COMP_OBS_V2(GAM_I(1:NP,1,I),FL, &
		GAM_I_STORE,GAMMA_NU_STORE,NF_GRID_PTS, &
		MU_AT_RMAX,HQW_AT_RMAX,NU_VEC_15,GAMMA_FLUX,NF_GRID_PTS, &
		V(1),R(1),'IPLUS','MON_INT',OBS_REL_FULL, &
		FIRST_OBS_COMP,NP)
	END DO
!
! The gamma-ray flux is now in Janskies for 1 kpc. I will output this
! the same as OBSFLUX does for CMFGEN
!
! Later, I will try add an option in plt_spec and like files to change the units to
! photons/cm^2/s/keV when using rd_mod in plt_spec. For now use
! cnvrt_gamflux.exe
!
	CALL GEN_ASCI_OPEN(LU_GAMOBS,'GAMFLUX','UNKNOWN',' ',' ',IZERO,IOS)
	WRITE(STRING,'(I10)')NF_GRID_PTS
	DO WHILE(STRING(1:1) .EQ. ' ')
	  STRING(1:)=STRING(2:)
	END DO
	STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	CALL WRITV_V2(NU_VEC_15,NF_GRID_PTS,ISIX,TRIM(STRING),LU_GAMOBS)
	CALL WRITV_V2(GAMMA_FLUX,NF_GRID_PTS,IFOUR,&
		'Observed intensity (Janskys)',LU_GAMOBS)
	CLOSE(UNIT=LU_GAMOBS)
!
!
!------------------------------------------------------------------------------
! Now I write some variables to check the outputs
!------------------------------------------------------------------------------
!
! Checking the optical depth to both x-rays and gamma-rays
!
	GAMMA_TAU=0.0_LDP
	XRAY_TAU=0.0_LDP
	DO K=1,NF_GRID_PTS
	  DO I=1,ND-1
	    T1=E_SCAT_ARRAY(K,I)
	    T2=E_SCAT_ARRAY(K,I+1)
	    T3=CHI_ARRAY_ABS(K,I)
	    T4=CHI_ARRAY_ABS(K,I+1)
	    GAMMA_TAU(K)=GAMMA_TAU(K)+0.5_LDP*(R(I)-R(I+1))*(T1+T2)
	    XRAY_TAU(K)=XRAY_TAU(K)+0.5_LDP*(R(I)-R(I+1))*(T3+T4)
	  END DO
	END DO
	GAMMA_TAU=GAMMA_TAU*1.0E10_LDP
	XRAY_TAU=XRAY_TAU*1.0E10_LDP
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/TAU_gam_xray.dat',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(7,'(A1,1X,A12,2X,A14,2X,A14)')'!','NU (keV)','TAU_GAM','TAU_XRAY'
	  DO I=1,NF_GRID_PTS
	    WRITE(7,'(2X,F14.6,ES16.6,ES16.6)')NU_GRID_VEC(I)*H_PL,GAMMA_TAU(I),XRAY_TAU(I)
	  END DO
	  CLOSE(UNIT=7)
	END IF
!
! Now we calculate the energy deposited from gamma-ray scattering.
! This calculates the difference between integral (chi*I-eta) dnu*dOmega
!
	CALL GAMMA_ENERGY_DEP_V7(GAM_I,GAM_ETA_SCAT,E_SCAT_ARRAY,CHI_ARRAY_ABS, &
		GAM_ENERGY_DEP,DECAY_KIN_E,NU_GRID_VEC,NF_GRID_PTS,V,R,NA,ND, &
	        GAMRAY_EMISS,SN_AGE_DAYS)
	RADIOACTIVE_DECAY_ENERGY=GAM_ENERGY_DEP
!
	K=0
	DO I=1,ND
	  K=MAX(R_MU(I)%MU_PTS,K)
	END DO
	ALLOCATE(TC_WRK(K))
!
	IF(WRITE_J_DATA)THEN
	  GAMMA_J=0D0
	  DO K=1,NF_GRID_PTS
	    DO I=1,ND
	      ID=R_MU(I)%MU_PTS
	      TC_WRK=0.0_LDP
	      DO J=1,ID
	        TC_WRK(J)=GAM_I(J,I,K)
	      END DO
	      CALL LUM_FROM_ETA(TC_WRK,R_MU(I)%MU_VECTOR,ID)
	      GAMMA_J(K,I)=SUM(TC_WRK(1:ID))
	    END DO
	  END DO
	  GAMMA_J=GAMMA_J/2.0_LDP
	  FILENAME='./data/GAMMA_J_'
	  OPTION='FREQ'
	  CALL WRITE_ARRAY_V2(GAMMA_J,ND,NF_GRID_PTS,NU_GRID_VEC,FILENAME,IZERO,OPTION)
	END IF
!
	IF(VERBOSE_GAMMA)THEN
	  GAM_ETA_MUAVG=0D0
	  DO K=1,NF_GRID_PTS
	    DO I=1,ND
	      ID=R_MU(I)%MU_PTS
	      DO J=1,ID-1
	        T1=R_MU(I)%MU_VECTOR(J)
	        T2=R_MU(I)%MU_VECTOR(J+1)
	        GAM_ETA_MUAVG(K,I)=GAM_ETA_MUAVG(K,I)+5.0E-1_LDP*(T1-T2)*( &
			GAM_ETA_SCAT(J,I,K)+GAM_ETA_SCAT(J+1,I,K))
	      END DO
	    END DO
	  END DO
	  FILENAME='./data/ETA_MUAVG_'
	  OPTION='FREQ'
	  CALL WRITE_ARRAY_V2(GAM_ETA_MUAVG,ND,NF_GRID_PTS,NU_GRID_VEC,FILENAME,IZERO,OPTION)
	END IF
!
	CALL TUNE(3,'')
	CLOSE(LU_GAM)
	WRITE(LUER,'(A)')' Finished gamma-ray scattering and deposition - Leaving gamray_sub_v3.f90'
!
	DEALLOCATE (TC_WRK)
	DEALLOCATE (NU_GRID_VEC)
	DEALLOCATE (NU_VEC)
	DEALLOCATE (NU_VEC_15)
	DEALLOCATE (GAMMA_TAU,XRAY_TAU)
	DEALLOCATE (NU_END)
	DEALLOCATE (GAM_E_ME)
	DEALLOCATE (GAM_NU_MAX)
	DEALLOCATE (KLEIN_ARRAY)
	DEALLOCATE (ETA_ISO)
	DEALLOCATE (GAMMA_FLUX)
	DEALLOCATE (DECAY_KIN_E)
	DEALLOCATE (ED_TOT)
	DEALLOCATE (GAMRAY_EMISS)
	DEALLOCATE (TOTAL_DECAY_LUM)
	DEALLOCATE (TA_WRK)
	DEALLOCATE (TB_WRK)
	DEALLOCATE (GAMMA_NU_STORE)
	DEALLOCATE (GAM_I_STORE)
!
	RETURN
	END SUBROUTINE
