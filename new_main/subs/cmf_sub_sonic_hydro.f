	MODULE CMF_SUB_SONIC_HYDRO
	USE SET_KIND_MODULE
!
	  REAL(KIND=LDP) GAM_EDD
	  REAL(KIND=LDP) GAM_LIM
	  REAL(KIND=LDP) GAM_FULL
	  REAL(KIND=LDP) AMU			!Atomic mass unit
	  REAL(KIND=LDP) MU_ATOM		!Mean atomic mass in amu
	  REAL(KIND=LDP) BC			!Boltzmann constant
	  REAL(KIND=LDP) C_CMS
	  REAL(KIND=LDP) SIGMA_TH		!Thompson cross-section
	  REAL(KIND=LDP) STEFAN_BC		!Stephan-Boltzmann constant
	  REAL(KIND=LDP) LOGG			!Log (surface gravity) (cgs units)
	  REAL(KIND=LDP) MDOT
	  REAL(KIND=LDP) GRAV_CON
	  REAL(KIND=LDP) REFERENCE_RADIUS
	  REAL(KIND=LDP) PREV_REF_RADIUS
	  REAL(KIND=LDP) RADIUS_AT_TAU_23
	  REAL(KIND=LDP) SOUND_SPEED
	  REAL(KIND=LDP) MOD_SOUND_SPEED
	  REAL(KIND=LDP) VTURB
!
! The following quantities are used when integrating the hydrostatic equation.
! KAP is used to refer to mass absorption opacities.
!
	  REAL(KIND=LDP) R_EST
	  REAL(KIND=LDP) T_EST
	  REAL(KIND=LDP) P_EST
	  REAL(KIND=LDP) V_EST
	  REAL(KIND=LDP) TAU_EST
	  REAL(KIND=LDP) ED_ON_NA_EST
	  REAL(KIND=LDP) ATOM_EST
	  REAL(KIND=LDP) TEFF
	  REAL(KIND=LDP) ROSS_ON_ES
	  REAL(KIND=LDP) KAP_ROSS
	  REAL(KIND=LDP) KAP_ES
	  REAL(KIND=LDP) F_EST
	  REAL(KIND=LDP) dASQdR_EST
	  REAL(KIND=LDP) ASQ_EST
	  REAL(KIND=LDP) dR
!
	  REAL(KIND=LDP) dVdR
	  REAL(KIND=LDP) dTAUdR
!
	  REAL(KIND=LDP) OLD_TEFF
	  REAL(KIND=LDP) OLD_REF_RADIUS
	  REAL(KIND=LDP) STEP_SIZE
!
	  LOGICAL PURE_LTE_EST
	  LOGICAL OLD_MODEL
	  LOGICAL WIND_PRESENT
	  LOGICAL PLANE_PARALLEL_MOD
	  LOGICAL RESET_REF_RADIUS
	  INTEGER, SAVE :: LAST_HYDRO_ITERATION=0
	  INTEGER DPTH_INDX
!
	  PRIVATE DO_CMF_DERIVS
	  PRIVATE STORE_OLD_GRID
	  PRIVATE CMF_HYDRO_NEW_EST
	  PRIVATE OUT_ESTIMATES
	  PRIVATE CHECK_HYD
!
	  INTEGER  LUV
	  INTEGER  LUV_KUTTA
	  INTEGER  LU_ERR
!
	END MODULE CMF_SUB_SONIC_HYDRO
!
CONTAINS
!
! Subroutine to calculate an updated hydrostatic structure. Routine is to be called from
! CMFGEN.
!
! Subroutine allows for:
!                      (a) A correction for the vdv/dr dynamical term below the sonic point.
!                      (b) A contributin by turbulent pressure.
!
	SUBROUTINE DO_CMF_HYDRO_V3(POPS,MOD_LUM,MOD_TEFF,MOD_LOGG,MOD_MASS,MOD_RSTAR,
	1		  MOD_RMAX,MOD_MDOT,MOD_VINF,MOD_BETA,
	1                 MOD_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NO_V,
	1                 MAIN_COUNTER,DONE_HYDRO,MOD_NC,MOD_ND,MOD_NP,NT)
	USE SET_KIND_MODULE
	USE CMF_SUB_SONIC_HYDRO
	USE OLD_GRID_MODULE
	USE UPDATE_KEYWORD_INTERFACE
	IMPLICIT NONE
!
! Created 10-Jan-2023 (based on do_cmf_hyro_v2.f)
!
	INTEGER MOD_NC			!Model numbr of core rays
	INTEGER MOD_ND			!Model number of depth points
	INTEGER MOD_NP			!Model number of impact papramters
	INTEGER NT			!Number of levels
!
	REAL(KIND=LDP) POPS(NT,MOD_ND)
!
	REAL(KIND=LDP) MOD_MDOT			!Mass loss rate (N=MOD_MDOT/V(kms)/r(10^10cm)^^2)
	REAL(KIND=LDP) MOD_LUM			!In Lsun
	REAL(KIND=LDP) MOD_TEFF			!In unitos of 10^4 K
	REAL(KIND=LDP) MOD_LOGG			!In cgs units
	REAL(KIND=LDP) MOD_MASS			!Mass of star in Msun (returned and output to VADAT)
	REAL(KIND=LDP) MOD_RSTAR		!Returned
	REAL(KIND=LDP) MOD_RMAX			!R(ND) on input
	REAL(KIND=LDP) MOD_VINF			!km/s
	REAL(KIND=LDP) MOD_BETA			!Classic velocity law exponent
	REAL(KIND=LDP) MOD_VTURB
!
	INTEGER MAIN_COUNTER
	LOGICAL PLANE_PARALLEL
	LOGICAL PLANE_PARALLEL_NO_V
	LOGICAL DONE_HYDRO
!
! The following vectors are used for the atmospheric structure resulting
! from the solution of the hydrostatic and tau equations.
!
	INTEGER, PARAMETER :: ND_MAX=4000
        REAL(KIND=LDP) R(ND_MAX)
        REAL(KIND=LDP) V(ND_MAX)
        REAL(KIND=LDP) SIGMA(ND_MAX)
        REAL(KIND=LDP) T(ND_MAX)
        REAL(KIND=LDP) ED(ND_MAX)
        REAL(KIND=LDP) TAU(ND_MAX)
        REAL(KIND=LDP) P(ND_MAX)
        REAL(KIND=LDP) ED_ON_NA(ND_MAX)
!
	REAL(KIND=LDP) dPdR_TERM(ND_MAX)
	REAL(KIND=LDP) GRAV_TERM(ND_MAX)
	REAL(KIND=LDP) RAD_TERM(ND_MAX)
	REAL(KIND=LDP) VdVdR_TERM(ND_MAX)
        REAL(KIND=LDP) ASQ(ND_MAX)
        REAL(KIND=LDP) dASQdR(ND_MAX)
        REAL(KIND=LDP) dASQdR_MON(ND_MAX)
!
        REAL(KIND=LDP) ROSS_MEAN(ND_MAX)
        REAL(KIND=LDP) FLUX_MEAN(ND_MAX)
        REAL(KIND=LDP) POP_ATOM(ND_MAX)
        REAL(KIND=LDP) MASS_DENSITY(ND_MAX)
        REAL(KIND=LDP) CLUMP_FAC(ND_MAX)
        REAL(KIND=LDP) POPION(ND_MAX)
        REAL(KIND=LDP) CHI_ROSS(ND_MAX)
        REAL(KIND=LDP) GAMMA_FULL(ND_MAX)
!
	REAL(KIND=LDP) TA(ND_MAX)
	REAL(KIND=LDP) TB(ND_MAX)
	REAL(KIND=LDP) TC(ND_MAX)
!
! These vectors are output, and contain the atmospheric structure to be used
! by CMFGEN.
!
	INTEGER NEW_ND
	REAL(KIND=LDP), ALLOCATABLE :: REV_TAU(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_R(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_V(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_SIGMA(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_POP_ATOM(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_ED(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_T(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_CHI_ROSS(:)
        REAL(KIND=LDP), ALLOCATABLE :: REV_GAMMA_FULL(:)
        REAL(KIND=LDP), ALLOCATABLE :: COEF(:,:)
!
! Parameters for for cumputing the final R grid.
!
	REAL(KIND=LDP) dLOG_TAU
	REAL(KIND=LDP) V_SCL_FAC
	REAL(KIND=LDP) OBND_PARS(20)
	INTEGER NUM_OBND_PARAMS
	CHARACTER(LEN=16) OUT_BND_OPT
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=20) HYDRO_OPT
!
! Wind parameters:
!
	REAL(KIND=LDP) VINF
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) BETA2
	REAL(KIND=LDP) RMAX
	REAL(KIND=LDP) CONNECTION_VEL
	REAL(KIND=LDP) CONNECTION_RADIUS
	REAL(KIND=LDP) OLD_CONNECTION_RADIUS
	REAL(KIND=LDP) RP2_ON_CON_RAD
	REAL(KIND=LDP) VEXT
	INTEGER CONNECTION_INDX
!
	INTEGER ND			!Number of points in RUnge-Kutta integration
	INTEGER I,J,IOS,IST
!
	REAL(KIND=LDP) PI			!
	REAL(KIND=LDP) NI_ZERO			!Density at outer boundary
	REAL(KIND=LDP) SCL_HT			!Atmopsheric scale height
	REAL(KIND=LDP) H			!Step size for Runge-Kutta integration
	REAL(KIND=LDP) PTURB_ON_NA		!Turbulent pressure due to turbulence/ (# of atoms)
	REAL(KIND=LDP) MASS_LOSS_SCALE_FACTOR	!Factor to convert MOD_MDOt to Msun/yr
	REAL(KIND=LDP) RBOUND
	REAL(KIND=LDP) TAU_MAX
	REAL(KIND=LDP) OLD_TAU_MAX
	REAL(KIND=LDP) BOLD,BNEW
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) GAM_LIM_STORE
	REAL(KIND=LDP) VC_ON_SS
	REAL(KIND=LDP) TAU_REF
	REAL(KIND=LDP) MAX_ED_ON_NA
	REAL(KIND=LDP) HYDRO_SUM
!
! Runge-Kutta estimates
!
	REAL(KIND=LDP) dV1,dTAU1
	REAL(KIND=LDP) dV2,dTAU2
	REAL(KIND=LDP) dV3,dTAU3
	REAL(KIND=LDP) dV4,dTAU4
!
	LOGICAL USE_OLD_VEL
	LOGICAL L_TEMP
	LOGICAL FILE_OPEN
	LOGICAL FILE_PRES
	LOGICAL VERBOSE_OUTPUT
	LOGICAL UPDATE_GREY_SCL
!
! These have cgs units.
!
	REAL(KIND=LDP) ATOMIC_MASS_UNIT,BOLTZMANN_CONSTANT,STEFAN_BOLTZ
	REAL(KIND=LDP) GRAVITATIONAL_CONSTANT,MASS_SUN,SPEED_OF_LIGHT,LUM_SUN
	INTEGER GET_INDX_DP,ERROR_LU
	EXTERNAL ATOMIC_MASS_UNIT,BOLTZMANN_CONSTANT,STEFAN_BOLTZ,GET_INDX_DP
	EXTERNAL GRAVITATIONAL_CONSTANT,MASS_SUN,SPEED_OF_LIGHT,ERROR_LU,LUM_SUN
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUSCR=8
	INTEGER  STRT_HYDRO_ITS
	INTEGER  FREQ_HYDRO_ITS
	INTEGER  NO_HYDRO_ITS
	INTEGER  NO_ITS_DONE
	INTEGER  ITERATION_COUNTER
	INTEGER  VEL_LAW
!
	INTEGER  LU
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	CHARACTER*2, PARAMETER :: FORMFEED=' '//CHAR(12)
!
! Set constants.
!
	LU_ERR=ERROR_LU()
	BC=1.0E+04_LDP*BOLTZMANN_CONSTANT()   			!erg/10^4 K
	AMU=ATOMIC_MASS_UNIT()          			!gm
	C_CMS=SPEED_OF_LIGHT()   	       			!cm/s
	GRAV_CON=1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*MASS_SUN()
	SIGMA_TH=6.65E-15_LDP					!cm^{-2} x 10^10
	STEFAN_BC=STEFAN_BOLTZ()
	MASS_LOSS_SCALE_FACTOR=3.02286E+23_LDP
	PI=ACOS(-1.0_LDP)
	STRT_HYDRO_ITS=20
	FREQ_HYDRO_ITS=8
	HYDRO_OPT='DEFAULT'
!
	CALL STORE_OLD_GRID(MU_ATOM,MOD_ND)
!
	TEFF=MOD_TEFF
	LOGG=MOD_LOGG
	VINF=MOD_VINF
	BETA=MOD_BETA
	MDOT=MOD_MDOT
	RMAX=MOD_RMAX/OLD_R(OLD_ND)
	VTURB=MOD_VTURB
	VEL_LAW=ITWO
	BETA2=0.0_LDP
	VEXT=-1.0_LDP
	RP2_ON_CON_RAD=-1.0_LDP
!
	VC_ON_SS=0.75_LDP
	TAU_REF=2.0_LDP/3.0_LDP
        dLOG_TAU=0.25_LDP
        V_SCL_FAC=0.75E00_LDP
        OBND_PARS(:)=0.0_LDP
        NUM_OBND_PARAMS=1
        OUT_BND_OPT='DEFAULT'
	UPDATE_GREY_SCL=.FALSE.
!
! Set defaults:
!
	IF(PLANE_PARALLEL)THEN
	  NI_ZERO=1.0E+06_LDP
	  PLANE_PARALLEL_MOD=.TRUE.
	  WIND_PRESENT=.TRUE.
	  RESET_REF_RADIUS=.FALSE.
	ELSE IF(PLANE_PARALLEL_NO_V)THEN
	  NI_ZERO=1.0E+06_LDP
	  PLANE_PARALLEL_MOD=.TRUE.
	  WIND_PRESENT=.FALSE.
	  RESET_REF_RADIUS=.FALSE.
	ELSE
	  PLANE_PARALLEL_MOD=.FALSE.
	  WIND_PRESENT=.TRUE.
	  RESET_REF_RADIUS=.TRUE.
	END IF
	USE_OLD_VEL=.FALSE.
	GAM_LIM=0.98_LDP
	GAM_LIM_STORE=GAM_LIM
!
! These are the default settings if no old model.
!
	PURE_LTE_EST=.FALSE.
	OLD_TAU_MAX=100.0_LDP
	RBOUND=0.0_LDP
	VERBOSE_OUTPUT=.TRUE.
!
! *************************************************************************
!
! Read in parameters describing the new model.
!
	CALL GEN_ASCI_OPEN(LUIN,'HYDRO_DEFAULTS','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LU_ERR,*)'Error opening HYDRO_DEFAULTS in WIND_HYD, IOS=',IOS
	  RETURN
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)

        CALL RD_STORE_INT(NO_HYDRO_ITS,'N_ITS',L_TRUE,'Number of hydro iterations remaining')
        CALL RD_STORE_INT(NO_ITS_DONE,'ITS_DONE',L_TRUE,'Number of hydro iterations completed')
        CALL RD_STORE_INT(STRT_HYDRO_ITS,'STRT_ITS',L_FALSE,'Iteration to start first hydro iteration')
        CALL RD_STORE_INT(FREQ_HYDRO_ITS,'FREQ_ITS',L_FALSE,'Frequency for hydro iterations')
	CALL RD_STORE_DBLE(NI_ZERO,'ATOM_DEN',L_FALSE,'Atom density at outer boundary (/cm^3)')
	CALL RD_STORE_LOG(USE_OLD_VEL,'OLD_V',L_FALSE,'Use old velocity law above connection point?')
	CALL RD_STORE_DBLE(RMAX,'MAX_R',L_FALSE,'Maximum radius in terms of Connection radius')
	CALL RD_STORE_LOG(RESET_REF_RADIUS,'RES_REF',L_FALSE,'Reset reference radius if using old velocity law')
	CALL RD_STORE_DBLE(GAM_LIM,'GAM_LIM',L_FALSE,'Limiting Eddington factor')
	GAM_LIM_STORE=GAM_LIM
	CALL RD_STORE_DBLE(VC_ON_SS,'VC_ON_SS',L_FALSE,'Connection velocity on sound speed')
	CALL RD_STORE_LOG(UPDATE_GREY_SCL,'UP_GREY_SCL',L_FALSE,'Update GREY_SCL_FAC_IN')
	CALL RD_STORE_DBLE(TAU_REF,'TAU_REF',L_FALSE,'Reference radius for g and Teff')
	CALL RD_STORE_INT(VEL_LAW,'VEL_LAW',L_FALSE,'Velocity law for wind region (2 or 3)')
	BETA2=BETA
	CALL RD_STORE_DBLE(BETA2,'BETA2',L_FALSE,'Second exponent for velocity law')
	CALL RD_STORE_DBLE(VEXT,'VEXT',L_FALSE,'Extension (km/s) to VINFr. Note: VINF1=VINF-VEXT')
	CALL RD_STORE_DBLE(RP2_ON_CON_RAD,'RP2_ON_RT',L_FALSE,'Ration of RP2 to transition radius')
	IF(VEL_LAW .EQ. 5)THEN
	  IF(BETA2 .EQ. 0.0_LDP)THEN
	    WRITE(LU_ERR,*)'Error in DO_CMF_HYDRO_V3 -- BETA2 cannot be zero for VEL_LAW 5'
	    STOP
	  END IF
	  IF(VEXT .LT. 0.0_LDP)THEN
	    WRITE(LU_ERR,*)'Error in DO_CMF_HYDRO_V3 -- VEXT cannot be less than zero for VEL_LAW 5'
	    STOP
	  END IF
	  IF(RP2_ON_CON_RAD .LT. 0.0_LDP)THEN
	    WRITE(LU_ERR,*)'Error in DO_CMF_HYDRO_V3 -- RP2_ON_CON_RAD cannot be less than zero for VEL_LAW 5'
	    STOP
	  END IF
	END IF
!
! Therse are the parameters used to define the new R grid to be output to RVSIG_COL.
!
	CALL RD_STORE_DBLE(dLOG_TAU,'dLOG_TAU',L_FALSE,'Logarithmic spacing in Tau for new R grid')
	CALL RD_STORE_DBLE(V_SCL_FAC,'VSCL_FAC',L_FALSE,'Maximum V(I-1)/V(I) for new R grid (<1)')
	I=10
	CALL RD_STORE_NCHAR(OUT_BND_OPT,'OB_OPT',I,L_FALSE,'Outer boundary option: POW, SPECIFY, DEFAULT, NONE')
	J=0; CALL RD_STORE_INT(J,'NOB_PARS',L_FALSE,'Number of outer boudary parameters')
	DO I=1,J
	  NUM_OBND_PARAMS=J
	  WRITE(STRING,'(I3)')I
	  STRING='OB_P'//ADJUSTL(STRING)
	  CALL RD_STORE_DBLE(OBND_PARS(I),TRIM(STRING),L_TRUE,'Paremeters for outer boundary condition')
	END DO
	CALL RD_STORE_CHAR(HYDRO_OPT,'HYDRO_OPT',L_FALSE,'FIXED_R or FIXED_V_FLUX or DEFAULT')
	IF(HYDRO_OPT .EQ. 'FIXED_V_FLUX')THEN
	  CALL RD_STORE_DBLE(OLD_TEFF,'OLD_TEFF',L_TRUE,'Effective temperatre of input model')
	END IF
	IF(TAU_REF .GT. 0.668_LDP .AND. HYDRO_OPT .EQ. 'DEFAULT')THEN
          HYDRO_OPT='FIXED_R_REF'
          WRITE(LU_ERR,*)'As TAU_REF > 2/3 HYDRO_DEFAULT is being set to FIXED_R_REF in DO_CMF_HYDRO_V2'
        END IF
	CALL CLEAN_RD_STORE()
!
        CLOSE(UNIT=LUIN)
        CLOSE(UNIT=LUSCR)
!
! Decide here whether we will do an iteration or not.
!
	DONE_HYDRO=.FALSE.
	IF(NO_HYDRO_ITS .EQ. 0)RETURN
	IF(MAIN_COUNTER .LT. STRT_HYDRO_ITS)RETURN
	IF( MOD( (MAIN_COUNTER-STRT_HYDRO_ITS),FREQ_HYDRO_ITS ) .NE. 0)RETURN
!
! Begin Hydro computation.
!
	IF(VERBOSE_OUTPUT)THEN
	  CALL GET_LU(LUV)
	  OPEN(UNIT=LUV,FILE='HYDRO_ITERATION_INFO',STATUS='UNKNOWN')
	  CALL SET_LINE_BUFFERING(LUV)
	  CALL GET_LU(LUV_KUTTA)
	  OPEN(UNIT=LUV_KUTTA,FILE='HYDRO_ITERATION_KUTTA_INFO',STATUS='UNKNOWN')
  	END IF
	CALL GET_LU(LU)				!For files open/shut immediately
!
! In TORSCL_V3, TA is TAU, TB is dTAU, and TC is used fro dCHIdR.
!
	WRITE(LU_ERR,'(/,A)')' Updating hydrostatic structure of the model'
	FLUSH(UNIT=6)
!
	IF(HYDRO_OPT .EQ. 'FIXED_R_REF')THEN
!
! This option keeps the radius, at a pre-specified TAU, fixed.
! To preserve the specifed effective temperature, the luminosity is
! updated.
!
	  WRITE(LU_ERR,'(A)')' Using FIXED_R_REF option in DO_CMF_HYDRO_V3'
	  CHI_ROSS(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	  I=7			!Use CHI(1) and CHI(I) to computed exponent.
	  CALL TORSCL_V3(TA,CHI_ROSS,OLD_R,TB,TC,MOD_ND,'LOGMON','PCOMP',I,L_FALSE)
	  DO I=1,MOD_ND
	    IF(TA(I) .GT. TAU_REF)THEN
	      T1=(TAU_REF-TA(I-1))/(TA(I)-TA(I-1))
	      REFERENCE_RADIUS=T1*OLD_R(I)+(1.0_LDP-T1)*OLD_R(I-1)
	      EXIT
	    END IF
	  END DO
!
	  T1=REFERENCE_RADIUS*TEFF*TEFF
	  MOD_LUM=4.0E+36_LDP*PI*STEFAN_BC*T1*T1/LUM_SUN()
	  CALL UPDATE_KEYWORD(MOD_LUM,'[LSTAR]','VADAT',L_TRUE,L_TRUE,LUIN)
	  CALL UPDATE_KEYWORD('DEFAULT','[HYDRO_OPT]','HYDRO_DEFAULTS',L_TRUE,L_TRUE,LUIN)
	  WRITE(LU_ERR,'(A)')' DO_CMF_HYDRO_V3 has adjusted LSTAR in VADAT'
!
! This option is useful for WR models where the key variable controlling the observed
! spectrum is the luminosity.
!
	ELSE IF(HYDRO_OPT .EQ. 'FIXED_LUM')THEN
	  WRITE(LU_ERR,'(A)')' Using FIXED_LUM option in DO_CMF_HYDRO_V3'
	  CHI_ROSS(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	  I=7			!Use CHI(1) and CHI(I) to computed exponent.
	  CALL TORSCL_V3(TA,CHI_ROSS,OLD_R,TB,TC,MOD_ND,'LOGMON','PCOMP',I,L_FALSE)
	  DO I=1,MOD_ND
	    IF(TA(I) .GT. TAU_REF)THEN
	      T1=(TAU_REF-TA(I-1))/(TA(I)-TA(I-1))
	      REFERENCE_RADIUS=T1*OLD_R(I)+(1.0_LDP-T1)*OLD_R(I-1)
	      EXIT
	    END IF
	  END DO
!
! Compute the revised effective temperature.
!
	  T1=4.0E+36_LDP*PI*STEFAN_BC*REFERENCE_RADIUS*REFERENCE_RADIUS
	  TEFF=(MOD_LUM*LUM_SUN()/T1)**0.25_LDP
	  CALL UPDATE_KEYWORD(TEFF,'[TEFF]','VADAT',L_TRUE,L_TRUE,LUIN)
	  WRITE(LU_ERR,'(A,ES14.4)')' DO_CMF_HYDRO_V3 has adjusted TEFF in VADAT: Teff=',TEFF
!
! This option (useful for O stars) attempts to preserve the V =-band flux.
!
	ELSE IF(HYDRO_OPT .EQ. 'FIXED_V_FLUX')THEN
	  WRITE(LU_ERR,'(A)')' Using FIXED_V_FLUX option in DO_CMF_HYDRO_V3'
	  I=7	!Use CHI(1) and CHI(I) to computed exponent.
	  CHI_ROSS(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	  CALL TORSCL_V3(TA,CHI_ROSS,OLD_R,TB,TC,MOD_ND,'LOGMON','PCOMP',I,L_FALSE)
	  DO I=1,MOD_ND
	    IF(TA(I) .GT. TAU_REF)THEN
	      T1=(TAU_REF-TA(I-1))/(TA(I)-TA(I-1))
	      REFERENCE_RADIUS=T1*OLD_R(I)+(1.0_LDP-T1)*OLD_R(I-1)
	      EXIT
	    END IF
	  END DO
          BOLD=1.0_LDP/(EXP(1.4388_LDP/0.55_LDP/OLD_TEFF)-1.0_LDP)
          BNEW=1.0_LDP/(EXP(1.4388_LDP/0.55_LDP/TEFF)-1.0_LDP)
	  REFERENCE_RADIUS=REFERENCE_RADIUS*SQRT(BOLD/BNEW)
	  T1=REFERENCE_RADIUS*TEFF*TEFF
	  MOD_LUM=4.0E+36_LDP*PI*STEFAN_BC*T1*T1/LUM_SUN()
	  CALL UPDATE_KEYWORD(MOD_LUM,'[LSTAR]','VADAT',L_TRUE,L_TRUE,LUIN)
	  CALL UPDATE_KEYWORD('DEFAULT','[HYDRO_OPT]','HYDRO_DEFAULTS',L_TRUE,L_FALSE,LUIN)
	  CALL UPDATE_KEYWORD(TEFF,'[OLD_TEFF]','HYDRO_DEFAULTS',L_FALSE,L_TRUE,LUIN)
	  WRITE(LU_ERR,'(A)')' DO_CMF_HYDRO_V3 has adjusted LSTAR in VADAT'
!
	ELSE IF(HYDRO_OPT .EQ. 'DEFAULT')THEN
	  IF( ABS(TAU_REF-2.0_LDP/3.0_LDP) .GT. 0.001_LDP)THEN
	    WRITE(LU_ERR,*)'Error -- for the DEFAULT HYDRO_OPTION, TAU_REF must be 2/3'
	    STOP
	  END IF
	  WRITE(LU_ERR,'(A)')' Reference radius: based on effective temperature and luminosity of star'
	  REFERENCE_RADIUS=1.0E-18_LDP*SQRT(MOD_LUM*LUM_SUN()/TEFF**4/STEFAN_BC/4.0_LDP/PI)
!
	ELSE
	  WRITE(LU_ERR,'(A)')' Error in do_cmf_hydro_v3.f: invlaid HYDRO_OPT option'
	  WRITE(LU_ERR,'(2A)')' HYDRO_OPT=',TRIM(HYDRO_OPT)
	  STOP
	END IF
	WRITE(LU_ERR,*)'Reference radius is',REFERENCE_RADIUS; FLUSH(UNIT=6)
!
!
!
! Compute the grey temperature distribution (returned in TC) and the Rosseland
! optical depth scale (returned in TA).
!
	TB(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	CALL COMP_GREY_V2(POPS,TC,TA,TB,LU_ERR,MOD_NC,MOD_ND,MOD_NP,NT)
	IF(MINVAL(OLD_ROSS_MEAN(1:MOD_ND)) .LE. 0.0_LDP)THEN
	   WRITE(LU_ERR,'(A)')' Bad Rosseland optical depth scale for T/TGREY output'
	   WRITE(LU_ERR,'(A)')' GREY_SCL_FAC_OUT not output'
	ELSE IF(UPDATE_GREY_SCL)THEN
	  OPEN(UNIT=LUIN,FILE='GREY_SCL_FAC_IN',STATUS='UNKNOWN')
	    WRITE(LUIN,'(A)')'!'
	    WRITE(LUIN,'(A,8X,A,7X,A,7X,A,6X,A)')'!','Log(Tau)','T/T(grey)','T(10^4 K)','L'
	    WRITE(LUIN,'(A)')'!'
	    WRITE(LUIN,*)MOD_ND
	    DO I=1,MOD_ND
	      WRITE(LUIN,'(2X,3ES16.6,4X,I3)')LOG10(TA(I)),OLD_T(I)/TC(I),OLD_T(I),I
	    END DO
	  CLOSE(LUIN)
	END IF
!
! Compute the Rosseland optical depth scale.
!
	TB(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	CALL TORSCL(OLD_TAU,TB,OLD_R,TA,TC,MOD_ND,'LOGMON','EXP')
	IF(PLANE_PARALLEL_MOD)THEN
	  OLD_REF_RADIUS=OLD_R(OLD_ND)
	ELSE
	  T1=TAU_REF
	  I=GET_INDX_DP(T1,OLD_TAU,MOD_ND)
	  T2=(LOG(T1)-LOG(OLD_TAU(I)))/(LOG(OLD_TAU(I+1))-LOG(OLD_TAU(I)))
	  OLD_REF_RADIUS=(1.0_LDP-T1)*OLD_R(I)+T1*OLD_R(I+1)
	  IF(VERBOSE_OUTPUT)THEN
	    WRITE(LUV,*)'OLD_TAU(1)=',OLD_TAU(1)
	    WRITE(LUV,*)'Reference radius (Tau=',TAU_REF,') of of old model is',OLD_REF_RADIUS
	  END IF
	END IF
!
	T1=OLD_REF_RADIUS/6.9599_LDP
	OLD_TEFF=TEFF
	OLD_SF(1:MOD_ND)=OLD_T(1:MOD_ND)/OLD_TEFF/(OLD_TAU(1:MOD_ND)+0.67_LDP)**0.25_LDP
	OLD_TAU_MAX=OLD_TAU(MOD_ND)
!
	IF(VERBOSE_OUTPUT)THEN
	  OPEN(UNIT=LU,FILE='HYDRO_OLD_MODEL',STATUS='UNKNOWN')
	  WRITE(LU,'(/,A,8(7X,A))')' Index','     R','   Tau','    Na',' Ne/Na','     T',
	1                  ' Kross','Kr/Kes','Kf/Kes'
	  DO I=1,MOD_ND
	    WRITE(LU,'(I6,8ES13.4)')I,OLD_R(I),OLD_TAU(I),OLD_POP_ATOM(I),OLD_ED(I)/OLD_POP_ATOM(I),
	1        OLD_T(I),OLD_KAP_ROSS(I),OLD_ROSS_MEAN(I)/OLD_ESEC(I),OLD_FLUX_MEAN(I)/OLD_ESEC(I)
	  END DO
!
	  WRITE(LU,'(/,A,/,A)')FORMFEED,' Old model mass absorption coefficients'
	  WRITE(LU,'(A,3(8X,A))')' Index',' Kross',' Kflux','  Kes'
	  DO I=1,MOD_ND
	    WRITE(LU,'(I6,3ES14.5)')I,OLD_KAP_ROSS(I),OLD_KAP_FLUX(I),OLD_KAP_ESEC(I)
	  END DO
	  CLOSE(UNIT=LU)
!
	  WRITE(LUV,*)' '
	  WRITE(LUV,'(A,ES14.4)')' Old TEFF is',OLD_TEFF
	  WRITE(LUV,'(A,ES14.4)')' Optical depth at inner boundary (om):',OLD_TAU(MOD_ND)
	  WRITE(LUV,*)' '
	END IF
!
	IF(WIND_PRESENT)THEN
	  T1=1.0E-10_LDP*BC/MU_ATOM/AMU
	  DO I=1,MOD_ND
	    ED_ON_NA_EST=OLD_ED(I)/OLD_POP_ATOM(I)
	    MOD_SOUND_SPEED=SQRT(T1*(1.0_LDP+ED_ON_NA_EST)*OLD_T(I)+0.5_LDP*VTURB*VTURB)
	    SOUND_SPEED=SQRT(T1*(1.0_LDP+ED_ON_NA_EST)*OLD_T(I))
	    IF(OLD_V(I) .LT. VC_ON_SS*MOD_SOUND_SPEED)THEN
	      CONNECTION_VEL=OLD_V(I)
	      CONNECTION_RADIUS=OLD_R(I)
	      CONNECTION_INDX=I
	      RMAX=RMAX*CONNECTION_RADIUS
	      WRITE(LU_ERR,*)'  Connection velocity is',CONNECTION_VEL
	      WRITE(LU_ERR,*)'    Connection radius is',CONNECTION_RADIUS
	      WRITE(LU_ERR,*)'       Maximum radius is',RMAX
	      WRITE(LU_ERR,*)'     Connection INDEX is',CONNECTION_INDX
	      WRITE(LU_ERR,*)'          Sound speed is',SOUND_SPEED
	      WRITE(LU_ERR,*)' Modified sound speed is',MOD_SOUND_SPEED
	      EXIT
	    END IF
	  END DO
	  OLD_CONNECTION_RADIUS=CONNECTION_RADIUS
	END IF
!
!
!
! Compute Eddington ratio, GAM_EDD. This formulae is set for one electron per ion.
! This formula holds at all radii, since g and Teff both scale as 1/r^2.
!
! Note: The factor of 10^6 in the expression occurs because:
!          (a) Teff is units of 10^4 K ==> 10^{16}
!          (b) SIGMA_TH * R is unitless, and R is in units of 10^10 cm.
!
	MAX_ED_ON_NA=0.0_LDP
	DO I=1,MOD_ND
	  MAX_ED_ON_NA=MAX(MAX_ED_ON_NA,OLD_ED(I)/OLD_POP_ATOM(I))
	END DO
	GAM_EDD=1.0E+06_LDP*SIGMA_TH*STEFAN_BC*(TEFF**4)/MU_ATOM/C_CMS/(10**LOGG)/AMU
	IF(GAM_EDD*MAX_ED_ON_NA .GT. 1.0_LDP)THEN
	  WRITE(LU_ERR,*)'An invalid Eddington parameter has been computed in DO_CMF_HYDRO_V3'
	  WRITE(LU_ERR,*)'Check the validity of Teff and Log G'
	  WRITE(LU_ERR,*)'The computed (maximum) Eddington parameter is ',GAM_EDD*MAX_ED_ON_NA
	  WRITE(LU_ERR,*)'                            Teff(K)/1.0+04 is',TEFF
	  WRITE(LU_ERR,*)'                         LOG_G (cgs units) is',LOGG
	  WRITE(LU_ERR,*)'                                MAX(Ne/Na) is',MAX_ED_ON_NA
	  STOP
	END IF
!
	IF(VERBOSE_OUTPUT)THEN
	  WRITE(LUV,'(A,ES14.6)')'          Surface gravity is:',LOGG
	  WRITE(LUV,'(A,ES14.6)')'             Mass of star is:',10**(LOGG)*(REFERENCE_RADIUS**2)/GRAV_CON
	  WRITE(LUV,'(A,ES14.6)')'         Mean atomic mass is:',MU_ATOM
	  WRITE(LUV,'(A,ES14.6)')'             Atom density is:',NI_ZERO
	  WRITE(LUV,'(A,ES14.6)')'New effective temperature is:',TEFF
	  WRITE(LUV,'(A,ES14.6,A)')'      Eddington parameter is:',GAM_EDD,'  (assuming Ne/NA=1)'
	  WRITE(LUV,'(A,ES14.6)')'               MAX(Ne/Na) is:',MAX_ED_ON_NA
	END IF
!
	PREV_REF_RADIUS=-1.0_LDP
	ITERATION_COUNTER=0
	DO WHILE(ABS(REFERENCE_RADIUS/PREV_REF_RADIUS-1.0_LDP) .GT. 1.0E-05_LDP)
	  ITERATION_COUNTER=ITERATION_COUNTER+1
	  IF(VERBOSE_OUTPUT)WRITE(LUV,*)' Beginning new hydro loop'
!
! The turbulent pressure is taken to be 0.5. roh . VTURB^2
!
	  PTURB_ON_NA=0.5E+10_LDP*VTURB*VTURB*MU_ATOM*AMU
!
!
! Set parameters/initial conditions at the outer boundary of the
! hydrostatic structure.
!
! The first section assumes we have a wind present. Given this we have 3 choices:
!   (1) We have an old model and will use its wind
!   (2) We have an old model, but will input a new wind.
!   (3) We don't have an old model.
!
	  WRITE(LU_ERR,*)'WIND_PRESENT=',WIND_PRESENT
	  IF(WIND_PRESENT)THEN
!
! In this case will use exactly the same grid as for the old model beyond
! the connection velocity (i.e., at larger V in the wind).
!
! NB: 1.0D-06 is  1.0E+04/(1.0E+05)**2. The first factor occurs since T is in
! units of 10^4K. The second is to convert V from cm/s to km/s, with allowance
! for the sqrt.
!
	      I=CONNECTION_INDX
	      ED_ON_NA_EST=OLD_ED(I)/OLD_POP_ATOM(I)
	      ROSS_ON_ES=OLD_ROSS_MEAN(I)/OLD_ESEC(I)
	      T2=-1.0_LDP
	      DO WHILE(T2 .LT. 0.0_LDP)
	        GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA_EST*(OLD_FLUX_MEAN(I)/OLD_ESEC(I)) )
	        T1=1.0E-10_LDP*(BC*TEFF*(1+ED_ON_NA_EST)+PTURB_ON_NA)/( (10.0_LDP**LOGG)*(1.0_LDP-GAM_FULL)*MU_ATOM*AMU )
	        T2=CONNECTION_VEL*(1.0_LDP/T1-2.0_LDP/CONNECTION_RADIUS)
	        WRITE(LU_ERR,*)'      Scale height is',T1
	        WRITE(LU_ERR,*)'      Connection dVdR',T2
	        WRITE(LU_ERR,*)'       Ne/Na estimate',ED_ON_NA_EST
	        WRITE(LU_ERR,*)'        Flux Op./esec',OLD_FLUX_MEAN(I)/OLD_ESEC(I)
                IF(T2 .LE. 0.0_LDP)THEN
	          WRITE(LU_ERR,*)'    Connection radius',CONNECTION_RADIUS
                  WRITE(LU_ERR,*)'Resetting GAM_LIM due to -ve velocity gradient'
                  WRITE(LU_ERR,*)'         GAM_FULL was',GAM_FULL
                  WRITE(LU_ERR,*)'          GAM_LIM was',GAM_LIM
	          GAM_LIM=GAM_LIM-0.01_LDP
                  WRITE(LU_ERR,*)'       New GAM_LIM is',GAM_LIM
	        END IF
	      END DO
	      FLUSH(UNIT=LU_ERR)
!
	      IF(USE_OLD_VEL .AND. HYDRO_OPT .EQ. 'FIXED_LUM')THEN !RESET_REF_RADIUS
	        J=CONNECTION_INDX
	        R(1:J)=OLD_R(1:J)
	        V(1:J)=OLD_V(1:J)
	        SIGMA(1:J)=OLD_SIGMA(1:J)
!	        CLUMP_FAC(1:J)=OLD_CLUMP_FAC(1:J)
	      ELSE IF(USE_OLD_VEL)THEN
	        J=CONNECTION_INDX
	        R(1:J)=OLD_R(1:J)*CONNECTION_RADIUS/OLD_CONNECTION_RADIUS
	        V(1:J)=OLD_V(1:J)
	        SIGMA(1:J)=OLD_SIGMA(1:J)
	      ELSE
	        FLUSH(UNIT=6)
	        CALL WIND_VEL_LAW_V3(R,V,SIGMA,VINF,BETA,BETA2,
	1            VEXT,RP2_ON_CON_RAD,RMAX,
	1            CONNECTION_RADIUS,CONNECTION_VEL,T2,VEL_LAW,J,ND_MAX)
	        WRITE(LUV,'(2X,A,ES15.6)')'Connection temperature is',OLD_T(J)
	     END IF
!
	     DO I=1,J
	       POP_ATOM(I)=MDOT/MU_ATOM/R(I)/R(I)/V(I)    !/CLUMP_FAC(I)
	       MASS_DENSITY(I)=MU_ATOM*AMU*POP_ATOM(I)
	     END DO
!
! To get other quantities we interpolate as a function of density.
! The atom density should be monotonic. At the outer boundary,
! we simply use the boundary value.
!
! POP_ATOM is not corrected for clumping.
!
	      TB(1:MOD_ND)=OLD_ED(1:MOD_ND)/ OLD_POP_ATOM(1:MOD_ND)
	      TA(1:MOD_ND)=LOG(OLD_POP_ATOM(1:MOD_ND)*OLD_CLUMP_FAC(1:MOD_ND))
	      TC(1:J)=LOG(POP_ATOM(1:J))
	      DO I=1,J
	        IF(TC(I) .LT. TA(1))TC(I)=TA(1)
	      END DO
	      CALL MON_INTERP(ED_ON_NA,J,IONE,TC,J,TB,MOD_ND,TA,MOD_ND)
	      TB(1:MOD_ND)=OLD_T(1:MOD_ND)*TEFF/OLD_TEFF
	      CALL MON_INTERP(T,J,IONE,TC,J,TB,MOD_ND,TA,MOD_ND)
	      TB(1:MOD_ND)=OLD_ROSS_MEAN(1:MOD_ND)/OLD_ESEC(1:MOD_ND)
	      CALL MON_INTERP(CHI_ROSS,J,IONE,TC,J,TB,MOD_ND,TA,MOD_ND)
	      DO I=1,J
	        ED(I)=ED_ON_NA(I)*POP_ATOM(I)
	        GAMMA_FULL(I)=GAM_FULL
	        CHI_ROSS(I)=ED_ON_NA(I)*SIGMA_TH*POP_ATOM(I)*CHI_ROSS(I)
	        P(I)=(BC*T(I)*(1.0_LDP+ED(I)/POP_ATOM(I))+PTURB_ON_NA)*POP_ATOM(I)
	        ASQ(I)=1.0E-10_LDP*BC*T(I)*(1+ED_ON_NA(I))/MU_ATOM/AMU               !km/s
	        IF(I .NE. 1)dASQdR(I)=(ASQ(I)-ASQ(I-1))/(R(I)-R(I-1))
	      END DO
	      ROSS_ON_ES=CHI_ROSS(J)/ED(J)/SIGMA_TH
	      dASQdR(1)=dASQdR(2)
	      CALL TORSCL(TAU,CHI_ROSS,R,TB,TC,J,'LOGMON',' ')
	      WRITE(LUV,'(10X,A,I7)')'Connection index is',J
	      WRITE(LUV,'(10X,A,ES15.6)')'Connection tau is',TAU(J)
	      I=J
!
! The following section is for the case of no wind.
!
	  ELSE
!
! This option is for a pure spherical or plane-parallel model. The depth variation
! of gravity is taken into account.
!
	      POP_ATOM(1)=NI_ZERO
	      T(1)=(TEFF/OLD_TEFF)*OLD_T(1)
	      P(1)=(BC*T(1)*(1.0_LDP+OLD_ED(1)/OLD_POP_ATOM(1))+PTURB_ON_NA)*POP_ATOM(1)
	      IF(RBOUND .EQ. 0.0_LDP)RBOUND=OLD_R(1)
	      R(1)=RBOUND
	      TAU(1)=OLD_TAU(1)
	      ED_ON_NA(1)=OLD_ED(1)/OLD_POP_ATOM(1)
	      ROSS_ON_ES=OLD_ROSS_MEAN(1)/OLD_ESEC(1)
	      GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA(1)*(OLD_FLUX_MEAN(1)/OLD_ESEC(1)) )
	      ED(1)=ED_ON_NA(1)*POP_ATOM(1)
	      CHI_ROSS(1)=ED_ON_NA(1)*SIGMA_TH*POP_ATOM(1)*ROSS_ON_ES
	      GAMMA_FULL(1)=GAM_FULL
	      V(1)=MDOT/MU_ATOM/POP_ATOM(1)/R(1)/R(1)
	      I=1
	  END IF
!
! Set sound speed at the point we match the wind to the hydrostatic structure.
! In no wind is present, we set the sound speed to a large number. This effectively
! set the correction to zero.
!
	  IF(WIND_PRESENT)THEN
	    SOUND_SPEED=(1.0_LDP+ED_ON_NA(I))*BC*T(I)/MU_ATOM/AMU
	    SOUND_SPEED=1.0E-05_LDP*SQRT(SOUND_SPEED)
	    WRITE(LU_ERR,'(A,3ES14.4)')'SOUND_SPEED',SOUND_SPEED
	  ELSE
	    IF(VTURB .EQ. 0.0_LDP)THEN
	      SOUND_SPEED=1.0E+30_LDP
	    ELSE
	      SOUND_SPEED=0.0E+00_LDP
	    END IF
	  END IF
!
! Compute the atmospheric pressure scale height. Close to the sonic point,
! Gamma can be very close to 1 which leads to a large scale height. But this
! will be come smaller as we move away from the sonic point. Thus to ensure
! that the step size is sufficiently small, we set the scale height using
! GAM_EDD (i.e., GAMMA computed the electron scattering opacity only) rather
! than GAM_FULL.
!
	    SCL_HT=(10.0_LDP**LOGG)*(1.0_LDP-GAM_EDD)*MU_ATOM*AMU/
	1            (BC*(1.0_LDP+ED_ON_NA(I))*T(I)+PTURB_ON_NA)
	    SCL_HT=1.0E-10_LDP/SCL_HT
!
! We set the step size to a fraction of the pressure scale height.
!
	    H=SCL_HT/40.0_LDP
!
!
!
! The boudary condition for the integration of the hydrostatic equation
! has been set, either at the outer boundary, or at the wind connection point.
! we can now perform the integration of the hydrostatic equation.
!
	  IST=I
	  DO WHILE( TAU(I) .LT. MAX(100.0_LDP,OLD_TAU_MAX) )
	    I=I+1
	    IF(I .GT. ND_MAX)THEN
	      WRITE(LU_ERR,*)'Error id DO_CMF_HYDRO'
	      WRITE(LU_ERR,*)'ND_MAX is too small'
	      STOP
	    END IF
!
	    IF(VERBOSE_OUTPUT)THEN
	      IF(MOD(I-IST-1,10) .EQ. 0)THEN
	        WRITE(LUV,*)' '
	        WRITE(LUV,'(A4,10A14)')'I','P(I-1)','Pgas','Pturb','R(I-1)','V(I-1)','T(I-1)',
	1                            'TAU(I-1)','ED/NA','NA(I-1)','GAM_FULL'
	      END IF
	      WRITE(LUV,'(I4,10ES14.5)')I-1,P(I-1),BC*T(I-1)*POP_ATOM(I-1)*(1+ED_ON_NA(I-1)),
	1                    0.5D+10*POP_ATOM(I-1)*MU_ATOM*AMU*VTURB*VTURB,R(I-1),V(I-1),T(I-1),
	1                    TAU(I-1),ED_ON_NA(I-1),POP_ATOM(I-1),GAM_FULL
	    END IF
!
! Set estimates at current location. Then integrate hydrostatic
! equation using 4th order Runge-Kutta.
!
	    CALL OUT_ESTIMATES('INIT',T1)
!
100	    DPTH_INDX=I
	    dR=H
	    STEP_SIZE=H
	    V_EST=V(I-1)
	    TAU_EST=TAU(I-1)
	    T_EST=T(I-1)
	    ED_ON_NA_EST=ED_ON_NA(I-1)
	    ATOM_EST=POP_ATOM(I-1)
	    R_EST=R(I-1)-H
	    dASQdR_EST=dASQdR(I-1)
!
	    CALL CMF_HYDRO_DERIVS( )
	    dV1=H*dVdR
	    dTAU1=H*dTAUdR
!
	    R_EST=R(I-1)-0.5_LDP*H
	    V_EST=V(I-1)-dV1/2
	    TAU_EST=TAU(I-1)-dTAU1/2
	    ATOM_EST=MDOT/R_EST/R_EST/V_EST/MU_ATOM
	    dR=0.5_LDP*H
	    CALL CMF_HYDRO_NEW_EST(ASQ(I-1))
	    CALL CMF_HYDRO_DERIVS( )
	    dV2=H*dVdR
	    dTAU2=H*dTAUdR
!
	    R_EST=R(I-1)-0.5_LDP*H
	    V_EST=V(I-1)-dV2/2
	    TAU_EST=TAU(I-1)-dTAU2/2
	    ATOM_EST=MDOT/R_EST/R_EST/V_EST/MU_ATOM
	    CALL CMF_HYDRO_NEW_EST(ASQ(I-1))
	    CALL CMF_HYDRO_DERIVS( )
	    dV3=H*dVdR
	    dTAU3=H*dTAUdR
!
	    dR=H
	    R_EST=R(I-1)-H
	    V_EST=V(I-1)-dV3
	    TAU_EST=TAU(I-1)-dTAU3
	    ATOM_EST=MDOT/R_EST/R_EST/V_EST/MU_ATOM
	    CALL CMF_HYDRO_NEW_EST(ASQ(I-1))
	    CALL CMF_HYDRO_DERIVS( )
	    dV4=H*dVdR
	    dTAU4=H*dTAUdR
!
! Update values at next grid point.
!
	    TAU(I)=TAU(I-1)-(dTAU1+2*dTAU2+2*dTAU3+dTAU4)/6.0_LDP
	    V(I)=V(I-1)-(dV1+2*dV2+2*dV3+dV4)/6.0_LDP
	    TAU_EST=TAU(I); V_EST=V(I)
!
! Check step size was not too large.
!
	    T1=100.0_LDP
            IF(V(I) .LT. 0.8*V(I-1))T1=V(I)/V(I-1)/0.8_LDP
            T2=LOG(TAU(I)/TAU(I-1))
            IF(T2 .GT. 0.2)T1=MIN(T1,0.2_LDP/T2)
	    IF(T1 .LT. 1.0_LDP)THEN
	       H=H*T1
	       GOTO 100
	    END IF
!
	    dR=H
	    R_EST=R(I-1)-H
	    R(I)=R(I-1)-H
	    POP_ATOM(I)=MDOT/R(I)/R(I)/V(I)/MU_ATOM
	    MASS_DENSITY(I)=MU_ATOM*AMU*POP_ATOM(I)
	    CALL CMF_HYDRO_NEW_EST(ASQ(I-1))
	    T(I)=T_EST
	    ED_ON_NA(I)=ED_ON_NA_EST
	    P(I)=POP_ATOM(I)*(1.0_LDP+ED_ON_NA(I))*BC*T(I)
	    ED(I)=ED_ON_NA(I)*POP_ATOM(I)
	    CHI_ROSS(I)=ED_ON_NA(I)*SIGMA_TH*POP_ATOM(I)*ROSS_ON_ES
	    GAMMA_FULL(I)=GAM_FULL
	    ASQ(I)=1.0E-10_LDP*BC*T(I)*(1.0_LDP+ED_ON_NA(I))/MU_ATOM/AMU
	    dASQdR(I)=dASQdR_EST
	   ND=I
!
	  END DO		!Loop over inner atmosphere
	  CLOSE(UNIT=75)
	  CLOSE(UNIT=76)
!
! Output estimates at the last depth.
!
	  I=ND
	  IF(VERBOSE_OUTPUT)THEN
	    WRITE(LUV,'(I4,10E14.5)')I,P(I-1),BC*T(I-1)*POP_ATOM(I-1)*(1+ED_ON_NA(I-1)),
	1                    0.5D+10*POP_ATOM(I-1)*MU_ATOM*AMU*VTURB*VTURB,R(I-1),V(I-1),T(I-1),
	1                    TAU(I-1),ED_ON_NA(I-1),POP_ATOM(I-1),GAM_FULL
	    WRITE(LUV,'(A4,9A10)')'I','P(I-1)','Pgas','Pturb','R(I-1)','V(I-1)','T(I-1)',
	1                            'TAU(I-1)','ED/NA','NA(I-1)','GAM_FULL'
	    WRITE(LUV,*)' '
	  END IF
!
! Output diagnostic files. These are on the calculate grid --- not the final
! grid.
!
	  ALLOCATE(COEF(ND,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ERR,*)'Error allocating COEF DO_CMF_HYDRO_V2'
	    WRITE(LU_ERR,*)'ND=',ND
	    STOP
	  END IF
	  CALL MON_INT_FUNS_V2(COEF,P,R,ND)
	  DO I=1,ND
	    dPdR_TERM(I)=1.0E-10_LDP*COEF(I,3)/POP_ATOM(I)/AMU/MU_ATOM
	  END DO
	  CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	  DO I=1,ND
	    VdVdR_TERM(I)=V(I)*COEF(I,3)
	    GRAV_TERM(I)=(10.0**LOGG*(REFERENCE_RADIUS/R(I))**2)
	    RAD_TERM(I)=GRAV_TERM(I)*GAMMA_FULL(I)
	  END DO
	  CALL MON_INT_FUNS_V2(COEF,ASQ,R,ND)
	  dASQdR_MON(1:ND)=COEF(1:ND,3)
	  DEALLOCATE(COEF)
!
	  IF(VERBOSE_OUTPUT)THEN
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)THEN
	      WRITE(LU,'(/,A,/)')FORMFEED
	    ELSE
	      OPEN(UNIT=LU,FILE='NEW_CALC_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	    END IF
	    WRITE(LU,'(A,I6,3I16,14I13)')'!',(I,I=1,18)
	    WRITE(LU,'(A,3A16,18A13)')'! Index','R','Vel','Tau','T','Pgas','Rho',
	1                  'Na','Ne/Na','a^2','da2dr','da2dRe','Kross','Kr/Kes','Gamma',
	1                  'dpdR/ROH','VdVdR_TERM','GRAV_TERM','RAD_TERM','SUM'
	    DO I=1,ND
	      T1=1.0E+10_LDP*MASS_DENSITY(I)
	      T2=SIGMA_TH*POP_ATOM(I)*ED_ON_NA(I)
	      T3=1.0D-10*(1.0_LDP+ED_ON_NA(I))*BC*T(I)/MU_ATOM/AMU
	      HYDRO_SUM=(dPdR_TERM(I)+VdVdR_TERM(I)+GRAV_TERM(I)-RAD_TERM(I))/GRAV_TERM(I)
	      WRITE(LU,'(I7,3ES16.7,18ES13.4)')I,R(I),V(I),TAU(I),T(I),P(I),MASS_DENSITY(I),
	1              POP_ATOM(I),ED_ON_NA(I),T3,dASQdR_MON(I),dASQdR(I),
	1              CHI_ROSS(I)/T1,CHI_ROSS(I)/T2,GAMMA_FULL(I),
	1              dPdR_TERM(I),VdVdR_TERM(I),GRAV_TERM(I),RAD_TERM(I),HYDRO_SUM
	    END DO
	    CLOSE(LU)
	  END IF
	  CALL CHECK_HYD(R,V,T,MASS_DENSITY,POP_ATOM,ED,P,GAMMA_FULL,GRAV_TERM,AMU,BC,MU_ATOM,ND,LU)

!
! Adjust the grid so that we get the correct reference radius,
! defined as Tau(Ross)=TAU_REF.
!
	  T1=TAU_REF
	  I=GET_INDX_DP(T1,TAU,ND)
	  T2=(LOG(T1)-LOG(TAU(I)))/(LOG(TAU(I+1))-LOG(TAU(I)))
	  T2=(1.0_LDP-T2)*R(I)+T2*R(I+1)
	  T1=T2-REFERENCE_RADIUS
	  RADIUS_AT_TAU_23=T2
	  IF(WIND_PRESENT .AND. USE_OLD_VEL .AND. HYDRO_OPT .EQ. 'FIXED_LUM')THEN
	    PREV_REF_RADIUS=REFERENCE_RADIUS
	  ELSE IF(WIND_PRESENT)THEN
	    PREV_REF_RADIUS=T2
	    IF(ITERATION_COUNTER .GT. 50)THEN
	      CONNECTION_RADIUS=CONNECTION_RADIUS-0.2_LDP*T1
	    ELSE
	      CONNECTION_RADIUS=CONNECTION_RADIUS-T1
	    END IF
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LUV,*)'    Old reference radius is',T2
	      WRITE(LUV,*)'Desired reference radius is',REFERENCE_RADIUS
	    END IF
	  ELSE IF(PLANE_PARALLEL_MOD)THEN
!
! For a plane-parallel model, the reference radius is the inner boundary,
! defined at TAU=100, not the radius at which tau=2/3.
!
	    T1=100.0_LDP
	    I=GET_INDX_DP(T1,TAU,ND)
	    T2=(LOG(T1)-LOG(TAU(I)))/(LOG(TAU(I+1))-LOG(TAU(I)))
	    T2=(1.0_LDP-T2)*R(I)+T2*R(I+1)
	    T1=T2-REFERENCE_RADIUS
	    PREV_REF_RADIUS=T2
	    RBOUND=RBOUND-T1
	    NI_ZERO=NI_ZERO*REFERENCE_RADIUS/PREV_REF_RADIUS
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LUV,*)'    Old reference radius is',T2
	      WRITE(LUV,*)'Desired reference radius is',REFERENCE_RADIUS
	    END IF
	  ELSE
	    PREV_REF_RADIUS=T2
	    RBOUND=RBOUND-T1
	    NI_ZERO=NI_ZERO*REFERENCE_RADIUS/PREV_REF_RADIUS
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LUV,*)'    Old reference radius is',T2
	      WRITE(LUV,*)'Desired reference radius is',REFERENCE_RADIUS
	    END IF
	  END IF
!
	  IF(ITERATION_COUNTER .GE. 100)THEN
	    WRITE(LU_ERR,*)'Exceed iteration count in DO_CMF_HYDRO_V3.'
	    WRITE(LU_ERR,*)'Aborting update of the hydro structure.'
	    WRITE(LU_ERR,*)'Iteration conunt =',ITERATION_COUNTER
	    IF(VERBOSE_OUTPUT)CLOSE(UNIT=LUV)
	    RETURN
	  END IF
	END DO			!Loop to set R(Tau=2/3)=REFERENCE_RADIUS
!
!	IF(WIND_PRESENT .AND. USE_OLD_VEL .AND. RESET_REF_RADIUS)THEN
!	  REFERENCE_RADIUS=RADIUS_AT_TAU_23
!	END IF
!
	IF(VERBOSE_OUTPUT)THEN
	  CLOSE(UNIT=LUV); CLOSE(LUV_KUTTA)
	END IF
!
!
! We now create the revised grid, At present it is equally spaced in Log(tau) with
! 2 extra points at either end of the grid.
!
	NEW_ND=MOD_ND
	ALLOCATE (REV_TAU(NEW_ND))
	ALLOCATE (REV_R(NEW_ND))
	ALLOCATE (REV_V(NEW_ND))
	ALLOCATE (REV_SIGMA(NEW_ND))
	ALLOCATE (REV_POP_ATOM(NEW_ND))
	ALLOCATE (REV_ED(NEW_ND))
	ALLOCATE (REV_T(NEW_ND))
	ALLOCATE (REV_CHI_ROSS(NEW_ND))
	ALLOCATE (REV_GAMMA_FULL(NEW_ND))
	ALLOCATE (COEF(NEW_ND,4))
!
	TAU_MAX=100.0_LDP
	IF(TAU_MAX .GT. TAU(ND))THEN
	  WRITE(LU_ERR,*)'Error --- TAU_MAX cannot be greater than calculated grid Tau'
	  WRITE(LU_ERR,*)'Setting TAU to maximum value in DO_CMF_HYDRO'
	  TAU_MAX=TAU(ND)
	END IF
!
	DO I=1,ND
	  V(I)=MDOT/MU_ATOM/POP_ATOM(I)/R(I)/R(I)
	END DO
	IF(WIND_PRESENT .AND. USE_OLD_VEL)THEN
	  J=CONNECTION_INDX
	  I=NEW_ND-J+1
	  CALL DET_R_GRID_V1(REV_TAU(J),I,ND_MAX,TAU_MAX,L_FALSE,
	1                    R(J),V(J),TAU(J),ND-J+1)
	  REV_TAU(1:J)=TAU(1:J)
	ELSE
	  CALL DET_R_GRID_V2(REV_TAU,NEW_ND,ND_MAX,TAU_MAX,
	1        dLOG_TAU,V_SCL_FAC,OUT_BND_OPT,OBND_PARS,NUM_OBND_PARAMS,
	1        R,V,TAU,ND)
	END IF
!
! We now compute the revised R grid. We then interplate on Log (r^2.rho) which
! is equivalent to interpolating on log V. This guarentees monotocity of V.
!
	TAU(1:ND)=LOG(TAU(1:ND))
	REV_TAU(1:NEW_ND)=LOG(REV_TAU(1:NEW_ND))
	CALL MON_INTERP(REV_R,NEW_ND,IONE,REV_TAU,NEW_ND,R,ND,TAU,ND)
	IF(PLANE_PARALLEL_MOD)THEN
	  POP_ATOM(1:ND)=LOG(POP_ATOM(1:ND))
	  CALL MON_INTERP(REV_POP_ATOM,NEW_ND,IONE,REV_R,NEW_ND,POP_ATOM,ND,R,ND)
	  REV_POP_ATOM(1:NEW_ND)=EXP(REV_POP_ATOM(1:NEW_ND))
	  POP_ATOM(1:ND)=EXP(POP_ATOM(1:ND))
	ELSE
	  POP_ATOM(1:ND)=LOG(POP_ATOM(1:ND)*R(1:ND)*R(1:ND))
	  CALL MON_INTERP(REV_POP_ATOM,NEW_ND,IONE,REV_R,NEW_ND,POP_ATOM,ND,R,ND)
	  REV_POP_ATOM(1:NEW_ND)=EXP(REV_POP_ATOM(1:NEW_ND))/REV_R(1:NEW_ND)/REV_R(1:NEW_ND)
	  POP_ATOM(1:ND)=EXP(POP_ATOM(1:ND))/R(1:ND)/R(1:ND)
	END IF
!
! This set the reference radius to R(ND), and keeps it the same as was read in.
!
	IF(PLANE_PARALLEL_MOD)THEN
	  DO I=1,ND
	    R(I)=R(I)+(REFERENCE_RADIUS-REV_R(NEW_ND))
	  END DO
	  DO I=1,NEW_ND
	    REV_R(I)=REV_R(I)+(REFERENCE_RADIUS-REV_R(NEW_ND))
	  END DO
	END IF
!
! Compute revised velocity.
!
	IF(WIND_PRESENT .AND. USE_OLD_VEL)THEN
	  REV_V(1:CONNECTION_INDX)=OLD_V(1:CONNECTION_INDX)
	  DO I=CONNECTION_INDX+1,NEW_ND
	    REV_V(I)=MDOT/MU_ATOM/REV_POP_ATOM(I)/REV_R(I)/REV_R(I)
	  END DO
	ELSE
	  DO I=1,NEW_ND
	    REV_V(I)=MDOT/MU_ATOM/REV_POP_ATOM(I)/REV_R(I)/REV_R(I)
	  END DO
	END IF
!
! Compute SIGMA by performing a monotonic cubic fit to V as a function of R.
!
	CALL MON_INT_FUNS_V2(COEF,REV_V,REV_R,NEW_ND)
	DO I=1,NEW_ND
	  REV_SIGMA(I)=REV_R(I)*COEF(I,3)/REV_V(I)-1.0_LDP
	END DO
!
!
!
! Saves the current RVSIG_FILE for recovery puposes. Except for the first
! iteration, this could be recovered from RVSIG_COL. For portability, we
! use only regular fortran commands.
!
	INQUIRE(FILE='RVSIG_COL',EXIST=FILE_PRES)
	IF(FILE_PRES)THEN
	  STRING=' '
	  WRITE(STRING,'(I4.4)')MAIN_COUNTER
	  STRING='RVSIG_COL_IT_'//TRIM(STRING)
   	  OPEN(UNIT=LUIN,FILE='RVSIG_COL',STATUS='OLD',ACTION='READ')
	  OPEN(UNIT=LU,FILE=TRIM(STRING),STATUS='UNKNOWN',ACTION='WRITE')
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    WRITE(LU,'(A)')TRIM(STRING)
	    IF(IOS .NE. 0)EXIT
	  END DO
	  CLOSE(LUIN)
	  CLOSE(LU)
	END IF
!
! Output revised hydrostatic structure. This can be used to restart the current
! model from scratch.
!
	OPEN(UNIT=LU,FILE='RVSIG_COL',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(A)')'! Note: The effective temperature and surface gravity are defined'
	  WRITE(LU,'(A)')'! at the reference radius, which (except when using the old'
	  WRITE(LU,'(A)')'! velocity) is the location where Tau=2/3.'
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(A,ES16.6)')'! Effective temperature (10^4 K) is:',TEFF
	  WRITE(LU,'(A,ES16.6)')'!      Log surface gravity (cgs) is:',LOGG
	  WRITE(LU,'(A,ES16.6)')'!         Core radius (10^10 cm) is:',REV_R(NEW_ND)
	  WRITE(LU,'(A,ES16.6)')'!    Reference radius (10^10 cm) is:',REFERENCE_RADIUS
	  WRITE(LU,'(A,ES16.6)')'!              Luminosity (Lsun) is:',( (TEFF/0.5770_LDP)**4 )*( (REFERENCE_RADIUS/6.9599_LDP)**2 )
	  WRITE(LU,'(A,ES16.6)')'!            Mass (Msun) of star is:',10.0_LDP**(LOGG)*REFERENCE_RADIUS*REFERENCE_RADIUS/GRAV_CON
	  WRITE(LU,'(A,ES16.6)')'!       Mass loss rate (Msun/yr) is:',MDOT/MASS_LOSS_SCALE_FACTOR
	  WRITE(LU,'(A,ES16.6)')'!         Mean atomic mass (amu) is:',MU_ATOM
	  WRITE(LU,'(A,ES16.6)')'!            Eddington parameter is:',GAM_EDD
	  WRITE(LU,'(A,ES16.6)')'!                   Atom density is:',NI_ZERO
	  WRITE(LU,'(A,F14.8)') '! Ratio of inner to outer radius is:',REV_R(1)/REV_R(NEW_ND)
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(3X,I5,10X,A)')NEW_ND,'!Number of depth points'
	  WRITE(LU,'(A)')'!'
	  IF(REV_R(1) .GT. 999999.0_LDP)THEN
	    WRITE(LU,'(A,4X,A,3(7X,A),3X,A)')'!','R(10^10cm)','V(km/s)','  Sigma','    Tau','  Index'
	    DO I=1,NEW_ND
	      WRITE(LU,'(F18.8,ES16.8,2ES14.6,6X,I4)')REV_R(I),REV_V(I),REV_SIGMA(I),EXP(REV_TAU(I)),I
	    END DO
	  ELSE IF( (REV_R(1)-REV_R(NEW_ND)) .LT. 1.0_LDP)THEN
	    WRITE(LU,'(A,10X,A,3(7X,A),3X,A)')'!','R(10^10cm)','V(km/s)','  Sigma','    Tau','  Index'
	    DO I=1,NEW_ND
	      WRITE(LU,'(F24.12,ES16.8,2ES14.6,6X,I4)')REV_R(I),REV_V(I),REV_SIGMA(I),EXP(REV_TAU(I)),I
	    END DO
	  ELSE
	    WRITE(LU,'(A,1X,A,3(7X,A),3X,A)')'!','R(10^10cm)','V(km/s)','  Sigma','    Tau','  Index'
	    DO I=1,NEW_ND
	      WRITE(LU,'(F15.8,ES16.8,2ES14.6,6X,I4)')REV_R(I),REV_V(I),REV_SIGMA(I),EXP(REV_TAU(I)),I
	    END DO
	  END IF
	CLOSE(UNIT=LU)
!
!
! Output estimate data for comparison with new model data. R is used as the
! dependent variable.
!
	CALL MON_INTERP(REV_T,NEW_ND,IONE,REV_R,NEW_ND,T,ND,R,ND)
	CALL MON_INTERP(REV_ED,NEW_ND,IONE,REV_R,NEW_ND,ED,ND,R,ND)
	CALL MON_INTERP(REV_CHI_ROSS,NEW_ND,IONE,REV_R,NEW_ND,CHI_ROSS,ND,R,ND)
	CALL MON_INTERP(REV_GAMMA_FULL,NEW_ND,IONE,REV_R,NEW_ND,GAMMA_FULL,ND,R,ND)
!
	OPEN(UNIT=LU,FILE='FIN_CAL_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(41,'(A,9(6X,A7))')'Index','     R','   TAU','    V','     T','    Na',
	1              ' Ne/Na',' Xross','Xr/Xes','   Gam'
	  DO I=1,NEW_ND
	    WRITE(LU,'(I5,9ES13.4)')I,REV_R(I),EXP(REV_TAU(I)),REV_V(I),REV_T(I),
	1             REV_POP_ATOM(I),REV_ED(I)/REV_POP_ATOM(I),REV_CHI_ROSS(I),
	1             REV_CHI_ROSS(I)/SIGMA_TH/REV_ED(I),REV_GAMMA_FULL(I)
	  END DO
	CLOSE(UNIT=LU)
!
	CALL SET_NEW_GRID_V2(REV_R,REV_T,REV_V,REV_SIGMA,REV_ED,REV_CHI_ROSS,NEW_ND)
	CALL SET_ABUND_CLUMP(T1,T2,LU_ERR,NEW_ND)
!	REV_ED(1:ND)=REV_ED(1:ND)/CLUMP_FAC(1:ND)
!	CALL SET_NEW_GRID_V2(REV_R,REV_T,REV_V,REV_SIGMA,REV_ED,REV_CHI_ROSS,NEW_ND)
	CALL ADJUST_POPS(POPS,LU_ERR,NEW_ND,NT)
!
	NO_HYDRO_ITS=NO_HYDRO_ITS-1
	NO_ITS_DONE=NO_ITS_DONE+1
	CALL UPDATE_KEYWORD(NO_HYDRO_ITS,'[N_ITS]','HYDRO_DEFAULTS',L_TRUE,L_FALSE,LUIN)
	CALL UPDATE_KEYWORD(NO_ITS_DONE,'[ITS_DONE]','HYDRO_DEFAULTS',L_FALSE,L_TRUE,LUIN)
	DONE_HYDRO=.TRUE.
!
! Make sure VADAT is consitent with revised RGRID & parameters in HYDRO_PARAMS.
! We also return the correct RSTAR, RMAX and stellar mass.
!
	MOD_RSTAR=REV_R(NEW_ND)
	MOD_RMAX=REV_R(1)
	MOD_MASS=10.0_LDP**(LOGG)*(REFERENCE_RADIUS**2)/GRAV_CON
	CALL UPDATE_KEYWORD(REV_R(NEW_ND),'[RSTAR]','VADAT',L_TRUE,L_FALSE,LUIN)
	T1=MOD_RMAX/REV_R(NEW_ND)
	CALL UPDATE_KEYWORD(T1,'[RMAX]','VADAT',L_FALSE,L_FALSE,LUIN)
	CALL UPDATE_KEYWORD(MOD_MASS,'[MASS]','VADAT',L_FALSE,L_TRUE,LUIN)
!	CALL UPDATE_KEYWORD(REV_V(1),'[VINF]','VADAT',L_FALSE,L_TRUE,LUIN)
	WRITE(LU_ERR,'(A)')' Revised hydrostatic structure and output new RVSIG_COL file'
	WRITE(LU_ERR,'(A)')' Updated RSTAR, RMAX and MASS in VADAT'
!
	RETURN
	END SUBROUTINE DO_CMF_HYDRO_V3
!	
! Subroutine to compute dPdR and dTAUdR for use with the
! Runge-Kutta integration.
!
	SUBROUTINE CMF_HYDRO_DERIVS( )
	USE CMF_SUB_SONIC_HYDRO
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	REAL(KIND=LDP) dlnTdR
	REAL(KIND=LDP) SURF_GRAV
	REAL(KIND=LDP) SOUND_SPEED_SQ
	REAL(KIND=LDP) MOD_SOUND_SPEED_SQ
	REAL(KIND=LDP) dGTAU
	REAL(KIND=LDP) T1
!
! The factor of 1.0E-10 is applied to get the sound speed in km/s.
! The constant BC already has a factor of 10^4 included to compensate for T being in units of 10^4 K.
! 
	T1=1.0E-10_LDP*BC/MU_ATOM/AMU
	SOUND_SPEED_SQ=T1*(1.0_LDP+ED_ON_NA_EST)*T_EST
	MOD_SOUND_SPEED_SQ=SOUND_SPEED_SQ+0.5_LDP*VTURB*VTURB

	SURF_GRAV=(10.0_LDP**LOGG)
	IF(.NOT.PLANE_PARALLEL_MOD)THEN
	  SURF_GRAV=SURF_GRAV*(REFERENCE_RADIUS/R_EST)**2
	END IF
!
! Notes: Since R is units of 10^10 cm and V in km/s, VdV/dR has cgs units (cm s^{-2})
!
	dTAUdR=-ED_ON_NA_EST*ATOM_EST*SIGMA_TH*ROSS_ON_ES
	T1=V_EST/(V_EST*V_EST-MOD_SOUND_SPEED_SQ)
	dVdR=T1*(2*MOD_SOUND_SPEED_SQ/R_EST - dASQdR_EST- SURF_GRAV*(1.0_LDP-GAM_FULL))
!
	CALL OUT_ESTIMATES('H',STEP_SIZE)
	CALL OUT_ESTIMATES('R_EST',R_EST)
	CALL OUT_ESTIMATES('V_EST',V_EST)
	CALL OUT_ESTIMATES('T_EST',T_EST)
	CALL OUT_ESTIMATES('TAU_EST',TAU_EST)
	T1=MU_ATOM*AMU*ATOM_EST;                           CALL OUT_ESTIMATES('DENS_EST',T1)
	CALL OUT_ESTIMATES('ATOM_EST',ATOM_EST)
	CALL OUT_ESTIMATES('ED_ON_NA',ED_ON_NA_EST)
	IF(VTURB .NE. 0)THEN
	  T1=0.5E+10_LDP*MU_ATOM*AMU*ATOM_EST*VTURB*VTURB; CALL OUT_ESTIMATES('    P(TURB)=',T1)
	END IF
	T1=BC*ATOM_EST*T_EST*(1+ED_ON_NA_EST);             CALL OUT_ESTIMATES('     P(GAS)=',T1)
	CALL OUT_ESTIMATES('dVdR',dVdR)
	CALL OUT_ESTIMATES('dTAUdR',dTAUdR)
	CALL OUT_ESTIMATES('ROSS_ON_ES',ROSS_ON_ES)
	CALL OUT_ESTIMATES('MC^2',MOD_SOUND_SPEED_SQ)
	CALL OUT_ESTIMATES('C^2',SOUND_SPEED_SQ)
	CALL OUT_ESTIMATES('2*MC^2/R',2*MOD_SOUND_SPEED_SQ/R_EST)
	CALL OUT_ESTIMATES('dASQdR',dASQdR_EST)
	CALL OUT_ESTIMATES('SURF_GRAV',SURF_GRAV)
	CALL OUT_ESTIMATES('GAM_FULL',GAM_FULL)
	CALL OUT_ESTIMATES('VTURB',VTURB)
	CALL OUT_ESTIMATES('LAST',VTURB)
	FLUSH(UNIT=LUV_KUTTA)
!
	RETURN
	END SUBROUTINE CMF_HYDRO_DERIVS
!
! Get new estimates of the model parameters. At present we
! can choose pure LTE estimates, or scaled LTE estimates.
!
	SUBROUTINE CMF_HYDRO_NEW_EST(ASQ_OLD)
	USE CMF_SUB_SONIC_HYDRO
	USE OLD_GRID_MODULE
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	REAL(KIND=LDP) ASQ_OLD
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) T2
	REAL(KIND=LDP) ED_EST
	REAL(KIND=LDP) SF
	REAL(KIND=LDP) FM
	REAL(KIND=LDP) ES
	REAL(KIND=LDP) RM
!
	REAL(KIND=LDP) OLD_ATOM
	REAL(KIND=LDP) LTE_OLD_ATOM
	REAL(KIND=LDP) OLD_TEMP
	REAL(KIND=LDP) LTE_ED
	REAL(KIND=LDP) LTE_OLD_ED
	REAL(KIND=LDP) KAP_ROSS_OLD
	REAL(KIND=LDP) KAP_ROSS_LTE
	REAL(KIND=LDP) KAP_ROSS_EST
	REAL(KIND=LDP) KAP_ES_OLD
	REAL(KIND=LDP) KAP_ES_LTE
!
	INTEGER, SAVE ::INDX=1
!
	IF(PURE_LTE_EST)THEN
	  T_EST=TEFF*(TAU_EST+0.67_LDP)**0.25_LDP
	  CALL GET_LTE_ROSS_V2(KAP_ROSS,KAP_ES,LTE_ED,ATOM_EST,T_EST)
	  ED_ON_NA_EST=LTE_ED/ATOM_EST
	  GAM_FULL=GAM_EDD*(KAP_ROSS/KAP_ES)*ED_ON_NA_EST
	  ROSS_ON_ES=KAP_ROSS/KAP_ES
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	  RETURN
	END IF
!
! Determine parameters at current grid point. Special allowance
! has to be made because of the boundaries. As we are access`ing
! the grid sequentially we can use INDX as the starting location
! of the search. We just need to check that we are not starting
! a new sequence.
!
	IF(TAU_EST .LT. OLD_TAU(INDX))INDX=1
	DO WHILE(TAU_EST .GT. OLD_TAU(INDX) .AND. INDX .LT. OLD_ND)
	  INDX=INDX+1
	END DO
	IF(INDX .EQ. 1)THEN
	  T_EST=T_EST     				!Use value passed
	  ED_ON_NA_EST=OLD_ED(1)/OLD_POP_ATOM(1)
	  ROSS_ON_ES=OLD_KAP_ROSS(1)/OLD_KAP_ESEC(1)
	  GAM_FULL=GAM_EDD*(OLD_KAP_FLUX(1)/OLD_KAP_ESEC(1))*ED_ON_NA_EST
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
!
! Fix to use LTE
!
	ELSE IF(TAU_EST .GE. OLD_TAU(OLD_ND))THEN
	  T_EST=TEFF*OLD_SF(OLD_ND)*(TAU_EST+0.67_LDP)**0.25_LDP
	  T2=OLD_T(OLD_ND)-OLD_T(OLD_ND-1)
	  ED_ON_NA_EST= OLD_ED(OLD_ND)/OLD_POP_ATOM(OLD_ND)
	  ROSS_ON_ES=OLD_KAP_ROSS(OLD_ND)/OLD_KAP_ESEC(OLD_ND)
	  GAM_FULL=GAM_EDD*(OLD_KAP_FLUX(OLD_ND)/OLD_KAP_ESEC(OLD_ND))*ED_ON_NA_EST
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	ELSE
!
! Determine parameters at current TAU in old model.
!
	  T1=(TAU_EST-OLD_TAU(INDX-1))/(OLD_TAU(INDX)-OLD_TAU(INDX-1))
	  SF=(1.0_LDP-T1)*OLD_SF(INDX-1)+T1*OLD_SF(INDX)
	  T_EST=SF*TEFF*(TAU_EST+0.67_LDP)**0.25_LDP
!
          ED_EST=(1.0_LDP-T1)*OLD_ED(INDX-1)+T1*OLD_ED(INDX)
          OLD_ATOM= (1.0_LDP-T1)*OLD_POP_ATOM(INDX-1)+T1*OLD_POP_ATOM(INDX)
          OLD_TEMP= (1.0_LDP-T1)*OLD_T(INDX-1)+T1*OLD_T(INDX)
          ED_ON_NA_EST=ED_EST/OLD_ATOM
	  IF(ATOM_EST .GT. 1.0E+09_LDP)THEN
	    CALL GET_LTE_ROSS_V2(KAP_ROSS_OLD,KAP_ES_OLD,LTE_OLD_ED,OLD_ATOM,OLD_TEMP)
	    CALL GET_LTE_ROSS_V2(KAP_ROSS_LTE,KAP_ES_LTE,LTE_ED,ATOM_EST,T_EST)
	    ED_EST=ED_EST*(LTE_ED/LTE_OLD_ED)/ATOM_EST
	  END IF
!
	  ASQ_EST=1.0E-10_LDP*BC*T_EST*(1.0_LDP+ED_ON_NA_EST)/MU_ATOM/AMU
	  dASQdR_EST=(ASQ_OLD-ASQ_EST)/dR
	  WRITE(40,*)T_EST,BC,MU_ATOM,AMU,1.0_LDP-ED_ON_NA_EST
	  WRITE(40,*)ASQ_EST,ASQ_OLD,dR,dASQdR_EST; FLUSH(UNIT=40)
!
	  FM=(1.0_LDP-T1)*OLD_KAP_FLUX(INDX-1)+T1*OLD_KAP_FLUX(INDX)
	  RM=(1.0_LDP-T1)*OLD_KAP_ROSS(INDX-1)+T1*OLD_KAP_ROSS(INDX)
	  ES=(1.0_LDP-T1)*OLD_KAP_ESEC(INDX-1)+T1*OLD_KAP_ESEC(INDX)
	  WRITE(6,*)SF,TEFF,T_EST; FLUSH(UNIT=6)
!
! Some fiddling may be required here to choose the optimal density for switching.
!
	  IF(ATOM_EST .GT. 1.0E+09_LDP)THEN
            KAP_ROSS_EST=(KAP_ROSS_LTE-KAP_ES_LTE)*( (RM-ES)/(KAP_ROSS_OLD-KAP_ES_OLD))+ES
	    KAP_ROSS=RM; KAP_ES=ES
	    GAM_FULL=GAM_EDD*(FM/RM)*(0.5_LDP*(KAP_ROSS_EST+RM)/KAP_ES)*ED_ON_NA_EST
	    ROSS_ON_ES=KAP_ROSS/KAP_ES
	    WRITE(76,'(ES16.8,3ES12.4,/,16X,3ES12.4,/,16X,ES12.4,/,16X,4ES12.4)')
	1                       R_EST,KAP_ROSS_OLD,KAP_ES_OLD,OLD_TEMP,
	1                       KAP_ROSS_LTE,KAP_ES_LTE,T_EST,
	1                       KAP_ROSS_EST,
	1                       RM,ES,FM,GAM_FULL
	  ELSE
	    GAM_FULL=GAM_EDD*(FM/ES)*ED_ON_NA_EST
	    ROSS_ON_ES=RM/ES
	  END IF
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	  WRITE(75,'(ES16.9,8ES12.3)')R_EST,GAM_FULL,ED_ON_NA_EST*GAM_EDD*FM/ES,TAU_EST,FM,RM,KAP_ROSS,KAP_ES,ED_ON_NA_EST
	END IF
!
	RETURN
	END SUBROUTINE CMF_HYDRO_NEW_EST

	SUBROUTINE STORE_OLD_GRID(MEAN_ATOMIC_MASS,ND)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE OLD_GRID_MODULE
	USE CMF_SUB_SONIC_HYDRO
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER IOS
	INTEGER ISPEC
	REAL(KIND=LDP) MEAN_ATOMIC_MASS
!
	MEAN_ATOMIC_MASS=0.0_LDP
	DO ISPEC=1,NUM_SPECIES
	    MEAN_ATOMIC_MASS=MEAN_ATOMIC_MASS+POP_SPECIES(ND,ISPEC)*AT_MASS(ISPEC)
	END DO
	MEAN_ATOMIC_MASS=MEAN_ATOMIC_MASS/POP_ATOM(ND)
!
	IF(.NOT. ALLOCATED(OLD_R))THEN
	  ALLOCATE (OLD_R(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_V(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_SIGMA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_T(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_ED(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_ROSS_MEAN(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_FLUX_MEAN(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_ESEC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_KAP_ROSS(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_KAP_FLUX(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_KAP_ESEC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_POP_ATOM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_MASS_DENSITY(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_POPION(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_CLUMP_FAC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_TAU(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_SF(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ERR,*)'Error in STORE_OLD_GRID -- error allocating atmospheric vectors'
	    WRITE(LU_ERR,*)'STATUS=',IOS
	    STOP
	  END IF
	END IF
!
	OLD_ND=ND
!
! Store vectors describing the old atmospheric structure.
!
	OLD_R=R
	OLD_T=T
	OLD_ED=ED
	OLD_SIGMA=SIGMA
	OLD_V=V
!
	OLD_POP_ATOM=POP_ATOM
	OLD_POPION=POPION
	OLD_CLUMP_FAC=CLUMP_FAC
	OLD_MASS_DENSITY=DENSITY
!
	OLD_ROSS_MEAN=ROSS_MEAN
	OLD_FLUX_MEAN=FLUX_MEAN
	OLD_ESEC=6.65E-15_LDP*OLD_ED
!
! Compute mass absorption coefficients in cgs units.
!
	OLD_KAP_ROSS=1.0E-10_LDP*OLD_ROSS_MEAN/OLD_MASS_DENSITY
	OLD_KAP_FLUX=1.0E-10_LDP*OLD_FLUX_MEAN/OLD_MASS_DENSITY
	OLD_KAP_ESEC=1.0E-10_LDP*OLD_ESEC/OLD_MASS_DENSITY
!
	RETURN
	END SUBROUTINE STORE_OLD_GRID
!
	SUBROUTINE OUT_ESTIMATES(ID,VALUE)
	USE CMF_SUB_SONIC_HYDRO
	USE SET_KIND_MODULE
	IMPLICIT NONE
	CHARACTER(LEN=*) ID
	REAL(KIND=LDP) VALUE
	INTEGER, SAVE :: CNT 
	CHARACTER(LEN=13) LOC_ID
	CHARACTER(LEN=300), SAVE :: HEADER_STRING
	CHARACTER(LEN=300), SAVE :: VALUE_STRING
!
	IF(ID .EQ. 'INIT')THEN
	  HEADER_STRING=' '
	  VALUE_STRING=' '
	  CNT=1 
	ELSE IF(ID .EQ. 'LAST')THEN
	  IF(MOD(CNT,4) .EQ. 1)THEN
	    WRITE(LUV_KUTTA,'(A)')' '
	    WRITE(LUV_KUTTA,'(A,4X,A,2X,A)')'!','I',TRIM(HEADER_STRING)
	    HEADER_STRING=' '
	  END IF
	  WRITE(LUV_KUTTA,'(I6,2X,A)')DPTH_INDX,TRIM(VALUE_STRING)
	  VALUE_STRING=' '
	  CNT=CNT+1
	ELSE
	  IF(INDEX(HEADER_STRING,' '//TRIM(ID)//' ') .EQ. 0)THEN
	    LOC_ID=ID; LOC_ID=ADJUSTR(LOC_ID)
	    HEADER_STRING=TRIM(HEADER_STRING)//LOC_ID
	  END IF
	  WRITE(LOC_ID,'(ES13.4)')VALUE
	  VALUE_STRING=TRIM(VALUE_STRING)//LOC_ID
	END IF
!
	RETURN
	END SUBROUTINE OUT_ESTIMATES
!
	SUBROUTINE CHECK_HYD(R,V,T,MASS_DENSITY,POP_ATOM,ED,P,GAMMA_FULL,GRAV_TERM,AMU,BC,MU_ATOM,ND,LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER ND
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) POP_ATOM(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) P(ND)
	REAL(KIND=LDP) GAMMA_FULL(ND)
	REAL(KIND=LDP) GRAV_TERM(ND)
	REAL(KIND=LDP) MASS_DENSITY(ND)
!
	REAL(KIND=LDP) BC
	REAL(KIND=LDP) AMU
	REAL(KIND=LDP) MU_ATOM
!
	INTEGER LU,I
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) ASQ(ND)
	REAL(KIND=LDP) dvdR(ND)
	REAL(KIND=LDP) dpdR_EST(ND)
	REAL(KIND=LDP) COEF(ND,4)
	REAL(KIND=LDP) dASQdr(ND)
	REAL(KIND=LDP) dpdR(ND)
	REAL(KIND=LDP) T1,T2
!
	OPEN(UNIT=LU,FILE='CHECK_HYD',STATUS='UNKNOWN',ACTION='WRITE')
!
	TA=ED/POP_ATOM
	ASQ=1.0E-10_LDP*BC*T*(1.0_LDP+TA)/MU_ATOM/AMU
	CALL MON_INT_FUNS_V2(COEF,ASQ,R,ND)
	dASQdR=COEF(:,3)
!
	CALL MON_INT_FUNS_V2(COEF,V,R,ND)
        dVdR=COEF(:,3)
!
	dPdR_EST=dASQdR  - 2*ASQ/R -ASQ*dVdR/V
!
	CALL MON_INT_FUNS_V2(COEF,P,R,ND)
        dpdR=1.0E-10_LDP*COEF(:,3)/MASS_DENSITY
!
	WRITE(LU,'(A5,A18,13A14)')'I','R','V','T','P','ED/NA','ASQ','dASQdR','dPdR_EST','dPdR','Gam','g(1-G)','LHS','RHS'
	DO I=1,ND
	  T1=(V(I)*V(I)-ASQ(I))*dVdR(I)/V(I)
	  T2=2*ASQ(I)/R(I)-dASQdR(I)-GRAV_TERM(I)*(1.0_LDP-GAMMA_FULL(I))
	  WRITE(LU,'(I5,ES18.8,12ES14.4)')I,R(I),V(I),T(I),P(I),TA(I),ASQ(I),dASQdR(I),dPdR_EST(I),dPdR(I),
	1                       GAMMA_FULL(I),GRAV_TERM(I)*(1.0_LDP-GAMMA_FULL(I)),T1,T2
	END DO
	CLOSE(UNIT=LU)
!
	RETURN
	END
