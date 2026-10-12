	SUBROUTINE RD_CONTROL_VARIABLES(LUIN,LUSCR,LUER,NUM_BNDS)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_log_unknown, jfnk_ptc_active,
	1     jfnk_ptc_dtau, jfnk_ptc_scale_ready
	IMPLICIT NONE
!
! Altered : 07-Oct-2026 : DO_JFNK is now the shared JFNK arm key for
!                           DO_NONMONO_KNOP and DO_NONMONO_SOB; USE_JFNK is a
!                           deprecated alias.  JFNK_COUPLED stays KNOP-only
!                           (Devin).
! Altered : 05-Oct-2026 : Reject legacy controls that use a multi-valued
!                           or non-positive V coordinate under non-mono
!                           solvers (Devin).
! Altered : 03-Oct-2026 : Non-monotonic mode prerequisites enforced (dev7):
!                            DO_NONMONO_SOB requires GLOBAL_LINE=SOB (or NONE
!                            with all TRANS_*=SOB) and USE_LAM_ES=T;
!                            DO_NONMONO_KNOP requires GLOBAL_LINE=BLANK (or NONE
!                            with all TRANS_*=BLANK) and LAM_SET=F.
!                            DO_BA_RES_ACC=F under KNOP warns.  Previously these
!                            were warnings only, so a mismatched VADAT ran with
!                            the flag inactive.
! Altered : 13-Jan-2024 : Added MAX_NO_GREY_ITERATIONS
! Altered : 19-Jun-2023 : Added SOL_ABUND_REF_SET
! Altered : 12-Aug-2022 : Added SN shock variables (following work by LUC.)
! Altered : 06-Jun-2022 : Variable COMP_STEQ_T_EHB ([COMP_EHB]) added.
! Altered : 15-Nov-2021 : Fixed read in calls for FIX_ALL_SPEC and BC_PAR1.
! Altered : 14-Jul-2019 : R_PNT_SRCE, TEFF_PNT_SRCE, NC_PNT_SRCE, LOGICAL PNT_SRCE_MOD, USE_ELEC_HEAT_BAL (from osiris)
!                            Added to ibis 14-Aug-2019
! Altered : 26-Apr-2019 : Added        J_CHK_OPTION (CHK_J)     - Done earlier on OSIRIS.
!                               AUTO_ADD_ION_STAGES (AUTO_ADD)
!                                  AUTO_SMOOTH_POPS (AUTO_SMOOTH)
! Altered : 01-Sep-2016 : Added variabe MINIMUM_ISO_POP (not required)
!                            TIME_SEQ_NO changed from integer to real.
! Altered : 15-Feb-2015 : Added INSTANTANEOUS_ENERGY_DEPOSITION option (12-Jan-2015 on OSPREY[cur_cmf_gam])
! Altered : 02/07-Nov-2011 : Changed location of COMP_GREY_LST_IT so that set for all model.
! Altered : 05-Apr-2011 : Changed the way REVISE_R_GRID is handled. Variable controlling
!                           R grid revision now read in by routine do the R-grid revision.
! Altered : 16-Jul-2010 : Added FIX_ALL_SPECIES variable, and assoicated options.
! Altered : 01-Feb-2010 : GAMRAY_TRANS read inserted
! Altered : 31-Jan-2010 : INNER_BND_METH and OUTER_BND_METH options installed.
! Altered : 23-Nov-2007 : Optional LAM_SCALE_OPT variable included.
! Altered : 29-Jan-2006 : Control variables for relativistic transfer and time
!                          dependent statistical equilibrium equations installed.
!
	INTEGER LUIN
	INTEGER LUSCR
	INTEGER LUER
	INTEGER NUM_BNDS
!
! Local variables.
!
	CHARACTER(LEN=80)  TMP_STRING
	CHARACTER(LEN=20)  TMP_KEY
	CHARACTER(LEN=132) TEMP_CHAR
!
	INTEGER I
	INTEGER ISPEC
	INTEGER ID
!
	INTEGER KN_M_MAXPASS_V,KN_RAY_MAXPASS_V,MOM_AA_M_V,
	1     MOM_AA_OUT_M_V
	REAL(KIND=LDP) KN_M_TOL_V,KN_RAY_TOL_V
	REAL(KIND=LDP) JNEXT_RELAX_V,INEXT_RELAX_V
!
	LOGICAL USE_JFNK_OLD
!
	REAL(KIND=LDP) ATOMIC_MASS_UNIT
	EXTERNAL ATOMIC_MASS_UNIT
!
	CALL GEN_ASCI_OPEN(LUIN,'VADAT','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening VADAT in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
! 
!
! Input model parameters and modelling specifications.
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(RP,'RSTAR',L_TRUE,'Stellar radius (in 10^10 cm)')
	  CALL RD_STORE_DBLE(RMAX,'RMAX',L_TRUE,'Maximum radius (in R*)')
	  RMAX=RMAX*RP
	  DO_HYDRO=.FALSE.
	  CALL RD_STORE_LOG(DO_HYDRO,'DO_HYDRO',L_FALSE,'Adjust hydrostatic structure')
!
	  CALL RD_STORE_INT(VELTYPE,'VEL_LAW',L_TRUE,'Velocity Law to be used')
	  IF(VELTYPE .EQ. 1 .OR. VELTYPE .EQ. 2)THEN
	    CALL RD_STORE_DBLE(VRP,'VRP',L_TRUE,'First velocity component')
	    CALL RD_STORE_DBLE(RN,'RN',L_TRUE,' ')
	    CALL RD_STORE_DBLE(VINF,'VINF',L_TRUE,' ')
	    CALL RD_STORE_DBLE(EPPS1,'EPSS1',L_TRUE,' ')
	    CALL RD_STORE_DBLE(GAMMA1,'GAMMA1',L_TRUE,' ')
!
	    CALL RD_STORE_DBLE(RP2,'RP2',L_TRUE,'Second velocity component')
	    CALL RD_STORE_DBLE(VRP2,'VRP2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(RN2,'RN2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(VINF2,'VINF2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(EPPS2,'EPPS2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(GAMMA2,'GAMMA2',L_TRUE,' ')
	    RN=RN*RP
	    RP2=RP*RP2
	    RN2=RN2*RP2
	  ELSE IF(VELTYPE .EQ. 3)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,
	1           'Core velocity (km/s)')
	    CALL RD_STORE_DBLE(VPHOT,'VPHOT',L_TRUE,
	1           'Photospheric velocity (km/s)')
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,
	1           'Terminal velocity (km/s)')
	    CALL RD_STORE_DBLE(SCL_HT,'SCL_HT',L_TRUE,
	1           'Scale Height (in R*) of photosphere')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA',L_TRUE,
	1           'Speed of velocity Law')
	    V_EPPS1=1.0_LDP
	    VINF=VINF1
	    VINF2=VINF1                !i.e. no 2nd component
	    V_BETA2=1.0_LDP
	    V_EPPS2=1.0_LDP
	    N_OB_INS=1                 !Old default
	    CONS_FOR_R_GRID=1.0_LDP
	    EXP_FOR_R_GRID=0.0_LDP
	  ELSE IF(VELTYPE .EQ. 4)THEN
!
! No parameters required
!
	  ELSE IF(VELTYPE .EQ. 5)THEN
	    WRITE(LUER,*)'Velocity law 5 not implemented in this version',
	1                 ' of CMFGEN'
	    STOP
	  ELSE IF(VELTYPE .EQ. 6)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,
	1           'Core velocity (km/s)')
	    CALL RD_STORE_DBLE(VPHOT,'VPHOT',L_TRUE,
	1           'Photospheric velocity (km/s)')
	    CALL RD_STORE_DBLE(SCL_HT,'SCL_HT',L_TRUE,
	1           'Scale Height (in R*) of photosphere')
	    CALL RD_STORE_DBLE(VINF1,'VINF1',L_TRUE,
	1            'Terminal velocity (km/s) if no 2nd comp.')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,
	1           'Speed of 1st Beta velocity Law')
	    CALL RD_STORE_DBLE(V_EPPS1,'EPPS1',L_TRUE,
	1           'Scale factor for 1s Beta velocity Law')
	    CALL RD_STORE_DBLE(VINF2,'VINF2',L_TRUE,
	1           'True terminal velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA2,'BETA2',L_TRUE,
	1           'Speed of 2nd Beta velocity Law')
	    CALL RD_STORE_DBLE(V_EPPS2,'EPPS2',L_TRUE,
	1           'Scale factor for 2nd Beta V law')
!
	    N_OB_INS=1          !Old default
	    N_IB_INS=2
	    CONS_FOR_R_GRID=-1.0_LDP
	    CALL RD_STORE_INT(N_OB_INS,'NBND_INS',L_FALSE,
	1           'Number of additional points to insert in radius grid at boundary')
	    CALL RD_STORE_DBLE(CONS_FOR_R_GRID,'C_R_GRID',L_FALSE,
	1           'Constant to allow improved shoice of R grid')
	    IF(CONS_FOR_R_GRID .GT. 0)THEN
	      CALL RD_STORE_DBLE(EXP_FOR_R_GRID,'E_R_GRID',L_TRUE,
	1           'Constant to allow improved shoice of R grid')
	    ELSE
	      CONS_FOR_R_GRID=1.0_LDP
	      EXP_FOR_R_GRID=0.0_LDP
	    END IF
!
! !Required by routines other than STARPCYG
!
	    VINF=VINF2	
	  ELSE IF(VELTYPE .EQ. 7)THEN
	    CALL RD_STORE_NCHAR(VEL_OPTION,'VEL_OPT',ITEN,L_TRUE,
	1                        'Velocity option: RVSIG_COL or deKOTER')
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,'Terminal velocity (km/s)')
	    VCORE=0.0_LDP		!Not used but initialized
	    VPHOT=0.0_LDP
	    SCL_HT=0.0_LDP
	    V_BETA1=0.0_LDP
	    V_EPPS1=1.0_LDP
	    VINF=VINF1
	    VINF2=VINF1		!i.e. no 2nd component
	    V_BETA2=1.0_LDP
	    V_EPPS2=1.0_LDP
	  ELSE IF(VELTYPE .EQ. 10)THEN
	    SN_MODEL=.TRUE.
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,'Initial velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,'Power of velocity Law')
	    RHO_ZERO=0.0_LDP
	    CALL RD_STORE_DBLE(RHO_ZERO,'RHO_ZERO',L_FALSE,'Initial density (gm/cm^3)')
	    IF(RHO_ZERO .EQ. 0.0_LDP)THEN
	      CALL RD_STORE_DBLE(RCUB_RHO_ZERO,'RCUB_RHO',L_TRUE,'r^3 . initial density 10^{-30} gm)')
	      RHO_ZERO=RCUB_RHO_ZERO/(RP**3)
	    END IF
	    RHO_ZERO=RHO_ZERO/ATOMIC_MASS_UNIT()
	    CALL RD_STORE_DBLE(N_RHO,'N_RHO',L_TRUE,'Density exponent (+ve)')
	    VINF=VCORE*(RMAX/RP)**V_BETA1
	  ELSE IF(VELTYPE .EQ. 11)THEN
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,'Terminal velocity (km/s)')
	    SN_MODEL=.TRUE.
	    SN_HYDRO_MODEL=.TRUE.
	  ELSE IF(VELTYPE .EQ. 12)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,'Initial velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,'Power of velocity Law')
	  ELSE IF(VELTYPE .EQ. 13 .OR. VELTYPE .EQ. 14)THEN
	    SN_MODEL=.TRUE.
	    SN_HYDRO_MODEL=.TRUE.
	  ELSE
	    WRITE(LUER,*)'Velocity law ',VELTYPE, ' not implemented',
	1                ' in this version of CMFGEN'
	    STOP
	  END IF
!
	  VAR_MDOT=.FALSE.
	  IF(VELTYPE .EQ. 7)THEN
	    CALL RD_STORE_LOG(VAR_MDOT,'VAR_MDOT',L_FALSE,'Variable mass-loss rate model?')
	    CALL RD_STORE_NCHAR(VAR_MDOT_FILE,'VM_FILE',ITEN,VAR_MDOT,
	1                        'File with density and clumping info')
	  END IF
!
	  IF(VAR_MDOT .OR. SN_MODEL)THEN
	    RMDOT=1.0E-20_LDP
	  ELSE
	    CALL RD_STORE_DBLE(RMDOT,'MDOT',L_TRUE,'Mass Loss rate (Msun/yr) ')
	  END IF
	  CALL RD_STORE_DBLE(LUM,'LSTAR',L_TRUE,'Stellar luminosity (Lsun)')
!
	  NC_PNT_SRCE=0.0_LDP
	  R_PNT_SRCE=0.0_LDP
	  TEFF_PNT_SRCE=0.0_LDP
	  PNT_SRCE_MOD=.FALSE.
	  CALL RD_STORE_LOG(PNT_SRCE_MOD,'PNT_SRCE',L_FALSE,'Include point source with hollow shell')
	  IF(PNT_SRCE_MOD)THEN
	    CALL RD_STORE_DBLE(R_PNT_SRCE,'R_PNT',L_TRUE,'Radius of point source in Rsun')
	    CALL RD_STORE_DBLE(TEFF_PNT_SRCE,'TEFF_PNT',L_TRUE,'Effective T of point source in 10^4K')
	    NC_PNT_SRCE=2
	  END IF
!
! TEFF and LOGG only need to be present if DO_HYDRO is TRUE. Values will still be read
! in when available if DO_HYDRO is FALSE.
!
	  TEFF=0.0_LDP; LOGG=0.0_LDP
	  IF(DO_HYDRO)THEN
	    CALL RD_STORE_DBLE(TEFF,'TEFF',DO_HYDRO,'Effective temperature (10^4 K)')
	    CALL RD_STORE_DBLE(LOGG,'LOGG',DO_HYDRO,'Log surface gravity (cgs units)')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA',DO_HYDRO,'Beta exponent for velocity law')
	    PRESSURE_VTURB=0.0_LDP
	    CALL RD_STORE_DBLE(PRESSURE_VTURB,'P_VTURB',L_FALSE,'Pressure trubulent velocity')
	  END IF
	  IF(TEFF .NE. 0.0_LDP .AND. .NOT. DO_HYDRO)THEN
	    WRITE(LUER,'(A)')' Possible Error in VADAT file.'
	    WRITE(LUER,'(A)')' You have set TEFF in VADAT but DO_HYDRO is false.'
            WRITE(LUER,'(A)')' TEFF will not be VALID unless hydro-iterations are performed.'
            WRITE(LUER,'(A)')' Remove TEFF if you are not going to do a hdyro iteration'
            WRITE(LUER,'(A)')' In this case TEFF is set by L, R and Mdot'
	    STOP
	  END IF

	  CALL RD_STORE_DBLE(STARS_MASS,'MASS',L_TRUE,'Stellar mass (Msun)')
!
! All clumping parameters are read in, even when CLUMPING is switched off.
!
	  CALL RD_STORE_LOG(DO_CLUMP_MODEL,'DO_CL',L_TRUE,
	1            'Calculate a model with clumping?')
	  CALL RD_STORE_NCHAR(CLUMP_LAW,'CL_LAW',ISIX,L_TRUE,
	1      'Which clumping law is being utilized?')
	  CALL SET_CASE_UP(CLUMP_LAW,IZERO,IZERO)
	  CALL RD_STORE_INT(N_CLUMP_PAR,'N_CL_PAR',L_TRUE,
	1          'Number of clumping parameters')
	  IF(N_CLUMP_PAR .GT. N_CLUMP_PAR_MAX)THEN
	    WRITE(LUER,*)'Error in CMFGEN'
	    WRITE(LUER,*)'N_CLUMP_PAR too large: N_CLUMP_PAR=',N_CLUMP_PAR
	    STOP
	  END IF
	  CLUMP_PAR(:)=0.0_LDP
	  DO I=1,N_CLUMP_PAR			!Should be less than 10
	    TEMP_CHAR='CL_PAR_'
	    WRITE(TEMP_CHAR(8:8),'(I1)')I
	    CALL RD_STORE_DBLE(CLUMP_PAR(I),TEMP_CHAR(1:8),L_TRUE,
	1             'Clumping parameters:')
	  END DO
!
	  REVISE_R_GRID=.FALSE.
	  CALL RD_STORE_LOG(REVISE_R_GRID,'REV_RGRID',L_FALSE,'Automatically revise R grid?')
!
! These calls read in parameters which automatically control the
! revision of the R grid automatically in SN models that have sharp
! ionization fronts.
!
	  JGREY_WITH_V_TERMS=.FALSE.
	  PURE_HUBBLE_FLOW=.FALSE.
	  COMP_GREY_LST_IT=.TRUE.
	  N_RG_PAR=0
	  TIME_SEQ_NO=0.0_LDP
	  IF(SN_MODEL)THEN
!
	    JGREY_WITH_V_TERMS=.TRUE.
	    CALL RD_STORE_LOG(JGREY_WITH_V_TERMS,'JG_W_V',L_FALSE,
	1            'Include V terms when computing Grey temperature?')
	    DO_CO_MOV_DDT=.FALSE.
	    CALL RD_STORE_LOG(DO_CO_MOV_DDT,'DO_DDT',L_FALSE,
	1            'Include comoving derivative in SE equations?')
	    CALL RD_STORE_DBLE(TIME_SEQ_NO,'TS_NO',DO_CO_MOV_DDT,
	1            'Time sequence number: 1 for inital model')
	    CALL RD_STORE_DBLE(SN_AGE_DAYS,'SN_AGE',DO_CO_MOV_DDT,'Age of SN in days')
	    CALL RD_STORE_LOG(PURE_HUBBLE_FLOW,'PURE_HUB',L_TRUE,
	1            'Force a pre hubble flow using age and radii of SN?')
	    CALL RD_STORE_LOG(INCL_RADIOACTIVE_DECAY,'INC_RAD_DECAYS',L_TRUE,
	1            'Allow for radiactive decays')
!
! Control parameters for handling non-thermal electrons.
!
	    CALL RD_STORE_LOG(TREAT_NON_THERMAL_ELECTRONS,'TRT_NON_TE',L_TRUE,
	1            'Treat non-thermal electrons')
	    SCL_NT_CROSEC=.FALSE.
	    SCL_NT_ION_CROSEC=.FALSE.
	    NT_OMIT_ION_SCALE=1.0E-03_LDP
	    NT_OMIT_LEV_SCALE=1.0E-04_LDP
	    NT_NKT=1000                 !Number of energy bins
	    NT_EMAX=1000.0_LDP            !eV
	    NT_EMIN=1.0_LDP               !ev
	    NON_THERMAL_IT_CNTRL=1
	    NT_SOURCE_TYPE='BELL_SHAPE'
	    READ_NON_THERM_SPEC=.FALSE.
	    IF(TREAT_NON_THERMAL_ELECTRONS)THEN
	      CALL RD_STORE_INT(NT_NKT,'NT_NKT',L_FALSE,'Number of energy bins')
	      CALL RD_STORE_DBLE(NT_EMIN,'NT_EMIN',L_FALSE,'Minimum energy of non-thermal electrons in eV')
	      CALL RD_STORE_DBLE(NT_EMAX,'NT_EMAX',L_FALSE,'Maximum energy of non-thermal electrons in eV')
	      CALL RD_STORE_LOG(SCL_NT_CROSEC,'SCL_NT_CROSEC',L_FALSE,
	1               'Scale the nonthermal excitation cross sections?')
	      CALL RD_STORE_LOG(SCL_NT_ION_CROSEC,'SCL_NT_ION_CROSEC',L_FALSE,
	1               'Scale the nonthermal ionization cross sections?')
	      CALL RD_STORE_INT(NON_THERMAL_IT_CNTRL,'NT_IT_CNTRL',L_FALSE,
	1               'Controls how often we update the non-thermal electron distiution')
	      CALL RD_STORE_DBLE(NT_OMIT_LEV_SCALE,'NT_OMIT_LEV_SCALE',L_FALSE,
	1               'Fractional populations below this level are excluded when computing the non-thermal electron spectrum')
	      CALL RD_STORE_DBLE(NT_OMIT_ION_SCALE,'NT_OMIT_ION_SCALE',L_FALSE,
	1               'Excludes ions with population NT_OMIT_SCALE_FRAC*(species of ion pop)')
	      CALL RD_STORE_CHAR(NT_SOURCE_TYPE,'NT_SOURCE',L_FALSE,'Non-thermal source type - INJECT_DIRAC, CONSTANT or BELL_SHAPE')
	      CALL RD_STORE_LOG(READ_NON_THERM_SPEC,'RD_NT_SPEC',L_FALSE,
	1               'Read non-thermal electron spectrum?')
	    ENDIF
!
	    ADD_DEC_NRG_SLOWLY=.FALSE.
	    CALL RD_STORE_LOG(ADD_DEC_NRG_SLOWLY,'GAMMA_SLOW',L_FALSE,
	1            'Add radioactivity decay energy slowly?')
	    IF(ADD_DEC_NRG_SLOWLY)THEN
	       CALL RD_STORE_DBLE(DEC_NRG_SCL_FAC_BEG,'DECNRG_SCLFAC_BEG',ADD_DEC_NRG_SLOWLY,
	1            'Initial Scale factor for adding decay energy')
            END IF
!
!       LUC: Parameters for SHOCK POWER
!
	    INC_SHOCK_POWER=.FALSE.
	    CALL RD_STORE_LOG(INC_SHOCK_POWER,'INC_SHOCK_POWER',L_FALSE,
	1             'Include the shock power from HYDRO simulation')
	    IF (INC_SHOCK_POWER) THEN
	      SCL_PWR_BY_FCL = .FALSE.
	      CALL RD_STORE_LOG(SCL_PWR_BY_FCL,'SCL_PWR_BY_FCL',L_TRUE,
	1             'Scale shock power by CLUMP_FAC')
	      ADD_SHOCK_POWER_SLOWLY = .FALSE.
	      CALL RD_STORE_LOG(ADD_SHOCK_POWER_SLOWLY,'ADD_SHOCK_POWER_SLOWLY',L_TRUE,
	1            'Introduce shock power slowly')
	     IF (ADD_SHOCK_POWER_SLOWLY) THEN
	       CALL RD_STORE_DBLE(SHOCK_POWER_FAC_BEG,'SHOCK_POWER_FAC_BEG',ADD_SHOCK_POWER_SLOWLY,
	1              'Initial Scale factor for adding shock power')
	     ENDIF
	     CALL RD_STORE_DBLE(PRESCRIBED_SHOCK_POWER,'PRESCRIBED_SHOCK_POWER',L_FALSE,'Luminosity for the shock')
	     CALL RD_STORE_DBLE(VLOC_SHOCK_POWER,'VLOC_SHOCK_POWER',L_FALSE,'Velocity location for shock power')
	     CALL RD_STORE_DBLE(DVLOC_SHOCK_POWER,'DVLOC_SHOCK_POWER',L_FALSE,'Gausian width for shock power')
	     LIN_INTERP_RD_SN_DATA=.FALSE.
	     CALL RD_STORE_LOG(LIN_INTERP_RD_SN_DATA,'LIN_RD_SN_DATA',L_FALSE,'Linear interpolaton in RD_SN_DATA?')
	   ENDIF
!
	    MINIMUM_ISO_POP=1.0E-20_LDP
	    CALL RD_STORE_DBLE(MINIMUM_ISO_POP,'MIN_ISO_POP',L_FALSE,'Minimum population for ant ISOTOPE')
!
	    CALL RD_STORE_LOG(COMP_GREY_LST_IT,'COMP_GREY_LST_IT',L_FALSE,'Compute grey solution on last iteration?')

	    CALL RD_STORE_NCHAR(SN_T_OPTION,'SN_T_OPT',ITEN,L_TRUE,
	1           'Method to get T with non-GRID option (USE_T_IN or USE_HYDRO)')
	    CALL RD_STORE_NCHAR(GAMRAY_TRANS,'GAMRAY_TRANS',ITEN,L_TRUE,
	1           'NONLOCAL or LOCAL gamma-ray enegry transport')
	    INSTANTANEOUS_ENERGY_DEPOSITION=.FALSE.
	    CALL RD_STORE_LOG(INSTANTANEOUS_ENERGY_DEPOSITION,'INS_DEP',L_FALSE,'Instantaneous gamma-ray energy deposition?')
!
	    N_IB_INS=2
	    N_OB_INS=3
	    RMAX_ON_RCORE=-1.0_LDP		!Implies use default.
	    CALL RD_STORE_INT(N_IB_INS,'N_IB_INS',L_FALSE,'# of points for fine grid at inner boundary')
	    CALL RD_STORE_INT(N_OB_INS,'N_OB_INS',L_FALSE,'# of points for fine grid at outer boundary')
	    CALL RD_STORE_DBLE(RMAX_ON_RCORE,'RMAX_ON_RCORE',L_FALSE,'RMAX/RCORE for SN if shrinking radius')
	  END IF	
	  DO_FULL_REL_OBS=.FALSE.
	  DO_FULL_REL_CMF=.FALSE.
	  CALL RD_STORE_LOG(DO_FULL_REL_OBS,'REL_OBS',L_FALSE,
	1          'Include all reltivistic terms in observer''s frame solution')
	  CALL RD_STORE_LOG(DO_FULL_REL_CMF,'REL_CMF',L_FALSE,
	1        'Include all reltivistic terms in CMF solution for observed intensity')
!
! These two options are installed only to retain consistency with earlier versions.
!
	  USE_OLD_MF_SCALING=.FALSE.
	  USE_OLD_MF_OUTPUT=.FALSE.
	  CALL RD_STORE_LOG(USE_OLD_MF_SCALING,'OLD_MFS',L_FALSE,
	1        'Use old mass scaling when reading SN_HYDRO_DATA')
	  CALL RD_STORE_LOG(USE_OLD_MF_OUTPUT,'OLD_MFO',L_FALSE,
	1        'Use old mass scaling when reading SN_HYDRO_DATA')
!
	  CALL RD_STORE_NCHAR(SOL_ABUND_REF_SET,'SOL_ABUND_SCL',ITEN,L_FALSE,'Solar abundance reference set')
!
! Read in the un-normalized fractional abundances.
!
	  DO ISPEC=1,NUM_SPECIES
	    TMP_KEY=TRIM(SPECIES(ISPEC))//'/X'
	    TMP_STRING=TRIM(SPECIES(ISPEC))//
	1            '/X fractional abundance by number (un-normalized)'
	    CALL RD_STORE_DBLE(AT_ABUND(ISPEC),TMP_KEY,SPECIES_PRES(ISPEC),
	1            TMP_STRING)
	  END DO
	  WRITE(LUSCR,'()')
!
	  CALL RD_STORE_LOG(RD_CONT_FREQ,'RD_CF_FILE',L_TRUE,
	1            'Read in continuum frequencies from file')
	  CALL RD_STORE_DBLE(MIN_CONT_FREQ,'MIN_CF',L_TRUE,
	1            'Minimum continuum frequency if calculating NU')
	  CALL RD_STORE_DBLE(MAX_CONT_FREQ,'MAX_CF',L_TRUE,
	1            'Maximum continuum frequency if calculating NU')
	  CALL RD_STORE_DBLE(SMALL_FREQ_RAT,'FRAC_SP',L_TRUE,
	1            'Fractional spacing for small frequencies')
	  CALL RD_STORE_DBLE(BIG_FREQ_AMP,'AMP_FAC',L_TRUE,
	1            'Amplification factor for large frequency ranges')
	  CALL RD_STORE_DBLE(dFREQ_bf_MAX,'MAX_BF',L_TRUE,
	1            'Maximum frequency spacing close to bf edge')
!
! Installed to allow earlier frequency grids to be used. 0 uses the
! the latest default grid.
!
	  FREQ_GRID_OPTION=0
	  CALL RD_STORE_INT(FREQ_GRID_OPTION,'FR_GRID',L_FALSE,
	1            'Which method to compute frequency grid?')
!
	  CALL RD_STORE_LOG(DO_LEV_DISSOLUTION,'DO_DIS',L_TRUE,
	1            'Allow for level dissolution of upper levels?')
	  CALL RD_STORE_DBLE(dV_LEV_DIS,'dV_LEV',L_TRUE,
	1             'Spacing (in km/s) on low side of bf edge for'//
	1             ' level dissolution')
	  CALL RD_STORE_DBLE(AMP_DIS,'AMP_DIS',L_TRUE,
	1            'Amplification factor on low side bf edge')
	  CALL RD_STORE_DBLE(MIN_FREQ_LEV_DIS,'MIN_DIS',L_TRUE,
	1            'Minimum frequency for level dissolution')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(COMPUTE_ALL_CROSS,'CROSS',L_TRUE,
	1            'Compute all photoionization cross-sections?')
	  CALL RD_STORE_DBLE(DELV_CONT,'V_CROSS',L_TRUE,
	1            'Max. vel. sep. (km/s) between evaluations of all'//
	1            '  phot. cross-sections?')
	  CALL RD_STORE_LOG(DIE_AS_LINE,'DIE_AS_LINE',L_TRUE,
	1        'Treat the dielectronic transitions as individual lines?')
	  CALL RD_STORE_DBLE(VSM_DIE_KMS,'VSM_DIE',L_TRUE,
	1        'Velocity (km/s) for smoothing dielectronic transitions')
	  CALL RD_STORE_DBLE(SIG_GAU_KMS,'SIG_GAU_KMS',L_TRUE,
	1        'Sigma of Gaussian used to smooth photoionization data')
	  FRAC_SIG_GAU=0.25_LDP
	  CALL RD_STORE_DBLE(FRAC_SIG_GAU,'FRAC_SIG',L_FALSE,
	1        'Fractional spacing a across smoothing Gauusian (use 0.25)')
	  CUT_ACCURACY=0.02_LDP
	  CALL RD_STORE_DBLE(CUT_ACCURACY,'CUT_ACC',L_FALSE,
	1        'Accuracy to retain data when omitting data points to save space (use 0.02)')
	  ABOVE_EDGE=.TRUE.
	  CALL RD_STORE_LOG(ABOVE_EDGE,'ABV_EDGE',L_FALSE,
	1        'Use only data above edge when smoothing (TRUE)')
!
	  CALL RD_STORE_DBLE(EXT_LINE_VAR,'EXT_LINE_VAR',L_TRUE,
	1        'Extent of line variation zone (V/INF) beyond'//
	1        'the resonancze zone')
          IF(EXT_LINE_VAR .LT. 0.0_LDP .OR. EXT_LINE_VAR .GT. 2.0_LDP)THEN
	    WRITE(LUER,*)'Error in CMFGEN --- invalid range for EXT_LINE_VAR'
	    STOP
	  END IF
	  NEW_VAR_STORAGE_METHOD=.TRUE.
	  CALL RD_STORE_LOG(NEW_VAR_STORAGE_METHOD,'NEW_VST_METH',L_FALSE,
	1         'Use new storage method for lines (different storage for each line)?')
	!
! NB: An ideal vale for ZNET_VAR_LIMIT is probably 0.01 or 0.001. If
! ZNET_VAR_LIMIT is zero, all depths will be included in the linearization,
! independent of ZNET. A very large value of ZNET (i.e. 10^4), will imply
! an interation on the NET_RATES, with no linearization.
!
	  CALL RD_STORE_DBLE(ZNET_VAR_LIMIT,'ZNET_VAR_LIM',L_TRUE,
	1            'Include lines in full varaition when '//
	1            ' ABS(ZNET-1) > ZNET_VAR_LIM')
	  CALL RD_STORE_LOG(WEAK_WITH_NET,'WNET',L_TRUE,'Use Lambda iteration for weak lines?')
	  CALL RD_STORE_DBLE(WEAK_LINE_LIMIT,'WK_LIM',L_TRUE,'Maximum opacity ratio for weak lines (0.01)?')
!
	  USE_WEAK_TAU_LIM=.FALSE.; WEAK_TAU_LINE_LIMIT=0.01_LDP
	  CALL RD_STORE_LOG(USE_WEAK_TAU_LIM,'WK_TAU',L_FALSE,'Use TAU(SOB) to decide weak lines?')
	  CALL RD_STORE_DBLE(WEAK_TAU_LINE_LIMIT,'WK_TAU_LIM',USE_WEAK_TAU_LIM,
	1                      'Maximum TAU for weak lines (0.01)?')
!
	  CALL RD_STORE_LOG(DIF,'DIF',L_TRUE,
	1            'Use Diffusion approximation at inner boundary ?')
	  IF(DIF)THEN
	    INNER_BND_METH='DIFFUSION'
	    CALL RD_STORE_CHAR(INNER_BND_METH,'IB_METH',L_FALSE,
	1           'Inner boundary method (DIFFUSION, HOLLOW, or ZERO_FLUX)')
	    IF(INNER_BND_METH .NE. 'DIFFUSION')THEN
	      WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES'
	      WRITE(LUER,*)'Inconsistent inner boundary condition'
	      WRITE(LUER,*)'DIF option and IB_METH are inconsistent'
	      STOP
	    END IF
	  ELSE
	    CALL RD_STORE_CHAR(INNER_BND_METH,'IB_METH',L_TRUE,
	1           'Inner boundary method (DIFUSION, HOLLOW, or ZERO_FLUX)')
	  END IF
!
! Centralized INNER_BND_METH validation: downstream code selects the
! inner boundary by string compare in several places, so a trailing
! blank or case difference would silently take the wrong branch.
! The whitelist covers every method the codebase implements; the
! Knop MOM solver applies its own stricter subset.
!
	  CALL SET_CASE_UP(INNER_BND_METH,IZERO,IZERO)
	  IF(INNER_BND_METH .NE. 'DIFFUSION' .AND.
	1     INNER_BND_METH .NE. 'ZERO_FLUX' .AND.
	2     INNER_BND_METH .NE. 'PNT_SRCE' .AND.
	3     INNER_BND_METH .NE. 'SCHUSTER' .AND.
	4     INNER_BND_METH .NE. 'HOLLOW' .AND.
	5     INNER_BND_METH .NE. 'JEQB')THEN
	    WRITE(LUER,*)'Invalid INNER_BND_METH=',INNER_BND_METH
	    STOP
	  END IF
	  IB_STAB_FACTOR=0.1_LDP
	  IF(INNER_BND_METH .EQ. 'DIFFUSION')IB_STAB_FACTOR=0.0_LDP
	  CALL RD_STORE_DBLE(IB_STAB_FACTOR,'IB_STAB',L_FALSE,'Inner boundary stabilization factor')
	  OUTER_BND_METH='HONJ'
	  CALL RD_STORE_CHAR(OUTER_BND_METH,'OB_METH',L_FALSE,'Outer boundary method (HONJ or HALF_MOM)')
!
	  CALL RD_STORE_LOG(RD_COHERENT_ES,'COH_ES',L_TRUE,
	1            'Assume coherent electron scattering? ')
	  CALL RD_STORE_LOG(USE_OLDJ_FOR_ES,'OLD_J',L_TRUE,
	1            'Use old file to provide initial estimate of J_ES?')
	  COHERENT_ES=RD_COHERENT_ES
	  CALL RD_STORE_LOG(MIXED_ES_VAR,'MIX_COH',L_TRUE,
	1            'Mix coherent/non-coherent e.s. in linearization?')
	  CALL RD_STORE_DBLE(ES_VAR_FAC,'ES_FAC',L_TRUE,
	1            'Fractional proximity of RJ and RJ_ES for coherent'//
	1            ' variation')
!
	  LTE_MODEL=.FALSE.
	  CALL RD_STORE_LOG(LTE_MODEL,'LTE_MOD',L_FALSE,
	1            'Force populations to LTE and iterate T only?')
!
	  CALL RD_STORE_NCHAR(METHOD,'METHOD',ISIX,L_TRUE,
	1         'Which method for continuum tau'//
	1         ' loglog, loglin, linear or zero ?')
	  LUM_FROM_ETA_METHOD='LINMON'
	  CALL RD_STORE_NCHAR(LUM_FROM_ETA_METHOD,'LUM_METH',ISIX,L_FALSE,
	1         'Which method for computing L from ETA '//
	1         ' loglog, loglin, linear or zero ?')
	  CALL RD_STORE_NCHAR(N_TYPE,'N_TYPE',ISIX,L_TRUE,
	1         'Method for to handle N for MOM_J_CMF -- '//
	1         'N_ON_J, MIXED, or G_ONLY')
!
	  J_CHK_OPTION='FORM_VAL'
	  H_CHK_OPTION='MAX_VAL'
	  XM_CHK_OPTION='NONE'
	  CALL RD_STORE_CHAR(J_CHK_OPTION,'CHK_J',L_FALSE,'FORM_VAL, ABS_VAL, or NONE')
	  CALL RD_STORE_CHAR(H_CHK_OPTION,'CHK_H',L_FALSE,'MAX_VAL, AV_VAL, or NONE')
	  CALL RD_STORE_CHAR(XM_CHK_OPTION,'CHK_XM',L_FALSE,'SET_POS, or NONE')
!
	  CALL RD_STORE_NCHAR(FG_SOL_OPTIONS,'FG_OPT',ITEN,L_TRUE,
	1         'Solution options for FG_J_CMF: DIFF/INS and INT/INS')
	  CALL RD_STORE_DBLE(DELV_FRAC_FG,'VFRAC_FG',L_TRUE,
	1         'Maximum velocity spacing (Doppler units) in FG_J_CMF_V10')
	  CALL RD_STORE_DBLE(DELV_FRAC_MOM,'VFRAC_MOM',L_TRUE,
	1         'Maximum velocity spacing (Doppler units) in MOM_J_CMF_V10')
!
	  CALL RD_STORE_LOG(RDTHK_CONT,'THK_CONT',L_TRUE,'Use thick boundary condition for continuum ? ')
	  RD_OUT_BC_TYPE=1
	  OUT_BC_TYPE=1
	  OUT_BC_PARAM_ONE=0.299794_LDP
	  CALL RD_STORE_INT(RD_OUT_BC_TYPE,'OBC_TYPE',L_FALSE,'Outer boundary condition type: 1=def=old')
	  CALL RD_STORE_DBLE(OUT_BC_PARAM_ONE,'BC_PAR1',L_FALSE,'Frequency to switch to new BC')
	  REXT_FAC=0.0_LDP
	  CALL RD_STORE_DBLE(REXT_FAC,'REXT_FAC',L_FALSE,'Factor ot extend R grid by for thick continuum')
!
	  INCL_INCID_RAD=.FALSE.
	  CALL RD_STORE_LOG(INCL_INCID_RAD,'INCID_RAD',L_FALSE,
	1           'Include incident radiation for plane-parellel mod with V?')
	  CALL RD_STORE_LOG(TRAPFORJ,'TRAP_J',L_TRUE,
	1           'Use trapazoidal weights to compute J? ')
!
! This section of the code should work with old VADAT  files when using a Doppler profile of
! fixed width (the default).
!
	  WRITE(LUSCR,'()')
	  FIX_DOP=.TRUE.
	  AMASS_DOP=1.0_LDP
	  CALL RD_STORE_LOG(FIX_DOP,'FIX_DOP',L_FALSE,
	1      'Use the same turbulent velocity for all species?')
	  CALL RD_STORE_DBLE(TDOP,'TDOP',L_TRUE,
	1      'Temperature to be used in Doppler profile (10^4K)')
	  IF(FIX_DOP)THEN
	    CALL RD_STORE_DBLE(AMASS_DOP,'AMASS_DOP',L_TRUE,
	1      'Atomic mass to be used in Doppler profile (amu''s)')
	    CALL RD_STORE_DBLE(VTURB,'VTURB',L_TRUE,
	1      'Turbulent velocity to be used in Doppler profile (km/s)')
	    VTURB_MIN=VTURB; VTURB_MAX=VTURB; VTURB_POW=1; VTURB_LAW='POW'
	    GLOBAL_LINE_PROF='DOP_FIX'
	  ELSE
	    VTURB_LAW='POW'; VTURB_POW=1
	    CALL RD_STORE_NCHAR(VTURB_LAW,'VTURB_LAW',ITEN,L_FALSE,
	1      'Law describing variation of turbulent velocity with V')
	    CALL RD_STORE_DBLE(VTURB_POW,'VTURB_POW',L_FALSE,
	1      'Power describing variation of turbulent velocity with V')
	    IF(VTURB_LAW .EQ. 'TRUNC_POW')THEN
	      CALL RD_STORE_DBLE(VTURB_VEND,'VTURB_VEND',L_TRUE,
	1        'V(km/s) at which lineare increase in Vturb ends')
	    END IF
	    CALL RD_STORE_DBLE(VTURB_MIN,'VTURB_MIN',L_TRUE,
	1      'Minimum turbulent velocity for Doppler profile (km/s)')
	    CALL RD_STORE_DBLE(VTURB_MAX,'VTURB_MAX',L_TRUE,
	1      'Maximum turbulent velocity for Doppler profile (km/s)')
	    VTURB=VTURB_MIN
!
	    WRITE(LUSCR,'()')
	    CALL RD_STORE_NCHAR(GLOBAL_LINE_PROF,'GLOBAL_PROF',ITEN,L_TRUE,
	1        'Global switch for intrinsic line absorption profile')
	    CALL SET_CASE_UP(GLOBAL_LINE_PROF,IZERO,IZERO)
	    IF( GLOBAL_LINE_PROF .NE. 'NONE' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOP_FIX' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOP_SPEC' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOPPLER' .AND.
	1       GLOBAL_LINE_PROF .NE. 'LIST' .AND.
	1       GLOBAL_LINE_PROF .NE. 'LIST_VGT' .AND.
	1       GLOBAL_LINE_PROF .NE. 'VOIGT' .AND.
	1       GLOBAL_LINE_PROF .NE. 'HZ_STARK')THEN
	      WRITE(LUER,*)'Invalid GLOBAL_LINE_PROF parameter'
	      STOP
	    END IF
	    CALL RD_STORE_LOG(SET_PROF_LIMS_BY_OPACITY,'OPAC_LIMS',L_TRUE,
	1           'Set prof limits by line to cont. ratio?')
	    CALL RD_STORE_DBLE(DOP_PROF_LIMIT,'DOP_LIM',L_TRUE,
	1           'Edge limits for Doppler line profile')
	    CALL RD_STORE_DBLE(VOIGT_PROF_LIMIT,'VOIGT_LIM',L_TRUE,
	1           'Edge limits for Voigt line profile')
	  END IF
!
! The following options do not apply to a Doppler profile.
!
	  V_PROF_LIMIT=3000.0_LDP
	  MAX_PROF_ED=1.0E+16_LDP
	  NORM_PROFILE=.TRUE.
	  CALL RD_STORE_DBLE(MAX_PROF_ED,'MAX_PROF_ED',L_FALSE,
	1           'Maximum electron density for Stark profile computation')
	  CALL RD_STORE_DBLE(V_PROF_LIMIT,'V_PROF_LIM',L_FALSE,
	1           'One-sided profile limit for Stark profiles (km/s)')
	  CALL RD_STORE_LOG(NORM_PROFILE,'NORM_PROF',L_FALSE,
	1           'When true, profiles are normalized to have unit area.')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(MAX_DOP,'MAX_DOP',L_TRUE,
	1      'Maximum half-width of resonance zone (in Doppler widths)')
	  CALL RD_STORE_DBLE(FRAC_DOP,'FRAC_DOP',L_TRUE,
	1      'Spacing in resonance zone (in Doppler widths)')
	  CALL RD_STORE_DBLE(dV_CMF_PROF,'dV_CMF_PROF',L_TRUE,
	1      'Spacing across cmf profile (in km/s)')
	  CALL RD_STORE_DBLE(dV_CMF_WING,'dV_CMF_WING',L_TRUE,
	1      'Spacing across e.s. wings of cmf profile(in km/s)')
	  CALL RD_STORE_DBLE(ES_WING_EXT,'ES_WING_EXT',L_TRUE,
	1      'Extent of BLUE e.s. wings from resonance core (in km/s)')
	  CALL RD_STORE_DBLE(R_CMF_WING_EXT,'R_CMF_WING_EXT',L_TRUE,
	1      'Extent of RED e.s. wings from RESONANCE core (in Vinf)')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(OBS_PRO_EXT_RAT,'OBS_EXT_RAT',L_TRUE,
	1      'Half width of profile in Vinf.')
	  CALL RD_STORE_DBLE(dV_OBS_PROF,'dV_OBS_PROF',L_TRUE,
	1      'Spacing across observed profile (in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_WING,'dV_OBS_WING',L_TRUE,
	1      'Spacing across e.s. wings of observed profile(in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_BIG,'dV_OBS_BIG',L_TRUE,
	1      'Frequency spacing between lines (in km/s)')
!
	  WRITE(LUSCR,'()')
	  USE_FIXED_J=.FALSE.
	  CALL RD_STORE_LOG(USE_FIXED_J,'USE_FIXED_J',L_FALSE,
	1           'Use previously computed J to evaluate ALL rates?')
	  IF(USE_FIXED_J .AND. .NOT. RD_LAMBDA)THEN
	    WRITE(LUER,'(A)')' Warning: RD_LAMBDA must be TRUE when USE_FIXED_J is TRUE.'
	    WRITE(LUER,'(A)')' Please change the setting in IN_ITS which is read after each iteration'
	    STOP
	  END IF
	  CALL RD_STORE_LOG(FLUX_CAL_ONLY,'FLUX_CAL_ONLY',L_TRUE,
	1           'Compute the observers frame flux only ?')
	  CALL RD_STORE_LOG(EXTEND_FRM_SOL,'EXT_FRM_SOL',L_TRUE,
	1           'Extrapolate the formal solution to larger radii?')
	  CALL RD_STORE_LOG(INSERT_FREQ_FRM_SOL,'INS_F_FRM_SOL',L_TRUE,
	1           'Insert extra frequencies for formal solution?')
	  CALL RD_STORE_NCHAR(CMF_FORM_OPTIONS,'FRM_OPT',ITEN,L_TRUE,
	1           'Solution options for CMF_FORM_SOL')
	  CALL RD_STORE_LOG(DO_SOBOLEV_LINES,'DO_SOB_LINES',L_TRUE,
	1        'Compute Sobolev rates and EWs for flux calculation?')
	  CALL RD_STORE_LOG(SOB_FREQ_IN_OBS,'SOB_FREQ_IN_OBS',L_TRUE,
	1        ' Allow for SOB & CMF lines in defining observers'//
	1        ' frequencies?')
!
	  VERBOSE_OUTPUT=.FALSE.
	  CALL RD_STORE_LOG(VERBOSE_OUTPUT,'VERBOSE_OUT',L_FALSE,
	1        'Switch on enhanced diagnostic output')
	  CALL SET_VERBOSE_INFO(VERBOSE_OUTPUT)
	  WRITE_RATES=.FALSE.
	  CALL RD_STORE_LOG(WRITE_RATES,'WRITE_RATES',L_FALSE,'Write out NETRATE, TOTRATE, etc')
	  WRITE_JH=.FALSE.
	  IF(SN_MODEL)WRITE_JH=.TRUE.
	  CALL RD_STORE_LOG(WRITE_JH,'WRITE_JH',L_FALSE,'Write out JH_AT_CURRENT_TIME')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_NCHAR(GLOBAL_LINE_SWITCH,'GLOBAL_LINE',ISIX,L_TRUE,
	1            'Global switch to indicate handeling of line')
	  CALL SET_CASE_UP(GLOBAL_LINE_SWITCH,IZERO,IZERO)
	  IF( GLOBAL_LINE_SWITCH(1:3) .NE. 'SOB' .AND.
	1       GLOBAL_LINE_SWITCH(1:3) .NE. 'CMF' .AND.
	1       GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE' .AND.
	1       GLOBAL_LINE_SWITCH(1:5) .NE. 'BLANK')THEN
	    WRITE(LUER,*)'Invalid GLOBAL_LINE SWITCH parameter'
	    STOP
	  END IF
	  CALL RD_STORE_LOG(SET_TRANS_TYPE_BY_LAM,'LAM_SET',L_TRUE,
	1         'Set long wavelengths to SOBOLEV approximation')
	  CALL RD_STORE_DBLE(FLUX_CAL_LAM_BEG,'F_LAM_BEG',L_TRUE,
	1         'Inital wavelength (A) for blanketed flux calculation')
	  CALL RD_STORE_DBLE(FLUX_CAL_LAM_END,'F_LAM_END',L_TRUE,
	1         'Final wavelength (A) for blanketed flux calculation')
	  CALL RD_STORE_DBLE(GF_CUT,'GF_CUT',L_TRUE,
	1          'gf value to omit transitions')
	  CALL RD_STORE_DBLE(AT_NO_GF_CUT,'AT_CUT',L_TRUE,
	1          'Only omit transitions if AT_NO >= AT_CUT')
	  CALL RD_STORE_INT(GF_LEV_CUT,'GF_LEV_CUT',L_TRUE,
	1          'Level above whit transitions omitted if gf < GF_CUT')
	  CALL RD_STORE_INT(MIN_NUM_TRANS,'MIN_TRANS',L_TRUE,
	1          'Minimum number of transitions from each level')
!
	  CALL RD_STORE_LOG(THK_LINE,'THK_LINE',L_TRUE,
	1           'Use thick boundary condition for lines?')
	  CALL RD_STORE_LOG(CHECK_LINE_OPAC,'CHK_L_POS',L_TRUE,
	1      'Ensure Line opacity is positive (SOB & CMF modes only)?')
	  CALL RD_STORE_NCHAR(NEG_OPAC_OPTION,'NEG_OPAC_OPT',ITEN,L_TRUE,
	1            'Method for negative opacities in BLANKETING mode')
	  CALL SET_CASE_UP(NEG_OPAC_OPTION,IZERO,IZERO)
	  IF(NEG_OPAC_OPTION .NE. 'SRCE_CHK' .AND.
	1                           NEG_OPAC_OPTION .NE. 'ESEC_CHK')THEN
	     WRITE(LUER,*)'Error in CMFGEN_SUB'
	     WRITE(LUER,*)'Invalid NEG_OPAC_OPTION'
	     WRITE(LUER,*)'Valid options are SRCE_CHK and ESEC_CHK'
	     STOP
	  END IF
!
! Non-monotonic velocity-field options (dev7).  Optional keywords; when
! FALSE (default) the original code path is executed unchanged.
!
!   DO_NONMONO_SOB : Dessart, Audit & Hillier (2015, arXiv:1503.05463)
!                    Sobolev treatment (unclamped sigma in SOB routines).
!   DO_NONMONO_KNOP: Knop et al. (2009) CMF formal solver
!                    (FG_J_CMF_NONMONO_V1 replaces FG_J_CMF_V13).
!
! The two modes are mutually exclusive.
!
	  DO_NONMONO_SOB=.FALSE.
	  CALL RD_STORE_LOG(DO_NONMONO_SOB,'DO_NONMONO_SOB',L_FALSE,
	1     'Use Dessart+2015 Sobolev treatment for non-monotonic V?')
!
	  DO_NONMONO_KNOP=.FALSE.
	  CALL RD_STORE_LOG(DO_NONMONO_KNOP,'DO_NONMONO_KNOP',L_FALSE,
	1     'Use Knop et al. 2009 non-monotonic formal solver?')
!
	  DO_NONLOC_SOB=.FALSE.
	  CALL RD_STORE_LOG(DO_NONLOC_SOB,'DO_NONLOC_SOB',L_FALSE,
	1     'Use Rybicki-Hummer non-local Sobolev coupling?')
!
	  IF(DO_NONLOC_SOB .AND. .NOT. DO_NONMONO_SOB)THEN
	    WRITE(LUER,*)'Warning: DO_NONLOC_SOB=T but DO_NONMONO_SOB=F'
	    WRITE(LUER,*)'The non-local coupling uses SIGMA_TRUE and only'
	    WRITE(LUER,*)'runs in the non-monotonic SOB path; it is ignored.'
	    DO_NONLOC_SOB=.FALSE.
	  END IF
!
	  IF(DO_NONMONO_KNOP .AND. DO_NONMONO_SOB)THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'DO_NONMONO_KNOP and DO_NONMONO_SOB are both .TRUE.'
	    WRITE(LUER,*)'These are mutually exclusive non-monotonic modes.'
	    WRITE(LUER,*)'Set exactly one to .TRUE., or both to .FALSE. for'
	    WRITE(LUER,*)'the original monotonic CMFGEN code path.'
	    STOP
	  END IF
!
! GET_POPS_AT_PREV_TIME_STEP_* maps populations with LOG(V) and
! interpolates R(V), both valid only for a positive single-valued
! Hubble-law velocity.  PURE_HUBBLE_FLOW likewise replaces the input
! field with V proportional to R.  Reject these expansion controls
! rather than discarding the non-monotonic structure.
!
	  IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1     (DO_CO_MOV_DDT .OR. PURE_HUBBLE_FLOW))THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'DO_DDT and PURE_HUB are incompatible with'
	    WRITE(LUER,*)'non-monotonic velocity fields.  The time-sequence'
	    WRITE(LUER,*)'restart interpolates on LOG(V), and PURE_HUB'
	    WRITE(LUER,*)'forces a monotonic V proportional to R.'
	    STOP
	  END IF
!
! Hydro adjustment can rewrite the prescribed velocity structure and
! must stay blocked under the non-monotonic solvers.  Automatic
! R-grid revision is now permitted: the only method whose R(V) map
! is multi-valued on a shocked field (VEL_SPACE) is rejected inside
! ADJUST_R_GRID_V4 itself, the interpolated-sigma STOP that assumed
! sigma >= -1 is bypassed under these flags, and
! RECOMPUTE_SIGMA_TRUE refreshes SIGMA_TRUE after every revision.
! TAU_SPACE, FULL_R and OB_ONLY interpolate LOG(V) on the monotone
! R coordinate, which is single-valued and non-monotonic safe.
!
! edited by Devin, 13-Oct-2026: was (DO_HYDRO .OR. REVISE_R_GRID).
!
	  IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1     DO_HYDRO)THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'DO_HYDRO is incompatible with non-monotonic'
	    WRITE(LUER,*)'velocity fields: it can rewrite the prescribed'
	    WRITE(LUER,*)'V(R).'
	    STOP
	  END IF
!
! N_LAM_RWARM (default 4): forced LAMBDA iterations at process start
! when DO_NONMONO_KNOP=T.  The Gauss-Seidel I_NEXT arrays in
! FG_J_CMF_NONMONO_V1 are process-local and are re-zeroed on every
! requeue; this forces LAMBDA for the first N_LAM_RWARM great
! iterations of each process so a BA step is never taken against a
! stale (zeroed) I_NEXT.  Set 0 to disable.
!
	  N_LAM_RWARM=4
	  CALL RD_STORE_INT(N_LAM_RWARM,'N_LAM_RWARM',L_FALSE,
	1     'Forced LAMBDA its at process start for KNOP I_NEXT re-init?')
	  IF(N_LAM_RWARM .LT. 0)N_LAM_RWARM=0
!
! MOM_TRIAL_MIN (default 30): frozen-iteration count that arms the
! marginal-release BA trial in SOLVE_FOR_POPS.  <=0 disables; the
! trial also requires DO_BA_RES_ACC so a bad step rolls back.
!
	  MOM_TRIAL_MIN=30
	  CALL RD_STORE_INT(MOM_TRIAL_MIN,'MOM_TRIAL_MIN',L_FALSE,
	1     'KNOP: frozen iters before a resacc-gated BA trial?')
	  IF(MOM_TRIAL_MIN .LT. 0)MOM_TRIAL_MIN=0
!
! GS_FREEZE_LIM (default 0 = disabled): hard stop on consecutive
! forced-LAMBDA iterations while the KNOP lagged field stays
! immature.  The submit script applies its own KNOP_STALL_MAX on
! the GS_STATE file, so a stalled run fails even with this key off.
!
	  GS_FREEZE_LIM=0
	  CALL RD_STORE_INT(GS_FREEZE_LIM,'GS_FREEZE_LIM',L_FALSE,
	1     'Stop model after this many consecutive GS-frozen its?')
	  IF(GS_FREEZE_LIM .LT. 0)GS_FREEZE_LIM=0
!
! GS_TRIAL_LIM (default 15): accepted resacc-gated BA trials the KNOP
! marginal-release path may take between raw-mature sweeps.  Accepted
! trials refund n_frozen (reset to MOM_TRIAL_MIN); rejected ones do
! not.  <=0 disables the refund.
!
	  GS_TRIAL_LIM=15
	  CALL RD_STORE_INT(GS_TRIAL_LIM,'GS_TRIAL_LIM',L_FALSE,
	1     'KNOP: accepted trial releases before freeze veto?')
	  IF(GS_TRIAL_LIM .LT. 0)GS_TRIAL_LIM=0
!
! GS_REL_STABLE (default 1): consecutive sweeps passing the maturity
! gate before the freeze releases the population update.  Guards
! against releasing on a single marginal sweep of a churning field.
!
	  GS_REL_STABLE=1
	  CALL RD_STORE_INT(GS_REL_STABLE,'GS_REL_STABLE',L_FALSE,
	1     'KNOP: mature sweeps required before freeze release?')
	  IF(GS_REL_STABLE .LT. 1)GS_REL_STABLE=1
!
! GS_REL_N (default 1): post-release LAMBDA iterations held at the
! near-mature correction cap while the field re-lags to the moved
! populations.
!
	  GS_REL_N=1
	  CALL RD_STORE_INT(GS_REL_N,'GS_REL_N',L_FALSE,
	1     'KNOP: damped LAMBDA its after a freeze release?')
	  IF(GS_REL_N .LT. 1)GS_REL_N=1
!
! GS_HOLD_POPS (default FALSE): on a GS-freeze iteration, roll the
! forced-LAMBDA population update back after SOLVEBA returns, so the
! great iteration converges the lagged radiation field at fixed
! populations.  The freeze already forces LAMBDA + fixed T; without
! the rollback the populations still move, which keeps the GS fixed
! point moving and defeats the nested iteration.  The resacc-gated
! MOM_TRIAL release path is excluded -- the trial exists to take a
! real BA step.  Default F preserves the previous behavior (a capped
! LAMBDA update still commits during the freeze).
!
	  GS_HOLD_POPS=.FALSE.
	  CALL RD_STORE_LOG(GS_HOLD_POPS,'GS_HOLD_POPS',L_FALSE,
	1     'Hold populations fixed during KNOP GS-freeze sweeps?')
	  GS_POPS_HELD=.FALSE.
	  GS_ATTEMPT_APP=0.0_LDP
!
! DO_KNOP_REFRESH (default DO_NONMONO_KNOP, or DO_NONMONO_SOB with
! USE_LAM_ES=F): per-ray re-solve of the lagged compressing-branch
! coupling after each frequency sweep (Knop Paper II Sect. 3.2
! ordering), replacing the one-bin-per-iteration lag that cannot
! mature on strongly non-monotonic models.  The SOB moment path needs
! it for the same reason: without a within-iteration refresh the
! lagged J_NEXT coupling injects a stale, dimmer upwind field and the
! emergent luminosity diverges secularly (m15 mom_probe, it441-453).
! Set FALSE for A/B comparison with the original lagged scheme.
!
	  DO_KNOP_REFRESH=DO_NONMONO_KNOP .OR.
	1    (DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES)
	  CALL RD_STORE_LOG(DO_KNOP_REFRESH,'DO_KNOP_REFRESH',L_FALSE,
	1     'Per-ray converge of KNOP lagged I_NEXT coupling?')
	  CALL SET_FG_KNOP_REFRESH(DO_KNOP_REFRESH)
!
! MOM refresh toggles (Oct-2026).  MOM_RF_GS_DESC stays FALSE
! (unvalidated).  MOM_RF_KEEPNEG and MOM_GATE_BULK default TRUE:
! storing repaired (non-fixed-point) values reseeds the lagged
! read with noise, and the strict any-negJ gate cannot clear a
! quasi-steady negJ residual (froze BA ~200 its on m15_2p82d;
! see cmfgen_obsidian/history/convergence_issues.md section 10).
!
	  MOM_RF_GS_DESC=.FALSE.
	  CALL RD_STORE_LOG(MOM_RF_GS_DESC,'MOM_RF_GS_DESC',L_FALSE,
	1     'Live J_NEXT read on descending MOM refresh passes?')
	  MOM_RF_KEEPNEG=.TRUE.
	  CALL RD_STORE_LOG(MOM_RF_KEEPNEG,'MOM_RF_KEEPNEG',L_FALSE,
	1     'Keep prior J_NEXT iterate at MOM-repaired cells?')
	  CALL SET_MOM_RF_OPTS(MOM_RF_GS_DESC,MOM_RF_KEEPNEG)
!
	  MOM_GATE_BULK=.TRUE.
	  CALL RD_STORE_LOG(MOM_GATE_BULK,'MOM_GATE_BULK',L_FALSE,
	1     'MOM maturity from bad-cell fraction, not max-norm?')
!
! MOM_RF_BULK_EXIT (default TRUE): exit the MOM refresh once the
! per-pass drift meets the same bulk bounds the gate applies to
! the sweep, instead of burning passes until the max-norm stall
! fires.  The retained-vs-sweep delta gate still applies to the
! converged field.
!
	  MOM_RF_BULK_EXIT=.TRUE.
	  CALL RD_STORE_LOG(MOM_RF_BULK_EXIT,'MOM_RF_BULK_EXIT',L_FALSE,
	1     'Exit MOM refresh on bulk-drift criteria?')
	  CALL SET_MOM_RF_BULK(MOM_RF_BULK_EXIT)
!
! MOM_RF_ACTIVE (default FALSE, experimental): skip refresh solves
! at frequencies whose residual and upwind read slots stayed below
! tolerance on the previous pass.
!
	  MOM_RF_ACTIVE=.FALSE.
	  CALL RD_STORE_LOG(MOM_RF_ACTIVE,'MOM_RF_ACTIVE',L_FALSE,
	1     'Active-set front-tracking in MOM refresh?')
	  CALL SET_MOM_RF_ACT(MOM_RF_ACTIVE)
!
! PRINT_KNOP_DEBUG (default FALSE): verbose Knop-solver diagnostics.
! Gates the O(NP*NI) I_NEXT depth-mapping table in OUTGEN; the mapping
! check itself always runs and reports anomalies regardless.
!
	  PRINT_KNOP_DEBUG=.FALSE.
	  CALL RD_STORE_LOG(PRINT_KNOP_DEBUG,'PRINT_KNOP_DEBUG',L_FALSE,
	1     'Print verbose KNOP solver diagnostics to OUTGEN?')
	  CALL SET_FG_KNOP_DEBUG(PRINT_KNOP_DEBUG)
!
! DO_NONMONO_FREQ_MATCH (default FALSE): rematch the lagged upwind
! frequency index at compressing points to the bin containing the
! photon's characteristic foot (per-cell drift in bins), and widen
! the frequency-derivative stencil to match.  Targets the
! oscillatory-eigenvalue instability identified by Ceulemans et
! al. 2025 when the per-cell Doppler shift exceeds the frequency
! spacing.  Only acts where the Knop solver already switches to
! the lagged redder upwind; expanding cells keep the implicit
! same-sweep coupling.
!
	  DO_NONMONO_FREQ_MATCH=.FALSE.
	  CALL RD_STORE_LOG(DO_NONMONO_FREQ_MATCH,
	1     'DO_NONMONO_FREQ_MATCH',L_FALSE,
	1     'Rematched upwind freq at compressing cells (KNOP only)?')
	  CALL SET_FG_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)
	  CALL SET_MOM_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)
!
! MOM refresh Anderson acceleration (default TRUE, depth 3).  For
! the affine-in-lagged-field refresh map AA(m) is the GMRES
! residual minimization on (I-A)x=b; expose it so a run can
! disable AA without a rebuild.
!
	  MOM_RF_AA=.TRUE.
	  CALL RD_STORE_LOG(MOM_RF_AA,'MOM_RF_AA',L_FALSE,
	1     'Anderson acceleration on MOM refresh?')
	  CALL SET_MOM_RF_AA(MOM_RF_AA)
!
! MOM_EPS_LIM (default TRUE): donor-cell positivity limiter on the
! lagged EPS midpoint transfer -- the lagged term may not debit a
! row past its positive RHS balance (removes the dominant negJ
! source on advection-dominated compressing grids).  MOM_RF_DIR
! (default 'AUTO'): refresh pass ordering; AUTO runs only the
! dominant coupling direction once one side holds >=90% of rows,
! and suppresses the stale-Jacobi in-sweep store write at
! compressing rows on a compression-dominant grid.  RF_EN passes
! DO_KNOP_REFRESH so suppression never runs without a refresh to
! own the store.
!
	  MOM_EPS_LIM=.TRUE.
	  CALL RD_STORE_LOG(MOM_EPS_LIM,'MOM_EPS_LIM',L_FALSE,
	1     'Positivity-limited lagged EPS transfer in MOM solve?')
	  MOM_RF_DIR='AUTO'
	  CALL RD_STORE_CHAR(MOM_RF_DIR,'MOM_RF_DIR',L_FALSE,
	1     'MOM refresh pass order: AUTO, ALT, ASC or DESC')
!
! MOM_RF_OWNW (default TRUE): inside a refresh pass, rows whose
! upwind read was stale this direction (compressing on ascending,
! expanding on descending) do not write the lagged store -- each row
! only takes fresh-read GS steps.  Removes the Jacobi re-injection
! floor of the alternating map; the fixed point is unchanged.
!
	  MOM_RF_OWNW=.TRUE.
	  CALL RD_STORE_LOG(MOM_RF_OWNW,'MOM_RF_OWNW',L_FALSE,
	1     'Pass-direction-owned lagged-store writes in MOM refresh?')
	  CALL SET_MOM_RF_DIR_OPT(MOM_RF_DIR,MOM_EPS_LIM,
	1     DO_KNOP_REFRESH,MOM_RF_OWNW)
!
! Refresh tuning knobs, previously compile-time parameters.  The
! defaults match the built-in values, so an unchanged VADAT is
! bit-identical in behavior; the setters clamp to safe ranges.
! MOM_AA_M is the Anderson history depth (memory per pair is
! ND_NM*NCF*8 bytes).  JNEXT_RELAX/INEXT_RELAX are the
! lagged-store under-relaxation factors.
!
	  KN_M_TOL_V=2.0E-3_LDP
	  CALL RD_STORE_DBLE(KN_M_TOL_V,'KN_M_TOL',L_FALSE,
	1     'MOM refresh max-norm convergence tolerance')
	  KN_M_MAXPASS_V=60
	  CALL RD_STORE_INT(KN_M_MAXPASS_V,'KN_M_MAXPASS',L_FALSE,
	1     'Max MOM refresh passes')
	  JNEXT_RELAX_V=0.5_LDP
	  CALL RD_STORE_DBLE(JNEXT_RELAX_V,'JNEXT_RELAX',L_FALSE,
	1     'MOM lagged-store under-relaxation factor')
	  MOM_AA_M_V=3
	  CALL RD_STORE_INT(MOM_AA_M_V,'MOM_AA_M',L_FALSE,
	1     'Anderson history depth for MOM refresh')
!
! Outer (call-level) AA on the committed lagged field: the residual
! is the successive committed-store difference, so the accelerated
! map includes the J<->Eddington coefficient drift the inner
! frozen-coefficient AA cannot see (the retained-field delta limit
! cycle).  MOM_AA_OUT_M is the pair count; memory per pair is
! 2*ND_NM*NCF*8 bytes.
!
	  MOM_AA_OUT=.TRUE.
	  CALL RD_STORE_LOG(MOM_AA_OUT,'MOM_AA_OUT',L_FALSE,
	1     'Call-level Anderson acceleration on MOM lagged field?')
	  MOM_AA_OUT_M_V=8
	  CALL RD_STORE_INT(MOM_AA_OUT_M_V,'MOM_AA_OUT_M',L_FALSE,
	1     'Anderson history depth for the call-level MOM map')
	  CALL SET_MOM_KN_OPTS(KN_M_TOL_V,KN_M_MAXPASS_V,
	1     JNEXT_RELAX_V,MOM_AA_M_V,MOM_AA_OUT,MOM_AA_OUT_M_V)
!
	  KN_RAY_TOL_V=2.0E-3_LDP
	  CALL RD_STORE_DBLE(KN_RAY_TOL_V,'KN_RAY_TOL',L_FALSE,
	1     'FG refresh per-ray convergence tolerance')
	  KN_RAY_MAXPASS_V=60
	  CALL RD_STORE_INT(KN_RAY_MAXPASS_V,'KN_RAY_MAXPASS',L_FALSE,
	1     'Max FG refresh passes per ray')
	  INEXT_RELAX_V=0.5_LDP
	  CALL RD_STORE_DBLE(INEXT_RELAX_V,'INEXT_RELAX',L_FALSE,
	1     'FG lagged-store under-relaxation factor')
	  CALL SET_FG_KN_OPTS(KN_RAY_TOL_V,KN_RAY_MAXPASS_V,
	1     INEXT_RELAX_V)
!
! KN_RAY_ORDERED (default TRUE): ordered characteristic march in
! the FG refresh (Knop Paper II ordering).  FALSE restores the
! alternating-direction passes.
!
	  KN_RAY_ORDERED=.TRUE.
	  CALL RD_STORE_LOG(KN_RAY_ORDERED,'KN_RAY_ORDERED',L_FALSE,
	1     'Ordered (Knop) march in the FG refresh?')
	  CALL SET_FG_ORD_SWEEP(KN_RAY_ORDERED)
!
! ------------------------------------------------------------------
! JFNK (Jacobian-Free Newton-Krylov) population solver.
! Replaces the SOLVEBA_V13 correction solve inside SOLVE_FOR_POPS
! when DO_JFNK=T.  Admitted under DO_NONMONO_KNOP or DO_NONMONO_SOB;
! admission additionally requires USE_LAM_ES=T (already mandatory for
! SOB) and, under KNOP only, a mature lagged field (GS_FREEZE=F).
!
! JFNK_COUPLED=T evaluates the residual with a fresh formal solve
! at each perturbed population state (JFNK_RAD_SWEEP_V1), so the
! FD Jacobian-vector product probes the full dJ/dn coupling --
! the Arramy et al. (2024, A&A 690, A12) formulation.  With
! JFNK_COUPLED=F the residual uses stored (fixed-J) radiative
! rates plus the CMAT_RAD linearization -- the same fixed-J system
! SOLVEBA_V13 solves directly, run through the globalized Newton
! machinery instead.
!
! JFNK_COUPLED is KNOP-only: JFNK_RAD_SWEEP_V1 replays the
! blanketing continuum chain and does not recompute the SOB line
! rates (SUB_SOB_LINE_V3), so under DO_NONMONO_SOB a coupled eval
! would solve a residual different from the equations the outer
! iteration assembles.  SOB therefore admits only the frozen-J tier.
!
! editted by Devin (V7-JFNK port, 02-Oct-2026)
! editted by Devin (DO_JFNK shared SOB/KNOP switch, 07-Oct-2026)
! ------------------------------------------------------------------
!
	  USE_JFNK_SOLVER=.FALSE.
	  CALL RD_STORE_LOG(USE_JFNK_SOLVER,'DO_JFNK',L_FALSE,
	1     'Use JFNK population solver for non-monotonic models?')
!
! USE_JFNK: pre-SOB name for the same switch, retained so existing
! KNOP VADAT files keep working.  Either key T arms the solver; an
! explicit F under one name cannot veto the other because the reader
! cannot distinguish an explicit F from an absent key.
!
	  USE_JFNK_OLD=.FALSE.
	  CALL RD_STORE_LOG(USE_JFNK_OLD,'USE_JFNK',L_FALSE,
	1     'Deprecated alias for DO_JFNK')
	  IF(USE_JFNK_OLD)THEN
	    WRITE(LUER,*)'Warning: USE_JFNK is deprecated; use DO_JFNK.'
	    USE_JFNK_SOLVER=.TRUE.
	  END IF
!
	  IF(USE_JFNK_SOLVER .AND.
	1     .NOT. (DO_NONMONO_KNOP .OR. DO_NONMONO_SOB))THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'DO_JFNK requires DO_NONMONO_KNOP or'
	    WRITE(LUER,*)'DO_NONMONO_SOB=.TRUE.'
	    STOP
	  END IF
!
	  JFNK_COUPLED=.FALSE.
	  CALL RD_STORE_LOG(JFNK_COUPLED,'JFNK_COUPLED',L_FALSE,
	1     'Evaluate JFNK residual with fresh formal solve (coupled)?')
!
	  IF(USE_JFNK_SOLVER .AND. JFNK_COUPLED .AND. DO_NONMONO_SOB)THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'JFNK_COUPLED=T is not available under'
	    WRITE(LUER,*)'DO_NONMONO_SOB.  The eval sweep callback'
	    WRITE(LUER,*)'(JFNK_RAD_SWEEP_V1) replays the KNOP continuum'
	    WRITE(LUER,*)'chain and does not recompute the SOB line rates'
	    WRITE(LUER,*)'(SUB_SOB_LINE_V3), so the coupled residual would'
	    WRITE(LUER,*)'not match the equations the outer iteration'
	    WRITE(LUER,*)'assembles.  Set JFNK_COUPLED=F: SOB admits only'
	    WRITE(LUER,*)'the frozen-J residual tier.'
	    STOP
	  END IF
!
	  JFNK_ANY_ITER=.FALSE.
	  CALL RD_STORE_LOG(JFNK_ANY_ITER,'JFNK_ANYITER',L_FALSE,
	1     'Debug: fire JFNK even on LAMBDA-designated calls?')
!
! Logarithmic unknown and pseudo-time live in the solver module so
! the Newton loop and the residual share them without a new entry
! in CONTROL_VARIABLE_MOD.  Defaults are the stable formulation
! (cmfgen_obsidian/solvers/JFNK_SOLVER.md).  Absent keys keep those defaults.
!
	  jfnk_log_unknown=.TRUE.
	  CALL RD_STORE_LOG(jfnk_log_unknown,'JFNK_LOGX',L_FALSE,
	1     'JFNK unknown is ln(n/n0) rather than a fraction?')
	  jfnk_ptc_active=.TRUE.
	  CALL RD_STORE_LOG(jfnk_ptc_active,'JFNK_PTC',L_FALSE,
	1     'Pseudo-time continuation on the JFNK residual?')
	  jfnk_ptc_dtau=1.0D0
	  jfnk_ptc_scale_ready=.FALSE.
!
! JFNK_IN_EVAL / JFNK_EVAL_BROKEN are internal state (no VADAT keys);
! re-init per read so a residual-eval abort can never leak a stale
! in-eval or broken-eval flag.
!
	  JFNK_IN_EVAL=.FALSE.
	  JFNK_EVAL_BROKEN=.FALSE.
!
	  JFNK_MAX_NEWTON=20
	  CALL RD_STORE_INT(JFNK_MAX_NEWTON,'JFNK_NITS',L_FALSE,
	1     'Max Newton iterations for JFNK solve')
!
	  JFNK_M_GMRES=10
	  CALL RD_STORE_INT(JFNK_M_GMRES,'JFNK_M',L_FALSE,
	1     'GMRES restart dimension for JFNK solve')
!
	  JFNK_MAX_GMRES=0
	  CALL RD_STORE_INT(JFNK_MAX_GMRES,'JFNK_MAXGMRES',L_FALSE,
	1     'Max GMRES its per Newton step (0=use restart dim m)')
!
	  JFNK_GTOL=1.0E-3_LDP
	  CALL RD_STORE_DBLE(JFNK_GTOL,'JFNK_GTOL',L_FALSE,
	1     'Relative GMRES tolerance (fraction of |F|)')
!
	  JFNK_TOL_REL=1.0E-4_LDP
	  CALL RD_STORE_DBLE(JFNK_TOL_REL,'JFNK_TOL',L_FALSE,
	1     'Relative Newton tolerance (target |F|/|F_0|)')
!
	  JFNK_EPS_JVP=SQRT(EPSILON(1.0_LDP))
	  CALL RD_STORE_DBLE(JFNK_EPS_JVP,'JFNK_EPS',L_FALSE,
	1     'Finite-difference step for Jacobian-vector products')
!
	  JFNK_DEL_INIT=-1.0_LDP
	  CALL RD_STORE_DBLE(JFNK_DEL_INIT,'JFNK_DEL',L_FALSE,
	1     'Initial trust-region radius (<0 = auto)')
!
	  JFNK_MNDL=1.0E-20_LDP
	  CALL RD_STORE_DBLE(JFNK_MNDL,'JFNK_MNDL',L_FALSE,
	1     'Minimum trust-region radius')
!
	  JFNK_MXDL=1.0E+20_LDP
	  CALL RD_STORE_DBLE(JFNK_MXDL,'JFNK_MXDL',L_FALSE,
	1     'Maximum trust-region radius')
!
	  JFNK_FIXED_T=.TRUE.
	  CALL RD_STORE_LOG(JFNK_FIXED_T,'JFNK_FIXT',L_FALSE,
	1     'Hold T fixed in JFNK? (aligned to FIX_T after read)')
!
	  JFNK_VERBOSE=.TRUE.
	  CALL RD_STORE_LOG(JFNK_VERBOSE,'JFNK_VERB',L_FALSE,
	1     'Print JFNK per-iteration diagnostics?')
!
	  JFNK_DIAGNOSTIC=.FALSE.
	  CALL RD_STORE_LOG(JFNK_DIAGNOSTIC,'JFNK_DIAG',L_FALSE,
	1     'Run one-shot Newton/GMRES diagnostic routines on first step?')
!
	  JFNK_WARMUP=6
	  CALL RD_STORE_INT(JFNK_WARMUP,'JFNK_WARMUP',L_FALSE,
	1     'Number of GS/LAMBDA warm-up iterations before JFNK (0=cold start)')
!
! MAX_JFNK: max population scale factor after a JFNK step (like MAX_LIN).
!
	  MAX_JFNK_COR=10.0_LDP
	  CALL RD_STORE_DBLE(MAX_JFNK_COR,'MAX_JFNK',L_FALSE,
	1     'Max population scale factor after a JFNK Newton step')
	  IF(MAX_JFNK_COR .LT. 1.01_LDP)THEN
	    WRITE(LUER,*)'Warning in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'MAX_JFNK=',MAX_JFNK_COR,' is too small; using 1.01'
	    MAX_JFNK_COR=1.01_LDP
	  END IF
!
	  JFNK_APPLY_BA_RIDGE=.FALSE.
	  JFNK_BA_RIDGE_FAC=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_BA_RIDGE_FAC,'JFNK_RIDGE',L_FALSE,
	1     'Seed/local ridge for JFNK BA (0=equilibrated LU first)')
	  IF(JFNK_BA_RIDGE_FAC .LT. 0.0_LDP)JFNK_BA_RIDGE_FAC=0.0_LDP
	  JFNK_RIDGE_PIV=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_RIDGE_PIV,'JFNK_RIDGE_PIV',L_FALSE,
	1     'Escalate ridge if min/max LU pivot below this (0=off)')
!
	  JFNK_PRECOND_REFINE=.TRUE.
	  CALL RD_STORE_LOG(JFNK_PRECOND_REFINE,'JFNK_REFINE',L_FALSE,
	1     'Iterative refinement (DGERFS) after JFNK precond solve?')
!
	  JFNK_LEFT_PRECOND=.TRUE.
	  CALL RD_STORE_LOG(JFNK_LEFT_PRECOND,'JFNK_LEFTPC',L_FALSE,
	1     'Use left GMRES preconditioning (Arramy 2024)?')
!
	  JFNK_PC_REFRESH_FRAC=0.8_LDP
	  CALL RD_STORE_DBLE(JFNK_PC_REFRESH_FRAC,'JFNK_PC_RF',L_FALSE,
	1     'Refresh LU cache when GMRES its > frac*m')
	  IF(JFNK_PC_REFRESH_FRAC .LT. 0.1_LDP)
	1     JFNK_PC_REFRESH_FRAC=0.1_LDP
	  IF(JFNK_PC_REFRESH_FRAC .GT. 1.0_LDP)
	1     JFNK_PC_REFRESH_FRAC=1.0_LDP
!
	  JFNK_LAM_AUTO_THRESH=50.0_LDP
	  CALL RD_STORE_DBLE(JFNK_LAM_AUTO_THRESH,'JFNK_LAM_THRESH',L_FALSE,
	1     'DO_LAM_AUTO Lambda->JFNK threshold (percent MAXCH)')
	  IF(JFNK_LAM_AUTO_THRESH .LT. 10.0_LDP)
	1     JFNK_LAM_AUTO_THRESH=10.0_LDP
	  IF(JFNK_LAM_AUTO_THRESH .GT. 100.0_LDP)
	1     JFNK_LAM_AUTO_THRESH=100.0_LDP
!
	  JFNK_COR_GATE=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_COR_GATE,'JFNK_COR_GATE',L_FALSE,
	1     'JFNK admission: last demanded MAXCH (%) below this to fire'//
	2     ' (0=off)')
!
	  JFNK_FRZ_GATE=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_FRZ_GATE,'JFNK_FRZ_GATE',L_FALSE,
	1     'Frozen-J JFNK admission gate (0=use JFNK_COR_GATE)')
!
	  JFNK_GATE_NREQ=1
	  CALL RD_STORE_INT(JFNK_GATE_NREQ,'JFNK_GATE_NREQ',L_FALSE,
	1     'Consecutive sub-gate solves required before JFNK fires')
	  IF(JFNK_GATE_NREQ .LT. 1)JFNK_GATE_NREQ=1
!
	  JFNK_DEEPDIAG=.FALSE.
	  CALL RD_STORE_LOG(JFNK_DEEPDIAG,'JFNK_DEEPDIAG',L_FALSE,
	1     'Run ~12-eval diagnostic suite on first JFNK activation?')
!
	  JFNK_MAX_EVALS=0
	  CALL RD_STORE_INT(JFNK_MAX_EVALS,'JFNK_MAX_EVALS',L_FALSE,
	1     'Max residual evaluations per JFNK activation (0=off)')
	  IF(JFNK_MAX_EVALS .LT. 0)JFNK_MAX_EVALS=0
!
	  JFNK_PRECOND_TYPE='LU'
	  CALL RD_STORE_CHAR(JFNK_PRECOND_TYPE,'JFNK_PRECOND',L_FALSE,
	1     'Preconditioner type: LU, LSQ, SCHUR or JACOBI')
	  CALL SET_CASE_UP(JFNK_PRECOND_TYPE,IZERO,IZERO)
	  IF(JFNK_PRECOND_TYPE(1:3) .NE. 'LU ' .AND.
	1    JFNK_PRECOND_TYPE(1:6) .NE. 'JACOBI' .AND.
	2    JFNK_PRECOND_TYPE(1:5) .NE. 'SCHUR' .AND.
	3    JFNK_PRECOND_TYPE(1:3) .NE. 'LSQ')THEN
	    JFNK_PRECOND_TYPE='LU'
	  END IF
!
! JFNK_LSQ_LAM (default 1D-2): relative Tikhonov damping for the
! LSQ preconditioner.  The LU solve (Dr*J*Dc)^{-1} returns huge
! corrections along near-null singular directions of the
! population Jacobian, which GMRES cannot use; LSQ instead solves
! min |J_eq y - v'|^2 + lam^2 |y|^2 per depth block via the
! normal equations (J_eq^T J_eq + lam^2 I) y = J_eq^T v', so
! unsolvable directions are damped to ~0 rather than amplified.
! lam is relative to max|J_eq| (which is ~1 after DGEEQU
! equilibration).  Verified against the algebraic B^T w apply in
! cmfgen_jfnk_unit_test (rel diff < 1D-13 vs the direct solve).
!
	  JFNK_LSQ_LAM=1.0D-2
	  CALL RD_STORE_DBLE(JFNK_LSQ_LAM,'JFNK_LSQ_LAM',L_FALSE,
	1     'Relative Tikhonov damping for JFNK_PRECOND=LSQ')
	  IF(JFNK_LSQ_LAM .LE. 0.0_LDP)JFNK_LSQ_LAM=1.0D-2
!
! JFNK_REACH_GATE (default 0 = off): admission precheck that skips
! the Newton solve when the preconditioned step direction carries no
! residual content.  z=M^-1 F and Jz=BA*z are evaluated before
! Newton; REACH=(2<F,Jz>-<Jz,Jz>)/<F,F> is the fraction of |F|^2 a
! linear step can remove.  On the SOB frozen-J residual the pinned
! radiative rows leave REACH~0, so every GMRES/JVP eval is wasted;
! the gate converts such admissions to ~one matvec.  Gate 0.1-0.5 is
! conservative: solvable test systems reach ~1.
!
	  JFNK_REACH_GATE=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_REACH_GATE,'JFNK_REACH_GATE',L_FALSE,
	1     'Skip Newton when residual reachable fraction below this')
	  IF(JFNK_REACH_GATE .LT. 0.0_LDP)JFNK_REACH_GATE=0.0_LDP
	  IF(JFNK_REACH_GATE .GT. 1.0_LDP)JFNK_REACH_GATE=1.0_LDP
!
! JFNK_REARM_EVERY (default 10): admitted calls skipped before a
! stagnated JFNK is retried.  Each failure cycle doubles the wait
! (cap x8), so a persistently failing Newton costs few retries while
! a recovering residual is re-attempted soon.  JFNK is never
! permanently disabled -- Lambda iterations carry the model forward
! between retries.
!
	  JFNK_REARM_EVERY=10
	  CALL RD_STORE_INT(JFNK_REARM_EVERY,'JFNK_REARM_EVERY',L_FALSE,
	1     'Admitted calls skipped before retrying a stagnated JFNK')
	  IF(JFNK_REARM_EVERY .LT. 1)JFNK_REARM_EVERY=10
!
! JFNK_ROW_DIAG (default 10 rows; 0 = off): at each admission, print
! the largest |F| rows labeled by ion/level/depth together with their
! collisional vs stored-radiation split, Jacobian row norm, and (when
! the reach precheck ran) the Jz_i the linear step assigns them.  The
! table distinguishes radiation-lag rows (|rad|>>|coll|, Jz~0) from
! population-driven ones, which decides whether a population-space
! residual (frozen or coupled) can ever reach them.
!
	  JFNK_ROW_DIAG=10
	  CALL RD_STORE_INT(JFNK_ROW_DIAG,'JFNK_ROW_DIAG',L_FALSE,
	1     'Top residual rows printed at each JFNK admission (0=off)')
	  IF(JFNK_ROW_DIAG .LT. 0)JFNK_ROW_DIAG=10
	  IF(JFNK_ROW_DIAG .GT. 50)JFNK_ROW_DIAG=50
!
! Freeze population-dependent occupation-probability inputs during a
! frozen-J activation.  This changes the residual operator, so keep it
! opt-in and report its status in the JFNK context.
	  JFNK_LAG_FACTORS=.FALSE.
	  CALL RD_STORE_LOG(JFNK_LAG_FACTORS,'JFNK_LAG_FACTORS',L_FALSE,
	1     'Hold POPION and dissolution factors during frozen-J evals?')
!
! FIX_T_AUTO damping ladder: the VARFIXT depth range becomes a
! damped-T domain instead of a hard mask.  RE is the T-row residual
! norm at admission; the rung advances when RE falls below
! JFNK_TAU_REL*RE0 after JFNK_TAU_NLOCK commits and reverts when the
! committed log-T step exceeds JFNK_TAU_EPS.
!
	  JFNK_TAU_A0=0.05_LDP
	  CALL RD_STORE_DBLE(JFNK_TAU_A0,'JFNK_TAU_A0',L_FALSE,
	1     'JFNK_TAUTO first damping rung')
	  JFNK_TAU_REL=0.3_LDP
	  CALL RD_STORE_DBLE(JFNK_TAU_REL,'JFNK_TAU_REL',L_FALSE,
	1     'JFNK_TAUTO RE/RE0 fraction to release a rung')
	  JFNK_TAU_EPS=0.1_LDP
	  CALL RD_STORE_DBLE(JFNK_TAU_EPS,'JFNK_TAU_EPS',L_FALSE,
	1     'JFNK_TAUTO max committed log-T step before revert')
	  JFNK_TAU_DMAX=0.5_LDP
	  CALL RD_STORE_DBLE(JFNK_TAU_DMAX,'JFNK_TAU_DMAX',L_FALSE,
	1     'JFNK_TAUTO hard log-T step cap per Newton step')
	  JFNK_TAU_NLOCK=2
	  CALL RD_STORE_INT(JFNK_TAU_NLOCK,'JFNK_TAU_NLOCK',L_FALSE,
	1     'JFNK_TAUTO commits before release allowed')
	  JFNK_TAU_PCW=.TRUE.
	  CALL RD_STORE_LOG(JFNK_TAU_PCW,'JFNK_TAU_PCW',L_FALSE,
	1     'Down-weight T rows in the JFNK preconditioned system')
	  IF(JFNK_TAU_A0 .LE. 0.0_LDP)JFNK_TAU_A0=0.05_LDP
	  IF(JFNK_TAU_REL .LE. 0.0_LDP)JFNK_TAU_REL=0.3_LDP
	  IF(JFNK_TAU_EPS .LE. 0.0_LDP)JFNK_TAU_EPS=0.1_LDP
	  IF(JFNK_TAU_DMAX .LE. 0.0_LDP)JFNK_TAU_DMAX=0.5_LDP
	  IF(JFNK_TAU_NLOCK .LT. 0)JFNK_TAU_NLOCK=2
!
	  JFNK_RSCALE_MODE='EUCLIDEAN'
	  CALL RD_STORE_CHAR(JFNK_RSCALE_MODE,'JFNK_RSCALE',L_FALSE,
	1     'JFNK residual scaling: EUCLIDEAN, SQRT_STEQ, DIAG_RDG, DIAG')
	  CALL SET_CASE_UP(JFNK_RSCALE_MODE,IZERO,IZERO)
	  IF(JFNK_RSCALE_MODE(1:9) .NE. 'EUCLIDEAN' .AND.
	1    JFNK_RSCALE_MODE(1:8) .NE. 'SQRT_STE' .AND.
	2    JFNK_RSCALE_MODE(1:8) .NE. 'DIAG_RDG' .AND.
	3    JFNK_RSCALE_MODE(1:4) .NE. 'DIAG')THEN
	    JFNK_RSCALE_MODE='EUCLIDEAN'
	  END IF
	  JFNK_MAXCH_CNV=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_MAXCH_CNV,'JFNK_MXCH_TOL',L_FALSE,
	1     'Converge JFNK when max |F_i|/|BA_ii| below this (0=off)')
	  JFNK_POP_FLOOR=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_POP_FLOOR,'JFNK_POP_FL',L_FALSE,
	1     'Freeze JFNK rows with n_i < JFNK_POP_FL*ion total (0=off)')
!
	  JFNK_MAX_RESTART=8
	  CALL RD_STORE_INT(JFNK_MAX_RESTART,'JFNK_MAX_RESTART',L_FALSE,
	1     'Max JFNK preconditioner rebuilds per great iteration')
	  IF(JFNK_MAX_RESTART .LT. 1)JFNK_MAX_RESTART=1
	  IF(JFNK_MAX_RESTART .GT. 50)JFNK_MAX_RESTART=50
!
	  JFNK_LINE_SEARCH=.TRUE.
	  CALL RD_STORE_LOG(JFNK_LINE_SEARCH,'JFNK_LS',L_FALSE,
	1     'Use Armijo line search on Newton direction?')
!
	  JFNK_LS_MAX=20
	  CALL RD_STORE_INT(JFNK_LS_MAX,'JFNK_LS_MAX',L_FALSE,
	1     'Max backtracking steps for line search')
!
	  JFNK_LS_BETA=0.5_LDP
	  CALL RD_STORE_DBLE(JFNK_LS_BETA,'JFNK_LS_BETA',L_FALSE,
	1     'Backtracking factor for line search')
!
	  JFNK_LS_C=1.0E-4_LDP
	  CALL RD_STORE_DBLE(JFNK_LS_C,'JFNK_LS_C',L_FALSE,
	1     'Armijo sufficient-decrease parameter')
!
	  JFNK_CAUCHY=.TRUE.
	  CALL RD_STORE_LOG(JFNK_CAUCHY,'JFNK_CAUCHY',L_FALSE,
	1     'Use Cauchy (M^-1 F) fallback direction?')
!
	  JFNK_USE_BANDED=.FALSE.
	  CALL RD_STORE_LOG(JFNK_USE_BANDED,'JFNK_USE_BANDED',L_FALSE,
	1     'Use banded BA_rad and banded LU preconditioner?')
!
	  JFNK_SELECTIVE_JVAR=.FALSE.
	  CALL RD_STORE_LOG(JFNK_SELECTIVE_JVAR,'JFNK_SELECTIVE_JVAR',
	1     L_FALSE,'Enable selective J-variation (cross-depth coupling)?')
!
! JFNK_ROW_TINY: rows of a per-depth preconditioner block with
! max|A_ij| < JFNK_ROW_TINY*blockmax are replaced by the identity row
! before equilibration.
!
	  JFNK_ROW_TINY=1.0D-30
	  CALL RD_STORE_DBLE(JFNK_ROW_TINY,'JFNK_ROW_TINY',L_FALSE,
	1     'Near-null row threshold for JFNK preconditioner blocks')
!
! KN_ORD_AITKEN (default FALSE): guarded Aitken delta-squared
! extrapolation on the ordered-march field iterate at pass
! boundaries -- accelerates the marginally contractive cells at
! strong compressions that otherwise stall the march short of the
! bulk-maturity exit.  Applies to the ordered march only.
!
	  KN_ORD_AITKEN=.FALSE.
	  CALL RD_STORE_LOG(KN_ORD_AITKEN,'KN_ORD_AITKEN',L_FALSE,
	1     'Aitken-accelerate the ordered FG refresh march?')
	  CALL SET_FG_AITKEN(KN_ORD_AITKEN)
!
! Residual-accepted BA update (default FALSE).  The BA correction
! is committed only when the SE residual after the previous BA
! step did not grow by more than BA_RES_GROW; per-depth
! corrections above BA_TRUST are zeroed before scaling.  Cold
! start is unaffected: LAMBDA iterations are not gated, and the
! first genuine BA step only arms the tracker.
!
	  DO_BA_RES_ACC=.FALSE.
	  CALL RD_STORE_LOG(DO_BA_RES_ACC,'DO_BA_RES_ACC',L_FALSE,
	1     'Accept BA steps only when SE residual decreases?')
	  BA_RES_GROW=1.5_LDP
	  CALL RD_STORE_DBLE(BA_RES_GROW,'BA_RES_GROW',L_FALSE,
	1     'Max residual growth factor before BA step is rejected')
	  IF(BA_RES_GROW .LT. 1.0_LDP)BA_RES_GROW=1.0_LDP
	  BA_RES_HARD=1.0E+02_LDP
	  CALL RD_STORE_DBLE(BA_RES_HARD,'BA_RES_HARD',L_FALSE,
	1     'Reject BA step if residual exceeds this x running min')
	  IF(BA_RES_HARD .LT. 1.0_LDP)BA_RES_HARD=1.0_LDP
	  BA_TRUST=1.0E+04_LDP
	  CALL RD_STORE_DBLE(BA_TRUST,'BA_TRUST',L_FALSE,
	1     'Max |STEQ| correction kept per depth (fractional units)')
	  IF(BA_TRUST .LT. 1.0_LDP)BA_TRUST=1.0_LDP
	  BA_NREJ_LIM=20
	  CALL RD_STORE_INT(BA_NREJ_LIM,'BA_NREJ_LIM',L_FALSE,
	1     'Stop after this many consecutive rejected BA steps')
	  BA_RES_SAT=1.0E+30_LDP
	  CALL RD_STORE_DBLE(BA_RES_SAT,'BA_RES_SAT',L_FALSE,
	1     'Res level above which the max channel is judged'//
	2     ' saturated and only the RES2 channel gates acceptance')
	  BA_REJ_FRAC=0.25_LDP
	  CALL RD_STORE_DBLE(BA_REJ_FRAC,'BA_REJ_FRAC',L_FALSE,
	1     'Depth-column failure fraction above which a BA step'//
	2     ' rolls back whole; <=0 disables masked acceptance')
	  IF(BA_REJ_FRAC .GT. 1.0_LDP)BA_REJ_FRAC=1.0_LDP
!
! Per-equation quarantine (default FALSE): parks the specific
! (equation,depth) pairs that repeatedly fail acceptance so they
! cannot veto the step or trip the whole-column trust mask.  See
! CONTROL_VARIABLE_MOD for the contract.  Inert without
! DO_BA_RES_ACC.
!
	  DO_BA_QUAR=.FALSE.
	  CALL RD_STORE_LOG(DO_BA_QUAR,'DO_BA_QUAR',L_FALSE,
	1     'Park persistently failing (eqn,depth) pairs out of BA'//
	2     ' acceptance and the trust mask')
	  BA_QUAR_MAX=25
	  CALL RD_STORE_INT(BA_QUAR_MAX,'BA_QUAR_MAX',L_FALSE,
	1     'Max quarantined population equations; <=0 disables'//
	2     ' quarantine growth')
	  BA_QUAR_COUNT=0
	  BA_SE_RES=0.0_LDP
	  BA_SE_RES2=0.0_LDP
	  BAACC_PENDING=.FALSE.
	  BA_STEP_REJECTED=.FALSE.
	  BA_SNAPSHOT_MISSING=.FALSE.
	  BAACC_RES_PREV=0.0_LDP
	  BAACC_RES_MIN=HUGE(1.0_LDP)
	  BAACC_RES2_PREV=0.0_LDP
	  BAACC_RES2_MIN=HUGE(1.0_LDP)
	  BA_TRUST_CUR=BA_TRUST
	  BA_REJ_MAXCH=0.0_LDP
	  BA_NREJ=0
	  BA_NREJ_CONSEC=0
	  BA_PART_THIS=.FALSE.
	  BA_SOFT_THIS=.FALSE.
	  BA_RAW_TRIAL_PENDING=.FALSE.
	  BA_RAW_REJECTED=.FALSE.
	  BA_RDLAM_PROBE=.FALSE.
	  BA_RDLAM_ACCEPT=.FALSE.
	  BAACC_RDLAM_PROBE_ARM=.FALSE.
	  BAACC_RDLAM_PROBE=.FALSE.
	  BA_NREJ_SOFT=0
	  BA_SAT_LIM=0
	  CALL RD_STORE_INT(BA_SAT_LIM,'BA_SAT_LIM',L_FALSE,
	1     'Demote to LAMBDA after this many saturated full steps')
	  BA_SAT_RAW=1.0E+06_LDP
	  CALL RD_STORE_DBLE(BA_SAT_RAW,'BA_SAT_RAW',L_FALSE,
	1     'Raw MAXCH level counted as a saturated BA demand')
	  BA_SAT_CNT=0
!
	  BA_DEM_PFLOOR=0.0_LDP
	  CALL RD_STORE_DBLE(BA_DEM_PFLOOR,'BA_DEM_PFLOOR',L_FALSE,
	1     'Population floor for the BA demanded-correction metric')
	  BA_DEM_RELFRAC=0.0_LDP
	  CALL RD_STORE_DBLE(BA_DEM_RELFRAC,'BA_DEM_RELFRAC',L_FALSE,
	1     'Demand floor as fraction of max population at the depth')
	  DO_FG13_ALO=.FALSE.
	  CALL RD_STORE_LOG(DO_FG13_ALO,'DO_FG13_ALO',L_FALSE,
	1     'Use FG_J_CMF_V13 diagonal ALO in the SOB continuum Jacobian')
!
! Legacy keys retained so existing VADAT files remain readable.  Shock
! deposition makes inner/outer luminosity agreement unsuitable for
! convergence control; both values are diagnostic-only.
!
	  BA_LUM_GATE=.FALSE.
	  CALL RD_STORE_LOG(BA_LUM_GATE,'BA_LUM_GATE',L_FALSE,
	1     'Legacy no-op; luminosity ratio is diagnostic only')
	  BA_LUM_DCHG=0.1_LDP
	  CALL RD_STORE_DBLE(BA_LUM_DCHG,'BA_LUM_DCHG',L_FALSE,
	1     'Legacy no-op luminosity-drift threshold')
!
! CONV_DIAG_LVL (default 2): verbosity of the dev convergence
! diagnostics.  2 reproduces the pre-existing output volume
! (MODCTL, MODE-OUT, BA-CTX, BA-SOLVE, BA13-TRACE plus the
! freeze dumps and the BA_RES_ACC trajectory); 1 drops the
! per-call trace lines and keeps ITER_SUM, ITER_MODE,
! CONV_PLATEAU, the freeze dumps and the trajectory line;
! 0 emits ITER_SUM and ITER_HISTORY only.
!
	  CONV_DIAG_LVL=2
	  CALL RD_STORE_INT(CONV_DIAG_LVL,'CONV_DIAG_LVL',L_FALSE,
	1     'Convergence diagnostic verbosity (0-2)?')
	  IF(CONV_DIAG_LVL .LT. 0)CONV_DIAG_LVL=0
	  IF(CONV_DIAG_LVL .GT. 2)CONV_DIAG_LVL=2
	  LUM_RATIO_OI=1.0_LDP
!
! The Knop chi_hat formalism is only implemented in the INTEGRAL
! (short-characteristics) path of FG_J_CMF_NONMONO_V1.  The DIFFERENCE
! path does not construct chi_hat = chi + |GAM|/dLOG_NU and can be
! unstable for SIGMA < -1.
!
	  IF(DO_NONMONO_KNOP .AND. FG_SOL_OPTIONS(1:4) .EQ. 'DIFF')THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'DO_NONMONO_KNOP requires the INTEGRAL method'
	    WRITE(LUER,*)'(FG_OPT=INT or INT/INS), but FG_OPT=',FG_SOL_OPTIONS
	    WRITE(LUER,*)'The DIFFERENCE method does not implement the'
	    WRITE(LUER,*)'Knop et al. 2009 chi_hat formalism and is unstable'
	    WRITE(LUER,*)'for non-monotonic velocity fields.'
	    STOP
	  END IF
!
	  CALL RD_STORE_LOG(SETZERO,'He2_RES=0',L_TRUE,
	1           'Set rates in He2 resonance lines to zero ?')
!
	  CALL RD_STORE_LOG(OVERLAP,'ALLOW_OL',L_TRUE,
	1           'Allow for overlap of close lines (SOB only) ?')
	  CALL RD_STORE_DBLE(OVER_FREQ_DIF,'OL_DIF',L_TRUE,
	1           'Max. difference (in km/s) for overlap')
	  OVER_FREQ_DIF=OVER_FREQ_DIF/2.998E+05_LDP
!
	  CALL RD_STORE_LOG(INCL_CHG_EXCH,'INC_CHG',L_TRUE,'Include charge exchange reactions?')
	  CALL RD_STORE_LOG(INCL_TWO_PHOT,'INC_TWO',L_TRUE,'Include two photon transitions?')
	  TWO_PHOTON_METHOD='USE_RAD'
	  CALL RD_STORE_CHAR(TWO_PHOTON_METHOD,'TWO_METH',L_FALSE,'USE_RAD, LTE, NOSTIM or OLD_DEFAULT')
	  INCL_PENNING_ION=.FALSE.
	  CALL RD_STORE_LOG(INCL_PENNING_ION,'INC_PEN',SN_MODEL,'Include Penning ionization?')
!
	  CALL RD_STORE_LOG(INCL_RAY_SCAT,'INC_RAY',L_TRUE,
	1           'Include opacity due to Rayleigh scattering?')
	  CALL RD_STORE_LOG(INCL_ADVECTION,'INC_ADV',L_TRUE,
	1           'Include advection terms in rate equations?')
	  CALL RD_STORE_LOG(INCL_ADIABATIC,'INC_AD',L_TRUE,
	1           'Include adiabatic cooling in energy equation')
	  SCL_SL_LINE_OPAC=.FALSE.
	  CALL RD_STORE_LOG(SCL_SL_LINE_OPAC,'SCL_SL_OPAC',L_FALSE,
	1            'Scale SL line opacities/emissivties for heating/cooling consistency?')
	  CALL RD_STORE_LOG(SCL_LINE_COOL_RATES,'SCL_LN',L_TRUE,
	1            'Scale line cooling rate for Rad. Eq. equation?')
	  IF(SCL_SL_LINE_OPAC .AND. SCL_LINE_COOL_RATES)THEN
	    WRITE(LUER,*)'Inconsitent use of SCL_SL_OPAC and SCL_LN'
	    WRITE(LUER,*)'Only one of the two options can be true'
	    STOP
	  END IF
	  CALL RD_STORE_DBLE(SCL_LINE_HT_FAC,'SCL_LN_FAC',L_TRUE,
	1            'Scale line cooling rate for Rad. Eq. equation?')
	  SCL_LINE_DENSITY_LIMIT=1.0E+30_LDP
	  CALL RD_STORE_DBLE(SCL_LINE_DENSITY_LIMIT,'SCL_DEN_LIM',L_FALSE,
	1            'Density beyond which line cooling scaling is switched off')
	  INCLUDE_dSLdT=.FALSE.
	  CALL RD_STORE_LOG(INCLUDE_dSLdT,'INCL_dSLdT',L_FALSE,
	1            'Include variation in distribution of level populations in a SL with T?')
	  NEW_LINE_BA=.FALSE.
	  IF(SN_MODEL)NEW_LINE_BA=.TRUE.
	  INDX_BA_METH_RD=45
	  CALL RD_STORE_LOG(NEW_LINE_BA,'NEW_LINE_BA',L_FALSE,
	1            'Use the new technique for comuting BA_T (i.e., dRE/dT)?')
	  CALL RD_STORE_INT(INDX_BA_METH_RD,'INDX_BA_METH',L_FALSE,'Depth location to begin NEW_LINE_BA technique')
!
	  LINEAR_ADV=.TRUE.
	  CALL RD_STORE_LOG(LINEAR_ADV,'LIN_ADV',L_FALSE,
	1           'Compute advection terms using derivatives in linear plane?')
	  ADVEC_RELAX_PARAM=1.0_LDP
	  CALL RD_STORE_DBLE(ADVEC_RELAX_PARAM,'ADV_RELAX',L_FALSE,
	1           'Parameter to allow advection terms to be included slowly')
!
! Except for the X-ray switch, the X-ray options are only needed if we
! are including X-rays.
!
	  VSMOOTH_XRAYS=3000.0_LDP
	  CALL RD_STORE_LOG(XRAYS,'INC_XRAYS',L_TRUE,
	1           'Include X-ray emission')
	  CALL RD_STORE_LOG(FF_XRAYS,'FF_XRAYS',XRAYS,
	1           'Use free-free processes to compute X-ray emission')
	  CALL RD_STORE_LOG(XRAY_SMOOTH_WIND,'X_SM_WIND',XRAYS,
	1           'Ignore clumping when computing X-ray emission')
	  CALL RD_STORE_DBLE(VSMOOTH_XRAYS,'VS_XRAYS',XRAYS,
	1           'X-ray smoothing width for SOB/CMF options')
	  CALL RD_STORE_DBLE(FILL_FAC_XRAYS_1,'FIL_FAC_1',XRAYS,
	1           'Filling factor for X-ray emission [1]')
	  CALL RD_STORE_DBLE(T_SHOCK_1,'T_SHOCK_1',XRAYS,
	1           'Shock T for X-ray emission [1]')
	  CALL RD_STORE_DBLE(V_SHOCK_1,'V_SHOCK_1',XRAYS,
	1           'Cut off velocity for X-ray emission [1]')
	  CALL RD_STORE_DBLE(FILL_FAC_XRAYS_2,'FIL_FAC_2',XRAYS,
	1           'Filling factor for X-ray emission [2]')
	  CALL RD_STORE_DBLE(T_SHOCK_2,'T_SHOCK_2',XRAYS,
	1           'Shock T for X-ray emission [2]')
	  CALL RD_STORE_DBLE(V_SHOCK_2,'V_SHOCK_2',XRAYS,
	1           'Cut off velocity for X-ray emission [2]')
	  CALL RD_STORE_LOG(ADD_XRAYS_SLOWLY,'XSLOW',XRAYS,
	1           'Add X-rays by slowly increasing filling factors?')
	  CALL RD_STORE_DBLE(FILL_FAC_X1_BEG,'XFI1_BEG',ADD_XRAYS_SLOWLY,
	1           'Initial filling factor for X-ray emission [1]')
	  CALL RD_STORE_DBLE(FILL_FAC_X2_BEG,'XFI2_BEG',ADD_XRAYS_SLOWLY,
	1           'Initial filling factor for X-ray emission [2]')
	  CALL RD_STORE_DBLE(SLOW_XRAY_SCL_FAC,'XSCL_FAC',ADD_XRAYS_SLOWLY,
	1           'Rate to increase X-ray filling factor')
!
	  SCALE_XRAY_LUM=.FALSE.
	  ALLOWED_XRAY_FLUX_ERROR=0.1_LDP
	  CALL RD_STORE_LOG(SCALE_XRAY_LUM,'SCL_XLUM',L_FALSE,
	1           'Scale X-ray emissivities to get a specified X-ray luminosiity')
	  CALL RD_STORE_DBLE(DESIRED_XRAY_LUM,'XRAY_LUM',SCALE_XRAY_LUM,
	1           'X-ray luminosity in units of LSTAR')
	  CALL RD_STORE_DBLE(ALLOWED_XRAY_FLUX_ERROR,'XRAY_ERR',L_FALSE,
	1           'Fractional error allowed in observed X-ray luminosity')
!
	  DELV_XRAY=0.5_LDP*VSMOOTH_XRAYS
	  CALL RD_STORE_DBLE(DELV_XRAY,'V_XRAY',L_FALSE,
	1            'Max. vel. sep. (km/s) between evaluations of '//
	1            '  phot. cross-sections in X-ray region?')
	  NU_XRAY_END=100.0_LDP
	  CALL RD_STORE_DBLE(NU_XRAY_END,'NU_XRAY',L_FALSE,
	1            'End of X-ray region for continuum definition')
!	
	  AUTO_ADD_ION_STAGES=.TRUE.
	  CALL RD_STORE_LOG(AUTO_ADD_ION_STAGES,'AUTO_ADD',L_FALSE,
	1       'Automatically add files with missing XzV_IN files?')
	  AUTO_SMOOTH_POPS=.FALSE.
	  CALL RD_STORE_LOG(AUTO_SMOOTH_POPS,'AUTO_SMOOTH',L_FALSE,
	1       'Automatically smooth populations as we iterate?')
!
	  WRITE(LUSCR,'()')
	  WRITE(LUSCR,'()')
	  INTERP_DC_SPH_TAU=.FALSE.
	  SET_LTE_AS_INIT_ESTIMATES=.FALSE.
	  CALL RD_STORE_LOG( RDINR,'RD_IN_R_GRID',L_TRUE,
	1        'Read in a predetermined R grid ?')
	  CALL RD_STORE_LOG(GRID,'LIN_INT',L_TRUE,
	1        'Use direct linear interpolation if  new model ?')
	  DC_INTERP_METHOD='ED'
	  IF(GRID)DC_INTERP_METHOD='R'
	  CALL RD_STORE_NCHAR(DC_INTERP_METHOD,'DC_METH',ITEN,L_FALSE,
	1        'Interpolation method: TR, R, ED, LTE, SPH_TAU')
	  CALL RD_STORE_LOG(INTERP_DC_SPH_TAU,'DC_SPH_TAU',L_FALSE,
	1        'Interpolate d.c''s on the spherical TAU scale?')
	  IF(INTERP_DC_SPH_TAU)DC_INTERP_METHOD='SPH_TAU'
	  CALL RD_STORE_LOG(SET_LTE_AS_INIT_ESTIMATES,'LTE_EST',L_FALSE,
	1        'Use LTE for the initial estimates')
!
	  IF(SET_LTE_AS_INIT_ESTIMATES)DC_INTERP_METHOD='LTE'
	  IF(DC_INTERP_METHOD .EQ. 'LTE')THEN
	    CALL RD_STORE_DBLE(T_EXCITE_MIN,'T_EXC',L_TRUE,'Minimum temp. when using LTE option')
	  END IF
	  CALL RD_STORE_LOG(DO_POP_SCALE,'POP_SCALE',L_TRUE,
	1        'Scale populations so that cons. Eq. satisfied ?')
	  CALL RD_STORE_DBLE(T_INIT_TAU,'T_INIT_TAU',L_TRUE,
	1        'Tau above which T is set exactly to T(spherical)')
	  CALL RD_STORE_LOG(ITERATE_INIT_T,'IT_ON_T',L_TRUE,
	1        'Improve initial T estimate by iteration ?')
	  INTERP_T_ON_R_GRID=.TRUE.
	  IF(DC_INTERP_METHOD .EQ. 'R')INTERP_T_ON_R_GRID=.FALSE.
	  CALL RD_STORE_LOG(INTERP_T_ON_R_GRID,'T_ON_R',L_FALSE,
	1       'Interpolate T on R grid - default is SPH tau grid')
	  CALL RD_STORE_DBLE(GREY_PAR,'GREY_TAU',L_TRUE,
	1        'Specifys Tau above which T is set to TGREY in iterative process')
	  MAX_NO_GREY_ITERATIONS=5
	  CALL RD_STORE_INT(MAX_NO_GREY_ITERATIONS,'NGREY_ITS',L_FALSE,
	1        'Number of TGREY iterations')
!
! We only need to specify the TRANSITION type if GLOBAL_LINE_SWITCH is NONE.
!
	  IF(GLOBAL_LINE_SWITCH(1:4) .EQ. 'NONE')THEN
	    WRITE(LUSCR,'()')
	    DO ID=1,NUM_IONS-1
	      TMP_KEY='TRANS_'//TRIM(ION_ID(ID))
	      TMP_STRING='Method for treating '//TRIM(ION_ID(ID))//' lines?'
	      CALL RD_STORE_NCHAR(ATM(ID)%XzV_TRANS_TYPE,TMP_KEY,
	1          ISIX,ATM(ID)%XZV_PRES,TMP_STRING)
	    END DO
	  END IF
!
	  WRITE(LUSCR,'()')
	  DO ID=1,NUM_IONS-1
	    TMP_KEY='DIE_'//TRIM(ION_ID(ID))
	    TMP_STRING='Include (?) LTDR AUOT, WI calc.s for '//TRIM(ION_ID(ID))
	    CALL RD_STORE_2LOG(ATM(ID)%DIE_AUTO_XzV,ATM(ID)%DIE_WI_XzV,
	1         TMP_KEY,ATM(ID)%XZV_PRES,TMP_STRING)
	  END DO
!
! When adding new species / ioization stages it is sometime useful to hold
! some levels fixed. This can be done in two ways. If we set FIX_ALL_SPEC true,
! all levels will be held fixed UNLESS we do an UNFIX_XzV command. This will be the quickest,
! and safest approach when adding a few ioization stages. We still need to set FIX_T and FIX_NE.
!
! Alternatively, we can set FIX_XzV for every ionization stage.
!
! NB: None of  the keywords need be present.
!
	  FIX_ALL_SPECIES=.FALSE.
	  CALL RD_STORE_LOG(FIX_ALL_SPECIES,'FIX_ALL_SPEC',L_FALSE,'Fix all species?')
!
! Since we don't use FIX_XzV very much, these have been changed to a hidden variable.
!
	  IF(FIX_ALL_SPECIES)THEN
	    DO ISPEC=1,NUM_SPECIES
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        ATM(ID)%FIX_NXzV=ATM(ID)%NXzV
	      END DO
	      FIX_SPECIES(ISPEC)=1
	    END DO
!
! We now allow the possibility of unfixing some levels.
! The use of UFIX_XzV avoids confusion and abiguity with the FIX_XzV command.
!
	    DO ISPEC=1,NUM_SPECIES
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        TMP_KEY='UNFIX_'//TRIM(ION_ID(ID))
	        TMP_STRING='Fix ? levels of '//TRIM(ION_ID(ID))//' : set to 0 to unfix'
	        I=-1; CALL RD_STORE_INT(I,TMP_KEY,L_FALSE,TMP_STRING)
	        IF(I .NE. -1)ATM(ID)%FIX_NXzV=I
	      END DO
	      TMP_KEY='FIX_'//TRIM(SPECIES(ISPEC))
	      TMP_STRING='Fix (?) highest ionization stage in '//TRIM(SPECIES(ISPEC))
	      I=-1; CALL RD_STORE_INT(I,TMP_KEY,L_FALSE,TMP_STRING)
	      IF(I .NE. -1)FIX_SPECIES(ISPEC)=I
	    END DO
!
	  ELSE
!
! Using the fix command to fix levels.
!
	    DO ISPEC=1,NUM_SPECIES
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ID .EQ. SPECIES_BEG_ID(ISPEC))WRITE(LUSCR,'()')
	        TMP_KEY='FIX_'//TRIM(ION_ID(ID))
	        TMP_STRING='Fix ? levels of '//TRIM(ION_ID(ID))
	        ATM(ID)%FIX_NXzV=0
	        CALL RD_STORE_INT(ATM(ID)%FIX_NXzV,TMP_KEY,L_FALSE,TMP_STRING)
	      END DO
	      TMP_KEY='FIX_'//TRIM(SPECIES(ISPEC))
	      TMP_STRING='Fix (?) highest ionization stage in '//TRIM(SPECIES(ISPEC))
	      FIX_SPECIES(ISPEC)=0
	      CALL RD_STORE_INT(FIX_SPECIES(ISPEC),TMP_KEY,L_FALSE,TMP_STRING)
	    END DO
	  END IF
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(FIXED_NE,'FIX_NE',L_TRUE,'Fix the electron density ?')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(RD_FIX_IMP,'FIX_IMP',L_TRUE,'Automatically fix impurity species?')
!
	  WRITE(LUSCR,'()')
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(RD_FIX_T,'FIX_T',L_TRUE,
	1            'Keep the Temperature fixed ?')
          FIX_IN_BOUND_T=.FALSE.
          CALL RD_STORE_LOG(FIX_IN_BOUND_T,'FIX_INB_T',L_FALSE,
	1            'Fix the Temperature at the inner boundary ?')
          FIX_LST_X_DPTHS=0
	  CALL RD_STORE_INT(FIX_LST_X_DPTHS,'FIX_X_DPTH',L_FALSE,
	1            'Fix the Temperature at the last ? depths')
	  CALL RD_STORE_LOG(VARFIXT,'FIX_T_AUTO',L_TRUE,
	1            'Fix the Temperature automatically ?')
!
! JFNK_TAUTO defaults to VARFIXT: FIX_T_AUTO=T switches the JFNK
! solver from hard-masking the variable-T range to adaptively
! damping it.  An explicit JFNK_TAUTO key overrides.
!
	  JFNK_TAUTO=VARFIXT
	  CALL RD_STORE_LOG(JFNK_TAUTO,'JFNK_TAUTO',JFNK_TAUTO,
	1     'Adaptive damped-T domain inside JFNK (default=FIX_T_AUTO)')
	  CALL RD_STORE_DBLE(TAU_SCL_T,'TAU_SCL_T',L_TRUE,
	1      'Electron scattering optical depth from which to fix T')
	  IF(TAU_SCL_T .EQ. 0.0_LDP)THEN
	    CON_SCL_T=0.0_LDP
	  ELSE
	    CON_SCL_T=1000.0_LDP
	  END IF
	  T_MIN=0.0_LDP
	  CALL RD_STORE_DBLE(T_MIN,'T_MIN',L_FALSE,'Minimum electron temperature')
	  DO_SRCE_VAR_ONLY=.FALSE.
	  CALL RD_STORE_LOG(DO_SRCE_VAR_ONLY,'SRCE_ONLY',L_FALSE,
	1            'Allow only ths source function to vary?')
	  ADD_ADDITIONAL_OPACITY=.FALSE.
	  ADD_OPAC_SCL_FAC=0.0_LDP
	  CALL RD_STORE_LOG(ADD_ADDITIONAL_OPACITY,'ADD_OPAC',L_FALSE,
	1            'Add additional pacity to help converge model?')
	  CALL RD_STORE_DBLE(ADD_OPAC_SCL_FAC,'OP_SCL_FAC',ADD_ADDITIONAL_OPACITY,
	1            'Scale factor for addition of extra opacity')
!
	  USE_ELEC_HEAT_BAL=.FALSE.
	  CALL RD_STORE_LOG(USE_ELEC_HEAT_BAL,'USE_EHB',L_FALSE,
	1            'Use electron heating/cooling balance as T constraint?')
	  COMP_STEQ_T_EHB=USE_ELEC_HEAT_BAL
	  CALL RD_STORE_LOG(COMP_STEQ_T_EHB,'COMP_EHB',L_FALSE,
	1            'Compute electron heating/cooling balance as T constraint check?')

!
	  CALL RD_STORE_NCHAR(METH_SOL,'SOL_METH',ISIX,L_TRUE,
	1            'Which Method To solve Matrix Equations'//
	1            ' DIAG, TRI, PEN, GSIT or MIN')
	  IF(NUM_BNDS .EQ. 1 .AND. METH_SOL .NE. 'DIAG')THEN
	    WRITE(LUER,*)'****************************************'
	    WRITE(LUER,*)'******WARNING in CMFGEN*****************'
	    WRITE(LUER,*)'Solution method inconsistent with NUM_BNDS'
	    WRITE(LUER,*)'METH_SOL=',METH_SOL,'NUM_BNDS=',NUM_BNDS
	    WRITE(LUER,'(X,A,/)')'Using diagnoal solution'
	    METH_SOL='DIAG'
	  END IF
	  CALL RD_STORE_NCHAR(SCALE_OPT,'SCALE_OPT',ISIX,L_TRUE,
	1           'Scale option (LOCAL, NONE or GLOBAL) ? ')
	  LAM_SCALE_OPT='LIMIT'
	  CALL RD_STORE_NCHAR(LAM_SCALE_OPT,'LAM_SCALE_OPT',ISIX,L_FALSE,
	1           'Only has an effect if set to LIMIT')
	  CALL RD_STORE_DBLE(EPS,'EPS_TERM',L_TRUE,
	1      'If maximum fractional % change < EPS terminate model ')
	  CALL RD_STORE_DBLE(MAX_LIN_COR,'MAX_LIN',L_TRUE,
	1      'Maximum fractional change for linearization ')
	  CALL RD_STORE_DBLE(MAX_LAM_COR,'MAX_LAM',L_TRUE,
	1      'Maximum fractional change for lambda iteration ')
	  MAX_dT_COR=0.2_LDP
	  CALL RD_STORE_DBLE(MAX_dT_COR,'MAX_dT',L_FALSE,
	1      'Maximum fractional change in the temperature')
	  IF(MAX_dT_COR .LE. 0.0_LDP .OR. MAX_dT_COR .GT. 0.201_LDP)THEN
	    WRITE(LUER,*)' Error: MAX_dT in VADAT has an invalid value of',MAX_dT_COR
	    WRITE(LUER,*)' Require 0 < MAX_dT_COR < 0.2'
	    STOP
	  END IF
	  CALL RD_STORE_DBLE(MAX_CHNG_LIM,'MAX_CHNG',L_TRUE,
	1      'If maximum % fractional change > MAX_CHNG terminate model ')
!
	  CALL RD_STORE_LOG(COMPUTE_BARDIN,'COMP_BA',L_TRUE,
	1            'Compute BA matrix ?')
	  CALL RD_STORE_LOG(WRBAMAT_RDIN,'STORE_BA',L_TRUE,
	1      'Store the BA matrix for/during each iteration ? ')
	  CALL RD_STORE_LOG(WR_BA_INV,'STORE_BA_INV',L_TRUE,
	1      'Store the INVERSE of the BA matrix on each iteration ? ')
	  CALL RD_STORE_LOG(WR_PART_OF_INV,'WR_PRT_INV',L_TRUE,
	1      'Store part of the INVERSE to reduce storage TRIDIAG only)?')
	  CALL RD_STORE_INT(N_ITS_TO_FIX_BA,'N_FIX_BA',L_TRUE,
	1      'Number of iterations to hold BA fixed')
	  CALL RD_STORE_DBLE(BA_CHK_FAC,'BA_CHK_FAC',L_TRUE,
	1      'If dJ < BA_CHK_FAC*RJ, ignore correction to BA')
	  CALL RD_STORE_DBLE(VAL_FIX_BA,'FIX_BA',L_TRUE,
	1      'Switch off BA computation if MAXCH< VAL_FIX_BA ')
	  CALL RD_STORE_DBLE(VAL_DO_LAM,'LAM_VAL',L_TRUE,
	1      'Do Lambda iterations if MAXCH > VAL_DO_LAM')
	  CALL RD_STORE_INT(RD_CNT_LAM,'NUM_LAM',L_TRUE,
	1      '# of Lambda iterations if MAXCH > VAL_DO_LAM')
	  CNT_LAM=0
!
! LAMBDA<->full hysteresis: LAM_ENTER is the MAXCH_CTRL level that
! forces (back) into LAMBDA, LAM_EXIT is the level a LAMBDA sequence
! must stay under for LAM_N_STABLE consecutive iterations before a
! full step is allowed again.  Defaulting both to LAM_VAL and the
! stable count to 1 reproduces the legacy single-threshold switch.
!
	  VAL_ENTER_LAM=VAL_DO_LAM
	  CALL RD_STORE_DBLE(VAL_ENTER_LAM,'LAM_ENTER',L_TRUE,
	1      'Re-enter LAMBDA if MAXCH_CTRL > LAM_ENTER')
	  VAL_EXIT_LAM=VAL_DO_LAM
	  CALL RD_STORE_DBLE(VAL_EXIT_LAM,'LAM_EXIT',L_TRUE,
	1      'Exit LAMBDA when MAXCH_CTRL < LAM_EXIT')
	  N_STABLE_LAM=1
	  CALL RD_STORE_INT(N_STABLE_LAM,'LAM_N_STABLE',L_TRUE,
	1      '# consecutive below-LAM_EXIT its needed to exit LAMBDA')
	  IF(N_STABLE_LAM .LT. 1)N_STABLE_LAM=1
	  CNT_STABLE_LAM=0
!
! BA_ENTER_RAW (optional, default 1E30 = gate off): while the raw
! demand of the last non-LAMBDA step is still at/above this level,
! a LAMBDA quiet window is not allowed to release a full step.
! Optional so VADAT files without the key are unaffected.
!
	  BA_ENTER_RAW=1.0E30_LDP
	  CALL RD_STORE_DBLE(BA_ENTER_RAW,'BA_ENTER_RAW',L_FALSE,
	1      'Hold LAMBDA->full exit while last raw BA demand >= this')
	  BA_LAST_RAW=0.0_LDP
!
! BA_PROBE_LAM (optional, default 0 = off): while the raw-demand
! gate above is holding, release one probe step after this many
! LAMBDA calls so BA_LAST_RAW can be re-measured.  Optional so
! VADAT files without the key get NUM_LAM-only probing.
!
	  BA_PROBE_LAM=0
	  CALL RD_STORE_INT(BA_PROBE_LAM,'BA_PROBE_LAM',L_FALSE,
	1      'Probe full step after N LAMBDA calls while gate holds')
!
! BA_DAMP_CAP (optional, default 0 = off): when positive, a
! saturated-demand full step stays in full mode and the next full
! call's correction factor is capped at MIN(MAX_LIN,BA_DAMP_CAP).
! Arming it gives BA_SAT_LIM a default leash of 8 so a saturated run
! is always bounded; set BA_SAT_LIM explicitly to change it.
!
	  BA_DAMP_CAP=0.0_LDP
	  CALL RD_STORE_DBLE(BA_DAMP_CAP,'BA_DAMP_CAP',L_FALSE,
	1      'Factor cap for damped full steps on saturated demand')
	  IF(BA_DAMP_CAP .GT. 0.0_LDP .AND. BA_SAT_LIM .LE. 0)
	1    BA_SAT_LIM=8
!
	  FRESH_BA_EXIT=.FALSE.
	  CALL RD_STORE_LOG(FRESH_BA_EXIT,'FRESH_BA_EXIT',L_TRUE,
	1      'Force fresh BA matrix on LAMBDA->full transition?')
!
	  N_FULL_TREL=0
	  CALL RD_STORE_INT(N_FULL_TREL,'N_FULL_TREL',L_TRUE,
	1      'Hold T for this many accepted full its after LAMBDA exit?')
	  IF(N_FULL_TREL .LT. 0)N_FULL_TREL=0
	  T_REL_CHG=3.0_LDP
	  CALL RD_STORE_DBLE(T_REL_CHG,'T_REL_CHG',L_TRUE,
	1      'MAXCH_CTRL level that also releases the post-LAMBDA T hold')
	  N_ACC_FULL=0
	  T_HELD_PHASE=.FALSE.
	  CALL RD_STORE_LOG(RDINSOL,'RD_SOL',L_TRUE,
	1            'RD in solution vector to update populations')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(EDD_CONT,'JC_W_EDD',L_TRUE,
	1        'Compute continuum intensity using Eddington factors')
	  NO_VEL_FOR_CONTINUUM=.FALSE.
	  CALL RD_STORE_LOG(NO_VEL_FOR_CONTINUUM,'NOV_CONT',L_FALSE,
	1        'In non-blanketed mode, ignore velocity terms when computing continuum')
	  CALL RD_STORE_LOG(EDD_LINECONT,'JBAR_W_EDD',L_TRUE,
	1    'Compute line continuum intensity using Eddington factors')
!
	  PLANE_PARALLEL_NO_V=.FALSE.
	  PLANE_PARALLEL=.FALSE.
	  CALL RD_STORE_LOG(PLANE_PARALLEL_NO_V,'PP_NOV',L_FALSE,'Plane-paralle geometry WITHOUT velocity field?')
	  CALL RD_STORE_LOG(PLANE_PARALLEL,'PP_MOD',L_FALSE,'Plane-paralle geometry with velocity field?')
	  IF(PLANE_PARALLEL_NO_V .AND. PLANE_PARALLEL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'PP_NOV and PP_MOD cannot both be true at the sdame time'
	    STOP
	  END IF
!
	  INCL_DJDT_TERMS=.FALSE.
	  CALL RD_STORE_LOG(INCL_DJDT_TERMS,'INCL_DJDT',L_FALSE,'DJDt terms in transfer equaton for SN models?')
	  IF(INCL_DJDT_TERMS)THEN
	    USE_DJDT_RTE=.TRUE.
	    USE_Dr4JDT=.TRUE.
	    CALL RD_STORE_LOG(USE_Dr4JDT,'USE_DR4JDT',L_FALSE,'Difference Dr4JDt')
	  ELSE
	    JGREY_WITH_V_TERMS=.FALSE.
	    USE_DJDT_RTE=.FALSE.
	    USE_Dr4JDT=.FALSE.
	    CALL RD_STORE_LOG(USE_DJDT_RTE,'USE_DJDT_RTE',L_FALSE,
	1    'Use solver which has DJDt terms in transfer equaton for SN models?')
	  END IF
	  DJDT_RELAX_PARAM=1.0_LDP
	  CALL RD_STORE_DBLE(DJDT_RELAX_PARAM,'DJDT_RELAX',L_FALSE,
	1          'Factor to scale DJDT terms to assist initial convergence')
!
	  USE_LAM_ES=.FALSE.
	  USE_J_REL=.FALSE.
	  INCL_REL_TERMS=.FALSE.
	  INCL_ADVEC_TERMS_IN_TRANS_EQ=.FALSE.
	  USE_FORMAL_REL=.FALSE.
	  CALL RD_STORE_LOG(USE_LAM_ES,'USE_LAM_ES',L_FALSE,'Use lambda iteration with ray solutions?')
	  CALL RD_STORE_LOG(USE_J_REL,'USE_J_REL',L_FALSE,'Use MOM_J_REL_VN to solve the moment equations?')
	  IF(USE_DJDT_RTE .OR. USE_J_REL)USE_FORMAL_REL=.TRUE.
	  CALL RD_STORE_LOG(USE_FORMAL_REL,'USE_FRM_REL',L_FALSE,'Use CMF_FORMAL_REL to compute F etc?')
	  CALL RD_STORE_LOG(INCL_REL_TERMS,'INCL_REL',USE_J_REL,'Include relativistic terms in the transfer equaton?')
	  IF(INCL_REL_TERMS .AND. .NOT. USE_J_REL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'Can only include relativistic terms if USE_J_REL is set to TRUE'
	    STOP
	  END IF
	  IF(USE_J_REL .AND. INCL_REL_TERMS)INCL_ADVEC_TERMS_IN_TRANS_EQ=.TRUE.
	  CALL RD_STORE_LOG(INCL_ADVEC_TERMS_IN_TRANS_EQ,'INCL_ADV_TRANS',USE_J_REL,
	1    'Include advection terms in the transfer equaton?')
	  IF(INCL_ADVEC_TERMS_IN_TRANS_EQ .AND. .NOT. USE_J_REL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'Can only include advection terms in transfer equation if USE_J_REL is set to TRUE'
	    STOP
	  END IF
	  IF(USE_DJDT_RTE .AND. USE_J_REL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'USE_DJDT_RTE and USE_J_REL cannot be TRUE at the same time'
	    STOP
	  END IF
	  IF( (PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V) .AND. (USE_DJDT_RTE .OR. USE_J_REL) )THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'USE_DJDT_RET and USE_J_REL cannot be TRUE at the same time'
	    STOP
	  END IF
!
! DO_NONMONO_SOB only steers the Sobolev line routines (SUB_SOB_LINE_V3,
! SOBEW*) through SIGMA_TRUE; if no line is SOB-treated the flag is dead
! code.  GLOBAL_LINE=SOB is therefore mandatory [GLOBAL_LINE].  The
! GLOBAL_LINE=NONE per-species form is accepted only when every ion
! present has TRANS_*=SOB (equivalent setup); a CMF or BLANK species
! would silently keep the stock clamped-sigma treatment the flag exists
! to replace.
!
	  IF(DO_NONMONO_SOB)THEN
	    WRITE(LUSCR,'(A)')' DO_NONMONO_SOB=T: using unclamped dlnV/dlnR-1'//
	1     ' (SIGMA_TRUE) with Dessart+2015 tau floor in Sobolev routines'
	    IF(GLOBAL_LINE_SWITCH(1:4) .EQ. 'NONE')THEN
	      DO ID=1,NUM_IONS-1
	        IF(ATM(ID)%XZV_PRES .AND.
	1          ATM(ID)%XzV_TRANS_TYPE(1:3) .NE. 'SOB')THEN
	          WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	          WRITE(LUER,*)'DO_NONMONO_SOB=T with GLOBAL_LINE=NONE'
	          WRITE(LUER,*)'requires TRANS_*=SOB for every ion present;',
	1            ' TRANS_'//TRIM(ION_ID(ID)),'=',
	2            ATM(ID)%XzV_TRANS_TYPE,' violates this.'
	          WRITE(LUER,*)'Set GLOBAL_LINE=SOB, or set every TRANS_*',
	1            ' key to SOB.'
	          STOP
	        END IF
	      END DO
	    ELSE IF(GLOBAL_LINE_SWITCH(1:3) .NE. 'SOB')THEN
	      WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	      WRITE(LUER,*)'DO_NONMONO_SOB=T requires GLOBAL_LINE=SOB'
	      WRITE(LUER,*)'(or GLOBAL_LINE=NONE with all TRANS_*=SOB).'
	      WRITE(LUER,*)'GLOBAL_LINE=',GLOBAL_LINE_SWITCH,' leaves no',
	1        ' line on the Sobolev routines the flag steers.'
	      STOP
	    END IF
!
! The verified non-monotonic SOB configuration takes J from the formal
! ray solution (USE_LAM_ES=T, the pipeline default).  USE_LAM_ES=F is
! permitted as an opt-in variant: continuum J then comes from
! MOM_J_CMF_NONMONO_V1 (12-Oct-2026; the stock MOM_J_CMF_V11 assumes
! expanding flow and diverges on this velocity field), with the
! fixed-J + GET_MOM_ALO_DIAG linearization in the BA matrix - the same
! pairing as KNOP+LAM_ES=F.  FG_J_CMF_V13 still runs each frequency,
! so the Eddington factors (N_ON_J etc.) come from its ray moments via
! GET_MOMS_NON_REL, and SIGMA_TRUE continues to drive the Sobolev
! lines unchanged.  Warn rather than STOP so the default
! (pipeline-forced LAM_ES=T) path is untouched.
!
	    IF(.NOT. USE_LAM_ES)THEN
	      WRITE(LUSCR,'(A)')' DO_NONMONO_SOB=T with USE_LAM_ES=F:'//
	1       ' continuum J from MOM_J_CMF_NONMONO_V1 (SIGMA_TRUE);'//
	2       ' fixed-J + MOM diagonal ALO linearization active'
	      WRITE(LUER,*)'Warning: DO_NONMONO_SOB with USE_LAM_ES=F uses'
	1       ,' the MOM_J_CMF_NONMONO_V1 moment solve (sign-split'
	2       ,' upwind on SIGMA_TRUE); Sobolev lines unchanged.'
	    END IF
!
! LAM_SET is redundant with GLOBAL_LINE=SOB (every line is already SOB)
! and resets GLOBAL_LINE_SWITCH to NONE at run time, which changes the
! OVERLAP handling in CMFGEN_SUB.  Warn only.
!
	    IF(SET_TRANS_TYPE_BY_LAM)THEN
	      WRITE(LUER,*)'Warning: LAM_SET=T is redundant with'
	      WRITE(LUER,*)'GLOBAL_LINE=SOB (all lines are already SOB) and'
	      WRITE(LUER,*)'resets GLOBAL_LINE_SWITCH to NONE at run time.',
	1        ' Recommend LAM_SET=F.'
	    END IF
	    IF(DO_NONLOC_SOB)THEN
	      WRITE(LUSCR,'(A)')' DO_NONLOC_SOB=T: Rybicki-Hummer 1978'//
	1         ' non-local coupling added to Sobolev net rates'
	    END IF
	  END IF
!
! DO_NONMONO_KNOP consistency checks (deferred until the solver-mode
! flags have been read).  The Knop formal solver replaces FG_J_CMF_V13
! on the standard (non-extended) grid.  Two J-update modes are valid:
!   USE_LAM_ES=T : J from the formal solution (Lambda iteration).
!   USE_LAM_ES=F : J from MOM_J_CMF_NONMONO_V1, the non-monotonic
!     variable-Eddington-factor moment solve (chi_hat generalized
!     opacity, sign-dependent frequency upwind, J_NEXT Gauss-Seidel) -
!     same equations as the formal solver, faster convergence.
! The stock MOM_J_CMF_V11 (clamped SIGMA, one-way PREV coupling) is
! never selected in KNOP mode, so no USE_LAM_ES guard is needed here.
!
	  IF(DO_NONMONO_KNOP)THEN
	    IF(USE_LAM_ES)THEN
	      WRITE(LUSCR,'(A)')' DO_NONMONO_KNOP=T: Knop et al. 2009'//
	1       ' formal solver (FG_J_CMF_NONMONO_V1, USE_LAM_ES)'
	    ELSE
	      WRITE(LUSCR,'(A)')' DO_NONMONO_KNOP=T: Knop et al. 2009'//
	1       ' formal + non-monotonic moment solver (SIGMA_TRUE)'
	      WRITE(LUER,*)'RD_CONTROL_VARIABLES: USE_LAM_ES=F with'
	      WRITE(LUER,*)'DO_NONMONO_KNOP=T - using MOM_J_CMF_NONMONO_V1'
	      WRITE(LUER,*)'(non-monotonic Knop moment solver).'
	    END IF
	    IF(USE_FIXED_J .OR. USE_FORMAL_REL .OR. USE_J_REL .OR.
	1      USE_DJDT_RTE .OR. PLANE_PARALLEL .OR.
	1      PLANE_PARALLEL_NO_V)THEN
	      WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	      WRITE(LUER,*)'DO_NONMONO_KNOP is incompatible with USE_FIXED_J,'
	      WRITE(LUER,*)'relativistic/DJDT solvers, and plane-parallel'
	      WRITE(LUER,*)'modes.  These paths bypass FG_J_CMF_NONMONO_V1.'
	      STOP
	    END IF
!
! The Knop solver only replaces the blanketed (continuum+line) transfer
! in COMP_J_BLANK; GLOBAL_LINE=BLANK is therefore mandatory.  The
! GLOBAL_LINE=NONE per-species form is accepted only when every ion
! present has TRANS_*=BLANK (equivalent setup); a SOB or CMF species
! would silently keep the stock clamped-sigma treatment the flag exists
! to replace.
!
	    IF(GLOBAL_LINE_SWITCH(1:4) .EQ. 'NONE')THEN
	      DO ID=1,NUM_IONS-1
	        IF(ATM(ID)%XZV_PRES .AND.
	1          ATM(ID)%XzV_TRANS_TYPE(1:5) .NE. 'BLANK')THEN
	          WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	          WRITE(LUER,*)'DO_NONMONO_KNOP=T with GLOBAL_LINE=NONE'
	          WRITE(LUER,*)'requires TRANS_*=BLANK for every ion',
	1            ' present; TRANS_'//TRIM(ION_ID(ID)),'=',
	2            ATM(ID)%XzV_TRANS_TYPE,' violates this.'
	          WRITE(LUER,*)'Set GLOBAL_LINE=BLANK, or set every TRANS_*',
	1            ' key to BLANK.'
	          STOP
	        END IF
	      END DO
	    ELSE IF(GLOBAL_LINE_SWITCH(1:5) .NE. 'BLANK')THEN
	      WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	      WRITE(LUER,*)'DO_NONMONO_KNOP=T requires GLOBAL_LINE=BLANK'
	      WRITE(LUER,*)'(or GLOBAL_LINE=NONE with all TRANS_*=BLANK).'
	      WRITE(LUER,*)'GLOBAL_LINE=',GLOBAL_LINE_SWITCH,' routes line',
	1        ' transfer outside COMP_J_BLANK, where the Knop solver',
	2        ' lives.'
	      STOP
	    END IF
!
! LAM_SET rewrites lines outside [F_LAM_BEG,F_LAM_END] to SOB at run
! time and resets GLOBAL_LINE_SWITCH to NONE: those lines would get the
! stock clamped-sigma Sobolev treatment, inconsistent with the
! non-monotonic field.  Not allowed under [DO_NONMONO_KNOP].
!
	    IF(SET_TRANS_TYPE_BY_LAM)THEN
	      WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	      WRITE(LUER,*)'DO_NONMONO_KNOP=T requires LAM_SET=F.'
	      WRITE(LUER,*)'LAM_SET would give edge lines the clamped-sigma'
	      WRITE(LUER,*)'Sobolev treatment, inconsistent with the Knop'
	      WRITE(LUER,*)'non-monotonic field.'
	      STOP
	    END IF
!
! DO_BA_RES_ACC is not a hard requirement, but it is the verified guard:
! a BA step taken against an immature lagged Knop field needs the
! residual gate to roll back.  Warn rather than abort so controlled
! comparisons (F baselines) still run.
!
	    IF(.NOT. DO_BA_RES_ACC)THEN
	      WRITE(LUER,*)'Warning: DO_BA_RES_ACC=F under DO_NONMONO_KNOP.'
	      WRITE(LUER,*)'BA steps are not residual-accepted; the verified'
	      WRITE(LUER,*)'KNOP configuration uses DO_BA_RES_ACC=T.'
	    END IF
	  END IF
!
! GET_NON_LOCAL_GAMMA_ENERGY_V2 maps NONLOCAL and cached RAD_TRANS
! deposition tables on V; that coordinate is multi-valued under a
! non-monotonic field.  An uncached RAD_TRANS call instead uses
! GAMRAY_SUB_V3/CMF_FORMAL_REL_V4_GAM, a signed-SIGMA legacy solver.
! LOCAL and ABS_TRANS do not use V as an independent coordinate.
!
	  IF(SN_MODEL .AND. (DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND.
	1     (GAMRAY_TRANS .EQ. 'RAD_TRANS' .OR.
	2      GAMRAY_TRANS .EQ. 'NONLOCAL'))THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'GAMRAY_TRANS=',GAMRAY_TRANS,' is incompatible'
	    WRITE(LUER,*)'with non-monotonic velocity fields.  Deposition'
	    WRITE(LUER,*)'tables are interpolated on V, and the RAD_TRANS'
	    WRITE(LUER,*)'solver is the legacy signed-SIGMA formal path.'
	    WRITE(LUER,*)'Use LOCAL, ABS_TRANS, or a radius-mapped table.'
	    STOP
	  END IF
!
! All available clumping laws are parameterized by V (or spline nodes
! on V), so repeated velocities on a non-monotonic V(R) profile force
! the same clumping factor onto physically different radii.  Several
! laws also normalize against V(1), which need not be the outer or
! maximum velocity scale in this mode.  Reject DO_CL rather than
! silently assigning ambiguous clumping.
!
	  IF((DO_NONMONO_SOB .OR. DO_NONMONO_KNOP) .AND. DO_CLUMP_MODEL)THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'DO_CL=T is incompatible with non-monotonic'
	    WRITE(LUER,*)'velocity fields.  The available CL_LAW options'
	    WRITE(LUER,*)'are functions of V and cannot distinguish repeated'
	    WRITE(LUER,*)'velocities at different radii.'
	    WRITE(LUER,*)'Set DO_CL=F.'
	    STOP
	  END IF
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(ACCURATE,'INC_GRID',L_TRUE,
	1          'Increase grid size to improve accuracy? ')
!
! DO_NONMONO_KNOP is incompatible with the enhanced-grid (ACCURATE)
! branch of COMP_J_BLANK, which calls MOM_J_CMF_V11 after the formal
! transfer and is not non-monotonic safe.
!
	  IF(DO_NONMONO_KNOP .AND. ACCURATE)THEN
	    WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'DO_NONMONO_KNOP=T requires INC_GRID=F'
	    WRITE(LUER,*)'The enhanced-grid branch calls MOM_J_CMF_V11 after'
	    WRITE(LUER,*)'formal transfer and is not non-monotonic safe.'
	    STOP
	  END IF
	  CALL RD_STORE_LOG(ALL_FREQ,'ALL_FREQ',L_TRUE,
	1          'Increase accuracy for all frequencies?')
	  CALL RD_STORE_DBLE(ACC_FREQ_END,'ACC_END',L_TRUE,
	1          'Increase accuracy for all frequencies < ACC_END?')
	  CALL RD_STORE_INT(NPINS,'N_INS',L_TRUE,
	1          'Number of points to be inserted in higher'//
	1          ' accuracy grid (1, 2 or 3) ')
	  CALL RD_STORE_INT(ST_INTERP_INDX,'ST_INT',L_TRUE,
	1          'Interpolate from ? ')
	  CALL RD_STORE_INT(END_INTERP_INDX,'END_INT',L_TRUE,
	1          'Interpolate to ? ')
	  CALL RD_STORE_INT(DEEP,'ND_QUAD',L_TRUE,
	1         'Quadratic interpolation from ND-? to ND')
	  CALL RD_STORE_NCHAR(INTERP_TYPE,'INTERP_TYPE',10,L_TRUE,
	1         'Perform interpolations in LOG or LIN plane')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_INT(N_PAR,'N_PAR',L_TRUE,
	1    'Rate of BA incrementation by BA_PAR in cont. loop (# of freq)')
!
! Next two variables apply for both ACCURATE and EDDINGTON.
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(COMPUTE_EDDFAC,'COMP_F',L_TRUE,
	1      'Compute new Eddington factors (f)')
	  CALL RD_STORE_DBLE(ACC_EDD_FAC,'ACC_F',L_TRUE,
	1      'Accuracy with which to compute the eddington factor f')
!
! SUP_FCONV_IT1 (dev7): summarize rather than print the per-frequency
! 'Possible error converging f' warnings during great iteration 1.  The
! iteration-1 EDDFACTOR build starts from FOLD=0 with unrelaxed
! populations, so thousands of marginal warnings are expected and flood
! OUTGEN.  Later-iteration warnings still print individually.  Default
! TRUE; set F in VADAT to restore per-frequency iteration-1 prints.
!
	  SUPPRESS_FCONV_IT1=.TRUE.
	  CALL RD_STORE_LOG(SUPPRESS_FCONV_IT1,'SUP_FCONV_IT1',L_FALSE,
	1      'Summarize iter-1 f-convergence warnings instead of printing?')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(NG_DO,'DO_NG',L_TRUE,
	1         'Perform NG acceleration when applicable ?')
	  CALL RD_STORE_LOG(DO_NG_VALIDITY_CHECK,'CHK_NG',L_TRUE,
	1         'Check if NG accleration has made a reasonable correction to T?')
	  CALL RD_STORE_DBLE(VAL_DO_NG,'BEG_NG',L_TRUE,
	1       'Percentage accuracy at which to begin NG acceleration')
	  CALL RD_STORE_INT(IT_TO_BEG_NG,'IBEG_NG',L_TRUE,
	1       'Iteration at which to begin NG acceleration')
	  CALL RD_STORE_INT(NG_BAND_WIDTH,'BW_NG',L_TRUE,
	1       'Depth band width for NG acceleration')
	  CALL RD_STORE_INT(ITS_PER_NG,'ITS/NG',L_TRUE,
	1         'Number of iterations between NG accelerations (>=4)')
	  IF(ITS_PER_NG .LT. 4)THEN
	     WRITE(LUER,*)'Error in CMFGEN - ITS_PER_NG too small'
	     STOP
	  END IF
!
! NG_ABSTAIN: see CONTROL_VARIABLE_MOD.  Optional (L_FALSE); absent
! key leaves the upstream all-or-nothing veto in place.
!
	  NG_ABSTAIN=.FALSE.
	  CALL RD_STORE_LOG(NG_ABSTAIN,'NG_ABSTAIN',L_FALSE,
	1   'Hold non-contractive/non-positive NG components per-component?')
!
	  WRITE(LUSCR,'()')
	  AVERAGE_DO=.FALSE.; NUM_OSC_AV=4; ITS_PER_AV=8
	  CALL RD_STORE_LOG(AVERAGE_DO,'DO_AV',L_FALSE,'Perform averaging of oscillating pops')
	  CALL RD_STORE_INT(NUM_OSC_AV,'NOSC_AV',L_FALSE,'# of consecquitive oscillations')
	  CALL RD_STORE_INT(ITS_PER_AV,'ITS/AV',L_FALSE,'# of iterations between averaging')
!
!
	  UNDO_LAST_IT=.FALSE.
	  CALL RD_STORE_LOG(UNDO_LAST_IT,'DO_UNDO',L_FALSE,'Undo corrections at last 5 depths')
!
	  STOP_IF_BAD_PARAM=.TRUE.
	  CALL RD_STORE_LOG(STOP_IF_BAD_PARAM,'STOP_IF_BP',L_FALSE,'Stop if bad/inconsistent parameters have been set')
	  STOP_IF_MAJOR_WARNING=.TRUE.
	  CALL RD_STORE_LOG(STOP_IF_MAJOR_WARNING,'STOP_IF_MW',L_FALSE,'Stop if code has potentially a major issue')
!
	  CALL CLEAN_RD_STORE()
!
	CLOSE(UNIT=LUIN)
!
! Check consistency of parameters.
!
	CALL CHECK_PARAM_CONSISTENCY()
!
	END
