	SUBROUTINE RD_CONTROL_VARIABLES(LUIN,LUSCR,LUER,NUM_BNDS)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE JFNK_SOLVER_V1_MOD, ONLY: jfnk_log_unknown, jfnk_ptc_active,
	1     jfnk_ptc_dtau, jfnk_ptc_scale_ready,
	2     jfnk_pc_set_options_v1
	IMPLICIT NONE
!
! Altered : 19-Oct-2026 : Force METH_SOL='DIAG' whenever NUM_BNDS<3:
!                            NUM_BNDS=2 with 'TRI' previously passed
!                            validation and STOPped in CMF_BLKBAND_V3
!                            after the BA was assembled (Devin).
! Altered : 21-Oct-2026 : USE_LAM_ES un-forced: SOB default is F
!                            (MOM_J_CMF_NONMONO_V1 + MOM ALO), explicit
!                            T honoured; KNOP absent-key default stays T
!                            but explicit F is honoured.  New SOB
!                            absent-key auto-sets: DO_NONMONO_FREQ_MATCH
!                            (LAM_ES=F), DO_FG13_ALO (LAM_ES=T),
!                            DO_BA_RES_ACC, DO_BA_RES_REUSE, NG_ABSTAIN.
!                            Explicit-value overrides print a banner
!                            (Devin).
! Altered : 16-Oct-2026 : Self-consistency enforcement: the
!                            non-monotonic/JFNK trigger flags now
!                            auto-configure their dependent switches
!                            (logged to LUER) instead of relying on
!                            VADAT being consistent.  USE_LAM_ES=F and
!                            DO_KNOP_REFRESH=F opt-ins are overridden;
!                            INC_GRID=F is auto-set under KNOP (Devin).
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
	LOGICAL L_FALSE_TMP
	EXTERNAL GET_KEY_STRING
!
	INTEGER I
	INTEGER ISPEC
	INTEGER ID
!
	INTEGER KN_M_MAXPASS_V,KN_RAY_MAXPASS_V,MOM_AA_M_V,
	1     MOM_AA_OUT_M_V,KN_ORD_AA_M_V
	INTEGER KN_RAY_HOTPASSES_V,KN_HOT_TOTAL_PASSES_V
	LOGICAL KN_ADAPT_RELAX_V
	LOGICAL FG_ADAPT_OMEGA_V2_V,FG_PLATEAU_DET_V
	LOGICAL MOM_ADAPT_OMEGA_V2_V
	REAL(KIND=LDP) KN_M_TOL_V,KN_RAY_TOL_V,KN_RAY_BULK_TOL_V
	REAL(KIND=LDP) JNEXT_RELAX_V,INEXT_RELAX_V
	REAL(KIND=LDP) KN_RAY_STALL_V
	REAL(KIND=LDP) FG_COMMIT_JUMP_V,FG_DIV_HEAL_CAP_V
	REAL(KIND=LDP) FG_COMMIT_NC_V
	LOGICAL FG_COMMIT_CONV_V,FG_BULK_CONTR_V
	REAL(KIND=LDP) MOM_DAMP_E
	EXTERNAL GET_MOM_DAMP_MIN
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
!
! Radius-space deposition (dev7, 13-Oct-2026): on a non-monotonic V
! field VLOC_SHOCK_POWER cannot identify a unique shock location, so
! RLOC_SHOCK_POWER (internal R units, 10^10 cm) selects a Gaussian in
! radius instead.  RLOC<=0 keeps the original V-space deposition.
! The strictly-decreasing-V guard in GET_EDEP_SHOCK_POWER applies only
! to the V-space path.
!
	     RLOC_SHOCK_POWER=0.0_LDP
	     CALL RD_STORE_DBLE(RLOC_SHOCK_POWER,'RLOC_SHOCK_POWER',L_FALSE,
	1        'Radius (1e10 cm) for shock power; <=0 uses VLOC')
	     DRLOC_SHOCK_POWER=0.0_LDP
	     CALL RD_STORE_DBLE(DRLOC_SHOCK_POWER,'DRLOC_SHOCK_POWER',L_FALSE,
	1        'Gaussian width (1e10 cm) for radius deposition')
	     IF(RLOC_SHOCK_POWER .GT. 0.0_LDP .AND.
	1       DRLOC_SHOCK_POWER .LE. 0.0_LDP)THEN
	       WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES:'
	       WRITE(LUER,*)'RLOC_SHOCK_POWER > 0 requires a positive',
	1         ' DRLOC_SHOCK_POWER.'
	       STOP
	     END IF
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
	  KN_NO_GAIN=.FALSE.
	  CALL RD_STORE_LOG(KN_NO_GAIN,'KN_NO_GAIN',L_FALSE,
	1     'KNOP: floor CHI at CHI_SCAT (no net photon gain)?')
!
	  DO_NONLOC_SOB=.FALSE.
	  CALL RD_STORE_LOG(DO_NONLOC_SOB,'DO_NONLOC_SOB',L_FALSE,
	1     'Use Rybicki-Hummer non-local Sobolev coupling?')
!
! ALLOW_SIGMA_UNDEF (default FALSE): RECOMPUTE_SIGMA_TRUE aborts
! when a cell has no usable velocity gradient -- on a LAW-14 grid
! that means a collapsed spacing, and continuing would mark the
! cell non-compressing without notice.  TRUE restores the
! warning-only behaviour for diagnosis.
!
	  ALLOW_SIGMA_UNDEF=.FALSE.
	  CALL RD_STORE_LOG(ALLOW_SIGMA_UNDEF,'ALLOW_SIGMA_UNDEF',L_FALSE,
	1     'Warn instead of abort on undefined SIGMA_TRUE cells?')
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
! FG_TRIAL_MIN (default 30, 0 disables): frozen-iteration count
! that arms the FG-only marginal-release trial -- the
! USE_LAM_ES=T analogue of MOM_TRIAL, which can never fire with
! MOM_RAN=F.  Requires the cured gate (no immaturity, no UNSAFE
! share) and the resacc snapshot; accepted trials refund N_FROZEN
! to this floor.
!
	  FG_TRIAL_MIN=30
	  CALL RD_STORE_INT(FG_TRIAL_MIN,'FG_TRIAL_MIN',L_FALSE,
	1     'KNOP: frozen iters before the FG-only LAMBDA trial?')
	  IF(FG_TRIAL_MIN .LT. 0)FG_TRIAL_MIN=0
!
! FG_TRIAL_CHURN (default 0.0 = off): store-churn multiplier under
! which a churn-only immature FG gate still arms the trial.  0 keeps
! the strict no-immaturity requirement; ~3 admits the m15-class
! asymptote where churn sits 2-3x the bulk bound with every other
! gate term clean.
!
	  FG_TRIAL_CHURN=0.0_LDP
	  CALL RD_STORE_DBLE(FG_TRIAL_CHURN,'FG_TRIAL_CHURN',L_FALSE,
	1     'KNOP: churn-bound factor arming the FG-only trial?')
	  IF(FG_TRIAL_CHURN .LT. 0.0_LDP)FG_TRIAL_CHURN=0.0_LDP
!
! FG_TRIAL_CHURN_MIN (default 3): consecutive churn-marginal calls
! before the churn arm may release a trial.  A single
! churn-marginal sweep could still be transient; N_FROZEN is the
! wrong clock because releases refund it under the same hold
! (m15_2p82d: streak ~5-6 vs FG_TRIAL_MIN=30).  0 disables the arm.
!
	  FG_TRIAL_CHURN_MIN=3
	  CALL RD_STORE_INT(FG_TRIAL_CHURN_MIN,'FG_TRIAL_CHURN_MIN',
	1     L_FALSE,'KNOP: churn-marginal calls before the FG trial?')
	  IF(FG_TRIAL_CHURN_MIN .LT. 0)FG_TRIAL_CHURN_MIN=0
!
! GS_SETTLE_TOL (default -1 = off): settle-asymptote tolerance on
! the lagged-field drift.  A 'cycle'-classified argmax floor at or
! below this value stops refreezing the field on store churn alone;
! an absent key keeps the GS_MATURE_TOL behaviour bit-identical.
!
	  GS_SETTLE_TOL=-1.0_LDP
	  CALL RD_STORE_DBLE(GS_SETTLE_TOL,'GS_SETTLE_TOL',L_FALSE,
	1     'KNOP: drift floor at which a cycle floor may settle?')
!
! MOM_TRIAL_LAM (default 0 = disabled): demand ceiling for the
! marginal-release trial.  With MAXCH above it the trial takes the
! damped GS_LAM_CAP LAMBDA step instead of the resacc BA step -- a
! 1E7-scale demand overshoots the BA linearization so far that the
! resacc rejects every release (m15 KNOP warm, Oct-2026).
!
	  MOM_TRIAL_LAM=0.0_LDP
	  CALL RD_STORE_DBLE(MOM_TRIAL_LAM,'MOM_TRIAL_LAM',L_FALSE,
	1     'KNOP: MAXCH above which a trial releases a LAMBDA step?')
	  IF(MOM_TRIAL_LAM .LT. 0.0_LDP)MOM_TRIAL_LAM=0.0_LDP
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
! GS_RELEASE_LIM (default 0 = disabled): rejected trial releases
! before the model stops -- the freeze/release cycle can otherwise
! refund the stall clock indefinitely.
!
	  GS_RELEASE_LIM=0
	  CALL RD_STORE_INT(GS_RELEASE_LIM,'GS_RELEASE_LIM',L_FALSE,
	1     'KNOP: rejected release attempts before stop (0=disabled)?')
	  IF(GS_RELEASE_LIM .LT. 0)GS_RELEASE_LIM=0
!
! GS_REL_NOPROG_LIM (default 0 = disabled): releases that made no
! demand progress or re-poisoned the field before the model stops.
! The count survives release refunds, so it catches the
! release/refreeze limit cycle the consecutive counters cannot.
!
	  GS_REL_NOPROG_LIM=0
	  CALL RD_STORE_INT(GS_REL_NOPROG_LIM,'GS_REL_NOPROG_LIM',
	1     L_FALSE,'KNOP: no-progress releases before stop?')
	  IF(GS_REL_NOPROG_LIM .LT. 0)GS_REL_NOPROG_LIM=0
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
! GS_REL_N (default 2): post-release LAMBDA iterations held at the
! near-mature correction cap while the field re-lags to the moved
! populations.  Default changed to 2 (edited by Devin, 19-Oct-2026):
! with N=1 only the first live step was damped and the it8-class
! undamped second step re-froze the field at factor-100 changes,
! feeding the release->refreeze alternation.
!
	  GS_REL_N=2
	  CALL RD_STORE_INT(GS_REL_N,'GS_REL_N',L_FALSE,
	1     'KNOP: damped LAMBDA its after a freeze release?')
	  IF(GS_REL_N .LT. 1)GS_REL_N=1
!
! GS_RAW_MATURE (default FALSE): when set, a sweep whose raw FG max
! change still exceeds GS_MATURE_TOL does not count toward the
! GS_REL_STABLE release streak, even if the bulk rule passed it.
!
	  GS_RAW_MATURE=L_FALSE
	  CALL RD_STORE_LOG(GS_RAW_MATURE,'GS_RAW_MATURE',L_FALSE,
	1     'KNOP: require raw max change for release maturity?')
!
! GS_PRECOND (default FALSE standalone; auto-set TRUE under
! DO_NONMONO_KNOP in the auto-set block below, explicit F
! respected): postpone the GS freeze while the demanded
! population correction entering the call is grossly
! unconverged.  Holding populations at ~1E10% demand parks the
! iteration 1E9 outside the solver's validity domain
! (m15_2p82d KNOP release/refreeze limit cycle).
!
	  GS_PRECOND=L_FALSE
	  CALL RD_STORE_LOG(GS_PRECOND,'GS_PRECOND',L_FALSE,
	1     'KNOP: postpone freeze while pop demand is unconverged?')
!
! GS_PRECOND_MAX (default 1.0E+03): percent demand threshold for
! the precondition gate.  Entry MAXCH is the previous call's
! demanded maximum correction.
!
	  GS_PRECOND_MAX=1.0E+03_LDP
	  CALL RD_STORE_DBLE(GS_PRECOND_MAX,'GS_PRECOND_MAX',L_FALSE,
	1     'KNOP: demand above which the GS freeze is postponed?')
	  IF(GS_PRECOND_MAX .LE. 0.0_LDP)GS_PRECOND_MAX=1.0E+03_LDP
!
! GS_PRECOND_LIM (default 0 = disabled): consecutive
! preconditioned calls before the model stops as non-convergent
! -- a demand that never falls is a structural failure, not a
! transient.
!
	  GS_PRECOND_LIM=0
	  CALL RD_STORE_INT(GS_PRECOND_LIM,'GS_PRECOND_LIM',L_FALSE,
	1     'KNOP: preconditioned calls before stop (0=disabled)?')
	  IF(GS_PRECOND_LIM .LT. 0)GS_PRECOND_LIM=0
!
! JFNK_CPL_RELWAIT (default FALSE): when set, the coupled JFNK tier
! is deferred while the call sits in the post-release damp window
! (N_REL_SINCE <= GS_REL_N).
!
	  JFNK_CPL_RELWAIT=L_FALSE
	  CALL RD_STORE_LOG(JFNK_CPL_RELWAIT,'JFNK_CPL_RELWAIT',L_FALSE,
	1     'KNOP: hold coupled JFNK during the GS release window?')
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
! DO_KNOP_REFRESH is read later, after USE_LAM_ES: its default
! expression depends on USE_LAM_ES, and Fortran does not guarantee
! short-circuit .OR. evaluation, so evaluating it here read an
! undefined variable on the SOB-only path.
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
! DO_NONMONO_FREQ_MATCH (default FALSE here; the KNOP auto-set
! below lifts it to T when the key is absent): rematch the lagged
! upwind frequency index at compressing points to the bin
! containing the photon's characteristic foot (per-cell drift in
! bins), and widen the frequency-derivative stencil to match.
! Targets the oscillatory-eigenvalue instability identified by
! Ceulemans et al. 2025 when the per-cell Doppler shift exceeds
! the frequency spacing.  Only acts where the Knop solver already
! switches to the lagged redder upwind; expanding cells keep the
! implicit same-sweep coupling.
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
! MOM_GAM_PIN (default TRUE): evaluate the lagged comoving-gradient
! factors in MOM_J_CMF_NONMONO_V1 (GAM_FAC_PREV, GAMH_FAC_PREV, the
! Eddington boundary PSIPREV, and the EPS_PREV angle factor) at the
! current frequency's Eddington factors instead of the bluer
! neighbour's.  Independent lagged evaluations let the explicit
! upwind weight exceed the implicit diagonal where 1+sigma*K drifts
! sign or magnitude across a frequency boundary -- the observed
! m15 lagged-map gain channel (12-Oct-2026 probes).  FALSE restores
! the independent lagged evaluation for A/B tests.
!
	  MOM_GAM_PIN=.TRUE.
	  CALL RD_STORE_LOG(MOM_GAM_PIN,'MOM_GAM_PIN',L_FALSE,
	1     'Pin lagged GAM factors to current-frequency values?')
!
! MOM_LAG_BND (default TRUE): per-row bound on the lagged frequency
! source in MOM_J_CMF_NONMONO_V1 -- where the row's lagged coupling
! weight (upwind J, HS*H_up, EPS_PREV transfers) exceeds the
! diagonal excess TB-|TA|-|TC|, the lagged contribution is scaled
! down so the lagged fixed-point map cannot amplify there.  FALSE
! restores the unbounded source for A/B tests.
!
	  MOM_LAG_BND=.TRUE.
	  CALL RD_STORE_LOG(MOM_LAG_BND,'MOM_LAG_BND',L_FALSE,
	1     'Bound lagged source by per-row diagonal excess?')
!
! MOM_ABS_FLR (default 0.0 = off): floor on the absorption fraction
! 1-COH in the MOM_J_CMF_NONMONO_V1 tridiagonal diagonal.  At
! ESEC>=CHI cells the clip at 0 removes the only J-independent
! diagonal weight; with PSI~0 near the 1+sigma*K=0 resonance the row
! then has no diagonal excess, MOM_LAG_BND nulls its whole lagged
! source, and the row sloshes between clipped and unclipped states
! rather than converging (KNOP probes, 13-Oct-2026).  A small floor
! (~1D-3) restores a physical absorption anchor.  0.0 keeps the
! exact clip for A/B tests.
!
	  MOM_ABS_FLR=0.0_LDP
	  CALL RD_STORE_DBLE(MOM_ABS_FLR,'MOM_ABS_FLR',L_FALSE,
	1     'Absorption-fraction floor in non-mono MOM diagonal')
!
! MOM_H_IMPL (default FALSE): add the lagged upwind-H coupling
! magnitude |VB|+|VC| to the MOM_J_CMF_NONMONO_V1 diagonal.  The
! lagged H divergence is then paid for as implicit weight -- the
! row is damped against both neighbour and lagged reads, and the
! MOM_LAG_BND excess grows by the same amount, so clipped shock
! rows stop cycling between clipped and unclipped fixed points.
!
	  MOM_H_IMPL=.FALSE.
	  CALL RD_STORE_LOG(MOM_H_IMPL,'MOM_H_IMPL',L_FALSE,
	1     'Implicit weight for lagged H coupling in MOM diagonal?')
	  CALL SET_MOM_LAG_OPTS(MOM_GAM_PIN,MOM_LAG_BND,MOM_ABS_FLR,
	1     MOM_H_IMPL)
!
! MOM_FAC_RLX (default 1.0 = off): under-relax the Eddington-factor
! update between great iterations in MOM_J_CMF_NONMONO_V1.  With
! pops held the formal-solve factor snap is the only channel that
! moves the lagged map's fixed point between iterations; on this
! model the snap displaces the map by >50% at ~2E5 cells per sweep
! while the refresh already converges at the old factors (KNOP cold
! probes, 13-Oct-2026).  1.0 keeps the raw snap for A/B tests.
!
	  MOM_FAC_RLX=1.0_LDP
	  CALL RD_STORE_DBLE(MOM_FAC_RLX,'MOM_FAC_RLX',L_FALSE,
	1     'Eddington-factor under-relaxation between MOM sweeps')
	  CALL SET_MOM_FAC_RLX(MOM_FAC_RLX)
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
!
! KN_RAY_BULK_TOL (default 6.0E-2 under DO_NONMONO_KNOP, 0 = off
! otherwise; edited by Devin, 20-Oct-2026): residual floor on the
! bulk exit (IXIT=2).  The bad-fraction test alone let rays exit
! converged at residual ~5E-2 -- one capped commit step from the
! fixed point -- pinning the refresh max just above the maturity
! gate (m15_2p82d: 5.19-5.24E-2) and ratcheting the store through
! WRITE_STORE_CAP.  The earlier 2.5E-2 default sat below that
! floor, so the bulk exit never fired on m15 and floored rays
! still hit the commit gate capped.  6.0E-2 clears the measured
! floor so those rays commit through WRITE_STORE.  The key is
! KNOP-scoped: SET_FG_KN_OPTS is global, so a shared nonzero
! default would silently floor bulk exits on non-KNOP runs --
! off unless the key is set explicitly.  0 disables the term.
!
	  KN_RAY_BULK_TOL_V=0.0_LDP
	  IF(DO_NONMONO_KNOP)KN_RAY_BULK_TOL_V=6.0E-2_LDP
	  CALL RD_STORE_DBLE(KN_RAY_BULK_TOL_V,'KN_RAY_BULK_TOL',
	1     L_FALSE,'Residual floor for the FG bulk exit (0=off)')
	  CALL SET_FG_KN_OPTS(KN_RAY_TOL_V,KN_RAY_MAXPASS_V,
	1     INEXT_RELAX_V,KN_RAY_BULK_TOL_V)
!
! KN_ADAPT_RELAX (default T): adaptive write-back weight for both
! lagged stores.  The relax factor becomes the seed value; each
! refresh call re-estimates the dominant-mode contraction rho_est
! and moves the weight toward 1/(1+rho_est).  Set FALSE to keep the
! fixed INEXT_RELAX/JNEXT_RELAX behaviour.
!
	  KN_ADAPT_RELAX_V=.TRUE.
	  CALL RD_STORE_LOG(KN_ADAPT_RELAX_V,'KN_ADAPT_RELAX',L_FALSE,
	1     'Adaptive lagged-store under-relaxation (FG+MOM)?')
	  CALL SET_FG_KN_ADAPT(KN_ADAPT_RELAX_V)
	  CALL SET_MOM_KN_ADAPT(KN_ADAPT_RELAX_V)
!
! MOM_ADAPT_OMEGA_V2 (default T, 19-Oct-2026): same-map omega
! estimator in the MOM refresh -- the relax weight is re-estimated
! from the median of the last ~4 per-pass residual ratios of the
! current refresh instead of the cross-call ratio, which miskeyed
! oscillation and marginal contraction to the same omega~0.5.
! F restores the cross-call rule (ports FG_ADAPT_OMEGA_V2, W5).
!
	  MOM_ADAPT_OMEGA_V2_V=.TRUE.
	  CALL RD_STORE_LOG(MOM_ADAPT_OMEGA_V2_V,'MOM_ADAPT_OMEGA_V2',
	1     L_FALSE,'Same-map omega estimator in the MOM refresh?')
	  CALL SET_MOM_OMEGA_V2(MOM_ADAPT_OMEGA_V2_V)
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
! admission additionally requires USE_LAM_ES=T (so JFNK is inert on
! the SOB default LAM_ES=F) and, under KNOP only, a mature lagged
! field (GS_FREEZE=F).
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
! Under DO_NONMONO_SOB the coupled tier is opt-in via
! [JFNK_COUPLED_SOB] (default F): the eval sweep then also replays
! the line-frequency section (SUB_SOB_LINE_V3 Sobolev rates against
! the frozen field), which is what makes the coupled residual the
! same equation the outer iteration assembles.  Without the flag,
! JFNK_COUPLED=T under SOB is reset to F with a warning -- never a
! stop -- and the SOB path stays frozen-only.  Even when enabled
! the coupled tier is not demand-admitted under SOB: it fires only
! as the reach-fail fallback inside JFNK_SOLVE_CORR_V1 once the
! frozen tier certifies the residual unreachable.
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
	  JFNK_COUPLED_SOB=.FALSE.
	  CALL RD_STORE_LOG(JFNK_COUPLED_SOB,'JFNK_COUPLED_SOB',L_FALSE,
	1     'Allow coupled JFNK tier under DO_NONMONO_SOB?')
	  JFNK_CPL_SOB_MAX_EVALS=0
	  CALL RD_STORE_INT(JFNK_CPL_SOB_MAX_EVALS,'JFNK_CPL_SOB_MAX_EVALS',
	1     L_FALSE,'Max residual evals per coupled-SOB activation')
	  IF(JFNK_CPL_SOB_MAX_EVALS .LT. 0)JFNK_CPL_SOB_MAX_EVALS=0
	  JFNK_CPL_SOB_EACH=.FALSE.
	  CALL RD_STORE_LOG(JFNK_CPL_SOB_EACH,'JFNK_CPL_SOB_EACH',L_FALSE,
	1     'Recheck coupled-SOB residual consistency each activation?')
	  JFNK_CPL_SOB_TOL=1.0E-5_LDP
	  CALL RD_STORE_DBLE(JFNK_CPL_SOB_TOL,'JFNK_CPL_SOB_TOL',
	1     L_FALSE,'Coupled-vs-production residual consistency tol')
	  IF(JFNK_CPL_SOB_TOL .LE. 0.0_LDP)JFNK_CPL_SOB_TOL=1.0E-5_LDP
!
! The coupled tier under SOB needs JFNK_COUPLED_SOB as the explicit
! opt-in: without it, a JFNK_COUPLED=T request used to STOP the run.
! Auto-F plus a loud warning keeps the run alive on the frozen tier
! -- the default path is unchanged.
!
	  IF(USE_JFNK_SOLVER .AND. JFNK_COUPLED .AND. DO_NONMONO_SOB .AND.
	1     .NOT. JFNK_COUPLED_SOB)THEN
	    WRITE(LUER,*)'Warning in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'JFNK_COUPLED=T under DO_NONMONO_SOB requires'
	    WRITE(LUER,*)'JFNK_COUPLED_SOB=T -- without it the coupled'
	    WRITE(LUER,*)'residual tier stays disabled and SOB admits'
	    WRITE(LUER,*)'only the frozen-J residual.  Forcing'
	    WRITE(LUER,*)'JFNK_COUPLED=F.'
	    JFNK_COUPLED=.FALSE.
	  END IF
	  IF(JFNK_COUPLED_SOB .AND. .NOT. JFNK_COUPLED)THEN
	    WRITE(LUER,*)'Warning in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'JFNK_COUPLED_SOB=T has no effect without'
	    WRITE(LUER,*)'JFNK_COUPLED=T -- forcing JFNK_COUPLED_SOB=F.'
	    JFNK_COUPLED_SOB=.FALSE.
	  END IF
	  IF(JFNK_COUPLED_SOB .AND. .NOT. DO_NONMONO_SOB)THEN
	    WRITE(LUER,*)'Warning in RD_CONTROL_VARIABLES:'
	    WRITE(LUER,*)'JFNK_COUPLED_SOB=T is only meaningful under'
	    WRITE(LUER,*)'DO_NONMONO_SOB -- forcing JFNK_COUPLED_SOB=F.'
	    JFNK_COUPLED_SOB=.FALSE.
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
	  JFNK_MAX_GMRES=50
	  CALL RD_STORE_INT(JFNK_MAX_GMRES,'JFNK_MAXGMRES',L_FALSE,
	1     'Max GMRES its per Newton step (0=use restart dim m)')
!
	  JFNK_GTOL=1.0E-3_LDP
	  CALL RD_STORE_DBLE(JFNK_GTOL,'JFNK_GTOL',L_FALSE,
	1     'Relative GMRES tolerance (fraction of |F|)')
!
	  JFNK_GMRES_STAG=1.0D-3
	  CALL RD_STORE_DBLE(JFNK_GMRES_STAG,'JFNK_GMRES_STAG',L_FALSE,
	1     'GMRES stagnation improvement threshold (0=never stagnate)')
!
	  JFNK_GMRES_RST=1
	  CALL RD_STORE_INT(JFNK_GMRES_RST,'JFNK_GMRES_RST',L_FALSE,
	1     'Max GMRES restarts per solve (safety cap; total its still'
	2     //' bounded by JFNK_MAXGMRES)')
!
	  JFNK_LGMRES=.FALSE.
	  CALL RD_STORE_LOG(JFNK_LGMRES,'JFNK_LGMRES',L_FALSE,
	1     'LGMRES: retain solution corrections across GMRES restarts?')
!
	  JFNK_LGMRES_K=2
	  CALL RD_STORE_INT(JFNK_LGMRES_K,'JFNK_LGMRES_K',L_FALSE,
	1     'Retained correction vectors for LGMRES (1..8)')
	  IF(JFNK_LGMRES_K .LT. 1)JFNK_LGMRES_K=1
	  IF(JFNK_LGMRES_K .GT. 8)JFNK_LGMRES_K=8
!
	  JFNK_GDEFL=.FALSE.
	  CALL RD_STORE_LOG(JFNK_GDEFL,'JFNK_GDEFL',L_FALSE,
	1     'GMRES-DR: deflate harmonic Ritz dirs at each restart?')
!
	  JFNK_GDEFL_K=4
	  CALL RD_STORE_INT(JFNK_GDEFL_K,'JFNK_GDEFL_K',L_FALSE,
	1     'Retained deflation directions per restart (1..8)')
	  IF(JFNK_GDEFL_K .LT. 1)JFNK_GDEFL_K=1
	  IF(JFNK_GDEFL_K .GT. 8)JFNK_GDEFL_K=8
!
	  JFNK_RECYC=.FALSE.
	  CALL RD_STORE_LOG(JFNK_RECYC,'JFNK_RECYC',L_FALSE,
	1     'Carry augmentation set across GMRES solves?')
!
	  JFNK_CYCSTALL=.FALSE.
	  CALL RD_STORE_LOG(JFNK_CYCSTALL,'JFNK_CYCSTALL',L_FALSE,
	1     'Exit GMRES after N flat restart cycles?')
!
	  JFNK_CYCSTALL_FAC=0.02_LDP
	  CALL RD_STORE_DBLE(JFNK_CYCSTALL_FAC,'JFNK_CYCSTALL_FAC',
	1     L_FALSE,
	2     'Min fractional LS-residual drop per cycle before stall')
!
	  JFNK_CYCSTALL_N=2
	  CALL RD_STORE_INT(JFNK_CYCSTALL_N,'JFNK_CYCSTALL_N',L_FALSE,
	1     'Consecutive flat cycles before stall exit')
	  IF(JFNK_CYCSTALL_N .LT. 1)JFNK_CYCSTALL_N=1
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
! JFNK_PRECOND_REFINE is parsed but has no implementation in this
! tree (10-Nov-2026 audit): nothing reads it past this point.  Default
! flipped to F so an explicit T warns instead of pretending the
! historical default did anything; JFNK_SOLVE_CORR_V1 prints the
! inactivity warning once per run.
!
	  JFNK_PRECOND_REFINE=.FALSE.
	  CALL RD_STORE_LOG(JFNK_PRECOND_REFINE,'JFNK_REFINE',L_FALSE,
	1     'Iterative refinement after JFNK precond solve? (INACTIVE)')
!
	  JFNK_LEFT_PRECOND=.TRUE.
	  CALL RD_STORE_LOG(JFNK_LEFT_PRECOND,'JFNK_LEFTPC',L_FALSE,
	1     'Use left GMRES preconditioning (Arramy 2024)?')
!
	  CALL RD_STORE_LOG(JFNK_FADAPT_V,'JFNK_FADAPT',L_FALSE,
	1     'Adaptive Eisenstat-Walker GMRES forcing term (F=fixed gtol)')
!
	  CALL RD_STORE_LOG(JFNK_FD_COMPLEX,'JFNK_FDCPLX',L_FALSE,
	1     'Complex-step JVP mode (NOT IMPLEMENTED - fails closed)')
!
	  CALL RD_STORE_LOG(JFNK_ENRICH_SFR,'JFNK_ESFR',L_FALSE,
	1     'SFR enrichment term (Arramy eqs 55-61; warns, plain path)')
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
! Max clip-crossing retries before keeping the last JVP attempt.
! Each retry shrinks h by 10x; the 21-Oct-2026 h-sweep showed frozen
! FD resolution requires h_eff ~1e-6+, so retries driven by
! degenerate floor crossings can push h below the cancellation
! floor.  JFNK_CLIP_RETRY=0 accepts the crossing (diagnostic use).
!
	  JFNK_CLIP_RETRY=2
	  CALL RD_STORE_INT(JFNK_CLIP_RETRY,'JFNK_CLIP_RETRY',L_FALSE,
	1     'Clip-crossing h retries per JVP (2=default, 0=accept)')
	  IF(JFNK_CLIP_RETRY .LT. 0)JFNK_CLIP_RETRY=0
	  JFNK_JVP_NDIR=4
	  CALL RD_STORE_INT(JFNK_JVP_NDIR,'JFNK_JVP_NDIR',L_FALSE,
	1     'JVPs reporting FD-vs-BA norm line (4=default, 0=off)')
	  IF(JFNK_JVP_NDIR .LT. 0)JFNK_JVP_NDIR=0
	  JFNK_BA_JVP=.FALSE.
	  CALL RD_STORE_LOG(JFNK_BA_JVP,'JFNK_BA_JVP',L_FALSE,
	1     'Analytic BA matvec in place of FD JVP (diagnostic)')
	  JFNK_LINTEST=.FALSE.
	  CALL RD_STORE_LOG(JFNK_LINTEST,'JFNK_LINTEST',L_FALSE,
	1     'JVP operator-consistency probes once per activation')
	  JFNK_HSCAN=.FALSE.
	  CALL RD_STORE_LOG(JFNK_HSCAN,'JFNK_HSCAN',L_FALSE,
	1     'Force the extended JVP h-sweep diagnostic once per run')
	  JFNK_HADAPT=.FALSE.
	  CALL RD_STORE_LOG(JFNK_HADAPT,'JFNK_HADAPT',L_FALSE,
	1     'Pick the FD step from the h-sweep plateau (implies h-sweep)')
	  JFNK_TRES_FAC=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_TRES_FAC,'JFNK_TRES_FAC',L_FALSE,
	1     'GMRES true-residual acceptance factor (0=off)')
	  IF(JFNK_TRES_FAC .LT. 0.0_LDP)JFNK_TRES_FAC=0.0_LDP
	  JFNK_MAX_SEC=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_MAX_SEC,'JFNK_MAX_SEC',L_FALSE,
	1     'Wall-clock seconds per JFNK activation (0=off)')
	  IF(JFNK_MAX_SEC .LT. 0.0_LDP)JFNK_MAX_SEC=0.0_LDP
!
! A zero budget used to warn here; the self-consistency pass below now
! raises JFNK_MAX_SEC=0 to the one-hour default under DO_JFNK and logs
! the override, so a stuck activation cannot consume the whole job.
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
! unit_tests (rel diff < 1D-13 vs the direct solve).
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
! conservative: solvable test systems reach ~1.  Frozen tier only
! (17-Oct-2026): REACH is priced against the stored outer BA, which
! is not the coupled operator -- the coupled tier is exempt.
!
	  JFNK_REACH_GATE=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_REACH_GATE,'JFNK_REACH_GATE',L_FALSE,
	1     'Skip Newton when residual reachable fraction below this')
!
! JFNK_REL_REACH_MAX (default 5): consecutive reach-gate skips that
! stand down the SOB stable-LAMBDA relief.  0 disables the stand-down.
!
	  JFNK_REL_REACH_MAX=5
	  CALL RD_STORE_INT(JFNK_REL_REACH_MAX,'JFNK_REL_REACH_MAX',
	1     L_FALSE,
	2     'SOB: reach skips before stable-LAMBDA relief stands down?')
	  IF(JFNK_REL_REACH_MAX .LT. 0)JFNK_REL_REACH_MAX=0
	  IF(JFNK_REACH_GATE .LT. 0.0_LDP)JFNK_REACH_GATE=0.0_LDP
	  IF(JFNK_REACH_GATE .GT. 1.0_LDP)JFNK_REACH_GATE=1.0_LDP
!
! JFNK direction gates (defaults 0 = off): admission defers to the
! caller's SOLVEBA mode when the preconditioned direction is not a
! usable descent direction.  All three metrics are computed in the
! reachability precheck from z=M^-1 F clipped to the commit bounds:
!   JFNK_COS_MIN       - reject cos(F,A z) below this.  Tested against
!                        both the stored-BA image (free) and the FD
!                        image (one extra residual eval); the FD test
!                        catches a stored operator that reports a
!                        descent the true Jacobian does not share.
!   JFNK_JZ_F_MAX      - reject |A z|/|F| above this (stored and FD).
!   JFNK_CLIP_FRAC_MAX - reject when the fraction of unknowns at the
!                        commit clip bound exceeds this (denominator
!                        JFNK_N, which includes regularized rows).
! A defer uses the same path as the reach gate: caller-mode BA runs
! once and the admission is not counted as stagnation.
!
	  JFNK_COS_MIN=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_COS_MIN,'JFNK_COS_MIN',L_FALSE,
	1     'Defer Newton when cos(F,Jz) below this (0=off)')
	  IF(JFNK_COS_MIN .LE. 0.0_LDP)JFNK_COS_MIN=0.0_LDP
	  IF(JFNK_COS_MIN .GT. 1.0_LDP)JFNK_COS_MIN=1.0_LDP
	  JFNK_JZ_F_MAX=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_JZ_F_MAX,'JFNK_JZ_F_MAX',L_FALSE,
	1     'Defer Newton when |Jz|/|F| exceeds this (0=off)')
	  IF(JFNK_JZ_F_MAX .LT. 0.0_LDP)JFNK_JZ_F_MAX=0.0_LDP
	  JFNK_CLIP_FRAC_MAX=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_CLIP_FRAC_MAX,'JFNK_CLIP_FRAC_MAX',
	1     L_FALSE,'Defer Newton when clipped-step fraction exceeds this')
	  IF(JFNK_CLIP_FRAC_MAX .LE. 0.0_LDP)JFNK_CLIP_FRAC_MAX=0.0_LDP
	  IF(JFNK_CLIP_FRAC_MAX .GT. 1.0_LDP)JFNK_CLIP_FRAC_MAX=1.0_LDP
!
! FD-image descent gates (17-Oct-2026): the stored-BA metrics above
! can pass while the true (finite-difference) operator maps the same
! step anti-parallel to F -- KNOP_jfnkrel it35 passed stored REACH
! 0.735 while cos(F,A_FD z) = -0.23.  JFNK_FD_COS_MIN gates cos on
! the FD image only; JFNK_FD_AUTODEF (default T) defers any
! anti-descent FD image even with all thresholds at 0.  Both arm the
! precheck block, which costs one extra residual eval per admission.
!
	  JFNK_FD_COS_MIN=0.05_LDP
	  CALL RD_STORE_DBLE(JFNK_FD_COS_MIN,'JFNK_FD_COS_MIN',L_FALSE,
	1     'Defer Newton when FD cos(F,Az) below this (0=off)')
	  IF(JFNK_FD_COS_MIN .LT. 0.0_LDP)JFNK_FD_COS_MIN=0.0_LDP
	  IF(JFNK_FD_COS_MIN .GT. 1.0_LDP)JFNK_FD_COS_MIN=1.0_LDP
	  JFNK_FD_AUTODEF=.TRUE.
	  CALL RD_STORE_LOG(JFNK_FD_AUTODEF,'JFNK_FD_AUTODEF',L_FALSE,
	1     'Defer Newton on anti-descent FD image (cos(F,Az)<0)?')
!
! W11 census gate (17-Oct-2026): price the frozen-tier residual's
! reachability from the F(0) vector and the stored BA row norms,
! BEFORE the LU factorization and the reach-precheck JVP evals.  On
! radiation-pinned residuals (deep-H rows carrying ~all of |F| with
! |BArow|/|F| ~ 1e8) this turns a ~78 s LU+CMAT admission into a
! ~1-eval defer.  Frozen tier only: the coupled operator re-solves J
! per eval so the stored-BA census cannot certify its reachability.
! All defaults 0 = off = pre-census behaviour.
!
	  JFNK_CENSUS_STIFF_CUT=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_CENSUS_STIFF_CUT,'JFNK_CENSUS_STIFF_CUT',
	1     L_FALSE,'Census: |BA row|/|F_i| stiffness cutoff (0=off)')
	  IF(JFNK_CENSUS_STIFF_CUT .LT. 0.0_LDP)JFNK_CENSUS_STIFF_CUT=0.0_LDP
	  JFNK_CENSUS_STIFF_FRAC=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_CENSUS_STIFF_FRAC,
	1     'JFNK_CENSUS_STIFF_FRAC',L_FALSE,
	2     'Defer Newton when stiff-row fraction exceeds this (0=off)')
	  IF(JFNK_CENSUS_STIFF_FRAC .LE. 0.0_LDP)
	1     JFNK_CENSUS_STIFF_FRAC=0.0_LDP
	  IF(JFNK_CENSUS_STIFF_FRAC .GT. 1.0_LDP)
	1     JFNK_CENSUS_STIFF_FRAC=1.0_LDP
	  JFNK_CENSUS_PIN_FRAC=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_CENSUS_PIN_FRAC,'JFNK_CENSUS_PIN_FRAC',
	1     L_FALSE,
	2     'Defer Newton when |F| share on regularized rows exceeds this')
	  IF(JFNK_CENSUS_PIN_FRAC .LE. 0.0_LDP)JFNK_CENSUS_PIN_FRAC=0.0_LDP
	  IF(JFNK_CENSUS_PIN_FRAC .GT. 1.0_LDP)JFNK_CENSUS_PIN_FRAC=1.0_LDP
!
! JFNK_DEFER_PROBE (default 0 = off): cap on deferral-path operator
! probes.  An FD-gate deferral normally skips the isolation probes
! entirely; this runs a bounded h-scanned FD-vs-BA comparison on the
! dominant residual rows for the first N deferrals so the cause of a
! deferred admission (stale stored operator vs FD cancellation vs
! genuinely unreachable residual) can be read from one log block.
!
	  JFNK_DEFER_PROBE=0
	  CALL RD_STORE_INT(JFNK_DEFER_PROBE,'JFNK_DEFER_PROBE',L_FALSE,
     1     'Probe FD-vs-BA on this many FD-gate deferrals (0=off)')
	  IF(JFNK_DEFER_PROBE .LT. 0)JFNK_DEFER_PROBE=0
!
! W12 (17-Oct-2026): low-rank field-response enrichment of the
! coupled-tier preconditioner (JVP-probed J - BA_A columns at the
! worst stored-operator anchors, Sherman-Morrison-Woodbury apply).
! Coupled tier only.  Default F = byte-identical to pre-enrichment.
!
	  JFNK_ENRICH=.FALSE.
	  CALL RD_STORE_LOG(JFNK_ENRICH,'JFNK_ENRICH',L_FALSE,
	1     'Enrich coupled-tier preconditioner with JVP columns?')
	  JFNK_ENRICH_K=16
	  CALL RD_STORE_INT(JFNK_ENRICH_K,'JFNK_ENRICH_K',L_FALSE,
	1     'Low-rank enrichment anchor count (clamped 2..64)')
	  IF(JFNK_ENRICH_K .LT. 2)JFNK_ENRICH_K=2
	  IF(JFNK_ENRICH_K .GT. 64)JFNK_ENRICH_K=64
	  JFNK_ENRICH_PERCENT=0.0_LDP
	  CALL RD_STORE_DBLE(JFNK_ENRICH_PERCENT,'JFNK_ENRICH_PERCENT',
	1     L_FALSE,'Cap on enrichment evals, percent of eval budget')
	  IF(JFNK_ENRICH_PERCENT .LE. 0.0_LDP)JFNK_ENRICH_PERCENT=0.0_LDP
	  JFNK_ENRICH_ARN=.FALSE.
	  CALL RD_STORE_LOG(JFNK_ENRICH_ARN,'JFNK_ENRICH_ARN',L_FALSE,
	1     'Arnoldi-basis enrichment instead of coordinate anchors?')
!
! DO_BA_DIR_CHK (default F): diagnostic-only JFNK admission on
! direct-BA iterations.  The demand gates are bypassed so the check
! fires on the first full step, where the direction quality matters
! most; the structural gates (warmup, non-LAMBDA, GS-freeze,
! USE_LAM_ES) still apply.  The activation runs the F(0) identity
! test, the stored-BA reach/direction metrics and the FD-vs-stored-
! operator check, then defers to caller-mode SOLVEBA through the
! info=-2 path -- the real correction still comes from SOLVEBA_V13.
! Requires DO_JFNK: the residual evaluator, stored-BA capture and
! LU/LSQ preconditioner only exist under JFNK arming.
!
	  DO_BA_DIR_CHK=.FALSE.
	  CALL RD_STORE_LOG(DO_BA_DIR_CHK,'DO_BA_DIR_CHK',L_FALSE,
	1     'Run JFNK direction diagnostics on BA steps, then defer')
	  BA_DIRCHK_MAX=1
	  JFNK_DIRCHK_N=0
	  CALL RD_STORE_INT(BA_DIRCHK_MAX,'BA_DIRCHK_MAX',L_FALSE,
	1     'Max dirchk defer events before JFNK admission resumes'//
	2     ' (default 1; -1 = every eligible call)')
	  IF(DO_BA_DIR_CHK .AND. .NOT. USE_JFNK_SOLVER)
	1    WRITE(LUER,*)'Warning: DO_BA_DIR_CHK=T has no effect'//
	2    ' without DO_JFNK=T.'
!
! JFNK_REARM_EVERY is retired (17-Oct-2026): the pause/re-arm
! machinery was replaced by the always-try budget ladder.  The key
! is still accepted so old VADAT files parse, but it does nothing.
!
	  JFNK_REARM_EVERY=10
	  CALL RD_STORE_INT(JFNK_REARM_EVERY,'JFNK_REARM_EVERY',L_FALSE,
	1     'RETIRED: pause/re-arm removed; JFNK budget ladder '//
	2     'trims costs instead')
	  CALL GET_KEY_STRING(TMP_STRING,'JFNK_REARM_EVERY',L_FALSE,
	1     L_FALSE_TMP)
	  IF(L_FALSE_TMP)
	1    WRITE(LUER,*)'Warning: JFNK_REARM_EVERY is retired; the '//
	2    'budget ladder (JFNK_RUNG_*) replaces pause/re-arm.'
!
! Budget ladder rungs (17-Oct-2026): JFNK is attempted on every
! admissible call.  Each counted failure raises JFNK_BUDGET_RUNG by
! one (max JFNK_RUNG_MAX); a committed outcome resets it.  Per rung
! the activation budgets halve (evals third) toward these floors:
! nits -> MAX(MIN_NITS,nits/2), maxgmres -> MAX(MIN_GMRES,m/2),
! evals -> MAX(MIN_EVALS,e/3), restarts -> 1 at rung >= 2.
! Defaults reproduce the pre-ladder per-activation clamps.
!
	  JFNK_RUNG_MAX=4
	  CALL RD_STORE_INT(JFNK_RUNG_MAX,'JFNK_RUNG_MAX',L_FALSE,
	1     'Top JFNK budget rung (budgets trimmed per rung; >=1)')
	  IF(JFNK_RUNG_MAX .LT. 1)JFNK_RUNG_MAX=1
	  JFNK_RUNG_MIN_NITS=4
	  CALL RD_STORE_INT(JFNK_RUNG_MIN_NITS,'JFNK_RUNG_MIN_NITS',
	1     L_FALSE,'Floor on Newton iterations per activation')
	  IF(JFNK_RUNG_MIN_NITS .LT. 1)JFNK_RUNG_MIN_NITS=1
	  JFNK_RUNG_MIN_GMRES=15
	  CALL RD_STORE_INT(JFNK_RUNG_MIN_GMRES,'JFNK_RUNG_MIN_GMRES',
	1     L_FALSE,'Floor on GMRES iterations per Newton step')
	  IF(JFNK_RUNG_MIN_GMRES .LT. 1)JFNK_RUNG_MIN_GMRES=1
	  JFNK_RUNG_MIN_EVALS=12
	  CALL RD_STORE_INT(JFNK_RUNG_MIN_EVALS,'JFNK_RUNG_MIN_EVALS',
	1     L_FALSE,'Floor on residual evals per activation')
	  IF(JFNK_RUNG_MIN_EVALS .LT. 1)JFNK_RUNG_MIN_EVALS=1
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
! Non-monotone (Li-Fukushima) reference for the JFNK line search:
! the sufficient-decrease test uses the max |F| over the last
! JFNK_LS_NMWIN accepted iterates instead of the current |F|, which
! keeps the search working when the coupled residual carries sweep
! noise.  Default T / window 4.
!
	  JFNK_LS_NM=.TRUE.
	  CALL RD_STORE_LOG(JFNK_LS_NM,'JFNK_LS_NM',L_FALSE,
	1     'Non-monotone (Li-Fukushima) JFNK line search?')
	  JFNK_LS_NMWIN=4
	  CALL RD_STORE_INT(JFNK_LS_NMWIN,'JFNK_LS_NMWIN',L_FALSE,
	1     'Non-monotone LS reference window (accepted steps)')
!
	  JFNK_CAUCHY=.TRUE.
	  CALL RD_STORE_LOG(JFNK_CAUCHY,'JFNK_CAUCHY',L_FALSE,
	1     'Use Cauchy (M^-1 F) fallback direction?')
!
	  JFNK_USE_BANDED=.FALSE.
	  CALL RD_STORE_LOG(JFNK_USE_BANDED,'JFNK_USE_BANDED',L_FALSE,
	1     'Use banded BA_rad and banded LU preconditioner?')
	  JFNK_BANDED_TDMA=.FALSE.
	  CALL RD_STORE_LOG(JFNK_BANDED_TDMA,'JFNK_BANDED_TDMA',L_FALSE,
	1     'Factor banded preconditioner by block-Thomas (not DGBTRF)?')
!
	  JFNK_SELECTIVE_JVAR=.FALSE.
	  CALL RD_STORE_LOG(JFNK_SELECTIVE_JVAR,'JFNK_SELECTIVE_JVAR',
	1     L_FALSE,'Enable selective J-variation (cross-depth coupling)?')
	  IF(JFNK_SELECTIVE_JVAR)
	1    WRITE(LUER,*)'Warning: JFNK_SELECTIVE_JVAR is read but not'//
	2    ' wired in this tree -- the selective-JVAR mask'//
	3    ' (V6 MASK_TX_TVX_NONMONO_V1) was never ported.  No effect.'
!
! JFNK_ROW_TINY: rows of a per-depth preconditioner block with
! max|A_ij| < JFNK_ROW_TINY*blockmax are replaced by the identity row
! before equilibration.
!
	  JFNK_ROW_TINY=1.0D-30
	  CALL RD_STORE_DBLE(JFNK_ROW_TINY,'JFNK_ROW_TINY',L_FALSE,
	1     'Near-null row threshold for JFNK preconditioner blocks')
!
! Preconditioner audit keys (jfnk-precond-fixes).  All default to
! historical behaviour; the jfnk_pc_* solver-module mirrors are set
! immediately below so the extracted block routines read them without
! a CONTROL_VARIABLE_MOD dependency.
!
	  JFNK_PC_SIGN_DIAG=.FALSE.
	  CALL RD_STORE_LOG(JFNK_PC_SIGN_DIAG,'JFNK_PC_SIGN_DIAG',
	1     L_FALSE,'Log PTC sign diagnostics at each factorization?')
	  JFNK_PTC_SIGNED=.FALSE.
	  CALL RD_STORE_LOG(JFNK_PTC_SIGNED,'JFNK_PTC_SIGNED',L_FALSE,
	1     'Sign-match the PTC diagonal shift to the BA diagonal?')
	  JFNK_REG_ZERO_COL=.FALSE.
	  CALL RD_STORE_LOG(JFNK_REG_ZERO_COL,'JFNK_REG_ZERO_COL',
	1     L_FALSE,'Also zero the column of regularized rows?')
	  JFNK_ROW_TINY_EQUIL=.FALSE.
	  CALL RD_STORE_LOG(JFNK_ROW_TINY_EQUIL,'JFNK_ROW_TINY_EQUIL',
	1     L_FALSE,'Classify near-null rows on the equilibrated block?')
	  JFNK_PC_RCOND_GATE=.FALSE.
	  CALL RD_STORE_LOG(JFNK_PC_RCOND_GATE,'JFNK_PC_RCOND_GATE',
	1     L_FALSE,'Gate factorizations on rcond instead of the x=1 test?')
	  JFNK_PC_RCOND_MIN=1.0D-13
	  CALL RD_STORE_DBLE(JFNK_PC_RCOND_MIN,'JFNK_PC_RCOND_MIN',
	1     L_FALSE,'rcond threshold for JFNK_PC_RCOND_GATE')
	  JFNK_PC_ITRY_CAP=3
	  CALL RD_STORE_INT(JFNK_PC_ITRY_CAP,'JFNK_PC_ITRY_CAP',
	1     L_FALSE,'Ridge-ladder cap per depth under RCOND_GATE')
	  JFNK_PC_PAR_FACTOR=.FALSE.
	  CALL RD_STORE_LOG(JFNK_PC_PAR_FACTOR,'JFNK_PC_PAR_FACTOR',
	1     L_FALSE,'Depth-parallel preconditioner factorization?'//
	2     ' (INACTIVE: parsed, not implemented)')
!
! The solver-module mirrors sync through a setter: their names
! collide with the CONTROL_VARIABLE_MOD declarations in this unit
! (ifort #6405).
	  CALL JFNK_PC_SET_OPTIONS_V1(JFNK_PTC_SIGNED,
	1     JFNK_REG_ZERO_COL, JFNK_ROW_TINY_EQUIL,
	2     JFNK_PC_RCOND_GATE, JFNK_PC_RCOND_MIN, JFNK_PC_ITRY_CAP)
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
! KN_ORD_AA (default FALSE): safeguarded Anderson acceleration of
! the ordered refresh march, replacing the KN_ORD_AITKEN
! experiment (depth-1 secants diverged on the seed transient).
! KN_ORD_AA_M is the secant depth; pairs cost 2*NI*NCF*8 B per
! slot and M+1 slots are kept per slow ray.
!
	  KN_ORD_AA=.FALSE.
	  CALL RD_STORE_LOG(KN_ORD_AA,'KN_ORD_AA',L_FALSE,
	1     'Anderson-accelerate the ordered FG refresh march?')
	  KN_ORD_AA_M_V=2
	  CALL RD_STORE_INT(KN_ORD_AA_M_V,'KN_ORD_AA_M',L_FALSE,
	1     'Anderson secant depth for the ordered FG march')
	  CALL SET_FG_ORD_AA(KN_ORD_AA,KN_ORD_AA_M_V)
!
! KN_ORD_P2DAMP (default FALSE): per-cell period-2 collapse in the
! ordered refresh march.  After KN_ORD_P2N consecutive update-sign
! flips at a significant cell, the commit becomes the bracket
! midpoint until KN_ORD_P2REL flip-free passes release it.  The
! fixed point is unchanged -- only the relaxation path to it.
!
	  KN_ORD_P2DAMP=.FALSE.
	  CALL RD_STORE_LOG(KN_ORD_P2DAMP,'KN_ORD_P2DAMP',L_FALSE,
	1     'Damp per-cell period-2 modes in the ordered FG march?')
	  KN_ORD_P2N=2
	  CALL RD_STORE_INT(KN_ORD_P2N,'KN_ORD_P2N',L_FALSE,
	1     'Consecutive sign flips before a cell is damped')
	  KN_ORD_P2REL=3
	  CALL RD_STORE_INT(KN_ORD_P2REL,'KN_ORD_P2REL',L_FALSE,
	1     'Flip-free passes before a damped cell is released')
	  CALL SET_FG_ORD_P2DAMP(KN_ORD_P2DAMP,KN_ORD_P2N,KN_ORD_P2REL)
!
! KN_GATE_DIVERGENCE (default FALSE, prototype): swap the absolute
! GS_UNSAFE_TOL unsafe tier in the SOLVE_FOR_POPS maturity gate for
! a divergence test on the refresh max residual -- unsafe only on
! non-finite values or growth above KN_GATE_DIV_FAC versus the
! value KN_GATE_DIV_M refresh sweeps back.  Dormant by design:
! enabling changes gate policy and needs owner sign-off.
!
	  KN_GATE_DIVERGENCE=.FALSE.
	  CALL RD_STORE_LOG(KN_GATE_DIVERGENCE,'KN_GATE_DIVERGENCE',
	1     L_FALSE,'Gate unsafe tier tests divergence, not size?')
	  KN_GATE_DIV_FAC=1.5_LDP
	  CALL RD_STORE_DBLE(KN_GATE_DIV_FAC,'KN_GATE_DIV_FAC',L_FALSE,
	1     'Growth factor over KN_GATE_DIV_M sweeps that is unsafe')
	  IF(KN_GATE_DIV_FAC .LT. 1.0_LDP)KN_GATE_DIV_FAC=1.0_LDP
	  KN_GATE_DIV_M=4
	  CALL RD_STORE_INT(KN_GATE_DIV_M,'KN_GATE_DIV_M',L_FALSE,
	1     'Sweeps of maxchg history in the divergence test')
	  KN_GATE_DIV_M=MAX(1,MIN(KN_GATE_DIV_M,8))
!
! KN_RAY_STALL (default 0.98): same-direction residual ratio above
! which a pass counts toward the NSTALL early exit.  Fires only on
! nearly flat passes so slowly contracting rays use more of
! KN_RAY_MAXPASS; 0.90 quit rays that were still contracting and
! held nearly all unconverged cells in the m15 2.82d validation.
!
	  KN_RAY_STALL_V=0.98_LDP
	  CALL RD_STORE_DBLE(KN_RAY_STALL_V,'KN_RAY_STALL',L_FALSE,
	1     'Same-direction residual ratio counting as a stall pass')
	  CALL SET_FG_STALL_LIM(KN_RAY_STALL_V)
!
! W1-W9 VADATs (17-Oct-2026, KNOP maturity-chain repair):
!
! FG_ADAPT_OMEGA_V2 (default T): same-map omega estimator -- the
! relaxation weight is re-estimated from the per-pass residual
! sequence of the current refresh instead of the cross-iteration
! ratio, which miskeyed oscillation and marginal contraction to
! the same omega~0.5.  F restores the cross-iteration rule.
!
	  FG_ADAPT_OMEGA_V2_V=.TRUE.
	  CALL RD_STORE_LOG(FG_ADAPT_OMEGA_V2_V,'FG_ADAPT_OMEGA_V2',
	1     L_FALSE,'Same-map omega estimator in the FG refresh?')
!
! FG_PLATEAU_DET (default T): plateau and combined-direction stall
! triggers in the refresh sweeps -- three flat same-direction
! ratios, or two straight >0.95 ratios, exit the ray at ~8-10
! passes instead of riding KN_RAY_MAXPASS.
!
	  FG_PLATEAU_DET_V=.TRUE.
	  CALL RD_STORE_LOG(FG_PLATEAU_DET_V,'FG_PLATEAU_DET',L_FALSE,
	1     'Plateau/combined-direction stall exits in FG refresh?')
	  CALL SET_FG_OMEGA_V2(FG_ADAPT_OMEGA_V2_V,FG_PLATEAU_DET_V)
!
! KN_RAY_HOTPASSES (default 120) / KN_HOT_TOTAL_PASSES (default
! 480): extra per-ray pass budget when a stall or cap exit leaves
! residual only inside a contiguous hot-frequency run; the run
! continues restricted to its window +/-2 bins.  The total bounds
! the refresh-wide spend so a hot-ray farm cannot dominate.
!
	  KN_RAY_HOTPASSES_V=120
	  CALL RD_STORE_INT(KN_RAY_HOTPASSES_V,'KN_RAY_HOTPASSES',
	1     L_FALSE,'Extra pass budget for hot-frequency runs per ray')
	  KN_HOT_TOTAL_PASSES_V=480
	  CALL RD_STORE_INT(KN_HOT_TOTAL_PASSES_V,'KN_HOT_TOTAL_PASSES',
	1     L_FALSE,'Refresh-wide cap on hot-run extra passes')
	  CALL SET_FG_HOTPASS(KN_RAY_HOTPASSES_V,KN_HOT_TOTAL_PASSES_V)
!
! FG_COMMIT_JUMP (default 50.0): whole-ray commit rollback -- a
! refreshed store whose max |V1-V2|/max(V2,1e-4*peak) exceeds this
! factor is rejected and the pre-refresh values stay.  Bounds the
! committed-near-blowup class (audit: retained delta 4.19E-1 ->
! 5.03E3) while allowing genuine factor-of-few rewrites.
! FG_DIV_HEAL_CAP (default 1e6): diverged-ray heal bound -- stored
! slots above this multiple of the per-depth median are zeroed so
! a poisoned slot cannot re-seed the next refresh.
!
	  FG_COMMIT_JUMP_V=50.0_LDP
	  CALL RD_STORE_DBLE(FG_COMMIT_JUMP_V,'FG_COMMIT_JUMP',L_FALSE,
	1     'Max committed-write relative jump before ray rollback')
	  FG_DIV_HEAL_CAP_V=1.0E+6_LDP
	  CALL RD_STORE_DBLE(FG_DIV_HEAL_CAP_V,'FG_DIV_HEAL_CAP',
	1     L_FALSE,'Heal bound on diverged-ray poisoned store slots')
!
! FG_COMMIT_CONV (17-Oct-2026, default T): fail-closed refresh
! commits.  When T, a ray whose refresh solve exits on stall,
! pass-cap or plateau (IXIT 3/4/6, residual still above KN_RAY_TOL)
! does not write its last iterate to the lagged store -- the
! pre-refresh generation stays as the seed instead of committing a
! non-converged field bounded only by FG_COMMIT_JUMP (w19 class:
! cap-exit commit at residual 0.465, retained delta 21).
! Default changed to T (edited by Devin, 19-Oct-2026): the guard
! only closed the w19 drift class when the VADAT key was set by
! hand, so fail-open was the wrong default.  F remains available
! as an opt-in diagnostic.
!
	  FG_COMMIT_CONV_V=.TRUE.
	  CALL RD_STORE_LOG(FG_COMMIT_CONV_V,'FG_COMMIT_CONV',L_FALSE,
	1     'Reject refresh commits on non-converged ray exits?')
!
! FG_COMMIT_NC (default 0.05): under FG_COMMIT_CONV a
! non-converged ray (bad exit class OR residual > KN_RAY_TOL --
! bulk exits are also caught, since IXIT=2 fires on the bad-cell
! fraction alone) takes a capped provisional commit, each cell
! bounded to NC times its retained value.  w19's bulk-exit commits
! of 8-27x at residual ~0.47 were the luminosity ratchet.
! Default changed to 0.05 (edited by Devin, 19-Oct-2026): the NC=0
! pure-reject default starved the store on models where every ray
! exits un-converged, re-arming the cold-reseed cycle the commit
! gate was built to break.
! FG_BULK_CONTR (default T) additionally requires same-direction
! contraction before the bulk exit may fire, so flat noncontracting
! rays are classed stall/plateau/cap instead of converged.
!
	  FG_COMMIT_NC_V=0.05_LDP
	  CALL RD_STORE_DBLE(FG_COMMIT_NC_V,'FG_COMMIT_NC',L_FALSE,
	1     'KNOP: per-cell cap factor on non-converged commits (0=reject)')
	  FG_BULK_CONTR_V=.TRUE.
	  CALL RD_STORE_LOG(FG_BULK_CONTR_V,'FG_BULK_CONTR',L_FALSE,
	1     'Require contraction before the bulk bad-fraction exit?')
	  CALL SET_FG_COMMIT(FG_COMMIT_JUMP_V,FG_DIV_HEAL_CAP_V,
	1     FG_COMMIT_CONV_V,FG_COMMIT_NC_V,FG_BULK_CONTR_V)
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
	  BA_LQ_DEM=1.0D+12
	  CALL RD_STORE_DBLE(BA_LQ_DEM,'BA_LQ_DEM',L_FALSE,
	1     'LAMBDA-phase quarantine demand bound (%) -- a row past'//
	2     ' this for BA_LQ_PAT consecutive LAMBDA calls is parked')
	  BA_LQ_PAT=8
	  CALL RD_STORE_INT(BA_LQ_PAT,'BA_LQ_PAT',L_FALSE,
	1     'Consecutive over-demand LAMBDA calls before an eqn parks')
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
	  BA_DEM_PFLOOR=-1.0_LDP
	  CALL RD_STORE_DBLE(BA_DEM_PFLOOR,'BA_DEM_PFLOOR',L_FALSE,
	1     'Population floor for the BA demanded-corr metric (-1=auto)')
	  BA_MAXCH_FAIL=1.0E+30_LDP
	  CALL RD_STORE_DBLE(BA_MAXCH_FAIL,'BA_MAXCH_FAIL',L_FALSE,
	1     'Reject and roll back a committed BA step above this change')
	  BA_TR_MAXSWP=3
	  CALL RD_STORE_INT(BA_TR_MAXSWP,'BA_TR_MAXSWP',L_FALSE,
	1     'Stage-B trace-polish sweeps after majors converge (0=off)')
	  BA_TR_ABSMAX=1.0E-08_LDP
	  CALL RD_STORE_DBLE(BA_TR_ABSMAX,'BA_TR_ABSMAX',L_FALSE,
	1     'Trace convergence: max |dn|/POP_ATOM over exempted rows')
	  BA_TRACE_PHASE=.FALSE.
	  BA_TR_NSWP=0
	  BA_TR_ABS=0.0_LDP
	  BA_DEM_EXM=0.0_LDP
	  BA_EX_RELFRAC=-1.0_LDP
	  CALL RD_STORE_DBLE(BA_EX_RELFRAC,'BA_EX_RELFRAC',L_FALSE,
	1     'Rel. floor for the converged-exit demand metric (-1=off)')
	  BA_EX_MAXCH=1.0E+30_LDP
	  BA_DEM_EXR=0.0_LDP
	  DO_FG13_ALO=.FALSE.
	  CALL RD_STORE_LOG(DO_FG13_ALO,'DO_FG13_ALO',L_FALSE,
	1     'Use FG_J_CMF_V13 diagonal ALO in the SOB continuum Jacobian')
!
! Residual-driven BA recall (default FALSE): while the judged
! residual still contracts, allow BAMAT recall past the FIX_BA
! demand window.  Inert without DO_BA_RES_ACC, which supplies the
! residual measurements; warn rather than alter the flag so the
! VADAT echo stays truthful.
!
	  DO_BA_RES_REUSE=.FALSE.
	  CALL RD_STORE_LOG(DO_BA_RES_REUSE,'DO_BA_RES_REUSE',L_FALSE,
	1     'Recall stored BA matrix while resacc steps still'//
	2     ' contract the SE residual')
	  BA_REUSE_TOL=0.9_LDP
	  CALL RD_STORE_DBLE(BA_REUSE_TOL,'BA_REUSE_TOL',L_FALSE,
	1     'RES2 ratio above which a judged recall step counts'//
	2     ' as stalled')
	  IF(BA_REUSE_TOL .LT. 0.0_LDP)BA_REUSE_TOL=0.0_LDP
	  IF(BA_REUSE_TOL .GT. 1.0_LDP)BA_REUSE_TOL=1.0_LDP
	  BA_REUSE_NSTALL=2
	  CALL RD_STORE_INT(BA_REUSE_NSTALL,'BA_REUSE_NSTALL',L_FALSE,
	1     'Consecutive stalled judged steps before forcing a'//
	2     ' fresh BA matrix')
	  IF(BA_REUSE_NSTALL .LT. 1)BA_REUSE_NSTALL=1
	  BA_JUDGED=.FALSE.
	  BA_JUDGE_RATIO=1.0_LDP
	  IF(DO_BA_RES_REUSE .AND. .NOT. DO_BA_RES_ACC)
	1    WRITE(LUER,*)'Warning: DO_BA_RES_REUSE=T has no effect'//
	2    ' without DO_BA_RES_ACC.'
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
	  CALL GET_KEY_STRING(TMP_STRING,'BA_LUM_GATE',L_FALSE,
	1     L_FALSE_TMP)
	  IF(L_FALSE_TMP)
	1    WRITE(LUER,*)'Warning: BA_LUM_GATE is a no-op; the '//
	2    'inner/outer luminosity gate is diagnostic only (shock '//
	3    'deposition).  Runaway detection is LUM_DIV_*.'
	  CALL GET_KEY_STRING(TMP_STRING,'BA_LUM_DCHG',L_FALSE,
	1     L_FALSE_TMP)
	  IF(L_FALSE_TMP)
	1    WRITE(LUER,*)'Warning: BA_LUM_DCHG is a no-op; see '//
	2    'LUM_DIV_FAC/LUM_DIV_N/LUM_DIV_CAP instead.'
!
! LUM_DIV_* (17-Oct-2026, DO_NONMONO_KNOP only): sustained outer-
! luminosity runaway detector.  Inner/outer agreement stays a
! diagnostic (shock deposition makes it physically ambiguous), but
! factor-of-LUM_DIV_FAC growth for LUM_DIV_N consecutive calls is
! unambiguous divergence; LUM_DIV_CAP (L_sun) is the nonphysical
! ceiling past which the run stops as non-convergent.
! LUM_DIV_N=0 disables the guard.
!
	  LUM_DIV_FAC=1.0D2
	  CALL RD_STORE_DBLE(LUM_DIV_FAC,'LUM_DIV_FAC',L_FALSE,
	1     'KNOP: per-call outer-lum growth factor counted as runaway')
	  IF(LUM_DIV_FAC .LE. 1.0_LDP)LUM_DIV_FAC=1.0D2
	  LUM_DIV_N=3
	  CALL RD_STORE_INT(LUM_DIV_N,'LUM_DIV_N',L_FALSE,
	1     'KNOP: consecutive runaway calls before LUM_BLOCKED (0=off)')
	  IF(LUM_DIV_N .LT. 0)LUM_DIV_N=3
	  LUM_DIV_CAP=1.0D15
	  CALL RD_STORE_DBLE(LUM_DIV_CAP,'LUM_DIV_CAP',L_FALSE,
	1     'KNOP: L_out (L_sun) past which the run stops divergent')
	  IF(LUM_DIV_CAP .LE. 0.0_LDP)LUM_DIV_CAP=1.0D15
!
! LUM_DIV_FAC2/LUM_DIV_N2 (17-Oct-2026): second, lower-factor streak.
! The FAC=1E2/N=3 detector misses a sustained geometric runaway of
! ~10x per call (w19: 9.5E9 -> 5.9E14 L_sun in 7 iterations) -- the
! slow tier catches LUM_DIV_N2 consecutive calls each growing by
! more than LUM_DIV_FAC2.  LUM_DIV_REC (default 1.0, 0 = never
! release) is the latch-recovery factor: once LUM_BLOCKED trips it
! holds until L_out falls back below LUM_DIV_REC times the
! pre-runaway level, so a single sub-threshold call cannot clear
! the hold at a still-runaway luminosity.
!
	  LUM_DIV_FAC2=5.0_LDP
	  CALL RD_STORE_DBLE(LUM_DIV_FAC2,'LUM_DIV_FAC2',L_FALSE,
	1     'KNOP: lower per-call growth factor for the slow streak')
	  IF(LUM_DIV_FAC2 .LE. 1.0_LDP)LUM_DIV_FAC2=5.0_LDP
	  LUM_DIV_N2=4
	  CALL RD_STORE_INT(LUM_DIV_N2,'LUM_DIV_N2',L_FALSE,
	1     'KNOP: consecutive slow-growth calls before LUM_BLOCKED')
	  IF(LUM_DIV_N2 .LT. 0)LUM_DIV_N2=4
	  LUM_DIV_REC=1.0_LDP
	  CALL RD_STORE_DBLE(LUM_DIV_REC,'LUM_DIV_REC',L_FALSE,
	1     'KNOP: release lum latch below REC*pre-runaway L_out (0=never)')
	  IF(LUM_DIV_REC .LT. 0.0_LDP)LUM_DIV_REC=1.0_LDP
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
	  IF(NUM_BNDS .LT. 3 .AND. METH_SOL .NE. 'DIAG')THEN
	    WRITE(LUER,*)'****************************************'
	    WRITE(LUER,*)'******WARNING in CMFGEN*****************'
	    WRITE(LUER,*)'Solution method requires NUM_BNDS of at least 3'
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
!
! edited by Devin, 17-Oct-2026: MAX_LAM/MAX_LIN sanity checks.  Values
! <= 1 forbid any increase and deadlock the solver; MAX_LAM ~ 1D10
! (used on the KNOP_jfnkrel run) makes the LAMBDA cap unbounded on
! every path GS_LAM_CAP does not clip - warn loudly rather than stop,
! since legacy pipelines set it high on purpose.  MAX_LAM < MAX_LIN
! inverts the documented intent (LAMBDA tolerates the larger step).
!
	  IF(MAX_LIN_COR .LE. 1.0_LDP .OR. MAX_LAM_COR .LE. 1.0_LDP)THEN
	    WRITE(LUER,*)' Error: MAX_LIN/MAX_LAM in VADAT must exceed 1.0'
	    WRITE(LUER,*)' MAX_LIN=',MAX_LIN_COR,' MAX_LAM=',MAX_LAM_COR
	    STOP
	  END IF
	  IF(MAX_LAM_COR .GT. 1.0D4)THEN
	    WRITE(LUER,*)'******** WARNING in CMFGEN *****************'
	    WRITE(LUER,*)'MAX_LAM=',MAX_LAM_COR,' is effectively '//
	1    'unbounded; a physically bounded value (1D2-1D3) is '//
	2    'recommended for production'
	  END IF
	  IF(MAX_LAM_COR .LT. MAX_LIN_COR)THEN
	    WRITE(LUER,*)'******** WARNING in CMFGEN *****************'
	    WRITE(LUER,*)'MAX_LAM=',MAX_LAM_COR,' < MAX_LIN=',
	1    MAX_LIN_COR,'; LAMBDA steps are meant to allow the '//
	2    'larger change'
	  END IF
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
	  LAM_STABLE_PROVEN=.FALSE.
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
!
! A user cap <= 1 lets MIN(T1,BA_DAMP_CAP) collapse every damped
! step to a contraction; fail at read time instead of letting the
! schedule reach a meaningless state (same contract as the MAX_dT
! validation above).
!
	  IF(BA_DAMP_CAP .GT. 0.0_LDP .AND. BA_DAMP_CAP .LE. 1.0_LDP)THEN
	    WRITE(LUER,*)' Error: BA_DAMP_CAP in VADAT has an invalid'//
	1      ' value of',BA_DAMP_CAP
	    WRITE(LUER,*)' Require BA_DAMP_CAP = 0 (off) or > 1.0'
	    STOP
	  END IF
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
!
! DO_KNOP_REFRESH (default DO_NONMONO_KNOP): per-ray re-solve of the
! lagged compressing-branch coupling after each frequency sweep (Knop
! Paper II Sect. 3.2 ordering), replacing the one-bin-per-iteration
! lag that cannot mature on strongly non-monotonic models.  The
! self-consistency pass below forces this flag T under
! DO_NONMONO_KNOP; under DO_NONMONO_SOB it serves the moment path
! (the SOB default USE_LAM_ES=F runs MOM_J_CMF_NONMONO_V1, whose
! lagged J_NEXT map needs the same within-sweep convergence).
! Read here, after USE_LAM_ES, so the default expression
! sees the final value.  The MOM-side wiring (MOM_RF_ENABLED) is
! re-synced on the post-read value so the earlier placeholder call
! does not leave the two sides disagreeing.
!
	  DO_KNOP_REFRESH=DO_NONMONO_KNOP .OR.
	1    (DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES)
	  CALL RD_STORE_LOG(DO_KNOP_REFRESH,'DO_KNOP_REFRESH',L_FALSE,
	1     'Per-ray converge of KNOP lagged I_NEXT coupling?')
	  CALL SET_FG_KNOP_REFRESH(DO_KNOP_REFRESH)
	  CALL SET_MOM_RF_DIR_OPT(MOM_RF_DIR,MOM_EPS_LIM,
	1     DO_KNOP_REFRESH,MOM_RF_OWNW)
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
! ------------------------------------------------------------------
! Non-monotonic / JFNK self-consistency enforcement (Devin,
! 16-Oct-2026).  The trigger flags DO_NONMONO_KNOP / DO_NONMONO_SOB /
! DO_JFNK auto-configure their dependent switches here so a model
! cannot run with an inconsistent combination; every override is
! logged to LUER.  This block sits after the last read of every
! dependent key it touches (GLOBAL_LINE/LAM_SET ~744, FG_OPT ~596,
! USE_LAM_ES ~2347, the JFNK block ~1236-1660, the BA trace keys
! ~1818) and before the validation blocks below, which then re-check
! the enforced values.  INC_GRID (ACCURATE) is read after this
! point, so its KNOP override lives with the INC_GRID read itself.
! NG_ABSTAIN is likewise read after this block: its SOB absent-key
! auto-set stays here, and its base default is set just below so the
! read site needs no re-default (a probe immediately after the read
! would miss the key -- the GET_KEY_STRING cursor brackets the index
! just consumed).
!
! One former opt-in override survives only as an absent-key default:
! under DO_NONMONO_KNOP an unset USE_LAM_ES defaults to T (the
! verified formal-solve pairing); an explicit F is honoured (the
! MOM_J_CMF_NONMONO_V1 moment solve).  Under DO_NONMONO_SOB the
! default is F -- the user-requested moment path -- with an explicit
! T selecting the V13 formal J.  DO_KNOP_REFRESH=F (the lagged-scheme
! A/B comparison) is still overridden under KNOP.
! ------------------------------------------------------------------
!
! NG_ABSTAIN has no module initializer; its read site carries no
! re-default (see below), so the base default must be set here --
! before the SOB auto-set -- on every path.
!
	  NG_ABSTAIN=.FALSE.
!
! Rule 1: DO_NONMONO_KNOP.  The Knop chi_hat formalism exists only in
! the BLANKeted INTEGRAL path of FG_J_CMF_NONMONO_V1, so GLOBAL_LINE
! must be BLANK (or NONE with all TRANS_*=BLANK, still validated
! below), LAM_SET must stay off (it would reset GLOBAL_LINE to NONE
! and hand edge lines to the clamped-sigma Sobolev routines), and
! FG_OPT=DIFF is rejected outright at its own read (no recovery).
! N_LAM_RWARM and GS_REL_STABLE are raised only when left at their
! floor values; higher explicit settings are respected.
!
	  IF(DO_NONMONO_KNOP)THEN
	    IF(GLOBAL_LINE_SWITCH(1:5) .NE. 'BLANK' .AND.
	1       GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE')THEN
	      GLOBAL_LINE_SWITCH='BLANK'
	      WRITE(LUER,'(A)')' KNOP auto-set: GLOBAL_LINE=BLANK'
	    END IF
	    IF(SET_TRANS_TYPE_BY_LAM)THEN
	      SET_TRANS_TYPE_BY_LAM=.FALSE.
	      WRITE(LUER,'(A)')' KNOP auto-set: LAM_SET=F'
	    END IF
	    IF(FG_SOL_OPTIONS(1:3) .NE. 'INT')THEN
	      FG_SOL_OPTIONS='INT/INS'
	      WRITE(LUER,'(A)')' KNOP auto-set: FG_OPT=INT/INS'
	    END IF
	    IF(.NOT. USE_LAM_ES)THEN
! edited by Devin, 21-Oct-2026: an explicit USE_LAM_ES=F is honoured
! under KNOP -- MOM_J_CMF_NONMONO_V1 is a self-consistent J-update in
! this mode (its ray moments come from FG_J_CMF_NONMONO_V1 on
! SIGMA_TRUE).  Only the absent-key default is lifted to T, the
! verified formal-solve pairing; explicit settings are never
! silently rewritten.
	      CALL GET_KEY_STRING(TMP_STRING,'USE_LAM_ES',L_FALSE,
	1     L_FALSE_TMP)
	      IF(.NOT. L_FALSE_TMP)THEN
	        USE_LAM_ES=.TRUE.
	        WRITE(LUER,'(A)')' KNOP auto-set: USE_LAM_ES=T'
	      ELSE
	        WRITE(LUER,'(A)')
	1       ' ============================================================'
	        WRITE(LUER,*)'KNOP field solver: MOM_J_CMF_NONMONO_V1'
	1       ,' moment solve -- explicit USE_LAM_ES=F honoured,'
	2       ,' no override applied.'
	        WRITE(LUER,'(A)')
	1       ' ============================================================'
	      END IF
	    END IF
	    IF(.NOT. DO_KNOP_REFRESH)THEN
	      DO_KNOP_REFRESH=.TRUE.
	      CALL SET_FG_KNOP_REFRESH(DO_KNOP_REFRESH)
! MOM_RF_ENABLED was wired from the pre-auto-set value; re-sync it
! here so the moment solver agrees with the FG side.
	      CALL SET_MOM_RF_DIR_OPT(MOM_RF_DIR,MOM_EPS_LIM,
	1       DO_KNOP_REFRESH,MOM_RF_OWNW)
	      WRITE(LUER,'(A)')' KNOP auto-set: DO_KNOP_REFRESH=T'
	    END IF
	    IF(NO_VEL_FOR_CONTINUUM)THEN
! NOV_CONT=T with GLOBAL_LINE other than BLANK (e.g. NONE) clears
! CONT_VEL in COMP_J_BLANK, which silently hands the continuum to a
! path with no KNOP coupling.  The non-monotonic run needs the
! comoving continuum, so the switch must stay off here.
	      NO_VEL_FOR_CONTINUUM=.FALSE.
	      WRITE(LUER,'(A)')' KNOP auto-set: NOV_CONT=F'
	    END IF
	    IF(N_LAM_RWARM .EQ. 0)THEN
	      N_LAM_RWARM=4
	      WRITE(LUER,'(A)')' KNOP auto-set: N_LAM_RWARM=4'
	    END IF
	    IF(GS_REL_STABLE .LT. 2)THEN
	      GS_REL_STABLE=2
	      WRITE(LUER,'(A)')' KNOP auto-set: GS_REL_STABLE=2'
	    END IF
!
! GS_PRECOND defaults T under KNOP (edited by Devin, 20-Oct-2026):
! the m15_2p82d KNOP run held populations at ~1E10% demand while
! the field cycled release/refreeze.  An explicit VADAT setting is
! respected -- only the absent-key default is lifted.
!
	    CALL GET_KEY_STRING(TMP_STRING,'GS_PRECOND',L_FALSE,
	1     L_FALSE_TMP)
	    IF(.NOT. L_FALSE_TMP)THEN
	      GS_PRECOND=.TRUE.
	      WRITE(LUER,'(A)')' KNOP auto-set: GS_PRECOND=T'
	    END IF
!
! GS_FREEZE_LIM defaults to 300 under KNOP (edited by Devin,
! 20-Oct-2026): a release/refreeze limit cycle resets N_FROZEN at
! every release, so with no in-model cap a stuck run burns the
! whole requeue budget doing forced-LAMBDA iterations that cannot
! converge (m15_2p82d: ~200 KNOP blocks with no progress).
! Explicit VADAT values, including 0 = off, are respected.
!
	    CALL GET_KEY_STRING(TMP_STRING,'GS_FREEZE_LIM',L_FALSE,
	1     L_FALSE_TMP)
	    IF(.NOT. L_FALSE_TMP)THEN
	      GS_FREEZE_LIM=300
	      WRITE(LUER,'(A)')' KNOP auto-set: GS_FREEZE_LIM=300'
	    END IF
!
! GS_PRECOND_LIM defaults to 50 under KNOP (edited by Devin,
! 20-Oct-2026): the precondition gate clears GS_FREEZE every call
! it fires on, so N_FROZEN never advances and GS_FREEZE_LIM is
! unreachable there.  A demand that stays above GS_PRECOND_MAX
! through 50 live updates is a structural failure, not a
! transient.  Explicit VADAT values, including 0 = off, are
! respected.
!
	    CALL GET_KEY_STRING(TMP_STRING,'GS_PRECOND_LIM',L_FALSE,
	1     L_FALSE_TMP)
	    IF(.NOT. L_FALSE_TMP)THEN
	      GS_PRECOND_LIM=50
	      WRITE(LUER,'(A)')' KNOP auto-set: GS_PRECOND_LIM=50'
	    END IF
!
! GS_REL_NOPROG_LIM defaults to 10 under KNOP (edited by Devin,
! 20-Oct-2026): the release/refreeze cycle refunds the freeze
! clock at every release, so the cumulative no-progress count is
! the earliest reachable stop on the m15_2p82d-class alternation
! (~10 bad releases inside ~60 iterations).  Explicit VADAT
! values, including 0 = off, are respected.
!
	    CALL GET_KEY_STRING(TMP_STRING,'GS_REL_NOPROG_LIM',
	1     L_FALSE,L_FALSE_TMP)
	    IF(.NOT. L_FALSE_TMP)THEN
	      GS_REL_NOPROG_LIM=10
	      WRITE(LUER,'(A)')' KNOP auto-set: GS_REL_NOPROG_LIM=10'
	    END IF
!
! DO_NONMONO_FREQ_MATCH defaults T under KNOP: with one-bin upwind
! coupling at compressing points the lagged map gain approaches 1
! when the per-cell Doppler shift exceeds the frequency spacing, so
! the refresh residual asymptotes above KN_RAY_TOL and the field
! never matures (Ceulemans et al. 2025, JQSRT, arXiv:2505.08309;
! Knop et al. 2009 II, A&A 501, 813).  An explicit VADAT setting is
! respected -- only the absent-key default is lifted.
!
	  CALL GET_KEY_STRING(TMP_STRING,'DO_NONMONO_FREQ_MATCH',L_FALSE,
	1     L_FALSE_TMP)
	  IF(DO_NONMONO_KNOP .AND. .NOT. L_FALSE_TMP)THEN
	    DO_NONMONO_FREQ_MATCH=.TRUE.
	    CALL SET_FG_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)
	    CALL SET_MOM_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)
	    WRITE(LUER,'(A)')' KNOP auto-set: DO_NONMONO_FREQ_MATCH=T'
	  END IF
!
! edited by Devin, 17-Oct-2026: the VADAT dump above echoes the
! pre-auto-set values, so it showed GS_REL_STABLE=1 while the solver
! ran with 2 (KNOP_jfnkrel triage read the stale number).  Re-print
! the effective post-auto-set knobs so the log reflects what runs.
!
	    WRITE(LUER,'(A,I4,A,I4,3A)')
	1   ' KNOP effective (post auto-set): N_LAM_RWARM=',N_LAM_RWARM,
	2   '  GS_REL_STABLE=',GS_REL_STABLE,'  GLOBAL_LINE=',
	3   GLOBAL_LINE_SWITCH(1:5),'  FG_OPT=',FG_SOL_OPTIONS(1:7)
!
! W1-W9 effective knobs (17-Oct-2026): echo the maturity-chain
! repair controls so the run log shows what the refresh/gate
! actually ran with.  commit_conv/commit_nc/bulk_contr print the
! effective values, so the post-19-Oct-2026 defaults (T/0.05/T)
! show up here whenever VADAT leaves them unset.
!
	    WRITE(LUER,'(A,L2,A,L2,A,ES10.3,A,I4,A,L2)')
	1   '   omega_v2=',FG_ADAPT_OMEGA_V2_V,
	2   '  plateau_det=',FG_PLATEAU_DET_V,
	3   '  commit_jump=',MAX(FG_COMMIT_JUMP_V,1.0_LDP),
	4   '  heal_cap_exp=',INT(LOG10(MAX(FG_DIV_HEAL_CAP_V,1.0_LDP))),
	5   '  commit_conv=',FG_COMMIT_CONV_V
	    WRITE(LUER,'(A,ES10.3,A,L2,A,L2)')
	1   '   commit_nc=',MAX(FG_COMMIT_NC_V,0.0_LDP),
	2   '  bulk_contr=',FG_BULK_CONTR_V,
	3   '  no_gain=',KN_NO_GAIN
	    WRITE(LUER,'(A,I4,A,I4,A,I4,A,ES10.3)')
	1   '   hot_passes=',KN_RAY_HOTPASSES_V,
	2   '  hot_total=',KN_HOT_TOTAL_PASSES_V,
	3   '  fg_trial_min=',FG_TRIAL_MIN,
	4   '  gs_settle=',GS_SETTLE_TOL
	    WRITE(LUER,'(A,L2,A,L2,A,L2,A,L2)')
	1   '   jfnk_fadapt=',JFNK_FADAPT_V,
	2   '  jfnk_fdcplx=',JFNK_FD_COMPLEX,
	3   '  jfnk_esfr=',JFNK_ENRICH_SFR,
	4   '  jfnk_earm=',JFNK_ENRICH_ARN
	    WRITE(LUER,'(A,L2,A,L2,A,L2)')
	1   '   freq_match=',DO_NONMONO_FREQ_MATCH,
	2   '  knop_refresh=',DO_KNOP_REFRESH,
	3   '  nov_cont=',NO_VEL_FOR_CONTINUUM
!
! MOM-side accelerators (echo only): Anderson on the refresh map
! (MOM_RF_AA, default T), the same-map omega estimator
! (MOM_ADAPT_OMEGA_V2, default T) and the per-cell cycle-damp floor
! MOM_DAMP_MIN.  These govern whether the ~1E-1 lagged-drift floor
! is attacked at all; MOM_DAMP_MIN is read back through a getter
! because it is a solver PARAMETER, not a VADAT key.
!
	    CALL GET_MOM_DAMP_MIN(MOM_DAMP_E)
	    WRITE(LUER,'(A,L2,A,L2,A,ES10.3)')
	1   '   mom_rf_aa=',MOM_RF_AA,
	2   '  mom_om_v2=',MOM_ADAPT_OMEGA_V2_V,
	3   '  mom_damp_min=',MOM_DAMP_E
	  END IF
!
! Rule 2: DO_NONMONO_SOB.  SIGMA_TRUE only reaches the Sobolev line
! routines, so GLOBAL_LINE must be SOB (or the validated all-SOB NONE
! form).  JFNK_COUPLED=T under SOB is rejected at its own read --
! the coupled residual replays the KNOP continuum chain and never
! recomputes the SOB line rates.  USE_LAM_ES is no longer forced
! (21-Oct-2026): the default is the moment path, so the within-sweep
! lag refresh (DO_KNOP_REFRESH, armed by default under LAM_ES=F) is
! live under SOB as well.
!
	  IF(DO_NONMONO_SOB)THEN
	    IF(GLOBAL_LINE_SWITCH(1:3) .NE. 'SOB' .AND.
	1       GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE')THEN
	      WRITE(LUER,'(A)')
	1     ' ============================================================'
	      WRITE(LUER,*)'Warning: DO_NONMONO_SOB overrides an explicit'
	1     ,' GLOBAL_LINE=',GLOBAL_LINE_SWITCH,' with SOB --'
	2     ,' SIGMA_TRUE only reaches the Sobolev line routines,'
	3     ,' so every line must be on the SOB treatment.'
	      WRITE(LUER,'(A)')
	1     ' ============================================================'
	      GLOBAL_LINE_SWITCH='SOB'
	      WRITE(LUER,'(A)')' SOB auto-set: GLOBAL_LINE=SOB'
	    END IF
!
! edited by Devin, 21-Oct-2026: USE_LAM_ES is no longer forced under
! SOB.  The user asked for the moment path as the SOB default
! J-update: MOM_J_CMF_NONMONO_V1 with the fixed-J + GET_MOM_ALO_DIAG
! + GET_MOM_CELL_QUAR linearization in do_var_cont.f -- the same
! recipe as KNOP+LAM_ES=F.  An explicit USE_LAM_ES=T keeps the V13
! formal J update.  Both modes are announced loudly so the solver
! in use is never silent.
!
	    IF(.NOT. USE_LAM_ES)THEN
	      WRITE(LUER,'(A)')
	1     ' ============================================================'
	      WRITE(LUER,*)'SOB field solver: MOM_J_CMF_NONMONO_V1 moment'
	1     ,' solve (USE_LAM_ES=F, the SOB default).  J satisfies'
	2     ,' the sign-split moment equation on SIGMA_TRUE and'
	3     ,' the continuum Jacobian carries the MOM diagonal ALO.'
	4     ,' FG_J_CMF_V13 still runs each frequency to feed'
	5     ,' K_ON_J/N_ON_J ray moments.  Known risk: bound-free'
	6     ,' edge doublets shrink dLOG_NU and inflate the'
	7     ,' advection coefficients.  USE_LAM_ES=T selects the'
	8     ,' V13 formal J update.'
	      WRITE(LUER,'(A)')
	1     ' ============================================================'
!
! DO_NONMONO_FREQ_MATCH defaults T under SOB+LAM_ES=F: the same
! one-bin upwind coupling argument as the KNOP auto-set applies --
! the moment solver's lagged J_NEXT map gains order unity when the
! per-cell Doppler shift exceeds the frequency spacing.  Explicit
! VADAT settings are respected.
!
	      CALL GET_KEY_STRING(TMP_STRING,'DO_NONMONO_FREQ_MATCH',
	1     L_FALSE,L_FALSE_TMP)
	      IF(.NOT. L_FALSE_TMP)THEN
	        DO_NONMONO_FREQ_MATCH=.TRUE.
	        CALL SET_FG_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)
	        CALL SET_MOM_FREQ_MATCH(DO_NONMONO_FREQ_MATCH)
	        WRITE(LUER,'(A)')' SOB auto-set: DO_NONMONO_FREQ_MATCH=T'
	      END IF
	    ELSE
	      WRITE(LUER,'(A)')
	1     ' ============================================================'
	      WRITE(LUER,*)'SOB field solver: FG_J_CMF_V13 formal-solve J'
	1     ,' update (explicit USE_LAM_ES=T; the SOB default is F).'
	2     ,' V13 is built on clamped SIGMA, so it is expansion-'
	3     ,' biased in compression zones; the verified'
	4     ,' non-monotonic formal solver is DO_NONMONO_KNOP.'
	      WRITE(LUER,'(A)')
	1     ' ============================================================'
!
! DO_FG13_ALO defaults T under SOB+LAM_ES=T: FG_J_CMF_V13 stores
! the diagonal OAB operator of the formal solve it just ran, so
! the continuum Jacobian gets real dJ/dCHI, dJ/dETA self-coupling
! instead of a pure Lambda step.  Explicit settings respected.
!
	      CALL GET_KEY_STRING(TMP_STRING,'DO_FG13_ALO',L_FALSE,
	1     L_FALSE_TMP)
	      IF(.NOT. L_FALSE_TMP)THEN
	        DO_FG13_ALO=.TRUE.
	        WRITE(LUER,'(A)')' SOB auto-set: DO_FG13_ALO=T'
	      END IF
	    END IF
!
! Accelerator defaults under SOB (absent-key only, 21-Oct-2026):
! DO_BA_RES_ACC is the rollback net -- the freeze/watchdog stack is
! KNOP-scoped, so a rejected-correction guard is the only safety on
! the unverified moment path; DO_BA_RES_REUSE recalls the stored BA
! matrix while the judged residual still contracts; NG_ABSTAIN keeps
! the contracting components of an NG correction rather than vetoing
! the whole step.  The NG_ABSTAIN probe must stay here: its read
! site sits ~300 lines later, after enough other reads that the
! GET_KEY_STRING cursor skip no longer applies (22-Oct-2026).
!
	    CALL GET_KEY_STRING(TMP_STRING,'DO_BA_RES_ACC',L_FALSE,
	1     L_FALSE_TMP)
	    IF(.NOT. L_FALSE_TMP)THEN
	      DO_BA_RES_ACC=.TRUE.
	      WRITE(LUER,'(A)')' SOB auto-set: DO_BA_RES_ACC=T'
	    END IF
	    IF(DO_BA_RES_ACC)THEN
	      CALL GET_KEY_STRING(TMP_STRING,'DO_BA_RES_REUSE',L_FALSE,
	1       L_FALSE_TMP)
	      IF(.NOT. L_FALSE_TMP)THEN
	        DO_BA_RES_REUSE=.TRUE.
	        WRITE(LUER,'(A)')' SOB auto-set: DO_BA_RES_REUSE=T'
	      END IF
	    END IF
	    CALL GET_KEY_STRING(TMP_STRING,'NG_ABSTAIN',L_FALSE,
	1     L_FALSE_TMP)
	    IF(.NOT. L_FALSE_TMP)THEN
	      NG_ABSTAIN=.TRUE.
	      WRITE(LUER,'(A)')' SOB auto-set: NG_ABSTAIN=T'
	    END IF
	  END IF
!
! Rule 3: DO_JFNK.  Arm the stable/globalized formulation when the
! solver is on: log unknown, pseudo-transient continuation, Armijo
! backtracking and the Cauchy fallback.  JFNK_MAX_SEC=0 leaves no
! wall-clock guard, so a stuck activation could consume the whole
! job; a zero value is raised to the one-hour budget.
!
	  IF(USE_JFNK_SOLVER)THEN
	    IF(JFNK_MAX_SEC .LE. 0.0_LDP)THEN
	      JFNK_MAX_SEC=3600.0_LDP
	      WRITE(LUER,'(A)')' JFNK auto-set: JFNK_MAX_SEC=3600'
	    END IF
	    IF(.NOT. jfnk_log_unknown)THEN
	      jfnk_log_unknown=.TRUE.
	      WRITE(LUER,'(A)')' JFNK auto-set: JFNK_LOGX=T'
	    END IF
	    IF(.NOT. jfnk_ptc_active)THEN
	      jfnk_ptc_active=.TRUE.
	      WRITE(LUER,'(A)')' JFNK auto-set: JFNK_PTC=T'
	    END IF
	    IF(.NOT. JFNK_LINE_SEARCH)THEN
	      JFNK_LINE_SEARCH=.TRUE.
	      WRITE(LUER,'(A)')' JFNK auto-set: JFNK_LS=T'
	    END IF
	    IF(.NOT. JFNK_CAUCHY)THEN
	      JFNK_CAUCHY=.TRUE.
	      WRITE(LUER,'(A)')' JFNK auto-set: JFNK_CAUCHY=T'
	    END IF
	    IF(JFNK_MAXCH_CNV .LE. 0.0_LDP)THEN
	      WRITE(LUER,*)'Warning: JFNK_MXCH_TOL=0 - the secondary'
	1    ,' convergence target is inactive; it needs a ready'
	2    ,' BA diagonal before it can gate admission.'
	    END IF
	  END IF
!
! Rule 4: trace-polish arming.  Only values left at their auto (-1)
! sentinels are set; explicit user settings, including 0 = off, are
! respected.  The demand floors exempt trace rows from the MAXCH
! metric, and BA_TR_MAXSWP bounds the Stage-B polish sweeps once the
! majors converge.  An explicit BA_TR_MAXSWP=0 stays 0 and disables
! the polish (BA_TRACE_NEEDED requires it > 0); the previous 0 -> 5
! rewrite made the documented disable unreachable.  Likewise an
! explicit BA_TR_ABSMAX=0 stays 0 and disables the absolute trace
! gate (both BA_TRACE_NEEDED legs require it > 0); the previous
! 0 -> 1E-08 rewrite made that off-state unreachable.
!
	  IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)THEN
	    IF(BA_DEM_PFLOOR .LT. 0.0_LDP)THEN
	      BA_DEM_PFLOOR=1.0E-30_LDP
	      WRITE(LUER,'(A)')' Trace auto-set: BA_DEM_PFLOOR=1.0E-30'//
	1     ' (-1 auto)'
	    END IF
	    IF(BA_TR_ABSMAX .LT. 0.0_LDP)THEN
	      BA_TR_ABSMAX=1.0E-08_LDP
	      WRITE(LUER,'(A)')' Trace auto-set: BA_TR_ABSMAX=1.0E-08'//
	1     ' (-1 auto)'
	    END IF
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
! The enforcement block above no longer forces USE_LAM_ES under SOB
! (21-Oct-2026): F is the default and takes continuum J from
! MOM_J_CMF_NONMONO_V1 with the fixed-J + GET_MOM_ALO_DIAG
! linearization; an explicit T keeps the V13 formal ray solution.
! The moment path's runtime status on the current build is
! unverified -- the m15 divergence probes predate the solver rework.
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
! on the standard (non-extended) grid, with J from the formal solution
! (Lambda iteration).  The enforcement block above keeps T as the
! absent-key default but honours an explicit USE_LAM_ES=F
! (MOM_J_CMF_NONMONO_V1); the stock MOM_J_CMF_V11 (clamped SIGMA,
! one-way PREV coupling) is never selected in KNOP mode either.
!
	  IF(DO_NONMONO_KNOP)THEN
	    IF(USE_LAM_ES)THEN
	      WRITE(LUSCR,'(A)')' DO_NONMONO_KNOP=T: Knop et al. 2009'//
	1       ' formal solver (FG_J_CMF_NONMONO_V1, USE_LAM_ES=T)'
	    ELSE
	      WRITE(LUSCR,'(A)')' DO_NONMONO_KNOP=T: Knop non-monotonic'//
	1       ' moment solver (MOM_J_CMF_NONMONO_V1, USE_LAM_ES=F)'
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
! residual gate to roll back, and the GS release trials are adjudicated
! through BAACC_PENDING -- unjudged without it.  Absent-key auto-set to
! T (22-Oct-2026): every validated KNOP config already sets it
! explicitly, and SOLVEBA_V13 handles the LAMBDA-phase calls and the
! MOM/FG lagged-field rollback under it.  An explicit F is still
! honoured for controlled F baselines (the warning below keeps
! covering that case).  The GET_KEY_STRING probe is safe here -- the
! key's read site sits ~1200 lines earlier.
!
	    CALL GET_KEY_STRING(TMP_STRING,'DO_BA_RES_ACC',L_FALSE,
	1     L_FALSE_TMP)
	    IF(.NOT. L_FALSE_TMP)THEN
	      DO_BA_RES_ACC=.TRUE.
	      WRITE(LUER,'(A)')' KNOP auto-set: DO_BA_RES_ACC=T'
	    END IF
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
! transfer and is not non-monotonic safe.  Auto-set rather than STOP
! (16-Oct-2026 self-consistency pass; INC_GRID is read after the main
! enforcement block above, so its override lives here).
!
	  IF(DO_NONMONO_KNOP .AND. ACCURATE)THEN
	    ACCURATE=.FALSE.
	    WRITE(LUER,'(A)')' KNOP auto-set: INC_GRID=F - the'//
	1   ' enhanced-grid branch calls MOM_J_CMF_V11 after the'//
	2   ' formal transfer and is not non-monotonic safe.'
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
! key leaves the upstream all-or-nothing veto in place.  The base
! default is set before the enforcement block above (NG_ABSTAIN
! carries no module initializer), so an SOB auto-set survives this
! read site -- the old unconditional .FALSE. here silently undid it
! (announced T, ran F; 22-Oct-2026).  NB: GET_KEY_STRING cannot
! re-probe this key here -- its I_DWN/I_UP cursor brackets the index
! just consumed by RD_STORE_LOG, so a probe immediately after the
! read reports not-found.  Presence checks must run earlier.
!
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
