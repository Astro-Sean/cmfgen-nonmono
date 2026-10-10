	MODULE CONTROL_VARIABLE_MOD
	USE SET_KIND_MODULE
!
! Altered : 13-Jan-2024 : Added MAX_NO_GREY_ITERATIONS
! Altered : 12-Aug-2022 : Added SN shock variables (following work by LUC.)
! Altered : 06-Jun-2022 : Added COMP_STEQ_T_EHB
! Altered : 14-Jul-2019 : R_PNT_SRCE, TEFF_PNT_SRCE, NC_PNT_SRCE, LOGICAL PNT_SRCE_MOD, USE_ELEC_HEAT_BAL (from osiris)
!                            Added to ibis 14-Aug-2019
! Altered : 17-Oct-2016 : H_CHK_OPTION replaces CHECK_H_ON_J.
! Altered : 01-Sep-2016 : TIME_SEQ_NO changed from integer to real.
! Altered : 15-Jun-2016 : NT_ITERATION_COUNTER added to this routine so that it is accessible by CMFGEN_SUB.
! Altered : 15-Feb-2015 : Added INSTANTANEOUS_ENERGY_DEPOSITION option (12-Jan-2015 on OSPREY[cur_cmf_gam])
! Incorporated: 02-Jan-2014: Changed to allow depth dependent profiles.
! Altered : 05-Apr-2011 : Added vriable R_GRD_REVISED (10-Feb-2011).
! Altered : 16-Jul-2010 : Added FIX_ALL_SPECIES variable.
! Altered : 23-Nov-2007 : LAM_SCALE_OPT included.
! Altered : 20-Feb-2006 : ABOVE_EDGE changed to LOGICAL from REAL(KIND=LDP)
! Altered : 29-Jan-2006 : Control variable fors relativistic transfer and time
!                          dependent statistical equilibrium equations installed.
!
! Set constants that are regularly used, and passed to subroutines.
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ITHREE=3
	INTEGER, PARAMETER :: IFOUR=4
	INTEGER, PARAMETER :: IFIVE=5
	INTEGER, PARAMETER :: ISIX=6
	INTEGER, PARAMETER :: ISEV=7
	INTEGER, PARAMETER :: ITEN=10
!
	REAL(KIND=LDP), PARAMETER :: RZERO=0.0_LDP
	REAL(KIND=LDP), PARAMETER :: RONE=1.0_LDP
	REAL(KIND=LDP), PARAMETER :: RTWO=2.0_LDP
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! Used to tack INPUT/OUPUT erros.
!
	INTEGER IOS
!
!*******************************************************************************
!
! Star and velocity specifications
!
	REAL(KIND=LDP) RP			!Stellar radius in units of 10^10 cm
	REAL(KIND=LDP) RMAX
	REAL(KIND=LDP) RMDOT
      	REAL(KIND=LDP) VRP,RN,VINF,EPPS1,GAMMA1
	REAL(KIND=LDP) LUM
	REAL(KIND=LDP) TEFF
	REAL(KIND=LDP) LOGG
	REAL(KIND=LDP) PRESSURE_VTURB
	REAL(KIND=LDP) RP2,VRP2,RN2,VINF2,EPPS2,GAMMA2
	REAL(KIND=LDP) CONS_FOR_R_GRID
	REAL(KIND=LDP) EXP_FOR_R_GRID
!
	INTEGER N_IB_INS			!Number of points for fine grid at inner boundary
	INTEGER N_OB_INS			!Number of points for fine grid at outer boundary
	INTEGER, PARAMETER :: NUM_V_OPTS=1
	CHARACTER*10 VEL_OPTION(NUM_V_OPTS)
!
! Parameters for VELTYPE=3 or VELTYPE=6 (STARPCYG_V2)
!
	REAL(KIND=LDP) SCL_HT,VCORE,VPHOT
	REAL(KIND=LDP) VINF1,V_BETA1,V_EPPS1
	REAL(KIND=LDP) V_BETA2,V_EPPS2		!VINF2 above
	INTEGER VELTYPE
!
! Variable Mdot parameters
!
	LOGICAL VAR_MDOT
	CHARACTER*80 VAR_MDOT_FILE
!
! Supernova variables
!
	REAL(KIND=LDP) RHO_ZERO			!Density at core in gm/cm^3
	REAL(KIND=LDP) RCUB_RHO_ZERO            !r^3 . density at core in (in 10^{-30} gm}
	REAL(KIND=LDP) N_RHO			!Exponent for density (+ve)
	REAL(KIND=LDP) DJDT_RELAX_PARAM         !Factor to assist inclusion of DJDT terms.
	REAL(KIND=LDP) SN_AGE_DAYS              !Age of SN in days.
	REAL(KIND=LDP) RMAX_ON_RCORE
	REAL(KIND=LDP) TIME_SEQ_NO              !Number of model in time sequence (need not be an integer).
!
	LOGICAL SN_MODEL
	LOGICAL SN_HYDRO_MODEL		!Use HYDRO model for SN input
	LOGICAL PURE_HUBBLE_FLOW        !Forces a pure Hubble flow
	LOGICAL INCL_RADIOACTIVE_DECAY	!Allow for radioactive decay.
	LOGICAL JGREY_WITH_V_TERMS 	!Include velocity terms when comuting T(grey).
	LOGICAL INCL_DJDT_TERMS         !Include DJDt terms in radiative transfer equation
	LOGICAL USE_DJDT_RTE            !Use the radiative transfer equation solver which includes DJDt terms
	LOGICAL USE_DR4JDT              !Explicitly difference Dr4JDt (rather than Dr3JDt)
	LOGICAL DO_CO_MOV_DDT		!Include comoving drivative in SE equations.
	LOGICAL DO_FULL_REL_OBS         !Include all relativistic terms in obs. frame computation.
	LOGICAL DO_FULL_REL_CMF         !Include all relativistic terms in CMF obs. frame computation.
!
! LUC: parameters to treat the SHOCK POWER in interacting supernova models
!
	LOGICAL INC_SHOCK_POWER        !Include shock power
	LOGICAL ADD_SHOCK_POWER_SLOWLY !Introduce the shock power slowly starting with scaling SHOCK_POWER_FAC_BEG
	LOGICAL SCL_PWR_BY_FCL         !Do we scale the shock_power by CLUMP_FAC?
	REAL(KIND=LDP)  SHOCK_POWER_FAC_BEG    !Initial scaling factor read from VADAT file
	REAL(KIND=LDP)  SHOCK_POWER_FAC        !Current scaling factor
	REAL(KIND=LDP)  PRESCRIBED_SHOCK_POWER !Value in erg/s of the desired power
	REAL(KIND=LDP)  VLOC_SHOCK_POWER       !Velocity at which we dump the shock power
	REAL(KIND=LDP)  DVLOC_SHOCK_POWER      !Gaussian width for the shock power deposition profile
!
! Non thermal parameters
!
	LOGICAL TREAT_NON_THERMAL_ELECTRONS     !
	LOGICAL ADD_DEC_NRG_SLOWLY
	LOGICAL READ_NON_THERM_SPEC
	REAL(KIND=LDP) DEC_NRG_SCL_FAC_BEG
	REAL(KIND=LDP) DEC_NRG_SCL_FAC
	REAL(KIND=LDP) NT_OMIT_ION_SCALE        !Ion omitted if < NT_OMIT_ION_SCL*(largest ion pop.).
	REAL(KIND=LDP) NT_OMIT_LEV_SCALE        !Fractional pops omitted if < NT_OMIT_SCL
        REAL(KIND=LDP) NT_EMAX			!Maximum energy of non-thermal electrons
        REAL(KIND=LDP) NT_EMIN			!Minimum energy of non-thermal electrons
	INTEGER NT_NKT
	INTEGER NON_THERMAL_IT_CNTRL    !Controls how often we update the nonthermal electron distribution.
	INTEGER NT_ITERATION_COUNTER
	LOGICAL COMP_GREY_LST_IT        !Comput J(GREY) on last iteration [DEFAULT is TRUE].
	LOGICAL SCL_NT_CROSEC
	LOGICAL SCL_NT_ION_CROSEC
	CHARACTER(LEN=12) NT_SOURCE_TYPE
!
	REAL(KIND=LDP) MINIMUM_ISO_POP          !Minimum isotop population for included isotope.
!
	LOGICAL USE_J_REL
	LOGICAL USE_FORMAL_REL
	LOGICAL USE_LAM_ES
	LOGICAL INCL_REL_TERMS
	LOGICAL INCL_ADVEC_TERMS_IN_TRANS_EQ
	LOGICAL USE_OLD_MF_SCALING
	LOGICAL USE_OLD_MF_OUTPUT
!
! Used to revise the R_GRID automatically after each iteration.
!
	INTEGER, PARAMETER :: N_RG_PAR_MAX=4
	INTEGER N_RG_PAR
	REAL(KIND=LDP) RG_PAR(N_RG_PAR_MAX)
	LOGICAL REVISE_R_GRID
	LOGICAL R_GRID_REVISED
        LOGICAL INSTANTANEOUS_ENERGY_DEPOSITION
	CHARACTER(LEN=10) NEW_RGRID_TYPE
	CHARACTER(LEN=10) SN_T_OPTION
	CHARACTER(LEN=10) GAMRAY_TRANS
!
! Variables for including clumping under the assumption that the clumping
! occurs on scales much smaller than any of the transfer scales. Valid
! for the continuum --- likely to be invalid for lines due to the
! shortness of the SObolev length.
!
! CLUMP_FAC represents the fractional volume occupied by material. Thus its
! value is always less than unity. The remaining volume is assumed to be a
! vacuum. In practice we are assuming dense spherical shells separated by
! a vacuum.
!
	INTEGER, PARAMETER ::  N_CLUMP_PAR_MAX=6
	REAL(KIND=LDP) CLUMP_PAR(N_CLUMP_PAR_MAX)
	INTEGER N_CLUMP_PAR
	LOGICAL DO_CLUMP_MODEL
	CHARACTER(LEN=6) CLUMP_LAW
!
!************************************************************************
!
! IF GRID is TRUE, variables are interpolated directly on old grid.
!
	LOGICAL GRID
	LOGICAL INTERP_DC_SPH_TAU
	LOGICAL SET_LTE_AS_INIT_ESTIMATES
	LOGICAL DO_HYDRO
	LOGICAL DONE_HYDRO_REVISION
        LOGICAL AUTO_ADD_ION_STAGES
	LOGICAL AUTO_SMOOTH_POPS
	CHARACTER(LEN=10) DC_INTERP_METHOD
!
! Used when constructing the Temperature distribution on the first
! iteration.
!
	REAL(KIND=LDP) GREY_PAR
	REAL(KIND=LDP) T_INIT_TAU
	INTEGER  MAX_NO_GREY_ITERATIONS
	LOGICAL ITERATE_INIT_T
	LOGICAL T_MIN_BA_EXTRAP
	LOGICAL INTERP_T_ON_R_GRID
!
! Used to limit the temperature while iterating on T.
!
	REAL(KIND=LDP) T_MIN
!
! When you use the LTE option to set the departure coefficients, T_EXCITE_MIN
! is used in the outer region to set the escitation coefficient.

	REAL(KIND=LDP) T_EXCITE_MIN
!
! Indicates that the R grid for the NEWMODEL should be rad in from a file --
! ist is not computed.
!
	LOGICAL RDINR
!
! Reads in a modified from of STEQ_VALS so at to apply corrections to the populations.
! For debuggung purposes only.
!
	LOGICAL RDINSOL
!
! Controls for dynamic smoothing of the photoioization cross-sections, and for
! treating dielectronic lines.
!
	REAL(KIND=LDP) SIG_GAU_KMS
	REAL(KIND=LDP) FRAC_SIG_GAU
	REAL(KIND=LDP) CUT_ACCURACY
	REAL(KIND=LDP) VSM_DIE_KMS
	LOGICAL DIE_AS_LINE
	LOGICAL ABOVE_EDGE
!
! Variables for controlling the frequency grid.
!
	REAL(KIND=LDP) MIN_CONT_FREQ 		!Minimum continuum frequency.
	REAL(KIND=LDP) MAX_CONT_FREQ    	!Maximum continuum frequency.
	REAL(KIND=LDP) SMALL_FREQ_RAT 		!Fractional spacing for small frequencies'
	REAL(KIND=LDP) dFREQ_bf_MAX		!Maximum spacing close to bf edge.
	REAL(KIND=LDP) BIG_FREQ_AMP		!Amplification factor
	REAL(KIND=LDP) dV_LEV_DIS		!dV on low side of bound-free edge.
	REAL(KIND=LDP) AMP_DIS			!Amplification factor
	REAL(KIND=LDP) MIN_FREQ_LEV_DIS		!Minimum frequency for lev dissolution.
	LOGICAL RD_CONT_FREQ		!Read in cont. frequencies from file.
	INTEGER FREQ_GRID_OPTION	!Which frequency grid
!
	REAL(KIND=LDP) R_PNT_SRCE
	REAL(KIND=LDP) TEFF_PNT_SRCE
	INTEGER NC_PNT_SRCE
	LOGICAL PNT_SRCE_MOD
!
! Used to describe the Doppler lines profile, which is currently assumed
! to be constant at all depths.
!
	REAL(KIND=LDP) TDOP
	REAL(KIND=LDP) VDOP_FIX
	REAL(KIND=LDP) VTURB
	REAL(KIND=LDP) AMASS_DOP
	REAL(KIND=LDP) VTURB_MIN,VTURB_MAX,VTURB_POW
	REAL(KIND=LDP) VTURB_VEND
	LOGICAL FIX_DOP
	CHARACTER(LEN=10) VTURB_LAW
!
	REAL(KIND=LDP) DOP_PROF_LIMIT
	REAL(KIND=LDP) VOIGT_PROF_LIMIT
	REAL(KIND=LDP) V_PROF_LIMIT
	REAL(KIND=LDP) MAX_PROF_ED
	LOGICAL SET_PROF_LIMS_BY_OPACITY
	LOGICAL RD_STARK_FILE
	LOGICAL NORM_PROFILE
	CHARACTER(LEN=10) GLOBAL_LINE_PROF
!
! Variables used to define what transitions are neglected. Cut is
! presently by the gf value, and the lower level of the transition.
!
	REAL(KIND=LDP) GF_CUT
	REAL(KIND=LDP) AT_NO_GF_CUT
	INTEGER GF_LEV_CUT
	INTEGER MIN_NUM_TRANS
!
! Variables for treating lines simultaneously with the continuum.
!
	REAL(KIND=LDP) V_DOP
	REAL(KIND=LDP) MAX_DOP
	REAL(KIND=LDP) FRAC_DOP
	REAL(KIND=LDP) dV_CMF_PROF
	REAL(KIND=LDP) dV_CMF_WING
	REAL(KIND=LDP) ES_WING_EXT
	REAL(KIND=LDP) R_CMF_WING_EXT
	REAL(KIND=LDP) EXT_LINE_VAR
	REAL(KIND=LDP) ZNET_VAR_LIMIT
!
! Variables for controlling the frequency grid in the Observer's frame.
!
	REAL(KIND=LDP) OBS_PRO_EXT_RAT
	REAL(KIND=LDP) dV_OBS_PROF
	REAL(KIND=LDP) dV_OBS_WING
	REAL(KIND=LDP) dV_OBS_BIG
!
! Compute all photoionization cross-sections at all frequencies. If not,
! we assume thei value from an earlier frequncy. Saves considerable computational
! effort. DELV_CONT gives the maximum velocity separation in km/s between
! evaluations of the phot. cross-sections.

	REAL(KIND=LDP) DELV_CONT
	LOGICAL COMPUTE_ALL_CROSS
	LOGICAL COMPUTE_NEW_CROSS
	LOGICAL USE_FIXED_J
!
!***************************************************************************************
!
! Method to handle N moment  (N_ON_J, MIXED, or G_ONLY)
!
	CHARACTER(LEN=6)  N_TYPE
	CHARACTER(LEN=10) J_CHK_OPTION
	CHARACTER(LEN=10) H_CHK_OPTION
	CHARACTER(LEN=10) XM_CHK_OPTION
!
! These insert extra points into the grid when solving for the radition field.
! These insertions are done internally, and do not directly effect the returned
! results. USed mainly with CMF_FLUX.
!
	REAL(KIND=LDP) DELV_FRAC_FG
	REAL(KIND=LDP) DELV_FRAC_MOM
!
! Indicate solution method used to ocmpute the Eddington factors by FG_J_CMF_V?.
! Options is usually  INT/INS (or DIFF/INS).
!
	CHARACTER(LEN=10) FG_SOL_OPTIONS
!
! Indicates type of BA matrix (DIAG, TRIDIAG etc).
!
	CHARACTER(LEN=6)  METH_SOL
	CHARACTER(LEN=10) CMF_FORM_OPTIONS		!Used for formal solution.
	CHARACTER(LEN=10) NEG_OPAC_OPTION
	CHARACTER(LEN=12) TWO_PHOTON_METHOD
!
! Specifies method used to compute optical depth.
!
	CHARACTER(LEN=6)  METHOD
	CHARACTER(LEN=6)  LUM_FROM_ETA_METHOD
!
        REAL(KIND=LDP) OVER_FREQ_DIF
        REAL(KIND=LDP) WEAK_LINE_LIMIT
	REAL(KIND=LDP) WEAK_TAU_LINE_LIMIT
!
	LOGICAL WEAK_WITH_NET
        LOGICAL USE_WEAK_TAU_LIM
	LOGICAL OVERLAP
	LOGICAL NEW_VAR_STORAGE_METHOD
	INTEGER NUM_OF_WEAK_LINES    	!Counter in CMFGEN_SUB only
!
! Indicates whether how ALL lines are to be treated.
!
	CHARACTER*6 GLOBAL_LINE_SWITCH
	LOGICAL SET_TRANS_TYPE_BY_LAM
	LOGICAL DO_SOBOLEV_LINES
!
	REAL(KIND=LDP) FLUX_CAL_LAM_BEG
	REAL(KIND=LDP) FLUX_CAL_LAM_END
!
! FLUX_CAL_ONLY provides a method for computing the continuous spectrum
! only (i.e. no linearization or population corrections):
!    To get a BLANKETED spectrum FLUX_CAL_ONLY should be set to
!       TRUE and GLOBAL_LINE_SWITCH to BLANK
!    To get a pure UNBLANKETED spectrum FLUX_CAL_ONLY should be set to
!       TRUE and GLOBAL_LINE_SWITCH to SOB
!
	LOGICAL FLUX_CAL_ONLY
!
! Indicates whether lines treated in SOB and CMF mode are allowed for when
! constructing the observers frame grid.
!
	LOGICAL SOB_FREQ_IN_OBS
!
! Variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	LOGICAL EDD_CONT
	LOGICAL COMPUTE_EDDFAC
	LOGICAL EDD_LINECONT
	LOGICAL NO_VEL_FOR_CONTINUUM
!
! Indicate whether the atmosphere is to be treated in plane-parallel mode
!
	LOGICAL PLANE_PARALLEL_NO_V
	LOGICAL PLANE_PARALLEL
!
! Variables required for adding additional points into the depth
! grid. Allows an increase in program accuracy to overcome
! rapid ionization changes.
!
	REAL(KIND=LDP) ACC_FREQ_END
	INTEGER NPINS			!Points inserted for error calc.
        INTEGER ST_INTERP_INDX        !Interp from ST_INT.. to END_INTERP..
        INTEGER END_INTERP_INDX
        LOGICAL ACCURATE
	LOGICAL ALL_FREQ
	CHARACTER*10 INTERP_TYPE
!
! These two variables used only in CMFGEN_SUB.
!
	LOGICAL THIS_FREQ_EXT
        LOGICAL INACCURATE
!
! ND-DEEP to DEEP we use a quadratic interpolation scheme so as to try
! and preserve "FLUX" in the diffusion approximation.
!
	INTEGER DEEP
	REAL(KIND=LDP) ACC_EDD_FAC
!
! SUPPRESS_FCONV_IT1 (dev7, diagnostic-volume control).  When TRUE, the
! per-frequency 'Possible error converging f - T1' warnings in
! COMP_J_BLANK are counted instead of written during great iteration 1,
! and CMFGEN_SUB prints a single count + max-residual line at the end of
! the iteration.  Iteration 1 builds EDDFACTOR from FOLD=0 with
! unrelaxed populations, so the f fixed-point can exceed its pass cap at
! thousands of frequencies and flood OUTGEN.  Warnings on later
! iterations (mid-run EDDFACTOR rebuilds) are unaffected.  VADAT keyword
! [SUP_FCONV_IT1]; default TRUE, set F in VADAT to restore upstream
! per-frequency printing on iteration 1.
!
	LOGICAL SUPPRESS_FCONV_IT1
	LOGICAL :: FCONV_IT1_ACTIVE=.FALSE.
	INTEGER :: FCONV_SUPPRESSED_CNT=0
	REAL(KIND=LDP) :: FCONV_SUPPRESSED_MAX=0.0_LDP
!
	LOGICAL EXTEND_FRM_SOL
	LOGICAL INSERT_FREQ_FRM_SOL
!
	REAL(KIND=LDP) ES_VAR_FAC
	LOGICAL MIXED_ES_VAR
!
!***********************************************************************************
! Variables for performing NG acceleration.
!
	INTEGER LAST_NG   		!Indicates iteration on which last
!                               	    NG acceleration occurred.
	INTEGER NEXT_NG   		!Indicates iteration on which next
!                                          NG acceleration is to occur.
	INTEGER IT_TO_BEG_NG       	!Iteration to beg NG acceleration
	INTEGER ITS_PER_NG  		!Iterations between NG accelerations.
	INTEGER NG_BAND_WIDTH      	!Number of depths to acclerate simultaneously
	REAL(KIND=LDP) VAL_DO_NG 		!Begin NG when MAXCH < VAL_DO_NG
	LOGICAL NG_DONE			!iIndicates successfull completion of NG
	LOGICAL NG_DO			!Switch on NG acceleration.
	LOGICAL DO_NG_VALIDITY_CHECK
!
! Variables for performing averaging of oscilating populations.
!
	INTEGER LAST_LAMBDA
	INTEGER LAST_AV
	INTEGER NEXT_AV
	INTEGER ITS_PER_AV
	INTEGER NUM_OSC_AV
	LOGICAL AVERAGE_DONE
	LOGICAL AVERAGE_DO
	LOGICAL UNDO_LAST_IT
!
! N_PAR is used to indicate how often the BA matrices should be incremented
! by BA_PAR. After the incrementation, the PAR matrices are zeroed.
!
	INTEGER N_PAR
	REAL(KIND=LDP) MAX_LAM_COR	!Maximum fractional change for Lambda iteration.
	REAL(KIND=LDP) MAX_LIN_COR	!Maximum fractional change for linearization.
	REAL(KIND=LDP) MAX_dT_COR       !Maximum allowed change in the temperature.
	REAL(KIND=LDP) MAX_CHNG_LIM
!
! Indicates how to scale the corections to the populations, so that
! populations don't go negative, and MAX_LIN_COR (or MAX_LAM_COR) is
! satisfied. Options are NONE, LOCAL, MAJOR, and GLOBAL.
! MAJOR is generally the preferred choice.
!
	CHARACTER(LEN=6)  SCALE_OPT
!
! Introduced Nov-2003 to limit LAMBDA changes prior to normal scaling of corrections.
! Option only has an effect if LAM_SCALE_OPT='LIMIT'
!
	CHARACTER(LEN=6)  LAM_SCALE_OPT
!
! Used to determine the accuracy with which the BA and BAION matrices are
! computed. 1.0D-10 < BA_CHK_FAC < 0.1. Smaller number means higher acuracy,
! but slower computation. Too large a value may affect convergence. Most
! models run with 1.0D-04. Basically updates to BA ignored if their effect
! on RJ is less than BA_CHK_FAC*RJ.
!
	REAL(KIND=LDP) BA_CHK_FAC
	LOGICAL INCLUDE_dSLdT
	LOGICAL NEW_LINE_BA
	INTEGER INDX_BA_METH_RD
!
! Variables for scaling the line cooling rates in oder that the radiative
! equilibrium equation is more consistent with the electron heating/cooling
! equation. The scaling is done when the line frequency is with a fraction
! of SCL_LINE_HT_FAC of the tmean frequency for the super-level under
! consideration. 0.5 is presently the prefered value.
!
! SCL_SL_LINE_OPAC effects the opacities to achieve consistency.
! SCL_LINE_COOL_RATES only modifies the line cooling rates to achieve consistency.
! Only one of these two options can be true at any one time.
!
	LOGICAL SCL_SL_LINE_OPAC
	LOGICAL SCL_LINE_COOL_RATES
	REAL(KIND=LDP) SCL_LINE_HT_FAC
	REAL(KIND=LDP) SCL_LINE_DENSITY_LIMIT
	REAL(KIND=LDP) EPS                    !If maximum fractional % change < EPS we terminate model
!
! Fix BA variation matrix if % change less than VAL_FIX_BA.
!
	REAL(KIND=LDP) VAL_FIX_BA
	INTEGER N_ITS_TO_FIX_BA
	INTEGER CNT_FIX_BA
	LOGICAL COMPUTE_BARDIN,COMPUTE_BA
	LOGICAL WRBAMAT,WRBAMAT_RDIN
!
! Indicates whether inverse of BA matrix should be output. This takes up
! disk space, but saves computaional effort, especially wen BA matrix is
! being held fixed.

	LOGICAL WR_BA_INV
	LOGICAL WR_PART_OF_INV
!
! Performs a lambda iteration if % change > VAL_DO_LAM
!
	REAL(KIND=LDP) VAL_DO_LAM
! MAXCH_CTRL is the current iteration's committed-change metric. Under
! DO_NONMONO_KNOP it is applied; otherwise it equals MAXCH, so the
! default path is unchanged. CMFGEN_SUB uses it for 50-percent gates.
	REAL(KIND=LDP) :: MAXCH_CTRL=0.0_LDP
	INTEGER CNT_LAM
	INTEGER RD_CNT_LAM
!
! VAL_ENTER_LAM / VAL_EXIT_LAM / N_STABLE_LAM (VADAT [LAM_ENTER],
! [LAM_EXIT], [LAM_N_STABLE]): hysteresis for the LAMBDA<->full mode
! switch.  A full step enters LAMBDA only when MAXCH_CTRL exceeds
! VAL_ENTER_LAM; a LAMBDA sequence exits to full only after
! MAXCH_CTRL has stayed below VAL_EXIT_LAM for N_STABLE_LAM
! consecutive (unforced) iterations.  CNT_STABLE_LAM counts the
! qualifying iterations.  Defaults reproduce the old single
! threshold: both thresholds = LAM_VAL, N_STABLE_LAM=1.
!
	REAL(KIND=LDP) VAL_ENTER_LAM
	REAL(KIND=LDP) VAL_EXIT_LAM
	INTEGER N_STABLE_LAM
	INTEGER CNT_STABLE_LAM
!
! FRESH_BA_EXIT (VADAT [FRESH_BA_EXIT], default F): when a LAMBDA
! sequence exits to a full step, force COMPUTE_BA and WRBAMAT for
! the next iteration so the full step runs against a Jacobian built
! from the current (post-LAMBDA) state instead of whatever matrix or
! inverse survived from before the LAMBDA window.
!
	LOGICAL FRESH_BA_EXIT
!
! BA_ENTER_RAW (VADAT [BA_ENTER_RAW], optional, default 1.0D30 =
! disabled) and BA_LAST_RAW: raw-demand gate on the LAMBDA->full
! merit exit.  BA_LAST_RAW records the raw (pre-limitation) MAXCH of
! the most recent non-LAMBDA SOLVEBA call.  On the non-monotonic
! paths (DO_NONMONO_SOB, DO_NONMONO_KNOP) the LAMBDA iterate can sit
! below LAM_EXIT while the BA
! direction is still non-contractive -- every admitted full step
! then demands a saturated ~1E5-1E7% decrease at the stiff inner
! cells and rebounds the SE residual several-fold (the
! lam_stable_exit/maxch_above_enter limit cycle).  When the gate is
! enabled a merit exit is held while the last observed BA demand is
! still >= BA_ENTER_RAW; the NUM_LAM escape still fires, so each
! escape doubles as a probe that refreshes BA_LAST_RAW and re-opens
! the gate once the Newton direction contracts.  Default leaves the
! schedule unchanged on every path.
!
! BA_PROBE_LAM (VADAT [BA_PROBE_LAM], optional, default 0 = off):
! while the raw-demand gate is holding, release one ordinary full
! step (bounded by BA_TRUST/resacc like any exit) once the hold has
! persisted for BA_PROBE_LAM LAMBDA calls.  BA_LAST_RAW only
! refreshes on a non-LAMBDA call, so without a probe the gate can
! never learn that the BA direction has begun to contract -- the
! NUM_LAM escape is the only re-open channel.  A failed probe
! demotes and resets CNT_LAM, so probing cannot spam; the probe
! cadence shortens NUM_LAM, it never extends it.
!
	REAL(KIND=LDP) BA_ENTER_RAW
	REAL(KIND=LDP) BA_LAST_RAW
	INTEGER BA_PROBE_LAM
!
! BA_DAMP_CAP (VADAT [BA_DAMP_CAP], optional, default 0 = off): when
! positive, a non-LAMBDA step whose raw correction demand saturates
! (MAXCH >= BA_SAT_RAW) is kept in full mode instead of demoted, and
! the next non-LAMBDA call's per-variable correction factor is
! clamped to MIN(MAX_LIN_COR,BA_DAMP_CAP) while BA_LAST_RAW stays
! saturated.  A capped step still moves the model along the BA
! direction; the damped cap bounds each excursion to ~(BA_DAMP_CAP-1)
! while the demand is re-measured every call.  Consecutive saturated
! steps are bounded by BA_SAT_LIM -- clamped to 8 in
! RD_CONTROL_VARIABLES when BA_DAMP_CAP is armed but BA_SAT_LIM is
! unset -- then ba_sat_demote fires.  Residual acceptance stays
! active on every step.  Optional so VADAT files without the key are
! unaffected.
!
	REAL(KIND=LDP) BA_DAMP_CAP
!
! Phased temperature release (VADAT [N_FULL_TREL] default 0=off,
! [T_REL_CHG] default 3.0): after a LAMBDA->full exit, keep FIXED_T
! set regardless of RD_FIX_T until N_ACC_FULL accepted non-LAMBDA,
! non-rejected, converged iterations have accumulated AND MAXCH_CTRL falls
! below T_REL_CHG, then restore RD_FIX_T.  T_HELD_PHASE marks the
! hold as active; both are reset on (re)entry to LAMBDA.  The hold
! lets the population Jacobian re-baseline after a LAMBDA window
! before temperature corrections resume.
!
	INTEGER N_FULL_TREL
	REAL(KIND=LDP) T_REL_CHG
	INTEGER N_ACC_FULL
	LOGICAL T_HELD_PHASE
	LOGICAL OLD_RD_LAMBDA
	LOGICAL RD_LAMBDA
	LOGICAL DO_LAMBDA_AUTO
	LOGICAL DO_GREY_T_AUTO
	LOGICAL DO_T_AUTO
	LOGICAL LAMBDA_ITERATION
	LOGICAL STOP_IF_BAD_PARAM
	LOGICAL STOP_IF_MAJOR_WARNING
!
! Variables for determining whether some populations are held fixed
! when the new populations are solved for.
!
	REAL(KIND=LDP) CON_SCL_T
	REAL(KIND=LDP) TAU_SCL_T
	LOGICAL RD_FIX_T,RD_FIX_IMP
	LOGICAL FIXED_T,FIXED_NE,FIX_IMPURITY
	LOGICAL VARFIXT
	LOGICAL DO_SRCE_VAR_ONLY
	LOGICAL FIX_ALL_SPECIES
	LOGICAL SET_POPS_D2_EQ_D1
	LOGICAL USE_ELEC_HEAT_BAL
	LOGICAL COMP_STEQ_T_EHB
!
	REAL(KIND=LDP) ADD_OPAC_SCL_FAC
 	LOGICAL ADD_ADDITIONAL_OPACITY
!
! Variable for specifying abundance set;
!
	CHARACTER(LEN=10) SOL_ABUND_REF_SET
!
! 
! Indicates whether level dissolution should be included in the model.
!
	LOGICAL DO_LEV_DISSOLUTION
!
	LOGICAL DIF		  	! Use diffusion approximation?
!
	LOGICAL INCL_ADIABATIC   	!Include adiabatic cooling
	LOGICAL INCL_CHG_EXCH		!Include charge exchange reactions.
	LOGICAL INCL_TWO_PHOT		!Include two-photon transitions
	LOGICAL INCL_RAY_SCAT           !Include Rayleigh scattering.
	LOGICAL INCL_PENNING_ION        !Include H/He Pening ionization
	LOGICAL LINEAR_ADV              !Comput derivatives using linear approximation.
!
	LOGICAL INCL_ADVECTION          !Include advection terms in rate equations.
	REAL(KIND=LDP)  ADVEC_RELAX_PARAM       !Allows advection terms to be added slowly
!
! Variables specifying computaion of e.s. souce function.
!
	LOGICAL COHERENT_ES
	LOGICAL RD_COHERENT_ES
	LOGICAL USE_OLDJ_FOR_ES
!
	REAL(KIND=LDP) TSTAR
	REAL(KIND=LDP) IC
!
! X-ray variables. XRAYS indicates whether or not to include X-ray emission
! from shocks in the wind? If FF_XRAYS is true, we assume the X-rays are
! generated by free-free proceses. If XRAY_SMOOTH_WIND is set, we factor out the
! effect of cluming when computing the effect of the X-ray emissivity.
!
	LOGICAL XRAYS
	LOGICAL ADD_XRAYS_SLOWLY
	LOGICAL FF_XRAYS
	LOGICAL XRAY_SMOOTH_WIND
!
! The shocks are characterized by up to 2 temperatures.
!
	REAL(KIND=LDP) FILL_FAC_XRAYS_1,T_SHOCK_1,V_SHOCK_1
	REAL(KIND=LDP) FILL_FAC_XRAYS_2,T_SHOCK_2,V_SHOCK_2
	REAL(KIND=LDP) FILL_FAC_X1_BEG,FILL_X1_SAV
	REAL(KIND=LDP) FILL_FAC_X2_BEG,FILL_X2_SAV
	REAL(KIND=LDP) SLOW_XRAY_SCL_FAC
	REAL(KIND=LDP) XRAY_EMISS_1,XRAY_EMISS_2
	REAL(KIND=LDP) VSMOOTH_XRAYS
!
	REAL(KIND=LDP) ALLOWED_XRAY_FLUX_ERROR       !Fractional error allowed before X-ray emissivties are scaled
        REAL(KIND=LDP) DESIRED_XRAY_LUM              !Desired X-ray luminosity (> 0.1 keV)
	LOGICAL SCALE_XRAY_LUM               !Indicates whether the X-ray emissivities will be scaled automatically.
!
	REAL(KIND=LDP) NU_XRAY_END
	REAL(KIND=LDP) DELV_XRAY
!
! Indicates number of time POPS array is to be written to scratch
! file per iteration.
!
        INTEGER RITE_N_TIMES
        PARAMETER (RITE_N_TIMES=1)
!
	INTEGER NLBEGIN
        LOGICAL NEWMOD
        LOGICAL WRITE_RVSIG
!
	CHARACTER(LEN=10) OUTER_BND_METH
	CHARACTER(LEN=10) INNER_BND_METH
!
	REAL(KIND=LDP) IB_STAB_FACTOR
	REAL(KIND=LDP) OUT_BC_PARAM_ONE
	REAL(KIND=LDP) REXT_FAC
	INTEGER RD_OUT_BC_TYPE
	INTEGER OUT_BC_TYPE
	LOGICAL RDTHK_CONT
	LOGICAL THK_CONT
	LOGICAL THK_LINE
	LOGICAL INCL_INCID_RAD
!
	LOGICAL LTE_MODEL
	LOGICAL SETZERO
	lOGICAL DO_POP_SCALE
	LOGICAL TRAPFORJ
	LOGICAL CHECK_LINE_OPAC
        LOGICAL SOBOLEV
        LOGICAL VERBOSE_OUTPUT
	LOGICAL WRITE_RATES
	LOGICAL WRITE_JH
!
! Non-monotonic velocity-field option (dev7).  When TRUE, the code uses
! the Dessart, Audit & Hillier (2015, arXiv:1503.05463) approach for
! non-monotonic velocity fields: standard CMF formal solver with a
! Sobolev line treatment.  VADAT keyword [DO_NONMONO_SOB]; default FALSE
! preserves the original code path.  Requires GLOBAL_LINE=SOB (or NONE
! with all TRANS_*=SOB) and USE_LAM_ES=T; RD_CONTROL_VARIABLES aborts
! otherwise (03-Oct-2026).
!
	LOGICAL DO_NONMONO_SOB
!
! Non-local Sobolev coupling (dev7).  When TRUE, and DO_NONMONO_SOB is
! also TRUE, the escape-probability line treatment adds the Rybicki &
! Hummer (1978, ApJ 219, 654) coupling between distant Sobolev points:
! photons emitted at one resonance are re-absorbed at a remote depth
! where the projected velocity along the ray matches.  Implemented in
! SOB_RH_COUPLE (subs/sob_rh_couple.f) as a correction to ZNET_SIM.
! VADAT keyword [DO_NONLOC_SOB]; default FALSE preserves the local
! escape-probability path.
!
	LOGICAL DO_NONLOC_SOB
!
! Non-monotonic velocity-field option (dev7).  When TRUE, COMP_J_BLANK
! calls FG_J_CMF_NONMONO_V1, the Knop et al. (2009) formal solver for
! non-monotonic velocity fields, in place of FG_J_CMF_V13.  The solver
! keeps |GAM| in the generalized opacity and the signed coupling in a
! lagged source term, so compression zones (sigma < -1) transfer
! correctly.  Both J-update modes are supported: USE_LAM_ES=T takes J
! from the formal solution (Lambda iteration); USE_LAM_ES=F calls
! MOM_J_CMF_NONMONO_V1, the non-monotonic VEF moment solver with the
! same chi_hat formalism.  The stock MOM_J_CMF_V11 is never selected:
! its clamped SIGMA is inconsistent with the non-monotonic ray/moment
! solution.  Mutually exclusive with DO_NONMONO_SOB.  Requires
! GLOBAL_LINE=BLANK (or NONE with all TRANS_*=BLANK), LAM_SET=F, and
! INC_GRID=F; RD_CONTROL_VARIABLES aborts otherwise (03-Oct-2026).
! VADAT keyword [DO_NONMONO_KNOP]; default FALSE preserves the
! original code path.
!
	LOGICAL DO_NONMONO_KNOP
!
! N_LAM_RWARM: forced LAMBDA iterations at process start when
! DO_NONMONO_KNOP=T.  The Gauss-Seidel I_NEXT arrays in
! FG_J_CMF_NONMONO_V1 are process-local and are re-zeroed on every
! requeue; this forces LAMBDA for the first N_LAM_RWARM great
! iterations of each process so a BA step is never taken against a
! stale (zeroed) I_NEXT.  Skipped when a mature FGINEXT_GS checkpoint
! restored the field.  VADAT keyword [N_LAM_RWARM]; default 4, 0
! disables.
!
	INTEGER N_LAM_RWARM
!
! MOM_TRIAL_MIN (default 30): first frozen iteration at which the
! marginal-release path may attempt a resacc-gated BA step.  A MOM
! lagged field that is frozen only on the marginal drift terms (no
! capped writes, no runaway cell, >50% drift confined to <1% of
! cells, negJ/quar fractions bounded) is mature enough to try a BA
! step: DO_BA_RES_ACC arms the residual acceptance on that call and
! rolls the step back if the SE residual grows, so the trial cannot
! silently damage the model.  N_FROZEN keeps counting through the
! trial, so KNOP_STALL_MAX still bounds a run that never matures.
! Set <=0 to disable marginal release.  Without DO_BA_RES_ACC the
! trial never fires -- an unchecked BA step on a marginal field is
! exactly the failure the gate exists to prevent.
! VADAT keyword [MOM_TRIAL_MIN].
!
	INTEGER MOM_TRIAL_MIN
!
! GS_FREEZE_LIM (default 0 = disabled): consecutive great iterations
! the KNOP forced-LAMBDA freeze may persist before SOLVE_FOR_POPS
! stops the model as non-convergent.  The freeze itself stays
! uncapped -- a field that cannot mature must not take BA steps --
! but a run that never matures should fail loudly rather than burn
! requeues and report as a completed model.  Independent of this
! key, SOLVE_FOR_POPS publishes the freeze state to the GS_STATE
! file each call so the job wrapper can apply its own limit.
! VADAT keyword [GS_FREEZE_LIM].
!
	INTEGER GS_FREEZE_LIM
!
! GS_TRIAL_LIM (default 15): accepted resacc-gated trial releases the
! marginal MOM path may take between raw-mature sweeps.  An accepted
! trial refunds the freeze clock -- N_FROZEN resets to MOM_TRIAL_MIN
! -- because the field demonstrated the marginal state and the resacc
! net committed the step: real progress the GS_FREEZE_LIM /
! KNOP_STALL watchdogs should not count against the run.  Rejected
! releases do not refund, so a churning run still hits the stall
! limit.  At GS_TRIAL_LIM accepted trials the marginal-release path
! is vetoed; the freeze resumes and the stall limit applies as
! before.  Set <=0 to disable the refund entirely (N_FROZEN then
! counts every frozen iteration, restoring the original behaviour).
! VADAT keyword [GS_TRIAL_LIM].
!
	INTEGER GS_TRIAL_LIM
!
! GS_REL_STABLE (default 1): consecutive sweeps that must pass the
! KNOP maturity gate before SOLVE_FOR_POPS releases the population
! update.  A single mature-looking sweep can be a coincidence of the
! churning field -- a release on one marginal sweep lets the update
! re-poison the field and the freeze/release cycle repeats without
! converging.  Requiring a stable streak confirms the lagged field
! is a fixed point before the populations move.  VADAT keyword
! [GS_REL_STABLE].
!
	INTEGER GS_REL_STABLE
!
! GS_REL_N (default 1): number of consecutive LAMBDA iterations
! after a GS-freeze release that keep the near-mature correction
! cap.  The first released step carries the deferred correction
! accumulated while pops were held; the following steps still see a
! field re-lagging to the moved populations, so they stay damped
! too.  Only LAMBDA steps are damped -- BA steps are covered by the
! residual-acceptance net.  VADAT keyword [GS_REL_N].
!
	INTEGER GS_REL_N
!
! DO_KNOP_REFRESH: per-ray Knop refresh (Oct-2026).  When TRUE and
! DO_NONMONO_KNOP=T, CMFGEN_SUB calls FG_KNOP_REFRESH_V1 after each
! frequency sweep to converge the compressing-branch wavelength
! coupling internally (alternating-direction per-ray re-solve on the
! sweep-captured opacities), replacing the one-frequency-bin-per-
! great-iteration lag that prevented GS maturity.  VADAT keyword
! [DO_KNOP_REFRESH]; default DO_NONMONO_KNOP (set FALSE to recover
! the lagged behaviour for A/B tests).
!
	LOGICAL DO_KNOP_REFRESH
!
! MOM_RF_GS_DESC (default FALSE): on descending MOM_KNOP_REFRESH
! passes, compressing depths read the live J_NEXT_FRQ(:,L+1) --
! updated this pass -- instead of the frozen pass-start snapshot.
! True Gauss-Seidel propagation along the red-to-blue compressing
! characteristic; the snapshot caps that propagation at one bin per
! pass.  VADAT keyword [MOM_RF_GS_DESC].
!
! MOM_RF_KEEPNEG (default TRUE): cells whose MOM solution hit a
! post-solve repair (significant negJ, non-finite or runaway
! rejection) keep the previous J_NEXT_FRQ iterate instead of storing
! the repaired value, so the lagged read is not reseeded with repair
! noise.  VADAT keyword [MOM_RF_KEEPNEG].
!
	LOGICAL MOM_RF_GS_DESC
	LOGICAL MOM_RF_KEEPNEG
!
! MOM_GATE_BULK (default TRUE): judge MOM lagged-field maturity by
! the fraction of cells still moving (nbad5/nbad50 vs NTOT, capped
! writes, max change vs GS_UNSAFE_TOL, negJ fraction <10%) instead
! of the max-norm + any-negJ test that permanently freezes BA on
! strongly compressed grids.  Mirrors the FG-side bulk criterion.
! VADAT keyword [MOM_GATE_BULK].
!
	LOGICAL MOM_GATE_BULK
!
! PRINT_KNOP_DEBUG (default FALSE): verbose Knop-solver diagnostics
! in OUTGEN -- currently the O(NP*NI) I_NEXT depth-mapping table.
! The mapping check itself always runs once per grid build and
! reports an ERROR on anomalies; the flag only controls the table
! print.
!
	LOGICAL PRINT_KNOP_DEBUG
!
! DO_NONMONO_FREQ_MATCH (default FALSE): at compressing points the
! lagged upwind frequency is rematched to the bin the photon
! actually occupies one cell upwind, n = |GAM|*ds/dln(nu) bins,
! instead of always the adjacent bin.  The frequency derivative
! then spans n*dln(nu), which is what the discretization needs
! when the per-cell Doppler shift exceeds one grid bin (the
! oscillatory-eigenvalue regime of Ceulemans et al. 2025, JQSRT).
! Implemented in FG_J_CMF_NONMONO_V1 and MOM_J_CMF_NONMONO_V1;
! expanding cells and the monotonic path are unchanged.  VADAT
! keyword [DO_NONMONO_FREQ_MATCH].
!
	LOGICAL DO_NONMONO_FREQ_MATCH
!
! KN_RAY_ORDERED (default TRUE): ordered (i,l) march inside
! FG_KNOP_REFRESH_V1 -- rows are marched along the characteristic
! and the per-row wavelength direction follows sign(GAM), the
! physically ordered Gauss-Seidel of Knop et al. (2009, Paper II
! Sec. 2).  Only the lagged slope corner and the downwind half of
! the compressing interpolant stay lagged, replacing the
! alternating-direction passes that leave half of each strong
! coupling stale every pass.  F restores them for A/B tests.
! Requires the refresh coefficient pack and NI>=3 per ray; the
! alternating scheme is the automatic fallback otherwise.  VADAT
! keyword [KN_RAY_ORDERED].
!
	LOGICAL KN_RAY_ORDERED
!
! KN_ORD_AITKEN (default FALSE): guarded Aitken delta-squared
! extrapolation on the ordered-march (KN_RAY_ORDERED=T) field
! iterate at pass boundaries inside SWEEP_FLD_ORD.  Targets the
! marginally contractive |Q|->1 compression cells that stall the
! march short of the bulk-maturity exit.  VADAT [KN_ORD_AITKEN].
!
	LOGICAL KN_ORD_AITKEN
!
! MOM_RF_AA (default TRUE): vector Anderson acceleration on the
! MOM refresh composite map, using the MOM_AA_M most recent
! iterate/residual pairs.  For the affine-in-lagged-field refresh
! map AA(m) is equivalent to the GMRES residual minimization on
! (I-A)x=b, so a separate GMRES pass is not currently warranted.
! VADAT keyword [MOM_RF_AA].
!
	LOGICAL MOM_RF_AA
!
! MOM_AA_OUT (default TRUE): call-level Anderson acceleration on the
! committed lagged J_NEXT_FRQ field.  The inner AA accelerates the
! frozen-coefficient refresh map; this one accelerates the outer
! sweep+refresh composition, whose iterates still carry the
! J<->Eddington coefficient drift.  Pairs cost 2*ND_NM*NCF*8 bytes
! each and MOM_AA_OUT_M of them are kept.  VADAT keywords
! [MOM_AA_OUT] and [MOM_AA_OUT_M].
!
	LOGICAL MOM_AA_OUT
!
! MOM_RF_BULK_EXIT (default TRUE): exit the MOM refresh once the
! per-pass bulk drift meets the same bounds the maturity gate
! applies to the sweep (<=0.1% of cells >50%, no capped writes,
! negJ/quarantined <=10%), instead of iterating to the max-norm
! stall.  The converged flag then routes the gate to the
! retained-vs-sweep delta, so a fixed point far from the field the
! populations saw is still rejected.  VADAT keyword
! [MOM_RF_BULK_EXIT].
!
	LOGICAL MOM_RF_BULK_EXIT
!
! MOM_RF_ACTIVE (default FALSE): active-set front-tracking over
! frequencies in the MOM refresh -- skip frequencies whose own
! residual and all upwind read slots stayed below KN_M_TOL on the
! previous pass.  VADAT keyword [MOM_RF_ACTIVE].
!
	LOGICAL MOM_RF_ACTIVE
!
! MOM_EPS_LIM (default TRUE): donor-cell positivity limiter on the
! lagged EPS midpoint transfer in the Knop moment solver -- the
! lagged auxiliary term may not debit a row past its accumulated
! positive RHS.  VADAT keyword [MOM_EPS_LIM].
!
! MOM_RF_DIR (default 'AUTO'): pass ordering for MOM_KNOP_REFRESH.
! AUTO measures the compressing-row fraction on the sweep and runs
! only the dominant coupling direction once one side holds >=90% of
! rows (the physically-guided ordering of Knop Paper II Sect. 3.2);
! ALT keeps the original alternating order, ASC/DESC force a
! direction.  VADAT keyword [MOM_RF_DIR].
!
! MOM_RF_OWNW (default TRUE): pass-direction ownership of
! lagged-store writes inside MOM_KNOP_REFRESH -- a row whose upwind
! read was stale this pass (compressing rows on an ascending pass,
! expanding on a descending one) does not write back, so each row
! only takes fresh-read Gauss-Seidel steps.  Removes the residual
! floor of the alternating GS*Jacobi map.  VADAT keyword
! [MOM_RF_OWNW].
!
	LOGICAL MOM_EPS_LIM
	LOGICAL MOM_RF_OWNW
	CHARACTER(LEN=4) MOM_RF_DIR
!
! DO_BA_RES_ACC (default FALSE): residual-accepted BA update.
! SOLVEBA_V13 measures the max statistical-equilibrium residual
! (the STEQ right-hand side at entry) into BA_SE_RES.  After the
! solve, a
! BA step is committed only when the residual measured at the
! post-step state does not exceed BA_RES_GROW times the residual
! at the state where the step was taken; otherwise the step is
! rolled back to the stored snapshot, the correction is
! discarded, and the BA matrix is recomputed next iteration.
! A second bound, BA_RES_HARD times the smallest residual seen
! since the tracker armed (BAACC_RES_MIN), closes the slow-growth
! loophole: a run of accepted 1.4x steps would otherwise diverge
! without ever tripping the per-step ratio.  Independently, any
! depth whose requested correction exceeds the current trust
! radius is zeroed, which blocks the stale-inverse
! 1e12-correction pathology.  The radius adapts like a textbook
! trust region: BA_TRUST_CUR starts at BA_TRUST, halves on each
! rejection, and recovers by 1.25x on each acceptance.
! VADAT keywords [DO_BA_RES_ACC], [BA_RES_GROW], [BA_RES_HARD],
! [BA_TRUST].
!
	LOGICAL DO_BA_RES_ACC
	REAL(KIND=LDP) BA_RES_GROW
	REAL(KIND=LDP) BA_RES_HARD
	REAL(KIND=LDP) BA_TRUST
!
! BA_NREJ_LIM (default 20): consecutive rejected BA steps at which
! SOLVEBA_V13 stops the model.  BA_TRUST_CUR floors at 1.0, so
! without a cap a model whose linearization cannot produce an
! acceptable step alternates LAMBDA and rejected-BA corrections
! forever and completes un-converged.  VADAT keyword
! [BA_NREJ_LIM]; <=0 disables.
!
	INTEGER BA_NREJ_LIM
!
! BA_RES_SAT (default 1E30): saturation level for the max-residual
! channel.  When the armed step's pre-step residual is already at or
! above this value, the channel is measuring a degenerate equation
! (the ~1E262 inner-boundary rows on the Dessart runs), not step
! quality -- a bounded correction cannot move it and any overflow
! above it trips the growth test spuriously.  In that regime only
! the RES2 rank channel judges the step.  VADAT keyword
! [BA_RES_SAT]; <=0 disables the exemption (legacy behaviour).
!
	REAL(KIND=LDP) BA_RES_SAT
!
! Internal state for the residual-accepted BA update (not VADAT
! keys).  The variables are maintained by SOLVEBA_V13 and the
! rejection flag is consumed by SOLVE_FOR_POPS.
!
	REAL(KIND=LDP) BA_SE_RES
! BA_SE_RES2: rank-based second residual channel.  It is the K-th
! largest per-depth column max of |STEQ|, so a persistent residual
! floor spanning fewer than K depth columns (the degenerate
! inner-boundary equations, ~1E262 on the Dessart run) cannot mask
! residual growth in the rest of the grid -- the single-metric
! acceptance test was blind to that (FreeNe LAM blow-up, Sep-2026).
	REAL(KIND=LDP) BA_SE_RES2
	LOGICAL BAACC_PENDING
	LOGICAL BA_STEP_REJECTED
! BA_SNAPSHOT_MISSING: set by SOLVEBA_V13 when the POPS_BA_SAVE arm
! allocation fails (or is otherwise unavailable) while DO_BA_RES_ACC
! is armed.  A BA step taken without a snapshot cannot be rolled
! back, so SOLVE_FOR_POPS demotes to a bounded LAMBDA update until
! the flag clears -- an unchecked BA step on a marginal lagged field
! is the failure the resacc net exists to prevent (LAM_dbg it39-48:
! alloc failed silently, unchecked steps wrecked the field).
	LOGICAL BA_SNAPSHOT_MISSING
	REAL(KIND=LDP) BAACC_RES_PREV
	REAL(KIND=LDP) BAACC_RES_MIN
	REAL(KIND=LDP) BAACC_RES2_PREV
	REAL(KIND=LDP) BAACC_RES2_MIN
	REAL(KIND=LDP) BA_TRUST_CUR
	REAL(KIND=LDP) BA_REJ_MAXCH
	INTEGER BA_NREJ
	INTEGER BA_NREJ_CONSEC
!
! BA_REJ_FRAC (VADAT [BA_REJ_FRAC], default 0.25): fraction of depth
! columns allowed to fail the per-depth residual test before an
! armed BA step is rolled back wholesale.  The observed failure is
! localized -- a handful of stiff cells at compression depths drive
! the rebound while the rest of the step contracted (m15: H2:lev1
! @d70, CV:1s_2p_1Po@d13) -- and damping alone does not fix it:
! the rebound is set by the direction, not the step size (SOB v2
! it509, a +50% capped step still rebounded the residual 120x).
! With 0 < frac <= 1 a step that fails at a minority of depths is
! partially accepted: the failing columns are restored from
! POPS_BA_SAVE, their fresh corrections are masked to zero, and the
! rest of the step stands.  The RES2 rank channel (K-th largest
! per-depth column residual) remains the bulk veto: a step that
! degrades the grid broadly still rolls back whole, as does one
! failing more than frac*ND columns.  FRAC<=0 restores the legacy
! whole-state decision on any failed channel.
!
	REAL(KIND=LDP) BA_REJ_FRAC
!
! Per-depth acceptance state (maintained by SOLVEBA_V13, not VADAT).
! BAACC_COL_PREV(I)/BAACC_COL_MIN(I) carry each depth column's
! pre-step and best-seen residual -- COL_RES(I) is already formed
! each call for BA_SE_RES/BA_SE_RES2, so the masked test reuses it.
! BA_RBK_MASK marks columns rolled back this call; SOLVEBA zeroes
! their fresh STEQ before scaling (those corrections were linearized
! at the rejected populations and cannot be trusted).  BA_PART_THIS
! flags a masked rollback to SOLVE_FOR_POPS: like BA_STEP_REJECTED
! it forces a BA recompute, but it does not demote to LAMBDA -- the
! accepted columns made progress and demoting is what recreates the
! FULL/LAMBDA alternation the mask exists to break.  BA_NREJ_SOFT
! counts armed steps whose residual grew but whose follow-up call
! was LAMBDA: no rollback there (undoing after a Lambda step
! reproduces the it33 cycle), but the trust radius still contracts
! and the event is logged rather than counted as acceptance.
!
	REAL(KIND=LDP), ALLOCATABLE :: BAACC_COL_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: BAACC_COL_MIN(:)
	LOGICAL, ALLOCATABLE :: BA_RBK_MASK(:)
	LOGICAL BA_PART_THIS
	INTEGER BA_NREJ_SOFT
!
! BA_SOFT_THIS flags a LAMBDA-follow-up soft rejection on this call:
! the armed step was judged degraded but not undone, so its record
! is a committed-but-bad iterate in the SCRTEMP history.  Consumed
! (and cleared) by SOLVE_FOR_POPS like BA_STEP_REJECTED -- its only
! current use is the NG window-cleanliness counter, since a
! divergent-but-committed iterate corrupts the secant fit exactly
! like a rolled-back one.
!
	LOGICAL BA_SOFT_THIS
!
! DO_BA_QUAR (default FALSE, VADAT [DO_BA_QUAR]): per-equation
! quarantine inside the residual-accepted BA step.  Masked
! acceptance rolls back a whole depth column when any one of its
! equations fails the residual test, and the trust mask zeroes a
! whole column when any variable's |STEQ| demand exceeds
! BA_TRUST_CUR; on the m15 2.82d SOB run a handful of stiff
! inner-ejecta equations (CV:1s_2s_1Se/1s_2p_1Po at d68-69, the NeV
! ground-term and HII ionization equations at d70-75) demand
! >=99.999% annihilation every armed step, which parked the trust
! radius at its floor, zeroed ~60 of 76 columns per armed solve,
! and vetoed every full step.  With DO_BA_QUAR the offending
! (equation,depth) pair is parked instead: its residual is excluded
! from COL_RES/BA_SE_RES/BA_SE_RES2 so it can no longer veto
! acceptance, and its correction is zeroed on non-LAMBDA solves
! before the trust test so it can no longer trip the whole column.
! LAMBDA calls still update it, so a parked variable keeps
! contracting under fixed-point iteration.  The set is populated by
! measured failure only: every over-trust population variable in a
! trust-masked column, and the worst residual equation of each
! rollback-masked column.  ED (J=NT-1) and T (J=NT) are never
! parked.  BA_QUAR_MAX (VADAT, default 25) bounds the set -- a
! parked equation is invisible to the acceptance test, so the blind
! spot stays bounded; <=0 disables growth.  Parked demand is logged
! in BA13-TRACE so the suppression is not silent, and it stays out
! of MAXCH: a parked equation is declared out of play for both
! acceptance and the demand metric, and its LAMBDA-call correction
! still counts toward MAXCH on Lambda iterations so no false
! converged exit is manufactured.  Inert unless DO_BA_RES_ACC; the
! default path is unchanged.
!
	LOGICAL DO_BA_QUAR
	INTEGER BA_QUAR_MAX
	LOGICAL, ALLOCATABLE :: BA_QUAR_MASK(:,:)
	INTEGER BA_QUAR_COUNT
! BA_RAW_TRIAL_PENDING tracks the residual check for the last full BA
! demand. With [DO_BA_RES_ACC] and either non-monotonic flag enabled,
! BA_RAW_REJECTED prevents an unsuccessful demand from reopening the
! LAMBDA or JFNK gates. The default path is unchanged.
	LOGICAL BA_RAW_TRIAL_PENDING
	LOGICAL BA_RAW_REJECTED
!
! BA_RDLAM_PROBE: armed by SOLVE_FOR_POPS when an [BA_PROBE_LAM]-cadence
! bounded trial should run inside the RD_LAMBDA phase.  Without it the
! DO_LAM_AUTO stable-exit gate has no escape when the LAMBDA attractor
! floors above LAM_EXIT -- the model parks in LAMBDA to NUM_ITS (m15
! SOB run: demand floored ~45% and rebounded while stable=0/3).
! CMFGEN_SUB reads the flag to skip the RD_LAMBDA -> LAMBDA_ITERATION
! re-imposition for exactly one call; SOLVE_FOR_POPS consumes it at
! entry.  BA_RDLAM_ACCEPT is armed when that probe survives residual
! acceptance on the next call; it pins CNT_STABLE_LAM at N_STABLE_LAM
! so the ordinary DO_LAM_AUTO switch performs the handoff (keeping the
! IN_ITS DO_LAM_IT write and the FIXED_T/FIX_IMPURITY restore). A clean
! probe is accepted on its Lambda follow-up; a degraded probe instead
! restores the full population/radiation snapshot and skips that call's
! stale Lambda solve. Ordinary BA/Lambda cycles retain soft accounting.
! These flags are inert unless DO_BA_RES_ACC and BA_PROBE_LAM/BA_DAMP_CAP
! are armed on a non-monotonic path; the default path is unchanged.
!
	LOGICAL BA_RDLAM_PROBE
	LOGICAL BA_RDLAM_ACCEPT
	LOGICAL :: BAACC_RDLAM_PROBE_ARM=.FALSE.
	LOGICAL :: BAACC_RDLAM_PROBE=.FALSE.
!
! NG_ABSTAIN (default FALSE, VADAT [NG_ABSTAIN]): when T,
! NG_MIT_BAND_OPT_V1 holds (rather than vetoes) any component whose
! newest stored difference is zero -- a resacc-restored or otherwise
! held column carries no contraction information -- and any
! component the secant would extrapolate non-positive.  Upstream the
! whole depth is vetoed for one bad component (>3 bad depths cancels
! the acceleration); on this model family a single oscillating level
! (OSIX et al., outer shell) vetoed ~50 of 70 depths every attempt.
! Abstention converts those vetoes into partial accelerations.
!
	LOGICAL NG_ABSTAIN
!
! BA saturation demotion (nonmono paths, disabled by default):
! counts consecutive non-LAMBDA calls whose raw correction demand
! stays at or above BA_SAT_RAW.  An applied correction pinned on
! the limitation cap while the raw demand never relaxes is
! marching, not converging -- LAM_dbg it36-47 committed 12 straight
! full steps at the 1E7 sentinel with SE_RES2 flat, until the
! lagged field crossed its basin edge.  Residual rejection can
! only fire after the cliff, so SOLVE_FOR_POPS demotes to LAMBDA
! and halves BA_TRUST_CUR when the count reaches BA_SAT_LIM.
! Extended KNOP -> SOB on 4-Oct-2026: the count also bounds
! [BA_DAMP_CAP] damped-full runs under SOB.  VADAT keywords
! [BA_SAT_LIM], [BA_SAT_RAW]; BA_SAT_LIM<=0 disables.
!
	INTEGER BA_SAT_LIM
	REAL(KIND=LDP) BA_SAT_RAW
	INTEGER BA_SAT_CNT
!
! BA_LUM_GATE and BA_LUM_DCHG are retained as accepted VADAT keys for
! restart compatibility.  They no longer affect convergence: shocks
! permit physical luminosity deposition between the inner and outer
! boundaries.  SOLVE_FOR_POPS reports the ratio and drift as diagnostic
! values only.
!
	LOGICAL BA_LUM_GATE
	REAL(KIND=LDP) BA_LUM_DCHG
!
! POPS_BA_SAVE: pre-step population snapshot taken at the start of
! each non-LAMBDA SOLVEBA call, so a residual-rejected BA step can
! be rolled back on the next call.  Allocated (NT,ND) on first use.
!
	REAL(KIND=LDP), ALLOCATABLE :: POPS_BA_SAVE(:,:)
!
! Convergence diagnostics + phased iteration control (dev7,
! 30-Sep-2026 session).  Everything here is either off by default
! or identical to the pre-existing path when left at its default.
!
! CONV_DIAG_LVL (VADAT [CONV_DIAG_LVL], default 2): verbosity of the
! per-iteration convergence diagnostics.
!   0 = ITER_SUM / ITER_HISTORY only
!   1 = + ITER_MODE transitions, CONV_PLATEAU, GS-freeze dumps, the
!       BA_RES_ACC trajectory line
!   2 = + MODCTL, MODE-OUT, BA-CTX, BA-SOLVE, BA13-TRACE (the output
!       volume produced before this control existed)
!
	INTEGER CONV_DIAG_LVL
!
! LUM_RATIO_OI: outer/inner radiative luminosity ratio
! RLUMST(1)/RLUMST(ND), set in CMFGEN_SUB just before the
! SOLVE_FOR_POPS call (RLUMST lives in RADIATION_MOD).
! Diagnostic only.
!
	REAL(KIND=LDP) LUM_RATIO_OI
!
! GS_HOLD_POPS (VADAT [GS_HOLD_POPS], default F): on a GS-freeze
! iteration, snapshot POPS before the forced-LAMBDA SOLVEBA update
! and restore them afterwards, so the sweep converges the lagged
! coupling field at fixed populations (Knop Paper II nested
! iteration ordering -- the field is a fixed point for the current
! populations, not a moving target).  GS_POPS_HELD records that the
! current call rolled its update back; GS_ATTEMPT_APP keeps the
! attempted committed change for the diagnostics; POPS_GS_SAVE is
! the pre-solve snapshot, allocated (NT,ND) on first use.  The T and
! ED rows are held as well -- the freeze is meant to be a
! transfer-only sweep.
!
	LOGICAL GS_HOLD_POPS
	LOGICAL GS_POPS_HELD
	REAL(KIND=LDP) GS_ATTEMPT_APP
	REAL(KIND=LDP), ALLOCATABLE :: POPS_GS_SAVE(:,:)
!
! ------------------------------------------------------------------
! JFNK (Jacobian-Free Newton-Krylov) population solver.
!
! USE_JFNK_SOLVER (VADAT [DO_JFNK], default F; [USE_JFNK] accepted
! as a deprecated alias): route the population correction solve
! inside SOLVE_FOR_POPS to the JFNK Newton-Krylov solver
! (JFNK_SOLVE_CORR_V1) instead of the direct SOLVEBA_V13
! factorization.  Requires DO_NONMONO_KNOP=T or DO_NONMONO_SOB=T;
! default preserves the original path bit-for-bit.
!
! JFNK_COUPLED (VADAT [JFNK_COUPLED], default F): the frontier
! option -- evaluate the statistical-equilibrium residual with a
! FRESH formal solution at the perturbed populations (each residual
! evaluation re-runs the opacity -> COMP_J_BLANK -> radiative-rate
! chain via JFNK_SWEEP_CB -> JFNK_RAD_SWEEP_V1), so the
! finite-difference JVP probes the full dJ/dn coupling.  This is
! the Arramy et al. (2024, A&A 690, A12) formulation: Newton on the
! coupled SE+RTE system.  KNOP only -- under DO_NONMONO_SOB the
! sweep callback does not recompute SOB line rates, so coupled mode
! is rejected at VADAT read and SOB admits only the frozen tier.
! When F, the residual uses the stored (fixed-J) radiative rates --
! the same fixed-J system SOLVEBA solves exactly, run through the
! globalized Newton machinery.
!
! The remaining keys mirror the V4/V6 JFNK implementation:
! solver tolerances, GMRES restart, trust-region (hookstep)
! controls, correction cap, preconditioner selection.
! ------------------------------------------------------------------
!
	LOGICAL USE_JFNK_SOLVER
	LOGICAL JFNK_COUPLED
! Internal runtime state (no VADAT key): true when [USE_JFNK]=T (default
! F) commits a best iterate, rebuilt baseline, or clipped Newton update.
! Convergence and phased-T gates must not count it as solved.
	LOGICAL :: JFNK_UNCONVERGED_COMMIT=.FALSE.
! JFNK_ANY_ITER (VADAT [JFNK_ANYITER], default F): debug/bring-up
! only -- lets the JFNK correction solve fire on a LAMBDA-designated
! call, so the machinery (residual/JVP/preconditioner path) can be
! exercised without waiting for the mode controller to release a full
! iteration.  The applied update still passes through the normal
! commit/hold/resacc logic, so a Lambda-hold iteration discards it.
	LOGICAL JFNK_ANY_ITER
! JFNK_IN_EVAL: internal runtime flag (no VADAT key).  Set only while
! JFNK_RAD_SWEEP_V1 re-runs the continuum sweep for a trial residual
! evaluation.  COMP_J_BLANK checks it to suppress its WRITE(LU_EDD)
! updates of the lagged-J store in EDDFACTOR -- evals must read the
! outer-iteration records but never overwrite them, else perturbed
! trial J pollutes the J initial guess for the next real sweep.
	LOGICAL JFNK_IN_EVAL
! JFNK_EVAL_BROKEN: internal runtime flag (no VADAT key), companion of
! JFNK_IN_EVAL.  COMP_J_BLANK sets it in place of a STOP when a
! data-integrity or blowup check fails inside an eval sweep -- a STOP
! there would kill the whole run mid-JVP instead of flagging the eval.
! JFNK_RAD_SWEEP_V1 consumes it in its epilogue by forcing a sentinel
! residual (SE%STEQ), which JFNK_RESIDUAL_V1's |Y|>1D30 cap turns into
! jfnk_resid_nbad>0, so the driver falls back to LAMBDA/BA.  Cleared at
! each sweep entry; meaningless outside JFNK_IN_EVAL=T regions.
	LOGICAL JFNK_EVAL_BROKEN
	INTEGER JFNK_MAX_NEWTON
	INTEGER JFNK_M_GMRES
	INTEGER JFNK_MAX_GMRES
	REAL(KIND=LDP) JFNK_GTOL
	REAL(KIND=LDP) JFNK_TOL_REL
	REAL(KIND=LDP) JFNK_EPS_JVP
	REAL(KIND=LDP) JFNK_DEL_INIT
	REAL(KIND=LDP) JFNK_MNDL
	REAL(KIND=LDP) JFNK_MXDL
	LOGICAL JFNK_FIXED_T
	LOGICAL JFNK_VERBOSE
	LOGICAL JFNK_DIAGNOSTIC
	INTEGER JFNK_WARMUP
	REAL(KIND=LDP) MAX_JFNK_COR
	LOGICAL JFNK_APPLY_BA_RIDGE
	REAL(KIND=LDP) JFNK_BA_RIDGE_FAC
	LOGICAL JFNK_PRECOND_REFINE
	LOGICAL JFNK_LEFT_PRECOND
	REAL(KIND=LDP) JFNK_PC_REFRESH_FRAC
	CHARACTER(LEN=8) JFNK_PRECOND_TYPE
	REAL(KIND=LDP) JFNK_LAM_AUTO_THRESH
	REAL(KIND=LDP) JFNK_COR_GATE
! JFNK_FRZ_GATE (VADAT [JFNK_FRZ_GATE], default 0=off): two-stage
! admission for the cheaper frozen-J Newton (JFNK_COUPLED=F residual
! semantics -- no radiation re-solve per eval).  When set, JFNK fires
! with the frozen residual once the demanded correction drops below
! FRZ_GATE; the coupled eval is admitted only when MAXCH also falls
! below JFNK_COR_GATE.  Default 0 reproduces the original single-gate
! behavior (one admission point, mode = [JFNK_COUPLED]).
	REAL(KIND=LDP) JFNK_FRZ_GATE
! JFNK_GATE_NREQ (VADAT [JFNK_GATE_NREQ], default 1): hysteresis on
! the admission gates -- Newton fires only after this many
! consecutive correction solves below the applicable gate, so a
! single lucky low-MAXCH Lambda step cannot admit JFNK while the
! field is still drifting.  Default 1 reproduces the instantaneous
! gate.
	INTEGER JFNK_GATE_NREQ
! JFNK_DEEPDIAG (VADAT [JFNK_DEEPDIAG], default F): run the
! compression-zone diagnostic suite (linearity probe, FD-eps sweep,
! FD-vs-BA audit, T-column audit) once per process at the first JFNK
! activation.  Costs ~12 residual evaluations -- a full formal sweep
! each under JFNK_COUPLED=T -- so it is opt-in only.
	LOGICAL JFNK_DEEPDIAG
! JFNK_MAX_EVALS (VADAT [JFNK_MAX_EVALS], default 0=off): per-
! activation residual-evaluation budget.  A coupled eval is a full
! formal sweep; an activation that exhausts the budget gets a
! sentinel residual so Newton rejects and the driver falls back to
! the best-iterate/Lambda path instead of burning sweeps.
	INTEGER JFNK_MAX_EVALS
	LOGICAL JFNK_LINE_SEARCH
	INTEGER JFNK_LS_MAX
	REAL(KIND=LDP) JFNK_LS_BETA
	REAL(KIND=LDP) JFNK_LS_C
	LOGICAL JFNK_CAUCHY
	CHARACTER(LEN=12) JFNK_RSCALE_MODE
	INTEGER JFNK_MAX_RESTART
	LOGICAL JFNK_USE_BANDED
	LOGICAL JFNK_SELECTIVE_JVAR
	REAL(KIND=LDP) JFNK_ROW_TINY
	REAL(KIND=LDP) JFNK_LSQ_LAM
	REAL(KIND=LDP) JFNK_REACH_GATE
	INTEGER JFNK_REARM_EVERY
	INTEGER JFNK_ROW_DIAG
!
! JFNK_TAUTO: adaptive damped-T domain inside the JFNK solve.  When
! active the VARFIXT depth range is not hard-masked: the T rows of F
! carry the energy residual, the JVP differentiates T, and every
! accepted Newton step has its T components scaled by the current
! ladder damping (JFNK_TAU_ALPHA in JFNK_CMFGEN_INT_MOD) and capped
! at JFNK_TAU_DMAX.  This is the STELLA-shock guard: a STELLA T
! profile is not in NLTE equilibrium, so the energy residual starts
! huge and a free Newton step would drive T to unphysical values.
! JFNK_TAUTO defaults to VARFIXT (VADAT FIX_T_AUTO) so the existing
! key switches the JFNK behaviour too; JFNK_TAUTO in VADAT overrides.
!
	LOGICAL JFNK_TAUTO
	REAL(KIND=LDP) JFNK_TAU_A0	! first damping rung
	REAL(KIND=LDP) JFNK_TAU_REL	! |RE|/|RE0| release fraction
	REAL(KIND=LDP) JFNK_TAU_EPS	! |dlogT|inf revert threshold
	REAL(KIND=LDP) JFNK_TAU_DMAX	! hard |dlogT| cap per step
	INTEGER JFNK_TAU_NLOCK	! commits before release allowed
	LOGICAL JFNK_TAU_PCW	! down-weight T rows in PC'd system
!
	END MODULE CONTROL_VARIABLE_MOD
