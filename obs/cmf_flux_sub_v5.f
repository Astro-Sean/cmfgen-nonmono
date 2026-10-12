!
! Subroutine to compute the emergent flux in the observer's frame using an
! Observer's frame flux calculation. The program first computes the
! emissivities, opacities, and radiation field in the CMF frame. COHERENT
! or INCOHERENT electron scattering can be assumed. This data is then passed
! to OBS_FRAME_SUB for computation of the OBSERVER's flux.
!
! This routine is a heavily stripped down and modified version of CMFGEN.
!
	SUBROUTINE CMF_FLUX_SUB_V5(ND,NC,NP,NDMAX,NPMAX,NT,NLINE_MAX)
	USE SET_KIND_MODULE
	USE DUST_MOD
	USE MOD_CMF_OBS
	USE MOD_FREQ_OBS
	USE CMF_FLUX_CNTRL_VAR_MOD
	USE MOD_LEV_DIS_BLK
	USE EDDFAC_REC_DEFS_MOD
	IMPLICIT NONE
!
! Altered: 12-Oct-2026 : Added DO_SOB_FLUX_OBS path: the Sobolev EWs
!                          computed in the DO_SOBOLEV_LINES block are
!                          deposited into OBS_FLUX (DEP_SOB_OBS_FLUX) so
!                          SOB lines reach the observed spectrum.
! Altered: 05-Oct-2026 : Validate V ordering against the selected
!                          non-monotonic mode before spectrum setup (Devin).
! Altered: 22-Jan-2024 : Fixed bug - dust and Rayleigh scattering were being treated like electron 
!                          scattering in the observer's frame calcualtion.
! Altered: 19-Mar-2024 : Added routines for dust treatment. All routines will treat isotropic 
!                          scattering. Only formal rel treats scattering with a HG phase function.
!                          Altered SRCE check for neg opacities (issue when dust added -- needs fixing)
! Altered: 28-Aug-2022 : Fixed bug - MAX_SIM now set for SOB computation with continuum calculation.
! Altered: 24-Aug-2022 : Added SOB_EW_LAM_BEG (and _END) to make Sobolev EW calculation more transparent
! Altered: 05-Jul-2022 : MAX_SIM is now determined by model, and vectors using
!                          this variable are now allocated.
!                        Now use INIT_PROFILE_MODULE_V2, SET_PROF_V6,
!                           and GET_PROFILE_STORAGE_LIMITS_V2.
! Altered: 27-Feb-2022 : Code no linger computes the opacities and emissivities on the
!                          2nd (and higher) electron scattering iterations. Instead it
!                          uses the data store in ETA_CMF_ST, CHI_CMF_ST
! Altered: 30-Apr-2020 : RAY_DATA created when Rayleigh scattering included.
! Altered: 16-Aug-2019 : PLANCKMEAN computed
! Altered: 26-Apr-2019 : Added variables/option to restrict range of OBS frame calcualtion.
!                           Added option to compute EW inc CMF (include ability to
!                               compute Int ABS(I/Ic-1) dv.
! Altered: 21-Sep-2016 : Error corrected -- I was writing out KAPPA from RVTJ rather than the from the CMF_FLUX
!                           calculaton. Note that the ROSSELAND mean is only correct (at depth) if we compute
!                           the spectrum over the full frequency range.
! Altered: 09-Sep-2015 : Changed to C4(line)= ABS(C4[upper]) + ABS(C4[lower) [I was just summing values]
! Altered: 18-May-2015 : Changed GAM2, GAM4 to C4 and C6 (quadratic and Van der Waals interacton constants).
!                           C4 is now utilized (read into VEC_C4). C6 is still not used (09-Jun-2015).
! Altered: 04-Apr-2015 : Changed SET_TWO_PHOT_V2 to _V3.
! Altered: 21-Jan-2013 : Change to vectors passed to SET_PRO_V3. Error probably affects IR HeI lines.
!                           No change top UV/optical spectrum was seen.
!                        Placed large vectors (dimension with NCF_MAX and NLINE_MAX) in module MOD_FREQ_OBS.
!                           These were being placed on the stack, and limiting the available stacks size.
!                           Error in STACK found using -g debug option with -Mchkstk.
!
! Altered:  1-Oct-2012 : Fixed bug: NC_OBS was not being set for the normal P-grid which caused
!                           a problem with the observer's frame calculation.
! Altered:  3-Aug-2012 : Changed to allow the computation of a revised grid for shell models.
! Altered: 17-Dec-2011 : Now call OPACITIES_V5.INC and EVAL_LTE_INC_V5.INC
!                          L_STAR_RATIO and U_STAR_RATIO now computed using XzVLTE_F_ON_S
!                          Done to allow lower wind temperatures.
! Altered: 20-Oct-2009 : Changed to call OBS_FRAME_SUB_V8 (WRITE_RTAU, TAU_REF added to call).
! Altered: 07-Jul-2008 : Substantial changes over a period of time to include
!                          relativistic radiative transfer. Lots of changes to
!                          COMP_JCONT.INC. NCEXT is now increased according to
!                          NPINS when ACCURATE option is used.
! Altered: 31-Jul-2007 : Call to MOM_JREL_V3 included in INCLUDE file.
!                          Variabled and control variables for rel. trans. installed.
! Altered: 31-Jul-2007 : Call to MOM_JREL_V3 included in INCLUDE file.
!                          Variabled and control variables for rel. trans. installed.
! Altered: 01-Feb-2006 : COMP_OBS_V2 installed. Frequency transformation between
!                          comoving and observing franes is now accuratle for all speeds.
!                          Frequency scaling of intensities currently switched off.
! Altered: 03-May-2004 : Options are now read into store in CMF_FLUX.
!                          STORE cleaned after use.
! Altered: 05-Jun-2002 : Min bug fix: Changed value of VDOP passed to
!                          DET_MAIN_CONT_FREQ. Reduces number of continuum frequencies.
! Altered: 16-Aug-2002 : Bug fix, NU_OBS_MAX was not always chosen correctly.
! Altered: 03-Jan-2001 : PROF_TYPE and VEC_CHAR_WRK now of length 12.
! Altered: 10-Mar-2000 : Variable type ATM installed to simplify access to
!                          different model atoms. Installation reduces size
!                          of executable, and allows for faster compilation.
!                          Changed to V4.
! Altered: 16-Jan-2000 : Changed to V2; Call to OBS_FRAME_V2 changed.
! Altered: 04-Jan-2000 : We now check that arrays can be allocated.
! Finalized: 5-Jan-1999
!
! Maximum number of lines whose profile overlap at a given frequency.
!
	INTEGER MAX_SIM
!
! Define vinteger variables to set aside for storage for the intrinsic line profiles.
! These are now computed by a subroutine (as of Nov-2013).
!
	INTEGER NLINES_PROF_STORE,NLINES_VOIGT_STORE
	INTEGER NFREQ_PROF_STORE,NFREQ_VOIGT_STORE
!
	INTEGER NCF
	INTEGER ND,NC,NP
	INTEGER NDMAX,NPMAX
	INTEGER NT,NLINE_MAX
!
	REAL(KIND=LDP) POPS(NT,ND)
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
! Internally used variables
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) OPLIN,EMLIN
	REAL(KIND=LDP) DTDR,BNUE,DBB,DDBBDT
	REAL(KIND=LDP) MEAN_ATOMIC_WEIGHT
	REAL(KIND=LDP) TSTAR,S1,IC,MAXCH
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) T1,T2,T3,T4
	REAL(KIND=LDP) CHI_FLR		!positive opacity floor for guarded divisions
	REAL(KIND=LDP) FL,FL_OLD
	REAL(KIND=LDP) FG_COUNT
!
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
 	INTEGER REC_SIZE
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER N_PER_REC
!
! 
! Logical Unit assignments. Those indicated with a # after the ! are open in
!  large sections of the code. Other units generally used temprarily.
!
	INTEGER               LUER        !Output/Error file.
	INTEGER, PARAMETER :: LUIN=7      !General input unit (closed after accesses).
	INTEGER, PARAMETER :: LUMOD=8     !Model Description file.
!
	INTEGER, PARAMETER :: LU_FLUX=10   	!Flux/Luminosity Data (OBSFLUX)
	INTEGER, PARAMETER :: LU_OPAC=18   	!Rosseland mean opacity etc.
	INTEGER, PARAMETER :: LU_EW=20     	!# EW data.
!
	INTEGER, PARAMETER :: LU_EDD=35       !Continuum Eddington factors.
	INTEGER, PARAMETER :: LU_JCOMP=37     !J_COMP
	INTEGER, PARAMETER :: LU_ES=38        !ES_J_CONV
!
! For listing of transitions with TOTAL negative opacity values at some depths.
!
	INTEGER, PARAMETER :: LU_NEG=75
!
! 
!
	INTEGER NL,NUP
	INTEGER MNL,MNUP
	INTEGER MNL_F,MNUP_F
	INTEGER I,J,K,L,ML,LS,LINE_INDX
	INTEGER ID,ID_SAV,ISPEC
	INTEGER ES_COUNTER
!
! Functions called
!
	INTEGER ICHRLEN,ERROR_LU
	REAL(KIND=LDP) LAMVACAIR
	REAL(KIND=LDP) ATOMIC_MASS_UNIT
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL JWEIGHT_V2,HWEIGHT_V2,KWEIGHT_V2,NWEIGHT_V2
	EXTERNAL JTRPWGT_V2,HTRPWGT_V2,KTRPWGT_V2,NTRPWGT_V2
	EXTERNAL ICHRLEN,ERROR_LU,SPEED_OF_LIGHT
	EXTERNAL LAMVACAIR,ATOMIC_MASS_UNIT
!
! 
	CHARACTER FMT*120
	CHARACTER SECTION*20
	CHARACTER TMP_KEY*20
	CHARACTER STRING*132
	CHARACTER EW_STRING*160
	CHARACTER TEMP_CHAR*132
!
! Global vectors:
!
	REAL(KIND=LDP) AMASS_ALL(NT)
!
! Arrays and variables for treating lines simultaneously.
!
	REAL(KIND=LDP), ALLOCATABLE :: EINA(:)
	REAL(KIND=LDP), ALLOCATABLE :: OSCIL(:)
	REAL(KIND=LDP), ALLOCATABLE :: GLDGU(:)
	REAL(KIND=LDP), ALLOCATABLE :: AMASS_SIM(:)
	REAL(KIND=LDP), ALLOCATABLE :: FL_SIM(:)
	INTEGER, ALLOCATABLE :: SIM_NL(:)
	INTEGER, ALLOCATABLE :: SIM_NUP(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: CHIL_MAT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETAL_MAT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: BB_COR(:,:)
!
	INTEGER NUM_SIM_LINES
	INTEGER SIM_INDX
	INTEGER TMP_MAX_SIM
	REAL(KIND=LDP) OVER_FREQ_DIF
	REAL(KIND=LDP) EW,ABS_EW
	REAL(KIND=LDP) CONT_INT
!
! Accumulates the Sobolev line flux deposited into the observed spectrum
! when DO_SOB_FLUX_OBS is true (see DEP_SOB_OBS_FLUX).
!
	REAL(KIND=LDP), ALLOCATABLE :: SOB_DEP(:)
	LOGICAL OVERLAP
	LOGICAL SOBOLEV
!
! L refers to the lower level, U to the upper level.
!
	REAL(KIND=LDP), ALLOCATABLE :: L_STAR_RATIO(:,:)   !ND,MAX_SIM)
	REAL(KIND=LDP), ALLOCATABLE :: U_STAR_RATIO(:,:)   !ND,MAX_SIM)
!
	REAL(KIND=LDP), ALLOCATABLE :: ETA_CMF_ST(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_CMF_ST(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RJ_CMF_ST(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RAY_CMF_ST(:,:)			!Rayleigh scattering opacity
	REAL(KIND=LDP), ALLOCATABLE :: DUST_CMF_ST(:,:)			!Dust scattering
	REAL(KIND=LDP), ALLOCATABLE :: ION_LINE_FORCE(:,:)
!
! Vectors for treating lines simultaneously with the continuum.
!
	REAL(KIND=LDP), ALLOCATABLE :: LINE_PROF_SIM(:,:)   !ND,MAX_SIM)
!
	REAL(KIND=LDP) NU_MAX_OBS
	REAL(KIND=LDP) NU_MIN_OBS
!
	REAL(KIND=LDP) CONT_FREQ
!
	CHARACTER*50, ALLOCATABLE :: TRANS_NAME_SIM(:)
!
	LOGICAL, ALLOCATABLE :: RESONANCE_ZONE(:)
	LOGICAL, ALLOCATABLE ::  END_RES_ZONE(:)
	LOGICAL, ALLOCATABLE ::  LINE_STORAGE_USED(:)
!
	INTEGER FREQ_INDX
	INTEGER FIRST_LINE
	INTEGER LAST_LINE
	INTEGER, ALLOCATABLE :: SIM_LINE_POINTER(:)
!
! 
!
! Opacity/emissivity
	REAL(KIND=LDP) CHI(ND)			!Continuum opacity (all sources)
	REAL(KIND=LDP) CHI_RAY(ND)
	REAL(KIND=LDP) CHI_SCAT(ND)
!
	REAL(KIND=LDP) ETA_DUST(ND)
	REAL(KIND=LDP) CHI_ABS_DUST(ND)
	REAL(KIND=LDP) CHI_SCAT_DUST(ND)
!
	REAL(KIND=LDP) ETA(ND)			!Continuum emissivity (all sources)
	REAL(KIND=LDP) CHIL(ND)			!Line opacity (without prof.)
	REAL(KIND=LDP) ETAL(ND)			!Line emissivity (without prof.)
	REAL(KIND=LDP) ESEC(ND)			!Continuum electron scattering coef.
	REAL(KIND=LDP) ZETA(ND)			!Source func. (all except elec. scat.)
	REAL(KIND=LDP) THETA(ND)		!Elec. scat. source coef.
	REAL(KIND=LDP) SOURCE(ND)		!Complete source function.
	REAL(KIND=LDP) DTAU(NDMAX)		!Optical depth (used in error calcs)
!		DTAU(I)=0.5*(CHI(I)+CHI(I+1))*(Z(I)-Z(I+1))
	REAL(KIND=LDP) dCHIdR(NDMAX) 		!Derivative of opacity.
!
	REAL(KIND=LDP) P(NP)
!
	REAL(KIND=LDP) CHI_CONT(ND)
	REAL(KIND=LDP) ETA_CONT(ND)
!
! These parameters are used when computing J and the variation of J.
!
	REAL(KIND=LDP) CHI_SCAT_DUST_CLUMP(ND)
	REAL(KIND=LDP) CHI_SCAT_CLUMP(ND)
	REAL(KIND=LDP) CHI_RAY_CLUMP(ND)
	REAL(KIND=LDP) CHI_CLUMP(ND)		!==CHI(I)*CLUMP_FAC(I)
	REAL(KIND=LDP) ETA_CLUMP(ND)		!==ETA(I)*CLUMP_FAC(I)
	REAL(KIND=LDP) ESEC_CLUMP(ND)		!==ESEC(I)*CLUMP_FAC(I)
!
! Variables to limit the computation of the continuum opacities and
! emissivities.
!
	REAL(KIND=LDP) EMHNUKT_CONT(ND)
	REAL(KIND=LDP) ETA_C_EVAL(ND)
	REAL(KIND=LDP) CHI_C_EVAL(ND)
!
! 
!
	REAL(KIND=LDP) Z_POP(NT)		!Ionic charge for each species
!
! Variables etc for computation of continuum in comoving frame.
!
	LOGICAL FIRST_FREQ
	LOGICAL NEW_FREQ
	REAL(KIND=LDP) dLOG_NU			!Step in frequency in Log plane
	REAL(KIND=LDP) FEDD_PREV(NDMAX)
	REAL(KIND=LDP) GEDD_PREV(NDMAX)
	REAL(KIND=LDP) H_ON_J(NDMAX)
	REAL(KIND=LDP) N_ON_J_NODE(NDMAX)
	REAL(KIND=LDP) dlnJdlnR(NDMAX)
	REAL(KIND=LDP) KMID_ON_J(NDMAX)
	REAL(KIND=LDP) N_ON_J(NDMAX)
	REAL(KIND=LDP) N_ON_J_PREV(NDMAX)
	REAL(KIND=LDP) JNU_PREV(NDMAX)
	REAL(KIND=LDP) RSQHNU_PREV(NDMAX)
	REAL(KIND=LDP) GEDD(NDMAX)
	REAL(KIND=LDP) RSQHNU(NDMAX)
	REAL(KIND=LDP) CHI_PREV(ND)
	REAL(KIND=LDP) ETA_PREV(ND)
	REAL(KIND=LDP) HBC_CMF(3),HBC_PREV(3)
	REAL(KIND=LDP) NBC_CMF(3),NBC_PREV(3),INBC_PREV
	REAL(KIND=LDP) HFLUX_AT_IB,HFLUX_AT_OB
!
! Quadrature weights.
	REAL(KIND=LDP) AQW(ND,NP)		!Angular quad. weights. (indep. of v)
	REAL(KIND=LDP) HQW(ND,NP)		!Angular quad. weights. (indep. of v) for flux integration.
	REAL(KIND=LDP) KQW(ND,NP)		!Angular quad. weights for K integration.
	REAL(KIND=LDP) NQW(ND,NP)		!Angular quad. weights. (indep. of v) for N integration.
	REAL(KIND=LDP) HMIDQW(ND,NP)		!Angular quad. weights. (indep. of v)
                                        !for flux integration. Defined at the
                                        !mid points of the radius mesh.
	REAL(KIND=LDP) NMIDQW(ND,NP)		!Angular quad. weights. (indep. of v)
                                        !for N integration. Defined at the
                                        !mid points of the radius mesh.
!
! Continuum matrices
	REAL(KIND=LDP) WM(ND,ND)		!Coef. matrix of J & %J vector
	REAL(KIND=LDP) FB(ND,ND)		!Coef. of J & %J vects in angular equ.
!
! Transfer equation vectors
	REAL(KIND=LDP) TA(NDMAX)
	REAL(KIND=LDP) TB(NDMAX)
	REAL(KIND=LDP) TC(NDMAX)
	REAL(KIND=LDP) XM(NDMAX)		!R.H.S. (SOURCE VECTOR)
	REAL(KIND=LDP) TCHI(ND)
!
! 
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	REAL(KIND=LDP) FEDD(NDMAX)
	REAL(KIND=LDP) QEDD(NDMAX)
!
	LOGICAL MID
	LOGICAL FIRST
	LOGICAL NEG_OPACITY(ND)
	LOGICAL FIRST_NEG
	LOGICAL LAMBDA_ITERATION
	LOGICAL LST_ITERATION
	LOGICAL LST_DEPTH_ONLY
!
!
!
! X-ray variables.
!
	REAL(KIND=LDP) FILL_VEC_SQ(ND)
	REAL(KIND=LDP) XRAY_LUM(ND)
	REAL(KIND=LDP) GFF,XCROSS_V2
	EXTERNAL GFF,XCROSS_V2
!
!
!
! ACESS_F is the current record we are writing in EDDFACTOR.
!
	INTEGER ACCESS_F
	CHARACTER(LEN=20) DA_FILE_DATE
!
	INTEGER NDEXT,NCEXT,NPEXT
	INTEGER INDX(NDMAX),POS_IN_NEW_GRID(ND)
	REAL(KIND=LDP) COEF(0:3,NDMAX)
	REAL(KIND=LDP) INBC,HBC_J,HBC_S			!Bound. Cond. for JFEAU
!
	REAL(KIND=LDP) REXT(NDMAX),PEXT(NPMAX),VEXT(NDMAX),LANG_COORDEXT(NDMAX)
	REAL(KIND=LDP) TEXT(NDMAX),SIGMAEXT(NDMAX)
	REAL(KIND=LDP) CHIEXT(NDMAX),ESECEXT(NDMAX),ETAEXT(NDMAX)
	REAL(KIND=LDP) CHI_RAY_EXT(NDMAX),CHI_SCAT_EXT(NDMAX)
	REAL(KIND=LDP) CHI_SCAT_DUST_EXT(NDMAX)
	REAL(KIND=LDP) ZETAEXT(NDMAX),THETAEXT(NDMAX)
	REAL(KIND=LDP) RJEXT(NDMAX),RJEXT_ES(NDMAX)
	REAL(KIND=LDP) FOLD(NDMAX),FEXT(NDMAX),QEXT(NDMAX),SOURCEEXT(NDMAX)
	REAL(KIND=LDP) VDOP_VEC_EXT(NDMAX)
!
!
	REAL(KIND=LDP) F2DAEXT(NDMAX,NDMAX)     !These arrays don't need to be
!
! If required, these arrays shoukd have size NDEXT*NPEXT
!
	REAL(KIND=LDP), ALLOCATABLE :: AQWEXT(:,:)	!Angular quad. weights. (indep. of v)
	REAL(KIND=LDP), ALLOCATABLE :: HQWEXT(:,:)	!Angular quad. weights for flux integration.
	REAL(KIND=LDP), ALLOCATABLE :: KQWEXT(:,:)	!Angular quad. weights for K integration.
	REAL(KIND=LDP), ALLOCATABLE :: NQWEXT(:,:)	!Angular quad. weights for N integration.
	REAL(KIND=LDP), ALLOCATABLE :: HMIDQWEXT(:,:)	!Angular quad. weights for flux integration.
	REAL(KIND=LDP), ALLOCATABLE :: NMIDQWEXT(:,:)	!Angular quad. weights for flux integration.
!
! Arrays for calculating mean opacities.
!
	REAL(KIND=LDP) FLUXMEAN(ND) 		!Flux mean opacity
	REAL(KIND=LDP) LINE_FLUXMEAN(ND) 	!Flux mean opacity due to lines.
	REAL(KIND=LDP) ROSSMEAN(ND)  		!Rosseland mean opacity
	REAL(KIND=LDP) PLANCKMEAN(ND)  	!Planck mean opacity
	REAL(KIND=LDP) INT_dBdT(ND)  		!Integral of dB/dT over nu (to calculate ROSSMEAN)
	REAL(KIND=LDP) FORCE_MULT(ND)
	REAL(KIND=LDP) NU_FORCE
	REAL(KIND=LDP) NU_FORCE_FAC
	INTEGER N_FORCE
	INTEGER ML_FORCE
	LOGICAL TMP_LOG
!
! Other arrays
	REAL(KIND=LDP) Z(NDMAX)			!Z displacement along a given array
	REAL(KIND=LDP) EMHNUKT(ND)		!EXP(-hv/kT)
	REAL(KIND=LDP) RLUMST(ND)		!Luminosity as a function of depth
	REAL(KIND=LDP) J_INT(ND)		!Frequency integrated J
	REAL(KIND=LDP) K_INT(ND)		!Frequency integrated K
	REAL(KIND=LDP) K_MOM(ND)		!Frequency dependent K moment
	REAL(KIND=LDP) SOB(ND)   	    	!Used in computing continuum flux
	REAL(KIND=LDP) RJ(ND)			!Mean intensity
	REAL(KIND=LDP) RJ_ES(ND)		!Convolution of RJ with e.s. R(v'v')
!
! Line variables.
!
	REAL(KIND=LDP) VAL_DO_NG
	REAL(KIND=LDP) RP
	REAL(KIND=LDP) VINF
	REAL(KIND=LDP) VTURB_VEC(ND)
	REAL(KIND=LDP) VDOP_VEC(ND)
	REAL(KIND=LDP) MAX_DEL_V_RES_ZONE(ND)
!
! Parameters, vectors, and arrays for computing the observed flux.
!
	INTEGER, PARAMETER :: NST_CMF=20000
	INTEGER NP_OBS_MAX
	INTEGER NP_OBS
	INTEGER NC_OBS
	REAL(KIND=LDP)  NU_STORE(NST_CMF)
	REAL(KIND=LDP) V_AT_RMAX		!Used if we extend the atmosphere.
	REAL(KIND=LDP) RMAX_OBS
	REAL(KIND=LDP) H_OUT,H_IN
!
! We allocate memory for the following vectors as we use them for the regular
! flux computation, and when extra depth points are inserted (ACCURATE=.TRUE.)
!
	REAL(KIND=LDP), ALLOCATABLE :: IPLUS_STORE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: P_OBS(:)
	REAL(KIND=LDP), ALLOCATABLE :: IPLUS(:)
	REAL(KIND=LDP), ALLOCATABLE :: MU_AT_RMAX(:)
	REAL(KIND=LDP), ALLOCATABLE :: dMU_AT_RMAX(:)
	REAL(KIND=LDP), ALLOCATABLE :: HQW_AT_RMAX(:)
!
! Supercedes OBS
!
	INTEGER N_OBS
	LOGICAL FIRST_OBS_COMP
!
! Indicates approximate frequencies for which TAU at outer boundary is written
! to OUTGEN on the last iteration.
!
! They are the He2 ege, NIII/CIII egde, HeI, HI, HI(N=2).
!
	INTEGER, PARAMETER :: N_TAU_EDGE=5
	REAL(KIND=LDP) TAU_EDGE(N_TAU_EDGE)
	DATA TAU_EDGE/13.16_LDP,11.60_LDP,5.95_LDP,3.29_LDP,0.83_LDP/
!
!
! Check whether EQUATION LABELLING is consistent. ' I ' is used as the
! number of the current equation. We also set the variable SPEC_PRES which
! indicates whether at least one ioization stage of a species is present.
! It is used to determine, foe example,  whether a number conservation
! equation is required.
!
!
!???????????????????????
!
	I=1
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    SPECIES_PRES(ISPEC)=.FALSE.
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      CALL CHK_EQ_NUM( ATM(ID)%XzV_PRES, ATM(ID)%EQXzV,
	1                      ATM(ID)%NXzV,
	1                      I,SPECIES_PRES(ISPEC),TRIM(ION_ID(ID)))
	    END DO
	    CALL CHK_EQ_NUM(SPECIES_PRES(ISPEC),EQ_SPECIES(ISPEC),IONE,I,
	1                          SPECIES_PRES(ISPEC),TRIM(SPECIES(ISPEC)))
	    END IF
	END DO
!
	LUER=ERROR_LU()
	IF(EQNE .NE. I)THEN
	  WRITE(LUER,*)'Error - EQNE has wrong value in CMF_FLUX_SUB_V5'
	  STOP
	END IF
	IF(NT .NE. I+1)THEN
	  WRITE(LUER,*)'Error - NT has wrong value in CMF_FLUX_SUB_V5'
	  STOP
	END IF
!
	CALL INIT_MOD_FREQ_OBS(NLINE_MAX)
! 
!
!	LAMBDA_ITERATION=.FALSE.
	LAMBDA_ITERATION=.TRUE.
	LST_ITERATION=.FALSE.
	LST_DEPTH_ONLY=.FALSE.
	RD_STARK_FILE=.FALSE.
	ACCESS_F=5
	RP=R(ND)
	VINF=V(1)
	TSTAR=T(ND)
	C_KMS=SPEED_OF_LIGHT()/1.0E+05_LDP
	DO_REL_IN_OBSFRAME=.FALSE.
	DA_FILE_DATE='20-Aug-2003'
!
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
! MAXCH and VAL_DO_NG are set so that they defined for the TEST in
! COMP_JCONT_V?.INC whether to do an accurate flux calculation. An
! accurate flux calculation can be avoided by doing a LAMBDA iteration.
!
	MAXCH=0.0_LDP
	VAL_DO_NG=5.0_LDP
!
! Set the vector Z_POP to contain the ionic charge for each species.
!
	Z_POP(1:NT)=0.0_LDP
!
	DO ID=1,NUM_IONS
	  CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1                       ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
	END DO
!
! Store atomic masses in vector of LENGTH NT for later use by line
! calculations. G_ALL and LEVEL_ID  are no longer used due to the use
! of super levels.
!
	AMASS_ALL(1:NT)=0.0_LDP
!
! We also set the mass of the ion corresponding to XzV to AT_MASS for each
! species . Thus the range is EQXzV:EQXzV+NXzV, NOT EQXzV:EQXzV+NXzV-1. Note
! that HeIII_PRES is always false, even when HeI and HeII are both present.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)AMASS_ALL(ATM(ID)%EQXzV:ATM(ID)%EQXzV+ATM(ID)%NXzV)=
	1              AT_MASS(SPECIES_LNK(ID))
	END DO
!
!
! Read in all options controlling the spectrum computation.
! As subroutine exits, it also cleans and releases the STORE.
!
	  CALL RD_CMF_FLUX_CONTROLS(ND,LUMOD,LUER)
!
! The stock CMF solvers are only validated for the legacy positive,
! decreasing V(R).  A dedicated non-monotonic flag selects the guarded
! path; without one, reject a repeated or reversed moving atmosphere
! rather than silently taking the wrong formal solution.  Static
! atmospheres remain allowed.
!
	  IF(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB)THEN
	    DO I=1,ND
	      IF(V(I) .LE. 0.0_LDP)THEN
	        WRITE(LUER,*)'Error in CMF_FLUX_SUB_V5:'
	        WRITE(LUER,*)'Non-monotonic modes require positive V.'
	        WRITE(LUER,*)'Depth, V=',I,V(I)
	        STOP
	      END IF
	    END DO
	  ELSE IF(MAXVAL(ABS(V(1:ND))) .GT. 0.0_LDP)THEN
	    DO I=1,ND-1
	      IF(V(I+1) .GE. V(I))THEN
	        WRITE(LUER,*)'Error in CMF_FLUX_SUB_V5:'
	        WRITE(LUER,*)'Non-monotonic V requires DO_NONMONO_SOB or'
	        WRITE(LUER,*)'DO_NONMONO_KNOP in CMF_FLUX_INPUT.'
	        WRITE(LUER,*)'Depth, V=',I,V(I),V(I+1)
	        STOP
	      END IF
	    END DO
	  END IF
!
! Stock CMFGEN assumes V(1) is the fastest material.  In a shocked
! CSM model V(1) is the slow outer CSM (~100 km/s) while the shell and
! ejecta interior move at several 1000 km/s.  VINF sets the CMF-to-
! observer frequency coverage margins (NU_MAX_OBS/NU_MIN_OBS) and the
! e.s. wing extents in INS_LINE_V5/INS_LINE_OBS_V4, so an undersized
! VINF leaves Doppler-shifted opacity outside the refined grid
! regions.  Use the full velocity-magnitude range on the non-monotonic paths only.
! Must sit after RD_CMF_FLUX_CONTROLS, which sets the DO_NONMONO flags.
!
	IF(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB)THEN
	  VINF=MAXVAL(ABS(V(1:ND)))
	  WRITE(LUER,'(A,F10.2,A,F10.2,A)')
	1   ' Non-monotonic V: VINF set to MAXVAL(ABS(V))=',VINF,
	2   ' km/s (V(1)=',V(1),')'
	END IF
!
! 
!
! Scale abunances of species if desired. This should only be done for exploratory
! spectral calculations, and only for IMPURITY species (i.e. not H or He).
!
	  K=0
	  DO ISPEC=1,NUM_SPECIES
	    IF(SPECIES_PRES(ISPEC) .AND. ABUND_SCALE_FAC(ISPEC) .NE. 1.0_LDP)THEN
	      IF(K .EQ. 0)THEN
	        WRITE(LUER,'(A)')' '
	        WRITE(LUER,'(A)')'******************Warning**********************'
	      END IF
	      K=1
	      WRITE(LUER,'(A,A,A)')'Abundance of species ',TRIM(SPECIES(ISPEC)),' scaled'
	      AT_ABUND(ISPEC)=AT_ABUND(ISPEC)*ABUND_SCALE_FAC(ISPEC)
	      POP_SPECIES(:,ISPEC)=POP_SPECIES(:,ISPEC)*ABUND_SCALE_FAC(ISPEC)
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        ATM(ID)%XzV_F=ATM(ID)%XzV_F*ABUND_SCALE_FAC(ISPEC)
	        ATM(ID)%DXzV_F=ATM(ID)%DXzV_F*ABUND_SCALE_FAC(ISPEC)
	      END DO
	    END IF
	  END DO
	  IF(K .NE. 0)THEN
	    WRITE(LUER,'(A)')'******************Warning**********************'
	    WRITE(LUER,'(A)')' '
	  END IF
!
! Evaluate constants used for level-dissolution calculations.
!
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! Read in file contaning line links to STARK tables.
!
	IF(RD_STARK_FILE .OR. GLOBAL_LINE_PROF(1:4) .EQ. 'LIST')THEN
	  CALL RD_STRK_LIST(LUIN)
	END IF
!
! Read in atomic data for 2-photon transitions.
!
	CALL RD_TWO_PHOT(LUIN,INCL_TWO_PHOT)
!
! Read in X-ray photoionization cross-sections.
!
	CALL RD_XRAY_FITS(LUIN)
!
	IF(XRAYS .AND. .NOT. FF_XRAYS)THEN
	  CALL RD_XRAY_SPEC(T_SHOCK_1,T_SHOCK_2,LUIN)
	END IF
!
	IF(INCL_DUST)THEN
	  CALL RD_DUST_CNTRL_FILE()
	  CALL SET_DUST_EJECTA_PROPERTIES(R,V,DENSITY,CLUMP_FAC,ND)
	END IF
!
! We now need to compute the populations for the model atom with Super-levels.
! We do this in reverse order (i.e. highest ionization stage first) in order
! that we the ion density for the lower ionization stage is available for
! the next call.
!
! For 1st call to FULL_TO_SUP, Last line contains FeX etc as FeXI not installed.
! We only do to NUM_IONS-1 to avoid access arror, and since the ion
! corresponding to NUM_IONS contains only 1 level and is done with NUM_IONS-1
!
	  DO ID=NUM_IONS-1,1,-1
	    CALL FULL_TO_SUP(
	1        ATM(ID)%XzV,        ATM(ID)%NXzV,   ATM(ID)%DXzV,
	1        ATM(ID)%XzV_PRES,   ATM(ID)%XzV_F,
	1        ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,
	1        ATM(ID+1)%XzV,      ATM(ID+1)%NXzV, ATM(ID+1)%XzV_PRES, ND)
	  END DO
!
! Store all quantities in POPS array.
!
	  DO ID=1,NUM_IONS
	    CALL IONTOPOP(POPS, ATM(ID)%XzV,        ATM(ID)%DXzV,  ED,T,
	1        ATM(ID)%EQXzV, ATM(ID)%NXzV,NT,ND, ATM(ID)%XzV_PRES)
	  END DO
!
! This routine not only evaluates the LTE populations of both model atoms, but
! it also evaluates the dln(LTE Super level Pop)/dT.
!
	INCLUDE 'EVAL_LTE_INC_V5.INC'
!
! Compute the turbulent velocity and MINIMUM Doppler velocity as a function of depth.
! For the later, the iron mas of 55.8amu is assumed.
!
!
! Bug fixes, not flag-gated (upstream code, 30-Sep-2026):
! LAW_V1 divided by V(1), which is zero for a finite base velocity
! offset, small at a resolved shock, or negative under infall -- any
! of these poisons VTURB_VEC (and hence VDOP_VEC).  Normalize by
! max|V| and use |V|, which reduces to the old law for monotonic
! outflow but stays defined for non-monotonic fields.
!
	IF(TURB_LAW .EQ. 'LAW_V1')THEN
	  T1=MAXVAL(ABS(V(1:ND)))
	  IF(T1 .GT. 1.0E-12_LDP)THEN
	    VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*
	1       MIN(1.0_LDP,ABS(V(1:ND))/T1)
	  ELSE
	    VTURB_VEC(1:ND)=VTURB_MIN
	  END IF
!
! Bug fix, not flag-gated (upstream code): the old loop read
! R(I+1)/ESEC(I+1) at I=ND (out of bounds) and built every TA(I)
! relative to TA(1) instead of accumulating, so the optical-depth
! weighting was wrong at all depths.  TA is now the cumulative
! electron-scattering depth integral inward from the outer boundary;
! the trapezoid assumes the CMFGEN convention R(1) > ... > R(ND).
!
	ELSE IF(TURB_LAW .EQ. 'LAW_TAU1')THEN
	  IF(ANY(R(1:ND-1) .LE. R(2:ND)))THEN
	    WRITE(LUER,*)'Error in CMF_FLUX_SUB_V5: LAW_TAU1 requires'
	    WRITE(LUER,*)'  R decreasing with depth index'
	    STOP
	  END IF
	  ESEC(1:ND)=6.65E-15_LDP*MAX(ED(1:ND),0.0_LDP)*
	1       MAX(CLUMP_FAC(1:ND),0.0_LDP)
	  TA(1)=1.0_LDP
	  DO I=2,ND
	    TA(I)=TA(I-1)+0.5_LDP*(R(I-1)-R(I))*(ESEC(I-1)+ESEC(I))
	  END DO
	  VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)/
	1       (1.0_LDP+MAX(TA(1:ND),0.0_LDP))
	ELSE
	  WRITE(LUER,'(A)')' Error --- TURBULENT velocity law not recognized'
	  WRITE(LUER,'(A)')TRIM(TURB_LAW)
	  STOP
	END IF
	VDOP_VEC(1:ND)=SQRT( VTURB_VEC(1:ND)**2 + 2.96_LDP*T(1:ND) )
!
	TA(1:ND)=ABS( CLUMP_FAC(1:ND)-1.0_LDP )
	T1=MAXVAL(TA)
	DO_CLUMP_MODEL=.FALSE.
	IF(T1 .GT. 1.0E-05_LDP)DO_CLUMP_MODEL=.TRUE.
!
!
! Compute profile frequencies such that for the adopted doppler
! velocity the profile ranges from 5 to -5 doppler widths.
! This section needs to be rewritten if we want the profile to
! vary with depth.
	  FIRST=.TRUE.		!Check cross section at edge is non-zero.
	  NCF=0 		!Initialize number of continuum frequencies.
!
	  DO ID=1,NUM_IONS
	    CALL SET_EDGE_FREQ_V3(ID,OBS,NCF,NCF_MAX,
	1           ATM(ID)%EDGEXzV_F,  ATM(ID)%NXzV_F, ATM(ID)%XzV_PRES,
	1           ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV,
	1           ATM(ID)%N_XzV_PHOT)
	  END DO

	  IF(XRAYS)THEN
	    DO ID=1,NUM_IONS-1
	      CALL SET_X_FREQ(OBS,NCF,NCF_MAX,AT_NO(SPECIES_LNK(ID)),
	1           ATM(ID)%ZXzV, ATM(ID)%XzV_PRES, ATM(ID+1)%XzV_PRES)
	    END DO
	  END IF
!
! Now insert addition points into frequency array. WSCI is used as a
! work array - okay since of length NCF_MAX, and zeroed in QUADSE.
! OBSF contains the bound-free edges - its contents are zero on
! subroutine exit. J is used as temporary variable for the number of
! frequencies transmitted to SET_CONT_FREQ. NCF is returned as the number
! of frequency points. FQW is used a an integer array for the sorting ---
! we know it has the correct length since it is the same size as NU.
! LUIN --- Used as temporary LU (opened and closed).
!
	  J=NCF
!	  CALL SET_CONT_FREQ(NU,OBS,FQW,
!	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
!	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
!	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
!	1                        J,NCF,NCF_MAX,LUIN)
	T1=5000.0_LDP
	  CALL SET_CONT_FREQ_V4(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_CONT,T1,
	1                        J,NCF,NCF_MAX,LUIN)
!
! 
!
! Set up lines that will be treated with the continuum calculation.
! This section of code is also used by the code treating purely lines
! (either single transition Sobolev or CMF, or overlapping Sobolev).
!
! To define the line transitions we need to operate on the FULL atom models.
! We thus perform separate loops for each species. VEV_TRANS_NAME is
! allocated temporaruly so that we can output the full transitions name
! to TRANS_INFO.
!
	ML=0			!Initialize line counter.
	ALLOCATE (VEC_TRANS_NAME(NLINE_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for VEC_TRANS_NAME'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
!
	ESEC(1:ND)=6.65E-15_LDP*ED(1:ND)
!
! The onlye species not present is the ion corresponding
! to the last ioization stage considered.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNL=1, ATM(ID)%NXzV_F-1
	      DO MNUP=MNL+1, ATM(ID)%NXzV_F
	        NL=ATM(ID)%F_TO_S_XzV(MNL)+ATM(ID)%EQXzV-1
	        NUP=ATM(ID)%F_TO_S_XzV(MNUP)+ATM(ID)%EQXzV-1
	        IF(ATM(ID)%AXzV_F(MNL,MNUP) .NE. 0)THEN
	          ML=ML+1
	          IF(ML .GT. NLINE_MAX)THEN
	            WRITE(LUER,*)'NLINE_MAX is too small in CMF_FLUX_SUB_V5'
	            STOP
	          END IF
	          VEC_FREQ(ML)=ATM(ID)%EDGEXzV_F(MNL)-ATM(ID)%EDGEXzV_F(MNUP)
	          VEC_SPEC(ML)=ION_ID(ID)
	          VEC_NL(ML)=NL
	          VEC_NUP(ML)=NUP
	          VEC_MNL_F(ML)=MNL
	          VEC_MNUP_F(ML)=MNUP
	          VEC_OSCIL(ML)=ATM(ID)%AXzV_F(MNL,MNUP)
	          VEC_EINA(ML)=ATM(ID)%AXzV_F(MNUP,MNL)
	          VEC_ARAD(ML)= ATM(ID)%ARAD(MNL)+ATM(ID)%ARAD(MNUP)
	          VEC_C4(ML)= ABS(ATM(ID)%C4(MNL)) + ABS(ATM(ID)%C4(MNUP))
	          VEC_TRANS_TYPE(ML)=ATM(ID)%XzV_TRANS_TYPE
	          VEC_TRANS_NAME(ML)=TRIM(VEC_SPEC(ML))//
	1             '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP))//'-'//
	1             TRIM(ATM(ID)%XzVLEVNAME_F(MNL))//')'
	          T1=VEC_OSCIL(ML)*OPLIN
	          T2=ATM(ID)%GXzV_F(MNL)/ATM(ID)%GXzV_F(MNUP)
	          DO I=1,ND
	            CHIL(I)=ABS(T1*(ATM(ID)%XzV_F(MNL,I)-T2*ATM(ID)%XzV_F(MNUP,I)))
	          END DO
	          PROF_TYPE(ML)=ATM(ID)%XzV_PROF_TYPE
	          IF(GLOBAL_LINE_PROF .NE. 'NONE')PROF_TYPE(ML)=GLOBAL_LINE_PROF
	          T2=0.0_LDP
	          CALL SET_PROF_LIMITS_V4(VEC_STRT_FREQ(ML),VEC_VDOP_MIN(ML),
	1             CHIL,ED,T,VTURB_VEC,ND,PROF_TYPE(ML),PROF_LIST_LOCATION(ML),
	1             VEC_FREQ(ML),MNL,MNUP,
	1             VEC_SPEC(ML),AT_MASS(SPECIES_LNK(ID)), ATM(ID)%ZXzV,
	1             VEC_ARAD(ML),VEC_C4(ML),TDOP,AMASS_DOP,VTURB_FIX,                !2: Garbage at present
	1             DOP_PROF_LIMIT,VOIGT_PROF_LIMIT,V_PROF_LIMIT,MAX_PROF_ED,
	1             SET_PROF_LIMS_BY_OPACITY)
	          IF(SPECIES(SPECIES_LNK(ID)) .EQ. 'HYD')THEN
	            IF(ATM(ID)%XzVLEVNAME_F(2) .NE. '2___' .AND. INDEX(PROF_TYPE(ML),'DOP') .EQ. 0)THEN
	              WRITE(LUER,*)'Error in CMF_FLUX -- only Doppler profiles are implemented for split H levels'
	              STOP
	            END IF
	          END IF
	        END IF
	      END DO
	    END DO
	  END IF
	END DO
	N_LINE_FREQ=ML
	CALL ADJUST_LINE_FREQ_V2(VEC_FREQ,VEC_STRT_FREQ,VEC_TRANS_NAME,N_LINE_FREQ,LUIN)
!
! 
!
! GLOBAL_LINE_SWITCH provides an option to handle all LINE by the same method.
! The local species setting only takes precedence when it is set to NONE.
!
	IF(GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE')THEN
	  DO I=1,N_LINE_FREQ
	    VEC_TRANS_TYPE(I)=GLOBAL_LINE_SWITCH
	  END DO
	ELSE
	  DO I=1,N_LINE_FREQ
	    CALL SET_CASE_UP(VEC_TRANS_TYPE(I),IZERO,IZERO)
	  END DO
	END IF
!
! If desired, we can set transitions with:
!      wavelengths > FLUX_CAL_LAM_END (in A) to the SOBOLEV option.
!      wavelengths < FLUX_CAL_LAM_BEG (in A) to the SOBOLEV option.
!
! The region defined by FLUX_CAL_LAM_BEG < LAM < FLUX_CAL_LAM_END will be computed using
! transition types determined by the earlier species and global options.
!
! Option has 2 uses:
!
! 1. Allows use of SOBOLEV approximation in IR where details of radiative
!    transfer is unimportant. In this case FLUX_CAL_LAM_BEG should be set to zero.
! 2. Allows a full flux calculation to be done in a limited wavelength region
!    as defined by FLUX_CAL_LAM_END and FLUX_CAL_LAM_BEG.
!
	IF(SET_TRANS_TYPE_BY_LAM)THEN
	  IF(FLUX_CAL_LAM_END .LT. FLUX_CAL_LAM_BEG)THEN
	    WRITE(LUER,*)'Error in CMF_FLUX_SUB_V5'
	    WRITE(LUER,*)'FLUX_CAL_LAM_END must be > FLUX_CAL_LAM_BEG'
	    STOP
	  END IF
	  GLOBAL_LINE_SWITCH='NONE'
	  T1=SPEED_OF_LIGHT()*1.0E-07_LDP
	  DO I=1,N_LINE_FREQ
	    IF(T1/VEC_FREQ(I) .GE. FLUX_CAL_LAM_END)THEN
	      VEC_TRANS_TYPE(I)='SOB'
	    END IF
	    IF(T1/VEC_FREQ(I) .LE. FLUX_CAL_LAM_BEG)THEN
	      VEC_TRANS_TYPE(I)='SOB'
	    END IF
	  END DO
	END IF
!
	DO ML=1,N_LINE_FREQ
	  IF(VEC_TRANS_TYPE(ML) .NE. 'BLANK')VEC_STRT_FREQ(ML)=VEC_FREQ(ML)
	END DO
!
! Sort lines into numerically decreaing frequency. This is used for
! outputing TRANS_INFO file. Need to sort all the VECTORS, as they
! are linked.
!
	CALL INDEXX(N_LINE_FREQ,VEC_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_C4,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_VDOP_MIN,VEC_INDX,VEC_DP_WRK)
!
	CALL SORTINT(N_LINE_FREQ,VEC_NL,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NUP,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,PROF_LIST_LOCATION,VEC_INDX,VEC_INT_WRK)
!
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,PROF_TYPE,VEC_INDX,VEC_CHAR_WRK)
!
	IF(WRITE_TRANS_INFO)THEN
	  I=160	!Record length - allow for long names
	  CALL GEN_ASCI_OPEN(LUIN,'TRANS_INFO','UNKNOWN',' ','WRITE',I,IOS)
	    WRITE(LUIN,*)
	1     '     I    NL_F  NUP_F        Nu',
	1     '       Lam(A)    /\V(km/s)    Transition'
	    WRITE(LUIN,
	1      '(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,16X,A)')
	1         IONE,VEC_MNL_F(1),VEC_MNUP_F(1),
	1         VEC_FREQ(1),LAMVACAIR(VEC_FREQ(1)),
	1         TRIM(VEC_TRANS_NAME(VEC_INDX(1)))
	    DO ML=2,N_LINE_FREQ
	      T1=LAMVACAIR(VEC_FREQ(ML))
	      T2=C_KMS*(VEC_FREQ(ML-1)-VEC_FREQ(ML))/VEC_FREQ(ML)
	      IF(T2 .GT. C_KMS)T2=C_KMS
	      IF(T1 .LT. 1.0E+04_LDP)THEN
	        WRITE(LUIN,'(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,2X,F10.2,4X,I7,2X,A12,3X,A)')
	1           ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1           VEC_FREQ(ML),T1,T2,PROF_LIST_LOCATION(ML),PROF_TYPE(ML),TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	      ELSE
	        WRITE(LUIN,'(1X,I6,2(1X,I6),2X,F10.6,2X,1P,E10.4,0P,2X,F10.2,4X,I7,2X,A12,3X,A)')
	1           ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1           VEC_FREQ(ML),T1,T2,PROF_LIST_LOCATION(ML),PROF_TYPE(ML),TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	      END IF
	    END DO
	  CLOSE(UNIT=LUIN)
	END IF
	DEALLOCATE (VEC_TRANS_NAME)
!
! Get lines and arrange in numerically decreasing frequency according to
! the START frequency of the line. This will allow us to consider line overlap,
! and to include lines with continuum frequencies so that the can be handled
! automatically.
!
	CALL INDEXX(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_C4,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_VDOP_MIN,VEC_INDX,VEC_DP_WRK)
!
	CALL SORTINT(N_LINE_FREQ,VEC_NL,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NUP,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,PROF_LIST_LOCATION,VEC_INDX,VEC_INT_WRK)
!
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,PROF_TYPE,VEC_INDX,VEC_CHAR_WRK)
!
!
!
! We have found all lines. If we are doing a blanketing calculation for this
! line we insert them into the continuum frequency set, otherwise the
! line is not included.
!
	DO ML=1,NCF
	  FQW(ML)=NU(ML)	!FQW has temporary storage of continuum freq.
	END DO
	V_DOP=MINVAL(VEC_VDOP_MIN)
	CALL INS_LINE_V5(  NU,LINES_THIS_FREQ,I,NCF_MAX,
	1		  VEC_FREQ,VEC_STRT_FREQ,VEC_VDOP_MIN,VEC_TRANS_TYPE,
	1                 LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,N_LINE_FREQ,
	1                 FQW,NCF,FRAC_DOP,VINF,dV_CMF_PROF,dV_CMF_WING,
	1                 ES_WING_EXT,R_CMF_WING_EXT,L_FALSE )
!
	K=NCF		!# of continuum frequencies: Need for DET_MAIN...
	NCF=I		!Revised
!
	WRITE(LUER,*)' '
	WRITE(LUER,'(A,1X,I7)')' Number of line frequencies is:',N_LINE_FREQ
	WRITE(LUER,'(A,6X,I7)')' Number of continuum frequencies is:',NCF
!
! Used inthis context, every edge must be with V_DOP km/s of a frequency in the NU
! array.
!
	V_DOP=FRAC_DOP*MINVAL(VEC_VDOP_MIN)
	CALL DET_MAIN_CONT_FREQ(NU,NCF,FQW,K,NU_EVAL_CONT,
	1             V_DOP,DELV_CONT,COMPUTE_ALL_CROSS)
!
! Redefine frequency quadrature weights.
!
	CALL SMPTRP(NU,FQW,NCF)
	DO ML=1,NCF
	  FQW(ML)=FQW(ML)*1.0E+15_LDP
	END DO
!
! Set observers frequencies. The slight fiddling in setting NU_MAX and NU_MIN
! is done so that the CMF frequencies encompass all observers frame
! frequencies. This allows computation of all observers fluxes allowing for
! velocity effects.
!
! We always insert lines into the observers frame frequencies. This allows
! a blanketed model spectrum to be divided by an unblanketed model
! spectrum.
!
	CALL TUNE(1,'INS_LINE_OBS')
	  WRITE(LUER,*)'Calling INS_LINE_OBS_V4'
	  T1=NU(1)/(1.0_LDP+2.0_LDP*VINF/C_KMS)
	  NU_MAX_OBS=MIN(T1,NU(3))
	  IF(RD_NU_MAX_OBS .GT. 0.0_LDP)NU_MAX_OBS=RD_NU_MAX_OBS
	  T1=NU(NCF)*(1.0_LDP+2.0_LDP*VINF/C_KMS)
	  NU_MIN_OBS=MAX(NU(NCF-3),T1)
	  IF(RD_NU_MIN_OBS .GT. 0.0_LDP)NU_MIN_OBS=RD_NU_MIN_OBS
	  CALL INS_LINE_OBS_V4(OBS_FREQ,N_OBS,NCF_MAX,
	1               VEC_FREQ,VEC_STRT_FREQ,VEC_VDOP_MIN,VEC_TRANS_TYPE,
	1               N_LINE_FREQ,SOB_FREQ_IN_OBS,
	1		NU_MAX_OBS,NU_MIN_OBS,VINF,
	1               FRAC_DOP_OBS,dV_OBS_PROF,dV_OBS_WING,dV_OBS_BIG,
	1               OBS_PRO_EXT_RAT,ES_WING_EXT,VTURB_MAX)
	  WRITE(LUER,*)'Finished calling INS_LINE_OBS_V4'
	CALL TUNE(2,'INS_LINE_OBS')
!
	CALL TUNE(1,'SET_PROF_STORE')
	WRITE(LUER,*)'Setting pofile storage'
	CALL GET_PROFILE_STORAGE_LIMITS_V2(NLINES_PROF_STORE,NFREQ_PROF_STORE,
	1         NLINES_VOIGT_STORE,NFREQ_VOIGT_STORE,MAX_SIM,
	1         LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,
	1         PROF_TYPE,VEC_TRANS_TYPE,N_LINE_FREQ,NCF)
!
	IF(MAX_SIM .EQ. 0)MAX_SIM=10			!Needs to be finite for EW computation
	ALLOCATE(EINA(MAX_SIM))
	ALLOCATE(OSCIL(MAX_SIM))
	ALLOCATE(GLDGU(MAX_SIM))
	ALLOCATE(AMASS_SIM(MAX_SIM))
	ALLOCATE(FL_SIM(MAX_SIM))
	ALLOCATE(SIM_NL(MAX_SIM))
	ALLOCATE(SIM_NUP(MAX_SIM))
	ALLOCATE(CHIL_MAT(ND,MAX_SIM))
	ALLOCATE(ETAL_MAT(ND,MAX_SIM))
	ALLOCATE(BB_COR(ND,MAX_SIM))
	ALLOCATE(L_STAR_RATIO(ND,MAX_SIM))
	ALLOCATE(U_STAR_RATIO(ND,MAX_SIM))
	ALLOCATE(LINE_PROF_SIM(ND,MAX_SIM))
	ALLOCATE(TRANS_NAME_SIM(MAX_SIM))
	ALLOCATE(RESONANCE_ZONE(MAX_SIM))
	ALLOCATE(END_RES_ZONE(MAX_SIM))
	ALLOCATE(LINE_STORAGE_USED(MAX_SIM))
	ALLOCATE(SIM_LINE_POINTER(MAX_SIM))
!
	CALL INIT_PROF_MODULE_V2(ND,NLINES_PROF_STORE,NFREQ_PROF_STORE,
	1         NLINES_VOIGT_STORE,MAX_SIM)
	WRITE(LUER,*)'Finished set profile storage'
	CALL TUNE(2,'SET_PROF_STORE')
	CALL TUNE(3,' ')
!
! 
! Need to calculate impact parameters, and angular quadrature weights here
! as these may be required when setting up the initial temperature
! distribution of the atmosphere (i.e. required by JGREY).
!
!
! Compute impact parameter values P
!
	CALL IMPAR(P,R,RP,NC,ND,NP)
!
! Compute the angular quadrature weights
!
	IF(TRAPFORJ)THEN
	  CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,JTRPWGT_V2)
	  CALL NORDANGQW_V2(HQW,R,P,NC,ND,NP,HTRPWGT_V2)
	  CALL NORDANGQW_V2(KQW,R,P,NC,ND,NP,KTRPWGT_V2)
	  CALL NORDANGQW_V2(NQW,R,P,NC,ND,NP,NTRPWGT_V2)
	  MID=.TRUE.
	  CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HTRPWGT_V2,MID)
	  CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NTRPWGT_V2,MID)
	ELSE
	  CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,JWEIGHT_V2)
	  CALL NORDANGQW_V2(HQW,R,P,NC,ND,NP,HWEIGHT_V2)
	  CALL NORDANGQW_V2(KQW,R,P,NC,ND,NP,KWEIGHT_V2)
	  CALL NORDANGQW_V2(NQW,R,P,NC,ND,NP,NWEIGHT_V2)
	  MID=.TRUE.
	  CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HWEIGHT_V2,MID)
	  CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NWEIGHT_V2,MID)
	END IF
!
	IF(ACCURATE)THEN
!
! We first verify that the interpolation range is valid.
!
	  IF(END_INTERP_INDX .GT. ND)END_INTERP_INDX=ND
	  IF(DEEP .GT. ND)DEEP=MIN(5,ND)
	  NDEXT=(END_INTERP_INDX-ST_INTERP_INDX)*NPINS+ND
	  IF(NDEXT .GT. NDMAX)THEN
	    WRITE(LUER,*)' Error - NDEXT larger than NDMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NDMAX in CMFGEN'
	    STOP
	  END IF
	  NCEXT=NC*(NPINS+1)
!
! NB: The following expression guarentees that NPEXT has the same relationship
! to NDEXT and NCEXT as does NP to ND and NC.
!
	  NPEXT=NDEXT+NCEXT+(NP-ND-NC)
	  IF(NPEXT .GT. NPMAX)THEN
	    WRITE(LUER,*)' Error - NPEXT larger than NPMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NPMAX in CMFGEN'
	    STOP
	  END IF
	  I=ND-DEEP
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1         ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1         V,T,SIGMA,ND)
	  CALL IMPAR(PEXT,REXT,RP,NCEXT,NDEXT,NPEXT)
	  LANG_COORDEXT=0.0_LDP
!
	  ALLOCATE (AQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (KQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	    WRITE(LUER,*)'Unable to allocate memory for AQWEXT'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
!	
! Compute the turbulent velocity and MINIMUM Doppler velocity as a function of depth.
! For the later, the iron mass of 55.8amu is assumed.
!
	  TA(1:NDEXT)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*VEXT(1:NDEXT)/VEXT(1)
	  VDOP_VEC_EXT(1:NDEXT)=SQRT( TA(1:NDEXT)**2 + 2.96_LDP*TEXT(1:NDEXT) )
!
! Note that the F2DAEXT vectors (here used as dummy variables) must be at least
! NPEXT long.
!
	  IF(TRAPFORJ)THEN
	    CALL NORDANGQW_V2(AQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,JTRPWGT_V2)
	    CALL NORDANGQW_V2(HQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HTRPWGT_V2)
	    CALL NORDANGQW_V2(KQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,KTRPWGT_V2)
	    CALL NORDANGQW_V2(NQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NTRPWGT_V2)
	    MID=.TRUE.
	    CALL GENANGQW_V2(HMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HTRPWGT_V2,MID)
	    CALL GENANGQW_V2(NMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NTRPWGT_V2,MID)
	  ELSE
	    CALL NORDANGQW(AQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,JWEIGHT_V2)
	    CALL NORDANGQW(HQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HWEIGHT_V2)
	    CALL NORDANGQW(KQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,KWEIGHT_V2)
	    CALL NORDANGQW(NQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NWEIGHT_V2)
	    MID=.TRUE.
	    CALL GENANGQW_V2(HMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HWEIGHT_V2,MID)
	    CALL GENANGQW_V2(NMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NWEIGHT_V2,MID)
	  END IF
	ELSE
	  NDEXT=ND ; NCEXT=NC; NPEXT=NP
	  TEXT(1:ND)=T(1:ND)
	END IF
        CALL SET_POP_FOR_TWOJ(POS_IN_NEW_GRID,EDD_CONT_REC,LU_EDD,NDEXT)
!
! Allocate arrays and vectors for computing observed fluxes.
!
	IF(ACCURATE)THEN
	  NP_OBS_MAX=NPEXT+12
	ELSE
	  NP_OBS_MAX=NP+12
	END IF
	ALLOCATE (IPLUS_STORE(NST_CMF,NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (P_OBS(NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (IPLUS(NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (MU_AT_RMAX(NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (dMU_AT_RMAX(NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (HQW_AT_RMAX(NP_OBS_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for IPLUS_STORE'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
!
! Used when computing the observed fluxes. These will get overwritten
! if we do an accurate comoving frame soluton using CMF_FORM_SOL.
!
! MAX() on the radicand guards the sqrt against a P exceeding R by a
! few ulps after P-grid interpolation or revision (tangential ray).
!
	IF(ACCURATE)THEN
	  DO LS=1,NPEXT
	    MU_AT_RMAX(LS)=SQRT( MAX(0.0_LDP,
	1        1.0_LDP -(PEXT(LS)/REXT(1))**2 ) )
	    HQW_AT_RMAX(LS)=HQWEXT(1,LS)
	  END DO
	  IF(PEXT(NPEXT) .EQ. REXT(1))MU_AT_RMAX(NPEXT)=0.0_LDP
	ELSE
	  DO LS=1,NP
	    MU_AT_RMAX(LS)=SQRT( MAX(0.0_LDP,
	1        1.0_LDP -(P(LS)/R(1))**2 ) )
	    HQW_AT_RMAX(LS)=HQW(1,LS)
	  END DO
	  IF(P(NP) .EQ. R(1))MU_AT_RMAX(NP)=0.0_LDP
	END IF
!
	ALLOCATE (ETA_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CHI_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RJ_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0 .AND. INCL_RAY_SCAT)ALLOCATE(RAY_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0 .AND. INCL_DUST)ALLOCATE(DUST_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for IPLUS_STORE'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
! 
!
! Set 2-photon data with current atomic models and populations.
!
	DO ID=1,NUM_IONS
	  ID_SAV=ID
	  CALL SET_TWO_PHOT_V3(TRIM(ION_ID(ID)),ID_SAV,
	1       ATM(ID)%XzVLTE,         ATM(ID)%NXzV,
	1       ATM(ID)%XzVLTE_F_ON_S,  ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,      ATM(ID)%GXzV_F,
	1       ATM(ID)%F_TO_S_XzV,     ATM(ID)%NXzV_F,     ND,
	1       ATM(ID)%ZXzV,           ATM(ID)%EQXzV,      ATM(ID)%XzV_PRES)
	END DO
!
	DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND))
! 
!
! Decide whether to use an file with old J values to provide an initial
! estimate of J with incoherent electron scattering. The options
!
!             .NOT. RD_COHERENT_ES .AND. COHERENT_ES
!
! indicate that no ES_J_CONV file exists already.
! Note that TEXT and NDEXT contain T and ND when ACCURATE is FALSE.
!
	IF(USE_OLDJ_FOR_ES .AND. .NOT. RD_COHERENT_ES .AND. COHERENT_ES)THEN
	  COHERENT_ES=RD_COHERENT_ES
	  I=ND*NCF
	  CALL COMP_J_CONV_V2(ETA_CMF_ST,I,NU,TEXT,NDEXT,NCF,LUIN,
	1        'OLD_J_FILE',
	1        EDD_CONT_REC,L_FALSE,L_TRUE,LU_ES,'ES_J_CONV')
	END IF
!
! Electron scattering iteration loop. We perform this loop to
! take into account incoherent electron scattering.
!
	IF(RD_COHERENT_ES)NUM_ES_ITERATIONS=1
	DO ES_COUNTER=1,NUM_ES_ITERATIONS
!
! This forces CMF_FORM_SOL_V? to be used for OBSFLUX computation.
!
	  IF(ES_COUNTER .EQ. NUM_ES_ITERATIONS)LST_ITERATION=.TRUE.
!
! NB: NDEXT contains ND if ACCURATE is false. The +1 arrises since
! we also write NU(ML) out.
!
	  IF(ACCURATE .OR. EDD_CONT .OR. EDD_LINECONT)THEN
	    INQUIRE(UNIT=LU_EDD,OPENED=TMP_LOG)
	    IF(.NOT. TMP_LOG)THEN
	      CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1       REXT,VEXT,LANG_COORDEXT,NDEXT,
	1       ACCESS_F,L_TRUE,COMPUTE_EDDFAC,USE_FIXED_J,'EDDFACTOR',LU_EDD)
	    END IF
	  END IF
!
! Open the file so it can be read in the CONTINUUM loop (read in  COMP_JCONT).
!
	  IF(.NOT. COHERENT_ES)THEN
	    I=WORD_SIZE*(ND+1)/UNIT_SIZE
	    IF(ACCURATE)I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
	    OPEN(UNIT=LU_ES,FILE='ES_J_CONV',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error opening ES_J_CONV - will compute new J'
	        COHERENT_ES=.TRUE.
	      END IF
	  END IF
!
!
!***************************************************************************
!***************************************************************************
!
!                         CONTINUUM LOOP
!
!***************************************************************************
!***************************************************************************
!
	EDDINGTON=EDD_CONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  IF(COMPUTE_EDDFAC)THEN
	    WRITE(LU_EDD,REC=EDD_CONT_REC)ACCESS_F,NCF,NDEXT
	  ELSE
	    READ(LU_EDD,REC=EDD_CONT_REC)ACCESS_F
	  END IF
	END IF
!
	FIRST_OBS_COMP=.TRUE.
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the BLANKETING portion of the code.
!
	LAST_LINE=0			!Updated as each line is done
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(MIN(LAST_LINE+1,N_LINE_FREQ))(1:4) .NE. 'BLAN')
	        LAST_LINE=LAST_LINE+1
	END DO
	DO SIM_INDX=1,MAX_SIM
	  LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	END DO
!
	CONT_FREQ=0.0_LDP
!
! Define parameters to allow the Cummulative force multipler to be output at
! a function of frequency. We presently output the force multiplier every 500km/s.
!
	N_FORCE=LOG(NU(NCF)/NU(1))/LOG(1.0_LDP-500.0_LDP/C_KMS)
	NU_FORCE=NU(1)
	NU_FORCE_FAC=(1.0_LDP-500.0_LDP/C_KMS)
	ML_FORCE=1
	FLUSH(LUER)
!
! Enter loop for each continuum frequency.
!
	CALL TUNE(IONE,'MLCF')
!
! DO_NONMONO_KNOP (Knop et al. 2009, Papers I/II): register the OBS
! frequency grid with the formal solver before the first call below.
! The FGINEXT_GS checkpoint was written on the model frequency grid
! (different NCF and spacing), so the solver init remaps the lagged
! I_NEXT coupling in frequency instead of restoring by position or
! cold-starting it.
!
	IF(DO_NONMONO_KNOP)CALL SET_FG_GS_OBS_GRID(NU,NCF)
!
	DO 10000 ML=1,NCF
	  FREQ_INDX=ML
	  FL=NU(ML)
	  IF(ML .EQ. 1)THEN
	    FIRST_FREQ=.TRUE.
	  ELSE
	    FIRST_FREQ=.FALSE.
	  END IF
	  SECTION='CONTINUUM'
!
	  IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	    COMPUTE_NEW_CROSS=.TRUE.
	    CONT_FREQ=NU_EVAL_CONT(ML)
	  ELSE
	    COMPUTE_NEW_CROSS=.FALSE.
	  END IF
!
! 
!
! Section to include lines automatically with the continuum.
!
!
!  LINES_THIS_FREQ --- Logical vector [NCF] indicating whether this frequency
!                        is part of the resonance zone (i.e. Doppler profile) of
!                        one (or more) lines.
!
! LINE_ST_INDX_IN_NU --- Integer vector [N_LINES] which specifies the starting
!                          frequency index for this lines resonance zone.
!
! LINE_END_INDX_IN_NU --- Integer vector [N_LINES] which specifies the final
!                          frequency index for this lines resonance zone.
!
! FIRST_LINE   ---- Integer specifying the index of the highest frequency
!                         line which we are taking into account in the
!                         transfer.
!
! LAST_LINE  ---- Integer specifying the index of the lowest frequency
!                         line which we are taking into account in the
!                         transfer.
!
! LINE_LOC   ---- Integer array. Used to locate location of a particular line
!                         in the SIM vectors/arrays.
!
! SIM_LINE_POINTER --- Integer array --- locates the line corresponding to
!                         the indicated storage location in the SIM vectors/
!                         arrays.
!
! Check whether we have to treat another line. We use a DO WHILE, rather
! than an IF statement, to handle lines which begin at the same (upper)
! frequency.
!
!
	CALL TUNE(IONE,'ADD_LINE')
	DO WHILE( LAST_LINE .LT. N_LINE_FREQ .AND.
	1                ML .EQ. LINE_ST_INDX_IN_NU(MIN(LAST_LINE+1,N_LINE_FREQ)) )
!
! Have another line --- need to find its storage location.
!
	  I=1
	  DO WHILE(LINE_STORAGE_USED(I))
	    I=I+1
	    IF(I .GT. MAX_SIM)THEN
	      FIRST_LINE=N_LINE_FREQ
	      DO SIM_INDX=1,MAX_SIM		!Not 0 as used!
	        FIRST_LINE=MIN(FIRST_LINE,SIM_LINE_POINTER(SIM_INDX))
	      END DO
	      IF( ML .GT. LINE_END_INDX_IN_NU(FIRST_LINE))THEN
!
! Free up storage location for line.
!
	        I=LINE_LOC(FIRST_LINE)
	        LINE_STORAGE_USED(I)=.FALSE.
	        SIM_LINE_POINTER(I)=0
	      ELSE
	        WRITE(LUER,*)'Too many lines have overlapping '//
	1                      'resonance zones'
	        WRITE(LUER,*)'Current frequency is:',FL
	        STOP
	      END IF
	    END IF
	  END DO
!
	  SIM_INDX=I
	  LAST_LINE=LAST_LINE+1
	  LINE_STORAGE_USED(SIM_INDX)=.TRUE.
	  LINE_LOC(LAST_LINE)=SIM_INDX
	  SIM_LINE_POINTER(SIM_INDX)=LAST_LINE
!
! Have located a storage location. Now must compute all relevant quantities
! necessary to include this line in the transfer calculations.
!
	  SIM_NL(SIM_INDX)=VEC_NL(LAST_LINE)
	  SIM_NUP(SIM_INDX)=VEC_NUP(LAST_LINE)
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
!
	  EINA(SIM_INDX)=VEC_EINA(LAST_LINE)
	  OSCIL(SIM_INDX)=VEC_OSCIL(LAST_LINE)
	  FL_SIM(SIM_INDX)=VEC_FREQ(LAST_LINE)
!
	  TRANS_NAME_SIM(SIM_INDX)=TRIM(VEC_SPEC(LAST_LINE))
	  AMASS_SIM(SIM_INDX)=AMASS_ALL(NL)
!
! 
!
! Compute U_STAR_RATIO and L_STAR_RATIO which are used to switch from
! the opacity/emissivity computed with a FULL_ATOM to an equivalent form
! but written in terms of the SUPER-LEVELS.
!
! L refers to the lower level of the transition.
! U refers to the upper level of the transition.
!
! At present we must treat each species separately (for those with both FULL
! and SUPER_LEVEL model atoms).
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
! MNL (MNUP) denotes the lower (upper) level in the super level model atom.
!
	MNL_F=VEC_MNL_F(LAST_LINE)
	MNUP_F=VEC_MNUP_F(LAST_LINE)
	DO K=1,ND
	  L_STAR_RATIO(K,SIM_INDX)=1.0_LDP
	  U_STAR_RATIO(K,SIM_INDX)=1.0_LDP
	END DO
!
! T1 is used to represent b(level)/b(super level). If no interpolation of
! the b values in a super level has been performed, this ratio will be unity .
! This ratio is NOT treated in the linearization.
!
	DO ID=1,NUM_IONS
	  IF(VEC_SPEC(LAST_LINE) .EQ. ION_ID(ID))THEN
	    MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	    MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	    DO K=1,ND
              T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzV(MNL,K))/ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
              L_STAR_RATIO(K,SIM_INDX)=T1*(ATM(ID)%W_XzV_F(MNUP_F,K)/ATM(ID)%W_XzV_F(MNL_F,K))*ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
              T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzV(MNUP,K))/ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	      U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	    END DO
	    GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	    TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1      '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1           TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	    EXIT
	  END IF
	END DO
! 
!
! Compute line opacity and emissivity for this line.
!
	  T1=OSCIL(SIM_INDX)*OPLIN
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=1,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    IF(CHIL_MAT(I,SIM_INDX) .EQ. 0.0_LDP)THEN
	      CHIL_MAT(I,SIM_INDX)=0.01_LDP*T1*POPS(NL,I)*L_STAR_RATIO(I,SIM_INDX)
	      WRITE(LUER,*)'Zero line opacity in CMFGEN_SUB'
	      WRITE(LUER,*)'This needs to be fixed'
	      J=ICHRLEN(TRANS_NAME_SIM(SIM_INDX))
	      WRITE(LUER,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)(1:J)
	    END IF
	  END DO
!
! Ensure that LAST_LINE points to the next LINE that is going to be handled
! in the BLANKETING portion of the code.
!
	  DO WHILE(LAST_LINE .LT. N_LINE_FREQ.AND.
	1            VEC_TRANS_TYPE(MIN(LAST_LINE+1,N_LINE_FREQ))(1:4) .NE. 'BLAN')
	       LAST_LINE=LAST_LINE+1
	  END DO
!	
	END DO	!Checking whether a  new line is being added.
	CALL TUNE(ITWO,'ADD_LINE')
!
! 
!
! Check whether current frequency is a resonance frequency for each line.
!
	DO SIM_INDX=1,MAX_SIM
	  RESONANCE_ZONE(SIM_INDX)=.FALSE.
	  END_RES_ZONE(SIM_INDX)=.FALSE.
	  IF(LINE_STORAGE_USED(SIM_INDX))THEN
	    L=SIM_LINE_POINTER(SIM_INDX)
	    IF( FREQ_INDX .GE. LINE_ST_INDX_IN_NU(L) .AND.
	1          FREQ_INDX .LT. LINE_END_INDX_IN_NU(L))THEN
	      RESONANCE_ZONE(SIM_INDX)=.TRUE.
	    ELSE IF(FREQ_INDX .EQ. LINE_END_INDX_IN_NU(L))THEN
 	      RESONANCE_ZONE(SIM_INDX)=.TRUE.
	      END_RES_ZONE(SIM_INDX)=.TRUE.
	    END IF
	  END IF
	END DO
!
	IF(ES_COUNTER .GT. 1)THEN
	  CALL ESOPAC(ESEC,ED,ND)		!Electron scattering emission factor.
	  CHI_SCAT=ESEC
	  IF(INCL_RAY_SCAT)THEN
	    CHI_RAY(1:ND)=RAY_CMF_ST(1:ND,ML)
	    CHI_SCAT=CHI_SCAT+CHI_RAY
	  END IF
	  IF(INCL_DUST)THEN
	    CHI_SCAT_DUST(1:ND)=DUST_CMF_ST(1:ND,ML)
	    CHI_SCAT=CHI_SCAT+CHI_SCAT_DUST
	  END IF
	  CHI(1:ND)=CHI_CMF_ST(1:ND,ML)
	  ETA(1:ND)=ETA_CMF_ST(1:ND,ML)
	  EMHNUKT(1:ND)=EXP(-HDKT*FL/T(1:ND))
	ELSE
!
! Compute profile: Doppler or Stark: T2 and T3 are presently garbage.
! TB is the proton density.
! TC is the He+ density.
!
	  CALL TUNE(IONE,'SET_PROF')
	  TB(1:ND)=0.0_LDP; TC(1:ND)=0.0_LDP
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. ION_ID(ID) .EQ. 'HI')THEN
	      TB(1:ND)=ATM(ID)%DxzV(1:ND)
	    ELSE IF(ATM(ID)%XzV_PRES .AND. ION_ID(ID) .EQ. 'HeI')THEN
	      TC(1:ND)=ATM(ID)%DxzV(1:ND)
	    END IF
	  END DO
!
	  I=ML
          CALL SET_PROF_V6(LINE_PROF_SIM,AMASS_SIM,
	1               RESONANCE_ZONE,END_RES_ZONE,
	1               SIM_LINE_POINTER,MAX_SIM,Z_POP,I,
	1               ED,TB,TC,T,VTURB_VEC,ND,NT,
	1               TDOP,AMASS_DOP,VTURB_FIX,MAX_PROF_ED,
	1               NORM_PROFILE,LUIN)
	  CALL TUNE(ITWO,'SET_PROF')
!
! 
!
! Determine which method will be used to compute continuum intensity.
!
	  IF(ACCURATE .AND. ALL_FREQ)THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE
	    THIS_FREQ_EXT=.FALSE.
	  END IF
!
! Compute opacity and emissivity.
!
	  CALL TUNE(IONE,'C_OPAC')
	  INCLUDE 'OPACITIES_V5.INC'
	  CALL TUNE(ITWO,'C_OPAC')
!
! Since resonance zones included, we must add the line opacity and
! emissivity to the raw continuum values. We first save the pure continuum
! opacity and emissivity. These are used in carrying the variation of J from
! one frequency to the next.
!
	  DO I=1,ND
	    CHI_CONT(I)=CHI(I)
	    ETA_CONT(I)=ETA(I)
	  END DO
!
	  CALL TUNE(IONE,'OP_TOT')
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      DO I=1,ND
	        CHI(I)=CHI(I)+CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	        ETA(I)=ETA(I)+ETAL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	      END DO
	    END IF
	  END DO
	  CALL TUNE(ITWO,'OP_TOT')
!
! CHECK for negative line opacities. We do not distinguish between lines.
!
	  AT_LEAST_ONE_NEG_OPAC=.FALSE.
	  NEG_OPACITY(1:ND)=.FALSE.
	  IF(NEG_OPAC_OPTION .EQ. 'SRCE_CHK')THEN
	    DO I=1,ND
!
! Bug fix, not flag-gated (upstream code): the SRCE_CHK target
! divides by ETA_CONT, which is zero or subnormal at frequencies
! with negligible continuum emission, producing Inf/NaN CHI that
! then poisons ZETA=ETA/CHI.  The target now falls back to the
! scattering floor, which itself takes a 1E-30 floor so the opacity
! cannot reach zero in electron-free zones.
!
	      T1=0.1_LDP*MAX(ESEC(I),CHI_SCAT(I),1.0E-30_LDP)
	      IF(ABS(ETA_CONT(I)) .GT. 1.0E-30_LDP)THEN
	        T2=0.1_LDP*ETA(I)*MAX(CHI_CONT(I)-ESEC(I),0.0_LDP)/
	1           ETA_CONT(I)
	      ELSE
	        T2=T1
	      END IF
	      IF(CHI(I) .LT. 0.1_LDP*CHI_CONT(I) .AND. CHI(I) .LT. T2)THEN
	        CHI(I)=T2
	        NEG_OPACITY(I)=.TRUE.
	        AT_LEAST_ONE_NEG_OPAC=.TRUE.
	        IF(CHI(I) .LE. 0.0_LDP)THEN
	          WRITE(6,*)'Possible error -- CHI still negative after correction'
	          WRITE(6,'(A,I4,5X,A,ES16.10)')'Depth=',I,'Freq=',FL
	          WRITE(6,'(5(8X,A))')'   ETA','ETA_CONT','    CHI','CHI_CONT','   ESEC'
	          WRITE(6,'(5ES16.10)')ETA(I),ETA_CONT(I),CHI(I),CHI_CONT(I),ESEC(I)
	          CHI(I)=T1
	         END IF
	      ELSE IF(CHI(I) .LT. T1)THEN
	        CHI(I)=T1
	        NEG_OPACITY(I)=.TRUE.
	        AT_LEAST_ONE_NEG_OPAC=.TRUE.
	      END IF
	    END DO
	  ELSE IF(NEG_OPAC_OPTION .EQ. 'ESEC_CHK')THEN
	    DO I=1,ND
	      T1=0.1_LDP*MAX(ESEC(I),CHI_SCAT(I),1.0E-30_LDP)
	      IF(CHI(I) .LT. T1)THEN
	        CHI(I)=T1
	        NEG_OPACITY(I)=.TRUE.
	        AT_LEAST_ONE_NEG_OPAC=.TRUE.
	      END IF
	    END DO
	  END IF
	  IF(LST_ITERATION .AND. AT_LEAST_ONE_NEG_OPAC)THEN
	    WRITE(LU_NEG,'(A,1P,E14.6)')
	1      ' Neg opacity for transition for frequency ',FL
	    DO SIM_INDX=1,MAX_SIM
	      IF(RESONANCE_ZONE(SIM_INDX))THEN
	        WRITE(LU_NEG,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)
	      END IF
	    END DO
	    J=0
	    K=0
	    DO I=1,ND
	     IF(NEG_OPACITY(I) .AND. K .EQ. 0)K=I
	     IF(NEG_OPACITY(I))J=I
	    END DO
	    WRITE(LU_NEG,'(A,2X,I3,5X,A,2XI3)')
	1        ' 1st depth',K,'Last depth',J
	  END IF
	END IF
!
!
! Bug fix, not flag-gated: a NEG_OPAC_OPTION of neither SRCE_CHK nor
! ESEC_CHK, or a zero scattering floor, could leave CHI <= 0 here, so
! ZETA/THETA took 0/0 or x/0.  Divide by a floored CHI instead; the
! floor tracks the continuum opacity scale so transparent zones get
! a large-but-finite source term rather than Inf.
!
	  CHI_FLR=MAX(1.0E-30_LDP,1.0E-12_LDP*MAXVAL(ABS(CHI_CONT(1:ND))))
	  DO I=1,ND
	    T2=MAX(CHI(I),CHI_FLR)
	    ZETA(I)=ETA(I)/T2
	    THETA(I)=CHI_SCAT(I)/T2
	  END DO
!
	  IF(LST_ITERATION .AND. ML .NE. NCF)THEN
	    T1=MAX( LOG(DENSITY(5)/DENSITY(1))/LOG(R(1)/R(5))-1.0_LDP,1.0_LDP )
	    DO I=1,N_TAU_EDGE
	      IF(NU(ML) .GE. TAU_EDGE(I) .AND. NU(ML+1) .LT. TAU_EDGE(I))THEN
	        IF(I .EQ. 1)WRITE(LUER,'(A)')' '
	        WRITE(LUER,'(A,1P,E10.4,A,E10.3)')' Tau(Nu=',NU(ML),
	1          ') at outer boundary is:',CHI_CONT(1)*R(1)/T1
	      END IF
	    END DO
	  END IF
!
! Compute continuum intensity.
!
	  IF(COMPUTE_J)THEN
	    CALL TUNE(IONE,'COMP_JCONT')
	    INCLUDE 'COMP_JCONT_V4.INC'	
	    CALL TUNE(ITWO,'COMP_JCONT')
	  END IF
!
!
! Free up LINE storage locations. As we are only computing the line flux,
! and not its variation, we can free up the memory space as soon as we
! exit the resonance zone. This procedure is much simpler than in CMFGEN,
!
	DO SIM_INDX=1,MAX_SIM
	  IF(END_RES_ZONE(SIM_INDX))THEN
	    LINE_LOC( SIM_LINE_POINTER(SIM_INDX) )=0
	    LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	    SIM_LINE_POINTER(SIM_INDX)=0
	  END IF
	END DO
!
!
! Compute flux distribution and luminosity (in L(sun)) of star. NB: For
! NORDFLUX we always assume coherent scattering.
!
!
	CALL TUNE(IONE,'FLUX_DIST')
	IF(.NOT. COMPUTE_J)THEN
	ELSE IF(THIS_FREQ_EXT .AND. .NOT. CONT_VEL)THEN
!
! Since ETAEXT is not required any more, it will be used
! flux.
!
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCEEXT,ZETAEXT,THETAEXT,RJEXT,NDEXT)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,
	1               SOURCEEXT,CHIEXT,dCHIdR,HQWEXT,ETAEXT,
	1               S1,THK_CONT,DIF,DBB,IC,NCEXT,NDEXT,NPEXT,METHOD)
	  CALL UNGRID(SOB,ND,ETAEXT,NDEXT,POS_IN_NEW_GRID)
	  SOB(2)=ETAEXT(2)				!Special case
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  OBS_FLUX(ML)=6.599341_LDP*SOB(1)*2.0_LDP		!2 DUE TO 0.5U
	ELSE IF(CONT_VEL)THEN
!
! TA is a work vector. TB initially used for extended SOB.
!
	   IF(ACCURATE)THEN
	     CALL REGRID_H(TB,REXT,RSQHNU,HFLUX_AT_OB,HFLUX_AT_IB,NDEXT,TA)
	     DO I=1,ND
	       SOB(I)=TB(POS_IN_NEW_GRID(I))
	     END DO
	   ELSE
	     CALL REGRID_H(SOB,R,RSQHNU,HFLUX_AT_OB,HFLUX_AT_IB,ND,TA)
	   END IF
	   IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V)THEN
	     SOB(1)=HFLUX_AT_OB; SOB(ND)=HFLUX_AT_IB
	     SOB(1:ND)=SOB(1:ND)*R(ND)*R(ND)
	   END IF
!
! Since the flux has been interpolated back onto the original size grid,
! we pass ND rather than NDEXT. The L_TRUE indicate NEW_MOD and
! that we will open a new file. J is used for ACCESS_F.
!
	   IF(WRITE_FLUX .AND. ES_COUNTER .EQ. NUM_ES_ITERATIONS)THEN
	     IF(ML .EQ. 1)THEN
	        CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1              REXT,VEXT,LANG_COORDEXT,ND,
	1              J,L_TRUE,L_TRUE,L_FALSE,'FLUX_FILE',LU_FLUX)
	     END IF
	     TA(1:ND)=13.1986_LDP*SOB(1:ND)
	     WRITE(LU_FLUX,REC=EDD_CONT_REC)INITIAL_ACCESS_REC,NCF,ND
	     WRITE(LU_FLUX,REC=INITIAL_ACCESS_REC+ML-1)(SOB(I),I=1,ND),FL
	   END IF
!
	   CALL TUNE(1,'COMP_OBS_CMF')
	   CALL COMP_OBS_V3(IPLUS,FL,
	1           IPLUS_STORE,NU_STORE,NST_CMF,
	1           P,MU_AT_RMAX,HQW_AT_RMAX,OBS_FREQ,OBS_FLUX,N_OBS,
	1           V_AT_RMAX,RMAX_OBS,'IPLUS','LIN_INT',DO_CMF_REL_OBS,
	1           FIRST_OBS_COMP,WRITE_CMF_IP,NP_OBS)
	   CALL TUNE(2,'COMP_OBS_CMF')
!
	ELSE
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCE,ZETA,THETA,RJ,ND)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,THETA,HQW,SOB,
	1               S1,THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  OBS_FLUX(ML)=6.599341_LDP*SOB(1)*2.0_LDP		!2 DUE TO 0.5U
	END IF
!
! The quantity output is the LINE FORCE multiplier M(t).
!
	IF(FL .LT. NU_FORCE .AND. LST_ITERATION .AND. WRITE_CMF_FORCE)THEN
	  INQUIRE(UNIT=82,OPENED=TMP_LOG)
	  J=82
	  IF(.NOT. TMP_LOG)THEN
	    CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1       REXT,VEXT,LANG_COORDEXT,ND,
	1       K,L_TRUE,L_TRUE,L_FALSE,'CMF_FORCE_DATA',J)
	    WRITE(82,REC=EDD_CONT_REC)K,N_FORCE,ND
	  END IF
	  K=INITIAL_ACCESS_REC
	  WRITE(82,REC=K-1+ML_FORCE)( LINE_FLUXMEAN(I)/STARS_LUM/ESEC(I) ,I=1,ND),NU_FORCE
	  ML_FORCE=ML_FORCE+1
	  NU_FORCE=NU_FORCE*NU_FORCE_FAC
	END IF
!
! The L_TRUE indicate NEW_MOD and that we will open a new file.
! J is used for ACCESS_F.
!
	   IF(WRITE_FLUX .AND. ES_COUNTER .EQ. NUM_ES_ITERATIONS)THEN
	     CALL TUNE(1,'WR_FLUX')
	     IF(ML .EQ. 1)THEN
	        CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1              REXT,VEXT,LANG_COORDEXT,ND,
	1              J,L_TRUE,L_TRUE,L_FALSE,'FLUX_FILE',LU_FLUX)
	     END IF
	     TA(1:ND)=13.1986_LDP*SOB(1:ND)
	     WRITE(LU_FLUX,REC=EDD_CONT_REC)INITIAL_ACCESS_REC,NCF,ND
	     WRITE(LU_FLUX,REC=INITIAL_ACCESS_REC+ML-1)(SOB(I),I=1,ND),FL
	     CALL TUNE(2,'WR_FLUX')
	   END IF
!
! Compute the luminosity, the FLUX mean opacity, and the ROSSELAND
! mean opacities.
!
	IF(ML .EQ. 1)THEN		!Need to move to main loop imit.
	  DO I=1,ND
	    RLUMST(I)=0.0_LDP
	    J_INT(I)=0.0_LDP
	    K_INT(I)=0.0_LDP
	    FLUXMEAN(I)=0.0_LDP
	    LINE_FLUXMEAN(I)=0.0_LDP
	    ROSSMEAN(I)=0.0_LDP
	    PLANCKMEAN(I)=0.0_LDP
	    INT_dBdT(I)=0.0_LDP
	  END DO
	END IF
!
	IF(COMPUTE_J)THEN
	  T1=TWOHCSQ*HDKT*FQW(ML)*(NU(ML)**4)
!
! dev7, not flag-gated: for x=HDKT*FL/T << 1 (IR/radio tail) the
! literal 1-EMHNUKT loses almost all digits to cancellation; the
! series x*(1 - x/2 + x^2/6) is accurate to ~1E-13 relative at the
! 1E-4 switchover.  The 1/CHI in the Rosseland weight is floored so
! a repaired-but-tiny opacity cannot inject Inf into the mean.
!
	  DO I=1,ND		              !(4*PI)**2*Dex(+20)/L(sun)
	    T2=SOB(I)*FQW(ML)*4.1274E-12_LDP
	    RLUMST(I)=RLUMST(I)+T2
	    J_INT(I)=J_INT(I)+RJ(I)*FQW(ML)*4.1274E-12_LDP
	    K_INT(I)=K_INT(I)+K_MOM(I)*FQW(ML)*4.1274E-12_LDP
	    FLUXMEAN(I)=FLUXMEAN(I)+T2*CHI(I)
	    LINE_FLUXMEAN(I)=LINE_FLUXMEAN(I)+T2*(CHI(I)-CHI_CONT(I))
	    T3=HDKT*FL/T(I)
	    T4=1.0_LDP-EMHNUKT(I)
	    IF(T3 .LT. 1.0E-4_LDP)
	1      T4=T3*(1.0_LDP-0.5_LDP*T3*(1.0_LDP-T3/3.0_LDP))
	    T2=T1*EMHNUKT(I)/(  ( T4*T(I) )**2  )
	    INT_dBdT(I)=INT_dBdT(I)+T2
	    ROSSMEAN(I)=ROSSMEAN(I)+T2/MAX(CHI(I),CHI_FLR)
	  END DO
	END IF
!
	T1=TWOHCSQ*FQW(ML)*(NU(ML)**3)
	DO I=1,ND
	  T3=HDKT*FL/T(I)
	  T4=1.0_LDP-EMHNUKT(I)
	  IF(T3 .LT. 1.0E-4_LDP)
	1    T4=T3*(1.0_LDP-0.5_LDP*T3*(1.0_LDP-T3/3.0_LDP))
	  T2=T1*EMHNUKT(I)/T4
	  PLANCKMEAN(I)=PLANCKMEAN(I)+T2*CHI(I)
	END DO
!
	CALL TUNE(ITWO,'FLUX_DIST')
!
! Compute and output line force contributed by each species.
!
	IF(WR_ION_LINE_FORCE)THEN
	  IF(.NOT. ALLOCATED(ION_LINE_FORCE))ALLOCATE (ION_LINE_FORCE(ND,NUM_IONS))
	  IF(ML .EQ. 1)ION_LINE_FORCE=0.0_LDP
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      K=SIM_LINE_POINTER(SIM_INDX)
	      DO ID=1,NUM_IONS
	        IF(VEC_SPEC(K) .EQ. ION_ID(ID))THEN
	          DO I=1,ND
	            ION_LINE_FORCE(I,ID)=ION_LINE_FORCE(I,ID)+4.1274E-12_LDP*FQW(ML)*SOB(I)*
	1                CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	          END DO
	          EXIT
	        END IF
	      END DO
	    END IF
	  END DO
	END IF
!
! The current opacities and emissivities are stored for the variation of the
! radiation field at the next frequency.
!
	DO I=1,ND
	  CHI_PREV(I)=CHI_CONT(I)
	  ETA_PREV(I)=ETA_CONT(I)
	END DO
!
! Store opacities and emissivities for use in observer's frame
! calculation.
!
	IF(ES_COUNTER .EQ. 1)THEN
	  IF(INCL_RAY_SCAT)RAY_CMF_ST(1:ND,ML)=CHI_RAY(1:ND)
	  IF(INCL_DUST)DUST_CMF_ST(1:ND,ML)=CHI_SCAT_DUST(1:ND)
	  ETA_CMF_ST(1:ND,ML)=ETA(1:ND)
	  CHI_CMF_ST(1:ND,ML)=CHI(1:ND)
	END IF
	RJ_CMF_ST(1:ND,ML)=RJ(1:ND)
!
10000	CONTINUE
	T1=1.0D0; WRITE(LU_EDD,REC=FINISH_REC)T1
	CALL TUNE(ITWO,'MLCF')
!
! NB: We use K here, rather than ACCESS_F, so that we don't corrupt EDDFACTOR if
!      evaluate EW is set to TRUE.
!
	IF(WRITE_ETA_AND_CHI .AND. (ES_COUNTER .EQ. NUM_ES_ITERATIONS .OR. .NOT. COMPUTE_J))THEN
	  CALL GET_LU(J,'ETA_DATA write')
	  CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1              REXT,VEXT,LANG_COORDEXT,ND,
	1              K,L_TRUE,L_TRUE,L_FALSE,'ETA_DATA',J)
	  K=INITIAL_ACCESS_REC
	  WRITE(J,REC=EDD_CONT_REC)K,NCF,ND
	  DO ML=1,NCF
	    WRITE(J,REC=K-1+ML)(ETA_CMF_ST(I,ML),I=1,ND),NU(ML)
	  END DO
	  T1=1.0D0; WRITE(J,REC=FINISH_REC)T1
	  CLOSE(UNIT=J)
!
	  CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1              REXT,VEXT,LANG_COORDEXT,ND,
	1              K,L_TRUE,L_TRUE,L_FALSE,'CHI_DATA',J)
	  K=INITIAL_ACCESS_REC
	  WRITE(J,REC=EDD_CONT_REC)K,NCF,ND
	  DO ML=1,NCF
	    WRITE(J,REC=K-1+ML)(CHI_CMF_ST(I,ML),I=1,ND),NU(ML)
	  END DO
	  T1=1.0D0; WRITE(J,REC=FINISH_REC)T1
	  CLOSE(UNIT=J)
!
	  IF(INCL_RAY_SCAT)THEN
	    CALL OPEN_RW_EDDFACTOR(R,V,LANG_COORD,ND,
	1              REXT,VEXT,LANG_COORDEXT,ND,
	1              K,L_TRUE,L_TRUE,L_FALSE,'RAY_DATA',J)
	    K=INITIAL_ACCESS_REC
	    WRITE(J,REC=EDD_CONT_REC)K,NCF,ND
	    DO ML=1,NCF
	      WRITE(J,REC=K-1+ML)(RAY_CMF_ST(I,ML),I=1,ND),NU(ML)
	    END DO
	    T1=1.0D0; WRITE(J,REC=FINISH_REC)T1
	    CLOSE(UNIT=J)
	  END IF
!
	END IF
!
! We output the PLANCKMEAN -- this was need for testing.
!
	IF(.NOT. COMPUTE_J)THEN
	  I=3; CALL TUNE(I,' ')
	  T1=7.218771E+11_LDP          ! T1=4 * [STEFAN BOLTZMAN CONS] * 1.0D+16 / pi
	  DO I=1,ND
	    PLANCKMEAN(I)=4.0_LDP*PLANCKMEAN(I)/T1/( T(I)**4 )
	  END DO
	  OPEN(UNIT=11,STATUS='UNKNOWN',FILE='PLANCK_KAPPA_MEAN')
	    WRITE(11,'(I5,2ES18.8,ES14.4)')(I,R(I),V(I),1.0D-10*PLANCKMEAN(I)/DENSITY(I),I=1,ND)
	  CLOSE(UNIT=11)
	  STOP
	END IF
!
	COMPUTE_EDDFAC=.FALSE.
!
! Compute the convolution of J with the electron redistribution function.
! If USE_FIXED_J is true, the frequencies in the EDDFACTOR file will not
! match those of this model. In addition, since we are not iterating, there
! is no need to create an ES_J_CONV file.
!
	COHERENT_ES=RD_COHERENT_ES
	IF(.NOT. COHERENT_ES .AND. .NOT. USE_FIXED_J)THEN
	  I=ND*ND
	  CALL COMP_J_CONV_V2(WM,I,NU,TEXT,NDEXT,NCF,LU_EDD,'EDDFACTOR',
	1          EDD_CONT_REC,L_FALSE,L_FALSE,LU_ES,'ES_J_CONV')
          CALL OUT_RV_TO_EDDFACTOR(R,V,LANG_COORD,ND,
	1        REXT,VEXT,LANG_COORDEXT,NDEXT,
	1        ACCESS_F,'ES_J_CONV',LU_ES)
	  CLOSE(LU_ES)	
!
	END IF
!
	END DO		!ES_COUNTER
!
!
!
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(STRING,'(I10)')N_OBS
	  DO WHILE(STRING(1:1) .EQ. ' ')
	       STRING(1:)=STRING(2:)
	  END DO
	  STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	  CALL WRITV_V2(OBS_FREQ,N_OBS,ISEV,TRIM(STRING),LU_FLUX)
	  CALL WRITV_V2(OBS_FLUX,N_OBS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	  CALL WRITV(RLUMST,ND,'Luminosity',LU_FLUX)
	CLOSE(UNIT=LU_FLUX)
!
! Make sure ALL frequencies have been output to CMF_FORCE_MULT
!
	IF(WRITE_CMF_FORCE)THEN
	  K=5   					!ACCESS_F
	  DO ML=ML_FORCE,N_FORCE
	    WRITE(82,REC=K-1+ML)( LINE_FLUXMEAN(I)/STARS_LUM/ESEC(I) ,I=1,ND),NU_FORCE
	    NU_FORCE=NU_FORCE*NU_FORCE_FAC
	  END DO
	  CLOSE(UNIT=82)
	END IF
!
! Output errors that have occurred in MOM_J_CMF
!
        I=1; CALL WRITE_J_CMF_ERR(I)
!
! Compute ROSSELAND and FLUX mean opacities. Compute the respective
! optical depth scales; TA for the FLUX mean optical depth scale,
! and TB for the ROSSELAND mean optical depth scale.
!
! T1=4 * [STEFAN BOLTZMAN CONS] * 1.0D+16 / pi
!
	T1=7.218771E+11_LDP
	DO I=1,ND
	  IF(ABS(RLUMST(I)) .LE. 1.0E-20_LDP)RLUMST(I)=1.0E-20_LDP
	  FLUXMEAN(I)=FLUXMEAN(I)/RLUMST(I)
!
! Bug fix, not flag-gated: ROSSMEAN is an accumulated 1/CHI-weighted
! integral that can vanish when every frequency's Planck weight
! underflows.  A zero denominator made both outputs Inf/NaN, and
! ROSSMEAN feeds TCHI and the derived tau scales below.  A dead cell
! now gets ROSSMEAN=0, which the existing MINVAL(TCHI)>0 test already
! traps as a diagnostic instead of propagating NaN.
!
	  IF(ROSSMEAN(I) .GT. 1.0E-30_LDP)THEN
	    INT_dBdT(I)=INT_dBdT(I)/ROSSMEAN(I)		!Program rosseland opac.
	    ROSSMEAN(I)=T1*( T(I)**3 )/ROSSMEAN(I)
	  ELSE
	    INT_dBdT(I)=0.0_LDP
	    ROSSMEAN(I)=0.0_LDP
	  END IF
	  PLANCKMEAN(I)=4.0_LDP*PLANCKMEAN(I)/T1/( T(I)**4 )
	END DO
!
	IF(WR_ION_LINE_FORCE)THEN
	    CALL OUT_LINE_FORCE(ION_LINE_FORCE,FLUXMEAN,ROSSMEAN,RLUMST,ESEC,
	1          R,V,DENSITY,ION_ID,ND,NUM_IONS)
	END IF
!
! Due to instabilities, the FLUX mean opacity can sometimes be
! negative. If so we can continue, but we note that the
! results may be in error, and need to be checked.
!
	TCHI(1:ND)=ROSSMEAN(1:ND)*CLUMP_FAC(1:ND)
	IF(MINVAL(TCHI) .GT. 0.0_LDP)THEN
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
        ELSE
	  dCHIDR(1:ND)=0.0_LDP
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB: Check MEANOPAC'
	  WRITE(LUER,*)'Rosseland mean opacity has negative values'
	END IF
	CALL NORDTAU(TA,TCHI,R,R,dCHIdR,ND)
!
	TCHI(1:ND)=FLUXMEAN(1:ND)*CLUMP_FAC(1:ND)
	IF(MINVAL(TCHI) .GT. 0.0_LDP)THEN
	  CALL DERIVCHI(dCHIdR,FLUXMEAN,R,ND,METHOD)
        ELSE
	  dCHIDR(1:ND)=0.0_LDP
	  WRITE(LUER,*)'Warning from CMF_FLUX_SUB: Check MEANOPAC'
	  WRITE(LUER,*)'Flux mean opacity has zero or negative values'
	END IF
        CALL NORDTAU(TB,TCHI,R,R,dCHIdR,ND)
!
	TCHI(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
        CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
!
	TA(ND)=0.0_LDP; TB(ND)=0.0_LDP; TC(ND)=0.0_LDP; DTAU(ND)=0.0_LDP
	CALL GEN_ASCI_OPEN(LU_OPAC,'MEANOPAC','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(LU_OPAC,
	1  '( ''         R          I   Tau(Ross)   /\Tau  Rat(Ross)'//
	1  ' Chi(Ross)  Chi(ross)  Chi(Flux)   Chi(es) '//
	1  '  Tau(Flux)  Tau(es)  Rat(Flux)  Rat(es)   Kappa(R) V(km/s)'' )' )
	IF(R(1) .GE. 1.0E+05_LDP)THEN
	  FMT='(ES17.10,I4,2ES10.3,ES10.2,4ES11.3,4ES10.2,2ES11.3)'
	ELSE
	  FMT='( F17.10,I4,2ES10.3,ES10.2,4ES11.3,4ES10.2,2ES11.3)'
	END IF
	  DO I=1,ND
	    IF(I .EQ. 1)THEN
	      T1=0.0_LDP		!Rosseland optical depth scale
	      T2=0.0_LDP		!Flux optical depth scale
	      T3=0.0_LDP		!Electron scattering optical depth scale.
	      TC(1:3)=0.0_LDP
	    ELSE
	      T1=T1+TA(I-1)
	      T2=T2+TB(I-1)
	      T3=T3+DTAU(I-1)
	      TC(1)=TA(I)/TA(I-1)
	      TC(2)=TB(I)/TB(I-1)
	      TC(3)=DTAU(I)/DTAU(I-1)
	    END IF
	    WRITE(LU_OPAC,FMT)R(I),I,T1,TA(I),TC(1),
	1      ROSSMEAN(I),INT_dBdT(I),FLUXMEAN(I),ESEC(I),
	1      T2,T3,TC(2),TC(3),1.0D-10*ROSSMEAN(I)/DENSITY(I),V(I)
	  END DO
	CLOSE(UNIT=LU_OPAC)
	WRITE(LU_OPAC,'(//,A,A)')
	1     'NB: Mean opacities do not include effect of clumping',
	1     'NB: Optical depth scale includes effect of clumping'
!
! Output hydrodynamical terms to allow check on radiation driving of the wind.
!
	I=18
	MEAN_ATOMIC_WEIGHT=DENSITY(ND)/POP_ATOM(ND)/ATOMIC_MASS_UNIT()
	CALL HYDRO_TERMS(POP_ATOM,R,V,T,SIGMA,ED,RLUMST,
	1                 STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1		  FLUXMEAN,ESEC,I,ND)
!
! 
!
! Update the emissivity for Rayleigh scattering which is is assumed to be coherent
! in the comoving frame.
!
	IF(INCL_RAY_SCAT)THEN
	  DO ML=1,NCF
	    ETA_CMF_ST(1:ND,ML)=ETA_CMF_ST(1:ND,ML)+RAY_CMF_ST(1:ND,ML)*RJ_CMF_ST(1:ND,ML)
	  END DO
	END IF
!
! At present we assume isotropic dust scattering in OBSFRAME.
!
	IF(INCL_DUST)THEN
	  DO ML=1,NCF
	    ETA_CMF_ST(1:ND,ML)=ETA_CMF_ST(1:ND,ML)+DUST_CMF_ST(1:ND,ML)*RJ_CMF_ST(1:ND,ML)
	  END DO
	END IF
!
! If requested, convolve J with the electron-scattering redistribution
! funtion. K is used for LUIN, and LUOUIT but is not accessed.
!
	IF(.NOT. COHERENT_ES)THEN
	  I=ND*NCF
	  CALL COMP_J_CONV_V2(RJ_CMF_ST,I,NU,T,ND,NCF,
	1         K,'J PASSED VIA CALL',EDD_CONT_REC,L_FALSE,L_FALSE,
	1         K,'RETURN J VIA CALL')
	END IF
!
	ESEC(1:ND)=6.65E-15_LDP*ED(1:ND)
	DO ML=1,NCF
	  ETA_CMF_ST(1:ND,ML)=ETA_CMF_ST(1:ND,ML) +
	1                        RJ_CMF_ST(1:ND,ML)*ESEC(1:ND)
	END DO
!
	DEALLOCATE (RJ_CMF_ST)
!
! 
! ***************************************************************************
! ***************************************************************************
!
! Determine the number of points inserted along each ray via
! MAX_DEL_V_RES_ZONE.
!
	DO I=1,ND
	  MAX_DEL_V_RES_ZONE(I)=VTURB_VEC(I)*FRAC_DOP*0.5_LDP
	END DO
!
	IF(DO_CLUMP_MODEL)THEN
	  DO ML=1,NCF
	    ETA_CMF_ST(:,ML)=ETA_CMF_ST(:,ML)*CLUMP_FAC(:)
	    CHI_CMF_ST(:,ML)=CHI_CMF_ST(:,ML)*CLUMP_FAC(:)
	  END DO
	END IF
!
!
! Ensure that MU and the quadrature weights are correctly defined.
! By definition, p * dp equals R**2 * mu * dmu. Integration over mu is
! more stable, and is to be preferred.
!
! For a plane-parallel model, we only have NC integration angles.
!
	IF(PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL)THEN
	  HQW_AT_RMAX(:)=0.0_LDP
	  CALL GAULEG(RZERO,RONE,MU_AT_RMAX,HQW_AT_RMAX,NC)
	  HQW_AT_RMAX(1:NC)=HQW_AT_RMAX(1:NC)*MU_AT_RMAX(1:NC)
	  HQW_AT_RMAX(1:NC)=HQW_AT_RMAX(1:NC)*R(ND)*R(ND)/R(1)/R(1)
	  P(1:NC)=R(ND)*SQRT( (1.0_LDP-MU_AT_RMAX(1:NC))*(1.0_LDP+MU_AT_RMAX(1:NC)) )
	  NC_OBS=NC; NP_OBS=NP
	  P_OBS(1:NP_OBS)=P(1:NP_OBS)
	ELSE IF(REVISE_P_GRID)THEN
	  DEALLOCATE (P_OBS)
	  I=NP*10; ALLOCATE (P_OBS(I))
	  CALL REVISE_OBS_P(P_OBS,NP_OBS,I,NC_OBS,NC,R,ND,LUIN,LUMOD)
!
	  DEALLOCATE (MU_AT_RMAX,dMU_AT_RMAX,HQW_AT_RMAX)
	  ALLOCATE (MU_AT_RMAX(NP_OBS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dMU_AT_RMAX(NP_OBS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HQW_AT_RMAX(NP_OBS),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocate MU_AT_RMAX in refine p grid section'
	    STOP
	  END IF
	  DO LS=1,NP_OBS
	    MU_AT_RMAX(LS)=SQRT( MAX(0.0_LDP,
	1        (R(1)-P_OBS(LS))*(R(1)+P_OBS(LS))) )/R(1)
	  END DO
	  CALL SET_ACC_dMU(MU_AT_RMAX,dMU_AT_RMAX,P_OBS,R(1),NP_OBS)
	  CALL HWEIGHT_V2(MU_AT_RMAX,dMU_AT_RMAX,HQW_AT_RMAX,NP_OBS)
!
	ELSE
	  DO LS=1,NP
	    MU_AT_RMAX(LS)=SQRT( MAX(0.0_LDP,
	1        (R(1)-P(LS))*(R(1)+P(LS))) )/R(1)
	  END DO
	  CALL SET_ACC_dMU(MU_AT_RMAX,dMU_AT_RMAX,P,R(1),NP)
	  CALL HWEIGHT_V2(MU_AT_RMAX,dMU_AT_RMAX,HQW_AT_RMAX,NP)
	  NC_OBS=NC; NP_OBS=NP
	  P_OBS(1:NP_OBS)=P(1:NP_OBS)
	END IF
!
! Output R, V and the Langrangian coordiante to the EDDACTOR file.
! If we adjust R, these values are consistent with the old grid --
! not the current grid. These only need to be output heew we are
! using an old EDDFACTOR file.
!
          CALL OUT_RV_TO_EDDFACTOR(R,V,LANG_COORD,ND,
	1        REXT,VEXT,LANG_COORDEXT,NDEXT,
	1        ACCESS_F,'EDDFACTOR',LU_EDD)
	  CLOSE(LU_EDD)	
!
! ***************************************************************************
! ***************************************************************************
!
! Compute the observer's frame fluxes. The fluxes are returned in Janskies.
! V6 of the observer's frame code can now handle a plane-parallel model atmosphere,
! with, or without, a velocity field.
!
! We use TA for V in the call to OBS_FRAME_SUB_V6.
!
	TMP_LOG=.FALSE.
	IF(PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL)TMP_LOG=.TRUE.
	IF(PLANE_PARALLEL_NO_V)THEN
	   TA(1:ND)=0.0_LDP
	ELSE
	   TA(1:ND)=V(1:ND)
	END IF
	CALL OBS_FRAME_SUB_V9(ETA_CMF_ST,CHI_CMF_ST,NU,
	1            R,TA,T,ED,ND,NCF,
	1            P_OBS,MU_AT_RMAX,HQW_AT_RMAX,NC_OBS,NP_OBS,
	1            OBS_FREQ,OBS_FLUX,N_OBS,
	1            MAX_DEL_V_RES_ZONE,OBS_TAU_MAX,OBS_ES_DTAU,
	1            N_INS_OBS,OBS_INT_METHOD,
	1            TAU_REF,WRITE_RTAU,WRITE_IP,WRITE_dFR,
	1            DO_REL_IN_OBSFRAME,TMP_LOG)
!
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFRAME','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(STRING,'(I10)')N_OBS
	  DO WHILE(STRING(1:1) .EQ. ' ')
	       STRING(1:)=STRING(2:)
	  END DO
	  STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	  CALL WRITV_V2(OBS_FREQ,N_OBS,ISEV,TRIM(STRING),LU_FLUX)
	  CALL WRITV_V2(OBS_FLUX,N_OBS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	CLOSE(UNIT=LU_FLUX)
!
! At present the line EW's are automatically computed when the
! CONTINUM Flux is evaluated.
!
	IF(.NOT. DO_SOBOLEV_LINES)RETURN
!
! This section of the program compute the SOBOLEV line EW's. These
! are only approximate, and can be used to determine which lines are
! imortant contributors to INDIVIDUAL emission features. NB: If a line
! has a P Cygni profile it may be an important contributor but give an
! EW of zero.
!
! Set access counter for Continuum Eddington file, ensuring that file
! is OPENED. FIRST is used as a temporary LOGICAL variable.
!
	COMPUTE_EDDFAC=.TRUE.
	EDDINGTON=EDD_LINECONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  INQUIRE(UNIT=LU_EDD,OPENED=FIRST)
	  IF(.NOT. FIRST)THEN
	    I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
	    OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I)
	  END IF
	  WRITE(LU_EDD,REC=4)ACCESS_F,N_LINE_FREQ,NDEXT
	END IF
!
! Open file to store EW data
!
	CALL GEN_ASCI_OPEN(LU_EW,'EWDATA','UNKNOWN',' ',' ',IZERO,IOS)
!
! SOB_DEP collects each line's deposited EW (in Jy per obs frequency)
! across the whole line loop. It is applied to OBS_FLUX after the loop,
! so the deposited spectrum is written with the same OBS_FREQ grid.
!
	IF(DO_SOB_FLUX_OBS)THEN
	  ALLOCATE(SOB_DEP(N_OBS))
	  SOB_DEP(1:N_OBS)=0.0_LDP
	END IF
!
! Zero the vector that will be used to store the force-multiplier computed
! for all lines using the SOBOLEV approximation.
!
	FORCE_MULT(1:ND)=0.0_LDP
	N_FORCE=LOG(NU(NCF)/NU(1))/LOG(1.0_LDP- 500.0_LDP/C_KMS)
	NU_FORCE=NU(1)
	NU_FORCE_FAC=(1.0_LDP-500.0_LDP/C_KMS)
	ML_FORCE=1
!
! Enter line loop.
!
	LINE_INDX=1
	NUM_SOB_LINES=N_LINE_FREQ
	IF(.NOT. DO_ALL_SOB_LINES)THEN
	  T1=SPEED_OF_LIGHT()*1.0E-07_LDP
	  DO I=1,N_LINE_FREQ
	    IF(T1/VEC_FREQ(I) .GE. SOB_EW_LAM_BEG)EXIT
	    LINE_INDX=I
	  END DO
	  DO I=N_LINE_FREQ,LINE_INDX,-1
	    IF(T1/VEC_FREQ(I) .LE. SOB_EW_LAM_END)EXIT
	    NUM_SOB_LINES=I
	  END DO
	END IF
!
	OVERLAP=.FALSE.
	CALL TUNE(1,'SOB_EW')
	DO WHILE (LINE_INDX .LE. NUM_SOB_LINES)      !N_LINE_FREQ)
!
	IF(OVERLAP)THEN
	  TMP_MAX_SIM=MAX_SIM
	ELSE
	  TMP_MAX_SIM=1
	  OVER_FREQ_DIF=0.0_LDP
	END IF
!
! Determine next line (or lines), and store line parameters (eg frequency,
! Einstein A coefficient etc) in appropriate locations for the subsequent
! treatment of line.
!
! The temporary variable J [=MAX(N_LINE_FREQ,LINE_INDX+SIM_INDX) ] is used so that
! VEC_FREQ(N_LINE_FREQ+1) is not accessed (remember that the arguments of
! a "DO WHILE" can be done in any order.
!
	SIM_INDX=0
	J=LINE_INDX
	SOBOLEV=.TRUE.
	DO WHILE( LINE_INDX+SIM_INDX .LE. NUM_SOB_LINES .AND.
	1      (VEC_FREQ(LINE_INDX)-VEC_FREQ(J))/VEC_FREQ(LINE_INDX)
	1          .LE. OVER_FREQ_DIF
	1          .AND. SIM_INDX .LT. TMP_MAX_SIM )
	  SIM_INDX=SIM_INDX+1
	  K=LINE_INDX+SIM_INDX-1
	  SIM_NL(SIM_INDX)=VEC_NL(K)
	  SIM_NUP(SIM_INDX)=VEC_NUP(K)
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  EINA(SIM_INDX)=VEC_EINA(K)
	  OSCIL(SIM_INDX)=VEC_OSCIL(K)
	  FL_SIM(SIM_INDX)=VEC_FREQ(K)
	  SIM_LINE_POINTER(SIM_INDX)=K
	  TRANS_NAME_SIM(SIM_INDX)=TRIM(VEC_SPEC(K))
	  J=MIN(LINE_INDX+SIM_INDX,N_LINE_FREQ)
	END DO
	NUM_SIM_LINES=SIM_INDX
!
! Determine center of blend. Used as continuum wavelength.
!
	FL=0.0_LDP
	DO SIM_INDX=1,NUM_SIM_LINES
	  FL=FL+FL_SIM(SIM_INDX)
	END DO
	FL=FL/NUM_SIM_LINES
!
	CONT_FREQ=FL
	COMPUTE_NEW_CROSS=.TRUE.
!
! Determine method to compute continuum intensity.
!
	IF(ACCURATE .AND. ALL_FREQ)THEN
	  THIS_FREQ_EXT=.TRUE.
	ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	  THIS_FREQ_EXT=.TRUE.
	ELSE
	  THIS_FREQ_EXT=.FALSE.
	END IF
! 
! Compute U_STAR_RATIO and L_STAR_RATIO which are used to switch from
! the opacity/emissivity computed with a FULL_ATOM to an equivalent form
! but written in terms of the SUPER-LEVELS.
!
! L refers to the lower level of the transition.
! U refers to the upper level of the transition.
!
! At present we must treat each species separately (for those with both FULL
! and SUPER_LEVEL model atoms).
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  I=SIM_LINE_POINTER(SIM_INDX)
	  MNL_F=VEC_MNL_F(I)
	  MNUP_F=VEC_MNUP_F(I)
	  L_STAR_RATIO(1:ND,SIM_INDX)=0.0_LDP
	  U_STAR_RATIO(1:ND,SIM_INDX)=0.0_LDP
!
	  DO ID=1,NUM_IONS
	    IF(VEC_SPEC(I) .EQ. ION_ID(ID))THEN
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
! T1 is used to represent b(level)/b(super level). If no interpolation of
! the b values in a super level has been performed, this ratio will be unity .
! This ratio is NOT treated in the linearization.
!
	      DO K=1,ND
	        T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzV(MNL,K))/ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
	        L_STAR_RATIO(K,SIM_INDX)=T1*(ATM(ID)%W_XzV_F(MNUP_F,K)/ATM(ID)%W_XzV_F(MNL_F,K))*ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
	        T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzV(MNUP,K))/ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	        U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
    	      END DO
	      GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	      TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1             '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1                  TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	      EXIT
	    END IF
	  END DO
	END DO
!
! Compute line opacity and emissivity.
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  T1=OSCIL(SIM_INDX)*OPLIN
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=1,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    NEG_OPACITY(I)=.FALSE.
	  END DO
	END DO
! 
!
	SECTION='LINE'
!
! Compute continuum opacity and emissivity at the line frequency.
!
	INCLUDE 'OPACITIES_V5.INC'
!
! Compute continuum intensity.
!
	INCLUDE 'COMP_JCONT_V4.INC'	
!
! SOURCE is used by SOBEW
!
	DO I=1,ND
	  SOURCE(I)=ZETA(I)+RJ(I)*THETA(I)
	END DO
! 
!
! Scale emissivity because of different frequencies.
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  DO I=1,ND
	    T1=FL**3/(EXP(HDKT*FL/T(I))-1.0_LDP)
	    T2=FL_SIM(SIM_INDX)**3/(EXP(HDKT*FL_SIM(SIM_INDX)/T(I))-1.0_LDP)
	    BB_COR(I,SIM_INDX)=T1/T2
	  END DO
	END DO
!
! AT present overlapping lines in the CMF frame is not installed. We still
! use CHIL and ETAL --- not CHIL_MAT etc.
!
	IF(SOBOLEV .OR. NUM_SIM_LINES .EQ. 1)THEN
	  DO I=1,ND
	    CHIL(I)=0.0_LDP
	    ETAL(I)=0.0_LDP
	    DO SIM_INDX=1,NUM_SIM_LINES
	      CHIL(I)=CHIL(I)+CHIL_MAT(I,SIM_INDX)
	      ETAL(I)=ETAL(I)+ETAL_MAT(I,SIM_INDX)*BB_COR(I,SIM_INDX)
	    END DO
	  END DO
	END IF
!
        IF(CHECK_LINE_OPAC .AND. SOBOLEV)THEN
	  T1=3.0E-10_LDP/FL
	  FIRST_NEG=.TRUE.
	  DO I=1,ND
	    T2=CHIL(I)*T1*R(I)/V(I)		!Tau(Sobolev) for dlnV/dlnr=1
	    IF(T2 .LT. -0.5_LDP)THEN
	      NEG_OPACITY(I)=.TRUE.
	      IF(LST_ITERATION)THEN
	        IF(FIRST_NEG)THEN
	          J=ICHRLEN(TRANS_NAME_SIM(1))
	          WRITE(LU_NEG,'(1X,A)')TRANS_NAME_SIM(1)(1:J)
	          FIRST_NEG=.FALSE.
	        END IF
	        WRITE(LU_NEG,2290)I,T2,ED(I)
2290	        FORMAT(1X,'I= ',I3,'  : TAU(Sob)= ',1P,E9.2,'  : Ne=',E8.2)
	        WRITE(LU_NEG,*)CHIL(I),POPS(NL,I),R(I),V(I),T1
	      END IF
	      CHIL(I)=1.0_LDP	!Reset after we output its value.
	    END IF
	  END DO
!
! Note that both CHIL and CHIL_MAT are used by SOBJBAR_SIM to compute ZNET.
! The following ensures consistency. Alternative is to make positive
! only those opacities that are negative --- then need  NEG_OPACITY indicator
! for each line?
!
	  DO SIM_INDX=1,NUM_SIM_LINES
	    DO I=1,ND
              IF(NEG_OPACITY(I))CHIL_MAT(I,SIM_INDX)=1.0_LDP/NUM_SIM_LINES
	    END DO
	  END DO
!
	END IF
!
! Estimate the line EW using a Modified Sobolev approximation.
!
	CHI_CLUMP(1:ND)=CHI(1:ND)*CLUMP_FAC(1:ND)
	ETA_CLUMP(1:ND)=ETA(1:ND)*CLUMP_FAC(1:ND)
	ESEC_CLUMP(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	CHI_SCAT_CLUMP(1:ND)=CHI_SCAT(1:ND)*CLUMP_FAC(1:ND)
!
	CHIL(1:ND)=CHIL(1:ND)*CLUMP_FAC(1:ND)
	ETAL(1:ND)=ETAL(1:ND)*CLUMP_FAC(1:ND)
!
	IF(FL .LT. NU_FORCE .AND. WRITE_SOB_FORCE)THEN
	  INQUIRE(UNIT=82,OPENED=TMP_LOG)
	  K=5   					!ACCESS_F
	  IF(.NOT. TMP_LOG)THEN
	    I=WORD_SIZE*(ND+1)/UNIT_SIZE; J=82
	    CALL OPEN_DIR_ACC_V1(ND,I,DA_FILE_DATE,'SOB_FORCE_DATA',J)
	    WRITE(82,REC=EDD_CONT_REC)K,N_FORCE,ND
	  END IF
	  WRITE(82,REC=K-1+ML_FORCE)(FORCE_MULT(I),I=1,ND),NU_FORCE
	  ML_FORCE=ML_FORCE+1
	  NU_FORCE=NU_FORCE*NU_FORCE_FAC
	END IF
!
! We use TA as a temporary vector which indicates the origin
! of the line emission. Not required in this code as used only for display purposes.
!
! Note: Using CMF_ABS_EW takes much longer and was added for determining strongest lines
!          to include in polarizatin calculations.
!
! With DO_NONMONO_SOB the Sobolev EW estimators use the unclamped
! SIGMA_TRUE (with the same max(|1+sig*mu^2|,-0.2(1+sig)) tau floor as
! the populations), matching the gradient the line rates used.  The
! default path passes the clamped SIGMA, identical to before.
!
	IF(DO_CMF_EW)THEN
	  INCLUDE 'CMF_ABS_EW.INC'
	ELSE IF(DO_NONMONO_SOB)THEN
	  CALL SOBEW_GRAD_V2(SOURCE,CHI_CLUMP,CHI_SCAT_CLUMP,CHIL,ETAL,
	1            V,SIGMA_TRUE,R,P,FORCE_MULT,STARS_LUM,AQW,HQW,TA,EW,
	1            CONT_INT,
	1            FL,INNER_BND_METH,DBB,IC,THK_CONT,L_FALSE,NC,NP,ND,METHOD)
	  ABS_EW=EW
	ELSE
	  CALL SOBEW_GRAD_V2(SOURCE,CHI_CLUMP,CHI_SCAT_CLUMP,CHIL,ETAL,
	1            V,SIGMA,R,P,FORCE_MULT,STARS_LUM,AQW,HQW,TA,EW,CONT_INT,
	1            FL,INNER_BND_METH,DBB,IC,THK_CONT,L_FALSE,NC,NP,ND,METHOD)
	  ABS_EW=EW
	END IF
!
! Deposit the Sobolev EW into OBS_FLUX. TRANS_TYPE='SOB' lines are
! excluded from the formal transfer above, so without this step the
! observed spectrum contains no line flux from SOB lines at all. The
! EW the populations converged under is spread over the observed-frame
! resonance span defined by the local velocity field, which keeps the
! render consistent with the population solve on a non-monotonic
! velocity field (where CMF/BLANK line transfer over-emits).
!
	IF(DO_SOB_FLUX_OBS .AND. ABS(ABS_EW) .GE. EW_CUT_OFF)THEN
	  CALL DEP_SOB_OBS_FLUX(SOB_DEP,OBS_FREQ,N_OBS,
	1        EW,CONT_INT,FL,ETAL,CHIL,V,R,ND,C_KMS)
	END IF
!
	IF(ABS(ABS_EW) .GE. EW_CUT_OFF)THEN
	  T1=0.01_LDP*C_KMS/FL_SIM(1) 		!Wavelength(Angstroms)
	  DO SIM_INDX=1,NUM_SIM_LINES
	    J=SIM_LINE_POINTER(SIM_INDX)
	    MNL=VEC_MNL_F(J)
	    MNUP=VEC_MNUP_F(J)
	    IF(SIM_INDX .EQ. 1)THEN
	      EW_STRING=TRANS_NAME_SIM(SIM_INDX)	!Used for transition name
	    ELSE
	      EW_STRING='## '//TRIM(TRANS_NAME_SIM(SIM_INDX))
	    END IF
	    CALL WRITE_OUT_EW(EW_STRING,T1,CONT_INT,EW,ABS_EW,SOBOLEV,MNL,MNUP,LU_EW)
	  END DO
	END IF
!
! Update line counter so that we move onto next line
!
	LINE_INDX=LINE_INDX+NUM_SIM_LINES
!
	END DO		!Line loop
!
! Apply the deposited Sobolev line flux and rewrite OBSFRAME/OBSFLUX.
! Both files were written before the EW loop (OBS_FLUX then carried
! continuum and BLANK lines only), so the deposit must be written out
! again here. The floor prevents stacked absorption lobes from driving
! a bin below zero.
!
	IF(DO_SOB_FLUX_OBS)THEN
	  DO J=1,N_OBS
	    OBS_FLUX(J)=MAX(OBS_FLUX(J)+SOB_DEP(J),
	1                  1.0D-6*OBS_FLUX(J))
	  END DO
	  DEALLOCATE(SOB_DEP)
!
	  CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(STRING,'(I10)')N_OBS
	    DO WHILE(STRING(1:1) .EQ. ' ')
	      STRING(1:)=STRING(2:)
	    END DO
	    STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	    CALL WRITV_V2(OBS_FREQ,N_OBS,ISEV,TRIM(STRING),LU_FLUX)
	    CALL WRITV_V2(OBS_FLUX,N_OBS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	  CLOSE(UNIT=LU_FLUX)
!
	  CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFRAME','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(STRING,'(I10)')N_OBS
	    DO WHILE(STRING(1:1) .EQ. ' ')
	      STRING(1:)=STRING(2:)
	    END DO
	    STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	    CALL WRITV_V2(OBS_FREQ,N_OBS,ISEV,TRIM(STRING),LU_FLUX)
	    CALL WRITV_V2(OBS_FLUX,N_OBS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	  CLOSE(UNIT=LU_FLUX)
	END IF
!
! Outout the Sobolev Line-Force multiplier
!
	IF(WRITE_SOB_FORCE)THEN
	  K=5   					!ACCESS_F
	  DO ML=ML_FORCE,N_FORCE
	    WRITE(82,REC=K-1+ML)(FORCE_MULT(I),I=1,ND),NU_FORCE
	    NU_FORCE=NU_FORCE*NU_FORCE_FAC
	  END DO
	  CLOSE(UNIT=82)
	  OPEN(UNIT=82,FILE='SOB_FORCE_MULT',STATUS='UNKNOWN')
	   DO I=1,ND
	      WRITE(82,'(1X,I3,2X,3ES14.6)')I,R(I),V(I),FORCE_MULT(I)
	    END DO
	  CLOSE(UNIT=82)
	END IF
!	
	CALL TUNE(2,'SOB_EW')
!
	RETURN
!
	END
!
!
!
! Deposit the Sobolev EW of one line (or overlapping-line group) into
! the observed-frame flux grid. TRANS_TYPE='SOB' lines are excluded
! from the formal transfer in CMF_FLUX_SUB_V5, so this routine is the
! only path for their flux into the observed spectrum.
!
! Each depth contributes its share of the line's net EW over the
! observed-frame span its projected velocities can reach:
!   - net emission (EW > 0): a shell at speed |V(I)| presents every
!     line-of-sight velocity in [-|V|, +|V|], so its share spreads
!     uniformly (flat box) over LAM0*(1-|V|/c) to LAM0*(1+|V|/c).
!     Depths are weighted by the line emissivity ETAL*R^2 so the
!     profile peaks where the line actually forms.
!   - net absorption (EW < 0): only gas in front of the continuum
!     photosphere absorbs, so the box is one-sided: [LAM0*(1-V/c),
!     LAM0] for outflowing gas (V > 0, blueshifted trough) and
!     [LAM0, LAM0*(1+|V|/c)] for infalling gas (V < 0), so velocity
!     reversals produce redshifted absorption. Weighted by the line
!     opacity CHIL*R^2.
!
! The deposited quantity is a flux density per Angstrom; its integral
! over the full span is EW*CONT_INT, the same line strength the
! population solve used (SIGMA_TRUE under DO_NONMONO_SOB). A minimum
! half-width of 100 km/s keeps depths near a velocity-reversal node
! depositing at line center instead of vanishing.
!
	SUBROUTINE DEP_SOB_OBS_FLUX(SOB_DEP,OBS_FREQ,N_OBS,
	1         EW,CONT_INT,FL,ETAL,CHIL,V,R,ND,C_KMS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: N_OBS,ND
	REAL(KIND=LDP), INTENT(INOUT) :: SOB_DEP(N_OBS)
	REAL(KIND=LDP), INTENT(IN) :: OBS_FREQ(N_OBS)
	REAL(KIND=LDP), INTENT(IN) :: EW,CONT_INT,FL,C_KMS
	REAL(KIND=LDP), INTENT(IN) :: ETAL(ND),CHIL(ND),V(ND),R(ND)
!
	INTEGER I,M,M1,M2,LO,HI,MID
	REAL(KIND=LDP) LAM0,LAM_LO,LAM_HI,NU_LO,NU_HI
	REAL(KIND=LDP) WT(ND),WT_SUM,BETA,DEP_DENSITY
	REAL(KIND=LDP), PARAMETER :: MIN_BETA=100.0_LDP/2.99792458D5
!
! Rest wavelength in Angstroms, matching the EWDATA convention.
!
	LAM0=0.01_LDP*C_KMS/FL
!
	IF(EW .GT. 0.0_LDP)THEN
	  WT(1:ND)=MAX(ETAL(1:ND),0.0_LDP)*R(1:ND)*R(1:ND)
	ELSE
	  WT(1:ND)=MAX(CHIL(1:ND),0.0_LDP)*R(1:ND)*R(1:ND)
	END IF
	WT_SUM=SUM(WT(1:ND))
	IF(WT_SUM .LE. 0.0_LDP)RETURN
	WT(1:ND)=WT(1:ND)/WT_SUM
!
	DO I=1,ND
	  IF(WT(I) .LE. 0.0_LDP)CYCLE
	  BETA=MAX(ABS(V(I))/C_KMS,MIN_BETA)
	  IF(EW .GT. 0.0_LDP)THEN
	    LAM_LO=LAM0*(1.0_LDP-BETA)
	    LAM_HI=LAM0*(1.0_LDP+BETA)
	  ELSE IF(V(I) .GE. 0.0_LDP)THEN
	    LAM_LO=LAM0*(1.0_LDP-BETA)
	    LAM_HI=LAM0
	  ELSE
	    LAM_LO=LAM0
	    LAM_HI=LAM0*(1.0_LDP+BETA)
	  END IF
	  IF(LAM_HI .LE. LAM_LO)CYCLE
!
! OBS_FREQ (1e15 Hz) is descending; NU_HI is the blue edge.
!
	  NU_HI=0.01_LDP*C_KMS/LAM_LO
	  NU_LO=0.01_LDP*C_KMS/LAM_HI
	  IF(NU_HI .LE. OBS_FREQ(N_OBS))CYCLE
	  IF(NU_LO .GE. OBS_FREQ(1))CYCLE
!
! M1 = first index with OBS_FREQ <= NU_HI; M2 = last with >= NU_LO.
!
	  IF(OBS_FREQ(1) .LE. NU_HI)THEN
	    M1=1
	  ELSE
	    LO=1; HI=N_OBS+1
	    DO WHILE(HI-LO .GT. 1)
	      MID=(LO+HI)/2
	      IF(OBS_FREQ(MID) .GT. NU_HI)THEN
	        LO=MID
	      ELSE
	        HI=MID
	      END IF
	    END DO
	    M1=HI
	  END IF
	  IF(OBS_FREQ(N_OBS) .GE. NU_LO)THEN
	    M2=N_OBS
	  ELSE
	    LO=0; HI=N_OBS
	    DO WHILE(HI-LO .GT. 1)
	      MID=(LO+HI)/2
	      IF(OBS_FREQ(MID) .GE. NU_LO)THEN
	        LO=MID
	      ELSE
	        HI=MID
	      END IF
	    END DO
	    M2=LO
	  END IF
	  IF(M2 .LT. M1)CYCLE
!
	  DEP_DENSITY=WT(I)*EW*CONT_INT/(LAM_HI-LAM_LO)
	  DO M=M1,M2
	    SOB_DEP(M)=SOB_DEP(M)+DEP_DENSITY
	  END DO
	END DO
!
	RETURN
	END
