!
! Main subroutine to examine model output from CMFGEN. Routine is called by
! DISPGEN. All Model and Atomic data is primarily read in by DISPGEN, and
! is contained in the module MOD_DISP
!
	SUBROUTINE MAINGEN(RMDOT,LUM,
	1                    ND,NP,NC,
	1                    N_MAX,ND_MAX,NC_MAX,NP_MAX,
	1                    N_LINE_MAX,N_PLT_MAX)
	USE SET_KIND_MODULE
	USE MOD_DISP
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE MOD_LEV_DIS_BLK
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered  05-Oct-2026 : Reject GRED axis inversion on non-monotonic ED.
! Altered  05-Oct-2026 : Validate log/power arguments in V diagnostics.
! Altered  05-Oct-2026 : Validate V and SIGMA in the SIGMA diagnostic.
! Altered  05-Oct-2026 : Reject GRV velocity inversion on non-monotonic V.
! Altered  29-Jul-2024 : Added TAUN opton (plots tau (Sobolev approximation) as a function of level #.
! Altered  03-Dec-2023 : GF option added (?)
! Altered  16-Jun-2023 : Changed angular quadrature weight routines to V2.
! Altered  19-May-2023 : Fixed bug with COLL option -- size of work vectors TA, TB, and TC increased.
! Altered  22-Feb-2023 : Minor changes to some plotting option. Variable CURVE_LAB added.
!          20-Jan-2023 : ABS_MEAN added.
! Altered  06-Sep-2022 : Changes to LINE_ID output (29-Aug-2022).
! Altered  09-Aug-2022 : Improved handling of FLUX_DEFICIT line id's (PLNID option).
! Altered  17-Nov-2021 : Copied from OSIRIS.
!                        Fixed call to SOBEW_GRAD_V2.
!                        Now call DP_CURVE_LAB (change done initially on OSIRIS).
! Altered  09-Sep-2021 : Now call FQCOMP_IBC_V2 and JEAU_IBC_V2.
! Altered  10-Dec-2019 : CLUMP is now allowed for when computing the recombination rate.
!                          XLOGRM1(LOG10(R(I)/R(ND)-1.0D0) added as option
!                          FLUX mean opacity added to Y options (FLUX) -- Also ES -- mass abs. coef.
!                          Option added to plot B(T) and dB/dT.
!                          Can now output a larger/fine grid with WRC and WRL options.
!                            (for use with MC calculations).
!
! Altered  01-Jul-2018 : Added species/ion options to TAUSOB.
! Altered  21-Mar-2018 : DERAD is now corrected for clumping.
! Altered  17-Mar-2018 : Adjusted to be compatible with OSIRIS version.
!                          On 14-MAr-2018 option to reverse M integration was added.
!                          WR_LINE keyword added to LTAU option.
! Altered  01-MAr-2016 : Adjusted PLNID -- only compute continuum opacity every 1009 km/s.
!                          (Results in considerable speed up [1-Feb-2016]).
! Altered  09-Jun-2015 : Added PLT_PRFS to call.
!                          Altered DTDP option.
! Altered  07-Jul-2011 : Included commands to look at f, and Mdot in "shell" models.
! Altered  15-Mar-2011 :  Section for plotting photoionization cross-sections
!                           removed to subroutine. Can now plot all ground-state
!                           cross section for a species.
!                         TCMF option installed to use GREY T from CMFGEN calculation.
! Altered  13-Jun-2010 :  Rayleigh scattering included in OPACITIES.INC.
!                         ESEC replaced by CHI_SCAT in may locations.
!                         Still some issues to be decided (i.e., Raleigh scattering to opacity file?)
!  	                  Changed to SOBEW_GRAD_V2 -- changed call. Adapted to work for very low opacities.
!                         Now include variable INNER_BND_METH.
!
! Altered  15-Aug-2003 :  QF option installed. Allows ion column densities to be compared.
!                         Option added to PLTPHOT to allow photoioization cross-section
!                           to be compted at a particular electron density. This allows the
!                           effects of level dissolution to be seen.
! Altered  28-Mar-2003 :  MAX_ION replaced bu NUM_IONS everywhere.
! Altered  01-Jub-1998 :  MOD_LEV_DIS_BLK replaces  include file.
!                         REXT_COEF changed to REXT_COEF_V2
! altered  12/16/96  DLM  Changed USRINP call to USR_OPTION.  USR_OPTION
!                         is a generic subroutine call which is defined in
!                         MOD_USR_OPTION.  The interface in MOD_USR_OPTION
!                         determines which usr_option subroutine to call by
!                         examining the types of variables passed to USR_OPTION.
!
!                         .sve files created by SVE_FILE routines
!
!                         Also cleaned up a few places:
!                           1) Changed all comparisons to input value X to compare
!                              to X(1:?) => needed by new .sve routines.
!                           2) Changed LOGR and LINR options to XLOGR and XLINR
!                              to be consistant => needed by new .sve routines.
!
! altered   2/27/97  DLM  Added box= option to create .box files
!                         Added hidden options for reading of .sve file
!
! 
!
	INTEGER NC,NP,ND
	INTEGER N_MAX,ND_MAX,NC_MAX,NP_MAX
	INTEGER N_LINE_MAX,N_PLT_MAX
!
	REAL(KIND=LDP) LUM,RMDOT
!
	INTEGER, PARAMETER :: NLF_MAX=1001
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
! 
!
	REAL(KIND=LDP) EINA
	REAL(KIND=LDP) OSCIL
	REAL(KIND=LDP) GUPDIE
	REAL(KIND=LDP) GLDGU
	REAL(KIND=LDP) EDGEDIE
	REAL(KIND=LDP) GLOW
	REAL(KIND=LDP) NUST(ND)
! 
!
! Arrays containing LINE frequencies in numerical order.
!
	REAL(KIND=LDP) VEC_FREQ(N_LINE_MAX)
	REAL(KIND=LDP) VEC_OSCIL(N_LINE_MAX)
	REAL(KIND=LDP) VEC_EINA(N_LINE_MAX)
	REAL(KIND=LDP) VEC_DP_WRK(N_LINE_MAX)
	INTEGER VEC_INDX(N_LINE_MAX)
	INTEGER VEC_ION_INDX(N_LINE_MAX)
	INTEGER VEC_MNL_F(N_LINE_MAX)
	INTEGER VEC_MNUP_F(N_LINE_MAX)
	INTEGER VEC_INT_WRK(N_LINE_MAX)
	CHARACTER*6 VEC_SPEC(N_LINE_MAX)
	CHARACTER*6 VEC_CHAR_WRK(N_LINE_MAX)
	CHARACTER*6 VEC_TRANS_TYPE(N_LINE_MAX)
	CHARACTER*80, ALLOCATABLE :: VEC_TRANS_NAME(:)
	CHARACTER*80, ALLOCATABLE :: TMP_TRANS_NAME(:)
	INTEGER N_LINE_FREQ
	INTEGER MNL_F
	INTEGER MNUP_F
!
! Angular quadrature weights.
!
	REAL(KIND=LDP) P(NP)				!Impact parameters
	REAL(KIND=LDP) JQW(ND,NP),HMIDQW(ND,NP)
	REAL(KIND=LDP) KQW(ND,NP),NMIDQW(ND,NP)
	REAL(KIND=LDP) HQW(ND,NP)
!
! Note that PFDOP contains the lineprofile in doppler units, whilst
! PF is use for the line profile in frequency units.
!
	REAL(KIND=LDP) LINE_PRO(NLF_MAX),LFQW(NLF_MAX),PF(NLF_MAX)
	REAL(KIND=LDP) ERF(NLF_MAX),PFDOP(NLF_MAX)
	REAL(KIND=LDP) JNU(ND,NLF_MAX+1),HNU(ND,NLF_MAX+1)
	REAL(KIND=LDP) FEDD(ND,NLF_MAX+1),GEDD(ND,NLF_MAX+1)
	REAL(KIND=LDP) HBC_VEC(3,NLF_MAX+1),NBC_VEC(3,NLF_MAX+1)
	REAL(KIND=LDP) INBC_VEC(NLF_MAX+1)
	REAL(KIND=LDP) TDOP,VTURB
!
	REAL(KIND=LDP) JBAR(ND),ZNET(ND)
	REAL(KIND=LDP) CHIROSS(ND),TAUROSS(ND)
!
! Generalized work vectors.
!
	REAL(KIND=LDP) TA(MAX(NP_MAX,N_MAX))
	REAL(KIND=LDP) TB(MAX(NP_MAX,N_MAX))
	REAL(KIND=LDP) TC(MAX(NP_MAX,N_MAX))
!
! Variables required to compute TGREY and for interpolations.
!
	REAL(KIND=LDP) JQWEXT(NP_MAX,NP_MAX),KQWEXT(NP_MAX,NP_MAX),PEXT(NP_MAX)
	REAL(KIND=LDP) Z(NP_MAX),DTAU(NP_MAX),XM(NP_MAX),RJ(NP_MAX)
	REAL(KIND=LDP) CHI(NP_MAX),REXT(NP_MAX),dCHIdr(NP_MAX)
	REAL(KIND=LDP) INBC,HBC,HBCNEW,NBC,FA(NP_MAX),GAM(NP_MAX),GAMH(NP_MAX)
	REAL(KIND=LDP) VEXT(NP_MAX),TEXT(NP_MAX),SIGMAEXT(NP_MAX)
	REAL(KIND=LDP) MASS_DENSITYEXT(NP_MAX),CLUMP_FACEXT(NP_MAX)
	REAL(KIND=LDP) HBC_J,HBC_S
!
	REAL(KIND=LDP) RJEXT(NP_MAX),FEXT(NP_MAX),Q(NP_MAX),FOLD(NP_MAX)
	REAL(KIND=LDP) CHIEXT(NP_MAX),ETAEXT(NP_MAX),ESECEXT(NP_MAX)
	REAL(KIND=LDP) CHILEXT(NP_MAX),ETALEXT(NP_MAX)
	REAL(KIND=LDP) SOURCEEXT(NP_MAX)
	REAL(KIND=LDP) ZETAEXT(NP_MAX),THETAEXT(NP_MAX)
	REAL(KIND=LDP) COEF(0:3,NP_MAX)
	REAL(KIND=LDP) TGREY(NP_MAX)
	INTEGER GRID(NP_MAX),INDX(NP_MAX)
	LOGICAL INACCURATE,REXT_COMPUTED,GREY_COMP
	LOGICAL GREY_WITH_V_TERMS
	LOGICAL GREY_REL
	LOGICAL PLANE_PARALLEL_NOV
	LOGICAL FLUSH_FILE
	LOGICAL OVERWRITE
!
	REAL(KIND=LDP) XV(N_PLT_MAX),XV_SAV(N_PLT_MAX),XNU(N_PLT_MAX)
	REAL(KIND=LDP) YV(N_PLT_MAX),ZV(N_PLT_MAX),WV(N_PLT_MAX)
!
! 
!
! Collisional matrices.
!
	REAL(KIND=LDP) OMEGA_F(N_MAX,N_MAX)
	REAL(KIND=LDP) dln_OMEGA_F_dlnT(N_MAX,N_MAX)
	REAL(KIND=LDP) OMEGA_S(N_MAX,N_MAX)
	REAL(KIND=LDP) dln_OMEGA_S_dlnT(N_MAX,N_MAX)
	EXTERNAL OMEGA_GEN_V3
!
	CHARACTER(LEN=120) COL_RATE_FILE
	REAL(KIND=LDP), ALLOCATABLE :: COL(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: dCOL(:,:,:)
!
	REAL(KIND=LDP) UNIT_VEC(N_MAX)
	REAL(KIND=LDP) ZERO_VEC(N_MAX)
!
	REAL(KIND=LDP) DOP_PRO
	REAL(KIND=LDP) S15ADF,XCROSS_V2
	EXTERNAL JTRPWGT_V2,HTRPWGT_V2,KTRPWGT_V2,NTRPWGT_V2
	EXTERNAL JWEIGHT_V2,HWEIGHT_V2,KWEIGHT_V2,NWEIGHT_V2,XCROSS
!
! Photoionization cross-section routines.
!
	EXTERNAL SUB_PHOT_GEN
!
	CHARACTER*30 UC
	EXTERNAL UC
	REAL(KIND=LDP) KEV_TO_HZ,ANG_TO_HZ
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) C_CMS
	REAL(KIND=LDP) C_KMS
!
	CHARACTER(LEN=6) METHOD
	CHARACTER(LEN=10) TYPE_ATM
	CHARACTER(LEN=10) INNER_BND_METH
!
	REAL(KIND=LDP), ALLOCATABLE :: CHI_PAR(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_PAR(:,:)
	REAL(KIND=LDP) YMAPV(MAX_ION)			!Defined so use instead of NUM_IONS
	REAL(KIND=LDP) XMAPV(ND)
	LOGICAL LAST_NON_ZERO
	CHARACTER*10 LOC_ION_ID
!
	REAL(KIND=LDP), ALLOCATABLE :: CHI_LAM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_LAM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_TOT_LAM(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_TOT_LAM(:)
	REAL(KIND=LDP), ALLOCATABLE :: LAM_VEC(:)
	INTEGER NLAM
	INTEGER NION
!
! Local variables
!
	REAL(KIND=LDP) FB(ND,ND),WM(ND,ND)
	REAL(KIND=LDP) GB(ND),U(ND),VB(ND),VC(ND),CHIL(ND),ETAL(ND)
	REAL(KIND=LDP) SOURCE(ND),TCHI(ND),ZETA(ND),THETA(ND)
	REAL(KIND=LDP) ETA(ND),ETA_WITH_ES(ND)
	REAL(KIND=LDP) ESEC(ND),EMHNUKT(ND),VT(ND)
	REAL(KIND=LDP) CHI_SCAT(ND),CHI_RAY(ND)
	REAL(KIND=LDP) AV(ND),CV(ND)
	REAL(KIND=LDP) FORCE_MULT(ND_MAX)
	LOGICAL DO_DPTH(ND)
!
! Used for computing distrbution of lines with respect to the Sobolev Line
! optical depth.
!
	INTEGER NBINS
	REAL(KIND=LDP) DELTA_TAU
	REAL(KIND=LDP) TAU_BEG
	REAL(KIND=LDP) TAU_MIN
	REAL(KIND=LDP) TAU_CONSTANT
	LOGICAL WEIGHT_NV
	LOGICAL MEAN_TAU
	LOGICAL RADIAL_TAU
	LOGICAL RAY_TAU
!
	INTEGER I,J,K,L
	INTEGER IMIN,IMAX
	INTEGER ISPEC,ID
	INTEGER NLF,ML
	INTEGER IOS,NFREQ,R_INDX
	INTEGER NL,NUP
	INTEGER NEXT_LOC
	INTEGER PHOT_ID
	INTEGER DPTH_INDX
	INTEGER LINE_INDX
!
	INTEGER, PARAMETER :: ITAU_GRT_LIM=6
	INTEGER TAU_GRT_LOGX(-ITAU_GRT_LIM:ITAU_GRT_LIM)
!
	REAL(KIND=LDP) DTDR,DBB,S1,IC
	REAL(KIND=LDP) HFLUX
	REAL(KIND=LDP) EXC_EN,EDGE_FREQ
	REAL(KIND=LDP) RVAL,TAU_VAL,ED_VAL
	REAL(KIND=LDP) XDIS,YDIS,DIS_CONST
	REAL(KIND=LDP) LAM_ST,LAM_EN,DEL_NU
	REAL(KIND=LDP) FREQ_ST,FREQ_EN
	REAL(KIND=LDP) NU_ST,NU_EN
	REAL(KIND=LDP) FREQ_RES,FREQ_MAX
	REAL(KIND=LDP) T1,T2,T3,T4,T5
	REAL(KIND=LDP) FLUX_DEFICIT
	REAL(KIND=LDP) VSHIFT,VDOP,DEL_V
	REAL(KIND=LDP) EK_EJECTA
	REAL(KIND=LDP) TMP_ED
	REAL(KIND=LDP) TAU_LIM
	REAL(KIND=LDP) DEPTH_LIM
	REAL(KIND=LDP) TEMP,TSTAR,NEW_RSTAR,NEW_VSTAR
	REAL(KIND=LDP) VSM_DIE_KMS
	REAL(KIND=LDP) DIST_KPC
!
	INTEGER NPINS,NCX,NDX,NPX
	INTEGER ND_TMP
	REAL(KIND=LDP) FREQ,FL,FL_SAVE
	REAL(KIND=LDP) TAU_SOB
	REAL(KIND=LDP)  TMP_GION
	EQUIVALENCE (FREQ,FL)
!
!
	REAL(KIND=LDP) SCLHT,VCORE,VPHOT,VINF1,V_BETA1,V_EPS1
	REAL(KIND=LDP) VINF2,V_BETA2,V_EPS2
!
	CHARACTER(LEN=20) FILE_FORMAT
	CHARACTER(LEN=20) MOD_NAME
	CHARACTER(LEN=30) TRANS_NAME
	CHARACTER(LEN=100), SAVE :: FILE_FOR_WRL=' '
	LOGICAL FILE_PRES
!
! Storage for constants used for evaluating the level disolution.
! Required by SUP_TO_FULL and LTE_POP_WLD.
!
	LOGICAL LEVEL_DISSOLUTION
!
	CHARACTER*80 NAME,XAXIS,YAXIS,XAXSAV,CURVE_LAB
	CHARACTER(LEN=200) TITLE(10)
!
	INTEGER, PARAMETER :: NSC=31
	COMMON/TOPBORD/ SCED(NSC),XED(NSC),NXED,TOPLABEL
	REAL(KIND=LDP) SCED,XED
	INTEGER NXED
	CHARACTER(LEN=30) TOPLABEL
	DATA SCED/2.0_LDP,2.5_LDP,3.0_LDP,3.5_LDP,4.0_LDP,4.5_LDP,5.0_LDP,5.5_LDP,6.0_LDP,6.5_LDP,7.0_LDP,7.5_LDP,8.0_LDP,8.5_LDP,9.0_LDP,9.5_LDP,10.0_LDP,
	1         10.5_LDP,11.0_LDP,11.5_LDP,12.0_LDP,12.5_LDP,13.0_LDP,13.5_LDP,14.0_LDP,14.5_LDP,15,15.5_LDP,
	1         16,16.5_LDP,17.0_LDP/
!
	INTEGER LEV(10)
	INTEGER IDEPTH(10)
	LOGICAL FLAG,LINV,TRAPFORJ,JONS,JONLY,IN_R_SUN,SONLY
	LOGICAL ELEC,DIF,SCALE,THICK,NORAD,ROSS,INC_RAY_SCAT
	LOGICAL CHECK_FOR_NAN,TMP_LOGICAL,DOWN_RATE
	LOGICAL SPEC_FRAC,RADIAL
	LOGICAL LST_DEPTH_ONLY
	LOGICAL DIE_REG,DIE_WI
	LOGICAL VALID_VALUE
	LOGICAL FOUND
	LOGICAL NEW_FORMAT,NEW_FILE
	LOGICAL NEW_FREQ
	LOGICAL DO_TAU
	LOGICAL DO_KAP
	LOGICAL DONE_LINE
	LOGICAL LINE_STRENGTH
	LOGICAL PLT_J,PLT_H,PLT_LF,PLT_FM
!
	REAL(KIND=LDP) RED_EXT
	REAL(KIND=LDP) VDOP_FG_FRAC
	REAL(KIND=LDP) VDOP_MOM_FRAC
	CHARACTER*30 TWO_PHOT_OPTION
	CHARACTER*30 FREQ_INPUT
	CHARACTER*10 SOL_OPT
	CHARACTER*10 FG_SOL_OPT
	CHARACTER*10 N_TYPE
	CHARACTER*10 LAST_COL_SPECIES
!
	LOGICAL KEV_INPUT,HZ_INPUT,ANG_INPUT
	LOGICAL LINE_BL,FULL_ES,EDDC,HAM,SKIPEW
	LOGICAL THK_CONT,THK_LINE,LIN_DC,LINX,LINY
	LOGICAL FIRST_RATE
	LOGICAL WR_LINE
	LOGICAL L_TRUE,L_FALSE
	DATA L_TRUE/.TRUE./
	DATA L_FALSE/.FALSE./
!
	LOGICAL XRAYS
	REAL(KIND=LDP) FILL_FAC_XRAYS,T_SHOCK,V_SHOCK
	REAL(KIND=LDP) FILL_VEC_SQ(ND)
!
! Equivalent width and flux variables.
!
	INTEGER CNT
	REAL(KIND=LDP) ERR(ND),EWACC
	REAL(KIND=LDP) EW,EWOLD,MOMEW,CONT_INT,MOMCONT_INT
!
	REAL(KIND=LDP) AMASS
!
	CHARACTER MAIN_OPT_STR*80
	CHARACTER X*20
	CHARACTER XOPT*10
	CHARACTER XSPEC*10
!
	CHARACTER TYPE*10
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=10) COL_OPT
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=200) PWD
!
	CHARACTER FMT*80
	CHARACTER*120 DEFAULT
	CHARACTER*120 TAUSOB_DEF
	CHARACTER*120 DESCRIPTION
	CHARACTER*10 XAXIS_OPT
	CHARACTER(LEN=1) FIRST_COLR
!
! External functions
!
	EXTERNAL WR_PWD
	CHARACTER*120 WR_PWD
!
	REAL(KIND=LDP) GFF,GBF,LAMVACAIR,SPEED_OF_LIGHT,BOLTZMANN_CONSTANT
	REAL(KIND=LDP) FUN_PI,SECS_IN_YEAR,MASS_SUN,LUM_SUN,ATOMIC_MASS_UNIT
	REAL(KIND=LDP) ASTRONOMICAL_UNIT,RAD_SUN,TEFF_SUN,STEFAN_BOLTZ
	INTEGER GET_INDX_DP
	EXTERNAL GFF,GBF,LAMVACAIR,SPEED_OF_LIGHT,GET_INDX_DP
	EXTERNAL FUN_PI,SECS_IN_YEAR,MASS_SUN,LUM_SUN,ATOMIC_MASS_UNIT
	EXTERNAL ASTRONOMICAL_UNIT,RAD_SUN,TEFF_SUN,BOLTZMANN_CONSTANT,STEFAN_BOLTZ
!
	LOGICAL PRESENT			!indicates whether species is linked
	LOGICAL EQUAL
!
	INTEGER, PARAMETER :: LU_IN=10		!File input (open/close)
	INTEGER, PARAMETER :: LU_OUT=11		!File input (open/close)
	INTEGER, PARAMETER :: LU_LOG=7
	INTEGER, PARAMETER :: LU_CROSS=15
	INTEGER, PARAMETER :: LU_NET=17
	INTEGER, PARAMETER :: LU_REC=18
	INTEGER, PARAMETER :: LU_COL=20		!Collisonal
	INTEGER, PARAMETER :: T_IN=5                 !Terminal input
	INTEGER, PARAMETER :: T_OUT=6                !Terminal output
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ITEN=10
	INTEGER, PARAMETER :: ITHREE=3
!
	REAL(KIND=LDP), PARAMETER :: RZERO=0.0_LDP
	REAL(KIND=LDP), PARAMETER :: RONE=1.0_LDP
	REAL(KIND=LDP), PARAMETER :: RTWO=2.0_LDP
	REAL(KIND=LDP), PARAMETER :: RTHREE=3.0_LDP
!
! 
!
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	LAST_COL_SPECIES=' '
	METHOD='LOGLOG'
	FIRST_COLR='T'
	DPTH_INDX=ND/2
	COL_OPT='CL'
	FILE_FORMAT='OLD'
	IC=RZERO
!
! NB: ZERO_FLUX condition is equivalent to DIF=.TRUE. and DBB=0.
!     Thus, for routine compatibility, we set DBB=0 when DIF is false.
!
	IF(V(1) .GT. 10000.0_LDP .OR. V(ND) .GT. 100.0_LDP)THEN
	  INNER_BND_METH='ZERO_FLUX'		!Default for SN model
	  DIF=.FALSE.
	  WRITE(6,*)'Inner boudary option zet to zero_flux'
	ELSE
	  INNER_BND_METH='DIFFUSION'		!Default for normal stars
	  DIF=.TRUE.
	  WRITE(6,*)'Inner boudary option zet to diffusion approximation'
	END IF
	TRAPFORJ=.TRUE.
	TSTAR=T(ND)
!
	TYPE_ATM='P'
	T1=LOG(MASS_DENSITY(5)/MASS_DENSITY(1))/LOG(R(1)/R(5))
	WRITE(6,'(A)')RED_PEN
	IF(T1 .LT. 2)THEN
	  WRITE(6,'(A)')' Power for density law variation at outer boudary is shallower than r^{-n} where n=2'
	  WRITE(6,'(A)')' Shown below is the power based on depths 1 and I'
	  WRITE(6,*)' '
	  DO I=2,MIN(15,ND)
	    WRITE(6,'(X,I4,ES12.3)')I,LOG(MASS_DENSITY(I)/MASS_DENSITY(1))/LOG(R(1)/R(I))
	  END DO
	  CALL USR_OPTION(T1,'POW','2','Power for rhj propto r^{-n}')
	  WRITE(6,*)' '
	END IF
!
	IF(T1 .LT. 100)THEN
	  WRITE(TYPE_ATM(2:6),'(F5.2)')T1
	  WRITE(6,'(A)')' For optical depth calculations we will assume a power density distribution'
	  WRITE(6,'(A,F5.2)')' at the outer boundary. Density exponent for atmosphere is ',T1
	ELSE
	  WRITE(TYPE_ATM(2:9),'(ES8.2)')T1
	  WRITE(6,'(A)')' For optical depth calculations we will assume a power density distribution'
	  WRITE(6,'(A,ES8.2)')' at the outer boundary. Density exponent for atmosphere is ',T1
	END IF
	WRITE(6,'(A)')DEF_PEN
!
! If the grey temperature is available in CMFGEN, we don't need to (although we can with the
! GREY option) compute it.
!
	GREY_COMP=.FALSE.
	IF(CMFGEN_TGREY(ND) .EQ. RZERO)THEN
	  GREY_COMP=.FALSE.
	ELSE
	  GREY_COMP=.TRUE.
	  TGREY(1:ND)=CMFGEN_TGREY(1:ND)
	END IF
	FIRST_RATE=.TRUE.
!
	XRAYS=.TRUE.                    !Changed def .FALSE.
	FILL_FAC_XRAYS=1.0E-100_LDP
	T_SHOCK=300.0_LDP
	V_SHOCK=200.0_LDP
	FILL_VEC_SQ(:)=RZERO
!
	LEVEL_DISSOLUTION=.TRUE.
	PI=FUN_PI()
	DIST_KPC=RONE
	VSM_DIE_KMS=3000.0_LDP
	REXT_COMPUTED=.FALSE.
!
	LST_DEPTH_ONLY=.FALSE.
!
	LAM_ST=RZERO
	LAM_EN=1.0E+06_LDP
	TAUSOB_DEF='1000'
!
	ZERO_VEC(:)=RZERO
	UNIT_VEC(:)=RONE
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	KEV_TO_HZ=0.241838E+03_LDP
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0E-07_LDP  	!10^8/10^15
	C_CMS=SPEED_OF_LIGHT()
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
!
	ANG_INPUT=.TRUE.
	HZ_INPUT=.FALSE.
	KEV_INPUT=.FALSE.
	FREQ_INPUT='Units are Angstroms'
!
! DTDR is used in the diffusion approximation inner boundary condition.
! In the MAIN code it is fixxd by the ROSSELAND mean opacity, and the
! LUMINOSITY.
!
	DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND))
!
	THK_CONT=.TRUE.
	THK_LINE=.TRUE.
	LINE_BL=.TRUE.
!
! Take log of R (default X axis).
!
	DO I=1,ND
	  XV(I)=LOG10(R(I)/R(ND))
	END DO
	XAXIS='Log(r/R\d*\u)'
	XAXSAV=XAXIS
	DEFAULT=WR_PWD()
	CALL RD_STRING(NAME,'Title',DEFAULT,
	1    'Title to be placed on all graphs')
!
! Allocate arrays:
!
	IOS=0
	ALLOCATE (CHI_PAR(ND,NUM_IONS),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ETA_PAR(ND,NUM_IONS),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error in MAINGEN'
	  WRITE(T_OUT,*)'Unable to allocate CHI_PAR or  ETA_PAR'
	  STOP
	END IF
!
! 
!
! Determine the population of the ground state, summing up over J levels that
! are within 2000 cm^-1 of each other. These are used to evaluate the
! free-free rates, and also provide a check on effective recombination
! rates.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    ATM(ID)%DXzV(1:ND)=ATM(ID)%DXzV_F(1:ND)
	    IF(ATM(ID+1)%XzV_PRES)THEN
	      DO J=2,ATM(ID+1)%NXzV_F
	        IF((ATM(ID+1)%EDGEXzV_F(J)-ATM(ID+1)%EDGEXzV_F(1))
	1          *1.0E+015_LDP/C_CMS .LT. 2000)THEN
	          DO I=1,ND
	            ATM(ID)%DXzV(I)=ATM(ID)%DXzV(I)+ATM(ID+1)%XzV_F(J,I)
	          END DO
	        END IF
	      END DO
	    END IF
	  END IF
	END DO
!
! 
!
! Determine impact parameters. Calculate all required Angular quadrature
! weights.
!
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
!
! The HMIDQW and NMIDQW qw are defined at the MIDPOINTS of the mesh, hence
! the option TRUE.
!
	IF(TRAPFORJ)THEN
	  CALL GENANGQW_V2(JQW,R,P,NC,ND,NP,JTRPWGT_V2,.FALSE.)
	  CALL GENANGQW_V2(HQW,R,P,NC,ND,NP,HTRPWGT_V2,.FALSE.)
	  CALL GENANGQW_V2(KQW,R,P,NC,ND,NP,KTRPWGT_V2,.FALSE.)
	  CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HTRPWGT_V2,.TRUE.)
	  CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NTRPWGT_V2,.TRUE.)
	ELSE
	  CALL GENANGQW_V2(JQW,R,P,NC,ND,NP,JWEIGHT_V2,.FALSE.)
	  CALL GENANGQW_V2(HQW,R,P,NC,ND,NP,HWEIGHT_V2,.FALSE.)
	  CALL GENANGQW_V2(KQW,R,P,NC,ND,NP,KWEIGHT_V2,.FALSE.)
	  CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HWEIGHT_V2,.TRUE.)
	  CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NWEIGHT_V2,.TRUE.)
	END IF
!
! 
! Compute vector constants for evaluating the level disolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD. NB: POPION was read in from
! the RVTJ file.
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,LEVEL_DISSOLUTION,ND)
!
! Set 2-photon data with current atomic models and populations.
!
        DO ID=1,NUM_IONS-1
          I=ID
          CALL SET_TWO_PHOT_DISP_V3(ION_ID(ID), ID,
	1       ATM(ID)%XzVLTE_F,        ATM(ID)%NXzV_F,
	1       RONE,                    ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,       ATM(ID)%GXzV_F,
	1       RONE,                    ATM(ID)%NXzV_F,  ND,
	1       ATM(ID)%ZXzV,            IONE,   ATM(ID)%XzV_PRES)
        END DO
!
!
! Compute LTE populations
!
	INCLUDE 'EVAL_LTE_FULL.INC'
!
! 
!
! Set up lines that will be treated with the continuum calculation.
! This section of code is also used by the code treating purely lines
! (either single transition Sobolev or CMF, or overlapping Sobolev).
!
! To define the line transitions we need to operate on the FULL atom models.
! We thus perform separate loops for each species.
!
	ML=0			!Initialize line counter.
	ALLOCATE (VEC_TRANS_NAME(N_LINE_MAX))
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNL_F=1,ATM(ID)%NXzV_F-1
	      DO MNUP_F=MNL_F+1,ATM(ID)%NXzV_F
	        IF(ATM(ID)%AXzV_F(MNL_F,MNUP_F) .NE. 0)THEN
	          ML=ML+1
	          IF(ML .GT. N_LINE_MAX)THEN
	            WRITE(T_OUT,*)'Error in MAINGEN'
	            WRITE(T_OUT,*)'N_LINE_MAX is too small'
	            STOP
	          END IF
	          VEC_FREQ(ML)=ATM(ID)%EDGEXzV_F(MNL_F)-ATM(ID)%EDGEXzV_F(MNUP_F)
	          VEC_SPEC(ML)=ION_ID(ID)
	          VEC_MNL_F(ML)=MNL_F
	          VEC_MNUP_F(ML)=MNUP_F
	          VEC_ION_INDX(ML)=ID
	          VEC_OSCIL(ML)=ATM(ID)%AXzV_F(MNL_F,MNUP_F)
	          VEC_EINA(ML)=ATM(ID)%AXzV_F(MNUP_F,MNL_F)
	          VEC_TRANS_NAME(ML)=TRIM(VEC_SPEC(ML))//
	1           '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1           TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	        END IF
	      END DO
	    END DO
	  END IF
	END DO
	WRITE(6,*)'Total number of line transition is',ML
!
!
! 
!
!
! Get lines and arange in numerically decreasing frequency. This will
! allow us to consider line overlap, and to include lines with continuum
! frequencies so that the can be handled automatically.
!
	N_LINE_FREQ=ML
!
	CALL INDEXX(N_LINE_FREQ,VEC_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_ION_INDX,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	ALLOCATE (TMP_TRANS_NAME(N_LINE_FREQ))
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_NAME,VEC_INDX,TMP_TRANS_NAME)
	DEALLOCATE (TMP_TRANS_NAME)
!
! 
!
! Entering main secion of the subroutine which controls which OPTIONS
! are executed. Above statements are initialization only.
!
! This message will only be printed once
!
	WRITE(T_OUT,*)
	WRITE(T_OUT,"(8X,A)")'(default is to write file '//
	1    'MAIN_OPT_STR.sve)'
	WRITE(T_OUT,"(8X,A)")'(append sve=filename to '//
	1    'write a new .sve file)'
	WRITE(T_OUT,"(8X,A)")'(box=filename to '//
	1    'write a .box file containing several .sve files)'
	WRITE(T_OUT,"(8X,A)")'(.filename to read .sve file)'
	WRITE(T_OUT,"(8X,A)")'(#filename to read .box file)'
	WRITE(T_OUT,*)
!
!
! This call resets the .sve algorithm.  Specifically it sets the next
! input answer to be a main option, and all subsequent inputs to be
! sub-options.
!
 3	CALL SVE_FILE('RESET')
!
	DESCRIPTION=' '				!As obvious.
	MAIN_OPT_STR='  '
	DEFAULT='GR'
	CALL USR_OPTION(MAIN_OPT_STR,'OPTION',DEFAULT,DESCRIPTION)
!
!   If the main option begins with a '.', a previously
!   written .sve file is read.
!
!   If the main option begins with a '#', a previously
!   written .box file is read.
!
!   If sve= is apended to the end of this main option, a new .sve file
!   is opened with the given name and the main option and all subsequent
!   sub-options are written to this file.
!
!   If box= is input then a .box file is created, which contains the name
!   of several .sve files to process.
!
!   If only a main option is given, the option and subsequent sub-options
!   are saved in a file called 'main option.sve'.  All following main
!   options are saved in separate files.
!
! Remove variable changes from main option.
!
	I=INDEX(MAIN_OPT_STR,'(')
	IF(I .EQ. 0)THEN
	  X=UC(TRIM(MAIN_OPT_STR))
	ELSE
	  X=UC(MAIN_OPT_STR(1:I-1))	!Remove line variables.
	END IF
!
! Remove possile file names etc.
!
	I=INDEX(X,' ')
	IF(I .EQ. 0)THEN
	  X=UC(TRIM(X))
	ELSE
	  X=UC(X(1:I-1))	!Remove file names
	END IF
!
! NB: X contains the full option
!     XOPT contains the main part of the option with the species identifier
!     XSPEC is the species identifier (eg. HE2)
!
! _ is presently used to separate the OPTION (e.g. DC) from the species
! (He2). This is the only location it needs ti specified.
!
	I=INDEX(X,'_')
	IF(I .EQ. 0)THEN
	  XSPEC=' '
	  XOPT=X
	ELSE
	  J=INDEX(X,' ')
	  XOPT=X(1:I-1)
	  XSPEC=X(I+1:J)
	END IF
!
! Check whether the species is
!                             (i)  In code
!                             (ii) Available
!
	IF(XSPEC .NE. ' ')THEN
	  PRESENT=.FALSE.
	  IF(XSPEC .EQ. 'ALL')THEN	!Required for SMOOTH option
	    PRESENT=.TRUE.
	  END IF
!
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. UC(SPECIES(ISPEC)))THEN
	      PRESENT=.TRUE.
	      IF(POPDUM(ND,ISPEC) .EQ. 0)THEN
	        WRITE(T_OUT,*)'Error --- this  species is unavailable'
	        GOTO 1
	      END IF
	    END IF
	   END DO
!
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      PRESENT=.TRUE.
	      IF(.NOT. ATM(ID)%XzV_PRES)THEN
	        WRITE(T_OUT,*)'Error --- this  ionization stage is unavailable'
	        GOTO 1
	      END IF
	    END IF
	  END DO
!
	  IF(.NOT. PRESENT)THEN
	    WRITE(T_OUT,*)'Error --- unknown species or ionization stage'
	    WRITE(T_OUT,*)'      --- ',XSPEC
	    WRITE(T_OUT,*)'Species may not be linked into DISPLAY package'
	    GOTO 1
	  END IF
	END IF
!
! 
!
! This IF statment is a preliminary to later work. It is
! used when line opacities are required. Only of use
! for options operating on lines. (First call is input NL and NUP).
!
! Line options.
!
	IF(XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR'
	1          .OR. XOPT .EQ. 'EW'
	1          .OR. XOPT .EQ. 'LAM'
	1          .OR. XOPT .EQ. 'SRCE'
	1          .OR. XOPT .EQ. 'SRCEBB'
	1          .OR. XOPT .EQ. 'SRCEJC'
	1          .OR. XOPT .EQ. 'EP'
	1          .OR. XOPT .EQ. 'CHIL'
	1          .OR. XOPT .EQ. 'TAUL'
	1          .OR. XOPT .EQ. 'TAULIP'
	1          .OR. XOPT .EQ. 'MTAULIP'
	1          .OR. XOPT .EQ. 'LINRC'
	1          .OR. XOPT .EQ. 'WRL'
	1          .OR. XOPT .EQ. 'BETA')THEN
!
	   CALL GET_LINE_INDICES(VEC_FREQ,VEC_MNL_F,VEC_MNUP_F,VEC_SPEC,VEC_TRANS_NAME,
	1                           N_LINE_FREQ,XSPEC,NL,NUP,FLAG)
	   LEV(1)=NL; LEV(2)=NUP
	   IF(.NOT. FLAG)GOTO 1
!
! Compute line opacity and emissivity.
!
	  FLAG=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)))THEN
	      IF(LEV(1) .LE. ATM(ID)%NXzV_F .AND. LEV(2) .LE. ATM(ID)%NXzV_F)THEN
	        FREQ=ATM(ID)%EDGEXzV_F(LEV(1))-ATM(ID)%EDGEXzV_F(LEV(2))
	        OSCIL=ATM(ID)%AXzV_F(LEV(1),LEV(2))
	        DO I=1,ND
	          T1=ATM(ID)%W_XzV_F(LEV(2),I)/ATM(ID)%W_XzV_F(LEV(1),I)
	          GLDGU=ATM(ID)%GXzV_F(LEV(1))/ATM(ID)%GXzV_F(LEV(2))
	          CHIL(I)=OPLIN*ATM(ID)%AXzV_F(LEV(1),LEV(2))*(
	1           T1*ATM(ID)%XzV_F(LEV(1),I)-GLDGU*ATM(ID)%XzV_F(LEV(2),I) )
	          ETAL(I)=EMLIN*FREQ*ATM(ID)%AXzV_F(LEV(2),LEV(1))*ATM(ID)%XzV_F(LEV(2),I)
	        END DO
	        CALL WRITE_LINE(LEV(1),LEV(2),FREQ,ION_ID(ID))
	        AMASS=AT_MASS(SPECIES_LNK(ID))
	        FLAG=.TRUE.
	        T1=ANG_TO_HZ/FREQ
	        IF(T1 .GT. 1000)THEN
	          I=NINT(T1)
	          WRITE(TRANS_NAME,'(I6)')I; TRANS_NAME=ADJUSTL(TRANS_NAME)
	        ELSE
	          WRITE(TRANS_NAME,'(F6.1)')T1; TRANS_NAME=ADJUSTL(TRANS_NAME)
	        END IF
	        TRANS_NAME=TRIM(PLT_ION_ID(ID))//' \gl'//TRIM(TRANS_NAME)
	      ELSE
	        WRITE(T_OUT,*)'Levels outside valid range'
	        WRITE(T_OUT,*)'Maximum level is',ATM(ID)%NXzV_F
	        GOTO 1
	      END IF
	    END IF
	  END DO
!
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)' Invalid population type or species unavailable.'
	    GOTO 1
	  END IF
!
	  IF(OSCIL .EQ. 0)THEN
	    WRITE(T_OUT,*)'Oscillator zero for this transition'
	    WRITE(T_OUT,*)'Rerturning to main loop for new option'
	    GOTO 1
	  END IF
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)'Levels outside permitted range'
	    WRITE(T_OUT,*)'Rerturning to main loop for new option'
	    GOTO 1
	  END IF
!
! Adjust the line opacities and emissivities for the effect of clumping.
!
	  DO I=1,ND
	    CHIL(I)=CHIL(I)*CLUMP_FAC(I)
	    ETAL(I)=ETAL(I)*CLUMP_FAC(I)
	  END DO
!
	END IF
!
! 
!
	IF(XOPT .EQ. 'DIE')THEN
	  CALL USR_OPTION(LEV,1,1,'LEVS',' ','Lower level ')
	  CALL USR_OPTION(EINA,'EINA',' ','Einstein A coefficient')
	  CALL USR_OPTION(GUPDIE,'G_UP',' ','G for Autoionizing level')
	  CALL USR_OPTION(FL,'FL',' ','Transition frequency')
!
	  FLAG=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      EDGEDIE=ATM(ID)%EDGEXzV_F(LEV(1))-FL
	      CALL LTEPOP(NUST,ED,ATM(ID)%DXzV_F,GUPDIE,EDGEDIE,T,
	1                 ATM(ID)%GIONXzV_F,1,ND)
	      GLOW=ATM(ID)%GXzV_F(LEV(1))
	      DO I=1,ND
	        TA(I)=ATM(ID)%XzV_F(LEV(1),I)*ATM(ID)%W_XzV_F(LEV(1),I)
	      END DO
	      FLAG=.TRUE.
	    END IF
	  END DO
!
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)' Invalid Population type'
	    GOTO 1
	  END IF
!
	  GLDGU=GLOW/GUPDIE
	  OSCIL=EINA*EMLIN/( GLDGU*OPLIN*TWOHCSQ*(FL**2) )
	  T1=OSCIL*OPLIN
	  T2=FL*EINA*EMLIN
!
! The effect of clumping is allowed for.
!
	  DO I=1,ND
	    CHIL(I)=T1*( TA(I)-GLDGU*NUST(I) )*CLUMP_FAC(I)
	    ETAL(I)=T2*NUST(I)*CLUMP_FAC(I)
	  END DO
	END IF
!
! 
!
! **************************************************************************
! **************************************************************************
!
! Continuum options requiring an INPUT frequency, and an indication of whether
! electron scattering is to be included in the opacity. Electron scattering
! opacity should be included for all radiation field options.
!
! **************************************************************************
! **************************************************************************
!
	IF( XOPT .EQ. 'OP'
	1          .OR. XOPT .EQ. 'KAPPA'
	1          .OR. XOPT .EQ. 'ETA'
	1          .OR. XOPT .EQ. 'XTAUC'
	1          .OR. XOPT .EQ. 'TAUC'
	1          .OR. XOPT .EQ. 'NEWRG'
	1          .OR. XOPT .EQ. 'DTAUC'
	1          .OR. XOPT .EQ. 'WRC'
	1          .OR. XOPT .EQ. 'JEXT'
	1          .OR. XOPT .EQ. 'CJBB')THEN
	  CALL USR_OPTION(FREQ,'LAM','0.0',FREQ_INPUT)
	  IF(FREQ .EQ. 0)THEN
	    CALL ESOPAC(CHI,ED,ND)		!Electron scattering
	    DO I=1,ND
	       CHI(I)=CHI(I)*CLUMP_FAC(I)
	    END DO
	    ELEC=.TRUE.
	  ELSE
	    IF(KEV_INPUT)THEN
	      FREQ=FREQ*KEV_TO_HZ
	    ELSE IF(ANG_INPUT)THEN
	      FREQ=ANG_TO_HZ/FREQ
	    END IF
	    CALL USR_OPTION(ELEC,'ELEC','T','Include electron scattering?')
	    CALL USR_OPTION(INC_RAY_SCAT,'RAY','F','Include Rayeigh scattering')
	  END IF
	END IF
!
! 
!
! **************************************************************************
! **************************************************************************
!
! Compute continuum opacity for options which require continuum
! opacity.
!
! **************************************************************************
! **************************************************************************
!
	IF(  (XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR'
	1          .OR. XOPT .EQ. 'EW'
	1          .OR. XOPT .EQ. 'SRCE'
	1          .OR. XOPT .EQ. 'EP'
	1          .OR. XOPT .EQ. 'BETA'
	1          .OR. XOPT .EQ. 'ETA'
	1          .OR. XOPT .EQ. 'OP'
	1          .OR. XOPT .EQ. 'KAPPA'
	1          .OR. XOPT .EQ. 'TAUC'
	1          .OR. XOPT .EQ. 'DTAUC'
	1          .OR. XOPT .EQ. 'XTAUC'
	1          .OR. XOPT .EQ. 'NEWRG'
	1          .OR. XOPT .EQ. 'WRC'
	1          .OR. XOPT .EQ. 'WRL'
	1          .OR. XOPT .EQ. 'JEXT'
	1          .OR. XOPT .EQ. 'CJBB'
	1          .OR. XOPT .EQ. 'DIE')
	1          .AND. FL .NE. 0  )THEN
!
! Compute continuum opacity and emissivity at the line frequency.
!
	  INCLUDE 'OPACITIES.INC'
!	  INCLUDE 'HOME:[jdh.cmf.carb.test]OPACITIES.INC'
!
! Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
! and DDBBDT= dB/dTR .
!
	  T1=HDKT*FL/T(ND)
	  T2=RONE-EMHNUKT(ND)
	  DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	  IF(.NOT. DIF)DBB=RZERO
!
! Solve for the continuum radiation field at the line frequency.
!
	  THICK=.TRUE.
	  S1=ETA(1)/CHI(1)
!
! Adjust the opacities and emissivities for the influence of clumping.
!
	  DO I=1,ND
	    ETA(I)=ETA(I)*CLUMP_FAC(I)
	    CHI(I)=CHI(I)*CLUMP_FAC(I)
	    ESEC(I)=ESEC(I)*CLUMP_FAC(I)
	    CHI_RAY(I)=CHI_RAY(I)*CLUMP_FAC(I)
	    CHI_SCAT(I)=CHI_SCAT(I)*CLUMP_FAC(I)
	  END DO
!
	END IF
! 
!
! **************************************************************************
! **************************************************************************
!
! Define reinded radius grid, and evaluate quadrature weights for angular
! integrations.
!
! **************************************************************************
! **************************************************************************
!
	IF(   (XOPT .EQ. 'GREY' .OR. XOPT .EQ. 'JEXT'
	1                         .OR. XOPT .EQ. 'EWEXT'
	1                         .OR. XOPT .EQ. 'INTERP')
	1    .AND. .NOT. REXT_COMPUTED)THEN
	  REXT_COMPUTED=.TRUE.
	  CALL USR_HIDDEN(NPINS,'NPINS','2','0, 1 or 2')
	  NDX=(ND-1)*NPINS+ND
	  I=ND-10		!Parabolic interp for > I
!
! NTERP has been set to ND
!
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDX,R,GRID,ND,
	1                     NPINS,.TRUE.,I,IONE,ND)
	  NCX=NC+10
	  NPX=NDX+NCX-2
	  CALL IMPAR(PEXT,REXT,R(ND),NCX,NDX,NPX)
!
! Note that the T? vectors (here used as dummy variables) must be at least
! NPX long.
!
	  IF(TRAPFORJ)THEN
	    CALL GENANGQW_V2(JQWEXT,REXT,PEXT,NCX,NDX,NPX,JTRPWGT_V2,.FALSE.)
	    CALL GENANGQW_V2(KQWEXT,REXT,PEXT,NCX,NDX,NPX,KTRPWGT_V2,.FALSE.)
	  ELSE
	    CALL GENANGQW_V2(JQWEXT,REXT,PEXT,NCX,NDX,NPX,JWEIGHT_V2,.FALSE.)
	    CALL GENANGQW_V2(KQWEXT,REXT,PEXT,NCX,NDX,NPX,KWEIGHT_V2,.FALSE.)
	  END IF
	END IF
!
! 
!
	IF(  (XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR') .AND. FIRST_RATE)THEN
	  FIRST_RATE=.FALSE.
	  CALL GEN_ASCI_OPEN(LU_NET,'RATES','UNKNOWN',' ',' ',IZERO,IOS)
	END IF
!
	IF(XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR')THEN
!
	  CALL USR_OPTION(NLF,'NLF','25','Number of requency points in line profile?')
	  IF(NLF .GT. NLF_MAX)THEN
	    WRITE(T_OUT,*)'Error - NLF is too large - Maximum values is',
	1	 NLF_MAX
	    GOTO 1
	  END IF
!
! Set frequencies equally spaced in frequency.
!
	  T1=12.0_LDP/(NLF-1)
	  DO I=1,NLF
	    PFDOP(I)=6.0_LDP-(I-1)*T1
	  END DO
!
! ERF is used in computin the Sobolev incident intensity at the
! outer boudary. ERF = int from -inf to x of e(-x^2)/sqrt(pi).
! Note that ERF is not the error function. ERF is related to the
! complementary error function by ERF =-0.5D0 . erfc(X).
! S15ADF is a NAG routine which returns erfc(x).
!
! The incident Sobolev intensity is S[ 1.0-exp(-tau(sob)*ERF) ]
!
	  J=0
	  DO I=1,NLF
	    ERF(I)=-0.5_LDP*S15ADF(PFDOP(I),J)
	  END DO
!
	  CALL USR_OPTION(TDOP,'TDOP','0.0',
	1    'Doppler temperature for Line profile (units 10^4K)?')
	  CALL USR_OPTION(VTURB,'VTURB','10.0',
	1      'Turbulent velcity for Line profile (units km/s)?')
!
	  T1=4.286299E-05_LDP*SQRT( TDOP/AMASS + (VTURB/12.85_LDP)**2 )
	  DO I=1,NLF
	    PF(I)=PFDOP(I)*T1
	  END DO
	  CALL TRAPUNEQ(PF,LFQW,NLF)
	  T1=0.0
	  DO ML=1,NLF
	    LFQW(ML)=LFQW(ML)*FL*1.0E+15_LDP
	    LINE_PRO(ML)=DOP_PRO(PF(ML),FL,TDOP,VTURB,AMASS)
	    T1=T1+LFQW(ML)*LINE_PRO(ML)
	  END DO
!
! Normalize frequency weights.
!
	  DO ML=1,NLF
	    LFQW(ML)=LFQW(ML)/T1
	  END DO
!
! Ensure total opacity is positive.
!
	  J=(NLF+1)/2.0_LDP
	  DO I=1,ND
	    T2=CHI(I)+CHIL(I)*LINE_PRO(J)
	    IF( T2 .LT. 0.0_LDP)THEN
	      T1=3.0E-10_LDP*CHIL(I)*R(I)/V(I)/FL
	      WRITE(T_OUT,
	1      "(' Negative total opacity at d= ',I3,' :  Tausob=',
	1      1P,E12.4,' :   X/Xc=',1P,E12.4)")I,T1,T2/CHI(I)
	      CHIL(I)=-0.8_LDP*CHI(I)/LINE_PRO(J)
	    END IF
	  END DO
	END IF
!
! 
!
! This section compute the CONTINUUM mean intensity using EDDINGTON factors, or
! a full ray-by ray solution.
!
	IF(XOPT .EQ. 'NETR'
	1          .OR. XOPT .EQ. 'MOMR'
	1          .OR. XOPT .EQ. 'SOBR'
	1          .OR. XOPT .EQ. 'SRCEJC'
	1          .OR. XOPT .EQ. 'CJBB'
	1          .OR. XOPT .EQ. 'DIE'
	1          .OR. ( XOPT .EQ. 'EW' .AND. XOPT .NE. 'EWEXT' ))THEN
!
	  EDDC=.TRUE.
	  DEFAULT='T'
	  IF(XOPT .EQ. 'MOMR')THEN
	    EDDC=.TRUE.
	    DEFAULT='T'
	  ENDIF
	  CALL USR_OPTION(EDDC,'EDDC',DEFAULT,'Use Eddington factors to compute Jc ?')
!
	  IF(EDDC)THEN
	    S1=ZETA(1)
	    DO I=1,ND
	      SOURCE(I)=ZETA(I)
	    END DO
	    INACCURATE=.TRUE.
	    HBC=0.99_LDP; HBC_S=RZERO
	    THK_CONT=.FALSE.
	    DO WHILE(INACCURATE)
	      CALL FQCOMP_IBC_V2(TA,TB,TC,XM,DTAU,R,Z,P,Q,FEXT,
	1            SOURCE,CHI,DCHIDR,JQW,KQW,DBB,HBC_J,HBC_S,
	1            INBC,IC,THK_CONT,INNER_BND_METH,NC,ND,NP,METHOD)
	      CALL JFEAU_IBC_V2(TA,TB,TC,DTAU,R,RJ,Q,FEXT,
	1            ZETA,THETA,CHI,DBB,IC,HBC_J,HBC_S,
	1            INBC,THK_CONT,INNER_BND_METH,ND,METHOD)
	      S1=ZETA(1)+THETA(1)*RJ(1)
	      INACCURATE=.FALSE.
	      T1=RZERO
	      DO I=1,ND
	       T1=MAX(ABS(FOLD(I)-FEXT(I)),T1)
	       FOLD(I)=FEXT(I)
	       SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	      END DO
	      IF(T1 .GT. 1.0E-05_LDP)INACCURATE=.TRUE.
	      WRITE(T_OUT,'('' Maximum fractional change is'',1PE11.4)')T1
	    END DO
	    WRITE(6,*)GAM(1),RJ(1)
	    DO I=1,ND
	      WRITE(26,'(ES18.8,6ES14.4)')R(I),ETA(I),CHI(I),THETA(I)*CHI(I),ZETA(I),THETA(I),RJ(I)
	    END DO
!
	  ELSE
!
! Compute J. DBB and S1 have been previously evaluated. VT is used
! for dCHI_dR.
!
	    CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1               ZETA,THETA,CHI,VT,JQW,
	1               THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
	  END IF
!
! Increment emissivity due to electron scatering from the continuum.
! it is assumed that electron scatering from the line does not
! contribute signifcantly to the emissivity.
!
	  DO I=1,ND
	    ETA_WITH_ES(I)=ETA(I)+RJ(I)*CHI_SCAT(I)
	  END DO
!
	END IF
! 
!
! ****************************************************************************
! ****************************************************************************
!
! Begin the individual option section.
!
! *****************************************************************************
! *****************************************************************************
!
! ****************************************************************************
! ****************************************************************************
!
! Set up X axis options -  Log(R/RP),R/RP,LOG(T(rossland)) or Ne .
!
! ****************************************************************************
! ****************************************************************************
!
	IF(XOPT .EQ. 'XLOGR')THEN
	  DO I=1,ND
	    XV(I)=LOG10(R(I)/R(ND))
	  END DO
	  XAXIS='Log(r/R\d*\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XLOGRM1')THEN
	  DO I=1,ND-1
	    XV(I)=LOG10(R(I)/R(ND)-RONE)
	  END DO
	  XV(ND)=XV(ND-1)-1
	  XAXIS='Log(r/R\d*\u)-1'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XLINR')THEN
	  CALL USR_HIDDEN(FLAG,'NORM','T','Normalized X axis (def=T)')
	  IF(FLAG)THEN
	    DO I=1,ND
	      XV(I)=R(I)/R(ND)
	    END DO
	    XAXIS='r/R\d*\u'
	  ELSE
	    T1=RONE
	    XAXIS='r(10\u10 \dcm)'
	    IF(R(ND) .GT. 1.0E+04_LDP)THEN
	      T1=1.0E+04_LDP
	      XAXIS='r(10\u14 \dcm)'
	    END IF
	    DO I=1,ND
	      XV(I)=R(I)/T1
	    END DO
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XAU')THEN
	  FLAG=.FALSE.
	  CALL USR_HIDDEN(FLAG,'LIN','F','Linear axis (def=F)')
	  T1=1.0E+10_LDP/ASTRONOMICAL_UNIT()
	  XV(1:ND)=T1*R(1:ND)
	  IF(FLAG)THEN
	    XAXIS='r(AU)'
	  ELSE
	    XV(1:ND)=LOG10(XV(1:ND))
	    XAXIS='Log r(AU)'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XRSUN')THEN
	  FLAG=.FALSE.
	  CALL USR_HIDDEN(FLAG,'LIN','F','Linear axis (def=F)')
	  T1=1.0E+10_LDP/RAD_SUN()
	  XV(1:ND)=T1*R(1:ND)
	  IF(FLAG)THEN
	    XAXIS='r(R\dsun\u)'
	  ELSE
	    XV(1:ND)=LOG10(XV(1:ND))
	    XAXIS='Log r(R\d\sun\u)'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XARC')THEN
!
! NB: 1 AU corresponds to exactly 1" at 1kpc.
!
	  FLAG=.FALSE.
	  DEFAULT=WR_STRING(DIST_KPC)
	  CALL USR_OPTION(DIST_KPC,'DIST',' ','Distance to star in kpc')
	  CALL USR_HIDDEN(FLAG,'LIN','F','Linear axis (def=F)')
	  T1=1.0E-03_LDP*(1.0E+10_LDP/ASTRONOMICAL_UNIT())/DIST_KPC
	  XV(1:ND)=T1*R(1:ND)
	  IF(FLAG)THEN
	    XAXIS='r(")'
	  ELSE
	    XV(1:ND)=LOG10(XV(1:ND))
	    XAXIS='Log r(")'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XED')THEN
	  DO I=1,ND
	    XV(I)=LOG10(ED(I))
	  END DO
	  XAXIS='Log(N\de\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XTEMP')THEN
	  DO I=1,ND
	    XV(I)=T(I)
	  END DO
	  XAXIS='T(10\u4\dK)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XN')THEN
	  DO I=1,ND
	    XV(I)=I
	  END DO
	  XAXIS='I'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XATOM')THEN
	  DO I=1,ND
	    XV(I)=LOG10(POP_ATOM(I))
	  END DO
	  XAXIS='Log(N\di\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XVEL')THEN
	  DO I=1,ND
	    XV(I)=V(I)
	  END DO
	  XAXIS='V(km s\u-1\d)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XLOGV')THEN
	  DO I=1,ND
	    IF(V(I) .LE. 0.0_LDP)THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: XLOGV requires positive V'
	      WRITE(T_OUT,*)'Depth, V=',I,V(I)
	      STOP
	    END IF
	  END DO
	  CALL LOGVEC(V,XV,ND)
	  XAXIS='Log V(kms\u-1\d)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XCOLD')THEN
!
	  DO I=1,ND
	    ZETA(I)=1.0E+10_LDP*MASS_DENSITY(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    XV(I)=LOG10(TA(I))
	  END DO
	  XAXIS='Log m(gm cm\u-2\d)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XMASS')THEN
!
	  CALL USR_OPTION(FLAG,'MDIR','T','Integrate inwards from outer boundary')
	  DO I=1,ND
	    ZETA(I)=4.0_LDP*PI*1.0E+30_LDP*MASS_DENSITY(I)*R(I)*R(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  IF(FLAG)THEN
	    WRITE(6,*)'Mass in integrated from outer boudary.'
	    WRITE(6,*)'Normalizing mass (in Msun) is',TA(ND)/1.989D+33
	    DO I=1,ND
	      XV(I)=LOG10(TA(I)/TA(ND))
	    END DO
	    XAXIS='Log Mass Fraction'
	  ELSE
	    WRITE(6,*)'Mass in integrated from inner boudary.'
	    WRITE(6,*)'Total mass (in Msun) is',TA(ND)/1.989D+33
	    DO I=1,ND
	      XV(I)=(TA(ND)-TA(I))/1.989E+33_LDP
	    END DO
	    XAXIS='Mass (Msun)'
	  END IF
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XFLUX')THEN
	    DO I=1,ND
	      TA(I)=ABS(CLUMP_FAC(I)*FLUX_MEAN(I))
	    END DO
	    CALL TORSCL(TB,TA,R,XM,TC,ND,METHOD,TYPE_ATM)
	    WRITE(T_OUT,*)'Flux mean opacity set positive'
	    WRITE(T_OUT,*)'Flux optical depth is : ',TB(ND)
	    XV(1:ND)=LOG10(TB(1:ND))
	    XAXIS='Log(\gt\dFlux\u)'
	    XAXSAV=XAXIS
!
	ELSE  IF(XOPT .EQ. 'XTAUC')THEN
	  IF(.NOT. ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
	  CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    XV(I)=LOG10(TA(I))
	  END DO
	  XAXIS='Log(\gt\dc\u)'
	  XAXSAV=XAXIS
!
	ELSE IF(XOPT .EQ. 'XCAKT')THEN
	  CALL USR_HIDDEN(T1,'VTH','10.0','Thermal doppler velocity (km/s)')
	  DO I=1,ND
	    IF(V(I) .LE. 0.0_LDP .OR. SIGMA(I)+RONE .LE. 0.0_LDP)THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: XCAKT requires V > 0 and SIGMA > -1'
	      WRITE(T_OUT,*)'Depth, V, SIGMA=',I,V(I),SIGMA(I)
	      STOP
	    END IF
	    T2=(SIGMA(I)+RONE)*V(I)/R(I)
!	    XV(I)=LOG10(6.65D-15*T1/T2)
	    XV(I)=LOG10(6.65E-15_LDP*ED(I)*CLUMP_FAC(I)*T1/T2)
	  END DO
	  XAXIS='Log(t)'
!
	ELSE IF(XOPT .EQ. 'SET-ATM')THEN
	  CALL USR_OPTION(TYPE_ATM,'ATM',' ','Type of atmosphere: EXP, WIND or POW')
	  TYPE_ATM=UC(TYPE_ATM)
	  IF(TYPE_ATM .EQ. 'EXP')THEN
	  ELSE IF(TYPE_ATM .EQ. 'POW')THEN
	    TYPE_ATM='P'
	    T1=LOG(MASS_DENSITY(5)/MASS_DENSITY(1))/LOG(R(1)/R(5))
	    IF(T1 .LT. 100)THEN
	      WRITE(TYPE_ATM(2:6),'(F5.2)')T1
	      WRITE(6,'(A,F5.2)')'Density exponent for atmospheres is ',T1
	    ELSE
	      WRITE(TYPE_ATM(2:9),'(ES8.2)')T1
	      WRITE(6,'(A,ES8.2)')'Density exponent for atmospheres is ',T1
	    END IF
	  ELSE
	    TYPE_ATM=' '
	    WRITE(T_OUT,*)'Atmosphere is assumed to have a wind'
	  ENDIF
!
	ELSE IF(XOPT .EQ. 'SET-METH')THEN
	  CALL USR_OPTION(METHOD,'METH',' ','Method for derivative evaluation: LOGLOG, LOGMON, ZERO')
	  METHOD=UC(METHOD)
	  IF(METHOD .NE. 'ZERO' .AND.
	1       METHOD .NE. 'LOGLOG' .AND.
	1       METHOD .NE. 'LOGMON')THEN
	    METHOD='LOGLOG'
	    WRITE(T_OUT,*)'Did not recognize option'
	    WRITE(T_OUT,*)'Setting METHOD=LOGLOG'
	  ENDIF
!
	ELSE IF(XOPT .EQ. 'SET-IBC')THEN
	  DEFAULT=INNER_BND_METH
	  CALL USR_OPTION(INNER_BND_METH,'IBC',DEFAULT,
	1      'Inner boundary condition: DIFFUSION, SCHUSTER, ZERO_FLUX, HOLLOW')
	  INNER_BND_METH=UC(INNER_BND_METH)
	  IF(INNER_BND_METH .NE. 'DIFFUSION' .AND. INNER_BND_METH .NE. 'SCHUSTER' .AND.
	1       INNER_BND_METH .NE. 'ZERO_FLUX' .AND. INNER_BND_METH .NE. 'HOLLOW')THEN
	    INNER_BND_METH='DIFFUSION'
	    WRITE(T_OUT,*)'Did not recognize option'
	    WRITE(T_OUT,*)'Setting INNER_BND_METH to DIFFUSION'
	  ENDIF
	  DIF=.TRUE.
	  IF(INNER_BND_METH .NE. 'DIFFUSION')DIF=.FALSE.
! 
! **************************************************************************
! **************************************************************************
!
! Information options. Provide a means of obtaining information about the
! model without using the debug option.
!
! **************************************************************************
! **************************************************************************
!

        ELSE IF(XOPT .EQ. 'WRFREQ')THEN
	  OPEN(UNIT=10,FILE='EDGE_FREQ',STATUS='NEW')
	    DO ID=1,NUM_IONS
	     IF(ATM(ID)%XzV_PRES)CALL EDGEWR(ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%NXzV_F,UC(ION_ID(ID)),10)
	    END DO
	  CLOSE(UNIT=10)
	  WRITE(T_OUT,*)'Edge frequencie written to unit 10'
! 
        ELSE IF(XOPT .EQ. 'WRN')THEN
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)WRITE(6,'(1X,A,T43,A,I4)')
	1        'Number of levels in FULL'//TRIM(ION_ID(ID)),
	1        ' model atom is:',ATM(ID)%NXzV_F
	  END DO
!
        ELSE IF(XOPT .EQ. 'WRA')THEN
	  CALL USR_HIDDEN(ELEC,'WRF','T','Write final ionization stage for each species?')
	  CALL USR_HIDDEN(FOUND,'DIAG','F','Output diagnostic info for name convesiont to latex format')
	  CALL WR_ATOM_SUM(ELEC,FOUND)
!
        ELSE IF(XOPT .EQ. 'WRID')THEN
!
	  CALL USR_HIDDEN(LEV,2,2,'LIMS','1,0','Imin, Imax')
!
	  FLAG=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)))THEN
	      IF(LEV(2) .EQ. 0)LEV(2)=ATM(ID)%NXzV_F
	      K=0
	      DO I=LEV(1),LEV(2)
	        K=MAX(K,LEN_TRIM(ATM(ID)%XzVLEVNAME_F(I)))
	      END DO
	      IF(K .LE. 16)THEN
	        WRITE(T_OUT,'(5(2X,A16,1X,I3))')(TRIM(ATM(ID)%XzVLEVNAME_F(I)),I,I=LEV(1),LEV(2))
	      ELSE
	        WRITE(T_OUT,'(3(2X,A30,1X,I3))')(TRIM(ATM(ID)%XzVLEVNAME_F(I)),I,I=LEV(1),LEV(2))
	      END IF
	      FLAG=.TRUE.
	    END IF
	  END DO
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)'Error --- invalid species identification'
	  END IF
!
	ELSE IF(XOPT .EQ. 'WRLST')THEN
	  OPEN(UNIT=LU_OUT,FILE='LST_LEVL_INFO',STATUS='UNKNOWN')
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      K=0
	      DO J=2,ATM(ID)%NXzV_F-1
	        DO L=1,J-1
	          IF(ATM(ID)%AXzV_F(J,L) .GE. 1.0E-20_LDP)K=K+1
	        END DO
	      END DO
	      I=ATM(ID)%NXzV_F
	      WRITE(LU_OUT,'(A7,3X,I4,I8,3X,A30)')ION_ID(ID),ATM(ID)%NXzV_F,K,ATM(ID)%XzVLEVNAME_F(I)
	    END IF
	  END DO
	  CLOSE(LU_OUT)
!
	ELSE IF(XOPT .EQ. 'WRTAB')THEN
	  OPEN(UNIT=10,FILE='NEW_TABLE',STATUS='UNKNOWN')
	  DO I=1,ND
	    TA(I)=CLUMP_FAC(I)*ROSS_MEAN(I)*(R(ND)/R(I))**2
	  END DO
	  CALL TORSCL(TAUROSS,TA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    WRITE(10,'(ES12.6,ES12.4,4X,F8.4,2X,ES12.4)')1.0E+10*R(I),POP_ATOM(I),T(I),TAUROSS(I)
	  END DO
!
	ELSE IF(XOPT .EQ. 'CHKA')THEN
!
	  WRITE(6,*)' '
          WRITE(6,*)BLUE_PEN//'Writing out levels with no decay routes'//DEF_PEN
	  WRITE(6,*)' '
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=2,ATM(ID)%NXZV_F
	        T1=SUM(ATM(ID)%AXzV_F(J,1:J-1))
	        IF(T1 .EQ. 0.0_LDP)THEN
	          WRITE(6,'(A6,3X,I4,4X,A)')ION_ID(ID),J,TRIM(ATM(ID)%XzVLEVNAME_F(J))
	          WRITE(30,'(A6,3X,I4,4X,A)')ION_ID(ID),J,TRIM(ATM(ID)%XzVLEVNAME_F(J))
	        END IF
	      END DO
	    END IF
	  END DO
!
! 
!
! Change units for input of frequencies/waevengths.
!
	ELSE IF(XOPT .EQ. 'XKEV')THEN
	  KEV_INPUT=.TRUE.
	  HZ_INPUT=.FALSE.
	  ANG_INPUT=.FALSE.
	  WRITE(T_OUT,*)'Frequencies are now input in keV'
	  FREQ_INPUT='Units are keV'
	ELSE IF(XOPT .EQ. 'XHZ')THEN
	  KEV_INPUT=.FALSE.
	  HZ_INPUT=.TRUE.
	  ANG_INPUT=.FALSE.
	  WRITE(T_OUT,*)'Frequencies are now input in Hz'
	  FREQ_INPUT='Units are 10^15 Hz'
	ELSE IF(XOPT .EQ. 'XANG')THEN
	  KEV_INPUT=.FALSE.
	  HZ_INPUT=.FALSE.
	  ANG_INPUT=.TRUE.
	  WRITE(T_OUT,*)'Frequencies are now input in Angstroms.'
	  FREQ_INPUT='Units are Angstroms'
!
! 
!
	ELSE IF(XOPT .EQ. 'DIE')THEN
	  DO I=1,ND
	    ZNET(I)=RONE-RJ(I)*CHIL(I)/ETAL(I)
	  END DO
	  WRITE(LU_NET,40001)'Znet for DIELECTRONIC transition'
	  WRITE(LU_NET,40002)LEV(1)
	  WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
	  WRITE(T_OUT,40003)(ZNET(I),I=1,ND)
!
! 
	ELSE IF(XOPT .EQ. 'NETR')THEN
!
	  CALL USR_OPTION(SOL_OPT,'FORM','FG','Method: FORM, FG, or HAM?')
!
	  IF(SOL_OPT(1:2) .EQ. 'FG')THEN
	    N_TYPE='N_ON_J'
	    CALL USR_OPTION(N_TYPE,'NTYPE','N_ON_J','N_ON_J .OR. G_ONLY?')
	    IF(UC(N_TYPE(1:1)) .EQ. 'N')N_TYPE='N_ON_J'
	    IF(UC(N_TYPE(1:1)) .EQ. 'G')N_TYPE='G_ONLY'
!
	    CALL USR_OPTION(IDEPTH,10,1,'DEPTHS','0','Depths to be plotted [0 plots JBAR] ')
	    CALL USR_OPTION(RED_EXT,'RED_EXT','14.0D0','Redward extension in Doppler widths' )
	    CALL USR_OPTION(VDOP_FG_FRAC,'FG_FRAC','0','Doppler fraction for FG[V10] ')
	    CALL USR_OPTION(VDOP_MOM_FRAC,'MOM_FRAC','0','Doppler fraction for MOM[V6] ')
	    CALL USR_HIDDEN(ELEC,'THK','T','Thick outer boundary in continuum & line')
	    CALL USR_OPTION(FG_SOL_OPT,'SOL_METH','INT/INS','FG solution method: INT/INS or DIFF/INS')
!
	    PLT_J=.FALSE.; PLT_H=.FALSE.; PLT_LF=.FALSE.; PLT_LF=.FALSE.; PLT_FM=.FALSE.
	    IF(IDEPTH(1) .EQ. 0)THEN
	    ELSE
	      CALL USR_OPTION(PLT_J,'PLT_J','F','Plot J')
	      CALL USR_OPTION(PLT_H,'PLT_H','F','Plot H')
	      CALL USR_OPTION(PLT_LF,'PLT_LF','F','Plot chi(L).H')
	      IF(PLT_J .OR. PLT_H .OR. PLT_FM)THEN
	      ELSE
	        CALL USR_OPTION(PLT_FM,'PLT_FM','T','Plot force multiplier')
	      END IF
	    END IF
!
! Note that CLUMP_FAC has already been included in ETA, CHI, ESEC, ETAL, and CHIL.
!
	    CALL COMP_JBAR(ETA,CHI,CHI_SCAT,ETAL,CHIL,
	1                T,V,SIGMA,R,P,
	1                JQW,HMIDQW,KQW,NMIDQW,
	1                JBAR,DIF,DTDR,IC,ELEC,
	1                VTURB,VDOP_FG_FRAC,VDOP_MOM_FRAC,
	1                RED_EXT,AMASS,FL,
	1                METHOD,FG_SOL_OPT,N_TYPE,
	1                LUM,PLT_J,PLT_H,PLT_LF,PLT_FM,
	1                LEV,10,NLF,NC,NP,ND)
	        DO I=1,ND
	          ZNET(I)=RONE-JBAR(I)*CHIL(I)/ETAL(I)
	        END DO
	      YAXIS=' '
	      IF(PLT_J)YAXIS='R\u2\dJ'
	      IF(PLT_H .AND. YAXIS .NE. ' ')YAXIS=TRIM(YAXIS)//'; R\u2\ddH'
	      IF(PLT_H .AND. YAXIS .EQ. ' ')YAXIS='R\u2\dH'
	      IF(PLT_LF .AND. YAXIS .NE. ' ')YAXIS=TRIM(YAXIS)//'; R\u2\dLF'
	      IF(PLT_LF .AND. YAXIS .EQ. ' ')YAXIS='R\u2\d^LF'
	      IF(PLT_LF .AND. YAXIS .EQ. ' ')YAXIS='MT'
	   ELSE
	    IDEPTH(1)=0
	    CALL USR_OPTION(FULL_ES,'ALLES','T',
	1        'Include line photons scatterd in resonace zone ?')
	    CALL USR_HIDDEN(SKIPEW,'SKIPEW','F','Skip accurate EW computation ?')
	    CALL USR_HIDDEN(EWACC,'EWACC','0.05','required % accuracy in EW')
!
	    INACCURATE=.TRUE.
	    CNT=0
	    EWOLD=1.0E+37_LDP
	    CALL DP_ZERO(GAM,ND)
	    DO WHILE(INACCURATE)
!
! The .FALSE. option is use to indicate that we do not require the
! approximate ANR operator to be computed. Use TA for APPROX_LAM
! which is not computed anyway.
!
	      IF(SOL_OPT(1:3) .EQ. 'HAM')THEN
	        CALL HAM_FORMSOL(ETA_WITH_ES,CHI,CHI_SCAT,CHIL,ETAL,V,SIGMA,R,P,
	1                  JBAR,ZNET,TA,TB,TC,.FALSE.,
	1                  RJ,AV,CV,GAM,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  JQW,HMIDQW,
	1                  PF,LINE_PRO,LFQW,ERF,FL,DIF,DBB,IC,
	1                  RONE,THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
	      ELSE
	        CALL FORMSOL(ETA_WITH_ES,CHI,CHI_SCAT,CHIL,ETAL,V,SIGMA,R,P,
	1                  JBAR,ZNET,TA,TB,TC,.FALSE.,
	1                  RJ,AV,CV,GAM,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  JQW,HMIDQW,
	1                  PF,LINE_PRO,LFQW,ERF,FL,DIF,DBB,IC,
	1                  RONE,THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
	      END IF
	      WRITE(T_OUT,180)EW,CONT_INT,LAMVACAIR(FL)
180	      FORMAT(1X,'EW=',1PE14.6,3X,'Ic=',E14.6,3X,'Lam=',E14.6)
!
! Use a NG acceleration every four iterations to speed up convergence.
! We use FEDD for JPREV since it is not used with FORMSOL.
!
	      CNT=CNT+1
	      L=5-MOD(CNT,4)
	      IF(L .EQ. 5)L=1
	      CALL SETFORNG(GAM,FEDD,L,ND)
	      IF(L .EQ. 1)CALL NGACCEL(GAM,FEDD,ND,.TRUE.)
!
	      ERR(CNT)=200.0_LDP*(EW-EWOLD)/(EW+EWOLD)
	      EWOLD=EW
	      IF( ABS(ERR(CNT)) .LT. EWACC )INACCURATE=.FALSE.
	      IF(CNT .GT. 20)THEN
	        WRITE(T_OUT,'(1X,'' Max changes for looop '',I3)')
	        WRITE(T_OUT,'(1P,(2X,5E12.3))')(ERR(I),I=1,CNT)
	        INACCURATE=.FALSE.
	      END IF
	      IF(SKIPEW)INACCURATE=.FALSE.
!
	    END DO
	  END IF
!
	  DO I=1,ND
	    YV(I)=LOG10(JBAR(I))
	  END DO
	  IF(IDEPTH(1) .EQ. 0)THEN
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(Jbar)'
	  END IF
!
	  WRITE(LU_NET,40001)'FORMSOL Transfer Solution'
	  WRITE(LU_NET,40002)LEV(1),LEV(2)
	  WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
	  FLUSH(UNIT=LU_NET)
! 
	ELSE IF(XOPT .EQ. 'MOMR')THEN
!
	  CALL USR_HIDDEN(HAM,'HAM','F','Use second order differencing ?')
	  CALL USR_OPTION(FULL_ES,'ALLES','T',
	1      'Include line photons scatterd in resonace zone ?')
	  CALL USR_HIDDEN(SKIPEW,'SKIPEW','F','Skip accurate EW computation ?')
!
	  INACCURATE=.TRUE.
	  CNT=0
	  EWOLD=1.0E+37_LDP
	  CALL DP_ZERO(JNU,ND*(NLF+1))
	  DO WHILE(INACCURATE)
!
! As ETAL and CHIL are known, the FORMAL solution gives FEDD and GEDD
! directly. We only need to iterate to obtian a reliable EW value in
! the presnece of electron scattering. We use ETAEDD for consistency
! purposes - ETAEDD contains the continuum electron scattering term
! computed using Eddington Factors.
!
	    IF(HAM)THEN
	      CALL FG_HAM(ETA_WITH_ES,CHI,CHI_SCAT,RJ,
	1                  CHIL,ETAL,V,SIGMA,R,P,
	1                  JNU,HNU,FEDD,GEDD,
	1                  JQW,HMIDQW,KQW,NMIDQW,
	1                  INBC_VEC,HBC_VEC,NBC_VEC,
	1                  PF,LINE_PRO,LFQW,ERF,FL,
	1                  EW,CONT_INT,LINE_BL,
	1                  DIF,DBB,IC,METHOD,
	1                  THK_CONT,THK_LINE,NLF,NC,NP,ND)
	      WRITE(T_OUT,*)'FGHAM_EW=',EW,'Ic=',CONT_INT
!
! We use TA for RADEQ, and TB for the FLUX vectors returned by the
! EW computation.
!
	      CALL MOMHAM(ETA_WITH_ES,CHI,CHI_SCAT,THETA,RJ,CHIL,ETAL,
	1                  V,SIGMA,R,JBAR,ZNET,
	1                  JNU,HNU,FEDD,GEDD,
	1                  HBC_VEC,INBC_VEC,NBC_VEC,TA,TB,
	1                  PF,LINE_PRO,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  MOMEW,MOMCONT_INT,LINE_BL,FULL_ES,
	1                  NLF,NC,NP,ND)
	      WRITE(T_OUT,180)MOMEW,MOMCONT_INT,LAMVACAIR(FL)
	    ELSE
	      CALL FG_COMP(ETA_WITH_ES,CHI,CHI_SCAT,RJ,
	1                  CHIL,ETAL,V,SIGMA,R,P,
	1                  JNU,HNU,FEDD,GEDD,
	1                  JQW,HMIDQW,KQW,NMIDQW,
	1                  INBC_VEC,HBC_VEC,NBC_VEC,
	1                  PF,LINE_PRO,LFQW,ERF,FL,
	1                  EW,CONT_INT,LINE_BL,
	1                  DIF,DBB,IC,METHOD,
	1                  THK_CONT,THK_LINE,NLF,NC,NP,ND)
	      WRITE(T_OUT,*)'FGCOMP_EW=',EW,'Ic=',CONT_INT
!
! We use TA for RADEQ, and TB for the FLUX vectors returned by the
! EW computation.
!
	      CALL MOMJBAR(ETA_WITH_ES,CHI,CHI_SCAT,THETA,RJ,CHIL,ETAL,
	1                  V,SIGMA,R,JBAR,ZNET,
	1                  JNU,HNU,FEDD,GEDD,
	1                  HBC_VEC,INBC_VEC,NBC_VEC,TA,TB,
	1                  PF,LINE_PRO,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  MOMEW,MOMCONT_INT,LINE_BL,FULL_ES,
	1                  NLF,NC,NP,ND)
	      WRITE(T_OUT,180)MOMEW,MOMCONT_INT,LAMVACAIR(FL)
	    END IF
!
	    CNT=CNT+1
	    ERR(CNT)=200.0_LDP*(MOMEW-EWOLD)/(MOMEW+EWOLD)
	    EWOLD=MOMEW
	    IF( ABS(ERR(CNT)) .LT. 0.001_LDP )INACCURATE=.FALSE.
	    IF(CNT .GT. 20)THEN
	      WRITE(T_OUT,'(1X,'' Max changes for looop '',I3)')
	      WRITE(T_OUT,'(1P,(2X,5E12.3))')(ERR(I),I=1,CNT)
	      INACCURATE=.FALSE.
	    END IF
	    IF(SKIPEW)INACCURATE=.FALSE.
!
	  END DO
!
	  WRITE(LU_NET,40001)'MOMHAM Transfer Solution'
	  WRITE(LU_NET,40002)LEV(1),LEV(2)
	  WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
!
! 
!
	ELSE IF( XOPT .EQ. 'SOBR' )THEN
!
	  CALL USR_HIDDEN(NORAD,'NORAD','F','Ignore radiation filed ?')
	  CALL USR_HIDDEN(DOWN_RATE,'DWNR','F','Compute new downward rate: Nu Aul Zul')
!
	  DO I=1,ND
	    SOURCE(I)=ETA_WITH_ES(I)/CHI(I)
	  END DO
!
! We use U and GB for BETA and BETAC respectively.
!
	  CALL SOBJBAR(SOURCE,CHI,CHIL,ETAL,
	1                V,SIGMA,R,P,JQW,
	1                JBAR,ZNET,VB,VC,U,GB,
	1                FL,DIF,DBB,IC,THK_CONT,NLF,NC,NP,ND,METHOD)
!
	  DO I=1,ND
	    YV(I)=LOG10(JBAR(I))
	  END DO
	  YAXIS='Log(Jbar)'
	  CALL DP_CURVE(ND,XV,YV)
!
	  IF(NORAD)THEN
	    WRITE(LU_NET,40001)
	1           'Sobolev Approx but radiation field not included'
	    WRITE(LU_NET,40002)LEV(1),LEV(2)
	    WRITE(LU_NET,40003)(U(I),I=1,ND)
	  ELSE IF(DOWN_RATE)THEN
	    WRITE(LU_NET,40001)'Ne Aul Zul'
	    WRITE(LU_NET,40002)LEV(1),LEV(2)
	    WRITE(LU_NET,40003)(ETAL(I)*ZNET(I)/EMLIN/FL,I=1,ND)
	  ELSE
	    WRITE(LU_NET,40001)'Sobolev Approximation'
	    WRITE(LU_NET,40002)LEV(1),LEV(2)
	    WRITE(LU_NET,40003)(ZNET(I),I=1,ND)
	  END IF
!
! 
!
! This section recomputes the LTE populations. Useful for adjusting
! the Temperature and hence the population levels to assist in obtaining
! a converged model.
!
	ELSE IF(XOPT .EQ. 'HMIN')THEN
	  DO I=1,ND
	    YV(I)=2.07E-22_LDP*ED(I)*ATM(1)%XzV_F(1,I)*EXP(HDKT*0.754_LDP/4.135667_LDP/T(I))
	    YV(I)=LOG(YV(I)/2.0_LDP/(T(I)**1.5_LDP))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'LTE')THEN
	  INCLUDE 'EVAL_LTE_FULL.INC'
!
! 
!
	ELSE IF(XOPT .EQ. 'EXTRAP')THEN
!
	  CALL USR_OPTION(LEV,ITWO,ITWO,'EWIN','1,1',
	1        'Depth band to be extrapolated over')
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO I=LEV(1),LEV(2)
	        T1=POP_ATOM(I)/POP_ATOM(LEV(2)+1)
	        DO J=1,ATM(ID)%NXzV_F
	          ATM(ID)%XzV_F(J,I)=T1*ATM(ID)%XzV_F(J,LEV(2)+1)
	        END DO
	        ATM(ID)%DXzV_F(I)=T1*ATM(ID)%DXzV_F(LEV(2)+1)
	      END DO
	    END IF
	  END DO
	  T(LEV(1):LEV(2))=T(LEV(2)+1)
!
	  INCLUDE 'EVAL_LTE_FULL.INC'
!
	ELSE IF(XOPT .EQ. 'WRDC' .OR. XOPT .EQ. 'WRPOP' .OR.
	1       XOPT .EQ. 'WRLTE' .OR. XOPT .EQ. 'WRTX' .OR. XOPT .EQ. 'WRREC')THEN
	  IF(XOPT .EQ. 'WRDC')THEN
	    TYPE='DC'
	  ELSE IF(XOPT .EQ. 'WRLTE')THEN
	    TYPE='LTE'
	  ELSE IF(XOPT .EQ. 'WRTX')THEN
	    TYPE='TX'
	  ELSE IF(XOPT .EQ. 'WRREC')THEN
	    TYPE='REC'
	  ELSE
	    TYPE='POP'
	  END IF
	  CALL USR_OPTION(OVERWRITE,'OVER','F','Overwrite existing files?')
	  CALL USR_HIDDEN(STRING,'EXT',TYPE,'File appendage')
	  IF(STRING .EQ. 'POP')STRING='POPS'
	  WRITE(T_OUT,*)'File names will have the form XzV'//TRIM(STRING)
!
! Set depth vector. The default allows will write all depths. We
! can omit ceartain depths using the OWIN option.
!
	  DO_DPTH(1:ND)=.TRUE.
	  CALL USR_HIDDEN(LEV,ITEN,ITWO,'OWIN','0,0,0,0,0,0,0,0,0,0',
	1        'Depth bands (pairs)NOT to be output [Max of 5]')
	  DO I=1,9,2
	    IF(LEV(I) .GT. 0 .AND. LEV(I+1) .GE. LEV(I) .AND.
	1        LEV(I+1) .LE. ND)THEN
	      DO_DPTH(LEV(I):LEV(I+1))=.FALSE.
	    END IF
	  END DO
	  T2=RZERO
	  CALL USR_HIDDEN(T2,'ROFFSET','0.0D0','Value to subtract from R')
!
! Write out departure coefficients to ASCI file.
! NB - 1 refers to dimension of DHYD (i.e. DHYD(1,nd)
!      1 refers to format for output.
!      1,NHY - For use with HEI.
!
	  TA(1:ND)=R(1:ND)-T2
	  IF(XOPT .EQ. 'WRLTE')THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        FILENAME=TRIM(ION_ID(ID))//TRIM(STRING)
	        CALL NEW_WRITEDC_V6(ATM(ID)%XzVLTE_F,ATM(ID)%LOG_XzVLTE_F,ATM(ID)%W_XzV_F,
	1             ATM(ID)%AXzV_F,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,ATM(ID)%NXzV_F,
	1             ATM(ID)%DXzV_F,ATM(ID)%GIONXzV_F,IONE,TA,T,ED,V,CLUMP_FAC,
	1             DO_DPTH,LUM,ND,FILENAME,TYPE,IONE,OVERWRITE,IOS)
	      END IF
	      IF(IOS .NE. 0)EXIT
	    END DO
	  ELSE
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        FILENAME=TRIM(ION_ID(ID))//TRIM(STRING)
	        CALL NEW_WRITEDC_V6(ATM(ID)%XzV_F,ATM(ID)%LOG_XzVLTE_F,ATM(ID)%W_XzV_F,
	1             ATM(ID)%AXzV_F,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,ATM(ID)%NXzV_F,
	1             ATM(ID)%DXzV_F,ATM(ID)%GIONXzV_F,IONE,TA,T,ED,V,CLUMP_FAC,
	1             DO_DPTH,LUM,ND,FILENAME,TYPE,IONE,OVERWRITE,IOS)
	      END IF
	      IF(IOS .NE. 0)EXIT
	    END DO
	  END IF
!
!
! 
!
! This page computes the Rosseland mean opacity from the temperature
! distribution and the population levels. An option allows the electron
! scattering opacity to be excluded from the calculations. The Rossland
! optical depth scale is given in TAUROSS, and TA is a working vector. The
! Rossland opacity is given in CHIROSS. Not written as a subroutine to allow
! easy inclusion of additional opacity sources.
!
	ELSE IF(XOPT .EQ. 'SXROSS' .OR. XOPT .EQ. 'SYROSS')THEN
	  DO I=1,ND
	    TA(I)=CLUMP_FAC(I)*ROSS_MEAN(I)*(R(ND)/R(I))**2
	  END DO
	  CALL TORSCL(TAUROSS,TA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(T_OUT,*)'Spherical rossland optical depth is : ',TAUROSS(ND)
	  IF(X .EQ. 'SXROSS')THEN
	    DO I=1,ND
	      XV(I)=LOG10(TAUROSS(I))
	    END DO
	    XAXIS='Log(\gt\dSRoss\u)'
	    XAXSAV=XAXIS
	  ELSE IF(XOPT .EQ. 'SYROSS')THEN
	    DO I=1,ND
	      YV(I)=LOG10(TAUROSS(I))
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\gt\dSRoss\u)'
	  END IF
!
	ELSE IF(XOPT .EQ. 'XROSS' .OR. XOPT .EQ. 'YROSS' .OR.
	1                                 XOPT .EQ. 'GREY')THEN
	  IF(ROSS_MEAN(1) .EQ. 0.0_LDP)THEN
	    WRITE(T_OUT,*)'Rosseland mean opacity unavailabe'
	    GOTO 1
	  END IF
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  ROSS=.TRUE.
	  WRITE(T_OUT,*)'Rossland optical depth is : ',TAUROSS(ND)
!
! 
	  IF(XOPT .EQ. 'XROSS')THEN
	    DO I=1,ND
	      XV(I)=LOG10(TAUROSS(I))
	    END DO
	    XAXIS='Log(\gt\dRoss\u)'
	    XAXSAV=XAXIS
	  ELSE IF(XOPT .EQ. 'YROSS')THEN
	    DO I=1,ND
	      YV(I)=LOG10(TAUROSS(I))
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\gt\dRoss\u)'
	  ELSE
!
! We have two different methods avaialable to compute the Grey Temperature
! distribution. One takes into account the velocity terms, and is needed for
! supernovae calculations.
!
	    IF(V(1) .GT. 5.0E+03_LDP)THEN
	      CALL USR_OPTION(GREY_REL,'GREY_REL','T','Use for relativistic solution')
	    ELSE
	      CALL USR_OPTION(GREY_REL,'GREY_REL','F','Use for relativistic solution')
	    END IF
	    CALL USR_HIDDEN(GREY_WITH_V_TERMS,'VT','F','File appendage')
	    CALL USR_HIDDEN(ELEC,'LOGT','F','Log of T?')
	    IF(R(1)/R(ND) .LT. RTWO .AND. V(1) .LT. RONE)THEN
	      CALL USR_OPTION(PLANE_PARALLEL_NOV,'PPNOV','T','Plane parallel model no V?')
	    ELSE
	      CALL USR_HIDDEN(PLANE_PARALLEL_NOV,'PPNOV','F','Plane parallel model no V?')
	    END IF
!
! Compute Grey temperature distribution.
!
! Will use FA for F, GAM for NEWRJ, and GAMH for NEWRK
!
	    YAXIS='T(10\u4 \dK)'
	    IF(GREY_REL)THEN
!
	      TA(1:ND)=0.3_LDP              !FEDD: Initial guess
	      TB(1:ND)=RZERO              !H_ON_J: Initial guess
	      TC(1:ND)=RZERO              !dlnJdlnR=0.0D0
	      GAMH(1:ND)=RZERO            !Old FEDD
	      XM(1:ND)=RZERO              !As grey solution, not needed (ETA)
	      NEW_FREQ=.TRUE.
	      HBC=0.7_LDP; INBC=0.1_LDP
	      WRITE(T_OUT,*)'Using MOM_JREL_GREY_V1 for grey solution'
!
! Note
!   HFLUX=LUM*Lsun/16/(PI*PI)/10**2/R**2 (10**2 for 1/R**2).
!   DBB = dBdR = 3.Chi.L/16(piR)**2
!   DBB is used for the lower boundary diffusion approximation.
!
	      HFLUX=3.826E+13_LDP*LUM/(4.0_LDP*PI*R(ND))**2
	      DBB=3.0_LDP*HFLUX*CHIROSS(ND)
	      IF(.NOT. DIF)DBB=RZERO
	      T1=RONE
	      DO WHILE(T1 .GT. 1.0E-05_LDP)
	        CALL MOM_JREL_GREY_V1(XM,CHIROSS,CHIROSS,V,SIGMA,R,
	1              TB,TA,TC,RJ,HNU,HBC,INBC,
	1              DIF,DBB,IC,METHOD,
	1              L_TRUE,L_TRUE,NEW_FREQ,ND)
!
	        CALL FGREY_NOREL_V1(TA,TB,RJ,CHIROSS,R,V,SIGMA,
	1              P,JQW,HMIDQW,KQW,LUM,IC,METHOD,
	1              HBC,INBC,DIF,ND,NC,NP)
	        T1=RZERO
	        DO I=1,ND
	          T1=MAX(ABS(TA(I)-GAMH(I)),T1)
	          GAMH(I)=TA(I)
	        END DO
	        NEW_FREQ=.FALSE.
	        WRITE(T_OUT,*)'Current grey iteration accuracy is',T1
	      END DO
	      CALL TORSCL(TA,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	      DO I=1,ND
	        ZV(I)=LOG10(TA(I))
	        TGREY(I)=((3.14159265_LDP/5.67E-05_LDP*RJ(I))**0.25_LDP)*1.0E-04_LDP
	        YV(I)=TGREY(I)
	      END DO
	      GREY_COMP=.TRUE.
	      IF(ELEC)THEN
	        YV(1:ND)=LOG10(YV(1:ND))
	        YAXIS='Log T(10\u4 \dK)'
	      END IF
	      CALL DP_CURVE(ND,ZV,YV)
	      XAXIS='Log(\gt\dRoss\u)'
	      XAXSAV=XAXIS
!

	  ELSE IF(GREY_WITH_V_TERMS)THEN
	      T2=1.0E-05_LDP          !Accuracy to converge f
              CALL JGREY_WITH_FVT(RJ,TB,CHIROSS,R,V,SIGMA,
	1                  P,JQW,HMIDQW,KQW,NMIDQW,
	1                  LUM,METHOD,DIF,IC,
	1                  T2,ND,NC,NP)
	      CALL TORSCL(TA,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	      DO I=1,ND
	        ZV(I)=LOG10(TA(I))
	        TGREY(I)=((3.14159265_LDP/5.67E-05_LDP*RJ(I))**0.25_LDP)*1.0E-04_LDP
	        YV(I)=TGREY(I)
	      END DO
	      GREY_COMP=.TRUE.
!
	      IF(ELEC)THEN
	        YV(1:ND)=LOG10(YV(1:ND))
	        YAXIS='Log T(10\u4 \dK)'
	      END IF
	      CALL DP_CURVE(ND,ZV,YV)
	      XAXIS='Log(\gt\dRoss\u)'
	      XAXSAV=XAXIS
	    ELSE IF(PLANE_PARALLEL_NOV)THEN
	      FOLD=0.3333_LDP
	      HBC=0.7_LDP; NBC=RZERO; FL=RONE; INBC=0.1_LDP
!
! Note HFLUX=LUM*Lsun/16/(PI*PI)/10**2 (10**2 for 1/R**2).
! DBB - dBdR = 3. Chi. L/16(piR)**2 and is used for the lower boundary
! diffusion approximation. Since we are dealing with a plane-parallel
! atmopshere, we divide HFLUX by R*^2.
!
              HFLUX=3.826E+13_LDP*LUM/16.0_LDP/PI**2/R(ND)/R(ND)
              DBB=3.0_LDP*CHIROSS(ND)*HFLUX
              IF(.NOT. DIF)DBB=RZERO
	      T1=1000.0_LDP
!
! Compute radial (vertical) optical depth increments.
!
	      CALL DERIVCHI(TB,CHIROSS,R,ND,METHOD)
	      CALL NORDTAU(DTAU,CHIROSS,R,R,TB,ND)
!
! Compute the solution vector. Note that the units need to be
! eventually included. The following follows direcly from d2K/d2Tau=0.
!
	      DO WHILE(T1 .GT. 1.0E-08_LDP)
	        T2=FOLD(1)/HBC
	        RJ(1)=HFLUX/HBC
	        DO I=2,ND
	          T2=T2+DTAU(I-1)
	          RJ(I)=T2*HFLUX/FOLD(I)
	        END DO
!
! Will use TA for IPLUS.
!
	        SOURCE(1:ND)=RJ
	        CALL FCOMP_PP_V2(R,TC,GAMH,SOURCE,CHIROSS,TA,HBC,
	1               NBC,INBC,DBB,IC,THK_CONT,DIF,ND,NC,METHOD)
	        T1=RZERO
                DO I=1,ND
                  T1=MAX(ABS(FOLD(I)-GAMH(I)),T1)
                  FOLD(I)=GAMH(I)
                END DO
              END DO
!
! Compute the temperature distribution, and the Rossland optical depth scale.
! Assumes LTE. NB sigma=5.67E-05 and the factor of 1.0E-04 is to convert
! T from units of K to units of 10^4 K.
!
	      GREY_COMP=.TRUE.
	      DO I=1,ND
	        TGREY(I)=((3.14159265_LDP/5.67E-05_LDP*RJ(I))**0.25_LDP)*1.0E-04_LDP
	      END DO
	      WRITE(T_OUT,*)'Use TGREY to plot T'
!
	    ELSE
!
! We use a finer grid here. The finer grid has been  previously defined.
!
! Interpolate CHIROSS onto a finer grid.
!
	      DO I=1,NDX
	        CHI(I)=RZERO
	        DO J=0,3
	          CHI(I)=CHI(I)+COEF(J,I)*LOG( CHIROSS(J+INDX(I)) )
	        END DO
	        CHI(I)=EXP(CHI(I))
	      END DO
!
	      DO I=1,NDX
	        FA(I)=RONE/RTHREE
	      END DO
	      HBC=RONE
	      T1=1000.0
	      DO WHILE(T1 .GT. 1.0E-05_LDP)
	        CALL JGREY(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,RJ,
	1          GAM,GAMH,Q,FA,CHI,dCHIdR,
	1          JQWEXT,KQWEXT,LUM,HBC,HBCNEW,NCX,NDX,NPX,METHOD)
	        T1=RZERO
	        DO I=1,NDX
	          T1=MAX(ABS(FA(I)-GAMH(I)),T1)
	          FA(I)=GAMH(I)
	        END DO
	        T1=MAX(ABS(HBC-HBCNEW),T1)
	        HBC=HBCNEW
	        WRITE(T_OUT,'('' Maximum change is'',1PE11.4)')T1
	      END DO
!
! Compute the temperature distribution, and the Rossland optical depth scale.
! Assumes LTE. NB sigma=5.67E-05 and the factor of 1.0E-04 is to convert
! T from units of K to units of 10^4 K. We use ZV for XV so as not to corrupt
! XV.
!
	      CALL TORSCL(TA,CHI,REXT,TB,TC,NDX,METHOD,TYPE_ATM)
	      DO I=1,NDX
	        ZV(I)=LOG10(TA(I))
	        YV(I)=((3.14159265_LDP/5.67E-05_LDP*RJ(I))**0.25_LDP)*1.0E-04_LDP
	      END DO
!
! Store grey temperature on Model Grid.
!
	      GREY_COMP=.TRUE.
	      DO I=1,ND
	        TGREY(I)=YV(GRID(I))
	      END DO
!
	      CALL USR_HIDDEN(ELEC,'LOGT','F','Log of T?')
	      IF(ELEC)THEN
	        DO I=1,NDX
	          YV(I)=LOG10(YV(I))
	        END DO
	        YAXIS='Log T(10\u4 \dK)'
	      END IF
	      CALL DP_CURVE(NDX,ZV,YV)
	      XAXIS='Log(\gt\dRoss\u)'
	      XAXSAV=XAXIS
	    END IF
!
! Set REXT(1)=0 to insure that REXT is recomputed for the interpolation
! section.
!
!	    REXT(1)=0.0
	  END IF
!
! 
!
! ****************************************************************************
! ****************************************************************************
!
! Y axis options.
!
! ****************************************************************************
! ****************************************************************************
!
	ELSE IF(XOPT .EQ. 'YFLUX')THEN
	  CALL USR_OPTION(ELEC,'dFLUX','F','Plot dTAU (flux mean opacity)')
	  DO I=1,ND
	    TA(I)=ABS(CLUMP_FAC(I)*FLUX_MEAN(I))
	  END DO
	  CALL TORSCL(TB,TA,R,XM,TC,ND,METHOD,TYPE_ATM)
	  WRITE(T_OUT,*)'Flux optical depth is : ',TB(ND)
!
	  IF(ELEC)THEN
	    DO I=1,ND-1
	      IF(TB(I) .GT. 0)THEN
	        YV(I)=LOG10(TB(I+1)/TB(I))
	      ELSE
	        YV(I)=10.0_LDP
	      END IF
	    END DO
	    CALL DP_CURVE(ND-1,XV,YV)
	    YAXIS='dLog \gt(Flux)'
	  ELSE
	    DO I=1,ND
	      IF(TB(I) .GT. 0)THEN
	        YV(I)=LOG10(TB(I))
	      ELSE
	        YV(I)=RZERO
	      END IF
	    END DO
	    YAXIS='Log \gt(Flux)'
	    CALL DP_CURVE(ND,XV,YV)
	  END IF
!
	ELSE IF(XOPT .EQ. 'ROSS' .OR. XOPT .EQ. 'FLUX' .OR. XOPT .EQ.  'ES' .OR.
	1          XOPT .EQ. 'PLANCK' .OR. XOPT .EQ. 'ABS')THEN
	  IF(XOPT .EQ. 'ROSS')THEN
	     TA(1:ND)=1.0E-10_LDP*ROSS_MEAN(1:ND)
	     CURVE_LAB='Rosseland'
	  ELSE IF(XOPT .EQ. 'FLUX')THEN
	     TA(1:ND)=1.0E-10_LDP*FLUX_MEAN(1:ND)
	     CURVE_LAB='Flux'
	  ELSE IF(XOPT .EQ. 'PLANCK')THEN
	     TA(1:ND)=1.0E-10_LDP*PLANCK_MEAN(1:ND)
	     CURVE_LAB='Planck'
	  ELSE IF(XOPT .EQ. 'ABS')THEN
	     TA(1:ND)=1.0E-10_LDP*ABS_MEAN(1:ND)
	     CURVE_LAB='Absorption'
	  ELSE
	     TA(1:ND)=6.65E-25_LDP*ED(1:ND)
	     CURVE_LAB='E.S.'
	  END IF
	  CALL USR_OPTION(ELEC,'KAPPA','T','Mass absorption coefficient?')
	  IF(TA(1) .NE. 0.0_LDP)THEN
	    IF(ELEC)THEN
	      DO I=1,ND
	        YV(I)=TA(I)/MASS_DENSITY(I)
	      END DO
	      YAXIS='\gk (cm\u2\d/g)'
	    ELSE
	      CALL USR_OPTION(ELEC,'ON_NE','T','Normalize by the electron scattering opacity?')
	      IF(ELEC)THEN
	        DO I=1,ND
	          YV(I)=TA(I)/(6.65E-25_LDP*ED(I))
	        END DO
	        YAXIS='\gx/\gsNe'           ! (cm\u-1\d)'
	      ELSE
	        DO I=1,ND
	          YV(I)=LOG10(TA(I))
	        END DO
	        YAXIS='\gx (cm\u-1\d)'
	        WRITE(T_OUT,*)'Volume filling factor not allowed for.'
	      END IF
	    END IF
	    CALL DP_CURVE_LAB(ND,XV,YV,CURVE_LAB)
	  ELSE
	    WRITE(T_OUT,*)'Opacity not available.'
	  END IF
!
	ELSE IF(XOPT .EQ. 'JINT')THEN
	  CURVE_LAB='Log (Integrated J moment)'
	  CALL DP_CURVE_LAB(ND,XV,J_INT,CURVE_LAB)
	ELSE IF(XOPT .EQ. 'HINT')THEN
	  CURVE_LAB='Integrated H moment'
	  CALL DP_CURVE_LAB(ND,XV,H_INT,CURVE_LAB)
	ELSE IF(XOPT .EQ. 'KINT')THEN
	  CURVE_LAB='Integrated K moment'
	  CALL DP_CURVE_LAB(ND,XV,K_INT,CURVE_LAB)
	ELSE IF(XOPT .EQ. 'BINT')THEN
	  TA(1:ND)=1.0E+16_LDP*STEFAN_BOLTZ()*(T(1:ND)**4)/PI
	  CURVE_LAB='Integrated Planck Function (sigma T^4 / pi)'
	  CALL DP_CURVE_LAB(ND,XV,TA,CURVE_LAB)
!
	ELSE IF(XOPT .EQ. 'YLOGR')THEN
	  DO I=1,ND
	    YV(I)=LOG10(R(I)/R(ND))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(r/R\d*\u)'
!
	ELSE IF(XOPT .EQ. 'ED')THEN
	  DO I=1,ND
	    YV(I)=LOG10(ED(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(N\de\u)'
!
	ELSE IF(XOPT .EQ. 'YDEN')THEN
	  DO I=1,ND
	    YV(I)=LOG10(MASS_DENSITY(I))
	  END DO
	  YAXIS='\gr(gm cm\u-3\d)'
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'YCOLD')THEN
	  IF(XSPEC .EQ. ' ')THEN
            CALL USR_OPTION(ELEC,'NA','F','Atom column density?')
	    IF(ELEC)THEN
	      DO I=1,ND
	        ZETA(I)=1.0E+10_LDP*POP_ATOM(I)*CLUMP_FAC(I)
	      END DO
	      YAXIS='N(atoms cm\u-2\d)'
	      DEFAULT='Atom column density'
	    ELSE
	      DO I=1,ND
	        ZETA(I)=1.0E+10_LDP*MASS_DENSITY(I)*CLUMP_FAC(I)
	      END DO
	      YAXIS='m(gm cm\u-2\d)'
	      DEFAULT=' '
	    END IF
	  ELSE
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)))THEN
	        ZETA(1:ND)=SUM(ATM(ID)%XzV_F,1)
	        ZETA(1:ND)=1.0E+10_LDP*ZETA(1:ND)*CLUMP_FAC(1:ND)
	        YAXIS='N(atoms cm\u-2\d)'
	        DEFAULT='N('//TRIM(ION_ID(ID))//')'
	        EXIT
	      END IF
	    END DO
	  END IF
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    YV(I)=LOG10(TA(I))
	  END DO
	  CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
!
! Used to verify that model has a constant mass-loss rate. It can also be used to
! see Mdot variation in time dependent model.
!
	ELSE IF(XOPT .EQ. 'YMDOT')THEN
	  T1=4.0E+25_LDP*PI*365.25_LDP*24.0_LDP*3600.0_LDP/MASS_SUN()
!	  WRITE(6,*)T1
!	  WRITE(6,*)MASS_DENSITY(1),CLUMP_FAC(1),R(1),V(1)
	  DO I=1,ND
	    WRITE(26,*)T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)*V(I)
	    YV(I)=LOG10( T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)*V(I) )
	  END DO
	  WRITE(6,'(A,ES14.4,A)')'Mass loss rate is ',10**(YV(1)),' Msun/yr'
	  WRITE(6,'(A,ES14.4,A)')'     Log(Mdot) is ',YV(1)
	  YAXIS='Mass Loss rate(Msun/yr)'
	  CALL DP_CURVE(ND,XV,YV)
!
! Designed to look at Mdot in shell models. Since these models are not smooth, we avearge
! (sum) data over the smoothing step size.
!
	ELSE IF(XOPT .EQ. 'YAVMDOT')THEN
	  CALL USR_OPTION(T4,'dLOGR','0.1','Smoothing step in log R space')
	  T1=4.0E+30_LDP*PI/MASS_SUN()
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(6,*)'Done mass determination'
	  DO I=1,ND
	    ZETA(I)=1.0E+05_LDP/V(I)
	  END DO
	  CALL TORSCL(ZV,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(6,*)'Done time determination'
	  DO I=1,ND-1
	    T1=365.25_LDP*24.0_LDP*3600.0_LDP
	    J=MIN(I+1,ND)
	    DO K=I+1,ND
	      IF( LOG10(R(I)/R(J)) .GT. T4)EXIT
	      J=K
	    END DO
	    YV(I)=LOG10( T1*(TA(J)-TA(I))/(ZV(J)-ZV(I)) )
	  END DO
	  T1=4.0E+25_LDP*PI*365.25_LDP*24.0_LDP*3600.0_LDP/MASS_SUN()
	  I=ND
	  YV(I)=LOG10( T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)*V(I) )
!
	  YAXIS='Mass Loss rate(Msun/yr)'
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'YF')THEN
	  YAXIS='Clumping factor'
	  CALL DP_CURVE(ND,XV,CLUMP_FAC)
!
! Designed to look at f, the clumping factor, in shell models. Since these models are not
! smooth, we avearge (sum) data over the smoothing step size.
!
	ELSE IF(XOPT .EQ. 'YAVF')THEN
	  CALL USR_OPTION(T4,'dLOGR','0.1','Smoothing step in log R space')
!
	  T1=4.0E+30_LDP*PI
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*R(I)*R(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
!
	  DO I=1,ND
	    ZETA(I)=ZETA(I)*MASS_DENSITY(I)
	  END DO
	  CALL TORSCL(ZV,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
!
! Compute dV
!
	  T1=4.0E+30_LDP*PI
	  DO I=1,ND
	    ZETA(I)=T1*R(I)*R(I)
	  END DO
	  CALL TORSCL(XM,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
!
	  CALL DETERM_CLUM_POS(MASS_DENSITY,R,ND,TC,K)
	  J=1
	  IF(K .GT. 2)THEN
	    L=(TC(2)-TC(1))/2
	  ELSE
	    L=1
	    WRITE(6,*)'L=',L
	  END IF
	  DO I=1,ND-5
	    IMIN=MAX(1,I-L)
	    IMAX=MIN(ND,I+L)
	    YV(I)=(TA(IMAX)-TA(IMIN))*(TA(IMAX)-TA(IMIN))/(ZV(IMAX)-ZV(IMIN))/(XM(IMAX)-XM(IMIN))
	    IF(I .GT. TC(K)+L/2)THEN
	      DO J=I,ND-1
	        IMAX=J+1; IMIN=J
	        YV(J)=(TA(IMAX)-TA(IMIN))*(TA(IMAX)-TA(IMIN))/(ZV(IMAX)-ZV(IMIN))/(XM(IMAX)-XM(IMIN))
	      END DO
	      YV(ND)=RONE
	      EXIT
	    END IF
	  END DO
!
!	  DO I=1,ND-1
!	    J=MIN(I+1,ND)
!	    DO K=I+1,ND
!	      IF(LOG10(R(I)/R(J)) .GT. T4)EXIT
!	      J=K
!	    END DO
!	    YV(I)=(TA(J)-TA(I))*(TA(J)-TA(I))/(ZV(J)-ZV(I))/(XM(J)-XM(I))
!	  END DO
!	  YV(ND)=RONE
!
	  YAXIS='Clumping factor'
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'YMASS')THEN
	  CALL USR_OPTION(TMP_LOGICAL,'OUTW','F','Intgerate mass outwards?')
	  T1=4.0E+30_LDP*PI/MASS_SUN()
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)
	  END DO
!
! Note that TB[i] contains dM[i] (R[i] to R[i+1]).
!
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  IF(TMP_LOGICAL)THEN
	    TA(ND)=RZERO
	    DO I=ND-1,1,-1
	      TA(I)=TA(I+1)+TB(I)
	    END DO
	  END IF
	  YAXIS='M(\dsun\u)'
	  WRITE(6,'(A,ES10.3,A)')'Mass of envelope (ejecta) is',TA(ND),' Msun'
	  CALL DP_CURVE(ND,XV,TA)
!
	ELSE IF(XOPT .EQ. 'IMASS')THEN
	  K=0
	  CALL USR_OPTION(TMP_LOGICAL,'OUTW','F','Intgerate mass outwards?')
	  WRITE(6,'(A)')' '
	  DO ISPEC=1,NSPEC
	    IF( (XSPEC .EQ. SPECIES(ISPEC) .OR. XSPEC .EQ. 'ALL') .AND.  POPDUM(ND,ISPEC) .GT. 0.0_LDP)THEN
	      T1=4.0E+30_LDP*PI*AT_MASS(ISPEC)*ATOMIC_MASS_UNIT()/MASS_SUN()
	      DO I=1,ND
	        ZETA(I)=T1*POPDUM(I,ISPEC)*CLUMP_FAC(I)*R(I)*R(I)
	      END DO
	      CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	      IF(TMP_LOGICAL)THEN
	        TA(ND)=RZERO
	        DO I=ND-1,1,-1
	          TA(I)=TA(I+1)+TB(I)
	        END DO
	      END IF
              DEFAULT=TRIM(SPECIES_ABR(ISPEC))
              J=INDEX(DEFAULT,'k'); IF(J .NE. 0)DEFAULT(J:J)='i'
              CALL DP_CURVE_LAB(ND,XV,TA,DEFAULT)
	      K=K+1
	      WRITE(6,'( A,A4,A,ES9.2,A)',ADVANCE='NO')' Mass of ',TRIM(SPECIES(ISPEC)),' is',TA(ND),' Msun'
	      IF(MOD(K,2) .NE. 0)WRITE(6,'(10X)',ADVANCE='NO')
	      IF(MOD(K,2) .EQ. 0)WRITE(6,'(A)')' '
	    END IF
	    IF(MOD(K,2) .NE. 0. .AND. XSPEC .EQ. 'ALL')WRITE(6,'(A)')' '
	  END DO
	  YAXIS='Mass(M\d\(9)\u)'
!
	ELSE IF(XOPT .EQ. 'LOGT')THEN
	  CALL LOGVEC(T,YV,ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(T[10\u4\dK])'
!
	ELSE IF(XOPT .EQ. 'TGREY')THEN
	  IF(GREY_COMP)THEN
	    CALL DP_CURVE(ND,XV,TGREY)
	    YAXIS='T(10\u4\dK)'
	  ELSE
	    WRITE(T_OUT,*)'Error -Grey Temperature distribution not available'
	    WRITE(T_OUT,*)'Call GREY option first --- ',
	1                'has higher spatial resolution.'
	  END IF
!
	ELSE IF(XOPT .EQ. 'ROP')THEN
	  WRITE(6,'(2X,A,2X,(A,3X))')'Depth','Log(T(K))',' LogR','Log K'
	  DO I=1,ND
	    YV(I)=LOG10(1.0E-10_LDP*ROSS_MEAN(I)/MASS_DENSITY(I))
	    ZV(I)=LOG10(1.0E+06_LDP*MASS_DENSITY(I)/T(I)**3)
	    WRITE(6,'(2X,I5,3X,3F8.3)')I,LOG10(T(I))+4.0D0,ZV(I),YV(I)
	  END DO
	  CALL DP_CURVE(ND,ZV,YV)
	  XAXIS='\gr/T\u3\d\d6\u'
	  YAXIS='\gk(ross)'
!
	ELSE IF(XOPT .EQ. 'THOP')THEN
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)*R(ND)*R(ND)/R(I)/R(I)
	  END DO
	  CALL TORSCL(TA,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  I=1
	  DO WHILE(TAUROSS(I) .LT. 0.1_LDP)
	    I=I+1
	  END DO
	  J=1
	  DO WHILE(TAUROSS(J) .LT. 1.0_LDP)
	    J=J+1
	  END DO
	  K=1
	  DO WHILE(TAUROSS(K) .LT. 3.0_LDP)
	    K=K+1
	  END DO
!	  Q0=1.0D0; QINF=0.5D0; GAM_HOPF=1.0
          T2=1.0E+10_LDP*R(ND)/RAD_SUN()
          T2=1.0E-04_LDP*TEFF_SUN()*(LUM/T2**2)**0.25_LDP                !Units of 10^4 K
	  DO I=1,ND
	    T1=1.3333333_LDP*( (T(I)/T2)**4 )*TAUROSS(I)/TA(I)-TAUROSS(I)
	    WRITE(6,*)I,TAUROSS(I),TA(I),T(I),T2
	    YV(I)=T1
	  END DO
	  CALL DP_CURVE(ND,TA,YV)
!
	ELSE IF(XOPT .EQ. 'TRAT')THEN
!
! Plots the ration of the ACTUAL T distribution to the GREY value.
!
	  IF(GREY_COMP)THEN
	    DO I=1,ND
	      YV(I)=T(I)/TGREY(I)
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='T(10\u4\dK)'
	  ELSE
	    WRITE(T_OUT,*)'Error -Grey Temperature distribution not available'
	    WRITE(T_OUT,*)'Call GREY option first --- '
	  END IF
!
	ELSE IF(XOPT .EQ. 'TCMF')THEN
	  DO I=1,ND
	    YV(I)=CMFGEN_TGREY(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='T(10\u4\dK)'
!
	ELSE IF(XOPT .EQ. 'T')THEN
	  DO I=1,ND
	    YV(I)=T(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='T(10\u4\dK)'
!
! Compute dlnT/dlnP for comparison with adiabatic temperature gradient.
!
	ELSE IF(XOPT .EQ. 'DTDP')THEN
	  T1=1.0E+04_LDP*BOLTZMANN_CONSTANT()
	  TA(1:ND)=T1*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)
	  ELEC=.TRUE.
	  CALL USR_OPTION(ELEC,'RP','T','Include radiation pressure (assumes black body)')
	  IF(ELEC)THEN
	    T1=1.0E+16_LDP*STEFAN_BOLTZ()/SPEED_OF_LIGHT()
	    TA(1:ND)=TA(1:ND)+T1*(T(1:ND)**4)
	  END IF
	  TA(1:ND)=LOG( TA(1:ND) )
	  TB(1:ND)=LOG( T(1:ND) )
	  CALL DERIVCHI(TC,TB,TA,ND,'LINMON')
	  YV(1:ND)=TC(1:ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='dlnT/dlnP'
!
	ELSE IF(XOPT .EQ. 'DERAD')THEN
	  CALL USR_OPTION(ELEC,'INTEG','F','Integrated luminosity')
	  IF(ELEC)THEN
	    YV(1:ND)=3.280E-03_LDP*dE_RAD_DECAY(1:ND)*CLUMP_FAC(1:ND)*R(1:ND)*R(1:ND)  !(4*PI*Dex(+30)/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='E(rad)(L\dsun\u)'
	  ELSE
	    YV(1:ND)=LOG10(dE_RAD_DECAY(1:ND))
	    YAXIS='Log \ge(ergs\u \dcm\u-3 \ds\u-1\d)'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'EK')THEN
	  CALL USR_OPTION(ELEC,'INTEG','T','Integrate thermal kinetic energy?')
	  YV(1:ND)=1.5E+04_LDP*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)*BOLTZMANN_CONSTANT()
	  IF(ELEC)THEN
	    YV(1:ND)=3.280E-03_LDP*YV(1:ND)*R(1:ND)*R(1:ND)  !(4*PI*Dex(+30)/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='E(kinetic)(s.L\dsun\u)'
	  ELSE
	    YV(1:ND)=LOG10(YV(1:ND))
	    YAXIS='Log EK(ergs\u \dcm\u-3)'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'ERAD')THEN
	  CALL USR_OPTION(ELEC,'INTEG','T','Integrate radiative energy?')
	  T1=4.0E+16_LDP*STEFAN_BOLTZ()/SPEED_OF_LIGHT()
	  YV(1:ND)=T1*(T(1:ND)**4)
	  IF(ELEC)THEN
	    YV(1:ND)=4.0E+30_LDP*PI*YV(1:ND)*R(1:ND)*R(1:ND)  !(*PI*Dex(+30) !/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='Log E(rad)'
	    WRITE(6,'(A,ES14.4,A)')'Total radiation density ',YV(1),' s.Lsun'
	    YV(1:ND-1)=LOG10(YV(1:ND-1))
	    CALL DP_CURVE(ND-1,XV,YV)
	  ELSE
	    YAXIS='Log E\drad\u(ergs\u \dcm\u-3)'
	    YV(1:ND)=LOG10(YV(1:ND))
	    CALL DP_CURVE(ND,XV,YV)
	  END IF
!
	ELSE IF(XOPT .EQ. 'EI')THEN
	  CALL USR_OPTION(ELEC,'INTEG','T','Integrate intenal energy?')
          YV(1:ND)=RZERO
!
	  DO ISPEC=1,NSPEC
	    T1=RZERO			!Total energy of state
	    T2=RZERO                    !Running ioization energy
!
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV_F
	        T1=T2+(ATM(ID)%EDGEXzV_F(1)-ATM(ID)%EDGEXzV_F(I))
	        DO J=1,ND
	          YV(J)=YV(J)+T1*ATM(ID)%XzV_F(I,J)
	        END DO
	      END DO
	      T2=T2+ATM(ID)%EDGEXzV_F(1)
	    END DO
	  END DO
	  YV(1:ND)=1.0E+04_LDP*HDKT*YV(1:ND)*BOLTZMANN_CONSTANT()
	  IF(ELEC)THEN
	    YV(1:ND)=3.280E-03_LDP*YV(1:ND)*R(1:ND)*R(1:ND)  !(4*PI*Dex(+30)/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    YAXIS='E(kinetic)(s.L\dsun\u)'
	  ELSE
	    YV(1:ND)=LOG10(YV(1:ND))
	    YAXIS='Log EI(ergs\u \dcm\u-3)'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'PGAS')THEN
	  TA(1:ND)=LOG10(1.0E+04_LDP*BOLTZMANN_CONSTANT()*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND))
	  CALL DP_CURVE(ND,XV,TA)
	  YAXIS='Pgas'
!
	ELSE IF(XOPT .EQ. 'PGONP')THEN
!
	  WRITE(6,*)'Assumes a black body spectrum to compute P(rad).'
	  TA(1:ND)=1.0E+04_LDP*BOLTZMANN_CONSTANT()*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)
	  T1=4.0E+16_LDP*STEFAN_BOLTZ()/SPEED_OF_LIGHT()/3.0_LDP
	  TB(1:ND)=T1*T(1:ND)**4
	  YV(1:ND)=TA(1:ND)/TB(1:ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Pgas/P'
!
	ELSE IF(XOPT .EQ. 'DPDR')THEN
	  TA(1:ND)=1.0E+04_LDP*BOLTZMANN_CONSTANT()*(ED(1:ND)+POP_ATOM(1:ND))*T(1:ND)
	  DO I=2,ND-1
	    YV(I)=1.0E-10_LDP*(TA(I+1)-TA(I-1))/(R(I-1)-R(I+1))/MASS_DENSITY(I)
	    TC(I)=1.0E-10_LDP*TA(I)*LOG(TA(I+1)/TA(I-1))/LOG(R(I-1)/R(I+1))/R(I)/MASS_DENSITY(I)
!	    TC(I)=LOG(V(I+1)/V(I-1))/LOG(R(I+1)/R(I-1))
!	    TB(I)=1.0D-03*SQRT(BOLTZMANN_CONSTANT()*(POP_ATOM(I)+ED(I))*T(I)/MASS_DENSITY(I))
!	    TC(I)=TB(I)*TB(I)*TC(I)/R(I)
	  END DO
	  YV(1)=YV(2); YV(ND)=YV(ND-1)
	  TC(1)=TC(2); TC(ND)=TC(ND-1)
	  DO I=1,ND
	    TB(I)=1.0E-03_LDP*SQRT(BOLTZMANN_CONSTANT()*(POP_ATOM(I)+ED(I))*T(I)/MASS_DENSITY(I))
	    TB(I)=TB(I)*TB(I)*(SIGMA(I)+RONE)/R(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DP_CURVE(ND,XV,TB)
	  CALL DP_CURVE(ND,XV,TC)
	  YAXIS='\gr\u-1\d|dPdR|'
!
	ELSE IF(XOPT .EQ. 'BEK')THEN
	  YV(1:ND)=0.5E+10_LDP*MASS_DENSITY(1:ND)*V(1:ND)*V(1:ND)
	  YV(1:ND)=4*PI*1.0E+30_LDP*YV(1:ND)*R(1:ND)*R(1:ND)
	  CALL LUM_FROM_ETA(YV,R,ND)
	  DO I=ND-1,1,-1
	    YV(I)=YV(I+1)+YV(I)
	  END DO	
	  EK_EJECTA=YV(1)
	  WRITE(6,*)RED_PEN
	  WRITE(6,'(A,ES10.3,A)')' Kinetic energy of ejecta is ',EK_EJECTA,' ergs'
	  YAXIS='Log EK(ergs\u \dcm\u-3)'
	  YV(1:ND-1)=LOG10(YV(1:ND-1))
	  CALL DP_CURVE(ND-1,XV,YV)
!
	  T1=4.0E+30_LDP*PI/MASS_SUN()
	  DO I=1,ND
	    ZETA(I)=T1*MASS_DENSITY(I)*CLUMP_FAC(I)*R(I)*R(I)
	  END DO
	  CALL TORSCL(TA,ZETA,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(6,'(A,F9.4,A)')' Mass of envelope (ejecta) is',TA(ND),' Msun'
	  T1=SQRT(2.0_LDP*EK_EJECTA/TA(ND)/MASS_SUN())/1.0E+05_LDP
	  WRITE(6,'(A,F10.2,A)')' SQRT(Mean square velocity) is ',T1,' km/s',DEF_PEN
!
!
!
	ELSE IF(XOPT .EQ. 'RONV')THEN
	  T1=1.0E+10_LDP/1.0E+05_LDP/24.0_LDP/3600.0_LDP
	  DO I=1,ND
	    YV(I)=T1*R(I)/V(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='t(days)'
!
	ELSE IF(XOPT .EQ. 'VEL')THEN
	  DO I=1,ND
	    YV(I)=V(I)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='V(kms\u-1\d)'
!
	ELSE IF(XOPT .EQ. 'LOGV') THEN
	  DO I=1,ND
	    IF(V(I) .LE. 0.0_LDP)THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: LOGV requires positive V'
	      WRITE(T_OUT,*)'Depth, V=',I,V(I)
	      STOP
	    END IF
	  END DO
	  CALL LOGVEC(V,YV,ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log V(kms\u-1\d) '
!
	ELSE IF(XOPT .EQ. 'VOL') THEN
	  CALL USR_OPTION(T1,'SN_AGE','1.0','SN age in days')
	  CALL USR_OPTION(T2,'dAGE','0.1','Fractional change in SN age')
	  T2=1.0E-05_LDP*T1*24.0_LDP*3600.0_LDP*T2
	  DO I=1,ND
	    TA(I)=(RONE+V(I)*T2/R(I))**2
	    TA(I)=TA(I)*(RONE+V(I)*T2/R(I)*(SIGMA(I)+RONE))
	    IF(TA(I) .LE. 0.0_LDP)THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: VOL requires a positive expansion factor'
	      WRITE(T_OUT,*)'Depth, V, SIGMA=',I,V(I),SIGMA(I)
	      STOP
	    END IF
	    TA(I)=TA(I)**(RONE/RTHREE)
	  END DO
	  CALL DP_CURVE(ND,XV,TA)
	  YAXIS='Vf/Vi**0.33'
!
	ELSE IF(XOPT .EQ. 'RZERO') THEN
	  CALL USR_OPTION(T1,'SN_AGE','1.0','SN age in days')
	  T2=1.0E-05_LDP*T1*24.0_LDP*3600.0_LDP
	  DO I=1,ND
	    TA(I)=R(I)-V(I)*T2
	  END DO
	  CALL DP_CURVE(ND,XV,TA)
	  YAXIS='Rzero'
!
! Allows various velocity parameters to be varied to test their
! effect on the velocity law.
!
	ELSE IF(XOPT .EQ. 'TSTV')THEN
	  CALL USR_OPTION(SCLHT,'SCLHT','0.01D0','DEL R/R*')
	  CALL USR_OPTION(VCORE,'VCORE','1.0D0',' ')
	  CALL USR_OPTION(VPHOT,'VPHOT','100.0D0',' ')
	  CALL USR_OPTION(V_BETA1,'BETA1','1.0D0',' ')
	  DEFAULT='1.0D0'
	  IF(V_BETA1 .LT. 1.0_LDP)DEFAULT='0.999D0'
	  CALL USR_OPTION(V_EPS1,'EPS1',DEFAULT,' ')
	  CALL USR_OPTION(VINF1,'VINF1','1.0D0',' ')
	  CALL USR_OPTION(V_BETA2,'BETA2','1.0D0',' ')
	  CALL USR_OPTION(V_EPS2,'EPS2','1.0D0',' ')
	  DEFAULT=WR_STRING(VINF1)
	  CALL USR_OPTION(VINF2,'VINF2',DEFAULT,' ')
!
	  CALL USR_OPTION(TYPE,'PLOT','VEL','VEL, SIG, LOGV, LOGS')
!
	  CALL STARPCYG_TST(Z,DTAU,ZV,R(1),R(ND),
	1                   SCLHT,VCORE,VPHOT,VINF1,V_BETA1,V_EPS1,
	1                   VINF2,V_BETA2,V_EPS2,ND,TA,TB,TC,L_FALSE,LU_IN)
!
	  IF(XAXIS .EQ. 'Log(r/R\d*\u)')THEN
	    Z(1:ND)=LOG10(Z(1:ND)/R(ND))
	  ELSE IF(XAXIS .EQ. 'Log(r/R\d*\u)-1')THEN
	    XAXIS='Log(r/R\d*\u)-1'
	    Z(1:ND-1)=LOG10(Z(1:ND-1)/R(ND)-1)
	    Z(ND)=Z(ND-1)-RONE
	  ELSE
	    Z(1:ND)=Z(1:ND)/R(ND)
	  END IF
	  TYPE=UC( TRIM(TYPE) )
	  IF(TYPE .EQ. 'LOGV')THEN
	    CALL LOGVEC(DTAU,YV,ND)
	    YAXIS='Log(V(kms\u-1\d))'
	  ELSE IF(TYPE .EQ. 'SIG')THEN
	    YV(1:ND)=ZV(1:ND)+RONE
	    YAXIS='\gs+1'
	  ELSE IF(TYPE .EQ. 'LOGS')THEN
	    YV(1:ND)=LOG10(ZV(1:ND)+RONE)
	    YAXIS='LOG(\gs+1)'
	  ELSE
	    YV(1:ND)=DTAU(1:ND)
	    YAXIS='V(kms\u-1\d)'
	  END IF
	  CALL DP_CURVE(ND,Z,YV)
	  XAXSAV=XAXIS
	  XAXIS='Log(r/R\d*\u)'
!
	ELSE IF(XOPT .EQ. 'TSTF')THEN
	  DO I=1,ND
	    ESEC(I)=6.65E-15_LDP*ED(I)*CLUMP_FAC(I)
	  END DO
	  CALL TORSCL(TA,ESEC,R,TB,TC,ND,METHOD,TYPE_ATM)
	  CALL TST_CLUMP_LAW(YV,R,V,TA,ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Volume filling factor'
!
! To convert from cm/s to km/s we cale the sound speed by a factor of 10^{-5}. However,
! T in units of 10^4 K, thus the scale factor is on 10^{-3}.
!
	ELSE IF(XOPT .EQ. 'ISOC')THEN
	  DO I=1,ND
	    YV(I)=1.0E-03_LDP*SQRT(BOLTZMANN_CONSTANT()*(POP_ATOM(I)+ED(I))*T(I)/MASS_DENSITY(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='c(km/s)'
!
	ELSE IF(XOPT .EQ. 'SIGMA')THEN
!
! LOG10(SIGMA+1), LOG(V), and the finite-difference estimate of
! dLOG(V)/dLOG(R) require positive V, SIGMA > -1, and the same local
! velocity ordering throughout the atmosphere.
!
	  DO I=1,ND
	    IF(V(I) .LE. 0.0_LDP .OR. SIGMA(I)+1.0_LDP .LE. 0.0_LDP)THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: SIGMA requires V > 0 and SIGMA > -1'
	      WRITE(T_OUT,*)'Depth, V, SIGMA=',I,V(I),SIGMA(I)
	      STOP
	    END IF
	  END DO
	  DO I=2,ND-1
	    IF(V(I+1) .GE. V(I-1))THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: SIGMA requires strictly decreasing V'
	      WRITE(T_OUT,*)'Repeated or non-monotonic velocity at depth',I
	      WRITE(T_OUT,*)'V(I-1), V(I+1)=',V(I-1),V(I+1)
	      STOP
	    END IF
	  END DO
!
	  WRITE(6,'(/,1X,A)')'Plotting 3 curves for comparison:'
	  WRITE(6,'(6X,A)')'log10(sigma+1)'
	  WRITE(6,'(6X,A)')'log10(ln(V(I+1)/V(I-1))/ln(R(I+1)/R(I-1)))'
	  WRITE(6,'(6X,A,/)')'log10(dlni/dlnr)'
!
	  DO I=1,ND
	    YV(I)=LOG10(SIGMA(I)+1.0_LDP)
	  END DO
	  DO I=2,ND-1
	    TC(I)=LOG10(LOG(V(I+1)/V(I-1))/LOG(R(I+1)/R(I-1)))
	  END DO
	  TC(1)=TC(2); TC(ND)=TC(ND-1)
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DP_CURVE(ND,XV,TC)
	  DO I=1,ND
	    TA(I)=LOG(R(I))
	    TB(I)=LOG(V(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DERIVCHI(V,TB,TA,ND,'LINMON')
	  YAXIS='Log(\gs+1)'
!
	ELSE IF(XOPT .EQ. 'FONR')THEN
	  IF(ROSS_MEAN(1) .NE. 0.0_LDP)THEN
	    DO I=1,ND
	      YV(I)=FLUX_MEAN(I)/ROSS_MEAN(I)
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='\gx(Flux)/\gx(Ross)'
	  ELSE
	    WRITE(6,*)'Error --- Rosseland mean opacity not defined'
	  END IF
!
	ELSE IF(XOPT .EQ. 'CAKT')THEN
	  CALL USR_HIDDEN(T1,'VTH','10.0','Thermal doppler velocity (km/s)')
	  DO I=1,ND
	    IF(V(I) .LE. 0.0_LDP .OR. SIGMA(I)+RONE .LE. 0.0_LDP)THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: CAKT requires V > 0 and SIGMA > -1'
	      WRITE(T_OUT,*)'Depth, V, SIGMA=',I,V(I),SIGMA(I)
	      STOP
	    END IF
	    T2=(SIGMA(I)+RONE)*V(I)/R(I)
	    YV(I)=LOG10(6.65E-15_LDP*ED(I)*CLUMP_FAC(I)*T1/T2)
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(t)'
!
! We subtract 1 since M(t), by definition, does not include the
! force on the elctrons. With the present definition, it does include
! bound-free and free contributins (only important at depth.)
!
	ELSE IF(XOPT .EQ. 'MT')THEN
	  T1=1.0E-30_LDP*LUM_SUN()*LUM/4.0_LDP/PI/C_CMS
	  T2=T1*6.65E-15_LDP
	  DO I=1,ND
	    YV(I)=FLUX_MEAN(I)/ED(I)/6.65E-15_LDP - RONE
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='M(t)'
!
! Sobolev radial optical depth scale. Assumes fg=1, POP_ATOM for the
! level population.
!
	ELSE IF(XOPT .EQ. 'TAUSOB')THEN
	  CALL USR_OPTION(FL,'LAM',TAUSOB_DEF,'Wavelength in Ang')
	  TAUSOB_DEF=WR_STRING(FL)
	  FL=ANG_TO_HZ/FL
	  TMP_LOGICAL=.FALSE.
	  IF(XSPEC .EQ. ' ')THEN
	    TA(1:ND)=POP_ATOM(1:ND)
	    TMP_LOGICAL=.TRUE.
	  ELSE
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. SPECIES(ISPEC))THEN
	        TA(1:ND)=POPDUM(1:ND,ISPEC)
	        TMP_LOGICAL=.TRUE.
	        EXIT
	      END IF
	    END DO
	    IF(.NOT. TMP_LOGICAL)THEN
	      DO ID=1,NUM_IONS
	        IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)))THEN
	          TA(1:ND)=ATM(ID)%XzV_F(1,1:ND)
	          K=INDEX(ATM(ID)%XzVLEVNAME_F(1),'[')
	          IF(K .EQ. 0)THEN
	            J=1
	          ELSE
	            DO L=1,ATM(ID)%NXzV_F
	              IF(INDEX(ATM(ID)%XzVLEVNAME_F(L),ATM(ID)%XzVLEVNAME_F(1)(1:K)) .EQ. 0)EXIT
	              J=L
	            END DO
	          END IF
	          WRITE(6,'(X,A,10F6.1)')'Stat. weights for lower term are:',(ATM(ID)%GXzV_F(L),L=1,J)
	          TMP_LOGICAL=.TRUE.
	          EXIT
	        END IF
	      END DO
	    END IF
	  END IF
	  IF(TMP_LOGICAL)THEN
	    DO I=1,ND
	      T1=3.0E-10_LDP*OPLIN*TA(I)*R(I)/V(I)/FL
	      YV(I)=LOG10(T1)
	      ZV(I)=LOG10(T1/(RONE+SIGMA(I)))
	    END DO
	    CALL DP_CURVE(ND,XV,ZV)
	    YAXIS='Log(\gt\dSob\u/fX)'
	  ELSE
	    WRITE(6,*)'Unrecognized species or ion: XSPEC=',TRIM(XSPEC)
	  END IF
!
	ELSE IF(XOPT .EQ. 'YATOM')THEN
	  DO I=1,ND
	    YV(I)=LOG10(POP_ATOM(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(N\di\u)'
	ELSE IF(XOPT .EQ. 'YSPEC')THEN
	  FOUND=.FALSE.
	  ELEC=.FALSE.
	  FLAG=.FALSE.
	  TITLE=' '
	  CNT=0
!
	  WRITE(6,*)'Option plots fractional abundance (N/N[atom]), mass fraction, or species density.'
!
	  CALL USR_OPTION(ELEC,'FRAC','T','Fractional abundance')
	  IF(.NOT. ELEC)CALL USR_OPTION(FLAG,'MF','T','Mass fraction')
	  IF(ELEC)THEN
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. SPECIES(ISPEC) .OR. (XSPEC .EQ. 'ALL' .AND.
	1              POPDUM(ND,ISPEC) .GT. 0.0_LDP))THEN
	        YV(1:ND)=LOG10(POPDUM(1:ND,ISPEC)/POP_ATOM(1:ND)+1.0E-100_LDP)
	        DEFAULT=TRIM(SPECIES_ABR(ISPEC))
	        J=INDEX(DEFAULT,'k'); IF(J .NE. 0)DEFAULT(J:J)='i'
	        CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
	        FOUND=.TRUE.
	        IF(XSPEC .NE. 'ALL')EXIT
	      END IF
	    END DO
	    YAXIS='Log Fractional abundance (N\dX\u/N\dA\u)'
	  ELSE IF(FLAG)THEN
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. SPECIES(ISPEC) .OR. (XSPEC .EQ. 'ALL' .AND.
	1              POPDUM(ND,ISPEC) .GT. 0.0_LDP))THEN
	        T1=AT_MASS(ISPEC)*ATOMIC_MASS_UNIT()
	        DO I=1,ND
	          YV(I)=T1*POPDUM(I,ISPEC)/MASS_DENSITY(I)
	        END DO
	        FOUND=.TRUE.
	        CALL DP_CURVE(ND,XV,YV)
	        CNT=CNT+1; WRITE(DEFAULT,'(I2)')MOD(CNT,14)+1; DEFAULT=ADJUSTL(DEFAULT)
	        DEFAULT=TRIM(SPECIES_ABR(ISPEC))//', \p'//TRIM(DEFAULT)//' '
	        J=INDEX(DEFAULT,'k'); IF(J .NE. 0)DEFAULT(J:J)='i'
	        TITLE(CNT/10+1)=TRIM(TITLE(CNT/10+1))//' '//TRIM(DEFAULT)//' '
	        IF(XSPEC .NE. 'ALL')EXIT
	      END IF
	    END DO
	    YAXIS='Mass fraction'
	  ELSE
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. SPECIES(ISPEC) .OR. (XSPEC .EQ. 'ALL' .AND.
	1              POPDUM(ND,ISPEC) .GT. 0.0_LDP))THEN
	        YV(1:ND)=LOG10(POPDUM(1:ND,ISPEC)+1.0E-100_LDP)
	        DEFAULT=TRIM(SPECIES_ABR(ISPEC))
	        J=INDEX(DEFAULT,'k'); IF(J .NE. 0)DEFAULT(J:J)='i'
	        CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
	        IF(XSPEC .NE. 'ALL')EXIT
	      END IF
	    END DO
	    YAXIS='Log Species density (cm\u-3\d)'
	  END IF
	  WRITE(6,'(A)')' '
!
	  IF(XSPEC .EQ. 'ALL')THEN
	    CALL WR_SPEC_SUM_V2(ELEC,FLAG,XV,ND)
	  ELSE IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)'Error --- unrecognized species'
	  END IF
!
	ELSE IF(XOPT .EQ. 'YION')THEN
	  DO I=1,ND
	    YV(I)=LOG10(POPION(I))
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(N\di\u)'
	  WRITE(T_OUT,*)'Warning: this option plots ionized species only'
!
	ELSE IF(XOPT .EQ. 'CLUMP')THEN
	
	  CALL USR_HIDDEN(ELEC,'RECIP','F','Plot 1 / fillingfactor')
	  IF(ELEC)THEN
	    YV(1:ND)=RONE/CLUMP_FAC(1:ND)
	    YAXIS='Recipricoal Filling Factor'
	  ELSE
	    YV(1:ND)=CLUMP_FAC(1:ND)
	    YAXIS='Filling Factor'
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'MOD-DEN')THEN
	  CALL DP_CURVE(ND,XV,MASS_DENSITY)
	  CALL MODIFY_DENSITY(MASS_DENSITY,R,V,TYPE_ATM,ND)
	  CALL DP_CURVE(ND,XV,MASS_DENSITY)
	  WRITE(6,*)'Mass density has been corrupted'
!
	ELSE IF(XOPT .EQ. 'SCL-DEN')THEN
	  CALL USR_OPTION(T1,'SF','1.0D0','Factor to mass scale density')
	  MASS_DENSITY=MASS_DENSITY*T1
	  WRITE(6,*)'Mass density has been corrupted'
!	
	ELSE IF(XOPT .EQ. 'PROF')THEN
	  CALL PLT_PROFS()
!
!
!
	ELSE IF(XOPT .EQ. 'LNID')THEN
!
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')' Create line ID file for stars with weak winds.'
	  WRITE(T_OUT,'(A)')' The list is created based on the line optical depth at a give Tau(e.s.)'
	  WRITE(T_OUT,'(A)')' Use XTAUC and TAUC option to find Tau(e.s) at Tau(cont)=2/3'
	  WRITE(T_OUT,'(A)')' Use PLNID for plane-parallel models'
	  WRITE(T_OUT,'(A)')' PLNID may also work better for O stars with weak winds'
	  WRITE(T_OUT,'(A)')' '
	  DO I=1,ND
	    ESEC(I)=6.65E-15_LDP*ED(I)
	  END DO
          CALL TORSCL(TA,ESEC,R,TB,TC,ND,METHOD,TYPE_ATM)
	  CALL USR_OPTION(TAU_VAL,'TAUES','0.67D0',
	1          'Electron scattering optical depth at whch line optical depth is evaluated')
	  CALL USR_OPTION(VTURB,'VTURB','10.0','Turbulent velcity for line profile (units km/s)?')
!
	  DO I=1,ND
	    J=I
	    IF(TA(I) .GT. TAU_VAL)THEN
	      EXIT
	    END IF
	  END DO
	  DPTH_INDX=J
!
	  DEFAULT='0.01'
	  CALL USR_OPTION(TAU_LIM,'TAU',DEFAULT,'We only output lines with Tau > ?')
	  CALL USR_OPTION(FLAG,'ADD','F','Add transition name to output file?')
!
	  WRITE(T_OUT,'(/,A,/)')' Values at depth just exceeeding requested TAU(e.s)'
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'     R(I)/R*=',R(DPTH_INDX)/R(ND)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'        V(I)=',V(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'        T(I)=',T(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       ED(I)=',ED(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'TAU(e.s.)(I)=',TA(DPTH_INDX)
!
	  WRITE(T_OUT,'(A)')' '
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  DEFAULT='LINE_ID'
	  CALL USR_OPTION(FILENAME,'FILE',DEFAULT,'Output file for line IDs')
	  CALL GEN_ASCI_OPEN(73,FILENAME,'UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(73,'(A)')'Species         Lambda      Tau            Lam        B_loc'//
	1                '   E_loc       log gf        hv/kT    NL   NUP                  Name'
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
!
	  CNT=0
!
! Determine index range encompassing desired wavelength range.
!
	  NU_ST=ANG_TO_HZ/LAM_ST
	  IF(NU_ST .LT. VEC_FREQ(1))THEN
	    NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NL) .GT. NU_ST)NL=NL+1
	  ELSE
	    NL=1
	  END IF
!
	  NU_EN=ANG_TO_HZ/LAM_EN
	  IF(NU_EN .GT. VEC_FREQ(N_LINE_FREQ))THEN
	    NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NUP) .LT. NU_EN)NUP=NUP-1
	  ELSE
	    NUP=N_LINE_FREQ
	  END IF
!
! We use T4 for the interpolation to the requested optical depth.
!
	  T4=(TAU_VAL-TA(DPTH_INDX-1))/(TA(DPTH_INDX)-TA(DPTH_INDX-1))
	  DO LINE_INDX=NL,NUP
	    MNL_F=VEC_MNL_F(LINE_INDX)
	    MNUP_F=VEC_MNUP_F(LINE_INDX)
	    FL=VEC_FREQ(LINE_INDX)
!
	    IF(ANG_TO_HZ/FL .GT. LAM_ST .AND. ANG_TO_HZ/FL .LT. LAM_EN)THEN
	      DO ID=1,NUM_IONS
	        TAU_SOB=RZERO
	        IF(VEC_SPEC(LINE_INDX) .EQ. ION_ID(ID) .AND.
	1           (XSPEC .EQ. ' ' .OR. XSPEC .EQ. UC(ION_ID(ID))) )THEN
	          DO I=1,DPTH_INDX
	            T1=ATM(ID)%W_XzV_F(MNUP_F,I)/ATM(ID)%W_XzV_F(MNL_F,I)
	            CHIL(I)=OPLIN*VEC_OSCIL(LINE_INDX)*( T1*ATM(ID)%XzV_F(MNL_F,I)-
	1               ATM(ID)%GXzV_F(MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)/ATM(ID)%GXzV_F(MNUP_F) )
	            T2=12.86_LDP*SQRT( T(I)/AT_MASS(SPECIES_LNK(ID))+(VTURB/12.86_LDP)**2 )/C_KMS
	            CHIL(I)=MAX(1.0E-15_LDP*CHIL(I)/FL/T2/SQRT(PI),1.0E-10_LDP)
	          END DO
	          T1=(ANG_TO_HZ/FL)/5000.0_LDP                               !Normalized to 5000 Ang.
                  T2=LOG10(VEC_OSCIL(LINE_INDX)*ATM(ID)%GXzV_F(MNL_F)*T1)
                  T3=HDKT*(ATM(ID)%EDGEXzV_F(1)-ATM(ID)%EDGEXzV_F(MNL_F))
!
! NB: The EW in ANg is proportional to tau.(lambda)^2. Hence EW/lambda
! is propotional to tau.lambda
!
                  CALL TORSCL_V2(TA,CHIL,R,TB,TC,DPTH_INDX,'ZERO',TYPE_ATM,L_FALSE)
	          TAU_SOB=T4*TA(DPTH_INDX)+(RONE-T4)*TA(DPTH_INDX-1)
	          T5=(10.0_LDP**T2)*EXP(-T3/2.75_LDP)/(TAU_SOB*T1)
!
	          STRING=' ';  IF(FLAG)STRING=VEC_TRANS_NAME(LINE_INDX)
	          IF(TAU_SOB .GT. TAU_LIM)THEN
	            WRITE(73,'(A,F12.3,ES15.5,F12.3,2F8.2,2(3X,ES10.3),2I6,3X,E10.3,5X,A)')
	1                        ION_ID(ID),ANG_TO_HZ/FL, TAU_SOB,ANG_TO_HZ/FL,
	1                        1.01,1.04,T2,T3,MNL_F,MNUP_F,T5,TRIM(STRING)
	          END IF
	        END IF
	      END DO
	    END IF
	  END DO
	  CLOSE(UNIT=73)
! 
!
! Designed for use with plane-parallel models. We compute the central intesnity of
! each line to deduce line-identifcations.
!
	ELSE IF(XOPT .EQ. 'PLNID')THEN
!
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')' Create line ID file for plane-paraillel models'
	  WRITE(T_OUT,'(A)')' The list is created based on an estimate of the central line intensity'
	  WRITE(T_OUT,'(A)')' The procedure may also work for most O stars -- especially those with weak winds'
	  WRITE(T_OUT,'(A)')' '
	  CALL USR_OPTION(VTURB,'VTURB','10.0','Turbulent velcity for line profile (units km/s)?')
!
	  DEFAULT='0.01'
	  CALL USR_OPTION(DEPTH_LIM,'DEPTH',DEFAULT,
	1        'Output lines whose central depth from the continuum is > LIMIT')
	  CALL USR_OPTION(FLAG,'ADD','F','Add transition name to output file?')
!
	  WRITE(T_OUT,'(A)')' '
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  DEFAULT='LINE_ID'
	  CALL USR_OPTION(FILENAME,'FILE',DEFAULT,'Output file for line IDs')
	  CALL GEN_ASCI_OPEN(73,FILENAME,'UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(73,'(A)')'Species         Lambda      Tau            Lam        B_loc'//
	1                '   E_loc       log gf        hv/kT    NL   NUP     Name'
	  CALL USR_HIDDEN(FLUSH_FILE,'FLUSH','F','Flush outout immediately?')
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
!
! Determine index range encompassing desired wavelength range.
!
	  NU_ST=ANG_TO_HZ/LAM_ST
	  IF(NU_ST .LT. VEC_FREQ(1))THEN
	    NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NL) .GT. NU_ST)NL=NL+1
	  ELSE
	    NL=1
	  END IF
!
	  NU_EN=ANG_TO_HZ/LAM_EN
	  IF(NU_EN .GT. VEC_FREQ(N_LINE_FREQ))THEN
	    NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NUP) .LT. NU_EN)NUP=NUP-1
	  ELSE
	    NUP=N_LINE_FREQ
	  END IF
	  WRITE(6,*)' '
	  WRITE(6,*)'Number of line transitions in the interval is',NUP-NL+1
!
	  FL_SAVE=RZERO
	  DO LINE_INDX=NL,NUP
	    MNL_F=VEC_MNL_F(LINE_INDX)
	    MNUP_F=VEC_MNUP_F(LINE_INDX)
	    FL=VEC_FREQ(LINE_INDX)
!
! Since the opacities are sequentially computed, we dont need to compute them
! at all frequencies. At present we simply compute them every 1000 km/s.
! This saves considerable computing time.
!
	    IF( 2.998E+05_LDP*ABS(FL_SAVE-FL)/FL .GT. 1000.0_LDP)THEN
	      FL_SAVE=FL
	      INCLUDE 'OPACITIES.INC'
!
! Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
! and DDBBDT= dB/dTR .
!
	      T1=HDKT*FL/T(ND)
	      T2=RONE-EMHNUKT(ND)
	      DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	      IF(.NOT. DIF)DBB=RZERO
!
! Adjust the opacities and emissivities for the influence of clumping.
!
	      DO I=1,ND
	        ETA(I)=ETA(I)*CLUMP_FAC(I)
	        CHI(I)=CHI(I)*CLUMP_FAC(I)
	        ESEC(I)=ESEC(I)*CLUMP_FAC(I)
	        CHI_RAY(I)=CHI_RAY(I)*CLUMP_FAC(I)
	        CHI_SCAT(I)=CHI_SCAT(I)*CLUMP_FAC(I)
	      END DO
	    END IF
!
! ID must be set here as it modified by OPACITIES.INC
!
	    ID=VEC_ION_INDX(LINE_INDX)
	    GLDGU=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	    DO I=1,ND
	      T1=ATM(ID)%W_XzV_F(MNUP_F,I)/ATM(ID)%W_XzV_F(MNL_F,I)
	      CHIL(I)=OPLIN*ATM(ID)%AXzV_F(MNL_F,MNUP_F)*(T1*ATM(ID)%XzV_F(MNL_F,I)-GLDGU*ATM(ID)%XzV_F(MNUP_F,I) )
	      ETAL(I)=EMLIN*FREQ*ATM(ID)%AXzV_F(MNUP_F,MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)
	      T2=12.86_LDP*SQRT( T(I)/AT_MASS(SPECIES_LNK(ID))+(VTURB/12.86_LDP)**2 )/C_KMS
	      T2=1.0E-15_LDP/FL/T2/SQRT(PI)
	      CHIL(I)=MAX(T2*CHIL(I)*CLUMP_FAC(I),1.0E-10_LDP)
	      ETAL(I)=ETAL(I)*T2*CLUMP_FAC(I)
	    END DO
	    CALL GET_FLUX_DEFICIT(FLUX_DEFICIT,R,ETA,CHI,CHI_RAY,CHI_SCAT,ESEC,ETAL,CHIL,FL,DBB,ND)
!
            T1=(ANG_TO_HZ/FL)/5000.0_LDP
            T2=LOG10(VEC_OSCIL(LINE_INDX)*ATM(ID)%GXzV_F(MNL_F)*T1)
            T3=HDKT*(ATM(ID)%EDGEXzV_F(1)-ATM(ID)%EDGEXzV_F(MNL_F))
!
	    IF(ABS(FLUX_DEFICIT).GT. DEPTH_LIM)THEN
	      STRING=' ';  IF(FLAG)STRING=VEC_TRANS_NAME(LINE_INDX)
	      T1=ANG_TO_HZ/FL
	      IF(.NOT. ATM(ID)%OBSERVED_LEVEL(MNL_F) .OR. .NOT. ATM(ID)%OBSERVED_LEVEL(MNUP_F))T1=-T1
	      WRITE(73,'(A,F12.3,ES15.5,F12.3,2F8.2,2(3X,ES10.3),2I6,5X,A)')ION_ID(ID),ANG_TO_HZ/FL,
	1                         FLUX_DEFICIT,ANG_TO_HZ/FL,1.01,1.04,T2,T3,MNL_F,MNUP_F,TRIM(STRING)
	      IF(FLUSH_FILE)FLUSH(UNIT=73)
	    END IF
	    IF(MOD(LINE_INDX-NL,MAX(10,(NUP-NL)/10)) .EQ. 0)THEN
	      WRITE(6,'(A,I7,A,I7,A)')' Done',LINE_INDX-NL+1,' of ',NUP-NL+1,' lines'
	    END IF
	  END DO
	  CLOSE(UNIT=73)
!
! Section to allow various line optical deph and wavelength distributions to be examined.
!
! 1. POW: Plots the # of lines in Log(Tau) or Log(Line strength) space. Each uniform logarithmic bin
!         is 0.25 wide in LOG10 space. For Log(tau) space, the line optical depth is normalized by
!         the elctron scatering opacity.
!
! 2. DIST: PLots the the ilogaritmic Sobolev optical depth (radial, or angle averaged) for every
!          line in a given wavelenth interval as a function of the line wavelngth. Lines with
!          a negative optical depth are set to -20.
!
! 3. NV:   Plots the  number of lines with optical depth > TAU_MIN in a given velocity interval.
!
	ELSE IF(XOPT .EQ. 'DIST'
	1        .OR. XOPT .EQ. 'NV'
	1        .OR. XOPT .EQ. 'POW')THEN
!
! Used to count the number lines as a function of optical depth over the entire
! wavelngth interval.
!
	  TAU_GRT_LOGX(:)=RZERO
!
	  IF(DPTH_INDX .LT. 1 .OR. DPTH_INDX .GT. ND)DPTH_INDX=ND/2
	  DEFAULT=WR_STRING(DPTH_INDX)
	  CALL USR_OPTION(DPTH_INDX,'DPTH',DEFAULT,'Depth for line plot')
!
! Ouput summary of model data at depth point
!
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'    R(I)/R*=',R(DPTH_INDX)/R(ND)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       V(I)=',V(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       T(I)=',T(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'      ED(I)=',ED(DPTH_INDX)
          TA(1:ND)=5.65E-15_LDP*ED(1:ND)
          CALL TORSCL(DTAU,TA,R,TB,TC,ND,METHOD,TYPE_ATM)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'  TAU_ES(I)=',DTAU(DPTH_INDX)
!
	  IF(XOPT .EQ. 'POW')THEN
	    CALL USR_OPTION(DO_TAU,'TAU','T','Plot versus TAU or Line strength')
	  ELSE IF(XOPT .EQ. 'DIST')THEN
	    CALL USR_OPTION(DO_TAU,'SRCE','F','Plot source function instead of tau')
	  END IF
	  CALL USR_HIDDEN(WEIGHT_NV,'WGT','T','Weight Tau by [1.0D0-E(-Tau)]?')
	  CALL USR_HIDDEN(MEAN_TAU,'MEAN','F','Use MEAN instead of radial Tau?')
!
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
!
! This makes sure LAM_ST,LAM_END are in A.
!
	  IF(KEV_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/(LAM_ST*KEV_TO_HZ)
	    LAM_EN=ANG_TO_HZ/(LAM_EN*KEV_TO_HZ)
	  ELSE IF(HZ_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/LAM_ST
	    LAM_EN=ANG_TO_HZ/LAM_EN
	  END IF
	  T1=MIN(LAM_ST,LAM_EN); LAM_EN=MAX(LAM_ST,LAM_EN); LAM_ST=T1
!
! TAU_CONSTANT is used to compute the optical depth. Two choices are
! availabe:
!    (1) Radial optical depth. Valid in outer region.
!    (2) Angle avearged optical depth. Might be more appropriated at depth.
!           For an O star it would be better to average over the disk of the
!           star. For a W-R star the disk is not well defined.
!
	  I=DPTH_INDX
	  TAU_CONSTANT=1.0E-15_LDP*C_KMS*R(I)/V(I)
	  IF(MEAN_TAU)THEN
	    T1=SIGMA(DPTH_INDX)
	    T2=SQRT(ABS(T1))
	    IF(T2 .LT. 0.01_LDP)THEN
	      TAU_CONSTANT=TAU_CONSTANT*(RONE-T1/3.0_LDP+T1*T1/5.0_LDP)
	    ELSE IF(T1 .LT. 0)THEN
	      TAU_CONSTANT=0.5_LDP*TAU_CONSTANT*LOG( (RONE+T2)/(RONE-T2) )/T2
	    ELSE
	      TAU_CONSTANT=TAU_CONSTANT*ATAN(T2)/T2
	    END IF
	  ELSE
	    TAU_CONSTANT=TAU_CONSTANT/(RONE+SIGMA(I))
	  END IF
!
! NB: DELTA_TAU is the spacing in LOG10(Tau)
!
	  IF(XOPT .EQ. 'POW')THEN
	    DELTA_TAU=0.25_LDP
	    TAU_BEG=-5.0_LDP
	    NBINS=(6.0_LDP-TAU_BEG)/DELTA_TAU+1
	    DO I=1,NBINS
 	      ZV(I)=TAU_BEG+(I-1)*DELTA_TAU
	    END DO
	    YV(:)=RZERO
	  ELSE IF(XOPT .EQ. 'NV')THEN
	    DELTA_TAU=V(DPTH_INDX)/C_KMS
	    NBINS=LOG(LAM_EN/LAM_ST)/LOG(RONE+DELTA_TAU)+1
	    DO I=1,NBINS
	      ZV(I)=LAM_ST*(RONE+DELTA_TAU)**(I-1)
	    END DO
	    YV(:)=RZERO
	    CALL USR_OPTION(TAU_MIN,'TAU_MIN','1.0D0','Minimum Tau')
	    WRITE(T_OUT,'(1X,A,1P,E11.4)')'Velocity is:',V(DPTH_INDX),' km/s'
	  END IF
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
!
	  CNT=0
!
! Determine index range encompassing desired wavelength range.
!
	  NU_ST=ANG_TO_HZ/LAM_ST
	  IF(NU_ST .LT. VEC_FREQ(1))THEN
	    NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NL) .GT. NU_ST)NL=NL+1
	  ELSE
	    NL=1
	  END IF
!
	  NU_EN=ANG_TO_HZ/LAM_EN
	  IF(NU_EN .GT. VEC_FREQ(N_LINE_FREQ))THEN
	    NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
	    IF(VEC_FREQ(NUP) .LT. NU_EN)NUP=NUP-1
	  ELSE
	    NUP=N_LINE_FREQ
	  END IF
!
	  DO LINE_INDX=NL,NUP
	    MNL_F=VEC_MNL_F(LINE_INDX)
	    MNUP_F=VEC_MNUP_F(LINE_INDX)
	    FL=VEC_FREQ(LINE_INDX)
	    I=DPTH_INDX
!
	    IF(ANG_TO_HZ/FL .GT. LAM_ST .AND. ANG_TO_HZ/FL .LT. LAM_EN)THEN
	      DONE_LINE=.FALSE.
	      IF(XSPEC .EQ. 'ALL')XSPEC=' '
	      DO ID=1,NUM_IONS
	        IF(VEC_SPEC(LINE_INDX) .EQ. ION_ID(ID) .AND.
	1           ( XSPEC .EQ. ' ' .OR. XSPEC .EQ. UC(ION_ID(ID)) .OR.
	1             XSPEC .EQ. SPECIES(SPECIES_LNK(ID)) ))THEN
	          T1=ATM(ID)%W_XzV_F(MNUP_F,I)/ATM(ID)%W_XzV_F(MNL_F,I)
	          CHIL(I)=OPLIN*VEC_OSCIL(LINE_INDX)*(
	1            T1*ATM(ID)%XzV_F(MNL_F,I)-
	1            ATM(ID)%GXzV_F(MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)/
	1            ATM(ID)%GXzV_F(MNUP_F) )
	          ETAL(I)=EMLIN*FL*ATM(ID)%AXzV_F(MNUP_F,MNL_F)*ATM(ID)%XzV_F(MNUP_F,I)
	          DONE_LINE=.TRUE.
	        END IF
	      END DO
!
	      IF(DONE_LINE)THEN
	        IF(XOPT .EQ. 'DIST' .OR. XOPT .EQ. 'NV' .OR. DO_TAU)THEN
	          TAU_SOB=TAU_CONSTANT*CHIL(I)/FL
	          DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	             IF(TAU_SOB .GT. 10.0_LDP**I)THEN
	               TAU_GRT_LOGX(I)=TAU_GRT_LOGX(I)+1
	               EXIT
	             END IF
	          END DO
	          IF(XOPT .EQ. 'POW' .AND. TAU_SOB .GT. 0)THEN
                    I=(LOG10(TAU_SOB)-TAU_BEG)/DELTA_TAU+1
	            IF(I .GT. 0 .AND. I .LE. NBINS)YV(I)=YV(I)+1
	          END IF
	        ELSE
!
! For simplicty, we choos a Doppler velocity of 10 km/s and ignore the
! factor of PI.
!
	          TAU_SOB=1.0E-15_LDP*CHIL(I)/(FL/2.998E+04_LDP)/(6.65E-15_LDP*ED(I))
	          IF(TAU_SOB .GT. 0)THEN
                    I=(LOG10(TAU_SOB)-TAU_BEG+0.5_LDP*DELTA_TAU)/DELTA_TAU+1
	            IF(I .GT. 0 .AND. I .LE. NBINS)YV(I)=YV(I)+1
	            IF(I .GE. NBINS)YV(NBINS)=YV(NBINS)+1
	          END IF
	          DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	             IF(TAU_SOB .GT. 10.0_LDP**I)THEN
	               TAU_GRT_LOGX(I)=TAU_GRT_LOGX(I)+1
	               EXIT
	             END IF
	          END DO
	        END IF
!
	        CNT=CNT+1
	        IF(XOPT .EQ. 'DIST' .AND. DO_TAU)THEN
	          ZV(2*CNT-1)=0.2998E+04_LDP/FL		!Angstroms
	          ZV(2*CNT)=0.2998E+04_LDP/FL
	          YV(2*CNT-1)=-30.0_LDP
	          T1=ETAL(DPTH_INDX)/CHIL(DPTH_INDX)
	          IF(T1 .GT. 1.0E-30_LDP)THEN
	            YV(2*CNT)=LOG10(T1)
	          ELSE IF(T1 .GE. 0.0_LDP)THEN
	            YV(2*CNT)=-30.0_LDP
	          ELSE IF(T1 .GT. -1.0E-30_LDP)THEN
	            YV(2*CNT)=-30.0_LDP
	          ELSE
	            YV(2*CNT)=-60.0_LDP-LOG10(-T1)
	          END IF
	          TAU_SOB=1.0E-15_LDP*CHIL(DPTH_INDX)/(FL/2.998E+04_LDP)/(6.65E-15_LDP*ED(DPTH_INDX))
	          WRITE(66,'(6ES14.4)')FL,0.29979D+04/FL,TAU_SOB,T1,YV(2*CNT-1),YV(2*CNT)
	        ELSE IF(XOPT .EQ. 'DIST')THEN
	          ZV(2*CNT-1)=0.2998E+04_LDP/FL		!Angstroms
	          ZV(2*CNT)=0.2998E+04_LDP/FL
	          YV(2*CNT-1)=-10.0_LDP
	          IF(TAU_SOB .GT. 1.0E-10_LDP)THEN
	            YV(2*CNT)=LOG10(TAU_SOB)
	          ELSE IF(TAU_SOB .GE. 0.0_LDP)THEN
	            YV(2*CNT)=-10.0_LDP
	          ELSE IF(TAU_SOB .GT. -1.0E-10_LDP)THEN
	            YV(2*CNT)=-10.0_LDP
	          ELSE
	            YV(2*CNT)=-10.0_LDP-LOG10(-1.0E+10_LDP*TAU_SOB)
	          END IF
	        ELSE IF(XOPT .EQ. 'NV')THEN
	          IF(TAU_SOB .GT. TAU_MIN)THEN
	            T1=C_KMS/FL/100.0_LDP			!Angstroms
	            T2=LAM_ST*(RONE-0.5_LDP*DELTA_TAU)
	             I=LOG(T1/T2)/LOG(RONE+DELTA_TAU)+1
	            T1=RZERO; T1=MAX(T1,TAU_SOB)
	            IF(WEIGHT_NV)THEN
	              YV(I)=YV(I)+(1.0_LDP-EXP(-T1))
	            ELSE
	              YV(I)=YV(I)+1
	            END IF
	          END IF
	        END IF			!Specifi option
	      END IF			!Line included (e.g. for species)
	    END IF			!Line between LAM_ST and LAM_END
	  END DO			!Loop overlines
!
	  WRITE(T_OUT,'(1X,A,1PE12.4,A)')'ED:   ',ED(DPTH_INDX),'/cm^3'
	  WRITE(T_OUT,'(1X,A,1PE12.4,A)')'Vel:  ',V(DPTH_INDX),'km/s'
	  IF(XOPT .EQ. 'POW' .AND. .NOT. DO_TAU)THEN
	    WRITE(T_OUT,'(1X,A,F9.4)')'/\Line_strength=',DELTA_TAU
	    DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	      T1=10.0_LDP**I
	      WRITE(T_OUT,'(1X,A,I2,A,I9)')'CHIL(Vo)/CHI_ES > 10**(',I,'):',TAU_GRT_LOGX(I)
	    END DO
	  ELSE
	    WRITE(T_OUT,'(1X,A,F9.4)')'/\Tau=',DELTA_TAU
	    WRITE(T_OUT,'(A)')' '
	    WRITE(T_OUT,'(1X,A,I7)')' Total number of lines in interval is',NUP-NL+1
	    WRITE(T_OUT,'(1X,A)')' Number of lines in each decade of optical depth (NOT cummualtive).'
	    WRITE(T_OUT,'(A)')' '
	    DO I=ITAU_GRT_LIM,-ITAU_GRT_LIM,-1
	      T1=10.0_LDP**I
	      WRITE(T_OUT,'(1X,A,I2,A,I9)')'Tau > 10**(',I,'):',TAU_GRT_LOGX(I)
	    END DO
	  END IF
!
	  IF(XOPT .EQ. 'DIST')THEN
	    I=2*CNT
	    IF(HZ_INPUT)THEN
	      DO J=1,I
	         ZV(J)=ANG_TO_HZ/ZV(J)
	      END DO
	    ELSE IF(KEV_INPUT)THEN
	      DO J=1,I
	         ZV(J)=ANG_TO_HZ/ZV(J)/KEV_TO_HZ
	      END DO
	    END IF
	    CALL DP_CURVE(I,ZV,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)','Log \gt\dSob\u',' ',' ')
	  ELSE IF(XOPT .EQ. 'NV')THEN
	    CALL DP_CURVE(NBINS,ZV,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)','N/V',' ',' ')
	  ELSE
!
! We divide by DELTA_TAU so that dN/dLOG_TAU * DELTA_TAU gives the number
! of lines in that bin. The bin is cenetered on the log(tau)=ZV(I).
!
	    DO I=1,NBINS
	      IF(YV(I) .NE. 0)THEN
	        YV(I)=LOG10(YV(I)/DELTA_TAU)
	      ELSE
	        YV(I)=-10.0_LDP
	      END IF
	    END DO
	    XAXIS='Log \gt'
	    YAXIS='Log(dN/dLog\gt)'
	    CALL DP_CURVE(NBINS,ZV,YV)
!
! Draw the canonical slope of 2/3 on the plot.
!
            I=-TAU_BEG/DELTA_TAU+1
	    T1=YV(I)
	    DO I=1,NBINS
	      YV(I)=T1-0.33333_LDP*ZV(I)
	    END DO
	    CALL DP_CURVE(NBINS,ZV,YV)
	  END IF
!
! 
!
	ELSE IF(XOPT .EQ. 'TAUN')THEN
!
	  WRITE(6,'(/,1X,A)')'This option computes the Sobolovev optical depth as a function of level'
	  WRITE(6,'(A)')'at a given depth. It assumes lambda = 1000Ang, f = 1, and ignores stimulated emission.'
	  WRITE(6,'(A,/)')'In general, Tau(Sob) scales as f.lam'
!	  
	  IF(DPTH_INDX .LE. 0 .OR. DPTH_INDX .GT. ND)DPTH_INDX=ND/2
	  WRITE(DEFAULT,*)DPTH_INDX; DEFAULT=ADJUSTL(DEFAULT)
	  CALL USR_OPTION(DPTH_INDX,'DPTH',DEFAULT,'Depth for line plot')
	  CALL USR_OPTION(XAXIS_OPT,'XAXIS','N','What X axis (Level number (def), of Lam (ion)=L)?')
!
! Output summary of model data at depth point
!
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'    R(I)/R*=',R(DPTH_INDX)/R(ND)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       V(I)=',V(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'       T(I)=',T(DPTH_INDX)
          WRITE(T_OUT,'(1X,A,1P,E14.6)')'      ED(I)=',ED(DPTH_INDX)
!
! Factor of 10^(-10) arises since I am assuming Lambda=1000Ang= 10^(-5) cm, since 
! OPLIN.R has the right units, and since I pick up a factor of 10^(-5}
! since V is in km/s.
!
	  I=DPTH_INDX
	  T1=1.00E-10_LDP*OPLIN*R(I)/V(I)
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (XSPEC .EQ. UC(ION_ID(ID)) .OR.
	1                       XSPEC .EQ. SPECIES(SPECIES_LNK(ID)) .OR.
	1                       XSPEC .EQ. 'ALL') )THEN
	        DO NL=1,ATM(ID)%NXzV_F
	           TA(NL)=NL
	           IF(XAXIS_OPT .NE.  'N')TA(NL)=ANG_TO_HZ/ATM(ID)%EDGEXzV_F(NL)
	           TB(NL)=LOG10(ATM(ID)%XzV_F(NL,DPTH_INDX)*T1)
	        END DO
	        CALL DP_CURVE(ATM(ID)%NXzV_F,TA,TB)
	     END IF
	  END DO
	  YAXIS='Log \gt\dSob\u.(1000\A/f/lam)'
	  XAXIS='Level index'
	  IF(XAXIS_OPT .NE.  'N')XAXIS='Ionization wavelength (Ang)'
!
	  IF(I .EQ. 0)THEN
	    WRITE(T_OUT,*)'Error no matching ionization stage'
	  END IF
!
! Plot gf as a function of wavelength.
!
	ELSE IF(XOPT .EQ. 'GF')THEN
!
	   I=0
	   DO LINE_INDX=1,N_LINE_FREQ
	     MNL_F=VEC_MNL_F(LINE_INDX)
	     MNUP_F=VEC_MNUP_F(LINE_INDX)
	     FL=VEC_FREQ(LINE_INDX)
!
	     T1=0
	     IF( XSPEC .EQ. UC(VEC_SPEC(LINE_INDX)) )THEN
	       DO ID=1,NUM_IONS
	         IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	           T1=ABS(VEC_OSCIL(LINE_INDX)*ATM(ID)%GXzV_F(MNL_F))
	         END IF
	       END DO
	     END IF
!
	     IF(T1 .GT. 0)THEN
	       I=I+1
	       ZV(2*I-1)=0.2998E+04_LDP/FL			!Angstroms
	       ZV(2*I)=0.2998E+04_LDP/FL
	       YV(2*I-1)=-10.0_LDP
	       YV(2*I)=LOG10(T1)
	     END IF
	   END DO
!
	   IF(I .EQ. 0)THEN
	     WRITE(T_OUT,*)'Error no matching ionization stage'
	   ELSE
             WRITE(T_OUT,*)'Number of transitions found is',I
	     I=2*I
	     CALL DP_CURVE(I,ZV,YV)
	   END IF
!
! 
	ELSE IF(XOPT .EQ. 'CJBB')THEN
!
! RJ has previously been computed.
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Option plots J/B (def), J/S or J versus X(def) or Log(Tau)'
!
	  CALL USR_OPTION(NORAD,'TPLT','F','Plot against Tau (alt is current X axis)?')
	  CALL USR_OPTION(JONS,'JONS','F','Plot J/S?')
	  JONLY=.FALSE.
	  IF(.NOT. JONS)THEN
	    CALL USR_OPTION(JONLY,'JONLY','F','Plot J only?')
	  END IF
	  IF(.NOT. JONS .AND. .NOT. JONLY)THEN
	    CALL USR_OPTION(SONLY,'SONLY','F','Plot S only?')
	  END IF
!
	  IF(JONS)THEN
	    DO I=1,ND
	      YV(I)=LOG10( RJ(I)*CHI(I)/ETA_WITH_ES(I) )
	    END DO
	    YAXIS='J\dc\u/S'
	  ELSE IF(JONLY)THEN
	    DO I=1,ND
	      YV(I)=LOG10(RJ(I))
	    END DO
	    YAXIS='J\dc\u'
	  ELSE IF(SONLY)THEN
	    DO I=1,ND
	      YV(I)=LOG10(ETA_WITH_ES(I)/CHI(I))
	    END DO
	    YAXIS='S'
	  ELSE
	    DO I=1,ND
	      YV(I)=LOG10( RJ(I)*
	1          ( EXP(HDKT*FREQ/T(I))-RONE )/TWOHCSQ/(FREQ**3)  )
	    END DO
	    YAXIS='J\dc\u/B(T\de\u)'
	  END IF
!
	  IF(NORAD)THEN
	    IF(.NOT. ELEC)THEN
	      DO I=1,ND
	       CHI(I)=CHI(I)-ESEC(I)
	      END DO
	    END IF
	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	    DO I=1,ND
	      ZV(I)=LOG10(TA(I))
	    END DO
	    CALL DP_CURVE(ND,ZV,YV)
	    XAXSAV=XAXIS
	    XAXIS='Log(\gt)'
	  ELSE
	    CALL DP_CURVE(ND,XV,YV)
	    WRITE(T_OUT,*)YV(ND-4),RJ(ND-4),CHI(ND-4),ETA(ND-4),ESEC(ND-4)
	    WRITE(T_OUT,*)YV(ND),RJ(ND),CHI(ND),ETA(ND),ESEC(ND)
	    OPEN(UNIT=37,FILE='RJ_COMP',STATUS='UNKNOWN')
	      WRITE(37,'(5X,A4,6X)')'  YV','  RJ',' CHI',' ETA','ESEC'
	      DO I=1,ND
	        WRITE(37,'(5E15.6)')YV(I),RJ(I),CHI(I),ETA(I),ESEC(I)
	      END DO
	    CLOSE(UNIT=37)
	  END IF
!
! 
! **************************************************************************
! **************************************************************************
!
! To enable a comparison of J (continuum) compared with that computed on a
! finer depth grid.
!
! **************************************************************************
! **************************************************************************
!
	ELSE IF(XOPT .EQ. 'JEXT'
	1  .OR. XOPT .EQ. 'EWEXT')THEN
!
! Compute J. DBB and S1 have been previously evaluated. VT is used
! for dCHI_dR.
!
	  CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1               ZETA,THETA,CHI,VT,JQW,
	1               THICK,DIF,DBB,IC,NC,ND,NP,METHOD)
!
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  DO I=1,ND
	    SOURCE(I)=(ETA(I)+CHI_SCAT(I)*RJ(I))/CHI(I)
	  END DO
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,dCHIdR,HQW,VT,
	1               S1,THICK,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  T1=6.599341_LDP*VT(1)*2.0_LDP		!2 DUE TO 0.5U
          WRITE(T_OUT,'('' Observed Flux is'',1Pe12.4,''Jy'')')T1
          WRITE(LU_LOG,'('' Observed Flux is'',1Pe12.4,''Jy'')')T1
!
	  CALL EXTEND_OPAC(CHIEXT,ETAEXT,ESECEXT,RJEXT,COEF,INDX,NDX
	1                     ,CHI,ETA,CHI_SCAT,RJ,ND)
!
	  DO I=1,NDX
	    ZETAEXT(I)=ETAEXT(I)/CHIEXT(I)
	    THETAEXT(I)=ESECEXT(I)/CHIEXT(I)
	    SOURCEEXT(I)=ZETAEXT(I)+THETAEXT(I)*RJEXT(I)
	    FOLD(I)=RONE/RTHREE
	  END DO
	  S1=SOURCE(1)
!
	  INACCURATE=.TRUE.
	  J=0				!Loop counter
	  DO WHILE(INACCURATE)
	    J=J+1
	    WRITE(T_OUT,'('' Beginning '',I3,''th loop'')')J
	    CALL FQCOMP_V2(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,Q,FEXT,
	1          SOURCEEXT,CHIEXT,ETAEXT,JQWEXT,KQWEXT,DBB,HBC_J,HBC_S,
	1          INBC,IC,THICK,INNER_BND_METH,NCX,NDX,NPX,METHOD)
	    CALL JFEAU_IBC_V2(TA,TB,TC,DTAU,REXT,RJEXT,Q,FEXT,
	1          ZETAEXT,THETAEXT,CHIEXT,DBB,IC,HBC_J,HBC_S,
	1          INBC,THICK,INNER_BND_METH,NDX,METHOD)
	    S1=ZETAEXT(1)+THETAEXT(1)*RJEXT(1)
	    INACCURATE=.FALSE.
	    T1=RZERO
	    DO I=1,NDX
	     T1=MAX(ABS(FOLD(I)-FEXT(I)),T1)
	     FOLD(I)=FEXT(I)
	     SOURCEEXT(I)=ZETAEXT(I)+THETAEXT(I)*RJEXT(I)
	    END DO
	    IF(T1 .GT. 1.0E-05_LDP)INACCURATE=.TRUE.
	    WRITE(T_OUT,'('' Maximum fractional change is'',1PE11.4)')T1
	  END DO
!
	  IF(XOPT .EQ. 'JEXT')THEN
!
! Will use AV for the new RJ
!
	    CALL UNGRID(AV,ND,RJEXT,NDX,GRID)
!
	    S1=(ETA(1)+AV(1)*CHI_SCAT(1))/CHI(1)
	    DO I=1,ND
	      SOURCE(I)=(ETA(I)+CHI_SCAT(I)*AV(I))/CHI(I)
	    END DO
	    CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,dCHIdR,HQW,VT,
	1                 S1,THICK,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	    T1=6.599341_LDP*VT(1)*2.0_LDP		!2 DUE TO 0.5U
            WRITE(T_OUT,'('' Observed Flux[Ext] is'',1Pe12.4,''Jy'')')T1
            WRITE(LU_LOG,'('' Observed Flux[Ext] is'',1Pe12.4,''Jy'')')T1
!
! With this definition for YV, can never get an error worse than 200%.
! An error of 100% implies that RJ and RJEXT differ by a factor of 3.
!
	    DO I=1,ND
	      YV(I)=200.0_LDP*(AV(I)-RJ(I))/(AV(I)+RJ(I))
	    END DO
	    YAXIS='% error'
!
	    CALL USR_HIDDEN(ELEC,'LOGE','F','Log of error?')
	    IF(ELEC)THEN
	      DO I=1,ND
	        YV(I)=LOG10(YV(I))
	      END DO
	      YAXIS='Log(% error)'
	    END IF
!
	    CALL DP_CURVE(ND,XV,YV)
	  ELSE IF(XOPT .EQ. 'EWEXT')THEN
!
	    DO I=1,ND
	      SOURCE(I)=ZETA(I)+THETA(I)*RJEXT(GRID(I))
	    END DO
!
! TA is used for the line flux. Integral of TA log(r) is
! the line EW.
!
	    CALL SOBEW_GRAD_V2(SOURCE,CHI,CHI_SCAT,CHIL,ETAL,
	1              V,SIGMA,R,P,FORCE_MULT,LUM,
	1              JQW,HQW,TA,T1,S1,
	1              FREQ,'HOLLOW',DBB,IC,THICK,.FALSE.,NC,NP,ND,METHOD)
!
	    CALL DP_CURVE(ND,XV,FORCE_MULT)
	    YAXIS='Force Multiplier'
!
	    OPEN(UNIT=18,FILE='DSOB_FORCE_MULT',STATUS='UNKNOWN')
	      WRITE(18,'(3X,A1,10X,A1,15X,A1,13X,A1)')'I','R','V','M'
              DO I=1,ND
                WRITE(18,'(1X,I3,3X,3ES14.6)')I, R(I),V(I),FORCE_MULT(I)
              END DO
            CLOSE(UNIT=18)
!
! S1 is the continuum flux in Jy for an object at 1kpc.
! T1 is the line equivalent width in Angstroms.
! T3 is the line flux in ergs/cm^2/s
!
	    T2=LAMVACAIR(FREQ)		!Wavelength(Angstroms)
	    T3=T1*S1*1.0E-23_LDP*FREQ*1.0E+15_LDP/T2
	    T4=T3*4*PI*(3.0856E+21_LDP)**2/LUM_SUN()
	    WRITE(LU_NET,40008)T1,S1,T2,T3,T4
	    WRITE(T_OUT,'(A)')RED_PEN
	    WRITE(T_OUT,40008)T1,S1,T2,T3,T4
40008	    FORMAT(1X,'EW =',ES10.3,' Ang',3X,'I =',ES10.3,' Jy',3X,
	1      'Lam =',ES11.4,' Ang',3X,'Line flux=',ES10.3,' ergs/cm^2/s =',ES10.3,' Lsun')
	    WRITE(T_OUT,'(A)')DEF_PEN
	  END IF
!
! 
! **************************************************************************
! **************************************************************************
!
! Comparison of line source function with Jc and B.
!
! **************************************************************************
! **************************************************************************
!
	ELSE IF(XOPT .EQ. 'SRCEBB')THEN
!
! Assumes line opacity and emissivity have been computed in the set up.
!
	  DO I=1,ND
	    IF(CHIL(I).LT. 1.0E-30_LDP)THEN
	      YV(I)=-30
	    ELSE
	      YV(I)=LOG10(  ETAL(I)/CHIL(I)*
	1          ( EXP(HDKT*FREQ/T(I))-RONE )/TWOHCSQ/(FREQ**3)  )
	    END IF
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='S/B(T\de\u)'
!
	ELSE IF(XOPT .EQ. 'SRCEJC')THEN
!
! RJ has previously been computed, as have the Line opacity and emissivity.
!
	  DO I=1,ND
	    IF(CHIL(I).LT. 1.0E-30_LDP)THEN
	      YV(I)=-30
	    ELSE
	      YV(I)=LOG10( ETAL(I)/CHIL(I)/RJ(I) )
	    END IF
	  END DO
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='S/J\dc\u'
!
	ELSE IF(XOPT .EQ. 'SRCE')THEN
!
! Assumes line opacity and emissivity have been computed in the set up.
!
	  CALL USR_OPTION(ELEC,'LIN','T','Plot S instead of Log S')
	  IF(ELEC)THEN
	    DO I=1,ND
	      IF( ABS(CHIL(I)) .LT. 1.0E-50_LDP)THEN
	         YV(I)=-30
	      ELSE
	        YV(I)=ETAL(I)/CHIL(I)
	      END IF
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='S'
	  ELSE
	    DO I=1,ND
	      IF(CHIL(I).LT. 1.0E-30_LDP)THEN
	         YV(I)=-30
	      ELSE
	        YV(I)=LOG10( ETAL(I)/CHIL(I) )
	      END IF
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log S'
	  END IF
!
	ELSE IF(XOPT .EQ. 'BB' .OR. XOPT .EQ. 'DBBDT')THEN
	  I=ND/2
	  DEFAULT=WR_STRING(I)
	  VALID_VALUE=.FALSE.
	  DO WHILE(.NOT. VALID_VALUE)
	    CALL USR_OPTION(I,'DEPTH',DEFAULT,'Depth for plotting BB')
	    IF(I .GE. 1 .AND. I .LE. ND)THEN
	       WRITE(T_OUT,'(A)')BLUE_PEN
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'     Radius at depth',I,' is',R(I)
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'   Velocity at depth',I,' is',V(I)
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'Temperature at depth',I,' is',T(I)
	       WRITE(T_OUT,'(A)')DEF_PEN
	       VALID_VALUE=.TRUE.
	       TEMP=T(I)
	    END IF
	  END DO
	  CALL USR_OPTION(ELEC,'BLAM','T','Compute Blam')
	  K=MIN(4000,N_PLT_MAX)
	  T1=0.01_LDP*C_KMS/10.0_LDP
	  T2=10**(4.0_LDP/(K-1))
	  XV(1)=T1
	  IF(XOPT .EQ. 'BB')THEN
	    DO I=1,K
	      IF(I .NE. 1)XV(I)=XV(I-1)/T2
	      T3=EXP(-HDKT*XV(I)/TEMP)
	      YV(I)=TWOHCSQ*(XV(I)**3)*T3/(RONE-T3)
	    END DO
	    YAXIS='B\d\gn\u(\gn)(ergs/cm\u2\d/Hz/str)'
	    IF(ELEC)THEN
	      YV(1:K)=1.0E+06_LDP*YV(1:K)*XV(1:K)*XV(1:K)/C_KMS
	      YAXIS='B\d\gl\u(\gl)(Gergs/cm\u2\d/\gA/str)'
	    END IF
	  ELSE
	    DO I=1,K
	      IF(I .NE. 1)XV(I)=XV(I-1)/T2
	      T4=HDKT*XV(I)/TEMP
	      T3=EXP(-T4)
	      YV(I)=1.0E-04_LDP*TWOHCSQ*T4*(XV(I)**3)*T3/((RONE-T3)**2)/TEMP
	    END DO
	    YAXIS='dB\d\gn\u(\gn)/dT(Mergs/cm\u2\d/\gHZ/str/K)'
	    IF(ELEC)THEN
	      YV(1:K)=1.0E+09_LDP*YV(1:K)*XV(1:K)*XV(1:K)/C_KMS
	      YAXIS='dB\d\gl\u(\gl)/dT(Mergs/cm\u2\d/\gA/str/K)'
	    END IF
	  END IF
	  XV(1:K)=0.01_LDP*C_KMS/XV(1:K)
	  CALL DP_CURVE(K,XV,YV)
	  XAXIS='\gl(\A)'
!
! 
!

	ELSE IF(XOPT .EQ. 'DION')THEN
	  TYPE=' '
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(TRIM(ION_ID(ID))) ) THEN
	      CALL LOGVEC(ATM(ID)%DXzV_F,YV,ND)
	      CALL DP_CURVE(ND,XV,YV)	
	      IF(ATM(ID)%ZXzV .LE. 9)WRITE(TYPE,'(I1)')NINT(ATM(ID)%ZXzV)
	      IF(ATM(ID)%ZXzV .GE. 10)WRITE(TYPE,'(I2)')NINT(ATM(ID)%ZXzV)
	      YAXIS='Log(H\u'//TRIM(TYPE)//'\d+)'
	      YV(1:ND)=ATM(ID)%DXzV(1:ND)
	    END IF
	  END DO
!
! 
!
	ELSE IF(XOPT .EQ. 'SCL')THEN
	  CALL USR_OPTION(T1,'SCL_FAC',' ',
	1   'Factor to multiply POPULATIONS by (<1)')
	  CALL USR_OPTION(K,'DEPTH',' ',
	1   'Depth index below which [for 1 to I] populations are reduced.')
!
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. SPECIES(ISPEC))THEN
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO I=1,K
	            DO J=1,ATM(ID)%NXzV_F
	              ATM(ID)%XzV_F(J,I)=ATM(ID)%XzV_F(J,I)*T1
	            END DO
	            ATM(ID)%DXzV_F(I)=ATM(ID)%DXzV_F(I)*T1
	            IF(.NOT. ATM(ID+1)%XzV_PRES)ATM(ID+1)%XzV_F(1,I)=
	1                ATM(ID)%XzV_F(1,I)*T1
	          END DO
	        END IF
	      END DO
	    END IF
	  END DO
! 
!
	ELSE IF(XOPT .EQ. 'AV')THEN
	  CALL USR_OPTION(K,'DEPTH',' ','Depth index for averaging')
	  DO ID=1,NUM_IONS
            IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%NXzV_F
	        ATM(ID)%XzV_F(J,K)=SQRT(ATM(ID)%XzV_F(J,K-1)*ATM(ID)%XzV_F(J,K+1) )
	      END DO
	      IF(ID .NE. 1)THEN
	        IF(ATM(ID-1)%XzV_PRES)THEN
                  ATM(ID-1)%DXzV_F(K)=SQRT(ATM(ID-1)%DXzV_F(K-1)*ATM(ID-1)%DXzV_F(K+1))
	        END IF
	      END IF
	      IF(ID .NE. NUM_IONS)THEN
	        IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
                  ATM(ID)%DXzV_F(K)=SQRT(ATM(ID)%DXzV_F(K-1)*ATM(ID)%DXzV_F(K+1))
	        END IF
	      END IF
	      EXIT
	    END IF
	  END DO
!
	ELSE IF(XOPT .EQ. 'REP')THEN
	  CALL USR_OPTION(K,'DEPTH',' ','Depth to replace')
	  I=K+1
	  CALL USR_OPTION(I,'DEPTH',' ','Replacement depth')
	  DO ID=1,NUM_IONS
            IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%NXzV_F
	        ATM(ID)%XzV_F(J,K)=ATM(ID)%XzV_F(J,I)*MASS_DENSITY(K)/MASS_DENSITY(I)
	      END DO
	      IF(ID .NE. 1)THEN
	        IF(ATM(ID-1)%XzV_PRES)THEN
                  ATM(ID-1)%DXzV_F(K)=ATM(ID-1)%DXzV_F(I)*MASS_DENSITY(K)/MASS_DENSITY(I)
	        END IF
	      END IF
	      IF(ID .NE. NUM_IONS)THEN
	        IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
                  ATM(ID)%DXzV_F(K)=ATM(ID)%DXzV_F(I)*MASS_DENSITY(K)/MASS_DENSITY(I)
	        END IF
	      END IF
	      EXIT
	    END IF
	  END DO
! 
!
	ELSE IF (XOPT .EQ. 'FIXT')THEN
!
	  CALL USR_OPTION(ELEC,'DT','F','Set T values at certain depths?')
	  IF(ELEC)THEN
	    CALL USR_OPTION(LEV,ITEN,ITWO,'OWIN','0,0,0,0,0,0,0,0,0,0',
	1        'Depth section to be changed [2: lower and upper limit]')
            DEFAULT='1.0'
	    DO I=LEV(1),LEV(2)
	      CALL USR_OPTION(T1,'TEMP',DEFAULT,'Revised temperature')
	      T(I)=T1
	      DEFAULT=WR_STRING(T1)
	    END DO
	  ELSE
	    CALL USR_HIDDEN(T1,'SCALE','1.0',
	1      'Optical Depth below which T is held approximately constant')
!
! As we don't require the old T, we can overite it straight away.
!
	    DO I=1,ND
	      IF(T1 .LE. 0)THEN
	        T1=RONE
	      ELSE
	        T1=TAUROSS(I)/T1
	        T1=RONE-EXP(-T1)
	      END IF
	      T(I)=T1*TGREY(I)+(RONE-T1)*T(I)
	    END DO
!

	    DO ID=1,NUM_IONS
	      CALL UPDATE_POPS(ATM(ID)%XzV_F,ATM(ID)%XzVLTE_F,
	1                    ATM(ID)%NXzV_F,ND,T,TA)
	    END DO
	    DO I=1,ND
	      T(I)=TA(I)
	    END DO
	  END IF
!
! Compute LTE populations
!
!            INCLUDE 'EVAL_LTE_FULL.INC'
!
	  ROSS=.FALSE.
	  GREY_COMP=.FALSE.
!
! 
	ELSE IF(XOPT .EQ. 'XCOOL')THEN
!
	   WRITE(6,*)' '
	   WRITE(6,*)' Compares X-ray cooling time and flow times assuming smooth flow'
	   WRITE(6,*)' '
!
! The factor 0.375 = 1.5/4 (3/2 kT and 4 for the shock compression).
!
	   T2=BOLTZMANN_CONSTANT()
	    DO I=1,ND
	      TA(I)=0.375_LDP*T2*(ED(I)+POP_ATOM(I))*1.0E+28_LDP/ED(I)/POP_ATOM(I)/CLUMP_FAC(I)
	      TB(I)=1.0E+05_LDP*R(I)/V(I)
	    END DO
	    CALL DP_CURVE(ND,XV,TA)
	    CALL DP_CURVE(ND,XV,TB)
	    YAXIS='t\dX\u(T/10\u6\d)(\gL/10\u22\d) s; t(flow)'
	ELSE IF(XOPT .EQ.'XRAY')THEN
	  XRAYS=.NOT. XRAYS
	  IF(XRAYS)THEN
	    WRITE(T_OUT,*)'Xray opacities now included'
	    CALL USR_HIDDEN(FILL_FAC_XRAYS,'FILL','1.0',
	1	 'X-ray filling factor')
	    CALL USR_HIDDEN(T_SHOCK,'TSHOCK','300',
	1	 'Shock temperature of region producing X-rays')
	    CALL USR_HIDDEN(V_SHOCK,'VSHOCK','200 ',
	1	 'Minimum velocity for X-ray production')
	  ELSE
	    WRITE(T_OUT,*)'Xray opacities NOT included'
	  END IF
!
! Switch on/of level dissolution.
!
	ELSE IF(XOPT .EQ. 'DIS')THEN
	  LEVEL_DISSOLUTION= .NOT. LEVEL_DISSOLUTION
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,LEVEL_DISSOLUTION,ND)
	  INCLUDE 'EVAL_LTE_FULL.INC'
	  IF(LEVEL_DISSOLUTION)THEN
	    WRITE(T_OUT,*)'Effect of level dissolution on opacities '//
	1	 'now included'
	  ELSE
	    WRITE(T_OUT,*)'Effect of level dissolution on opacities '//
	1	 'NOT included'
	  END IF
	  WRITE(T_OUT,*)'Warning --- level dissolution option should match' ,
	1    'that of main code to avoid inconsisteancies.'
!
! 
! Computes the optical depth of lines of a given ionization stage and species at a
! given depth in the atmosphere. Either the radial or tangential Sobolev optical
! depth is used.
!
	ELSE IF(XOPT .EQ. 'ARAT')THEN
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')'Plotting line optical depth versus wavelength'
	  WRITE(6,'(A)')'Lines with different A/SUM ranges ar eplotted in different colors'
	  WRITE(6,'(A)')BLACK_PEN
!
	  XAXSAV=XAXIS
          CALL BRANCH_RAT(OMEGA_F,XV,YV,XAXIS,YAXIS,XSPEC,N_MAX,N_PLT_MAX,ND)
!
! Computes the optical depth of lines of a given ionization stage and species at a
! given depth in the atmosphere. Either the radial or tangential Sobolev optical
! depth is used.
!
	ELSE IF(XOPT .EQ. 'LTAU')THEN
!
	  XV_SAV=XV
	  I=ND/2
	  DEFAULT=WR_STRING(I)
	  VALID_VALUE=.FALSE.
	  DO WHILE(.NOT. VALID_VALUE)
	    CALL USR_OPTION(I,'DEPTH',DEFAULT,'Depth for plotting TAUL')
	    IF(I .GE. 1 .AND. I .LE. ND)THEN
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'     Radius at depth',I,' is',R(I)
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'   Velocity at depth',I,' is',V(I)
	       WRITE(T_OUT,'(A,I4,A,ES12.4)')'Temperature at depth',I,' is',T(I)
	       VALID_VALUE=.TRUE.
	    END IF
	  END DO
	  CALL USR_HIDDEN(FLAG,'LOGX','F','Logarithmic in Wavelength?')
	  CALL USR_OPTION(RADIAL_TAU,'RD_TAU','TRUE',
	1      'Use radial (alt. is TANGENTIAL) direction to evaluate the Sobolev optical depth')
	  CALL USR_OPTION(LINE_STRENGTH,'LS','F','Plot line strength')
	  TAU_MIN=1.0E-04_LDP
	  DEFAULT=WR_STRING(TAU_MIN)
	  CALL USR_HIDDEN(TAU_MIN,'TAU_MIN',DEFAULT,'Minimum Tau')
	  TAU_MIN=LOG10(TAU_MIN)
!
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  CALL USR_OPTION(WR_LINE,'WR_LINE','F','Write line?')
!
! This makes sure LAM_ST,LAM_END are in A.
!
	  IF(KEV_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/(LAM_ST*KEV_TO_HZ)
	    LAM_EN=ANG_TO_HZ/(LAM_EN*KEV_TO_HZ)
	  ELSE IF(HZ_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/LAM_ST
	    LAM_EN=ANG_TO_HZ/LAM_EN
	  END IF
	  T1=MIN(LAM_ST,LAM_EN); LAM_EN=MAX(LAM_ST,LAM_EN); LAM_ST=T1
	  FREQ_ST=ANG_TO_HZ/LAM_ST; FREQ_EN=ANG_TO_HZ/LAM_EN
!
	  IF(LINE_STRENGTH)THEN
!
! Vdop=10 km/s and ignoring SQRT(PI)
!
	    WRITE(6,*)'Opacity at line center(for Vdop=10km/s) * SQRT(PI) normalized be e.s. opacity'
	    TAU_CONSTANT=1.0E-15_LDP*OPLIN*2.998E+04_LDP/(6.65E-15_LDP*ED(I))
	  ELSE IF(V(I) .LE. 20.0_LDP)THEN
	    WRITE(6,*)'Using modified static optical depth'
	    TAU_CONSTANT=OPLIN*1.6914E-11_LDP		!Assumes Vdop=10 km/s
	  ELSE
	    TAU_CONSTANT=OPLIN*R(I)*2.998E-10_LDP/V(I)
	    IF(RADIAL_TAU)TAU_CONSTANT=TAU_CONSTANT/(RONE+SIGMA(I))
	  END IF
!
! Compute line opacity and emissivity.
!
	  J=0
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (XSPEC .EQ. UC(ION_ID(ID)) .OR.
	1                       XSPEC .EQ. SPECIES(SPECIES_LNK(ID)) .OR.
	1                       XSPEC .EQ. 'ALL') )THEN
	      DO NL=1,ATM(ID)%NXzV_F
	        DO NUP=NL+1,ATM(ID)%NXzV_F
	          WRITE(25,*)NL,NUP,ATM(ID)%AXzV_F(NL,NUP)
	          IF(ATM(ID)%AXzV_F(NL,NUP) .NE. 0)THEN
	            J=J+1
	            T1=ATM(ID)%W_XzV_F(NUP,I)/ATM(ID)%W_XzV_F(NL,I)
	            FREQ=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	            GLDGU=ATM(ID)%GXzV_F(NL)/ATM(ID)%GXzV_F(NUP)
	            IF(FLAG)THEN
                      XV(J)=LOG10( ANG_TO_HZ/FREQ )
	            ELSE
                      XV(J)=ANG_TO_HZ/FREQ
	            END IF
	            T2=ATM(ID)%AXzV_F(NL,NUP)*(T1*ATM(ID)%XzV_F(NL,I)-GLDGU*ATM(ID)%XzV_F(NUP,I))
	            WRITE(25,*)NL,NUP,T1,T2
	            WRITE(25,*)J,XV(J),T2
	            IF(T2 .NE. 0 .AND. FREQ .LE. FREQ_ST .AND. FREQ .GE.  FREQ_EN)THEN
	              IF(V(I) .LE. 20.0_LDP .AND. .NOT. LINE_STRENGTH)THEN
	                YV(J)=RZERO
	                T3=T2
	                DO K=I,2,-1
	                  T2=T3
	                  T1=ATM(ID)%W_XzV_F(NUP,K)/ATM(ID)%W_XzV_F(NL,K)
	                  T3=T1*ATM(ID)%XzV_F(NL,K)-GLDGU*ATM(ID)%XzV_F(NUP,K)
	                  YV(J)=YV(J)+(R(K-1)-R(K))*(T2+T3)
	                  IF(V(K+1) .GT. 20.0_LDP)EXIT
	                END DO
	                YV(J)=LOG10( 0.5_LDP*ATM(ID)%AXzV_F(NL,NUP)*ABS(YV(J))*TAU_CONSTANT/FREQ)
	              ELSE
	                YV(J)=LOG10( ABS(T2)*TAU_CONSTANT/FREQ )
	              END IF
	              IF(WR_LINE)THEN
	                WRITE(40,'(F12.4,ES12.4,2X,A10,2I6)')XV(J),YV(J),ION_ID(ID), NL, NUP
	              END IF
	              IF(YV(J) .LE. TAU_MIN)THEN
	                J=J-1
	              ELSE
	               FOUND=.TRUE.
	              END IF
	            ELSE
	              J=J-1
	            END IF
	            IF(J .GT. 0)WRITE(25,*)J,XV(J),T2
	          END IF
	        END DO
	      END DO
	    END IF
	  END DO
	  IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)' Invalid population type or species unavailable.'
	    GOTO 1
	  END IF
!
	  XAXSAV=XAXIS
	  IF(FLAG)THEN
	    XAXIS='Log(\gl(\A))'
	  ELSE
	    XAXIS='\gl(\gA)'
	  END IF
	  YAXIS='Log(\gt)'
!
	  IF(J .NE. 0)WRITE(T_OUT,*)J,' lines plotted'
	  IF(J  .NE. 0)CALL DP_CURVE(J,XV,YV)
	  XV=XV_SAV
!
!`
	ELSE IF(XOPT .EQ. 'VLTAU')THEN
!
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  NU_ST=ANG_TO_HZ/LAM_ST
	  NU_EN=ANG_TO_HZ/LAM_EN
	  T1=NU_ST
	  NU_ST=MIN(NU_ST,NU_EN)
	  NU_EN=MAX(T1,NU_EN)
!
	  DO I=1,ND
	    TA(I)=OPLIN*R(I)*2.998E-10_LDP/V(I)
	  END DO
!
	  J=0
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (XSPEC .EQ. UC(ION_ID(ID)) .OR. XSPEC .EQ. 'ALL') )THEN
	      DO NL=1,ATM(ID)%NXzV_F
	        DO NUP=NL+1,ATM(ID)%NXzV_F
	          FREQ=ATM(ID)%EDGEXzV_F(NL)-ATM(ID)%EDGEXzV_F(NUP)
	          GLDGU=ATM(ID)%GXzV_F(NL)/ATM(ID)%GXzV_F(NUP)
	          IF(ATM(ID)%AXzV_F(NL,NUP) .NE. 0 .AND. FREQ .GT. NU_ST .AND. FREQ .LT. NU_EN)THEN
	            J=J+1
	            IF(FLAG)THEN
                      XV(J)=LOG10( ANG_TO_HZ/FREQ )
	            ELSE
                      XV(J)=ANG_TO_HZ/FREQ
	            END IF
	            ELEC=.FALSE.
	            DO I=1,ND
	              T1=ATM(ID)%W_XzV_F(NUP,I)/ATM(ID)%W_XzV_F(NL,I)
	              T2=ATM(ID)%AXzV_F(NL,NUP)*(T1*ATM(ID)%XzV_F(NL,I)-GLDGU*ATM(ID)%XzV_F(NUP,I))
	              T2=T2*TA(I)/FREQ
	              IF(ELEC .AND. T2 .LT. RONE)THEN
	                ZV(J)=(YV(J)-V(I))/2
	                YV(J)=(YV(J)+V(I))/2
	                EXIT
	              ELSE IF(T2 .GT. RONE .AND. .NOT. ELEC)THEN
	                ELEC=.TRUE.
	                YV(J)=V(I)
	              ELSE IF(ELEC .AND. I .EQ. ND)THEN
	                ZV(J)=(YV(J)-V(I))/2
	                YV(J)=(YV(J)+V(I))/2
	              ELSE IF(I .EQ. ND)THEN
	                J=J-1
	              END IF
	            END DO
	          END IF
	        END DO
	      END DO
	    END IF
	  END DO
!
	  XAXSAV=XAXIS
	  IF(FLAG)THEN
	    XAXIS='Log(\gl(\A))'
	  ELSE
	    XAXIS='\gl(\gV)'
	  END IF
	  YAXIS='V(km/s)'
	  IF(J .NE. 0)WRITE(T_OUT,*)J,' lines plotted'
	  IF(J  .NE. 0)CALL DP_CURVE_AND_ER(J,XV,YV,ZV,' ')
!
!
! This section creates an image of the partial opacities (i.e., the opacities due
! to each ionizaion stage) at a given frequency as a function of depth.
!
	ELSE IF(XOPT .EQ. 'DCHI')THEN
	  CALL USR_OPTION(FREQ,'LAM','0.0',FREQ_INPUT)
	  IF(FREQ .EQ. 0)GOTO 1
	  IF(KEV_INPUT)THEN
	    FREQ=FREQ*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    FREQ=ANG_TO_HZ/FREQ
	  END IF
	  FL=FREQ
!
	  INCLUDE 'PAR_OPACITIES.INC'
	  CALL USR_OPTION(ELEC,'ELEC','F','Include electron scattering?')
	  IF(ELEC)THEN
	    CALL ESOPAC(CHI,ED,ND)
	  ELSE
	    CHI(1:ND)=RZERO
	  END IF
!
	  DO ID=1,NUM_IONS
	    CHI(1:ND)=CHI(1:ND)+CHI_PAR(1:ND,ID)
	  END DO
	  WRITE(6,*)'Computed the total opacity'
!
	  DO ID=1,NUM_IONS
	    CHI_PAR(:,ID)=CHI_PAR(:,ID)/CHI(1:ND)
            YMAPV(ID)=ID
	  END DO
	  CALL USR_OPTION(ELEC,'Depth','T','Plot against depth index')
	  IF(ELEC)THEN
	    DO I=1,ND
	      XMAPV(I)=I
	    END DO
	  ELSE
	    DO I=1,ND
	      WRITE(6,*)'I=',I,ED(I)
	      XMAPV(I)=LOG10(ED(I))
	    END DO
	  END IF
!
! Remove rows containing all zero's.
!
	  ETA_PAR(:,:)=CHI_PAR(:,:)
	  J=0
	  LAST_NON_ZERO=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(SUM(ETA_PAR(:,ID)) .NE. 0)THEN
	      J=J+1
              CHI_PAR(:,J)=ETA_PAR(:,ID)
	      LAST_NON_ZERO=.TRUE.
	    ELSE IF(LAST_NON_ZERO)THEN
	      J=J+1
              CHI_PAR(:,J)=ETA_PAR(:,ID)
	      LAST_NON_ZERO=.FALSE.
	      WRITE(6,*)J,ION_ID(ID)
	    ELSE
	    END IF
	  END DO
!
! ETA_PAR is passed, but should not be imaged.
!
	  CALL MAP_PLOT(CHI_PAR,ETA_PAR,XMAPV,YMAPV,ND,J,
	1               " "," "," "," ")
!
!
!
! This section does one of two things:
!           1. Creates an image of the fractional contribution to the opacity
!              and emissivity of a single species as a function of depth and
!              wavelength.
!           2. When no species is specified, the maximum fractional contribution
!               of each species at ANY depth is ouput to the terminal
!           3. Plot species dependent contributions at a given depth.
!
	ELSE IF(XOPT .EQ. 'MCHI')THEN
	  CALL USR_OPTION(LAM_ST,'LAM_ST','50.0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAM_END','10000.0',FREQ_INPUT)
	  CALL USR_OPTION(NLAM,'NLAM','1000','# of wavlengths')
	  CALL USR_OPTION(ELEC,'ELEC','F','Include electron scattering?')
	  CALL USR_OPTION(DO_KAP,'DO_KAP','F','Plot species dependent contributions?')
	  IF(DO_KAP)THEN
	    CALL USR_OPTION(DPTH_INDX,'DPTH','30','Depth index?')
	  END IF
!
	  LOC_ION_ID='$$'
	  DO WHILE(LOC_ION_ID .EQ. '$$' .AND. .NOT. DO_KAP)
	    CALL USR_OPTION(LOC_ION_ID,'ION',' ','Ion ID (eg HI)')
	    DO ID=1,NUM_IONS
	      IF(LOC_ION_ID .EQ. ION_ID(ID) .OR. LOC_ION_ID .EQ. ' ')THEN
	         EXIT
	      END IF
	      IF(ID .EQ. NUM_IONS)THEN
	        WRITE(6,*)'Invalid ION_ID',LOC_ION_ID
	        LOC_ION_ID='$$'
	      END IF
	    END DO
	  END DO
!
	  IF(KEV_INPUT)THEN
	    NU_ST=LAM_ST*KEV_TO_HZ
	    NU_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    NU_EN=ANG_TO_HZ/LAM_ST
	    NU_ST=ANG_TO_HZ/LAM_EN
	  ELSE
	    NU_ST=LAM_ST
	    NU_EN=LAM_EN
	  END IF
	  DEL_NU=10.0_LDP**( LOG10(NU_EN/NU_ST)/(NLAM-1) )
!
	  IF(LOC_ION_ID .NE. ' ')THEN
	    ALLOCATE (CHI_LAM(NLAM,ND))
	    ALLOCATE (ETA_LAM(NLAM,ND))
	    ALLOCATE (LAM_VEC(NLAM));     LAM_VEC(:)=RZERO
	  ELSE IF(DO_KAP)THEN
	    ALLOCATE (CHI_LAM(NLAM,NUM_IONS))
	    ALLOCATE (ETA_LAM(NLAM,NUM_IONS))
	    ALLOCATE (LAM_VEC(NLAM));     LAM_VEC(:)=RZERO
	  ELSE
	    ALLOCATE (CHI_TOT_LAM(NUM_IONS)); CHI_TOT_LAM(:)=RZERO
	    ALLOCATE (LAM_VEC(NUM_IONS));     LAM_VEC(:)=RZERO
	  END IF
!
	  LAST_NON_ZERO=.FALSE.
	  NION=0
	  DO ML=1,NLAM
	    WRITE(6,*)'ML=',ML
	    FREQ=NU_ST*(DEL_NU**(ML-1))
	    FL=FREQ
!
	    INCLUDE 'PAR_OPACITIES.INC'
	    IF(ELEC)THEN
	      CALL ESOPAC(CHI,ED,ND)
	    ELSE
	      CHI(1:ND)=RZERO
	    END IF
	    ETA(1:ND)=RZERO
!
	    DO ID=1,NUM_IONS
	      CHI(1:ND)=CHI(1:ND)+CHI_PAR(1:ND,ID)
	      ETA(1:ND)=ETA(1:ND)+ETA_PAR(1:ND,ID)
	    END DO
!
! NB: First dimension of CHI_PAR is ND.
!
	    IF(.NOT. DO_KAP)THEN
	      DO ID=1,NUM_IONS
	        CHI_PAR(:,ID)=CHI_PAR(:,ID)/CHI(1:ND)
	        ETA_PAR(:,ID)=ETA_PAR(:,ID)/ETA(1:ND)
                YMAPV(ID)=ID
	      END DO
	    END IF
!
	    IF(LOC_ION_ID .EQ. ' ')THEN
	      DO ID=1,NUM_IONS
	        T1=MAXVAL(CHI_PAR(:,ID))
	        T2=MAXVAL(ETA_PAR(:,ID))
	        CHI_TOT_LAM(ID)=MAX(T1,CHI_TOT_LAM(ID))
	        LAM_VEC(ID)=MAX(T2,LAM_VEC(ID))
	      END DO
	    ELSE IF(DO_KAP)THEN
	       LAM_VEC(ML)=LOG10(ANG_TO_HZ/FREQ)
	       DO ID=1,NUM_IONS
	         CHI_LAM(ML,ID)=CHI_PAR(DPTH_INDX,ID)
	         ETA_LAM(ML,ID)=ETA_PAR(DPTH_INDX,ID)
	       END DO
	    ELSE
	       LAM_VEC(ML)=LOG10(ANG_TO_HZ/FREQ)
	       DO ID=1,NUM_IONS
	         IF(LOC_ION_ID .EQ. ION_ID(ID))THEN
	           CHI_LAM(ML,:)=CHI_PAR(:,ID)
	           ETA_LAM(ML,:)=ETA_PAR(:,ID)
	         END IF
	      END DO
	    END IF
	  END DO
!
	  IF(LOC_ION_ID .EQ. ' ')THEN
	    WRITE(6,'(2X,A3,2X,A7,A10,A10)')'ID','  ION  ','  % CHI   ','  % ETA'
	    DO ID=1,NUM_IONS
	      IF(CHI_TOT_LAM(ID) .NE. 0)THEN
	         ISPEC=SPECIES_LNK(ID)
	         WRITE(6,'(2X,I3,2X,A7,2F10.3)')
	1          ID,TRIM(ION_ID(ID)),100.0D0*CHI_TOT_LAM(ID),100.0D0*LAM_VEC(ID)
	      END IF
	    END DO
	    DEALLOCATE (CHI_TOT_LAM)
	    DEALLOCATE (LAM_VEC)
	  ELSE IF(DO_KAP)THEN
	    DO ID=1,NUM_IONS
	      CHI_LAM(:,ID)=1.0E-10_LDP*CHI_LAM(:,ID)/MASS_DENSITY(DPTH_INDX)
	    END DO
	    DO ISPEC=1,NSPEC
	      YV(1:NLAM)=RZERO
	      ELEC=.FALSE.
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        YV(1:NLAM)=YV(1:NLAM)+CHI_LAM(1:NLAM,ID)
	        ELEC=.TRUE.
	      END DO
	      IF(ELEC)THEN
	        CALL DP_CURVE(NLAM,LAM_VEC,YV)
	        WRITE(6,*)ISPEC,SPECIES(ISPEC)
	      END IF
	    END DO
	  ELSE
	    DO I=1,ND
              XMAPV(I)=I
	    END DO
	    CALL MAP_PLOT(CHI_LAM,ETA_LAM,LAM_VEC,XMAPV,NLAM,ND,
	1               " "," "," "," ")
!
	    DEALLOCATE (CHI_LAM)
	    DEALLOCATE (ETA_LAM)
	    DEALLOCATE (LAM_VEC)
	  END IF
!
!
! This section creates an image of the fractional contribution to the opacity
!              as a function of species and wavelength.
!
	ELSE IF(XOPT .EQ. 'LCHI' .OR. XOPT .EQ. 'PPHOT')THEN
	  CALL USR_OPTION(LAM_ST,'LAM_ST','50.0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAM_END','10000.0',FREQ_INPUT)
	  CALL USR_OPTION(K,'DEPTH','20.0','Depth for opacities')
	  CALL USR_OPTION(NLAM,'NLAM','1000.0','# of wavlengths')
	  CALL USR_OPTION(ELEC,'ELEC','F','Include electron scattering?')
!
	  IF(KEV_INPUT)THEN
	    NU_ST=LAM_ST*KEV_TO_HZ
	    NU_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    NU_EN=ANG_TO_HZ/LAM_ST
	    NU_ST=ANG_TO_HZ/LAM_EN
	  ELSE
	   NU_ST=LAM_ST
	   NU_EN=LAM_EN
	  END IF
!
	  LAST_NON_ZERO=.FALSE.
	  NION=0
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      NION=NION+1
	      LAST_NON_ZERO=.TRUE.
	    ELSE IF(LAST_NON_ZERO)THEN
	      NION=NION+1
	      LAST_NON_ZERO=.FALSE.
	    ELSE
	    END IF
	  END DO
!
	  ALLOCATE (CHI_LAM(NLAM,NION))
	  ALLOCATE (ETA_LAM(NLAM,NION))
	  ALLOCATE (CHI_TOT_LAM(NLAM))
	  ALLOCATE (ETA_TOT_LAM(NLAM))
	  ALLOCATE (LAM_VEC(NLAM))
	  DEL_NU=EXP(LOG(NU_EN/NU_ST)/(NLAM-1))
!
	  ETA_LAM(:,:)=RZERO
	  CHI_LAM(:,:)=RZERO
	  CHI_TOT_LAM(:)=RZERO
	  ETA_TOT_LAM(:)=RZERO
	  WRITE(6,'(A3)',ADVANCE='NO')'ML='
	  DO ML=1,NLAM
	    IF(MOD(ML,200) .EQ. 0)THEN
	        WRITE(6,*)' '
	        WRITE(6,'(A3)',ADVANCE='NO')'ML='
	    END IF
	    IF(ML .EQ. NLAM)WRITE(6,'(A)')' '
	    IF(MOD(ML,20) .EQ. 0)WRITE(6,'(1X,I5)',ADVANCE='NO')ML
	    FREQ=NU_ST*(DEL_NU**(ML-1))
	    LAM_VEC(ML)=LOG10(ANG_TO_HZ/FREQ)
	    INCLUDE 'PAR_OPACITIES.INC'
	    J=0; CNT=0
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        J=J+1
                CHI_LAM(ML,J)=CHI_PAR(K,ID)
                ETA_LAM(ML,J)=ETA_PAR(K,ID)
	        LAST_NON_ZERO=.TRUE.
	        CHI_TOT_LAM(ML)=CHI_TOT_LAM(ML)+CHI_LAM(ML,J)
	        ETA_TOT_LAM(ML)=ETA_TOT_LAM(ML)+ETA_LAM(ML,J)
	        IF(XOPT .NE. 'PPHOT' .AND. ML .EQ. NLAM)THEN
	          CNT=CNT+1
	          IF(MOD(CNT-1,5) .EQ. 0)THEN
	            WRITE(6,'(4X,A5,1X,I3)')TRIM(ION_ID(ID)),J
	          ELSE
	            WRITE(6,'(4X,A5,1X,I3)',ADVANCE='NO')TRIM(ION_ID(ID)),J
	          END IF
	        END IF
	      ELSE IF(LAST_NON_ZERO)THEN
	        J=J+1
                CHI_LAM(ML,J)=RZERO
                ETA_LAM(ML,J)=RZERO
	        LAST_NON_ZERO=.FALSE.
	      ELSE
	      END IF
	    END DO
	  END DO
	  IF(XOPT .NE. 'PPHOT')WRITE(6,'(A,/)')' '
	  NION=J
!
	  CALL ESOPAC(ESEC,ED,ND)
	  IF(ELEC)THEN
	    CHI_TOT_LAM(1:NLAM)=CHI_TOT_LAM(1:NLAM)+ESEC(K)
	  END IF
!
	  IF(XOPT .EQ. 'PPHOT')THEN
	    WRITE(6,*)'Total opacity is plotted in red'
	    XV(1:NLAM)=10**(LAM_VEC(1:NLAM))
	    YV(1:NLAM)=CHI_TOT_LAM(1:NLAM)/ESEC(K)
	    CALL DP_CURVE(NLAM,XV,YV)
	    YAXIS='Opacity/(e.s. opacity)'
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES .AND. XSPEC .EQ. UC(ION_ID(ID)) )THEN
	        YV(1:NLAM)=CHI_LAM(1:NLAM,ID)
	        CALL DP_CURVE(NLAM,XV,YV)
	        EXIT
	      END IF
	    END DO
	    DO ISPEC=1,NSPEC
	      IF(XSPEC .EQ. UC(SPECIES(ISPEC)) )THEN
	        DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	          YV(1:NLAM)=CHI_LAM(1:NLAM,ID)
	          CALL DP_CURVE(NLAM,XV,YV)
	        END DO
	      END IF
	    END DO
	  ELSE
	    DO ID=1,NION
	      DO ML=1,NLAM
	        CHI_LAM(ML,ID)=CHI_LAM(ML,ID)/CHI_TOT_LAM(ML)
	        ETA_LAM(ML,ID)=ETA_LAM(ML,ID)/ETA_TOT_LAM(ML)
	      END DO
	    END DO
!
	    DO ID=1,NION
              YMAPV(ID)=ID
	    END DO
!
	    CALL MAP_PLOT(CHI_LAM,ETA_LAM,LAM_VEC,YMAPV,NLAM,NION,
	1               " "," "," "," ")
	  END IF
!
	  DEALLOCATE (CHI_LAM)
	  DEALLOCATE (ETA_LAM)
	  DEALLOCATE (LAM_VEC)
	  DEALLOCATE (CHI_TOT_LAM)
	  DEALLOCATE (ETA_TOT_LAM)
! 
!
! Graphical options.
!
	ELSE IF(XOPT .EQ. 'GR') THEN
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,' ')
	  XAXIS=XAXSAV
!
	ELSE IF(XOPT .EQ. 'GRNL') THEN
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  XAXIS=XAXSAV
!
	ELSE IF(XOPT .EQ. 'GRED')THEN
!
! GRED labels XV by inverting XV(ED).  A repeated electron density has
! two depths, so the mapping is not defined unless LOG10(ED) is
! strictly increasing with depth.
!
	  DO I=1,ND-1
	    IF(ED(I) .LE. 0.0_LDP .OR. ED(I+1) .LE. ED(I))THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: GRED requires positive,'
	      WRITE(T_OUT,*)'strictly increasing electron density'
	      WRITE(T_OUT,*)'ED(I), ED(I+1)=',ED(I),ED(I+1)
	      STOP
	    END IF
	  END DO
	  IF(ED(ND) .LE. 0.0_LDP)THEN
	    WRITE(T_OUT,*)'Error in MAINGEN: GRED requires positive ED'
	    STOP
	  END IF
!
	  TA(1:ND)=XV(1:ND)
!
	  I=LOG10(ED(1));J=LOG10(ED(ND))
	  IF(I .LT. ED(1))I=I+1
	  NXED=(J-I)*2+1
	  NXED=MIN(NXED,NSC)
	  SCED(1)=I
	  DO I=2,NXED
	    SCED(I)=SCED(1)+(I-1)*0.5_LDP
	  END DO
	  DO I=1,ND
	    TB(I)=LOG10(ED(I))
	  END DO
	  CALL MON_INTERP(XED,NXED,IONE,SCED,NXED,TA,ND,TB,ND)
	  TOPLABEL='Log N\de\u(cm\u-3\d)'
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,'TOPLAB')
	  XAXIS=XAXSAV
!
	ELSE IF(XOPT .EQ. 'GRV')THEN
!
! GRV labels XV by inverting XV(V).  A repeated V has two depths, so the
! mapping is not defined unless V is strictly decreasing with depth.
!
	  DO I=1,ND-1
	    IF(V(I+1) .GE. V(I))THEN
	      WRITE(T_OUT,*)'Error in MAINGEN: GRV requires strictly decreasing V'
	      WRITE(T_OUT,*)'Repeated or non-monotonic velocity at depth',I
	      WRITE(T_OUT,*)'V(I), V(I+1)=',V(I),V(I+1)
	      STOP
	    END IF
	  END DO
!
	  TA(1:ND)=XV(1:ND)
!
	  SCED(1)=100*INT(V(1)/100)
	  SCED(2)=2000
	  SCED(3)=1500
	  SCED(4)=1250
	  SCED(5)=1000
	  SCED(6)=750
	  SCED(7)=500
	  SCED(8)=300
	  SCED(9)=100
	  SCED(10)=20
	  SCED(11)=10
	  SCED(12)=1
	  NXED=12
	  TB(1:ND)=V(1:ND)
	  CALL MON_INTERP(XED,NXED,IONE,SCED,NXED,TA,ND,TB,ND)
	  TOPLABEL='V(km/s)'
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,'TOPLAB/ALL')
	  XAXIS=XAXSAV
!
	ELSE IF(XOPT .EQ. 'GROSS')THEN
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL_V2(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM,L_FALSE)
	  WRITE(T_OUT,*)'Rossland optical depth is : ',TAUROSS(ND)
!
	  TB(1:ND)=LOG10(TAUROSS(1:ND))
	  TA(1:ND)=XV(1:ND)
!
	  I=TB(1); J=TB(ND)
	  IF(I .LT. TB(1))I=I+1
	  IF(J .GT. TB(ND))J=J-1
	  SCED(1)=I
	  NXED= (J-I)*2+1
	  NXED=MIN(NXED,NSC)
	  DO I=2,NXED
	    SCED(I)=SCED(I-1)+0.5_LDP
	  END DO
	  CALL MON_INTERP(XED,NXED,IONE,SCED,NXED,TA,ND,TB,ND)
	  TOPLABEL='Log(\gt\dRoss\u)'
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,'TOPLAB')
	  XAXIS=XAXSAV
!
! 
	ELSE IF(XOPT .EQ. 'SM')THEN
	  FLAG=.FALSE.
	  CALL USR_OPTION(T1,'VAL','3.0','Fractional change allowed.')
!
	  DEFAULT=WR_STRING(ND/2)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Smooth from DEPTH to 1')
	  CALL USR_HIDDEN(TYPE,'TYPE','POP','Smooth in POP (def) or DC')
!
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL SMOOTH_V2(ATM(ID)%XzV_F,ATM(ID)%XzVLTE_F,T1,
	1         ATM(ID)%NXzV_F,TYPE,K,
	1         ND,XSPEC,ION_ID(ID),FLAG)
	      CALL SMOOTH_V2(ATM(ID)%DXzV_F,POPION,T1,1,TYPE,K,
	1         ND,XSPEC,ION_ID(ID),FLAG)
	      FOUND=.TRUE.
	    END IF
	  END DO
!
	  IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)' Invalid population type or species unavailable.'
	    GOTO 1
	  END IF
	  IF(.NOT. FLAG)THEN
	    WRITE(T_OUT,*)'Warning - invalid SM extension.'
	  END IF
!
! 
!
! Allow the OCCUPATION PROBAILITIES to be plotted.
!
	ELSE IF(XOPT .EQ. 'WLD')THEN
	  CALL USR_OPTION(LEV,ITEN,IONE,'LEVS','1,2,3,4,5,6,7,8,9,10',
	1      'Levels to plot (10 max)')
	  J=1
	  FOUND=.FALSE.
	  DO WHILE(LEV(J) .GT. 0 .AND. LEV(J) .LE. 200)
	    DO ID=1,NUM_IONS
	      IF(XSPEC .EQ. ION_ID(ID))THEN
	        DO I=1,ND
	          YV(I)=ATM(ID)%W_XzV_F(LEV(J),I)
	        END DO
	        FOUND=.TRUE.
	      END IF
	    END DO
	    IF(.NOT. FOUND)THEN
	      WRITE(T_OUT,*)'Invalid species descritor ', X
	    END IF
	    J=J+1
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='WLD'
	  END DO
!
! 
!
! ***********************************************************************
! ***********************************************************************
!
! DC  --- Plot of Log(departure coefficient).
!           IF(LIN_DC is TRUE, a linear plot is preseneted)
!
! POP --- LOG plot of the population of an individual level.
!
! RAT --- Plot of the ionization fraction of a species. Fraction can be
!           relative to the ion density, or the total species density.
!
! ***********************************************************************
! ***********************************************************************
!
	ELSE IF(XOPT .EQ.'DC' .OR. XOPT .EQ. 'DCRGS' .OR. XOPT .EQ. 'POP'
	1             .OR. XOPT .EQ. 'RAT' .OR. XOPT .EQ. 'TX') THEN
	  DO I=1,10
	    LEV(I)=0
	  END DO
	  IF(XOPT(1:2) .EQ. 'DC')CALL USR_HIDDEN(LIN_DC,'LIN','F','Linear D.C. plots?')
	  CALL USR_HIDDEN(SPEC_FRAC,'SPEC_FRAC','F','Species fraction?')
	  IF(XOPT .EQ. 'RAT')THEN
	   HAM=.TRUE.
	  ELSE
	   HAM=.FALSE.
	  END IF
	  DEFAULT=WR_STRING(HAM)
	  CALL USR_HIDDEN(HAM,'SUM',DEFAULT,'Sum all ion levels?')
	  IF(HAM)THEN
	    LEV(1)=1000
	    LEV(2)=0
	  ELSE
	    CALL USR_OPTION(LEV,ITEN,IONE,'LEVS','1,2,3,4,5,6,7,8,9,10',
	1	 'Levels to plot (10 max)')
	  END IF
!
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      DO I=1,10
	        IF(LEV(I) .EQ. 0 .OR. LEV(I) .GT. 20000)GOTO 1
	          FLAG=.FALSE.
	          CALL SET_DC_OR_POP_OR_TX_V2(YV,LEV(I),ATM(ID)%XzV_F,ATM(ID)%LOG_XzVLTE_F,
	1            ATM(ID)%EDGEXzV_F,ATM(ID)%NXzV_F,ND,T,X,UC(ION_ID(ID)),FLAG)
	          IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
	            T1=RZERO
	            CALL SET_DC_OR_POP_OR_TX(YV,LEV(I),ATM(ID)%DXzV_F,ATM(ID)%DXzV_F,
	1                   T1,IONE,ND,T,X,'D'//TRIM(UC(ION_ID(ID))),FLAG)
	          END IF
!
	        IF(XOPT .EQ. 'DCRGS')THEN
	          YAXIS='Log(b/b\d1\u)'
	          IF(LIN_DC)THEN		!Saves altering SETDC_OR_POP
	            DO J=1,ND
	              YV(J)=10.0_LDP**(YV(J))
	            END DO
	            YAXIS='b/b\d1\u'
	          END IF
	        ELSE IF(XOPT .EQ. 'DC')THEN
	          YAXIS='Log(b)'
	          IF(LIN_DC)THEN		!Saves altering SETDC_OR_POP
	            DO J=1,ND
	              YV(J)=10.0_LDP**(YV(J))
	            END DO
	            YAXIS='b'
	          END IF
	        ELSE IF(XOPT .EQ. 'TX')THEN
	          YAXIS='T\dX\u(10\u4\ \dK)'
	        ELSE IF(XOPT .EQ. 'RAT')THEN
	          IF(SPEC_FRAC)THEN
	          ELSE
	            DO J=1,ND
	              YV(J)=YV(J)-LOG10(POP_ATOM(J))
	            END DO
	            YAXIS='Log(N\dx\u/N\di\u)'
	          END IF
	        ELSE
	          YAXIS='Log(N)'
	          IF(SPEC_FRAC)THEN
	            ISPEC=SPECIES_LNK(ID)
	            DO J=1,ND
	              YV(J)=YV(J)-LOG10(POPDUM(J,ISPEC))
	            END DO
	            YAXIS='Log(N/Nsp)'
	          END IF
	        END IF
!
	        WRITE(DEFAULT,'(I4)')LEV(I); DEFAULT=ADJUSTL(DEFAULT)
	        CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
	      END DO
	    END IF
	  END DO
!
! Plots departure coeeficient for each level at each depth.
!
	ELSE IF(XOPT .EQ. 'DCS')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)BLUE_PEN//'Use B option in PLT_SPEC for plotting'//DEF_PEN
	  WRITE(6,*)' '
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF( (XSPEC .EQ. UC(ION_ID(ID)) .OR. XSPEC .EQ. 'ALL') .AND. ATM(ID)%XzV_PRES)THEN
	      K=0
	      DO L=1,ND
	        DO J=1,ATM(ID)%NXzV_F
	          K=K+1
	          ZV(K)=XV(L)
	          YV(K)=LOG10(ATM(ID)%XzV_F(J,L)/ATM(ID)%XzVLTE_F(J,L))
	        END DO
	      END DO
	      CALL DP_CURVE(K,ZV,YV)
	      FOUND=.TRUE.
	      YAXIS='Log b'
	    END IF
	  END DO
!
	ELSE IF(XOPT .EQ.'TPOP' .OR. XOPT .EQ. 'LPOP' .OR. XOPT .EQ. 'LDC')THEN
!
	  CALL USR_OPTION(K,'DEPTH','20','Depth at which quantities are to be plotted.')
	  FOUND=.FALSE.
	  XAXIS='Level'
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      L=ATM(ID)%NXzV_F
	      IF(XOPT .EQ. 'LDC')THEN
	        DO I=1,ATM(ID)%NXzV_F
	          YV(I)=LOG10( ATM(ID)%XzV_F(I,K) / ATM(ID)%XzVLTE_F(I,K) )
	          ZV(I)=I
	        END DO
	        YAXIS='log b'
	      ELSE
	        DO I=1,ATM(ID)%NXzV_F
	          YV(I)=LOG10(ATM(ID)%XzV_F(I,K))
	          ZV(I)=I
	        END DO
	        YAXIS='log N'
	      END IF
	      FOUND=.TRUE.
	    END IF
	  END DO
	  IF(XOPT .EQ. 'TPOP' .AND. FOUND)THEN
	    WRITE(6,*)'Wavelength set to 5000 Angstroms'
	    FL=ANG_TO_HZ/5000.0_LDP
	    T1=LOG10( 3.0E-10_LDP*OPLIN*R(K)/V(K)/(RONE+SIGMA(K))/FL )
	    DO I=1,L
	      YV(I)=YV(I)+T1
	    END DO
	    YAXIS='Tau/f'
	  END IF
	  IF(FOUND)THEN
	    CALL DP_CURVE(L,ZV,YV)
	  ELSE
	    WRITE(6,*)'Species not found'
	  END IF
!
! 
!
! Compute the ionization ratio for the total population of one ionization
! state to the total population of the next ionization state.
!
	ELSE IF(XOPT .EQ.'ION') THEN
!
	  FOUND=.FALSE.
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. ION_ID(ID))THEN
	      IF(ATM(ID+1)%XzV_PRES)THEN
	        CALL COMPION(ATM(ID)%XzV_F,ATM(ID)%NXzV_F,
	1                ATM(ID+1)%XzV_F,ATM(ID+1)%NXzV_F,YV,ND)
	      ELSE
	        CALL COMPION(ATM(ID)%XzV_F,ATM(ID)%NXzV_F,ATM(ID)%DXzV_F,IONE,YV,ND)
	      END IF
	      FOUND=.TRUE.
	    END IF
	  END DO
	  IF(.NOT. FOUND)THEN
	    WRITE(T_OUT,*)' Invalid species'
	    GOTO 1
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='Log(X\u(n-1)+\d/X\un+\d)'
! 
! ^L
!
! Computes the column density for 2 successive ionization stages, and for
! the atom. If XV is the column density for species XzV, the program returns
! XV/X and XSIX/X where X is the species column density. Inserted for comparison
! with SETI results.
!
	ELSE IF(XOPT .EQ. 'QF') THEN
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. SPECIES(ISPEC))THEN
	      ZV(ND)=RZERO
	      YAXIS='q'
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO I=1,ND
	            T2=RZERO
	            DO J=1,ATM(ID)%NXzV_F
	              T2=T2+ATM(ID)%XzV_F(J,I)
	            END DO
	            TA(I)=POPDUM(I,ISPEC)
	            TB(I)=T2
	            TC(I)=ATM(ID)%DXzV_F(I)
	          END DO
	          CALL TORSCL(DTAU,TA,R,Z,XM,ND,METHOD,TYPE_ATM)
	          TA(1:ND)=DTAU(1:ND)
	          CALL TORSCL(DTAU,TB,R,Z,XM,ND,METHOD,TYPE_ATM)
	          TB(1:ND)=DTAU(1:ND)
	          CALL TORSCL(DTAU,TC,R,Z,XM,ND,METHOD,TYPE_ATM)
	          TC(1:ND)=DTAU(1:ND)
!
	          DO I=1,ND
	            YV(I)=LOG10(TB(I)/TA(I))
	            ZV(I)=LOG10(TC(I)/TA(I))
	          END DO
	          CALL DP_CURVE(ND,XV,YV)
	        END IF
	      END DO
	      IF(ZV(ND) .NE. 0.0_LDP)CALL DP_CURVE(ND,XV,ZV)
	      YAXIS='CD(XV)/CD(X); CD(XSIX)/CD(X)'
	    END IF
	  END DO
! 
!
! Compute the ionization ratio for the total population of one ionization
! state to the total population of the next ionization state.
!
	ELSE IF(XOPT .EQ. 'IF') THEN
	  CALL USR_OPTION(SPEC_FRAC,'SPEC_FRAC','T','Species fraction?')
	  DO ISPEC=1,NSPEC
	    IF(XSPEC .EQ. SPECIES(ISPEC))THEN
	      IF(SPEC_FRAC)THEN
	        DO I=1,ND
	          TA(I)=POPDUM(I,ISPEC)
	        END DO
	        YAXIS='Log N('//TRIM(SPECIES_ABR(ISPEC))//
	1                       '\un+\d)/N('//TRIM(SPECIES_ABR(ISPEC))//')'
	      ELSE
	        DO I=1,ND
	          TA(I)=POP_ATOM(I)
	        END DO
	        YAXIS='Log N('//TRIM(SPECIES_ABR(ISPEC))//'\un+\d)/N(total)'
	      END IF
	      TITLE=' '
	      CNT=0
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ATM(ID)%XzV_PRES)THEN
	          CNT=CNT+1
	          DO I=1,ND
	            T1=RZERO			!Using T1 avoids underflow
	            DO J=1,ATM(ID)%NXzV_F
	              T1=T1+ATM(ID)%XzV_F(J,I)
	            END DO
	            YV(I)=LOG10(T1/TA(I))
	            ZV(I)=LOG10(ATM(ID)%DXzV_F(I)/TA(I))
	          END DO
	          DEFAULT=TRIM(PLT_ION_ID(ID))
	          J=INDEX(DEFAULT,'k'); IF(J .NE. 0)DEFAULT(J:J)='i'
	          CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
	          WRITE(6,'(A,I2,A,A)')' Curve ',CNT,' is due to: ',TRIM(ION_ID(ID))
	        END IF
	      END DO
	      ID=SPECIES_END_ID(ISPEC)
              DEFAULT=TRIM(PLT_ION_ID(ID))
              J=INDEX(DEFAULT,'k'); IF(J .NE. 0)DEFAULT(J:J)='i'
	      CALL DP_CURVE_LAB(ND,XV,ZV,DEFAULT)
	      CNT=CNT+1
!	      WRITE(6,'(A,I2,A,A)')' Curve ',CNT,' is due to: ',TRIM(ION_ID(ID))
	      WRITE(6,'(A,I2,A,A,F6.2)')' Curve ',CNT,' is due to: ',TRIM(ION_ID(ID)),ATM(ID)%ZXzV
	    END IF
	  END DO
! 
	ELSE IF(XOPT .EQ. 'MODSUM')THEN
!
! This option was installed to revise the incorrect vaules of R(Tau) and
! Tau output by CMFGEN to MOD_SUM when clumping is present. Problem potentially
! applies to all MOD_SUM files with a CMFGEN program date earlier than
! 01-Sep-199.
!
! Compute the Radius and Velocity at Tau=10, and at Tau=2/3.
!
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	  TA(1:ND)=RZERO ; DO I=2,ND ; TA(I) = TA(I-1)+DTAU(I-1) ; END DO
	  TB(1)=2.0_LDP/3.0_LDP  ; TB(2)=10.0_LDP
	  CALL MON_INTERP(TC,ITWO,IONE,TB,ITWO,R,ND,TA,ND)
	  CALL MON_INTERP(AV,ITWO,IONE,TB,ITWO,V,ND,TA,ND)
!
	  OPEN(UNIT=LU_IN,FILE='MOD_SUM',STATUS='OLD',ACTION='READ')
	  OPEN(UNIT=LU_OUT,FILE='REV_MOD_SUM',STATUS='UNKNOWN')
	   READ(LU_IN,'(A)')STRING
	   DO WHILE( INDEX(STRING,'Tau=') .EQ. 0)
	     WRITE(LU_OUT,'(A)')TRIM(STRING)
	     READ(LU_IN,'(A)')STRING
	   END DO
!
! Output summary of Teff, R, and V at RSTAR, Tau=10, and TAU=2/3.
!
	   NEXT_LOC=1  ;   STRING=' '
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TA(ND))
	   T1=1.0E+10_LDP*R(ND)/RAD_SUN(); CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*/Rsun',T1)
	   T1=TEFF_SUN()*(LUM/T1**2)**0.25_LDP
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'T*  ',T1)
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',V(ND))
	   WRITE(LU_OUT,'(A)')TRIM(STRING)
!
	   NEXT_LOC=1  ;   STRING=' '
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(2))		!10.0D0
	   T1=1.0E+10_LDP*TC(2)/RAD_SUN()
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	   T1=TEFF_SUN()*(LUM/T1**2)**0.25_LDP
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(2))
	   WRITE(LU_OUT,'(A)')TRIM(STRING)
!
	   NEXT_LOC=1  ;   STRING=' '
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(1))		!0.67D0
	   T1=1.0E+10_LDP*TC(1)/RAD_SUN()
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	   T1=TEFF_SUN()*(LUM/T1**2)**0.25_LDP
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	   CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(1))
	   WRITE(LU_OUT,'(A)')TRIM(STRING)
!
	   READ(LU_IN,'(A)')STRING		!Tau=10 string
	   READ(LU_IN,'(A)')STRING		!Tau=1 string
!
	   IOS=0
	   READ(LU_IN,'(A)',IOSTAT=IOS)STRING
	   DO WHILE(IOS .EQ. 0)
	     WRITE(LU_OUT,'(A)')TRIM(STRING)
	     READ(LU_IN,'(A)',IOSTAT=IOS)STRING
	   END DO
	   WRITE(LU_OUT,'(A)')'R(Tau) values have been revised.'
	   CLOSE(LU_OUT); CLOSE(LU_IN)
!
	ELSE IF(XOPT .EQ. 'MEANOPAC')THEN
!
! Output Opacities, Optical depth scales and optical depth increments to file
! in an identical format to MEANOPAC created by CMFGEN.
!
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(TA,TCHI,R,R,dCHIdR,ND)
	  TCHI(1:ND)=FLUX_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(TB,TCHI,R,R,dCHIdR,ND)
	  CALL ESOPAC(ESEC,ED,ND)		!Elec. scattering opacity.
	  TCHI(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	  TA(ND)=RZERO; TB(ND)=RZERO; DTAU(ND)=RZERO
!
	  CALL GEN_ASCI_OPEN(LU_OUT,'MEANOPAC','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(LU_OUT,
	1    '( ''     R        I   Tau(Ross)   /\Tau   Rat(Ross)'//
	1    '  Chi(Ross)  Chi(ross)  Chi(Flux)   Chi(es) '//
	1    '  Tau(Flux)  Tau(es)  Rat(Flux)  Rat(es)'' )' )
	    IF(R(1) .GE. 1.0E+05_LDP)THEN
	      FMT='( 1X,1P,E10.4,2X,I3,1X,1P,E9.3,2(2X,E8.2),1X,'//
	1          '4(2X,E9.3),2(2X,E8.2),2(2X,E8.2) )'
	    ELSE
	      FMT='( 1X,F10.4,2X,I3,1X,1P,E9.3,2(2X,E8.2),1X,'//
	1          '4(2X,E9.3),2(2X,E8.2),2(2X,E8.2) )'
	    END IF
	    DO I=1,ND
	      IF(I .EQ. 1)THEN
	        T1=RZERO		!Rosseland optical depth scale
	        T2=RZERO		!Flux optical depth scale
	        T3=RZERO		!Elec. scattering optical depth scale.
	        TC(1:3)=RZERO
	      ELSE
	        T1=T1+TA(I-1)
	        T2=T2+TB(I-1)
	        T3=T3+DTAU(I-1)
	        TC(1)=TA(I)/TA(I-1)
	        TC(2)=TB(I)/TB(I-1)
	        TC(3)=DTAU(I)/DTAU(I-1)
	      END IF
	      WRITE(LU_OUT,FMT)R(I),I,T1,TA(I),TC(1),
	1        ROSS_MEAN(I),ROSS_MEAN(I),FLUX_MEAN(I),ESEC(I),
	1        T2,T3,TC(1),TC(2)
	    END DO
	    WRITE(LU_OUT,'(//,A,/,A,/,A)')
	1       'NB: Mean opacities do not include effect of clumping.',
	1       'NB: Optical depth scale includes effect of clumping.',
	1       'INT_dBdT(opacity) set to ROSS_MEAN.'
	  CLOSE(UNIT=LU_OUT)
!
! Option to determine the velocity and radius of a star at a give
! Rosseland optical depth.
!
	ELSE IF(XOPT .EQ. 'CHKV')THEN
	  CALL USR_HIDDEN(T1,'TAU','100.0',' ')
!
! Determine radius where TAUROSS=T1
!
	  IF(.NOT. ROSS)THEN
	    WRITE(T_OUT,*)'Error - Roesseland optical depth not computed'
	  ELSE
	    DO I=1,ND
	      IF(V(I) .LE. 0.0_LDP)THEN
	        WRITE(T_OUT,*)'Error in MAINGEN: CHKV requires positive V'
	        WRITE(T_OUT,*)'Depth, V=',I,V(I)
	        STOP
	      END IF
	    END DO
	    IF(TAUROSS(ND) .LT. T1)THEN
	      T2=LOG(T1/TAUROSS(ND-2)) / LOG(TAUROSS(ND-1)/TAUROSS(ND-2))
	      NEW_RSTAR=EXP( LOG(R(ND-2)) - T2*LOG(R(ND-2)/R(ND)) )
	      T2=LOG(R(ND-2)/NEW_RSTAR) / LOG(R(ND-2)/R(ND))
	      NEW_VSTAR=EXP( LOG(V(ND-2)) + T2*LOG( V(ND)/V(ND-2) ) )
	    ELSE
	      I=ND-1
	      DO WHILE(TAUROSS(I) .GT. T1)
	        I=I-1
	      END DO
	      T2=LOG(T1/TAUROSS(I)) / LOG(TAUROSS(I+1)/TAUROSS(I))
	      NEW_RSTAR=EXP( LOG(R(I)) - T2*LOG(R(I)/R(I+1)) )
	      T2=LOG(R(I)/NEW_RSTAR) / LOG(R(I)/R(I+1))
	      NEW_VSTAR=EXP(  LOG(V(I)) + T2*LOG( V(I+1)/V(I) )  )
	    END IF
	    WRITE(T_OUT,730)T1,NEW_RSTAR,NEW_VSTAR
730	    FORMAT(1X,'Estimate of R* where Rosseland Optical depth is',
	1     F8.2,' is ',F10.4,/,
	1     1X,'Velocity at R* is ',1PE14.5,'km/s')
	  END IF
!
	ELSE IF(XOPT .EQ. 'COLR')THEN
	  IF(DPTH_INDX .LT. 1 .OR. DPTH_INDX .GT. ND)DPTH_INDX=ND/2
	  DEFAULT=WR_STRING(DPTH_INDX)
	  CALL USR_OPTION(DPTH_INDX,'Depth',DEFAULT,'Input depth to check collision quantities')
!
	  WRITE(6,'(/,A,3X,ES14.4)')' Radius (10^10 cm)', ED(DPTH_INDX)
	  WRITE(6,'(A,3X,ES14.4)')' Electron density', ED(DPTH_INDX)
	  WRITE(6,'(A,3X,ES14.4)')' Velocity (km/s) ', ED(DPTH_INDX)
!
	  STRING=' '
	  WRITE(6,*)' '
	  WRITE(6,*)'Options are:'
	  WRITE(6,*)'   NR: Output net rates;       DR: Downward rates'
	  WRITE(6,*)'   CR: Output cooling terms;   CL: Cooling lines'
	  WRITE(6,*)'   SL: Ouput rates for a single level (or multiple levels in SL)'
	  WRITE(6,*)' '
	  DEFAULT=COL_OPT
	  CALL USR_OPTION(COL_OPT,'TYPE',DEFAULT,'NR, DR, CR, CL,  SL')
	  IF(UC(COL_OPT(1:2)) .EQ. 'NR')THEN
	    STRING='NET_RATES'
	  ELSE IF(UC(COL_OPT(1:2)) .EQ. 'CR')THEN
	    STRING='COOL_RATES'
	  ELSE IF(UC(COL_OPT(1:2)) .EQ. 'CL')THEN
	    STRING='COOL_LINES'
	  ELSE IF(UC(COL_OPT(1:2)) .EQ. 'SL')THEN
	    STRING='SUPER_LEVEL'
	  ELSE
	    STRING=' '
	  END IF
	  I=DPTH_INDX
	  TMP_ED=ED(I)
	  T1=T(I)
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      CALL SUBCOL_MULTI_V4(
	1         OMEGA_F,dln_OMEGA_F_dlnT,
	1         OMEGA_S,dln_OMEGA_S_dlnT,
	1         ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),UNIT_VEC,ATM(ID)%NXzV_F,
	1         ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),ATM(ID)%W_XzV_F(1,I),
	1         ATM(ID)%EDGEXzV_F,ATM(ID)%AXzV_F,ATM(ID)%GXzV_F,
	1         ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,ATM(ID)%ZXzV,
	1         ID,TRIM(ION_ID(ID))//'_COL_DATA',OMEGA_GEN_V3,
	1         ATM(ID)%F_TO_S_XzV,TEMP,T1,TMP_ED,IONE)
	      IF(STRING .EQ. 'COOL_LINES')THEN
	  	CALL USR_OPTION(FLAG,'NEW_FILE',FIRST_COLR,'New file?')
	        FIRST_COLR='F'
	        J=10
	        CALL WR_COL_LINES(OMEGA_S,ATM(ID)%AXzV_F,ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),
	1               ATM(ID)%W_XzV_F(1,I),ATM(ID)%GXzV_F,
	1               ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,
	1               TRIM(ION_ID(ID)),R(I),V(I),SIGMA(I),I,J,LU_COL,FLAG)
	      ELSE IF(STRING .EQ. 'SUPER_LEVEL')THEN
!
	        WRITE(6,'(/,A)')' You may input a single level or multiple levels that are linked into a super level.'
	  	WRITE(6,'(A,/)')' At this stage you mus use the level number as input.'
	        CALL USR_OPTION(LEV,ITEN,IONE,'LEV','1','Levels of interest')
	        CALL WR_COL_SL(OMEGA_S,ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),
	1                ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,
	1                TRIM(XSPEC),LU_COL,LEV,ITEN,DPTH_INDX)
	      ELSE
	        CALL WR_COL_RATES(OMEGA_S,ATM(ID)%XzV_F(1,I),ATM(ID)%XzVLTE_F(1,I),
	1               ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,
	1               TRIM(XSPEC)//'R',LU_COL,STRING)
	      END IF
	    END IF
	  END DO
!
	ELSE IF(XOPT .EQ. 'COLL')THEN
!
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'Options are:'
	  WRITE(6,*)' NSR -- net rate TO level of interest summed of all levels'
	  WRITE(6,*)' NR  -- net rate TO level of interest'
	  WRITE(6,'(A,A,/)')' NI  -- indicidual rates',DEF_PEN
!
	  IF(COL_OPT .EQ. 'CL')COL_OPT='IND_RATES'
	  DEFAULT=COL_OPT
	  CALL USR_OPTION(COL_OPT,'TYPE',DEFAULT,'NR (NET_RATES), IN (IND_RATES), or TR')
	  IF(UC(COL_OPT(1:2)) .EQ. 'NR' .OR. UC(COL_OPT) .EQ. 'NET_RATES')THEN
	    COL_OPT='NET_RATES'
	  ELSE IF(UC(COL_OPT(1:2)) .EQ. 'IN' .OR. UC(COL_OPT) .EQ.  'IND_RATES')THEN
	    COL_OPT='IND_RATES'
	  ELSE IF(UC(COL_OPT(1:2)) .EQ. 'TR' .OR. UC(COL_OPT) .EQ.  'TOT_RATE')THEN
	    COL_OPT='TOTAL_RATE'
	  ELSE
	    WRITE(6,*)'Col option not recognized'
	    GOTO 1
	  END IF
!
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      COL_RATE_FILE=TRIM(ION_ID(ID))//'_'//TRIM(COL_OPT)
	      INQUIRE(FILE=COL_RATE_FILE,EXIST=FLAG)
	      IF(FLAG)THEN
	        CALL USR_OPTION(ELEC,'OVER','F','Overwite existing file?')
	        IF(ELEC)THEN
	          OPEN(UNIT=LU_COL,FILE=TRIM(COL_RATE_FILE),STATUS='OLD',ACTION='WRITE')
	        ELSE
	          OPEN(UNIT=LU_COL,FILE=TRIM(COL_RATE_FILE),STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	        END IF
	      ELSE
	        OPEN(UNIT=LU_COL,FILE=TRIM(COL_RATE_FILE),STATUS='NEW',ACTION='WRITE')
	      END IF
	      EXIT
	    END IF
	  END DO
!
	  CALL USR_OPTION(NL,'NL','1','Level of interest (1 to NF)')
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      IF(LAST_COL_SPECIES .NE. UC(ION_ID(ID)))THEN
	        IF(ALLOCATED(COL))DEALLOCATE(COL,dCOL)
	        J=ATM(ID)%NXzV_F
	        ALLOCATE(COL(J,J,ND),dCOL(J,J,ND))
	        CALL SUBCOL_MULTI_V4(
	1           OMEGA_F,dln_OMEGA_F_dlnT,
	1           COL,dCOL,
	1           ATM(ID)%XzV_F,ATM(ID)%XzVLTE_F,UNIT_VEC,ATM(ID)%NXzV_F,
	1           ATM(ID)%XzV_F,ATM(ID)%XzVLTE_F,ATM(ID)%W_XzV_F,
	1           ATM(ID)%EDGEXzV_F,ATM(ID)%AXzV_F,ATM(ID)%GXzV_F,
	1           ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,ATM(ID)%ZXzV,
	1           ID,TRIM(ION_ID(ID))//'_COL_DATA',OMEGA_GEN_V3,
	1           ATM(ID)%F_TO_S_XzV,TA,T,ED,ND)
	         LAST_COL_SPECIES=UC(ION_ID(ID))
	      END IF
!
	      IF(COL_OPT .EQ. 'TOTAL_RATE')THEN
	        STRING=TRIM(ION_ID(ID))//'('//TRIM(ATM(ID)%XzVLEVNAME_F(NL))//')'
	        TC=RZERO
	        DO K=1,ATM(ID)%NXzV_F
	          TA(1:ND)=COL(K,NL,1:ND)*ATM(ID)%XzV_F(K,1:ND)			!Col rates into NL
	          TB(1:ND)=COL(NL,K,1:ND)*ATM(ID)%XzV_F(NL,1:ND)		!Col rates out of NL
	          TC(1:ND)=TC(1:ND)+(TA(1:ND)-TB(1:ND))
	        END DO
	        CALL WRITE_VEC(TC,ND,STRING,LU_COL)
!
	      ELSE IF(COL_OPT .EQ. 'NET_RATE')THEN
	        CALL  USR_OPTION(L,'Depth','30','Depth to order rates')
	        TA=RZERO; TB=RZERO; K=ATM(ID)%NXzV_F
	        TA(1:K)=COL(:,NL,L)*ATM(ID)%XzV_F(:,L)		!Col rates into NL
	        TA(NL)=COL(NL,NL,L)*ATM(ID)%XzVLTE_F(NL,L)	!Col rates into NL
	        TB(1:K)=COL(NL,:,L)*ATM(ID)%XzV_F(NL,L)		!Col rates out of NL
	        TA(1:K)=TB(1:K)-TA(1:K)
	        DO ISPEC=1,MIN(30,ATM(ID)%NXzV_F)
	          LEV(1:1)=MAXLOC(TA); K=LEV(1)
	          STRING=TRIM(ION_ID(ID))//'('//TRIM(ATM(ID)%XzVLEVNAME_F(K))//
	1                  '-'//TRIM(ATM(ID)%XzVLEVNAME_F(NL))//')'
	          TC(1:ND)=COL(K,NL,1:ND)*ATM(ID)%XzV_F(K,1:ND)	-
	1                  COL(NL,K,1:ND)*ATM(ID)%XzV_F(NL,1:ND)		!Col rates out of NL
	          CALL WRITE_VEC(TC,ND,STRING,LU_COL)
	          TA(K)=RZERO
	        END DO
!
	      ELSE
!
	        CALL  USR_OPTION(L,'Depth','30','Depth to order rates')
	        TA=RZERO; TB=RZERO; K=ATM(ID)%NXzV_F
	        TA(1:K)=COL(:,NL,L)*ATM(ID)%XzV_F(:,L)		!Col rates into NL
	        TA(NL)=COL(NL,NL,L)*ATM(ID)%XzVLTE_F(NL,L)	!Col rates into NL
	        TB(1:K)=COL(NL,:,L)*ATM(ID)%XzV_F(NL,L)		!Col rates out of NL
!
	        DO ISPEC=1,MIN(30,ATM(ID)%NXzV_F)
	          LEV(1:1)=MAXLOC(TA); K=LEV(1)
	          IF(K .EQ. NL)THEN
	            STRING=TRIM(ION_ID(ID))//'(ION-'//TRIM(ATM(ID)%XzVLEVNAME_F(NL))//')'
	            TC(1:ND)=COL(NL,NL,1:ND)*ATM(ID)%XzVLTE_F(NL,1:ND)		!Col rates into NL
	          ELSE
	            STRING=TRIM(ION_ID(ID))//'('//TRIM(ATM(ID)%XzVLEVNAME_F(K))//
	1                  '-'//TRIM(ATM(ID)%XzVLEVNAME_F(NL))//')'
	            TC(1:ND)=COL(K,NL,1:ND)*ATM(ID)%XzV_F(K,1:ND)		!Col rates into NL
	          END IF
	          CALL WRITE_VEC(TC,ND,STRING,LU_COL)
	          TA(K)=RZERO
	        END DO
!
	        DO ISPEC=1,MIN(30,ATM(ID)%NXzV_F)
	          LEV(1:1)=MAXLOC(TB); K=LEV(1)
	          IF(K .EQ. NL)THEN
	            STRING=TRIM(ION_ID(ID))//'('//TRIM(ATM(ID)%XzVLEVNAME_F(NL))//'-ION)'
	          ELSE
	            STRING=TRIM(ION_ID(ID))//'('//TRIM(ATM(ID)%XzVLEVNAME_F(NL))//
	1                  '-'//TRIM(ATM(ID)%XzVLEVNAME_F(K))//')'
	    	  END IF
	          TC(1:ND)=COL(NL,K,1:ND)*ATM(ID)%XzV_F(NL,1:ND)
	          CALL WRITE_VEC(TC,ND,STRING,LU_COL)
	          TB(K)=RZERO
	        END DO
	      END IF
	      CLOSE(LU_COL)
	      EXIT
	    END IF
	  END DO
!
!	ELSE IF(XOPT .EQ. 'NONT')THEN
!          DO ID=1,NUM_IONS
!            IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
!	       CALL DISP_BETHE_APPROX_V5(Q,NL,NUP,XKT,dXKT_ON_XKT,NKT,ID,DPTH_INDX)
!	    END IF
!	  END DO
!
	ELSE IF(XOPT .EQ. 'CSUM')THEN
	  TMP_ED=RONE
	  CALL USR_OPTION(TEMP,'T','1.0','Input T')
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)) .OR. (XSPEC .EQ. 'ALL' .AND.  ATM(ID)%XzV_PRES))THEN
	      CALL GET_COL_SUMMARY_V1(OMEGA_F,
	1         ATM(ID)%EDGEXzV_F,ATM(ID)%AXzV_F,ATM(ID)%GXzV_F,
	1         ATM(ID)%XzVLEVNAME_F,ATM(ID)%ZXzV,ATM(ID)%NXzV_F,
	1         TEMP,TRIM(ION_ID(ID))//'_COL_DATA')
	     END IF
	   END DO
!
	ELSE IF(XOPT .EQ. 'COL' .OR. XOPT .EQ. 'CRIT')THEN
	  TMP_ED=RONE
	  CALL USR_OPTION(T1,'T','1.0','Input T')
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      CALL SUBCOL_MULTI_V4(
	1         OMEGA_F,dln_OMEGA_F_dlnT,
	1         OMEGA_S,dln_OMEGA_S_dlnT,
	1         UNIT_VEC,UNIT_VEC,ZERO_VEC,ATM(ID)%NXzV_F,
	1         UNIT_VEC,UNIT_VEC,UNIT_VEC,
	1         ATM(ID)%EDGEXzV_F,ATM(ID)%AXzV_F,ATM(ID)%GXzV_F,
	1         ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,ATM(ID)%ZXzV,
	1         ID,TRIM(ION_ID(ID))//'_COL_DATA',OMEGA_GEN_V3,
	1         ATM(ID)%F_TO_S_XzV,TA,T,TMP_ED,IONE)
	          IF(XOPT .EQ. 'COL')THEN
	            CALL WR_COL(OMEGA_F,ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,XSPEC,LU_COL,' ')
	          ELSE
	            CALL WR_CRIT(OMEGA_F,T1,HDKT,ATM(ID)%AXzV_F,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,
	1                    ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,XSPEC,LU_COL)
	         END IF
	    END IF
	  END DO
!
! 
!
	ELSE IF(XOPT .EQ. 'RNT')THEN
	  DEFAULT='NON_THERM_SPEC_INFO'
	  CALL USR_OPTION(FILENAME,'FILE',DEFAULT,'File name with fractional energies')
	  CALL RD_NON_THERM_SPEC(YV,ZV,WV,ND,FILENAME,LU_IN,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error unable to open file. IOS=',IOS
	    GOTO 1
	  END IF
!
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')'Curve ordering is Felec, Fion, Fexec '
	  CALL DP_CURVE(ND,XV,YV)
	  CALL DP_CURVE(ND,XV,ZV)
	  CALL DP_CURVE(ND,XV,WV)
	  YAXIS='Channel fraction'
!
	ELSE IF(XOPT .EQ. 'PHOT')THEN
	  CALL USR_OPTION(TEMP,'Nu','0.0','Input frequency (10^15)Hz')
	  CALL USR_OPTION(PHOT_ID,'Nu','1','Photoionization route')
	  FLAG=.FALSE.				!Don't return edge value.
	  IF(TEMP .EQ. 0.0_LDP)FLAG=.TRUE.
	  I=1
!
	  K=0
	  OMEGA_F(:,:)=RZERO
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      K=ATM(ID)%NXzV_F
	      CALL SUB_PHOT_GEN(ID,OMEGA_F,TEMP,ATM(ID)%EDGEXzV_F,
	1            ATM(ID)% NXzV_F,PHOT_ID,FLAG)
	      EXIT
	    END IF
	 END DO
	 IF(K .EQ. 0)THEN
	    WRITE(T_OUT,*)'Species not recognized'
	    GOTO 1
	  END IF
!
	  WRITE(T_OUT,*)'Cross sections (Mbarns) for frequency ',TEMP
	  DO I=1,K
	    WRITE(T_OUT,'(1X,I4,3X,1P,E12.5)')I,OMEGA_F(I,1)/1.0D-08
	  END DO
	  WRITE(LU_CROSS,*)'Cross sections (Mbarns) for frequency ',TEMP
	  DO I=1,K
	    WRITE(LU_CROSS,'(1X,I4,3X,1P,E12.5)')I,OMEGA_F(I,1)/1.0D-08
	  END DO
! 
!
	ELSE IF(XOPT .EQ. 'PLTPHOT')THEN
	  CALL PLTPHOT_SUB(XSPEC,XV,YV,WV,ZV,N_PLT_MAX,OMEGA_F,N_MAX,XAXIS,YAXIS,XRAYS,VSM_DIE_KMS,ND)
!
	ELSE IF(XOPT .EQ. 'PLTARN')THEN
	  CALL PLT_ARN(XSPEC,ND,XV,YV,N_PLT_MAX)
	  XAXIS='eV'; YAXIS='\gs(cm\u-2\d)'	
!
	ELSE IF(XOPT .EQ. 'RDDIE')THEN
	  CALL USR_OPTION(DIE_REG,'REG','F','Include dielectronic lines')
	  CALL USR_OPTION(DIE_WI,'WI','F','Include WI dielectronic lines')
	  CALL USR_OPTION(VSM_DIE_KMS,'VSM','3000.0D0',
	1                                   'Smoothing width in km/s')
	  DO ID=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	      IF(DIE_REG .OR. DIE_WI)THEN
	        CALL RD_PHOT_DIE_V1(ID,ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,
	1             ATM(ID)%NXzV_F,ATM(ID)%GIONXzV_F,
	1             VSM_DIE_KMS,DIE_REG,DIE_WI,
	1             ION_ID(ID),LU_IN,LU_LOG,'DIE'//TRIM(ION_ID(ID)))
	      END IF
	    END IF
	  END DO
!
c 
!
! Section to compute the recombination rates for various species.
! Default is for ionizatio to the ground state. The photionzation route
! is a hidden parameters. Both GION and EXC_EN must be specified. These
! are input to accomodate TERM spitting.
!
	ELSE IF(XOPT .EQ. 'NRR')THEN
	  CALL USR_OPTION(T1,'T','1.0','Input T (in 10^4 K)')
	  CALL USR_HIDDEN(PHOT_ID,'PHOT_ID','1',
	1      'Photoionization route')
	  EXC_EN=RZERO
	  IF(PHOT_ID .NE. 1)THEN
	    CALL USR_OPTION(EXC_EN,'EXC_EN',' ',
	1     'Excitation Energy (cm^-1) of final state')
	  END IF
	  EXC_EN=1.0E-15_LDP*C_CMS*EXC_EN
	  CALL USR_OPTION(TMP_GION,'GION',' ',
	1      'G for ION (No def)')
	  J=0		!Set in case no identification.
	  I=INDEX(X,' ')
!
	  DO ID=1,NUM_IONS
	    I=ID
	    IF(XSPEC .EQ. UC(ION_ID(ID)) .AND. ATM(ID)%XzV_PRES)THEN
	      CALL RECOM_CHK_V2(TA,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,TMP_GION,
	1            ATM(ID)% NXzV_F,EXC_EN,PHOT_ID,SUB_PHOT_GEN,I,T1)
	      STRING=ION_ID(ID)
              J=ATM(ID)%NXzV_F
	      EXIT
	    END IF
	  END DO
!
	  IF(J .NE. 0)THEN
	    WRITE(LU_REC,*)' '
	    WRITE(LU_REC,*)'Temperatre (10^4 K) is',T1
	    WRITE(T_OUT,*)'Recombination rates for ',TRIM(STRING)
	    WRITE(LU_REC,*)'Recombination rates for ',TRIM(STRING)
	    WRITE(T_OUT,'(1X,A,1X,I3)')'Photoionization route number is',PHOT_ID
	    WRITE(LU_REC,'(1X,A,1X,I3)')'Photoionization route number is',PHOT_ID
	    WRITE(T_OUT,*)' '
	    WRITE(LU_REC,*)' '
	    T1=RZERO
	    DO I=1,J
	      WRITE(T_OUT, '(I4,ES12.4,2X,A)')I,TA(I),TRIM(ATM(ID)%XzVLEVNAME_F(I))
	      WRITE(LU_REC,'(I4,ES12.4,2X,A)')I,TA(I),TRIM(ATM(ID)%XzVLEVNAME_F(I))
	      T1=T1+TA(I)
	    END DO
	    WRITE(T_OUT,'(1X,A,1X,1PE12.4)')'Total recombination rate is:',T1
	    WRITE(T_OUT,'(1X,A,1X,1PE10.2)')'Value for GION used was:',TMP_GION
	    WRITE(LU_REC,'(1X,A,1X,1PE12.4)')'Total recombination rate is:',T1
	    WRITE(LU_REC,'(1X,A,1X,1PE10.2)')'Value for GION used was:',TMP_GION
	  END IF
!
! 
!
	ELSE IF(XOPT .EQ. 'GNTFF')THEN
	    CALL USR_OPTION(CHI,3,3,'GNT',' ','LAM(um) , T & Z')
	    CHI(1)=1.0E-04_LDP*ANG_TO_HZ/CHI(1)
	    CHI(4)=GFF(CHI(1),CHI(2),CHI(3))
	    WRITE(T_OUT,2000)CHI(4)
2000	    FORMAT(3X,'The Free-free gaunt factor is',1X,1PE10.3)
!
	ELSE IF(XOPT .EQ. 'GNTBF')THEN
	  CALL USR_OPTION(CHI,3,3,'GNT',' ','Lam(um) ,Level & Z')
	    CHI(1)=1.0E-04_LDP*ANG_TO_HZ/CHI(1)
	    I=(CHI(2)+0.000002_LDP)
	    CHI(4)=GBF(CHI(1),I,CHI(3))
	    WRITE(T_OUT,2001)CHI(4)
2001	    FORMAT(3X,'The Bound-free gaunt factor is',1X,1PE10.3)
!
	ELSE IF(XOPT .EQ. 'LAM')THEN
!
! NB. FL and FREQ have been equivalenced.
!
	  DEFAULT=WR_STRING(FREQ)
	  CALL USR_HIDDEN(T1,'FREQ',DEFAULT,' ')
	  T2=LAMVACAIR(T1)		!Wavelength(Angstroms)
	  WRITE(T_OUT,'(1X,''Lambda(in air)='',1PE14.6)')T2
	  WRITE(T_OUT,'(1X,''Lambda(in vac)='',1PE14.6)')ANG_TO_HZ/T1
	ELSE IF(XOPT .EQ. 'TIT')THEN
	  CALL USR_OPTION(NAME,'Title',' ','Title for all graphs')
	ELSE IF(XOPT .EQ. 'METHOD')THEN
	  CALL USR_OPTION(METHOD,'LOGLOG',' ','Tau option (LOGLOG,LOGMON, ZERO')
!
! 
!
! SETREC sets TA to DCI if DCI is non-zero, otherwise TA is set to C2(1, ).
! Thus can compute C2-CI recombination rate even if CI is not present.
!
	ELSE IF(XOPT .EQ. 'RR')THEN
	   FOUND=.FALSE.
	   CALL USR_HIDDEN(LEV,2,2,'LIMS','1,1','Imin, Imax')
	   DO ID=1,NUM_IONS
	     IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	       IF(LEV(2) .EQ. 1)THEN
	         CALL SETREC(TA,ATM(ID)%DXzV,ATM(ID+1)%XzV_F,ATM(ID+1)%NXzV_F,ND)
	       ELSE
	         TA=RZERO
	         DO J=1,ND
	           DO I=LEV(1),LEV(2)
	             TA(J)=TA(J)+ATM(ID+1)%XzV_F(I,J)
	           END DO
	         END DO
	       END IF
	       FOUND=.TRUE.
	     END IF
	   END DO
	   IF(.NOT. FOUND)THEN
	     WRITE(T_OUT,*)'Invalid Recombination Request'
	     GOTO 3
	   END IF
!
! Compute the recombination rate. Two recombination rates will be
! computed. The first will ignore the temperature variation, whilst the
! second will normalize the number of recombinations to 10^4K by assuming
! that the recombintion rate goes as T^{-0.8}. For d=1kpc, alpha=10^{-12}
! the recombination rate will have the units ergs/cm^2/s . (need to
! divide bu the transition wavelength in mum)
!
	   CALL USR_HIDDEN(ELEC,'TVAR','T','Include T variation?')
	   CALL USR_HIDDEN(T2,'EXP','0.8',
	1	'Exponent for T variation?')
!
	   T1=-26.68059_LDP
!
	   DO I=1,ND
	     TA(I)=CLUMP_FAC(I)*(ED(I)/1.0E+10_LDP)*TA(I)*( R(I)**3 )
	     TB(I)=TA(I)*( T(I)**(-T2) )
	   END DO
	   IF(ELEC)THEN
!
! Need to use TC to compute YV as can get floating point underflow
! at outer boundary.
!
	     TC(1)=TB(1)
	     DO I=2,ND
	       TC(I)=TC(I-1)+(TB(I)+TB(I-1))*
	1            ( LOG(R(I-1))-LOG(R(I)) )*0.5
	       YV(I-1)=LOG10(TC(I-1))+T1
	     END DO
	     YV(ND)=LOG10(TC(ND))+T1
	     CALL DP_CURVE(ND,XV,YV)
	      YAXIS='Log(c.Ne.X\uN+\d.T\u0.8\d)'
	   ELSE
	     TC(1)=TA(1)
	     DO I=2,ND
	       TC(I)=TC(I-1)+(TA(I)+TA(I-1))*
	1            ( LOG(R(I-1))-LOG(R(I)) )*0.5
	       YV(I-1)=LOG10(TC(I-1))+T1
	     END DO
	     YV(ND)=LOG10(TC(ND))+T1
	     CALL DP_CURVE(ND,XV,YV)
	     YAXIS='Log(c.Ne.X\uN+\d)'
	   END IF
!
! 
!
! Plots the Rayleigh scattering opacity as a function of depth.
!
	ELSE IF(XOPT .EQ. 'RAY')THEN
	  CALL USR_OPTION(FREQ,'LAM','0.0',FREQ_INPUT)
	  IF(FREQ .EQ. 0)THEN
	  ELSE IF(KEV_INPUT)THEN
	      FREQ=FREQ*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	      FREQ=ANG_TO_HZ/FREQ
	  END IF
          IF(ATM(1)%XzV_PRES)THEN
	    CHI_RAY(1:ND)=RZERO
            CALL RAYLEIGH_SCAT(CHI_RAY,ATM(1)%XzV_F,ATM(1)%AXzV_F,ATM(1)%EDGEXZV_F,
	1             ATM(1)%NXzV_F,FREQ,ND)
	    CHI_RAY(1:ND)=1.0E-10_LDP*CHI_RAY(1:ND)*CLUMP_FAC(1:ND)
          END IF
	  CALL DP_CURVE(ND,XV,CHI_RAY)
	  YAXIS='\gx(cm\u-1\d)'
!
! Plots the Rayleigh scattering cross-section as a function of opacity.
! In units of the Thompson-cross section.
!
	ELSE IF(XOPT .EQ. 'RAYL')THEN
	  CALL USR_OPTION(LAM_ST,'LAM_ST','100.0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAM_EN','10000.0',FREQ_INPUT)
	  CALL USR_OPTION(K,'NPNTS=','1000','Number of points')
	  K=MIN(K,N_PLT_MAX)
	  T2=EXP(LOG(LAM_EN/LAM_ST)/(K-1))
	  T1=LAM_ST
	  I=0
	  DO WHILE(I .LT. K)
	    I=I+1
	    FREQ=ANG_TO_HZ/T1
	    WV(I)=T1
            IF(ATM(1)%XzV_PRES)THEN
	      YV(I)=RZERO
              CALL RAYLEIGH_SCAT(YV(I),ATM(1)%XzV_F,ATM(1)%AXzV_F,ATM(1)%EDGEXZV_F,
	1             ATM(1)%NXzV_F,FREQ,IONE)
	    ELSE
	      GOTO 1
            END IF
	    T1=T1*T2
	  END DO
	  YV(1:K)=YV(1:K)/ATM(1)%XzV_F(1,1)/6.65E-15_LDP
	  CALL DP_CURVE(K,WV,YV)
	  XAXIS='\gl(\A)'
	  YAXIS='\gs/\gs\dT\u'
!
	ELSE IF(XOPT .EQ. 'ETA')THEN
	  ELEC=.FALSE.
	  CALL USR_OPTION(ELEC,'RCUBE','F','Multiply ETA by R^3?')
	  IF(ELEC)THEN
	    DO I=1,ND
	      YV(I)=LOG10(R(I)*R(I)*R(I)*ETA(I)+1.0E-250_LDP)+20.0_LDP
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(r\u3\d\ge(ergs/s/Hz)'
	  ELSE
	    DO I=1,ND
	      YV(I)=LOG10(ETA(I)+1.0E-250_LDP)-10.0_LDP
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\ge(ergs/cm\u3\d/s/Hz)'
	  END IF
!
	ELSE IF(XOPT .EQ. 'DIFFT')THEN
	  ELEC=.FALSE.
	  CALL USR_OPTION(ELEC,'DB','T','Diffusion time to outer boundary?')
	  DO I=1,ND
	    CHIROSS(I)=CLUMP_FAC(I)*ROSS_MEAN(I)
	  END DO
	  CALL TORSCL(TAUROSS,CHIROSS,R,TB,TC,ND,METHOD,TYPE_ATM)
	  WRITE(T_OUT,*)'Rossland optical depth is : ',TAUROSS(ND)
	  DO I=1,ND-1
	    YV(I)=(RONE+TAUROSS(I+1)-TAUROSS(I))*(R(I)-R(I+1))*1.0E+05_LDP/C_KMS/24.0_LDP/3600.0_LDP
	  END DO
	  J=ND-1
	  IF(ELEC)THEN
	    ZV(1)=RZERO
	    DO I=2,ND
	      ZV(I)=ZV(I-1)+YV(I-1)
	    END DO
	    CALL DP_CURVE(J,XV,ZV)
	    YAXIS='Diffusion time to outer boundary (days)'
	  ELSE
	    WRITE(6,*)'Computed diffusion time between grid points'
	    CALL DP_CURVE(J,XV,YV)
	    YAXIS='Diffusion time (days)'
	  END IF
!
	ELSE IF(XOPT .EQ.'OP' .OR.
	1       XOPT .EQ. 'KAPPA' .OR.
	1       XOPT .EQ. 'TAUC' .OR.
	1       XOPT .EQ. 'DTAUC') THEN
!
	  IF(.NOT. ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
!
	  IF(XOPT .EQ. 'TAUC')THEN
	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	    DO I=1,ND
	      YV(I)=LOG10(TA(I))
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log(\gt)'
	  ELSE IF(XOPT .EQ. 'DTAUC')THEN
	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	    DO I=1,ND-1
	      IF(TA(I+1) .GT. TA(I))THEN
	        YV(I)=LOG10(TA(I+1)-TA(I))
	      ELSE
	         YV(I)=10
	      END IF
	    END DO
	    I=ND-1
	    CALL DP_CURVE(I,XV,YV)
	    YAXIS='Log(\gD\gt)'
	  ELSE IF(XOPT .EQ. 'KAPPA')THEN
	     YV(1:ND)=1.0E-10_LDP*CHI(1:ND)/MASS_DENSITY(1:ND)/CLUMP_FAC(1:ND)
	     CALL DP_CURVE(ND,XV,YV)
	     YAXIS='\gk(cm\u3 \d/g)'
	  ELSE
!
! We subtract 10 to put CHI in units of cm-1
!
	    WRITE(6,*)'Warning: The plotted CHI has units of cm^-1 - not program units'
	    DO I=1,ND
	      YV(I)=LOG10(CHI(I))-10.0
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    YAXIS='Log \gx(cm\u-1\d)'
	  END IF
!
	ELSE IF(XOPT .EQ. 'TWOOP')THEN
	  CALL USR_OPTION(LAM_ST,'LAMST','911.0D0',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAMEN','9000.0D0',FREQ_INPUT)
	  IF(KEV_INPUT)THEN
	    LAM_ST=LAM_ST*KEV_TO_HZ
	    LAM_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/LAM_ST
	    LAM_EN=ANG_TO_HZ/LAM_EN
	  END IF
	  CALL USR_OPTION(NFREQ,'NPTS','-100','+ve lin. spacing, -ve log')
!
	  IF(NFREQ .LT. 0)THEN
	    LINX=.FALSE.
	    NFREQ=-NFREQ
	    T1=LOG(LAM_EN/LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      XNU(I)=EXP( LOG(LAM_ST)+(I-1)*T1 )
	    END DO
	  ELSE
	    LINX=.TRUE.
	    T1=(LAM_EN-LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      XNU(I)=LAM_ST+(I-1)*T1
	    END DO
	  END IF
!
	  DO I=1,NUM_IONS
	    IF(XSPEC .EQ. UC(ION_ID(I)))THEN
	      ID=I
	      EXIT
	    END IF
	  END DO
	  DEFAULT='LTE'
	  CALL USR_OPTION(TWO_PHOT_OPTION,'OPT','LTE','Two photon option: LTE, NOSTIM, RAD, and OLD_DEF')
	  TWO_PHOT_OPTION=UC(TWO_PHOT_OPTION)
!
	  XAXSAV=XAXIS
	  CALL USR_OPTION(ELEC,'OPAC','F','Plot opacity instead of emissivity?')
	  DO ML=1,NFREQ
	    ETA(1:ND)=RZERO; CHI(1:ND)=RZERO
	    CALL TWO_PHOT_OPAC_DISP_V3(ETA,CHI,ATM(ID)%XzV_F,T,XNU(ML),TWO_PHOT_OPTION,ID,ND,ATM(ID)%NXzV_F)
	    IF(ELEC)THEN
	      YV(ML)=CHI(1)
	      WV(ML)=CHI(ND/2)
	      ZV(ML)=CHI(ND)
	    ELSE
	      YV(ML)=ETA(1)
	      WV(ML)=ETA(ND/2)
	      ZV(ML)=ETA(ND)
	    END IF
	  END DO
	  IF(ELEC)THEN
	    YAXIS='Opacity'
	  ELSE
	    YAXIS='Emissivity'
	  END IF
!
	  CALL USR_OPTION(ELEC,'PHOT','F','Plot photon emission rate instead of emissivity?')
	  IF(ELEC)THEN
	    DO ML=1,NFREQ
	      T1=4.0_LDP*PI/6.626E-27_LDP/1.0E+15_LDP/XNU(ML)/1.0E+10_LDP
	      YV(ML)=YV(ML)*T1/ATM(ID)%XzV_F(2,1)
	      WV(ML)=WV(ML)*T1/ATM(ID)%XzV_F(2,ND/2)
	      ZV(ML)=ZV(ML)*T1/ATM(ID)%XzV_F(2,ND)
	    END DO
	    YAXIS='# of photons/s/Hz'
	  END IF
!
	  XAXIS='\gv(10\u15\d Hz)'
	  CALL DP_CURVE(NFREQ,XNU,YV)
	  CALL DP_CURVE(NFREQ,XNU,WV)
	  CALL DP_CURVE(NFREQ,XNU,ZV)
!
! Estimate the effective absorbative optical depth scale for Gamma-rays.
!
	ELSE IF(XOPT .EQ. 'TAUGAM' .OR. XOPT .EQ. 'GAMABS')THEN
	  TA(1:ND)=0.5_LDP			!Number of electrons per baryon
	  DO ISPEC=1,NSPEC
	    IF('HYD' .EQ. SPECIES(ISPEC))THEN
	      TA(1:ND)=0.5_LDP*(RONE+POPDUM(1:ND,ISPEC)/POP_ATOM(1:ND))
	    END IF
	  END DO
	  XM(1:ND)=0.06_LDP*TA(1:ND)*MASS_DENSITY(1:ND)*1.0E+10_LDP
	  WRITE(6,*)XM(ND)*R(ND),XM(ND)*(R(ND-1)-R(ND))
!
	  IF(XOPT .EQ. 'GAMABS')THEN
	    CALL GAM_ABS(R,V,XM,ND)
	    YAXIS='Rad. Dec. Energy (ergs/cm\u2\d/sec)'
!
	  ELSE
	    CALL TORSCL(TA,XM,R,TB,TC,ND,METHOD,TYPE_ATM)
	    CALL USR_OPTION(ELEC,'dTAU','F','Plot dTAU)')
	    IF(ELEC)THEN
	      DO I=1,ND-1
	        YV(I)=LOG10(TA(I+1)/TA(I))
	      END DO
	      CALL DP_CURVE(ND-1,XV,YV)
	      YAXIS='dLog \gt\d\g\u'
	    ELSE
	      DO I=1,ND
	        YV(I)=LOG10(TA(I))
	      END DO
	      CALL DP_CURVE(ND,XV,YV)
	      YAXIS='Log(\gt\d\g\u)'
	    END IF
	  END IF
!
! These options allow you to plot tau at a particular R (TAUR), or
! alternatively, R at a given value of Tau (RTAU).
!
! We can also write out the continuum opacities as a function of
! wavelength.
!
! To be read TAU_at_R and R_at_TAU respectively.
!
	ELSE IF(XOPT .EQ. 'RTAU'  .OR. XOPT .EQ. 'VTAU' .OR. XOPT .EQ. 'TAUR' .OR.
	1       XOPT .EQ. 'EDTAU' .OR.
	1       XOPT .EQ. 'KAPR'  .OR. XOPT .EQ. 'CHIR' .OR. XOPT .EQ. 'ALBEDO' .OR.
	1       XOPT .EQ. 'ETAR'  .OR. XOPT .EQ. 'WROPAC')THEN
	  IF(XRAYS)WRITE(T_OUT,*)'Xray opacities (i.e. K shell) are included'
	  IF(.NOT. XRAYS)WRITE(T_OUT,*)'Xray opacities (i.e. K shell) are NOT included'
	  CALL USR_OPTION(LAM_ST,'LAMST',' ',FREQ_INPUT)
	  CALL USR_OPTION(LAM_EN,'LAMEN',' ',FREQ_INPUT)
	  IF(KEV_INPUT)THEN
	    LAM_ST=LAM_ST*KEV_TO_HZ
	    LAM_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    LAM_ST=ANG_TO_HZ/LAM_ST
	    LAM_EN=ANG_TO_HZ/LAM_EN
	  END IF
	  CALL USR_OPTION(NFREQ,'NPTS','10','+ve lin. spacing, -ve log')
!
	  IF(NFREQ .LT. 0)THEN
	    LINX=.FALSE.
	    NFREQ=-NFREQ
	    T1=LOG(LAM_EN/LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      Q(I)=EXP( LOG(LAM_ST)+(I-1)*T1 )
	      ZV(I)=Q(I) 			!Use ZV .: don't corrupt XV
	    END DO
	  ELSE
	    LINX=.TRUE.
	    T1=(LAM_EN-LAM_ST)/(NFREQ-1)
	    DO I=1,NFREQ
	      Q(I)=LAM_ST+(I-1)*T1
	      ZV(I)=Q(I) 			!Use ZV .: don't corrupt XV
	    END DO
	  END IF
	  CALL USR_OPTION(ELEC,'ELEC','T','Include elec?')
!
	  CALL USR_HIDDEN(LINX,'LINX','F','Linear X spacing?')
	  IF(KEV_INPUT)THEN
	    DO I=1,NFREQ
	      ZV(I)=ZV(I)/KEV_TO_HZ
	    END DO
	    XAXIS='E(keV)'
	  ELSE IF(ANG_INPUT)THEN
	    DO I=1,NFREQ
	      ZV(I)=ANG_TO_HZ/ZV(I)
	    END DO
	    XAXIS='\gl(\A)'
	  ELSE
	    XAXIS='\gn(10\u15\d Hz)'
	  END IF
	  IF(.NOT. LINX)THEN
	    XAXIS='Log '//XAXIS
	    DO I=1,NFREQ
	      ZV(I)=LOG10(ZV(I))
	    END DO
	  END IF
	  CALL USR_HIDDEN(LINY,'LINY','F','Linear Y Axis')
!
	  IF(XOPT .EQ. 'WROPAC')THEN
	    WRITE(6,*)' '
	    WRITE(6,*)'Writing contiinuum opacities to CONT_OPAC '
	    WRITE(6,*)' '
	    OPEN(LU_OUT,FILE='CONT_OPAC',STATUS='UNKNOWN')
	      WRITE(LU_OUT,'(I5,T30,A)')ND,'!Number of depth points'
	      WRITE(LU_OUT,'(I5,T30,A)')NFREQ,'!Number of frequency points'
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Radius grid'
	      WRITE(LU_OUT,'(6ES14.6)')R(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Velocity(km/s)'
	      WRITE(LU_OUT,'(6ES14.6)')V(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Density (gm/cm^3)'
	      WRITE(LU_OUT,'(6ES14.6)')MASS_DENSITY(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Clumping factor'
	      WRITE(LU_OUT,'(6ES14.6)')CLUMP_FAC(1:ND)
	      WRITE(LU_OUT,'(A)')' '
	      WRITE(LU_OUT,'(A)')' Kappa Table (cm^2/gm)'
	  ELSE IF(XOPT .EQ. 'TAUR' .OR. XOPT .EQ. 'CHIR' .OR. XOPT .EQ. 'KAPR' .OR. XOPT .EQ. 'ETAR')THEN
	    CALL USR_OPTION(RVAL,'RAD',' ','Radius in R* (-ve for depth index)')
	    IF(RVAL .GT. 0)THEN
	      RVAL=RVAL*R(ND)
	      IF(RVAL .GT. R(1))RVAL=R(1)
	      IF(RVAL .LT. R(ND))RVAL=R(ND)
	      R_INDX=ND
	      DO WHILE(RVAL .GT. R(R_INDX))
	        R_INDX=R_INDX-1
	      END DO
	      R_INDX=MIN(R_INDX,ND-1)
	    ELSE
	      R_INDX=NINT(ABS(RVAL))
	      R_INDX=MAX(1,R_INDX)
	      R_INDX=MIN(R_INDX,ND)
	      RVAL=R(R_INDX)
	      R_INDX=MIN(R_INDX,ND-1)
	    END IF
	    WRITE(6,*)'   RVAL is',RVAL
	    WRITE(6,*)'RVAL/R* is',RVAL/R(ND)
	    YAXIS='Log(\gt[R])'
	    IF(LINY)YAXIS='\gt[R]'
	    IF(XOPT .EQ. 'CHIR')YAXIS='\gx (cm\u-1\d)'
	    IF(XOPT .EQ. 'KAPR')YAXIS='\gk (cm\u2\d/g)'
	    IF(XOPT .EQ. 'ETAR')YAXIS='\ge (ergs cm\u3 \ds\u-1\d)'
	  ELSE IF(XOPT .EQ. 'VTAU')THEN
	    CALL USR_OPTION(TAU_VAL,'TAU',' ','Tau value for which V is to be determined')
	    YAXIS='V(km/s)'
	  ELSE IF(XOPT .EQ. 'EDTAU')THEN
	    CALL USR_OPTION(TAU_VAL,'TAU',' ','Tau value for which Ne is to be determined')
	    YAXIS='Ne(/cm\u3\d)'
	  ELSE
	    CALL USR_OPTION(TAU_VAL,'TAU',' ','Tau value for which R is to be determined')
	    YAXIS='Log(R[\gt]/R\d*\u)'
	    IF(LINY)YAXIS='R[\gt]/R\d*\u'
	  END IF
!
	  CALL USR_HIDDEN(IN_R_SUN,'IN_R_SUN','F','Y In R_SUN?')
	  IF(IN_R_SUN)THEN
	    IF(LINY)YAXIS='R/R\dO\u'
	    IF(.NOT. LINY)YAXIS='Log(R/R\dO\u)'
	  END IF
!
	  DO ML=1,NFREQ
!
! Compute continuum opacity and emissivity at the line frequency.
!
	    FL=Q(ML)
	    INCLUDE 'OPACITIES.INC'
!	    INCLUDE 'HOME:[jdh.cmf.carb.test]OPACITIES.INC'
!
	    IF(.NOT. ELEC)THEN
	      DO I=1,ND
	        CHI(I)=CHI(I)-ESEC(I)
	      END DO
	    END IF
!
	    IF(XOPT .EQ. 'WROPAC')THEN
	      WRITE(LU_OUT,'(/,ES16.6)')ANG_TO_HZ/FL
	      WRITE(LU_OUT,'(6ES14.4)')1.0D-10*CHI(1:ND)/MASS_DENSITY(1:ND)
	    ELSE IF(XOPT .EQ. 'ETAR')THEN
	      T1=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	      T2=ETA(R_INDX+1)*CLUMP_FAC(R_INDX+1)
	      T3=ETA(R_INDX)*CLUMP_FAC(R_INDX)
	      YV(ML)=1.0E-10_LDP*( T1*T2 + (1.0_LDP-T1)*T3 )
	      IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	    ELSE IF(XOPT .EQ. 'CHIR')THEN
	      T1=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	      T2=CHI(R_INDX+1)*CLUMP_FAC(R_INDX+1)
	      T3=CHI(R_INDX)*CLUMP_FAC(R_INDX)
	      YV(ML)=1.0E-10_LDP*( T1*T2 + (1.0_LDP-T1)*T3 )
	      IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	    ELSE IF(XOPT .EQ. 'KAPR')THEN
	      T1=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	      T2=CHI(R_INDX+1)/MASS_DENSITY(R_INDX+1)
	      T3=CHI(R_INDX)/MASS_DENSITY(R_INDX)
	      YV(ML)=1.0E-10_LDP*( T1*T2 + (1.0_LDP-T1)*T3 )
	      IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	    ELSE
!
! Adjust opacities for the effect of clumping.
!
	      DO I=1,ND
	        ETA(I)=ETA(I)*CLUMP_FAC(I)
	        CHI(I)=CHI(I)*CLUMP_FAC(I)
	        ESEC(I)=ESEC(I)*CLUMP_FAC(I)
	        CHI_RAY(I)=CHI_RAY(I)*CLUMP_FAC(I)
	        CHI_SCAT(I)=CHI_SCAT(I)*CLUMP_FAC(I)
	      END DO
!
	      CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	      IF(XOPT .EQ. 'ALBEDO')THEN
	        I=1
	        DO WHILE(TAU_VAL .GT. TA(I) .AND. I .LT. ND)
	          I=I+1
	        END DO
	        IF(TAU_VAL .GT. TA(ND-1))THEN
	          YV(ML)=ESEC(ND)/CHI(ND)
	        ELSE IF(TAU_VAL .LE. TA(1))THEN
	          YV(ML)=ESEC(1)/CHI(1)
	        ELSE
	          T2=(TA(I)-TAU_VAL)/(TA(I)-TA(I-1))
	          YV(ML)=( (1.0_LDP-T2)*ESEC(I)+T2*ESEC(I-1) )/( (1.0_LDP-T2)*CHI(I)+T2*CHI(I-1) )
	        END IF
	        YAXIS='Albedo'
	      ELSE
	        IF(XOPT .EQ. 'TAUR')THEN
	          T2=(R(R_INDX)-RVAL)/(R(R_INDX)-R(R_INDX+1))
	          YV(ML)=T2*TA(R_INDX+1) + (1.0_LDP-T2)*TA(R_INDX)
	          IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	        ELSE IF(XOPT .EQ. 'VTAU' .OR. XOPT .EQ. 'EDTAU')THEN
	          TB(1:ND)=V(1:ND)
	          IF(XOPT .EQ. 'EDTAU')TB(1:ND)=ED(1:ND)
	          I=1
	          DO WHILE(TAU_VAL .GT. TA(I) .AND. I .LT. ND)
	            I=I+1
	          END DO
	          IF(TAU_VAL .GT. TA(ND-1))THEN
	            YV(ML)=TB(ND)
	          ELSE IF(TAU_VAL .LE. TA(1))THEN
	            YV(ML)=TB(1)*TA(1)/TAU_VAL
	          ELSE
	            T2=(TA(I)-TAU_VAL)/(TA(I)-TA(I-1))
                    YV(ML)=( (1.0_LDP-T2)*TB(I)+T2*TB(I-1) )
	          END IF
	        ELSE
	          I=1
	          DO WHILE(TAU_VAL .GT. TA(I) .AND. I .LT. ND)
	            I=I+1
	          END DO
	          IF(TAU_VAL .GT. TA(ND-1))THEN
	            YV(ML)=1.0
	          ELSE IF(TAU_VAL .LE. TA(1))THEN
	            YV(ML)=R(1)*TA(1)/TAU_VAL/R(ND)
	          ELSE
	            T2=(TA(I)-TAU_VAL)/(TA(I)-TA(I-1))
                    YV(ML)=( (1.0_LDP-T2)*R(I)+T2*R(I-1) )/R(ND)
	          END IF
	          IF(IN_R_SUN)YV(ML)=YV(ML)*R(ND)/6.96_LDP
	          IF(.NOT. LINY)YV(ML)=LOG10(YV(ML))
	        END IF
	      END IF
	    END IF
	  END DO
	  IF(XOPT .NE. 'WROPAC')CALL DP_CURVE(NFREQ,ZV,YV)
!
	ELSE IF(XOPT .EQ. 'INTERP')THEN
!
! Routine interpolates the last varible plotted. Variable can only
! be plotted against R. Crude, but effective.
! Routine assumes YV is in the correct form for the interpolation.
!
	  DO I=1,ND
	    TA(I)=YV(I)
	  END DO
!
! Do the interpolation - Based on INTERPTHREE.
!
	  DO I=1,NDX
	    TB(I)=RZERO
	    DO J=0,3
	      TB(I)=TB(I)+COEF(J,I)*TA(J+INDX(I))
	    END DO
	    YV(I)=TB(I)
	  END DO
!
	  DO I=1,NDX
	    XNU(I)=LOG10( REXT(I)/REXT(NDX) )
	  END DO
	  CALL DP_CURVE(NDX,XNU,YV)
!
! Option computes a new radius grid, equally spaced in Log(Tau).
!
	ELSE IF(XOPT .EQ. 'NEWRG')THEN
!
! Compute the optical depth scale, which is stored, in LOG form, in TA.
!
	  IF(.NOT. ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
	  CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
	  DO I=1,ND
	    TA(I)=LOG10(TA(I))
	  END DO
!
	  DEFAULT=WR_STRING(ND)
	  CALL USR_OPTION(ND_TMP,'ND=',DEFAULT,'Number of radal depths')
	  IF(ND_TMP .GT. NP_MAX)THEN
	    WRITE(T_OUT,*)'Error in MAINGEN'
	    WRITE(T_OUT,*)'Value for ND is too large'
	    WRITE(T_OUT,*)'Maximum value for ND is',ND
	    GOTO 1
	  END IF
!
! Compute the new optical depth scale, equally spaced in Log(TAU), with
! a finer spacing near both boundaies.
!
	  DTAU(1)=(TA(ND)-TA(1))/(ND_TMP-5)
	  TB(1)=TA(1)
	  TB(2)=TA(1)+DTAU(1)/10
	  TB(3)=TA(1)+DTAU(1)/2
	  DO I=4,ND_TMP-3
	    TB(I)=TA(1)+(I-3)*DTAU(1)
	  END DO
	  TB(ND_TMP)=TA(ND)
	  TB(ND_TMP-1)=TA(ND)-DTAU(1)/10
	  TB(ND_TMP-2)=TA(ND)-DTAU(1)/2
!
	  CALL MON_INTERP(TC,ND_TMP,IONE,TB,ND,R,ND,TA,ND)
	  CALL MON_INTERP(XM,ND_TMP,IONE,TB,ND,V,ND,TA,ND)
	  CALL MON_INTERP(RJ,ND_TMP,IONE,TB,ND,SIGMA,ND,TA,ND)
!
	  WRITE(T_OUT,*)'Plotting Tau, V, Sigma versus Log(r/RP).'
	  WRITE(T_OUT,*)'Plotting Tau, V, Sigma versus depth index.'
	  WRITE(T_OUT,*)'Choose plot in plot package.'

	  T1=TC(ND)
	  DO I=1,ND
	    XNU(I)=LOG10(TC(I)/T1)
	  END DO
	  CALL DP_CURVE(NDX,XNU,TC)
	  CALL DP_CURVE(NDX,XNU,XM)
	  CALL DP_CURVE(NDX,XNU,RJ)
!
	  DO I=1,ND
	    XNU(I)=I
	  END DO
	  CALL DP_CURVE(NDX,XNU,TC)
	  CALL DP_CURVE(NDX,XNU,XM)
	  CALL DP_CURVE(NDX,XNU,RJ)
!
	  CALL GEN_ASCI_OPEN(27,'NEW_R_GRID','UNKNOWN',' ',' ',IZERO,IOS)
	    WRITE(27,*)ND_TMP
	    WRITE(27,'(3X,1P,3E17.8)')(TC(I),XM(I),RJ(I),I=1,ND_TMP)
	  CLOSE(UNIT=27)
! 
!
!
	ELSE IF(XOPT .EQ. 'EW')THEN
!
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	  END DO
!
! TA is used for the line flux. Integral of TA log(r) is
! the line EW.
!
	  FORCE_MULT(1:ND)=RZERO
	  CALL SOBEW_GRAD_V2(SOURCE,CHI,CHI_SCAT,CHIL,ETAL,
	1             V,SIGMA,R,P,FORCE_MULT,LUM,
	1             JQW,HQW,TA,T1,S1,
	1             FREQ,'HOLLOW',DBB,IC,THICK,.FALSE.,NC,NP,ND,METHOD)
!
          OPEN(UNIT=18,FILE='DSOB_FORCE_MULT',STATUS='UNKNOWN')
	    WRITE(18,'(3X,A1,10X,A1,15X,A1,13X,A1)')'I','R','V','M'
            DO I=1,ND
              WRITE(18,'(1X,I3,3X,3ES14.6)')I, R(I),V(I),FORCE_MULT(I)
            END DO
          CLOSE(UNIT=18)
!
! S1 is the continuum flux in Jy for an object at 1kpc.
! T1 is the line equivalent width in Angstroms.
! T3 is the line flux in ergs/cm^2/s
!
	  T2=LAMVACAIR(FREQ)		!Wavelength(Angstroms)
	  T3=T1*S1*1.0E-23_LDP*FREQ*1.0E+15_LDP/T2
	  T4=T3*4*PI*(3.0856E+21_LDP)**2/LUM_SUN()
	  WRITE(LU_NET,40008)T1,S1,T2,T3,T4
	  WRITE(T_OUT,'(A)')RED_PEN
	  WRITE(T_OUT,40008)T1,S1,T2,T3,T4
	  WRITE(T_OUT,'(A)')DEF_PEN
!
	  CALL USR_OPTION(ELEC,'PLOT','T','PLot Line Origin?')
!
! We can now use the current X-axis for the X-axis. Zeta(x)
! is defined similarly to before: Zeta(x) dx gives the
! emission about x in an interval dx.
!
! We use JNU for COEF which must be dimensioned JNU(ND,*) [with * > 3]
! NB: T1 contains the line EW.
!
	  IF(ELEC)THEN
	    WRITE(T_OUT,*)' '
	    WRITE(T_OUT,*)BLUE_PEN//'Area under curve is normalized to unity'//DEF_PEN
	    WRITE(T_OUT,*)' '
	    DO I=1,ND-2
	      IF( (XV(I+1)-XV(I))*(XV(I+2)-XV(I+1)) .LE. 0.0_LDP)THEN
	        WRITE(T_OUT,*)RED_PEN//'Error -- can only plot line origin against a mononotonix X-axis'//DEF_PEN
	        GOTO 1
	      END IF
	    END DO
	    CALL MON_INT_FUNS_V2(JNU,XV,R,ND)
	    T4=RONE; IF(XV(1) .LT. XV(4))T4=-RONE
	    DO I=1,ND
	      YV(I)=T4*TA(I)/R(I)/JNU(I,3)/T1
	    END DO
	    CALL DP_CURVE_LAB(ND,XV,YV,TRANS_NAME)
	    YAXIS='\gz'
	  ELSE
	    CALL DP_CURVE_LAB(ND,XV,FORCE_MULT,TRANS_NAME)
	    YAXIS='Force Multiplier'
	  END IF
!
	ELSE IF(XOPT .EQ. 'EP' .OR. XOPT .EQ. 'BETA')THEN
!
! Used to indicate where the line emission is coming from.
! It should be plotted against log(R) . Area under curve is
! emission. Takes continuous opacity (not .e.s.) into account.
! This section requires tha the line and continuous opacity have
! been previously computed.
!
	  CALL USR_HIDDEN(ELEC,'ELEC','T',
	1      'Ignore Electron Scattering?')
	  IF(ELEC)THEN
	    DO I=1,ND
	      CHI(I)=CHI(I)-ESEC(I)
	    END DO
	  END IF
	  CALL USR_HIDDEN(ELEC,'NEW','T','Continuum included?')
	  CALL USR_HIDDEN(LINV,'LINV','F','Linear R Scale?')
!
! We are using ETA for BETA, the escape probability.
!
	  IF(ELEC)THEN
	    CALL BETANEW(CHI,CHIL,SIGMA,ETA,R,Z,V,FREQ,ND)
	  ELSE
!	    CALL TORSCL(TA,CHI,R,TB,TC,ND,METHOD,TYPE_ATM)
 	    CALL WRBETA(CHIL,SIGMA,ETA,R,V,FREQ,ND)
	  END IF
!
	  IF(XOPT .EQ. 'EP')THEN
	    YAXIS='\gc'
	    IF(LINV)THEN
	      DO I=1,ND
	        ZV(I)=R(I)/R(ND)
	        YV(I)=ETA(I)*R(I)*R(I)*ETAL(I)
	      END DO
	      XAXIS='r/R\d*\u'
	    ELSE
	      DO I=1,ND
	        ZV(I)=LOG10(R(I)/R(ND))
!	        YV(I)=ETA(I)*R(I)*R(I)*R(I)*ETAL(I)*EXP(-TA(I))
	        YV(I)=ETA(I)*R(I)*R(I)*R(I)*ETAL(I)
	      END DO
	      XAXIS='Log(r/R\d*\u)'
	    END if
	    T1=0.0
	    DO I=1,ND-1
	      T1=T1+0.5_LDP*(YV(I)+YV(I+1))*(ZV(I)-ZV(I+1))
	    END DO
!
! Normalize to unit area.
!
	    CALL USR_HIDDEN(SCALE,'SCL','T','Your scalinge for EP?')
	    IF(.NOT. SCALE)THEN
	      WRITE(T_OUT,*)'SCALE PARAM=',T1
	      CALL USR_OPTION(T1,'T1',' ','Scale parameter')
	    END IF
	    DO I=1,ND
	      YV(I)=YV(I)/T1
	    END DO
!
	  ELSE
	    Yaxis='\gb'
	    XAXIS='Log(r/R\d*\u)'
	    DO I=1,ND
	      ZV(I)=LOG10(R(I)/R(ND))
	      YV(I)=ETA(I)
	    END DO
	  END IF
	  CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'LINRC')THEN
	  DO I=1,ND
	    ZV(I)=LOG10(R(I)/R(ND))
	    YV(I)=R(I)*R(I)*R(I)*ETAL(I)
	  END DO
	  WV=RZERO
	  WV(1)=YV(1)
	  DO I=2,ND
	    WV(I)=WV(I-1)+0.5_LDP*(YV(I-1)+YV(I))*(ZV(I-1)-ZV(I))
	  END DO
	  CALL DP_CURVE(ND,XV,WV)
!
! 
!
! Require CHIL to have been computed in setup.
!
	ELSE IF(XOPT .EQ. 'TAUL')THEN
	  CALL USR_HIDDEN(ELEC,'STAT','F','Stationary opactical depth?')
	  CALL USR_HIDDEN(RADIAL,'RADS','F','Radial Sobolev optical depth?')
	  CALL USR_HIDDEN(RAY_TAU,'RAY_TAU','F','Ray Sobolev optical depth?')
	  IF(ELEC)THEN
	    CALL TORSCL(TA,CHIL,R,TB,TC,ND,METHOD,TYPE_ATM)
!
! Assumes V_D=10kms.
!
	    IF(MINVAL(TA(1:ND)) .LE. 0.0_LDP)THEN
	      WRITE(6,*)'Negative optical depth encountered using linear plot'
	      T1=1.6914E-11_LDP/FREQ
	      DO I=1,ND
	         YV(I)=T1*TA(I)
		 WRITE(6,*)I,TA(I)
	      END DO
	      YAXIS='Log(\gt\dstat\u)'
	    ELSE
	      T1=LOG10(1.6914E-11_LDP/FREQ)
	      DO I=1,ND
	        IF(TA(I) .GT. 0)THEN
	          YV(I)=T1+LOG10(TA(I))
	        ELSE IF(TA(I) .LT. 0)THEN
	          YV(I)=T1+LOG10(-TA(I))-20.0_LDP
	        ELSE
	          YV(I)=-30.0
	        END IF
	      END DO
	      YAXIS='Log(\gt\dstat\u)'
	    END IF
	    DEFAULT=TRIM(TRANS_NAME)//'(stat)'
	    CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
	  ELSE IF(RAY_TAU)THEN
	    CALL IMPAR(P,R,R(ND),NC,ND,NP)
	    WRITE(6,*)'NC,ND=',NC,ND
	    CALL USR_OPTION(K,'RAY_INDX','20','Depth indx (-ve for rays 1 to NC')
	    IF(K .LT. 0)THEN
	      T1=P(ABS(K))*P(ABS(K))
	      K=ABS(K)
	    ELSE
	      T1=R(K)*R(K)
	    END IF
	    DO I=1,K
	      ZV(I)=SQRT(R(I)*R(I)-T1)
	      T2=ZV(I)/R(I)
	      YV(I)=CHIL(I)*R(I)*2.998E-10_LDP/FREQ/V(I)
	      YV(I)=YV(I)/(RONE+T2*T2*SIGMA(I))
	    END DO
	    DEFAULT=TRIM(TRANS_NAME)
	    CALL DP_CURVE_LAB(K,ZV,YV,DEFAULT)
	  ELSE
	    DO I=1,ND
	      YV(I)=CHIL(I)*R(I)*2.998E-10_LDP/FREQ/V(I)
	      IF(RADIAL)YV(I)=YV(I)/(RONE+SIGMA(I))
	    END DO
	    IF(MINVAL(YV(1:ND)) .LE. 0)THEN
	      WRITE(6,*)'Use C option in PGPLT with LG to plot on log scale'
	      YAXIS='\gt\dSob\u'
	    ELSE
	      DO I=1,ND
	        YV(I)=LOG10(YV(I))
	      END DO
	      YAXIS='Log(\gt\dSob\u)'
	    END IF
	    DEFAULT=TRIM(TRANS_NAME)
	    IF(RADIAL)DEFAULT=TRIM(TRANS_NAME)//'(radial)'
	    CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
	  END IF
	
!
	ELSE IF(XOPT .EQ. 'TAULIP')THEN
	  CALL IMPAR(P,R,R(ND),NC,ND,NP)
	  WRITE(6,*)'NC,ND=',NC,ND
	  CALL USR_OPTION(I,'RAY_INDX','20','Depth indx (-ve for rays 1 to NC')
	  CALL USR_OPTION(VSHIFT,'VSHIFT','0','Velocity for line center in km/s')
	  CALL USR_OPTION(XAXIS_OPT,'XAXIS','NZ','What X axisi (R,NR,Z,NZ,uV,I)?')
	  CALL USR_OPTION(VDOP,'VDOP','10.0D0','Doppler veolicity')
	  DEL_V=VDOP/2.0_LDP
	  XAXIS_OPT=UC(XAXIS_OPT)
	  IF(I .LT. 0)THEN
	    I=ABS(I)
	    T1=V(ND)
	  ELSE
	    T1=V(I)
	    I=NC+(ND-I)+1
	  END IF
	  WRITE(6,*)'P, V=',P(I),T1
	  CALL SET_FINE_RAY_GRID(ETAL,CHIL,
	1     R,V,SIGMA,MASS_DENSITY,P(I),DEL_V,ND)	
	  CALL NON_SOB_TAUL(VSHIFT,VDOP,FREQ,METHOD,TYPE_ATM,XAXIS,XAXIS_OPT)
!
	ELSE IF(XOPT .EQ. 'MTAULIP')THEN
	  CALL IMPAR(P,R,R(ND),NC,ND,NP)
	  CALL USR_OPTION(VDOP,'VDOP','10.0D0','Doppler veolicity')
	  DEL_V=VDOP/2.0_LDP
	  CALL USR_OPTION(VSHIFT,'VSHIFT','0','Velocity for line center in km/s')
	  DO I=1,NP-2
	    T1=P(I)
	    IF(I .EQ. 1)T1=0.1_LDP*P(2)
	    CALL SET_FINE_RAY_GRID(ETAL,CHIL,
	1       R,V,SIGMA,MASS_DENSITY,T1,DEL_V,ND)	
	    CALL MAX_NON_SOB_TAUL(YV(I),VSHIFT,VDOP,FREQ,METHOD,TYPE_ATM)
	    WV(I)=T1/R(ND)
	  END DO
	  XAXIS='P/R(ND)'
	  CALL DP_CURVE(NP-2,WV,YV)
!
! Require CHIL to have been computed in setup.
!
	ELSE IF(XOPT .EQ. 'CHIL')THEN
!
! Assumes V_D=10kms.
!
	  CALL USR_OPTION(ELEC,'KAPPA','T','Plot kappa insted of chi?')
	  WRITE(6,*)'A Doppler velocity of 10 km/s is assumed'
	  WRITE(6,*)'Warning: CHIL now has units of cm^-1 - not program units'
	  IF(ELEC)THEN
	    T1=1.6914E-21_LDP/FREQ
	    YV(1:ND)=T1*CHIL(1:ND)/MASS_DENSITY(1:ND)/CLUMP_FAC(1:ND)
	    YAXIS='\gk(cm\u3 \d/g)'
	    CALL DP_CURVE(ND,XV,YV)
	  ELSE
	    VALID_VALUE=.TRUE.
	    ZV(1:ND)=RZERO; YV(1:ND)=RZERO
	    T1=LOG10(1.6914E-11_LDP/FREQ)-10.0_LDP
	    DO I=1,ND
	      IF(CHIL(I) .GT. 0)THEN
	        YV(I)=T1+LOG10(CHIL(I))
	      ELSE IF(CHIL(I) .LT. 0)THEN
	        ZV(I)=T1+LOG10(-CHIL(I))
	        VALID_VALUE=.FALSE.
	      ELSE
	        YV(I)=-30.0
	      END IF
	    END DO
	    YAXIS='Log \gx\dL\u(cm\u-1\d)'
	    CALL DP_CURVE(ND,XV,YV)
	    IF(.NOT. VALID_VALUE)THEN
	      WRITE(6,*)'Warning -- negative line opacities'
	      WRITE(6,*)'Sending to data plots to DURVE - use MARK to illustrate'
	      CALL DP_CURVE(ND,XV,ZV)
	    END IF
	  END IF
! 
!
! Write out line/continuum opacities and emissivities for use with
! the polarization codes, or profile codes.
!
	ELSE IF(XOPT .EQ. 'WRC' .OR. XOPT .EQ. 'WRL')THEN
	  CONT_INT=TWOHCSQ*(FREQ**3)/( EXP(HDKT*FREQ/T(ND))-RONE )
	  CALL USR_OPTION(FILE_FORMAT,'FORMAT',FILE_FORMAT,'Format: OLD, NEW, MULTI')
!
	  NEW_FILE=.TRUE.
	  IF(FILE_FORMAT .NE. 'OLD')THEN
	    CALL USR_OPTION(NEW_FILE,'NEW_FILE','F','Open new file?')
	  END IF
!
	  IF(NEW_FILE)THEN
	    FILE_FOR_WRL='LINEDATA'
	    FILE_PRES=.TRUE.
	    DO WHILE(FILE_PRES)
	      CALL USR_OPTION(FILE_FOR_WRL,'FILE',FILE_FOR_WRL,'Output file')
	      INQUIRE(FILE=FILE_FOR_WRL,EXIST=FILE_PRES)
	      IF(FILE_PRES)THEN
	        WRITE(6,*)'Error - file exists - enter alternate name'
	      END IF
	    END DO
	    OPEN(UNIT=LU_OUT,FILE=FILE_FOR_WRL,STATUS='NEW')
	  ELSE
	    INQUIRE(FILE=FILE_FOR_WRL,EXIST=FILE_PRES)
	    IF(.NOT. FILE_PRES)THEN
	      WRITE(6,*)'Error - file does not exist'
	      CALL USR_OPTION(FILE_FOR_WRL,'FILE',FILE_FOR_WRL,'Output file')
	    END IF
	    OPEN(UNIT=LU_OUT,FILE=FILE_FOR_WRL,STATUS='OLD',POSITION='APPEND')
	  END IF
!
	  IF(XOPT .EQ. 'WRL')THEN
	    WRITE(TRANS_NAME,'(1X,A,A,I3,A,I3,A,T25,A,T40,A)')
	1            TRIM(XSPEC),'(',LEV(2),'-',LEV(1),')'
	  ELSE
	    TRANS_NAME='Continuum'
	  END IF
!
	  IF(NEW_FILE)THEN
	    CALL GETCWD(PWD)
	    WRITE(LU_OUT,'(A)')'!'
	    WRITE(LU_OUT,'(A,A)')'! Model directory is ',TRIM(PWD)
	    WRITE(LU_OUT,'(A)')'!'
	    MOD_NAME=PWD
	    IF(LEN_TRIM(PWD) .GT. 20)THEN
	      I=LEN_TRIM(PWD)
	      DO J=I-1,1,-1
	        IF(PWD(J:J) .EQ. '/')THEN
	          MOD_NAME=PWD(MAX(I-19,J+1):I)
	          EXIT
	        END IF
	      END DO
	    END IF
	  END IF
!
	  CALL USR_OPTION(TMP_LOGICAL,'BGRID','F','Use a bigger grid')
	  IF(XOPT .EQ. 'WRC')THEN
	    ETAL=1.0E-10_LDP; CHIL=1.0E-10_LDP
	  END IF
	  IF(TMP_LOGICAL)THEN
	    CALL USR_OPTION(NPINS,'NPINS','1','0, 1 or 2')
	    NDX=(ND-1)*NPINS+ND
	    I=ND-10		!Parabolic interp for > I
	    CALL REXT_COEF_V2(REXT,COEF,INDX,NDX,R,GRID,ND,NPINS,.TRUE.,I,IONE,ND)
	    WRITE(6,*)'New R grid set'
	    CALL LOG_MON_INTERP(ETAEXT,NDX,IONE,REXT,NDX,ETA,ND,R,ND,L_FALSE,L_TRUE)
	    WRITE(6,*)'Eta defined'
	    CALL LOG_MON_INTERP(CHIEXT,NDX,IONE,REXT,NDX,CHI,ND,R,ND,L_FALSE,L_TRUE)
	    WRITE(6,*)'CHI defined'
	    CALL LOG_MON_INTERP(ESECEXT,NDX,IONE,REXT,NDX,ESEC,ND,R,ND,L_FALSE,L_TRUE)
	    WRITE(6,*)'ESEC defined'
	    CALL LOG_MON_INTERP(ETALEXT,NDX,IONE,REXT,NDX,ETAL,ND,R,ND,L_FALSE,L_TRUE)
	    WRITE(6,*)'Etal defined'
	    CALL LOG_MON_INTERP(CHILEXT,NDX,IONE,REXT,NDX,CHIL,ND,R,ND,L_FALSE,L_TRUE)
	    WRITE(6,*)'Chil defined'
	    CALL LOG_MON_INTERP(VEXT,NDX,IONE,REXT,NDX,V,ND,R,ND,L_FALSE,L_TRUE)
	    WRITE(6,*)'VEXT defined'
	    CALL LOG_MON_INTERP(SIGMAEXT,NDX,IONE,REXT,NDX,SIGMA,ND,R,ND,L_FALSE,L_FALSE)
	    WRITE(6,*)'SIGEXT defined'
	    CALL LOG_MON_INTERP(TEXT,NDX,IONE,REXT,NDX,T,ND,R,ND,L_FALSE,L_TRUE)
	    CALL LOG_MON_INTERP(MASS_DENSITYEXT,NDX,IONE,REXT,NDX,MASS_DENSITY,ND,R,ND,L_FALSE,L_TRUE)
	    CALL LOG_MON_INTERP(CLUMP_FACEXT,NDX,IONE,REXT,NDX,CLUMP_FAC,ND,R,ND,L_FALSE,L_TRUE)
	    ND_TMP=NDX
	  ELSE
	    ND_TMP=ND
	    REXT(1:ND)=R(1:ND); VEXT(1:ND)=V(1:ND); SIGMAEXT(1:ND)=SIGMA(1:ND);  TEXT(1:ND)=T(1:ND)
	    ETAEXT(1:ND)=ETA(1:ND); CHIEXT(1:ND)=CHI(1:ND); ESECEXT(1:ND)=ESEC(1:ND)
	    ETALEXT(1:ND)=ETAL(1:ND); CHILEXT(1:ND)=CHIL(1:ND)
	    MASS_DENSITYEXT(1:ND)=MASS_DENSITY(1:ND); CLUMP_FACEXT(1:ND)=CLUMP_FAC(1:ND)
	  END IF	

	  IF(FILE_FORMAT .EQ. 'NEW')THEN
	    CALL WRITE_LINE_12MAY98(TRANS_NAME,MOD_NAME,
	1          DIF,IC,FREQ,AMASS,REXT,VEXT,SIGMAEXT,TEXT,
	1          ETAEXT,CHIEXT,ESECEXT,CHILEXT,ETALEXT,ND_TMP,LU_OUT)
	  ELSE IF(FILE_FORMAT .EQ. 'MULTI')THEN
	    CALL WRITE_LINE_MULTI(TRANS_NAME,MOD_NAME,NEW_FILE,
	1          DIF,IC,FREQ,AMASS,REXT,VEXT,SIGMAEXT,TEXT,
	1          MASS_DENSITYEXT,CLUMP_FACEXT,
	1          ETAEXT,CHIEXT,ESECEXT,CHILEXT,ETALEXT,ND_TMP,LU_OUT)
	  ELSE
	    CALL WRITE_LINE_OLD(TRANS_NAME,MOD_NAME,
	1          DIF,IC,FREQ,AMASS,REXT,VEXT,SIGMAEXT,TEXT,
	1          MASS_DENSITYEXT,CLUMP_FACEXT,
	1          ETAEXT,CHIEXT,ESECEXT,CHILEXT,ETALEXT,ND_TMP,LU_OUT,L_FALSE)
	  END IF
!
! Option to output information concerning bound-bound transitions.
! Output file has same format as that generated by CMFGEN. Only a section
! of wavelength space needs to be output.
!
	ELSE IF(XOPT .EQ. 'WRTRANS')THEN
!
	  DEFAULT=WR_STRING(LAM_ST)
	  CALL USR_OPTION(LAM_ST,'LAMST',DEFAULT,FREQ_INPUT)
	  LAM_EN=LAM_ST*1.005_LDP
	  DEFAULT=WR_STRING(LAM_EN)
	  CALL USR_OPTION(LAM_EN,'LAMEN',DEFAULT,FREQ_INPUT)
	  IF(KEV_INPUT)THEN
	    NU_ST=LAM_ST*KEV_TO_HZ
	    NU_EN=LAM_EN*KEV_TO_HZ
	  ELSE IF(ANG_INPUT)THEN
	    NU_ST=ANG_TO_HZ/LAM_ST
	    NU_EN=ANG_TO_HZ/LAM_EN
	  ELSE
	    NU_ST=LAM_ST
	    NU_EN=LAM_EN
	  END IF
!
	  NL=GET_INDX_DP(NU_ST,VEC_FREQ,N_LINE_FREQ)
	  NUP=GET_INDX_DP(NU_EN,VEC_FREQ,N_LINE_FREQ)
!	
	  I=160	!Record length - allow for long names
	  CALL GEN_ASCI_OPEN(LU_OUT,'TRANS_INFO','UNKNOWN',' ','WRITE',I,IOS)
	  WRITE(LU_OUT,*)
	1     '     I  NL_F  NUP_F      Nu',
	1     '       Lam(A)      gf      /\V(km/s)    Transition'
	  IF(NUP-NL .LT. 50)WRITE(T_OUT,*)
	1     '     I  NL_F  NUP_F      Nu',
	1     '       Lam(A)      gf      /\V(km/s)    Transition'
	  DO ML=NL,NUP
	    T1=LAMVACAIR(VEC_FREQ(ML))
	    T2=2.998E+05_LDP*(VEC_FREQ(MAX(1,ML-1))-VEC_FREQ(ML))/VEC_FREQ(ML)
	    IF(T2 .GT. 2.998E+05_LDP)T2=2.998E+05_LDP
	    ID=VEC_ION_INDX(ML)
	    T3=VEC_OSCIL(ML)*ATM(ID)%GXzV_F(VEC_MNL_F(ML))	!gf
	    IF(T1 .LT. 1.0E+04_LDP)THEN
	      WRITE(LU_OUT,
	1      '(1X,I6,2I6,F10.6,2X,F10.3,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	      IF(NUP-NL .LT. 50)WRITE(T_OUT,
	1      '(1X,I6,2I6,F10.6,2X,F10.3,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	    ELSE
	      WRITE(LU_OUT,
	1      '(1X,I6,2(1X,I6),F10.6,1X,ES11.4,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	      IF(NUP-NL .LT. 50)WRITE(T_OUT,
	1      '(1X,I6,2(1X,I6),F10.6,1X,ES11.4,ES10.2,1X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T3,T2,TRIM(VEC_TRANS_NAME(ML))
	    END IF
	  END DO
	CLOSE(UNIT=LU_OUT)

	ELSE IF(XOPT .EQ. 'WRRTK')THEN
	  DO I=1,ND
	    T1=1.0E-10_LDP*ROSS_MEAN(I)/MASS_DENSITY(I)
	    WRITE(25,'(I3,ES15.5,6ES14.4)')I,R(I)*1.0D+10,T(I)*1.0D+04,
	1          MASS_DENSITY(I),T1,
	1          LOG10(T1),LOG10(T(I))+4.0D0,LOG10(MASS_DENSITY(I)/T(I)**3)+6.0D0
	  END DO
!
	ELSE IF(XOPT .EQ. 'WRRVSIG')THEN
	  CALL USR_OPTION(ELEC,'DCF','F','Include the density and clumping factor in file?')
	  CALL GEN_ASCI_OPEN(LU_OUT,'NEW_RVSIG','UNKNOWN',' ','WRITE',I,IOS)
	  WRITE(LU_OUT,'(A)')'!'
	  IF(ELEC)THEN
	    WRITE(LU_OUT,'(A,10X,A,9X,10X,A,12X,A,2(7X,A),3X,A)')'!','R','V(km/s)','Sigma',
	1                   '    Density','Clump. Fac.','Depth'
	  ELSE
	    WRITE(LU_OUT,'(A,7X,A,9X,10X,A,12X,A,3X,A)')'!','R','V(km/s)','Sigma','Depth'
	  END IF
	  WRITE(LU_OUT,'(A)')'!'
	  WRITE(LU_OUT,'(A)')' '
	  WRITE(LU_OUT,'(I4,20X,A)')ND,'!Number of depth points`'
	  WRITE(LU_OUT,'(A)')' '
	  IF(ELEC)THEN
	    DO I=1,ND
	      WRITE(LU_OUT,'(F18.8,ES20.10,F17.7,2ES18.8,4X,I4)')R(I),V(I),SIGMA(I),MASS_DENSITY(I),CLUMP_FAC(I),I
	    END DO
	  ELSE
	    DO I=1,ND
	      WRITE(LU_OUT,'(F18.8,ES17.7,F17.7,4X,I4)')R(I),V(I),SIGMA(I),I
	    END DO
	  END IF
	  CLOSE(LU_OUT)
	  WRITE(T_OUT,'(A,ES17.10)')'     R(ND)=',R(ND)
	  WRITE(T_OUT,'(A,ES17.10)')'RMAX/R(ND)=',R(1)/R(ND)
	  WRITE(T_OUT,'(A)')'RVSIG data output to NEW_RVSIG'
	ELSE IF(XOPT .EQ. 'WRPLOT')THEN
	  WRITE(T_OUT,*)' OPTION `WRPLOT` NOT AVAILABLE'
	  WRITE(T_OUT,*)' Use WP option in plot package'
	ELSE IF(XOPT .EQ. 'RDPLOT')THEN
	  WRITE(T_OUT,*)' OPTION `RDPLOT` NOT AVAILABLE'
	  WRITE(T_OUT,*)' Use RP option in plot package'
!
! 
!
! Output suumary of commands, or a brief command summary.
! The files containing the descriptions need to be linked via DEFINE
! statements (VMS), or through soft links (UNIX).
!
	ELSE IF(XOPT .EQ. 'LI' .OR. XOPT .EQ. 'LIST' .OR.
	1       XOPT .EQ. 'HE' .OR. XOPT .EQ. 'HELP')THEN
	  IF(XOPT(1:2) .EQ. 'LI')THEN
	    CALL GEN_ASCI_OPEN(LU_IN,'MAINGEN_OPT_DESC','OLD',' ','READ',
	1                  IZERO,IOS)
	  ELSE
	    CALL GEN_ASCI_OPEN(LU_IN,'MAINGEN_OPTIONS','OLD',' ','READ',
	1                  IZERO,IOS)
	  END IF
!
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening description file'
	    GOTO 1
	  END IF
	  READ(LU_IN,*)I,K			!For page formating (I=22,K=12)
	  DO WHILE(1.EQ. 1)
	    DO J=1,I
	      READ(LU_IN,'(A)',END=700)STRING
	      L=LEN_TRIM(STRING)
	      IF(L .EQ. 0)THEN
	        WRITE(T_OUT,'(1X)')
	      ELSE
	        WRITE(T_OUT,'(1X,A)')STRING(1:L)
	      END IF
	    END DO
	    READ(T_IN,'(A)')STRING
	    IF(STRING(1:1) .EQ. 'E' .OR. STRING(1:1) .EQ. 'e' .OR.
	1       STRING(1:1) .EQ. 'Q' .OR.
	1       STRING(1:1) .EQ. 'q')GOTO 700		!Exit from listing.
	    I=K
	  END DO
700	  CONTINUE
	  CLOSE(UNIT=LU_IN)
!
	ELSE IF(XOPT .EQ. 'EX' .OR. XOPT .EQ. 'EXIT') THEN
	  STOP
	ELSE IF(X(1:4) .EQ. 'BOX=') THEN
	  CALL WR_BOX_FILE(MAIN_OPT_STR)
	ELSE
	  PRINT*,'OPTION REQUESTED DOES NOT EXIST'
	END IF
!
 1	  GO TO 3
!
! Formats used for writing out headers, net rates tec.
!
40001	  FORMAT(/,1X,A70)
40002	  FORMAT(1X,'NL=',I3,5X,'NUP=',I3)
40003	  FORMAT(3X,1P5E16.5)
!
	  END
