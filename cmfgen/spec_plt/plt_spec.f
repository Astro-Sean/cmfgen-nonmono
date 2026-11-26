C
C General routine for plotting and comparing observational and model
C spectra. Model input spectra are read from OBSLFLUX.
C
C Various options are available to redden and normalize the model spectra.
C Several different units can be used for the X and Y axes.
C
	PROGRAM PLT_SPEC
	USE SET_KIND_MODULE
C
	USE FILT_PASS_BAND
C
C Interface routines for IO routines.
C
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE GEN_IN_INTERFACE
C
	IMPLICIT NONE
C
C Altered 14-Oct-2023 : Disablled FFT option
C Altered S4-Sep-2023 : Modifiex RXY option so as to take logs if logy is set -- LONG ver -- 25-Oct-2023)
C Altered 22-Jan-2023 : MKGL option added (make Gaussian line).
C Altered 19-Apr-2020 : Added CLIP option to RD_CONT.
C Altered 16-May-2019 : NU, OBSF etc are now allocatable arrays.
C Altered 18-Jan-2015 : Added REM_BP (automatic removal of bad pixel).
C Altered 04-Nov-2015 : Fixed SMC redenning law in optical. Previous formula only
C                         valid for UV. Joined UV smoothly (at 2950A) to CCM law with R=2.74.
C Altered 15-Mar-2011 : SMC reddening law added (done by Kathryn Neugent).
C                       RED option installed -- crude method to get redenning.
C                       Crude procedure to remove cosmic rays (or spikes) form observed
C                         data implemented with rd_obs option.
C                       Fixed comvolution option to make more transparent.
C Altered 17-Jun-1996 : Bug fixed with Wavelength normalization for BB option.
C Altered 16-mar-1997 : Cleaned: USR_OPTION installation finalized.
C Altered 26-Nov-1996 : Norm option fixed so that entire MOD spectrum is
C                        plotted.
C
C
C Determines largest single plot that can read in.
C
	INTEGER, PARAMETER :: NCF_MAX=30000000
C
C Used to indicate number of data points in BB spectrum
c
	INTEGER, PARAMETER :: NBB=2000
C
	INTEGER NCF		!Number of data points in default data
	INTEGER NCF_MOD		!Used when plotting another model data set
c
	REAL(KIND=LDP), ALLOCATABLE :: NU(:)
	REAL(KIND=LDP), ALLOCATABLE :: OBSF(:)
	REAL(KIND=LDP), ALLOCATABLE :: FQW(:)
	REAL(KIND=LDP), ALLOCATABLE :: AL_D_EBmV(:)
C
	INTEGER NCF_CONT
	REAL(KIND=LDP), ALLOCATABLE :: NU_CONT(:)
	REAL(KIND=LDP), ALLOCATABLE :: OBSF_CONT(:)
!
	INTEGER NOBS
	REAL(KIND=LDP), ALLOCATABLE ::  NU_OBS(:)
	REAL(KIND=LDP), ALLOCATABLE ::  OBSF_OBS(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: XVEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: YVEC(:)
C
C Indicates which columns the observatoinal data is in.
C
	INTEGER OBS_COLS(2)
C
C Vectors for passing data to plot package via calls to CURVE.
C
	REAL*4, ALLOCATABLE ::  XV(:)
	REAL*4, ALLOCATABLE ::  YV(:)
	REAL*4, ALLOCATABLE ::  ZV(:)
	REAL*4 SP_VAL
C
	CHARACTER*80 NAME		!Default title for plot
	CHARACTER*80 XAXIS,XAXSAV	!Label for Absisca
	CHARACTER*80 YAXIS		!Label for Ordinate
	CHARACTER*80 TITLE
C
	REAL(KIND=LDP) LUM
	REAL(KIND=LDP) TOT_LUM
	REAL(KIND=LDP) N_PHOT
	REAL(KIND=LDP) ZAN_FREQ(10)
C
	REAL(KIND=LDP) ANG_TO_HZ
	REAL(KIND=LDP) KEV_TO_HZ
	REAL(KIND=LDP) NORM_WAVE
	REAL(KIND=LDP) NORM_FREQ
	REAL(KIND=LDP) NORM_FLUX
	REAL(KIND=LDP) DNU
	REAL(KIND=LDP) BB_FLUX
	REAL(KIND=LDP) SCALE_FAC
	REAL(KIND=LDP) XFAC
	REAL(KIND=LDP) ADD_FAC
	REAL(KIND=LDP) LAMC
	REAL(KIND=LDP) RAD_VEL			!Radial velcity in km/s
	REAL(KIND=LDP) C_CMS
	REAL(KIND=LDP) WT(30)			!Used when smoothing observational data.
	REAL(KIND=LDP) LAM_RANGE(2)
C
	LOGICAL NON_MONOTONIC
	LOGICAL SMOOTH			!Smooth observational data?
	LOGICAL CLEAN			!Remove IUE bad pixels?
	LOGICAL REMOVE_BAD_PIX
	LOGICAL CLN_CR			!Remove cosmic-ray spikes
	LOGICAL TREAT_AS_MOD		
	LOGICAL READ_OBS
	LOGICAL AIR_LAM
	LOGICAL CLIP_FEATURES
	LOGICAL FIRST
C
	LOGICAL WR_PLT,OVER,ABS_VALUE
C
	LOGICAL LIN_INT
	LOGICAL UNEQUAL
	LOGICAL LOG_X,LOG_Y
	CHARACTER*10 Y_PLT_OPT,X_UNIT
	CHARACTER*80 IS_FILE
	CHARACTER*80 DIRECTORY
	CHARACTER*80 XKEY,YKEY
	CHARACTER*200 FILENAME
!
! Variable for applying interstellar absorption to model spectrum.
!
	LOGICAL HI_ABS                   ! correct for HI absorption
	LOGICAL H2_ABS                   ! correct for H2 absorption
	REAL(KIND=LDP) T_IN_K                    ! temp in K of intersteallar H&HII
	REAL(KIND=LDP) V_TURB                    ! turbulent velocity of "     "
	REAL(KIND=LDP) V_R                       ! radial v of star w.r.t. ISM
	REAL(KIND=LDP) LOG_NTOT                  ! log of H column density
	REAL(KIND=LDP) LOG_H2_NTOT               ! log of H2 column density
	LOGICAL FFT_CONVOLVE             ! use FFT method to convolve data
	REAL(KIND=LDP) INST_RES                  ! desired instrument resolution (dl)
	REAL(KIND=LDP) MIN_RES_KMS               ! minimum resolution for model data
	REAL(KIND=LDP) NUM_RES                   ! number of res. elements to cons.
	REAL(KIND=LDP) RESOLUTION                ! desired resolution (R=l/dl)
	REAL(KIND=LDP) WAVE_MAX                  ! max wavelength to convolve over
	REAL(KIND=LDP) WAVE_MIN                  ! min wavelength to convolve over
	REAL(KIND=LDP) VSINI
	REAL(KIND=LDP) EPSILON
C
C Miscellaneous variables.
C
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,K,L,ML,MLST
	INTEGER IST,IEND
	INTEGER CNT
	INTEGER NHAN
!
	REAL(KIND=LDP) LAM_CENT
	REAL(KIND=LDP) HEIGHT
	REAL(KIND=LDP) SIGMA
!
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) SUM
	REAL(KIND=LDP) TEMP
	REAL(KIND=LDP) TMP_FREQ
C
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_IN=5		!For file I/O
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LU_IN=10	!For file I/O
	INTEGER, PARAMETER :: LU_OUT=11
C
	REAL(KIND=LDP), PARAMETER :: EDGE_HYD=3.28808662499619_LDP
	REAL(KIND=LDP), PARAMETER :: EDGE_HEI=5.94520701882481_LDP
	REAL(KIND=LDP), PARAMETER :: EDGE_HE2=13.1581564178623_LDP
C
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
	LOGICAL ADD_NOISE
	INTEGER IRAN
	REAL(KIND=LDP) R_SEED
	REAL(KIND=LDP) COUNTS
	REAL(KIND=LDP) POIDEV
C
C USR_OPTION variables
C
	CHARACTER MAIN_OPT_STR*80	!Used for input of the main otion
	CHARACTER X*10			!Used for the idividual option
	CHARACTER STRING*80
	CHARACTER*120 DEFAULT
	CHARACTER*120 DESCRIPTION
C
	REAL(KIND=LDP) SPEED_OF_LIGHT,FAC,LAM_VAC
	INTEGER GET_INDX_SP
	LOGICAL EQUAL
	CHARACTER*30 UC
	CHARACTER*30 FILTER_SET
	EXTERNAL SPEED_OF_LIGHT,EQUAL,FAC,UC,LAM_VAC,GET_INDX_DP
C
C 
C Filter and extinction data:
C
	REAL(KIND=LDP) NORM_LUM
	REAL(KIND=LDP) DIST
	REAL(KIND=LDP) R_EXT
	REAL(KIND=LDP) EBMV_GAL
	REAL(KIND=LDP) EBMV_LMC
	REAL(KIND=LDP) EBMV_CCM
	REAL(KIND=LDP) EBMV_SMC         !KN - SMC law
	REAL(KIND=LDP) RAX,RBX
	REAL(KIND=LDP) FILTLAM(8),FILTZP(8),FLAM(24),ZERO_POINT
	CHARACTER*1 FILT(8)
	REAL(KIND=LDP) C1,C2,C3,C4,D,F  !KN For SMC
!
	REAL(KIND=LDP) RESPONSE1
	REAL(KIND=LDP) RESPONSE2
	REAL(KIND=LDP) FILT_INT_BEG
	INTEGER IF
	INTEGER ML_ST
	INTEGER ML_END

C
	DATA FILT/' ','u','b','v','r','J','H','K'/
	DATA FILTLAM/0.3644_LDP,0.3650_LDP,0.427_LDP,0.5160_LDP,0.60_LDP,1.25_LDP,1.65_LDP,2.20_LDP/
	DATA FLAM/0.1_LDP,.125_LDP,0.15_LDP,0.175_LDP,0.20_LDP,0.225_LDP,0.25_LDP,0.275_LDP,0.300_LDP,0.325_LDP
	1, 0.350_LDP,0.4_LDP,0.5_LDP,0.6_LDP,0.7_LDP,0.8_LDP,0.9_LDP,1.25_LDP,1.65_LDP,2.2_LDP,3.6_LDP,4.8_LDP,10.0_LDP,20.0_LDP/
C
C Zero point is Based on Flagstaff calibration which gives F_nu(Alpha Lyrae)=3560Jy. Note the
C V magnitude of (Alpha Lyrae) is 0.03.
C
	DATA ZERO_POINT/3560/
	DATA FILTZP/3560.0_LDP,3560.0_LDP,3560.0_LDP,3560.0_LDP,3560.0_LDP,1564.0_LDP,1008.0_LDP,628.0_LDP/
C
C 
C
	IOS=0
	ALLOCATE (NU(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (OBSF(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (FQW(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (AL_D_EBmV(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (NU_CONT(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (OBSF_CONT(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (NU_OBS(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (OBSF_OBS(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (XV(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (YV(NCF_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ZV(NCF_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error allocating OBSF, FQW etc in PLT_SPEC'
	  WRITE(6,*)'Error is ',IOS
	  STOP
	END IF
C
C Set constants.
C
	CHIBF=2.815E-06_LDP
	CHIFF=3.69E-29_LDP
	HDKT=4.7994145_LDP
	TWOHCSQ=0.0147452575_LDP
	OPLIN=2.6540081E+08_LDP
	EMLIN=5.27296E-03_LDP
C
	C_CMS=SPEED_OF_LIGHT()
C
C Set defaults.
C
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	NAME=' '
	TITLE=' '
!	PI=FUN_PI()
C
	LOG_X=.FALSE.
	LOG_Y=.FALSE.
	X_UNIT='ANG'
	Y_PLT_OPT='FNU'
!
	NCF=0
	NOBS=0
C
C Conversion factor from Kev to units of 10^15 Hz.
C Conversion factor from Angstroms to units of 10^15 Hz.
C
	KEV_TO_HZ=0.241838E+03_LDP
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0E-07_LDP  	!10^8/10^15
C
C  Read in default model.
C
	NCF=0			!In case no model read
	FILENAME='OBSFLUX'
	CALL GEN_IN(FILENAME,'Model file')
	CALL RD_MOD(NU,OBSF,NCF_MAX,NCF,FILENAME,IOS)
C
C 
C
C This message will only be printed once
C
	WRITE(T_OUT,*)
	WRITE(T_OUT,"(8X,A)")'(default is to write file '//
	1    'main_option.sve)'
	WRITE(T_OUT,"(8X,A)")'(append sve=filename to '//
	1    'write a new .sve file)'
	WRITE(T_OUT,"(8X,A)")'(box=filename to '//
	1    'write a .box file containing several .sve files)'
	WRITE(T_OUT,"(8X,A)")'(.filename to read .sve file)'
	WRITE(T_OUT,"(8X,A)")'(#filename to read .box file)'
	WRITE(T_OUT,*)
C
C This call resets the .sve algorithm.  Specifically it sets the next
C input answer to be a main option, and all subsequent inputs to be
C sub-options.
C
 3	CALL SVE_FILE('RESET')
C
	MAIN_OPT_STR='  '
	DEFAULT='GR'
	DESCRIPTION=' '					!Obvious main option
	CALL USR_OPTION(MAIN_OPT_STR,'OPTION',DEFAULT,DESCRIPTION)
	IF(MAIN_OPT_STR .EQ. ' ')GOTO 3
C
C   If the main option begins with a '.', a previously
C   written .sve file is read.
C
C   If the main option begins with a '#', a previously
C   written .box file is read.
C
C   If sve= is apended to the end of this main option, a new .sve file
C   is opened with the given name and the main option and all subsequent
C   sub-options are written to this file.
C
C   If box= is input then a .box file is created, which contains the name
C   of several .sve files to process.
C
C   If only a main option is given, the option and subsequent sub-options
C   are saved in a file called 'main option.sve'.  All following main
C   options are saved in separate files.
C
	X=UC(MAIN_OPT_STR)
	I=INDEX(X,'(')
	IF(I .NE. 0)X=X(1:I-1)		!Remove line variables
C
C 
C
	IF(X(1:3) .EQ. 'TIT')THEN
	  CALL USR_OPTION(NAME,'Title',' ',' ')
C
C Set X-Ais plotting options.
C
	ELSE IF(X(1:2) .EQ.'LX' .OR. X(1:4) .EQ. 'LOGX' .OR.
	1                            X(1:4) .EQ. 'LINX')THEN
	  LOG_X=.NOT. LOG_X
	  IF(LOG_X)WRITE(T_OUT,*)'Now using Logarithmic X axis'
	  IF(.NOT. LOG_X)WRITE(T_OUT,*)'Now using Linear X axis'
	ELSE IF(X(1:2) .EQ.'XU' .OR. X(1:6) .EQ. 'XUNITS')THEN
	  CALL USR_OPTION(X_UNIT,'X_UNIT','Ang',
	1                  'Ang, AA[Air Ang], nm, um, eV, keV, Hz, Mm/s, km/s')
	  CALL SET_CASE_UP(X_UNIT,IZERO,IZERO)
	  IF(X_UNIT .NE. 'ANG' .AND.
	1        X_UNIT .NE. 'AA' .AND.
	1        X_UNIT .NE. 'UM' .AND.
	1        X_UNIT .NE. 'NM' .AND.
	1        X_UNIT .NE. 'EV' .AND.
	1        X_UNIT .NE. 'KEV' .AND.
	1        X_UNIT .NE. 'HZ' .AND.
	1        X_UNIT .NE. 'MM/S' .AND.
	1        X_UNIT .NE. 'KM/S')THEN
	     WRITE(T_OUT,*)'Invalid X unit: Try again'
	   END IF
C
C NB: We offer the option to use the central frequency to avoid
C air/vacuum confusions. Model data is in vacuum wavelengths, which
C we use in plotting at all wavelngths.
C
	   IF(X_UNIT .EQ. 'MM/S' .OR. X_UNIT .EQ. 'KM/S')THEN
	     CALL USR_OPTION(LAMC,'LAMC','0.0',
	1             'Central Lambda(Ang) [-ve for frequency (10^15 Hz)]')
	     IF(LAMC .LT. 0)THEN
	       LAMC=1.0E-07_LDP*C_CMS/ABS(LAMC)
	     ELSE
	       IF(LAMC .GT. 2000)THEN
                 CALL USR_OPTION(AIR_LAM,'AIR','T',
	1                'Air wavelength [only for Lam > 2000A]?')
	       ELSE
	         AIR_LAM=.FALSE.
	       END IF
	     END IF
	     IF(AIR_LAM)LAMC=LAM_VAC(LAMC)
	   END IF
C
C Set Y axis plotting options.
C
	ELSE IF(X(1:2) .EQ. 'LY' .OR. X(1:4) .EQ. 'LOGY' .OR.
	1                             X(1:4) .EQ. 'LINY')THEN
	  LOG_Y=.NOT. LOG_Y
	  IF(LOG_Y)WRITE(T_OUT,*)'Now using Logarithmic Y axis'
	  IF(.NOT. LOG_Y)WRITE(T_OUT,*)'Now using Linear Y axis'
	ELSE IF(X(1:2) .EQ.'YU' .OR. X(1:6) .EQ. 'YUNITS')THEN
	  CALL USR_OPTION(Y_PLT_OPT,'Y_UNIT',' ',
	1          'FNU, NU_FNU, FLAM, LAM_FLAM')
	  CALL SET_CASE_UP(Y_PLT_OPT,IZERO,IZERO)
	  IF(Y_PLT_OPT .NE. 'FNU' .AND.
	1        Y_PLT_OPT .NE. 'NU_FNU' .AND.
	1        Y_PLT_OPT .NE. 'LAM_FLAM' .AND.
	1        Y_PLT_OPT .NE. 'FLAM')THEN
	     WRITE(T_OUT,*)'Invalid Y Plot option: Try again'
	  END IF
C
C 
C
	ELSE IF(X(1:6) .EQ. 'RD_MOD')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Model file')
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  CALL USR_HIDDEN(XFAC,'XFAC','1.0D0',' ')
	  CALL USR_HIDDEN(TITLE,'TIT',' ',' ')
	  CALL USR_HIDDEN(RAD_VEL,'RAD_VEL','0.0D0','Radial velocity(km/s) of star')
	  IF(XFAC .NE. 1.0_LDP .AND. RAD_VEL .NE. 0.0_LDP)THEN
	    WRITE(6,*)'Only one of XFAC and RAD_VEL can be changed from their default values of 1 and 0'
	    GOTO 1
	  ELSE IF(RAD_VEL .NE. 0.0_LDP)THEN
	    XFAC=(1.0_LDP+1.0E+05_LDP*RAD_VEL/C_CMS)
	  END IF
	  IF(OVER)THEN
C
C This option allows all normal model options to be done on the data
C (e.g. redenning).
C
	    CALL RD_MOD(NU,OBSF,NCF_MAX,NCF,FILENAME,IOS)
	    IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error reading model data -- no data read'
	       GOTO 1		!Get another option
	    END IF
	    OVER=.FALSE.
!
	    IF(SCALE_FAC .NE. 1.0_LDP .OR. XFAC .NE. 1.0_LDP)THEN
	      OBSF(1:NCF)=OBSF(1:NCF)*SCALE_FAC
	      NU(1:NCF)=NU(1:NCF)*XFAC
	      WRITE(T_OUT,*)'Model has been scaled!'
	    ELSE
	      WRITE(T_OUT,*)'No scaling done with new model data'
	    END IF
!
	    WRITE(T_OUT,*)'New model data replaces old data'
	    WRITE(T_OUT,*)'No plots done with new model data'
	  ELSE
C
C This option is now similar to RD_CONT
C
	    CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_MOD,FILENAME,IOS)
	    IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error reading model data -- no data read'
	       GOTO 1		!Get another option
	    END IF
	    DO I=1,NCF_MOD
	      XV(I)=NU_CONT(I)*XFAC
	      YV(I)=OBSF_CONT(I)*SCALE_FAC
	    END DO
	    CALL CNVRT(XV,YV,NCF_MOD,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL CURVE_LAB(NCF_MOD,XV,YV,TITLE)
	  END IF
C
	ELSE IF(X(1:7) .EQ. 'RD_CONT')THEN
	  FILENAME=' '
	  CLIP_FEATURES=.FALSE.
	  CALL USR_OPTION(FILENAME,'File',' ','Continuum file')
	  CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_CONT,FILENAME,IOS)
	  CALL USR_HIDDEN(CLIP_FEATURES,'CLIP','F','Clipe out continuum features')
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  CALL USR_HIDDEN(TITLE,'TIT',' ',' ')
	  IF(IOS .NE. 0)GOTO 1		!Get another option
	  SCALE_FAC=1.0_LDP
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  IF(CLIP_FEATURES)THEN
	    WRITE(6,*)'Clipping fetaures from spectrum'
	    CALL CLIP(NU_CONT,OBSF_CONT,NCF_CONT)
	  END IF
!
	  IF(OVER)THEN
	    NCF=NCF_CONT
	    DO I=1,NCF
	      NU(I)=NU_CONT(I)
	      OBSF(I)=OBSF_CONT(I)*SCALE_FAC
	    END DO
	  ELSE
	    DO I=1,NCF_CONT
	      XV(I)=NU_CONT(I)
	      YV(I)=OBSF_CONT(I)*SCALE_FAC
	    END DO
	    CALL CNVRT(XV,YV,NCF_CONT,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	    IF(WR_PLT)THEN
	      DO I=1,NCF_CONT
	        WRITE(50,*)XV(I),YV(I)
	      END DO
	    ELSE
	      CALL CURVE_LAB(NCF_CONT,XV,YV,TITLE)
	    END IF
	  END IF
!
	ELSE IF(X(1:6) .EQ. 'RD_POL')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Polarization file')
	  CALL USR_OPTION(IST,'IREC','1','Angle (polarization record)')
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  CALL USR_HIDDEN(TITLE,'TIT',' ',' ')
	  IF(OVER)THEN
	    CALL RD_SING_POL_I(OBSF,NU,NCF,NCF_MAX,IST,FILENAME,LU_IN,IOS)
	    IF(IOS .NE. 1)GOTO 1
	  ELSE
	    CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	    CALL RD_SING_POL_I(OBSF_CONT,NU_CONT,NCF_CONT,NCF_MAX,IST,FILENAME,LU_IN,IOS)
	    IF(IOS .NE. 0)GOTO 1
	    DO I=1,NCF_CONT
	      XV(I)=NU_CONT(I)
	      YV(I)=OBSF_CONT(I)
	    END DO
	    CALL CNVRT(XV,YV,NCF_CONT,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	    IF(WR_PLT)THEN
	      DO I=1,NCF_CONT
	        WRITE(50,*)XV(I),YV(I)
	      END DO
	    ELSE
	      CALL CURVE_LAB(NCF_CONT,XV,YV,TITLE)
	    END IF
	  END IF
!
! This option simply reads in data in XY format. No conversion is done to the data.
! Comments (bginning with !) and blank lines are ignored in the data file.
!
! If logy is set, the logarithm of the y data is taken.
!
	ELSE IF(X(1:3) .EQ. 'RXY')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Model file')
	  CALL RD_XY_DATA_USR(NU_CONT,OBSF_CONT,NCF_CONT,NCF_MAX,FILENAME,LU_IN,IOS)
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  CALL USR_HIDDEN(ABS_VALUE,'ABSR','F','Use absolute values when taking log')
	  CALL USR_HIDDEN(TITLE,'TIT',' ',' ')
	  IF(OVER .AND. IOS .EQ. 0)THEN
	    NCF=NCF_CONT
	    NU(1:NCF)=NU_CONT(1:NCF)
	    OBSF(1:NCF)=OBSF_CONT(1:NCF)
	    WRITE(6,*)'As read in to PLT_SPEC buffer, X will assumed to be NU'
	  ELSE IF(IOS .EQ. 0)THEN
	    T1=MAXVAL(OBSF_CONT(1:NCF_CONT))
	    IF(T1 .GT. 1.0E+38_LDP)THEN
	      WRITE(6,*)'Data exceeds single precision range: Maximum=',T1
	      WRITE(6,*)'Necessary to scale data for plotting'
	      T1=1.0_LDP
	      CALL USR_OPTION(T1,'SCL_FAC','1.0D+40','Factor to divide data by')
	      OBSF_CONT(1:NCF_CONT)=OBSF_CONT(1:NCF_CONT)/T1
	    END IF
	    IF(LOG_Y .AND. ABS_VALUE)THEN
	      OBSF_CONT(1:NCF_CONT)=LOG10(ABS(OBSF_CONT(1:NCF_CONT)))
	    ELSE IF(LOG_Y)THEN
	      FIRST=.TRUE.
	      DO I=1,NCF_CONT
	        IF(OBSF_CONT(I) .GT. 0)THEN
	          OBSF_CONT(I)=LOG10(OBSF_CONT(I))
	        ELSE
	          OBSF_CONT(I)=-200
	          IF(FIRST)THEN
	             FIRST=.FALSE.
	             WRITE(6,*)'Warning -- zero or neagtive values encountered'
	 	     WRITE(6,*)'Log set to -200'
	          END IF
	        END IF
	      END DO
	    END IF
	    CALL DP_CURVE_LAB(NCF_CONT,NU_CONT,OBSF_CONT,TITLE)
	  END IF
!
	ELSE IF(X(1:7) .EQ. 'RROW')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Model file')
	  CALL USR_OPTION(XKEY,'XKEY',' ','Key with X data')
	  CALL USR_OPTION(YKEY,'YKEY',' ','Key with Y data')
	  CALL USR_OPTION(NCF_CONT,'N',' ','Number of data points')
	  CALL RD_ROW_DATA(NU_CONT,OBSF_CONT,NCF_CONT,XKEY,YKEY,FILENAME,LU_IN,IOS)
	  IF(IOS .EQ. 0)THEN
	    CALL DP_CURVE_LAB(NCF_CONT,NU_CONT,OBSF_CONT,TITLE)
	  END IF
C
	ELSE IF(X(1:4) .EQ. 'NORM')THEN
C
	  READ_OBS=.FALSE.
	  CALL USR_HIDDEN(READ_OBS,'RD_OBS','F',' ')
	  CALL USR_HIDDEN(TITLE,'TIT',' ',' ')
C
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Continuum file')
	  IF(READ_OBS)THEN
	    CALL USR_HIDDEN(OBS_COLS,2,2,'COLS','1,2','Columns with data')
	    CALL RD_OBS_DATA_V2(XV,YV,NCF_MAX,J,FILENAME,OBS_COLS,IOS)
	    IF(IOS .NE. 0)GOTO 1		!Get another option
	    DO I=1,J
	      NU_CONT(I)=ANG_TO_HZ/XV(I)
	      OBSF_CONT(I)=YV(I)
	    END DO
	    NCF_CONT=J
	  ELSE
	    CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_CONT,FILENAME,IOS)
	    IF(IOS .NE. 0)GOTO 1		!Get another option
	  END IF
C
	  T1=1.0E-08_LDP
	  I=1
	  UNEQUAL=.FALSE.
	  IF(NCF_CONT .NE. NCF)UNEQUAL=.TRUE.
	  DO WHILE(.NOT. UNEQUAL .AND. I .LE. NCF_CONT)
	    IF( EQUAL(NU_CONT(I),NU(I),T1) )THEN
	      XV(I)=NU_CONT(I)
	      I=I+1
	    ELSE
	      UNEQUAL=.TRUE.
	    END IF
	  END DO
	  IF(UNEQUAL)CALL USR_HIDDEN(LIN_INT,'LIN','F',
	1               'Use linear interpolation')
	  IF(UNEQUAL .AND. LIN_INT)THEN
	    L=1
	    DO I=1,NCF
  	      XV(I)=NU(I)
	      IF(NU(I) .GT. NU_CONT(1))THEN
	        YV(I)=0.0
	      ELSE IF(NU(I) .LT. NU_CONT(NCF_CONT))THEN
	        YV(I)=0.0
	      ELSE
	        DO WHILE (NU(I) .LT. NU_CONT(L+1))
	          L=L+1
	        END DO
	        T1=(NU(I)-NU_CONT(L+1))/(NU_CONT(L)-NU_CONT(L+1))
	        T2=(1.0_LDP-T1)*OBSF_CONT(L+1)+T1*OBSF_CONT(L)
	        YV(I)=0.0
	        IF(T2 .NE. 0)THEN
	          T2=OBSF(I)/T2
	          IF(T2 .LT. 1.0E+020_LDP)YV(I)=T2
	        END IF
	      END IF
	    END DO
	  ELSE IF(UNEQUAL)THEN
C
C We will use monotonic cubic interpolation. We first verify the range.
C I & J are temporary variables for the callt o MON_INTERP. I denotes the
C first element. Initially J denotes the last element, then the numer of
C elements that can be interpolated.
C
	    I=1
	    DO WHILE(NU(I) .GT. NU_CONT(1))
	      I=I+1
	    END DO
	    J=NCF
	    DO WHILE(NU(J) .LE. NU_CONT(NCF_CONT))
	      J=J-1
	    END DO
	    J=J-I+1
C
	    FQW(1:NCF)=0.0_LDP				!Temporary usage
	    CALL MON_INTERP(FQW(I),J,IONE,NU(I),J,
	1            OBSF_CONT,NCF_CONT,NU_CONT,NCF_CONT)
  	    XV(1:NCF)=NU(1:NCF)
	    DO I=1,NCF
	      IF(FQW(I) .GT. 0)THEN
	         T2=OBSF(I)/FQW(I)
	         IF(T2 .LT. 1.0E+20_LDP)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  ELSE
	    DO I=1,NCF
	      YV(I)=0
	      IF(OBSF_CONT(I) .GT. 0)THEN
	         T2=OBSF(I)/OBSF_CONT(I)
	         IF(T2 .LT. 1.0E+20_LDP)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  END IF
!
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	  CALL USR_HIDDEN(ADD_NOISE,'ADD','F','Add poisonian noise?')
	  OVER=.FALSE.
	  CALL USR_OPTION(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  IF(ADD_NOISE)THEN
	    CALL USR_OPTION(COUNTS,'CNTS','100','Counts in continuum')
	    CALL USR_OPTION(R_SEED,'R_SEED','0.2',
	1              'Random number seed (0 to 1)')
	    CALL USR_OPTION(LAM_RANGE,2,2,'LAM=','3000.0,7000.0',
	1              'Wavelength Range (A)')
	    IRAN=-R_SEED*1234567
	    DO I=1,NCF
	      T1=YV(I)*COUNTS
	      IF(T1 .GE. 1 .AND. XV(I) .GT. LAM_RANGE(1) .AND. XV(I)
	1                                       .LT. LAM_RANGE(2))
	1           YV(I)=POIDEV(T1,IRAN)/COUNTS
	    END DO
	  END IF
	  IF(OVER)THEN
	     OBSF(1:NCF)=YV(1:NCF)
	  ELSE IF(WR_PLT)THEN
	    DO I=1,NCF
	      WRITE(50,*)XV(I),YV(I)
	    END DO
	  ELSE
	    CALL CURVE_LAB(NCF,XV,YV,TITLE)
	  END IF
	  YAXIS='F\d\gn\u/F\dc\u'
C
C 
C
	ELSE IF(X(1:6) .EQ. 'RD_EW')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','EW file')
C
	    CALL RD_EW(XV,YV,NCF_MAX,J,FILENAME,IOS)
	    IF(IOS .NE. 0)GOTO 1		!Get another option
	    CALL CURVE_LAB(J,XV,YV,TITLE)
C
C 
C
	ELSE IF(X(1:6) .EQ. 'RD_OBS')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ',' ')
C
	  SCALE_FAC=1.0_LDP
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  ADD_FAC=0.0_LDP
	  CALL USR_HIDDEN(ADD_FAC,'ADD','0.0D0',' ')
	  CALL USR_HIDDEN(TITLE,'TIT',' ',' ')
C
	  RAD_VEL=0.0_LDP
	  CALL USR_HIDDEN(RAD_VEL,'RAD_VEL','0.0D0','Radial velocity of star (+ve if away)')
C
	  CLEAN=.FALSE.
	  CALL USR_HIDDEN(CLEAN,'CLEAN','F',' ')
C
	  REMOVE_BAD_PIX=.FALSE.
	  CALL USR_HIDDEN(REMOVE_BAD_PIX,'REM_BP','F','Remove bad pixels? ')
C
	  CLN_CR=.FALSE.
	  CALL USR_HIDDEN(CLN_CR,'CLN_CR','F','Remove cosmic ray spikes? ')
C
	  SMOOTH=.FALSE.
	  CALL USR_HIDDEN(SMOOTH,'SMOOTH','F',' ')
C
C Option allows the observational data to be treated as though it were
C a model. Thus can reden the data etc. No plot is done, and the model
C is overwritten.
C
	  TREAT_AS_MOD=.FALSE.
	  CALL USR_HIDDEN(TREAT_AS_MOD,'OVER','F',' ')
C
	  IF(SMOOTH)THEN
	    NHAN=5
	    CALL USR_OPTION(NHAN,'HAN','5','Number of points for HAN [ODD]')
	    NHAN=2*(NHAN/2)+1	!Ensures odd.
	  END IF
C
	  CALL USR_HIDDEN(OBS_COLS,2,2,'COLS','1,2','Columns with data')
	  CALL RD_OBS_DATA_V2(XV,YV,NCF_MAX,J,FILENAME,OBS_COLS,IOS)
	  IF(IOS .NE. 0)GOTO 1		!Get another option
!
! Convert from Ang (Vacuum) to Hz.
!
	  DO I=1,J
	    XV(I)=ANG_TO_HZ/XV(I)
	    YV(I)=YV(I)*SCALE_FAC+ADD_FAC
	  END DO
C
	  IF(RAD_VEL .NE. 0)THEN
	    DO I=1,J
	     XV(I)=XV(I)*(1.0_LDP+1.0E+05_LDP*RAD_VEL/C_CMS)
	    END DO
	  END IF
C
C Procdure to remove single pizels that are zero due quirks with IUE.
C
	  IF(CLEAN)THEN
	    DO I=1,J
	      ZV(I)=YV(I)
	    END DO
	    DO I=2,J-1
	      IF(YV(I) .EQ. 0)THEN
	        YV(I)=0.5_LDP*(ZV(I-1)+ZV(I+1))
	      ELSE IF(YV(I) .LT. -1.0E+10_LDP)THEN
	        YV(I)=0.0_LDP
	      END IF
	    END DO
	  END IF
!
	  IF(REMOVE_BAD_PIX)THEN
	    DO L=3,J-50,90
	      T1=0.0_LDP; T2=0.0_LDP; T3=0.0_LDP
	      DO K=L,MIN(L+99,J)
	        T1=T1+YV(K)
	        T2=T2+YV(K)*YV(K)
	        T3=T3+1
	      END DO
	      T1=T1/T3
	      T2=SQRT( (T2-T3*T1*T1)/(T3-1) )
	      DO K=L+1,MIN(L+98,J-1)
	        IF( ABS(YV(K)-T1) .GT. 5.0_LDP*T2 .AND.
	1           ABS(YV(K-1)-T1) .LT. 3.0_LDP*T2 .AND.
	1           ABS(YV(K+1)-T1) .LT. 3.0_LDP*T2)THEN
	           YV(K)=0.5_LDP*(YV(K-1)+YV(K+1))
	        END IF
	      END DO
	    END DO
	  END IF
!
	  IF(CLN_CR)THEN
	    ZV(1:J)=YV(1:J)
	    DO L=3,J-10
	      T1=MAXVAL(YV(L:L+9))
	      DO I=L,MIN(J-2,L+9)
	         IF(YV(I) .EQ. T1)THEN
	           K=I
	           EXIT
	         END IF
	      END DO
	      T1=0.0_LDP; T2=0.0_LDP; CNT=0
	      DO I=MAX(1,K-10),MIN(K+10,J)
	        IF(I .LT. K-1 .OR. I .GT. K+1)THEN
	          T1=T1+YV(I)
	          T2=T2+YV(I)*YV(I)
	          CNT=CNT+1
	        END IF
	      END DO
	      IF(CNT .GE. 4)THEN
	        T1=T1/CNT
	        T2=T2-T1*T1*CNT
	        IF(T2 .GT. 0.0_LDP)THEN
	          T2=SQRT(T2/(CNT-1.0_LDP))
	        ELSE
	          T1=YV(K)
	          T2=1.0E+20_LDP
	        END IF
	      ELSE
	        T1=YV(K)
	        T2=1.0E+20_LDP
	      END IF
	      IF(YV(K) .GT. T1+4.0_LDP*T2)THEN
	        YV(K)=YV(K-1)
	        IF(YV(K-1) .GT. T1+4.0_LDP*T2)YV(K)=YV(K-2)
		T3=YV(K+1)
	        IF(T3 .GT. T1+4.0_LDP*T2)T3=YV(K+2)
	        YV(K)=0.5_LDP*(YV(K)+T3)
	      END IF
	    END DO
	  END IF
C
C We check whether the X axis is monotonic. If not, we smooth each section
C separately. Designed for non-merged overlapping ECHELLE orders.
C
	  NON_MONOTONIC=.FALSE.
	  IF(SMOOTH)THEN
	    T2=XV(2)-XV(1)
	    DO I=1,J-1
	      IF( (XV(I)-XV(I+1))*T2 .LT. 0)THEN
	        NON_MONOTONIC=.TRUE.
	        EXIT
	      END IF
	    END DO
	  END IF
C
	  K=NHAN
	  IF(SMOOTH .AND. NON_MONOTONIC)THEN
	    ZV(1:J)=YV(1:J)
	    DO I=1,K
	      WT(I)=FAC(K-1)/FAC(I-1)/FAC(K-I)
	    END DO
	    IST=1
	    IEND=0
	    T2=XV(2)-XV(1)
	    DO WHILE(IEND .LT. J)
	      IEND=J
	      DO I=IST,J-1
	        IF( (XV(I+1)-XV(I))*T2 .LT. 0)THEN
	          IEND=I
	          EXIT
	        END IF
	      END DO
	      WRITE(6,*)IST,IEND
	      DO I=IST,IEND
	        T1=0.0_LDP
	        YV(I)=0.0_LDP
	        DO L=MAX(IST,I-K/2),MIN(IEND,I+k/2)
                  ML=L-I+K/2+1
	          T1=T1+WT(ML)
	          YV(I)=YV(I)+ZV(L)*WT(ML)
	        END DO
	        YV(I)=YV(I)/T1
	      END DO
	      IST=IEND+1
	    END DO
C
	  ELSE IF(SMOOTH)THEN
	    DO I=1,J
	      ZV(I)=YV(I)
	    END DO
	    DO I=1,K
	      WT(I)=FAC(K-1)/FAC(I-1)/FAC(K-I)
	    END DO
	    DO I=1,J
	      T1=0.0_LDP
	      YV(I)=0.0_LDP
	      DO L=MAX(1,I-K/2),MIN(J,I+k/2)
                ML=L-I+K/2+1
	        T1=T1+WT(ML)
	        YV(I)=YV(I)+ZV(L)*WT(ML)
	      END DO
	      YV(I)=YV(I)/T1
	    END DO	
	  END IF

C
	  IF(TREAT_AS_MOD)THEN
	    NCF=J
	    DO I=1,NCF
	      NU(I)=XV(I)
	      OBSF(I)=YV(I)
	    END DO
	    WRITE(T_OUT,*)'Observational data replaces model data'
	  ELSE
	    NOBS=J
	    NU_OBS(1:NOBS)=XV(1:NOBS)
	    OBSF_OBS(1:NOBS)=YV(1:NOBS)
	    CALL CNVRT(XV,YV,J,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	    IF(WR_PLT)THEN
	      DO I=1,J
	        WRITE(50,*)XV(I),YV(I)
	      END DO
	    ELSE
	      CALL CURVE_LAB(J,XV,YV,TITLE)
	    END IF
	  END IF
!
! 
!
	ELSE IF(X(1:8) .EQ. 'RD_MONTE')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ',' ')
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  CALL USR_HIDDEN(XFAC,'XFAC','1.0D0',' ')
	  CALL USR_HIDDEN(RAD_VEL,'RAD_VEL','0.0D0','Radial velocity(km/s) of star')
	  IF(XFAC .NE. 1.0_LDP .AND. RAD_VEL .NE. 0.0_LDP)THEN
	    WRITE(6,*)'Only one of XFAC and RAD_VEL can be changed from their default values of 1 and 0'
	    GOTO 1
	  ELSE IF(RAD_VEL .NE. 0.0_LDP)THEN
	    XFAC=(1.0_LDP+1.0E+05_LDP*RAD_VEL/C_CMS)
	  END IF
!
	 IF(OVER)THEN
!
! This option allows all normal model options to be done on the data
! (e.g. redenning).
!
	    CALL RD_MONTE_LINE(NU,OBSF,NCF_MAX,NCF,FILENAME,IOS)
	    IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error reading model data -- no data read'
	       GOTO 1		!Get another option
	    END IF
	    OVER=.FALSE.
!
	    IF(SCALE_FAC .NE. 1.0_LDP .OR. XFAC .NE. 1.0_LDP)THEN
	      OBSF(1:NCF)=OBSF(1:NCF)*SCALE_FAC
	      NU(1:NCF)=NU(1:NCF)*XFAC
	      WRITE(T_OUT,*)'Model has been scaled!'
	    ELSE
	      WRITE(T_OUT,*)'No scaling done with new model data'
	    END IF
!
	    WRITE(T_OUT,*)'New model data replaces old data'
	    WRITE(T_OUT,*)'No plots done with new model data'
	  ELSE
!
! This option is now similar to RD_CONT
!
	    CALL RD_MONTE_LINE(NU_CONT,OBSF_CONT,NCF_MAX,NCF_MOD,FILENAME,IOS)
	    IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error reading model data -- no data read'
	       GOTO 1		!Get another option
	    END IF
	    DO I=1,NCF_MOD
	      XV(I)=NU_CONT(I)*XFAC
	      YV(I)=OBSF_CONT(I)*SCALE_FAC
	    END DO
	    CALL CNVRT(XV,YV,NCF_MOD,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL CURVE_LAB(NCF_MOD,XV,YV,TITLE)
	  END IF
!
	ELSE IF(X(1:5) .EQ. 'ISABS') THEN
	  CALL USR_OPTION(T_IN_K,'T_IN_K','100d0','Temp. in Kelvin')
	  CALL USR_OPTION(V_TURB,'V_TURB','10d0','Turbulent Velocity (km/s)')
	  CALL USR_OPTION(LOG_NTOT,'LOG_NTOT','20d0','Log of H column density')
	  CALL USR_OPTION(LOG_H2_NTOT,'LOG_H2_NTOT','20d0','Log H2 column dens')
!
	  CALL USR_HIDDEN(HI_ABS,'HI_ABS','T','Correct for HI absorption')
	  CALL USR_HIDDEN(H2_ABS,'H2_ABS','T','Correct for HII absorption')
	  CALL USR_HIDDEN(V_R,'V_R','0.0D0','Radial Velocity (km/s)')
	  CALL USR_HIDDEN(WAVE_MIN,'WAVE_MIN','900d0','Minimum Wavelength')
	  CALL USR_HIDDEN(WAVE_MAX,'WAVE_MAX','3000d0','Maximum Wavelength')
	  CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','2.0D0','Minimum Model Resolution')
	  CALL USR_HIDDEN(IS_FILE,'IS_FILE','IS_LINE_LIST','File wth list and strengths of IS lines')
!
	  CALL UVABS_V2(NU,OBSF,NCF,NCF_MAX,
     1              T_IN_K,V_TURB,LOG_NTOT,
     1	            LOG_H2_NTOT,V_R,MIN_RES_KMS,WAVE_MAX,WAVE_MIN,
     1	            HI_ABS,H2_ABS,IS_FILE)
!
!
!
	ELSE IF(X(1:5) .EQ. 'CNVLV') THEN
!
	 IF(NCF .EQ. 0)THEN
	    WRITE(T_OUT,*)'Error: no data in buffer to operate on.'
	    WRITE(T_OUT,*)'Use rd_mod(OVER=T) with rd_mod option, or'
	    WRITE(T_OUT,*)'Use rd_obs(MOD=T) wit rd_obs option.'
	    GOTO 1		!Get another option
	 END IF
!
! Instrumental profiles is assumed to be Gaussian.
!
	 WRITE(T_OUT,*)' '
	 WRITE(T_OUT,*)' Two choices are possible: '
	 WRITE(T_OUT,*)'   (1) Fixed resolution (INST_RES) in angstroms'
	 WRITE(T_OUT,*)'   (2) Fixed velocity resoluton (Lam/dLam) '
	 WRITE(T_OUT,*)' The non-zero value is used'
	 WRITE(T_OUT,*)' '
!
! If the RESOLUTION is non zero, it gets used instead of INST_RES.
! Resolution allows a spectrum to be smoothed to a constant velocity
! resolution at all wavelengths.
!
100	 CALL USR_OPTION(INST_RES,'INST_RES','0.0D0','Instrumental Resolution in Angstroms [dLam - FWHM]')
	 CALL USR_OPTION(RESOLUTION,'RES','0.0D0','Resolution [Lam/dLam(FWHM)] (km/s if -ve)')
	 IF(RESOLUTION .LT. 0.0_LDP)THEN
	   RESOLUTION=1.0E-05_LDP*C_CMS/ABS(RESOLUTION)
	 ELSE IF(RESOLUTION .EQ. 0.0_LDP .AND. INST_RES .EQ. 0.0_LDP)THEN
	   WRITE(T_OUT,*)'Only one INST_RES and RES can be zero'
	   GOTO 100
	 ELSE IF(RESOLUTION .NE. 0.0_LDP .AND. INST_RES .NE. 0.0_LDP)THEN
	   WRITE(T_OUT,*)'Only one INST_RES and RES can be non-zero'
	   GOTO 100
	 END IF
!
! Defaults are those of HUT.
!
	 CALL USR_OPTION(WAVE_MIN,'WAVE_MIN','900d0','Minimum Wavelength')
	 CALL USR_OPTION(WAVE_MAX,'WAVE_MAX','10000d0','Maximum Wavelength')
	 CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','1.0d0',
	1           'Minimum Model Resolution (km/s)')
	 CALL USR_HIDDEN(NUM_RES,'NUM_RES','5.0d0',
     1              'Number of resolution elements condidered outside bandpass')
	 FFT_CONVOLVE=.FALSE.
!	 CALL USR_HIDDEN(FFT_CONVOLVE,'FFT','F',
!     1              'Use FFT methods for convolution')
!
	 VSINI=0.0_LDP		!For rotational broadening, so set to zero
	 EPSILON=0.0_LDP
	 CALL SMEAR_V2(NU,OBSF,NCF,
	1	      WAVE_MAX,WAVE_MIN,
	1             INST_RES,RESOLUTION,VSINI,EPSILON,
	1             MIN_RES_KMS,NUM_RES,FFT_CONVOLVE)
!
	ELSE IF(X(1:3) .EQ. 'ROT') THEN
!
	  IF(NCF .EQ. 0)THEN
	    WRITE(T_OUT,*)'Error: no data in buffer to operate on.'
	    WRITE(T_OUT,*)'Use rd_mod(OVER=T) with rd_mod option, or'
	    WRITE(T_OUT,*)'Use rd_obs(MOD=T) wit rd_obs option.'
	    GOTO 1		!Get another option
	  END IF
!
! Perform a crude rotational broadening. This is not meant to be rigorous.
!
	 CALL USR_OPTION(VSINI,'VSINI','100.0D0','Vsini')
	 CALL USR_OPTION(EPSILON,'EPS','0.5D0',
	1            'I(mu)/I(mu=1) = 1-eps + eps*mu')
!
! Defaults are observable spectral region.
!
	 CALL USR_OPTION(WAVE_MIN,'WAVE_MIN','900.0d0','Minimum Wavelength')
	 CALL USR_OPTION(WAVE_MAX,'WAVE_MAX','7000.0d0','Maximum Wavelength')
!
! If the RESOLUTION is non zero, it gets used instead of INST_RES.
! Resolution allows a spectrum to be smoothed to a constant velocity
! resolution at all wavelengths.
!
	 CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','1.0d0',
	1                      'Minimum Model Resolution (km/s)')
	 CALL USR_HIDDEN(NUM_RES,'NUM_RES','5.0d0',
     1                   'Number of resolution elements condidered')
!
! Rotational broadening cannot use FFT option at present.
!
	 FFT_CONVOLVE=.FALSE.
	 CALL SMEAR_V2(NU,OBSF,NCF,
	1	         WAVE_MAX,WAVE_MIN,
	1                INST_RES,RESOLUTION,VSINI,EPSILON,
	1                MIN_RES_KMS,NUM_RES,FFT_CONVOLVE)
!
	ELSE IF(X(1:3) .EQ. 'EXT') THEN
!
! Extract spectrum at a given pixel resolution.
!
	  CALL USR_OPTION(RESOLUTION,'RES','1000.0D0','R / /\R ')
	  I=1
	  T1=NU(1)
	  ML=1
	  DO WHILE(T1 .GT. NU(NCF))
	    DO WHILE( T1 .LT. NU(ML+1) )
	      ML=ML+1
	    END DO
	    T2=(T1-NU(ML))/(NU(ML+1)-NU(ML))
	    YV(I)=(1.0_LDP-T2)*OBSF(ML)+T2*OBSF(ML+1)
	    XV(I)=T1
	    T1=T1*RESOLUTION/(1.0_LDP+RESOLUTION)
	    I=I+1
	  END DO
	  NCF=I-1
	  NU(1:NCF)=XV(1:NCF)
	  OBSF(1:NCF)=YV(1:NCF)
C 
C
C This option s a simpliied combination of RD_MOD, ROT, and NORM. It is specifically
C design for examining normalized model spectra. The separate options are more general.
C
	ELSE IF(X(1:3) .EQ. 'GEN')THEN
	  DIRECTORY=' '
	  CALL USR_OPTION(DIRECTORY,'Dir',' ','Model directory')
	  CALL USR_HIDDEN(FILENAME,'File','obs_fin','File (def=obs_fin)')
C
	  IF(INDEX(DIRECTORY,']') .EQ. 0)THEN
	    FILENAME=TRIM(DIRECTORY)//'/obs/'//TRIM(FILENAME)
	  ELSE
	    I=LEN_TRIM(DIRECTORY)
	    CALL SET_CASE_UP(FILENAME,IZERO,IZERO)
	    FILENAME=DIRECTORY(1:I-1)//'.OBS]'//TRIM(FILENAME)
	  END IF
	  CALL RD_MOD(NU,OBSF,NCF_MAX,NCF,FILENAME,IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error reading model data -- no data read'
	     GOTO 1		!Get another option
	  END IF
	  WRITE(T_OUT,*)'New model data replaces old data'
!
! Perform a crude rotational broadening. This is not meant to be rigorous.
!
	  CALL USR_OPTION(VSINI,'VSINI','100.0D0','Vsini')
	  CALL USR_OPTION(EPSILON,'EPS','0.5D0',
	1            'I(mu)/I(mu=1) = 1-eps + eps*mu')
!
! Defaults are observable spectral region.
!
	  CALL USR_OPTION(WAVE_MIN,'WAVE_MIN','900.0d0','Minimum Wavelength')
	  CALL USR_OPTION(WAVE_MAX,'WAVE_MAX','7000.0d0','Maximum Wavelength')
!
! If the RESOLUTION is non zero, it gets used instead of INST_RES.
! Resolution allows a spectrum to be smoothed to a constant velocity
! resolution at all wavelengths.
!
	  CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','1.0d0',
	1                     'Minimum Model Resolution (km/s)')
	  CALL USR_HIDDEN(NUM_RES,'NUM_RES','5.0d0',
     1                   'Number of resolution elements condidered')
!
! Rotational broadening cannot use FFT option at present.
!
	  FFT_CONVOLVE=.FALSE.
	  CALL SMEAR_V2(NU,OBSF,NCF,
	1	         WAVE_MAX,WAVE_MIN,
	1                INST_RES,RESOLUTION,VSINI,EPSILON,
	1                MIN_RES_KMS,NUM_RES,FFT_CONVOLVE)
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Rotational smoothing has been performed'
!
	  CALL USR_OPTION(T_IN_K,'T_IN_K','100d0','Temp. in Kelvin')
	  CALL USR_OPTION(V_TURB,'V_TURB','10d0','Turbulent Velocity (km/s)')
	  CALL USR_OPTION(LOG_NTOT,'LOG_NTOT','20d0','Log of H column density')
	  CALL USR_OPTION(LOG_H2_NTOT,'LOG_H2_NTOT','20d0','Log H2 column dens')
!
	  CALL USR_HIDDEN(HI_ABS,'HI_ABS','T','Correct for HI absorption')
	  CALL USR_HIDDEN(H2_ABS,'H2_ABS','T','Correct for HII absorption')
	  CALL USR_HIDDEN(V_R,'V_R','0.0d0','Radial Velocity (km/s)')
	  CALL USR_HIDDEN(WAVE_MIN,'WAVE_MINI','900d0','Minimum Wavelength')
	  CALL USR_HIDDEN(WAVE_MAX,'WAVE_MAXI','1300d0','Maximum Wavelength')
	  CALL USR_HIDDEN(IS_FILE,'IS_FILE','IS_LINE_LIST','File wth list and strengths of IS lines')
!
	  CALL UVABS_V2(NU,OBSF,NCF,NCF_MAX,
     1              T_IN_K,V_TURB,LOG_NTOT,
     1	            LOG_H2_NTOT,V_R,MIN_RES_KMS,WAVE_MAX,WAVE_MIN,
     1	            HI_ABS,H2_ABS,IS_FILE)
!
	  WRITE(T_OUT,*)'ISABS correction has been performed'
!

	  IF(INDEX(DIRECTORY,']') .EQ. 0)THEN
	    FILENAME=TRIM(DIRECTORY)//'/obs/obs_cont'
	  ELSE
	    I=LEN_TRIM(DIRECTORY)
	    FILENAME=DIRECTORY(1:I-1)//'.OBS]OBS_cont'
	  END IF
	  CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_CONT,FILENAME,IOS)
	  IF(IOS .NE. 0)GOTO 1		!Get another option
!
	  T1=1.0E-08_LDP
	  I=1
	  UNEQUAL=.FALSE.
	  IF(NCF_CONT .NE. NCF)UNEQUAL=.TRUE.
	  DO WHILE(.NOT. UNEQUAL .AND. I .LE. NCF_CONT)
	    IF( EQUAL(NU_CONT(I),NU(I),T1) )THEN
	      XV(I)=NU_CONT(I)
	      I=I+1
	    ELSE
	      UNEQUAL=.TRUE.
	    END IF
	  END DO
	  IF(UNEQUAL)CALL USR_HIDDEN(LIN_INT,'LIN','F',
	1               'Use linear interpolation')
	  IF(UNEQUAL .AND. LIN_INT)THEN
	    L=1
	    DO I=1,NCF
  	      XV(I)=NU(I)
	      IF(NU(I) .GT. NU_CONT(1))THEN
	        YV(I)=0.0
	      ELSE IF(NU(I) .LT. NU_CONT(NCF_CONT))THEN
	        YV(I)=0.0
	      ELSE
	        DO WHILE (NU(I) .LT. NU_CONT(L+1))
	          L=L+1
	        END DO
	        T1=(NU(I)-NU_CONT(L+1))/(NU_CONT(L)-NU_CONT(L+1))
	        T2=(1.0_LDP-T1)*OBSF_CONT(L+1)+T1*OBSF_CONT(L)
	        YV(I)=0.0
	        IF(T2 .NE. 0)THEN
	          T2=OBSF(I)/T2
	          IF(T2 .LT. 1.0E+020_LDP)YV(I)=T2
	        END IF
	      END IF
	    END DO
	  ELSE IF(UNEQUAL)THEN
C
C We will use monotonic cubic interpolation. We first verify the range.
C I & J are temporary variables for the callt o MON_INTERP. I denotes the
C first element. Initially J denotes the last element, then the numer of
C elements that can be interpolated.
C
	    I=1
	    DO WHILE(NU(I) .GT. NU_CONT(1))
	      I=I+1
	    END DO
	    J=NCF
	    DO WHILE(NU(J) .LE. NU_CONT(NCF_CONT))
	      J=J-1
	    END DO
	    J=J-I+1
C
	    FQW(1:NCF)=0.0_LDP				!Temporary usage
	    CALL MON_INTERP(FQW(I),J,IONE,NU(I),J,
	1            OBSF_CONT,NCF_CONT,NU_CONT,NCF_CONT)
  	    XV(1:NCF)=NU(1:NCF)
	    DO I=1,NCF
	      IF(FQW(I) .GT. 0)THEN
	         T2=OBSF(I)/FQW(I)
	         IF(T2 .LT. 1.0E+20_LDP)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  ELSE
	    DO I=1,NCF
	      YV(I)=0
	      IF(OBSF_CONT(I) .GT. 0)THEN
	         T2=OBSF(I)/OBSF_CONT(I)
	         IF(T2 .LT. 1.0E+20_LDP)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  END IF
!
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL CURVE_LAB(NCF,XV,YV,TITLE)
	  YAXIS='F\d\gn\u/F\dc\u'
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Data has been normalized and plotted'
!
!
!
! Option to do a least fit CCM redenning Law. The observational data
! should have been read in using RD_OBS, and must be contained in one file.
!
	ELSE IF(X(1:3) .EQ. 'RED')THEN
	  IF(NOBS .EQ. 0)THEN
	    WRITE(6,*)'Error -- observational data has not been read in'
	  ELSE IF(NCF .EQ. 0)THEN
	    WRITE(6,*)'Error -- model data has not been strored in the buffer'
	  ELSE
	    T1=2.5_LDP; T2=5.5
	    CALL DETERM_REDDENING(OBSF_OBS,NU_OBS,NOBS,OBSF,NU,NCF,T1,T2)
	  END IF
!
	ELSE IF(X(1:4) .EQ. 'FLAM' .OR.
	1     X(1:4) .EQ. 'WRFL' .OR.  X(1:3) .EQ. 'FNU' .OR.
	1                X(1:4) .EQ. 'EBMV') THEN
!
! If NCF is defined, we use that freuqency grid when computing the extinction curve.
! Otherwise, we define the grid.
!
	  IF(X(1:4) .EQ. 'EBMV' .AND. NCF .EQ. 0)THEN
	    T1=ANG_TO_HZ/900.0_LDP; T2=ANG_TO_HZ/5.0E+04_LDP
	    NCF=1000
	    T2=EXP(LOG(T1/T2)/(NCF-1))
	    NU(1)=T1
	    DO I=2,NCF
	      NU(I)=NU(I-1)/T2
	    END DO
	  END IF
!
	  CALL USR_OPTION(EBMV_CCM,'EBMV_CCM','0.0D0',
	1             'CCM E(B-V) to correct for I.S. extinction')
	  IF(EBMV_CCM .NE. 0)THEN
	    CALL USR_OPTION(r_ext,'R_EXT','3.1D0',
	1         'R_EXT for Cardelli, Clayton, Mathis extinction law')
	  END IF
!
	  CALL USR_OPTION(EBMV_GAL,'EBMV_GAL','0.0',
	1             'Galactic E(B-V) to correct for I.S. extinction')
	  CALL USR_OPTION(EBMV_LMC,'EBMV_LMC','0.0',
	1             'LMC E(B-V) to correct for I.S. extinction')
	  CALL USR_OPTION(EBMV_SMC,'EBMV_SMC','0.0',
	1      'SMC E(B-V) to correct for I.S. extinction') !KN SMC
!
	  DIST=1.0_LDP
	  CALL USR_OPTION(DIST,'DIST','1.0D0',' (in kpc) ')
!
	  OVER=.FALSE.
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite store')
!
	  DO I=1,NCF
	    XV(I)=NU(I)
	    YV(I)=OBSF(I)/DIST/DIST
	  END DO
C
	  IF(EBMV_CCM .NE. 0.0_LDP)THEN
	    DO I=1,NCF
	      T1=ANG_TO_HZ/NU(I)
	      T1=(10000.0_LDP/T1)				!1/Lambda(um)
	      IF(T1 .LT. 1.1_LDP)THEN
	        RAX=0.574_LDP*(T1**1.61_LDP)
	        RBX=-0.527_LDP*(T1**1.61_LDP)
	      ELSE IF(T1. LT. 3.3_LDP)THEN
	        T2=T1-1.82_LDP
	        RAX=1+T2*(0.17699_LDP-T2*(0.50447_LDP+T2*(0.02427_LDP-T2*(0.72085_LDP
	1                +T2*(0.01979_LDP-T2*(0.77530_LDP-0.32999_LDP*T2))))))
	        RBX=T2*(1.41338_LDP+T2*(2.28305_LDP+T2*(1.07233_LDP-T2*(5.38434_LDP
	1                +T2*(0.62251_LDP-T2*(5.30260_LDP-2.09002_LDP*T2))))))
	      ELSE IF(T1 .lT. 5.9_LDP)THEN
	        RAX=1.752_LDP-0.316_LDP*T1-0.104_LDP/((T1-4.67_LDP)**2+0.341_LDP)
	        RBX=-3.090_LDP+1.825_LDP*T1+1.206_LDP/((T1-4.62_LDP)**2+0.263_LDP)
	      ELSE IF(T1 .LT. 8.0_LDP)THEN
  	        T2=T1-5.9
	        RAX=1.752_LDP-0.316_LDP*T1-0.104_LDP/((T1-4.67_LDP)**2+0.341_LDP) -
	1                       T2*T2*(0.04773_LDP+0.009779_LDP*T2)
	        RBX=-3.090_LDP+1.825_LDP*T1+1.206_LDP/((T1-4.62_LDP)**2+0.263_LDP)+
	1                       T2*T2*(0.2130_LDP+0.1207_LDP*T2)
	      ELSE IF(T1 .LT. 10)THEN
	        T2=T1-8
	        RAX=-1.073_LDP-T2*(0.628_LDP-T2*(0.137_LDP-0.070_LDP*T2))
	        RBX=13.670_LDP+T2*(4.257_LDP-T2*(0.420_LDP-0.374_LDP*T2))
	      ELSE
	        T1=10
	        T2=T1-8
	        RAX=-1.073_LDP-T2*(0.628_LDP-T2*(0.137_LDP-0.070_LDP*T2))
	        RBX=13.670_LDP+T2*(4.257_LDP-T2*(0.420_LDP-0.374_LDP*T2))
	      END IF
              AL_D_EBmV(I)=R_EXT*(RAX+RBX/R_EXT)
	    END DO
	    IF(X(1:4) .EQ. 'EBMV')THEN
	      DO I=1,NCF
	        XV(I)=NU(I)
	        YV(I)=EBMV_CCM*AL_D_EBmV(I)
	      END DO
 	      XAXIS='\gl(\A)'
	      YAXIS='A\d\gl\u'
	    ELSE
	      DO I=1,NCF
	        YV(I)=YV(I)*( 10.0_LDP**(-0.4_LDP*EBMV_CCM*AL_D_EBmV(I)) )
	      END DO
	    END IF
	  END IF
C
C Set galactic interstellar extinction curve. Curve is from Howarth
C (1983, MNRAS, 203, 301) and Seaton (1979, MNRAS, 187, 73P).
C
	  IF(EBMV_GAL .NE. 0)THEN
	    R_EXT=3.1
	    DO I=1,NCF
	      T1=ANG_TO_HZ/NU(I)
	      T1=(10000.0_LDP/T1)				!1/Lambda(um)
	      IF( T1 .LT. 1.83_LDP)THEN
	        AL_D_EBmV(I)=((1.86_LDP-0.48_LDP*T1)*T1-0.1_LDP)*T1
	      ELSE IF( T1 .LT. 2.75_LDP)THEN
	        AL_D_EBmV(I)=R_EXT+2.56_LDP*(T1-1.83_LDP)-0.993_LDP*(T1-1.83_LDP)**2
	      ELSE IF( T1 .LT. 3.65_LDP)THEN
	        AL_D_EBmV(I)=(R_EXT-1.64_LDP)+1.048_LDP*T1+1.01_LDP/( (T1-4.60_LDP)**2+0.28_LDP )
	      ELSE IF( T1 .LT. 7.14_LDP)THEN
	        AL_D_EBmV(I)=(R_EXT-0.91_LDP)+0.848_LDP*T1+1.01_LDP/( (T1-4.60_LDP)**2+0.28_LDP )
	      ELSE IF( T1 .LT. 11)THEN
                AL_D_EBmV(I)=(R_EXT+12.97_LDP)-3.20_LDP*T1+0.2975_LDP*T1*T1
	      ELSE
	        T1=11.0
                AL_D_EBmV(I)=(R_EXT+12.97_LDP)-3.20_LDP*T1+0.2975_LDP*T1*T1
	      END IF
	    END DO
	  END IF
	  IF(X(1:4) .EQ. 'EBMV' .AND. EBMV_GAL .NE. 0)THEN
	    DO I=1,NCF
	      XV(I)=NU(I)
	      YV(I)=EBMV_GAL*AL_D_EBmV(I)
	    END DO
	    XAXIS='\gl(\A)'
	    YAXIS='A\d\gl\u'
	  ELSE IF(X(1:4) .NE. 'EBMV' .AND. EBMV_GAL .NE. 0)THEN
	    DO I=1,NCF
	      YV(I)=YV(I)*( 10.0_LDP**(-0.4_LDP*EBMV_GAL*AL_D_EBmV(I)) )
	    END DO
	  END IF
C
C Set LMC interstellar extinction curve. Curve is from Howarth
C (1983, MNRAS, 203, 301).
C
	  IF(EBMV_LMC .NE. 0)THEN
	    R_EXT=3.1
	    DO I=1,NCF
	      T1=ANG_TO_HZ/NU(I)
	      T1=(10000.0_LDP/T1)				!1/Lambda(um)
	      IF( T1 .LT. 1.83_LDP)THEN
    	        AL_D_EBmV(I)=((1.86_LDP-0.48_LDP*T1)*T1-0.1_LDP)*T1
	      ELSE IF( T1 .LT. 2.75_LDP)THEN
	        AL_D_EBmV(I)=R_EXT+2.04_LDP*(T1-1.83_LDP)+0.094_LDP*(T1-1.83_LDP)**2
	      ELSE IF( T1 .LT. 11.0_LDP)THEN	!Strictly only below 9
	        AL_D_EBmV(I)=(R_EXT-0.236_LDP)+0.462_LDP*T1+0.105_LDP*T1*T1+
	1                       0.454_LDP/( (T1-4.557_LDP)**2+0.293_LDP )
	      ELSE
	        T1=11.0
	        AL_D_EBmV(I)=(R_EXT-0.236_LDP)+0.462_LDP*T1+0.105_LDP*T1*T1+
	1                       0.454_LDP/( (T1-4.557_LDP)**2+0.293_LDP )
	      END IF
	    END DO
	  END IF
C
	  IF(X(1:4) .EQ. 'EBMV' .AND. EBMV_LMC .NE. 0)THEN
	    DO I=1,NCF
	      XV(I)=NU(I)
!	      XV(I)=ANG_TO_HZ/NU(I)
	      YV(I)=EBMV_LMC*AL_D_EBmV(I)
	    END DO
	    XAXIS='\gl(\A)'
	    YAXIS='A\d\gl\u'
	  ELSE IF(X(1:4) .NE. 'EBMV' .AND. EBMV_LMC .NE. 0)THEN
	    DO I=1,NCF
	      YV(I)=YV(I)*( 10.0_LDP**(-0.4_LDP*EBMV_LMC*AL_D_EBmV(I)) )
	    END DO
	  END IF
C
C       KN: Set SMC interstellar extinction curve.
C
	  IF(EBMV_SMC .NE. 0)THEN
	     R_EXT=2.74_LDP
	     DO I=1,NCF
	        T1=ANG_TO_HZ/NU(I)
	        T1=(10000.0_LDP/T1) !1/Lambda(um)
	        IF(T1 .LT. 1.1_LDP)THEN
	          RAX=0.574_LDP*(T1**1.61_LDP)
	          RBX=-0.527_LDP*(T1**1.61_LDP)
	          AL_D_EBmV(I)=R_EXT*(RAX+RBX/R_EXT)
	        ELSE IF(T1. LT. 3.39_LDP)THEN
	          T2=T1-1.82_LDP
	          RAX=1+T2*(0.17699_LDP-T2*(0.50447_LDP+T2*(0.02427_LDP-T2*(0.72085_LDP
	1                +T2*(0.01979_LDP-T2*(0.77530_LDP-0.32999_LDP*T2))))))
	          RBX=T2*(1.41338_LDP+T2*(2.28305_LDP+T2*(1.07233_LDP-T2*(5.38434_LDP
	1                 +T2*(0.62251_LDP-T2*(5.30260_LDP-2.09002_LDP*T2))))))
	          AL_D_EBmV(I)=R_EXT*(RAX+RBX/R_EXT)
	        ELSE
	          C1=-4.959_LDP
	          C2=2.264_LDP*T1
	          D=(T1**2)/(((T1**2-4.6_LDP**2)**2)+(T1**2))
	          C3=0.389_LDP*D
	          IF(T1 .LT. 5.9_LDP)THEN
	             F=0
	          ELSE
	             F=0.5392_LDP*((T1-5.9_LDP)**2)+0.05644_LDP*((T1-5.9_LDP)**3)
	          END IF
	          C4=0.461_LDP*F
                  AL_D_EBmV(I)=C1+C2+C3+C4+R_EXT
               END IF
	     END DO
          END IF
	  IF(X(1:4) .EQ. 'EBMV' .AND. EBMV_SMC .NE. 0)THEN
	     DO I=1,NCF
	        XV(I)=NU(I)
	        YV(I)=EBMV_SMC*AL_D_EBmV(I)
	     END DO
	     XAXIS='\gl(\A)'
	     YAXIS='A\d\gl\u'
	  ELSE IF(X(1:4) .NE. 'EBMV' .AND. EBMV_SMC .NE. 0)THEN
	     DO I=1,NCF
	        YV(I)=YV(I)*( 10.0_LDP**(-0.4_LDP*EBMV_SMC*AL_D_EBmV(I)) )
	     END DO
	  END IF
C
	  IF(X(1:4) .NE. 'EBMV')THEN
	    CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	  ELSE
	    CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  END IF
 	  IF(X(1:4) .Eq. 'WRFL')THEN
	    DO I=1,NCF
	      WRITE(50,*)XV(I),YV(I)
	    END DO
	  ELSE IF(OVER)THEN
	    OBSF(1:NCF)=YV(1:NCF)
	  ELSE
	    CALL CURVE_LAB(NCF,XV,YV,TITLE)
	  END IF
	ELSE IF(X(1:2) .EQ.'RV')THEN
	  RAD_VEL=0.0_LDP
	  CALL USR_OPTION(RAD_VEL,'RAD_VEL','0.0D0','Radial velcoity (+ve if away)(km/s)')
	  RAD_VEL=1.0E+05_LDP*RAD_VEL
	  DO I=1,NCF
	    NU(I)=NU(I)*(1.0_LDP-RAD_VEL/C_CMS)
	  END DO
C
C The follwing option cumputes the luminosity below the 3 main H/He edges,
C and the number of photons emitted.
C
	ELSE IF(X(1:3) .EQ.'ZAN')THEN
	  CALL TRAPUNEQ(NU,FQW,NCF)
	  CALL USR_OPTION(ZAN_FREQ,5,1,'LEVS','0,0,0,0,0,0',
	1      'Edge frequencies to compute photon flux (not H/He)')
C
	  ZAN_FREQ(4:8)=ZAN_FREQ(1:5)
	  ZAN_FREQ(1)=EDGE_HYD
	  ZAN_FREQ(2)=EDGE_HEI
	  ZAN_FREQ(3)=EDGE_HE2
C
	  TOT_LUM=0.0_LDP
	  DO I=1,NCF
	    TOT_LUM=TOT_LUM+OBSF(I)*FQW(I)
	  END DO
	  TOT_LUM=TOT_LUM*312.7_LDP			!4pi*(1kpc)**2*(1E+15)*(1D-23)/Lsun
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Total luminosity is:               ',TOT_LUM
C
	  ML=1
	  DO WHILE(ZAN_FREQ(ML) .GT. 0.0_LDP)
	    TMP_FREQ=ZAN_FREQ(ML)
	    N_PHOT=0.0_LDP
	    LUM=0.0_LDP
	    J=1
	    DO WHILE(NU(J+1) .GT. TMP_FREQ)
	      J=J+1
	    END DO
	    CALL TRAPUNEQ(NU,FQW,J)
	    DO I=1,J
              N_PHOT=N_PHOT+FQW(I)*OBSF(I)/NU(I)
	      LUM=LUM+FQW(I)*OBSF(I)
	    END DO
	    LUM=LUM*312.7_LDP			!4pi*(1kpc)**2*(1E+15)*(1D-23)/Lsun
	    N_PHOT=47.2566_LDP+LOG10(N_PHOT) !LOG10(4pi*(1kpc)**2*(1D-23)/h)
	    T1=ANG_TO_HZ/TMP_FREQ
	    WRITE(T_OUT,*)' '
	    WRITE(T_OUT,'(1X,A,F4.0,A,1PE11.4)')
	1        'Luminosity shortward of        ',T1,'A is:   ',LUM
	    WRITE(T_OUT,'(1X,A,F4.0,A,1PE11.4)')
	1        'Log(#) of photons shortward of ',T1,'A is:   ',N_PHOT
	    ML=ML+1
	  END DO
C
C
C The follwing option plot the cummulative phopton flux as a function of
C frequency.
C
	ELSE IF(X(1:4) .EQ. 'PZAN')THEN
C
	  XV(1:NCF)=NU(1:NCF)
	  YV(1)=0.0_LDP
	  DO ML=2,NCF
	    YV(ML)=YV(ML-1)+ 0.5_LDP*(NU(ML-1)-NU(ML))*
	1              (OBSF(ML-1)/NU(ML-1)+OBSF(ML)/NU(ML))
	  END DO
	  DO ML=2,NCF
	    IF(YV(ML) .LT. 1.0E-30_LDP)YV(ML)=1.0E-30_LDP
	    YV(ML)=47.2566_LDP+LOG10(YV(ML)) 		!LOG10(4pi*(1kpc)**2*(1E-23)/h)
	  END DO
C
	  I=NCF-1
	  CALL CNVRT(XV(2),YV(2),I,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL CURVE_LAB(I,XV(2),YV(2),TITLE)
	  YAXIS='Log(N)'
C
	ELSE IF(X(1:3) .EQ.'CUM')THEN
	  NORM_LUM=0.0_LDP
	  CALL USR_OPTION(NORM_LUM,'NORM_LUM',' ',' (in Lsun) ')
	  YV(1)=0.0_LDP
	  DO I=2,NCF
	    YV(I)=YV(I-1)+(OBSF(I)+OBSF(I-1))*(NU(I-1)-NU(I))*0.5_LDP
	  END DO
	  T2=YV(NCF)*312.7_LDP		!4pi*(1kpc)**2*(1D+15)*(1D-23)/Lsun
	  IF(NORM_LUM .EQ. 0.0_LDP)THEN
	    NORM_LUM=YV(NCF)
	  ELSE
	    NORM_LUM=NORM_LUM/312.7_LDP
	  END IF
	  DO I=1,NCF
	    XV(I)=NU(I)
	    YV(I)=YV(I)/NORM_LUM
	  END DO
	  WRITE(T_OUT,*)'Total Luminosity is',T2
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL CURVE_LAB(NCF,XV,YV,TITLE)
C
	ELSE IF(X(1:4) .EQ. 'FILT')THEN
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=OBSF(1:NCF)
	  CALL CNVRT(XV,YV,NCF,L_FALSE,L_FALSE,'ANG','FLAM',
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	  DO IF=1,NF
!
	    SP_VAL=FILT_ST_LAM(IF); ML_ST=GET_INDX_SP(SP_VAL,XV,NCF)-1
	    SP_VAL=FILT_ST_LAM(IF)+24.0_LDP*FILT_DEL_LAM(IF)
	    ML_END=GET_INDX_SP(SP_VAL,XV,NCF)+1
!
	    WRITE(T_OUT,*)'Filter ',IF
	    WRITE(T_OUT,*)FILT_ST_LAM(IF),XV(ML_ST),XV(ML_END),YV(ML_ST),YV(ML_END)
	    WRITE(T_OUT,*)ML_ST,ML_END
!
	    ZV(IF)=0.0_LDP
	    SUM=0.0_LDP
	    RESPONSE2=0.0_LDP
	    I=1
	    DO ML=ML_ST,ML_END
	      IF(XV(ML) .GE. FILT_ST_LAM(IF)+24.0_LDP*FILT_DEL_LAM(IF))EXIT
	      DO WHILE(XV(ML) .GT. FILT_ST_LAM(IF)+I*FILT_DEL_LAM(IF))
	        I=I+1
	      END DO
	      FILT_INT_BEG =FILT_ST_LAM(IF)+(I-1)*FILT_DEL_LAM(IF)
	      RESPONSE1=RESPONSE2
	      T1=(XV(ML+1)-FILT_INT_BEG)/FILT_DEL_LAM(IF)
	      RESPONSE2=(1.0_LDP-T1)*NORM_PASS(I,IF)+T1*NORM_PASS(I+1,IF)
	      IF(ABS(T1) .GT. 1)RESPONSE2=0.0_LDP
	      ZV(IF)=ZV(IF)+0.5_LDP*(RESPONSE1*YV(ML)+RESPONSE2*YV(ML+1))*(XV(ML+1)-XV(ML))
	      SUM=SUM+0.5_LDP*(RESPONSE1+RESPONSE2)*(XV(ML+1)-XV(ML))
	    END DO
	    WRITE(T_OUT,*)ZV(IF),SUM,ZV(IF)/SUM
	    ZV(IF)=-FILT_ZP(IF)-2.5_LDP*LOG10(ZV(IF)/SUM)
	    WRITE(T_OUT,'(A,1PE10.3)')TRIM(FILT_NAME(IF)),ZV(IF)
	  END DO
	  DO IF=1,NF
            WRITE(T_OUT,'(A,F12.3)')TRIM(FILT_NAME(IF)),ZV(IF)
	  END DO

	ELSE IF(X(1:3) .EQ. 'MAG')THEN
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'The magnitudes computed here are NOT accurate.'
	  WRITE(T_OUT,*)'Filter convolutions need to be implemented.'
	  WRITE(T_OUT,*)'Originally designed for crude estimates using continuum fluxe.s'
	  WRITE(T_OUT,*)' '
	  DIST=1.0
	  CALL USR_OPTION(DIST,'DIST','1.0D0',' (in kpc) ')
	  CALL USR_OPTION(FILTER_SET,'FSET','uvby','Filter set (case sensitive)')
	  CALL GEN_ASCI_OPEN(LU_OUT,'MAG','UNKNOWN',' ',' ',IZERO,IOS)
!
	  IF( UC(FILTER_SET) .EQ. 'ALL')THEN
	    FILTER_SET='ubvy'
	    CALL GET_MAG(NU,OBSF,NCF,DIST,FILTER_SET,LU_OUT)
	    FILTER_SET='UBV'
	    CALL GET_MAG(NU,OBSF,NCF,DIST,FILTER_SET,LU_OUT)
	  ELSE
	    CALL GET_MAG(NU,OBSF,NCF,DIST,FILTER_SET,LU_OUT)
	  END IF
!
!	  WRITE(LU_OUT,104)DIST	
!	  DO L=1,8
!	    MLST=1
!	    DO ML=MLST,NCF
!	      IF(NU(ML) .LT. 0.2998/FILTLAM(L))THEN
!	        T1=LOG10(OBSF(ML-1))-LOG10(NU(ML-1)*FILTLAM(L)/0.2998)
!	1       *(LOG10(OBSF(ML)/OBSF(ML-1))/(LOG10(NU(ML)/NU(ML-1))))
!	        T1=5.0*LOG10(DIST)-2.5*T1+2.5*LOG10(FILTZP(L))
!	        MLST=ML
!	        GOTO 102
!	       END IF
!	    END DO
!102	  CONTINUE
!	    WRITE(LU_OUT,103)FILTLAM(L),FILTZP(L),FILT(L),T1
!103	    FORMAT(2X,F7.4,5X,F7.1,8X,A1,5X,F6.2)
!104	    FORMAT(2X,' Assumed Distance is',F5.1,' kpc',/
!	1  ,/,2x,'    Lam          ZP     Filt        Mag'//)
!	  END DO
!	  WRITE(LU_OUT,107)
!	  DO L=1,24
!	    MLST=1
!	    DO ML=MLST,NCF
!	      IF(NU(ML) .LT. 0.2998/FLAM(L))THEN
!	        T1=LOG10(OBSF(ML-1))-LOG10(NU(ML-1)*FLAM(L)/0.2998)
!	1       *(LOG10(OBSF(ML)/OBSF(ML-1))/(LOG10(NU(ML)/NU(ML-1))))
!	        T1=5.0*LOG10(DIST)-2.5*T1+2.5*LOG10(ZERO_POINT)
!	        MLST=ML
!	        GOTO 105
!	       END IF
!	    END DO
!105	  CONTINUE
!	    WRITE(LU_OUT,106)FLAM(L),ZERO_POINT,T1
!106	    FORMAT(2X,F7.4,5X,F7.1,5X,F6.2)
!107	    FORMAT(//,'     Lam          ZP     Mag//')
!	  END DO
!	  CLOSE(UNIT=LU_OUT)
!
C
	ELSE IF(X .EQ. 'BB')THEN
	
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model data')
	  CALL USR_OPTION(TEMP,'TEMP','3.0',' ')
!
	  IF(NCF .EQ. 0)THEN
	     T2=EXP(LOG(1.0E+4_LDP)/(NBB-1))
	     T1=3.289_LDP*912.0_LDP/10.0_LDP
	     CALL USR_HIDDEN(LIN_INT,'NORM','T','Normalize peak to unity?')
	     DO I=1,NBB
	       NU(I)=T1/(T2**I)
	       T3=HDKT*NU(I)/TEMP
	       IF(T3 .GT. 1.0_LDP)THEN
	         YV(I)=TWOHCSQ*(NU(I)**3)*EXP(-T3)/(1.0_LDP-EXP(-T3))
	       ELSE
	         YV(I)=TWOHCSQ*(NU(I)**3)/(EXP(T3)-1.0_LDP)
	       END IF
	       XV(I)=NU(I)
	     END DO
	     IF(LIN_INT)THEN
	       T1=MAXVAL(YV(1:NBB))
	       YV(1:NBB)=YV(1:NBB)/T1
	     END IF
	  ELSE
	    CALL USR_OPTION(NORM_WAVE,'NW',' ',
	1       'Norm Wave (Angstroms) or Radius (<0 [in Rsun])')
C
C Rather than choose the full NCF points we adopt a set that extends from
C NU_MAX to NU_MIN but with a larger spacing.
C
	    DNU=LOG10(NU(NCF)/NU(1))/(NBB-1)
	    DO I=1,NBB
	      T1=NU(1)*10.0_LDP**(DNU*(I-1))
	       T3=HDKT*T1/TEMP
	       IF(T3 .GT. 1.0_LDP)THEN
	         YV(I)=TWOHCSQ*(T1**3)*EXP(-T3)/(1.0_LDP-EXP(-T3))
	       ELSE
	         YV(I)=TWOHCSQ*(T1**3)/(EXP(T3)-1.0_LDP)
	       END IF
	      XV(I)=T1
	    END DO
	    IF(NORM_WAVE .EQ. 0.0_LDP)THEN
	      T1=4.0_LDP*ATAN(1.0_LDP)    !PI
	      DO I=1,NBB
	        YV(I)=YV(I)*T1
	      END DO
	    ELSE IF(NORM_WAVE .GT. 0.0_LDP)THEN
	      NORM_FREQ=2.998E+03_LDP/NORM_WAVE
	      DO I=2,NCF
	        IF(NORM_FREQ .LE. NU(I-1) .AND. NORM_FREQ .GT. NU(I))K=I
	      END DO
	      T1=LOG(NU(K-1)/NORM_FREQ)/LOG(NU(K-1)/NU(K))
	      NORM_FLUX=EXP( (1.0_LDP-T1)*LOG(OBSF(K-1))+T1*LOG(OBSF(K)) )
	      BB_FLUX=TWOHCSQ*(NORM_FREQ**3)/(EXP(HDKT*NORM_FREQ/TEMP)-1.0_LDP)
	      SCALE_FAC=NORM_FLUX/BB_FLUX
	      DO I=1,NBB
	        YV(I)=YV(I)*SCALE_FAC
	      END DO
	    ELSE
	      SCALE_FAC=NORM_WAVE*NORM_WAVE*159.8413_LDP     	 !1E+23*pi*(Rsun/1kpc)**2
	      DO I=1,NBB
	        YV(I)=YV(I)*SCALE_FAC
	      END DO
	    END IF
	  END IF
	  IF(OVER)THEN
	    OBSF(1:NBB)=YV(1:NBB)
	    NU(1:NBB)=XV(1:NBB)
	    NCF=NBB
	    WRITE(T_OUT,*)'New model data replaces old data'
	    WRITE(T_OUT,*)'No plots done with new model data'
	    WRITE(T_OUT,*)'No scaling done with new model data'
	  ELSE
	    CALL CNVRT(XV,YV,NBB,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    IF(LIN_INT .AND. NCF .EQ. 0)THEN
	       T1=MAXVAL(YV(1:NBB))
	       YV(1:NBB)=YV(1:NBB)/T1
	    END IF
	    CALL CURVE_LAB(NBB,XV,YV,TITLE)
	  END IF
!
	ELSE IF(X(1:5) .EQ. 'MKGL')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)' Note:  FWHM - 2.35482*SIGMA '
	  WRITE(6,*)' '
	  CALL USR_OPTION(LAM_CENT,'LAMC','6562','Central wavelngth')
	  CALL USR_OPTION(SIGMA,'SIGMA','1.0','Sigma (-ve for km/s)')
	  IF(SIGMA .LE. 0.0_LDP)THEN
	     SIGMA=LAM_CENT*ABS(SIGMA)/(1.0E-05_LDP*C_CMS)
	  END IF
	  CALL USR_OPTION(HEIGHT,'HEIGHT','1.0','Height (-ve for absorption line)')
	  WRITE(DEFAULT,'(ES10.4)')0.2D0*SIGMA
	  CALL USR_OPTION(T1,'dLAM',DEFAULT,'Spacing (def=SIGMA/5.0)')
	  K=20.0_LDP*SIGMA/T1
	  IF(ALLOCATED(XVEC))DEALLOCATE(XVEC,YVEC)
	  ALLOCATE (XVEC(K),YVEC(K))
	  CALL CREATE_GAUSS(XVEC,YVEC,LAM_CENT,SIGMA,HEIGHT,T1,K)
	  CALL DP_CURVE(K,XVEC,YVEC)
!
	ELSE IF(X(1:5) .EQ. 'AV_EN')THEN
	  T1=0.0_LDP
	  T2=0.0_LDP
	  DO I=1,NCF-1
	    T1=T1+0.5_LDP*(NU(I)-NU(I+1))*(OBSF(I)+OBSF(I+1))
	    T2=T2+0.5_LDP*(NU(I)-NU(I+1))*(OBSF(I)/NU(I)+OBSF(I+1)/NU(I+1))
	  END DO
	  WRITE(6,*)'Average photon frequency is (in 10^15 Hz)',T1/T2
	  WRITE(6,*)'Effective photon energy is',6.626075D-12*T1/T2/2.7D0/1.380658D-16,'K'
C 
C
C Plot section:
C
	ELSE IF(X(1:2) .EQ. 'GR') THEN
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,' ')
	  XAXIS=XAXSAV
C
	ELSE IF(X(1:4) .EQ.'GRNL') THEN
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  XAXIS=XAXSAV
C
C 
C
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:4) .EQ. 'LIST' .OR.
	1                             X(1:2) .EQ. 'HE'
	1          .OR. X(1:4) .EQ. 'HELP')THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Please see the HTML web pages in $CMFDIST/txt_files for help'
	  WRITE(T_OUT,*)' '
	ELSE IF(X(1:2) .EQ. 'EX') THEN
	  STOP
	ELSE IF(X(1:3) .EQ. 'BOX') THEN
	  CALL WR_BOX_FILE(MAIN_OPT_STR)
	ELSE
	  PRINT*,'OPTION REQUESTED DOES NOT EXIST'
	END IF
C
1	CONTINUE
	GO TO 3
C
	END
C 
C
	FUNCTION FAC(N)
	USE SET_KIND_MODULE
	REAL(KIND=LDP) FAC
	INTEGER N
	INTEGER, PARAMETER :: T_OUT=5
C
	IF(N .EQ. 0)THEN
	  FAC=1
	ELSE IF(N .LT. 0)THEN
	  WRITE(T_OUT,*)'Error in FAC --- invalid argument'
	  STOP
	ELSE
	  FAC=1
	  DO I=2,N
	    FAC=FAC*I
	  END DO
	END IF
C
	RETURN
   	END
