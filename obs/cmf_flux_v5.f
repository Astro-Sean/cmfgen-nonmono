!
! Program to compute the emergent flux in the observer's frame using an
! Observer's frame flux calculation. The program first computes the
! emissivities, opacities, and radiation field in the CMF frame. COHERENT
! or INCOHERENT electron scattering can be assumed. This data is then passed
! to OBS_FRAME_SUB for computation of the OBSERVER's flux.
!
! CMF_FLUX calling program is partially based on DISPGEN.
!
	PROGRAM CMF_FLUX_V5
	USE SET_KIND_MODULE
	USE CMF_FLUX_CNTRL_VAR_MOD
	USE MOD_CMF_OBS
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Altered: 27-Jul-2024 : Added Zn and Mn.
! Altered: 30-Dec-2023 : Added dE_OPTION for reading in SL assignments.
! Altered: 18-May-2015 : Changed GAM2, GAM4 to C4 and C6 (quadratic and Van der Waals
!                          interacton constants)(09-Jun-2105).
! Altered: 17-Mar-2003: SCRAT & SCRATREC are now initialized.
! Altered: 03-Mar-2000: Variable type ATM installed to simplify handling
!	                   of multiple species.
!
! Created:  5-Jan-1998=9 (Progran began late Dec, 1998)
!
	INTEGER ND		!Actula number of depth points in atmosphere
	INTEGER NC		!Actual number of core rays
	INTEGER NP		!Total number of rays (ND+NC)
!
	INTEGER ND_MAX,NP_MAX,NC_MAX
	INTEGER N_MAX
	INTEGER N_LINE_MAX
!
	CHARACTER(LEN=20) TIME
	CHARACTER(LEN=80) FILENAME,BLANK,DIR_NAME,STRING
	CHARACTER(LEN=20) FILE_EXTENT
	CHARACTER(LEN=11) FORMAT_DATE
	CHARACTER(LEN=10) NAME_CONVENTION
	CHARACTER(LEN=20) SL_OPTION		!Option to fudge SL assignments
	CHARACTER(LEN=20) dE_OPTION
!
	LOGICAL ASK           			!Ask of filenames or uset defaults.
	LOGICAL FILE_PRES
	LOGICAL SCRAT
	LOGICAL F_TO_S_RD_ERROR
	INTEGER I,J,NT,SCRATREC
	INTEGER N_SIG_UNDEF
	INTEGER LEN_DIR
	INTEGER EQ_TEMP
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) RMDOT
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER, PARAMETER :: T_IN=5		!Terminal IO
	INTEGER, PARAMETER :: T_OUT=6		!Terminal IO
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUMOD=8
	INTEGER, PARAMETER :: LU_TMP=10
	INTEGER, PARAMETER :: LU=30			!Used for pop file io
	INTEGER, PARAMETER :: LUSCRAT=33
!
	INTEGER NF
	INTEGER NS
	INTEGER ID
	INTEGER ISPEC
	LOGICAL PREV_STAGE_PRES
!
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=12), PARAMETER :: PRODATE='18-Jun-2025'
!
	DATA BLANK/' '/
!
! Set constants.
!
	CHIBF=2.815E-06_LDP
	CHIFF=3.69E-29_LDP
	HDKT=4.7994145_LDP
	TWOHCSQ=0.0147452575_LDP
	OPLIN=2.6540081E+08_LDP
	EMLIN=5.27296E-03_LDP
	OPLIN=2.6540081E+08_LDP
	EMLIN=5.27296E-03_LDP
	LUER=ERROR_LU()
	FILE_EXTENT=' '
!
! Open output file for all errors and comments.
!
        LUER=ERROR_LU(); IOS=0
        CALL GEN_ASCI_OPEN(LUER,'OUT_FLUX','UNKNOWN','APPEND',' ',IZERO,IOS)
        IF(IOS .NE. 0)THEN
          WRITE(LUER,*)'Error opening OUT_FLUX in CMF_FLUX_V5, IOS=',IOS
          STOP
        END IF
        CALL SET_LINE_BUFFERING(LUER)
!
!
! Set all atomic data. New species can be simple added by insertion.
! Try to add species in order of atomic number. Hydrogen should ALWAYS
! be species 1, Helium should ALWAYS be species 2.
!
! While this tabulation is Verbose, it is simple to change.
! Note that the Solar abundances are only used for reference in
! the output table MOD_SUM.
!
	ID=1
	AT_NO(1)=1.0_LDP;		    AT_MASS(ID)=1.0_LDP
	SPECIES(ID)='HYD';	    SPECIES_ABR(ID)='H'
	SOL_ABUND_HSCL(ID)=12.0_LDP
!
	ID=ID+1
	AT_NO(ID)=2.0_LDP; 	    AT_MASS(ID)=4.0_LDP		!Helium
	SPECIES(ID)='HE';	    SPECIES_ABR(ID)='He'
	SOL_ABUND_HSCL(ID)=11.0_LDP
!
	ID=ID+1
	AT_NO(ID)=6.0_LDP;	    AT_MASS(ID)=12.0_LDP		!Carbon
	SPECIES(ID)='CARB';	    SPECIES_ABR(ID)='C'
	SOL_ABUND_HSCL(ID)=8.56_LDP
!
	ID=ID+1
	AT_NO(ID)=7.0_LDP;	    AT_MASS(ID)=14.0_LDP		!Nitrogen
	SPECIES(ID)='NIT';	    SPECIES_ABR(ID)='N'
	SOL_ABUND_HSCL(ID)=8.05_LDP
!
	ID=ID+1
	AT_NO(ID)=8.0_LDP; 	    AT_MASS(ID)=16.0_LDP		!Oxygen
	SPECIES(ID)='OXY';	    SPECIES_ABR(ID)='O'
	SOL_ABUND_HSCL(ID)=8.93_LDP
!
	ID=ID+1
	AT_NO(ID)=9.0_LDP;	    AT_MASS(ID)=19.00_LDP		!Fluorine
	SPECIES(ID)='FLU';	    SPECIES_ABR(ID)='F'
	SOL_ABUND_HSCL(ID)=4.56_LDP
!
	ID=ID+1
	AT_NO(ID)=10.0_LDP;	    AT_MASS(ID)=20.2_LDP		!Neon
	SPECIES(ID)='NEON';	    SPECIES_ABR(ID)='Ne'
	SOL_ABUND_HSCL(ID)=8.09_LDP
!
	ID=ID+1
	AT_NO(ID)=11.0_LDP;	    AT_MASS(ID)=23.0_LDP		!Sodium
	SPECIES(ID)='SOD';	    SPECIES_ABR(ID)='Na'
	SOL_ABUND_HSCL(ID)=6.33_LDP
!
	ID=ID+1
	AT_NO(ID)=12.0_LDP;	    AT_MASS(ID)=24.3_LDP		!Magnesium
	SPECIES(ID)='MAG';	    SPECIES_ABR(ID)='Mg'
	SOL_ABUND_HSCL(ID)=7.58_LDP
!
	ID=ID+1
	AT_NO(ID)=13.0_LDP;	    AT_MASS(ID)=27.0_LDP		!Aluminium
	SPECIES(ID)='ALUM';	    SPECIES_ABR(ID)='Al'
	SOL_ABUND_HSCL(ID)=6.47_LDP
!
	ID=ID+1
	AT_NO(ID)=14.0_LDP;	    AT_MASS(ID)=28.1_LDP		!Silicon
	SPECIES(ID)='SIL';	    SPECIES_ABR(ID)='Sk'
	SOL_ABUND_HSCL(ID)=7.55_LDP
!
	ID=ID+1
	AT_NO(ID)=15.0_LDP;	    AT_MASS(ID)=31.0_LDP		!Phosphorous
	SPECIES(ID)='PHOS';	    SPECIES_ABR(ID)='P'
	SOL_ABUND_HSCL(ID)=5.45_LDP
!
	ID=ID+1
	AT_NO(ID)=16.0_LDP;	    AT_MASS(ID)=32.1_LDP		!Sulpher
	SPECIES(ID)='SUL';	    SPECIES_ABR(ID)='S'
	SOL_ABUND_HSCL(ID)=7.21_LDP
!
	ID=ID+1
	AT_NO(ID)=17.0_LDP;	    AT_MASS(ID)=35.5_LDP		!Chlorine
	SPECIES(ID)='CHL';	    SPECIES_ABR(ID)='Cl'
	SOL_ABUND_HSCL(ID)=5.5_LDP
!
	ID=ID+1
	AT_NO(ID)=18.0_LDP;	    AT_MASS(ID)=39.9_LDP		!Argon
	SPECIES(ID)='ARG';	    SPECIES_ABR(ID)='Ar'
	SOL_ABUND_HSCL(ID)=6.56_LDP
!
	ID=ID+1
	AT_NO(ID)=19.0_LDP;	    AT_MASS(ID)=39.1_LDP		!Potassium
	SPECIES(ID)='POT';	    SPECIES_ABR(ID)='K'
	SOL_ABUND_HSCL(ID)=5.12_LDP
!
	ID=ID+1
	AT_NO(ID)=20.0_LDP;	    AT_MASS(ID)=40.1_LDP		!Calcium
	SPECIES(ID)='CAL';	    SPECIES_ABR(ID)='Ca'
	SOL_ABUND_HSCL(ID)=6.36_LDP
!
        ID=ID+1
        AT_NO(ID)=21.0_LDP;           AT_MASS(ID)=44.96_LDP         !Scandium
        SPECIES(ID)='SCAN';         SPECIES_ABR(ID)='Sc'
        SOL_ABUND_HSCL(ID)=3.10_LDP
!
	ID=ID+1
	AT_NO(ID)=22.0_LDP;	    AT_MASS(ID)=47.88_LDP		!Titanium
	SPECIES(ID)='TIT';	    SPECIES_ABR(ID)='Tk'	!Actual symbol is Ti
	SOL_ABUND_HSCL(ID)=4.99_LDP
!
	ID=ID+1
	AT_NO(ID)=23.0_LDP;           AT_MASS(ID)=50.94_LDP         !Vandium
	SPECIES(ID)='VAN';          SPECIES_ABR(ID)='V'         !Actual symbol is V
	SOL_ABUND_HSCL(ID)=4.00_LDP
!
	ID=ID+1
	AT_NO(ID)=24.0_LDP;	    AT_MASS(ID)=52.0_LDP		!Chromium
	SPECIES(ID)='CHRO';	    SPECIES_ABR(ID)='Cr'
	SOL_ABUND_HSCL(ID)=5.67_LDP
!
	ID=ID+1
	AT_NO(ID)=25.0_LDP;	    AT_MASS(ID)=54.9_LDP		!Maganese
	SPECIES(ID)='MAN';	    SPECIES_ABR(ID)='Mn'
	SOL_ABUND_HSCL(ID)=5.39_LDP
!
	ID=ID+1
	AT_NO(ID)=26.0_LDP;	    AT_MASS(ID)=55.8_LDP		!Iron
	SPECIES(ID)='IRON';	    SPECIES_ABR(ID)='Fe'
	SOL_ABUND_HSCL(ID)=7.54_LDP
!
	ID=ID+1
	AT_NO(ID)=27.0_LDP;	    AT_MASS(ID)=58.9_LDP		!Cobalt
	SPECIES(ID)='COB';	    SPECIES_ABR(ID)='Co'
	SOL_ABUND_HSCL(ID)=4.92_LDP
!
	ID=ID+1
	AT_NO(ID)=28.0_LDP;	    AT_MASS(ID)=58.7_LDP		!Nickel
	SPECIES(ID)='NICK';	    SPECIES_ABR(ID)='Nk'
	SOL_ABUND_HSCL(ID)=6.25_LDP
!
	ID=ID+1
	AT_NO(ID)=29.0_LDP;         AT_MASS(ID)=63.6_LDP            !Copper
	SPECIES(ID)='COP';          SPECIES_ABR(ID)='Cu'
	SOL_ABUND_HSCL(ID)=4.20_LDP
!
	ID=ID+1
	AT_NO(ID)=30.0_LDP;         AT_MASS(ID)=65.40_LDP           !Zinz
	SPECIES(ID)='ZINC';         SPECIES_ABR(ID)='Zn'
	SOL_ABUND_HSCL(ID)=4.60_LDP
!
	ID=ID+1
	AT_NO(ID)=56.0_LDP;           AT_MASS(ID)=137.33_LDP        !Barium
	SPECIES(ID)='BAR';          SPECIES_ABR(ID)='Ba'
	SOL_ABUND_HSCL(ID)=2.13_LDP
!
	IF(ID .NE. NUM_SPECIES)THEN
	  WRITE(LUER,*)'Error in CMFGEN: Invalid species setup'
	  STOP
	END IF
!
! Convert from the abundance on a logarithmic scale with H=12.0 dex, to
! mass-fractions.
!
	T1=0.0_LDP
	DO ID=1,NUM_SPECIES
	  SOL_MASS_FRAC(ID)=10.0_LDP**(SOL_ABUND_HSCL(ID)-12.0_LDP)
	  SOL_MASS_FRAC(ID)=AT_MASS(ID)*SOL_MASS_FRAC(ID)
	  T1=T1+SOL_MASS_FRAC(ID)
	END DO
	SOL_MASS_FRAC(:)=SOL_MASS_FRAC(:)/T1
!
! Initilaization: These parameters will remain as they are when a species
! is not present.
!
	EQ_SPECIES(:)=0
	SPECIES_PRES(:)=.FALSE.
	AT_ABUND(:)=0.0_LDP
	ABUND_SCALE_FAC(:)=1.0_LDP
	DO ID=1,MAX_NUM_IONS
	  ATM(ID)%NXzV_F=1; 	ATM(ID)%NXzV=0
	  ATM(ID)%XzV_PRES=.FALSE.
	END DO
!
! Read in the gaunt factors for individual l states of hydrogen.
!
	CALL RD_HYD_BF_DATA(LUIN,LUMOD,T_OUT)
!
! *************************************************************************
!
! Read in basic model. This includes scaler quantities which describe the
! model (ie MDOT, ND, NC) as well as important vectors [i.e. Density
! structure, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
10	FILENAME=' '
	DIR_NAME=' '
	ASK=.FALSE.
	SCRAT=.FALSE.
	SCRATREC=1
	IOS=2			!Filename has to exits, blank not allowed.
	WRITE(T_OUT,*)'Program date is: ',PRODATE
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'Append (ask) to so that subsequent FILE names are not defaulted'
	WRITE(T_OUT,*)'Append (scrat) to get scratch output.'
	WRITE(T_OUT,*)' '
	I=80
	CALL RD_NCHAR(FILENAME,'RVTJ',I,T_IN,LUER,'Input RVTJ filename')
	STRING=FILENAME
	CALL SET_CASE_UP(STRING,0,0)
	IF( INDEX(STRING,'(SCRAT)') .NE. 0)SCRAT=.TRUE.
	IF( INDEX(STRING,'(ASK)') .NE. 0)ASK=.TRUE.
	I= INDEX(STRING,'(')
	IF(I .NE. 0)FILENAME=FILENAME(1:I-1)//' '
	INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	IF(.NOT. FILE_PRES)THEN
	  WRITE(LUER,*)'Unable to find file'
	  WRITE(LUER,*)'File is ',TRIM(FILENAME)
	  GOTO 10
	END IF
!
! Get default directory.
!
	DIR_NAME=' '		!Valid DIR_NAME if not present.
	LEN_DIR=0
	J=LEN(FILENAME)
	DO WHILE(J .GT. 0)
	  IF( FILENAME(J:J) .EQ. ']' .OR.
	1     FILENAME(J:J) .EQ. ':' .OR.
	1     FILENAME(J:J) .EQ. '/'        )THEN
	    DIR_NAME=FILENAME(1:J)
	    LEN_DIR=J
	    J=0
	  END IF
	  J=J-1
	END DO
!
! Get default extension.
!
	FILE_EXTENT=' '
	J=LEN(FILENAME)
	DO WHILE(J .GT. LEN_DIR)
	  IF( FILENAME(J:J) .EQ. '.' )THEN
	    FILE_EXTENT=FILENAME(J:)
	    J=0
	  END IF
	  J=J-1
	END DO
!
	CALL RD_RVTJ_PARAMS_V4(RMDOT,STARS_LUM,AT_ABUND(1),TIME,NAME_CONVENTION,
	1                    ND,NC,NP,FORMAT_DATE,FILENAME,LUIN)
	ALLOCATE (R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SIGMA(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SIGMA_TRUE(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (T(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ED(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (LANG_COORD(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CMFGEN_TGREY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (dE_RAD_DECAY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ROSS_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (FLUX_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (PLANCK_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POP_ATOM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (DENSITY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POPION(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CLUMP_FAC(ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX'
	  WRITE(LUER,*)'Unable to allocate Atmosphere arrays'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
	CALL RD_RVTJ_VEC_V4(R,V,SIGMA,ED,T,CMFGEN_TGREY,dE_RAD_DECAY,
	1       ROSS_MEAN,FLUX_MEAN,PLANCK_MEAN,
	1       POP_ATOM,POPION,DENSITY,CLUMP_FAC,FORMAT_DATE,ND,LUIN)
!
! RVTJ stores the infrastructure-clamped SIGMA (>= -0.999).  Re-derive
! the unclamped SIGMA_TRUE through the shared worker (SIGMA_TRUE_CORE,
! new_main/subs/recompute_sigma_true.f) so the EW estimators see the
! same gradient -- including the stagnation-cell regularization and
! degenerate-cell count -- as the transfer solve.  Called
! unconditionally: CONTROL_VARIABLE_MOD's nonmono flags are not set in
! this executable, only the CMF_FLUX_CNTRL_VAR_MOD copies are.
!
	CALL SIGMA_TRUE_CORE(SIGMA_TRUE,R,V,ND,N_SIG_UNDEF)
	IF(DO_NONMONO_SOB)THEN
	  WRITE(T_OUT,'(A,2ES12.4)')
	1    'SIGMA_TRUE range:',MINVAL(SIGMA_TRUE(1:ND)),
	2                        MAXVAL(SIGMA_TRUE(1:ND))
	  IF(N_SIG_UNDEF .GT. 0)WRITE(T_OUT,'(A,I6,A)')
	1    '  SIGMA_TRUE undefined at',N_SIG_UNDEF,' cells (set to 0)'
	END IF
!
	SN_MODEL=.FALSE.
	IF(V(ND) .GT. 100.0_LDP)SN_MODEL=.TRUE.
!
	CLOSE(UNIT=LUIN)
	WRITE(T_OUT,3)TIME
3	FORMAT(1X,'Model completed on ',A20)
!
! Convert to old naming convention if necessary.
!
	CALL SET_CASE_UP(NAME_CONVENTION,IZERO,IZERO)
	IF(NAME_CONVENTION .EQ. 'K_FOR_I')THEN
	ELSE IF(NAME_CONVENTION .EQ. 'X_FOR_I')THEN
	  DO I=1,NUM_SPECIES
	    IF(SPECIES_ABR(I)(2:2) .EQ. 'k')SPECIES_ABR(I)(2:2)='x'
	  END DO
	ELSE
	  WRITE(T_OUT,*)'Don''t recognize naming convention in DISPGEN'
	  WRITE(T_OUT,*)'NAME_CONVENTION= ',NAME_CONVENTION
	  STOP
	END IF
!
! 
! Open MODEL data file to get N_S, DIE_AUTO, and DIE_WI for each species.
!
	IF(ASK)THEN
	  IOS=1  		!File has to exist if filename input.
	  CALL GET_FILENAME(FILENAME,1,2,'MODEL','!',IOS)
	ELSE
	  FILENAME=' '
	  FILENAME=DIR_NAME(1:LEN_DIR)//'MODEL'//FILE_EXTENT
	  INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	  IF(.NOT. FILE_PRES)IOS=-1
	END IF
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODEL file'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
	CALL RD_MODEL_FILE(FILENAME,LU,IOS)
	CALL RD_DBLE(STARS_MASS,'MASS',T_IN,LUER,'Stars mass')
	IF(STARS_MASS .EQ. 0)CALL RD_MODEL_DBLE(STARS_MASS,'MASS')
	CALL RD_MODEL_LOG(DIE_AS_LINE,'DIE_AS_LINE')
	CLOSE(LU)
!
! Read in other parameters from batch file. These must be in order.
! Using a STRIN to read ONLY_OBS_LINES allows me to preserve backward
! compatibility, while adding the new variable ONLY_UNOBS)_LINES.
!
	CALL RD_CHAR(STRING,'ONLY_OBS_LINES',T_IN,LUER,'Observed lines only?')
	STRING=ADJUSTL(STRING)
	ONLY_OBS_LINES=.FALSE.; ONLY_UNOBS_LINES=.FALSE.
	IF(STRING(1:1) .EQ. 'T' .OR. STRING(1:2) .EQ. '.T')ONLY_OBS_LINES=.TRUE.
	IF(STRING(1:3) .EQ. 'OBS')ONLY_OBS_LINES=.TRUE.
	IF(STRING(1:5) .EQ. 'UNOBS')ONLY_UNOBS_LINES=.TRUE.
!
	FILENAME=TRIM(DIR_NAME)//'MODEL_SPEC'
	CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODEL_SPEC in CMF_FLUX_V5, IOS=',IOS
	  STOP
	END IF
	WRITE(LUER,'(A)')' '
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUER)
	SL_OPTION=' '; dE_OPTION=' ';  I=20
        CALL RD_STORE_NCHAR(SL_OPTION,'SL_OPT',I,L_FALSE,'Auto SL splitting option')
	CALL RD_STORE_NCHAR(dE_OPTION,'dE_OPT',I,L_FALSE,'Option to split non-LS SLs by enegry')
	WRITE(6,*)'CMF_FLUX_V5 is adopting the following SL_OPTION: ',TRIM(SL_OPTION)
	WRITE(6,*)'CMF_FLUX_V5 is adopting the following dE_OPTION: ',TRIM(dE_OPTION)
	WRITE(LUER,'(A)')' '
	CALL CLEAN_RD_STORE()
!
! 
!
!
! Allocate population vectors.
!
	ALLOCATE (POP_SPECIES(ND,NUM_SPECIES)); POP_SPECIES=0.0_LDP
!
	EQ_TEMP=1
	NT=2
	ID=1
	DO ISPEC=1,NUM_SPECIES
!
! Open input file for each species: This file is an ASCI (new format) of a direct
! access file (old format) and contains the ouput for the entrie atom.
! The file need not exist. We do this for all species including H and He.
!
	  IOS=1
	  DO WHILE(IOS .NE. 0)
	    IOS=0
	    IF(ASK .AND. IOS .EQ. 0)THEN
	      IOS=1  		!File has to exist if filename input.
	      CALL GET_FILENAME(FILENAME,1,2,TRIM(SPECIES(ISPEC)),'!',IOS)
	    ELSE
	      FILENAME=DIR_NAME(1:LEN_DIR)//'POP'//TRIM(SPECIES(ISPEC))//FILE_EXTENT
	      INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	      IF(.NOT. FILE_PRES)FILENAME=BLANK
	    END IF
!
! Read in basic model data (TIMECHK, and ABUNDC)
!
	    IF(FILENAME .NE. BLANK)THEN
	      CALL OP_SPEC_FILE_V2(FILENAME,LU,AT_ABUND(ISPEC),POP_SPECIES(1,ISPEC),
	1          ND,FORMAT_DATE,IOS,TIME,TRIM(SPECIES(ISPEC)))
	    END IF
	  END DO
!
	  PREV_STAGE_PRES=.FALSE.
	  IF(FILENAME .NE. BLANK)THEN
!
	    DO J=1,MAX_IONS_PER_SPECIES
	      ION_ID(ID)=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))
	      CALL RD_POP_DIM(ATM(ID)%NXzV_F, ATM(ID)%XzV_PRES,
	1                     TRIM(ION_ID(ID)),FORMAT_DATE,LU)
!
! Only allocate memory if ion is available. We need, however, to allocate
! for the ion also.
!
	      IF(ATM(ID)%XzV_PRES .OR. PREV_STAGE_PRES)THEN
	        CALL RD_MODEL_SPEC_INFO( TRIM(ION_ID(ID)),
	1          ATM(ID)%NXzV_F,       ATM(ID)%XzV_PRES,  ATM(ID)%NXzV,
	1          ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV)
                SPECIES_PRES(ISPEC)=.TRUE.
	        IF(SPECIES_BEG_ID(ISPEC) .EQ. 0)SPECIES_BEG_ID(ISPEC)=ID
	        IF(SL_OPTION .EQ. 'FULL')THEN
	          ATM(ID)%NXzV=ATM(ID)%NXzV_F
	        END IF
	        NF=ATM(ID)%NXzV_F
	        NS=ATM(ID)%NXzV
!
	                      ALLOCATE (ATM(ID)%XzV(NS,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE(NS,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%LOG_XzVLTE(NS,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%dlnXzVLTE_dlnT(NS,ND),STAT=IOS)
!
	  	IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzV_F(NF,ND),STAT=IOS)
	  	IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F(NF,ND),STAT=IOS)
	  	IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%LOG_XzVLTE_F(NF,ND),STAT=IOS)
	  	IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F_ON_S(NF,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%W_XzV_F(NF,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%AXzV_F(NF,NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%EDGEXzV_F(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%F_TO_S_XzV(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%INT_SEQ_XzV(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLEVNAME_F(NF),STAT=IOS)
!
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV_F(ND),STAT=IOS)       ; ATM(ID)%DXzV_F(:)=0.0_LDP
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV(ND),STAT=IOS)         ; ATM(ID)%DXzV(:)=0.0_LDP
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GXzV_F(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%ARAD(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%C4(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%C6(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%OBSERVED_LEVEL(NF),STAT=IOS)
	        IF(IOS .NE. 0)THEN
	          WRITE(LUER,*)'Error in CMF_FLUX'
	          WRITE(LUER,*)'Unable to allocate arrays for species XzV'
	          WRITE(LUER,*)'STATUS=',IOS
	          STOP
	        END IF
!
	        DO I=1,ATM(ID)%NXzV_F ; ATM(ID)%F_TO_S_XzV(I)=I ; END DO
	        CALL RD_ION_POP_V3(ATM(ID)%XzV_F, ATM(ID)%DXzV_F,
	1             ATM(ID)%XzV_PRES,    ATM(ID)%NXzV_F,
	1             ATM(ID)%XzV_OSCDATE, TRIM(ION_ID(ID)),
	1             FORMAT_DATE,LU,ND,SCRAT,LUSCRAT,SCRATREC)
	        SPECIES_LNK(ID)=ISPEC
	        PREV_STAGE_PRES=ATM(ID)%XzV_PRES
	        IF(.NOT. ATM(ID)%XZV_PRES)SPECIES_END_ID(ISPEC)=ID
	        ATM(ID)%EQXzV=EQ_TEMP
	        IF(.NOT. ATM(ID)%XzV_PRES)THEN
	          NT=NT+1
	          EQ_SPECIES(ISPEC)=EQ_TEMP
	          EQ_TEMP=EQ_TEMP+1
	        ELSE
	          NT=NT+ATM(ID)%NXzV
	          EQ_TEMP=EQ_TEMP+ATM(ID)%NXzV
	        END IF
	        ID=ID+1
	      END IF
	    END DO
	  END IF
	END DO			!Over species
!
	NUM_IONS=ID-1
!
	CLOSE(UNIT=LU)
	EQNE=NT-1
!
	CALL CLEAN_MODEL_STORE
!
! 
!
! Read options into store. Most of these will be read in CMF_FLUX_SUB_V5, but a few
! are needed in order to handle the photoioization data correct.
!
	CALL GEN_ASCI_OPEN(LUIN,'CMF_FLUX_PARAM','OLD',' ','READ',IZERO,IOS)
        IF(IOS .NE. 0)THEN
          WRITE(LUER,*)'Error opening CMF_FLUX_PARAM in CMFFLUX_SUB_V5, IOS=',IOS
          STOP
        END IF
        CALL GEN_ASCI_OPEN(LUMOD,'OUT_PARAMS','UNKNOWN',' ','WRITE',IZERO,IOS)
        CALL RD_OPTIONS_INTO_STORE(LUIN,LUMOD)
	CLOSE(LUIN)
!
! Ideally, VSM_DIE_KMS and SIG_GAU_KMS will have the same value.
!
	CALL RD_STORE_DBLE(VSM_DIE_KMS,'VSM_DIE',L_TRUE,
	1         'Sigma of Gaussian used for dielectronic lines that are read in')
	CALL RD_STORE_DBLE(SIG_GAU_KMS,'SIG_GAU_KMS',L_TRUE,
	1         'Sigma of Gaussian used to smooth photoionization data')
	FRAC_SIG_GAU=0.25_LDP
	CALL RD_STORE_DBLE(FRAC_SIG_GAU,'FRAC_SIG',L_FALSE,
	1         'Fractional spacing a across smoothing Gauusian (use 0.25)')
	CUT_ACCURACY=0.02_LDP
	CALL RD_STORE_DBLE(CUT_ACCURACY,'CUT_ACC',L_FALSE,
	1         'Accuracy to retain data when omitting data points to save space (use 0.02)')
	ABOVE_EDGE=.TRUE.
	CALL RD_STORE_LOG(ABOVE_EDGE,'ABV_EDGE',L_FALSE,
	1         'Use only data above edge when smoothing (TRUE)')
!
	GF_CUT=0.0_LDP
	GF_LEV_CUT=1000
	GF_ACTION=' '
	MIN_NUM_TRANS=1000
	XRAYS=.FALSE.
	T2=0.0_LDP
!
	F_TO_S_RD_ERROR=.FALSE.
	DO ID=NUM_IONS,1,-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    FILENAME=TRIM(ION_ID(ID))//'_F_OSCDAT'
	    CALL GENOSC_V9(ATM(ID)%AXzV_F,  ATM(ID)%EDGEXzV_F,
	1          ATM(ID)%GXzV_F, ATM(ID)%XzVLEVNAME_F,
	1          ATM(ID)%ARAD,ATM(ID)%C4,ATM(ID)%C6,
	1          ATM(ID)%OBSERVED_LEVEL,T1,ATM(ID)%ZXzV,
	1          ATM(ID)%NEW_XzV_OSCDATE, ATM(ID)%NXzV_F, I,
	1          GF_ACTION,GF_CUT,GF_LEV_CUT,MIN_NUM_TRANS,
	1          ONLY_OBS_LINES,ONLY_UNOBS_LINES,
	1          LUIN,LU_TMP,FILENAME)
	    IF(ATM(ID)%XzV_OSCDATE .NE. ATM(ID)%NEW_XzV_OSCDATE)THEN
	       WRITE(T_OUT,*)'Warning --- invalid date for ',ION_ID(ID)
	       WRITE(T_OUT,*)'Old oscilator date:',ATM(ID)%XzV_OSCDATE
	       WRITE(T_OUT,*)'New oscilator date:',ATM(ID)%NEW_XzV_OSCDATE
	    END IF
	    FILENAME=TRIM(ION_ID(ID))//'_AUTO_DATA'
	    CALL ADD_AUTO_RATES(ATM(ID)%ARAD,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,
	1              ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,FILENAME)
!
	    FILENAME=TRIM(ION_ID(ID))//'_F_TO_S'
	    CALL RD_F_TO_S_IDS_V4(ATM(ID)%F_TO_S_XzV,ATM(ID)%INT_SEQ_XzV,
	1           ATM(ID)%XzVLEVNAME_F,ATM(ID)%NXzV_F,ATM(ID)%NXzV,
	1           LUIN,FILENAME,SL_OPTION,dE_OPTION,F_TO_S_RD_ERROR)
	    CALL RDPHOT_GEN_V2(ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1            ATM(ID)%GIONXzV_F,      AT_NO(SPECIES_LNK(ID)),
	1            ATM(ID)%ZXzV,           ATM(ID)%NXzV_F,
	1            ATM(ID)%XzV_ION_LEV_ID, ATM(ID)%N_XzV_PHOT, NPHOT_MAX,
	1            ATM(ID+1)%XzV_PRES,     ATM(ID+1)%EDGEXzV_F,
	1            ATM(ID+1)%GXzV_F,       ATM(ID+1)%F_TO_S_XzV,
	1            ATM(ID+1)%XzVLEVNAME_F, ATM(ID)%NXzV_F,
	1            SIG_GAU_KMS,FRAC_SIG_GAU,CUT_ACCURACY,ABOVE_EDGE,
	1            XRAYS,ID,ION_ID(ID),LUIN,LU_TMP)
            IF(ATM(ID+1)%XzV_PRES)ATM(ID)%GIONXzV_F=ATM(ID+1)%GXzV_F(1)
 	    IF(.NOT. DIE_AS_LINE .AND. (ATM(ID)%DIE_AUTO_XzV .OR. ATM(ID)%DIE_WI_XzV) )THEN
	      FILENAME='DIE'//TRIM(ION_ID(ID))
	      CALL RD_PHOT_DIE_V1(ID,  ATM(ID)%EDGEXzV_F,    ATM(ID)%XzVLEVNAME_F,
	1             ATM(ID)%NXzV_F,  ATM(ID)%GIONXzV_F,
	1             VSM_DIE_KMS,     ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LU_TMP,FILENAME)
 	    ELSE IF(DIE_AS_LINE .AND. (ATM(ID)%DIE_AUTO_XzV .OR. ATM(ID)%DIE_WI_XzV) )THEN
	      WRITE(LUER,*)'Warning in CMF_FLUX_V5'
	      WRITE(LUER,*)'Dielectronic lines from files not included'//
	1                    ' in spectrum calculation.'
	    END IF
	  END IF
	END DO		!Over NUM_SPECIES
        IF(F_TO_S_RD_ERROR)THEN
          WRITE(6,*)' '
          WRITE(6,*)'There are errors reading in the super level links.'
          WRITE(6,*)'These need to be fixed before the code can run.'
          WRITE(6,*)'See F_TO_S_RD_ERRORS for details.'
          WRITE(6,*)' '
          STOP
        END IF
!
! This surbroutine allow forbidden lines to be omitted -- either for individual
! ionization stages, or for all species.
!
        CALL SET_FORBID_ZERO(LUIN)
!
! 
! Determine the total number of bound-bound transitions in the model.
!
	N_MAX=0
	N_LINE_MAX=0
	DO ID=1,NUM_IONS
	  N_MAX=MAX(N_MAX,ATM(ID)%NXzV_F)
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO I=1,ATM(ID)%NXzV_F-1
	      DO J=I+1,ATM(ID)%NXzV_F
	        IF(ATM(ID)%AXzV_F(I,J) .NE. 0)N_LINE_MAX=N_LINE_MAX+1
	      END DO
	    END DO
	  END IF
	END DO
!
	ND_MAX=4*ND-1
	NC_MAX=4*NC
	NP_MAX=ND_MAX+NC_MAX
	FLUSH(LUER)
!
!*****************************************************************************
!*****************************************************************************
!
! Call program to compute emissivities, opacities, and the mean intensities
! in the CMF. The routine then passes (calls) OBS_FRAM_SUb to compute the
! model spectrum
!
	CALL TUNE(1,'CMF_FLUX_SUB')
	CALL CMF_FLUX_SUB_V5(ND,NC,NP,ND_MAX,NP_MAX,NT,N_LINE_MAX)
	CALL TUNE(2,'CMF_FLUX_SUB')
!
	CALL TUNE(3,' ')
!
	WRITE(LUER,'(//,A,//)')' CMF_FLUX has finished'
	STOP
!
	END
