!
! Program to compute the emergent flux in the observer's frame using an
! Observer's frame flux calculation. The program first computes the
! emissivities, opacities, and radiation field in the CMF frame. COHERENT
! or INCOHERENT electron scattering can be assumed. This data is then passed
! to OBS_FRAME_SUB for computation of the OBSERVER's flux.
!
! CMF_FLUX calling program is partially based on DISPGEN.
!
	PROGRAM CMFGEN
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Altered 25-Jul-2024: Added new species -- Zn and Mn.
! Altered 02-Sep-2023: Now zero ABS_MEAN
! Altered 24-Sep-2023: Updated phsical and opacity constants.
! Altered 20-Jan-2023: ABS_MEAN added. The code now checks that star models go to
!                           sufficently high optical depth. It also issues a
!                           warning when there is poor flux convergence.
! Altered 15-Jul-2019 (osiris - added ATM(ID)%dWCRXzVdT allocation).
!                        Incoprated IBIS -- 17-Aug-2019.
! Icorporated 02-Jun-2014: Changes to allow depth dependent profiles.
! Altered: 29-Nov-2011: Memory allocation for OLD_LEV_POP_AVAIL added.
! Altered: 25-Sep-2011: ION_ID now set for last ioization stage (LAST_ION).
!                          LOG_XzVLTE, LOG_XzVLTE_F and XzVLTE_F_ON_S allocated.
! Altered: 19-Jan-2008: SL_OPTION installed. Allows number of SL's to be
!                         automatically changed.
! Altered: 20-Mar-2005: FLUX_MEAN & ROSS_MEAN now zeroed. These vectors now
!                         used in CMFGEN_SUB.
! Altered: 03-Mar-2000: Variable type ATM installed to simplify handling
!	                   of multiple species.
!
! Created:  5-Jan-1998=9 (Program began late Dec, 1998)
!
	INTEGER ND		!Actual number of depth points in atmosphere
	INTEGER NC		!Actual number of core rays
	INTEGER NP		!Total number of rays (ND+NC)
	INTEGER NUM_BNDS	!Number of bans in linearization matrix
	INTEGER MAX_SIM	!Maximum number of lines that can be treated sim.
	INTEGER NCF_MAX	!Maximum number of frequencies that can be treated.
!
	INTEGER ND_MAX,NP_MAX
	INTEGER N_LINE_MAX
!
	INTEGER NM
	INTEGER NLF
	INTEGER NM_KI
	INTEGER TX_OFFSET
	INTEGER NION
	INTEGER DIAG_INDX
!
	CHARACTER*20 TEMP_KEY
!
	REAL(KIND=LDP) T1		!Temporary variable
	INTEGER I,J,IOS,NT
	INTEGER EQ_TEMP
	INTEGER LAST_ION_J
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LU_IN=7
	INTEGER, PARAMETER :: LU_OUT=8
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	INTEGER NF
	INTEGER NS
	INTEGER NV
	INTEGER ID
	INTEGER ISPEC
	INTEGER NUM_IONS_RD
!
	INTEGER  LUER,LUWARN
	INTEGER  ERROR_LU,WARNING_LU
	REAL(KIND=LDP)   PLANCKS_CONSTANT, FUN_PI, SPEED_OF_LIGHT, BOLTZMANN_CONSTANT
	REAL(KIND=LDP)   ELECTRON_MASS, ELECTRON_CHARGE
	EXTERNAL ERROR_LU,WARNING_LU, PLANCKS_CONSTANT, FUN_PI
	EXTERNAL SPEED_OF_LIGHT, BOLTZMANN_CONSTANT, ELECTRON_CHARGE, ELECTRON_MASS
!
	LOGICAL AT_LEAST_ONE_ION_PRES
	LOGICAL FND_END_OF_IONS
	LOGICAL DO_TERM_OUT
	CHARACTER(LEN=20) TIME
!
! Open output files for all errors and comments. Change DO_TERM_OUT to
! have the output go to the terminal/batch log file. NB: The WARNINGS
! file is overwritten.
!
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
	DO_TERM_OUT=.FALSE.                                   !TRUE.
	IF(.NOT. DO_TERM_OUT)THEN
	  CALL GEN_ASCI_OPEN(LUER,'OUTGEN','UNKNOWN','APPEND',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening OUTGEN in CMFGEN, IOS=',IOS
	    STOP
	  END IF
	  CALL SET_LINE_BUFFERING(LUER)
	END IF
	CALL GEN_ASCI_OPEN(LUWARN,'WARNINGS','UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUWARN,*)'Error opening WARNINGS in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL SET_LINE_BUFFERING(LUWARN)
        CALL DATE_TIME(TIME)
        WRITE(LUER,'(//,'' Model started on:'',15X,(A))')TIME
!
! Set constants.
!
	CHIBF=2.815E-06_LDP
	CHIFF=3.69E-29_LDP
	HDKT=1.0E+11_LDP*PLANCKS_CONSTANT()/BOLTZMANN_CONSTANT()                                  !Old value: HDKT=4.7994145D0
	TWOHCSQ=2.0E45_LDP*PLANCKS_CONSTANT()/(SPEED_OF_LIGHT())**2                               !Old value: TWOHCSQ=0.0147452575D0
	OPLIN=1.0E+10_LDP*FUN_PI()*(ELECTRON_CHARGE())**2/ELECTRON_MASS()/SPEED_OF_LIGHT()        !Old value  OPLIN=2.6540081D+08
	EMLIN=1.0E+25_LDP*PLANCKS_CONSTANT()/4.0_LDP/FUN_PI()                                       !Old value: EMLIN=5.27296D-03
!
	WRITE(LUER,*)' '
	WRITE(LUER,'(1X,A)')'Opacity/excitaion parameters adopted in CMFGEN are:'
	WRITE(LUER,'(30X,A,ES16.8)')'     CHIBF=',CHIBF
	WRITE(LUER,'(30X,A,ES16.8)')'     CHIFF=',CHIFF
	WRITE(LUER,'(30X,A,ES16.8)')'      HDKT=',HDKT
	WRITE(LUER,'(30X,A,ES16.8)')' TWOHONCSQ=',TWOHCSQ
	WRITE(LUER,'(30X,A,ES16.8)')'     OPLIN=',OPLIN
	WRITE(LUER,'(30X,A,ES16.8)')'     EMLIN=',EMLIN
	WRITE(LUER,*)' '
	WRITE(LUER,'(A70,I5)')'The KIND of the selected floating'//
	1    ' point format is:',KIND(HDKT)
	WRITE(LUER,'(A70,I5)')'The number of storage bits in the selected floating'//
	1    ' point format is:',STORAGE_SIZE(HDKT)
	WRITE(LUER,'(A70,I5)')'The number of significant digits (precision) is:',PRECISION(HDKT)
	WRITE(LUER,'(A70,I5)')'The maximum exponent range is:',RANGE(HDKT)
	WRITE(LUER,*)' '
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
	AT_NO(ID)=9.0_LDP;            AT_MASS(ID)=19.00_LDP         !Fluorine
	SPECIES(ID)='FLU';          SPECIES_ABR(ID)='F'
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
	AT_NO(ID)=21.0_LDP;	    AT_MASS(ID)=44.96_LDP		!Scandium
	SPECIES(ID)='SCAN';	    SPECIES_ABR(ID)='Sc'
	SOL_ABUND_HSCL(ID)=3.10_LDP
!
	ID=ID+1
	AT_NO(ID)=22.0_LDP;	    AT_MASS(ID)=47.88_LDP		!Titanium
	SPECIES(ID)='TIT';	    SPECIES_ABR(ID)='Tk'	!Actual symbol is Ti
	SOL_ABUND_HSCL(ID)=4.99_LDP
!
	ID=ID+1
	AT_NO(ID)=23.0_LDP;	    AT_MASS(ID)=50.94_LDP		!Vandium
	SPECIES(ID)='VAN';	    SPECIES_ABR(ID)='V'		!Actual symbol is V
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
	AT_NO(ID)=29.0_LDP;         AT_MASS(ID)=63.6_LDP                 !Copper
	SPECIES(ID)='COP';          SPECIES_ABR(ID)='Cu'
	SOL_ABUND_HSCL(ID)=4.20_LDP
!
	ID=ID+1
	AT_NO(ID)=30.0_LDP;         AT_MASS(ID)=65.40_LDP       	 !Zinc
	SPECIES(ID)='ZINC';         SPECIES_ABR(ID)='Zn'
	SOL_ABUND_HSCL(ID)=4.60_LDP
!
	ID=ID+1
	AT_NO(ID)=56.0_LDP;	    AT_MASS(ID)=137.33_LDP	!Barium
	SPECIES(ID)='BAR';	    SPECIES_ABR(ID)='Ba'
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
	SPECIES_BEG_ID(:)=0
	SPECIES_END_ID(:)=-1
!
! Although these will be set later, we set them here to allow us to check
! that the included ionization stages are sequential (for each species).
!
	DO ID=1,MAX_NUM_IONS
	  ATM(ID)%NXzV_F=1; 	ATM(ID)%NXzV=1
	  ATM(ID)%XzV_PRES=.FALSE.
	END DO
!
! Get data describing number of depth points, number of atomic levels
! etc.
!
	OPEN(UNIT=LU_IN,FILE='MODEL_SPEC',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to open MODEL_SPEC in cmfgen.f'
	  STOP
	END IF
	CALL RD_OPTIONS_INTO_STORE(LU_IN,LU_OUT)
	CALL RD_STORE_INT(ND,'ND',L_TRUE,'Number of depth points')
	CALL RD_STORE_INT(NC,'NC',L_TRUE,'Number of core rays')
	CALL RD_STORE_INT(NP,'NP',L_TRUE,'Number of impact parameters')
	CALL RD_STORE_INT(NUM_BNDS,'NUM_BNDS',L_TRUE,
	1        'Number of bands in linearization matrix (BA)')
	IF(NUM_BNDS .NE. 1 .AND. NUM_BNDS .NE. 3)THEN
	  WRITE(LUER,*)'Error: NUM_BNDS in MODEL_SPEC must be 1 or 3. Value is',NUM_BNDS
	  STOP
	END IF
	CALL RD_STORE_INT(MAX_SIM,'MAX_SIM',L_TRUE,
	1        'Maximum # of lines that cab treated simultaneously')
	CALL RD_STORE_INT(N_LINE_MAX,'NLINE_MAX',L_TRUE,
	1        'Maximum # of lines that can be treated')
	CALL RD_STORE_INT(NCF_MAX,'NCF_MAX',L_TRUE,
	1        'Maximum # of frequencies that can be treated')
	CALL RD_STORE_INT(NLF,'NLF',L_TRUE,
	1        'Number of frequencies per Doppler profile in CMF mode (21)')
	FL_OPTION=' '; I=20
	CALL RD_STORE_NCHAR(FL_OPTION,'FL_OPT',I,L_FALSE,
	1        'Option to automatically fudge number of FLs in all species')
	SL_OPTION=' '; I=20
	CALL RD_STORE_NCHAR(SL_OPTION,'SL_OPT',I,L_FALSE,
	1        'Option to automatically fudge number of SLs in all species')
	dE_OPTION=' '; I=20
	CALL RD_STORE_NCHAR(dE_OPTION,'dE_OPT',I,L_FALSE,
	1        'Option to split non-LS SLs by enegry')
	IL_OPTION=' '; I=20
	CALL RD_STORE_NCHAR(IL_OPTION,'IL_OPT',I,L_FALSE,
	1        'Option to automatically fudge number of ILs in all species')
!
! We now get the number of atomic levels. Old MODEL_SPEC files, with NSF in the
! keyword can be read. ISF take's precident over NSF, and there is no check that
! there is not both an effectively identical NSF and ISF keyowrd. In this
! case, the number of important variables is assumed to be the same as NS.
!
	ID=0
	NUM_IONS_RD=0
	EQ_TEMP=1
	DO ISPEC=1,NUM_SPECIES
	  AT_LEAST_ONE_ION_PRES=.FALSE.
	  FND_END_OF_IONS=.FALSE.
	  DO J=1,MIN(MAX_IONS_PER_SPECIES,NINT(AT_NO(ISPEC))+1)
	    TEMP_KEY=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))//'_ISF'
	    NS=0; NV=0; NF=0
	    CALL RD_STORE_3INT(NV,NS,NF,TEMP_KEY,L_FALSE,' ')
	    IF(NS .EQ .0)THEN
	      TEMP_KEY=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))//'_NSF'
	      CALL RD_STORE_2INT(NS,NF,TEMP_KEY,L_FALSE,' ')
	      NV=NS
	    END IF
	    IF(NV .GT. NS)THEN
	      NV=NS
	      WRITE(LUER,*)'Warning NV should be less than, or equal to, NS -- resetting NV'
	      WRITE(LUER,'(1X,A,I5,4X,A,I5)')'NV=',NV,'NS=',NS
	      WRITE(LUER,*)'Species is ',TRIM(TEMP_KEY)
	    END IF
	    IF(NS .GT. NF)THEN
	      WRITE(LUER,*)'Error NS should be less than, or equal to, NF'
	      WRITE(LUER,'(1X,A,I5,4X,A,I5)')'NS=',NS,'NF=',NF
	      WRITE(LUER,*)'Species is ',TRIM(TEMP_KEY)
	      STOP
	    END IF
!
	    IF(NS .NE. 0)THEN
	      TEMP_KEY=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))//'_F_TO_S'
	      CALL FDG_F_TO_S_NS_V2(NF,NS,NV,FL_OPTION,SL_OPTION,dE_OPTION,IL_OPTION,LU_IN,TEMP_KEY)
	      ID=ID+1
	      ION_ID(ID)=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))
	      ATM(ID)%XzV_PRES=.TRUE.
	      ATM(ID)%NXzV_F=NF
	      ATM(ID)%NXzV=NS
	      ATM(ID)%NXzV_IV=NV
	      AT_LEAST_ONE_ION_PRES=.TRUE.
	      IF(SPECIES_BEG_ID(ISPEC) .EQ. 0)SPECIES_BEG_ID(ISPEC)=ID
	      ATM(ID)%EQXzV=EQ_TEMP
	      EQ_TEMP=EQ_TEMP+ATM(ID)%NXzV
	      SPECIES_LNK(ID)=ISPEC
	      IF(FND_END_OF_IONS)THEN
	        WRITE(LUER,*)'Error in CMFGEN'
	        WRITE(LUER,*)'Ionization stages for species ',
	1                      TRIM(SPECIES_ABR(ISPEC)),' are not consecutive'
	        WRITE(LUER,*)'Check your MODEL_SPEC file for typo''s'
	        STOP
	      END IF
	      LAST_ION_J=J+1
	    ELSE IF(AT_LEAST_ONE_ION_PRES)THEN
	      FND_END_OF_IONS=.TRUE.
	    END IF
	  END DO
!
! We now set the variables for the SINGLE level of the HIGHEST ionization stage.
! XzV_PRES must be false.
!
	  IF(AT_LEAST_ONE_ION_PRES)THEN
	    ID=ID+1
            SPECIES_PRES(ISPEC)=.TRUE.
	    SPECIES_END_ID(ISPEC)=ID
	    SPECIES_LNK(ID)=ISPEC
	    ION_ID(ID)=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(LAST_ION_J))
	    ATM(ID)%XzV_PRES=.FALSE.
	    ATM(ID)%EQXzV=EQ_TEMP
	    EQ_SPECIES(ISPEC)=EQ_TEMP
	    EQ_TEMP=EQ_TEMP+1
	    NUM_IONS_RD=NUM_IONS_RD+1		!i.e. number of species
	  END IF
	END DO
	NUM_IONS=ID
	NUM_IONS_RD=NUM_IONS-NUM_IONS_RD	!
	EQNE=EQ_TEMP
	NT=EQNE+1
!
! Check whether any ION type'os in MODEL_SPEC file --- i.e. check that all ions
! specified by [XzV_NSF] are valid.
!
	CALL CNT_NUM_KEYS(I,'_ISF]')
	CALL CNT_NUM_KEYS(J,'_NSF]')
	IF(I+J .NE. NUM_IONS_RD)THEN
	  WRITE(LUER,*)'Error in CMFGEN'
	  WRITE(LUER,*)'You have an invalid [XzV_ISF] or [XZV_NsF] key in MODEL_SPEC'
	  WRITE(LUER,*)'NUM_IONS=',NUM_IONS
	  WRITE(LUER,*)'NUM_IONS_RD=',NUM_IONS_RD
	  WRITE(LUER,*)'No IONS in MODEL_SPEC',I
	  STOP
	END IF

!
	CALL CLEAN_RD_STORE()
	CLOSE(LU_IN)
!
!
	IF(NP .NE. ND+NC .AND. NP .NE. ND+NC-2)THEN
	  WRITE(LUER,*)'Error in CMFGEN'
	  WRITE(LUER,*)'Invalid NP: should be ND+NC or ND+NC-2'
	  STOP
	END IF
!
! Check whether EQUATION LABELLING is consistent. ' I ' is used as the
! number of the current equation. We also set the variable SPEC_PRES which
! indicates whether at least one ionization stage of a species is present.
! It is used to determine, for example,  whether a number conservation
! equation is required.
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
	IF(EQNE .NE. I)THEN
	  WRITE(LUER,*)'Error - EQNE has wrong value in CMFGEN'
	  STOP
	END IF
	IF(NT .NE. I+1)THEN
	  WRITE(LUER,*)'Error - NT has wrong value in CMFGEN'
	  STOP
	END IF
!
!
	ALLOCATE (R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SIGMA(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SIGMA_TRUE(ND),STAT=IOS); SIGMA_TRUE=0.0_LDP
	IF(IOS .EQ. 0)ALLOCATE (LANG_COORD(ND),STAT=IOS); LANG_COORD=0.0_LDP
	IF(IOS .EQ. 0)ALLOCATE (T(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ED(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ROSS_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (PLANCK_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ABS_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (FLUX_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POP_ATOM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (DENSITY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POPION(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CLUMP_FAC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (VOL_EXP_FAC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (VTURB_VEC(ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX'
	  WRITE(LUER,*)'Unable to allocate Atmosphere arrays'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
	FLUX_MEAN(:)=0.0_LDP
	ROSS_MEAN(:)=0.0_LDP
	PLANCK_MEAN(:)=0.0_LDP
	ABS_MEAN(:)=0.0_LDP
	VOL_EXP_FAC(:)=0.0_LDP
!
! Allocate population vectors.
!
	ALLOCATE (POP_SPECIES(ND,NUM_SPECIES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (GAM_SPECIES(ND,NUM_SPECIES),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMFGEN'
	  WRITE(LUER,*)'Unable to allocate NUM_SPECIES and GAM_SPECIES arrays'
	  WRITE(LUER,*)'STATUS=',IOS
	  WRITE(LUER,*)ND,NUM_SPECIES
	  STOP
	END IF
	POP_SPECIES=0.0_LDP; GAM_SPECIES=0.0_LDP
!
! Allocate ATM memory
!
	DO ID=1,NUM_IONS
!
	  NF=ATM(ID)%NXzV_F
	  NS=ATM(ID)%NXzV
	                ALLOCATE (ATM(ID)%XzV(NS,ND),STAT=IOS);            ATM(ID)%XzV=0.0_LDP
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%LOG_XzVLTE(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%dlnXzVLTE_dlnT(NS,ND),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzV_F(NF,ND),STAT=IOS);          ATM(ID)%XzV_F=0.0_LDP
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F(NF,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%LOG_XzVLTE_F(NF,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F_ON_S(NF,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%W_XzV_F(NF,ND),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%AXzV_F(NF,NF),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%EDGEXzV_F(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%F_TO_S_XzV(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%INT_SEQ_XzV(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLEVNAME_F(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GXzV_F(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%ARAD(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM2(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM4(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%OBSERVED_LEVEL(NF),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV_F(ND),STAT=IOS)       ; ATM(ID)%DXzV_F(:)=0.0_LDP
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV(ND),STAT=IOS)         ; ATM(ID)%DXzV(:)=0.0_LDP
!
! These 4 arrays are now allocated by a sparate call after the
! photoioizations cross-sections have been read in.
!
!	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WSXzV(NS,ND,NPHOT_MAX),STAT=IOS)
!	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WCRXzV(NS,ND,NPHOT_MAX),STAT=IOS)
!	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%dWSXzVdT(NS,ND,NPHOT_MAX),STAT=IOS)
!	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%dWCRXzVdT(NS,ND,NPHOT_MAX),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WSE_X_XzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WCR_X_XzV(NS,ND),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%APRXzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%ARRXzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%BFCRXzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%FFXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CPRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CRRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CHG_PRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CHG_RRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%COOLXzV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%NTCXzV(ND),STAT=IOS)       ; ATM(ID)%NTCXzV(:)=0.0_LDP
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%NTIXzV(ND),STAT=IOS)       ; ATM(ID)%NTIXzV(:)=0.0_LDP
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%NTIXzV_2E(ND),STAT=IOS)    ; ATM(ID)%NTIXzV_2E(:)=0.0_LDP
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%NT_ION_CXzV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%NT_EXC_CXzV(ND),STAT=IOS)
!
!
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in CMF_FLUX'
	    WRITE(LUER,*)'Unable to allocate arrays for species XzV'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
	END DO
!
	I=0
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    I=I+1
	    ATM(ID)%INDX_XzV=I
	    ATM(ID)%EQXzVBAL=I
	  ELSE
	    ATM(ID)%INDX_XzV=0
	    ATM(ID)%EQXzVBAL=0
	  END IF
	END DO
!
! Set non-thermal scale-factors to unity.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    ATM(ID)%CROSEC_NTFAC=1.0_LDP
	    ATM(ID)%ION_CROSEC_NTFAC=1.0_LDP
	  END IF
	END DO
!
! Create a vector, IMP_VAR, to inidicate which variables are considered as
! important for the atmospheric structure.
!
	ALLOCATE (IMP_VAR(NT))
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	     CALL SET_IMP_VEC(IMP_VAR,ATM(ID)%NXzV,ATM(ID)%NXzV_IV,
	1                        ATM(ID)%EQXzV,NT,ATM(ID)%XzV_PRES)
	  END IF
	END DO
!
! Allocated logical vector to idicate whether a level was in available
! in the earlier model of a time sequence.
!
	ALLOCATE (OLD_LEV_POP_AVAIL(NT))
!
! Used when computing dCHI and dETA for CMF mode. Required dimension is 4.
!
	NM_KI=4
!
! Indicates what storage locations in dzNET can be used for lines.
! The first 2 locations are for the current total opacity and emissivity.
! The next 3 locations are for the continuum opacity, emissivity, and
! electron scattering opacity.
! levels can be used for
!
	TX_OFFSET=6
!
! Total number of storage locations to be set aside. We multiply MAX_SIM by
! 2 to account for both the upper and lower levels.
!
	NM=TX_OFFSET+2*MAX_SIM
!
! NION and NUM_IONS can be used interchangably. NION is passed to CMFGEN_SUB
! for dynamic array allocation.
!
	NION=NUM_IONS
	DIAG_INDX=NUM_BNDS/2+1
!
! Set arrays asside for temoray storage, and in an ACCURATE J calculation
! is to be performed.
!
	ND_MAX=MAX(NT,2*ND)
	NP_MAX=ND_MAX+2*NC
!
!       CALL INIT_PROF_MODULE(ND,NLINES_PROF_STORE,NFREQ_PROF_STORE)
!
	CALL CMFGEN_SUB(ND,NC,NP,NT,
	1            NUM_BNDS,NION,DIAG_INDX,
	1            ND_MAX,NP_MAX,NCF_MAX,N_LINE_MAX,
	1            TX_OFFSET,MAX_SIM,NM,NM_KI,NLF)
!
	CALL TUNE(3,' ')
	STOP
!
	END
