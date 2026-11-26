	MODULE MOD_GUESS_DC
	USE SET_KIND_MODULE
	INTEGER NUM_FILES
	INTEGER ID
!
! Altered 18-Oct-2021 - Added depth index to DC ouput. More output  precision for R.
! Altered 20-May-2019 - Bug fixes
! Created 02-Apr-2019
!
! Needed when reading EDDFACTOR
!
	INTEGER ND,NCF			
	REAL(KIND=LDP), ALLOCATABLE :: RJ(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NU(:)
!
! Needed when reading RVTJ.
!
	INTEGER ND_ATM,NC_ATM,NP_ATM
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: T(:)
	REAL(KIND=LDP), ALLOCATABLE :: ION_POP(:)
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)
	REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: POP_ATOM(:)
	REAL(KIND=LDP), ALLOCATABLE :: MASS_DENSITY(:)
	REAL(KIND=LDP), ALLOCATABLE :: POPION(:)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: PHOT_SUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: RECOM_SUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: GS_DC(:)		!GS denotes GS.
	REAL(KIND=LDP), ALLOCATABLE :: DC(:,:)		!Departure coefficents
	REAL(KIND=LDP), ALLOCATABLE :: T_EXC(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: GS_ION_POP(:)	!Ground-state ion population.
	REAL(KIND=LDP), ALLOCATABLE :: DC_RUB(:)
	REAL(KIND=LDP), ALLOCATABLE :: WORK(:)
	SAVE

	END MODULE MOD_GUESS_DC
!
! General subroutine to guess the departure coefficients for a new species or
! ionization stage. The following files are required:
!
!     EDDFACTOR
!     RVTJ
!     XzV_F_OSCDAT
!     XzIV_IN          (needed for added high ionization species only)
!
! This routine will only be used when the *_IN file does not exist.
!
	SUBROUTINE SUB_GUESS_DC(SPECIES,REF_SPECIES,GION_SPEC,GION_REF)
	USE SET_KIND_MODULE
	USE MOD_GUESS_DC
	IMPLICIT NONE
!
! Altered 20-May-2019 - Bug fixes
! Created 02-Apr-2019
!
	REAL(KIND=LDP) GION_SPEC
	REAL(KIND=LDP) GION_REF
	CHARACTER(LEN=6) SPECIES		!Ionization stage we're adding
	CHARACTER(LEN=6) REF_SPECIES		!Previous ionization staged when adding new species.
	LOGICAL, SAVE :: FIRST=.TRUE.
!
! Used when reading file containing energy levels, oscillator strengths, etc
! (e.g. CV_F_OSCDAT)
!
	INTEGER, PARAMETER :: N_MAX=5000
	CHARACTER(LEN=30) NAME(N_MAX)
	REAL(KIND=LDP) FEDGE(N_MAX)
	REAL(KIND=LDP) ENERGY(N_MAX)
	REAL(KIND=LDP) G(N_MAX)
	REAL(KIND=LDP) ION_EN
	REAL(KIND=LDP) ZION
	CHARACTER(LEN=30) EN_DATE
	INTEGER NLEV
!
        INTEGER ACCESS_F
        INTEGER, PARAMETER :: EDD_CONT_REC=3
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
! These can come from a similar model but the frequency range must be sufficient.
! For high ionization species, the model should be identical, since the departure
! coefficients are very sensitive to the electron temperature.  Program estimates
! the ground state departure coefficient, and use the same excitation temperature
! for all other levels. If adding a whole new species, start with lowest ionization
! species.
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_IN=5		!For terminal input
	INTEGER, PARAMETER :: T_OUT=6           !For terminal output
	INTEGER, PARAMETER :: LU_IN=10		!For file I/O
	INTEGER, PARAMETER :: LU_OUT=11
	INTEGER, PARAMETER :: LU_HEAD=12
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	REAL(KIND=LDP) RMDOT
	REAL(KIND=LDP) RLUM
	REAL(KIND=LDP) ABUND_HYD
	INTEGER GET_INDX_DP
!
	CHARACTER(LEN=11)  FORMAT_DATE
	CHARACTER(LEN=21)  TIME
	CHARACTER(LEN=10)  DATA_TYPE
	CHARACTER(LEN=40)  FILE_DATE
	CHARACTER(LEN=80)  FILENAME
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=10)  NAME_CONVENTION
!
	CHARACTER(LEN=80) RVTJ_FILE_NAME
	CHARACTER(LEN=12) IN_FILE
	CHARACTER(LEN=10) OUT_FILE
!
	REAL(KIND=LDP) SPEED_OF_LIGHT
!
! Miscellaneous variables.
!
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,ML
	INTEGER EDGE_ML
	INTEGER ST_REC
	INTEGER REC_LENGTH
	INTEGER ND_RD,NLEV_RD
	REAL(KIND=LDP) NU1,NU2
	REAL(KIND=LDP) RJ1,RJ2
	REAL(KIND=LDP) T1,T2,T3,T4
!
! 
!
	IF(FIRST)THEN
          CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
!  Read in EDDFACTOR file. This is used to compute the photoionization and
!         recombination rates.
!
	  FILENAME='EDDFACTOR'
	  CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening/reading INFO file: check format'
	    WRITE(T_OUT,*)'Also check error file or fort.2'
	    STOP
	  END IF
	  OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error opening ',TRIM(FILENAME)
	       WRITE(T_OUT,*)'IOS=',IOS
	      STOP
	    END IF
	    READ(LU_IN,REC=3)ST_REC,NCF,ND
	    ND=ND; NCF=NCF
	    ALLOCATE (RJ(ND,NCF))
	    ALLOCATE (NU(NCF))
	    DO ML=1,NCF
	      READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(RJ(I,ML),I=1,ND),NU(ML)
	      IF(IOS .NE. 0)THEN
	        WRITE(T_OUT,*)'Error reading all frequencies'
	        NCF=ML-1
	        EXIT
	      END IF
	    END DO
	  CLOSE(LU_IN)
	  WRITE(T_OUT,*)'Successfully read in ',TRIM(FILENAME),' file as MODEL A (default)'
	  WRITE(T_OUT,*)'Number of depth points is',ND
	  WRITE(T_OUT,*)'Number of frequencies is ',NCF
	  WRITE(T_OUT,*)' '
	  FLUSH(T_OUT)
!
!
!
! *************************************************************************
!
! Read in basic model [i.e. R, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
10	  RVTJ_FILE_NAME='RVTJ'
	  OPEN(UNIT=LU_IN,FILE=RVTJ_FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to open RVTJ in new_main/subs/sub_guess_dc.f: IOS=',IOS
	   STOP
	  END IF
	  CLOSE(LU_IN)
	  CALL RD_RVTJ_PARAMS_V4(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1             ND_ATM,NC_ATM,NP_ATM,FORMAT_DATE,RVTJ_FILE_NAME,LU_IN)
	  WRITE(T_OUT,*)'Successfully read RVTJ params'
	  FLUSH(T_OUT)
!
	  IF(ND .NE. ND_ATM)THEN
	    WRITE(T_OUT,*)'Error in SUB_GUESS_DC -- ND and ND_ATM must be the same'
	    WRITE(T_OUT,*)'The files MUST be constructed by the same model atmosphere'
	    WRITE(T_OUT,*)'ND in EDDFACTOR is', ND
	    WRITE(T_OUT,*)' ND_ATM in RVTJ is', ND_ATM
	  END IF
!
	  ALLOCATE (R(ND_ATM))
	  ALLOCATE (V(ND_ATM))
	  ALLOCATE (SIGMA(ND_ATM))
	  ALLOCATE (T(ND_ATM))
	  ALLOCATE (ION_POP(ND_ATM))
	  ALLOCATE (ED(ND_ATM))
	  ALLOCATE (WORK(ND_ATM))
	  ALLOCATE (ROSS_MEAN(ND_ATM))
	  ALLOCATE (FLUX_MEAN(ND_ATM))
	  ALLOCATE (POP_ATOM(ND_ATM))
	  ALLOCATE (MASS_DENSITY(ND_ATM))
	  ALLOCATE (POPION(ND_ATM))
	  ALLOCATE (CLUMP_FAC(ND_ATM))
	  CALL RD_RVTJ_VEC_V4(R,V,SIGMA,ED,T,
	1        WORK,WORK,ROSS_MEAN,FLUX_MEAN,WORK,
	1        POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,FORMAT_DATE,ND_ATM,LU_IN)
	  CLOSE(LU_IN)
	  WRITE(T_OUT,*)'Successfully read RVTJ'
	  FLUSH(T_OUT)!
	  FIRST=.FALSE.
	END IF
! 
!
! Open Oscillator file, and read in ENEGRY names, levels, and statistical weights..
!
	FILENAME=TRIM(SPECIES)//'_F_OSCDAT'
	CALL GEN_ASCI_OPEN(LU_HEAD,'HEAD_INFO','UNKNOWN',' ','WRITE',IZERO,IOS)
	CALL RD_ENERGY(NAME,G,ENERGY,FEDGE,NLEV,N_MAX,
	1       ION_EN,ZION,EN_DATE,FILENAME,LU_IN,LU_HEAD,IOS)
	IF(IOS .NE. 0)THEN
	   WRITE(T_OUT,*)'Error occurred reading Oscillator file: try again'
	   STOP
	END IF
	CLOSE(LU_HEAD); CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read oscillator file for ',TRIM(SPECIES)
	FLUSH(T_OUT)
!
! Check EDDFACTOR file extends to high enough frequencies.
!
	IF(FEDGE(1) .GT. NU(1))THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Error --- maximum frequency in EDDFACTOR too small '
	  WRITE(T_OUT,*)'EDDFACTOR NU_MAX=',NU(1)
	  WRITE(T_OUT,*)'Required  NU_MAX>',FEDGE(1)
	  WRITE(T_OUT,*)' '
	  STOP
	END IF
!
	ML=1
	DO WHILE(FEDGE(1) .LT. NU(ML))
	  ML=ML+1
	END DO
	EDGE_ML=ML-1
!
! Compute the recombination and photoionization rate to the ground state.
! For simplicity we assume the ground state population is set by a balance
! between photoionizations and recombinations. We ignore the frequency
! dependence of the photoionization cross-section, and note that its numerical
! value does not matter.
!
	IF(.NOT. ALLOCATED(PHOT_SUM))THEN
	  ALLOCATE (PHOT_SUM(1:ND))
	  ALLOCATE (RECOM_SUM(1:ND))
	END IF
	PHOT_SUM(1:ND)=0.0_LDP
	RECOM_SUM(1:ND)=0.0_LDP
	DO ML=1,EDGE_ML-1
	  DO I=1,ND
	    PHOT_SUM(I)=PHOT_SUM(I)+0.5_LDP*(NU(ML)-NU(ML+1))*(RJ(I,ML)+RJ(I,ML+1))
	    T1=(TWOHCSQ*(NU(ML)**3)+RJ(I,ML))*EXP(-HDKT*NU(ML)/T(I))
	    T2=(TWOHCSQ*(NU(ML+1)**3)+RJ(I,ML+1))*EXP(-HDKT*NU(ML+1)/T(I))
	    RECOM_SUM(I)=RECOM_SUM(I)+0.5_LDP*(NU(ML)-NU(ML+1))*(T1+T2)
	  END DO
	END DO
!
! Add contribution where NU(EDGE_ML) is very different from FEDGE(1)
!
	IF(NU(EDGE_ML) .GT. 1.00000001_LDP*FEDGE(1))THEN
	  DO ML=1,10
	    NU1=NU(EDGE_ML)-(ML-1)*(NU(EDGE_ML)-FEDGE(1))/10
	    NU2=NU(EDGE_ML)-ML*(NU(EDGE_ML)-FEDGE(1))/10
	    T1=(NU(EDGE_ML)-NU1)/(NU(EDGE_ML)-NU(EDGE_ML+1))
	    T2=(NU(EDGE_ML)-NU2)/(NU(EDGE_ML)-NU(EDGE_ML+1))
	    DO I=1,ND
	      RJ1=(1.0_LDP-T1)*RJ(I,EDGE_ML)+T1*RJ(I,EDGE_ML+1)
	      RJ2=(1.0_LDP-T2)*RJ(I,EDGE_ML)+T2*RJ(I,EDGE_ML+1)
	      PHOT_SUM(I)=PHOT_SUM(I)+0.5_LDP*(NU1-NU2)*(RJ1+RJ2)
	      T3=(TWOHCSQ*(NU1**3)+RJ1)*EXP(-HDKT*NU1/T(I))
	      T4=(TWOHCSQ*(NU2**3)+RJ2)*EXP(-HDKT*NU2/T(I))
	      RECOM_SUM(I)=RECOM_SUM(I)+0.5_LDP*(NU1-NU2)*(T3+T4)
	    END DO
	  END DO
	END IF
!
! Compute the ground-state departure coefficient, and the excitation temperature
! of the ground state.
!
	IF(.NOT. ALLOCATED(GS_DC))THEN
	  ALLOCATE (GS_DC(1:ND))
	  ALLOCATE (T_EXC(1:ND))
	END IF
	GS_DC(1:ND)=RECOM_SUM(1:ND)/PHOT_SUM(1:ND)
	T1=1.0E-02_LDP      !Accuracy parameter
	CALL GET_EXCITE_TEMP(T,T_EXC,GS_DC,FEDGE,HDKT,T1,ND)
!
! Compute departure coefficients of all levels. We assume that they have the same
! excitation temeprature as the ground state. In many cases, this will
! be a poor approximation.
!
	IF(ALLOCATED(DC))DEALLOCATE(DC)
	ALLOCATE (DC(NLEV,ND))
	DO I=1,ND
	 DO J=1,NLEV
	   DC(J,I)=((T(I)/T_EXC(I))**1.5_LDP )*EXP(HDKT*FEDGE(J)*(1.0_LDP/T_EXC(I)-1.0_LDP/T(I)))
	 END DO
	END DO
!
! Read in ION file to get ground state population. For a lower ionization species,
! ion population does not matter, since the actual ion population gets used when
! species is read into CMFGEN (ie. level 1 of CIV will get used as the ion population
! for CIII, not DCIII. This case occurs when we pass REF_SPECIES as the same of species.
!
! Setting ION_POP to 10^{-50} should avoid (we hope) floating overflow if adding
! all species.
!
	IF(SPECIES .EQ. REF_SPECIES)THEN
	  ION_POP(1:ND)=1.0E-50_LDP
	ELSE
	  FILENAME=TRIM(REF_SPECIES)//'_IN'
          CALL GEN_ASCI_OPEN(LU_IN,FILENAME,'OLD',' ','READ',IZERO,IOS)
          I=0
          STRING=' '
          DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
            I=I+1
            READ(LU_IN,'(A)')STRING
          END DO
          IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU_IN)
          READ(LU_IN,*)T1,T2,NLEV_RD,ND_RD
	  IF(ALLOCATED(DC_RUB))DEALLOCATE(DC_RUB)
	  ALLOCATE (DC_RUB(NLEV_RD))
          IF(ND_RD .NE. ND)THEN
	    WRITE(T_OUT,*)'ND in ion file must be same as current ND'
	    STOP
	  END IF
	  IF(.NOT. ALLOCATED(GS_ION_POP))ALLOCATE (GS_ION_POP(ND))
          DO J=1,ND_RD
            READ(LU_IN,*)T1,GS_ION_POP(J)
            READ(LU_IN,*)(DC_RUB(I),I=1,NLEV_RD)
	  END DO
	  CLOSE(LU_IN)
	END IF
!
! Here we acturally compute the ION_POP. ION_POP has already been set for the
! case SPECIES=REF_SPECIES.
!
	IF(SPECIES .NE. REF_SPECIES)THEN
	  DO I=1,ND
	    T1=LOG(2.07078E-22_LDP*ED(I)*DC(1,I))
	    T2=T1+HDKT*FEDGE(1)/T(I)
	    WRITE(6,'(4ES16.4)')ED(I),DC(1,I),T1,HDKT*FEDGE(1)/T(I)
	    T1=GION_REF/(T(I)**1.5_LDP)/GION_SPEC
	    ION_POP(I)=EXP(LOG(GS_ION_POP(I)/T1)-T2)
	  END DO
!
! Limit ionization ratio to a factor of 10^10
!
	  DO I=1,ND
	    IF(ION_POP(I) .LT. 1.0E-10_LDP*GS_ION_POP(I))THEN
	      T1=1.0E-10_LDP*GS_ION_POP(I)/ION_POP(I)
	      DO J=1,NLEV
	        DC(J,I)=DC(J,I)/T1
	      END DO	
	      ION_POP(I)=T1*ION_POP(I)
	    END IF
	  END DO
	END IF
!
! Output population estimates to XzV_IN.
!
	FILENAME=TRIM(SPECIES)//'_IN'
	CALL GEN_ASCI_OPEN(LU_OUT,FILENAME,'NEW',' ','WRITE',IZERO,IOS)
	  WRITE(LU_OUT,'(/,1X,A,T40,A)')'07-Jul-1997','!Format date'
	  WRITE(LU_OUT,2120)R(ND),RLUM,NLEV,ND
	  DO I=1,ND
	    WRITE(LU_OUT,2122)R(I),ION_POP(I),ED(I),T(I),0.0,V(I),CLUMP_FAC(I),I
	    WRITE(LU_OUT,'(1X,1P,5E17.7)')(DC(J,I),J=1,NLEV)
	  END DO
	CLOSE(LU_OUT)
2120	FORMAT(/,1X,ES14.8,5X,1PE12.6,5X,0P,I4,5X,I4)
2122	FORMAT(/,1X,ES17.10,6ES17.8,2X,I4,1X)
!
	RETURN
	END
!
	SUBROUTINE DEALLOCATE_MOD_GUESS_DC
	USE SET_KIND_MODULE
	USE MOD_GUESS_DC
	IMPLICIT NONE
!
	IF(.NOT. ALLOCATED(RJ))RETURN
	DEALLOCATE (RJ)
	DEALLOCATE (NU)
	DEALLOCATE (R)
	DEALLOCATE (V)
	DEALLOCATE (SIGMA)
	DEALLOCATE (T)
	DEALLOCATE (ION_POP)
	DEALLOCATE (ED)
	DEALLOCATE (ROSS_MEAN)
	DEALLOCATE (FLUX_MEAN)
	DEALLOCATE (POP_ATOM)
	DEALLOCATE (MASS_DENSITY)
	DEALLOCATE (POPION)
	DEALLOCATE (CLUMP_FAC)
	IF(ALLOCATED(DC))DEALLOCATE(DC)
	IF(ALLOCATED(GS_ION_POP))DEALLOCATE(GS_ION_POP)
	IF(ALLOCATED(GS_DC))DEALLOCATE(GS_DC)
	IF(ALLOCATED(T_EXC))DEALLOCATE(T_EXC)
!
	RETURN
	END
