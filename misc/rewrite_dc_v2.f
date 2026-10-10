!
! Program to rewrite a new Departure Coefficient file to accomodate an
! increase in the number of levels due to level splitting.
!
! Required: OLD Departure coefficient file.
! Required: OLD oscilator file
!
! Required: NEW oscilator file
!
! The number of levels for the OLD case is  determined by the number of
! levels in the DC file.
!
! The number of levels for the NEW case is  determined by the last level
! corresponding to the HIGHEST level in the OLD data set.
!
! The name of any new levels can be individually reset if there is a
! known mismatch. For example, 4z2Z renamed to 4f2Fo
! (provided not last level).
!
	PROGRAM REWRITE_DC_V2
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 01-Jun-1998 Big fix --- STRING length changed from 80 to 132
!                         so that clump factor output.
! Altered; 15-Mar-1997 !Cleaned; GEN_IN installed etc
!                        T_OUT installed ect.
! Altered: 22-OCT-1996 !Cleaned.
! Created: 16-Sep-1996
!
! Storage for oscilators etc from old oscilator file.
! A_O is needed for the call, but is not used in the main routine.
!
	INTEGER N_O
	REAL(KIND=LDP), ALLOCATABLE :: A_O(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: EDGE_O(:)
	REAL(KIND=LDP), ALLOCATABLE :: G_O(:)
	CHARACTER*50, ALLOCATABLE ::LEVNAME_O(:)
	CHARACTER*132 OLD_OSC_FILE
!
! Storage for oscilators etc from new oscilator file.
!
	INTEGER N_N
	INTEGER, ALLOCATABLE :: INDX(:)
	REAL(KIND=LDP), ALLOCATABLE :: A_N(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: EDGE_N(:)
	REAL(KIND=LDP), ALLOCATABLE :: G_N(:)
	CHARACTER*50, ALLOCATABLE ::LEVNAME_N(:)
	CHARACTER*132 NEW_OSC_FILE
!
	INTEGER, PARAMETER :: T_IN=5		!Terminal input
	INTEGER, PARAMETER :: T_OUT=6		!Terminal output
	INTEGER, PARAMETER :: LUIN=10		!File input
	INTEGER, PARAMETER :: LUSCR=11	!
	INTEGER, PARAMETER :: DCIN=15
	INTEGER, PARAMETER :: DCOUT=16
!
	REAL(KIND=LDP), ALLOCATABLE :: DC(:)
!
	INTEGER ND
	REAL(KIND=LDP) LSTAR
	REAL(KIND=LDP) RSTAR
	CHARACTER*132 DC_FILE
	CHARACTER*30 TMP_NAME
!
	REAL(KIND=LDP) GF_LEV_CUT
	REAL(KIND=LDP) EN_LEV_CUT
	REAL(KIND=LDP) Z
	REAL(KIND=LDP) T1
	CHARACTER*30 OSCDATE
	CHARACTER*132 STRING,STRING_N
	CHARACTER*132 LNK_FILE
	CHARACTER*132 DCOUT_FILE
	LOGICAL OLD_OSC_EXIST,CREATE_LNK_FILE
!
	INTEGER I,J,K,KJ
	INTEGER IZERO,IOS
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
	CHIBF=2.815E-06_LDP
	CHIFF=3.69E-29_LDP
	HDKT=4.7994145_LDP
	TWOHCSQ=0.0147452575_LDP
	OPLIN=2.6540081E+08_LDP		!pi*e*e/m/c*1.0E+10
	EMLIN=5.27296E-03_LDP		!pc*1.0E+025/4.0/pi
	DC_FILE=' '
!
	IOS=1
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(DC_FILE,'Name of file with old depart. coef.')
	  CALL GEN_ASCI_OPEN(DCIN,DC_FILE,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)WRITE(T_OUT,*)' Error opening DC file: Try again'
	END DO
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file. If it has, we save it to output. Note that
! the header to each depth (i.e. that contianing R, Ne, etc) output to the
! final file has EXACTLY the same format as the main input file.
!
	I=0
	STRING_N=' '
	DO WHILE(INDEX(STRING_N,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(DCIN,'(A)')STRING_N
	END DO
	IF( INDEX(STRING_N,'!Format date') .EQ. 0)THEN
	   REWIND(DCIN)
	   STRING_N=' '
	END IF
	READ(DCIN,*)RSTAR,LSTAR,N_O,ND
!
	ALLOCATE (A_O(N_O,N_O))
	ALLOCATE (EDGE_O(N_O))
	ALLOCATE (G_O(N_O))
	ALLOCATE (LEVNAME_O(N_O))
!
	OLD_OSC_EXIST=.FALSE.
	OLD_OSC_FILE=' '
	DO WHILE(.NOT. OLD_OSC_EXIST)
	  CALL GEN_IN(OLD_OSC_FILE,
	1       'Oscillator file assoc. with old D.C. file')
	  INQUIRE(FILE=OLD_OSC_FILE,EXIST=OLD_OSC_EXIST)
	  IF(.NOT. OLD_OSC_EXIST)
	1       WRITE(T_OUT,*)'Error opening OLD_OSC_FILE'
	END DO
	CALL GENOSC_V5(A_O,EDGE_O,G_O,
	1                   LEVNAME_O,T1,Z,
	1                   OSCDATE,N_O,I,EN_LEV_CUT,GF_LEV_CUT,
	1                   LUIN,LUSCR,OLD_OSC_FILE)
!
	DO I=1,N_O
	  LEVNAME_O(I)=ADJUSTL(LEVNAME_O(I))
	END DO
!
200	WRITE(T_OUT,'(A)',ADVANCE='NO')
	IOS=1
	NEW_OSC_FILE=' '
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(NEW_OSC_FILE,
	1    'Oscillator file to be assoc. with NEW D.C. file')
!
! Do an initial open to get the total number of available energy levels.
!
	  CALL GEN_ASCI_OPEN(LUIN,NEW_OSC_FILE,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)
	1    WRITE(T_OUT,*)'Error opening NEW OSC. file: Try again'
	END DO
!
	J=0
	DO WHILE(J .EQ. 0)
	  READ(LUIN,'(A)')STRING
	  J=INDEX(STRING,'!Number of energy levels')
	END DO
	READ(STRING,*)N_N
	CLOSE(UNIT=LUIN)
!
	ALLOCATE (DC(N_N))
	ALLOCATE (INDX(N_N))
	ALLOCATE (A_N(N_N,N_N))
	ALLOCATE (EDGE_N(N_N))
	ALLOCATE (G_N(N_N))
	ALLOCATE (LEVNAME_N(N_N))
!
! Now re-open to use standard routine to get levelnames etc.
!
	CALL GENOSC_V5(A_N,EDGE_N,G_N,
	1                   LEVNAME_N,T1,Z,
	1                   OSCDATE,N_N,I,EN_LEV_CUT,GF_LEV_CUT,
	1                   LUIN,LUSCR,NEW_OSC_FILE)
!
	DO I=1,N_N
	  LEVNAME_N(I)=ADJUSTL(LEVNAME_N(I))
	END DO
!
! Find level in NEW data set corresponding to last level in OLD DATA set.
! It is assumed that the level ordering is the same.
!
	INDX(1:N_N)=0
	DO I=1,N_N
	  DO J=1,N_O
	    IF(LEVNAME_N(I) .EQ. LEVNAME_O(J))THEN
	      INDX(I)=J
	      EXIT
	    END IF
	  END DO
	END DO
!
	DO I=1,N_N
	  K=INDEX(LEVNAME_N(I),'[')
	  IF(INDX(I) .NE. 0)THEN
	  ELSE IF(K .EQ. 0)THEN
	    DO J=1,N_O
	      IF(LEVNAME_N(I) .EQ. LEVNAME_O(J))THEN
	        INDX(I)=J
	        EXIT
	      END IF
	    END DO
	  ELSE
	    DO J=1,N_O
	      IF(LEVNAME_N(I)(1:K-1) .EQ. LEVNAME_O(J))THEN
	        INDX(I)=J
	        EXIT
	      END IF
	    END DO
	  END IF
!
	  IF(INDX(I) .EQ. 0)THEN
	    DO J=1,N_O
	      IF(EDGE_N(I) .LE. EDGE_O(N_O))THEN
	        INDX(I)=N_O
	        EXIT
	      ELSE
	        IF(EDGE_N(I) .LE. EDGE_O(J) .AND. EDGE_N(I) .GT. EDGE_O(J+1))THEN
	          INDX(I)=J
	          EXIT
	        END IF
	      END IF
	    END DO
	  END IF
100	  CONTINUE
	END DO
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file. If it has, we save it to output. Note that
! the header to each depth (i.e. that contianing R, Ne, etc) output to the
! final file has EXACTLY the same format as the main input file.
!
	IF(STRING_N .NE. ' ')THEN
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')TRIM(STRING_N)
	END IF
	WRITE(DCOUT,'(A)')' '
	WRITE(DCOUT,'(F9.4,3X,1P,E12.4,5X,0P,I4,5X,I4)')RSTAR,LSTAR,N_N,ND
!
! Can now read in DC file an rewrite out.
!
	DO K=1,ND
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ')
	    READ(DCIN,'(A)')STRING
	  END DO
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')STRING(1:LEN(STRING))
	  READ(DCIN,*)(DC(J),J=1,N_O)
	  WRITE(DCOUT,'(1X,1P,5E15.5)')(DC(INDX(I)),I=1,N_N)
	END DO
!
! Create a one-to-one link file which may be simpler to use.
!
	CREATE_LNK_FILE=.TRUE.
	CALL GEN_IN(CREATE_LNK_FILE,'Create a link file')
	IF(CREATE_LNK_FILE)THEN
	  LNK_FILE='LINK'
	  CALL GEN_IN(LNK_FILE,'Name for link file')
	  CALL GEN_ASCI_OPEN(DCOUT,LNK_FILE,
	1              'UNKNOWN',' ','WRITE',IZERO,IOS)
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A,A)')'New oscillator file was ',
	1         TRIM(NEW_OSC_FILE)
	  WRITE(DCOUT,'(A,A)')'Old oscillator file was ',
	1         TRIM(OLD_OSC_FILE)
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(2X,A,6X,A,T25,A,T55,A)')'New','Old',
	1              'New name','Old name'
	  WRITE(DCOUT,'(A)')'!'		!Signifies beginning of data.
	  DO I=1,N_N
	    WRITE(DCOUT,'(1X,I4,5X,I4,T25,A,T55,A)')I,INDX(I),
	1       TRIM(LEVNAME_N(I)),TRIM(LEVNAME_O(INDX(I)))
	  END DO
	END IF
!
	STOP
	END
