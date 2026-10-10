!
! Subroutine to read in the levale names and energies for use by
! WR_F_TO_S.
!
	SUBROUTINE RD_ENERGY_V2(LEVNAME,STAT_WT,ENERGY,FEDGE,N,NMAX,
	1                   RD_ONLY_N_LEVELS,IONIZATION_EN,ZION,
	1                   OSCDATE,FILNAME,LUIN,LUOUT,IOS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered: 06-Jul-2021  Changed to allow only a finite number of levels to be read.
!
	INTEGER N,NMAX
	INTEGER LUIN,LUOUT,IOS
	REAL(KIND=LDP) STAT_WT(NMAX)
	REAL(KIND=LDP) ENERGY(NMAX)
	REAL(KIND=LDP) FEDGE(NMAX)
	REAL(KIND=LDP) IONIZATION_EN,ZION
	LOGICAL RD_ONLY_N_LEVELS
	CHARACTER*(*) LEVNAME(NMAX)
	CHARACTER*(*) FILNAME,OSCDATE
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
! External functions.
!
	EXTERNAL SPEED_OF_LIGHT,ICHRLEN,RD_FREE_VAL,ERROR_LU
	REAL(KIND=LDP) SPEED_OF_LIGHT,RD_FREE_VAL
	INTEGER ICHRLEN,ERROR_LU
!
! Local variables
!
	INTEGER, PARAMETER :: IZERO=0
	REAL(KIND=LDP) SPEED_LIGHT
	INTEGER N_LOC
	INTEGER I,J,L1,BIGLEN,MAXLEN,LUER
	CHARACTER*132 STRING
	CHARACTER*40 LOCNAME(NMAX)
!
! Variables for free-format internal reads.
!
	INTEGER NEXT,STR_LEN
	CHARACTER*80 DESC
	DATA STR_LEN/80/
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
	SPEED_LIGHT=SPEED_OF_LIGHT()			!cm/s^-1
	LUER=ERROR_LU()
!
! Initialize arrays.
!
	LEVNAME(:)=' '
	FEDGE(:)=0
	ENERGY(:)=0
	STAT_WT(:)=0
!
	CALL GEN_ASCI_OPEN(LUIN,FILNAME,'OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening '//FILNAME//' in RD_ENERGY_V2'
	  IOS=2
	  RETURN
	END IF
!
! We keep reading the file until we come upon the date. From
! then on the file has to have a fixed format. All file header
! information is written to unit LUOUT. This can be the MODEL file,
! and hence will contain relevant model data.
!
50	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Date')
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading in Oscilator Date from '//
	1                  FILNAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    RETURN
	  END IF
	  IF(L1 .NE. 0)THEN
	    OSCDATE(1:11)=STRING(1:11)
	  ELSE
	    GOTO 50
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Number of energy levels')
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='N Read in RD_ENERGY_V2-'//FILNAME
	    N_LOC=RD_FREE_VAL(STRING,1,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)'Error reading in # of energy levels from '//TRIM(FILNAME)
	    IOS=3
	    RETURN
	  END IF
	  IF(RD_ONLY_N_LEVELS)THEN
	    IF(N_LOC .LT. N)THEN
	      WRITE(LUER,*)'Error reading '//TRIM(FILNAME)//' in RD_ENERGY_V2 '
	      WRITE(LUER,*)'Insufficient levels to read: N=, N(des)=',N_LOC,N
	      STOP
	    ELSE IF(N .GT. NMAX)THEN
	      WRITE(LUER,*)'Error '//TRIM(FILNAME)//' in RD_ENERGY_V2 '
	      WRITE(LUER,*)'NMAX=',NMAX,'N=',N
	      STOP
	    END IF
	  ELSE IF(N_LOC .GT. NMAX)THEN
	    WRITE(LUER,*)'Error reading '//TRIM(FILNAME)//' in RD_ENERGY_V2 - insufficient storage'
	    WRITE(LUER,*)'N read, NMAX',N_LOC,NMAX
	    STOP
	  ELSE
	    N=N_LOC
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Ionization energy')
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='Ionization energy read in RD_ENERGY_V2-'//TRIM(FILNAME)
	    IONIZATION_EN=RD_FREE_VAL(STRING,1,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1    'Error reading in Ionization Energy from '//TRIM(FILNAME)
	    IOS=4
	    RETURN
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Screened nuclear charge')
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='NW Read in RD_ENERGY_V2-'//TRIM(FILNAME)
	    ZION=RD_FREE_VAL(STRING,1,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1     'Error reading in Screened Charge from '//TRIM(FILNAME)
	    IOS=5
	    RETURN
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Number of transitions')
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)
	1     'Error reading in # of energy levels from '//TRIM(FILNAME)
	    IOS=6
	    RETURN
	  END IF
!
! Skip blank record.
!
	  READ(LUIN,'(A)')STRING
	  IF(STRING .NE. ' ')THEN
	    WRITE(LUER,*)'Error reading blank(1) from '//FILNAME
	    IOS=7
	    RETURN
	  END IF	
!
! We first read the record into a character string so that we can do
! an unformatted read on the real variables. LEVNAME name (or transition)
! must be separated by at least 2 spaces from the real data.
! Level names need not be the same length. Note that FEDGE is initially
! the excitation energy in cm^-1.
!
	  BIGLEN=0
	  MAXLEN=MIN( LEN(LEVNAME(1)),LEN(LOCNAME(1)) )
	  DO I=1,N
	    READ(LUIN,'(A)')STRING
	    L1=INDEX(STRING,'  ')-1
	    IF( L1 .LE. 0 .OR. L1 .GT. MAXLEN )THEN
	      WRITE(LUER,*)
	1            'Error reading in Level Names from '//FILNAME
	      WRITE(LUER,*)'Invalid length of name'
	      STOP
	    END IF
	    LEVNAME(I)=STRING(1:L1)
	    DO WHILE(LEVNAME(I)(1:1) .EQ. ' ')		!Strip leading blanks.
	      LEVNAME(I)(1:)=LEVNAME(I)(2:)
	    END DO
	    LOCNAME(I)=LEVNAME(I)
	    DESC='STAT_WT read in RD_ENERGY_V2-'//FILNAME
	    STAT_WT(I)=RD_FREE_VAL(STRING,L1+1,STR_LEN,NEXT,DESC)
	    DESC='FEDGE read in RD_ENERGY_V2-'//FILNAME
	    ENERGY(I)=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,DESC)
	    FEDGE(I)=(IONIZATION_EN-ENERGY(I))*SPEED_LIGHT*1.0E-15_LDP
	    BIGLEN=MAX(BIGLEN,L1)
	  END DO
!
	RETURN
	END
