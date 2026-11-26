	SUBROUTINE RD_MOD(OBS_FREQ,OBS_FLUX,NCF_MAX,NCF,FILENAME,IOS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 09-Aug-2022: Changed to detect possible problem with type of input file.
C Altered 04-Jan-1999: Can now read in an arbitrary number of frequencies
C                       on a given line.
C Altered 02-May-1995: IOS installed in call.
C Altered 02-Aug-1995: Bug fix --- OBS_FREQ needed to be zeroed for correct
C                         input operations.
C
	INTEGER NCF,NCF_MAX
	REAL(KIND=LDP) OBS_FREQ(NCF_MAX)
	REAL(KIND=LDP) OBS_FLUX(NCF_MAX)
	CHARACTER*(*) FILENAME
C
	INTEGER IOS
	INTEGER, PARAMETER :: T_OUT=6
C
	CHARACTER*132 STRING
	CHARACTER(LEN=10) ANS
	INTEGER I,CNT
C
	OPEN(UNIT=10,FILE=FILENAME,ACTION='READ',STATUS='OLD',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to open file in RD_MOD'
	  RETURN
	END IF
C
C The initialization required for when I find out how many numbers I have
C read in.
C
	OBS_FREQ(1:NCF_MAX)=0.0_LDP
	OBS_FLUX(1:NCF_MAX)=0.0_LDP
C
C Continue reading observed frequencies until we hit FLUX header string.
C
	STRING=' '
	CNT=0
	DO WHILE (INDEX(STRING,'Frequencies') .EQ. 0)
	  READ(10,'(A)')STRING
	  CNT=CNT+1
	  IF(CNT .EQ. 30)THEN
	    WRITE(6,'(/,/,A)')'30 lines have been read -- does this file have the correct OBSFLUX format?'
	    WRITE(6,'(A)',ADVANCE='NO')'Do you wish to continue (return continues): '
	    READ(5,'(A)')ANS
	    IF(ANS(1:1) .NE. ' ')THEN
	      CLOSE(UNIT=30)
	      RETURN
	    END IF
	  END IF
	END DO
	READ(10,*,ERR=50)( OBS_FREQ(I),I=1,NCF_MAX)
50	CONTINUE
C
C Check how many frequencies have been read in.
C
	NCF=NCF_MAX
	DO WHILE (OBS_FREQ(NCF) .EQ. 0)
	  NCF=NCF-1
	END DO
	IF(NCF .EQ. NCF_MAX)THEN
	   WRITE(T_OUT,*)'Error --- NCF_MAX limit in RD_MOD is too small'
	   STOP
	END IF
C
C Now read in the observed fluxes.
C
	BACKSPACE(10)
	STRING=' '
	DO WHILE (INDEX(STRING,'intensity') .EQ. 0)
	  READ(10,'(A)')STRING
	END DO
	READ(10,*)( OBS_FLUX(I),I=1,NCF)
	CLOSE(UNIT=10)
C
	WRITE(T_OUT,*)NCF,' Data values read in'
C
	RETURN
	END
