C
C Set and returns those quantities that are computer dependent when
C outputing 'direct access', unformatted files.
C
C REC_SIZE_LIM is the maximum record length in bytes,
C UNIT_SIZE    is the nuber of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the numer of bytes used to represent the number.
C MAX_NUM_REC  is the maximum number of numbers that can be written to
C                 a single record.
C
	SUBROUTINE DIR_ACC_PARS(REC_SIZE_LIM,UNIT_SIZE,
	1                        WORD_SIZE,MAX_NUM_REC)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) X
C
C Created 24-Jun-1998 : Machine independent version
C
	INTEGER REC_SIZE_LIM,UNIT_SIZE,WORD_SIZE,MAX_NUM_REC
C
	REC_SIZE_LIM=4094*4
	IF(STORAGE_SIZE(X)/8 .EQ. 10)THEN
	  WORD_SIZE=8                        !To get into bytes
        ELSE
	  WORD_SIZE=STORAGE_SIZE(X)/8        !To get into bytes
	END IF
	INQUIRE (IOLENGTH=UNIT_SIZE)X
        UNIT_SIZE=WORD_SIZE/UNIT_SIZE
	MAX_NUM_REC=REC_SIZE_LIM/WORD_SIZE
C
	RETURN
	END

