!
! Subroutine to output informatio about model atoms/ions used in a CMFGEN model.
! Species, number for full levels, and number of super levels are output.
!
! Routine also outputs the last level of each ion in LATEX format. Some
! levels will need to be revised.
!
! The raw level name is also output. This can be stripped from the file
! using, for example, visual mode in VIM.
!
	SUBROUTINE WR_ATOM_SUM(DO_LAST_ION,VERBOSE)
	USE SET_KIND_MODULE
	USE MOD_DISP
	IMPLICIT NONE
!
! Created: 26-Sep-2021
!
!Ouputs diagnostic information related to conversion of level name to LATE fomrat.
!
	LOGICAL VERBOSE
!
! Indicates if highest ioization state (for which there is only 1 level)
! is output.
!
	LOGICAL DO_LAST_ION
!
	INTEGER NF
	INTEGER I,K,ID
	INTEGER II,IS
	INTEGER NL,NUP,CNT
	INTEGER LU
	INTEGER LU_DIAG
	INTEGER NS(NUM_IONS)
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=80) TMP_STR
	CHARACTER(LEN=80) MOD_NAME
!
	CALL GET_LU(LU,'In WR_ATOM_SUM')
	OPEN(UNIT=LU,FILE='MODEL_SPEC',STATUS='OLD')
	  DO WHILE( INDEX(STRING,'_ISF') .EQ. 0)
	    READ(LU,'(A)')STRING
	  END DO
!
	  NS=0
	  DO WHILE(1 .EQ. 1)
	    IS=INDEX(STRING,'[')+1
	    II=INDEX(STRING,'_ISF')-1
	    DO ID=1,NUM_IONS
	      IF(ION_ID(ID) .EQ. STRING(IS:II))THEN
	        READ(STRING,*)I,NS(ID)
	        EXIT
	      END IF
	    END DO
	    READ(LU,'(A)',END=100)STRING
	  END DO
100	  CONTINUE
	CLOSE(LU)
	WRITE(6,'(/,A)')' Successfully read in MODEL_SPEC file to get number of super levels for each ion.'
!
	IF(VERBOSE)THEN
	  CALL GET_LU(LU_DIAG,'In WR_ATOM_SUM')
	  OPEN(UNIT=LU_DIAG,FILE='LATEX_CONV_DIAGNOSTICS',STATUS='UNKNOWN',ACTION='WRITE')
	END IF
!
	OPEN(UNIT=LU,FILE='MOD_ION_SUMMARY',STATUS='UNKNOWN',ACTION='WRITE')
!
! For each ion we need to create the latex version of the name.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    TMP_STR=SPECIES(SPECIES_LNK(ID))
	    IS=LEN_TRIM(SPECIES_ABR(SPECIES_LNK(ID)))
	    NF=ATM(ID)%NXzV_F
	    DO I=1,NION_MAX
	      II=LEN_TRIM(ION_ID(ID))
	      IF(ION_ID(ID)(IS+1:II) .EQ. GEN_ION_ID(I))THEN
	        WRITE(TMP_STR,'(I3)')I-1
	        TMP_STR=ADJUSTL(TMP_STR)
	        TMP_STR='\ion{'//TRIM(SPECIES_ABR(SPECIES_LNK(ID)))//'}{'//TRIM(TMP_STR)//'}'
	        TMP_STR(16:16)='&'
	        K=INDEX(TMP_STR,'k')
	        IF(K .NE. 0)TMP_STR(K:K)='i'
	        MOD_NAME=ATM(ID)%XzVLEVNAME_F(NF)
	        CNT=0
	        DO NL=1,NF-1
	          DO NUP=2,NF
	            IF(ATM(ID)%AXzV_F(NUP,NL) .GT. 0.0_LDP)CNT=CNT+1
	          END DO
	        END DO
	        CALL LATEX_NAME(MOD_NAME,VERBOSE,LU_DIAG)
	        WRITE(LU,'(A,T20,I4,A,I7,A,T40,A,T90,A,I10,3A)')TRIM(TMP_STR),NS(ID),'  &',NF,'  &',
	1                  TRIM(MOD_NAME),'&',CNT,'     \\','%',TRIM(ATM(ID)%XzVLEVNAME_F(NF))
	        EXIT
	      END IF
	    END DO
!
! This section adds the descripter for the final ionization stage of each
! atomic species.
!
	  ELSE IF(DO_LAST_ION .AND. ATM(MAX(ID-1,1))%XzV_PRES .AND. .NOT. ATM(ID)%XzV_PRES)THEN
	    DO I=1,NION_MAX
	      II=LEN_TRIM(ION_ID(ID-1))
	      IF(ION_ID(ID-1)(IS+1:II) .EQ. GEN_ION_ID(I))THEN
	        WRITE(TMP_STR,'(I3)')I
	        TMP_STR=ADJUSTL(TMP_STR)
	        TMP_STR='\ion{'//TRIM(SPECIES_ABR(SPECIES_LNK(ID)))//'}{'//TRIM(TMP_STR)//'}'
	        TMP_STR(16:16)='&'
	        WRITE(LU,'(A,T20,I4,A,I7,A,T90,A)')TRIM(TMP_STR),1,'  &',1,'  &','     \\'
	        EXIT
	     END IF
	    END DO
	  END IF
	END DO	
	IF(VERBOSE)THEN	
	  CLOSE(LU_DIAG)
	  WRITE(6,*)'LATEX diagnostic information written to LATEX_CONV_DIAGNOSTICS'
	END IF
	CLOSE(UNIT=LU)
!
	WRITE(6,*)' '
	WRITE(6,*)'Model ion summary has been written to MOD_ION_SUMMARY'
	WRITE(6,*)'Conversion of names to LATEX format will need to be checked'
	WRITE(6,*)' '
!
	RETURN
	END
!
! Subroutine designed to modify the level name in latex format. This
! will work better for species with newer atomic data. Due to some
! vageries in the atomic data, some latex names will need to be modified.
!
	SUBROUTINE LATEX_NAME(MOD_NAME,VERBOSE,LU_DIAG)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER LU_DIAG
	LOGICAL VERBOSE
	CHARACTER(LEN=*) MOD_NAME
!
	INTEGER I,J,K,L
!
	INTEGER, PARAMETER :: NBANG=9
	INTEGER, PARAMETER :: NSANG=8
	CHARACTER(LEN=1) BANG(NBANG)
	CHARACTER(LEN=1) SANG(NSANG)
!
	DATA BANG/'S','P','D','F','G','H','I','W','Z'/
	DATA SANG/'s','p','d','f','g','h','i','w'/
!
! Remove j designations.
!
	J=INDEX(MOD_NAME,'[')
	IF(J .NE. 0)MOD_NAME(J:)=' '
	IF(VERBOSE)WRITE(LU_DIAG,'(A)')TRIM(MOD_NAME)
	IF(VERBOSE)FLUSH(UNIT=LU_DIAG)
!
! Put S in correct format.
!
	DO L=1,NBANG
	  J=INDEX(MOD_NAME,BANG(L))
	  DO WHILE(J .NE. 0)
	     IF(MOD_NAME(J-1:J-1) .GT. '0' .AND. MOD_NAME(J-1:J-1) .LE. '9')THEN
	       MOD_NAME(J-1:)='$^'//MOD_NAME(J-1:J-1)//'$'//MOD_NAME(J:)
	     END IF
	     I=INDEX(MOD_NAME(J+3:),BANG(L))		!In case multiple
	     IF(I .NE. 0)THEN
	        J=J+I+2
	     ELSE
	        J=0
	     END IF
	  END DO
	END DO
	IF(VERBOSE)WRITE(LU_DIAG,'(A)')TRIM(MOD_NAME)
	IF(VERBOSE)FLUSH(UNIT=LU_DIAG)
!
! Now do the l states: e.g. p2 should be p$^2$.
! We do the loop several times as there can be multiple l designations.
!
	DO K=1,3
	  DO L=1,NSANG
	    J=INDEX(MOD_NAME,SANG(L))
	    IF(J .NE. 0)THEN
	      J=J+1
	      IF(MOD_NAME(J:J) .GT. '0' .AND. MOD_NAME(J:J) .LE. '9')THEN
	         DO I=1,NSANG
	           IF(MOD_NAME(J+1:J+1) .EQ. SANG(I))THEN
	              MOD_NAME(J:)=' '//MOD_NAME(J:)
	              GOTO 20
	           END IF
	         END  DO
	         MOD_NAME(J:)='$^'//MOD_NAME(J:J)//'$'//MOD_NAME(J+1:)
	      END IF
	    END IF
	  END DO
	END DO
20	CONTINUE
	IF(VERBOSE)WRITE(LU_DIAG,'(A)')TRIM(MOD_NAME)
	IF(VERBOSE)FLUSH(UNIT=LU_DIAG)
!
! Fix the parity. We remove all even designations, consistent with
! nomral atomic nomenclature.
!
	J=INDEX(MOD_NAME,'e')
	IF(J .NE. 0)MOD_NAME(J:)=MOD_NAME(J+1:)
	J=INDEX(MOD_NAME,'e')
	IF(J .NE. 0)MOD_NAME(J:)=MOD_NAME(J+1:)
	IF(VERBOSE)WRITE(LU_DIAG,'(A)')TRIM(MOD_NAME)
	IF(VERBOSE)FLUSH(UNIT=LU_DIAG)
!
! Put odd designations in correct format. The trick with the % allows
! multipe odd designations to be fixed. Recall that the core can also
! contain a parity designation.
!	
	J=INDEX(MOD_NAME,'o')
	IF(J .NE. 0)MOD_NAME(J:)='$^{\sr %}$'//MOD_NAME(J+1:)
	J=INDEX(MOD_NAME,'o')
	IF(J .NE. 0)MOD_NAME(J:)='$^{\sr o}$'//MOD_NAME(J+1:)
	J=INDEX(MOD_NAME,'%')
	IF(J .NE. 0)MOD_NAME(J:J)='o'
	IF(VERBOSE)WRITE(LU_DIAG,'(A)')TRIM(MOD_NAME)
	IF(VERBOSE)FLUSH(UNIT=LU_DIAG)
!
! Now replace _ by blank spaces.
!
	DO K=1,3
	  J=INDEX(MOD_NAME,'_')
	  IF(J .NE. 0)MOD_NAME(J:J)=' '
	END DO
	IF(VERBOSE)WRITE(LU_DIAG,'(A)')TRIM(MOD_NAME)
	IF(VERBOSE)FLUSH(UNIT=LU_DIAG)
!
	RETURN
	END
