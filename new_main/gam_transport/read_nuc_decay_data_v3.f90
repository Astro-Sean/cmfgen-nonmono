!	This subroutine is going to read the bottom half of NUC_DECAY_DATA
!	which is the part containing the nuclear decay gamma energies
!	It reads them in and stores them into vectors, so each line is linked
!	to its species, atomic mass, and probability
!
!       CREATED JULY 2013
!	Edited 8/28/13 to include storing into structure
!	Edited July 15, 2014		Changed structure type so that it looks more
!					similar to John Hillier's format for his structure
!------------------------------------------------------------------------------------------
!
	SUBROUTINE GAM_NUC_DECAY_DATA_SUB_V3(VERBOSE_GAMMA)
	USE SET_KIND_MODULE
	USE GAMMA_NUC_DECAY_V2
	IMPLICIT NONE
!
! Altered 19-Nov-2021 - Added VERBOSE_GAMMA to call
!
	CHARACTER(LEN=140) :: STRING
	CHARACTER(LEN=30) :: FILENAME
	CHARACTER(LEN=4) :: SPECIES
	CHARACTER(LEN=12) :: ATOM_MASS, E_KIN
	CHARACTER(LEN=6) :: GAMMA_E, PROB
	INTEGER  :: LINE_INDEX, INDX
	INTEGER :: COUNTER, IOS
	INTEGER :: I,J,K,L,R,S
	INTEGER :: LU
	LOGICAL VERBOSE_GAMMA
!
! Opening and skipping to the probabilities and energies part since
! rd_nuc_decay_data.f already reads in the top part of the file
!
	WRITE(6,*)'Running GAM_NUC_DECAY_DATA_SUB_V3'
	LU=340
	OPEN(UNIT=LU, FILE='NUC_DECAY_DATA',STATUS="OLD",ACTION="READ",IOSTAT=IOS)
	IF (IOS .NE. 0) STOP "*** CANNOT OPEN NUC_DECAY_DATA FILE***"
!
	STRING=' '
	DO WHILE (INDEX(STRING, 'fully at decay location') .EQ. 0)
	  READ(LU, '(A)') STRING
	END DO
	DO WHILE (INDEX(STRING, 'NICK') .EQ. 0)
	  READ(LU, '(A)') STRING
	END DO
!
	BACKSPACE(LU)
!
! NOW LOOP OVER ALL THE SPECIES FIRST IN THE EXTERIOR DO LOOP
! STORING THE SPECIES AND ATOMIC MASS INTO VECTORS
! IN THE SECOND (INTERIOR) DO LOOP WE ASSIGN THE EMISSION ENERGIES
! AND PROBABILITIES FOR EACH SPECIES STORING THEM INTO VECTORS
!
	INDX=0
	LINE_INDEX=0
	DO I=1,SPECIES_MAX
	  READ(LU, '(A)', END=200) STRING
	  K=INDEX(STRING,' ')
	  INDX=INDX+1
	  IF (INDX .EQ. SPECIES_MAX) STOP "**(FOR READ_ATOMIC_DATA) NEED TO CHANGE SPECIES_MAX TO INCLUDE MORE SPECIES **"
	  SPECIES=STRING(1:K-1)
	  GAM_ISO(INDX)%SPECIES=SPECIES
	  STRING=STRING(K:)
	  STRING=ADJUSTL(STRING)
	  STRING=TRIM(STRING)
	  ATOM_MASS=STRING
	  READ(ATOM_MASS,*)GAM_ISO(INDX)%ATOMIC_MASS
!
! Skipping line "E (MeV) Probability" since no information to extract
!
	  READ(LU,'(A)') STRING
!
! Reading in the different photon energies and probabilities and storing them in vectors and structures
!
	  COUNTER=0	! I will eventually set NUM_GAM=COUNTER so I need to reset it to zero for each species
	  DO J=1,LINE_MAX
	    READ(LU,'(A)') STRING
	    IF ((INDEX(STRING, '!Mean') .NE. 0) .OR. (INDEX(STRING, '! tot') .NE. 0)) THEN
	      K=INDEX(STRING, ' ')
	      E_KIN=STRING(1:K-1)
	      READ(E_KIN,*) E_KIN_VEC(INDX)
	      GAM_ISO(INDX)%KIN_ENERGY=E_KIN_VEC(INDX)
!	WRITE(6,*)"E_KIN:",GAM_ISO(INDX)%KIN_ENERGY
	      GAM_ISO(INDX)%NUM_GAMMA=COUNTER
!
! This next read statement was added because an extra blank line was added to the format of the NUC_DECAY_DATA file
!
	      READ(LU,'(A)',END=200) STRING
	      EXIT
	    ELSE
	      COUNTER=COUNTER+1
	      LINE_INDEX=LINE_INDEX+1
	      IF (LINE_INDEX .EQ. LINE_MAX) STOP "**(FOR READ_ATOMIC_DATA) INCREASE LINE_MAX TO INCLUDE MORE GAMMA LINES**"
	      K=INDEX(STRING,' ')
	      GAMMA_E=STRING(1:K-1)
	      STRING=STRING(K:)
	      STRING=ADJUSTL(STRING)
	      STRING=TRIM(STRING)
	      PROB=STRING
	      READ(PROB,*) PROB_VEC(LINE_INDEX)
	      READ(GAMMA_E,*) GAMMA_E_VEC(LINE_INDEX)
	      READ(ATOM_MASS,*) ATOM_MASS_VEC(LINE_INDEX)
	      SPECIES_VEC(LINE_INDEX)=SPECIES
	      GAM_ISO(INDX)%E_GAMMA(COUNTER)=GAMMA_E_VEC(LINE_INDEX)
	      GAM_ISO(INDX)%PROB(COUNTER)=PROB_VEC(LINE_INDEX)
	    END IF
	  END DO
	END DO
        200 CONTINUE
        N_GAMMA=LINE_INDEX
        N_SPECIES=INDX
	WRITE(6,'(A,1X,I3)')' # of Gamma-Ray decay lines read:',N_GAMMA
	WRITE(6,'(A,1X,I3)')' # of isotopes read:',N_SPECIES
	CLOSE(LU)
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/GAMMA_RAY_LINE_INFO',STATUS='UNKNOWN',ACTION='WRITE')
	  DO I=1,N_SPECIES
	    WRITE(7,'(A10,2X,F11.6)')GAM_ISO(I)%SPECIES,GAM_ISO(I)%ATOMIC_MASS
	    WRITE(7,'(A10,4X,A10)')'ENERGY','PROB'
	    DO J=1,GAM_ISO(I)%NUM_GAMMA
	      WRITE(7,'(F10.4,4X,F10.4)')GAM_ISO(I)%E_GAMMA(J),GAM_ISO(I)%PROB(J)
	    END DO
	    WRITE(7,*)' '
	  END DO
	  CLOSE(UNIT=7)
	END IF
!
	RETURN
        END SUBROUTINE GAM_NUC_DECAY_DATA_SUB_V3
