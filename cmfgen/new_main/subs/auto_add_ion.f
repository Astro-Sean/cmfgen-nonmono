!
! Subroutine does the following:
!     (a) Determine which ionization stages have input files.
!     (b) Call a routine to create XZV_IN files for species not present.
!
	SUBROUTINE AUTO_ADD_ION()
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered: 23-Oct-2026 - Use ATM(NID)%NXzV_F in the high-stage and no-stage
!                          generation calls. ATM(ID) held the stale index of
!                          the last existing ion (or one below the species
!                          range), writing the wrong level count into the
!                          generated XzV_IN file. Only affects the
!                          AUTO_ADD_ION_STAGES=T path.
! Altered: 26-Aug-2024 - Fixed bug in call to SUB_GUESS_DC_V2 -- wrong N being passed.
! Altered: 25-Jul-2024 - Changed call to SUB_GUESS_DC_V2
!                           DCs estimates no only computed for levesl in the model.
! Altered: 12-Aug-2019 - Better error reporting.
! Altered: 20-May-2019 - Can do entire species
!
	INTEGER ISPEC
	INTEGER ID,NID
	INTEGER FST_ID,LST_ID
	LOGICAL FILE_PRES
	CHARACTER(LEN=12) TMP_STR
!
	IF(.NOT. AUTO_ADD_ION_STAGES)RETURN
!
	WRITE(6,*)'Entered AUTO_ADD_ION'
	WRITE(6,*)'To generate input files, this routine requires: RVTJ'
	WRITE(6,*)'                                                EDDFACTOR_INFO'
	WRITE(6,*)'                                                EDDFACTOR'
	WRITE(6,*)'from a previous (similar) model'
	FLUSH(UNIT=6)
!
! Determine which ions have the necessary input files.
! We first look at the lowest ionization stages.
!
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_BEG_ID(ISPEC) .EQ. 0)GOTO 1000
	  FST_ID=0
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    TMP_STR=TRIM(ION_ID(ID))//'_IN'
	    INQUIRE(FILE=TMP_STR,EXIST=FILE_PRES)
	    IF(FILE_PRES)THEN
	      FST_ID=ID
	      EXIT
	    END IF
	  END DO
!
! Note: In CMFGEN we start with the highest population when initializing the code.

! Add species below the first ionizaton limit.
! Note: In this case the second argument in SUB_GUESS_DC must be the same as the first.
!       The fourth argument will not be used.
!
!	  WRITE(6,'(L,3X,3I4)')FILE_PRES,FST_ID
	  DO NID=FST_ID-1,SPECIES_BEG_ID(ISPEC),-1
	    WRITE(6,*)'Generating input file for ',SPECIES(ISPEC),ION_ID(NID), ATM(NID)%NXzV_F
	    CALL SUB_GUESS_DC_V2(ATM(NID)%NXzV_F,ION_ID(NID),ION_ID(NID),ATM(NID)%GXZV_F(1),ATM(NID+1)%GXZV_F(1))
	  END DO
!
! Now look at the high ionization stages.
!
	  LST_ID=SPECIES_END_ID(ISPEC)
	  DO ID=SPECIES_END_ID(ISPEC)-1,SPECIES_BEG_ID(ISPEC),-1
	    TMP_STR=TRIM(ION_ID(ID))//'_IN'
	    INQUIRE(FILE=TMP_STR,EXIST=FILE_PRES)
	    IF(FILE_PRES)THEN
	      LST_ID=ID
	      EXIT
	    END IF
	  END DO
!
! Add the additonal high ionization species. This is case is distint from the low ioizaton
! case, as we need the ground stat population of the last opulation. For example, if
! adding OV we need to DION(OIV) which is the same as the ground state OV population.
!
	  DO NID=LST_ID+1,SPECIES_END_ID(ISPEC)-1
	    WRITE(6,*)'Generating input file for ',SPECIES(ISPEC),ION_ID(NID)
	    CALL SUB_GUESS_DC_V2(ATM(NID)%NXzV_F,ION_ID(NID),ION_ID(NID-1),ATM(NID)%GXZV_F(1),ATM(NID-1)%GXZV_F(1))
	  END DO
!
! When no ionization stages are present, we need to add all ionization stages.
!
	  IF(FST_ID .EQ. 0)THEN
	    WRITE(6,*)SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1; FLUSH(UNIT=6)
	    DO NID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      WRITE(6,*)'Generating input file for ',SPECIES(ISPEC),ION_ID(NID)
	      CALL SUB_GUESS_DC_V2(ATM(NID)%NXzV_F,ION_ID(NID),ION_ID(NID),ATM(NID)%GXZV_F(1),ATM(NID)%GXZV_F(1))
	    END DO
	  END IF
!
1000	  CONTINUE
	END DO
	CALL DEALLOCATE_MOD_GUESS_DC()
!
	RETURN
	END
