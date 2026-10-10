      SUBROUTINE RD_DUST_FILE
      USE DUST_MOD
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
! Created 12-Feb-2023
!
      CHARACTER*80 :: STRING
      INTEGER I,IOS
      INTEGER UNIT_DUST
      INTEGER, PARAMETER :: IZERO=0

      HZ_TO_MICRON	= 2.99792458D10 * 1D4 / 1D15
      
      CALL GET_LU(UNIT_DUST,'RD_DUST_FILE')
      OPEN(UNIT=UNIT_DUST, FILE='DUST_FILE', STATUS='OLD', IOSTAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(6,*)'Unable to open dust fule in RD_DUST_FILE'
          WRITE(6,*)'IOS=',IOS
          STOP
        END IF
!
      STRING = ' '
      DO WHILE (1 .EQ. 1)
         READ(UNIT_DUST,'(A)')STRING
         IF(INDEX(STRING,'# File:').NE. 0)THEN
            I=INDEX(STRING,'File:')+5
            READ(STRING(I:),*)DUST_FILE_NAME
         ENDIF
         IF(INDEX(STRING,'# Number of entries:') .NE. 0)THEN
            I=INDEX(STRING,'entries:')+8
            READ(STRING(I:),*)N_DUST_FILE
            EXIT
         ENDIF
      ENDDO
!
! Allocate storage for dust
!
      ALLOCATE (LAM_DUST_IN(N_DUST_FILE))
      ALLOCATE (KAP_ABS_DUST_IN(N_DUST_FILE)) ; KAP_ABS_DUST_IN(1:N_DUST_FILE) = 0.D0
      ALLOCATE (KAP_SCAT_DUST_IN(N_DUST_FILE)) ; KAP_SCAT_DUST_IN(1:N_DUST_FILE) = 0.D0
!
      IF (INDEX(DUST_FILE_NAME,'abs_sca').NE.0) THEN
         DUST_TYPE = 'MIX'
         PURE_ABS_DUST = .FALSE.
         DO I=1,N_DUST_FILE
            READ(UNIT_DUST,*) LAM_DUST_IN(I),KAP_ABS_DUST_IN(I),KAP_SCAT_DUST_IN(I)
         ENDDO
      ELSE IF (INDEX(DUST_FILE_NAME,'sca').NE.0) THEN
         DUST_TYPE = 'SCA'
         PURE_ABS_DUST = .FALSE.
         DO I=1,N_DUST_FILE
             READ(UNIT_DUST,*) LAM_DUST_IN(I),KAP_SCAT_DUST_IN(I)
         ENDDO
      ELSE IF (INDEX(DUST_FILE_NAME,'abs').NE.0) THEN
         DUST_TYPE = 'ABS'
         PURE_ABS_DUST = .TRUE.
         DO I=1,N_DUST_FILE
            READ(UNIT_DUST,*) LAM_DUST_IN(I),KAP_ABS_DUST_IN(I)
         ENDDO
      ELSE
         WRITE(6,*) 'Dust file name not recognized: ',DUST_FILE_NAME
         STOP
      ENDIF
      CLOSE(UNIT_DUST)
!
      WRITE(6,*) 'DUST_TYPE :', DUST_TYPE, N_DUST_FILE, DUST_FILE_NAME
      CALL FLUSH(6)

      RETURN
      END
