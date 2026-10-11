! Driver for NG_MIT_BAND_OPT_V1 with the NG_ABSTAIN switch.
!
!   drv_ng <case_file>
!
! Case file:
!   line1: NT ND NBAND
!   line2: NG_ABSTAIN (T/F)
!   next NT*ND lines: K L RDPOPS(K,L,1..4)   (J=1 newest -> J=4 oldest)
!
! POPS is initialized to RDPOPS(:,:,1) (newest estimate).  The driver
! calls NG_MIT_BAND_OPT_V1 and prints NG_DONE, MAXINC, MAXDEC and the
! final POPS.
!
      PROGRAM DRV_NG
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: LDP=KIND(1.0D0)
      INTEGER, PARAMETER :: NT=4, ND=8
!
      REAL(KIND=LDP) POPS(NT,ND), RDPOPS(NT,ND,4)
      REAL(KIND=LDP) MAXINC, MAXDEC
      LOGICAL NG_DONE, NG_ABSTAIN
      CHARACTER(LEN=256) :: CFILE
      CHARACTER(LEN=4) :: SAB
      INTEGER NT_IN, ND_IN, NBAND, I, K, L, J, LUER, IOS
!
      CALL GET_COMMAND_ARGUMENT(1,CFILE)
      OPEN(10,FILE=TRIM(CFILE),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,*)'cannot open ',TRIM(CFILE)
        STOP
      END IF
      READ(10,*)NT_IN,ND_IN,NBAND
      READ(10,*)SAB
      NG_ABSTAIN=(SAB(1:1) .EQ. 'T')
      DO I=1,NT_IN*ND_IN
        READ(10,*)K,L,(RDPOPS(K,L,J),J=1,4)
      END DO
      CLOSE(10)
!
      POPS=RDPOPS(:,:,1)
      LUER=6
      CALL NG_WRAP(POPS,RDPOPS,ND_IN,NT_IN,NBAND, &
                  NG_DONE,MAXDEC,MAXINC,NG_ABSTAIN,LUER)
!
      WRITE(*,'(a,l2)')'RESULT_OK NG_DONE=',NG_DONE
      WRITE(*,'(a,2(1x,es24.16e3))')'MAXID ', &
           MAXINC,MAXDEC
      DO L=1,ND_IN
        DO K=1,NT_IN
          WRITE(*,'(a,2i4,1x,es24.16e3)')'POP ',K,L,POPS(K,L)
        END DO
      END DO
!
      END
