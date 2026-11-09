! Driver for OPEN_RW_EDDFACTOR + OUT_RV_TO_EDDFACTOR
! (new_main/subs/open_rw_eddfactor.f).
!
! Two alterations under test:
!   * filename gate on the FINISH_REC check: upstream compared
!     FILENAME to itself (always true), so the "not computed ->
!     recompute" branch fired for every file name. The fix gates on
!     FILENAME='EDDFACTOR' or 'ES_J_CONV'.
!   * IOSTAT guard on the RV_REC pointer read in
!     OUT_RV_TO_EDDFACTOR: a truncated file reports IREC=0 and is
!     treated as empty instead of aborting on the bare READ.
!
! Cases (run in a scratch dir):
!   fresh    no files -> creates EDDFACTOR + _INFO, ACCESS_F=15,
!            RV_REC record holds INITIAL_RV_REC=10
!   gate     call twice per name; EDDFACTOR re-warns+recomputes on
!            the second call (ACCESS_F=15) while NOT_EDD returns
!            silently with COMPUTE_EDDFAC=F, ACCESS_F=0
!   trunc    NOT_EDD created then truncated to one record;
!            OUT_RV_TO_EDDFACTOR reports the empty pointer and still
!            writes the RV block (ACCESS_F bumps by 3)
!
      PROGRAM DRV_EDDFACTOR
      USE SET_KIND_MODULE
      USE EDDFAC_REC_DEFS_MOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=4,ND_EXT=4,LU_EDD=30
      REAL(KIND=LDP) R(ND),V(ND),LC(ND)
      REAL(KIND=LDP) R_EXT(ND_EXT),V_EXT(ND_EXT),LC_EXT(ND_EXT)
      INTEGER ACCESS_F,IOS,IREC,I
      LOGICAL NEWMOD,COMPUTE_EDDFAC,USE_FIXED_J,TMP
      CHARACTER(LEN=32) CASE
      CHARACTER(LEN=32) FN
      INTEGER NARG
      INTEGER REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC,RECSZ
      EXTERNAL OPEN_RW_EDDFACTOR,OUT_RV_TO_EDDFACTOR,DIR_ACC_PARS
!
      NARG=COMMAND_ARGUMENT_COUNT()
      IF(NARG .NE. 1)STOP 'usage: drv_eddfactor <case>'
      CALL GET_COMMAND_ARGUMENT(1,CASE)
!
      DO I=1,ND
        R(I)=100.0_LDP-10.0_LDP*I
        V(I)=500.0_LDP*I
        LC(I)=0.1_LDP*I
      END DO
      R_EXT=R; V_EXT=V; LC_EXT=LC
      NEWMOD=.FALSE.; USE_FIXED_J=.FALSE.
! Record length must match the production formula exactly: the
! routine derives it from DIR_ACC_PARS (ifort direct-access RECL is
! in 4-byte units, so this is NOT 8*(ND_EXT+1) bytes).
      CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
      RECSZ=WORD_SIZE*(ND_EXT+1)/UNIT_SIZE
!
      SELECT CASE(TRIM(CASE))
      CASE('fresh')
        COMPUTE_EDDFAC=.FALSE.
        FN='EDDFACTOR'
        CALL OPEN_RW_EDDFACTOR(R,V,LC,ND,R_EXT,V_EXT,LC_EXT,ND_EXT, &
     &    ACCESS_F,NEWMOD,COMPUTE_EDDFAC,USE_FIXED_J,FN,LU_EDD)
        WRITE(*,'(A,1X,I5,1X,L1)')'AF',ACCESS_F,COMPUTE_EDDFAC
        INQUIRE(FILE='EDDFACTOR',EXIST=TMP)
        WRITE(*,'(A,1X,L1)')'FEX',TMP
        INQUIRE(FILE='EDDFACTOR_INFO',EXIST=TMP)
        WRITE(*,'(A,1X,L1)')'IEX',TMP
        CLOSE(LU_EDD)
        OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED', &
     &    ACCESS='DIRECT',STATUS='OLD',RECL=RECSZ,IOSTAT=IOS)
        IF(IOS .NE. 0)STOP 'reopen failed'
        READ(LU_EDD,REC=RV_REC)IREC
        WRITE(*,'(A,1X,I5)')'RVP',IREC
        CLOSE(LU_EDD)
      CASE('gate')
! First call per name creates the file (fresh-create path always
! leaves FINISH_REC = 0). Second call reads it back: the name gate
! decides whether the "not computed" warning and recompute fire.
        DO I=1,2
          COMPUTE_EDDFAC=.FALSE.
          FN='EDDFACTOR'
          CALL OPEN_RW_EDDFACTOR(R,V,LC,ND,R_EXT,V_EXT,LC_EXT,ND_EXT, &
     &      ACCESS_F,NEWMOD,COMPUTE_EDDFAC,USE_FIXED_J,FN,LU_EDD)
          CLOSE(LU_EDD)
        END DO
        WRITE(*,'(A,1X,I5,1X,L1)')'EDDAF',ACCESS_F,COMPUTE_EDDFAC
        DO I=1,2
          COMPUTE_EDDFAC=.FALSE.
          FN='NOT_EDD'
          CALL OPEN_RW_EDDFACTOR(R,V,LC,ND,R_EXT,V_EXT,LC_EXT,ND_EXT, &
     &      ACCESS_F,NEWMOD,COMPUTE_EDDFAC,USE_FIXED_J,FN,LU_EDD)
          CLOSE(LU_EDD)
        END DO
        WRITE(*,'(A,1X,I5,1X,L1)')'NOTED',ACCESS_F,COMPUTE_EDDFAC
      CASE('trunc')
        COMPUTE_EDDFAC=.FALSE.
        FN='NOT_EDD'
        CALL OPEN_RW_EDDFACTOR(R,V,LC,ND,R_EXT,V_EXT,LC_EXT,ND_EXT, &
     &    ACCESS_F,NEWMOD,COMPUTE_EDDFAC,USE_FIXED_J,FN,LU_EDD)
        CLOSE(LU_EDD)
! Truncate the data file to a single record so REC=RV_REC (2) is
! beyond EOF.
        OPEN(UNIT=LU_EDD,FILE='NOT_EDD',FORM='UNFORMATTED', &
     &    ACCESS='DIRECT',STATUS='REPLACE',RECL=RECSZ)
        WRITE(LU_EDD,REC=1)0
        CLOSE(LU_EDD)
        ACCESS_F=INITIAL_ACCESS_REC
        CALL OUT_RV_TO_EDDFACTOR(R,V,LC,ND,R_EXT,V_EXT,LC_EXT,ND_EXT, &
     &    ACCESS_F,FN,LU_EDD)
        WRITE(*,'(A,1X,I5)')'AF',ACCESS_F
        CLOSE(LU_EDD)
        OPEN(UNIT=LU_EDD,FILE='NOT_EDD',FORM='UNFORMATTED', &
     &    ACCESS='DIRECT',STATUS='OLD',RECL=RECSZ,IOSTAT=IOS)
        READ(LU_EDD,REC=RV_REC)IREC
        WRITE(*,'(A,1X,I5)')'RVP',IREC
        CLOSE(LU_EDD)
      CASE DEFAULT
        STOP 'bad case'
      END SELECT
      WRITE(*,'(A)')'DONE'
      END
