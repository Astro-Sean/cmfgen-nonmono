! Driver for the Rybicki-Hummer non-local Sobolev coupling.
!
!   drv_sob_rh <case_file>
!
! Case file layout:
!   line1: ND NSIM FL
!   line2: NEG_OPACITY(1..ND) as T/F
!   next ND lines: R V SIGMA
!   next ND lines: CHIL ETAL BBCOR ZNET   (per SIM block, NSIM blocks)
!
! The driver calls SOB_RH_GEOM, dumps the cached pair tables (via
! USE SOB_RH_MOD), then calls SOB_RH_COUPLE and prints DRC plus the
! diagnostic counters.  Geometry+module access makes the test a real
! check of the production pair selection, not a reimplementation.
!
      PROGRAM DRV_SOB_RH
      USE SET_KIND_MODULE
      USE SOB_RH_MOD
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: NDX=64, MS=4
      REAL(KIND=LDP) R(NDX), V(NDX), SIGMA(NDX)
      REAL(KIND=LDP), ALLOCATABLE :: CHIL(:,:), ETAL(:,:), BBC(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: ZNET(:,:), DRC(:,:)
      LOGICAL NEG(NDX)
      CHARACTER(LEN=256) :: CFILE
      INTEGER ND, NSIM, I, J, K, L, LU, IOS
      REAL(KIND=LDP) FL
!
      CALL GET_COMMAND_ARGUMENT(1,CFILE)
      OPEN(10,FILE=TRIM(CFILE),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,*)'cannot open ',TRIM(CFILE)
        STOP
      END IF
      READ(10,*)ND,NSIM,FL
      ALLOCATE(CHIL(ND,NSIM),ETAL(ND,NSIM),BBC(ND,NSIM), &
               ZNET(ND,NSIM),DRC(ND,NSIM))
      READ(10,*)(NEG(I),I=1,ND)
      DO I=1,ND
        READ(10,*)R(I),V(I),SIGMA(I)
      END DO
      DO L=1,NSIM
        DO I=1,ND
          READ(10,*)CHIL(I,L),ETAL(I,L),BBC(I,L),ZNET(I,L)
        END DO
      END DO
      CLOSE(10)
!
      LU=6
      CALL SOB_RH_GEOM(R,V,SIGMA,ND,LU)
!
      WRITE(*,'(a,i8)')'TOTPAIRS ',RH_TOTPAIRS
      DO I=1,ND
        DO K=1,RH_NP(I)
          J=RH_ID(I,K)
          WRITE(*,'(a,3i6,3(1x,es24.16e3))')'PAIR ',I,K,J, &
        RH_T1(I,K),RH_T1P(I,K),RH_PW(I,K)
        END DO
      END DO
!
      CALL SOB_RH_COUPLE(DRC,ZNET,CHIL,ETAL,BBC,NEG, &
        V,SIGMA,R,FL,ND,NSIM,LU)
!
      WRITE(*,'(a,2i8)')'NCALL_ND ',RH_NCALL,RH_ND
      WRITE(*,'(a,2(1x,es24.16e3))')'MAXDJZ ',RH_MAXDJ,RH_MAXDZ
      DO L=1,NSIM
        DO I=1,ND
          WRITE(*,'(a,2i6,1x,es24.16e3)')'DRC ',L,I,DRC(I,L)
        END DO
      END DO
!
      END
