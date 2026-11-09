! Driver for GET_JH_AT_PREV_TIME_STEP (plane/get_jh_at_prev_time_step.f).
!
! Alterations under test:
!   - IOSTAT on the header/R-V/J-H direct-access reads so a truncated
!     JH_AT_OLD_TIME fails loudly instead of interpolating garbage.
!   - OLD_V must be positive and strictly decreasing before LOG(V) is
!     used as the MON_INTERP coordinate; the same guard runs on the
!     new V grid when INIT is true.
!
! File layout manufactured by the driver (matches the production
! writer): REC 3 = (ST_IREC,NCF_OLD,ND_OLD), REC ST_IREC = OLD_R/OLD_V,
! REC ST_IREC+1 = grey (ND_OLD J + ND_OLD H + 2 BC scalars), then
! NCF_OLD frequency records (ND_OLD J + ND_OLD-1 H + 2 BC + NU).
! Frequency records are always read in NSM=100 blocks, so NCF_OLD is
! kept at 100 to avoid the record-rewind path.
!
! Usage: drv CASE   (runs inside a scratch cwd holding the two files)
!
      PROGRAM DRV_JH_PREV
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=4, ND_OLD=4, NCF=100, ST_IREC=4
      INTEGER, PARAMETER :: LUW=20, LUI=21
      REAL(KIND=LDP) R(ND),V(ND),NEW_J(ND),NEW_H(ND-1)
      REAL(KIND=LDP) OLD_R(ND_OLD),OLD_V(ND_OLD)
      REAL(KIND=LDP) JWRK(ND_OLD),HWRK(ND_OLD-1)
      REAL(KIND=LDP) NU,H_INBC,H_OUTBC,DT,NUW
      CHARACTER(LEN=32) CASE
      CHARACTER*30 FILE_DATE
      INTEGER I,ML,RECU,IOS
      INTEGER REC_SIZE_LIM,UNIT_SIZE,WORD_SIZE,MAX_NUM_REC
      LOGICAL WRITE_HDR,WRITE_RV,WRITE_DATA
      CALL GETARG(1,CASE)
!
      OLD_R=(/10.0_LDP,8.0_LDP,6.0_LDP,4.0_LDP/)
      OLD_V=(/40.0_LDP,30.0_LDP,20.0_LDP,10.0_LDP/)
      R=(/11.0_LDP,9.0_LDP,7.0_LDP,5.0_LDP/)
      V=(/35.0_LDP,25.0_LDP,15.0_LDP,10.0_LDP/)
      NU=50.0_LDP
      WRITE_HDR=.TRUE.; WRITE_RV=.TRUE.; WRITE_DATA=.TRUE.
      IF(CASE .EQ. 'oldv_neg')OLD_V(2)=-1.0_LDP
      IF(CASE .EQ. 'oldv_flat')OLD_V=(/40.0_LDP,30.0_LDP,30.0_LDP,10.0_LDP/)
      IF(CASE .EQ. 'v_neg')V(2)=0.0_LDP
      IF(CASE .EQ. 'v_flat')V=(/35.0_LDP,25.0_LDP,25.0_LDP,10.0_LDP/)
      IF(CASE .EQ. 'inner_mm')V(4)=9.0_LDP
      IF(CASE .EQ. 'outer_big')V(1)=45.0_LDP
      IF(CASE .EQ. 'dt_zero')R(4)=3.0_LDP
      IF(CASE .EQ. 'nu_high')NU=200.0_LDP
      IF(CASE .EQ. 'nu_low')NU=0.5_LDP
      IF(CASE .EQ. 'bad_hdr')WRITE_HDR=.FALSE.
      IF(CASE .EQ. 'bad_rv')WRITE_DATA=.FALSE.
!
      CALL DIR_ACC_PARS(REC_SIZE_LIM,UNIT_SIZE,WORD_SIZE,MAX_NUM_REC)
      RECU=(2*ND_OLD+2)*WORD_SIZE/UNIT_SIZE
      OPEN(UNIT=LUW,FILE='JH_AT_OLD_TIME',STATUS='REPLACE',&
     &     ACCESS='DIRECT',FORM='UNFORMATTED',RECL=RECU)
      WRITE(LUW,REC=1)1,2
      WRITE(LUW,REC=2)3,4
      IF(WRITE_HDR)WRITE(LUW,REC=3)ST_IREC,NCF,ND_OLD
      IF(WRITE_DATA)THEN
        IF(WRITE_RV)WRITE(LUW,REC=ST_IREC)OLD_R,OLD_V
        WRITE(LUW,REC=ST_IREC+1)(1.5_LDP,I=1,ND_OLD),&
     &       (0.2_LDP,I=1,ND_OLD),0.3_LDP,0.4_LDP
        DO ML=1,NCF
          JWRK=REAL(ML,KIND=LDP)
          HWRK=0.1_LDP*ML
          NUW=REAL(101-ML,KIND=LDP)
          WRITE(LUW,REC=ST_IREC+1+ML)(JWRK(I),I=1,ND_OLD),&
     &         (HWRK(I),I=1,ND_OLD-1),2.0_LDP*ML,3.0_LDP*ML,NUW
        END DO
      END IF
      CLOSE(LUW)
      FILE_DATE='01-Jan-2000'
      CALL WRITE_DIRECT_INFO_V3(2*ND_OLD+2,RECU,FILE_DATE,&
     &                          'JH_AT_OLD_TIME',LUI)
!
      NEW_J=-1.0_LDP; NEW_H=-1.0_LDP
      H_INBC=-1.0_LDP; H_OUTBC=-1.0_LDP; DT=-1.0_LDP
      IF(CASE .EQ. 'ogrey')THEN
        CALL GET_JH_AT_PREV_TIME_STEP(NEW_J,NEW_H,&
     &       H_INBC,H_OUTBC,DT,NU,R,V,ND,.TRUE.,'OGREY')
      ELSE
        CALL GET_JH_AT_PREV_TIME_STEP(NEW_J,NEW_H,&
     &       H_INBC,H_OUTBC,DT,NU,R,V,ND,.TRUE.,'NORMAL')
      END IF
      DO I=1,ND
        WRITE(*,'(A,I2,1X,ES16.9)')'NJ ',I,NEW_J(I)
      END DO
      DO I=1,ND-1
        WRITE(*,'(A,I2,1X,ES16.9)')'NH ',I,NEW_H(I)
      END DO
      WRITE(*,'(A,ES16.9)')'HIBC ',H_INBC
      WRITE(*,'(A,ES16.9)')'HOBC ',H_OUTBC
      WRITE(*,'(A,ES16.9)')'DT   ',DT
      END
