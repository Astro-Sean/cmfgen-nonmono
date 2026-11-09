! Driver for INS_LINE_V6 (subs/ins_line_v6.f) and
! INS_LINE_OBS_V5 (subs/ins_line_obs_v5.f).
!
! Shared alteration under test: the DO WHILE scans over TRANS_TYPE
! that locate the next blanketed line.  Upstream combined the bounds
! test and the element test with .AND., which Fortran does not
! short-circuit -- an all-SOB line deck read TRANS_TYPE(N_LINES+1)
! out of bounds.  The split DO WHILE/EXIT form ends LN_INDX at
! N_LINES+1 when no BLANK transition exists.
!
! Cases:
!   v6_sob : all-SOB deck -> early return, FREQ=NU_CONT, NFREQ=NCF
!   v6_bla : one BLANK + one SOB -> merged grid, line indices set
!            only for the BLANK line
!   v5_sob : obs grid, INCL_ALL=F, all-SOB -> continuum spacing only
!   v5_bla : obs grid, INCL_ALL=F, one BLANK -> extra freqs inserted
!   v5_all : obs grid, INCL_ALL=T, all-SOB -> all lines included
!
      PROGRAM DRV_INS_LINE
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER, PARAMETER :: NCF=15, N_LINES=2, NFREQ_MAX=20000
      REAL(KIND=LDP) FREQ(NFREQ_MAX),NU_CONT(NCF)
      REAL(KIND=LDP) NU_LINE(N_LINES),NU_STRT_LINE(N_LINES)
      REAL(KIND=LDP) VDOP(N_LINES)
      INTEGER LINES_THIS_FREQ(NFREQ_MAX)
      INTEGER LINE_ST_INDX(N_LINES),LINE_END_INDX(N_LINES)
      CHARACTER*6 TRANS_TYPE(N_LINES)
      CHARACTER(LEN=32) CASE
      INTEGER NFREQ,I,NZ
      REAL(KIND=LDP) FRAC_DOP,VINF,DVP,DVW,ESX,RXC
      CALL GETARG(1,CASE)
!
      NU_CONT=(/40.0_LDP,30.0_LDP,25.0_LDP,20.0_LDP,18.0_LDP,&
     &          16.0_LDP,14.0_LDP,12.0_LDP,10.0_LDP,8.0_LDP,&
     &          6.0_LDP,4.0_LDP,3.0_LDP,2.0_LDP,1.0_LDP/)
      NU_LINE=(/17.0_LDP,7.0_LDP/)
      NU_STRT_LINE=(/19.0_LDP,8.5_LDP/)
      IF(CASE .EQ. 'v5_all')THEN
! The profile K-loop above has no K<=N_LINES bound: under
! INCL_ALL_LINES=T it reads TRANS_TYPE(K) out of bounds as soon as the
! grid falls below the last line's blue profile edge (a latent
! upstream defect outside the patched scans).  Keeping NU_MIN inside
! the last line's red edge exercises the flag on a terminating grid.
        NU_LINE=(/39.9_LDP,39.8_LDP/)
        NU_STRT_LINE=(/39.95_LDP,39.9_LDP/)
      END IF
      VDOP=(/30.0_LDP,30.0_LDP/)
      TRANS_TYPE=(/'SOB   ','SOB   '/)
      IF(CASE .EQ. 'v6_bla' .OR. CASE .EQ. 'v5_bla')&
     &     TRANS_TYPE(1)='BLANK '
      NFREQ=0
!
      IF(CASE(1:3) .EQ. 'v6_')THEN
        FRAC_DOP=0.5_LDP; VINF=3000.0_LDP
        DVP=30.0_LDP; DVW=300.0_LDP; ESX=1000.0_LDP; RXC=2.0_LDP
        CALL INS_LINE_V6(&
     &       FREQ,LINES_THIS_FREQ,NFREQ,NFREQ_MAX,&
     &       NU_LINE,NU_STRT_LINE,VDOP,TRANS_TYPE,&
     &       LINE_ST_INDX,LINE_END_INDX,N_LINES,&
     &       NU_CONT,NCF,FRAC_DOP,VINF,DVP,&
     &       DVW,ESX,RXC,.TRUE.)
      ELSE IF(CASE .EQ. 'v5_all')THEN
        CALL INS_LINE_OBS_V5(&
     &       FREQ,NFREQ,NFREQ_MAX,&
     &       NU_LINE,NU_STRT_LINE,VDOP,TRANS_TYPE,&
     &       N_LINES,.FALSE.,&
     &       40.0_LDP,36.0_LDP,3000.0_LDP,&
     &       0.5_LDP,30.0_LDP,300.0_LDP,2000.0_LDP,&
     &       1.5_LDP,1000.0_LDP,30.0_LDP)
        WRITE(*,'(A,I6)')'NF0',NFREQ
        CALL INS_LINE_OBS_V5(&
     &       FREQ,NFREQ,NFREQ_MAX,&
     &       NU_LINE,NU_STRT_LINE,VDOP,TRANS_TYPE,&
     &       N_LINES,.TRUE.,&
     &       40.0_LDP,36.0_LDP,3000.0_LDP,&
     &       0.5_LDP,30.0_LDP,300.0_LDP,2000.0_LDP,&
     &       1.5_LDP,1000.0_LDP,30.0_LDP)
      ELSE
        CALL INS_LINE_OBS_V5(&
     &       FREQ,NFREQ,NFREQ_MAX,&
     &       NU_LINE,NU_STRT_LINE,VDOP,TRANS_TYPE,&
     &       N_LINES,.FALSE.,&
     &       40.0_LDP,1.0_LDP,3000.0_LDP,&
     &       0.5_LDP,30.0_LDP,300.0_LDP,2000.0_LDP,&
     &       1.5_LDP,1000.0_LDP,30.0_LDP)
      END IF
!
      WRITE(*,'(A,I6)')'NF ',NFREQ
      DO I=1,NFREQ
        WRITE(*,'(A,I6,1X,ES16.9)')'FR ',I,FREQ(I)
      END DO
      IF(CASE(1:3) .EQ. 'v6_')THEN
        DO I=1,N_LINES
          WRITE(*,'(A,I2,1X,I6,1X,I6)')'LIX',I,LINE_ST_INDX(I),&
     &         LINE_END_INDX(I)
        END DO
        NZ=COUNT(LINES_THIS_FREQ(1:NFREQ) .NE. 0)
        WRITE(*,'(A,I6)')'NZ ',NZ
      END IF
      END
