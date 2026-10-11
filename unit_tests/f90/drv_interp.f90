! Driver for the interpolation-input guards added to the JFNK tree.
!
!   drv_interp <routine> <case>
!
! routine:  lin   -> LIN_INTERP
!           mon   -> MON_INTERP
!           monf  -> MON_INTERP_FAST
!           monf2 -> MON_INTERP_FAST_V2
!           mif2  -> MON_INT_FUNS_V2
!
! case (shared grid scenario):
!   valid_dec   strictly decreasing grids, valid ranges
!   valid_inc   strictly increasing grids, valid ranges
!   rnonmono    R has a local reversal
!   rrepeat     R has a repeated point
!   ureverse    U/QZR ordered opposite to R
!   unonmono    U/QZR has a repeated point
!   urange      U/QZR does not cover R (pre-existing check)
!   qzrrepeat   repeated QZR (allowed by MON_INTERP only)
!   ndsmall     ND too small
!   ninsmall    NIN/NX too small
!
! On success the driver prints RESULT_OK followed by data rows.  If the
! production routine STOPs on a guard, RESULT_OK is absent.
!
      PROGRAM DRV_INTERP
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: LDP=KIND(1.0D0)
      INTEGER, PARAMETER :: MX=64
      INTEGER, PARAMETER :: NDM=20
!
      CHARACTER(LEN=32) :: RTN, CNAME
      INTEGER ND, NIN, NX, LIN_END, NQ, NV
      REAL(KIND=LDP) R(MX), V(MX), U(MX), W(MX)
      REAL(KIND=LDP) QZR(MX), QZ(MX,4), VARRAY(MX,4)
      REAL(KIND=LDP) COEF(NDM,4), CHI(NDM)
      REAL(KIND=LDP) T1, XX, YY
      INTEGER I, J, NARGS
      LOGICAL DEC
!
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS .LT. 2)THEN
        WRITE(*,*)'usage: drv_interp <routine> <case>'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,RTN)
      CALL GET_COMMAND_ARGUMENT(2,CNAME)
!
      DEC=(CNAME .EQ. 'valid_dec' .OR. CNAME .EQ. 'rnonmono' .OR. &
           CNAME .EQ. 'rrepeat' .OR. CNAME .EQ. 'ureverse' .OR. &
           CNAME .EQ. 'unonmono' .OR. CNAME .EQ. 'urange' .OR. &
           CNAME .EQ. 'qzrrepeat')
!
      ND=NDM
      NIN=30
      NX=30
      IF(CNAME .EQ. 'ndsmall') ND=1
      IF(CNAME .EQ. 'ninsmall')THEN
        NIN=1
        NX=0
      END IF
!
! Default grids: R decreasing 1e5 -> 1e3 (CMFGEN depth order).
! U decreasing, spanning R, on a finer mesh.
!
      DO I=1,MX
        R(I)=1.0E5_LDP*(1.0E3_LDP/1.0E5_LDP)**(REAL(I-1,LDP)/(NDM-1))
        U(I)=1.2E5_LDP*(0.8E3_LDP/1.2E5_LDP)**(REAL(I-1,LDP)/(NIN-1))
        W(I)=3.0_LDP*U(I)+0.5_LDP           ! linear -> exact for LIN
      END DO
      DO I=1,NDM
        V(I)=0.0_LDP
        VARRAY(I,1)=EXP(-R(I)/3.0E4_LDP)
        VARRAY(I,2)=R(I)*R(I)
        CHI(I)=VARRAY(I,1)
      END DO
      IF(.NOT. DEC)THEN
        DO I=1,MX
          R(I)=1.0E3_LDP*(1.0E5_LDP/1.0E3_LDP)**(REAL(I-1,LDP)/(NDM-1))
          U(I)=0.8E3_LDP*(1.2E5_LDP/0.8E3_LDP)**(REAL(I-1,LDP)/(NIN-1))
          W(I)=3.0_LDP*U(I)+0.5_LDP
        END DO
        DO I=1,NDM
          VARRAY(I,1)=EXP(-R(I)/3.0E4_LDP)
          VARRAY(I,2)=R(I)*R(I)
          CHI(I)=VARRAY(I,1)
        END DO
      END IF
!
! For the MON* routines the "source" is R/VARRAY (ND) and the "target"
! is QZR (NX), which must lie inside R at both ends.
!
      DO I=1,NX
        QZR(I)=R(1)+(R(NDM)-R(1))*(0.02_LDP+0.96_LDP* &
               REAL(I-1,LDP)/MAX(NX-1,1))
      END DO
!
! Case perturbations.
!
      IF(CNAME .EQ. 'rnonmono') R(10)=R(8)
      IF(CNAME .EQ. 'rrepeat') R(6)=R(5)
      IF(CNAME .EQ. 'ureverse')THEN
        DO I=1,NIN/2
          T1=U(I)
          U(I)=U(NIN+1-I)
          U(NIN+1-I)=T1
        END DO
        DO I=1,NX/2
          T1=QZR(I)
          QZR(I)=QZR(NX+1-I)
          QZR(NX+1-I)=T1
        END DO
      END IF
      IF(CNAME .EQ. 'unonmono' .OR. CNAME .EQ. 'qzrrepeat')THEN
        U(8)=U(7)
        QZR(8)=QZR(7)
      END IF
      IF(CNAME .EQ. 'urange')THEN
        U(1)=R(1)*0.5_LDP                 ! no longer covers R(1)
        QZR(1)=R(1)*1.5_LDP               ! outside R for MON*
      END IF
!
! NQ/NV are the array leading dimensions (strides), not logical sizes.
!
      NQ=MX
      NV=MX
      LIN_END=2
      QZ=0.0_LDP
!
      IF(RTN .EQ. 'lin')THEN
        CALL LIN_INTERP(R,V,ND,U,W,NIN)
        WRITE(*,'(a)')'RESULT_OK'
        DO I=1,ND
          WRITE(*,'(i6,2(1x,es24.16e3))') I,R(I),V(I)
        END DO
      ELSE IF(RTN .EQ. 'mon')THEN
        CALL MON_INTERP(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND)
        WRITE(*,'(a)')'RESULT_OK'
        DO I=1,NX
          WRITE(*,'(i6,3(1x,es24.16e3))') I,QZR(I),QZ(I,1),QZ(I,2)
        END DO
      ELSE IF(RTN .EQ. 'monf')THEN
        CALL MON_INTERP_FAST(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND)
        WRITE(*,'(a)')'RESULT_OK'
        DO I=1,NX
          WRITE(*,'(i6,3(1x,es24.16e3))') I,QZR(I),QZ(I,1),QZ(I,2)
        END DO
      ELSE IF(RTN .EQ. 'monf2')THEN
        CALL MON_INTERP_FAST_V2(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV, &
                                R,ND,'drv_interp')
        WRITE(*,'(a)')'RESULT_OK'
        DO I=1,NX
          WRITE(*,'(i6,3(1x,es24.16e3))') I,QZR(I),QZ(I,1),QZ(I,2)
        END DO
      ELSE IF(RTN .EQ. 'mif2')THEN
        CALL MON_INT_FUNS_V2(COEF,CHI,R,ND)
        WRITE(*,'(a)')'RESULT_OK'
        DO I=1,MIN(ND-1,NDM-1)
          DO J=0,4
            XX=R(I)+(R(I+1)-R(I))*REAL(J,LDP)/4.0_LDP
            T1=XX-R(I)
            YY=COEF(I,1)*T1*T1*T1+COEF(I,2)*T1*T1+ &
               COEF(I,3)*T1+COEF(I,4)
            WRITE(*,'(2(1x,es24.16e3))') XX,YY
          END DO
        END DO
      ELSE
        WRITE(*,*)'unknown routine ',TRIM(RTN)
        STOP
      END IF
!
      END
