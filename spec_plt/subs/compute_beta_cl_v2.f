!
	SUBROUTINE COMPUTE_BETA_CL_V2(BETA_VEC,NBETA,TOP_BOT_SYM,
	1               BETA_LAW,PARAMS,NPAR)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 04-Sep-2021 : Changed to V2; and BETA_LAW, PARAMS, NPAR
! Created 10-Aug-2021 (Extracted from write_line_data_jh.f)
!
! Input values
!
	INTEGER NBETA
	LOGICAL TOP_BOT_SYM
	LOGICAL UNIFORM_BETA
!
	INTEGER NPAR
	REAL(KIND=LDP) PARAMS(NPAR)
	CHARACTER(LEN=*) BETA_LAW
!
! Returned.
!
	REAL(KIND=LDP)  BETA_VEC(NBETA)
!
! Local variables
!
	REAL(KIND=LDP) PI,T1,T2
	REAL(KIND=LDP) dBETA_MIN, dBETA_MAX
	INTEGER I,K
!
	IF(MOD(NBETA,2) .NE. 1)THEN
	  WRITE(6,*)'Error in COMPUTE_BETA_CL_V2'
	  WRITE(6,*)'NBETA must be ODD'
	  STOP
	END IF
!
! Two options are present to specify the azimuth angle Beta.
! (i) BETA uniformly spaced from 0 to PI/2 (0 to PI if not TOP_BOT_SYM).
! (i) COS_BETA uniformly spaced from 1 to 0 (-1 to 1 if not TOP_BOT_SYM).
!
! Interpolation section in LONG_CHAR may need to be changed, if BETA selction
! is changed.
!
! NB. Make sure this is a double precision calculation.
!
	PI=4.0*ATAN(1.0D0)
	IF(BETA_LAW .EQ. 'UNIFORM_BETA')THEN
          IF(TOP_BOT_SYM)THEN
            T1=0.5D0*PI/(NBETA-1.0D0)
            DO I=1,NBETA
              BETA_VEC(I)=(I-1)*T1
            END DO
          ELSE
            T1=PI/(NBETA-1.0D0)
            DO I=1,NBETA
              BETA_VEC(I)=(I-1)*T1
            END DO
          END IF
        ELSE IF(BETA_LAW .EQ. 'UNIFORM_COSB')THEN
          IF(TOP_BOT_SYM)THEN
            T1=1.0D0/(NBETA-1.0D0)
            DO I=1,NBETA
              T2=1.0D0-(I-1)*T1
              BETA_VEC(I)=ACOS(T2)
            END DO
          ELSE
            T1=2.0D0/(NBETA-1.0D0)
            DO I=1,NBETA
              T2=1.0D0-(I-1)*T1
              BETA_VEC(I)=ACOS( T2 )
            END DO
          END IF
        ELSE IF(BETA_LAW .EQ. 'POWER')THEN
          IF(TOP_BOT_SYM)THEN
            T1=1.0D0/(NBETA-1.0D0)
            DO I=1,NBETA
              T2=((I-1)*T1)**PARAMS(1)
              BETA_VEC(I)=0.5*PI*T2
            END DO
          ELSE
            T1=2.0D0/(NBETA-1.0D0)
            DO I=1,NBETA/2+1
              T2=((I-1)*T1)**PARAMS(1)
              BETA_VEC(I)=0.5*PI*T2
	      BETA_VEC(NBETA-I+1)=PI-BETA_VEC(I)
            END DO
          END IF
        ELSE IF(BETA_LAW .EQ. 'BMIN')THEN
          IF(TOP_BOT_SYM)THEN
            T1=1.0D0/(NBETA-1.0D0)
	    K=NINT(PARAMS(1))
            DO I=1,K
              BETA_VEC(I)=(I-1)*PARAMS(2)
	    END DO
	    T2=(0.5*PI-BETA_VEC(K))/(NBETA-K)
            DO I=K+1,NBETA
              T1=(I-K)*T2
              BETA_VEC(I)=BETA_VEC(K)+T1
            END DO
          ELSE
	  END IF
        ELSE 
	  WRITE(6,*)'Unrecognized law for choosing Beta in COMPUTE_BETA_CL_V2'
	  WRITE(6,*)'BETA_LAW=',BETA_LAW
          STOP
	END IF
!
	dBETA_MIN=100.0D0
	dBETA_MAX=-100.0D0
	DO I=1,NBETA-1
	  WRITE(45,'(I5,3ES14.4)')I,BETA_VEC(I),BETA_VEC(I+1)-BETA_VEC(I)
	  dBETA_MIN=MIN(dBETA_MIN,BETA_VEC(I+1)-BETA_VEC(I))
	  dBETA_MAX=MAX(dBETA_MAX,BETA_VEC(I+1)-BETA_VEC(I))
	END DO
	WRITE(45,'(I5,3ES14.4)')NBETA,BETA_VEC(NBETA)
	CLOSE(UNIT=45)
	WRITE(6,*)'dBETA (min and max) are', dBETA_MIN, dBETA_MAX
	WRITE(6,*)'BETA(1), BETA(NMAX) are', BETA_VEC(1), BETA_VEC(NBETA)
!
	RETURN
	END
