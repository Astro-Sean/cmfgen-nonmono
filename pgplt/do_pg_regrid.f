!
! Subroutine to interpolate a data vector in PGPLOT to a uniform grid in
! X space or LOG(X) space.
!
	  SUBROUTINE DO_PG_REGRID(IN1,OUT,XMIN,XMAX,REG_OPT,VAR)
	  USE MOD_CURVE_DATA
	  IMPLICIT NONE
!
! Created 16-Nov-2020
!
	  INTEGER IN1
	  INTEGER OUT
	  REAL*4 XMIN,XMAX
	  REAL*4 VAR
	  CHARACTER(LEN=*) REG_OPT
!
	  INTEGER NINS
	  REAL*4 dLAM
	  REAL*4 RES
!
	  INTEGER GET_INDX_SP
	  EXTERNAL GET_INDX_SP
!
! Local variables
!
	  REAL*4, ALLOCATABLE :: XV_OLD(:)
	  REAL*4, ALLOCATABLE :: YV_OLD(:)
	  REAL*4, ALLOCATABLE :: XV_NEW(:)
	  REAL*4, ALLOCATABLE :: YV_NEW(:)
!
	  INTEGER N1
	  INTEGER N2
!
	  REAL*4 T1
	  INTEGER I,J,L
	  INTEGER IL,IU
	  INTEGER LOW_LIM,UP_LIM
!
	  INTEGER, PARAMETER :: IONE=1
	  INTEGER, PARAMETER :: T_OUT=6
!
! Check validity of maps:
!
	  IF(NPTS(IN1) .EQ. 0)THEN
	    WRITE(T_OUT,*)'Invalid Plot ID''s in DO_VEC_OP'
	    WRITE(T_OUT,*)'Plot ID=',IN1,'  NPTS=',NPTS(IN1)
	    RETURN
	  END IF
!
	  IF(REG_OPT(1:2) .EQ. 'DX')THEN
	    dLAM=VAR
	  ELSE IF(REG_OPT(1:1) .EQ. 'R')THEN
	    RES=VAR
	  ELSE IF(REG_OPT(1:1) .EQ. 'N')THEN
	    NINS=NINT(VAR); NINS=MAX(NINS,1)
	  ELSE IF(REG_OPT(1:2) .EQ. 'UG')THEN
	  ELSE
	    WRITE(T_OUT,*)'Error in DO_PG_REGRID'
	    WRITE(T_OUT,*)'Unrecoginized optio:',TRIM(REG_OPT)
	    WRITE(T_OUT,*)'Available options are: DX, R(ES), N(NIS)'
	    RETURN
	  END IF
!
	  N1=NPTS(IN1)
	  IF(XMIN .EQ. XMAX)THEN
	    LOW_LIM=1; UP_LIM=N1
	  ELSE
	    LOW_LIM=GET_INDX_SP(XMIN,CD(IN1)%XVEC,N1)
	    UP_LIM=GET_INDX_SP(XMAX,CD(IN1)%XVEC,N1)
	    J=LOW_LIM
	    LOW_LIM=MIN(LOW_LIM,UP_LIM)
	    IF(J .NE. LOW_LIM)UP_LIM=J
	  END IF
	  WRITE(T_OUT,'(/,A,T18,2ES14.4)')' XMIN,XMAX',XMIN,XMAX
	  WRITE(T_OUT,'(A,T18,2I14)')' LOW_LIM,UP_LIM',LOW_LIM,UP_LIM
!
	  ALLOCATE(XV_OLD(NPTS(IN1)))
	  ALLOCATE(YV_OLD(NPTS(IN1)))
	  N1=UP_LIM-LOW_LIM+1
	  XV_OLD(1:N1)=CD(IN1)%XVEC(LOW_LIM:UP_LIM)
	  YV_OLD(1:N1)=CD(IN1)%DATA(LOW_LIM:UP_LIM)
	  WRITE(T_OUT,'(/,A)')' Set old vectors'
!
!   [X*(1+t1)-X]/X=T1;   R=1/T1
!
! R=Lam/dL  => dL=LAM/R =? LOG(dL)=Log(LAM)-LOG(R)
!
	  IF(REG_OPT(1:2) .EQ. 'DX')THEN
	    N2=ABS(XV_OLD(1)-XV_OLD(N1))/dLAM
	    ALLOCATE(XV_NEW(N2))
	    ALLOCATE(YV_NEW(N2))
	    J=1; IF(XV_OLD(1) .GT. XV_OLD(N1))J=-1
	    DO I=1,N2
	      XV_NEW(I)=XV_OLD(1)+(I-1)*ABS(dLAM)
	    END DO
	    WRITE(6,'(A,ES14.4)')'Uniform grid with dX=',dLAM
!
	  ELSE IF(REG_OPT(1:1) .EQ. 'R')THEN
	    T1=LOG(1.0+1.0/RES)
	    N2=ABS(LOG(XV_OLD(1)/XV_OLD(N1))/T1)
	    T1=EXP(LOG(XV_OLD(N1)/XV_OLD(1))/(N2-1))
	    ALLOCATE(XV_NEW(N2))
	    ALLOCATE(YV_NEW(N2))
	    WRITE(T_OUT,*)T1,N2
	    XV_NEW(1)=XV_OLD(1); XV_NEW(N2)=XV_OLD(N1)
	    DO I=2,N2-1
	      XV_NEW(I)=XV_OLD(1)*(T1**(I-1))
	    END DO
!
	  ELSE IF(REG_OPT(1:1) .EQ. 'N')THEN
	    N2=(N1-1)*NINS+N1
	    ALLOCATE(XV_NEW(N2))
	    ALLOCATE(YV_NEW(N2))
	    L=0
	    DO I=1,N1-1
	      L=L+1; XV_NEW(L)=XV_OLD(I)
	      T1=(XV_OLD(I+1)-XV_OLD(I))/(NINS+1)
	      DO J=1,NINS
	        L=L+1; XV_NEW(L)=XV_NEW(L-1)+T1
	      END DO
	     END DO
             XV_NEW(N2)=XV_OLD(N1)
!
	  ELSE IF(REG_OPT(1:2) .EQ. 'UG')THEN
	    N2=N1
	    ALLOCATE(XV_NEW(N2))
	    ALLOCATE(YV_NEW(N2))
	    T1=(XV_OLD(N1)-XV_OLD(1))/(N2-1)
	    XV_NEW(1)=XV_OLD(1); XV_NEW(N2)=XV_OLD(N1)
	    DO I=2,N2-1
	      XV_NEW(I)=XV_NEW(1)+(I-1)*T1
	    END DO
	    WRITE(6,'(A,ES14.4)')'Uniform grid option used with dX=',T1
!
	  ELSE
	    WRITE(T_OUT,*)'Error in DO_PG_REGRID - XV_NEW set up'
	    WRITE(T_OUT,*)'Unrecoginized optio:',TRIM(REG_OPT)
	    WRITE(T_OUT,*)'Available options are: DX, R(ES), N(NIS)'
	    RETURN
	  END IF
!
	  WRITE(T_OUT,*)'Number of data points in regridded data: ',N2
	  WRITE(6,'(A)')'Check on new X vector'
	  WRITE(T_OUT,'(A,T10,4ES14.4)')'1:4',XV_NEW(1:MIN(4,N2))
	  WRITE(T_OUT,'(A,T10,4ES14.4)')'N2-3:N2',XV_NEW(MAX(1,N2-3):N2)
!
	  NPTS(OUT)=N2
	  ERR(OUT)=.FALSE.
	  IF(OUT .GT. NPLTS)NPLTS=OUT
!
 	  CALL SP_MON_INTERP(YV_NEW,N2,IONE,XV_NEW,N2,YV_OLD,N1,XV_OLD,N1)
!
	  IF(ALLOCATED(CD(OUT)%XVEC))THEN
	    DEALLOCATE(CD(OUT)%DATA,CD(OUT)%XVEC)
	  END IF
	  ALLOCATE(CD(OUT)%DATA(N2),CD(OUT)%XVEC(N2))
!
	  CD(OUT)%XVEC(1:N2)=XV_NEW(1:N2)
	  CD(OUT)%DATA(1:N2)=YV_NEW(1:N2)
	  DEALLOCATE(XV_NEW,XV_OLD,YV_OLD,YV_NEW)
!
	  RETURN
	  END
