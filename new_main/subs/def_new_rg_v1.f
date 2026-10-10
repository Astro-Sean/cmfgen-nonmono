!
! Designed to create a new R-grid so that the variations in certain populations (as specified by LEVELS)
! are better resolved. As present, thw whole R gird is adjusted. The step size in R and scale factor
! can be passed, but are computed if pased as zero. Originally developed for SN shock model.
!
	SUBROUTINE DEF_NEW_RG_V1(RNEW,R,POPS,LEVELS,PASS_MAX_RATIO,PASS_DR_STEP,PASS_DR_FAC,
	1                           DONE_R_REV,NLEVS,NT,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 28-Aug-2024:
!
	INTEGER NLEVS
	INTEGER NT
	INTEGER ND
	REAL(KIND=LDP) RNEW(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) POPS(NT,ND)
!
! Control variable for controlling the new R grid
!
	REAL(KIND=LDP) PASS_MAX_RATIO		!Input
	REAL(KIND=LDP) PASS_DR_STEP		!Computed if 0
	REAL(KIND=LDP) PASS_DR_FAC		!Computed if 0
	INTEGER LEVELS(NLEVS)
	LOGICAL DONE_R_REV
!
! Local variables and work vectors.
!
! Used to store populations that are being used to define the new R grid.
!
	REAL(KIND=LDP), ALLOCATABLE :: OLD_POPS(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RF(:)		!Fine R grid
	REAL(KIND=LDP), ALLOCATABLE :: RTMP(:)		!Temporary stroage for the new R grid.
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)		!Used for plotting, interpolation etc
	REAL(KIND=LDP), ALLOCATABLE :: TB(:)
!
	REAL(KIND=LDP) CUR_POPS(NLEVS)
	REAL(KIND=LDP) NEXT_POPS(NLEVS)
	REAL(KIND=LDP) RATIO(NLEVS)
!
	REAL(KIND=LDP) MAX_RATIO
	REAL(KIND=LDP) DR_STEP
	REAL(KIND=LDP) DR_FAC
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) dR
	REAL(KIND=LDP) dR_STEP_RATIO
	REAL(KIND=LDP) PREV_STEP_SIZE
!
	INTEGER, PARAMETER :: NINS=30
!
	INTEGER NF		!Number of points in fine R grid
	INTEGER ND_MAX		!Uset to set temporary R grid length (which can be > ND)
	INTEGER CNT		!Interation count
	INTEGER LST_NX
	INTEGER LUOUT
	INTEGER I,J,K,L
	INTEGER IADD
!
! IOLD refer to the location in the OLD R grid.
!
	INTEGER IOLD,IOLD_SAVE
!
! Inner and outer boundary index
!
	INTEGER IB_INDX,OB_INDX
!
	CALL GET_LU(LUOUT,'In DEF_NEW_RG_F')
!
! We make a fine grid to better resolve complex stuctures.
! We simply use the data from this grid to defined the new R grid.
!
! Deine the R gid, with an extra NINS points per ineterval.
!
	NF=(NINS+1)*(ND-1)+1
	IF(ALLOCATED(RF))DEALLOCATE(RF,TA,TB)
	ALLOCATE(RF(NF),TA(NF),TB(NF))
	DO I=1,ND-1
	  K=I+NINS*(I-1)
	  T1=(R(I+1)-R(I))/(NINS+1)
	  DO J=0,NINS
	    RF(K+J)=R(I)+J*T1
	  END DO
	END DO
	RF(NF)=R(ND)
!
! Put the old pops onto the fine grid.
!
	IF(ALLOCATED(OLD_POPS))DEALLOCATE(OLD_POPS)
	ALLOCATE(OLD_POPS(NLEVS,NF))
	DO I=1,NLEVS
	  K=LEVELS(I)
	  TA(1:ND)=LOG(POPS(K,1:ND)); J=1
	  CALL MON_INTERP(TB,NF,J,RF,NF,TA,ND,R,ND)
	  OLD_POPS(I,1:NF)=EXP(TB(1:NF))
	END DO
!
! Allocate a long array so we can iterate on the new R grid until we get
! the correct number of grid points.
!
	ND_MAX=10*ND
	IF(ALLOCATED(RTMP))DEALLOCATE(RTMP)
	ALLOCATE(RTMP(ND_MAX))
!
! We leave the grid at the inner and outer boundaries intact.
!
	IB_INDX=ND-2
	OB_INDX=3
	RTMP=0.0_LDP
	RTMP(1:OB_INDX)=R(1:OB_INDX)
!
! Set control parmeters if not set. We use local values so as not tochange the passed values.
!
	DR_FAC=PASS_DR_FAC; DR_STEP=PASS_DR_STEP; MAX_RATIO=PASS_MAX_RATIO
	IF(DR_FAC .EQ. 0.0_LDP)THEN
	  DR_FAC=EXP(LOG(R(OB_INDX)/R(IB_INDX))/0.75_LDP/(IB_INDX-OB_INDX-1))
	END IF
	IF(DR_STEP .EQ. 0.0_LDP)THEN
	  DR_STEP=1.2_LDP*(R(OB_INDX)-R(IB_INDX))/(IB_INDX-OB_INDX-1)
	END IF
	CNT=1
	LST_NX=0
!
! Iterate until we get the correct number of data points
!
	WRITE(6,'(A,9X,A,3X,A,3X,A,3X,A)')' Count','dR_STEP','MAX_RATIO','Need NX',' Cur NX'
	DO WHILE(CNT .LE. 30)
!
	  IOLD=(OB_INDX-1)*(NINS+1)+1
	  K=OB_INDX+1
	  RTMP(K:)=0.0_LDP
	  CUR_POPS=OLD_POPS(:,IOLD)
	  PREV_STEP_SIZE=R(OB_INDX-1)-R(OB_INDX)
	  WRITE(6,'(I6,F16.6,F12.3,2I10)')CNT,DR_STEP,MAX_RATIO,IB_INDX-1,LST_NX
!
	  DO WHILE(RTMP(K) .EQ. 0)
	    IOLD=IOLD+1
	    RATIO(:)=OLD_POPS(:,IOLD)/CUR_POPS(:)
	    WHERE(RATIO .LT. 1)RATIO=1/RATIO
	    T1=MAXVAL(RATIO)
	    T2=RTMP(K-1)-R(IB_INDX)
	    dR=RTMP(K-1)-RF(IOLD)
	    dR_STEP_RATIO=dR/PREV_STEP_SIZE
	    WRITE(20,*)K,RTMP(K-1),dR,dR_STEP_RATIO,PREV_STEP_SIZE,DR_STEP,MAX_RATIO; FLUSH(UNIT=20)
	    IF(T2 .LT. 1.9_LDP*DR_STEP)THEN
	      IF(T2 .GT. 1.0_LDP*DR_STEP)THEN
	        RTMP(K)=RTMP(K-1)-0.5_LDP*T2
	        K=K+1
	      END IF
	      EXIT
	    ELSE IF(dR  .GT. DR_STEP .OR. T1 .GT. MAX_RATIO .OR. dR_STEP_RATIO .GT. 3.0_LDP)THEN
	      RTMP(K)=RF(IOLD)
	      PREV_STEP_SIZE=RTMP(K-1)-RTMP(K)
	      CUR_POPS=OLD_POPS(:,IOLD)
	      IOLD=IOLD+1
	      K=K+1
	    END IF
	  END DO
!
	  K=K-1;   L=0
	  DO WHILE(L .LT. 2)
	    L=L+1; IADD=0
	    TA(1:K)=RTMP(1:K)
	    J=OB_INDX+1
	    DO I=OB_INDX+1,K-1
	      WRITE(22,*)(TA(I)-TA(I+1))/(TA(I-1)-TA(I))
	      IF( (TA(I)-TA(I+1))/(TA(I-1)-TA(I)) .LT. 0.33_LDP)THEN
	        RTMP(J)=TA(I-1)-0.67_LDP*(TA(I-1)-TA(I))
	        J=J+1
	        IADD=IADD+1
	      END IF
	      RTMP(J)=TA(I)
	      J=J+1	
	    END DO
	    J=K; K=K+IADD	
	    RTMP(K-1:K)=TA(J-1:J)
	  END DO
!
	  LST_NX=K
	  IF(K .EQ. IB_INDX-1)THEN
	     EXIT
	  ELSE IF(K .EQ. IB_INDX-2)THEN
	    RATIO=0.0_LDP
	    J=1; T2=0.0_LDP
	    DO I=OB_INDX,IB_INDX-3
	      T1=RTMP(I)/RTMP(I+1)
	      IF(T1 .GT. T2)THEN
	        J=I; T2=T1
	      END IF
	    END DO
	    TA(1:ND)=RTMP(1:ND)
	    RTMP(J+1)=0.5_LDP*(TA(J)+TA(J+1))
	    RTMP(J+2:IB_INDX-1)=TA(J+1:IB_INDX-2)
	    EXIT
	  ELSE
	     T1=FLOAT(IB_INDX-OB_INDX)/FLOAT(K+1-OB_INDX)
	     DR_STEP=DR_STEP/T1**0.75_LDP
	     MAX_RATIO=EXP(LOG(MAX_RATIO)/T1)
	  END IF
	  CNT=CNT+1
	END DO
!
	RNEW(1:IB_INDX-1)=RTMP(1:IB_INDX-1)
	RNEW(IB_INDX:ND)=R(IB_INDX:ND)
!
	IF(CNT .GE. 30)THEN
	  WRITE(6,*)'Insufficient iterations to get desired number of grid pointss'
	  DONE_R_REV=.FALSE.
	ELSE
	  DONE_R_REV=.TRUE.
	  DO I=1,ND-1
	    WRITE(16,*)I,RNEW(I)
	    IF(RNEW(I) .LE. RNEW(I+1))THEN
	      WRITE(6,'(/,A)')' Error -- R is not monotonic'
	      WRITE(6,'(I6,ES20.10)')(J,RNEW(J),J=MAX(1,I-5),MIN(I+5,ND))
	      WRITE(6,*)'See file BAD_R_GRID'
	      OPEN(LUOUT,FILE='BAD_R_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	        DO J=1,ND-1
	           WRITE(16,*)J,RNEW(J)
	        END DO
	      CLOSE(LUOUT)
	      DONE_R_REV=.FALSE.
	    END IF
	  END DO
	END IF
!
	IF(DONE_R_REV)THEN
!
	  OPEN(UNIT=LUOUT,FILE='NEW_RDINR',STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(LUOUT,'(/,1X,A,/)')' 24-FEB-2004                   !Format date'
	    WRITE(LUOUT,'(ES17.10,ES14.4,5X,I4,5X,I4)')RNEW(ND),1.0D0,1,ND
	    DO I=1,ND
	      WRITE(LUOUT,'(/,ES17.10,6ES17.7,3X,I5)')RNEW(I),(1.0D0, J=1,6),I
	      WRITE(LUOUT,'(4X,ES14.7)')1.0D0
	  END DO
	  CLOSE(LUOUT)	
!
	  DO I=1,NLEVS
	    K=LEVELS(I)
	    TA(1:ND)=LOG(POPS(K,1:ND)); J=1
	    CALL MON_INTERP(TB,ND,J,RNEW,ND,TA,ND,R,ND)
	    TB(1:ND)=EXP(TB(1:ND))
	  END DO
	END IF
!
	RETURN
	END
