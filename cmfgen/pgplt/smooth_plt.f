!
! Simple rountine to smooth data already in GRAMON.
! Simple HAN procedure is used.
!
	SUBROUTINE SMOOTH_PLT(SMOOTH_PLOT,BOX_FILTER,NHAN)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
! Altered 28-Aug-2017 - Added box filter option.
! Created 04-Aug-2017 - Based on smothing in PLT_SPEC.
!
	LOGICAL SMOOTH_PLOT(MAX_PLTS)
	LOGICAL BOX_FILTER
	INTEGER NHAN
	INTEGER, PARAMETER :: NHAN_MAX=50
!
	REAL(KIND=LDP), ALLOCATABLE :: ZV(:)
	REAL(KIND=LDP) WT(NHAN_MAX)
	REAL(KIND=LDP) T1,T2
	INTEGER I,J,K,IP,ML
	INTEGER NX
	INTEGER IST,IEND
	LOGICAL NON_MONOTONIC
!
	WRITE(6,*)'Entered smooth plot'
	WRITE(6,*)'NHAN=',NHAN
	IF(NHAN .GT. NHAN_MAX)THEN
	  WRITE(6,*)'NHAN too large -- maximum value is ',NHAN_MAX
	  WRITE(6,*)'No smoothing done'
	  RETURN
	END IF
	WRITE(6,*)BOX_FILTER
!
	NX=1
	DO IP=1,NPLTS
	  IF(SMOOTH_PLOT(IP))NX=MAX(NX,NPTS(IP))
	END DO
	IF(ALLOCATED(ZV))DEALLOCATE(ZV)
	ALLOCATE(ZV(NX))
!
! Compute weights. We either use a box-filter or HAN which
! is equivalent to binomial smothing and is a true low pass filter.
!
	IF(BOX_FILTER)THEN
	  DO J=1,NHAN
	    WT(J)=1.0D0/NHAN
	  END DO
	ELSE
	  WT=0.0D0
	  WT(1)=1.0D0
	  WT(2)=1.0D0
	  DO I=3,NHAN
	    T2=WT(1)
	    DO J=2,NHAN-1
	      T1=T2
	      T2=WT(J)
	      WT(J)=WT(J)+T1
	    END DO
	    WT(I)=1
	  END DO
          WT=WT/SUM(WT)
	END IF
!
! We check whether the X axis is monotonic. If not, we smooth each section
!
	DO IP=1,NPLTS
	  IF(SMOOTH_PLOT(IP))THEN
	    NX=NPTS(IP)
!
	    NON_MONOTONIC=.FALSE.
	    T2=CD(IP)%XVEC(2)-CD(IP)%XVEC(1)
	    DO I=1,NX-1
	      IF( (CD(IP)%XVEC(I)-CD(IP)%XVEC(I+1))*T2 .LT. 0)THEN
	        NON_MONOTONIC=.TRUE.
	        EXIT
	      END IF
	    END DO
!
	    K=NHAN
	    IF(NON_MONOTONIC)THEN
	      ZV(1:NX)=CD(IP)%DATA(1:NX)
	      IST=1
	      IEND=0
	      T2=CD(IP)%XVEC(2)-CD(IP)%XVEC(1)
	      DO WHILE(IEND .LT. NX)
	        IEND=NX
	        DO I=IST,NX-1
	          IF( (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))*T2 .LT. 0)THEN
	            IEND=I
	            EXIT
	          END IF
	        END DO
	        DO I=IST,IEND
	          T1=0.0D0
	          CD(IP)%DATA(I)=0.0D0
	          DO ML=1,NHAN
                    J=I-1-NHAN/2+ML
	            IF(J .LT. IST)J=IST             !To get correct weighting at ends
	            IF(J .GT. IEND)J=IEND
	            T1=T1+WT(ML)
	            CD(IP)%DATA(I)=CD(IP)%DATA(I)+ZV(J)*WT(ML)
	          END DO
	          CD(IP)%DATA(I)=CD(IP)%DATA(I)/T1        !T1 should be unity
	        END DO
	        IST=IEND+1
	      END DO
!
	    ELSE
	      DO I=1,NX
	        ZV(I)=CD(IP)%DATA(I)
	       END DO
	       DO I=1,NX
	         T1=0.0D0
	         CD(IP)%DATA(I)=0.0D0
	         DO ML=1,NHAN
                   J=I-1-NHAN/2+ML
	           IF(J .LT. IST)J=IST             !To get correct weighting at ends
	           IF(J .GT. IEND)J=IEND
	           T1=T1+WT(ML)
	           CD(IP)%DATA(I)=CD(IP)%DATA(I)+ZV(J)*WT(ML)
	         END DO
	         CD(IP)%DATA(I)=CD(IP)%DATA(I)/T1
	       END DO	
	    END IF
	  END IF
	END DO
!
	RETURN
	END
