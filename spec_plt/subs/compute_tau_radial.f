!
! Subroutine designed to compute the radial optical depth done to a
! given velocity (VMIN). If VMIN is zero, the core optical depth is
! returned.
!
	SUBROUTINE COMPUTE_TAU_RADIAL(TAU, CHI, VEL, R, NU, dV, VMIN, ND, NF)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created: 05-Apr-2020
! Altered: 05-Oct-2026 -- Reject non-monotonic VEL.
!
	INTEGER NF		!Number of frequency points
	INTEGER ND		!Numer of depth points
	REAL(KIND=LDP) dV		!Spacing in velcity space -- set to V(DOP)/2
	REAL(KIND=LDP) VMIN		!Lower limit for integration
!
	REAL(KIND=LDP) NU(NF)		!Freuqency (10^15 Hz)
	REAL(KIND=LDP) TAU(NF)
	REAL(KIND=LDP) CHI(ND,NF)	!Opacity (scaled by 10^10) in comoving frame
	REAL(KIND=LDP) VEL(ND)		!Velocity (km/s)
	REAL(KIND=LDP) R(ND)		!Radius (in units of 10^10 cm).
!
	REAL(KIND=LDP), ALLOCATABLE :: NEW_CHI(:,:)		!New opacity (interpolated in space and frequency).
	REAL(KIND=LDP), ALLOCATABLE :: NEW_VEL(:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_R(:)
!
! Work vectors
!
	REAL(KIND=LDP), ALLOCATABLE :: TMP_CHI(:)
	REAL(KIND=LDP), ALLOCATABLE :: TMP_FREQ(:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_VEC(:)
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) T2
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER NF_INT
	INTEGER ML_END
	INTEGER ND_LIM		!Use when not integrating all the way to the core.
	INTEGER NX		!New number of depth points
	INTEGER NJ		!Used when inserying points in a depth interval.
	INTEGER I,J,ML,K
	INTEGER N_THREAD
	INTEGER N_THREAD_MAX
!
	REAL(KIND=LDP) SPEED_OF_LIGHT,C_KMS
	EXTERNAL SPEED_OF_LIGHT
!
! The grid refinement uses the velocity drop per radial interval, and the
! VMIN cutoff assumes the first crossing is the inner boundary.  Repeated
! or increasing velocities make both assumptions invalid.
!
	DO I=2,ND
	  IF(VEL(I) .GE. VEL(I-1))THEN
	    WRITE(6,*)'Error in COMPUTE_TAU_RADIAL:'
	    WRITE(6,*)'VEL must be strictly decreasing'
	    WRITE(6,*)'Depth, VEL(I-1), VEL(I)=',I,VEL(I-1),VEL(I)
	    STOP
	  END IF
	END DO
!
! Determine new grid size.
!
	NX=1
	ND_LIM=ND
	DO I=2,ND
	  NJ=1+(VEL(I-1)-VEL(I))/dV
	  NX=NX+NJ
	  IF(VEL(I) .LT. VMIN)THEN
	    ND_LIM=I
	    EXIT
	  END IF
	END DO
	WRITE(6,*)'Number of data points is ',NX
!
	ALLOCATE (NEW_CHI(NX,NF))
	ALLOCATE (NEW_R(NX))
	ALLOCATE (NEW_VEL(NX))
!
	ALLOCATE (TMP_CHI(NF))
	ALLOCATE (TMP_FREQ(NF))
	ALLOCATE (WRK_VEC(NF))
!
! Determine the new spatial grid -- maximum spacing in velcoity space is dV km/s.
!
	NEW_R(1)=R(1)
	NEW_VEL(1)=VEL(1)
	J=1
	DO I=2,ND_LIM
	  NJ=1+(VEL(I-1)-VEL(I))/dV
	  T1=(VEL(I)-VEL(I-1))/NJ
	  T2=(R(I)-R(I-1))/NJ
	  WRITE(6,*)J,NJ
	  DO K=1,NJ-1
	    J=J+1
	    NEW_VEL(J)=NEW_VEL(J-1)+T1
	    NEW_R(J)=NEW_R(J-1)+T2
	  END DO
	  J=J+1
	  NEW_VEL(J)=VEL(I)
	  NEW_R(J)=R(I)
	END DO
!
	WRITE(6,*)'Starting radial interpolation'
	CALL OMP_GET_MAX_THREADS(N_THREAD_MAX)
	CALL OMP_GET_NUM_THREADS(N_THREAD)
	J=MAX(8,N_THREAD_MAX); CALL OMP_SET_NUM_THREADS(J)
	CALL MON_INTERP_FAST(NEW_CHI,NX,NF,NEW_R,NX,CHI,ND,R,ND)
	WRITE(6,*)'Finished the radial interpolation'
!
! We now need to put verything on the same frequency grid.
!
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
!
!$OMP PARALLEL DO PRIVATE(TMP_FREQ, TMP_CHI, WRK_VEC, T1, I, ML, ML_END, NF_INT)
	DO I=1,NX
	  T1=SQRT( (1.0_LDP+NEW_VEL(I)/C_KMS)/(1.0_LDP-NEW_VEL(I)/C_KMS) )
	  TMP_FREQ(1:NF)=NU(1:NF)*T1
	  DO ML=NF,1,-1
	   IF(TMP_FREQ(NF) .LE. NU(ML))THEN
	     ML_END=ML
	     EXIT
	   END IF
	  END DO
	  NF_INT=ML_END
!	  WRITE(6,'(3I12,4ES14.4)')I,NF_INT,NF,T1,NEW_VEL(I),TMP_FREQ(NF),NU(NF_INT)
	  TMP_CHI(1:NF)=NEW_CHI(I,1:NF)
	  CALL MON_INTERP(WRK_VEC,NF_INT,IONE,NU,NF_INT,TMP_CHI,NF,TMP_FREQ,NF)
	  NEW_CHI(I,1:ML_END)=WRK_VEC(1:ML_END)
	END DO
!$OMP END PARALLEL DO
!
	TAU(1:NF)=NEW_CHI(1,1:NF)*R(1)/10.0_LDP
!$OMP PARALLEL DO PRIVATE(ML,I)
	DO ML=1,NF
	  DO I=2,NX
	    TAU(ML)=TAU(ML)+0.5_LDP*(NEW_R(I-1)-NEW_R(I))*(NEW_CHI(I-1,ML)+NEW_CHI(I,ML))
	    IF(NEW_VEL(I) .LT. VMIN)EXIT
	  END DO
	END DO
!$OMP END PARALLEL DO
!
	T1=0.0_LDP	
	DO I=2,NX
	  IF(NEW_VEL(I) .LT. VMIN)THEN
	    T1=NEW_VEL(I)
	    EXIT
	  END IF
	END DO
!
	IF(T1 .EQ. 0.0_LDP)THEN
	  WRITE(6,*)'Computed optical depth to the core'
	ELSE
	  WRITE(6,'(A,ES14.4,A)')'Computed optical depth down to a velcoity of ',T1,' kms'
	END IF
!
	DEALLOCATE (NEW_R,NEW_CHI,NEW_VEL)
	DEALLOCATE (TMP_FREQ,TMP_CHI,WRK_VEC)
!
! Reset number of threads
!
	CALL OMP_SET_NUM_THREADS(N_THREAD)
!
	RETURN
	END
