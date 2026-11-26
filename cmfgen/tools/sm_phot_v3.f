!============================================================================
!
	SUBROUTINE SM_PHOT_V3(NU,CROSS,NCROSS,NSM_MAX,
	1                        SIG_GAU_KMS,FRAC_SIG_GAU,
	1                        CUT_ACCURACY,ABOVE)
	USE SET_KIND_MODULE
!
! Altered: 17-Sep-2018 : Bug fix -- caused crash with certain opacity cross-sections.
! Altered:  8-Oct-2008 : Bug fixed
!                        Ensure cross section falls off at least as 1/nu
!                         in extrapolation region.
!
!============================================================================
!
! Routine to smooth an arbitrary photionization cross-section. Cross-section
! is smoothed with a Gaussian whose width is equal to SIG_GAU times the
! current smoothed frequency. Thus cross-section area is not conserved
! exactly. Allows wider freqeuncy spacing as the frequency is increased.
!
! On entry
!         NU            Should contain the frequency in units of nu_zero
!         CROSS         The cross section.
!         NCROSS        Number of cross sections.
!         SIG_GAU       The sigma of the gaussian smoothing routine.
!         DEL_NU        Frequency spacing for smoothed cross-section.
!                         (should be fraction of SIG_GAU).
!         CUT_ACCURACY  Accuracy to limit omission of data points.
!                          Assumes inear interpolation.
!         ABOVE         Smooth from above threshold only?
!
! On exit
!        NU_SM       Contains the evenly spaced frequencies at which
!                    the cross-section is tabulated.
!        CROSS_SM    Smoothed Cross-Section
!        NSM         Number of frequency points in smooth cross-section.
!
! Altered 24-May-2005 : Fixed to handle case where mirroring of cross-section
!                        about NU=1 gives negative frequencies.
! Created 19-Apr-2004 : Based on SM_PHOT_V2
!                       Designed to be incorporated into CMFGEN
!                       Some cleaning.
!			See SM_PHOT_V2 for earlier changes.
!
	IMPLICIT NONE
	INTEGER NCROSS
!
	INTEGER NSM_MAX
	REAL(KIND=LDP) NU(NSM_MAX)
	REAL(KIND=LDP) CROSS(NSM_MAX)
!
	REAL(KIND=LDP) SIG_GAU_KMS
	REAL(KIND=LDP) FRAC_SIG_GAU
	REAL(KIND=LDP) CUT_ACCURACY
	LOGICAL ABOVE
!
! Work arrays.
!
	REAL(KIND=LDP) NU_SM(NSM_MAX)
	REAL(KIND=LDP) CROSS_SM(NSM_MAX)
!
	LOGICAL FLAG
	EXTERNAL EQUAL
	LOGICAL EQUAL
!
! Local variables
!
	INTEGER, PARAMETER :: NZ=10000
	INTEGER, PARAMETER :: NLOC=1000000
	REAL(KIND=LDP) Z(NZ)
	REAL(KIND=LDP) NU_FINE(NLOC)
	REAL(KIND=LDP) CROSS_FINE(NLOC)
	INTEGER NFINE
!
	REAL(KIND=LDP) T1,T2,DIFF,SUM,SIG_NEW
	REAL(KIND=LDP) LAST_FREQ
	REAL(KIND=LDP) DNU
	REAL(KIND=LDP) dCROSS
	REAL(KIND=LDP) SIG_GAU
	INTEGER I,J,K,ML
	INTEGER IST,IEND,IBEG,N_INS
	INTEGER NSM
!
	INTEGER NU_INDX_ST
	REAL(KIND=LDP) NU_ST
!
	REAL(KIND=LDP) MAX_DEL
	REAL(KIND=LDP) PTS_PER_SIG
!
	PTS_PER_SIG=4.0_LDP			!Used for the convolution.
	SIG_GAU=SIG_GAU_KMS/2.998E+05_LDP		!C does not need to be accurate
	DNU=FRAC_SIG_GAU*SIG_GAU
!
! Check grid is mononotnic.
!
	DO I=2,NCROSS
	  IF(NU(I) .LT. NU(I-1) .OR. NU(I-1) .LT. 0)THEN
	    WRITE(6,*)'Original (but modified) frequency grid is not monotonic or has -ve values'
	    WRITE(6,*)'See unit 30 for information on the unmodified coarse grid'
	    WRITE(6,*)I,NU(I-1),NU(I)
	    WRITE(6,*)IBEG,NCROSS
	    DO J=1,NCROSS-1
	      WRITE(30,'(I6,3ES16.6)')J,NU(J),3.0D+05*(NU(J+1)/NU(J)-1.0D0),CROSS(J)
	    END DO
	    J=NCROSS;  WRITE(30,'(I6,3ES16.6)')J,NU(J),0.0D0,CROSS(J)
!	    STOP
	  END IF
	END DO
!
! If ABOVE is TRUE, we ignore cross-section (except for point immediately
! proceeding threshold) when smoothing cross-section.  We reflect the cross
! section about the threshold frequency. This ensures that we conserve area.
! If ABOVE is FALSE, we smooth across the threshold.
!
	IF(ABOVE)THEN
!
! Use the FINE arrays as temporary vectors as we shift the
! coarse frequency grid. We add extra points near the origin
! to ensure that the rfelection about NU(1) does extend
! to far below NU=1 (and to -ve values).
!
	  K=1
	  NU_FINE(1)=NU(1)
	  CROSS_FINE(1)=CROSS(1)
	  DO I=2,NCROSS
	    T1=(NU(I)-NU(I-1))/0.05_LDP
	    IF(K+1+T1 .GT. NLOC)THEN
	      WRITE(6,*)'Error in SM_PHOT_V3'
	      WRITE(6,*)'NLOC is too small'
	      WRITE(6,'(1X,A,I7,6X,A,I7)')'NLOC-',NLOC,'NCROSS=',NCROSS
	    STOP
	    END IF
	    IF(T1 .GT. 1.0_LDP .AND. NU_FINE(MAX(1,K-1)) .LT. 2.0_LDP)THEN
	      J=T1
	      T1=(NU(I)-NU(I-1))/(J+1)
	      dCROSS=(CROSS(I)-CROSS(I-1))/(J+1)
	      DO ML=1,J
	        K=K+1
	        NU_FINE(K)=NU_FINE(K-1)+T1
	        CROSS_FINE(K)=CROSS_FINE(K-1)+dCROSS
	      END DO
	    END IF
	    K=K+1
	    NU_FINE(K)=NU(I)
	    CROSS_FINE(K)=CROSS(I)
	  END DO
	  NCROSS=K
	  NFINE=NCROSS
!
! How we proceed depends on whether NU=1 is present.
!
	  I=1
	  DO WHILE(NU_FINE(I) .LT. 1.0_LDP)
	    I=I+1
	  END DO
	  IBEG=I
	  IF(IBEG .NE. 1)THEN
	    K=IBEG-1
	    IF(ABS(NU_FINE(K)-1.0_LDP) .LT. 1.0E-07_LDP)THEN
	      NU_FINE(K)=1.0_LDP
	      IBEG=K
	    END IF
	  END IF
	  IF(ABS(NU_FINE(IBEG)-1.0_LDP) .LT. 1.0E-07_LDP)NU_FINE(IBEG)=1.0_LDP
!
! Find range of Gaussian.
!
	  DO WHILE(NU_FINE(I) .LT. 1.0_LDP+6.0_LDP*SIG_GAU)
	    I=I+1
	  END DO
	  N_INS=I-IBEG+1
!
	  IF(NU_FINE(IBEG) .NE. 1.0_LDP)THEN
	    IF(IBEG .EQ. 1)THEN
	      T2=0.0_LDP
	    ELSE
	      T1=(NU_FINE(IBEG-1)-1.0_LDP)/(NU_FINE(IBEG-1)-NU_FINE(IBEG))
	      T2=(1.0_LDP-T1)*CROSS_FINE(IBEG-1)+T1*CROSS_FINE(IBEG)
	    END IF
	    N_INS=N_INS+1
	    DO J=NCROSS,IBEG,-1
	      K=J+N_INS-IBEG+1
	      NU(K)=NU_FINE(J)
	      CROSS(K)=CROSS_FINE(J)
	    END DO
	    NU(N_INS)=1.0_LDP
	    CROSS(N_INS)=T2
	    DO J=N_INS-1,1,-1
	      K=IBEG+(N_INS-1-J)
	      NU(J)=2.0_LDP-NU_FINE(K)
	      CROSS(J)=CROSS_FINE(K)
	    END DO
	    NCROSS=NCROSS-IBEG+1+N_INS
!
	  ELSE
	    DO J=NCROSS,IBEG,-1
	      K=J+N_INS-IBEG+1
	      NU(K)=NU_FINE(J)
	      CROSS(K)=CROSS_FINE(J)
	    END DO
	    DO J=1,N_INS
	      K=IBEG+(N_INS+1-J)
	      NU(J)=2.0_LDP-NU_FINE(K)
	      CROSS(J)=CROSS_FINE(K)
	    END DO
	    NCROSS=NCROSS-IBEG+1+ N_INS
	  END IF
	ELSE
!
! Usually we will have enough points below threshold. For simplicty,
! we assume a constant cross-section below these points.
!
	 IF(NU(1) .GT. 1.0_LDP-6.0_LDP*SIG_GAU)THEN
	   DO J=NCROSS,1,-1
	     K=J+1
	     NU(K)=NU(J)
	     CROSS(K)=CROSS(J)
	   END DO
	   NU(1)=1.0_LDP-6.0_LDP*SIG_GAU
	   CROSS(1)=CROSS(2)
	   NCROSS=NCROSS+1
	  END IF
	END IF
!
	IF(NCROSS+1 .GT. NSM_MAX)THEN
	  WRITE(6,*)'Error in SM_PHOT_V3'
	  WRITE(6,*)'NSM_MAX is too small'
	  STOP
	END IF
!
! Extend cross-sections at high frequencies. We ensure cross-section
! falls off as least as 1/nu.
!
	NCROSS=NCROSS+1
	NU(NCROSS)=NU(NCROSS-1)*(1.0_LDP+6.0_LDP*SIG_GAU)
	T1=LOG(NU(NCROSS-1)/NU(NCROSS-2))
	T2=LOG(CROSS(NCROSS-1)/CROSS(NCROSS-2))/T1
	IF(T2 .GT. -1.0_LDP)T2=-1.0_LDP
	CROSS(NCROSS)=CROSS(NCROSS-1)*(NU(NCROSS)/NU(NCROSS-1))**T2
!
! Check grid is mononotnic.
!
	DO I=2,NCROSS
	  IF(NU(I) .LT. NU(I-1) .OR. NU(I-1) .LT. 0)THEN
	    WRITE(6,*)'Original (but modified) frequency grid is not monotonic or has -ve values'
	    WRITE(6,*)'See unit 30 for information on the unmodified coarse grid'
	    WRITE(6,*)'See unit 31 for information on the (modified) coarse grid'
	    WRITE(6,*)I,NU(I-1),NU(I)
	    WRITE(6,*)IBEG,N_INS
	    DO J=1,NFINE-1
	      WRITE(30,'(I6,3ES16.6)')J,NU_FINE(J),
	1          3.0D+05*(NU_FINE(J+1)/NU_FINE(J)-1.0D0),CROSS_FINE(J)
	    END DO
	    DO J=1,NCROSS-1
	      WRITE(31,'(I6,3ES16.6)')J,NU(J),3.0D+05*(NU(J+1)/NU(J)-1.0D0),CROSS(J)
	    END DO
!	    STOP
	  END IF
	END DO
!
!	CALL DP_CURVE(NCROSS,NU,CROSS)
!
! Interpolate cross-section onto a sufficiently fine grid so that the Gaussian
! correctly samples the cross-section data. At present we assume linear
! interpolation. This is satisfactory given the accuracy of the cross-sections.
!
! NB: The Gaussian is assumed to have constant width in velocity space.
!     For this reason SIG_GAU is mutipled by NU(ML). For frequencies < 1,
!     we use SIG_GAU. This prevents DIFF from becoming arbitrarily small if
!     NU should approach zero.
!
	I=1
	NU_FINE(1)=NU(1)
	CROSS_FINE(1)=CROSS(1)
	DO ML=2,NCROSS
	  MAX_DEL=SIG_GAU*NU(ML)/PTS_PER_SIG
	  IF(NU(ML) .LE. 1.0_LDP)MAX_DEL=SIG_GAU/PTS_PER_SIG
	  DIFF=NU(ML)-NU(ML-1)
	  IF( DIFF .GT. MAX_DEL)THEN
            J=DIFF/MAX_DEL+1
	    T2=DIFF/J
	    DO K=1,J-1
	      IF(I+K .GT. NLOC)GOTO 999
	      NU_FINE(I+K)=NU_FINE(I)+K*T2
	      T1=(NU_FINE(I+K)-NU(ML-1))/(NU(ML)-NU(ML-1))
	      CROSS_FINE(I+K)=(1.0_LDP-T1)*CROSS(ML-1)+T1*CROSS(ML)
	    END DO
	    I=I+J-1
	  END IF
	  I=I+1
	  IF(I .GT. NLOC)GOTO 999
	  NU_FINE(I)=NU(ML)
	  CROSS_FINE(I)=CROSS(ML)
	END DO
	NFINE=I
!
!	CALL DP_CURVE(NFINE,NU_FINE,CROSS_FINE)
!
! Define the smooth mesh. DNU is the frequency spacing for the
! smoothed cross-section output.
!
	FLAG=.TRUE.
	NSM=1
	NU_SM(1)=1.0_LDP
	LAST_FREQ=NU(NCROSS)
	DO WHILE(FLAG)
	  IF(NSM+1 .GT. NSM_MAX)THEN
	    WRITE(6,10)NSM,NSM_MAX
 10         FORMAT(' Error NSM_MAX TOO small : NSM =',I6,': NSM_MAX =',I6)
	    WRITE(6,*)'DNU=',DNU
	    STOP
	  END IF
	  NU_SM(NSM+1)=NU_SM(1)*(1.0_LDP+DNU)**NSM
	  IF(NU_SM(NSM+1) .GT. LAST_FREQ)EXIT
	  NSM=NSM+1
	END DO
C 
C Now determine smooth cross section for each (smooth) frequency.
C
C We use Z to generate for gaussian over the desired freqency band.
C We also integrate Z to get the proper normailization.
C
C IST  is index for ML at which left side of Gaussian begins.
C IEND is index for ML at which right side of Gaussian ends.
C
	IST=1
	IEND=2
	DO ML=1,NSM
	  SIG_NEW=SIG_GAU*NU_SM(ML)
	  DO WHILE (NU_FINE(IST) .LT. NU_SM(ML)-5.0_LDP*SIG_NEW)
	    IST=IST+1
	  END DO
	  DO WHILE ( (NU_FINE(IEND)-NU_SM(ML)) .LT. 5.0_LDP*SIG_NEW)
	    IF(IEND .EQ. NFINE)GOTO 100
	    IEND=IEND+1
	  END DO
100	CONTINUE
!
! Determine the smoothing Gaussian.
!
	  IF(IEND-IST+1 .GT. NZ)THEN
	    WRITE(6,*)'Error in SM_PHOT_V3: NZ too small'
	    WRITE(6,*)'NZ =',NZ
	    WRITE(6,*)'IEND-IST+1 =',IST,IEND,IEND-IST+1
	    WRITE(6,*)NU_FINE(IST),NU_FINE(IEND),SIG_NEW
	    WRITE(6,*)SIG_GAU_KMS,FRAC_SIG_GAU
	    WRITE(6,*)'See unit 30 for information on the fine grid'
	    WRITE(6,*)'See unit 31 for information on the (modified) coarse grid'
	    DO I=1,NFINE-1
	      WRITE(30,*)I,NU_FINE(I),3.0D+05*(NU_FINE(I+1)/NU_FINE(I)-1.0D0)
	    END DO
	    DO I=1,NCROSS-1
	      WRITE(31,'(I6,3ES16.6)')I,NU(I),3.0D+05*(NU(I+1)/NU(I)-1.0D0),CROSS(I)
	    END DO
	    STOP
	  END IF
	  DO I=IST,IEND
	    J=I-IST+1
            Z(J)=EXP( -0.5_LDP*( (NU_SM(ML)-NU_FINE(I))/SIG_NEW )**2 )
	  END DO
!
! Perform the integration for the current frequency.
!
	  CROSS_SM(ML)=0.0_LDP
	  SUM=0.0_LDP
	  DO I=IST,IEND-1
	    J=I-IST+1
	    CROSS_SM(ML)=CROSS_SM(ML)+ (NU_FINE(I+1)-NU_FINE(I))*
	1         (CROSS_FINE(I)*Z(J)+CROSS_FINE(I+1)*Z(J+1))
	    SUM=SUM + (NU_FINE(I+1)-NU_FINE(I))*(Z(J)+Z(J+1))
	  END DO
	  IF(SUM .NE. 0.0_LDP)THEN
	    CROSS_SM(ML)=CROSS_SM(ML)/SUM
	  ELSE
	    WRITE(6,*)'Zero SUM in SM_PHOT_V3'
	    WRITE(6,*)IST,IEND,ML,NSM
	    WRITE(6,*)'NU_SM(ML)=',NU_SM(ML)
	    WRITE(6,*)'NU_FINE(IST)=',NU_FINE(IST)
	    WRITE(6,*)'NU_FINE(IEND)=',NU_FINE(IEND)
	    WRITE(6,*)'SIG_GAU=',SIG_GAU
	    WRITE(6,*)'See unit 30 for information on the fine grid'
	    WRITE(6,*)'See unit 31 for information on the (modified) coarse grid'
	    WRITE(30,*)NFINE
	    DO J=1,NFINE
	      WRITE(30,'(I6,2ES16.6)')J,NU_FINE(J),CROSS_FINE(J)
	    END DO
	    DO J=1,NCROSS
	      WRITE(31,'(I6,2ES16.6)')J,NU(J),CROSS(J)
	    END DO
	    STOP
	  END IF
	END DO				!Frequency loop
!
! Omit pints that are uneessary to maintain an accuracy of CUT_ACCURACY in the
! cross-section assuming linear interpolation.
!
        CALL CUT_POINTS_V3(NU,CROSS,NCROSS,NU_SM,CROSS_SM,NSM,CUT_ACCURACY)
!
	RETURN
!
999	WRITE(6,*)'Error --- NLOC too small in SM_PHOT_V3'
	WRITE(6,*)'NU(1)=',NU(1)
	WRITE(6,*)'NU(NCROSS)=',NU(NCROSS)
	WRITE(6,*)'DNU=',DNU
!
	STOP
	END
