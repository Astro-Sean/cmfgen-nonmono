!
! Subroutine to read in non-local rdiactive energy deposition.
! The data should be computed by a separate calculation, and is
! read in from:
!                current_nonlocal_decay_energy.dat
!
	SUBROUTINE GET_NON_LOCAL_GAMMA_ENERGY_V2(R,V,ND,LU)
	USE SET_KIND_MODULE
	USE CONTROL_VARIABLE_MOD
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
! Altered 20-May-2023 - This version had not been updated in main source. Made two version compatible.
! Altered 21-Nov-2021 - Added option to read in GAMRAY_ENERGY_DEP if present and when
!                          GAMRAY_TRANS='RAD_TRANS'
! Altered 29-Jan-2014 - Call DO_GAM_ABS_APPROX_V2
!                         Compute emitted and absorbed decau luminosities.
! Altered 05-Jan-2014 - R added to call.
!                         Now output energy emitted & absorbed to check_edep.
!                         Changed to V2.
!
	INTEGER ND
	INTEGER LU
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
!
! Variables for reading in the non-local energy deposition from outside file
!
	INTEGER NDTMP
	REAL(KIND=LDP), ALLOCATABLE :: VTMP(:),EDEPTMP(:),EDEPNEW(:)
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) CONV_CONST
	REAL(KIND=LDP) WRK(ND)
	INTEGER LUER
	INTEGER I
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=80) STRING
	LOGICAL FILE_EXISTS
!!
	LUER=ERROR_LU()
        IF (GAMRAY_TRANS .EQ. 'LOCAL') THEN
!
        ELSE IF (GAMRAY_TRANS.EQ.'ABS_TRANS') THEN
!
! Before the call, WRK must contain the locally EMITTED energy. After the call,
! RADIOACTIVE_DECAY_ENERGY will contain the locally ABSORBED energy.
!
	   WRK=RADIOACTIVE_DECAY_ENERGY
	   CALL DO_GAM_ABS_APPROX_V2(RADIOACTIVE_DECAY_ENERGY,WRK,KINETIC_DECAY_ENERGY,ND)
!
	ELSE IF (GAMRAY_TRANS .EQ. 'NONLOCAL')THEN
!
	   CLOSE(LU)
	   OPEN(LU,FILE='current_nonlocal_decay_energy.dat',STATUS='UNKNOWN')
!
! Can now have comments at top of file.
!
	   DO WHILE(INDEX(STRING,'Number of depth points') .EQ. 0)
	     READ(LU,'(A)',IOSTAT=IOS)STRING
	     IF(IOS .NE. 0)THEN
	       WRITE(6,*)'Error reading number of depth points from current_nonlocal_decay_energy.dat'
	       STOP
	     END IF
	   END DO
!
	   READ(STRING,*) NDTMP
	   READ(LU,*) T1
	   IF (T1 .NE. SN_AGE_DAYS) THEN
	      WRITE(LUER,'(A80)') 'New time for Gamma-ray transport calculation is incompatible'
	      WRITE(LUER,'(A50)') 'Monte Carlo new time [days]',T1
	      WRITE(LUER,'(A50)') 'CMFGEN new time [days]',SN_AGE_DAYS
	      STOP
	   ENDIF
!
	   ALLOCATE (VTMP(NDTMP),EDEPTMP(NDTMP),EDEPNEW(ND))
	   DO I=1,NDTMP
	      READ(LU,*) VTMP(I),EDEPTMP(I)
	   ENDDO
	   CLOSE(LU)
!
! Allow for rounding error at the outer boundary.
!
	   IF(V(1) .GT. VTMP(1))THEN
	     VTMP(1)=VTMP(1)+0.01_LDP*(VTMP(1)-VTMP(2))
	     IF(V(1) .GT. VTMP(1))THEN
	       WRITE(6,*)'Error in interpolating energy deposition in get_non_local_gamma_energy_v2'
	       WRITE(6,*)'VTMP(1) in file is to small: VTMP(1)=',VTMP(1)
	       WRITE(6,*)'V(1) in model is',V(1)
	       STOP
	     END IF
	   END IF
	   VTMP(NDTMP) = V(ND)
	   CALL LIN_INTERP(V,EDEPNEW,ND,VTMP,EDEPTMP,NDTMP)
!
	   OPEN(LU,FILE='check_edep.dat',STATUS='UNKNOWN')
	   WRITE(LU,'(I5,A50)') ND,' !Number of depth points'
!
!
! The conversion factor:
!   (a) 4pi r^2 dr --> 4pi .10^30 (as R is in units of 10^10 cm).
!   (b) The extra factor of 4 arises as ATAN(1.0D0) is pi/4.
!   (c) /Lsun to convert to solar luminosities.
!
	   CONV_CONST=16.0_LDP*ATAN(1.0_LDP)*1.0E+30_LDP/3.826E+33_LDP
	   DO I=1,ND
	     WRK(I)=RADIOACTIVE_DECAY_ENERGY(I)*R(I)*R(I)
	   END DO
	   CALL LUM_FROM_ETA(WRK,R,ND)
	   T2=SUM(WRK)*CONV_CONST
	   WRITE(LU,'(A)')'!'
	   WRITE(LU,'(A,ES13.5,A)')'!            Radioactive energy emitted is: ',T2,' Lsun'
	   DO I=1,ND
	     WRK(I)=EDEPNEW(I)*R(I)*R(I)
	   END DO
	   CALL LUM_FROM_ETA(WRK,R,ND)
	   T1=SUM(WRK)*CONV_CONST
	   WRITE(LU,'(A,ES13.5,A)')'!            Radioactive energy absorbed is:',T1,' Lsun'
	   WRITE(LU,'(A,ES13.5,A)')'!Fraction of radioactive energy absorbed is:',T1/T2
	   WRITE(LU,'(A)')'!'
!
! Switched ordering of EDEPNEW and RADIOACTIVE_DECAY_ENERGY to be consistent
! with that in 'current_nonlocal_decay_energy.dat'
!
	   WRITE(LU,'(A,12X,A,5X,A,5X,A)')'!','Velocity','Non-Local Decay','    local Decay'
	   WRITE(LU,'(A,12X,A,5X,A,5X,A)')'!','  km/s  ','    ergs/cm^3/s','   ergs/cm^3/s'
	   DO I=1,ND
	      WRITE(LU,'(3ES20.8)')V(I),EDEPNEW(I),RADIOACTIVE_DECAY_ENERGY(I)
	   ENDDO
!
	   WRITE(LU,'(/,/,I5,A50)') NDTMP,' !Number of depth points in MC computation'
	   DO I=1,NDTMP
	      WRITE(LU,'(3ES20.8)') VTMP(I),EDEPTMP(I)
	   ENDDO
	   WRITE(LU,'(A)')'!'
!
	   CLOSE(LU)
!
! We overwrite the edep computed before assuming local deposition
!
	   RADIOACTIVE_DECAY_ENERGY = EDEPNEW
	   DEALLOCATE (VTMP,EDEPTMP,EDEPNEW)
!
	ELSE IF (GAMRAY_TRANS .EQ. 'RAD_TRANS') THEN
!
! Check to see if we have already computed the GAMMA deposition profile as a function
! of radius. Do not need to change the keyword in VADAT -- works even if model is
! restarted.
!
	  INQUIRE(FILE='GAMRAY_ENERGY_DEP',EXIST=FILE_EXISTS)	
	  IF(FILE_EXISTS)THEN
	    OPEN(UNIT=LU,FILE='GAMRAY_ENERGY_DEP',STATUS='OLD',ACTION='READ')
	    T1=-1.0_LDP; NDTMP=-1
	    DO WHILE(1 .EQ. 1)
	      READ(LU,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Error reading number of depth points or age from GAMRAY_ENERGY_DEP'
	        STOP
	      END IF
	      IF(INDEX(STRING,'Number of depth points') .NE. 0)THEN
	        READ(STRING,*) NDTMP
	      ELSE IF(INDEX(STRING,'Current time after explosion') .NE. 0)THEN
	        READ(STRING,*)T1
	      END IF
	      IF(T1 .GT. 0 .AND. NDTMP .GT. 0)EXIT
	    END DO
	    IF(T1 .LT. 0 .OR. NDTMP .LT. 0)THEN
	      WRITE(6,*)'Error reading number of depth points or age from GAMRAY_ENERGY_DEP'
              WRITE(6,*)'NDTEMP=',NDTMP; WRITE(6,*)'AGE=',T1
	      STOP
	    END IF
!
	    IF (ABS( (T1-SN_AGE_DAYS)/(T1+SN_AGE_DAYS) ) .GT. 1.0E-07_LDP) THEN
	      WRITE(LUER,'(A)') 'New time for Gamma-ray transport calculation is incompatible'
	      WRITE(LUER,'(A,ES14.8)')'Monte Carlo new time [days]',T1
	      WRITE(LUER,'(A,ES14.8)')'CMFGEN new time [days]',SN_AGE_DAYS
	      STOP
	    ENDIF
!
! Skip additional header info.
!
	    STRING='!'
	    DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	      READ(LU,'(A)',IOSTAT=IOS)STRING
	    END DO
	    BACKSPACE(LU)
!
! We may need to inteprolate because of automatic changes to the grid.
! We interpolate in V, but could interpolate in R.
!
	    ALLOCATE (VTMP(NDTMP),EDEPTMP(NDTMP),EDEPNEW(ND))
	    DO I=1,NDTMP
	       READ(LU,*)T1,VTMP(I),EDEPTMP(I)
	    ENDDO
	    IF( (VTMP(1)-V(1))/V(1) .LT. 1.0E-06_LDP)VTMP(1)=V(1)
	    IF( (V(ND)-VTMP(NDTMP))/V(ND) .LT. 1.0E-06_LDP)VTMP(NDTMP)=V(ND)
	    CALL LIN_INTERP(V,RADIOACTIVE_DECAY_ENERGY,ND,VTMP,EDEPTMP,NDTMP)
	    DEALLOCATE (VTMP,EDEPTMP,EDEPNEW)
	    CLOSE(LU)
	  ELSE
	    RADIOACTIVE_DECAY_ENERGY=0.0_LDP
	  END IF
!
	ELSE
	   WRITE(LUER,*)'Error in GET_NON_LOCAL_GAMMA_ENERGY'
	   WRITE(LUER,*)'Unrecognized gamma-ray transport option'
	   WRITE(LUER,*)GAMRAY_TRANS
	   STOP
	ENDIF
!
	RETURN
	END
