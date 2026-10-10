!
! Simple program to create an INCIDENT INTENSITY file for use with cmfgen and
! the plane-parallel subroutine:
!                              pp_form_cmf_v2.f
! The incident intensity is outout in table format, and thus arbitrary
! distributions are possible. At present, the blackbody flux is described
! a blackbody distribution.
!
	PROGRAM CREATE_INCID_INTEN
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	REAL(KIND=LDP), ALLOCATABLE :: MU(:)
	REAL(KIND=LDP), ALLOCATABLE :: DIST(:)
	REAL(KIND=LDP) NU
	REAL(KIND=LDP) BNU
!
	REAL(KIND=LDP) TEFF_STAR
	REAL(KIND=LDP) TEFF_INCID
	REAL(KIND=LDP) REL_LUM
	REAL(KIND=LDP) DILUTION_FACTOR
!
	REAL(KIND=LDP), PARAMETER :: HDKT=4.7994145_LDP
	REAL(KIND=LDP), PARAMETER :: TWOHCSQ=0.0147452575_LDP
	REAL(KIND=LDP), PARAMETER :: NU_MAX=100.0_LDP
	REAL(KIND=LDP), PARAMETER :: NU_MIN=0.0001_LDP
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) FLUX
	REAL(KIND=LDP) X
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) NORMALIZED_ANGLE_INTEGRAL
!
	INTEGER I
	INTEGER NANG
	INTEGER NCF
	INTEGER NPTS_PER_DECADE
!
	CHARACTER*132 FILENAME
!
	NANG=11; NPTS_PER_DECADE=20; REL_LUM=0.1_LDP
	CALL GEN_IN(TEFF_STAR,'Teff for star (in 10^4K)')
	CALL GEN_IN(TEFF_INCID,'Teff describing incident radiation (in 10^4 K)')
	CALL GEN_IN(REL_LUM,'Relative fluxes (incident/stellar)')
	CALL GEN_IN(NPTS_PER_DECADE,'Number of points per decade in frequency')
	CALL GEN_IN(NANG,'Number of points in angle')
!
	ALLOCATE (MU(NANG))
	ALLOCATE (DIST(NANG))
!
! Set angles used to describe the incident radiation field.
!
	MU(1)=0.0_LDP
	MU(NANG)=1.0_LDP
	DO I=2,NANG-1
	  MU(I)=(I-1.0_LDP)/(NANG-1.0_LDP)
	END DO
!
! Set the distribution of flux as a function of angle.
!
	DO I=1,NANG
	  DIST(I)=1.0_LDP-4.0_LDP*(MU(I)-0.5_LDP)**2
	END DO
!
! Estimate the integral I(mu).mu dmu. For isotropic radiation, the integral
! is 0.5. We thus normalize the integral by 0.5.
!
	T1=0.0_LDP
	DO I=1,NANG-1
	  T1=T1+0.5_LDP*(MU(I+1)-MU(I))*(MU(I)*DIST(I)+MU(I+1)*DIST(I+1))
	END DO
	NORMALIZED_ANGLE_INTEGRAL=2.0_LDP*T1
        WRITE(6,*)'Normalized angle weighting integral is',NORMALIZED_ANGLE_INTEGRAL
!
! Determine the dilution factor to give the requested relative luminosity fir
! the requested incident Teff (fluxes) and the stellar model Teff.
!
	DILUTION_FACTOR=(REL_LUM/NORMALIZED_ANGLE_INTEGRAL)*(TEFF_STAR/TEFF_INCID)**4
!
! Determine number of frequency points.
!
	NCF=LOG10(NU_MAX/NU_MIN)*NPTS_PER_DECADE+1
!
! INCID_INTEN is the name of the output file to contain the incident fluxes as a
! function of angle and frequency.
!
	FILENAME='INCID_INTEN'
	CALL GEN_IN(FILENAME,'Output file')
	OPEN(UNIT=15,STATUS='UNKNOWN',FILE=TRIM(FILENAME))
!
	WRITE(15,'(A)')'! '
	WRITE(15,'(A)')'! The effective temperature of the incident radiation is defined as the'
	WRITE(15,'(A)')'! temperature of the blackbody that gives the equivalent radiation flux.'
	WRITE(15,'(A)')'! Teff(star) is used to set the relative luminosity and dilution factor.'
	WRITE(15,'(A)')'! Thus the dilution factor (W) is defined by:'
	WRITE(15,'(A)')'! '
	WRITE(15,'(A)')'!         Rel. Lum = W . NAI . (Teff/Teff[*])*0.25'
	WRITE(15,'(A)')'! '
	WRITE(15,'(A)')'! where NIA is the normalized angle integral. The actual dilution factor'
	WRITE(15,'(A)')'! used to scale the blackbody is output just prior to the angle grid.'
	WRITE(15,'(A)')'! This can be easily changed without recreating a new incident flux file.'
	WRITE(15,'(A)')'! Frequency is in units of 10^15 Hz, and I is in CGS units'
	WRITE(15,'(A)')'! '
	WRITE(15,'(A,ES14.6)')'! Teff(K) of incident radiation      :',TEFF_INCID*1.0D+04
	WRITE(15,'(A,ES14.6)')'! Normalized angle weighting integral:',NORMALIZED_ANGLE_INTEGRAL
	WRITE(15,'(A,ES14.6)')'! '
	WRITE(15,'(A,ES14.6)')'! Assumed Teff(K) of star is:         ',TEFF_STAR*1.0D+04
	WRITE(15,'(A,ES14.6)')'! Relative luminosity:                ',REL_LUM
	WRITE(15,'(A,ES14.6)')'! Dilution factor is:                 ',DILUTION_FACTOR
	WRITE(15,*)' '
!
	WRITE(15,'(A,T30,A)')'   10-Apr-2006','!Format date'
	WRITE(15,'(I14,T30,A)')NCF,'!Number of continuum frequencies'
	WRITE(15,'(I14,T30,A)')NANG,'!Number of angles'
	WRITE(15,'(ES14.6,T30,A)')DILUTION_FACTOR,'!Dilution factor used for scaling'
	WRITE(15,*)' '
	WRITE(15,*)'Angle grid'
	WRITE(15,*)(MU(I),I=1,NANG)
	WRITE(15,*)' '
	WRITE(15,*)'Intensity variation with angle'
	WRITE(15,*)(DIST(I),I=1,NANG)
!
	WRITE(15,*)' '
	T1=EXP(LOG(10.0_LDP)/NPTS_PER_DECADE)
	NU=NU_MAX*T1
	FLUX=0.0_LDP
	DO I=1,NCF
	  NU=NU/T1
	  X=EXP(-HDKT*NU/TEFF_INCID)
	  BNU=TWOHCSQ*NU*NU*NU*X/(1.0_LDP-X)
	  FLUX=FLUX+BNU*(NU*T1-NU/T1)/2.0_LDP
	  WRITE(15,'(X,2ES14.6)')NU,BNU
	END DO
	PI=4.0_LDP*ATAN(1.0_LDP)
	T1=(FLUX*PI*1.0E+15_LDP/5.6705E-05_LDP)**(0.25_LDP)
	WRITE(6,'(A,ES14.6,A)')'Equivalent effective temperature is',T1,' K'
!
	STOP
	END
