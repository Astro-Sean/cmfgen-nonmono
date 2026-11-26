!
! Subroutine to compare terms in the Hydronamical equations. Routine is
! meant for illustration purposes only.
!
	SUBROUTINE HYDRO_TERMS_V2(POP_ATOM,R,V,T,SIGMA,ED,LUM_STAR,
	1              STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1              FLUX_MEAN_OPAC,ELEC_MEAN_OPAC,
	1              LAM_FLUXMEAN_BAND_END,BAND_FLUXMEAN,
	1              BAND_FLUX,N_FLUXMEAN_BANDS,
	1              LU_OUT,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 16-Jan-2005 : Accuracy of dPdR/ROH improved. Now use monotonic
!                         interpolation. Depth index added output.
! Altered 17-Mar-2001 : Altered V output format
! Altered 27-Oct-2000 : Percentage error improved.
!                       STARS mass now output.
! Altered 24-May-1996 : L GEN_ASCI_OPEN and ERROR_LU installed.
! Altered 06-Mar-1996 : Description of terms added to output.
! Created 28-Mar-1994
!
	INTEGER ND,LU_OUT
!
	REAL(KIND=LDP) POP_ATOM(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) LUM_STAR(ND)
!
	REAL(KIND=LDP) FLUX_MEAN_OPAC(ND)
	REAL(KIND=LDP) ELEC_MEAN_OPAC(ND)
!
	INTEGER N_FLUXMEAN_BANDS
	REAL(KIND=LDP) LAM_FLUXMEAN_BAND_END(ND)
	REAL(KIND=LDP) BAND_FLUXMEAN(ND,N_FLUXMEAN_BANDS)
	REAL(KIND=LDP) BAND_FLUX(ND,N_FLUXMEAN_BANDS)
!
	REAL(KIND=LDP) STARS_MASS
	REAL(KIND=LDP) MEAN_ATOMIC_WEIGHT
!
! External functions
!
	REAL(KIND=LDP) BOLTZMANN_CONSTANT
	REAL(KIND=LDP) GRAVITATIONAL_CONSTANT
	REAL(KIND=LDP) MASS_SUN
	REAL(KIND=LDP) LUM_SUN
	REAL(KIND=LDP) FUN_PI
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) ATOMIC_MASS_UNIT
!
	EXTERNAL BOLTZMANN_CONSTANT,GRAVITATIONAL_CONSTANT,MASS_SUN
	EXTERNAL LUM_SUN,FUN_PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT
!
	REAL(KIND=LDP) MOD_PRESSURE(ND)
	REAL(KIND=LDP) COEF(ND,4)
!
! Local variables.
!
	REAL(KIND=LDP) VdVdR
	REAL(KIND=LDP) dPdR_ON_ROH
	REAL(KIND=LDP) g_TOT
	REAL(KIND=LDP) g_ELEC
	REAL(KIND=LDP) g_GRAV
	REAL(KIND=LDP) g_RAD
	REAL(KIND=LDP) ERROR
	REAL(KIND=LDP) RLUMST(ND)
!
	REAL(KIND=LDP) GRAV_CON
	REAL(KIND=LDP) dP_CON
	REAL(KIND=LDP) RAD_CON
	REAL(KIND=LDP) T1
!
	CHARACTER*80 FMT
	INTEGER I,J,IOS,ERROR_LU
	EXTERNAL ERROR_LU
!
	I=MAX(132,(N_FLUXMEAN_BANDS+1)*12+16)
	CALL GEN_ASCI_OPEN(LU_OUT,'HYDRO','UNKNOWN',' ',' ',I,IOS)
	IF(IOS .NE. 0)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error opening HYDRO in HYDRO_TERMS'
	  WRITE(I,*)'IOS=',IOS
	END IF
!
! Output file header.
!
	WRITE(LU_OUT,'(1X,6X,A,6X, 8X,A,3X, 2X,A, 5(5X,A,1X), 4X,A,2X,A)')
	1 'R','V','% Error','   VdVdR',
	1                   'dPdR/ROH',
	1                   '   g_TOT',
	1                   '   g_RAD',
	1                    '  g_ELEC','Gamma','Depth'
!
! 1.0D-06 = [ 10^(-10) {from 1/r) * 10^4 {from T}
	dP_CON=1.0E-06_LDP*BOLTZMANN_CONSTANT()/MEAN_ATOMIC_WEIGHT/
	1               ATOMIC_MASS_UNIT()
!
	GRAV_CON=1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()
!
!       1.0E-10 (CHI) / 1.0E-20 (1/r^2)
	RAD_CON=1.0E-30_LDP*LUM_SUN()/ATOMIC_MASS_UNIT()/4.0_LDP/FUN_PI()/
	1                  SPEED_OF_LIGHT()/MEAN_ATOMIC_WEIGHT
!
	DO I=1,ND
	  MOD_PRESSURE(I)=(POP_ATOM(I)+ED(I))*T(I)
	END DO
	CALL MON_INT_FUNS_V2(COEF,MOD_PRESSURE,R,ND)
!	
	DO I=1,ND
	  VdVdR=V(I)*V(I)*(SIGMA(I)+1.0_LDP)/R(I)
	  dPdR_ON_ROH=dP_CON*COEF(I,3)/POP_ATOM(I)
!
! Gravitational force per unit mass (in cgs units).
!
	  g_grav=GRAV_CON/R(I)/R(I)
!
! Radiation pressure forces
!
	  g_rad=RAD_CON*LUM_STAR(I)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	  g_elec=RAD_CON*LUM_STAR(I)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
!
	  g_TOT= g_RAD-g_GRAV
	  ERROR=200.0_LDP*(VdVdR+dPdR_ON_ROH-g_TOT)/
	1        ( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(g_TOT) )
!
	  IF(R(I) .GT. 9.99E+04_LDP)THEN
	    FMT='(1X,ES12.6,ES13.4,F9.2,5(ES14.4),F9.2,I7)'
	  ELSE
	    FMT='(1X,F12.6,ES13.4,F9.2,5(ES14.4),F9.2,I7)'
	  END IF
	  WRITE(LU_OUT,FMT)
	1             R(I),V(I),ERROR,VdVdR,dPdR_ON_ROH,
	1             g_TOT,g_RAD,g_ELEC,g_RAD/g_GRAV,I
!
	END DO
!
! Write out a summary of momentum euqation etc to remind the user of the
! meaning of individual terms.
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Momentum equation is:',
	1                    ' VdV/dr = - dPdR/ROH - g + g_RAD'
	WRITE(LU_OUT,'(1X,A,A)')'        or          :',
	1                    ' VdV/dr = - dPdR/ROH + g_tot'
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Error is 200.0D0*(VdVdR+dPdR_ON_ROH-g_TOT)/',
	1        '( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(g_TOT) )'
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Gamma = g_rad/g [g=g_GRAV] '
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,1PE14.4)')
	1          'Surface gravity is: ',GRAV_CON/R(ND)/R(ND)
	WRITE(LU_OUT,'(1X,A,F8.2,A)')
	1          'Stars mass is: ',STARS_MASS,' Msun'
!
	RLUMST(:)=BAND_FLUX(:,N_FLUXMEAN_BANDS)
        DO J=N_FLUXMEAN_BANDS,2,-1
          BAND_FLUXMEAN(:,J)=BAND_FLUXMEAN(:,J)-BAND_FLUXMEAN(:,J-1)
          BAND_FLUX(:,J)=BAND_FLUX(:,J)-BAND_FLUX(:,J-1)
        END DO
        DO J=1,N_FLUXMEAN_BANDS
          BAND_FLUXMEAN(:,J)=BAND_FLUXMEAN(:,J)/FLUX_MEAN_OPAC(:)
          BAND_FLUX(:,J)=BAND_FLUX(:,J)/RLUMST(:)
        END DO
!
	WRITE(LU_OUT,'(A)')CHAR(12)                         !Formfeed
        WRITE(LU_OUT,'(3X,A,5X,A,6X,14F12.2))')'I','V',
	1                      (LAM_FLUXMEAN_BAND_END(J),J=1,N_FLUXMEAN_BANDS)
        DO I=1,ND
          WRITE(LU_OUT,'(X,I3,15ES12.3))')I,V(I),(BAND_FLUXMEAN(I,J),J=1,N_FLUXMEAN_BANDS)
          WRITE(LU_OUT,'(16X,15ES12.3))')(BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
	  T1=FLUX_MEAN_OPAC(I)/ELEC_MEAN_OPAC(I)
          WRITE(LU_OUT,'(16X,15ES12.3))')(T1*BAND_FLUXMEAN(I,J)/BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
        END DO
!
	CLOSE(LU_OUT)
	RETURN
	END
