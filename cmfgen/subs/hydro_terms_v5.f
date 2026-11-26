!
! Subroutine to compare terms in the Hydronamical equations. Routine is
! meant for illustration purposes only.
!
	SUBROUTINE HYDRO_TERMS_V5(POP_ATOM,R,V,T,SIGMA,ED,CLUMP_FAC,LUM_STAR,
	1              LOGG,STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1              FLUX_MEAN_OPAC,ROSS_MEAN_OPAC,ELEC_MEAN_OPAC,
	1              LAM_FLUXMEAN_BAND_END,BAND_FLUXMEAN,
	1              PRESSURE_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NOV,
	1              LST_ITERATION,BAND_FLUX,N_FLUXMEAN_BANDS,LU_OUT,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 21-Jun-2023 : Added log g consistency check.
! Altered 24-Aug-2022 : Added descrition for flux and flux mean opacity output.
! Altered 11-Jan-2021 : Increased minimum size of record in HYDRO_TERMS
! Altered 20-Oct-2020 : Output new error (100 E/g_rad) in final column.
! Altered 01-Jan-2020 : Increased R precision. Increased Tau precision.
! Altered 13-Feb-2019 : Tau (Rosseland) now output.
! Altered 08-Nov-2016 : Now output Gamma (e.s.) at photosphere.
! Altered 01-Jan-2015 : Turbulent pressure term added, and also output when non-zero.
! Altered 26-Jan-2014 : Changed ouput to TYPE_ATM and increased its length.
! Altered 06-Jan-2014 : Changed to V5 (inserted LOGG, ROSS_MEAN_OPAC and CLUMP_FAC in call).
! Altered 31-Mar-2010 : Changed from V3 to V4: Inserted LST_ITERATION in call.
! Altered 11-Mar-2008 : Changed from V2 to V3
!                       PRESSURE_VTURB,PLANE_PARALLEL & PLANE_PARALLEL_NOV inserted into call.
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
	REAL(KIND=LDP) CLUMP_FAC(ND)
	REAL(KIND=LDP) LUM_STAR(ND)
!
	REAL(KIND=LDP) FLUX_MEAN_OPAC(ND)
	REAL(KIND=LDP) ROSS_MEAN_OPAC(ND)
	REAL(KIND=LDP) ELEC_MEAN_OPAC(ND)
!
	INTEGER N_FLUXMEAN_BANDS
	REAL(KIND=LDP) LAM_FLUXMEAN_BAND_END(ND)
	REAL(KIND=LDP) BAND_FLUXMEAN(ND,N_FLUXMEAN_BANDS)
	REAL(KIND=LDP) BAND_FLUX(ND,N_FLUXMEAN_BANDS)
!
	REAL(KIND=LDP) LOGG
	REAL(KIND=LDP) STARS_MASS
	REAL(KIND=LDP) MEAN_ATOMIC_WEIGHT
	REAL(KIND=LDP) PRESSURE_VTURB
	LOGICAL PLANE_PARALLEL
	LOGICAL PLANE_PARALLEL_NOV
	LOGICAL LST_ITERATION
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
	REAL(KIND=LDP) TAU(ND)
	REAL(KIND=LDP) OPAC(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) TC(ND)
!
! Local variables.
!
	REAL(KIND=LDP) VdVdR
	REAL(KIND=LDP) dPdR_ON_ROH
	REAL(KIND=LDP) dPTURBdR_ON_ROH
	REAL(KIND=LDP) g_TOT
	REAL(KIND=LDP) g_ELEC
	REAL(KIND=LDP) g_GRAV
	REAL(KIND=LDP) g_RAD
	REAL(KIND=LDP) ERROR,REV_ERROR
	REAL(KIND=LDP) RLUMST(ND)
!
	REAL(KIND=LDP) GRAV_CON
	REAL(KIND=LDP) dP_CON
	REAL(KIND=LDP) PTURB_CON
	REAL(KIND=LDP) RAD_CON
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) SOUND_SPEED
	REAL(KIND=LDP) RPHOT
	REAL(KIND=LDP) GAM_ES_PHOT
	REAL(KIND=LDP) GPHOT
	REAL(KIND=LDP) ERROR_MAX
	REAL(KIND=LDP) ERROR_SUM
	REAL(KIND=LDP) ERROR_SQ
	REAL(KIND=LDP) PHOT_LOGG
!
	CHARACTER(LEN=12) TYPE_ATM
	CHARACTER(LEN=100) FMT
	CHARACTER(LEN=20) TMP_FMT
	INTEGER I,J,IOS,ERROR_LU
	INTEGER ERROR_CNT
	EXTERNAL ERROR_LU
	INTEGER, PARAMETER :: IONE=1
!
	I=MAX(200,(N_FLUXMEAN_BANDS+1)*12+16)
	CALL GEN_ASCI_OPEN(LU_OUT,'HYDRO','UNKNOWN',' ',' ',I,IOS)
	IF(IOS .NE. 0)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error opening HYDRO in HYDRO_TERMS'
	  WRITE(I,*)'IOS=',IOS
	END IF
!
	ERROR_SUM=0.0_LDP
	ERROR_SQ=0.0_LDP
	ERROR_MAX=0.0_LDP
	ERROR_CNT=0
!
! 1.0D-06 = [ 10^(-10) {from 1/r) * 10^4 {from T}
!
	dP_CON=1.0E-06_LDP*BOLTZMANN_CONSTANT()/MEAN_ATOMIC_WEIGHT/
	1               ATOMIC_MASS_UNIT()
	PTURB_CON=1.0E+10_LDP*0.5_LDP*PRESSURE_VTURB*PRESSURE_VTURB
!
	GRAV_CON=1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()
!
!       1.0E-10 (CHI) / 1.0E-20 (1/r^2)
!
	RAD_CON=1.0E-30_LDP*LUM_SUN()/ATOMIC_MASS_UNIT()/4.0_LDP/FUN_PI()/
	1                  SPEED_OF_LIGHT()/MEAN_ATOMIC_WEIGHT
!
	TYPE_ATM='P'
        T1=LOG(POP_ATOM(5)/POP_ATOM(1))/LOG(R(1)/R(5))
        WRITE(TYPE_ATM(2:),'(ES10.3)')T1
!
! Compute photospheric radius and photospheric surface gravity.
!
	DO I=1,ND
	  OPAC(I)=MAX(ROSS_MEAN_OPAC(I),ELEC_MEAN_OPAC(I))*CLUMP_FAC(I)
	END DO
	CALL TORSCL(TAU,OPAC,R,TB,TC,ND,'LOGMON',TYPE_ATM)
	T1=0.6667_LDP
	CALL MON_INTERP(RPHOT,IONE,IONE,T1,IONE,R,ND,TAU,ND)
	IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NOV)RPHOT=R(ND)
	GPHOT=GRAV_CON/RPHOT/RPHOT
!
! Get Gamma due to e.s. at photosphere. TB is used as a temporary work vector
! for G_ELEC.
!
	IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NOV)THEN
	  DO I=1,ND
	    TB(I)=RAD_CON*LUM_STAR(ND)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(ND)/R(ND)
	  END DO
	ELSE
	  DO I=1,ND
	    TB(I)=RAD_CON*LUM_STAR(ND)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	  END DO
	END IF
	T1=0.6667_LDP
	CALL MON_INTERP(GAM_ES_PHOT,IONE,IONE,T1,IONE,TB,ND,TAU,ND)
	GAM_ES_PHOT=GAM_ES_PHOT/GPHOT
!
! Output file header.
!
	IF(PRESSURE_VTURB .EQ. 0.0_LDP)THEN
	  WRITE(LU_OUT,'(1X,6X,A,6X, 12X,A,3X, 2X,A, 5(4X,A,1X),4X,A,2X,A,9X,A,3X,A9,3X,A)')
	1     'R','V','% Error','    VdVdR',
	1                       ' dPdR/ROH',
	1                       '    g_TOT',
	1                       '    g_RAD',
	1                       '   g_ELEC','Gamma','Depth','Tau','Vsound','E/g_rad%'
	ELSE
	  WRITE(LU_OUT,'(1X,6X,A,6X, 12X,A,3X, 2X,A, 6(4X,A,1X), 4X,A,2X,A,9X,A,3X,A,3X,A)')
	1     'R','V','% Error','    VdVdR',
	1                       ' dPdR/ROH',
	1                       'dTPdR/ROH',
	1                       '    g_TOT',
	1                       '    g_RAD',
	1                       '   g_ELEC','Gamma','Depth','Tau','Vsound','E/g_rad)%'
	END IF
!
	DO I=1,ND
	  MOD_PRESSURE(I)=(POP_ATOM(I)+ED(I))*T(I)
	END DO
	CALL MON_INT_FUNS_V2(COEF,MOD_PRESSURE,R,ND)
!
! Note: g_elec (for I=ND) is used late in the code.
!
	DO I=1,ND
	  VdVdR=V(I)*V(I)*(SIGMA(I)+1.0_LDP)/R(I)
	  dPdR_ON_ROH=dP_CON*COEF(I,3)/POP_ATOM(I)
	  dPTURBdR_ON_ROH=PTURB_CON*1.0E-10_LDP*(-2.0_LDP/R(I)-(SIGMA(I)+1.0_LDP)/R(I))
!
! Gravitational force per unit mass (in cgs units) and radiation pressure forces.
!
	  IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NOV)THEN
	    g_grav=GRAV_CON/R(ND)/R(ND)
	    g_rad=RAD_CON*LUM_STAR(ND)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(ND)/R(ND)
	    g_elec=RAD_CON*LUM_STAR(ND)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(ND)/R(ND)
	  ELSE
	    g_grav=GRAV_CON/R(I)/R(I)
	    g_rad=RAD_CON*LUM_STAR(ND)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	    g_elec=RAD_CON*LUM_STAR(ND)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
!	    g_rad=RAD_CON*LUM_STAR(I)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
!	    g_elec=RAD_CON*LUM_STAR(I)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	  END IF
!
	  g_TOT= g_RAD-g_GRAV
	  ERROR=200.0_LDP*(VdVdR+dPdR_ON_ROH+dPTURBdR_ON_ROH-g_TOT)/
	1        ( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(dPTURBdR_ON_ROH)+ ABS(g_TOT) )
	  REV_ERROR=100.0_LDP*(VdVdR+dPdR_ON_ROH+dPTURBdR_ON_ROH-g_TOT)/ABS(g_RAD)
!
	  IF(V(I) .LT. 5.0_LDP)THEN
	    ERROR_SUM=ERROR_SUM+ERROR
	    ERROR_SQ=ERROR_SQ+ERROR*ERROR
	    ERROR_MAX=MAX(ERROR_MAX,ABS(ERROR))
	    ERROR_CNT=ERROR_CNT+1
	  END IF
!
! Compute the sound speed in km/s. Note: The constant is 1.0D-06 because T is in unts of 10^4 K,
! and a factor of 10^{-10} is used to convert to km/s.
!
	  T1=1.0E-06_LDP*BOLTZMANN_CONSTANT()/MEAN_ATOMIC_WEIGHT/ATOMIC_MASS_UNIT()
	  T2=PRESSURE_VTURB
	  SOUND_SPEED=SQRT( T1*T(I)*(1.0_LDP+ED(I)/POP_ATOM(I)) + 0.5_LDP*T2*T2)
!
	  TMP_FMT='F10.2,F10.2,F12.2)'
	  IF(TAU(I) .LT. 1.0_LDP)TMP_FMT='F10.3,F10.2,F12.2)'
	  IF(TAU(I) .LT. 0.1_LDP)TMP_FMT='ES10.2,F10.2,F12.2)'
	  IF(PRESSURE_VTURB .EQ. 0.0_LDP)THEN
	    IF(R(I) .GT. 9.99E+04_LDP)THEN
	      FMT='(1X,ES12.6,ES17.8,F9.2,5(ES14.4),F9.2,I7,2X,'//TMP_FMT
	    ELSE
	      FMT='(1X,F12.6,ES17.8,F9.2,5(ES14.4),F9.2,I7,2X,'//TMP_FMT
	    END IF
	    WRITE(LU_OUT,FMT)
	1             R(I),V(I),ERROR,VdVdR,dPdR_ON_ROH,
	1             g_TOT,g_RAD,g_ELEC,g_RAD/g_GRAV,I,TAU(I),SOUND_SPEED,REV_ERROR
	  ELSE
	    IF(R(I) .GT. 9.99E+04_LDP)THEN
	      FMT='(1X,ES15.9,ES17.8,F9.2,6(ES14.4),F9.2,I7,2X,'//TMP_FMT
	    ELSE
	      FMT='(1X,F15.9,ES17.8,F9.2,6(ES14.4),F9.2,I7,2X,'//TMP_FMT
	    END IF
	    WRITE(LU_OUT,FMT)
	1             R(I),V(I),ERROR,VdVdR,dPdR_ON_ROH,dPTURBdR_ON_ROH,
	1             g_TOT,g_RAD,g_ELEC,g_RAD/g_GRAV,I,TAU(I),SOUND_SPEED,REV_ERROR
	  END IF
!
	END DO
!
! Write out a summary of momentum euqation etc to remind the user of the
! meaning of individual terms.
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Momentum equation is:',
	1                    ' VdV/dr = - dPdR/ROH - dTPdR/ROH - g + g_RAD'
	WRITE(LU_OUT,'(1X,A,A)')'        or          :',
	1                    ' VdV/dr = - dPdR/ROH - dTPdR/ROH + g_tot'
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Error is 200.0D0*(VdVdR+dPdR_ON_ROH+dTPdR/ROH-g_TOT)/',
	1        '( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(dTPdR/ROH)+ ABS(g_TOT) )'
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Gamma = g_rad/g [g=g_GRAV] '
!
	WRITE(LU_OUT,'(1X,A)')' '
	PHOT_LOGG=GRAV_CON/RPHOT/RPHOT
	WRITE(LU_OUT,'(1X,A,ES14.4,A,4X,ES10.4,A)')
	1          '         Photospheric radius is: ',RPHOT,'(10^10 cm)',RPHOT/6.96D0,'(Rsun)'
	WRITE(LU_OUT,'(1X,A,ES14.4,A,F7.4,A)')
	1          'Photospheric surface gravity is: ',PHOT_LOGG,' (',LOG10(PHOT_LOGG),')'
	WRITE(LU_OUT,'(1X,A,ES14.4,A,F7.4,A)')
	1          '   Specified surface gravity is: ',10**LOGG,' (',LOGG,')'
	WRITE(LU_OUT,'(1X,A,F10.4,A)')
	1          '                  Stars mass is: ',STARS_MASS,' Msun'
	WRITE(LU_OUT,'(1X,A,F8.5,A)')
	1          'Edd. fac. (e.s.) at photosphere: ',GAM_ES_PHOT
	WRITE(LU_OUT,'(1X,A,F8.5,A)')
	1          'Edd. fac. (e.s.) at inner bound: ',g_elec*R(ND)*R(ND)/GRAV_CON
	WRITE(LU_OUT,'(A,A)')' Power law exponent at outer boundary (TYPE_ATM)=',TRIM(TYPE_ATM(2:))
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
	WRITE(LU_OUT,'(/,A)')'Row 1: Fraction of flux mean opacity contained in the band.'
	WRITE(LU_OUT,  '(A)')'Row 2: Fraction of flux contained in the band.'
	WRITE(LU_OUT,'(A,/)')'Row 3: Flux mean opacity normalized by the electron scatering opacity in the band.'
!
        WRITE(LU_OUT,'(3X,A,5X,A,6X,13F12.2,ES12.2)')'I','V',
	1                      (LAM_FLUXMEAN_BAND_END(J),J=1,N_FLUXMEAN_BANDS)
        DO I=1,ND
          WRITE(LU_OUT,'(X,I3,15ES12.3)')I,V(I),(BAND_FLUXMEAN(I,J),J=1,N_FLUXMEAN_BANDS)
          WRITE(LU_OUT,'(16X,15ES12.3)')(BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
	  T1=FLUX_MEAN_OPAC(I)/ELEC_MEAN_OPAC(I)
	  DO J=1,N_FLUXMEAN_BANDS
	    IF(ABS(BAND_FLUX(I,J)) .LT. 1.0E-80_LDP)BAND_FLUX(I,J)=1.0E-80_LDP
          END DO
	  WRITE(LU_OUT,'(16X,15ES12.3)')(T1*BAND_FLUXMEAN(I,J)/BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
        END DO
!
	CLOSE(LU_OUT)
!
	IF(LST_ITERATION)THEN
!
	  IF( ABS(LOG10(PHOT_LOGG)-LOGG) .GT. 0.02_LDP )THEN
	    I=ERROR_LU()
	    WRITE(I,*)'*****************************************************************************************'
	    WRITE(I,*)'Warning inconsistent surface gravities'
	    WRITE(I,*)'This error has occured because you forgot to set HYDRO_DEFAULTS or'
	    WRITE(I,*)'       or you have done insufficient interations'
	    WRITE(I,*)'You can ignore this error if not doing Hydro iterations.'
	    WRITE(I,*)'*****************************************************************************************'
	  END IF

	  ERROR_SQ=SQRT(ERROR_SQ/MAX(1,ERROR_CNT))
	  ERROR_SUM=ERROR_SUM/MAX(1,ERROR_CNT)
	  T1=ABS(LOGG-LOG10(GPHOT))
	  IF(ERROR_SQ .GT. 5.0_LDP .OR. ERROR_MAX .GT. 20.0_LDP .OR. T1 .GT. 0.02_LDP)THEN
	    I=ERROR_LU()
	    WRITE(I,*)' '
	    WRITE(I,*)'*****************************************************************************************'
	    WRITE(I,'(A)')' Possible error with hydrostatic structure --- large error in photosphere'
	    WRITE(I,'(A,ES10.2)')'              Mean error is',ERROR_SUM
	    WRITE(I,'(A,ES10.2)')' Root mean squared error is',ERROR_SQ
	    WRITE(I,'(A,ES10.2)')'           Maximum error is',ERROR_MAX
	    WRITE(I,'(A,ES10.2)')'       Stars mass (Msun) is',STARS_MASS
	    WRITE(I,'(A,ES10.2)')'        R(phot, 10^10cm) is',RPHOT
	    WRITE(I,'(A,ES10.2)')'             log G(phot) is',LOG10(GPHOT)
	    WRITE(I,'(A,ES10.2)')'         Specified log g is',LOGG
	    WRITE(I,*)'*****************************************************************************************'
	    WRITE(I,*)' '
	  END IF
	END IF
!
	RETURN
	END
