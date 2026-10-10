!
! Subroutine to compute CLUMP_FAC, and the abundance vectors for
! each species. R and V must be defined. This subroutine can be
! called again if the R grid is altered.
!
! Notes:
!     AT_ABUND can be a number abundance or a mass fraction (if -ve).
!              On exit, it is always a number fraction.
!     MEAN_ATOMIC_WEIGHT is the mean atomic/ionic weight in amu. Electrons
!              are not included.
!     ABUND_SUM is the total number of atoms present on the scale where some
!         species (usually H or He) is one.
!
	SUBROUTINE SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
	USE SET_KIND_MODULE
	USE CONTROL_VARIABLE_MOD
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 05-Oct-2026 : Use ABS(V)/MAXVAL(ABS(V)) for non-mono VTURB (Devin).
! Altered 29-Jul-2023 : Update VTURB_VEC -- now allow for POW law option (LONG ver -- 15-Oct-2023).
! Altered 16-Apr-2023 : Added RUN2 and REX2 from Paco (Paul now uses RUNA which is
!	                      the same option).
! Altered 03-MAr-2023 : Added clumping law from Paul Crowther
! Altered 15-May-2012 : Modications to clumping laws (imported from OSIRIS: 21-Jul-2021).
! Altered 09-May-2021 : Added RPOW option.
! Altered 14-Dec-2020 : Altered EXPO option to inhibit clumping below a certain to velocity.
! Altered 01-Oct-2018 : Added MEXP option.
! Altered 14-Sep-2018 : Added POW option.
! Altered 21-Mar-2018 : Added SNCL clumping option.
! Altered 25-Jun-2010 : Correctly compute MEAN_ATOMIC_WEIGHT and ABUND_SUM for
!                         SN model. Abunances are surface abundances.
! Altered 08-Feb-2006 : New clumping law installed (used by PACO for P Cygni).
! Altered 05-Sep-2005 : VAR_MDOT section included.
! Created 19-Dec-2004
!
	INTEGER ND
	INTEGER LUER
	REAL(KIND=LDP) MEAN_ATOMIC_WEIGHT	!In amu (does not include electrons)
	REAL(KIND=LDP) ABUND_SUM
!
! Local variables.
!
	REAL(KIND=LDP) T1,T2,T3,T4,PI
!
	INTEGER I
	INTEGER ISPEC
	INTEGER K
!
	INTEGER LUIN,NHYDRO
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP), ALLOCATABLE :: R_HYDRO(:),RHO_HYDRO(:),T_HYDRO(:),LOG_R_HYDRO(:)
	CHARACTER*80 STRING
	LOGICAL ADD_ACC_ZONE
!
        REAL(KIND=LDP) ATOMIC_MASS_UNIT
	EXTERNAL ATOMIC_MASS_UNIT
!
!
! Allow for the possibility that the material in the wind is clumped. The
! assumption that goes into our treatment are not rigorous --- rather this
! treatment should be used as a gauge to discern the errors arising by treating
! O and W-R winds as homogeneous.
!
	IF(DO_CLUMP_MODEL)THEN
	  IF(GLOBAL_LINE_SWITCH .NE. 'BLANK' .AND.
	1                        .NOT. FLUX_CAL_ONLY)THEN
	    WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	    WRITE(LUER,*)'Clumping is only treated for a fully blanketed model'
	    STOP
	  ELSE IF(GLOBAL_LINE_SWITCH .NE. 'BLANK')THEN
	    WRITE(LUER,*)'Warning in SET_ABUND_CLUMP'
	    WRITE(LUER,*)'Clumping is only treated for a fully blanketed model'
	    WRITE(LUER,*)'EW in EWDATA file are incorrect'
	    WRITE(LUER,*)'Continuum fluxes will be okay'
	  END IF
	  IF(CLUMP_LAW(1:4) .EQ. 'EXPO')THEN
!
! CLUMP_PAR(1) is the clumping factor at infinity.
! CLUMP_PAR(2) is a velocity, and determines how fast the clumping factor
! approach CLUMP_PAR(1).
!
	    IF(CLUMP_PAR(3) .EQ. 0.0_LDP)CLUMP_PAR(4)=1.0_LDP
	    T2=CLUMP_PAR(6)
	    DO K=1,ND
	      T1=V(K)
	      IF(N_CLUMP_PAR .EQ. 6)T1=LOG(1.0_LDP+EXP(T2*(V(K)-CLUMP_PAR(5))))/T2
	      CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1)-CLUMP_PAR(3))*
	1                     EXP(-T1/CLUMP_PAR(2))+
	1                     CLUMP_PAR(3)*EXP(-T1/CLUMP_PAR(4))
	    END DO
!
	  ELSE IF(CLUMP_LAW(1:4) .EQ. 'MEXP')THEN
	    IF(N_CLUMP_PAR .NE. 5)THEN
	      WRITE(LUER,*)'Error in SET_ABUND_CLUMP for MEXP N_CLUMP_PAR should be 5'
	      WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	      STOP
	    END IF
	    DO K=1,ND
	      T1=1.0_LDP-(MAX(1.0_LDP,CLUMP_PAR(4)*R(1)/R(K)))**CLUMP_PAR(5)
	      CLUMP_FAC(K)=(CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1))*EXP(-V(K)/CLUMP_PAR(2)))/
	1                     (1.0_LDP+CLUMP_PAR(3)*EXP(T1))
	    END DO

	  ELSE IF(CLUMP_LAW(1:4) .EQ. 'DBLE')THEN
!
! CLUMP_PAR(1) is the clumping factor at infinity.
! CLUMP_PAR(2)
!
	    DO K=1,ND
	      T1=10.0_LDP/CLUMP_PAR(4)
	      T2=10.0_LDP/CLUMP_PAR(3)
	      CLUMP_FAC(K)= CLUMP_PAR(1)
	1             + (CLUMP_PAR(2)-CLUMP_PAR(1))/(1.0_LDP+EXP( T1*(V(K)-CLUMP_PAR(4)) ) )
	1             + (1.0_LDP-CLUMP_PAR(2)-CLUMP_PAR(1))/(1.0_LDP+EXP( T2*(V(K)-CLUMP_PAR(3)) ) )
	    END DO
!
	  ELSE IF(CLUMP_LAW(1:4) .EQ. 'REXP')THEN
	    IF(N_CLUMP_PAR .NE. 3)THEN
	      WRITE(LUER,*)'Error in SET_ABUND_CLUMP - for MEXP N_CLUMP_PAR should be 3'
	      WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	      STOP
	    END IF
	    DO K=1,ND
	      CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1))*
	1          EXP(-V(K)/CLUMP_PAR(2))+
	1	   (1.0_LDP-CLUMP_PAR(1))*EXP( (V(K)-V(1))/CLUMP_PAR(3))
	    END DO
!
	  ELSE IF(CLUMP_LAW(1:3) .EQ. 'POW')THEN
	    IF(N_CLUMP_PAR .NE. 2)THEN
	      WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	      WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	      STOP
	    END IF
	    DO K=1,ND
	     CLUMP_FAC(K)=1.0_LDP-(1.0_LDP-CLUMP_PAR(1))*(V(K)/V(1))**CLUMP_PAR(2)
	    END DO
!
	  ELSE IF(CLUMP_LAW(1:4) .EQ. 'RPOW')THEN
	    IF(N_CLUMP_PAR .NE. 2 .AND. N_CLUMP_PAR .NE. 4)THEN
	      WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	      WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	      STOP
	    END IF
	    IF(N_CLUMP_PAR .EQ. 2)THEN
	      CLUMP_PAR(3)=0.0_LDP; CLUMP_PAR(4)=1.0_LDP
	    END IF
	    DO K=1,ND
	      T1=1.0_LDP/CLUMP_PAR(1)-1.0_LDP
	      CLUMP_FAC(K)=1.0_LDP/(1.0_LDP+T1*(V(K)/V(1))**CLUMP_PAR(2))/
	1               (1.0_LDP+CLUMP_PAR(3)*(V(K)/V(1))**CLUMP_PAR(4))
	    END DO
!
	  ELSE IF(CLUMP_LAW(1:4) .EQ. 'SNCL')THEN
	    IF (CLUMP_PAR(3) .EQ. 0.0_LDP) THEN
	       DO K=1,ND
		  T1 = (V(K) - V(ND)) / CLUMP_PAR(2)
		  CLUMP_FAC(K)=1.0_LDP + (CLUMP_PAR(1)-1.0_LDP)*EXP(-T1*T1)
	       END DO
	    ELSE
	       DO K=1,ND
		  IF (V(K).LE.CLUMP_PAR(3)) THEN
		     CLUMP_FAC(K) = CLUMP_PAR(1)
		  ELSE
		     T1 = (V(K) - CLUMP_PAR(3)) / CLUMP_PAR(2)
		     CLUMP_FAC(K)=1.0_LDP + (CLUMP_PAR(1)-1.0_LDP)*EXP(-T1*T1)
		  ENDIF
	       END DO
	    ENDIF
!
! For consistency with bith Paul and Paco nomenclature.
!
	  ELSE IF(CLUMP_LAW(1:4) .EQ. 'RUNA' .OR. CLUMP_LAW(1:4) .EQ. 'PEXP')THEN
	    IF(N_CLUMP_PAR .NE. 4)THEN
	      WRITE(LUER,*)'Error in CMFGEN'
	      WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	      STOP
	    END IF
	    DO K=1,ND
	      CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1))*
	1     EXP(-V(K)/CLUMP_PAR(2))+
	1     (CLUMP_PAR(4)-CLUMP_PAR(1))*
	1     EXP( (V(K)-V(1))/CLUMP_PAR(3))
	      IF(CLUMP_FAC(K).GT. 1._LDP) CLUMP_FAC(K)=1.0_LDP
	    END DO
!
          ELSE IF(CLUMP_LAW(1:4) .EQ. 'REX2')THEN
	    IF(N_CLUMP_PAR .NE. 3)THEN
	      WRITE(LUER,*)'Error in CMFGEN'
	      WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	      STOP
	    END IF
	    DO K=1,ND
	      CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1))*
	1          EXP(-(V(K)/CLUMP_PAR(2))**2._LDP)+
	1	   (1.0_LDP-CLUMP_PAR(1))*EXP(-((V(K)-V(1))/CLUMP_PAR(3))**2.)
              IF(CLUMP_FAC(K).GT. 1._LDP) CLUMP_FAC(K)=1.0_LDP
	    END DO
!
          ELSE IF(CLUMP_LAW(1:4) .EQ. 'RUN2')THEN
            IF(N_CLUMP_PAR .NE. 4)THEN
	     WRITE(LUER,*)'Error in CMFGEN'
             WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	     STOP
	    END IF
	    DO K=1,ND
              CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1))*
     &          EXP(-(V(K)/CLUMP_PAR(2))**2._LDP)+
     &          (CLUMP_PAR(4)-CLUMP_PAR(1))*
     &          EXP(-((V(K)-V(1))/CLUMP_PAR(3))**2.)
              IF(CLUMP_FAC(K).GT. 1._LDP) CLUMP_FAC(K)=1.0_LDP
	    END DO
!
	  ELSE IF(CLUMP_LAW(1:3) .EQ. 'SIN')THEN
	    IF(N_CLUMP_PAR .NE. 2)THEN
	      WRITE(LUER,*)'Error in CMFGEN'
	      WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',N_CLUMP_PAR
	      STOP
	    END IF
	    PI=ACOS(-1.0_LDP)
	    DO K=1,ND
	      T1=MOD(CLUMP_PAR(1)*LOG(R(K)/R(ND)),LOG(R(1)/R(ND)))/LOG(R(1)/R(ND))
	      T2=V(K)*(SIN(2*PI*T1)-CLUMP_PAR(2))/V(1)
	      CLUMP_FAC(K)=10**T2
	      IF(CLUMP_FAC(K).GT. 1._LDP) CLUMP_FAC(K)=1.0_LDP
	    END DO
	  ELSE IF(CLUMP_LAW(1:6) .EQ. 'SPLINE')THEN
	    CALL SPL_CLUMP(CLUMP_FAC,V,ND)
!
	 ELSE
	    WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	    WRITE(LUER,*)'Clump law=',TRIM(CLUMP_LAW)
	    WRITE(LUER,*)'Invalid law for computing clumping factor'
	    STOP
	 END IF
	ELSE
	  DO K=1,ND
	    CLUMP_FAC(K)=1.0_LDP
	  END DO
	END IF
!
	IF(SN_HYDRO_MODEL .AND. VELTYPE .NE. 12)THEN
	  ABUND_SUM=0.0_LDP
	  CALL RD_SN_DATA(ND,NEWMOD,7)
	  DO ISPEC=1,NUM_SPECIES
	    IF(AT_ABUND(ISPEC) .GT. 0)THEN
	      ABUND_SUM=ABUND_SUM+AT_ABUND(ISPEC)
              MEAN_ATOMIC_WEIGHT=MEAN_ATOMIC_WEIGHT+AT_MASS(ISPEC)*AT_ABUND(ISPEC)
	    END IF
	  END DO
	  MEAN_ATOMIC_WEIGHT=MEAN_ATOMIC_WEIGHT/ABUND_SUM
	  CALL UPDATE_VTURB()
	  RETURN
	END IF
! 
! Compute the atomic number density as function of radius. The abundances
! can be specified in two ways:
!     (1) As a fractional abundance relative to some arbitrary species X.
!     (2) As a mass fraction: Assumed when the abundance is negative.
!
!      N#= SUM[ Xi/X ]            ; Xi/X > 0
!      mu#= SUM[ Ai Xi/X ] / N#   ; Xi/X > 0
!      F#= SUM[Fi]=SUM[ Xi/X ]    ; Xi/X < 0
!      1/mu = f#/mu# + SUM[Fi/Ai]
!
	T1=0.0_LDP		!Sum of abundance fractions (-> N#)
	T2=0.0_LDP		!Sum of A_i * abundance fractions
	T3=0.0_LDP		!Sum of mass fractions      (-> F#)
	T4=0.0_LDP                !Sum of f/A
!
	DO ISPEC=1,NUM_SPECIES
	  IF(AT_ABUND(ISPEC) .GT. 0)THEN
	    T1=T1+AT_ABUND(ISPEC)
            T2=T2+AT_MASS(ISPEC)*AT_ABUND(ISPEC)
	  ELSE
	    T3=T3+ABS(AT_ABUND(ISPEC))
	    T4=T4+ABS(AT_ABUND(ISPEC))/AT_MASS(ISPEC)
	  END IF
	END DO
	IF(T1 .NE. 0.0_LDP .AND. T3 .GT. 1.0_LDP)THEN
	  WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	  WRITE(LUER,*)'Sum of mass fractions should be less than unity'
	  WRITE(LUER,*)'Mass fraction=',T3
	  STOP
	END IF
	IF(T1 .EQ. 0.0_LDP .AND. ABS(T3-1.0_LDP) .GT. 0.001_LDP)THEN
	  WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	  WRITE(LUER,*)'Sum of mass fractions should be unity'
	  WRITE(LUER,*)'Mass fraction=',T3
	  STOP
	END IF
!
!NB: 1-T3 is the mass fraction of species whose abundances were specified in
!          terms of fractional number densities.
!
	IF(T2 .NE. 0.0_LDP)T2=(1.0_LDP-T3)*T1/T2		!f#/mu#
	MEAN_ATOMIC_WEIGHT=1.0_LDP/(T2+T4)
!
! Convert any mass fractions to number fractions. If no number fractions
! are present, the first species with non-zero abundance is arbitrarily set
! to have an abundance of unity.
!
	ABUND_SUM=0.0_LDP
	IF(T1 .EQ. 0.0_LDP)T1=1.0_LDP
	DO ISPEC=1,NUM_SPECIES
	  IF(T2 .EQ. 0.0_LDP)T2=ABS(AT_ABUND(ISPEC))/AT_MASS(ISPEC)
	  IF(AT_ABUND(ISPEC) .LT. 0.0_LDP)AT_ABUND(ISPEC)=
	1         T1*ABS(AT_ABUND(ISPEC))/T2/AT_MASS(ISPEC)
	  ABUND_SUM=ABUND_SUM+AT_ABUND(ISPEC)
	END DO
!
	IF (VELTYPE .EQ. 12) THEN
          LUIN = 7
          OPEN(UNIT=LUIN,STATUS='OLD',FILE='input_hydro.dat',IOSTAT=IOS)
          READ(LUIN,*) ADD_ACC_ZONE
          IF (ADD_ACC_ZONE)READ(LUIN,*)(T1, I=1,5)
          READ(LUIN,*) STRING
          READ(LUIN,*) NHYDRO
          READ(LUIN,*) STRING  ! reads the bogus line from python script
          ALLOCATE (R_HYDRO(NHYDRO),T_HYDRO(NHYDRO),RHO_HYDRO(NHYDRO),LOG_R_HYDRO(NHYDRO))
          DO I=1,NHYDRO
             READ(LUIN,*) R_HYDRO(NHYDRO-I+1),T1,RHO_HYDRO(NHYDRO-I+1),T_HYDRO(NHYDRO-I+1)
          ENDDO
          CLOSE(LUIN)
!
          LOG_R=LOG(R)
          LOG_R_HYDRO = LOG(R_HYDRO)
          RHO_HYDRO = LOG(RHO_HYDRO)
          CALL MON_INTERP(DENSITY,ND,IONE,LOG_R,ND,RHO_HYDRO,NHYDRO,LOG_R_HYDRO,NHYDRO)
          DENSITY=EXP(DENSITY)
!
          DO K=1,ND
            POP_ATOM(K)=DENSITY(K)/MEAN_ATOMIC_WEIGHT/ATOMIC_MASS_UNIT()
            POP_ATOM(K)=POP_ATOM(K)/CLUMP_FAC(K)
          END DO
!
	ELSE IF(SN_MODEL)THEN
	  DO K=1,ND
	    POP_ATOM(K)=(RHO_ZERO/MEAN_ATOMIC_WEIGHT)*(R(ND)/R(K))**N_RHO
	    POP_ATOM(K)=POP_ATOM(K)/CLUMP_FAC(K)
	  END DO
!
	ELSE IF(VAR_MDOT)THEN
!
! Treat a model where the mass-loss might be variable as a function of time.
! e.g, For a model of AG Car. We also read in the clumping factor.
! R, V, and SIGMA should be read in with VEL_LAW option 7.
!
	  CALL RD_MOD_DENSITY(DENSITY,CLUMP_FAC,R,ND,VAR_MDOT_FILE)
	  DO K=1,ND
	    POP_ATOM(K)=DENSITY(K)/MEAN_ATOMIC_WEIGHT/ATOMIC_MASS_UNIT()
	  END DO
	ELSE
	  DO K=1,ND
	    POP_ATOM(K)=RMDOT/MEAN_ATOMIC_WEIGHT/(V(K)*R(K)*R(K)*CLUMP_FAC(K))
	  END DO
	END IF
!
	DO ISPEC=1,NUM_SPECIES
	  DO K=1,ND
	    POP_SPECIES(K,ISPEC)=AT_ABUND(ISPEC)*POP_ATOM(K)/ABUND_SUM
	  END DO
	END DO
	T1=ATOMIC_MASS_UNIT()*MEAN_ATOMIC_WEIGHT
!
	T1=ATOMIC_MASS_UNIT()*MEAN_ATOMIC_WEIGHT
	DO I=1,ND
	  DENSITY(I)=POP_ATOM(I)*T1			!gm/cm^3
	END DO
!
	CALL UPDATE_VTURB()
!	
	RETURN
!
	CONTAINS
	SUBROUTINE UPDATE_VTURB()
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	IF(VTURB_LAW .EQ. 'POW')THEN
	  IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)THEN
!
! Turbulent width is a speed, not a signed radial velocity.  Normalize
! by the full velocity range so an interior maximum is used and
! repeated/non-monotonic V values cannot make V/V(1) non-positive.
!
	    T2=MAXVAL(ABS(V(1:ND)))
	    IF(T2 .GT. 1.0E-12_LDP)THEN
	      VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*
	1        (ABS(V(1:ND))/T2)**VTURB_POW
	    ELSE
	      VTURB_VEC(1:ND)=VTURB_MIN
	    END IF
	  ELSE
	    VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*(V(1:ND)/V(1))**VTURB_POW
	  END IF
	ELSE IF(VTURB_LAW .EQ. 'TRUNC_POW')THEN
	  DO I=1,ND
	    IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)THEN
	      T1=MIN(ABS(V(I)),VTURB_VEND)/VTURB_VEND
	    ELSE
	      T1=MIN(V(I),VTURB_VEND)/VTURB_VEND
	    END IF
	    VTURB_VEC(I)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*(T1**VTURB_POW)
	  END DO
	ELSE
	  WRITE(6,*)'Error VTURB law not recognized in SET_ABUND_CLUMP'
	  STOP
	END IF
	RETURN
	END SUBROUTINE UPDATE_VTURB
	END SUBROUTINE SET_ABUND_CLUMP
