!
! Sets the R grid, V and SIGMA from SN_HYDRO_DATA using the
! RV_SN_MODEL_SNIIN (VEL_LAW=12) grid recipe: uniform spacing in
! log(tau_es) with NBND_INS extra points packed at each boundary.
! VEL_LAW=13 differs from VEL_LAW=12 only in the data source --
! SN_HYDRO_DATA instead of input_hydro.dat -- so the SN_MODEL and
! SN_HYDRO_MODEL paths (density/abundance via RD_SN_DATA, VADAT
! RSTAR/RMAX update, requeue expansion bookkeeping) stay active.
!
! The V12 reader interpolates R(V) and passes an uninitialized
! coordinate array, which the strict LIN_INTERP checks now reject.
! Here V is instead interpolated in R by the same monotone cubic
! used in SET_RV_HYDRO_MODEL_V3, and SIGMA comes from its analytic
! derivative; this preserves the sign-changing gradient regions that
! the non-monotonic solvers need and avoids the two-sided difference
! bias near the inserted boundary cells.
!
	SUBROUTINE RV_SN_MODEL_SNHYD(R,V,SIGMA,RMAX,RCORE,SN_AGE_DAYS,ND,LU,
	1                           TAU_OPT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER LU
	INTEGER TAU_OPT
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) RMAX
	REAL(KIND=LDP) RCORE
	REAL(KIND=LDP) SN_AGE_DAYS
!
	REAL(KIND=LDP), ALLOCATABLE :: R_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: T_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: DENSITY_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: KAPPA_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: TAU_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_TAU(:)
	REAL(KIND=LDP), ALLOCATABLE :: U_HYDRO(:)
	REAL(KIND=LDP), ALLOCATABLE :: COEF(:,:)
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TAU(ND)
	REAL(KIND=LDP) KAPPA(ND)
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP), ALLOCATABLE :: LOG_R_HYDRO(:)
	LOGICAL, ALLOCATABLE :: USED(:)
	LOGICAL INS
	REAL(KIND=LDP) OLD_SN_AGE_DAYS
	REAL(KIND=LDP) DEN_SCL_FAC
	REAL(KIND=LDP) TAU_LO,DLOG_TAU,DLOG_U
	REAL(KIND=LDP) T1,T2
!
! TAU_OPT selects the opacity used for the grid coordinate:
!   0: constant electron-scattering kappa (KAPPA_ES, the VEL_LAW=12
!      input_hydro.dat convention).  For a grid uniform in LOG(tau) a
!      constant rescales all depths identically, so the value only
!      matters through depth variations it ignores (ionization state).
!   1: the Kappa (Rosseland) column already in SN_HYDRO_DATA --
!      composition- and ionization-aware per depth, which is the
!      relevant behaviour for hydrogen-poor models.  This option also
!      applies two literature-backed refinements (Dessart & Hillier
!      2010; Dessart, Audit & Hillier 2015):
!      a) a blended grid coordinate U = integral of
!         dlog_tau + ALPHA_R*|dlogR| along the hydro profile, with
!         log_tau clamped below at TAU_LO (no more than LOG_DEC_CAP
!         decades under tau_max).  A grid uniform in U is uniform in
!         log tau where the flow is optically thick and uniform in
!         log R wherever tau is nearly flat -- the thin outer
!         atmosphere, or an internal gap in a detached shell -- and
!         gets both where tau and R are jointly steep, so no fixed
!         outer zone or zone boundary is needed.  TA(1) is pinned to
!         the hydro edge so the geometry is unchanged.
!      b) structure-aware insertion: nearest uniform targets are moved
!         onto velocity-reversal entries/extrema (sigma<-1) and density
!         spikes, concentrating resolution where the solver is stiff.
!         A U_GAP_FAC floor rejects insertions that would leave the
!         cell under 5% of a uniform-cell width (adjacent features
!         bracketing a cell boundary); the output grid is then checked
!         to be strictly decreasing and the tightest cell reported.
!
	REAL(KIND=LDP), PARAMETER :: KAPPA_ES=0.34_LDP
	REAL(KIND=LDP), PARAMETER :: LOG_DEC_CAP=4.0_LDP
	REAL(KIND=LDP), PARAMETER :: RHO_SPIKE_FAC=1.5_LDP
	REAL(KIND=LDP), PARAMETER :: ALPHA_R=1.5_LDP
	REAL(KIND=LDP), PARAMETER :: U_GAP_FAC=5.0E-2_LDP
	INTEGER, PARAMETER :: NBND_INS=3
!
	INTEGER NX,NSP,MND
	INTEGER I,J,I0,IOS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=200) STRING
!
	LUER=ERROR_LU()
	IF(ND .LT. 2*NBND_INS+4)THEN
	  WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD --- ND too small for boundary insertions'
	  WRITE(LUER,*)'ND=',ND,' NBND_INS=',NBND_INS
	  STOP
	END IF
!
	OPEN(UNIT=LU,FILE='SN_HYDRO_DATA',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD --- SN_HYDRO_DATA file not found'
	    STOP
	  END IF
!
! Get the number of data points in the HYDRO model, and the number
! of species.
!
	NX=0; NSP=0; OLD_SN_AGE_DAYS=0.0_LDP
	DO WHILE (1 .EQ. 1)
	  STRING=' '
	  DO WHILE (STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING
	  END DO
	  IF(INDEX(STRING,'Number of data points:') .NE. 0)THEN
	    I=INDEX(STRING,'points:')+7
	    READ(STRING(I:),*)NX
	  ELSE IF(INDEX(STRING,'Number of mass fractions:') .NE. 0)THEN
	    I=INDEX(STRING,'fractions:')+10
	    READ(STRING(I:),*)NSP
	  ELSE IF(INDEX(STRING,'Time(days) since explosion:') .NE. 0)THEN
	    I=INDEX(STRING,'explosion:')+10
	    READ(STRING(I:),*)OLD_SN_AGE_DAYS
	  ELSE IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	     IF(NSP .EQ. 0 .OR. NX .EQ. 0 .OR. OLD_SN_AGE_DAYS .EQ. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	      WRITE(LUER,*)'NSP or NX undefined'
	      STOP
	    ELSE
	      EXIT
	    END IF
	  END IF
	END DO
!
	ALLOCATE (V_HYDRO(NX));          V_HYDRO=0.0_LDP
	ALLOCATE (SIGMA_HYDRO(NX));      SIGMA_HYDRO=0.0_LDP
	ALLOCATE (DENSITY_HYDRO(NX));    DENSITY_HYDRO=0.0_LDP
	ALLOCATE (KAPPA_HYDRO(NX));      KAPPA_HYDRO=0.0_LDP
	ALLOCATE (R_HYDRO(NX));          R_HYDRO=0.0_LDP
	ALLOCATE (LOG_R_HYDRO(NX));      LOG_R_HYDRO=0.0_LDP
	ALLOCATE (TAU_HYDRO(NX));        TAU_HYDRO=0.0_LDP
	ALLOCATE (T_HYDRO(NX));          T_HYDRO=0.0_LDP
!
! Get basic HYDRO grid vectors.
!
	DO WHILE (1 .EQ. 1)
	  DO WHILE (STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING
	  END DO
	  IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	    READ(LU,*)R_HYDRO
	  ELSE IF(INDEX(STRING,'Velocity') .NE. 0)THEN
	    READ(LU,*)V_HYDRO
	  ELSE IF(INDEX(STRING,'Sigma') .NE. 0)THEN
	    READ(LU,*)SIGMA_HYDRO
	  ELSE IF(INDEX(STRING,'Density') .NE. 0)THEN
	    READ(LU,*)DENSITY_HYDRO
	  ELSE IF(INDEX(STRING,'Kappa') .NE. 0)THEN
	    READ(LU,*)KAPPA_HYDRO
	  ELSE IF(INDEX(STRING,'Temperature') .NE. 0)THEN
	    READ(LU,*)T_HYDRO
	  ELSE IF(INDEX(STRING,'ass fraction') .NE. 0)THEN
	    IF(R_HYDRO(1) .EQ. 0.0_LDP .OR.
	1             V_HYDRO(1) .EQ. 0.0_LDP .OR.
	1             DENSITY_HYDRO(1) .EQ. 0.0_LDP)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error reading SN data in RV_SN_MODEL_SNHYD'
	      WRITE(LUER,*)'R, V, the DENSITY is zero'
	      STOP
	    ELSE
	      EXIT
	    END IF
	  END IF
	  STRING=' '
	END DO
	WRITE(LUER,*)'Successfuly read SN data in RV_SN_MODEL_SNHYD'
	CLOSE(LU)
!
! Adjust R grid of SN for expansion, exactly as in
! SET_RV_HYDRO_MODEL_V3: R=R0+V.dt and the matching density scale.
! Valid for any positive V(R); a non-monotonic velocity can violate
! the ordering, so the expanded grid is validated.
!
	T1=24.0_LDP*3600.0_LDP*1.0E+05_LDP*(SN_AGE_DAYS-OLD_SN_AGE_DAYS)/1.0E+10_LDP
	DO I=1,NX-1
	  IF(R_HYDRO(I) .LE. R_HYDRO(I+1))THEN
	    WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	    WRITE(LUER,*)'Input hydro R grid is not strictly decreasing'
	    WRITE(LUER,*)'I,R(I),R(I+1)=',I,R_HYDRO(I),R_HYDRO(I+1)
	    STOP
	  END IF
	END DO
	DO I=1,NX
	  IF(V_HYDRO(I) .LE. 0.0_LDP)THEN
	    WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	    WRITE(LUER,*)'Hydro velocity must be positive'
	    WRITE(LUER,*)'I,V(I)=',I,V_HYDRO(I)
	    STOP
	  END IF
	  T2=1.0_LDP+T1*V_HYDRO(I)/R_HYDRO(I)*(SIGMA_HYDRO(I)+1.0_LDP)
	  IF(T2 .LE. 0.0_LDP)THEN
	    WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	    WRITE(LUER,*)'Hydro density expansion factor is non-positive'
	    WRITE(LUER,*)'I,R,V,SIGMA,factor=',I,R_HYDRO(I),V_HYDRO(I),
	1                   SIGMA_HYDRO(I),T2
	    STOP
	  END IF
	  DEN_SCL_FAC=1.0_LDP/T2/
	1                   (1.0_LDP+T1*V_HYDRO(I)/R_HYDRO(I))**2
	  DENSITY_HYDRO(I)=DEN_SCL_FAC*DENSITY_HYDRO(I)
	  R_HYDRO(I)=R_HYDRO(I)+T1*V_HYDRO(I)
	END DO
	DO I=1,NX-1
	  IF(R_HYDRO(I) .LE. R_HYDRO(I+1))THEN
	    WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	    WRITE(LUER,*)'Expanded hydro R grid is not strictly decreasing'
	    WRITE(LUER,*)'I,R(I),R(I+1)=',I,R_HYDRO(I),R_HYDRO(I+1)
	    STOP
	  END IF
	END DO
	LOG_R_HYDRO(1:NX)=LOG(R_HYDRO(1:NX))
!
! Electron-scattering optical depth measured inwards from the outer
! boundary.  The trapezoid gives zero at the outer edge; the log
! gridding below needs a positive first tau, so seed it at 1% of
! the second value, matching the guard used when input_hydro.dat
! is written for VEL_LAW=12.
!
	IF(TAU_OPT .EQ. 1)THEN
! File Kappa column is in cm^2/g; floor it so a zero-opacity zone
! cannot produce duplicate tau values (LIN_INTERP needs strict
! monotonicity) while leaving the relative depth variation intact.
	  T2=MAXVAL(KAPPA_HYDRO)
	  IF(T2 .LE. 0.0_LDP)THEN
	    WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	    WRITE(LUER,*)'TAU_OPT=1 requires a positive Kappa section',
	1               ' in SN_HYDRO_DATA'
	    STOP
	  END IF
	  T2=1.0E-06_LDP*T2
	  DO I=1,NX
	    KAPPA_HYDRO(I)=MAX(KAPPA_HYDRO(I),T2)*DENSITY_HYDRO(I)*1.0E+10_LDP
	  END DO
	ELSE
	  DO I=1,NX
	    KAPPA_HYDRO(I)=KAPPA_ES*DENSITY_HYDRO(I)*1.0E+10_LDP
	  END DO
	END IF
	TAU_HYDRO(1)=0.0_LDP
	DO I=2,NX
	  TAU_HYDRO(I)=TAU_HYDRO(I-1)+0.5_LDP*(KAPPA_HYDRO(I-1)+KAPPA_HYDRO(I))*
	1                  (R_HYDRO(I-1)-R_HYDRO(I))
	END DO
	IF(TAU_HYDRO(2) .GT. 0.0_LDP)THEN
	  TAU_HYDRO(1)=0.01_LDP*TAU_HYDRO(2)
	ELSE
	  WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	  WRITE(LUER,*)'Outer electron-scattering tau is not positive'
	  STOP
	END IF
	TAU_HYDRO(1:NX)=LOG(TAU_HYDRO(1:NX))
!
! New interior grid: MND points uniform in LOG(tau_es).  TA comes
! out ordered outer to inner (decreasing R), the CMFGEN depth
! convention, so no reordering is needed.
!
	MND=ND-2*NBND_INS
	ALLOCATE (NEW_TAU(MND))
	IF(TAU_OPT .EQ. 1)THEN
! Blended-coordinate construction.  Where tau is nearly flat in R
! (the thin outer atmosphere, or an internal gap in a detached
! shell) a uniform log-tau grid cannot resolve the region: U adds
! an |dlogR| term to each cell's weight, so points spread
! geometrically through any thin zone wherever it sits.  TAU_LO
! bounds the tau contribution so the artificial outer seed cannot
! dominate U; TA(1) is pinned to the hydro edge, keeping RMAX exact.
	  TAU_LO=MAX(TAU_HYDRO(2),TAU_HYDRO(NX)-LOG_DEC_CAP)
	  ALLOCATE (U_HYDRO(NX))
	  U_HYDRO(1)=0.0_LDP
	  DO I=2,NX
	    U_HYDRO(I)=U_HYDRO(I-1)+
	1     MAX(TAU_HYDRO(I),TAU_LO)-MAX(TAU_HYDRO(I-1),TAU_LO)+
	2     ALPHA_R*(LOG_R_HYDRO(I-1)-LOG_R_HYDRO(I))
	  END DO
	  DLOG_U=U_HYDRO(NX)/REAL(MND-1,KIND=LDP)
	  DO I=1,MND
	    NEW_TAU(I)=DLOG_U*(I-1)
	  END DO
	  NEW_TAU(1)=0.0_LDP
	  NEW_TAU(MND)=U_HYDRO(NX)
!
! Structure-aware insertion: overwrite the nearest uniform target
! with the U of each velocity-reversal entry, each local sigma
! minimum inside a reversal, and each sharp density maximum.  A
! feature U always lies within half a cell of the overwritten
! target, so strict monotonicity is preserved.
	  ALLOCATE (USED(MND)); USED=.FALSE.
	  DO J=2,NX-1
	    INS=.FALSE.
	    IF(SIGMA_HYDRO(J) .LT. -1.0_LDP .AND.
	1      SIGMA_HYDRO(J-1) .GE. -1.0_LDP)INS=.TRUE.
	    IF(SIGMA_HYDRO(J) .LT. -1.0_LDP .AND.
	1      SIGMA_HYDRO(J) .LE. SIGMA_HYDRO(J-1) .AND.
	1      SIGMA_HYDRO(J) .LE. SIGMA_HYDRO(J+1))INS=.TRUE.
	    IF(DENSITY_HYDRO(J) .GT. RHO_SPIKE_FAC*DENSITY_HYDRO(J-1) .AND.
	1      DENSITY_HYDRO(J) .GT. RHO_SPIKE_FAC*DENSITY_HYDRO(J+1))INS=.TRUE.
	    IF(INS .AND. U_HYDRO(J) .GT. U_HYDRO(1) .AND.
	1      U_HYDRO(J) .LT. U_HYDRO(NX))THEN
	      I0=NINT(U_HYDRO(J)/DLOG_U)+1
	      I0=MAX(2,MIN(MND-1,I0))
! The strictness check guards the rare case where an adjacent cell
! was already overwritten with a nearly equal feature coordinate.
! The U_GAP_FAC floor bounds the collapse the ordering check alone
! allows: two features bracketing the same uniform-cell boundary
! could otherwise sit arbitrarily close, leaving a near-zero
! dlogR cell whose centred-difference denominators
! (R(I-1)-R(I+1), Euler-Maclaurin DTAU corrections) stiffen the
! exact shock feature the insertion exists to resolve.  A rejected
! feature keeps the uniform target; the cell stays ~1 cell wide.
	      IF(.NOT. USED(I0) .AND.
	1        U_HYDRO(J)-NEW_TAU(I0-1) .GT. U_GAP_FAC*DLOG_U .AND.
	2        NEW_TAU(I0+1)-U_HYDRO(J) .GT. U_GAP_FAC*DLOG_U)THEN
	        NEW_TAU(I0)=U_HYDRO(J)
	        USED(I0)=.TRUE.
	      END IF
	    END IF
	  END DO
	  DEALLOCATE (USED)
! Interpolate log R in U: uniform-U points come out uniform in log R
! through thin zones and uniform in log tau through thick ones.
	  CALL LIN_INTERP(NEW_TAU,TA(1:MND),MND,U_HYDRO,LOG_R_HYDRO,NX)
	  DEALLOCATE (U_HYDRO)
	  TA(1:MND)=EXP(TA(1:MND))
	  TA(1)=R_HYDRO(1)
	  TA(MND)=R_HYDRO(NX)
	ELSE
	  TAU_LO=TAU_HYDRO(1)
	  T1=TAU_HYDRO(NX)-TAU_LO
	  DLOG_TAU=T1/REAL(MND-1,KIND=LDP)
	  DO I=1,MND
	    NEW_TAU(I)=TAU_LO+DLOG_TAU*(I-1)
	  END DO
	  NEW_TAU(1)=TAU_LO
	  NEW_TAU(MND)=TAU_HYDRO(NX)
	  CALL LIN_INTERP(NEW_TAU,TA(1:MND),MND,TAU_HYDRO,R_HYDRO,NX)
	END IF
	DEALLOCATE (NEW_TAU)
!
	RMAX=TA(1)
	RCORE=TA(MND)
!
! Pack NBND_INS points inside the outermost and innermost interior
! cells, using the same fractions as RV_SN_MODEL_SNIIN.
!
	R(1)=TA(1)
	R(ND)=TA(MND)
	DO I=2,MND-1
	  R(I+NBND_INS)=TA(I)
	END DO
	R(2)=TA(1)-(TA(1)-TA(2))/20.0_LDP
	R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/20.0_LDP
	R(3)=TA(1)-(TA(1)-TA(2))/8.0_LDP
	R(ND-2)=R(ND)+(TA(MND-1)-TA(MND))/8.0_LDP
	R(4)=TA(1)-(TA(1)-TA(2))/3.0_LDP
	R(ND-3)=R(ND)+(TA(MND-1)-TA(MND))/3.0_LDP
!
! Grid integrity check.  LIN_INTERP is monotone in its coordinate,
! so a flat U(logR) shelf or a boundary cell degenerate in TA can
! still collapse an output spacing; centred-difference denominators
! further down (SIGMA stencil, DTAU corrections) assume strict
! ordering.  Refuse a non-decreasing grid outright and report the
! tightest surviving cell so stiff regions are visible in the log.
!
	T1=HUGE(T1)
	I0=0
	DO I=1,ND-1
	  IF(R(I) .LE. R(I+1))THEN
	    WRITE(LUER,*)'Error in RV_SN_MODEL_SNHYD'
	    WRITE(LUER,*)'Output R grid is not strictly decreasing'
	    WRITE(LUER,*)'I,R(I),R(I+1)=',I,R(I),R(I+1)
	    STOP
	  END IF
	  IF(R(I)/R(I+1)-1.0_LDP .LT. T1)THEN
	    T1=R(I)/R(I+1)-1.0_LDP
	    I0=I
	  END IF
	END DO
	IF(T1 .LT. 1.0E-4_LDP)THEN
	  WRITE(LUER,'(A,ES12.4,A,I6,A,ES12.4)')
	1   'RV_SN_MODEL_SNHYD: tightest cell dR/R-1=',T1,
	2   ' at I=',I0,' (feature insertion floor is',
	3   U_GAP_FAC,' of a uniform cell)'
	END IF
!
! Monotone-cubic interpolation of V in R with SIGMA from the cubic
! derivative (same evaluation as SET_RV_HYDRO_MODEL_V3).  Both
! grids are strictly decreasing.
!
	ALLOCATE (COEF(NX,4))
	CALL MON_INT_FUNS_V2(COEF,V_HYDRO,R_HYDRO,NX)
	J=1
	DO I=1,ND
! LIN_INTERP can return a last-ulp excursion below R_HYDRO(NX) at the
! pinned endpoint; J is capped so the final segment extrapolates
! instead of reading R_HYDRO(NX+1).
	  DO WHILE (R(I) .LT. R_HYDRO(J+1) .AND. J .LT. NX-1)
	    J=J+1
	  END DO
	  T1=R(I)-R_HYDRO(J)
	  V(I)=((COEF(J,1)*T1+COEF(J,2))*T1+COEF(J,3))*T1+COEF(J,4)
	  SIGMA(I)=((3.0_LDP*COEF(J,1)*T1+2.0_LDP*COEF(J,2))*T1+COEF(J,3))*R(I)/V(I)-1.0_LDP
	END DO
!
! If very close to Hubble law, we set SIGMA to zero.
!
	IF(MAXVAL(ABS(SIGMA)) .LT. 1.0E-04_LDP)SIGMA(1:ND)=0.0_LDP
!
! Diagnostic tau scale on the new grid (kappa_interp * dR).
!
	LOG_R(1:ND)=LOG(R(1:ND))
	CALL MON_INTERP(KAPPA,ND,1,LOG_R,ND,KAPPA_HYDRO,NX,LOG_R_HYDRO,NX)
	TAU(1)=KAPPA(1)*R(1)
	T1=LOG(KAPPA(2)/KAPPA(1))/LOG(R(1)/R(2))
	IF(T1 .GT. 1.5_LDP)TAU(1)=TAU(1)/(T1-1)
	DO I=2,ND
	  TAU(I)=TAU(I-1)+0.5_LDP*(KAPPA(I-1)+KAPPA(I))*(R(I-1)-R(I))
	END DO
!
	OPEN(UNIT=LU,FILE='NEW_SN_R_GRID',STATUS='UNKNOWN')
	  WRITE(LU,'(A,15X,A,3(5X,A))')'    I','R','V(km/s)','  SIGMA','    TAU'
	  DO I=1,ND-1
	    WRITE(LU,'(I5,ES16.8,7ES12.4)')I,R(I),V(I),SIGMA(I),TAU(I),R(I)/R(I+1),LOG10(TAU(I+1)/TAU(I))
	  END DO
	  WRITE(LU,'(I5,ES16.8,3ES12.4)')ND,R(ND),V(ND),SIGMA(ND),TAU(ND)
	CLOSE(LU)
!
	IF(TAU_OPT .EQ. 1)THEN
	  WRITE(LUER,*)'Computed blended tau-R/log-R grid in RV_SN_MODEL_SNHYD'
	ELSE
	  WRITE(LUER,*)'Computed tau-ES R grid in RV_SN_MODEL_SNHYD'
	END IF
	WRITE(LUER,'(A,ES14.6,A,ES14.6)')'  RMAX=',RMAX,' RCORE=',RCORE
!
	DEALLOCATE (R_HYDRO,LOG_R_HYDRO,V_HYDRO,SIGMA_HYDRO)
	DEALLOCATE (DENSITY_HYDRO,KAPPA_HYDRO,TAU_HYDRO,T_HYDRO,COEF)
!
	RETURN
	END
