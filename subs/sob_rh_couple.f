C
C Non-local Sobolev escape-probability coupling for non-monotonic
C velocity fields (dev7).
C
C Implements the Rybicki & Hummer (1978, ApJ 219, 654) coupling between
C distant Sobolev points, in the geometry-free form summarized by
C Feldmeier & Nikutta (2006, A&A 446, 661).  In a non-monotonic flow a
C photon emitted at the Sobolev point r' can return to line centre at a
C distant point r where the projected gas velocity along the same ray
C matches.  The local escape probability then underestimates the true
C line mean intensity; each line receives an extra contribution
C
C    dJ(r) = 0.5 * INT dr' r' v(r') GAMH_l(r') S_l(r') K(r,r')
C
C    K(r,r') = [(1-e^-t)/t] [(1-e^-t')/t'] / SQRT[(r^2-r'^2)(r^2v'^2-r'^2v^2)]
C
C with t = GAMH_l(r)/|1+mu^2 sigma(r)|, t' the same quantity at r', and
C
C    mu^2  = (r'^2-r^2) v'^2 / (r'^2 v^2 - r^2 v'^2)
C    mu'^2 = (r'^2-r^2) v^2  / (r'^2 v^2 - r^2 v'^2)
C
C Only pairs with 0 < mu^2 <= 1 and 0 < mu'^2 <= 1 couple; the mu^2<=1
C cut also bounds the kernel away from its (integrable) edge
C singularity, since D2 >= D1 v'^2 is then guaranteed.
C
C The correction is returned as a rate-space absorption term
C
C    DRC_l = CHIL_l*dJ_l/BB_COR_l      (same units as ETAL_l*ZNET_l)
C
C which SUB_SOB_LINE_V3 subtracts inside T2 = ETAL*ZNET - DRC.
C The remote term is deliberately kept OUT of ZNET_SIM and JBAR: the
C BA Jacobian uses ZNET_SIM/POPS(NUP) directly, and under
C LAMBDA_ITERATION builds VB_SIM=-JBAR/ETAL, so both arrays feed local
C derivatives.  A lagged remote-pumping term has no local-population
C derivative; folding it in there injected an O(dJ/S) error into the
C Jacobian and destabilized SOLVEBA (verified on TypeIIn 10d).
C
C Notes and deliberate limitations:
C  - Two-point (RH78) kernel only.  Directions crossing a third
C    resonance need the extra exp(-t'') extinction of Feldmeier &
C    Nikutta; without it the coupling there is a bounded overestimate.
C  - The kernel uses the UNCLAMPED |1+mu^2 sigma|: at the coupling
C    angles t -> infinity merely damps the term, which is the correct
C    non-monotonic limit.  The Dessart, Audit & Hillier floor is a
C    patch for the local beta integral and is not applied here.
C  - S_l is evaluated with the current populations; the term is a
C    lagged (lambda-iteration) correction, so no Jacobian (VB/VC)
C    derivatives are added.
C  - GAMH_l(I) = CHIL_MAT(I,l)*3.0D-10*R(I)/V(I)/FL, matching the
C    SOBJBAR_SIM normalization, so that kappa*v_th*rho maps to
C    GAMH_l*v/r identically in code units (R in cm, V in km/s).
C
	MODULE SOB_RH_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
	SAVE
C
C Cached coupling geometry.  The tables depend only on R, V and SIGMA,
C which are fixed during an iteration sweep, so they are rebuilt only
C when any of the three input arrays changes.  The signature is an
C element-wise copy of the inputs: a scalar digest (the old SUM check)
C cannot distinguish profiles that share a sum, and it did not cover
C SIGMA at all even though RH_T1/RH_T1P cache |1+mu^2*SIGMA|.
C
	INTEGER :: RH_ND=0
	INTEGER :: RH_TOTPAIRS=0
	INTEGER :: RH_NCALL=0
	INTEGER :: RH_LASTREP=0
	INTEGER :: RH_NDOM=0
	REAL(KIND=LDP) :: RH_MAXDJ=0.0_LDP
	REAL(KIND=LDP) :: RH_MAXDJ1=0.0_LDP
	REAL(KIND=LDP) :: RH_MAXDZ=0.0_LDP
	REAL(KIND=LDP), ALLOCATABLE :: RH_RSAV(:)	!R at last build
	REAL(KIND=LDP), ALLOCATABLE :: RH_VSAV(:)	!V at last build
	REAL(KIND=LDP), ALLOCATABLE :: RH_SSAV(:)	!SIGMA at last build
	INTEGER, ALLOCATABLE :: RH_NP(:)	!# partners per depth
	INTEGER, ALLOCATABLE :: RH_ID(:,:)	!partner depth index
	REAL(KIND=LDP), ALLOCATABLE :: RH_T1(:,:)	!|1+mu^2 sigma| at I
	REAL(KIND=LDP), ALLOCATABLE :: RH_T1P(:,:)	!|1+mu'^2 sigma| at J
	REAL(KIND=LDP), ALLOCATABLE :: RH_PW(:,:)	!w_J*R_J*V_J*Upsilon1
C
	END MODULE SOB_RH_MOD
C
C
C
	SUBROUTINE SOB_RH_GEOM(R,V,SIGMA,ND,LU)
	USE SET_KIND_MODULE
	USE SOB_RH_MOD
	IMPLICIT NONE
C
C Build the coupling tables: for every depth pair (I,J) decide whether
C a photon emitted at J can return to line centre at I.  The pair
C geometry (mu^2, mu'^2, Upsilon1, integration weight) is line
C independent; only the kernel strengths depend on the line opacity.
C
	INTEGER ND,LU
	REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND)
C
	INTEGER I,J,K
	INTEGER RH_NDROP
	REAL(KIND=LDP) D1,D2,MU2,MUP2,UPS
	REAL(KIND=LDP) WJ(ND)
C
C Reuse the cached tables only when every input is bit-identical to the
C call that built them (BUG-002).  The old 1D-25-sum fingerprint could
C not see a profile change whose sum is preserved, and it never looked
C at SIGMA -- a stale T1/T1P then feeds wrong resonance taus into DRC.
C Element-wise equality is O(ND), far cheaper than the O(ND^2) rebuild.
C The RH_ND test gates the comparison: ALL() on an unallocated array is
C undefined, and the saved arrays exist only after a build.
C
	IF(RH_ND .EQ. ND)THEN
	  IF(ALL(R .EQ. RH_RSAV) .AND. ALL(V .EQ. RH_VSAV) .AND.
	1     ALL(SIGMA .EQ. RH_SSAV))RETURN
	END IF
C
	IF(ALLOCATED(RH_NP))THEN
	  DEALLOCATE(RH_NP,RH_ID,RH_T1,RH_T1P,RH_PW,
	1             RH_RSAV,RH_VSAV,RH_SSAV)
	END IF
	ALLOCATE(RH_NP(ND),RH_ID(ND,ND),RH_T1(ND,ND),
	1         RH_T1P(ND,ND),RH_PW(ND,ND))
	ALLOCATE(RH_RSAV(ND),RH_VSAV(ND),RH_SSAV(ND))
	RH_ND=ND
	RH_RSAV=R
	RH_VSAV=V
	RH_SSAV=SIGMA
	RH_NP=0
	RH_TOTPAIRS=0
	RH_NDROP=0
C
C Trapezoid weights for the r' integral on the descending grid.
C
	WJ(1)=ABS(R(2)-R(1))
	DO J=2,ND-1
	  WJ(J)=0.5_LDP*ABS(R(J+1)-R(J-1))
	END DO
	WJ(ND)=ABS(R(ND)-R(ND-1))
C
	DO I=1,ND
	  DO J=1,ND
	    IF(J .EQ. I)CYCLE
C
C RH78 assumes a monotone expanding flow: the mu^2 formulas derive
C from the Sobolev resonance condition on a same-sign profile, and
C V(J) enters the pair weight as a velocity-space measure.  A pair
C straddling a sign flip or touching a stagnation depth (V<=0) is
C outside the kernel's derivation -- drop it rather than couple an
C expanding point to an infalling one.  Dropped pairs are counted
C and reported below.  [DO_NONLOC_SOB]
C
	    IF(V(I) .LE. 0.0_LDP .OR. V(J) .LE. 0.0_LDP)THEN
	      RH_NDROP=RH_NDROP+1
	      CYCLE
	    END IF
	    D1=R(J)*R(J)-R(I)*R(I)
	    D2=R(J)*R(J)*V(I)*V(I)-R(I)*R(I)*V(J)*V(J)
	    IF(D1*D2 .LE. 0.0_LDP)CYCLE
	    MU2=D1*V(J)*V(J)/D2
	    MUP2=D1*V(I)*V(I)/D2
	    IF(MU2 .LE. 0.0_LDP .OR. MU2 .GT. 1.0_LDP)CYCLE
	    IF(MUP2 .LE. 0.0_LDP .OR. MUP2 .GT. 1.0_LDP)CYCLE
	    UPS=1.0_LDP/SQRT(D1*D2)
	    K=RH_NP(I)+1
	    RH_NP(I)=K
	    RH_ID(I,K)=J
	    RH_T1(I,K)=MAX(ABS(1.0_LDP+MU2*SIGMA(I)),1.0D-30)
	    RH_T1P(I,K)=MAX(ABS(1.0_LDP+MUP2*SIGMA(J)),1.0D-30)
	    RH_PW(I,K)=WJ(J)*R(J)*V(J)*UPS
	    RH_TOTPAIRS=RH_TOTPAIRS+1
	  END DO
	END DO
C
	WRITE(LU,'(A,I8,A,I5,A)')' SOB_RH: ',RH_TOTPAIRS,
	1  ' non-local coupling pairs over',ND,' depths'
	IF(RH_NDROP .GT. 0)WRITE(LU,'(A,I8,A)')
	1  ' SOB_RH: ',RH_NDROP,' pairs dropped at V<=0 depths',
	2  ' (non-expanding regions are not coupled)'
C
	RETURN
	END
C
C
C
	SUBROUTINE SOB_RH_COUPLE(DRC,ZNET_SIM,CHIL_MAT,ETAL_MAT,
	1                BB_COR,NEG_OPACITY,V,SIGMA,R,FL,ND,NSIM,LU)
	USE SET_KIND_MODULE
	USE SOB_RH_MOD
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
C
C Apply the RH78 coupling correction as a separate rate-space term.
C For each line l the remote-resonance contribution
C
C    dJ_l(I) = 0.5*SUM_K PW*GAMH_l(J)*S_l(J)*EXF(t_I)*EXF(t_J)
C
C is computed over the partner list of depth I, then returned as
C DRC_l = CHIL_l*dJ_l/BB_COR_l (the extra absorption rate, same units
C as ETAL_l*ZNET_l).  ZNET_SIM is read only for diagnostics and is
C not modified.
C
	INTEGER ND,NSIM,LU
	REAL(KIND=LDP) DRC(ND,NSIM)
	REAL(KIND=LDP) ZNET_SIM(ND,NSIM)
	REAL(KIND=LDP) CHIL_MAT(ND,NSIM)
	REAL(KIND=LDP) ETAL_MAT(ND,NSIM)
	REAL(KIND=LDP) BB_COR(ND,NSIM)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND),FL
	LOGICAL NEG_OPACITY(ND)
C
	INTEGER I,J,L,K
	REAL(KIND=LDP) GAMH_L(ND),S_L(ND)
	REAL(KIND=LDP) CFAC,TI,TJ,DJ,DJ1,DRATE,RH_JRAT
	REAL(KIND=LDP) EXPONX
	EXTERNAL EXPONX
C
	CALL SOB_RH_GEOM(R,V,SIGMA,ND,LU)
	DRC=0.0_LDP
	IF(RH_TOTPAIRS .EQ. 0)RETURN
C
	CFAC=3.0E-10_LDP/FL
	DO L=1,NSIM
	  DO I=1,ND
C
C |V| with a 1 km/s floor: the local Sobolev strength is a magnitude,
C and a stagnation depth (V->0) would otherwise leave Inf in the
C table.  Rows with V<=0 hold no pairs anyway (SOB_RH_GEOM drops
C them), so the floor only removes a NaN/Inf hazard.
C
	    GAMH_L(I)=CHIL_MAT(I,L)*CFAC*R(I)/MAX(ABS(V(I)),1.0_LDP)
	    IF(CHIL_MAT(I,L) .GT. 0.0_LDP .AND.
	1       ETAL_MAT(I,L) .GT. 0.0_LDP)THEN
	      S_L(I)=ETAL_MAT(I,L)/CHIL_MAT(I,L)
	    ELSE
	      S_L(I)=0.0_LDP
	    END IF
	  END DO
	  DO I=1,ND
	    IF(RH_NP(I) .EQ. 0)CYCLE
	    IF(NEG_OPACITY(I))CYCLE
	    IF(GAMH_L(I) .LE. 0.0_LDP)CYCLE
	    DJ=0.0_LDP
	    DJ1=0.0_LDP
	    DO K=1,RH_NP(I)
	      J=RH_ID(I,K)
	      IF(NEG_OPACITY(J) .OR. S_L(J) .LE. 0.0_LDP)CYCLE
	      TI=GAMH_L(I)/RH_T1(I,K)
	      TJ=GAMH_L(J)/RH_T1P(I,K)
	      DJ=DJ+RH_PW(I,K)*GAMH_L(J)*S_L(J)*EXPONX(TI)*EXPONX(TJ)
	      IF(SOB_RH_TEST_S1)THEN
	        DJ1=DJ1+RH_PW(I,K)*GAMH_L(J)*EXPONX(TI)*EXPONX(TJ)
	      END IF
	    END DO
	    DJ=0.5_LDP*DJ
	    IF(SOB_RH_TEST_S1)THEN
	      DJ1=0.5_LDP*DJ1
	      RH_MAXDJ1=MAX(RH_MAXDJ1,DJ1)
	    END IF
	    IF(DJ .LE. 0.0_LDP)CYCLE
	    IF(BB_COR(I,L) .EQ. 0.0_LDP)CYCLE
	    DRATE=DJ*CHIL_MAT(I,L)/BB_COR(I,L)
C A non-finite correction must not enter DRC silently: it poisons the
C STEQ residual for the whole depth (m15_2p82d_VAL_sob_nonloc died this
C way -- an all-NaN STEQ column at depth 11).  Report the culprit
C values and stop so the producer can be traced, instead of letting a
C NaN sail through the .LE./.GT. style guards into POPS.
	    IF(.NOT. (ABS(DRATE) .LE. HUGE(1.0_LDP)))THEN
	      WRITE(LU,'(A,2I8,A,1P3E14.6,A,I6)')
	1     ' SOB_RH: non-finite DRATE at I,L=',I,L,
	2     ' DJ,CHIL,BB_COR=',DJ,CHIL_MAT(I,L),BB_COR(I,L),
	3     ' NPART=',RH_NP(I)
	      STOP 'SOB_RH: non-finite DRATE'
	    END IF
	    DRC(I,L)=DRC(I,L)+DRATE
	    RH_NCALL=RH_NCALL+1
	    RH_MAXDZ=MAX(RH_MAXDZ,ABS(DRATE))
C Dominance census in intensity units: count cells where the remote
C dJ exceeds the local mean intensity, i.e. the lagged coupling is
C larger than what it corrects -- the artefact signature of the
C two-point kernel overestimating at multi-resonance crossings.
C The rate-space version |DRC|>|ETAL*ZNET| was dropped: ZNET~0 at
C detailed balance makes it trivially true on thick lines.
	    IF(S_L(I) .GT. 0.0_LDP .AND.
	1       ZNET_SIM(I,L) .LT. 1.0_LDP)THEN
	      RH_JRAT=DJ/((1.0_LDP-ZNET_SIM(I,L))*S_L(I)*BB_COR(I,L))
	      RH_MAXDJ=MAX(RH_MAXDJ,RH_JRAT)
	      IF(RH_JRAT .GT. 1.0_LDP)RH_NDOM=RH_NDOM+1
	    END IF
	  END DO
	END DO
C
C Self-paced strength report: the geometry rebuild trigger never
C fires when R,V are frozen by the hydro input, so report every
C fixed block of applied corrections instead.
C
	IF(RH_NCALL-RH_LASTREP .GE. 500000)THEN
	  IF(SOB_RH_TEST_S1)THEN
	    WRITE(LU,'(A,1PD10.3,A,1PD10.3,A,1PD10.3,A,I8,A,I8,A)')
	1    ' SOB_RH: max dJ/J=',RH_MAXDJ,' max dJ(S=1)=',RH_MAXDJ1,
	2    ' max dRATE=',RH_MAXDZ,' over',RH_NCALL-RH_LASTREP,
	3    ' corrections;',RH_NDOM,' dominant (dJ>J_local)'
	    RH_MAXDJ1=0.0_LDP
	  ELSE
	    WRITE(LU,'(A,1PD10.3,A,1PD10.3,A,I8,A,I8,A)')
	1    ' SOB_RH: max dJ/J=',RH_MAXDJ,' max dRATE=',RH_MAXDZ,
	2    ' over',RH_NCALL-RH_LASTREP,' corrections;',
	3    RH_NDOM,' dominant (dJ>J_local)'
	  END IF
	  RH_LASTREP=RH_NCALL
	  RH_MAXDJ=0.0_LDP
	  RH_MAXDZ=0.0_LDP
	  RH_NDOM=0
	END IF
C
	RETURN
	END
