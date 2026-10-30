	SUBROUTINE RECOMPUTE_SIGMA_TRUE(SIGMA_TRUE,SIGMA,R,V,ND)
	USE SET_KIND_MODULE
	USE CONTROL_VARIABLE_MOD, ONLY: DO_NONMONO_SOB,DO_NONMONO_KNOP,
	1                              ALLOW_SIGMA_UNDEF
	IMPLICIT NONE
!
! Created 18-Sep-2026 (dev7) : Fill SIGMA_TRUE with the unclamped
!                       velocity-gradient parameter
!                           SIGMA = dlnV/dlnR - 1 = (R/V) dV/dR - 1
!                       evaluated by centred differences on the current
!                       R and V arrays.
!
! Required by the DO_NONMONO_SOB path (Dessart, Audit & Hillier 2015,
! arXiv:1503.05463): the Sobolev optical depth must see the true
! velocity gradient, including dv/dr < 0 at shocks and cold dense
! shells.  DO_NONMONO_KNOP needs the same array: the Knop chi_hat and
! sign-dependent upwind are built from the unclamped gradient.
!
! Under DO_NONMONO_SOB the work array SIGMA is re-clamped to >= -0.999:
! the model readers (SET_RV_HYDRO_MODEL_V3 et al.) do not clamp sigma,
! so sigma < -1 leaks into the stock infrastructure (MOM_J_CMF_V11 builds
! W=GAMH*(1+sigma*NMID_ON_HMID) and EPS~1/(1+W) -- sign flips/poles),
! while the Sobolev line routines take SIGMA_TRUE with their own
! Dessart+2015 tau floor.  Clamping here, at the only refresh point,
! restores the sigma >= -1 invariant those routines assume.
! The clamp is not applied under DO_NONMONO_KNOP: the Knop solvers read
! SIGMA_TRUE for all sigma-dependent physics and the unclamped SIGMA
! behaviour is part of the verified KNOP baseline.
!
! Called at model setup and whenever R or V are revised (hydro or wind
! velocity update, NG grid revision) so SIGMA_TRUE cannot go stale.
!
! CMFGEN indexes depth outside-in: R(1) > R(2) > ... > R(ND).
!
	INTEGER ND
	REAL(KIND=LDP) SIGMA_TRUE(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
!
	INTEGER I
	INTEGER N_UNDEF
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	IF(.NOT. DO_NONMONO_SOB .AND. .NOT. DO_NONMONO_KNOP)RETURN
!
	IF(DO_NONMONO_SOB)THEN
	  DO I=1,ND
	    IF(SIGMA(I) .LT. -0.999_LDP)SIGMA(I)=-0.999_LDP
	  END DO
	END IF
!
	CALL SIGMA_TRUE_CORE(SIGMA_TRUE,R,V,ND,N_UNDEF)
!
	IF(N_UNDEF .GT. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'RECOMPUTE_SIGMA_TRUE:',N_UNDEF,
	1   ' of',ND,'cells have no usable velocity gradient;',
	2   ' SIGMA_TRUE=0 marks them non-compressing'
!
! Fail closed by default (Oct-2026 audit): an undefined cell on a
! hydro-driven grid means a collapsed spacing, and continuing would
! drop it from the compressing census silently.  VADAT key
! ALLOW_SIGMA_UNDEF=T downgrades to the former warning-only path.
!
	  IF(.NOT. ALLOW_SIGMA_UNDEF)THEN
	    WRITE(LUER,*)'Aborting: undefined SIGMA_TRUE cells make the'
	1   //' compressing census unreliable.  Set ALLOW_SIGMA_UNDEF=T'
	2   //' in VADAT to downgrade this to a warning.'
	    STOP
	  END IF
	END IF
!
	RETURN
	END
!
! 
!
	SUBROUTINE SIGMA_TRUE_CORE(SIGMA_TRUE,R,V,ND,N_UNDEF)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Shared worker: unclamped sigma for every cell, independent of the
! CONTROL_VARIABLE_MOD flags.  RECOMPUTE_SIGMA_TRUE gates on the
! nonmono flags before calling this; CMF_FLUX_V5 calls it directly so
! its EW estimators cannot drift out of sync with the transfer solve.
!
! Derivative reconstruction (30-Oct-2026 audit): sigma is computed as
! d ln|V|/d lnR - 1 = (R/V)dV/dR - 1 via the quadratic interpolant in
! lnR through the three-point stencil.  The log-space form is EXACT
! for any power law V ~ R^b on any grid (homologous SN profiles give
! sigma=0 to roundoff) and equals (R/V)dV/dR for either sign of V as
! long as the stencil carries one sign.  Where the stencil contains a
! sign change (zero crossing / shock), ln|V| is singular there and the
! code falls back to the linear-space quadratic derivative in R.
! The old two-point secant (V(I-1)-V(I+1))/(R(I-1)-R(I+1)) is the
! slope at the linear midpoint of the stencil, not at R(I); on a
! non-uniform grid that is an O(dx) bias for both forms.
!
! |V|>tol is tested (not V>tol): the algebraic form (R/V)dV/dR-1 is
! defined for V<0 and propagates the infall sign correctly into GAM.
!
! Stagnation cells |V(i)|<=VTOL keep the SIGNED V(i) in the 1/V
! factor: sigma = (R/V)dV/dR - 1 diverges with sign =
! sign(V)*sign(dV/dR), so a near-zero cell still reports the side it
! sits on correctly.  The previous T1=MIN(|V(I-1)|,|V(I+1)|)
! substitution discarded sign(V(I)) and flipped sigma's sign for
! V<0 stagnation cells.  An exactly-zero V(I) is the one case the
! sigma parametrization cannot express (the GAM product stays
! finite, ~mu^2 dV/dR/c); it is reported N_UNDEF rather than being
! assigned an arbitrary classification.
!
	INTEGER ND
	INTEGER N_UNDEF
	REAL(KIND=LDP) SIGMA_TRUE(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
!
	INTEGER I,K,N_CLAMP
	REAL(KIND=LDP) VTOL,RTOL
	REAL(KIND=LDP) HP,HM,DVDR
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL, SAVE :: VNEG_WARNED=.FALSE.
!
	N_UNDEF=0
	N_CLAMP=0
!
! ND<3: the three-point stencil needs I-1/I+1 and the boundary
! blocks read V(3)/R(ND-2), so guard the degenerate grid.
!
	IF(ND .LT. 3)THEN
	  SIGMA_TRUE(1:ND)=0.0_LDP
	  N_UNDEF=ND
	  RETURN
	END IF
!
	IF(.NOT. VNEG_WARNED .AND. MINVAL(V(1:ND)) .LT. 0.0_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'SIGMA_TRUE_CORE: negative velocities present;',
	1   '  infall path is signed-consistent but untested'
	  VNEG_WARNED=.TRUE.
	END IF
!
	VTOL=1.0E-6_LDP*MAX(MAXVAL(ABS(V(1:ND))),1.0E-30_LDP)
!
! Relative radius floor for the spacing guards: the old absolute
! 1E-10 test was unit-dependent (never fires on cm grids ~1E15).
!
	RTOL=1.0E-12_LDP
!
	DO I=2,ND-1
	  HP=R(I-1)-R(I)
	  HM=R(I)-R(I+1)
	  IF(MIN(HP,HM) .GT. RTOL*R(I))THEN
	    IF(V(I-1)*V(I+1) .GT. 0.0_LDP .AND.
	1       V(I)*V(I-1) .GT. 0.0_LDP .AND.
	1       MIN(ABS(V(I-1)),ABS(V(I)),ABS(V(I+1))) .GT.
	1       0.1_LDP*MAX(ABS(V(I-1)),ABS(V(I)),ABS(V(I+1))))THEN
!
! Stencil has one sign AND |V| varies by less than a decade across
! it: sigma = d ln|V|/d lnR - 1 via the non-uniform quadratic in
! lnR.  Exact for V ~ R^b on any grid; HP/HM reuse the log spacings.
! The decade gate keeps the log form off cells adjacent to a sign
! crossing or a sharp velocity dip, where ln|V| has a kink the
! quadratic cannot follow (the linear form is safer there and is
! exact when V is locally linear in R).
!
	      HP=LOG(R(I-1)/R(I))
	      HM=LOG(R(I)/R(I+1))
	      DVDR=LOG(ABS(V(I-1)))*HM/(HP*(HP+HM)) +
	1           LOG(ABS(V(I)))*(HP-HM)/(HP*HM) -
	1           LOG(ABS(V(I+1)))*HP/(HM*(HP+HM))
	      SIGMA_TRUE(I)=DVDR-1.0_LDP
	    ELSE IF(V(I) .NE. 0.0_LDP)THEN
!
! Sign change inside the stencil (shock or zero crossing): ln|V| is
! singular here, so use the linear-space quadratic for dV/dR.
!
	      DVDR=V(I-1)*HM/(HP*(HP+HM)) +
	1           V(I)*(HP-HM)/(HP*HM) -
	1           V(I+1)*HP/(HM*(HP+HM))
	      SIGMA_TRUE(I)=(R(I)/V(I))*DVDR-1.0_LDP
	    ELSE
	      SIGMA_TRUE(I)=0.0_LDP
	      N_UNDEF=N_UNDEF+1
	      CYCLE
	    END IF
	    SIGMA_TRUE(I)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1     SIGMA_TRUE(I)))
	    IF(ABS(SIGMA_TRUE(I)) .GE. 1.0E+6_LDP)N_CLAMP=N_CLAMP+1
	  ELSE
	    SIGMA_TRUE(I)=0.0_LDP
	    N_UNDEF=N_UNDEF+1
	  END IF
	END DO
!
! Boundary depths: second-order one-sided quadratic differences.
!
	HP=R(1)-R(2)
	HM=R(2)-R(3)
	IF(MIN(HP,HM) .GT. RTOL*R(1))THEN
	  IF(V(1)*V(2) .GT. 0.0_LDP .AND.
	1     V(2)*V(3) .GT. 0.0_LDP .AND.
	1     MIN(ABS(V(1)),ABS(V(2)),ABS(V(3))) .GT.
	1     0.1_LDP*MAX(ABS(V(1)),ABS(V(2)),ABS(V(3))))THEN
	    HP=LOG(R(1)/R(2))
	    HM=LOG(R(2)/R(3))
	    DVDR=LOG(ABS(V(1)))*(2.0_LDP*HP+HM)/(HP*(HP+HM)) -
	1         LOG(ABS(V(2)))*(HP+HM)/(HP*HM) +
	1         LOG(ABS(V(3)))*HP/(HM*(HP+HM))
	    SIGMA_TRUE(1)=DVDR-1.0_LDP
	    SIGMA_TRUE(1)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1     SIGMA_TRUE(1)))
	    IF(ABS(SIGMA_TRUE(1)) .GE. 1.0E+6_LDP)N_CLAMP=N_CLAMP+1
	  ELSE IF(V(1) .NE. 0.0_LDP)THEN
	    DVDR=V(1)*(2.0_LDP*HP+HM)/(HP*(HP+HM)) -
	1         V(2)*(HP+HM)/(HP*HM) +
	1         V(3)*HP/(HM*(HP+HM))
	    SIGMA_TRUE(1)=(R(1)/V(1))*DVDR-1.0_LDP
	    SIGMA_TRUE(1)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1     SIGMA_TRUE(1)))
	    IF(ABS(SIGMA_TRUE(1)) .GE. 1.0E+6_LDP)N_CLAMP=N_CLAMP+1
	  ELSE
	    SIGMA_TRUE(1)=0.0_LDP
	    N_UNDEF=N_UNDEF+1
	  END IF
	ELSE
	  SIGMA_TRUE(1)=0.0_LDP
	  N_UNDEF=N_UNDEF+1
	END IF
!
	HP=R(ND-1)-R(ND)
	HM=R(ND-2)-R(ND-1)
	IF(MIN(HP,HM) .GT. RTOL*R(ND))THEN
	  IF(V(ND)*V(ND-1) .GT. 0.0_LDP .AND.
	1     V(ND-1)*V(ND-2) .GT. 0.0_LDP .AND.
	1     MIN(ABS(V(ND)),ABS(V(ND-1)),ABS(V(ND-2))) .GT.
	1     0.1_LDP*MAX(ABS(V(ND)),ABS(V(ND-1)),ABS(V(ND-2))))THEN
	    HP=LOG(R(ND-1)/R(ND))
	    HM=LOG(R(ND-2)/R(ND-1))
	    DVDR=LOG(ABS(V(ND-2)))*HP/(HM*(HP+HM)) -
	1         LOG(ABS(V(ND-1)))*(HP+HM)/(HP*HM) +
	1         LOG(ABS(V(ND)))*(2.0_LDP*HP+HM)/(HP*(HP+HM))
	    SIGMA_TRUE(ND)=-DVDR-1.0_LDP
	    SIGMA_TRUE(ND)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1     SIGMA_TRUE(ND)))
	    IF(ABS(SIGMA_TRUE(ND)) .GE. 1.0E+6_LDP)N_CLAMP=N_CLAMP+1
	  ELSE IF(V(ND) .NE. 0.0_LDP)THEN
!
! The linear-space form needs the mirror of the I=1 stencil; the
! weights carry the opposite overall sign at the inner edge (the
! 30-Oct-2026 T1 failure caught a missing negation here).
!
	    DVDR=V(ND-2)*HP/(HM*(HP+HM)) -
	1         V(ND-1)*(HP+HM)/(HP*HM) +
	1         V(ND)*(2.0_LDP*HP+HM)/(HP*(HP+HM))
	    SIGMA_TRUE(ND)=(R(ND)/V(ND))*(-DVDR)-1.0_LDP
	    SIGMA_TRUE(ND)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1     SIGMA_TRUE(ND)))
	    IF(ABS(SIGMA_TRUE(ND)) .GE. 1.0E+6_LDP)N_CLAMP=N_CLAMP+1
	  ELSE
	    SIGMA_TRUE(ND)=0.0_LDP
	    N_UNDEF=N_UNDEF+1
	  END IF
	ELSE
	  SIGMA_TRUE(ND)=0.0_LDP
	  N_UNDEF=N_UNDEF+1
	END IF
!
	IF(N_CLAMP .GT. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'SIGMA_TRUE_CORE:',N_CLAMP,' of',ND,
	1   'cells clamped at |sigma|=1E6 (1/V amplification)'
	  K=0
	  DO I=1,ND
	    IF(ABS(SIGMA_TRUE(I)) .GE. 0.999999E+6_LDP)THEN
	      K=K+1
	      IF(K .LE. 10)WRITE(LUER,'(A,I6,3ES16.7)')
	1        '   clip at I,R,V,sigma=',I,R(I),V(I),SIGMA_TRUE(I)
	    END IF
	  END DO
	END IF
!
	RETURN
	END
