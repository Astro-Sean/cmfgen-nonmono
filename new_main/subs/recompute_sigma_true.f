	SUBROUTINE RECOMPUTE_SIGMA_TRUE(SIGMA_TRUE,SIGMA,R,V,ND)
	USE SET_KIND_MODULE
	USE CONTROL_VARIABLE_MOD, ONLY: DO_NONMONO_SOB,DO_NONMONO_KNOP
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
! |V|>tol is tested (not V>tol): the algebraic form (R/V)dV/dR-1 is
! defined for V<0 and propagates the infall sign correctly into GAM.
! Only |V|<=tol cells take the stagnation branch.
!
	INTEGER ND
	INTEGER N_UNDEF
	REAL(KIND=LDP) SIGMA_TRUE(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
!
	INTEGER I,K,N_CLAMP
	REAL(KIND=LDP) T1,VTOL
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL, SAVE :: VNEG_WARNED=.FALSE.
!
	N_UNDEF=0
	N_CLAMP=0
!
! ND<3: the interior centred difference needs I-1/I+1 and the
! boundary blocks read V(2)/R(ND-1), so guard the degenerate grid.
! A 1-2 cell model has no usable velocity gradient; report the cells
! as undefined rather than reading out of bounds.
!
	IF(ND .LT. 3)THEN
	  SIGMA_TRUE(1:ND)=0.0_LDP
	  N_UNDEF=ND
	  RETURN
	END IF
!
! V<0 (true infall) is algebraically handled here and flagged
! compressing by the solvers' masks, but no validation run has ever
! exercised it; warn once per process rather than silently proceed.
!
	IF(.NOT. VNEG_WARNED .AND. MINVAL(V(1:ND)) .LT. 0.0_LDP)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'SIGMA_TRUE_CORE: negative velocities present;',
	1   '  infall path is signed-consistent but untested'
	  VNEG_WARNED=.TRUE.
	END IF
!
! Relative velocity floor: an absolute 1E-10 never engages for km/s
! fields (a cell at 1E-9 km/s would take the algebraic branch and
! return sigma ~ 1E9).  1E-6 of the peak speed marks cells whose
! 1/V amplification is discretisation noise rather than physics.
!
	VTOL=1.0E-6_LDP*MAX(MAXVAL(ABS(V(1:ND))),1.0E-30_LDP)
!
	DO I=2,ND-1
	  IF(ABS(V(I)) .GT. VTOL .AND.
	1     ABS(R(I-1)-R(I+1)) .GT. 1.0E-10_LDP)THEN
	    SIGMA_TRUE(I)=(R(I)/V(I))*
	1     (V(I-1)-V(I+1))/(R(I-1)-R(I+1))-1.0_LDP
!
! edited by Devin, 28-Sep-2026: the same 1E6 bound as the stagnation
! branch.  Cells just above VTOL (V ~ 1E-6*Vmax) can return sigma ~
! 1E9 on a smooth profile, which poisons GAM and chi_hat; such cells
! carry no physical gradient information anyway.
!
	    IF(ABS(SIGMA_TRUE(I)) .GT. 1.0E+6_LDP)THEN
	      SIGMA_TRUE(I)=SIGN(1.0E+6_LDP,SIGMA_TRUE(I))
	      N_CLAMP=N_CLAMP+1
	    END IF
	  ELSE IF(ABS(R(I-1)-R(I+1)) .GT. 1.0E-10_LDP)THEN
! |V|~0 (stagnation cell): the true sigma diverges as 1/V and its
! sign flips discontinuously across the zero crossing, so no choice
! is exact; scaling by the smaller neighbour velocity adopts the
! sign of dV/dR (dV/dR<0 -> compressing), matching the algebraic
! branch's classification on the V>0 side of the crossing.
! The 1E6 cap bounds the magnitude: a stagnation cell needs only to
! read "strongly compressing", not carry an unbounded value into
! GAM and the chi_hat opacity.
	    T1=MIN(ABS(V(I-1)),ABS(V(I+1)))
	    IF(T1 .GT. VTOL)THEN
	      SIGMA_TRUE(I)=(R(I)/T1)*
	1       (V(I-1)-V(I+1))/(R(I-1)-R(I+1))-1.0_LDP
	      SIGMA_TRUE(I)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1       SIGMA_TRUE(I)))
	    ELSE
	      SIGMA_TRUE(I)=0.0_LDP
	      N_UNDEF=N_UNDEF+1
	    END IF
	  ELSE
	    SIGMA_TRUE(I)=0.0_LDP
	    N_UNDEF=N_UNDEF+1
	  END IF
	END DO
!
! Boundary depths: one-sided differences.
!
	IF(ABS(V(1)) .GT. VTOL .AND.
	1   ABS(R(1)-R(2)) .GT. 1.0E-10_LDP)THEN
	  SIGMA_TRUE(1)=(R(1)/V(1))*
	1   (V(1)-V(2))/(R(1)-R(2))-1.0_LDP
	  IF(ABS(SIGMA_TRUE(1)) .GT. 1.0E+6_LDP)THEN
	    SIGMA_TRUE(1)=SIGN(1.0E+6_LDP,SIGMA_TRUE(1))
	    N_CLAMP=N_CLAMP+1
	  END IF
	ELSE IF(ABS(R(1)-R(2)) .GT. 1.0E-10_LDP .AND.
	1        ABS(V(2)) .GT. VTOL)THEN
	  SIGMA_TRUE(1)=(R(1)/ABS(V(2)))*
	1   (V(1)-V(2))/(R(1)-R(2))-1.0_LDP
	  SIGMA_TRUE(1)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1   SIGMA_TRUE(1)))
	ELSE
	  SIGMA_TRUE(1)=0.0_LDP
	  N_UNDEF=N_UNDEF+1
	END IF
!
	IF(ABS(V(ND)) .GT. VTOL .AND.
	1   ABS(R(ND-1)-R(ND)) .GT. 1.0E-10_LDP)THEN
	  SIGMA_TRUE(ND)=(R(ND)/V(ND))*
	1   (V(ND-1)-V(ND))/(R(ND-1)-R(ND))-1.0_LDP
	  IF(ABS(SIGMA_TRUE(ND)) .GT. 1.0E+6_LDP)THEN
	    SIGMA_TRUE(ND)=SIGN(1.0E+6_LDP,SIGMA_TRUE(ND))
	    N_CLAMP=N_CLAMP+1
	  END IF
	ELSE IF(ABS(R(ND-1)-R(ND)) .GT. 1.0E-10_LDP .AND.
	1        ABS(V(ND-1)) .GT. VTOL)THEN
	  SIGMA_TRUE(ND)=(R(ND)/ABS(V(ND-1)))*
	1   (V(ND-1)-V(ND))/(R(ND-1)-R(ND))-1.0_LDP
	  SIGMA_TRUE(ND)=MAX(-1.0E+6_LDP,MIN(1.0E+6_LDP,
	1   SIGMA_TRUE(ND)))
	ELSE
	  SIGMA_TRUE(ND)=0.0_LDP
	  N_UNDEF=N_UNDEF+1
	END IF
!
	IF(N_CLAMP .GT. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'SIGMA_TRUE_CORE:',N_CLAMP,' of',ND,
	1   'cells clamped at |sigma|=1E6 (1/V amplification)'
!
! Clipped cells are where the discretised gradient carried no
! physical information; knowing WHERE they sit (a shock front vs a
! stagnation point) matters for interpreting the solve, so report
! the first few instead of only the count.
!
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
