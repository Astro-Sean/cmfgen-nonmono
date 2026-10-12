	SUBROUTINE GET_MOMS_NON_REL(R, V, FREQ, N_TYPE, ND)
	USE SET_KIND_MODULE
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) FREQ
	CHARACTER(LEN=*) N_TYPE
!
! Local vectors
!
	REAL(KIND=LDP) TMP_VEC(1:ND)
	REAL(KIND=LDP) RSQ_JNU(1:ND)
	REAL(KIND=LDP) RSQ_JNU_MID(1:ND)
	REAL(KIND=LDP) RSQ_HNU_MID(1:ND)
	REAL(KIND=LDP) RSQ_NNU_MID(1:ND)
!
	REAL(KIND=LDP) T1,T2
	INTEGER, PARAMETER :: IONE=1
	INTEGER I
	INTEGER IOS
	INTEGER NS
	INTEGER NDM1
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	NS=ND_STORE
!
! The module arrays below persist between calls but were sized by the
! first caller.  The moment solver's expanded grid ND_NM changes on a
! NEW_R_GRID rebuild, and differs from the model ND passed by other
! callers (the USE_LAM_ES fill in COMP_J_BLANK); a stale size makes
! the RMID and MON_INTERP writes below run out of bounds.  On a size
! change, drop the whole block so the allocation re-runs at the
! current ND.  The lagged _PREV/_SAVE Eddington state is lost -- the
! same cold start as the first call -- benign next to an OOB write.
! Not flag-gated: upstream callers pass a fixed ND, so the stock path
! is unchanged.
!
	IF(ALLOCATED(K_ON_J) .AND. SIZE(K_ON_J) .NE. ND)THEN
	   DEALLOCATE(RMID)
	   DEALLOCATE(K_ON_J,NMID_ON_HMID,NMID_ON_J)
	   DEALLOCATE(K_ON_J_PREV,NMID_ON_HMID_PREV,NMID_ON_J_PREV)
	   DEALLOCATE(K_ON_J_SAVE,NMID_ON_HMID_SAVE,NMID_ON_J_SAVE)
	 END IF
	IF(.NOT. ALLOCATED(K_ON_J))THEN
	   ALLOCATE (RMID(ND),STAT=IOS)
!
	   IF(IOS .EQ. 0)ALLOCATE (K_ON_J(ND),STAT=IOS)
	   IF(IOS .EQ. 0)ALLOCATE (NMID_ON_HMID(ND),STAT=IOS)
	   IF(IOS .EQ. 0)ALLOCATE (NMID_ON_J(ND),STAT=IOS)
!
	   IF(IOS .EQ. 0)ALLOCATE (K_ON_J_PREV(ND),STAT=IOS)
	   IF(IOS .EQ. 0)ALLOCATE (NMID_ON_HMID_PREV(ND),STAT=IOS)
	   IF(IOS .EQ. 0)ALLOCATE (NMID_ON_J_PREV(ND),STAT=IOS)
!
	   IF(IOS .EQ. 0)ALLOCATE (K_ON_J_SAVE(ND),STAT=IOS)
	   IF(IOS .EQ. 0)ALLOCATE (NMID_ON_HMID_SAVE(ND),STAT=IOS)
	   IF(IOS .EQ. 0)ALLOCATE (NMID_ON_J_SAVE(ND),STAT=IOS)
	   IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Error allocating memory in GET_MOMS_NON_REL: IOS=',IOS
	     STOP
	   END IF
	   NMID_ON_HMID=0.0_LDP; NMID_ON_J=0.0_LDP; NMID_ON_HMID_PREV=0.0_LDP; NMID_ON_J_PREV=0.0_LDP
	   K_ON_J_SAVE(:)=0.0_LDP;  NMID_ON_HMID_SAVE(1:ND)=0.0_LDP; NMID_ON_J_SAVE(1:ND)=0.0_LDP
	 END IF
	 DO I=1,ND-1
	   RMID(I)=0.5_LDP*(R(I)+R(I+1))
	 END DO
!
! Compute FEDD=K/J. Since J & K are always computed at the nodes,
! this calculaton is independent of HN_DEF_AT_NODES.
!
! We only need RSQ_JNU if we are computing N(mid)/J(node)
!
	IF(ND .NE. NS)THEN
	  IF(N_TYPE .NE. 'G_ONLY')THEN
	    TMP_VEC(1:NS)=R_STORE(1:NS)*R_STORE(1:NS)*JNU_STORE(1:NS)
	    CALL MON_INTERP(RSQ_JNU,ND,IONE,R,ND,TMP_VEC,NS,R_STORE,NS)
	  END IF
!
! Degenerate-J guard (dev7, mirrors the FEDD guard in
! FG_J_CMF_NONMONO_V1): where the stored J is sub-tiny the K/J ratio
! is undefined and produced garbage Eddington factors (negative or
! ~1e-3 instead of ~1/3).  QFROMF then integrated the bad F into a
! degenerate QF (bounded at 1e30) over most of the grid, which made
! the moment tridiagonal near-singular and drove exponential J growth.
! Use the isotropic limit F=1/3 at cells where JNU_STORE is sub-tiny.
!
	  DO I=1,NS
	    IF(JNU_STORE(I) .GT. 1.0E-30_LDP)THEN
	      TMP_VEC(I)=KNU_STORE(I)/JNU_STORE(I)
	    ELSE
	      TMP_VEC(I)=1.0_LDP/3.0_LDP
	    END IF
	  END DO
	  CALL MON_INTERP(K_ON_J,ND,IONE,R,ND,TMP_VEC,NS,R_STORE,NS)
	ELSE
	  IF(N_TYPE .NE. 'G_ONLY')THEN
	    RSQ_JNU(1:NS)=JNU_STORE(1:NS)*R_STORE(1:NS)*R_STORE(1:NS)
	  END IF
	  DO I=1,NS
	    IF(JNU_STORE(I) .GT. 1.0E-30_LDP)THEN
	      K_ON_J(I)=KNU_STORE(I)/JNU_STORE(I)
	    ELSE
	      K_ON_J(I)=1.0_LDP/3.0_LDP
	    END IF
	  END DO
	END IF
!
	IF(HN_DEF_ON_NODES)THEN
!
! Compute r^2.{H,N} at the mid points. These are then used to compute
! Eddington-like factor at the grid mid points.
!
!	  NDM1=ND-1
!	  TMP_VEC(1:NS)=R_STORE(1:NS)*R_STORE(1:NS)*HNU_STORE(1:NS)
!	  CALL MON_INTERP(RSQ_HNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
!	  TMP_VEC(1:NS)=R_STORE(1:NS)*R_STORE(1:NS)*NNU_STORE(1:NS)
!	  CALL MON_INTERP(RSQ_NNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
!
	  NDM1=ND-1
!
! Same degenerate-J guard for the H/J and N/J node ratios: where
! JNU_STORE is sub-tiny the ratios are undefined and would poison the
! mid-point Eddington factors through the interpolation below.  The
! physical limit with no field is zero directional information.
!
	  DO I=1,NS
	    IF(JNU_STORE(I) .GT. 1.0E-30_LDP)THEN
	      TMP_VEC(I)=HNU_STORE(I)/JNU_STORE(I)
	    ELSE
	      TMP_VEC(I)=0.0_LDP
	    END IF
	  END DO
	  CALL MON_INTERP(RSQ_HNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
	  DO I=1,NS
	    IF(JNU_STORE(I) .GT. 1.0E-30_LDP)THEN
	      TMP_VEC(I)=NNU_STORE(I)/JNU_STORE(I)
	    ELSE
	      TMP_VEC(I)=0.0_LDP
	    END IF
	  END DO
	  CALL MON_INTERP(RSQ_NNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
	  TMP_VEC(1:NS)=JNU_STORE(1:NS)*R_STORE(1:NS)*R_STORE(1:NS)
	  CALL MON_INTERP(RSQ_JNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,NS,R_STORE,NS)
	  RSQ_HNU_MID(1:NDM1)=RSQ_HNU_MID(1:NDM1)*RSQ_JNU_MID(1:NDM1)
	  RSQ_NNU_MID(1:NDM1)=RSQ_NNU_MID(1:NDM1)*RSQ_JNU_MID(1:NDM1)
!
	ELSE IF(ND .EQ. NS)THEN

	  DO I=1,NS-1
	    RSQ_HNU_MID(I)=HNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	    RSQ_NNU_MID(I)=NNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	  END DO

	ELSE 			!ND .NE. NS
!
! Compute H/J at the nodes. We interpolate in r^2 H.
!
	  NDM1=ND-1
	  TMP_VEC(1)=HNU_AT_OB*R(1)*R(1); TMP_VEC(NS+1)=HNU_AT_IB*R(ND)*R(ND)
	  DO I=1,NS-1
	    TMP_VEC(I+1)=HNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	  END DO
	  I=NS+1
	  CALL MON_INTERP(RSQ_HNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,I,EXT_RMID_STORE,I)
!
	  TMP_VEC(1)=NNU_AT_OB*R(1)*R(1); TMP_VEC(NS+1)=NNU_AT_IB*R(ND)*R(ND)
	  DO I=1,NS-1
	    TMP_VEC(I+1)=NNU_STORE(I)*RMID_STORE(I)*RMID_STORE(I)
	  END DO
	  I=NS+1
	  CALL MON_INTERP(RSQ_NNU_MID,NDM1,IONE,RMID,NDM1,TMP_VEC,I,EXT_RMID_STORE,I)
!
	END IF
!
! Now compute the "G" Eddington factors. Because H may be zero we have to
! be careful. Depending on N_TYPE we can specify N in terms of J, N in
! terms of H, or N in terms of J and H.
!
	NDM1=ND-1
	IF(N_TYPE .EQ. 'N_ON_J')THEN
	  NMID_ON_HMID(1:NDM1)=0.0_LDP
	  DO I=1,NDM1
!
! Degenerate-J guard: where the local J is sub-tiny the N/J ratio is
! undefined; use 0 (no field -> no coupling) rather than a garbage
! quotient.
!
	    IF(RSQ_JNU(I)+RSQ_JNU(I+1) .GT. 1.0E-30_LDP)THEN
	      NMID_ON_J(I)=RSQ_NNU_MID(I)/(RSQ_JNU(I)+RSQ_JNU(I+1))
	    ELSE
	      NMID_ON_J(I)=0.0_LDP
	    END IF
	  END DO
!
! Opt-in midpoint g = N/H for the sign-split moment solver
! (GMNR_NMJ_G).  With g=0 the factor 1+sigma*g is always positive,
! so the compressing branches of MOM_J_CMF_NONMONO_V1 never engage
! and the sigma*n*J term is differenced downwind wherever sigma*n<0.
! Where |H| is resolvable use the measured ratio; where it is not
! (|H| < 1E-3|N|), the sigma*n*J piece dominates the characteristic
! and a large signed g gives 1+sigma*g the sign of sigma*N -- the
! same direction the EPS term takes.
!
	  IF(GMNR_NMJ_G)THEN
	    DO I=1,NDM1
	      IF(ABS(RSQ_HNU_MID(I)) .GT. 1.0E-3_LDP*
	1         MAX(ABS(RSQ_NNU_MID(I)),1.0E-290_LDP))THEN
	        NMID_ON_HMID(I)=RSQ_NNU_MID(I)/RSQ_HNU_MID(I)
	      ELSE
	        NMID_ON_HMID(I)=SIGN(1.0E+10_LDP,RSQ_NNU_MID(I))
	      END IF
	    END DO
	  END IF
	ELSE IF(N_TYPE .EQ. 'MIXED')THEN
	  NMID_ON_J(1:NDM1)=0.0_LDP
	  DO I=1,NDM1
	    IF(RSQ_HNU_MID(I) .NE. 0)THEN
	      T1=RSQ_NNU_MID(I)/RSQ_HNU_MID(I)
	    ELSE
	      T1=100.0_LDP
	    END IF
	    IF(T1 .GT. 1.1_LDP .OR. T1 .LT. 0.05_LDP)THEN
! Same degenerate-J guard as the N_ON_J branch above.
	      IF(RSQ_JNU(I)+RSQ_JNU(I+1) .GT. 1.0E-30_LDP)THEN
	        NMID_ON_J(I)=RSQ_NNU_MID(I)/(RSQ_JNU(I)+RSQ_JNU(I+1))
	      ELSE
	        NMID_ON_J(I)=0.0_LDP
	      END IF
	      NMID_ON_HMID(I)=0.0_LDP
	    ELSE
	      NMID_ON_HMID(I)=T1
	    END IF
	  END DO
	ELSE IF(N_TYPE .EQ. 'G_ONLY')THEN
!
! Compute G Eddington factor storing in N. We also check the validity of
! the Eddington factor in case strange values are occurring because H
! is near zero (switching sign?).
!
	  NMID_ON_J(1:NDM1)=0.0_LDP
	  DO I=1,NDM1
	    IF(RSQ_HNU_MID(I) .NE. 0)THEN
	       NMID_ON_HMID(I)=RSQ_NNU_MID(I)/RSQ_HNU_MID(I)
	       IF(NMID_ON_HMID(I) .GT. 1.0_LDP)THEN
	         NMID_ON_HMID(I)=1.0_LDP
	       ELSE IF(NMID_ON_HMID(I) .LT. 0.01_LDP)THEN
	         NMID_ON_HMID(I)=0.01_LDP
	       END IF
	    ELSE
	      NMID_ON_HMID(I)=0.0_LDP              !Later replaced by average.
	      LUER=ERROR_LU()
	      WRITE(LUER,'(1X,A,1PE16.8)')'HNU zero for frequency:',FREQ
	      WRITE(LUER,'(1X,A,I4)')'Error occurred at depth:',I
	     END IF
	  END DO
	END IF
!
	RETURN
	END
