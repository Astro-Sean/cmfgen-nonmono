!
! Subroutine to perform an NG accleration for a Comoving-Frame Model. Progam
! uses the last 4 iterations which are stored in the last "4 records"
! (effectively) of SCRTEMP.
!
! The NG acceleration is perfomed separately on each depth, or over a band
! of depths.
!
! Output is to the last record of SCRTEMP.
!
! Input files required:
!                    POINT1(.DAT)
!                    SCRTEMP(.DAT)
! Output files:
!                    POINT1(.DAT)
!                    POINT2(.DAT)
!                    SCRTEMP(.DAT)
!
	SUBROUTINE DO_NG_BAND_ACCEL_V4(POPS,R,V,SIGMA,ROLD,
	1             NT,ND,NG_BAND,NG_DONE,MAXDEC,MAXINC,
	1             DO_NG_VALIDITY_CHECK,NG_ABSTAIN,LUSCR,LUER)
	USE SET_KIND_MODULE
	USE CONTROL_VARIABLE_MOD, ONLY: DO_NONMONO_SOB,DO_NONMONO_KNOP,
	1                               NG_SHADOW_DIAG
	IMPLICIT NONE
!
! Altered 05-Oct-2026 : Freeze prescribed R,V,SIGMA under non-mono.
! Altered 15-Jan-2009 : Minor bug fixed.
!                       Populations were reset to IREC-3 iteration if NG acceleration failed.
! Altered 02-May-2004 : Changed to handle new SCRTEMP file format.
!                       Allows for models in which R, V, and SIGMA are also output to
!                          SCRTEMP on every iteration.
!                       Renamed from V1: R, V, SIGMA, and ROLD added to call.
	INTEGER NT
	INTEGER ND
	INTEGER NG_BAND
	INTEGER LUSCR
	INTEGER LUER
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) ROLD(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
!
	REAL(KIND=LDP) MAXDEC
	REAL(KIND=LDP) MAXINC
	LOGICAL NG_DONE
	LOGICAL DO_NG_VALIDITY_CHECK
!
! edited 05-Oct-2026: [NG_ABSTAIN] pass-through to
! NG_MIT_BAND_OPT_V1 -- per-component hold instead of per-depth veto.
! Default F preserves the upstream behavior exactly.
!
	LOGICAL NG_ABSTAIN
!
! Local arrays. We will acelerate R, V, and SIGMA as well. If R is not changing,
! the accleration wiull keep R constant since acceleration is linear.
!
	REAL(KIND=LDP) BIG_POPS(NT+3,ND)
	REAL(KIND=LDP) RDPOPS(NT+3,ND,4)
!
	REAL(KIND=LDP) RATIO
	REAL(KIND=LDP) MEAN
	REAL(KIND=LDP) MIN_RATIO, MAX_RATIO
	REAL(KIND=LDP) T1
!
! Local variables which are adjusted to match the particular model under
! consideration.
!
	INTEGER IOS
	INTEGER IREC
	INTEGER NITSF
	INTEGER IT_STEP
	INTEGER LST_NG
	INTEGER IFLAG
!
	INTEGER NT_MOD
	INTEGER I,J
	INTEGER ICOUNT
	INTEGER, PARAMETER :: RITE_N_TIMES=1
!
! Shadow-window diagnostics (NG_SHADOW_DIAG only). Extended
! history is read into RDX through scratch arrays so POPS, R, V
! and SIGMA keep the newest-record values the caller expects.
!
	REAL(KIND=LDP), ALLOCATABLE :: RDX(:,:,:)
	REAL(KIND=LDP) POPX(NT,ND)
	REAL(KIND=LDP) RX(ND),VX(ND),SX(ND)
	LOGICAL NEWMOD_S,WRVSIG_S
	INTEGER NREC_AVA,NITSF_S,LST_NG_S
!
	LOGICAL NEWMOD
	LOGICAL WRITE_RVSIG
!
	IT_STEP=1 		!Use last 4 consecuitive iterations.
	IFLAG=0
!
! Read POPULATIONS that were output on last iteration. This is primarily
! done to get NITSF etc.
!
	NEWMOD=.FALSE.
	IREC=0
	CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                    WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	IF(NEWMOD)THEN
	  WRITE(LUER,*)'Unable to read last iteration in scratch file in DO_NG_BAND_ACCEL_V4'
	  WRITE(LUER,*)'Unrecoverable error (1) as these should have been written recently'
	  STOP
	  NG_DONE=.FALSE.
	  RETURN
	END IF
	BIG_POPS(1:NT,:)=POPS(1:NT,:)
	BIG_POPS(NT+1,:)=R(:)
	BIG_POPS(NT+2,:)=V(:)
	BIG_POPS(NT+3,:)=SIGMA(:)+1.0_LDP
!
! Read in the last 4 estimates of the poplations, as output to SCRTEMP.
! We read in the data backwards so that POPS, R, V and SIGMA contain
! the current iteration values on the last read.
!
! If the read fails, we try to recover the current POPS, R, V and SIGMA.
! We should not get a failure here.
!
	NEWMOD=.FALSE.
	DO I=4,1,-1
	  J=IREC-(I-1)*IT_STEP
	  CALL SCR_READ_V2(R,V,SIGMA,POPS,J,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  RDPOPS(1:NT,:,I)=POPS(1:NT,1:ND)
	  RDPOPS(NT+1,:,I)=R(:)
	  RDPOPS(NT+2,:,I)=V(:)
	  RDPOPS(NT+3,:,I)=SIGMA(:)+1.0_LDP
	  IF(NEWMOD)THEN
	    WRITE(LUER,*)'Unable to read scratch file in DO_NG_BAND_ACCEL_V4'
	    WRITE(LUER,*)'Trying to read IREC ',J
	    CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	    IF(NEWMOD)THEN
	      WRITE(LUER,*)'Error is not recoverable'
	      STOP
	    END IF
            NG_DONE=.FALSE.
	    RETURN
	  END IF
	END DO
!
! [DO_NONMONO_SOB/DO_NONMONO_KNOP] R, V, and SIGMA are prescribed
! structure coordinates in these modes, not unknowns for the population
! map.  Make their NG history constant so the secant fit cannot move the
! supplied non-monotonic field or reuse a stale grid from SCRTEMP.
!
	IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)THEN
	  DO I=1,4
	    RDPOPS(NT+1,:,I)=R(:)
	    RDPOPS(NT+2,:,I)=V(:)
	    RDPOPS(NT+3,:,I)=SIGMA(:)+1.0_LDP
	  END DO
	END IF
!
! [NG_SHADOW_DIAG, default F] diagnostics only (bounded
! investigation NG1, docs/CHECKLIST.md): extend the history
! window to 16 records and report, per band and window size,
! the affine-fit coefficient norm, implied prediction move,
! predicted-residual ratio, and the successive-prediction
! settle ratio that gates adaptive Ng (Ceulemans et al. 2024,
! arXiv:2411.03040, Sec 4.2). Nothing here mutates POPS,
! BIG_POPS, flags or scheduling; extra reads use scratch
! arrays. The non-mono structure freeze is applied to the
! extended records so the diagnostic sees the same history the
! fit sees. A failed extra read just shortens the window.
!
	IF(NG_SHADOW_DIAG)THEN
	  ALLOCATE(RDX(NT+3,ND,16))
	  RDX(:,:,1:4)=RDPOPS(:,:,:)
	  NREC_AVA=4
	  DO I=5,16
	    J=IREC-(I-1)*IT_STEP
	    IF(J.LT.1)EXIT
	    CALL SCR_READ_V2(RX,VX,SX,POPX,J,NITSF_S,RITE_N_TIMES,
	1                LST_NG_S,WRVSIG_S,NT,ND,LUSCR,NEWMOD_S)
	    IF(NEWMOD_S)EXIT
	    NREC_AVA=I
	    RDX(1:NT,:,I)=POPX(1:NT,:)
	    RDX(NT+1,:,I)=RX(:)
	    RDX(NT+2,:,I)=VX(:)
	    RDX(NT+3,:,I)=SX(:)+1.0_LDP
	  END DO
	  IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)THEN
	    DO I=5,NREC_AVA
	      RDX(NT+1,:,I)=R(:)
	      RDX(NT+2,:,I)=V(:)
	      RDX(NT+3,:,I)=SIGMA(:)+1.0_LDP
	    END DO
	  END IF
	  CALL NG_SHADOW_REPORT_V1(RDX,NREC_AVA,NT+3,ND,NG_BAND,
	1                LST_NG,IREC,LUER)
	  DEALLOCATE(RDX)
	END IF
!
! Now perform the NG acceleration. The accelerated estimates for the populations
! are returned in BIG_POPS.
!
	NT_MOD=NT+3
	CALL NG_MIT_BAND_OPT_V1(BIG_POPS,RDPOPS,ND,NT_MOD,NG_BAND,
	1                NG_DONE,MAXDEC,MAXINC,NG_ABSTAIN,LUER)
!
! Check to see is accleration if reasonable.
!
	IF(DO_NG_VALIDITY_CHECK .AND. NG_DONE)THEN
	  MIN_RATIO=1.0E+30_LDP; MAX_RATIO=-1.0E+30_LDP
	  IF(MAXVAL(BIG_POPS(NT,:)-POPS(NT,:)) .GT. 0.0_LDP)THEN
	    MEAN=0.0_LDP
	    ICOUNT=0
	    DO I=1,ND
	      T1=(RDPOPS(NT,I,1)-RDPOPS(NT,I,2))
! Skip zero iteration-to-iteration differences: 0/0 produced a NaN
! RATIO that silently corrupted MEAN and could mask a bad NG update.
	      IF(T1 .NE. 0.0_LDP)THEN
	        RATIO=(BIG_POPS(NT,I)-POPS(NT,I))/T1
	        MEAN=MEAN+RATIO
	        MIN_RATIO=MIN(MIN_RATIO,RATIO)
	        MAX_RATIO=MAX(MAX_RATIO,RATIO)
	        ICOUNT=ICOUNT+1
	      END IF
	    END DO
	    IF(ICOUNT .GT. 0)MEAN=MEAN/ICOUNT
	    WRITE(6,*)' '
	    IF(MEAN .GT. 5.0_LDP .AND. MIN_RATIO .GT. 0.5_LDP*MEAN .AND.
	1             MAX_RATIO .LT. 2.0_LDP*MEAN)THEN
	      WRITE(6,*)'NG acceleration will be performed'
	    ELSE
	      NG_DONE=.FALSE.
	      WRITE(6,*)'No NG acceleration will be performed'
	      WRITE(6,*)'Eiher:'
	      WRITE(6,*)'  (a) mean correction to T is less than a factor of 5 or'
	      WRITE(6,*)'  (b) variation in T corection is more than a factor or 2'
	    END IF
	    WRITE(6,'(3(3X,A))')'MIN_RATIO','     MEAN','MAX_RATIO'
	    WRITE(6,'(3ES12.2)')MIN_RATIO,MEAN,MAX_RATIO
	    WRITE(6,*)' '
	  END IF
	END IF
!
! We don't set R(1), R(ND) etc, so as to make sure these remain identical
! to the initial values., These should not change, even if we are changing the
! R grid.
!
	IF(NG_DONE)THEN
	  ROLD(1:ND)=R(1:ND)
          POPS(1:NT,:)=BIG_POPS(1:NT,:)
          R(2:ND-1)=BIG_POPS(NT+1,2:ND-1)
          V(2:ND-1)=BIG_POPS(NT+2,2:ND-1)
          SIGMA(2:ND-1)=BIG_POPS(NT+3,2:ND-1)-1.0_LDP
	END IF
!
	RETURN	
	END
!
!
!
! Routine is deisgned to allow NG acceleration over NBAND depths
! simultaneously. For use with CMFGEN.
!
	SUBROUTINE NG_MIT_BAND_OPT_V1(POPS,RDPOPS,ND,NT,NBAND,
	1                NG_DONE,MAXDEC,MAXINC,NG_ABSTAIN,LUER)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 01-Jun-2003: NUM_BAD_NG now correctly initialized to zero.
!
	INTEGER ND
	INTEGER NT
	INTEGER NBAND
	INTEGER LUER
	REAL(KIND=LDP) POPS(NT,ND)              !Initial/corrected populations (I/O)
	REAL(KIND=LDP) RDPOPS(NT,ND,4)		!Populations to be accelerated
!
! Return maximum % increase and decrease, and whether the NG acceleration
! was succesful.
!
	REAL(KIND=LDP) MAXINC
	REAL(KIND=LDP) MAXDEC
	LOGICAL NG_DONE
!
! edited 05-Oct-2026: [NG_ABSTAIN] -- when T, components
! with no contraction information or a non-physical extrapolation
! are held at their current value instead of vetoing the depth.
!
	LOGICAL NG_ABSTAIN
!
! Local arrays
!
	REAL(KIND=LDP) NEWPOP(NT,ND)
	REAL(KIND=LDP) VEC_INC(ND)
	REAL(KIND=LDP) VEC_DEC(ND)
	REAL(KIND=LDP) INT_ARRAY(ND)
!
! Local variables
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) LOCINC,LOCDEC
!
	INTEGER NS
	INTEGER NUM_BAD_NG
	INTEGER DEC_LOC,INC_LOC
	INTEGER I,K,L
	INTEGER LST,LEND
	INTEGER LOC_NBAND
	INTEGER, PARAMETER :: IONE=1
!
	NUM_BAD_NG=0
	IF(NBAND .LE. 0 .OR. NBAND .GE. ND)THEN
          LOC_NBAND=ND
	  NS=NT*ND
	  CALL DO_TRANS_AND_NG_V1(NEWPOP,RDPOPS,IONE,ND,NT,ND,NS)
	ELSE
	  DO K=ND,1,-NBAND
	    LST=MAX(K-NBAND+1,1)
	    LEND=LST+NBAND-1
	    NS=(LEND-LST+1)*NT
	    CALL DO_TRANS_AND_NG_V1(NEWPOP(1,LST),RDPOPS,LST,LEND,NT,ND,NS)
	  END DO
          LOC_NBAND=NBAND
	END IF
!
! edited 05-Oct-2026: [NG_ABSTAIN] per-component hold.
! The banded secant fit produces a single (A,B) pair per band; any
! component that does not share the band's contraction gets a bad
! extrapolation.  Upstream marks the whole depth bad when one of
! its ~NT components is bad, and >3 bad depths cancels everything.
! Two component classes are held at their current value instead:
!   - zero newest difference (RDPOPS(.,.,1).EQ.RDPOPS(.,.,2)): a
!     resacc-restored/held column or a static row (R,V,T,SIGMA+1)
!     carries no contraction information; the fitted A,B would move
!     it backwards along the record sequence;
!   - non-positive extrapolation (or non-positive current value,
!     where the acceptance ratio is undefined): the secant crossed
!     zero, which is never physical for a population.
! Held components report ratio 1, so they pass the LOCINC/LOCDEC
! test below without vetoing the depth.  Default F: unchanged.
!
	IF(NG_ABSTAIN)THEN
	  DO L=1,ND
	    DO K=1,NT
	      IF(NEWPOP(K,L) .LE. 0.0_LDP .OR.
	1        POPS(K,L) .LE. 0.0_LDP .OR.
	2        RDPOPS(K,L,1) .EQ. RDPOPS(K,L,2))
	3        NEWPOP(K,L)=POPS(K,L)
	    END DO
	  END DO
	END IF
!
! Now check whether the NG acceleation has been reasonable.
! If it has, store the new estimates in POPS.
!
	MAXINC=-1000.0_LDP
	MAXDEC=1000.0_LDP
	DO L=1,ND
	  LOCINC=-1000.0_LDP
	  LOCDEC=1000.0_LDP
	  DO K=1,NT
	    T1=NEWPOP(K,L)/POPS(K,L)
	    LOCINC=MAX(LOCINC,T1)
	    LOCDEC=MIN(LOCDEC,T1)
	  END DO
	  VEC_INC(L)=100.0_LDP*(LOCINC-1.0_LDP)
	  VEC_DEC(L)=100.0_LDP*(1.0_LDP/LOCDEC-1.0_LDP)
!
! Before storing the NG acceleration at this depth, we check to
! see whether the predicted corrections are "reasonable".
!
	  T1=NEWPOP(NT,L)-RDPOPS(NT,L,1)
	  IF(T1 .NE. 0.0_LDP)T2=(RDPOPS(NT,L,1)-RDPOPS(NT,L,2))/T1
	  IF(T1 .NE. 0 .AND. T2 .LT. 0.0_LDP)THEN
	    NUM_BAD_NG=NUM_BAD_NG+1
	    WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	    WRITE(LUER,'(A,I4)')'No accelleration perfromed at depth ',L
	    WRITE(LUER,'(A,I4)')'Temperature correction has wrong sign'
	  ELSE IF(LOCINC .GT. 10.1_LDP .OR. LOCDEC .LT. 0.09_LDP)THEN
	    NUM_BAD_NG=NUM_BAD_NG+1
	    WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	    WRITE(LUER,9000)L,LOCINC,LOCDEC
9000	    FORMAT(1X,'No NG acceleration performed at depth ',I3,/,
	1          1X,'Biggest increase was ',1PE10.2,/,
	1          1X,'Biggest decrease was ',E10.2)
	
	  ELSE
	    DO K=1,NT
	      POPS(K,L)=NEWPOP(K,L)
	    END DO
	    IF(LOCINC .GT. MAXINC)THEN
	      MAXINC=LOCINC
	      INC_LOC=L
	    END IF
	    IF(LOCDEC .LT. MAXDEC)THEN
	      MAXDEC=LOCDEC
	      DEC_LOC=L
	    END IF
	  END IF
	END DO
!
	WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	IF(NUM_BAD_NG .GT. 3)THEN
	  WRITE(LUER,*)'Too many bad NG accelerations - '//
	1              'NG acceleration cancelled'
	  NG_DONE=.FALSE.
!
! Restore old population estimates for all depths.
!
	  DO L=1,ND
	    DO K=1,NT
	      POPS(K,L)=RDPOPS(K,L,1)
	    END DO
	  END DO
	  RETURN
	END IF
!
! [NG_ABSTAIN] vacuous-acceleration guard: if every component at
! every accepted depth was held, MAXINC/MAXDEC remain at their
! no-change values (both exactly 1.0 -- a real acceleration cannot
! leave every ratio bit-identical).  Committing would append an
! identical SCRTEMP record at MAIN_COUNTER+1, polluting the next
! window with a zero-difference pair.  Report as a failed attempt;
! nothing was applied, so no restore is needed and the retry
! cadence (+8) is unchanged.
!
	IF(NG_ABSTAIN .AND. MAXINC .EQ. 1.0_LDP .AND.
	1   MAXDEC .EQ. 1.0_LDP)THEN
	  WRITE(LUER,*)'NG_ABSTAIN: all components held, '//
	1             'no acceleration applied'
	  NG_DONE=.FALSE.
	  RETURN
	END IF
!
! By definition MAXINC and MININC are both positive. These are only
! defined for successful NG accelerations (does not include any depths
! at which NG acceleration did not work).
!
	MAXINC=100.0_LDP*(MAXINC-1.0_LDP)
	MAXDEC=100.0_LDP*(1.0_LDP/MAXDEC-1.0_LDP)
	WRITE(LUER,'(A,I3)')' NG Acceleration performed using NG_BAND=',LOC_NBAND
	WRITE(LUER,'(A,I3,A,ES10.2)')
	1  ' Max NG % increase at depth ',INC_LOC,' is',MAXINC
	WRITE(LUER,'(A,I3,A,ES10.2)')
	1  ' Max NG % decrease at depth ',DEC_LOC,' is',MAXDEC
!
	NG_DONE=.TRUE.  		!Successful acceleration
!
	RETURN
	END
!
! 
!
	SUBROUTINE DO_TRANS_AND_NG_V1(NEWPOP,RDPOPS,LST,LEND,NT,ND,NS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER LST            !Start depth
	INTEGER LEND           !Final depth
	INTEGER NT
	INTEGER ND
!
! Total number of populations to be accelerated.
! This is simplye (LEND-LST+1)*NT
!
	INTEGER NS
!
	REAL(KIND=LDP) NEWPOP(NS)             !Returned with new estimates.
	REAL(KIND=LDP) RDPOPS(NT,ND,4)
!
	REAL(KIND=LDP) TEMP(4,NS)
	INTEGER I,J,K,L
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	LOGICAL WEIGHT
!
	IF(NS .NE. NT*(LEND-LST+1))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in DO_TRANS_AND_NG_V1'
	  WRITE(LUER,*)'Invalid NS dimension'
	  WRITE(LUER,*)'      NS=',NS
	  WRITE(LUER,*)'Check NS=',NT*(LEND-LST+1)
	  STOP
	END IF
!
! Use the NG acceleration method to improve the population estimates.
! We can perform the NG acceleration at each depth individually, or
! over a range of depths. Weight is used to indicate that we
! are to minimize the percentage errors - not the absolute magnitudec
! of the errors. The absolute maximum percentage change is also
! determined.
!
! Rewrite the relevant poulations in a form suitable for NGACCEL.
!
	DO J=1,4
	  DO L=LST,LEND
	    DO I=1,NT
	      K=I+NT*(L-LST)
	      TEMP(J,K)=RDPOPS(I,L,J)
	    END DO
	  END DO
	END DO
	WEIGHT=.TRUE.
!
	CALL NGACCEL(NEWPOP,TEMP,NS,WEIGHT)
!
	RETURN
	END
!
!
! Diagnostics-only shadow of the NG history fit, called from
! DO_NG_BAND_ACCEL_V4 when [NG_SHADOW_DIAG] is T (default F -
! the production path is unchanged when absent). For each
! depth band - partitioned exactly as NG_MIT_BAND_OPT_V1 - and
! each window size N=4..min(NREC,16, records since LST_NG), a
! weighted affine fit is formed:
!
!   p = sum_i c_i x_(i+1),   sum_i c_i = 1,
!
! with c the constrained least-squares solution minimizing the
! weighted norm of the implied residual sum_i c_i r_i, where
! r_i(k) = x_i(k) - x_(i+1)(k) and w(k) = 1/x_1(k), the same
! convention as NGACCEL. The Gram matrix of the residuals is
! accumulated once per band; each window size solves a leading
! (N-1)x(N-1) sub-system through LAPACK DGESV. The "previous
! running prediction" is the same fit on records 2..N - what
! the window contained one solver call earlier - so the settle
! ratio
!
!   max_rel(p_cur, p_prev) / max_rel(x_1, x_2)
!
! mirrors the Magritte adaptive-Ng acceptance test (Ceulemans
! et al. 2024, arXiv:2411.03040, Sec 4.2): settle < 1 means
! the prediction is converging faster than the ordinary
! iteration and a kick would be accepted. Also reported are
! the affine coefficient L1 norm (error amplification of the
! fit) and the relative move the prediction would apply.
! This routine prints to LUER and returns; it mutates nothing.
!
! edited 24-Nov-2026: bounded investigation NG1.
!
	SUBROUTINE NG_SHADOW_REPORT_V1(RDX,NREC,NT_MOD,ND,NBAND,
	1            LST_NG,IREC,LUER)
	USE SET_KIND_MODULE
	USE ITER_LOG_FMT_MOD, ONLY: ILF_KVB,ILF_VI
	IMPLICIT NONE
!
	INTEGER NREC            !records available in RDX (>=4)
	INTEGER NT_MOD          !rows per depth (NT+3)
	INTEGER ND
	INTEGER NBAND
	INTEGER LST_NG          !record index of last NG commit
	INTEGER IREC            !newest record index
	INTEGER LUER
!
	REAL(KIND=LDP) RDX(NT_MOD,ND,NREC)
!
	INTEGER, PARAMETER :: NMAX=16
	INTEGER LST,LEND,NW
	INTEGER K
	INTEGER IPIV(NMAX-1)
	REAL(KIND=LDP) G(NMAX-1,NMAX-1)
	REAL(KIND=LDP) GA(NMAX-1,NMAX-1)
	REAL(KIND=LDP) RHS(NMAX-1),C(NMAX-1),CP(NMAX-1)
	REAL(KIND=LDP) RN(NMAX-1)
!
! Records older than the last injection belong to a different
! iteration epoch - the adaptive window resets on a kick, so
! the explored window is capped at IREC-LST_NG as well as the
! storage bound.
!
	NW=MIN(NREC,NMAX)
	IF(LST_NG .GT. 0)NW=MIN(NW,IREC-LST_NG)
	CALL ILF_KVB(LUER,'ng shadow:',3,
	1 (/'irec   ','lst_ng','nwin  '/),
	2 (/ILF_VI(IREC),ILF_VI(LST_NG),ILF_VI(NW)/))
	IF(NW .LT. 4)RETURN
!
	IF(NBAND .LE. 0 .OR. NBAND .GE. ND)THEN
	  LST=1
	  LEND=ND
	  CALL NG_SHD_BAND_V1(RDX,NREC,NT_MOD,ND,LST,LEND,NW,
	1         G,GA,RHS,C,CP,RN,IPIV,NMAX,LUER)
	ELSE
	  DO K=ND,1,-NBAND
	    LST=MAX(K-NBAND+1,1)
	    LEND=LST+NBAND-1
	    CALL NG_SHD_BAND_V1(RDX,NREC,NT_MOD,ND,LST,LEND,NW,
	1         G,GA,RHS,C,CP,RN,IPIV,NMAX,LUER)
	  END DO
	END IF
!
	RETURN
	END
!
!
! Per-band worker for NG_SHADOW_REPORT_V1. Builds the residual
! Gram matrix over band depths LST..LEND and reports one line
! per window size. Work arrays are passed in to keep the
! caller-side band loop free of local storage.
!
! edited 24-Nov-2026: bounded investigation NG1.
!
	SUBROUTINE NG_SHD_BAND_V1(RDX,NREC,NT_MOD,ND,LST,LEND,NW,
	1         G,GA,RHS,C,CP,RN,IPIV,NMAX,LUER)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NREC,NT_MOD,ND,LST,LEND,NW,NMAX,LUER
	INTEGER IPIV(NMAX-1)
	REAL(KIND=LDP) RDX(NT_MOD,ND,NREC)
	REAL(KIND=LDP) G(NMAX-1,NMAX-1)
	REAL(KIND=LDP) GA(NMAX-1,NMAX-1)
	REAL(KIND=LDP) RHS(NMAX-1),C(NMAX-1),CP(NMAX-1)
	REAL(KIND=LDP) RN(NMAX-1)
!
	INTEGER I,J,K,L,N,NR,NR2,INFO
	REAL(KIND=LDP) W,T1,T2,DEN,CL1,RESP,DENP
	REAL(KIND=LDP) MOVE,SETN,SETD,SETTLE
	EXTERNAL DGESV
!
! Accumulate the weighted residual Gram once; leading
! sub-matrices serve every window size.
!
	G=0.0_LDP
	DO L=LST,LEND
	  DO K=1,NT_MOD
	    W=1.0_LDP
	    IF(RDX(K,L,1) .NE. 0.0_LDP)W=1.0_LDP/RDX(K,L,1)
	    DO I=1,NW-1
	      RN(I)=RDX(K,L,I)-RDX(K,L,I+1)
	    END DO
	    DO I=1,NW-1
	      DO J=I,NW-1
	        G(I,J)=G(I,J)+W*RN(I)*RN(J)
	      END DO
	    END DO
	  END DO
	END DO
	DO I=1,NW-1
	  DO J=1,I-1
	    G(I,J)=G(J,I)
	  END DO
	END DO
!
	DO N=4,NW
	  NR=N-1
	  GA=0.0_LDP
	  GA(1:NR,1:NR)=G(1:NR,1:NR)
	  RHS=0.0_LDP
	  RHS(1:NR)=1.0_LDP
	  CALL DGESV(NR,1,GA,NMAX-1,IPIV,RHS,NMAX-1,INFO)
	  IF(INFO .NE. 0)THEN
	    WRITE(LUER,'(A,I3,A,I3,A,I2,A)')
	1    'ng shadow: band ',LST,'-',LEND,' N=',N,' singular'
	    CYCLE
	  END IF
	  DEN=SUM(RHS(1:NR))
	  IF(ABS(DEN) .LT. 1.0E-60_LDP)THEN
	    WRITE(LUER,'(A,I3,A,I3,A,I2,A)')
	1    'ng shadow: band ',LST,'-',LEND,' N=',N,' den~0'
	    CYCLE
	  END IF
	  C(1:NR)=RHS(1:NR)/DEN
	  CL1=SUM(ABS(C(1:NR)))
!
! Predicted residual contraction: c^T G c relative to the
! norm of the newest residual r_1.
!
	  RESP=0.0_LDP
	  DO I=1,NR
	    T2=0.0_LDP
	    DO J=1,NR
	      T2=T2+C(J)*G(J,I)
	    END DO
	    RESP=RESP+C(I)*T2
	  END DO
	  IF(RESP .GT. 0.0_LDP .AND. G(1,1) .GT. 0.0_LDP)THEN
	    RESP=SQRT(RESP/G(1,1))
	  ELSE
	    RESP=-1.0_LDP
	  END IF
!
! Previous running prediction: same constrained fit on records
! 2..N (residuals r_2..r_(N-1), coefficients on x_3..x_N).
!
	  NR2=N-2
	  GA=0.0_LDP
	  GA(1:NR2,1:NR2)=G(2:NR2+1,2:NR2+1)
	  RHS=0.0_LDP
	  RHS(1:NR2)=1.0_LDP
	  CALL DGESV(NR2,1,GA,NMAX-1,IPIV,RHS,NMAX-1,INFO)
	  DENP=0.0_LDP
	  IF(INFO .EQ. 0)THEN
	    DENP=SUM(RHS(1:NR2))
	    CP(1:NR2)=RHS(1:NR2)/DENP
	  END IF
!
	  MOVE=0.0_LDP
	  SETN=0.0_LDP
	  SETD=0.0_LDP
	  DO L=LST,LEND
	    DO K=1,NT_MOD
	      T1=0.0_LDP
	      DO I=1,NR
	        T1=T1+C(I)*RDX(K,L,I+1)
	      END DO
	      T2=0.0_LDP
	      IF(INFO .EQ. 0 .AND. ABS(DENP) .GE. 1.0E-60_LDP)THEN
	        DO I=1,NR2
	          T2=T2+CP(I)*RDX(K,L,I+2)
	        END DO
	      END IF
	      IF(ABS(RDX(K,L,1)) .GT. 1.0E-300_LDP)
	1       MOVE=MAX(MOVE,ABS(T1-RDX(K,L,1))/ABS(RDX(K,L,1)))
	      IF(ABS(T2) .GT. 1.0E-300_LDP)
	1       SETN=MAX(SETN,ABS(T1-T2)/ABS(T2))
	      IF(ABS(RDX(K,L,2)) .GT. 1.0E-300_LDP)
	1       SETD=MAX(SETD,
	2       ABS(RDX(K,L,1)-RDX(K,L,2))/ABS(RDX(K,L,2)))
	    END DO
	  END DO
	  SETTLE=-1.0_LDP
	  IF(SETD .GT. 0.0_LDP .AND. INFO .EQ. 0 .AND.
	1   ABS(DENP) .GE. 1.0E-60_LDP)SETTLE=SETN/SETD
	  WRITE(LUER,'(A,I3,A,I3,A,I2,4(A,ES10.2))')
	1  'ng shadow: band ',LST,'-',LEND,' N=',N,
	2  ' |c|1=',CL1,' move=',MOVE,' predres=',RESP,
	3  ' settle=',SETTLE
	END DO
!
	RETURN
	END
