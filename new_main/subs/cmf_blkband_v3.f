!
! Subroutine to solve a "Block-Banded" system of simultaneous equations
!                 BA. X. = STEQ
! The solution X is returned in STEQ. BA is corrupted.
!
! The "BANDED" matrix can either be 'Diagonal' or 'Tridiagonal'
!
! Routine is currently designed to operate on the SMALL variation matrix
! which is of dimension N.NIV.NUM_BNDS.ND where NUM_BNDS refers to the number
! of bands.
!
! A integer matrix, LNK_IV_TO_F indicates how BA is expanded into BA_BIG,
! which has dimension N,N,NUM_BNDS,ND.
!
! This routine was desined specifically to handle BA_SM, while conserving
! memory.
!
	SUBROUTINE CMF_BLKBAND_V3(STEQ,POPS,SOL_TYPE,FLAG,
	1                           DIAG_INDX,N,NION,NUM_BNDS,ND,
	1                           BA_COMPUTED,WR_BA_INV,WR_PRT_INV)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 01-Mar-2023: Removed superfulous IFAIL variable from calls to MAT5PEN.
! Altered 15-Nov-2021: Changed RUB to matrix to avoid scale/matrix conflicts in calls.
! Altered 28-Jan-2001: Based on CMF_BLKBAND_V2
!                      BA_COMPUTED and WR_BA_INV included in call.
!                      Routine will write out INVERSE of BA matrix if requested.
!                      Routine will also run on LINUX system with limited memory.
! Created 23-Mar-2001: Based on BLKBAND (See BLKBAND for earlier alterations.
!
! 
! The description here is from BLKBAND, which also allowed for a PENTDIAGONAL
! matrix. FOr simplicity,w e retained the smae notation.
!
! Let A[k], B[k], C[k], D[k], and E[k] be the sub-matrices (dimension N*N) of
! the large Block-Pentadiagonal matrix. Then
!                               A[k]=BA( , ,DI-2,K)
!                               B[k]=BA( , ,DI-1,K)
!                               C[k]=BA( , ,DI,K)
!                               D[k]=BA( , ,DI+1,K)
!                               E[k]=BA( , ,DI+2,K)
!
! where DI=K if ND=NUM_BNDS, else DI=DIAG_INDX.
!
! Likewise we have              L[k] =STEQ( , K)
!
! Elimination defined by:
!                  phi[1]= C^{-1) . L[1]
!                  alpha[1]= C^{-1) . D[1]
!                  gamma[1]= C^{-1) . E[1]
!                  SOL[1] = phi[1] - alpha[1].SOL[2] - gamma[1].SOL[3]
! Then
!                  m[k]=C[k]-B[k].alpha[k-1] - A[k].( gamma[k-2]
!                                            - alpha[k-2].alpha[k-1] )
!                  phi[k]=m[k]^{-1}(L[k] - B[k].phi[k-1]- A[k].( phi[k-2]
!                                            - alpha[k-2].phi(k-1) )
!                  alpha[k]=m[k]^{-1}(D[k]-B[k].gamma[k-1] +
!                                            + A[k].alpha[k-2].gamma[k-1]
!                  gamma[k]=m[k]^{-1} . E[k]
! with
!                  SOL[k]= phi[k] - alpha[k]. SOL[k+1] - gamma[k].SOL[k+2]
!
! Note that for k=2, we take A[k]==0.
!           for k=ND-1, we hace E[ND-1]==0 and hence gamma[ND-1]==0
!           for k=ND, we have D[ND]=E[ND]=0 and hence alpha[ND]=gamma[ND]==0.
! Thus
!                  m[ND]=C[ND]-B[ND].alpha[ND-1] - A[ND].( gamma[ND-2]
!                                            - alpha[ND-2].alpha[ND-1] )
!                  phi[ND]=m[ND]^{-1}(L[ND] - B[ND].phi[ND-1]
!                                            - A[ND].( phi[ND-2]
!                                            - alpha[ND-2].phi(ND-1) )
! and              SOL(ND)=phi[ND]
!
! The back substitution is simly obtained from the exprssion for SOL[k]
! since SOL[ND] is now known.
!
! The tridiagonal case can be obtained from these expressinons by setting
! A[k]=E[k]=0, and A[1]=D[ND]=0. When the tridiagonal option is indicated,
! these matrices are ASSUMED to be zero.
! 
!**************************************************************************
! SUBROUTINES CALLED
! ******************
!
!          DGEMV(TRANS,M,N,ALPHA,A,IDA,X,INCX,BETA,Y,INCY)
!
! BLAS2 routine.
!
! Compute Y <--- alpha A.X + beta Y
!
! Where
!
!  	A = Double preciosion A(IDA,N)
!       M = Number of rows in A
!	N = Number of columns in A
!       IDA = First dimension of A [ > max(1,m) ]
!
!	ALPHA, BETA are REAL
!
!       X=REAL   X(1+(N-1}INCX) if TRANS='N'
!       X=REAL   X(1+{M-1}INCY) if TRANS='T'
!	INCX = Stride for X (>0, generally 1)
!
!       Y=REAL   Y(1+(M-1}INCX) if TRANS='N'
!       Y=REAL   Y(1+(N-1}INCX) if TRANS='T'
!	INCY = Stride for Y (>0, generally 1)
! 
! ****************************************************************************
!
!          DGETRF(M,N,A,LDA,IPRIV,INFO)
!
! LAPACK routine.
!
! Perform the LU decomposition of the general matrix A where:
!
!		M   = Number of rows in A
!		N   = Number of columns in A
!       	IDA = First dimension of A [ > max(1,m) ]
!  		A   = Double precision A(IDA,N)
!
!	      IPIV  = Integer work vector with pivot
!                   = IPIV(min[m,n])
!
!             Integer INFO = 0 (successful exit)
!                          = -i (i th argument has illeagal value)
!                          = i (Pivot A(i,i) is exactly zero)
!
! LU decomposition is stored in A.
!
!******************
!
!                   DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
!
! LAPACK routine.
!
! Solves for the solution of A^{-1} . B . DGETRF must have
! previously been called, and A and IPIV must be passed to DGETRS unchanged.
! The solution is returned in B.
!
!		N      = Order of A
!		NRHS   = Number of RHS
!       	IDA    = First dimension of A [ > max(1,m) ]
!  		A      = Double precision A(IDA,N)
!
!	      IPIV  = Integer work vector with pivot
!                   = IPIV(min[m,n])
!
!  		B      = Double precision B(IDB,NRHS)
!
!             Integer INFO = 0 (successful exit)
!                          = -i (i th argument has illeagal value)
! 
!*******************
!
!	    DGEEQU(M,N,A,LDA,ROW_SF,COL_SF,ROW_CND,COL_CND,MAX_VAL,IFAIL)
!
! LAPACK routine
!
! Compute the Equibriation factors for the M by N matrix A whose first
! dimension is LDA. These factors help improve the condition of the matrix.
!
!
! The equilibrated matrix is found from
!
!         E(I,J)=ROW_SF(I)*A(I,J)*COL_SF(J)
!
! To solve  a set of simultaneous equtions multiply RHS by ROW_SF before
! solution, and by COL_SF after solution.
!
!
!*******************
!                  MAT5PEN(A,B,C,D,E,F,WRKMAT,VEC,N,NS,ENONZERO,DNONZERO)
!
! Evaluates A = A - B. C + D. (F . C - E)
!
! where:
!               A = real A[N,NS] - Retuned with solution
!               B = real B[N,N]  - Unchanged
!               C = real C[N,NS] - Unchanged
!               D = real D[N,N]  - Unchanged
!               E = real E[N,NS] - Unchanged
!               F = real F[N,N] -  Unchanged
!               WRKMAT = real A[N,NS] - Working matrix (corrupted)
!               VEC = real VEC[N] - Working vector (passed but no longer used)
!               ENONZERO - If FALSE, E is assumed to be zero, and E is not
!                             accessed.
!               DNONZERO - If FALSE, D is assumed to be zero, and D, F and E
!                             are not accessed.
!
!
!******************************************************************************
! 
!
	INTEGER N,NION,ND,NUM_BNDS,DIAG_INDX
	CHARACTER*(*) SOL_TYPE
	REAL(KIND=LDP) STEQ_STORE(N,ND)
        REAL(KIND=LDP) STEQ(N,ND)
        REAL(KIND=LDP) POPS(N,ND)
	LOGICAL FLAG
!
	LOGICAL BA_COMPUTED
	LOGICAL WR_BA_INV
	LOGICAL WR_PRT_INV
!
! Local variables
!
        REAL(KIND=LDP), ALLOCATABLE :: B_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: C_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: D_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: ORIG_POPS(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: PREV_D_MAT(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: D_MAT_STORE(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: RUB(:,:)	      	      !Not accessed when passed.
!
        REAL(KIND=LDP) ROW_SF(N)
        REAL(KIND=LDP) COL_SF(N)
        REAL(KIND=LDP) ROW_CND,COL_CND,MAX_VAL
	INTEGER IPIVOT(N)
	INTEGER RUB_VEC(N)
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL REPLACE_EQ(NION)
	LOGICAL ZERO_STEQ(N)
	LOGICAL USE_PASSED_REP
!
	INTEGER DEPTH_INDX
	INTEGER BAND_INDX
!
        INTEGER I,J,K,JJ
	INTEGER N_DIV_FLOOR
	REAL(KIND=LDP) POP_DIV
        INTEGER IOS,IFAIL
	CHARACTER*10 DESC
	LOGICAL ANONZERO,ENONZERO
	LOGICAL FIRST_MATRIX,LAST_MATRIX
	LOGICAL WR_D_MAT
!
        REAL(KIND=LDP),      PARAMETER :: DP_NEG_ONE=-1.0_LDP
        REAL(KIND=LDP),      PARAMETER :: DP_ONE=1.0_LDP
        INTEGER,   PARAMETER :: INT_ONE=1
        INTEGER,   PARAMETER :: NSNG=1
        CHARACTER*1, PARAMETER :: NO_TRANS='N'
	CHARACTER(LEN=50) TMP_STR
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	LOGICAL BCD_OK,BCD_TRI_OK
	CHARACTER(LEN=200) BCD_DIR,BCD_F1,BCD_F2
!
	LUER=ERROR_LU()
	N_DIV_FLOOR=0
!
! The DIAG recall path reads the cached LU factors (CSCRATCH<depth>)
! written by WRITE_BCD_MAT.  A stored BAMAT can outlive the cache:
! warm-restart staging, directory refreshes, and solves that factor
! via the JFNK path all write BAMAT/BAMATPNT without the BCD files.
! READ_BCD_MAT STOPs on a missing file, so verify the boundary
! depths first and fall through to the compute path, which rebuilds
! the factors from the loaded SE%BA.  The check mirrors the startup
! guard in cmfgen_sub.f but covers mid-run cache loss as well.
!
	BCD_OK=.FALSE.
	IF(SOL_TYPE(1:4) .EQ. 'DIAG' .AND. .NOT. BA_COMPUTED .AND.
	1   WR_BA_INV)THEN
	  CALL GET_ENVIRONMENT_VARIABLE('SCRATCH_DIR',BCD_DIR)
	  IF(LEN_TRIM(BCD_DIR) .EQ. 0)BCD_DIR='SCRATCH/'
	  WRITE(BCD_F1,'(A,A,I4.4)')TRIM(BCD_DIR),'CSCRATCH',1
	  WRITE(BCD_F2,'(A,A,I4.4)')TRIM(BCD_DIR),'CSCRATCH',ND
	  INQUIRE(FILE=TRIM(BCD_F1),EXIST=BCD_OK)
	  IF(BCD_OK)INQUIRE(FILE=TRIM(BCD_F2),EXIST=BCD_OK)
	  IF(.NOT. BCD_OK)THEN
	    WRITE(LUER,*)'CMF_BLKBAND_V3: BAMAT recall requested but'//
	1      ' CSCRATCH factor cache is absent; recomputing LU'//
	2      ' from stored BA'
	  END IF
	END IF
!
	BCD_TRI_OK=.FALSE.
	IF(SOL_TYPE(1:3) .EQ. 'TRI' .AND. .NOT. BA_COMPUTED .AND.
	1   WR_BA_INV)THEN
	  CALL GET_ENVIRONMENT_VARIABLE('SCRATCH_DIR',BCD_DIR)
	  IF(LEN_TRIM(BCD_DIR) .EQ. 0)BCD_DIR='SCRATCH/'
	  WRITE(BCD_F1,'(A,A,I4.4)')TRIM(BCD_DIR),'BCSCRATCH',1
	  WRITE(BCD_F2,'(A,A,I4.4)')TRIM(BCD_DIR),'DSCRATCH',ND
	  INQUIRE(FILE=TRIM(BCD_F1),EXIST=BCD_TRI_OK)
	  IF(BCD_TRI_OK)INQUIRE(FILE=TRIM(BCD_F2),EXIST=BCD_TRI_OK)
	  IF(.NOT. BCD_TRI_OK)THEN
	    WRITE(LUER,*)'CMF_BLKBAND_V3: BAMAT recall requested but'//
	1      ' BCSCRATCH/DSCRATCH factor cache is absent;'//
	2      ' recomputing LU from stored BA'
	  END IF
	END IF
!
!
!
! If the BA matrix was not computed, it means that its inverse must have been previously
! computed. IF WR_BA_INV is true, we can use the previous inverse, which is stored on disk.
!
	IF(SOL_TYPE(1:4) .EQ. 'DIAG' .AND. .NOT. BA_COMPUTED .AND.
	1   WR_BA_INV .AND. BCD_OK)THEN
	  FIRST_MATRIX=.TRUE.
	  LAST_MATRIX=.FALSE.
	  USE_PASSED_REP=.TRUE.
          ALLOCATE (C_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (D_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RUB(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ORIG_POPS(N,ND),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_BLKBAND'
            WRITE(LUER,*)'Unable to allocate C_MAT etc'
            WRITE(LUER,*)'STAT=',IOS
            STOP
          END IF
	  C_MAT=0.0_LDP; D_MAT=0.0_LDP; ORIG_POPS=0.0_LDP
!
	  DO K=1,ND
!
! Read in  LU decomposition, and other necessary data.
!
	    DEPTH_INDX=K
	    CALL READ_BCD_MAT(RUB,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1          ORIG_POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'C')
!
! Map the small BA rray onto the full BA array (one depth at a time).
! Really only need to do this for STEQ. We pass D_MAT, since we don't
! want to corrupt C_MAT. D_MAT is not used.
!
	    IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	    CALL GENERATE_FULL_MATRIX_V3(
	1         D_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1         N,ND,NION,NUM_BNDS,
	1         DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    FIRST_MATRIX=.FALSE.
!
	    STEQ_STORE(:,K)=STEQ(:,K)
!
! Correct STEQ for the eqilibrization of C
!
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	    END DO
!
! Now perform the solution.
!
	    CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ(1,K),N,IFAIL)
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	    END DO
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_DIAG'
	      GOTO 9999
	    END IF
!
! Correct for the use of old POPS populations in doing the scaling.
!
	   DO I=1,N
!
! A transiently non-positive POPS (rejected pre-step residue)
! divides to Inf/NaN here; floor it at a fraction of the saved
! pre-step population so the correction stays finite, and count
! the intervention.
!
	     POP_DIV=MAX(POPS(I,K),
	1        1.0E-40_LDP*MAX(ABS(ORIG_POPS(I,K)),1.0E-300_LDP))
	     IF(POP_DIV .NE. POPS(I,K))THEN
	       N_DIV_FLOOR=N_DIV_FLOOR+1
	       FLAG=.TRUE.
	     END IF
	     STEQ(I,K)=STEQ(I,K)*(ORIG_POPS(I,K)/POP_DIV)
	   END DO
!
	  END DO
          FLAG=.TRUE.
          DEALLOCATE (C_MAT,D_MAT,ORIG_POPS,RUB)
	  CALL WR2D_V2(STEQ_STORE,N,ND,'STEQ_ARRAY','*',L_TRUE,16)
!
	  RETURN
	END IF
!
!
! Check to see if we only require the diagonal solution. For this special
! case, the solutions at each depth are independent. We must compute the
! inverse, which will be saved on DISK if WR_BA_INV is true.
!
	IF(SOL_TYPE(1:4) .EQ. 'DIAG')THEN
!
	  FIRST_MATRIX=.TRUE.
	  LAST_MATRIX=.FALSE.
	  USE_PASSED_REP=.FALSE.
          ALLOCATE (C_MAT(N,N),STAT=IOS)
          ALLOCATE (RUB(N,N),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_BLKBAND'
            WRITE(LUER,*)'Unable to allocate C_MAT etc'
            WRITE(LUER,*)'STAT=',IOS
            STOP
	  END IF
	  C_MAT=0.0_LDP
!
          DO K=1,ND
!
! Map the small BA rray onto the full BA array (one depth at a time).
!
	    DEPTH_INDX=K
	    IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	    CALL GENERATE_FULL_MATRIX_V3(
	1         C_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1         N,ND,NION,NUM_BNDS,
	1         DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    FIRST_MATRIX=.FALSE.
!
	    STEQ_STORE(:,K)=STEQ(:,K)
!
!	  WRITE(TMP_STR,'(I2)')K; TMP_STR='BA_ASCI_N_D'//ADJUSTL(TMP_STR)
!	  OPEN(UNIT=96,FILE=TMP_STR,STATUS='UNKNOWN')
!           CALL WR2D_MA(POPS(1,K),N,1,'POPS_D1',96)
!           CALL WR2D_MA(STEQ(1,K),N,1,'STEQ_VEC_D1',96)
!            CALL WR2D_MA(C_MAT,N,N,'C_MAT_D1',96)
!         CLOSE(UNIT=96)
!
! Perform the LU decomposition using DGETRF. We first equilibrize the matrix
! using DGEEQU sot the the maximum row and column values are approximately
! unity.
!
	    CALL TUNE(1,'DGEEQU')
	    CALL DGEEQU(N,N,C_MAT,N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	    CALL TUNE(2,'DGEEQU')
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
              DO I=1,N
	        C_MAT(I,J)=C_MAT(I,J)*ROW_SF(I)*COL_SF(J)
	      END DO
	    END DO
	    CALL TUNE(1,'DGERTF')
	    CALL DGETRF(N,N,C_MAT,N,IPIVOT,IFAIL)
	    CALL TUNE(2,'DGERTF')
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRF_DIAG'
	      WRITE(LUER,*)'Error in CMF_BLKBAND_V3'
	      WRITE(LUER,*)'Unable to get solution at depth',K
	      WRITE(LUER,*)'Setting fractional corrections to zero'
	      WRITE(LUER,*)'IFAIL=',IFAIL
	      STEQ(:,K)=0.0_LDP
	      GOTO 500			!9999
	    END IF
	    IF(WR_BA_INV)THEN
	      CALL WRITE_BCD_MAT(RUB,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1              POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'C')
	    END IF
!
! Now perform the solution.
!
	    CALL TUNE(1,'DGETRS')
	    CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ(1,K),N,IFAIL)
	    CALL TUNE(2,'DGETRS')
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	    END DO
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_DIAG'
	      WRITE(LUER,*)'Error in CMF_BLKBAND_V3'
	      WRITE(LUER,*)'Unable to get solution at depth',K
	      WRITE(LUER,*)'Setting fractional corrections to zero'
	      WRITE(LUER,*)'IFAIL=',IFAIL
	      STEQ(:,K)=0.0_LDP
	    END IF
500	    CONTINUE
	  END DO
          FLAG=.TRUE.
          DEALLOCATE (C_MAT,RUB)
	  CALL WR2D_V2(STEQ_STORE,N,ND,'STEQ_ARRAY','*',L_TRUE,16)
	  RETURN
	END IF
!
! 
!
! If we get here, we must be performing a TRI diagonal solution.
!
	IF(SOL_TYPE(1:3) .EQ. 'TRI')THEN
	  IF(NUM_BNDS .LT. 3)THEN
	    WRITE(LUER,*)'Error in CMF_BLKBAND : NUM_BNDS too small ',
	1             'for solution type'
	    FLAG=.FALSE.
	    STOP
	  END IF
	ELSE
	  WRITE(LUER,*)'Error in CMF_BLKBAND - invalid SOL_TYPE - ',
	1          'SOLTYPE= ',SOL_TYPE
	  STOP
        END IF
!
! Solve the TRIDIAGONAL system of equations assuming that we have already performed
! an LU decomposition on the BA matrix.
!
	IF(SOL_TYPE(1:3) .EQ. 'TRI' .AND. .NOT. BA_COMPUTED .AND.
	1   WR_BA_INV .AND. BCD_TRI_OK)THEN
!
! Allocate needed work arrays. These have been saved, 1 depth at a time.
!
          ALLOCATE (B_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (C_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (D_MAT(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RUB(N,N),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ORIG_POPS(N,ND),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_BLKBAND'
            WRITE(LUER,*)'Unable to allocate B_MAT, C_MAT & D_MAT'
            WRITE(LUER,*)'STAT=',IOS
            STOP
          END IF
	  B_MAT=0.0_LDP; C_MAT=0.0_LDP; D_MAT=0.0_LDP; ORIG_POPS=0.0_LDP
!
! Needed for MAT5PEN and GENERATE_FULL_MATRIX.
!
          ANONZERO=.FALSE.
          ENONZERO=.FALSE.
	  FIRST_MATRIX=.TRUE.
	  LAST_MATRIX=.FALSE.
	  USE_PASSED_REP=.TRUE.
! 
!
! Do the forward elimination etc.
!
	  DO K=1,ND
	    DEPTH_INDX=K
!
	    CALL READ_BCD_MAT(B_MAT,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1           ORIG_POPS(1,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'BC')
!
! Map the small BA array onto the full BA array (one depth at a time).
! We just need STEQ as a modified C_MAT will be read in. To avoid
! corruptin C_MAT, we use D_MAT.
!
	    IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	    CALL GENERATE_FULL_MATRIX_V3(
	1         D_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1         N,ND,NION,NUM_BNDS,
	1         DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	    FIRST_MATRIX=.FALSE.
!
	    STEQ_STORE(:,K)=STEQ(:,K)
!
! Computes phi[k]' (i.e. ' denotes not yet multiplied by m[k]^{-1} )
! Stored in L[k].
!
	    IF(K .NE. 1)THEN
	      ENONZERO=.FALSE.
	      ANONZERO=.FALSE.
	      CALL MAT5PEN(STEQ(1,K),B_MAT,
	1                 STEQ(1,K-1),RUB,RUB,RUB,RUB,RUB,
	1                 N,NSNG,ENONZERO,ANONZERO)
	    END IF
!
! Computes phi[k] (Stored in L[k])
!
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	    END DO
	    CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ(1,K),N,IFAIL)
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_4'
	      GOTO 9999
	    END IF
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	    END DO
!
	  END DO
!
!**************************************************************************
!
! Now we start the back substitution.
!
	  DO K=ND-1,1,-1
!
	    DEPTH_INDX=K
	    CALL READ_BCD_MAT(RUB,RUB,D_MAT,RUB,RUB,RUB_VEC,
	1          ORIG_POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'D')
	    CALL DGEMV(NO_TRANS,N,N,DP_NEG_ONE,D_MAT,N,STEQ(1,K+1),
	1                 INT_ONE,DP_ONE,STEQ(1,K),INT_ONE)
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGEMV_2'
	      GOTO 9999
	    END IF
!
	  END DO
!
! Correct for the use of old POPS populations in doing the scaling.
!
	  DO K=1,ND
	    DO I=1,N
	      POP_DIV=MAX(POPS(I,K),
	1        1.0E-40_LDP*MAX(ABS(ORIG_POPS(I,K)),1.0E-300_LDP))
	      IF(POP_DIV .NE. POPS(I,K))THEN
	        N_DIV_FLOOR=N_DIV_FLOOR+1
	        FLAG=.TRUE.
	      END IF
	      STEQ(I,K)=STEQ(I,K)*(ORIG_POPS(I,K)/POP_DIV)
	    END DO
	  END DO
!
! Successfull solution obtained.
!
	  FLAG=.TRUE.
!
	  DEALLOCATE (B_MAT,C_MAT,D_MAT)
	  DEALLOCATE (ORIG_POPS,RUB)
	  CALL WR2D_V2(STEQ_STORE,N,ND,'STEQ_ARRAY','*',L_TRUE,16)
!
	  RETURN
	END IF
!
! 
!
! Perform the TRIDIAGONAL solution. If we reach here, we need to redo
! the LU decomposition of BA.
!
	WRITE(6,*)'Beginning TRI solution in CMF_BLKBAND'
!
	ALLOCATE (B_MAT(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (C_MAT(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (D_MAT(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (RUB(N,N),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (PREV_D_MAT(N,N),STAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(LUER,*)'Error in CMF_BLKBAND'
          WRITE(LUER,*)'Unable to allocate D_MAT etc'
          WRITE(LUER,*)'STAT=',IOS
          STOP
        END IF
	B_MAT=0.0_LDP; C_MAT=0.0_LDP; D_MAT=0.0_LDP; PREV_D_MAT=0.0_LDP
!
! WR_BA_INV is false, we perform the LU decomosition using dynamic memory
! allocation only. If insufficent storage can be allocated to store D_MAT
! (in D_MAT_STORE), we output it to disk, one depth at a time, for later
! use. In this case, D_MAT will be output, independent of WR_BA_INV.
!
	WR_D_MAT=WR_PRT_INV
	IF(.NOT. WR_BA_INV .AND. .NOT. WR_PRT_INV)THEN
          ALLOCATE (D_MAT_STORE(N,N,ND-1),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'Error in CMF_BLKBAND'
            WRITE(LUER,*)'Unable to allocate D_MAT_STORE etc'
            WRITE(LUER,*)'STAT=',IOS
            WR_D_MAT=.TRUE.
          END IF
	END IF
!
! 
!
! Solve the tridiagonal equations.
!
        ANONZERO=.FALSE.
        ENONZERO=.FALSE.
	FIRST_MATRIX=.TRUE.
	LAST_MATRIX=.FALSE.
	USE_PASSED_REP=.FALSE.
!
	DO K=1,ND
!
! Map the small BA rray onto the full BA array (one depth at a time).
!
	  DEPTH_INDX=K
	  CALL TUNE(1,'TRI_GEN')
	  CALL GENERATE_FULL_MATRIX_V3(
	1         C_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1         N,ND,NION,NUM_BNDS,
	1         DIAG_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	  FIRST_MATRIX=.FALSE.
	  IF(K .NE. 1)THEN
	    IF(K .EQ. ND)LAST_MATRIX=.TRUE.
	    BAND_INDX=DIAG_INDX-1
	    CALL GENERATE_FULL_MATRIX_V3(
	1           B_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1           N,ND,NION,NUM_BNDS,
	1           BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1         FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	  END IF
	  IF(K .NE. ND)THEN
	    BAND_INDX=DIAG_INDX+1
	    CALL GENERATE_FULL_MATRIX_V3(
	1           D_MAT,STEQ(1,K),POPS,REPLACE_EQ,ZERO_STEQ,
	1           N,ND,NION,NUM_BNDS,
	1           BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1           FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	  END IF
	  CALL TUNE(2,'TRI_GEN')
!
	  STEQ_STORE(:,K)=STEQ(:,K)
C
C Computes phi[k]' (i.e. ' denotes not yet multiplied by m[k]^{-1} )
C Stored in L[k].
C
	  CALL TUNE(1,'TRI_MAT5')
	  IF(K .NE. 1)THEN
	    CALL MAT5PEN(STEQ(1,K),B_MAT,
	1                  STEQ(1,K-1),RUB,RUB,RUB,RUB,RUB,
	1                  N,NSNG,ENONZERO,ANONZERO)
C
C Computes m[k] - Stored in C[k]
C
            CALL MAT5PEN(C_MAT,B_MAT,PREV_D_MAT,
	1                 RUB,RUB,RUB,RUB,RUB,
	1                 N,N,ENONZERO,ANONZERO)
	  END IF
	  CALL TUNE(2,'TRI_MAT5')
C
C Do LU decompostion of m[k]. We first equilibrilze C_MAT.
C
	  CALL TUNE(1,'TRI_DGEEQU')
	  CALL DGEEQU(N,N,C_MAT,N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	  CALL TUNE(2,'TRI_DGEEQU')
!
	  CALL TUNE(1,'TRI_DGETRF')
	  DO J=1,N
            DO I=1,N
	      C_MAT(I,J)=C_MAT(I,J)*ROW_SF(I)*COL_SF(J)
	    END DO
	  END DO
	  CALL DGETRF(N,N,C_MAT,N,IPIVOT,IFAIL)
	  CALL TUNE(2,'TRI_DGETRF')
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGETRF_2'
	    GOTO 9999
	  END IF
C
C Computes phi[k] (Stored in L[k])
C
	  CALL TUNE(1,'TRI_DGETRS')
	  DO J=1,N
	    STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	  END DO
	  CALL DGETRS(NO_TRANS,N,NSNG,C_MAT,N,IPIVOT,STEQ(1,K),N,IFAIL)
	  CALL TUNE(2,'TRI_DGETRS')
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGETRS_4'
	    GOTO 9999
	  END IF
	  DO J=1,N
	    STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	  END DO
C
C Computes alpha[k] (Stored in D[k)).
C
	  IF (K .NE. ND)THEN
	    DO J=1,N
              DO I=1,N
	        D_MAT(I,J)=D_MAT(I,J)*ROW_SF(I)
	      END DO
	    END DO
	    CALL TUNE(1,'TRI_DGETRS')
	    CALL DGETRS(NO_TRANS,N,N,C_MAT,N,IPIVOT,D_MAT,N,IFAIL)
	    CALL TUNE(2,'TRI_DGETRS')
	    DO J=1,N
              DO I=1,N
	        D_MAT(I,J)=D_MAT(I,J)*COL_SF(I)
	      END DO
	    END DO
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_5'
	      GOTO 9999
	    END IF
!
	    IF(WR_BA_INV .OR. WR_D_MAT)THEN
	      CALL WRITE_BCD_MAT(RUB,RUB,D_MAT,RUB,RUB,RUB_VEC,
	1            POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'D')
	    ELSE
	      D_MAT_STORE(:,:,K)=D_MAT
	    END IF
	    PREV_D_MAT=D_MAT
	  END IF
!
	  IF(WR_BA_INV)THEN
	    CALL WRITE_BCD_MAT(B_MAT,C_MAT,RUB,ROW_SF,COL_SF,IPIVOT,
	1           POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'BC')
	  END IF
C
	END DO
C
!
	CALL WR2D_V2(STEQ_STORE,N,ND,'STEQ_ARRAY','*',L_TRUE,16)
! 
!**************************************************************************
!
! Now we start the back substitution.
!
	DO K=ND-1,1,-1
!
	  DEPTH_INDX=K
	  IF(WR_BA_INV .OR. WR_D_MAT)THEN
	    CALL READ_BCD_MAT(RUB,RUB,D_MAT,RUB,RUB,RUB_VEC,
	1         POPS(:,K),REPLACE_EQ,ZERO_STEQ,N,NION,DEPTH_INDX,'D')
	  ELSE
	    D_MAT=D_MAT_STORE(:,:,K)
	  END IF
	  CALL DGEMV(NO_TRANS,N,N,DP_NEG_ONE,D_MAT,N,STEQ(1,K+1),
	1                 INT_ONE,DP_ONE,STEQ(1,K),INT_ONE)
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGEMV_2'
	    GOTO 9999
	  END IF
C
	END DO
C
C Successfull solution obtained.
C
	FLAG=.TRUE.
!
	DEALLOCATE (B_MAT,C_MAT,D_MAT,PREV_D_MAT,RUB)
	IF(ALLOCATED(D_MAT_STORE))DEALLOCATE (D_MAT_STORE)
C
	RETURN
C
C Error handling routine.
C
9999	CONTINUE
	WRITE(LUER,*)'Error in LINPAC (or BLAS) routine',
	1               DESC,' in CMF_BLKBAND'
	WRITE(LUER,100)K,IFAIL
100	FORMAT(1x,'depth=',I3,10x,'IFAIL=',I3)
	FLAG=.FALSE.
	RETURN
C
	END
