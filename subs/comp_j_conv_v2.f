C
C Routine to compute the convolution of J with the electron scattering
C redistribution function.
C
C Method uses a "Power law [in v] approximation" for the electron scattering
C function (Exponential in ln v].
C
C Basic method used is due to Rybcki and Hummer (A\&A,  ). This has
C been modified so that electron scattering preserves the Planck-function at
C depth. This scaling is necessary to recover LTE depth, and arises since
C the RH formalism neglects both Compton scattering, and stimulated electron
C scattering. The preservation of the Planck function is achieved by a
C wavelength shift on the Wien side, and a scaling everywhere else.
C
	SUBROUTINE COMP_J_CONV_V2(J_STORE,J_SIZE,NU,TEMP,ND,NCF,
	1             LU_IN,FILE_IN,CONT_REC,RD_NU,ALLOW_UNEQUAL_FREQ,
	1             LU_OUT,FILE_OUT)
	USE SET_KIND_MODULE
	USE, INTRINSIC :: IEEE_ARITHMETIC
	IMPLICIT NONE
C
C Altered 14-Nov-2026 : Planck-ratio rescale suppressed where the shifted
C                        PLANCK_ES evaluation is not genuine (upper endpoint
C                        and nodes whose bracket contains it).  B_es/B there
C                        is not a calibration and crushed live non-thermal
C                        fields (P8.1 audit sec 8.5.1).  Not flag-gated.
C                        Same edit: PLANCK_NU match loops bounded so a
C                        target outside the PLANCK_FN range, a non-finite
C                        value, or a degenerate interval cannot run ML or I
C                        off the arrays; on failure the matched region ends
C                        (or the whole correction is bypassed when no valid
C                        node exists) and the uncorrected J_ES stands.
C                        Non-finite J input is skipped with a diagnostic.
C                        Same date, follow-up: all three fields checked
C                        finite before the LOGs (MINVAL cannot spot a
C                        NaN), every continuation and interpolation
C                        denominator is validated (ES_PLANCK_NU_CONT,
C                        ES_SHIFTED_RESCALE), and a derived PLANCK_NU
C                        that is not finite-positive fails the map.  A
C                        bad node takes the same uncorrected treatment
C                        and never contaminates its neighbours.
C                        Same date, final bounds pass: the shifted
C                        walk is now written without any reliance on
C                        short-circuit evaluation (PNU(I+1) is read
C                        only after I is proven inside the written
C                        map and the node finite), the bracket
C                        arithmetic is validated strictly before use,
C                        and the interpolation weight must lie inside
C                        [0,1] up to a documented round-off band -
C                        an out-of-range weight is extrapolation
C                        through a broken map and takes the
C                        uncorrected treatment like any bad bracket.
C Altered 20-Oct-2003 : No scaling peformed when PLANCK_FN has zero''s.
C Altered 18-Aug-2003 : FILE_DATE inserted into INFO call.
C Altered 17-Oct-2002 : Computation of BB changed to avoid floating overflow.
C Altered 16-Jun-2000 : Record length can no be obtained form INFO files.
C Altered 21-Dec-1999 : Allocation of additional memory checked.
C Altered 03-Jan-1999 : J can now be passed in J_STORE in the CALL.
C                       If passed, FILE_IN must contain
C                                        'J PASSED VIA CALL'
C                      Similarly, the convolved J can be returned in
C                      call.  To do this, FILE_OUT must contain
C                                       'RETURN J VIA CALL'
C
C Altered 03-Dec-1998 : Call altered: Changed to version V2
C                         NDEXT can now be passes as ND.
C                         Better handling of dimensions of C_STORE.
C                         ALLOW_UNEQ_FREQ inserted in call. This will
C                           allow a smaller EDDFACTOR file to be used.
C			  NT, NUM_BNDS removed from call
C
C Altered 10-Feb-1998 : Section introduced to preserve Planck Function at
C                          depth.
C Altered 05-Feb-1998 : Bug fix. Routine was erroneously using T instead of
c                          SQRT(T) to compute BETA. Effect would have been
C                          for P Cygni type models, but large at depth where
C                          T is large.
C
	INTEGER J_SIZE	!Total number of elements in J_STORE
	INTEGER ND		!Number of depth points.
	INTEGER NCF		!Number of continuum frequencies
	INTEGER LU_IN
	INTEGER LU_OUT
	INTEGER CONT_REC	!Record in file that points to beginning of Jc
	LOGICAL RD_NU
	LOGICAL ALLOW_UNEQUAL_FREQ
C
C It is assumed that BA can be passed for J_STORE. J_SIZE is the maximum
C size of the passed vector.
C
	REAL(KIND=LDP) J_STORE(ND,J_SIZE/ND)		!BA matrix in main program.
	REAL(KIND=LDP) NU(NCF)				!Frequency (10^15 Hz)
	REAL(KIND=LDP) TEMP(ND)				!Temperature (10^4 K)
C
	CHARACTER*(*) FILE_IN
	CHARACTER*(*) FILE_OUT
C
C Simultaneous equations are assumed to have the form
C
C    A(i).X(i-1) - [H(i)+A(i)+C(i)].X(i) - C(i).X(i+1) = D(i)
C
	REAL(KIND=LDP) A(NCF)
	REAL(KIND=LDP) H(NCF)
	REAL(KIND=LDP) C(NCF)
	REAL(KIND=LDP) D(NCF)
C
	REAL(KIND=LDP) A_STORE(NCF)
	REAL(KIND=LDP) C_STORE(NCF)
	REAL(KIND=LDP) J_ES(NCF)
C
	REAL(KIND=LDP) PLANCK_FN(NCF)
	REAL(KIND=LDP) PLANCK_ES(NCF)
	REAL(KIND=LDP) PLANCK_NU(NCF)
C
C Provide some extra storage in case J_STORE dimensions are not large
C enough.
C
	REAL(KIND=LDP), ALLOCATABLE :: EXTRA_J_ST(:,:)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER WARNING_LU
	EXTERNAL WARNING_LU
	LOGICAL FILE_OPEN
	LOGICAL UNEQUAL_FREQ
C
	INTEGER IONE
	PARAMETER (IONE=1)
C
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) D1,D2,DH
	INTEGER I,K,L,ML,INIT_REC,IOS
	INTEGER IREC_LEN
	INTEGER J_DIM
	INTEGER NCF_RD
	INTEGER ND_RD
	INTEGER FAIL,LUE,NPNU
	LOGICAL MAP_OK,J_OK,FLD_OK
	CHARACTER*20 FILE_DATE
C
C Constants for opacity etc [Set in CMFGEN].
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) OPLIN,EMLIN
C
C Parameters for fit to Electrons Scattering redistribution function
C (dipole form). From Rybicki and Hummer (A&A, 290,553)
C
	INTEGER NCOEF
	PARAMETER (NCOEF=2)
	REAL(KIND=LDP) ACOEF(2),BCOEF(2)
	DATA ACOEF/1.690703717290_LDP,-0.690703717290_LDP/
	DATA BCOEF/1.614249968779_LDP,2.154326524957_LDP/
C
	IF(FILE_IN .EQ. 'J PASSED VIA CALL')THEN
	  IF(J_SIZE .NE. NCF*ND)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in COMP_J_COV_V2: When J is passed in call'
	    WRITE(I,*)'J_SIZE must be the same as NCF*ND'
	    WRITE(I,*)'J_SIZE=',J_SIZE,'NCF*ND=',NCF*ND
	    STOP
	  END IF
	  J_DIM=NCF		!Must be set
	  GOTO 1000		!Jump the section reading in the data
	END IF
!
	IF(FILE_OUT .EQ. 'RETURN J VIA CALL')THEN
	  IF(J_SIZE .NE. NCF*ND)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in COMP_J_COV_V2: When J is returned in call'
	    WRITE(I,*)'J_SIZE must be the same as NCF*ND'
	    WRITE(I,*)'J_SIZE=',J_SIZE,'NCF*ND=',NCF*ND
	    STOP
	  END IF
	END IF
C
C It is assumed that we can use BA to store J. If this is not the case we
C allocate some extra storage.
C
	J_DIM=J_SIZE/ND
	IF(J_DIM .LT. NCF)THEN
	  ALLOCATE (EXTRA_J_ST(ND,J_DIM+1:NCF),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    I=ERROR_LU()
	    WRITE(I,'(A)')'Error in COMP_J_CONV_V2'
	    WRITE(I,'(A)')'Unable to allocate memeory'
	    STOP
	  END IF
	ELSE
	  J_DIM=NCF
	END IF
C
C Open file with mean intensities. Read in J for all depths at all frequencies.
C
	INQUIRE(UNIT=LU_IN,OPENED=FILE_OPEN)
        IF(FILE_OPEN)INQUIRE(FILE=FILE_IN,OPENED=FILE_OPEN)
	IF(.NOT. FILE_OPEN)THEN
	  CALL READ_DIRECT_INFO_V3(I,IREC_LEN,FILE_DATE,FILE_IN,LU_IN,IOS)
	  IF(IOS .NE. 0)STOP
	  OPEN(UNIT=LU_IN,FILE=FILE_IN,STATUS='OLD',RECL=IREC_LEN,
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
	END IF
C
	INQUIRE(UNIT=LU_IN,RECL=IREC_LEN)
	READ(LU_IN,REC=CONT_REC)INIT_REC,NCF_RD,ND_RD
	IF(ND_RD .NE. ND)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in COMP_J_COV_V2'
	  WRITE(I,*)'ND and ND_RD must agree'
	  WRITE(I,*)'ND=',ND,'ND_RD=',ND_RD
	  STOP
	END IF
	IF(NCF_RD .GT. NCF)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in COMP_J_COV_V2'
	  WRITE(I,*)'NCF must be greater than NCF_RD'
	  WRITE(I,*)'NCF=',NCF,'NCF_RD=',NCF_RD
	  STOP
	END IF
C
	IF(RD_NU)THEN
	  DO ML=1,J_DIM
	    READ(LU_IN,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),NU(ML)
	  END DO
	  DO ML=J_DIM+1,NCF
	    READ(LU_IN,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),NU(ML)
	  END DO
	ELSE IF(ALLOW_UNEQUAL_FREQ)THEN
	  UNEQUAL_FREQ=.FALSE.
	  DO ML=1,NCF_RD
	    IF(ML .LE. J_DIM)THEN
	      READ(LU_IN,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),A(ML)
	    ELSE
	      READ(LU_IN,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),A(ML)
	    END IF
	    IF(A(ML) .NE. NU(ML))UNEQUAL_FREQ=.TRUE.
	  END DO
	  IF(NCF .NE. NCF_RD)UNEQUAL_FREQ=.TRUE.
	  IF(UNEQUAL_FREQ)THEN
	    IF(NU(1) .GT. A(1) .OR. NU(NCF) .LT. A(NCF_RD))THEN
	      I=ERROR_LU()
	      WRITE(I,*)'Insufficient range of frequencies in EDFFACTOR'
	      WRITE(I,*)'Routine is COMP_J_CONV'
	      STOP
	    END IF
	    DO K=1,ND
	      H(1:J_DIM)=J_STORE(K,1:J_DIM)
	      IF(NCF_RD .GT. J_DIM)H(J_DIM+1:NCF_RD)=J_STORE(K,J_DIM+1:NCF_RD)
	      CALL MON_INTERP(C,NCF,IONE,NU,NCF,H,NCF_RD,A,NCF_RD)
	      J_STORE(K,1:J_DIM)=C(1:J_DIM)
	      IF(NCF .GT. J_DIM)EXTRA_J_ST(K,J_DIM+1:NCF)=C(J_DIM+1:NCF)
	    END DO
	  END IF
	ELSE
	  DO ML=1,NCF
	    IF(ML .LE. J_DIM)THEN
	      READ(LU_IN,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),T1
	    ELSE
	      READ(LU_IN,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),T1
	    END IF
	    IF(T1 .NE. NU(ML))THEN
	      I=ERROR_LU()
	      WRITE(I,*)'Invalid frequency in ',FILE_OUT
	      WRITE(I,*)'Routine is COMP_J_CONV'
	      STOP
	    END IF
	  END DO
	END IF
	CLOSE(UNIT=LU_IN)
C
C*****************************************************************************
C*****************************************************************************
C
C All data has now been read in. We can now perform the convolutions.
c
1000	CONTINUE
C
C Compute those parts of the TRIDIAGONAL vectors which are independent of
C depth, and the fitting parameters.
C
	A_STORE(1)=0.0_LDP
	C_STORE(1)=-2.0_LDP/( LOG(NU(1)/NU(2)) )**2
	DO ML=2,NCF-1
	  D1=LOG(NU(ML-1)/NU(ML))
	  D2=LOG(NU(ML)/NU(ML+1))
	  DH=0.5_LDP*(D1+D2)
	  A_STORE(ML)=-1.0_LDP/D1/DH
	  C_STORE(ML)=-1.0_LDP/D2/DH
	END DO
	A_STORE(NCF)=-2.0_LDP/( LOG(NU(NCF-1)/NU(NCF)) )**2
	C_STORE(NCF)=0.0_LDP
C
C Compute the triadiagonal quantities for performing the convolution, and
C perform the convolution. Due to the depth dependence of BETA, the vectors
C are depth dependent. We convolve both J and the Planck Function.
C
!
! We can scale the Placnk function to avoid overflow.
!
	DO K=1,ND
C
C A non-finite input field cannot be convolved or logged - decline the
C depth and leave the invalid input visible downstream rather than
C writing solver garbage over it.  The EXTRA_J_ST check is nested so
C the unallocated array is never referenced.
C
	  J_OK=ALL(IEEE_IS_FINITE(J_STORE(K,1:J_DIM)))
	  IF(J_OK .AND. J_DIM .LT. NCF)THEN
	    J_OK=ALL(IEEE_IS_FINITE(EXTRA_J_ST(K,J_DIM+1:NCF)))
	  END IF
	  IF(.NOT. J_OK)THEN
	    LUE=WARNING_LU()
	    WRITE(LUE,'(A,I5,A)')'*** Warning in COMP_J_CONV_V2: '//
     1	    'non-finite J input at depth',K,' - convolution skipped'
	    CYCLE
	  END IF
C
	  BETA=1.84E-03_LDP*SQRT(TEMP(K))
	  T3=0.5_LDP*HDKT*NU(1)/TEMP(K)
	  IF(T3 .LT. 1)T3=0.0_LDP
	  IF(T3 .GT. 700.0_LDP)T3=700.0_LDP
	  DO ML=1,NCF
	    T1=EXP( -HDKT*NU(ML)/TEMP(K) )
	    T2=EXP( T3-HDKT*NU(ML)/TEMP(K) )
	    PLANCK_FN(ML)=TWOHCSQ*(NU(ML)**3)*T2/ (1.0_LDP-T1)
	  END DO
C
	  J_ES(:)=0.0_LDP
	  PLANCK_ES(:)=0.0_LDP
C
	  DO L=1,NCOEF
	    T1=BETA*BETA/BCOEF(L)/BCOEF(L)
	    A(:)=T1*A_STORE(:)			!Over frequency
	    H(:)=-1.0_LDP
	    C(:)=T1*C_STORE(:)
	    D(1:J_DIM)=J_STORE(K,1:J_DIM)
	    IF(J_DIM .LT. NCF)D(J_DIM+1:NCF)=EXTRA_J_ST(K,J_DIM+1:NCF)
	    CALL THOMAS_RH(A,H,C,D,NCF,IONE)
	    J_ES(:)=J_ES(:)+ACOEF(L)*D(:)
	    D(:)=PLANCK_FN(:)
	    CALL SIMPTH_RH(A,H,C,D,NCF,IONE)
	    PLANCK_ES(:)=PLANCK_ES(:)+ACOEF(L)*D(:)
	  END DO				!Fit parameter loop.
C
C Derive the wavelength shift so that the electron scattered planck function
C gives the Planck function. This is only done on the Wien side of the
C BB curve. At longer wavelengths we will use a simple scaling.
C
C Because of the steep variation of B on the Wien side, we operate on Log(B).
C
C MINVAL alone cannot be trusted to spot a non-finite field, so the
C finiteness of all three logarithm arguments is checked explicitly
C before any LOG is taken.  A non-finite convolver output (J_ES,
C PLANCK_ES or PLANCK_FN) is diagnosed per field and bypassed, and
C the invalid result is left in J_ES for the caller to see - the same
C policy as the non-finite input check above, not a silent repair.
C Finite but non-positive values still take the stock MINVAL gate:
C that branch is an expected physical state (a zero Planck node on
C the coldest depths), not a numerical defect.
C
	  FLD_OK=ALL(IEEE_IS_FINITE(PLANCK_FN)) .AND.
     1	       ALL(IEEE_IS_FINITE(PLANCK_ES)) .AND.
     2	       ALL(IEEE_IS_FINITE(J_ES))
	  IF(.NOT. FLD_OK)THEN
	    LUE=WARNING_LU()
	    WRITE(LUE,'(A,I5,A)')'*** Warning in COMP_J_CONV_V2: '//
     1	    'non-finite convolver output at depth',K,
     2	    ' - Planck correction bypassed'
	    IF(.NOT. ALL(IEEE_IS_FINITE(J_ES)))WRITE(LUE,'(A)')
     1	    '    J_ES (convolver output) is non-finite'
	    IF(.NOT. ALL(IEEE_IS_FINITE(PLANCK_FN)))WRITE(LUE,'(A)')
     1	    '    PLANCK_FN is non-finite'
	    IF(.NOT. ALL(IEEE_IS_FINITE(PLANCK_ES)))WRITE(LUE,'(A)')
     1	    '    PLANCK_ES is non-finite'
	  END IF
	  IF( FLD_OK .AND. MINVAL(PLANCK_FN) .GT. 0.0_LDP .AND.
     1	    MINVAL(PLANCK_ES) .GT. 0.0_LDP .AND.
     2	    MINVAL(J_ES) .GT. 0.0_LDP)THEN
	    PLANCK_ES(:)=LOG(PLANCK_ES(:))
	    PLANCK_FN(:)=LOG(PLANCK_FN(:))
C
	    PLANCK_NU(1)=NU(1)
	    I=2
	    ML=1
	    MAP_OK=.TRUE.
	    DO WHILE(I .LE. NCF)
	      IF(NU(I) .LE. 1.5_LDP*TEMP(K))EXIT
C
C The match itself lives in ES_PLANCK_NU_MATCH (end of this file):
C a bounded form of the stock two-phase walk.  NU decreases with
C index; the NU(I) > 1.5*TEMP gate keeps the entered region on the
C Wien side of the Planck maximum (nu_peak ~ 0.59*TEMP), where
C PLANCK_FN increases with ML, so the intended bracket is the first
C ML in [2,NCF] with PLANCK_FN(ML-1) <= PLANCK_ES(I) <= PLANCK_FN(ML).
C The stock loops walked ML with no bound: a target outside the
C PLANCK_FN range ran off the end of the array, and the ML-1 lookup
C at ML=1 read PLANCK_FN(0).  The subroutine keeps every PLANCK_FN
C access in bounds without relying on short-circuit evaluation.
C
	      CALL ES_PLANCK_NU_MATCH(PLANCK_ES(I),PLANCK_FN,ML,NCF,
     1	                            FAIL,T1)
	      IF(FAIL .NE. 0)THEN
C
C No genuine bracket: a target outside the PLANCK_FN range has no
C on-grid shifted position, and extrapolating would invent a
C PLANCK_NU that violates the ln B_es(NU) = ln B(PLANCK_NU) relation.
C With at least one derived node the matched region simply ends here
C and the continuation below supplies the last valid shift, the same
C treatment the stock code applies below the match region.  A
C first-node failure means no valid map exists at all, so the whole
C Planck correction is bypassed and the raw convolved J_ES stands
C (the uncorrected treatment of the 14-Nov-2026 edit).
C
	        LUE=WARNING_LU()
	        WRITE(LUE,'(A,I5,A,I5,A,I2)')'*** Warning in '//
     1	        'COMP_J_CONV_V2: PLANCK_NU match failed, K=',K,
     2	        ' I=',I,' FAIL=',FAIL
	        WRITE(LUE,'(A,4ES15.7)')'    NU(I), PLANCK_ES(I), '//
     1	        'PLANCK_FN(1), PLANCK_FN(NCF)=',NU(I),PLANCK_ES(I),
     2        PLANCK_FN(1),PLANCK_FN(NCF)
	        IF(I .GT. 2)THEN
	          WRITE(LUE,'(A)')'    matched region truncated; '//
     1	          'last valid shift continues the map'
	        ELSE
	          WRITE(LUE,'(A)')'    no valid node: Planck '//
     1	          'correction bypassed, uncorrected J_ES returned'
	          MAP_OK=.FALSE.
	        END IF
	        EXIT
	      END IF
              PLANCK_NU(I)=T1*NU(ML-1)+(1.0_LDP-T1)*NU(ML)
C
C A matched weight guarantees PLANCK_NU between NU(ML) and NU(ML-1)
C only while NU itself is positive and finite.  The map must stay a
C positive frequency - anything else would poison the continuation
C anchors and the interpolation below.
C
	      IF(.NOT. IEEE_IS_FINITE(PLANCK_NU(I)) .OR.
     1	         PLANCK_NU(I) .LE. 0.0_LDP)THEN
	        FAIL=7
	        LUE=WARNING_LU()
	        WRITE(LUE,'(A,I5,A,I5,A,I2)')'*** Warning in '//
     1	        'COMP_J_CONV_V2: PLANCK_NU match failed, K=',K,
     2	        ' I=',I,' FAIL=',FAIL
	        WRITE(LUE,'(A,4ES15.7)')'    NU(I), PLANCK_ES(I), '//
     1	        'PLANCK_FN(1), PLANCK_FN(NCF)=',NU(I),PLANCK_ES(I),
     2        PLANCK_FN(1),PLANCK_FN(NCF)
	        IF(I .GT. 2)THEN
	          WRITE(LUE,'(A)')'    matched region truncated; '//
     1	          'last valid shift continues the map'
	        ELSE
	          WRITE(LUE,'(A)')'    no valid node: Planck '//
     1	          'correction bypassed, uncorrected J_ES returned'
	          MAP_OK=.FALSE.
	        END IF
	        EXIT
	      END IF
	      I=I+1
	    END DO
C
	    IF(MAP_OK)THEN
C
C The matched region ends at I-1; ES_PLANCK_NU_CONT continues the
C map with the same wavelength shift as the last derived node, the
C stock treatment below the matched region.  Every reciprocal is
C guarded inside the subroutine: a bad anchor or a non-positive,
C non-finite or vanishing denominator ends the map at NPNU rather
C than producing an infinite or negative frequency.  Only the
C nodes in 1:NPNU are genuine shifted positions.
C
	    CALL ES_PLANCK_NU_CONT(NU,PLANCK_NU,I,NCF,NPNU,FAIL)
	    IF(FAIL .NE. 0)THEN
	      LUE=WARNING_LU()
	      WRITE(LUE,'(A,I5,A,I5,A,I2)')'*** Warning in '//
     1	      'COMP_J_CONV_V2: PLANCK_NU continuation failed, K=',
     2	      K,' at ML=',NPNU+1,' FAIL=',FAIL
	      WRITE(LUE,'(A,I5,A)')'    map truncated after node',
     1	      NPNU,': below-map nodes take the uncorrected treatment'
	    END IF
C
C ES_SHIFTED_RESCALE interpolates the scattered Planck function and
C J_ES back onto the frequency grid through the shifted map, then
C applies the B/B_es rescale.  Every interpolation bracket is
C validated inside: a truncated, non-finite, non-positive or
C degenerate bracket leaves that node in the established
C uncorrected state instead of dividing through a bad interval.
C
	    CALL ES_SHIFTED_RESCALE(NU,PLANCK_NU,PLANCK_FN,PLANCK_ES,
     1	    J_ES,A,C,NCF,NPNU,K)
C
	    END IF
C
C When MAP_OK is false no valid shift map exists: J_ES still holds the
C raw convolution, which is what the correction would have reduced to
C anyway where no genuine shifted PLANCK_ES exists.
C
	END IF
C
C Store the computed  J to be used to treat the electron scattering. We
C overwrite J_STORE since the mean intensity is no longer required in this
C routine.
C
	  J_STORE(K,1:J_DIM)=J_ES(1:J_DIM)
	  IF(J_DIM .LT. NCF)EXTRA_J_ST(K,J_DIM+1:NCF)=J_ES(J_DIM+1:NCF)
C
	END DO					!Depth
C
C 
C ****************************************************************************
C ****************************************************************************
C
	IF(FILE_OUT .NE. 'RETURN J VIA CALL')THEN
C
C Now output J convolution to data file. Data file will have exactly the same
C format as EDDFACTOR file (file that contains RJ).
C
	  INQUIRE(UNIT=LU_OUT,OPENED=FILE_OPEN)
          IF(FILE_OPEN)INQUIRE(FILE=FILE_OUT,OPENED=FILE_OPEN)
	  IF(.NOT. FILE_OPEN)THEN
	    CALL WRITE_DIRECT_INFO_V3(ND,IREC_LEN,'20-Aug-2000',FILE_OUT,LU_OUT)
	    OPEN(UNIT=LU_OUT,FILE=FILE_OUT,FORM='UNFORMATTED',
	1      ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IREC_LEN)
	  END IF
	  WRITE(LU_OUT,REC=CONT_REC)INIT_REC,NCF,ND
	  DO ML=1,J_DIM
	    WRITE(LU_OUT,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),NU(ML)
	  END DO
	  DO ML=J_DIM+1,NCF
	    WRITE(LU_OUT,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),NU(ML)
	  END DO
	  CLOSE(UNIT=LU_OUT)
C
	  J_STORE(:,1:J_DIM)=0.0_LDP
	  IF(J_DIM .NE. NCF)DEALLOCATE (EXTRA_J_ST)
	END IF
C
	RETURN
	END
C
C
C
	SUBROUTINE ES_PLANCK_NU_MATCH(PES,PFN,ML,NCF,FAIL,T1)
	USE SET_KIND_MODULE
	USE, INTRINSIC :: IEEE_ARITHMETIC
	IMPLICIT NONE
C
C Bounded bracket search for the COMP_J_CONV_V2 Planck wavelength
C shift, split out 14-Nov-2026 so the same production code can be
C exercised directly by the manufactured driver.  PES is a log-space
C target PLANCK_ES(I); PFN is the full log-B grid.  ML carries the
C walk position between calls and only ever moves upward, as in the
C stock two-loop search.
C
C The caller gates on NU(I) > 1.5*TEMP, which keeps the entered
C region above the Planck maximum (nu_peak ~ 0.59*TEMP); there PFN
C increases with index and the intended match is the first ML in
C [2,NCF] with PFN(ML-1) <= PES <= PFN(ML).  Because NU decreases
C with index, T1 = (PES-PFN(ML))/(PFN(ML-1)-PFN(ML)) gives the
C shifted position T1*NU(ML-1) + (1-T1)*NU(ML).
C
C FAIL codes:
C   0  bracket found; T1 set.
C   1  PES above every PFN(ML): the smeared target is larger than
C      any on-grid B, so no shifted position exists.
C   2  PES below the reachable floor: either PES <= PFN(1) (the
C      stock code read PFN(0) here) or PFN(ML-1) still exceeds PES
C      when ML has walked off the top of the array.
C   3  degenerate interval PFN(ML-1) = PFN(ML): the denominator
C      is zero or non-finite and no weight can be formed.
C   4  PES is not finite.
C In every failure case no PLANCK_NU is produced: inventing one
C would violate the ln B_es(NU) = ln B(PLANCK_NU) relation, so the
C caller falls back to the established uncorrected-convolution
C treatment rather than extrapolating.
C
	INTEGER NCF,ML,FAIL
	REAL(LDP) PES,T1
	REAL(LDP) PFN(NCF)
	REAL(LDP) T2
C
	FAIL=0
	T1=0.0_LDP
C
	IF(.NOT. IEEE_IS_FINITE(PES))THEN
	  FAIL=4
	  RETURN
	END IF
C
	DO WHILE(ML .LE. NCF)
	  IF(PES .LE. PFN(ML))EXIT
	  ML=ML+1
	END DO
	IF(ML .GT. NCF)THEN
	  FAIL=1
	  RETURN
	END IF
	IF(ML .EQ. 1)THEN
	  FAIL=2
	  RETURN
	END IF
C
	DO WHILE(ML .LE. NCF)
	  IF(PES .GE. PFN(ML-1))EXIT
	  ML=ML+1
	END DO
	IF(ML .GT. NCF)THEN
	  FAIL=2
	  RETURN
	END IF
C
	T2=PFN(ML-1)-PFN(ML)
	IF(T2 .EQ. 0.0_LDP .OR. .NOT. IEEE_IS_FINITE(T2))THEN
	  FAIL=3
	  RETURN
	END IF
C
	T1=(PES-PFN(ML))/T2
	RETURN
	END
C
C
C
	SUBROUTINE ES_PLANCK_NU_CONT(NU,PNU,I,NCF,NPNU,FAIL)
	USE SET_KIND_MODULE
	USE, INTRINSIC :: IEEE_ARITHMETIC
	IMPLICIT NONE
C
C Wavelength-shift continuation for the COMP_J_CONV_V2 Planck map,
C split out 14-Nov-2026 with the same purpose as
C ES_PLANCK_NU_MATCH: the same production code is driven directly
C by the manufactured driver.  PNU(1:I-1) already holds genuine
C derived nodes (endpoint plus match results); nodes I:NCF take the
C constant wavelength shift of the last derived node, exactly the
C stock continuation:
C
C   T1 = 1/NU(I-1) - 1/PNU(I-1)
C   PNU(ML) = 1/(1/NU(ML) - T1)
C
C Every reciprocal is guarded: the anchors and each denominator
C must be finite, and the denominator must be strictly positive so
C the continued frequency stays a positive physical frequency.
C With NU strictly decreasing and a positive PNU(I-1) anchor the
C denominator is > 1/PNU(I-1) by construction, so the denominator
C failure is unreachable through production input - it is a
C defensive bound exercised directly by the driver.
C
C On failure no extrapolated shift is invented: the map keeps only
C its genuine nodes and NPNU reports the last valid index.  Nodes
C beyond NPNU are left untouched (never written) and the caller's
C interpolation treats the region below PNU(NPNU) as off-map.
C
C FAIL codes (continuing the ES_PLANCK_NU_MATCH sequence):
C   5  anchor NU(I-1) or PNU(I-1) not finite-positive, or the
C      shift T1 they produce is non-finite; NPNU=I-1.
C   6  denominator at ML not finite or not positive (covers a
C      non-positive or non-finite NU(ML)); NPNU=ML-1.
C
	INTEGER NCF,I,NPNU,FAIL
	REAL(LDP) NU(NCF),PNU(NCF)
	REAL(LDP) T1,DEN
	INTEGER ML
C
	FAIL=0
	NPNU=NCF
C
C Sequential checks again: a NaN operand inside a chained .OR.
C compare is an unordered comparison and can trip an
C invalid-operation trap, so each value is proven finite before
C it is compared.
C
	IF(FAIL .EQ. 0)THEN
	  IF(.NOT. IEEE_IS_FINITE(NU(I-1)))FAIL=5
	END IF
	IF(FAIL .EQ. 0)THEN
	  IF(NU(I-1) .LE. 0.0_LDP)FAIL=5
	END IF
	IF(FAIL .EQ. 0)THEN
	  IF(.NOT. IEEE_IS_FINITE(PNU(I-1)))FAIL=5
	END IF
	IF(FAIL .EQ. 0)THEN
	  IF(PNU(I-1) .LE. 0.0_LDP)FAIL=5
	END IF
	IF(FAIL .NE. 0)THEN
	  NPNU=I-1
	  RETURN
	END IF
	T1=1.0_LDP/NU(I-1)-1.0_LDP/PNU(I-1)
	IF(.NOT. IEEE_IS_FINITE(T1))THEN
	  FAIL=5
	  NPNU=I-1
	  RETURN
	END IF
C
	DO ML=I,NCF
	  IF(.NOT. IEEE_IS_FINITE(NU(ML)))THEN
	    FAIL=6
	    NPNU=ML-1
	    RETURN
	  END IF
	  IF(NU(ML) .LE. 0.0_LDP)THEN
	    FAIL=6
	    NPNU=ML-1
	    RETURN
	  END IF
	  DEN=1.0_LDP/NU(ML)-T1
	  IF(.NOT. IEEE_IS_FINITE(DEN))THEN
	    FAIL=6
	    NPNU=ML-1
	    RETURN
	  END IF
	  IF(DEN .LE. 0.0_LDP)THEN
	    FAIL=6
	    NPNU=ML-1
	    RETURN
	  END IF
	  PNU(ML)=1.0_LDP/DEN
	END DO
	RETURN
	END
C
C
C
	SUBROUTINE ES_SHIFTED_RESCALE(NU,PNU,PFN,PES,JES,A,C,NCF,
     1	                            NPNU,K)
	USE SET_KIND_MODULE
	USE, INTRINSIC :: IEEE_ARITHMETIC
	IMPLICIT NONE
C
C Shifted-grid interpolation and the Planck-preserving rescale of
C COMP_J_CONV_V2, split out 14-Nov-2026 so the manufactured driver
C can feed synthetic maps through the same object code.  PNU(1:NPNU)
C is the validated shifted-frequency map from ES_PLANCK_NU_MATCH /
C ES_PLANCK_NU_CONT; PFN, PES and JES are the log-space Planck,
C scattered Planck and convolved fields on the NU grid.
C
C For each interior output node the bracket
C PNU(I+1) <= NU(ML) <= PNU(I) is located by the stock walk, then
C checked before it is used: both endpoints must be finite positive
C frequencies, the interval finite and nonzero, the weight finite
C and inside [0,1] - a genuine shifted evaluation is an
C interpolation, never an extrapolation through a broken map.
C A mathematically monotone map can still carry duplicate,
C degenerate or locally non-monotone nodes numerically.
C
C Where no genuine shifted evaluation exists the node takes the
C established invalid-map treatment: PLANCK_ES(ML)=PLANCK_FN(ML)
C (unity rescale) and J_ES(ML)=A(ML) (raw convolved field).  That
C covers the I=1 bracket containing the forced endpoint
C (P8.1 Wien-end correction), an output node below the last valid
C map frequency PNU(NPNU) (map truncated by
C ES_PLANCK_NU_CONT), any degenerate bracket, and any out-of-range
C weight.  Nodes whose bracket is valid keep their interpolated
C values - an invalid node never contaminates a valid one.
C
	INTEGER NCF,NPNU,K
	REAL(LDP) NU(NCF),PNU(NCF),PFN(NCF),PES(NCF),JES(NCF)
	REAL(LDP) A(NCF),C(NCF)
C
C Round-off band on the interpolation weight: the bracket
C arithmetic can round a genuinely valid node a few ulp past the
C endpoint, so W_TOL ulps are snapped back to [0,1].  Anything
C beyond the band is a real extrapolation and is left invalid.
C
	REAL(LDP), PARAMETER :: W_TOL=64.0_LDP*EPSILON(1.0_LDP)
C
	INTEGER ML,I,NBAD,LUE
	REAL(LDP) T1,T2
	LOGICAL BRK_OK,NU_OK
	INTEGER WARNING_LU
	EXTERNAL WARNING_LU
C
	A(:)=LOG(JES(:))
	JES(1)=A(1)		!Must be in LOG form: not modified by interp.
	JES(NCF)=A(NCF)
	C(:)=PES(:)		!Already taken LOG
	I=1
	NBAD=0
	DO ML=2,NCF-1
C
C Locate the bracket by the stock upward walk, but bound it in
C separate control flow: Fortran need not short-circuit .AND.,
C so the bound on I is tested in its own IF before PNU(I+1) is
C touched, and PNU(I+1) is checked finite before the compare -
C an unordered NaN compare would stop the walk anyway and could
C trip an invalid-operation trap.  For NPNU < 2 the walk body
C is skipped entirely and PNU(2) is never read.
C
	  NU_OK=IEEE_IS_FINITE(NU(ML)) .AND. NU(ML) .GT. 0.0_LDP
	  IF(NU_OK .AND. NPNU .GE. 2)THEN
	    DO
	      IF(I .GE. NPNU-1)EXIT
	      IF(.NOT. IEEE_IS_FINITE(PNU(I+1)))EXIT
	      IF(NU(ML) .LT. PNU(I+1))THEN
	        I=I+1
	      ELSE
	        EXIT
	      END IF
	    END DO
	  END IF
C
C A truncated map ends at PNU(NPNU): an output node below that
C frequency has no on-grid bracket at all.  The stock clamp would
C only extrapolate inside the last bracket, so the node takes the
C uncorrected treatment instead.  When NPNU=NCF the bound
C reproduces the stock NCF-1 cap exactly.
C
	  BRK_OK=NU_OK .AND. NPNU .GE. 2
	  IF(BRK_OK .AND. NPNU .LT. NCF)THEN
	    IF(I .EQ. NPNU-1)THEN
	      IF(IEEE_IS_FINITE(PNU(NPNU)))THEN
	        IF(NU(ML) .LT. PNU(NPNU))BRK_OK=.FALSE.
	      END IF
	    END IF
	  END IF
	  IF(BRK_OK)BRK_OK=I .GE. 1 .AND. I+1 .LE. NPNU
	  IF(BRK_OK)BRK_OK=IEEE_IS_FINITE(PNU(I))
	  IF(BRK_OK)BRK_OK=PNU(I) .GT. 0.0_LDP
	  IF(BRK_OK)BRK_OK=IEEE_IS_FINITE(PNU(I+1))
	  IF(BRK_OK)BRK_OK=PNU(I+1) .GT. 0.0_LDP
	  IF(BRK_OK)THEN
	    T2=PNU(I+1)-PNU(I)
	    BRK_OK=IEEE_IS_FINITE(T2) .AND. T2 .NE. 0.0_LDP
	  END IF
	  IF(BRK_OK)THEN
	    T1=(NU(ML)-PNU(I))/T2
	    BRK_OK=IEEE_IS_FINITE(T1)
	  END IF
	  IF(BRK_OK)BRK_OK=T1 .GE. -W_TOL .AND.
     1	                 T1 .LE. 1.0_LDP+W_TOL
	  IF(BRK_OK)THEN
	    T1=MAX(0.0_LDP,MIN(1.0_LDP,T1))
	    PES(ML)=T1*C(I+1)+(1.0_LDP-T1)*C(I)
	    JES(ML)=T1*A(I+1)+(1.0_LDP-T1)*A(I)
	  END IF
C
C The I=1 bracket containing the forced PLANCK_NU(1)=NU(1) anchor
C and every invalid bracket share the same treatment: no genuine
C shifted evaluation exists, so the B/B_es factor is unity and the
C raw convolved J_ES stands (the 14-Nov-2026 Wien-end correction).
C
	  IF(.NOT. BRK_OK .OR. I .EQ. 1)THEN
	    IF(.NOT. BRK_OK)NBAD=NBAD+1
	    PES(ML)=PFN(ML)
	    JES(ML)=A(ML)
	  END IF
	END DO
C
	IF(NBAD .GT. 0)THEN
	  LUE=WARNING_LU()
	  WRITE(LUE,'(A,I5,A,I5,A)')'*** Warning in '//
     1	  'COMP_J_CONV_V2: ',NBAD,' node(s) at depth',K,
     2	  ' have no genuine shifted bracket - uncorrected J_ES kept'
	END IF
C
C The same invalid-map treatment applies at the upper endpoint:
C its shifted position lies off the grid (the smeared B_es is
C smallest at the boundary itself), so no genuine shifted
C PLANCK_ES exists and the rescale is set to unity.
C
	PES(1)=PFN(1)
C
C To remove any residual variations (primarily at low frequencies) we now
C scale J according the the ratio of the PLANCK functions.
C
	JES(:)=EXP( JES(:)+ (PFN(:)-PES(:)) )
	RETURN
	END
