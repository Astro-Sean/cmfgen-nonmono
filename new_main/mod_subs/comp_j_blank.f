!
! Subroutine to compute the mean intensity. This routine replaces
! COMP_J_CONT.INC and was developed to faciltate inclusion of additional
! options without making CMFGEN_SUB continually larger.
!
! Available options for computing J: These should be the same as for DO_VAR_CONT.
! NB: DO_VAR CONT depends on results from this routine, however this routine is
! independent of DO_VAR_CONT.
!
!       CONT_VEL=.TRUE.			!Velocity field taken into account.
!            ACURATE=.TRUE.		!Enhanced spatial grid
!            ACCURATE=.FALSE.		!Regular spatial grid
!
!If the following options hold, we use Eddington factors and the enhanced grid.
!
!       CONT_VEL=.FALSE.  & THIS_FREQ_EXT=.TRUE.
!
!If the following options hold, we use Eddington factors and the regular grid.
!
!	CONT_VEL=.FALSE., EDDINGTON=.TRUE. & THIS_FREQ_EXT=.FALSE.
!
!If the following options hold we use ray-by ray solution for the full computaion.
!Spherical model, and no velocity terms.
!
!	CONT_VEL=.FALSE., EDDINGTON=.FALSE. & THIS_FREQ_EXT=.FALSE.
!
	SUBROUTINE COMP_J_BLANK(SECTION,EDDINGTON,FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,MAXCH,
	1                              LUER,LU_ES,LU_JCOMP,LU_EDD,ACCESS_F,
	1                              ND,NC,NP,NCF,NDEXT,NCEXT,NPEXT)
	USE SET_KIND_MODULE
!
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE RADIATION_MOD
	USE CONTROL_VARIABLE_MOD
	USE MOD_RAY_MOM_STORE, ONLY: HNU_AT_OB
	IMPLICIT NONE
!
! Altered : 05-Oct-2026 : Keep stock CMF_FORM_SOL off the Knop path (Devin).
! Altered : 05-Oct-2026 : Demoted solver banners to iteration subheaders (edited by Devin).
! Altered : 05-Oct-2026 : JFNK eval hygiene: softened STOPs inside
!                            JFNK_IN_EVAL eval sweeps; cold FOLD seed.
! Altered : 19-Sep-2026 : dev7 Knop integration (DO_NONMONO_KNOP): dispatch to
!                            FG_J_CMF_NONMONO_V1 (USE_LAM_ES=T) or MOM_J_CMF_NONMONO_V1
!                            (USE_LAM_ES=F), one-sided dLOG_NU at the bluest frequency,
!                            GET_MOMS_NON_REL Eddington fill, and the upstream
!                            'ZERO FLUX' -> 'ZERO_FLUX' comparison fix.
! Altered : 09-Jun-2019 : Added POINT srce options - osiris (added to IBIS 17-Aug-2019)
! Altered : 30-Apr-2019 : Changed to MOM_JREL_V9 (add XM_CHK_OPTION and J_CHK_OPTION).
! Altered : 29-Apr-2019 : XM_CHK_OPTION added. Changed to MOM_J_DDT_V6.
! Altered : 20-Apr-2019 : J_CHK_OPTION added. Changed to MOM_J_DDT_V5.
! Altered : 17-Oct-2016 : H_CHK_OPTION added to moment routines.
!                            Moment routines are: MOM_J_CMF_V11.F, MOM_J_DDT_V4.F, and MOM_JREL_V8.F.
! Altered : 17-Feb-2015 : r^2.J and r^2.H now ouput on last iteration when using USE_LAM_ES option.
! Altered : 14-Dec-2014 : RSQHNU etc now set to one when not computing J. This avoids issues with
!                             possible NaNs.
! Altered : 16-Dec-2013 : CMF_FORM_SOL_V2 (non EXT option) no longer called when ND > 199.
!                             CMF_FORM_SOL_V2 is not parallelized and slows down large clumped models.
! Altered : 16-Feb-2006 : CMF_FORM_SOL_V2 used for last iteration when MAXCH<100, and
!                            not LAMBDA iteration. Sometimes it might be useful to
!                            change so that CMF_FORM_SOL_V2 is also called when LMABDA
!                            iteration used. FG_COUNt was not being initialized.
! Finalized: 17-Dec-2004
!
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	REAL(KIND=LDP) FL
	REAL(KIND=LDP) MAXCH
!
	INTEGER ACCESS_F
	INTEGER LU_EDD
	INTEGER LUER
	INTEGER LU_ES
	INTEGER LU_JCOMP
	INTEGER ND,NC,NP
	INTEGER NDEXT,NCEXT,NPEXT
	INTEGER NCF
	INTEGER FREQ_INDX
	CHARACTER(LEN=*) SECTION
!
	LOGICAL LST_ITERATION
	LOGICAL FIRST_FREQ
	LOGICAL EDDINGTON
!
! Local variables and arrays.
!
! Constants for opacity etc.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL(KIND=LDP), SAVE :: FL_OLD
	REAL(KIND=LDP) BNUE
	REAL(KIND=LDP) S1
	REAL(KIND=LDP) T1,T2
!
! FG_COUNT is used to determine the average number of calls to FG_J_CMF_V10 per
! frequency.
!
	INTEGER, SAVE :: FG_COUNT
	INTEGER I
	INTEGER J_IT_COUNTER
	LOGICAL NEW_FREQ
!
! edited by Devin, 27-Oct-2026: ES-ALI state [DO_ES_ALI].  ES_ALI_ON
! is the per-call gate (evaluated just before the iteration loop);
! ES_DJ is the max relative RJ update of the current pass; ES_ALI_FB
! marks a pass that fell back to plain Lambda.  The SAVE'd counters
! accumulate per-frequency results for the once-per-iteration print.
! DJDCHI_TMP/DJDETA_TMP hold the V13 stored diagonal.
!
	LOGICAL ES_ALI_ON
	LOGICAL ES_ALI_FB
	LOGICAL ALO_OK
	REAL(KIND=LDP) ES_DJ
	REAL(KIND=LDP) RJNEW
	REAL(KIND=LDP) DEN
	REAL(KIND=LDP) DJDCHI_TMP(ND)
	REAL(KIND=LDP) DJDETA_TMP(ND)
	INTEGER, SAVE :: ES_ALI_NUNC
	INTEGER, SAVE :: ES_ALI_NFALL
	REAL(KIND=LDP), SAVE :: ES_ALI_DJMAX
!
! See COMP_J_CONT.INC to see earlier changes to this routine.
!
! 
!
! Determine outer boundary confition. For CONT_VEL (i.e. full blanketing)
! THK_CONT is always set to RDRHK_CONT, and cannot change during a model run.
!
	  IF(ATM(1)%XzV_PRES)THEN		!Hydrogen
	    T1=ATM(1)%EDGEXzV_F(1)
	  ELSE IF(ATM(4)%XzV_PRES)THEN		!Helium II
	    T1=ATM(4)%EDGEXzV_F(2)
	  ELSE
	    T1=ATM(3)%EDGEXzV_F(5)		!Helium I
	  END IF
	  IF(RDTHK_CONT .AND. FL .GT. T1)THEN
	    THK_CONT=.TRUE.
	  ELSE
	    THK_CONT=.FALSE.
	  END IF
	  C_KMS=SPEED_OF_LIGHT(I)
	  IF(FIRST_FREQ)WRITE(LUER,*)' '
!
! Option only avaialble for non-reativistic soultion, and no additinal
! points inserted.
!
	  OUT_BC_TYPE=1
	  IF(FL .LE. OUT_BC_PARAM_ONE)OUT_BC_TYPE=RD_OUT_BC_TYPE
!
       	IF(SECTION .EQ. 'CONTINUUM')THEN
	  CONT_VEL=.TRUE.
	  THK_CONT=RDTHK_CONT
	  IF(GLOBAL_LINE_SWITCH(1:5) .NE. 'BLANK' .AND. NO_VEL_FOR_CONTINUUM)THEN
	    CONT_VEL=.FALSE.
	  END IF
	ELSE
	  CONT_VEL=.FALSE.
	END IF
!
! Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
! and DDBBDT= dB/dTR .
!
	T1=HDKT*FL/T(ND)
	T2=1.0_LDP-EMHNUKT(ND)
	BNUE=TWOHCSQ*( FL**3 )*EMHNUKT(ND)/T2
	DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	DDBBDT=DBB*(T1*(1.0_LDP+EMHNUKT(ND))/T2-2.0_LDP)/T(ND)
	HFLUX_AT_IB=DBB/CHI(ND)/3.0_LDP
!
	IF(PNT_SRCE_MOD)THEN
	  T1=-HDKT*FL/TEFF_PNT_SRCE
	  IC=TWOHCSQ*( FL**3 )*EXP(T1)/(1.0_LDP-EXP(T1))
	END IF
!
! Switch to using CHI_CLUMP, ETA_CLUMP, and ESEC_CLUMP in case the model
! has clumping.
!
	CHI_CLUMP(1:ND)=CHI(1:ND)*CLUMP_FAC(1:ND)
	ETA_CLUMP(1:ND)=ETA(1:ND)*CLUMP_FAC(1:ND)
	ESEC_CLUMP(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	CHI_SCAT_CLUMP(1:ND)=CHI_SCAT(1:ND)*CLUMP_FAC(1:ND)
!
! 
!
	IF(CONT_VEL .AND. USE_FIXED_J)THEN
	  CALL RD_CONT_J(FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,
	1          ACCURATE,LUER,LU_EDD,ACCESS_F,ND,NP)
	  RSQHNU=1.0_LDP; HFLUX_AT_OB=0.001_LDP; HFLUX_AT_IB=0.0_LDP
!
	ELSE IF(.NOT. CONT_VEL .AND. THIS_FREQ_EXT)THEN
C
C Solve for the mean intensity J . We can either solve for J with or without
C Eddington factors. Generally use Eddington factors when there is many
C grid points.
C
	  CALL TUNE(IONE,'JFEAUEXT')
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,
	1                      NDEXT,CHI_CLUMP,ETA_CLUMP,CHI_SCAT_CLUMP,ND)
C
	  DO I=1,NDEXT
	    ZETAEXT(I)=ETAEXT(I)/CHIEXT(I)
	    THETAEXT(I)=ESECEXT(I)/CHIEXT(I)
	  END DO
C
	  IF(SECTION .EQ. 'CONTINUUM' .AND. FREQ_INDX .EQ. 1)FEDD=0.0_LDP
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,NDEXT
	      RJEXT(I)=0.0_LDP
	      FOLD(I)=0.0_LDP
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),T1
	    IF(ABS(T1/FL-1.0_LDP) .GT. 1.0E-10_LDP)THEN
	      WRITE(LUER,*)'Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
C
C edited by Devin, 05-Oct-2026: inside a JFNK eval sweep
C (JFNK_IN_EVAL, CONTROL_VARIABLE_MOD -- no VADAT key) a STOP would
C kill the whole run mid-JVP.  Instead flag the eval
C (JFNK_EVAL_BROKEN, consumed by JFNK_RAD_SWEEP_V1's epilogue, which
C forces a sentinel residual so JFNK_RESIDUAL_V1's |Y|>1D30 cap rejects
C it) and RETURN.  Skipping the f-iteration and the ACCESS_F advance
C is safe here: the record mismatch means the lagged-J file does not
C match this grid, and the epilogue restores ACCESS_F to
C JFNK_EDD_START_REC regardless of how far the cursor walked.
C Production path (JFNK_IN_EVAL=F) is unchanged -- still STOPs.
C
	      IF(JFNK_IN_EVAL)THEN
	        JFNK_EVAL_BROKEN=.TRUE.
	        RETURN
	      END IF
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
	  DO WHILE(INACCURATE)
	    DO I=1,NDEXT
	      SOURCEEXT(I)=ZETAEXT(I)+THETAEXT(I)*RJEXT(I)
	    END DO
	    S1=SOURCEEXT(1)
	    CALL FQCOMP_IBC(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,QEXT,FEXT,
	1            SOURCEEXT,CHIEXT,dCHIdR,AQWEXT,KQWEXT,
	1            DBB,HBC_J,HBC_S,INBC,IC,
	1            THK_CONT,DIF,NCEXT,NDEXT,NPEXT,METHOD)
	    CALL JFEAU_IBC(TA,TB,TC,DTAU,REXT,RJEXT,QEXT,FEXT,
	1          ZETAEXT,THETAEXT,CHIEXT,DBB,IC,HBC_J,HBC_S,
	1          INBC,THK_CONT,DIF,NDEXT,METHOD)
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 3 .OR. COMPUTE_EDDFAC)THEN	!Changed 8-Feb-95
	        T1=0.0_LDP
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEXT(I)),T1)
	          FOLD(I)=FEXT(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(J_IT_COUNTER .GT. 15)THEN
	         IF(FCONV_IT1_ACTIVE)THEN
	            FCONV_SUPPRESSED_CNT=FCONV_SUPPRESSED_CNT+1
	            FCONV_SUPPRESSED_MAX=MAX(FCONV_SUPPRESSED_MAX,T1)
	         ELSE
	            WRITE(LUER,*)'Possible error converging f - T1 is',T1
	            WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	         END IF
	      	 INACCURATE=.FALSE.
	      END IF
	
	    END DO
C
C Put accurate calculation of J on old grid.
C
	    CALL UNGRID(RJ,ND,RJEXT,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(K_MOM,ND,FEXT,NDEXT,POS_IN_NEW_GRID)
	    DO I=1,ND
	      K_MOM(I)=K_MOM(I)*RJ(I)
	    END DO
C
C Optput Mean intensity for subsequent iterations.
C JFNK_IN_EVAL (CONTROL_VARIABLE_MOD, no VADAT key): suppress during
C JFNK residual eval sweeps so perturbed trial J does not overwrite
C the outer-iteration records read as the ES source-function guess.
C
	    IF(.NOT. JFNK_IN_EVAL)
	1     WRITE(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),FL
C
C Update record for next frequency
	    ACCESS_F=ACCESS_F+1
C
	  CALL TUNE(ITWO,'JFEAUEXT')
C
C
C 
C
	ELSE IF(CONT_VEL .AND. ACCURATE)THEN
	  CALL TUNE(IONE,'JCONT_ACC')
C
C Interpolate the opacity and emissivity using a LINEAR interpolation
C law. CHIEXT etc. will contain the opacities etc. on the transfer grid
C with the clumping corrections. CHI_CLUMP etc refer to the appropriate
C quantities on the population grid.
C
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,NDEXT,
	1              CHI_CLUMP,ETA_CLUMP,CHI_SCAT_CLUMP,ND)
C
C NB: CHI_PREV is used to refer to the continuum opacity at the previous
C frequency. Is does not need to be multiplied by CLUMP_FAC, as it is
C compared directly to CHI_CONT. Since it is used for describing the
C variation in chi from one frequency to the next, we also do not need to
C use the extended vectors.
C
C For HBC and NBC only the first vector element is used.
C
	  CALL TUNE(IONE,'CONT_VEL')
	  NEW_FREQ=.TRUE.
	  IF(FIRST_FREQ)THEN
	    CHI_PREV(1:ND)=CHI(1:ND)
	    ETA_PREV(1:ND)=ETA(1:ND)
C
	    FEDD_PREV(1:NDEXT)=0.0_LDP		!Not required.
	    GEDD_PREV(1:NDEXT)=0.0_LDP
	    JNU_PREV(1:NDEXT)=0.0_LDP
	    N_ON_J_PREV(1:NDEXT)=0.0_LDP
	    RSQHNU_PREV(1:NDEXT)=0.0_LDP
C
	    HBC_PREV(:)=0.0_LDP		!1:3
	    NBC_PREV(:)=0.0_LDP		!1:3
	    HBC_CMF(:)=0.0_LDP		!1:3
	    NBC_CMF(:)=0.0_LDP		!1:3
	    FG_COUNT=0
	  ELSE
	    dLOG_NU=LOG(FL_OLD/FL)
	    FEDD_PREV(1:NDEXT)=FEDD(1:NDEXT)
	    GEDD_PREV(1:NDEXT)=GEDD(1:NDEXT)
	    N_ON_J_PREV(1:NDEXT)=N_ON_J(1:NDEXT)
	    JNU_PREV(1:NDEXT)=RJEXT(1:NDEXT)
	    RSQHNU_PREV(1:NDEXT)=RSQHNU(1:NDEXT)
C
	    HBC_PREV(:)=HBC_CMF(:)
	    NBC_PREV(:)=NBC_CMF(:)
	  END IF
C
	  IF(COMPUTE_EDDFAC)THEN
	    IF(FIRST_FREQ)THEN
	      RJEXT(1:NDEXT)=0.0_LDP
	      RJEXT_ES(1:NDEXT)=0.0_LDP
	      FOLD(1:NDEXT)=0.0_LDP
	    END IF
	  ELSE
C
C Same JFNK_IN_EVAL cold start as the CONT_VEL branch below: FOLD is
C SAVE'd RADIATION_MOD state carried between sweeps, so without the
C seed an eval's iter-1 f convergence test compares fresh FEDD against
C whatever the previous eval/production sweep left -- eval-history
C dependence.  Production path (JFNK_IN_EVAL=F) unchanged.
C
	    IF(JFNK_IN_EVAL .AND. FIRST_FREQ)FOLD(1:NDEXT)=0.0_LDP
	    READ(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),T1
	    IF(ABS(T1/FL-1.0_LDP) .GT. 1.0E-10_LDP)THEN
	      WRITE(LUER,*)'Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
C
C JFNK_IN_EVAL soften -- see the first EDDFACTOR read site above:
C flag the eval broken and RETURN so JFNK_RAD_SWEEP_V1's epilogue can
C reject it; the cursor is restored there, so skipping the record walk
C is safe.  Production path unchanged (still STOPs).
C
	      IF(JFNK_IN_EVAL)THEN
	        JFNK_EVAL_BROKEN=.TRUE.
	        RETURN
	      END IF
	      STOP
	    END IF
	  END IF
C
C If we are using incoherent electron scattering, RJEXT_ES must be available.
C
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJEXT_ES(I),I=1,NDEXT),T1
	    IF(ABS(T1/FL-1.0_LDP) .GT. 1.0E-10_LDP)THEN
	      WRITE(LUER,*)'Error - incorrect reading of ES_J_CONV in COMO_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete ES_J_CONV'
C
C JFNK_IN_EVAL soften -- same pattern as the EDDFACTOR sites above:
C RJEXT_ES is read (not iterated) so a wrong record would silently
C corrupt the eval; flag it and RETURN.  Production path unchanged.
C
	      IF(JFNK_IN_EVAL)THEN
	        JFNK_EVAL_BROKEN=.TRUE.
	        RETURN
	      END IF
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
	  DO WHILE(INACCURATE)
C
	     IF(COHERENT_ES)THEN
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1           ESECEXT(1:NDEXT)*RJEXT(1:NDEXT)
	     ELSE
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1           ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	     END IF
C
C NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
C
	     CALL TUNE(IONE,'FG_J_CMF_ACC')
	     CALL FG_J_CMF_V13(TA,CHIEXT,ESECEXT,
	1            VEXT,SIGMAEXT,REXT,PEXT,TC,FEDD,
	1            AQWEXT,HQWEXT,KQWEXT,NQWEXT,HMIDQWEXT,NMIDQWEXT,
	1            INBC,HBC_CMF(1),IPLUS,FL,dLOG_NU,
	1            INNER_BND_METH,DBB,IC,
	1            VDOP_VEC_EXT,DELV_FRAC_FG,REXT_FAC,
	1            METHOD,FG_SOL_OPTIONS,THK_CONT,
	1            FIRST_FREQ,NEW_FREQ,NCEXT,NPEXT,NDEXT)
	     CALL TUNE(ITWO,'FG_J_CMF_ACC')
	     FG_COUNT=FG_COUNT+1
C
	     IF(COHERENT_ES)THEN
	       TA(1:NDEXT)=ETAEXT(1:NDEXT)
	     ELSE
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1            ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF_ACC')
	     CALL MOM_J_CMF_V11(TA,CHIEXT,ESECEXT,VEXT,SIGMAEXT,REXT,
	1  	       RJEXT,RSQHNU,VDOP_VEC_EXT,DELV_FRAC_MOM,
	1              FL,dLOG_NU,INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
	1              N_TYPE,H_CHK_OPTION,METHOD,COHERENT_ES,OUT_BC_TYPE,
	1              FIRST_FREQ,NEW_FREQ,NCEXT,NPEXT,NDEXT)
	     CALL TUNE(ITWO,'MOM_J_CMF_ACC')
             IF(.NOT. DIF)HFLUX_AT_IB=0.5_LDP*IC*(0.5_LDP+INBC)-INBC*RJEXT(NDEXT)
             HFLUX_AT_OB=HBC_CMF(1)*RJEXT(1)
C
C We set NEW_FREQ to false so that FG_J_CMF continues to use the same
C AV_PREV and CV_PREV. NEW_FREQ must be set to true again outside the
C F iteration loop.
C
	     NEW_FREQ=.FALSE.
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0_LDP
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(J_IT_COUNTER .GT. 10)THEN
	         IF(FCONV_IT1_ACTIVE)THEN
	            FCONV_SUPPRESSED_CNT=FCONV_SUPPRESSED_CNT+1
	            FCONV_SUPPRESSED_MAX=MAX(FCONV_SUPPRESSED_MAX,T1)
	         ELSE
	            WRITE(LUER,*)'Possible error converging f - T1 is',T1
	            WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	         END IF
	      	 INACCURATE=.FALSE.
	      END IF
	    END DO
C
C Output RJ for subsequent iterations.
C JFNK_IN_EVAL (CONTROL_VARIABLE_MOD): suppress during JFNK residual
C eval sweeps -- see note at the first LU_EDD write site.
C
	    IF(.NOT. JFNK_IN_EVAL)
	1     WRITE(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),FL
C
C Store J on the normal mesh. No interpolation is involved here as
C it is assumed that the fine grid was created by the addition of extra
C points only.
C
	    CALL UNGRID(RJ,ND,RJEXT,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(RJ_ES,ND,RJEXT_ES,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(K_MOM,ND,FEDD,NDEXT,POS_IN_NEW_GRID)
C
C Compute K for use in computing mechanical energy loss.
C
	    K_MOM(1:ND)=RJ(1:ND)*K_MOM(1:ND)
C
	    IF(COHERENT_ES)THEN
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ(1:ND)
	    ELSE
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ_ES(1:ND)
	    END IF
C
C Update record for next frequency
C
	    ACCESS_F=ACCESS_F+1
	    FL_OLD=FL
!
! Note that TC is one the EXTENDED grid, hence we access its value at the
! inner boundary using NDEXT.
!
	  IF(LST_ITERATION)THEN
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=200.0_LDP*(RJ(1)-TC(1))/T1
	    T2=ABS(RJ(ND))+ABS(TC(NDEXT))
	    IF(T2 .NE. 0)T2=200.0_LDP*(RJ(ND)-TC(NDEXT))/T2
	    IF(FIRST_FREQ)THEN
	      OPEN(UNIT=LU_JCOMP,STATUS='UNKNOWN',FILE='J_COMP')
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,'(A)')'Comparison of J at Outer and Inner',
	1       ' boundaries computed using Moments and Ray techniques.'
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,
	1       '(3X,A,7X,A,7X,A,6X,A,5X,A,6X,A,5X,A,6X,A,5X,A)')
	1       'Indx','Nu','J(mom)','J(ray)','%Diff','HBC_CMF',
	1       'J(mom)','J(ray)','%Diff'
	    END IF
	    WRITE(LU_JCOMP,'(I7,ES16.6,2ES12.4,F10.2,3ES12.4,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(NDEXT),T2
	    IF(FREQ_INDX .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(FREQ_INDX .EQ. NCF)THEN
	    WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FLOAT(FG_COUNT)/FLOAT(NCF)
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=ABS(200.0_LDP*(RJ(1)-TC(1))/T1)
	    IF(T1 .GT. 100.0_LDP)THEN
	      WRITE(LUER,'(A)')'***************************************************************************'
	      WRITE(LUER,'(A)')' Error --- J(mom) and J(ray) differ by more than 100% for last frequency'
	      WRITE(LUER,'(A)')' It is STRONGLY suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error and/or use plt_jh'
	      WRITE(LUER,'(A)')'***************************************************************************'
	    ELSE IF(T1 .GT. 50.0_LDP)THEN
	      WRITE(LUER,'(A)')' Error --- J(mom) and J(ray) differ by more than 50% for last frequency'
	      WRITE(LUER,'(A)')' It is strongly suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    ELSE IF(T1 .GT. 20.0_LDP)THEN
	      WRITE(LUER,'(A)')' Error --- J(mom) and J(ray) differ by more than 20% for last frequency'
	      WRITE(LUER,'(A)')' Although this is nlikely to effect the colution, it is suggested that'
	      WRITE(LUER,'(A)')' you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    END IF
	  END IF
	  CALL TUNE(ITWO,'CONT_VEL')
C
C Set up for the compuation of the observes flux. LST_ITERATION is
C TRUE if FLUX_CAL_ONLY is true (single iteration with coherent,
C last iteration if non-coherent).
C
	    IF( (LST_ITERATION .AND. .NOT. LAMBDA_ITERATION .AND.
	1         MAXCH .LT. 100.0_LDP .AND. .NOT. DO_NONMONO_KNOP .AND.
	1         .NOT. SN_MODEL) )THEN
C
C Quick and dirty method to ge an extended DENSITY vector. Will use TB in
C the call to CMF_FORM_SOL.
C
	       CALL EXTEND3OPAC(TA,TB,TC,COEF,INDX,NDEXT,
	1              DENSITY,DENSITY,DENSITY,ND)
C
	      IF(COHERENT_ES)THEN
	        TA(1:NDEXT)=ETAEXT(1:NDEXT)+ESECEXT(1:NDEXT)*RJEXT(1:NDEXT)
	      ELSE
	        TA(1:NDEXT)=ETAEXT(1:NDEXT)+ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	      END IF
C
C NB Using TA for ETA, U for P_OBS (temporay measure), I for NP_OBS.
C
	      CALL TUNE(IONE,'CMF_FORM_SOL')
	      CALL CMF_FORM_SOL_V2(TA,CHIEXT,ESECEXT,
	1                 TB,VEXT,SIGMAEXT,REXT,PEXT,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,CMF_FORM_OPTIONS,
	1                 FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING,
	1                 FIRST_FREQ,NCEXT,NPEXT,NDEXT)
	      CALL TUNE(ITWO,'CMF_FORM_SOL')
	    ELSE IF(FIRST_FREQ)THEN
C
C So as defined for normal OBSFLUX calculation.
C
	      NP_OBS=NPEXT
	      P_OBS(1:NPEXT)=PEXT(1:NPEXT)
	      RMAX_OBS=R(1)
	      V_AT_RMAX=V(1)
	    END IF
!
	  CALL TUNE(ITWO,'JCONT_ACC')
C
C
C 
C
	ELSE IF(CONT_VEL)THEN
C
C NB: CHI_PREV is used to refer to the continuum opacity at the previous
C frequency. Is does not need to be multiplied by CLUMP_FAC, as it is compared
C directly to CHI_CONT.
C
C For HBC and NBC only the first vector element is used.
C
	  CALL TUNE(IONE,'CONT_VEL')
	  NEW_FREQ=.TRUE.
	  IF(FIRST_FREQ)THEN
	    DO I=1,ND
	      CHI_PREV(I)=CHI(I)
	      ETA_PREV(I)=ETA(I)
	      FEDD_PREV(I)=0.0_LDP		!Not required.
	      GEDD_PREV(I)=0.0_LDP
	      JNU_PREV(I)=0.0_LDP
	      N_ON_J_PREV(I)=0.0_LDP
	      RSQHNU_PREV(I)=0.0_LDP
	    END DO
!
! At the bluest frequency FL_OLD is undefined, so dLOG_NU cannot be
! evaluated as LOG(FL_OLD/FL).  The correct one-sided spacing is
! LOG(FL/NU_NEXT_FREQ) (NU_NEXT_FREQ set by CMFGEN_SUB).  Without this
! the non-monotonic chi_hat = chi + |GAM|/dLOG_NU sees a garbage/zero
! dLOG_NU at frequency 1 and produces a non-finite moment solution.
! Gated on DO_NONMONO_KNOP, and on DO_NONMONO_SOB with USE_LAM_ES=F
! (MOM_J_CMF_NONMONO_V1 is the field solver there too, 22-Oct-2026),
! so the stock path is unchanged.
!
	    IF((DO_NONMONO_KNOP .OR.
	1       (DO_NONMONO_SOB .AND. .NOT. USE_LAM_ES) .OR.
	1       (DO_NONMONO_SOB .AND. DO_SOB_FG_KNOP)) .AND.
	2       NU_NEXT_FREQ .GT. 0.0_LDP .AND.
	3       NU_NEXT_FREQ .LT. FL)THEN
	      dLOG_NU=LOG(FL/NU_NEXT_FREQ)
	    END IF
	    HBC_PREV(:)=0.0_LDP		!1:3
	    NBC_PREV(:)=0.0_LDP		!1:3
	    HBC_CMF(:)=0.0_LDP		!1:3
	    NBC_CMF(:)=0.0_LDP		!1:3
	    FG_COUNT=0
!
! ES-ALI per-iteration counters [DO_ES_ALI]; reset with FG_COUNT so
! the NCF print covers exactly one great iteration.  Inert when the
! feature is off.
!
	    ES_ALI_NUNC=0
	    ES_ALI_NFALL=0
	    ES_ALI_DJMAX=0.0_LDP
	  ELSE
	    dLOG_NU=LOG(FL_OLD/FL)
	    DO I=1,ND
	      FEDD_PREV(I)=FEDD(I)
	      GEDD_PREV(I)=GEDD(I)
	      N_ON_J_PREV(I)=N_ON_J(I)
	      JNU_PREV(I)=RJ(I)
	      RSQHNU_PREV(I)=RSQHNU(I)
	    END DO
	    HBC_PREV(:)=HBC_CMF(:)
	    NBC_PREV(:)=NBC_CMF(:)
	  END IF
C
	  IF(COMPUTE_EDDFAC)THEN
	    IF(FIRST_FREQ)THEN
	      DO I=1,ND
	        RJ(I)=0.0_LDP
	        RJ_ES(I)=0.0_LDP
	        FOLD(I)=0.0_LDP
	      END DO
	    END IF
	  ELSE
C
C edited by Devin, 05-Oct-2026: under COMPUTE_EDDFAC the FIRST_FREQ
C block above cold-starts the f-iteration with FOLD=0; with stored
C Eddington factors (COMPUTE_EDDFAC=F) FOLD is SAVE'd module state in
C RADIATION_MOD and is never reset, so an eval sweep's iter-1
C convergence test would compare the fresh FEDD against whatever the
C previous eval or production sweep left -- eval-history dependence.
C Seed the cold start the same way the COMPUTE_EDDFAC branch does.
C Gated on JFNK_IN_EVAL so production warm starts are bit-identical.
C
	    IF(JFNK_IN_EVAL .AND. FIRST_FREQ)FOLD(1:ND)=0.0_LDP
	    READ(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),T1
	    IF(ABS(T1/FL-1.0_LDP) .GT. 1.0E-10_LDP)THEN
	      WRITE(LUER,*)'Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
C
C JFNK_IN_EVAL soften -- see the first EDDFACTOR read site above:
C flag the eval broken and RETURN; the sweep epilogue restores the
C ACCESS_F cursor and forces the sentinel residual.  Production path
C unchanged (still STOPs).
C
	      IF(JFNK_IN_EVAL)THEN
	        JFNK_EVAL_BROKEN=.TRUE.
	        RETURN
	      END IF
	      STOP
	    END IF
	  END IF
C
C If we are using incoherent electron scattering, RJEXT_ES must be available.
C
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJ_ES(I),I=1,ND),T1
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of ES_J_CONV in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete ES_J_CONV'
C
C JFNK_IN_EVAL soften -- same pattern as the EDDFACTOR sites above:
C RJ_ES is a read (not iterated) input to the source function, so a
C wrong record would silently corrupt the eval; flag it and RETURN.
C Production path unchanged (still STOPs).
C
	      IF(JFNK_IN_EVAL)THEN
	        JFNK_EVAL_BROKEN=.TRUE.
	        RETURN
	      END IF
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
C
C edited by Devin, 27-Oct-2026: ES-ALI gate [DO_ES_ALI].  The
C Olson-Auer-Buchler update below uses the diagonal stored by
C FG_J_CMF_V13, so it is restricted to the V13 path: DO_NONMONO_SOB +
C USE_LAM_ES + COHERENT_ES and not DO_SOB_FG_KNOP (which runs
C FG_J_CMF_NONMONO_V1 instead).  RJ entering the loop is the
C EDDFACTOR-read previous-iteration J (warm seed), and TA is rebuilt
C from the updated RJ at the top of every pass.  DO_ES_ALI=F leaves
C the loop bit-identical to stock.
C
	  ES_ALI_ON=DO_ES_ALI .AND. DO_NONMONO_SOB .AND. USE_LAM_ES .AND.
	1     COHERENT_ES .AND. .NOT. DO_SOB_FG_KNOP
	  ES_ALI_FB=.FALSE.
	  ES_DJ=0.0_LDP
	  DO WHILE(INACCURATE)
C
	     IF(ES_ALI_ON)ES_DJ=0.0_LDP
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ_ES(1:ND)
	     END IF
C
C NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
C
	     CALL TUNE(IONE,'FG_J_CMF')
	     IF(PLANE_PARALLEL_NO_V)THEN
	        IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling FCOMP_PP'
	        SOURCE(1:ND)=TA(1:ND)/CHI_CLUMP(1:ND)
	        CALL FCOMP_PP_V2(R,TC,FEDD,SOURCE,CHI_CLUMP,IPLUS,HBC_CMF,
	1               NBC_CMF,INBC,DBB,IC,THK_CONT,DIF,ND,NC,METHOD)
	     ELSE IF(PLANE_PARALLEL)THEN
	        IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling PP_FORM_CMF_V2'
	        CALL PP_FORM_CMF_V2(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1               TC,TB,FEDD,GEDD,N_ON_J,INBC,
	1               HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1               IPLUS,FL,dLOG_NU,DIF,DBB,IC,VDOP_VEC,DELV_FRAC_FG,
	1               METHOD,FG_SOL_OPTIONS,THK_CONT,INCL_INCID_RAD,
	1               FIRST_FREQ,NEW_FREQ,N_TYPE,NC,ND)
!
	     ELSE IF(USE_FORMAL_REL)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)
	1            WRITE(LUER,*)'Calling CMF_FORMAL_REL_V4 in COMP_J_BLANK'
	       CALL CMF_FORMAL_REL_V4
	1                 (TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P,
	1                  TC,FEDD,HFLUX_AT_IB,HFLUX_AT_OB,IPLUS,
	1                  FL,dLOG_NU,BNUE,DBB,
	1                  INNER_BND_METH,THK_CONT,
	1                  VDOP_VEC,DELV_FRAC_FG,REXT_FAC,
	1                  METHOD,FIRST_FREQ,NEW_FREQ,NC,NP,ND)
!
	     ELSE
	       IF(DO_NONMONO_KNOP .OR.
	1         (DO_NONMONO_SOB .AND. USE_LAM_ES .AND. DO_SOB_FG_KNOP))THEN
!	         Knop et al. (2009, Paper I) CMF formal solver for
!	         non-monotonic velocity fields: chi_hat = chi + |GAM|/dLOG_NU stays
!	         positive in compression zones (sigma < -1), and the
!	         signed frequency coupling is carried in the source term
!	         with the upwind direction chosen by the sign of GAM.
!	         Uses the unclamped SIGMA_TRUE; SIGMA stays clamped for
!	         the rest of the infrastructure.
!	         DO_SOB_FG_KNOP (opt-in) brings the same solver to the
!	         SOB+LAM_ES=T path -- V13 is expansion-biased at sigma<-1.
	         IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)THEN
	           WRITE(LUER,'(/,A)')' -------------------------------- KNOP SOLVER --------------------------------'
	           WRITE(LUER,*)'Calling FG_J_CMF_NONMONO_V1 (Knop et al. 2009)'
	         END IF
!
! Pass nu(k+1) to the solver so the compressing branch uses the
! forward spacing ln(nu_k/nu_{k+1}) (asymmetric-spacing fix,
! Oct-2026).  NU_NEXT_FREQ is set by CMFGEN_SUB each frequency.
!
	         CALL SET_FG_NU_NEXT(NU_NEXT_FREQ)
	         CALL FG_J_CMF_NONMONO_V1(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P,
	1                  TC,FEDD,AQW,HQW,KQW,NQW,HMIDQW,NMIDQW,
	1                  INBC,HBC_CMF(1),IPLUS,FL,dLOG_NU,
	1                  INNER_BND_METH,DBB,IC,
	1                  VDOP_VEC,DELV_FRAC_FG,REXT_FAC,
	1                  METHOD,FG_SOL_OPTIONS,THK_CONT,
	1                  FIRST_FREQ,NEW_FREQ,NC,NP,ND,SIGMA_TRUE,NCF)
	       ELSE
	         IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)THEN
	           IF(DO_NONMONO_SOB)WRITE(LUER,'(/,A)')
	1            ' --------------------------------- SOB SOLVER ---------------------------------'
	           WRITE(LUER,*)'Calling FG_J_CMF_V13 in COMP_J_BLANK'
	         END IF
	         CALL FG_J_CMF_V13(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P,
	1                  TC,FEDD,AQW,HQW,KQW,NQW,HMIDQW,NMIDQW,
	1                  INBC,HBC_CMF(1),IPLUS,FL,dLOG_NU,
	1                  INNER_BND_METH,DBB,IC,
	1                  VDOP_VEC,DELV_FRAC_FG,REXT_FAC,
	1                  METHOD,FG_SOL_OPTIONS,THK_CONT,
	1                  FIRST_FREQ,NEW_FREQ,NC,NP,ND)
	       END IF
!
	     END IF
	     CALL TUNE(ITWO,'FG_J_CMF')
	     FG_COUNT=FG_COUNT+1
C
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ_ES(1:ND)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF')
	     IF(PLANE_PARALLEL_NO_V)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling MOM_J_PP_V1'
	       CALL MOM_J_PP_V1(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1                  R,FEDD,RJ,RSQHNU,HBC_CMF,NBC_CMF,INBC,
	1                  FL,DIF,DBB,IC,METHOD,COHERENT_ES,
	1                  IZERO,FIRST_FREQ,NEW_FREQ,ND)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
	       IF(.NOT. DIF)HFLUX_AT_IB=0.5_LDP*IC*(0.5_LDP+INBC)-INBC*RJ(ND)
	       IF(LST_ITERATION)THEN
	         DO I=1,ND
	           TA(I)=RJ(I)*R(I)*R(I)
	         END DO
	         T1=HFLUX_AT_IB*R(ND)*R(ND)
	         T2=HFLUX_AT_OB/RJ(1)
	         CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,R,V,ND,FIRST_FREQ,'NORMAL')
	       END IF
	     ELSE IF(PLANE_PARALLEL)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling PP_MOM_CMF_V1'
	       CALL PP_MOM_CMF_V1(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1                  FEDD,GEDD,N_ON_J,RJ,RSQHNU,
	1                  VDOP_VEC,DELV_FRAC_MOM,
	1                  INBC,HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1                  FL,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,COHERENT_ES,
	1                  FIRST_FREQ,NEW_FREQ,ND)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
	       IF(.NOT. DIF)HFLUX_AT_OB=0.5_LDP*IC*(0.5_LDP+INBC)-INBC*RJ(ND)
	     ELSE IF(USE_DJDT_RTE)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling MOM_J_DDT_V6'
	       CALL MOM_J_DDT_V6(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1              V,R,FEDD,RJ,RSQHNU,DJDt_TERM,
	1              HFLUX_AT_IB,HFLUX_AT_OB,
	1              VDOP_VEC,DELV_FRAC_MOM,FL,dLOG_NU,DBB,
	1              XM_CHK_OPTION,J_CHK_OPTION,H_CHK_OPTION,
	1              INNER_BND_METH,OUTER_BND_METH,
	1              METHOD,COHERENT_ES,FIRST_FREQ,NEW_FREQ,
	1              INCL_DJDT_TERMS,USE_DR4JDT,DJDT_RELAX_PARAM,NC,NP,ND,NCF)
	       HBC_CMF(1)=HFLUX_AT_OB/RJ(1)
	     ELSE IF(USE_J_REL)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling MOM_JREL_V9'
	       CALL MOM_JREL_V9(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1             RJ,RSQHNU,HFLUX_AT_IB,HFLUX_AT_OB,
	1             VDOP_VEC,DELV_FRAC_MOM,
	1             FL,dLOG_NU,DBB,
	1             XM_CHK_OPTION,J_CHK_OPTION,H_CHK_OPTION,IB_STAB_FACTOR,
	1             INNER_BND_METH,OUTER_BND_METH,
	1             METHOD,COHERENT_ES,N_TYPE,
	1             INCL_ADVEC_TERMS_IN_TRANS_EQ,INCL_REL_TERMS,FIRST_FREQ,ND)
	       IF(LST_ITERATION)THEN
	         DO I=1,ND
	           TA(I)=RJ(I)*R(I)*R(I)
	         END DO
	         T1=HFLUX_AT_IB*R(ND)*R(ND)
	         T2=HFLUX_AT_OB/RJ(1)
	         CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,R,V,ND,FIRST_FREQ,'NORMAL')
	       END IF
	     ELSE IF(USE_LAM_ES)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)
	1         'Field solver: ray solutions with lambda'//
	2         ' iteration (USE_LAM_ES); the moment equations'//
	3         ' are not solved'
	       IF(ES_ALI_ON)THEN
!
! edited by Devin, 27-Oct-2026: Olson-Auer-Buchler Jacobi update of
! the coherent-ES iterate [DO_ES_ALI].  TA was built from the
! previous RJ, so TC is Lambda[RJ]; the V13 stored diagonal
! DJDETA = dJ/dETA converts the residual into an approximate solve
! of J = Lambda[ETA + chi_scat J]:
! RJNEW = RJ + (TC-RJ)/DEN with DEN = 1 - DJDETA*chi_scat, floored
! at 0.01; a non-positive RJNEW falls back to the Lambda value TC.
! ES_DJ records the max relative update for the exit test below.  If
! the diagonal is unavailable the pass degrades to plain RJ=TC, the
! frequency is counted as a fallback, and no extra passes are forced.
!
	         CALL GET_FG13_ALO_DIAG(DJDCHI_TMP,DJDETA_TMP,ND,ALO_OK)
	         IF(ALO_OK)THEN
	           DO I=1,ND
	             DEN=MAX(1.0_LDP-DJDETA_TMP(I)*CHI_SCAT_CLUMP(I),
	1                0.01_LDP)
	             RJNEW=RJ(I)+(TC(I)-RJ(I))/DEN
	             IF(RJNEW .LE. 0.0_LDP)RJNEW=TC(I)
	             ES_DJ=MAX(ES_DJ,
	1                ABS(RJNEW-RJ(I))/MAX(ABS(RJNEW),1.0E-300_LDP))
	             RJ(I)=RJNEW
	           END DO
	         ELSE
	           RJ(1:ND)=TC(1:ND)
	           ES_ALI_FB=.TRUE.
	         END IF
	       ELSE
	         RJ(1:ND)=TC(1:ND)
	       END IF
!
! Fill Eddington factors (K_ON_J, NMID_ON_J, NMID_ON_HMID) from the
! ray moments stored by FG_J_CMF_NONMONO_V1.  Normally these are
! filled by MOM_J_CMF_V11, which is NOT called when USE_LAM_ES=T.
! Gated on DO_NONMONO_KNOP: USE_LAM_ES is a standard CMFGEN option
! and the stock path is unchanged without the flag.  DO_SOB_FG_KNOP
! runs the same ray solver, so its stored moments feed here too.
!
	       IF(DO_NONMONO_KNOP .OR.
	1         (DO_NONMONO_SOB .AND. DO_SOB_FG_KNOP))THEN
	         CALL GET_MOMS_NON_REL(R,V,FL,N_TYPE,ND)
	       END IF
	       IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	          HFLUX_AT_IB=0.0_LDP
	          DO I=1,NC_PNT_SRCE
	            HFLUX_AT_IB=HFLUX_AT_IB+0.5_LDP*HQW(ND,I)*IC
	          END DO
! dev7: upstream typo fixed - stock tested 'ZERO FLUX' (space), which
! never matches the VADAT value 'ZERO_FLUX', so this branch was dead
! code and the diffusion-flux ELSE always ran (spurious HFLUX_AT_IB).
	       ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	          HFLUX_AT_IB=0.0_LDP
! With ZERO_FLUX now live, the generic branch that sets HFLUX_AT_OB is
! skipped; under USE_LAM_ES the Eddington factor HBC_CMF(1) is stale
! because MOM_J_CMF is never called.  The Knop ray solve measures the
! emergent flux directly, so use HNU_AT_OB for the SOB/RLUMST outer
! boundary.  Stock path (DO_NONMONO_KNOP=F) is unchanged.
	          IF(DO_NONMONO_KNOP .OR.
	1            (DO_NONMONO_SOB .AND. DO_SOB_FG_KNOP))
	1            HFLUX_AT_OB=HNU_AT_OB
	       ELSE IF(.NOT. USE_FORMAL_REL)THEN
	         IF(.NOT. DIF)HFLUX_AT_IB=0.5_LDP*IC*(0.5_LDP+INBC)-INBC*RJ(ND)
                 HFLUX_AT_OB=HBC_CMF(1)*RJ(1)
	       END IF
	       CALL GET_RSQH_REL(RSQHNU,R,V,FL,ND)
	       IF(LST_ITERATION)THEN
	         DO I=1,ND
	           TA(I)=RJ(I)*R(I)*R(I)
	         END DO
	         T1=HFLUX_AT_IB*R(ND)*R(ND)
	         T2=HFLUX_AT_OB/RJ(1)
	         CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,R,V,ND,FIRST_FREQ,'NORMAL')
	       END IF
	     ELSE
	       IF(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB)THEN
!	        Non-monotonic moment solver (Knop et al. 2009 Paper I
!	        generalized opacity split).  Uses SIGMA_TRUE (unclamped) and sign-
!	        dependent upwind (J_PREV for GAM>=0, J_NEXT for GAM<0).
!	        Called under DO_NONMONO_KNOP, and under DO_NONMONO_SOB when
!	        the opt-in USE_LAM_ES=F moment path is active (12-Oct-2026):
!	        MOM_J_CMF_V11 assumes expanding flow and diverges on this
!	        velocity field.  V13 still runs each frequency above, so
!	        GET_MOMS_NON_REL inside the solver reads fresh N_ON_J ray
!	        moments.  Default (no non-mono flag) keeps MOM_J_CMF_V11.
	         IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)
	1             WRITE(LUER,*)'Field solver: moment equations'//
	2             ' (MOM_J_CMF_NONMONO_V1)'
!
! Pass nu(k+1) to the moment solver so compressing depths use the
! forward spacing ln(nu_k/nu_{k+1}) (asymmetric-spacing fix,
! Oct-2026; same convention as SET_FG_NU_NEXT above).
!
	         CALL SET_MOM_NU_NEXT(NU_NEXT_FREQ)
	         CALL MOM_J_CMF_NONMONO_V1(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1              V,SIGMA,SIGMA_TRUE,R,
	1              RJ,RSQHNU,VDOP_VEC,DELV_FRAC_MOM,
	1              FL,dLOG_NU,INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
	1              N_TYPE,H_CHK_OPTION,METHOD,COHERENT_ES,OUT_BC_TYPE,
	1              FIRST_FREQ,NEW_FREQ,NC,NP,ND,FREQ_INDX,NCF)
	       ELSE
	         IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)
	1             WRITE(LUER,*)'Calling MOM_J_CMF_V11'
	         CALL MOM_J_CMF_V11(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1              RJ,RSQHNU,VDOP_VEC,DELV_FRAC_MOM,
	1              FL,dLOG_NU,INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
	1              N_TYPE,H_CHK_OPTION,METHOD,COHERENT_ES,OUT_BC_TYPE,
	1              FIRST_FREQ,NEW_FREQ,NC,NP,ND)
	       END IF
	       IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	          HFLUX_AT_IB=0.0_LDP
	          DO I=1,NC_PNT_SRCE
	            HFLUX_AT_IB=HFLUX_AT_IB+0.5_LDP*HQW(ND,I)*IC
	          END DO
! dev7: same upstream 'ZERO FLUX' -> 'ZERO_FLUX' typo fix as above;
! the dead branch had let the generic ELSE overwrite HFLUX_AT_IB.
	       ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	          HFLUX_AT_IB=0.0_LDP
	       ELSE IF(INNER_BND_METH(1:3) .EQ. 'DIF')THEN
	       ELSE
	          HFLUX_AT_IB=0.5_LDP*IC*(0.5_LDP+INBC)-INBC*RJ(ND)
	       END IF
!	       IF(.NOT. DIF)HFLUX_AT_IB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
               HFLUX_AT_OB=HBC_CMF(1)*RJ(1)
	       IF(LST_ITERATION .AND. WRITE_JH)THEN
	         DO I=1,ND
	           TA(I)=RJ(I)*R(I)*R(I)
	         END DO
	         T1=HFLUX_AT_IB*R(ND)*R(ND)
	         T2=HFLUX_AT_OB/RJ(1)
	         CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,R,V,ND,FIRST_FREQ,'NORMAL')
	       END IF
	     END IF
	     CALL TUNE(ITWO,'MOM_J_CMF')
C
C We set NEW_FREQ to false so that FG_J_CMF continues to use the same
C AV_PREV and CV_PREV. NEW_FREQ must be set to true again outside the
C F iteration loop.
C
	     NEW_FREQ=.FALSE.
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0_LDP
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
C
C edited by Devin, 27-Oct-2026: under DO_ES_ALI the ES iterate keeps
C passing until the last update is below ES_ALI_TOL or the pass cap
C ES_ALI_MAX is reached (clamped to 20 at read time, so the existing
C forced exit still bounds the loop).  Fallback frequencies are not
C forced to extra passes.
C
	        IF(ES_ALI_ON .AND. .NOT. ES_ALI_FB .AND.
	1          ES_DJ .GT. ES_ALI_TOL .AND.
	2          J_IT_COUNTER .LT. ES_ALI_MAX)INACCURATE=.TRUE.
	      END IF
C
	      IF(J_IT_COUNTER .GT. 20)THEN
	         IF(FCONV_IT1_ACTIVE)THEN
	            FCONV_SUPPRESSED_CNT=FCONV_SUPPRESSED_CNT+1
	            FCONV_SUPPRESSED_MAX=MAX(FCONV_SUPPRESSED_MAX,T1)
	         ELSE
	            WRITE(LUER,*)'Possible error converging f - T1 is',T1
	            WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	         END IF
	      	 INACCURATE=.FALSE.
	      END IF
	    END DO
C
C edited by Devin, 27-Oct-2026: accumulate the per-frequency ES-ALI
C result for the once-per-iteration print [DO_ES_ALI].
C
	    IF(ES_ALI_ON)THEN
	      ES_ALI_DJMAX=MAX(ES_ALI_DJMAX,ES_DJ)
	      IF(ES_ALI_FB)THEN
	        ES_ALI_NFALL=ES_ALI_NFALL+1
	      ELSE IF(ES_DJ .GT. ES_ALI_TOL)THEN
	        ES_ALI_NUNC=ES_ALI_NUNC+1
	      END IF
	    END IF
C
	    IF(COHERENT_ES)THEN
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ(1:ND)
	    ELSE
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ_ES(1:ND)
	    END IF
C
C Output RJ for subsequent iterations.
C JFNK_IN_EVAL (CONTROL_VARIABLE_MOD): suppress during JFNK residual
C eval sweeps -- see note at the first LU_EDD write site.
C
	    IF(.NOT. JFNK_IN_EVAL)
	1     WRITE(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),FL
C
C Compute K for use in computing mechanical energy loss.
C
	    DO I=1,ND
	      K_MOM(I)=RJ(I)*FEDD(I)
	    END DO
C
C Update record for next frequency
C
	    ACCESS_F=ACCESS_F+1
	    FL_OLD=FL
C
C Set up for the compuation of the observes flux. LST_ITERATION is
C TRUE if FLUX_CAL_ONLY is true (single iteration with coherent,
C last iteration if non-coherent).
C
	    IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V)THEN
!
! So as defined for normal OBSFLUX calculation. HQW_AT_RMAX is initially set
! to JQW. Thus we need to multiply by MU to get the actual H weights at the
! outer boundary. For a plane-parallel atmosphere, RMAX_OBS is only scaling
! constant. Setting its value to ND means that the observed luminosity should
! correspond to the luminosity in VADAT (in absence of significant velocity
! effects).
!
	      IF(FIRST_FREQ)THEN
	        NP_OBS=NC
	        CALL GAULEG(RZERO,RONE,MU_AT_RMAX,HQW_AT_RMAX,NC)
	        HQW_AT_RMAX(1:NC)=HQW_AT_RMAX(1:NC)*MU_AT_RMAX(1:NC)
	        RMAX_OBS=R(ND)
	        V_AT_RMAX=V(1)
	        IF(PLANE_PARALLEL_NO_V)V_AT_RMAX=0.0_LDP
	      END IF
	    ELSE IF( (LST_ITERATION .AND. .NOT. LAMBDA_ITERATION .AND.
	1         MAXCH .LT.  100.0_LDP .AND. .NOT. DO_NONMONO_KNOP .AND.
	1         .NOT. SN_MODEL .AND. ND .LT. 200) )THEN
	      IF(COHERENT_ES)THEN
     	        TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ(1:ND)
	      ELSE
     	        TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ_ES(1:ND)
	      END IF
C
C NB Using TA for ETA, U for P_OBS (temporay measure), I for NP_OBS.
C
	      CALL TUNE(IONE,'CMF_FORM_SOL')
	      CALL CMF_FORM_SOL_V2(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1                 DENSITY,V,SIGMA,R,P,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,CMF_FORM_OPTIONS,
	1                 FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING,
	1                 FIRST_FREQ,NC,NP,ND)
	      CALL TUNE(ITWO,'CMF_FORM_SOL')
	    ELSE IF(FIRST_FREQ)THEN
!
! So as defined for normal OBSFLUX calculation.
!
	      NP_OBS=NP
	      P_OBS(1:NP)=P(1:NP)
	      RMAX_OBS=R(1)
	      V_AT_RMAX=V(1)
	    END IF
C
	  IF(LST_ITERATION)THEN
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=200.0_LDP*(RJ(1)-TC(1))/T1
	    T2=ABS(RJ(ND))+ABS(TC(ND))
	    IF(T2 .NE. 0)T2=200.0_LDP*(RJ(ND)-TC(ND))/T2
	    IF(FIRST_FREQ)THEN
	      OPEN(UNIT=LU_JCOMP,STATUS='UNKNOWN',FILE='J_COMP')
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,'(A)')'Comparison of J at Outer and Inner',
	1       ' boundaries computed using Moments and Ray techniques.'
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,
	1       '(3X,A,14X,A,6X,A,6X,A,5X,A,5X,A,5X,A,7X,A,5X,A)')
	1       'Indx','Nu','J(mom)','J(ray)','%Diff','HBC_CMF',
	1       'J(mom)','J(ray)','%Diff'
	    END IF
	    WRITE(LU_JCOMP,'(I7,ES16.6,2ES12.4,F10.2,3ES12.4,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(ND),T2
	    IF(FREQ_INDX .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(FREQ_INDX .EQ. NCF)THEN
	    WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FLOAT(FG_COUNT)/FLOAT(NCF)
!
! edited by Devin, 27-Oct-2026: one-line ES-ALI summary per great
! iteration [DO_ES_ALI]: frequencies that exited with the last RJ
! update above ES_ALI_TOL, the worst final relative update, and the
! count of passes that fell back to plain Lambda.
!
	    IF(ES_ALI_ON)WRITE(LUER,'(A,I6,A,ES10.2,A,I6)')
	1     ' ES_ALI: freqs unconverged=',ES_ALI_NUNC,
	2     ', max final dJ/J=',ES_ALI_DJMAX,
	3     ', ALO fallbacks=',ES_ALI_NFALL
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=ABS(200.0_LDP*(RJ(1)-TC(1))/T1)
	    IF(T1 .GT. 100.0_LDP)THEN
	      WRITE(LUER,'(/,A)')' Error --- J(mom) and J(ray) differ by more than 1000% for last frequency'
	      WRITE(LUER,'(A)')' It is STRONGLY suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error and/or use plt_jh'
	    ELSE IF(T1 .GT. 50.0_LDP)THEN
	      WRITE(LUER,'(/,A)')' Error --- J(mom) and J(ray) differ by more than 50% for last frequency'
	      WRITE(LUER,'(A)')' It is strongly suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    ELSE IF(T1 .GT. 20.0_LDP)THEN
	      WRITE(LUER,'(/,A)')' Error --- J(mom) and J(ray) differ by more than 20% for last frequency'
	      WRITE(LUER,'(A)')' Although this is ulikely to effect the colution, it is suggested that'
	      WRITE(LUER,'(A)')' you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    END IF
	  END IF
	  CALL TUNE(ITWO,'CONT_VEL')
C
C 
C
	ELSE IF(EDDINGTON)THEN
C
C Calculation of "static" J in the continuum using Edington factors.
C
	  CALL TUNE(IONE,'JFEAU')
	  IF(SECTION .EQ. 'CONTINUUM' .AND. FREQ_INDX .EQ. 1)FEDD=0.0_LDP
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,ND
	      RJ(I)=0.0_LDP
              FOLD(I)=FEDD(I)
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),T1
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,'(/,A)')' Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
C
C JFNK_IN_EVAL soften -- same pattern as the other EDDFACTOR read
C sites: flag the eval broken and RETURN; the sweep epilogue restores
C the ACCESS_F cursor and forces the sentinel residual.  Production
C path unchanged (still STOPs).
C
	      IF(JFNK_IN_EVAL)THEN
	        JFNK_EVAL_BROKEN=.TRUE.
	        RETURN
	      END IF
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
	  DO WHILE(INACCURATE)
	    DO I=1,ND
	      SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	    END DO
	    S1=SOURCE(1)
	    CALL FQCOMP_IBC(TA,TB,TC,XM,DTAU,R,Z,P,QEDD,FEDD,
	1            SOURCE,CHI_CLUMP,dCHIdR,AQW,KQW,DBB,HBC_J,HBC_S,
	1            INBC,IC,THK_CONT,DIF,NC,ND,NP,METHOD)
	    CALL JFEAU_IBC(TA,TB,TC,DTAU,R,RJ,QEDD,FEDD,
	1          ZETA,THETA,CHI_CLUMP,DBB,IC,HBC_J,HBC_S,
	1          INBC,THK_CONT,DIF,ND,METHOD)
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 3 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0_LDP
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(J_IT_COUNTER .GT. 15)THEN
	         IF(FCONV_IT1_ACTIVE)THEN
	            FCONV_SUPPRESSED_CNT=FCONV_SUPPRESSED_CNT+1
	            FCONV_SUPPRESSED_MAX=MAX(FCONV_SUPPRESSED_MAX,T1)
	         ELSE
	            WRITE(LUER,*)'Possible error converging f - T1 is',T1
	            WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	         END IF
	      	 INACCURATE=.FALSE.
	      END IF
	
	    END DO
C
	    DO I=1,ND
	      K_MOM(I)=RJ(I)*FEDD(I)
	    END DO
C
C Output mean intensity for subsequent iterations.
C JFNK_IN_EVAL (CONTROL_VARIABLE_MOD): suppress during JFNK residual
C eval sweeps -- see note at the first LU_EDD write site.
C
	    IF(.NOT. JFNK_IN_EVAL)
	1     WRITE(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),FL
C
C Update record for next frequency
	    ACCESS_F=ACCESS_F+1
C
	  CALL TUNE(ITWO,'JFEAU')
	ELSE
C
C Calculation of "static" J in the continuum using Rybick method.
C
	  CALL TUNE(IONE,'JSOL')
	  CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1       ZETA,THETA,CHI_CLUMP,dCHIdR,AQW,
	1       THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
	  CALL TUNE(ITWO,'JSOL')
C
C Compute K_MOM.
C
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	  END DO
	  S1=SOURCE(1)
	  CALL FQCOMP_IBC(TA,TB,TC,XM,DTAU,R,Z,P,QEDD,FEDD,
	1            SOURCE,CHI_CLUMP,dCHIdR,AQW,KQW,DBB,HBC_J,HBC_S,
	1            INBC,IC,THK_CONT,DIF,NC,ND,NP,METHOD)
	  DO I=1,ND
	    K_MOM(I)=K_MOM(I)*FEDD(I)
	  END DO
	END IF
!
	IF(ABS(RJ(1)) .GT. 1.0E+30_LDP)THEN
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)'Error in comp_j_blank.f'
	  WRITE(LUER,*)'Mean intensity blowing up, which is due to an instabilty'
	  WRITE(LUER,*)'Try a finer grid at the outer boudary (preferred choice?)'
	  WRITE(LUER,*)'Alternatively try a different N_TYPE option'
!
! edited by Devin, 05-Oct-2026: under JFNK_IN_EVAL a blown-up RJ is a
! failed eval, not a fatal model error -- a STOP here kills the run
! mid-JVP.  Flag the eval (JFNK_EVAL_BROKEN) and return; the huge RJ
! was already computed and feeds the residual, and the sweep epilogue
! forces the sentinel so JFNK_RESIDUAL_V1's |Y|>1D30 cap rejects it.
! All per-frequency bookkeeping (ACCESS_F advance, suppressed LU_EDD
! write) already ran in the branch above, so RETURN loses nothing.
! Production path (JFNK_IN_EVAL=F) still STOPs -- bit-identical.
!
	  IF(JFNK_IN_EVAL)THEN
	    JFNK_EVAL_BROKEN=.TRUE.
	    RETURN
	  END IF
	  STOP
	END IF
!
	RETURN
	END
