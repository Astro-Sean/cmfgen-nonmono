C
C This subroutine computes the variation of J as a function of the emissivity
C ETA and the opacity CHI.
C
C TX has dimensions (ND,ND_SM,NM) and TVX has dimension (ND-1,ND_SM,NM).
C
C ND refelcts the number of depth points at which J is computed (
C      ==NDEXT in main code).
C ND_SM reflects the number of nodes at which the populations (and hence
C    opacities) are known (==ND in main code). This notation was adopted
C    to avoid unnecessary editing.
C NM reflects the total number of matrices required to compute the variation
C    of J. For this routine only the firts 2 are important and are used to
C    compute the variation of J with respect to the CURRENT opacity and the
C    CURRENT emissivity.
C
C The others are modified, but do not need to have any specific order.
C
C Typicall they will be CHI_C, and ETA_C, CHIL, AND ETAL of othe lines.
C
C NB TX(i,m, ) =dJ(i)/d(CHI(m),ETA(m),...)
C NB TVX(i,m, ) =dRSQH(i)/d(CHI(m),ETA(m),...)
C
C TX_DIF_d_T and TX_DIFF_d_dTdR are used to describe the variations in J
C caused by the DIFFUSION approxination at the inner boundary.
C
C The particular choice of the outer boundary condition adopted is irrelevant
C for this routine. Such information is incorporated by the outer boundary
C Eddington factors HBC, and NBC.
C
C Note TA, TBC, TC, HU etc have been defined so that the program computes
C JNU and r^2 HNU.
C
C This routine assumes that we are computing the radiation field at ND
C grid points. The radiation field is a function of the Emissivity and
C Opacity at these points. The opacity and emissivity at the transfer
C nodes have been derived from their values at ND_SM nodes. The
C populations are only computed on the ND_SM nodes. INDX and COEF are
C used to specify how the opacity and emissivities have been interpolated.
C
C NB:Compared to the calling routine:
C                                    ND=ND_EXT (calling)
C                                    ND_SM=ND (calling)
C
	SUBROUTINE VAR_MOM_JEXT_CMF_V3(ETA,CHI,ESEC,V,SIGMA,R,
	1                  ETA_SM,CHI_SM,ESEC_SM,
	1                  INDX,COEF,INTERP_TYPE,ND_SM,
	1                  TX,TVX,TX_DIF_d_T,TX_DIF_d_dTdR,
	1                  TVX_DIF_d_T,TVX_DIF_d_dTdR,
	1                  KI,RHS_dHdCHI,
	1                  INIT,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1                  DO_THIS_TX_MATRIX,METHOD,COHERENT,ND,NM,NM_KI)
	USE SET_KIND_MODULE
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
! Created 28-Dec-2008: Based on VAR_MOM_JEXT_CM_V1
!                      Inserted module MOD_RAY_MOM_STORE
!
	INTEGER ND
	INTEGER ND_SM
	INTEGER NM
	INTEGER NM_KI
C
C Emissivities etc on the transfer grid.
C
	REAL(KIND=LDP) ETA(ND),CHI(ND),ESEC(ND)
	REAL(KIND=LDP) V(ND),SIGMA(ND),R(ND)
C
C Emissivities etc on the (small) population grid.
C
	REAL(KIND=LDP) ETA_SM(ND_SM),CHI_SM(ND_SM),ESEC_SM(ND_SM)
	INTEGER INDX(ND)
	REAL(KIND=LDP) COEF(0:3,ND)
	CHARACTER*(*) INTERP_TYPE
C
C Variation arrays and vectors. NB: The variation is defined as the
C variation of J (or H) on the transfer grid with respect to the
C opacity/emissivity on the population (small) grid.
C
	REAL(KIND=LDP) TX(ND,ND_SM,NM)
	REAL(KIND=LDP) TVX(ND-1,ND_SM,NM)
C
C KI and RHS_dHdCHI are defined on the full grid, and then reduced to
C the smaller grid.
C
	REAL(KIND=LDP) KI(ND,ND,NM_KI)
	REAL(KIND=LDP) RHS_dHdCHI(ND-1,ND)
	REAL(KIND=LDP) TX_DIF_d_T(ND),TX_DIF_d_dTdR(ND)
	REAL(KIND=LDP) TVX_DIF_d_T(ND),TVX_DIF_d_dTdR(ND)
C
	LOGICAL DO_THIS_TX_MATRIX(NM)
C
	REAL(KIND=LDP) dLOG_NU,dTdR,DBB,dDBBdT,IC
	CHARACTER*6 METHOD
C
C INIT is used to indicate that there is no coupling to the previous frequency.
C We are thus solving the normal continuum transfer equation (i.e. the absence
C of velocity fields)
C
C COHERENT indicates whether the scattering is coherent. If it is, we
C explicitly take it into account. If COHERENT is FALSE, any electron
C scattering term should be included directly in the ETA that is passed
C to the routine.
C
	LOGICAL DIF,INIT,COHERENT
C
C Vectors required by future calls to VAR_MOM_J_CMF.
C
	INTEGER NV
	PARAMETER (NV=200)
	REAL(KIND=LDP) JNUM1(NV),RSQ_HNUM1(NV)
	SAVE JNUM1,RSQ_HNUM1
C
	REAL(KIND=LDP) KI_SM(ND,ND_SM,NM_KI)
	REAL(KIND=LDP) RHS_dHdCHI_SM(ND-1,ND_SM)
	REAL(KIND=LDP) WORKMAT(ND,ND)
C
C Work vectors.
C
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,DTAU,RSQ_DTAUONQ,
	1                   XM,SOURCE,Q,JNU,RSQ_HNU,
	1                   VB,VC,HU,HL,HS,THETA,
	1                   EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                   TX_OLD_d_T,TX_OLD_d_dTdR,
	1                   GAM,GAMH,W,WPREV,PSI,PSIPREV_MOD,PSIPREV
C
	REAL(KIND=LDP) TA(NV),TB(NV),TC(NV),DTAU(NV),RSQ_DTAUONQ(NV)
	REAL(KIND=LDP) XM(NV),SOURCE(NV),Q(NV),JNU(NV),RSQ_HNU(NV)
	REAL(KIND=LDP) VB(NV),VC(NV),HU(NV),HL(NV),HS(NV),THETA(NV)
	REAL(KIND=LDP) GAM(NV),GAMH(NV),W(NV),WPREV(NV)
	REAL(KIND=LDP) PSI(NV),PSIPREV_MOD(NV),PSIPREV(NV)
	REAL(KIND=LDP) EPS_A(NV),EPS_B(NV)
	REAL(KIND=LDP) EPS_PREV_A(NV),EPS_PREV_B(NV)
	REAL(KIND=LDP) TX_OLD_d_T(NV),TX_OLD_d_dTdR(NV)
C
	REAL(KIND=LDP) PROGDESC	
	REAL(KIND=LDP), PARAMETER :: PROG_ID=3.4281463E+08_LDP  !Must be unique (VAR_MOM_)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER I
	REAL(KIND=LDP) AV_SIGMA
C
C PROGDESC is a variable used to confirm that the scratch block is not
C being used by some other routine.
C
	PROGDESC=PROG_ID
	IF(ND .GT. NV)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in VAR_MOM_J_CMF - NV smaller than ND'
	  WRITE(I,*)'ND=',ND,'NV',NV
	  STOP
	END IF
C
C Zero common block. There are currently 29 vectors in the common block.
C TA must be the first vector, and PSIPREV the last.
C
	PSIPREV(NV-1)=1.0_LDP
	PSIPREV(NV)=1.0_LDP
	I=(NV*29)-1
	CALL DP_ZERO(TA,I)
	IF(PSIPREV(NV-1) .NE. 0.0_LDP .OR. PSIPREV(NV) .NE. 1.0_LDP)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in zeroing SCRATCH block in VAR_MOM_J_CMF'
	  STOP
	ELSE
	  PSIPREV(NV)=0.0_LDP
	END IF
C
C 
C
C Zero relevant vectors and matrices.
C
	JNU(1:ND)=0.0_LDP
	RSQ_HNU(1:ND)=0.0_LDP
C
C Allows us to indicate whether the electron scattering is assumed to
C to be COHERENT. If it is not coherent, we assume the ETA has already
C been updated (presumably assuming that the emission from scattering is
C from continuum photons only).
C
	IF(COHERENT)THEN
	  DO I=1,ND
	    THETA(I)=ESEC(I)/CHI(I)
	  END DO
	ELSE
	  THETA(1:ND)=0.0_LDP
	END IF
C
	IF(INIT)THEN
	  TX(:,:,:)=0.0_LDP			!ND*ND_SM*NM
	  TVX(:,:,:)=0.0_LDP			!(ND-1)*ND_SM*NM
	  DO I=1,ND
	    JNUM1(I)=0.0_LDP
	    RSQ_HNUM1(I)=0.0_LDP
	    GAMH(I)=0.0_LDP
	    GAM(I)=0.0_LDP
	    W(I)=0.0_LDP
	    WPREV(I)=0.0_LDP
	    PSI(I)=0.0_LDP
	    PSIPREV(I)=0.0_LDP
	    TX_DIF_d_T(I)=0.0_LDP
	    TX_DIF_d_dTdR(I)=0.0_LDP
	    EPS_A(I)=0.0_LDP
	    EPS_B(I)=0.0_LDP
	    EPS_PREV_A(I)=0.0_LDP
	    EPS_PREV_B(I)=0.0_LDP
 	  END DO
	END IF
C
	SOURCE(1:ND)=ETA(1:ND)/CHI(1:ND)
C
C Compute the Q factors from F. Then compute optical depth scale.
C
	CALL QFROMF(K_ON_J,Q,R,TA,TB,ND)	!TA work vector
	DO I=1,ND
	  TA(I)=CHI(I)*Q(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
C
C We need to call d_DERIVCHI_dCHI to set the TRAP derivatives.
C
	CALL d_DERIVCHI_dCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
	DO I=2,ND
	  RSQ_DTAUONQ(I)=0.5_LDP*R(I)*R(I)*(DTAU(I)+DTAU(I-1))/Q(I)
	END DO
C
C 
C
C Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
C NB - By definition, G is defined at the mesh midpoints.
C
C We evaluate and store the constant terms in the coputation of GAMH and
C GAM, since numer of operations only proprtional to ND. Late on scaing is
C proportional to NM*ND*ND.
C
	IF(.NOT. INIT)THEN
	  DO I=1,ND-1
	    AV_SIGMA=0.5_LDP*(SIGMA(I)+SIGMA(I+1))
	    GAMH(I)=2.0_LDP*3.33564E-06_LDP*(V(I)+V(I+1))/(R(I)+R(I+1))
	1         /dLOG_NU/( CHI(I)+CHI(I+1) )
	    W(I)=GAMH(I)*( 1.0_LDP+AV_SIGMA*NMID_ON_HMID(I) )
	    WPREV(I)=GAMH(I)*( 1.0_LDP+AV_SIGMA*NMID_ON_HMID(I) )
	    EPS_A(I)=GAMH(I)*AV_SIGMA*NMID_ON_J(I)/(1.0_LDP+W(I))
	    EPS_B(I)=EPS_A(I)*R(I+1)*R(I+1)
	    EPS_A(I)=EPS_A(I)*R(I)*R(I)
	    EPS_PREV_A(I)=GAMH(I)*AV_SIGMA*NMID_ON_J_PREV(I)/(1.0_LDP+W(I))
	    EPS_PREV_B(I)=EPS_PREV_A(I)*R(I+1)*R(I+1)
	    EPS_PREV_A(I)=EPS_PREV_A(I)*R(I)*R(I)
	  END DO
	  DO I=1,ND
	    GAM(I)=3.33564E-06_LDP*V(I)/R(I)/CHI(I)/dLOG_NU
	  END DO
C
C PSIPREV is equivalent to the U vector of FORMSOL.
C
	  PSI(1)=R(1)*R(1)*GAM(1)*( HBC+NBC*SIGMA(1) )
	  PSIPREV(1)=R(1)*R(1)*GAM(1)*( HBC_PREV+NBC_PREV*SIGMA(1) )
	  DO I=2,ND
	    PSI(I)=RSQ_DTAUONQ(I)*GAM(I)*( 1.0_LDP+SIGMA(I)*K_ON_J(I) )
	    PSIPREV(I)=RSQ_DTAUONQ(I)*GAM(I)*
	1                   ( 1.0_LDP+SIGMA(I)*K_ON_J_PREV(I) )
	  END DO
	END IF
C
C 
C
C Compute vectors used to compute the flux vector H.
C
	DO I=1,ND-1
	  HU(I)=R(I+1)*R(I+1)*K_ON_J(I+1)*Q(I+1)/(1.0_LDP+W(I))/DTAU(I)
	  HL(I)=R(I)*R(I)*K_ON_J(I)*Q(I)/(1.0_LDP+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0_LDP+W(I))
	END DO
C
C Compute the TRIDIAGONAL operators, and the RHS source vector.
C
	DO I=2,ND-1
	  TA(I)=-HL(I-1)-EPS_A(I-1)
	  TC(I)=-HU(I)+EPS_B(I)
	  TB(I)=RSQ_DTAUONQ(I)*(1.0_LDP-THETA(I)) + PSI(I) +HU(I-1) +HL(I)
	1             -EPS_B(I-1)+EPS_A(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=RSQ_DTAUONQ(I)*SOURCE(I)
	END DO
C
C Evaluate TA,TB,TC for boudary conditions
C
	TC(1)=-R(2)*R(2)*K_ON_J(2)*Q(2)/DTAU(1)
	TB(1)=R(1)*R(1)*( K_ON_J(1)*Q(1)/DTAU(1) + HBC ) + PSI(1)
	XM(1)=0.0_LDP
	TA(1)=0.0_LDP
	VB(1)=0.0_LDP
	VC(1)=0.0_LDP
C
	TA(ND)=-R(ND-1)*R(ND-1)*K_ON_J(ND-1)*Q(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  TB(ND)=R(ND)*R(ND)*K_ON_J(ND)/DTAU(ND-1)
	  XM(ND)=DBB*R(ND)*R(ND)/3.0_LDP/CHI(ND)
	ELSE
	  TB(ND)=R(ND)*R(ND)*K_ON_J(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=R(ND)*R(ND)*IC*(0.25_LDP+0.5_LDP*IN_HBC)
	END IF
	TC(ND)=0.0_LDP
	VB(ND)=0.0_LDP
	VC(ND)=0.0_LDP
	PSIPREV(ND)=0.0_LDP
C
C We create PSIPREV_MOD to save multiplications in the UP_TX_TVX routine/
C It is only different from PSIPREV when N_ON_J is non zero.
C
	PSIPREV_MOD(1)=PSIPREV(1)
	PSIPREV_MOD(ND)=PSIPREV(ND)
	DO I=2,ND-1
	  PSIPREV_MOD(I)=(EPS_PREV_A(I)-EPS_PREV_B(I-1)) + PSIPREV(I)
	END DO
C
	XM(1)=XM(1) + PSIPREV_MOD(1)*JNUM1(1)
	DO I=2,ND-1
	  XM(I)=XM(I) + VB(I)*RSQ_HNUM1(I-1) + VC(I)*RSQ_HNUM1(I)
	1          + PSIPREV_MOD(I)*JNUM1(I)
	1          + ( EPS_PREV_B(I)*JNUM1(I+1)
	1               - EPS_PREV_A(I-1)*JNUM1(I-1) )
	END DO
	XM(ND)=XM(ND) + PSIPREV_MOD(ND)*JNUM1(ND)
C
C Solve for the radiation field along ray for this frequency.
C
	CALL THOMAS(TA,TB,TC,XM,ND,1)
C
	DO I=1,ND
	  JNU(I)=XM(I)
	END DO
C
	DO I=1,ND-1
	  RSQ_HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQ_HNUM1(I) +
	1              (EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*XM(I)) +
	1              (EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*XM(I+1))
	END DO
C 
C
C The J and H components of the radiation field have been found. We can thus
C begin the variation computation.
C
C Compute d{non-radiation field}/dchi matrix.
C
	CALL TUNE(1,'MOM_EDD')
	I=2
	CALL EDD_J_VAR_V4(KI,RHS_dHdCHI,WORKMAT,
	1                SOURCE,CHI,ESEC,DTAU,R,SIGMA,
	1                K_ON_J,Q,HU,HL,HS,RSQ_DTAUONQ,
	1                W,WPREV,PSI,PSIPREV,
	1                EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                JNU,JNUM1,RSQ_HNUM1,
	1                DBB,DIF,COHERENT,ND,I)
	CALL TUNE(2,'MOM_EDD')
C
C Now need to remove variation on the values of chi, eta at the
C inserted nodes.
C
	CALL FIX_dCHI(KI_SM,RHS_dHdCHI_SM,CHI_SM,ETA_SM,ND_SM,
	1                 KI,RHS_dHdCHI,CHI,ETA,ND,
	1                 COEF,INDX,INTERP_TYPE)
C
C Evaluate the intensity variations.
C                                  TX=dJ/d(chi,eta,....)
C                          and     TVX=dRSQH/d(chi,eta,...)
C
C WORKMAT is dimension (ND,ND_SM) is is used to temporaily save
c TX( , ,K) for each K.
C
	CALL TUNE(1,'UP_TX')
	CALL UP_TX_TVX_EXT_V1(TX,TVX,KI_SM,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI_SM,
	1                       EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                       WORKMAT,ND,ND_SM,NM,NM_KI,
	1                       INIT,DO_THIS_TX_MATRIX)
	CALL TUNE(2,'UP_TX')
C
C 
C
C Evaluate diffusion variation. TXD is initially the value from the
C previous frquency.
C
	DO I=1,ND
	  TX_OLD_d_T(I)=TX_DIF_d_T(I)
	  TX_OLD_d_dTdR(I)=TX_DIF_d_dTDR(I)
	END DO
	IF(DIF)THEN
	  TX_DIF_d_T(1)=PSIPREV(1)*TX_DIF_d_T(1)
	  TX_DIF_d_dTdR(1)=PSIPREV(1)*TX_DIF_d_dTdR(1)
	  DO I=2,ND-1
	    TX_DIF_d_T(I)= PSIPREV_MOD(I)*TX_DIF_d_T(I)
	1             + VB(I)*TVX_DIF_d_T(I-1) + VC(I)*TVX_DIF_d_T(I)
	1          + ( EPS_PREV_B(I)*TX_OLD_d_T(I+1)
	1               - EPS_PREV_A(I-1)*TX_OLD_d_T(I-1) )
	    TX_DIF_d_dTdR(I)= PSIPREV_MOD(I)*TX_DIF_d_dTdR(I)
	1             + VB(I)*TVX_DIF_d_dTdR(I-1) + VC(I)*TVX_DIF_d_dTdR(I)
	1          + ( EPS_PREV_B(I)*TX_OLD_d_dTdR(I+1)
	1               - EPS_PREV_A(I-1)*TX_OLD_d_dTdR(I-1) )
	  END DO
	  TX_DIF_d_T(ND)=dDBBdT*R(ND)*R(ND)/3.0_LDP/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_T(ND)
	  TX_DIF_d_dTdR(ND)=DBB/dTdR*R(ND)*R(ND)/3.0_LDP/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_dTdR(ND)
C
C Solve for the radiation field along ray for this frequency.
C
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_T,ND,1)
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_dTdR,ND,1)
C
	  DO I=1,ND-1
	    TVX_DIF_d_T(I)=HU(I)*TX_DIF_d_T(I+1) - HL(I)*TX_DIF_d_T(I) +
	1        HS(I)*TVX_DIF_d_T(I) +
	1        (EPS_PREV_A(I)*TX_OLD_d_T(I)-EPS_A(I)*TX_DIF_d_T(I)) +
	1        (EPS_PREV_B(I)*TX_OLD_d_T(I+1)-EPS_B(I)*TX_DIF_d_T(I+1))
	    TVX_DIF_d_dTdR(I)=HU(I)*TX_DIF_d_dTdR(I+1) -
	1     HL(I)*TX_DIF_d_dTdR(I) +
	1     HS(I)*TVX_DIF_d_dTdR(I) +
	1     (EPS_PREV_A(I)*TX_OLD_d_dTdR(I)-EPS_A(I)*TX_DIF_d_dTdR(I)) +
	1     (EPS_PREV_B(I)*TX_OLD_d_dTdR(I+1)-EPS_B(I)*TX_DIF_d_dTdR(I+1))
	  END DO
C
	END IF	    	    !DIF END
C
C 
C
C Save JNU and RSQ_HNU for next frequency integration.
C
	JNUM1(1:ND)=JNU(1:ND)
	RSQ_HNUM1(1:ND)=RSQ_HNU(1:ND)
C
C 
C
	IF(PROGDESC .NE. PROG_ID)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error - SCRATCH block corrupted in VAR_MOM_J_CMF'
	  STOP
	END IF
C
	RETURN
	END
