
! Data module for PP_MOM_CMF_V1. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE PP_MOM_CMF_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! To be dimenensioned ND_SM where ND_SM is the size of the R grid
! as passed to MOM_J_CMF.
!
	REAL(KIND=LDP), ALLOCATABLE :: LOG_R_SM(:)
!
! Dimensioned ND_SM,4
!
	REAL(KIND=LDP), ALLOCATABLE :: V_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: F_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: G_COEF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: N_ON_J_COEF(:,:)
!
! All the following rays have dimension ND, where ND >= ND_SM.
! Some of the data in the arrays is need in subsequent calls.
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: JNU(:)
	REAL(KIND=LDP), ALLOCATABLE :: HNU(:)
	REAL(KIND=LDP), ALLOCATABLE :: F(:)
	REAL(KIND=LDP), ALLOCATABLE :: G(:)
	REAL(KIND=LDP), ALLOCATABLE :: N_ON_J(:)
	REAL(KIND=LDP), ALLOCATABLE :: F_SAV(:)
	REAL(KIND=LDP), ALLOCATABLE :: G_SAV(:)
	REAL(KIND=LDP), ALLOCATABLE :: N_ON_J_SAV(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: JNU_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: HNU_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: F_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: G_PREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: N_ON_J_PREV(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: CON_GAM(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_GAMH(:)
	REAL(KIND=LDP), ALLOCATABLE :: AV_SIGMA(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: DTAU(:)
	REAL(KIND=LDP), ALLOCATABLE :: DTAU_MID(:)
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
	REAL(KIND=LDP), ALLOCATABLE :: TB(:)
	REAL(KIND=LDP), ALLOCATABLE :: TC(:)
	REAL(KIND=LDP), ALLOCATABLE :: XM(:)
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE(:)
	REAL(KIND=LDP), ALLOCATABLE :: VB(:)
	REAL(KIND=LDP), ALLOCATABLE :: VC(:)
	REAL(KIND=LDP), ALLOCATABLE :: HU(:)
	REAL(KIND=LDP), ALLOCATABLE :: HL(:)
	REAL(KIND=LDP), ALLOCATABLE :: HS(:)
	REAL(KIND=LDP), ALLOCATABLE :: COH_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAMH(:)
	REAL(KIND=LDP), ALLOCATABLE :: W(:)
	REAL(KIND=LDP), ALLOCATABLE :: WPREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: PSI(:)
	REAL(KIND=LDP), ALLOCATABLE :: PSIPREV(:)
	REAL(KIND=LDP), ALLOCATABLE :: EPS(:)
	REAL(KIND=LDP), ALLOCATABLE :: EPS_PREV(:)
!
! Boundary conditions
!
        REAL(KIND=LDP) HBC_PREV
	REAL(KIND=LDP) IN_HBC_PREV
	REAL(KIND=LDP) NBC_PREV
	REAL(KIND=LDP) NBC_INCID_PREV
!
	REAL(KIND=LDP) HBC_SAV
	REAL(KIND=LDP) IN_HBC_SAV
	REAL(KIND=LDP) NBC_SAV
	REAL(KIND=LDP) NBC_INCID_SAV
	REAL(KIND=LDP) VDOP_FRAC_SAV
!
	INTEGER ND
!
! R_PNT(K) defines the interpolation for the variable at depth K.
!
	INTEGER, ALLOCATABLE :: R_PNT(:)
!
! ?_INDX are used to indicate the location of J and H on the small grid
! in the larger array.
!
	INTEGER, ALLOCATABLE :: J_INDX(:)
	INTEGER, ALLOCATABLE :: H_INDX(:)
!
	LOGICAL FIRST_TIME
	DATA FIRST_TIME/.TRUE./
	DATA VDOP_FRAC_SAV/-10001.1_LDP/    !Absurd value
!
	END MODULE PP_MOM_CMF_MOD
!
!
!
!
! Routine to compute the mean intensity J at a single frequency in the
! Comoving-Frame. The computed intensity thus depends on the intensity
! computed for the previous (bluer) frequency. Routine is designed for
! a plane-parallel atomoshere.
!
! The F, G, and N_ON_J Eddingto factors must be supplied.
!
! NB:
!	F = K / J
!	G=  N / H
!	N_ON_J(I) = N(I)/(J(I)+J(I+1))
!
! NB: Only one of G and N_ON_J is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
!       N_ON_J=0 at all depths.
!     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) N_ON_J is defined at all
!       depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or N_ON_J is
!       non-zero, and is the value to be used in MOM_J_CMF
!
	SUBROUTINE PP_MOM_CMF_V1(ETA_SM,CHI_SM,ESEC_SM,
	1                  V_SM,SIGMA_SM,R_SM,
	1		   F_SM,G_SM,N_ON_J_SM,
	1                  JNU_SM,HNU_SM,
	1                  VDOP_VEC,VDOP_FRAC,
	1                  IN_HBC,HBC,HBC_INCID,NBC,NBC_INCID,
	1                  FREQ,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,COHERENT,
	1                  INIT,NEW_FREQ,ND_SM)
	USE SET_KIND_MODULE
	USE PP_MOM_CMF_MOD
	IMPLICIT NONE
!
! Altered:   17-Feb-2007 : dNBC_INCID was not zeroed on when INIT set.
! Created:   02-Mar-2006 : Based on MOM_J_CMF_V6
!
	INTEGER ND_SM
	REAL(KIND=LDP) ETA_SM(ND_SM)
	REAL(KIND=LDP) CHI_SM(ND_SM)
	REAL(KIND=LDP) ESEC_SM(ND_SM)
	REAL(KIND=LDP) V_SM(ND_SM)
	REAL(KIND=LDP) SIGMA_SM(ND_SM)
	REAL(KIND=LDP) R_SM(ND_SM)
!
! Radiation field variables. F, G, JNU_PREV, and HNU_PREV must be supplied.
! JNU and HNU recomputed.
!
	REAL(KIND=LDP) F_SM(ND_SM)
	REAL(KIND=LDP) G_SM(ND_SM)
	REAL(KIND=LDP) N_ON_J_SM(ND_SM)
	REAL(KIND=LDP) JNU_SM(ND_SM)
	REAL(KIND=LDP) HNU_SM(ND_SM)
	REAL(KIND=LDP) VDOP_VEC(ND_SM)
	REAL(KIND=LDP) VDOP_FRAC
!
	INTEGER N_ERR_MAX,MOM_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL(KIND=LDP) MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
! Boundary conditions.
!
	REAL(KIND=LDP) HBC,HBC_INCID
	REAL(KIND=LDP) NBC,NBC_INCID
	REAL(KIND=LDP) IN_HBC
!
	REAL(KIND=LDP) DBB,IC,FREQ,dLOG_NU
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
! NEW_FREQ is used to indicatae that we are computing J for a new frequency.
! If we were iterating between computing J and the Eddington factors, NEW_FREQ
! would be set to false.
!
! COHERENT indicates whether the scattering is coherent. If it is, we
! explicitly take it into account. If COHERENT is FALSE, any electron
! scattering term should be included directly in the ETA that is passed
! to the routine.
!
	LOGICAL DIF,INIT,COHERENT,NEW_FREQ
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP)  DELTA_R
	REAL(KIND=LDP)  dNBC_INCID
	INTEGER IT1,I,J,K
	INTEGER IOS
	LOGICAL NEW_R_GRID
!
!
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. ALLOCATED(R))THEN
	  DO I=1,ND_SM
	    IF(LOG(R_SM(I)) .NE. LOG_R_SM(I))THEN
	      NEW_R_GRID=.TRUE.
	      WRITE(171,*)'Updating RGRID in MOM_J_CMF_V6'
	      EXIT
	    END IF
	  END DO
          IF(VDOP_FRAC .NE. VDOP_FRAC_SAV)NEW_R_GRID=.TRUE.
	END IF
	IF(FIRST_TIME)THEN
          OPEN(UNIT=47,FILE='MOM_J_ERRORS',STATUS='UNKNOWN')
        ELSE IF(INIT)THEN
	  REWIND(47)
	END IF
!
! Deallocate all arrayes if we have changed VDOP_FRAC. This will only
! be done in testing this routine (e.g., using DISPGEN).
!
	IF(ALLOCATED(R) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R, R_PNT, LOG_R_SM )
	  DEALLOCATE ( V_COEF, SIGMA_COEF, ETA_COEF, ESEC_COEF, CHI_COEF )
	  DEALLOCATE ( F_COEF, G_COEF, N_ON_J_COEF )
	  DEALLOCATE ( V, SIGMA, ETA, ESEC, CHI )
	  DEALLOCATE ( JNU, HNU, F, G, N_ON_J )
	  DEALLOCATE ( F_SAV, G_SAV, N_ON_J_SAV )
	  DEALLOCATE ( JNU_PREV, HNU_PREV, F_PREV, G_PREV, N_ON_J_PREV )
	  DEALLOCATE ( CON_GAM, CON_GAMH, AV_SIGMA )
	  DEALLOCATE ( DTAU, DTAU_MID, TA, TB, TC, XM, SOURCE )
	  DEALLOCATE ( VB, VC, HU, HL, HS, COH_VEC )
	  DEALLOCATE ( GAM, GAMH, W, WPREV, PSI, PSIPREV )
	  DEALLOCATE ( EPS, EPS_PREV, J_INDX, H_INDX )
	END IF
	VDOP_FRAC_SAV=VDOP_FRAC
!
! On the very first entry, we define the improved R grid, and allocate all
! data arrays.
!
	IF( FIRST_TIME .OR. .NOT. ALLOCATED(R) )THEN
!
! Determine the number of points for the expanded R grid.
! We always insert an EVEN number of points. This guarentees that
! H_SM (defined at the midpoints of the pass grid) has an exact correspondence
! with H defined on the extended gid.
!
          K=1
          T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:ND_SM))
          DO I=1,ND_SM-1
            IT1=INT( (V_SM(I)-V_SM(I+1))/T2 )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
            IF(IT1 .GT. 0)K=K+IT1
            K=K+1
          END DO
          ND=K
!
	  ALLOCATE ( R(ND), STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( R_PNT(ND), STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Unable to allocate R & R_PNT in PP_MOM_CMF_V1'
	    WRITE(LUER,*)'Status=',IOS
	  END IF
!
          K=1
	  R(1)=R_SM(1)
          R_PNT(1)=1
	  DO I=1,ND_SM-1
            IT1=INT( (V_SM(I)-V_SM(I+1))/T2 )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
            IF(IT1 .GT. 0)THEN
              DELTA_R=(R_SM(I+1)-R_SM(I))/(IT1+1)
              DO J=1,IT1
                K=K+1
                R(K)=R(K-1)+DELTA_R
                R_PNT(K)=I
              END DO
            END IF
            K=K+1
            R(K)=R_SM(I+1)
            R_PNT(K)=I
	  END DO
!
	  ALLOCATE ( LOG_R_SM(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (V_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (SIGMA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ETA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ESEC_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (G_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (N_ON_J_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(lUER,*)'Unable to allocate COEF memory in PP_MOM_CMF_V1'
	     STOP
	  END IF
!
	  ALLOCATE ( V(ND), SIGMA(ND), ETA(ND), ESEC(ND), CHI(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (JNU(ND), HNU(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F(ND), G(ND), N_ON_J(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F_SAV(ND), G_SAV(ND), N_ON_J_SAV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (JNU_PREV(ND), HNU_PREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F_PREV(ND), G_PREV(ND), N_ON_J_PREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_GAM(ND), CON_GAMH(ND), AV_SIGMA(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAU(ND), DTAU_MID(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TA(ND),  TB(ND),  TC(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (XM(ND), SOURCE(ND), COH_VEC(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VB(ND), VC(ND), HU(ND), HL(ND), HS(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM(ND), GAMH(ND), W(ND), WPREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PSI(ND), PSIPREV(ND), EPS(ND), EPS_PREV(ND), STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(lUER,*)'Unable to allocate vectors in PP_MOM_CMF_V1'
	     STOP
	  END IF
!
!
	  LOG_R_SM(1:ND_SM)=LOG(R_SM(1:ND_SM))
          CALL MON_INT_FUNS_V2(V_COEF,V_SM,LOG_R_SM,ND_SM)
          CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_SM,LOG_R_SM,ND_SM)
          DO I=1,ND
            K=R_PNT(I)
            T1=LOG(R(I)/R_SM(K))
            V(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
            SIGMA(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
          END DO
!
	  ALLOCATE ( J_INDX(ND_SM),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( H_INDX(ND_SM),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Unable to allocate J_INDX & H_INDX in PP_MOM_CMF_V1'
	    WRITE(LUER,*)'Status=',IOS
	  END IF
	  J_INDX(1:ND_SM)=0;    H_INDX(1:ND_SM)=0
!
	  K=1
	  DO I=1,ND_SM
	    DO WHILE(J_INDX(I) .EQ. 0)
	      IF(R_SM(I) .LE. R(K) .AND. R_SM(I) .GE. R(K+1))THEN
	        IF( (R(K)-R_SM(I)) .LT. (R_SM(I)-R(K+1)) )THEN
	          J_INDX(I)=K
	        ELSE
	          J_INDX(I)=K+1
	        END IF
	      ELSE
	        K=K+1
	      END IF
	    END DO
	  END DO
!
	  K=1
	  DO I=1,ND_SM-1
	    T1=0.5_LDP*(R_SM(I)+R_SM(I+1))
	    DO WHILE(H_INDX(I) .EQ. 0)
	      IF(T1 .LT. R(K) .AND. T1 .GT. R(K+1))THEN
	        H_INDX(I)=K
	      ELSE
	        K=K+1
	      END IF
	    END DO
	  END DO
!
	  IF(N_TYPE .NE. 'G_ONLY' .AND. N_TYPE .NE. 'N_ON_J')THEN
	    IF(ND .NE. ND_SM)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in MOM_J_CMF_V6'
	      WRITE(LUER,*)'Cannot use N_TYPE=''MIXED'' when inserting extra points'
	    END IF
	  END IF
!
	  FIRST_TIME=.FALSE.
	END IF
!
!
	IF(ND .GT. ND_SM)THEN
!
! Interpolate quantities onto revised grid.
!
	  TA(1:ND_SM)=LOG(CHI_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(CHI_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ESEC_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ESEC_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ETA_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ETA_COEF,TA,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(F_COEF,F_SM,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(G_COEF,G_SM,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(N_ON_J_COEF,N_ON_J_SM,LOG_R_SM,ND_SM)

	  DO I=1,ND
	    K=R_PNT(I)
	    T1=LOG(R(I)/R_SM(K))
	    T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	    CHI(I)=EXP(T2)
	    T2=((ESEC_COEF(K,1)*T1+ESEC_COEF(K,2))*T1+ESEC_COEF(K,3))*T1+ESEC_COEF(K,4)
	    ESEC(I)=EXP(T2)
	    T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	    ETA(I)=EXP(T2)
	    F(I)=((F_COEF(K,1)*T1+F_COEF(K,2))*T1+F_COEF(K,3))*T1+F_COEF(K,4)
	    G(I)=((G_COEF(K,1)*T1+G_COEF(K,2))*T1+G_COEF(K,3))*T1+G_COEF(K,4)
	    N_ON_J(I)=((N_ON_J_COEF(K,1)*T1+N_ON_J_COEF(K,2))*T1+
	1                               N_ON_J_COEF(K,3))*T1+N_ON_J_COEF(K,4)
	  END DO
!
	ELSE
	  ESEC(1:ND)=ESEC_SM(1:ND)
	  CHI(1:ND)=CHI_SM(1:ND)
	  ETA(1:ND)=ETA_SM(1:ND)
	  F(1:ND)=F_SM(1:ND)
	  G(1:ND)=G_SM(1:ND)
	  N_ON_J(1:ND)=N_ON_J_SM(1:ND)
	END IF
!
! 
!

	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0_LDP
	  END DO
	  MOM_ERR_CNT=0
          JNU=0.0_LDP; HNU=0.0_LDP
	  JNU_PREV=0.0_LDP; HNU_PREV=0.0_LDP
	  N_ON_J_SAV=0.0_LDP; G_SAV=0.0_LDP; F_SAV=0.0_LDP
	END IF
!
!*****************************************************************************
!
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
	IF(COHERENT)THEN
	  DO I=1,ND
	    COH_VEC(I)=ESEC(I)/CHI(I)
	  END DO
	ELSE
	  DO I=1,ND
	    COH_VEC(I)=0.0_LDP
	  END DO
	END IF
!
	DO I=1,ND
	  IF(G(I) .GT. 1.0_LDP)G(I)=1.0_LDP
	  IF(G(I) .LT. -1.0_LDP)G(I)=-1.0_LDP
	END DO
!
! NB: We solve for J.
!
!
	CALL DERIVCHI(TB,CHI,R,ND,METHOD)
	CALL NORDTAU(DTAU,CHI,R,R,TB,ND)
!
	IF(NEW_FREQ)THEN
	  F_PREV(1:ND)=F_SAV(1:ND)
	  G_PREV(1:ND)=G_SAV(1:ND)
	  N_ON_J_PREV(1:ND)=N_ON_J_SAV(1:ND)
	  JNU_PREV(1:ND)=JNU(1:ND)
	  HNU_PREV(1:ND)=HNU(1:ND)
	  HBC_PREV=HBC_SAV;  IN_HBC_PREV=IN_HBC_SAV
	  NBC_PREV=NBC_SAV
	  NBC_INCID_PREV=NBC_INCID_SAV
	END IF
!
!
	IF(INIT)THEN
	  DO I=1,ND
	    GAMH(I)=0.0_LDP
	    GAM(I)=0.0_LDP
	    W(I)=0.0_LDP
	    WPREV(I)=0.0_LDP
	    PSI(I)=0.0_LDP
	    PSIPREV(I)=0.0_LDP
	    JNU_PREV(I)=0.0_LDP
	    HNU_PREV(I)=0.0_LDP
	    EPS(I)=0.0_LDP
	    EPS_PREV(I)=0.0_LDP
	  END DO
	  HBC_PREV=0.0_LDP;  IN_HBC_PREV=0.0_LDP; NBC_PREV=0.0_LDP
	  dNBC_INCID=0.0_LDP
	ELSE
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	  DO I=1,ND-1
	    CON_GAMH(I)=2.0_LDP*3.33564E-06_LDP*(V(I)+V(I+1))/(R(I)+R(I+1))
	    AV_SIGMA(I)=0.5_LDP*(SIGMA(I)+SIGMA(I+1))
	    CON_GAM(I)=3.33564E-06_LDP*V(I)/R(I)
	  END DO
	  CON_GAM(ND)=3.33564E-06_LDP*V(ND)/R(ND)
!
! Since we are intgerating from blue to red, FL_PREV is always larger than
! FL. dLOG_NU is define as vd / dv which is the same as d / d ln v.
!
! EPS is used if we define N in terms of J rather than H, This is sometimes
! useful as H can approach zero, and hence N/H is undefined.
!
	  IF(N_TYPE .EQ. 'G_ONLY')THEN
	    DO I=1,ND-1
	      GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	      W(I)=GAMH(I)*(1.0_LDP+AV_SIGMA(I))*G(I)
	      WPREV(I)=GAMH(I)*(1.0_LDP+AV_SIGMA(I))*G_PREV(I)
	    END DO
	  ELSE
	    DO I=1,ND-1
	      GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	      W(I)=GAMH(I)*(1.0_LDP+AV_SIGMA(I))*G(I)
	      WPREV(I)=GAMH(I)*(1.0_LDP+AV_SIGMA(I))*G_PREV(I)
	      EPS(I)=GAMH(I)*(1.0_LDP+AV_SIGMA(I))*N_ON_J(I)/(1.0_LDP+W(I))
	      EPS_PREV(I)=GAMH(I)*(1.0_LDP+AV_SIGMA(I))*N_ON_J_PREV(I)/(1.0_LDP+W(I))
	    END DO
	  END IF
!
	  DO I=1,ND
	    GAM(I)=CON_GAM(I)/CHI(I)/dLOG_NU
	  END DO
!
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	  PSI(1)=GAM(1)*NBC*(1.0_LDP+SIGMA(1))
	  PSIPREV(1)=GAM(1)*NBC_PREV*(1.0_LDP+SIGMA(1))
	  dNBC_INCID=GAM(1)*(1.0_LDP+SIGMA(1))*(NBC_INCID-NBC_INCID_PREV)
	END IF
!
! 
!
	DO I=2,ND-1
	  DTAU_MID(I)=0.5_LDP*(DTAU(I)+DTAU(I-1))
	  PSI(I)=DTAU_MID(I)*GAM(I)*(1.0_LDP+SIGMA(I))*F(I)
	  PSIPREV(I)=DTAU_MID(I)*GAM(I)*(1.0_LDP+SIGMA(I))*F_PREV(I)
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=F(I+1)/(1.0_LDP+W(I))/DTAU(I)
	  HL(I)=F(I)/(1.0_LDP+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0_LDP+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)
	    TC(I)=-HU(I)
	    TB(I)=DTAU_MID(I)*(1.0_LDP-COH_VEC(I)) + PSI(I) + HL(I) + HU(I-1)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAU_MID(I)*SOURCE(I)
	  END DO
	ELSE
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)-EPS(I-1)
	    TC(I)=-HU(I)+EPS(I)
	    TB(I)=DTAU_MID(I)*(1.0_LDP-COH_VEC(I)) + PSI(I) + HL(I) + HU(I-1)
	1               -EPS(I-1)+EPS(I)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAU_MID(I)*SOURCE(I)
	  END DO
	END IF
!
! Evaluate TA,TB,TC for boudary conditions
!
	TC(1)=-F(2)/DTAU(1)
	TB(1)=F(1)/DTAU(1) + PSI(1) + HBC
	XM(1)=HBC_INCID+dNBC_INCID
	TA(1)=0.0_LDP
	VB(1)=0.0_LDP
	VC(1)=0.0_LDP
!
	TA(ND)=-F(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  TB(ND)=F(ND)/DTAU(ND-1)
	  XM(ND)=DBB/3.0_LDP/CHI(ND)
	ELSE
	  TB(ND)=F(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=IC*(0.25_LDP+0.5_LDP*IN_HBC)
	END IF
	TC(ND)=0.0_LDP
	VB(ND)=0.0_LDP
	VC(ND)=0.0_LDP
	PSIPREV(ND)=0.0_LDP
!
! Note that EPS and EPS_PREV will be identically zero hence when N_TYPE is
! G_ONLY.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  XM(1)=XM(1) + PSIPREV(1)*JNU_PREV(1)
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*HNU_PREV(I-1) + VC(I)*HNU_PREV(I)
	1          + PSIPREV(I)*JNU_PREV(I)
	  END DO
	  XM(ND)=XM(ND)
	ELSE
	  XM(1)=XM(1) + PSIPREV(1)*JNU_PREV(1)
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*HNU_PREV(I-1) + VC(I)*HNU_PREV(I)
	1          + PSIPREV(I)*JNU_PREV(I)
	1          - EPS_PREV(I-1)*(JNU_PREV(I-1)+JNU_PREV(I))
	1          + EPS_PREV(I)*(JNU_PREV(I)+JNU_PREV(I+1))
	  END DO
	  XM(ND)=XM(ND)
	END IF
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
! Check that no negative mean intensities have been computed.
!
	IF(MINVAL(XM(1:ND)) .LE. 0.0_LDP)THEN
	   WRITE(47,*)'Freq=',FREQ
	   TA(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
	   CALL WRITV(TA,ND,'XM Vec',47)
	   CALL WRITV(F,ND,'F Vec',47)
	   CALL WRITV(G,ND,'G Vec',47)
	   CALL WRITV(ETA,ND,'ETA Vec',47)
	   CALL WRITV(ESEC,ND,'ESEC Vec',47)
	   CALL WRITV(CHI,ND,'CHI Vec',47)
	END IF
!
	RECORDED_ERROR=.FALSE.
	DO I=1,ND
	  IF(XM(I) .LT. 0.0_LDP)THEN
	    XM(I)=ABS(XM(I))/10.0_LDP
	  END IF
	  IF(.NOT. RECORDED_ERROR)THEN
	    IF(MOM_ERR_CNT .GT. N_ERR_MAX)THEN
	      MOM_ERR_CNT=MOM_ERR_CNT+1
	    ELSE IF(MOM_ERR_ON_FREQ(MOM_ERR_CNT) .NE. FREQ)THEN
	      MOM_ERR_CNT=MOM_ERR_CNT+1
	      IF(MOM_ERR_CNT .LT. N_ERR_MAX)MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	    END IF
	    RECORDED_ERROR=.TRUE.
	  END IF
	END DO
!
! Save R^2 J for next iteration.
!
	JNU(1:ND)=XM(1:ND)
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=1,ND-1
	    HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*HNU_PREV(I)
	  END DO
	ELSE
	  DO I=1,ND-1
	    HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*HNU_PREV(I)+
	1              ( EPS_PREV(I)*(JNU_PREV(I)+JNU_PREV(I+1)) -
	1                  EPS(I)*(XM(I)+XM(I+1)) )
	  END DO
	END IF
!
! Regrid derived J and H values onto small grid. We devide J by R^2 so that
! we return J.
!
	DO I=1,ND_SM
	  K=J_INDX(I)
	  JNU_SM(I)=JNU(K)
	END DO
!
	DO I=1,ND_SM-1
	  K=H_INDX(I)
	  HNU_SM(I)=HNU(K)
	END DO
!
	F_SAV(1:ND)=F(1:ND)
	G_SAV(1:ND)=G(1:ND)
	N_ON_J_SAV(1:ND)=N_ON_J(1:ND)
	HBC_SAV=HBC
	IN_HBC_SAV=IN_HBC
	NBC_SAV=NBC
	NBC_INCID_SAV=NBC_INCID
!
	RETURN
	END
