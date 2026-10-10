!
! Subroutine to increment the Linearization of Electron Energy Balance equation
! for terms that depend directly on the intensity J. The Radiative equilibrium equation
! is not altered.
!
	SUBROUTINE VEHB_BYJ_V1(ID,
	1             WSE,dWSEdT,WCR,dWCRdT,
	1             HN,HNST,dlnHNST_dlnT,NLEV,EQGS,
	1             DI,LOG_DIST,dlnDIST_dlnT,N_DI,ION_LEV,
	1             ED,T,
	1             JREC,dJRECdT,JPHOT,
	1             JREC_CR,dJREC_CRdT,JPHOT_CR,
	1             FIXED_T,ND,DST,DEND,NT)
	USE SET_KIND_MODULE
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Created 13-Jul-2019
!
	INTEGER ID		!Number of ionization stage
	INTEGER NLEV		!Numer of levls in HN
        INTEGER N_DI		!Number of levels in target ion
	INTEGER ION_LEV		!Super level target in ION
	INTEGER ND		!Number of depth points
        INTEGER NION		!Numer of Eqns. in ionization matrix.
	INTEGER EQGS            !Equation number of GS in BA_T matrix.
!
! NB --- NION is the total number of ionic species i.e. for
! HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and CV).
!
	INTEGER DST,DEND
!
	REAL(KIND=LDP) WSE(NLEV,ND),dWSEdT(NLEV,ND)
	REAL(KIND=LDP) WCR(NLEV,ND),dWCRdT(NLEV,ND)
!
! Populations of species undergoing photoionization.
!
	REAL(KIND=LDP) HN(NLEV,ND),HNST(NLEV,ND),dlnHNST_dlnT(NLEV,ND)
!
! Ion populations.
!
	REAL(KIND=LDP) DI(N_DI,ND)
	REAL(KIND=LDP) LOG_DIST(N_DI,ND)
	REAL(KIND=LDP) dlnDIST_dlnT(N_DI,ND)
!
	REAL(KIND=LDP) ED(ND),T(ND)
	REAL(KIND=LDP) JREC(ND)
	REAL(KIND=LDP) dJRECdT(ND)
	REAL(KIND=LDP) JPHOT(ND)
	REAL(KIND=LDP) JREC_CR(ND)
	REAL(KIND=LDP) dJREC_CRdT(ND)
	REAL(KIND=LDP) JPHOT_CR(ND)
	LOGICAL FIXED_T
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables
!
	INTEGER J,K
	INTEGER JF
	INTEGER NT
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) B_RAT
	REAL(KIND=LDP) LOG_B_RAT
!
	REAL(KIND=LDP) PLANCKS_CONSTANT,PC
	EXTERNAL PLANCKS_CONSTANT
!
! REV_HNST referes to the LTE population  of the level defined with respect
! to the actual destination (target) level.
!
	REAL(KIND=LDP) REV_HNST
	REAL(KIND=LDP) WSE_BY_RJ,DI_FAC,ED_FAC,T_FAC
!
	PC=1.0E+15_LDP*PLANCKS_CONSTANT()
!
!
!
!$OMP PARALLEL DO PRIVATE(K,J,JF,LOG_B_RAT,B_RAT,WSE_BY_RJ,REV_HNST,T1,T2,T3,DI_FAC,ED_FAC,T_FAC)
!
	DO K=DST,DEND			!Which depth point.
	  IF(ION_LEV .NE. 1)THEN
	    LOG_B_RAT=LOG(DI(ION_LEV,K)/DI(1,K))+LOG_DIST(1,K)-LOG_DIST(ION_LEV,K)
	    B_RAT=0.0_LDP
	    IF(LOG_B_RAT .LT. 780.0_LDP)B_RAT=EXP(LOG_B_RAT)
	  ELSE
	    B_RAT=1.0_LDP
	    LOG_B_RAT=0.0_LDP
	  END IF
!
	  DO J=1,NLEV			!Which equation (for S.E. only)
	    IF(WSE(J,K) .NE. 0.0_LDP)THEN
	      WSE_BY_RJ=PC*(WSE(J,K)*JPHOT_CR(K)+WCR(J,K)*JPHOT(K))
	      JF=EQGS+J-1
	      BA_T_PAR_EHB(JF,K)=BA_T_PAR_EHB(JF,K)+WSE_BY_RJ
!
	      REV_HNST=HNST(J,K)*B_RAT
	      T3=PC*REV_HNST*(WCR(J,K)*JREC(K)+WSE(J,K)*JREC_CR(K))
	      DI_FAC=T3/DI(ION_LEV,K)
	      ED_FAC=T3/ED(K)
	      JF=EQGS+NLEV+(ION_LEV-1)
	      BA_T_PAR_EHB(JF,K)=BA_T_PAR_EHB(JF,K)  - DI_FAC
	      BA_T_PAR_EHB(NT-1,K) =BA_T_PAR_EHB(NT-1,K) - ED_FAC
!
	      T1=dWSEdT(J,K)*(REV_HNST*JREC_CR(K)-HN(J,K)*JPHOT_CR(K)) + REV_HNST*WSE(J,K)*dJREC_CRdT(K)
	      T2=dWCRdT(J,K)*(REV_HNST*JREC(K)-HN(J,K)*JPHOT(K)) + REV_HNST*WCR(J,K)*dJRECdT(K)
	      T_FAC=T3*( dlnHNST_dlnT(J,K) +(dlnDIST_dlnT(1,K)-dlnDIST_dlnt(ION_LEV,K)) )/T(K)
	      BA_T_PAR_EHB(NT,K)   =BA_T_PAR_EHB(NT,K)  - (T_FAC+PC*(T1+T2))
!
	    END IF		!WSE(J,K) .NE. 0
	  END DO
	END DO
!$OMP END PARALLEL DO
!
	RETURN
	END
