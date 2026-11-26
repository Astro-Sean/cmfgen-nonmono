C
C Routine to compute the collisional recobination and ionization rates, and
C the collisional cooling rate for ions with super levels.
C
	SUBROUTINE COLCOOL_SL_V4(CPR,CRR,COOL,CNM,DCNM,
	1             HN_S,HNST_S,dlnHNST_S_dlnT,N_S,
	1             HN_F,HNST_F,A_F,W_F,EDGE_F,G_F,LEVNAME_F,
	1             F_TO_S_MAP,N_F,
	1             ZION,ID,COL_FILE,OMEGA_COL,ED,T,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 24-Dec-1996 : SUB_PHOT replaces PHOT_FUN (superficial).
C Altered 24-May-1995 : N_F_MAX removed using F90
C Altered 01-Feb-1996 : Routine now calls SUBCOL_MULTI_V3 (not _V2). The _V3
C                        routine was introdued when interpolationg between
C                        levels in a given super level.
C                        Changed to _V3 for consistency with SUBCOL.
C Created 07-JUn-1995 :Based on COLGENCOOL
C
	EXTERNAL OMEGA_COL
	EXTERNAL ERROR_LU
	INTEGER ERROR_LU	
C
	INTEGER ID
	INTEGER N_S,N_F,ND
	REAL(KIND=LDP) COOL(ND),CRR(ND),CPR(ND)
	REAL(KIND=LDP) CNM(N_S,N_S),DCNM(N_S,N_S)
	REAL(KIND=LDP) ZION,T(ND),ED(ND)
C
	REAL(KIND=LDP) HN_S(N_S,ND),HNST_S(N_S,ND)
	REAL(KIND=LDP) dlnHNST_S_dlnT(N_S,ND)
C
	REAL(KIND=LDP) HN_F(N_F,ND),HNST_F(N_F,ND)	
	REAL(KIND=LDP) W_F(N_F,ND)
	REAL(KIND=LDP) A_F(N_F,N_F),EDGE_F(N_F),G_F(N_F)
	INTEGER F_TO_S_MAP(N_F)
	CHARACTER*(*) COL_FILE,LEVNAME_F(N_F)
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables.
C
	REAL(KIND=LDP) H,TMP_ED
	INTEGER I,J,IONE
	PARAMETER (IONE=1)
C
	REAL(KIND=LDP) OMEGA_F(N_F,N_F)
	REAL(KIND=LDP) dln_OMEGA_F_dlnT(N_F,N_F)
C
	H=6.6261965E-12_LDP		!H*1.0E+15  (1.0E+15 due to times frequency)
	TMP_ED=1.0_LDP
C
	DO I=1,ND			!Which depth
C
C Compute collisional cross-sections (and their T derivatives)
C
	  COOL(I)=0.0_LDP
	  CALL SUBCOL_MULTI_V4(OMEGA_F,dln_OMEGA_F_dlNT,
	1          CNM,DCNM,
	1          HN_S(1,I),HNST_S(1,I),dlnHNST_S_dlnT(1,I),N_S,
	1          HN_F(1,I),HNST_F(1,I),W_F(1,I),EDGE_F,
	1          A_F,G_F,LEVNAME_F,N_F,
	1          ZION,ID,COL_FILE,OMEGA_COL,
	1          F_TO_S_MAP,COOL(I),T(I),TMP_ED,IONE)
C
	  CPR(I)=0.0_LDP
	  CRR(I)=0.0_LDP
	  COOL(I)=COOL(I)*ED(I)*H
C
	  DO J=1,N_S				!Level
	    CPR(I)=CPR(I)+HN_S(J,I)*ED(I)*CNM(J,J)
	    CRR(I)=CRR(I)+HNST_S(J,I)*ED(I)*CNM(J,J)
	  END DO
C
	END DO
C
	RETURN
	END
