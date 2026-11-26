!
! Routine to compute the collisional recobination and ionization rates, and
! the collisional cooling rate for ions with super levels.
!
	SUBROUTINE COLCOOL_SL_V5(CPR,CRR,COOL,CNM,DCNM,
	1             HN_S,HNST_S,dlnHNST_S_dlnT,N_S,
	1             HN_F,HNST_F_ON_S,A_F,W_F,EDGE_F,G_F,LEVNAME_F,
	1             F_TO_S_MAP,N_F,
	1             ZION,ID,COL_FILE,OMEGA_COL,ED,T,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 20-Sec-2023 : Updated to use standard value for Planck's constant (LONG -- 15-Oct-2023).
! Altered 04-Oct-2016 : SUBCOL_MULTI updated to V6.
! Altered 05-Apr-2011 : Updated from V4
!                          Most changes 02-Dec-10.
!                          HNST replaced in call by HNST_F_ON_S.
!	                   Call to SUBCOL_MULTI_V5 updated (from V4)
! Altered 24-Dec-1996 : SUB_PHOT replaces PHOT_FUN (superficial).
! Altered 24-May-1995 : N_F_MAX removed using F90
! Altered 01-Feb-1996 : Routine now calls SUBCOL_MULTI_V3 (not _V2). The _V3
!                          routine was introdued when interpolationg between
!                          levels in a given super level.
!                          Changed to _V3 for consistency with SUBCOL.
! Created 07-JUn-1995 :Based on COLGENCOOL
!
	REAL(KIND=LDP) PLANCKS_CONSTANT
	EXTERNAL OMEGA_COL, PLANCKS_CONSTANT
	EXTERNAL ERROR_LU
	INTEGER ERROR_LU	
!
	INTEGER ID
	INTEGER N_S,N_F,ND
	REAL(KIND=LDP) COOL(ND),CRR(ND),CPR(ND)
	REAL(KIND=LDP) CNM(N_S,N_S),DCNM(N_S,N_S)
	REAL(KIND=LDP) ZION,T(ND),ED(ND)
!
	REAL(KIND=LDP) HN_S(N_S,ND),HNST_S(N_S,ND)
	REAL(KIND=LDP) dlnHNST_S_dlnT(N_S,ND)
!
	REAL(KIND=LDP) HN_F(N_F,ND)
	REAL(KIND=LDP) HNST_F_ON_S(N_F,ND)	
	REAL(KIND=LDP) W_F(N_F,ND)
	REAL(KIND=LDP) A_F(N_F,N_F),EDGE_F(N_F),G_F(N_F)
	INTEGER F_TO_S_MAP(N_F)
	CHARACTER*(*) COL_FILE,LEVNAME_F(N_F)
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables.
!
	REAL(KIND=LDP) H,TMP_ED
	INTEGER I,J
	INTEGER, PARAMETER :: IONE=1
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	REAL(KIND=LDP) OMEGA_F(N_F,N_F)
	REAL(KIND=LDP) dln_OMEGA_F_dlnT(N_F,N_F)
!
	H=PLANCKS_CONSTANT()*1.0E+15_LDP             !ergs/s (*1.0E+15 due to *nu)
	TMP_ED=1.0_LDP
!
	DO I=1,ND			!Which depth
!
! Compute collisional cross-sections (and their T derivatives)
! The last line is COMPUTE_BA, FIXED_T, LST_ITERATION.
! With the adopted settings we do not compute dln_OMEGA_F_dlNT.
!
	  COOL(I)=0.0_LDP
	  CALL SUBCOL_MULTI_V6(OMEGA_F,dln_OMEGA_F_dlNT,
	1          CNM,DCNM,
	1          HN_S(1,I),HNST_S(1,I),dlnHNST_S_dlnT(1,I),N_S,
	1          HN_F(1,I),HNST_F_ON_S(1,I),W_F(1,I),EDGE_F,
	1          A_F,G_F,LEVNAME_F,N_F,
	1          ZION,ID,COL_FILE,OMEGA_COL,
	1          F_TO_S_MAP,COOL(I),T(I),TMP_ED,IONE,
	1          L_FALSE,L_TRUE,L_TRUE)
!
	  CPR(I)=0.0_LDP
	  CRR(I)=0.0_LDP
	  COOL(I)=COOL(I)*ED(I)*H
!
	  DO J=1,N_S				!Level
	    CPR(I)=CPR(I)+HN_S(J,I)*(ED(I)*CNM(J,J))
	    CRR(I)=CRR(I)+HNST_S(J,I)*(ED(I)*CNM(J,J))
	  END DO
!
	END DO
!
	RETURN
	END
