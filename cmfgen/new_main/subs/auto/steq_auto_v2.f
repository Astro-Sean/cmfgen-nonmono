!
! Subroutine to compute the value of the statistical equilibrium
! equations and the variation of the statistical equilibrium matrix for
! AUTOINIZATION terms.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
! NB - ZION is the charge on the ion - thus ZHYD=1.0D0
!
! Routine does not work for NUM_BNDS=ND.
!
	SUBROUTINE STEQ_AUTO_V2(ED,T,
	1       HN_S,N_S,DI_S,
	1       HN_F,HNST_F,FEDGE_F,G_F,LEVNAME_F,N_F,
	1       F_TO_S_MAPPING,ID,
	1       DIERECOM,DIECOOL,AUTO_FILE,
	1       NUM_BNDS,ND,COMPUTE_BA,DST,DEND)
	USE SET_KIND_MODULE
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 28-Aug-2015 : Altered to allow levels VERY close (but below isolated atom threshold)
!                          to autoionize.
! Altered 23-Sep-2005 : Bug fix for dT term in computation of BA.
! Altered 13-Sep-2002 : Bug fix & cooling term added.
!                       DIECOOL included in call hence changed to V2
!                       STEQ(EQION, ) was not being computed.
!
	INTEGER NT
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER DST,DEND
!
	INTEGER ID
	INTEGER N_S,N_F
!
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) ED(ND)
!
	REAL(KIND=LDP) HN_S(N_S,ND)
	REAL(KIND=LDP) DI_S(ND)
!
	REAL(KIND=LDP) HN_F(N_F,ND)
	REAL(KIND=LDP) HNST_F(N_F,ND)
	REAL(KIND=LDP) FEDGE_F(N_F)
	REAL(KIND=LDP) G_F(N_F)
	CHARACTER*(*) LEVNAME_F(N_F)
	CHARACTER*(*) AUTO_FILE
	INTEGER F_TO_S_MAPPING(N_F)
	LOGICAL COMPUTE_BA
!
	REAL(KIND=LDP) DIERECOM(ND)
	REAL(KIND=LDP) DIECOOL(ND)
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) PLNCKS_CONST
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL(KIND=LDP) AUTO(N_F)
!
	INTEGER EQION
	INTEGER EQ_NUM_CONV
	INTEGER VION,VED,VT
	INTEGER I,J,K,M,IBEG_AUTO
	LOGICAL AUTO_FILE_EXISTS
!
	PLNCKS_CONST=6.6261965E-12_LDP         !H*1.0E+15  (1.0E+15 due to times frequency)
!
! The autoionization probabilities are depth independent, and hence can be
! read in at the beginning.
!
        INQUIRE(FILE=AUTO_FILE,EXIST=AUTO_FILE_EXISTS)
        IF(AUTO_FILE_EXISTS)THEN
          CALL RD_AUTO_V1(AUTO,FEDGE_F,G_F,LEVNAME_F,N_F,AUTO_FILE)
	  IBEG_AUTO=0
	  DO I=1,N_F
	    IF(AUTO(I) .GT. 0.0_LDP)THEN
	      IBEG_AUTO=I
	      EXIT
	    END IF
	  END DO
	  IF(IBEG_AUTO .EQ. 0)RETURN
        ELSE IF(FEDGE_F(N_F) .LT. 0.0_LDP)THEN
          LUER=ERROR_LU()
          WRITE(LUER,*)' '
          WRITE(LUER,*)'Warning: possible error in STEQ_AUTO_V1'
          WRITE(LUER,*)'No autoionization probabilities available'
          WRITE(LUER,*)'Model has states above the ionization limit'
          WRITE(LUER,'(A,A)')' AUTO_FILE is ',TRIM(AUTO_FILE)
          WRITE(LUER,*)' '
          RETURN
	ELSE
	  RETURN
        END IF
!
        EQION=N_S+1			!Ion equation
        VION=N_S+1
        VED=SE(ID)%N_IV-1
        VT=SE(ID)%N_IV
	EQ_NUM_CONV=SE(ID)%NUMBER_BAL_EQ
	M=(NUM_BNDS/2)+1
!
! 
!
	DO I=DST,DEND			!Which depth
!
	  DO K=IBEG_AUTO,N_F
	    J=F_TO_S_MAPPING(K)
	    T1=AUTO(K)*(HNST_F(K,I)-HN_F(K,I))
	    SE(ID)%STEQ(J,I)=SE(ID)%STEQ(J,I)+T1
	    SE(ID)%STEQ(EQION,I)=SE(ID)%STEQ(EQION,I)-T1
	    DIERECOM(I)=DIERECOM(I)+T1
	    DIECOOL(I)=DIECOOL(I)-PLNCKS_CONST*FEDGE_F(K)*T1
	  END DO
!
	  IF(COMPUTE_BA)THEN
	    DO K=IBEG_AUTO,N_F			!Which level
	      J=F_TO_S_MAPPING(K)
!
	      SE(ID)%BA(J,J,M,I)       =SE(ID)%BA(J,J,M,I)        -AUTO(K)*HN_F(K,I)/HN_S(J,I)
	      SE(ID)%BA(EQION,J,M,I)   =SE(ID)%BA(EQION,J,M,I)    +AUTO(K)*HN_F(K,I)/HN_S(J,I)
!
	      SE(ID)%BA(J,VION,M,I)    =SE(ID)%BA(J,VION,M,I)     +AUTO(K)*HNST_F(K,I)/DI_S(I)
	      SE(ID)%BA(EQION,VION,M,I)=SE(ID)%BA(EQION,VION,M,I) -AUTO(K)*HNST_F(K,I)/DI_S(I)
!
	      SE(ID)%BA(J,VED,M,I)     =SE(ID)%BA(J,VED,M,I)      +AUTO(K)*HNST_F(K,I)/ED(I)
	      SE(ID)%BA(EQION,VED,M,I) =SE(ID)%BA(EQION,VED,M,I)  -AUTO(K)*HNST_F(K,I)/ED(I)
!
	      SE(ID)%BA(J,VT,M,I)      =SE(ID)%BA(J,VT,M,I)       -AUTO(K)*HNST_F(K,I)*
	1                                                               (1.5_LDP+HDKT*FEDGE_F(K)/T(I))/T(I)
	      SE(ID)%BA(EQION,VT,M,I)  =SE(ID)%BA(EQION,VT,M,I)   +AUTO(K)*HNST_F(K,I)*
	1                                                               (1.5_LDP+HDKT*FEDGE_F(K)/T(I))/T(I)
!
	     END DO          		!Over level variable
	  END IF
!
	END DO		!Loop over depth (I)
!
	IF(DST .EQ. 1 .AND. DEND .EQ. ND)THEN
	  OPEN(UNIT=10,FILE=TRIM(AUTO_FILE)//'RATE',STATUS='UNKNOWN')
	    WRITE(10,'(A)')'!'
	    WRITE(10,'(A)')'! Anti-autoionization ratse, Autoionization rates '
	    WRITE(10,'(A)')'!'
	    WRITE(10,'(I5,T20,A)')N_F-IBEG_AUTO+1,'!Number of transitions '
	    WRITE(10,'(I5,T20,A)')ND,'!Number of depths'
	    WRITE(10,'(A)')'!'
	    DO K=IBEG_AUTO,N_F
	      J=F_TO_S_MAPPING(K)
	      WRITE(10,'(A,3X,I5)')LEVNAME_F(K),K
	      WRITE(10,'(8ES16.8)')(AUTO(K)*HNST_F(K,I),AUTO(K)*HN_F(K,I),I=1,ND)
	      WRITE(10,'(A)')' '
	    END DO
	  CLOSE(UNIT=10)
	END IF
!	
	RETURN
	END
