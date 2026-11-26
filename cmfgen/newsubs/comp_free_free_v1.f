!
! Subroutine to compute the contribution to the opacity AND emissivity
! by FREE-FREE and BOUND-FREE processes for a general ion. The
! contribution is added directly to the opacity CHI and emissivity ETA.
!
	SUBROUTINE COMP_FREE_FREE_V1(CHI_FF,ETA_FF,VCHI_FF,VETA_FF,CONT_FREQ,FREQ,DO_VAR,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 04-May-2022 : CHange HMI to H0
! Created 11-Jul-2019 : Based on COMP_OPAC and GENOPAETA_V10
!
	INTEGER ND
	INTEGER NT
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) CONT_FREQ
	REAL(KIND=LDP) CHI_FF(ND)			!Opacity
	REAL(KIND=LDP) ETA_FF(ND)			!Emissivity
	REAL(KIND=LDP) VCHI_FF(NT,ND)			!Opacity
	REAL(KIND=LDP) VETA_FF(NT,ND)			!Emissivity
	LOGICAL DO_VAR
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Vectors to save computational effort.
!
	REAL(KIND=LDP) EMHNUKT(ND)		!EXP(-hv/kT)
	REAL(KIND=LDP) GFF_VAL(ND)		!g(ff) as a function of depth
	REAL(KIND=LDP) POP_SUM(ND)		!Factor to convert HNST for ION_LEV
!
! Local constants.
!
	INTEGER ID
	INTEGER EQION
	INTEGER I,K,L
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) NU
	REAL(KIND=LDP) ALPHA,TCHI1,TETA1
	REAL(KIND=LDP) EMIS
	REAL(KIND=LDP) HNUONKT
!
	NU=FREQ
        T1=-HDKT*NU
        DO I=1,ND
          EMHNUKT(I)=EXP(T1/T(I))
          CHI_FF(I)=0.0_LDP
          ETA_FF(I)=0.0_LDP
        END DO
	VCHI_FF=0.0_LDP
	VETA_FF=0.0_LDP
!
! Add in free-free contribution. Because SN can be dominated by elements other
! than H and He, we now sum over all levels. To make sure that we only do this
! once, we only include the FREE-FREE contribution for the ion when PHOT_ID is one.
!
	DO ID=1,NUM_IONS
	  IF(.NOT. ATM(ID)%XzV_PRES)THEN
	  ELSE IF(ATM(ID)%ZXzV .EQ. 0.0_LDP)THEN
	    POP_SUM=ATM(ID)%XzV_F(1,1:ND)
	    I=7                                 !Used to read in data on first entry (will not be used here.)
	    CALL DO_H0_FF(ETA_FF,CHI_FF,POP_SUM,ED,T,EMHNUKT,NU,I,ND)
	  ELSE
!
! Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
!
	    CALL GFF_VEC(GFF_VAL,NU,T,ATM(ID)%ZXzV,ND)
	    CALL FF_RES_GAUNT(GFF_VAL,NU,T,ID,ATM(ID)%GIONXzV_F,ATM(ID)%ZXzV,ND)
!
! We use POP_SUM as a temporary vector containing the sum of all level populations in
! the ion at each depth.
!
	    POP_SUM(1:ND)=SUM(ATM(ID+1)%XzV,1)
	    TCHI1=CHIFF*ATM(ID)%ZXzV*ATM(ID)%ZXzV/(FREQ*FREQ*FREQ)
	    TETA1=CHIFF*ATM(ID)%ZXzV*ATM(ID)%ZXzV*TWOHCSQ
	    DO I=1,ND
	      ALPHA=ED(I)*POP_SUM(I)*GFF_VAL(I)/SQRT(T(I))
	      CHI_FF(I)=CHI_FF(I)+TCHI1*ALPHA*(1.0_LDP-EMHNUKT(I))
	      ETA_FF(I)=ETA_FF(I)+TETA1*ALPHA*EMHNUKT(I)
	    END DO
!
	    IF(DO_VAR)THEN
	      EQION=ATM(ID+1)%EQXzV
	      DO I=1,ND
	        ALPHA=POP_SUM(I)*GFF_VAL(I)/SQRT(T(I))
	        VCHI_FF(NT-1,I)=VCHI_FF(NT-1,I)+TCHI1*ALPHA*(1.0_LDP-EMHNUKT(I))
	        VETA_FF(NT-1,I)=VETA_FF(NT-1,I)+TETA1*ALPHA*EMHNUKT(I)
!
	        ALPHA=TCHI1*ED(I)*GFF_VAL(I)*(1.0_LDP-EMHNUKT(I))/SQRT(T(I))
	        EMIS=TETA1*ED(I)*GFF_VAL(I)*EMHNUKT(I)/SQRT(T(I))
	        DO L=1,ATM(ID+1)%NXzV
	           K=EQION+L-1
	           VCHI_FF(K,I)=VCHI_FF(K,I)+ALPHA
	           VETA_FF(K,I)=VETA_FF(K,I)+EMIS
	        END DO
!
	        HNUONKT=HDKT*NU/T(I)
	        ALPHA=ED(I)*POP_SUM(I)*GFF_VAL(I)/SQRT(T(I))/T(I)
		VCHI_FF(NT,I)=VCHI_FF(NT,I)-TCHI1*ALPHA*(0.5_LDP+(HNUONKT-0.5_LDP)*EMHNUKT(I))
		VETA_FF(NT,I)=VETA_FF(NT,I)+TETA1*ALPHA*(HNUONKT-0.5_LDP)*EMHNUKT(I)
	      END DO
	    END IF
	  END IF
	END DO
!
	RETURN
	END
