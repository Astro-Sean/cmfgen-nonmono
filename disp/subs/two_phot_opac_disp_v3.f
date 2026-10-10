!
! Routine to increment the Emissivity and Opacity for 2-photon processes.
!
! This routine should be the same at TWO_PHOT_DISP_V2 except it works directly on the
! atom rather than pops, and ION_ID is included in the call.
!
	SUBROUTINE TWO_PHOT_OPAC_DISP_V3(ETA,CHI,POPS,T,FREQ,TWO_PHOTON_METHOD,ION_ID,ND,NT)
	USE SET_KIND_MODULE
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Altered 14-Jul-2015: Added two options to improve two photon absorption.
!                        USE_TWO  - Uses actual radiation field from EDDFACTOR file.
!                        USE_W    - Uses dilution factor.
!	
! Created 26-Jun-1998
!
	INTEGER NT,ND,ION_ID
	REAL(KIND=LDP) ETA(ND)		!Emisivity
	REAL(KIND=LDP) CHI(ND)		!Opacity [in (10^10 cm)^-1]
	REAL(KIND=LDP) POPS(NT,ND)	!Atomic populations
	REAL(KIND=LDP) T(ND)		!in 10^4 K
	REAL(KIND=LDP) FREQ		!Current frequency (10^15 Hz)
!
	CHARACTER(LEN=*) TWO_PHOTON_METHOD
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local vectors and constants.
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) PLANKS_CONSTANT	!cgs units
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) CONST
	REAL(KIND=LDP) ETA_CONST	!Used to evaluate ETA
	REAL(KIND=LDP) CHI_CONST	!Used to evaluate CHI
	REAL(KIND=LDP) FREQ_B		!Frequency of other photon
	REAL(KIND=LDP) T1,T2
!
! The 2-photon distribution functions, AY, are usually in wrtten in terms
! of the variable y=FREQ/MAX_FREQ, and which extends from 0 to 1.
! As the distribution, AY, is symmetric about y/2 we have defined the
! variables  U=Y(1-Y) and FU=4*U
!
	REAL(KIND=LDP) AY,Y,U,FU
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER J,L,ML
	INTEGER NL,NUP
!
	PLANKS_CONSTANT=6.626E-27_LDP			!cgs units
	PI=4.0_LDP*ATAN(1.0_LDP)
	LUER=ERROR_LU()
	IF(TWO_METHOD .NE. TWO_PHOTON_METHOD)THEN
	  TWO_METHOD=TWO_PHOTON_METHOD
	  WRITE(LUER,*)'Using ',TRIM(TWO_METHOD),' method for two-photon decay'
	END IF
!
! The factor 10^10 arises since R is units of 10^10cm, and we
! scale CHI (and ETA) so that R.CHI is dimensionless, and
! ETA/CHI has the dimensions of the Planck Function.
! We don't have to worry about the frequency units, as it is
! multilplied by a ratio of 2 frequencies.
!
	CONST=1.0E+10_LDP*PLANKS_CONSTANT/4.0_LDP/PI
!
	DO J=1,N_TWO
	  IF(ION_ID_TWO(J) .EQ. ION_ID .AND. TWO_PHOT_AVAILABLE(J) .AND. FREQ .LT. FREQ_TWO(J))THEN
!
	    FREQ_B=FREQ_TWO(J)-FREQ
	    NL=LOW_LEV_TWO(J)
	    NUP=UP_LEV_TWO(J)
!
	    IF(TYPE_TWO(J) .EQ. 1)THEN
	      Y=FREQ/FREQ_TWO(J)
	      U=Y*(1.0_LDP-Y)
	      FU=4.0_LDP*U
	      AY=24.56_LDP*COEF_TWO(J,1)*( U*(1.0_LDP-FU**0.8_LDP) +
	1                 0.88_LDP*(U**1.53_LDP)*(FU**0.8_LDP) )
	    ELSE
	      WRITE(LUER,*)'Error in TWO_PHOT_OPAC'
	      WRITE(LUER,*)'Unrecognized two photon transition'
	      STOP
	    END IF
!
! The exact expressions for the two photon emissivity and opacity are:
!
! Eta=  (hv/4pi) (nu/nu_o) A(v) N_u [1 + Y(nu_B)]
! Chi=  (c^2/2hv^3) (hv/4pi) (nu/nu_o) A(v) .
!	    [ N_l (g_u/g_l) Y(nu_B) - N_u {1+Y(nu_B)} ]
! where
!       Y(v)=J(v)/[2hv^3/c^2]	
!
	    IF(TWO_METHOD .EQ. 'USE_RAD')THEN
	      CALL GET_J_FOR_TWO_PHOT(PHOT_OC_TWO(1,J),FREQ_B,LST_FREQ_INDX_TWO(J),ND)
	      T1=TWOHCSQ*(FREQ_B**3)
	      PHOT_OC_TWO(:,J)=PHOT_OC_TWO(:,J)/T1
	    ELSE IF(TWO_METHOD .EQ. 'LTE')THEN
	      DO L=1,ND
	        T1=EXP(-HDKT*FREQ_B/T(L))
	        PHOT_OC_TWO(L,J)=T1/(1.0_LDP-T1)
	      END DO
	    ELSE IF(TWO_METHOD .EQ. 'NOSTIM')THEN
	      PHOT_OC_TWO(:,J)=0.0_LDP
	    ELSE IF(TWO_METHOD .EQ. 'OLD_DEFAULT')THEN
	    ELSE
	      WRITE(LUER,*)'Unrecognized option in TWO_PHOT_OPAC_V3'
	      WRITE(LUER,*)'Options is: ',TRIM(TWO_METHOD)
	      WRITE(LUER,*)'Available options are USE_RAD, LTE, NOSTIM, and OLD_DEFAULT'
	      STOP
	    END IF
!
	    ETA_CONST=CONST*FREQ/FREQ_TWO(J)
	    CHI_CONST=G_UP_TWO(J)*ETA_CONST/TWOHCSQ/FREQ**3
!
	    IF(TWO_METHOD .EQ. 'OLD_DEFAULT')THEN
!
! The expressions for CHI and ETA are approximate. Their ratio gives the
! Planck function at depth, and the expression for ETA is exact at low
! densities and for small dilution factors.
!
	      DO L=1,ND
	        ETA(L)=ETA(L) + ETA_CONST*AY*POPS(NUP,L)*FS_RAT_UP(L,J)
	        T1=EXP(-HDKT*FREQ_B/T(L))
	        CHI(L)=CHI(L) + CHI_CONST*AY*(
	1                 POPS(NL,L)*FS_RAT_LOW(L,J)*T1/G_LOW_TWO(J)-
	1                 POPS(NUP,L)*FS_RAT_UP(L,J)/G_UP_TWO(J) )
	      END DO
!
	    ELSE
!
! When PHOT_OC_TWO is correctly set, these expressions are exact.
!
	      DO L=1,ND
	        ETA(L)=ETA(L) + ETA_CONST*AY*POPS(NUP,L)*FS_RAT_UP(L,J)*(1.0_LDP+PHOT_OC_TWO(L,J))
	        CHI(L)=CHI(L) + CHI_CONST*AY*(
	1                 POPS(NL,L)*FS_RAT_LOW(L,J)*PHOT_OC_TWO(L,J)/G_LOW_TWO(J)-
	1                 POPS(NUP,L)*FS_RAT_UP(L,J)/G_UP_TWO(J)*(1.0_LDP+PHOT_OC_TWO(L,J)) )
	      END DO
!
	    END IF
	  END IF
	END DO
!
	RETURN
	END
