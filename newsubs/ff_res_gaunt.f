!
! Subroutine to incorporate a free-free line resonance into the
! normal free-free opacity. This is done by modifying the free-free
! gaunt factor for the species under consideration.
!
	SUBROUTINE FF_RES_GAUNT(GFF,FREQ,T,ID,GION,ZION,ND)
	USE SET_KIND_MODULE
	USE PHOT_DATA_MOD
	IMPLICIT NONE
!
	INTEGER ID		!Species identifier
	INTEGER ND		!Number of depth points
	REAL(KIND=LDP) FREQ		!Frequency in units of 10^15 Hz
	REAL(KIND=LDP) GION		!Ion statistical weight
	REAL(KIND=LDP) ZION		!Charge on ion.
	REAL(KIND=LDP) GFF(ND)		!Free-free Gaunt factor
	REAL(KIND=LDP) T(ND)		!Temperature in unitsof 10^4 K
!
	REAL(KIND=LDP) VOIGT
	EXTERNAL VOIGT
!
! Common block with opacity/emissivity constants.
!
        REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
!
! Local variables
!
	REAL(KIND=LDP) DOP_NU
	REAL(KIND=LDP) A_VOIGT
	REAL(KIND=LDP) V_VOIGT
	REAL(KIND=LDP) PHI_VOIGT
	REAL(KIND=LDP) CONST
	REAL(KIND=LDP) T1
!
	INTEGER I,K
!
! Check if any free-free resonances for this species.
!
	IF(PD(ID)%NUM_FF .EQ. 0)RETURN
!
! Compute resonance contribution. At present a fixed Doppler width is
! assumed. The Doppler width should be chosen to avoid undersampling of
! the intrinsic line profile with the adopted frequency grid. The intrinsic
! width of the line, set by the autoionization probabilities of the lower and
! upper levels, is taken into account.
!
	DO I=1,PD(ID)%NUM_FF
	  IF(FREQ .LE. PD(ID)%FF_NU_MAX(I) .AND. FREQ .GE. PD(ID)%FF_NU_MIN(I))THEN
            DOP_NU=PD(ID)%NU_ZERO(I)*PD(ID)%VSM_KMS/2.998E+05_LDP
            A_VOIGT=PD(ID)%GAMMA(I)/DOP_NU
            V_VOIGT=(FREQ-PD(ID)%FF_NU_ZERO(I))/DOP_NU
	    CONST=1.489E+15_LDP*PD(ID)%FF_GF(I)*(PD(ID)%FF_NU_ZERO(I)**3)/GION/ZION/ZION
            PHI_VOIGT=1.0E-15_LDP*VOIGT(A_VOIGT,V_VOIGT)/DOP_NU
            T1=HDKT*PD(ID)%FF_NU_EXCITE(I)
!	    WRITE(2,*)FREQ,PD(ID)%FF_NU_ZERO(I)
!	    WRITE(2,*)DOP_NU,A_VOIGT,V_VOIGT,CONST,PHI_VOIGT,T1
	    DO K=1,ND
	      GFF(K)=GFF(K)+CONST*PHI_VOIGT*EXP(-T1/T(K))/T(K)
	    END DO
	  END IF
	END DO
!
	RETURN
	END
