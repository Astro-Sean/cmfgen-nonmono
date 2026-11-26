C
C Module to store data necssary to compute occupation probabilities
C etc associated with level dissolution.
C
C The vectors should be defined as follows:
C
C  B_LEV_DIS=( 8.589E+14*(POPION**0.333333D0)/ED)**1.5D0
C  A_LEV_DIS=0.0009*(ED**0.1667)/SQRT(T)
C  X_LEV_DIS=(1+A_LEV_DIS)**3.15
C
C NB: B_LEV_DIS is undefined for ED=0, but this should be of no
C       practical concern, if we specifically handle the case when
C       level dissoultion is switched off (corresponds to ED=0)
C
C Based on the include file LEV_DIS_BLK.INC. Has dynamic allocation
C and variable indicating whether level dissolution is switched on
C in now explicitly included in the data block.
C
	MODULE MOD_LEV_DIS_BLK
	REAL(KIND=LDP), ALLOCATABLE :: B_LEV_DIS(:)
	REAL(KIND=LDP), ALLOCATABLE :: A_LEV_DIS(:)
	REAL(KIND=LDP), ALLOCATABLE :: X_LEV_DIS(:)
	LOGICAL MOD_DO_LEV_DIS
C
	END MODULE MOD_LEV_DIS_BLK
C
C Subroutine to:
C          (1) Allocate desired memory if not done already.
C          (2) Compute level dissolution constants.
C
	SUBROUTINE COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DIS,ND)
	USE SET_KIND_MODULE
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
C
	INTEGER ND
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) POPION(ND)
	REAL(KIND=LDP) T(ND)
	LOGICAL DO_LEV_DIS
C
	IF(.NOT. ALLOCATED(B_LEV_DIS))THEN
	  ALLOCATE( B_LEV_DIS(ND) )
	  ALLOCATE( A_LEV_DIS(ND) )
	  ALLOCATE( X_LEV_DIS(ND) )
	END IF
C
	MOD_DO_LEV_DIS=DO_LEV_DIS
	IF(DO_LEV_DIS)THEN
	  B_LEV_DIS=( 8.589E+14_LDP*(POPION**0.333333_LDP)/ED)**1.5_LDP
	  A_LEV_DIS=0.0009_LDP*(ED**0.1667_LDP)/SQRT(T)
	  X_LEV_DIS=(1.0_LDP+A_LEV_DIS)**3.15_LDP
	ELSE
	  B_LEV_DIS=1.0E+30_LDP		!Large number
	  A_LEV_DIS=0.0_LDP
	  X_LEV_DIS=1.0_LDP
	END IF
C
	RETURN
	END
