!
! Module for vectors required when computing the continuum radiation field,
! with blanketing. Some are simply scratch vectors. Others
! are required in CMFGEN_SUB and in the variation routines.
! Vectors required for dJ are in VAR_RAD_MOD.
!
	MODULE RADIATION_MOD
	USE SET_KIND_MODULE
!
	REAL(KIND=LDP), ALLOCATABLE :: RJ(:)      	!ND - Mean intensity  (ND)
	REAL(KIND=LDP), ALLOCATABLE :: RJ_ES(:)   	!ND - Convolution of RJ with e.s. R(v'v')
!
	REAL(KIND=LDP), ALLOCATABLE :: J_INT(:)   	!ND - Frequency integrated J
	REAL(KIND=LDP), ALLOCATABLE :: H_INT(:)   	!ND - Frequency integrated K
	REAL(KIND=LDP), ALLOCATABLE :: K_INT(:)   	!ND - Frequency integrated K
	REAL(KIND=LDP), ALLOCATABLE :: H_MOM(:)   	!ND - Frequency dependent H moment
	REAL(KIND=LDP), ALLOCATABLE :: K_MOM(:)   	!ND - Frequency dependent K moment
!
! Arrays for calculating mean opacities.
!
	REAL(KIND=LDP), ALLOCATABLE :: INT_dBdT(:)     	!ND - Int. of dB/dT dv (to calculate ROSSMEAN)
!
	REAL(KIND=LDP), ALLOCATABLE :: RLUMST(:)      	!ND - Luminosity as a function of depth
	REAL(KIND=LDP), ALLOCATABLE :: MECH_LUM(:)     	!ND - Mechanical luminosity
	REAL(KIND=LDP), ALLOCATABLE :: SOB(:)      	!ND - Used in computing continuum flux
	REAL(KIND=LDP), ALLOCATABLE :: LLUMST(:)      	!ND - Line luminosity.
	REAL(KIND=LDP), ALLOCATABLE :: DIELUM(:)      	!ND - Dielectronic line emission luminosity.
	REAL(KIND=LDP), ALLOCATABLE :: DEP_RAD_EQ(:)    !ND - Integrated departure from radiative equilibrium
	REAL(KIND=LDP), ALLOCATABLE :: DJDt_FLUX(:)     !ND - DJDT correction to integrated flux.
	REAL(KIND=LDP), ALLOCATABLE :: DJDt_TERM(:)     !ND -
!
! Vector giving the MINIMUM Doppler width at each depth.
!
	REAL(KIND=LDP), ALLOCATABLE :: VDOP_VEC(:)      !ND
!
! Transfer equation vectors
!
	REAL(KIND=LDP), ALLOCATABLE :: Z(:)        	!NDMAX - Z displacement along a given array
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: TB(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: TC(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: XM(:)        	!NDMAX - R.H.S. (SOURCE VECTOR)
	REAL(KIND=LDP), ALLOCATABLE :: DTAU(:)        	!NDMAX - Optical depth (used in error calcs)
	REAL(KIND=LDP), ALLOCATABLE :: dCHIdR(:)        !NDMAX - Derivative of opacity.
!
! Continuum matrices
!
	REAL(KIND=LDP), ALLOCATABLE :: WM(:,:)        	!ND,ND - Coef. matrix of J & %J vector
	REAL(KIND=LDP), ALLOCATABLE :: FB(:,:)        	!ND,ND - Coef. of J & %J vects in angular equ.
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	REAL(KIND=LDP), ALLOCATABLE :: FEDD(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: GEDD(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: QEDD(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: N_ON_J(:)        !NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: RSQHNU(:)        !NDMAX -
!
! For MOM_J_REL_V2 --- i.e., the inclusion of relatvistic, but not time
! dependence, for the computation of J.
!
	REAL(KIND=LDP), ALLOCATABLE :: N_ON_J_NODE(:)		!
	REAL(KIND=LDP), ALLOCATABLE :: H_ON_J(:)		!
	REAL(KIND=LDP), ALLOCATABLE :: KMID_ON_J(:)		!
	REAL(KIND=LDP), ALLOCATABLE :: dlnJdlNR(:)		!
!
! Boundary conditions.
!
	REAL(KIND=LDP) HBC_CMF(3)
	REAL(KIND=LDP) NBC_CMF(3)
	REAL(KIND=LDP) INBC
	REAL(KIND=LDP) HBC_J
	REAL(KIND=LDP) HBC_S			!Bound. Cond. for JFEAU
	REAL(KIND=LDP) HBC_PREV(3)
	REAL(KIND=LDP) NBC_PREV(3)
	REAL(KIND=LDP) INBC_PREV
!
	REAL(KIND=LDP), ALLOCATABLE :: FEDD_PREV(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: GEDD_PREV(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: N_ON_J_PREV(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: JNU_PREV(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: RSQHNU_PREV(:)        	!NDMAX -
	REAL(KIND=LDP), ALLOCATABLE :: FOLD(:)        		!NDMAX
!
! 
!
! Use for interpolating opacities tec onto the foine grid.
!
	INTEGER, ALLOCATABLE :: INDX(:)                 !NDMAX
	INTEGER, ALLOCATABLE :: POS_IN_NEW_GRID(:)      !ND
	REAL(KIND=LDP), ALLOCATABLE :: COEF(:,:)        	!0:3,NDMAX
!
! Variables and arrays required on the fine grid.
!
	REAL(KIND=LDP), ALLOCATABLE :: REXT(:)        		!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: VEXT(:)        		!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: LANG_COORDEXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: TEXT(:)        		!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: SIGMAEXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: VDOP_VEC_EXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: CHIEXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: ESECEXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: ETAEXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: ZETAEXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: THETAEXT(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: RJEXT(:)        		!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: RJEXT_ES(:)        	!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: FEXT(:)        		!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: QEXT(:)        		!NDMAX
	REAL(KIND=LDP), ALLOCATABLE :: SOURCEEXT(:)        	!NDMAX
!
! Diffusion approximation variables
!
	REAL(KIND=LDP) DTDR
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) DDBBDT
!
	REAL(KIND=LDP) dLOG_NU
	LOGICAL CONT_VEL
!
	END MODULE RADIATION_MOD
!
	SUBROUTINE SET_RADIATION_MOD(ND,NDMAX,NPMAX)
	USE SET_KIND_MODULE
	USE RADIATION_MOD
	IMPLICIT NONE
	INTEGER ND,NDMAX,NPMAX
!
	INTEGER IOS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	IOS=0
	IF(IOS .EQ. 0)ALLOCATE ( RJ(ND),STAT=IOS )		!Mean intensity
	IF(IOS .EQ. 0)ALLOCATE ( RJ_ES(ND),STAT=IOS )		!Convolution of RJ with e.s. R(v'v')
	IF(IOS .EQ. 0)ALLOCATE ( J_INT(ND),STAT=IOS )		!Frequency integrated J
	IF(IOS .EQ. 0)ALLOCATE ( H_INT(ND),STAT=IOS )		!Frequency integrated H
	IF(IOS .EQ. 0)ALLOCATE ( K_INT(ND),STAT=IOS )		!Frequency integrated K
	IF(IOS .EQ. 0)ALLOCATE ( H_MOM(ND),STAT=IOS )		!Frequency dependent H moment
	IF(IOS .EQ. 0)ALLOCATE ( K_MOM(ND),STAT=IOS )		!Frequency dependent K moment
!
! Arrays for calculating mean opacities.
!
	IF(IOS .EQ. 0)ALLOCATE ( INT_dBdT(ND),STAT=IOS )  	!Integral of dB/dT over nu (to calculate ROSSMEAN)
!
	IF(IOS .EQ. 0)ALLOCATE ( RLUMST(ND),STAT=IOS )		!Luminosity as a function of depth
	IF(IOS .EQ. 0)ALLOCATE ( MECH_LUM(ND),STAT=IOS )	!Mechanical luminosity
	IF(IOS .EQ. 0)ALLOCATE ( SOB(ND),STAT=IOS )   		!Used in computing continuum flux
	IF(IOS .EQ. 0)ALLOCATE ( LLUMST(ND),STAT=IOS )    	!Line luminosity.
	IF(IOS .EQ. 0)ALLOCATE ( DIELUM(ND),STAT=IOS )    	!Dielectronic line emission luminosity.
	IF(IOS .EQ. 0)ALLOCATE ( DEP_RAD_EQ(ND),STAT=IOS )    	!Depature from radiative equilibrium.
	IF(IOS .EQ. 0)ALLOCATE ( DJDt_TERM(ND),STAT=IOS )	!
	IF(IOS .EQ. 0)ALLOCATE ( DJDt_FLUX(ND),STAT=IOS )    	!DJDt correction to integrated flux.
	IF(IOS .EQ. 0)ALLOCATE ( VDOP_VEC(ND),STAT=IOS )
!
! Transfer equation vectors
!
	IF(IOS .EQ. 0)ALLOCATE ( Z(NDMAX),STAT=IOS )		!Z displacement along a given array
	IF(IOS .EQ. 0)ALLOCATE ( TA(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( TB(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( TC(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( XM(NDMAX),STAT=IOS )		!R.H.S. (SOURCE VECTOR)
	IF(IOS .EQ. 0)ALLOCATE ( DTAU(NDMAX),STAT=IOS )       	!Optical depth (used in error calcs)
	IF(IOS .EQ. 0)ALLOCATE ( dCHIdR(NDMAX),STAT=IOS ) 	!Derivative of opacity.
!
! Continuum matrices
!
	IF(IOS .EQ. 0)ALLOCATE ( WM(ND,ND),STAT=IOS )		!Coef. matrix of J & %J vector
	IF(IOS .EQ. 0)ALLOCATE ( FB(ND,ND),STAT=IOS )		!Coef. of J & %J vects in angular equ.
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	IF(IOS .EQ. 0)ALLOCATE ( FEDD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( GEDD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( QEDD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( N_ON_J(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RSQHNU(NDMAX),STAT=IOS )
!
	IF(IOS .EQ. 0)ALLOCATE ( N_ON_J_NODE(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( H_ON_J(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( KMID_ON_J(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dlnJdlnR(NDMAX),STAT=IOS )
!
	IF(IOS .EQ. 0)ALLOCATE ( FEDD_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( GEDD_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( N_ON_J_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( JNU_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RSQHNU_PREV(NDMAX),STAT=IOS )
!
! 
!
	IF(IOS .EQ. 0)ALLOCATE ( INDX(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( POS_IN_NEW_GRID(ND),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( COEF(0:3,NDMAX),STAT=IOS )
!
	IF(IOS .EQ. 0)ALLOCATE ( REXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( VEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( LANG_COORDEXT(NDMAX),STAT=IOS ); LANG_COORDEXT=0.0_LDP
	IF(IOS .EQ. 0)ALLOCATE ( TEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( SIGMAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( VDOP_VEC_EXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( CHIEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( ESECEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( ETAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( ZETAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( THETAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RJEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RJEXT_ES(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( FOLD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( FEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( QEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( SOURCEEXT(NDMAX),STAT=IOS )
!
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error allocating memory in SET_RADIATION_MOD'
	  WRITE(LUER,*)'STAT=',IOS
	  STOP
	END IF
!
! The following initializations are required when we do a model with 0 iterations.
!
	J_INT=0.0_LDP; H_INT=0.0_LDP; K_INT=0.0_LDP
!
	RETURN
	END
