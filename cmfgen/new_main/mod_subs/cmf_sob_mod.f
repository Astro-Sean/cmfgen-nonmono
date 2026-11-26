!
! Module for making declarations required when using the Sobolev approximation, or CMF
! section to compute JBAR for individual lines. This data is NOT required for the blanketing
! section. Storage is allocated by a call to SET_CMF_SOB_MOD.
!
	MODULE CMF_SOB_MOD
	USE SET_KIND_MODULE
!
! Arrays for computation of JBAR using moment equations.
!
        REAL(KIND=LDP), ALLOCATABLE :: JNU(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: HNU(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: F_LINE(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: G_LINE(:,:)               !ND,NLF+1
        REAL(KIND=LDP), ALLOCATABLE :: HBC_LINE(:,:)             !3,NLF+1
        REAL(KIND=LDP), ALLOCATABLE :: NBC_LINE(:,:)             !3,NLF+1
        REAL(KIND=LDP), ALLOCATABLE :: IN_HBC_LINE(:)            !NLF+1
!
! Variables and Vectors for EW's and LINE blanketing.
!
	REAL(KIND=LDP), ALLOCATABLE :: JBLANK(:)                !ND
	REAL(KIND=LDP), ALLOCATABLE :: HBLANK(:)                !ND
	REAL(KIND=LDP), ALLOCATABLE :: JEW(:)                   !ND : Used to compute line Ew's.
	REAL(KIND=LDP), ALLOCATABLE :: JBAR(:)                  !ND : Mean line intensity.
	REAL(KIND=LDP), ALLOCATABLE :: ZNET(:)                  !ND : Net radiative rate
!
! Line profile arrays and variables.
!
	REAL(KIND=LDP), ALLOCATABLE :: PF(:)                    !NLF : Prof. freq. for line computations
	REAL(KIND=LDP), ALLOCATABLE :: PROF(:)                  !NLF : Line profile
	REAL(KIND=LDP), ALLOCATABLE :: LFQW(:)                  !NLF : Quad. weights assoc. with line prof
	REAL(KIND=LDP), ALLOCATABLE :: ERF(:)
!
! Used only in LINEGEN.INC
!
! Variation arrays
! Variable,depth of variable,depth of J. If NUM_BNDS .ne. ND the
! the variable depth is given by [ VJ(I1,I2,I3) ] I3+I2-NDIAG.
!
	REAL(KIND=LDP), ALLOCATABLE :: VZNET(:,:,:)             !NT,NUM_BNDS,ND
	REAL(KIND=LDP), ALLOCATABLE :: FQAF(:,:,:)              !ND,ND,NM_KI
	REAL(KIND=LDP), ALLOCATABLE :: FQAFD(:)                 !ND
!
	END MODULE CMF_SOB_MOD
!
! Allocate storage for Sobolev and CMF sections.
!
	SUBROUTINE SET_CMF_SOB_MOD(ND,NUM_BNDS,NT,NM_KI,NLF,LUER)
	USE SET_KIND_MODULE
	USE CMF_SOB_MOD
	IMPLICIT NONE
	INTEGER ND
	INTEGER NUM_BNDS
	INTEGER NT
	INTEGER NM_KI
	INTEGER NLF
	INTEGER LUER
	INTEGER IOS
!
	IOS=0
        IF(IOS .EQ.0)ALLOCATE(JNU(ND,NLF+1))
        IF(IOS .EQ.0)ALLOCATE(HNU(ND,NLF+1))
        IF(IOS .EQ.0)ALLOCATE(F_LINE(ND,NLF+1))
        IF(IOS .EQ.0)ALLOCATE(G_LINE(ND,NLF+1))
        IF(IOS .EQ.0)ALLOCATE(HBC_LINE(3,NLF+1))
        IF(IOS .EQ.0)ALLOCATE(NBC_LINE(3,NLF+1))
        IF(IOS .EQ.0)ALLOCATE(IN_HBC_LINE(NLF+1))
!
! Variables and Vectors for EW's and LINE blanketing.
!
	IF(IOS .EQ. 0)ALLOCATE( JBLANK(ND), STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( HBLANK(ND), STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( JEW(ND), STAT=IOS)    			!Used to compute line Ew's.
	IF(IOS .EQ. 0)ALLOCATE( JBAR(ND), STAT=IOS)			!Mean line intensity.
	IF(IOS .EQ. 0)ALLOCATE( ZNET(ND), STAT=IOS)			!Net radiative rate
!
! Line profile arrays and variables.
!
	IF(IOS .EQ. 0)ALLOCATE( PF(NLF), STAT=IOS)                  !Prof. freq. for line computations
	IF(IOS .EQ. 0)ALLOCATE( PROF(NLF), STAT=IOS)                !Line profile
	IF(IOS .EQ. 0)ALLOCATE( LFQW(NLF), STAT=IOS)                !Quad. weights assoc. with line prof
	IF(IOS .EQ. 0)ALLOCATE( ERF(NLF), STAT=IOS)
!
! Variation arrays
! Variable,depth of variable,depth of J. If NUM_BNDS .ne. ND the
! the variable depth is given by [ VJ(I1,I2,I3), STAT=IOS) ] I3+I2-NDIAG.
!
	IF(IOS .EQ. 0)ALLOCATE( VZNET(NT,NUM_BNDS,ND), STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( FQAF(ND,ND,NM_KI), STAT=IOS)	!Used only in LINEGEN.INC
	IF(IOS .EQ. 0)ALLOCATE( FQAFD(ND), STAT=IOS)
!
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in SET_CMF_SOB_MOD'
	  WRITE(LUER,*)'Unable to allcoate required memory'
	  WRITE(LUER,*)'IOS=',IOS
	END IF
!
	RETURN
	END
