!
! Module to contain the data vectors with d(dCHIdR)/dCHI. Replaces common block
! TRAPDERIVATIVES.
!
! A(I)= d[dCHIdR(I)]/dCHI(I-1)
! B(I)= d[dCHIdR(I)]/dCHI(I)
! C(I)= d[dCHIdR(I)]/dCHI(I+1)
!
! Created: 02-Mar-1998. Allows arbitrarily large dimensions.
!
	MODULE MOD_TRAP_DERIVATIVES
	USE SET_KIND_MODULE
	  INTEGER ND_TRAP
	  REAL(KIND=LDP), ALLOCATABLE :: A(:)
	  REAL(KIND=LDP), ALLOCATABLE :: B(:)
	  REAL(KIND=LDP), ALLOCATABLE :: C(:)
          DATA ND_TRAP/0/
	END MODULE MOD_TRAP_DERIVATIVES
