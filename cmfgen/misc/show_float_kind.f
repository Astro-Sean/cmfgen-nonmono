	PROGRAM SHOW_FLOAT_KIND
	USE SET_KIND_MODULE
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	REAL(KIND=LDP) T1
	WRITE(6,*)'The number of storage bits in the selected floating'//
	1    ' point format is:',STORAGE_SIZE(T1)
	WRITE(6,*)'The KIND of the selected floating'//
	1    ' point format is:',KIND(T1)
!
	STOP
	END
