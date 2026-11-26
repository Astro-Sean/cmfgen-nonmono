	LOGICAL FUNCTION JFNK_T_FIXED_MODE_V1(FIXED_T_MODE,RD_FIX_T,
	1                                      JFNK_FIXT)
	IMPLICIT NONE
	LOGICAL FIXED_T_MODE,RD_FIX_T,JFNK_FIXT
!
! RD_FIX_T is the model policy and must override temporary mode changes.
!
	JFNK_T_FIXED_MODE_V1=FIXED_T_MODE .OR. RD_FIX_T .OR. JFNK_FIXT
	RETURN
	END
