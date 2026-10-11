C Stub for SCR_READ_V2 -- required only to satisfy the linker for
C DO_NG_BAND_ACCEL_V4's unused caller path.  The unit tests call
C NG_MIT_BAND_OPT_V1 directly and never reach this routine.
	SUBROUTINE SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,
     1                       LST_NG_REC_CHANGE,ND,NT,LUSCR,LUER)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER ND,NT,IREC,NITSF,RITE_N_TIMES,LST_NG_REC_CHANGE
	INTEGER LUSCR,LUER
	REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND),POPS(NT,ND)
	WRITE(LUER,*)'SCR_READ_V2 stub: not supported in unit tests'
	STOP
	END
