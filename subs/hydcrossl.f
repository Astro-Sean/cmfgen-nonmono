C
C Function to return the hydrogenic cross section. Cross-section has
C been scaled by 10^10 (Thus at threshold they should be of
C the order of 10^{-8}. N, LST and LEND must be integer valued.
C All other variables must be double precision.
C
	FUNCTION HYDCROSSL(NU,NUION,Z,NEF,N,LST,LEND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 24-May-1996 ERROR_LU etc installed
C Created 28-Aug-1989
C
	REAL(KIND=LDP) HYDCROSSL
	REAL(KIND=LDP) NU,NUION,Z,NEF
	INTEGER N,LST,LEND
C
	INTEGER ERROR_LU,LUER
	LOGICAL EQUAL
	EXTERNAL ERROR_LU,EQUAL
C
	INTEGER L,NMAX
	PARAMETER (NMAX=100)
	REAL(KIND=LDP) G(0:NMAX,0:NMAX),SIGMA(0:NMAX)
	REAL(KIND=LDP) NSAV,USAV,GSUM,U
	SAVE NSAV,USAV,SIGMA
C
	IF(N .GT. 100 .OR. N .LT. 1)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid N in HYDCROSSL - LIMIT IS 100'
	  WRITE(LUER,*)'N=',N,' NUION=',NUION,' Z=',Z
	  STOP
	END IF
C
	IF(LST .GT. LEND .OR. LST .LT. 0 .OR. LEND .GT. N-1)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid LST and LEND in HYDCROSSL'
	  WRITE(LUER,*)'LST=',LST,' LEND=',LEND
	  WRITE(LUER,*)'N=',N,' NUION=',NUION,' Z=',Z
	  STOP
	END IF
C
C Compute scaled energy of ejected electron in rydbergs.
C
	U=(NU-NUION)/3.2882_LDP/Z/Z
C
C Check whether cross-sections need to be computed. Limit on
C on equality of U and USAV is not very stringent, but cross-section
C returned should have accuracy better than 1% and much better
C at threshold.
C
	IF( N .NE. NSAV .OR. .NOT. EQUAL(U,USAV,1.0E-03_LDP) )THEN
	  CALL GNTHYDL(SIGMA,G,N,U)
	  USAV=U
	  NSAV=N
	END IF
C
C Compute effective gaunt factor if state corresponds to more than
C one L value.
C
	HYDCROSSL=0.0_LDP
	GSUM=0.0_LDP
	DO L=LST,LEND,1
	  HYDCROSSL=HYDCROSSL+(2.0_LDP*L+1.0_LDP)*SIGMA(L)
	  GSUM=GSUM+(2.0_LDP*L+1.0_LDP)
	END DO
	HYDCROSSL=HYDCROSSL/GSUM
C
C Convert from Gaunt factor to cross section. Cross-section has
C been scaled by 10^10 (Thus at threshold they should be of
C the order of 10^{-8}.
C
	HYDCROSSL=2.815E-06_LDP*HYDCROSSL*( (Z/NU/NEF)**4 )*NU/N
C
	RETURN
	END
