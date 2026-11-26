C
C Function to return the wavelength of a transition which is specified
C by a frequency (10^15 Hz). The wavelength is returned in Angstroms.
C
C If the wavelenth is less than 2000Ang, it is a vacuum wavelength,
C otherwise it is in air.
C
	FUNCTION LAMVACAIR(FREQ)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) LAMVACAIR,FREQ,T1
C
	LAMVACAIR=2997.92458_LDP/FREQ
C
	IF(LAMVACAIR .GT. 1999.352_LDP)THEN
	  T1=1.0E+08_LDP/(LAMVACAIR*LAMVACAIR)
	  LAMVACAIR=LAMVACAIR/(  1.0_LDP+
	1               1.0E-07_LDP*( 643.28_LDP+
	1              294981.0_LDP/(146.0_LDP-T1)+2554.0_LDP/(41.0_LDP-T1) )  )
	END IF
C
	RETURN
	END
C
	FUNCTION LAM_AIR(LAM)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) LAM_AIR,LAM,T1
C
	LAM_AIR=LAM
C
	IF(LAM_AIR .GT. 1999.352_LDP)THEN
	  T1=1.0E+08_LDP/(LAM_AIR*LAM_AIR)
	  LAM_AIR=LAM_AIR/(  1.0_LDP+
	1               1.0E-07_LDP*( 643.28_LDP+
	1              294981.0_LDP/(146.0_LDP-T1)+2554.0_LDP/(41.0_LDP-T1) )  )
	END IF
C
	RETURN
	END
