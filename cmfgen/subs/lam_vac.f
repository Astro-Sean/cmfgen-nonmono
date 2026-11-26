C
C Routine to compute the vacuum wavelenth from a wavelength in AIR. The
C input and output units are Angstroms.
C
C Formulae is only valid for LAM(AIR) > than 2000Ang.
C
C Results checked against table in Allen: Astrophysical quantities.
C
	FUNCTION LAM_VAC(LAM_AIR)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Created 11-Apr-1997 : Based on LAMVACAIR
C
	REAL(KIND=LDP) LAM_VAC			!Vacuum wavelength in Ang.
	REAL(KIND=LDP) LAM_AIR			!Air wavelength in Ang.
	REAL(KIND=LDP) T1
	INTEGER I
C
	IF(LAM_AIR .LT. 1999.352_LDP)THEN
	  WRITE(6,*)' ERROR --- in LAM_VAC'
	  WRITE(6,*)' Formulae for conversion of air to vacuum wavelenths',
	1                ' invalid below 2000Ang'
	  STOP
	END IF
C
C Need to iterate, as T1 is unknown vacuum wavelength. Because difference
C between LAM_VAC and LAM_AIR is small, convergence is VERY rapid.
C
	LAM_VAC=LAM_AIR
	DO I=1,3
	  T1=1.0E+08_LDP/(LAM_VAC*LAM_VAC)
	  LAM_VAC=LAM_AIR*(  1.0_LDP+
	1               1.0E-07_LDP*( 643.28_LDP+
	1              294981.0_LDP/(146.0_LDP-T1)+2554.0_LDP/(41.0_LDP-T1) )  )
	END DO
C
	RETURN
	END
