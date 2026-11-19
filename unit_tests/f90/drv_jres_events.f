! drv_jres_events.f
!
! Manufactured driver for the J_RESCUE_EVENTS writer (obs/
! jres_validity.f).  Emits one NES_SAN event, one NRESC event in a
! second section, and the end-of-render JRESC_TOT, into a scratch
! J_RESCUE_EVENTS file in the run directory.  The Python test
! parses the records for positional-field sanity.
!
	PROGRAM DRV_JRES_EVENTS
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER LUER
	EXTERNAL JRESC_EV,JRESC_TOT
	LUER=6
	CALL JRESC_EV('NES_SAN',3,7,2.5E+01_LDP,'CONTINUUM',.FALSE.,
	1   'MOM_J_CMF_V8','MOM_J_CMF_V8',0,.FALSE.,3,0,12,
	2   1.0E-02_LDP,3.0E-01_LDP,2.0E+05_LDP,1.0E-03_LDP,
	3   1.0E+10_LDP,.TRUE.,4.0E+04_LDP,0.0_LDP,1.0E+01_LDP,
	4   5.0E+02_LDP,-5.0E-01_LDP,-1.5E+00_LDP,
	5   0.0_LDP,0.0_LDP,7.0E-05_LDP,9.0E+00_LDP,2.0E-05_LDP,LUER)
	CALL JRESC_EV('NRESC',9,2,3.0E+01_LDP,'LINE',.FALSE.,
	1   'FG_J_CMF_V12','LAM_ES_TC',14,.TRUE.,3,0,5,
	2   8.0E+06_LDP,1.1E+00_LDP,0.0_LDP,1.0E-03_LDP,1.0E+10_LDP,
	3   .TRUE.,8.0E+06_LDP,3.0E-01_LDP,2.0E+01_LDP,6.0E+02_LDP,
	4   2.0E+00_LDP,4.0E+00_LDP,
	5   7.0E+06_LDP,4.0E-01_LDP,1.0E-04_LDP,8.0E+00_LDP,3.0E-05_LDP,
	6   LUER)
	CALL JRESC_TOT()
	STOP
	END
