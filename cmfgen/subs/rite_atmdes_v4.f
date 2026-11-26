C
C Routine to output basic descriptors for each ionic species. These include
C number of levels, core charge, and statistical weight of ion. Species
C not include in model are indicated by *****.
C
	SUBROUTINE RITE_ATMDES_V4(CI_PRES,
	1                      NCI_S,ZCI,EQCI,
	1                      CI_LEVNAME_F,NCI_F,GIONCI_F,
	1                      N_CI_PHOT,AT_NO_CARB,MASS_CARB,
	1                      ID,ISPEC,LUMOD,DESC)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 10-Nov-2023 --- Changed from V3. Add ID and ISPEC to call.
C Altered 13-Oct-2021 --- Changed format so masses > 100 can be written correctly.
C Altered 29-May-1996 --- <> deleted from FORMAT 100.
C Altered 22-Sep-1990 --- STRING installed, varaible format <> deleted for
C                         CRAY compatibility.
C Altered 20-Sep-1990 --- EQCI installed in call. Note that this routine is
C                         also called by OBSCIV and PROCMF4TH. In these two
C                         routines EQCI is undefined.
C Created 18-Aug-1989.
C
	LOGICAL CI_PRES
	INTEGER NCI_S,NCI_F
	INTEGER LUMOD,EQCI
	INTEGER N_CI_PHOT
	INTEGER ID
	INTEGER ISPEC
	REAL(KIND=LDP) ZCI,GIONCI_F
	REAL(KIND=LDP) MASS_CARB
	REAL(KIND=LDP) AT_NO_CARB
	CHARACTER(LEN=*) CI_LEVNAME_F(NCI_F)
	CHARACTER(LEN=*) DESC
C
	REAL(KIND=LDP) AT_NO_SAVE
	SAVE AT_NO_SAVE
	DATA AT_NO_SAVE/0/
C
	INTEGER L,ICHRLEN
	CHARACTER STAR*4
C
	L=ICHRLEN(DESC)
	IF(AT_NO_SAVE .NE. AT_NO_CARB)WRITE(LUMOD,'(A)')' '
	IF(CI_PRES)THEN
	  STAR='    '
	  WRITE(LUMOD,100)DESC(1:L),STAR,NCI_F,NCI_S,
	1                ZCI,GIONCI_F,MASS_CARB,AT_NO_CARB,
	1                N_CI_PHOT,EQCI,ISPEC,ID,TRIM(CI_LEVNAME_F(NCI_F))
	ELSE
	  STAR='****'
	  WRITE(LUMOD,100)DESC(1:L),STAR
	END IF
	AT_NO_SAVE=AT_NO_CARB
C
100	FORMAT(1X,A,T12,A4,:2(2X,I4):,4(3X,F5.1,1X),4X,I2,4X,I4,3X,I4,2X,I4,4X,A)
C
	RETURN
	END

	SUBROUTINE RITE_ATMHD_V4(LUMOD)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER LUMOD
C
	WRITE(LUMOD,100)
100	FORMAT(1X,'Species',T12,6X,' N_F',2X,' N_S ',4X,' Z ',
	1          5X,'Gion',5X,'Mass',6X,'At.',
	1          4X,'No.',4X,'Eqn.',3X,'ISPEC',4X,'ID',4X,'Last level')
	WRITE(LUMOD,110)
110	FORMAT(40X,'(gs)',15X,'No.',3X,'Phot.',4X,'No.')
C
	RETURN
	END
