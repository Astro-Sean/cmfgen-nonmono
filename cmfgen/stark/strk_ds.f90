!
!-----Isolated and multi component HeI and NIT line profiles
!     Data are input from file HEI_OTHERSTARK
!     See below for required data.
!
	SUBROUTINE STRK_DS(PRO,FREQ,NF,ED,TE,VTURB,ND,NU_ZERO,    &
	                         SPECIES,NL,NUP,AMASS,            &
	                         FILENAME,LU_STRK)
	USE SET_KIND_MODULE
!
! Altered 09-Oct-2002 : Bug fixed (effect compiler dependent).
!                         STKTB_PRES was not being initialized.
! Created 20-Sep-2000 : Based on STRK_DS and STRK_DIMITR for F. Najarro.
!
	IMPLICIT NONE
!
	INTEGER*4 NF		!Number of frequencies
	INTEGER*4 ND		!Number of depth points
	INTEGER*4 LU_STRK	!Input logical unit
	REAL(KIND=LDP) FREQ(NF)		!Frequncy in units of 10^15 Hz
	REAL(KIND=LDP) PRO(ND,NF)	!Profile
	REAL(KIND=LDP) TE(ND)		!Electron temperature (10^4 K)
	REAL(KIND=LDP) ED(ND)		!Electron density (cgs units)
	REAL(KIND=LDP) VTURB(ND)	!Turbulent velocity km/s
	REAL(KIND=LDP) NU_ZERO		!Line frequency (10^15 Hz)
	REAL(KIND=LDP) AMASS		!Species mass in AMU
!
	INTEGER*4 NL			!Lower transition level
	INTEGER*4 NUP			!Upper transition level
	CHARACTER(LEN=*) SPECIES	!Species (e.g. HeI, etc)
	CHARACTER(LEN=*) FILENAME	!
!
! Profile parameters:
!
	INTEGER*4, PARAMETER :: MXS=10
	INTEGER*4 STKTB_NTS
	INTEGER*4 STKTB_NL
	INTEGER*4 STKTB_NU
	REAL(KIND=LDP) STKTB_TS(MXS)
	REAL(KIND=LDP) STKTB_WS(MXS)
	REAL(KIND=LDP) STKTB_DS(MXS)
	REAL(KIND=LDP) STKTB_AS(MXS)
	LOGICAL STKTB_PRESS
	REAL(KIND=LDP) STKTB_DLP
	REAL(KIND=LDP) STKTB_ELE
	INTEGER*4 NL_RD
	INTEGER*4 NUP_RD
	INTEGER*4 NLSTR
!
	INTEGER*4 ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables:
!
	INTEGER*4 J,IB,IA,IDE,IOS
	REAL(KIND=LDP) DLS,FOS,FT,RFT,WT,RWT,RBHZ,RBA,TT,EE,VMOT,WF,Y,PP
	REAL(KIND=LDP) VM,RHOM,ALF,SIGMA,X,A,CON,DB,P,VA,VB,D,W
	REAL(KIND=LDP) TMP
	REAL(KIND=LDP) CLIGHT,SRT
	CHARACTER*80 STRING
	DATA CLIGHT/2.997925E18_LDP/
	DATA SRT/1.414213562_LDP/
!
	REAL(KIND=LDP) VOIGTN
	EXTERNAL VOIGTN
!
	FT=NU_ZERO*1.E15_LDP
	RFT=1._LDP/FT
!
! Obtain profile information: At present this is slow, since it may
! necessitate reading the entire file. We skipe over comments, indicated
! by a " ! " in the file.
!
	STKTB_PRESS=.FALSE.
	OPEN (UNIT=LU_STRK,FILE=FILENAME,FORM='FORMATTED',     &
	           STATUS='OLD',ACTION='READ')
	  DO WHILE( .NOT. STKTB_PRESS )
	    STRING='!'
	    DO WHILE(STRING(1:1) .EQ. '!')
	      READ(LU_STRK,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        WRITE(ERROR_LU(),*)'Error reading ',TRIM(FILENAME),'-- End of file'
	        STOP
	      END IF
	    END DO
	    STRING=ADJUSTL(STRING)
	    NLSTR=INDEX(STRING,' ')
	    READ(STRING(NLSTR:),*)NL_RD,NUP_RD,STKTB_NL,STKTB_NU
	    IF ( STRING(1:NLSTR-1) .EQ. SPECIES .AND.          &
                   NL_RD .EQ. NL .AND. NUP_RD .EQ. NUP ) THEN
	      STKTB_PRESS=.TRUE.
	    END IF
	    READ(LU_STRK,*) TMP,STKTB_DLP,STKTB_ELE,STKTB_NTS
	    READ(LU_STRK,*) ( STKTB_TS(J), J=1,STKTB_NTS )
	    READ(LU_STRK,*) ( STKTB_WS(J), J=1,STKTB_NTS )
	    READ(LU_STRK,*) ( STKTB_DS(J), J=1,STKTB_NTS )
	    READ(LU_STRK,*) ( STKTB_AS(J), J=1,STKTB_NTS )
	  END DO
	CLOSE(UNIT=LU_STRK)
!
! Now compute the profile.
!
!-----Doppler quantities					
!
	WT=CLIGHT*RFT
	RWT=FT/CLIGHT
!
DEPTH:	DO IDE=1,ND
	  TT=TE(IDE)*1.E4_LDP
	  EE=ED(IDE)
	  VMOT=12.85_LDP*SQRT(TE(IDE)/AMASS+(VTURB(IDE)/12.85_LDP)**2._LDP)/2.997925E5_LDP
	  VMOT=1._LDP/VMOT	! (CLIGHT/(Vth+Vtur)) with Amass
	  RBHZ=VMOT*RFT							
	  RBA=VMOT*RWT								
!
!-----Set up interpolation in T						
!
	  DO J=2,STKTB_NTS
	   IA=J
	   IF(STKTB_TS(J).GE.TT)EXIT
	  END DO
	  IB=IA-1
	  WF=(TT-STKTB_TS(IB))/(STKTB_TS(IA)-STKTB_TS(IB))
	  IF(TT.LT.STKTB_TS(IB)) THEN
	    WF = 0.0
	  END IF
	  IF(TT.GT.STKTB_TS(IA)) THEN
	    WF = 1.0
	  ENDIF
!
!-----Perturber quantities
!
	  Y=EE/STKTB_ELE
	  PP=EE
	  VM=8.78_LDP*(2.0_LDP*EE-PP)/(EE*SRT)
	  RHOM=0.62029_LDP/EE**0.333333_LDP
!
!-----Impact width width (A)
!
	  W=(WF*(STKTB_WS(IA)-STKTB_WS(IB))+STKTB_WS(IB))*Y
!
!-----Ratio impact shift/width
!
	  D=WF*(STKTB_DS(IA)-STKTB_DS(IB))+STKTB_DS(IB)
!
!-----Ion broadening parameters
!
	  ALF=(WF*(STKTB_AS(IA)-STKTB_AS(IB))+STKTB_AS(IB)) *Y**0.25_LDP
	  SIGMA=1.885E14_LDP*W*RHOM*VM*RWT*RWT
	  X=ALF*0.888889_LDP/SIGMA**0.666667_LDP
!
!-----Total width in Doppler units
!
	  A=W*(1.0_LDP+1.36_LDP*X)*RBA
!
!-----Total shift in Doppler units
!
	  DLS=W*D*(1.0_LDP+2.36_LDP*X/ABS(D))*RBA
	  FOS=1.0
!
!-----Satellite components
!							
	  X=FOS
	  CON=6.268773E-2_LDP*RBHZ/X
!
!-----Compute profile
!
	  DO J=1,NF
	    DB=-(FREQ(J)-NU_ZERO)*VMOT/NU_ZERO
	    VA=DB-DLS
	    VB=VA-STKTB_DLP
	    P=FOS*(8.0_LDP*VOIGTN(A,VA)+VOIGTN(A,VB))
	    PRO(IDE,J)=CON*P
	  END DO
  	END DO DEPTH
!
	RETURN
	END
