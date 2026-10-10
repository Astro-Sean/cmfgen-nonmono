!
! Routine to determine the most importan collision cooling lines. The flux in the
! line is due to collisons only (i.e., other processes like fluoresence are ignored) and
! we assume that all electrons eventually cascade to the ground state.
!
! As a check we output the total collsion cooling rate (ignoring ioizations).
!
	SUBROUTINE WR_COL_LINES(OMEGA,AXzV,XzV,XZVLTE,W_XzV,GXzV,EDGE,LEV_NAME,N,DESC,
	1                       R,V,SIGMA,DEPTH_INDX,NSTRONG_MAX,LU,NEW_FILE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created  29-Jan-2004 : Based on WR_COL_FRACTIONS
!
	INTEGER N
	INTEGER DEPTH_INDX
	INTEGER LU
	INTEGER NSTRONG_MAX
!
	REAL(KIND=LDP) XzV(N)
	REAL(KIND=LDP) XzVLTE(N)
	REAL(KIND=LDP) W_XzV(N)
	REAL(KIND=LDP) GXzV(N)
	REAL(KIND=LDP) EDGE(N)
	REAL(KIND=LDP) OMEGA(N,N)
	REAL(KIND=LDP) AXzV(N,N)
	CHARACTER(LEN=*) LEV_NAME(N),DESC
!
	REAL(KIND=LDP) R
	REAL(KIND=LDP) V
	REAL(KIND=LDP) SIGMA
	LOGICAL NEW_FILE
!
	REAL(KIND=LDP), PARAMETER :: OPLIN=2.6540081E+08_LDP
!
! Internal variables.
!
	REAL(KIND=LDP) FRACTION(N)
	REAL(KIND=LDP) COOL(N,N)
	REAL(KIND=LDP) NET(N)
!
	REAL(KIND=LDP) COOLING_CONSTANT
	REAL(KIND=LDP) COLLISION_COOLING
	REAL(KIND=LDP) TOTAL_COOLING
!
	REAL(KIND=LDP) ESCAPE_PROB
	REAL(KIND=LDP) TOTAL
	REAL(KIND=LDP) TAUL
	REAL(KIND=LDP) FREQ
	REAL(KIND=LDP) CHIL
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) GLDGU
!
	INTEGER LOCATION(2)
	INTEGER LEVEL
	INTEGER I
	INTEGER NL
	INTEGER NUP
	INTEGER IOS
	INTEGER LMAX
	CHARACTER(LEN=5) C1,C2
	CHARACTER(LEN=80) OUT_FORM
	CHARACTER(LEN=80) HEAD_FORM
	CHARACTER(LEN=60) TRANS_NAME
!
	INTEGER, PARAMETER :: IZERO=0
!
	OUT_FORM='MAIN_COOLING_LINES'
	IF(NEW_FILE)THEN
	  CALL GEN_ASCI_OPEN(LU,OUT_FORM,'UNKNOWN',' ','WRITE',IZERO,IOS)
	ELSE
	  CALL GEN_ASCI_OPEN(LU,OUT_FORM,'UNKNOWN','APPEND','WRITE',IZERO,IOS)
	END IF
	WRITE(LU,'(A)')' '
!
! Determine maximum name length. Used for formatting. Then determine the
! maximum numbers to output in one line, allowing for the level name.
!
	LMAX=0
	DO I=1,N
	  LMAX=MAX(LEN_TRIM(LEV_NAME(I)),LMAX)
	END DO
	WRITE(C1,'(I2)')LMAX+3;      C1=ADJUSTL(C1)
	WRITE(C2,'(I3)')2*LMAX+6;   C2=ADJUSTL(C2)
	OUT_FORM='(1X,A,T'//TRIM(C1)//',A,T'//TRIM(C2)//',2I5,2ES14.4)'
	HEAD_FORM='(1X,A,T'//TRIM(C1)//',A,T'//TRIM(C2)//',2A,7X,A,6X,A)'
	COOLING_CONSTANT=6.6261965E-12_LDP
!
	COOL=0.0_LDP
	NET=0.0_LDP
!
! Determine net number of collisional transitions into each level.
!
	DO LEVEL=2,N
	  DO NL=1,LEVEL-1
	    NET(LEVEL)=NET(LEVEL) + (OMEGA(NL,LEVEL)*XzV(NL)-OMEGA(LEVEL,NL)*XzV(LEVEL) )
	  END DO
	  DO NUP=LEVEL+1,N
	    NET(LEVEL)=NET(LEVEL) + (OMEGA(NUP,LEVEL)*XzV(NUP)-OMEGA(LEVEL,NUP)*XzV(LEVEL) )
	  END DO
	END DO
!
! Compute collisional cooling assuming ionization cooling is unimportant.
!
	T1=0.0_LDP
	DO LEVEL=1,N
	  DO NL=1,N
	    T1 = T1 + (EDGE(NL)-EDGE(LEVEL))*OMEGA(NL,LEVEL)*XzV(NL)
	  END DO
	END DO
	COLLISION_COOLING=T1*COOLING_CONSTANT
!
! We do this loop backwards as we need to allow for the cascades to the
! lower levels.
!
	DO NUP=N,2,-1
!
! Compute the fraction that decay to each lower level. To do this will first
! need to compute the total number of decays.
!
	  TOTAL=0.0_LDP
	  FRACTION=0.0_LDP
	  DO NL=1,NUP-1
	    FREQ=EDGE(NL)-EDGE(NUP)
	    IF(FREQ .NE. 0.0_LDP .AND. AXzV(NL,NUP) .NE. 0.0_LDP)THEN
              T1=W_XzV(NUP)/W_XzV(NL)
              GLDGU=GXzV(NL)/GXzV(NUP)
              CHIL=OPLIN*AXzV(NL,NUP)*(T1*XzV(NL)-GLDGU*XzV(NUP))
	      TAUL=CHIL*R*2.998E-10_LDP/FREQ/V
	      ESCAPE_PROB=(1.0_LDP-EXP(-TAUL))/TAUL
	      FRACTION(NL)=AXzV(NUP,NL)*ESCAPE_PROB
	      TOTAL=TOTAL+FRACTION(NL)
	    END IF
	  END DO
	  IF(TOTAL .EQ. 0.0_LDP)TOTAL=1.0_LDP
	  FRACTION=FRACTION/TOTAL
!
! Work out cascades.
!
	  DO NL=1,NUP-1
	    FREQ=EDGE(NL)-EDGE(NUP)
	    IF(FREQ .NE. 0.0_LDP .AND. AXzV(NL,NUP) .NE. 0.0_LDP)THEN
	      NET(NL)=NET(NL)+NET(NUP)*FRACTION(NL)
	    END IF
	  END DO
!
! Can now work out cooling lines:
!
	  DO NL=1,NUP-1
	    FREQ=EDGE(NL)-EDGE(NUP)
	    IF(FREQ .NE. 0.0_LDP .AND. AXzV(NL,NUP) .NE. 0.0_LDP)THEN
	      COOL(NUP,NL)=FRACTION(NL)*(EDGE(NL)-EDGE(NUP))*COOLING_CONSTANT*NET(NUP)
	    END IF
	  END DO
!
	END DO
!
	T1=0.0_LDP
	DO I=1,N
	  T1=T1+COOLING_CONSTANT*EDGE(NL)*(XzV(NL)-XzVLTE(NL))*OMEGA(NL,NL)
	END DO
!
! Find 10 strongest cooling lines
!
	TOTAL_COOLING=SUM(COOL)
	WRITE(LU,'(A)')' '
	WRITE(LU,'(3A,I3)')' Strongest cooling lines for ',TRIM(DESC),' at ',DEPTH_INDX
	WRITE(LU,'(/,A,ES14.4,A)')' The total collisonal line cooling rate is',TOTAL_COOLING,' ergs/cm^3/s'
	WRITE(LU,'(A,ES14.4,A)')  '         Collisional ionization cooling is',T1,' ergs/cm^3/s'
	WRITE(LU,'(A,ES14.4,A)')  '              Total collisional cooling is',COLLISION_COOLING,' ergs/cm^3/s'
	WRITE(LU,'(A)')' '
!
	WRITE(6,'(A)')' '
	WRITE(6,'(/,A,ES14.4,A)')' The total collisonal line cooling rate is',TOTAL_COOLING,' ergs/cm^3/s'
	WRITE(6,'(A,ES14.4,A)')  '         Collisional ionization cooling is',T1,' ergs/cm^3/s'
	WRITE(6,'(A,ES14.4,A)')  '              Total collisional cooling is',COLLISION_COOLING,' ergs/cm^3/s'
	WRITE(6,'(A)')' '
!
	WRITE(LU,HEAD_FORM)'Upper','Lower','  NUP','   NL','Cooling','Fraction'
	DO I=1,NSTRONG_MAX
	  LOCATION=MAXLOC(COOL)
	  NL=LOCATION(2); NUP=LOCATION(1)
	  WRITE(LU,OUT_FORM)TRIM(LEV_NAME(NUP)),TRIM(LEV_NAME(NL)),NUP,NL,
	1                      COOL(NUP,NL),COOL(NUP,NL)/TOTAL_COOLING
	  IF(COOL(NUP,NL)/TOTAL_COOLING .LT. 0.01_LDP)EXIT
	  COOL(NUP,NL)=0.0_LDP			!or use MASK option.
	END DO
!
	RETURN
	END
