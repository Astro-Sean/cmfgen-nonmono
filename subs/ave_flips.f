!
! Subroutine to create an averaged set of populations from the data in
! SCRTEMP. Only populations whose direction of changes has oscilated
! NUM_OSC times are effected. NUM_OSC can be 1 or larger.
!
	SUBROUTINE AVE_FLIPS(NEW_POPS,POPS,NT,ND,CUR_IT_NO,NUM_OSC,SIMPLE,AV_DONE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Aletered 06-Jul-2017: Routine returns AV_DONE=.FALSE is ICOUNT (number of flips) is zero.
!
	INTEGER NT
	INTEGER ND
	REAL(KIND=LDP) NEW_POPS(NT*ND)
	REAL(KIND=LDP) POPS(NT*ND)
!
	INTEGER CUR_IT_NO
	INTEGER NUM_OSC
	INTEGER IFLAG
	LOGICAL SIMPLE
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables
!
	REAL(KIND=LDP) OLD_POPS(NT*ND,NUM_OSC+2)
	REAL(KIND=LDP) R_TMP(ND)
	REAL(KIND=LDP) V_TMP(ND)
	REAL(KIND=LDP) SIGMA_TMP(ND)
	REAL(KIND=LDP) R
	REAL(KIND=LDP) T1,T2
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER NTD
	INTEGER ICOUNT
	INTEGER I,J
	INTEGER IST
	INTEGER NITSF
	INTEGER IREC
	INTEGER LST_NG
	INTEGER, PARAMETER :: RITE_N_TIMES=1
	INTEGER, PARAMETER :: LUSCR=10
!
	LOGICAL NEWMOD
	LOGICAL AV_DONE
	LOGICAL WRITE_RVSIG
!
	AV_DONE=.FALSE.
	NEW_POPS=POPS
	IF(CUR_IT_NO .LE. 3+NUM_OSC)RETURN
	NTD=NT*ND
	LUER=ERROR_LU()
!
! Read in the last NSTEP estimates of the populations, as output to SCRTEMP.
! NEW_POPS is used for temporary storage.
!
        NEWMOD=.FALSE.
	IST=CUR_IT_NO-NUM_OSC-1
	DO I=IST,CUR_IT_NO
	  IREC=I
          CALL SCR_READ_V2(R_TMP,V_TMP,SIGMA_TMP,NEW_POPS,
	1            IREC,NITSF,RITE_N_TIMES,LST_NG,
	1            WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  IF(I .NE. IREC .OR. NEWMOD)THEN
	    WRITE(LUER,*)'Error reading correct record in SCRTEMP'
	    RETURN
	  END IF
	  OLD_POPS(:,1+(CUR_IT_NO-I))=NEW_POPS
	END DO
!
! Average those populations which have oscillated. We use either a geometric
! series to get the best estimate, or simple average the result.
!
	NEW_POPS=POPS
	R=1.0_LDP
	ICOUNT=0
	DO I=1,NT*ND
	  DO J=2,NUM_OSC
	    T1=OLD_POPS(I,J)-OLD_POPS(I,J+1)
	    T2=OLD_POPS(I,J+1)-OLD_POPS(I,J+2)
	    IF(T1*T2 .GE. 0.0_LDP)GOTO 100
	  END DO
	  T1=OLD_POPS(I,1)-OLD_POPS(I,2)
	  T2=OLD_POPS(I,2)-OLD_POPS(I,3)
	  IF(T1*T2 .GE. 0.0_LDP .AND. ABS(T1) .LT. ABS(T2) )GOTO 100
	  IF(.NOT. SIMPLE)R=ABS(T1/T2)
	  NEW_POPS(I)=POPS(I)-T1*R/(1.0_LDP+R)
	  ICOUNT=ICOUNT+1
100	  CONTINUE
	END DO
	IF(ICOUNT .GT. 0)AV_DONE=.TRUE.
	WRITE(LUER,*)'Number of averages performed in AVE_FLIPS is',ICOUNT
	WRITE(LUER,*)'      Maximum number of possible averages is',NT*ND
!
	RETURN
	END
