!
! Routine to output collisional data in 132 column format. Only 3 significant
! digits are output. OPTION is presently not utilized, while DESC is
! used to detrmine the file name.
!
	SUBROUTINE WR_COL_SL(OMEGA,XzV,XZVLTE,EDGE,LEV_NAME,N,DESC,LU,LEVEL,NL,DEPTH)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered  24-May-2020 : When OPTION=NET_RATES we write out to XzV_TOTR_COL the net rate,
!                          from all levels, into each leval
! Created  29-Jan-2004 : Based on WR_COL
!
	INTEGER N,NL,LU,DEPTH
	REAL(KIND=LDP) XzV(N)
	REAL(KIND=LDP) XzVLTE(N)
	REAL(KIND=LDP) EDGE(N)
	REAL(KIND=LDP) OMEGA(N,N)
	CHARACTER*(*) LEV_NAME(N),DESC
	INTEGER LEVEL(NL)
!
! Internal variables.
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) SUM
	REAL(KIND=LDP) RATE_IN,RATE_OUT
	REAL(KIND=LDP) TOT_IN,TOT_OUT
	REAL(KIND=LDP) NET_IN
	INTEGER I,J,K,L,M,IOS
	INTEGER N_PER_LINE,LMAX,LIM,ST_POS
	INTEGER LEV
	LOGICAL MEM_SAME_SL(NL)
	CHARACTER*80 FILENAME
!
	INTEGER, PARAMETER :: IZERO=0
!
	FILENAME=TRIM(DESC)//'_LEV_COL'
	CALL GEN_ASCI_OPEN(LU,FILENAME,'UNKNOWN','APPEND',' ',IZERO,IOS)
!
! Determine maximum name length. Used for formatting. Then determine the
! maximum numbers to output in one line, allowing for the level name.
!
! Rates out of level
!
	DO J=1,NL
	  IF(LEVEL(J) .LE. 0)EXIT
	  WRITE(LU,'(/,A)')' Level of interest is '//TRIM(LEV_NAME(LEVEL(J)))
	END DO
	WRITE(LU,'(A,3X,I5,/)')' Depth is:',DEPTH
!
	MEM_SAME_SL=.FALSE.
	DO I=1,N
	  DO J=1,NL
	    IF(I .EQ. LEVEL(J))THEN
	      MEM_SAME_SL(I)=.TRUE.
              EXIT
	    END IF
	  END DO
	END DO
!
! We exclude couplings between levels in the SL.
!
	WRITE(LU,'(8X,A,3X,A,2(6X,A),T60,A)')'Net in','% in',' Rate in','Rate out','Level'
	TOT_IN=0.0_LDP; TOT_OUT=0.0_LDP
	DO I=1,N
	  RATE_IN=0.0_LDP; RATE_OUT=0.0_LDP
	  DO J=1,NL
	    LEV=LEVEL(J)
	    IF(LEV .LE. 0)EXIT
	    IF(I .EQ. LEV)THEN
	      RATE_IN=RATE_IN+OMEGA(I,I)*XzVLTE(I)
	      RATE_OUT=RATE_OUT+OMEGA(LEV,I)*XzV(LEV)
	    ELSE IF(.NOT. MEM_SAME_SL(I))THEN
	      RATE_IN=RATE_IN+OMEGA(I,LEV)*XzV(I)
	      RATE_OUT=RATE_OUT+OMEGA(LEV,I)*XzV(LEV)
	    END IF
	  END DO
	  NET_IN=RATE_IN-RATE_OUT
	  TOT_IN=TOT_IN+RATE_IN
	  TOT_OUT=TOT_OUT+RATE_OUT
	  T1=(RATE_IN+RATE_OUT)/100.0_LDP
!	  WRITE(LU,'(ES14.4,F8.3,2ES14.4,T60,A)')NET_IN,NET_IN/T1,RATE_IN,RATE_OUT,TRIM(LEV_NAME(I))
	  WRITE(LU,'(ES14.4,F8.3,2ES14.4,T60,A)')NET_IN,T1,RATE_IN,RATE_OUT,TRIM(LEV_NAME(I))
	END DO
	FLUSH(LU)
!
	DO I=6,LU,LU-6
	  WRITE(I,'(A)')' '
	  WRITE(I,'(ES14.4,T20,A)')TOT_IN,'Total col rate  in'
	  WRITE(I,'(ES14.4,T20,A)')TOT_OUT,'Total col rate out'
	  WRITE(I,'(ES14.4,T20,A)')TOT_IN-TOT_OUT,'Net col rate in'
	  T1=100.0_LDP*(TOT_IN-TOT_OUT)/(TOT_IN+TOT_OUT)
	  WRITE(I,'(ES14.4,T20,A)')T1,'% Net col rate in (i.e, 100*NET/(TOT_IN+TOT_OUT)'
	END DO
	WRITE(6,'(/,2A)')' Collision rates written (appended) to: ',TRIM(FILENAME)
	CLOSE(LU)
!
	RETURN
	END
