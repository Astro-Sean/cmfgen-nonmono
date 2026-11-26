!
! Routine to output collisional data in 132 column format. Only 3 significant
! digits are output. OPTION is presently not utilized, while DESC is
! used to detrmine the file name.
!
	SUBROUTINE WR_COL_RATES(OMEGA,XzV,XZVLTE,EDGE,LEV_NAME,N,DESC,LU,OPTION)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered  24-May-2020 : When OPTION=NET_RATES we write out to XzV_TOTR_COL the net rate,
!                          from all levels, into each leval
! Created  29-Jan-2004 : Based on WR_COL
!
	INTEGER N,LU
	REAL(KIND=LDP) XzV(N)
	REAL(KIND=LDP) XzVLTE(N)
	REAL(KIND=LDP) EDGE(N)
	REAL(KIND=LDP) OMEGA(N,N)
	CHARACTER*(*) LEV_NAME(N),OPTION,DESC
!
! Internal variables.
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) SUM
	INTEGER I,J,K,L,M,IOS
	INTEGER N_PER_LINE,LMAX,LIM,ST_POS
	CHARACTER*80 FORM
	CHARACTER*132 STRING
	CHARACTER*30 TMP_NAME
!
! Number of digits to output OMEGA. Must include spaces used to separate
! number from proceeding number. Must be the same as the number X in
! EX.3
!
	INTEGER, PARAMETER :: DIG_PER_NUM=11
	INTEGER, PARAMETER :: IZERO=0
!
	WRITE(LU,'(A)')' '
	IF( TRIM(OPTION) .EQ. 'NET_RATES')THEN
	  WRITE(LU,'(A)')'  Net collision rate bewteen upper and lower levels'
	  WRITE(LU,'(A)')'  Rate is +ve when net change into lower level'
	  WRITE(LU,'(A)')'  I,I is the net collisional recombination rate.'
	  FORM=TRIM(DESC)//'_NR_COL'
	ELSE IF( TRIM(OPTION) .EQ. 'COOL_RATES')THEN
	  WRITE(LU,'(A)')'  Net electron cooling rate (cgs units) bewteen upper and lower levels'
	  WRITE(LU,'(A)')'  Rate is +ve when cooling thermal electrons'
	  WRITE(LU,'(A)')'  I,I is the net ionzation/recombination cooling rate.'
	  FORM=TRIM(DESC)//'_CR_COL'
	ELSE
	  WRITE(LU,'(A)')'  collision rate down to lower level from upper leve'
	  WRITE(LU,'(A)')'  I,I is the collisional ionization rate.'
	  FORM=TRIM(DESC)//'_DR_COL'
	END IF
	WRITE(LU,'(A)')' '
!
	CALL GEN_ASCI_OPEN(LU,FORM,'UNKNOWN',' ',' ',IZERO,IOS)
!
! Determine maximum name length. Used for formatting. Then determine the
! maximum numbers to output in one line, allowing for the level name.
!
	LMAX=0
	DO I=1,N
	  LMAX=MAX(LEN_TRIM(LEV_NAME(I)),LMAX)
	END DO
	N_PER_LINE=MIN(N,(132-LMAX-2)/DIG_PER_NUM)
!
	DO M=1,(N+N_PER_LINE-1)/N_PER_LINE
	  LIM=MIN(M*N_PER_LINE,N)
!
! Determine header string. Only the last DIG_PER_NUM-2 characters in the
! name are output.
!
	  ST_POS=LMAX+1
	  STRING=' '
	  DO J=(N_PER_LINE*(M-1)+1),LIM
            K=LEN_TRIM(LEV_NAME(J))
	    IF(K .LT. DIG_PER_NUM-2)THEN
	      TMP_NAME=' '
	      TMP_NAME(DIG_PER_NUM+1-K:DIG_PER_NUM)=LEV_NAME(J)
	    ELSE
	      TMP_NAME='  '//LEV_NAME(J)(K-DIG_PER_NUM+3:K)
	    END IF
	    STRING(ST_POS:)=TMP_NAME
	    ST_POS=ST_POS+DIG_PER_NUM
	  END DO
	  WRITE(LU,*)' '
	  WRITE(LU,'(A)')TRIM(STRING)
!
! We only ouput OMEGA(I,J) for J>=I. K and L help set up the appropriate
! format to do this.  For the INTEL compiler, K must be at least 1
! (altered 24-Jun-2003)
!
	  K=1
	  L=N_PER_LINE
	  IF( TRIM(OPTION) .EQ. 'NET_RATES')THEN
	    DO I=1,MIN(N_PER_LINE*M,N)
	      IF(N_PER_LINE*M-I .LT. N_PER_LINE-1)THEN
                K=K+DIG_PER_NUM
	        L=L-1
	      END IF
	      WRITE(FORM,'(A,I4.4,A,I4.4,A,I4.4,A)')
	1                   '(A',LMAX,',',K,'X,1P,',L,'E11.3)'
	      IF(MAX((N_PER_LINE*(M-1)+1),I) .EQ. I)THEN
	        WRITE(LU,FORM)
	1           LEV_NAME(I),OMEGA(I,I)*(XzVLTE(I)-XzV(I)),
	1              ( (OMEGA(J,I)*XzV(J)-OMEGA(I,J)*XzV(I)), J=MAX((N_PER_LINE*(M-1)+1),I)+1,LIM)
	      ELSE
	        WRITE(LU,FORM)
	1         LEV_NAME(I),( (OMEGA(J,I)*XzV(J)-OMEGA(I,J)*XzV(I)),
	1                J=MAX((N_PER_LINE*(M-1)+1),I),LIM )
	      END IF
	    END DO
!	
	  ELSE IF( TRIM(OPTION) .EQ. 'COOL_RATES')THEN
	    T1=6.626E-12_LDP
	    DO I=1,MIN(N_PER_LINE*M,N)
	      IF(N_PER_LINE*M-I .LT. N_PER_LINE-1)THEN
                K=K+DIG_PER_NUM
	        L=L-1
	      END IF
	      WRITE(FORM,'(A,I4.4,A,I4.4,A,I4.4,A)')
	1                   '(A',LMAX,',',K,'X,1P,',L,'E11.3)'
	      IF(MAX((N_PER_LINE*(M-1)+1),I) .EQ. I)THEN
	        WRITE(LU,FORM)
	1           LEV_NAME(I),OMEGA(I,I)*(XzV(I)-XzVLTE(I))*EDGE(I)*T1,
	1              ( (OMEGA(I,J)*XzV(I)-OMEGA(J,I)*XzV(J))*(EDGE(I)-EDGE(J))*T1,
	1               J=MAX((N_PER_LINE*(M-1)+1),I)+1,LIM)
	      ELSE
	        WRITE(LU,FORM)
	1         LEV_NAME(I),( (OMEGA(I,J)*XzV(I)-OMEGA(J,I)*XzV(J))*(EDGE(I)-EDGE(J))*T1,
	1                J=MAX((N_PER_LINE*(M-1)+1),I),LIM )
	      END IF
	    END DO
	  ELSE
	    DO I=1,MIN(N_PER_LINE*M,N)
	      IF(N_PER_LINE*M-I .LT. N_PER_LINE-1)THEN
                K=K+DIG_PER_NUM
	        L=L-1
	      END IF
	      WRITE(FORM,'(A,I4.4,A,I4.4,A,I4.4,A)')
	1                   '(A',LMAX,',',K,'X,1P,',L,'E11.3)'
	      WRITE(LU,FORM)LEV_NAME(I),(OMEGA(J,I)*XzV(J),J=MAX((N_PER_LINE*(M-1)+1),I),LIM)
	    END DO
	  END IF
	END DO
!
	IF(TRIM(OPTION) .EQ. 'COOL_RATES')THEN
	  T1=6.626E-12_LDP
	  SUM=0.0_LDP
	  DO I=1,N
	    SUM=SUM+OMEGA(I,I)*(XzV(I)-XzVLTE(I))*EDGE(I)*T1
	    DO J=I+1,N
	      SUM=SUM+(OMEGA(I,J)*XzV(I)-OMEGA(J,I)*XzV(J))*(EDGE(I)-EDGE(J))*T1
	    END DO
	 END DO
	 WRITE(LU,'(A)')' '
	 WRITE(LU,'(A,ES14.4)')'Total colissional cooling rate is:',SUM
!
	ELSE IF(TRIM(OPTION) .EQ. 'NET_RATES')THEN
	  SUM=0.0_LDP
	  DO I=1,N
	    SUM=SUM+OMEGA(I,I)*(XzV(I)-XzVLTE(I))
	   END DO
	 WRITE(LU,'(A)')' '
	 WRITE(LU,'(A,ES14.4)')'Net ionization rate is:',SUM
	END IF
!
	CLOSE(LU)
!
! Determine the NET collision rate into each level. We can simply sum
! when I=J as the rate will be zero.
!
	IF( TRIM(OPTION) .EQ. 'NET_RATES')THEN
	  FORM=TRIM(DESC)//'_TOTR_COL'
	  CALL GEN_ASCI_OPEN(LU,FORM,'UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(FORM,'(I2.2)')LMAX+5
	  FORM='(A,T'//FORM(1:2)//',ES16.8)'
	  DO I=1,N
	    T1=OMEGA(I,I)*(XzVLTE(I)-XzV(I))
	    DO J=1,N
	      T1=T1+(OMEGA(J,I)*XzV(J)-OMEGA(I,J)*XzV(I))
	    END DO
	    WRITE(LU,TRIM(FORM))LEV_NAME(I),T1
	  END DO
	  CLOSE(LU)
	END IF
!
	RETURN
	END
