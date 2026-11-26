!
! This option intepolates new populations between two depths (not
! inclusive). Interpolation is only peformed for those population where
! very the large corrections were made.
! 
! Preferred option is to read in STEQ_VALS.
!
          FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID1=2; ID2=ND/2-1; MAX_CHANGE=900.0D0
	  CALL GEN_IN(IT,'Iteration # (zero to exit)')
	  CALL GEN_IN(ID1,'Start of interpolating rangir (exclusive)')
	  CALL GEN_IN(ID2,'End of interpolating range (exclusive)')
	  CALL GEN_IN(MAX_CHANGE,'Interpolate values with FRACTIONAL correction > >')
!
	  CALL READ_CORRECTIONS(CORRECTIONS,POPS,ND,NT,NIT,LU_IN)
!
!	  OPEN(UNIT=LU_OUT,FILE='NEW_COR_SUM',STATUS='UNKNOWN',ACTION='WRITE')	
!	  WRITE(LU_OUT,'(A,4X,A,2X,A,6X,A,2(2X,A),4(8X,A))')' Depth','Ion','Ion','N(levs)',
!	1                 'Neg.','Pos.','Ave(-ve)','Min(-ve)','Ave(+ve)','Max(+ve)'
!	  WRITE(LU_OUT,'(34X,A,3X,A)')'Cnt','Cnt'
!
	  T2=(1.0D0-1.0D0/MAX_CHANGE)
	  ID=ID1
	  ID1=MIN(ID1,ID2); ID2=MAX(ID,ID2)
	  DO ID=ID1+1,ID2-1
	    DO I=1,NUM_IONS
	      DO IVAR=ION_INDEX(I),ION_INDEX(I)+NION(I)-1
	        IF(CORRECTIONS(IVAR,ID) .LT. -MAX_CHANGE
	1           .OR. CORRECTIONS(IVAR,ID) .GT. T2)THEN
	          LOW_ID=ID1
	          DO K=ID1+1,ID-1
	            IF(CORRECTIONS(IVAR,K) .GT. MAX_CHANGE .AND. 
	1              CORRECTIONS(IVAR,K) .LT. T2)LOW_ID=ID
	          END DO
	          UP_ID=ID2
	          DO K=ID+1,ID2-1
	            IF(CORRECTIONS(IVAR,K) .GT. MAX_CHANGE .AND. 
	1              CORRECTIONS(IVAR,K) .LT. T2)UP_ID=ID
	          END DO
	          T1=LOG(R(ID)/R(LOW_ID))/LOG(R(UP_ID)/R(LOW_ID))
	          POPS(IVAR,ID,IT)=EXP( T1*LOG(POPS(IVAR,UP_ID,IT)) +
	1                     (1.0D0-T1)*LOG(POPS(IVAR,LOW_ID,IT)) )
	        END IF
	      END DO
	    END DO
	  END DO
	  CLOSE(LU_OUT)
	  IREC=NIT			!IREC is updated on write
          IF(FDG_COUNTER .EQ. 1)NITSF=NITSF+1
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)' '
	  WRITE(6,*)' Broad summary of changes as function of ION written to NEW_COR_SUM'
	  WRITE(6,*)' Values listed only refer to changes exceeding MAX_CHANGE'
	  WRITE(6,*)' '
	  WRITE(6,*)'Corrections written to SCRTEMP as new iteration.'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'//RED_PEN
	  TMP_STR=' '; CALL GEN_IN(TMP_STR,'Hit any character to continue')
	  WRITE(6,*)DEF_PEN
	  GOTO 200
