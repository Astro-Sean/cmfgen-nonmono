!
! Subroutine to measure the CHI^2 of an observation relative to a model.
! The slope across the band defining the lines is chosen to minimize chi^2.
! This is done using a direct numerical integration of the data.
! Routine is designed to be called in GRAMON_PGPLOT.
!
	SUBROUTINE DO_FILE_BALMER_V1(FILE_WITH_LINE_LIMS,IP_OBS)
	USE SET_KIND_MODULE
	USE MOD_CURVE_DATA
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	INTEGER IP_OBS
	CHARACTER(LEN=*) FILE_WITH_LINE_LIMS
!
	INTEGER, PARAMETER :: IONE=1
!
! Local parameters.
!
	INTEGER N_OMIT
	INTEGER NC_OMIT
	INTEGER OMIT_ST(20), OMIT_END(20)
	REAL*4  OLAM_ST(20), OLAM_END(20)
!
        REAL*4  C_OLAM_ST(20), C_OLAM_END(20)
        INTEGER C_OMIT_ST(20), C_OMIT_END(20)
!
	INTEGER NPIX
	INTEGER IST,IEND		!Line limits in pixel space
	REAL*4 LAM_ST,LAM_END           !Line limits in lambda space
	REAL(KIND=LDP) LAM_CENT			!Line centroid
	REAL(KIND=LDP) MEAN 			!Used computing line centroid
	REAL(KIND=LDP) EW_MOD,EW_OBS
        REAL(KIND=LDP) EW_OBS_OMIT,EW_MOD_OMIT
!
! Work variables
!
	REAL*4 T1,T2,T3,T4	!Work variable
	REAL*4 X1,X2		!Work variable
	REAL*4 D1,D2		!Work variable
	REAL*4 LAM
	REAL*4 O_DATA
	REAL*4 SLOPE,INTER,CONT
!
! Work variable used to compute slope of line that yields the lowest chi^2.
!
	REAL(KIND=LDP) SUM_OSQ, SUM_LOSQ, SUM_LSQ_OSQ
	REAL(KIND=LDP) SUM_MO, SUM_LMO
	REAL(KIND=LDP) A, B
	REAL(KIND=LDP) DET, DA, DB
	REAL(KIND=LDP) CHISQ,RAW_CHISQ,RED_CHISQ
!
! Mask is to omit regions containg weak lines in the Balmer wings.
!
	REAL*4, ALLOCATABLE :: OBS_DATA(:)
	REAL*4, ALLOCATABLE :: MOD_DATA(:)
	REAL*4, ALLOCATABLE :: MASK(:)
	REAL*4, ALLOCATABLE :: C_MASK(:)
!
	INTEGER, SAVE :: LU_LIMS=7
	INTEGER, SAVE :: LU_OUT=10
!
	INTEGER IP
	INTEGER IP_MOD
	INTEGER IP_OUT
	INTEGER I,J,K
	INTEGER NP
!
! External function.
!
	INTEGER GET_INDX_SP
	CHARACTER(LEN=30) UC
	EXTERNAL GET_INDX_SP
	EXTERNAL UC
!
	LOGICAL END_FILE
	LOGICAL FILE_PRES
!
	CHARACTER(LEN=80) OUT_FILE
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=80) TRANS_NAME
!
! Open output file. Data is appended if it alread exists.
!
	OUT_FILE='BALMER_DATA'
	CALL GEN_IN(OUT_FILE,'File to OUTPUT chi^2 and EW etc')
	INQUIRE(FILE=OUT_FILE,EXIST=FILE_PRES)
	IF(FILE_PRES)THEN
	  WRITE(6,*)'File already exists -- appending new data'
	  OPEN(UNIT=LU_OUT,FILE=OUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	ELSE
	  OPEN(UNIT=LU_OUT,FILE=OUT_FILE,STATUS='NEW',ACTION='WRITE')
	  CALL WRITE_BALMER_HEADER(LU_OUT)
	END IF
!
	OPEN(UNIT=LU_LIMS,STATUS='OLD',ACTION='READ',FILE=FILE_WITH_LINE_LIMS)
	DO WHILE(1 .EQ. 1)
	  DO WHILE(1 .EQ. 1)
	    READ(LU_LIMS,'(A)',END=1000)STRING
	    IF(STRING .NE. ' ' .AND. STRING(1:1) .NE. '!')EXIT
	  END DO
!
	  READ(STRING,*,IOSTAT=IOS)T1                 !This is the results.
	  IF(IOS .NE. 0)THEN
	     WRITE(6,*)'Error reading balmer line record'
	     WRITE(6,*)'STRING follows:'
	     WRITE(6,*)TRIM(STRING)
	  END IF
	  READ(LU_LIMS,*)LAM_ST,LAM_END,N_OMIT,NC_OMIT
	  IF(N_OMIT .NE. 0)READ(LU_LIMS,*)(OLAM_ST(I),OLAM_END(I),I=1,N_OMIT)
	  IF(NC_OMIT .NE. 0)READ(LU_LIMS,*)(C_OLAM_ST(I),C_OLAM_END(I),I=1,NC_OMIT)
	  DO IP_MOD=1,NPLTS
	    WRITE(6,*)'IP_MOD=',IP_MOD
	    IF(IP_MOD .EQ. IP_OBS)THEN
	    ELSE
!
! Interpolate data onto OBSERVATIONAL grid.
!
	      IF(ALLOCATED(MOD_DATA))THEN
	        DEALLOCATE(MOD_DATA,OBS_DATA,MASK,C_MASK)
	      END IF
	      NP=NPTS(IP_OBS)
	      ALLOCATE(MOD_DATA(NP),MASK(NP),C_MASK(NP),OBS_DATA(NP))
	      CALL MON_INTERP_SP(MOD_DATA,NP,IONE,CD(IP_OBS)%XVEC,NP,CD(IP_MOD)%DATA,NPTS(IP_MOD),CD(IP_MOD)%XVEC,NPTS(IP_MOD))
	      OBS_DATA=CD(IP_OBS)%DATA
	      WRITE(6,*)'Interpolated data'
!
! This defines the full extent of the line.
!
	      IST=GET_INDX_SP(LAM_ST,CD(IP_OBS)%XVEC,NPTS(IP_OBS))
	      IEND=GET_INDX_SP(LAM_END,CD(IP_OBS)%XVEC,NPTS(IP_OBS))
!
	      NPIX=0; MASK=1.0D0; C_MASK=1.0D0
	      DO J=1,N_OMIT
	        OMIT_ST(J)=GET_INDX_SP(OLAM_ST(J),CD(IP_OBS)%XVEC,NPTS(IP_OBS))
	        OMIT_END(J)=GET_INDX_SP(OLAM_END(J),CD(IP_OBS)%XVEC,NPTS(IP_OBS))
	        DO I=OMIT_ST(J)+1,OMIT_END(J)-1
	          MASK(I)=0.0D0
	        END DO
	      END DO
	      WRITE(6,*)'Aq'
!
	      DO J=1,NC_OMIT
	        C_OMIT_ST(J)=GET_INDX_SP(OLAM_ST(J),CD(IP_OBS)%XVEC,NPTS(IP_OBS))
	        C_OMIT_END(J)=GET_INDX_SP(OLAM_END(J),CD(IP_OBS)%XVEC,NPTS(IP_OBS))
	        DO I=C_OMIT_ST(J)+1,C_OMIT_END(J)-1
	          MASK(I)=0.0D0
	        END DO
	      END DO
!
	      SUM_OSQ=0.0D0;    SUM_LOSQ=0.0D0;   SUM_LSQ_OSQ=0.0D0
	      SUM_MO=0.0D0;     SUM_LMO=0.0D0;    NPIX=0.0D0
	      DO I=IST,IEND
	        NPIX=NPIX+NINT(MASK(I)*1.0D0)
	        LAM=CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST)
	        O_DATA=CD(IP_OBS)%DATA(I)*MASK(I)*C_MASK(I)
!
	        SUM_OSQ=SUM_OSQ+O_DATA*O_DATA
	        SUM_LOSQ=SUM_LOSQ+LAM*O_DATA*O_DATA
	        SUM_LSQ_OSQ=SUM_LSQ_OSQ+LAM*LAM*O_DATA*O_DATA
!
	        SUM_MO=SUM_MO+O_DATA*MOD_DATA(I)
	        SUM_LMO=SUM_LMO+LAM*O_DATA*MOD_DATA(I)
	      END DO
	      DET=SUM_OSQ*SUM_LSQ_OSQ-SUM_LOSQ*SUM_LOSQ
	      DA=SUM_MO*SUM_LSQ_OSQ-SUM_LMO*SUM_LOSQ
	      DB=SUM_OSQ*SUM_LMO-SUM_LOSQ*SUM_MO
	      A=DA/DET; B=DB/DET
!
! Compute CHI^2
!
	      RAW_CHISQ=0.0D0; CHISQ=0.0D0
	      DO I=IST,IEND
	        T1=CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST)
	        RAW_CHISQ=RAW_CHISQ+MASK(I)*(MOD_DATA(I)-OBS_DATA(I))**2/MOD_DATA(I)
	        OBS_DATA(I)=(A+B*T1)*OBS_DATA(I)
	        CHISQ=CHISQ+MASK(I)*(MOD_DATA(I)-OBS_DATA(I))**2/MOD_DATA(I)
	      END DO
	      WRITE(6,*)'A'
!
! Simple trapzoidal rule integration.
! Compute model continuum level assuming it is defined close to the line bounds.
!
	     T1=0.0D0; T2=0.0D0
	     DO I=IST,IST+2
	       T1=T1+MOD_DATA(I)
	     END DO
	     DO I=IEND-2,IEND
	       T2=T2+MOD_DATA(I)
	     END DO
	     T1=T1/3; T2=T2/3
	     SLOPE=(T2-T1)/(CD(IP_OBS)%XVEC(IEND-1)-CD(IP_OBS)%XVEC(IST+1))
	     INTER=T1
	    WRITE(6,*)'B'
!
! Since the obersevations have been "normalized" we assume that the
! continuum is the same for the observations.
!
	     EW_OBS=0.0D0; EW_MOD=0.0D0
	     DO I=IST,IEND
	       CONT=INTER+SLOPE*(CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST+1))
	       T1=CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST)
	       T2=0.5D0*(CD(IP_OBS)%XVEC(MAX(IST,I-1))-CD(IP_OBS)%XVEC(MIN(I+1,IEND)))
	       EW_OBS=EW_OBS+(1.0D0-OBS_DATA(I)/CONT)*ABS(T2)
	       EW_MOD=EW_MOD+(1.0D0-MOD_DATA(I)/CONT)*ABS(T2)
	     END DO
	    WRITE(6,*)'C'
!
! Simple trapzoidal rule integration. These EWs exclude the omitted regions.
!
	     EW_OBS_OMIT=0.0D0; EW_MOD_OMIT=0.0D0
	     DO I=IST,IEND
	       CONT=INTER+SLOPE*(CD(IP_OBS)%XVEC(I)-CD(IP_OBS)%XVEC(IST+1))
	       IF(MASK(I) .GT. 0.01)THEN
	         J=MAX(I-1,IST)
	         IF(MASK(J) .LT. 0.01)J=I
	         K=MIN(I+1,IEND)
	        IF(MASK(K) .LT. 0.01)K=I
	        T1=0.5D0*(CD(IP_OBS)%XVEC(J)-CD(IP_OBS)%XVEC(K))
	        EW_OBS_OMIT=EW_OBS_OMIT+(1.0D0-OBS_DATA(I)/CONT)*ABS(T1)
	        EW_MOD_OMIT=EW_MOD_OMIT+(1.0D0-MOD_DATA(I)/CONT)*ABS(T1)
	      END IF
	    END DO
!
! We ignore the clipped regions for computing the mean wavelength.
!
	      MEAN=0.0D0; LAM_CENT=0.0D0	
	      DO I=IST,IEND-1
	        X1=CD(IP_OBS)%XVEC(I); X2=CD(IP_OBS)%XVEC(I+1)
	        D1=MOD_DATA(I)-1.0D0
	        D2=MOD_DATA(I+1)-1.0D0
	        MEAN=MEAN+(D1+D2)*(X2-X1)
	        LAM_CENT=LAM_CENT+(X1*D1+X2*D2)*(X2-X1)
	     END DO
	     LAM_CENT=LAM_CENT/MEAN
!
	     RED_CHISQ=CHISQ*1.0D+04/(NPIX-2)
             T2=EW_MOD; T3=LAM_CENT; T4=50.0   !Rough FWHN in km/s to get closest line
             CALL GET_LINE_ID_PG(TRANS_NAME,T1,T2,T3,T4)
             WRITE(LU_OUT,'(I4,F12.4,6F12.3,6X,A)')IP_MOD,LAM_CENT,EW_MOD,EW_OBS,
	1                          EW_MOD_OMIT,EW_OBS_OMIT,RED_CHISQ,T1,TRIM(TRANS_NAME)
	     FLUSH(LU_OUT)
!
	    END IF		!IP NE IP_OPS
	  END DO		!Over plots
	  WRITE(LU_OUT,*)' '; FLUSH(LU_OUT)
	END DO			!Over lines
1000	CONTINUE
	CLOSE(LU_LIMS)
!
	RETURN
!
	CONTAINS
!
	SUBROUTINE WRITE_BALMER_HEADER(LU)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER LU,IP
!
 	WRITE(LU,'(A)')'! '
 	WRITE(LU,'(A)')'! Reduced Chi^2 value assumes a SN to 100'
 	WRITE(LU,'(A)')'! Reduced Chi^2 at another SN =  Chi^2 * (SN/100)^2'
 	WRITE(LU,'(A)')'! Lam(min) refers to the wavelength with the minimum data values in the defined band'
	WRITE(LU,'(A)')'! '
	DO IP=1,NPLTS
	  WRITE(LU,'(A,I3,5X,A,2X,A)')'! Plot #:',IP,'Plot title:',TRIM(CD(IP)%CURVE_ID)
        END DO
        WRITE(LU,'(A)')'!'
	WRITE(LU,'(A,1X,A,9X,A,8(5X,A))')'!','IP','Lam','EW(mod)','EW(obs)',
	1	         'EW(mod)','EW(obs)','Chi^2','Lam(ID)','Transition'
 	WRITE(LU,'(A,7X)')'!',' Omit window pairs'
 	WRITE(LU,'(A)')'! '
	FLUSH(LU)
!	
	END SUBROUTINE WRITE_BALMER_HEADER
	END SUBROUTINE DO_FILE_BALMER_V1
