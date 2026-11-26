	MODULE VEL_LAW_PARAMS
	USE SET_KIND_MODULE
	INTEGER VEL_TYPE
	INTEGER TRANS_I
	REAL(KIND=LDP) R_TRANS,V_TRANS
	REAL(KIND=LDP) dVdR_TRANS
	REAL(KIND=LDP) VINF,BETA
	REAL(KIND=LDP) BETA2
	REAL(KIND=LDP) VEXT
	REAL(KIND=LDP) RP2
	REAL(KIND=LDP) SCALE_HEIGHT
	REAL(KIND=LDP) RO
	REAL(KIND=LDP) dVdR
	REAL(KIND=LDP) ALPHA
	END MODULE VEL_LAW_PARAMS
!
CONTAINS
!
! Program to modify RVSIG_COL. Various options are available.
! Ideal for revising grid etc.
!
	PROGRAM REVISE_RVSIG
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	USE VEL_LAW_PARAMS
	IMPLICIT NONE
!
! Altered 05-Oct-2026: Reject velocity-coordinate revisions for non-monotonic V.
! Altered 05-Jun-2025: Increased precision of V output.
! Altered 27-Apr-2025: Fixed output bug with CUR/mv option.
! Altered 08-Jan-2023: Option 6 no longer askes for VINF and BETA.
! Altered 06-Jul-2021: Added SSCLR option (transferred from OSIRIS: 21-Jul-2021).
! Altered 07-Dec-2020: Updated SCLR from IBIS version -- to used old V(r) law'
! Altered 08-May-2020: Can use log axis for R or V with CUR option.
! Altered 01-Feb-2020: Added 'm' optoion to CUR option.
! Altered 23-Jan-2020: Some cleaning done.
! Altered 02-Jan-2020: Improved CUR option (R is normalized), added SIG option.
! Altered 19-Aug-2019: Added VEL_LAW_PARAMS, and the two subroutines.
!                        Done to make velocity law calculation for SCLR
!                        and MDOT full consistent.
! Altered 10-Jun-2019: Changed to allow reding of density and clumping fator (necessary, for variable MDOT,
!                        and point source models.
! Altered 05-Mar-2018: Added option to make a plane-parallel models out of an RVSIG_COL
!                        file that was used for a spherical model. Option is 'STOP'
!                        i.e., Spherical To Parallel.
! Altered 06-Nov-2017: Made compatible with version of osiris.
!                        TAU option cleaned, options added.
! Altered 24-Oct-2013: Added VEL_TYPE=3 so that can do a velocity law with 2 components.
! Altered 14-Mar-2011: Improved header output to RVSIG_COL.
! Altered 15-Oct-2010: Fixed bug with SIGMA computation for the extra
!                        points added with the EXTR option.
! Altered 27-Aug-2007: Revised file read so as all ! (1st character)
! comments are ignored.
!
	INTEGER, PARAMETER :: NMAX=5000
	INTEGER, PARAMETER :: IONE=1
!
	REAL(KIND=LDP) OLD_R(NMAX)
	REAL(KIND=LDP) OLD_V(NMAX)
	REAL(KIND=LDP) OLD_SIGMA(NMAX)
	REAL(KIND=LDP) OLD_DENSITY(NMAX)
	REAL(KIND=LDP) OLD_CLUMP_FAC(NMAX)
	REAL(KIND=LDP) OLD_LOG_R(NMAX)
!
	REAL(KIND=LDP) R_PLT(NMAX)
	REAL(KIND=LDP) V_PLT(NMAX)
!
	REAL(KIND=LDP) R(NMAX)
	REAL(KIND=LDP) V(NMAX)
	REAL(KIND=LDP) SIGMA(NMAX)
	REAL(KIND=LDP) DENSITY(NMAX)
	REAL(KIND=LDP) CLUMP_FAC(NMAX)
	REAL(KIND=LDP) LOG_R(NMAX)
!
	REAL(KIND=LDP) RTMP(NMAX)
	REAL(KIND=LDP) OLD_TAU(NMAX)
	REAL(KIND=LDP) TAU_SAV(NMAX)
	REAL(KIND=LDP) TAU(NMAX)
	REAL(KIND=LDP) CHI_ROSS(NMAX)
!
	REAL(KIND=LDP) TMP_R(NMAX)
	REAL(KIND=LDP) X1(NMAX)
	REAL(KIND=LDP) X2(NMAX)	
	REAL(KIND=LDP) WRK_VEC(NMAX)	
!
	REAL(KIND=LDP), ALLOCATABLE :: COEF(:,:)
!
	INTEGER ND_OLD
	INTEGER ND
!
	REAL(KIND=LDP) RX
	REAL(KIND=LDP) NEW_RSTAR
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) FAC
	REAL(KIND=LDP) V_MAX
	REAL(KIND=LDP) V_MIN
!
	REAL*4 XVAL,YVAL,SYMB_EXP_FAC
	REAL(KIND=LDP) MDOT
	REAL(KIND=LDP) OLD_MDOT
	REAL(KIND=LDP) LSTAR
	REAL(KIND=LDP) OLD_LSTAR
	REAL(KIND=LDP) TOP,BOT
	REAL(KIND=LDP) dTOPdR,dBOTdR
!
	REAL(KIND=LDP) V_CON,V_RAT_MAX,IB_RAT,OB_RAT,DTAU2_ON_DTAU1
	INTEGER N_IB_INS,N_OB_INS
!
	INTEGER NN
	INTEGER NX
	INTEGER N_ADD
	INTEGER NX_IN
	INTEGER NX_OUT
	INTEGER NINS
!
	INTEGER I,J,K
	INTEGER I_ST,I_END
	INTEGER N_HEAD
	INTEGER IOS
	INTEGER PGCURS
	INTEGER CURSERR
!
	LOGICAL ROUND_ERROR
	LOGICAL RD_MEANOPAC
	LOGICAL REPLOT
	LOGICAL V_MONO_DEC
	LOGICAL RD_DENSITY
	LOGICAL LOG_X,LOG_Y
	LOGICAL EQUAL; EXTERNAL EQUAL
!
	INTEGER GET_INDX_DP
        CHARACTER*30 UC
        EXTERNAL UC,GET_INDX_DP
!
	CHARACTER(LEN=1) CURSVAL
	CHARACTER(LEN=10) OPTION
	CHARACTER(LEN=80) OLD_RVSIG_FILE
	CHARACTER(LEN=80) NEW_RVSIG_FILE
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=80) OLD_HEADER(30)
	CHARACTER(LEN=20) XLAB,YLAB	
!
	OLD_RVSIG_FILE='RVSIG_COL_OLD'
	CALL GEN_IN(OLD_RVSIG_FILE,'File containing old R, V and sigma values')
	RD_DENSITY=.FALSE.
	CALL GEN_IN(RD_DENSITY,'Does the file contain the density and clumping factors')
	OPEN(UNIT=10,FILE=OLD_RVSIG_FILE,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Old RVSIG_COL file: ',TRIM(OLD_RVSIG_FILE),' does not exist'
	    STOP
	  END IF
	  STRING=' '
	  N_HEAD=0
	  DO WHILE (INDEX(STRING,'!Number of depth points') .EQ. 0)
	    READ(10,'(A)')STRING
	    N_HEAD=N_HEAD+1
	    N_HEAD=MIN(30,N_HEAD)
	    OLD_HEADER(N_HEAD)=STRING
	  END DO
	  N_HEAD=N_HEAD-1
	  READ(STRING,*)ND_OLD
	  STRING=' '
	  DO WHILE (STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(10,'(A)')STRING
	  END DO
	  IF(RD_DENSITY)THEN
	    READ(STRING,*)OLD_R(1),OLD_V(1),OLD_SIGMA(1),OLD_DENSITY(1),OLD_CLUMP_FAC(1)
	    DO I=2,ND_OLD
	      READ(10,*)OLD_R(I),OLD_V(I),OLD_SIGMA(I),OLD_DENSITY(I),OLD_CLUMP_FAC(I)
	    END DO
	  ELSE
	    READ(STRING,*)OLD_R(1),OLD_V(1),OLD_SIGMA(1)
	    DO I=2,ND_OLD
	      READ(10,*)OLD_R(I),OLD_V(I),OLD_SIGMA(I)
	    END DO
	  END IF
	CLOSE(UNIT=10)
!
	OPTION='NEW_ND'
	WRITE(6,'(A)')
	WRITE(6,'(A)')'Current options are:'
	WRITE(6,'(A)')'         ADDR: explicitly add extra grid points'
	WRITE(6,'(A)')'          CUR: Revise V (and/or add extra R values) using a cursor.'
	WRITE(6,'(A)')'         EXTR: extend grid to larger radii'
	WRITE(6,'(A)')'         FGOB: insert N points at outer boundary (to make finer grid)'
	WRITE(6,'(A)')'           FG: insert N points ito make a fine grid (uniform spacing)'
	WRITE(6,'(A)')'          FGT: insert N points ito make a fine grid (uniform spacing in TAU)'
	WRITE(6,'(A)')'         MDOT: change the mass-loss rate or velocity law'
	WRITE(6,'(A)')'         NEWG: revise grid between two velocities or depth indicies'
	WRITE(6,'(A)')'       NEW_ND: revise number of depth points by simple scaling'
	WRITE(6,'(A)')'         NINS: insert points in specified region'
	WRITE(6,'(A)')'         PLOT: plot V and SIGMA from old RVSIG file'
	WRITE(6,'(A)')'          RDR: read in new R grid (DC format)'
	WRITE(6,'(A)')'         SCLR: scale radius of star to new value'
	WRITE(6,'(A)')'        SSCLR: simple scale radius of star to new value (all radii changed).'
	WRITE(6,'(A)')'         SCLV: scale velocity law to new value'
        WRITE(6,'(A)')'         SIGU: update sigma using V read from RVSIG file'
	WRITE(6,'(A)')'          SPP: convert plane-parallel model to a spherical -paralell model'
	WRITE(6,'(A)')'         STOP: convert spherical model to a plane-paralell model'
	WRITE(6,'(A)')'         SUBO: subtract constant offset from R grid (for a plane-paralell model)'
	WRITE(6,'(A)')'          TAU: Create entrie grid equally spaced in log[TAU] with constraints on dV'
	WRITE(6,'(A)')
	WRITE(6,'(A)')'R_GRID_CHK (if MEANOPAC readd) will contain dTAU and dV information.'
	WRITE(6,'(A)')

	CALL GEN_IN(OPTION,'Enter option for revised RVSIG file')
	OPTION=UC(TRIM(OPTION))
!
! These revisions locate a connection or refinement interval by V, so a
! repeated velocity selects an arbitrary radius branch.  Keep the
! radius-coordinate edits available, but reject velocity-coordinate
! edits on fields that are not positive and strictly decreasing.
!
	V_MONO_DEC=.TRUE.
	DO I=1,ND_OLD-1
	  IF(OLD_V(I) .LE. 0.0_LDP .OR.
	1     OLD_V(I+1) .GE. OLD_V(I))V_MONO_DEC=.FALSE.
	END DO
	IF(.NOT. V_MONO_DEC .AND.
	1   (OPTION .EQ. 'TAU' .OR. OPTION .EQ. 'NEWG' .OR.
	2    OPTION .EQ. 'SCLR' .OR. OPTION .EQ. 'MDOT' .OR.
	3    OPTION .EQ. 'STOP'))THEN
	  WRITE(6,*)'Error in REVISE_RVSIG:'
	  WRITE(6,*)'Option ',TRIM(OPTION),' requires a positive,'
	  WRITE(6,*)'strictly decreasing velocity field because it uses'
	  WRITE(6,*)'V as an independent coordinate.'
	  STOP
	END IF
!
! Read in optical depth scale. Needed for RTAU and TAU options. Also used
! for checking purposes (if available) for some other options.
!
! We use TAU_SAV for dTAU, and is used to increase the precision of the TAU
! scale (as insufficient digits may be print out).
!
	ROUND_ERROR=.FALSE.
	OPEN(UNIT=20,FILE='MEANOPAC',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(20,'(A)')STRING
	    DO I=1,ND_OLD
	     READ(20,*)RTMP(I),J,OLD_TAU(I),TAU_SAV(I),T1,CHI_ROSS(I)
	      J=MAX(I,2)
	      T1=OLD_R(J-1)-OLD_R(J)
	      IF( ABS(RTMP(I)-R(I))/T1 .GT. 2.0E-03_LDP .AND. .NOT. ROUND_ERROR)THEN
	        WRITE(6,*)' '
	        WRITE(6,*)'Possible eror with MEANOPAC -- inconsistent R grid'
	        WRITE(6,*)'Error could simply be a lack of sig. digits in MEANOPAC'
	        WRITE(6,*)' RMO(I)=',RTMP(I)
	        WRITE(6,*)'   R(I)=',OLD_R(I)
	        WRITE(6,*)' R(I+1)=',OLD_R(I+1)
	        ROUND_ERROR=.TRUE.
	        CALL GEN_IN(ROUND_ERROR,'Continue as only rounding error?')
	        IF(.NOT. ROUND_ERROR)STOP
	      END IF
	    END DO
	    DO I=8,1,-1
              OLD_TAU(I)=OLD_TAU(I+1)-TAU_SAV(I)
            END DO
	    CLOSE(UNIT=20)
	    RD_MEANOPAC=.TRUE.
	    IF(ROUND_ERROR)RTMP(1:ND_OLD)=OLD_R(1:ND_OLD)
	    WRITE(6,*)'Successfully read MEANOPAC'
	  ELSE
	    RD_MEANOPAC=.FALSE.
	    WRITE(6,*)'MEANOPAC has not been read'
	    IF(OPTION .EQ. 'SPP' .OR. OPTION .EQ. 'TAU')THEN
	      WRITE(6,*)'MENAOPAC  is equired for the TAU option'
	      STOP
	    END IF
	  END IF
!
	IF(OPTION .EQ. 'SPP')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option takes a plane-parallel model and oututs a spherical model'
	  WRITE(6,'(A)')'The MEANOPAC from the plane-parallel model is required'
	  WRITE(6,'(A)')'The VADAT file is also required'
	  WRITE(6,'(A)')'Option not yet implemented'
	  STOP
	END IF
!
	IF(OPTION .EQ. 'TAU')THEN
!
	  ND=ND_OLD
	  V_TRANS=30.0_LDP; V_RAT_MAX=1.5
	  IB_RAT=2.0_LDP; OB_RAT=4.0_LDP; DTAU2_ON_DTAU1=100.0_LDP
	  N_IB_INS=2; N_OB_INS=1
!
	  CALL GEN_IN(ND,'Number of depth points')
	  CALL GEN_IN(V_CON,'Connection velocity in km/s')
	  CALL GEN_IN(IB_RAT,'Ratio of spacing at the inner boundary')
	  CALL GEN_IN(OB_RAT,'Ratio of spacing at the outer boundary')
	  CALL GEN_IN(N_IB_INS,'Number of extra points to insert at the inner boudary')
	  CALL GEN_IN(N_OB_INS,'Number of extra points to insert at the outer boudary')
	  CALL GEN_IN(DTAU2_ON_DTAU1,'DTAU(2)/DTAU(1)')
!
	  CALL ADJUST_ATM_R_GRID(R,OLD_R,OLD_V,OLD_TAU,
	1         V_CON,V_RAT_MAX,IB_RAT,OB_RAT,
	1         DTAU2_ON_DTAU1,N_IB_INS,N_OB_INS,ND,ND_OLD)
!
! Now compute the revised SIGMA. V has already been computed.
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
	  WRITE(6,*)'Determined R and V on the new grid'
!
	ELSE IF(OPTION .EQ. 'RDR')THEN
	  OLD_RVSIG_FILE='RDINR'
	  CALL GEN_IN(OLD_RVSIG_FILE,'File with new R values -- DC format')
	  OPEN(UNIT=9,FILE=OLD_RVSIG_FILE,STATUS='OLD',ACTION='READ')
	  DO I=1,3
	    READ(9,'(A)')STRING
	  END DO
	  READ(9,'(A)')STRING
	  READ(STRING,*)T1,T1,J,ND
	  READ(9,'(A)')STRING                 !Final blank line
!
	  DO I=1,ND
	    READ(9,'(A)')STRING
	    READ(STRING,*)R(I)
	    DO WHILE(STRING .NE. ' ')
	      READ(9,'(A)',END=100)STRING
	    END DO
	  END DO
100	  CONTINUE
	  CLOSE(UNIT=9)
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
!
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)

	ELSE IF(OPTION .EQ. 'NEW_ND')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows a new R grid to be output'
	  WRITE(6,'(A)')'Grid spacing is similar to input grid.'
	  WRITE(6,'(A)')' '
	  ND=70
	  CALL GEN_IN(ND,'Number of depth points')
	  DO I=1,ND_OLD
	      X1(I)=I
	  END DO
	  T1=DFLOAT(ND_OLD-1)/DFLOAT(ND-1)
	  DO I=1,ND
	    X2(I)=(I-1)*T1+1
	  END DO
	  X2(1)=X1(1)
	  X2(ND)=X1(ND_OLD)
	
	  CALL MON_INTERP(R,ND,IONE,X2,ND,OLD_R,ND_OLD,X1,ND_OLD)
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
!
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
!
	ELSE IF(OPTION .EQ. 'FIXV')THEN
	  ND=ND_OLD
	  R(1:ND)=OLD_R(1:ND)
	  V(1:ND)=OLD_V(1:ND)
	  SIGMA(1:ND)=OLD_SIGMA(1:ND)
	  T1=(SIGMA(3)+1.0D0)*V(3)/R(3)
	  V(3)=V(3)+T1*(R(2)-R(3))
	  T1=(SIGMA(2)+1.0D0)*V(2)/R(2)
	  V(1)=V(2)+T1*(R(1)-R(2))
!
	ELSE IF(OPTION .EQ. 'NINS')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows grid to be increaed by iserting points in old grid'
	  WRITE(6,'(A)')'The old grid points are retained'
	  WRITE(6,'(A)')' '
	  NINS=1
	  CALL GEN_IN(NINS,'Number of points to insert')
	  CALL GEN_IN(I_ST,'Miminum index for grid reginement')
	  CALL GEN_IN(I_END,'Maximumindex for grid reginement')
!
	  R(1:I_ST)=OLD_R(1:I_ST); J=I_ST
	  DO I=I_ST,I_END-1
	    T1=(OLD_R(I+1)-OLD_R(I))/(NINS+1)
	    DO K=1,NINS
	      J=J+1
	      R(J)=R(J-1)+T1
	    END DO
	    J=J+1
	    R(J)=OLD_R(I+1)
	  END DO
	  DO I=I_END+1,ND_OLD
	    J=J+1
	    R(J)=OLD_R(I)
	  END DO
	  ND=J
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
!
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'SCLV')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option simply scales V'
	  WRITE(6,'(A)')'Sigma will not be changed by this operation.'
	  WRITE(6,'(A)')' '
	  T1=1.1_LDP
	  ND=ND_OLD
	  CALL GEN_IN(T1,'Factor to sclae velocity')
	  T2=0.0_LDP
	  CALL GEN_IN(T2,'Velocity below which scaling is skipped')
	  DO I=1,ND
	    IF(OLD_V(I) .GT. T2)THEN
	      V(I)=T1*OLD_V(I)
	    ELSE
	      V(I)=OLD_V(I)
	    END IF
	  END DO
	  R(1:ND)=OLD_R(1:ND)
	  SIGMA(1:ND)=OLD_SIGMA(1:ND)
!
	ELSE IF(OPTION .EQ. 'FG' .OR. OPTION .EQ. 'FGT')THEN
          I_ST=1; I_END=ND
          CALL GEN_IN(I_ST,'Start index for fine grid')
          CALL GEN_IN(I_END,'End index for fine grid')
          WRITE(6,*)'Number of points in requeted interval is',I_END-I_ST-1
          NX=I_END-I_ST
          CALL GEN_IN(NX,'New number of grid points for this interval')
          ND=I_ST+NX+(ND_OLD-I_END)+1
!
	  IF(OPTION .EQ. 'FG')THEN
            DO I=1,I_ST
              X2(I)=I
            END DO
            T1=(I_END-I_ST)/(NX+1.0_LDP)
            DO I=1,NX
              X2(I_ST+I)=I_ST+I*T1
            END DO
            DO I=I_END,ND_OLD
              X2(I_ST+NX+I+1-I_END)=I
            END DO
            DO I=1,ND_OLD
              X1(I)=I
            END DO
	  ELSE
            DO I=1,I_ST
              X2(I)=OLD_TAU(I)
            END DO
            T1=EXP(LOG(OLD_TAU(I_END)/OLD_TAU(I_ST))/(NX+1.0_LDP))
            DO I=1,NX
              X2(I_ST+I)=X2(I_ST+I-1)*T1
	      WRITE(6,*)X2(I_ST+I),X2(I_ST+I-1)
            END DO
            DO I=I_END,ND_OLD
              X2(I_ST+NX+I+1-I_END)=OLD_TAU(I)
	      WRITE(6,*)X2(I_ST+NX+I+1-I_END)
            END DO
            DO I=1,ND_OLD
              X1(I)=OLD_TAU(I)
            END DO
          END IF
	  CALL MON_INTERP(R,ND,IONE,X2,ND,OLD_R,ND_OLD,X1,ND_OLD)
	  DO I=1,ND_OLD
	     WRITE(33,*)I,OLD_R(I),X1(I),OLD_TAU(I)
	   END DO
	  DO I=1,ND
	     WRITE(34,*)I,R(I),X2(I)
	   END DO
	  WRITE(6,*)'Done interpolation'
!
! Now compute the revised V and SIGMA.
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'ADDR')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows extra grid points to be added'
	  WRITE(6,'(A)')'The grid points are specfied by the user as requested'
	  WRITE(6,'(A)')'Please insert the grid points largest to smallest'
	  N_ADD=1
	  CALL GEN_IN(N_ADD,'Number of grid points to be added')
	  DO I=1,N_ADD
	    CALL GEN_IN(X1(I),'New grid point:')
	  END DO
	  DO I=1,N_ADD-1
	    IF(X1(I) .LE. X1(I+1))THEN
	      WRITE(6,*)'Error: new grid points are invalid -- equal or wrong order'
	      STOP
	    END IF
	  END DO
	  IF(X1(1) .GE. OLD_R(1))THEN
	    WRITE(6,*)'Error: can''t insert grid points outside defined grid'
	    WRITE(6,*)'OLD_R(1)=',OLD_R(1)
	    WRITE(6,*)'R_INS=',X1(1)
	    STOP
	  END IF
	  IF(X1(N_ADD) .LE. OLD_R(ND_OLD))THEN
	    WRITE(6,*)'Error: can''t insert grid points outside defined grid'
	    WRITE(6,*)'OLD_R(ND_OLD)=',OLD_R(ND_OLD)
	    WRITE(6,*)'R_INS=',X1(N_ADD)
	    STOP
	  END IF
!
	  ND=N_ADD+ND_OLD
	  J=1
	  I=2
	  K=1
	  R(J)=OLD_R(1)
	  DO WHILE (K .LE. N_ADD)
	    IF(X1(K) .LE. R(J) .AND. X1(K) .GT. OLD_R(I))THEN
	      J=J+1
	      R(J)=X1(K)
	      K=K+1
	    ELSE
	      J=J+1
	      R(J)=OLD_R(I)
	      I=I+1
	    END IF
	  END DO
	  DO K=J+1,ND
	    R(K)=OLD_R(I)
	    I=I+1
	  END DO
!
! Check ordering etc.
!
	DO I=1,ND-1
	  IF(R(I) .LE. R(I+1))THEN
	    WRITE(6,*)'Error: R values not monotonic'
	    WRITE(6,*)'I=',I
	    WRITE(6,*)'R(I)=',R(I)
	    WRITE(6,*)'R(I+1)=',R(I+1)
	    STOP
	  END IF
	END DO
!
! Now compute the revised V and SIGMA.
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'NEWG')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows the grid to be redefined over a specified interval'
	  WRITE(6,'(A)')'Originally defined for small ranges'
	  WRITE(6,'(A)')'Data points in the requested interval are equally spaced on log r^2. V'
	  WRITE(6,'(A)')' '
!
	  CALL GEN_IN(V_MAX,'Maximum velocity for grid refinement (-ve for depth index)')
	  CALL GEN_IN(V_MIN,'Initial velocity for grid refinement (-ve for depth index)')
	  IF(V_MAX*V_MIN .LT. 0.0_LDP)THEN
	    WRITE(6,*)'Error: either velocities of depth index but not mixture'
	    STOP
	  END IF
!
	  IF(V_MAX .LT. 0.0_LDP)THEN
	    V_MAX=ABS(V_MAX); V_MIN=ABS(V_MIN)
	    I_ST=NINT(MIN(V_MIN,V_MAX))
	    I_END=NINT(MAX(V_MIN,V_MAX))
	  ELSE
	    IF(V_MAX .GE. OLD_V(1) .OR. V_MAX .LE. OLD_V(ND_OLD))THEN
	      WRITE(6,*)'Error: can''t insert grid points outside velociy grid'
	      STOP
	    END IF
	    IF(V_MIN .GE. OLD_V(1) .OR. V_MIN .LE. OLD_V(ND_OLD))THEN
	      WRITE(6,*)'Error: can''t insert grid points outside velociy grid'
	      STOP
	    END IF
!
! Find interval
!
	    DO I=1,ND_OLD-1
	      IF(V_MAX .GT. OLD_V(I+1))THEN
	        I_ST=I
	        EXIT
	      END IF
	    END DO
	    IF(V_MAX-OLD_V(I+1) .LT. OLD_V(I)-V_MAX)I_ST=I+1
!
	    DO I=1,ND_OLD-1
	      IF(V_MIN .GT. OLD_V(I+1))THEN
	        I_END=I
	        EXIT
	      END IF
	    END DO
	    IF(V_MIN-OLD_V(I_END+1) .LT. OLD_V(I_END)-V_MIN)I_END=I_END+1
	  END IF
!
	  WRITE(6,*)' I_ST=',I_ST,OLD_V(I_ST)
	  WRITE(6,*)'I_END=',I_END,OLD_V(I_END)
	  V_MAX=OLD_V(I_ST); V_MIN=OLD_V(I_END)
	  WRITE(6,*)'Current number of points in interval is',I_END-I_ST-1
	  CALL GEN_IN(N_ADD,'Number of points in interval (exclusive)')
!
	  X1(1:ND_OLD)=OLD_V(1:ND_OLD)*OLD_R(1:ND_OLD)*OLD_R(1:ND_OLD)
	  ND=N_ADD+ND_OLD-(I_END-I_ST-1)
	  X2(1:I_ST)=X1(1:I_ST)
	  T1=EXP(LOG(V_MAX*OLD_R(I_ST)*OLD_R(I_ST)/V_MIN/OLD_R(I_END)/OLD_R(I_END))/(N_ADD+1))
	  DO I=I_ST+1,I_ST+N_ADD
	   X2(I)=X2(I-1)/T1
	  END DO
	  X2(I_ST+N_ADD+1:ND)=X1(I_END:ND_OLD)
!
	  DO I=1,ND; WRITE(6,*)I,X2(I); END DO; FLUSH(UNIT=6)
	  WRITE(6,*)'Call ing MON_INTERP'
	  CALL MON_INTERP(R,ND,IONE,X2,ND,OLD_R,ND_OLD,X1,ND_OLD)
!
! Now compute the revised SIGMA. V has already been computed.
!
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  J=1
	  I=1
	  DO WHILE (I .LE. ND)
	    IF(R(I) .GE. OLD_R(J+1))THEN
	      T1=R(I)-OLD_R(J)
	      V(I)=COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1)))
	      SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	      I=I+1
	    ELSE
	      J=J+1
	    END IF
	  END DO
	  DEALLOCATE (COEF)
	
	ELSE IF(OPTION .EQ. 'FGOB')THEN
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'This option allows the grid to be redefined at the outer boundary'
	  WRITE(6,'(A)')'Adding N points decreases the spacing at the outer boundary by 3^N'
!
	  N_ADD=1
	  CALL GEN_IN(N_ADD,'Number of grid points to be added at outer boundary')
	  ND=ND_OLD+N_ADD
	  R(N_ADD+2:ND)=OLD_R(2:ND_OLD)
	  R(1)=OLD_R(1)
	  DO I=N_ADD,1,-1
	    R(I+1)=R(1)-0.3333_LDP*(R(1)-R(I+2))
	    WRITE(6,*)I,R(1),R(I+1)
	  END DO
!
! Now compute the revised V & SIGMA. The new points lie in the first interval.
!
	  V(1)=OLD_V(1); V(N_ADD+2:ND)=OLD_V(2:ND_OLD)
	  SIGMA(1)=OLD_SIGMA(1); SIGMA(N_ADD+2:ND)=OLD_SIGMA(2:ND_OLD)
	  ALLOCATE (COEF(ND_OLD,4))
	  CALL MON_INT_FUNS_V2(COEF,OLD_V,OLD_R,ND_OLD)
	  DO I=2,N_ADD+1
	    T1=R(I)-OLD_R(1)
	    V(I)=COEF(1,4)+T1*(COEF(1,3)+T1*(COEF(1,2)+T1*COEF(1,1)))
	    SIGMA(I)=COEF(1,3)+T1*(2.0_LDP*COEF(1,2)+3.0_LDP*T1*COEF(1,1))
	    SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	  END DO
	  DEALLOCATE (COEF)
	
	ELSE IF(OPTION .EQ. 'EXTR')THEN
	  FAC=2.0_LDP
	  CALL GEN_IN(FAC,'Factor to extend RMAX by')
          NX_OUT=2
	  VINF=1000.0_LDP
	  CALL GEN_IN(VINF,'Velocity at infinity in km/s')
	  BETA=1.0_LDP
	  CALL GEN_IN(BETA,'Beta for velocity law')
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Outputing old grid near outer boudary for NX estimate.'
	  WRITE(6,*)'NB: The grid ratio will differ for NX+1 values.'
	  WRITE(6,*)' '
	  DO I=1,7
	    WRITE(6,*)I,OLD_R(I)/OLD_R(I+1)
	  END DO
	  NX_OUT=2
	  CALL GEN_IN(NX_OUT,'Number of points used at outer boundary to refine grid')
!
	  RX=OLD_R(1)*(1.0_LDP-(OLD_V(1)/VINF)**(1.0_LDP/BETA))
	  WRITE(6,*)'RX=',RX
!
! Set up a rough grid so we can define a better grid, equally spaced
! in log density.
!
	  TMP_R(1)=FAC*OLD_R(1)
	  NX=20
	  T1=EXP(LOG(FAC)/NX)
	  TMP_R(1)=FAC*OLD_R(1)
	  DO I=2,NX
	    TMP_R(I)=TMP_R(I-1)/T1
	  END DO
	  TMP_R(NX+1)=OLD_R(1)
	  DO I=1,NX+1
	    V(I)=VINF*(1.0_LDP-RX/TMP_R(I))**BETA
	    X1(I)=TMP_R(I)*TMP_R(I)*V(I)
	  END DO
	  NX=NX+1
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Outputing old V.r^2 near outer boudary for GRID estimate.'
	  WRITE(6,*)' '
	  DO I=1,7
	    WRITE(6,*)I,(OLD_V(I)/OLD_V(I+1))*(OLD_R(I)/OLD_R(I+1))**2
	  END DO
	  FAC=SQRT(OLD_R(5)*OLD_R(5)*OLD_V(5)/OLD_R(7)/OLD_R(7)/OLD_V(7))
	  CALL GEN_IN(FAC,'Grid spacing ratio')
	  T1=LOG(X1(1)/X1(NX))
	  NN=NINT(T1/LOG(FAC))
	  T1=EXP(T1/NN)
	  WRITE(6,*)'Number of points to be add is',NN
	  WRITE(6,*)'Grid spacing factor is ',T1
!
	  X2(1)=X1(1)
	  DO I=2,NN
	    X2(I)=X2(I-1)/T1
	  END DO
	  CALL MON_INTERP(R,NN,IONE,X2,NN,TMP_R,NX,X1,NX)
!
	  DO I=NN,2,-1
	    R(I+NX_OUT)=R(I)
	  END DO
	  IF(NX_OUT .EQ. 1)THEN
	    R(2)=R(1)-0.2_LDP*(R(1)-R(3))
	  ELSE IF(NX_OUT .EQ. 2)THEN
	    R(2)=R(1)-0.1_LDP*(R(1)-R(4))
	    R(3)=R(1)-0.4_LDP*(R(1)-R(4))
	  ELSE IF(NX_OUT .EQ. 3)THEN
	    R(2)=R(1)-0.05_LDP*(R(1)-R(5))
	    R(3)=R(1)-0.15_LDP*(R(1)-R(5))
	    R(4)=R(1)-0.40_LDP*(R(1)-R(5))
	  ELSE IF(NX_OUT .EQ. 3)THEN
	    R(2)=R(1)-0.015_LDP*(R(1)-R(6))
	    R(3)=R(1)-0.05_LDP*(R(1)-R(6))
	    R(4)=R(1)-0.15_LDP*(R(1)-R(6))
	    R(5)=R(1)-0.40_LDP*(R(1)-R(6))
	  END IF
!
	  DO I=1,NN+NX_OUT
	    V(I)=VINF*(1.0_LDP-RX/R(I))**BETA
	    SIGMA(I)=BETA*RX/R(I)/(1.0_LDP-RX/R(I))-1.0_LDP
	  END DO
!
! We remove the fine grid in the old model.
!
	  J=NN+NX_OUT+1
	  R(J)=OLD_R(1)
	  V(J)=OLD_V(1)
	  SIGMA(J)=OLD_SIGMA(1)
!
	  DO I=NX_OUT+2,ND_OLD
	    J=NN+I
	    R(J)=OLD_R(I)
	    V(J)=OLD_V(I)
	    SIGMA(J)=OLD_SIGMA(I)
	  END DO
	  ND=NN+ND_OLD
!
	ELSE IF(OPTION .EQ. 'SSCLR')THEN
!
	  ND=ND_OLD
	  NEW_RSTAR=OLD_R(ND_OLD)
	  CALL GEN_IN(NEW_RSTAR,'New radius')
	  R(1:ND)=OLD_R(1:ND)*(NEW_RSTAR/OLD_R(ND_OLD))
	  V(1:ND)=OLD_V(1:ND)
	  SIGMA(1:ND)=OLD_SIGMA(1:ND)
!
	ELSE IF(OPTION .EQ. 'SCLR')THEN
!
	  ND=ND_OLD
	  NEW_RSTAR=OLD_R(ND_OLD)
	  V_TRANS=4.0_LDP
	  CALL GEN_IN(NEW_RSTAR,'New radius')
	  CALL GEN_IN(V_TRANS,'Connection velocity in km/s')
	  CALL GEN_IN(OLD_MDOT,'Old mass-loss rate in Msun/yr')
	  MDOT=OLD_MDOT
!
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'The defaut mass assume Mdot/R^1.5 is to be preserved'
	  WRITE(6,*)DEF_PEN
	  MDOT=OLD_MDOT*(NEW_RSTAR/OLD_R(ND_OLD))**1.5_LDP
	  CALL GEN_IN(MDOT,'New mass-loss rate in Msun/yr')
!
	  CALL GEN_IN(OLD_LSTAR,'Old stellar luminosity in Lsun')
	  T1=OLD_LSTAR*(NEW_RSTAR/OLD_R(ND_OLD))**2
	  WRITE(6,*)RED_PEN
	  WRITE(6,'(A,ES14.4)')'New luminosity if Teff is to be preserved',T1
	  WRITE(6,*)DEF_PEN
!
	  CALL DESCRIBE_VEL_LAWS()
!
	  VEL_TYPE=1
	  CALL GEN_IN(VEL_TYPE,'Velocity law to be used: 1 or 2')
	  IF(VEL_TYPE .NE. 6)THEN
	    VINF=1000.0_LDP
	    CALL GEN_IN(VINF,'Velocity at infinity in km/s')
	    BETA=1.0_LDP
	    CALL GEN_IN(BETA,'Beta for velocity law')
	  END IF
!
! Find conection velocity and index.
!
	  TRANS_I=0
	  DO I=1,ND_OLD
	    IF(V_TRANS .LE. OLD_V(I) .AND. V_TRANS .GE. OLD_V(I+1))THEN
	      TRANS_I=I
	      EXIT
	    END IF
	  END DO
	  IF(TRANS_I .EQ. 0)THEN
	    WRITE(6,'(/,1X,A)')'Error V_TRANS is outside range'
	    WRITE(6,'(1X,3(A,ES15.8,3X))')'V_TRANS=',V_TRANS,'OLD_V(1)=',OLD_V(1),'V(ND_OLD)=',OLD_V(ND_OLD)
	  END IF

	  IF( OLD_V(TRANS_I)-V_TRANS .GT. V_TRANS-OLD_V(TRANS_I+1))TRANS_I=TRANS_I+1
	  V_TRANS=OLD_V(TRANS_I)
	  R(1:ND)=OLD_R(1:ND_OLD)+(NEW_RSTAR-OLD_R(ND_OLD))
	  R(ND)=NEW_RSTAR
!
! Now do the new wind law, keeping the same radius grid.
!
	  IF(VEL_TYPE .EQ. 6)THEN
	    T1=MDOT/OLD_MDOT
	    DO I=1,ND
	      T3=(OLD_R(I)/R(I))**2
	      T2=(OLD_V(I)/V_TRANS)**2
	      V(I)=T3*OLD_V(I)*(T1+(1.0_LDP-T1)*T2/(1.0_LDP+T2))
	    END DO
	    ALLOCATE (COEF(ND,4))
	    CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	    DO I=1,ND
	      SIGMA(I)=COEF(I,3)
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	    END DO
	    DEALLOCATE (COEF)
	  ELSE
!
! In the hydrostatic region, the velocity is simply scaled so as to keep the
! density constant.
!
	  T1=R(ND)/OLD_R(ND_OLD)
	  DO I=TRANS_I,ND
	    V(I)=MDOT*OLD_V(I)/T1/T1/OLD_MDOT
	    SIGMA(I)=(OLD_SIGMA(I)+1.0_LDP)*R(I)/OLD_R(I)-1.0_LDP
	  END DO
	  V_TRANS=V(TRANS_I)
!
!
!
! Now compute the velocity law beyond the sonic point.
!
	  R_TRANS=R(TRANS_I)
	  dVdR_TRANS=(SIGMA(TRANS_I)+1.0_LDP)*V_TRANS/R_TRANS
	  CALL CALCULATE_VEL(R,V,SIGMA,ND)
	  END IF
!
	ELSE IF(OPTION .EQ. 'MDOT')THEN
!
	  ND=ND_OLD
	  OLD_MDOT=0.0_LDP; VINF=0.0_LDP; BETA=1.0_LDP; V_TRANS=10.0_LDP
	  CALL GEN_IN(OLD_MDOT,'Old mass-loss rate in Msun/yr')
	  MDOT=OLD_MDOT
	  CALL GEN_IN(MDOT,'New mass-loss rate in Msun/yr')
!
	  WRITE(6,*)' '
	  WRITE(6,*)' If increasing Mdot, simply enter a number slightly smaller than the sound speed.'
	  WRITE(6,*)' If decreasing Mdot, you may need to enter a smaller number,'
	  WRITE(6,*)' since as you decrease Mdot, the extent of the photosphere increase'
!
	  CALL DESCRIBE_VEL_LAWS()
!
	  VEL_TYPE=2
	  CALL GEN_IN(VEL_TYPE,'Velocity law to be used: 1, 2, 3, 4, 5, or 6')
	  IF(VEL_TYPE .NE. 6)THEN
	    CALL GEN_IN(VINF,'Velocity at infinity in km/s')
	    CALL GEN_IN(BETA,'Beta for velocity law')
	  END IF
	  CALL GEN_IN(V_TRANS,'Connection velocity in km/s')
	  V_TRANS=V_TRANS*OLD_MDOT/MDOT
!
! Find conection velocity and index.
!
	  DO I=1,ND_OLD
	    IF(V_TRANS .LE. OLD_V(I) .AND. V_TRANS .GE. OLD_V(I+1))THEN
	      TRANS_I=I
	      EXIT
	    END IF
	  END DO
	  IF( OLD_V(TRANS_I)-V_TRANS .GT. V_TRANS-OLD_V(TRANS_I+1))TRANS_I=TRANS_I+1
	  R(1:ND)=OLD_R(1:ND_OLD)
!
! Now adjust the velocity law.
!
! In the hydrostatic region, the velocity is simply scaled by the change in
! mass-loss rate. This preserves the density. Only valid if wind does not have
! a significant optical depth.
!
!
	  IF(VEL_TYPE .EQ. 6 .OR. VEL_TYPE .EQ. 7)THEN
	    T1=MDOT/OLD_MDOT
	    DO I=1,ND
	      T2=(OLD_V(I)/V_TRANS)**2
	      V(I)=OLD_V(I)*(T1+(1.0_LDP-T1)*T2/(1.0_LDP+T2))
	    END DO
	    IF(VEL_TYPE .EQ. 7)THEN
	     DO I=1,ND
	        V(I)=V(I)*( 1.3D0-0.3*(1.0-SQRT(R(I)/R(1)))**4)
	      END DO
	    END IF
	    ALLOCATE (COEF(ND,4))
	    CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	    DO I=1,ND
	      SIGMA(I)=COEF(I,3)
	      SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	    END DO
	    DEALLOCATE (COEF)
!
	  ELSE
	    V_TRANS=OLD_V(TRANS_I)
	    DO I=TRANS_I,ND
	      V(I)=MDOT*OLD_V(I)/OLD_MDOT
	      SIGMA(I)=OLD_SIGMA(I)
	    END DO
!
	    R_TRANS=R(TRANS_I)
	    V_TRANS=MDOT*V_TRANS/OLD_MDOT
	    dVdR_TRANS=(SIGMA(TRANS_I)+1.0_LDP)*V_TRANS/R_TRANS
	    CALL CALCULATE_VEL(R,V,SIGMA,ND)
	  END IF
!
	ELSE IF(OPTION .EQ. 'STOP')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)' Create a new RVSIG_COL file for a plane-parallel model from a spherical modl'
	  WRITE(6,*)' '
!
	  ND=ND_OLD
	  CALL GEN_IN(ND,'Increase number depth points (larger?)')
	  CALL GEN_IN(V_TRANS,'Velocity at which to extend plane-parallel atmosphere (km/s)')
	  IF(ND .EQ. ND_OLD)THEN
	    R=OLD_R; V=OLD_V
	  ELSE IF(ND .LT. ND_OLD)THEN
	    J=ND_OLD-ND
	    R(1:ND)=OLD_R(J+1:ND)
	    V(1:ND)=OLD_V(J+1:ND)
	  ELSE
	    J=ND-ND_OLD
	    R(J+1:ND)=OLD_R(1:ND)
	    V(J+1:ND)=OLD_V(1:ND)
	    DO I=J,1,-1
	      V(I)=V(I+1)*1.01_LDP		!Value irrelevant
	    END DO
	  END IF
	  TRANS_I=GET_INDX_DP(V_TRANS,V,ND)
	  SCALE_HEIGHT=(R(TRANS_I)-R(TRANS_I+1))/LOG(V(TRANS_I)/V(TRANS_I+1))
!
	  WRITE(6,*)'     TRANS_I=',TRANS_I
	  WRITE(6,*)'  R(TRANS_I)=',R(TRANS_I)
	  WRITE(6,*)'  V(TRANS_I)=',V(TRANS_I)
	  WRITE(6,*)'SCALE_HEIGHT=',SCALE_HEIGHT
!
	  T1=3.0_LDP; CALL GEN_IN(T1,'Number of depth points per scale height')
	  IF(T1 .GT. 1)T1=1.0_LDP/T1
	  DO I=TRANS_I-1,1,-1
	    R(I)=R(I+1)+T1*SCALE_HEIGHT
	    V(I)=V(TRANS_I)*EXP( (R(I)-R(TRANS_I))/SCALE_HEIGHT)
	  END DO
	  V(1:ND)=1.0E-15_LDP*V(1:ND)
	  WRITE(6,*)'Decrease the mass-loss rate in the VADAT file by 10^15'
!
	  ALLOCATE (COEF(ND,4))
	  CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	  DO I=1,ND
	    SIGMA(I)=COEF(I,3)
	    SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'SUBO')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)' Create a new RVSIG_COL file for a plane-parallel model'
	  WRITE(6,*)' Subtract a constant offset off - designed to facilitate WD modeling.'
	  WRITE(6,*)' Gives a more reasonable scale height compard to R'
	  WRITE(6,*)' '
!
	  CALL GEN_IN(T1,'Radius value to subtract from file')
!
! The scaling of V preserves the density for the aopted mass-loss rate.
!
	  ND=ND_OLD
	  R=OLD_R; V=OLD_V
	  DO I=1,ND
	    R(I)=R(I)-T1
	  END DO
	  V(1:ND)=V(1:ND)*(OLD_R(1:ND)/R(1:ND))**2
!
	  ALLOCATE (COEF(ND,4))
	  CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	  DO I=1,ND
	    SIGMA(I)=COEF(I,3)
	    SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'SIGU')THEN
!
! Update SIGMA
!
	  ND=ND_OLD
	  ALLOCATE (COEF(ND,4))
	  R(1:ND)=OLD_R(1:ND)
	  V(1:ND)=OLD_V(1:ND)
	  CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	  DO I=1,ND
	    SIGMA(I)=COEF(I,3)
	    SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'CUR')THEN
!
          WRITE(6,'(A)')BLUE_PEN
          WRITE(6,*)' '
          WRITE(6,*)' First create desired plot, amd use M to mark data point'
          WRITE(6,*)' Do NOT change type of axes (i.e. take log etc)'
          WRITE(6,*)' Then exit plot package -- cursor input done outside plot routine'
          WRITE(6,'(A)')DEF_PEN
!
	  LOG_X=.FALSE.; LOG_Y=.FALSE.
	  CALL GEN_IN(LOG_X,'Are we using LOG axes for R')	
	  CALL GEN_IN(LOG_Y,'Are we using LOG axes for V')	
	  ND=ND_OLD;    T1=OLD_R(ND)
	  IF(LOG_X)THEN
	    R_PLT(1:ND)=LOG10(OLD_R(1:ND)/T1)
	    XLAB='Log R/R\d*\u'
	  ELSE
	    R_PLT(1:ND)=OLD_R(1:ND)/T1
	    XLAB='R/R\d*\u'
	  END IF
	  IF(LOG_Y)THEN
	    V_PLT(1:ND)=LOG10(OLD_V(1:ND))
	    YLAB='Log V(km/s)'
	  ELSE
	    V_PLT(1:ND)=OLD_V(1:ND)
	    YLAB='V(km/s)'
	  END IF
	  CALL DP_CURVE(ND_OLD,R_PLT,V_PLT)
 	  CALL GRAMON_PGPLOT(XLAB,YLAB,' ',' ')
!
	  R(1:ND)=R_PLT(1:ND); V(1:ND)=V_PLT(1:ND)
	  REPLOT=.TRUE.
!
	  DO WHILE(1 .EQ. 1)				!Multiple plotting
!
	    WRITE(6,'(A)')BLUE_PEN
	    WRITE(6,'(A)')' '
	    WRITE(6,*)'Click on top bar of plot window to activate cursor'
	    WRITE(6,'(A)')' '
	    WRITE(6,*)'Use ''r'' to replace a data point'
	    WRITE(6,*)'Use ''a'' to add a data point'
	    WRITE(6,*)'Use ''d'' to delete a data point'
	    WRITE(6,*)'Use ''m'' move next, and points further out/inward by a fixed amount'
	    WRITE(6,*)'Use ''v'' move next, and points further  up/down by a fixed amount in V'
	    WRITE(6,*)'Use ''e'' to exit'
	    WRITE(6,'(A)')DEF_PEN
!
	    XVAL=5.0_LDP; CALL PGSCH(XVAL)
	    DO WHILE(1 .EQ. 1)				!Multiple cursor entries
	      CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	      WRITE(6,*)'Cursor values are:',XVAL,YVAL
	      IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')EXIT
	      IF(CURSVAL .EQ. 'r' .OR. CURSVAL .EQ. 'R')THEN
	        T1=XVAL
	        TMP_R(1:ND)=ABS(R(1:ND)-XVAL)
	        I=MINLOC(TMP_R(1:ND),IONE)
	        V(I)=YVAL
	        J=1; CALL PGSCI(J)
!	        SYMB_EXP_FAC=5.0; CALL PGSCH(SYMB_EXP_FAC)
	        CALL PGPT(IONE,XVAL,YVAL,IONE)
	        WRITE(6,*)'Replaced V for R=',R(I)
	      ELSE IF(CURSVAL .EQ. 'd' .OR. CURSVAL .EQ. 'D')THEN
	        T1=XVAL
	        TMP_R(1:ND)=ABS(R(1:ND)-T1)
	        I=MINLOC(TMP_R(1:ND),IONE)
	        R(I:ND-1)=R(I+1:ND)
	        V(I:ND-1)=V(I+1:ND)
	        ND=ND-1
	        WRITE(6,*)'Deleted R=',R(I),' from grid'
	      ELSE IF(CURSVAL .EQ. 'a' .OR. CURSVAL .EQ. 'A')THEN
	        DO I=1,ND-1
	          IF( (R(I)-XVAL)*(R(I+1)-XVAL) .LT. 0 )THEN
	            DO J=ND,I+1,-1
                      R(J+1)=R(J)
                      V(J+1)=V(J)
	            END DO
	            R(I+1)=XVAL; V(I+1)=YVAL
	            ND=ND+1
	            J=1; CALL PGSCI(J)
	            CALL PGPT(IONE,XVAL,YVAL,IONE)
	            WRITE(6,*)'Add R=',R(I),' to grid'
	            EXIT
	          END IF
	        END DO
	      ELSE IF(CURSVAL .EQ. 'm' .OR. CURSVAL .EQ. 'M')THEN
	        WRK_VEC(1:ND)=ABS(V(1:ND)-YVAL)
	        I=MINLOC(WRK_VEC(1:ND),IONE)
	        J=1; CALL PGSCI(J)
	        T1=XVAL-R(I)
	        DO J=1,I
	          R(J)=R(J)+T1
	          XVAL=R(J); YVAL=V(J)
	          CALL PGPT(IONE,XVAL,YVAL,IONE)
	        END DO	
	        WRITE(6,*)'Shited R grid by',T1
	      ELSE IF(CURSVAL .EQ. 'v' .OR. CURSVAL .EQ. 'V')THEN
	        WRK_VEC(1:ND)=ABS(R(1:ND)-XVAL)
	        I=MINLOC(WRK_VEC(1:ND),IONE)
	        J=1; CALL PGSCI(J)
	        T1=YVAL-V(I)
	        DO J=1,I
	          V(J)=V(J)+T1
	          XVAL=R(J); YVAL=V(J)
	          CALL PGPT(IONE,XVAL,YVAL,IONE)
	        END DO	
	        WRITE(6,*)'Shited R grid by',T1
	      ELSE
	        WRITE(6,*)RED_PEN
	        WRITE(6,*)'Error - use r(eplace), a(dd), d(elete), e'
	        WRITE(6,*)DEF_PEN
	      END IF
	    END DO
	    CALL GEN_IN(REPLOT,'Replot to see revision to V')
	    IF(REPLOT)THEN
	      CALL DP_CURVE(ND_OLD,R_PLT,V_PLT)
	      CALL DP_CURVE(ND,R,V)
 	      CALL GRAMON_PGPLOT(XLAB,YLAB,' ',' ')
	    ELSE
	      IF(LOG_X)THEN
	        R(1:ND)=(10.0_LDP**R(1:ND))*OLD_R(ND_OLD)
	      ELSE
	        R(1:ND)=R(1:ND)*OLD_R(ND_OLD)
	      END IF
	      IF(LOG_Y)THEN
	        V(1:ND)=(10.0_LDP**V(1:ND))
	      END IF
	      R(ND)=OLD_R(ND_OLD)
	      EXIT
	    END IF
	  END DO
!
! Compute SIGMA
!
	  ALLOCATE (COEF(ND,4))
	  CALL MON_INT_FUNS_V2(COEF,V,R,ND)
	  DO I=1,ND
	    SIGMA(I)=COEF(I,3)
	    SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
	  END DO
	  DEALLOCATE (COEF)
!
	ELSE IF(OPTION .EQ. 'PLOT')THEN
	  ND=ND_OLD
	  DO I=1,ND
	    R(I)=OLD_R(I)
	    V(I)=OLD_V(I)
	    SIGMA(I)=OLD_SIGMA(I)
	  END DO
	ELSE
	  WRITE(6,'(A)')'Option not recognixed'
	  STOP
	END IF
!
	IF(OPTION .NE. 'PLOT')THEN
	  NEW_RVSIG_FILE='RVSIG_COL_NEW'
	  CALL GEN_IN(NEW_RVSIG_FILE,'File for new R, V and sigma values')
	  OPEN(UNIT=10,FILE=NEW_RVSIG_FILE,STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(10,'(A)')'!'
	    IF(OPTION .EQ. 'MDOT')THEN
	      WRITE(10,'(A,ES14.6)')'! Old mass-loss rate in Msun/yr=',OLD_MDOT
	      WRITE(10,'(A,ES14.6)')'! New mass-loss rate in Msun/yr=',MDOT
	      WRITE(10,'(A,ES14.6)')'! Velocity at infinity in km/s =',VINF
	      WRITE(10,'(A,ES14.6)')'! Beta for velocity law        =',BETA
	      WRITE(10,'(A,I3)'    )'! Velocity law (type)          =',VEL_TYPE
	      WRITE(10,'(A,ES14.6)')'! Transition velocity is       =',V_TRANS
	      IF(VEL_TYPE .EQ. 5)THEN
	        WRITE(10,'(A,ES14.6)')'! Beta for 2nd vel. component  =',BETA2
	        WRITE(10,'(A,ES14.6)')'! VEXT for 2nd vel. component  =',VEXT
	        WRITE(10,'(A,ES14.6)')'! RP2  for 2nd vel. component  =',RP2/R(ND)
	        ELSE IF(VEL_TYPE .EQ. 3 .OR. VEL_TYPE .EQ. 3)THEN
	      WRITE(10,'(A,ES14.6)')'! Beta in inner wind           =',BETA2
	        WRITE(10,'(A,ES14.6)')'! ALPHA                        =',ALPHA
	        END IF
	    ELSE IF(OPTION .EQ. 'SCLR')THEN
	      WRITE(10,'(A,ES14.6)')'! Old stellar radius           =',OLD_R(ND_OLD)
	      WRITE(10,'(A,ES14.6)')'! New stellar rdaius           =',R(ND_OLD)
	      WRITE(10,'(A,ES14.6)')'! Velocity at infinity in km/s =',VINF
	      WRITE(10,'(A,ES14.6)')'! Beta for velocity law        =',BETA
	      WRITE(10,'(A,I3)'    )'! Velocity law (type)          =',VEL_TYPE
	      WRITE(10,'(A,ES14.6)')'! Transition velocity is       =',V_TRANS
	    END IF
	    WRITE(10,'(A,ES14.6)')'! R(1)/R(ND)                   =',R(1)/R(ND)
	    T1=R(1)/R(ND); T2=OLD_R(1)/OLD_R(ND_OLD); T3=1.0E-10_LDP
	    IF(.NOT. EQUAL(T1,T2,T3))THEN
	      WRITE(6,*)RED_PEN
	      WRITE(6,*)'Warning: update VADAT as R(1)/R(ND) has changed'
	      WRITE(6,*)DEF_PEN
	    END IF
	    T3=1.0E-02_LDP
	    IF(.NOT. EQUAL(V(1),OLD_V(1),T3))THEN
	      WRITE(6,*)RED_PEN
	      WRITE(6,*)'Warning: check if V(1) in VADAT needs to be updated'
	      WRITE(6,*)DEF_PEN
	    END IF
!
	    IF(RD_DENSITY)THEN
	      WRK_VEC(1:ND_OLD)=LOG(OLD_DENSITY(1:ND_OLD))
	      LOG_R(1:ND)=LOG(R(1:ND)); OLD_LOG_R(1:ND)=LOG(OLD_R(1:ND))
	      CALL MON_INTERP(DENSITY,ND,IONE,LOG_R,ND,WRK_VEC,ND_OLD,OLD_LOG_R,ND_OLD)
	      DENSITY(1:ND)=EXP(DENSITY(1:ND))
	      CALL MON_INTERP(CLUMP_FAC,ND,IONE,LOG_R,ND,OLD_CLUMP_FAC,ND_OLD,OLD_LOG_R,ND_OLD)
	    END IF
!
! Use a !! to inidicate from old file.
!
	    DO I=1,N_HEAD
	      WRITE(10,'(A,A)')'!',TRIM(OLD_HEADER(I))
	    END DO
	    WRITE(10,'(A)')'!'
	    WRITE(10,'(A,9X,A,19X,A,15X,A,5X,A)')'!','R','V(km/s)','Sigma','Depth'
	    WRITE(10,'(A)')'!'
	    WRITE(10,'(A)')' '
	    WRITE(10,'(I4,20X,A)')ND,'!Number of depth points`'
	    WRITE(10,'(A)')' '
	    IF(RD_DENSITY)THEN
	      DO I=1,ND
	        WRITE(10,'(F22.12,ES20.12,F17.7,2ES18.8,I4)')R(I),V(I),SIGMA(I),DENSITY(I),CLUMP_FAC(I),I
	      END DO
	    ELSE
	      DO I=1,ND
	        WRITE(10,'(F22.12,ES20.12,F17.7,4X,I4)')R(I),V(I),SIGMA(I),I
	      END DO
	    END IF
	  CLOSE(UNIT=10)
!
	  IF(RD_MEANOPAC)THEN
	    OPEN(UNIT=30,FILE='R_GRID_CHK',STATUS='UNKNOWN',ACTION='WRITE')
	    CALL MON_INTERP(TAU,ND,IONE,R,ND,OLD_TAU,ND_OLD,OLD_R,ND_OLD)
	    WRITE(6,*)'Determined TAU on the new grid'
	    I=1
	    WRITE(30,'(3X,A,6X,A,T28,A,T40,A,T53,A,T66,A,T78,A,T93,A,T101,A)')
	1           'I','R','V','SIGMA','TAU','dTAU','dTRAT','dV(I)','V(I)/V(I+1)'
	    FLUSH(UNIT=30)
	    T2=0.00_LDP
	    WRITE(30,'(1X,I3,F16.9,7ES13.4)')I,R(I),V(I),SIGMA(I),TAU(I),TAU(I+1)-TAU(I),T2,V(I)-V(I+1),V(I)/V(I+1)
	    DO I=2,ND-1
	      T1=TAU(I+1)-TAU(I)
	      T2=T1/(TAU(I)-TAU(I-1))
	      WRITE(30,'(1X,I3,F16.9,7ES13.4)')I,R(I),V(I),SIGMA(I),TAU(I),T1,T2,V(I)-V(I+1),V(I)/V(I+1)
	      FLUSH(UNIT=30)
	    END DO
	    I=ND
	    WRITE(30,'(1X,I3,F16.9,6ES13.4)')I,R(I),V(I),SIGMA(I),TAU(I)
	    CLOSE(UNIT=30)
	  END IF
	END IF
!
	T1=R(ND)
	R(1:ND)=R(1:ND)/T1
	T1=OLD_R(ND_OLD)
	OLD_R(1:ND_OLD)=OLD_R(1:ND_OLD)/T1
	WRITE(6,*)'Plotting V versus R/R*'
	CALL DP_CURVE(ND_OLD,OLD_R,OLD_V)
	CALL DP_CURVE(ND,R,V)
	CALL GRAMON_PGPLOT('R/R\d*\u','V(km/s)',' ',' ')
	WRITE(6,*)'Plotting Sigma versus R/R*'
	CALL DP_CURVE(ND_OLD,OLD_R,OLD_SIGMA)
	CALL DP_CURVE(ND,R,SIGMA)
	CALL GRAMON_PGPLOT('R/R\d*\u','SIGMA',' ',' ')
!
	STOP
	END
!
CONTAINS
	SUBROUTINE DESCRIBE_VEL_LAWS()	
	USE SET_KIND_MODULE
	USE MOD_COLOR_PEN_DEF
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')'Type 1: W(r).V(r) with'
	WRITE(6,'(A)')'            W(r) = VINF*(1-rt/r)**BETA'
	WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')'Type 2: W(r).V(r) with'
	WRITE(6,'(A)')'            W(r) = 2*VTRANS + (VINF-2*VTRANS)*(1-rt/r)**BETA'
	WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')'Type 3: W(r).V(r) with'
        WRITE(6,'(A)')'              X  = 1-rt/r'
	WRITE(6,'(A)')'            W(r) = 2*VTRANS + (VINF-2*VTRANS)*X**[BETA+(BETA2-BETA)*X]'
	WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')'Type 4: W(r).V(r) with'
        WRITE(6,'(A)')'              X  = 1-rt/r'
	WRITE(6,'(A)')'            W(r) = 3*VTRANS + (VINF-2*VTRANS)*X**[BETA+(BETA2-BETA)*X]'
	WRITE(6,'(A)')'        and V(r) = 1.0D0+exp( (rt-r)/h )'
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')'Type 5: (W(r)+E(r))/V(r) with'
	WRITE(6,'(A)')'            W(r) = 2*VTRANS + (VINF-2*VTRANS)*(1-rt/r)**BETA'
	WRITE(6,'(A)')'            E(r)= VEXT*(1-RP2/r)**BETA2 '
	WRITE(6,'(A)')'        and V(r) = 1.0D0+2*exp( (rt-r)/h )'
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')'Type 6: Attach existing existing velocity law at connection depth'
	WRITE(6,'(A)')DEF_PEN
!
	RETURN
	END

	SUBROUTINE CALCULATE_VEL(R,V,SIGMA,ND)
	USE SET_KIND_MODULE
	USE VEL_LAW_PARAMS
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
!
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) TOP,BOT
	REAL(KIND=LDP) dTOPdR,dBOTdR
!
	INTEGER I
	INTEGER V_TYPE
!
	IF(VEL_TYPE .EQ. 1)THEN
	  RO = R_TRANS * (1.0_LDP - (2.0_LDP*V_TRANS/VINF)**(1.0_LDP/BETA) )
	  T1= R_TRANS * dVdR_TRANS / V_TRANS
	  SCALE_HEIGHT =  0.5_LDP*R_TRANS / (T1 - BETA*RO/(R_TRANS-RO) )
!
	  WRITE(6,*)'  Transition radius is',R_TRANS
	  WRITE(6,*)'Transition velocity is',V_TRANS
	  WRITE(6,*)'                 R0 is',RO
	  WRITE(6,*)'       Scale height is',SCALE_HEIGHT
!
	  DO I=1,TRANS_I-1
            T1=RO/R(I)
            T2=1.0_LDP-T1
            TOP = VINF* (T2**BETA)
            BOT = 1.0_LDP + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
            V(I) = TOP/BOT

!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.

            dTOPdR = VINF * BETA * T1 / R(I) * T2**(BETA - 1.0_LDP)
            dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT )  / SCALE_HEIGHT
            dVdR = dTOPdR / BOT  + V(I)*dBOTdR/BOT
            SIGMA(I)=R(I)*dVdR/V(I)-1.0_LDP
	  END DO
!
	ELSE IF(VEL_TYPE .EQ. 2)THEN
	  SCALE_HEIGHT = V_TRANS / (2.0_LDP * DVDR_TRANS)
	  WRITE(6,*)'  Transition radius is',R_TRANS
	  WRITE(6,*)'Transition velocity is',V_TRANS
	  WRITE(6,*)'       Scale height is',SCALE_HEIGHT
	  DO I=1,TRANS_I-1
	    T1=R_TRANS/R(I)
	    T2=1.0_LDP-T1
	    TOP = 2.0_LDP*V_TRANS + (VINF-2.0_LDP*V_TRANS) * T2**BETA
	    BOT = 1.0_LDP + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
	    V(I) = TOP/BOT

!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.

	    dTOPdR = (VINF - 2.0_LDP*V_TRANS) * BETA * T1 / R(I) * T2**(BETA - 1.0_LDP)
	    dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
	    dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
            SIGMA(I)=R(I)*dVdR/V(I)-1.0_LDP
	  END DO
!
	ELSE IF(VEL_TYPE .EQ. 3 .OR. VEL_TYPE .EQ. 4)THEN
	  SCALE_HEIGHT = V_TRANS / (2.0_LDP * DVDR_TRANS)
	  WRITE(6,*)'  Transition radius is',R_TRANS
	  WRITE(6,*)'Transition velocity is',V_TRANS
	  WRITE(6,*)'       Scale height is',SCALE_HEIGHT
	  ALPHA=2.0_LDP
	  IF(VEL_TYPE .EQ. 4)ALPHA=3.0_LDP
	  CALL GEN_IN(BETA2,'Beta2 for velocity law')
	  DO I=1,TRANS_I-1
	    T1=R_TRANS/R(I)
	    T2=1.0_LDP-T1
	    T3=BETA+(BETA2-BETA)*T2
	    TOP = (VINF-ALPHA*V_TRANS) * T2**T3
	    BOT = 1.0_LDP + (ALPHA-1.0_LDP)*exp( (R_TRANS-R(I))/SCALE_HEIGHT )

!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.

	    dTOPdR = (VINF - ALPHA*V_TRANS) * BETA * T1 / R(I) * T2**(T3-1.0_LDP) +
	1                  T1*TOP*(BETA2-BETA)*(1.0_LDP+LOG(T2))/R(I)
	    dBOTdR=  (ALPHA-1.0_LDP)*exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
!
	    TOP = ALPHA*V_TRANS + TOP
	    dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
	    V(I) = TOP/BOT
            SIGMA(I)=R(I)*dVdR/V(I)-1.0_LDP
	  END DO
!
	ELSE IF(VEL_TYPE .EQ. 5)THEN
	  BETA2=BETA; VEXT=0.1_LDP*VINF; RP2=2.0_LDP
	  CALL GEN_IN(VEXT,'Additonal V component (add to VINF)- VEXT')
	  CALL GEN_IN(BETA2,'Beta2 for velocity law')
	  CALL GEN_IN(RP2,'RP2 for velocity law (in terms of R(ND)')
	  RP2=RP2*R(ND)
	  SCALE_HEIGHT = V_TRANS / (2.0_LDP * DVDR_TRANS)
	  WRITE(6,*)'  '
	  WRITE(6,*)'  Transition radius is',R_TRANS
	  WRITE(6,*)'Transition velocity is',V_TRANS
	  WRITE(6,*)'       Scale height is',SCALE_HEIGHT
	  WRITE(6,*)'  '
	  DO I=1,TRANS_I-1
	    T1=R_TRANS/R(I)
	    T2=1.0_LDP-T1
	    TOP = 2.0_LDP*V_TRANS + (VINF-2.0_LDP*V_TRANS) * T2**BETA
	    BOT = 1.0_LDP + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
            IF(RP2/R(I) .LT. 1.0_LDP)THEN
               TOP=TOP+ VEXT*(1.0_LDP-RP2/R(I))**BETA2
               dTOPdR=(RP2/R(I)/R(I))*BETA2*VEXT*(1.0_LDP-RP2/R(I))**(BETA2-1)
            ELSE
               dTOPdR=0.0_LDP
	    END IF
	    V(I) = TOP/BOT
!
!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.
!
	    dTOPdR = dTOPdR + (VINF - 2.0_LDP*V_TRANS) * BETA * T1 / R(I) * T2**(BETA - 1.0_LDP)
	    dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
	    dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
            SIGMA(I)=R(I)*dVdR/V(I)-1.0_LDP
	  END DO
	ELSE
	  WRITE(6,*)'Unrecognized veloity type: your type is',VEL_TYPE
	END IF
!
	RETURN
	END
