C
C Routine to plot J from EDDFACTOR file. This J is convolved with the
C electron redistribution function using the 1-parameter formulation of
C Hummer and Rybicki. For comparison, RJ_ES from the ES_J_CONV file,
C may also be plotted.
C
	PROGRAM PLT_RJ
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
C
	REAL(KIND=LDP), ALLOCATABLE :: NU(:)
	REAL(KIND=LDP), ALLOCATABLE :: RJ(:,:)
C
	REAL(KIND=LDP), ALLOCATABLE :: RJ_ES_RD(:,:)
C
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_RJ(:)
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_ES(:)
C
	REAL(KIND=LDP), ALLOCATABLE :: PLANCK_FN(:)
C
	REAL(KIND=LDP), ALLOCATABLE :: A(:)
	REAL(KIND=LDP), ALLOCATABLE :: B(:)
	REAL(KIND=LDP), ALLOCATABLE :: C(:)
	REAL(KIND=LDP), ALLOCATABLE :: D(:)
C
	REAL*4, ALLOCATABLE :: XV(:)
	REAL*4, ALLOCATABLE :: YV(:)
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_OUT=6
C
	REAL(KIND=LDP) RJ_FLUX,ES_FLUX,ES2_FLUX
C
	REAL(KIND=LDP) T_ELEC
	REAL(KIND=LDP) BETA
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) D1,D2,DH
	REAL(KIND=LDP) SPEED_OF_LIGHT,C_KMS,NU_0
	EXTERNAL SPEED_OF_LIGHT
C
	INTEGER I,J,K,ML
	INTEGER ND,NCF
	INTEGER ND2,NCF2
	INTEGER LU_IN
	INTEGER IOS
	INTEGER REC_LENGTH
C
	CHARACTER*80 FILENAME
	CHARACTER*80 FILE_DATE
C
	C_KMS=SPEED_OF_LIGHT()/1.0E+05_LDP
C
	WRITE(T_OUT,*)
	1  ' Routine to plot J from EDDFACTOR file. This J is convolved'
	WRITE(T_OUT,*)
	1  ' with the electron redistrbution function using the 1-parameter'
	WRITE(T_OUT,*)
	1  ' formulation of Hummer and Rybicki. For comparison, RJ_ES from'
	WRITE(T_OUT,*)
	1  ' the ES_J_CONV file, may also be plotted.'
!
!	CALL GEN_IN(ND,'Number of model depth points')
!	CALL GEN_IN(NCF,'Number of model frequency points')
!
	LU_IN=35
	T_elec=1.0_LDP		!10^4 K
C
C Open EDDFACTOR file.
C
1	FILENAME='EDDFACTOR'
	CALL GEN_IN(FILENAME,'Filename with RJ')
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to open INFO file for J data'
	  STOP
	END IF
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',RECL=REC_LENGTH,
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ',ERR=1)
	  READ(LU_IN,REC=3)I,NCF,ND
	  ALLOCATE ( NU(NCF) )
	  ALLOCATE ( RJ(ND,NCF) )
	  DO ML=1,NCF
	    READ(LU_IN,REC=I+ML-1)(RJ(K,ML),K=1,ND),NU(ML)
	  END DO
	CLOSE(LU_IN)
!
! Allocate other varrays.
!
	ALLOCATE ( RJ_ES_RD(ND,NCF) )
C
	ALLOCATE ( FLUX_RJ(ND) )
	ALLOCATE ( FLUX_ES(ND) )
C
	ALLOCATE ( PLANCK_FN(NCF) )
	ALLOCATE ( A(NCF) )
	ALLOCATE ( B(NCF) )
	ALLOCATE ( C(NCF) )
	ALLOCATE ( D(NCF) )
C
	ALLOCATE ( XV(NCF) )
	ALLOCATE ( YV(NCF) )
!
! 
!
	WRITE(T_OUT,*)'If ES_J_CONV unavailable, enter EDDFACTOR again'
2	FILENAME='ES_J_CONV'
	CALL GEN_IN(FILENAME,'Filename with RJ_ES')
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',RECL=REC_LENGTH,
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ',ERR=2)
	  READ(LU_IN,REC=3)I,NCF2,ND2
	  IF(ND2 .NE. ND .OR. NCF2 .NE. NCF)THEN
	     WRITE(T_OUT,*)'Error: ND and NCF not the same'
	     WRITE(T_OUT,*)'ND=',ND,'NCF=',NCF
	     WRITE(T_OUT,*)'ND2=',ND2,'NCF2=',NCF2
	     STOP
	   END IF
	   DO ML=1,NCF
	     READ(LU_IN,REC=I+ML-1)(RJ_ES_RD(K,ML),K=1,ND),T1
	   END DO
	CLOSE(LU_IN)
C
C Compute integrals as a function of depth to check flux conservation.
C
	FLUX_RJ(1:ND)=+0.0_LDP
	FLUX_ES(1:ND)=+0.0_LDP
	DO ML=2,NCF-1
	  DO I=1,ND
	    FLUX_RJ(I)=FLUX_RJ(I)+(NU(ML)-NU(ML+1))*(RJ(I,ML)+
	1                  RJ(I,ML+1))
	    FLUX_ES(I)=FLUX_ES(I)+(NU(ML)-NU(ML+1))*(RJ_ES_RD(I,ML)+
	1                  RJ_ES_RD(I,ML+1))
	  END DO
	END DO
	WRITE(9,*)'% Flux error '
	WRITE(9,'(1P5E12.5)')200.0D0*(FLUX_RJ(1:ND)-FLUX_ES(1:ND))/
	1                    (FLUX_RJ(1:ND)+FLUX_ES(1:ND))
C
	FLUX_RJ(1:ND)=+0.0_LDP
	FLUX_ES(1:ND)=+0.0_LDP
	DO ML=2,NCF-1
	  DO I=1,ND
	    T1=LOG(NU(ML)/NU(ML+1))
	    FLUX_RJ(I)=FLUX_RJ(I)+T1*(RJ(I,ML)+RJ(I,ML+1))
	    FLUX_ES(I)=FLUX_ES(I)+T1*(RJ_ES_RD(I,ML)+RJ_ES_RD(I,ML+1))
	  END DO
	END DO
	WRITE(9,*)' '
	WRITE(9,*)' '
	WRITE(9,*)'% Photon error '
	WRITE(9,'(1P5E12.5)')200.0D0*(FLUX_RJ(1:ND)-FLUX_ES(1:ND))/
	1                    (FLUX_RJ(1:ND)+FLUX_ES(1:ND))
C
	K=ND
200	CALL GEN_IN(K,'Input depth for plotting (0 to stop)')
	IF(K .EQ. 0)STOP
C
	NU_0=0
	IF(K .LT. 0)THEN
	  CALL GEN_IN(NU_0,'Input central frequecy for velocity plot')
	  K=ABS(K)
	END IF
	DO ML=1,NCF
	  XV(ML)=NU(ML)
	  YV(ML)=LOG10(RJ(K,ML))
	END DO
	IF(NU_0 .NE. 0)XV(1:NCF)=C_KMS*(NU(1:NCF)-NU_0)/NU_0
	CALL CURVE(NCF,XV,YV)
C
	DO ML=1,NCF
	  YV(ML)=LOG10(RJ_ES_RD(K,ML))
	END DO
	CALL CURVE(NCF,XV,YV)
C
C Convolve RJ with electron-scattering redistribution function.
C
	CALL GEN_IN(T_elec,'Electron temperature (in units 10^4 K)')
	BETA=1.84E-03_LDP*SQRT(T_elec)
	T1=0.5_LDP*BETA*BETA
	C(1)=0.0_LDP
	A(1)=0.0_LDP
	B(1)=-1.0_LDP
	D(1)=RJ(K,1)
	DO ML=2,NCF-1
	  D1=LOG(NU(ML-1)/NU(ML))
	  D2=LOG(NU(ML)/NU(ML+1))
	  DH=0.5_LDP*(D1+D2)
	  A(ML)=-T1/D1/DH
	  B(ML)=-1.0_LDP
	  C(ML)=-T1/D2/DH
	  D(ML)=RJ(K,ML)
	END DO
	A(NCF)=0.0
	C(NCF)=0.0
	B(NCF)=-1.0_LDP
	D(NCF)=RJ(K,NCF)
C
	CALL THOMAS_RH(A,B,C,D,NCF,IONE)
	DO I=1,NCF
	  YV(I)=LOG10(D(I))
	END DO
	CALL CURVE(NCF,XV,YV)
!
! Compute integrals as a function of depth to check flux conservation.
!
	RJ_FLUX=+0.0_LDP
	ES_FLUX=+0.0_LDP
	ES2_FLUX=+0.0_LDP
	DO ML=2,NCF-1
	    RJ_FLUX=RJ_FLUX+(NU(ML)-NU(ML+1))*(RJ(K,ML)+RJ(K,ML+1))
	    ES_FLUX=ES_FLUX+(NU(ML)-NU(ML+1))*(RJ_ES_RD(K,ML)+RJ_ES_RD(K,ML+1))
	    ES2_FLUX=ES2_FLUX+(NU(ML)-NU(ML+1))*(D(ML)+D(ML+1))
	END DO
	ES_FLUX=200.0_LDP*(RJ_FLUX-ES_FLUX)/(RJ_FLUX+ES_FLUX)
	ES2_FLUX=200.0_LDP*(RJ_FLUX-ES2_FLUX)/(RJ_FLUX+ES2_FLUX)
	WRITE(5,'(A,1X,1P,2E13.5)')'  %Flux errors:  ',ES_FLUX,ES2_FLUX
C
	RJ_FLUX=+0.0_LDP
	ES_FLUX=+0.0_LDP
	ES2_FLUX=+0.0_LDP
	DO ML=2,NCF-1
	  DO I=1,ND
	    T1=LOG(NU(ML)/NU(ML+1))
	    RJ_FLUX=RJ_FLUX+T1*(RJ(K,ML)+RJ(K,ML+1))
	    ES_FLUX=ES_FLUX+T1*(RJ_ES_RD(K,ML)+RJ_ES_RD(K,ML+1))
	    ES2_FLUX=ES2_FLUX+T1*(D(ML)+D(ML+1))
	  END DO
	END DO
	ES_FLUX=200.0_LDP*(RJ_FLUX-ES_FLUX)/(RJ_FLUX+ES_FLUX)
	ES2_FLUX=200.0_LDP*(RJ_FLUX-ES2_FLUX)/(RJ_FLUX+ES2_FLUX)
	WRITE(5,'(A,1X,1P,2E13.5)')'  %Photon errors:',ES_FLUX,ES2_FLUX
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'Plotting J, Jes, Jes'
	CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)',
	1          'Log J(ergs cm\u-2\d s\u-1\d Hz\u-1\d)',' ',' ')
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'Plotting normalized % difference between J_RD, Jes_RD'
	DO ML=1,NCF
	  YV(ML)=100.0_LDP*(RJ(K,ML)-RJ_ES_RD(K,ML))/RJ(K,ML)
	END DO
	CALL CURVE(NCF,XV,YV)
	CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','%Difference',' ',' ')
!
! Covolve Planck function and plot
!
	BETA=1.84E-03_LDP*SQRT(T_elec)
	T1=0.5_LDP*BETA*BETA
	C(1)=0.0_LDP
	A(1)=0.0_LDP
	B(1)=-1.0_LDP
	D(1)=NU(1)**3/(EXP(4.7994_LDP*NU(1)/T_ELEC)-1)
	DO ML=2,NCF-1
	  D1=LOG(NU(ML-1)/NU(ML))
	  D2=LOG(NU(ML)/NU(ML+1))
	  DH=0.5_LDP*(D1+D2)
	  A(ML)=-T1/D1/DH
	  B(ML)=-1.0_LDP
	  C(ML)=-T1/D2/DH
	  D(ML)=NU(ML)**3/(EXP(4.7994_LDP*NU(ML)/T_ELEC)-1)
	END DO
	A(NCF)=0.0
	C(NCF)=0.0
	B(NCF)=-1.0_LDP
	D(NCF)=NU(NCF)**3/(EXP(4.7994_LDP*NU(NCF)/T_ELEC)-1)
	PLANCK_FN(1:NCF)=D(1:NCF)
C
	CALL THOMAS_RH(A,B,C,D,NCF,IONE)
C
	PLANCK_FN=LOG10(PLANCK_FN)
	D=LOG10(D)
C
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'Plotting Log(B_es/B)'
	DO I=1,NCF
	  YV(I)=D(I)-PLANCK_FN(I)
	END DO
	CALL CURVE(NCF,XV,YV)
	CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','Log(B_es/B)',' ',' ')
!
	CALL CURVE(NCF,XV,YV)
!
! Shift and scale e.s. funtion so still have Placnk function
!
	A(1)=NU(1)
	I=2
	ML=1
	DO WHILE(NU(I) .GT. 1.5_LDP*T_ELEC)
	  DO WHILE(D(I) .GT. PLANCK_FN(ML))
	    ML=ML+1
	  END DO
	  DO WHILE(D(I) .LT. PLANCK_FN(ML-1))
	    ML=ML+1
	  END DO
	  T1=(D(I)-PLANCK_FN(ML))/(PLANCK_FN(ML-1)-PLANCK_FN(ML))
          A(I)=T1*NU(ML-1)+(1.0_LDP-T1)*NU(ML)
	  I=I+1
	END DO
C
	T1=1/NU(I-1)-1/A(I-1)		!Wavelength shift
	DO J=I,NCF
	  A(J)=1.0_LDP/( 1.0_LDP/NU(J) - T1)
	END DO
C
C Perfrom a simple linear interpolation back onto the old frequency grid.
C
	B(1)=D(1)
	I=1
	DO ML=2,NCF-2
	  DO WHILE (NU(ML) .LT. A(I+1))
	     I=I+1
	  END DO
	  T1=(NU(ML)-A(I))/(A(I+1)-A(I))
	  B(ML)=T1*D(I+1)+(1.0_LDP-T1)*D(I)
	END DO
	B(NCF-1:NCF)=D(NCF-1:NCF)
C
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'Plotting scaled/shifted Log(B_es/B)'
	DO I=1,NCF
	  YV(I)=B(I)-PLANCK_FN(I)
	END DO
	CALL CURVE(NCF,XV,YV)
C
	CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)',' ',' ',' ')
	GOTO 200
C
	END
C
C Subroutine to solve a tridiagonal system of N1 simultaneous
C Equations which are tridiagonal in nature for N2 R.H.sides.
C
C The equations to be solved are assumed to have the form
C
C    A(i).X(i-1) - [H(i)+A(i)+B(i)].X(i+1) - C(i).X(i+1)
C
	SUBROUTINE THOMAS_RH(A,H,C,D,N1,N2)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	INTEGER N1,N2
	REAL(KIND=LDP) A(N1),H(N1),C(N1),D(N1,N2)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
	INTEGER I,J
	REAL(KIND=LDP) DIV(N1)
C
C
C Compute quantities that will be used repeatedly if the same tridiagonal
C system is used for many R.H. Sides.
C
	C(N1)=0				!As used.
	DIV(1)=1.0_LDP/(C(1)+H(1))
	C(1)=C(1)*DIV(1)
	H(1)=H(1)*DIV(1)
	DO I=2,N1
	  DIV(I)=1.0_LDP/(A(I)*H(I-1)+H(I)+C(I))
	  H(I)=(A(I)*H(I-1)+H(I))*DIV(I)
	  C(I)=C(I)*DIV(I)
	END DO
C
C Entry for Thomas algorithim when H,C have been previously modified.
C
	ENTRY SIMPTH_RH(A,H,C,D,N1,N2)
C
	DO J=1,N2
	  D(1,J)=D(1,J)*DIV(1)
	  DO I=2,N1
	    D(I,J)=(D(I,J)+A(I)*D(I-1,J))*DIV(I)
	  END DO
	END DO
C
	DO J=1,N2
	  D(N1,J)=-D(N1,J)
	  DO I=N1-1,1,-1
	    D(I,J)=C(I)*D(I+1,J)-D(I,J)
	  END DO
	END DO
C
	RETURN
	END
