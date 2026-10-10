!
	SUBROUTINE WRITE_LINE_DATA_V3(R,VEL,T,SIGMA,ESEC,BETA_VEC,
	1                             ES_J,CRSE_FREQ,NF_CRSE,
	1                             ND,NC,NP,NPHI,NBETA)
        USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 01-May-2022 -- Based on write_line_data_cl.f
!
! Input values
!
	INTEGER NPHI
	INTEGER NBETA
	INTEGER ND
	INTEGER NC
	INTEGER NP
	INTEGER NF_CRSE   
!
	REAL(KIND=LDP)  R(ND)
	REAL(KIND=LDP)  VEL(ND)
	REAL(KIND=LDP)  T(ND)
	REAL(KIND=LDP)  SIGMA(ND)
	REAL(KIND=LDP)  ESEC(ND)
	REAL(KIND=LDP)  BETA_VEC(NBETA)
	REAL(KIND=LDP)  ES_J(NF_CRSE,NBETA,ND)
!
	LOGICAL RAND_IN_BETA

	REAL(KIND=LDP)  TA(ND)
	REAL(KIND=LDP)  P(NP)
!
	REAL(KIND=LDP)  PHI_VEC(NPHI)
!
	REAL(KIND=LDP)  CRSE_FREQ(NF_CRSE)
	REAL(KIND=LDP)  DEL_NU
!
	INTEGER NPNT_CRSE
	INTEGER NC_CRSE
	INTEGER NP_CRSE
	INTEGER NP_INS
!
	REAL(KIND=LDP)  PI
	REAL(KIND=LDP)  T1, T2, T3
	INTEGER I, J, K, ML
!
	CHARACTER FORMFEED*2
	CHARACTER TIME*20,FMT*132
	CHARACTER*12 PRODATE
	PARAMETER (PRODATE='13-Nov-2018') !Must be changed after alterations
!
	INTEGER, PARAMETER :: LU   =300
	INTEGER, PARAMETER :: LUOUT=9
        INTEGER, PARAMETER :: LUMOM=10
        INTEGER, PARAMETER :: LUMOD=20
        INTEGER, PARAMETER :: LUOBS=23
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	  I=12
          FORMFEED=' '//CHAR(I)
	  TA(1:ND) = 0.0E+00
!
          PI=ACOS(-1.0D0)
          DO I=1,NPHI
            PHI_VEC(I)=(I-1)*PI/(NPHI-1)
          END DO
!
! Create the file LINE_MOM_DATA here
!
	  CALL DATE_TIME(TIME)
          OPEN(UNIT=LUMOD,FILE='LINE_MOM_DATA',STATUS='UNKNOWN')
          WRITE(LUMOD,'(/,'' Model Started on:'',15X,(A))')TIME
          WRITE(LUMOD,
	1       '('' Main program last changed on:'',3X,(A))')PRODATE
          WRITE(LUMOD,'('' Main Format Date:'',15X,(A))')'05-Feb-1994'
!
! We will use a blank line in the outputfile to indicate the end of a section.
!
! *****************************
!
! Spatial Co-ordinates.
!
          WRITE(LUMOD,'()')
          FMT='(7X,I8,3X,''!Number of depth points [ND]'')'
          WRITE(LUMOD,FMT)ND
          FMT='(7X,I8,3X,''!Number of polar angles [NBETA]'')'
          WRITE(LUMOD,FMT)NBETA
!
! Frequency dependance.
!
          WRITE(LUMOD,'()')
          FMT='(7X,I8,3X,''!Number of CMF frequencies [NF_CRSE]'')'
          WRITE(LUMOD,FMT)NF_CRSE
!
! Co-ordiantes for angular description of radiation field (at a given
! spatial depth).
!
          WRITE(LUMOD,'(7X,I8,3X,A)')NPHI,'!Number of azimuthal angles [NPHI]'
          WRITE(LUMOD,'(7X,I8,3X,A)')0,'!Number of point source rays [NPNT_CRSE]'
          WRITE(LUMOD,'(7X,I8,3X,A)')NC,'!Number of core rays [NC_CRSE]'
          WRITE(LUMOD,'(7X,I8,3X,A)')NP,'!Total number of rays [NP_CRSE]'
          WRITE(LUMOD,'(7X,I8,3X,A)')0,'!Number of rays inserted [NP_INS]'
          WRITE(LUMOD,'(7X,I8,3X,A)')ND,'!Number of ray points [NR_MAX]'
!
! *****************************
!
! Compute coarse impact grid.
!
	  CALL IMPAR(P,R,R(ND),NC,ND,NP)
!
          CALL RITE_2DJ_ASC_JH(R,ND,1,'Radius Coordinates - [R]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(BETA_VEC,NBETA,1,'Polar Coordinates - [BETA]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(CRSE_FREQ,NF_CRSE,1,'Coarse CMF frequency grid - [FREQ]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(PHI_VEC,NPHI,1,'Azimuthal Coordinates - [PHI]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(P,NP,1,'Impact Parameters - [P]:',LUMOD)
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A)')FORMFEED
!
! *****************************
!
          WRITE(LUMOD,'()')
          FMT='('' Boundary Parameters Format Date:'',14X,(A))'
          WRITE(LUMOD,FMT)'06-Dec-1993'
          FMT='(14X,L1,3X,''!Optically thin core?'')'
          WRITE(LUMOD,FMT).FALSE.   !THIN_CORE
          FMT='(14X,L1,3X,''!Schuster boundary condition at core?'')'
          WRITE(LUMOD,FMT).FALSE.   !SCHUS_BC
          FMT='(X,1PE14.6,3X,''!Schuster intensity'')'
          WRITE(LUMOD,FMT)1.0E+00   !I_SCHUS
!
          FMT='(14X,L1,3X,''!Diffusion boundary condition at core?'')'
          WRITE(LUMOD,FMT).TRUE.    !DIFF
          FMT='(X,1PE14.6,3X,''!Diffusion flux'')'
          WRITE(LUMOD,FMT)0.0E+00   !DBB
!
          FMT='(14X,L1,3X,''!Point source?'')'
          WRITE(LUMOD,FMT).FALSE.   !PNT_SRCE
          FMT='(X,1PE14.6,3X,''!Radius of point source'')'
          WRITE(LUMOD,FMT)0.0E+00   !R_PNT_SRCE
          FMT='(X,1PE14.6,3X,''!B Point source'')'
          WRITE(LUMOD,FMT)0.0E+00   !B_PNT_SRCE
!
          FMT='(14X,L1,3X,''!Optically thick at outer boundary?'')'
          WRITE(LUMOD,FMT).FALSE.   !THK
!
! *****************************
!
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A,15X,A)')' Angular Parameters Format Date:','30-Dec-1993'
          WRITE(LUMOD,'(A,15X,A)')' Angular Law:','POW_E3' !DEN_FN_NAME
          FMT='(7X,I8,3X,''!Number of density parameters'')'
          WRITE(LUMOD,FMT) 3               !N_FN_PRMS
          DO I=1,3
            FMT='(X,1PE14.6,3X,''!A'',I1)'
            WRITE(LUMOD,FMT)0.0E+0,I       !DEN_FN_PRMS(I),I
          END DO
!
          WRITE(LUMOD,'(X,1PE14.6,3X,A)')0.0,'!TH_ANG_EXP'
          WRITE(LUMOD,'(X,1PE14.6,3X,A)')0.0,'!ETA_ANG_EXP'
          WRITE(LUMOD,'(X,1PE14.6,3X,A)')0.0,'!ESEC_ANG_EXP'
          WRITE(LUMOD,'(X,1PE14.6,3X,A)')0.0,'!CHIL_ANG_EXP'
          WRITE(LUMOD,'(X,1PE14.6,3X,A)')0.0,'!ETAL_ANG_EXP'
!
! *****************************
!
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A,15X,A)')'Opacity Parameters Format Date:','06-Dec-1993'
          WRITE(LUMOD,'(1X,A,15X,A)')'Opacity Law:','File'
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Constant for thermal Opacity'		!MCHI(1)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Exponent for thermal Opacity'		!MCHI(2)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Constant for thermal Emissivity'		!META(1)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Exponent for thermal Emissivity'		!META(2)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Constant for electron scattering'		!MESEC(1)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Exponent for electron scattering'		!MESEC(2)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Constant for line Opacity'			!MCHIL(1)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Exponent for line Opacity'			!MCHIL(2)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Constant for line Emissivity'		!MCHIL(1)
          WRITE(LUMOD,'(1X,ES14.6,3X,A)')0.0E+00,'!Exponent for line Emissivity'		!MCHIL(2)
!
! *****************************
!
	  TA(1:ND) = 0.0E+00
          CALL RITE_2DJ_ASC_JH(TA,ND,1,'Thermal Opacity.',LUMOD)
          CALL RITE_2DJ_ASC_JH(TA,ND,1,'Thermal Emissivity.',LUMOD)
          CALL RITE_2DJ_ASC_JH(ESEC,ND,1,'Electron Scattering Opacity.',LUMOD)
          CALL RITE_2DJ_ASC_JH(TA,ND,1,'Line Opacity.',LUMOD)
          CALL RITE_2DJ_ASC_JH(TA,ND,1,'Line Emissivity.',LUMOD)
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A)')FORMFEED
!
! *****************************
!
          WRITE(LUMOD,'()')
          FMT='('' Structure Parameters Format Date:'',15X,(A))'
          WRITE(LUMOD,FMT)'06-Dec-1993'
!
          CALL RITE_2DJ_ASC_JH(T,ND,1,'Temperature (10^4K)',LUMOD)
          CALL RITE_2DJ_ASC_JH(VEL,ND,1,'Velocity (km/s)',LUMOD)
          CALL RITE_2DJ_ASC_JH(SIGMA,ND,1,'Velocity gradient [dlnV/Dlnr-1]',LUMOD)
!
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A)')' Moments of Radiation field:'
          FMT='(5X,I8,5X,''!Number of radiation moments'')'
          WRITE(LUMOD,FMT)5
          I=NBETA*NF_CRSE       !Temporary Variable
          CALL RITE_2DJ_ASC_JH(ES_J,I,ND,'J_L moment.',LUMOD)
          CALL RITE_2DJ_ASC_JH(ES_J,I,ND,'J_R moment.',LUMOD)
	  WRITE(6,*)'ABOUT TO WRITE OUT K20'; FLUSH(UNIT=6)
!
! Compute K moment assuming J/3
!
	  ES_J=ES_J/3.0_LDP 

          I=NBETA*NF_CRSE       !Temporary Variable
          CALL RITE_2DJ_ASC_JH(ES_J,I,ND,'K20_L moment.',LUMOD)
	  WRITE(6,*)'DONE'; FLUSH(UNIT=6)
!
	  ES_J=0.0_LDP
          CALL RITE_2DJ_ASC_JH(ES_J,I,ND,'GT moment.',LUMOD)
          CALL RITE_2DJ_ASC_JH(ES_J,I,ND,'KT moment.',LUMOD)
!
	  RETURN
!
	END SUBROUTINE
