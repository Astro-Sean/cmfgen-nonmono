!
      SUBROUTINE DEFINE_GRID_V2_GAM(R_EXT,V_EXT,VDOP_VEC,VDOP_FRAC,ND_EXT,R,P,ND,NC,NP,ANG_MULT)
	USE SET_KIND_MODULE
      USE MOD_SPACE_GRID_V2_GAM
      USE GAM_MU_MOD
      USE MOD_RD_GAMRAY_CNTRL_VARIABLES, ONLY : VERBOSE_GAMMA
      IMPLICIT NONE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Determine radial grid from given input parameters
! Written   11-94 DLM
! Altered 2/21/96 DLM Updated to F90 standard
! Altered 6/12/97 DLM Changed name from rpz.f
! Altered 26/Nov/98 DJH: Changed to V2
!                        Altered to allow an extended R_GRID (R_EXT).
!--------------------------------------------------------------------
!
! Grid size variables: passed from calling routine
!
      INTEGER NP
      INTEGER ND
      INTEGER NC
!
      INTEGER ND_EXT
      REAL(KIND=LDP) R_EXT(ND_EXT)
      REAL(KIND=LDP) V_EXT(ND_EXT)
      REAL(KIND=LDP) VDOP_VEC(ND_EXT)
      REAL(KIND=LDP) VDOP_FRAC
!
! Grid variables
!
      REAL(KIND=LDP) R(ND)
      REAL(KIND=LDP) V(ND)
      REAL(KIND=LDP) P(NP)
!
! Local variables
!
      REAL(KIND=LDP) MU(NP)
      REAL(KIND=LDP) TEMP_VEC(NP)
!
      INTEGER ID,IP
      INTEGER I,J
      INTEGER NRAY
      INTEGER IOS
      INTEGER LUER, ERROR_LU
      LOGICAL NEW_GRID
      EXTERNAL ERROR_LU
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	INTEGER DI
	INTEGER ANG_MULT
	INTEGER NA
!
!=======================================================================
!	V_EXT=0.0D0
!
! If necessary, calculate characteristic rays for moving atmospheres
!
!------------------------------------------------------------------------
!
      LUER=ERROR_LU()
      IF(ALLOCATED(R_EXT_SAV) .AND. ND_EXT .EQ. ND_EXT_SAV)THEN
	NEW_GRID=.FALSE.
	DO I=1,ND_EXT
	  IF(R_EXT(I) .NE. R_EXT_SAV(I))THEN
	    NEW_GRID=.TRUE.
	    EXIT
	  END IF
	END DO
      ELSE
	NEW_GRID=.TRUE.
      END IF
!
      IF(NEW_GRID)THEN
        WRITE(LUER,*)'   Calculating characteristic rays.....'
        WRITE(LUER,*)'   From CHARACTERISTICS_V2_GAM.....'
!	I=ND_EXT/2
!	WRITE(LUER,'(A7,I3,1X,3(A5,I3,1X))')'ND_EXT:',ND_EXT,'ND:',ND,'NC:',NC,'NP:',NP
!	WRITE(LUER,'(A,1X,F6.2,2X,A9,1X,3(ES12.3))')'VDOP_FRAC:',VDOP_FRAC,'VDOP_VEC:',
!	1	VDOP_VEC(1),VDOP_VEC(I),VDOP_VEC(ND_EXT)
        CALL CHARACTERISTICS_V2_GAM(R_EXT,P,V_EXT,VDOP_VEC,VDOP_FRAC,ND_EXT,NC,NP)
	WRITE(LUER,*)'   Exited CHARACTERISTICS_V2_GAM.....'
      END IF
!--------------------------------------------------------------------
!
      IF(NEW_GRID .AND. ALLOCATED(TAU))THEN
	DEALLOCATE (TAU, DTAU)
	DO IP=1,NP
          DEALLOCATE (RAY(IP)%I_P, RAY(IP)%I_M, RAY(IP)%LNK)
          DEALLOCATE (RAY(IP)%I_P_PREV, RAY(IP)%I_M_PREV)
          DEALLOCATE (RAY(IP)%I_P_SAVE, RAY(IP)%I_M_SAVE)
	END DO
      END IF
!
      IF(NEW_GRID)THEN
        J=0
	DO IP=1,NP
          NRAY=RAY(IP)%NZ
          J=MAX(J,NRAY)
	  ALLOCATE (RAY(IP)%I_P(NRAY))
          ALLOCATE (RAY(IP)%I_M(NRAY))
	  ALLOCATE (RAY(IP)%I_P_PREV(NRAY))
          ALLOCATE (RAY(IP)%I_M_PREV(NRAY))
	  ALLOCATE (RAY(IP)%I_P_SAVE(NRAY))
          ALLOCATE (RAY(IP)%I_M_SAVE(NRAY))
          ALLOCATE (RAY(IP)%LNK(ND)); RAY(IP)%LNK(:)=0
	  ALLOCATE (RAY(IP)%ETA_M(NRAY))
	  ALLOCATE (RAY(IP)%ETA_P(NRAY))
        END DO
        ALLOCATE (TAU(J))
        ALLOCATE (DTAU(J))
      END IF
!
! Determine the links betweem the extended-fine grid and the original R grid.
!
	DO IP=1,NP
	  RAY(IP)%LNK(1:ND)=0
	  ID=1
	  DO J=1,RAY(IP)%NZ
	    IF(RAY(IP)%R_RAY(J) .EQ. R(ID))THEN
	      RAY(IP)%LNK(ID)=J
	      WRITE(160,*)IP,J,RAY(IP)%LNK(ID)
	      ID=ID+1
	    END IF
	  END DO
	END DO
!
! Check all lnks allocated
!
	DO IP=1,NP
	  DO ID=1,MIN(ND,ND-(IP-NC-1))
	    IF(RAY(IP)%LNK(ID) .EQ. 0)THEN
	      WRITE(LUER,*)'Error in DEFNE_GRID_V2'
	      WRITE(LUER,*)'LNK not defined: IP=',IP,'ID=',ID
	      STOP
	    END IF
	  END DO
	END DO
!
	IF( ALLOCATED(JQW_P) .AND. (ND .NE. ND_SAV .OR. NP .NE. NP_SAV) )THEN
          DEALLOCATE(JQW_P, HQW_P, KQW_P, NQW_P)
          DEALLOCATE(JQW_M, HQW_M, KQW_M, NQW_M)
          DEALLOCATE(I_P_GRID, I_M_GRID)
	END IF
!
	IF(.NOT. ALLOCATED(JQW_P))THEN
	  ALLOCATE (JQW_P(ND,NP),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (HQW_P(ND,NP),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (KQW_P(ND,NP),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (NQW_P(ND,NP),STAT=IOS)
!
          IF(IOS .EQ. 0)ALLOCATE (JQW_M(ND,NP),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (HQW_M(ND,NP),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (KQW_M(ND,NP),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (NQW_M(ND,NP),STAT=IOS)
!
          IF(IOS .EQ. 0)ALLOCATE (I_P_GRID(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (I_M_GRID(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in DEFINE_GRID_V2_GAM'
	    WRITE(LUER,*)'Unable to allocate JQW_P etc'
	    WRITE(LUER,*)'Status =',IOS
	    STOP
	  END IF
	END IF
!
! Calculate J, H, K and N integration weights for mu+ and mu-
!
	DO ID=1,ND
	  J=NP+1-ID
	  DO IP=1,J
            MU(IP)=RAY(IP)%MU_P(RAY(IP)%LNK(ID))
	  END DO
	  CALL J_WEIGHT(J,MU,TEMP_VEC); JQW_P(ID,1:J)=TEMP_VEC(1:J)
	  CALL H_WEIGHT(J,MU,TEMP_VEC); HQW_P(ID,1:J)=TEMP_VEC(1:J)
	  CALL K_WEIGHT(J,MU,TEMP_VEC); KQW_P(ID,1:J)=TEMP_VEC(1:J)
	  CALL N_WEIGHT(J,MU,TEMP_VEC); NQW_P(ID,1:J)=TEMP_VEC(1:J)
	END DO
!
	DO ID=1,ND
	  J=NP+1-ID
	  DO IP=1,J
            MU(IP)=RAY(IP)%MU_M(RAY(IP)%LNK(ID))
	  END DO
	  CALL J_WEIGHT(J,MU,TEMP_VEC); JQW_M(ID,1:J)=TEMP_VEC(1:J)
	  CALL H_WEIGHT(J,MU,TEMP_VEC); HQW_M(ID,1:J)=TEMP_VEC(1:J)
	  CALL K_WEIGHT(J,MU,TEMP_VEC); KQW_M(ID,1:J)=TEMP_VEC(1:J)
	  CALL N_WEIGHT(J,MU,TEMP_VEC); NQW_M(ID,1:J)=TEMP_VEC(1:J)
	END DO
!
! Quadrature weights for mu_m are negative.  Change them to positive.
!
        DO IP=1,NP
          DO J=1,ND
            JQW_M(J,IP)=-JQW_M(J,IP)
            HQW_M(J,IP)=-HQW_M(J,IP)
            KQW_M(J,IP)=-KQW_M(J,IP)
            NQW_M(J,IP)=-NQW_M(J,IP)
          END DO
        END DO
!
        IF(NEW_GRID)THEN
          IF(ALLOCATED(R_EXT_SAV))DEALLOCATE(R_EXT_SAV)
          ALLOCATE(R_EXT_SAV(ND_EXT))
	  R_EXT_SAV(:)=R_EXT
	  ND_EXT_SAV=ND_EXT
          ND_SAV=ND
          NP_SAV=NP
	END IF
!
! saving the inverted array for the gamma-ray transfer routines
!
	IF(.NOT. ALLOCATED(R_MU) .AND. VERBOSE_GAMMA)THEN
	  ALLOCATE(R_MU(ND))
	  WRITE(LUER,*)'ALLOCATED R_MU'
	  WRITE(LUER,*)'Size of R_MU:',SIZE(R_MU)
	END IF
	DO ID=1,ND
	  J=2*(NP-ID)+1
	  ALLOCATE(R_MU(ID)%MU_VECTOR(J))
	  R_MU(ID)%MU_PTS=J
	  R_MU(ID)%MU_VECTOR=0.0_LDP
	  J=NP-ID+1
	  DO IP=1,J
	    R_MU(ID)%MU_VECTOR(IP)=RAY(IP)%MU_P(ID)
	  END DO
	  I=NP-ID+1
	  DO IP=J-1,1,-1 !this will step down since MU_M starts at -1
	    I=I+1
	    R_MU(ID)%MU_VECTOR(I)=RAY(IP)%MU_M(ID)
	  END DO
	END DO
!
	IF(VERBOSE_GAMMA)THEN
	  OPEN(UNIT=7,FILE='./data/GAM_MU_GRID',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  IF(IOS .NE. 0) STOP '****Cannot open GAM_MU_GRID****'
	  WRITE(7,'(A,I4)')'Num Depth pts:',ND
	  CALL SET_LINE_BUFFERING(7)
	  DO ID=1,ND
	    NA=R_MU(ID)%MU_PTS
	    WRITE(7,'(A11,1X,I3,4X,A10,1X,I3)')'# Angles:',NA,'Should be:',2*(NP-ID+1)-1
	    WRITE(7,'(A6,1X,I3)')'Depth:',ID
	    DO I=1,NA
	      WRITE(7,'(I3,1X,F14.8)')I,R_MU(ID)%MU_VECTOR(I)
	    END DO
	    WRITE(7,*) ' '
	  END DO
	  CLOSE(UNIT=7)
	END IF
!
	DO ID=1,ND
	  IF(.NOT. ALLOCATED(R_MU(ID)%MON_MU))THEN
	    R_MU(ID)%MON_MU_PTS=R_MU(ID)%MU_PTS*ANG_MULT
	    NA=R_MU(ID)%MON_MU_PTS
	    ALLOCATE(R_MU(ID)%MON_MU(NA))
	    R_MU(ID)%dMU=2.0_LDP/(NA-1)
	    R_MU(ID)%MON_MU(1)=1.0_LDP
	    DO I=2,NA-1
	      R_MU(ID)%MON_MU(I)=1.0_LDP-(I-1)*R_MU(ID)%dMU
	    END DO
	    R_MU(ID)%MON_MU(NA)=-1.0_LDP
	  END IF
	END DO
!
      RETURN
      END
