	SUBROUTINE DET_LTE_ED(TMIN,ND,DO_LEV_DIS)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER ND
	LOGICAL DO_LEV_DIS
!
	REAL(KIND=LDP) ION_POPS(ND,NUM_IONS)
	REAL(KIND=LDP) FSAHA(NUM_IONS)
	REAL(KIND=LDP) XZ(NUM_IONS)
	REAL(KIND=LDP) XZW(NUM_IONS)
	REAL(KIND=LDP) ED_OLD(ND)
	REAL(KIND=LDP) POPION_OLD(ND)
	REAL(KIND=LDP) TEMP(ND)
	REAL(KIND=LDP) TMIN
!
	REAL(KIND=LDP) XEDW
	REAL(KIND=LDP) XED_OLD
	REAL(KIND=LDP) XERR
	REAL(KIND=LDP) XKBT
	REAL(KIND=LDP) XG0
	REAL(KIND=LDP) XG1
	REAL(KIND=LDP) XGE
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) HSQR
	REAL(KIND=LDP), PARAMETER :: HDKT=4.7994145_LDP
	REAL(KIND=LDP), PARAMETER :: EPS=1.0E-05_LDP
!
	REAL(KIND=LDP) PLANCKS_CONSTANT
	REAL(KIND=LDP) BOLTZMANN_CONSTANT
	REAL(KIND=LDP) ELECTRON_MASS
	EXTERNAL PLANCKS_CONSTANT
	EXTERNAL BOLTZMANN_CONSTANT
	EXTERNAL ELECTRON_MASS
!
	INTEGER ISTART,IEND
	INTEGER L
	INTEGER K
	INTEGER I
	INTEGER IW
	INTEGER ID
	INTEGER ISPEC
	LOGICAL CONVERGED
!
	PI=ACOS(-1.0_LDP)
        HSQR = PLANCKS_CONSTANT()*PLANCKS_CONSTANT()
!
! Perform initializations.
!
	TEMP(1:ND)=T(1:ND)
	DO I=1,ND
	  IF(T(I) .LT. TMIN)TEMP(I)=T(I)
	END DO
	ED_OLD=0.0_LDP
	ED=POP_ATOM
	POPION=POP_ATOM
	POPION_OLD=0.0_LDP
	FSAHA=0.0_LDP
	CONVERGED=.FALSE.
!
! Compute the effective statistical weight for all levels in all ions.
! The effective statistical weight is the product of the actual statistical
! weight, and the Boltzmann excitation factor. We store it in the FULL
! level population vector.
!
	DO L=1,ND
          XKBT = BOLTZMANN_CONSTANT() * TEMP(L) * 1.0E4_LDP
          T1 = (2.0_LDP*PI*ELECTRON_MASS()*XKBT/HSQR)**1.5_LDP
          T2 = HDKT/TEMP(L)
          DO ID=1,NUM_IONS
            IF (ATM(ID)%XzV_PRES) THEN
              ATM(ID)%XzV_F(1,L)= ATM(ID)%GXzV_F(1)
	      DO I=2,ATM(ID)%NXzV_F
	        ATM(ID)%XzV_F(I,L)=ATM(ID)%GXzV_F(I)*EXP(T2*(ATM(ID)%EDGEXzV_F(I)-ATM(ID)%EDGEXzV_F(1)))
	      END DO
            END IF
          END DO
	END DO
!
	DO WHILE(.NOT. CONVERGED)
	   CONVERGED=.TRUE.
!
! Compute the occupation probabilities, first updating the dissolution constants.
!
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DIS,ND)
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL OCCUPATION_PROB(ATM(ID)%W_XZV_F,ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%ZXzV,ATM(ID)%NXzV_F,ND)
	    END IF
	  END DO
!
! Solve the Saha equation for the mixture considered: we obtain the electron density
! and the relative population of all ions.
!
	  DO L=1,ND
	    DO WHILE( ABS(1.0_LDP-ED_OLD(L)/ED(L)) .GT. EPS .OR.
	1             ABS(1.0_LDP-POPION_OLD(L)/POPION(L)) .GT. EPS )
	      WRITE(266,*)L,ED(L)
	      ED_OLD(L)=ED(L)
	      POPION_OLD(L)=POPION(L)
	      CONVERGED=.FALSE.
!
              XKBT = BOLTZMANN_CONSTANT() * TEMP(L) * 1.0E4_LDP
              T1 = (2.0_LDP*PI*ELECTRON_MASS()*XKBT/HSQR)**1.5_LDP
              T2 = HDKT/TEMP(L)
              DO ID=1,NUM_IONS
                IF (ATM(ID)%XzV_PRES) THEN
                  XGE = 2.0_LDP
                  XG0 = ATM(ID)%GXzV_F(1)
	          DO I=2,ATM(ID)%NXzV_F
	            XG0=XG0+ATM(ID)%W_XzV_F(I,L)*ATM(ID)%XzV_F(I,L)
	          END DO
                  XG1 = ATM(ID)%GIONXzV_F
                  IF(ATM(ID+1)%XzV_PRES) THEN
	            DO I=2,ATM(ID+1)%NXzV_F
	              XG1=XG1+ATM(ID+1)%W_XzV_F(I,L)*ATM(ID+1)%XzV_F(I,L)
                    END DO
	          END IF
	          FSAHA(ID) = (XG1*XGE/XG0) * T1 * EXP(-ATM(ID)%EDGEXzV_F(1)*T2)
                END IF
              END DO	!ion loop.
!
	      XEDW = ED(L)
              XED_OLD = XEDW
              XERR = 2.0_LDP * EPS
	      DO WHILE (XERR .GT. EPS)
                T3 = 0.0_LDP
                DO ISPEC=1,NUM_SPECIES
                  ISTART = SPECIES_BEG_ID(ISPEC)
                  IF (ISTART.NE.0) THEN
                    IEND = SPECIES_END_ID(ISPEC)
                    T1 = 1.0_LDP
                    T2 = 1.0_LDP
                    DO IW=ISTART+1,IEND
                      T1 = (FSAHA(IW-1)/XEDW) * T1
                      T2 = T2 + T1
                      XZW(IW) = T1
                    END DO
                    XZ(ISTART) = 1.0_LDP / T2
                    DO IW=ISTART+1,IEND
                      XZ(IW) = XZW(IW) * XZ(ISTART)
                      T3 = T3 + XZ(IW) * POP_SPECIES(L,ISPEC) * ATM(IW-1)%ZXzV
                    END DO
                    XZ(ISTART:IEND) = XZ(ISTART:IEND) * POP_SPECIES(L,ISPEC)
                  END IF
                END DO
                XED_OLD = XEDW
                XEDW = T3
                XERR = ABS(1.0_LDP-XED_OLD/XEDW)
	        XEDW=0.5_LDP*(XEDW+XED_OLD)
	        WRITE(6,*)T(L),XED_OLD,XEDW
              END DO
              ED(L)=XEDW
              ION_POPS(L,:) = XZ(:)  ! each entry corresponds to ions ID=1,NUM_IONS
!
	    END DO	!Only do depths that are not finished
	  END DO	!Loop over L index
!
! Get total ion population for level dissolution. Note that ATM(ID)$ZxzV
! is the charge on the ion AFTER the valence electron is removed
! (i.e., it =1 for HI).
!
	  POPION=0.0_LDP
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      IF(ATM(ID)%ZXzV .GT. 1.01_LDP)THEN
	        DO L=1,ND
	          POPION(L)=POPION(L)+ION_POPS(L,ID)
	        END DO
	      END IF
	    END DO
	    IF(SPECIES_BEG_ID(ISPEC) .NE. 0)THEN
	      DO L=1,ND
	        POPION(L)=POPION(L)+ION_POPS(L,SPECIES_END_ID(ISPEC))
	      END DO
	    END IF
	  END DO
!
	END DO		!All Ne are not accurate
!
! Compute the occupation probabilities, first updating the dissolution constants
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DIS,ND)
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL OCCUPATION_PROB(ATM(ID)%W_XZV_F,ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%ZXzV,ATM(ID)%NXzV_F,ND)
	  END IF
	END DO
!
	RETURN
	END
