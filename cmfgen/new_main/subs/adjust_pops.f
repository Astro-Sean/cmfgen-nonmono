!
	SUBROUTINE ADJUST_POPS(POPS,LUER,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE OLD_GRID_MODULE
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered 04-Dec-2023: Moved initialization of Z_POP so it was still being set.
! Altered 25-Nov-2023: Now initialize Z_POP (NT and NT-1 wer not set), and added LU_VB.
! Altered 27-Mar-2023: Added VERBOSE_OUTPUT to most of the output options.
!                        File 175 now has name (HYDRO_GRID_SUMMARIES).
! Altered 04-Jan-2020: Changed LTEPOP_WLD and CNVT_FR_DC to V2.
! Altered 20-Jun-2008: Z_POP now set (was not set)a
!                      Improved handling at boundaries.
! Created
!
	INTEGER ND
	INTEGER NT
	INTEGER LUER				!Unit for error messages
	INTEGER LU_GRID
	INTEGER LU_VB
!
	REAL(KIND=LDP) LOG_OLD_TAU(OLD_ND)
	REAL(KIND=LDP) TB(OLD_ND)
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) Z_POP(NT)			!Vector containing Z of atom/ion (not core)
!
	REAL(KIND=LDP) TAU(ND)
	REAL(KIND=LDP) LOG_TAU(ND)
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) T1
!
	INTEGER NXST,NX
	INTEGER I,J,L
	INTEGER ISPEC
	INTEGER ID
	INTEGER N_INT
	LOGICAL FIRST
!
! We assume a power law density distribution at the outer boundary.
!
	TB(1:OLD_ND)=OLD_ROSS_MEAN(1:OLD_ND)*OLD_CLUMP_FAC(1:OLD_ND)
	T1=LOG(TB(5)/TB(1))/LOG(OLD_R(1)/OLD_R(5))
	OLD_TAU(1)=TB(1)*OLD_R(1)/T1
	DO I=2,OLD_ND
	  OLD_TAU(I)=OLD_TAU(I-1)+0.5_LDP*(OLD_R(I-1)-OLD_R(I))*(TB(I-1)+TB(I))
	END DO
!
	VERBOSE_OUTPUT=.FALSE.
	IF(VERBOSE_OUTPUT)THEN
	  CALL GET_LU(LU_VB,'Verbose information in ADJUST_POPS')
	  OPEN(UNIT=LU_VB,FILE='RGRID_INFO',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
	  WRITE(LU_VB,*)ND
	  J=5; CALL WRITV_V2(OLD_TAU,ND,J,'OLD_TAU',LU_VB)
	  FLUSH(UNIT=LU_VB)
	END IF
!
	TA(1:ND)=ROSS_MEAN(1:ND)			!*CLUMP_FAC(1:ND)
	T1=LOG(TA(5)/TA(1))/LOG(R(1)/R(5))
	IF(T1 .LT. 2.0_LDP)THEN
	  IF(R(1) .GT. 4*R(ND))THEN
	    TAU(1)=TA(1)*R(1)
	  ELSE
	    TAU(1)=TA(1)*(R(1)-R(ND))/8.0_LDP
	  END IF
	ELSE
	  TAU(1)=TA(1)*R(1)/T1
	END IF
	DO I=2,ND
	  TAU(I)=TAU(I-1)+0.5_LDP*(R(I-1)-R(I))*(TA(I-1)+TA(I))
	END DO
!
	IF(VERBOSE_OUTPUT)THEN
	  WRITE(LU_VB,*)ND
	  CALL WRITV_V2(R,ND,J,'R',LU_VB)
	  CALL WRITV_V2(ROSS_MEAN,ND,J,'ROSS_MEAN',LU_VB)
	  CALL WRITV_V2(TAU,ND,J,'TAU',LU_VB)
	  FLUSH(UNIT=LU_VB)
	END IF
!
! Find indices for new tau grid contained in old tau grid.
!
	NXST=1
	DO WHILE(TAU(NXST) .LT. OLD_TAU(1))
	  NXST=NXST+1
	END DO
	NX=ND
	DO WHILE(TAU(NX) .GT. OLD_TAU(OLD_ND))
	  NX=NX-1
	END DO
	LOG_TAU=LOG(TAU)
	LOG_OLD_TAU=LOG(OLD_TAU)
	N_INT=NX-NXST+1
!
	TB(1:OLD_ND)=LOG(OLD_T(1:OLD_ND))
	CALL LINPOP(LOG_TAU(NXST),TA(NXST),N_INT,LOG_OLD_TAU,TB,OLD_ND)
	T(NXST:NX)=EXP(TA(NXST:NX))
	T(1:NXST-1)=OLD_T(1)
	DO I=NX+1,ND
	  T(I)=OLD_T(OLD_ND)*(TAU(I)+0.67_LDP)/(OLD_TAU(OLD_ND)+0.67_LDP)
	END DO
!
	IF(VERBOSE_OUTPUT)THEN
	  WRITE(LU_VB,*)ND
	  J=5;  CALL WRITV_V2(T,ND,J,'T',LU_VB)
	  J=10; CALL WRITV_V2(CLUMP_FAC,ND,J,'CLUMP_FAC',LU_VB)
	  FLUSH(UNIT=LU_VB)
	END IF
!
	ED=ED/CLUMP_FAC
!
	IF(VERBOSE_OUTPUT)THEN
	  CALL GET_LU(LU_GRID,'Grid summaries in ADJUST_POPS')
	  OPEN(UNIT=LU_GRID,FILE='HYDRO_GRID_SUMMARIES',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
	    WRITE(LU_GRID,'(4X,10(5X,A))')
	1       '     OLD_R','   OLD_T','  OLD_ED',' OLD_TAU','OLD_ROSS',
	1       '     NEW_R','   NEW_T','  NEW_ED',' NEW_TAU','NEW_ROSS'
	    DO I=1,ND
	      WRITE(LU_GRID,'(I4,2(ES15.6,4ES13.3))')I,
	1                OLD_R(I),OLD_T(I),OLD_ED(I),OLD_TAU(I),OLD_ROSS_MEAN(I),
	1                R(I),T(I),ED(I),TAU(I),ROSS_MEAN(I)/CLUMP_FAC(I)
	    END DO
	    FLUSH(UNIT=LU_GRID)
	  CLOSE(LU_GRID)
	END IF
!
! Loop over all species and ionization stages
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    IF(VERBOSE_OUTPUT)WRITE(6,*)'In loop',ISPEC,ID
!
! Compute departure coefficients
!
	    ATM(ID)%XzV_F=LOG(ATM(ID)%XzV_F)-ATM(ID)%LOG_XzVLTE_F
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LU_VB,'(/,1X,A)')'Ion info'
	      WRITE(LU_VB,'(A,3X,3ES14.4)')TRIM(ION_ID(ID)),ATM(ID)%XzV_F(1,1),
	1                  ATM(ID)%XzVLTE_F(1,1),ATM(ID)%DXzV_F(1)
	      FLUSH(UNIT=LU_VB)
	    END IF
!
! Put the departure coefficients on the new grid.
!
	    N_INT=NX-NXST+1
	    DO I=1,ATM(ID)%NXZV_F
	       TB(1:OLD_ND)=ATM(ID)%XzV_F(I,1:OLD_ND)
	       CALL LINPOP(LOG_TAU(NXST),TA(NXST),N_INT,LOG_OLD_TAU,TB,OLD_ND)
	       ATM(ID)%XzV_F(I,NXST:NX)=TA(NXST:NX)
	    END DO
!
! Now do the boudaries which may extend outside the old grid.
!
	    J=ATM(ID)%NXzV_F
	    DO L=2,NXST-1
	      ATM(ID)%XzV_F(1:J,L)=ATM(ID)%XzV_F(1:J,1)
	    END DO
	    DO L=NX+1,ND-1
	      ATM(ID)%XzV_F(1:J,L)=ATM(ID)%XzV_F(1:J,ND)
	    END DO
!
	  END DO	!Ion loop
!
! Interplate the ion density.
!
	  IF(SPECIES_PRES(ISPEC))THEN
	    ID=SPECIES_END_ID(ISPEC)-1
	    TB(1:OLD_ND)=LOG(ATM(ID)%DXzV_F(1:OLD_ND))
	    CALL LINPOP(LOG_TAU(NXST),TA(NXST),N_INT,LOG_OLD_TAU,TB,OLD_ND)
	    ATM(ID)%DXzV_F(NXST:NX)=EXP(TA(NXST:NX))
	    T1=ATM(ID)%DXzV_F(1)/OLD_MASS_DENSITY(1)
	    ATM(ID)%DXzV_F(1:NXST-1)=T1*DENSITY(1:NXST-1)
	    T1=ATM(ID)%DXzV_F(OLD_ND)/OLD_MASS_DENSITY(OLD_ND)
	    ATM(ID)%DXzV_F(NX+1:ND)=T1*DENSITY(NX+1:ND)
	  END IF
!
	END DO		!Species loop
!
!
! Compute vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD.
!
	TB(1:OLD_ND)=OLD_POPION/OLD_POP_ATOM
	CALL LINPOP(LOG_TAU(NXST),TA(NXST),N_INT,LOG_OLD_TAU,TB,OLD_ND)
	TA(1:NXST-1)=TB(1)
	TA(NX+1:ND)=TB(OLD_ND)
	POPION(1:ND)=TA(1:ND)*POP_ATOM(1:ND)
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! Compute the LTE populations and convert from departure coefficients to
! populations. The populations are scaled to ensure number conservation -
! during the scaling the ionization remains fixed. On each call to CNVT_FR_DC,
! TA is incremented by the population. IF FIRST is .TRUE., TA is zeroed first.
! The second flag indicates whether to add in the ION contribution, which
! will also be added as the ground state population of the next species.
! If it is FALSE, the higher ionization stage is assumed not to be present,
! and DION is added in.
!
! TB is used as a dummy vector when we are dealing with the lowest ionization
! stage. It is returned with the ground state population.
!
! CNVT_FR_DC_V2 requires the we pass LOG(DCs) to it.
!
	WRITE(6,*)'Beginning DC interpolation in ADJUST_POPS'; FLUSH(UNIT=6)
	DO ISPEC=1,NUM_SPECIES
	  FIRST=.TRUE.
	  DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL LTEPOP_WLD_V2(ATM(ID)%XzVLTE_F, ATM(ID)%LOG_XzVLTE_F, ATM(ID)%W_XzV_F,
	1              ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,
	1              ATM(ID)%ZXzV,      ATM(ID)%GIONXzV_F,
	1              ATM(ID)%NXzV_F,    ATM(ID)%DXzV_F,     ED,T,ND)
!	      ATM(ID)%XzV_F=LOG(ATM(ID)%XzV_F)
	      IF(VERBOSE_OUTPUT)THEN
	        WRITE(LU_VB,'(/,1X,A)')'Ion info'
	        WRITE(LU_VB,'(A,3X,3ES14.4)')TRIM(ION_ID(ID)),ATM(ID)%XzV_F(1,1),
	1                                ATM(ID)%LOG_XzVLTE_F(1,1),ATM(ID)%DXzV_F(1)
	        FLUSH(LU_VB)
	      END IF
	      CALL CNVT_FR_DC_V2(ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
	1              ATM(ID)%DXzV_F,    ATM(ID)%NXzV_F,
	1              TB,                TA,ND,
	1              FIRST,             ATM(ID+1)%XzV_PRES)
!	    CALL WRITEDC_V3( ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
!	1          ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,IONE,
!	1          R,T,ED,V,CLUMP_FAC,LUM,ND,
!	1          TRIM(ION_ID(ID))//'BFM','DC',IONE)
	      IF(ID .NE. SPECIES_BEG_ID(ISPEC))ATM(ID-1)%DXzV_F(1:ND)=TB(1:ND)
	    END IF
	  END DO
!
! Scale the populatons to match the actual density.
!
	  IF(VERBOSE_OUTPUT)THEN
	    J=5; CALL WRITV_V2(TA,ND,J,SPECIES(ISPEC),LU_VB); FLUSH(UNIT=LU_VB)
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      CALL WRITEDC_V3( ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
	1            ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,IONE,
	1            R,T,ED,V,CLUMP_FAC,LUM,ND,
	1            TRIM(ION_ID(ID))//'OUTDC','DC',IONE)
	      CALL WRITEDC_V3( ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
	1            ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,IONE,
	1            R,T,ED,V,CLUMP_FAC,LUM,ND,
	1            TRIM(ION_ID(ID))//'POPS','POP',IONE)
	    END DO
	  END IF
!
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    CALL SCALE_POPS(ATM(ID)%XzV_F,ATM(ID)%DXzV_F,
	1              POP_SPECIES(1,ISPEC),TA,ATM(ID)%NXzV_F,ND)
	  END DO
	END DO		!Loop over species
	IF(VERBOSE_OUTPUT)CLOSE(LU_VB)
!
! We now need to compute the populations for the model atom with Super-levels.
! We do this in reverse order (i.e. highest ionization stage first) in order
! that we the ion density for the lower ionization stage is available for
! the next call.
!
! For 1st call to FULL_TO_SUP, Last line contains FeX etc as FeXI not installed.
!
	DO ID=NUM_IONS-1,1,-1
	  CALL FULL_TO_SUP(
	1      ATM(ID)%XzV,   ATM(ID)%NXzV,       ATM(ID)%DXzV,   ATM(ID)%XzV_PRES,
	1      ATM(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,
	1      ATM(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES,  ND)
	END DO
!
! Store all quantities in POPS array. This is done here as it enables POPION
! to be readily computed. It also ensures that POS is correct if we don't
! iterate on T.
!
	DO ID=1,NUM_IONS-1
	  CALL IONTOPOP(POPS,  ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1          ATM(ID)%EQXzV, ATM(ID)%NXzV, NT,ND, ATM(ID)%XzV_PRES)
	END DO
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	Z_POP=0.0_LDP
        DO ID=1,NUM_IONS-1
	  CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1              ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
        END DO
!
	DO J=1,ND
	  POPION(J)=0.0_LDP
	  DO I=1,NT
	    IF(Z_POP(I) .GT. 0.01_LDP)POPION(J)=POPION(J)+POPS(I,J)
	  END DO
	END DO
!
! Evaluates LTE populations for both the FULL atom, and super levels.
!
	CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
!
	RETURN
	END
