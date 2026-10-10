!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
	SUBROUTINE ADJUST_R_GRID(POPS,P,
	1             JQW,HQW,KQW,HMIDQW,NMIDQW,
	1             MU_AT_RMAX,HQW_AT_RMAX,TRAPFORJ,
	1             FLUXMEAN,ESEC,GRID_TYPE,RG_PARS,N_PARS,
	1             ND,NT,NC,NP)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered : 18-Mar-2004 FLUXMEAN passed insted of dTAU_OLD
!                       ESEC passed in call.
!                       As still developing, subroutine name NOT changed.
! Created : 25-Feb-2004
!
	INTEGER ND,NT,NC,NP
!
	REAL(KIND=LDP) POPS(NT,ND)
!
! NB: dTAU_OLD(I) = Optical depth increment between depth I & I+1.
!                   It is NOT the optical depth scale.
!                   Generally computed using the FLUX_MEAN opacity.
!
	LOGICAL TRAPFORJ
	REAL(KIND=LDP) P(NP)
	REAL(KIND=LDP) JQW(ND,NP)
	REAL(KIND=LDP) HQW(ND,NP)
	REAL(KIND=LDP) KQW(ND,NP)
	REAL(KIND=LDP) HMIDQW(ND,NP)
	REAL(KIND=LDP) NMIDQW(ND,NP)
	REAL(KIND=LDP) MU_AT_RMAX(NP)
	REAL(KIND=LDP) HQW_AT_RMAX(NP)
!
	REAL(KIND=LDP) FLUXMEAN(ND)
	REAL(KIND=LDP) ESEC(ND)
!
	INTEGER N_PARS
	CHARACTER(LEN=*) GRID_TYPE
	REAL(KIND=LDP) RG_PARS(N_PARS)
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	EXTERNAL JWEIGHT,HWEIGHT,KWEIGHT,NWEIGHT
        EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
!
! Local variables.
!
	REAL(KIND=LDP) dTAU_OLD(ND)
	REAL(KIND=LDP) R_OLD(ND)
	REAL(KIND=LDP) LOG_R_OLD(ND)
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP) TAU_OLD(ND)
	REAL(KIND=LDP) TAU(ND)
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
!
! The fine grid (FG) is chosen to cover the ioization front. The default valuse are
! -3.0 to 1.5D0 in log(TAU) space.
!
	REAL(KIND=LDP) FG_MIN                   !Default=-3.0D0
	REAL(KIND=LDP) FG_MAX			!Default=1.5D0
	REAL(KIND=LDP) FG_RANGE			!Default=4.5D0
!
	REAL(KIND=LDP) W1(NP),W2(NP),W3(NP)
!
	REAL(KIND=LDP) DTAU
	REAL(KIND=LDP) T1
	INTEGER, PARAMETER :: IONE=1
	INTEGER LS,I,I1,I2
	INTEGER NX
	INTEGER ISPEC
	LOGICAL MID
!
	WRITE(139,*)'Entering ADJUST_R_GRID'
!
! Compute opitcal depth scale.
!
	DO I=1,ND
	  TA(I)=FLUXMEAN(I)
	  IF(TA(I) .LE. 0.0_LDP)TA(I)=ESEC(I)
	  TA(I)=TA(I)*CLUMP_FAC(I)
	END DO
        TB(1:ND)=0.0_LDP                              !Used for dCHIdR
        CALL NORDTAU(dTAU_OLD,TA,R,R,TB,ND)
!
! Save existing grid, which will be used for the interplations.
!
	R_OLD(1:ND)=R(1:ND)
!
! Compute optical depth scale. Note that we are passed the optical detph
! increments, not the optical depth scale.
!
	TAU_OLD(1:ND)=0.0_LDP
	DO I=2,ND
	  TAU_OLD(I)=TAU_OLD(I-1)+dTAU_OLD(I-1)
	END DO
!
! Revise TAU near the outer boundary to allow for the fact that TAU_OLD(1)
! is zero. These expressions are correct for a power law opacity.
!
	TAU_OLD(1)=TAU_OLD(3)*R_OLD(3)/R_OLD(1)
	TAU_OLD(2)=TAU_OLD(3)*R_OLD(2)/R_OLD(1)
	TAU_OLD(1:ND)=LOG10(TAU_OLD(1:ND))
!
	FG_MIN=-3.0_LDP
	FG_MAX=1.50_LDP
	FG_RANGE=FG_MAX-FG_MIN
!
	IF( TRIM(GRID_TYPE) .EQ. 'UNIFORM')THEN
	  DTAU=(TAU_OLD(ND)-TAU_OLD(1))/(ND-5)
	  IF(N_PARS .EQ. 2)THEN
	    FG_MIN=RG_PARS(1)
	    FG_MAX=RG_PARS(2)
	    FG_RANGE=FG_MAX-FG_MIN
	  END IF
	  NX=(FG_RANGE+0.7_LDP*DTAU)/DTAU
	ELSE IF( TRIM(GRID_TYPE) .EQ. 'FIX_NX')THEN
	  NX=RG_PARS(1)
	  IF(N_PARS .EQ. 3)THEN
	    FG_MIN=RG_PARS(2)
	    FG_MAX=RG_PARS(3)
	    FG_RANGE=FG_MAX-FG_MIN
	  END IF
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid GRID_TYPE in ADJUST_R_GRID'
	  WRITE(LUER,*)'GRID_TYPE=',TRIM(GRID_TYPE)
	  STOP
	END IF
!
! Compute the new grid. In the present version, NX points are
! inserted between log(TAU)=-3.0 and 1.5. We also have 2 special
! points at either boundary, which are placed separately.
!
	T1=TAU_OLD(ND)-TAU_OLD(1)-FG_RANGE
        DTAU=T1/(ND-4-NX)
	I1=(FG_MIN-TAU_OLD(1))/DTAU
	DTAU=(FG_MIN-TAU_OLD(1))/I1
        I1=I1+2
	TAU(1)=TAU_OLD(1)
	DO I=4,I1
	  TAU(I)=TAU_OLD(1)+DTAU*(I-3)
	END DO
        TAU(2)=MIN(TAU_OLD(2),TAU(1)+0.1_LDP*DTAU)
        TAU(3)=MIN(TAU_OLD(3),TAU(1)+0.3_LDP*DTAU)
!
! Do the insertion in the crtical section.
!
	DTAU=4.5_LDP/(NX-1)
	DO I=I1+1,I1+NX
	  TAU(I)=FG_MIN+DTAU*(I-I1-1)
	END DO
!
! Now do the last section, towards the inner boundary.
!
	I2=ND-(NX+I1+2)
	T1=(TAU_OLD(ND)-FG_MAX)
	DTAU=T1/I2
	DO I=I1+NX+1,ND-3
	  TAU(I)=TAU(I-1)+DTAU
	END DO
	TAU(ND)=TAU_OLD(ND)
	TAU(ND-1)=MAX(TAU_OLD(ND-1),TAU(ND)-0.1E00_LDP*DTAU)
	TAU(ND-2)=TAU(ND)-0.5E00_LDP*DTAU
!
! Compute the new radius grid. Linear interpolation is
! more than adequate, since we're just defining a new grid.
!
	CALL LININT(TAU,R,ND,TAU_OLD,R_OLD,ND)
!
! We now need to regrid all the populations. All interpolations (except
! sigma) are performed in the LOG-LOG plane. For SN this is ideal, since
! the density and velocity are power laws in r. For SIGMA, we do not take
! the log.
!
! We do not need to interpolate T, and ED directly, since these are part
! of POPS.
!
	LOG_R_OLD=LOG(R_OLD)
	LOG_R=LOG(R)
!
	TA(1:ND)=LOG(V(1:ND))
	CALL MON_INTERP(V,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	V(1:ND)=EXP(V(1:ND))
!
	TA(1:ND)=SIGMA(1:ND)
	CALL MON_INTERP(SIGMA,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
!
	DO I=1,NT
	  TA(1:ND)=LOG(POPS(I,1:ND))
	  CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	  POPS(I,1:ND)=EXP(TB(1:ND))
	END DO
!
! Now need to intepolate densities, and the cliumping factor, which all
! depend on the adopted R grid.
!
	TA(1:ND)=LOG(CLUMP_FAC(1:ND))
	CALL MON_INTERP(CLUMP_FAC,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	CLUMP_FAC(1:ND)=EXP(CLUMP_FAC(1:ND))
!
	TA(1:ND)=LOG(POP_ATOM(1:ND))
	CALL MON_INTERP(POP_ATOM,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	POP_ATOM(1:ND)=EXP(POP_ATOM(1:ND))
!
	TA(1:ND)=LOG(DENSITY(1:ND))
	CALL MON_INTERP(DENSITY,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	DENSITY(1:ND)=EXP(DENSITY(1:ND))
!
	DO ISPEC=1,NUM_SPECIES
	  IF(POP_SPECIES(1,ISPEC) .NE. 0.0_LDP)THEN
	    TA=LOG(POP_SPECIES(:,ISPEC))
	    CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	    POP_SPECIES(:,ISPEC)=EXP(TB)
	  END IF
	END DO
!
! Finally we compute quadrature weights etc which are grid dependent.
!
! Compute impact parameter values P
!
	CALL IMPAR(P,R,R_OLD(ND),NC,ND,NP)
!
! Compute the angular quadrature weights
!
	IF(TRAPFORJ)THEN
	  CALL NORDANGQW(JQW,R,P,W1,W2,W3,NC,ND,NP,JTRPWGT)
	  CALL NORDANGQW(HQW,R,P,W1,W2,W3,NC,ND,NP,HTRPWGT)
	  CALL NORDANGQW(KQW,R,P,W1,W2,W3,NC,ND,NP,KTRPWGT)
	  MID=.TRUE.
	  CALL GENANGQW(HMIDQW,R,P,W1,W2,W3,NC,ND,NP,HTRPWGT,MID)
	  CALL GENANGQW(NMIDQW,R,P,W1,W2,W3,NC,ND,NP,NTRPWGT,MID)
	ELSE
	  CALL NORDANGQW(JQW,R,P,W1,W2,W3,NC,ND,NP,JWEIGHT)
	  CALL NORDANGQW(HQW,R,P,W1,W2,W3,NC,ND,NP,HWEIGHT)
	  CALL NORDANGQW(KQW,R,P,W1,W2,W3,NC,ND,NP,KWEIGHT)
	  MID=.TRUE.
	  CALL GENANGQW(HMIDQW,R,P,W1,W2,W3,NC,ND,NP,HWEIGHT,MID)
	  CALL GENANGQW(NMIDQW,R,P,W1,W2,W3,NC,ND,NP,NWEIGHT,MID)
	END IF
!
	DO LS=1,NP
	  MU_AT_RMAX(LS)=SQRT( 1.0_LDP -(P(LS)/R(1))**2 )
	  HQW_AT_RMAX(LS)=HQW(1,LS)
	END DO
!
	RETURN
	END
