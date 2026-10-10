!
! Simple subroutine to estimate amount of gamma-ray energy absorbed. This routine
! assumes the pure absorption approximation.
!
! It will be replaced by a more sophisticated version.
!
	SUBROUTINE DO_GAM_ABS_APPROX_V2(LOCAL_ABS_ENERGY,TOTAL_DECAY_ENERGY,KINETIC_DECAY_ENERGY,ND)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered: 26-Aug-2023 -- Updated to use new JTRPWGT_V2 and KTRPWGT_V2.
! Altered: 29-Oct-2018 -- Fixed up issues with clumping. Clumping can now vary with depth.
! Altered: 27-May-2016 -- CHI is now multiplied by CLUMP_FAC so that clumping is correctly accounted for.
! Altered: 29-Jan-2015 -- LOCAL_ABS_ENEGRY added to call. Fixed bug in luminosity calculation.
! Altered: 06-Jan-2015 -- KINETIC_DECAY_ENERGY included in the call.
!                         Changed to V2. Kinetic energy is assumed to be
!                         absorbed locally.
! Created: 29-May-2009
!
	INTEGER ND
	INTEGER ND_EXT,NC_EXT,NP_EXT
!
	REAL(KIND=LDP) LOCAL_ABS_ENERGY(ND)     	!Returned - energy absorbed locally
	REAL(KIND=LDP) TOTAL_DECAY_ENERGY(ND)           !Passed   - total energey EMITTED locally
	REAL(KIND=LDP) KINETIC_DECAY_ENERGY(ND)         !Passed   - energy local emitted as kinetic energy
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) CHI(ND)
!
	REAL(KIND=LDP) R_EXT(3*ND-6)
!
	INTEGER I,J,ISPEC
!
	ND_EXT=3*ND-6
	NC_EXT=50
	NP_EXT=ND_EXT+NC_EXT
!
! Determine number of electron per baryon.
!
	TA(1:ND)=0.5_LDP                        !Number of electrons per baryon
	DO ISPEC=1,NUM_SPECIES
	  IF('HYD' .EQ. SPECIES(ISPEC))THEN
	   TA(1:ND)=0.5_LDP*(1.0_LDP+POP_SPECIES(1:ND,ISPEC)/POP_ATOM(1:ND))
	  END IF
	END DO
!
! Compute the absorbative opacity. This will later be corrected for clumping.
!
	CHI(1:ND)=0.06_LDP*TA(1:ND)*DENSITY(1:ND)*1.0E+10_LDP
!
! Calculate a larger R grid. We attempt to keep the outer grid spacing small.
!
	R_EXT(1)=R(1)
	R_EXT(2)=R(1)+(R(2)-R(1))/3.0_LDP
	R_EXT(3)=R(1)+(R(3)-R(1))/3.0_LDP
	R_EXT(4)=R(1)+(R(3)-R(1))/1.5_LDP
	I=3
	J=5
	DO I=3,ND-2
	  R_EXT(J)=R(I)
	  R_EXT(J+1)=R(I)+(R(I+1)-R(I))/3.0_LDP
	  R_EXT(J+2)=R(I)+(R(I+1)-R(I))/1.5_LDP
	  J=J+3
	END DO
	R_EXT(ND_EXT-1)=R(ND)+(R(ND-1)-R(ND))/3.0_LDP
	R_EXT(ND_EXT)=R(ND)
!
	DO I=1,ND_EXT-1
	  IF(R_EXT(I+1) .GE. R_EXT(I))THEN
	    WRITE(6,*)'Error -- invalid extended R grid in DO_GAM_ABS_APPROX'
	    WRITE(6,*)'Depth is',I
	    WRITE(6,*)R_EXT
	    STOP
	  END IF
	END DO
!
! Now do the transfer.
!
	CALL SUB_GAM_ABS_APPROX_V2(R_EXT,ND_EXT,NC_EXT,NP_EXT,R,V,CHI,CLUMP_FAC,
	1          LOCAL_ABS_ENERGY,TOTAL_DECAY_ENERGY,KINETIC_DECAY_ENERGY,ND)
!
	RETURN
	END
!
	SUBROUTINE SUB_GAM_ABS_APPROX_V2(R,ND,NC,NP,SM_R,SM_V,SM_CHI,SM_CLUMP_FAC,
	1              LOCAL_ABS_ENERGY,TOTAL_DECAY_ENERGY,KINETIC_DECAY_ENERGY,SM_ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! ALtered 20-Nov-2022 : Fixed bug affecting the computation of the diagnostics for the energy
!                          emitted and absorbed (I was summing to ND instead of ND_SM).
!
	INTEGER SM_ND
	INTEGER ND,NC,NP
	REAL(KIND=LDP) R(ND)
!
	REAL(KIND=LDP) SM_R(SM_ND)
	REAL(KIND=LDP) SM_V(SM_ND)
	REAL(KIND=LDP) SM_CHI(SM_ND)
	REAL(KIND=LDP) SM_CLUMP_FAC(SM_ND)
	REAL(KIND=LDP) LOCAL_ABS_ENERGY(SM_ND)
	REAL(KIND=LDP) TOTAL_DECAY_ENERGY(SM_ND)
	REAL(KIND=LDP) KINETIC_DECAY_ENERGY(SM_ND)
!
	INTEGER, PARAMETER :: IONE=1
!
	REAL(KIND=LDP) P(NP)
	REAL(KIND=LDP) V(ND)
!
	REAL(KIND=LDP) Z(ND),RJ(ND),Q(ND),F(ND)
	REAL(KIND=LDP) DTAU(ND),XM(ND),ETA(ND)
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND)
	REAL(KIND=LDP) SOURCE(ND),CHI(ND),DCHIDR(ND),THETA(ND)
	REAL(KIND=LDP) JQW(ND,NP),KQW(ND,NP)
	REAL(KIND=LDP) WM(ND,ND),FB(ND,ND)
!
	REAL(KIND=LDP) ABS_OPT_DEPTH
	REAL(KIND=LDP) IC,T1,T2,T3,CONV_FAC,DBB,HBC_J,HBC_S,INBC
	LOGICAL THK_CONT
	CHARACTER(LEN=6) METHOD
        CHARACTER(LEN=9) INNER_BND_METH
!
	INTEGER I
	INTEGER IOS
	CHARACTER(LEN=132) STRING
        EXTERNAL JTRPWGT_V2,HTRPWGT_V2,KTRPWGT_V2,NTRPWGT_V2
!
! NB: When using FQCOM_INC_V2, DBB=0 and the diffusion approximation
! will give the same result as a ZERO_FLUX option, and this will be
! essentially the same as using the HOLLOW_CORE option. The later
! will potentially be out because of the small velcoity shifts induced
! by the non-zero velocity of the inner boudary.
!
	METHOD='LOGLOG'
	INNER_BND_METH='DIFFUSION'
	THK_CONT=.FALSE.
	DBB=0.0_LDP
!
! Estimate effective X-ray absorbative optical depth to the inner core.
!
	TA(1:SM_ND)=SM_CLUMP_FAC(1:SM_ND)*SM_CHI(1:SM_ND)
	T1=0.0_LDP
	DO I=1,SM_ND-1
	  T1=T1+(SM_R(I)-SM_R(I+1))*(TA(I)+TA(I+1))
	END DO
	ABS_OPT_DEPTH=0.5_LDP*T1
!
	CALL MON_INTERP(  V,ND,IONE,R,ND,SM_V,        SM_ND,SM_R,SM_ND)
!
! Correct the opacity for clumping, and interpolate it onto finer grid.
!
	TA(1:SM_ND)=SM_CLUMP_FAC(1:SM_ND)*SM_CHI(1:SM_ND)
	CALL MON_INTERP(CHI,ND,IONE,R,ND,TA,      SM_ND,SM_R,SM_ND)
!
! Correct the emissivity for clumping, and interpolate it onto finer grid.
!
	TA(1:SM_ND)=SM_CLUMP_FAC(1:SM_ND)*(TOTAL_DECAY_ENERGY(1:SM_ND)-KINETIC_DECAY_ENERGY(1:SM_ND))
	CALL MON_INTERP(ETA,ND,IONE,R,ND,TA,SM_ND,SM_R,SM_ND)
!
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
        CALL GENANGQW_V2(JQW,R,P,NC,ND,NP,JTRPWGT_V2,.FALSE.)
        CALL GENANGQW_V2(KQW,R,P,NC,ND,NP,KTRPWGT_V2,.FALSE.)
!
! The emissivity is actually ETA/4PI, but we would need to multiply again to get
! the absorbed energy.
!
	SOURCE=ETA/CHI
	CALL FQCOMP_IBC_V2(TA,TB,TC,XM,DTAU,R,Z,P,Q,F,
	1      SOURCE,CHI,DCHIDR,JQW,KQW,DBB,HBC_J,HBC_S,
	1      INBC,IC,THK_CONT,INNER_BND_METH,NC,ND,NP,METHOD)
!
! RJ is the gamma-ray mean intensity.
!
	RJ=XM
!
! Interpolate absorbed energy on to the CMFGEN grid, and add in the locally deposited
! kinetic energy.
!
! We mulitply by SM_CHI since we want the energy absorbed in the CLUMPS (not a volume
! avearged absorption).
!
	CALL MON_INTERP(TA,SM_ND,IONE,SM_R,SM_ND,RJ,ND,R,ND)
	TA(1:SM_ND)=TA(1:SM_ND)*SM_CHI(1:SM_ND)
	LOCAL_ABS_ENERGY(1:SM_ND)=TA(1:SM_ND)+KINETIC_DECAY_ENERGY(1:SM_ND)
!
! To get the total energy absorbed we need to include the clumping factor.
!
	DO I=1,SM_ND
	  TA(I)=TOTAL_DECAY_ENERGY(I)*SM_CLUMP_FAC(I)*SM_R(I)*SM_R(I)
	  TB(I)=LOCAL_ABS_ENERGY(I)*SM_CLUMP_FAC(I)*SM_R(I)*SM_R(I)
	  TC(I)=KINETIC_DECAY_ENERGY(I)*SM_CLUMP_FAC(I)*SM_R(I)*SM_R(I)
	END DO
	CALL LUM_FROM_ETA(TA,SM_R,SM_ND)
	CALL LUM_FROM_ETA(TB,SM_R,SM_ND)
	CALL LUM_FROM_ETA(TC,SM_R,SM_ND)
!
! The conversion factor:
!   (a) 4pi r^2 dr  --> 4pi .10^30 (as R is in units of 10^10 cm).
!   (b) The extra factor of 4 arises as ATAN(1.0D0) is pi/4.
!   (c) /Lsun to convert to solar luminosities
!
! NB: TA etc are decarled to be of length ND, thus I need to explictly give the
!       upper limit in the SUM.

	CONV_FAC=16.0_LDP*ATAN(1.0_LDP)*1.0E+30_LDP/3.826E+33_LDP
	T1=SUM(TA(1:SM_ND))*CONV_FAC
	T2=SUM(TB(1:SM_ND))*CONV_FAC
	T3=SUM(TC(1:SM_ND))*CONV_FAC
!
	OPEN(UNIT=10,FILE='check_edep.dat',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A,ES13.5,A)')'!                Radioactive energy emitted is :',T1,' Lsun'
	  WRITE(10,'(A,ES13.5,A)')'!                Radioactive energy absorbed is:',T2,' Lsun'
	  WRITE(10,'(A,ES13.5,A)')'!    Fraction of radioactive energy absorbed is:',T2/T1
	  WRITE(10,'(A,ES13.5,A)')'!          Fraction absorbed that is kinetic is:',T3/T2
	  WRITE(10,'(A,ES13.5,A)')'! Effec. absorbative gamma-ray optical depth is:',ABS_OPT_DEPTH
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A)')'! The energies have been averaged over volume, and thus can be'
	  WRITE(10,'(A)')'! directly compard with an identical unclumped model.'
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A)')'! V(kms)   E(non local)   E(local)   Eg(local)'
	  WRITE(10,'(A)')'!'
!
	  DO I=1,SM_ND
	    T1=SM_CLUMP_FAC(I)
	    WRITE(10,'(5ES14.6)')SM_V(I),LOCAL_ABS_ENERGY(I)*T1,TOTAL_DECAY_ENERGY(I)*T1,
	1                          KINETIC_DECAY_ENERGY(I)*T1
	  END DO
!
	 IF(MINVAL(SM_CLUMP_FAC) .LT. 0.99999_LDP)THEN
	   WRITE(10,'(A)')'!'
	   WRITE(10,'(A)')'!'
	   WRITE(10,'(A)')'! These energies have not been averaged over volume'
	   WRITE(10,'(A)')'!'
	   WRITE(10,'(A)')'! V(kms)   E(non local)   E(local)   Eg(local)'
	   WRITE(10,'(A)')'!'
	   DO I=1,SM_ND
	     WRITE(10,'(5ES14.6)')SM_V(I),LOCAL_ABS_ENERGY(I),TOTAL_DECAY_ENERGY(I),
	1                          KINETIC_DECAY_ENERGY(I)
	   END DO
	 END IF
!
! Outout the results on the fine grid -- useful for checking.
! Get total decay energy on grid -- not just that due to gamma-rays.
!
	  TA(1:SM_ND)=SM_CLUMP_FAC(1:SM_ND)*LOCAL_ABS_ENERGY(1:SM_ND)
	  CALL MON_INTERP(RJ,ND,IONE,R,ND,TA,SM_ND,SM_R,SM_ND)
!
	  TA(1:SM_ND)=SM_CLUMP_FAC(1:SM_ND)*TOTAL_DECAY_ENERGY(1:SM_ND)
	  CALL MON_INTERP(TB,ND,IONE,R,ND,TA,SM_ND,SM_R,SM_ND)
	  WRITE(10,'(A)')'!'
	  WRITE(10,'(A)')'! Large grid'
	  WRITE(10,'(A)')'! V(kms)   E(non local)   E(local)   Eg(local)'
	  WRITE(10,'(A)')'!'
	  DO I=1,ND
	    WRITE(10,'(5ES14.6)')V(I),RJ(I),TB(I),ETA(I)
	  END DO
!
	CLOSE(UNIT=10)
!
	RETURN
	END
