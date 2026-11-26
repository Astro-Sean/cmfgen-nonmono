!
! Subroutine to compute the opacity due to Rayleigh scattering at
! ND depth points. This routine has good accuracy (5%) redward of the
! Lyman lines. Within 1000 km/s of Lyman alpha, cros-section is continuous
! and constant.
!
! References:
!         Gavrila (Phys Rev., 1967, 163, p147)i
!         Hee-Won Lee and Hee Il Kim, 2004, MNRAS, 347, 802
!         Nussbaumer, H., Schmid, H. M., Vogel, M., 1989, A&A, 211, L27
!
	SUBROUTINE RAYLEIGH_SCAT(RKI,HYD,AHYD,EDGE_HYD,NHYD,FREQ,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 11-Apr-2005 : Major bug fix.
! Created 17-Mar-2004
!
	INTEGER NHYD
	INTEGER ND
!
	REAL(KIND=LDP) AHYD(NHYD,NHYD)		!No longer used
	REAL(KIND=LDP) HYD(NHYD,ND)
	REAL(KIND=LDP) EDGE_HYD(NHYD)		!Ground state value used only.
	REAL(KIND=LDP) RKI(ND)
	REAL(KIND=LDP) FREQ
!
! We include the oscillator strength in the file. That way we don't
! have to worry about whether the l states are split.
!
	REAL(KIND=LDP), SAVE ::  FOSC(20)
	DATA FOSC(2:20)/4.162E-01_LDP,7.910E-02_LDP,2.899E-02_LDP,1.394E-02_LDP,7.799E-03_LDP,
	1         4.814E-03_LDP,3.183E-03_LDP,2.216E-03_LDP,1.605E-03_LDP,1.201E-03_LDP,
	1         9.214E-04_LDP,7.227E-04_LDP,5.774E-04_LDP,4.686E-04_LDP,3.856E-04_LDP,
	1         3.211E-04_LDP,2.702E-04_LDP,2.296E-04_LDP,1.967E-04_LDP/
!
	REAL(KIND=LDP) EDGE_FREQ
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) DELF
	INTEGER I,IRES
!
! Use the results of Gavrila (Phys Rev., 1967, 163, p147) below
! the Lyman limit. Below 0.647 we use the fitting formula due to
! Ferland (Cloudy). Blueward of Lyman edge did a simple fit to
! results of Gavrila.
!
	IF(FREQ .GT. EDGE_HYD(1))THEN
	  T1=FREQ/EDGE_HYD(1)
	  RKI(1:ND)=RKI(1:ND)+6.65E-15_LDP*HYD(1,1:ND)*(1.0_LDP+1.66_LDP/SQRT(T1))
	  RETURN
	ELSE IF(FREQ .LT. 0.647_LDP*EDGE_HYD(1))THEN
	  T1=(FREQ/EDGE_HYD(1))**2
	  RKI(1:ND)=RKI(1:ND)+6.65E-15_LDP*HYD(1,1:ND)*
	1               T1*T1*(1.26_LDP+5.068_LDP*T1+708.3_LDP*(T1**5))
	  RETURN
	END IF
!
! Now sum up over all the resonances.
! See: Nussbaumer, H., Schmid, H. M., Vogel, M., 1989, A&A, 211, L27
!
	T1=0.0_LDP
	IRES=1
	DELF=100
	DO I=2,20
	  EDGE_FREQ=EDGE_HYD(1)*(1.0_LDP-1.0_LDP/I/I)
	  T2=(EDGE_FREQ/FREQ)**2-1.0_LDP
	  IF(DELF .GT. ABS(T2))THEN
	    DELF=ABS(T2)
	    IRES=I
	  END IF
	  T2=MAX(ABS(T2),0.0001_LDP)		!Prevent overflow
	  IF(EDGE_FREQ .LT. FREQ)T2=-T2
	  T1=T1+FOSC(I)/T2
	END DO
!
! The 0.35 gives a match to fitting formula above, and prevents resonances
! from going to zero.
!
	T1=T1*T1+0.35_LDP
	IF(IRES .GT. 4)IRES=4
	T1=MIN(1.0E4_LDP*FOSC(IRES)*FOSC(IRES),T1)
	RKI(1:ND)=RKI(1:ND)+6.65E-15_LDP*T1*HYD(1,1:ND)
!
	RETURN
	END
