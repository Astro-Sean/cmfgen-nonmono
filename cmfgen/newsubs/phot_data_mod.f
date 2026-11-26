!
! Data Module to store cross-sectional date for computing the photoionization
! cross-sections.
!
	MODULE PHOT_DATA_MOD
	USE SET_KIND_MODULE
!!
	  TYPE PHOTO_DATA
!
! Altered 28-Mar-2003 : NPD_MAX variable installed.
! Altered Sep-18-1996 : Dimension of LST_U, LST_CROSS, LST_LOC changed.
!                       (LST_U renamed to LST_FREQ).
!
	   INTEGER MAX_TERMS		!# of terms (cross-sections).
	   INTEGER MAX_CROSS		!# of cross sections/term.
	   INTEGER N_FILES		!Number of files read in
	   INTEGER NUM_PHOT_ROUTES	!# of photoionization routes
!
! For clarity define NLOC as the maximum number of levels in ion.
!                    NT as the number tables read in.
!                    NS as the numer of photoionization paths (NUM_PHOT_ROUTES).
!
! Data block for photoionization cross section data. The arrays are allocated
! as the photoionization data is read in.
!
	  REAL(KIND=LDP), POINTER :: NU_NORM(:)		!MAX_CROSS
	  REAL(KIND=LDP), POINTER :: CROSS_A(:)         !MAX_CROSS
	  REAL(KIND=LDP), POINTER :: EXC_FREQ(:)        !NT
	  REAL(KIND=LDP), POINTER :: GION(:)           	!NT
!
	  INTEGER, POINTER :: A_ID(:,:)    		!NLOC,NS
          REAL(KIND=LDP), POINTER :: NEF(:,:)           !NLOC,NS
!
	  INTEGER, POINTER :: ST_LOC(:,:)    		!MAX_TERMS,NT
	  INTEGER, POINTER :: END_LOC(:,:)    		!MAX_TERMS,NT
	  INTEGER, POINTER :: CROSS_TYPE(:,:)    	!MAX_TERMS,NT
	  LOGICAL*4, POINTER :: DO_PHOT(:,:)    	!NT,NT
!
	  REAL(KIND=LDP), POINTER :: LST_CROSS(:,:)	!NLOC,NS
	  REAL(KIND=LDP), POINTER :: LST_FREQ(:,:)	!NLOC,NS
	  REAL(KIND=LDP), POINTER :: LST_LOC(:,:)	!NLOC,NS
	  REAL(KIND=LDP), POINTER :: SCALE_FAC(:)      	!NT
!
	  INTEGER, POINTER :: PHOT_GRID(:)		!NT
	  INTEGER, POINTER :: ION_LEV_ID(:)		!NT

	  REAL(KIND=LDP) AT_NO
	  REAL(KIND=LDP) ZION
	  REAL(KIND=LDP) ALPHA_BF
!
	  LOGICAL*4 DO_KSHELL_W_GS
!
! Data block for dielectronic transitions. Arrays are allocated when the
! dielectronic data is read in.
!
	  INTEGER NDIE_MAX
	  REAL(KIND=LDP), POINTER :: OSC(:)			!NDIE_MAX
	  REAL(KIND=LDP), POINTER :: GAMMA(:)			!NDIE_MAX
!
! NU_ZERO is the central frequency in 10^15 Hz.
! NU_MIN is the minimum frequency for which the line has to be considered.
! NU_MAX is the maximum frequency for which the line has to be considered.
!
	  REAL(KIND=LDP), POINTER :: NU_ZERO(:)		!NDIE_MAX
	  REAL(KIND=LDP), POINTER :: NU_MAX(:)		!NDIE_MAX
	  REAL(KIND=LDP), POINTER :: NU_MIN(:)		!NDIE_MAX
!
! The Dielectronic lines are arranged according to lower levels of the
! recombination transitions:
!
! ST_INDEX refers to the first dielectronic line associated with level I.
! END_INDX refers to the final dielectronic line associated with level I.
!
	  INTEGER, POINTER :: ST_INDEX(:)		!NLOC
	  INTEGER, POINTER :: END_INDEX(:)		!NLOC
!
	  REAL(KIND=LDP) VSM_KMS
	  INTEGER NUM_DIE
!
! Data block for free-free dielectronic transitions. Arrays are allocated when the
! free-free dielectronic data is read in.
!
	  INTEGER NUM_FF
!
! FF_GF    is the gf value for the line.
! FF_GAMMA is the Lorentz paramter, normalized by 10^{-15}
!              It thus should be number < 10, since autoionization
!              widths are typically 10^14 or less.
!
	  REAL(KIND=LDP), POINTER :: FF_GF(:)
	  REAL(KIND=LDP), POINTER :: FF_GAMMA(:)
!
! FF_NU_ZERO   is the central frequency of the line in 10^15 Hz.
! FF_NU_EXCITE is the excitation frequency of the lower level in 10^15 Hz,
!                 realtive to the ion. It is positive.
! FF_NU_MIN    is the minimum frequency for which the line has to be considered.
! F_NU_MAX     is the maximum frequency for which the line has to be considered.
!
	  REAL(KIND=LDP), POINTER :: FF_NU_ZERO(:)
	  REAL(KIND=LDP), POINTER :: FF_NU_EXCITE(:)
	  REAL(KIND=LDP), POINTER :: FF_NU_MAX(:)
	  REAL(KIND=LDP), POINTER :: FF_NU_MIN(:)
!
	  CHARACTER(LEN=20), ALLOCATABLE :: FINAL_ION_NAME(:)
!
	END TYPE PHOTO_DATA
!
	INTEGER, PARAMETER :: NPD_MAX=200
	TYPE (PHOTO_DATA) PD(NPD_MAX)
!
! CONV_FAC is the factor required to convert from Megabarns to program units
! such that R is in units of 10^10 cm, and CHI.R is dimensionless.
!
	REAL(KIND=LDP) CONV_FAC
	REAL(KIND=LDP) LG10_CONV_FAC
	PARAMETER (CONV_FAC=1.0E-08_LDP)
	PARAMETER (LG10_CONV_FAC=-8.0_LDP)
!
	END MODULE PHOT_DATA_MOD
