	SUBROUTINE OUT_LINE_FORCE(ION_LINE_FORCE,FLUXMEAN,ROSSMEAN,RLUMST,
	1              ESEC,R,V,DENSITY,ION_ID,ND,NUM_IONS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NUM_IONS
!
	REAL(KIND=LDP) ION_LINE_FORCE(ND,NUM_IONS)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) DENSITY(ND)
	REAL(KIND=LDP) FLUXMEAN(ND)
	REAL(KIND=LDP) ROSSMEAN(ND)
	REAL(KIND=LDP) RLUMST(ND)
	REAL(KIND=LDP) ESEC(ND)
	CHARACTER(LEN=*) ION_ID(NUM_IONS)
!
	REAL(KIND=LDP) TA(ND)
	INTEGER I,L,ID,LU
!
	CALL GET_LU(LU,'OUT_LINE_FORCE')
!
	DO I=1,ND
	  ION_LINE_FORCE(I,:)=ION_LINE_FORCE(I,:)/ESEC(I)/RLUMST(I)
	  TA(I)=FLUXMEAN(I)/ESEC(I)
	END DO
!
	OPEN(UNIT=LU,FILE='ION_LINE_FORCE_TABLE',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(A)')' Summary of line force contributions by individual ions.'
	  WRITE(LU,'(A)')' Ion contributions are expressed as a % of total radiation force.'
	  WRITE(LU,'(A)')' At depth, continuum opacities will also be important.'
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(3X,A,500(A8))')'d',' V(km/s)',' M(t)',(TRIM(ION_ID(ID)),ID=1,NUM_IONS)
	  DO I=1,ND
	    WRITE(LU,'(I4,1X,F8.3,500(F8.2))')I,V(I),TA(I),
	1         (100.0_LDP*ION_LINE_FORCE(I,ID)/TA(I),ID=1,NUM_IONS)
	  END DO
	  WRITE(LU,'(3X,A,500(A8))')'d',' V(km/s)',' M(t)',(TRIM(ION_ID(ID)),ID=1,NUM_IONS)
        CLOSE(LU)
!
	OPEN(UNIT=LU,FILE='ION_FLUX_MEAN_OPAC',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(A)')' Summary of line force contributions by individual ions.'
	  WRITE(LU,'(A)')' Contributions are expressed as M(t)'
	  WRITE(LU,'(A)')' At depth, continuum opacities will also be important.'
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(I4,T30,A)')ND,'!Number of depth points'
	  WRITE(LU,'(I4,T30,A)')NUM_IONS,'!Number of ions'
	  WRITE(LU,'(A)')'R(10^10cm)'
	  WRITE(LU,'(10ES18.8)')(R(L),L=1,ND)
	  WRITE(LU,'(A)')'V(km/s)'
	  WRITE(LU,'(10ES14.4)')(V(L),L=1,ND)
	  WRITE(LU,'(A)')'Kappa electron scattering'
	  WRITE(LU,'(10ES14.4)')(1.0D-10*ESEC(L)/DENSITY(L),L=1,ND)
	  WRITE(LU,'(A)')'Normalized Rosseland mean opacity -- M(t)'
	  WRITE(LU,'(10ES14.4)')(ROSSMEAN(L)/ESEC(L),L=1,ND)
	  WRITE(LU,'(A)')'Normalized flux mean opacity -- M(t) '
	  WRITE(LU,'(10ES14.4)')(FLUXMEAN(L)/ESEC(L),L=1,ND)
	  DO ID=1,NUM_IONS
	    WRITE(LU,'(A)')TRIM(ION_ID(ID))
	    WRITE(LU,'(10ES14.4)')(ION_LINE_FORCE(L,ID),L=1,ND)
	  END DO
	CLOSE(UNIT=LU)
!
	RETURN
	END
