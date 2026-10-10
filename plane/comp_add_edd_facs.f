!
! Subroutine to compute additional Eddington factors required by the relativistic
! transfer routines. These additional moments can be obtained from the Eddington
! factors and moments already computed by subroutines such as FG_J_CMF_V10.
!
	SUBROUTINE COMPUTE_ADD_EDD_FACTS(
	1              H_ON_J,N_ON_J_NODE,KMID_ON_J,
	1              RJ,HNU,FEDD,GEDD,N_ON_J,HBC_CMF,NBC_CMF,INBC,
	1              IC,DBB,DIF,R,V,CHI,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) H_ON_J(ND)		!H/J at note
	REAL(KIND=LDP) N_ON_J_NODE(ND)		!N/J at node
	REAL(KIND=LDP) KMID_ON_J(ND)		!K/J at center of interval
	REAL(KIND=LDP) RJ(ND)			!Mean intensity (on node)
	REAL(KIND=LDP) HNU(ND)			!Eddington flux
	REAL(KIND=LDP) FEDD(ND)			!K/J at node
	REAL(KIND=LDP) GEDD(ND)			!N/H at center of interval
	REAL(KIND=LDP) N_ON_J(ND)		!Ni/(Ji+Jk)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) HBC_CMF			!H/J at outerboundar
	REAL(KIND=LDP) NBC_CMF			!N/J at outer boundary
	REAL(KIND=LDP) IC
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) INBC
	LOGICAL DIF			!Use Diffusion approximation?
!
! Local work arrays
!
	REAL(KIND=LDP) RN_MOM(ND)		!N moment at mid points
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) MIDR(ND)
	REAL(KIND=LDP) GAM_REL(ND)
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
	REAL(KIND=LDP) T1,T2
	INTEGER NDM1,NDM2
	INTEGER, PARAMETER :: IONE=1
	INTEGER I,J
!
	NDM1=ND-1
	NDM2=ND-2
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
!
	DO I=1,ND
	  GAM_REL(I)=SQRT(1.0_LDP/(1.0_LDP-(V(I)/C_KMS)**2))
	END DO
	DO I=1,NDM1
	  MIDR(I)=0.5_LDP*(R(I)+R(I+1))
	END DO
!
! Compute H/J at the nodes. We interpolate in r^2 H. Note that H_IN and
! H_OUT are multiplied by r^2 in REGRID_H.
!
	T1=HBC_CMF*RJ(1)
	TA(1:NDM1)=HNU(1:NDM1)*MIDR(1:NDM1)*MIDR(1:NDM1)
	T2=DBB/CHI(ND)/3.0_LDP
	IF(.NOT. DIF)T2=0.5_LDP*IC*(0.5_LDP+INBC)-INBC*RJ(ND)
	CALL REGRID_H(H_ON_J,R,TA,T1,T2,ND,TB)			!TB is work vector.
	H_ON_J(1:ND)=H_ON_J(1:ND)/RJ(1:ND)/R(1:ND)/R(1:ND)
!
!Compute N on J at the nodes. We first get N at the mid-points.
!NB: RN_MOM will contain r^2 N.
!
	DO I=1,NDM1
	  IF(GEDD(I) .EQ. 0.0_LDP)THEN
!	    RN_MOM(I)=N_ON_J(I)*(R(I)*R(I)*RJ(I)+R(I+1)*R(I+1)*RJ(I+1))
	    RN_MOM(I)=N_ON_J(I)
	  ELSE
!           RN_MOM(I)=HNU(I)*GEDD(I)*MIDR(I)*MIDR(I)
            RN_MOM(I)=HNU(I)*GEDD(I)*MIDR(I)*MIDR(I)/
	1                (R(I)*R(I)*RJ(I)+R(I+1)*R(I+1)*RJ(I+1))
	  END IF
	END DO
!
	CALL MON_INTERP(N_ON_J_NODE(2),NDM2,IONE,R(2),NDM2,RN_MOM,NDM1,MIDR,NDM1)
	DO I=2,NDM1
	  N_ON_J_NODE(I)=2.0_LDP*N_ON_J_NODE(I)        !/RJ(I)/R(I)/R(I)
	END DO
        N_ON_J_NODE(1)=NBC_CMF
        N_ON_J_NODE(ND)=N_ON_J_NODE(NDM1)
!
	DO I=1,ND
	  IF(ABS(N_ON_J_NODE(I)) .GE. 1.0_LDP)THEN
	    DO J=1,ND
	      WRITE(146,'(8ES14.4)')R(J),RJ(J),N_ON_J_NODE(J)
	    END DO
	    DO J=1,ND-1
	      WRITE(146,'(8ES14.4)')MIDR(J),N_ON_J(J),GEDD(J),RJ(J),RJ(J+1)
	    END DO
	    STOP
	  END IF
	END DO
!
!Recompute RSQN_ON_RSQJ (defined at the midpoints) to take into account GAMMA.
!
	DO I=1,NDM1
	  IF(N_ON_J(I) .NE. 0)THEN
	    T1=SQRT(1.0_LDP/(1.0_LDP-0.25_LDP*(V(I)/C_KMS+V(I+1)/C_KMS)**2))
	    N_ON_J(I) =T1*RN_MOM(I)*( R(I)*R(I)*RJ(I)+R(I+1)*R(I+1)*RJ(I+1) )/
	1         ( GAM_REL(I)*R(I)*R(I)*RJ(I)+
	1           GAM_REL(I+1)*R(I+1)*R(I+1)*RJ(I+1) )
	  END IF
	END DO
!
!Compute r^2 Gam K at midpoints.
!
	TA(1:ND)=RJ(1:ND)*FEDD(1:ND)*R(1:ND)*R(1:ND)*GAM_REL(1:ND)
	CALL MON_INTERP(KMID_ON_J,NDM1,IONE,MIDR,NDM1,TA,ND,R,ND)
	TA(1:ND)=RJ(1:ND)*R(1:ND)*R(1:ND)*GAM_REL(1:ND)
	KMID_ON_J(1:NDM1)=KMID_ON_J(1:NDM1)/(TA(1:NDM1)+TA(2:ND))
!
	RETURN
	END
