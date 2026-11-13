      PROGRAM DRV_OBS_FRAME_COV
      USE SET_KIND_MODULE
      IMPLICIT NONE
C
C Frequency-coverage driver for OBS_FRAME_SUB_V9.  Same compact model
C skeleton as drv_obs_frame_ib: a small spherical grid whose ray
C directions exercise both the front and mirrored far side.  The cases
C probe the P5 coverage work: monotonic outflow, several reversals,
C negative velocities, relativistic mapping and its non-relativistic
C limit, NUM_BANDS partitioning, grid-extension stability, a line
C feature inside the Doppler reach, and (in UNDER modes) a global-grid
C undercoverage that must abort with a named diagnostic.
C
      INTEGER, PARAMETER :: ND=5
      INTEGER, PARAMETER :: NC=3
      INTEGER, PARAMETER :: NP=7
      INTEGER, PARAMETER :: NCF=4
      INTEGER, PARAMETER :: NCFX=6
      INTEGER, PARAMETER :: NOS=2
      INTEGER, PARAMETER :: NOS_B=2000
      INTEGER, PARAMETER :: LU_OUT=21
C
      REAL(KIND=LDP) :: ETA_CMF(ND,NCFX)
      REAL(KIND=LDP) :: CHI_CMF(ND,NCFX)
      REAL(KIND=LDP) :: ETA4(ND,NCF)
      REAL(KIND=LDP) :: CHI4(ND,NCF)
      REAL(KIND=LDP) :: FREQ_CMF(NCFX)
      REAL(KIND=LDP) :: FREQ_X(NCFX)
      REAL(KIND=LDP) :: R(ND)
      REAL(KIND=LDP) :: V(ND)
      REAL(KIND=LDP) :: T(ND)
      REAL(KIND=LDP) :: ED(ND)
      REAL(KIND=LDP) :: P(NP)
      REAL(KIND=LDP) :: MU_AT_RMAX(NP)
      REAL(KIND=LDP) :: HQW_AT_RMAX(NP)
      REAL(KIND=LDP) :: OBS_FREQ(NOS)
      REAL(KIND=LDP) :: OBS_FLUX(NOS)
      REAL(KIND=LDP), ALLOCATABLE :: OBS_FREQ_B(:)
      REAL(KIND=LDP), ALLOCATABLE :: OBS_FLUX_B(:)
      REAL(KIND=LDP) :: MAX_DEL_V_RES_ZONE(ND)
      REAL(KIND=LDP) :: TAU_MAX
      REAL(KIND=LDP) :: ES_DTAU
      REAL(KIND=LDP) :: TAU_REF
      REAL(KIND=LDP) :: IB_J_CMF(NCFX)
      REAL(KIND=LDP) :: IB_H_CMF(NCFX)
      INTEGER :: I,J
      INTEGER :: N_INS_OBS
      LOGICAL :: WRITE_RTAU
      LOGICAL :: WRITE_IP
      LOGICAL :: WRITE_DFR
      LOGICAL :: DO_REL
      LOGICAL :: PLANE_PARALLEL
      CHARACTER(LEN=256) :: OUT_DIR
      CHARACTER(LEN=16) :: MODE
C
      CALL GET_COMMAND_ARGUMENT(1,OUT_DIR)
      CALL GET_COMMAND_ARGUMENT(2,MODE)
      IF(LEN_TRIM(OUT_DIR) .EQ. 0)THEN
        WRITE(*,*)'usage: drv_obs_frame_cov OUT_DIR MODE'
        STOP 1
      END IF
C
      R=(/10.0_LDP,8.0_LDP,6.0_LDP,4.0_LDP,2.0_LDP/)
      T=1.0_LDP
      ED=1.0_LDP
      P=(/0.0_LDP,1.0_LDP,1.9_LDP,3.0_LDP,5.0_LDP,7.0_LDP,9.0_LDP/)
      MU_AT_RMAX=SQRT(1.0_LDP-(P/R(1))**2)
      HQW_AT_RMAX=1.0_LDP
      MAX_DEL_V_RES_ZONE=50.0_LDP
      TAU_MAX=20.0_LDP
      ES_DTAU=1.0_LDP
      N_INS_OBS=1
      TAU_REF=1.0_LDP
      WRITE_RTAU=.FALSE.
      WRITE_IP=.FALSE.
      WRITE_DFR=.FALSE.
      DO_REL=.FALSE.
      PLANE_PARALLEL=.FALSE.
      IB_J_CMF=0.0_LDP
      IB_H_CMF=0.0_LDP
      OBS_FREQ=(/4.45_LDP,4.05_LDP/)
C
C Frequency-dependent opacity so the T1 bracket actually selects
C different opacities; the line near the grid top sits just outside
C the Doppler reach of the base model.
C
      FREQ_CMF=0.0_LDP
      FREQ_CMF(1:NCF)=(/5.0_LDP,4.5_LDP,4.0_LDP,3.5_LDP/)
      FREQ_X=(/5.5_LDP,5.0_LDP,4.5_LDP,4.0_LDP,3.5_LDP,3.0_LDP/)
      DO J=1,NCFX
        DO I=1,ND
          CHI_CMF(I,J)=1.0E-8_LDP*(FREQ_X(J)/4.0_LDP)
          ETA_CMF(I,J)=1.0E-10_LDP*(FREQ_X(J)/4.0_LDP)
        END DO
      END DO
      DO J=1,NCF
        DO I=1,ND
          CHI4(I,J)=1.0E-8_LDP*(FREQ_CMF(J)/4.0_LDP)
          ETA4(I,J)=1.0E-10_LDP*(FREQ_CMF(J)/4.0_LDP)
        END DO
      END DO
C
      IF(MODE .EQ. 'under_blue')THEN
C
C Global-grid undercoverage, blue edge.  FREQ_CMF tops out below the
C required nu_obs*(1+MAX_VMU); the serial assertion must abort with
C GLOBAL_CMF_GRID_UNDERCOVERAGE, not emit a truncated spectrum.
C
        FREQ_CMF(1:NCF)=(/4.46_LDP,4.44_LDP,4.42_LDP,4.40_LDP/)
        V=(/1500.0_LDP,1200.0_LDP,900.0_LDP,600.0_LDP,300.0_LDP/)
        OBS_FLUX=0.0_LDP
        CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1      FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2      NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3      TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4      WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5      IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
        WRITE(*,*)'UNEXPECTED: undercovered run completed'
        STOP 2
C
      ELSE IF(MODE .EQ. 'under_red')THEN
C
C Same failure mode at the red edge: the grid bottom sits above the
C required red Doppler reach.
C
        FREQ_CMF(1:NCF)=(/4.60_LDP,4.30_LDP,4.20_LDP,4.10_LDP/)
        V=(/1500.0_LDP,1200.0_LDP,900.0_LDP,600.0_LDP,300.0_LDP/)
        OBS_FLUX=0.0_LDP
        CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1      FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2      NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3      TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4      WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5      IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
        WRITE(*,*)'UNEXPECTED: undercovered run completed'
        STOP 2
      END IF
C
C TEST A: monotonic outflow, V increasing outward.
C
      V=(/1500.0_LDP,1200.0_LDP,900.0_LDP,600.0_LDP,300.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_mono.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C TEST B: several velocity reversals along the profile.
C
      V=(/5000.0_LDP,1000.0_LDP,4500.0_LDP,500.0_LDP,3000.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_revers.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C TEST C: negative velocities (infall cells) in the profile.
C
      V=(/200.0_LDP,-300.0_LDP,400.0_LDP,-100.0_LDP,600.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_negv.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C TEST H: relativistic mapping.  V ~ 0.03c makes (1-VMU)*GAM differ
C from the non-rel 1+-MAX_VMU by ~1%; the exact extrema keep the
C window covered.  A second call at V=1 km/s confirms the
C non-relativistic limit.
C
      V=(/10000.0_LDP,9000.0_LDP,9500.0_LDP,7000.0_LDP,8000.0_LDP/)
      DO_REL=.TRUE.
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA_CMF,CHI_CMF,FREQ_X,
     1    R,V,T,ED,ND,NCFX,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_rel.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=1.0E-3_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_rel_lo.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      DO_REL=.FALSE.
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_norel_lo.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C TEST F: NUM_BANDS partitioning.  NOS_B=2000 gives NUM_BANDS=2; the
C 1999-point grid gives NUM_BANDS=1.  Identical observer frequencies
C in the overlap must give identical flux.
C
      ALLOCATE(OBS_FREQ_B(NOS_B),OBS_FLUX_B(NOS_B))
      DO I=1,NOS_B
        OBS_FREQ_B(I)=4.45_LDP-0.40_LDP*REAL(I-1,LDP)/REAL(NOS_B-1,LDP)
      END DO
      V=(/1500.0_LDP,1200.0_LDP,900.0_LDP,600.0_LDP,300.0_LDP/)
      OBS_FLUX_B=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ_B,OBS_FLUX_B,NOS_B,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_b2.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(ES25.17)')OBS_FLUX_B
      CLOSE(LU_OUT)
C
      OBS_FLUX_B=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ_B(1:NOS_B-1),OBS_FLUX_B(1:NOS_B-1),NOS_B-1,
     3    MAX_DEL_V_RES_ZONE,
     4    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     5    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     6    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_b1.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(ES25.17)')OBS_FLUX_B(1:NOS_B-1)
      CLOSE(LU_OUT)
      DEALLOCATE(OBS_FREQ_B,OBS_FLUX_B)
C
C TEST P5.7: convergence under grid extension.  The extended grid
C adds points only outside the required Doppler reach, so every T1
C bracket is unchanged and the flux must be bitwise identical.
C
      V=(/1500.0_LDP,1200.0_LDP,900.0_LDP,600.0_LDP,300.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_base.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA_CMF,CHI_CMF,FREQ_X,
     1    R,V,T,ED,ND,NCFX,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_ext.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C Line feature inside the Doppler reach: a CHI spike at FREQ_CMF(2)
C anchors the bracket every required T1 interpolates through, so the
C emergent flux must change - the proof the mapping samples the
C feature rather than a substituted edge value.  Only CHI is
C perturbed; ETA is unchanged from the base run.
C
      DO I=1,ND
        CHI4(I,2)=1.0E-2_LDP
      END DO
      OBS_FLUX=0.0_LDP
      CALL OBS_FRAME_SUB_V9(ETA4,CHI4,
     1    FREQ_CMF(1:NCF),R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,'STAU',TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF(1:NCF),IB_H_CMF(1:NCF),.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_feat.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      STOP
      END
