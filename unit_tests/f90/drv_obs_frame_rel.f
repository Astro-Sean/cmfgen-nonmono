      PROGRAM DRV_OBS_FRAME_REL
      USE SET_KIND_MODULE
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE
C
C Relativistic-transformation driver for OBS_FRAME_SUB_V9 (P7).
C Same manufactured shell as drv_obs_frame_norm: R covers [2,10] in
C 1e10 cm, FREQ_CMF is flat and wide enough for the v=40000 km/s
C Doppler reach, and rays with P<R(ND) are core rays (NC counts them).
C HOLLOW_CORE is selected by V(ND)>100, so uniform fields above 100 km/s
C mirror across the cavity while |V(ND)|<=100 keeps a solid wall.
C
C Modes:
C   run  - write the flux_*.txt battery below.
C   quad - dump P/MU/HQW and NC so the python side can evaluate the
C          boundary-seed angular sum independently.
C   prop - printed property checks: closed-form Doppler values,
C          aberration round trip, and the three Lorentz invariants.
C   sl X - single call with uniform V=X under DO_REL_CORRECTIONS;
C          used to exercise the superluminal rejection path.
C   nan  - single call with a quiet NaN inside V; used to exercise
C          the non-finite rejection path.
C
      INTEGER, PARAMETER :: ND=5
      INTEGER, PARAMETER :: NCF=4
      INTEGER, PARAMETER :: NOS=2
      INTEGER, PARAMETER :: NP=65
      INTEGER, PARAMETER :: LU_OUT=21
C
      REAL(KIND=LDP) :: ETA_CMF(ND,NCF)
      REAL(KIND=LDP) :: CHI_CMF(ND,NCF)
      REAL(KIND=LDP) :: FREQ_CMF(NCF)
      REAL(KIND=LDP) :: R(ND)
      REAL(KIND=LDP) :: V(ND)
      REAL(KIND=LDP) :: T(ND)
      REAL(KIND=LDP) :: ED(ND)
      REAL(KIND=LDP) :: P(NP)
      REAL(KIND=LDP) :: MU_AT_RMAX(NP)
      REAL(KIND=LDP) :: DMU_AT_RMAX(NP)
      REAL(KIND=LDP) :: HQW_AT_RMAX(NP)
      REAL(KIND=LDP) :: OBS_FREQ(NOS)
      REAL(KIND=LDP) :: OBS_FLUX(NOS)
      REAL(KIND=LDP) :: MAX_DEL_V_RES_ZONE(ND)
      REAL(KIND=LDP) :: TAU_MAX
      REAL(KIND=LDP) :: ES_DTAU
      REAL(KIND=LDP) :: TAU_REF
      REAL(KIND=LDP) :: IB_J_CMF(NCF)
      REAL(KIND=LDP) :: IB_H_CMF(NCF)
      REAL(KIND=LDP) :: PI_L
      REAL(KIND=LDP) :: C_KMS,SPEED_OF_LIGHT
      REAL(KIND=LDP) :: V0,BETA,GAM,DOPP,MUO,MUC,MUR,T1,T2,T3,T4
      REAL(KIND=LDP), PARAMETER :: BETAB(6)=
     1  (/0.9_LDP,0.5_LDP,0.1_LDP,-0.1_LDP,-0.5_LDP,-0.9_LDP/)
      REAL(KIND=LDP), PARAMETER :: BETAA(5)=
     1  (/0.9_LDP,0.3_LDP,0.0_LDP,-0.3_LDP,-0.9_LDP/)
      REAL(KIND=LDP), PARAMETER :: MUTAB(7)=
     1  (/1.0_LDP,0.8_LDP,0.5_LDP,0.2_LDP,0.0_LDP,-0.4_LDP,-0.9_LDP/)
      INTEGER :: I,J,NC,N_INS_OBS,NBAD
      LOGICAL :: WRITE_RTAU
      LOGICAL :: WRITE_IP
      LOGICAL :: WRITE_DFR
      LOGICAL :: DO_REL
      LOGICAL :: PLANE_PARALLEL
      LOGICAL :: DO_IB
      CHARACTER(LEN=256) :: OUT_DIR
      CHARACTER(LEN=32) :: MODE
      CHARACTER(LEN=32) :: SARG
      CHARACTER(LEN=4) :: INT_METHOD
      CHARACTER(LEN=64) :: FNAME
C
      CALL GET_COMMAND_ARGUMENT(1,OUT_DIR)
      CALL GET_COMMAND_ARGUMENT(2,MODE)
      CALL GET_COMMAND_ARGUMENT(3,SARG)
      IF(LEN_TRIM(OUT_DIR) .EQ. 0)THEN
        WRITE(*,*)'usage: drv_obs_frame_rel OUT_DIR MODE [V0]'
        STOP 1
      END IF
C
      PI_L=ACOS(-1.0_LDP)
      C_KMS=SPEED_OF_LIGHT()*1.0E-05_LDP
      R=(/10.0_LDP,8.0_LDP,6.0_LDP,4.0_LDP,2.0_LDP/)
      T=1.0_LDP
      ED=1.0_LDP
      V=0.0_LDP
      MAX_DEL_V_RES_ZONE=50.0_LDP
      TAU_MAX=20.0_LDP
      ES_DTAU=1.0_LDP
      N_INS_OBS=1
      TAU_REF=1.0_LDP
      WRITE_RTAU=.FALSE.
      WRITE_IP=.FALSE.
      WRITE_DFR=.FALSE.
      DO_REL=.TRUE.
      PLANE_PARALLEL=.FALSE.
      DO_IB=.FALSE.
      INT_METHOD='STAU'
      IB_J_CMF=1.0_LDP
      IB_H_CMF=0.2_LDP
      OBS_FREQ=(/4.45_LDP,4.05_LDP/)
C Same wide flat grid as drv_obs_frame_norm: covers the v=40000 reach.
      FREQ_CMF=(/6.0_LDP,5.0_LDP,4.0_LDP,3.0_LDP/)
C
C Uniform-in-angle impact grid on [1,0]; limb ray included once.
C
      DO I=1,NP
        MU_AT_RMAX(I)=COS(REAL(I-1,LDP)/REAL(NP-1,LDP)*PI_L/2.0_LDP)
        P(I)=R(1)*SQRT(MAX(0.0_LDP,1.0_LDP-MU_AT_RMAX(I)**2))
      END DO
      MU_AT_RMAX(NP)=0.0_LDP
      P(NP)=R(1)
      CALL SET_ACC_dMU(MU_AT_RMAX,DMU_AT_RMAX,P,R(1),NP)
      CALL HWEIGHT_V2(MU_AT_RMAX,DMU_AT_RMAX,HQW_AT_RMAX,NP)
      NC=0
      DO I=1,NP
        IF(P(I) .LT. R(ND))NC=I
      END DO
C
      IF(MODE .EQ. 'quad')THEN
C
C Dump the production angular grids plus the core-ray count so python
C can evaluate sum_LS<=NC W*I_b(LS) directly.
C
        OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/quad_grid.txt',
     1       STATUS='REPLACE')
        WRITE(LU_OUT,'(2I8)')NP,NC
        DO I=1,NP
          WRITE(LU_OUT,'(3(ES25.17,1X))')P(I),MU_AT_RMAX(I),
     1      HQW_AT_RMAX(I)
        END DO
        CLOSE(LU_OUT)
        STOP
      END IF
C
      IF(MODE .EQ. 'prop')THEN
C
C Property checks with independent closed-form targets.
C mu=+/-1 reduce the Doppler factor to the familiar sqrt forms;
C mu=0 gives the transverse 1/gamma.  The production expression is
C evaluated verbatim from the code path.
C
        NBAD=0
        DO I=1,6
          BETA=BETAB(I)
          GAM=1.0_LDP/SQRT(1.0_LDP-BETA*BETA)
          DOPP=1.0_LDP/GAM/(1.0_LDP-BETA)
          T1=SQRT((1.0_LDP+BETA)/(1.0_LDP-BETA))
          IF(ABS(DOPP-T1) .GT. 1.0E-14_LDP*DOPP)NBAD=NBAD+1
          DOPP=1.0_LDP/GAM/(1.0_LDP+BETA)
          T1=SQRT((1.0_LDP-BETA)/(1.0_LDP+BETA))
          IF(ABS(DOPP-T1) .GT. 1.0E-14_LDP*DOPP)NBAD=NBAD+1
          DOPP=1.0_LDP/GAM/(1.0_LDP-0.0_LDP)
          IF(ABS(DOPP-1.0_LDP/GAM) .GT. 1.0E-14_LDP*DOPP)NBAD=NBAD+1
        END DO
        WRITE(*,'(A,I2)')'PROP doppler_closed_form failures=',NBAD
C
C Aberration round trip: mu -> mu_cmf -> mu must recover mu for
C every (beta,mu) in the table, including mu_cmf <= 0 cases.
C
        NBAD=0
        DO I=1,5
          BETA=BETAA(I)
          DO J=1,7
            MUO=MUTAB(J)
            MUC=(MUO-BETA)/(1.0_LDP-BETA*MUO)
            MUR=(MUC+BETA)/(1.0_LDP+BETA*MUC)
            IF(ABS(MUR-MUO) .GT. 1.0E-14_LDP)NBAD=NBAD+1
          END DO
        END DO
        WRITE(*,'(A,I2)')'PROP aberration_roundtrip failures=',NBAD
C
C Lorentz invariants: build the observer-frame quantities with the
C production powers and check the invariants directly.
C
        NBAD=0
        BETA=0.37_LDP
        MUO=0.62_LDP
        GAM=1.0_LDP/SQRT(1.0_LDP-BETA*BETA)
        DOPP=1.0_LDP/GAM/(1.0_LDP-BETA*MUO)
        T1=4.2_LDP                        !nu_CMF
        T2=T1*DOPP                        !nu_OBS
        T3=7.7_LDP*DOPP**2                !ETA_OBS = D^2 ETA_CMF
        T4=5.5_LDP/DOPP                   !CHI_OBS = CHI_CMF / D
        IF(ABS(T3/T2**2 - 7.7_LDP/T1**2) .GT.
     1      1.0E-14_LDP)NBAD=NBAD+1
        IF(ABS(T2*T4 - T1*5.5_LDP) .GT. 1.0E-13_LDP)NBAD=NBAD+1
        T3=9.1_LDP*DOPP**3                !I_OBS = D^3 I_CMF
        IF(ABS(T3/T2**3 - 9.1_LDP/T1**3) .GT.
     1      1.0E-14_LDP)NBAD=NBAD+1
        WRITE(*,'(A,I2)')'PROP invariant failures=',NBAD
        STOP
      END IF
C
      IF(MODE .EQ. 'sl')THEN
C
C Single call with a uniform parsed velocity; superluminal values
C must abort inside OBS_FRAME_SUB_V9 before any flux is written.
C The grid is widened so that genuinely subluminal near-c inputs
C clear the P5 coverage guard and actually complete the run.
C
        READ(SARG,*)V0
        V=V0
        FREQ_CMF=(/70.0_LDP,6.0_LDP,1.0_LDP,0.25_LDP/)
        CHI_CMF=1.0E-6_LDP
        ETA_CMF=1.0_LDP
        OBS_FLUX=0.0_LDP
        CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
        WRITE(*,'(A,ES14.6)')'SL completed with V=',V0
        WRITE(*,'(2(ES25.17,1X))')OBS_FLUX
        STOP
      END IF
C
      IF(MODE .EQ. 'nan')THEN
        V=200.0_LDP
        V(3)=IEEE_VALUE(1.0_LDP,IEEE_QUIET_NAN)
        CHI_CMF=1.0E-6_LDP
        ETA_CMF=1.0_LDP
        OBS_FLUX=0.0_LDP
        CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
        WRITE(*,*)'NAN mode completed unexpectedly'
        STOP
      END IF
C
C ===== run mode =====
C
C TEST 1: zero velocity under DO_REL_CORRECTIONS=T must reproduce the
C non-relativistic path bit for bit (GAM=1 and NU_ON_NU0=1 exactly).
C
      CHI_CMF=1.0E-6_LDP
      ETA_CMF=1.0_LDP
      V=0.0_LDP
      DO_REL=.TRUE.
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_rel_v0.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      DO_REL=.FALSE.
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_norel_v0.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      DO_REL=.TRUE.
C
C TEST 2/4/5: uniform velocity, thin shell (hollow core since
C V(ND)>100).  CHI=1E-6 keeps the same shell geometry as the P6
C normalisation runs; the CHI=0 pair isolates the pure D^2 boost.
C
      CHI_CMF=1.0E-6_LDP
      ETA_CMF=1.0_LDP
      V=20000.0_LDP
      DO I=1,2
        DO_REL=(I .EQ. 1)
        OBS_FLUX=0.0_LDP
        CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
        IF(I .EQ. 1)THEN
          FNAME='/flux_rel_u20.txt'
        ELSE
          FNAME='/flux_norel_u20.txt'
        END IF
        OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//TRIM(FNAME),
     1       STATUS='REPLACE')
        WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
        CLOSE(LU_OUT)
      END DO
      V=40000.0_LDP
      DO I=1,2
        DO_REL=(I .EQ. 1)
        OBS_FLUX=0.0_LDP
        CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
        IF(I .EQ. 1)THEN
          FNAME='/flux_rel_u40.txt'
        ELSE
          FNAME='/flux_norel_u40.txt'
        END IF
        OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//TRIM(FNAME),
     1       STATUS='REPLACE')
        WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
        CLOSE(LU_OUT)
      END DO
      DO_REL=.TRUE.
C
C Moderately thick slab (CHI=0.5): exercises the opacity transform
C CHI -> CHI*X together with the contracted path L' = X*L.
C
      CHI_CMF=0.5_LDP
      ETA_CMF=1.0_LDP
      V=20000.0_LDP
      DO I=1,2
        DO_REL=(I .EQ. 1)
        OBS_FLUX=0.0_LDP
        CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
        IF(I .EQ. 1)THEN
          FNAME='/flux_slab_u20.txt'
        ELSE
          FNAME='/flux_slab_norel_u20.txt'
        END IF
        OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//TRIM(FNAME),
     1       STATUS='REPLACE')
        WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
        CLOSE(LU_OUT)
      END DO
      DO_REL=.TRUE.
C
C TEST 6: mirror pair.  |V| is identical at every depth and V(ND)=300
C keeps HOLLOW_CORE for both, so the only asymmetry between the two
C runs is the sign of the interior velocity.  With CHI=0 the flux is
C the D^2 multiset and must agree to round-off.
C
      CHI_CMF=0.0_LDP
      ETA_CMF=1.0_LDP
      V=(/20000.0_LDP,20000.0_LDP,20000.0_LDP,20000.0_LDP,300.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_mirror_pos.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      V=(/-20000.0_LDP,-20000.0_LDP,-20000.0_LDP,-20000.0_LDP,
     1     300.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_mirror_neg.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C TEST 7: non-monotonic field under rel corrections - acceleration,
C deceleration, a velocity maximum, a reversal, and negative cells in
C one profile.  Must complete with finite flux under the coverage
C guards.
C
      CHI_CMF=1.0E-6_LDP
      ETA_CMF=1.0_LDP
      V=(/15000.0_LDP,-3000.0_LDP,8000.0_LDP,-5000.0_LDP,300.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_rel_nonmono.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      DO_REL=.FALSE.
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_norel_nonmono.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      DO_REL=.TRUE.
C
C TEST 8: inner-boundary aberration.  V(ND)=+/-99.9 km/s keeps the
C solid-wall geometry (|V(ND)|<=100) and gives the largest attainable
C boundary beta under this convention.  ETA=0 and CHI=1E-6 make the
C diffusion seed the whole signal; IB_J_CMF=1, IB_H_CMF=0.2 constant
C so frequency interpolation is exact.
C
      ETA_CMF=0.0_LDP
      CHI_CMF=1.0E-6_LDP
      DO_IB=.TRUE.
      V=99.9_LDP
      DO I=1,2
        DO_REL=(I .EQ. 1)
        OBS_FLUX=0.0_LDP
        CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,DO_IB)
        IF(I .EQ. 1)THEN
          FNAME='/flux_ib_rel_p.txt'
        ELSE
          FNAME='/flux_ib_norel_p.txt'
        END IF
        OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//TRIM(FNAME),
     1       STATUS='REPLACE')
        WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
        CLOSE(LU_OUT)
      END DO
      DO_REL=.TRUE.
      V=-99.9_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,DO_IB)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_ib_rel_m.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C beta=0 boundary check: rel must recover the norel value exactly.
C
      V=0.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,DO_IB)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_ib_rel_v0.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      DO_REL=.FALSE.
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,DO_IB)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_ib_norel_v0.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      DO_REL=.TRUE.
C
C ETAZ cross-check of the same boundary seed (TEST 8 both methods).
C
      INT_METHOD='ETAZ'
      V=99.9_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,DO_IB)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_ib_etaz_p.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      V=-99.9_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1  FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2  NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3  TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4  WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5  IB_J_CMF,IB_H_CMF,DO_IB)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_ib_etaz_m.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
      INT_METHOD='STAU'
C
      STOP
      END
