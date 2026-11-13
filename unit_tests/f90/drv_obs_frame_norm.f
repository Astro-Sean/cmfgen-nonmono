      PROGRAM DRV_OBS_FRAME_NORM
      USE SET_KIND_MODULE
      IMPLICIT NONE
C
C Flux-normalisation driver for OBS_FRAME_SUB_V9 (P6).  R covers
C [2,10] in the usual 1e10 cm units; rays with P<R(ND) are core rays
C (NC counts them) and terminate at the solid inner boundary, while
C grazing rays traverse the full mirrored chord.  With a constant
C emissivity j the emergent intensity is analytic:
C     I(p) = j (sqrt(R1^2-p^2)-sqrt(RND^2-p^2))   for p <= RND
C     I(p) = 2 j sqrt(R1^2-p^2)                    for p >  RND
C (the solid core blocks the back-side half of core-ray chords).
C For a thick constant-source model I(p) -> S and int I mu dmu -> 0.5 S.
C HQW_AT_RMAX is built with the production SET_ACC_dMU + HWEIGHT_V2 so
C the run exercises the same quadrature the caller uses.
C
      INTEGER, PARAMETER :: ND=5
      INTEGER, PARAMETER :: NCF=4
      INTEGER, PARAMETER :: NOS=2
      INTEGER, PARAMETER :: NPN=5
      INTEGER, PARAMETER :: NPMAX=129
      INTEGER, PARAMETER :: LU_OUT=21
C
      REAL(KIND=LDP) :: ETA_CMF(ND,NCF)
      REAL(KIND=LDP) :: CHI_CMF(ND,NCF)
      REAL(KIND=LDP) :: FREQ_CMF(NCF)
      REAL(KIND=LDP) :: R(ND)
      REAL(KIND=LDP) :: V(ND)
      REAL(KIND=LDP) :: T(ND)
      REAL(KIND=LDP) :: ED(ND)
      REAL(KIND=LDP) :: P(NPMAX)
      REAL(KIND=LDP) :: MU_AT_RMAX(NPMAX)
      REAL(KIND=LDP) :: DMU_AT_RMAX(NPMAX)
      REAL(KIND=LDP) :: HQW_AT_RMAX(NPMAX)
      REAL(KIND=LDP) :: OBS_FREQ(NOS)
      REAL(KIND=LDP) :: OBS_FLUX(NOS)
      REAL(KIND=LDP) :: MAX_DEL_V_RES_ZONE(ND)
      REAL(KIND=LDP) :: TAU_MAX
      REAL(KIND=LDP) :: ES_DTAU
      REAL(KIND=LDP) :: TAU_REF
      REAL(KIND=LDP) :: IB_J_CMF(NCF)
      REAL(KIND=LDP) :: IB_H_CMF(NCF)
      REAL(KIND=LDP) :: PI_L
      INTEGER :: I,J
      INTEGER :: N_INS_OBS
      INTEGER :: NPL(NPN)
      INTEGER :: NP
      INTEGER :: NC
      LOGICAL :: WRITE_RTAU
      LOGICAL :: WRITE_IP
      LOGICAL :: WRITE_DFR
      LOGICAL :: DO_REL
      LOGICAL :: PLANE_PARALLEL
      CHARACTER(LEN=256) :: OUT_DIR
      CHARACTER(LEN=16) :: MODE
      CHARACTER(LEN=4) :: INT_METHOD
      CHARACTER(LEN=64) :: FNAME
C
      CALL GET_COMMAND_ARGUMENT(1,OUT_DIR)
      CALL GET_COMMAND_ARGUMENT(2,MODE)
      IF(LEN_TRIM(OUT_DIR) .EQ. 0)THEN
        WRITE(*,*)'usage: drv_obs_frame_norm OUT_DIR MODE'
        STOP 1
      END IF
C
      PI_L=ACOS(-1.0_LDP)
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
      DO_REL=.FALSE.
      PLANE_PARALLEL=.FALSE.
      IB_J_CMF=0.0_LDP
      IB_H_CMF=0.0_LDP
      OBS_FREQ=(/4.45_LDP,4.05_LDP/)
C FREQ_CMF must cover the full Doppler reach at v0=40000 km/s
C (max T1 = OBS_FREQ(1)*(1+0.134)*GAM ~ 5.09); the flat ETA/CHI arrays
C make the extra range free.
      FREQ_CMF=(/6.0_LDP,5.0_LDP,4.0_LDP,3.0_LDP/)
      NPL=(/9,17,33,65,129/)
C
C Uniform-in-angle impact grid: mu_i = cos((i-1)/(NP-1)*pi/2) spans
C [1,0], so P spans [0,R1] and the limb ray is included exactly once.
C
      NP=NPL(NPN)
      DO I=1,NP
        MU_AT_RMAX(I)=COS(REAL(I-1,LDP)/REAL(NP-1,LDP)*PI_L/2.0_LDP)
        P(I)=R(1)*SQRT(MAX(0.0_LDP,1.0_LDP-MU_AT_RMAX(I)**2))
      END DO
      MU_AT_RMAX(NP)=0.0_LDP
      P(NP)=R(1)
      NC=0
      DO I=1,NP
        IF(P(I) .LT. R(ND))NC=I
      END DO
      CALL SET_ACC_dMU(MU_AT_RMAX,DMU_AT_RMAX,P,R(1),NP)
      CALL HWEIGHT_V2(MU_AT_RMAX,DMU_AT_RMAX,HQW_AT_RMAX,NP)
C
      IF(MODE .EQ. 'quad')THEN
C
C Quadrature-only mode: dump the production MU and W grids for every
C NP in the family; the python side checks sum W*V(mu) against the
C analytic first moments.  This isolates SET_ACC_dMU + HWEIGHT_V2.
C
        DO J=1,NPN
          NP=NPL(J)
          DO I=1,NP
            MU_AT_RMAX(I)=COS(REAL(I-1,LDP)/REAL(NP-1,LDP)*PI_L/2.0_LDP)
            P(I)=R(1)*SQRT(MAX(0.0_LDP,1.0_LDP-MU_AT_RMAX(I)**2))
          END DO
          MU_AT_RMAX(NP)=0.0_LDP
          P(NP)=R(1)
          CALL SET_ACC_dMU(MU_AT_RMAX,DMU_AT_RMAX,P,R(1),NP)
          CALL HWEIGHT_V2(MU_AT_RMAX,DMU_AT_RMAX,HQW_AT_RMAX,NP)
          WRITE(FNAME,'(A,I3.3,A)')'/quad_w_np',NP,'.txt'
          OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//TRIM(FNAME),
     1         STATUS='REPLACE')
          WRITE(LU_OUT,'(I8)')NP
          DO I=1,NP
            WRITE(LU_OUT,'(3(ES25.17,1X))')MU_AT_RMAX(I),
     1        DMU_AT_RMAX(I),HQW_AT_RMAX(I)
          END DO
          CLOSE(LU_OUT)
        END DO
        STOP
      END IF
C
C CASE 1: optically thick shell, constant source S = ETA/CHI = 1.
C Emergent intensity saturates to I0 = 1 everywhere the chord is
C optically thick; the limb ray itself carries zero quadrature weight
C by construction (V(mu=0)=0 convention in HWEIGHT_V2).
C
      CHI_CMF=100.0_LDP
      ETA_CMF=100.0_LDP
C
      INT_METHOD='STAU'
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_thick_stau.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      INT_METHOD='ETAZ'
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_thick_etaz.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C CASE 2: optically thin shell, constant emissivity j = 1.  V(ND)>100
C forces HOLLOW_CORE so core rays mirror across the cavity: the
C emergent intensity is then the analytic shell chord
C     I(p) = 2j(sqrt(R1^2-p^2)-sqrt(RND^2-p^2))   p <= RND
C     I(p) = 2j sqrt(R1^2-p^2)                     p >  RND
C Run the NP series for the angular-quadrature convergence study,
C both formal methods at the finest grid.
C
      CHI_CMF=1.0E-6_LDP
      ETA_CMF=1.0_LDP
      V=(/400.0_LDP,350.0_LDP,300.0_LDP,250.0_LDP,200.0_LDP/)
      DO J=1,NPN
        NP=NPL(J)
        DO I=1,NP
          MU_AT_RMAX(I)=COS(REAL(I-1,LDP)/REAL(NP-1,LDP)*PI_L/2.0_LDP)
          P(I)=R(1)*SQRT(MAX(0.0_LDP,1.0_LDP-MU_AT_RMAX(I)**2))
        END DO
        MU_AT_RMAX(NP)=0.0_LDP
        P(NP)=R(1)
        NC=0
        DO I=1,NP
          IF(P(I) .LT. R(ND))NC=I
        END DO
        CALL SET_ACC_dMU(MU_AT_RMAX,DMU_AT_RMAX,P,R(1),NP)
        CALL HWEIGHT_V2(MU_AT_RMAX,DMU_AT_RMAX,HQW_AT_RMAX,NP)
        INT_METHOD='STAU'
        OBS_FLUX=0.0_LDP
        CALL OBS_V9(ETA_CMF,CHI_CMF,
     1      FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2      NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3      TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4      WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5      IB_J_CMF,IB_H_CMF,.FALSE.)
        WRITE(FNAME,'(A,I3.3,A)')'/flux_thin_np',NP,'.txt'
        OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//TRIM(FNAME),
     1       STATUS='REPLACE')
        WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
        CLOSE(LU_OUT)
      END DO
C
      NP=NPL(NPN)
      DO I=1,NP
        MU_AT_RMAX(I)=COS(REAL(I-1,LDP)/REAL(NP-1,LDP)*PI_L/2.0_LDP)
        P(I)=R(1)*SQRT(MAX(0.0_LDP,1.0_LDP-MU_AT_RMAX(I)**2))
      END DO
      MU_AT_RMAX(NP)=0.0_LDP
      P(NP)=R(1)
      NC=0
      DO I=1,NP
        IF(P(I) .LT. R(ND))NC=I
      END DO
      CALL SET_ACC_dMU(MU_AT_RMAX,DMU_AT_RMAX,P,R(1),NP)
      CALL HWEIGHT_V2(MU_AT_RMAX,DMU_AT_RMAX,HQW_AT_RMAX,NP)
      INT_METHOD='ETAZ'
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_thin_etaz.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C CASE 3: flat-spectrum thin shell under velocity laws that keep the
C same core geometry.  With ETA/CHI constant in frequency the Doppler
C mapping is invisible, so any flux difference between fields of the
C same geometry class is ray miscounting.  V(ND)>100 keeps
C HOLLOW_CORE=T throughout this block; the solid-boundary pair at the
C end must match for the same reason.  NP=65 keeps the runtime down;
C the quadrature itself is verified by the quad-mode weights and the
C thin-shell NP series.
C
      NP=65
      DO I=1,NP
        MU_AT_RMAX(I)=COS(REAL(I-1,LDP)/REAL(NP-1,LDP)*PI_L/2.0_LDP)
        P(I)=R(1)*SQRT(MAX(0.0_LDP,1.0_LDP-MU_AT_RMAX(I)**2))
      END DO
      MU_AT_RMAX(NP)=0.0_LDP
      P(NP)=R(1)
      NC=0
      DO I=1,NP
        IF(P(I) .LT. R(ND))NC=I
      END DO
      CALL SET_ACC_dMU(MU_AT_RMAX,DMU_AT_RMAX,P,R(1),NP)
      CALL HWEIGHT_V2(MU_AT_RMAX,DMU_AT_RMAX,HQW_AT_RMAX,NP)
      INT_METHOD='STAU'
      V=200.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_flat_v200.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=(/1500.0_LDP,1200.0_LDP,900.0_LDP,600.0_LDP,300.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_flat_mono.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=(/5000.0_LDP,1000.0_LDP,4500.0_LDP,500.0_LDP,3000.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_flat_rev.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=(/200.0_LDP,-300.0_LDP,400.0_LDP,-100.0_LDP,600.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_flat_neg.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C Solid-boundary pair: V(ND)<=100 keeps HOLLOW_CORE=F, so the two
C fields below must still agree with each other.
C
      V=50.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_flat_s50.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=(/60.0_LDP,40.0_LDP,55.0_LDP,30.0_LDP,50.0_LDP/)
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_flat_s50b.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C CASE 4: relativistic boost on the thin shell.  Uniform V = v0; the
C code applies j -> j D^2 with D = 1/GAM/(1 - v z/(r c)) pointwise
C along each mirrored chord, and the python side predicts the same
C disk integral independently.  NP=33 is enough for the comparison
C tolerance.
C
      NP=33
      DO I=1,NP
        MU_AT_RMAX(I)=COS(REAL(I-1,LDP)/REAL(NP-1,LDP)*PI_L/2.0_LDP)
        P(I)=R(1)*SQRT(MAX(0.0_LDP,1.0_LDP-MU_AT_RMAX(I)**2))
      END DO
      MU_AT_RMAX(NP)=0.0_LDP
      P(NP)=R(1)
      NC=0
      DO I=1,NP
        IF(P(I) .LT. R(ND))NC=I
      END DO
      CALL SET_ACC_dMU(MU_AT_RMAX,DMU_AT_RMAX,P,R(1),NP)
      CALL HWEIGHT_V2(MU_AT_RMAX,DMU_AT_RMAX,HQW_AT_RMAX,NP)
      DO_REL=.TRUE.
      V=20000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_rel_v20k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=40000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_rel_v40k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      DO_REL=.FALSE.
      V=20000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_norel_v20k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=40000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_norel_v40k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C ETAZ rel pair at the same velocities: the direct emissivity
C integral has a different discretisation error than the STAU
C source-function cubic, so agreement between the two rel diffs
C isolates the formal-method contribution to the boost residual.
C
      DO_REL=.TRUE.
      INT_METHOD='ETAZ'
      V=20000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_etaz_v20k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=40000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_etaz_v40k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
C Pure-boost probe: CHI=0 makes the ETAZ emissivity integral exactly
C int j*D^2 dz along each chord, isolating the applied emissivity
C boost from every attenuation effect.
C
      CHI_CMF=0.0_LDP
      INT_METHOD='ETAZ'
      DO_REL=.TRUE.
      V=20000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_boost_v20k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=40000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_boost_v40k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      DO_REL=.FALSE.
      V=20000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_boost_norel_v20k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      V=40000.0_LDP
      OBS_FLUX=0.0_LDP
      CALL OBS_V9(ETA_CMF,CHI_CMF,
     1    FREQ_CMF,R,V,T,ED,ND,NCF,P,MU_AT_RMAX,HQW_AT_RMAX,
     2    NC,NP,OBS_FREQ,OBS_FLUX,NOS,MAX_DEL_V_RES_ZONE,
     3    TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,TAU_REF,
     4    WRITE_RTAU,WRITE_IP,WRITE_DFR,DO_REL,PLANE_PARALLEL,
     5    IB_J_CMF,IB_H_CMF,.FALSE.)
      OPEN(UNIT=LU_OUT,FILE=TRIM(OUT_DIR)//'/flux_boost_norel_v40k.txt',
     1     STATUS='REPLACE')
      WRITE(LU_OUT,'(2(ES25.17,1X))')OBS_FLUX
      CLOSE(LU_OUT)
C
      STOP
      END
