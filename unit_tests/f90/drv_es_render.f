C
C File-input render driver for the P8.1 ES emissivity convention
C comparison.  Reads a manufactured ETA/CHI field set written by
C the test harness (or drv_es_conv qray) and runs the production
C observer-frame transfer OBS_V9 on it, so two emissivity
C conventions can be compared through the production render.
C
C Usage: drv_es_render <out_dir> <tag>
C   input : <out_dir>/esrend_in_<tag>.txt
C   output: <out_dir>/esrend_flux_<tag>.txt
C
C Input layout (whitespace free-format):
C   ND NCF NOS
C   NU(1..NCF)              decreasing, 10^15 Hz
C   R(1..ND)                decreasing, 1e10 cm
C   V(1..ND)                km/s
C   T(1..ND), ED(1..ND)
C   OBS_FREQ(1..NOS)
C   ETA(I,ML) rows          I=1..ND outer blocks of NCF values
C   CHI(I,ML) rows          same order
C
      PROGRAM DRV_ES_RENDER
      USE SET_KIND_MODULE
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: NDMAX=8
      INTEGER, PARAMETER :: NCFMX=1024
      INTEGER, PARAMETER :: NOSMX=16
      INTEGER, PARAMETER :: NP=65
      INTEGER, PARAMETER :: LU_IN=22
      INTEGER, PARAMETER :: LU_OUT=21
C
      REAL(KIND=LDP), ALLOCATABLE :: ETA(:,:),CHI(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: FREQ(:),R(:),V(:),T(:),ED(:)
      REAL(KIND=LDP), ALLOCATABLE :: OBS_FREQ(:),OBS_FLUX(:)
      REAL(KIND=LDP) :: P(NP),MU_AT_RMAX(NP),DMU_AT_RMAX(NP)
      REAL(KIND=LDP) :: HQW_AT_RMAX(NP)
      REAL(KIND=LDP) :: MAX_DEL_V_RES_ZONE(NDMAX)
      REAL(KIND=LDP) :: TAU_MAX,ES_DTAU,TAU_REF
      REAL(KIND=LDP) :: PI_L
      REAL(KIND=LDP), ALLOCATABLE :: IB_J(:),IB_H(:)
C
      CHARACTER(LEN=256) OUT_DIR
      CHARACTER(LEN=64) TAG
      INTEGER ND,NCF,NOS,NC,I,ML,N_INS_OBS
      LOGICAL DO_REL,PLANE_PARALLEL,DO_IB
      LOGICAL WRITE_RTAU,WRITE_IP,WRITE_DFR
      CHARACTER(LEN=8) INT_METHOD
C
      CALL GETARG(1,OUT_DIR)
      CALL GETARG(2,TAG)
      OPEN(LU_IN,FILE=TRIM(OUT_DIR)//'/esrend_in_'//
     1     TRIM(TAG)//'.txt',STATUS='OLD')
      READ(LU_IN,*)ND,NCF,NOS
      IF(ND .GT. NDMAX .OR. NCF .GT. NCFMX .OR. NOS .GT. NOSMX)
     1  STOP 'input dims exceed driver limits'
      ALLOCATE(ETA(ND,NCF),CHI(ND,NCF),FREQ(NCF),R(ND),V(ND),
     1  T(ND),ED(ND),OBS_FREQ(NOS),OBS_FLUX(NOS),
     2  IB_J(NCF),IB_H(NCF))
      READ(LU_IN,*)(FREQ(ML),ML=1,NCF)
      READ(LU_IN,*)(R(I),I=1,ND)
      READ(LU_IN,*)(V(I),I=1,ND)
      READ(LU_IN,*)(T(I),I=1,ND)
      READ(LU_IN,*)(ED(I),I=1,ND)
      READ(LU_IN,*)(OBS_FREQ(I),I=1,NOS)
      DO I=1,ND
        READ(LU_IN,*)(ETA(I,ML),ML=1,NCF)
      END DO
      DO I=1,ND
        READ(LU_IN,*)(CHI(I,ML),ML=1,NCF)
      END DO
      CLOSE(LU_IN)
C
      PI_L=4.0_LDP*ATAN(1.0_LDP)
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
      MAX_DEL_V_RES_ZONE=1.0_LDP
      N_INS_OBS=1
      TAU_REF=1.0_LDP
      TAU_MAX=15.0_LDP
      ES_DTAU=0.1_LDP
      WRITE_RTAU=.FALSE.
      WRITE_IP=.FALSE.
      WRITE_DFR=.FALSE.
      DO_REL=.FALSE.
      PLANE_PARALLEL=.FALSE.
      DO_IB=.FALSE.
      INT_METHOD='STAU'
      IB_J=1.0_LDP
      IB_H=0.2_LDP
      OBS_FLUX=0.0_LDP
C
      CALL OBS_V9(ETA,CHI,FREQ,R,V,T,ED,ND,NCF,
     1  P,MU_AT_RMAX,HQW_AT_RMAX,NC,NP,
     2  OBS_FREQ,OBS_FLUX,NOS,
     3  MAX_DEL_V_RES_ZONE,TAU_MAX,ES_DTAU,N_INS_OBS,
     4  INT_METHOD,TAU_REF,WRITE_RTAU,WRITE_IP,WRITE_DFR,
     5  DO_REL,PLANE_PARALLEL,IB_J,IB_H,DO_IB)
C
      OPEN(LU_OUT,FILE=TRIM(OUT_DIR)//'/esrend_flux_'//
     1     TRIM(TAG)//'.txt',STATUS='REPLACE')
      DO I=1,NOS
        WRITE(LU_OUT,'(2(ES27.17E3,1X))')OBS_FREQ(I),OBS_FLUX(I)
      END DO
      CLOSE(LU_OUT)
      STOP
      END
