      SUBROUTINE OBS_V9(ETA_CMF,CHI_CMF,FREQ_CMF,
     1            R,V,T,ED,ND,NCF,
     2            P,MU_AT_RMAX,HQW_AT_RMAX,NC,NP,
     3            OBS_FREQ,OBS_FLUX,NOS,
     4            MAX_DEL_V_RES_ZONE,
     5            TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,
     6            TAU_REF,WRITE_RTAU,WRITE_IP,WRITE_dFR,
     7            DO_REL_CORRECTIONS,PLANE_PARALLEL,
     8            IB_J_CMF,IB_H_CMF,DO_IB_DIFF_FLUX)
      USE SET_KIND_MODULE
      IMPLICIT NONE
C
C Adapter that keeps one call signature in drv_obs_frame_norm so the
C driver can link against either this tree's OBS_FRAME_SUB_V9 (this
C file) or the upstream reference (obs_v9_adapter_ref.f).
C
      INTEGER ND,NCF,NC,NP,NOS,N_INS_OBS
      REAL(KIND=LDP) :: ETA_CMF(ND,NCF),CHI_CMF(ND,NCF),FREQ_CMF(NCF)
      REAL(KIND=LDP) :: R(ND),V(ND),T(ND),ED(ND)
      REAL(KIND=LDP) :: P(NP),MU_AT_RMAX(NP),HQW_AT_RMAX(NP)
      REAL(KIND=LDP) :: OBS_FREQ(NOS),OBS_FLUX(NOS)
      REAL(KIND=LDP) :: MAX_DEL_V_RES_ZONE(ND)
      REAL(KIND=LDP) :: TAU_MAX,ES_DTAU,TAU_REF
      REAL(KIND=LDP) :: IB_J_CMF(NCF),IB_H_CMF(NCF)
      LOGICAL WRITE_RTAU,WRITE_IP,WRITE_dFR
      LOGICAL DO_REL_CORRECTIONS,PLANE_PARALLEL,DO_IB_DIFF_FLUX
      CHARACTER(LEN=*) :: INT_METHOD
C
      CALL OBS_FRAME_SUB_V9(ETA_CMF,CHI_CMF,FREQ_CMF,
     1            R,V,T,ED,ND,NCF,
     2            P,MU_AT_RMAX,HQW_AT_RMAX,NC,NP,
     3            OBS_FREQ,OBS_FLUX,NOS,
     4            MAX_DEL_V_RES_ZONE,
     5            TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,
     6            TAU_REF,WRITE_RTAU,WRITE_IP,WRITE_dFR,
     7            DO_REL_CORRECTIONS,PLANE_PARALLEL,
     8            IB_J_CMF,IB_H_CMF,DO_IB_DIFF_FLUX)
      RETURN
      END
