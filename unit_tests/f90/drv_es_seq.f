C
C P8.1 T5 sequencing harness. Mirrors the production ES-iteration
C ordering of CMF_FLUX_SUB_V5 on the OBS path, with a scripted
C solver standing in for the COMP_JCONT + formal-solve block and
C the real COMP_J_CONV_V2 for every convolution.
C
C Production operations represented (same order, same array
C roles):
C   sweep k: consume JBAR_{k-1} (ES_J_CONV role)
C            -> scripted solve produces J_k
C            -> store J_k into ST (RJ_CMF_ST role)
C   inter-sweep: refresh slot (counter only, mirrors the
C            FG_KNOP_REFRESH_V1 / MOM_KNOP_REFRESH_V1 call site)
C            -> convolve stored ST -> JBAR_k (ES_J_CONV write)
C   render: convolve stored ST -> JBAR_N (in-place RJ_CMF_ST
C            convolve role), TA = ETA + ESEC*JBAR_N.
C
C Modes:
C   correct  - production ordering.
C   stale    - sweep 2 is fed JBAR_0 (simulates a stale read).
C   preconv  - convolution runs before J_1 is stored, so JBAR_1
C              is built from the pre-sweep ST content.
C   oldrend  - the render folds the previous file JBAR_1 instead
C              of the in-place convolve of the stored J_2.
C
C The scripted solver is J_k = LAM*JBAR_consumed + NEW_k, so the
C genealogy of every field is measurable: perturbing NEW_1 or
C NEW_2 (args) propagates through exactly the stages that consume
C it in the correct mode, and is blocked in the broken modes.
C This harness does not test CMF_FLUX internals; it exercises the
C sequencing contract on the real convolver.
C
C Usage: drv_es_seq <out_dir> <mode> [pert1] [pert2]
C Writes esseq_fields.txt: one header line then per frequency
C NU JBAR0 J1 JBAR1 J2 JBAR2 TA_RENDER  plus REFRESH_COUNT.
C
      PROGRAM DRV_ES_SEQ
      USE SET_KIND_MODULE
      IMPLICIT NONE
C
      INTEGER ND,NCF,ML,K,LU_OUT
      PARAMETER (ND=2,NCF=64)
      CHARACTER(LEN=200) OUT_DIR,MODE,SARG
      REAL(KIND=LDP) NU(NCF),T(ND)
      REAL(KIND=LDP) JBAR0(ND,NCF),JBAR1(ND,NCF),JBAR2(ND,NCF)
      REAL(KIND=LDP) J1(ND,NCF),J2(ND,NCF),ST(ND,NCF)
      REAL(KIND=LDP) ETA(ND,NCF),ESEC(ND),TA(ND,NCF)
      REAL(KIND=LDP) NEWK(ND,NCF),PERT1,PERT2,LAM,NU_HI,NU_LO
      REAL(KIND=LDP) BB,T1,T2,TWOHCSQ,HDKT
      INTEGER LUDUM,CONT_REC,REFRESH_COUNT
      LOGICAL L_FALSE
      EXTERNAL COMP_J_CONV_V2
C
      COMMON/CONSTANTS/ DUM1,DUM2,HDKT,TWOHCSQ
      REAL(KIND=LDP) DUM1,DUM2
      COMMON/LINE/ OPLIN,EMLIN
      REAL(KIND=LDP) OPLIN,EMLIN
C
      BB(T1,T2)=TWOHCSQ*(T2**3)/( EXP(HDKT*T2/T1)-1.0_LDP )
C
      HDKT=4.7994145_LDP
      TWOHCSQ=0.0147452575_LDP
      OPLIN=0.0_LDP
      EMLIN=0.0_LDP
      L_FALSE=.FALSE.
      LUDUM=7
      CONT_REC=0
      LU_OUT=41
      LAM=0.85_LDP
      REFRESH_COUNT=0
C
      CALL GETARG(1,OUT_DIR)
      CALL GETARG(2,MODE)
      PERT1=0.0_LDP
      PERT2=0.0_LDP
      CALL GETARG(3,SARG)
      IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)PERT1
      CALL GETARG(4,SARG)
      IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)PERT2
C
      T(1)=2.0_LDP
      T(2)=1.0_LDP
      NU_HI=120.0_LDP
      NU_LO=0.05_LDP
      DO ML=1,NCF
        NU(ML)=EXP( LOG(NU_HI)+(LOG(NU_LO)-LOG(NU_HI))*
     1            (ML-1.0_LDP)/(NCF-1.0_LDP) )
      END DO
C
C Manufactured fields. JBAR_0, the initial convolved field, is a
C diluted Planck; NEW_K gives each scripted sweep output its own
C spectral signature so the generations are distinguishable.
C
      DO K=1,ND
        ESEC(K)=0.05_LDP*K
        DO ML=1,NCF
          JBAR0(K,ML)=(0.3_LDP+0.4_LDP*(K-1.0_LDP))*BB(T(K),NU(ML))
          NEWK(K,ML)=0.5_LDP*JBAR0(K,ML)*
     1      EXP(-0.5_LDP*((NU(ML)-1.8_LDP)/0.6_LDP)**2)
          ETA(K,ML)=1.0D-2*ESEC(K)*JBAR0(K,ML)
        END DO
      END DO
C
C ---------------- sweep 1 ----------------
C Consume JBAR_0, produce J_1, store ST=J_1.
      DO K=1,ND
        DO ML=1,NCF
          J1(K,ML)=LAM*JBAR0(K,ML) +
     1      NEWK(K,ML)*(1.0_LDP+PERT1)
        END DO
      END DO
      ST(:,:)=J1(:,:)
C
C ---------------- inter-sweep block ----------------
C Refresh slot (FG_KNOP_REFRESH_V1 / MOM_KNOP_REFRESH_V1 site).
      REFRESH_COUNT=REFRESH_COUNT+1
C
C Convolve stored J_1 -> JBAR_1 (ES_J_CONV write). In mode
C 'preconv' the convolution runs before the store, i.e. on the
C pre-sweep content of ST, which here is JBAR_0.
      IF(MODE .EQ. 'preconv')ST(:,:)=JBAR0(:,:)
      CALL COMP_J_CONV_V2(ST,NCF*ND,NU,T,ND,NCF,
     1     LUDUM,'J PASSED VIA CALL',CONT_REC,L_FALSE,L_FALSE,
     2     LUDUM,'RETURN J VIA CALL')
      JBAR1(:,:)=ST(:,:)
C
C ---------------- sweep 2 ----------------
C Consume JBAR_1 (or JBAR_0 in 'stale'), produce J_2, store.
      IF(MODE .EQ. 'stale')JBAR1(:,:)=JBAR0(:,:)
      DO K=1,ND
        DO ML=1,NCF
          J2(K,ML)=LAM*JBAR1(K,ML) +
     1      NEWK(K,ML)*(1.0_LDP+PERT2)
        END DO
      END DO
      ST(:,:)=J2(:,:)
C
C ---------------- render ----------------
C In-place convolve of the stored J -> JBAR_2, then fold once:
C TA = ETA + ESEC*JBAR_2 (RJ_CMF_ST convolve + ETA_RND_ST role).
      CALL COMP_J_CONV_V2(ST,NCF*ND,NU,T,ND,NCF,
     1     LUDUM,'J PASSED VIA CALL',CONT_REC,L_FALSE,L_FALSE,
     2     LUDUM,'RETURN J VIA CALL')
      JBAR2(:,:)=ST(:,:)
      IF(MODE .EQ. 'oldrend')JBAR2(:,:)=JBAR1(:,:)
      DO K=1,ND
        DO ML=1,NCF
          TA(K,ML)=ETA(K,ML)+ESEC(K)*JBAR2(K,ML)
        END DO
      END DO
C
      OPEN(LU_OUT,FILE=TRIM(OUT_DIR)//'/esseq_fields.txt',
     1     STATUS='UNKNOWN')
      WRITE(LU_OUT,'(A,I3)')'# REFRESH_COUNT=',REFRESH_COUNT
      WRITE(LU_OUT,'(A)')'# ML NU JBAR0 J1 JBAR1 J2 JBAR2 TAM'
      DO ML=1,NCF
        WRITE(LU_OUT,'(I5,7(ES27.17E3))')
     1    ML,NU(ML),JBAR0(ND,ML),J1(ND,ML),JBAR1(ND,ML),
     2    J2(ND,ML),JBAR2(ND,ML),TA(ND,ML)
      END DO
      CLOSE(LU_OUT)
      STOP
      END
