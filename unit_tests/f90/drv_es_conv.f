C
C Manufactured driver for the production electron-scattering
C convolution COMP_J_CONV_V2 (P8.1 acceptance).
C
C The routine is exercised in its 'J PASSED VIA CALL' /
C 'RETURN J VIA CALL' mode so no EDDFACTOR/ES_J_CONV files are
C needed: the manufactured J goes in, the convolved J comes back.
C
C Usage: drv_es_conv <out_dir> <mode> [arg3..arg6]
C   modes:
C     planck   J = B_nu(T(K))            - Planck preservation
C     const    J = const                 - constant-field response
C     line     J = B + narrow Gaussian   - spectral feature smearing
C     edge     narrow features near grid ends
C     tempu    J = B(TTST) uniform T     - pure T dependence of BETA
C     pos      J = B with zeroed blocks  - positivity
C     grid     J = B(T) at arg3 = NCF    - grid-refinement series
C     linepair cont and cont+feature convolved separately on the
C              production grid (arg3=NCF, arg4=NU0, arg5=DNU,
C              arg6=AMP) - baseline-subtracted smearing
C     bndp     live features at both production bounds, production
C              grid (arg3=NCF)
C     bnd      same field on the appended-node extended grid
C              (arg3=NCF_P; NCF=NCF_P+NCF_P/2 internally)
C     linext   linepair on the same extended grid - boundary
C              sensitivity of the smeared feature
C     wcp      smooth non-Planck continuum alive through the
C              Wien-dead upper end, production grid (arg3=NCF)
C     wce      same continuum on the extended grid (arg3=NCF_P)
C     swp      absolute-amplitude Gaussian at arg3=NU0 - upper
C              boundary feature-position sweep (arg4=DNU,arg5=AMP)
C     lte      J = B(T(K)) LTE trade-off diagnostic: emits the
C              replica convolved Planck field, the replicated
C              bounded PLANCK_NU map, a per-node zone flag, and
C              JBAR (arg3=NCF; arg4/arg5 = NU_HI/NU_LO grid bounds
C              for pathological-grid match-failure probes)
C     nanj     J = B with a NaN injected at depth 2 - non-finite
C              input handling
C     neg      J = B with a negative mid-grid block - exercises
C              the MINVAL(J_ES) gate ahead of LOG(J_ES)
C     brk      synthetic ES_PLANCK_NU_MATCH battery: valid, above-
C              range, below-range, floor-walk, degenerate, NaN, Inf,
C              huge-negative, boundary-exact, and persisted-ML
C              targets through the production object code
C     qray     INC_RAY manufactured atmosphere; emits both
C              emissivity conventions + all fields for the
C              CHI_RAY*(RJ_ES-RJ) measurement
C
C   arg3 = NCF for grid/linepair/linext/bnd/bndp, NU0 for
C          line/edge/tempu, AMP for qray
C   arg4 = Gaussian width (line) / NU0 (linepair,linext)
C   arg5 = peak amplitude (line) / DNU (linepair,linext)
C   arg6 = AMP (linepair,linext)
C
      PROGRAM DRV_ES_CONV
      USE SET_KIND_MODULE
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: ND=4
      INTEGER, PARAMETER :: NCFMX=1024
      INTEGER, PARAMETER :: LU_OUT=31
C
      REAL(KIND=LDP), ALLOCATABLE :: J_STORE(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: J_IN(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: RJ_CONT(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: J_LINE(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: NU(:)
      REAL(KIND=LDP) T(ND)
      REAL(KIND=LDP) NU0,DNU,AMP,NU_HI,NU_LO,TTST,DLN
      CHARACTER(LEN=20) FNAME
C
C COMMON blocks required by COMP_J_CONV_V2.
C
      COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
      COMMON/LINE/ OPLIN,EMLIN
      REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
      REAL(KIND=LDP) OPLIN,EMLIN
C
      CHARACTER(LEN=256) OUT_DIR,SARG
      CHARACTER(LEN=32) MODE
      INTEGER I,K,ML,NCF,NCF_P,NEXT,CONT_REC,LUDUM
      LOGICAL L_FALSE
      REAL(KIND=LDP) BB,T1
C
      BB(T1,NU0)=TWOHCSQ*(NU0**3)/
     1   (EXP(MIN(700.0_LDP,HDKT*NU0/T1))-1.0_LDP)
C
      CALL GETARG(1,OUT_DIR)
      CALL GETARG(2,MODE)
      NCF=256
      NU0=20.0_LDP
      DNU=0.5_LDP
      AMP=100.0_LDP
      CALL GETARG(3,SARG)
      IF(LEN_TRIM(SARG) .GT. 0)THEN
        IF(MODE .EQ. 'grid' .OR. MODE .EQ. 'linepair' .OR.
     1     MODE .EQ. 'linext' .OR. MODE .EQ. 'wcp' .OR.
     1     MODE .EQ. 'wce' .OR. MODE .EQ. 'lte' .OR.
     1     MODE .EQ. 'nanj' .OR. MODE .EQ. 'neg' .OR.
     1     MODE .EQ. 'bnd' .OR. MODE .EQ. 'bndp')THEN
          READ(SARG,*)NCF
        ELSE IF(MODE .EQ. 'qray')THEN
          READ(SARG,*)AMP
        ELSE
          READ(SARG,*)NU0
        END IF
      ELSE IF(MODE .EQ. 'qray')THEN
        AMP=0.025_LDP
      END IF
C tempu passes the uniform temperature through NU0.
      TTST=NU0
C
C Extended-domain modes: resolve the true NCF BEFORE allocating -
C the arrays must span the extended grid or the convolver walks
C out of bounds into the neighbouring allocation.
      NCF_P=NCF
      NEXT=0
      IF(MODE .EQ. 'bnd' .OR. MODE .EQ. 'linext' .OR.
     1   MODE .EQ. 'wce')THEN
        NEXT=NCF_P/4
        NCF=NCF_P+2*NEXT
      END IF
      IF(MODE .NE. 'linepair' .AND. MODE .NE. 'linext')THEN
        CALL GETARG(4,SARG)
        IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)DNU
        CALL GETARG(5,SARG)
        IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)AMP
      ELSE
        CALL GETARG(4,SARG)
        IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)NU0
        CALL GETARG(5,SARG)
        IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)DNU
        CALL GETARG(6,SARG)
        IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)AMP
      END IF
C
      IF(NCF .GT. NCFMX)STOP 'NCF>NCFMX'
      ALLOCATE(J_STORE(ND,NCF),J_IN(ND,NCF),RJ_CONT(ND,NCF),
     1     J_LINE(ND,NCF),NU(NCF))
C
      CHIBF=0.0_LDP
      CHIFF=0.0_LDP
      HDKT=4.7994145_LDP
      TWOHCSQ=0.0147452575_LDP
      OPLIN=0.0_LDP
      EMLIN=0.0_LDP
      L_FALSE=.FALSE.
      LUDUM=7
      CONT_REC=0
C
      T(1)=5.0_LDP
      T(2)=2.0_LDP
      T(3)=1.0_LDP
      T(4)=0.5_LDP
C
      IF(MODE .EQ. 'brk')THEN
C Bracket battery: exercises the production ES_PLANCK_NU_MATCH
C directly on synthetic log-space targets - the Planck path cannot
C reach the failure codes (see test header), so they are driven
C here through the same object code the convolver calls.
        CALL BRK_MATCH_BLOCK(OUT_DIR,LU_OUT)
        STOP
      END IF
C
C Decreasing log-spaced frequency grid (10^15 Hz), spanning the
C production range where e.s. redistribution acts.
C
      NU_HI=300.0_LDP
      NU_LO=0.02_LDP
      IF(MODE .EQ. 'lte')THEN
C Pathological grids for the match-failure probes: arg4/arg5
C override the frequency bounds (e.g. a grid inside the Planck
C peak triggers the below-range failure at the first Wien node).
        CALL GETARG(4,SARG)
        IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)NU_HI
        CALL GETARG(5,SARG)
        IF(LEN_TRIM(SARG) .GT. 0)READ(SARG,*)NU_LO
      END IF
      IF(MODE .EQ. 'bnd' .OR. MODE .EQ. 'linext' .OR.
     1   MODE .EQ. 'wce')THEN
C Extended domain: the production log spacing is continued NEXT
C points beyond BOTH ends, so every production node coincides
C exactly with an extended node.  The truncated/extended
C comparison then isolates the boundary treatment of the
C operator; interpolating between different spacings would
C measure resolution error on steep slopes instead.
        DLN=LOG(NU_HI/NU_LO)/(NCF_P-1.0_LDP)
        DO ML=1,NCF
          NU(ML)=NU_HI*EXP((NEXT-(ML-1.0_LDP))*DLN)
        END DO
      ELSE
        DO ML=1,NCF
          NU(ML)=EXP( LOG(NU_HI)+(LOG(NU_LO)-LOG(NU_HI))*
     1            (ML-1.0_LDP)/(NCF-1.0_LDP) )
        END DO
      END IF
C
      J_STORE(:,:)=0.0_LDP
      IF(MODE .EQ. 'planck' .OR. MODE .EQ. 'grid' .OR.
     1   MODE .EQ. 'lte' .OR. MODE .EQ. 'nanj')THEN
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=BB(T(K),NU(ML))
          END DO
        END DO
        IF(MODE .EQ. 'nanj')THEN
C Non-finite input probe: one NaN mid-grid at depth 2 leaves the
C other depths untouched, so the same run shows the diagnosed skip
C and the normal path side by side.
          J_STORE(2,NCF/2)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
        END IF
      ELSE IF(MODE .EQ. 'const')THEN
        J_STORE(:,:)=1.0_LDP
      ELSE IF(MODE .EQ. 'line' .OR. MODE .EQ. 'edge')THEN
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=BB(T(K),NU(ML))
          END DO
        END DO
        IF(MODE .EQ. 'line')THEN
          DO K=1,ND
            DO ML=1,NCF
              J_STORE(K,ML)=J_STORE(K,ML)+
     1          AMP*BB(T(K),NU0)*
     2          EXP(-0.5_LDP*((NU(ML)-NU0)/DNU)**2)
            END DO
          END DO
        ELSE
          DO K=1,ND
            DO ML=1,NCF
              J_STORE(K,ML)=J_STORE(K,ML)+
     1          AMP*BB(T(K),NU(2))*
     2          EXP(-0.5_LDP*((NU(ML)-NU(2))/DNU)**2)+
     3          AMP*BB(T(K),NU(NCF-1))*
     4          EXP(-0.5_LDP*((NU(ML)-NU(NCF-1))/DNU)**2)
            END DO
          END DO
        END IF
      ELSE IF(MODE .EQ. 'qray')THEN
C Non-Planckian J: a diluted blackbody plus a narrow emission
C feature at the CHI_RAY band centre gives a real RJ_ES-RJ
C contrast - the convention difference CHI_RAY*(RJ_ES-RJ)
C vanishes identically when J is already Planckian, and a
C feature narrower than the kernel width is what e.s.
C redistribution actually smears (the SOB-wing motivation).
C RJ_CONT keeps the unbumped continuum for the ETA definition.
        DO K=1,ND
          DO ML=1,NCF
            RJ_CONT(K,ML)=(0.15_LDP+0.85_LDP*(K-1.0_LDP)/
     1        (ND-1.0_LDP))*BB(T(3),NU(ML))
            J_STORE(K,ML)=RJ_CONT(K,ML) +
     2        50.0_LDP*BB(T(3),6.0_LDP)*
     3        EXP(-0.5_LDP*((NU(ML)-6.0_LDP)/0.5_LDP)**2)*
     4        (1.0_LDP-0.5_LDP*(K-1.0_LDP)/(ND-1.0_LDP))
          END DO
        END DO
      ELSE IF(MODE .EQ. 'tempu')THEN
C Uniform temperature at every depth: isolates the local
C beta(T)=1.84E-3*sqrt(T) width dependence from the
C Planck-preserving correction (the old 'temp' mode mixed
C different T per depth with a single B(T3) field).
        DO K=1,ND
          T(K)=TTST
          DO ML=1,NCF
            J_STORE(K,ML)=BB(TTST,NU(ML))
          END DO
        END DO
      ELSE IF(MODE .EQ. 'bnd' .OR. MODE .EQ. 'bndp')THEN
C Diluted continuum plus features that stay ALIVE at the
C production frequency bounds (NU=0.02 and NU=300): each
C Gaussian sits ~1.7 sigma inside the end node so the field
C and its ln(nu) curvature are significant where the
C one-sided endpoint stencils act. Absolute amplitudes are
C used - the Planck continuum is Wien-dead at both ends for
C T=1e4 K, so scaling to B there would zero the features.
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=(0.3_LDP+0.5_LDP*(K-1.0_LDP)/
     1        (ND-1.0_LDP))*BB(T(3),NU(ML)) +
     2        0.5_LDP*
     3        EXP(-0.5_LDP*((NU(ML)-0.045_LDP)/0.015_LDP)**2)+
     4        0.5_LDP*
     5        EXP(-0.5_LDP*((NU(ML)-260.0_LDP)/25.0_LDP)**2)
          END DO
        END DO
      ELSE IF(MODE .EQ. 'wcp' .OR. MODE .EQ. 'wce')THEN
C Smooth non-Planck continuum alive across the whole Wien-dead
C upper end: a NU**0.5 ramp with depth-varying normalisation on a
C weak Planck floor.  Tests whether the upper-boundary
C suppression is specific to the narrow bnd feature or applies
C to the entire non-thermal field.
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=(0.2_LDP+0.6_LDP*(K-1.0_LDP)/
     1        (ND-1.0_LDP))*(NU(ML)/NU_HI)**0.5_LDP+
     2        0.01_LDP*BB(T(K),NU(ML))
          END DO
        END DO
      ELSE IF(MODE .EQ. 'swp')THEN
C Feature-position sweep: one absolute-amplitude Gaussian moved
C progressively inward from the upper endpoint, on a weak Planck
C floor.  Same field at every depth - the per-depth T rows give
C the temperature response for free.
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=0.05_LDP*BB(T(3),NU(ML))+
     1        AMP*EXP(-0.5_LDP*((NU(ML)-NU0)/DNU)**2)
          END DO
        END DO
      ELSE IF(MODE .EQ. 'linepair' .OR. MODE .EQ. 'linext')THEN
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=BB(T(K),NU(ML))
          END DO
        END DO
      ELSE IF(MODE .EQ. 'pos')THEN
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=BB(T(K),NU(ML))
          END DO
        END DO
        J_STORE(:,NCF/3:2*NCF/5)=0.0_LDP
      ELSE IF(MODE .EQ. 'neg')THEN
C Signed input probe: the convolved J_ES picks up negative nodes,
C so the MINVAL(J_ES) gate must skip the Planck block before the
C A(:)=LOG(J_ES) statement can produce NaN.
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=BB(T(K),NU(ML))
          END DO
        END DO
        J_STORE(:,NCF/2-2:NCF/2+2)=-1.0_LDP
      ELSE IF(MODE .EQ. 'temp')THEN
        DO K=1,ND
          DO ML=1,NCF
            J_STORE(K,ML)=BB(T(3),NU(ML))
          END DO
        END DO
      ELSE
        STOP 'unknown mode'
      END IF
C
      IF(MODE .EQ. 'grid' .OR. MODE .EQ. 'linepair' .OR.
     1   MODE .EQ. 'wcp' .OR. MODE .EQ. 'bndp')THEN
        WRITE(FNAME,'(A,I5.5)')TRIM(MODE),NCF
      ELSE IF(MODE .EQ. 'bnd' .OR. MODE .EQ. 'linext' .OR.
     1   MODE .EQ. 'wce')THEN
        WRITE(FNAME,'(A,I5.5)')TRIM(MODE),NCF_P
      ELSE IF(MODE .EQ. 'lte')THEN
C Grid-tagged name so the pathological-grid probes do not collide
C with the default-grid diagnostic in the same work directory.
        WRITE(FNAME,'(A,I4.4,A,I4.4)')'lte',NINT(NU_HI),'x',
     1    NINT(100.0_LDP*NU_LO)
      ELSE IF(MODE .EQ. 'swp')THEN
        WRITE(FNAME,'(A,I5.5)')TRIM(MODE),NINT(10.0_LDP*NU0)
      ELSE
        FNAME=MODE
      END IF
C
      IF(MODE .EQ. 'linepair' .OR. MODE .EQ. 'linext')THEN
C Baseline-subtracted smearing measurement: convolve the plain
C continuum and the continuum+feature separately through the
C real operator so the excess (CONV(J_line)-CONV(J_cont)) can be
C measured against the input excess.
        J_IN(:,:)=J_STORE(:,:)
        CALL COMP_J_CONV_V2(J_STORE,NCF*ND,NU,T,ND,NCF,
     1       LUDUM,'J PASSED VIA CALL',CONT_REC,L_FALSE,L_FALSE,
     2       LUDUM,'RETURN J VIA CALL')
        RJ_CONT(:,:)=J_STORE(:,:)
        DO K=1,ND
          DO ML=1,NCF
            J_LINE(K,ML)=J_IN(K,ML)+
     1        AMP*BB(T(3),NU0)*
     2        EXP(-0.5_LDP*((NU(ML)-NU0)/DNU)**2)
          END DO
        END DO
        J_STORE(:,:)=J_LINE(:,:)
        CALL COMP_J_CONV_V2(J_STORE,NCF*ND,NU,T,ND,NCF,
     1       LUDUM,'J PASSED VIA CALL',CONT_REC,L_FALSE,L_FALSE,
     2       LUDUM,'RETURN J VIA CALL')
        OPEN(LU_OUT,FILE=TRIM(OUT_DIR)//'/esconv_'//
     1       TRIM(FNAME)//'.txt',STATUS='UNKNOWN')
        WRITE(LU_OUT,'(A,I3,A,I5)')'# ND=',ND,' NCF=',NCF
        WRITE(LU_OUT,'(A)')'# ML NU T J_CONT JBAR_CONT'//
     1    ' J_LINE JBAR_LINE'
        DO K=1,ND
          DO ML=1,NCF
            WRITE(LU_OUT,'(I5,ES27.17E3,ES15.7,4(ES27.17E3))')
     1        ML,NU(ML),T(K),J_IN(K,ML),RJ_CONT(K,ML),
     2        J_LINE(K,ML),J_STORE(K,ML)
          END DO
        END DO
        CLOSE(LU_OUT)
        STOP
      END IF
C
      J_IN(:,:)=J_STORE(:,:)
C
      CALL COMP_J_CONV_V2(J_STORE,NCF*ND,NU,T,ND,NCF,
     1     LUDUM,'J PASSED VIA CALL',CONT_REC,L_FALSE,L_FALSE,
     2     LUDUM,'RETURN J VIA CALL')
C
      IF(MODE .EQ. 'qray')THEN
        CALL QRAY_BLOCK(OUT_DIR,NU,T,J_IN,J_STORE,RJ_CONT,
     1     ND,NCF,LU_OUT,AMP)
      ELSE IF(MODE .EQ. 'lte')THEN
        CALL LTE_MAP_BLOCK(OUT_DIR,NU,T,J_IN,J_STORE,
     1     ND,NCF,LU_OUT,FNAME)
      ELSE
        OPEN(LU_OUT,FILE=TRIM(OUT_DIR)//'/esconv_'//
     1       TRIM(FNAME)//'.txt',STATUS='UNKNOWN')
        WRITE(LU_OUT,'(A,I3,A,I5)')'# ND=',ND,' NCF=',NCF
        WRITE(LU_OUT,'(A)')'# ML NU T J_IN JBAR'
        DO K=1,ND
          DO ML=1,NCF
            WRITE(LU_OUT,'(I5,ES27.17E3,ES15.7,2(ES27.17E3))')
     1        ML,NU(ML),T(K),J_IN(K,ML),J_STORE(K,ML)
          END DO
        END DO
        CLOSE(LU_OUT)
      END IF
C
      STOP
      END
C
C*******************************************************************
C
C Q-RAY block: manufacture an INC_RAY=T atmosphere and emit the two
C production emissivity conventions elementwise.
C
C   TA_model = ETA_CLUMP + CHI_SCAT_CLUMP*RJ_ES   (comp_j_blank:843)
C   TA_obs   = ETA_CLUMP + ESEC_CLUMP*RJ_ES
C             + CHI_RAY_CLUMP*RJ                  (COMP_JCONT:751)
C
C CHI_RAY uses the production Rayleigh shape - neutral-H Rayleigh
C scattering rises steeply with frequency; we take chi_ray ∝ NU**4
C normalized per depth (the standard H I Rayleigh form, as produced
C by RAYLEIGH_SCAT below the Lyman edges). ESEC is uniform; ETA is
C a thermal (absorptive) emissivity; CLUMP_FAC is applied to all
C opacity/emissivity factors as in production.
C
      SUBROUTINE QRAY_BLOCK(OUT_DIR,NU,T,RJ,RJ_ES,RJ_CONT,
     1     ND,NCF,LU,RRAY0)
      USE SET_KIND_MODULE
      IMPLICIT NONE
      INTEGER ND,NCF,LU
      REAL(KIND=LDP) NU(NCF),T(ND),RJ(ND,NCF),RJ_ES(ND,NCF)
      REAL(KIND=LDP) RJ_CONT(ND,NCF),RRAY0
      CHARACTER(LEN=*) OUT_DIR
C
      REAL(KIND=LDP) CHI_RAY(ND,NCF),ESEC(ND),ETA(ND,NCF)
      REAL(KIND=LDP) CHI(ND,NCF),CLUMP(ND)
      REAL(KIND=LDP) TA_MODEL(ND,NCF),TA_OBS(ND,NCF)
      INTEGER I,ML
C
C Rayleigh opacity rising as NU**4, weighted down with depth so it
C resembles a recombined outer atmosphere. Amplitude chosen so the
C max chi_ray/esec ratio is O(1) at the outer edge - a deliberately
C strong but plausible case for bounding the convention difference.
C
      DO I=1,ND
        ESEC(I)=0.05_LDP*I                    ! rises inward
        CLUMP(I)=1.0_LDP+0.05_LDP*(ND-I)      ! mild clumping outward
        DO ML=1,NCF
          CHI_RAY(I,ML)=RRAY0*ESEC(I)*
     1      EXP(-0.5_LDP*((NU(ML)-6.0_LDP)/4.0_LDP)**2)*
     2      (1.0_LDP+2.0_LDP*(I-1)/(ND-1.0_LDP))   ! stronger inward
          ETA(I,ML)=0.01_LDP*ESEC(I)*RJ_CONT(I,ML)
          CHI(I,ML)=ESEC(I)+CHI_RAY(I,ML)+0.01_LDP*I
        END DO
      END DO
C
      DO I=1,ND
        DO ML=1,NCF
          TA_MODEL(I,ML)=ETA(I,ML)*CLUMP(I)+
     1      (ESEC(I)+CHI_RAY(I,ML))*CLUMP(I)*RJ_ES(I,ML)
          TA_OBS(I,ML)=ETA(I,ML)*CLUMP(I)+
     1      ESEC(I)*CLUMP(I)*RJ_ES(I,ML)+
     2      CHI_RAY(I,ML)*CLUMP(I)*RJ(I,ML)
        END DO
      END DO
C
      OPEN(LU,FILE=TRIM(OUT_DIR)//'/esconv_qray_fields.txt',
     1     STATUS='UNKNOWN')
      WRITE(LU,'(A)')'# I ML NU RJ RJ_ES CHI_RAY ESEC ETA CHI'//
     1  ' TA_MODEL TA_OBS CLUMP'
      DO I=1,ND
        DO ML=1,NCF
          WRITE(LU,'(I3,I5,10(ES27.17E3))')
     1      I,ML,NU(ML),RJ(I,ML),RJ_ES(I,ML),CHI_RAY(I,ML),
     2      ESEC(I),ETA(I,ML),CHI(I,ML),TA_MODEL(I,ML),
     3      TA_OBS(I,ML),CLUMP(I)
        END DO
      END DO
      CLOSE(LU)
      RETURN
      END
C
C*******************************************************************
C
C LTE endpoint diagnostic (P8.1 robustness audit): replicates the
C convolved Planck field B_es through the production THOMAS_RH and
C the bounded PLANCK_NU map so the per-node treatment zone can be
C measured against the production JBAR.
C
C The replication is diagnostic only: it mirrors the production
C formulas (same BB scaling with the T3 cap, same stencil
C coefficients, same bounded bracket search and continuation) so the
C zone flags describe what the production routine actually did.
C
C ZONE codes: 0 corrected interior, 1 upper endpoint, 2 invalid-map
C bracket (I=1 during the shift interpolation), 3 gate-skip
C (MINVAL block not entered), 4 continuation tail (bracket beyond
C the last derived node), 5 no-valid-map bypass, 9 lower endpoint.
C
      SUBROUTINE LTE_MAP_BLOCK(OUT_DIR,NU,T,JIN,JBAR,ND,NCF,LU,FNAME)
      USE SET_KIND_MODULE
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE
      INTEGER ND,NCF,LU
      REAL(KIND=LDP) NU(NCF),T(ND),JIN(ND,NCF),JBAR(ND,NCF)
      CHARACTER(LEN=*) OUT_DIR,FNAME
C
      COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
      REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
C
      INTEGER NCOEF,IONE
      PARAMETER (NCOEF=2,IONE=1)
      REAL(KIND=LDP) ACOEF(2),BCOEF(2)
      DATA ACOEF/1.690703717290_LDP,-0.690703717290_LDP/
      DATA BCOEF/1.614249968779_LDP,2.154326524957_LDP/
C
      REAL(KIND=LDP) A_ST(NCF),C_ST(NCF)
      REAL(KIND=LDP) A(NCF),H(NCF),C(NCF),D(NCF)
      REAL(KIND=LDP) B(NCF),BES(NCF)
      REAL(KIND=LDP) PES(NCF),PFN(NCF),PNU(NCF)
      INTEGER ZN(NCF)
      REAL(KIND=LDP) BETA,T1,T2,T3,D1,D2,DH
      INTEGER I,K,L,ML,FAIL,LASTOK
      LOGICAL MAPOK
C
      A_ST(1)=0.0_LDP
      C_ST(1)=-2.0_LDP/( LOG(NU(1)/NU(2)) )**2
      DO ML=2,NCF-1
        D1=LOG(NU(ML-1)/NU(ML))
        D2=LOG(NU(ML)/NU(ML+1))
        DH=0.5_LDP*(D1+D2)
        A_ST(ML)=-1.0_LDP/D1/DH
        C_ST(ML)=-1.0_LDP/D2/DH
      END DO
      A_ST(NCF)=-2.0_LDP/( LOG(NU(NCF-1)/NU(NCF)) )**2
      C_ST(NCF)=0.0_LDP
C
      OPEN(LU,FILE=TRIM(OUT_DIR)//'/esconv_'//TRIM(FNAME)//'.txt',
     1     STATUS='UNKNOWN')
      WRITE(LU,'(A,I3,A,I5)')'# ND=',ND,' NCF=',NCF
      WRITE(LU,'(A)')'# ML NU T B B_ES PNU ZONE FAIL JBAR'
      DO K=1,ND
        BETA=1.84E-03_LDP*SQRT(T(K))
        T3=0.5_LDP*HDKT*NU(1)/T(K)
        IF(T3 .LT. 1.0_LDP)T3=0.0_LDP
        IF(T3 .GT. 700.0_LDP)T3=700.0_LDP
        DO ML=1,NCF
          T1=EXP( -HDKT*NU(ML)/T(K) )
          T2=EXP( T3-HDKT*NU(ML)/T(K) )
          B(ML)=TWOHCSQ*(NU(ML)**3)*T2/(1.0_LDP-T1)
        END DO
        BES(:)=0.0_LDP
        DO L=1,NCOEF
          T1=BETA*BETA/BCOEF(L)/BCOEF(L)
          A(:)=T1*A_ST(:)
          H(:)=-1.0_LDP
          C(:)=T1*C_ST(:)
          D(:)=B(:)
          CALL THOMAS_RH(A,H,C,D,NCF,IONE)
          BES(:)=BES(:)+ACOEF(L)*D(:)
        END DO
C
        PNU(:)=0.0_LDP
        ZN(:)=3
        FAIL=-1
        IF(MINVAL(B) .GT. 0.0_LDP .AND. MINVAL(BES) .GT. 0.0_LDP)THEN
          PES(:)=LOG(BES(:))
          PFN(:)=LOG(B(:))
          PNU(1)=NU(1)
          I=2
          ML=1
          LASTOK=1
          MAPOK=.TRUE.
          DO WHILE(I .LE. NCF)
            IF(NU(I) .LE. 1.5_LDP*T(K))EXIT
C The replica calls the same production match the convolver does,
C so the zone map can never drift from the shipped search.
            CALL ES_PLANCK_NU_MATCH(PES(I),PFN,ML,NCF,FAIL,T1)
            IF(FAIL .NE. 0)THEN
              IF(I .LE. 2)MAPOK=.FALSE.
              EXIT
            END IF
            PNU(I)=T1*NU(ML-1)+(1.0_LDP-T1)*NU(ML)
            LASTOK=I
            FAIL=-1
            I=I+1
          END DO
          IF(MAPOK)THEN
            T1=1.0_LDP/NU(I-1)-1.0_LDP/PNU(I-1)
            DO ML=I,NCF
              PNU(ML)=1.0_LDP/(1.0_LDP/NU(ML)-T1)
            END DO
            ZN(1)=1
            ZN(NCF)=9
            I=1
            DO ML=2,NCF-1
              DO WHILE(NU(ML) .LT. PNU(I+1) .AND. I .LT. NCF-1)
                I=I+1
              END DO
              IF(I .EQ. 1)THEN
                ZN(ML)=2
              ELSE IF(I+1 .GT. LASTOK)THEN
                ZN(ML)=4
              ELSE
                ZN(ML)=0
              END IF
            END DO
          ELSE
            ZN(:)=5
          END IF
        END IF
        DO ML=1,NCF
          WRITE(LU,'(I5,ES27.17E3,ES15.7,3(ES27.17E3),2I3,
     1      ES27.17E3)')ML,NU(ML),T(K),JIN(K,ML),BES(ML),
     2      PNU(ML),ZN(ML),FAIL,JBAR(K,ML)
        END DO
      END DO
      CLOSE(LU)
      RETURN
      END
C
C
C
      SUBROUTINE BRK_MATCH_BLOCK(OUT_DIR,LU)
      USE SET_KIND_MODULE
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE
C
C Synthetic battery for ES_PLANCK_NU_MATCH.  The convolver feeds the
C match only log-space Planck targets, and on the Wien-gated region
C PLANCK_FN is monotone and PLANCK_ES is a kernel average, so the
C failure codes cannot be produced through the public call - they
C are driven here on the same compiled code the convolver uses.
C
C Columns: ICASE, ML_IN, PES, FAIL, ML_OUT, T1, PNU
C ICASE list:
C  1  valid interior bracket, fresh walk (ML=1)
C  2  target above every PFN (fresh walk)         -> FAIL=1, ML=NCF+1
C  3  target below PFN(1) at first node           -> FAIL=2, ML=1
C  4  persistent ML floor walk off the top        -> FAIL=2, ML=NCF+1
C  5  degenerate PFN(ML-1)=PFN(ML) bracket        -> FAIL=3
C  6  NaN target                                 -> FAIL=4
C  7  +Inf target                                -> FAIL=4
C  8  huge-negative target (log of a tiny Planck) -> FAIL=2
C  9  target exactly on a PFN node               -> FAIL=0, T1=0
C 10  continuation from a persisted ML            -> FAIL=0
C
      CHARACTER(LEN=*) OUT_DIR
      INTEGER LU
      INTEGER, PARAMETER :: NB=8
      REAL(KIND=LDP) NU_B(NB),PFN_B(NB),PFN_D(NB)
      REAL(KIND=LDP) PES,T1,PNU
      INTEGER ML,FAIL,I
C
      DO I=1,NB
        NU_B(I)=NB+1.0_LDP-I
        PFN_B(I)=2.0_LDP*I-10.0_LDP
        PFN_D(I)=PFN_B(I)
      END DO
      PFN_D(4)=PFN_D(3)
C
      OPEN(LU,FILE=TRIM(OUT_DIR)//'/esconv_brk.txt',
     1     STATUS='UNKNOWN')
      WRITE(LU,'(A)')'# ICASE ML_IN PES FAIL ML_OUT T1 PNU'
C
      ML=1
      PES=-5.0_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      PNU=0.0_LDP
      IF(FAIL .EQ. 0)PNU=T1*NU_B(ML-1)+(1.0_LDP-T1)*NU_B(ML)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')1,1,PES,FAIL,ML,T1,PNU
C
      ML=1
      PES=10.0_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')2,1,PES,FAIL,ML,T1,0.0_LDP
C
      ML=1
      PES=-20.0_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')3,1,PES,FAIL,ML,T1,0.0_LDP
C
      ML=5
      PES=-20.0_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')4,5,PES,FAIL,ML,T1,0.0_LDP
C
      ML=4
      PES=-4.0_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_D,ML,NB,FAIL,T1)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')5,4,PES,FAIL,ML,T1,0.0_LDP
C
      ML=1
      PES=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')6,1,PES,FAIL,ML,T1,0.0_LDP
C
      ML=1
      PES=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')7,1,PES,FAIL,ML,T1,0.0_LDP
C
      ML=1
      PES=-1.0E300_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')8,1,PES,FAIL,ML,T1,0.0_LDP
C
      ML=1
      PES=0.0_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      PNU=0.0_LDP
      IF(FAIL .EQ. 0)PNU=T1*NU_B(ML-1)+(1.0_LDP-T1)*NU_B(ML)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')9,1,PES,FAIL,ML,T1,PNU
C
      ML=3
      PES=-1.0_LDP
      CALL ES_PLANCK_NU_MATCH(PES,PFN_B,ML,NB,FAIL,T1)
      PNU=0.0_LDP
      IF(FAIL .EQ. 0)PNU=T1*NU_B(ML-1)+(1.0_LDP-T1)*NU_B(ML)
      WRITE(LU,'(I3,I4,ES23.14E3,I4,I4,2ES15.7)')10,3,PES,FAIL,ML,T1,PNU
C
      CLOSE(LU)
      RETURN
      END
