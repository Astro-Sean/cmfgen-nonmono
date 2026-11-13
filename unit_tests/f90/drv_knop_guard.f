C
C Driver for the production RD_CMF_FLUX_CONTROLS validation path
C (P8.1 T6). Reads a manufactured control file into the option
C store exactly as cmf_flux does, then calls the real
C RD_CMF_FLUX_CONTROLS so the DO_NONMONO_KNOP ES-prerequisite
C guard at :654 fires (or not) on the production code path.
C
C Usage: drv_knop_guard <control_file>
C Prints 'GUARD_PASS' to stdout when the reader returns; a guard
C STOP exits non-zero with the named diagnostic on LUER.
C
C SET_FG_KNOP_DEBUG/SET_FG_FREQ_MATCH/SET_MOM_FREQ_MATCH are
C local stubs: they forward solver-module flags unrelated to the
C ES guard. Linking the full Knop and moment solver files to
C reach three one-line setters is not justified; the guarded
C path itself is entirely production code.
C
      PROGRAM DRV_KNOP_GUARD
      USE SET_KIND_MODULE
      USE MOD_CMF_OBS
      IMPLICIT NONE
      INTEGER ND,LUMOD,LUER,LU_IN,LU_OUT,ERROR_LU
      CHARACTER(LEN=200) CFILE
      EXTERNAL RD_OPTIONS_INTO_STORE,RD_CMF_FLUX_CONTROLS,CLEAN_RD_STORE
      EXTERNAL ERROR_LU
C
      CALL GETARG(1,CFILE)
      ND=4
      LUMOD=6
      LUER=ERROR_LU()
      NUM_IONS=0
C
      LU_IN=31
      OPEN(LU_IN,FILE=TRIM(CFILE),STATUS='OLD',ACTION='READ')
      LU_OUT=32
      OPEN(LU_OUT,STATUS='SCRATCH',FORM='FORMATTED')
      CALL RD_OPTIONS_INTO_STORE(LU_IN,LU_OUT)
C
      CALL RD_CMF_FLUX_CONTROLS(ND,LUMOD,LUER)
      WRITE(*,'(A)')'GUARD_PASS'
C CLEAN_RD_STORE is skipped: the verdict is already emitted and
C process teardown frees the store.
      STOP
      END
C
      SUBROUTINE SET_FG_KNOP_DEBUG(DO_DBG)
      LOGICAL DO_DBG
      RETURN
      END
C
      SUBROUTINE SET_FG_FREQ_MATCH(DO_MATCH)
      LOGICAL DO_MATCH
      RETURN
      END
C
      SUBROUTINE SET_MOM_FREQ_MATCH(DO_MATCH)
      LOGICAL DO_MATCH
      RETURN
      END
