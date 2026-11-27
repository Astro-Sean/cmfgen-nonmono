! Stubs for externals referenced by fg_j_cmf_nonmono_v1.f and
! mom_j_cmf_nonmono_v1.f that the RD_CONTROL_VARIABLES parse path
! never executes.  Only the SET_* entry points in those files are
! reached by the driver.
      SUBROUTINE KN_ORD_SPV()
      END
      SUBROUTINE KN_ORD_SPW_FILL()
      END
      SUBROUTINE SWEEP_FLD_ORD()
      END
      SUBROUTINE TUNE()
      END
      SUBROUTINE TUVGHD_RH()
      END
      SUBROUTINE WRITE_DIRECT_INFO_V3()
      END
C WP-2: JFNK_FACTOR_PRECOND_V1 is referenced (flag-gated) by
C newton_hook_v1 in jfnk_solver_v1.f90, which the vadat driver
C links; the parse path never calls it.
      SUBROUTINE JFNK_FACTOR_PRECOND_V1(DO_FACTOR)
      LOGICAL DO_FACTOR
      END
