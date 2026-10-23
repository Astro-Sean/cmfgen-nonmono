!----------------------------------------------------------------------
! JFNK_INT_STUB: stand-in for jfnk_cmfgen_int_mod for the standalone
! Newton-core test.  Provides exactly the symbols that
! new_main/mod_subs/jfnk_solver_v1.f90 imports; the production
! interface module pulls in the whole CMFGEN model state, which the
! test does not need.
! edited by Devin, 06-Oct-2026
!----------------------------------------------------------------------
module jfnk_cmfgen_int_mod
   implicit none
   save
   integer :: JFNK_NT = 3
   integer :: JFNK_ND = 0
   logical :: JFNK_COUPLED_LOC = .false.
   double precision, allocatable :: JFNK_RSCALE(:)
   logical :: JFNK_RSCALE_READY = .false.
   ! Only referenced by the diagnostic routines (never called by the
   ! test); kept so the `use ... only:` lists resolve.
   logical :: JFNK_BA_BUILT = .false.
   double precision, allocatable :: JFNK_BA_DIAG(:,:)
   logical :: JFNK_BA_DIAG_READY = .false.
   double precision, allocatable :: JFNK_BA_A(:,:,:)
   logical :: JFNK_LU_READY = .false.
   double precision, allocatable :: JFNK_BA_RIDGE_ADD(:,:)
   ! Budget-sentinel plumbing added by the commit-verification fix
   ! (21b52ca): the stub mirrors the production module so the G11
   ! contract test exercises the same symbols.
   logical :: JFNK_EVAL_BYPASS = .false.
   integer :: JFNK_MAX_EVALS = 0
   integer :: JFNK_N_RESIDUAL_EVALS = 0
   ! Activation-local budget the driver trims per rung; the residual
   ! sentinel enforces the tighter of it and JFNK_MAX_EVALS.
   integer :: JFNK_EVAL_LIM = 0
end module jfnk_cmfgen_int_mod
