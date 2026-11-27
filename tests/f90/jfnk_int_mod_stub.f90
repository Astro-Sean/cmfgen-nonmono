! Stub for jfnk_cmfgen_int_mod: provides only the symbols that
! jfnk_solver_v1.f90 imports, so the production Newton/GMRES code
! can be unit-tested without the full CMFGEN atmosphere state.
!
module jfnk_cmfgen_int_mod
   implicit none
   save
   integer, save :: JFNK_NT = 0
   integer, save :: JFNK_ND = 0
   logical, save :: JFNK_COUPLED_LOC = .false.
   double precision, allocatable, save :: JFNK_RSCALE(:)
   logical, save :: JFNK_RSCALE_READY = .false.
   logical, save :: JFNK_BA_BUILT = .false.
   double precision, allocatable, save :: JFNK_BA_DIAG(:,:)
   logical, save :: JFNK_BA_DIAG_READY = .false.
   double precision, allocatable, save :: JFNK_BA_DIAG_RAW(:,:)
   logical, save :: JFNK_BA_DIAG_RAW_READY = .false.
   double precision, allocatable, save :: JFNK_BA_A(:,:,:)
   logical, save :: JFNK_LU_READY = .false.
   double precision, allocatable, save :: JFNK_BA_RIDGE_ADD(:,:)
   ! WP-0: ledger/accounting routines report the interface-side eval
   ! total; the stub counter stays 0 unless a driver's f() bumps it.
   integer, save :: JFNK_N_RESIDUAL_EVALS = 0
   ! WP-3: mirrors production -- solver checks the flag after each
   ! eval; tests set it directly to exercise the graceful stop.
   logical, save :: JFNK_EVAL_BUDGET_HIT = .false.
   integer, save :: JFNK_EVAL_CAP = 0
   ! WP-2: counts graded-stagnation refresh calls so tests can see
   ! the consumption on the rejection path.
   integer, save :: JFNK_STUB_NFACTOR = 0
end module jfnk_cmfgen_int_mod

! Production JFNK_FACTOR_PRECOND_V1 is an external routine (not a
! module member), so the stub mirrors that shape.
subroutine JFNK_FACTOR_PRECOND_V1(do_factor)
   use jfnk_cmfgen_int_mod, only: JFNK_STUB_NFACTOR
   implicit none
   logical, intent(in) :: do_factor
   if(do_factor) JFNK_STUB_NFACTOR = JFNK_STUB_NFACTOR + 1
end subroutine JFNK_FACTOR_PRECOND_V1
