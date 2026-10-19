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
   double precision, allocatable, save :: JFNK_BA_A(:,:,:)
   logical, save :: JFNK_LU_READY = .false.
   double precision, allocatable, save :: JFNK_BA_RIDGE_ADD(:,:)
end module jfnk_cmfgen_int_mod
