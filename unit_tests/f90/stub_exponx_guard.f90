module exponx_guard_mod
  ! Stub matching the module defined at the top of the jfnk-tree
  ! cmfgen/subs/exponx.f.  The reference exponx.f has no X < -700
  ! bound and no counters; this stub lets drv_exponx link against the
  ! reference object so identical inputs give a direct comparison.
  implicit none
  integer :: n_exponx_guard = 0
  integer :: n_dexponx_guard = 0
end module exponx_guard_mod
