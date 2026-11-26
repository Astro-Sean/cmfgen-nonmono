! Driver for the X < -700 overflow bound in
!   subs/exponx.f  (EXPONX, d_EXPONX_dX)   [SOB_check WP3]
!
! Usage: drv_exponx <casefile>
!   casefile: line 1 = N; then N rows, each a single X value.
! Output: one row per X as "i  X  EXPONX(X)  d_EXPONX_dX(X)", then a
! terminal "GUARD n_f n_d" row with the module activation counters.
! The reference build links stub_exponx_guard.f90 (reference exponx.f
! carries no bound and no module), so the same case file gives a
! direct new-vs-reference comparison.
!
program drv_exponx
  use set_kind_module
  use exponx_guard_mod
  implicit none

  character(len=512) :: cfile
  integer :: lu, ios, n, i
  real(kind=ldp) :: x, f, d
  real(kind=ldp) :: exponx, d_exponx_dx
  external exponx, d_exponx_dx

  call get_command_argument(1, cfile)
  open(newunit=lu, file=trim(cfile), status='old', action='read', &
       iostat=ios)
  if (ios /= 0) then
     write(*,*) 'cannot open ', trim(cfile)
     stop 2
  end if
  read(lu,*) n
  do i = 1, n
     read(lu,*) x
     f = exponx(x)
     d = d_exponx_dx(x)
     write(*,'(i6,3(1x,es24.16e3))') i, x, f, d
  end do
  close(lu)
  write(*,'(a,2(1x,i12))') 'GUARD ', n_exponx_guard, n_dexponx_guard

end program drv_exponx
