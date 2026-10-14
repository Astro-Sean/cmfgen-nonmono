! Driver for RECOMPUTE_SIGMA_TRUE / SIGMA_TRUE_CORE
! (new_main/subs/recompute_sigma_true.f).
!
! Usage:
!   drv_sigma_true core   <rvfile>
!   drv_sigma_true recomp <rvfile> <off|sob|knop> [allow]
!
! rvfile: first line ND, then ND rows "R V SIGMA" (R outer -> inner).
! For "recomp" the SIGMA column seeds the work array that the SOB path
! re-clamps; for "core" the column is ignored but must be present.
! The optional 4th arg "allow" sets ALLOW_SIGMA_UNDEF=T, downgrading
! the N_UNDEF>0 abort to a warning.
!
program drv_sigma_true
  use set_kind_module
  use control_variable_mod, only: do_nonmono_sob, do_nonmono_knop, &
                                allow_sigma_undef
  implicit none

  character(len=64)  :: cmd, mode
  character(len=512) :: rfile
  integer :: nd, i, n_undef, lu, ios
  real(kind=ldp), allocatable :: r(:), v(:), sig(:), st(:)

  call get_command_argument(1, cmd)
  call get_command_argument(2, rfile)
  open(newunit=lu, file=trim(rfile), status='old', action='read', &
       iostat=ios)
  if (ios /= 0) then
     write(*,*) 'cannot open ', trim(rfile)
     stop 2
  end if
  read(lu,*) nd
  allocate(r(nd), v(nd), sig(nd), st(nd))
  do i = 1, nd
     read(lu,*) r(i), v(i), sig(i)
  end do
  close(lu)

  ! Sentinel: entries left untouched by the callee must show up in the
  ! output so the test can distinguish "computed 0" from "not written".
  st = -999.0_ldp

  if (cmd == 'core') then
     call sigma_true_core(st, r, v, nd, n_undef)
     write(*,'(a,i6)') 'N_UNDEF ', n_undef
     do i = 1, nd
        write(*,'(i6,3(1x,es24.16e3))') i, r(i), v(i), st(i)
     end do
  else if (cmd == 'recomp') then
     call get_command_argument(3, mode)
     do_nonmono_sob  = (mode == 'sob')
     do_nonmono_knop = (mode == 'knop')
     allow_sigma_undef = .false.
     if (command_argument_count() >= 4) then
        call get_command_argument(4, mode)
        allow_sigma_undef = (mode == 'allow')
     end if
     call recompute_sigma_true(st, sig, r, v, nd)
     do i = 1, nd
        write(*,'(i6,4(1x,es24.16e3))') i, r(i), v(i), sig(i), st(i)
     end do
  else
     write(*,*) 'unknown cmd ', trim(cmd)
     stop 2
  end if

end program drv_sigma_true
