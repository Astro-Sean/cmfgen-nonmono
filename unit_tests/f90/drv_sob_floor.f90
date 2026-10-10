! Driver for the Dessart+2015 tau-floor edits in
!   subs/sobjbar_sim.f   (SOBJBAR_SIM)
!   subs/sobew.f         (SOBEW)
!   subs/sobew_grad_v2.f (SOBEW_GRAD_V2)
!   subs/sobew_grad_v3.f (SOBEW_GRAD_V3)
!
! The same driver links against either the jfnk-tree source or the
! pristine reference source, so identical inputs give a direct
! new-vs-reference comparison.
!
! Usage: drv_sob_floor <routine> <casefile>
!   routine: sobjsim | sobew | sobewg2 | sobewg3
!
! casefile:
!   line1: ND NC NP NSIM
!   line2: FL DIF THICK DIE DBB IC RLUM
!   line3: METHOD INNER_BND_METH
!   ND rows:  R V SIGMA CHI ESEC CHIL ETAL SOURCE
!   NP rows:  P AQW HAQW
!   NSIM rows (sobjsim only): CHIL_MAT ETAL_MAT BBCOR  (uniform cols)
!
! AQW(I,LS) = w(LS) for all I (ray weights constant with depth), and
! only rays with LS <= NC are used (NC=NP -> every ray sees all ND
! depths), so BETA(I) = SUM_LS w(LS)*EXF(GAM(I,LS)) is exactly
! predictable from the floored tau formula.
!
program drv_sob_floor
  use set_kind_module
  implicit none

  integer, parameter :: NDMAX = 400, NPMAX = 40, NSMAX = 8

  character(len=32)  :: routine
  character(len=512) :: cfile
  character(len=16)  :: method, ibnd
  integer :: lu, ios, nd, nc, np, nsim, i, l
  real(kind=ldp) :: fl, dbb, xic, rlum
  logical :: dif, thick, die

  real(kind=ldp) :: r(ndmax), v(ndmax), sig(ndmax)
  real(kind=ldp) :: chi(ndmax), esec(ndmax), chil(ndmax), etal(ndmax)
  real(kind=ldp) :: src(ndmax), p(npmax), aqw(ndmax, npmax)
  real(kind=ldp) :: haqw(ndmax, npmax)
  real(kind=ldp) :: chilm(ndmax, nsmax), etalm(ndmax, nsmax)
  real(kind=ldp) :: bbcor(ndmax, nsmax)
  real(kind=ldp) :: jbar(ndmax), beta(ndmax), znet(ndmax, nsmax)
  real(kind=ldp) :: vb(ndmax, nsmax), vc(ndmax, nsmax)
  real(kind=ldp) :: vb2(ndmax), vc2(ndmax), betac(ndmax, nsmax)
  real(kind=ldp) :: line_flux(ndmax), fm(ndmax), ew, emiss_ew, cflux

  call get_command_argument(1, routine)
  call get_command_argument(2, cfile)
  open(newunit=lu, file=trim(cfile), status='old', action='read', &
       iostat=ios)
  if (ios /= 0) then
     write(*,*) 'cannot open ', trim(cfile)
     stop 2
  end if
  read(lu,*) nd, nc, np, nsim
  read(lu,*) fl, dif, thick, die, dbb, xic, rlum
  read(lu,*) method, ibnd
  do i = 1, nd
     read(lu,*) r(i), v(i), sig(i), chi(i), esec(i), chil(i), &
                etal(i), src(i)
  end do
  aqw = 0.0_ldp
  haqw = 0.0_ldp
  do l = 1, np
     read(lu,*) p(l), aqw(1, l), haqw(1, l)
     aqw(:, l) = aqw(1, l)
     haqw(:, l) = haqw(1, l)
  end do
  chilm = 0.0_ldp
  etalm = 1.0_ldp
  bbcor = 1.0_ldp
  do l = 1, nsim
     read(lu,*,iostat=ios) chilm(1, l), etalm(1, l), bbcor(1, l)
     chilm(:, l) = chilm(1, l)
     etalm(:, l) = etalm(1, l)
     bbcor(:, l) = bbcor(1, l)
  end do
  close(lu)

  if (routine == 'sobjsim') then
     call sobjbar_sim(src, chi, chil, etal, v, sig, r, p, aqw, &
                      jbar, beta, znet, vb, vc, vb2, vc2, betac, &
                      chilm, etalm, bbcor, &
                      fl, dif, dbb, xic, thick, &
                      nd, nc, np, nd, nsim, method)
     do i = 1, nd
        write(*,'(i6,4(1x,es24.16e3))') i, beta(i), jbar(i), znet(i,1), vb2(i)
     end do
  else if (routine == 'sobew') then
     call sobew(src, chi, esec, chil, etal, v, sig, r, p, &
                aqw, haqw, line_flux, ew, cflux, &
                fl, dif, dbb, xic, thick, die, nc, np, nd, method)
     write(*,'(a,2(1x,es24.16e3))') 'EW_CF ', ew, cflux
     do i = 1, nd
        write(*,'(i6,1x,es24.16e3)') i, line_flux(i)
     end do
  else if (routine == 'sobewg2') then
     call sobew_grad_v2(src, chi, esec, chil, etal, v, sig, r, p, &
                fm, rlum, aqw, haqw, line_flux, ew, cflux, &
                fl, ibnd, dbb, xic, thick, die, nc, np, nd, method)
     write(*,'(a,2(1x,es24.16e3))') 'EW_CF ', ew, cflux
     do i = 1, nd
        write(*,'(i6,2(1x,es24.16e3))') i, fm(i), line_flux(i)
     end do
  else if (routine == 'sobewg3') then
     call sobew_grad_v3(src, chi, esec, chil, etal, v, sig, r, p, &
                fm, rlum, aqw, haqw, line_flux, ew, emiss_ew, cflux, &
                fl, ibnd, dbb, xic, thick, die, nc, np, nd, method)
     write(*,'(a,3(1x,es24.16e3))') 'EW_CF ', ew, emiss_ew, cflux
     do i = 1, nd
        write(*,'(i6,2(1x,es24.16e3))') i, fm(i), line_flux(i)
     end do
  else
     write(*,*) 'unknown routine ', trim(routine)
     stop 2
  end if

end program drv_sob_floor
