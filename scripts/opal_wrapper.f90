! Thin F90 wrapper around the SNEC OPAL opacity interpolator.
! Exposes opal_grid() callable from Python via f2py.
!
! opal_grid: evaluate log10(kappa_Ross) on a flat array of (t6, r) points
!            for a single fixed composition.
!
! Inputs (all real*4 scalars unless noted):
!   z           metallicity (e.g. 0.02)
!   xh          hydrogen mass fraction
!   xxc         enhanced carbon mass fraction (delta Xc)
!   xxo         enhanced oxygen mass fraction (delta Xo)
!   t6arr(n)    temperatures in units of 1e6 K
!   rarr(n)     OPAL R = rho / t6^3
!   n           number of grid points
!   table_dir   path to directory CONTAINING the tables/ subdirectory (len=512)
!
! Output:
!   log_kappa_arr(n)  log10(Rosseland mean opacity [cm^2/g])

subroutine opal_grid(z, xh, xxc, xxo, t6arr, rarr, n, log_kappa_arr, table_dir)
  implicit none
  integer,       intent(in)  :: n
  real*4,        intent(in)  :: z, xh, xxc, xxo
  real*4,        intent(in)  :: t6arr(n), rarr(n)
  real*4,        intent(out) :: log_kappa_arr(n)
  character*512, intent(in)  :: table_dir

  ! f2py intent declarations (must repeat for f2py parser)
  !f2py intent(in)  :: n, z, xh, xxc, xxo, t6arr, rarr, table_dir
  !f2py intent(out) :: log_kappa_arr

  ! OPAL common /e/ (must match opal_opacity.F90)
  real*4 opact, dopact, dopacr, dopactd
  common/e/ opact, dopact, dopacr, dopactd

  integer :: i
  real*4  :: xxc_i, xxo_i

  ! Change to table_dir so OPAL readco finds 'tables/codataa' etc.
  call chdir(trim(table_dir))

  do i = 1, n
    xxc_i = xxc
    xxo_i = xxo
    call opac(z, xh, xxc_i, xxo_i, t6arr(i), rarr(i))
    log_kappa_arr(i) = opact
  end do

end subroutine opal_grid
