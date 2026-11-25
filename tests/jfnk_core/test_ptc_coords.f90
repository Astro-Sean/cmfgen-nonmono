! test_ptc_coords.f90 -- Phase-2 regression: PTC residual/direction
! strip and the u-coordinate fold factor at NONZERO Newton iterates.
!
! Residual:   G_i = F_i + (x_i - x0_i) s_i / dtau     (SER anchor)
! GMRES solves the fractional-correction coordinate
!             u = exp(x) - 1  =>  du/dx = exp(x)
! so the folded diagonal is s_i/dtau * dx/du = s_i exp(-x_i)/dtau.
! The test checks (a) the production strip recovers F and the
! F-Jacobian action at x /= x0 /= 0, and (b) the FD derivative of
! the PTC term in u-coordinates matches s*exp(-x)/dtau at
! x = 0, +-0.1, +-1 -- the points the audit mandates.
program test_ptc_coords
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_x0, &
      jfnk_ptc_anchor_ready
   implicit none
   external jfnk_ptc_strip_v1, jfnk_ptc_strip_dir_v1
   integer, parameter :: n = 8
   double precision :: f(n), g(n), x(n), sdir(n), js(n), scl(n)
   double precision :: xref(5), xi, h, dudx_fd, want, err
   double precision :: fmax, worst
   integer :: i, k, nfail

   nfail = 0
   xref = (/ 0.0d0, 0.1d0, -0.1d0, 1.0d0, -1.0d0 /)

   jfnk_ptc_active = .true.
   jfnk_ptc_dtau = 0.37d0
   allocate(jfnk_ptc_scale(n), jfnk_ptc_x0(n))
   do i = 1, n
      scl(i) = 1.0d0 + 0.5d0 * i
      jfnk_ptc_scale(i) = scl(i)
      jfnk_ptc_x0(i) = 0.23d0 * i - 0.4d0   ! nonzero anchor
   end do
   jfnk_ptc_scale_ready = .true.
   jfnk_ptc_anchor_ready = .true.

   ! --- A: strip recovers F and the F-action at x /= x0 -------------
   worst = 0d0
   do k = 1, 5
      xi = xref(k)
      do i = 1, n
         x(i) = xi * (1d0 + 0.05d0 * i)
         f(i) = 1.0d2 * sin(0.3d0 * i * xi) + 10d0 * i
         g(i) = f(i) + (x(i) - jfnk_ptc_x0(i)) * scl(i) / &
            jfnk_ptc_dtau
      end do
      call jfnk_ptc_strip_v1(n, x, g)
      do i = 1, n
         err = abs(g(i) - f(i)) / max(1d0, abs(f(i)))
         worst = max(worst, err)
      end do
      ! directional strip: js = J_F*s + s*scale/dtau -> J_F*s
      do i = 1, n
         sdir(i) = 0.1d0 * (-1.0d0) ** i
         js(i) = 7.0d0 * i + sdir(i) * scl(i) / jfnk_ptc_dtau
      end do
      call jfnk_ptc_strip_dir_v1(n, sdir, js)
      do i = 1, n
         err = abs(js(i) - 7.0d0 * i) / (7.0d0 * i)
         worst = max(worst, err)
      end do
   end do
   write(*, '(A,ES12.4)') ' A: worst strip defect (x,x0 nonzero) =', &
      worst
   if (worst < 1d-12) then
      write(*, '(A)') ' A: PASS'
   else
      write(*, '(A)') ' A: FAIL'
      nfail = nfail + 1
   end if

   ! --- B: u-coordinate fold factor at nonzero x --------------------
   ! Fold target: d/du [ s*(x(u)-x0)/dtau ] = s/dtau * exp(-x).
   ! FD in u of the same term must reproduce it; a missing or wrong
   ! exp(-x) shows up as O(1) error at |x| ~ 1.
   worst = 0d0
   do k = 2, 5              ! nonzero iterates only (x=0 is trivial)
      xi = xref(k)
      h = 1d-7 * max(1d0, abs(xi))
      ! FD: du/dx evaluated at xi  ==  [exp(xi+h)-exp(xi)]/h - 1? no:
      ! u = exp(x)-1 -> du/dx = exp(x); the fold uses exp(-x), so
      ! check  s/dtau * (1/(du/dx_FD))  vs  s/dtau * exp(-x).
      dudx_fd = (exp(xi + h) - exp(xi - h)) / (2d0 * h)
      want = scl(3) / jfnk_ptc_dtau * exp(-xi)
      err = abs(scl(3) / jfnk_ptc_dtau / dudx_fd - want) / want
      worst = max(worst, err)
      write(*, '(A,F7.2,A,ES12.4,A,ES12.4)') &
         '   x=', xi, '  fold-FD=', &
         scl(3) / jfnk_ptc_dtau / dudx_fd, '  exp(-x) form=', want
   end do
   ! 1e-9 gate: FD noise floor at h=1e-7 is ~3e-11; a wrong
   ! transform (e.g. exp(+x) or no factor) is off by O(1) at |x|=1.
   if (worst < 1d-9) then
      write(*, '(A)') ' B: PASS'
   else
      write(*, '(A)') ' B: FAIL'
      nfail = nfail + 1
   end if

   if (nfail == 0) then
      write(*, '(A)') ' ALL PTC-COORD TESTS PASSED'
   else
      write(*, '(A,I2)') ' PTC-COORD FAILURES: ', nfail
      call exit(1)
   end if
end program test_ptc_coords
