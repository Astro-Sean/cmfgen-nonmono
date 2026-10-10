! CLI driver for the production JFNK solver (newton_hook_v1 +
! gmres_m_v1 in new_main/mod_subs/jfnk_solver_v1.f90), intended for
! subprocess-based Python testing: the harness writes a tiny command
! set on argv/stdin and reads results from stdout.  Compiled against
! the stub jfnk_cmfgen_int_mod so the real solver path is exercised
! without CMFGEN atmosphere state.
!
! Usage:
!   jfnk_cli fdjvp  <eps> <x1> <x2> <v1> <v2>    -> JV <jv1> <jv2>
!     Forward-difference J*v on the Rosenbrock residual using the
!     production JFNK_JVP_V1 step formula
!     h = eps * max(1,|x|inf) / |v|inf.
!   jfnk_cli newton <x1> <x2>                   -> INFO/NHIST/X/HIST
!     Jacobian-free Newton on Rosenbrock F = (1-x, 100(y-x^2)):
!     FD-JVP (machine-eps rule), exact Jacobi preconditioner,
!     residual history recorded via the sub() callback.
!   jfnk_cli gmres <n> <m> <itmax> <pc> <afile> <bfile>
!     gmres_m_v1 on dense A (afile: n*n values, bfile: n values).
!     pc=0 identity psolve, pc=1 Jacobi (v/diag(A)).  For each
!     k=1..itmax a fresh solve capped at k Krylov iterations gives
!     the true ||A x_k - b|| curve: "<k> <its> <res>".

module jfnk_cli_prob
   implicit none
   double precision, allocatable :: amat(:,:), pdiag(:)
   double precision, allocatable :: hist_res(:), hist_dtau(:)
   integer :: hist_n = 0
   logical :: hist_on = .false.
end module jfnk_cli_prob


subroutine rosen_f(n, x, fx)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: fx(n)
   fx(1) = 1d0 - x(1)
   fx(2) = 1d2 * (x(2) - x(1)**2)
end subroutine rosen_f


subroutine rosen_sub()
   use jfnk_cli_prob
   use jfnk_solver_v1_mod, only: jfnk_tol, jfnk_ptc_dtau
   implicit none
   if (hist_on) then
      hist_n = hist_n + 1
      if (hist_n > size(hist_res)) then
         hist_res = [hist_res, 0d0 * hist_res]
         hist_dtau = [hist_dtau, 0d0 * hist_dtau]
      end if
      hist_res(hist_n) = jfnk_tol
      hist_dtau(hist_n) = jfnk_ptc_dtau
   end if
end subroutine rosen_sub


subroutine jvp_fd_rosen(n, v, jv)
   ! Jacobian-free matvec: forward difference of the Rosenbrock
   ! residual with the production JFNK_JVP_V1 step rule
   ! (h = eps*max(1,|x|inf)/|v|inf, eps = sqrt(eps_mach)).
   use jfnk_solver_v1_mod, only: jfnk_x
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   call fd_jvp_rosen(n, jfnk_x, v, sqrt(epsilon(1d0)), jv)
end subroutine jvp_fd_rosen


subroutine fd_jvp_rosen(n, x, v, eps, jv)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n), v(n), eps
   double precision, intent(out) :: jv(n)
   double precision :: h, xp(n), fp(n), f0(n)
   if (maxval(abs(v)) < 1d-300) then
      jv = 0d0
      return
   end if
   h = eps * max(1d0, maxval(abs(x))) / maxval(abs(v))
   xp = x + h * v
   call rosen_f(n, x, f0)
   call rosen_f(n, xp, fp)
   jv = (fp - f0) / h
end subroutine fd_jvp_rosen


subroutine rosen_dfp(n, v)
   ! Jacobi preconditioner for Rosenbrock: J = [[-1,0],[-200x,100]],
   ! so M^-1 v = [-v1, v2/100].  Production semantics: psolve solves
   ! M z = v and returns z in v.
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   v(1) = -v(1)
   v(2) = v(2) / 1d2
end subroutine rosen_dfp


subroutine mv_dense(n, v, jv)
   use jfnk_cli_prob
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   jv = matmul(amat, v)
end subroutine mv_dense


subroutine pc_ident(n, v)
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
end subroutine pc_ident


subroutine pc_jacobi(n, v)
   use jfnk_cli_prob
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   v = v / pdiag
end subroutine pc_jacobi


double precision function dp_ident(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_ident = dot_product(a, b)
end function dp_ident


subroutine cli_reset(n)
   use jfnk_cli_prob
   use jfnk_solver_v1_mod
   implicit none
   integer, intent(in) :: n
   if (allocated(jfnk_x)) deallocate(jfnk_x)
   if (allocated(jfnk_fx)) deallocate(jfnk_fx)
   allocate(jfnk_x(n), jfnk_fx(n))
   jfnk_x = 0d0
   jfnk_fx = 0d0
   jfnk_fx_fresh = .false.
   jfnk_tol = 0d0
   jfnk_tol_prev = 0d0
   jfnk_nits = 0
   jfnk_gits = 0
   jfnk_del = 0d0
   jfnk_use_left_pc = .true.
   jfnk_max_gmres_its = 0
   jfnk_ptc_active = .false.
   jfnk_ptc_scale_ready = .false.
   jfnk_ptc_anchor_ready = .false.
   jfnk_ptc_dtau = 1d0
   jfnk_jvp_failed = .false.
   jfnk_pc_clip_event = .false.
   jfnk_need_precond_refresh = 0
end subroutine cli_reset


program jfnk_cli_driver
   use jfnk_cli_prob
   use jfnk_solver_v1_mod
   implicit none
   character(len=32) :: cmd, arg
   character(len=256) :: afile, bfile
   integer :: n, m, itmax, pc_mode, k, info, its, lu, ios
   double precision :: x0(2), v2(2), eps, jv2(2)
   double precision, allocatable :: a(:,:), b(:), x(:)
   double precision, allocatable :: h(:,:), v(:,:)
   double precision :: res, del, pred_res, tru_res
   external :: rosen_f, jvp_fd_rosen, rosen_dfp, rosen_sub
   external :: mv_dense, pc_ident, pc_jacobi
   double precision, external :: dp_ident

   call get_command_argument(1, cmd)
   select case (trim(cmd))

   case ('fdjvp')
      call get_command_argument(2, arg); read(arg, *) eps
      call get_command_argument(3, arg); read(arg, *) x0(1)
      call get_command_argument(4, arg); read(arg, *) x0(2)
      call get_command_argument(5, arg); read(arg, *) v2(1)
      call get_command_argument(6, arg); read(arg, *) v2(2)
      call fd_jvp_rosen(2, x0, v2, eps, jv2)
      write(*, '(A,2(1X,ES24.16))') 'JV ', jv2

   case ('newton')
      call get_command_argument(2, arg); read(arg, *) x0(1)
      call get_command_argument(3, arg); read(arg, *) x0(2)
      call cli_reset(2)
      if (allocated(hist_res)) deallocate(hist_res, hist_dtau)
      allocate(hist_res(64), hist_dtau(64))
      hist_on = .true.
      hist_n = 0
      jfnk_x = x0
      info = 0
      call newton_hook_v1(rosen_f, jvp_fd_rosen, rosen_dfp, rosen_sub, &
           dp_ident, 10, 2, 1d-3, 1d-10, 0d0, 1d-20, 1d20, 30, info, 6)
      hist_on = .false.
      write(*, '(A,I0)') 'INFO ', info
      write(*, '(A,I0)') 'NHIST ', hist_n
      write(*, '(A,2(1X,ES24.16))') 'X ', jfnk_x
      write(*, '(A)') 'HIST'
      do k = 1, hist_n
         write(*, '(I0,1X,ES24.16,1X,ES24.16)') k - 1, &
            hist_res(k), hist_dtau(k)
      end do

   case ('gmres')
      call get_command_argument(2, arg); read(arg, *) n
      call get_command_argument(3, arg); read(arg, *) m
      call get_command_argument(4, arg); read(arg, *) itmax
      call get_command_argument(5, arg); read(arg, *) pc_mode
      call get_command_argument(6, afile)
      call get_command_argument(7, bfile)
      allocate(a(n, n), b(n), x(n), h(m + 1, m), v(n, m + 1))
      open(newunit=lu, file=trim(afile), status='old', iostat=ios)
      if (ios /= 0) stop 2
      read(lu, *) a
      close(lu)
      open(newunit=lu, file=trim(bfile), status='old', iostat=ios)
      if (ios /= 0) stop 2
      read(lu, *) b
      close(lu)
      if (allocated(amat)) deallocate(amat)
      allocate(amat(n, n)); amat = a
      if (allocated(pdiag)) deallocate(pdiag)
      allocate(pdiag(n))
      do k = 1, n
         pdiag(k) = a(k, k)
      end do
      call cli_reset(n)
      del = 0d0   ! no trust region: raw GMRES least-squares step
      h = 0d0
      v = 0d0
      do k = 1, itmax
         x = 0d0
         its = k
         info = 0
         res = 1d-14   ! entry value of res is the target tolerance
         h = 0d0
         v = 0d0
         if (pc_mode == 1) then
            call gmres_m_v1(m, n, x, b, mv_dense, pc_jacobi, &
                 dp_ident, h, v, res, del, pred_res, its, info, 6)
         else
            call gmres_m_v1(m, n, x, b, mv_dense, pc_ident, &
                 dp_ident, h, v, res, del, pred_res, its, info, 6)
         end if
         ! True unpreconditioned residual of the returned iterate,
         ! plus the solver's own LS residual (preconditioned norm when
         ! pc=1) so the harness can compare each method's minimized
         ! quantity.
         tru_res = sqrt(sum((matmul(a, x) - b)**2))
         write(*, '(I0,1X,I0,2(1X,ES24.16),1X,I0)') &
            k, its, res, tru_res, info
      end do
      do k = 1, n
         write(*, '(A,1X,ES24.16)') 'X', x(k)
      end do

   case default
      write(*, '(A)') 'unknown command: ' // trim(cmd)
      stop 64
   end select
end program jfnk_cli_driver
