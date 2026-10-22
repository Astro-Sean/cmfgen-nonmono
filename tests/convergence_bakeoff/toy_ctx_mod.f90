!----------------------------------------------------------------------
! TOY_CTX_MOD: bake-off driver context.  Holds the per-run method
! configuration that the external newton_hook_v1 callbacks read:
!
!   res_mode     : eval mode for the JFNK residual (1=frozen, 0=coupled)
!   use_wdot     : candidate B -- weighted dot with w_i = 1/G_i so
!                  |F|_w is the gross-flux-normalized norm
!   red_space    : candidate A -- Krylov/Newton on majors+ED+T only;
!                  traces solved exactly inside each F eval
!   theta_c1     : candidate C1 -- secant FD step h*||v|| = theta*del
!                  (del = trust radius), latched for the whole solve
!   smooth_d_run : candidate D -- freeze Jacobi counts for the solve
!   tr_in_eval   : (module-level in toy_se_mod) run trace_solve inside
!                  eval_F after the J determination
!
! edited by Devin, Phase 1 convergence bake-off
!----------------------------------------------------------------------
module toy_ctx_mod
   use toy_se_mod
   implicit none

   integer  :: res_mode = 2
   logical  :: use_wdot = .false.
   ! set by f_cb when an eval returns the budget sentinel; df_cb turns
   ! it into jfnk_jvp_failed so GMRES aborts the same way production
   ! does when JFNK_RESIDUAL_V1 flags jfnk_resid_nbad
   logical  :: f_sent_last = .false.
   integer  :: n_cheap_entry = 0  ! n_cheap at activation entry
   integer  :: eval_cap = 2000    ! per-activation cheap-eval cap
   logical  :: red_space = .false.
   real(dp) :: theta_c1 = 0.0d0
   real(dp) :: h_latch = 0.0d0
   logical  :: latch_h = .false.
   logical  :: smooth_d_run = .false.

   real(dp) :: step_cap = 0.5d0    ! stock update clamp (log units)
   real(dp) :: tol_r = 1.0d-8      ! success: max|r| < max(tol, floor)
   real(dp) :: tol_frac = 1.0d-4   ! inner solver tol = frac*|F0|_w
   real(dp) :: budget = 2.0d2      ! sweep-equivalents (at c_ba=1)
   integer  :: max_outer = 400
   integer  :: stag_win = 25       ! stagnation window (outer iters)

   ! row weights for the dot product (neq or nred entries)
   real(dp), allocatable :: wrow(:)

   ! reduced-space map: full-row index for each reduced unknown
   integer  :: nred = 0
   integer,  allocatable :: red_map(:)
   real(dp), allocatable :: ba_red(:,:,:)
   integer,  allocatable :: ip_red(:,:)
   logical  :: ba_red_ready = .false.

   ! per-run method bookkeeping
   integer  :: n_fallback = 0      ! inner-solve failures -> stock step
   integer  :: n_outer = 0
   integer  :: n_and_rej = 0       ! Anderson safeguard rejections
   real(dp) :: r_hist(0:10000)

   contains

   ! reduced-map rows per depth: 1..nl (majors), nt-1 (ED), nt (T)
   subroutine build_red_map()
      implicit none
      integer :: k, m, nmj, idx
      nmj = nl + 2
      nred = nmj*nd
      if(allocated(red_map)) deallocate(red_map, ba_red, ip_red)
      allocate(red_map(nred), ba_red(nmj,nmj,nd), ip_red(nmj,nd))
      idx = 0
      do k = 1, nd
         do m = 1, nl
            idx = idx + 1
            red_map(idx) = (k-1)*nt + m
         end do
         idx = idx + 1
         red_map(idx) = (k-1)*nt + (nt-1)
         idx = idx + 1
         red_map(idx) = (k-1)*nt + nt
      end do
   end subroutine build_red_map

   ! factor the major subblocks of the assembled full BA.  One-way
   ! coupling (dF_major/dn_trace = 0) makes the subblock the exact
   ! reduced Jacobian model -- no Schur complement needed.
   subroutine factor_ba_red()
      implicit none
      integer :: k, i, j, info, nmj
      external dgetrf
      nmj = nl + 2
      do k = 1, nd
         do i = 1, nmj
            do j = 1, nmj
               ! red_map holds global row indices; extract from the
               ! unfactored ba_raw (ba_blk holds LU factors by now)
               ba_red(j,i,k) = ba_raw( &
                  mod(red_map((k-1)*nmj+j)-1,nt)+1, &
                  mod(red_map((k-1)*nmj+i)-1,nt)+1, k)
            end do
         end do
         call dgetrf(nmj, nmj, ba_red(:,:,k), nmj, ip_red(:,k), info)
         if(info /= 0) then
            ba_red_ready = .false.
            return
         end if
      end do
      ba_red_ready = .true.
   end subroutine factor_ba_red

   subroutine apply_ba_red(dv)
      implicit none
      real(dp), intent(inout) :: dv(nred)
      integer :: k, info, nmj
      external dgetrs
      nmj = nl + 2
      do k = 1, nd
         call dgetrs('N', nmj, 1, ba_red(:,:,k), nmj, ip_red(:,k), &
                     dv((k-1)*nmj+1:k*nmj), nmj, info)
      end do
      dv = -dv
   end subroutine apply_ba_red

   ! reduced x -> full pops (majors from exp map; traces by the caller)
   subroutine red_to_pops(xr, n_cur)
      implicit none
      real(dp), intent(in)  :: xr(nred)
      real(dp), intent(out) :: n_cur(nt,nd)
      integer :: k, m, nmj, idx
      real(dp) :: xt
      nmj = nl + 2
      n_cur = 0d0
      do k = 1, nd
         do m = 1, nmj
            idx = red_map((k-1)*nmj + m)
            xt = min(max(xr((k-1)*nmj+m), -x_clip), x_clip)
            ! idx is a global row index; local row = mod
            n_cur(mod(idx-1,nt)+1, k) = n0(mod(idx-1,nt)+1, k)*exp(xt)
         end do
      end do
      ! trace block: start from baseline; caller runs trace_solve
      do k = 1, nd
         do m = 1, ntr
            n_cur(nl+m,k) = n0(nl+m,k)
         end do
      end do
   end subroutine red_to_pops

   ! extract the major rows of a full residual into reduced space
   subroutine rows_to_red(F, Fr)
      implicit none
      real(dp), intent(in)  :: F(nt,nd)
      real(dp), intent(out) :: Fr(nred)
      integer :: i, nmj
      nmj = nl + 2
      do i = 1, nred
         Fr(i) = F(mod(red_map(i)-1,nt)+1, (red_map(i)-1)/nt + 1)
      end do
   end subroutine rows_to_red

end module toy_ctx_mod
