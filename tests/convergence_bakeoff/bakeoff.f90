!----------------------------------------------------------------------
! BAKEOFF: convergence-strategy driver for the toy NLTE problem.
!
! Usage:
!   bakeoff sanity  <param-file>          deterministic sanity checks
!   bakeoff <method> <param-file>         run one method on one instance
!   bakeoff noise   <param-file>          noise characterization probe
!
! Method strings (+-combined variants):
!   s0            stock Picard (BA step per lagged sweep)
!   andM          Anderson/Ng-equivalent on stock map, depth M (S1, C2)
!   s2f / s2c     JFNK via newton_hook_v1, frozen / coupled residual
!   +b            gross-flux-normalized metric (weighted dot)
!   +c1           secant JVP h = theta*del/|v|, latched per solve
!   +a            reduced majors-only space, traces solved in eval
!   +d            frozen Jacobi counts during the Newton solve
!
! Output: one RESULT CSV line.
!
! edited by Devin, Phase 1 convergence bake-off
!----------------------------------------------------------------------
program bakeoff
   use toy_se_mod
   use toy_ctx_mod
   use jfnk_solver_v1_mod
   implicit none
   external :: f_cb, df_cb, dfp_cb, sub_cb
   double precision, external :: dp_cb
   character(len=256) :: meth, pfile
   character(len=64)  :: w
   integer :: iarg, ios, seed, stag_cnt, outer, info, mres
   real(dp) :: rmax, rrms, rmax_true, r_prev_best, r_tr
   logical :: success
   integer :: nml_seed, nml_nd, nml_nl, nml_ntr
   real(dp) :: nml_stiff, nml_tau, nml_g0, nml_alstar, nml_rho
   real(dp) :: nml_trscale, nml_jtol, nml_jit, nml_start, nml_tfrac
   integer :: nml_jsolve, nml_solveT
   real(dp) :: nml_theta, nml_budget, nml_tol
   namelist /toy/ nml_seed, nml_nd, nml_nl, nml_ntr, nml_solveT, &
      nml_stiff, nml_tau, nml_g0, nml_alstar, nml_rho, nml_trscale, &
      nml_jsolve, nml_jtol, nml_jit, nml_start, nml_theta, nml_budget, &
      nml_tfrac, nml_tol
   !
   nml_seed=1; nml_nd=8; nml_nl=4; nml_ntr=3; nml_solveT=1
   nml_stiff=1d8; nml_tau=1d0; nml_g0=0.45d0; nml_alstar=0.5d0
   nml_rho=0d0; nml_trscale=1d-15; nml_jsolve=0; nml_jtol=1d-8
   nml_jit=0d0; nml_start=0.3d0; nml_theta=0.1d0; nml_budget=200d0
   nml_tfrac=1d-4; nml_tol=1d-8
   !
   iarg = command_argument_count()
   if(iarg < 2) then
      write(*,'(A)') 'usage: bakeoff <method|sanity|noise> <param-file>'
      call exit(2)
   end if
   call get_command_argument(1, meth)
   call get_command_argument(2, pfile)
   open(10, file=trim(pfile), status='old', iostat=ios)
   if(ios /= 0) then
      write(*,'(A)') 'cannot open param file '//trim(pfile)
      call exit(2)
   end if
   read(10, nml=toy)
   close(10)
   !
   seed       = nml_seed
   nd         = nml_nd
   nl         = nml_nl
   ntr        = nml_ntr
   solve_t    = (nml_solveT /= 0)
   stiff      = nml_stiff
   tau_scale  = nml_tau
   g0         = nml_g0
   al_star    = nml_alstar
   rho_lag    = nml_rho
   tr_scale   = nml_trscale
   jsolve_mode= nml_jsolve
   j_tol      = nml_jtol
   jit_frac   = nml_jit
   budget     = nml_budget
   tol_frac   = nml_tfrac
   tol_r      = nml_tol
   !
   call toy_alloc()
   call toy_build_rates(seed)
   call init_start(nml_start, seed)
   !
   if(trim(meth) == 'sanity') then
      call run_sanity()
      call exit(0)
   else if(trim(meth) == 'noise') then
      call run_noise()
      call exit(0)
   end if
   !
   call parse_method(meth)
   !
   call counters_reset()
   n_fallback = 0; n_outer = 0; n_and_rej = 0
   mres = merge(nred, neq, red_space)
   if(allocated(wrow)) deallocate(wrow)
   allocate(wrow(mres)); wrow = 1d0
   call refresh_L(pops)             ! initial formal solve -> j_store
   call max_r(pops, 2, rmax, rrms)  ! cheap baseline residual measure
   r_hist(0) = rmax
   r_prev_best = rmax
   stag_cnt = 0
   success = .false.; info = -1
   w = 'budget'
   !
   do outer = 1, max_outer
      n_outer = outer
      call one_outer_step(meth, info)
      call max_r(pops, 2, rmax, rrms)
      if(len_trim(dbg_env()) > 0) then
         block
            integer :: iw, kw
            real(dp) :: rw
            call worst_row(iw, kw, rw)
            write(0,'(A,I4,3ES13.5,A,2I4,A,I3,A,2I6)') &
               '  it', outer, rmax, dble(n_sweep), dble(n_ba), &
               '  worst(i,k)=', iw, kw, '  info=', info, &
               '  cheap/jvp=', n_cheap, n_jvp
         end block
      end if
      if(rmax /= rmax .or. rmax > 1.0d3) then
         w = 'diverged'; exit
      end if
      if(rmax < r_prev_best*0.95d0) then
         r_prev_best = rmax; stag_cnt = 0
      else
         stag_cnt = stag_cnt + 1
      end if
      r_hist(outer) = rmax
      if(rmax < max(tol_r, 1.0d-12)*1.0d3) then
         call true_r(r_tr)
         if(r_tr < tol_r) then
            success = .true.; w = 'converged'; exit
         end if
      end if
      if(dble(n_sweep) + dble(n_ba) > budget) then
         w = 'budget'; exit
      end if
      if(stag_cnt >= stag_win) then
         w = 'stagnation'; exit
      end if
   end do
   !
   call true_r(rmax_true)
   write(*,'(A)') 'RESULT,'//trim(meth)//','//itoa(n_sweep)//','// &
      itoa(n_ba)//','//itoa(n_cheap)//','//itoa(n_jvp)//','// &
      itoa(n_outer)//','//itoa(n_fallback)//','//itoa(n_and_rej)//','// &
      trim(w)//','//es(rmax_true)//','//es(rrms)
   call exit(0)

contains

   ! initial state: the (ed,J)-linearized major equilibrium plus a
   ! bounded jitter.  At fixed ed and J the major rows are linear in
   ! n, so one small solve per depth gives the local equilibrium;
   ! without this the random start sits many dex away and every
   ! method degenerates into a clamp-rate march before the coupling
   ! physics is even exercised.  ed/T/traces filled after the solve.
   subroutine init_start(off, isd)
      implicit none
      real(dp), intent(in) :: off
      integer, intent(in)  :: isd
      integer :: k, i, j, is, ipiv(64), info
      real(dp) :: u, Ml(64,64), rhs(64), ed_g, Jg(nd), S(nd), chi(nd)
      real(dp) :: dS(nt), r_ij, r_ji
      external dgesv
      is = isd*7919 + 17
      ! first-guess field: ground-state populations, J from that guess
      pops = 0d0
      do k = 1, nd
         pops(1,k)   = ntot
         pops(nt-1,k)= ed_ext
         pops(nt,k)  = 1.0d4
         do i = 1, ntr
            pops(nl+i,k) = tr_scale*ntot
         end do
      end do
      do k = 1, nd
         call source_S(pops(:,k), pops(nt,k), S(k), dS)
         call opac_chi(pops(:,k), chi(k))
      end do
      call solve_J(S, chi, Jg, i)
      do k = 1, nd
         ! linear major equilibrium at fixed ed_g, Jg(k)
         ed_g = ed_ext
         Ml = 0d0; rhs = 0d0
         do i = 1, nl-1
            do j = 1, nl
               if(j == i) cycle
               r_ij = cmat(i,j)*ed_g + rmat(i,j)*Jg(k)
               r_ji = cmat(j,i)*ed_g + rmat(j,i)*Jg(k) + &
                      rmat(j,i)*1.0d-3
               Ml(i,j) = r_ji
               Ml(i,i) = Ml(i,i) - r_ij
            end do
         end do
         do j = 1, nl
            Ml(nl,j) = 1d0
         end do
         rhs(nl) = ntot
         call dgesv(nl, 1, Ml, 64, ipiv, rhs, 64, info)
         if(info /= 0) then
            do i = 1, nl
               pops(i,k) = ntot/nl
            end do
         else
            do i = 1, nl
               pops(i,k) = max(rhs(i), tr_scale*ntot*1d-6)
            end do
         end if
         ! ED consistent with the solved ion stage
         pops(nt-1,k) = pops(nl,k) + ed_ext
         ! jitter the majors around equilibrium
         do i = 1, nl
            u = toy_rand(is)
            pops(i,k) = pops(i,k) * exp(off*(u - 0.5d0)*2.0d0)
         end do
         u = toy_rand(is)
         pops(nt-1,k) = pops(nt-1,k)*exp(0.1d0*off*(u-0.5d0))
         u = toy_rand(is)
         pops(nt,k)   = pops(nt,k)*exp(0.1d0*off*(u-0.5d0))
         do i = 1, ntr
            u = toy_rand(is)
            pops(nl+i,k) = tr_scale*ntot * 10.0d0**(-dble(i)) * &
                           exp(off*(u - 0.5d0)*2.0d0)
         end do
         ! renormalize majors to keep conservation exact at start
         pops(1:nl,k) = pops(1:nl,k) * ntot/sum(pops(1:nl,k))
      end do
      ! traces at their equilibrium for the first-guess field
      j_store = Jg
      call trace_solve(pops, Jg)
      n0 = pops
   end subroutine init_start

   ! parse method string: base[+variant...]
   subroutine parse_method(m)
      implicit none
      character(len=*), intent(in) :: m
      character(len=256) :: rest, tok
      integer :: p
      rest = m
      p = index(rest, '+')
      if(p > 0) then
         tok = rest(1:p-1); rest = rest(p+1:)
      else
         tok = rest; rest = ''
      end if
      if(tok == 's0' .or. tok(1:3) == 'and') then
         res_mode = 2
      else if(tok == 's2f') then
         res_mode = 1
      else if(tok == 's2c') then
         res_mode = 0
      else
         write(*,'(A)') 'unknown base method '//trim(tok)
         call exit(2)
      end if
      do while(len_trim(rest) > 0)
         p = index(rest, '+')
         if(p > 0) then
            tok = rest(1:p-1); rest = rest(p+1:)
         else
            tok = rest; rest = ''
         end if
         select case(trim(tok))
         case('b');  use_wdot = .true.
         case('c1'); latch_h = .true.; theta_c1 = nml_theta
         case('a');  red_space = .true.; tr_in_eval = .true.
         case('d');  smooth_d_run = .true.
         case default
            write(*,'(A)') 'unknown variant '//trim(tok)
            call exit(2)
         end select
      end do
      if(red_space) call build_red_map()
   end subroutine parse_method

   subroutine one_outer_step(m, info)
      implicit none
      character(len=*), intent(in) :: m
      integer, intent(inout) :: info
      character(len=256) :: base
      integer :: p
      base = m
      p = index(base, '+')
      if(p > 0) base = base(1:p-1)
      if(base(1:3) == 'and') then
         call step_anderson(base(4:4))
      else if(base == 's0') then
         call step_stock()
      else
         call step_jfnk(info)
      end if
      call refresh_L(pops)
   end subroutine one_outer_step

   ! ---- S0: one stock BA step --------------------------------------
   subroutine step_stock()
      implicit none
      real(dp) :: Fv(nt,nd), Gd(nt,nd), dn(nt,nd)
      real(dp) :: Fr(nred), dnr(nred)
      call eval_F(pops, Fv, Gd, res_mode)
      call build_BA(pops)
      if(.not. ba_ready) return
      if(red_space) then
         call rows_to_red(Fv, Fr)
         call factor_ba_red()
         if(.not. ba_red_ready) return
         dnr = Fr
         call apply_ba_red(dnr)
         call apply_step_red(dnr)
      else
         dn = Fv
         call apply_BA(dn)
         call apply_step(dn)
      end if
   end subroutine step_stock

   ! clamped log-space update (production MAX_CHNG analog)
   subroutine apply_step(dn)
      implicit none
      real(dp), intent(in) :: dn(nt,nd)
      real(dp) :: st
      integer :: k, i
      do k = 1, nd
         do i = 1, nt
            if(.not. solve_t .and. i == nt) cycle
            st = dn(i,k)/max(pops(i,k), p_floor)
            st = min(max(st, -step_cap), step_cap)
            pops(i,k) = pops(i,k)*exp(st)
            if(abs(st) >= step_cap) n_clip = n_clip + 1
         end do
      end do
   end subroutine apply_step

   subroutine apply_step_red(dnr)
      implicit none
      real(dp), intent(in) :: dnr(nred)
      real(dp) :: st
      integer :: k, m, nmj, irow, iloc
      nmj = nl + 2
      do k = 1, nd
         do m = 1, nmj
            irow = red_map((k-1)*nmj + m)
            iloc = mod(irow-1, nt) + 1
            if(.not. solve_t .and. iloc == nt) cycle
            st = dnr((k-1)*nmj+m)/max(pops(iloc,k), p_floor)
            st = min(max(st, -step_cap), step_cap)
            pops(iloc,k) = pops(iloc,k)*exp(st)
            if(abs(st) >= step_cap) n_clip = n_clip + 1
         end do
      end do
      call trace_solve(pops, j_store)
   end subroutine apply_step_red

   ! ---- S1/C2: Anderson(m) on the stock Picard map -------------------
   ! x = ln(pops); g(x) = x + clamped BA step.  Normal-system solve
   ! for the mixing coefficients, residual-increase safeguard.
   ! Ng-equivalent: stock NGACCEL is a secant fit on the same map;
   ! this Anderson implementation is the documented substitute.
   subroutine step_anderson(mchar)
      implicit none
      character(len=1), intent(in) :: mchar
      integer, save :: and_m = 2, and_n = 0, and_seed = -1
      real(dp), save, allocatable :: xh(:,:), fh(:,:)
      real(dp), allocatable :: xv(:), gv(:), xnew(:), alpha(:)
      real(dp) :: Fv(nt,nd), Gd(nt,nd), dn(nt,nd), Fr(nred), dnr(nred)
      real(dp) :: Am(8,8), rhs(8), r_old, r_new, rr, wgt
      real(dp) :: pops_sav(nt,nd)
      integer :: mm, i, j, info, ipiv(8), mres, irow, iloc, k
      integer :: nmj
      external dgesv
      mres = merge(nred, neq, red_space)
      nmj = nl + 2
      if(allocated(wrow)) then
         if(size(wrow) /= mres) then
            deallocate(wrow); allocate(wrow(mres)); wrow = 1d0
         end if
      else
         allocate(wrow(mres)); wrow = 1d0
      end if
      if(and_seed /= seed .or. .not. allocated(xh)) then
         and_seed = seed
         and_m = ichar(mchar) - ichar('0')
         and_n = 0
         if(allocated(xh)) deallocate(xh, fh)
         allocate(xh(mres,and_m), fh(mres,and_m))
      end if
      allocate(xv(mres), gv(mres), xnew(mres), alpha(and_m))
      ! Picard step on the active (possibly reduced) space
      call eval_F(pops, Fv, Gd, res_mode)
      call build_BA(pops)
      if(.not. ba_ready) return
      if(use_wdot) then
         if(red_space) then
            call rows_to_red(Gd, dnr)
            wrow = 1d0/max(abs(dnr), tiny(1d0))
         else
            wrow = 1d0/max(abs(reshape(Gd,(/mres/))), tiny(1d0))
         end if
      end if
      if(red_space) then
         call rows_to_red(Fv, Fr)
         call factor_ba_red()
         if(.not. ba_red_ready) return
         dnr = Fr
         call apply_ba_red(dnr)
         ! dnr is now the reduced n-space correction
         do i = 1, nred
            irow = red_map(i); iloc = mod(irow-1,nt)+1
            k = (irow-1)/nt + 1
            xv(i) = log(max(pops(iloc,k),p_floor))
            dnr(i) = min(max(dnr(i)/max(pops(iloc,k),p_floor), &
               -step_cap), step_cap)
            gv(i) = xv(i) + dnr(i)
         end do
      else
         dn = Fv
         call apply_BA(dn)
         do k = 1, nd
            do i = 1, nt
               xv((k-1)*nt+i) = log(max(pops(i,k),p_floor))
               dn(i,k) = min(max(dn(i,k)/max(pops(i,k),p_floor), &
                  -step_cap), step_cap)
               gv((k-1)*nt+i) = xv((k-1)*nt+i) + dn(i,k)
            end do
         end do
      end if
      ! push history (newest at column 1)
      do i = and_m, 2, -1
         xh(:,i) = xh(:,i-1); fh(:,i) = fh(:,i-1)
      end do
      xh(:,1) = xv
      fh(:,1) = gv - xv
      and_n = and_n + 1
      mm = min(and_n, and_m)
      if(mm < 2) then
         xnew = gv
      else
         Am = 0d0; rhs = 0d0
         do i = 1, mm-1
            do j = 1, mm-1
               Am(i,j) = dot_product(wrow(1:mres)*(fh(:,i)-fh(:,mm)), &
                                     wrow(1:mres)*(fh(:,j)-fh(:,mm)))
            end do
            rhs(i) = dot_product(wrow(1:mres)*(fh(:,i)-fh(:,mm)), &
                                 wrow(1:mres)*fh(:,mm))
         end do
         do i = 1, mm-1
            Am(i,i) = Am(i,i) + 1.0d-12*max(1d0,Am(i,i))
         end do
         call dgesv(mm-1, 1, Am, 8, ipiv, rhs, 8, info)
         if(info /= 0) then
            xnew = gv
         else
            alpha(mm) = 1d0
            do i = 1, mm-1
               alpha(i) = -rhs(i)
               alpha(mm) = alpha(mm) - alpha(i)
            end do
            xnew = 0d0
            do i = 1, mm
               xnew = xnew + alpha(i)*(xh(:,i) + fh(:,i))
            end do
         end if
      end if
      ! safeguard: reject if the lagged residual worsens
      pops_sav = pops
      r_old = r_hist(max(0,n_outer-1))
      do i = 1, mres
         if(red_space) then
            irow = red_map(i)
         else
            irow = i
         end if
         iloc = mod(irow-1,nt)+1; k = (irow-1)/nt + 1
         pops(iloc,k) = exp(xnew(i))
      end do
      if(red_space) call trace_solve(pops, j_store)
      call max_r(pops, 2, r_new, rr)
      ! strict non-increase: the production resacc gate rejects any
      ! step that raises the residual; a 2% tolerance absorbs jitter
      if(r_new > 1.02d0*max(r_old,1d-30)) then
         n_and_rej = n_and_rej + 1
         pops = pops_sav
         do i = 1, mres
            if(red_space) then
               irow = red_map(i)
            else
               irow = i
            end if
            iloc = mod(irow-1,nt)+1; k = (irow-1)/nt + 1
            pops(iloc,k) = exp(gv(i))
         end do
         if(red_space) call trace_solve(pops, j_store)
      end if
      deallocate(xv, gv, xnew, alpha)
   end subroutine step_anderson

   ! ---- S2: JFNK activation ------------------------------------------
   subroutine step_jfnk(info)
      implicit none
      integer, intent(inout) :: info
      real(dp) :: Fv(nt,nd), Gd(nt,nd), Fr(nred), Gr(nred)
      real(dp) :: gtol_l, tol_l, del_l, fnrm0
      integer :: mres, ni
      mres = merge(nred, neq, red_space)
      n0 = pops
      if(allocated(wrow)) then
         if(size(wrow) /= mres) then
            deallocate(wrow); allocate(wrow(mres)); wrow = 1d0
         end if
      else
         allocate(wrow(mres)); wrow = 1d0
      end if
      ! baseline residual + weights (frozen at admission, like RSCALE)
      call eval_F(pops, Fv, Gd, res_mode)
      if(red_space) then
         call rows_to_red(Fv, Fr)
         call rows_to_red(Gd, Gr)
         if(allocated(jfnk_fx)) deallocate(jfnk_fx)
         allocate(jfnk_fx(mres))
         jfnk_fx = Fr
         if(use_wdot) wrow = 1d0/max(abs(Gr), tiny(1d0))
         fnrm0 = sqrt(dp_cb(mres, jfnk_fx, jfnk_fx))
      else
         if(allocated(jfnk_fx)) deallocate(jfnk_fx)
         allocate(jfnk_fx(neq))
         jfnk_fx = reshape(Fv, (/neq/))
         if(use_wdot) wrow = 1d0/max(abs(reshape(Gd,(/neq/))), &
                                    tiny(1d0))
         fnrm0 = sqrt(dp_cb(neq, jfnk_fx, jfnk_fx))
      end if
      jfnk_fx_fresh = .true.
      call set_ptc(mres, jfnk_fx)
      jfnk_ptc_active = .true.
      jfnk_ptc_dtau = 1.0d0
      call build_BA(pops)
      if(red_space) call factor_ba_red()
      if(allocated(jfnk_x)) deallocate(jfnk_x)
      allocate(jfnk_x(mres))
      jfnk_x = 0d0
      gtol_l = 1.0d-2
      tol_l  = tol_frac*fnrm0
      del_l  = -1d0
      if(latch_h) h_latch = -1d0
      if(smooth_d_run) jcount_latch = .true.
      ni = 60
      info = 0
      n_cheap_entry = n_cheap
      ! mndl floor matches production: the driver enforces
      ! JFNK_MNDL_LOC >= 1d-8 in jfnk_solve_corr_v1.f
      call newton_hook_v1(f_cb, df_cb, dfp_cb, sub_cb, dp_cb, &
                          20, mres, gtol_l, tol_l, del_l, &
                          1.0d-8, 1.0d99, ni, info, 0)
      jcount_latch = .false.
      if(len_trim(dbg_env()) > 0) write(0,'(A,6I8)') &
         '     jfnk nits/gits/feval/jvp/pcalls/rej =', jfnk_nits, &
         jfnk_gits, jfnk_total_f_evals, jfnk_total_jvp_evals, &
         jfnk_precond_calls, jfnk_nan_rejects
      if(red_space) then
         call red_commit()
      else
         call x_to_pops(jfnk_x, pops)
      end if
      ! production driver falls back to a stock step on inner failure
      if(info /= 0) then
         n_fallback = n_fallback + 1
         call step_stock()
      end if
   end subroutine step_jfnk

   subroutine set_ptc(mres, f0)
      implicit none
      integer, intent(in) :: mres
      real(dp), intent(in) :: f0(mres)
      real(dp) :: fmax, floor_use
      if(allocated(jfnk_ptc_scale)) deallocate(jfnk_ptc_scale, jfnk_ptc_x0)
      allocate(jfnk_ptc_scale(mres), jfnk_ptc_x0(mres))
      fmax = maxval(abs(f0))
      floor_use = max(1.0d-30, 1.0d-8*fmax)
      jfnk_ptc_scale = max(abs(f0), floor_use)
      jfnk_ptc_x0 = 0d0
      jfnk_ptc_scale_ready = .true.
      jfnk_ptc_anchor_ready = .true.
   end subroutine set_ptc

   subroutine red_commit()
      implicit none
      real(dp) :: np(nt,nd)
      call red_to_pops(jfnk_x, np)
      pops = np
      call trace_solve(pops, j_store)
   end subroutine red_commit

   ! locate the worst row (diagnostic)
   subroutine worst_row(i_out, k_out, rr_out)
      implicit none
      integer, intent(out) :: i_out, k_out
      real(dp), intent(out) :: rr_out
      real(dp) :: Fv(nt,nd), Gd(nt,nd), rr
      integer :: i, k
      call eval_F(pops, Fv, Gd, 2)
      rr_out = 0d0; i_out = 0; k_out = 0
      do k = 1, nd
         do i = 1, nt
            if(Gd(i,k) > 1d-200) then
               rr = abs(Fv(i,k))/Gd(i,k)
               if(rr > rr_out) then
                  rr_out = rr; i_out = i; k_out = k
               end if
            end if
         end do
      end do
   end subroutine worst_row

   ! true residual: exact J, no lag blend (one sweep)
   subroutine true_r(r)
      implicit none
      real(dp), intent(out) :: r
      real(dp) :: rho_sav, rr
      rho_sav = rho_lag
      rho_lag = 0d0
      call max_r(pops, 0, r, rr)
      rho_lag = rho_sav
   end subroutine true_r

   ! ---- noise probe --------------------------------------------------
   ! Reproduces the production FD-window experiment on the toy: for a
   ! fixed direction v (fractional perturbation), |F(n(1+hv))-F(n)|/h
   ! is compared with |BA v| (the approximate-Jacobian action).  A
   ! window where the ratio is ~1 = trustworthy FD JVPs; divergence
   ! at small h = evaluation non-smoothness (Jacobi-count switching),
   ! divergence at large h = truncation.  Run on instances with
   ! nml_jsolve=0 vs 1 to separate deterministic vs iterative noise.
   subroutine run_noise()
      implicit none
      real(dp) :: F1(nt,nd), F2(nt,nd), F3(nt,nd), Gd(nt,nd)
      real(dp) :: vv(nt,nd), np(nt,nd), fd(nt,nd), bav(nt,nd)
      real(dp) :: hrel, nfd, nbv, num
      integer  :: ih, k, i, is
      call refresh_L(pops)
      jcount_base = -1
      call eval_F(pops, F1, Gd, 0)
      call eval_F(pops, F2, Gd, 0)
      write(*,'(A,ES12.4)') 'NOISE same-x |F1-F2|/max|F1| = ', &
         maxval(abs(F1-F2))/max(1d-300,maxval(abs(F1)))
      ! deterministic unit direction in fractional-population space
      is = seed*2654435761 + 97
      do k = 1, nd
         do i = 1, nt
            vv(i,k) = 2d0*toy_rand(is) - 1d0
         end do
      end do
      vv = vv/maxval(abs(vv))
      call build_BA(pops)          ! ba_raw holds the unfactored block
      do k = 1, nd
         bav(:,k) = matmul(ba_raw(:,:,k), vv(:,k)*pops(:,k))
      end do
      nbv = sqrt(sum(bav*bav))
      write(*,'(A)') 'NOISE h_rel      |FD|/|BAv|      jac(0) jac(h)'
      do ih = 0, 10
         hrel = 10.0d0**(-dble(ih))
         np = pops*(1d0 + hrel*vv)
         call eval_F(np, F3, Gd, 0)
         fd = (F3 - F1)/hrel
         nfd = sqrt(sum(fd*fd))
         write(*,'(A,ES10.2,ES16.6,2I7)') 'NOISE ', hrel, &
            nfd/max(nbv,1d-300), jit_hist(1), n_jac_tot
      end do
      ! The cancellation noise floor only binds where |F| << |G|:
      ! near the fixed point.  March to convergence with stock steps
      ! first (the state Newton actually polishes), then rescan.
      block
         integer :: it
         do it = 1, 300
            call step_stock(); call refresh_L(pops)
            call max_r(pops, 2, rmax_true, rrms)
            if(rmax_true < 1.0d-9) exit
            if(dble(n_sweep)+dble(n_ba) > budget) exit
         end do
      end block
      call eval_F(pops, F1, Gd, 0)
      call build_BA(pops)
      do k = 1, nd
         bav(:,k) = matmul(ba_raw(:,:,k), vv(:,k)*pops(:,k))
      end do
      nbv = sqrt(sum(bav*bav))
      write(*,'(A,ES12.4)') 'NOISE post-stock max|r| mode2 = ', &
         rmax_true
      write(*,'(A)') 'NOISE h_rel      |FD|/|BAv| (near fixed point)'
      do ih = 0, 12
         hrel = 10.0d0**(-dble(ih))
         np = pops*(1d0 + hrel*vv)
         call eval_F(np, F3, Gd, 0)
         fd = (F3 - F1)/hrel
         nfd = sqrt(sum(fd*fd))
         write(*,'(A,ES10.2,ES16.6)') 'NOISE ', hrel, &
            nfd/max(nbv,1d-300)
      end do
   end subroutine run_noise

   ! ---- sanity: same fixed point from every converging method --------
   subroutine run_sanity()
      implicit none
      real(dp) :: ra, rb, rc
      real(dp) :: nst(nt,nd), n_s0(nt,nd), n_s2f(nt,nd), n_s2c(nt,nd)
      nst = pops
      ! S0
      pops = nst; n0 = nst
      res_mode = 2; red_space=.false.; tr_in_eval=.false.
      use_wdot=.false.; latch_h=.false.; smooth_d_run=.false.
      call counters_reset(); call refresh_L(pops)
      block
         integer :: it
         do it = 1, 300
            call step_stock(); call refresh_L(pops)
            call true_r(ra)
            if(ra < tol_r) exit
            if(dble(n_sweep)+dble(n_ba) > budget) exit
         end do
      end block
      n_s0 = pops
      write(*,'(A,ES12.4,A,I5,A,I5)') 'SANITY s0   maxr=', ra, &
         '  sweeps=', n_sweep, '  ba=', n_ba
      ! S2 frozen
      pops = nst
      res_mode = 1
      call counters_reset(); call refresh_L(pops)
      block
         integer :: it, inf
         do it = 1, 60
            call step_jfnk(inf); call refresh_L(pops)
            call true_r(rb)
            if(rb < tol_r) exit
            if(dble(n_sweep)+dble(n_ba) > budget) exit
         end do
      end block
      n_s2f = pops
      write(*,'(A,ES12.4,A,I5,A,I5,A,I6)') 'SANITY s2f  maxr=', rb, &
         '  sweeps=', n_sweep, '  ba=', n_ba, '  cheap=', n_cheap
      ! S2 coupled
      pops = nst
      res_mode = 0
      call counters_reset(); call refresh_L(pops)
      block
         integer :: it, inf
         do it = 1, 60
            call step_jfnk(inf); call refresh_L(pops)
            call true_r(rc)
            if(rc < tol_r) exit
            if(dble(n_sweep)+dble(n_ba) > budget) exit
         end do
      end block
      n_s2c = pops
      write(*,'(A,ES12.4,A,I5,A,I5,A,I6)') 'SANITY s2c  maxr=', rc, &
         '  sweeps=', n_sweep, '  ba=', n_ba, '  cheap=', n_cheap
      write(*,'(A,ES12.4)') 'SANITY |s0-s2f|rel = ', &
         maxval(abs(n_s0-n_s2f)/max(1d-300,maxval(abs(n_s0))))
      write(*,'(A,ES12.4)') 'SANITY |s0-s2c|rel = ', &
         maxval(abs(n_s0-n_s2c)/max(1d-300,maxval(abs(n_s0))))
   end subroutine run_sanity

   function dbg_env()
      character(len=8) :: dbg_env
      call get_environment_variable('DBG', dbg_env)
   end function dbg_env

   character(len=12) function itoa(i)
      integer, intent(in) :: i
      write(itoa,'(I12)') i
      itoa = adjustl(itoa)
   end function itoa

   character(len=14) function es(x)
      real(dp), intent(in) :: x
      write(es,'(ES13.6)') x
      es = adjustl(es)
   end function es

end program bakeoff


!======================================================================
! newton_hook_v1 callbacks (external, mirroring the production
! interface routines).  Context comes from toy_ctx_mod.
!======================================================================

subroutine f_cb(n, x, y)
   use toy_se_mod
   use toy_ctx_mod
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_x0, &
      jfnk_ptc_anchor_ready
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   real(dp) :: n_cur(nt,nd), Fv(nt,nd), Gd(nt,nd), Fr(nred)
   integer :: i
   f_sent_last = .false.
   ! eval budget: mirror production's sentinel (Y huge -> step
   ! rejected, trust region collapses, driver falls back).  Without
   ! it one coupled activation can burn the whole budget in a single
   ! newton_hook_v1 call.
   if(res_mode == 0 .and. dble(n_sweep) + dble(n_ba) > budget) then
      y = 1.0d20
      f_sent_last = .true.
      return
   end if
   ! cheap-eval cap: per-activation like production's JFNK_EVAL_LIM
   ! (JFNK_N_RESIDUAL_EVALS resets in JFNK_INIT_V1).  A global cap
   ! poisons every later activation once crossed, and the sentinel
   ! must propagate via jfnk_jvp_failed (set in df_cb) or the finite
   ! 1e20 JVP is "valid" and the trust-region recompute cycles
   ! forever: the ared<0/pred_negligible/rho<0.75 rejects never
   ! increment n_reject and del parks at 0.5*|s| above mndl.
   if(n_cheap - n_cheap_entry > eval_cap) then
      y = 1.0d20
      f_sent_last = .true.
      return
   end if
   if(red_space) then
      call red_to_pops(real(x,dp), n_cur)
      call eval_F(n_cur, Fv, Gd, res_mode)
      call rows_to_red(Fv, Fr)
      y = Fr
   else
      call x_to_pops(real(x,dp), n_cur)
      call eval_F(n_cur, Fv, Gd, res_mode)
      y = reshape(Fv, (/n/))
   end if
   if(jfnk_ptc_active .and. jfnk_ptc_scale_ready .and. &
      jfnk_ptc_anchor_ready) then
      if(allocated(jfnk_ptc_x0) .and. allocated(jfnk_ptc_scale)) then
         if(size(jfnk_ptc_x0) == n) then
            do i = 1, n
               y(i) = y(i) + (x(i) - jfnk_ptc_x0(i))* &
                  jfnk_ptc_scale(i)/max(jfnk_ptc_dtau,1d-8)
            end do
         end if
      end if
   end if
end subroutine f_cb


subroutine df_cb(n, v, jv)
   use toy_se_mod
   use toy_ctx_mod
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx, jfnk_del, &
      jfnk_jvp_failed
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   double precision :: h, vinf, xpert(n), fpert(n)
   n_jvp = n_jvp + 1
   vinf = maxval(abs(v))
   if(vinf == 0d0) then
      jv = 0d0; return
   end if
   if(latch_h) then
      ! |h*v| = theta*del latched at the solve's first JVP -> a fixed
      ! linear operator for the whole GMRES cycle
      if(h_latch <= 0d0) h_latch = theta_c1*max(jfnk_del,1d-8)
      h = h_latch/vinf
   else
      h = sqrt(epsilon(1d0))*max(1d0,maxval(abs(jfnk_x)))/vinf
   end if
   xpert = jfnk_x + h*v
   call f_cb(n, xpert, fpert)
   ! production JFNK_JVP_V1 sets jfnk_jvp_failed when the perturbed
   ! eval comes back flagged (jfnk_resid_nbad > 0) -- the solver then
   ! aborts GMRES instead of interpreting a sentinel JVP as valid.
   if(f_sent_last) then
      jfnk_jvp_failed = .true.
      jv = 0d0
      return
   end if
   jv = (fpert - jfnk_fx)/h
end subroutine df_cb


subroutine dfp_cb(n, v)
   use toy_se_mod
   use toy_ctx_mod
   use jfnk_solver_v1_mod, only: jfnk_x
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   real(dp) :: w(nt,nd), wr(nred), npop(nt,nd)
   integer :: i, k, iloc
   if(red_space) then
      call red_to_pops(real(jfnk_x,dp), npop)
      wr = real(v,dp)
      call apply_ba_red(wr)
      do i = 1, nred
         iloc = mod(red_map(i)-1, nt) + 1
         k = (red_map(i)-1)/nt + 1
         v(i) = wr(i)/max(npop(iloc,k), p_floor)
      end do
   else
      call x_to_pops(real(jfnk_x,dp), npop)
      w = reshape(real(v,dp), (/nt,nd/))
      call apply_BA(w)
      w = w/max(npop, p_floor)
      v = reshape(w, (/n/))
   end if
end subroutine dfp_cb


subroutine sub_cb()
   implicit none
end subroutine sub_cb


double precision function dp_cb(n, a, b)
   use toy_ctx_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_cb = dot_product(a, b)
   if(use_wdot) then
      if(allocated(wrow)) then
         if(size(wrow) == n) dp_cb = dot_product(wrow*a, wrow*b)
      end if
   end if
end function dp_cb
