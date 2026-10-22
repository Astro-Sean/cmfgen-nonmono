!----------------------------------------------------------------------
! TOY_SE_MOD: toy NLTE statistical-equilibrium problem for the
! convergence bake-off.  Reproduces the structural features of the
! production population solve that decide which convergence strategy
! wins, without needing the full CMFGEN model state:
!
!   1. Depth-local nonlinear SE blocks: NL major levels per depth with
!      collisional rates (bilinear in n and ED), radiative rates
!      (linear in n at fixed J), a number-conservation row, an ED row
!      and a T row (energy balance).  Per-depth block NT = NL+NTR+2.
!      With solve_t=F the T row is masked (dummy identity row), the
!      production FIX_T analog.
!   2. Non-local radiation: J(k) solves (I - g L) J = S(n) on the
!      depth grid (tridiagonal).  The BA preconditioner sees only the
!      diagonal Lambda* response (al_star fraction of the local
!      dS/dn), so BA is a genuinely imperfect Jacobian model.
!   3. Gross-flux stiffness: rate rows scaled by STIFF, so the
!      inflow/outflow terms are STIFF times larger than the net
!      residual that conservation allows.
!   4. Trace levels: NTR rows per depth, populations ~TR_SCALE of
!      majors, one-way coupled (tr rates may depend on majors and J;
!      majors never see traces -> block lower-triangular, linear in
!      n_tr at fixed majors+J).
!   5. Lagged store J_STORE: the production I_NEXT analog.  Frozen
!      inner evals see J_STORE plus the linearized BA_rad increment;
!      coupled evals see J_eff = (1-rho)*J_exact(n) + rho*J_STORE.
!      Refreshed once per outer iteration by a real transfer solve.
!   6. Eval noise/non-smoothness: J may be solved by tol-terminated
!      Jacobi iterations (iteration count varies with input -> a
!      deterministic non-smooth residual), an optional jitter emulates
!      non-reproducible noise, and populations carry a floor clip.
!
! Cost accounting (sweep equivalents):
!   n_sweep    : evals that solve the transfer operator (coupled
!                evals, coupled JVPs, stock outer refreshes)
!   n_ba       : BA assembly + factorization count (cost c_ba each)
!   n_cheap    : frozen/lagged evals (no transfer solve)
!   sweep_eq   : n_sweep + c_ba*n_ba, reported for c_ba in {0.1,0.3,1}
!
! edited by Devin, Phase 1 convergence bake-off
!----------------------------------------------------------------------
module toy_se_mod
   implicit none
   save

   integer, parameter :: dp = kind(1.0d0)

   ! ---- problem geometry -------------------------------------------
   integer :: nd = 8            ! depths
   integer :: nl = 4            ! major levels per depth
   integer :: ntr = 3           ! trace levels per depth
   logical :: solve_t = .true.  ! include T row as a real unknown
   integer :: nt = 0            ! nl + ntr + 2 (ED at nt-1, T at nt)
   integer :: neq = 0           ! nt*nd

   ! ---- physics knobs ----------------------------------------------
   real(dp) :: stiff = 1.0d8    ! gross-flux / net scale ratio
   real(dp) :: tau_scale = 1.0d0 ! opacity scale (thin->thick)
   real(dp) :: g0 = 0.45d0      ! J coupling strength (<0.5 for stability)
   real(dp) :: al_star = 0.5d0  ! Lambda* fidelity in BA (0=lambda,1=full local)
   real(dp) :: rho_lag = 0.0d0  ! lagged-store fraction in coupled evals
   real(dp) :: tr_scale = 1.0d-15 ! trace population scale vs majors
   real(dp) :: j_inc = 1.0d0    ! incident J at the inner boundary
   real(dp) :: tr_stiff = 1.0d2 ! extra rate stiffness inside trace block

   ! ---- noise knobs -------------------------------------------------
   ! jsolve_mode: 0=exact tridiag, 1=tol-terminated Jacobi
   integer  :: jsolve_mode = 0
   real(dp) :: j_tol = 1.0d-8   ! Jacobi tolerance (iter mode)
   integer  :: j_itmax = 200    ! Jacobi iteration cap
   real(dp) :: jit_frac = 0.0d0 ! non-reproducible jitter on J
   integer  :: jit_seed = 12345
   ! when smooth_d=T, freeze the Jacobi count to the baseline count
   ! during JVPs (candidate D); baseline evals (latch off) record it.
   logical  :: smooth_d = .false.
   integer  :: jcount_base = -1
   logical  :: jcount_latch = .false.

   ! ---- state --------------------------------------------------------
   real(dp), allocatable :: pops(:,:)   ! (nt,nd) current populations
   real(dp), allocatable :: n0(:,:)     ! baseline populations for log-x
   real(dp), allocatable :: j_store(:)  ! (nd) lagged field
   real(dp), allocatable :: j_now(:)    ! (nd) last computed J
   real(dp), allocatable :: x_base(:)   ! (neq) baseline unknowns
   integer,  allocatable :: jit_hist(:) ! last Jacobi count per eval

   ! rate tables (deterministic, seed-built)
   real(dp), allocatable :: cmat(:,:)   ! (nl,nl) collisional strength
   real(dp), allocatable :: rmat(:,:)   ! (nl,nl) radiative strength
   real(dp), allocatable :: atr(:,:)    ! (ntr,ntr) trace self-rates
   real(dp), allocatable :: bsrc(:,:)   ! (ntr,nl) trace source strengths
   real(dp) :: ntot, ed_ext             ! gas number, external charge
   real(dp) :: s0, chi0                 ! source/opacity scales
   real(dp) :: heat_c, cool_c           ! T-row strengths

   ! ---- counters -----------------------------------------------------
   integer :: n_sweep = 0, n_ba = 0, n_cheap = 0, n_jvp = 0
   integer :: n_clip = 0, n_jac_tot = 0, n_jac_evals = 0

   ! ---- BA workspace --------------------------------------------------
   real(dp), allocatable :: ba_blk(:,:,:) ! (nt,nt,nd) assembled BA
   integer,  allocatable :: ba_ipiv(:,:)  ! (nt,nd) LU pivots
   logical  :: ba_ready = .false.

   ! candidate A: solve traces exactly inside each F eval, using the
   ! same J the eval sees (set by the driver context)
   logical  :: tr_in_eval = .false.

   ! floor for populations and log-unknown clip (production bounds)
   real(dp), parameter :: p_floor = 1.0d-300
   real(dp), parameter :: x_clip  = 40.0d0

contains

   !------------------------------------------------------------------
   subroutine toy_alloc()
      implicit none
      nt = nl + ntr + 2
      neq = nt*nd
      if(allocated(pops)) deallocate(pops, n0, j_store, j_now, x_base, &
         jit_hist, cmat, rmat, atr, bsrc, ba_blk, ba_ipiv)
      allocate(pops(nt,nd), n0(nt,nd), j_store(nd), j_now(nd), &
         x_base(neq), jit_hist(nd), cmat(nl,nl), rmat(nl,nl), &
         atr(ntr,ntr), bsrc(ntr,nl), ba_blk(nt,nt,nd), ba_ipiv(nt,nd))
   end subroutine toy_alloc

   ! simple LCG; deterministic across compilers
   real(dp) function toy_rand(iseed)
      implicit none
      integer, intent(inout) :: iseed
      iseed = mod(iseed*1103515245 + 12345, 2147483647)
      toy_rand = dble(iseed)/2147483647.0d0
   end function toy_rand

   ! log-spaced jittered rate tables -> stiff, badly-scaled SE matrix
   subroutine toy_build_rates(iseed)
      implicit none
      integer, intent(inout) :: iseed
      integer :: i, j
      real(dp) :: u
      do i = 1, nl
         do j = 1, nl
            if(i == j) then
               cmat(i,j) = 0d0; rmat(i,j) = 0d0
            else
               u = toy_rand(iseed)
               cmat(i,j) = stiff * 10.0d0**(-4.0d0*u - dble(abs(i-j)))
               u = toy_rand(iseed)
               rmat(i,j) = stiff * 10.0d0**(-3.0d0*u - 2.0d0*dble(abs(i-j)))
            end if
         end do
      end do
      do i = 1, ntr
         do j = 1, ntr
            if(i == j) then
               atr(i,j) = 0d0
            else
               u = toy_rand(iseed)
               atr(i,j) = tr_stiff*stiff*10.0d0**(-3.0d0*u - dble(abs(i-j)))
            end if
         end do
      end do
      do i = 1, ntr
         do j = 1, nl
            u = toy_rand(iseed)
            bsrc(i,j) = tr_scale*tr_stiff*stiff * &
                        10.0d0**(3.0d0*u - 0.5d0*dble(j) - 1.0d0)
         end do
      end do
      ntot   = 1.0d12
      ed_ext = 1.0d10
      s0     = 1.0d0
      chi0   = 1.0d-12
      heat_c = 1.0d-10
      cool_c = 1.0d-6
   end subroutine toy_build_rates

   !------------------------------------------------------------------
   ! J solve: (I - g(k) L) J = S.  g(k) = g0/(1 + tau_scale*chi(k)).
   ! Mode 0: exact tridiagonal Thomas.  Mode 1: tol-terminated Jacobi
   ! (iteration count depends on the input -> non-smooth residual).
   ! With jcount_latch=T the count is frozen at jcount_base (candidate
   ! D: makes F smooth during a Newton step).
   !------------------------------------------------------------------
   subroutine solve_J(S, chi, J, iters_used)
      implicit none
      real(dp), intent(in)  :: S(nd), chi(nd)
      real(dp), intent(out) :: J(nd)
      integer,  intent(out) :: iters_used
      real(dp) :: g(nd), al(nd), be(nd), ga(nd), rhs(nd)
      real(dp) :: Jold(nd), err, den
      integer  :: k, it, n_it_want
      do k = 1, nd
         g(k) = g0 / (1.0d0 + tau_scale*chi(k))
      end do
      if(jsolve_mode == 0) then
         al = 0d0; be = 0d0; ga = 0d0; rhs = S
         do k = 2, nd-1
            al(k) = -g(k); be(k) = 1d0 + 2d0*g(k); ga(k) = -g(k)
         end do
         be(1) = 1d0 + g(1); ga(1) = -g(1)
         al(nd) = -g(nd); be(nd) = 1d0 + g(nd)
         rhs(nd) = S(nd) + g(nd)*j_inc
         do k = 2, nd
            be(k) = be(k) - al(k)*ga(k-1)/be(k-1)
            rhs(k) = rhs(k) - al(k)*rhs(k-1)/be(k-1)
         end do
         J(nd) = rhs(nd)/be(nd)
         do k = nd-1, 1, -1
            J(k) = (rhs(k) - ga(k)*J(k+1))/be(k)
         end do
         iters_used = 1
      else
         J = S
         iters_used = 0
         n_it_want = j_itmax
         if(jcount_latch .and. jcount_base > 0) n_it_want = jcount_base
         do it = 1, n_it_want
            Jold = J
            J(1) = (S(1) + g(1)*Jold(2)) / (1d0 + g(1))
            do k = 2, nd-1
               J(k) = (S(k) + g(k)*(Jold(k-1)+Jold(k+1))) / &
                      (1d0 + 2d0*g(k))
            end do
            J(nd) = (S(nd) + g(nd)*(Jold(nd-1)+j_inc)) / &
                    (1d0 + g(nd))
            err = maxval(abs(J - Jold)/max(1d0,maxval(abs(J))))
            iters_used = it
            if(.not. jcount_latch .and. err < j_tol) exit
         end do
         if(.not. jcount_latch) jcount_base = iters_used
      end if
      if(jit_frac > 0d0) then
         ! seeded per call so repeated evals of the SAME n differ:
         ! the OMP-reduction / non-reproducible-noise analog
         do k = 1, nd
            den = toy_rand(jit_seed)
            J(k) = J(k) * (1d0 + jit_frac*(den - 0.5d0))
         end do
      end if
   end subroutine solve_J

   !------------------------------------------------------------------
   ! Source function (nonlinear: two-level ratio, weak T dependence)
   ! and opacity (linear in the ground population).
   !------------------------------------------------------------------
   subroutine source_S(n_loc, T_loc, S, dS_dn)
      implicit none
      real(dp), intent(in)  :: n_loc(nt), T_loc
      real(dp), intent(out) :: S
      real(dp), intent(out) :: dS_dn(nt)
      real(dp) :: ratio, nlo
      nlo = max(n_loc(1), p_floor)
      ratio = n_loc(2)/nlo
      S = s0 * ratio * (1d0 + 0.1d0*T_loc/1.0d4)
      dS_dn = 0d0
      dS_dn(1) = -s0 * n_loc(2)/nlo**2 * (1d0 + 0.1d0*T_loc/1.0d4)
      dS_dn(2) =  s0/nlo * (1d0 + 0.1d0*T_loc/1.0d4)
      dS_dn(nt) = s0*ratio*0.1d0/1.0d4
   end subroutine source_S

   subroutine opac_chi(n_loc, chi)
      implicit none
      real(dp), intent(in)  :: n_loc(nt)
      real(dp), intent(out) :: chi
      chi = chi0 * max(n_loc(1), p_floor)
   end subroutine opac_chi

   !------------------------------------------------------------------
   ! Residual.  mode selects the radiation treatment:
   !   0 = TRUE coupled: J = (1-rho)*J_exact(n) + rho*j_store
   !       (one transfer solve per eval -> counts a sweep)
   !   1 = frozen: coll(n) + rad(n0; j_store) + BA_rad*(n - n0)
   !   2 = stock: rad(n; j_store) with the stored field (cheap)
   ! G(i,k) = gross flux: sum of |inflow| + |outflow| terms per row.
   !------------------------------------------------------------------
   subroutine eval_F(n_cur, F, G, mode)
      implicit none
      real(dp), intent(in)  :: n_cur(nt,nd)
      real(dp), intent(out) :: F(nt,nd), G(nt,nd)
      integer,  intent(in)  :: mode
      real(dp) :: Jvv(nd), S(nd), chi(nd), dS(nt), n_eval(nt,nd)
      real(dp) :: ed, Tl, nloc(nt), infl, outfl, r_ij, r_ji
      real(dp) :: src, chi_k
      integer  :: k, i, j, it_used
      n_eval = n_cur
      if(mode == 0) then
         do k = 1, nd
            call source_S(n_eval(:,k), n_eval(nt,k), S(k), dS)
            call opac_chi(n_eval(:,k), chi(k))
         end do
         call solve_J(S, chi, Jvv, it_used)
         n_sweep = n_sweep + 1
         n_jac_tot = n_jac_tot + it_used
         n_jac_evals = n_jac_evals + 1
         Jvv = (1d0 - rho_lag)*Jvv + rho_lag*j_store
         j_now = Jvv
      else
         Jvv = j_store
         n_cheap = n_cheap + 1
      end if
      ! candidate A hook: traces eliminated exactly at the eval's J
      if(tr_in_eval) call trace_solve(n_eval, Jvv)
      F = 0d0; G = 0d0
      do k = 1, nd
         nloc = n_eval(:,k)
         ed = nloc(nt-1)
         Tl = nloc(nt)
         if(.not. solve_t) Tl = 1.0d4
         ! ---- major level rows i=1..nl-1: rate balance ----------
         do i = 1, nl-1
            if(mode == 1) then
               infl = 0d0; outfl = 0d0
               do j = 1, nl
                  if(j == i) cycle
                  infl  = infl  + nloc(j)*cmat(j,i)*ed
                  outfl = outfl + nloc(i)*cmat(i,j)*ed
               end do
               call frozen_rad_row(k, i, n_eval, infl, outfl)
            else
               infl = 0d0; outfl = 0d0
               do j = 1, nl
                  if(j == i) cycle
                  r_ij = cmat(i,j)*ed + rmat(i,j)*Jvv(k)
                  r_ji = cmat(j,i)*ed + rmat(j,i)*Jvv(k) + &
                         rmat(j,i)*1.0d-3
                  infl  = infl  + nloc(j)*r_ji
                  outfl = outfl + nloc(i)*r_ij
               end do
            end if
            F(i,k) = infl - outfl
            G(i,k) = infl + outfl
         end do
         ! ---- row nl: number conservation ------------------------
         F(nl,k) = stiff*(sum(nloc(1:nl)) - ntot)
         G(nl,k) = stiff*ntot
         ! ---- trace rows: one-way, linear in n_tr ---------------
         do i = 1, ntr
            infl = 0d0; outfl = 0d0; src = 0d0
            do j = 1, ntr
               if(j == i) cycle
               infl  = infl  + nloc(nl+j)*atr(j,i)
               outfl = outfl + nloc(nl+i)*atr(i,j)
            end do
            do j = 1, nl
               src = src + nloc(j)*bsrc(i,j)*(Jvv(k) + 0.1d0)
            end do
            outfl = outfl + nloc(nl+i)*tr_stiff*stiff
            F(nl+i,k) = src + infl - outfl
            G(nl+i,k) = src + infl + outfl
         end do
         ! ---- ED row ---------------------------------------------
         F(nt-1,k) = (nloc(nl) + ed_ext - ed)*stiff
         G(nt-1,k) = (nloc(nl) + ed_ext + ed)*stiff
         ! ---- T row ----------------------------------------------
         call opac_chi(nloc, chi_k)
         if(solve_t) then
            F(nt,k) = heat_c*chi_k*(Jvv(k) - 1.0d-16*Tl**4) - &
                      cool_c*ed*(Tl - 3.0d3)
            G(nt,k) = heat_c*chi_k*(Jvv(k) + 1.0d-16*Tl**4) + &
                      cool_c*ed*(Tl + 3.0d3)
         else
            ! masked T row (FIX_T analog): F=0, BA row = identity
            F(nt,k) = 0d0
            G(nt,k) = 1.0d-300
         end if
      end do
   end subroutine eval_F

   !------------------------------------------------------------------
   ! Frozen-mode linearized radiative increment for level row i at
   ! depth k:  rad(n0; J_store) + BA_rad*(n - n0), with BA_rad built
   ! from the same Lambda* as the preconditioner.  This is exactly
   ! the production frozen residual (STEQ_rad(n0) + CMAT_RAD*x).
   !------------------------------------------------------------------
   subroutine frozen_rad_row(k, i, n_cur, infl, outfl)
      implicit none
      integer,  intent(in) :: k, i
      real(dp), intent(in) :: n_cur(nt,nd)
      real(dp), intent(inout) :: infl, outfl
      real(dp) :: dS(nt), dn(nt), dJ_dn(nt), S_b
      integer  :: j, m
      ! baseline radiative contribution (exact at n0 with stored J)
      do j = 1, nl
         if(j == i) cycle
         infl  = infl  + n0(j,k)*(rmat(j,i)*j_store(k) + rmat(j,i)*1.0d-3)
         outfl = outfl + n0(i,k)*rmat(i,j)*j_store(k)
      end do
      ! linearized increment; dJ/dn_m = al_star*dS/dn_m (local only)
      dn = n_cur(:,k) - n0(:,k)
      call source_S(n0(:,k), n0(nt,k), S_b, dS)
      dJ_dn = al_star * dS
      do m = 1, nt
         if(dn(m) == 0d0) cycle
         do j = 1, nl
            if(j == i) cycle
            if(m == j) infl = infl + dn(m)*(rmat(j,i)*j_store(k) + &
                                           rmat(j,i)*1.0d-3)
            if(m == i) outfl = outfl + dn(m)*rmat(i,j)*j_store(k)
            infl  = infl  + n0(j,k)*rmat(j,i)*dJ_dn(m)*dn(m)
            outfl = outfl + n0(i,k)*rmat(i,j)*dJ_dn(m)*dn(m)
         end do
      end do
   end subroutine frozen_rad_row

   !------------------------------------------------------------------
   ! BA preconditioner: per-depth dense block, LU-factored.  Radiation
   ! response restricted to the Lambda* local model -> genuinely
   ! imperfect vs the true non-local Jacobian.  The masked-T row is
   ! the identity so dn_T = 0.
   !------------------------------------------------------------------
   subroutine build_BA(n_cur)
      implicit none
      real(dp), intent(in) :: n_cur(nt,nd)
      real(dp) :: dS(nt), S_b, nloc(nt), ed, Tl
      real(dp) :: r_ij, r_ji, dJ(nt)
      integer  :: k, i, j, m, info, ii
      external dgetrf
      ba_blk = 0d0
      do k = 1, nd
         nloc = n_cur(:,k)
         ed = nloc(nt-1)
         Tl = nloc(nt)
         if(.not. solve_t) Tl = 1.0d4
         call source_S(nloc, Tl, S_b, dS)
         if(.not. solve_t) dS(nt) = 0d0
         dJ = al_star*dS
         ! ---- major rows ----------------------------------------
         do i = 1, nl-1
            do j = 1, nl
               if(j == i) cycle
               r_ij = cmat(i,j)*ed + rmat(i,j)*j_store(k)
               r_ji = cmat(j,i)*ed + rmat(j,i)*j_store(k) + &
                      rmat(j,i)*1.0d-3
               ba_blk(i,j,k) = ba_blk(i,j,k) + r_ji
               ba_blk(i,i,k) = ba_blk(i,i,k) - r_ij
               ba_blk(i,nt-1,k) = ba_blk(i,nt-1,k) + &
                  nloc(j)*cmat(j,i) - nloc(i)*cmat(i,j)
               do m = 1, nt
                  ba_blk(i,m,k) = ba_blk(i,m,k) + &
                     (nloc(j)*rmat(j,i) - nloc(i)*rmat(i,j))*dJ(m)
               end do
            end do
         end do
         ! ---- conservation row -----------------------------------
         do j = 1, nl
            ba_blk(nl,j,k) = stiff
         end do
         ! ---- trace rows -----------------------------------------
         do i = 1, ntr
            ii = nl + i
            do j = 1, ntr
               if(j == i) cycle
               ba_blk(ii,nl+j,k) = atr(j,i)
               ba_blk(ii,ii,k)   = ba_blk(ii,ii,k) - atr(i,j)
            end do
            ba_blk(ii,ii,k) = ba_blk(ii,ii,k) - tr_stiff*stiff
            do j = 1, nl
               ba_blk(ii,j,k) = bsrc(i,j)*(j_store(k) + 0.1d0)
            end do
            do m = 1, nt
               do j = 1, nl
                  ba_blk(ii,m,k) = ba_blk(ii,m,k) + &
                     nloc(j)*bsrc(i,j)*dJ(m)
               end do
            end do

         end do
         ! ---- ED row ----------------------------------------------
         ba_blk(nt-1,nl,k)   = stiff
         ba_blk(nt-1,nt-1,k) = -stiff
         ! ---- T row ------------------------------------------------
         if(solve_t) then
            ba_blk(nt,1,k)  = heat_c*chi0*(j_store(k) - 1.0d-16*Tl**4)
            do m = 1, nt
               ba_blk(nt,m,k) = ba_blk(nt,m,k) + &
                  heat_c*chi0*nloc(1)*dJ(m)
            end do
            ba_blk(nt,nt,k) = ba_blk(nt,nt,k) &
               - heat_c*chi0*nloc(1)*4.0d-16*Tl**3 - cool_c*ed
            ba_blk(nt,nt-1,k) = -cool_c*(Tl - 3.0d3)
         else
            ba_blk(nt,nt,k) = 1d0
         end if
         call dgetrf(nt, nt, ba_blk(:,:,k), nt, ba_ipiv(:,k), info)
         if(info /= 0) then
            ba_ready = .false.
            return
         end if
      end do
      n_ba = n_ba + 1
      ba_ready = .true.
   end subroutine build_BA

   ! dv := -BA^{-1} Fv (in place on dv, which enters holding Fv)
   subroutine apply_BA(dv)
      implicit none
      real(dp), intent(inout) :: dv(nt,nd)
      integer :: k, info
      external dgetrs
      do k = 1, nd
         call dgetrs('N', nt, 1, ba_blk(:,:,k), nt, ba_ipiv(:,k), &
                     dv(:,k), nt, info)
      end do
      dv = -dv
   end subroutine apply_BA

   !------------------------------------------------------------------
   ! Exact trace solve at fixed majors+J: the trace block is linear
   ! in n_tr -> one small LU per depth.  One-way coupling means the
   ! majors/J are unaffected, so no extra sweep is needed.  Uses the
   ! same J the caller's F eval used (passed in).
   !------------------------------------------------------------------
   subroutine trace_solve(n_cur, J_eval)
      implicit none
      real(dp), intent(inout) :: n_cur(nt,nd)
      real(dp), intent(in)    :: J_eval(nd)
      real(dp) :: Am(ntr,ntr), rhs(ntr), ed
      integer :: k, i, j, ipiv(ntr), info
      external dgesv
      do k = 1, nd
         ed = n_cur(nt-1,k)
         Am = 0d0; rhs = 0d0
         do i = 1, ntr
            do j = 1, ntr
               if(j /= i) Am(i,j) = atr(j,i)
            end do
            Am(i,i) = -sum(atr(i,:)) - tr_stiff*stiff
            do j = 1, nl
               rhs(i) = rhs(i) - n_cur(j,k)*bsrc(i,j)*(J_eval(k) + &
                        0.1d0)
            end do
         end do
         call dgesv(ntr, 1, Am, ntr, ipiv, rhs, ntr, info)
         if(info == 0) then
            do i = 1, ntr
               n_cur(nl+i,k) = max(rhs(i), tr_scale*ntot*1.0d-8)
            end do
         end if
      end do
   end subroutine trace_solve

   !------------------------------------------------------------------
   ! Outer refresh: the real sweep.  J_store <- J_exact(n).
   !------------------------------------------------------------------
   subroutine refresh_L(n_cur)
      implicit none
      real(dp), intent(in) :: n_cur(nt,nd)
      real(dp) :: S(nd), chi(nd), dS(nt), Jn(nd)
      integer :: k, it_used
      do k = 1, nd
         call source_S(n_cur(:,k), n_cur(nt,k), S(k), dS)
         call opac_chi(n_cur(:,k), chi(k))
      end do
      call solve_J(S, chi, Jn, it_used)
      n_sweep = n_sweep + 1
      j_store = Jn
      j_now = Jn
   end subroutine refresh_L

   !------------------------------------------------------------------
   ! Utilities: pops <-> log unknowns, metrics, counters.
   !------------------------------------------------------------------
   subroutine pops_to_x(n_cur, x)
      implicit none
      real(dp), intent(in)  :: n_cur(nt,nd)
      real(dp), intent(out) :: x(neq)
      integer :: k, i
      do k = 1, nd
         do i = 1, nt
            x((k-1)*nt + i) = log(max(n_cur(i,k),p_floor)/ &
                                  max(n0(i,k),p_floor))
         end do
      end do
   end subroutine pops_to_x

   subroutine x_to_pops(x, n_cur)
      implicit none
      real(dp), intent(in)  :: x(neq)
      real(dp), intent(out) :: n_cur(nt,nd)
      integer :: k, i
      real(dp) :: xt
      do k = 1, nd
         do i = 1, nt
            xt = min(max(x((k-1)*nt+i), -x_clip), x_clip)
            n_cur(i,k) = n0(i,k)*exp(xt)
         end do
      end do
   end subroutine x_to_pops

   ! max and RMS of r_i = F_i/G_i (gross-flux normalized residual)
   subroutine max_r(n_cur, mode_true, rmax, rrms)
      implicit none
      real(dp), intent(in) :: n_cur(nt,nd)
      integer,  intent(in) :: mode_true
      real(dp), intent(out) :: rmax, rrms
      real(dp) :: F(nt,nd), G(nt,nd)
      integer :: i, k, ncnt
      call eval_F(n_cur, F, G, mode_true)
      rmax = 0d0; rrms = 0d0; ncnt = 0
      do k = 1, nd
         do i = 1, nt
            if(G(i,k) > 1.0d-200) then
               rmax = max(rmax, abs(F(i,k))/G(i,k))
               rrms = rrms + (F(i,k)/G(i,k))**2
               ncnt = ncnt + 1
            end if
         end do
      end do
      rrms = sqrt(rrms/max(1,ncnt))
   end subroutine max_r

   subroutine counters_reset()
      implicit none
      n_sweep = 0; n_ba = 0; n_cheap = 0; n_jvp = 0
      n_clip = 0; n_jac_tot = 0; n_jac_evals = 0
   end subroutine counters_reset

end module toy_se_mod
