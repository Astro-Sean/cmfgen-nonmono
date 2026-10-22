! se_prob_mod: surrogate statistical-equilibrium problem for the
! cmfgen_jfnk unit tests.  Mirrors the frozen-J residual used by
! JFNK_RESIDUAL_V1 (new_main/mod_subs/jfnk_cmfgen_interface_v1.f):
!
!   F_i(x) = sum_j A_ij n0_j exp(x_j) + sum_j C_ij x_j + r_i
!            [+ (x_i - x0_i) ptc_scale_i / dtau]
!
! plus the production safeguards that matter to solver behavior:
! exp argument clipping, +/-1D20 residual cap, fixed-row masking
! (T/Ne), pseudo-transient term, and a faithful replica of the
! JFNK_PRECOND_V1 depth-block LU apply chain (DGEEQU equilibration,
! row regularization + z=0 projection, exp(-x) post-rescale, scalar
! PTC diagonal correction).  The exp(-x) and PTC-fix stages are
! switchable so the tests can A/B them.
!
module se_prob_mod
   implicit none
   save
   integer :: se_n = 0, se_nt = 0, se_nd = 0
   double precision, allocatable :: a_se(:,:), c_se(:,:)
   double precision, allocatable :: r_se(:), n0_se(:)
   integer, allocatable :: mask_se(:)

   ! switches / parameters (set from CLI args)
   logical          :: se_expneg = .true.    ! z *= exp(-x) post-solve
   logical          :: se_ptcfix = .true.    ! z *= d/(d+sig) post-solve
   logical          :: se_eqscale = .false.  ! row-equilibrate equations
   double precision :: se_rowtiny = 0d0      ! row-reg threshold
   double precision :: se_ridge = 0d0        ! base LM ridge
   double precision :: se_eps = 0d0          ! FD JVP step scale

   ! FIX_T_AUTO replica state.  In the surrogate the T row is flat
   ! index i with (i-1) mod se_nt = se_trow-1 (row 35, 0-based, in
   ! gen_se_system.py); the Ne row (36) stays masked under tauto.
   ! When se_tauto_on the T rows stay live in F/JVP and the
   ! preconditioned output T components get scaled by se_tau_alpha
   ! (the JFNK_TAU_PCW analogue); step damping itself lives in the
   ! production solver module via jfnk_tauto_*.
   logical          :: se_tauto_on = .false.
   logical          :: se_tau_pcw = .true.
   double precision :: se_tau_alpha = 1d0
   integer          :: se_trow = 0
   ! weighted inner product (production JFNK_RSCALE analog): se_dp
   ! returns sum (w*a)*(w*b).  Emulates SQRT_STE with w = 1/sqrt|F_base|
   ! loaded from wscale.txt -- metric only, equations untouched.
   logical          :: se_wscale_on = .false.
   double precision, allocatable :: se_wscale(:)

   ! structure-flip row: surrogate for GENERATE_FULL_MATRIX_V3's
   ! REPLACE decision.  Row se_flip_row carries the rate equation
   ! while x(se_flip_col) <= se_flip_thr and a conservation form
   ! sum(w_j n_j) - tgt above it.  se_flip_freeze pins the branch to
   ! se_flip_branch0, captured at the baseline x: perturbed evals
   ! replay the baseline equation set, the USE_PASSED_REP=T analogue,
   ! so an FD step cannot cross the flip.
   logical :: se_flip_on = .false., se_flip_freeze = .false.
   logical :: se_flip_branch0 = .false.
   integer :: se_flip_row = 0, se_flip_col = 0, se_flip_ng = 0
   double precision :: se_flip_thr = 0d0, se_flip_tgt = 0d0
   integer, allocatable :: se_flip_grp(:)
   double precision, allocatable :: se_flip_w(:)

   ! Aggregate-factor analogue: live mode multiplies the rate residual
   ! by a factor of the group population sum; lag mode pins that factor
   ! at its baseline value so its Jacobian matches the stored BA.
   logical :: se_factor_on = .false., se_factor_lag = .false.
   integer :: se_factor_ng = 0
   double precision :: se_factor_alpha = 0d0, se_factor_target = 0d0
   integer, allocatable :: se_factor_grp(:)

   ! preconditioner state (per-depth LU, production replica)
   double precision, allocatable :: lu_fac(:,:,:), row_sf(:,:), col_sf(:,:)
   double precision, allocatable :: ba_diag(:,:), pc_diag(:,:), m0_blk(:,:,:)
   ! raw (pre-ridge) block and applied ridge per row, for the
   ! algebraic least-squares apply that never rebuilds B
   double precision, allocatable :: a_raw_blk(:,:,:), ridge_save(:,:)
   double precision, allocatable :: vsave(:)
   integer, allocatable :: lu_ipiv(:,:), lu_info_k(:)
   ! schur pc state: per-depth pop-block Cholesky (nt x nt storage,
   ! leading np x np used) + dense S factor + equilibration
   double precision, allocatable :: se_schur_fac(:,:,:), &
      se_schur_s(:,:), se_rsf2(:,:), se_csf2(:,:)
   integer, allocatable :: se_schur_ipiv(:)
   logical :: se_schur_ok = .false.
   logical, allocatable :: row_reg(:,:), lu_ok(:)
   integer :: se_nreg = 0, se_nridge_depths = 0

contains

   ! True when flat index i is a T row AND tauto keeps it live.
   logical function se_t_live(i)
      integer, intent(in) :: i
      se_t_live = se_tauto_on .and. se_trow > 0 .and. &
                  mod(i - 1, se_nt) == se_trow - 1
   end function se_t_live

end module se_prob_mod


subroutine se_load(dir, stat)
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx
   implicit none
   character(len=*), intent(in) :: dir
   integer, intent(out) :: stat
   integer :: lu, ios, k, i
   character(len=1024) :: fn
   double precision, allocatable :: esc(:)

   stat = 0
   fn = trim(dir) // '/params.txt'
   open(newunit=lu, file=trim(fn), status='old', iostat=ios)
   if (ios /= 0) then
      stat = 2; return
   end if
   read(lu, *) se_nt, se_nd, se_n
   close(lu)
   ! T row is the second-to-last row of each depth block (T_ROW=35,
   ! 0-based, in gen_se_system.py -> Fortran row j = se_nt-1).
   se_trow = se_nt - 1

   allocate(a_se(se_n,se_n), c_se(se_n,se_n), r_se(se_n), n0_se(se_n), &
            mask_se(se_n), vsave(se_n), esc(se_n), stat=ios)
   if (ios /= 0) then
      stat = 3; return
   end if

   fn = trim(dir) // '/n0.txt'
   open(newunit=lu, file=trim(fn), status='old'); read(lu,*) n0_se; close(lu)
   fn = trim(dir) // '/acoll.txt'
   open(newunit=lu, file=trim(fn), status='old'); read(lu,*) a_se; close(lu)
   fn = trim(dir) // '/cmat.txt'
   open(newunit=lu, file=trim(fn), status='old'); read(lu,*) c_se; close(lu)
   fn = trim(dir) // '/rvec.txt'
   open(newunit=lu, file=trim(fn), status='old'); read(lu,*) r_se; close(lu)
   fn = trim(dir) // '/mask.txt'
   open(newunit=lu, file=trim(fn), status='old'); read(lu,*) mask_se; close(lu)

   ! eqscale: row-equilibrate the equations, F' = D^-1 F.  Scaling A,
   ! C and r before any downstream use makes the residual, the JVP
   ! (an FD of F') and the preconditioner matrix (built from scaled
   ! rows) all consistent -- equivalent to scaling inside
   ! JFNK_RESIDUAL_V1.  The root x* is unchanged.
   if (se_eqscale) then
      fn = trim(dir) // '/escale.txt'
      open(newunit=lu, file=trim(fn), status='old')
      read(lu,*) esc
      close(lu)
      do k = 1, se_n
         if (esc(k) > 0d0) then
            a_se(k,:) = a_se(k,:) / esc(k)
            c_se(k,:) = c_se(k,:) / esc(k)
            r_se(k)   = r_se(k)   / esc(k)
         end if
      end do
   end if
   deallocate(esc)

   if (se_wscale_on) then
      allocate(se_wscale(se_n), stat=ios)
      if (ios /= 0) then
         stat = 4; return
      end if
      se_wscale = 1d0
      fn = trim(dir) // '/wscale.txt'
      open(newunit=lu, file=trim(fn), status='old', iostat=ios)
      if (ios == 0) then
         read(lu, *) se_wscale
         close(lu)
      end if
   end if

   se_flip_on = .false.
   se_factor_on = .false.
   se_factor_lag = .false.

   ! optional structure-flip row spec, emitted only by the 'flip'
   ! system spec; absent elsewhere so the residual stays branch-A.
   fn = trim(dir) // '/flip.txt'
   open(newunit=lu, file=trim(fn), status='old', iostat=ios)
   if (ios == 0) then
      read(lu, *) se_flip_row, se_flip_col, se_flip_thr, se_flip_tgt
      read(lu, *) se_flip_ng
      if (allocated(se_flip_grp)) deallocate(se_flip_grp, se_flip_w)
      allocate(se_flip_grp(se_flip_ng), se_flip_w(se_flip_ng))
      do i = 1, se_flip_ng
         read(lu, *) se_flip_grp(i), se_flip_w(i)
      end do
      close(lu)
      se_flip_on = .true.
   end if

   fn = trim(dir) // '/factor.txt'
   open(newunit=lu, file=trim(fn), status='old', iostat=ios)
   if (ios == 0) then
      read(lu, *) se_factor_alpha, se_factor_target, se_factor_ng
      if (allocated(se_factor_grp)) deallocate(se_factor_grp)
      allocate(se_factor_grp(se_factor_ng))
      read(lu, *) se_factor_grp
      close(lu)
      se_factor_on = .true.
   end if

   if (allocated(jfnk_x)) deallocate(jfnk_x)
   if (allocated(jfnk_fx)) deallocate(jfnk_fx)
   allocate(jfnk_x(se_n), jfnk_fx(se_n))
   jfnk_x = 0d0
   jfnk_fx = 0d0
   i = 0
end subroutine se_load


subroutine se_f_raw(n, x, fx)
   ! Raw residual: F = A (n0 e^x) + C x + r, exp arg clipped at +-80,
   ! components capped at +-1D20 (production NaN/Inf policy), masked
   ! rows forced to 0.  No PTC term -- se_f adds it when armed.
   use se_prob_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: fx(n)
   double precision :: npop(n), an(n), factor
   integer :: i
   logical :: fb
   npop = n0_se * exp(min(80d0, max(-80d0, x)))
   an = matmul(a_se, npop)
   factor = 1d0
   if (se_factor_on .and. .not. se_factor_lag) then
      factor = 1d0 + se_factor_alpha * &
         (sum(npop(se_factor_grp)) - se_factor_target)
   end if
   fx = factor * an + matmul(c_se, x) + r_se
   if (se_flip_on .and. se_flip_row > 0) then
      ! branch B is the conservation equation; it need not share the
      ! rate equation's value at the flip point, which is what makes
      ! a live (unfrozen) branch discontinuous in the FD sense.
      if (se_flip_freeze) then
         fb = se_flip_branch0
      else
         fb = (x(se_flip_col) > se_flip_thr)
      end if
      if (fb) fx(se_flip_row) = sum(se_flip_w * npop(se_flip_grp)) &
                              - se_flip_tgt
   end if
   do i = 1, n
      if (mask_se(i) == 1 .and. .not. se_t_live(i)) fx(i) = 0d0
      if (fx(i) /= fx(i)) then
         fx(i) = 1d20
      else if (abs(fx(i)) > 1d20) then
         fx(i) = sign(1d20, fx(i))
      end if
   end do
end subroutine se_f_raw


subroutine se_f(n, x, fx)
   ! Residual callback for newton_hook_v1: raw F plus the production
   ! pseudo-time term (x - x0)*s/dtau read from the solver module, so
   ! re-anchoring inside newton_hook_v1 works exactly as in CMFGEN.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_scale, &
      jfnk_ptc_scale_ready, jfnk_ptc_x0, jfnk_ptc_anchor_ready, &
      jfnk_ptc_dtau
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: fx(n)
   integer :: i
   call se_f_raw(n, x, fx)
   if (jfnk_ptc_active .and. jfnk_ptc_scale_ready) then
      if (jfnk_ptc_anchor_ready .and. allocated(jfnk_ptc_x0)) then
         do i = 1, n
            fx(i) = fx(i) + (x(i) - jfnk_ptc_x0(i)) * jfnk_ptc_scale(i) &
                  / max(jfnk_ptc_dtau, 1d-8)
         end do
      else
         do i = 1, n
            fx(i) = fx(i) + x(i) * jfnk_ptc_scale(i) / &
                    max(jfnk_ptc_dtau, 1d-8)
         end do
      end if
   end if
end subroutine se_f


subroutine se_jexact(n, x, v, jv)
   ! Analytic J(x) v with J = A diag(n0 e^x) + C (+ ptc diag if armed).
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_scale, &
      jfnk_ptc_scale_ready, jfnk_ptc_dtau
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n), v(n)
   double precision, intent(out) :: jv(n)
   double precision :: npop(n), an(n), factor, dfactor
   integer :: i
   logical :: fb
   npop = n0_se * exp(min(80d0, max(-80d0, x)))
   an = matmul(a_se, npop)
   factor = 1d0
   dfactor = 0d0
   if (se_factor_on .and. .not. se_factor_lag) then
      factor = 1d0 + se_factor_alpha * &
         (sum(npop(se_factor_grp)) - se_factor_target)
      dfactor = se_factor_alpha * &
         sum(npop(se_factor_grp) * v(se_factor_grp))
   end if
   jv = factor * matmul(a_se, npop * v) + matmul(c_se, v) &
      + dfactor * an
   if (se_flip_on .and. se_flip_row > 0) then
      ! the exact J is taken w.r.t. whichever branch the residual is
      ! actually evaluating, so live-vs-frozen comparisons measure
      ! the FD artifact, not a model mismatch.
      if (se_flip_freeze) then
         fb = se_flip_branch0
      else
         fb = (x(se_flip_col) > se_flip_thr)
      end if
      if (fb) jv(se_flip_row) = sum(se_flip_w * n0_se(se_flip_grp) * &
         exp(min(80d0, max(-80d0, x(se_flip_grp)))) * &
         v(se_flip_grp))
   end if
   if (jfnk_ptc_active .and. jfnk_ptc_scale_ready) then
      jv = jv + jfnk_ptc_scale * v / max(jfnk_ptc_dtau, 1d-8)
   end if
   do i = 1, n
      if (mask_se(i) == 1 .and. .not. se_t_live(i)) jv(i) = 0d0
   end do
end subroutine se_jexact


subroutine se_jvp_fd(n, v, jv)
   ! Production FD rule (JFNK_JVP_V1): h = eps*max(1,|x|inf)/|v|inf,
   ! forward difference of se_f at the module state jfnk_x.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_x
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   double precision :: h, xp(n), fp(n), f0(n)
   if (maxval(abs(v)) < 1d-300) then
      jv = 0d0
      return
   end if
   h = se_eps * max(1d0, maxval(abs(jfnk_x))) / maxval(abs(v))
   xp = jfnk_x + h * v
   call se_f(n, jfnk_x, f0)
   call se_f(n, xp, fp)
   jv = (fp - f0) / h
end subroutine se_jvp_fd


subroutine se_jvp_fd_stale(n, v, jv)
   ! Pre-fix production semantics for the regression test: the FD
   ! subtracts the stored baseline jfnk_fx instead of a fresh eval at
   ! jfnk_x, so any offset between the stored and evaluated residual
   ! enters every matvec amplified by 1/h.  Kept deliberately wrong.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   double precision :: h, xp(n), fp(n)
   if (maxval(abs(v)) < 1d-300) then
      jv = 0d0
      return
   end if
   h = se_eps * max(1d0, maxval(abs(jfnk_x))) / maxval(abs(v))
   xp = jfnk_x + h * v
   call se_f(n, xp, fp)
   jv = (fp - jfnk_fx) / h
end subroutine se_jvp_fd_stale


subroutine se_factor_lu()
   ! Replica of JFNK_FACTOR_PRECOND_V1 (LU path): per-depth block of
   ! the baseline Jacobian J0 = A diag(n0) + C, ridge retry ladder,
   ! near-null row -> identity (rowmax < rowtiny*blockmax), then
   ! DGEEQU + DGETRF.  Also stores the unaugmented block and the
   ! effective M block (with ridge + identity rows) for diagnostics.
   use se_prob_mod
   implicit none
   integer :: k, i, j, info, itry
   double precision :: blkmax, rowmax, ridadd, rowcnd, colcnd, amax
   double precision, allocatable :: wrk(:,:)
   double precision :: tr(5) = (/0d0, 1d-12, 1d-8, 1d-4, 1d-2/)

   if (allocated(lu_fac)) deallocate(lu_fac, row_sf, col_sf, ba_diag, &
      pc_diag, m0_blk, lu_ipiv, lu_info_k, row_reg, lu_ok, &
      a_raw_blk, ridge_save)
   allocate(lu_fac(se_nt,se_nt,se_nd), row_sf(se_nt,se_nd), &
            col_sf(se_nt,se_nd), ba_diag(se_nt,se_nd), &
            pc_diag(se_nt,se_nd), m0_blk(se_nt,se_nt,se_nd), &
            a_raw_blk(se_nt,se_nt,se_nd), ridge_save(se_nt,se_nd), &
            lu_ipiv(se_nt,se_nd), lu_info_k(se_nd), &
            row_reg(se_nt,se_nd), lu_ok(se_nd), wrk(se_nt,se_nt))
   se_nreg = 0
   se_nridge_depths = 0

   do k = 1, se_nd
      ! baseline x-space Jacobian block (POPS multiplication is built
      ! into BA via GENERATE_FULL_MATRIX_V3 C_MAT(I,J)*POPS(J))
      do i = 1, se_nt
         do j = 1, se_nt
            m0_blk(i,j,k) = a_se((k-1)*se_nt+i, (k-1)*se_nt+j) * &
                            n0_se((k-1)*se_nt+j) + &
                            c_se((k-1)*se_nt+i, (k-1)*se_nt+j)
         end do
         ba_diag(i,k) = m0_blk(i,i,k)
         pc_diag(i,k) = m0_blk(i,i,k)
      end do
      a_raw_blk(:,:,k) = m0_blk(:,:,k)
      ridge_save(:,k) = 0d0
      blkmax = maxval(abs(m0_blk(:,:,k)))
      row_reg(:,k) = .false.
      do i = 1, se_nt
         rowmax = maxval(abs(m0_blk(i,:,k)))
         if (blkmax > 0d0 .and. rowmax < se_rowtiny*blkmax) then
            row_reg(i,k) = .true.
            se_nreg = se_nreg + 1
         end if
      end do

      lu_ok(k) = .false.
      do itry = 1, 5
         ridadd = max(tr(itry), se_ridge)
         wrk = m0_blk(:,:,k)
         do i = 1, se_nt
            wrk(i,i) = wrk(i,i) + ridadd * &
                       max(abs(m0_blk(i,i,k)), 1d-30*blkmax)
         end do
         do i = 1, se_nt
            if (row_reg(i,k)) then
               wrk(i,:) = 0d0
               wrk(i,i) = 1d0
            end if
         end do
         lu_fac(:,:,k) = wrk
         call dgeequ(se_nt, se_nt, lu_fac(1,1,k), se_nt, &
                     row_sf(1,k), col_sf(1,k), rowcnd, colcnd, amax, info)
         if (info /= 0) cycle
         do j = 1, se_nt
            do i = 1, se_nt
               lu_fac(i,j,k) = lu_fac(i,j,k) * row_sf(i,k) * col_sf(j,k)
            end do
         end do
         call dgetrf(se_nt, se_nt, lu_fac(1,1,k), se_nt, &
                     lu_ipiv(1,k), info)
         if (info == 0) then
            lu_ok(k) = .true.
            lu_info_k(k) = itry - 1
            do i = 1, se_nt
               pc_diag(i,k) = m0_blk(i,i,k) + ridadd * &
                              max(abs(m0_blk(i,i,k)), 1d-30*blkmax)
               ridge_save(i,k) = ridadd * &
                                 max(abs(m0_blk(i,i,k)), 1d-30*blkmax)
               if (row_reg(i,k)) pc_diag(i,k) = 1d0
            end do
            if (ridadd > 0d0) se_nridge_depths = se_nridge_depths + 1
            ! rebuild effective M for diagnostics: ridged + reg rows,
            ! unequilibrated
            m0_blk(:,:,k) = wrk
            exit
         end if
      end do
   end do
   deallocate(wrk)
end subroutine se_factor_lu


subroutine se_psolve(n, v)
   ! Replica of JFNK_PRECOND_V1 LU apply chain:
   !   1. v' = ROW_SF v;  2. DGETRS;  3. z = COL_SF z';
   !   4. regularized rows -> z = 0
   !   5. [switch] z *= exp(-x): the production log-unknown rescale
   !   6. [switch] scalar PTC fix: z *= d'/(d'+sig), sig = scale/dtau,
   !      with the null-diagonal branch z = v/sig
   !   7. mask fixed unknowns, nonfinite -> 0 + clip event
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_ptc_active, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_dtau, &
      jfnk_pc_clip_event
   implicit none
   external :: jfnk_schur_apply_v1
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: k, i, idx, info
   double precision :: zba, xpc, dpc, sigpc
   vsave = v
   do k = 1, se_nd
      if (lu_ok(k)) then
         do i = 1, se_nt
            idx = (k-1)*se_nt + i
            v(idx) = v(idx) * row_sf(i,k)
         end do
         call dgetrs('N', se_nt, 1, lu_fac(1,1,k), se_nt, &
                     lu_ipiv(1,k), v((k-1)*se_nt+1), se_nt, info)
         do i = 1, se_nt
            idx = (k-1)*se_nt + i
            v(idx) = v(idx) * col_sf(i,k)
            if (row_reg(i,k)) v(idx) = 0d0
         end do
      else
         do i = 1, se_nt
            idx = (k-1)*se_nt + i
            if (row_reg(i,k)) then
               v(idx) = 0d0
            else if (abs(pc_diag(i,k)) > 1d-30) then
               v(idx) = v(idx) / pc_diag(i,k)
            else
               v(idx) = 0d0
            end if
         end do
      end if
   end do
   do idx = 1, n
      k = (idx-1)/se_nt + 1
      i = idx - (k-1)*se_nt
      if (allocated(row_reg) .and. row_reg(i,k)) then
         v(idx) = 0d0
         cycle
      end if
      zba = v(idx)
      if (se_expneg) then
         xpc = min(max(jfnk_x(idx), -40d0), 40d0)
         zba = zba * exp(-xpc)
      end if
      if (se_ptcfix .and. jfnk_ptc_active .and. jfnk_ptc_scale_ready) then
         dpc = pc_diag(i,k)
         sigpc = jfnk_ptc_scale(idx) / max(jfnk_ptc_dtau, 1d-8)
         if (abs(dpc) > 1d-30) then
            if (abs(dpc+sigpc) > 1d-30*max(1d0, abs(dpc))) then
               v(idx) = zba * dpc / (dpc + sigpc)
            else
               v(idx) = zba
            end if
         else if (abs(sigpc) > 1d-30) then
            v(idx) = vsave(idx) / sigpc
         else
            v(idx) = zba
         end if
      else
         v(idx) = zba
      end if
      if (mask_se(idx) == 1) then
         if (se_t_live(idx)) then
            if (se_tau_pcw) v(idx) = v(idx) * se_tau_alpha
         else
            v(idx) = 0d0
         end if
      end if
      if (v(idx) /= v(idx)) then
         v(idx) = 0d0
         jfnk_pc_clip_event = .true.
      end if
   end do
end subroutine se_psolve


subroutine se_psolve_jacobi(n, v)
   ! Jacobi path: v/diag(BA+ridge) then the same switchable post
   ! chain as the LU path.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_ptc_active, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_dtau, &
      jfnk_pc_clip_event
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: k, i, idx
   double precision :: zba, xpc, dpc, sigpc
   vsave = v
   do k = 1, se_nd
      do i = 1, se_nt
         idx = (k-1)*se_nt + i
         if (row_reg(i,k)) then
            v(idx) = 0d0
         else if (abs(pc_diag(i,k)) > 1d-30) then
            v(idx) = v(idx) / pc_diag(i,k)
         else
            v(idx) = 0d0
         end if
      end do
   end do
   do idx = 1, n
      k = (idx-1)/se_nt + 1
      i = idx - (k-1)*se_nt
      if (row_reg(i,k)) then
         v(idx) = 0d0
         cycle
      end if
      zba = v(idx)
      if (se_expneg) then
         xpc = min(max(jfnk_x(idx), -40d0), 40d0)
         zba = zba * exp(-xpc)
      end if
      if (se_ptcfix .and. jfnk_ptc_active .and. &
          jfnk_ptc_scale_ready) then
         dpc = pc_diag(i,k)
         sigpc = jfnk_ptc_scale(idx) / max(jfnk_ptc_dtau, 1d-8)
         if (abs(dpc) > 1d-30) then
            if (abs(dpc+sigpc) > 1d-30*max(1d0, abs(dpc))) then
               v(idx) = zba * dpc / (dpc + sigpc)
            else
               v(idx) = zba
            end if
         else if (abs(sigpc) > 1d-30) then
            v(idx) = vsave(idx) / sigpc
         else
            v(idx) = zba
         end if
      else
         v(idx) = zba
      end if
      if (mask_se(idx) == 1) then
         if (se_t_live(idx)) then
            if (se_tau_pcw) v(idx) = v(idx) * se_tau_alpha
         else
            v(idx) = 0d0
         end if
      end if
      if (v(idx) /= v(idx)) then
         v(idx) = 0d0
         jfnk_pc_clip_event = .true.
      end if
   end do
end subroutine se_psolve_jacobi


subroutine se_psolve_tsvd(n, v)
   ! Truncated-SVD (Tikhonov-type) preconditioner candidate: per
   ! depth block solve M z = v by DGELSD with relative cutoff
   ! se_ridge.  Singular directions get the minimum-norm answer
   ! instead of the LU blow-up -> z stays finite and points along the
   ! solvable part of the residual.  Reg-row projection and the fixed
   ! mask still apply.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_pc_clip_event
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: k, i, idx, info, lwork, liwork, rank
   double precision, allocatable :: u(:,:), s(:), svv(:,:), wk(:)
   integer, allocatable :: iwk(:)
   double precision :: rcond, wkq(1)
   integer :: iwkq(1)
   vsave = v
   ! workspace query (lwork = -1): dgelsd returns optimal sizes.
   ! WORK and IWORK are rank-1 arrays in LAPACK even on a query call;
   ! passing scalars is a rank mismatch under gfortran.
   allocate(u(se_nt,se_nt), s(se_nt), svv(se_nt,se_nt), stat=info)
   if (info /= 0) then
      v = 0d0
      return
   end if
   u = m0_blk(:,:,1)
   wkq = 0d0
   iwkq = 0
   call dgelsd(se_nt, se_nt, 1, u, se_nt, v(1), se_nt, s, 1d-8, &
               rank, wkq, -1, iwkq, info)
   ! fall back to documented bounds if the query failed or returned
   ! unusable sizes, so a bad query cannot size the real call
   lwork = int(wkq(1)) + 100
   liwork = iwkq(1) + 100
   if (info /= 0 .or. lwork < 12*se_nt+2*se_nt*se_nt .or. &
       liwork < 12*se_nt) then
      lwork = 12*se_nt + 4*se_nt*se_nt + 512
      liwork = 3*se_nt*max(1, int(log(dble(se_nt))/log(2d0)) + 1) + &
               11*se_nt
   end if
   allocate(wk(lwork), iwk(liwork), stat=info)
   if (info /= 0) then
      deallocate(u, s, svv)
      v = 0d0
      return
   end if
   rcond = max(se_ridge, 1d-14)
   do k = 1, se_nd
      u = m0_blk(:,:,k)
      call dgelsd(se_nt, se_nt, 1, u, se_nt, &
                  v((k-1)*se_nt+1), se_nt, s, rcond, rank, &
                  wk, lwork, iwk, info)
      if (info /= 0) then
         do i = 1, se_nt
            idx = (k-1)*se_nt + i
            v(idx) = 0d0
         end do
         jfnk_pc_clip_event = .true.
         cycle
      end if
      do i = 1, se_nt
         idx = (k-1)*se_nt + i
         if (row_reg(i,k)) v(idx) = 0d0
         if (mask_se(idx) == 1) then
         if (se_t_live(idx)) then
            if (se_tau_pcw) v(idx) = v(idx) * se_tau_alpha
         else
            v(idx) = 0d0
         end if
      end if
         if (v(idx) /= v(idx)) then
            v(idx) = 0d0
            jfnk_pc_clip_event = .true.
         end if
      end do
   end do
   deallocate(u, s, svv, wk, iwk)
end subroutine se_psolve_tsvd


subroutine se_psolve_ptclu(n, v)
   ! PTC-consistent preconditioner candidate: factor
   ! M = J0_blk + diag(scale/dtau) and solve -- the damped operator
   ! GMRES actually sees when PTC is armed, instead of correcting the
   ! BA solve diagonally after the fact.  Reg rows keep their z=0
   ! projection; fixed mask applies.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_scale, &
      jfnk_ptc_scale_ready, jfnk_ptc_dtau, jfnk_pc_clip_event
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: k, i, idx, info
   double precision, allocatable :: wm(:,:)
   integer, allocatable :: ip(:)
   double precision :: sigpc
   vsave = v
   allocate(wm(se_nt,se_nt), ip(se_nt))
   do k = 1, se_nd
      wm = m0_blk(:,:,k)
      if (jfnk_ptc_active .and. jfnk_ptc_scale_ready) then
         do i = 1, se_nt
            idx = (k-1)*se_nt + i
            sigpc = jfnk_ptc_scale(idx) / max(jfnk_ptc_dtau, 1d-8)
            wm(i,i) = wm(i,i) + sigpc
         end do
      end if
      do i = 1, se_nt
         if (row_reg(i,k)) then
            wm(i,:) = 0d0
            wm(i,i) = 1d0
         end if
      end do
      call dgesv(se_nt, 1, wm, se_nt, ip, v((k-1)*se_nt+1), se_nt, info)
      do i = 1, se_nt
         idx = (k-1)*se_nt + i
         if (row_reg(i,k)) v(idx) = 0d0
         if (mask_se(idx) == 1) then
         if (se_t_live(idx)) then
            if (se_tau_pcw) v(idx) = v(idx) * se_tau_alpha
         else
            v(idx) = 0d0
         end if
      end if
         if (v(idx) /= v(idx)) then
            v(idx) = 0d0
            jfnk_pc_clip_event = .true.
         end if
      end do
   end do
   deallocate(wm, ip)
end subroutine se_psolve_ptclu


subroutine se_psolve_lsq(n, v)
   ! Ridge least-squares candidate: solve z = (J^T J + lam^2 I)^{-1}
   ! J^T v per depth block via DPOTRF.  Near-null directions are
   ! damped by lam instead of exploding as in the LU solve, and the
   ! normal equations factor needs the same storage as the LU it
   ! would replace -- unlike TSVD it keeps no U/VT.  lam =
   ! se_ridge * max|J| so se_ridge plays the rcond role.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_pc_clip_event
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: k, i, idx, info
   double precision, allocatable :: gn(:,:), rhs(:)
   double precision :: lam, blkmax
   allocate(gn(se_nt,se_nt), rhs(se_nt))
   do k = 1, se_nd
      blkmax = maxval(abs(m0_blk(:,:,k)))
      lam = max(se_ridge, 1d-14) * max(blkmax, 1d-300)
      call dgemm('T', 'N', se_nt, se_nt, se_nt, 1d0, &
                 m0_blk(:,:,k), se_nt, m0_blk(:,:,k), se_nt, &
                 0d0, gn, se_nt)
      call dgemv('T', se_nt, se_nt, 1d0, m0_blk(:,:,k), se_nt, &
                 v((k-1)*se_nt+1), 1, 0d0, rhs, 1)
      do i = 1, se_nt
         gn(i,i) = gn(i,i) + lam*lam
      end do
      call dpotrf('L', se_nt, gn, se_nt, info)
      if (info == 0) then
         call dpotrs('L', se_nt, 1, gn, se_nt, rhs, se_nt, info)
      end if
      do i = 1, se_nt
         idx = (k-1)*se_nt + i
         if (info /= 0) then
            v(idx) = 0d0
            jfnk_pc_clip_event = .true.
         else
            v(idx) = rhs(i)
         end if
         if (row_reg(i,k)) v(idx) = 0d0
         if (mask_se(idx) == 1) then
         if (se_t_live(idx)) then
            if (se_tau_pcw) v(idx) = v(idx) * se_tau_alpha
         else
            v(idx) = 0d0
         end if
      end if
         if (v(idx) /= v(idx)) then
            v(idx) = 0d0
            jfnk_pc_clip_event = .true.
         end if
      end do
   end do
   deallocate(gn, rhs)
end subroutine se_psolve_lsq


subroutine se_psolve_lsq2(n, v)
   ! Algebraic equivalent of se_psolve_lsq, written the way the
   ! production apply will be: B = A_raw + D_ridge + reg-row
   ! substitution is never rebuilt.  Instead
   !   B^T w = A_raw^T w + D_r w + sum_{i reg} w_i (e_i - A_raw(i,:)^T)
   ! uses DGEMV('T') on the stored raw block plus rank-1 corrections.
   ! The test suite cross-checks |z_lsq2 - z_lsq| to confirm the
   ! formula before it is mirrored into the interface.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_pc_clip_event
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: k, i, j, idx, info
   double precision, allocatable :: gn(:,:), rhs(:), w(:)
   double precision :: lam, blkmax, wi
   allocate(gn(se_nt,se_nt), rhs(se_nt), w(se_nt))
   do k = 1, se_nd
      blkmax = maxval(abs(m0_blk(:,:,k)))
      lam = max(se_ridge, 1d-14) * max(blkmax, 1d-300)
      call dgemm('T', 'N', se_nt, se_nt, se_nt, 1d0, &
                 m0_blk(:,:,k), se_nt, m0_blk(:,:,k), se_nt, &
                 0d0, gn, se_nt)
      ! w = v for this block; rhs = B^T w via the algebraic split.
      ! m0_blk already holds B (ridged + reg-substituted), so the
      ! formula is exercised on the RAW block a0 stored in a_raw_blk.
      w = v((k-1)*se_nt+1:k*se_nt)
      call dgemv('T', se_nt, se_nt, 1d0, a_raw_blk(:,:,k), se_nt, &
                 w, 1, 0d0, rhs, 1)
      do i = 1, se_nt
         rhs(i) = rhs(i) + ridge_save(i,k) * w(i)
      end do
      do i = 1, se_nt
         if (row_reg(i,k)) then
            wi = w(i)
            do j = 1, se_nt
               rhs(j) = rhs(j) - wi * a_raw_blk(i,j,k)
            end do
            ! reg row i in B is e_i^T, so its ridge and raw-row
            ! contributions are replaced outright
            rhs(i) = rhs(i) + wi - wi * ridge_save(i,k)
         end if
      end do
      do i = 1, se_nt
         gn(i,i) = gn(i,i) + lam*lam
      end do
      call dpotrf('L', se_nt, gn, se_nt, info)
      if (info == 0) then
         call dpotrs('L', se_nt, 1, gn, se_nt, rhs, se_nt, info)
      end if
      do i = 1, se_nt
         idx = (k-1)*se_nt + i
         if (info /= 0) then
            v(idx) = 0d0
            jfnk_pc_clip_event = .true.
         else
            v(idx) = rhs(i)
         end if
         if (row_reg(i,k)) v(idx) = 0d0
         if (mask_se(idx) == 1) then
         if (se_t_live(idx)) then
            if (se_tau_pcw) v(idx) = v(idx) * se_tau_alpha
         else
            v(idx) = 0d0
         end if
      end if
         if (v(idx) /= v(idx)) then
            v(idx) = 0d0
            jfnk_pc_clip_event = .true.
         end if
      end do
   end do
   deallocate(gn, rhs, w)
end subroutine se_psolve_lsq2


subroutine se_psolve_ident(n, v)
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   v = v
end subroutine se_psolve_ident


double precision function se_dp(n, a, b)
   use se_prob_mod, only: se_wscale, se_wscale_on
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   if (se_wscale_on .and. allocated(se_wscale)) then
      se_dp = dot_product(se_wscale * a, se_wscale * b)
   else
      se_dp = dot_product(a, b)
   end if
end function se_dp


subroutine se_sub()
   use jfnk_solver_v1_mod, only: jfnk_tol, jfnk_ptc_dtau
   implicit none
   write(*, '(A,2(1X,ES16.8))') 'SUB ', jfnk_tol, jfnk_ptc_dtau
end subroutine se_sub


subroutine se_factor_schur()
   ! Surrogate replica of the SCHUR factor path: per-depth population
   ! sub-block (rows 1..nt-2) with ridge + reg-row substitution,
   ! DGEEQU equilibration, normal-equations DPOTRF stored in
   ! se_schur_fac.  S assembled by the PRODUCTION routine
   ! jfnk_schur_form_v1 (jfnk_solver_v1.f90).  Surrogate row map:
   ! T at se_trow = nt-1, Ne at nt (production has them swapped --
   ! the jt/jne arguments absorb the difference).
   use se_prob_mod
   implicit none
   external :: jfnk_schur_form_v1
   integer :: k, i, j, np, info, ierr
   double precision, allocatable :: b0(:,:), gn(:,:)
   double precision :: lam, blkmax, rcnd, ccnd, amax
   np = se_nt - 2
   call se_factor_lu()
   if (allocated(se_schur_fac)) deallocate(se_schur_fac, se_schur_s, &
      se_schur_ipiv, se_rsf2, se_csf2)
   allocate(se_schur_fac(se_nt,se_nt,se_nd), se_schur_s(se_nd,se_nd), &
            se_schur_ipiv(se_nd), se_rsf2(se_nt,se_nd), &
            se_csf2(se_nt,se_nd))
   se_rsf2 = 1d0
   se_csf2 = 1d0
   se_schur_fac = 0d0
   do k = 1, se_nd
      allocate(b0(np,np), gn(np,np))
      b0 = a_raw_blk(1:np,1:np,k)
      do i = 1, np
         b0(i,i) = b0(i,i) + ridge_save(i,k)
      end do
      do i = 1, np
         if (row_reg(i,k)) then
            b0(i,:) = 0d0
            b0(i,i) = 1d0
         end if
      end do

      blkmax = max(maxval(abs(b0)), 1d-300)
      lam = max(se_ridge, 1d-14) * blkmax
      call dsyrk('L', 'T', np, np, 1d0, b0, np, 0d0, gn, np)
      do i = 1, np
         gn(i,i) = gn(i,i) + lam*lam
      end do
      call dpotrf('L', np, gn, np, info)
      if (info /= 0) then
         ! fall back to a diagonal block so the apply stays finite
         gn = 0d0
         do i = 1, np
            gn(i,i) = max(abs(b0(i,i)), lam*lam)
         end do
         call dpotrf('L', np, gn, np, info)
      end if
      se_schur_fac(1:np,1:np,k) = gn
      deallocate(b0, gn)
   end do
   call jfnk_schur_form_v1(se_nd, se_nt, np, se_trow, &
        se_schur_fac, a_raw_blk, ridge_save, row_reg, &
        se_rsf2, se_csf2, se_schur_s, se_schur_ipiv, ierr)
   se_schur_ok = (ierr == 0)
end subroutine se_factor_schur


subroutine se_psolve_schur(n, v)
   ! Two-level (n,T) Schur apply via the PRODUCTION routine
   ! jfnk_schur_apply_v1, then the same post chain as se_psolve:
   ! reg projection, exp(-x) rescale, PTC fix, mask, clip flag.
   use se_prob_mod
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_ptc_active, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_dtau, &
      jfnk_pc_clip_event
   implicit none
   external :: jfnk_schur_apply_v1
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: k, i, idx, info
   double precision :: zba, xpc, dpc, sigpc
   vsave = v
   if (.not. se_schur_ok) then
      call se_factor_schur()
   end if
   call jfnk_schur_apply_v1(se_nd, se_nt, se_nt-2, se_trow, se_nt, &
        se_schur_fac, a_raw_blk, ridge_save, row_reg, &
        se_rsf2, se_csf2, se_schur_s, se_schur_ipiv, v, info)
   if (info /= 0) then
      v = 0d0
      jfnk_pc_clip_event = .true.
      return
   end if
   do idx = 1, n
      k = (idx-1)/se_nt + 1
      i = idx - (k-1)*se_nt
      if (allocated(row_reg) .and. row_reg(i,k)) then
         v(idx) = 0d0
         cycle
      end if
      zba = v(idx)
      if (se_expneg) then
         xpc = min(max(jfnk_x(idx), -40d0), 40d0)
         zba = zba * exp(-xpc)
      end if
      if (se_ptcfix .and. jfnk_ptc_active .and. &
          jfnk_ptc_scale_ready) then
         dpc = pc_diag(i,k)
         sigpc = jfnk_ptc_scale(idx) / max(jfnk_ptc_dtau, 1d-8)
         if (abs(dpc) > 1d-30) then
            if (abs(dpc+sigpc) > 1d-30*max(1d0, abs(dpc))) then
               v(idx) = zba * dpc / (dpc + sigpc)
            else
               v(idx) = zba
            end if
         else if (abs(sigpc) > 1d-30) then
            v(idx) = vsave(idx) / sigpc
         else
            v(idx) = zba
         end if
      else
         v(idx) = zba
      end if
      if (mask_se(idx) == 1) then
         if (se_t_live(idx)) then
            if (se_tau_pcw) v(idx) = v(idx) * se_tau_alpha
         else
            v(idx) = 0d0
         end if
      end if
      if (v(idx) /= v(idx)) then
         v(idx) = 0d0
         jfnk_pc_clip_event = .true.
      end if
   end do
end subroutine se_psolve_schur
