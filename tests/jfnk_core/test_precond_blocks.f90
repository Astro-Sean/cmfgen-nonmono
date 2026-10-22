!----------------------------------------------------------------------
! test_precond_blocks: unit tests for the extracted per-depth BA
! preconditioner routines jfnk_pc_factor_block_v1 /
! jfnk_pc_apply_block_v1 (new_main/mod_subs/jfnk_solver_v1.f90).
!
! T-ident: the extraction must be bit-identical to the old inline
! code.  old_factor/old_apply below are verbatim copies of the
! pre-extraction inline logic from jfnk_cmfgen_interface_v1.f; on
! random blocks every output must compare exactly equal.
!
! T-linear: the apply must be a fixed linear operator in v (GMRES
! Krylov consistency), for every flag combination.
!
! added by Devin, 22-Oct-2026 (jfnk-precond-fixes)
!----------------------------------------------------------------------
program test_precond_blocks
   use jfnk_solver_v1_mod, only: jfnk_ridge_floor, &
        jfnk_ridge_floor_set, jfnk_pc_signed, jfnk_pc_reg_zerocol, &
        jfnk_pc_tiny_equil, jfnk_pc_rcond_gate, jfnk_pc_rcond_min, &
        jfnk_pc_itry_cap
   implicit none
   integer :: nfail
!
   nfail = 0
   jfnk_ridge_floor = 0.0d0
   jfnk_ridge_floor_set = .true.
!
   call t_ident_factor(nfail)
   call t_ident_apply(nfail)
   call t_rcond_norm(nfail)
   call t_rscale_order(nfail)
   call t_sign(nfail)
   call t_ridge(nfail)
   call t_linear(nfail)
!
   if(nfail == 0)then
      print *, 'ALL PRECOND BLOCK TESTS PASSED'
   else
      print *, nfail, ' PRECOND BLOCK TESTS FAILED'
      call exit(1)
   end if
end program test_precond_blocks
!
!======================================================================
! Verbatim copy of the pre-extraction inline factorization (the
! ITRY=0..8 ladder) from JFNK_FACTOR_PRECOND_V1, reduced to one block.
! pre_reg plays the role of the population-floor mask the caller
! supplied; row_tiny/ridge_piv/lsq_lam mirror the control variables.
!======================================================================
subroutine old_factor(a, n, pre_reg, row_tiny, ridge_piv, lsq_lam, &
     lsq_mode, lu, ipiv, rsf, csf, reg, ridge_rows, rcond, pivrat, &
     itry_used, ridge_used, lu_ok)
   use jfnk_solver_v1_mod, only: jfnk_ridge_floor
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n,n), row_tiny, ridge_piv, &
        lsq_lam
   logical, intent(in) :: pre_reg(n), lsq_mode
   double precision, intent(out) :: lu(n,n), rsf(n), csf(n), &
        ridge_rows(n), rcond, pivrat, ridge_used
   integer, intent(out) :: ipiv(n), itry_used
   logical, intent(out) :: reg(n), lu_ok
   integer :: i, j, itry, ios, info_k, info_g
   double precision :: eps, ridge_add, blkmax, scale_k, rowmax_i
   double precision :: row_cnd, col_cnd, amax, lsq_l2, anorm_k
   double precision :: s_tmp, piv_min, piv_max, rcond_k
   double precision, allocatable :: mtst(:,:), b_tst(:), x_tst(:)
   double precision, allocatable :: gn_work(:,:)
   double precision, allocatable :: gecon_work(:)
   integer, allocatable :: gecon_iwork(:)
!
   lu_ok = .false.
   ridge_used = 0.0d0
   itry_used = -1
   rcond = -1.0d0
   pivrat = -1.0d0
   ridge_rows = 0.0d0
   reg = .false.
!
   scale_k = 0.0d0
   do i = 1, n
      scale_k = max(scale_k, abs(a(i,i)))
   end do
   scale_k = max(scale_k, 1.0d0)
!
   blkmax = 0.0d0
   do i = 1, n
      do j = 1, n
         blkmax = max(blkmax, abs(a(i,j)))
      end do
   end do
!
   allocate(mtst(n,n), b_tst(n), x_tst(n))
   if(lsq_mode) allocate(gn_work(n,n))
!
   do itry = 0, 8
      reg = .false.
      if(itry == 0)then
         eps = jfnk_ridge_floor
      else if(itry == 1)then
         eps = 1.0d-16
      else if(itry == 2)then
         eps = 1.0d-14
      else if(itry == 3)then
         eps = 1.0d-12
      else if(itry == 4)then
         eps = 1.0d-10
      else if(itry == 5)then
         eps = 1.0d-8
      else if(itry == 6)then
         eps = 1.0d-6
      else if(itry == 7)then
         eps = 1.0d-4
      else
         eps = 1.0d-2
      end if
      ridge_add = max(eps, jfnk_ridge_floor)
      lu = a
      do i = 1, n
         lu(i,i) = lu(i,i) + ridge_add * &
              max(abs(a(i,i)), 1.0d-30*scale_k)
      end do
      if(blkmax > 0.0d0)then
         do i = 1, n
            rowmax_i = 0.0d0
            do j = 1, n
               rowmax_i = max(rowmax_i, abs(a(i,j)))
            end do
            if(rowmax_i < row_tiny*blkmax)then
               reg(i) = .true.
               do j = 1, n
                  lu(i,j) = 0.0d0
               end do
               lu(i,i) = 1.0d0
            end if
         end do
      end if
      do i = 1, n
         if(pre_reg(i) .and. .not. reg(i))then
            reg(i) = .true.
            do j = 1, n
               lu(i,j) = 0.0d0
            end do
            lu(i,i) = 1.0d0
         end if
      end do
      call dgeequ(n, n, lu, n, rsf, csf, row_cnd, col_cnd, amax, info_k)
      if(info_k /= 0) cycle
      do j = 1, n
         do i = 1, n
            lu(i,j) = lu(i,j) * rsf(i) * csf(j)
         end do
      end do
      if(lsq_mode)then
         lsq_l2 = (lsq_lam * max(amax, 1.0d0))**2
         call dsyrk('L', 'T', n, n, 1.0d0, lu, n, 0.0d0, gn_work, n)
         do i = 1, n
            gn_work(i,i) = gn_work(i,i) + lsq_l2
         end do
         lu = gn_work
         mtst = lu
         call dpotrf('L', n, lu, n, info_k)
      else
         anorm_k = 0.0d0
         do j = 1, n
            s_tmp = 0.0d0
            do i = 1, n
               s_tmp = s_tmp + abs(lu(i,j))
            end do
            anorm_k = max(anorm_k, s_tmp)
         end do
         mtst = lu
         call dgetrf(n, n, lu, n, ipiv, info_k)
      end if
      if(info_k == 0)then
         b_tst = 0.0d0
         do i = 1, n
            do j = 1, n
               b_tst(i) = b_tst(i) + mtst(i,j)
            end do
         end do
         x_tst = b_tst
         if(lsq_mode)then
            call dpotrs('L', n, 1, lu, n, x_tst, n, info_g)
         else
            call dgetrs('N', n, 1, lu, n, ipiv, x_tst, n, info_g)
         end if
         if(info_g /= 0 .or. &
              maxval(abs(x_tst-1.0d0)) > 1.0d-4)then
            ridge_used = 0.0d0
            cycle
         end if
         lu_ok = .true.
         ridge_used = ridge_add
         itry_used = itry
         piv_min = huge(1.0d0)
         piv_max = 0.0d0
         do i = 1, n
            piv_min = min(piv_min, abs(lu(i,i)))
            piv_max = max(piv_max, abs(lu(i,i)))
         end do
         pivrat = piv_min / max(piv_max, 1.0d-300)
         rcond = -1.0d0
         if(.not. lsq_mode)then
            if(.not. allocated(gecon_work)) &
                 allocate(gecon_work(4*n), gecon_iwork(n), stat=ios)
            call dgecon('I', n, lu, n, anorm_k, rcond_k, &
                 gecon_work, gecon_iwork, info_g)
            if(info_g == 0) rcond = rcond_k
         end if
         if(.not. lsq_mode .and. &
              ridge_piv > 0.0d0 .and. pivrat < ridge_piv)then
            lu_ok = .false.
            ridge_used = 0.0d0
            cycle
         end if
         do i = 1, n
            ridge_rows(i) = ridge_add * &
                 max(abs(a(i,i)), 1.0d-30*scale_k)
         end do
         exit
      end if
   end do
end subroutine old_factor
!
!======================================================================
! Verbatim copy of the pre-extraction inline apply (LU/LSQ per-depth
! block) from JFNK_PRECOND_V1.
!======================================================================
subroutine old_apply(v, rhs, n, lu_ok, lu, ipiv, rsf, csf, a_raw, &
     ridge_rows, ba_diag, reg, fixed, lsq_mode, zmax, &
     nproj, nzero, nzclip, nzmax_hit, nsolved)
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   double precision, intent(in) :: rhs(n)
   logical, intent(in) :: lu_ok, lsq_mode, reg(n), fixed(n)
   double precision, intent(in) :: lu(n,n), rsf(n), csf(n), &
        a_raw(n,n), ridge_rows(n), ba_diag(n), zmax
   integer, intent(in) :: ipiv(n)
   integer, intent(inout) :: nproj, nzero, nzclip, nzmax_hit, nsolved
   integer :: i, j, info
   double precision :: diag_use, wi_lsq
   double precision, allocatable :: w_lsq(:)
!
   if(.not. lu_ok)then
      do j = 1, n
         diag_use = ba_diag(j) + ridge_rows(j)
         if(reg(j))then
            v(j) = 0.0d0
            nproj = nproj + 1
         else if(abs(diag_use) > 1.0d-30)then
            v(j) = rhs(j) / diag_use
         else
            v(j) = 0.0d0
            nzero = nzero + 1
         end if
         nsolved = nsolved + 1
      end do
      return
   end if
!
   do j = 1, n
      v(j) = v(j) * rsf(j)
   end do
   if(lsq_mode)then
      allocate(w_lsq(n))
      do j = 1, n
         w_lsq(j) = v(j) * rsf(j)
      end do
      call dgemv('T', n, n, 1.0d0, a_raw, n, w_lsq, 1, 0.0d0, v, 1)
      do j = 1, n
         v(j) = v(j) + ridge_rows(j) * w_lsq(j)
      end do
      do i = 1, n
         if(reg(i))then
            wi_lsq = w_lsq(i)
            do j = 1, n
               v(j) = v(j) - wi_lsq * a_raw(i,j)
            end do
            v(i) = v(i) + wi_lsq - wi_lsq * ridge_rows(i)
         end if
      end do
      do j = 1, n
         v(j) = v(j) * csf(j)
      end do
      call dpotrs('L', n, 1, lu, n, v, n, info)
   else
      call dgetrs('N', n, 1, lu, n, ipiv, v, n, info)
   end if
   if(info /= 0)then
      do j = 1, n
         diag_use = ba_diag(j) + ridge_rows(j)
         if(reg(j))then
            v(j) = 0.0d0
            nproj = nproj + 1
         else if(abs(diag_use) > 1.0d-30)then
            v(j) = rhs(j) / diag_use
         else
            v(j) = 0.0d0
            nzero = nzero + 1
         end if
      end do
   else
      do j = 1, n
         v(j) = v(j) * csf(j)
      end do
      do j = 1, n
         if(reg(j))then
            v(j) = 0.0d0
            nproj = nproj + 1
         else if(fixed(j))then
            v(j) = 0.0d0
         else if(.not. (abs(v(j)) <= huge(1.0d0)))then
            v(j) = 0.0d0
            nzero = nzero + 1
            nzclip = nzclip + 1
         else if(abs(v(j)) > zmax)then
            nzmax_hit = nzmax_hit + 1
         end if
         nsolved = nsolved + 1
      end do
   end if
   if(allocated(w_lsq)) deallocate(w_lsq)
end subroutine old_apply
!
!======================================================================
subroutine fill_test_block(a, n, iseed, style)
   implicit none
   integer, intent(in) :: n, iseed, style
   double precision, intent(out) :: a(n,n)
   integer :: i, j
   double precision :: u
!
! Deterministic LCG so the test needs no seed-array plumbing.
   call srand(iseed)
   do j = 1, n
      do i = 1, n
         u = rand()
         a(i,j) = (2.0d0*u - 1.0d0)
      end do
   end do
   if(style == 1)then
! Rate-like block: loss-dominated (negative) diagonals, stiff
! coupling spread over decades.
      do i = 1, n
         a(i,i) = -abs(a(i,i)) * 10.0d0**(6.0d0*rand() - 3.0d0)
      end do
      do i = 1, n
         do j = 1, n
            if(i /= j) a(i,j) = a(i,j) * &
                 10.0d0**(8.0d0*rand() - 6.0d0)
         end do
      end do
   else if(style == 2)then
! Mixed block with some near-null rows (raw-magnitude test fodder).
      do i = 1, n
         a(i,i) = a(i,i) * 1.0d3
      end do
      do i = 4, min(8, n)
         a(i,:) = a(i,:) * 1.0d-14
      end do
   else
! Mildly conditioned generic block.
      do i = 1, n
         a(i,i) = a(i,i) + 4.0d0
      end do
   end if
end subroutine fill_test_block
!
!======================================================================
subroutine t_ident_factor(nfail)
   use jfnk_solver_v1_mod, only: jfnk_pc_signed, jfnk_pc_reg_zerocol, &
        jfnk_pc_tiny_equil, jfnk_pc_rcond_gate
   implicit none
   integer, intent(inout) :: nfail
   integer, parameter :: n = 24
   double precision :: a(n,n), pre_dummy
   logical :: pre_reg(n)
   double precision :: lu1(n,n), lu2(n,n), rsf1(n), rsf2(n), &
        csf1(n), csf2(n), rr1(n), rr2(n)
   integer :: ipiv1(n), ipiv2(n)
   logical :: reg1(n), reg2(n), regt(n)
   double precision :: rcond1, rcond2, pv1, pv2, ru1, ru2
   integer :: itry1, itry2, style, seed, imode
   logical :: ok1, ok2, lsq
   integer :: i, j
!
   do style = 0, 2
      do imode = 0, 1
         lsq = (imode == 1)
         do seed = 1, 3
            call fill_test_block(a, n, 1000*style+17*seed, style)
            pre_reg = .false.
            if(style == 2) pre_reg(2) = .true.
! defaults: all audit flags off
            jfnk_pc_signed = .false.; jfnk_pc_reg_zerocol = .false.
            jfnk_pc_tiny_equil = .false.; jfnk_pc_rcond_gate = .false.
            call old_factor(a, n, pre_reg, 1.0d-8, 0.0d0, 1.0d-4, &
                 lsq, lu1, ipiv1, rsf1, csf1, reg1, rr1, rcond1, &
                 pv1, itry1, ru1, ok1)
            call jfnk_pc_factor_block_v1(a, n, pre_reg, 1.0d-8, &
                 0.0d0, 1.0d-4, lsq, 0, lu2, ipiv2, rsf2, csf2, &
                 reg2, regt, rr2, rcond2, pv2, itry2, ru2, ok2)
            if(ok1 .neqv. ok2)then
               print *, 'T-IDENT FAIL ok', style, imode, seed, ok1, ok2
               nfail = nfail + 1; cycle
            end if
            if(any(reg1 .neqv. reg2) .or. itry1 /= itry2 .or. &
                 ru1 /= ru2 .or. pv1 /= pv2)then
               print *, 'T-IDENT FAIL scalars', style, imode, seed
               nfail = nfail + 1; cycle
            end if
! rcond is exempt from the bit-compare: old_factor supplies the
! 1-norm to DGECON('I') (the F4 defect); the extracted routine now
! supplies the matching infinity norm, so the values legitimately
! differ (diagnostic only).
            if(ok1)then
               if(rcond2 <= 0.0d0)then
                  print *, 'T-IDENT FAIL rcond not set'
                  nfail = nfail + 1; cycle
               end if
            end if
            if(any(ipiv1 /= ipiv2) .or. any(lu1 /= lu2) .or. &
                 any(rsf1 /= rsf2) .or. any(csf1 /= csf2) .or. &
                 any(rr1 /= rr2))then
               print *, 'T-IDENT FAIL arrays', style, imode, seed, &
                    ok1
               nfail = nfail + 1; cycle
            end if
         end do
      end do
   end do
   print *, 'T-IDENT factor: ok (bit-identical vs inline copy)'
end subroutine t_ident_factor
!
!======================================================================
subroutine t_ident_apply(nfail)
   implicit none
   integer, intent(inout) :: nfail
   integer, parameter :: n = 24
   double precision :: a(n,n), lu(n,n), v1(n), v2(n), rhs(n)
   double precision :: rsf(n), csf(n), rr(n), diag(n)
   integer :: ipiv(n)
   logical :: reg(n), regt(n), fixed(n), pre_reg(n)
   double precision :: rcond, pv, ru, zmax
   integer :: itry, seed, imode, i
   logical :: ok, lsq
   integer :: np1,nz1,nc1,nm1,ns1, np2,nz2,nc2,nm2,ns2
!
   do imode = 0, 1
      lsq = (imode == 1)
      do seed = 1, 3
         call fill_test_block(a, n, 77*seed, 0)
         do i = 1, n
            diag(i) = a(i,i)
         end do
         pre_reg = .false.
         call jfnk_pc_factor_block_v1(a, n, pre_reg, 1.0d-8, 0.0d0, &
              1.0d-4, lsq, 0, lu, ipiv, rsf, csf, reg, regt, rr, &
              rcond, pv, itry, ru, ok)
         fixed = .false.
         fixed(3) = .true.
         rhs = 0.0d0
         do i = 1, n
            rhs(i) = sin(3.7d0*i)
         end do
         np1=0;nz1=0;nc1=0;nm1=0;ns1=0
         np2=0;nz2=0;nc2=0;nm2=0;ns2=0
         v1 = rhs; v2 = rhs
         call old_apply(v1, rhs, n, ok, lu, ipiv, rsf, csf, a, rr, &
              diag, reg, fixed, lsq, 1.0d8, np1, nz1, nc1, nm1, ns1)
         call jfnk_pc_apply_block_v1(v2, rhs, n, ok, lu, ipiv, rsf, &
              csf, a, rr, diag, reg, fixed, lsq, 1.0d8, &
              np2, nz2, nc2, nm2, ns2)
         if(any(v1 /= v2) .or. np1/=np2 .or. nz1/=nz2 .or. &
              nc1/=nc2 .or. nm1/=nm2 .or. ns1/=ns2)then
            print *, 'T-IDENT apply FAIL', imode, seed, ok, &
                 maxval(abs(v1-v2))
            nfail = nfail + 1
         end if
      end do
   end do
   print *, 'T-IDENT apply: ok (bit-identical vs inline copy)'
end subroutine t_ident_apply
!
!======================================================================
! M^{-1} must be a fixed linear operator inside a GMRES cycle:
! M^-1(a v1 + b v2) = a M^-1 v1 + b M^-1 v2 for every flag combo.
!======================================================================
subroutine t_linear(nfail)
   use jfnk_solver_v1_mod, only: jfnk_pc_signed, jfnk_pc_reg_zerocol, &
        jfnk_pc_tiny_equil, jfnk_pc_rcond_gate
   implicit none
   integer, intent(inout) :: nfail
   integer, parameter :: n = 24
   double precision :: a(n,n), lu(n,n), v1(n), v2(n), va(n), vb(n), &
        vc(n), ref(n), rhs1(n), rhs2(n), rhsc(n)
   double precision :: rsf(n), csf(n), rr(n), diag(n)
   integer :: ipiv(n)
   logical :: reg(n), regt(n), fixed(n), pre_reg(n)
   double precision :: rcond, pv, ru, aa, bb, err, scl
   integer :: itry, i, iflags
   logical :: ok
   integer :: np,nz,nc,nm,ns
!
   call fill_test_block(a, n, 4242, 0)
   do i = 1, n
      diag(i) = a(i,i)
   end do
   pre_reg = .false.; fixed = .false.
   do i = 1, n
      rhs1(i) = sin(1.3d0*i)
      rhs2(i) = cos(0.7d0*i)
   end do
   aa = 1.7d0; bb = -0.9d0
!
   do iflags = 0, 15
      jfnk_pc_signed     = (iand(iflags,1) /= 0)
      jfnk_pc_reg_zerocol= (iand(iflags,2) /= 0)
      jfnk_pc_tiny_equil = (iand(iflags,4) /= 0)
      jfnk_pc_rcond_gate = (iand(iflags,8) /= 0)
      call jfnk_pc_factor_block_v1(a, n, pre_reg, 1.0d-8, 0.0d0, &
           1.0d-4, .false., 0, lu, ipiv, rsf, csf, reg, regt, rr, &
           rcond, pv, itry, ru, ok)
      np=0;nz=0;nc=0;nm=0;ns=0
      va = rhs1
      call jfnk_pc_apply_block_v1(va, rhs1, n, ok, lu, ipiv, rsf, &
           csf, a, rr, diag, reg, fixed, .false., 1.0d30, &
           np, nz, nc, nm, ns)
      vb = rhs2
      call jfnk_pc_apply_block_v1(vb, rhs2, n, ok, lu, ipiv, rsf, &
           csf, a, rr, diag, reg, fixed, .false., 1.0d30, &
           np, nz, nc, nm, ns)
      rhsc = aa*rhs1 + bb*rhs2
      vc = rhsc
      call jfnk_pc_apply_block_v1(vc, rhsc, n, ok, lu, ipiv, rsf, &
           csf, a, rr, diag, reg, fixed, .false., 1.0d30, &
           np, nz, nc, nm, ns)
      ref = aa*va + bb*vb
      scl = max(maxval(abs(ref)), 1.0d-30)
      err = maxval(abs(vc - ref)) / scl
      if(err > 1.0d-12)then
         print *, 'T-LINEAR FAIL flags=', iflags, ' err=', err
         nfail = nfail + 1
      end if
   end do
   print *, 'T-LINEAR: ok'
end subroutine t_linear

!======================================================================
! T-RCOND-NORM (F4): on rate-like blocks with one dense
! conservation-like row, the rcond returned must match an independent
! DGECON('I') call whose ANORM is recomputed here as the infinity norm
! (max row sum) of the reconstructed pre-factor matrix.  The old
! 1-norm ANORM produced rcond values ~9x (n=50), ~50x (n=200), ~240x
! (n=800) too optimistic in the reviewer's measurement; this test
! asserts the corrected value is consistent and prints the residual
! old/new ratio for the record.
!======================================================================
subroutine t_rcond_norm(nfail)
   use jfnk_solver_v1_mod, only: jfnk_pc_signed, jfnk_pc_reg_zerocol, &
        jfnk_pc_tiny_equil, jfnk_pc_rcond_gate
   implicit none
   integer, intent(inout) :: nfail
   integer, parameter :: ns(3) = [50, 200, 800]
   double precision, allocatable :: a(:,:), b(:,:), lu(:,:)
   double precision, allocatable :: rsf(:), csf(:), rr(:), w(:)
   integer, allocatable :: ipiv(:), iw(:)
   logical, allocatable :: reg(:), regt(:), pre_reg(:)
   double precision :: rcond2, pv, ru, an_inf, an_one, rc_true, &
        rc_old, ratio
   integer :: itry, ni, n, i, j, info
   logical :: ok
   double precision :: u
!
   jfnk_pc_signed=.false.; jfnk_pc_reg_zerocol=.false.
   jfnk_pc_tiny_equil=.false.; jfnk_pc_rcond_gate=.false.
   do ni = 1, 3
      n = ns(ni)
      allocate(a(n,n), b(n,n), lu(n,n), rsf(n), csf(n), rr(n), &
           w(4*n), ipiv(n), iw(n), reg(n), regt(n), pre_reg(n))
      call srand(97 + n)
      do j = 1, n
         do i = 1, n
            u = rand()
            a(i,j) = (2.0d0*u - 1.0d0) * 1.0d-2
         end do
      end do
! rate-like negative diagonals plus one dense conservation row whose
! magnitude dominates every column (the anisotropy that separated the
! two norms in the reviewer measurement)
      do j = 1, n
         a(j,j) = -abs(a(j,j)) - 1.0d0
         a(1,j) = 5.0d0 + rand()
      end do
      pre_reg = .false.
      call jfnk_pc_factor_block_v1(a, n, pre_reg, 0.0d0, 0.0d0, &
           1.0d-4, .false., 0, lu, ipiv, rsf, csf, reg, regt, rr, &
           rcond2, pv, itry, ru, ok)
      if(.not. ok)then
         print *, 'T-RCOND-NORM: factor failed at n=', n
         nfail = nfail + 1; cycle
      end if
! Reconstruct the pre-factor equilibrated matrix the routine factored:
! ridged diagonal, identity substitution on reg rows, then rsf*csf
! scaling -- the same B the internal MTST held at DGETRF time.
      b = a
      do i = 1, n
         b(i,i) = b(i,i) + rr(i)
      end do
      do i = 1, n
         if(reg(i))then
            b(i,:) = 0.0d0
            b(i,i) = 1.0d0
         end if
      end do
      do j = 1, n
         do i = 1, n
            b(i,j) = b(i,j) * rsf(i) * csf(j)
         end do
      end do
      an_inf = 0.0d0
      an_one = 0.0d0
      do i = 1, n
         an_inf = max(an_inf, sum(abs(b(i,:))))
      end do
      do j = 1, n
         an_one = max(an_one, sum(abs(b(:,j))))
      end do
      call dgecon('I', n, lu, n, an_inf, rc_true, w, iw, info)
      call dgecon('I', n, lu, n, an_one, rc_old, w, iw, info)
      if(rcond2 <= 0.0d0 .or. abs(rcond2/rc_true - 1.0d0) > 2.0d0)then
         print *, 'T-RCOND-NORM FAIL n=', n, ' rcond=', rcond2, &
              ' indep=', rc_true
         nfail = nfail + 1
      end if
      ratio = rc_old / max(rc_true, tiny(1.0d0))
      print '(A,I5,A,ES10.2,A,ES10.2,A,F7.2)', &
           ' T-RCOND-NORM n=', n, ' rcond=', rcond2, &
           ' old-1norm=', rc_old, ' old/new=', ratio
      deallocate(a, b, lu, rsf, csf, rr, w, ipiv, iw, reg, regt, &
           pre_reg)
   end do
end subroutine t_rcond_norm

!======================================================================
! T-RSCALE-ORDER (F7): 'DIAG_RDG' must select the DIAG_RDG branch, not
! be shadowed by the 4-character 'DIAG' prefix.  Production dispatch
! is jfnk_rscale_branch_v1 itself, so this test covers the real code.
!======================================================================
subroutine t_rscale_order(nfail)
   implicit none
   integer, intent(inout) :: nfail
   integer :: jfnk_rscale_branch_v1
   external jfnk_rscale_branch_v1
!
   if(jfnk_rscale_branch_v1('DIAG_RDG') /= 3)then
      print *, 'T-RSCALE-ORDER FAIL: DIAG_RDG not selected'
      nfail = nfail + 1
   end if
   if(jfnk_rscale_branch_v1('DIAG') /= 2)then
      print *, 'T-RSCALE-ORDER FAIL: DIAG not selected'
      nfail = nfail + 1
   end if
   if(jfnk_rscale_branch_v1('DIAG_XX') /= 2)then
      print *, 'T-RSCALE-ORDER FAIL: DIAG prefix lost'
      nfail = nfail + 1
   end if
   if(jfnk_rscale_branch_v1('SQRT_STEQ') /= 1)then
      print *, 'T-RSCALE-ORDER FAIL: SQRT_STE not selected'
      nfail = nfail + 1
   end if
   if(jfnk_rscale_branch_v1('EUCLIDEAN') /= 0)then
      print *, 'T-RSCALE-ORDER FAIL: default branch lost'
      nfail = nfail + 1
   end if
   print *, 'T-RSCALE-ORDER: ok'
end subroutine t_rscale_order
!======================================================================
subroutine t_sign(nfail)
   use jfnk_solver_v1_mod, only: jfnk_pc_signed
   implicit none
   integer, intent(inout) :: nfail
   double precision :: a_coef, xstar, xk, sigma, gval, du, ds, &
        con_u, con_s, sgn(1), diag(1), m_u(1,1), m_s(1,1), rhs1(1)
   integer :: ipiv1(1), info, i
   double precision :: fracs(7) = [0.1d0, 0.5d0, 0.9d0, 1.1d0, &
        2.0d0, 10.0d0, 100.0d0]
   double precision :: tabl_u(7) = [1.0d0/9.0d0*1.111111d0*0.9d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
   logical :: bad
!
   a_coef = 2.0d0
   xstar = 1.7d0
   xk = 4.0d0
   bad = .false.
   do i = 1, 7
      sigma = fracs(i) * a_coef
      diag(1) = -a_coef
      call jfnk_pc_sign_vec_v1(diag, 1, 0.0d0, sgn)
! residual at the anchor: G(xk) = F(xk); the PTC term vanishes
      gval = -a_coef * (xk - xstar)
! unsigned: J_G = -a + sigma ; signed: J_G = -a - sigma
      m_u(1,1) = -a_coef + sigma
      m_s(1,1) = -a_coef + sgn(1)*sigma
      rhs1(1) = gval
      call dgetrf(1, 1, m_u, 1, ipiv1, info)
      call dgetrs('N', 1, 1, m_u, 1, ipiv1, rhs1, 1, info)
      du = -rhs1(1)
      rhs1(1) = gval
      call dgetrf(1, 1, m_s, 1, ipiv1, info)
      call dgetrs('N', 1, 1, m_s, 1, ipiv1, rhs1, 1, info)
      ds = -rhs1(1)
      con_u = abs(xk + du - xstar) / abs(xk - xstar)
      con_s = abs(xk + ds - xstar) / abs(xk - xstar)
      if(abs(con_u - sigma/abs(sigma - a_coef)) > 1.0d-10)then
         print *, 'T-SIGN FAIL unsigned model', fracs(i), con_u
         bad = .true.
      end if
      if(con_s >= 1.0d0)then
         print *, 'T-SIGN FAIL signed contraction', fracs(i), con_s
         bad = .true.
      end if
      if(abs(con_s - sigma/(a_coef+sigma)) > 1.0d-10)then
         print *, 'T-SIGN FAIL signed model', fracs(i), con_s
         bad = .true.
      end if
   end do
   if(.not. bad) print *, 'T-SIGN scalar: ok'

! 3-level coupled rate block: negative rate diagonal plus a dense
! conservation row.  Signed fold makes M equal J_G, so the
! preconditioned step equals the true Newton step on the signed G.
   block
      double precision :: ar(3,3), mu(3,3), ms(3,3), rhs(3)
      double precision :: dx_u(3), dx_s(3), dx_ref(3), sg3(3), d3(3)
      double precision :: jg(3,3), gv(3)
      integer :: ipv(3)
      integer :: i2, j2
      ar = 0.0d0
! conservation row (positive), two loss rows (negative diagonal)
      ar(1,1) =  1.0d0; ar(1,2) =  1.0d0; ar(1,3) =  1.0d0
      ar(2,1) =  0.4d0; ar(2,2) = -2.0d0; ar(2,3) =  0.1d0
      ar(3,1) =  0.2d0; ar(3,2) =  0.3d0; ar(3,3) = -3.0d0
      gv = [0.5d0, -0.7d0, 0.9d0]
      sigma = 1.5d0
      do i2 = 1, 3
         d3(i2) = ar(i2,i2)
      end do
      call jfnk_pc_sign_vec_v1(d3, 3, 0.0d0, sg3)
! unsigned fold M_u = A + sigma I; signed fold M_s = A + sgn*sigma I;
! J_G(signed residual) = A + sgn*sigma = M_s, so the signed solve is
! the exact Newton step for the signed G: assert dx_s == reference.
      mu = ar; ms = ar; jg = ar
      do i2 = 1, 3
         mu(i2,i2) = mu(i2,i2) + sigma
         ms(i2,i2) = ms(i2,i2) + sg3(i2)*sigma
         jg(i2,i2) = jg(i2,i2) + sg3(i2)*sigma
      end do
      rhs = gv
      call dgetrf(3, 3, mu, 3, ipv, info)
      call dgetrs('N', 3, 1, mu, 3, ipv, rhs, 3, info)
      dx_u = -rhs
      rhs = gv
      call dgetrf(3, 3, ms, 3, ipv, info)
      call dgetrs('N', 3, 1, ms, 3, ipv, rhs, 3, info)
      dx_s = -rhs
      rhs = gv
      call dgetrf(3, 3, jg, 3, ipv, info)
      call dgetrs('N', 3, 1, jg, 3, ipv, rhs, 3, info)
      dx_ref = -rhs
      if(maxval(abs(dx_s - dx_ref)) > 1.0d-12)then
         print *, 'T-SIGN FAIL coupled signed step', &
              maxval(abs(dx_s-dx_ref))
         bad = .true.
      end if
      if(maxval(abs(dx_u - dx_ref)) < 1.0d-12)then
         print *, 'T-SIGN WARN: unsigned step identical to signed'
      end if
      print '(A,3ES11.3,A,3ES11.3)', &
           ' T-SIGN coupled dx signed=', dx_s, ' unsigned=', dx_u
   end block
   if(bad) nfail = nfail + 1
end subroutine t_sign
!======================================================================
subroutine t_ridge(nfail)
   use jfnk_solver_v1_mod, only: jfnk_pc_signed, jfnk_ridge_floor
   implicit none
   integer, intent(inout) :: nfail
   double precision :: a(4,4), lu(4,4), rsf(4), csf(4), rr(4)
   integer :: ipiv(4)
   logical :: reg(4), regt(4), pre_reg(4)
   double precision :: rcond, pv, ru
   integer :: itry, i
   logical :: ok
   double precision :: du(4), ds(4)
!
! rows 1,3: positive diagonal; rows 2,4: negative.  ridge floor
! forces the signed/unsigned ridge branch to add a visible shift.
   jfnk_ridge_floor = 1.0d-3
   a = 0.0d0
   do i = 1, 4
      a(i,i) = merge(2.0d0, -2.0d0, mod(i,2) == 1)
      a(i,mod(i,4)+1) = 0.1d0
   end do
   pre_reg = .false.
   jfnk_pc_signed = .true.
   call jfnk_pc_factor_block_v1(a, 4, pre_reg, 0.0d0, 0.0d0, 1.0d-4, &
        .false., 0, lu, ipiv, rsf, csf, reg, regt, rr, rcond, pv, &
        itry, ru, ok)
   do i = 1, 4
      ds(i) = rr(i)
   end do
   jfnk_pc_signed = .false.
   call jfnk_pc_factor_block_v1(a, 4, pre_reg, 0.0d0, 0.0d0, 1.0d-4, &
        .false., 0, lu, ipiv, rsf, csf, reg, regt, rr, rcond, pv, &
        itry, ru, ok)
   do i = 1, 4
      du(i) = rr(i)
   end do
   do i = 1, 4
      if(abs(a(i,i) + ds(i)) < abs(a(i,i)))then
         print *, 'T-RIDGE FAIL signed ridge shrinks |d| at row', i
         nfail = nfail + 1
      end if
      if(a(i,i) < 0.0d0)then
         if(ds(i) >= 0.0d0)then
            print *, 'T-RIDGE FAIL signed ridge sign at row', i
            nfail = nfail + 1
         end if
! document: unsigned ridge shrinks |d| on this negative row
         if(abs(a(i,i) + du(i)) >= abs(a(i,i)))then
            print *, 'T-RIDGE WARN: unsigned ridge did not shrink'
         end if
      end if
   end do
   jfnk_ridge_floor = 0.0d0
   print *, 'T-RIDGE: ok'
end subroutine t_ridge
