!----------------------------------------------------------------------
! TEST_PC_QUALITY: numerical quality of the per-depth preconditioner
! block factorization/apply (jfnk_pc_factor_block_v1 /
! jfnk_pc_apply_block_v1), linked from the production object.
!
! The audit question (JT2): is M^{-1} a good inverse of the block J?
! The production acceptance tests check M's own invertibility (x=1
! micro-test or rcond(M)) -- never ||I - M^{-1} A||, which is what
! GMRES actually sees. These tests measure it over all canonical
! directions.
!
!   Q1  well-conditioned block -> lu_ok, small ||I - M^{-1}A||.
!   Q2  acceptance blindspot: A = I - (1-d)*u u^T, u .|. 1. The x=1
!       micro-test passes exactly while lambda_u(A) = d = 1e-9 makes
!       M^{-1}A annihilate u (~1e-7). Reports what the rcond gate sees
!       (rcond of M, not of M^{-1}A) -- neither gate catches it.
!   Q3  exact-zero row -> reg substitution; apply projects z_i = 0.
!   Q4  LSQ mode quality on the same well-conditioned block.
!   Q5  row_tiny raw vs equilibrated classification on a unit-scaled
!       full-rank row ([JFNK_ROW_TINY_EQUIL]).
!   Q6  ridge_piv gate rejects a successful-but-ill-conditioned LU.
!
! edited 09-Nov-2026
!----------------------------------------------------------------------
program test_pc_quality
   use jfnk_solver_v1_mod, only: jfnk_pc_set_options_v1, &
      jfnk_pc_tiny_equil, jfnk_rcond_gate, jfnk_rcond_min, &
      jfnk_itry_cap, jfnk_ridge_floor
   implicit none
   integer, parameter :: N = 12, LUER = 6
   double precision :: a(N,N), lu(N,N), rsf(N), csf(N), rrg(N)
   double precision :: ba_diag(N), ridge_rows(N), rcond, pivrat
   double precision :: ridge_used, v(N), rhs(N), z(N,N), err_m, lam
   integer :: ipiv(N), itry_used, nproj, nzero, nzclip, nzmax_hit
   integer :: nsolved, i, j, info, nfail
   logical :: pre_reg(N), reg(N), reg_tiny(N), fixed(N), lu_ok
   logical :: lsq_mode
   double precision :: u(N), av(N)
   !
   nfail = 0
   fixed = .false.
   pre_reg = .false.
   lsq_mode = .false.
   !
   ! restore production defaults for every case
   call jfnk_pc_set_options_v1(.false., .false., .false., .false., &
                             1.0d-13, 3)
   !
   ! ================= Q1: well-conditioned block =================
   !
   call gen_wellcond(a, N)
   call factor_and_apply_err(a, N, lsq_mode, lu_ok, itry_used, &
                             ridge_used, rcond, err_m)
   write(luer,'(A,L3,A,I3,A,ES10.2,A,ES10.2)') &
      'Q1: lu_ok=', lu_ok, ' itry=', itry_used, &
      ' ridge=', ridge_used, '  ||I-MinvA||inf=', err_m
   if(lu_ok .and. itry_used == 0 .and. err_m < 0.15d0) then
      write(luer,'(A)') 'Q1: PASS'
   else
      write(luer,'(A)') 'Q1: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ======== Q2: acceptance blindspot (u .perp. 1, lam=1e-9) ========
   !
   u = 0d0
   u(1) = 1d0 / sqrt(2d0)
   u(2) = -1d0 / sqrt(2d0)
   lam = 1d-9
   do j = 1, N
      do i = 1, N
         a(i,j) = merge(1d0, 0d0, i == j) - (1d0 - lam) * u(i) * u(j)
      end do
   end do
   ! A*1 = 1 exactly, so the default micro-test sees a perfect solve.
   call factor_and_apply_err(a, N, lsq_mode, lu_ok, itry_used, &
                             ridge_used, rcond, err_m)
   write(luer,'(A,L3,A,I3,A,ES10.2,A,ES10.2)') &
      'Q2: lu_ok=', lu_ok, ' itry=', itry_used, &
      ' rcond(M)=', rcond, '  ||I-MinvA||inf=', err_m
   if(lu_ok .and. err_m > 0.5d0) then
      write(luer,'(A)') 'Q2: BLINDSPOT CONFIRMED -- accepted at ' // &
         'itry=0 while M^{-1}A deviates by >0.5'
   else
      write(luer,'(A)') 'Q2: blindspot not reproduced'
   end if
   ! same matrix under the F5 rcond gate: gate measures rcond(M),
   ! which is ~ridge/(1+ridge) ~ 1e-2 -- still far above 1e-13, so
   ! the gate accepts the same bad preconditioner.
   call jfnk_pc_set_options_v1(.false., .false., .false., .true., &
                             1.0d-13, 3)
   call factor_and_apply_err(a, N, lsq_mode, lu_ok, itry_used, &
                             ridge_used, rcond, err_m)
   write(luer,'(A,L3,A,ES10.2,A,ES10.2)') &
      'Q2b: gate-on lu_ok=', lu_ok, ' rcond(M)=', rcond, &
      ' ||I-MinvA||=', err_m
   call jfnk_pc_set_options_v1(.false., .false., .false., .false., &
                             1.0d-13, 3)
   !
   ! ============ Q3: exact-zero row -> reg projection ============
   !
   call gen_wellcond(a, N)
   a(5,:) = 0d0            ! row 5 identically zero
   call factor_full(a, N, .false., lu_ok, itry_used, ridge_used, &
                    rcond, lu, ipiv, rsf, csf, reg, reg_tiny, &
                    ridge_rows, ba_diag)
   write(luer,'(A,L3,A,L3,A,I3)') 'Q3: lu_ok=', lu_ok, &
      ' reg(5)=', reg(5), ' itry=', itry_used
   if(lu_ok .and. reg(5)) then
      ! apply on e_7: the reg row must return exactly 0 while the
      ! free component z_7 ~ (A^-1)_77 ~ 1/a77 is the honest solve.
      rhs = 0d0; rhs(7) = 1d0
      v = rhs
      call apply_block(v, rhs, N, lu_ok, lu, ipiv, rsf, csf, a, &
                       ridge_rows, ba_diag, reg, fixed, .false.)
      if(v(5) == 0d0 .and. abs(v(7) - 1d0/a(7,7)) < 0.1d0) then
         write(luer,'(A,2ES10.2)') 'Q3: PASS  z_5, z_7=', v(5), v(7)
      else
         write(luer,'(A,2ES12.4)') 'Q3: FAIL  v5,v7=', v(5), v(7)
         nfail = nfail + 1
      end if
   else
      write(luer,'(A)') 'Q3: FAIL (reg row not detected)'
      nfail = nfail + 1
   end if
   !
   ! ============== Q4: LSQ mode quality ==============
   ! Known defect this asserts the fix for: DSYRK('L','T') fills only
   ! the lower triangle of gn_work, so mtst's upper triangle is
   ! garbage and the b=mtst*1 micro-test can never pass -> lsq_mode
   ! always returns lu_ok=F and production silently falls back to
   ! Jacobi when JFNK_PRECOND=LSQ.
   call gen_wellcond(a, N)
   lsq_mode = .true.
   call factor_and_apply_err(a, N, lsq_mode, lu_ok, itry_used, &
                             ridge_used, rcond, err_m)
   write(luer,'(A,L3,A,I3,A,ES10.2)') &
      'Q4: lsq lu_ok=', lu_ok, ' itry=', itry_used, &
      '  ||I-MinvA||inf=', err_m
   if(lu_ok .and. err_m < 0.2d0) then
      write(luer,'(A)') 'Q4: PASS'
   else
      write(luer,'(A)') 'Q4: FAIL -- lsq factorization never ' // &
         'accepted (mtst upper triangle of gn_work uninitialized)'
      nfail = nfail + 1
   end if
   lsq_mode = .false.
   !
   ! ===== Q5: raw vs equilibrated near-null row classification =====
   !
   a = 0d0
   do i = 1, N
      a(i,i) = 1d0
   end do
   a(1,1) = 1d-9         ! row 1: only entry 1e-9, owns its column
   ! raw test:  rowmax=1e-9 < row_tiny*blkmax = 1e-4 -> tiny
   ! equil test: share = 1e-9/csf(1)=1e-9/1e-9 = 1 -> NOT tiny
   call jfnk_pc_set_options_v1(.false., .false., .false., .false., &
                             1.0d-13, 3)
   call factor_reg_only(a, N, reg_tiny)
   write(luer,'(A,L3)') 'Q5: raw  reg_tiny(1)=', reg_tiny(1)
   if(.not. reg_tiny(1)) then
      write(luer,'(A)') 'Q5: FAIL (raw test should flag row 1)'
      nfail = nfail + 1
   end if
   call jfnk_pc_set_options_v1(.false., .false., .true., .false., &
                             1.0d-13, 3)
   call factor_reg_only(a, N, reg_tiny)
   write(luer,'(A,L3)') 'Q5: equil reg_tiny(1)=', reg_tiny(1)
   if(reg_tiny(1)) then
      write(luer,'(A)') 'Q5: FAIL (equil test should clear row 1)'
      nfail = nfail + 1
   else
      write(luer,'(A)') 'Q5: PASS'
   end if
   call jfnk_pc_set_options_v1(.false., .false., .false., .false., &
                             1.0d-13, 3)
   !
   ! ===== Q6: ridge_piv gate reachability =====
   ! pivrat is measured on the POST-equilibration LU diagonals:
   ! DGEEQU already normalized rows/columns, so partial pivoting
   ! returns O(1) pivots and pivrat ~1 on typical blocks -- the gate
   ! is nearly untriggerable by construction. Documented here on a
   ! block whose raw diagonal spans 1e-14..1.
   a = 0d0
   do i = 1, N
      a(i,i) = 1d0
   end do
   a(3,3) = 1d-13
   a(3,4) = 2d-4         ! row max 2e-4 > 1e-4*blkmax: survives reg
   a(4,3) = 1d0          ! pins csf(3)=1 so equil cannot lift the diag
   a(4,4) = 0d0
   call factor_piv(a, N, 0.5d0, lu_ok, itry_used, pivrat)
   write(luer,'(A,L3,A,I3,A,ES10.2)') &
      'Q6: piv-gate lu_ok=', lu_ok, ' itry=', itry_used, &
      ' pivrat=', pivrat
   if(pivrat > 0.9d0) then
      write(luer,'(A)') 'Q6: FINDING -- equil+pivot self-normalize;' // &
         ' ridge_piv never sees the raw ill-conditioning'
   else
      write(luer,'(A)') 'Q6: gate can fire on this matrix'
   end if
   !
   if(nfail == 0) then
      print *, 'ALL PC-QUALITY TESTS PASSED'
   else
      print *, nfail, ' PC-QUALITY TESTS FAILED'
   end if
   call exit(nfail)
contains

   subroutine gen_wellcond(m, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(out) :: m(n,n)
      integer :: i, j
      do j = 1, n
         do i = 1, n
            m(i,j) = 0.05d0 * sin(1.7d0 * i + 0.9d0 * j) / n
         end do
      end do
      do i = 1, n
         m(i,i) = 1d0 + 0.1d0 * i
      end do
   end subroutine gen_wellcond

   ! factor, then measure max over basis e_j of ||apply(A*e_j)-e_j||
   subroutine factor_and_apply_err(m, n, lsq, ok, itry, ru, rc, err)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: m(n,n)
      logical, intent(in) :: lsq
      logical, intent(out) :: ok
      integer, intent(out) :: itry
      double precision, intent(out) :: ru, rc, err
      double precision :: l_lu(n,n), l_rsf(n), l_csf(n), l_rr(n)
      double precision :: l_bd(n), l_pv, l_v(n), l_rhs(n)
      integer :: l_ip(n)
      logical :: l_reg(n), l_rt(n), l_fx(n)
      integer :: j
      call jfnk_pc_factor_block_v1(m, n, pre_reg, 1.0d-4, 0.0d0, &
           1.0d-2, lsq, luer, l_lu, l_ip, l_rsf, l_csf, l_reg, &
           l_rt, l_rr, rc, l_pv, itry, ru, ok)
      if(.not. ok) then
         err = -1d0
         return
      end if
      l_bd = 0d0
      do j = 1, n
         l_bd(j) = m(j,j)
      end do
      l_fx = .false.
      err = 0d0
      do j = 1, n
         l_rhs = m(:,j)
         l_v = l_rhs
         call apply_block(l_v, l_rhs, n, ok, l_lu, l_ip, l_rsf, &
                          l_csf, m, l_rr, l_bd, l_reg, l_fx, lsq)
         l_v(j) = l_v(j) - 1d0
         err = max(err, maxval(abs(l_v)))
      end do
   end subroutine factor_and_apply_err

   subroutine apply_block(vv, rr, n, ok, l_lu, l_ip, l_rsf, l_csf, &
                          m, l_rr, l_bd, l_reg, l_fx, lsq)
      implicit none
      integer, intent(in) :: n
      double precision, intent(inout) :: vv(n)
      double precision, intent(in) :: rr(n), l_lu(n,n), l_rsf(n), &
           l_csf(n), m(n,n), l_rr(n), l_bd(n)
      integer, intent(in) :: l_ip(n)
      logical, intent(in) :: ok, l_reg(n), l_fx(n), lsq
      integer :: np_, nz_, nc_, nm_, ns_
      np_ = 0; nz_ = 0; nc_ = 0; nm_ = 0; ns_ = 0
      call jfnk_pc_apply_block_v1(vv, rr, n, ok, l_lu, l_ip, l_rsf, &
           l_csf, m, l_rr, l_bd, l_reg, l_fx, lsq, 1.0d8, &
           np_, nz_, nc_, nm_, ns_)
   end subroutine apply_block

   ! full factor returning every output needed downstream
   subroutine factor_full(m, n, lsq, ok, itry, ru, rc, l_lu, l_ip, &
                          l_rsf, l_csf, l_reg, l_rt, l_rr, l_bd)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: m(n,n)
      logical, intent(in) :: lsq
      logical, intent(out) :: ok, l_reg(n), l_rt(n)
      integer, intent(out) :: itry, l_ip(n)
      double precision, intent(out) :: ru, rc, l_lu(n,n), &
           l_rsf(n), l_csf(n), l_rr(n), l_bd(n)
      double precision :: l_pv
      integer :: i
      call jfnk_pc_factor_block_v1(m, n, pre_reg, 1.0d-4, 0.0d0, &
           1.0d-2, lsq, luer, l_lu, l_ip, l_rsf, l_csf, l_reg, &
           l_rt, l_rr, rc, l_pv, itry, ru, ok)
      l_bd = 0d0
      l_bd(1:n) = [(m(i,i), i=1,n)]
   end subroutine factor_full

   subroutine factor_reg_only(m, n, l_rt)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: m(n,n)
      logical, intent(out) :: l_rt(n)
      double precision :: l_lu(n,n), l_rsf(n), l_csf(n), l_rr(n)
      double precision :: rc, l_pv, ru
      integer :: l_ip(n), itry
      logical :: ok, l_reg(n)
      call jfnk_pc_factor_block_v1(m, n, pre_reg, 1.0d-4, 0.0d0, &
           1.0d-2, .false., luer, l_lu, l_ip, l_rsf, l_csf, l_reg, &
           l_rt, l_rr, rc, l_pv, itry, ru, ok)
   end subroutine factor_reg_only

   subroutine factor_piv(m, n, rp, ok, itry, pv)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: m(n,n), rp
      logical, intent(out) :: ok
      integer, intent(out) :: itry
      double precision, intent(out) :: pv
      double precision :: l_lu(n,n), l_rsf(n), l_csf(n), l_rr(n)
      double precision :: rc, ru
      integer :: l_ip(n)
      logical :: l_reg(n), l_rt(n)
      call jfnk_pc_factor_block_v1(m, n, pre_reg, 1.0d-4, rp, &
           1.0d-2, .false., luer, l_lu, l_ip, l_rsf, l_csf, l_reg, &
           l_rt, l_rr, rc, pv, itry, ru, ok)
   end subroutine factor_piv

end program test_pc_quality
