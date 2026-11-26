program drv_cmf_blkband_pin
   use set_kind_module
   implicit none
   integer, parameter :: n = 3
   real(kind=ldp) :: c_mat(n,n), rhs(n), saved_mat(n,n), saved_rhs(n)
   integer :: n_pin, failures

   external :: cmf_blkband_pin_null_v1
   external :: cmf_blkband_pin_sym_v1

   failures = 0
   c_mat = 0.0_ldp
   c_mat(1,1) = 2.0_ldp
   c_mat(3,3) = 3.0_ldp
   rhs = [4.0_ldp, 0.0_ldp, 9.0_ldp]
   call cmf_blkband_pin_null_v1(c_mat, rhs, n, n_pin)
   call check(n_pin == 1 .and. all(c_mat == reshape( &
      [2.0_ldp, 0.0_ldp, 0.0_ldp, 0.0_ldp, 1.0_ldp, 0.0_ldp, &
       0.0_ldp, 0.0_ldp, 3.0_ldp], [n,n])), &
      'only the isolated zero equation is pinned')
   call check(all(rhs == [4.0_ldp, 0.0_ldp, 9.0_ldp]), &
      'pinning preserves the right-hand side')

   c_mat = 0.0_ldp
   c_mat(1,1) = 1.0_ldp
   c_mat(2,1) = 2.0_ldp
   c_mat(3,3) = 3.0_ldp
   rhs = [1.0_ldp, 5.0_ldp, 0.0_ldp]
   saved_mat = c_mat
   saved_rhs = rhs
   call cmf_blkband_pin_null_v1(c_mat, rhs, n, n_pin)
   call check(n_pin == 0 .and. all(c_mat == saved_mat) .and. &
      all(rhs == saved_rhs), 'zero column keeps its nonzero equation')

   c_mat = 0.0_ldp
   c_mat(2,2) = 1.0_ldp
   c_mat(3,3) = 3.0_ldp
   rhs = [1.0_ldp, 0.0_ldp, 0.0_ldp]
   saved_mat = c_mat
   saved_rhs = rhs
   call cmf_blkband_pin_null_v1(c_mat, rhs, n, n_pin)
   call check(n_pin == 0 .and. all(c_mat == saved_mat) .and. &
      all(rhs == saved_rhs), 'zero row keeps its nonzero residual')

   c_mat = 0.0_ldp
   c_mat(2,1) = 1.0_ldp
   c_mat(2,2) = 2.0_ldp
   c_mat(3,3) = 3.0_ldp
   rhs = 0.0_ldp
   saved_mat = c_mat
   call cmf_blkband_pin_null_v1(c_mat, rhs, n, n_pin)
   call check(n_pin == 0 .and. all(c_mat == saved_mat), &
      'zero row with a live column is unchanged')

   c_mat = 0.0_ldp
   c_mat(1,1) = 1.0e-90_ldp
   c_mat(2,2) = 2.0_ldp
   c_mat(3,3) = 3.0_ldp
   rhs = 0.0_ldp
   saved_mat = c_mat
   call cmf_blkband_pin_null_v1(c_mat, rhs, n, n_pin)
   call check(n_pin == 0 .and. all(c_mat == saved_mat), &
      'near-zero live diagonal is unchanged')

! --- WP-7c symmetric elimination (cmf_blkband_pin_sym_v1) ---------
! Manufactured near-null trace column: a fourth unknown couples to
! the 3x3 major block at ~1e-43 (a floor-pinned trace level,
! POPS ~ 1e-31 scaled).  After symmetric elimination with XBAR=0
! the block solve must stay finite, the three major corrections
! must match a reference solve of the 3x3 block alone, and the NC
! row residual (sum of level corrections minus the conserved
! source) must hold to 1e-12 relative.

   block
      integer, parameter :: m = 4, mm = 3
      real(kind=ldp) :: blk(m,m), r(m), xbar(m), ref(mm,mm), rr(mm)
      logical :: pinmask(m)
      integer :: ipv(m), info, np
      real(kind=ldp) :: trace_c

      trace_c = 1.0e-43_ldp
      blk = 0.0_ldp
      blk(1,1) = 2.0_ldp; blk(1,2) = -0.5_ldp
      blk(2,1) = -0.5_ldp; blk(2,2) = 2.0_ldp; blk(2,3) = -0.25_ldp
      blk(3,1) = 1.0_ldp; blk(3,2) = 1.0_ldp; blk(3,3) = 1.0_ldp
      blk(4,4) = trace_c
      blk(1,4) = trace_c; blk(4,1) = trace_c
      r = [1.0_ldp, 0.5_ldp, 2.5_ldp, 0.0_ldp]
      xbar = 0.0_ldp
      pinmask = .false.; pinmask(4) = .true.

      ref = blk(1:mm,1:mm)
      rr = r(1:mm)
      call dgesv(mm, 1, ref, mm, ipv, rr, mm, info)
      call check(info == 0, 'sym-pin: reference 3x3 solve succeeds')

      call cmf_blkband_pin_sym_v1(blk, r, xbar, pinmask, m, np)
      call check(np == 1, 'sym-pin: one unknown eliminated')
      call check(blk(4,4) == 1.0_ldp .and. r(4) == 0.0_ldp, &
         'sym-pin: pinned row is the identity prescription')
      call check(all(blk(1:3,4) == 0.0_ldp), &
         'sym-pin: trace column removed from majors')

      call dgesv(m, 1, blk, m, ipv, r, m, info)
      call check(info == 0, 'sym-pin: pinned block solve succeeds')
      call check(all(r == r), 'sym-pin: solution finite')
      call check(maxval(abs(r(1:mm) - rr)) < 1.0e-12_ldp, &
         'sym-pin: majors match the no-trace reference')
      call check(r(4) == 0.0_ldp, &
         'sym-pin: pinned correction is the prescribed zero')

! Nonzero XBAR: pinned value must carry its column into the RHS.
! x4=2 adds trace_c*2 to row residuals 1 and 4 before the row is
! rewritten, so majors match the reference with rhs1 shifted.
      blk = 0.0_ldp
      blk(1,1) = 2.0_ldp; blk(1,2) = -0.5_ldp
      blk(2,1) = -0.5_ldp; blk(2,2) = 2.0_ldp; blk(2,3) = -0.25_ldp
      blk(3,1) = 1.0_ldp; blk(3,2) = 1.0_ldp; blk(3,3) = 1.0_ldp
      blk(4,4) = 1.0_ldp
      blk(1,4) = 0.5_ldp; blk(4,1) = 0.5_ldp
      r = [1.0_ldp, 0.5_ldp, 2.5_ldp, 0.0_ldp]
      xbar = 0.0_ldp; xbar(4) = 2.0_ldp

      ref = blk(1:mm,1:mm)
      rr = r(1:mm) - blk(1:mm,4)*xbar(4)
      call dgesv(mm, 1, ref, mm, ipv, rr, mm, info)
      call cmf_blkband_pin_sym_v1(blk, r, xbar, pinmask, m, np)
      call dgesv(m, 1, blk, m, ipv, r, m, info)
      call check(info == 0, 'sym-pin: carried block solve succeeds')
      call check(maxval(abs(r(1:mm) - rr)) < 1.0e-12_ldp, &
         'sym-pin: column carry keeps majors exact')
      call check(r(4) == 2.0_ldp, &
         'sym-pin: pinned correction is the prescribed value')
   end block

   if (failures /= 0) stop 1
   write(*,'(A)') 'CMF_BLKBAND_PIN_NULL_V1: PASS'
   write(*,'(A)') 'CMF_BLKBAND_PIN_SYM_V1: PASS'

contains

   subroutine check(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label

      if (ok) then
         write(*,'(A)') 'PASS: '//trim(label)
      else
         failures = failures + 1
         write(*,'(A)') 'FAIL: '//trim(label)
      end if
   end subroutine check

end program drv_cmf_blkband_pin
