program drv_cmf_blkband_pin
   use set_kind_module
   implicit none
   integer, parameter :: n = 3
   real(kind=ldp) :: c_mat(n,n), rhs(n), saved_mat(n,n), saved_rhs(n)
   integer :: n_pin, failures

   external :: cmf_blkband_pin_null_v1

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

   if (failures /= 0) stop 1
   write(*,'(A)') 'CMF_BLKBAND_PIN_NULL_V1: PASS'

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
