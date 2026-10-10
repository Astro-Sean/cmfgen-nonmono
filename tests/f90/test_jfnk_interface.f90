! Integration test for the JFNK-CMFGEN interface layer.  Links the
! real jfnk_cmfgen_interface_v1.o / jfnk_cmfgen_int_mod.o from the
! production library and drives the state-free routines directly:
!
!   JFNK_BA_RAD_MATVEC_V1  -- the CMAT_RAD * x_frac operator piece
!   JFNK_MASK_FIXED_UNKNOWNS -- fixed T / Ne masking
!   JFNK_DOT_V1            -- RSCALE-weighted inner product
!
! These routines depend only on JFNK_CMFGEN_INT_MOD state (set here
! directly), so no CMFGEN atmosphere is required.

program test_jfnk_interface
   use jfnk_cmfgen_int_mod
   use jfnk_solver_v1_mod, only: jfnk_log_unknown
   implicit none
   external :: JFNK_MASK_FIXED_UNKNOWNS, JFNK_BA_RAD_MATVEC_V1
   double precision, external :: JFNK_DOT_V1
   integer :: luer, npass, nfail, n, i, k
   double precision, allocatable :: v(:), y(:), yref(:), a(:), b(:)
   double precision :: denom

   luer = 6
   npass = 0
   nfail = 0

   write(luer,'(A)') '================================================'
   write(luer,'(A)') ' JFNK interface unit tests (real library code)'
   write(luer,'(A)') '================================================'

   ! Layout: NT=4 (levels 1..NT-2, Ne=NT-1, T=NT), ND=3.
   JFNK_NT = 4
   JFNK_ND = 3
   JFNK_NION = 2
   JFNK_NUM_BNDS = 1
   JFNK_DIAG_INDX = 1
   JFNK_USE_BANDED_LOC = .false.
   n = JFNK_NT * JFNK_ND
   allocate(v(n), y(n), yref(n), a(n), b(n))

   !----------------------------------------------------------------
   ! I1: matvec in fractional mode: y_k = CMAT_RAD(:,:,k) * x_k
   !----------------------------------------------------------------
   jfnk_log_unknown = .false.
   JFNK_FIXED_T_LOC = .false.
   JFNK_FIXED_NE_LOC = .false.
   JFNK_FIX_T_D_ST = 0
   JFNK_FIX_T_D_END = 0
   allocate(JFNK_CMAT_RAD(JFNK_NT, JFNK_NT, JFNK_ND))
   JFNK_CMAT_RAD = 0d0
   do k = 1, JFNK_ND
      do i = 1, JFNK_NT
         JFNK_CMAT_RAD(i, i, k) = 10d0*k + i          ! diagonal
      end do
      JFNK_CMAT_RAD(1, 2, k) = 0.5d0*k               ! one off-diag
   end do
   JFNK_BA_RAD_READY = .true.
   v = 0d0
   do i = 1, n
      v(i) = 0.01d0 * i
   end do
   y = -99d0
   call JFNK_BA_RAD_MATVEC_V1(v, y)
   yref = 0d0
   do k = 1, JFNK_ND
      yref((k-1)*JFNK_NT+1:k*JFNK_NT) = v((k-1)*JFNK_NT+1:k*JFNK_NT)
      ! T is always log-mapped (residual applies T = T0*exp(x_T)),
      ! so its frac is exp(x_T)-1 even under the fractional pop map.
      yref(k*JFNK_NT) = exp(v(k*JFNK_NT)) - 1d0
      yref((k-1)*JFNK_NT+1:k*JFNK_NT) = &
         matmul(JFNK_CMAT_RAD(:,:,k), yref((k-1)*JFNK_NT+1:k*JFNK_NT))
   end do
   call check('I1 frac-mode matvec = CMAT*x (T log-mapped)', &
              maxval(abs(y - yref)) < 1d-14)

   !----------------------------------------------------------------
   ! I2: log mode: x_frac = exp(x)-1, so y = CMAT*(e^x - 1)
   !----------------------------------------------------------------
   jfnk_log_unknown = .true.
   y = -99d0
   call JFNK_BA_RAD_MATVEC_V1(v, y)
   do k = 1, JFNK_ND
      yref((k-1)*JFNK_NT+1:k*JFNK_NT) = &
         matmul(JFNK_CMAT_RAD(:,:,k), &
                exp(v((k-1)*JFNK_NT+1:k*JFNK_NT)) - 1d0)
   end do
   call check('I2 log-mode matvec = CMAT*(e^x-1)', &
              maxval(abs(y - yref)) < 1d-14)

   !----------------------------------------------------------------
   ! I3: matvec masks fixed unknowns
   !----------------------------------------------------------------
   jfnk_log_unknown = .false.
   JFNK_FIXED_T_LOC = .true.
   JFNK_FIXED_NE_LOC = .true.
   y = -99d0
   call JFNK_BA_RAD_MATVEC_V1(v, y)
   ! fixed rows contribute nothing: reference uses x with fixed comps 0
   yref = 0d0
   do k = 1, JFNK_ND
      yref((k-1)*JFNK_NT+1:k*JFNK_NT) = v((k-1)*JFNK_NT+1:k*JFNK_NT)
      yref(k*JFNK_NT-1) = 0d0
      yref(k*JFNK_NT)   = 0d0
      yref((k-1)*JFNK_NT+1:k*JFNK_NT) = &
         matmul(JFNK_CMAT_RAD(:,:,k), yref((k-1)*JFNK_NT+1:k*JFNK_NT))
   end do
   call check('I3 matvec zeroes fixed-unknown rows', &
              maxval(abs(y - yref)) < 1d-14)
   JFNK_FIXED_T_LOC = .false.
   JFNK_FIXED_NE_LOC = .false.

   !----------------------------------------------------------------
   ! I4: MASK_FIXED_UNKNOWNS: full FIX_T + FIX_NE
   !----------------------------------------------------------------
   JFNK_FIXED_T_LOC = .true.
   JFNK_FIXED_NE_LOC = .true.
   v = 1d0
   call JFNK_MASK_FIXED_UNKNOWNS(n, v)
   call check('I4 Ne components masked', &
              all(v((/(k*JFNK_NT-1, k=1,JFNK_ND)/)) == 0d0))
   call check('I4 T components masked', &
              all(v((/(k*JFNK_NT,   k=1,JFNK_ND)/)) == 0d0))
   call check('I4 pop components kept', &
              all(v(1:2) == 1d0))

   !----------------------------------------------------------------
   ! I5: partial T-fix range (FIX_T_D_ST..END), FIX_T off
   !----------------------------------------------------------------
   JFNK_FIXED_T_LOC = .false.
   JFNK_FIX_T_D_ST = 2
   JFNK_FIX_T_D_END = 2
   v = 1d0
   call JFNK_MASK_FIXED_UNKNOWNS(n, v)
   call check('I5 T masked only in range (K=2)', &
              v(2*JFNK_NT) == 0d0 .and. v(JFNK_NT) == 1d0 &
              .and. v(3*JFNK_NT) == 1d0)
   JFNK_FIX_T_D_ST = 0
   JFNK_FIX_T_D_END = 0
   JFNK_FIXED_NE_LOC = .false.

   !----------------------------------------------------------------
   ! I6: weighted dot: sum (w*a)(w*b); unweighted when not ready
   !----------------------------------------------------------------
   allocate(JFNK_RSCALE(n))
   JFNK_RSCALE = 2d0
   JFNK_RSCALE_READY = .true.
   a = 1d0
   b = 3d0
   call check('I6 weighted dot = 4*sum(a*b)', &
              abs(JFNK_DOT_V1(n, a, b) - 4d0*3d0*n) < 1d-10)
   JFNK_RSCALE_READY = .false.
   call check('I6 unweighted fallback = sum(a*b)', &
              abs(JFNK_DOT_V1(n, a, b) - 3d0*n) < 1d-10)
   JFNK_RSCALE_READY = .true.

   write(luer,'(A)') '================================================'
   write(luer,'(A,I3,A,I3)') ' RESULTS: ', npass, ' passed, ', nfail, ' failed'
   write(luer,'(A)') '================================================'
   if (nfail > 0) stop 1

contains

   subroutine check(name, ok)
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      if (ok) then
         npass = npass + 1
         write(luer,'(A,A)') '  PASS  ', name
      else
         nfail = nfail + 1
         write(luer,'(A,A)') '  FAIL  ', name
      end if
   end subroutine check

end program test_jfnk_interface
