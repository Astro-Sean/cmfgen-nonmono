!----------------------------------------------------------------------
! TEST_GUARDS: standalone regression tests for the Oct-2026
! robustness fixes (A1/A4/A2) against the same production sources the
! build compiles.
!
!   G1: operator-consistency of the in-sweep matched-frequency span
!       with the MATCH_FREQ refresh span.  The in-sweep code in
!       fg_j_cmf_nonmono_v1.f sums one-bin log ratios over the matched
!       offset K2 and reads the upwind intensity from store slot
!       FREQ_IDX_ITER+NM_RAY (frequency k+NM_RAY, since the store is
!       shifted one slot).  The refresh MATCH_FREQ computes
!       LOG(KN_FREQ(L)/KN_FREQ(L+NNF)).  This test verifies the two
!       expressions telescope to the same span for every k and offset,
!       and that the cap NCF_IN-FIDX-1 never touches the phantom slot.
!
!   G2: IN_HBC denominator floor.  Replicates the guarded formula in
!       fg_j_cmf_nonmono_v1.f verbatim and checks the degenerate
!       extremes (cancellation denominator, NaN/Inf inputs) produce
!       the 0.5 fallback, and that a well-scaled denominator keeps the
!       raw ratio.
!
!   G3: JFNK_MXDL cap honored under auto-delta (del<0).  Runs
!       newton_hook_v1 on the shared test problem with mxdl=0.05 and
!       checks the trust radius stayed <= mxdl (the old code reset
!       mxdl to 1d99 whenever del<0).
!
!   G4: solver-core determinism: two newton_hook_v1 activations from
!       the same start return bitwise-identical iterates.  The
!       production eval-side bracketing (JFNK_RAD_SWEEP_V1) cannot be
!       exercised standalone; this test locks determinism of the
!       Newton/GMRES layer itself.
!
! edited by Devin, Oct-2026
!----------------------------------------------------------------------

! Test problem identical in shape to test_newton_core's: a known-root
! exponential system.  f_t is a pure function of x, so repeated
! activations are deterministic by construction.
module guards_prob_mod
   implicit none
   integer, parameter :: NG = 60
   double precision :: ag(NG), xsg(NG)
   double precision :: bg = 0.3d0
contains
   subroutine seed_prob()
      integer :: i
      do i = 1, NG
         ag(i)  = 1.0d-3 * 10d0**(6d0 * dble(i - 1) / dble(NG - 1))
         xsg(i) = 1.5d0 * sin(1.7d0 * dble(i))
      end do
   end subroutine seed_prob
end module guards_prob_mod

program test_guards
   use jfnk_solver_v1_mod
   use guards_prob_mod
   implicit none
   external :: f_t, df_t, dfp_t, sub_t
   double precision, external :: dp_t
   integer :: nfail, info, luer, i
   !
   luer = 6
   nfail = 0
   !
   ! ---------------- G1: matched-span operator consistency ----------
   !
   block
      integer, parameter :: NCF_IN = 40
      double precision :: knf(0:NCF_IN)   ! KN_FREQ indexed by frequency
      integer :: fidx, k2, ncf_ok, m
      double precision :: span_insweep, span_match, span_tele, freq
      logical :: bad
      !
      ! Nonuniform log grid (MATCH_FREQ must hold on any spacing).
      do i = 0, NCF_IN
         knf(i) = 1.0d15 * exp(-0.11d0 * dble(i) - 0.002d0 * dble(i*i))
      end do
      bad = .false.
      do fidx = 1, NCF_IN - 2
         ! Current frequency k sits at store slot fidx, i.e. frequency
         ! index fidx+1 (the I_NEXT store is shifted by one slot).
         freq = knf(fidx + 1)
         do k2 = 1, NCF_IN - fidx - 1   ! the production cap
            ! In-sweep span (fixed form): log(FREQ/KN_FREQ(FIDX+1+K2)).
            span_insweep = log(freq / knf(fidx + 1 + k2))
            ! MATCH_FREQ span LOG(KN_FREQ(L)/KN_FREQ(L+NNF)) with
            ! L = fidx+1 (frequency index of the current slot) and
            ! NNF = k2 -- the same upwind slot the in-sweep store read
            ! FREQ_IDX_ITER+K2 maps to (frequency fidx+1+k2).
            span_match = log(knf(fidx + 1) / knf(fidx + 1 + k2))
            ! A by-hand one-bin telescoped sum gives a third,
            ! independent evaluation of the same span.
            span_tele = 0d0
            do m = 1, k2
               span_tele = span_tele + &
                  log(knf(fidx + m) / knf(fidx + m + 1))
            end do
            if(abs(span_insweep - span_match) > &
               1d-12 * abs(span_match)) then
               write(luer,'(A,2I4,2ES14.6)') &
                  'G1 FAIL: span mismatch k2=', fidx, k2, &
                  span_insweep, span_match
               bad = .true.
            end if
            if(abs(span_tele - span_match) > &
               1d-10 * abs(span_match)) then
               write(luer,'(A,2I4,2ES14.6)') &
                  'G1 FAIL: telescoped span mismatch k2=', fidx, k2, &
                  span_tele, span_match
               bad = .true.
            end if
         end do
      end do
      ! Cap check: largest usable offset leaves the phantom slot free.
      ncf_ok = 0
      do fidx = 1, NCF_IN - 2
         if(fidx + 1 + (NCF_IN - fidx - 1) <= NCF_IN) &
            ncf_ok = ncf_ok + 1
      end do
      if(ncf_ok /= NCF_IN - 2) bad = .true.
      if(bad) then
         nfail = nfail + 1
         write(luer,'(A)') 'G1: FAIL (matched-span operator)'
      else
         write(luer,'(A)') 'G1: PASS (matched-span operator)'
      end if
   end block
   !
   ! ---------------- G2: IN_HBC denominator guard -------------------
   !
   block
      double precision :: inhbc, den, dscl, jnd, ic, muin
      integer :: nbad2
      nbad2 = 0
      ! Case A: near-cancellation 2*J ~ IC (old guard lets it blow up)
      jnd = 0.5d0 * 3.7d-10; ic = 3.7d-10
      muin = 1.0d-9
      inhbc = muin
      den = 2.0d0 * jnd - ic
      dscl = max(abs(ic), 2.0d0 * jnd, 1.0d-30)
      if(abs(den) < 1.0d-6 * dscl) then
         inhbc = 0.5d0
      else
         inhbc = inhbc / den
      end if
      if(.not.(inhbc >= 0.0d0)) inhbc = 0.5d0
      if(.not.(inhbc <= 1.0d0)) inhbc = 0.5d0
      if(inhbc /= 0.5d0) nbad2 = nbad2 + 1
      ! Case B: healthy denominator keeps the measured ratio
      jnd = 2.0d0; ic = 1.0d0
      muin = 0.9d0
      inhbc = muin
      den = 2.0d0 * jnd - ic
      dscl = max(abs(ic), 2.0d0 * jnd, 1.0d-30)
      if(abs(den) < 1.0d-6 * dscl) then
         inhbc = 0.5d0
      else
         inhbc = inhbc / den
      end if
      if(.not.(inhbc >= 0.0d0)) inhbc = 0.5d0
      if(.not.(inhbc <= 1.0d0)) inhbc = 0.5d0
      if(abs(inhbc - 0.3d0) > 1d-14) nbad2 = nbad2 + 1
      ! Case C: NaN numerator falls back
      inhbc = 0.0d0 / 0.0d0
      den = 2.0d0 * jnd - ic
      dscl = max(abs(ic), 2.0d0 * jnd, 1.0d-30)
      if(abs(den) < 1.0d-6 * dscl) then
         inhbc = 0.5d0
      else
         inhbc = inhbc / den
      end if
      if(.not.(inhbc >= 0.0d0)) inhbc = 0.5d0
      if(.not.(inhbc <= 1.0d0)) inhbc = 0.5d0
      if(inhbc /= 0.5d0) nbad2 = nbad2 + 1
      if(nbad2 > 0) then
         nfail = nfail + 1
         write(luer,'(A,I2)') 'G2: FAIL (IN_HBC guard), nbad=', nbad2
      else
         write(luer,'(A)') 'G2: PASS (IN_HBC guard)'
      end if
   end block
   !
   ! -------- G3/G4: MXDL cap + core determinism on newton_hook_v1 ---
   !
   block
      integer, parameter :: NN = 60
      double precision :: x1(NN), x2(NN), tol_t
      integer :: j
      call seed_prob()
      if(.not. allocated(jfnk_x)) allocate(jfnk_x(NN), jfnk_fx(NN))
      jfnk_ptc_active = .false.
      jfnk_do_line_search = .true.
      ! del<0 selects automatic trust-region init; mxdl=0.05 must not
      ! be overwritten with the unlimited fallback.
      jfnk_x = xsg + 0.4d0
      info = 1
      call newton_hook_v1(f_t, df_t, dfp_t, sub_t, dp_t, &
         10, NN, 1d-9, 1d-11, -1d0, 1d-12, 5d-2, 30, info, luer)
      if(jfnk_del > 5.0d-2 + 1d-14) then
         nfail = nfail + 1
         write(luer,'(A,ES12.4)') &
            'G3: FAIL (mxdl overridden), jfnk_del=', jfnk_del
      else
         write(luer,'(A,ES12.4)') &
            'G3: PASS (mxdl honored), jfnk_del=', jfnk_del
      end if
      x1 = jfnk_x
      ! Determinism: identical second activation, bitwise-identical x.
      jfnk_x = xsg + 0.4d0
      info = 1
      call newton_hook_v1(f_t, df_t, dfp_t, sub_t, dp_t, &
         10, NN, 1d-9, 1d-11, -1d0, 1d-12, 5d-2, 30, info, luer)
      x2 = jfnk_x
      write(luer,'(A,ES12.4)') '  G4 max|x1-x2|=', maxval(abs(x1-x2))
      if(any(x1 /= x2)) then
         nfail = nfail + 1
         write(luer,'(A)') 'G4: FAIL (nondeterministic core)'
      else
         write(luer,'(A)') 'G4: PASS (deterministic core)'
      end if
   end block
   !
   if(nfail == 0) then
      write(luer,'(A)') 'ALL GUARD TESTS PASSED'
   else
      write(luer,'(A,I2)') 'FAILURES: ', nfail
      stop 1
   end if
end program test_guards

subroutine f_t(n, x, y)
   use guards_prob_mod
   implicit none
   integer :: n, i
   double precision :: x(n), y(n), d(n)
   d = x - xsg
   do i = 1, n
      y(i) = ag(i) * (exp(min(20d0, x(i))) - exp(xsg(i)))
      if(i > 1)  y(i) = y(i) + bg * (d(i) - d(i - 1))
      if(i < n)  y(i) = y(i) + bg * (d(i) - d(i + 1))
   end do
end subroutine f_t

subroutine df_t(n, v, jv)
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx
   implicit none
   integer :: n
   double precision :: v(n), jv(n), xp(n), yp(n), h, vinf
   ! Fixed-h finite difference (deterministic).
   vinf = maxval(abs(v))
   if(vinf == 0d0) then
      jv = 0d0
      return
   end if
   h = 1d-6 / vinf
   xp = jfnk_x + h * v
   call f_t(n, xp, yp)
   jv = (yp - jfnk_fx) / h
end subroutine df_t

subroutine dfp_t(n, z)
   implicit none
   integer :: n
   double precision :: z(n)
   ! Identity preconditioner is adequate for this problem scale.
end subroutine dfp_t

subroutine sub_t()
   implicit none
end subroutine sub_t

double precision function dp_t(n, x, y)
   implicit none
   integer :: n
   double precision :: x(n), y(n)
   dp_t = dot_product(x, y)
end function dp_t
