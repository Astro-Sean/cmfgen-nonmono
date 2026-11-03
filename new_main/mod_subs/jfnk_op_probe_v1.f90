! jfnk_op_probe_v1.f90 -- JFNK operator-consistency battery.
!
! Runs the production-interface linearity probes against whatever
! residual/JVP callbacks the caller supplies: the production
! JFNK_RESIDUAL_V1/JFNK_JVP_V1 pair at a live admission (via the
! JFNK_HSCAN hook in JFNK_SOLVE_CORR_V1), or the surrogate
! se_prob_mod pair in the unit-test driver.  Keeping the battery
! callback-based means the same measurement code is exercised in
! both places; only the operators differ.
!
! Emitted lines all begin ' OPPROBE ' and are k=v parseable:
!   BEGIN / REP (baseline repeatability) / REPP (perturbed)
!   DIR (per-direction FD vs reference + policy metadata)
!   HS  (per-direction FD h-sweep rows)
!   ADD / HOM / END (eval accounting)
!
! Unknown layout assumed depth-major: idx=(k-1)*nt+j, j=nt-1 -> ED,
! j=nt -> T (matches JFNK_NT ordering and the surrogate se layout).

module jfnk_op_probe_mod
   implicit none
   abstract interface
      subroutine op_resid(n, x, f)
         integer, intent(in) :: n
         double precision, intent(in) :: x(n)
         double precision, intent(out) :: f(n)
      end subroutine op_resid
      ! jvpa wraps the production JVP call and reports the effective
      ! step and policy events so the battery can separate operator
      ! error from h-selection noise.
      subroutine op_jvpa(n, v, jv, h_eff, capped, nretr, nclip, ok)
         integer, intent(in) :: n
         double precision, intent(in) :: v(n)
         double precision, intent(out) :: jv(n)
         double precision, intent(out) :: h_eff
         logical, intent(out) :: capped, ok
         integer, intent(out) :: nretr, nclip
      end subroutine op_jvpa
      ! bamv is the independent reference derivative action: the BA
      ! diagonal-band matvec on production (approximate), or an exact
      ! assembled Jacobian on the surrogate.
      subroutine op_bamv(n, v, jv)
         integer, intent(in) :: n
         double precision, intent(in) :: v(n)
         double precision, intent(out) :: jv(n)
      end subroutine op_bamv
      function op_dotp(n, a, b) result(d)
         integer, intent(in) :: n
         double precision, intent(in) :: a(n), b(n)
         double precision :: d
      end function op_dotp
   end interface
contains

   ! Deterministic LCG so probe directions are reproducible across
   ! runs and compilers.  Returns U(0,1); state carried in s.
   double precision function op_rand(s)
      implicit none
      integer(kind=8), intent(inout) :: s
      s = mod(6364136223846793005_8*s + 1442695040888963407_8, &
              9223372036854775807_8)
      if (s <= 0_8) s = s + 9223372036854775806_8
      op_rand = dble(s) / 9223372036854775807d0
   end function op_rand

   double precision function op_wnorm(n, v, dotp)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: v(n)
      procedure(op_dotp) :: dotp
      op_wnorm = sqrt(max(0d0, dotp(n, v, v)))
   end function op_wnorm

   ! Build probe direction isel in v, unit inf-norm.  Slot j<nt-1 are
   ! populations, j=nt-1 ED, j=nt T; depth k owns slots
   ! (k-1)*nt+1 .. k*nt.
   subroutine op_mkdir(n, nt, nd, f0, dotp, isel, v, tag)
      implicit none
      integer, intent(in) :: n, nt, nd, isel
      double precision, intent(in) :: f0(n)
      procedure(op_dotp) :: dotp
      double precision, intent(out) :: v(n)
      character(len=*), intent(out) :: tag
      integer(kind=8) :: s
      integer :: i, j, k, kk, kmax, ncomp
      double precision :: fmax, sc
      v = 0d0
      tag = 'unk'
      s = 1103515245_8 + 97_8*isel
      select case (isel)
      case (1, 2)                       ! full random dirs
         do i = 1, n
            v(i) = 2d0*op_rand(s) - 1d0
         end do
         write(tag, '(a,i1)') 'rnd', isel
      case (3)                          ! populations only
         do i = 1, n
            j = mod(i - 1, nt) + 1
            if (j <= nt - 2) v(i) = 2d0*op_rand(s) - 1d0
         end do
         tag = 'pops'
      case (4)                          ! ED only
         do i = 1, n
            j = mod(i - 1, nt) + 1
            if (j == nt - 1) v(i) = 2d0*op_rand(s) - 1d0
         end do
         tag = 'ed'
      case (5)                          ! T only
         do i = 1, n
            j = mod(i - 1, nt) + 1
            if (j == nt) v(i) = 2d0*op_rand(s) - 1d0
         end do
         tag = 't'
      case (6)                          ! single depth 1
         kk = 1
         do i = (kk - 1)*nt + 1, kk*nt
            v(i) = 2d0*op_rand(s) - 1d0
         end do
         tag = 'k1'
      case (7)                          ! single mid depth
         kk = max(1, (nd + 1)/2)
         do i = (kk - 1)*nt + 1, kk*nt
            v(i) = 2d0*op_rand(s) - 1d0
         end do
         tag = 'kmid'
      case (8)                          ! worst-residual depth
         fmax = 0d0
         kmax = 1
         do i = 1, n
            if (abs(f0(i)) > fmax) then
               fmax = abs(f0(i))
               kmax = (i - 1)/nt + 1
            end if
         end do
         do i = (kmax - 1)*nt + 1, min(kmax*nt, n)
            v(i) = 2d0*op_rand(s) - 1d0
         end do
         tag = 'kworst'
      case (9)                          ! three depths spread
         do kk = 1, 3
            k = 1 + (kk - 1)*max(1, (nd - 1)/2)
            k = min(k, nd)
            do i = (k - 1)*nt + 1, k*nt
               v(i) = 2d0*op_rand(s) - 1d0
            end do
         end do
         tag = 'blk3'
      case (10)                         ! single basis vector
         i = 1 + int(op_rand(s)*dble(n - 1))
         v(max(1, i)) = 1d0
         tag = 'basis'
      case (11)                         ! ED + T together
         do i = 1, n
            j = mod(i - 1, nt) + 1
            if (j >= nt - 1) v(i) = 2d0*op_rand(s) - 1d0
         end do
         tag = 'edt'
      case (12)                         ! spike: forces tiny h so the
         ! perturbed eval likely crosses a clip margin or cap
         i = 1 + int(op_rand(s)*dble(n - 1))
         v(max(1, i)) = 1d3
         tag = 'spike'
      case (13)                         ! cross-depth checkerboard
         ! on the first population slot of each depth
         do k = 1, nd
            i = (k - 1)*nt + 1
            if (i <= n) v(i) = (-1d0)**k
         end do
         tag = 'xdep'
      end select
      ncomp = count(abs(v) > 0d0)
      if (ncomp == 0) then
         v(1) = 1d0
         tag = 'fb1'
      end if
      sc = maxval(abs(v))
      if (sc > 0d0) v = v/sc
   end subroutine op_mkdir

   ! Forward-difference JVP at an explicit step (independent of the
   ! production h policy -- the battery's own measurement).
   subroutine op_fd(n, x0, f0, v, h, resid, jv)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: x0(n), f0(n), v(n), h
      double precision, intent(out) :: jv(n)
      procedure(op_resid) :: resid
      double precision :: xp(n), fp(n)
      xp = x0 + h*v
      call resid(n, xp, fp)
      jv = (fp - f0)/h
   end subroutine op_fd

   subroutine jfnk_op_probe_v1(luer, n, x0, f0, nt, nd, resid, &
                               jvpa, bamv, dotp, mode)
      implicit none
      integer, intent(in) :: luer, n, nt, nd
      double precision, intent(in) :: x0(n), f0(n)
      procedure(op_resid) :: resid
      procedure(op_jvpa) :: jvpa
      procedure(op_bamv) :: bamv
      procedure(op_dotp) :: dotp
      character(len=*), intent(in) :: mode

      integer, parameter :: NDIR = 13, NREP = 4, NSCAN = 5
      double precision, parameter :: SCAN_M(NSCAN) = &
         (/1d-2, 1d-1, 1d0, 1d1, 1d2/)
      integer :: i, isel, irep, isc, nretr, nclip
      integer :: n_ev, n_jvp, n_fail, n_cap, n_retr, n_clip
      double precision :: v(n), u(n), w(n), jv(n), jv2(n), jvu(n)
      double precision :: jvv(n), bv(n), fref(n), frep(n), xp(n)
      double precision :: fp(n), fm(n), jvc(n), fdu(n), fdv(n)
      double precision :: fduv(n), heff, h_fd, dif, dn, nv, nbv
      double precision :: efd, efc, emin, ef, d_pol, d_op, d2_pol
      double precision :: d2_op, umin, hu, hv, hw
      logical :: capped, ok
      character(len=8) :: tag, t2

      n_ev = 0
      n_jvp = 0
      n_fail = 0
      n_cap = 0
      n_retr = 0
      n_clip = 0
      if (luer <= 0) return

      write(luer, '(a,a,a,i8,a,i4,a,i4)') ' OPPROBE BEGIN mode=', &
         trim(mode), ' n=', n, ' nt=', nt, ' nd=', nd

      !--- Repeatability at the baseline point ---------------------
      fref = f0
      do irep = 1, NREP
         call resid(n, x0, frep)
         n_ev = n_ev + 1
         dif = maxval(abs(frep - fref))
         dn = max(1d-300, maxval(abs(fref)))
         write(luer, '(a,i2,a,es12.4,a,es12.4)') ' OPPROBE REP k=', &
            irep, ' df_inf=', dif, ' df_rel=', dif/dn
      end do

      !--- Direction battery ---------------------------------------
      do isel = 1, NDIR
         call op_mkdir(n, nt, nd, f0, dotp, isel, v, tag)
         call jvpa(n, v, jv, heff, capped, nretr, nclip, ok)
         n_jvp = n_jvp + 1
         n_ev = n_ev + 1
         if (.not. ok) then
            n_fail = n_fail + 1
            write(luer, '(a,a,a)') ' OPPROBE DIR tag=', tag, &
               ' fail=1'
            cycle
         end if
         if (capped) n_cap = n_cap + 1
         n_retr = n_retr + nretr
         n_clip = n_clip + nclip
         nv = max(1d-300, op_wnorm(n, jv, dotp))

         ! Perturbed-point repeatability on the first direction only.
         if (isel == 1) then
            xp = x0 + heff*v
            call resid(n, xp, fref)
            call resid(n, xp, frep)
            n_ev = n_ev + 2
            dif = maxval(abs(frep - fref))
            dn = max(1d-300, maxval(abs(fref)))
            write(luer, '(a,es12.4,a,es12.4)') &
               ' OPPROBE REPP df_inf=', dif, ' df_rel=', dif/dn
         end if

         ! Reference action and central-difference check at h_eff.
         call bamv(n, v, bv)
         nbv = max(1d-300, op_wnorm(n, bv, dotp))
         xp = x0 + heff*v
         call resid(n, xp, fp)
         xp = x0 - heff*v
         call resid(n, xp, fm)
         n_ev = n_ev + 2
         jvc = (fp - fm)/(2d0*heff)
         efd = op_wnorm(n, jv - bv, dotp)/nbv
         efc = op_wnorm(n, jvc - bv, dotp)/nbv

         ! h-sweep: forward-FD vs reference across heff*SCAN_M.
         ! The sweep separates truncation (linear in h) from the
         ! roundoff/noise floor (diverging as h -> 0) and exposes a
         ! plateau where the FD resolves the operator.
         emin = huge(emin)
         umin = 1d0
         do isc = 1, NSCAN
            h_fd = heff*SCAN_M(isc)
            call op_fd(n, x0, f0, v, h_fd, resid, w)
            n_ev = n_ev + 1
            ef = op_wnorm(n, w - bv, dotp)/nbv
            write(luer, '(a,a,a,es10.2,a,es12.4)') &
               ' OPPROBE HS tag=', tag, ' mult=', SCAN_M(isc), &
               ' ef=', ef
            if (ef < emin) then
               emin = ef
               umin = SCAN_M(isc)
            end if
         end do

         write(luer, '(a,a,a,es12.4,a,l1,a,i2,a,i3,3(a,es12.4))') &
            ' OPPROBE DIR tag=', tag, ' heff=', heff, ' capped=', &
            capped, ' retr=', nretr, ' clip=', nclip, ' njv=', nv, &
            ' eFD=', efd, ' eFC=', efc
         write(luer, '(a,a,2(a,es12.4))') ' OPPROBE DM tag=', tag, &
            ' emin=', emin, ' hstar=', umin
      end do

      !--- Additivity ----------------------------------------------
      ! d_pol measures the production JVP's observed nonlinearity
      ! (per-direction h policy included); d_op repeats the test with
      ! all three forward differences evaluated at one shared h so
      ! the operator nonlinearity is isolated from the h policy.
      call op_mkdir(n, nt, nd, f0, dotp, 1, u, t2)
      call op_mkdir(n, nt, nd, f0, dotp, 3, v, t2)
      call jvpa(n, u, jvu, hu, capped, nretr, nclip, ok)
      n_jvp = n_jvp + 1
      if (.not. ok) n_fail = n_fail + 1
      call jvpa(n, v, jvv, hv, capped, nretr, nclip, ok)
      n_jvp = n_jvp + 1
      if (.not. ok) n_fail = n_fail + 1
      ! w = u + v left unnormalized: additivity is J(u+v) = Ju + Jv,
      ! and rescaling w would only inject a spurious (1/m - 1) term
      ! into the defect.
      w = u + v
      call jvpa(n, w, jv2, hw, capped, nretr, nclip, ok)
      n_jvp = n_jvp + 1
      if (.not. ok) n_fail = n_fail + 1
      ! Defects are normalized by the operand norms: an absolute
      ! defect is meaningless against images of size ~1e14.
      dn = max(1d-300, op_wnorm(n, jvu, dotp) + &
               op_wnorm(n, jvv, dotp) + op_wnorm(n, jv2, dotp))
      d_pol = op_wnorm(n, jv2 - jvu - jvv, dotp)/dn
      ! Shared-h forward differences at the smallest returned step
      ! so no direction uses a step its own policy rejected; this
      ! isolates operator nonlinearity from the per-direction h
      ! policy baked into the JVP images above.
      h_fd = min(hu, min(hv, hw))
      call op_fd(n, x0, f0, u, h_fd, resid, fdu)
      call op_fd(n, x0, f0, v, h_fd, resid, fdv)
      call op_fd(n, x0, f0, w, h_fd, resid, fduv)
      n_ev = n_ev + 3
      dn = max(1d-300, op_wnorm(n, fdu, dotp) + &
               op_wnorm(n, fdv, dotp) + op_wnorm(n, fduv, dotp))
      d_op = op_wnorm(n, fduv - fdu - fdv, dotp)/dn
      write(luer, '(a,2(a,es12.4),a,es12.4)') ' OPPROBE ADD', &
         ' d_pol=', d_pol, ' d_op=', d_op, ' h=', h_fd

      !--- Homogeneity ---------------------------------------------
      w = 2d0*u
      call jvpa(n, w, jv2, hw, capped, nretr, nclip, ok)
      n_jvp = n_jvp + 1
      if (.not. ok) n_fail = n_fail + 1
      dn = max(1d-300, op_wnorm(n, jv2, dotp) + &
               2d0*op_wnorm(n, jvu, dotp))
      d2_pol = op_wnorm(n, jv2 - 2d0*jvu, dotp)/dn
      call op_fd(n, x0, f0, w, h_fd, resid, fduv)
      n_ev = n_ev + 1
      dn = max(1d-300, op_wnorm(n, fduv, dotp) + &
               2d0*op_wnorm(n, fdu, dotp))
      d2_op = op_wnorm(n, fduv - 2d0*fdu, dotp)/dn
      write(luer, '(a,2(a,es12.4))') ' OPPROBE HOM a=2', &
         ' d_pol=', d2_pol, ' d_op=', d2_op

      write(luer, '(a,6(a,i6))') ' OPPROBE END', ' evals=', n_ev, &
         ' jvps=', n_jvp, ' fail=', n_fail, ' caps=', n_cap, &
         ' retr=', n_retr, ' clip=', n_clip
      flush(luer)
   end subroutine jfnk_op_probe_v1

end module jfnk_op_probe_mod
