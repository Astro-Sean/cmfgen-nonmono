!----------------------------------------------------------------------
! TEST_SCHUR: the two-level (n,T) Schur-complement preconditioner
! jfnk_schur_form_v1 / jfnk_schur_apply_v1 / jfnk_blk_lsq_solve,
! linked from the production object.
!
! Inputs are built exactly as the production SCHUR path does
! (jfnk_cmfgen_interface_v1.f ~4580-4670): per depth, the population
! block B = Jnn + ridge*diag with reg rows identity-substituted, the
! normal matrix G = B^T B + lsq_lam^2 * blkmax^2 * I factored by
! DPOTRF into chol; braw = the RAW block; rsf=csf=1.
!
!   S1  exactness: v = A_full x_true -> apply recovers x on pops+T
!       (small lsq_lam), Ne row zeroed. Verifies the Schur algebra.
!   S2  rank-deficient Jnn + lsq_lam>0: LSQ silently regularizes the
!       rank deficiency -- apply completes with a pseudo-inverse.
!   S3  same matrix, lsq_lam=0: DPOTRF fails -> chol holds a partial
!       factor. DPOTRS does not detect non-PD input (info=0), so
!       schur_form's "if(ierr/=0) w=0" guard never fires and the
!       failure propagates silently into S and the apply output.
!   S4  reg-row pass-through: a null Jnn row (regularized to I) makes
!       blk_lsq return z_i ~ v_i on that row.
!
! edited 09-Nov-2026
!----------------------------------------------------------------------
program test_schur
   implicit none
   integer, parameter :: NT = 6, NP = 4, ND = 3, LUER = 6
   double precision :: a(NT,NT,ND), chol(NT,NT,ND), braw(NT,NT,ND)
   double precision :: ridge(NT,ND), rsf(NT,ND), csf(NT,ND)
   double precision :: s(ND,ND), v(ND*NT), x_true(ND*NT)
   integer :: sipiv(ND), info, nfail, k, i, ierr, nbad
   logical :: reg(NT,ND), lsq_ok(ND)
   double precision :: err, lsq_lam, zl(ND*NT)
   !
   nfail = 0
   rsf = 1d0
   csf = 1d0
   reg = .false.
   !
   ! ============ S1a: production ridge -> ~1% solve bias ============
   ! build_inputs uses ridge_add=1e-2 (the production ITRY=0 ridge);
   ! the Schur solve is then M^-1 A on a ridged Jnn, so ~1% deviation
   ! is the regularization floor, not an algebra error.
   do k = 1, ND
      call gen_block(a(1,1,k), NT, NP, k)
   end do
   lsq_lam = 1d-8
   call build_inputs(a, chol, braw, ridge, reg, rsf, csf, &
                     lsq_lam, 1d-2, lsq_ok, NT, NP, ND)
   call jfnk_schur_form_v1(ND, NT, NP, NT, chol, braw, ridge, &
        reg, rsf, csf, s, sipiv, info)
   write(luer,'(A,I3)') 'S1a: schur_form info=', info
   call seed_x(x_true, ND*NT)
   v = 0d0
   do k = 1, ND
      do i = 1, NT
         v((k-1)*NT+i) = dot_product(a(i,1:NT,k), &
              x_true((k-1)*NT+1:k*NT))
      end do
   end do
   call jfnk_schur_apply_v1(ND, NT, NP, NT, NT-1, chol, braw, &
        ridge, reg, rsf, csf, s, sipiv, v, ierr)
   err = 0d0
   do k = 1, ND
      do i = 1, NP
         err = max(err, abs(v((k-1)*NT+i) - x_true((k-1)*NT+i)))
      end do
      err = max(err, abs(v((k-1)*NT+NT) - x_true((k-1)*NT+NT)))
   end do
   write(luer,'(A,I3,A,ES12.4,A,ES12.4)') &
      'S1a: apply ierr=', ierr, '  max|z-x|=', err, &
      '  Ne out=', v(NT-1)
   if(info == 0 .and. ierr == 0 .and. err < 0.1d0 .and. &
      v(NT-1) == 0d0) then
      write(luer,'(A)') 'S1a: PASS (deviation = ridge bias)'
   else
      write(luer,'(A)') 'S1a: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ===== S1b: ridge=0, tiny lam -> algebraic exactness =====
   ! The Schur solve drops the Ne equation (z_Ne := 0), so it can
   ! only invert the (n,T) subsystem exactly -- x_true must have a
   ! zero Ne component or the leftover JnNe*x_Ne term masquerades
   ! as solve error (~1e-2 with these couplings).
   call build_inputs(a, chol, braw, ridge, reg, rsf, csf, &
                     1d-10, 0d0, lsq_ok, NT, NP, ND)
   call jfnk_schur_form_v1(ND, NT, NP, NT, chol, braw, ridge, &
        reg, rsf, csf, s, sipiv, info)
   call seed_x(x_true, ND*NT)
   do k = 1, ND
      x_true((k-1)*NT + NT-1) = 0d0
   end do
   v = 0d0
   do k = 1, ND
      do i = 1, NT
         v((k-1)*NT+i) = dot_product(a(i,1:NT,k), &
              x_true((k-1)*NT+1:k*NT))
      end do
   end do
   call jfnk_schur_apply_v1(ND, NT, NP, NT, NT-1, chol, braw, &
        ridge, reg, rsf, csf, s, sipiv, v, ierr)
   err = 0d0
   do k = 1, ND
      do i = 1, NP
         err = max(err, abs(v((k-1)*NT+i) - x_true((k-1)*NT+i)))
      end do
      err = max(err, abs(v((k-1)*NT+NT) - x_true((k-1)*NT+NT)))
   end do
   write(luer,'(A,I3,A,ES12.4)') &
      'S1b: apply ierr=', ierr, '  max|z-x|=', err
   if(info == 0 .and. ierr == 0 .and. err < 1d-6) then
      write(luer,'(A)') 'S1b: PASS (algebra exact)'
   else
      write(luer,'(A)') 'S1b: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ===== S2: rank-deficient Jnn under lsq_lam > 0 =====
   !
   do k = 1, ND
      call gen_block(a(1,1,k), NT, NP, k)
   end do
   a(2,:,2) = a(1,:,2)      ! depth 2: Jnn row 2 = row 1 -> rank NP-1
   reg = .false.
   lsq_lam = 1d-2
   call build_inputs(a, chol, braw, ridge, reg, rsf, csf, &
                     lsq_lam, 1d-2, lsq_ok, NT, NP, ND)
   write(luer,'(A,3L3)') 'S2: per-depth potrf ok=', lsq_ok
   call jfnk_schur_form_v1(ND, NT, NP, NT, chol, braw, ridge, &
        reg, rsf, csf, s, sipiv, info)
   call seed_x(x_true, ND*NT)
   v = 0d0
   do k = 1, ND
      do i = 1, NT
         v((k-1)*NT+i) = dot_product(a(i,1:NT,k), &
              x_true((k-1)*NT+1:k*NT))
      end do
   end do
   call jfnk_schur_apply_v1(ND, NT, NP, NT, NT-1, chol, braw, &
        ridge, reg, rsf, csf, s, sipiv, v, ierr)
   nbad = 0
   do i = 1, ND*NT
      if(.not. (abs(v(i)) <= huge(1d0))) nbad = nbad + 1
   end do
   err = 0d0
   do k = 1, ND
      do i = 1, NP
         err = max(err, abs(v((k-1)*NT+i) - x_true((k-1)*NT+i)))
      end do
   end do
   write(luer,'(A,I3,A,I3,A,ES12.4)') &
      'S2: apply ierr=', ierr, ' nonfinite=', nbad, &
      '  max|z-x| pops=', err
   if(ierr == 0 .and. nbad == 0) then
      write(luer,'(A)') 'S2: LSQ regularizes rank deficiency ' // &
         'silently (pseudo-inverse, no flag)'
   else
      write(luer,'(A)') 'S2: surfaced failure'
   end if
   !
   ! ===== S3: potrf-failed chol -> silent propagation =====
   ! Jnn zero COLUMN at depth 2: B^T B is exactly singular (rank
   ! NP-1), and the row-based reg cannot see it -> DPOTRF fails and
   ! chol keeps a partial factor. DPOTRS then returns info=0 while
   ! dividing by the failed pivot, so the defect propagates into S
   ! and the apply output without any ierr flag.
   do k = 1, ND
      call gen_block(a(1,1,k), NT, NP, k)
   end do
   a(1:NP,3,2) = 0d0
   reg = .false.
   lsq_lam = 0d0
   call build_inputs(a, chol, braw, ridge, reg, rsf, csf, &
                     lsq_lam, 0d0, lsq_ok, NT, NP, ND)
   write(luer,'(A,3L3)') 'S3: per-depth potrf ok=', lsq_ok
   call jfnk_schur_form_v1(ND, NT, NP, NT, chol, braw, ridge, &
        reg, rsf, csf, s, sipiv, info)
   v = 1d0
   call jfnk_schur_apply_v1(ND, NT, NP, NT, NT-1, chol, braw, &
        ridge, reg, rsf, csf, s, sipiv, v, ierr)
   nbad = 0
   do i = 1, ND*NT
      if(.not. (abs(v(i)) <= huge(1d0))) nbad = nbad + 1
   end do
   write(luer,'(A,I3,A,I3,A,I3,A,L3)') &
      'S3: form info=', info, ' apply ierr=', ierr, &
      ' nonfinite=', nbad, ' S finite=', all(abs(s) <= huge(1d0))
   write(luer,'(A)') 'S3: potrf-failed chol propagates through ' // &
      'dpotrs (info=0) -- schur_form ierr guard never fires'
   !
   ! ===== S4: reg row pass-through in blk_lsq =====
   !
   do k = 1, ND
      call gen_block(a(1,1,k), NT, NP, k)
   end do
   a(3,1:NP,1) = 0d0       ! depth 1 pop row 3 null -> regularized
   reg = .false.
   reg(3,1) = .true.
   lsq_lam = 1d-2
   call build_inputs(a, chol, braw, ridge, reg, rsf, csf, &
                     lsq_lam, 1d-2, lsq_ok, NT, NP, ND)
   ! blk_lsq on depth 1, rhs = e_3: reg row should pass z_3 ~ v_3
   block
      double precision :: w(NP), zz(NP)
      integer :: io
      w = 0d0
      w(3) = 1d0
      call jfnk_blk_lsq_solve(NP, chol(1:NP,1:NP,1), &
           braw(1:NP,1:NP,1), ridge(1:NP,1), reg(1:NP,1), &
           rsf(1:NP,1), csf(1:NP,1), w, zz, io)
      write(luer,'(A,I3,A,ES12.4)') &
         'S4: blk_lsq info=', io, ' z_3/v_3=', zz(3)
      if(io == 0 .and. abs(zz(3) - 1d0) < 0.2d0) then
         write(luer,'(A)') 'S4: PASS (reg row ~ identity pass-through)'
      else
         write(luer,'(A)') 'S4: FAIL'
         nfail = nfail + 1
      end if
   end block
   !
   if(nfail == 0) then
      print *, 'ALL SCHUR TESTS PASSED'
   else
      print *, nfail, ' SCHUR TESTS FAILED'
   end if
   call exit(nfail)
contains

   subroutine gen_block(m, nt, np, k)
      implicit none
      integer, intent(in) :: nt, np, k
      double precision, intent(out) :: m(nt,nt)
      integer :: i, j
      m = 0d0
      ! well-conditioned Jnn (diag dominant)
      do j = 1, np
         do i = 1, np
            m(i,j) = 0.05d0 * sin(1.3d0 * i + 0.7d0 * j + k) / np
         end do
      end do
      do i = 1, np
         m(i,i) = 1d0 + 0.2d0 * i
      end do
      ! Ne row/col (nt-1): small but nonzero couplings
      do i = 1, np
         m(nt-1,i) = 0.01d0 * sin(2.1d0 * i + k)
         m(i,nt-1) = 0.01d0 * cos(1.9d0 * i + k)
      end do
      m(nt-1,nt-1) = 1d0
      ! T row/column (nt): JTn, JnT, JTT
      do i = 1, np
         m(nt,i) = 0.03d0 * cos(1.1d0 * i + 0.5d0 * k)
         m(i,nt) = 0.03d0 * sin(0.9d0 * i + 0.3d0 * k)
      end do
      m(nt-1,nt) = 0.01d0 * sin(k + 0.4d0)
      m(nt,nt-1) = 0.01d0 * cos(k + 0.8d0)
      m(nt,nt) = 0.8d0 + 0.1d0 * k
   end subroutine gen_block

   subroutine seed_x(x, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(out) :: x(n)
      integer :: i
      do i = 1, n
         x(i) = sin(0.37d0 * i) + 0.3d0 * cos(0.11d0 * i)
      end do
   end subroutine seed_x

   ! Replicates the production per-depth chain (jfnk_cmfgen_interface_
   ! v1.f ~4590-4670): B = Jnn + ridge*diag with reg rows -> e_i,
   ! G = B^T B + (lam*max(blkmax,1))^2 I, chol = potrf(G).
   subroutine build_inputs(m, ch, br, rg, rg_l, rs, cs, lam, &
                           ridge_in, ok, nt, np, nd)
      implicit none
      integer, intent(in) :: nt, np, nd
      double precision, intent(in) :: m(nt,nt,nd), lam, ridge_in
      double precision, intent(out) :: ch(nt,nt,nd), br(nt,nt,nd), &
           rg(nt,nd), rs(nt,nd), cs(nt,nd)
      logical, intent(in) :: rg_l(nt,nd)
      logical, intent(out) :: ok(nd)
      double precision :: b(np,np), g(np,np), blkmax, ridge_add, l2
      integer :: k, i, j, io
      ch = 0d0
      br = m
      rs = 1d0
      cs = 1d0
      rg = 0d0
      ridge_add = ridge_in          ! 1e-2 = production ITRY=0 ridge
      do k = 1, nd
         blkmax = maxval(abs(m(1:np,1:np,k)))
         b = m(1:np,1:np,k)
         do i = 1, np
            rg(i,k) = ridge_add * max(abs(m(i,i,k)), 1d-30)
            b(i,i) = b(i,i) + rg(i,k)
            if(rg_l(i,k)) then
               b(i,:) = 0d0
               b(i,i) = 1d0
            end if
         end do
         l2 = (lam * max(blkmax, 1d0))**2
         g = matmul(transpose(b), b)
         do i = 1, np
            g(i,i) = g(i,i) + l2
         end do
         ch(1:np,1:np,k) = g
         call dpotrf('L', np, ch(1,1,k), nt, io)
         ok(k) = (io == 0)
      end do
   end subroutine build_inputs

end program test_schur
