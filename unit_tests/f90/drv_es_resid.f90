! Driver for the P8.2 ES redistribution convergence diagnostic.
! Manufactures field pairs and calls the production routines in
! obs/es_resid_diag.f (linked as the production .o, never a copy).
!
!   drv_es_resid compare <case> [q]
!       Runs ES_RESID_COMPARE on the manufactured pair and prints
!       R_GLOB/K_GLOB/RREL_GLOB/D_GLOB/TOT_*/VALID plus per-depth
!       rows.  Cases:
!         A  identical fields          R=0
!         B  NEW = 1.01*OLD            R ~ 0.01
!         C  single-node perturbation  R at (k0,m0)
!         D  zero fields               R=0
!         E  zero -> nonzero           R=1
!         F  weak-node perturbation    node masked
!         G1 NaN in NEW                VALID=F, BAD_WHICH=1
!         G2 +Inf in OLD               VALID=F, BAD_WHICH=-1
!         H  negative values           counted, not clamped
!         I  contraction q             R_k ~ q^k (q arg, default .8)
!         J  oscillating               R stays O(1)
!   drv_es_resid seq
!       Full production chain on manufactured ES_J_CONV files:
!       seed file -> SEED -> update file -> UPDATE -> update file
!       -> UPDATE -> RENDER(latest) -> RENDER(stale) -> FINAL.
!       Exercises the generation bookkeeping and the render
!       consumption check end to end (TEST L).
!   drv_es_resid omp
!       Same manufactured pair evaluated once; the Python side runs
!       this under OMP_NUM_THREADS=1 and =4 and diffs output (TEST K).
!
      program drv_es_resid
      use set_kind_module
      use es_resid_diag_mod
      use, intrinsic :: ieee_arithmetic
      implicit none
!
      integer, parameter :: nd = 8, ncf = 12
      real(kind=ldp) :: rj_new(nd, ncf), rj_old(nd, ncf)
      real(kind=ldp) :: r_k(nd), d_k(nd), s_k(nd), rrel_k(nd)
      integer :: nact_k(nd), nmask_k(nd), nnonf_k(nd), nneg_k(nd)
      real(kind=ldp) :: afloor, r_glob, rrel_glob, d_glob, q
      integer :: k_glob, tot_act, tot_inval, tot_neg
      integer :: bad_k, bad_m, bad_which
      logical :: valid
      character(len=32) :: mode, cname, carg
      integer :: k, m, ngen
!
      call get_command_argument(1, mode)
      if (len_trim(mode) == 0) then
        write(6, *) 'usage: drv_es_resid compare <case> [q] | seq | omp'
        stop 1
      end if
!
      afloor = 1.0e-12_ldp
      rj_old = 0.0_ldp
      rj_new = 0.0_ldp
      rj_old = base_field()
      rj_new = rj_old
!
      if (trim(mode) == 'compare') then
        call get_command_argument(2, cname)
        call get_command_argument(3, carg)
        q = 0.8_ldp
        if (len_trim(carg) > 0) read(carg, *) q
        select case (trim(cname))
        case ('A')
        case ('B')
          rj_new = 1.01_ldp * rj_old
        case ('C')
          rj_new(4, 7) = rj_old(4, 7) + 3.0e-6_ldp
        case ('D')
          rj_old = 0.0_ldp
          rj_new = 0.0_ldp
        case ('E')
          rj_old = 0.0_ldp
        case ('F')
! The weak node must sit below AFLOOR*S_K - 12+ orders under the
! depth scale - to land on the masked side of the active mask.
          rj_old(:, :) = base_field()
          rj_old(2, 5) = 1.0e-25_ldp
          rj_new(2, 5) = 3.0e-25_ldp
        case ('G1')
          rj_new(3, 6) = ieee_value(0.0_ldp, ieee_quiet_nan)
        case ('G2')
          rj_old(5, 9) = ieee_value(0.0_ldp, ieee_positive_inf)
        case ('H')
          rj_old(2, 3) = -rj_old(2, 3)
          rj_new(6, 10) = -rj_new(6, 10)
        case ('I', 'J')
          call seq_series(q, cname == 'J')
          stop
        case default
          write(6, *) 'unknown case ', trim(cname)
          stop 1
        end select
!
        bad_k = 0
        bad_m = 0
        bad_which = 0
        call es_resid_compare(rj_new, rj_old, nd, ncf, afloor,       &
          r_k, d_k, s_k, rrel_k, nact_k, nmask_k, nnonf_k, nneg_k,  &
          r_glob, k_glob, rrel_glob, d_glob, tot_act, tot_inval,    &
          tot_neg, bad_k, bad_m, bad_which, valid)
        call print_result('SINGLE', r_glob, k_glob, rrel_glob,       &
          d_glob, tot_act, tot_inval, tot_neg, bad_k, bad_m,         &
          bad_which, valid)
!
      else if (trim(mode) == 'seq') then
        call run_seq()
      else if (trim(mode) == 'seqstale') then
        call run_seq_stale()
      else if (trim(mode) == 'seqcold') then
        call run_seq_cold()
      else if (trim(mode) == 'omp') then
        bad_k = 0
        bad_m = 0
        bad_which = 0
        rj_new = 1.01_ldp * rj_old
        call es_resid_compare(rj_new, rj_old, nd, ncf, afloor,       &
          r_k, d_k, s_k, rrel_k, nact_k, nmask_k, nnonf_k, nneg_k,  &
          r_glob, k_glob, rrel_glob, d_glob, tot_act, tot_inval,    &
          tot_neg, bad_k, bad_m, bad_which, valid)
        call print_result('OMP', r_glob, k_glob, rrel_glob,          &
          d_glob, tot_act, tot_inval, tot_neg, bad_k, bad_m,         &
          bad_which, valid)
      else
        write(6, *) 'unknown mode ', trim(mode)
        stop 1
      end if
!
      contains
!
!
!
      function base_field() result(f)
      real(kind=ldp) :: f(nd, ncf)
      integer :: kk, mm
      do kk = 1, nd
        do mm = 1, ncf
          f(kk, mm) = 1.0e-8_ldp * kk * exp(-0.3_ldp * mm)
        end do
      end do
      end function base_field
!
!
!
      subroutine print_result(tag, rg, kg, rrg, dg, ta, ti, tn,     &
                              bk, bm, bw, ok)
      character(len=*) :: tag
      real(kind=ldp) :: rg, rrg, dg
      integer :: kg, ta, ti, tn, bk, bm, bw
      logical :: ok
      integer :: kk
      write(6, '(A)') 'RESULT ' // trim(tag)
      write(6, '(A,ES16.8)') 'R_GLOB=    ', rg
      write(6, '(A,I0)') 'K_GLOB=', kg
      write(6, '(A,ES16.8)') 'RREL_GLOB= ', rrg
      write(6, '(A,ES16.8)') 'D_GLOB=    ', dg
      write(6, '(A,I0)') 'TOT_ACT=', ta
      write(6, '(A,I0)') 'TOT_INVAL=', ti
      write(6, '(A,I0)') 'TOT_NEG=', tn
      write(6, '(A,I0,A,I0,A,I0)') 'BAD=', bk, ',', bm, ',', bw
      write(6, '(A,L1)') 'VALID=', ok
      do kk = 1, nd
        write(6, '(A,I0,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,' //    &
          'A,I0,A,I0,A,I0,A,I0)')                                  &
          'DROW K=', kk, ' R=', r_k(kk), ' RREL=', rrel_k(kk),     &
          ' D=', d_k(kk), ' S=', s_k(kk), ' NACT=', nact_k(kk),    &
          ' NMASK=', nmask_k(kk), ' NNONFIN=', nnonf_k(kk),        &
          ' NNEG=', nneg_k(kk)
      end do
      end subroutine print_result
!
!
!
      subroutine seq_series(qv, osc)
! Manufacture a generation sequence and compare consecutive pairs:
! contraction J_k = J_fix + q^k*D (R should fall by ~q each step),
! oscillation J_k = J_fix + (-1)^k*D (R stays large).
      real(kind=ldp), intent(in) :: qv
      logical, intent(in) :: osc
      real(kind=ldp) :: fj(nd, ncf), delta(nd, ncf), sgn
      real(kind=ldp) :: rg, rrg, dg
      integer :: kg, ta, ti, tn, bk, bm, bw, kk, n
      logical :: ok
! Contraction needs delta small against the field scale, else the
! depth normalizer S_K shrinks with the decaying delta and the
! measured ratio only asymptotes to q.  Oscillation wants the
! opposite so the residual stays large.
      if (osc) then
        delta = 1.0e-7_ldp
      else
        delta = 1.0e-10_ldp
      end if
      fj = base_field()
! Generation 0 is identical to its predecessor so the first reported
! residual is 0 and the geometric decay shows up from gen 1 on.
      rj_old = fj + delta
      do n = 0, 5
        if (osc) then
          sgn = merge(-1.0_ldp, 1.0_ldp, mod(n, 2) == 1)
          rj_new = fj + sgn * delta
        else
          rj_new = fj + qv**n * delta
        end if
        bk = 0
        bm = 0
        bw = 0
        call es_resid_compare(rj_new, rj_old, nd, ncf, afloor,      &
          r_k, d_k, s_k, rrel_k, nact_k, nmask_k, nnonf_k, nneg_k, &
          rg, kg, rrg, dg, ta, ti, tn, bk, bm, bw, ok)
        write(6, '(A,I0,A,ES16.8)') 'SERIES GEN=', n, ' R=', rg
        rj_old = rj_new
      end do
      end subroutine seq_series
!
!
!
      subroutine write_esj(fname, field, fnd, fncf)
! Write a manufactured ES_J_CONV in the exact EDDFACTOR layout
! COMP_J_CONV_V2 produces: sidecar INFO with ND/RECL, record 3
! header (INIT_REC, NCF, ND), data records INIT_REC+ML-1.
      character(len=*) :: fname
      integer :: fnd, fncf
      real(kind=ldp) :: field(fnd, fncf)
      integer :: recl_lim, usize, wsize, maxrec, irecl, lu, ml, kk
      real(kind=ldp) :: nuv
      call dir_acc_pars(recl_lim, usize, wsize, maxrec)
      irecl = wsize * (fnd + 1) / usize
      call get_lu(lu, 'drv write_esj')
      call write_direct_info_v3(fnd, irecl, '15-Nov-2026', fname, lu)
      open(unit=lu, file=fname, form='unformatted', access='direct', &
           status='replace', recl=irecl)
      write(lu, rec=3) 15, fncf, fnd
      do ml = 1, fncf
        nuv = 200.0_ldp - ml
        write(lu, rec=15 + ml - 1) (field(kk, ml), kk = 1, fnd), nuv
      end do
      close(lu)
      end subroutine write_esj
!
!
!
      subroutine run_seq()
! Seed with generation-0 file, two update generations, then the
! render consumption check against both the latest and a stale
! field (TEST L: the check must distinguish them).
      real(kind=ldp) :: g0(nd, ncf), g1(nd, ncf), g2(nd, ncf)
      real(kind=ldp) :: render(nd, ncf)
      g0 = base_field() * 1.10_ldp
      g1 = base_field() * 1.05_ldp
      g2 = base_field() * 1.01_ldp
      call write_esj('ES_J_CONV', g0, nd, ncf)
      call es_resid_seed(nd, ncf, .true., afloor, .true., 3, 6)
      call write_esj('ES_J_CONV', g1, nd, ncf)
      call es_resid_update(nd, ncf, 1, 3, 6)
      call write_esj('ES_J_CONV', g2, nd, ncf)
      call es_resid_update(nd, ncf, 2, 3, 6)
!
      render = g2
      call es_resid_render(render, nd, ncf, 6)
      call es_resid_final(3, 6)
      write(6, '(A)') 'SEQ_DONE'
      end subroutine run_seq
!
!
!
      subroutine run_seq_cold()
! Cold start: no ES_J_CONV at SEED time, so the first UPDATE must
! report a first generation rather than inventing a comparison
! (spec case 13).
      real(kind=ldp) :: g1(nd, ncf)
      g1 = base_field() * 1.05_ldp
      call es_resid_seed(nd, ncf, .false., afloor, .true., 2, 6)
      call write_esj('ES_J_CONV', g1, nd, ncf)
      call es_resid_update(nd, ncf, 1, 2, 6)
      call es_resid_final(2, 6)
      write(6, '(A)') 'SEQ_COLD_DONE'
      end subroutine run_seq_cold
!
!
!
      subroutine run_seq_stale()
      real(kind=ldp) :: g1(nd, ncf), g2(nd, ncf), render(nd, ncf)
      g1 = base_field() * 1.05_ldp
      g2 = base_field() * 1.01_ldp
      call write_esj('ES_J_CONV', g1, nd, ncf)
      call es_resid_seed(nd, ncf, .false., afloor, .true., 2, 6)
      call write_esj('ES_J_CONV', g2, nd, ncf)
      call es_resid_update(nd, ncf, 1, 2, 6)
      render = g1
      call es_resid_render(render, nd, ncf, 6)
      call es_resid_final(2, 6)
      write(6, '(A)') 'SEQ_STALE_DONE'
      end subroutine run_seq_stale
!
!
!
      end program drv_es_resid
