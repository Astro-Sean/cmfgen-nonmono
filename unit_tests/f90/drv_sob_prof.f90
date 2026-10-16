! Driver for ADD_SOB_LINE_PROF (embedded at the end of
! obs/cmf_flux_sub_v5.f).
!
! Usage: drv_sob_prof <scenario>
!   mono     : strictly decreasing NU, one SOB line inside the grid
!   offgrid  : SOB line rest frequency outside [NU(NCF),NU(1)]
!   maser    : CHIL_MAT<=0 at 60% of depths -> per-line drop warning
!   nonmono  : NU has one non-decreasing pair -> assertion STOP
!
! Arrays are built procedurally; ND=12, NCF=200, one blend group of
! NUM_SIM=2 lines (one SOB, one CMF) drawn from N_LINES=2.
!
program drv_sob_prof
  use set_kind_module
  implicit none

  character(len=64) :: scen
  integer, parameter :: nd = 12, ncf = 200, max_sim = 2, n_lines = 2
  integer, parameter :: num_sim = 2
  real(kind=ldp) :: nu(ncf), chil(nd,max_sim), etal(nd,max_sim)
  real(kind=ldp) :: clump(nd), fl_sim(max_sim), vdop_min(n_lines)
  real(kind=ldp) :: eta_st(nd,ncf), chi_st(nd,ncf)
  real(kind=ldp) :: c_kms, drop_sum, inj_sum, integ, nu_prev
  integer :: sim_ptr(max_sim), n_drop
  logical :: neg_opacity(nd)
  character(len=10) :: trans_type(n_lines)
  integer :: i, k, idep

  call get_command_argument(1, scen)

  ! NU strictly decreasing (units of 1e15 Hz): coarse 100->14, a
  ! fine segment 14->8 around the test line at F0=10, coarse 8->0.5.
  ! The fine segment must resolve the Doppler core (W=3*DNU_D) or
  ! XNORM=0 and the line is skipped by design.
  nu(1:60) = 100.0_ldp + (14.0_ldp-100.0_ldp)* &
       (/(real(k-1,ldp),k=1,60)/)/59.0_ldp
  nu(61:160) = 14.0_ldp + (8.0_ldp-14.0_ldp)* &
       (/(real(k-1,ldp),k=1,100)/)/99.0_ldp
  nu(161:200) = 8.0_ldp + (0.5_ldp-8.0_ldp)* &
       (/(real(k-1,ldp),k=1,40)/)/39.0_ldp
  ! Guard: keep the grid strictly decreasing across the seams.
  nu(60) = nu(59) - 0.5_ldp*(nu(59)-nu(61))
  nu(160) = nu(159) - 0.5_ldp*(nu(159)-nu(161))
  if (scen == 'nonmono') nu(100) = nu(99)   ! break strict decrease

  chil = 1.0e-3_ldp
  etal = 1.0e-6_ldp
  clump = 1.0_ldp
  fl_sim(1) = 10.0_ldp
  fl_sim(2) = 20.0_ldp
  vdop_min = 5000.0_ldp   !ejecta-scale width so W spans several cells
  c_kms = 2.99792458e5_ldp
  sim_ptr(1) = 1
  sim_ptr(2) = 2
  trans_type(1) = 'SOB'
  trans_type(2) = 'CMF'
  neg_opacity = .false.
  eta_st = 0.0_ldp
  chi_st = 0.0_ldp
  n_drop = 0
  drop_sum = 0.0_ldp
  inj_sum = 0.0_ldp

  if (scen == 'offgrid') fl_sim(1) = 200.0_ldp   ! above NU(1)=100
  if (scen == 'maser') chil(1:8,1) = -1.0e-3_ldp ! 8/12 depths masing

  call add_sob_line_prof(eta_st, chi_st, nu, ncf, chil, etal, clump, &
       fl_sim, sim_ptr, vdop_min, trans_type, neg_opacity, &
       num_sim, n_lines, max_sim, nd, c_kms, n_drop, drop_sum, inj_sum)

  write(*,'(a,2(1x,i6))') 'N_DROP', n_drop, nd
  write(*,'(a,2(1x,es24.16e3))') 'DROP_INJ ', drop_sum, inj_sum

  ! Frequency integral of the stored CHI_ST column at the first
  ! non-masing depth.  With the cell-width normalization it should
  ! equal CHIL*CLUMP*1e-15 (the per-Hz integral over the 1e15-Hz grid
  ! carries the explicit 1e-15 factor).
  idep = 1
  do i = 1, nd
     if (chil(i,1) > 0.0_ldp .and. .not. neg_opacity(i)) then
        idep = i
        exit
     end if
  end do
  integ = 0.0_ldp
  do k = 1, ncf
     if (k == 1) then
        nu_prev = nu(1) - nu(2)
     else if (k == ncf) then
        nu_prev = nu(ncf-1) - nu(ncf)
     else
        nu_prev = 0.5_ldp*(nu(k-1)-nu(k+1))
     end if
     integ = integ + chi_st(idep,k)*nu_prev
  end do
  write(*,'(a,2(1x,es24.16e3))') 'CHI_INTEG_EXPECT ', integ, &
       chil(idep,1)*clump(idep)*1.0e-15_ldp
  write(*,'(a,es24.16e3)') 'ETA_SUM ', sum(eta_st)
  write(*,'(a,es24.16e3)') 'CHI_SUM ', sum(chi_st)

end program drv_sob_prof
