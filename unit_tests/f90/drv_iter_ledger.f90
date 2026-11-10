! drv_iter_ledger -- CLI driver for the T-06 ITER_LEDGER module
! (audit F-22/F-24). Runs in a scratch dir, exercises the reset /
! counter / timer / CSV-write path, and prints markers the Python
! test parses. Arguments:
!
!   argv1: 'on'  -> arm LEDG_ITER_RESET then emit two rows
!          'off' -> call LEDG_WRITE_ITER without arming; the internal
!                   LEDG_ACTIVE guard must leave no file behind (the
!                   stock-path acceptance: zero ledger output)
!
! Rows are written with distinct caller args so the test can verify
! every column against the values passed in.
!
program drv_iter_ledger
  use set_kind_module
  use iter_ledger_mod
  implicit none
!
  character(len=16) mode
  real(kind=LDP) wt1,wt2,rss
  real(kind=LDP) omp_get_wtime
  external omp_get_wtime
!
  call getarg(1,mode)
!
  if (mode .eq. 'on') then
    call ledg_iter_reset(100.0_LDP, 40, 12)
    ledg_n_fg   = 1300
    ledg_n_mom  = 1300
    ledg_n_rf   = 1
    ledg_np_rf  = 7
    ledg_n_frz  = 1
    ledg_t_rad  = 55.25_LDP
    ledg_t_rf   = 12.5_LDP
    ledg_t_p0   = 200.0_LDP
    wt1 = 300.0_LDP
    call ledg_write_iter(wt1, &
         12,'LAMBDA_TFIX','FULL','LAMBDA', &
         'gs_immature_demote+hold_pops', &
         .true.,.true., 1,0, .false.,.false., 52,24, 0, &
         3,2,1_8, &
         4.5E2_LDP,1.1E2_LDP,9.9E1_LDP, &
         3.2E4_LDP,1.7E3_LDP,0.5_LDP, &
         1,1,0,4,0, &
         3_8,1_8,0_8,2_8, &
         2_8,1_8,0_8,0,1)
    wt2 = 360.0_LDP
    call ledg_write_iter(wt2, &
         13,'FULL','FULL','HELD', &
         'untracked', &
         .false.,.false., 0,1, .true.,.true., 60,30, 1, &
         0,4,2_8, &
         3.0E-1_LDP,2.5E-1_LDP,3.0E-1_LDP, &
         8.0E1_LDP,2.0E1_LDP,1.0_LDP, &
         1,0,0,3,1, &
         0_8,0_8,0_8,0_8, &
         0_8,0_8,0_8,1,0)
    write(*,'(A)')'ROWS_WRITTEN=2'
    call ledg_rss_mb(rss)
    write(*,'(A,F12.1)')'RSS_MB=',rss
  else
    wt1 = omp_get_wtime()
! unarmed: LEDG_WRITE_ITER must return without creating ITER_LEDGER
    call ledg_write_iter(wt1, &
         1,'FULL','FULL','HELD', &
         ' ', &
         .false.,.false., 0,0, .false.,.false., 0,0, 0, &
         0,0,0_8, &
         1.0_LDP,1.0_LDP,1.0_LDP, &
         1.0_LDP,1.0_LDP,1.0_LDP, &
         0,0,0,0,0, &
         0_8,0_8,0_8,0_8, &
         0_8,0_8,0_8,0,0)
    write(*,'(A)')'DONE'
  end if
!
end program drv_iter_ledger
