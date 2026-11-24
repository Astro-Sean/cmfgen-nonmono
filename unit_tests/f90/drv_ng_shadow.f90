! Driver for NG_SHADOW_DIAG tests: loads a synthetic SCRTEMP
! history through the serving stub and calls the production
! DO_NG_BAND_ACCEL_V4, with the diagnostic flag set from the
! command line. Prints NG_DONE and the resulting populations
! so the test can verify the diagnostics change nothing.
!
      program drv_ng_shadow
      use control_variable_mod, only: ng_shadow_diag,             &
     &                              do_nonmono_sob,               &
     &                              do_nonmono_knop
      use stub_ng_mod
      use set_kind_module
      implicit none

      integer nt, nd, nband
      real(kind=ldp), allocatable :: pops(:,:)
      real(kind=ldp), allocatable :: r(:), v(:), sigma(:), rold(:)
      real(kind=ldp) maxdec, maxinc
      logical ng_done
      character(len=512) cfile, flagarg
      integer k, l

      call getarg(1, cfile)
      call getarg(2, flagarg)
      call stub_ng_load(trim(cfile))

      nt = nt_s
      nd = nd_s
      allocate(pops(nt, nd), r(nd), v(nd), sigma(nd), rold(nd))

      ! Non-mono flags stay F: the freeze path is exercised by
      ! the production runs, and keeping it off lets the
      ! structure rows carry record-to-record variation.
      do_nonmono_sob = .false.
      do_nonmono_knop = .false.
      ng_shadow_diag = (flagarg(1:1) .eq. 'T')

      r = rs
      v = vs
      sigma = ss
      rold = rs
      pops = recs(:, :, nrec_s)

      nband = 4
      call do_ng_band_accel_v4(pops, r, v, sigma, rold,           &
     &    nt, nd, nband, ng_done, maxdec, maxinc,                 &
     &    .false., .true., 6, 6)

      write(6, '(a,l2)') 'NG_DONE=', ng_done
      write(6, '(a,2es14.6)') 'MAXINC/DEC=', maxinc, maxdec
      do l = 1, nd
        do k = 1, nt
          write(6, '(a,i4,i4,es24.16)') 'POP', k, l, pops(k, l)
        end do
      end do
      end program drv_ng_shadow
