! Driver for the ilf_kvb / ilf_v* block-format tests: prints each
! edge case behind a 'CASE <name>' marker so the python side can
! split stdout per scenario.  Only formatting is exercised; no
! CMFGEN state is needed.
!
      program drv_ilf_kvb
      use iter_log_fmt_mod
      use set_kind_module
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan,        &
     &    ieee_positive_inf, ieee_negative_inf
      implicit none

      integer, parameter :: ml = 60
      character(len=ml) :: lab(9)
      character(len=ilf_val_len) :: val(9)
      character(len=100) :: wlab(3)
      character(len=ilf_val_len) :: wval(3)
      real(kind=ldp) :: nanv, pinf, ninf

      nanv = ieee_value(0.0_ldp, ieee_quiet_nan)
      pinf = ieee_value(0.0_ldp, ieee_positive_inf)
      ninf = ieee_value(0.0_ldp, ieee_negative_inf)

! The spec example: pad set by the two 17-char labels, mixed
! integer and ES11.2E2 fields.
      write(6,'(a)') 'CASE spec'
      lab(1) = 'freqs unconverged'
      lab(2) = 'max final dJ/J'
      lab(3) = 'ALO fallbacks'
      lab(4) = 'J-floor hits run'
      lab(5) = 'BA-floor hits run'
      val(1) = ilf_vi(98)
      val(2) = ilf_vr(1.96e-3_ldp, 2)
      val(3) = ilf_vi(0)
      val(4) = ilf_vi(0)
      val(5) = ilf_vi(0)
      call ilf_kvb(6, ' ES_ALI:', 5, lab, val)

! Real edges: 1E+99 and 1E-99 stay on two-digit exponents, zero and
! negatives print normally, |x|>=1E100 and |x|<1E-99 widen to E3,
! NaN and +/-Inf print as text.  No case may emit '*****'.
      write(6,'(a)') 'CASE edges'
      lab(1) = 'one e plus 99'
      lab(2) = 'one e minus 99'
      lab(3) = 'exactly zero'
      lab(4) = 'negative'
      lab(5) = 'exponent over 99'
      lab(6) = 'exponent under -99'
      lab(7) = 'a NaN'
      lab(8) = 'plus infinity'
      lab(9) = 'minus infinity'
      val(1) = ilf_vr(1.0e99_ldp, 2)
      val(2) = ilf_vr(1.0e-99_ldp, 2)
      val(3) = ilf_vr(0.0_ldp, 2)
      val(4) = ilf_vr(-4.5e3_ldp, 3)
      val(5) = ilf_vr(1.23e100_ldp, 2)
      val(6) = ilf_vr(1.2e-101_ldp, 2)
      val(7) = ilf_vr(nanv, 2)
      val(8) = ilf_vr(pinf, 2)
      val(9) = ilf_vr(ninf, 2)
      call ilf_kvb(6, 'EDGES:', 9, lab, val)

! A 10-digit integer fills I10 exactly; logicals take L1-style
! right alignment; the string field keeps its text.
      write(6,'(a)') 'CASE intlog'
      lab(1) = 'big counter'
      lab(2) = 'neg counter'
      lab(3) = 'flag on'
      lab(4) = 'flag off'
      lab(5) = 'status'
      val(1) = ilf_vi(1234567890)
      val(2) = ilf_vi(-42)
      val(3) = ilf_vl(.true.)
      val(4) = ilf_vl(.false.)
      val(5) = ilf_vs('(BA computed)')
      call ilf_kvb(6, 'INTLOG:', 5, lab, val)

! Empty field list: the header prints alone.
      write(6,'(a)') 'CASE empty'
      call ilf_kvb(6, 'EMPTY:', 0, lab, val)

! Labels past the 30-char ilf_kv width align at their own max.
      write(6,'(a)') 'CASE longlab'
      lab(1) = 'short'
      lab(2) = 'a label that is well past thirty characters'
      val(1) = ilf_vi(7)
      val(2) = ilf_vi(8)
      call ilf_kvb(6, 'LONGLAB:', 2, lab, val)

! A label too long for the shared '=' column wraps its value to a
! continuation line; the shorter labels still align under the
! shrunken padding cap and every line stays <= 100 columns.
      write(6,'(a)') 'CASE wide'
      wlab(1) = 'label a'
      wlab(2) = 'label b'
      wlab(3) = 'a very long label '//repeat('x', 70)//' end'
      wval(1) = ilf_vi(1)
      wval(2) = ilf_vi(2)
      wval(3) = ilf_vi(3)
      call ilf_kvb(6, 'WIDE:', 3, wlab, wval)

! LOG_COMPACT_KV collapses the block to the legacy one-liner.
      write(6,'(a)') 'CASE compact'
      log_compact_kv = .true.
      lab(1) = 'freqs unconverged'
      lab(2) = 'max final dJ/J'
      lab(3) = 'ALO fallbacks'
      lab(4) = 'J-floor hits run'
      lab(5) = 'BA-floor hits run'
      val(1) = ilf_vi(98)
      val(2) = ilf_vr(1.96e-3_ldp, 2)
      val(3) = ilf_vi(0)
      val(4) = ilf_vi(0)
      val(5) = ilf_vi(0)
      call ilf_kvb(6, ' ES_ALI:', 5, lab, val)

      end program drv_ilf_kvb
