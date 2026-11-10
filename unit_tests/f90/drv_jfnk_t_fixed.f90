program drv_jfnk_t_fixed
  implicit none
  logical, external :: jfnk_t_fixed_mode_v1
  integer :: failures

  failures = 0
  call check(.false., .true., .false., .true., &
      "VADAT FIX_T overrides free JFNK mode")
  call check(.false., .false., .true., .true., &
      "JFNK_FIXT holds temperature for Newton")
  call check(.true., .false., .false., .true., &
      "transient fixed-T mode holds temperature")
  call check(.false., .false., .false., .false., &
      "variable-T mode stays free")
  if (failures > 0) stop 1

contains

  subroutine check(fixed_t_mode, rd_fix_t, jfnk_fixt, expected, label)
    implicit none
    logical, intent(in) :: fixed_t_mode, rd_fix_t, jfnk_fixt, expected
    character(len=*), intent(in) :: label
    logical :: actual

    actual = jfnk_t_fixed_mode_v1(fixed_t_mode, rd_fix_t, jfnk_fixt)
    if (actual .eqv. expected) then
      write(*, '(A)') "  [PASS] " // label
    else
      write(*, '(A)') "  [FAIL] " // label
      failures = failures + 1
    end if
  end subroutine check

end program drv_jfnk_t_fixed
