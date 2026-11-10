    program drv_check_ions
    use set_kind_module
    use mod_cmfgen
    implicit none

    integer, parameter :: nd_in = 4
    character(len=32) :: case_name

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    num_ions = 1
    species_pres = .false.
    at_no = 0.0_ldp
    species_beg_id = 0
    species_end_id = 0
    allocate (t(nd_in), ed(nd_in))
    ed = [1.0_ldp, 2.0_ldp, 3.0_ldp, 4.0_ldp]

    select case (case_name)
    case ('zero_t')
        t = 0.0_ldp
    case ('pos_t')
        t = [1.0_ldp, 3.0_ldp, 2.0_ldp, 0.5_ldp]
    case ('neg_t')
        t = -1.0_ldp
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    call check_ions_present(nd_in)
    write (*, '(a)') 'DONE_' // case_name

    end program drv_check_ions
