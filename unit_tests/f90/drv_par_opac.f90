    program drv_par_opac
    use set_kind_module
    use mod_par_opac
    implicit none

    integer, parameter :: nd = 4
    character(len=32) :: case_name
    integer :: i

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    allocate (r(nd), v(nd))
    r = [(20.0_ldp - 2.0_ldp*(i - 1), i = 1, nd)]
!   All species absent -> the equation-labelling loop is a no-op,
!   leaving I=1 so EQNE=1 and NT=I+1=2 satisfy the pre-checks.
    species_pres = .false.
    eqne = 1

    select case (case_name)
    case ('flat')
        v = [1000.0_ldp, 800.0_ldp, 800.0_ldp, 500.0_ldp]
    case ('reversal')
        v = [1000.0_ldp, 800.0_ldp, 1200.0_ldp, 500.0_ldp]
    case ('neg')
        v = [1000.0_ldp, 800.0_ldp, -100.0_ldp, 500.0_ldp]
    case ('good')
        v = [1000.0_ldp, 800.0_ldp, 600.0_ldp, 500.0_ldp]
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    call par_opac_sub(nd, 2, 6, nd, 6, 2, 10)
    write (*, '(a)') 'DONE_' // case_name

    end program drv_par_opac
