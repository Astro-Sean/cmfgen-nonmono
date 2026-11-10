    program drv_vturb
    use set_kind_module
    use control_variable_mod
    use mod_cmfgen
    implicit none

    integer, parameter :: nd_in = 5
    real(kind=ldp) :: maw, asum
    character(len=32) :: case_name
    integer :: i

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    do_clump_model = .false.
    sn_hydro_model = .false.
    sn_model = .false.
    var_mdot = .false.
    veltype = 1
    rmdot = 1.0e-5_ldp
    do_nonmono_sob = .false.
    do_nonmono_knop = .false.

    vturb_min = 5.0_ldp
    vturb_max = 25.0_ldp
    vturb_pow = 0.5_ldp
    vturb_vend = 15.0_ldp
    vturb_law = 'POW'

    at_abund = 0.0_ldp
    at_abund(1) = 1.0_ldp
    at_mass = 1.0_ldp

    allocate (r(nd_in), v(nd_in), clump_fac(nd_in))
    allocate (pop_atom(nd_in), pop_species(nd_in, num_species))
    allocate (density(nd_in), vturb_vec(nd_in))
    r = [12.0_ldp, 10.0_ldp, 8.0_ldp, 6.0_ldp, 3.0_ldp]
    clump_fac = 1.0_ldp

    select case (case_name)
    case ('pow_nm')
        do_nonmono_sob = .true.
        v = [40.0_ldp, 10.0_ldp, -5.0_ldp, 20.0_ldp, 30.0_ldp]
    case ('pow_nm_knop')
        do_nonmono_knop = .true.
        v = [40.0_ldp, 10.0_ldp, -5.0_ldp, 20.0_ldp, 30.0_ldp]
    case ('pow_off')
        v = [40.0_ldp, 30.0_ldp, 20.0_ldp, 15.0_ldp, 10.0_ldp]
    case ('pow_off_nm')
        v = [40.0_ldp, 10.0_ldp, -5.0_ldp, 20.0_ldp, 30.0_ldp]
    case ('pow_zero')
        do_nonmono_sob = .true.
        v = 0.0_ldp
    case ('tpow_nm')
        do_nonmono_sob = .true.
        vturb_law = 'TRUNC_POW'
        v = [40.0_ldp, 10.0_ldp, -5.0_ldp, 20.0_ldp, 30.0_ldp]
    case ('tpow_off')
        vturb_law = 'TRUNC_POW'
        v = [40.0_ldp, 30.0_ldp, 20.0_ldp, 15.0_ldp, 10.0_ldp]
    case ('bad_law')
        vturb_law = 'BOGUS'
        v = [40.0_ldp, 30.0_ldp, 20.0_ldp, 15.0_ldp, 10.0_ldp]
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    call set_abund_clump(maw, asum, 6, nd_in)

    write (*, '(a)') 'DONE_' // case_name
    write (*, '(a,5es14.6)') 'VT=', vturb_vec
    write (*, '(a,l1)') 'VT_FIN=', &
        all(vturb_vec == vturb_vec) .and. &
        all(vturb_vec < huge(1.0_ldp))
    write (*, '(a,l1)') 'VT_NAN=', any(vturb_vec /= vturb_vec)
    write (*, '(a,2es14.6)') 'MAW_ASUM=', maw, asum

    end program drv_vturb
