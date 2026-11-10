    program drv_adiab
    use set_kind_module
    use mod_cmfgen
    use steq_data_mod
    use control_variable_mod, only: use_elec_heat_bal, &
        comp_steq_t_ehb, do_nonmono_sob, do_nonmono_knop
    implicit none

    integer, parameter :: nd = 5, nt = 6, nbnds = 3, didx = 2
    real(kind=ldp) :: ad_cr_v(nd), ad_cr_dt(nd)
    real(kind=ldp) :: pops(nt, nd), ave_energy(nt)
    character(len=32) :: case_name
    integer :: i

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    allocate (r(nd), v(nd), t(nd), ed(nd), pop_atom(nd))
    allocate (sigma(nd))
    allocate (steq_t(nd), ba_t(nt, nbnds, nd))
    allocate (steq_t_ehb(nd), ba_t_ehb(nt, nbnds, nd))

    r = [20.0_ldp, 15.0_ldp, 10.0_ldp, 6.0_ldp, 3.0_ldp]
    v = [900.0_ldp, 700.0_ldp, 500.0_ldp, 300.0_ldp, 150.0_ldp]
    t = [3.0_ldp, 4.0_ldp, 5.0_ldp, 6.5_ldp, 8.0_ldp]
    ed = [1.0e5_ldp, 3.0e5_ldp, 8.0e5_ldp, 2.0e6_ldp, &
        5.0e6_ldp]
    pop_atom = 1.0e10_ldp
    sigma = [-0.5_ldp, -0.5_ldp, -0.5_ldp, -0.5_ldp, -0.5_ldp]
    species_beg_id = 0
    species_end_id = 0
    steq_t = 0.0_ldp
    ba_t = 0.0_ldp
    steq_t_ehb = 0.0_ldp
    ba_t_ehb = 0.0_ldp
    use_elec_heat_bal = .false.
    comp_steq_t_ehb = .false.
    pops = 0.0_ldp
    ave_energy = [(1.0_ldp*i, i = 1, nt)]

    do_nonmono_sob = .false.
    do_nonmono_knop = .false.
    select case (case_name)
    case ('sob_true')
        do_nonmono_sob = .true.
        allocate (sigma_true(nd))
        sigma_true = [0.3_ldp, -0.8_ldp, 0.5_ldp, -0.2_ldp, &
            0.7_ldp]
    case ('sob_unalloc')
        do_nonmono_sob = .true.
    case ('knop_true')
        do_nonmono_knop = .true.
        allocate (sigma_true(nd))
        sigma_true = [0.3_ldp, -0.8_ldp, 0.5_ldp, -0.2_ldp, &
            0.7_ldp]
    case ('off')
        allocate (sigma_true(nd))
        sigma_true = [0.3_ldp, -0.8_ldp, 0.5_ldp, -0.2_ldp, &
            0.7_ldp]
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    call eval_adiabatic_v3(ad_cr_v, ad_cr_dt, pops, ave_energy, &
        1.0e-2_ldp, .true., .true., didx, nbnds, nt, nd)

    write (*, '(a)') 'DONE_' // case_name
    write (*, '(a,5es14.6)') 'ADCRV=', ad_cr_v
    write (*, '(a,5es14.6)') 'ADCRDT=', ad_cr_dt

    end program drv_adiab
