    program drv_tau_regrid
    use set_kind_module
    use mod_cmfgen
    implicit none

    integer, parameter :: nd = 8, nt = 4
    real(kind=ldp) :: pops(nt, nd), esec(nd)
    logical :: done_r_rev
    character(len=32) :: case_name
    integer :: i, lu_opts

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    allocate (r(nd), v(nd), sigma(nd), flux_mean(nd), &
        clump_fac(nd))
    r = [20.0_ldp, 15.0_ldp, 11.0_ldp, 8.0_ldp, 5.8_ldp, &
        4.2_ldp, 3.0_ldp, 2.0_ldp]
    v = [900.0_ldp, 700.0_ldp, 540.0_ldp, 400.0_ldp, &
        290.0_ldp, 210.0_ldp, 150.0_ldp, 100.0_ldp]
    sigma = -0.5_ldp
    esec = 1.0e-3_ldp
    clump_fac = 1.0_ldp
    pops = 0.0_ldp

    select case (case_name)
    case ('zero_head')
        flux_mean = [0.0_ldp, 1.0e-2_ldp, 3.0e-2_ldp, &
            6.0e-2_ldp, 1.2e-1_ldp, 3.0e-1_ldp, 8.0e-1_ldp, &
            2.0e0_ldp]
    case ('normal')
        flux_mean = [1.0e-2_ldp, 1.4e-2_ldp, 3.0e-2_ldp, &
            6.0e-2_ldp, 1.2e-1_ldp, 3.0e-1_ldp, 8.0e-1_ldp, &
            2.0e0_ldp]
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    open (newunit=lu_opts, file='tau_opts.txt', status='old', &
        action='read')
    call rd_options_into_store(lu_opts, 6)
    close (lu_opts)

    call do_tau_regrid_v2(pops, esec, done_r_rev, nd, nt, 6)

    write (*, '(a)') 'DONE_' // case_name
    write (*, '(a,l1)') 'REV=', done_r_rev
    write (*, '(a,8es14.6)') 'RNEW=', r

    end program drv_tau_regrid
