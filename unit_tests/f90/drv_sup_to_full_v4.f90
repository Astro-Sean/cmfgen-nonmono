    program drv_sup_to_full_v4
    use set_kind_module
    use mod_cmfgen
    use mod_lev_dis_blk, only: mod_hold_lev_dis, b_lev_dis, &
        a_lev_dis, x_lev_dis
    implicit none

    integer, parameter :: nd = 2, nf = 5, ns = 2, nt_in = 4
    real(kind=ldp) :: pops(nt_in, nd), z_pop(nt_in)
    real(kind=ldp) :: save_s(ns, nd), sum_l, err, errmax
    character(len=32) :: case_name
    integer :: i, k, l

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    num_ions = 2
    allocate (popion(nd), t(nd), ed(nd))
    t = [3.0_ldp, 6.0_ldp]
    ed = [1.0e6_ldp, 1.0e8_ldp]

    call fill_atom()

    pops(1, :) = [1.0e-3_ldp, 2.0e-4_ldp]
    pops(2, :) = [4.0e-5_ldp, 8.0e-3_ldp]
    pops(3, :) = [9.9e9_ldp, 9.9e9_ldp]
    pops(4, :) = [7.0e-2_ldp, 3.0e-2_ldp]
    z_pop = [1.0_ldp, 1.0_ldp, 0.005_ldp, 1.0_ldp]
    save_s(1, :) = atm(1)%xzv(1, :)
    save_s(2, :) = atm(1)%xzv(2, :)

    select case (case_name)
    case ('hold_f')
        mod_hold_lev_dis = .false.
    case ('hold_t')
        mod_hold_lev_dis = .true.
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    popion = 777.0_ldp
    allocate (b_lev_dis(nd), a_lev_dis(nd), x_lev_dis(nd))
    b_lev_dis = 888.0_ldp
    a_lev_dis = 888.0_ldp
    x_lev_dis = 888.0_ldp

    call sup_to_full_v4(pops, z_pop, .false., nd, nt_in)

    write (*, '(a)') 'DONE_' // case_name
    write (*, '(a,2es14.6)') 'POPION=', popion
    write (*, '(a,2es14.6)') 'BLD=', b_lev_dis

    errmax = 0.0_ldp
    do k = 1, nd
        do l = 1, ns
            sum_l = 0.0_ldp
            do i = 1, nf
                if (atm(1)%f_to_s_xzv(i) == l) &
                    sum_l = sum_l + atm(1)%xzv_f(i, k)
            end do
            err = abs(sum_l - save_s(l, k))/save_s(l, k)
            if (err > errmax) errmax = err
        end do
    end do
    write (*, '(a,es12.4)') 'CONS_ERR=', errmax

contains

    subroutine fill_atom()
        integer :: i
        allocate (atm(1)%xzv_f(nf, nd), atm(1)%xzvlte_f(nf, nd))
        allocate (atm(1)%log_xzvlte_f(nf, nd))
        allocate (atm(1)%xzvlte_f_on_s(nf, nd))
        allocate (atm(1)%w_xzv_f(nf, nd), atm(1)%edgexzv_f(nf))
        allocate (atm(1)%gxzv_f(nf), atm(1)%f_to_s_xzv(nf))
        allocate (atm(1)%int_seq_xzv(nf), atm(1)%dxzv_f(nd))
        allocate (atm(1)%xzv(ns, nd), atm(1)%xzvlte(ns, nd))
        allocate (atm(1)%log_xzvlte(ns, nd), atm(1)%dxzv(nd))
        allocate (atm(2)%edgexzv_f(1), atm(2)%gxzv_f(1))
        allocate (atm(2)%f_to_s_xzv(1))
        atm(1)%edgexzv_f = [3.0_ldp, 2.6_ldp, 1.9_ldp, 1.1_ldp, &
            0.4_ldp]
        atm(1)%gxzv_f = [2.0_ldp, 6.0_ldp, 4.0_ldp, 2.0_ldp, &
            10.0_ldp]
        atm(1)%f_to_s_xzv = [1, 1, 2, 2, 2]
        atm(1)%int_seq_xzv = [0, 0, 0, 0, 0]
        atm(1)%nxzv_f = nf
        atm(1)%nxzv = ns
        atm(1)%gionxzv_f = 1.0_ldp
        atm(1)%zxzv = 1.0_ldp
        atm(1)%xzv_pres = .true.
        atm(1)%xzv(1, :) = [1.0e-3_ldp, 2.0e-4_ldp]
        atm(1)%xzv(2, :) = [4.0e-5_ldp, 8.0e-3_ldp]
        atm(1)%xzv_f = 0.0_ldp
        atm(1)%w_xzv_f = 1.0_ldp
        atm(1)%dxzv_f = 1.0e-4_ldp
        atm(1)%dxzv = 1.0e-4_ldp
        do i = 1, nf
            atm(1)%xzvlte_f(i, :) = atm(1)%gxzv_f(i) &
                *exp(-atm(1)%edgexzv_f(i)*4.7994145_ldp/t) &
                *1.0e-3_ldp
        end do
        atm(1)%log_xzvlte_f = log(atm(1)%xzvlte_f)
        atm(1)%xzvlte_f_on_s = atm(1)%xzvlte_f
        atm(1)%xzvlte = 1.0e-3_ldp
        atm(1)%log_xzvlte = log(atm(1)%xzvlte)
        atm(2)%edgexzv_f = [0.0_ldp]
        atm(2)%gxzv_f = [1.0_ldp]
        atm(2)%f_to_s_xzv = [1]
        atm(2)%nxzv_f = 1
        atm(2)%xzv_pres = .false.
    end subroutine fill_atom

    end program drv_sup_to_full_v4
