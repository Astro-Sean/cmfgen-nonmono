    program drv_add_dust
    use set_kind_module
    use mod_space_grid_v2
    implicit none

    integer, parameter :: nd = 4, nray = 5, np = 8, ip = 2
    real(kind=ldp) :: dust_scat_opac(nd), r(nd), eta(nray)
    real(kind=ldp) :: eta_p_ray(nray), eta_m_ray(nray)
    real(kind=ldp) :: eta_p(nd), eta_m(nd)
    character(len=32) :: case_name
    integer :: i, k, nmu

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    r = [2.0_ldp, 4.0_ldp, 6.0_ldp, 8.0_ldp]
    dust_scat_opac = [1.0e-2_ldp, 2.0e-2_ldp, 3.0e-2_ldp, 4.0e-2_ldp]
    eta = [1.0e-4_ldp, 2.0e-4_ldp, 3.0e-4_ldp, 4.0e-4_ldp, 5.0e-4_ldp]

    do k = 1, np
        allocate (ray(k)%lnk(nd), ray(k)%i_p_save(nd), &
            ray(k)%i_m_save(nd))
        ray(k)%lnk = [1, 2, 3, 4]
        ray(k)%i_p_save = [(10.0_ldp*k + i, i = 1, nd)]
        ray(k)%i_m_save = [(20.0_ldp*k + i, i = 1, nd)]
    end do

    allocate (dw(ip)%wgt_p_to_p(np, nd), dw(ip)%wgt_p_to_m(np, nd))
    allocate (dw(ip)%wgt_m_to_p(np, nd), dw(ip)%wgt_m_to_m(np, nd))
    do k = 1, np
        do i = 1, nd
            dw(ip)%wgt_p_to_p(k, i) = 1.0e-5_ldp*(k + i)
            dw(ip)%wgt_m_to_p(k, i) = 2.0e-5_ldp*(k - i + 0.5_ldp)
            dw(ip)%wgt_p_to_m(k, i) = 3.0e-5_ldp*(k + 2*i)
            dw(ip)%wgt_m_to_m(k, i) = 4.0e-5_ldp*(k - 0.25_ldp*i)
        end do
    end do

    ray(ip)%nz = nray
    allocate (ray(ip)%r_ray(nray))
    select case (case_name)
    case ('outside')
        ray(ip)%r_ray = [2.5_ldp, 3.0_ldp, 4.0_ldp, 6.0_ldp, 8.0_ldp]
    case ('crossing')
        ray(ip)%r_ray = [2.0_ldp, 3.0_ldp, 4.0_ldp, 6.0_ldp, 8.0_ldp]
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    eta_p_ray = -9.0_ldp
    eta_m_ray = -9.0_ldp

    call add_eta_dust(dust_scat_opac, eta, eta_m_ray, eta_p_ray, &
        ip, nray, r, nd, np)

    write (*, '(a)') 'DONE_' // case_name
    write (*, '(a,5es16.8)') 'EPRAY=', eta_p_ray
    write (*, '(a,5es16.8)') 'EMRAY=', eta_m_ray

    end program drv_add_dust
