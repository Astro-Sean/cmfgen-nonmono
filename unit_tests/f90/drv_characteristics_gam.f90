    program drv_characteristics_gam
    use set_kind_module
    use mod_space_grid_v2_gam
    implicit none

    integer, parameter :: nd = 20, nc = 6, np = nd + nc
    real(kind=ldp) :: r_grid(nd), p(np), v_grid(nd), vdop_vec(nd)
    real(kind=ldp) :: vdop_frac, t1
    character(len=32) :: case_name
    integer :: i, k, bad, nz_tot
    real(kind=ldp) :: mumin, mumax, rmin

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

!   Same grid convention as misc/tst_char_v2.f: R(1) is the outer
!   boundary, R(ND) the inner; core rays fill P < R(ND) and each
!   tangent ray lands exactly on a grid shell.
    r_grid(1) = 20.0_ldp
    r_grid(nd) = 2.0_ldp
    t1 = exp(log(r_grid(1)/r_grid(nd))/(nd - 1))
    do i = 2, nd
        r_grid(i) = r_grid(nd)*(t1**(nd - i))
    end do
    vdop_vec = 20.0_ldp
    vdop_frac = 0.5_ldp

    t1 = r_grid(nd)/nc
    do i = 1, nc
        p(i) = t1*(i - 1)
    end do
    do i = nc + 1, np
        p(i) = r_grid(nd + nc + 1 - i)
    end do

    select case (case_name)
    case ('mono')
        v_grid = 50.0_ldp*(r_grid/r_grid(nd))
    case ('nudge')
!       A tangent ray whose impact radius sits a hair above its
!       shell: upstream sqrt goes NaN, the clamp keeps it real.
        v_grid = 50.0_ldp*(r_grid/r_grid(nd))
        p(nc + 5) = r_grid(nd - 4)*(1.0_ldp + 2.0e-4_ldp)
        p(np - 1) = r_grid(2)*(1.0_ldp + 2.0e-4_ldp)
    case ('dip')
        v_grid = 50.0_ldp*(r_grid/r_grid(nd))
        v_grid(6) = 0.25_ldp*v_grid(6)
        v_grid(7) = 0.5_ldp*v_grid(7)
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    call characteristics_v2_gam(r_grid, p, v_grid, vdop_vec, vdop_frac, &
        nd, nc, np)

    nz_tot = 0
    write (*, '(a,l1)') 'INS=', ray_points_inserted
    do k = 1, np
        if (allocated(ray(k)%s_p)) then
            bad = 0
            mumin = 2.0_ldp
            mumax = -2.0_ldp
            rmin = 1.0e30_ldp
            do i = 1, ray(k)%nz
                if (ray(k)%s_p(i) /= ray(k)%s_p(i)) bad = bad + 1
                if (ray(k)%mu_p(i) /= ray(k)%mu_p(i)) bad = bad + 1
                if (ray(k)%b_p(i) /= ray(k)%b_p(i)) bad = bad + 1
                if (ray(k)%r_ray(i) /= ray(k)%r_ray(i)) &
                    bad = bad + 1
                mumin = min(mumin, ray(k)%mu_p(i))
                mumax = max(mumax, ray(k)%mu_p(i))
                rmin = min(rmin, ray(k)%r_ray(i))
            end do
            nz_tot = nz_tot + ray(k)%nz
            write (*, '(a,i2,a,i4,a,i2,a,2es12.4,a,es12.4)') &
                'RAY', k, ' NZ=', ray(k)%nz, ' BAD=', bad, &
                ' MU=', mumin, mumax, ' RMIN=', rmin
        else
            write (*, '(a,i2,a)') 'RAY', k, ' NOT_ALLOC'
        end if
    end do
    write (*, '(a,i6)') 'NZTOT=', nz_tot

    write (*, '(a)') 'DONE_' // case_name

    end program drv_characteristics_gam
