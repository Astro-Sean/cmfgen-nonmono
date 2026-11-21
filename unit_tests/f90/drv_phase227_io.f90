program drv_phase227_io
    use fg_j_cmf_nonmono_mod_v1
    implicit none
    character(len=4096) :: mode, path
    integer :: unit, ray
    integer :: classes(4) = [0, 1, 3, 5]
    real(kind=ldp) :: residuals(2), omegas(2)

    call get_command_argument(1, mode)
    call get_command_argument(2, path)
    residuals = [0.001_ldp, 0.02_ldp]
    omegas = [1.0_ldp, 0.51412345678901234_ldp]
    select case (trim(mode))
    case ('omega')
        call kn_replay_load(trim(path), 81, 0)
        if (any(kn_replay_count /= [2, 1])) error stop 'unexpected counts'
    case ('census')
        open(newunit=unit, file=trim(path), status='new')
        do ray = 1, 4
            call kn_census_row(unit, 0_8, ray, classes(ray), 99, -1, &
                huge(0_8), huge(1.0_ldp), 1.0e-100_ldp, 4, &
                residuals, [1, 3], [2, 100000], omegas, .true.)
        end do
        close(unit)
    case ('roundtrip')
        kn_stage_enabled = .true.
        kn_stage_ls = 81
        call kn_omega_applied(81, 1, 1, omegas(2), omegas(2))
    case default
        error stop 'unknown test mode'
    end select
end program drv_phase227_io
