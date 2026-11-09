    program drv_get_pops_v5
    use set_kind_module
    use mod_cmfgen
    implicit none

    integer, parameter :: nd_in = 5, nd_old = 6, nt_in = 4, lu = 30
    real(kind=ldp) :: pops(nt_in, nd_in), old_pops(nt_in, nd_in)
    real(kind=ldp) :: old_r(nd_in)
    real(kind=ldp) :: fr(nd_old), fv(nd_old), fsg(nd_old)
    real(kind=ldp) :: fpa(nd_old), fdens(nd_old), ft(nd_old), fed(nd_old)
    real(kind=ldp) :: sn_age
    character(len=32) :: case_name
    integer :: i

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    num_ions = 1
    species_beg_id = 0
    species_end_id = 0
    species = 'ZZZ'
    allocate (v(nd_in), t(nd_in))
    allocate (old_lev_pop_avail(nt_in), vol_exp_fac(nd_in))
    v = [40.0_ldp, 32.5_ldp, 25.0_ldp, 17.5_ldp, 10.0_ldp]
    t = 3.0_ldp
    vol_exp_fac = 1.0_ldp

    fv = [40.0_ldp, 35.0_ldp, 30.0_ldp, 25.0_ldp, 20.0_ldp, 10.0_ldp]
    fr = [12.0_ldp, 11.0_ldp, 10.0_ldp, 9.0_ldp, 8.0_ldp, 3.0_ldp]
    fsg = 0.5_ldp
    ft = [5.0_ldp, 4.6_ldp, 4.2_ldp, 3.8_ldp, 3.4_ldp, 3.0_ldp]
    fed = [2.0_ldp, 1.8_ldp, 1.6_ldp, 1.4_ldp, 1.2_ldp, 1.0_ldp]
    fpa = 1.0_ldp
    fdens = 0.1_ldp
    sn_age = 4.0_ldp

    select case (case_name)
    case ('good')
    case ('v_neg')
        v(2) = -1.0_ldp
    case ('v_flat')
        v(3) = v(2)
    case ('oldv_neg')
        fv(2) = -1.0_ldp
    case ('oldv_flat')
        fv(3) = fv(2)
    case ('v_ib_mismatch')
        fv(nd_old) = 9.9_ldp
    case ('outer_small')
        fv(1) = v(1)*0.9_ldp
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    call write_model_file()

    pops = 1.0_ldp
    old_pops = -1.0_ldp
    old_r = -1.0_ldp

    call get_pops_at_prev_time_step_v5(pops, old_pops, old_r, &
        .false., .false., .false., 2.0_ldp, nd_in, nt_in, lu)

    write (*, '(a)') 'COMPLETED'
    write (*, '(a,2es16.8)') 'OLD_R_END=', old_r(1), old_r(nd_in)
    write (*, '(a,l1)') 'OLD_R_MONO=', &
        all(old_r(2:nd_in) < old_r(1:nd_in - 1))
    write (*, '(a,l1)') 'T_FIN=', &
        all(old_pops(nt_in, :) > 0.0_ldp .and. &
            old_pops(nt_in, :) < huge(1.0_ldp))
    write (*, '(a,2es16.8)') 'T_END=', &
        old_pops(nt_in, 1), old_pops(nt_in, nd_in)
    write (*, '(a,l1)') 'ED_ZERO=', all(old_pops(nt_in - 1, :) == 0.0_ldp)

contains

    subroutine write_model_file()
        open (unit=lu, file='OLD_MODEL_DATA', form='unformatted', &
            status='replace')
        write (lu) '21-Jul-2008'
        write (lu) nd_old, 0
        write (lu) sn_age
        write (lu) fr
        write (lu) fv
        write (lu) fsg
        write (lu) ft
        write (lu) fed
        write (lu) fpa
        write (lu) fdens
        close (lu)
    end subroutine write_model_file

    end program drv_get_pops_v5
