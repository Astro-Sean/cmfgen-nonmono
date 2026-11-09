    program drv_rd_cont_j
    use set_kind_module
    use mod_cmfgen, only: r, v
    use opac_mod
    use radiation_mod
    use ang_qw_mod
    implicit none

    integer, parameter :: nd_in = 4, np_in = 3, luer = 6, luedd = 31
    real(kind=ldp) :: arr(nd_in)
    character(len=32) :: case_name
    integer :: access_f, i, ios

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    allocate (r(nd_in), v(nd_in))
    r = [12.0_ldp, 9.0_ldp, 6.0_ldp, 3.0_ldp]
    v = [40.0_ldp, 30.0_ldp, 20.0_ldp, 10.0_ldp]
    allocate (p(np_in), p_obs(np_in), iplus(np_in))
    p = [1.0_ldp, 5.0_ldp, 12.0_ldp]
    allocate (rj(nd_in), k_mom(nd_in), rsqhnu(nd_in))
    allocate (source(nd_in), zeta(nd_in), theta(nd_in))
    rj = -1.0_ldp
    source = -1.0_ldp
    zeta = 0.5_ldp
    theta = 0.1_ldp

    open (unit=luedd, file='EDDF_TST', form='unformatted', &
        access='direct', recl=(nd_in + 1)*2, status='replace', &
        iostat=ios)
    if (ios /= 0) then
        write (*, '(a,i5)') 'OPEN_FAIL=', ios
        stop
    end if

    select case (case_name)
    case ('torn_hi')
        access_f = 1
        call rd_cont_j(5.0_ldp, 1, .true., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call dump_state('TORN_HI')

    case ('torn_lo')
        arr = 8.0_ldp
        write (luedd, rec=1) arr, 10.0_ldp
        access_f = 1
        call rd_cont_j(5.0_ldp, 1, .true., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call dump_state('TORN_LO')

    case ('first_ok')
        arr = 8.0_ldp
        write (luedd, rec=1) arr, 10.0_ldp
        arr = 4.0_ldp
        write (luedd, rec=2) arr, 5.0_ldp
        access_f = 1
        call rd_cont_j(7.0_ldp, 1, .true., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call dump_state('FIRST_OK')
        write (*, '(a,i3)') 'ACC=', access_f
        write (*, '(a,2es14.6)') 'P_OBS=', p_obs(1), p_obs(np_in)

    case ('hi_clamp')
        arr = 8.0_ldp
        write (luedd, rec=1) arr, 10.0_ldp
        arr = 4.0_ldp
        write (luedd, rec=2) arr, 5.0_ldp
        access_f = 1
        call rd_cont_j(12.0_ldp, 1, .true., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call dump_state('HI_CLAMP')

    case ('walk')
        arr = 8.0_ldp
        write (luedd, rec=1) arr, 10.0_ldp
        arr = 4.0_ldp
        write (luedd, rec=2) arr, 5.0_ldp
        arr = 2.0_ldp
        write (luedd, rec=3) arr, 2.0_ldp
        access_f = 1
        call rd_cont_j(7.0_ldp, 1, .true., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call rd_cont_j(4.0_ldp, 2, .false., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call dump_state('WALK')
        write (*, '(a,i3)') 'ACC=', access_f

    case ('walk_torn')
        arr = 8.0_ldp
        write (luedd, rec=1) arr, 10.0_ldp
        arr = 4.0_ldp
        write (luedd, rec=2) arr, 5.0_ldp
        access_f = 1
        call rd_cont_j(7.0_ldp, 1, .true., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call rd_cont_j(4.0_ldp, 2, .false., .false., .false., &
            luer, luedd, access_f, nd_in, np_in)
        call dump_state('WALK_TORN')

    case ('accurate')
        access_f = 1
        call rd_cont_j(5.0_ldp, 1, .true., .false., .true., &
            luer, luedd, access_f, nd_in, np_in)
        write (*, '(a)') 'ACC_NO_STOP'

    case default
        write (*, '(a)') 'UNKNOWN_CASE'
    end select

contains

    subroutine dump_state(tag)
        character(len=*), intent(in) :: tag
        write (*, '(a)') 'DONE_' // tag
        write (*, '(a,4es14.6)') 'RJ=', rj
        write (*, '(a,l1)') 'KMOM_OK=', all(k_mom == rj/3.0_ldp)
        write (*, '(a,l1)') 'SRC_OK=', &
            all(source == zeta + theta*rj)
        write (*, '(a,l1)') 'RSQ_OK=', all(rsqhnu == 1.0e-20_ldp)
        write (*, '(a,l1)') 'HBC_OK=', all(hbc_cmf == 1.0_ldp)
        write (*, '(a,l1)') 'IPLUS_OK=', all(iplus == 1.0e-10_ldp)
    end subroutine dump_state

    end program drv_rd_cont_j
