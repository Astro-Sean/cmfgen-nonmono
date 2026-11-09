    program drv_mom_gs
    use set_kind_module
    use mod_mom_j_nonmono_v1
    use mod_cmfgen, only: r, v
    use control_variable_mod, only: conv_diag_lvl
    implicit none

    integer, parameter :: nd_in = 4, ndnm_in = 6, ncf_in = 4, lu = 30
    real(kind=ldp) :: sigma(nd_in)
    real(kind=ldp) :: r_alt(nd_in), v_alt(nd_in)
    real(kind=ldp), allocatable :: j_orig(:, :), h_orig(:, :)
    integer(kind=2), allocatable :: q_orig(:, :)
    character(len=32) :: case_name
    logical :: ok_flag

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    call setup_state()

    select case (case_name)
    case ('dump_only')
        mom_gs_tag_exp = 7
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        write (*, '(a)') 'DUMP_ONLY_DONE'

    case ('restore_only')
        mom_gs_tag_exp = 7
        j_next_frq = -99.0_ldp
        call restore_mom_jnext_gs(r, v, nd_in)
        write (*, '(a,l1,a,l1)') 'REST_OK=', mom_j_gs_ok, &
            ' JVALID=', j_next_valid
        write (*, '(a,es12.5)') 'J11=', j_next_frq(1, 1)
        write (*, '(a,es12.5)') 'CHG_FILE=', mom_gs_chg_file

    case ('rt')
        mom_gs_tag_exp = 7
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        inquire (file='MOMJNXT_GS', exist=ok_flag)
        write (*, '(a,l1)') 'FILE_OK=', ok_flag
        inquire (file='MOMJNXT_GS.tmp', exist=ok_flag)
        write (*, '(a,l1)') 'TMP_GONE=', .not. ok_flag
        j_next_frq = -77.0_ldp
        h_next_frq = -55.0_ldp
        mom_quar_age = 0_2
        mom_j_gs_ok = .false.
        j_next_valid = .false.
        mom_gs_chg_file = -1.0_ldp
        mom_gs_cap_file = -1
        call restore_mom_jnext_gs(r, v, nd_in)
        write (*, '(a,l1,a,l1)') 'REST_OK=', mom_j_gs_ok, &
            ' JVALID=', j_next_valid
        write (*, '(a,l1)') 'JBIT=', all(j_next_frq == j_orig)
        write (*, '(a,l1)') 'HBIT=', all(h_next_frq == h_orig)
        write (*, '(a,l1)') 'QBIT=', all(mom_quar_age == q_orig)
        write (*, '(a,es12.5)') 'CHG_FILE=', mom_gs_chg_file
        write (*, '(a,i4)') 'CAP_FILE=', mom_gs_cap_file

    case ('no_file')
        mom_gs_tag_exp = 7
        call restore_mom_jnext_gs(r, v, nd_in)
        write (*, '(a,l1)') 'REST_OK=', mom_j_gs_ok

    case ('no_tag')
        mom_gs_tag_exp = -1
        call restore_mom_jnext_gs(r, v, nd_in)
        write (*, '(a,l1)') 'REST_OK=', mom_j_gs_ok

    case ('stale')
        mom_gs_tag_exp = 7
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        mom_gs_tag_exp = 99
        call restore_mom_jnext_gs(r, v, nd_in)
        write (*, '(a,l1)') 'REST_OK=', mom_j_gs_ok

    case ('bad_rv')
        mom_gs_tag_exp = 7
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        r_alt = r
        v_alt = v
        r_alt(1) = r(1)*1.01_ldp
        call restore_mom_jnext_gs(r_alt, v_alt, nd_in)
        write (*, '(a,l1)') 'REST_OK=', mom_j_gs_ok

    case ('grid_mm')
        mom_gs_tag_exp = 7
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        r_nm(1) = r_nm(1)*(1.0_ldp + 1.0e-10_ldp)
        call restore_mom_jnext_gs(r, v, nd_in)
        write (*, '(a,l1)') 'REST_OK=', mom_j_gs_ok

    case ('neg_j')
        mom_gs_tag_exp = 7
        j_next_frq(2, 1) = -5.0_ldp
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        j_next_frq = 3.0_ldp
        call restore_mom_jnext_gs(r, v, nd_in)
        write (*, '(a,l1)') 'REST_OK=', mom_j_gs_ok
        write (*, '(a,l1)') 'J_ZEROED=', all(j_next_frq == 0.0_ldp)

    case ('latest')
        mom_gs_tag_exp = 7
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        j_next_frq = -66.0_ldp
        h_next_frq = -44.0_ldp
        j_next_valid = .false.
        mom_gs_tag_exp = 99
        call restore_mom_jnext_latest(ok_flag)
        write (*, '(a,l1)') 'LATEST_OK=', ok_flag
        write (*, '(a,l1)') 'JBIT=', all(j_next_frq == j_orig)
        write (*, '(a,l1)') 'HBIT=', all(h_next_frq == h_orig)

    case ('latest_mm')
        mom_gs_tag_exp = 7
        call dump_mom_jnext_gs(7, r, v, sigma, nd_in, lu)
        j_next_frq = -33.0_ldp
        r_nm(1) = r_nm(1)*(1.0_ldp + 1.0e-10_ldp)
        call restore_mom_jnext_latest(ok_flag)
        write (*, '(a,l1)') 'LATEST_OK=', ok_flag
        write (*, '(a,l1)') 'J_UNTOUCHED=', all(j_next_frq == -33.0_ldp)

    case default
        write (*, '(a)') 'UNKNOWN_CASE'
    end select

contains

    subroutine setup_state()
        integer :: ii, kk
        allocate (r(nd_in), v(nd_in))
        r = [12.0_ldp, 9.0_ldp, 6.0_ldp, 3.0_ldp]
        v = [40.0_ldp, 30.0_ldp, 20.0_ldp, 10.0_ldp]
        sigma = 0.5_ldp
        allocate (j_next_frq(ndnm_in, ncf_in), h_next_frq(ndnm_in, ncf_in))
        allocate (r_nm(ndnm_in), r_pnt_nm(ndnm_in))
        allocate (mom_quar_age(ndnm_in, ncf_in))
        allocate (j_orig(ndnm_in, ncf_in), h_orig(ndnm_in, ncf_in))
        allocate (q_orig(ndnm_in, ncf_in))
        r_nm = [12.0_ldp, 10.5_ldp, 9.0_ldp, 7.5_ldp, 6.0_ldp, 3.0_ldp]
        r_pnt_nm = [1, 1, 2, 3, 3, 4]
        do ii = 1, ndnm_in
            do kk = 1, ncf_in
                j_orig(ii, kk) = (ii - 1)*10.0_ldp + kk
                h_orig(ii, kk) = (ii - 1)*10.0_ldp - kk*0.5_ldp
                q_orig(ii, kk) = int(ii + kk, kind=2)
            end do
        end do
        j_next_frq = j_orig
        h_next_frq = h_orig
        mom_quar_age = q_orig
        mom_gs_ncf = ncf_in
        mom_gs_maxchg = 1.0e-3_ldp
        mom_inext_capped = 2_8
        mom_rf_ran = .false.
        conv_diag_lvl = 0
    end subroutine setup_state

    end program drv_mom_gs
