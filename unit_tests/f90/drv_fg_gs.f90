    program drv_fg_gs
    use set_kind_module
    use fg_j_cmf_nonmono_mod_v1
    use mod_cmfgen, only: r, v, sigma
    use control_variable_mod, only: conv_diag_lvl, gs_foreign_seed
    implicit none

    integer, parameter :: nd_in = 4, ncomp_in = 3, np_in = 2, &
        ncf_in = 5, lu = 30
    real(kind=ldp) :: r_l(nd_in), v_l(nd_in), sg_l(nd_in)
    real(kind=ldp) :: r_alt(nd_in)
    real(kind=ldp), allocatable :: ip_orig(:, :, :), im_orig(:, :, :)
    character(len=32) :: case_name
    logical :: ok_flag

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))

    call setup_state()

    select case (case_name)
    case ('dump_only')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        write (*, '(a)') 'DUMP_ONLY_DONE'

    case ('restore_only')
        fg_gs_tag_exp = 7
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok

    case ('rt')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        inquire (file='FGINEXT_GS', exist=ok_flag)
        write (*, '(a,l1)') 'FILE_OK=', ok_flag
        inquire (file='FGINEXT_GS.tmp', exist=ok_flag)
        write (*, '(a,l1)') 'TMP_GONE=', .not. ok_flag
        i_p_depth_next = -77.0_ldp
        i_m_depth_next = -55.0_ldp
        i_next_gs_ok = .false.
        fg_gs_chg_file = -1.0_ldp
        fg_gs_cap_file = -1
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok
        write (*, '(a,l1)') 'IPBIT=', all(i_p_depth_next == ip_orig)
        write (*, '(a,l1)') 'IMBIT=', all(i_m_depth_next == im_orig)
        write (*, '(a,l1)') 'SEEDED=', all(fg_ray_seeded == 1)
        write (*, '(a,es12.5)') 'CHG_FILE=', fg_gs_chg_file
        write (*, '(a,i4)') 'CAP_FILE=', fg_gs_cap_file

    case ('no_file')
        fg_gs_tag_exp = 7
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok

    case ('no_tag')
        fg_gs_tag_exp = -1
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok

    case ('stale')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        fg_gs_tag_exp = 99
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok

    case ('foreign')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        fg_gs_tag_exp = 99
        gs_foreign_seed = .true.
        i_p_depth_next = -9.0_ldp
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok
        write (*, '(a,l1)') 'IPBIT=', all(i_p_depth_next == ip_orig)

    case ('bad_rv')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        r_alt = r_l
        r_alt(2) = r_l(2)*(1.0_ldp + 1.0e-12_ldp)
        call restore_fg_inext_gs(r_alt, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok

    case ('bad_cmp')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        comp_depth_idx(1) = comp_depth_idx(1) + 1
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok

    case ('neg_i')
        fg_gs_tag_exp = 7
        i_p_depth_next(1, 1, 1) = -5.0_ldp
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        i_p_depth_next = 3.0_ldp
        call restore_fg_inext_gs(r_l, v_l, sg_l, nd_in, np_in, ncf_in)
        write (*, '(a,l1)') 'REST_OK=', i_next_gs_ok
        write (*, '(a,l1)') 'IP_ZEROED=', all(i_p_depth_next == 0.0_ldp)

    case ('latest')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        i_p_depth_next = -66.0_ldp
        i_m_depth_next = -44.0_ldp
        fg_gs_tag_exp = 99
        call restore_fg_inext_latest(ok_flag)
        write (*, '(a,l1)') 'LATEST_OK=', ok_flag
        write (*, '(a,l1)') 'IPBIT=', all(i_p_depth_next == ip_orig)
        write (*, '(a,l1)') 'IMBIT=', all(i_m_depth_next == im_orig)

    case ('latest_mm')
        fg_gs_tag_exp = 7
        call dump_fg_inext_gs(7, r_l, v_l, sg_l, nd_in, lu)
        i_p_depth_next = -33.0_ldp
        sigma(1) = sigma(1)*(1.0_ldp + 1.0e-10_ldp)
        call restore_fg_inext_latest(ok_flag)
        write (*, '(a,l1)') 'LATEST_OK=', ok_flag
        write (*, '(a,l1)') 'IP_UNTOUCHED=', &
            all(i_p_depth_next == -33.0_ldp)

    case default
        write (*, '(a)') 'UNKNOWN_CASE'
    end select

contains

    subroutine setup_state()
        integer :: ii, kk, ll
        allocate (r(nd_in), v(nd_in), sigma(nd_in))
        r = [12.0_ldp, 9.0_ldp, 6.0_ldp, 3.0_ldp]
        v = [40.0_ldp, 30.0_ldp, 20.0_ldp, 10.0_ldp]
        sigma = [0.5_ldp, -1.5_ldp, 0.2_ldp, 0.8_ldp]
        r_l = r
        v_l = v
        sg_l = sigma
        allocate (i_p_depth_next(ncomp_in, np_in, ncf_in))
        allocate (i_m_depth_next(ncomp_in, np_in, ncf_in))
        allocate (ip_orig(ncomp_in, np_in, ncf_in))
        allocate (im_orig(ncomp_in, np_in, ncf_in))
        allocate (kn_freq(ncf_in), comp_depth_idx(nd_in))
        allocate (fg_ray_seeded(np_in))
        kn_freq = [50.0_ldp, 40.0_ldp, 30.0_ldp, 20.0_ldp, 10.0_ldp]
        comp_depth_idx = [1, 0, 2, 0]
        fg_ray_seeded = 0
        do ii = 1, ncomp_in
            do kk = 1, np_in
                do ll = 1, ncf_in
                    ip_orig(ii, kk, ll) = (ii - 1)*100.0_ldp &
                        + (kk - 1)*10.0_ldp + ll
                    im_orig(ii, kk, ll) = (ii - 1)*100.0_ldp &
                        + (kk - 1)*10.0_ldp + ll*0.25_ldp
                end do
            end do
        end do
        i_p_depth_next = ip_orig
        i_m_depth_next = im_orig
        n_comp = ncomp_in
        kn_nd = nd_in
        fg_gs_np = np_in
        fg_gs_ncf = ncf_in
        fg_gs_maxchg = 2.0e-3_ldp
        fg_inext_capped = 3_8
        fg_gs_ndiv = 1
        fg_commit_rej = 0
        fg_div_heal = 1
        conv_diag_lvl = 0
        gs_foreign_seed = .false.
    end subroutine setup_state

    end program drv_fg_gs
