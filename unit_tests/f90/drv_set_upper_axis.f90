    program drv_set_upper_axis
    use set_kind_module
    implicit none

    integer, parameter :: nd = 5, nsc = 31
    real(kind=ldp) :: ed(nd), xv(nd)
    character(len=30) :: top_lab
    character(len=32) :: case_name
    common /topbord/ sced(nsc), xed(nsc), nxed, toplabel
    real(kind=ldp) :: sced, xed
    integer :: nxed
    character(len=30) :: toplabel

    call get_command_argument(1, case_name)
    case_name = trim(adjustl(case_name))
    xv = [1.0_ldp, 2.0_ldp, 3.0_ldp, 4.0_ldp, 5.0_ldp]
    top_lab = 'TEST'

    select case (case_name)
    case ('good')
        ed = [1.0e4_ldp, 1.0e5_ldp, 1.0e6_ldp, 1.0e7_ldp, 1.0e8_ldp]
    case ('flat')
        ed = [1.0e4_ldp, 1.0e5_ldp, 1.0e5_ldp, 1.0e7_ldp, 1.0e8_ldp]
    case ('drop')
        ed = [1.0e4_ldp, 1.0e6_ldp, 1.0e5_ldp, 1.0e7_ldp, 1.0e8_ldp]
    case ('neg')
        ed = [-1.0e4_ldp, 1.0e5_ldp, 1.0e6_ldp, 1.0e7_ldp, 1.0e8_ldp]
    case ('zero')
        ed = [0.0_ldp, 1.0e5_ldp, 1.0e6_ldp, 1.0e7_ldp, 1.0e8_ldp]
    case default
        write (*, '(a)') 'UNKNOWN_CASE'
        stop
    end select

    call set_upper_axis(ed, xv, nd, top_lab)
    write (*, '(a)') 'DONE_' // case_name
    write (*, '(a,i3)') 'NXED=', nxed
    write (*, '(a,a)') 'LAB=', toplabel

    end program drv_set_upper_axis
