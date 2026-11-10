! Stubs for dead-tool dependencies of par_opac_sub.o.
! These routines sit in the standalone par_opac.exe link graph
! but are not archived anywhere in the current build; they are
! only reachable past the V monotonicity guard, which the test
! never crosses.

    subroutine set_prof_limits_v2()
    stop 'stub set_prof_limits_v2 reached'
    end

    subroutine set_prof_v3()
    stop 'stub set_prof_v3 reached'
    end

    subroutine set_two_phot()
    stop 'stub set_two_phot reached'
    end

    subroutine write_direct_info()
    stop 'stub write_direct_info reached'
    end
