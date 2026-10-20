module knop_status_mod
!
! Shared status sentences for the non-monotonic lagged-field
! machinery (KNOP and the SOB moment path -- both drive
! MOM_J_CMF_NONMONO_V1).  Each pure function returns one complete,
! human-readable sentence of at most 200 characters; the caller
! writes it with a plain (A) format.  Keeping the wording in one
! module means the sweep report, the refresh narration, the
! checkpoint-restore messages and the population-gate sentences all
! say the same thing everywhere they print.
!
! No machine records are produced here: GS_STATE, MOMJNXT_GS and
! the ITER_SUM/ITER_HISTORY rows keep their own formats.
!
    use set_kind_module, only: ldp
    implicit none
    private

    public :: mom_state_sentence
    public :: mom_sweep_sentence
    public :: mom_refresh_sentence
    public :: mom_pass_sentence
    public :: mom_restart_sentence
    public :: mom_direction_sentence

    character(len=*), parameter :: DIR_CMP = 'compressing'
    character(len=*), parameter :: DIR_EXP = 'expanding'

contains

    ! Left-adjusted ES-format for embedding a real in a sentence.
    ! len 13 / es13.4: a negative 3-digit exponent needs 12 chars
    ! (e.g. -1.2345E-105), which overran len 10 and printed stars.
    pure function es(x) result(str)
        real(kind=ldp), intent(in) :: x
        character(len=13) :: str
        write (str, '(es13.4)') x
        str = adjustl(str)
    end function es

    ! len 24: kind-8 cell counts (GS_NTOT ~ 6.5E8 now, larger on
    ! expanded grids) plus sign fit without overflow.
    pure function i8(x) result(str)
        integer(kind=8), intent(in) :: x
        character(len=24) :: str
        write (str, '(i0)') x
    end function i8

    ! Population-gate narration for the moment field.  Selects the
    ! single clause describing why the call held, released or
    ! postponed the population update; 'it' is the great-iteration
    ! counter used in the held sentence.
    pure function mom_state_sentence(freeze, held, unsafe, broad, &
            trial, precond, rf_conv, it) result(str)
        logical, intent(in) :: freeze, held, unsafe, broad
        logical, intent(in) :: trial, precond, rf_conv
        integer, intent(in) :: it
        character(len=200) :: str
        character(len=12) :: num
        write (num, '(i0)') it
        if (unsafe) then
            ! 'holding' is only true when GS_HOLD_POPS is on; with the
            ! hold off the capped LAMBDA update still commits.
            if (held) then
                str = 'Moment field unsafe - catastrophic drift; ' // &
                      'holding populations.'
            else
                str = 'Moment field unsafe - catastrophic drift; ' // &
                      'hold disabled.'
            end if
        else if (precond) then
            str = 'Populations far from equilibrium - postponing the ' // &
                  'lagged-field hold and taking a live LAMBDA step.'
        else if (trial) then
            str = 'Moment field marginal - trying a full step with ' // &
                  'the rollback net armed.'
        else if (broad) then
            str = 'Moment field drifting over most cells - holding ' // &
                  'populations.'
        else if (freeze .and. held) then
            str = 'Moment field not stable - holding populations ' // &
                  'while the lagged J_NEXT field converges ' // &
                  '(iteration ' // trim(num) // ').'
        else if (freeze) then
            str = 'Moment field not stable - populations may move ' // &
                  'this step (hold disabled).'
        else if (rf_conv) then
            str = 'Moment refresh converged.'
        else
            str = 'Moment field stable.'
        end if
    end function mom_state_sentence

    ! Per-sweep lagged-field summary.  view='pops' describes the
    ! field the population solve saw (no drift number: the census is
    ! against the stored field); view='sweep' describes the fresh
    ! solve versus the stored J_NEXT, led by the largest drift.
    pure function mom_sweep_sentence(view, maxchg, nbad5, nbad50, &
            negj, nquar, ntot, npquar) result(str)
        character(len=*), intent(in) :: view
        real(kind=ldp), intent(in) :: maxchg
        integer(kind=8), intent(in) :: nbad5, nbad50, negj
        integer(kind=8), intent(in) :: nquar, ntot, npquar
        character(len=200) :: str
        character(len=200) :: tail
        tail = trim(i8(nbad5)) // ' cells moved more than 5%, ' // &
            trim(i8(nbad50)) // ' more than 50%; ' // &
            trim(i8(negj)) // ' cells with negative J_NEXT ' // &
            '(repaired); ' // trim(i8(nquar)) // ' cells quarantined'
        if (view .eq. 'pops') then
            str = 'Moment field, as seen by the populations: ' // &
                trim(tail) // '; ' // trim(i8(ntot)) // &
                ' significant cells; ' // trim(i8(npquar)) // &
                ' physically quarantined.'
        else
            str = 'Moment field, fresh solve versus stored ' // &
                'J_NEXT: largest drift ' // trim(es(maxchg)) // &
                '; ' // trim(tail) // '.'
        end if
    end function mom_sweep_sentence

    ! Refresh-outcome narration.  exit_kind selects the sentence:
    ! 'converged'  converged at a pass
    ! 'bulk'       bulk-converged exit (n50 carries the >50% count)
    ! 'stall'      stopped on a stall
    ! 'diverge'    divergence abort (best carries the best residual)
    ! 'keep_best'  unconverged, best-pass field retained
    ! 'keep_prior' unconverged, pre-refresh field retained
    ! 'delta'      retained-field delta written back to the store
    ! 'done'       pass count and wall time (best carries seconds)
    pure function mom_refresh_sentence(conv, dchg, n5, n50, &
            exit_kind, ipass, best) result(str)
        logical, intent(in) :: conv
        real(kind=ldp), intent(in) :: dchg, best
        integer(kind=8), intent(in) :: n5, n50
        character(len=*), intent(in) :: exit_kind
        integer, intent(in) :: ipass
        character(len=200) :: str
        character(len=12) :: num
        write (num, '(i0)') ipass
        select case (exit_kind)
        case ('converged')
            str = 'Moment refresh converged at pass ' // &
                trim(num) // ' (residual ' // trim(es(dchg)) // ').'
        case ('bulk')
            str = 'Moment refresh bulk-converged at pass ' // &
                trim(num) // ' (' // trim(i8(n50)) // &
                ' cells above 50%).'
        case ('stall')
            str = 'Moment refresh stalled at pass ' // trim(num) // &
                ' (residual ' // trim(es(dchg)) // ').'
        case ('diverge')
            str = 'Moment refresh diverged at pass ' // trim(num) // &
                ' (residual ' // trim(es(dchg)) // &
                ' versus the best ' // trim(es(best)) // '); aborted.'
        case ('keep_best')
            str = 'Moment refresh did not converge - retained the ' // &
                'best pass (residual ' // trim(es(dchg)) // '; ' // &
                trim(i8(n5)) // ' cells above 5%, ' // &
                trim(i8(n50)) // ' above 50%).'
        case ('keep_prior')
            str = 'Moment refresh did not converge - retained the ' // &
                'pre-refresh field (best residual ' // &
                trim(es(dchg)) // ').'
        case ('delta')
            if (conv) then
                str = 'Moment refresh wrote the stored field: ' // &
                    'converged; largest change ' // &
                    trim(es(dchg)) // '; ' // trim(i8(n5)) // &
                    ' cells moved more than 5%, ' // &
                    trim(i8(n50)) // ' more than 50%.'
            else
                str = 'Moment refresh wrote the stored field: ' // &
                    'not converged; largest change ' // &
                    trim(es(dchg)) // '; ' // trim(i8(n5)) // &
                    ' cells moved more than 5%, ' // &
                    trim(i8(n50)) // ' more than 50%.'
            end if
        case ('done')
            write (str, '(a,i0,a,f0.2,a)') &
                'Moment refresh finished: ', ipass, ' passes in ', &
                best, ' s.'
        case default
            str = 'Moment refresh: residual ' // trim(es(dchg)) // '.'
        end select
    end function mom_refresh_sentence

    ! One refresh pass.  idir picks the direction name: 2 is the
    ! compressing direction (descending frequency order, upwind
    ! reads from J_NEXT), anything else is expanding.
    pure function mom_pass_sentence(ipass, idir, chg, ratio, n5, &
            n50) result(str)
        integer, intent(in) :: ipass, idir
        real(kind=ldp), intent(in) :: chg, ratio
        integer(kind=8), intent(in) :: n5, n50
        character(len=200) :: str
        character(len=12) :: dname, num
        if (idir .eq. 2) then
            dname = DIR_CMP
        else
            dname = DIR_EXP
        end if
        write (num, '(i0)') ipass
        str = 'Moment pass ' // trim(num) // ' (' // trim(dname) // &
            '): largest change ' // trim(es(chg)) // &
            '; same-direction contraction ratio ' // &
            trim(es(ratio)) // '; ' // trim(i8(n5)) // &
            ' cells above 5%, ' // trim(i8(n50)) // ' above 50%.'
    end function mom_pass_sentence

    ! Lagged-field checkpoint restore narration.  'reason' selects
    ! the rejection clause when valid is false: 'stale' (tag or
    ! foreign checkpoint), 'dims' (dimension or expanded-grid
    ! mismatch), 'data' (non-finite or negative values), 'read'
    ! (record/format failure).
    pure function mom_restart_sentence(present, valid, reason, &
            tag, expected) result(str)
        logical, intent(in) :: present, valid
        character(len=*), intent(in) :: reason
        integer, intent(in) :: tag, expected
        character(len=200) :: str
        character(len=12) :: tnum, enum
        write (tnum, '(i0)') tag
        write (enum, '(i0)') expected
        if (.not. present) then
            str = 'J_NEXT checkpoint: absent - cold start, ' // &
                  'rebuilding from zero.'
        else if (valid) then
            str = 'J_NEXT checkpoint: present, accepted ' // &
                '(iteration tag ' // trim(tnum) // ').'
        else
            select case (reason)
            case ('stale')
                str = 'J_NEXT checkpoint: rejected because it is ' // &
                    'stale or from another model (tag ' // &
                    trim(tnum) // ' versus expected ' // &
                    trim(enum) // '); rebuilding from zero.'
            case ('data')
                str = 'J_NEXT checkpoint: rejected because the ' // &
                      'data are non-finite or negative; ' // &
                      'rebuilding from zero.'
            case ('read')
                str = 'J_NEXT checkpoint: rejected because the ' // &
                      'file could not be read or parsed; ' // &
                      'rebuilding from zero.'
            case default
                str = 'J_NEXT checkpoint: rejected because the ' // &
                      'dimensions or the expanded grid do not ' // &
                      'match; rebuilding from zero.'
            end select
        end if
    end function mom_restart_sentence

    ! Resolved refresh pass ordering.
    pure function mom_direction_sentence(dom_cmp, dom_exp) result(str)
        logical, intent(in) :: dom_cmp, dom_exp
        character(len=200) :: str
        if (dom_cmp) then
            str = 'Pass order: compressing direction only ' // &
                  '(dominant-direction schedule).'
        else if (dom_exp) then
            str = 'Pass order: expanding direction only ' // &
                  '(dominant-direction schedule).'
        else
            str = 'Pass order: alternating compressing and ' // &
                  'expanding directions.'
        end if
    end function mom_direction_sentence

end module knop_status_mod
