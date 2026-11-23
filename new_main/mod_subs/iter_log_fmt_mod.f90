module iter_log_fmt_mod
!
! Shared layout for the per-iteration OUTGEN report. Section
! banners, key/value rows and the end-of-iteration status line are
! emitted through this module so the widths stay identical no
! matter which routine does the printing. Formatting only: no
! numerical state is stored or altered here.
!
! Verbosity rides on the existing CONV_DIAG_LVL VADAT key
! (default 2, clamped 0..2 by RD_CONTROL_VARIABLES):
!   0 - iteration framing, the ITER_SUM block, Model Convergence,
!       timing, and the IT status banner; the summary survives at
!       level 0 by design (see SOLVE_FOR_POPS).
!   1 - the above plus section banners and the per-call status
!       lines (solve type, max corrections, luminosity, RUNAWAY
!       table, ES_ALI/JCOMP census, JFNK verdicts).
!   2 - everything, including the trace noise (setup chatter,
!       non-thermal species dump, DEC/INC vectors, the
!       BA-CTX/MODCTL/MODE-OUT/BA13-TRACE dev blocks).
! Lines that parsers consume keep their exact wording at every
! level that prints them; the level only decides whether they
! print, never how they are spelled.
!
    use set_kind_module, only: ldp
    implicit none
    private

    public :: ilf_rule, ilf_section
    public :: ilf_kv, ilf_kvd, ilf_kvi
    public :: ilf_set_banner, ilf_banner

    integer, parameter :: ilf_lw = 80          ! section/rule width

    ! Status-banner state exported by SOLVE_FOR_POPS once per great
    ! iteration and read by CMFGEN_SUB at the timing line. ILF_IT
    ! guards staleness: a call that never ran the mode controller
    ! (stock path) prints no banner at all.
    character(len=16) :: ilf_mode_a = ' '
    character(len=16) :: ilf_mode_b = ' '
    character(len=8)  :: ilf_ba     = ' '
    integer           :: ilf_npct   = -1
    integer           :: ilf_it     = 0

contains

    ! One 80-column rule of the fill character.
    subroutine ilf_rule(lu, fill)
        integer, intent(in) :: lu
        character(len=1), intent(in) :: fill
        write (lu, '(a)') repeat(fill, ilf_lw)
    end subroutine ilf_rule

    ! '--- TITLE ---' padded with dashes to 80 columns, with a
    ! blank line before and after so sections stay visually
    ! separated regardless of what precedes them.
    subroutine ilf_section(lu, title)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: title
        character(len=ilf_lw) :: line
        line = ' '
        if (len_trim(title) + 5 <= ilf_lw) then
            line = '--- '//trim(title)//' '//repeat('-', ilf_lw)
        else
            line = '--- '//trim(title)
        end if
        write (lu, '(a)') ' '
        write (lu, '(a)') line
        write (lu, '(a)') ' '
    end subroutine ilf_section

    ! '  key                          : value' -- label field 30
    ! wide so the colon column is fixed for every key/value row.
    ! The key is buffered to 30 chars first: an A30 edit on a short
    ! item right-justifies, which would pull the label away from
    ! the left margin.
    subroutine ilf_kv(lu, key, val)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: key, val
        character(len=30) :: kb
        kb = key
        write (lu, '(2x,a30,a,a)') kb, ': ', trim(val)
    end subroutine ilf_kv

    subroutine ilf_kvd(lu, key, x)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: key
        real(kind=ldp), intent(in) :: x
        character(len=30) :: kb
        kb = key
        write (lu, '(2x,a30,a,es12.4)') kb, ': ', x
    end subroutine ilf_kvd

    subroutine ilf_kvi(lu, key, n)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: key
        integer, intent(in) :: n
        character(len=30) :: kb
        kb = key
        write (lu, '(2x,a30,a,i0)') kb, ': ', n
    end subroutine ilf_kvi

    ! Record this call's status words for the banner. Called from
    ! SOLVE_FOR_POPS once CONV_NPCT exists, inside the same
    ! dev-path gate as ITER_SUM.
    subroutine ilf_set_banner(mode_a, mode_b, ba, npct, it)
        character(len=*), intent(in) :: mode_a, mode_b, ba
        integer, intent(in) :: npct, it
        ilf_mode_a = mode_a
        ilf_mode_b = mode_b
        ilf_ba     = ba
        ilf_npct   = npct
        ilf_it     = it
    end subroutine ilf_set_banner

    ! End-of-iteration status line, printed next to 'Time taken
    ! for iteration'. Skipped when no banner state was exported
    ! for this counter (stock path or a call that stopped early).
    subroutine ilf_banner(lu, it, tsec)
        integer, intent(in) :: lu, it
        real(kind=ldp), intent(in) :: tsec
        character(len=12) :: tbuf
        if (ilf_it /= it) return
        write (tbuf, '(f8.1)') tsec
        write (lu, '(a)') ' '
        call ilf_rule(lu, '-')
        if (ilf_npct >= 0) then
            write (lu, '(a,i0,7a,i0,3a)') ' IT ', it, ' | ', &
                trim(ilf_mode_a), ' -> ', trim(ilf_mode_b), &
                ' | BA ', trim(ilf_ba), ' | Conv ', ilf_npct, &
                '% | ', trim(adjustl(tbuf)), ' s'
        else
            write (lu, '(a,i0,9a)') ' IT ', it, ' | ', &
                trim(ilf_mode_a), ' -> ', trim(ilf_mode_b), &
                ' | BA ', trim(ilf_ba), ' | ', &
                trim(adjustl(tbuf)), ' s'
        end if
        call ilf_rule(lu, '-')
        write (lu, '(a)') ' '
    end subroutine ilf_banner

end module iter_log_fmt_mod
