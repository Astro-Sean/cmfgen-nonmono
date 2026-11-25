module iter_log_fmt_mod
!
! Shared layout for the per-iteration OUTGEN report. Section
! banners, key/value rows, tables, tagged diagnostics and the
! end-of-iteration status line are emitted through this module so
! the widths stay identical no matter which routine does the
! printing. Formatting only: no numerical state is stored or
! altered here.
!
! Two key/value styles exist. ilf_kv* emits the fixed 'key : value'
! rows used by the iteration framing. ilf_kvb emits the vertical
! 'label = value' block that converts multi-field one-liners: the
! header text is kept, labels pad to the longest in the block, and
! the ilf_v* builders right-align each value field (I10 integers,
! ES reals at the caller's decimals, T/F logicals). LOG_COMPACT_KV
! collapses a block back to the legacy one-line form for diffs.
!
! Verbosity rides on OUT_LEVEL (preferred VADAT spelling) or the
! legacy CONV_DIAG_LVL alias, clamped 0..3 by RD_CONTROL_VARIABLES:
!   0 - iteration framing, ITER_SUM, Model Convergence, timing,
!       the IT status banner, and the ISSUES block.
!   1 - the above plus section banners, the moment-refresh table
!       and summary, the RUNAWAY table, and per-call status lines.
!   2 - the former full output volume in the new layout, including
!       one compact detail block for the last moment-refresh pass.
!   3 - all per-pass moment-refresh detail blocks.
! Lines that parsers consume keep their exact wording at every
! level that prints them; the level only decides whether they
! print, never how they are spelled.
!
    use set_kind_module, only: ldp
    implicit none
    private

    public :: ilf_rule, ilf_section, ilf_subsection
    public :: ilf_kv, ilf_kvd, ilf_kvi
    public :: ilf_kvb, ilf_vi, ilf_vr, ilf_vl, ilf_vs
    public :: ilf_val_len
    public :: log_compact_kv
    public :: ilf_set_banner, ilf_banner
    public :: ilf_es, ilf_es_trim, ilf_tag
    public :: ilf_table_header, ilf_table_row, ilf_metric_row
    public :: ilf_metric_header
    public :: ilf_issue_reset, ilf_issue_report
    public :: ilf_repeat_reset, ilf_repeat_write, ilf_repeat_flush

    integer, parameter :: ilf_lw = 80          ! section/rule width
    integer, parameter :: ilf_maxline = 100    ! output line budget
    integer, parameter :: ilf_val_len = 40     ! ilf_v* value buffer
    integer, parameter :: ilf_maxissue = 128
    integer, parameter :: ilf_msglen = 180

    ! Layout switch for ilf_kvb blocks.  .false. prints the vertical
    ! 'label = value' block; .true. prints the legacy one-line
    ! 'label=val, label=val' style so a converted log can be diffed
    ! against a pre-conversion run.  Compact lines keep the legacy
    ! length and may exceed ilf_maxline; the 100-column budget only
    ! binds the block layout.  Flipped by the driver/test code, not
    ! by VADAT: no key exists for it yet.
    logical :: log_compact_kv = .false.

    ! Status-banner state exported by SOLVE_FOR_POPS once per great
    ! iteration and read by CMFGEN_SUB at the timing line. ILF_IT
    ! guards staleness: a call that never ran the mode controller
    ! (stock path) prints no banner at all.
    character(len=16) :: ilf_mode_a = ' '
    character(len=16) :: ilf_mode_b = ' '
    character(len=8)  :: ilf_ba     = ' '
    integer           :: ilf_npct   = -1
    integer           :: ilf_it     = 0

    ! Per-iteration issue census. Exact level/source/message triples
    ! are counted; values embedded in the message therefore make
    ! distinct entries, while unchanged repeated diagnostics fold.
    integer :: ilf_nissue = 0
    integer :: ilf_issue_overflow = 0
    integer :: ilf_issue_count(ilf_maxissue) = 0
    character(len=8) :: ilf_issue_level(ilf_maxissue) = ' '
    character(len=16) :: ilf_issue_source(ilf_maxissue) = ' '
    character(len=ilf_msglen) :: ilf_issue_msg(ilf_maxissue) = ' '

    ! One-slot run-length filter used by hot-loop diagnostics.
    logical :: ilf_rep_active = .false.
    integer :: ilf_rep_count = 0
    integer :: ilf_rep_lu = -1
    character(len=8) :: ilf_rep_level = ' '
    character(len=16) :: ilf_rep_source = ' '
    character(len=ilf_msglen) :: ilf_rep_msg = ' '

contains

    ! Fixed-width scientific notation. Two exponent digits are used
    ! whenever the exponent fits; three digits are used beyond that,
    ! still in a fixed 11-character field. This avoids both width
    ! drift and '*****' overflow at the values CMFGEN prints.
    pure function ilf_es(x) result(str)
        real(kind=ldp), intent(in) :: x
        character(len=11) :: str
        if (x /= x) then
            str = '        NaN'
        else if (abs(x) >= 1.0e100_ldp .or. &
                 (x /= 0.0_ldp .and. abs(x) < 1.0e-99_ldp)) then
            write (str, '(es11.3e3)') x
        else
            write (str, '(es11.4e2)') x
        end if
    end function ilf_es

    pure function ilf_es_trim(x) result(str)
        real(kind=ldp), intent(in) :: x
        character(len=11) :: str
        str = adjustl(ilf_es(x))
    end function ilf_es_trim

    ! One 80-column rule of the fill character.
    subroutine ilf_rule(lu, fill)
        integer, intent(in) :: lu
        character(len=1), intent(in) :: fill
        write (lu, '(a)') repeat(fill, ilf_lw)
    end subroutine ilf_rule

    ! '--- TITLE ---' padded with dashes to 80 columns, with a
    ! blank line before and after. Blank lines are zero length, not
    ! whitespace-only records.
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
        write (lu, '(a)') ''
        write (lu, '(a)') trim(line)
        write (lu, '(a)') ''
    end subroutine ilf_section

    ! Single subsection style: '-- TITLE ' followed by dashes.
    subroutine ilf_subsection(lu, title)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: title
        character(len=ilf_lw) :: line
        line = ' '
        if (len_trim(title) + 4 <= ilf_lw) then
            line = '-- '//trim(title)//' '//repeat('-', ilf_lw)
        else
            line = '-- '//trim(title)
        end if
        write (lu, '(a)') trim(line)
    end subroutine ilf_subsection

    ! '  key                          : value' -- label field 30
    ! wide so the colon column is fixed for every key/value row.
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
        write (lu, '(2x,a30,a,a)') kb, ': ', trim(ilf_es_trim(x))
    end subroutine ilf_kvd

    subroutine ilf_kvi(lu, key, n)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: key
        integer, intent(in) :: n
        character(len=30) :: kb
        kb = key
        write (lu, '(2x,a30,a,i0)') kb, ': ', n
    end subroutine ilf_kvi

    ! ------------------------------------------------------------------
    ! Vertical key/value blocks: a header line followed by one
    ! 'label = value' row per field.  The ilf_v* builders return the
    ! value field already right-aligned inside an ilf_val_len buffer;
    ! leading blanks are part of the field so the '=' column only
    ! depends on the label padding.
    ! ------------------------------------------------------------------

    ! Integer field, right-aligned in I10.  Values that do not fit
    ! fall back to the minimal I0 form rather than '*****'.
    pure function ilf_vi(n) result(str)
        integer, intent(in) :: n
        character(len=ilf_val_len) :: str
        character(len=16) :: tmp
        str = ' '
        write (tmp, '(i10)') n
        if (index(tmp, '*') > 0) then
            write (tmp, '(i0)') n
            str = ' '//adjustl(tmp)
        else
            str(1:10) = tmp
        end if
    end function ilf_vi

    ! Real field, right-aligned in an ES field of width
    ! max(10, ndec+8) with a two-digit exponent.  NDEC is the
    ! decimals the old format carried so the printed precision is
    ! unchanged.  Three exponent digits are used only when two do
    ! not fit (|x|>=1E100 or underflow below 1E-99), so '*****'
    ! never appears.  NaN and +/-Inf print as 'NaN'/'Inf'/'-Inf',
    ! right-aligned in the same field.
    pure function ilf_vr(x, ndec) result(str)
        real(kind=ldp), intent(in) :: x
        integer, intent(in) :: ndec
        character(len=ilf_val_len) :: str
        character(len=16) :: fmt
        integer :: w, d
        str = ' '
        d = max(0, ndec)
        w = max(10, d + 8)
        if (x /= x) then
            str(w - 2:w) = 'NaN'
            return
        end if
        if (abs(x) > huge(x)) then
            if (x > 0.0_ldp) then
                str(w - 2:w) = 'Inf'
            else
                str(w - 3:w) = '-Inf'
            end if
            return
        end if
        ! Three-digit exponents need the E3 form: inside an E2
        ! format the processor drops the E ('1.23+100'), which is
        ! less readable than the widened field.
        if (abs(x) >= 1.0e100_ldp .or. &
            (x /= 0.0_ldp .and. abs(x) < 1.0e-99_ldp)) then
            write (fmt, '("(es",i0,".",i0,"e3)")') &
                min(w + 1, ilf_val_len), d
        else
            write (fmt, '("(es",i0,".",i0,"e2)")') w, d
        end if
        write (str, fmt) x
        if (index(str, '*') > 0) then
            str = ' '
            str(w - 2:w) = 'ovf'
        end if
    end function ilf_vr

    ! Logical field, 'T'/'F' right-aligned in 10 columns like I10 so
    ! mixed blocks keep a tidy value edge.
    pure function ilf_vl(l) result(str)
        logical, intent(in) :: l
        character(len=ilf_val_len) :: str
        str = ' '
        if (l) then
            str(10:10) = 'T'
        else
            str(10:10) = 'F'
        end if
    end function ilf_vl

    ! Short string field: leading blanks are stripped and one blank
    ! separates '=' from the text.  Content past the buffer is
    ! truncated; string fields longer than ilf_val_len-1 need a
    ! hand-written row, not this helper.
    pure function ilf_vs(s) result(str)
        character(len=*), intent(in) :: s
        character(len=ilf_val_len) :: str
        str = ' '//adjustl(s)
    end function ilf_vs

    ! 'header' then N rows '    label  = value'.  Labels pad to the
    ! longest label in the block so every '=' sits in one column.
    ! Rows are indented 4 columns from an unindented header (a
    ! header's own leading indent carries over).  If the natural
    ! label width would push a row past ilf_maxline the padding is
    ! shrunk; a single label too long to fit even at zero padding
    ! puts its value on an extra-indented continuation line instead
    ! of overflowing the budget.  With LOG_COMPACT_KV=.true. the
    ! whole block collapses to the legacy 'header label=val, ...'
    ! one-liner.
    subroutine ilf_kvb(lu, header, n, labels, vals)
        integer, intent(in) :: lu
        integer, intent(in) :: n
        character(len=*), intent(in) :: header
        character(len=*), intent(in) :: labels(:), vals(:)
        integer :: i, m, ind0, ind, pad, padcap, vmax, lt, lv, pos
        character(len=320) :: line
        character(len=ilf_val_len + 64) :: frag
        m = min(n, size(labels), size(vals))
        ind0 = 0
        do while (ind0 < len(header))
            if (header(ind0 + 1:ind0 + 1) /= ' ') exit
            ind0 = ind0 + 1
        end do
        ind = max(4, ind0 + 3)
        if (log_compact_kv) then
            line = ' '
            if (len_trim(header) > 0) then
                line = header(1:len_trim(header))
            end if
            pos = len_trim(line) + 1
            do i = 1, m
                frag = ' '//trim(adjustl(labels(i)))//'='// &
                    trim(adjustl(vals(i)))
                if (i > 1) frag = ','//frag
                if (pos + len_trim(frag) > len(line)) exit
                line(pos:) = trim(frag)
                pos = len_trim(line) + 1
            end do
            write (lu, '(a)') trim(line)
            return
        end if
        if (len_trim(header) > 0) then
            write (lu, '(a)') header(1:len_trim(header))
        end if
        if (m <= 0) return
        pad = 0
        vmax = 0
        do i = 1, m
            pad = max(pad, len_trim(labels(i)))
            vmax = max(vmax, len_trim(vals(i)))
        end do
        padcap = ilf_maxline - ind - 3 - vmax
        pad = min(pad, max(0, padcap))
        do i = 1, m
            lt = len_trim(labels(i))
            lv = len_trim(vals(i))
            if (lt <= pad) then
                line = repeat(' ', ind)//labels(i)(1:pad)//'  ='// &
                    vals(i)(1:lv)
            else if (ind + lt + 3 + lv <= ilf_maxline) then
                line = repeat(' ', ind)//labels(i)(1:lt)//'  ='// &
                    vals(i)(1:lv)
            else
                ! Label too long to share the row with its value:
                ! split rather than truncate the field text.
                call ilf_issue_add('WARN', 'ILF', &
                    'kv row over line budget; value wrapped')
                write (lu, '(a)') repeat(' ', ind)// &
                    labels(i)(1:lt)//'  ='
                line = repeat(' ', ind + 4)//vals(i)(1:lv)
            end if
            write (lu, '(a)') trim(line)
        end do
    end subroutine ilf_kvb

    ! Tagged diagnostic line. WARN/BAD lines also enter the
    ! per-iteration issue census; printing is never conditional here,
    ! so the caller retains ownership of when the diagnostic exists.
    subroutine ilf_tag(lu, level, source, msg)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: level, source, msg
        character(len=8) :: lev
        character(len=16) :: src
        lev = adjustl(level)
        src = adjustl(source)
        write (lu, '(2x,a)') '['//trim(lev)//']['//trim(src)//'] '//trim(msg)
        if (trim(lev) == 'WARN' .or. trim(lev) == 'BAD') then
            call ilf_issue_add(trim(lev), trim(src), trim(msg))
        end if
    end subroutine ilf_tag

    subroutine ilf_issue_add(level, source, msg)
        character(len=*), intent(in) :: level, source, msg
        integer :: j
        do j = 1, ilf_nissue
            if (ilf_issue_level(j) == level .and. &
                ilf_issue_source(j) == source) then
                ilf_issue_count(j) = ilf_issue_count(j) + 1
                return
            end if
        end do
        if (ilf_nissue < ilf_maxissue) then
            ilf_nissue = ilf_nissue + 1
            j = ilf_nissue
            ilf_issue_level(j) = level
            ilf_issue_source(j) = source
            ilf_issue_msg(j) = msg
            ilf_issue_count(j) = 1
        else
            ilf_issue_overflow = ilf_issue_overflow + 1
        end if
    end subroutine ilf_issue_add

    subroutine ilf_issue_reset()
        ilf_nissue = 0
        ilf_issue_overflow = 0
        ilf_issue_count = 0
        ilf_issue_level = ' '
        ilf_issue_source = ' '
        ilf_issue_msg = ' '
    end subroutine ilf_issue_reset

    ! End-of-iteration issue census. N counts distinct level/source
    ! tags; the first message is retained as an example. The block
    ! intentionally exists even when empty so a
    ! reader never has to infer whether warnings were lost.
    subroutine ilf_issue_report(lu)
        integer, intent(in) :: lu
        integer :: j
        write (lu, '(a,i0,a)') 'ISSUES (', ilf_nissue, ')'
        do j = 1, ilf_nissue
            if (ilf_issue_count(j) > 1) then
                write (lu, '(2x,a,i0,a)') '['//trim(ilf_issue_level(j))// &
                    ']['//trim(ilf_issue_source(j))//'] '// &
                    trim(ilf_issue_msg(j))//' (x', ilf_issue_count(j), ')'
            else
                write (lu, '(2x,a)') '['//trim(ilf_issue_level(j))// &
                    ']['//trim(ilf_issue_source(j))//'] '// &
                    trim(ilf_issue_msg(j))
            end if
        end do
        if (ilf_issue_overflow > 0) then
            write (lu, '(2x,a,i0,a)') '[WARN][LOG] issue census full; ', &
                ilf_issue_overflow, ' further occurrences omitted'
        end if
        call ilf_issue_reset()
    end subroutine ilf_issue_report

    ! Generic fixed-width table row emitter. WIDTHS includes the
    ! visible cell width; one space separates cells. RIGHT selects
    ! right alignment (useful for numbers). Long cells are truncated
    ! rather than allowed to break the column boundary.
    subroutine ilf_table_row(lu, cols, widths, right, indent)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: cols(:)
        integer, intent(in) :: widths(:)
        logical, intent(in), optional :: right(:)
        integer, intent(in), optional :: indent
        integer :: j, n, pos, w, lenc, ind
        character(len=256) :: line
        character(len=64) :: cell
        logical :: ralign
        n = min(size(cols), size(widths))
        ind = 0
        if (present(indent)) ind = max(0, indent)
        line = ''
        if (ind > 0) line(1:ind) = repeat(' ', ind)
        pos = ind + 1
        do j = 1, n
            w = widths(j)
            ralign = .false.
            if (present(right)) then
                if (j <= size(right)) ralign = right(j)
            end if
            cell = adjustl(cols(j))
            lenc = len_trim(cell)
            if (pos + w - 1 > len(line)) exit
            if (lenc == 0) then
                line(pos:pos+w-1) = repeat(' ', w)
            else if (ralign) then
                if (lenc < w) then
                    line(pos:pos+w-lenc-1) = repeat(' ', w-lenc)
                    line(pos+w-lenc:pos+w-1) = cell(1:lenc)
                else
                    line(pos:pos+w-1) = cell(1:w)
                end if
            else
                line(pos:pos+w-1) = cell(1:min(w, lenc))
            end if
            pos = pos + w + 1
        end do
        write (lu, '(a)') trim(line)
    end subroutine ilf_table_row

    subroutine ilf_table_header(lu, cols, widths, indent)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: cols(:)
        integer, intent(in) :: widths(:)
        integer, intent(in), optional :: indent
        integer :: j, n, w, ind, pos
        character(len=256) :: sep
        call ilf_table_row(lu, cols, widths, indent=indent)
        ind = 0
        if (present(indent)) ind = max(0, indent)
        n = min(size(cols), size(widths))
        sep = ''
        if (ind > 0) sep(1:ind) = repeat(' ', ind)
        pos = ind + 1
        do j = 1, n
            w = widths(j)
            if (pos + w - 1 > len(sep)) exit
            sep(pos:pos+w-1) = repeat('-', w)
            pos = pos + w + 1
        end do
        write (lu, '(a)') trim(sep)
    end subroutine ilf_table_header

    ! Standard metric/A/B/status/note row used by ITER_SUM and the
    ! solver censuses. Widths keep the whole row below 100 columns.
    subroutine ilf_metric_header(lu)
        integer, intent(in) :: lu
        character(len=12) :: cols(5)
        integer :: widths(5)
        cols = [character(len=12) :: 'metric', 'A', 'B', 'status', 'note']
        widths = [14, 11, 11, 7, 38]
        call ilf_table_header(lu, cols, widths, indent=2)
    end subroutine ilf_metric_header

    subroutine ilf_metric_row(lu, metric, a, b, status, note)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: metric, a, b, status, note
        character(len=40) :: cols(5)
        integer :: widths(5)
        logical :: right(5)
        cols(1) = metric
        cols(2) = a
        cols(3) = b
        cols(4) = status
        cols(5) = note
        widths = [14, 11, 11, 7, 38]
        right = [.false., .true., .true., .false., .false.]
        call ilf_table_row(lu, cols, widths, right=right, indent=2)
    end subroutine ilf_metric_row

    ! Run-length compression for consecutive identical tagged
    ! messages. The first instance is buffered until the next unlike
    ! message or an explicit flush, so ordering with unrelated writes
    ! is the caller's responsibility: flush before leaving the hot
    ! loop that uses the filter.
    subroutine ilf_repeat_reset()
        ilf_rep_active = .false.
        ilf_rep_count = 0
        ilf_rep_lu = -1
        ilf_rep_level = ' '
        ilf_rep_source = ' '
        ilf_rep_msg = ' '
    end subroutine ilf_repeat_reset

    subroutine ilf_repeat_write(lu, level, source, msg)
        integer, intent(in) :: lu
        character(len=*), intent(in) :: level, source, msg
        if (ilf_rep_active .and. ilf_rep_lu == lu .and. &
            ilf_rep_level == level .and. ilf_rep_source == source .and. &
            ilf_rep_msg == msg) then
            ilf_rep_count = ilf_rep_count + 1
            return
        end if
        call ilf_repeat_flush(lu)
        ilf_rep_active = .true.
        ilf_rep_count = 1
        ilf_rep_lu = lu
        ilf_rep_level = level
        ilf_rep_source = source
        ilf_rep_msg = msg
    end subroutine ilf_repeat_write

    subroutine ilf_repeat_flush(lu)
        integer, intent(in) :: lu
        character(len=ilf_msglen+16) :: msg
        if (.not. ilf_rep_active) return
        if (lu /= ilf_rep_lu) return
        if (ilf_rep_count > 1) then
            write (msg, '(a,a,i0,a)') trim(ilf_rep_msg), ' (x', &
                ilf_rep_count, ')'
            call ilf_tag(lu, trim(ilf_rep_level), trim(ilf_rep_source), &
                trim(msg))
        else
            call ilf_tag(lu, trim(ilf_rep_level), trim(ilf_rep_source), &
                trim(ilf_rep_msg))
        end if
        call ilf_repeat_reset()
    end subroutine ilf_repeat_flush

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
        write (lu, '(a)') ''
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
        write (lu, '(a)') ''
    end subroutine ilf_banner

end module iter_log_fmt_mod
