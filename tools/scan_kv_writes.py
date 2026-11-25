#!/usr/bin/env python3
"""Inventory multi-field label=value WRITE/PRINT sites in CMFGEN.

Scans every .f/.f90/.for/.inc file, assembles logical statements
(fixed-form column-6 and free-form & continuations, tab-aware),
resolves FORMAT labels, and reports every output statement that emits
two or more label=value, label:value or 'label N' pairs on one output
line.

Classification (kind):
  kv-line    several labelled scalars on one line (convert)
  table-row  inside a DO loop; one row of a repeated table / loop line
  data-file  written to a unit opened on a data/restart/model file
  parsed     emits a token an external script greps (tokens file)
  namelist   NML= write; already vertical key=value output
  unclear    format or unit could not be resolved statically

CSV columns: file,routine,line,unit,header_text,n_fields,est_width,
kind,in_loop,conditional,fmt_ref,notes

Usage: scan_kv_writes.py [root] [--tokens parser_tokens.txt]
                              [--csv out.csv]
"""

import csv
import re
import sys
from collections import Counter
from pathlib import Path

FIXED_EXTS = {".f", ".for"}
FREE_EXTS = {".f90"}
SKIP_DIRS = {".git", "cmfgen_obsidian", "docs", "web"}

# Opened-file names that are data products read back by code or
# scripts (model/restart, machine-readable tables, unformatted dumps).
# Writes to these are never converted.
DATA_FILE_RE = re.compile(
    r"^(RVTJ.*|MODEL.*|MODEL_SCR|SCRTEMP|SIGE.*|RVSIGE.*|STEQ_VALS.*|"
    r"OBSFLUX.*|OBS.*|MEANOPAC.*|EWDATA.*|NETRATE.*|TOTRATE.*|"
    r"LINEHEAT.*|NEG_OPAC.*|MOD_SUM.*|NON_THERM.*|GREY_SCL.*|"
    r"IMPURITYJ|ES_J_CONV|JEW|IN_ITS|OUT_DIRECT|EDDCONT|EDDFACTOR|"
    r"BA.*DATA|STORE_BA.*|DD_.*|NUC_.*|GAMMA.*|SPEC_.*|FLUX.*|"
    r"POPDAT|POPUL.*|HYDRO.*|VEL.*DATA|DEPTH.*|GRID.*|XRAY.*|"
    r".*\.(DAT|BIN|NPY|NPZ|H5|HDF|NC))$",
    re.IGNORECASE,
)

DESC_NUM = re.compile(r"^(I|F|E|ES|EN|D|G|B|O|Z)(\d+)(\.\d+)?(E\d+)?$",
                      re.IGNORECASE)
DESC_A = re.compile(r"^A(\d+)?$", re.IGNORECASE)
DESC_L = re.compile(r"^L(\d+)$", re.IGNORECASE)
DESC_CTRL = re.compile(r"^(\d+X|TR?\d*|TL\d+|\d*P|S|SP|SS|BN|BZ|Q)$",
                       re.IGNORECASE)
LABELISH_RE = re.compile(r"[A-Za-z0-9]")

ROUTINE_RE = re.compile(
    r"^\s*(?:[\w*(),=\s:]*?\s)?(subroutine|program|module)\s+(\w+)",
    re.IGNORECASE)
FUNC_RE = re.compile(r"\bfunction\s+(\w+)", re.IGNORECASE)
END_RE = re.compile(r"^\s*end(\s+(program|subroutine|function|module|"
                    r"block\s+data)\b.*)?$", re.IGNORECASE)


def strip_bang(text):
    """Remove a trailing ! comment outside quotes."""
    out = []
    q = None
    for c in text:
        if q is not None:
            out.append(c)
            if c == q:
                q = None
        elif c in "'\"":
            q = c
            out.append(c)
        elif c == "!":
            break
        else:
            out.append(c)
    return "".join(out)


def is_fixed_form(path, lines):
    ext = path.suffix.lower()
    if ext in FREE_EXTS:
        return False
    if ext in FIXED_EXTS:
        return True
    for raw in lines[:300]:
        if strip_bang(raw).rstrip().endswith("&"):
            return False
    return True


def fixed_fields(raw):
    """Split a fixed-form line -> (label, cont, code) or None.

    A tab in the leading whitespace jumps to column 7; a digit right
    after such a tab is the continuation indicator (ifort fixed-form
    rule).
    """
    if not raw.strip():
        return None
    if raw[0] in "Cc*!Dd#":
        return None
    # A tab in columns 1-6 jumps to column 7 (ifort fixed-form rule):
    # '3\tFORMAT' is a labelled statement, '\t1 X' marks X as a
    # continuation, and a non-blank at column 6 before the tab is the
    # continuation indicator.
    tidx = raw.find("\t", 0, 7)
    if tidx == 0:
        rest = raw[1:]
        if re.match(r"^\d(\s|$)", rest):
            return "", True, strip_bang(rest[1:]).rstrip()
        lm = re.match(r"^(\d+)", rest)
        if lm:
            return lm.group(1), False, \
                strip_bang(rest[lm.end():]).rstrip()
        return "", False, strip_bang(rest).rstrip()
    if 0 < tidx <= 6:
        before = raw[:tidx]
        after = strip_bang(raw[tidx + 1:]).rstrip()
        if len(before) > 5 and before[5] not in (" ", "0"):
            return before[:5].strip(), True, after
        return before.strip(), False, after
    cont = len(raw) > 5 and raw[5] not in (" ", "0")
    lab = raw[:5].strip()
    code = strip_bang(raw[6:72]).rstrip() if len(raw) > 6 else ""
    return lab, cont, code


def assemble(lines, fixed):
    """Return list of (lineno, label, text) logical statements."""
    stmts = []
    cur = None
    start = 0
    label = None
    if fixed:
        for i, raw in enumerate(lines, 1):
            parsed = fixed_fields(raw)
            if parsed is None:
                continue
            lab, cont, code = parsed
            if cont and cur is not None:
                cur += " " + code.strip()
                continue
            if cur is not None:
                stmts.append((start, label, cur.strip()))
            cur = code
            start = i
            label = lab if lab.isdigit() else None
    else:
        for i, raw in enumerate(lines, 1):
            code = strip_bang(raw).strip()
            if not code:
                continue
            lab = None
            m = re.match(r"^(\d+)\s+(.+)", code)
            if m and not code.startswith("&"):
                lab, code = m.group(1), m.group(2)
            if code.startswith("&"):
                code = code[1:].lstrip()
            cont_next = code.rstrip().endswith("&")
            if cont_next:
                code = code.rstrip()[:-1].rstrip()
            if cur is None:
                cur, start = "", i
            if code:
                cur += (" " if cur else "") + code
            if lab and label is None:
                label = lab
            if not cont_next:
                stmts.append((start, label, cur.strip()))
                cur = None
                label = None
    if cur is not None:
        stmts.append((start, label, cur.strip()))
    return stmts


def split_top(s, sep=","):
    """Split on top-level commas, respecting parens and quotes."""
    parts, depth, q, last = [], 0, None, 0
    for i, c in enumerate(s):
        if q is not None:
            if c == q:
                q = None
        elif c in "'\"":
            q = c
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == sep and depth == 0:
            parts.append(s[last:i])
            last = i + 1
    parts.append(s[last:])
    return parts


def close_paren(s, i):
    """Index of the ')' matching '(' at s[i], respecting quotes."""
    depth, q = 0, None
    for j in range(i, len(s)):
        c = s[j]
        if q is not None:
            if c == q:
                q = None
        elif c in "'\"":
            q = c
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return j
    return -1


def unquote(s):
    s = s.strip()
    if len(s) >= 2 and s[0] in "'\"" and s[-1] == s[0]:
        return s[1:-1].replace(s[0] * 2, s[0])
    return None


def parse_ctl(ctl):
    """Parse a WRITE control list -> dict(unit, fmt, rec, nml)."""
    out = {"unit": None, "fmt": None, "rec": False, "nml": False}
    pos = 0
    for item in split_top(ctl):
        item = item.strip()
        m = re.match(r"^(\w+)\s*=\s*(.*)$", item, re.DOTALL)
        key = m.group(1).upper() if m else None
        val = m.group(2).strip() if m else item
        if key == "UNIT":
            out["unit"] = val
        elif key == "FMT":
            out["fmt"] = val
        elif key == "NML":
            out["nml"] = True
        elif key == "REC":
            out["rec"] = True
        elif key is None:
            if pos == 0:
                out["unit"] = val
            elif pos == 1:
                out["fmt"] = val
            pos += 1
    return out


def parse_format_body(body):
    """Tokenise a format body into (kind, payload) tokens.

    kind in {'lit','desc','ctrl','nl'}; parenthesised repeat groups are
    expanded (cap 99).
    """
    toks = []
    n = len(body)

    def parse_group(j):
        """Body[j]=='('; return (toks, index after matching ')')."""
        depth, q = 0, None
        k = j + 1
        while k < n:
            c = body[k]
            if q is not None:
                if c == q:
                    q = None
            elif c in "'\"":
                q = c
            elif c == "(":
                depth += 1
            elif c == ")":
                if depth == 0:
                    break
                depth -= 1
            k += 1
        return parse_format_body(body[j + 1:k]), k + 1

    i = 0
    while i < n:
        c = body[i]
        if c in " ,\t":
            i += 1
            continue
        if c == "/":
            toks.append(("nl", "/"))
            i += 1
            continue
        if c == ":":
            toks.append(("ctrl", ":"))
            i += 1
            continue
        if c in "'\"":
            j = i + 1
            txt = []
            while j < n:
                if body[j] == c:
                    if j + 1 < n and body[j + 1] == c:
                        txt.append(c)
                        j += 2
                        continue
                    break
                txt.append(body[j])
                j += 1
            toks.append(("lit", "".join(txt)))
            i = j + 1
            continue
        m = re.match(r"^(\d+)[Hh]", body[i:])
        if m:
            cnt = int(m.group(1))
            s = i + m.end()
            toks.append(("lit", body[s:s + cnt]))
            i = s + cnt
            continue
        m = re.match(r"^(\d+)\s*\(", body[i:])
        if m:
            rep = min(int(m.group(1)), 99)
            inner, i = parse_group(i + m.end() - 1)
            for _ in range(rep):
                toks.extend(inner)
            continue
        m = re.match(r"^[A-Za-z]+[0-9]*(\.[0-9]+)?([Ee][0-9]+)?",
                     body[i:])
        if m:
            d = m.group(0).upper()
            if DESC_CTRL.match(d) or d == "X":
                toks.append(("ctrl", d))
            else:
                toks.append(("desc", d))
            i += m.end()
            continue
        toks.append(("ctrl", c))
        i += 1
    return toks


def desc_width(d):
    m = DESC_NUM.match(d)
    if m:
        return int(m.group(2))
    m = DESC_A.match(d)
    if m:
        return int(m.group(1)) if m.group(1) else 8
    m = DESC_L.match(d)
    if m:
        return int(m.group(1))
    m = re.match(r"^(\d+)X$", d)
    if m:
        return int(m.group(1))
    return 0


def labelish(text):
    return bool(LABELISH_RE.search(text))


def header_of(first_lit):
    """Header = leading text of the first literal up to its first
    ':' or '=' (the part that is not itself a field label)."""
    if not first_lit:
        return ""
    m = re.search(r"[:=]", first_lit)
    if m:
        head = first_lit[:m.end()]
        return head if labelish(head) else ""
    return ""


def fields_in_toks(toks, n_lit_items):
    """Count label->value fields on the busiest /-separated line.

    A literal followed by one or more value descriptors is one field.
    Literal I/O items also count when the format has A descriptors.
    Returns (max_fields_per_line, first_lit, est_width).
    """
    per_line = []
    cur = 0
    pending = None
    consumed = False
    width = 0
    first_lit = ""
    seen_lit = False
    has_a = False
    for kind, pay in toks:
        if kind == "nl":
            per_line.append(cur)
            cur = 0
            pending = None
            consumed = False
            continue
        if kind == "lit":
            if not seen_lit:
                first_lit = pay
                seen_lit = True
            pending = pay
            consumed = False
            width += len(pay)
            continue
        if kind == "desc":
            width += desc_width(pay)
            if DESC_A.match(pay):
                has_a = True
            if pending is not None and not consumed and \
                    labelish(pending):
                cur += 1
                consumed = True
            continue
    per_line.append(cur)
    best = max(per_line) if per_line else 0
    if has_a and n_lit_items > best:
        best = n_lit_items
    return best, first_lit, width


def parse_io_items(list_txt):
    """Split an I/O list into ('lit',text)/('val',raw)/('ido',raw)."""
    items = []
    for part in split_top(list_txt):
        p = part.strip()
        if not p:
            continue
        if p.startswith("(") and re.search(r",\s*\w+\s*=", p):
            items.append(("ido", p))
            continue
        q = unquote(p)
        if q is not None:
            items.append(("lit", q))
            continue
        lits = [lm.group(1) if lm.group(1) is not None else lm.group(2)
                for lm in re.finditer(
                    r"'((?:[^']|'')*)'|\"((?:[^\"]|\"\")*)\"", p)]
        if lits:
            items.append(("lit", "".join(lits)))
        else:
            items.append(("val", p))
    return items


def fields_in_items(items):
    """Alternating literal->value rule for list-directed output."""
    per_line = []
    cur = 0
    pending = None
    consumed = False
    width = 0
    first_lit = ""
    seen_lit = False
    for kind, pay in items:
        if kind == "lit":
            if not seen_lit:
                first_lit = pay
                seen_lit = True
            pending = pay
            consumed = False
            width += len(pay)
            continue
        if kind == "val":
            width += 13
            if pending is not None and not consumed and \
                    labelish(pending):
                cur += 1
                consumed = True
            continue
        if kind == "ido":
            width += 13
            if pending is not None and not consumed and \
                    labelish(pending):
                cur += 1
                consumed = True
    per_line.append(cur)
    return max(per_line) if per_line else 0, first_lit, width


def collect_opens(stmts):
    """Map unit literals / var names -> opened file info."""
    lit_map = {}
    var_map = {}
    for _, _, st in stmts:
        m = re.match(r"(?i)call\s+gen_asci_open\s*\(", st)
        if m:
            j = close_paren(st, m.end() - 1)
            args = split_top(st[m.end():j]) if j > 0 else []
            if len(args) >= 2:
                unit = args[0].strip()
                name = unquote(args[1]) or args[1].strip()
                if unit.isdigit():
                    lit_map.setdefault(unit.lstrip("0") or "0",
                                       []).append(name.upper())
                elif re.match(r"^\w+$", unit):
                    var_map.setdefault(unit.upper(), set()).add(
                        name.upper())
            continue
        m = re.match(r"(?i)open\s*\(", st)
        if m:
            j = close_paren(st, m.end() - 1)
            if j < 0:
                continue
            unit = None
            name = None
            unformatted = False
            first = True
            for item in split_top(st[m.end():j]):
                item = item.strip()
                kv = re.match(r"^(\w+)\s*=\s*(.*)$", item, re.DOTALL)
                if kv:
                    k = kv.group(1).upper()
                    v = kv.group(2).strip()
                    if k == "UNIT":
                        unit = v
                    elif k == "FILE":
                        name = unquote(v) or v
                    elif k == "FORM" and "UNFORM" in v.upper():
                        unformatted = True
                elif first and unit is None:
                    unit = item
                first = False
            if unit:
                nm = ("U:" if unformatted else "") + (name or "?").upper()
                if unit.isdigit():
                    lit_map.setdefault(unit.lstrip("0") or "0",
                                       []).append(nm)
                elif re.match(r"^\w+$", unit):
                    var_map.setdefault(unit.upper(), set()).add(nm)
    return lit_map, var_map


def unit_kind(unit, assigns, lit_map, var_map):
    """Resolve a unit expr -> (unit_desc, kind_hint, note)."""
    u = unit.strip().upper()
    if u in ("*", "6", "LUER", "LU_ER", "LUERR", "OUTPUT_UNIT",
             "STDOUT_LU"):
        return "stdout", "console", ""
    if u in ("0", "ERROR_UNIT"):
        return "stderr", "console", ""
    if u.isdigit():
        names = lit_map.get(u.lstrip("0") or "0")
        if names:
            nm = names[0]
            if nm.startswith("U:"):
                return "file:" + nm[2:], "data", ""
            return "file:" + nm, "file", ""
        return "unit:" + u, "unknown", "unopened-literal"
    if re.match(r"^\w+$", u):
        seen = set()
        v = u
        while v in assigns and v not in seen:
            seen.add(v)
            v = assigns[v].upper()
        if v != u:
            if v.isdigit():
                return unit_kind(v, assigns, lit_map, var_map)
            if v in ("ERROR_LU()", "OUTPUT_UNIT()"):
                return "stdout", "console", "via " + v
            if v in ("LUER", "LU_ER", "LUERR"):
                return "stdout", "console", "via " + v
            if v in var_map:
                names = var_map[v]
                if len(names) == 1:
                    nm = next(iter(names))
                    if nm.startswith("U:"):
                        return "file:" + nm[2:], "data", ""
                    return "file:" + nm, "file", ""
                return "var:" + u, "unknown", \
                    "multi-file:" + ";".join(sorted(names))
            return "var:" + u, "unknown", "unresolved-expr:" + v
        names = var_map.get(u)
        if names and len(names) == 1:
            nm = next(iter(names))
            if nm.startswith("U:"):
                return "file:" + nm[2:], "data", ""
            return "file:" + nm, "file", ""
        if names:
            return "var:" + u, "unknown", \
                "multi-file:" + ";".join(sorted(names))
        return "var:" + u, "unknown", "unresolved-var"
    return "expr:" + u[:20], "unknown", "unresolved-expr"


def is_data_unit(udesc):
    return udesc.startswith("file:") and bool(
        DATA_FILE_RE.match(udesc[5:]))


def scan_file(path, tokens):
    try:
        # split on \n only: splitlines() also breaks on form feeds,
        # which old Fortran sources use as page separators and which
        # would desync reported line numbers from the file's own.
        lines = [l.replace("\x0c", " ").replace("\x0b", " ")
                 for l in path.read_text(errors="replace").split("\n")]
    except OSError:
        return [], {}, {}
    fixed = is_fixed_form(path, lines)
    stmts = assemble(lines, fixed)
    lit_map, var_map = collect_opens(stmts)

    fmt_map = {}
    for _, lab, st in stmts:
        if lab:
            fm = re.match(r"(?i)^format\s*\(", st)
            if fm:
                j = close_paren(st, fm.end() - 1)
                if j > 0:
                    fmt_map[lab] = st[fm.end():j]

    assigns = {}
    char_vars = set()
    char_assigns = {}
    for _, _, st in stmts:
        m = re.match(
            r"^(\w+)\s*=\s*(ERROR_LU\(\)|GET_LU\(\)|OUT_LU\(\)|\w+|\d+)"
            r"\s*$", st, re.IGNORECASE)
        if m:
            assigns[m.group(1).upper()] = m.group(2).upper()
            continue
        # integer constant declarations: INTEGER, PARAMETER :: T_OUT=6
        if re.match(r"(?i)^\s*integer\b", st):
            for km in re.finditer(r"(?i)\b(\w+)\s*=\s*(\d+)\b", st):
                assigns.setdefault(km.group(1).upper(),
                                   km.group(2))
        # CHARACTER declarations: a WRITE to a character variable is
        # an internal-file string build, not console output
        cm = re.match(r"(?i)^\s*character\b[^=]*?(::)?\s*(.*)$", st)
        if cm:
            for nm in re.finditer(r"(?i)\b(\w+)\s*(?:\(|\*|\d|,|$)", st):
                w = nm.group(1).upper()
                if w not in ("CHARACTER", "LEN", "KIND", "FUNCTION"):
                    char_vars.add(w)
        # character fmt variables: FMT = '(A,I5)' possibly // joined
        fm = re.match(r"(?i)^(\w+)\s*=\s*(('([^']|'')*'|\"([^\"]|\"\")*\")"
                      r"(\s*//\s*('([^']|'')*'|\"([^\"]|\"\")*\"))*)"
                      r"\s*$", st)
        if fm:
            lits = re.findall(
                r"'((?:[^']|'')*)'|\"((?:[^\"]|\"\")*)\"",
                fm.group(2))
            val = "".join(a or b for a, b in lits)
            char_assigns.setdefault(fm.group(1).upper(), set()).add(val)

    routine = "(outside)"
    loop_stack = []
    rows = []
    char_frag = {}
    for lineno, lab, st in stmts:
        rm = ROUTINE_RE.match(st)
        fm = None if rm else FUNC_RE.search(st)
        if rm:
            routine = rm.group(2)
        elif fm and not END_RE.match(st):
            routine = fm.group(1)
        elif END_RE.match(st):
            routine = "(outside)"

        work = st
        cond = False
        im = re.match(r"(?i)^if\s*\(", work)
        if im:
            j = close_paren(work, work.find("("))
            if j > 0 and not re.match(
                    r"(?i)^\s*then\b", work[j + 1:]):
                cond = True
                work = work[j + 1:].strip()

        dm = re.match(r"(?i)^do\s+(\d+)\s+\w+\s*=", work)
        dbm = re.match(r"(?i)^do\b\s*(\w+\s*=|while|$)", work)
        if dm:
            loop_stack.append(("label", dm.group(1)))
        elif dbm:
            loop_stack.append(("block", None))
        elif re.match(r"(?i)^end\s*do\b|^enddo\b", work):
            for t in range(len(loop_stack) - 1, -1, -1):
                if loop_stack[t][0] == "block":
                    loop_stack.pop(t)
                    break

        in_loop = bool(loop_stack)

        # collect string fragments assigned to variables so parser
        # tokens carried by e.g. TRIM(STRING) still flag the WRITE
        am = re.match(r"(?i)^(\w+)\s*=\s*(.*)$", work)
        if am and ("'" in am.group(2) or '"' in am.group(2)):
            lits = [lm.group(0) for lm in re.finditer(
                r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"", am.group(2))]
            if lits:
                frag = " ".join(unquote(l) or l[1:-1] for l in lits)
                rk = routine if routine != "(outside)" else ""
                key = (rk, am.group(1).upper())
                cur = char_frag.get(key, "")
                if len(cur) < 400:
                    char_frag[key] = (cur + " " + frag)[:400]

        wm = re.match(r"(?i)^write\s*\(", work)
        pm = re.match(r"(?i)^print\b(.*)", work)
        if not (wm or pm):
            if lab:
                while loop_stack and loop_stack[-1][0] == "label" and \
                        loop_stack[-1][1] == lab:
                    loop_stack.pop()
            continue

        ctl = None
        io_list = ""
        if wm:
            j = close_paren(work, work.find("("))
            if j < 0:
                if lab:
                    while loop_stack and \
                            loop_stack[-1][0] == "label" and \
                            loop_stack[-1][1] == lab:
                        loop_stack.pop()
                continue
            ctl = parse_ctl(work[work.find("(") + 1:j])
            io_list = work[j + 1:]
        else:
            rest = pm.group(1).strip()
            if not rest:
                continue
            fmt_spec, _, io_list = rest.partition(",")
            if fmt_spec.strip().startswith("("):
                jj = close_paren(rest, 0)
                fmt_spec = rest[:jj + 1]
                io_list = rest[jj + 1:].lstrip(",")
            elif fmt_spec.strip().startswith(("'", '"')):
                qq = unquote(rest)
                if qq is not None:
                    fmt_spec, io_list = rest, ""
                else:
                    fmt_spec, _, io_list = rest.partition(",")
            ctl = {"unit": "*", "fmt": fmt_spec.strip(), "rec": False,
                   "nml": False}

        if lab:
            while loop_stack and loop_stack[-1][0] == "label" and \
                    loop_stack[-1][1] == lab:
                loop_stack.pop()

        base = {"file": str(path), "routine": routine, "line": lineno,
                "in_loop": in_loop, "cond": cond,
                "src": re.sub(r"\s+", " ", work)[:140]}

        # internal-file write into a CHARACTER variable: a string
        # build, not console/log output
        if (ctl.get("unit") or "").strip().upper() in char_vars:
            continue

        if ctl.get("nml"):
            unit_info = unit_kind(ctl["unit"] or "*", assigns,
                                  lit_map, var_map)
            rows.append({**base, "unit": unit_info[0],
                         "header": "(namelist)", "n_fields": -1,
                         "width": -1, "kind": "namelist",
                         "fmt": "NML", "notes": unit_info[2]})
            continue

        if ctl.get("rec"):
            unit_info = unit_kind(ctl["unit"] or "*", assigns,
                                  lit_map, var_map)
            rows.append({**base, "unit": unit_info[0],
                         "header": "(direct-access)", "n_fields": -1,
                         "width": -1, "kind": "data-file",
                         "fmt": "-", "notes": "direct-access;" +
                         unit_info[2]})
            continue

        fmt = ctl.get("fmt")
        n_fields = 0
        width = 0
        first_lit = ""
        note = ""
        fmt_ref = ""
        unresolved = False
        io_items = parse_io_items(io_list)
        n_lit_items = sum(1 for k, p in io_items
                          if k == "lit" and labelish(p))
        if fmt is None:
            if io_items:
                unresolved = True
                note = "no-fmt-with-items"
            else:
                continue
        elif fmt.strip() == "*":
            n_fields, first_lit, width = fields_in_items(io_items)
            fmt_ref = "*"
        else:
            q = unquote(fmt)
            body = None
            if q is not None:
                body = q
                fmt_ref = "inline"
            elif fmt.strip().isdigit():
                fmt_ref = "FMT " + fmt.strip()
                if fmt.strip().lstrip("0") in fmt_map or \
                        fmt.strip() in fmt_map:
                    body = fmt_map.get(fmt.strip(),
                                       fmt_map[fmt.strip().lstrip("0")])
                else:
                    unresolved = True
                    note = "fmt-label-missing:" + fmt.strip()
            elif fmt.strip().startswith("("):
                body = unquote(fmt.strip()) or fmt.strip()[1:-1]
                fmt_ref = "inline"
            elif fmt.strip().upper() in assigns and \
                    assigns[fmt.strip().upper()].isdigit():
                fmt_ref = "FMT " + assigns[fmt.strip().upper()]
                lab = assigns[fmt.strip().upper()].lstrip("0") or "0"
                if lab in fmt_map:
                    body = fmt_map[lab]
                else:
                    unresolved = True
                    note = "fmt-label-missing:" + lab
            elif fmt.strip().upper() in char_assigns and \
                    len(char_assigns[fmt.strip().upper()]) == 1:
                body = next(iter(char_assigns[fmt.strip().upper()]))
                fmt_ref = "fmt-char:" + fmt.strip()[:14]
            else:
                unresolved = True
                note = "fmt-var:" + fmt.strip()[:20]
                fmt_ref = fmt.strip()[:20]
            if body is not None:
                b = body.strip()
                if b.startswith("(") and \
                        close_paren(b, 0) == len(b) - 1:
                    b = b[1:-1]
                toks = parse_format_body(b)
                n_fields, first_lit, width = fields_in_toks(
                    toks, n_lit_items)

        if n_fields < 2 and not unresolved:
            continue

        unit_info = unit_kind(ctl["unit"] or "*", assigns, lit_map,
                              var_map)
        hit = ""
        if tokens:
            probe_parts = [first_lit] + \
                [p for k, p in io_items if k == "lit"]
            for k, p in io_items:
                if k != "lit":
                    for vn in re.findall(r"\b(\w+)\b", p):
                        probe_parts.append(
                            char_frag.get((routine, vn.upper()), ""))
                        probe_parts.append(
                            char_frag.get(("", vn.upper()), ""))
            probe = " ".join(probe_parts)
            src_probe = re.sub(r"\s+", " ", work)
            for t in tokens:
                if t.lower() in probe.lower() or \
                        t.lower() in src_probe.lower():
                    hit = t
                    break

        udesc, uk, unote = unit_info
        if unresolved and note == "no-fmt-with-items" and \
                uk in ("file", "data"):
            kind = "data-file"
        elif uk == "data" or is_data_unit(udesc):
            kind = "data-file"
        elif hit:
            kind = "parsed"
            unote = (unote + ";" if unote else "") + "tok:" + hit
        elif unresolved:
            kind = "unclear"
            unote = (note + ";" + unote) if unote else note
        elif uk == "unknown":
            kind = "unclear"
        elif in_loop:
            kind = "table-row"
        else:
            kind = "kv-line"

        rows.append({**base, "unit": unit_info[0],
                     "header": header_of(first_lit),
                     "n_fields": n_fields, "width": width,
                     "kind": kind, "fmt": fmt_ref, "notes": unote})
    return rows, lit_map, var_map, assigns


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 and
                not sys.argv[1].startswith("--") else ".")
    csv_path = "kv_writes_inventory.csv"
    tok_path = None
    args = sys.argv[1:]
    for i, a in enumerate(args):
        if a == "--csv" and i + 1 < len(args):
            csv_path = args[i + 1]
        if a == "--tokens" and i + 1 < len(args):
            tok_path = args[i + 1]
    tokens = []
    if tok_path:
        tokens = [t.strip() for t in
                  Path(tok_path).read_text().splitlines()
                  if t.strip() and not t.startswith("#")]

    var_map_g = {}
    lit_map_g = {}
    all_rows = []
    files = sorted(p for p in root.rglob("*")
                   if p.suffix.lower() in
                   (FIXED_EXTS | FREE_EXTS | {".inc"})
                   and not any(d in SKIP_DIRS for d in p.parts))
    for p in files:
        rows, lit_map, var_map, _ = scan_file(p, tokens)
        for k, v in var_map.items():
            var_map_g.setdefault(k, set()).update(v)
        for k, v in lit_map.items():
            lit_map_g.setdefault(k, set()).update(v)
        all_rows.extend(rows)

    # second pass: resolve var:/unit: units against the global OPEN
    # maps and re-classify
    def reclass(r):
        if r["kind"] == "unclear" and \
                "unopened-literal" in r["notes"]:
            r["kind"] = "table-row" if r["in_loop"] else "kv-line"
            r["notes"] = r["notes"].replace("unopened-literal", "").\
                strip(";")
        if is_data_unit(r["unit"]) and r["kind"] != "parsed":
            r["kind"] = "data-file"

    for r in all_rows:
        if r["unit"].startswith("var:"):
            name = r["unit"][4:]
            names = var_map_g.get(name.upper())
            if names and len(names) == 1:
                nm = next(iter(names))
                if nm.startswith("U:"):
                    r["unit"] = "file:" + nm[2:]
                    r["kind"] = "data-file"
                else:
                    r["unit"] = "file:" + nm
                    reclass(r)
            elif names:
                if all(n.startswith("U:") or
                       DATA_FILE_RE.match(n) for n in names):
                    r["unit"] = "file:" + sorted(names)[0].lstrip("U:")
                    r["kind"] = "data-file"
                else:
                    r["notes"] = (r["notes"] + ";" if r["notes"]
                                  else "") + \
                        "multi-file:" + ";".join(sorted(names))
        elif r["unit"].startswith("unit:"):
            num = r["unit"][5:]
            names = lit_map_g.get(num)
            if names and len(names) == 1:
                nm = next(iter(names))
                if nm.startswith("U:"):
                    r["unit"] = "file:" + nm[2:]
                    r["kind"] = "data-file"
                else:
                    r["unit"] = "file:" + nm
                    reclass(r)
            elif names:
                if all(n.startswith("U:") or
                       DATA_FILE_RE.match(n) for n in names):
                    r["unit"] = "file:" + sorted(names)[0].lstrip("U:")
                    r["kind"] = "data-file"
                else:
                    r["notes"] = (r["notes"] + ";" if r["notes"]
                                  else "") + \
                        "multi-open:" + ";".join(sorted(names))

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "routine", "line", "unit", "header_text",
                    "n_fields", "est_width", "kind", "in_loop",
                    "conditional", "fmt_ref", "notes", "stmt"])
        for r in sorted(all_rows, key=lambda x: (x["file"], x["line"])):
            w.writerow([r["file"], r["routine"], r["line"], r["unit"],
                        r["header"], r["n_fields"], r["width"],
                        r["kind"], r["in_loop"], r["cond"], r["fmt"],
                        r["notes"], r["src"]])

    by_kind = Counter(r["kind"] for r in all_rows)
    by_dir = Counter()
    for r in all_rows:
        rel = Path(r["file"]).resolve().relative_to(root.resolve())
        by_dir[(rel.parts[0] if len(rel.parts) > 1 else ".",
                r["kind"])] += 1
    print(f"files scanned: {len(files)}")
    print(f"sites found:   {len(all_rows)}")
    print("\nby kind:")
    for k, c in by_kind.most_common():
        print(f"  {k:10s} {c}")
    print("\nby directory x kind:")
    kinds = [k for k, _ in by_kind.most_common()]
    print("  dir".ljust(22) + "".join(k.rjust(11) for k in kinds))
    for d in sorted({d for d, _ in by_dir}):
        line = "  " + d.ljust(20)
        for k in kinds:
            line += str(by_dir.get((d, k), 0)).rjust(11)
        print(line)
    print("\ntop 20 sites by field count:")
    for r in sorted(all_rows, key=lambda x: -x["n_fields"])[:20]:
        rel = Path(r["file"]).resolve().relative_to(root.resolve())
        print(f"  {r['n_fields']:3d}  {rel}:{r['line']:5d} "
              f"{r['routine'][:20]:20s} {r['kind']:9s} "
              f"{r['src'][:70]}")


if __name__ == "__main__":
    main()
