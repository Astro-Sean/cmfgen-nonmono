#!/usr/bin/env python3
"""
End-of-iteration trace-species handling in this CMFGEN fork, drawn as
a flowchart.

The BA_DEM floors exempt dead trace rows from the convergence metric
while their corrections still commit.  When the majors have converged
but the exempted rows are still moving, the run does not exit: it
arms BA_TRACE_PHASE and loops for up to BA_TR_MAXSWP Stage-B sweeps
in which only the exempted rows are updated.

Box colours follow plot_iteration_cycle.py: grey = unchanged CMFGEN,
amber = existing routine with new switches, blue = new code.

Output: stage_b_flow.png (next to this script).
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch

STYLE = {
    "stock": dict(fc="#f2f2f2", ec="#777777"),
    "mod":   dict(fc="#fdf1dc", ec="#c9821e"),
    "new":   dict(fc="#e2eefb", ec="#1f6fb4"),
}

FS_TITLE = 11
FS_BODY = 9.2
LH_TITLE = 0.34
LH_BODY = 0.265
PAD_V = 0.16
TXT_PAD_L = 0.28
TXT_PAD_R = 0.25
GAP = 0.62
DPI = 160
FIG_W_IN = 12.5
AX_W_FRAC = 0.775


def box_height(lines):
    return PAD_V * 2 + LH_TITLE + LH_BODY * len(lines)


def measure_px(items):
    figm = plt.figure(figsize=(8, 3), dpi=DPI)
    axm = figm.add_axes([0, 0, 1, 1])
    artists = [axm.text(0, 0, t, fontsize=fs, fontweight=w)
               for t, fs, w in items]
    figm.canvas.draw()
    out = {a.get_text(): a.get_window_extent().width for a in artists}
    plt.close(figm)
    return out


def draw_box(ax, xc, y_top, title, lines, kind, width):
    h = box_height(lines)
    st = STYLE[kind]
    ax.add_patch(FancyBboxPatch(
        (xc - width / 2, y_top - h), width, h,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        fc=st["fc"], ec=st["ec"], lw=1.6, zorder=2))
    ax.text(xc, y_top - PAD_V - LH_TITLE / 2, title,
            ha="center", va="center", fontsize=FS_TITLE,
            fontweight="bold", zorder=3)
    y = y_top - PAD_V - LH_TITLE - LH_BODY / 2
    for ln in lines:
        ax.text(xc - width / 2 + TXT_PAD_L, y, ln,
                ha="left", va="center", fontsize=FS_BODY, zorder=3)
        y -= LH_BODY
    return y_top - h, h


def draw_diamond(ax, xc, yc, w, h, text):
    ax.add_patch(plt.Polygon(
        [(xc, yc + h / 2), (xc + w / 2, yc), (xc, yc - h / 2),
         (xc - w / 2, yc)],
        fc="#f2f2f2", ec="#777777", lw=1.6, zorder=2))
    ax.text(xc, yc, text, ha="center", va="center",
            fontsize=FS_TITLE, fontweight="bold", zorder=3)


def arrow(ax, x0, y0, x1, y1, label=None, lx=None, ly=None,
          colour="#444444"):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16,
        lw=1.4, color=colour, zorder=1))
    if label:
        ax.text(lx if lx is not None else x0 + 0.15,
                ly if ly is not None else (y0 + y1) / 2,
                label, fontsize=8.6, color=colour, va="center")


MAIN = [
    dict(kind="stock", title="Great iteration",
         lines=[
             "radiation field, populations, checkpoint",
         ]),
    dict(kind="mod", title="Measure demand   [SOLVEBA_V13]",
         lines=[
             "MAXCH: largest correction on non-exempt rows",
             "BA_DEM_EXM: largest correction the floors hid",
             "BA_TR_ABS: largest |dn|/gas on exempted rows",
         ]),
]

DEC1 = "majors\nconverged?"
SIDE = dict(kind="stock", title="Normal exit or continue",
            lines=[
                "stock logic decides:",
                "converged -> write model, return",
                "else      -> next great iteration",
            ])

ARM = dict(kind="new", title="Arm Stage-B sweep   [BA_STAGEB_ARM]",
           lines=[
               "BA_TRACE_PHASE=T; sweep count +1",
               "reassemble BA matrix (COMPUTE_BA)",
               "another great iteration runs",
           ])

SOLVE = dict(kind="new", title="Stage-B solve   [BA_TRACE_PHASE]",
             lines=[
                 "direct BA step; JFNK/LAMBDA off",
                 "only floor-exempted rows commit;",
                 "majors, ED and T restored from POPS_OLD",
             ])

DEC2 = "trace settled\nor sweeps spent?"

EXIT = dict(kind="stock", title="Exit",
            lines=[
                "converged write-out; a spent sweep",
                "budget never holds the run",
            ])

TRIGGER = [
    "Stage-B arms only when all of these hold:",
    "DO_NONMONO_SOB or DO_NONMONO_KNOP on",
    "BA_DEM_PFLOOR or BA_DEM_RELFRAC set",
    "MAXCH < EPS (majors converged)",
    "BA_TR_ABS > BA_TR_ABSMAX (traces moving)",
    "sweeps used < BA_TR_MAXSWP",
]

# ---------------------------------------------------------------
# Measure text so the canvas fits the widest box.
# ---------------------------------------------------------------
items = []
for b in MAIN + [SIDE, ARM, SOLVE, EXIT]:
    items.append((b["title"], FS_TITLE, "bold"))
    items += [(ln, FS_BODY, "normal") for ln in b["lines"]]
for t in [DEC1, DEC2] + TRIGGER:
    for ln in t.split("\n"):
        items.append((ln, FS_TITLE, "bold"))

px = measure_px(items)
aw_px = FIG_W_IN * DPI * AX_W_FRAC
wbox_px = max(px.values())

S = (0.6 + TXT_PAD_L + TXT_PAD_R + 0.3) / (1.0 - wbox_px / aw_px)
BOX_W = wbox_px * S / aw_px + TXT_PAD_L + TXT_PAD_R
# Keep the side column clear of the main column.
min_sep = BOX_W + 1.2
if S * 0.46 < min_sep:
    S = min_sep / 0.46
X_MAIN = S * 0.33
X_SIDE = S * 0.79
XLIM = S

fig_h = 17.0
fig, ax = plt.subplots(figsize=(FIG_W_IN, fig_h * 0.58))
ax.set_xlim(0, XLIM)
ax.set_ylim(0, fig_h)
ax.axis("off")

fig.suptitle("Trace-species handling at the converged exit",
             fontsize=13, fontweight="bold", y=0.985)

handles = [Patch(fc=STYLE["new"]["fc"], ec=STYLE["new"]["ec"],
                 label="new code"),
           Patch(fc=STYLE["mod"]["fc"], ec=STYLE["mod"]["ec"],
                 label="existing code, new switches"),
           Patch(fc=STYLE["stock"]["fc"], ec=STYLE["stock"]["ec"],
                 label="unchanged CMFGEN")]
fig.legend(handles=handles, loc="upper center", ncol=3,
           bbox_to_anchor=(0.5, 0.955), frameon=False, fontsize=8.8,
           handlelength=1.6)

y = fig_h - 1.7

bounds = []
for b in MAIN:
    yb, h = draw_box(ax, X_MAIN, y, b["title"], b["lines"],
                     b["kind"], BOX_W)
    bounds.append((yb, h))
    if len(bounds) > 1:
        arrow(ax, X_MAIN, bounds[-2][0], X_MAIN, y)
    y = yb - GAP

dec_h = 1.7
dec_w = BOX_W * 0.62
d1_yc = y - dec_h / 2
draw_diamond(ax, X_MAIN, d1_yc, dec_w, dec_h, DEC1)
arrow(ax, X_MAIN, bounds[-1][0], X_MAIN, d1_yc + dec_h / 2)

# Side box: normal exit/continue for the "no" branch.
yb_side, h_side = draw_box(ax, X_SIDE, y + dec_h / 2, SIDE["title"],
                           SIDE["lines"], SIDE["kind"], BOX_W)
arrow(ax, X_MAIN + dec_w / 2, d1_yc, X_SIDE - BOX_W / 2, d1_yc)
ax.text((X_MAIN + dec_w / 2 + X_SIDE - BOX_W / 2) / 2, d1_yc + 0.3,
        "no", fontsize=8.6, color="#444444", ha="center")

y = d1_yc - dec_h / 2 - GAP
yb_arm, h_arm = draw_box(ax, X_MAIN, y, ARM["title"], ARM["lines"],
                         ARM["kind"], BOX_W)
arrow(ax, X_MAIN, d1_yc - dec_h / 2, X_MAIN, y,
      label="yes + trace demand", lx=X_MAIN + 0.15)

y = yb_arm - GAP
yb_sol, h_sol = draw_box(ax, X_MAIN, y, SOLVE["title"],
                         SOLVE["lines"], SOLVE["kind"], BOX_W)
arrow(ax, X_MAIN, yb_arm, X_MAIN, y)

d2_yc = yb_sol - GAP - dec_h / 2
draw_diamond(ax, X_MAIN, d2_yc, dec_w, dec_h, DEC2)
arrow(ax, X_MAIN, yb_sol, X_MAIN, d2_yc + dec_h / 2)

y = d2_yc - dec_h / 2 - GAP
yb_exit, h_exit = draw_box(ax, X_MAIN, y, EXIT["title"],
                           EXIT["lines"], EXIT["kind"], BOX_W)
arrow(ax, X_MAIN, d2_yc - dec_h / 2, X_MAIN, y,
      label="yes", lx=X_MAIN + 0.15)

# Loop-back for "no" on the second decision: re-arm Stage-B.
x_loop = X_MAIN - BOX_W / 2 - 0.9
arm_mid = yb_arm - h_arm / 2
ax.add_patch(FancyArrowPatch(
    (X_MAIN - dec_w / 2, d2_yc), (x_loop, d2_yc),
    arrowstyle="-", lw=1.4, color="#444444"))
ax.add_patch(FancyArrowPatch(
    (x_loop, d2_yc), (x_loop, arm_mid),
    arrowstyle="-", lw=1.4, color="#444444"))
ax.add_patch(FancyArrowPatch(
    (x_loop, arm_mid), (X_MAIN - BOX_W / 2, arm_mid),
    arrowstyle="-|>", mutation_scale=16, lw=1.4, color="#444444"))
ax.text(X_MAIN - dec_w / 2 - 0.15, d2_yc + 0.25,
        "no: arm another sweep", fontsize=8.6, color="#444444",
        ha="right")

# Trigger-condition note in the right column, beside the Stage-B
# boxes.
note_x = X_SIDE - BOX_W / 2
note_y = yb_arm + 0.2
ax.text(note_x, note_y, TRIGGER[0],
        fontsize=9.0, fontweight="bold", va="top")
for ln in TRIGGER[1:]:
    note_y -= LH_BODY
    ax.text(note_x, note_y, ln, fontsize=9.0, va="top")

out = Path(__file__).with_name("stage_b_flow.png")
fig.savefig(out, dpi=DPI, bbox_inches="tight")
print(f"wrote {out}")
