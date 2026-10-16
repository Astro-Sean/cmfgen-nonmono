#!/usr/bin/env python3
"""
End-of-iteration trace-species handling in this CMFGEN fork, drawn as
a flowchart.

Demand floors exempt dead trace rows from the convergence metric
while their corrections still commit.  When the major species have
converged but the exempted rows are still moving, the run does not
exit: it arms a trace-polish phase and loops for a bounded number of
sweeps in which only the exempted rows are updated.

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
TXT_PAD_L = 0.30
TXT_PAD_R = 0.30
GAP = 0.45
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
             "Radiation field, populations, checkpoint",
         ]),
    dict(kind="mod", title="Measure the demand",
         lines=[
             "Largest correction on ordinary rows",
             "Largest correction the floors hid",
             "Largest trace change per gas density",
         ]),
]

DEC1 = "Majors\nconverged?"
SIDE = dict(kind="stock", title="Normal exit or continue",
            lines=[
                "Stock logic decides:",
                "Converged -> write model, return",
                "Else      -> next great iteration",
            ])

ARM = dict(kind="new", title="Arm a polish sweep",
           lines=[
               "Mark the next iteration as a polish sweep;",
               "Sweep count +1",
               "Rebuild the update matrix",
               "Another great iteration runs",
           ])

SOLVE = dict(kind="new", title="Polish solve",
             lines=[
                 "Direct matrix step; Krylov and",
                 "fixed-point solvers stay off",
                 "Only floor-exempted rows commit;",
                 "Majors, density and temperature held",
             ])

DEC2 = "Trace settled\nor sweeps spent?"

EXIT = dict(kind="stock", title="Exit",
            lines=[
                "Converged write-out; a spent sweep",
                "Budget never holds the run",
            ])

TRIGGER = [
    "A polish sweep runs only when all of these hold:",
    "A non-monotonic solver is active",
    "A demand floor is set",
    "Ordinary rows have converged",
    "Traces are still moving",
    "Sweeps remain in the budget",
]

# ---------------------------------------------------------------
# Measure text so the canvas fits the widest box.
# ---------------------------------------------------------------
# Only box text sizes the boxes; notes and diamonds have their own
# layout.  TRIGGER lines render at 9 pt, so they are measured at 9 pt.
items = []
for b in MAIN + [SIDE, ARM, SOLVE, EXIT]:
    items.append((b["title"], FS_TITLE, "bold"))
    items += [(ln, FS_BODY, "normal") for ln in b["lines"]]

px = measure_px(items)
px.update(measure_px(
    [(ln, FS_TITLE, "bold") for t in [DEC1, DEC2]
     for ln in t.split("\n")]
    + [(TRIGGER[0], 9.0, "bold")]
    + [(ln, 9.0, "normal") for ln in TRIGGER[1:]]))

aw_px = FIG_W_IN * DPI * AX_W_FRAC
wbox_px = max(px[ln] for b in MAIN + [SIDE, ARM, SOLVE, EXIT]
              for ln in [b["title"]] + b["lines"])

# Layout: main column at X_MAIN, note column at X_SIDE.  BOX_W grows
# with XLIM (pixels convert at XLIM/aw_px) and XLIM must contain
# BOX_W twice plus the column gap and margins, so solve directly.
COL_GAP = 1.0
A = wbox_px * 1.06 / aw_px
B = TXT_PAD_L + TXT_PAD_R
# XLIM = 2*BOX_W + COL_GAP + left loop margin (1.1) + right margin (0.4)
XLIM = (2 * B + COL_GAP + 1.5) / (1.0 - 2.0 * A)
BOX_W = A * XLIM + B
X_MAIN = BOX_W / 2 + 1.1
X_SIDE = X_MAIN + BOX_W + COL_GAP

fig_h = 14.6
fig, ax = plt.subplots(figsize=(FIG_W_IN, fig_h * 0.58))
ax.set_xlim(0, XLIM)
ax.set_ylim(0, fig_h)
ax.axis("off")

fig.suptitle("Trace-species polish at the converged exit",
             fontsize=13, fontweight="bold", y=0.985)

handles = [Patch(fc=STYLE["new"]["fc"], ec=STYLE["new"]["ec"],
                 label="New code"),
           Patch(fc=STYLE["mod"]["fc"], ec=STYLE["mod"]["ec"],
                 label="Existing code, new switches"),
           Patch(fc=STYLE["stock"]["fc"], ec=STYLE["stock"]["ec"],
                 label="Unchanged CMFGEN")]
fig.legend(handles=handles, loc="upper center", ncol=3,
           bbox_to_anchor=(0.5, 0.955), frameon=False, fontsize=8.8,
           handlelength=1.6)

y = fig_h - 1.0

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
        "No", fontsize=8.6, color="#444444", ha="center")

y = d1_yc - dec_h / 2 - GAP
yb_arm, h_arm = draw_box(ax, X_MAIN, y, ARM["title"], ARM["lines"],
                         ARM["kind"], BOX_W)
arrow(ax, X_MAIN, d1_yc - dec_h / 2, X_MAIN, y,
      label="Yes + trace demand", lx=X_MAIN + 0.15)

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
      label="Yes", lx=X_MAIN + 0.15)

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
        "No: arm another sweep", fontsize=8.6, color="#444444",
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
