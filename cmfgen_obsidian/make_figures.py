#!/usr/bin/env python3
"""Regenerate the schematic figures embedded in the obsidian vault.

Usage:  python3 obsidian/make_figures.py
Output: obsidian/attachments/fig_*.png

All figures are schematic teaching aids (not model data).  They are
kept faithful to the real test models where noted:
  - m15 2.82d has sigma < -1 at essentially all depths and |dV| up to
    ~1650 km/s per cell (~330 frequency bins at 5 km/s spacing).
  - TypeIIn 10d: slow shell ~160 km/s, fast interior ~5800 km/s,
    sigma_min ~ -98.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "attachments/"
plt.rcParams.update({
    "font.size": 11,
    "figure.dpi": 150,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ----------------------------------------------------------------------
# fig 1: velocity law + sigma, expanding vs compressing regions
# ----------------------------------------------------------------------
def fig_sigma_velocity():
    # Schematic Type-IIn-like profile: homologous ejecta (V prop. R,
    # sigma = 0) inside, steep drop through the stalled shell
    # (sigma << -1), slow gently-rising CSM outside.
    R = np.logspace(-0.6, 0.3, 2000)                     # normalized
    V_hom = 2200.0 * R                                   # km/s
    drop = 1.0 / (1.0 + np.exp(-(R - 0.85) / 0.012))     # 0 in -> 1 out
    V_csm = 60.0 + 60.0 * (R - 0.85)                     # slow CSM
    V = V_hom * (1.0 - drop) + V_csm * drop
    sigma = np.gradient(np.log(V), np.log(R)) - 1.0      # sigma_true
    comp = sigma < -1.5   # shade only the genuinely compressing zone

    fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True,
                           gridspec_kw={"hspace": 0.09})

    ax[0].plot(R, V, lw=2, color="C0")
    ax[0].axvspan(R[comp][0], R[comp][-1], color="C3", alpha=0.12)
    ax[0].set_ylabel("V(R)  [km/s]")
    ax[0].set_yscale("log")
    ax[0].set_ylim(40, 8000)
    ax[0].annotate("homologous ejecta\nV prop. R, sigma = 0",
                   xy=(0.33, 320), ha="center", fontsize=10)
    ax[0].annotate("stalled shell:\nV collapses outward\n"
                   "(reverse-shocked CDS)",
                   xy=(0.90, 700), xytext=(0.235, 2500), fontsize=10,
                   ha="left",
                   arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
    ax[0].annotate("slow dense CSM", xy=(1.35, 160), ha="center",
                   fontsize=10)
    ax[0].set_title("Non-monotonic velocity field (schematic, "
                    "m15 / TypeIIn-like)", fontsize=11)

    ax[1].plot(R, sigma, lw=2, color="C1")
    ax[1].axhline(-1.0, color="k", ls="--", lw=1)
    ax[1].axvspan(R[comp][0], R[comp][-1], color="C3", alpha=0.12)
    ax[1].set_ylabel(r"$\sigma = d\ln V/d\ln R - 1$")
    ax[1].set_xlabel("R  (normalized)")
    ax[1].set_ylim(-140, 30)
    ax[1].annotate("sigma = -1 : GAM changes sign\n"
                   "(wavelength derivative term)",
                   xy=(0.75, -1), xytext=(0.25, -60),
                   arrowprops=dict(arrowstyle="->", color="k", lw=0.8),
                   fontsize=10)
    ax[1].annotate("expansion, GAM>0:\nupwind frequency already\n"
                   "computed (live Gauss-Seidel)",
                   xy=(0.35, -120), ha="center", fontsize=9)
    ax[1].annotate("compression, GAM<0:\nupwind frequency not yet\n"
                   "computed -> lagged I_NEXT store",
                   xy=(0.88, -125), xytext=(1.02, -125), fontsize=9)
    fig.savefig(OUT + "fig_sigma_velocity.png", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# fig 2: chi_eff vs chi_hat through a compression zone
# ----------------------------------------------------------------------
def fig_chihat():
    # one ray, one frequency pair; schematic cell scan through a
    # compressing region where GAM/dlnnu is large and negative.
    x = np.linspace(0, 1, 400)
    chi = 1.0 + 0.05 * x                              # baseline opacity
    gam = -8.0 * np.exp(-((x - 0.55) / 0.12) ** 2)    # GAM/dlnnu < 0
    chi_eff = chi + gam                               # stock: can go < 0
    chi_hat = chi + np.abs(gam)                       # Knop: always >= chi

    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    ax.plot(x, chi_eff, lw=2, color="C3",
            label=r"stock: $\chi + \Gamma/d\ln\nu$")
    ax.plot(x, chi_hat, lw=2, color="C0",
            label=r"Knop: $\hat\chi = \chi + |\Gamma|/d\ln\nu$")
    ax.plot(x, chi, lw=1.2, color="gray", ls=":",
            label=r"$\chi$ (baseline)")
    ax.axhline(0, color="k", lw=1)
    ax.axvspan(0.43, 0.67, color="C3", alpha=0.10)
    ax.fill_between(x, chi_eff, 0, where=chi_eff < 0, color="C3",
                    alpha=0.35)
    ax.annotate("negative generalized opacity\n"
                "-> negative / NaN intensities",
                xy=(0.55, -3.4), xytext=(0.02, -5.6), fontsize=10,
                arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
    ax.annotate("advection into the cell acts like\n"
                "extra positive absorption; its upwind\n"
                "source |Q|*I_next stays explicit",
                xy=(0.66, 6.8), xytext=(0.60, 8.9), fontsize=9,
                arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
    ax.set_xlabel("position along ray (through compression zone)")
    ax.set_ylabel("effective opacity  (arb. units)")
    ax.set_ylim(-8.5, 12.5)
    ax.legend(loc="upper left", fontsize=9, framealpha=1.0)
    ax.set_title("Knop positive generalized opacity "
                 "(Paper I, Eq. 7)", fontsize=11)
    fig.savefig(OUT + "fig_chihat.png", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# fig 3: Sobolev resonance points, monotonic vs non-monotonic
# ----------------------------------------------------------------------
def fig_sobolev():
    s = np.linspace(-1, 1, 600)
    v_mono = s * 1.0                                  # homologous LOS
    v_nonm = s - 1.6 * np.exp(-((s - 0.25) / 0.18) ** 2)  # kink/dip
    v_res = -0.25

    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.6), sharey=True)
    for a, v, ttl in [(ax[0], v_mono, "monotonic: one Sobolev point"),
                      (ax[1], v_nonm, "non-monotonic: three Sobolev points")]:
        a.plot(s, v, lw=2, color="C0")
        a.axhline(v_res, color="k", ls="--", lw=1)
        a.set_xlabel("position s along ray")
        a.set_title(ttl, fontsize=11)
        a.set_xlim(-1, 1)
        a.set_ylim(-1.8, 1.2)
    ax[0].set_ylabel(r"$v_{\rm los}(s)$")
    ax[0].plot(v_res, v_res, "o", color="C3", ms=8, zorder=5)
    ax[0].annotate("P1", xy=(v_res + 0.05, v_res + 0.12), fontsize=10,
                   color="C3")

    # crossings of the non-monotonic curve
    idx = np.where(np.diff(np.sign(v_nonm - v_res)) != 0)[0]
    for k, i in enumerate(idx):
        ax[1].plot(s[i], v_res, "o", color="C3", ms=8, zorder=5)
        ax[1].annotate(f"P{k+1}", xy=(s[i] - 0.04, v_res + 0.10),
                       fontsize=10, color="C3")
    ax[1].annotate("comoving resonance\n(same line redshifted into\n"
                   "resonance at 3 places)",
                   xy=(-0.9, 0.45), fontsize=10)
    fig.savefig(OUT + "fig_sobolev_points.png", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# fig 4: three-phase convergence schematic (Lambda -> frz -> cpl JFNK)
# ----------------------------------------------------------------------
def fig_phases():
    it = np.arange(0, 90)
    # stylized demanded-correction decay with the observed signatures:
    # slow Lambda grind, then Newton plunge.
    demand = 1e8 * np.exp(-it / 16.0)
    demand += 1e4 * np.exp(-((it - 30) / 6.0) ** 2)   # mid-run bump
    demand[it > 62] = 8e1 * np.exp(-(it[it > 62] - 62) / 3.0)

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    ax.semilogy(it, demand, lw=2, color="C0")
    ax.axhline(1e3, color="C2", ls="--", lw=1.2)
    ax.axhline(1e2, color="C4", ls="--", lw=1.2)
    ax.text(1, 1.6e3, "JFNK_FRZ_GATE = 1e3", color="C2", fontsize=10)
    ax.text(1, 1.4e2, "JFNK_COR_GATE = 1e2", color="C4", fontsize=10)

    ax.axvspan(0, 55, color="0.5", alpha=0.10)
    ax.axvspan(55, 70, color="C2", alpha=0.12)
    ax.axvspan(70, 90, color="C4", alpha=0.12)
    ax.text(27, 4e4, "Phase 0: Lambda + Knop sweeps\n"
            "field matures, demand decays\n(GS_FREEZE until the store is\n"
            "a fixed point; capped steps)",
            ha="center", fontsize=9)
    ax.text(62.5, 5e5, "Phase 1:\nfrozen-J\nJFNK\n(cheap)",
            ha="center", va="center", fontsize=9, color="C2")
    ax.text(80, 5e5, "Phase 2:\ncoupled\nJFNK\n(endgame)",
            ha="center", va="center", fontsize=9, color="C4")
    ax.set_xlabel("iteration")
    ax.set_ylabel("demanded correction  MAXCH  [%]")
    ax.set_title("Three-phase convergence (schematic)", fontsize=11)
    ax.set_xlim(0, 90)
    ax.set_ylim(0.5, 3e8)
    fig.savefig(OUT + "fig_three_phases.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_sigma_velocity()
    fig_chihat()
    fig_sobolev()
    fig_phases()
    print("wrote fig_sigma_velocity.png fig_chihat.png",
          "fig_sobolev_points.png fig_three_phases.png ->", OUT)
