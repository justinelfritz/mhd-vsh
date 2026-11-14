#!/usr/bin/env python3
"""Plots f_H(r), B_r(r), and Rmag(r) for the old ad hoc seed IC vs the
new l=1 force-free dipole IC, on the same production grid, as written by
app/mhdvsh_ic_rmag_compare.f90. Tests the hypothesis that the new IC's
much smaller Hall-CFL dt comes from its field being concentrated near
the surface -- exactly where f_H/Rmag are most extreme -- unlike the old
seed, which is engineered to vanish at both boundaries.

Usage: plot_ic_rmag_compare.py <data_file> <output_png>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png>")
    data_path, out_path = sys.argv[1], sys.argv[2]

    data = np.loadtxt(data_path, comments="#")
    r, f_hall = data[:, 0], data[:, 2]
    old_br, old_rmag = data[:, 5], data[:, 6]
    new_br, new_rmag = data[:, 8], data[:, 9]

    fig, axes = plt.subplots(3, 1, figsize=(7, 8.5), dpi=150, sharex=True)
    fig.patch.set_facecolor(SURFACE)

    ax = axes[0]
    ax.set_facecolor(SURFACE)
    ax.plot(r, f_hall, color=PALETTE[3], linewidth=2.2, solid_capstyle="round")
    ax.set_yscale("log")
    ax.set_ylabel(r"$f_H(r)$  [code units]", color=TEXT_SECONDARY, fontsize=10)

    ax = axes[1]
    ax.set_facecolor(SURFACE)
    ax.plot(r, old_br, color=PALETTE[0], linewidth=2.2, solid_capstyle="round", label="old ad hoc seed")
    ax.plot(r, new_br, color=PALETTE[1], linewidth=2.2, solid_capstyle="round", label="new force-free dipole")
    ax.set_ylabel(r"$B_r$ at the pole  $[10^{12}\mathrm{G}]$", color=TEXT_SECONDARY, fontsize=10)
    ax.legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=9, loc="upper left")

    ax = axes[2]
    ax.set_facecolor(SURFACE)
    ax.plot(r, np.abs(old_rmag), color=PALETTE[0], linewidth=2.2, solid_capstyle="round", label="old ad hoc seed")
    ax.plot(r, np.abs(new_rmag), color=PALETTE[1], linewidth=2.2, solid_capstyle="round", label="new force-free dipole")
    ax.set_yscale("log")
    ax.set_ylabel(r"$R_{\rm mag}=f_H B_r/\eta$", color=TEXT_SECONDARY, fontsize=10)
    ax.legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=9, loc="upper left")

    for ax in axes:
        ax.grid(True, color=GRID, linewidth=0.8, which="both")
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(GRID)
        ax.tick_params(colors=TEXT_SECONDARY)

    axes[-1].set_xlabel("r [km]", color=TEXT_SECONDARY)
    fig.suptitle("Old seed vs new force-free dipole: where is the field strong?",
                 color=TEXT_PRIMARY, fontsize=12, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
