#!/usr/bin/env python3
"""Plots Rmag(r) two ways on the SAME fine (387-point, surface-clustered)
radial grid already used for results/eos_comparison/new_eos_profile.dat:
the original flat-B_max=1.912 reference convention vs the real,
radially-varying B_r(r) of the new l=1 force-free dipole IC, as written
by app/mhdvsh_ic_rmag_eos_compare.f90.

Usage: plot_ic_rmag_eos_compare.py <data_file> <output_png>
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
    r = data[:, 0]
    br_new = data[:, 3]
    rmag_const = data[:, 4]
    rmag_real = data[:, 5]

    fig, axes = plt.subplots(2, 1, figsize=(7, 6.5), dpi=150, sharex=True)
    fig.patch.set_facecolor(SURFACE)

    ax = axes[0]
    ax.set_facecolor(SURFACE)
    ax.plot(r, br_new, color=PALETTE[1], linewidth=2.2, solid_capstyle="round")
    ax.set_ylabel(r"$B_r$, new dipole  $[10^{12}\mathrm{G}]$", color=TEXT_SECONDARY, fontsize=10)

    ax = axes[1]
    ax.set_facecolor(SURFACE)
    ax.plot(r, np.abs(rmag_const), color=PALETTE[0], linewidth=2.2, solid_capstyle="round",
            label=r"const $B_{max}=1.912$ (eos_comparison's own reference)")
    ax.plot(r, np.abs(rmag_real), color=PALETTE[1], linewidth=2.2, solid_capstyle="round",
            label="real new force-free dipole $B_r(r)$")
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
    fig.suptitle("Rmag on the true surface-resolved grid: constant field vs real dipole",
                 color=TEXT_PRIMARY, fontsize=11, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
