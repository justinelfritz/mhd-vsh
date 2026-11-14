#!/usr/bin/env python3
"""Plots B_r(r) at the magnetic pole for the l=1 force-free dipole
production run, at the two full-field snapshots actually available:
t=0 (the IC, reconstructed exactly from the namelist) and the latest
checkpoint. As written by app/mhdvsh_br_snapshot_compare.f90.

Usage: plot_br_snapshots.py <data_file> <output_png>
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

    with open(data_path) as f:
        header = f.readline()
    t_yr = float(header.split("t_yr=")[1])

    data = np.loadtxt(data_path, comments="#")
    r, br_t0, br_ck = data[:, 0], data[:, 1], data[:, 2]

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.plot(r, br_t0, color=PALETTE[0], linewidth=2.2, solid_capstyle="round", label="t=0 (IC)")
    ax.plot(r, br_ck, color=PALETTE[1], linewidth=2.2, solid_capstyle="round", label=f"t={t_yr:.3f} yr (checkpoint)")
    ax.axhline(0.0, color=GRID, linewidth=1.0, zorder=0)

    ax.set_xlabel("r [km]", color=TEXT_SECONDARY)
    ax.set_ylabel(r"$B_r$ at the pole  $[10^{12}\mathrm{G}]$", color=TEXT_SECONDARY)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY)
    ax.legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=10)

    fig.suptitle("l=1 force-free dipole: B_r(r) at the pole, t=0 vs latest checkpoint",
                 color=TEXT_PRIMARY, fontsize=12, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
