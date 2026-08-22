#!/usr/bin/env python3
"""Plots a comprehensive view of the crust EOS/conductivity profile
produced by app/mhdvsh_tov.f90: number density, mass density, eta,
f_Hall, and the derived magnetic Reynolds number Rmag=f_H*Bmax/eta,
all against radius r, on one figure with five independently scaled
y-axes sharing the same x-axis.

Bmax (the max magnetic field magnitude at t=0, code units of 10**12 G)
is a fixed constant, not read from the data file -- this profile has no
attached field simulation, so it must be supplied. Default 1.912
matches the max Br at t=0 for the standard single-mode Phi(l=1,m=0)
seed IC used by app/mhdvsh_hall.f90/mhdvsh_hall_adaptive.f90; override
with a 3rd CLI argument for a different IC/simulation.

Reads the whitespace-delimited data file mhdvsh_tov writes (header:
"# r_km  rhocgs  n_e_cm-3  eta_km2_per_yr  f_hall_km2_per_1e12G_per_yr")
and writes a PNG.

Usage: plot_eos_profile.py <data_file> <output_png> [Bmax]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# High-contrast categorical palette -- deliberately NOT the shared
# plot_energy_budget.py PALETTE (its 4th slot is a yellow/amber too easy
# to lose against the light SURFACE background, flagged by the user
# 2026-08-21 specifically for f_H). No yellow anywhere in this set.
BLUE   = "#1c6dd0"
RED    = "#d6393a"
GREEN  = "#178a5a"
PURPLE = "#7b4fb0"
MAGENTA = "#c2185b"

DEFAULT_BMAX = 1.912  # code units, 10**12 G -- see module docstring

TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"


def read_columns(path):
    with open(path) as f:
        header_line = next(line for line in f if line.startswith("#"))
    names = header_line.lstrip("#").split()
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return names, data


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png> [Bmax]")
    data_path, out_path = sys.argv[1], sys.argv[2]
    bmax = float(sys.argv[3]) if len(sys.argv) == 4 else DEFAULT_BMAX

    _, data = read_columns(data_path)
    r, rhocgs, n_e, eta, f_hall = (data[:, i] for i in range(5))
    rmag = f_hall * bmax / eta

    # (label, values, color, outward spine offset in points)
    series = [
        (r"$n_e$ [cm$^{-3}$]", n_e, BLUE, 0),
        (r"$\rho$ [g cm$^{-3}$]", rhocgs, RED, 52),
        (r"$\eta$ [km$^2$ yr$^{-1}$]", eta, GREEN, 104),
        (r"$f_H$ [km$^2$ (10$^{12}$G)$^{-1}$ yr$^{-1}$]", f_hall, PURPLE, 156),
        (rf"$R_{{mag}}=f_H B_{{max}}/\eta$, $B_{{max}}$={bmax:g}", rmag, MAGENTA, 208),
    ]

    fig, host = plt.subplots(figsize=(10.5, 5.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    host.set_facecolor(SURFACE)

    axes = [host]
    for _ in range(len(series) - 1):
        axes.append(host.twinx())

    lines = []
    for ax, (label, values, color, offset) in zip(axes, series):
        (line,) = ax.semilogy(r, values, color=color, linewidth=2.2,
                               solid_capstyle="round", label=label)
        lines.append(line)
        ax.set_ylabel(label, color=color, fontsize=10)
        ax.tick_params(axis="y", colors=color)
        ax.spines["right"].set_visible(ax is not host)
        if ax is not host:
            ax.spines["right"].set_position(("outward", offset))
            ax.spines["right"].set_color(color)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)
        else:
            ax.spines["left"].set_color(BLUE)
            ax.spines["top"].set_visible(False)

    host.set_xlabel("r [km]", color=TEXT_SECONDARY)
    host.set_title("Neutron star crust EOS profile", color=TEXT_PRIMARY,
                    fontsize=13, fontweight="bold", loc="left")
    host.grid(True, color=GRID, linewidth=0.8)
    host.set_axisbelow(True)
    host.tick_params(axis="x", colors=TEXT_SECONDARY)

    legend = host.legend(lines, [l.get_label() for l in lines],
                          frameon=False, labelcolor=TEXT_SECONDARY,
                          loc="upper center", bbox_to_anchor=(0.5, -0.12),
                          ncol=3)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE, bbox_extra_artists=(legend,),
                bbox_inches="tight")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
