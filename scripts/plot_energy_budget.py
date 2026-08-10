#!/usr/bin/env python3
"""Plots a regime's time-resolved energy-budget series (total magnetic
energy today; Joule dissipation, Poynting flux, etc. once that physics
is derived and logged as additional columns -- see
app/mhdvsh_diffusion.f90 and src/core/field_diagnostics.f90).

Reads a whitespace-delimited data file with a header comment line
naming the columns (as written by the Fortran driver), plots every
column after the first two (step, t) against t, and writes a PNG.
Written generically against the column set so that adding a term in
Fortran (a new column) makes it appear here with no script changes.

Usage: plot_energy_budget.py <data_file> <output_png>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Validated categorical palette (light mode), fixed assignment order --
# see the dataviz skill's references/palette.md. Column 3 (first
# plotted series) always gets slot 1, column 4 slot 2, and so on, so
# colors stay stable as terms are added rather than being reassigned.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]

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
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png>")
    data_path, out_path = sys.argv[1], sys.argv[2]

    names, data = read_columns(data_path)
    t = data[:, 1]
    series_names = names[2:]
    series = data[:, 2:]

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    for i, name in enumerate(series_names):
        ax.plot(t, series[:, i], color=PALETTE[i % len(PALETTE)],
                 linewidth=2, solid_capstyle="round",
                 label=name.replace("_", " "))

    ax.set_xlabel("t", color=TEXT_SECONDARY)
    ax.set_ylabel("energy", color=TEXT_SECONDARY)
    ax.set_title("Diffusion regime: energy budget vs time", color=TEXT_PRIMARY,
                 fontsize=12, fontweight="bold", loc="left")

    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY)

    if len(series_names) > 1:
        legend = ax.legend(frameon=False, labelcolor=TEXT_SECONDARY)
    else:
        # A single series needs no legend box -- the title/axis label
        # already names it (dataviz skill: legend only for >=2 series).
        pass

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
