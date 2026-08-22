#!/usr/bin/env python3
"""Plots the time-dependent growth of every Phi_l0 or Psi_l0 mode from a
Hall stability experiment run (app/mhdvsh_hall_stability_experiment.f90).

Reads the per-mode amplitude series it logs (step, t, Phi_l1..Phi_lLMAX,
Psi_l1..Psi_lLMAX, max_imag_fraction) and plots Re(<field>_l0) at the
domain midpoint against t, one line per l, for whichever field is named
on the command line. l is an ORDERED index (mode degree), not a
category, so it's colored with a single-hue light->dark ordinal ramp
(dataviz skill's documented blue ramp, references/palette.md) and a
colorbar rather than a 30-entry legend -- consistent with how the skill
treats position-in-a-sequence encodings, not identity. Both fields share
one style (same ramp, same layout) so the two plots read as a matched
pair, not as unrelated charts.

Usage: plot_hall_mode_amplitudes.py <data_file> <output_png> <Phi|Psi>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np

# Documented ordinal blue ramp, dataviz skill references/palette.md,
# steps 250->700 (light mode: lightest step no lighter than step 250).
BLUE_RAMP = [
    "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
    "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]
BLUE_CMAP = LinearSegmentedColormap.from_list("ordinal_blue", BLUE_RAMP)

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
    if len(sys.argv) != 4 or sys.argv[3] not in ("Phi", "Psi"):
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png> <Phi|Psi>")
    data_path, out_path, field = sys.argv[1], sys.argv[2], sys.argv[3]

    names, data = read_columns(data_path)
    t = data[:, 1]
    # Columns are step, t, Phi_l1..Phi_lLMAX, Psi_l1..Psi_lLMAX,
    # max_imag_fraction -- every "<field>_l*" column is one mode of the
    # requested field; the last column is a realness sanity check shared
    # across both fields (whichever drifted furthest from real), not a
    # mode to plot.
    prefix = f"{field}_l"
    mode_cols = [i for i, n in enumerate(names) if n.startswith(prefix)]
    l_values = [int(names[i].replace(prefix, "")) for i in mode_cols]
    max_imag_fraction = data[:, -1]

    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    norm = Normalize(vmin=min(l_values), vmax=max(l_values))
    for i, l in zip(mode_cols, l_values):
        ax.plot(t, data[:, i], color=BLUE_CMAP(norm(l)), linewidth=1.5,
                 solid_capstyle="round")

    # Modes span ~69 orders of magnitude (l=1 ~1, l=30 ~1e-69) -- most of
    # that range is below double-precision noise (~1e-16 relative to the
    # largest mode), not real cascade content. symlog (not a magnitude-
    # only log) keeps the sign information a linear axis would hide
    # entirely, while the linear threshold at the actual noise floor
    # stops the plot from pretending to resolve meaningless digits.
    NOISE_FLOOR = 1.0E-15
    ax.set_yscale("symlog", linthresh=NOISE_FLOOR)
    ax.axhspan(-NOISE_FLOOR, NOISE_FLOOR, color=GRID, alpha=0.5, zorder=0)

    symbol = r"\Phi" if field == "Phi" else r"\Psi"
    ax.set_xlabel("t [yr]", color=TEXT_SECONDARY)
    ax.set_ylabel(rf"Re(${symbol}_{{l,0}}$) at midpoint (symlog)", color=TEXT_SECONDARY)
    ax.set_title(f"Hall stability experiment: {field} mode growth vs time",
                 color=TEXT_PRIMARY, fontsize=12, fontweight="bold", loc="left")

    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY)

    sm = ScalarMappable(norm=norm, cmap=BLUE_CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("l (mode degree)", color=TEXT_SECONDARY)
    cbar.ax.tick_params(colors=TEXT_SECONDARY)
    cbar.outline.set_edgecolor(GRID)

    worst_imag = float(np.max(max_imag_fraction)) if len(max_imag_fraction) else 0.0
    ax.text(0.02, 0.02,
            f"max |Im|/|.| across all Phi/Psi modes/steps: {worst_imag:.2e}",
            transform=ax.transAxes, fontsize=7.5, color=TEXT_SECONDARY,
            va="bottom", ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")
    print(f"max |Im|/|.| across all Phi/Psi modes/steps: {worst_imag:.6e}")


if __name__ == "__main__":
    main()
