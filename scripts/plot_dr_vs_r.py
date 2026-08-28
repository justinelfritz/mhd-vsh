#!/usr/bin/env python3
"""Plots the TOV profile's own radial grid spacing dr(r) = r[i+1]-r[i]
against r, for a data file in mhdvsh_tov's own format (same file
plot_eos_profile.py reads). Diagnostic for how densely the log(P)-
uniform TOV resampling clusters points near the surface, independent
of (and much denser than) the actual simulation's own uniform-r grid.

Optionally marks a "knee" radius (e.g. found via a log-log slope scan
of dr vs rho, since the r column's own print precision makes dr read
as exactly 0 -- and the slope diagnostic meaningless -- for many rows
right at the surface) with a vertical line and annotation.

Usage: plot_dr_vs_r.py <data_file> <output_png> [knee_r_km]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE = "#1c6dd0"
RED = "#d6393a"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png> [knee_r_km]")
    data_path, out_path = sys.argv[1], sys.argv[2]
    knee_r = float(sys.argv[3]) if len(sys.argv) == 4 else None

    r = np.loadtxt(data_path, comments="#", usecols=0)
    dr = np.diff(r)
    r_mid = 0.5 * (r[:-1] + r[1:])

    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.semilogy(r_mid, dr, color=BLUE, linewidth=2.0, marker="o",
                markersize=3, solid_capstyle="round", zorder=2)

    if knee_r is not None:
        ax.axvline(knee_r, color=RED, linewidth=1.6, linestyle="--", zorder=1)
        ymin, ymax = ax.get_ylim()
        y_label = ymax / 2.2
        r_max = r[-1]
        ax.annotate(
            f"knee: r={knee_r:.4f} km\n({(r_max-knee_r)*1000:.1f} m below surface)",
            xy=(knee_r, y_label), xytext=(8, 0), textcoords="offset points",
            color=RED, fontsize=9, va="center", ha="left",
        )

    ax.set_xlabel("r [km]", color=TEXT_SECONDARY)
    ax.set_ylabel(r"$\Delta r$ [km]  (TOV profile grid spacing)", color=TEXT_SECONDARY)
    ax.set_title("TOV profile radial grid spacing vs r", color=TEXT_PRIMARY,
                 fontsize=13, fontweight="bold", loc="left")
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=TEXT_SECONDARY)
    for spine in ax.spines.values():
        spine.set_color(GRID)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
