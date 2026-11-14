#!/usr/bin/env python3
"""Plots the radial profile of a force-free Bessel-Riccati IC (Phi, Psi,
and B_r at the magnetic pole) both inside the star (the real eigenmode
solution) and outside it (the analytic vacuum continuation), as written
by app/mhdvsh_force_free_profile.f90.

B_r is reported at the magnetic pole (theta=0), where the l=1 angular
factor is exactly 1 -- see that driver's own header for why (no (Phi,Psi)
-> physical B converter exists yet in this codebase, ROADMAP.md
programming task 2). Psi is plotted directly as the toroidal potential,
not converted to a B_theta/B_phi amplitude, for the same reason.

Usage: plot_force_free_profile.py <data_file> <output_png>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Same validated categorical palette as scripts/plot_energy_budget.py
# (dataviz skill's reference/palette.md 8-slot default), for consistency.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"
SURFACE_LINE_COLOR = "#8a8880"


def read_metadata(path):
    """Pulls R_min/R_max/mu out of the Fortran driver's own '#' header
    lines (written as 'key = value ...'), rather than hardcoding them."""
    meta = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            s = line.lstrip("#").strip()
            for key, tag in [("R_min", "R_min"), ("R_max", "R_max"), ("mu", "mu    =")]:
                if s.startswith(tag):
                    meta[key] = float(s.split("=")[1].split()[0])
    return meta


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png>")
    data_path, out_path = sys.argv[1], sys.argv[2]

    meta = read_metadata(data_path)
    r_min, r_max = meta["R_min"], meta["R_max"]

    rows = np.genfromtxt(data_path, comments="#", dtype=None, encoding="utf-8")
    region = np.array([row[0] for row in rows])
    r = np.array([row[1] for row in rows], dtype=float)
    phi = np.array([row[2] for row in rows], dtype=float)
    psi = np.array([row[3] for row in rows], dtype=float)
    b_r = np.array([row[4] for row in rows], dtype=float)

    is_int = region == "interior"
    is_ext = region == "exterior"

    fig, axes = plt.subplots(3, 1, figsize=(7, 8.5), dpi=150, sharex=True)
    fig.patch.set_facecolor(SURFACE)

    panels = [
        (axes[0], phi, r"$\Phi_{l=1,m=0}(r)$", PALETTE[0]),
        (axes[1], psi, r"$\Psi_{l=1,m=0}(r)$", PALETTE[1]),
        (axes[2], b_r, r"$B_r$ at the pole  $[10^{12}\ \mathrm{G}]$", PALETTE[2]),
    ]
    for ax, y, ylabel, color in panels:
        ax.set_facecolor(SURFACE)
        ax.plot(r[is_int], y[is_int], color=color, linewidth=2.2,
                 solid_capstyle="round", label="interior (crust)")
        ax.plot(r[is_ext], y[is_ext], color=color, linewidth=2.2, linestyle="--",
                 dash_capstyle="round", label="exterior (vacuum)")
        ax.axvline(r_min, color=SURFACE_LINE_COLOR, linewidth=1.0, linestyle=":", zorder=0)
        ax.axvline(r_max, color=SURFACE_LINE_COLOR, linewidth=1.2, zorder=0)
        ax.set_ylabel(ylabel, color=TEXT_SECONDARY, fontsize=10)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.axhline(0.0, color=GRID, linewidth=1.0, zorder=0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(GRID)
        ax.tick_params(colors=TEXT_SECONDARY)

    ylo, yhi = axes[-1].get_ylim()
    axes[-1].text(r_max, ylo + 0.06*(yhi-ylo), " surface (R_max)", color=TEXT_SECONDARY,
                  fontsize=8, va="bottom", ha="left")
    axes[-1].text(r_min, ylo + 0.06*(yhi-ylo), "core/crust (R_min) ", color=TEXT_SECONDARY,
                  fontsize=8, va="bottom", ha="right")
    axes[0].legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=9, loc="lower right")

    axes[-1].set_xlabel("r [km]", color=TEXT_SECONDARY)
    fig.suptitle(f"Force-free Bessel-Riccati IC, l=1 ($\\mu$={meta['mu']:.4f} km$^{{-1}}$)",
                 color=TEXT_PRIMARY, fontsize=12, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
