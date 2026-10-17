#!/usr/bin/env python3
"""Energy per mode (E_pol_l, E_tor_l vs l) at the LAST available
timestep of a production_run_energy_by_l.dat file -- a spectrum
snapshot, complementary to plot_energy_by_l_heatmap.py's time-resolved
view: this answers "how much energy is actually in each l, right now"
directly, rather than "when did each l become active".

Usage: plot_energy_spectrum_snapshot.py <data_file> <output_png>
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
        header_line = next(line for line in f if line.startswith("#"))
    names = header_line.lstrip("#").split()
    data = np.loadtxt(data_path, comments="#")

    t = data[-1, names.index("t")]
    pol_cols = [i for i, n in enumerate(names) if n.startswith("Epol_l")]
    tor_cols = [i for i, n in enumerate(names) if n.startswith("Etor_l")]
    l_pol = np.array([int(names[i][len("Epol_l"):]) for i in pol_cols])
    l_tor = np.array([int(names[i][len("Etor_l"):]) for i in tor_cols])
    order_p = np.argsort(l_pol); order_t = np.argsort(l_tor)
    l_pol, e_pol = l_pol[order_p], data[-1, [pol_cols[i] for i in order_p]]
    l_tor, e_tor = l_tor[order_t], data[-1, [tor_cols[i] for i in order_t]]

    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.plot(l_pol, e_pol, "o-", color=PALETTE[0], linewidth=1.8, markersize=4, label="Poloidal")
    ax.plot(l_tor, e_tor, "o-", color=PALETTE[1], linewidth=1.8, markersize=4, label="Toroidal")
    ax.set_yscale("log")
    ax.set_xlabel("l (mode degree)", color=TEXT_SECONDARY)
    ax.set_ylabel(r"$E_l$  [erg]", color=TEXT_SECONDARY)
    ax.grid(True, color=GRID, linewidth=0.8, which="both")
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY)
    ax.legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=10)

    run_name = Path(data_path).parent.name
    fig.suptitle(f"Energy per mode at t={t:.4f} yr: {run_name}", color=TEXT_PRIMARY,
                 fontsize=12, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
