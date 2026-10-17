#!/usr/bin/env python3
"""Spectrogram-style view of the per-l energy cascade: l (mode degree)
on the y-axis, t on the x-axis, log10(E_l) as color -- replaces plotting
all 30 modes as 30 overlapping lines spanning ~30 orders of magnitude
(unreadable) with a single compact 2D view, plus a derived "cascade
front" line (the highest l whose energy is within FRONT_REL_THRESHOLD
of the dominant l=1 mode at that instant) as the one number that most
directly answers "how far has the cascade reached, and how fast is it
advancing toward LMAX".

Usage: plot_energy_by_l_heatmap.py <data_file> <output_png> <pol|tor>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"
FRONT_COLOR = "#eb6834"

ENERGY_UNIT = 1.0E40  # erg, matches plot_hall_energy_by_l.py's own convention
FRONT_REL_THRESHOLD = 1.0E-6  # "active" := within 1e-6 of the dominant l's own energy


def main():
    if len(sys.argv) != 4 or sys.argv[3] not in ("pol", "tor"):
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png> <pol|tor>")
    data_path, out_path, field = sys.argv[1], sys.argv[2], sys.argv[3]

    with open(data_path) as f:
        header_line = next(line for line in f if line.startswith("#"))
    names = header_line.lstrip("#").split()
    data = np.loadtxt(data_path, comments="#")

    prefix = "Epol_l" if field == "pol" else "Etor_l"
    l_cols = [i for i, n in enumerate(names) if n.startswith(prefix)]
    lvals = np.array([int(names[i][len(prefix):]) for i in l_cols])
    order = np.argsort(lvals)
    lvals = lvals[order]
    l_cols = [l_cols[i] for i in order]

    t = data[:, names.index("t")]
    E = data[:, l_cols] / ENERGY_UNIT   # (n_steps, n_l), units of 1e40 erg

    # Cascade front: highest l with E_l within FRONT_REL_THRESHOLD of the
    # dominant (l=1) mode's own energy at that instant -- a single,
    # scale-invariant trace instead of guessing an absolute erg cutoff.
    e_dominant = E[:, 0]
    front = np.full(len(t), np.nan)
    for i in range(len(t)):
        active = np.where(E[i, :] >= FRONT_REL_THRESHOLD * e_dominant[i])[0]
        if len(active) > 0:
            front[i] = lvals[active.max()]

    logE = np.log10(np.clip(E, 1e-300, None))
    # Floor the color scale at the noise level just below the smallest
    # genuinely-resolved value (not at the data's own arbitrary minimum,
    # which is often double-precision garbage many more decades down).
    vmax = np.nanmax(logE)
    vmin = max(np.nanmin(logE[np.isfinite(logE)]), vmax - 30)

    fig, ax = plt.subplots(figsize=(7.5, 5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    mesh = ax.pcolormesh(t, lvals, logE.T, shading="nearest", cmap="viridis",
                          vmin=vmin, vmax=vmax)
    ax.plot(t, front, color=FRONT_COLOR, linewidth=2.2, label=f"cascade front ($E_l \\geq 10^{{-6}} E_{{l=1}}$)")

    ax.set_xscale("log")
    ax.set_xlabel("t [yr]", color=TEXT_SECONDARY)
    ax.set_ylabel("l (mode degree)", color=TEXT_SECONDARY)
    ax.tick_params(colors=TEXT_SECONDARY)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=9, loc="upper left")

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10} E_l$  $[10^{40}\,\mathrm{erg}]$", color=TEXT_SECONDARY)
    cbar.ax.tick_params(colors=TEXT_SECONDARY)

    field_label = "Poloidal" if field == "pol" else "Toroidal"
    run_name = Path(data_path).parent.name
    fig.suptitle(f"{field_label} energy cascade: {run_name}", color=TEXT_PRIMARY,
                 fontsize=12, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
