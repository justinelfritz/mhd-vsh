#!/usr/bin/env python3
"""Plots the time-dependent growth of magnetic energy stored in each
poloidal or toroidal mode, for drivers that log the FIELD_DIAGNOSTICS::
POLOIDAL_MAGNETIC_ENERGY_BY_L/TOROIDAL_MAGNETIC_ENERGY_BY_L column
convention in physical erg (UNITS::ENERGY_UNIT_ERG) -- app/mhdvsh_hall.f90,
app/mhdvsh_hall_crust_profile.f90.

Reads the per-degree energy series (step, t, Epol_l1..Epol_lLMAX,
Etor_l1..Etor_lLMAX -- the mhd-vsh-relations.pdf E_B,pol(t)/E_B,tor(t)
volume-integrated energy formulas broken out by degree l instead of
summed over it) and plots one line per l against t, for whichever field
is named on the command line. Same styling as plot_hall_mode_amplitudes.py
(ordinal blue ramp by l, colorbar not a 30-entry legend). Plain log, not
symlog: energy is non-negative by construction, no sign to preserve.

Log-x, and y in units of 10^40 erg -- per user direction (2026-08-29).
Values plotted in erg directly (no unit conversion applied elsewhere in
this script's own history) span the full range from ~1e40 down through
the double-precision noise floor several times over, so a FIXED-unit
y-axis (rather than the prior "deliberately unit-agnostic" design, which
served a now-unused code-units harness -- app/mhdvsh_hall_stability_
experiment.f90 was never wired to this script via any build target) is
the right choice for the erg-only driver family this script now targets.
log-x spreads the early, fast-evolving transient (the first ~10-100yr,
where these per-l curves show real dynamical structure -- dips/wobbles
matching ROADMAP.md's still-open "first ~10yr energy-balance violation"
item) out against the full multi-century run the same way
plot_energy_budget.py's own log-x does.

Usage: plot_hall_energy_by_l.py <data_file> <output_png> <pol|tor>
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

ENERGY_UNIT = 1.0E40  # erg -- fixed unit for this script's y-axis


def read_columns(path):
    with open(path) as f:
        header_line = next(line for line in f if line.startswith("#"))
    names = header_line.lstrip("#").split()
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return names, data


def main():
    if len(sys.argv) != 4 or sys.argv[3] not in ("pol", "tor"):
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png> <pol|tor>")
    data_path, out_path, field = sys.argv[1], sys.argv[2], sys.argv[3]

    names, data = read_columns(data_path)
    t = data[:, 1]
    prefix = f"E{field}_l"
    mode_cols = [i for i, n in enumerate(names) if n.startswith(prefix)]
    l_values = [int(names[i].replace(prefix, "")) for i in mode_cols]

    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    norm = Normalize(vmin=min(l_values), vmax=max(l_values))
    for i, l in zip(mode_cols, l_values):
        ax.plot(t, data[:, i] / ENERGY_UNIT, color=BLUE_CMAP(norm(l)), linewidth=1.5,
                 solid_capstyle="round")

    # Energy is non-negative by construction (it's |.|**2-built), so a
    # plain log axis is the right tool -- no sign to lose, unlike the
    # Phi/Psi amplitude plots which need symlog. Values that are exactly
    # zero (a mode never populated, e.g. every even l in the poloidal
    # cascade) are simply omitted by the log scale rather than plotted
    # at -inf; matplotlib does this automatically for zero/negative
    # values on a log axis.
    ax.set_yscale("log")
    # log-x: the seed IC's fast initial relaxation (first ~10-100yr) and
    # the slow multi-century decay afterward both need visible space --
    # same rationale as plot_energy_budget.py's own log-x. t=0 (the IC's
    # own row) is dropped automatically, same as log(0) would be.
    ax.set_xscale("log")
    # Floor set RELATIVE to this run's own peak (not a fixed absolute
    # number -- unaffected by the ENERGY_UNIT division above, since
    # dividing every series by the same constant leaves peak/value
    # ratios unchanged). 1e-30 relative crops the double-precision noise
    # floor (energy is amplitude-squared, so a ~1e-16 relative amplitude
    # floor becomes ~1e-32 in energy) while keeping real cascade content
    # visible.
    peak = np.nanmax(data[:, mode_cols]) if mode_cols else 1.0
    if peak > 0:
        ax.set_ylim(bottom=peak * 1.0E-30 / ENERGY_UNIT)

    # Regime name inferred from the data path's own artifacts/<regime>/...
    # convention, same as plot_energy_budget.py -- this script is shared
    # across driver families, not specific to one experiment.
    path_parts = Path(data_path).parts
    if "artifacts" in path_parts:
        title = f"{path_parts[path_parts.index('artifacts') + 1].capitalize()} regime: " \
                f"{'poloidal' if field == 'pol' else 'toroidal'} energy by mode vs time"
    else:
        title = f"{'Poloidal' if field == 'pol' else 'Toroidal'} energy by mode vs time: {Path(data_path).parent.name}"

    ax.set_xlabel("t [yr]", color=TEXT_SECONDARY)
    ax.set_ylabel(rf"$E_{{\mathrm{{{field}}}}}(l)$  $\left[10^{{40}}\ \mathrm{{erg}}\right]$", color=TEXT_SECONDARY)
    ax.set_title(title, color=TEXT_PRIMARY, fontsize=12, fontweight="bold", loc="left")

    ax.grid(True, color=GRID, linewidth=0.8, which="both")
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

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
