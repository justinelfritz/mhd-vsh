#!/usr/bin/env python3
"""Plots a regime's energy-budget rate terms vs time: dE_pol/dt, dE_tor/dt,
Joule dissipation rate, the (diffusive + Hall) Poynting flux rates, and
their total balance E_tot (per user direction, 2026-08-29) -- the signed
sum of all five, which should sit near zero if the run's energy budget
closes; growth away from zero flags a real (or numerical) imbalance.

Reads a whitespace-delimited data file with a header comment line naming
the columns (as written by the Fortran driver). Only rate-type quantities
are plotted -- raw energies (erg) and step-size bookkeeping (dt_yr,
tc_raw_yr) are deliberately excluded, per user direction (2026-08-25):
"dt and tc do not belong in this plot. The only quantities we want are
the time derivate of the B-field components (pol and tor), and joule
rate, and the poynting fluxes."

Plotted in units of 10^40 erg/Myr (1 Myr = 10^6 yr), per user direction
(2026-08-28) -- the Fortran driver logs erg/s, which put the y-axis on a
symlog scale spanning ~60 orders of magnitude peak-to-trough. Converting
units alone doesn't shrink that span (a linear rescale shifts every
value by the same constant, so log(x*c) - log(y*c) = log(x) - log(y));
what actually matters is that 10^40 erg/Myr is close to this problem's
own natural scale (total crustal magnetic energy ~10^38-10^39 erg,
evolving over ~Myr), so the bulk of the run's own quasi-steady-state
values land within a normal O(1)-O(100) linear range. The very first
~10s of years (the initial relaxation transient, before the profile
settles) still swings far outside that band -- left on-scale rather
than clipped, since silently cropping real data would misrepresent the
run; see the module docstring in the driver script's own output/report
for that caveat.

If the data file already has EDOT_poloidal_erg_per_s/EDOT_toroidal_erg_per_s
columns (app/mhdvsh_hall_crust_profile.f90, 2026-08-25 onward), those are
used directly. Older data files that only logged raw E_poloidal_erg/
E_toroidal_erg have dE/dt derived here via the same backward-difference
the Fortran driver itself uses internally for its residual term, so
existing runs don't need to be re-executed just to get this plot.

Usage: plot_energy_budget.py <data_file> <output_png>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# src/core/units.f90's TIME_UNIT_S -- needed because dt_yr is in years but
# the EDOT_* columns (when derived here rather than read directly from
# the Fortran driver) must be in erg/s to match joule_dissipation_rate_
# erg_per_s etc. Dividing raw-energy deltas by dt_yr directly (years, not
# seconds) was an early bug in this fallback: it under-converted by
# ~3.16e7x, making dE/dt appear ~7-8 orders of magnitude larger than the
# Joule/Poynting rates it should balance against.
SECONDS_PER_YEAR = 3.15576E7
YEARS_PER_MYR = 1.0E6

# erg/s -> (10^40 erg)/Myr: multiply by seconds/Myr, then divide by 10^40.
RATE_UNIT_LABEL = r"rate  $\left[10^{40}\ \mathrm{erg}\,\mathrm{Myr}^{-1}\right]$"
RATE_CONVERT = (SECONDS_PER_YEAR * YEARS_PER_MYR) / 1.0E40

# Validated categorical palette (light mode), fixed assignment order --
# the dataviz skill's own reference/palette.md 8-slot default. Assignment
# is fixed per named quantity (not per column index) so a given series
# keeps the same color/style regardless of which columns happen to be
# present in a given data file.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
DASH_STYLES = ["-", "--", ":", "-."]

# (column name, LaTeX legend label, color, style) -- checked/plotted in
# this fixed order. This is the exact, and only, quantity list the user
# asked for: dE_pol/dt, dE_tor/dt, Joule rate, Poynting rate, Hall-
# Poynting rate. Nothing else in the data file is plotted.
SERIES = [
    ("EDOT_poloidal_erg_per_s",           r"$\dot{E}_{B,\mathrm{pol}}$",   PALETTE[0], DASH_STYLES[0]),
    ("EDOT_toroidal_erg_per_s",           r"$\dot{E}_{B,\mathrm{tor}}$",   PALETTE[1], DASH_STYLES[1]),
    ("joule_dissipation_rate_erg_per_s",  r"$\dot{E}_{J}$",                PALETTE[2], DASH_STYLES[2]),
    ("poynting_flux_rate_erg_per_s",      r"$\dot{E}_{S}$",                PALETTE[3], DASH_STYLES[3]),
    ("hall_poynting_flux_rate_erg_per_s", r"$\dot{E}_{S,\mathrm{Hall}}$",  PALETTE[4], DASH_STYLES[0]),
]

# Total energy-rate balance -- literally the sum of all five contributions
# above, signed so a perfectly balanced budget gives ~0: the field's own
# rate of change (EDOT_poloidal + EDOT_toroidal) should equal everything
# draining/feeding it (Joule dissipation + the two Poynting flux terms).
# This is the exact same quantity the Fortran driver itself already
# computes as its residual check (app/mhdvsh_hall_crust_profile.f90's
# RESIDUAL = (EDOT_POL+EDOT_TOR) - (EDOT_J+EDOT_S+EDOT_S_HALL)) --
# recomputed here from the already-plotted columns rather than read from
# a separate "energy_balance_residual_erg_per_s" column, so the plot
# script stays self-documenting about what E_tot actually is, and still
# works against a data file that predates that column. No absolute value
# is taken anywhere in this sum -- it is a signed net balance, not a
# magnitude -- so no |...| notation is needed on its legend label.
EDOT_TOT_COLUMNS = [col for col, *_ in SERIES]
EDOT_TOT_SIGNS = [+1, +1, -1, -1, -1]
EDOT_TOT_LABEL = r"$\dot{E}_{\mathrm{tot}}$"
EDOT_TOT_COLOR = PALETTE[5]
EDOT_TOT_STYLE = DASH_STYLES[1]

TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID = "#e3e2dd"


def read_columns(path):
    # The column-header line is the LAST "#" line immediately before the
    # data starts, not the first -- a driver may write extra metadata/
    # provenance comment lines above it (e.g. app/mhdvsh_hall_crust_
    # profile.f90's own "run inputs:" block), and taking the first "#"
    # line in that case grabs a metadata line instead of the real header.
    header_line = None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                header_line = line
            else:
                break
    names = header_line.lstrip("#").split()
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return names, data


def derive_edot(names, data):
    """Backward-difference dE/dt from raw E_poloidal_erg/E_toroidal_erg,
    for data files logged before the Fortran driver exposed EDOT_* columns
    directly. Matches the driver's own EDOT_POL = (E_POL-PREV_E_POL)/DT_NOW
    formula exactly (same t_yr/dt_yr pairing) when a dt_yr column exists
    (adaptive-dt drivers). Fixed-dt drivers (e.g. the diffusion regime's
    own energy_budget.dat, which never logs a per-row dt_yr at all) fall
    back to np.diff(t_yr) -- equivalent since dt is constant there."""
    t_yr = data[:, names.index("t_yr")]
    if "dt_yr" in names:
        dt_s = data[:, names.index("dt_yr")] * SECONDS_PER_YEAR
    else:
        dt_s = np.empty_like(t_yr)
        dt_s[1:] = np.diff(t_yr) * SECONDS_PER_YEAR
    e_pol = data[:, names.index("E_poloidal_erg")]
    e_tor = data[:, names.index("E_toroidal_erg")]
    edot_pol = np.full_like(e_pol, np.nan)
    edot_tor = np.full_like(e_tor, np.nan)
    edot_pol[1:] = (e_pol[1:] - e_pol[:-1]) / dt_s[1:]
    edot_tor[1:] = (e_tor[1:] - e_tor[:-1]) / dt_s[1:]
    names = names + ["EDOT_poloidal_erg_per_s", "EDOT_toroidal_erg_per_s"]
    data = np.column_stack([data, edot_pol, edot_tor])
    return names, data


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png>")
    data_path, out_path = sys.argv[1], sys.argv[2]

    names, data = read_columns(data_path)
    if "EDOT_poloidal_erg_per_s" not in names:
        names, data = derive_edot(names, data)

    t = data[:, names.index("t_yr")]

    present = [(col, label, color, style) for col, label, color, style in SERIES if col in names]
    missing = [col for col, *_ in SERIES if col not in names]
    if missing:
        print(f"warning: columns not found, skipping: {missing}", file=sys.stderr)

    path_parts = Path(data_path).parts
    if "artifacts" in path_parts:
        suptitle = f"{path_parts[path_parts.index('artifacts') + 1].capitalize()} regime: energy budget vs time"
    else:
        suptitle = f"Energy budget vs time: {Path(data_path).stem}"

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    all_y = []
    for col, label, color, style in present:
        y = data[:, names.index(col)] * RATE_CONVERT
        ax.plot(t, y, color=color, linestyle=style, linewidth=2,
                 solid_capstyle="round", dash_capstyle="round", label=label)
        all_y.append(y)

    if all(col in names for col in EDOT_TOT_COLUMNS):
        edot_tot = sum(sign * data[:, names.index(col)] for sign, col in zip(EDOT_TOT_SIGNS, EDOT_TOT_COLUMNS))
        edot_tot = edot_tot * RATE_CONVERT
        ax.plot(t, edot_tot, color=EDOT_TOT_COLOR, linestyle=EDOT_TOT_STYLE, linewidth=2,
                 solid_capstyle="round", dash_capstyle="round", label=EDOT_TOT_LABEL)
        all_y.append(edot_tot)
    else:
        print("warning: can't compute E_tot, missing one of: "
              f"{[c for c in EDOT_TOT_COLUMNS if c not in names]}", file=sys.stderr)

    # Log-x (age), linear-y in the new 10^40 erg/Myr units, per user
    # direction (2026-08-28): spreads the fast initial relaxation
    # transient (order-years) out against the slow ~1000yr plateau decay
    # on the time axis, standard convention for NS field-decay plots.
    # The very first row (t=0) has no defined EDOT_* (nothing to
    # backward-difference against) and is dropped automatically by the
    # NaN there; log(0) would be undefined regardless.
    ax.set_xscale("log")

    # Clip y-limits rather than the data itself, per user direction
    # (2026-08-28): every point is still plotted (a brief initial-
    # relaxation transient, ~10^2-10^3x the eventual plateau, runs off
    # the top/bottom of the visible frame) but the axis range is set from
    # the 2nd/98th percentile of the plotted values with 15% padding, so
    # the quasi-steady-state structure the rest of the run actually shows
    # stays readable instead of being flattened by that one transient.
    all_y = np.concatenate(all_y)
    all_y = all_y[np.isfinite(all_y)]
    ylo, yhi = np.percentile(all_y, [2, 98])
    pad = 0.15 * (yhi - ylo)
    ax.set_ylim(ylo - pad, yhi + pad)

    ax.set_ylabel(RATE_UNIT_LABEL, color=TEXT_SECONDARY)
    ax.set_xlabel("t [yr]", color=TEXT_SECONDARY)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.axhline(0.0, color=GRID, linewidth=1.0, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY)
    ax.legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=10)

    fig.suptitle(suptitle, color=TEXT_PRIMARY, fontsize=12, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
