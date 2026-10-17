#!/usr/bin/env python3
"""Every energy channel from production_run.dat, each on its OWN panel
with its own y-scale -- unlike plot_energy_budget.py (one shared,
percentile-clipped axis, which hides any channel much smaller than the
largest one, e.g. Joule/diffusive-Poynting next to a large Hall-Poynting
term). Built specifically to diagnose WHERE an energy-conservation
violation originates: E_pol/E_tor themselves, their own EDOT rates, each
dissipation/flux term individually, and the residual that should sum
them to zero.

Every panel has an explicit y-axis label stating exactly what quantity
and units are plotted -- do not rely on legend text alone (a single-
series panel never shows a legend, so that text would otherwise appear
nowhere on the figure at all).

The top panel (E_tot percent deviation from t=0) uses a LINEAR x-axis,
not log: t=0 is a load-bearing reference point there (exactly 0% by
construction) and a log x-axis cannot show t=0 at all (log(0) is
undefined), which silently hides that anchor and makes the first
VISIBLE point look like an unexplained jump. The remaining panels use
log-x (their own t=0 rows are NaN anyway, by construction -- EDOT/
residual need a prior sample to backward-difference against).

Usage: plot_energy_channels.py <data_file> <output_png>
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


def read_columns(path):
    header_line = None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                header_line = line
            else:
                break
    names = header_line.lstrip("#").split()
    data = np.loadtxt(path, comments="#")
    return names, data


def style_axis(ax):
    ax.set_facecolor(SURFACE)
    ax.axhline(0.0, color=GRID, linewidth=1.0, zorder=0)
    ax.grid(True, color=GRID, linewidth=0.8, which="both")
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png>")
    data_path, out_path = sys.argv[1], sys.argv[2]

    names, data = read_columns(data_path)
    t = data[:, names.index("t_yr")]
    e_pol = data[:, names.index("E_poloidal_erg")]
    e_tor = data[:, names.index("E_toroidal_erg")]
    e_tot = e_pol + e_tor
    e0 = e_tot[0]
    assert t[0] == 0.0, "expected the first logged row to be t=0"

    edot_pol = data[:, names.index("EDOT_poloidal_erg_per_s")]
    edot_tor = data[:, names.index("EDOT_toroidal_erg_per_s")]
    joule = data[:, names.index("joule_dissipation_rate_erg_per_s")]
    poynting = data[:, names.index("poynting_flux_rate_erg_per_s")]
    hall_poynting = data[:, names.index("hall_poynting_flux_rate_erg_per_s")]
    residual = data[:, names.index("energy_balance_residual_erg_per_s")]

    # (title, ylabel, [(y, legend_label, color), ...], symlog?, log-x?)
    panels = [
        (r"$E_{\rm tot}$ deviation from $t{=}0$ ($E_{\rm tot}(0)=%.4e$ erg)" % e0,
         "percent deviation  [%]",
         [(100*(e_tot-e0)/e0, None, PALETTE[5])], False, False),
        ("E_poloidal, E_toroidal", "energy  [erg]",
         [(e_pol, "E_pol", PALETTE[0]), (e_tor, "E_tor", PALETTE[1])], False, True),
        (r"$\dot E_{pol}$, $\dot E_{tor}$", "rate  [erg/s]",
         [(edot_pol, "EDOT_pol", PALETTE[0]), (edot_tor, "EDOT_tor", PALETTE[1])], True, True),
        (r"Joule dissipation $\dot E_J$", "rate  [erg/s]",
         [(joule, None, PALETTE[2])], True, True),
        (r"Diffusive Poynting $\dot E_S$", "rate  [erg/s]",
         [(poynting, None, PALETTE[3])], True, True),
        (r"Hall Poynting $\dot E_{S,Hall}$", "rate  [erg/s]",
         [(hall_poynting, None, PALETTE[4])], True, True),
        (r"Residual: $(\dot E_{pol}+\dot E_{tor}) - (\dot E_J+\dot E_S+\dot E_{S,Hall})$", "rate  [erg/s]",
         [(residual, None, PALETTE[7])], True, True),
    ]

    fig, axes = plt.subplots(len(panels), 1, figsize=(7.5, 3.0*len(panels)), dpi=150)
    fig.patch.set_facecolor(SURFACE)

    for ax, (title, ylabel, series, symlog, logx) in zip(axes, panels):
        for y, label, color in series:
            ax.plot(t, y, color=color, linewidth=1.8, label=label)
        if symlog:
            finite = series[0][0][np.isfinite(series[0][0])]
            linthresh = max(1.0, np.nanmax(np.abs(finite))*1e-6) if len(finite) else 1.0
            ax.set_yscale("symlog", linthresh=linthresh)
        if logx:
            ax.set_xscale("log")
        style_axis(ax)
        ax.set_title(title, color=TEXT_PRIMARY, fontsize=10, loc="left")
        ax.set_ylabel(ylabel, color=TEXT_SECONDARY, fontsize=9)
        ax.set_xlabel("t [yr]", color=TEXT_SECONDARY, fontsize=8)
        if any(label is not None for _, label, _ in series):
            ax.legend(frameon=False, labelcolor=TEXT_SECONDARY, fontsize=8, loc="upper left")

    run_name = Path(data_path).parent.name
    fig.suptitle(f"Energy channels (individually scaled): {run_name}", color=TEXT_PRIMARY,
                 fontsize=13, fontweight="bold", x=0.02, ha="left")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, facecolor=SURFACE)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
