#!/usr/bin/env python3
"""EOS-derived quantities and grid spacing on the ACTUAL simulation grid
(N_R=50, theta_max=0.2 truncation, RADIAL_GEOMETRIC_OUTER Q=0.90 -- see
app/mhdvsh_hall_crust_profile.f90), evaluated at every grid node (not a
dense reference profile) -- one panel per quantity, each with its own
y-scale, grid nodes marked explicitly as points so the reader can see
directly how the non-uniform grid resolves (or doesn't resolve) each
quantity's radial variation.

Usage: plot_eos_profile_truncated.py <data_file> <output_png>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"
GRID_COLOR = "#e3e2dd"


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
    ax.grid(True, color=GRID_COLOR, linewidth=0.8, which="both")
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID_COLOR)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <data_file> <output_png>")
    data_path, out_path = sys.argv[1], sys.argv[2]

    names, data = read_columns(data_path)
    r = data[:, names.index("r_km")]
    dr = data[:, names.index("dr_km")]
    rho = data[:, names.index("rho_gcc")]
    press = data[:, names.index("pressure_dyncm2")]
    nb = data[:, names.index("n_baryon_fm-3")]
    ne = data[:, names.index("n_e_cgs")]
    eta = data[:, names.index("eta_codeunits_km2yr")]
    f_hall = data[:, names.index("f_hall_codeunits")]
    theta = data[:, names.index("theta")]
    br = data[:, names.index("B_r_1e12G")]
    rmag = data[:, names.index("Rmag")]

    r_min, r_max = r[0], r[-1]

    # (title, ylabel, values, color)
    panels = [
        (r"Grid spacing $dr(r)$ -- the non-uniform grid itself", "dr  [km]", dr, PALETTE[6]),
        (r"Mass density $\rho(r)$", r"$\rho$  [g/cm$^3$]", rho, PALETTE[0]),
        (r"Pressure $P(r)$", r"$P$  [dyn/cm$^2$]", press, PALETTE[1]),
        (r"Baryon number density $n_b(r)$", r"$n_b$  [fm$^{-3}$]", nb, PALETTE[2]),
        (r"Electron number density $n_e(r)$", r"$n_e$  [cm$^{-3}$]", ne, PALETTE[3]),
        (r"Resistivity $\eta(r)$ (code units)", r"$\eta$  [km$^2$/yr]", eta, PALETTE[4]),
        (r"Hall prefactor $f_H(r)$ (code units)", r"$f_H$", f_hall, PALETTE[5]),
        (r"Degeneracy parameter $\theta(r)=k_BT/E_F$ -- truncation criterion, "
         r"cut at $\theta_{\max}=0.2$ (the $r_{\max}$ edge)", r"$\theta$", theta, PALETTE[7]),
        (r"Radial field at the pole $B_r(r)$ (l=1 force-free IC)", r"$B_r$  [$10^{12}$ G]", br, PALETTE[6]),
        (r"Local Hall magnetic Reynolds number $R_{\rm mag}(r)=f_H(r)B_r(r)/\eta(r)$",
         r"$R_{\rm mag}$", rmag, PALETTE[0]),
    ]

    fig, axes = plt.subplots(len(panels), 1, figsize=(7.5, 2.8 * len(panels)), dpi=150)
    fig.patch.set_facecolor(SURFACE)

    for ax, (title, ylabel, y, color) in zip(axes, panels):
        ax.plot(r, y, color=color, linewidth=1.4, zorder=2)
        ax.scatter(r, y, color=color, s=14, zorder=3, edgecolors=SURFACE, linewidths=0.6,
                   label="grid node")
        ax.set_yscale("log")
        style_axis(ax)
        ax.set_title(title, color=TEXT_PRIMARY, fontsize=10, loc="left")
        ax.set_ylabel(ylabel, color=TEXT_SECONDARY, fontsize=9)
        ax.set_xlabel("r  [km]", color=TEXT_SECONDARY, fontsize=8)
        ax.set_xlim(r_min, r_max)
        ax.legend(loc="best", fontsize=7, framealpha=0.9, facecolor=SURFACE, edgecolor=GRID_COLOR,
                  labelcolor=TEXT_SECONDARY)

    theta_ax, rmag_ax = axes[7], axes[9]  # indices into `panels` above
    theta_ax.axhline(0.2, color=TEXT_SECONDARY, linewidth=0.9, linestyle="--", zorder=1)
    theta_ax.text(r_min, 0.2, r"  $\theta_{\max}=0.2$", color=TEXT_SECONDARY, fontsize=7,
                  va="bottom", ha="left")
    rmag_ax.text(r[-1], rmag[-1], f"  $R_{{\\rm mag}}(r_{{\\max}})$={rmag[-1]:.3g}",
                 color=TEXT_SECONDARY, fontsize=7, va="center", ha="left")

    fig.suptitle(
        f"EOS parameters on the truncated simulation grid "
        f"(N_R=50, $\\theta_{{\\max}}$=0.2, $r$=[{r_min:.4f}, {r_max:.4f}] km)",
        color=TEXT_PRIMARY, fontsize=12, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(out_path, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
