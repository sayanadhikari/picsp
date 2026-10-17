"""Consistent figure style for PICSP3D plots."""
import matplotlib as mpl
import matplotlib.pyplot as plt

# categorical series colours, assigned in fixed order (validated palette)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"

# colormaps by the job the colour does
CMAP_SEQ = "viridis"     # magnitudes (density, |E|)
CMAP_HOT = "inferno"     # temperatures, energies
CMAP_DIV = "RdBu_r"      # signed quantities (potential, charge density), centred at 0


def apply():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 10.5,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=SERIES),
        "lines.linewidth": 1.6,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "font.family": "DejaVu Sans",
        "figure.dpi": 110,
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
        "image.cmap": CMAP_SEQ,
    })


def series_color(i):
    return SERIES[i % len(SERIES)]


apply()
