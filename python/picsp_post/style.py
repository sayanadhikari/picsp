"""Consistent figure style for PICSP3D plots.

PICSP_STYLE=paper switches to a journal style (Nature-like): figures resized to the printed
column width, 7 pt Helvetica-class type, bold panel letters in place of axes titles and thin
axes. `journal_layout(fig)` applies it to a finished figure and returns the removed titles.
"""
import os
import re

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
        "savefig.dpi": int(os.environ.get("PICSP_FIG_DPI", 160)),
        "savefig.bbox": "tight",
        "image.cmap": CMAP_SEQ,
    })


def series_color(i):
    return SERIES[i % len(SERIES)]


PAPER = os.environ.get("PICSP_STYLE", "") == "paper"
MM = 1 / 25.4
COLUMN = 89 * MM          # single column
DOUBLE = 183 * MM         # full page width
FONT = ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"]   # Helvetica and Arial metrics
SIZE, SMALL, LETTER = 7, 6.5, 8


def apply_paper():
    mpl.rcParams.update({
        "font.family": FONT,     # a list: glyphs missing in the first font come from the next
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT[0],
        "mathtext.it": f"{FONT[0]}:italic",
        "mathtext.bf": f"{FONT[0]}:bold",
        "mathtext.sf": FONT[0],
        "mathtext.fallback": "stixsans",
        "font.size": SIZE,
        "axes.titlesize": SIZE,
        "axes.labelsize": SIZE,
        "xtick.labelsize": SMALL,
        "ytick.labelsize": SMALL,
        "legend.fontsize": SMALL,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.edgecolor": INK,
        "axes.linewidth": 0.6,
        "axes.grid": False,
        "axes.labelpad": 2.5,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.pad": 2,
        "ytick.major.pad": 2,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "legend.handlelength": 1.6,
        "legend.borderaxespad": 0.3,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.pad_inches": 0.02,
        "axes.formatter.use_mathtext": True,
        "axes.formatter.limits": (-3, 4),
    })


# "v_x", "n_O+", "ω_pe", "B_max": subscripts written with an underscore in plain-text labels
_SUB = re.compile(r"(?<![\w$\\])([A-Za-z\u0370-\u03ff]{1,3})_([A-Za-z0-9+]+)")


_SUPER = dict(zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+"))
_SUPER_RUN = re.compile("[" + "".join(_SUPER) + "]+")
_MATH_SYMBOLS = {"⟨": r"$\langle $", "⟩": r"$\rangle $"}


def _subscripts(text):
    """Plain-text label -> mathtext: "v_x" becomes v with subscript x. Once a label holds
    mathtext, unicode superscripts and brackets are converted too, since mathtext takes its
    glyphs from one font only."""
    if "$" not in text and "_" in text:
        text = _SUB.sub(lambda m: rf"${m.group(1)}_{{\mathrm{{{m.group(2)}}}}}$", text)
    if "$" in text:
        text = _SUPER_RUN.sub(lambda m: "$^{" + "".join(_SUPER[c] for c in m.group(0)) + "}$", text)
        for k, v in _MATH_SYMBOLS.items():
            text = text.replace(k, v)
        text = text.replace("$$", "")
    return text


def _is_colorbar(ax):
    return hasattr(ax, "_colorbar") or ax.get_label() == "<colorbar>"


def journal_layout(fig, width=None):
    """Resize a finished figure to a journal column, replace axes titles by panel letters and
    restyle its text. Returns the removed titles as [(letter, title)]."""
    w, h = fig.get_size_inches()
    fixed = getattr(fig, "_journal_width", None)    # laid out at its printed size already
    width = width or fixed or (DOUBLE if w > 8.5 else COLUMN)
    if not fixed:
        fig.set_size_inches(width, max(1.6, h * width / w * 1.08))
    main, seen = [], []
    for ax in fig.axes:
        pos = tuple(round(v, 3) for v in ax.get_position().bounds)
        if _is_colorbar(ax) or pos in seen or not ax.get_visible():
            continue
        seen.append(pos)
        main.append(ax)
    if fig._suptitle is not None:
        fig._suptitle.set_visible(False)
    titles = []
    for ax in fig.axes:
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_linewidth(0.6)
            sp.set_color(INK)
        ax.tick_params(which="both", labelsize=SMALL, colors=INK, width=0.6, length=2.5, pad=2)
        for lab in (ax.xaxis.label, ax.yaxis.label):
            lab.set_text(_subscripts(lab.get_text()))
            lab.set_fontsize(SIZE)
            lab.set_color(INK)
        for t in ax.texts:
            t.set_text(_subscripts(t.get_text()))
            t.set_fontsize(min(t.get_fontsize(), SIZE))
        leg = ax.get_legend()
        if leg is not None:
            leg.get_frame().set_visible(False)
            for t in leg.get_texts():
                t.set_text(_subscripts(t.get_text()))
                t.set_fontsize(SMALL)
        for line in ax.get_lines():
            line.set_linewidth(min(line.get_linewidth(), 1.2))
    multi = len(main) > 1
    for i, ax in enumerate(main):
        title = ax.get_title()
        ax.set_title("")
        titles.append((chr(97 + i) if multi else "", title))
    if not fixed:
        try:
            fig.tight_layout(pad=0.3, w_pad=1.2, h_pad=1.0)
        except Exception:
            pass
    for ax in fig.axes:
        cb = getattr(ax, "_colorbar", None)
        if cb is not None and getattr(cb.mappable, "axes", None) is not None:
            ax.yaxis.label.set_text(_subscripts(ax.yaxis.label.get_text()))
            parent = cb.mappable.axes.get_position()
            pos = ax.get_position()
            ax.set_position([pos.x0, parent.y0, pos.width, parent.height])
    if multi:
        r = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        for (letter, _), ax in zip(titles, main):
            x0 = inv.transform(ax.get_tightbbox(r).p0)[0]
            y1 = ax.get_position().y1
            fig.text(x0, y1, letter, fontsize=LETTER, fontweight="bold", va="bottom", ha="left", color=INK)
    return titles


apply()
if PAPER:
    apply_paper()
