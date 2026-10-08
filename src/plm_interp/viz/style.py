"""Shared figure style: one palette, thin marks, recessive axes."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Categorical slots in fixed order (validated reference palette, light mode)
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
VIOLET = "#4a3aa7"
GRAY = "#8a8985"
LIGHT_GRAY = "#c9c8c3"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
SURFACE = "#fcfcfb"

MODEL_COLORS = {"8M": ORANGE, "35M": BLUE}

SEQUENTIAL = LinearSegmentedColormap.from_list(
    "seq_blue", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]
)
DIVERGING = LinearSegmentedColormap.from_list(
    "div_blue_red", ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f3a5a4", "#e34948", "#a8211f"]
)


def apply() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "axes.labelcolor": TEXT_2,
            "axes.edgecolor": LIGHT_GRAY,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#ecebe7",
            "grid.linewidth": 0.6,
            "xtick.color": TEXT_2,
            "ytick.color": TEXT_2,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.frameon": False,
            "legend.fontsize": 8,
            "lines.linewidth": 2.0,
            "lines.markersize": 4,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
        }
    )


def save(fig, path_stem) -> None:
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_stem.with_suffix(".png"))
    fig.savefig(path_stem.with_suffix(".pdf"))
    plt.close(fig)
