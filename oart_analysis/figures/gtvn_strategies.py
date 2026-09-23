"""Figure S4b–d: GTVn strategy comparison, all patient pairs, no inferential tests."""

from __future__ import annotations

import argparse

import csv

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import numpy as np

from matplotlib.lines import Line2D

from PIL import Image

plt.rcParams["font.family"] = "sans-serif"

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]

plt.rcParams["svg.fonttype"] = "none"

plt.rcParams["pdf.fonttype"] = 42

plt.rcParams.update(
    {
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.3,
        "axes.linewidth": 0.65,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    }
)

STRATEGIES = (
    "Image-only baseline",
    "Direct DIR propagation",
    "PCG-UNet",
)

METRICS = ("DSC", "HD95", "ASD")

SHORT_METHODS = ("Baseline", "DIR", "PCG-UNet")

BOX_COLORS = ("#D8DCE4", "#ECDCC1", "#B9D0E6")

POINT_COLORS = ("#53606F", "#A47232", "#15518F")

METRIC_SPECS = {
    "DSC": {"title": "DSC", "ylim": (0.55, 0.86), "yticks": (0.60, 0.70, 0.80)},
    "HD95": {
        "title": "HD95 (mm)",
        "ylim": (3.75, 8.60),
        "yticks": (4, 5, 6, 7, 8),
    },
    "ASD": {"title": "ASD (mm)", "ylim": (0.72, 3.17), "yticks": (1, 1.5, 2, 2.5, 3)},
}

def draw_metric(ax, metric: str, patients: list[str], values: dict) -> None:
    """Draw boxplots, all individual points, and subtle patient-pairing lines."""
    data = [np.array([values[(p, s, metric)] for p in patients]) for s in STRATEGIES]
    # Fixed offsets make all runs deterministic and preserve the same patient's
    # horizontal position in every strategy and metric panel.
    offsets = np.linspace(-0.11, 0.11, len(patients))
    for i, patient in enumerate(patients):
        ax.plot(
            np.arange(1, 4) + offsets[i],
            [values[(patient, s, metric)] for s in STRATEGIES],
            color="#B8BEC6",
            lw=0.50,
            alpha=0.40,
            zorder=0.5,
        )

    boxes = ax.boxplot(
        data,
        positions=(1, 2, 3),
        widths=0.52,
        whis=1.5,
        patch_artist=True,
        showfliers=False,  # All observations are overlaid as points instead.
        boxprops={"color": "#59616B", "linewidth": 0.75},
        whiskerprops={"color": "#59616B", "linewidth": 0.75},
        capprops={"color": "#59616B", "linewidth": 0.75},
        medianprops={"color": "#20272F", "linewidth": 1.15},
    )
    for patch, color in zip(boxes["boxes"], BOX_COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.76)
        patch.set_zorder(1)

    for category, (series, color) in enumerate(zip(data, POINT_COLORS), start=1):
        ax.scatter(
            category + offsets,
            series,
            s=15,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=0.94,
            zorder=3,
        )

    spec = METRIC_SPECS[metric]
    ax.set_title(spec["title"], pad=7, fontweight="semibold", color="#252B32")
    ax.set_xlim(0.57, 3.43)
    ax.set_ylim(*spec["ylim"])
    ax.set_yticks(spec["yticks"])
    ax.set_xticks((1, 2, 3))
    ax.set_xticklabels(SHORT_METHODS, rotation=24, ha="right")
    ax.tick_params(axis="both", direction="out", length=2.6, width=0.65, color="#626973")
    ax.spines["left"].set_color("#626973")
    ax.spines["bottom"].set_color("#626973")
    ax.grid(axis="y", color="#E8EBEE", lw=0.45, zorder=0)

def create_lower_figure(patients: list[str], values: dict):
    """Create the lower quantitative row as a standalone vector figure."""
    fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 78 / 25.4))
    fig.subplots_adjust(left=0.075, right=0.985, top=0.83, bottom=0.20, wspace=0.37)
    for letter, ax, metric in zip("bcd", axes, METRICS):
        draw_metric(ax, metric, patients, values)
        ax.text(-0.20, 1.10, letter, transform=ax.transAxes, fontsize=9, fontweight="bold")
    fig.text(0.075, 0.965, "12 paired cases", ha="left", va="top", fontsize=7.5, color="#3D4650")
    return fig

def save_figure(fig, stem: Path, raster_dpi: int) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.04)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.04)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(stem.with_suffix(".tiff"), dpi=raster_dpi, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
