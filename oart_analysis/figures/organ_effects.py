"""Figure S7: organ-level mean paired changes and pointwise 95% CIs."""

import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]

CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]

DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]

SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]

ACCENT_RED  = "#B2182B"

GREY        = "#999999"

BLACK       = "#222222"

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Export exact-size vectors plus 300-dpi PNG and 600-dpi TIFF."""
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(f"{filename}.pdf", bbox_inches=None, dpi=300)
        fig.savefig(f"{filename}.svg", bbox_inches=None, dpi=300)
        fig.savefig(f"{filename}.png", bbox_inches=None, dpi=300)
        fig.savefig(
            f"{filename}.tiff", bbox_inches=None, dpi=600,
            pil_kwargs={"compression": "tiff_lzw"},
        )

import argparse

import csv

import itertools

import json

import math

import platform

from datetime import datetime, timezone

from pathlib import Path

import warnings

mpl.use("Agg")

import matplotlib.pyplot as plt

from matplotlib import font_manager

from matplotlib.lines import Line2D

FIGURE_NAME = "FigureS7_OAR_change_forest"

CENTRES = ("Main centre", "Centre I", "Centre II", "Centre III", "Centre IV")

CENTRE_SHORT = ("Main centre", "Centre I", "Centre II", "Centre III", "Centre IV")

SAMPLE_SIZES = {"Main centre": 29, "Centre I": 14, "Centre II": 21, "Centre III": 14, "Centre IV": 25}

COLOURS = {
    "Main centre": BLACK, "Centre I": CATEGORICAL[0],
    "Centre II": CATEGORICAL[2], "Centre III": CATEGORICAL[4],
    "Centre IV": CATEGORICAL_EXTENDED[9],
}

MARKERS = {"Main centre": "s", "Centre I": "o", "Centre II": "o", "Centre III": "o", "Centre IV": "o"}

METRICS = ("DSC", "HD95", "ASD")

METRIC_TITLES = {
    "DSC": "DSC improvement",
    "HD95": "HD95 reduction (mm)",
    "ASD": "ASD reduction (mm)",
}

PANEL_LETTERS = {"DSC": "a", "HD95": "b", "ASD": "c"}

SMALL_OARS = (
    "Optic chiasm", "Cochlea", "Eustachian tube bone", "Internal auditory canal",
    "Lens", "Optic nerve", "Pituitary", "Temporomandibular joint",
    "Tympanic cavity", "Vestibule/semicircular canals",
)

LARGE_OARS = (
    "Eye", "Glottic larynx", "Supraglottic larynx", "Mandible", "Oral cavity",
    "Parotid gland", "Brain stem", "Thyroid", "Temporal lobe",
    "Submandibular gland", "Spinal cord",
)

OARS = SMALL_OARS + LARGE_OARS

SCALES = {
    "DSC": {
        "main": (-0.01, 0.06, (0.00, 0.03, 0.06)),
        "Small OARs": (-0.10, 0.55, (-0.10, 0.00, 0.50)),
        "Large OARs": (-0.10, 0.25, (-0.10, 0.00, 0.20)),
    },
    "HD95": {
        "main": (-0.10, 1.10, (0.0, 0.5, 1.0)),
        "Small OARs": (-3.0, 9.0, (-3.0, 0.0, 9.0)),
        "Large OARs": (-2.0, 2.5, (-2.0, 0.0, 2.0)),
    },
    "ASD": {
        "main": (-0.02, 0.23, (0.0, 0.1, 0.2)),
        "Small OARs": (-2.0, 8.0, (-2.0, 0.0, 8.0)),
        "Large OARs": (-1.0, 1.0, (-1.0, 0.0, 1.0)),
    },
}

SOURCE_FIELDS = (
    "centre", "roi_group", "roi", "metric", "unit", "n", "before_model",
    "before_mean", "before_sd", "before_min", "before_max", "before_q1",
    "before_median", "before_q3", "after_model", "after_mean", "after_sd",
    "after_min", "after_max", "after_q1", "after_median", "after_q3",
    "mean_difference", "sd_paired_difference", "ci_lower", "ci_upper",
    "p_value", "q_value", "source_sheet", "source_row", "difference_cell", "q_cell",
)

NUMERIC_FIELDS = (
    "before_mean", "before_sd", "before_min", "before_max", "before_q1",
    "before_median", "before_q3", "after_mean", "after_sd", "after_min",
    "after_max", "after_q1", "after_median", "after_q3", "mean_difference",
    "sd_paired_difference", "ci_lower", "ci_upper",
)

def parse_probability(value):
    text = str(value).strip().lstrip("<").strip()
    number = float(text)
    if not 0 <= number <= 1:
        raise ValueError(f"Probability outside [0,1]: {value}")
    return number

def validate_data(rows):
    expected = set(itertools.product(CENTRES, OARS, METRICS))
    keys = [(r["centre"], r["roi"], r["metric"]) for r in rows]
    if len(keys) != 315 or len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError("Table S5 coverage is incomplete or duplicated")
    for row in rows:
        key = (row["centre"], row["roi"], row["metric"])
        expected_group = "Small OARs" if row["roi"] in SMALL_OARS else "Large OARs"
        if row["roi_group"] != expected_group:
            raise ValueError(f"Wrong OAR group: {key}")
        if row["n"] != SAMPLE_SIZES[row["centre"]]:
            raise ValueError(f"Wrong paired n: {key}")
        if not all(math.isfinite(row[field]) for field in NUMERIC_FIELDS):
            raise ValueError(f"Non-finite number: {key}")
        if not row["ci_lower"] <= row["mean_difference"] <= row["ci_upper"]:
            raise ValueError(f"CI excludes estimate: {key}")
    return {"rows": len(rows), "unique_comparisons": len(set(keys)), "omitted": 0}

def oriented_effect(row):
    """Return estimate and CI with positive values consistently meaning better."""
    estimate = row["mean_difference"]
    lower, upper = row["ci_lower"], row["ci_upper"]
    if row["metric"] in ("HD95", "ASD"):
        estimate, lower, upper = -estimate, -upper, -lower
    return estimate, lower, upper

def select_font(requested=None):
    candidates = ([requested] if requested else []) + ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
    for candidate in candidates:
        try:
            font_manager.findfont(candidate, fallback_to_default=False)
        except ValueError:
            continue
        mpl.rcParams["font.sans-serif"] = [candidate]
        if requested and requested != candidate:
            warnings.warn(f"{requested!r} unavailable; using {candidate!r}")
        return candidate
    raise RuntimeError("No supported sans-serif font is available")

def tick_text(value, metric, scale_kind):
    if metric == "DSC":
        if abs(value) < 1e-12:
            return "0"
        return f"{value:.2f}"
    if metric == "ASD" and scale_kind == "main":
        return "0" if abs(value) < 1e-12 else f"{value:.1f}"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.1f}"

def axis_for(fig, left_mm, bottom_mm, width_mm, height_mm, total_width, total_height):
    return fig.add_axes([
        left_mm / total_width, bottom_mm / total_height,
        width_mm / total_width, height_mm / total_height,
    ])

def draw_facet(ax, rows, centre, metric, group, roi_order, show_labels):
    selected = {
        row["roi"]: row for row in rows
        if row["centre"] == centre and row["metric"] == metric and row["roi_group"] == group
    }
    if set(selected) != set(roi_order):
        raise ValueError(f"Incomplete facet: {centre}, {metric}, {group}")
    scale_kind = "main" if centre == "Main centre" else group
    lo, hi, ticks = SCALES[metric][scale_kind]
    y = list(range(len(roi_order)))

    for yi in y:
        ax.hlines(yi, lo, hi, color="#D8DDE1", linewidth=0.35, zorder=0)
    ax.axvline(0, color="#6F767C", linewidth=0.7, zorder=1)

    for yi, roi in zip(y, roi_order):
        row = selected[roi]
        estimate, lower, upper = oriented_effect(row)
        if lower < lo or upper > hi:
            raise ValueError(f"CI clipped by axis: {centre}, {metric}, {group}, {roi}, {(lower, upper)}")
        colour = COLOURS[centre]
        ax.plot([lower, upper], [yi, yi], color="#BFC4C8", linewidth=1.05, solid_capstyle="round", zorder=2)
        significant = row["q_numeric"] < 0.05
        ax.plot(
            estimate, yi, marker=MARKERS[centre], markersize=3.5,
            markerfacecolor=colour if significant else "white",
            markeredgecolor=colour, markeredgewidth=0.85,
            linestyle="none", zorder=3,
        )

    ax.set_xlim(lo, hi)
    ax.set_ylim(len(roi_order) - 0.5, -0.5)
    ax.set_xticks(ticks)
    ax.set_xticklabels([tick_text(v, metric, scale_kind) for v in ticks], fontsize=5.0)
    ax.tick_params(axis="x", length=2.0, pad=1.0, width=0.5, color="#555555")
    ax.set_yticks(y)
    if show_labels:
        ax.set_yticklabels(roi_order, fontsize=5.1)
        ax.tick_params(axis="y", length=0, pad=1.5)
        for label in ax.get_yticklabels():
            label.set_rotation(15)
            label.set_rotation_mode("anchor")
            label.set_horizontalalignment("right")
            label.set_verticalalignment("center")
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#6F767C")
    ax.spines["bottom"].set_linewidth(0.55)
    return {
        "centre": centre, "metric": metric, "group": group,
        "axis_min": lo, "axis_max": hi, "ticks": list(ticks),
        "n_rows": len(roi_order), "significant": sum(selected[r]["q_numeric"] < 0.05 for r in roi_order),
    }

def export_plotted_intervals(rows, path):
    fields = (
        "panel", "centre", "roi_group", "roi", "metric", "unit", "estimate_positive_is_improvement",
        "ci_lower_positive_is_improvement", "ci_upper_positive_is_improvement",
        "bh_q_original", "bh_q_numeric", "significant_bh_q_lt_0_05",
        "source_row", "difference_cell", "q_cell",
    )
    with Path(path).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            estimate, lower, upper = oriented_effect(row)
            writer.writerow({
                "panel": PANEL_LETTERS[row["metric"]], "centre": row["centre"],
                "roi_group": row["roi_group"], "roi": row["roi"], "metric": row["metric"],
                "unit": row["unit"], "estimate_positive_is_improvement": estimate,
                "ci_lower_positive_is_improvement": lower, "ci_upper_positive_is_improvement": upper,
                "bh_q_original": row["q_value"], "bh_q_numeric": row["q_numeric"],
                "significant_bh_q_lt_0_05": row["q_numeric"] < 0.05,
                "source_row": row["source_row"], "difference_cell": row["difference_cell"],
                "q_cell": row["q_cell"],
            })

def draw_figure(rows, output_dir):
    width_mm, height_mm = 183.0, 240.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")
    lefts = (41.0, 66.0, 93.0, 120.0, 147.0)
    widths = (20.0, 24.0, 24.0, 24.0, 24.0)
    section_bottoms = {"DSC": 163.0, "HD95": 90.0, "ASD": 17.0}
    axes_report = []

    for metric in METRICS:
        bottom = section_bottoms[metric]
        # Panel title and direction statement.
        fig.text(4 / width_mm, (bottom + 71.0) / height_mm, PANEL_LETTERS[metric],
                 fontsize=8.5, fontweight="bold", ha="left", va="center", color=BLACK)
        fig.text(9 / width_mm, (bottom + 71.0) / height_mm, METRIC_TITLES[metric],
                 fontsize=8.0, fontweight="bold", ha="left", va="center", color=BLACK)
        fig.text(179 / width_mm, (bottom + 71.0) / height_mm,
                 "Mean paired difference (95% CI); positive favours locked/adapted  →",
                 fontsize=5.2, ha="right", va="center", color="#555555")

        # Vertical group labels preserve space for the 15° ROI labels and avoid
        # collisions at the boundary between the two OAR-size blocks.
        fig.text(3.2 / width_mm, (bottom + 49.5) / height_mm, "Small OARs (≤4 cm³)",
                 rotation=90, fontsize=5.7, fontweight="bold", ha="center", va="center", color="#555555")
        fig.text(3.2 / width_mm, (bottom + 18.0) / height_mm, "Large OARs (>4 cm³)",
                 rotation=90, fontsize=5.7, fontweight="bold", ha="center", va="center", color="#555555")
        # Centre headings are repeated for rapid scanning.
        for centre, label, left, width in zip(CENTRES, CENTRE_SHORT, lefts, widths):
            fig.text((left + width / 2) / width_mm, (bottom + 64.5) / height_mm,
                     label, fontsize=5.8, fontweight="bold", ha="center", va="center", color=COLOURS[centre])
            fig.text((left + width / 2) / width_mm, (bottom + 61.7) / height_mm,
                     f"n={SAMPLE_SIZES[centre]}", fontsize=5.1, ha="center", va="center", color="#555555")

        # Main-centre divider emphasises the different comparison while external centres remain comparable.
        fig.add_artist(Line2D(
            [63.5 / width_mm, 63.5 / width_mm],
            [(bottom + 2.0) / height_mm, (bottom + 61.5) / height_mm],
            transform=fig.transFigure, color="#C4C9CD", linewidth=0.55,
        ))

        for centre_index, (centre, left, width) in enumerate(zip(CENTRES, lefts, widths)):
            small_ax = axis_for(fig, left, bottom + 38.0, width, 23.0, width_mm, height_mm)
            large_ax = axis_for(fig, left, bottom + 4.0, width, 28.0, width_mm, height_mm)
            axes_report.append(draw_facet(
                small_ax, rows, centre, metric, "Small OARs", SMALL_OARS, centre_index == 0
            ))
            axes_report.append(draw_facet(
                large_ax, rows, centre, metric, "Large OARs", LARGE_OARS, centre_index == 0
            ))

        if metric != "ASD":
            fig.add_artist(Line2D(
                [4 / width_mm, 179 / width_mm],
                [(bottom - 1.5) / height_mm, (bottom - 1.5) / height_mm],
                transform=fig.transFigure, color="#D6DADD", linewidth=0.45,
            ))

    # Direct statistical key and scale disclosure.
    fig.add_artist(Line2D([7 / width_mm], [10.4 / height_mm], marker="o", markersize=3.6,
                          markerfacecolor=CATEGORICAL[0], markeredgecolor=CATEGORICAL[0],
                          linestyle="none", transform=fig.transFigure))
    fig.text(10 / width_mm, 10.4 / height_mm, "BH-adjusted q < 0.05", fontsize=5.4,
             ha="left", va="center", color=BLACK)
    fig.add_artist(Line2D([45 / width_mm], [10.4 / height_mm], marker="o", markersize=3.6,
                          markerfacecolor="white", markeredgecolor=CATEGORICAL[0],
                          markeredgewidth=0.85, linestyle="none", transform=fig.transFigure))
    fig.text(48 / width_mm, 10.4 / height_mm, "q ≥ 0.05", fontsize=5.4,
             ha="left", va="center", color=BLACK)
    fig.add_artist(Line2D([68 / width_mm, 75 / width_mm], [10.4 / height_mm, 10.4 / height_mm],
                          color="#BFC4C8", linewidth=1.05, transform=fig.transFigure))
    fig.text(77 / width_mm, 10.4 / height_mm, "95% CI", fontsize=5.4,
             ha="left", va="center", color=BLACK)
    fig.text(179 / width_mm, 10.4 / height_mm,
             "Main centre: square; external centres: circle", fontsize=5.2,
             ha="right", va="center", color="#555555")
    fig.text(4 / width_mm, 5.5 / height_mm,
             "Main centre uses one zoomed scale per metric; Centres I–IV share a scale within each OAR-size block. All limits cover the complete 95% CIs.",
             fontsize=5.1, ha="left", va="center", color="#555555")
    fig.text(4 / width_mm, 2.3 / height_mm,
             "Main: locked source model − conventional U-Net. External centres: locally adapted − pre-adaptation (distance reductions are sign-reversed).",
             fontsize=5.1, ha="left", va="center", color="#555555")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    export_plotted_intervals(rows, output_dir / "plotted_intervals.csv")
    save_cns_figure(fig, output_dir / FIGURE_NAME)
    plt.close(fig)
    return axes_report
