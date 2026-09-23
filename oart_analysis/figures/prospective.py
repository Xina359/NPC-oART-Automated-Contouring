"""Figure 6: prospective concordance and ranked recorded time intervals."""

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
    """Export exact-size PDF/SVG plus 300-dpi PNG and 600-dpi TIFF."""
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

import hashlib

import itertools

import json

import math

import platform

from collections import defaultdict

from datetime import datetime, timezone

from pathlib import Path

import warnings

mpl.use("Agg")

import matplotlib.pyplot as plt

import numpy as np

from matplotlib import font_manager

from matplotlib.lines import Line2D

from matplotlib.patches import Patch

FIGURE_NAME = "Figure6_prospective_quality_and_workflow"

METRICS = ("DSC", "HD95", "ASD")

CATEGORIES = ("GTVp", "CTV1", "CTV2", "GTVn", "Large OARs", "Small OARs")

TARGETS = ("GTVp", "CTV1", "CTV2", "GTVn")

OAR_GROUPS = ("Large OARs", "Small OARs")

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

GROUP_ROIS = {"Small OARs": SMALL_OARS, "Large OARs": LARGE_OARS}

TARGET_BLUE = "#3B73B9"

OAR_GREEN = "#66BD7A"

ALGORITHM_TEAL = "#4C9FAD"

REVIEW_ORANGE = "#F2A65A"

ROW_Y = {"GTVp": 5.5, "CTV1": 4.5, "CTV2": 3.5, "GTVn": 2.5, "Large OARs": 1.0, "Small OARs": 0.0}

AXIS_LIMITS = {"DSC": (0.68, 1.00), "HD95": (0.0, 8.0), "ASD": (0.0, 5.3)}

AXIS_TICKS = {"DSC": (0.70, 0.80, 0.90, 1.00), "HD95": (0, 2, 4, 6, 8), "ASD": (0, 1, 2, 3, 4, 5)}

METRIC_TITLES = {"DSC": "DSC", "HD95": "HD95 (mm)", "ASD": "ASD (mm)"}

PERFORMANCE_FIELDS = (
    "study_id", "category", "category_type", "metric", "value", "unit",
    "n_components", "source_sheet", "source_rows", "source_cells", "derivation",
)

TIMING_FIELDS = (
    "study_id", "algorithm_s", "review_s", "combined_s", "source_sheet",
    "source_row", "algorithm_cell", "review_cell", "combined_cell",
)

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def read_source_workbook(path):
    """Read the two prospective sheets and derive patient-level OAR-group means."""
    from openpyxl import load_workbook

    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb["Prospective patient metrics"]
        values = list(ws.values)
        header = list(values[0])
        required = {"Study ID", "ROI group", "ROI", "Metric", "Value", "Unit"}
        if not required <= set(header):
            raise ValueError(f"Missing prospective metric columns: {required-set(header)}")
        raw = []
        for excel_row, cells in enumerate(values[1:], 2):
            if not any(v is not None for v in cells):
                continue
            record = dict(zip(header, cells))
            record["_row"] = excel_row
            raw.append(record)

        study_ids = sorted({str(r["Study ID"]) for r in raw})
        performance = []
        lookup = defaultdict(list)
        for row in raw:
            lookup[(str(row["Study ID"]), str(row["Metric"]))].append(row)
        for study_id in study_ids:
            for metric in METRICS:
                patient_rows = lookup[(study_id, metric)]
                if len(patient_rows) != 25:
                    raise ValueError(f"{study_id}/{metric} must contain 25 ROIs, found {len(patient_rows)}")
                for target in TARGETS:
                    matches = [r for r in patient_rows if r["ROI"] == target]
                    if len(matches) != 1:
                        raise ValueError(f"Expected one {target} row for {study_id}/{metric}")
                    row = matches[0]
                    performance.append({
                        "study_id": study_id, "category": target, "category_type": "Target",
                        "metric": metric, "value": float(row["Value"]), "unit": row["Unit"],
                        "n_components": 1, "source_sheet": "Prospective patient metrics",
                        "source_rows": str(row["_row"]), "source_cells": f"H{row['_row']}",
                        "derivation": "direct ROI value",
                    })
                for group in OAR_GROUPS:
                    matches = [r for r in patient_rows if r["ROI group"] == group]
                    expected_n = len(GROUP_ROIS[group])
                    if len(matches) != expected_n or {r["ROI"] for r in matches} != set(GROUP_ROIS[group]):
                        raise ValueError(f"Incomplete {group} for {study_id}/{metric}")
                    performance.append({
                        "study_id": study_id, "category": group, "category_type": "Grouped OARs",
                        "metric": metric, "value": float(np.mean([float(r["Value"]) for r in matches])),
                        "unit": matches[0]["Unit"], "n_components": expected_n,
                        "source_sheet": "Prospective patient metrics",
                        "source_rows": ";".join(str(r["_row"]) for r in matches),
                        "source_cells": ";".join(f"H{r['_row']}" for r in matches),
                        "derivation": f"within-patient mean across {expected_n} OARs",
                    })

        ws = wb["Prospective timing"]
        values = list(ws.values)
        header = list(values[0])
        timing = []
        for excel_row, cells in enumerate(values[1:], 2):
            if not any(v is not None for v in cells):
                continue
            record = dict(zip(header, cells))
            timing.append({
                "study_id": str(record["Study ID"]),
                "algorithm_s": float(record["Algorithm-processing time (s)"]),
                "review_s": float(record["Review-and-modification interval (s)"]),
                "combined_s": float(record["Derived combined time (s)"]),
                "source_sheet": "Prospective timing", "source_row": excel_row,
                "algorithm_cell": f"B{excel_row}", "review_cell": f"C{excel_row}",
                "combined_cell": f"D{excel_row}",
            })
    finally:
        wb.close()
    return performance, timing

def write_csv(rows, path, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def validate_data(performance, timing):
    study_ids = sorted({row["study_id"] for row in performance})
    expected = set(itertools.product(study_ids, CATEGORIES, METRICS))
    keys = [(r["study_id"], r["category"], r["metric"]) for r in performance]
    if len(study_ids) != 135 or len(keys) != 2430 or len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError("Performance data must contain 135 × 6 × 3 unique values")
    for row in performance:
        value = row["value"]
        if not math.isfinite(value) or value < 0 or (row["metric"] == "DSC" and value > 1):
            raise ValueError(f"Invalid performance value: {row}")
        expected_components = 1 if row["category"] in TARGETS else len(GROUP_ROIS[row["category"]])
        if row["n_components"] != expected_components:
            raise ValueError(f"Wrong component count: {row}")
    timing_ids = [row["study_id"] for row in timing]
    if len(timing) != 135 or len(timing_ids) != len(set(timing_ids)) or set(timing_ids) != set(study_ids):
        raise ValueError("Timing and performance cohorts must contain the same 135 patients")
    residual = 0.0
    for row in timing:
        values = [row[k] for k in ("algorithm_s", "review_s", "combined_s")]
        if not all(math.isfinite(v) and v >= 0 for v in values):
            raise ValueError(f"Invalid timing value: {row}")
        residual = max(residual, abs(row["combined_s"] - row["algorithm_s"] - row["review_s"]))
    if residual > 1e-9:
        raise ValueError(f"Combined time is not the exact sum; max residual={residual}")
    return {
        "patients": 135, "performance_values": len(performance), "timing_rows": len(timing),
        "missing": 0, "maximum_combined_time_residual_s": residual,
    }

def select_font(requested=None):
    candidates = ([requested] if requested else []) + ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
    for candidate in candidates:
        try:
            font_manager.findfont(candidate, fallback_to_default=False)
        except ValueError:
            continue
        mpl.rcParams["font.sans-serif"] = [candidate]
        if requested and candidate != requested:
            warnings.warn(f"{requested!r} unavailable; using {candidate!r}")
        return candidate
    raise RuntimeError("No supported sans-serif font available")

def summarise_performance(performance):
    summary = []
    for metric in METRICS:
        for category in CATEGORIES:
            values = np.array([
                row["value"] for row in performance
                if row["metric"] == metric and row["category"] == category
            ], dtype=float)
            p05, q1, median, q3, p95 = np.percentile(values, [5, 25, 50, 75, 95])
            summary.append({
                "metric": metric, "category": category, "n": len(values),
                "mean": float(np.mean(values)), "sd": float(np.std(values, ddof=1)),
                "p05": float(p05), "q1": float(q1), "median": float(median),
                "q3": float(q3), "p95": float(p95), "minimum": float(np.min(values)),
                "maximum": float(np.max(values)),
            })
    return summary

def draw_metric_axis(ax, metric, summary, show_labels):
    by_category = {row["category"]: row for row in summary if row["metric"] == metric}
    for category in CATEGORIES:
        row = by_category[category]
        y = ROW_Y[category]
        colour = TARGET_BLUE if category in TARGETS else OAR_GREEN
        ax.hlines(y, row["p05"], row["p95"], color=colour, linewidth=0.75, zorder=2)
        ax.plot([row["q1"], row["q3"]], [y, y], color=colour, linewidth=5.0,
                solid_capstyle="round", zorder=3)
        ax.plot(row["median"], y, marker="o", markersize=3.7, markerfacecolor="white",
                markeredgecolor=colour, markeredgewidth=0.8, linestyle="none", zorder=4)
    ax.axhline(1.75, color="#AEB5BA", linewidth=0.55)
    ax.set_xlim(*AXIS_LIMITS[metric])
    ax.set_xticks(AXIS_TICKS[metric])
    # Keep the adjacent DSC and HD95 endpoint labels visually separate.
    if metric == "DSC":
        ax.get_xticklabels()[-1].set_horizontalalignment("right")
    elif metric == "HD95":
        ax.get_xticklabels()[0].set_horizontalalignment("left")
    ax.set_ylim(-0.6, 6.15)
    ax.set_yticks([ROW_Y[c] for c in CATEGORIES])
    if show_labels:
        ax.set_yticklabels(CATEGORIES, fontsize=6.4)
        ax.tick_params(axis="y", length=0, pad=3)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
    ax.set_title(METRIC_TITLES[metric], fontsize=8.0, fontweight="bold", pad=6)
    ax.tick_params(axis="x", labelsize=6.0, length=2.5, pad=2)
    ax.grid(axis="x", color="#D9DEE2", linewidth=0.45)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#6F767C")
    ax.spines["bottom"].set_linewidth(0.6)

def draw_figure(performance, timing, output_dir):
    summary = summarise_performance(performance)
    width_mm, height_mm = 183.0, 175.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")

    fig.text(4 / width_mm, 169 / height_mm, "a", fontsize=9, fontweight="bold", ha="left", va="center")
    fig.text(10 / width_mm, 169 / height_mm, "Prospective contour concordance",
             fontsize=8.5, fontweight="bold", ha="left", va="center")
    fig.text(179 / width_mm, 169 / height_mm, "n = 135", fontsize=6.2,
             ha="right", va="center", color="#555555")

    legend_handles = [
        Line2D([], [], marker="o", markersize=4.4, linestyle="none", markerfacecolor=TARGET_BLUE,
               markeredgecolor=TARGET_BLUE, label="Target"),
        Line2D([], [], marker="s", markersize=4.4, linestyle="none", markerfacecolor=OAR_GREEN,
               markeredgecolor=OAR_GREEN, label="Grouped OARs"),
        Line2D([], [], color=BLACK, linewidth=0.75, label="5th–95th percentiles"),
        Line2D([], [], color=BLACK, linewidth=4.0, solid_capstyle="round", label="IQR"),
        Line2D([], [], marker="o", markersize=4.0, linestyle="none", markerfacecolor="white",
               markeredgecolor=BLACK, label="Median"),
    ]
    fig.legend(handles=legend_handles, loc="center", bbox_to_anchor=(0.53, 160.8 / height_mm),
               ncol=5, frameon=False, handlelength=1.8, columnspacing=1.5, fontsize=6.0)

    axes_specs = ((35.0, 48.0, "DSC"), (85.0, 42.0, "HD95"), (132.0, 43.0, "ASD"))
    for index, (left, axis_width, metric) in enumerate(axes_specs):
        ax = fig.add_axes([left / width_mm, 96 / height_mm, axis_width / width_mm, 56 / height_mm])
        draw_metric_axis(ax, metric, summary, show_labels=index == 0)

    fig.text(4 / width_mm, 82 / height_mm, "b", fontsize=9, fontweight="bold", ha="left", va="center")
    fig.text(10 / width_mm, 82 / height_mm, "Recorded intervals and their derived sum",
             fontsize=8.5, fontweight="bold", ha="left", va="center")

    timing_sorted = sorted(timing, key=lambda row: (row["combined_s"], row["study_id"]))
    algorithm = np.array([row["algorithm_s"] for row in timing_sorted], dtype=float)
    review = np.array([row["review_s"] for row in timing_sorted], dtype=float)
    combined = np.array([row["combined_s"] for row in timing_sorted], dtype=float)
    algorithm_mean, algorithm_sd = float(np.mean(algorithm)), float(np.std(algorithm, ddof=1))
    review_mean, review_sd = float(np.mean(review)), float(np.std(review, ddof=1))
    combined_mean, combined_sd = float(np.mean(combined)), float(np.std(combined, ddof=1))

    fig.text(70 / width_mm, 77.0 / height_mm,
             f"Algorithm {algorithm_mean:.0f} ± {algorithm_sd:.0f} s", fontsize=5.8,
             fontweight="bold", color=ALGORITHM_TEAL, ha="center", va="center")
    fig.text(116 / width_mm, 77.0 / height_mm,
             f"Review/modification {review_mean:.0f} ± {review_sd:.0f} s", fontsize=5.8,
             fontweight="bold", color="#C87516", ha="center", va="center")
    fig.text(163 / width_mm, 77.0 / height_mm,
             f"Derived sum {combined_mean:.0f} ± {combined_sd:.0f} s", fontsize=5.8,
             fontweight="bold", color=BLACK, ha="center", va="center")

    # Use the same 35–175 mm plot boundary as panel a. A stacked rank ribbon
    # retains every patient-level composition without the barcode appearance of
    # 135 separate bars.
    ax = fig.add_axes([35 / width_mm, 18 / height_mm, 140 / width_mm, 51 / height_mm])
    # Show the actual patient rank so the horizontal extent directly reflects
    # the prospective cohort size (n = 135), rather than a 0–100 percentile
    # transformation that can be mistaken for the number of patients.
    patient_rank = np.arange(1, len(timing_sorted) + 1)
    ax.fill_between(
        patient_rank, 0, algorithm, step="mid", color=ALGORITHM_TEAL,
        linewidth=0, alpha=0.96, zorder=2,
    )
    ax.fill_between(
        patient_rank, algorithm, combined, step="mid", color=REVIEW_ORANGE,
        linewidth=0, alpha=0.96, zorder=2,
    )
    ax.plot(patient_rank, algorithm, color="#287A88", linewidth=0.65,
            drawstyle="steps-mid", zorder=3)
    ax.plot(patient_rank, combined, color=BLACK, linewidth=0.8,
            drawstyle="steps-mid", zorder=4)
    for percentile in (5, 50, 95):
        value = float(np.percentile(combined, percentile))
        rank_position = 1 + (len(timing_sorted) - 1) * percentile / 100
        ax.plot(rank_position, value, marker="o", markersize=3.0, markerfacecolor="white",
                markeredgecolor=BLACK, markeredgewidth=0.7, zorder=5)
        ax.annotate(
            f"P{percentile}  {value:.0f} s", xy=(rank_position, value),
            xytext=(0, 5 if percentile != 95 else -10), textcoords="offset points",
            ha="center", va="bottom" if percentile != 95 else "top",
            fontsize=5.2, color=BLACK,
        )
    ax.set_xlim(1, len(timing_sorted))
    ax.set_ylim(0, 420)
    ax.set_xticks([1, 34, 68, 101, 135])
    ax.set_yticks([0, 100, 200, 300, 400])
    ax.set_xlabel("Patient rank by derived sum (n = 135)", fontsize=6.5, labelpad=3)
    ax.set_ylabel("Recorded interval / derived sum (s)", fontsize=6.5, labelpad=4)
    ax.tick_params(axis="both", labelsize=5.8, length=2.5)
    ax.grid(axis="y", color="#D9DEE2", linewidth=0.45)
    ax.set_axisbelow(True)
    ax.spines["left"].set_color("#6F767C")
    ax.spines["bottom"].set_color("#6F767C")
    ax.legend(
        handles=[
            Patch(facecolor=ALGORITHM_TEAL, edgecolor="none", label="Algorithm processing"),
            Patch(facecolor=REVIEW_ORANGE, edgecolor="none", label="Review and modification"),
            Line2D([], [], color=BLACK, linewidth=0.8, label="Derived sum"),
        ],
        loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=3, frameon=False,
        fontsize=5.5, handlelength=1.5, columnspacing=1.5,
    )

    fig.text(91.5 / width_mm, 6.5 / height_mm,
             "Derived sum = algorithm processing + review/modification; the recorded intervals overlap.",
             fontsize=5.3, ha="center", va="center", color="#555555")
    fig.text(91.5 / width_mm, 3.2 / height_mm,
             "The sum is not end-to-end elapsed time. All 135 patients are shown in order of their derived sum.",
             fontsize=5.3, ha="center", va="center", color="#555555")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(summary, output_dir / "performance_summary.csv", tuple(summary[0].keys()))
    ranked = []
    for rank, row in enumerate(timing_sorted, 1):
        ranked.append({
            "rank": rank, "study_id": row["study_id"], "algorithm_s": row["algorithm_s"],
            "review_s": row["review_s"], "combined_s": row["combined_s"],
            "source_row": row["source_row"],
        })
    write_csv(ranked, output_dir / "timing_ranked.csv", tuple(ranked[0].keys()))
    save_cns_figure(fig, output_dir / FIGURE_NAME)
    plt.close(fig)
    return summary, ranked, {
        "algorithm_mean": algorithm_mean, "algorithm_sd": algorithm_sd,
        "review_mean": review_mean, "review_sd": review_sd,
        "combined_mean": combined_mean, "combined_sd": combined_sd,
        "algorithm_min": float(np.min(algorithm)), "algorithm_max": float(np.max(algorithm)),
        "review_min": float(np.min(review)), "review_max": float(np.max(review)),
        "combined_min": float(np.min(combined)), "combined_max": float(np.max(combined)),
    }
