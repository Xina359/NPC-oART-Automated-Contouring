"""Figure S5: within-case OAR differences for the selected representative case."""

from __future__ import annotations

import argparse

import csv

import json

from pathlib import Path

from statistics import mean

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

MM = 1 / 25.4

WIDTH_MM = 183

HEIGHT_MM = 205

INK = "#263746"

MUTED = "#71808A"

BASE = "#7A858D"

MODEL = "#176B8C"

MODEL_LIGHT = "#D1E4EC"

ACCENT = "#2D8D83"

GRID = "#DDE6EA"

GROUP_A = "#F5F8F9"

GROUP_B = "#FCFDFD"

WHITE = "#FFFFFF"

DISPLAY_NAMES = {
    "Eustachian tube bone": "ET bone",
    "Temporomandibular joint": "TM joint",
    "Vestibule/semicircular canals": "Vestibular canals",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6,
    "axes.linewidth": 0.55,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})

def layout_rows(rows):
    positions, headers, ranges = [], {}, {}
    y, current, start = 0.0, None, None
    for row in rows:
        group = row["group"]
        if group != current:
            if current is not None:
                ranges[current] = (start, positions[-1])
                y += .70
            current, start = group, y
            headers[group] = y - .54
        positions.append(y)
        y += 1.0
    ranges[current] = (start, positions[-1])
    return positions, headers, ranges, y - 1

def add_group_structure(ax, ranges, groups):
    for i, group in enumerate(groups):
        lo, hi = ranges[group]
        ax.axhspan(lo - .48, hi + .48,
                   color=GROUP_A if i % 2 == 0 else GROUP_B,
                   linewidth=0, zorder=-20)
        if i:
            ax.axhline(lo - .71, color=GRID, lw=.65, zorder=-3)

def style_change_axis(ax, xlim, ticks, title, subtitle, ymax):
    ax.set_xlim(*xlim)
    ax.set_ylim(ymax + .60, -1.08)
    ax.set_yticks([])
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", colors=MUTED, labelsize=5.3,
                   length=2.5, width=.55, pad=2)
    for side in ("top", "left", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#98A5AD")
    ax.grid(axis="x", color=GRID, lw=.50, zorder=-10)
    ax.axvline(0, color="#A9B5BC", lw=.75, zorder=-2)
    ax.set_xlabel(title, fontsize=6.5, fontweight="bold", color=INK, labelpad=3)
    ax.set_title(subtitle, fontsize=4.9, color=MUTED, pad=4)

def draw_baseline_centred_pair(ax, y, change, baseline, model,
                               absolute_decimals, change_decimals, xmax):
    # The pale under-stroke creates a soft, rounded pairing ribbon.
    ax.plot([0, change], [y, y], color=MODEL_LIGHT, lw=4.5,
            solid_capstyle="round", zorder=2)
    ax.plot([0, change], [y, y], color=MODEL, lw=1.15,
            alpha=.9, solid_capstyle="round", zorder=3)
    ax.scatter([0], [y], s=29, marker="o", facecolor=WHITE,
               edgecolor=BASE, linewidth=1.05, zorder=6)
    ax.scatter([change], [y], s=35, marker="D", facecolor=MODEL,
               edgecolor=WHITE, linewidth=.75, zorder=7)

    abs_fmt = f"{{:.{absolute_decimals}f}}"
    ax.annotate(abs_fmt.format(baseline), (0, y), xytext=(0, 4.2),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=4.35, color=BASE, zorder=8)
    ax.annotate(abs_fmt.format(model), (change, y), xytext=(0, -4.2),
                textcoords="offset points", ha="center", va="top",
                fontsize=4.45, color=MODEL, fontweight="bold", zorder=8)

    delta = f"{change:+.{change_decimals}f}"
    ax.annotate(delta, (change, y),
                xytext=(6, 0), textcoords="offset points",
                ha="left", va="center",
                fontsize=4.15, color=ACCENT, fontweight="bold", zorder=8)

def draw_header(fig, rows):
    ax = fig.add_axes([.025, .908, .95, .075])
    ax.set_axis_off()
    ax.text(0, .86, "MC_R023 · representative OAR case",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=8.6, fontweight="bold", color=INK)
    ax.text(0, .47, f"{len(rows)} OARs · {sum(r[f'beneficial_change_{m}'] > 0 for r in rows for m in ('dsc', 'hd95', 'asd'))} of {len(rows)*3} changes favour VAG-UNet",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=5.55, color=MUTED)

    summaries = [
        ("DSC", mean(r["baseline_dsc"] for r in rows),
         mean(r["vag_dsc"] for r in rows), 3, f"{mean(r['vag_dsc']-r['baseline_dsc'] for r in rows):+.3f}"),
        ("HD95", mean(r["baseline_hd95"] for r in rows),
         mean(r["vag_hd95"] for r in rows), 2, f"{mean(r['vag_hd95']-r['baseline_hd95'] for r in rows):+.2f} mm"),
        ("ASD", mean(r["baseline_asd"] for r in rows),
         mean(r["vag_asd"] for r in rows), 2, f"{mean(r['vag_asd']-r['baseline_asd'] for r in rows):+.2f} mm"),
    ]
    x0, width = .47, .172
    for i, (metric, before, after, dec, delta) in enumerate(summaries):
        x = x0 + i * width
        if i:
            ax.plot([x - .018, x - .018], [.18, .82], color=GRID, lw=.8,
                    transform=ax.transAxes, clip_on=False)
        ax.text(x, .78, metric, transform=ax.transAxes,
                fontsize=5.3, color=MUTED, fontweight="bold", va="top")
        ax.text(x, .48, f"{before:.{dec}f} → {after:.{dec}f}",
                transform=ax.transAxes, fontsize=6.0, color=INK,
                fontweight="bold", va="top")
        ax.text(x, .18, delta, transform=ax.transAxes,
                fontsize=5.0, color=ACCENT, fontweight="bold", va="top")

    ax.plot([0, 1], [.01, .01], color=GRID, lw=.8,
            transform=ax.transAxes, clip_on=False)

def render(rows, output_dir: Path):
    ypos, headers, ranges, ymax = layout_rows(rows)
    groups = list(dict.fromkeys(row["group"] for row in rows))

    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), facecolor=WHITE)
    draw_header(fig, rows)
    gs = fig.add_gridspec(
        1, 4, left=.025, right=.982, top=.874, bottom=.065,
        width_ratios=[1.50, 1.18, 1.18, 1.18], wspace=.20,
    )
    ax_label = fig.add_subplot(gs[0, 0])
    axes = [fig.add_subplot(gs[0, i]) for i in range(1, 4)]
    for ax in [ax_label, *axes]:
        add_group_structure(ax, ranges, groups)

    ax_label.set_xlim(0, 1)
    ax_label.set_ylim(ymax + .60, -1.08)
    ax_label.set_axis_off()
    for row, y in zip(rows, ypos):
        ax_label.text(.055, y, DISPLAY_NAMES.get(row["roi"], row["roi"]),
                      ha="left", va="center", fontsize=5.65, color=INK)
    for group in groups:
        ax_label.text(0, headers[group], group, ha="left", va="center",
                      fontsize=5.95, color=MODEL, fontweight="bold")

    configs = [
        ("dsc", (-.006, .096), [0, .02, .04, .06, .08],
         "DSC gain", "baseline centred at 0 · absolute values label points",
         3, 3, .096),
        ("hd95", (-.08, 1.31), [0, .3, .6, .9, 1.2],
         "HD95 reduction (mm)", "rightward distance = improvement",
         2, 2, 1.31),
        ("asd", (-.012, .205), [0, .04, .08, .12, .16, .20],
         "ASD reduction (mm)", "rightward distance = improvement",
         2, 2, .205),
    ]
    for cfg, ax in zip(configs, axes):
        metric, xlim, ticks, title, subtitle, abs_dec, change_dec, xmax = cfg
        observed = [r[f"beneficial_change_{metric}"] for r in rows]
        if min(observed) < xlim[0] or max(observed) > xlim[1] * .88:
            raise ValueError(f"Selected-case values exceed the configured {metric} display range; adjust the range.")
        style_change_axis(ax, xlim, ticks, title, subtitle, ymax)
        for row, y in zip(rows, ypos):
            draw_baseline_centred_pair(
                ax, y, row[f"beneficial_change_{metric}"],
                row[f"baseline_{metric}"], row[f"vag_{metric}"],
                abs_dec, change_dec, xmax,
            )

    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=WHITE,
               markeredgecolor=BASE, markeredgewidth=1.0,
               markersize=5.2, label="Baseline (zero change)"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=MODEL,
               markeredgecolor=WHITE, markeredgewidth=.7,
               markersize=5.4, label="VAG-UNet"),
    ]
    fig.legend(handles=legend, loc="lower left", bbox_to_anchor=(.025, .886),
               ncol=2, frameon=False, fontsize=5.25,
               handletextpad=.45, columnspacing=1.25)
    fig.text(.025, .022,
             "Horizontal spans encode within-OAR gain or reduction; numeric point labels are the absolute metric values.",
             fontsize=4.75, color=MUTED, ha="left", va="bottom")

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "FigureS5_MC_R023_baseline_centered"
    metadata = {"Title": "MC_R023 baseline-centred paired improvements"}
    fig.savefig(stem.with_suffix(".svg"), metadata=metadata)
    fig.savefig(stem.with_suffix(".pdf"), dpi=300,
                metadata={**metadata, "Creator": "Python matplotlib",
                          "CreationDate": None, "ModDate": None})
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})

    svg = stem.with_suffix(".svg").read_text(encoding="utf-8")
    report = {
        "study_id": "MC_R023",
        "n_oars": len(rows),
        "n_improving_changes": sum(
            row[f"beneficial_change_{m}"] > 0
            for row in rows for m in ("dsc", "hd95", "asd")
        ),
        "coordinate_system": "baseline-centred within-OAR improvement",
        "absolute_values_retained_as_labels": True,
        "width_mm": WIDTH_MM,
        "height_mm": HEIGHT_MM,
        "svg_text_nodes": svg.count("<text"),
        "svg_embeds_raster": "<image" in svg,
    }
    (output_dir / "render_report_baseline_centered.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    if report["svg_embeds_raster"]:
        raise RuntimeError("SVG unexpectedly contains an embedded raster image")
    plt.close(fig)
