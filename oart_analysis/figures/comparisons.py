"""Figure 5b/c and Figure S6: paired means and direction-aligned changes."""

import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
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
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

import argparse

import csv

import hashlib

import itertools

import json

import math

from pathlib import Path

import platform

import warnings

mpl.use("Agg")

import matplotlib.pyplot as plt

from matplotlib import font_manager

from matplotlib.lines import Line2D

from matplotlib.transforms import Bbox

CENTRES = ("Main centre", "Centre I", "Centre II", "Centre III", "Centre IV")

ROIS = ("Small OARs (≤ 4 cm³)", "Large OARs (> 4 cm³)", "GTVp", "CTV1", "CTV2", "GTVn")

ROI_LABELS = ("Small OARs", "Large OARs", "GTVp", "CTV1", "CTV2", "GTVn")

METRICS = ("DSC", "HD95", "ASD")

SAMPLE_SIZES = {"Centre I": 14, "Centre II": 21, "Centre III": 14, "Centre IV": 25}

COLOURS = (BLACK, CATEGORICAL[0], CATEGORICAL[2], CATEGORICAL[4], CATEGORICAL_EXTENDED[9])

Y = (5.5, 4.5, 3.0, 2.0, 1.0, 0.0)

LIMITS = {"DSC": (0.5, 1.0), "HD95": (0.0, 9.0), "ASD": (0.0, 5.0)}

TICKS = {"DSC": [0.5, 0.75, 1.0], "HD95": [0, 4, 8], "ASD": [0, 2, 4]}

SOURCE_COLUMNS = {
    "centre": "Centre", "roi": "ROI", "metric": "Metric", "n": "Paired n",
    "unit": "Unit", "before_model": "Before model", "after_model": "After model",
    "before_mean": "Before mean", "after_mean": "After mean",
    "before_sd": "Before SD", "after_sd": "After SD",
    "mean_difference": "Mean difference (after-before)",
    "ci_lower": "95% CI lower", "ci_upper": "95% CI upper",
    "p_value": "Two-sided paired t-test P", "q_value": "BH-adjusted q",
}

NUMERIC = ("before_mean", "after_mean", "before_sd", "after_sd", "mean_difference", "ci_lower", "ci_upper")

def load_summary(path):
    """Read all rows, without sorting or rounding any measurement values."""
    path = Path(path)
    if path.suffix.lower() == ".xlsx":
        from openpyxl import load_workbook
        wb = load_workbook(path, read_only=True, data_only=True)
        try:
            ws = wb["Table S4"]
            data = list(ws.values)
            head = list(data[0])
            absent = set(SOURCE_COLUMNS.values()) - set(head)
            if absent:
                raise ValueError(f"Missing Table S4 columns: {sorted(absent)}")
            rows = []
            for excel_row, cells in enumerate(data[1:], 2):
                if not any(c is not None for c in cells):
                    continue
                record = dict(zip(head, cells))
                row = {key: record[col] for key, col in SOURCE_COLUMNS.items()}
                row.update(source_sheet="Table S4", source_row=excel_row,
                           before_cell=f"H{excel_row}", after_cell=f"P{excel_row}")
                rows.append(row)
        finally:
            wb.close()
    elif path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
    else:
        raise ValueError("Source must be a .csv or .xlsx file.")
    for row in rows:
        for key in NUMERIC:
            row[key] = float(row[key])
        count = float(row["n"])
        if not count.is_integer():
            raise ValueError(f"Non-integer sample size: {row}")
        row["n"] = int(count)
        row["source_row"] = int(row["source_row"])
    return rows

def validate_summary(rows):
    """Reject omitted/duplicated comparisons and wrong model/sample definitions."""
    wanted = set(itertools.product(CENTRES, ROIS, METRICS))
    keys = [(r["centre"], r["roi"], r["metric"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate centre/ROI/metric combination.")
    if set(keys) != wanted:
        raise ValueError(f"Coverage mismatch: missing={wanted-set(keys)}, extra={set(keys)-wanted}")
    for r in rows:
        key = (r["centre"], r["roi"], r["metric"])
        if not all(math.isfinite(r[c]) for c in NUMERIC):
            raise ValueError(f"Non-finite numeric value: {key}")
        expected_n = (12 if r["roi"] in ("GTVp", "GTVn") else 29) if r["centre"] == "Main centre" else SAMPLE_SIZES[r["centre"]]
        if r["n"] != expected_n:
            raise ValueError(f"Unexpected n for {key}: {r['n']} != {expected_n}")
        models = ("Conventional U-Net comparator", "Locked source model") if r["centre"] == "Main centre" else ("Pre-adaptation model", "Locally adapted model")
        if (r["before_model"], r["after_model"]) != models:
            raise ValueError(f"Model labels do not match the comparison: {key}")
        expected_unit = "unitless" if r["metric"] == "DSC" else "mm"
        if r["unit"] != expected_unit:
            raise ValueError(f"Unexpected unit: {key}")
        for col in ("before_mean", "after_mean"):
            value = r[col]
            if value < 0 or (r["metric"] == "DSC" and value > 1):
                raise ValueError(f"Out-of-range metric mean: {key}, {col}")
            lo, hi = LIMITS[r["metric"]]
            if not lo <= value <= hi:
                raise ValueError(f"Mean outside the figure axis limits: {key}; update LIMITS.")
        if min(r["before_sd"], r["after_sd"]) < 0:
            raise ValueError(f"Negative SD: {key}")
        # Stored means and the paired difference were each rounded to 4 decimals.
        if abs(r["after_mean"]-r["before_mean"]-r["mean_difference"]) > 0.00015001:
            raise ValueError(f"Inconsistent stored mean difference: {key}")
    return {"comparisons": len(rows), "mean_points": 2*len(rows),
            "centres": 5, "endpoints": 6, "metrics": 3, "omitted_rows": 0,
            "unit_of_analysis": "patient", "summary_statistic": "mean across patients"}

def select_font(requested):
    for candidate in ([requested] if requested else []) + ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]:
        try:
            font_manager.findfont(candidate, fallback_to_default=False)
        except ValueError:
            continue
        mpl.rcParams["font.sans-serif"] = [candidate]
        if requested and candidate != requested:
            warnings.warn(f"{requested!r} unavailable; using {candidate!r}.")
        return candidate
    raise RuntimeError("No supported sans-serif font available.")

def facet_scale(rows, centre, metric):
    """Set a linear range around ALL 12 means within this centre/metric facet.

    Six percent of the observed span is added at each end before rounding
    outward to 0.02 DSC, 0.2 mm HD95 or 0.1 mm ASD. The same deterministic
    rule applies to every centre. Theoretical bounds are 0--1 for DSC and
    nonnegative for distance metrics. No individual endpoint gets its own scale.
    """
    values=[r[k] for r in rows if r['centre']==centre and r['metric']==metric
            for k in ('before_mean','after_mean')]
    if len(values)!=12:raise ValueError('Each facet must contain all 12 means.')
    data_min,data_max=min(values),max(values)
    step={'DSC':0.02,'HD95':0.2,'ASD':0.1}[metric]
    pad=max((data_max-data_min)*0.06,step*0.5)
    lo=round(max(0,math.floor((data_min-pad)/step)*step),8)
    hi=round(math.ceil((data_max+pad)/step)*step,8)
    if metric=='DSC':hi=min(1.0,hi)
    if not lo<=data_min<=data_max<=hi or hi<=lo:
        raise ValueError('Invalid facet range.')
    ticks=[lo,round((lo+hi)/2,8),hi]
    old_lo,old_hi=LIMITS[metric]
    return dict(centre=centre,metric=metric,data_min=data_min,data_max=data_max,
                axis_min=lo,axis_max=hi,ticks=ticks,scale='linear',
                width_gain_vs_original=(old_hi-old_lo)/(hi-lo))

from matplotlib.patches import Rectangle, FancyArrowPatch

LEFTS=(38.,70.,98.,126.,154.)

GAIN_LIMITS={'DSC':(-0.06,0.35),'HD95':(-2.8,5.3),'ASD':(-0.65,3.15)}

GAIN_TICKS={'DSC':[0,0.30],'HD95':[-2,0,4],'ASD':[0,3]}

ROW_GUIDE='#D5D9DC'

CI_GREY='#C7C9CC'

ZERO_GREY='#777D82'

def gain_values(r):
    if r['metric']=='DSC':return r['mean_difference'],r['ci_lower'],r['ci_upper']
    return -r['mean_difference'],-r['ci_upper'],-r['ci_lower']

def new_canvas(height):
    fig=plt.figure(figsize=(183/25.4,height/25.4),facecolor='white')
    def label(x,y,text,**kw):
        return fig.text(x/183,y/height,text,color=BLACK,**kw)
    return fig,label

def draw_centres(fig,label,height,group_y,name_y,n_y):
    label(49,group_y,'Internal testing',fontsize=7,ha='center',va='center',fontweight='bold')
    label(123,group_y,'External-centre adaptation',fontsize=7,ha='center',va='center',fontweight='bold')
    for j,c in enumerate(CENTRES):
        label(LEFTS[j]+11,name_y,c,fontsize=7.2,ha='center',va='center',fontweight='bold')
        label(LEFTS[j]+11,n_y,'n = 29 / 12*' if j==0 else f'n = {SAMPLE_SIZES[c]}',
              fontsize=6.5,ha='center',va='center')

def draw_row(fig,label,height,rows,metric,panel,bottom,plot_height,kind,logs,ranges):
    lookup={(r['centre'],r['roi'],r['metric']):r for r in rows}
    axes=[]
    for j,c in enumerate(CENTRES):
        ax=fig.add_axes([LEFTS[j]/183,bottom/height,22/183,plot_height/height])
        axes.append(ax)
        ax.set_ylim(-0.6,6.1)
        if kind=='mean':
            s=facet_scale(rows,c,metric);ranges.append(dict(panel=panel,kind=kind,**s))
            ax.set_xlim(s['axis_min'],s['axis_max']);ax.set_xticks(s['ticks'])
            ax.set_xticklabels([f'{t:.2f}' for t in s['ticks']])
        else:
            ax.set_xlim(*GAIN_LIMITS[metric]);ax.set_xticks(GAIN_TICKS[metric])
            ax.set_xticklabels([f'{t:.2f}' if metric=='DSC' else f'{t:g}' for t in GAIN_TICKS[metric]])
            ax.axvline(0,color=ZERO_GREY,lw=0.65,zorder=1)
        ax.set_yticks(Y);ax.set_yticklabels(ROI_LABELS if j==0 else ['']*6)
        ax.tick_params(axis='y',length=0,pad=6,labelsize=7)
        ax.tick_params(axis='x',length=2.2,pad=3,labelsize=6.5,colors=BLACK)
        ax.spines['left'].set_visible(False);ax.spines['bottom'].set_color('#7D858D')
        for guide_y in Y:
            ax.axhline(guide_y,color=ROW_GUIDE,lw=0.45,zorder=0)
        ax.axhline(3.75,color='#BFC5C9',lw=0.7,zorder=0)
        for roi,y in zip(ROIS,Y):
            r=lookup[(c,roi,metric)]
            common=dict(panel=panel,centre=c,roi=roi,metric=metric,n=r['n'],source_sheet=r['source_sheet'])
            if kind=='mean':
                x0,x1=r['before_mean'],r['after_mean'];marker='s' if j==0 else 'o'
                ax.plot([x0,x1],[y,y],color=CI_GREY,lw=1.7,solid_capstyle='round',zorder=2)
                ax.plot(x0,y,marker=marker,ms=4.6,mfc='white',mec=COLOURS[j],mew=0.9,ls='none',zorder=3)
                ax.plot(x1,y,marker=marker,ms=3.3,mfc=COLOURS[j],mec='white',mew=0.35,ls='none',zorder=4)
                for state,x in [('before',x0),('after',x1)]:
                    logs['means'].append(dict(**common,state=state,mean=x,source_cell=r[state+'_cell']))
            else:
                x,lo,hi=gain_values(r)
                if not GAIN_LIMITS[metric][0]<=lo<=x<=hi<=GAIN_LIMITS[metric][1]:
                    raise ValueError(f'Gain CI outside axis: {c}, {roi}, {metric}')
                marker='s' if j==0 else 'o'
                ax.plot([lo,hi],[y,y],color=CI_GREY,lw=3.0,solid_capstyle='round',zorder=2)
                ax.plot([lo,hi],[y,y],color=COLOURS[j],lw=0.65,solid_capstyle='round',alpha=0.75,zorder=3)
                ax.plot(x,y,marker=marker,ms=3.8,mfc=COLOURS[j],mec='white',mew=0.45,ls='none',zorder=4)
                logs['gains'].append(dict(**common,gain=x,ci_lower=lo,ci_upper=hi,
                                         definition=(
                                             'locked-source-minus-conventional-U-Net' if c=='Main centre' and metric=='DSC'
                                             else 'conventional-U-Net-minus-locked-source' if c=='Main centre'
                                             else 'locally-adapted-minus-pre-adaptation' if metric=='DSC'
                                             else 'pre-adaptation-minus-locally-adapted'),
                                         difference_cell=f"W{r['source_row']}",
                                         ci_source_cells=f"Y{r['source_row']}:Z{r['source_row']}"))
    return axes

def symbol(fig,x,y,height,shape,filled):
    fig.add_artist(Line2D([x/183],[y/height],transform=fig.transFigure,marker=shape,
                         ms=2.7 if filled else 4.2,mfc=BLACK if filled else 'white',
                         mec=BLACK,mew=0.7,ls='none'))

def model_legend(fig,label,height,y):
    label(2.5,y,'Main centre:',fontsize=6.2,va='center')
    symbol(fig,29,y,height,'s',False);label(32,y,'Conventional U-Net',fontsize=6.2,va='center')
    symbol(fig,70,y,height,'s',True);label(73,y,'Locked source model',fontsize=6.2,va='center')
    label(2.5,y-5,'External centres:',fontsize=6.2,va='center')
    symbol(fig,29,y-5,height,'o',False);label(32,y-5,'Pre-adaptation',fontsize=6.2,va='center')
    symbol(fig,70,y-5,height,'o',True);label(73,y-5,'Locally adapted',fontsize=6.2,va='center')

def check_layout(fig,axes):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();frame=fig.bbox
    for t in fig.findobj(mpl.text.Text):
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0<-1 or b.x1>frame.x1+1 or b.y0<-1 or b.y1>frame.y1+1:
            raise ValueError(f'Text outside canvas: {t.get_text()!r}')
    all_ticks=[t for ax in axes for t in ax.get_xticklabels() if t.get_visible()]
    for i,a in enumerate(all_ticks):
        for b in all_ticks[i+1:]:
            if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):
                raise ValueError('Overlapping tick labels.')

def compose_figure5(rows):
    height=158.;fig,label=new_canvas(height);logs=dict(means=[],gains=[]);ranges=[]
    label(2.5,151,'b',fontsize=10,fontweight='bold',va='center')
    label(10,151,'DSC',fontsize=8.5,fontweight='bold',va='center')
    draw_centres(fig,label,height,151,144,139)
    axes=draw_row(fig,label,height,rows,'DSC','5b',94,38,'mean',logs,ranges)
    fig.add_artist(Line2D([65.5/183]*2,[90/height,154/height],transform=fig.transFigure,color='#D5DADF',lw=0.6))
    model_legend(fig,label,height,84)
    label(2.5,70,'c',fontsize=10,fontweight='bold',va='center')
    label(10,70,'DSC improvement',fontsize=8.5,fontweight='bold',va='center')
    label(178,70,'Mean paired difference (95% CI); positive favours locked/adapted →',fontsize=6.5,ha='right',va='center')
    for j,c in enumerate(CENTRES):
        label(LEFTS[j]+11,63,c,fontsize=7,ha='center',va='center',fontweight='bold')
    axes+=draw_row(fig,label,height,rows,'DSC','5c',21,37,'gain',logs,ranges)
    label(2.5,12,'b: centre-specific DSC ranges. c: common scale; main = locked − conventional, external = adapted − pre-adaptation.',fontsize=6,va='center')
    label(2.5,7.5,'*Main centre: OARs/CTVs n = 29; GTVs n = 12. Small/large OARs: ≤4/>4 cm³. Points in b: patient means.',fontsize=5.8,va='center')
    label(2.5,3,'CIs in c are pointwise, not multiplicity-adjusted. Distance metrics: Fig. S6. Full statistics: Table S4.',fontsize=6,va='center')
    check_layout(fig,axes)
    assert len(logs['means'])==60 and len(logs['gains'])==30
    return fig,logs,dict(width_mm=183,height_mm=height,quantitative_axes=10,axis_ranges=ranges,
                         gain_axis=GAIN_LIMITS['DSC'],gain_centres=list(CENTRES),
                         mean_points=60,gain_points=30,clipped_values=0)

def compose_figureS6(rows):
    height=226.;fig,label=new_canvas(height);logs=dict(means=[],gains=[]);ranges=[];axes=[]
    draw_centres(fig,label,height,222,214,209)
    specs=[('HD95','a',166,'mean'),('HD95','b',119,'gain'),
           ('ASD','c',72,'mean'),('ASD','d',25,'gain')]
    for metric,p,bottom,kind in specs:
        label(2.5,bottom+36,p,fontsize=10,fontweight='bold',va='center')
        title=f'{metric} (mm)' if kind=='mean' else f'{metric} reduction (mm)'
        label(10,bottom+36,title,fontsize=8.2,fontweight='bold',va='center')
        label(178,bottom+36,'← Lower is better' if kind=='mean' else 'Paired change (95% CI); improvement →',
              fontsize=6.3,ha='right',va='center')
        axes+=draw_row(fig,label,height,rows,metric,'S6'+p,bottom,30,kind,logs,ranges)
    model_legend(fig,label,height,15)
    label(2.5,5,'Absolute values: centre-specific ranges. Reductions: main conventional − locked; external pre − adapted. CIs: pointwise 95%.',fontsize=5.7,va='center')
    label(2.5,1.3,'*Main centre n = 29 for OARs/CTVs; n = 12 for GTVs. Small/large OARs: ≤4/>4 cm³. Full statistics: Table S4.',fontsize=5.6,va='center')
    check_layout(fig,axes)
    assert len(logs['means'])==120 and len(logs['gains'])==60
    return fig,logs,dict(width_mm=183,height_mm=height,quantitative_axes=20,axis_ranges=ranges,
                         mean_points=120,gain_points=60,clipped_values=0)
