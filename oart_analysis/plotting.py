"""Render the same figures from normalized bundled or workbook plotting records."""
from pathlib import Path
from .io import write_csv, write_json

FIGURES = ("5", "6", "S4", "S5", "S6", "S7")


def save_figure(fig, stem):
    import matplotlib.pyplot as plt
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "pdf", "png", "tiff"):
        options = {"dpi": 600 if ext == "tiff" else 300, "bbox_inches": None}
        if ext == "tiff": options["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(stem.with_suffix('.'+ext), **options)
    plt.close(fig)


def summary_rows(frame):
    names = ("centre", "roi_group", "roi", "metric", "unit", "n", "before_model", "before_mean", "before_sd", "before_min",
             "before_max", "before_q1", "before_median", "before_q3", "after_model", "after_mean", "after_sd", "after_min",
             "after_max", "after_q1", "after_median", "after_q3", "mean_difference", "sd_paired_difference", "ci_lower", "ci_upper", "p_value", "q_value")
    columns = ("Centre", "ROI group", "ROI", "Metric", "Unit", "Paired n", "Before model", "Before mean", "Before SD", "Before minimum",
               "Before maximum", "Before Q1", "Before median", "Before Q3", "After model", "After mean", "After SD", "After minimum",
               "After maximum", "After Q1", "After median", "After Q3", "Mean difference (after-before)", "SD of paired differences",
               "95% CI lower", "95% CI upper", "Two-sided paired t-test P", "BH-adjusted q")
    rows=[]
    for raw in frame.to_dict('records'):
        row=dict(zip(names, [raw[c] for c in columns]))
        n=int(raw['source_row'])
        row.update(source_sheet='Table S5',source_row=n,difference_cell=f'W{n}',q_cell=f'AB{n}')
        row['p_numeric']=float(str(row['p_value']).lstrip('<'))
        row['q_numeric']=float(str(row['q_value']).lstrip('<'))
        rows.append(row)
    return rows


def render(data, figure, output):
    """Render one figure from validated plotting records."""
    output = Path(output)/f"Figure{figure}"
    output.mkdir(parents=True, exist_ok=True)
    if figure in ("5", "S6"):
        from .figures import comparisons as plot
        plot.select_font(None)
        rows=data['summary_s4']
        plot.validate_summary(rows)
        composer=plot.compose_figure5 if figure=='5' else plot.compose_figureS6
        fig,logs,layout=composer(rows)
        save_figure(fig,output/("Figure5_bc" if figure=='5' else 'FigureS6_abcd'))
        for kind,records in logs.items(): write_csv(records,output/f'plotted_{kind}.csv')
        write_json(layout,output/'layout.json')
    elif figure=='6':
        from .figures import prospective as plot
        plot.select_font(None)
        performance,timing=data['performance'],data['timing']
        plot.validate_data(performance,timing)
        summary,ranked,stats=plot.draw_figure(performance,timing,output)
        write_csv(performance,output/'plotted_patient_values.csv')
        write_json(stats,output/'timing_summary.json')
    elif figure=='S4':
        from .figures import gtvn_strategies as plot
        frame=data['gtvn']
        values={(r['Study ID'],r['Strategy'],r['Metric']):r['Value'] for r in frame}
        patients=sorted({r['Study ID'] for r in frame})
        for (patient,strategy,metric),value in values.items():
            lo,hi=plot.METRIC_SPECS[metric]['ylim']
            if not lo<=value<=hi: raise ValueError(f'GTVn {metric} outside plot axis')
        fig=plot.create_lower_figure(patients,values)
        save_figure(fig,output/'FigureS4_bcd')
        write_csv(frame,output/'plotted_patient_values.csv')
    elif figure=='S5':
        from .figures import representative_case as plot
        rows=data['case']
        plot.render(rows,output)
        write_csv(rows,output/'plotted_patient_values.csv')
    elif figure=='S7':
        from .figures import organ_effects as plot
        rows=data['summary_s5']
        plot.validate_data(rows)
        plot.select_font(None)
        axes=plot.draw_figure(rows,output)
        write_json(axes,output/'axes.json')
    else:
        raise ValueError(f'Unknown figure {figure}. Choose {FIGURES}')
    return output
