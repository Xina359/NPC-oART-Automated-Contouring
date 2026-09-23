"""Recompute the study tables from patient records and validate archived summaries."""
from pathlib import Path
import math
import numpy as np
import pandas as pd
from .constants import (CENTRES, METRICS, TARGETS, SMALL, LARGE, GROUPS, GROUP_LABELS,
                        ENDPOINTS, STRATEGIES, TIMES, CASE_ORDER, models, expected_n)
from .io import require_unique, require_numeric, write_csv, write_json
from .statistics import describe, paired_test, benjamini_hochberg, orient_effect, precision_sample_size

KEY = ["Centre", "Model state", "Study ID", "ROI", "Metric"]
COMPOSITE_KEY = ["Centre", "Model state", "Study ID", "ROI group", "Metric"]


def validate_metrics(frame, prospective=False):
    context = "prospective" if prospective else "retrospective"
    require_unique(frame, KEY, context)
    require_numeric(frame, ["Value"], context)
    if set(frame["Metric"]) != set(METRICS) or set(frame["ROI"]) != set((*SMALL, *LARGE, *TARGETS)):
        raise ValueError(f"Unexpected metric/ROI definitions: {context}")
    for row in frame.to_dict("records"):
        group = next((g for g, rois in GROUPS.items() if row["ROI"] in rois), "Targets")
        unit = "unitless" if row["Metric"] == "DSC" else "mm"
        if row["ROI group"] != group or row["Unit"] != unit:
            raise ValueError(f"ROI group or unit mismatch at source row {row['source_row']}")
        if row["Value"] < 0 or (row["Metric"] == "DSC" and row["Value"] > 1):
            raise ValueError(f"Metric outside logical range at source row {row['source_row']}")
    roi_codes = frame[["ROI", "ROI code"]].drop_duplicates()
    if roi_codes["ROI"].duplicated().any() or roi_codes["ROI code"].duplicated().any() or roi_codes.isna().any().any():
        raise ValueError("ROI labels and codes are not a one-to-one mapping")
    if prospective:
        ids = set(frame["Study ID"])
        keys = set(map(tuple, frame[["Study ID", "ROI", "Metric"]].to_numpy()))
        expected = {(p, roi, m) for p in ids for roi in (*SMALL, *LARGE, *TARGETS) for m in METRICS}
        if len(ids) != 135 or len(frame) != len(expected) or keys != expected:
            raise ValueError("Prospective cohort must be complete: 135 patients × 25 ROIs × 3 metrics")
        if set(frame['Centre']) != {'Main centre prospective cohort'} or set(frame['Model state']) != {'Final clinician-approved contours versus automatic contours'}:
            raise ValueError('Unexpected prospective centre/model labels')
    else:
        if set(frame["Centre"]) != set(CENTRES):
            raise ValueError("Expected all five centres")
        for centre in CENTRES:
            c = frame[frame["Centre"] == centre]
            if set(c["Model state"]) != set(models(centre)):
                raise ValueError(f"Unexpected model states: {centre}")
            cohorts = {}
            for roi in (*SMALL, *LARGE, *TARGETS):
                for model in models(centre):
                    for metric in METRICS:
                        g = c[(c.ROI == roi) & (c["Model state"] == model) & (c.Metric == metric)]
                        ids = set(g["Study ID"])
                        if len(g) != expected_n(centre, roi):
                            raise ValueError(f"Incorrect count: {centre}/{model}/{roi}/{metric}")
                        category = "GTV" if centre == "Main centre" and roi in ("GTVp", "GTVn") else "CTV/OAR"
                        if category in cohorts and ids != cohorts[category]:
                            raise ValueError(f"Patient IDs differ across paired conditions: {centre}/{roi}")
                        cohorts[category] = ids
            if centre == "Main centre" and cohorts["GTV"] & cohorts["CTV/OAR"]:
                raise ValueError("The two main-centre test cohorts must be distinct")


def compute_composites(raw, round_digits=None):
    """Equal weight per ROI within each patient, requiring all 10 or all 11 ROIs."""
    rows = []
    oars = raw[raw["ROI group"].isin(GROUPS)]
    require_unique(oars, KEY, "OAR components")
    for key, frame in oars.groupby(COMPOSITE_KEY, sort=True):
        expected = GROUPS[key[3]]
        if len(frame) != len(expected) or set(frame.ROI) != set(expected):
            raise ValueError(f"Incomplete OAR composite: {key}")
        value = math.fsum(frame.Value)/len(frame)
        if round_digits is not None:
            value = round(value, round_digits)
        row = dict(zip(COMPOSITE_KEY, key))
        row.update({"ROI": GROUP_LABELS[row["ROI group"]], "ROI group": GROUP_LABELS[row["ROI group"]],
                    "Value": value, "Unit": "unitless" if row["Metric"] == "DSC" else "mm",
                    "Number of OARs averaged": len(frame),
                    "source_rows": ";".join(map(str, frame.source_row))})
        rows.append(row)
    return pd.DataFrame(rows)


def check_stored_composites(raw, stored):
    require_unique(stored, COMPOSITE_KEY, "stored OAR means")
    require_numeric(stored, ["Within-patient OAR mean", "Number of OARs averaged"], "stored OAR means")
    calculated = compute_composites(raw)
    joined = calculated.merge(stored, on=COMPOSITE_KEY, how="outer", validate="one_to_one", indicator=True, suffixes=("", "_stored"))
    if not joined._merge.eq("both").all():
        raise ValueError("Stored composites are missing or have unmatched patients")
    joined["absolute_difference"] = abs(joined.Value-joined["Within-patient OAR mean"])
    if (joined.absolute_difference > 0.00005000001).any() or not joined["Number of OARs averaged"].eq(joined["Number of OARs averaged_stored"]).all() or not joined.Unit.eq(joined.Unit_stored).all():
        raise ValueError("Stored OAR means disagree with constituent ROIs or their unit/count")
    return joined.drop(columns="_merge")


def paired_values(frame, centre, roi, metric):
    """Exact outer-join check, then deterministic patient ordering."""
    group = frame[(frame.Centre == centre) & (frame.ROI == roi) & (frame.Metric == metric)]
    require_unique(group, KEY, "paired input")
    before, after = models(centre)
    x = group[group["Model state"] == before].set_index("Study ID").Value
    y = group[group["Model state"] == after].set_index("Study ID").Value
    if set(x.index) != set(y.index) or len(x) != expected_n(centre, roi):
        raise ValueError(f"Unmatched patient IDs or incorrect paired n: {centre}/{roi}/{metric}")
    ids = sorted(x.index)
    return ids, x.loc[ids].to_numpy(float), y.loc[ids].to_numpy(float)


def calculate_table(frame, endpoints):
    result = []
    for centre in CENTRES:
        family = []
        for roi in endpoints:
            for metric in METRICS:
                ids, x, y = paired_values(frame, centre, roi, metric)
                values = paired_test(x, y)
                before, after = models(centre)
                group = "Small OARs" if roi in SMALL else "Large OARs" if roi in LARGE else roi if roi in GROUP_LABELS.values() else "Targets"
                row = {"Centre": centre, "ROI group": group, "ROI": roi, "Metric": metric,
                       "Unit": "unitless" if metric == "DSC" else "mm", "Paired n": len(ids),
                       "Before model": before}
                for state in ("Before", "After"):
                    if state == "After": row["After model"] = after
                    for label, name in (("mean", "mean"), ("SD", "sd"), ("minimum", "minimum"), ("maximum", "maximum"),
                                        ("Q1", "q1"), ("median", "median"), ("Q3", "q3")):
                        row[f"{state} {label}"] = values[state.lower()][name]
                d = values["difference"]
                row.update({"Mean difference (after-before)": d["mean"], "SD of paired differences": d["sd"],
                            "95% CI lower": d["ci_lower"], "95% CI upper": d["ci_upper"],
                            "Two-sided paired t-test P": values["p"]})
                family.append(row)
        for row, q in zip(family, benjamini_hochberg([r["Two-sided paired t-test P"] for r in family])):
            row["BH-adjusted q"] = float(q)
        result.extend(family)
    return pd.DataFrame(result)


def compare_summary(calculated, stored, sheet):
    keys = ["Centre", "ROI", "Metric"]
    require_unique(stored, keys, sheet)
    got = stored.set_index(keys)
    if set(got.index) != set(calculated.set_index(keys).index):
        raise ValueError(f"{sheet} endpoint coverage mismatch")
    checks = []
    numeric = [c for c in calculated.columns if c not in (*keys, "ROI group", "Unit", "Before model", "After model")]
    for row in calculated.to_dict("records"):
        key = tuple(row[k] for k in keys)
        archive = got.loc[key]
        for col in ("Unit", "Before model", "After model"):
            if archive[col] != row[col]: raise ValueError(f"{sheet} {col} mismatch at {key}")
        for col in numeric:
            expected, reported = float(row[col]), archive[col]
            if isinstance(reported, str) and reported.startswith("<") and col in ("Two-sided paired t-test P", "BH-adjusted q"):
                passed = 0 <= expected < float(reported[1:])
            else:
                tolerance = 0 if col == "Paired n" else 0.0000500001
                passed = abs(expected-float(reported)) <= tolerance
            checks.append(dict(sheet=sheet, **dict(zip(keys, key)), field=col, reported=reported,
                               recalculated=expected, within_storage_precision=bool(passed)))
    return pd.DataFrame(checks)


def check_external_summary(table, stored):
    """Recreate the legacy wide external-centre descriptive summary."""
    expected={f'{centre} — {model}' for centre in CENTRES[1:] for model in models(centre)}
    stats={'Mean':'mean','Max':'maximum','Min':'minimum','Q1':'Q1','Q2':'median','Q3':'Q3','STD':'SD'}
    require_unique(stored,['Centre/model state','Metric/statistic'],'External-centre summary')
    keys=set(map(tuple,stored[['Centre/model state','Metric/statistic']].values))
    if keys!={(c,f'{m} {s}') for c in expected for m in METRICS for s in stats}:
        raise ValueError('External summary is incomplete or contains undefined rows')
    mapping={'Small OAR':GROUP_LABELS['Small OARs'],'Large OAR':GROUP_LABELS['Large OARs'],**{r:r for r in TARGETS}}
    records=[]
    calculated=[]
    for row in stored.to_dict('records'):
        centre,model=row['Centre/model state'].split(' — ')
        metric,stat=row['Metric/statistic'].split(' ')
        state='Before' if model==models(centre)[0] else 'After'
        output={'Centre/model state':row['Centre/model state'],'Metric/statistic':row['Metric/statistic']}
        for column,roi in mapping.items():
            g=table[(table.Centre==centre)&(table.ROI==roi)&(table.Metric==metric)]
            value=float(g[f'{state} {stats[stat]}'].iloc[0])
            output[column]=value
            records.append(dict(centre=centre,model=model,roi=roi,metric=metric,statistic=stat,
                                reported=row[column],recalculated=value,
                                within_storage_precision=abs(value-float(row[column]))<=.0000500001))
        calculated.append(output)
    return pd.DataFrame(calculated),pd.DataFrame(records)


def case_values(raw, study_id="MC_R023"):
    rows = []
    for group, roi in CASE_ORDER:
        g = raw[(raw.Centre == "Main centre") & (raw["Study ID"] == study_id) & (raw.ROI == roi)]
        require_unique(g, ["Model state", "Metric"], "representative case")
        lookup = g.set_index(["Model state", "Metric"])
        if len(lookup) != 6: raise ValueError(f"Incomplete representative case: {roi}")
        record = dict(study_id=study_id, group=group, roi=roi)
        for metric in METRICS:
            before = lookup.loc[(models("Main centre")[0], metric)]
            after = lookup.loc[(models("Main centre")[1], metric)]
            m = metric.lower()
            record.update({f"baseline_{m}": float(before.Value), f"vag_{m}": float(after.Value),
                           f"beneficial_change_{m}": float((after.Value-before.Value) * (1 if metric == "DSC" else -1)),
                           f"baseline_{m}_cell": f"H{int(before.source_row)}", f"vag_{m}_cell": f"H{int(after.source_row)}"})
        rows.append(record)
    return pd.DataFrame(rows)


def analyse(source, output):
    output = Path(output)
    r, p, a, timing, strategy = [source[s] for s in ("Retrospective patient metrics", "Prospective patient metrics",
                                "Patient-level OAR averages", "Prospective timing", "Fig S4 GTVn comparison")]
    validate_metrics(r)
    validate_metrics(p, prospective=True)
    if set(map(tuple,r[['ROI','ROI code']].drop_duplicates().values)) != set(map(tuple,p[['ROI','ROI code']].drop_duplicates().values)):
        raise ValueError('ROI code/name mapping differs across cohorts')
    composite_checks = check_stored_composites(r, a)
    stored = a.rename(columns={"Within-patient OAR mean": "Value"}).copy()
    stored["ROI"] = stored["ROI group"]
    s4_input = pd.concat([r[r.ROI.isin(TARGETS)], stored], ignore_index=True)
    s4 = calculate_table(s4_input, ENDPOINTS)
    s5 = calculate_table(r, (*SMALL, *LARGE))
    checks = pd.concat([compare_summary(s4, source["Table S4"], "Table S4"),
                        compare_summary(s5, source["Table S5"], "Table S5")], ignore_index=True)
    write_csv(checks, output/"table_recalculation_checks.csv")
    if not checks.within_storage_precision.all():
        raise ValueError("Archived tables do not match patient-level recalculation; see table_recalculation_checks.csv")
    for name, frame in [("Table_S4", s4), ("Table_S5", s5), ("OAR_composite_checks", composite_checks)]:
        write_csv(frame, output/f"{name}.csv")
    # Sensitivity to using unrounded retrospective composite means, without replacing S4.
    full = calculate_table(pd.concat([r[r.ROI.isin(TARGETS)], compute_composites(r)], ignore_index=True), ENDPOINTS)
    write_csv(full, output/"Table_S4_unrounded_composites_sensitivity.csv")
    ext,ext_checks=check_external_summary(full,source['External-centre summary'])
    write_csv(ext,output/'External_centre_summary.csv')
    write_csv(ext_checks,output/'External_summary_checks.csv')
    if not ext_checks.within_storage_precision.all():
        raise ValueError('External summary does not match patient-level calculation; see External_summary_checks.csv')
    prospective = pd.concat([p, compute_composites(p)], ignore_index=True)
    descriptive = []
    for (roi, metric), frame in prospective.groupby(["ROI", "Metric"], sort=True):
        descriptive.append(dict(ROI=roi, Metric=metric, Unit="unitless" if metric == "DSC" else "mm", **describe(frame.Value)))
    descriptive = pd.DataFrame(descriptive)
    write_csv(descriptive, output/"Prospective_all_ROI_descriptives.csv")
    write_csv(descriptive[descriptive.ROI.isin(ENDPOINTS)], output/"Table_S6.csv")
    require_unique(timing, ["Study ID"], "timing")
    require_numeric(timing, TIMES, "timing")
    if set(timing["Study ID"]) != set(p["Study ID"]) or (timing[list(TIMES)] < 0).any().any():
        raise ValueError("Timing IDs or nonnegative constraints failed")
    if not np.allclose(timing[TIMES[0]]+timing[TIMES[1]], timing[TIMES[2]], atol=1e-9, rtol=0):
        raise ValueError("Derived combined time is not the sum of the two intervals")
    time_stats = pd.DataFrame([dict(endpoint=c, unit="s", **describe(timing[c])) for c in TIMES])
    write_csv(time_stats, output/"Timing_summary.csv")
    require_unique(strategy, ["Study ID", "Strategy", "Metric"], "GTVn strategies")
    require_numeric(strategy, ["Value"], "GTVn strategies")
    ids = set(r[(r.Centre == "Main centre") & (r.ROI == "GTVn")]["Study ID"])
    expected = {(i,s,m) for i in ids for s in STRATEGIES for m in METRICS}
    if set(map(tuple, strategy[["Study ID", "Strategy", "Metric"]].values)) != expected:
        raise ValueError("GTVn strategy pairs are incomplete or do not match the main GTVn cohort")
    if not strategy.Centre.eq('Main centre').all() or not strategy.ROI.eq('GTVn').all():
        raise ValueError('Unexpected GTVn strategy centre or ROI')
    reused = 0
    for row in strategy.to_dict("records"):
        if row['Unit'] != ('unitless' if row['Metric']=='DSC' else 'mm') or row['Value']<0 or (row['Metric']=='DSC' and row['Value']>1):
            raise ValueError('Invalid GTVn strategy unit or value')
        if row["Strategy"] == "Direct DIR propagation": continue
        model = models("Main centre")[0 if row["Strategy"] == "Image-only baseline" else 1]
        match = r[(r["Study ID"] == row["Study ID"]) & (r.ROI == "GTVn") & (r.Metric == row["Metric"]) & (r["Model state"] == model)]
        if len(match) != 1 or float(match.Value.iloc[0]) != row["Value"]:
            raise ValueError("GTVn strategy sheet does not match retrospective data")
        reused += 1
    gs = pd.DataFrame([dict(strategy=s, metric=m, **describe(f.Value)) for (s,m),f in strategy.groupby(["Strategy","Metric"])])
    write_csv(gs, output/"GTVn_strategy_descriptives.csv")
    for row in source['GTVn stored summary'].to_dict('records'):
        calc=gs[(gs.strategy==row['Strategy']) & (gs.metric==row['Metric'])]
        if len(calc)!=1 or int(calc.n.iloc[0])!=row['n'] or abs(float(calc['mean'].iloc[0])-row['Mean'])>5.00001e-5 or abs(float(calc.sd.iloc[0])-row['Sample SD'])>5.00001e-5:
            raise ValueError('GTVn stored summary mismatch')
    write_csv(case_values(r), output/"FigureS5_case_values.csv")
    families=[]
    for name,table in [('Table S4',s4),('Table S5',s5)]:
        for centre,f in table.groupby('Centre',sort=False):
            effects=np.array([orient_effect(row['Mean difference (after-before)'],row['95% CI lower'],row['95% CI upper'],row['Metric'])[0] for row in f.to_dict('records')])
            sig=f['BH-adjusted q'].to_numpy()<.05
            families.append(dict(table=name,centre=centre,family_size=len(f),significant_improvements=int(np.sum(sig & (effects>0))),significant_deteriorations=int(np.sum(sig & (effects<0))),not_significant=int(np.sum(~sig))))
    write_csv(families, output/'Multiplicity_families.csv')
    summary=dict(**source["metadata"], records={s:len(source[s]) for s in source if isinstance(source[s],pd.DataFrame)},
                 table_cells_verified=len(checks), composite_values_verified=len(composite_checks),
                 external_summary_values_verified=len(ext_checks),
                 GTVn_reused_values_verified=reused, patients_prospective=len(timing),
                 sample_size=precision_sample_size(),
                 method="Patient pairing; two-sided paired t; pointwise t CI; BH within centre and table (18/63 tests)",
                 raw_imaging_verified=False)
    write_json(summary, output/"analysis_validation.json")
    return summary
