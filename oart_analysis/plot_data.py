"""One plotting-data interface for packaged CSVs and an explicitly supplied workbook."""
from collections import defaultdict
from importlib.resources import files
import csv
import hashlib
import io
import itertools
import json
import math
import numbers
import numpy as np
from .analysis import case_values
from .constants import CENTRES, METRICS, ENDPOINTS, SMALL, LARGE, GROUPS, TARGETS, STRATEGIES, models, expected_n
from .plotting import summary_rows

DATASETS = {
    'summary_s4': 'table_s4_plot.csv', 'summary_s5': 'table_s5_plot.csv',
    'performance': 'prospective_endpoints.csv', 'timing': 'prospective_timing.csv',
    'gtvn': 'gtvn_strategies.csv', 'case': 'representative_oar_case.csv',
}


def from_workbook(source):
    """Convert a validated workbook to plotting records without rounding values."""
    mapping = {
        'centre': 'Centre', 'roi': 'ROI', 'metric': 'Metric', 'n': 'Paired n', 'unit': 'Unit',
        'before_model': 'Before model', 'after_model': 'After model', 'before_mean': 'Before mean',
        'after_mean': 'After mean', 'before_sd': 'Before SD', 'after_sd': 'After SD',
        'mean_difference': 'Mean difference (after-before)', 'ci_lower': '95% CI lower',
        'ci_upper': '95% CI upper', 'p_value': 'Two-sided paired t-test P', 'q_value': 'BH-adjusted q',
    }
    s4 = []
    for record in source['Table S4'].to_dict('records'):
        row = {key: record[column] for key, column in mapping.items()}
        n = int(record['source_row'])
        row.update(source_sheet='Table S4', source_row=n, before_cell=f'H{n}', after_cell=f'P{n}')
        s4.append(row)
    raw = source['Prospective patient metrics'].to_dict('records')
    lookup = defaultdict(list)
    for row in raw:
        lookup[(row['Study ID'], row['Metric'])].append(row)
    performance = []
    for study_id in sorted({row['Study ID'] for row in raw}):
        for metric in METRICS:
            rows = lookup[(study_id, metric)]
            for category in (*TARGETS, 'Large OARs', 'Small OARs'):
                group = category in GROUPS
                matches = [r for r in rows if r['ROI group' if group else 'ROI'] == category]
                count = len(GROUPS[category]) if group else 1
                if len(matches) != count or (group and {r['ROI'] for r in matches} != set(GROUPS[category])):
                    raise ValueError(f'Incomplete prospective endpoint: {study_id}/{category}/{metric}')
                performance.append({
                    'study_id': study_id, 'category': category,
                    'category_type': 'Grouped OARs' if group else 'Target', 'metric': metric,
                    'value': float(np.mean([float(r['Value']) for r in matches])),
                    'unit': matches[0]['Unit'], 'n_components': count,
                    'source_sheet': 'Prospective patient metrics',
                    'source_rows': ';'.join(str(r['source_row']) for r in matches),
                    'source_cells': ';'.join(f"H{r['source_row']}" for r in matches),
                    'derivation': f'within-patient mean across {count} OARs' if group else 'direct ROI value',
                })
    timing = []
    for record in source['Prospective timing'].to_dict('records'):
        n = int(record['source_row'])
        timing.append(dict(study_id=record['Study ID'],
                           algorithm_s=float(record['Algorithm-processing time (s)']),
                           review_s=float(record['Review-and-modification interval (s)']),
                           combined_s=float(record['Derived combined time (s)']),
                           source_sheet='Prospective timing', source_row=n,
                           algorithm_cell=f'B{n}', review_cell=f'C{n}', combined_cell=f'D{n}'))
    data = dict(summary_s4=s4, summary_s5=summary_rows(source['Table S5']), performance=performance,
                timing=timing, gtvn=source['Fig S4 GTVn comparison'].to_dict('records'),
                case=case_values(source['Retrospective patient metrics']).to_dict('records'))
    validate(data)
    return data, dict(input_mode='workbook', **source['metadata'], full_patient_level_analysis=True)


def _unique(rows, keys, label):
    actual = [tuple(r[k] for k in keys) for r in rows]
    if len(actual) != len(set(actual)) or any(v is None or str(v).strip() == '' for key in actual for v in key):
        raise ValueError(f'Missing or duplicate keys: {label}')
    return set(actual)


def _finite(value):
    return isinstance(value, numbers.Real) and not isinstance(value, (bool, np.bool_)) and math.isfinite(value)


def _metric(value, metric, unit):
    if metric not in METRICS or not _finite(value) or value < 0 or (metric == 'DSC' and value > 1):
        raise ValueError('Invalid plotting metric value')
    if unit != ('unitless' if metric == 'DSC' else 'mm'):
        raise ValueError('Invalid plotting metric unit')


def _probability(value):
    number = float(str(value).lstrip('<'))
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise ValueError('Invalid P/q value in plotting summary')
    return number


def validate(data):
    """Check the plotting subset; this does not reconstruct omitted patient pairs."""
    if set(data) != set(DATASETS): raise ValueError('Expected exactly six plotting datasets')
    for name, rois in (('summary_s4', ENDPOINTS), ('summary_s5', (*SMALL, *LARGE))):
        rows = data[name]
        keys = _unique(rows, ('centre', 'roi', 'metric'), name)
        if keys != set(itertools.product(CENTRES, rois, METRICS)):
            raise ValueError(f'Incomplete endpoint coverage: {name}')
        for row in rows:
            if row['n'] != expected_n(row['centre'], row['roi']): raise ValueError('Unexpected paired sample size')
            if (row['before_model'], row['after_model']) != models(row['centre']): raise ValueError('Unexpected model comparison')
            for state in ('before', 'after'):
                _metric(row[state+'_mean'], row['metric'], row['unit'])
                if not _finite(row[state+'_sd']) or row[state+'_sd'] < 0: raise ValueError('Invalid sample SD')
            if not all(_finite(row[k]) for k in ('mean_difference', 'ci_lower', 'ci_upper')): raise ValueError('Invalid effect or CI')
            if not row['ci_lower'] <= row['mean_difference'] <= row['ci_upper']: raise ValueError('CI does not bracket its estimate')
            if abs(row['after_mean']-row['before_mean']-row['mean_difference']) > .00015001:
                raise ValueError('Mean difference is inconsistent with stored means')
            for key in ('p_value','q_value'): _probability(row[key])
    performance, timing = data['performance'], data['timing']
    ids = {r['study_id'] for r in performance}
    keys = _unique(performance, ('study_id','category','metric'), 'prospective endpoints')
    if len(ids) != 135 or keys != set(itertools.product(ids, (*TARGETS, *GROUPS), METRICS)):
        raise ValueError('Expected 135 complete prospective patients')
    for row in performance:
        _metric(row['value'], row['metric'], row['unit'])
        if row['n_components'] != (len(GROUPS[row['category']]) if row['category'] in GROUPS else 1):
            raise ValueError('Incorrect prospective OAR component count')
    if _unique(timing, ('study_id',), 'timing') != {(p,) for p in ids}: raise ValueError('Performance and timing patient IDs differ')
    for row in timing:
        if not all(_finite(row[k]) and row[k] >= 0 for k in ('algorithm_s','review_s','combined_s')): raise ValueError('Invalid timing value')
        if abs(row['algorithm_s']+row['review_s']-row['combined_s']) > 1e-9: raise ValueError('Incorrect patient timing sum')
    strategies = data['gtvn']
    gtvn_ids = {r['Study ID'] for r in strategies}
    if len(gtvn_ids) != 12 or _unique(strategies, ('Study ID','Strategy','Metric'), 'GTVn') != set(itertools.product(gtvn_ids, STRATEGIES, METRICS)):
        raise ValueError('Incomplete GTVn strategy pairs')
    for row in strategies:
        _metric(row['Value'], row['Metric'], row['Unit'])
        if row['Centre'] != 'Main centre' or row['ROI'] != 'GTVn': raise ValueError('Unexpected GTVn cohort')
    case = data['case']
    if _unique(case, ('study_id','roi'), 'representative case') != {('MC_R023', r) for r in (*SMALL, *LARGE)}:
        raise ValueError('Expected the 21 OARs of MC_R023')
    for row in case:
        for metric in METRICS:
            m = metric.lower()
            for state in ('baseline','vag'): _metric(row[state+'_'+m], metric, 'unitless' if metric == 'DSC' else 'mm')
            calculated = (row['vag_'+m]-row['baseline_'+m]) * (1 if metric == 'DSC' else -1)
            if not _finite(row['beneficial_change_'+m]) or abs(calculated-row['beneficial_change_'+m]) > 1e-12:
                raise ValueError('Inconsistent representative-case change')
    return {'dataset_rows': {key:len(rows) for key,rows in data.items()}, 'checks': 'passed',
            'scope': 'plotting subset integrity; not a full patient-level retrospective reanalysis'}


def load_bundled():
    """Read installed package resources, independent of the current directory."""
    directory = files('oart_analysis').joinpath('data')
    manifest_bytes = directory.joinpath('manifest.json').read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest['schema_version'] != 1 or set(manifest['datasets']) != set(DATASETS): raise ValueError('Unsupported bundled-data manifest')
    result = {}
    for name, filename in DATASETS.items():
        spec = manifest['datasets'][name]
        if spec['file'] != filename: raise ValueError('Unexpected bundled dataset filename')
        raw = directory.joinpath(filename).read_bytes()
        if hashlib.sha256(raw).hexdigest() != spec['sha256']: raise ValueError(f'Bundled dataset checksum mismatch: {filename}')
        reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
        if reader.fieldnames != list(spec['types']): raise ValueError(f'CSV schema mismatch: {filename}')
        rows = []
        for record in reader:
            row = {}
            for key, kind in spec['types'].items():
                value = record[key]
                if kind == 'float': value = float(value)
                elif kind == 'int':
                    number = float(value)
                    if not math.isfinite(number) or not number.is_integer(): raise ValueError('Noninteger count or source row')
                    value = int(number)
                elif kind != 'str': raise ValueError(f'Unsupported CSV field type: {kind}')
                row[key] = value
            rows.append(row)
        if len(rows) != spec['rows']: raise ValueError(f'CSV row count mismatch: {filename}')
        result[name] = rows
    validate(result)
    return result, dict(input_mode='bundled', dataset_version=manifest['dataset_version'],
                        manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                        source_workbook_sha256=manifest['source_workbook_sha256'], full_patient_level_analysis=False)
