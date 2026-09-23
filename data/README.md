# Plotting data and workbook input

The real plotting datasets are packaged inside `oart_analysis/data/` so they are
available both in a cloned repository and after installation. The complete
`Source_data.xlsx` workbook and patient imaging files are not included.

| CSV file | Rows | Observation unit | Used by |
|---|---:|---|---|
| `table_s4_plot.csv` | 90 | Centre × group/target × metric summary | Fig. 5b–c, S6a–d |
| `table_s5_plot.csv` | 315 | Centre × individual OAR × metric summary | Fig. S7a–c |
| `prospective_endpoints.csv` | 2430 | Patient × 6 endpoints × 3 metrics | Fig. 6a |
| `prospective_timing.csv` | 135 | One patient's recorded intervals and derived sum | Fig. 6b |
| `gtvn_strategies.csv` | 108 | Patient × 3 strategies × 3 metrics | Fig. S4b–d |
| `representative_oar_case.csv` | 21 | One OAR in case MC_R023, with both states and three metrics | Fig. S5 |

These files contain study measurements and derived plotting summaries, not
synthetic demonstration data. Study IDs are the coded identifiers used by the
study workbook. No names, dates of birth, patient imaging or clinical records
are included.

## Field definitions

- `centre`, `roi`, `metric`, `before_model`, `after_model`: the comparison being
  displayed. DSC is unitless; HD95 and ASD are in millimetres. `n` is the number
  of paired patients, not the number of structures.
- `before_mean`, `after_mean`, `before_sd`, `after_sd`: patient-level means and
  sample SDs. Other Table S5 summary columns retain min/max/quartiles.
- `mean_difference`, `ci_lower`, `ci_upper`: after minus before and its
  pointwise 95% confidence interval. Distance improvement plots reverse these
  signs and exchange interval endpoints during rendering.
- `p_value`, `q_value`: archived paired t-test P and centre-specific BH q. A
  leading `<` is a bound, not an exact probability. S5's `p_numeric` and
  `q_numeric` fields use that bound numerically for display classification;
  they are not substituted for exact values in new statistical testing.
- `study_id`, `category`, `value`, `n_components`: prospective patient endpoint
  records. `n_components` is 1 for a target and 10 or 11 for an OAR composite.
  Composite values retain the full precision used by the figure. Their
  constituent organ records are available through workbook mode.
- `algorithm_s`, `review_s`, `combined_s`: seconds. The last is the patient-level
  sum of the first two and is not end-to-end elapsed time.
- GTVn strategy rows retain the workbook's English column names (`Study ID`,
  `Strategy`, `Metric`, `Value`, etc.). The same patients appear in all strategies.
- Representative-case `baseline_*` and `vag_*` are absolute values;
  `beneficial_change_*` is after minus before for DSC and before minus after for
  HD95/ASD. Anatomical `group` labels organise the display; they are distinct
  from the fixed small/large OAR statistical groups.
- `source_sheet`, `source_row`, `source_rows` and `*_cell` fields identify the
  worksheet locations or constituent cells used to form a plotting value.

## Precision and integrity

`manifest.json` records file names, row counts, field types, dataset version,
SHA-256 checksums and the source-workbook digest. CSV serialization preserves
the underlying binary64 values without an extra display-rounding step. Source
summary values retain the precision originally stored in the workbook.

The loader checks each file's checksum and validates endpoint coverage,
sample counts, units, unique keys, interval ordering, timing sums and
representative-case changes. The bundled subset cannot independently rebuild
all retrospective P/q calculations, because most underlying retrospective
patient pairs are not included.

## Choosing input data

`python run_all.py` uses these bundled CSVs. It does not inspect a local
workbook automatically. To use your workbook, explicitly run:

```bash
python run_all.py --source-data "/path/to/Source_data.xlsx" --output outputs/workbook
```

The selected workbook must retain the required sheet names and headers. A
missing or invalid selected workbook produces an error with no silent fallback.
Full statistical reanalysis also requires this explicit workbook path.
