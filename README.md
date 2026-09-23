# NPC oART Automated Contouring

Analysis and figure reproduction for automated contouring in nasopharyngeal carcinoma (NPC) during online adaptive radiotherapy (oART).

**Version 1.1.0** · [Data dictionary](data/README.md) · [Statistical methods](docs/METHODS.md) · [Imaging metrics](docs/IMAGING.md) · [Figure mapping](docs/FIGURE_MAPPING.md)

## Quick start: no workbook required

The repository includes the **real numerical data used by the included figures**. They are stored as six CSV files in `oart_analysis/data/`, with original storage precision and worksheet/cell references where applicable. They are not simulated examples. The complete `Source_data.xlsx` workbook and patient images are not included.

Python 3.11 or later is required; Python 3.12 was used for local validation. From the repository root:

```bash
python -m venv .venv
# Activate .venv using the command for your operating system.
python -m pip install -e ".[test]"
python run_all.py
```

This generates all six included quantitative figure sets using the bundled CSVs. No workbook, DICOM images, GPU or uTPS installation is required for this plotting workflow. Dependencies must be installed first; the scripts do not download or upload data.

To use the versions tested locally, install `requirements-tested.txt` before the editable package. It records one tested Python 3.12 environment, including optional imaging/test dependencies.

## Two data modes

| Mode | How to select it | What it does |
|---|---|---|
| **Bundled plotting data (default)** | Omit `--source-data` | Validates bundled data and reproduces every included figure |
| **Local workbook** | Pass `--source-data "/path/to/Source_data.xlsx"` | Performs complete patient-level workbook analysis and validation, then reproduces the figures |

```bash
# Default: use bundled data, even if a workbook is present locally
python run_all.py --output outputs/bundled

# Explicitly select a local workbook
python run_all.py --source-data "/path/to/Source_data.xlsx" --output outputs/workbook

# Reproduce one figure from bundled data
python scripts/plot_figure5.py
python -m oart_analysis.cli plot --figure S7

# Reproduce one figure from a workbook
python scripts/plot_figure5.py --source-data "/path/to/Source_data.xlsx"

# Complete statistical reanalysis requires the local workbook
python -m oart_analysis.cli analyse --source-data "/path/to/Source_data.xlsx"
```

The default mode does not search for `Source_data.xlsx`. If a workbook is explicitly selected but is missing or fails validation, the command reports an error; it does not silently use bundled data. Workbooks are opened read-only.

Each run prints its data mode and writes `data_source.json`, including input hashes, validation scope and whether full patient-level analysis was performed. Each figure folder contains the same mode record. Use separate output directories when comparing modes.

**Scope of verification:** the bundled subset contains S4/S5 plotting summaries, prospective patient endpoints and timing, GTVn strategy records, and the selected OAR case. It reproduces the plots, their displayed intervals and significance markers, and the descriptive summaries calculated from those inputs. It does not contain all retrospective patient pairs or all constituent prospective organ measurements. Recomputing the complete retrospective P/q analyses, verifying all OAR constituent means and regenerating all statistical tables therefore requires the local workbook. The default mode does not present stored P/q values as newly calculated tests.

## Included figures and outputs

| Folder under `outputs/figures/` | Included panels |
|---|---|
| `Figure5/` | 5b–c: DSC means and paired improvements across all five centres |
| `Figure6/` | 6a–b: prospective geometric agreement and all 135 ranked patient time intervals |
| `FigureS4/` | S4b–d: three-strategy GTVn comparison in 12 paired patients |
| `FigureS5/` | S5: 21 OARs in the selected MC_R023 case |
| `FigureS6/` | S6a–d: HD95/ASD absolute means and paired reductions |
| `FigureS7/` | S7a–c: organ-level paired changes and confidence intervals |

Every figure is exported as editable SVG, PDF, PNG preview and TIFF, with the plotted numerical values saved beside it. The complete workbook mode additionally writes S4/S5/S6 calculations, OAR composite checks, the external-centre summary, timing/strategy summaries, multiplicity families and validation reports to `outputs/analysis/`.

Fig. 2, Fig. 3 and Fig. 4g are excluded. Non-data workflow, mechanism and network diagrams and clinical image plates are also excluded, including Fig. 5a and Fig. S4a. Fig. 5 and Fig. S4 outputs are quantitative panel sets for assembly with separate artwork. Final manuscript panel lettering is retained.

## Statistical and display conventions

- Patients are the analysis unit. Workbook pairing uses study identifiers, not row positions or metric ranks. Missing pairs and duplicate keys cause errors.
- Small and large OAR composites are equal-weight means of fixed groups of 10 and 11 organs. S4 uses the patient means stored to four decimals. Prospective composites and the auxiliary external-centre summary use unrounded constituent means. See [METHODS.md](docs/METHODS.md).
- Workbook analyses use two-sided paired t-tests, sample SD, pointwise 95% t confidence intervals and BH correction within each centre: 18 tests for S4 and 63 for S5. Bundled figure summaries retain these reported results.
- Table differences are after minus before. DSC gain keeps that sign; HD95/ASD reduction reverses it and negates/exchanges the CI endpoints.
- Fig. 6a shows the median, IQR and 5th–95th percentiles. These ranges are not confidence intervals.
- Fig. 6b ranks patients 1–135 by their derived timing sum. The recorded intervals can overlap, so the sum is not end-to-end elapsed time or active editing time.
- Fig. S5 is a selected single case and carries no cohort inference. Its header means average the 21 organs within that case. GTVn strategy plots are descriptive.
- No automatic outlier removal, imputation or data alteration is performed.

## DSC, HD95 and ASD from local imaging

The optional imaging module evaluates already exported masks; it does not run the proprietary segmentation networks.

```bash
python -m pip install -e ".[imaging]"
# Copy config/imaging.example.json to config/imaging.local.json,
# then replace the input paths and ROI names/labels with your own values.
python -m oart_analysis.cli metrics --config config/imaging.local.json --output outputs/imaging
```

Inputs can be geometrically aligned scalar labelmaps or RTSTRUCTs on one CT series. Bilateral components are unioned before evaluation. The surface, physical spacing, pooled bidirectional distance and empty-mask conventions are documented in [IMAGING.md](docs/IMAGING.md). Original patient images are required to investigate agreement with the study's archived metrics; the bundled numerical data alone cannot verify the original masks.

## Algorithm availability

The automated contouring algorithms, including PCG-UNet, PAC-UNet and VAG-UNet, are available within the licensed automated contouring module in uTPS (United Imaging Healthcare, Shanghai, China). Due to commercial restrictions, the algorithm implementations and model weights are not publicly distributed through this repository. Requests for access for academic research and independent verification should be directed to the corresponding author. Following review of the request, the corresponding author will seek approval from United Imaging Healthcare. The scope and terms of any access granted will be subject to approval by United Imaging Healthcare and the applicable licensing conditions.

Please use the corresponding-author contact details listed in the associated manuscript.

## Tests and installation checks

```bash
python -m pytest
python -m oart_analysis.cli sample-size
```

Tests cover statistics, patient pairing, bundled-data integrity, mode selection and synthetic mask geometry. Optional imaging tests require the `[imaging]` extra. The GitHub Actions workflow tests the installed package, including its packaged CSV resources, without the complete workbook.

## Repository layout

```text
oart_analysis/         workbook analysis, plotting and segmentation metrics
oart_analysis/data/    six real plotting CSVs and their integrity manifest
scripts/              individual figure, analysis and mask-evaluation entry points
config/               local imaging path/ROI example
data/README.md        dataset and field descriptions
docs/                 methods, figure mapping, terminology and validation scope
tests/                numerical, input-mode and synthetic imaging tests
```

Clinical system and workflow descriptions use **automated contouring**; network/mask operations and metric functions use **segmentation**. See [TERMINOLOGY.md](docs/TERMINOLOGY.md).

The software is released under [Apache License 2.0](LICENSE). The software licence does not grant rights in proprietary uTPS software or model weights. The bundled study data are supplied to support reproduction of the included figures.

When using this repository, cite the associated study and specify the code version. No provisional article DOI is assigned in this repository.
