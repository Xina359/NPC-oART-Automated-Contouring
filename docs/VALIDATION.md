# Validation scope

Local verification used the current `Source_data.xlsx` supplied for the study.
The workbook was opened read-only. Its SHA-256 digest is recorded in the local
runtime report rather than embedded as an input restriction in the public code.

- All 90 S4 and 315 S5 comparisons were recomputed from patient-level values;
  8,505 numerical table fields agreed within storage precision, including the
  specified P/q bounds.
- All 1,236 stored patient OAR composites agreed with their complete constituent
  organs within four-decimal storage precision.
- All 1,008 external-centre descriptive values were verified with unrounded
  composites, as used by that auxiliary worksheet.
- GTVn strategy coverage and the 72 deliberately shared baseline/PCG-UNet
  measurements were checked across sheets.
- The prospective performance and timing sheets contained the same 135
  patients. Derived sums were verified per patient.
- All six quantitative figure sets generated successfully in SVG, PDF, PNG and
  TIFF. The local preview was visually inspected. Source positions/intervals
  are exported alongside plots for numerical checking.
- Synthetic tests exercise statistics, pairing, input validation and mask
  geometry without study data. Run `python -m pytest` for the current test count.

The release package includes the six plotting datasets. It excludes the complete
workbook, patient imaging and generated validation/figure output folders. The
bundled data preserve the numeric values used to draw the included figures;
complete workbook-based analyses are available through `--source-data`.

Version 1.1.0 validation additionally checks the CSV precision round trip,
bundled-versus-workbook figure equivalence, explicit input-mode selection,
checksum failures and package-data availability after installation. Both input
modes use the same rendering functions. Bundled-mode reports do not claim that
all patient-level retrospective tests were repeated.

- All 13,867 numeric fields in the exported plotting datasets survived the CSV
  round trip with exact numeric equality.
- All six PNG previews were pixel-identical between bundled and workbook modes
  in the same local environment. All six SVGs retained editable text and vector
  elements without embedded raster images.
- All 25 automated tests passed, including explicit-input failures and bundled
  dataset integrity checks.
- A built wheel was installed separately and successfully loaded all six
  datasets and plotted Fig. S4 from an empty working directory.

These checks establish executable numerical consistency, not correctness of
unavailable original imaging or clinical records. There was no raw-patient-image
validation. The RTSTRUCT round-trip test uses a simple synthetic CT and ROI.
The shipped GitHub Actions workflow has not been run on a remote repository as
part of this local delivery.
