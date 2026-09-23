# Analysis specification

## Input modes

Bundled mode is the plotting default. It validates and displays the included
numerical subset and retains the stored S4/S5 statistical results. It does not
reconstruct retrospective tests from summary values. The complete patient-level
checks and statistical recalculations described below run when a local workbook
is explicitly supplied with `--source-data`.

The same rendering functions receive normalized records from either input mode.
No automatic workbook discovery or fallback occurs. Output metadata identifies
the selected mode and whether complete patient-level analysis was performed.

## Observation unit and pairing

Each observation is identified by centre, model state, study ID, ROI and metric.
Raw ROI codes and names must correspond one-to-one. Paired testing matches study
IDs within centre/ROI/metric; it never pairs observations by row position or
metric rank. Missing records, missing values, duplicate identifiers, undefined
ROIs/models/units and inconsistent counts cause explicit errors. There is no
imputation, outlier deletion, interpolation between patients or pooling of
organs as if they were independent patients.

The main-centre GTV test cohort contains 12 patients; the CTV/OAR test cohort
contains 29 distinct patients. External held-out cohorts have 14, 21, 14 and 25
patients. The prospective cohort has 135 patients with 25 ROIs × 3 metrics each.
Main-centre paired states are conventional U-Net and locked source model;
external paired states are pre-adaptation and locally adapted model.

## OAR composites and stored precision

For patient i and a fixed group of K organs, the composite is
`mean_i = sum(metric_i,k, k=1..K) / K`.
K is 10 for small OARs and 11 for large OARs. All constituent organs are required.
Categories are fixed, not assigned anew from a patient's predicted volume.
The full membership list is in `oart_analysis/constants.py`.

S4 uses the archived patient composites in `Patient-level OAR averages`, which
were stored to four decimal places. These are checked against the constituent
means with a half-unit-in-the-last-place tolerance of approximately 0.00005.
The primary S4 calculation does not substitute unrounded means. An unrounded
alternative is exported separately and labelled as sensitivity analysis.
The `External-centre summary` auxiliary worksheet uses unrounded patient OAR
means; its 1,008 descriptive values are checked using that convention, and a
recalculated wide summary is exported. It is not used to replace S4's primary
values for the figures.
Prospective group means use full stored ROI precision without intermediate
rounding. In both cohorts, group summaries use the distribution of patient
composites, not a pooled patient-by-organ distribution.

## Descriptive and inferential statistics

For n observations, mean is the arithmetic mean; sample variance divides by
`n−1`. `SEM = SD/sqrt(n)`; `CV = SD/mean` is dimensionless and is left undefined
when mean is zero. Quartiles and percentiles use linear interpolation
(`numpy.quantile(..., method='linear')`). All observations are retained.

For paired differences `d_i = after_i − before_i`,
`t = mean(d) / (SD(d)/sqrt(n))`, with n−1 degrees of freedom. The two-sided P
value is `2 * t.sf(abs(t), n−1)`. A pointwise 95% CI is
`mean(d) ± t.ppf(0.975,n−1)*SD(d)/sqrt(n)`.
If every difference is exactly zero, P=1. If every difference is the same
nonzero number, the limiting P=0 is recorded and the t statistic is undefined;
the data should be reviewed for the reason for the degenerate distribution.

BH correction sorts m P values and uses
`q_(i) = min(1, min_{j≥i}(m*p_(j)/j))`, then restores the original endpoint order.
Families are prespecified **within each centre**: 18 group/target/metric
comparisons in S4; 63 organ/metric comparisons in S5. The criterion is q<0.05.
Pointwise CIs are not multiplicity-adjusted. Their overlap with zero need not
give the same decision as q<0.05.

All table effects are after−before. For HD95/ASD improvement plots,
`effect = −difference`, `lower = −original_upper`, `upper = −original_lower`.
DSC improvement keeps the original sign and endpoints. Fixed plotting limits
are checked so complete displayed ranges/CIs are retained.

Archived table means, SDs, differences and CIs are checked within 0.00005
storage rounding. P/q strings such as `<0.0001` are interpreted as bounds, not
as exact P values. Calculated P/q values are exported without display rounding.

## Prospective and case reporting

Prospective summaries include n, mean, SD, SEM, CV, quantiles and pointwise mean
CI. These are descriptive; no unplanned group hypothesis tests are added.
Fig. 6a shows P5–P95, IQR and the median. Table_S6.csv contains full numerical
summaries for the six group/target endpoints and three metrics; the manuscript
can select the relevant reporting columns.

Timing is checked per patient: derived sum = algorithm interval + review and
modification interval. Means, SDs and quantiles are calculated from each
patient's values, including the paired sum. SDs are not added. Ranked display
positions are 1–135, not percentiles, chronological order or treatment fractions.
The two recorded intervals can overlap, so their sum is not the actual
end-to-end elapsed time or active-editing-only time.

The prospective sample-size calculation uses SD=0.0653, 95% normal confidence,
half-width 0.012 and 15% non-evaluable allowance. It rounds the required
evaluable count up first, then inflates and rounds up the recruitment count.
The output is 114 evaluable and 135 recruited patients.

Fig. S4 includes all 12 patients under the three GTVn strategies. Baseline and
PCG-UNet values are checked against the retrospective sheet. Its summaries are
descriptive and no additional P/q values are created.

Fig. S5 uses MC_R023. Its 21 organ pairs are read from the selected input source;
the figure annotations are calculated from those pairs. The three header means
average the 21 organs within this case and are not cohort results. No case-level
P values or CIs are shown.

## Limits of verification

The programme verifies the numerical relationships available in the workbook;
it does not certify that all source measurements are correct. Training,
inference, image registration, original masks and unrecorded clinical metadata
cannot be reconstructed from metric-level tables. The independent mask
evaluator needs separate local images and an explicitly documented surface
convention; see [IMAGING.md](IMAGING.md).

These analyses do not infer clinical outcomes, dosimetry or causal treatment
benefit from geometric contour concordance.
