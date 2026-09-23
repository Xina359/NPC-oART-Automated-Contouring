# Final manuscript figure and table mapping

Default plotting inputs are the six bundled CSVs documented in
[`data/README.md`](../data/README.md). The workbook column in the table below
describes the explicit `--source-data` mode. Both modes use the same rendering
functions and preserve the final panel letters. Complete table recomputation
requires workbook mode.

| Manuscript item | Workbook sheet / derivation | Entry point |
|---|---|---|
| Fig. 5b | Table S4, DSC before/after means | `scripts/plot_figure5.py` |
| Fig. 5c | Table S4, DSC mean paired differences and pointwise CI, all five centres | same |
| Fig. 6a | Prospective patient metrics; 4 targets and equal-weight within-patient 10/11-OAR means | `scripts/plot_figure6.py` |
| Fig. 6b | Prospective timing; all 135 patients ordered by derived sum then ID | same |
| Fig. S4b–d | Fig S4 GTVn comparison, columns A:G; all three strategies for each patient | `scripts/plot_figureS4.py` |
| Fig. S5 | Retrospective patient metrics; MC_R023, 21 OARs, 3 metrics, both model states | `scripts/plot_figureS5.py` |
| Fig. S6a–d | Table S4, HD95/ASD means and sign-reversed paired differences/CI | `scripts/plot_figureS6.py` |
| Fig. S7a–c | Table S5, individual OAR differences/CI and centre-specific BH q | `scripts/plot_figureS7.py` |
| Table S4 | Retrospective patient metrics + Patient-level OAR averages | `analyse` command |
| Table S5 | Retrospective patient metrics; individual organs | `analyse` command |
| Table S6 | Prospective patient metrics; full-precision composites | `analyse` command |

The archived S4/S5 values are used for figure positions. Workbook mode verifies
them against patient-level recalculation; bundled mode checks the stored plotting
subset without repeating those tests. Calculated tables in workbook mode are
exported at full numerical precision; apply final display rounding when
assembling manuscript tables.

The code excludes Fig. 2, Fig. 3 and Fig. 4g. It also excludes non-data diagrams
and image plates (including Fig. 5a and Fig. S4a). These are not missing data
panels in this package. Final panel letters are retained in partial exports.

Small/large OAR grouping is fixed by the study's organ categories. The five
anatomical display groups in Fig. S5 are a layout classification and do not
replace the small/large OAR statistical groups.
