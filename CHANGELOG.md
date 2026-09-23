# Changelog

## 1.1.0

- Include six real plotting datasets with field descriptions and integrity
  checks, enabling figure reproduction without a workbook.
- Use bundled data by default; select a local workbook explicitly with
  `--source-data`. Invalid selected workbooks do not trigger a fallback.
- Keep complete patient-level reanalysis in workbook mode and report the data
  source and validation scope for every run.
- Package CSV resources for use after installation, with tests for input modes,
  data integrity and installed plotting.
- Use the project name `npc-oart-automated-contouring` and consistent clinical
  and technical terminology.
- Document commercial licensing and corresponding-author access requests for
  PCG-UNet, PAC-UNet and VAG-UNet.

The complete workbook, patient imaging and proprietary networks are excluded.
