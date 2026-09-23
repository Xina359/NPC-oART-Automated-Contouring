# Calculate DSC, HD95 and ASD from local masks

Install `python -m pip install -e ".[imaging]"`. Copy
`config/imaging.example.json` to `config/imaging.local.json` and replace the
`REPLACE_...` values with paths and ROI definitions. Keep the input-format
example you need; delete the unused example. Add more cases/ROIs as required.

```bash
python -m oart_analysis.cli metrics --config config/imaging.local.json --output outputs/imaging
```

## Supported inputs

- **Labelmaps:** scalar 3-D NIfTI, NRRD or other formats supported by SimpleITK.
  Reference and prediction must have identical size, spacing, origin and
  direction. Supply a list of positive integer label IDs for each ROI. Listing
  left and right IDs unions the masks before calculation. A binary mask uses
  `[1]`. Probability maps must be thresholded explicitly before use.
- **DICOM RTSTRUCT:** provide a directory containing one CT series, the
  reference and prediction RTSTRUCT paths, and exact ROI names. Keep RTSTRUCTs
  outside the CT-only directory. Both must refer to the same CT frame of
  reference. The adapter uses rt-utils to rasterize contours onto that CT grid.
  Bilateral ROI name lists are unioned first. Irregular slice positions,
  ambiguous frames, sheared grids and missing/empty structures are rejected.
  DICOM SEG and multiple-series directories are not supported by this adapter;
  export a geometrically verified labelmap instead.

No automatic registration or resampling is performed. If images occupy
different grids, prepare and review the alignment before running evaluation.
For RTSTRUCTs, inspect the rasterized masks locally against the CT, especially
for contour holes, small structures and unusual contour geometry. The package
does not replace clinical contour review.

## Numerical definitions

DSC is `2 * |prediction ∩ reference| / (|prediction| + |reference|)`.

The surface is the set of foreground voxel centres removed by a single
six-neighbour binary erosion, with the region outside the array considered
background. Nearest-neighbour Euclidean distances are calculated in physical
millimetres in both directions. If A and B are the two directed distance
arrays:

- HD95 = the 95th percentile of the concatenated array `[A, B]`, using linear
  interpolation.
- ASD = `(sum(A) + sum(B)) / (len(A) + len(B))`.

The latter weights the two directions by their surface-point counts. It is not
the unweighted average of the two directional means when those counts differ.
Each surface voxel is weighted equally; no mesh-area/surfel weighting is used.
The maximum of two separate directional 95th percentiles is also a different
HD95 definition and is not used here.

Array axis order matters: SimpleITK labelmaps are z,y,x, so x,y,z image spacing
is reversed; rt-utils masks are row,column,slice, so spacing follows that order.
Both mask grids and their physical metadata are checked before comparison.
Volume is reported as voxel count × voxel volume / 1000 in cm³.

Empty reference or prediction masks raise an error. The code does not assign
an arbitrary distance, remove the patient, or treat a missing ROI as zero.
All evaluated masks in the reported study were nonempty; another dataset may
need a separate, prespecified empty-mask policy.

## Interpretation and testing

This module provides an explicit evaluation convention and configurable file
reading. Metric values alone in `Source_data.xlsx` do not identify the original
surface sampling/rasterization implementation. Therefore the repository does
not claim an exact raw-image reproduction of the paper values in the absence
of those images. Differences in contour rasterization, surface extraction,
resampling, weighting or percentile definition can change distance metrics.

Synthetic tests cover identical/disjoint/shifted masks, anisotropic spacing,
unequal bidirectional surface counts, bilateral union, empty-mask rejection,
labelmap geometry and a simple RTSTRUCT round trip. They contain no patient
images.

API references: [SciPy distance transform](https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html),
[SimpleITK DICOM series reading](https://simpleitk.readthedocs.io/en/master/link_DicomSeriesReader_docs.html),
[rt-utils documentation](https://github.com/qurit/rt-utils).
