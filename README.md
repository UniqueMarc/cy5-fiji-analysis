# Cy5-mRNA confocal image analysis in Fiji

Local release candidate **1.0.0**, prepared 2026-10-01 from the pipeline used
for the 20260930 U2OS dataset. Public upload has **not** been performed.

The pipeline projects Hoechst (405 nm) and Cy5 (640 nm) ND2 Z-stacks, detects
Cy5-positive connected regions, exports nucleus-centred crops, and measures
area and raw integrated fluorescence on exact foreground pixels.

**Use the exact-pixel results, not the older polygon/per-cell tables.**
The final dataset contains 115 included fields, 7,152 Cy5 regions and 537
automatically segmented nuclear regions. These nuclear regions are not
independently verified cell counts. The release contains only the valid
dataset used for the final analysis.

## Contents

| Item | Purpose |
| --- | --- |
| `Run_in_Fiji.groovy` | Graphical entry point; runs detection and exact-pixel measurement |
| `scripts/` | Portable Fiji pipeline, CSV export and statistical functions |
| `parameters/batch_20260930.json` | Actual settings for the latest dataset |
| `data_manifest.json` / `.csv` | Stable image IDs, current names/groups and raw SHA-256 for the valid dataset |
| `results/final_20260930/` | Frozen final measurements, wide tables and descriptive results |
| `examples/` | Small real raw-ZMAX/mask example with expected corrected labels |
| `tests/` | Numerical regression tests and the pixel-example test |
| `docs/` | Methods, parameters, validation and submission checklist |
| `provenance/` | Hashes of the historical sources used to prepare this release |

## Requirements

- Fiji with ImageJ1, Bio-Formats and Groovy scripting support.
- Tested component versions are recorded in `docs/VALIDATION.md`.
- Python 3.9 or later for CSV summaries/tests; only the standard library is
  required. Python is not needed for the Fiji image-processing steps.
- A graphical environment is required by the legacy ImageJ ROI Manager;
  this release does not claim support for pure headless Fiji.
- Start with 4 GB Java heap. Raw ND2 files are external and are not bundled.

Do not update Fiji between reproductions without checking the supplied tests.
Fiji/Bio-Formats/Java are dependencies, not redistributed in this repository.

## Run in Fiji

1. Open `Run_in_Fiji.groovy` in Fiji's Script Editor and select **Groovy**.
2. Click **Run**. Select this release folder, the input parent folder, an
   output parent, and `parameters/batch_20260930.json`.
3. Keep the manifest option checked to reproduce the final 20260930 data.
   The input parent must contain the eight current group folders and the
   included ND2 filenames listed in `data_manifest.json`.
4. Results go into a new timestamped folder; input files are never renamed
   or overwritten. Read failures abort completion after recording status.

For a **new dataset**, copy/edit the configuration, set its group folder
names and channel indices, and uncheck the final-manifest option. Files are
then discovered recursively inside those group folders. Do not reuse the
20260930 image manifest for a new experiment. The configuration's display
range is batch-specific; it is not a recommended universal exposure range.

The input contract is one series, one time point, uint16 images, valid
micrometre calibration, and channel names `405` and `640` at the configured
indices. Multi-timepoint/multi-series files are rejected, not silently
interpreted as extra Z sections.

## Outputs and final tables

- `405_Hoechst_ZMAX/`, `640_Cy5_ZMAX/`: raw calibrated TIFFs and display PNGs.
- `Merged_QC/`: full-field merged/QC images.
- `Single_Cell_Crops/`: 405, 640, unmarked merged, and green-outline QC crops.
- `Cy5_Puncta_Masks/`, `Nucleus_Masks/`, `ROI_Files/`: detection audit data.
- `Exact_pixel_measurements/Geometry_QC_results.json`: **final measurements**.
- `Exact_pixel_measurements/Labels/`: canonical, non-overlapping pixel labels.
- `Legacy_NOT_for_final_quantitation/`: provisional historical per-cell and
  polygon summaries; these are not the source of the final violin tables.

To produce Excel/Prism-readable CSVs, use:

```text
python scripts/summarize.py "PATH_TO_RUN/Exact_pixel_measurements/Geometry_QC_results.json" "NEW_TABLE_FOLDER"
```

For the bundled final data:

```text
python scripts/summarize.py results/final_20260930/frozen_exact_measurements.json reproduced_tables
python tests/test_release.py
```

`Area_um2_wide.csv` and `Integrated_density_raw_wide.csv` contain **one group
per column** with blank padding, not zero padding. Row positions across
groups do not represent matched observations. The already-generated copies
are under `results/final_20260930/tables/`.

## Measurement and display are separate

For this dataset, all Cy5 display PNGs use **131–220 raw intensity units**,
linearly mapped to magenta with clipping. Hoechst uses per-image display
limits and blue colour. Display pixels, enlarged crops, and RGB PNGs are
never used for measurement. No scale bar is drawn automatically.

Detection uses raw Cy5 ZMAX intensities to set an image-specific robust
threshold, then thresholds a Gaussian-smoothed detection copy. Fixed display
limits do **not** imply a fixed detection threshold. The final area is the
number of pixels in each accepted 8-connected component times native XY
pixel area. Integrated density is the sum of **unsmoothed, unscaled raw
ZMAX** pixel values in that component, without background subtraction.

The latest release measures **all accepted regions in the full field**,
including unassigned/extracellular candidates. Perinuclear assignment and
35 × 35 µm crops are supplementary views, not the population used for the
final distributions. Green QC outlines show all accepted intersecting
regions, not every coloured pixel. See `docs/METHODS_AND_LIMITATIONS.md`.

## Statistical scope

The bundled two-sided, equal-variance Student tests reproduce the historical
**pooled-region exploratory comparisons**. Raw and Holm-adjusted P values
are separate. Puncta within a cell/image are not independent biological
replicates; these tests must not be presented as tests of three independent
experiments. The authoritative image-to-biological-replicate mapping is not
available in this release, so that field is blank rather than invented.
Per-image summaries are provided for a later appropriately nested analysis.
The manually supplied three-row normalized tables are not mixed into this
image-derived dataset.

## Raw-data archive and publication status

Full ND2 stacks are not included. `data_manifest.csv` identifies the 115 valid
raw files by stable ID and SHA-256. Paths are relative to the raw input parent.
The files should later be deposited in a suitable data archive; add its
verified accession/DOI here before publication. No accession is fabricated.

Before public upload, confirm author order/contact, code/data licenses,
biological-replicate metadata and the data deposit. See
`docs/PUBLICATION_CHECKLIST.md`. No open-source license has yet been selected;
`LICENSE_PENDING.txt` is intentionally not a license grant.

Chinese quick start: `使用说明.md`.
