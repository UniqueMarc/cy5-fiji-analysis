# Cy5-mRNA confocal image analysis in Fiji

Process Hoechst (405 nm) and Cy5 (640 nm) ND2 Z-stacks, measure exact Cy5 components, and calculate mean punctum area and mean integrated fluorescence for each nucleus-associated analysis region.

## Requirements

- Fiji with ImageJ1, Bio-Formats and Groovy scripting support.
- A graphical desktop session and approximately 4 GB Java heap memory.
- Python 3.9 or later for CSV export; no additional Python packages are required.

Full ND2 stacks are external inputs. Small projection and ROI examples are supplied. The complete measured dataset and ready-to-use per-cell tables are in `results/per_cell/`.

## Run in Fiji

1. Open `Run_in_Fiji.groovy` in Fiji's Script Editor, select **Groovy**, and click **Run**.
2. Select the project folder, input parent folder, output parent folder and `parameters/batch_20260930.json`.
3. Keep the manifest option selected for the supplied dataset. Place each ND2 at the relative path recorded in `data_manifest.json`.
4. Wait for all three stages: image processing, exact-pixel measurement, and punctum-to-territory assignment.

Input layout is `INPUT_PARENT/<Group>/<filename>.nd2`. For new inputs, use matching group folder names and uncheck the manifest option. Inputs must contain one series, one time point, 16-bit unsigned pixels and valid micrometre calibration. Channel metadata must identify 405 and 640 at the configured indices.

## Export per-cell tables

From the project folder, using a new destination folder:

```sh
python scripts/summarize.py "PATH_TO_RUN/Per_cell_measurements/per_cell_measurements.json" "NEW_TABLE_FOLDER"
```

To reproduce the supplied tables without ND2 files:

```sh
python scripts/summarize.py results/per_cell/per_cell_measurements.json reproduced_tables
```

| Output | Contents |
| --- | --- |
| `Per_cell_results.csv` | One row per nuclear territory, including punctum count, means and sums |
| `Mean_punctum_area_um2_wide.csv` | Per-territory mean punctum areas, one group per column |
| `Mean_punctum_integrated_intensity_raw_wide.csv` | Per-territory mean integrated punctum fluorescence |
| `Puncta_assignments.csv` | Every measured punctum, its cell ID, geometry and raw intensity |
| `Group_summary.csv` | Valid/empty territory counts and descriptive summaries of per-cell means |
| `Images.csv` | Image identifiers, groups and population counts |

Wide columns are independent groups; rows across columns are not matched observations. Blank padding and missing punctum means are not zeros. Open the CSV files in Excel or GraphPad Prism.

## Measurement definitions

The nuclear ROI is expanded outward by 20 micrometres, including the nucleus. Overlap belongs to the nearest nuclear contour centroid. These are **nucleus-associated territories**, not segmented cell membranes. Each whole punctum is assigned once, using the rounded centroid of its original detection contour.

Within each territory:

- Mean punctum area = sum of assigned exact component areas / assigned punctum count.
- Mean integrated punctum fluorescence = sum of assigned raw component intensities / assigned punctum count.
- A zero-punctum territory retains count and sums of zero; its punctum means are missing and excluded from the mean-punctum distributions.
- Unassigned puncta retain `Cell_ID = 0` in the detail table and do not enter cell means.

The supplied dataset contains 115 images, 537 nuclear territories and 7,152 puncta: 6,507 assigned and 645 unassigned. There are 505 nonempty territories and 32 empty territories. Mean-punctum measurements retain full precision; an average can be fractional even though raw pixel sums are integers.

No background subtraction is applied to Cy5 raw integrated intensity. Display ranges, 35 × 35 micrometre crops and crop resizing do not affect measurements or inclusion. Crops are not used to select the quantitative population. No significance tests are automatically exported. See [methods](docs/METHODS_AND_LIMITATIONS.md) and [parameters](docs/PARAMETERS.md).

## Image and intermediate outputs

- Raw ZMAX projections, display images, detection masks and ROI ZIP files are written by the first stage.
- `Exact_pixel_measurements/` contains exact component labels and per-punctum measurements.
- `Per_cell_measurements/` contains assignment data and non-overlapping `Territory_labels/` TIFFs.
- `results/final_20260930/frozen_exact_measurements.json` contains the complete per-punctum geometry reference; use `results/per_cell/` for the per-cell figure distributions.

## Check the installation

```sh
python -m unittest discover -s tests -v
```

See [tests and examples](docs/TESTING.md) to reproduce pixel measurements and nuclear-territory assignments with the bundled example.
