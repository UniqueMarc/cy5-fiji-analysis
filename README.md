# Cy5-mRNA confocal image analysis in Fiji

Process Hoechst (405 nm) and Cy5 (640 nm) ND2 Z-stacks, detect Cy5 regions, export nucleus-centred crops, and measure region area and raw integrated fluorescence.

## Requirements

- Fiji with ImageJ1, Bio-Formats and Groovy scripting support.
- A graphical desktop session and approximately 4 GB Java heap memory.
- Python 3.9 or later for CSV export; no additional Python packages are required.

Full ND2 stacks are not included. Small example images and an ROI ZIP are in `examples/`. Supplied CSV tables can be opened directly without running Fiji.

## Run in Fiji

1. Open `Run_in_Fiji.groovy` in Fiji's Script Editor.
2. Select **Groovy** and click **Run**.
3. Select the project folder, input parent folder, output parent folder and `parameters/batch_20260930.json`.
4. Keep the manifest option selected to process the dataset listed in `data_manifest.json`. The input parent must contain the group folders and ND2 filenames specified by `Raw_relative_path`.
5. Wait for the completion message. Results are saved in a new timestamped folder under the output parent.

Input layout:

```text
input_parent/
  <Group>/
    <filename>.nd2
```

For a new dataset, copy the parameter file and edit `group_order`, channel indices and other settings. Use matching group folder names and uncheck the manifest option. ND2 files are then discovered within those folders.

Inputs must contain one series, one time point, 16-bit unsigned pixels and valid micrometre calibration. Channel metadata must identify `405` and `640` at the configured indices. See [parameters](docs/PARAMETERS.md).

## Export tables

From the project folder:

```sh
python scripts/summarize.py "PATH_TO_RUN/Exact_pixel_measurements/Geometry_QC_results.json" "NEW_TABLE_FOLDER"
```

To export tables from the supplied measurements:

```sh
python scripts/summarize.py results/final_20260930/frozen_exact_measurements.json reproduced_tables
```

Choose a new destination folder. Supplied tables are in `results/final_20260930/tables/`.

| Output | Contents |
| --- | --- |
| `Area_um2_wide.csv` | Cy5 region area, one group per column |
| `Integrated_density_raw_wide.csv` | Raw integrated intensity, one group per column |
| `Exact_regions_long.csv` | Individual region measurements |
| `Group_summary.csv` | Group counts and descriptive statistics |
| `Per_image_summary.csv` | Per-image descriptive statistics |
| `Images.csv` | Image-level metadata |
| `Exploratory_Student_tests.csv` | Exploratory pooled-region comparisons |

Wide tables use blank padding. Rows across group columns are not matched observations. Open CSV files in Excel or import them into GraphPad Prism.

## Image outputs

- `405_Hoechst_ZMAX/` and `640_Cy5_ZMAX/`: maximum-intensity projections and display images.
- `Merged_QC/`: full-field merged images.
- `Single_Cell_Crops/`: nucleus-centred crops and outlined views.
- `Cy5_Puncta_Masks/`, `Nucleus_Masks/` and `ROI_Files/`: masks and ROI files.
- `Exact_pixel_measurements/Geometry_QC_results.json`: region measurements used for CSV export.
- `Exact_pixel_measurements/Labels/`: non-overlapping region label images.

## Use the measurements

Area is the accepted component's pixel count multiplied by native pixel area. Integrated density is the sum of unsmoothed raw Cy5 projection intensities, without background subtraction. All accepted full-field Cy5 regions are included, including regions without a nuclear assignment.

The supplied settings display Cy5 in magenta over raw intensities 131-220 and Hoechst in blue with per-image limits. Display ranges and crop resizing do not affect measurements. Crops are 35 x 35 micrometres and exported at 512 x 512 pixels. Scale bars are not added automatically.

Nuclear counts are automatically segmented regions. Region-level Student tests are exploratory; regions within an image are not independent biological replicates. Use experiment-level metadata for replicate-based or hierarchical inference. See [measurement definitions](docs/METHODS_AND_LIMITATIONS.md).

## Check the installation

```sh
python tests/test_release.py
```

See [example and command-line usage](docs/TESTING.md) to run the supplied pixel example or selected ND2 inputs.
