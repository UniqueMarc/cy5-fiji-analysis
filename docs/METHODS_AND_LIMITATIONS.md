# Methods and interpretation

## Processing order

1. Import ND2 with Bio-Formats; verify channel identity, dimensions, type and
   spatial calibration. The recorded preset uses channel **3 = 405**,
   **4 = 640** (one-based), with other channels ignored.
2. Independently take the maximum across every Z section in each channel.
   These are 2-D maximum-intensity projections, not 3-D object measurements.
3. Segment nuclei on a separate Hoechst projection: Gaussian sigma 0.5 µm,
   Otsu threshold × 1.0, binary close, fill holes, no watershed, no padding;
   accept 30–1000 µm², circularity 0.1–1.0, excluding edge nuclei.
4. Compute Cy5 raw-ZMAX median and median absolute deviation (MAD).
   Automatic threshold = median + 4 × 1.4826 × max(1, MAD). A positive
   configured manual threshold overrides it. Do not clamp a threshold above
   the image maximum down to the maximum, which would force a false positive.
5. Gaussian-smooth a **detection copy** by 0.6 native pixels, threshold it,
   and accept ImageJ connected components using 0.02–400 µm² and no peak
   splitting. Edge Cy5 regions are not excluded. Area acceptance is applied
   during particle analysis.
6. Flood-fill each accepted component from the ImageJ particle seed using
   8-connectivity. Validate its size against ImageJ's internal pixel count.
   Exclude internal holes and disconnected internal islands from that
   component. Verify no pixel is counted in more than one final label.
7. Measure raw unsmoothed ZMAX on those native pixels; export label TIFFs and
   composite ROI vectors. Label TIFFs are authoritative: rerasterizing an ROI
   vector can differ at thin boundaries.

## Measurements

| Field | Definition |
| --- | --- |
| `Area_pixels` | Exact accepted connected-component pixel count |
| `Area_um2` | `Area_pixels × pixelWidth_um × pixelHeight_um` |
| `Integrated_density_raw` | Sum of unscaled raw Cy5 ZMAX pixel intensities |
| `Mean_raw_intensity` | Raw sum divided by exact component pixel count |
| `Max_raw_intensity` | Maximum raw value in that component |
| `X_pixel`, `Y_pixel` | Mean native pixel coordinates of component members |
| `Unique_ID` | Stable field ID plus component scan-order ID |

Raw intensity units are detector units, not calibrated photon/molecule counts.
No background subtraction, size trimming, normalization or post hoc outlier
deletion is applied to these final measurements. Cy5 objects can contain
aggregates and overlapping spots; the no-split configuration does not resolve
individual molecules. Area reflects projected segmentation, not a particle's
physical diameter or 3-D volume.

## Assignment and per-cell aggregation

Nuclear ROIs are expanded outward by 20 micrometres. Nuclear interiors are included. Overlap belongs to the nearest nuclear contour centroid; equal distances retain the earlier nuclear ROI ID. Each territory has a stable cell ID within its image.

Assign each entire punctum at its original detection-contour centroid, rounded using Java Math.round. The exact component pixel centroid is a separate geometry measurement and is not used for assignment. No punctum is assigned twice. Unassigned puncta have Cell_ID 0, remain in the detailed output, and are excluded from per-cell aggregation.

For each territory, sum the exact component areas and raw integrated intensities of its assigned puncta. Divide each sum by the assigned punctum count to obtain mean punctum area and mean integrated punctum fluorescence. Each nonempty territory contributes one value to each distribution. Group summaries weight territories equally, regardless of punctum count. An intensity sum over all puncta divided by all puncta in a group is a different quantity and is not the reported per-cell mean.

Zero-punctum territories retain zero counts and sums. Their punctum means are null in JSON and blank in CSV; they do not contribute zeros to mean-punctum distributions. Numerical means retain full precision.

## Display crops

Nucleus-centred crops are nominally 35 × 35 micrometres and resampled to 512 × 512 for display. Crops crossing the field edge are skipped; minimum neighbouring nuclear-centroid distance is 20 micrometres. Crop eligibility does not determine inclusion in per-cell quantification. All nuclear territories are retained.

Hoechst display uses per-image percentile limits and Cy5 uses fixed raw limits 131–220. Display adjustment and resizing do not change measured values. These territories are not membrane segmentations and cannot establish that every assigned spot is intracellular.

## Interpretation

One value per nuclear territory describes the cell-associated sampling unit. Cells and fields from the same experiment remain nested; cell counts do not establish the number of independent biological experiments. This exporter provides descriptive values and does not generate significance annotations or assume paired observations from row order.
