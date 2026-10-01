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

## Assignment and crops

Nuclear ROIs are expanded outward by 20 µm. Overlapping candidate territories
are assigned to the nearest nuclear centroid; nuclear interiors are included.
These assignments feed only auxiliary per-cell context. The final quantitative
population is **all accepted full-field Cy5 regions**, independent of
assignment. Bright puncta elsewhere in a displayed crop can belong to another
cell or be extracellular; no membrane was segmented.

Eligible nucleus-centred crops are nominally 35 × 35 µm. Native crop dimensions
are rounded from calibration, so actual extent differs by at most approximately
half a native pixel per dimension. Crops crossing the field edge are skipped;
the minimum neighbouring nuclear-centroid distance is 20 µm. This does not
guarantee that a second nucleus is absent or that a merged nucleus was split.
The crop is bicubically resampled to 512 × 512 for display only. Cyan outlines
mark the target nucleus; green outlines mark accepted regions intersecting
the crop. No number labels, yellow territory ring or scale bars are added.

Hoechst display uses per-image 0.175th/99.825th-percentile clipping; Cy5 uses
fixed 131–220 raw limits. RGB channels are blue for Hoechst and magenta for
Cy5. Changing a display range alone cannot change the quantitative outputs.

## Statistical limitations

Student tests are independent, two-sided, equal-variance tests on pooled region
measurements, provided for exploratory comparisons.
Holm corrections are supplied per metric (10 planned comparisons) and across
both metrics (20). Raw and adjusted stars must not be mixed. Thresholds are
P ≤ 0.05 (*), ≤ 0.01 (**), ≤ 0.001 (***), ≤ 0.0001 (****); otherwise ns.

Regions/cells/fields from the same independent experiment are nested. Neither
the number of regions nor automatically segmented nuclei is the biological
replicate n. This dataset manifest does not yet establish independent
experiment identities or paired designs. A replicate-level or hierarchical
analysis requires that mapping. Region-level P values alone cannot establish
biological reproducibility. Representative crops are illustrations, not the
data selection rule for the distribution or tests.
