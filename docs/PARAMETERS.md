# Parameters for the latest dataset

Edit `parameters/batch_20260930.json`; retain a separate copy for new batches.
The complete configuration is saved with every run.

| Parameter | Recorded setting | Effect |
| --- | --- | --- |
| `hoechstChannel`, `cy5Channel` | 3, 4 | One-based channels; metadata must say 405 and 640 |
| `fixedCy5DisplayMinRaw`, `fixedCy5DisplayMaxRaw` | 131, 220 | Display only; shared across the batch |
| `nucleusThresholdMethod`, `nucleusThresholdScale` | Otsu, 1.0 | Nuclear mask threshold |
| `nucleusSigmaUm` | 0.5 | Nuclear detection smoothing |
| `nucleusBackgroundUm`, `nucleusPaddingUm` | 0, 0 | Disabled subtraction/padding |
| `minNucleusAreaUm2`, `maxNucleusAreaUm2` | 30, 1000 | Accepted nuclear ROI size |
| `separateTouchingNuclei` | false | No nuclear watershed |
| `cy5ThresholdK`, `cy5ManualThreshold` | 4, 0 | Robust auto threshold; manual override off |
| `cy5DetectionSigmaPixels` | 0.6 | Smoothing on detection copy only |
| `minPunctumAreaUm2`, `maxPunctumAreaUm2` | 0.02, 400 | Particle acceptance in native calibrated image |
| `splitConnectedCy5` | false | Connected objects remain connected; true is rejected |
| `cy5PeakTolerance` | 12 | Inactive because peak splitting is off |
| `expansionUm`, `excludeNuclearInterior` | 20, false | Auxiliary assignment only; not final field population |
| `singleCellCropUm`, `singleCellOutputPixels` | 35, 512 | Physical crop size and display resampling |
| `minimumNeighborDistanceUm` | 20 | Crop eligibility by nuclear centroid distance |
| `includeSubfolders` | true | Recursive discovery within configured group folders |

Manifest mode reads the explicitly listed included files rather than applying
recursive discovery. Use a new config and discover mode for new experiments.
Recorded display limits came from an initial mRNA reference (131–1125), then
a display-only batch revision to 131–220. The current package uses the final
display limits directly; detection and exact raw measurements are unchanged.
