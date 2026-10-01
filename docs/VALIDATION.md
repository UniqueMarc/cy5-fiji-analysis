# Local validation on 2026-10-01

## Environment

Windows 10; Java 1.8.0_322; ImageJ 1.53t; Bio-Formats 6.11.1;
Groovy 3.0.4. Python regression tests use only the standard library.
Equivalent behaviour on newer Fiji versions or other operating systems has
not been established.

## Tests actually performed

| Check | Result |
| --- | --- |
| ND2 end-to-end reprocessing | Stable image IDs 4 and 17; 76 total Cy5 regions, 11 nuclear regions, 3 crop sets |
| Raw 405/640 maximum projections | Both fields pixel-identical to archived projections |
| Cy5 and nuclear detection masks | Both fields pixel-identical to archived masks |
| Final label TIFFs | Both fields pixel-identical to archived exact-pixel labels |
| Final region measurements | All 76 IDs, groups, areas, sums, means, maxima and centroids match |
| Fixed Cy5 display | Full-field 131–220 PNG pixels match the final display revision |
| Packaged small example | 58 regions match expected values and label pixels; 2 geometries corrected |
| Final table regression | 5 tests pass: population, means/medians/SD, P values/stars, geometry definitions, group-wide columns |
| Raw integrity inventory | The 115 valid ND2 files have verified SHA-256 values; previously recorded moved-file hashes agree |
| Example display inspection | Blue Hoechst, magenta Cy5, cyan target nucleus, green detected-region outlines; no scale bar |

This validation did **not** reprocess all 115 included raw images. The full
final tables were reproduced from the frozen corrected measurements, while
two raw images were reprocessed end-to-end. The GUI dialog was not operated
in this validation; its shared processing entry point was exercised via CLI.
`validation_results.json` records the tested IDs and component versions.

## Re-run tests

From this release folder:

```text
python tests/test_release.py
```

Small pixel-correction example (Windows PowerShell; supply your Fiji and Java
locations and a new destination):

```powershell
.\tests\Run_PixelExample.ps1 -FijiFolder "PATH_TO_FIJI" -JavaExe "PATH_TO_JAVA_EXE" -NewOutputFolder "NEW_EXAMPLE_OUTPUT"
```

End-to-end smoke test requires the external included ND2 dataset:

```powershell
.\scripts\Run_CLI.ps1 -FijiFolder "PATH_TO_FIJI" -JavaExe "PATH_TO_JAVA_EXE" -InputRoot "INCLUDED_RAW_ROOT" -OutputParent "OUTPUT_PARENT" -StableIds "4,17"
```

Omit `-StableIds` to process every included manifest file. Raw images are
never renamed or modified by these commands.

## Packaging changes, not scientific re-tuning

The executed historical wrapper was expanded into a portable source file;
local drive paths and historical label-reassignment plumbing were replaced
by explicit configuration and the final manifest. Final display limits are
applied directly. All accepted ROI vectors are now exported in the first
stage so exact-pixel measurement runs automatically rather than relying on a
separate historical ROI-export step. Provisional per-cell/polygon tables
were moved under a clearly labelled legacy output directory. An unused
hand-built XLSX writer was removed; canonical tables use portable CSV.

The detection/acceptance algorithms were not retuned to obtain a different
group ranking. Source hashes and frozen final results are included for
traceability. This is a local release candidate, not yet a public archive.
