# Tests and examples

## Check the supplied per-cell tables

```sh
python -m unittest discover -s tests -v
```

Checks cover all 537 per-cell records against an independent expected table, 7,152 puncta, 115 image IDs, raw intensity conservation, zero-punctum handling, invalid/duplicate assignments, and every exported wide-table value.

## Run the pixel and assignment example

The example includes a raw Cy5 projection, detection mask, punctum ROIs, nuclear ROIs and expanded territory ROIs. It verifies exact connected-component measurements, label pixels, punctum ownership and territory pixel counts.

```powershell
.\tests\Run_PixelExample.ps1 -FijiFolder "PATH_TO_FIJI" -JavaExe "PATH_TO_JAVA_EXE" -NewOutputFolder "NEW_EXAMPLE_OUTPUT"
```

## Process selected ND2 inputs

Arrange raw files according to `data_manifest.json`:

```powershell
.\scripts\Run_CLI.ps1 -FijiFolder "PATH_TO_FIJI" -JavaExe "PATH_TO_JAVA_EXE" -InputRoot "INPUT_PARENT" -OutputParent "OUTPUT_PARENT" -StableIds "4,17"
```

Omit `-StableIds` to process all included inputs. A complete run produces `Per_cell_measurements/per_cell_measurements.json`, which is the input to `scripts/summarize.py`. All commands require new output folders and refuse to overwrite existing analyses.
