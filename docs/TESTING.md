# Run checks and examples

## Check the supplied tables

From the project folder:

```sh
python tests/test_release.py
```

## Run the pixel example

The example uses the bundled Cy5 projection, mask and ROI files. In Windows PowerShell, provide your Fiji folder, Java executable and a new output folder:

```powershell
.\tests\Run_PixelExample.ps1 -FijiFolder "PATH_TO_FIJI" -JavaExe "PATH_TO_JAVA_EXE" -NewOutputFolder "NEW_EXAMPLE_OUTPUT"
```

## Process selected ND2 inputs

Provide the external ND2 dataset arranged according to `data_manifest.json`:

```powershell
.\scripts\Run_CLI.ps1 -FijiFolder "PATH_TO_FIJI" -JavaExe "PATH_TO_JAVA_EXE" -InputRoot "INPUT_PARENT" -OutputParent "OUTPUT_PARENT" -StableIds "4,17"
```

Omit `-StableIds` to process every included manifest entry. Run Fiji with a graphical desktop session.
