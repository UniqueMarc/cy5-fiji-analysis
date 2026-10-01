param(
    [Parameter(Mandatory=$true)][string]$FijiFolder,
    [Parameter(Mandatory=$true)][string]$JavaExe,
    [Parameter(Mandatory=$true)][string]$NewOutputFolder
)
$ErrorActionPreference='Stop'
$releaseRoot=Split-Path $PSScriptRoot -Parent
$jarRoots=@((Join-Path $FijiFolder 'jars'),(Join-Path $FijiFolder 'plugins'))
$releaseClasspath=((Get-ChildItem -LiteralPath $jarRoots -Filter '*.jar' -File -Recurse |
    Select-Object -ExpandProperty DirectoryName -Unique | ForEach-Object { Join-Path $_ '*' }) -join ';')
& $JavaExe '-Xmx2g' '-cp' $releaseClasspath 'groovy.ui.GroovyMain' (Join-Path $PSScriptRoot 'run_pixel_example.groovy') $releaseRoot $NewOutputFolder
if($LASTEXITCODE -ne 0){throw "Pixel example failed (exit $LASTEXITCODE)"}
