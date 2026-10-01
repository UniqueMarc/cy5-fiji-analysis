param(
    [Parameter(Mandatory=$true)][string]$FijiFolder,
    [Parameter(Mandatory=$true)][string]$JavaExe,
    [Parameter(Mandatory=$true)][string]$InputRoot,
    [Parameter(Mandatory=$true)][string]$OutputParent,
    [ValidateSet('manifest','discover')][string]$Mode='manifest',
    [string]$StableIds=''
)
$ErrorActionPreference='Stop'
$releaseRoot=Split-Path $PSScriptRoot -Parent
$jarRoots=@((Join-Path $FijiFolder 'jars'),(Join-Path $FijiFolder 'plugins'))
$releaseClasspath=((Get-ChildItem -LiteralPath $jarRoots -Filter '*.jar' -File -Recurse |
    Select-Object -ExpandProperty DirectoryName -Unique | ForEach-Object { Join-Path $_ '*' }) -join ';')
$arguments=@('-Xmx4g',"-Dplugins.dir=$(Join-Path $FijiFolder 'plugins')","-Dij.dir=$FijiFolder",'-cp',$releaseClasspath,
    'groovy.ui.GroovyMain',(Join-Path $PSScriptRoot 'run.groovy'),$releaseRoot,$InputRoot,$OutputParent,
    (Join-Path $releaseRoot 'parameters/batch_20260930.json'),$Mode)
if($StableIds){$arguments+=$StableIds}
& $JavaExe @arguments
if($LASTEXITCODE -ne 0){throw "Fiji analysis failed (exit $LASTEXITCODE)"}
