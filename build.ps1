$ErrorActionPreference = 'Stop'
$paperRoot = [System.IO.Path]::GetFullPath($PSScriptRoot)
$buildRoot = Join-Path $paperRoot '.build'
New-Item -ItemType Directory -Path $buildRoot -Force | Out-Null
$env:TECTONIC_CACHE_DIR = Join-Path $buildRoot 'tectonic-cache'
$env:TEMP = $buildRoot
$env:TMP = $buildRoot
$tectonicCommand = Get-Command tectonic -ErrorAction SilentlyContinue
if (-not $tectonicCommand) {
    throw 'Install Tectonic, add it to PATH, and rerun this script; or run pdflatex main.tex twice.'
}
Push-Location -LiteralPath $paperRoot
try {
    & $tectonicCommand.Source --keep-logs --outdir $paperRoot (Join-Path $paperRoot 'main.tex')
    if ($LASTEXITCODE -ne 0) { throw 'TeX compilation failed.' }
} finally {
    Pop-Location
}
