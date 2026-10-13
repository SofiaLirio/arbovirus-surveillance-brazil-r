# Full rebuild from raw data: packages -> raw->db ETL -> all analyses (R).
# Usage (from any folder):  powershell -ExecutionPolicy Bypass -File run_all.ps1 [-SkipEtl] [-SkipInstall]
param([switch]$SkipEtl, [switch]$SkipInstall)
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
if (-not $Rscript) {
    $cand = Get-ChildItem "C:\Program Files\R\R-*\bin\Rscript.exe" -ErrorAction SilentlyContinue | Sort-Object FullName -Descending | Select-Object -First 1
    if ($cand) { $Rscript = $cand.FullName } else { throw "Rscript not found: install R 4.5+ from https://cloud.r-project.org" }
}
$args_ = @()
if ($SkipEtl) { $args_ += "--skip-etl" }
if ($SkipInstall) { $args_ += "--skip-install" }
& $Rscript (Join-Path $Root "run_all.R") @args_
if ($LASTEXITCODE) { throw "run_all.R failed (exit $LASTEXITCODE)" }
