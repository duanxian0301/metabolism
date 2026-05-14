$rootDir = "D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion"
$stepDir = Join-Path $rootDir "step21_native_wsl_factor_gwas_inputs_nonlipid8"
$stdoutLog = Join-Path $stepDir "prepare_nonlipid8_wsl_stdout.log"
$stderrLog = Join-Path $stepDir "prepare_nonlipid8_wsl_stderr.log"
$chunkManifest = Join-Path $stepDir "chunk_manifest.tsv"

Write-Output "**Files**"
Get-ChildItem -Force $stepDir | Select-Object Name, Length, LastWriteTime

Write-Output "`n**WSL R Processes**"
wsl.exe -e bash -lc "ps -ef | grep '11_prepare_wsl_factor_gwas_inputs_nonlipid8.R' | grep -v grep" 2>$null

Write-Output "`n**chunk_manifest exists**"
Test-Path $chunkManifest

if (Test-Path $stdoutLog) {
  Write-Output "`n**stdout tail**"
  Get-Content $stdoutLog -Tail 40
}

if (Test-Path $stderrLog) {
  Write-Output "`n**stderr tail**"
  Get-Content $stderrLog -Tail 40
}
