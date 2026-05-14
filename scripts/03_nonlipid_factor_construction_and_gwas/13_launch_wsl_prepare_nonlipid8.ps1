$rootDir = "D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion"
$workDir = "D:\codex\GenomicSEM\metabolic"
$logDir = Join-Path $rootDir "step21_native_wsl_factor_gwas_inputs_nonlipid8"
$stdoutLog = Join-Path $logDir "prepare_nonlipid8_wsl_stdout.log"
$stderrLog = Join-Path $logDir "prepare_nonlipid8_wsl_stderr.log"

New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$wslCommand = @"
cd /mnt/d/codex/GenomicSEM/metabolic &&
/usr/bin/Rscript /mnt/d/codex/GenomicSEM/metabolic/nonlipid_module/scripts/11_prepare_wsl_factor_gwas_inputs_nonlipid8.R
"@

$proc = Start-Process `
  -FilePath "wsl.exe" `
  -ArgumentList @("bash", "-lc", $wslCommand) `
  -WorkingDirectory $workDir `
  -RedirectStandardOutput $stdoutLog `
  -RedirectStandardError $stderrLog `
  -PassThru

Write-Output ("Started WSL prepare_nonlipid8 job. PID={0}" -f $proc.Id)
Write-Output ("stdout: {0}" -f $stdoutLog)
Write-Output ("stderr: {0}" -f $stderrLog)
