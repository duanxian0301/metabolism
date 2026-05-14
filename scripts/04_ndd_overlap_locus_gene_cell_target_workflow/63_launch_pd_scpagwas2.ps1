$ErrorActionPreference = "Stop"

$root = "D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd"
$logDir = Join-Path $root "results\logs_long_tasks"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$pd20k = "D:\scPagwas\MSSM_PD_20k.rds"
$rscript = "D:\R\R-4.4.0\bin\Rscript.exe"

$run1 = Join-Path $root "scripts\61_run_scpagwas2_nonlipid8_F1_MSSM_PD.R"
$run2 = Join-Path $root "scripts\62_run_scpagwas2_lipid8_F1_MSSM_PD.R"

$out1 = Join-Path $logDir "scpagwas2_nonlipid8_F1_MSSM_PD.stdout.log"
$err1 = Join-Path $logDir "scpagwas2_nonlipid8_F1_MSSM_PD.stderr.log"
$out2 = Join-Path $logDir "scpagwas2_lipid8_F1_MSSM_PD.stdout.log"
$err2 = Join-Path $logDir "scpagwas2_lipid8_F1_MSSM_PD.stderr.log"
$launcherLog = Join-Path $logDir "scpagwas2_pd_launcher.log"

function Write-LauncherLog($msg) {
  $ts = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
  Add-Content -Path $launcherLog -Value "[$ts] $msg"
}

Write-LauncherLog "Launcher started."

while (-not (Test-Path $pd20k)) {
  Write-LauncherLog "Waiting for MSSM_PD_20k.rds ..."
  Start-Sleep -Seconds 60
}

Write-LauncherLog "Detected MSSM_PD_20k.rds. Starting nonlipid8_F1 × PD scPagwas2."
$p1 = Start-Process -FilePath $rscript -ArgumentList $run1 -RedirectStandardOutput $out1 -RedirectStandardError $err1 -WindowStyle Hidden -PassThru
$p1.WaitForExit()
Write-LauncherLog "nonlipid8_F1 × PD finished with exit code $($p1.ExitCode)."

Write-LauncherLog "Starting lipid8_F1 × PD scPagwas2."
$p2 = Start-Process -FilePath $rscript -ArgumentList $run2 -RedirectStandardOutput $out2 -RedirectStandardError $err2 -WindowStyle Hidden -PassThru
$p2.WaitForExit()
Write-LauncherLog "lipid8_F1 × PD finished with exit code $($p2.ExitCode)."

Write-LauncherLog "Launcher completed."
