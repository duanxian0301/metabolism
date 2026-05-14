param(
    [Parameter(Mandatory = $true)]
    [string]$ConfigPath
)

$ErrorActionPreference = "Stop"

$PleioRoot = "D:\pleioFDR\pleiofdr-master"
$ConfigTarget = Join-Path $PleioRoot "config.txt"
$ConfigName = [System.IO.Path]::GetFileNameWithoutExtension($ConfigPath)
$LogDir = "D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd\results\logs_long_tasks"
New-Item -ItemType Directory -Force $LogDir | Out-Null
$StdoutLog = Join-Path $LogDir ($ConfigName + ".stdout.log")
$StderrLog = Join-Path $LogDir ($ConfigName + ".stderr.log")

Copy-Item -LiteralPath $ConfigPath -Destination $ConfigTarget -Force

$MatlabCmd = "try, run('runme.m'); catch ME, disp(getReport(ME)); exit(1); end; exit(0);"

Push-Location $PleioRoot
try {
    & "C:\matlab\R2024b\bin\matlab.exe" -nodisplay -nosplash -nodesktop -r $MatlabCmd 1>> $StdoutLog 2>> $StderrLog
}
finally {
    Pop-Location
}
