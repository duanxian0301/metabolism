$rootDir = "D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion"
$inputDir = Join-Path $rootDir "step21_native_wsl_factor_gwas_inputs_nonlipid8"
$outputDir = Join-Path $rootDir "step22_native_wsl_usergwas_nonlipid8_results"
$manifestPath = Join-Path $inputDir "chunk_manifest.tsv"
$workDir = "D:\codex\GenomicSEM\metabolic"

if (-not (Test-Path $manifestPath)) {
  throw "Missing chunk manifest: $manifestPath"
}

$manifest = Import-Csv -Delimiter "`t" $manifestPath
$n = $manifest.Count
if ($n -lt 1) {
  throw "chunk_manifest.tsv is empty: $manifestPath"
}

$chunk1End = [int][Math]::Ceiling($n / 3.0)
$chunk2End = [int][Math]::Ceiling(2 * $n / 3.0)

$groups = @(
  @{ Name = "grp1"; Start = 1; End = $chunk1End; Cores = 8 },
  @{ Name = "grp2"; Start = ($chunk1End + 1); End = $chunk2End; Cores = 8 },
  @{ Name = "grp3"; Start = ($chunk2End + 1); End = $n; Cores = 8 }
)

New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

foreach ($g in $groups) {
  if ($g.Start -gt $g.End) {
    continue
  }

  $stdoutLog = Join-Path $outputDir ("usergwas_nonlipid8_{0}_stdout.log" -f $g.Name)
  $stderrLog = Join-Path $outputDir ("usergwas_nonlipid8_{0}_stderr.log" -f $g.Name)

  $wslCommand = @"
cd /mnt/d/codex/GenomicSEM/metabolic &&
START_CHUNK=$($g.Start) END_CHUNK=$($g.End) USE_CORES=$($g.Cores) /usr/bin/Rscript /mnt/d/codex/GenomicSEM/metabolic/nonlipid_module/scripts/12_run_wsl_usergwas_nonlipid8.R
"@

  $proc = Start-Process `
    -FilePath "wsl.exe" `
    -ArgumentList @("bash", "-lc", $wslCommand) `
    -WorkingDirectory $workDir `
    -RedirectStandardOutput $stdoutLog `
    -RedirectStandardError $stderrLog `
    -PassThru

  Write-Output ("Started {0}: PID={1}, chunks {2}-{3}, cores={4}" -f $g.Name, $proc.Id, $g.Start, $g.End, $g.Cores)
  Write-Output ("stdout: {0}" -f $stdoutLog)
  Write-Output ("stderr: {0}" -f $stderrLog)
}
