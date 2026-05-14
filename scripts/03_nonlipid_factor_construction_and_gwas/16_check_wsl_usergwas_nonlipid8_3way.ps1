$rootDir = "D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion"
$outputDir = Join-Path $rootDir "step22_native_wsl_usergwas_nonlipid8_results"

Write-Output "**WSL userGWAS processes**"
wsl.exe -e bash -lc "ps -ef | grep '12_run_wsl_usergwas_nonlipid8.R' | grep -v grep" 2>$null

Write-Output "`n**Current output files**"
Get-ChildItem -Recurse -Force $outputDir | Select-Object FullName, Length, LastWriteTime

Write-Output "`n**Recent logs**"
Get-ChildItem $outputDir -Filter "*.log" -ErrorAction SilentlyContinue | ForEach-Object {
  Write-Output ("`n== {0} ==" -f $_.Name)
  Get-Content $_.FullName -Tail 30
}
