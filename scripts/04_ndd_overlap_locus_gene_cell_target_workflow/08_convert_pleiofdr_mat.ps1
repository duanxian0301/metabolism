$ErrorActionPreference = "Stop"

$PleioRoot = "D:\pleioFDR\pleiofdr-master"
$Python = "python"

Push-Location $PleioRoot
try {
  & $Python "python_convert\sumstats.py" mat --sumstats "data\lipid8_F2_fdr.txt" --ref "9545380.ref" --out "data\lipid8_F2_fdr.mat" --force
  & $Python "python_convert\sumstats.py" mat --sumstats "data\AD_metabolic_fdr.txt" --ref "9545380.ref" --out "data\AD_metabolic_fdr.mat" --force
}
finally {
  Pop-Location
}
