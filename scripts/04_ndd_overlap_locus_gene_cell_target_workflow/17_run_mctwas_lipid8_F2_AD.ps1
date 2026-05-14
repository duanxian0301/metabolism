$ErrorActionPreference = "Stop"

$Script = "/mnt/d/codex/GenomicSEM/scripts/04_postgwas/21_run_mctwas_single_trait.R"
$OutputRoot = "/home/shenjing/ctwas_paths/metabolic_mctwas_work"
$RefRoot = "/home/shenjing/ctwas_paths/mctwas_work"

# Reuse the reference manifests and LD matrices prepared for the ALPS manuscript.
wsl -e bash -lc "mkdir -p $OutputRoot && cp -r $RefRoot/manifests $OutputRoot/ && cp -r $RefRoot/refs $OutputRoot/"

wsl -e bash -lc "Rscript $Script --trait lipid8_F2 --gwas-file /mnt/d/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/work/clean_factor_inputs/lipid8_F2_clean.txt --output-root $OutputRoot --ncore 4 --niter-prefit 3 --niter 30 --L 5 --seed 20260421"
wsl -e bash -lc "Rscript $Script --trait AD --gwas-file '/mnt/d/文章/4NDD/NDDGWAS/AD.txt' --output-root $OutputRoot --ncore 4 --niter-prefit 3 --niter 30 --L 5 --seed 20260421"
