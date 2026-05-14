$ErrorActionPreference = "Stop"

$Script = "/mnt/d/codex/GenomicSEM/scripts/04_postgwas/21_run_mctwas_single_trait.R"
$OutputRoot = "/home/shenjing/ctwas_paths/metabolic_mctwas_work_lipid8_F1_PD"
$RefRoot = "/home/shenjing/ctwas_paths/mctwas_work"

wsl -e bash -lc "mkdir -p $OutputRoot && cp -r $RefRoot/manifests $OutputRoot/ && cp -r $RefRoot/refs $OutputRoot/"

wsl -e bash -lc "Rscript $Script --trait lipid8_F1 --gwas-file /mnt/d/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/work/clean_factor_inputs/lipid8_F1_clean.txt --output-root $OutputRoot --ncore 4 --niter-prefit 3 --niter 30 --L 5 --seed 20260429"
wsl -e bash -lc "Rscript $Script --trait PD --gwas-file /mnt/d/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/work/ctwas_inputs/PD.txt --output-root $OutputRoot --ncore 4 --niter-prefit 3 --niter 30 --L 5 --seed 20260429"
