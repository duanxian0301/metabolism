$ErrorActionPreference = "Stop"

$ProjectRoot = "D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd"
$MixerRoot = Join-Path $ProjectRoot "work\04_mixer_inputs"
$RunDir = Join-Path $MixerRoot "runs\rep1"
$RefDir = "D:\mixer\mixer\reference\1000G_EUR_Phase3_plink"
$Image = "ghcr.io/precimed/gsa-mixer:2.2.1"
$Threads = "12"
$Rep = "1"

New-Item -ItemType Directory -Force $RunDir | Out-Null

function Invoke-Mixer {
    param(
        [Parameter(Mandatory=$true)]
        [string[]]$ArgsList
    )

    docker run --rm `
        -v "${MixerRoot}:/work" `
        -v "${RefDir}:/ref" `
        -w "/work/runs/rep1" `
        $Image `
        python /tools/mixer/precimed/mixer.py @ArgsList
}

function Run-Fit1 {
    param([string]$Trait)
    $JsonPath = Join-Path $RunDir "$Trait.fit.rep$Rep.json"
    if (Test-Path $JsonPath) {
        Write-Host "[SKIP] fit1 $Trait"
        return
    }
    Write-Host "[RUN] fit1 $Trait"
    Invoke-Mixer @(
        "fit1",
        "--trait1-file", "/work/inputs/$Trait.sumstats.gz",
        "--out", "$Trait.fit.rep$Rep",
        "--extract", "/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$Rep.snps",
        "--bim-file", "/ref/1000G.EUR.QC.@.bim",
        "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
        "--threads", $Threads,
        "--exclude-ranges", "MHC"
    )
}

function Run-Fit2Test2 {
    param([string]$Trait1, [string]$Trait2)
    $Prefix = "${Trait1}_vs_${Trait2}"
    $FitJson = Join-Path $RunDir "$Prefix.fit.rep$Rep.json"
    $TestJson = Join-Path $RunDir "$Prefix.test.rep$Rep.json"

    if (-not (Test-Path $FitJson)) {
        Write-Host "[RUN] fit2 $Prefix"
        Invoke-Mixer @(
            "fit2",
            "--trait1-file", "/work/inputs/$Trait1.sumstats.gz",
            "--trait2-file", "/work/inputs/$Trait2.sumstats.gz",
            "--trait1-params", "$Trait1.fit.rep$Rep.json",
            "--trait2-params", "$Trait2.fit.rep$Rep.json",
            "--out", "$Prefix.fit.rep$Rep",
            "--extract", "/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$Rep.snps",
            "--bim-file", "/ref/1000G.EUR.QC.@.bim",
            "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
            "--threads", $Threads,
            "--exclude-ranges", "MHC"
        )
    } else {
        Write-Host "[SKIP] fit2 $Prefix"
    }

    if (-not (Test-Path $TestJson)) {
        Write-Host "[RUN] test2 $Prefix"
        Invoke-Mixer @(
            "test2",
            "--trait1-file", "/work/inputs/$Trait1.sumstats.gz",
            "--trait2-file", "/work/inputs/$Trait2.sumstats.gz",
            "--load-params", "$Prefix.fit.rep$Rep.json",
            "--out", "$Prefix.test.rep$Rep",
            "--bim-file", "/ref/1000G.EUR.QC.@.bim",
            "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
            "--threads", $Threads,
            "--exclude-ranges", "MHC"
        )
    } else {
        Write-Host "[SKIP] test2 $Prefix"
    }
}

$Traits = @("lipid8_F2", "AD", "nonlipid8_F1", "lipid8_F1", "PD")
foreach ($Trait in $Traits) {
    Run-Fit1 $Trait
}

Run-Fit2Test2 "lipid8_F2" "AD"
Run-Fit2Test2 "nonlipid8_F1" "PD"
Run-Fit2Test2 "lipid8_F1" "PD"

Write-Host "MiXeR rep1 key analyses complete."
