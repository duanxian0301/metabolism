$ErrorActionPreference = "Stop"

$ProjectRoot = "D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd"
$MixerRoot = Join-Path $ProjectRoot "work\04_mixer_inputs"
$RunDir = Join-Path $MixerRoot "runs\fast_rep1"
$RefDir = "D:\mixer\mixer\reference\1000G_EUR_Phase3_plink"
$Image = "ghcr.io/precimed/gsa-mixer:2.2.1"
$Threads = "12"
$Rep = "1"
$FastArgs = @("--fit-sequence", "diffevo-fast", "neldermead-fast", "--diffevo-fast-repeats", "2")

New-Item -ItemType Directory -Force $RunDir | Out-Null

function Invoke-Mixer {
    param([Parameter(Mandatory=$true)][string[]]$ArgsList)
    docker run --rm `
        -v "${MixerRoot}:/work" `
        -v "${RefDir}:/ref" `
        -w "/work/runs/fast_rep1" `
        $Image `
        python /tools/mixer/precimed/mixer.py @ArgsList
}

function Run-Fit1Fast {
    param([string]$Trait)
    $JsonPath = Join-Path $RunDir "$Trait.fit.fast.rep$Rep.json"
    if (Test-Path $JsonPath) {
        Write-Host "[SKIP] fit1 fast $Trait"
        return
    }
    Write-Host "[RUN] fit1 fast $Trait"
    $ArgsList = @(
        "fit1"
    ) + $FastArgs + @(
        "--trait1-file", "/work/inputs/$Trait.sumstats.gz",
        "--out", "$Trait.fit.fast.rep$Rep",
        "--extract", "/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$Rep.snps",
        "--bim-file", "/ref/1000G.EUR.QC.@.bim",
        "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
        "--threads", $Threads,
        "--exclude-ranges", "MHC"
    )
    Invoke-Mixer $ArgsList
}

function Run-Fit2Test2Fast {
    param([string]$Trait1, [string]$Trait2)
    $Prefix = "${Trait1}_vs_${Trait2}"
    $FitJson = Join-Path $RunDir "$Prefix.fit.fast.rep$Rep.json"
    $TestJson = Join-Path $RunDir "$Prefix.test.fast.rep$Rep.json"

    if (-not (Test-Path $FitJson)) {
        Write-Host "[RUN] fit2 fast $Prefix"
        $ArgsList = @(
            "fit2"
        ) + $FastArgs + @(
            "--trait1-file", "/work/inputs/$Trait1.sumstats.gz",
            "--trait2-file", "/work/inputs/$Trait2.sumstats.gz",
            "--trait1-params-file", "$Trait1.fit.fast.rep$Rep.json",
            "--trait2-params-file", "$Trait2.fit.fast.rep$Rep.json",
            "--out", "$Prefix.fit.fast.rep$Rep",
            "--extract", "/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$Rep.snps",
            "--bim-file", "/ref/1000G.EUR.QC.@.bim",
            "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
            "--threads", $Threads,
            "--exclude-ranges", "MHC"
        )
        Invoke-Mixer $ArgsList
    } else {
        Write-Host "[SKIP] fit2 fast $Prefix"
    }

    if (-not (Test-Path $TestJson)) {
        Write-Host "[RUN] test2 fast $Prefix"
        Invoke-Mixer @(
            "test2",
            "--trait1-file", "/work/inputs/$Trait1.sumstats.gz",
            "--trait2-file", "/work/inputs/$Trait2.sumstats.gz",
            "--load-params", "$Prefix.fit.fast.rep$Rep.json",
            "--out", "$Prefix.test.fast.rep$Rep",
            "--bim-file", "/ref/1000G.EUR.QC.@.bim",
            "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
            "--threads", $Threads,
            "--exclude-ranges", "MHC"
        )
    } else {
        Write-Host "[SKIP] test2 fast $Prefix"
    }
}

Run-Fit1Fast "lipid8_F2"
Run-Fit1Fast "AD"
Run-Fit2Test2Fast "lipid8_F2" "AD"

Write-Host "Fast MiXeR lipid8_F2 vs AD complete."
