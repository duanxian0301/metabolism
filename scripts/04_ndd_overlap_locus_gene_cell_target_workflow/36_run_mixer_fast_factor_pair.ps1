param(
    [Parameter(Mandatory = $true)]
    [string]$Trait1,
    [Parameter(Mandatory = $true)]
    [string]$Trait2,
    [string]$RunName = "fast_rep1",
    [string]$Threads = "12",
    [string]$Rep = "1"
)

$ErrorActionPreference = "Stop"

$ProjectRoot = "D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd"
$MixerRoot = Join-Path $ProjectRoot "work\04_mixer_inputs"
$RunDir = Join-Path $MixerRoot "runs\$RunName"
$RefDir = "D:\mixer\mixer\reference\1000G_EUR_Phase3_plink"
$Image = "ghcr.io/precimed/gsa-mixer:2.2.1"
$FastArgs = @("--fit-sequence", "diffevo-fast", "neldermead-fast", "--diffevo-fast-repeats", "2")

New-Item -ItemType Directory -Force $RunDir | Out-Null

function Invoke-Mixer {
    param([Parameter(Mandatory = $true)][string[]]$ArgsList)
    docker run --rm `
        -v "${MixerRoot}:/work" `
        -v "${RefDir}:/ref" `
        -w "/work/runs/$RunName" `
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
    param([string]$T1, [string]$T2)
    $Prefix = "${T1}_vs_${T2}"
    $FitJson = Join-Path $RunDir "$Prefix.fit.fast.rep$Rep.json"
    $TestJson = Join-Path $RunDir "$Prefix.test.fast.rep$Rep.json"

    if (-not (Test-Path $FitJson)) {
        Write-Host "[RUN] fit2 fast $Prefix"
        $ArgsList = @(
            "fit2"
        ) + $FastArgs + @(
            "--trait1-file", "/work/inputs/$T1.sumstats.gz",
            "--trait2-file", "/work/inputs/$T2.sumstats.gz",
            "--trait1-params-file", "$T1.fit.fast.rep$Rep.json",
            "--trait2-params-file", "$T2.fit.fast.rep$Rep.json",
            "--out", "$Prefix.fit.fast.rep$Rep",
            "--extract", "/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$Rep.snps",
            "--bim-file", "/ref/1000G.EUR.QC.@.bim",
            "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
            "--threads", $Threads,
            "--exclude-ranges", "MHC"
        )
        Invoke-Mixer $ArgsList
    }
    else {
        Write-Host "[SKIP] fit2 fast $Prefix"
    }

    if (-not (Test-Path $TestJson)) {
        Write-Host "[RUN] test2 fast $Prefix"
        Invoke-Mixer @(
            "test2",
            "--trait1-file", "/work/inputs/$T1.sumstats.gz",
            "--trait2-file", "/work/inputs/$T2.sumstats.gz",
            "--load-params", "$Prefix.fit.fast.rep$Rep.json",
            "--out", "$Prefix.test.fast.rep$Rep",
            "--bim-file", "/ref/1000G.EUR.QC.@.bim",
            "--ld-file", "/ref/1000G.EUR.QC.@.run4.ld",
            "--threads", $Threads,
            "--exclude-ranges", "MHC"
        )
    }
    else {
        Write-Host "[SKIP] test2 fast $Prefix"
    }
}

Run-Fit1Fast $Trait1
Run-Fit1Fast $Trait2
Run-Fit2Test2Fast $Trait1 $Trait2

Write-Host "Fast MiXeR complete for ${Trait1}_vs_${Trait2}."
