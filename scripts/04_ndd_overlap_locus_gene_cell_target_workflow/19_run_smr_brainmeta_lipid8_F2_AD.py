from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
SMR_EXE = Path(r"D:\SMR\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe")
LD_BFILE = Path(r"D:\SMR\g1000\g1000_eur")
SUMSTAT_DIR = ROOT / "results" / "11_smr_lipid8_F2_AD" / "input" / "smr_sumstats"
BRAINMETA_PREFIX = Path(r"D:\SMR\Brain-mMeta\Brain-mMeta")
OUT_ROOT = ROOT / "results" / "11_smr_lipid8_F2_AD"


def run_trait(trait: str) -> None:
    gwas_file = SUMSTAT_DIR / f"{trait}.smr.txt"
    if not gwas_file.exists():
        raise FileNotFoundError(gwas_file)
    out_dir = OUT_ROOT / "single_trait" / trait / "brainmeta"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / f"{trait}_BrainMeta"
    cmd = [
        str(SMR_EXE),
        "--bfile",
        str(LD_BFILE),
        "--gwas-summary",
        str(gwas_file),
        "--beqtl-summary",
        str(BRAINMETA_PREFIX),
        "--out",
        str(out_prefix),
    ]
    if trait == "lipid8_F2":
        cmd.append("--disable-freq-ck")
    log_path = Path(str(out_prefix) + ".run.log")
    with log_path.open("w", encoding="utf-8") as log_f:
        proc = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, text=True)
    print(f"{trait}\texit_code={proc.returncode}\tlog={log_path}\tout={out_prefix}.smr")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--traits", nargs="+", default=["lipid8_F2", "AD"])
    args = parser.parse_args()
    for trait in args.traits:
        run_trait(trait)


if __name__ == "__main__":
    main()
