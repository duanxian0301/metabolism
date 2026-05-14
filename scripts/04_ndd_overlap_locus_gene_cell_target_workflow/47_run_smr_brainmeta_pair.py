from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
SMR_EXE = Path(r"D:\SMR\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe")
LD_BFILE = Path(r"D:\SMR\g1000\g1000_eur")
BRAINMETA_PREFIX = Path(r"D:\SMR\Brain-mMeta\Brain-mMeta")


def run_trait(pair: str, result_dir_name: str, sumstat_dir_name: str, trait: str, disable_freq_ck: bool) -> None:
    gwas_file = ROOT / "results" / result_dir_name / "input" / "smr_sumstats" / f"{trait}.smr.txt"
    if not gwas_file.exists():
        raise FileNotFoundError(gwas_file)
    out_dir = ROOT / "results" / result_dir_name / "single_trait" / trait / "brainmeta"
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
    if disable_freq_ck:
        cmd.append("--disable-freq-ck")
    log_path = Path(str(out_prefix) + ".run.log")
    with log_path.open("w", encoding="utf-8") as log_f:
        proc = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, text=True)
    print(f"{pair}\t{trait}\texit_code={proc.returncode}\tlog={log_path}\tout={out_prefix}.smr")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair", required=True)
    parser.add_argument("--result-dir-name", required=True)
    parser.add_argument("--traits", nargs="+", required=True)
    parser.add_argument("--disable-freq-ck-traits", nargs="*", default=[])
    args = parser.parse_args()
    disable_set = set(args.disable_freq_ck_traits)
    for trait in args.traits:
        run_trait(args.pair, args.result_dir_name, args.result_dir_name, trait, trait in disable_set)


if __name__ == "__main__":
    main()
