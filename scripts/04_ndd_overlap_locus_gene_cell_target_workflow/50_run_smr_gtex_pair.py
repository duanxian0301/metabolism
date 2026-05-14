from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import psutil


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
SMR_EXE = Path(r"D:\SMR\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe")
LD_BFILE = Path(r"D:\SMR\g1000\g1000_eur")
GTEX_BRAIN_DIR = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\smr_gtex_v8_brain_hg19\brain_lite\eQTL_besd_lite")

GTEX_BRAIN_TISSUES = [
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
]


def count_active_smr() -> int:
    total = 0
    for proc in psutil.process_iter(["name"]):
        try:
            name = (proc.info.get("name") or "").lower()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
        if "smr" in name:
            total += 1
    return total


def out_prefix(result_dir_name: str, trait: str, tissue: str) -> Path:
    return ROOT / "results" / result_dir_name / "single_trait" / trait / "gtex_v8_brain" / f"{trait}_{tissue}"


def output_done(result_dir_name: str, trait: str, tissue: str) -> bool:
    out = Path(str(out_prefix(result_dir_name, trait, tissue)) + ".smr")
    log = Path(str(out_prefix(result_dir_name, trait, tissue)) + ".run.log")
    if not out.exists() or out.stat().st_size == 0:
        return False
    if log.exists() and "Analysis completed:" in log.read_text(encoding="utf-8", errors="ignore"):
        return True
    return out.stat().st_size > 0


def build_cmd(result_dir_name: str, trait: str, tissue: str, disable_freq_ck: bool) -> list[str]:
    pair = result_dir_name.replace("14_smr_gtex_", "", 1)
    gwas = ROOT / "results" / f"11_smr_{pair}" / "input" / "smr_sumstats" / f"{trait}.smr.txt"
    ref = GTEX_BRAIN_DIR / f"{tissue}.lite"
    out = out_prefix(result_dir_name, trait, tissue)
    out.parent.mkdir(parents=True, exist_ok=True)
    for suffix in [".besd", ".epi", ".esi"]:
        if not Path(str(ref) + suffix).exists():
            raise FileNotFoundError(f"Missing GTEx reference: {ref}{suffix}")
    cmd = [
        str(SMR_EXE),
        "--bfile",
        str(LD_BFILE),
        "--gwas-summary",
        str(gwas),
        "--beqtl-summary",
        str(ref),
        "--out",
        str(out),
    ]
    if disable_freq_ck:
        cmd.append("--disable-freq-ck")
    return cmd


def run_one(result_dir_name: str, trait: str, tissue: str, disable_freq_ck: bool) -> int:
    out = out_prefix(result_dir_name, trait, tissue)
    log = Path(str(out) + ".run.log")
    cmd = build_cmd(result_dir_name, trait, tissue, disable_freq_ck)
    with log.open("w", encoding="utf-8") as log_f:
        proc = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, text=True)
    return proc.returncode


def queue(args: argparse.Namespace) -> None:
    out_root = ROOT / "results" / args.result_dir_name
    out_root.mkdir(parents=True, exist_ok=True)
    status_path = out_root / "gtex_queue_status.tsv"
    disable_set = set(args.disable_freq_ck_traits)
    pending = [(trait, tissue) for trait in args.traits for tissue in args.tissues if not output_done(args.result_dir_name, trait, tissue)]
    while pending:
        launched = 0
        for trait, tissue in pending:
            if output_done(args.result_dir_name, trait, tissue):
                continue
            if count_active_smr() >= args.max_total_smr:
                continue
            child_cmd = [
                sys.executable,
                __file__,
                "single",
                "--result-dir-name",
                args.result_dir_name,
                "--trait",
                trait,
                "--tissue",
                tissue,
            ]
            if trait in disable_set:
                child_cmd.append("--disable-freq-ck")
            subprocess.Popen(child_cmd, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
            launched += 1
            time.sleep(args.launch_gap_seconds)
        rows = []
        for trait in args.traits:
            for tissue in args.tissues:
                rows.append(
                    {
                        "trait": trait,
                        "tissue": tissue,
                        "done": output_done(args.result_dir_name, trait, tissue),
                        "output": str(Path(str(out_prefix(args.result_dir_name, trait, tissue)) + ".smr")),
                    }
                )
        import pandas as pd

        pd.DataFrame(rows).to_csv(status_path, sep="\t", index=False)
        pending = [(trait, tissue) for trait in args.traits for tissue in args.tissues if not output_done(args.result_dir_name, trait, tissue)]
        print(f"pending\t{len(pending)}\tlaunched\t{launched}\tactive_smr\t{count_active_smr()}", flush=True)
        if pending:
            time.sleep(args.poll_seconds)
    print("status\tcomplete", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)

    p_single = sub.add_parser("single")
    p_single.add_argument("--result-dir-name", required=True)
    p_single.add_argument("--trait", required=True)
    p_single.add_argument("--tissue", required=True, choices=GTEX_BRAIN_TISSUES)
    p_single.add_argument("--disable-freq-ck", action="store_true")

    p_queue = sub.add_parser("queue")
    p_queue.add_argument("--result-dir-name", required=True)
    p_queue.add_argument("--traits", nargs="+", required=True)
    p_queue.add_argument("--disable-freq-ck-traits", nargs="*", default=[])
    p_queue.add_argument("--tissues", nargs="+", default=GTEX_BRAIN_TISSUES, choices=GTEX_BRAIN_TISSUES)
    p_queue.add_argument("--max-total-smr", type=int, default=2)
    p_queue.add_argument("--poll-seconds", type=int, default=60)
    p_queue.add_argument("--launch-gap-seconds", type=int, default=5)
    args = parser.parse_args()
    if args.mode == "single":
        code = run_one(args.result_dir_name, args.trait, args.tissue, args.disable_freq_ck)
        print(f"{args.trait}\t{args.tissue}\texit_code={code}")
        raise SystemExit(code)
    queue(args)


if __name__ == "__main__":
    main()
