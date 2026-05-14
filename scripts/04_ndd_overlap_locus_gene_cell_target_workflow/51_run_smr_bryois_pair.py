from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import psutil


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
SMR_EXE = Path(r"D:\SMR\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe")
LD_BFILE = Path(r"D:\SMR\g1000\g1000_eur")
BRYOIS_PREPARED = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\bryois2022_celltype_eqtl\prepared")

CELLTYPES = [
    "Astrocytes",
    "Microglia",
    "Endothelial.cells",
    "Excitatory.neurons",
    "Inhibitory.neurons",
    "OPCs...COPs",
    "Oligodendrocytes",
    "Pericytes",
]
CHROMS = list(range(1, 23))


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


def out_prefix(result_dir_name: str, trait: str, celltype: str, chrom: int) -> Path:
    return ROOT / "results" / result_dir_name / "single_trait" / trait / "bryois2022_celltype" / celltype / f"{trait}_{celltype}_chr{chrom}"


def output_done(result_dir_name: str, trait: str, celltype: str, chrom: int) -> bool:
    out = Path(str(out_prefix(result_dir_name, trait, celltype, chrom)) + ".smr")
    log = Path(str(out_prefix(result_dir_name, trait, celltype, chrom)) + ".run.log")
    if not out.exists() or out.stat().st_size == 0:
        return False
    if log.exists() and "Analysis completed:" in log.read_text(encoding="utf-8", errors="ignore"):
        return True
    return out.stat().st_size > 0


def ref_prefix(celltype: str, chrom: int) -> Path:
    return BRYOIS_PREPARED / celltype / f"chr{chrom}" / f"{celltype}_chr{chrom}_filtered"


def build_cmd(result_dir_name: str, trait: str, celltype: str, chrom: int, disable_freq_ck: bool) -> list[str]:
    pair = result_dir_name.replace("15_smr_bryois_", "", 1)
    gwas = ROOT / "results" / f"11_smr_{pair}" / "input" / "smr_sumstats" / f"{trait}.smr.txt"
    ref = ref_prefix(celltype, chrom)
    for suffix in [".besd", ".epi", ".esi"]:
        if not Path(str(ref) + suffix).exists():
            raise FileNotFoundError(f"Missing Bryois reference: {ref}{suffix}")
    out = out_prefix(result_dir_name, trait, celltype, chrom)
    out.parent.mkdir(parents=True, exist_ok=True)
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


def run_one(result_dir_name: str, trait: str, celltype: str, chrom: int, disable_freq_ck: bool) -> int:
    out = out_prefix(result_dir_name, trait, celltype, chrom)
    log = Path(str(out) + ".run.log")
    cmd = build_cmd(result_dir_name, trait, celltype, chrom, disable_freq_ck)
    with log.open("w", encoding="utf-8") as log_f:
        proc = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, text=True)
    return proc.returncode


def all_tasks(traits: list[str], celltypes: list[str], chroms: list[int]) -> list[tuple[str, str, int]]:
    return [(trait, cell, chrom) for cell in celltypes for chrom in chroms for trait in traits]


def write_status(result_dir_name: str, tasks: list[tuple[str, str, int]]) -> None:
    out_root = ROOT / "results" / result_dir_name
    out_root.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "trait": trait,
            "celltype": cell,
            "chrom": chrom,
            "done": output_done(result_dir_name, trait, cell, chrom),
            "output": str(Path(str(out_prefix(result_dir_name, trait, cell, chrom)) + ".smr")),
        }
        for trait, cell, chrom in tasks
    ]
    pd.DataFrame(rows).to_csv(out_root / "bryois_queue_status.tsv", sep="\t", index=False)


def queue(args: argparse.Namespace) -> None:
    tasks = all_tasks(args.traits, args.celltypes, args.chroms)
    disable_set = set(args.disable_freq_ck_traits)
    while True:
        pending = [t for t in tasks if not output_done(args.result_dir_name, *t)]
        write_status(args.result_dir_name, tasks)
        if not pending:
            print("status\tcomplete", flush=True)
            break
        launched = 0
        for trait, cell, chrom in pending:
            if count_active_smr() >= args.max_total_smr:
                break
            child_cmd = [
                sys.executable,
                __file__,
                "single",
                "--result-dir-name",
                args.result_dir_name,
                "--trait",
                trait,
                "--celltype",
                cell,
                "--chrom",
                str(chrom),
            ]
            if trait in disable_set:
                child_cmd.append("--disable-freq-ck")
            subprocess.Popen(child_cmd, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
            launched += 1
            time.sleep(args.launch_gap_seconds)
        print(f"pending\t{len(pending)}\tlaunched\t{launched}\tactive_smr\t{count_active_smr()}", flush=True)
        time.sleep(args.poll_seconds)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)

    p_single = sub.add_parser("single")
    p_single.add_argument("--result-dir-name", required=True)
    p_single.add_argument("--trait", required=True)
    p_single.add_argument("--celltype", required=True, choices=CELLTYPES)
    p_single.add_argument("--chrom", required=True, type=int, choices=CHROMS)
    p_single.add_argument("--disable-freq-ck", action="store_true")

    p_queue = sub.add_parser("queue")
    p_queue.add_argument("--result-dir-name", required=True)
    p_queue.add_argument("--traits", nargs="+", required=True)
    p_queue.add_argument("--disable-freq-ck-traits", nargs="*", default=[])
    p_queue.add_argument("--celltypes", nargs="+", default=CELLTYPES, choices=CELLTYPES)
    p_queue.add_argument("--chroms", nargs="+", type=int, default=CHROMS, choices=CHROMS)
    p_queue.add_argument("--max-total-smr", type=int, default=6)
    p_queue.add_argument("--poll-seconds", type=int, default=90)
    p_queue.add_argument("--launch-gap-seconds", type=int, default=2)
    args = parser.parse_args()
    if args.mode == "single":
        code = run_one(args.result_dir_name, args.trait, args.celltype, args.chrom, args.disable_freq_ck)
        print(f"{args.trait}\t{args.celltype}\tchr{args.chrom}\texit_code={code}")
        raise SystemExit(code)
    queue(args)


if __name__ == "__main__":
    main()
