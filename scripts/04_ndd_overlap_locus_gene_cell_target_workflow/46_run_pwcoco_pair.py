from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
PWC = Path(r"D:\PWCoCo\build\Release\pwcoco.exe")


def run_one(base: Path, pair: str, row: pd.Series) -> dict[str, str | int]:
    region_id = row["region_id"]
    run_dir = base / "runs" / region_id
    run_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = run_dir / f"{pair}_{region_id}"
    log_path = run_dir / f"{pair}_{region_id}.log.txt"
    coloc_path = Path(str(out_prefix) + ".coloc")
    if coloc_path.exists():
        coloc_path.unlink()

    cmd = [
        str(PWC),
        "--bfile",
        str(row["bfile"]),
        "--sum_stats1",
        str(row["sum_stats1"]),
        "--sum_stats2",
        str(row["sum_stats2"]),
        "--chr",
        str(int(row["chrnum"])),
        "--out",
        str(out_prefix),
        "--log",
        str(log_path),
        "--threads",
        "8",
        "--init_h4",
        "100",
        "--verbose",
    ]
    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (run_dir / f"{pair}_{region_id}.stdout.txt").write_text(proc.stdout, encoding="utf-8", errors="replace")
    return {
        "region_id": region_id,
        "exit_code": proc.returncode,
        "coloc_file": str(coloc_path),
        "log_file": str(log_path),
    }


def summarize(base: Path, pair: str, run_rows: list[dict[str, str | int]]) -> None:
    frames = []
    for rr in run_rows:
        path = Path(str(rr["coloc_file"]))
        if not path.exists() or path.stat().st_size == 0:
            continue
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "region_id", rr["region_id"])
        frames.append(df)
    if not frames:
        return
    out = pd.concat(frames, ignore_index=True)
    out["result_type"] = out.apply(
        lambda r: "unconditioned"
        if r["SNP1"] == "unconditioned" and r["SNP2"] == "unconditioned"
        else (
            "conditioned_both"
            if r["SNP1"] != "unconditioned" and r["SNP2"] != "unconditioned"
            else ("conditioned_trait" if r["SNP1"] != "unconditioned" else "conditioned_disease")
        ),
        axis=1,
    )
    out.to_csv(base / f"pwcoco_{pair}_summary.tsv", sep="\t", index=False)
    best = out.sort_values(["region_id", "H4"], ascending=[True, False]).groupby("region_id", as_index=False).head(1)
    best.to_csv(base / f"pwcoco_{pair}_best_h4.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair", required=True)
    parser.add_argument("--trait1", required=True)
    parser.add_argument("--trait2", required=True)
    args = parser.parse_args()

    base = ROOT / "results" / f"10_pwcoco_{args.pair}"
    manifest = pd.read_csv(base / "pwcoco_region_input_manifest.tsv", sep="\t")
    trait1_n_col = "nsnps_trait1" if "nsnps_trait1" in manifest.columns else f"nsnps_{args.trait1}"
    trait2_n_col = "nsnps_trait2" if "nsnps_trait2" in manifest.columns else f"nsnps_{args.trait2}"
    manifest = manifest[(manifest[trait1_n_col] >= 50) & (manifest[trait2_n_col] >= 50)].copy()
    run_rows = [run_one(base, args.pair, row) for _, row in manifest.iterrows()]
    if run_rows:
        with (base / "pwcoco_run_manifest.tsv").open("w", newline="", encoding="utf-8") as fout:
            writer = csv.DictWriter(fout, fieldnames=list(run_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(run_rows)
        summarize(base, args.pair, run_rows)
    print(f"Wrote PWCoCo run outputs to {base}")


if __name__ == "__main__":
    main()
