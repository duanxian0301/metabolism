from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
PWC = Path(r"D:\PWCoCo\build\Release\pwcoco.exe")
BASE = ROOT / "results" / "10_pwcoco_lipid8_F2_AD"
MANIFEST = BASE / "pwcoco_region_input_manifest.tsv"
RUN_ROOT = BASE / "runs"


def run_one(row: pd.Series) -> dict[str, str | int]:
    region_id = row["region_id"]
    run_dir = RUN_ROOT / region_id
    run_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = run_dir / f"lipid8_F2_AD_{region_id}"
    log_path = run_dir / f"lipid8_F2_AD_{region_id}.log.txt"
    coloc_path = Path(str(out_prefix) + ".coloc")
    if coloc_path.exists():
      coloc_path.unlink()

    cmd = [
        str(PWC),
        "--bfile",
        row["bfile"],
        "--sum_stats1",
        row["sum_stats1"],
        "--sum_stats2",
        row["sum_stats2"],
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
    (run_dir / f"lipid8_F2_AD_{region_id}.stdout.txt").write_text(proc.stdout, encoding="utf-8", errors="replace")
    return {
        "region_id": region_id,
        "exit_code": proc.returncode,
        "coloc_file": str(coloc_path),
        "log_file": str(log_path),
    }


def summarize(run_rows: list[dict[str, str | int]]) -> None:
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
    out.to_csv(BASE / "pwcoco_lipid8_F2_AD_summary.tsv", sep="\t", index=False)
    best = out.sort_values(["region_id", "H4"], ascending=[True, False]).groupby("region_id", as_index=False).head(1)
    best.to_csv(BASE / "pwcoco_lipid8_F2_AD_best_h4.tsv", sep="\t", index=False)


def main() -> None:
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(MANIFEST, sep="\t")
    manifest = manifest[(manifest["nsnps_lipid8_F2"] >= 50) & (manifest["nsnps_AD"] >= 50)].copy()
    run_rows = [run_one(row) for _, row in manifest.iterrows()]
    with (BASE / "pwcoco_run_manifest.tsv").open("w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=list(run_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(run_rows)
    summarize(run_rows)
    print(f"Wrote PWCoCo run outputs to {BASE}")


if __name__ == "__main__":
    main()
