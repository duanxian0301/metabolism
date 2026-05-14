from __future__ import annotations

import csv
import math
from pathlib import Path


PROJECT_ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
CLEAN_REPORT = PROJECT_ROOT / "results" / "01_clean_factor_inputs" / "clean_factor_inputs_report.tsv"
RESULTS_DIR = PROJECT_ROOT / "results" / "02_clean_input_qc"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DISEASE_TRAITS = {
    "AD": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
    "PD": Path(r"D:\文章\4NDD\NDDGWAS\PD.txt"),
    "LBD": Path(r"D:\文章\4NDD\NDDGWAS\LBD.txt"),
}

REQUIRED_COLUMNS = ["SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"]


def to_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        out = float(value)
    except ValueError:
        return None
    if not math.isfinite(out):
        return None
    return out


def inspect(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        missing = [col for col in REQUIRED_COLUMNS if col not in header]
        rows = 0
        missing_chr_bp = 0
        bad_se = 0
        bad_p = 0
        n_values = set()
        min_chr = None
        max_chr = None
        for row in reader:
            rows += 1
            if len(n_values) < 10 and row.get("N"):
                n_values.add(str(row["N"]))
            chrom = to_float(row.get("CHR"))
            bp = to_float(row.get("BP"))
            if chrom is None or bp is None:
                missing_chr_bp += 1
            else:
                c = int(chrom)
                min_chr = c if min_chr is None else min(min_chr, c)
                max_chr = c if max_chr is None else max(max_chr, c)
            se = to_float(row.get("SE"))
            p = to_float(row.get("P"))
            if se is None or se <= 0:
                bad_se += 1
            if p is None or p <= 0 or p > 1:
                bad_p += 1
    return {
        "file": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size,
        "row_count": rows,
        "columns": ",".join(header),
        "missing_required": ",".join(missing),
        "n_unique_values_sample": ",".join(sorted(n_values)),
        "min_chr": "" if min_chr is None else min_chr,
        "max_chr": "" if max_chr is None else max_chr,
        "missing_chr_bp": missing_chr_bp,
        "bad_se": bad_se,
        "bad_p": bad_p,
    }


def main() -> None:
    with CLEAN_REPORT.open("r", encoding="utf-8", newline="") as handle:
        clean_rows = list(csv.DictReader(handle, delimiter="\t"))

    rows = []
    for row in clean_rows:
        qc = inspect(Path(row["clean_file"]))
        rows.append({"trait": row["trait"], "group": "metabolic_factor", **qc})
    for trait, path in DISEASE_TRAITS.items():
        qc = inspect(path)
        rows.append({"trait": trait, "group": "ndd_trait", **qc})

    manifest = RESULTS_DIR / "clean_input_qc_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    blockers = []
    for row in rows:
        for key in ["missing_required", "missing_chr_bp", "bad_se", "bad_p"]:
            value = row[key]
            if value not in ("", 0, "0"):
                blockers.append(f"{row['trait']}: {key}={value}")

    blocker_path = RESULTS_DIR / "clean_input_qc_blockers.txt"
    blocker_path.write_text("\n".join(blockers) + ("\n" if blockers else ""), encoding="utf-8")
    print(f"Wrote QC manifest: {manifest}")
    if blockers:
        print(f"Found blockers: {blocker_path}")
        for item in blockers:
            print("BLOCKER:", item)
        raise SystemExit(1)
    print("No clean input blockers found.")


if __name__ == "__main__":
    main()
