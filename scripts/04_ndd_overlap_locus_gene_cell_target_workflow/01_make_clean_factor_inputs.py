from __future__ import annotations

import csv
import math
from pathlib import Path


PROJECT_ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
INPUT_LOCK = PROJECT_ROOT / "results" / "00_input_lock" / "metabolic_ad_pdlbd_input_manifest.tsv"
CLEAN_DIR = PROJECT_ROOT / "work" / "clean_factor_inputs"
CLEAN_DIR.mkdir(parents=True, exist_ok=True)

REPORT_DIR = PROJECT_ROOT / "results" / "01_clean_factor_inputs"
REPORT_DIR.mkdir(parents=True, exist_ok=True)


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


def normal_p_from_z(z: float) -> float:
    # Two-sided normal tail via erfc is stable until very extreme z; clamp to
    # a positive floor so downstream tools never see P=0.
    p = math.erfc(abs(z) / math.sqrt(2.0))
    if not math.isfinite(p) or p <= 0.0:
        return 1e-300
    return max(p, 1e-300)


def clean_one(trait: str, source: Path) -> dict[str, object]:
    out_path = CLEAN_DIR / f"{trait}_clean.txt"
    total = 0
    written = 0
    dropped_missing_effect = 0
    dropped_bad_se = 0
    dropped_bad_p = 0
    p_zero_recomputed = 0
    p_clamped_min = 0

    with source.open("r", encoding="utf-8-sig", errors="replace", newline="") as fin, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "FRQ" in fieldnames and "FREQ" not in fieldnames:
            fieldnames = ["FREQ" if col == "FRQ" else col for col in fieldnames]
        required = ["SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"]
        missing = [col for col in required if col not in fieldnames]
        if missing:
            raise ValueError(f"{trait}: missing columns after FRQ/FREQ normalization: {missing}")

        writer = csv.DictWriter(fout, fieldnames=required, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            total += 1
            if "FRQ" in row and "FREQ" not in row:
                row["FREQ"] = row.get("FRQ", "")

            beta = to_float(row.get("BETA"))
            se = to_float(row.get("SE"))
            p = to_float(row.get("P"))

            if beta is None:
                dropped_missing_effect += 1
                continue
            if se is None or se <= 0:
                dropped_bad_se += 1
                continue

            if p is None or p > 1:
                dropped_bad_p += 1
                continue
            if p <= 0:
                p = normal_p_from_z(beta / se)
                p_zero_recomputed += 1
                if p == 1e-300:
                    p_clamped_min += 1

            out_row = {col: row.get(col, "") for col in required}
            out_row["P"] = f"{p:.12g}"
            writer.writerow(out_row)
            written += 1

    return {
        "trait": trait,
        "source_file": str(source),
        "clean_file": str(out_path),
        "total_rows": total,
        "written_rows": written,
        "dropped_rows": total - written,
        "dropped_missing_beta": dropped_missing_effect,
        "dropped_bad_se": dropped_bad_se,
        "dropped_bad_p": dropped_bad_p,
        "p_zero_recomputed": p_zero_recomputed,
        "p_clamped_to_1e_300": p_clamped_min,
    }


def main() -> None:
    with INPUT_LOCK.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    factor_rows = [row for row in rows if row["group"] == "metabolic_factor"]
    reports = []
    for row in factor_rows:
        reports.append(clean_one(row["trait"], Path(row["standard_file"])))

    report_path = REPORT_DIR / "clean_factor_inputs_report.tsv"
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(reports[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(reports)

    print(f"Wrote clean files to: {CLEAN_DIR}")
    print(f"Wrote report: {report_path}")
    for row in reports:
        print(
            row["trait"],
            "written=", row["written_rows"],
            "dropped=", row["dropped_rows"],
            "p0_recomputed=", row["p_zero_recomputed"],
        )


if __name__ == "__main__":
    main()
