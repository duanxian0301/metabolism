import csv
import math
from pathlib import Path
from statistics import mean


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")

FACTOR_FILES = {
    "lipid8_F1": ROOT / "step14_native_wsl_usergwas_final8_results" / "merged_lipid_final8" / "standard_txt" / "lipid_final8_F1_standard.txt",
    "lipid8_F2": ROOT / "step14_native_wsl_usergwas_final8_results" / "merged_lipid_final8" / "standard_txt" / "lipid_final8_F2_standard.txt",
    "lipid8_F3": ROOT / "step14_native_wsl_usergwas_final8_results" / "merged_lipid_final8" / "standard_txt" / "lipid_final8_F3_standard.txt",
    "nonlipid8_F1": ROOT / "step22_native_wsl_usergwas_nonlipid8_results" / "merged_nonlipid_final8" / "standard_txt" / "nonlipid_final8_F1_standard.txt",
    "nonlipid8_F2": ROOT / "step22_native_wsl_usergwas_nonlipid8_results" / "merged_nonlipid_final8" / "standard_txt" / "nonlipid_final8_F2_standard.txt",
    "nonlipid8_F3": ROOT / "step22_native_wsl_usergwas_nonlipid8_results" / "merged_nonlipid_final8" / "standard_txt" / "nonlipid_final8_F3_standard.txt",
}

METADATA_FILES = [
    ROOT / "step14_native_wsl_usergwas_final8_results" / "merged_lipid_final8" / "standard_txt" / "lipid_final8_standard_txt_metadata.tsv",
    ROOT / "step22_native_wsl_usergwas_nonlipid8_results" / "merged_nonlipid_final8" / "standard_txt" / "nonlipid_final8_standard_txt_metadata.tsv",
]

SUMMARY_FILE = ROOT / "factor_standard_txt_neff_summary.tsv"


def estimate_neff(path: Path):
    values = []
    total_rows = 0
    used_rows = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total_rows += 1
            frq = row.get("FRQ", "")
            se = row.get("SE", "")
            if not frq or not se:
                continue
            maf = float(frq)
            se_val = float(se)
            if 0.1 <= maf <= 0.4 and se_val > 0:
                values.append(1.0 / (2.0 * maf * (1.0 - maf) * (se_val ** 2)))
                used_rows += 1

    if not values:
        raise RuntimeError(f"No valid SNPs for Neff estimation in {path}")

    neff_mean = int(round(mean(values)))
    values.sort()
    neff_median = int(round(values[len(values) // 2]))
    return {
        "total_rows": total_rows,
        "used_rows": used_rows,
        "neff_mean": neff_mean,
        "neff_median": neff_median,
    }


def rewrite_n_column(path: Path, n_value: int):
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with path.open("r", encoding="utf-8-sig", newline="") as src, tmp_path.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = reader.fieldnames
        if fieldnames is None or "N" not in fieldnames:
            raise RuntimeError(f"Missing N column in {path}")
        writer = csv.DictWriter(dst, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            row["N"] = str(n_value)
            writer.writerow(row)
    tmp_path.replace(path)


def update_metadata(path: Path, neff_map):
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if fieldnames is None:
            raise RuntimeError(f"Missing header in {path}")
        name = path.name.lower()
        if name.startswith("nonlipid_final8"):
            module = "nonlipid8"
        elif name.startswith("lipid_final8"):
            module = "lipid8"
        else:
            raise RuntimeError(f"Could not infer module from {path}")
        for row in reader:
            factor = row["factor"]
            key = f"{module}_{factor}"
            if key not in neff_map:
                raise RuntimeError(f"Could not map {key} for {path}")
            row["n_used"] = str(neff_map[key]["neff_mean"])
            note = row.get("notes", "")
            extra = " N replaced with factor-specific Neff estimated as mean[1/(2*MAF*(1-MAF)*SE^2)] across SNPs with 0.1<=MAF<=0.4."
            if extra not in note:
                row["notes"] = note + extra
            rows.append(row)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    summary_rows = []
    neff_map = {}

    for key, path in FACTOR_FILES.items():
        stats = estimate_neff(path)
        neff_map[key] = stats
        rewrite_n_column(path, stats["neff_mean"])
        summary_rows.append(
            {
                "factor_id": key,
                "standard_file": str(path),
                "n_original_label": 599249,
                "neff_mean": stats["neff_mean"],
                "neff_median": stats["neff_median"],
                "total_rows": stats["total_rows"],
                "rows_used_for_neff": stats["used_rows"],
                "maf_window": "0.1-0.4",
                "formula": "mean(1/(2*MAF*(1-MAF)*SE^2))",
            }
        )

    for metadata_path in METADATA_FILES:
        update_metadata(metadata_path, neff_map)

    with SUMMARY_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "factor_id",
                "standard_file",
                "n_original_label",
                "neff_mean",
                "neff_median",
                "total_rows",
                "rows_used_for_neff",
                "maf_window",
                "formula",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
