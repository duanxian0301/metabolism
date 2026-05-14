from __future__ import annotations

import csv
import math
from pathlib import Path


PROJECT_ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
RESULTS_DIR = PROJECT_ROOT / "results" / "00_input_lock"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

METABOLIC_ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
NEFF_SUMMARY = METABOLIC_ROOT / "factor_standard_txt_neff_summary.tsv"
FINALIZATION_SUMMARY = METABOLIC_ROOT / "factor_txt_gz_finalization_summary.tsv"
CHR_BP_SUMMARY = METABOLIC_ROOT / "factor_standard_txt_chr_bp_finalfill_summary.tsv"

DISEASE_TRAITS = {
    "AD": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
    "PD": Path(r"D:\文章\4NDD\NDDGWAS\PD.txt"),
    "LBD": Path(r"D:\文章\4NDD\NDDGWAS\LBD.txt"),
}

REQUIRED_COLUMNS = ["SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_path(value: str) -> Path:
    return Path(value.replace("/", "\\"))


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


def inspect_tsv(path: Path) -> dict[str, object]:
    size_bytes = path.stat().st_size if path.exists() else None
    header: list[str] = []
    row_count = 0
    missing_required: list[str] = []
    n_values: set[str] = set()
    min_chr: int | None = None
    max_chr: int | None = None
    missing_chr_bp = 0
    bad_se = 0
    bad_p = 0

    if not path.exists():
        return {
            "exists": False,
            "size_bytes": size_bytes,
            "row_count": "",
            "columns": "",
            "missing_required": ",".join(REQUIRED_COLUMNS),
            "n_unique_values_sample": "",
            "min_chr": "",
            "max_chr": "",
            "missing_chr_bp": "",
            "bad_se": "",
            "bad_p": "",
        }

    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        missing_required = [col for col in REQUIRED_COLUMNS if col not in header]

        for row in reader:
            row_count += 1
            if len(n_values) < 10 and row.get("N"):
                n_values.add(str(row["N"]))

            chrom = to_float(row.get("CHR"))
            bp = to_float(row.get("BP"))
            if chrom is None or bp is None:
                missing_chr_bp += 1
            else:
                chrom_i = int(chrom)
                min_chr = chrom_i if min_chr is None else min(min_chr, chrom_i)
                max_chr = chrom_i if max_chr is None else max(max_chr, chrom_i)

            se = to_float(row.get("SE"))
            if se is None or se <= 0:
                bad_se += 1

            p = to_float(row.get("P"))
            if p is None or p <= 0 or p > 1:
                bad_p += 1

    return {
        "exists": True,
        "size_bytes": size_bytes,
        "row_count": row_count,
        "columns": ",".join(header),
        "missing_required": ",".join(missing_required),
        "n_unique_values_sample": ",".join(sorted(n_values)),
        "min_chr": "" if min_chr is None else min_chr,
        "max_chr": "" if max_chr is None else max_chr,
        "missing_chr_bp": missing_chr_bp,
        "bad_se": bad_se,
        "bad_p": bad_p,
    }


def main() -> None:
    neff_rows = read_tsv(NEFF_SUMMARY)
    finalization_by_file = {
        str(normalize_path(row["standard_file"])): row for row in read_tsv(FINALIZATION_SUMMARY)
    }
    chrbp_by_file = {
        str(normalize_path(row["file"])): row for row in read_tsv(CHR_BP_SUMMARY)
    }

    factor_manifest: list[dict[str, object]] = []
    for row in neff_rows:
        factor_id = row["factor_id"]
        path = normalize_path(row["standard_file"])
        path_key = str(path)
        inspect = inspect_tsv(path)
        final_row = finalization_by_file.get(path_key, {})
        chrbp_row = chrbp_by_file.get(path_key, {})

        factor_manifest.append(
            {
                "trait": factor_id,
                "group": "metabolic_factor",
                "standard_file": str(path),
                "n_original_label": row.get("n_original_label", ""),
                "neff_mean": row.get("neff_mean", ""),
                "neff_median": row.get("neff_median", ""),
                "neff_formula": row.get("formula", ""),
                "rows_after_finalization": final_row.get("txt_rows_after", ""),
                "rows_removed_missing_chr_bp": final_row.get("rows_removed", ""),
                "pct_chr_bp_filled_total": chrbp_row.get("pct_filled_total", ""),
                **inspect,
            }
        )

    disease_manifest: list[dict[str, object]] = []
    for trait, path in DISEASE_TRAITS.items():
        inspect = inspect_tsv(path)
        disease_manifest.append(
            {
                "trait": trait,
                "group": "ndd_trait",
                "standard_file": str(path),
                "n_original_label": "",
                "neff_mean": "",
                "neff_median": "",
                "neff_formula": "",
                "rows_after_finalization": "",
                "rows_removed_missing_chr_bp": "",
                "pct_chr_bp_filled_total": "",
                **inspect,
            }
        )

    all_rows = factor_manifest + disease_manifest
    manifest_path = RESULTS_DIR / "metabolic_ad_pdlbd_input_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)

    blockers = []
    for row in all_rows:
        if not row["exists"]:
            blockers.append(f"{row['trait']}: missing file")
        if row["missing_required"]:
            blockers.append(f"{row['trait']}: missing required columns {row['missing_required']}")
        if row["missing_chr_bp"] not in ("", 0):
            blockers.append(f"{row['trait']}: {row['missing_chr_bp']} rows missing CHR/BP")
        if row["bad_se"] not in ("", 0):
            blockers.append(f"{row['trait']}: {row['bad_se']} rows with invalid SE")
        if row["bad_p"] not in ("", 0):
            blockers.append(f"{row['trait']}: {row['bad_p']} rows with invalid P")

    blockers_path = RESULTS_DIR / "metabolic_ad_pdlbd_input_blockers.txt"
    blockers_path.write_text("\n".join(blockers) + ("\n" if blockers else ""), encoding="utf-8")

    print(f"Wrote manifest: {manifest_path}")
    if blockers:
        print(f"Found {len(blockers)} blockers. See: {blockers_path}")
        for item in blockers[:20]:
            print("BLOCKER:", item)
        raise SystemExit(1)
    print("No input blockers found.")


if __name__ == "__main__":
    main()
