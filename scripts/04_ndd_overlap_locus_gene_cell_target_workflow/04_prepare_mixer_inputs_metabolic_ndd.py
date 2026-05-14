from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path


PROJECT_ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
OUT_DIR = PROJECT_ROOT / "work" / "04_mixer_inputs" / "inputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR = PROJECT_ROOT / "results" / "04_mixer_inputs"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CLEAN_DIR = PROJECT_ROOT / "work" / "clean_factor_inputs"

TRAITS = [
    ("lipid8_F1", CLEAN_DIR / "lipid8_F1_clean.txt", 276215.0),
    ("lipid8_F2", CLEAN_DIR / "lipid8_F2_clean.txt", 532584.0),
    ("lipid8_F3", CLEAN_DIR / "lipid8_F3_clean.txt", 236386.0),
    ("nonlipid8_F1", CLEAN_DIR / "nonlipid8_F1_clean.txt", 101905.0),
    ("nonlipid8_F2", CLEAN_DIR / "nonlipid8_F2_clean.txt", 123638.0),
    ("nonlipid8_F3", CLEAN_DIR / "nonlipid8_F3_clean.txt", 34657.0),
    ("AD", Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"), 4.0 / (1.0 / (39106.0 + 46828.0) + 1.0 / 401577.0)),
    ("PD", Path(r"D:\文章\4NDD\NDDGWAS\PD.txt"), 4.0 / (1.0 / (63555.0 + 17700.0) + 1.0 / 1746386.0)),
    ("LBD", Path(r"D:\文章\4NDD\NDDGWAS\LBD.txt"), 4.0 / (1.0 / 2591.0 + 1.0 / 4027.0)),
]


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


def parse_chr(value: str | None) -> int | None:
    if value is None:
        return None
    text = value.strip().upper().replace("CHR", "")
    try:
        chrom = int(float(text))
    except ValueError:
        return None
    return chrom if 1 <= chrom <= 22 else None


def transform_trait(trait: str, source: Path, n_value: float) -> dict[str, object]:
    target = OUT_DIR / f"{trait}.sumstats.gz"
    kept = 0
    skipped = 0

    with source.open("r", encoding="utf-8-sig", errors="replace", newline="") as fin, gzip.open(
        target, "wt", encoding="utf-8", newline=""
    ) as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "N", "Z"])

        for row in reader:
            snp = (row.get("SNP") or "").strip()
            chrom = parse_chr(row.get("CHR"))
            bp = to_float(row.get("BP"))
            a1 = (row.get("A1") or "").strip().upper()
            a2 = (row.get("A2") or "").strip().upper()
            beta = to_float(row.get("BETA"))
            se = to_float(row.get("SE"))
            if not snp or chrom is None or bp is None or not a1 or not a2 or beta is None or se is None or se <= 0:
                skipped += 1
                continue
            z = beta / se
            if not math.isfinite(z):
                skipped += 1
                continue
            writer.writerow([snp, chrom, int(bp), a1, a2, f"{n_value:.6f}", f"{z:.10f}"])
            kept += 1

    return {
        "trait": trait,
        "source_file": str(source),
        "mixer_input": str(target),
        "N_used": f"{n_value:.6f}",
        "rows_kept": kept,
        "rows_skipped": skipped,
    }


def main() -> None:
    rows = [transform_trait(*trait_cfg) for trait_cfg in TRAITS]
    manifest = RESULTS_DIR / "mixer_input_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote MiXeR inputs to: {OUT_DIR}")
    print(f"Wrote manifest: {manifest}")
    for row in rows:
        print(row["trait"], "kept=", row["rows_kept"], "skipped=", row["rows_skipped"])


if __name__ == "__main__":
    main()
