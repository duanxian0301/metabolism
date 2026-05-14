from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
PLEIO_DATA = Path(r"D:\pleioFDR\pleiofdr-master\data")

INPUTS = {
    "lipid8_F2": ROOT / "work" / "clean_factor_inputs" / "lipid8_F2_clean.txt",
    "AD_metabolic": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
}


def as_float(value: str) -> float:
    if value is None or value == "":
        return float("nan")
    return float(value)


def write_fdr_txt(name: str, src: Path, out_dir: Path) -> dict[str, int | str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    dst = out_dir / f"{name}_fdr.txt"

    seen = set()
    total = written = dropped = duplicate = 0

    with src.open("r", newline="") as fin, dst.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        required = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{src} missing required columns: {sorted(missing)}")

        writer = csv.DictWriter(
            fout,
            fieldnames=["SNP", "CHR", "BP", "A1", "A2", "Z", "PVAL", "N"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            total += 1
            snp = row["SNP"]
            if not snp or snp in seen:
                duplicate += 1
                continue
            seen.add(snp)

            try:
                beta = as_float(row["BETA"])
                se = as_float(row["SE"])
                p = as_float(row["P"])
                n = as_float(row["N"])
            except ValueError:
                dropped += 1
                continue

            if not (
                math.isfinite(beta)
                and math.isfinite(se)
                and se > 0
                and math.isfinite(p)
                and 0 < p <= 1
                and math.isfinite(n)
                and n > 0
            ):
                dropped += 1
                continue

            writer.writerow(
                {
                    "SNP": snp,
                    "CHR": row["CHR"],
                    "BP": row["BP"],
                    "A1": row["A1"],
                    "A2": row["A2"],
                    "Z": beta / se,
                    "PVAL": p,
                    "N": n,
                }
            )
            written += 1

    return {
        "trait": name,
        "source": str(src),
        "output": str(dst),
        "total": total,
        "written": written,
        "dropped": dropped,
        "duplicate": duplicate,
    }


def main() -> None:
    reports = [write_fdr_txt(name, src, PLEIO_DATA) for name, src in INPUTS.items()]

    report_path = ROOT / "results" / "05_pleiofdr_inputs" / "pleiofdr_input_report.tsv"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=list(reports[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(reports)

    print(f"Wrote report: {report_path}")


if __name__ == "__main__":
    main()
