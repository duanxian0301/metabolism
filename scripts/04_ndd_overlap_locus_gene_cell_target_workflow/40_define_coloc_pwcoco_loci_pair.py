from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
WINDOW_BP = 500_000


def build_sentinel_table(loci: pd.DataFrame, trait: str, disease: str, pair: str) -> pd.DataFrame:
    out = loci.loc[:, ["locusnum", "snpid", "chrnum", "chrpos", "min_conjfdr"]].copy()
    out = out.sort_values(["chrnum", "chrpos", "min_conjfdr"]).reset_index(drop=True)
    out.insert(0, "sentinel_id", [f"S{i:03d}" for i in range(1, len(out) + 1)])
    out["trait"] = trait
    out["disease"] = disease
    out["pair"] = pair
    out["region_start_500kb"] = (out["chrpos"].astype(int) - WINDOW_BP).clip(lower=1)
    out["region_end_500kb"] = out["chrpos"].astype(int) + WINDOW_BP
    return out


def build_region_table(sentinels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for chrnum, sub in sentinels.groupby("chrnum", sort=True):
        sub = sub.sort_values("chrpos").reset_index(drop=True)
        current = None
        for _, row in sub.iterrows():
            start = int(row["region_start_500kb"])
            end = int(row["region_end_500kb"])
            if current is None or start > current["region_end"]:
                if current is not None:
                    rows.append(current)
                current = {
                    "chrnum": int(chrnum),
                    "region_start": start,
                    "region_end": end,
                    "sentinel_ids": [row["sentinel_id"]],
                    "sentinel_snps": [row["snpid"]],
                    "best_conjfdr": float(row["min_conjfdr"]),
                }
            else:
                current["region_end"] = max(current["region_end"], end)
                current["sentinel_ids"].append(row["sentinel_id"])
                current["sentinel_snps"].append(row["snpid"])
                current["best_conjfdr"] = min(current["best_conjfdr"], float(row["min_conjfdr"]))
        if current is not None:
            rows.append(current)
    out = pd.DataFrame(rows)
    out.insert(0, "region_id", [f"R{i:03d}" for i in range(1, len(out) + 1)])
    out["n_sentinels"] = out["sentinel_ids"].apply(len)
    out["sentinel_ids"] = out["sentinel_ids"].apply(lambda x: ";".join(x))
    out["sentinel_snps"] = out["sentinel_snps"].apply(lambda x: ";".join(x))
    out["region_width_bp"] = out["region_end"] - out["region_start"] + 1
    out["priority"] = out["n_sentinels"].apply(lambda n: "PWCoCo_priority" if n >= 2 else "coloc_first")
    return out.loc[
        :,
        [
            "region_id",
            "chrnum",
            "region_start",
            "region_end",
            "region_width_bp",
            "n_sentinels",
            "sentinel_ids",
            "sentinel_snps",
            "best_conjfdr",
            "priority",
        ],
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", required=True)
    parser.add_argument("--disease", required=True)
    parser.add_argument("--pleio-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    pair = f"{args.trait}_{args.disease}"
    pleio_loci = Path(args.pleio_dir) / f"{pair}_conjfdr_0.05_loci.csv"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    loci = pd.read_csv(pleio_loci)
    sentinels = build_sentinel_table(loci, args.trait, args.disease, pair)
    regions = build_region_table(sentinels)
    sentinels.to_csv(out_dir / f"{pair}_conjfdr_sentinels.tsv", sep="\t", index=False)
    regions.to_csv(out_dir / f"{pair}_regions_500kb.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "pair": pair,
                "n_conjfdr_locus_rows": len(loci),
                "n_sentinels": len(sentinels),
                "n_regions_500kb": len(regions),
                "n_pwcoco_priority_regions": int((regions["n_sentinels"] >= 2).sum()),
            }
        ]
    ).to_csv(out_dir / f"{pair}_locus_region_summary.tsv", sep="\t", index=False)
    print(out_dir)


if __name__ == "__main__":
    main()
