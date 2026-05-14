from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
PLEIO_LOCI = ROOT / "results" / "06_pleiofdr_lipid8_F2_AD" / "lipid8_F2_AD_conjfdr_0.05_loci.csv"
OUT_DIR = ROOT / "results" / "08_coloc_pwcoco_loci_lipid8_F2_AD"
WINDOW_BP = 500_000


def build_sentinel_table(loci: pd.DataFrame) -> pd.DataFrame:
    out = loci.loc[:, ["locusnum", "snpid", "chrnum", "chrpos", "min_conjfdr"]].copy()
    out = out.sort_values(["chrnum", "chrpos", "min_conjfdr"]).reset_index(drop=True)
    out.insert(0, "sentinel_id", [f"S{i:03d}" for i in range(1, len(out) + 1)])
    out["trait"] = "lipid8_F2"
    out["disease"] = "AD"
    out["pair"] = "lipid8_F2_AD"
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
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    loci = pd.read_csv(PLEIO_LOCI)
    sentinels = build_sentinel_table(loci)
    regions = build_region_table(sentinels)
    sentinels.to_csv(OUT_DIR / "lipid8_F2_AD_conjfdr_sentinels.tsv", sep="\t", index=False)
    regions.to_csv(OUT_DIR / "lipid8_F2_AD_regions_500kb.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "pair": "lipid8_F2_AD",
                "n_conjfdr_locus_rows": len(loci),
                "n_sentinels": len(sentinels),
                "n_regions_500kb": len(regions),
                "n_pwcoco_priority_regions": int((regions["n_sentinels"] >= 2).sum()),
            }
        ]
    ).to_csv(OUT_DIR / "lipid8_F2_AD_locus_region_summary.tsv", sep="\t", index=False)
    print(f"Wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
