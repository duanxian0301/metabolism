from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
REGION_FILE = ROOT / "results" / "08_coloc_pwcoco_loci_lipid8_F2_AD" / "lipid8_F2_AD_regions_500kb.tsv"
COLOC_FILE = ROOT / "results" / "09_coloc_lipid8_F2_AD" / "coloc_lipid8_F2_AD_regions.tsv"
PWC_FILE = ROOT / "results" / "10_pwcoco_lipid8_F2_AD" / "pwcoco_lipid8_F2_AD_summary.tsv"
OUT_DIR = ROOT / "results" / "10_pwcoco_lipid8_F2_AD"


def h4_class(h4: float | None) -> str:
    if pd.isna(h4):
        return "missing"
    if h4 >= 0.8:
        return "strong_H4"
    if h4 >= 0.5:
        return "moderate_H4"
    return "low_H4"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    regions = pd.read_csv(REGION_FILE, sep="\t")
    coloc = pd.read_csv(COLOC_FILE, sep="\t")
    pwc = pd.read_csv(PWC_FILE, sep="\t")

    pwc["H4"] = pd.to_numeric(pwc["H4"], errors="coerce")
    pwc["H3"] = pd.to_numeric(pwc["H3"], errors="coerce")
    pwc["conditional"] = pwc["result_type"] != "unconditioned"
    pwc["h4_class"] = pwc["H4"].apply(h4_class)

    counts = (
        pwc.groupby(["result_type", "h4_class"], dropna=False)
        .size()
        .reset_index(name="n_rows")
        .sort_values(["result_type", "h4_class"])
    )
    counts.to_csv(OUT_DIR / "pwcoco_lipid8_F2_AD_result_counts.tsv", sep="\t", index=False)

    best_all = pwc.sort_values(["region_id", "H4"], ascending=[True, False]).groupby("region_id", as_index=False).head(1)
    best_cond = (
        pwc[pwc["conditional"]]
        .sort_values(["region_id", "H4"], ascending=[True, False])
        .groupby("region_id", as_index=False)
        .head(1)
    )
    best_uncond = pwc[pwc["result_type"] == "unconditioned"].copy()

    keep = ["region_id", "SNP1", "SNP2", "H3", "H4", "result_type", "h4_class"]
    best_all = best_all[keep].rename(
        columns={c: f"pwcoco_best_{c}" for c in keep if c != "region_id"}
    )
    best_cond = best_cond[keep].rename(
        columns={c: f"pwcoco_best_cond_{c}" for c in keep if c != "region_id"}
    )
    best_uncond = best_uncond[keep].rename(
        columns={c: f"pwcoco_uncond_{c}" for c in keep if c != "region_id"}
    )

    coloc_keep = [
        "region_id",
        "nsnps",
        "PP.H3",
        "PP.H4",
        "status",
        "coloc_class",
    ]
    if "coloc_class" not in coloc.columns:
        coloc["coloc_class"] = coloc["PP.H4"].apply(h4_class)
    merged = (
        regions.merge(coloc[coloc_keep], on="region_id", how="left")
        .merge(best_uncond, on="region_id", how="left")
        .merge(best_cond, on="region_id", how="left")
        .merge(best_all, on="region_id", how="left")
    )

    merged["pwcoco_region_class"] = merged["pwcoco_best_cond_H4"].apply(h4_class)
    merged["priority_shared_signal"] = (
        merged["PP.H4"].fillna(0).ge(0.5) | merged["pwcoco_best_cond_H4"].fillna(0).ge(0.8)
    )
    merged.to_csv(OUT_DIR / "coloc_pwcoco_lipid8_F2_AD_region_integrated.tsv", sep="\t", index=False)

    priority = merged[merged["priority_shared_signal"]].copy()
    conj_col = "best_conjfdr" if "best_conjfdr" in priority.columns else None
    sort_cols = ["pwcoco_best_cond_H4", "PP.H4"]
    sort_ascending = [False, False]
    if conj_col:
        sort_cols.append(conj_col)
        sort_ascending.append(True)
    priority = priority.sort_values(sort_cols, ascending=sort_ascending)
    priority.to_csv(OUT_DIR / "coloc_pwcoco_lipid8_F2_AD_priority_regions.tsv", sep="\t", index=False)

    print(
        {
            "n_regions": int(len(merged)),
            "n_coloc_h4_ge_0_5": int(merged["PP.H4"].fillna(0).ge(0.5).sum()),
            "n_pwcoco_cond_h4_ge_0_8": int(merged["pwcoco_best_cond_H4"].fillna(0).ge(0.8).sum()),
            "n_priority_regions": int(len(priority)),
        }
    )


if __name__ == "__main__":
    main()
