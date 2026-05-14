from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
COMPACT_REVIEW_DIR = ROOT / "nonlipid_module_from_full_manifest" / "compact_panel_review"
COMPACT_ESEM_DIR = ROOT / "step17_efa_esem_nonlipid_module_compact15"
OUT_DIR = COMPACT_REVIEW_DIR / "ultra_pure_3factor_review"
LOCAL_SELECTION_DIR = Path(__file__).resolve().parents[1] / "inputs" / "selection"

KEEP_TRAITS = {
    "Val",
    "Leu",
    "Phe",
    "Acetoacetate",
    "Acetate",
    "bOHbutyrate",
    "Lactate",
    "Glucose",
}

ULTRA_PURE_FACTOR = {
    "Val": "F2_amino_acid_axis",
    "Leu": "F2_amino_acid_axis",
    "Phe": "F2_amino_acid_axis",
    "Acetoacetate": "F1_ketone_axis",
    "bOHbutyrate": "F1_ketone_axis",
    "Acetate": "F3_energy_bridge_axis",
    "Lactate": "F3_energy_bridge_axis",
    "Glucose": "F3_energy_bridge_axis",
}

SELECTION_REASON = {
    "Val": "Retained as the strongest and most stable BCAA anchor across ALL and ODD with a very large loading gap on the amino-acid factor.",
    "Leu": "Retained as the second core BCAA anchor with highly stable primary loading on the amino-acid factor.",
    "Phe": "Retained as the cleanest aromatic amino-acid marker and the best third indicator for the amino-acid factor in the high-fit panels.",
    "Tyr": "Dropped because tyrosine produced slightly weaker fit than phenylalanine in the best-performing compact-to-ultrapure search panels.",
    "His": "Dropped because histidine cross-loaded on the ketone and bridge factors and did not improve fit in the strongest ultrapure candidates.",
    "Gln": "Dropped because glutamine showed unstable primary-factor assignment between ALL and ODD and reduced model purity.",
    "Acetoacetate": "Retained as the clearest ketone-body anchor with high and stable primary loading on the ketone factor.",
    "bOHbutyrate": "Retained as the complementary ketone-body anchor with strong and stable ketone-factor loading.",
    "Acetate": "Retained as the defining marker for the energy/inflammation bridge factor and as the key separator from the core ketone factor.",
    "Lactate": "Retained as the strongest glycolytic bridge marker and a stable negative pole of the energy-bridge factor.",
    "Glucose": "Retained because the best-performing 8-trait ultrapure panels consistently included glucose as the third bridge-factor indicator.",
    "Citrate": "Dropped because citrate had weak loadings and was not needed in the best-fitting ultrapure panels.",
    "GlycA": "Dropped because GlycA had unstable primary-factor assignment across ALL and ODD and worsened fit in the best ultrapure panels.",
    "Albumin": "Dropped because albumin was weakly loading and less effective than glucose in the top-scoring compact subset search.",
    "Creatinine": "Dropped because creatinine had very small loadings and did not help define a clean ultrapure factor.",
}


def summarize_primary_factors() -> pd.DataFrame:
    all_df = pd.read_csv(COMPACT_ESEM_DIR / "ALL_3factor_loadings.tsv", sep="\t")
    odd_df = pd.read_csv(COMPACT_ESEM_DIR / "ODD_3factor_loadings.tsv", sep="\t")

    def build(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
        work = df.copy()
        work["abs_loading"] = work["std.all"].abs()
        rows = []
        for trait, sub in work.groupby("trait"):
            sub = sub.sort_values("abs_loading", ascending=False)
            rows.append(
                {
                    "trait_code": trait,
                    f"primary_factor_{suffix}": sub.iloc[0]["factor"],
                    f"top_abs_{suffix}": sub.iloc[0]["abs_loading"],
                    f"second_abs_{suffix}": sub.iloc[1]["abs_loading"],
                    f"gap_{suffix}": sub.iloc[0]["abs_loading"] - sub.iloc[1]["abs_loading"],
                    f"f1_{suffix}": float(sub.loc[sub["factor"] == "f1", "std.all"].iloc[0]),
                    f"f2_{suffix}": float(sub.loc[sub["factor"] == "f2", "std.all"].iloc[0]),
                    f"f3_{suffix}": float(sub.loc[sub["factor"] == "f3", "std.all"].iloc[0]),
                }
            )
        return pd.DataFrame(rows)

    merged = build(all_df, "all").merge(build(odd_df, "odd"), on="trait_code")
    merged["same_primary"] = merged["primary_factor_all"] == merged["primary_factor_odd"]
    merged["mean_top_abs"] = merged[["top_abs_all", "top_abs_odd"]].mean(axis=1)
    merged["mean_gap"] = merged[["gap_all", "gap_odd"]].mean(axis=1)
    return merged


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_SELECTION_DIR.mkdir(parents=True, exist_ok=True)

    compact = pd.read_csv(COMPACT_REVIEW_DIR / "nonlipid_module_compact_kept.tsv", sep="\t")
    metrics = summarize_primary_factors()
    review = compact.merge(metrics, on="trait_code", how="left")

    review["selection_status"] = review["trait_code"].map(
        lambda x: "kept" if x in KEEP_TRAITS else "dropped"
    )
    review["ultra_pure_factor"] = review["trait_code"].map(
        lambda x: ULTRA_PURE_FACTOR.get(x, "not_retained")
    )
    review["selection_reason"] = review["trait_code"].map(SELECTION_REASON)

    review = review[
        [
            "study_accession",
            "trait_code",
            "biomarker_name",
            "group",
            "proposed_module",
            "marker_role",
            "selection_status",
            "ultra_pure_factor",
            "same_primary",
            "mean_top_abs",
            "mean_gap",
            "primary_factor_all",
            "primary_factor_odd",
            "f1_all",
            "f2_all",
            "f3_all",
            "f1_odd",
            "f2_odd",
            "f3_odd",
            "selection_reason",
            "sumstats_file",
        ]
    ].sort_values(["selection_status", "ultra_pure_factor", "trait_code"], ascending=[False, True, True])

    kept = review[review["selection_status"] == "kept"].copy()
    dropped = review[review["selection_status"] == "dropped"].copy()

    for target_dir in [OUT_DIR, LOCAL_SELECTION_DIR]:
        review.to_csv(target_dir / "nonlipid_module_ultrapure3_review.tsv", sep="\t", index=False)
        kept.to_csv(target_dir / "nonlipid_module_ultrapure3_kept.tsv", sep="\t", index=False)
        dropped.to_csv(target_dir / "nonlipid_module_ultrapure3_dropped.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
