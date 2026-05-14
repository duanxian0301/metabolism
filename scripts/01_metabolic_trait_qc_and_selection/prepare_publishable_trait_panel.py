from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
ANNOT_PATH = BASE / "balanced_selection_3to4f" / "traits_balanced_annotated.tsv"
LOADINGS_PATH = BASE / "step4_esem_target_balanced_56" / "ALL_4factor_loadings.tsv"
RG_PATH = (
    BASE
    / "step2_efa_balanced_3to4f_minres"
    / "balanced_56trait_ODD_rg_matrix_smoothed.csv"
)
OUT_DIR = BASE / "publishable_compact_trait_panel_review"


KEEP_SPECS = {
    "VLDL_CE": {
        "module": "F1_VLDL_CE_structural_lipids",
        "role": "core_marker",
        "reason": "strongest VLDL/CE transport marker with very high factor specificity",
    },
    "XS_VLDL_FC": {
        "module": "F1_VLDL_CE_structural_lipids",
        "role": "core_marker",
        "reason": "captures very-small-VLDL free cholesterol while replacing the highly redundant XS_VLDL_L",
    },
    "S_VLDL_L": {
        "module": "F1_VLDL_CE_structural_lipids",
        "role": "core_marker",
        "reason": "retained as a compact total-lipid marker for the small-VLDL compartment",
    },
    "Total_L": {
        "module": "F1_VLDL_CE_structural_lipids",
        "role": "core_marker",
        "reason": "broad total-lipid transport marker retained for global structural lipid burden",
    },
    "Total_CE": {
        "module": "F1_VLDL_CE_structural_lipids",
        "role": "core_marker",
        "reason": "broad esterified-cholesterol marker retained for non-subclass-specific CE signal",
    },
    "LA": {
        "module": "F1_VLDL_CE_structural_lipids",
        "role": "core_marker",
        "reason": "chosen as the most interpretable omega-6/PUFA structural fatty-acid representative",
    },
    "ApoA1": {
        "module": "F2_HDL_CE_phospholipids",
        "role": "core_marker",
        "reason": "canonical HDL/apolipoprotein marker with very strong loading and good separation",
    },
    "M_HDL_CE": {
        "module": "F2_HDL_CE_phospholipids",
        "role": "core_marker",
        "reason": "retained as the clearest HDL cholesteryl-ester subclass marker",
    },
    "S_HDL_C": {
        "module": "F2_HDL_CE_phospholipids",
        "role": "core_marker",
        "reason": "selected over S_HDL_CE to reduce redundancy while keeping small-HDL cholesterol content",
    },
    "Phosphatidylc": {
        "module": "F2_HDL_CE_phospholipids",
        "role": "core_marker",
        "reason": "represents the phospholipid/choline cluster after collapsing the near-duplicate phospholipid trio",
    },
    "HDL_CE": {
        "module": "F2_HDL_CE_phospholipids",
        "role": "core_marker",
        "reason": "global HDL cholesteryl-ester marker retained for broader HDL esterification signal",
    },
    "S_HDL_L": {
        "module": "F2_HDL_CE_phospholipids",
        "role": "core_marker",
        "reason": "retained as a total small-HDL lipid marker complementary to HDL cholesterol composition",
    },
    "VLDL_size": {
        "module": "F3_TG_particle_size_BCAA",
        "role": "core_marker",
        "reason": "strong particle-size marker anchoring the TG/size axis",
    },
    "LDL_size": {
        "module": "F3_TG_particle_size_BCAA",
        "role": "core_marker",
        "reason": "kept to preserve the opposite size direction on the particle-size dimension",
    },
    "M_HDL_TG": {
        "module": "F3_TG_particle_size_BCAA",
        "role": "core_marker",
        "reason": "strong HDL-TG marker for the atherogenic TG enrichment axis",
    },
    "S_HDL_TG": {
        "module": "F3_TG_particle_size_BCAA",
        "role": "core_marker",
        "reason": "small-HDL TG marker retained to stabilize the TG-rich HDL subdimension",
    },
    "S_VLDL_TG": {
        "module": "F3_TG_particle_size_BCAA",
        "role": "core_marker",
        "reason": "chosen as the compact VLDL-TG representative after collapsing the highly correlated TG trio",
    },
    "MUFA": {
        "module": "F3_TG_particle_size_BCAA",
        "role": "core_marker",
        "reason": "retained as the clearest fatty-acid marker aligned with the TG/IR-like axis",
    },
    "Val": {
        "module": "F3_TG_particle_size_BCAA",
        "role": "core_marker",
        "reason": "selected as the BCAA representative after collapsing the Ile/Leu/Val redundancy cluster",
    },
    "GlycA": {
        "module": "F4_HDL_size_inflammation_ketone_support",
        "role": "support_marker",
        "reason": "retained to preserve the inflammatory component that repeatedly appears on the fourth dimension",
    },
    "bOHbutyrate": {
        "module": "F4_HDL_size_inflammation_ketone_support",
        "role": "support_marker",
        "reason": "best ketone-body marker for the fourth dimension despite moderate loading magnitude",
    },
    "Acetate": {
        "module": "F4_HDL_size_inflammation_ketone_support",
        "role": "support_marker",
        "reason": "kept as a secondary ketone/acetate support trait for the weaker fourth factor",
    },
    "Gln": {
        "module": "F4_HDL_size_inflammation_ketone_support",
        "role": "support_marker",
        "reason": "retained as a non-lipoprotein support marker linked to the weaker fourth dimension",
    },
    "HDL_size": {
        "module": "F4_HDL_size_inflammation_ketone_support",
        "role": "support_marker",
        "reason": "anchors the HDL-size component that differentiates the fourth dimension from the main HDL factor",
    },
}


DROPPED_REASON_OVERRIDES = {
    "Cholines": "dropped as part of the phospholipid/choline redundancy cluster; Phosphatidylc retained",
    "Phosphoglyc": "dropped as part of the phospholipid/choline redundancy cluster; Phosphatidylc retained",
    "XS_VLDL_L": "dropped for extreme redundancy with XS_VLDL_FC",
    "S_VLDL_FC": "dropped for extreme redundancy with VLDL_CE/S_VLDL_CE and to avoid repeated small-VLDL FC markers",
    "XS_VLDL_CE": "dropped because VLDL_CE and XS_VLDL_FC already capture the same structural-lipid block with less redundancy",
    "S_VLDL_CE": "dropped because VLDL_CE retained as the cleaner CE anchor for this block",
    "L_HDL_C": "dropped for extreme redundancy with HDL_size and overlap with the main HDL composition factor",
    "S_HDL_CE": "dropped for strong redundancy with S_HDL_C and heavy cross-loading onto factor 4",
    "Ile": "dropped as part of the BCAA redundancy cluster; Val retained",
    "Leu": "dropped as part of the BCAA redundancy cluster; Val retained",
    "M_LDL_TG": "dropped as part of the TG-rich lipoprotein redundancy cluster; S_VLDL_TG retained",
    "S_LDL_TG": "dropped as part of the TG-rich lipoprotein redundancy cluster; S_VLDL_TG retained",
    "XS_VLDL_TG": "dropped as part of the TG-rich lipoprotein redundancy cluster; S_VLDL_TG retained",
    "SFA": "dropped for near-redundancy with Total_FA and broad cross-loading",
    "Omega_6": "dropped for near-redundancy with LA",
    "PUFA": "dropped for near-redundancy with LA/Omega_6",
    "Acetoacetate": "dropped because of weak loading and minimal factor separation",
    "Sphingomyelins": "dropped for strong cross-loading between factors 1 and 2",
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    annot = pd.read_csv(ANNOT_PATH, sep="\t")
    loadings = pd.read_csv(LOADINGS_PATH, sep="\t")
    rg = pd.read_csv(RG_PATH, index_col=0)

    wide = loadings.pivot_table(
        index=["trait", "biomarker_name", "group"],
        columns="factor",
        values="std.all",
    )
    abswide = wide.abs()

    stats_rows = []
    for idx, row in abswide.iterrows():
        ordered = row.sort_values(ascending=False)
        stats_rows.append(
            {
                "trait_code": idx[0],
                "biomarker_name": idx[1],
                "group": idx[2],
                "primary_factor": ordered.index[0],
                "top_loading_abs": float(ordered.iloc[0]),
                "second_loading_abs": float(ordered.iloc[1]),
                "loading_gap": float(ordered.iloc[0] - ordered.iloc[1]),
            }
        )
    stats = pd.DataFrame(stats_rows)

    merged = annot.merge(stats, on=["trait_code", "biomarker_name", "group"], how="left")
    merged["selection_status"] = np.where(
        merged["trait_code"].isin(KEEP_SPECS), "keep", "drop"
    )
    merged["proposed_module"] = merged["trait_code"].map(
        {k: v["module"] for k, v in KEEP_SPECS.items()}
    )
    merged["marker_role"] = merged["trait_code"].map(
        {k: v["role"] for k, v in KEEP_SPECS.items()}
    )
    merged["selection_reason"] = merged["trait_code"].map(
        {k: v["reason"] for k, v in KEEP_SPECS.items()}
    )

    merged["selection_reason"] = merged["selection_reason"].fillna(
        merged["trait_code"].map(DROPPED_REASON_OVERRIDES)
    )

    default_drop_reason = (
        "dropped because of weaker factor purity, lower interpretability, or redundancy relative to retained markers"
    )
    merged["selection_reason"] = merged["selection_reason"].fillna(default_drop_reason)

    high_rg_links = []
    trait_codes = merged["trait_code"].tolist()
    for trait in trait_codes:
        if trait not in rg.index:
            high_rg_links.append("")
            continue
        vals = rg.loc[trait, trait_codes].drop(labels=[trait], errors="ignore")
        vals = vals[vals.abs() > 0.95].sort_values(key=lambda s: s.abs(), ascending=False)
        if vals.empty:
            high_rg_links.append("")
        else:
            links = [f"{idx}({val:.3f})" for idx, val in vals.items()]
            high_rg_links.append("; ".join(links))
    merged["high_rg_links_gt_0p95"] = high_rg_links

    kept = merged[merged["selection_status"] == "keep"].copy()
    dropped = merged[merged["selection_status"] == "drop"].copy()

    kept = kept.sort_values(
        ["proposed_module", "marker_role", "top_loading_abs"],
        ascending=[True, True, False],
    )
    dropped = dropped.sort_values(["group", "top_loading_abs"], ascending=[True, False])

    summary_lines = [
        "Proposed compact publishable trait panel",
        f"Starting balanced panel: {len(merged)} traits",
        f"Retained for next-step module-aware modeling: {len(kept)} traits",
        f"Dropped at review stage: {len(dropped)} traits",
        "",
        "Retention principles:",
        "1. Preserve 3-4 interpretable metabolic dimensions from the 4-factor target-rotation ESEM.",
        "2. Prefer traits with stronger primary loading and wider loading gap.",
        "3. Collapse extreme |rg| > 0.95 redundancy clusters to one representative marker.",
        "4. Keep a small fourth-dimension support block even though factor 4 is weaker than factors 1-3.",
    ]

    summary_path = OUT_DIR / "selection_summary.txt"
    summary_path.write_text("\n".join(summary_lines), encoding="utf-8")

    review_cols = [
        "study_accession",
        "trait_code",
        "biomarker_name",
        "group",
        "Z",
        "h2_obs",
        "n_independent_leads",
        "primary_factor",
        "top_loading_abs",
        "second_loading_abs",
        "loading_gap",
        "selection_status",
        "proposed_module",
        "marker_role",
        "selection_reason",
        "high_rg_links_gt_0p95",
    ]
    merged[review_cols].sort_values(
        ["selection_status", "proposed_module", "group", "top_loading_abs"],
        ascending=[False, True, True, False],
    ).to_csv(
        OUT_DIR / "publishable_compact_trait_panel_review.tsv",
        sep="\t",
        index=False,
    )
    kept[review_cols].to_csv(
        OUT_DIR / "publishable_compact_trait_panel_kept.tsv",
        sep="\t",
        index=False,
    )
    dropped[review_cols].to_csv(
        OUT_DIR / "publishable_compact_trait_panel_dropped.tsv",
        sep="\t",
        index=False,
    )

    module_counts = (
        kept.groupby(["proposed_module", "marker_role"])
        .size()
        .reset_index(name="n_traits")
        .sort_values(["proposed_module", "marker_role"])
    )
    module_counts.to_csv(
        OUT_DIR / "publishable_compact_trait_panel_module_counts.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
