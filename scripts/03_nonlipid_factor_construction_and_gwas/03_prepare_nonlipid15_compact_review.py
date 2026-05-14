from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
CANDIDATE_DIR = ROOT / "nonlipid_module_from_full_manifest"
OUT_DIR = CANDIDATE_DIR / "compact_panel_review"
LOCAL_SELECTION_DIR = Path(__file__).resolve().parents[1] / "inputs" / "selection"
RG_PATH = ROOT / "step1_ldsc_results" / "Main_Zgt4_nonproportion_rg_matrix.csv"

KEEP_TRAITS = {
    "Gln",
    "His",
    "Leu",
    "Phe",
    "Tyr",
    "Val",
    "Albumin",
    "Creatinine",
    "Citrate",
    "Glucose",
    "Lactate",
    "GlycA",
    "Acetate",
    "Acetoacetate",
    "bOHbutyrate",
}

COMPACT_MODULE_MAP = {
    "Gln": "amino_acid_core",
    "His": "amino_acid_core",
    "Leu": "amino_acid_core",
    "Phe": "amino_acid_core",
    "Tyr": "amino_acid_core",
    "Val": "amino_acid_core",
    "Albumin": "renal_inflammation_bridge",
    "Creatinine": "renal_inflammation_bridge",
    "GlycA": "renal_inflammation_bridge",
    "Citrate": "energy_metabolism_core",
    "Glucose": "energy_metabolism_core",
    "Lactate": "energy_metabolism_core",
    "Acetate": "energy_metabolism_core",
    "Acetoacetate": "energy_metabolism_core",
    "bOHbutyrate": "energy_metabolism_core",
}

COMPACT_MARKER_ROLE = {
    "Gln": "support_marker",
    "His": "support_marker",
    "Leu": "core_marker",
    "Phe": "core_marker",
    "Tyr": "core_marker",
    "Val": "core_marker",
    "Albumin": "sentinel_marker",
    "Creatinine": "sentinel_marker",
    "GlycA": "bridge_marker",
    "Citrate": "support_marker",
    "Glucose": "core_marker",
    "Lactate": "core_marker",
    "Acetate": "core_marker",
    "Acetoacetate": "core_marker",
    "bOHbutyrate": "core_marker",
}

SELECTION_REASON = {
    "Ala": "Dropped at compact stage because the amino-acid block is already well covered and alanine is less specific than the retained BCAA, aromatic, and glutamine/histidine representatives.",
    "Gln": "Retained as a high-h2 non-BCAA amino-acid marker that broadens the amino-acid module beyond BCAA and aromatic signals.",
    "His": "Retained as a lower-redundancy amino-acid marker that may capture unique non-BCAA variation during EFA/ESEM.",
    "Ile": "Dropped at compact stage because it is nearly redundant with both leucine and valine (|rg| > 0.95) and was not needed once the BCAA cluster remained represented by Leu and Val.",
    "Leu": "Retained as a core BCAA representative with strong h2 and slightly less extreme redundancy than isoleucine.",
    "Phe": "Retained as an aromatic amino-acid core marker and likely anchor for any aromatic/BCAA separation.",
    "Tyr": "Retained as the complementary aromatic amino-acid marker and a possible bridge to the broader amino-acid factor.",
    "Val": "Retained as a core BCAA representative with very high h2 Z and strong signal across the amino-acid cluster.",
    "Albumin": "Retained as a fluid-balance / serum abundance sentinel needed to test whether renal-systemic variation forms a bridge factor.",
    "Creatinine": "Retained as a renal-handling sentinel with high h2 and strong interpretability for a nonlipid bridge factor.",
    "Citrate": "Retained as a low-redundancy central-carbon marker that may contribute unique structure despite modest pairwise correlations.",
    "Glucose": "Retained as the canonical glycemic marker for the energy-metabolism block.",
    "Lactate": "Retained as a glycolytic marker and plausible bridge to inflammatory or fasting-related biology.",
    "GlycA": "Retained as the sole inflammation trait because it may bridge amino-acid, glucose, and renal-systemic variation.",
    "Acetate": "Retained as a ketone/fasting-related marker with distinctive negative correlation to lactate and GlycA.",
    "Acetoacetate": "Retained as a canonical ketone marker and likely core anchor for any ketone-body factor.",
    "Acetone": "Dropped at compact stage because the ketone block is already represented by acetate, acetoacetate, and beta-hydroxybutyrate, while acetone is the most dispensable support marker.",
    "bOHbutyrate": "Retained as a canonical beta-hydroxybutyrate marker with cleaner complementarity to acetate and acetoacetate than acetone.",
}


def load_rg_subset(traits: list[str]) -> pd.DataFrame:
    rg = pd.read_csv(RG_PATH)
    rg.index = rg.columns[1:]
    rg = rg.drop(columns=[rg.columns[0]])
    return rg.loc[traits, traits].astype(float)


def high_rg_links(rg_sub: pd.DataFrame, trait: str, threshold: float = 0.95) -> str:
    hits = []
    for other in rg_sub.columns:
        if other == trait:
            continue
        val = float(rg_sub.loc[trait, other])
        if abs(val) > threshold:
            hits.append(f"{other}:{val:.3f}")
    return "; ".join(hits)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_SELECTION_DIR.mkdir(parents=True, exist_ok=True)

    candidate = pd.read_csv(CANDIDATE_DIR / "nonlipid_module_candidates.tsv", sep="\t")
    rg_sub = load_rg_subset(candidate["trait_code"].tolist())

    review = candidate.copy()
    review["selection_status"] = review["trait_code"].map(
        lambda x: "kept" if x in KEEP_TRAITS else "dropped"
    )
    review["proposed_module"] = review["trait_code"].map(
        lambda x: COMPACT_MODULE_MAP.get(x, "not_retained")
    )
    review["marker_role"] = review["trait_code"].map(
        lambda x: COMPACT_MARKER_ROLE.get(x, review.loc[review["trait_code"] == x, "marker_role"].iloc[0])
    )
    review["selection_reason"] = review["trait_code"].map(SELECTION_REASON)
    review["high_rg_links_gt_0p95"] = review["trait_code"].map(
        lambda x: high_rg_links(rg_sub, x, threshold=0.95)
    )

    review = review[
        [
            "study_accession",
            "trait_code",
            "biomarker_name",
            "group",
            "module_name",
            "proposed_submodule",
            "candidate_bucket",
            "marker_role",
            "priority_for_compact_review",
            "h2",
            "intercept",
            "Z",
            "selected_for_efa",
            "max_abs_rg_with_any_trait",
            "n_independent_leads",
            "sumstats_file",
            "selection_status",
            "proposed_module",
            "selection_reason",
            "high_rg_links_gt_0p95",
        ]
    ].sort_values(["selection_status", "proposed_module", "group", "trait_code"], ascending=[False, True, True, True])

    kept = review[review["selection_status"] == "kept"].copy()
    dropped = review[review["selection_status"] == "dropped"].copy()

    group_counts = (
        kept.groupby(["proposed_module", "group"])
        .size()
        .reset_index(name="n_traits")
        .sort_values(["proposed_module", "n_traits", "group"], ascending=[True, False, True])
    )

    summary_lines = [
        "Nonlipid compact review",
        f"Candidate traits reviewed: {len(review)}",
        f"Compact traits retained: {len(kept)}",
        f"Compact traits dropped: {len(dropped)}",
        "Compact strategy: preserve one publishable nonlipid panel while trimming only obvious redundancy and lower-priority support markers.",
        "Dropped traits: " + ", ".join(dropped["trait_code"].tolist()),
    ]

    for target_dir in [OUT_DIR, LOCAL_SELECTION_DIR]:
        review.to_csv(target_dir / "nonlipid_module_compact_review.tsv", sep="\t", index=False)
        kept.to_csv(target_dir / "nonlipid_module_compact_kept.tsv", sep="\t", index=False)
        dropped.to_csv(target_dir / "nonlipid_module_compact_dropped.tsv", sep="\t", index=False)
        group_counts.to_csv(target_dir / "nonlipid_module_compact_group_counts.tsv", sep="\t", index=False)
        (target_dir / "selection_summary.txt").write_text("\n".join(summary_lines), encoding="utf-8")


if __name__ == "__main__":
    main()
