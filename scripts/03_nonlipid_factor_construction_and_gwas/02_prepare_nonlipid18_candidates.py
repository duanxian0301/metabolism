from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
STEP1 = ROOT / "step1_ldsc_results"
OUT_DIR = ROOT / "nonlipid_module_from_full_manifest"
LOCAL_SELECTION_DIR = Path(__file__).resolve().parents[1] / "inputs" / "selection"

ALL_GROUP_METADATA = {
    "Amino acids": {
        "include_in_nonlipid": True,
        "proposed_submodule": "amino_acid_axis",
        "biological_logic": "Canonical non-lipid small molecules, including BCAA and aromatic amino acids.",
    },
    "Glycolysis related metabolites": {
        "include_in_nonlipid": True,
        "proposed_submodule": "energy_metabolism_axis",
        "biological_logic": "Central carbon and glycolysis-related intermediates that should be evaluated jointly with ketone markers.",
    },
    "Ketone bodies": {
        "include_in_nonlipid": True,
        "proposed_submodule": "energy_metabolism_axis",
        "biological_logic": "Non-lipid ketone body markers expected to inform an energy/fasting axis.",
    },
    "Inflammation": {
        "include_in_nonlipid": True,
        "proposed_submodule": "renal_inflammation_bridge",
        "biological_logic": "Inflammatory bridge indicator retained for joint non-lipid model search.",
    },
    "Fluid balance": {
        "include_in_nonlipid": True,
        "proposed_submodule": "renal_inflammation_bridge",
        "biological_logic": "Albumin and creatinine are non-lipid sentinels that may bridge hydration, renal handling, and systemic metabolism.",
    },
    "Apolipoproteins": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because apolipoproteins are core lipid/lipoprotein measures.",
    },
    "Cholesterol": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because cholesterol traits are part of the lipid/lipoprotein module.",
    },
    "Cholesteryl esters": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because cholesteryl ester traits are lipid-core measures.",
    },
    "Fatty acids": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because fatty acid traits belong to the lipid module definition for this project.",
    },
    "Free cholesterol": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because free cholesterol traits are lipid-core measures.",
    },
    "Lipoprotein particle sizes": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because lipoprotein particle-size measures are part of the lipid module.",
    },
    "Lipoprotein subclasses": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because subclass-resolved lipoprotein measures define the lipid module.",
    },
    "Other lipids": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because the retained traits in this group are lipid-structure measurements rather than non-lipid metabolites.",
    },
    "Phospholipids": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because phospholipid traits are lipid-core measures.",
    },
    "Total lipids": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because total-lipid traits belong to the lipid module.",
    },
    "Triglycerides": {
        "include_in_nonlipid": False,
        "proposed_submodule": "excluded_lipid",
        "biological_logic": "Excluded because triglyceride traits belong to the lipid module.",
    },
}

TRAIT_METADATA = {
    "Ala": {
        "marker_role": "support_marker",
        "candidate_bucket": "general_amino_acid",
        "inclusion_rationale": "Retain as a broad amino-acid marker to test whether an amino-acid factor extends beyond BCAA/aromatic signals.",
    },
    "Gln": {
        "marker_role": "support_marker",
        "candidate_bucket": "general_amino_acid",
        "inclusion_rationale": "Retain as a high-h2 non-BCAA amino-acid marker that may anchor a distinct glutamine/glutamate-adjacent signal.",
    },
    "His": {
        "marker_role": "support_marker",
        "candidate_bucket": "general_amino_acid",
        "inclusion_rationale": "Retain as a lower-correlation amino-acid indicator that may improve breadth if it loads cleanly in EFA/ESEM.",
    },
    "Ile": {
        "marker_role": "core_marker",
        "candidate_bucket": "bcaa_axis",
        "inclusion_rationale": "Retain as a core BCAA marker with very high heritability Z and strong shared structure with leucine and valine.",
    },
    "Leu": {
        "marker_role": "core_marker",
        "candidate_bucket": "bcaa_axis",
        "inclusion_rationale": "Retain as a core BCAA marker; redundancy will be reviewed during compact reduction rather than excluded upfront.",
    },
    "Phe": {
        "marker_role": "core_marker",
        "candidate_bucket": "aromatic_axis",
        "inclusion_rationale": "Retain as an aromatic amino-acid marker with stronger pairing to tyrosine than to the ketone or fluid traits.",
    },
    "Tyr": {
        "marker_role": "core_marker",
        "candidate_bucket": "aromatic_axis",
        "inclusion_rationale": "Retain as an aromatic amino-acid marker and potential bridge between aromatic and broader amino-acid variation.",
    },
    "Val": {
        "marker_role": "core_marker",
        "candidate_bucket": "bcaa_axis",
        "inclusion_rationale": "Retain as a core BCAA marker with the strongest observed rg links in the amino-acid cluster.",
    },
    "Albumin": {
        "marker_role": "sentinel_marker",
        "candidate_bucket": "renal_fluid_axis",
        "inclusion_rationale": "Retain as a non-lipid serum abundance / fluid-balance sentinel for testing whether it contributes a distinct systemic factor.",
    },
    "Creatinine": {
        "marker_role": "sentinel_marker",
        "candidate_bucket": "renal_fluid_axis",
        "inclusion_rationale": "Retain as a renal-handling sentinel with strong h2 and clear non-lipid interpretation.",
    },
    "Citrate": {
        "marker_role": "support_marker",
        "candidate_bucket": "central_carbon_axis",
        "inclusion_rationale": "Retain as a low-redundancy TCA/central-carbon marker that may survive compact reduction if it contributes unique structure.",
    },
    "Glucose": {
        "marker_role": "core_marker",
        "candidate_bucket": "glycolysis_glucose_axis",
        "inclusion_rationale": "Retain as the canonical glycemic marker for the energy-metabolism block.",
    },
    "Lactate": {
        "marker_role": "core_marker",
        "candidate_bucket": "glycolysis_glucose_axis",
        "inclusion_rationale": "Retain as a glycolytic/anaerobic metabolism marker and candidate bridge to inflammation or ketone biology.",
    },
    "GlycA": {
        "marker_role": "bridge_marker",
        "candidate_bucket": "inflammation_axis",
        "inclusion_rationale": "Retain as the only inflammation trait because it may act as a bridge across amino-acid and glucose-related variation.",
    },
    "Acetate": {
        "marker_role": "core_marker",
        "candidate_bucket": "ketone_fasting_axis",
        "inclusion_rationale": "Retain as a ketone/fasting-related marker with distinctive negative correlation to lactate and GlycA.",
    },
    "Acetoacetate": {
        "marker_role": "core_marker",
        "candidate_bucket": "ketone_fasting_axis",
        "inclusion_rationale": "Retain as a canonical ketone-body indicator and likely core anchor for any ketone factor.",
    },
    "Acetone": {
        "marker_role": "support_marker",
        "candidate_bucket": "ketone_fasting_axis",
        "inclusion_rationale": "Retain as a ketone-body support marker; redundancy with acetoacetate will be assessed at compact reduction.",
    },
    "bOHbutyrate": {
        "marker_role": "core_marker",
        "candidate_bucket": "ketone_fasting_axis",
        "inclusion_rationale": "Retain as a canonical beta-hydroxybutyrate marker and likely ketone-axis anchor.",
    },
}


def build_group_table(manifest: pd.DataFrame) -> pd.DataFrame:
    counts = manifest.groupby("group").size().rename("n_traits").reset_index()
    group_df = pd.DataFrame(
        [
            {
                "group": group,
                "include_in_nonlipid": meta["include_in_nonlipid"],
                "proposed_submodule": meta["proposed_submodule"],
                "biological_logic": meta["biological_logic"],
            }
            for group, meta in ALL_GROUP_METADATA.items()
        ]
    )
    group_df = counts.merge(group_df, on="group", how="left").sort_values(
        ["include_in_nonlipid", "n_traits", "group"], ascending=[False, False, True]
    )
    return group_df


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_SELECTION_DIR.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(STEP1 / "trait_manifest.tsv", sep="\t")
    summary = pd.read_csv(STEP1 / "Main_Zgt4_nonproportion_ldsc_summary.tsv", sep="\t")
    qc = pd.read_csv(STEP1 / "efa_selection" / "trait_qc_and_selection.tsv", sep="\t")

    all_df = (
        manifest.merge(summary, left_on="trait_code", right_on="trait", how="left")
        .merge(
            qc[
                [
                    "trait_code",
                    "Z",
                    "selected_for_efa",
                    "max_abs_rg_with_any_trait",
                    "n_independent_leads",
                ]
            ],
            on="trait_code",
            how="left",
        )
        .drop(columns=["trait"])
    )

    group_table = build_group_table(manifest)
    included_groups = set(group_table.loc[group_table["include_in_nonlipid"], "group"])

    candidate = all_df[all_df["group"].isin(included_groups)].copy()
    candidate["module_name"] = "nonlipid_main"
    candidate["proposed_submodule"] = candidate["group"].map(
        lambda g: ALL_GROUP_METADATA[g]["proposed_submodule"]
    )
    candidate["marker_role"] = candidate["trait_code"].map(
        lambda t: TRAIT_METADATA[t]["marker_role"]
    )
    candidate["candidate_bucket"] = candidate["trait_code"].map(
        lambda t: TRAIT_METADATA[t]["candidate_bucket"]
    )
    candidate["inclusion_rationale"] = candidate["trait_code"].map(
        lambda t: TRAIT_METADATA[t]["inclusion_rationale"]
    )
    candidate["priority_for_compact_review"] = candidate["marker_role"].map(
        {
            "core_marker": "high",
            "bridge_marker": "high",
            "sentinel_marker": "medium",
            "support_marker": "medium",
        }
    )

    candidate = candidate[
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
            "inclusion_rationale",
            "sumstats_file",
        ]
    ].sort_values(["proposed_submodule", "group", "trait_code"])

    candidate_counts = (
        candidate.groupby(["module_name", "proposed_submodule", "group"])
        .size()
        .reset_index(name="n_traits")
        .sort_values(["module_name", "n_traits", "group"], ascending=[True, False, True])
    )

    module_grouping = (
        candidate[
            [
                "group",
                "proposed_submodule",
                "candidate_bucket",
                "marker_role",
                "trait_code",
                "biomarker_name",
            ]
        ]
        .sort_values(["proposed_submodule", "candidate_bucket", "trait_code"])
        .reset_index(drop=True)
    )

    summary_text = "\n".join(
        [
            "Nonlipid module candidate panel",
            f"Total traits in full Main_Zgt4_nonproportion manifest: {len(all_df)}",
            f"Nonlipid candidate traits: {len(candidate)}",
            "Included groups: "
            + ", ".join(sorted(included_groups)),
            "Excluded groups: "
            + ", ".join(sorted(set(manifest['group']) - included_groups)),
            "Recommended top-level strategy: keep one nonlipid main module, then use staged reduction within the module.",
        ]
    )

    for target_dir in [OUT_DIR, LOCAL_SELECTION_DIR]:
        candidate.to_csv(target_dir / "nonlipid_module_candidates.tsv", sep="\t", index=False)
        candidate_counts.to_csv(target_dir / "nonlipid_module_group_counts.tsv", sep="\t", index=False)
        group_table.to_csv(target_dir / "nonlipid_group_inclusion_table.tsv", sep="\t", index=False)
        module_grouping.to_csv(target_dir / "nonlipid_module_grouping.tsv", sep="\t", index=False)
        (target_dir / "summary.txt").write_text(summary_text, encoding="utf-8")


if __name__ == "__main__":
    main()
