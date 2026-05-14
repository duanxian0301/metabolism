from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
MASTER = ROOT / "results" / "17_candidate_gene_integration_lipid8_F2_AD" / "lipid8_F2_AD_cross_evidence_master.tsv"
SC_CELL = ROOT / "results" / "12_scpagwas2_lipid8_F2_MSSM_AD" / "lipid8_F2_MSSM_AD_Merged_celltype_pvalue_withFDR.csv"
OUTDIR = ROOT / "results" / "18_virtual_knockout_lipid8_F2_AD"


def evidence_text(row: pd.Series) -> str:
    bits = []
    if bool(row.get("fuma_supported")):
        bits.append("FUMA")
    if bool(row.get("coloc_supported")):
        bits.append("coloc")
    if bool(row.get("pwcoco_supported")):
        bits.append("PWCoCo")
    if bool(row.get("smr_bulk_supported")):
        bits.append("bulk SMR")
    if bool(row.get("smr_gtex_supported")):
        bits.append("GTEx SMR")
    if bool(row.get("smr_bryois_supported")):
        bits.append("Bryois SMR")
    if bool(row.get("ctwas_supported")):
        bits.append("cTWAS")
    return "; ".join(bits)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(MASTER, sep="\t", low_memory=False)
    sc = pd.read_csv(SC_CELL)

    sc = sc.sort_values(["celltype_FDR", "pvalue"])
    top_sc = sc[["celltype", "pvalue", "celltype_FDR"]].copy()
    top_sc["scpagwas_priority"] = [
        "primary" if fdr < 0.05 else "secondary" if fdr < 0.10 else "background"
        for fdr in top_sc["celltype_FDR"]
    ]

    design_rows = [
        {
            "analysis_branch": "core",
            "cell_branch": "astrocyte",
            "gene": "CAB39L",
            "group_id": "astrocyte__CAB39L",
            "priority_round": 1,
            "selection_basis": "A-tier convergent gene with FUMA+coloc+PWCoCo+GTEx SMR+Bryois astrocyte SMR+cTWAS",
        },
        {
            "analysis_branch": "core",
            "cell_branch": "astrocyte",
            "gene": "LRRC37A2",
            "group_id": "astrocyte__LRRC37A2",
            "priority_round": 1,
            "selection_basis": "A-tier shared-convergent locus gene with Bryois astrocyte SMR+GTEx SMR+cTWAS+PWCoCo",
        },
        {
            "analysis_branch": "core",
            "cell_branch": "astrocyte",
            "gene": "KNOP1",
            "group_id": "astrocyte__KNOP1",
            "priority_round": 1,
            "selection_basis": "B-tier AD-weighted gene with bulk SMR+GTEx SMR+Bryois astrocyte SMR+cTWAS",
        },
        {
            "analysis_branch": "core",
            "cell_branch": "excitatory_neuron",
            "gene": "LRRC37A",
            "group_id": "excitatory_neuron__LRRC37A",
            "priority_round": 1,
            "selection_basis": "A-tier shared-convergent gene with bulk SMR+GTEx SMR+Bryois excitatory neuron SMR+PWCoCo+cTWAS",
        },
        {
            "analysis_branch": "core",
            "cell_branch": "excitatory_neuron",
            "gene": "MAP1LC3A",
            "group_id": "excitatory_neuron__MAP1LC3A",
            "priority_round": 1,
            "selection_basis": "B-tier gene with GTEx SMR+Bryois excitatory neuron SMR+cTWAS and neuronal autophagy relevance",
        },
        {
            "analysis_branch": "core",
            "cell_branch": "oligodendrocyte",
            "gene": "ARL17B",
            "group_id": "oligodendrocyte__ARL17B",
            "priority_round": 1,
            "selection_basis": "A-tier 17q21 convergent gene with Bryois oligodendrocyte SMR+GTEx SMR+PWCoCo+cTWAS",
        },
        {
            "analysis_branch": "core",
            "cell_branch": "oligodendrocyte",
            "gene": "PEX6",
            "group_id": "oligodendrocyte__PEX6",
            "priority_round": 1,
            "selection_basis": "B-tier gene with bulk SMR+GTEx SMR+Bryois oligodendrocyte SMR+cTWAS and strong PIP",
        },
        {
            "analysis_branch": "core",
            "cell_branch": "opc_cop",
            "gene": "ANKRD36B",
            "group_id": "opc_cop__ANKRD36B",
            "priority_round": 1,
            "selection_basis": "B-tier gene with bulk SMR+GTEx SMR+Bryois OPC/COP SMR+cTWAS",
        },
        {
            "analysis_branch": "exploratory",
            "cell_branch": "pericyte",
            "gene": "TSPAN14",
            "group_id": "pericyte__TSPAN14",
            "priority_round": 2,
            "selection_basis": "Exploratory pericyte branch because scPagwas2 pericyte is the only FDR-significant cell type; gene chosen for AD-relevant bulk+GTEx SMR support",
        },
        {
            "analysis_branch": "exploratory",
            "cell_branch": "pericyte",
            "gene": "ACE",
            "group_id": "pericyte__ACE",
            "priority_round": 2,
            "selection_basis": "Exploratory pericyte branch because scPagwas2 pericyte is the only FDR-significant cell type; gene chosen for vascular/AD interpretability and GTEx support",
        },
    ]

    design = pd.DataFrame(design_rows)
    merged = design.merge(df, how="left", left_on="gene", right_on="Gene")

    merged["evidence_summary"] = merged.apply(evidence_text, axis=1)
    merged["scpagwas_direct_support"] = merged["cell_branch"].map(
        {
            "pericyte": "Yes: scPagwas2 FDR<0.05",
            "astrocyte": "No direct scPagwas2 signal in current run; selected from cross-method gene evidence",
            "excitatory_neuron": "No direct scPagwas2 signal in current run; selected from cross-method gene evidence",
            "oligodendrocyte": "Secondary scPagwas2 support nearby in OPC branch (FDR<0.10)",
            "opc_cop": "Secondary scPagwas2 support nearby in OPC branch (FDR<0.10)",
        }
    )
    merged["recommended_use"] = merged["priority_round"].map(
        {1: "first-round KO", 2: "second-round exploratory KO"}
    )

    keep_cols = [
        "analysis_branch",
        "priority_round",
        "recommended_use",
        "cell_branch",
        "gene",
        "group_id",
        "priority_tier",
        "priority_score",
        "n_evidence_layers",
        "evidence_summary",
        "selection_basis",
        "scpagwas_direct_support",
        "smr_bryois_best_celltype",
        "smr_bulk_supported",
        "smr_gtex_supported",
        "smr_bryois_supported",
        "ctwas_supported",
        "coloc_supported",
        "pwcoco_supported",
        "max_pip",
    ]
    merged = merged[keep_cols].sort_values(["priority_round", "cell_branch", "priority_score"], ascending=[True, True, False])

    merged.to_csv(OUTDIR / "lipid8_F2_AD_virtual_ko_design.tsv", sep="\t", index=False)
    top_sc.to_csv(OUTDIR / "lipid8_F2_AD_scpagwas2_celltype_priority.tsv", sep="\t", index=False)

    counts = pd.DataFrame(
        [
            {"set_name": "core_first_round_groups", "n_groups": int((merged["priority_round"] == 1).sum())},
            {"set_name": "exploratory_second_round_groups", "n_groups": int((merged["priority_round"] == 2).sum())},
            {"set_name": "all_planned_groups", "n_groups": int(len(merged))},
        ]
    )
    counts.to_csv(OUTDIR / "lipid8_F2_AD_virtual_ko_group_counts.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
