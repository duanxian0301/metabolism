from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
MASTER = ROOT / "results" / "17_candidate_gene_integration_lipid8_F2_AD" / "lipid8_F2_AD_cross_evidence_master.tsv"
OUTDIR = ROOT / "results" / "18_virtual_knockout_lipid8_F2_AD"


def bool_col(df: pd.DataFrame, col: str) -> pd.Series:
    return df[col].fillna(False).astype(bool)


def evidence_summary(row: pd.Series) -> str:
    bits = []
    if bool(row.get("fuma_posMap_supported")):
        bits.append("FUMA-posMap")
    if bool(row.get("fuma_eqtlMap_supported")):
        bits.append("FUMA-eQTL")
    if bool(row.get("fuma_ciMap_supported")):
        bits.append("FUMA-ciMap")
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

    fuma_any = (
        bool_col(df, "fuma_posMap_supported")
        | bool_col(df, "fuma_eqtlMap_supported")
        | bool_col(df, "fuma_ciMap_supported")
        | bool_col(df, "fuma_supported")
    )
    snp_region_support = bool_col(df, "coloc_supported") | bool_col(df, "pwcoco_supported") | bool_col(df, "fuma_supported")
    gene_layer_support = (
        bool_col(df, "smr_bulk_supported")
        | bool_col(df, "smr_gtex_supported")
        | bool_col(df, "smr_bryois_supported")
        | bool_col(df, "ctwas_supported")
    )

    strict = df[fuma_any & snp_region_support & gene_layer_support].copy()
    strict["evidence_summary"] = strict.apply(evidence_summary, axis=1)

    strict_keep = [
        "Gene",
        "priority_tier",
        "priority_score",
        "n_evidence_layers",
        "fuma_posMap_supported",
        "fuma_eqtlMap_supported",
        "fuma_ciMap_supported",
        "coloc_supported",
        "pwcoco_supported",
        "smr_bulk_supported",
        "smr_gtex_supported",
        "smr_bryois_supported",
        "smr_bryois_best_celltype",
        "ctwas_supported",
        "max_pip",
        "evidence_summary",
    ]
    strict = strict[strict_keep].sort_values(
        ["priority_score", "n_evidence_layers", "max_pip"], ascending=[False, False, False]
    )
    strict.to_csv(OUTDIR / "lipid8_F2_AD_virtual_ko_strict_candidate_set.tsv", sep="\t", index=False)

    strict_cell = strict[strict["smr_bryois_best_celltype"].isin(["Pericytes", "OPCs...COPs"])].copy()
    strict_cell.to_csv(OUTDIR / "lipid8_F2_AD_virtual_ko_strict_pe_opc_only.tsv", sep="\t", index=False)

    # If PE/OPC is empty or weak, relax only the cell-type consistency requirement but retain SNP/FUMA support.
    relaxed_genes = ["CAB39L", "LRRC37A2", "LRRC37A", "ARL17B", "KNOP1", "MAP1LC3A", "YPEL3", "INO80E"]
    relaxed = strict[strict["Gene"].isin(relaxed_genes)].copy()
    relaxed["recommended_cell_branch"] = relaxed["smr_bryois_best_celltype"].replace(
        {
            "Astrocytes": "astrocyte",
            "Excitatory.neurons": "excitatory_neuron",
            "Oligodendrocytes": "oligodendrocyte",
        }
    )
    relaxed["recommended_reason"] = [
        "strict SNP/FUMA-supported candidate; chosen without forcing scPagwas cell consistency"
        for _ in range(len(relaxed))
    ]
    relaxed.to_csv(OUTDIR / "lipid8_F2_AD_virtual_ko_relaxed_no_cell_consistency.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"rule_set": "strict_snp_fuma_gene", "n_total": int(len(strict))},
            {"rule_set": "strict_snp_fuma_gene_pe_opc_only", "n_total": int(len(strict_cell))},
            {"rule_set": "relaxed_no_cell_consistency_recommended", "n_total": int(len(relaxed))},
        ]
    )
    summary.to_csv(OUTDIR / "lipid8_F2_AD_virtual_ko_rule_set_counts.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
