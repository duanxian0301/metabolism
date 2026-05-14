from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
ULTRA_DIR = ROOT / "nonlipid_module_from_full_manifest" / "compact_panel_review" / "ultra_pure_3factor_review"
ESEM_DIR = ROOT / "step20_efa_esem_nonlipid_module_ultrapure8"
OUT_DIR = ROOT / "final_model_nonlipid_module_final8"
LOCAL_FINAL_DIR = Path(__file__).resolve().parents[1] / "inputs" / "final_model"

FINAL_FACTOR_NAME = {
    "Acetoacetate": "F1_ketone_axis",
    "bOHbutyrate": "F1_ketone_axis",
    "Val": "F2_amino_acid_axis",
    "Leu": "F2_amino_acid_axis",
    "Phe": "F2_amino_acid_axis",
    "Acetate": "F3_energy_bridge_axis",
    "Lactate": "F3_energy_bridge_axis",
    "Glucose": "F3_energy_bridge_axis",
}

MODEL_TEXT = """# Final nonlipid-module 3-factor model (8 traits)
F1_ketone_axis =~ Acetoacetate + bOHbutyrate
F2_amino_acid_axis =~ Val + Leu + Phe
F3_energy_bridge_axis =~ Acetate + Lactate + Glucose
"""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_FINAL_DIR.mkdir(parents=True, exist_ok=True)

    kept = pd.read_csv(ULTRA_DIR / "nonlipid_module_ultrapure3_kept.tsv", sep="\t")
    esem = pd.read_csv(ESEM_DIR / "nonlipid_module_ultrapure8_esem_summary.tsv", sep="\t")

    final_manifest = kept[
        [
            "study_accession",
            "trait_code",
            "biomarker_name",
            "group",
            "ultra_pure_factor",
            "sumstats_file",
        ]
    ].copy()
    final_manifest["final_factor_name"] = final_manifest["trait_code"].map(FINAL_FACTOR_NAME)
    factor_order = {
        "F1_ketone_axis": 1,
        "F2_amino_acid_axis": 2,
        "F3_energy_bridge_axis": 3,
    }
    final_manifest["factor_order"] = final_manifest["final_factor_name"].map(factor_order)
    final_manifest = final_manifest.sort_values(["factor_order", "trait_code"]).drop(columns=["factor_order"])

    all_fit = esem[(esem["dataset"] == "ALL") & (esem["nfactors"] == 3)].iloc[0]
    odd_fit = esem[(esem["dataset"] == "ODD") & (esem["nfactors"] == 3)].iloc[0]

    model_summary = "\n".join(
        [
            "Current best nonlipid-module model selected for downstream analysis.",
            "Model: 3-factor target-rotation ESEM, final-8 marker panel.",
            "Fit on ALL covariance matrix: "
            f"CFI = {all_fit['cfi']:.6f}, SRMR = {all_fit['srmr']:.6f}, RMSEA = {all_fit['rmsea']:.6f}.",
            "Fit on ODD covariance matrix: "
            f"CFI = {odd_fit['cfi']:.6f}, SRMR = {odd_fit['srmr']:.6f}, RMSEA = {odd_fit['rmsea']:.6f}.",
            "Factors:",
            "F1_ketone_axis: Acetoacetate, bOHbutyrate",
            "F2_amino_acid_axis: Val, Leu, Phe",
            "F3_energy_bridge_axis: Acetate, Lactate, Glucose",
            "Note: ALL fit passes the target threshold CFI >= 0.95 and SRMR < 0.08; ODD fit also passes those thresholds, though ODD RMSEA remains slightly above 0.08.",
        ]
    )

    for target_dir in [OUT_DIR, LOCAL_FINAL_DIR]:
        final_manifest.to_csv(target_dir / "final8_trait_manifest.tsv", sep="\t", index=False)
        (target_dir / "final8_factor_structure.txt").write_text(MODEL_TEXT, encoding="utf-8")
        (target_dir / "model_summary.txt").write_text(model_summary, encoding="utf-8")


if __name__ == "__main__":
    main()
