import pandas as pd
from pathlib import Path


EXCEL_PATH = Path(r"D:\metabolic\metabolite_FGWAS_selection_lists.xlsx")
SHEET_NAME = "Main_Zgt4_nonproportion"
OUT_DIR = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\balanced_selection_3to4f")

SELECTED_TRAITS = [
    # Amino acids
    "Ala", "Gln", "His", "Ile", "Leu", "Phe", "Tyr", "Val",
    # Fatty acids
    "DHA", "LA", "MUFA", "Omega_3", "Omega_6", "PUFA", "SFA", "Total_FA",
    # Fluid balance / glycolysis / inflammation / ketone bodies
    "Albumin", "Creatinine", "Citrate", "Glucose", "Lactate", "GlycA",
    "Acetate", "Acetoacetate", "Acetone", "bOHbutyrate",
    # Broad lipid composition
    "ApoA1", "HDL_CE", "Total_CE", "VLDL_CE",
    "Cholines", "Phosphatidylc", "Phosphoglyc", "Sphingomyelins", "Total_L",
    # Particle sizes
    "HDL_size", "LDL_size", "VLDL_size",
    # HDL / LDL / VLDL subclass anchors
    "L_HDL_C", "M_HDL_CE", "M_HDL_TG",
    "S_HDL_C", "S_HDL_CE", "S_HDL_L", "S_HDL_TG",
    "M_LDL_TG", "S_LDL_TG",
    "S_VLDL_CE", "S_VLDL_FC", "S_VLDL_L", "S_VLDL_TG",
    "XL_HDL_TG",
    "XS_VLDL_CE", "XS_VLDL_FC", "XS_VLDL_L", "XS_VLDL_TG",
]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_NAME)
    df = df.rename(
        columns={
            "trait": "trait_code",
            "biomarker name": "biomarker_name",
        }
    )
    df["trait_code"] = df["trait_code"].astype(str)

    missing = sorted(set(SELECTED_TRAITS) - set(df["trait_code"]))
    if missing:
        raise ValueError(f"Traits missing from sheet {SHEET_NAME}: {missing}")

    out = df[df["trait_code"].isin(SELECTED_TRAITS)].copy()
    out["selection_scheme"] = "balanced_manual_3to4factor_v1"
    out["selection_rationale"] = "retain biological breadth across metabolic domains and lipoprotein subdimensions"
    out["selection_order"] = out["trait_code"].apply(SELECTED_TRAITS.index)
    out = out.sort_values("selection_order")

    out.to_csv(OUT_DIR / "traits_balanced_annotated.tsv", sep="\t", index=False)
    out[["trait_code"]].to_csv(OUT_DIR / "traits_balanced.tsv", sep="\t", index=False)
    (
        out.groupby("group")["trait_code"]
        .count()
        .reset_index(name="n_kept")
        .sort_values(["n_kept", "group"], ascending=[False, True])
        .to_csv(OUT_DIR / "traits_balanced_group_counts.tsv", sep="\t", index=False)
    )

    print(f"balanced_n\t{len(out)}")
    print(out[["trait_code", "biomarker_name", "group"]].to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
