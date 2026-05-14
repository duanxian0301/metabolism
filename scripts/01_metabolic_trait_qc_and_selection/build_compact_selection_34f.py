import pandas as pd
from pathlib import Path


EXCEL_PATH = Path(r"D:\metabolic\metabolite_FGWAS_selection_lists.xlsx")
SHEET_NAME = "Main_Zgt4_nonproportion"
OUT_DIR = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\compact_selection_34f")

SELECTED_TRAITS = [
    # F1-like: VLDL / CE / TG transport
    "XS_VLDL_FC", "XS_VLDL_L", "VLDL_CE", "S_VLDL_FC", "S_VLDL_CE",
    "XS_VLDL_CE", "S_VLDL_L", "XS_VLDL_TG", "M_LDL_TG", "S_LDL_TG", "Total_L",
    # F2-like: HDL / phospholipid / membrane-lipid composition
    "ApoA1", "M_HDL_CE", "S_HDL_C", "S_HDL_CE", "S_HDL_L", "HDL_CE",
    "Phosphatidylc", "Cholines", "Phosphoglyc", "Sphingomyelins", "DHA", "Omega_3",
    # F3-like: BCAA / particle size / TG-related metabolism
    "Ile", "Leu", "Val", "Phe", "LDL_size", "VLDL_size", "HDL_size", "M_HDL_TG", "S_HDL_TG",
    # F4-like: ketone-body axis
    "bOHbutyrate", "Acetoacetate", "Acetate",
]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_NAME)
    df = df.rename(columns={"trait": "trait_code", "biomarker name": "biomarker_name"})
    df["trait_code"] = df["trait_code"].astype(str)

    missing = sorted(set(SELECTED_TRAITS) - set(df["trait_code"]))
    if missing:
        raise ValueError(f"Traits missing from sheet {SHEET_NAME}: {missing}")

    out = df[df["trait_code"].isin(SELECTED_TRAITS)].copy()
    out["selection_scheme"] = "balanced_compact_34factor_v1"
    out["selection_rationale"] = "retain strongest and cleanest indicators for 3-4-factor identification and better CFA fit"
    out["selection_order"] = out["trait_code"].apply(SELECTED_TRAITS.index)
    out = out.sort_values("selection_order")

    out.to_csv(OUT_DIR / "traits_compact_annotated.tsv", sep="\t", index=False)
    out[["trait_code"]].to_csv(OUT_DIR / "traits_compact.tsv", sep="\t", index=False)
    (
        out.groupby("group")["trait_code"]
        .count()
        .reset_index(name="n_kept")
        .sort_values(["n_kept", "group"], ascending=[False, True])
        .to_csv(OUT_DIR / "traits_compact_group_counts.tsv", sep="\t", index=False)
    )

    print(f"compact_n\t{len(out)}")
    print(out[["trait_code", "biomarker_name", "group"]].to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
