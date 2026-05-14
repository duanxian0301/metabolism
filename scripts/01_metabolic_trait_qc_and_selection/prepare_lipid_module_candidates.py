from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
STEP1 = ROOT / "step1_ldsc_results"
OUT_DIR = ROOT / "lipid_module_from_full_manifest"

LIPID_GROUPS = {
    "Apolipoproteins",
    "Cholesterol",
    "Cholesteryl esters",
    "Fatty acids",
    "Free cholesterol",
    "Lipoprotein particle sizes",
    "Lipoprotein subclasses",
    "Other lipids",
    "Phospholipids",
    "Total lipids",
    "Triglycerides",
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(STEP1 / "trait_manifest.tsv", sep="\t")
    summary = pd.read_csv(STEP1 / "Main_Zgt4_nonproportion_ldsc_summary.tsv", sep="\t")

    df = manifest.merge(summary, left_on="trait_code", right_on="trait", how="left")
    df["is_lipid_module"] = df["group"].isin(LIPID_GROUPS)
    lipid = df[df["is_lipid_module"]].copy()

    lipid = lipid[
        [
            "study_accession",
            "trait_code",
            "biomarker_name",
            "group",
            "h2",
            "intercept",
            "sumstats_file",
        ]
    ].sort_values(["group", "trait_code"])

    lipid.to_csv(OUT_DIR / "lipid_module_candidates.tsv", sep="\t", index=False)

    counts = (
        lipid.groupby("group")
        .size()
        .reset_index(name="n_traits")
        .sort_values(["n_traits", "group"], ascending=[False, True])
    )
    counts.to_csv(OUT_DIR / "lipid_module_group_counts.tsv", sep="\t", index=False)

    summary_text = "\n".join(
        [
            "Lipid module candidate panel",
            f"Total traits in full Main_Zgt4_nonproportion manifest: {len(df)}",
            f"Lipid-module candidate traits: {len(lipid)}",
            "Included groups: "
            + ", ".join(sorted(LIPID_GROUPS)),
            "Rationale: restart module-specific modeling from multivariate LDSC using a lipid-only trait universe.",
        ]
    )
    (OUT_DIR / "summary.txt").write_text(summary_text, encoding="utf-8")


if __name__ == "__main__":
    main()
