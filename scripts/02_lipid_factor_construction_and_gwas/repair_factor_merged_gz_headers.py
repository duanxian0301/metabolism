from pathlib import Path

import pandas as pd


TASKS = [
    (
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\lipid_final8_F1_userGWAS_merged.tsv.gz"),
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F1_standard.txt"),
    ),
    (
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\lipid_final8_F2_userGWAS_merged.tsv.gz"),
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F2_standard.txt"),
    ),
    (
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\lipid_final8_F3_userGWAS_merged.tsv.gz"),
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F3_standard.txt"),
    ),
    (
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\nonlipid_final8_F1_userGWAS_merged.tsv.gz"),
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F1_standard.txt"),
    ),
    (
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\nonlipid_final8_F2_userGWAS_merged.tsv.gz"),
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F2_standard.txt"),
    ),
    (
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\nonlipid_final8_F3_userGWAS_merged.tsv.gz"),
        Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F3_standard.txt"),
    ),
]

SUMMARY = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\factor_merged_gz_header_repair_summary.tsv")


def repair_one(merged_path: Path, standard_path: Path):
    merged = pd.read_csv(merged_path, sep="\t", compression="gzip")
    lookup = pd.read_csv(standard_path, sep="\t", usecols=["SNP", "CHR", "BP", "A1", "A2", "FRQ"])
    lookup = lookup.rename(columns={"FRQ": "MAF"})

    wrong_prefix = [c for c in ["CHR", "BP", "MAF", "A1", "A2"] if c in merged.columns]
    suffix_cols = [c for c in merged.columns if c not in ["SNP"] + wrong_prefix]

    repaired = merged[["SNP"] + suffix_cols].merge(lookup, on="SNP", how="left")
    ordered = ["SNP", "CHR", "BP", "MAF", "A1", "A2"] + suffix_cols
    repaired = repaired[ordered]
    repaired.to_csv(merged_path, sep="\t", index=False, compression="gzip")

    matched = repaired["CHR"].notna().sum()
    return {
        "merged_file": str(merged_path),
        "standard_lookup": str(standard_path),
        "n_rows": len(repaired),
        "n_chr_bp_filled": int(matched),
        "pct_chr_bp_filled": matched / len(repaired),
        "n_chr_bp_missing": int(len(repaired) - matched),
    }


def main():
    rows = [repair_one(merged_path, standard_path) for merged_path, standard_path in TASKS]
    pd.DataFrame(rows).to_csv(SUMMARY, sep="\t", index=False)
    print(SUMMARY)


if __name__ == "__main__":
    main()
