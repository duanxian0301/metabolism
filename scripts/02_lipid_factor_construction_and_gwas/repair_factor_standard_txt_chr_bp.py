from pathlib import Path

import pandas as pd


BIM = Path(r"D:\SMR\g1000\g1000_eur.bim")
TARGETS = [
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F1_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F2_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\standard_txt\lipid_final8_F3_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F1_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F2_standard.txt"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\standard_txt\nonlipid_final8_F3_standard.txt"),
]
SUMMARY = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\factor_standard_txt_chr_bp_repair_summary.tsv")


def main():
    coord = pd.read_csv(
        BIM,
        sep="\t",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "BIM_A1", "BIM_A2"],
        usecols=["CHR", "SNP", "BP"],
    ).drop_duplicates("SNP")

    summary_rows = []
    for path in TARGETS:
        dt = pd.read_csv(path, sep="\t")
        had_chr_bp = {"CHR", "BP"}.issubset(dt.columns)
        if had_chr_bp:
            dt = dt.drop(columns=["CHR", "BP"])
        dt = dt.merge(coord, on="SNP", how="left")
        ordered = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
        dt = dt[ordered]
        matched = dt["CHR"].notna().sum()
        dt.to_csv(path, sep="\t", index=False)
        summary_rows.append(
            {
                "file": str(path),
                "n_rows": len(dt),
                "n_chr_bp_filled": int(matched),
                "pct_chr_bp_filled": matched / len(dt),
                "n_chr_bp_missing": int(len(dt) - matched),
                "coord_source": str(BIM),
            }
        )

    pd.DataFrame(summary_rows).to_csv(SUMMARY, sep="\t", index=False)
    print(SUMMARY)


if __name__ == "__main__":
    main()
