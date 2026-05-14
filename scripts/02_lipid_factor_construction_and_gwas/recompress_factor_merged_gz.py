from pathlib import Path
import gzip
import shutil


FILES = [
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\lipid_final8_F1_userGWAS_merged.tsv.gz"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\lipid_final8_F2_userGWAS_merged.tsv.gz"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step14_native_wsl_usergwas_final8_results\merged_lipid_final8\lipid_final8_F3_userGWAS_merged.tsv.gz"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\nonlipid_final8_F1_userGWAS_merged.tsv.gz"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\nonlipid_final8_F2_userGWAS_merged.tsv.gz"),
    Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step22_native_wsl_usergwas_nonlipid8_results\merged_nonlipid_final8\nonlipid_final8_F3_userGWAS_merged.tsv.gz"),
]


def recompress(path: Path):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(path, "rb") as src, gzip.open(tmp, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    path.unlink()
    tmp.rename(path)


def main():
    for path in FILES:
        recompress(path)
        print(path)


if __name__ == "__main__":
    main()
