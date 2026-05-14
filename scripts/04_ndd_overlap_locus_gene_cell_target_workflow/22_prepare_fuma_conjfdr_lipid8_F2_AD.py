from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
PLEIO_DIR = ROOT / "results" / "06_pleiofdr_lipid8_F2_AD"
OUT_DIR = ROOT / "results" / "13_fuma_lipid8_F2_AD"

ALL_IN = PLEIO_DIR / "lipid8_F2_AD_conjfdr_0.05_all.csv"
LOCI_IN = PLEIO_DIR / "lipid8_F2_AD_conjfdr_0.05_loci.csv"


def load_conjfdr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(
        columns={
            "snpid": "SNP",
            "chrnum": "CHR",
            "chrpos": "BP",
            "pval_lipid8_F2": "P_lipid8_F2",
            "fdr_lipid8_F2": "FDR_lipid8_F2",
            "conjfdr_lipid8_F2_AD": "conjFDR_lipid8_F2_AD",
        }
    )
    required = ["SNP", "CHR", "BP", "P_lipid8_F2", "conjFDR_lipid8_F2_AD"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")
    df = df.dropna(subset=["SNP", "CHR", "BP", "conjFDR_lipid8_F2_AD"]).copy()
    df["SNP"] = df["SNP"].astype(str)
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce").astype("Int64")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce").astype("Int64")
    df["P"] = pd.to_numeric(df["conjFDR_lipid8_F2_AD"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP", "P"])
    df = df[(df["P"] > 0) & (df["P"] <= 0.05)]
    df = df.drop_duplicates(subset=["SNP"], keep="first")
    df = df.sort_values(["CHR", "BP", "SNP"], kind="stable")
    return df


def write_fuma(df: pd.DataFrame, prefix: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fuma = df[["SNP", "CHR", "BP", "P"]].copy()
    fuma.to_csv(OUT_DIR / f"{prefix}_for_FUMA_SNP2GENE.txt", sep="\t", index=False)

    # Optional richer file for manual checking and later linking back to pleioFDR evidence.
    trace_cols = [
        c
        for c in [
            "locusnum",
            "SNP",
            "geneid",
            "CHR",
            "BP",
            "P_lipid8_F2",
            "FDR_lipid8_F2",
            "conjFDR_lipid8_F2_AD",
            "min_conjfdr",
            "prune_lipid8_F2_AD",
        ]
        if c in df.columns
    ]
    df[trace_cols].to_csv(OUT_DIR / f"{prefix}_traceback.tsv", sep="\t", index=False)


def main() -> None:
    all_df = load_conjfdr(ALL_IN)
    loci_df = load_conjfdr(LOCI_IN)

    write_fuma(all_df, "lipid8_F2_AD_conjFDR0.05_all_positive_SNPs")
    write_fuma(loci_df, "lipid8_F2_AD_conjFDR0.05_sentinel_loci")

    readme = f"""FUMA SNP2GENE upload files for lipid8_F2 x AD conjFDR-positive SNPs.

Primary upload file:
  lipid8_F2_AD_conjFDR0.05_all_positive_SNPs_for_FUMA_SNP2GENE.txt

Columns:
  SNP: rsID
  CHR: GRCh37 chromosome
  BP: GRCh37 base-pair position
  P: conjFDR_lipid8_F2_AD, used as the FUMA p-value column for annotation/prioritization.

Important interpretation note:
  The P column is a conjFDR value, not a conventional single-trait GWAS P value. Use this as a
  target-SNP annotation input and describe it as conjFDR-prioritized SNP mapping.

Backup/sentinel-only upload file:
  lipid8_F2_AD_conjFDR0.05_sentinel_loci_for_FUMA_SNP2GENE.txt

Traceback tables retain the original pleioFDR columns for manuscript linkage.

Counts:
  all positive SNPs: {len(all_df)}
  sentinel/locus rows: {len(loci_df)}
"""
    (OUT_DIR / "README_FUMA_upload_lipid8_F2_AD.txt").write_text(readme, encoding="utf-8")
    print(f"out_dir\t{OUT_DIR}")
    print(f"n_all_positive_snps\t{len(all_df)}")
    print(f"n_sentinel_loci\t{len(loci_df)}")


if __name__ == "__main__":
    main()
