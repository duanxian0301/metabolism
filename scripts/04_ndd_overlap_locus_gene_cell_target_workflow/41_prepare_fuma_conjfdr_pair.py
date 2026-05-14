from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_conjfdr(path: Path, trait1: str, disease: str) -> pd.DataFrame:
    pair = f"{trait1}_{disease}"
    df = pd.read_csv(path)
    df = df.rename(
        columns={
            "snpid": "SNP",
            "chrnum": "CHR",
            "chrpos": "BP",
            f"pval_{trait1}": f"P_{trait1}",
            f"fdr_{trait1}": f"FDR_{trait1}",
            f"conjfdr_{pair}": f"conjFDR_{pair}",
        }
    )
    required = ["SNP", "CHR", "BP", f"P_{trait1}", f"conjFDR_{pair}"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")
    df = df.dropna(subset=["SNP", "CHR", "BP", f"conjFDR_{pair}"]).copy()
    df["SNP"] = df["SNP"].astype(str)
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce").astype("Int64")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce").astype("Int64")
    df["P"] = pd.to_numeric(df[f"conjFDR_{pair}"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP", "P"])
    df = df[(df["P"] > 0) & (df["P"] <= 0.05)]
    df = df.drop_duplicates(subset=["SNP"], keep="first")
    df = df.sort_values(["CHR", "BP", "SNP"], kind="stable")
    return df


def write_fuma(df: pd.DataFrame, prefix: str, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fuma = df[["SNP", "CHR", "BP", "P"]].copy()
    fuma.to_csv(out_dir / f"{prefix}_for_FUMA_SNP2GENE.txt", sep="\t", index=False)
    trace_cols = [
        c
        for c in [
            "locusnum",
            "SNP",
            "geneid",
            "CHR",
            "BP",
            "min_conjfdr",
            "P",
        ]
        if c in df.columns
    ]
    df[trace_cols].to_csv(out_dir / f"{prefix}_traceback.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", required=True)
    parser.add_argument("--disease", required=True)
    parser.add_argument("--pleio-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    pair = f"{args.trait}_{args.disease}"
    pleio_dir = Path(args.pleio_dir)
    out_dir = Path(args.out_dir)
    all_df = load_conjfdr(pleio_dir / f"{pair}_conjfdr_0.05_all.csv", args.trait, args.disease)
    loci_df = load_conjfdr(pleio_dir / f"{pair}_conjfdr_0.05_loci.csv", args.trait, args.disease)
    write_fuma(all_df, f"{pair}_conjFDR0.05_all_positive_SNPs", out_dir)
    write_fuma(loci_df, f"{pair}_conjFDR0.05_sentinel_loci", out_dir)
    print(f"all_positive_snps\t{len(all_df)}")
    print(f"sentinel_loci\t{len(loci_df)}")


if __name__ == "__main__":
    main()
