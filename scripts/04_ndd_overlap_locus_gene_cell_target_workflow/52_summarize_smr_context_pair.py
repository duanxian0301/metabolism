from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")


def read_table(path: Path, extra: dict[str, object]) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(path, sep=r"\s+")
    for key, value in extra.items():
        df.insert(len(df.columns) if key in df.columns else 0, key, value)
    return df


def summarize_gtex(pair: str, result_dir_name: str, traits: list[str], tissues: list[str]) -> None:
    base = ROOT / "results" / result_dir_name
    out_dir = base / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for trait in traits:
        for tissue in tissues:
            path = base / "single_trait" / trait / "gtex_v8_brain" / f"{trait}_{tissue}.smr"
            frames.append(read_table(path, {"trait": trait, "tissue": tissue}))
    frames = [x for x in frames if not x.empty]
    if not frames:
        print(f"No GTEx outputs found for {pair}")
        return
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_dir / f"smr_gtex_{pair}_combined.tsv", sep="\t", index=False)
    if "p_SMR" in combined.columns:
        top = combined.sort_values("p_SMR").groupby("trait", as_index=False).head(50)
        top.to_csv(out_dir / f"smr_gtex_{pair}_top50_by_trait.tsv", sep="\t", index=False)
        counts = (
            combined.assign(sig=pd.to_numeric(combined["p_SMR"], errors="coerce") < 5e-8)
            .groupby(["trait", "tissue"], as_index=False)
            .agg(n_genes=("Gene", "count"), n_sig=("sig", "sum"))
        )
        counts.to_csv(out_dir / f"smr_gtex_{pair}_counts.tsv", sep="\t", index=False)
    if "Gene" in combined.columns:
        shared = combined.groupby("Gene").filter(lambda x: x["trait"].nunique() == len(traits))
        shared.to_csv(out_dir / f"smr_gtex_{pair}_shared_genes.tsv", sep="\t", index=False)
    print(f"Wrote GTEx summary to {out_dir}")


def summarize_bryois(pair: str, result_dir_name: str, traits: list[str], celltypes: list[str], chroms: list[int]) -> None:
    base = ROOT / "results" / result_dir_name
    out_dir = base / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for trait in traits:
        for celltype in celltypes:
            for chrom in chroms:
                path = base / "single_trait" / trait / "bryois2022_celltype" / celltype / f"{trait}_{celltype}_chr{chrom}.smr"
                frames.append(read_table(path, {"trait": trait, "celltype": celltype, "chrom": chrom}))
    frames = [x for x in frames if not x.empty]
    if not frames:
        print(f"No Bryois outputs found for {pair}")
        return
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_dir / f"smr_bryois_{pair}_combined.tsv", sep="\t", index=False)
    if "p_SMR" in combined.columns:
        top = combined.sort_values("p_SMR").groupby("trait", as_index=False).head(100)
        top.to_csv(out_dir / f"smr_bryois_{pair}_top100_by_trait.tsv", sep="\t", index=False)
        counts = (
            combined.assign(sig=pd.to_numeric(combined["p_SMR"], errors="coerce") < 5e-8)
            .groupby(["trait", "celltype"], as_index=False)
            .agg(n_genes=("Gene", "count"), n_sig=("sig", "sum"))
        )
        counts.to_csv(out_dir / f"smr_bryois_{pair}_counts.tsv", sep="\t", index=False)
    if "Gene" in combined.columns:
        shared = combined.groupby(["Gene", "celltype"]).filter(lambda x: x["trait"].nunique() == len(traits))
        shared.to_csv(out_dir / f"smr_bryois_{pair}_shared_genes.tsv", sep="\t", index=False)
    print(f"Wrote Bryois summary to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)

    p_gtex = sub.add_parser("gtex")
    p_gtex.add_argument("--pair", required=True)
    p_gtex.add_argument("--result-dir-name", required=True)
    p_gtex.add_argument("--traits", nargs="+", required=True)
    p_gtex.add_argument("--tissues", nargs="+", required=True)

    p_bry = sub.add_parser("bryois")
    p_bry.add_argument("--pair", required=True)
    p_bry.add_argument("--result-dir-name", required=True)
    p_bry.add_argument("--traits", nargs="+", required=True)
    p_bry.add_argument("--celltypes", nargs="+", required=True)
    p_bry.add_argument("--chroms", nargs="+", type=int, required=True)

    args = parser.parse_args()
    if args.mode == "gtex":
        summarize_gtex(args.pair, args.result_dir_name, args.traits, args.tissues)
    else:
        summarize_bryois(args.pair, args.result_dir_name, args.traits, args.celltypes, args.chroms)


if __name__ == "__main__":
    main()
