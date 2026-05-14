from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")


def read_smr(base: Path, trait: str) -> pd.DataFrame:
    path = base / "single_trait" / trait / "brainmeta" / f"{trait}_BrainMeta.smr"
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(path, sep=r"\s+")
    df.insert(0, "trait", trait)
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair", required=True)
    parser.add_argument("--result-dir-name", required=True)
    parser.add_argument("--traits", nargs="+", required=True)
    args = parser.parse_args()

    base = ROOT / "results" / args.result_dir_name
    frames = [read_smr(base, trait) for trait in args.traits]
    frames = [x for x in frames if not x.empty]
    out_dir = base / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    if not frames:
        print("No SMR outputs found")
        return
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_dir / f"smr_brainmeta_{args.pair}_combined.tsv", sep="\t", index=False)
    p_col = "p_SMR" if "p_SMR" in combined.columns else None
    if p_col:
        top = combined.sort_values(p_col).groupby("trait", as_index=False).head(50)
        top.to_csv(out_dir / f"smr_brainmeta_{args.pair}_top50_by_trait.tsv", sep="\t", index=False)
    if "Gene" in combined.columns:
        shared = combined.groupby("Gene").filter(lambda x: x["trait"].nunique() == len(args.traits))
        shared.to_csv(out_dir / f"smr_brainmeta_{args.pair}_shared_genes.tsv", sep="\t", index=False)
    print(f"Wrote SMR summary to {out_dir}")
