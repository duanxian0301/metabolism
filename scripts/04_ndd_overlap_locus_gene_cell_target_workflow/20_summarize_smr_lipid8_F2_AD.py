from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
BASE = ROOT / "results" / "11_smr_lipid8_F2_AD"


def read_smr(trait: str) -> pd.DataFrame:
    path = BASE / "single_trait" / trait / "brainmeta" / f"{trait}_BrainMeta.smr"
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(path, sep=r"\s+")
    df.insert(0, "trait", trait)
    return df


def main() -> None:
    frames = [read_smr("lipid8_F2"), read_smr("AD")]
    frames = [x for x in frames if not x.empty]
    out_dir = BASE / "summary"
    out_dir.mkdir(parents=True, exist_ok=True)
    if not frames:
        print("No SMR outputs found")
        return
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_dir / "smr_brainmeta_lipid8_F2_AD_combined.tsv", sep="\t", index=False)
    p_col = "p_SMR" if "p_SMR" in combined.columns else None
    if p_col:
        top = combined.sort_values(p_col).groupby("trait", as_index=False).head(50)
        top.to_csv(out_dir / "smr_brainmeta_lipid8_F2_AD_top50_by_trait.tsv", sep="\t", index=False)
    if "Gene" in combined.columns:
        shared = combined.groupby("Gene").filter(lambda x: x["trait"].nunique() == 2)
        shared.to_csv(out_dir / "smr_brainmeta_lipid8_F2_AD_shared_genes.tsv", sep="\t", index=False)
    print(f"Wrote SMR summary to {out_dir}")


if __name__ == "__main__":
    main()
