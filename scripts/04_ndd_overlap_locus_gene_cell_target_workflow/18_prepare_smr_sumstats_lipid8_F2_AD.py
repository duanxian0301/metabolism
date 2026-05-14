from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
OUT = ROOT / "results" / "11_smr_lipid8_F2_AD" / "input" / "smr_sumstats"

TRAIT_INPUTS = {
    "lipid8_F2": ROOT / "work" / "clean_factor_inputs" / "lipid8_F2_clean.txt",
    "AD": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
}


def prepare_trait(trait: str, path: Path) -> Path:
    df = pd.read_csv(path, sep="\t", dtype={"CHR": str})
    required = ["SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{trait} missing columns: {missing}")
    out = df.loc[:, required].copy()
    out = out.rename(columns={"FREQ": "freq", "BETA": "b", "SE": "se", "P": "p", "N": "n"})
    out = out.dropna(subset=["SNP", "A1", "A2", "freq", "b", "se", "p", "n"])
    out = out.loc[out["SNP"].astype(str).str.startswith("rs")]
    out = out.loc[(out["freq"] > 0) & (out["freq"] < 1) & (out["se"] > 0) & (out["p"] > 0) & (out["p"] <= 1) & (out["n"] > 0)]
    out = out.drop_duplicates(subset=["SNP"], keep="first")
    out = out.loc[:, ["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]]
    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / f"{trait}.smr.txt"
    out.to_csv(out_path, sep="\t", index=False, float_format="%.10g")
    return out_path


def main() -> None:
    for trait, path in TRAIT_INPUTS.items():
        out = prepare_trait(trait, path)
        print(f"{trait}\t{out}")


if __name__ == "__main__":
    main()
