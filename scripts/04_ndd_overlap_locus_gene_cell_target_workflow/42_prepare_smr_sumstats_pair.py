from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def prepare_trait(trait: str, path: Path, out_dir: Path) -> Path:
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
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{trait}.smr.txt"
    out.to_csv(out_path, sep="\t", index=False, float_format="%.10g")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait1", required=True)
    parser.add_argument("--trait1-input", required=True)
    parser.add_argument("--trait2", required=True)
    parser.add_argument("--trait2-input", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    for trait, path in [
        (args.trait1, Path(args.trait1_input)),
        (args.trait2, Path(args.trait2_input)),
    ]:
        out = prepare_trait(trait, path, out_dir)
        print(f"{trait}\t{out}")


if __name__ == "__main__":
    main()
