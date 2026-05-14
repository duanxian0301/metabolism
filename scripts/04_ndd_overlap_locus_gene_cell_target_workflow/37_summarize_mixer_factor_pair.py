from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
FIELDS = [
    "rg",
    "rho_beta",
    "pi1",
    "pi2",
    "pi12",
    "pi1u",
    "pi2u",
    "dice",
    "nc1",
    "nc2",
    "nc12",
    "nc1u",
    "nc2u",
    "pi12_over_pi1u",
    "pi12_over_pi2u",
    "pi12_over_totalpi",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait1", required=True)
    parser.add_argument("--trait2", required=True)
    parser.add_argument("--run-name", default="fast_rep1")
    parser.add_argument("--rep", default="1")
    args = parser.parse_args()

    json_path = (
        ROOT
        / "work"
        / "04_mixer_inputs"
        / "runs"
        / args.run_name
        / f"{args.trait1}_vs_{args.trait2}.fit.fast.rep{args.rep}.json"
    )
    if not json_path.exists():
        alt = json_path.with_name(f"{args.trait1}_vs_{args.trait2}.fit.fast2.rep{args.rep}.json")
        if alt.exists():
            json_path = alt
        else:
            raise FileNotFoundError(json_path)

    out_dir = ROOT / "results" / f"05_mixer_{args.trait1}_{args.trait2}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"mixer_{args.trait1}_{args.trait2}_summary.tsv"

    data = json.loads(json_path.read_text())
    ci = data["ci"]
    row = {"trait1": args.trait1, "trait2": args.trait2, "source": str(json_path)}
    for field in FIELDS:
        row[field] = ci[field]["point_estimate"]

    with out_path.open("w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=list(row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    print(out_path)


if __name__ == "__main__":
    main()
