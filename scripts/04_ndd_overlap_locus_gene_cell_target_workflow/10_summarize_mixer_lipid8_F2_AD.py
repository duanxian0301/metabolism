from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
MIXER_JSON = (
    ROOT
    / "work"
    / "04_mixer_inputs"
    / "runs"
    / "fast_rep1"
    / "lipid8_F2_vs_AD.fit.fast2.rep1.json"
)
OUT = ROOT / "results" / "05_mixer_lipid8_F2_AD" / "mixer_lipid8_F2_AD_summary.tsv"

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
    data = json.loads(MIXER_JSON.read_text())
    ci = data["ci"]
    row = {"trait1": "lipid8_F2", "trait2": "AD", "source": str(MIXER_JSON)}
    for field in FIELDS:
        row[field] = ci[field]["point_estimate"]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=list(row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
