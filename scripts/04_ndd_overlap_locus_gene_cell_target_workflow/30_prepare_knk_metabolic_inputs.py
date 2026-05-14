import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


SEED = 123
np.random.seed(SEED)

H5AD_PATH = Path("D:/codex/GenomicSEM/data_100k/MSSM_AD_noPD_noLBD.h5ad")
OUT_ROOT = Path("D:/codex/GenomicSEM/metabolic/postgwas_ad_pdlbd/work/knk_inputs")
DATASET = "MSSM_AD"

TARGETS = {
    "Pericyte": {
        "cell_types": ["pericyte"],
        "max_cells": 4000,
    },
    "Oligodendrocyte_precursor_cell": {
        "cell_types": ["oligodendrocyte precursor cell"],
        "max_cells": 4000,
    },
}


def write_lines_gz(path: Path, values) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for val in values:
            fh.write(f"{val}\n")


def write_mtx_gz(path: Path, matrix) -> None:
    tmp = path.with_suffix("")
    mmwrite(tmp.as_posix(), matrix)
    with open(tmp, "rb") as src, gzip.open(path, "wb") as dst:
        dst.writelines(src)
    tmp.unlink()


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"[prepare] reading {H5AD_PATH}", flush=True)
    adata = ad.read_h5ad(H5AD_PATH, backed="r")
    obs = adata.obs.copy()
    obs["cell_type"] = obs["cell_type"].astype(str)

    feature_col = "feature_name" if "feature_name" in adata.var.columns else None
    feature_symbols = (
        adata.var[feature_col].astype(str).tolist()
        if feature_col is not None
        else adata.var_names.astype(str).tolist()
    )
    feature_ids = adata.var_names.astype(str).tolist()

    manifest_rows = []
    for group_name, spec in TARGETS.items():
        cell_types = spec["cell_types"]
        max_cells = spec["max_cells"]
        group_dir = OUT_ROOT / DATASET / group_name
        group_dir.mkdir(parents=True, exist_ok=True)

        idx = np.where(obs["cell_type"].isin(cell_types).values)[0]
        available_n = int(len(idx))
        row = {
            "dataset": DATASET,
            "group": group_name,
            "cell_types": "|".join(cell_types),
            "available_cells": available_n,
            "sampled_cells": 0,
            "status": "PENDING",
            "note": "",
        }
        if available_n == 0:
            row["status"] = "SKIP_NO_CELLTYPE"
            row["note"] = "No matching cells in source h5ad"
            manifest_rows.append(row)
            print(f"[prepare] {group_name}: no matching cells", flush=True)
            continue

        if available_n > max_cells:
            idx = np.random.choice(idx, size=max_cells, replace=False)
        sampled_n = int(len(idx))
        row["sampled_cells"] = sampled_n
        print(
            f"[prepare] {group_name}: available={available_n}, sampled={sampled_n}",
            flush=True,
        )

        sub = adata[idx, :].to_memory()
        X = sub.X
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        X = X.tocsc()

        metadata = sub.obs.copy()
        metadata.index = metadata.index.astype(str)

        write_mtx_gz(group_dir / "matrix.mtx.gz", X.T)
        write_lines_gz(group_dir / "barcodes.tsv.gz", metadata.index.tolist())
        with gzip.open(group_dir / "features.tsv.gz", "wt", encoding="utf-8") as fh:
            for gid, gsym in zip(feature_ids, feature_symbols):
                fh.write(f"{gid}\t{gsym}\n")
        metadata.to_csv(group_dir / "metadata.csv.gz", compression="gzip")
        with open(group_dir / "manifest.json", "w", encoding="utf-8") as fh:
            json.dump(
                {
                    "dataset": DATASET,
                    "group": group_name,
                    "cell_types": cell_types,
                    "available_cells": available_n,
                    "sampled_cells": sampled_n,
                    "feature_col": feature_col,
                },
                fh,
                ensure_ascii=False,
                indent=2,
            )

        row["status"] = "OK"
        manifest_rows.append(row)

    manifest_path = OUT_ROOT / "metabolic_knk_input_manifest.csv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
    print(f"[prepare] wrote manifest: {manifest_path}", flush=True)


if __name__ == "__main__":
    main()
