import numpy as np
import pandas as pd
from pathlib import Path


BASE_DIR = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\step1_ldsc_results")
EXCEL_PATH = Path(r"D:\metabolic\metabolite_FGWAS_selection_lists.xlsx")
SHEET_NAME = "Main_Zgt4_nonproportion"
PRUNE_THRESHOLD = 0.99


def load_inputs():
    s = pd.read_csv(BASE_DIR / "Main_Zgt4_nonproportion_S_matrix.csv", index_col=0)
    rg = pd.read_csv(BASE_DIR / "Main_Zgt4_nonproportion_rg_matrix.csv", index_col=0)
    i_mat = pd.read_csv(BASE_DIR / "Main_Zgt4_nonproportion_I_matrix.csv", index_col=0)
    trait_order = list(rg.columns)

    s.index = trait_order
    s.columns = trait_order
    rg.index = trait_order
    i_mat.index = trait_order
    i_mat.columns = trait_order

    manifest = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_NAME)
    manifest = manifest.rename(
        columns={
            "trait": "trait_code",
            "biomarker name": "biomarker_name",
            "study_accession": "study_accession",
        }
    )
    manifest["trait_code"] = manifest["trait_code"].astype(str)
    manifest = manifest.drop_duplicates(subset=["trait_code"])
    return s, rg, i_mat, manifest, trait_order


def greedy_prune(rg, h2, trait_order, threshold):
    keep = set(trait_order)
    steps = []

    while True:
        keep_list = [t for t in trait_order if t in keep]
        cur = rg.loc[keep_list, keep_list]
        upper = cur.where(np.triu(np.ones(cur.shape), 1).astype(bool)).stack()
        upper = upper[upper.abs() > threshold]
        if upper.empty:
            break

        t1, t2 = upper.abs().sort_values(ascending=False).index[0]
        rg_value = rg.loc[t1, t2]
        drop = t1 if h2[t1] < h2[t2] else t2
        keep.remove(drop)
        steps.append(
            {
                "trait1": t1,
                "trait2": t2,
                "rg": rg_value,
                "threshold": threshold,
                "dropped_trait": drop,
                "kept_trait": t2 if drop == t1 else t1,
                "h2_trait1": h2[t1],
                "h2_trait2": h2[t2],
                "drop_reason": f"|rg| > {threshold} and lower h2 within the pair",
            }
        )

    return [t for t in trait_order if t in keep], pd.DataFrame(steps)


def main():
    s, rg, i_mat, manifest, trait_order = load_inputs()
    h2 = pd.Series(np.diag(s.values), index=trait_order)
    intercept = pd.Series(np.diag(i_mat.values), index=trait_order)

    kept_traits, pruning_steps = greedy_prune(rg, h2, trait_order, PRUNE_THRESHOLD)
    dropped_traits = [t for t in trait_order if t not in kept_traits]

    summary = pd.DataFrame(
        {
            "trait_code": trait_order,
            "h2": h2[trait_order].values,
            "intercept": intercept[trait_order].values,
            "max_abs_rg_with_any_trait": rg.abs().where(
                np.triu(np.ones(rg.shape), 1).astype(bool)
            ).max(axis=1).reindex(trait_order).values,
            "selected_for_efa": [t in kept_traits for t in trait_order],
        }
    )
    summary = summary.merge(manifest, on="trait_code", how="left")

    kept_df = summary[summary["selected_for_efa"]].copy()
    dropped_df = summary[~summary["selected_for_efa"]].copy()

    if not pruning_steps.empty:
        dropped_df = dropped_df.merge(
            pruning_steps[["dropped_trait", "kept_trait", "rg", "drop_reason"]],
            left_on="trait_code",
            right_on="dropped_trait",
            how="left",
        )

    upper = rg.where(np.triu(np.ones(rg.shape), 1).astype(bool)).stack()
    matrix_diag = pd.DataFrame(
        [
            {
                "metric": "n_input_traits",
                "value": len(trait_order),
            },
            {
                "metric": "n_kept_for_efa",
                "value": len(kept_traits),
            },
            {
                "metric": "n_dropped_by_redundancy",
                "value": len(dropped_traits),
            },
            {
                "metric": "input_min_h2",
                "value": float(h2.min()),
            },
            {
                "metric": "input_max_h2",
                "value": float(h2.max()),
            },
            {
                "metric": "input_min_intercept",
                "value": float(intercept.min()),
            },
            {
                "metric": "input_max_intercept",
                "value": float(intercept.max()),
            },
            {
                "metric": "input_pairs_abs_rg_gt_0_99",
                "value": int((upper.abs() > 0.99).sum()),
            },
            {
                "metric": "input_pairs_abs_rg_gt_1",
                "value": int((upper.abs() > 1).sum()),
            },
            {
                "metric": "input_min_eigenvalue_S",
                "value": float(np.linalg.eigvalsh(s.values).min()),
            },
            {
                "metric": "kept_min_eigenvalue_S",
                "value": float(np.linalg.eigvalsh(s.loc[kept_traits, kept_traits].values).min()),
            },
        ]
    )

    out_dir = BASE_DIR / "efa_selection"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary.to_csv(out_dir / "trait_qc_and_selection.tsv", sep="\t", index=False)
    kept_df.to_csv(out_dir / "traits_kept_for_efa.tsv", sep="\t", index=False)
    dropped_df.to_csv(out_dir / "traits_dropped_for_efa.tsv", sep="\t", index=False)
    pruning_steps.to_csv(out_dir / "pruning_steps.tsv", sep="\t", index=False)
    matrix_diag.to_csv(out_dir / "selection_diagnostics.tsv", sep="\t", index=False)

    group_counts = (
        kept_df.groupby("group", dropna=False)["trait_code"]
        .count()
        .reset_index(name="n_kept")
        .sort_values(["n_kept", "group"], ascending=[False, True])
    )
    group_counts.to_csv(out_dir / "kept_group_counts.tsv", sep="\t", index=False)

    print(f"Output directory: {out_dir}")
    print(f"Input traits: {len(trait_order)}")
    print(f"Kept for EFA: {len(kept_traits)}")
    print(f"Dropped for redundancy: {len(dropped_traits)}")


if __name__ == "__main__":
    main()
