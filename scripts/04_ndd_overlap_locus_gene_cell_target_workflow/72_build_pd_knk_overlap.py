from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
SCP_ROOT = Path(r"D:\scPagwas\metabolic_scpagwas2")


ANALYSES = [
    {
        "name": "nonlipid8_F1_PD",
        "knk_root": ROOT / "results" / "24_knk_nonlipid8_F1_PD",
        "summary_csv": ROOT
        / "results"
        / "24_knk_nonlipid8_F1_PD"
        / "summary"
        / "knk_summary_nonlipid8_F1_PD.csv",
        "scp_root": SCP_ROOT / "nonlipid8_F1_MSSM_PD",
        "out_root": ROOT / "results" / "26_knk_scpagwas_overlap_nonlipid8_F1_PD",
        "cell_file_map": {
            "Pericyte": "pericyte",
            "Oligodendrocyte_precursor_cell": "oligodendrocyte_precursor_cell",
            "VIP_GABAergic_cortical_interneuron": "VIP_GABAergic_cortical_interneuron",
        },
    },
    {
        "name": "lipid8_F1_PD",
        "knk_root": ROOT / "results" / "25_knk_lipid8_F1_PD",
        "summary_csv": ROOT
        / "results"
        / "25_knk_lipid8_F1_PD"
        / "summary"
        / "knk_summary_lipid8_F1_PD.csv",
        "scp_root": SCP_ROOT / "lipid8_F1_MSSM_PD",
        "out_root": ROOT / "results" / "27_knk_scpagwas_overlap_lipid8_F1_PD",
        "cell_file_map": {
            "Pericyte": "pericyte",
            # lipid8_F1-PD scPagwas2 only yielded pericyte pathway TRS outputs.
        },
    },
]


def _read_pathway_ids(path: Path, top_n: int) -> list[str]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    df = pd.read_csv(path)
    if df.empty:
        return []
    for col in ("ID", "pathway_id"):
        if col in df.columns:
            ids = df[col].dropna().astype(str).tolist()
            return ids[:top_n]
    return []


def _jaccard(a: Iterable[str], b: Iterable[str]) -> float:
    set_a = set(a)
    set_b = set(b)
    union = set_a | set_b
    if not union:
        return math.nan
    return len(set_a & set_b) / len(union)


def _find_knk_pathway_file(knk_root: Path, cell_type: str, gene: str) -> Path | None:
    matches = list(knk_root.glob(f"{cell_type}/**/*{gene}*_pathway_enrichment.csv"))
    if not matches:
        return None
    return matches[0]


def _refresh_summary_top_pathways(summary_csv: Path, knk_root: Path) -> pd.DataFrame:
    summary_df = pd.read_csv(summary_csv)
    summary_df.columns = [c.strip() for c in summary_df.columns]
    updated = False
    for idx, row in summary_df.iterrows():
        if str(row.get("ko_success", "")).strip().lower() not in {"true", "1"}:
            continue
        current_raw = row.get("top_pathways", "")
        current = "" if pd.isna(current_raw) else str(current_raw).strip()
        if current and current.lower() != "nan":
            continue
        pathway_file = _find_knk_pathway_file(knk_root, str(row["cell_type"]), str(row["gene"]))
        top_ids = _read_pathway_ids(pathway_file, 5) if pathway_file else []
        if top_ids:
            summary_df.at[idx, "top_pathways"] = "; ".join(top_ids)
            updated = True
    if updated:
        summary_df.to_csv(summary_csv, index=False)
    return summary_df


def _write_scp_top_tables(pathway_trs_dir: Path, out_root: Path, file_stem: str) -> tuple[list[str], list[str]]:
    all_file = pathway_trs_dir / f"Result_{file_stem}_Pathway_vs_TRS_all.csv"
    top20 = _read_pathway_ids(all_file, 20)
    top50 = _read_pathway_ids(all_file, 50)
    if all_file.exists():
        df = pd.read_csv(all_file)
        df.head(20).to_csv(out_root / f"scpagwas_{file_stem}_top20.tsv", sep="\t", index=False)
        df.head(50).to_csv(out_root / f"scpagwas_{file_stem}_top50.tsv", sep="\t", index=False)
    return top20, top50


def build_one(meta: dict) -> None:
    out_root = meta["out_root"]
    out_root.mkdir(parents=True, exist_ok=True)

    summary_df = _refresh_summary_top_pathways(meta["summary_csv"], meta["knk_root"])

    scp_top_sets: dict[str, tuple[list[str], list[str]]] = {}
    pathway_trs_dir = meta["scp_root"] / "Pathway_TRS"
    for cell_type, file_stem in meta["cell_file_map"].items():
        scp_top_sets[cell_type] = _write_scp_top_tables(pathway_trs_dir, out_root, file_stem)

    rows = []
    for _, row in summary_df.iterrows():
        if str(row.get("ko_success", "")).strip().lower() not in {"true", "1"}:
            continue
        if pd.isna(row.get("perturbed_gene_n")) or float(row.get("perturbed_gene_n", 0)) <= 0:
            continue

        cell_type = str(row["cell_type"])
        gene = str(row["gene"])
        pathway_file = _find_knk_pathway_file(meta["knk_root"], cell_type, gene)
        knk_top20 = _read_pathway_ids(pathway_file, 20) if pathway_file else []
        knk_top50 = _read_pathway_ids(pathway_file, 50) if pathway_file else []

        record = {
            "cell_type": cell_type,
            "gene": gene,
            "knk_n_pathways": len(_read_pathway_ids(pathway_file, 100000)) if pathway_file else 0,
        }

        if cell_type not in scp_top_sets:
            record.update(
                {
                    "scpagwas_match_available": False,
                    "scpagwas_reference_cell": "",
                    "knk_top20_vs_cell_trs_top20_overlap_n": math.nan,
                    "knk_top20_vs_cell_trs_top20_overlap_ids": "",
                    "knk_top20_vs_cell_trs_top20_jaccard": math.nan,
                    "knk_top50_vs_cell_trs_top50_overlap_n": math.nan,
                    "knk_top50_vs_cell_trs_top50_overlap_ids": "",
                    "knk_top50_vs_cell_trs_top50_jaccard": math.nan,
                    "remarks": "No matching scPagwas2 Pathway_TRS file for this cell type",
                }
            )
            rows.append(record)
            continue

        scp20, scp50 = scp_top_sets[cell_type]
        overlap20 = [x for x in knk_top20 if x in set(scp20)]
        overlap50 = [x for x in knk_top50 if x in set(scp50)]

        record.update(
            {
                "scpagwas_match_available": True,
                "scpagwas_reference_cell": meta["cell_file_map"][cell_type],
                "knk_top20_vs_cell_trs_top20_overlap_n": len(overlap20),
                "knk_top20_vs_cell_trs_top20_overlap_ids": ";".join(overlap20),
                "knk_top20_vs_cell_trs_top20_jaccard": _jaccard(knk_top20, scp20),
                "knk_top50_vs_cell_trs_top50_overlap_n": len(overlap50),
                "knk_top50_vs_cell_trs_top50_overlap_ids": ";".join(overlap50),
                "knk_top50_vs_cell_trs_top50_jaccard": _jaccard(knk_top50, scp50),
                "remarks": "",
            }
        )
        rows.append(record)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_root / f"knk_vs_scpagwas_{meta['name']}_overlap.tsv", sep="\t", index=False)


def main() -> None:
    for meta in ANALYSES:
        build_one(meta)


if __name__ == "__main__":
    main()
