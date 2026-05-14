from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
RESULTS = ROOT / "results"


def as_bool01(value) -> bool:
    try:
        return int(float(value)) > 0
    except Exception:
        return False


def join_unique(values) -> str | None:
    vals = sorted({str(v) for v in values if pd.notna(v) and str(v) not in {"", "nan", "NA"}})
    return ";".join(vals) if vals else None


def genes_for_locus(fuma_genes: pd.DataFrame, genomic_locus, lead_snp: str) -> dict[str, str | int | None]:
    sub = pd.DataFrame()
    if pd.notna(genomic_locus):
        key = str(genomic_locus)
        if key.replace(".0", "").isdigit():
            key = str(int(float(key)))
        sub = fuma_genes[fuma_genes["GenomicLocus"].astype(str) == key].copy()
    if sub.empty:
        sub = fuma_genes[fuma_genes["IndSigSNPs"].fillna("").astype(str).str.contains(lead_snp, regex=False)].copy()
    if sub.empty:
        return {
            "fuma_mapped_gene_n": 0,
            "fuma_mapped_genes_any": None,
            "fuma_mapped_genes_pos": None,
            "fuma_mapped_genes_eqtl": None,
            "fuma_mapped_genes_ci": None,
        }
    return {
        "fuma_mapped_gene_n": int(sub["symbol"].nunique()),
        "fuma_mapped_genes_any": join_unique(sub["symbol"]),
        "fuma_mapped_genes_pos": join_unique(sub.loc[pd.to_numeric(sub["posMapSNPs"], errors="coerce").fillna(0) > 0, "symbol"]),
        "fuma_mapped_genes_eqtl": join_unique(sub.loc[pd.to_numeric(sub["eqtlMapSNPs"], errors="coerce").fillna(0) > 0, "symbol"]),
        "fuma_mapped_genes_ci": join_unique(sub.loc[sub["ciMap"].fillna("").astype(str).str.upper().eq("YES"), "symbol"]),
    }


def pick_best_region(df: pd.DataFrame, score_col: str) -> pd.Series | None:
    if df.empty:
        return None
    tmp = df.copy()
    tmp[score_col] = pd.to_numeric(tmp[score_col], errors="coerce")
    tmp = tmp.sort_values(score_col, ascending=False)
    return tmp.iloc[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", required=True)
    ap.add_argument("--factor-trait", required=True)
    ap.add_argument("--fuma-job-dir", required=True)
    args = ap.parse_args()

    pair = args.pair
    factor = args.factor_trait
    outdir = RESULTS / f"20_snp_evidence_{pair}"
    outdir.mkdir(parents=True, exist_ok=True)

    conj = pd.read_csv(RESULTS / f"06_pleiofdr_{pair}" / f"{pair}_conjfdr_0.05_loci.csv", low_memory=False)
    coloc = pd.read_csv(RESULTS / f"09_coloc_{pair}" / f"coloc_{pair}_regions.tsv", sep="\t", low_memory=False)
    pwcoco = pd.read_csv(RESULTS / f"10_pwcoco_{pair}" / f"coloc_pwcoco_{pair}_region_integrated.tsv", sep="\t", low_memory=False)
    fuma_snps = pd.read_csv(Path(args.fuma_job_dir) / "snps.txt", sep="\t", low_memory=False)
    fuma_genes = pd.read_csv(Path(args.fuma_job_dir) / "genes.txt", sep="\t", low_memory=False)

    for df, chr_col in [(conj, "chrnum"), (coloc, "chrnum"), (pwcoco, "chrnum"), (fuma_snps, "chr")]:
        df[chr_col] = pd.to_numeric(df[chr_col], errors="coerce")
    conj["chrpos"] = pd.to_numeric(conj["chrpos"], errors="coerce")
    coloc["region_start"] = pd.to_numeric(coloc["region_start"], errors="coerce")
    coloc["region_end"] = pd.to_numeric(coloc["region_end"], errors="coerce")
    pwcoco["region_start"] = pd.to_numeric(pwcoco["region_start"], errors="coerce")
    pwcoco["region_end"] = pd.to_numeric(pwcoco["region_end"], errors="coerce")
    fuma_snps["pos"] = pd.to_numeric(fuma_snps["pos"], errors="coerce")
    fuma_snps["CADD"] = pd.to_numeric(fuma_snps["CADD"], errors="coerce")
    fuma_snps["r2"] = pd.to_numeric(fuma_snps["r2"], errors="coerce")

    conj_col = f"conjfdr_{pair}"
    pval_col = f"pval_{factor}"
    fdr_col = f"fdr_{factor}"

    rows: list[dict] = []
    for _, r in conj.iterrows():
        chrnum = r["chrnum"]
        chrpos = r["chrpos"]
        snpid = str(r["snpid"])

        fuma_match = fuma_snps[fuma_snps["rsID"].astype(str) == snpid].copy()
        if fuma_match.empty:
            fuma_match = fuma_snps[(fuma_snps["chr"] == chrnum) & (fuma_snps["pos"] == chrpos)].copy()
        fuma_row = fuma_match.iloc[0] if not fuma_match.empty else None

        coloc_hits = coloc[(coloc["chrnum"] == chrnum) & (coloc["region_start"] <= chrpos) & (coloc["region_end"] >= chrpos)].copy()
        coloc_best = pick_best_region(coloc_hits, "PP.H4")

        pwcoco_hits = pwcoco[(pwcoco["chrnum"] == chrnum) & (pwcoco["region_start"] <= chrpos) & (pwcoco["region_end"] >= chrpos)].copy()
        pwcoco_best = pick_best_region(pwcoco_hits, "pwcoco_best_H4")

        genomic_locus = fuma_row["GenomicLocus"] if fuma_row is not None and "GenomicLocus" in fuma_row else None
        gene_info = genes_for_locus(fuma_genes, genomic_locus, snpid)

        out = {
            "pair": pair,
            "locusnum": r["locusnum"],
            "lead_snp": snpid,
            "chrnum": chrnum,
            "chrpos": chrpos,
            conj_col: r.get(conj_col),
            pval_col: r.get(pval_col),
            fdr_col: r.get(fdr_col),
            "coloc_region_id": None if coloc_best is None else coloc_best["region_id"],
            "coloc_class": None if coloc_best is None else coloc_best["coloc_class"],
            "coloc_PP.H4": None if coloc_best is None else coloc_best["PP.H4"],
            "coloc_status": None if coloc_best is None else coloc_best["status"],
            "coloc_sentinel_snps": None if coloc_best is None else coloc_best["sentinel_snps"],
            "coloc_n_overlapping_regions": int(len(coloc_hits)),
            "coloc_supported": False if coloc_best is None else str(coloc_best["coloc_class"]) in {"strong_H4", "moderate_H4"},
            "pwcoco_region_id": None if pwcoco_best is None else pwcoco_best["region_id"],
            "pwcoco_best_h4_class": None if pwcoco_best is None else pwcoco_best["pwcoco_best_h4_class"],
            "pwcoco_best_H4": None if pwcoco_best is None else pwcoco_best["pwcoco_best_H4"],
            "pwcoco_best_result_type": None if pwcoco_best is None else pwcoco_best["pwcoco_best_result_type"],
            "pwcoco_region_class": None if pwcoco_best is None else pwcoco_best["pwcoco_region_class"],
            "pwcoco_priority_shared_signal": None if pwcoco_best is None else pwcoco_best["priority_shared_signal"],
            "pwcoco_sentinel_snps": None if pwcoco_best is None else pwcoco_best["sentinel_snps"],
            "pwcoco_n_overlapping_regions": int(len(pwcoco_hits)),
            "pwcoco_supported": False if pwcoco_best is None else str(pwcoco_best["pwcoco_best_h4_class"]) in {"strong_H4", "moderate_H4"},
            "fuma_rsID": None if fuma_row is None else fuma_row["rsID"],
            "fuma_uniqID": None if fuma_row is None else fuma_row["uniqID"],
            "fuma_GenomicLocus": None if fuma_row is None else fuma_row["GenomicLocus"],
            "fuma_IndSigSNP": None if fuma_row is None else fuma_row["IndSigSNP"],
            "fuma_r2": None if fuma_row is None else fuma_row["r2"],
            "nearestGene": None if fuma_row is None else fuma_row["nearestGene"],
            "dist": None if fuma_row is None else fuma_row["dist"],
            "func": None if fuma_row is None else fuma_row["func"],
            "CADD": None if fuma_row is None else fuma_row["CADD"],
            "RDB": None if fuma_row is None else fuma_row["RDB"],
            "posMapFilt": False if fuma_row is None else as_bool01(fuma_row["posMapFilt"]),
            "eqtlMapFilt": False if fuma_row is None else as_bool01(fuma_row["eqtlMapFilt"]),
            "ciMapFilt": False if fuma_row is None else as_bool01(fuma_row["ciMapFilt"]),
            **gene_info,
        }
        out["fuma_any_map_supported"] = out["posMapFilt"] or out["eqtlMapFilt"] or out["ciMapFilt"]
        out["support_count_core"] = int(out["coloc_supported"]) + int(out["pwcoco_supported"]) + int(out["posMapFilt"]) + int(out["eqtlMapFilt"]) + int(out["ciMapFilt"])
        rows.append(out)

    evidence = pd.DataFrame(rows).sort_values(["support_count_core", conj_col, "chrnum", "chrpos"], ascending=[False, True, True, True])
    evidence.to_csv(outdir / f"{pair}_lead_snp_evidence_table.tsv", sep="\t", index=False)

    counts = pd.DataFrame(
        {
            "metric": [
                "n_lead_snps",
                "n_coloc_supported",
                "n_pwcoco_supported",
                "n_fuma_any_map_supported",
                "n_posMapFilt",
                "n_eqtlMapFilt",
                "n_ciMapFilt",
                "n_support_count_core_ge_3",
            ],
            "value": [
                len(evidence),
                int(evidence["coloc_supported"].sum()),
                int(evidence["pwcoco_supported"].sum()),
                int(evidence["fuma_any_map_supported"].sum()),
                int(evidence["posMapFilt"].sum()),
                int(evidence["eqtlMapFilt"].sum()),
                int(evidence["ciMapFilt"].sum()),
                int((evidence["support_count_core"] >= 3).sum()),
            ],
        }
    )
    counts.to_csv(outdir / f"{pair}_lead_snp_evidence_counts.tsv", sep="\t", index=False)

    with pd.ExcelWriter(outdir / f"{pair}_lead_snp_evidence_tables.xlsx") as writer:
        evidence.to_excel(writer, sheet_name="lead_snp_evidence", index=False)
        counts.to_excel(writer, sheet_name="counts", index=False)


if __name__ == "__main__":
    main()
