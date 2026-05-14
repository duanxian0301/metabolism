from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd")
RESULTS = ROOT / "results"
OUTDIR = RESULTS / "17_candidate_gene_integration_lipid8_F2_AD"

FUMA_DIR = RESULTS / "13_fuma_lipid8_F2_AD" / "fuma" / "FUMA_job729202"
PLEIO_DIR = RESULTS / "06_pleiofdr_lipid8_F2_AD"
COLOC_DIR = RESULTS / "09_coloc_lipid8_F2_AD"
PWCOCO_DIR = RESULTS / "10_pwcoco_lipid8_F2_AD"
SMR_BRAINMETA_DIR = RESULTS / "11_smr_lipid8_F2_AD" / "summary"
SMR_GTEX_DIR = RESULTS / "14_smr_gtex_lipid8_F2_AD" / "summary"
SMR_BRYOIS_DIR = RESULTS / "15_smr_bryois_lipid8_F2_AD" / "summary"
CTWAS_DIR = RESULTS / "16_ctwas_lipid8_F2_AD" / "summary"


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce")
    out = pd.Series(index=p.index, dtype=float)
    valid = p.dropna().sort_values()
    n = len(valid)
    if n == 0:
        return out
    ranks = pd.Series(range(1, n + 1), index=valid.index, dtype=float)
    q = valid * n / ranks
    q = q.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    out.loc[q.index] = q
    return out


def norm_sf_twosided_from_z(z: pd.Series) -> pd.Series:
    znum = pd.to_numeric(z, errors="coerce").abs()
    return znum.map(lambda x: math.erfc(x / math.sqrt(2.0)) if pd.notna(x) else math.nan)


def parse_ctwas_gene(feature: str) -> str | None:
    if pd.isna(feature):
        return None
    feature = str(feature)
    if feature.startswith("rs"):
        return None
    left = feature.split("|", 1)[0]
    if "." not in left:
        return None
    return left.split(".", 1)[1]


def parse_ctwas_feature_context(feature: str) -> tuple[str | None, str | None]:
    if pd.isna(feature):
        return None, None
    feature = str(feature)
    if feature.startswith("rs"):
        return None, None
    left = feature.split("|", 1)[0]
    if "." not in left:
        return None, None
    context, gene = left.split(".", 1)
    return gene, context


def best_smr_by_gene(path: Path, context_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["p_SMR"] = pd.to_numeric(df["p_SMR"], errors="coerce")
    if "p_HEIDI" in df.columns:
        df["p_HEIDI"] = pd.to_numeric(df["p_HEIDI"], errors="coerce")
    df["fdr_SMR"] = bh_fdr(df["p_SMR"])
    heidi_ok = df["p_HEIDI"].isna() | (df["p_HEIDI"] > 0.01)
    df["supported"] = (df["fdr_SMR"] < 0.05) & heidi_ok
    best = (
        df.sort_values(["Gene", "supported", "p_SMR"], ascending=[True, False, True])
        .groupby("Gene", as_index=False)
        .first()
    )
    best = best.rename(
        columns={
            context_col: "best_context",
            "topSNP": "best_topSNP",
            "p_SMR": "best_p_SMR",
            "p_HEIDI": "best_p_HEIDI",
            "fdr_SMR": "best_fdr_SMR",
            "supported": "supported",
        }
    )
    return best[
        [
            "Gene",
            "best_context",
            "best_topSNP",
            "best_p_SMR",
            "best_p_HEIDI",
            "best_fdr_SMR",
            "supported",
        ]
    ]


def summarize_bryois_shared_by_gene(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["p_SMR"] = pd.to_numeric(df["p_SMR"], errors="coerce")
    if "p_HEIDI" in df.columns:
        df["p_HEIDI"] = pd.to_numeric(df["p_HEIDI"], errors="coerce")
    df["fdr_SMR"] = bh_fdr(df["p_SMR"])
    heidi_ok = df["p_HEIDI"].isna() | (df["p_HEIDI"] > 0.01)
    df["supported"] = (df["fdr_SMR"] < 0.05) & heidi_ok

    any_best = (
        df.sort_values(["Gene", "supported", "p_SMR"], ascending=[True, False, True])
        .groupby("Gene", as_index=False)
        .first()[
            ["Gene", "celltype", "topSNP", "p_SMR", "p_HEIDI", "fdr_SMR", "supported"]
        ]
        .rename(
            columns={
                "celltype": "any_best_context",
                "topSNP": "any_best_topSNP",
                "p_SMR": "any_best_p_SMR",
                "p_HEIDI": "any_best_p_HEIDI",
                "fdr_SMR": "any_best_fdr_SMR",
                "supported": "any_supported",
            }
        )
    )

    shared = (
        df[df["supported"]]
        .pivot_table(
            index=["Gene", "celltype"],
            columns="trait",
            values="p_SMR",
            aggfunc="min",
        )
        .reset_index()
    )
    for col in ["AD", "lipid8_F2"]:
        if col not in shared.columns:
            shared[col] = pd.NA
    shared["shared_supported"] = shared["AD"].notna() & shared["lipid8_F2"].notna()
    shared = shared[shared["shared_supported"]].copy()

    if shared.empty:
        shared_best = pd.DataFrame(
            columns=[
                "Gene",
                "shared_celltypes",
                "shared_n_celltypes",
                "shared_best_context",
                "shared_best_p_SMR_AD",
                "shared_best_p_SMR_lipid8_F2",
                "shared_best_joint_score",
                "shared_supported",
            ]
        )
    else:
        shared["joint_score"] = shared[["AD", "lipid8_F2"]].max(axis=1)
        shared_celltypes = (
            shared.groupby("Gene", as_index=False)
            .agg(
                shared_celltypes=("celltype", lambda x: ";".join(sorted(set(map(str, x))))),
                shared_n_celltypes=("celltype", "nunique"),
            )
        )
        shared_best = (
            shared.sort_values(["Gene", "joint_score", "AD", "lipid8_F2", "celltype"], ascending=[True, True, True, True, True])
            .groupby("Gene", as_index=False)
            .first()[
                ["Gene", "celltype", "AD", "lipid8_F2", "joint_score"]
            ]
            .rename(
                columns={
                    "celltype": "shared_best_context",
                    "AD": "shared_best_p_SMR_AD",
                    "lipid8_F2": "shared_best_p_SMR_lipid8_F2",
                    "joint_score": "shared_best_joint_score",
                }
            )
        )
        shared_best = shared_celltypes.merge(shared_best, on="Gene", how="left")
        shared_best["shared_supported"] = True

    return any_best.merge(shared_best, on="Gene", how="left")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fuma_genes = pd.read_csv(FUMA_DIR / "genes.txt", sep="\t")
    fuma_gene_table = pd.read_csv(
        RESULTS / "13_fuma_lipid8_F2_AD" / "fuma" / "FUMA_gene2func729262" / "geneTable.txt",
        sep="\t",
    )
    magma = pd.read_csv(FUMA_DIR / "magma.genes.out", sep=r"\s+")
    coloc = pd.read_csv(COLOC_DIR / "coloc_lipid8_F2_AD_regions.tsv", sep="\t")
    pwcoco_regions = pd.read_csv(PWCOCO_DIR / "coloc_pwcoco_lipid8_F2_AD_region_integrated.tsv", sep="\t")
    conj_loci = pd.read_csv(PLEIO_DIR / "lipid8_F2_AD_conjfdr_0.05_loci.csv")

    smr_brainmeta = best_smr_by_gene(
        SMR_BRAINMETA_DIR / "smr_brainmeta_lipid8_F2_AD_combined.tsv", "trait"
    )
    # brainmeta combined has no context column; use probe tissue placeholder from probe source
    if "best_context" in smr_brainmeta.columns:
        smr_brainmeta["best_context"] = "BrainMeta_bulk"

    smr_gtex = best_smr_by_gene(
        SMR_GTEX_DIR / "smr_gtex_lipid8_F2_AD_combined.tsv", "tissue"
    )
    smr_bryois = summarize_bryois_shared_by_gene(
        SMR_BRYOIS_DIR / "smr_bryois_lipid8_F2_AD_combined.tsv"
    )

    ctwas_z = []
    ctwas_pip = []
    for trait in ["AD", "lipid8_F2"]:
        zdf = pd.read_csv(CTWAS_DIR / f"{trait}_z_gene.tsv", sep="\t")
        zdf["Gene"] = zdf["id"].map(parse_ctwas_gene)
        zdf["context"] = zdf["id"].map(lambda x: parse_ctwas_feature_context(x)[1])
        zdf = zdf.dropna(subset=["Gene"]).copy()
        zdf["trait"] = trait
        zdf["abs_z"] = pd.to_numeric(zdf["z"], errors="coerce").abs()
        zdf["p_ctwas"] = norm_sf_twosided_from_z(zdf["z"])
        zdf["fdr_ctwas"] = bh_fdr(zdf["p_ctwas"])
        ctwas_z.append(zdf)

        pdf = pd.read_csv(CTWAS_DIR / f"{trait}_finemap_res.tsv", sep="\t")
        pdf = pdf[pdf["type"] != "SNP"].copy()
        pdf["Gene"] = pdf["molecular_id"].map(parse_ctwas_gene)
        pdf["context"] = pdf["molecular_id"].map(lambda x: str(x).split(".", 1)[0] if pd.notna(x) and "." in str(x) else None)
        pdf["trait"] = trait
        pdf["susie_pip"] = pd.to_numeric(pdf["susie_pip"], errors="coerce")
        ctwas_pip.append(pdf)
    ctwas_z = pd.concat(ctwas_z, ignore_index=True)
    ctwas_pip = pd.concat(ctwas_pip, ignore_index=True)

    ctw_z_best = (
        ctwas_z.sort_values(["Gene", "fdr_ctwas", "abs_z"], ascending=[True, True, False])
        .groupby("Gene", as_index=False)
        .first()[["Gene", "id", "context", "abs_z", "p_ctwas", "fdr_ctwas"]]
        .rename(columns={"id": "best_feature", "context": "best_context"})
    )
    ctw_pip_best = (
        ctwas_pip.sort_values(["Gene", "susie_pip"], ascending=[True, False])
        .groupby("Gene", as_index=False)
        .first()[["Gene", "molecular_id", "context", "susie_pip", "region_id"]]
        .rename(
            columns={
                "molecular_id": "top_pip_feature",
                "context": "top_pip_context",
                "susie_pip": "max_pip",
                "region_id": "top_pip_region_id",
            }
        )
    )
    ctwas_best = ctw_z_best.merge(ctw_pip_best, on="Gene", how="outer")
    ctwas_best["ctwas_supported"] = (ctwas_best["fdr_ctwas"] < 0.05) | (ctwas_best["max_pip"] >= 0.5)

    magma["magma_gene_supported"] = magma["P"] < (0.05 / len(magma))
    magma = magma.rename(columns={"SYMBOL": "Gene", "P": "magma_gene_p", "ZSTAT": "magma_gene_z"})

    # Build coordinate backbone from FUMA and MAGMA
    coords = pd.concat(
        [
            fuma_genes[["symbol", "chr", "start", "end"]].rename(
                columns={"symbol": "Gene", "chr": "chrnum", "start": "gene_start", "end": "gene_end"}
            ),
            magma[["Gene", "CHR", "START", "STOP"]].rename(
                columns={"CHR": "chrnum", "START": "gene_start", "STOP": "gene_end"}
            ),
        ],
        ignore_index=True,
    )
    coords = (
        coords.dropna(subset=["Gene", "chrnum", "gene_start", "gene_end"])
        .groupby("Gene", as_index=False)
        .first()
    )

    gene_universe = sorted(
        set(fuma_genes["symbol"].dropna())
        | set(magma["Gene"].dropna())
        | set(smr_brainmeta["Gene"].dropna())
        | set(smr_gtex["Gene"].dropna())
        | set(smr_bryois["Gene"].dropna())
        | set(ctwas_best["Gene"].dropna())
    )
    master = pd.DataFrame({"Gene": gene_universe})
    master["trait_pair"] = "lipid8_F2_AD"
    master = master.merge(coords, on="Gene", how="left")

    # FUMA aggregation
    fuma_genes["brain_eqtl_supported"] = fuma_genes["eqtlMapts"].fillna("").str.contains("Brain_|PsychENCODE|BRAINEAC|CMC", regex=True)
    fuma_agg = (
        fuma_genes.assign(
            fuma_supported=True,
            fuma_posMap_supported=fuma_genes["posMapSNPs"].fillna(0).astype(float) > 0,
            fuma_eqtlMap_supported=fuma_genes["eqtlMapSNPs"].fillna(0).astype(float) > 0,
            fuma_ciMap_supported=fuma_genes["ciMap"].fillna("").eq("Yes"),
        )
        .groupby("symbol", as_index=False)
        .agg(
            fuma_supported=("fuma_supported", "max"),
            fuma_posMap_supported=("fuma_posMap_supported", "max"),
            fuma_eqtlMap_supported=("fuma_eqtlMap_supported", "max"),
            fuma_ciMap_supported=("fuma_ciMap_supported", "max"),
            fuma_brain_eqtl_supported=("brain_eqtl_supported", "max"),
            fuma_eqtl_tissues=("eqtlMapts", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v) and str(v) != "NA"}))[:2000]),
            fuma_indsig_snps=("IndSigSNPs", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v)}))),
            fuma_genomic_loci=("GenomicLocus", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v) and str(v) != "NA"}))),
            fuma_min_gwas_p=("minGwasP", "min"),
        )
        .rename(columns={"symbol": "Gene"})
    )
    master = master.merge(fuma_agg, on="Gene", how="left")

    master = master.merge(magma[["Gene", "magma_gene_p", "magma_gene_z", "magma_gene_supported"]], on="Gene", how="left")

    # Overlap coloc/pwcoco regions with gene coordinates
    region_cols = ["region_id", "chrnum", "region_start", "region_end"]
    region_gene = master[["Gene", "chrnum", "gene_start", "gene_end"]].dropna().copy()
    region_gene["chrnum"] = pd.to_numeric(region_gene["chrnum"], errors="coerce")
    coloc["chrnum"] = pd.to_numeric(coloc["chrnum"], errors="coerce")
    pwcoco_regions["chrnum"] = pd.to_numeric(pwcoco_regions["chrnum"], errors="coerce")

    region_gene = region_gene.dropna(subset=["chrnum", "gene_start", "gene_end"])

    coloc_links = region_gene.merge(coloc, on="chrnum", how="inner")
    coloc_links = coloc_links[
        (coloc_links["gene_start"] <= coloc_links["region_end"]) & (coloc_links["gene_end"] >= coloc_links["region_start"])
    ].copy()
    coloc_links["coloc_supported"] = coloc_links["coloc_class"].isin(["strong_H4", "moderate_H4"])
    coloc_best = (
        coloc_links.sort_values(["Gene", "coloc_supported", "PP.H4"], ascending=[True, False, False])
        .groupby("Gene", as_index=False)
        .first()[["Gene", "region_id", "coloc_class", "PP.H4", "coloc_supported"]]
        .rename(columns={"region_id": "coloc_best_region", "PP.H4": "coloc_best_PP4"})
    )
    master = master.merge(coloc_best, on="Gene", how="left")

    pw_links = region_gene.merge(pwcoco_regions, on="chrnum", how="inner")
    pw_links = pw_links[
        (pw_links["gene_start"] <= pw_links["region_end"]) & (pw_links["gene_end"] >= pw_links["region_start"])
    ].copy()
    pw_links["pwcoco_supported"] = pw_links["pwcoco_best_h4_class"].isin(["strong_H4", "moderate_H4"])
    pw_best = (
        pw_links.sort_values(["Gene", "pwcoco_supported", "pwcoco_best_H4"], ascending=[True, False, False])
        .groupby("Gene", as_index=False)
        .first()[["Gene", "region_id", "pwcoco_best_h4_class", "pwcoco_best_H4", "pwcoco_supported", "priority_shared_signal"]]
        .rename(columns={"region_id": "pwcoco_best_region"})
    )
    master = master.merge(pw_best, on="Gene", how="left")

    # Merge SMR and cTWAS
    master = master.merge(
        smr_brainmeta.rename(
            columns={
                "best_context": "smr_bulk_best_context",
                "best_topSNP": "smr_bulk_best_topSNP",
                "best_p_SMR": "smr_bulk_best_p_SMR",
                "best_p_HEIDI": "smr_bulk_best_p_HEIDI",
                "best_fdr_SMR": "smr_bulk_best_fdr_SMR",
                "supported": "smr_bulk_supported",
            }
        ),
        on="Gene",
        how="left",
    )
    master = master.merge(
        smr_gtex.rename(
            columns={
                "best_context": "smr_gtex_best_tissue",
                "best_topSNP": "smr_gtex_best_topSNP",
                "best_p_SMR": "smr_gtex_best_p_SMR",
                "best_p_HEIDI": "smr_gtex_best_p_HEIDI",
                "best_fdr_SMR": "smr_gtex_best_fdr_SMR",
                "supported": "smr_gtex_supported",
            }
        ),
        on="Gene",
        how="left",
    )
    master = master.merge(
        smr_bryois.rename(
            columns={
                "shared_celltypes": "smr_bryois_shared_celltypes",
                "shared_n_celltypes": "smr_bryois_shared_n_celltypes",
                "shared_best_context": "smr_bryois_best_celltype",
                "shared_best_p_SMR_AD": "smr_bryois_best_p_SMR_AD",
                "shared_best_p_SMR_lipid8_F2": "smr_bryois_best_p_SMR_lipid8_F2",
                "shared_best_joint_score": "smr_bryois_best_joint_p_SMR",
                "shared_supported": "smr_bryois_supported",
                "any_best_context": "smr_bryois_any_best_celltype",
                "any_best_topSNP": "smr_bryois_any_best_topSNP",
                "any_best_p_SMR": "smr_bryois_any_best_p_SMR",
                "any_best_p_HEIDI": "smr_bryois_any_best_p_HEIDI",
                "any_best_fdr_SMR": "smr_bryois_any_best_fdr_SMR",
                "any_supported": "smr_bryois_any_supported",
            }
        ),
        on="Gene",
        how="left",
    )
    master = master.merge(ctwas_best, on="Gene", how="left")

    bool_cols = [
        "fuma_supported",
        "fuma_posMap_supported",
        "fuma_eqtlMap_supported",
        "fuma_ciMap_supported",
        "fuma_brain_eqtl_supported",
        "magma_gene_supported",
        "coloc_supported",
        "pwcoco_supported",
        "priority_shared_signal",
        "smr_bulk_supported",
        "smr_gtex_supported",
        "smr_bryois_supported",
        "ctwas_supported",
    ]
    for col in bool_cols:
        if col in master.columns:
            master[col] = master[col].fillna(False).astype(bool)

    master["evidence_locus"] = master["fuma_supported"] | master["magma_gene_supported"]
    master["evidence_smr_bulk"] = master["smr_bulk_supported"]
    master["evidence_smr_gtex"] = master["smr_gtex_supported"]
    master["evidence_smr_bryois"] = master["smr_bryois_supported"]
    master["evidence_ctwas"] = master["ctwas_supported"]
    master["evidence_coloc"] = master["coloc_supported"]
    master["evidence_pwcoco"] = master["pwcoco_supported"]

    evidence_cols = [
        "evidence_locus",
        "evidence_coloc",
        "evidence_pwcoco",
        "evidence_smr_bulk",
        "evidence_smr_gtex",
        "evidence_smr_bryois",
        "evidence_ctwas",
    ]
    master["n_evidence_layers"] = master[evidence_cols].sum(axis=1)

    def assign_tier(row: pd.Series) -> str:
        if row["evidence_ctwas"] and (row["evidence_smr_bulk"] or row["evidence_smr_gtex"] or row["evidence_smr_bryois"]) and row["evidence_locus"] and (row["evidence_coloc"] or row["evidence_pwcoco"]):
            return "A_high_convergent"
        if row["n_evidence_layers"] >= 4:
            return "B_multi_source"
        if row["n_evidence_layers"] >= 3:
            return "C_supportive"
        if row["n_evidence_layers"] >= 2:
            return "D_two_layer"
        return "E_single_layer"

    master["priority_tier"] = master.apply(assign_tier, axis=1)
    master["priority_score"] = (
        master["evidence_ctwas"].astype(int) * 3
        + (master["evidence_smr_bulk"] | master["evidence_smr_gtex"] | master["evidence_smr_bryois"]).astype(int) * 2
        + master["evidence_coloc"].astype(int)
        + master["evidence_pwcoco"].astype(int)
        + master["evidence_locus"].astype(int)
    )

    bryois_sort_col = "smr_bryois_best_joint_p_SMR" if "smr_bryois_best_joint_p_SMR" in master.columns else "smr_bryois_any_best_p_SMR"
    master = master.sort_values(
        ["priority_score", "n_evidence_layers", "max_pip", "smr_bulk_best_p_SMR", "smr_gtex_best_p_SMR", bryois_sort_col, "Gene"],
        ascending=[False, False, False, True, True, True, True],
    )

    shortlist = master[master["priority_tier"].isin(["A_high_convergent", "B_multi_source", "C_supportive"])].copy()

    master_path = OUTDIR / "lipid8_F2_AD_cross_evidence_master.tsv"
    shortlist_path = OUTDIR / "lipid8_F2_AD_cross_evidence_shortlist.tsv"
    master.to_csv(master_path, sep="\t", index=False)
    shortlist.to_csv(shortlist_path, sep="\t", index=False)

    counts = (
        master.groupby("priority_tier")
        .size()
        .reset_index(name="n_genes")
        .sort_values("priority_tier")
    )
    counts.to_csv(OUTDIR / "lipid8_F2_AD_cross_evidence_tier_counts.tsv", sep="\t", index=False)

    try:
        with pd.ExcelWriter(OUTDIR / "lipid8_F2_AD_cross_evidence_tables.xlsx") as writer:
            master.to_excel(writer, sheet_name="master", index=False)
            shortlist.to_excel(writer, sheet_name="shortlist", index=False)
            counts.to_excel(writer, sheet_name="tier_counts", index=False)
    except Exception as exc:
        (OUTDIR / "xlsx_write_error.txt").write_text(str(exc), encoding="utf-8")


if __name__ == "__main__":
    main()
