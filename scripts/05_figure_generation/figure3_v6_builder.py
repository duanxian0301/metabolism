import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

workbook = r"D:\codex\GenomicSEM\metabolic\postgwas_ad_pdlbd\results\22_supplement_tables_lipid8_F2_AD\metabolic_factor_triplet_supplementary_tables_v12_submission_clean.xlsx"
out_png = r"figures/figure3_shared_loci_genes_v6_final_polished.png"
out_pdf = r"figures/figure3_shared_loci_genes_v6_final_polished.pdf"
out_svg = r"figures/figure3_shared_loci_genes_v6_final_polished.svg"

pair_order = ["lipid8_F2_AD", "nonlipid8_F1_PD", "lipid8_F1_PD"]
pair_labels = {
    "lipid8_F2_AD": "HDL-core - AD",
    "nonlipid8_F1_PD": "Ketone-body - PD",
    "lipid8_F1_PD": "TG/VLDL - PD",
}
pair_colors = {
    "lipid8_F2_AD": "#5B86B2",
    "nonlipid8_F1_PD": "#64A99A",
    "lipid8_F1_PD": "#C77C5F",
}
pair_conjfdr_col = {
    "lipid8_F2_AD": "conjfdr_lipid8_F2_AD",
    "nonlipid8_F1_PD": "conjfdr_nonlipid8_F1_PD",
    "lipid8_F1_PD": "conjfdr_lipid8_F1_PD",
}

summary = pd.read_excel(workbook, sheet_name="S26_PleioFDR_summary", header=6)
loci = pd.read_excel(workbook, sheet_name="S27_PleioFDR_loci", header=6)
snp = pd.read_excel(workbook, sheet_name="S30_SNP_evidence", header=6)
genes = pd.read_excel(workbook, sheet_name="S35_Core_genes", header=6)

for df in [summary, loci, snp, genes]:
    if "pair" in df.columns:
        df["pair"] = df["pair"].astype(str)

def as_bool(x):
    if isinstance(x, bool):
        return x
    if pd.isna(x):
        return False
    if isinstance(x, (int, float)):
        return x != 0
    return str(x).strip().lower() in {"true", "1", "yes", "y"}

def neglog10(x, cap=30):
    if pd.isna(x):
        return np.nan
    x = float(x)
    if x <= 0:
        return cap
    return min(-math.log10(x), cap)

summary = summary.set_index("pair").loc[pair_order].reset_index()
summary["label"] = summary["pair"].map(pair_labels)
summary["color"] = summary["pair"].map(pair_colors)

# Approximate hg19 chromosome lengths for chromosomal strip plotting.
chr_lengths = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
    6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
    11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
    16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
    21: 48129895, 22: 51304566,
}
chr_starts = {}
offset = 0
for chrom in range(1, 23):
    chr_starts[chrom] = offset
    offset += chr_lengths[chrom]
chr_centers = {c: chr_starts[c] + chr_lengths[c] / 2 for c in range(1, 23)}

loci = loci[loci["pair"].isin(pair_order)].copy()
loci["chrnum"] = pd.to_numeric(loci["chrnum"], errors="coerce").astype("Int64")
loci["chrpos"] = pd.to_numeric(loci["chrpos"], errors="coerce")
loci = loci[loci["chrnum"].between(1, 22) & loci["chrpos"].notna()].copy()
loci["genome_pos"] = loci.apply(lambda r: chr_starts[int(r["chrnum"])] + float(r["chrpos"]), axis=1)
loci["chr_equal_pos"] = loci.apply(lambda r: int(r["chrnum"]) - 1 + float(r["chrpos"]) / chr_lengths[int(r["chrnum"])], axis=1)
loci["neglog_conjfdr"] = np.nan
for pair, col in pair_conjfdr_col.items():
    idx = loci["pair"].eq(pair)
    loci.loc[idx, "neglog_conjfdr"] = loci.loc[idx, col].map(lambda v: neglog10(v, 12))

snp = snp[snp["pair"].isin(pair_order)].copy()
for c in ["support_count_core", "coloc_PP.H4", "pwcoco_best_H4"]:
    snp[c] = pd.to_numeric(snp[c], errors="coerce")
snp["coloc_supported_b"] = snp["coloc_supported"].map(as_bool)
snp["pwcoco_supported_b"] = snp["pwcoco_supported"].map(as_bool)
snp["fuma_supported_b"] = snp["fuma_any_map_supported"].map(as_bool)
snp["sort_h4"] = snp[["coloc_PP.H4", "pwcoco_best_H4"]].max(axis=1)
top_loci = (
    snp.sort_values(["pair", "support_count_core", "sort_h4"], ascending=[True, False, False])
       .groupby("pair", sort=False)
       .head(4)
       .copy()
)
top_loci["row_label"] = top_loci.apply(
    lambda r: f"{str(r.get('nearestGene', 'NA'))[:12]} ({str(r.get('lead_snp', ''))})", axis=1
)
top_loci["pair"] = pd.Categorical(top_loci["pair"], pair_order, ordered=True)
top_loci = top_loci.sort_values(["pair", "support_count_core", "sort_h4"], ascending=[True, False, False]).reset_index(drop=True)
top_loci["shown_evidence_n"] = (
    1
    + top_loci["coloc_supported_b"].map(as_bool).astype(int)
    + top_loci["pwcoco_supported_b"].map(as_bool).astype(int)
    + top_loci["fuma_supported_b"].map(as_bool).astype(int)
)

genes = genes[genes["pair"].isin(pair_order)].copy()
for c in ["priority_score", "n_evidence_layers"]:
    genes[c] = pd.to_numeric(genes[c], errors="coerce")
smr_source_cols = ["evidence_smr_gtex", "evidence_smr_bulk", "evidence_smr_bryois"]
genes["evidence_smr_any"] = genes[smr_source_cols].apply(lambda row: any(as_bool(v) for v in row), axis=1)
evidence_cols = [
    ("Locus", "evidence_locus"),
    ("Shared\nregion", "evidence_shared_signal_region"),
    ("FUMA", "fuma_supported"),
    ("SMR", "evidence_smr_any"),
    ("cTWAS", "evidence_ctwas"),
]
genes["evidence_shared_signal_region"] = genes[["evidence_coloc", "evidence_pwcoco"]].apply(lambda row: any(as_bool(v) for v in row), axis=1)
for _, c in evidence_cols:
    genes[c + "_b"] = genes[c].map(as_bool)
top_genes = (
    genes.sort_values(["pair", "priority_score", "n_evidence_layers", "Gene"], ascending=[True, False, False, True])
         .groupby("pair", sort=False)
         .head(5)
         .copy()
)
top_genes["pair"] = pd.Categorical(top_genes["pair"], pair_order, ordered=True)
top_genes = top_genes.sort_values(["pair", "priority_score", "n_evidence_layers", "Gene"], ascending=[True, False, False, True]).reset_index(drop=True)
top_genes["shown_evidence_n"] = 0
for _, col in evidence_cols:
    top_genes["shown_evidence_n"] += top_genes[col + "_b"].map(as_bool).astype(int)

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.68,
})

fig = plt.figure(figsize=(11.8, 6.75), dpi=300)
gs = GridSpec(2, 6, figure=fig, height_ratios=[0.88, 1.18], width_ratios=[1, 1, 1, 1, 1, 1], wspace=0.90, hspace=0.40)

# A. Shared locus counts.
axA = fig.add_subplot(gs[0, 0:2])
y = np.arange(len(summary))
axA.barh(y, summary["unique_locusnum"], color=summary["color"], height=0.44)
for i, v in enumerate(summary["unique_locusnum"]):
    axA.text(v + 4, i, f"{int(v)}", va="center", ha="left", fontsize=7.0)
axA.set_yticks(y)
axA.set_yticklabels(summary["label"])
axA.invert_yaxis()
axA.set_xlabel("Shared loci, conjFDR < 0.05")
axA.set_xlim(0, max(summary["unique_locusnum"]) * 1.23)
axA.set_title("Shared-locus discovery", loc="left", fontweight="bold")
axA.tick_params(axis="y", length=0)
for s in ["top", "right"]:
    axA.spines[s].set_visible(False)
axA.text(-0.22, 1.13, "A", transform=axA.transAxes, fontsize=13.5, fontweight="bold", va="top")

# B. Chromosomal distribution.
axB = fig.add_subplot(gs[0, 2:6])
for i, pair in enumerate(pair_order):
    d = loci[loci["pair"].eq(pair)]
    sizes = 6 + 4.7 * d["neglog_conjfdr"].fillna(1).clip(0, 12)
    axB.scatter(d["chr_equal_pos"], np.full(len(d), i), s=sizes, color=pair_colors[pair], alpha=0.72, edgecolor="white", linewidth=0.10)
for c in range(23):
    axB.axvline(c, color="#E9E9E9", lw=0.36, zorder=0)
axB.set_yticks(range(len(pair_order)))
axB.set_yticklabels([pair_labels[p] for p in pair_order])
axB.invert_yaxis()
axB.set_xlim(0, 22)
axB.set_xticks(np.arange(22) + 0.5)
axB.set_xticklabels([str(c) for c in range(1, 23)], fontsize=5.8)
axB.set_xlabel("Chromosome")
axB.set_title("Chromosomal map of shared loci", loc="left", fontweight="bold")
axB.tick_params(axis="y", length=0)
for s in ["top", "right", "left"]:
    axB.spines[s].set_visible(False)
axB.text(-0.08, 1.13, "B", transform=axB.transAxes, fontsize=13.5, fontweight="bold", va="top")

# C. Top locus-level evidence.
axC = fig.add_subplot(gs[1, 0:3])
c_cols = ["conjFDR", "coloc", "PWCoCo", "FUMA"]
c_x = np.arange(len(c_cols))
c_y = np.arange(len(top_loci))
for yi, r in top_loci.iterrows():
    pair = str(r["pair"])
    vals = [
        1.0,
        float(r["coloc_PP.H4"]) if as_bool(r["coloc_supported_b"]) and pd.notna(r["coloc_PP.H4"]) else np.nan,
        float(r["pwcoco_best_H4"]) if as_bool(r["pwcoco_supported_b"]) and pd.notna(r["pwcoco_best_H4"]) else np.nan,
        1.0 if as_bool(r["fuma_supported_b"]) else np.nan,
    ]
    for xi, v in enumerate(vals):
        if pd.notna(v):
            size = 26 + 95 * min(max(v, 0), 1)
            axC.scatter(xi, yi, s=size, color=pair_colors[pair], edgecolor="#2B2B2B", lw=0.30, alpha=0.92)
axC.set_xticks(c_x)
axC.set_xticklabels(c_cols, rotation=24, ha="right")
axC.set_yticks(c_y)
axC.set_yticklabels(top_loci["row_label"])
axC.invert_yaxis()
axC.set_xlim(-0.55, len(c_cols) - 0.45)
axC.set_title("Locus-level shared-signal evidence", loc="left", fontweight="bold")
axC.tick_params(axis="both", length=0)
for yi in range(len(top_loci)):
    axC.axhline(yi + 0.5, color="#EFEFEF", lw=0.32, zorder=0)
for xi in c_x:
    axC.axvline(xi, color="#F3F3F3", lw=0.35, zorder=0)
for s in axC.spines.values():
    s.set_visible(False)
axC.text(-0.17, 1.08, "C", transform=axC.transAxes, fontsize=13.5, fontweight="bold", va="top")

# D. Core-gene evidence matrix.
axD = fig.add_subplot(gs[1, 3:6])
d_cols = [x[0] for x in evidence_cols]
d_x = np.arange(len(d_cols))
d_y = np.arange(len(top_genes))
for yi, r in top_genes.iterrows():
    pair = str(r["pair"])
    for xi, (_, col) in enumerate(evidence_cols):
        if as_bool(r[col + "_b"]):
            axD.scatter(xi, yi, s=58, color=pair_colors[pair], edgecolor="#2B2B2B", lw=0.30)
axD.set_xticks(d_x)
axD.set_xticklabels(d_cols, rotation=24, ha="right")
axD.set_yticks(d_y)
axD.set_yticklabels(top_genes["Gene"])
axD.invert_yaxis()
axD.set_xlim(-0.55, len(d_cols) - 0.45)
axD.set_title("Prioritized genes linked to shared-signal loci", loc="left", fontweight="bold")
axD.tick_params(axis="both", length=0)
for yi in range(len(top_genes)):
    axD.axhline(yi + 0.5, color="#EFEFEF", lw=0.32, zorder=0)
for xi in d_x:
    axD.axvline(xi, color="#F3F3F3", lw=0.35, zorder=0)
for s in axD.spines.values():
    s.set_visible(False)
axD.text(-0.15, 1.08, "D", transform=axD.transAxes, fontsize=13.5, fontweight="bold", va="top")

legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=pair_colors[p], markeredgecolor="#2B2B2B", markersize=6.0, label=pair_labels[p])
    for p in pair_order
]
fig.legend(handles=legend_handles, frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.52, 0.006), fontsize=6.9, handletextpad=0.45, columnspacing=1.35)

plt.savefig(out_png, bbox_inches="tight", dpi=300)
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
print(out_png)
print(out_pdf)
print(out_svg)
