from pathlib import Path
import ast
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(r"D:\codex\GenomicSEM\metabolic")
FIG = ROOT / "figures"
SUPP = ROOT / "postgwas_ad_pdlbd/results/22_supplement_tables_lipid8_F2_AD/metabolic_factor_triplet_supplementary_tables_v12_submission_clean.xlsx"
MOD = ROOT / "manuscript/target_validation_modules"
EXT = MOD / "routeA_external_validation/routeA_external_disease_state_validation_summary.csv"

OUT_PNG = FIG / "figure5_translational_prioritization_v2_smr_union.png"
OUT_PDF = FIG / "figure5_translational_prioritization_v2_smr_union.pdf"
OUT_SVG = FIG / "figure5_translational_prioritization_v2_smr_union.svg"

PAIR_ORDER = ["lipid8_F2_AD", "nonlipid8_F1_PD", "lipid8_F1_PD"]
PAIR_LABELS = {
    "lipid8_F2_AD": "HDL-core / AD",
    "nonlipid8_F1_PD": "Ketone-core / PD",
    "lipid8_F1_PD": "TG/VLDL-core / PD",
}
PAIR_COLORS = {
    "lipid8_F2_AD": "#4E7FA8",
    "nonlipid8_F1_PD": "#4C9A86",
    "lipid8_F1_PD": "#C66F55",
}
TIER_COLORS = {
    "A_binding_or_drug_traceable": "#9C2F2F",
    "B_modality_supported": "#D08A3C",
    "C_structure_or_network_supported": "#7A90A4",
}
EVIDENCE_COLORS = {
    "genetic": "#4E7FA8",
    "regulatory": "#6C8F63",
    "cell": "#C66F55",
    "external": "#8B6BB8",
    "translation": "#D08A3C",
}


def as_bool(x):
    if isinstance(x, bool):
        return x
    if pd.isna(x):
        return False
    if isinstance(x, (int, float)):
        return bool(x)
    return str(x).strip().lower() in {"true", "1", "yes", "y"}


def parse_pd_score(txt):
    if pd.isna(txt) or not str(txt).strip():
        return {}
    out = {}
    for part in str(txt).split(";"):
        if "=" in part:
            k, v = part.strip().split("=", 1)
            try:
                out[k.strip()] = float(v)
            except ValueError:
                out[k.strip()] = v.strip()
    return out


def parse_ad_effect(txt):
    if pd.isna(txt) or not str(txt).strip():
        return np.nan, np.nan
    try:
        d = ast.literal_eval(str(txt))
        return float(d.get("Effect", np.nan)), float(d.get("p-value FDR", np.nan))
    except Exception:
        return np.nan, np.nan


def panel_label(ax, label, x=-0.10, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13.5, fontweight="bold", va="top")


core = pd.read_csv(MOD / "knk_prioritized_core_genes.csv")
target = pd.read_csv(MOD / "module1_3_4_consolidated_target_summary.csv")
tract = pd.read_csv(MOD / "module4_binding_modality_summary.csv")
brain = pd.read_csv(MOD / "module3_brain_expression_disease_support_summary.csv")
ext = pd.read_csv(EXT)
s35 = pd.read_excel(SUPP, sheet_name="S35_Core_genes", header=6)

target_genes = target["gene"].dropna().astype(str).tolist()
priority_gene_order = (
    tract.set_index("gene")
    .reindex(target_genes)
    .sort_values(["binding_modality_score", "gene"], ascending=[False, True])
    .index.tolist()
)

# Evidence convergence matrix.
s35 = s35[s35["Gene"].isin(target_genes)].copy()
evidence_layers = [
    ("Shared locus", "genetic"),
    ("coloc/PWCoCo", "genetic"),
    ("SMR", "regulatory"),
    ("cTWAS", "regulatory"),
    ("scPagwas cell", "cell"),
    ("KNK perturb.", "cell"),
    ("Disease-state", "external"),
    ("Brain protein", "external"),
    ("Tractability", "translation"),
]
matrix = pd.DataFrame(0, index=priority_gene_order, columns=[x[0] for x in evidence_layers])
for gene in priority_gene_order:
    d = s35[s35["Gene"].eq(gene)]
    if len(d):
        matrix.loc[gene, "Shared locus"] = int(d["evidence_locus"].map(as_bool).any())
        matrix.loc[gene, "coloc/PWCoCo"] = int(d["evidence_coloc"].map(as_bool).any() or d["evidence_pwcoco"].map(as_bool).any())
        matrix.loc[gene, "SMR"] = int(
            d["evidence_smr_gtex"].map(as_bool).any()
            or d["evidence_smr_bulk"].map(as_bool).any()
            or d["evidence_smr_bryois"].map(as_bool).any()
        )
        matrix.loc[gene, "cTWAS"] = int(d["evidence_ctwas"].map(as_bool).any())
matrix.loc[matrix.index.isin(core["gene"].astype(str)), "KNK perturb."] = 1
matrix.loc[matrix.index.isin(core.loc[core["remarks"].astype(str).str.contains("scPagwas", case=False, na=False), "gene"].astype(str)), "scPagwas cell"] = 1
matched_genes = set(ext.loc[ext["matched"].astype(str).str.lower().eq("yes"), "gene"].astype(str))
matrix.loc[matrix.index.isin(matched_genes), "Disease-state"] = 1
brain_prot = set(brain.loc[brain["has_brain_protein_evidence"].map(as_bool), "gene"].astype(str))
matrix.loc[matrix.index.isin(brain_prot), "Brain protein"] = 1
tractable = set(tract.loc[tract["binding_modality_tier"].astype(str).str.startswith(("A_", "B_")), "gene"].astype(str))
matrix.loc[matrix.index.isin(tractable), "Tractability"] = 1

gene_pairs = (
    core.groupby("gene")["pair"]
    .apply(lambda x: sorted(set(map(str, x))))
    .reindex(priority_gene_order)
    .to_dict()
)

# Disease-state summaries.
pd_rows = []
for gene in priority_gene_order:
    row = ext[(ext["gene"].eq(gene)) & (ext["disease"].eq("PD"))]
    info = parse_pd_score(row["effect_or_score"].iloc[0]) if len(row) else {}
    pd_rows.append({
        "gene": gene,
        "D_overall": info.get("D.overall", np.nan),
        "rank": info.get("Rank.overall", np.nan),
        "directionality": info.get("directionality_score", np.nan),
        "midbrain": str(info.get("expressed_Midbrain", "NO")).upper() == "YES",
    })
pd_support = pd.DataFrame(pd_rows)
ad_rows = []
for gene in priority_gene_order:
    rows = ext[(ext["gene"].eq(gene)) & (ext["disease"].eq("AD")) & (ext["matched"].astype(str).str.lower().eq("yes"))]
    effs, fdrs = [], []
    for _, r in rows.iterrows():
        e, f = parse_ad_effect(r["effect_or_score"])
        if not np.isnan(e):
            effs.append(e)
            fdrs.append(f)
    ad_rows.append({"gene": gene, "matched": len(effs) > 0, "effect": np.nanmean(effs) if effs else np.nan, "fdr": np.nanmin(fdrs) if fdrs else np.nan})
ad_support = pd.DataFrame(ad_rows)

# Axis-gene map.
axis_edges = []
for gene in priority_gene_order:
    pairs = gene_pairs.get(gene, [])
    for pair in pairs:
        axis_edges.append((pair, gene))

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7.4,
    "axes.titlesize": 8.7,
    "axes.labelsize": 7.6,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 7.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.65,
})

fig = plt.figure(figsize=(12.8, 7.7), dpi=300)
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.48, 1.0], height_ratios=[1.15, 1.0], wspace=0.36, hspace=0.42)

# A. Evidence matrix.
axA = fig.add_subplot(gs[0, 0])
axA.set_title("Convergent evidence supporting prioritized genes", loc="left", fontweight="bold")
vals = matrix.to_numpy()
for i, gene in enumerate(matrix.index):
    for j, (layer, cat) in enumerate(evidence_layers):
        color = EVIDENCE_COLORS[cat] if vals[i, j] else "#EFEFEF"
        axA.scatter(j, i, s=62 if vals[i, j] else 22, marker="s", color=color, edgecolor="white", linewidth=0.45)
for j in range(len(evidence_layers)):
    axA.axvline(j, color="#F4F4F4", lw=0.35, zorder=0)
axA.set_xticks(range(len(evidence_layers)))
axA.set_xticklabels([x[0] for x in evidence_layers], rotation=35, ha="right")
axA.set_yticks(range(len(matrix.index)))
axA.set_yticklabels(matrix.index)
axA.invert_yaxis()
axA.tick_params(length=0)
for s in axA.spines.values():
    s.set_visible(False)
legend_e = [Patch(facecolor=v, edgecolor="none", label=k.capitalize()) for k, v in EVIDENCE_COLORS.items()]
axA.legend(handles=legend_e, frameon=False, ncol=5, fontsize=6.5, loc="upper left", bbox_to_anchor=(0.00, -0.23))
panel_label(axA, "A")

# B. External disease-state support.
axB = fig.add_subplot(gs[0, 1])
axB.set_title("Supportive public disease-state summaries", loc="left", fontweight="bold")
y = np.arange(len(priority_gene_order))
pd_support = pd_support.set_index("gene").reindex(priority_gene_order)
ad_support = ad_support.set_index("gene").reindex(priority_gene_order)
size = 35 + 360 * (pd_support["D_overall"].fillna(0) / max(0.001, pd_support["D_overall"].max()))
colors = np.where(pd_support["midbrain"], "#6F3E8E", "#B7A6C7")
axB.scatter(pd_support["D_overall"], y, s=size, color=colors, edgecolor="white", linewidth=0.65, label="PD transcriptomic meta-signature")
ad_y = [i for i, g in enumerate(priority_gene_order) if bool(ad_support.loc[g, "matched"])]
if ad_y:
    ad_x = [0.112] * len(ad_y)
    axB.scatter(ad_x, ad_y, s=68, marker="D", color="#4E7FA8", edgecolor="white", linewidth=0.65, label="AD CSF proteomics present")
axB.set_yticks(y)
axB.set_yticklabels(priority_gene_order)
axB.invert_yaxis()
axB.set_xlabel("PD D.overall score")
axB.set_xlim(-0.004, 0.122)
axB.grid(axis="x", color="#E8E8E8", lw=0.45)
axB.set_axisbelow(True)
for s in ["top", "right"]:
    axB.spines[s].set_visible(False)
axB.legend(frameon=False, fontsize=6.3, loc="lower right")
axB.text(0.02, 0.02, "Purple fill: midbrain expressed", transform=axB.transAxes, fontsize=6.1, color="#666666")
panel_label(axB, "B")

# C. Tractability and modality.
axC = fig.add_subplot(gs[1, 0])
tr = tract.set_index("gene").reindex(priority_gene_order).reset_index()
tr = tr.sort_values("binding_modality_score", ascending=True)
y2 = np.arange(len(tr))
bar_cols = [TIER_COLORS.get(t, "#AAAAAA") for t in tr["binding_modality_tier"]]
axC.barh(y2, tr["binding_modality_score"], color=bar_cols, height=0.62)
axC.set_yticks(y2)
axC.set_yticklabels(tr["gene"])
axC.set_xlabel("Binding/modality support score")
axC.set_title("Target tractability and modality annotation", loc="left", fontweight="bold")
for yi, (_, r) in enumerate(tr.iterrows()):
    tags = []
    if as_bool(r.get("direct_approved_drug_gene_evidence")):
        tags.append("approved drug")
    if as_bool(r.get("small_molecule_binding_or_pocket_support")):
        tags.append("small molecule")
    if as_bool(r.get("chemical_probe_available")):
        tags.append("probe")
    if as_bool(r.get("antibody_or_surface_modality_support")):
        tags.append("surface/antibody")
    axC.text(r["binding_modality_score"] + 0.13, yi, ", ".join(tags) if tags else "network/structure", va="center", fontsize=6.4, color="#555555")
axC.set_xlim(0, 9.7)
axC.grid(axis="x", color="#E8E8E8", lw=0.45)
axC.set_axisbelow(True)
for s in ["top", "right"]:
    axC.spines[s].set_visible(False)
tier_handles = [Patch(facecolor=v, edgecolor="none", label=k.split("_", 1)[0]) for k, v in TIER_COLORS.items()]
axC.legend(handles=tier_handles, title="Tier", frameon=False, fontsize=6.4, title_fontsize=6.6, loc="lower right")
panel_label(axC, "C")

# D. Axis-to-target convergence map.
axD = fig.add_subplot(gs[1, 1])
axD.set_title("Metabolic axes converging on target candidates", loc="left", fontweight="bold")
axis_y = {p: len(priority_gene_order) - 1 - i * (len(priority_gene_order) - 1) / 2 for i, p in enumerate(PAIR_ORDER)}
gene_y = {g: i for i, g in enumerate(priority_gene_order)}
for pair in PAIR_ORDER:
    axD.scatter(0, axis_y[pair], s=120, color=PAIR_COLORS[pair], edgecolor="white", linewidth=0.8, zorder=3)
    axD.text(-0.04, axis_y[pair], PAIR_LABELS[pair], va="center", ha="right", fontsize=7.0, color=PAIR_COLORS[pair], fontweight="bold")
score_map = tract.set_index("gene")["binding_modality_score"].to_dict()
for gene in priority_gene_order:
    sc = 38 + 13 * score_map.get(gene, 1)
    tier = tract.set_index("gene").loc[gene, "binding_modality_tier"]
    axD.scatter(1, gene_y[gene], s=sc, color=TIER_COLORS.get(tier, "#AAAAAA"), edgecolor="white", linewidth=0.7, zorder=3)
    axD.text(1.055, gene_y[gene], gene, va="center", ha="left", fontsize=7.0)
for pair, gene in axis_edges:
    axD.plot([0.04, 0.96], [axis_y[pair], gene_y[gene]], color=PAIR_COLORS[pair], alpha=0.48, lw=1.1)
axD.set_xlim(-0.36, 1.38)
axD.set_ylim(-0.8, len(priority_gene_order) - 0.2)
axD.axis("off")
panel_label(axD, "D", x=-0.05)

fig.subplots_adjust(left=0.070, right=0.985, top=0.945, bottom=0.110)
fig.savefig(OUT_PNG, dpi=450)
fig.savefig(OUT_PDF)
fig.savefig(OUT_SVG)
print(OUT_PNG)
