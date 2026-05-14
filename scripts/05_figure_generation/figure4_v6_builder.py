from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe


ROOT = Path(r"D:\codex\GenomicSEM\metabolic")
FIG = ROOT / "figures"
UMAP_DIR = FIG / "figure4_scpagwas_umap_data"
WORKBOOK = ROOT / "postgwas_ad_pdlbd/results/22_supplement_tables_lipid8_F2_AD/metabolic_factor_triplet_supplementary_tables_v12_submission_clean.xlsx"

OUT_PNG = FIG / "figure4_singlecell_pathway_v6_clean_labels.png"
OUT_PDF = FIG / "figure4_singlecell_pathway_v6_clean_labels.pdf"
OUT_SVG = FIG / "figure4_singlecell_pathway_v6_clean_labels.svg"

PAIR_ORDER = ["lipid8_F2_AD", "nonlipid8_F1_PD", "lipid8_F1_PD"]
PAIR_LABELS = {
    "lipid8_F2_AD": "HDL-core / AD",
    "nonlipid8_F1_PD": "Ketone-core / PD",
    "lipid8_F1_PD": "TG/VLDL-core / PD",
}
KEY_TO_PAIR = {
    "HDL_core_AD": "lipid8_F2_AD",
    "Ketone_core_PD": "nonlipid8_F1_PD",
    "TG_VLDL_core_PD": "lipid8_F1_PD",
}
KEY_TITLES = {
    "HDL_core_AD": "HDL-core / AD TRS score",
    "Ketone_core_PD": "Ketone-core / PD TRS score",
    "TG_VLDL_core_PD": "TG/VLDL-core / PD TRS score",
}
PAIR_COLORS = {
    "lipid8_F2_AD": "#4E7FA8",
    "nonlipid8_F1_PD": "#4C9A86",
    "lipid8_F1_PD": "#C66F55",
}
METHOD_COLORS = {"scPagwas": "#3F8F91", "KNK": "#D28A34"}
TRS_CMAP = LinearSegmentedColormap.from_list("trs_gray_red", ["#D6D6D6", "#F2B8A5", "#D85A4A", "#6D0015"])

CELL_COLORS = {
    "Astrocyte": "#7A6BB7",
    "Endothelial": "#D67B83",
    "GABA neuron": "#E2A54E",
    "L2/3 IT": "#8A9E3D",
    "L2/3-6 IT": "#B6A038",
    "L5/6 NP": "#C18A3D",
    "L6 CT": "#9C7E53",
    "L6 IT": "#BA8EA8",
    "L6b Glu": "#B97B4A",
    "LAMP5 GABA": "#E0BD62",
    "Microglia": "#548BBE",
    "NK cell": "#7B7B7B",
    "OPC": "#65A963",
    "Oligodendrocyte": "#4E69B2",
    "PVALB GABA": "#E3A25B",
    "PVM": "#7AA6A1",
    "Pericyte": "#D35F49",
    "SMC": "#9E8E77",
    "SST GABA": "#D2A34B",
    "T cell": "#9B9B9B",
    "VIP GABA": "#D9B24D",
    "VLMC": "#B56E7E",
}

KEGG_NAMES = {
    "hsa03010": "Ribosome",
    "hsa03040": "Spliceosome",
    "hsa03420": "Nucleotide excision repair",
    "hsa04015": "Rap1 signaling",
    "hsa04022": "cGMP-PKG signaling",
    "hsa04062": "Chemokine signaling",
    "hsa04066": "HIF-1 signaling",
    "hsa04145": "Phagosome",
    "hsa04261": "Adrenergic signaling",
    "hsa04330": "Notch signaling",
    "hsa04514": "Cell adhesion molecules",
    "hsa04520": "Adherens junction",
    "hsa04611": "Platelet activation",
    "hsa04621": "NOD-like receptor signaling",
    "hsa04662": "B-cell receptor signaling",
    "hsa04723": "Endocannabinoid signaling",
    "hsa04724": "Glutamatergic synapse",
    "hsa04728": "Dopaminergic synapse",
    "hsa04912": "GnRH signaling",
    "hsa04925": "Aldosterone synthesis/secretion",
    "hsa04926": "Relaxin signaling",
    "hsa04932": "NAFLD",
    "hsa04935": "Growth hormone synthesis/secretion",
    "hsa04960": "Aldosterone-regulated sodium reabsorption",
    "hsa05131": "Shigellosis",
    "hsa05167": "KSHV infection",
    "hsa05202": "Transcriptional misregulation",
    "hsa05230": "Central carbon metabolism",
    "hsa05410": "Hypertrophic cardiomyopathy",
    "hsa05414": "Dilated cardiomyopathy",
}


def short_label(text, n=28):
    text = str(text)
    return text if len(text) <= n else text[: n - 1] + "..."


def clean_cell(s):
    return str(s).replace("_", " ")


def star(fdr):
    if pd.isna(fdr):
        return ""
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    if fdr < 0.05:
        return "*"
    return ""


def panel_label(ax, label, x=-0.12, y=1.055):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13.5, fontweight="bold", va="top", ha="left")


def load_umap(key):
    df = pd.read_csv(UMAP_DIR / f"{key}_umap_trs.csv")
    df["pair"] = KEY_TO_PAIR[key]
    return df


def load_celltype_fdr():
    paths = {
        "HDL_core_AD": r"D:\scPagwas\metabolic_scpagwas2\lipid8_F2_MSSM_AD\lipid8_F2_MSSM_AD_Merged_celltype_pvalue_withFDR.csv",
        "Ketone_core_PD": r"D:\scPagwas\metabolic_scpagwas2\nonlipid8_F1_MSSM_PD\nonlipid8_F1_MSSM_PD_Merged_celltype_pvalue_withFDR.csv",
        "TG_VLDL_core_PD": r"D:\scPagwas\metabolic_scpagwas2\lipid8_F1_MSSM_PD\lipid8_F1_MSSM_PD_Merged_celltype_pvalue_withFDR.csv",
    }
    out = []
    for key, p in paths.items():
        d = pd.read_csv(p)
        d["analysis_key"] = key
        d["pair"] = KEY_TO_PAIR[key]
        out.append(d)
    return pd.concat(out, ignore_index=True)


def regulatory_support():
    genes = pd.read_excel(WORKBOOK, sheet_name="S35_Core_genes", header=6)
    genes = genes[genes["pair"].isin(PAIR_ORDER)].copy()
    cols = {
        "evidence_smr_gtex": "GTEx SMR",
        "evidence_smr_bulk": "Bulk SMR",
        "evidence_smr_bryois": "Cell SMR",
        "evidence_ctwas": "cTWAS",
    }
    rows = []
    for pair in PAIR_ORDER:
        d = genes[genes["pair"].eq(pair)]
        for col, lab in cols.items():
            vals = d[col].astype(str).str.lower().isin(["true", "1", "yes"])
            rows.append({"pair": pair, "layer": lab, "n": int(vals.sum())})
    return pd.DataFrame(rows)


def pathway_name_from_scpagwas():
    mp = dict(KEGG_NAMES)
    files = [
        ROOT / "postgwas_ad_pdlbd/results/21_knk_scpagwas_overlap_lipid8_F2_AD/scpagwas_pericyte_top50.tsv",
        ROOT / "postgwas_ad_pdlbd/results/26_knk_scpagwas_overlap_nonlipid8_F1_PD/scpagwas_pericyte_top50.tsv",
        ROOT / "postgwas_ad_pdlbd/results/26_knk_scpagwas_overlap_nonlipid8_F1_PD/scpagwas_oligodendrocyte_precursor_cell_top50.tsv",
        ROOT / "postgwas_ad_pdlbd/results/26_knk_scpagwas_overlap_nonlipid8_F1_PD/scpagwas_VIP_GABAergic_cortical_interneuron_top50.tsv",
        ROOT / "postgwas_ad_pdlbd/results/27_knk_scpagwas_overlap_lipid8_F1_PD/scpagwas_pericyte_top50.tsv",
    ]
    for p in files:
        if p.exists():
            d = pd.read_csv(p, sep="\t")
            if "pathway_id" in d and "pathway_name" in d:
                mp.update(dict(zip(d["pathway_id"], d["pathway_name"])))
    return mp


def sc_rank_file(pair, cell):
    if pair == "lipid8_F2_AD":
        return ROOT / "postgwas_ad_pdlbd/results/21_knk_scpagwas_overlap_lipid8_F2_AD/scpagwas_pericyte_top50.tsv"
    if pair == "nonlipid8_F1_PD":
        return ROOT / f"postgwas_ad_pdlbd/results/26_knk_scpagwas_overlap_nonlipid8_F1_PD/scpagwas_{cell}_top50.tsv"
    return ROOT / "postgwas_ad_pdlbd/results/27_knk_scpagwas_overlap_lipid8_F1_PD/scpagwas_pericyte_top50.tsv"


def knk_file(pair, cell, gene):
    if pair == "lipid8_F2_AD":
        base = ROOT / "postgwas_ad_pdlbd/results/19_knk_lipid8_F2_AD_core4"
        prefix = "MSSM_AD"
    elif pair == "nonlipid8_F1_PD":
        base = ROOT / "postgwas_ad_pdlbd/results/24_knk_nonlipid8_F1_PD"
        prefix = "MSSM_PD"
    else:
        base = ROOT / "postgwas_ad_pdlbd/results/25_knk_lipid8_F1_PD"
        prefix = "MSSM_PD"
    return base / cell / f"{prefix}_{cell}_{gene}" / f"{prefix}_{cell}_{gene}_pathway_enrichment.csv"


def pathway_overlap_records():
    names = pathway_name_from_scpagwas()
    configs = [
        ("AD", "lipid8_F2_AD", ROOT / "postgwas_ad_pdlbd/results/21_knk_scpagwas_overlap_lipid8_F2_AD/knk_vs_scpagwas_pericyte_trs_overlap.tsv"),
        ("PD", "nonlipid8_F1_PD", ROOT / "postgwas_ad_pdlbd/results/26_knk_scpagwas_overlap_nonlipid8_F1_PD/knk_vs_scpagwas_nonlipid8_F1_PD_overlap.tsv"),
        ("PD", "lipid8_F1_PD", ROOT / "postgwas_ad_pdlbd/results/27_knk_scpagwas_overlap_lipid8_F1_PD/knk_vs_scpagwas_lipid8_F1_PD_overlap.tsv"),
    ]
    rows = []
    for disease, pair, path in configs:
        ov = pd.read_csv(path, sep="\t")
        id_col = [c for c in ov.columns if c.endswith("top50_overlap_ids")][0]
        ref_col = "scpagwas_reference_cell" if "scpagwas_reference_cell" in ov.columns else None
        for _, r in ov.iterrows():
            ids = [x for x in str(r.get(id_col, "")).split(";") if x and x != "nan"]
            if not ids:
                continue
            cell = str(r["cell_type"])
            ref_cell = str(r[ref_col]) if ref_col else "pericyte"
            gene = str(r["gene"])
            sc_path = sc_rank_file(pair, ref_cell)
            sc = pd.read_csv(sc_path, sep="\t") if sc_path.exists() else pd.DataFrame()
            if len(sc):
                sc["rank_sc"] = np.arange(1, len(sc) + 1)
            kfile = knk_file(pair, cell, gene)
            knk = pd.read_csv(kfile) if kfile.exists() else pd.DataFrame()
            if len(knk):
                knk = knk.sort_values(["p.adjust", "pvalue"]).reset_index(drop=True)
                knk["rank_knk"] = np.arange(1, len(knk) + 1)
            for pid in ids:
                sc_rank = float(sc.loc[sc["pathway_id"].eq(pid), "rank_sc"].min()) if len(sc) and pid in set(sc["pathway_id"]) else np.nan
                knk_rank = float(knk.loc[knk["ID"].eq(pid), "rank_knk"].min()) if len(knk) and pid in set(knk["ID"]) else np.nan
                rows.append({
                    "disease": disease,
                    "pair": pair,
                    "pathway_id": pid,
                    "pathway": names.get(pid, KEGG_NAMES.get(pid, pid)),
                    "sc_score": 51 - sc_rank if not np.isnan(sc_rank) else np.nan,
                    "knk_score": 51 - knk_rank if not np.isnan(knk_rank) else np.nan,
                    "gene": gene,
                })
    df = pd.DataFrame(rows)
    return df.groupby(["disease", "pair", "pathway_id", "pathway"], as_index=False).agg(
        sc_score=("sc_score", "max"),
        knk_score=("knk_score", "max"),
        genes=("gene", lambda x: ", ".join(sorted(set(map(str, x)))[:4])),
    )


def draw_umap_celltypes(ax, data):
    for ct, d in data.groupby("cell_type_short", sort=False):
        ax.scatter(d["umap_1"], d["umap_2"], s=1.25, color=CELL_COLORS.get(ct, "#999999"), alpha=0.78, linewidth=0)
    ax.set_title("MSSM cell-type UMAP", loc="left", fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    med = data.groupby("cell_type_short")[["umap_1", "umap_2"]].median()
    label_pos = {
        "L2/3-6 IT": (-2.0, 11.0),
        "L2/3 IT": (-9.8, 7.8),
        "L6 IT": (-4.6, 5.7),
        "Oligodendrocyte": (7.1, 3.2),
        "L5/6 NP": (-15.8, 3.9),
        "PVALB GABA": (-11.4, 1.7),
        "SST GABA": (-8.0, 1.0),
        "L6 CT": (-16.7, 0.0),
        "L6b Glu": (-16.8, -1.5),
        "VIP GABA": (-2.1, -1.7),
        "Endothelial": (1.4, -1.0),
        "GABA neuron": (-4.7, -3.3),
        "LAMP5 GABA": (-8.7, -4.4),
        "Pericyte": (0.8, -4.65),
        "SMC": (1.3, -6.0),
        "VLMC": (-1.9, -7.0),
        "PVM": (5.4, -6.7),
        "Microglia": (6.7, -7.9),
        "NK cell": (2.6, -9.3),
        "T cell": (4.2, -9.45),
        "OPC": (-11.1, -8.9),
        "Astrocyte": (-3.3, -13.9),
    }
    for ct, r in med.iterrows():
        tx, ty = label_pos.get(ct, (r["umap_1"], r["umap_2"]))
        moved = abs(tx - r["umap_1"]) + abs(ty - r["umap_2"]) > 0.65
        ax.annotate(
            ct,
            xy=(r["umap_1"], r["umap_2"]),
            xytext=(tx, ty),
            fontsize=4.25,
            ha="center",
            va="center",
            color="#202020",
            fontweight="bold" if ct == "Pericyte" else "normal",
            arrowprops=dict(arrowstyle="-", lw=0.28, color="#777777", alpha=0.62, shrinkA=1, shrinkB=1) if moved else None,
            path_effects=[pe.withStroke(linewidth=1.05, foreground="white", alpha=0.92)],
        )


def draw_umap_trs(ax, data, title):
    score = data["scPagwas.TRS.Score"].to_numpy()
    lo, hi = np.nanpercentile(score, [3, 99])
    vals = (np.clip(score, lo, hi) - lo) / max(1e-8, hi - lo)
    sc = ax.scatter(data["umap_1"], data["umap_2"], c=vals, s=1.25, cmap=TRS_CMAP, alpha=0.86, linewidth=0, vmin=0, vmax=1)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return sc


def plot_overlap(ax, overlap, disease, title, nmax):
    d = overlap[overlap["disease"].eq(disease)].copy()
    d = d.sort_values(["sc_score", "knk_score"], ascending=False).head(nmax)
    d["label"] = d["pathway"].map(lambda x: short_label(x, 24))
    d = d.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(d))
    ax.barh(y + 0.18, d["sc_score"], height=0.30, color=METHOD_COLORS["scPagwas"], label="scPagwas")
    ax.barh(y - 0.18, d["knk_score"], height=0.30, color=METHOD_COLORS["KNK"], label="KNK")
    for yi, (_, r) in enumerate(d.iterrows()):
        ax.scatter(0.8, yi, s=20, color=PAIR_COLORS[r["pair"]], edgecolor="white", linewidth=0.35, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(d["label"])
    ax.set_xlim(-0.2, 51)
    ax.set_xlabel("Top-50 rank score")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(axis="x", color="#E8E8E8", lw=0.45)
    ax.set_axisbelow(True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)


mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7.3,
    "axes.titlesize": 8.4,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.8,
    "ytick.labelsize": 6.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.62,
})

reg = regulatory_support()
all_umap = {key: load_umap(key) for key in KEY_TO_PAIR}
fdr = load_celltype_fdr()
overlap = pathway_overlap_records()

fig = plt.figure(figsize=(13.4, 9.4), dpi=300)
gs = GridSpec(4, 8, figure=fig, height_ratios=[1.0, 1.0, 0.92, 1.36], hspace=0.54, wspace=0.68)

axA = fig.add_subplot(gs[0:2, 0:2])
layers = ["GTEx SMR", "Bulk SMR", "Cell SMR", "cTWAS"]
x = np.arange(len(layers))
w = 0.24
for j, pair in enumerate(PAIR_ORDER):
    d = reg[reg["pair"].eq(pair)].set_index("layer").loc[layers]
    axA.bar(x + (j - 1) * w, d["n"], width=w, color=PAIR_COLORS[pair], label=PAIR_LABELS[pair])
axA.set_xticks(x)
axA.set_xticklabels(layers, rotation=30, ha="right")
axA.set_ylabel("Genes with support")
axA.set_title("Regulatory evidence", loc="left", fontweight="bold")
axA.legend(frameon=False, fontsize=6.2, loc="upper center", bbox_to_anchor=(0.50, 1.01), ncol=1, handlelength=1.8)
for s in ["top", "right"]:
    axA.spines[s].set_visible(False)
panel_label(axA, "A", x=-0.18, y=1.06)

axB = fig.add_subplot(gs[0, 2:5])
draw_umap_celltypes(axB, all_umap["HDL_core_AD"])
panel_label(axB, "B", x=-0.10, y=1.055)

axC = fig.add_subplot(gs[0, 5:8])
sc = draw_umap_trs(axC, all_umap["HDL_core_AD"], KEY_TITLES["HDL_core_AD"])
panel_label(axC, "C", x=-0.10, y=1.055)
cbar = fig.colorbar(sc, ax=axC, fraction=0.040, pad=0.010)
cbar.set_label("TRS percentile", fontsize=6.6)
cbar.set_ticks([0, 0.5, 1.0])
cbar.set_ticklabels(["low", "mid", "high"])
cbar.ax.tick_params(labelsize=5.8)

axD = fig.add_subplot(gs[1, 2:5])
draw_umap_trs(axD, all_umap["Ketone_core_PD"], KEY_TITLES["Ketone_core_PD"])
panel_label(axD, "D", x=-0.10, y=1.055)

axE = fig.add_subplot(gs[1, 5:8])
draw_umap_trs(axE, all_umap["TG_VLDL_core_PD"], KEY_TITLES["TG_VLDL_core_PD"])
panel_label(axE, "E", x=-0.10, y=1.055)

subV = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2, :], wspace=0.30)
selected = {
    "HDL_core_AD": ["Pericyte", "OPC", "Astrocyte", "Microglia", "Oligodendrocyte"],
    "Ketone_core_PD": ["GABA neuron", "OPC", "Astrocyte", "VIP GABA", "Pericyte"],
    "TG_VLDL_core_PD": ["Pericyte", "OPC", "Astrocyte", "Oligodendrocyte", "VIP GABA"],
}
raw_map = {
    "Pericyte": "pericyte",
    "OPC": "oligodendrocyte precursor cell",
    "Astrocyte": "astrocyte",
    "Microglia": "microglial cell",
    "Oligodendrocyte": "oligodendrocyte",
    "VIP GABA": "VIP GABAergic cortical interneuron",
    "GABA neuron": "GABAergic neuron",
    "Endothelial": "endothelial cell",
}
for i, key in enumerate(["HDL_core_AD", "Ketone_core_PD", "TG_VLDL_core_PD"]):
    ax = fig.add_subplot(subV[0, i])
    d = all_umap[key]
    cts = selected[key]
    vals = [d.loc[d["cell_type_short"].eq(ct), "scPagwas.TRS.Score"].dropna().to_numpy() for ct in cts]
    parts = ax.violinplot(vals, positions=np.arange(len(cts)), widths=0.78, showmeans=False, showextrema=False, showmedians=True)
    color = PAIR_COLORS[KEY_TO_PAIR[key]]
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.68)
    parts["cmedians"].set_color("#262626")
    parts["cmedians"].set_linewidth(0.7)
    ax.set_title(PAIR_LABELS[KEY_TO_PAIR[key]], color=color, fontweight="bold", loc="left")
    ax.set_xticks(np.arange(len(cts)))
    ax.set_xticklabels(cts, rotation=21, ha="right")
    ax.axhline(0, color="#CFCFCF", lw=0.55, zorder=0)
    ax.set_ylabel("TRS score" if i == 0 else "")
    ff = fdr[fdr["analysis_key"].eq(key)]
    ytop = max([np.nanpercentile(v, 98) if len(v) else 0 for v in vals])
    yrange = max(0.01, ytop - min([np.nanpercentile(v, 2) if len(v) else 0 for v in vals]))
    for j, ct in enumerate(cts):
        hit = ff[ff["celltype"].astype(str).eq(raw_map.get(ct, ct))]
        txt = star(float(hit["celltype_FDR"].iloc[0])) if len(hit) else ""
        if txt:
            ax.text(j, ytop + 0.05 * yrange, txt, ha="center", va="bottom", fontsize=8.0, fontweight="bold")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    panel_label(ax, chr(ord("F") + i), x=-0.11, y=1.055)

axI = fig.add_subplot(gs[3, 0:4])
plot_overlap(axI, overlap, "AD", "AD: overlapping scPagwas-KNK pathways", 8)
panel_label(axI, "I", x=-0.12, y=1.055)

axJ = fig.add_subplot(gs[3, 4:8])
plot_overlap(axJ, overlap, "PD", "PD: overlapping scPagwas-KNK pathways", 10)
panel_label(axJ, "J", x=-0.12, y=1.055)

handles = [
    Line2D([0], [0], marker="s", color="none", markerfacecolor=METHOD_COLORS["scPagwas"], markersize=6, label="scPagwas"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=METHOD_COLORS["KNK"], markersize=6, label="KNK"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=PAIR_COLORS["lipid8_F2_AD"], markersize=5.5, label="HDL-core / AD"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=PAIR_COLORS["nonlipid8_F1_PD"], markersize=5.5, label="Ketone-core / PD"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=PAIR_COLORS["lipid8_F1_PD"], markersize=5.5, label="TG/VLDL-core / PD"),
]
fig.legend(handles=handles, frameon=False, ncol=5, loc="lower center", bbox_to_anchor=(0.55, 0.010), fontsize=6.9)

fig.subplots_adjust(left=0.132, right=0.985, top=0.958, bottom=0.082)
fig.savefig(OUT_PNG, dpi=450)
fig.savefig(OUT_PDF)
fig.savefig(OUT_SVG)
print(OUT_PNG)
