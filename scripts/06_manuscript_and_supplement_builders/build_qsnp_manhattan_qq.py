from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
REF_DIR = Path(r"D:\LDSC\ldsc-master\eur_w_ld_chr")
OUT_DIR = ROOT / "q_snp_figures_final8"
OUT_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="white", context="talk")
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["figure.dpi"] = 180

MODULE_COLORS = {
    "lipid": ("#B24745", "#E7B6B2"),
    "nonlipid": ("#2F6C8F", "#B7D0E1"),
}


def save_fig(fig, stem):
    fig.savefig(OUT_DIR / f"{stem}.png", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def load_reference_map():
    parts = []
    for chrom in range(1, 23):
        path = REF_DIR / f"{chrom}.l2.ldscore.gz"
        df = pd.read_csv(path, sep="\t", compression="gzip", usecols=["CHR", "SNP", "BP"])
        parts.append(df)
    ref = pd.concat(parts, ignore_index=True).drop_duplicates("SNP")
    ref["CHR"] = ref["CHR"].astype(np.int64)
    ref["BP"] = ref["BP"].astype(np.int64)
    ref = ref.sort_values(["CHR", "BP"]).reset_index(drop=True)

    chr_max = ref.groupby("CHR", as_index=False)["BP"].max().rename(columns={"BP": "chr_len"})
    chr_max["offset"] = chr_max["chr_len"].cumsum().shift(fill_value=0).astype(np.int64)
    ref = ref.merge(chr_max[["CHR", "offset"]], on="CHR", how="left")
    ref["pos"] = (ref["BP"] + ref["offset"]).astype(np.int64)

    centers = ref.groupby("CHR", as_index=False).agg(start=("pos", "min"), end=("pos", "max"))
    centers["center"] = (centers["start"] + centers["end"]) / 2
    return ref[["SNP", "CHR", "BP", "pos"]], centers[["CHR", "center"]]


def load_qsnp_for_plot(path, ref_map):
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=["SNP", "Q_SNP_pval"],
        low_memory=False,
    )
    df = df.rename(columns={"Q_SNP_pval": "P"})
    df = df.merge(ref_map, on="SNP", how="inner")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df = df.dropna(subset=["P", "CHR", "BP"])
    df = df[(df["P"] > 0) & (df["P"] <= 1)].copy()
    df["minus_log10_p"] = -np.log10(df["P"])
    return df.sort_values(["CHR", "BP"])


def qq_data(pvals):
    pvals = np.asarray(pvals, dtype=float)
    pvals = pvals[np.isfinite(pvals) & (pvals > 0) & (pvals <= 1)]
    pvals = np.sort(pvals)
    n = len(pvals)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(pvals)
    return exp, obs


def add_panel_label(ax, label):
    ax.text(-0.10, 1.05, label, transform=ax.transAxes, fontsize=17, fontweight="bold", va="bottom")


def make_manhattan_panel(module, files, ref_map, chr_centers):
    dark, light = MODULE_COLORS[module]
    fig, axes = plt.subplots(3, 1, figsize=(18, 10), sharex=True, gridspec_kw={"hspace": 0.12})
    chr_colors = {chrom: (dark if chrom % 2 else light) for chrom in range(1, 23)}
    labels = ["A", "B", "C"]

    ymax = 0
    plot_data = []
    for _, _, path in files:
        df = load_qsnp_for_plot(path, ref_map)
        ymax = max(ymax, np.nanpercentile(df["minus_log10_p"], 99.95), 8)
        plot_data.append(df)
    ymax = max(ymax + 0.8, 8)

    for ax, label, (factor, title, _), df in zip(axes, labels, files, plot_data):
        for chrom in range(1, 23):
            sub = df[df["CHR"] == chrom]
            if len(sub) == 0:
                continue
            ax.scatter(sub["pos"], sub["minus_log10_p"], s=4, c=chr_colors[chrom], alpha=0.75, linewidths=0, rasterized=True)
        ax.axhline(-np.log10(5e-8), color="black", linestyle="--", linewidth=1)
        top = df.nsmallest(1, "P").iloc[0]
        ax.text(0.995, 0.92, f"{factor}: min Q_SNP P = {top['P']:.2e}", transform=ax.transAxes, ha="right", va="top", fontsize=10, color=dark)
        ax.set_ylim(0, ymax)
        ax.set_ylabel("-log10(Q_SNP P)")
        ax.set_title(title, loc="left", color=dark, fontsize=15)
        ax.grid(axis="y", alpha=0.15)
        sns.despine(ax=ax)
        add_panel_label(ax, label)

    odd_chr = chr_centers[chr_centers["CHR"] % 2 == 1]
    axes[-1].set_xticks(odd_chr["center"])
    axes[-1].set_xticklabels(odd_chr["CHR"])
    axes[-1].set_xlabel("Chromosome")
    fig.suptitle(f"{module.capitalize()} Q_SNP Manhattan plots", y=1.01, fontsize=20, fontweight="bold")
    save_fig(fig, f"Figure_{module}_qsnp_manhattan")


def make_qq_panel(module, files, ref_map):
    dark, light = MODULE_COLORS[module]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), gridspec_kw={"wspace": 0.28})
    labels = ["A", "B", "C"]

    max_axis = 0
    qq_series = []
    for _, _, path in files:
        df = load_qsnp_for_plot(path, ref_map)
        exp, obs = qq_data(df["P"].values)
        qq_series.append((df, exp, obs))
        max_axis = max(max_axis, np.nanpercentile(obs, 99.9), np.nanpercentile(exp, 99.9))
    max_axis = max(max_axis, 8)

    for ax, panel_label, (factor, title, _), (df, exp, obs) in zip(axes, labels, files, qq_series):
        ax.scatter(exp, obs, s=5, c=dark, alpha=0.6, linewidths=0, rasterized=True)
        ax.plot([0, max_axis], [0, max_axis], linestyle="--", color=light, linewidth=2)
        top = df.nsmallest(1, "P").iloc[0]
        ax.text(0.98, 0.06, f"n = {len(df):,}\nmin P = {top['P']:.2e}", transform=ax.transAxes, ha="right", va="bottom", fontsize=10, color=dark)
        ax.set_xlim(0, max_axis)
        ax.set_ylim(0, max_axis)
        ax.set_xlabel("Expected -log10(P)")
        ax.set_ylabel("Observed -log10(P)")
        ax.set_title(title, color=dark, fontsize=15)
        sns.despine(ax=ax)
        add_panel_label(ax, panel_label)

    fig.suptitle(f"{module.capitalize()} Q_SNP QQ plots", y=1.03, fontsize=20, fontweight="bold")
    save_fig(fig, f"Figure_{module}_qsnp_qq")


def write_manifest():
    manifest = pd.DataFrame(
        [
            ("Figure_lipid_qsnp_manhattan.png", "Three-panel Manhattan plot for lipid Q_SNP results across F1-F3."),
            ("Figure_lipid_qsnp_qq.png", "Three-panel QQ plot for lipid Q_SNP results across F1-F3."),
            ("Figure_nonlipid_qsnp_manhattan.png", "Three-panel Manhattan plot for nonlipid Q_SNP results across F1-F3."),
            ("Figure_nonlipid_qsnp_qq.png", "Three-panel QQ plot for nonlipid Q_SNP results across F1-F3."),
        ],
        columns=["file", "description"],
    )
    manifest.to_csv(OUT_DIR / "q_snp_figure_manifest.tsv", sep="\t", index=False)


def main():
    ref_map, chr_centers = load_reference_map()
    lipid_files = [
        ("F1", "Lipid F1 Q_SNP", ROOT / "step14_native_wsl_usergwas_final8_results" / "merged_lipid_final8" / "lipid_final8_F1_userGWAS_merged.tsv.gz"),
        ("F2", "Lipid F2 Q_SNP", ROOT / "step14_native_wsl_usergwas_final8_results" / "merged_lipid_final8" / "lipid_final8_F2_userGWAS_merged.tsv.gz"),
        ("F3", "Lipid F3 Q_SNP", ROOT / "step14_native_wsl_usergwas_final8_results" / "merged_lipid_final8" / "lipid_final8_F3_userGWAS_merged.tsv.gz"),
    ]
    nonlipid_files = [
        ("F1", "Nonlipid F1 Q_SNP", ROOT / "step22_native_wsl_usergwas_nonlipid8_results" / "merged_nonlipid_final8" / "nonlipid_final8_F1_userGWAS_merged.tsv.gz"),
        ("F2", "Nonlipid F2 Q_SNP", ROOT / "step22_native_wsl_usergwas_nonlipid8_results" / "merged_nonlipid_final8" / "nonlipid_final8_F2_userGWAS_merged.tsv.gz"),
        ("F3", "Nonlipid F3 Q_SNP", ROOT / "step22_native_wsl_usergwas_nonlipid8_results" / "merged_nonlipid_final8" / "nonlipid_final8_F3_userGWAS_merged.tsv.gz"),
    ]

    make_manhattan_panel("lipid", lipid_files, ref_map, chr_centers)
    make_qq_panel("lipid", lipid_files, ref_map)
    make_manhattan_panel("nonlipid", nonlipid_files, ref_map, chr_centers)
    make_qq_panel("nonlipid", nonlipid_files, ref_map)
    write_manifest()
    print(f"Wrote Q_SNP Manhattan/QQ figures to: {OUT_DIR}")


if __name__ == "__main__":
    main()
