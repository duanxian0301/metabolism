from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


MATCHED = Path(r"D:\metabolic\233\ldsc_univariate\external_bivariate_final_summary\external_bivariate_matched_pairs_round1_round2.tsv")
SUMMARY = Path(r"D:\metabolic\233\ldsc_univariate\external_bivariate_final_summary\external_bivariate_factor_summary_round1_round2.tsv")
OUT_DIR = Path(r"D:\codex\GenomicSEM\metabolic\manuscript\figures_external_validation")
OUT_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["figure.dpi"] = 180

MODULE_COLORS = {"lipid": "#B24745", "nonlipid": "#2F6C8F"}
PAIRSET_MARKERS = {"primary_match": "o", "support_match": "s"}

FACTOR_DISPLAY = {
    "lipid8_F1": "Lipid F1",
    "lipid8_F2": "Lipid F2",
    "lipid8_F3": "Lipid F3",
    "nonlipid8_F1": "Nonlipid F1",
    "nonlipid8_F2": "Nonlipid F2",
    "nonlipid8_F3": "Nonlipid F3",
}


def add_panel_label(ax, label):
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=17, fontweight="bold", va="bottom")


def save_fig(fig, stem):
    fig.savefig(OUT_DIR / f"{stem}.png", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def prepare_pairs():
    df = pd.read_csv(MATCHED, sep="\t")
    df["display_factor"] = df["factor"].map(FACTOR_DISPLAY)
    df["label"] = df["display_factor"] + " | " + df["internal_trait"] + " <- " + df["external_trait"]
    df["ci_low"] = df["aligned_rg"] - 1.96 * df["rg_se"]
    df["ci_high"] = df["aligned_rg"] + 1.96 * df["rg_se"]
    factor_order = ["lipid8_F3", "lipid8_F2", "lipid8_F1", "nonlipid8_F3", "nonlipid8_F2", "nonlipid8_F1"]
    df["factor_rank"] = df["factor"].map({f: i for i, f in enumerate(factor_order)})
    df["pair_rank"] = df["pair_set"].map({"support_match": 0, "primary_match": 1})
    df = df.sort_values(["module", "factor_rank", "pair_rank", "aligned_rg"], ascending=[True, True, True, False]).reset_index(drop=True)
    return df


def prepare_summary():
    df = pd.read_csv(SUMMARY, sep="\t")
    df = df[df["pair_set"].isin(["primary_match", "support_match"])].copy()
    df["display_factor"] = df["factor"].map(FACTOR_DISPLAY)
    df["pair_label"] = df["pair_set"].map({"primary_match": "Primary", "support_match": "Support"})
    order = ["Lipid F1", "Lipid F2", "Lipid F3", "Nonlipid F1", "Nonlipid F2", "Nonlipid F3"]
    pair_order = ["Primary", "Support"]
    mat = df.pivot(index="display_factor", columns="pair_label", values="mean_aligned_rg")
    return mat.reindex(index=order, columns=pair_order)


def plot_forest(ax, sub, title, module):
    y = np.arange(len(sub))[::-1]
    for yi, (_, row) in zip(y, sub.iterrows()):
        ax.hlines(yi, row["ci_low"], row["ci_high"], color=MODULE_COLORS[module], linewidth=2.2, alpha=0.85)
        ax.scatter(
            row["aligned_rg"],
            yi,
            s=90,
            color=MODULE_COLORS[module],
            marker=PAIRSET_MARKERS[row["pair_set"]],
            edgecolor="black",
            linewidth=0.6,
            zorder=3,
        )
    ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["label"], fontsize=10)
    ax.set_xlim(-0.2, 1.05)
    ax.set_xlabel("Orientation-aligned genetic correlation")
    ax.set_title(title, color=MODULE_COLORS[module])
    sns.despine(ax=ax, left=False, bottom=False)


def write_legend():
    text = """Figure 5. External validation of lipid and nonlipid factor GWAS.

Panel A and Panel B show orientation-aligned genetic correlations (aligned rg) from bivariate LDSC between each final factor GWAS and its matched external metabolite GWAS from Karjalainen et al. Primary matches are shown as circles and same-domain support matches as squares. Horizontal lines indicate 95% confidence intervals.

Panel C summarizes factor-level mean aligned rg across primary and support matches. Lipid factors showed strong external correspondence across primary matches, with additional support for the CE-structural axis. Nonlipid F1 and F2 showed clear external support, whereas nonlipid F3 was supported mainly by glucose and remained weaker across secondary support traits. For nonlipid F3, direction was aligned before summary because factor sign is arbitrary.
"""
    (OUT_DIR / "Figure5_external_validation_legend_zh.md").write_text(text, encoding="utf-8")


def main():
    pairs = prepare_pairs()
    summary = prepare_summary()

    lipid = pairs[pairs["module"] == "lipid"].copy()
    nonlipid = pairs[pairs["module"] == "nonlipid"].copy()

    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.25, 0.9], hspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])

    plot_forest(ax1, lipid, "Lipid external matched-trait validation", "lipid")
    plot_forest(ax2, nonlipid, "Nonlipid external matched-trait validation", "nonlipid")

    sns.heatmap(
        summary,
        ax=ax3,
        cmap=sns.light_palette("#2F6C8F", as_cmap=True),
        vmin=0,
        vmax=0.95,
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={"label": "Mean aligned rg", "shrink": 0.8},
    )
    ax3.set_title("Factor-level summary across external primary and support matches")
    ax3.set_xlabel("")
    ax3.set_ylabel("")
    ax3.tick_params(axis="x", rotation=0)
    ax3.tick_params(axis="y", rotation=0)

    add_panel_label(ax1, "A")
    add_panel_label(ax2, "B")
    add_panel_label(ax3, "C")
    fig.suptitle("External validation of lipid and nonlipid factor GWAS", y=1.01, fontsize=20, fontweight="bold")
    save_fig(fig, "Figure5_external_validation")
    write_legend()
    print(f"Wrote figure outputs to: {OUT_DIR}")


if __name__ == "__main__":
    main()
