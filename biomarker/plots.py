"""Figures -> base64 PNG strings (so the HTML report is one self-contained file).
Styled to match the report: warm paper background, one warm and one cool accent, no chart junk."""
import base64
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA

PAPER, INK, MUTED, RULE = "#fffdf9", "#1f2430", "#6b6f7a", "#e4dfd2"
UP, DOWN, NS, GOLD = "#b5412b", "#2b6a8c", "#cfc9bb", "#b7791f"
DIV = LinearSegmentedColormap.from_list("div", [DOWN, "#f6f2e8", UP])
SEQ = LinearSegmentedColormap.from_list("seq", ["#fbf8f0", "#9fc0cf", DOWN])
PALETTE = ["#2b6a8c", "#b5412b", "#6b8e4e", "#8a5a9e", "#c28a1e", "#3f8f8f", "#7d7d7d", "#a05a6a", "#4f6fa8", "#9a7b4f"] + \
    list(plt.get_cmap("tab20").colors)

plt.rcParams.update({
    "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
    "axes.edgecolor": "#a39d8c", "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titleweight": "bold", "axes.titlesize": 10, "axes.labelsize": 9, "font.size": 9,
    "font.family": "sans-serif", "grid.color": RULE, "legend.frameon": False,
})


def _b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def color_map(classes):
    return {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(classes)}


def pca_plot(Xl, y):
    top = Xl.var().sort_values(ascending=False).head(2000).index
    Z = Xl[top]
    Z = (Z - Z.mean()) / (Z.std() + 1e-9)
    p = PCA(2, random_state=0)
    c = p.fit_transform(Z.values)
    cm = color_map(sorted(y.unique()))
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    for k in sorted(y.unique()):
        m = (y == k).values
        ax.scatter(c[m, 0], c[m, 1], s=22, alpha=.8, color=cm[k], label=k)
    ax.set_xlabel(f"PC1 ({p.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({p.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.legend(fontsize=8, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_title("Sample map (PCA, 2,000 most variable genes)")
    return _b64(fig)


def volcano(df, title, cfg):
    fig, ax = plt.subplots(figsize=(6, 4.6))
    ylog = -np.log10(df.padj.clip(lower=1e-300))
    col = df.direction.map({"up": UP, "down": DOWN, "ns": NS})
    ax.scatter(df.log2FC, ylog, s=7, c=col, alpha=.7, linewidths=0)
    ax.axhline(-np.log10(cfg.padj_cutoff), ls="--", lw=.7, c="k")
    for x in (-cfg.lfc_cutoff, cfg.lfc_cutoff):
        ax.axvline(x, ls="--", lw=.7, c="k")
    top = df[df.significant].assign(y=ylog[df.significant]).sort_values("y", ascending=False).head(8)
    for _, r in top.iterrows():
        ax.annotate(r.gene, (r.log2FC, r.y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("log2 fold-change"); ax.set_ylabel("-log10 adjusted p")
    ax.set_title(title, fontsize=10)
    return _b64(fig)


def marker_heatmap(Xl, y, panels):
    genes = list(dict.fromkeys(g for p in panels.values() for g in p.gene if g in Xl.columns))
    if not genes:
        return None
    classes = sorted(y.unique())
    M = pd.DataFrame({c: Xl.loc[(y == c).values, genes].mean() for c in classes})
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1) + 1e-9, axis=0)
    fig, ax = plt.subplots(figsize=(max(4, .9 * len(classes) + 2), max(4, .17 * len(genes) + 1.5)))
    im = ax.imshow(Z.values, aspect="auto", cmap=DIV, vmin=-2, vmax=2)
    ax.set_xticks(range(len(classes))); ax.set_xticklabels(classes, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
    fig.colorbar(im, ax=ax, shrink=.5, label="relative expression (z-score)")
    ax.set_title("Biomarker panel: mean expression per disease")
    return _b64(fig)


def enrichment_bar(df, title, n=10):
    d = df.head(n).iloc[::-1]
    if d.empty:
        return None
    fig, ax = plt.subplots(figsize=(7, max(2.2, .38 * len(d) + 1)))
    ax.barh([t[:60] for t in d.term], -np.log10(d.padj.astype(float).clip(lower=1e-300)), color=DOWN)
    ax.axvline(-np.log10(0.05), ls="--", c="k", lw=.7)
    ax.set_xlabel("-log10 adjusted p"); ax.set_title(title, fontsize=10)
    ax.tick_params(axis="y", labelsize=7)
    return _b64(fig)


def confusion_plot(cm, classes):
    fig, ax = plt.subplots(figsize=(max(4, .8 * len(classes) + 1.5),) * 2)
    ax.imshow(cm, cmap=SEQ)
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, cm[i, j], ha="center", va="center", color="white" if cm[i, j] > cm.max() / 2 else "black")
    ax.set_xticks(range(len(classes))); ax.set_xticklabels(classes, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(classes))); ax.set_yticklabels(classes, fontsize=8)
    ax.set_xlabel("Predicted"); ax.set_ylabel("True"); ax.set_title("Held-out test samples")
    return _b64(fig)


def prob_bars(pred_df, classes, max_n=40):
    d = pred_df.head(max_n).iloc[::-1]
    cm = color_map(classes)
    fig, ax = plt.subplots(figsize=(7, max(2.5, .32 * len(d) + 1)))
    left = np.zeros(len(d))
    for c in classes:
        v = d[f"P({c})"].values
        ax.barh(d["sample"], v, left=left, color=cm[c], label=c)
        left += v
    ax.set_xlim(0, 1); ax.set_xlabel("probability")
    ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.tick_params(axis="y", labelsize=7)
    return _b64(fig)


# ---------------------------------------------------------------- comprehensive-analysis figures
def corr_heatmap(C, title, vmin=-1.0, cmap=DIV, annotate=None):
    n = len(C)
    size = max(4.5, .42 * n + 2.5)
    fig, ax = plt.subplots(figsize=(size, size * .9))
    im = ax.imshow(C.values.astype(float), cmap=cmap, vmin=vmin, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(C.columns, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(C.index))); ax.set_yticklabels(C.index, fontsize=7)
    if annotate or (annotate is None and n <= 14):
        for i in range(len(C.index)):
            for j in range(n):
                ax.text(j, i, f"{C.values[i, j]:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, shrink=.6)
    ax.set_title(title, fontsize=10)
    return _b64(fig)


def rect_heatmap(M, title, cmap=DIV, vmin=-1, vmax=1, label="rho"):
    fig, ax = plt.subplots(figsize=(max(4, .55 * M.shape[1] + 3), max(3, .26 * M.shape[0] + 1.5)))
    im = ax.imshow(M.values.astype(float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=7)
    fig.colorbar(im, ax=ax, shrink=.6, label=label)
    ax.set_title(title, fontsize=10)
    return _b64(fig)


def pca_labeled(pc_df, var, y):
    cm = color_map(sorted(y.unique()))
    fig, ax = plt.subplots(figsize=(6, 4.6))
    for k in sorted(y.unique()):
        m = (y == k).values
        ax.scatter(pc_df.iloc[m, 0], pc_df.iloc[m, 1], s=45, color=cm[k], label=k)
    for s_, r in pc_df.iterrows():
        ax.annotate(str(s_).replace("LLS-RNA-seq-", "S"), (r.iloc[0], r.iloc[1]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)"); ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Samples (labelled)", fontsize=10)
    return _b64(fig)


def de_heatmap(Xl, y, genes, n=40):
    genes = [g for g in genes if g in Xl.columns][:n]
    if len(genes) < 3:
        return None
    from scipy.cluster.hierarchy import leaves_list, linkage
    order_s = list(y.sort_values(kind="stable").index)
    Z = Xl.loc[order_s, genes]
    Z = (Z - Z.mean()) / (Z.std(ddof=0) + 1e-9)
    lk = leaves_list(linkage(Z.T.values, "average"))
    Z = Z.iloc[:, lk]
    fig, ax = plt.subplots(figsize=(max(4, .6 * len(order_s) + 2.5), max(4, .2 * len(genes) + 1.5)))
    im = ax.imshow(Z.T.values, aspect="auto", cmap=DIV, vmin=-2, vmax=2)
    ax.set_xticks(range(len(order_s)))
    ax.set_xticklabels([f"{str(s).replace('LLS-RNA-seq-', 'S')}\n{y[s]}" for s in order_s], fontsize=7)
    ax.set_yticks(range(len(Z.columns))); ax.set_yticklabels(Z.columns, fontsize=7)
    fig.colorbar(im, ax=ax, shrink=.5, label="z-score across samples")
    ax.set_title("Top differentially expressed genes", fontsize=10)
    return _b64(fig)


def qc_bars(qc_df):
    cols = [c for c in ("library_size_M", "genes_detected", "mito_pct") if c in qc_df and qc_df[c].notna().any()]
    cm = color_map(sorted(qc_df.group.unique()))
    fig, axs = plt.subplots(1, len(cols), figsize=(4.2 * len(cols), 3.4), squeeze=False)
    for ax, c in zip(axs[0], cols):
        ax.bar(range(len(qc_df)), qc_df[c], color=[cm[g] for g in qc_df.group])
        ax.set_xticks(range(len(qc_df)))
        ax.set_xticklabels([str(s).replace("LLS-RNA-seq-", "S") for s in qc_df.index], rotation=60, fontsize=7)
        ax.set_title(c, fontsize=9)
    from matplotlib.patches import Patch
    axs[0][0].legend(handles=[Patch(color=cm[g], label=g) for g in sorted(cm)], fontsize=7, loc="lower left")
    return _b64(fig)


def chrom_plot(merged, tab):
    chroms = [c for c in tab.chromosome]
    data = [merged.loc[merged.chromosome == c, "log2FC"].clip(-5, 5).values for c in chroms]
    fig, ax = plt.subplots(figsize=(9, 3.8))
    bp = ax.boxplot(data, showfliers=False, patch_artist=True)      # works on every matplotlib version
    ax.set_xticks(range(1, len(chroms) + 1))
    ax.set_xticklabels(chroms)
    flag = set(tab.loc[tab.padj_enrichment < 0.05, "chromosome"])
    for patch, c in zip(bp["boxes"], chroms):
        patch.set_facecolor(UP if c in flag else "#b9c9d3")
    ax.axhline(0, c="k", lw=.6)
    ax.set_ylabel("log2 fold-change"); ax.set_xlabel("chromosome")
    ax.set_title("Fold-change by chromosome (red = significant genes over-represented)", fontsize=10)
    return _b64(fig)


def gsea_plot(df, title, n=8, fdr=0.25):
    if df is None or len(df) == 0:
        return None
    d = df.head(60)
    pos = d[d.NES > 0].head(n)
    neg = d[d.NES < 0].head(n)
    d = pd.concat([pos, neg]).sort_values("NES")
    fig, ax = plt.subplots(figsize=(8, max(2.5, .33 * len(d) + 1)))
    ax.barh([t[:60] + (" *" if p < fdr else "") for t, p in zip(d.term, d.padj)], d.NES,
            color=[UP if v > 0 else DOWN for v in d.NES])
    ax.axvline(0, c="k", lw=.6)
    ax.set_xlabel("normalised enrichment score (red = higher in case; * = FDR passes)")
    ax.tick_params(axis="y", labelsize=7)
    ax.set_title(title, fontsize=10)
    return _b64(fig)


def network_plot(edges, nodes, de):
    dirs = dict(zip(de.gene.astype(str).str.upper(), de.direction))
    names = list(nodes.gene)
    try:
        import networkx as nx
        g = nx.from_pandas_edgelist(edges, "a", "b")
        pos = nx.spring_layout(g, seed=1, k=1.2 / np.sqrt(max(len(g), 1)))
        pts = np.array([pos[n] for n in names])
    except ImportError:
        ang = np.linspace(0, 2 * np.pi, len(names), endpoint=False)
        pts = np.c_[np.cos(ang), np.sin(ang)]
    idx = {n: i for i, n in enumerate(names)}
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for a, b in zip(edges.a, edges.b):
        if a in idx and b in idx:
            ax.plot(*zip(pts[idx[a]], pts[idx[b]]), c="#cccccc", lw=.5, zorder=1)
    col = [{"up": UP, "down": DOWN}.get(dirs.get(n, "ns"), "#a39d8c") for n in names]
    ax.scatter(pts[:, 0], pts[:, 1], s=20 + 14 * nodes.degree.values, c=col, zorder=2, edgecolors="white", linewidths=.5)
    for i in np.argsort(-nodes.degree.values)[:20]:
        ax.annotate(names[i], pts[i], fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axis("off")
    ax.set_title("Protein-interaction network (red = up, blue = down; size = connections)", fontsize=10)
    return _b64(fig)


def rna_ribo_scatter(m, st):
    cm = {"RNA and translation agree": NS, "translation UP beyond RNA": UP, "translation DOWN beyond RNA": DOWN}
    fig, ax = plt.subplots(figsize=(5.6, 5))
    for k, c in cm.items():
        s = m[m.translation_class == k]
        ax.scatter(s.log2FC_rna, s.log2FC_ribo, s=6, c=c, alpha=.6, label=f"{k} ({len(s)})", linewidths=0)
    lim = 6
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=.6)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("RNA-seq log2 fold-change"); ax.set_ylabel("Ribo-seq log2 fold-change")
    ax.set_title(f"Transcription vs translation (Spearman r = {st['spearman']:.2f})", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    return _b64(fig)
