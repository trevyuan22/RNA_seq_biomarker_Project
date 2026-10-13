"""Sample quality control, sex check, cell-type marker panels and the per-sample factor table."""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

SEX_Y = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "ZFY", "EIF1AY", "USP9Y", "TXLNGY", "NLGN4Y"]
MARKERS = {
    "Pluripotency": ["POU5F1", "NANOG", "LIN28A", "DPPA4", "TDGF1"],
    "Neural progenitor": ["SOX2", "NES", "PAX6", "SOX1", "HES1", "VIM"],
    "Neuron": ["MAP2", "TUBB3", "RBFOX3", "SYT1", "STMN2", "DCX", "SNAP25"],
    "Astrocyte": ["GFAP", "S100B", "AQP4", "ALDH1L1", "SLC1A2"],
    "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "MCM2"],
}


def sample_qc(X_raw, Xl, y, kind):
    d = pd.DataFrame(index=Xl.index)
    d["group"] = y.values
    if kind == "counts":
        tot = X_raw.sum(axis=1)
        d["library_size_M"] = tot / 1e6
        mt = [c for c in X_raw.columns if str(c).upper().startswith("MT-")]
        d["mito_pct"] = X_raw[mt].sum(axis=1) / tot * 100 if mt else np.nan
        top = np.sort(X_raw.values, axis=1)[:, -10:].sum(axis=1)
        d["top10_genes_pct"] = top / tot.values * 100
    d["genes_detected"] = (Xl >= 1).sum(axis=1)
    C = Xl.T.corr()
    same = []
    other = []
    for s in Xl.index:
        m = (y == y[s]).values
        same.append(C.loc[s, Xl.index[m & (Xl.index != s)]].mean())
        other.append(C.loc[s, Xl.index[~m]].mean())
    d["corr_own_group"], d["corr_other_group"] = same, other
    return d


def sex_check(Xl):
    cols = {c.upper(): c for c in Xl.columns}
    y_genes = [cols[g] for g in SEX_Y if g in cols]
    y_score = Xl[y_genes].mean(axis=1) if y_genes else pd.Series(np.nan, index=Xl.index)
    xist = Xl[cols["XIST"]] if "XIST" in cols else pd.Series(np.nan, index=Xl.index)
    ym = y_score.mean()
    if not y_genes:
        call = "No Y-chromosome genes found in the data, so sex could not be checked."
    elif ym >= 3:
        call = "Y-chromosome genes are clearly expressed: samples look male-derived."
    elif ym < 1.5:
        call = ("Y-chromosome genes are barely detected: samples look female-derived (or Y-gene expression is lost). "
                + ("XIST is " + ("expressed." if xist.mean() >= 2 else "low/absent, which in stem-cell-derived lines can indicate X-inactivation erosion.")
                   if "XIST" in cols else "XIST is not in the data table, so X-inactivation status could not be assessed."))
    else:
        call = "Y-chromosome gene expression is intermediate/ambiguous; sex could not be called reliably."
    return pd.DataFrame({"y_gene_score": y_score, "XIST": xist}), call, len(y_genes)


def marker_scores(Xl):
    cols = {c.upper(): c for c in Xl.columns}
    out, cover = {}, {}
    for panel, genes in MARKERS.items():
        have = [cols[g] for g in genes if g in cols]
        cover[panel] = f"{len(have)}/{len(genes)}"
        if have:
            out[panel] = Xl[have].mean(axis=1)
    return pd.DataFrame(out), cover


def pcs(Xl, n=3):
    top = Xl.var().sort_values(ascending=False).head(2000).index
    Z = Xl[top]
    Z = (Z - Z.mean()) / (Z.std() + 1e-9)
    k = min(n, len(Xl) - 1)
    p = PCA(k, random_state=0)
    c = p.fit_transform(Z.values)
    return (pd.DataFrame(c, index=Xl.index, columns=[f"PC{i + 1}" for i in range(k)]),
            p.explained_variance_ratio_)


def build_factors(y, case, qc_df, sex_df, markers, pc_df, extra=None):
    F = pd.DataFrame(index=y.index)
    F[f"group ({case}=1)"] = (y == case).astype(int).values
    for c in ("library_size_M", "genes_detected", "mito_pct", "top10_genes_pct"):
        if c in qc_df and qc_df[c].notna().any():
            F[c] = qc_df[c]
    F["Y-gene score"] = sex_df["y_gene_score"]
    if sex_df["XIST"].notna().any():
        F["XIST"] = sex_df["XIST"]
    for c in markers.columns:
        F[f"marker: {c}"] = markers[c]
    if extra is not None:
        for c in extra.columns:
            F[c] = extra[c]
    for c in pc_df.columns:
        F[c] = pc_df[c]
    return F
