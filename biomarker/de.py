"""Differential expression. DESeq2 (via PyDESeq2) for raw counts when feasible,
otherwise a fast Welch t-test on log2 values. Always BH-adjusted."""
import warnings
import numpy as np
import pandas as pd
from scipy import stats


def bh(p):
    """Benjamini-Hochberg FDR."""
    p = np.nan_to_num(np.asarray(p, dtype=float), nan=1.0)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def make_contrasts(y, control=None):
    classes = sorted(y.unique())
    out = []
    if control:
        if control not in classes:
            raise ValueError(f"--control '{control}' is not one of the labels: {classes}")
        for c in classes:
            if c != control:
                out.append(dict(name=f"{c}_vs_{control}", cls=c,
                                case=(y == c).values, ref=(y == control).values))
    else:
        for c in classes:
            out.append(dict(name=f"{c}_vs_rest", cls=c, case=(y == c).values, ref=(y != c).values))
    return out


def _welch(Xl, case, ref):
    a, b = Xl.values[case], Xl.values[ref]
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, p = stats.ttest_ind(a, b, equal_var=False, axis=0)
    with np.errstate(all="ignore"):
        se = np.sqrt(a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b))
    return a.mean(0) - b.mean(0), np.nan_to_num(p, nan=1.0), bh(p), np.nan_to_num(se, nan=1.0)


def _moderated(Xl, case, ref, d0=4.0):
    """Simplified limma-style moderated t-test: shrinks each gene's variance toward the typical
    variance. Much more stable than a plain t-test when groups have only 2-4 samples."""
    a, b = Xl.values[case], Xl.values[ref]
    n1, n2 = len(a), len(b)
    df = n1 + n2 - 2
    s2 = ((n1 - 1) * a.var(0, ddof=1) + (n2 - 1) * b.var(0, ddof=1)) / df
    s0 = np.median(s2[s2 > 0])
    s2m = (d0 * s0 + df * s2) / (d0 + df)
    d = a.mean(0) - b.mean(0)
    t = d / np.sqrt(s2m * (1 / n1 + 1 / n2))
    pv = 2 * stats.t.sf(np.abs(t), df + d0)
    return d, pv, bh(pv), np.sqrt(s2m * (1 / n1 + 1 / n2))


def _deseq2(counts, case, ref):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    keep = case | ref
    c = counts.loc[keep].round().astype(int)
    meta = pd.DataFrame({"condition": np.where(case[keep], "case", "ref")}, index=c.index)
    try:
        dds = DeseqDataSet(counts=c, metadata=meta, design="~condition", quiet=True)
    except TypeError:                       # older PyDESeq2
        dds = DeseqDataSet(counts=c, metadata=meta, design_factors="condition", quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "case", "ref"], quiet=True)
    st.summary()
    r = st.results_df.reindex(counts.columns)
    return (r["log2FoldChange"].values, r["pvalue"].values, r["padj"].fillna(1).values,
            r["lfcSE"].fillna(1).values)


def run_de(Xl, y, contrasts, cfg, counts=None):
    """Return {contrast_name: DataFrame}, and a note on the method used."""
    results, methods, warned = {}, set(), False
    for c in contrasts:
        case, ref = c["case"], c["ref"]
        lfc = p = padj = se = None
        if counts is not None and (case | ref).sum() <= cfg.max_samples_deseq2:
            try:
                lfc, p, padj, se = _deseq2(counts, case, ref)
                methods.add("DESeq2 (PyDESeq2)")
            except Exception as e:           # missing package, convergence problem, etc.
                if not warned:
                    print(f"  [note] DESeq2 not used ({type(e).__name__}: install pydeseq2 to enable); using Welch t-test")
                    warned = True
        if lfc is None and min(case.sum(), ref.sum()) < 5:
            lfc, p, padj, se = _moderated(Xl, case, ref)
            methods.add("moderated t-test on log2 values (simplified limma-style; groups are small)")
        if lfc is None:
            lfc, p, padj, se = _welch(Xl, case, ref)
            methods.add("Welch t-test on log2 values")
        df = pd.DataFrame({
            "gene": Xl.columns, "log2FC": lfc, "lfcSE": se, "pvalue": p, "padj": padj,
            "mean_case": Xl.values[case].mean(0), "mean_ref": Xl.values[ref].mean(0)})
        df["log2FC"] = df["log2FC"].fillna(0.0)
        df["padj"] = df["padj"].fillna(1.0)
        df["significant"] = (df.padj < cfg.padj_cutoff) & (df.log2FC.abs() >= cfg.lfc_cutoff)
        df["direction"] = np.where(~df.significant, "ns", np.where(df.log2FC > 0, "up", "down"))
        df["_o"] = df.log2FC.abs()
        results[c["name"]] = df.sort_values(["padj", "_o"], ascending=[True, False]).drop(columns="_o").reset_index(drop=True)
    return results, " + ".join(sorted(methods))
