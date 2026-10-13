"""Integrate Ribo-seq (translation) with RNA-seq (transcription) for the same comparison.
Genes whose protein-synthesis change differs from their RNA change point to translational regulation."""
import numpy as np
import pandas as pd
from scipy import stats
from . import de, io_utils as io, preprocess
from .de import bh


def _map_groups(y, case, control):
    """Map whatever group names the Ribo-seq samples carry onto the case/control labels."""
    kc, kt = io.norm_key(case), io.norm_key(control)
    out = {}
    for s, g in y.items():
        k = io.norm_key(g)
        if k == kc or (kc in k and kt not in k):
            out[s] = case
        elif k == kt or (kt in k and kc not in k):
            out[s] = control
    return pd.Series(out, dtype=object)


def run_ribo(ribo_expression, ribo_labels, case, control, rna_de, cfg, input_type="auto", group_column=None):
    X = io.load_expression(ribo_expression, sample_hint=io.sheet_hint(ribo_labels))
    y = None
    if ribo_labels is not None:
        try:
            y = _map_groups(io.labels_for(list(X.index), ribo_labels, group_column), case, control)
        except ValueError:
            y = None
    if y is None or len(y) < 4:                       # fall back: group names written inside the sample names
        direct = pd.Series({s: (case if io.norm_key(case) in io.norm_key(s) and io.norm_key(control) not in io.norm_key(s)
                                else control if io.norm_key(control) in io.norm_key(s) and io.norm_key(case) not in io.norm_key(s)
                                else None) for s in X.index}, dtype=object).dropna()
        if len(direct) >= 4:
            y = direct
        else:
            guess, _ = io.infer_groups_from_names(list(X.index))
            y = _map_groups(guess, case, control) if guess is not None else pd.Series(dtype=object)
    if len(y) < 4:
        raise ValueError(f"Could not tell which Ribo-seq samples belong to '{case}' and '{control}'. "
                         "Include a sample sheet that lists the Ribo-seq samples too, or put the group names in the sample names.")
    X = X.loc[y.index]
    gc = y.value_counts().to_dict()
    if min(gc.get(case, 0), gc.get(control, 0)) < 2:
        raise ValueError(f"Ribo-seq needs at least 2 samples per group (found {gc}).")
    kind = preprocess.detect_scale(X) if input_type == "auto" else input_type
    Xf = preprocess.filter_genes(X, kind, cfg)
    Xl = preprocess.to_log2(Xf, kind, X.sum(axis=1) if kind == "counts" else None)
    con = [dict(name="ribo", cls=case, case=(y == case).values, ref=(y == control).values)]
    res, method = de.run_de(Xl, y, con, cfg, counts=Xf if kind == "counts" else None)
    rb = res["ribo"]
    m = rna_de[["gene", "log2FC", "lfcSE", "padj"]].merge(rb[["gene", "log2FC", "lfcSE", "padj"]], on="gene", suffixes=("_rna", "_ribo"))
    pear = stats.pearsonr(m.log2FC_rna, m.log2FC_ribo)
    spear = stats.spearmanr(m.log2FC_rna, m.log2FC_ribo)
    delta = m.log2FC_ribo - m.log2FC_rna
    m["delta_TE"] = delta - delta.median()        # remove any global library-composition shift
    se = np.sqrt(m.lfcSE_rna ** 2 + m.lfcSE_ribo ** 2).replace(0, np.nan)
    z = (m.delta_TE / se).fillna(0)
    m["p_TE"] = 2 * stats.norm.sf(np.abs(z))
    m["padj_TE"] = bh(m.p_TE.values)
    m["translation_class"] = np.where(m.padj_TE >= cfg.padj_cutoff, "RNA and translation agree",
                                      np.where(m.delta_TE > 0, "translation UP beyond RNA", "translation DOWN beyond RNA"))
    stats_ = dict(n_genes=len(m), pearson=float(pear[0]), spearman=float(spear[0]), method=method,
                  n_ribo_samples=len(y), n_te_sig=int((m.padj_TE < cfg.padj_cutoff).sum()))
    return m.sort_values("padj_TE").reset_index(drop=True), rb, stats_
