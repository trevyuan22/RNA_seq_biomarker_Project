"""Quality filtering and normalisation."""
import numpy as np
import pandas as pd


def detect_scale(X):
    """'counts' (raw integers), 'linear' (TPM/FPKM-like), or 'log2' (already logged)."""
    v = X.values.ravel()
    rng = np.random.default_rng(0)
    v = v[rng.choice(len(v), min(len(v), 200_000), replace=False)]
    if v.min() < 0:
        return "log2"
    intlike = np.mean(np.isclose(v, np.round(v))) > 0.99
    if intlike and v.max() > 100:
        return "counts"
    if v.max() > 100:
        return "linear"
    return "log2"


def filter_genes(X, kind, cfg):
    if kind == "counts":
        cpm = X.div(X.sum(axis=1), axis=0) * 1e6
        keep = (cpm >= cfg.min_cpm).mean(axis=0) >= cfg.min_sample_frac
    elif kind == "linear":
        keep = (X >= 1).mean(axis=0) >= cfg.min_sample_frac
    else:
        keep = X.std(axis=0) > 1e-6
        if X.values.min() >= 0:
            keep &= X.mean(axis=0) > 0.5
    return X.loc[:, keep]


def to_log2(X, kind, totals=None):
    """Return log2-scale values. For counts: log2(CPM+1) using library sizes `totals`."""
    if kind == "counts":
        totals = X.sum(axis=1) if totals is None else totals
        return np.log2(X.div(totals, axis=0) * 1e6 + 1)
    if kind == "linear":
        return np.log2(X + 1)
    return X.copy()
