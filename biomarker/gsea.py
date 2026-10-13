"""Pre-ranked gene-set enrichment analysis (GSEA-style), implemented from scratch so it needs no extra
packages. Uses ALL genes ranked by evidence, so it can detect coordinated small shifts that a
significance cut-off misses. Simplified: gene-set permutation null (binned by set size), BH-adjusted p."""
import numpy as np
import pandas as pd
from .de import bh


def rank_metric(de):
    """Signed -log10(p) (capped) with a tiny fold-change tie-breaker; index = upper-case gene symbol."""
    m = np.sign(de.log2FC) * np.clip(-np.log10(de.pvalue.clip(lower=1e-300)), 0, 20) + 1e-3 * de.log2FC
    s = pd.Series(m.values, index=de.gene.astype(str).str.upper())
    s = s[~s.index.duplicated()]
    return s.sort_values(ascending=False)


def _arrays(idx, r, N):
    k = len(idx)
    w = r[idx]
    W = w.sum()
    if W <= 0:
        return None
    cum = np.cumsum(w) / W
    miss = (idx - np.arange(k)) / (N - k)
    return cum - miss, (cum - w / W) - miss


def _es(idx, r, N):
    a = _arrays(idx, r, N)
    if a is None:
        return 0.0
    up, down = a
    mx, mn = up.max(), down.min()
    return mx if mx >= -mn else mn


def _leading_edge(idx, r, N, es):
    up, down = _arrays(idx, r, N)
    return idx[: int(np.argmax(up)) + 1] if es > 0 else idx[int(np.argmin(down)):]


def preranked_gsea(rank, collections, n_perm=1000, min_size=10, max_size=500, seed=0):
    genes = rank.index.to_numpy()
    N = len(genes)
    r = np.abs(rank.values)
    pos = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(seed)
    grid = np.unique(np.round(np.geomspace(min_size, max_size, 40)).astype(int))
    nulls = {}
    rows = []
    for lib, sets in collections.items():
        for term, members in sets.items():
            idx = np.array(sorted(pos[g] for g in members if g in pos))
            k = len(idx)
            if k < min_size or k > max_size or N - k < 2:
                continue
            es = _es(idx, r, N)
            kb = int(grid[np.argmin(np.abs(grid - k))])
            if kb not in nulls:
                nulls[kb] = np.array([_es(np.sort(rng.choice(N, kb, replace=False)), r, N) for _ in range(n_perm)])
            null = nulls[kb]
            same = null[null > 0] if es > 0 else null[null < 0]
            if len(same) == 0:
                continue
            p = (1 + np.sum(np.abs(same) >= abs(es))) / (1 + len(same))
            nes = es / np.abs(same).mean()
            led = _leading_edge(idx, r, N, es)
            rows.append(dict(library=lib, term=term, size=k, ES=es, NES=nes, pvalue=p,
                             leading_edge=", ".join(genes[led][:15])))
    if not rows:
        return pd.DataFrame(columns=["library", "term", "size", "ES", "NES", "pvalue", "padj", "leading_edge"])
    df = pd.DataFrame(rows)
    df["padj"] = bh(df.pvalue.values)
    return df.assign(_a=df.NES.abs()).sort_values(["padj", "_a"], ascending=[True, False]).drop(columns="_a").reset_index(drop=True)
