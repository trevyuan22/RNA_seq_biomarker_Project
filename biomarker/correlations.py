"""Correlations between factors, genes, modules and pathways."""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA


def spearman_matrix(F):
    F = F.loc[:, F.nunique() > 1].astype(float)
    cols = list(F.columns)
    n = len(cols)
    r, p = np.eye(n), np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            rho, pv = stats.spearmanr(F.iloc[:, i], F.iloc[:, j])
            r[i, j] = r[j, i] = rho
            p[i, j] = p[j, i] = pv
    return pd.DataFrame(r, cols, cols), pd.DataFrame(p, cols, cols)


def top_pairs(r, p, min_abs=0.8, skip_pc_pc=True):
    rows = []
    cols = list(r.columns)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = cols[i], cols[j]
            if skip_pc_pc and a.startswith("PC") and b.startswith("PC"):
                continue
            if abs(r.iloc[i, j]) >= min_abs:
                rows.append(dict(factor_a=a, factor_b=b, rho=r.iloc[i, j], pvalue=p.iloc[i, j]))
    return pd.DataFrame(rows, columns=["factor_a", "factor_b", "rho", "pvalue"]).sort_values("rho", key=abs, ascending=False)


def confounders(r, p, group_col, min_abs):
    """Factors that move together with the group (technical/biological confounders)."""
    if group_col not in r.columns:
        return pd.DataFrame()
    d = pd.DataFrame({"rho_with_group": r[group_col], "pvalue": p[group_col]}).drop(index=group_col)
    d = d[~d.index.str.startswith("PC")]
    return d[d.rho_with_group.abs() >= min_abs].sort_values("rho_with_group", key=abs, ascending=False)


def gene_factor_corr(Xl, genes, F):
    genes = [g for g in genes if g in Xl.columns and Xl[g].std() > 0]
    F = F.loc[:, F.nunique() > 1]
    big = pd.concat([Xl[genes], F], axis=1).corr(method="spearman")
    return big.loc[genes, F.columns]


def coexpression_modules(Xl, genes, y, case, n_modules=None):
    """Hierarchical clustering of genes by correlation across samples (WGCNA-lite).
    Returns gene table (module, kME), module eigengenes, module-vs-group correlation, gene order."""
    genes = [g for g in genes if g in Xl.columns and Xl[g].std() > 0]
    if len(genes) < 6:
        return None
    Z = Xl[genes]
    Z = (Z - Z.mean()) / (Z.std(ddof=0) + 1e-9)
    C = Z.corr().fillna(0).values
    D = np.clip(1 - C, 0, 2)
    np.fill_diagonal(D, 0)
    link = linkage(squareform((D + D.T) / 2, checks=False), "average")
    k = n_modules or int(min(4, max(2, len(genes) // 15)))
    mod = fcluster(link, k, criterion="maxclust")
    eig, rows = {}, []
    grp = (y == case).astype(int).values
    for m in sorted(set(mod)):
        mg = [g for g, mm in zip(genes, mod) if mm == m]
        pc = PCA(1).fit_transform(Z[mg].values)[:, 0]
        if np.corrcoef(pc, Z[mg].mean(axis=1))[0, 1] < 0:
            pc = -pc
        eig[f"M{m}"] = pc
        for g in mg:
            rows.append(dict(gene=g, module=f"M{m}", kME=float(np.corrcoef(Z[g], pc)[0, 1])))
    E = pd.DataFrame(eig, index=Xl.index)
    mg = pd.DataFrame({"module": E.columns,
                       "n_genes": [sum(1 for r_ in rows if r_["module"] == c) for c in E.columns],
                       "rho_with_group": [stats.spearmanr(E[c], grp)[0] for c in E.columns]})
    order = [genes[i] for i in leaves_list(link)]
    return dict(genes=pd.DataFrame(rows), eigengenes=E, module_group=mg, order=order, corr=pd.DataFrame(C, genes, genes))


def pathway_scores(Xl, sets, terms):
    """Per-sample activity of gene sets: mean z-score of member genes."""
    up = {c.upper(): c for c in Xl.columns}
    Z = (Xl - Xl.mean()) / (Xl.std(ddof=0) + 1e-9)
    out = {}
    for t in terms:
        have = [up[g] for g in sets.get(t, ()) if g in up]
        if len(have) >= 3:
            out[t[:55]] = Z[have].mean(axis=1)
    return pd.DataFrame(out)
