"""Pathway enrichment: offline over-representation analysis (your own GMT file) or
online Enrichr (via gseapy). Fails soft: returns None with a reason."""
import numpy as np
import pandas as pd
from scipy import stats
from .de import bh


def read_gmt(path):
    sets = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = set(g for g in parts[2:] if g)
    return sets


def ora(genes, universe, gmt, min_overlap=2, min_size=5, library="user GMT"):
    U, G = set(universe), set(genes) & set(universe)
    N, n = len(U), len(G)
    rows = []
    for term, members in gmt.items():
        S = set(members) & U
        K = len(S)
        k = len(G & S)
        if K < min_size or k < min_overlap:
            continue
        p = stats.hypergeom.sf(k - 1, N, K, n)
        rows.append(dict(term=term, library=library, overlap=f"{k}/{K}", n_overlap=k,
                         fold_enrichment=(k / n) / (K / N), pvalue=p,
                         genes=", ".join(sorted(G & S)[:25])))
    if not rows:
        return pd.DataFrame(columns=["term", "library", "overlap", "n_overlap", "fold_enrichment", "pvalue", "padj", "genes"])
    df = pd.DataFrame(rows)
    df["padj"] = bh(df.pvalue.values)
    return df.sort_values("padj").reset_index(drop=True)


def run_enrichment(de_results, universe, cfg, gmt_path=None, offline=False):
    """Over-representation analysis of up/down genes for every contrast.
    Returns ({contrast: {'up': df, 'down': df}} | None, note)."""
    from . import genesets
    cols, notes = genesets.load_collections(genesets.PATHWAYS, gmt_path, offline)
    if not cols:
        return None, (notes[0] if notes else "Pathway enrichment skipped.")
    U = [u.upper() for u in universe]
    out = {}
    for name, df in de_results.items():
        out[name] = {}
        for d in ("up", "down"):
            genes = [g.upper() for g in df[df.direction == d].gene.head(cfg.max_genes_for_enrichment)]
            if len(genes) < 5:
                out[name][d] = pd.DataFrame()
                continue
            frames = [f for f in (ora(genes, U, sets, library=lib) for lib, sets in cols.items()) if len(f)]
            out[name][d] = pd.concat(frames).sort_values("padj").reset_index(drop=True) if frames else pd.DataFrame()
    note = ("Over-representation analysis (hypergeometric test, BH-adjusted) against " + ", ".join(cols)
            + "; background = expressed genes.")
    return out, note + ((" " + " ".join(notes)) if notes else "")
