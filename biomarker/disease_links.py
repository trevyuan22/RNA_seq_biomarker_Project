"""Link differentially expressed genes to known diseases (DisGeNET, GWAS Catalog, OMIM, ...)."""
import pandas as pd
from .enrichment import ora


def disease_enrichment(genes, universe, collections):
    out = [ora([g.upper() for g in genes], [u.upper() for u in universe], sets, library=lib)
           for lib, sets in collections.items()]
    out = [o for o in out if len(o)]
    if not out:
        return pd.DataFrame()
    return pd.concat(out).sort_values("padj").reset_index(drop=True)


def gene_disease_table(genes, collections, keyword=None, max_terms=3):
    """For each gene: how many disease terms list it, and which terms match the keyword."""
    gset = {g.upper() for g in genes}
    terms = {g: [] for g in gset}
    for lib, sets in collections.items():
        for term, members in sets.items():
            for g in members & gset:
                terms[g].append(term)
    kw = keyword.lower() if keyword else None
    rows = []
    for g, ts in terms.items():
        hits = [t for t in ts if kw and kw in t.lower()]
        rows.append(dict(gene=g, n_disease_terms=len(ts), keyword_hits=len(hits),
                         example_terms="; ".join((hits or ts)[:max_terms])))
    return pd.DataFrame(rows).sort_values(["keyword_hits", "n_disease_terms"], ascending=False).reset_index(drop=True)
