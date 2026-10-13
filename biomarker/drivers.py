"""Candidate-driver prioritisation: combine every line of evidence into one transparent table.
This ranks genes worth testing in the lab. It does NOT prove a gene causes disease."""
import numpy as np
import pandas as pd


def candidate_drivers(de, pathway_n=None, tf_map=None, ppi_degree=None, modules=None, disease=None,
                      ribo_ok=None, chrom=None, flagged_chroms=(), top=40):
    sig = de[de.significant].copy()
    if sig.empty:
        return pd.DataFrame()
    sig["gene_u"] = sig.gene.astype(str).str.upper()
    pathway_n, tf_map, ppi_degree = pathway_n or {}, tf_map or {}, ppi_degree or {}
    modules, disease, ribo_ok, chrom = modules or {}, disease or {}, ribo_ok or {}, chrom or {}
    hub_cut = max(3, np.percentile(list(ppi_degree.values()), 75)) if ppi_degree else np.inf
    rows = []
    for _, r in sig.iterrows():
        g = r.gene_u
        n_path = pathway_n.get(g, 0)
        tfs = tf_map.get(g, [])
        deg = ppi_degree.get(g, 0)
        kme = modules.get(g, np.nan)
        d_n, d_kw = disease.get(g, (0, 0))
        ribo = ribo_ok.get(g)
        ch = chrom.get(g, "")
        flags = [n_path > 0, len(tfs) > 0, deg >= hub_cut, (not np.isnan(kme)) and kme >= 0.8, d_n > 0, bool(ribo)]
        de_score = min(-np.log10(max(r.padj, 1e-30)), 30) / 30 * min(abs(r.log2FC), 5) / 5
        caution = []
        if ch and ch in flagged_chroms:
            caution.append(f"on chromosome {ch}, which is over-represented among DE genes (check dosage/sex confounding)")
        if max(r.mean_case, r.mean_ref) < 2:
            caution.append("low expression: fold-change unreliable")
        rows.append(dict(gene=r.gene, direction=r.direction, log2FC=r.log2FC, padj=r.padj, chromosome=ch,
                         pathways=n_path, regulators=", ".join(tfs[:4]), ppi_degree=deg,
                         coexpr_kME=kme, disease_terms=d_n, keyword_hits=d_kw,
                         translation_supported=("yes" if ribo else ("no" if ribo is not None else "")),
                         evidence_lines=int(sum(flags)) + (1 if d_kw else 0),
                         priority_score=round(de_score + sum(flags) + (1 if d_kw else 0), 3),
                         caution="; ".join(caution)))
    out = pd.DataFrame(rows).sort_values(["priority_score", "padj"], ascending=[False, True])
    return out.head(top).reset_index(drop=True)
