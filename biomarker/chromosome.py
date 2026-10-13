"""Chromosome-level check: are differentially expressed genes concentrated on one chromosome
(e.g. X-chromosome dosage / X-inactivation differences between cell lines or donors)?"""
import urllib.parse
import urllib.request
from io import StringIO
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from .de import bh

CACHE = Path("data/annotation/gene_chromosomes.csv")
ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
_XML = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" '
        'header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" '
        'interface="default"><Attribute name="external_gene_name"/><Attribute name="chromosome_name"/></Dataset></Query>')


def parse_biomart(text):
    df = pd.read_csv(StringIO(text), sep="\t", header=None, names=["gene", "chromosome"]).dropna()
    df["gene"] = df.gene.astype(str).str.upper()
    df["chromosome"] = df.chromosome.astype(str)
    df = df[df.chromosome.isin(ORDER)]
    return df.drop_duplicates("gene").reset_index(drop=True)


def load_annotation(path=None, offline=False):
    """Returns (DataFrame gene,chromosome | None, note)."""
    if path:
        df = pd.read_csv(path)
        df.columns = [c.lower() for c in df.columns]
        df = df.rename(columns={df.columns[0]: "gene", df.columns[1]: "chromosome"})
        df["gene"] = df.gene.astype(str).str.upper()
        df["chromosome"] = df.chromosome.astype(str).str.replace("chr", "", case=False)
        return df[["gene", "chromosome"]], f"Gene locations from {path}."
    if CACHE.exists():
        return pd.read_csv(CACHE), "Gene locations from cached Ensembl BioMart download."
    if offline:
        return None, "Chromosome check skipped (offline). Supply --annotation genes.csv (columns: gene,chromosome)."
    last = None
    for host in ("https://www.ensembl.org", "https://useast.ensembl.org"):
        try:
            url = f"{host}/biomart/martservice?query={urllib.parse.quote(_XML)}"
            with urllib.request.urlopen(url, timeout=180) as r:
                df = parse_biomart(r.read().decode("utf-8"))
            if len(df) > 1000:
                CACHE.parent.mkdir(parents=True, exist_ok=True)
                df.to_csv(CACHE, index=False)
                return df, "Gene locations from Ensembl BioMart."
        except Exception as e:
            last = e
    return None, f"Chromosome check skipped: could not download gene locations ({type(last).__name__})."


def chromosome_summary(de, ann, cfg):
    d = de.assign(gene_u=de.gene.astype(str).str.upper()).merge(ann, left_on="gene_u", right_on="gene", suffixes=("", "_a"))
    coverage = len(d) / max(len(de), 1)
    tot_sig = int(d.significant.sum())
    rows = []
    for c in ORDER:
        s = d[d.chromosome == c]
        if len(s) < 10:
            continue
        other = d[d.chromosome != c]
        k, n = int(s.significant.sum()), len(s)
        p_enr = stats.hypergeom.sf(k - 1, len(d), tot_sig, n) if tot_sig else 1.0
        p_lfc = stats.mannwhitneyu(s.log2FC, other.log2FC).pvalue if n > 20 else np.nan
        rows.append(dict(chromosome=c, genes_tested=n, significant=k, pct_significant=100 * k / n,
                         up=int((s.direction == "up").sum()), down=int((s.direction == "down").sum()),
                         median_log2FC=s.log2FC.median(), p_enrichment=p_enr, p_shift=p_lfc))
    t = pd.DataFrame(rows)
    if len(t):
        t["padj_enrichment"] = bh(t.p_enrichment.values)
    return t, d, coverage


def mean_expression(Xl, ann, chrom):
    up = {c.upper(): c for c in Xl.columns}
    g = [up[x] for x in ann[ann.chromosome == chrom].gene if x in up]
    return Xl[g].mean(axis=1) if len(g) >= 20 else None
