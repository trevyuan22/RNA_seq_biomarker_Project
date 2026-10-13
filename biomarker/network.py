"""Protein-protein interaction network among differentially expressed genes (STRING database).
Hub genes (many interaction partners) are good candidates for key regulators."""
import urllib.parse
import urllib.request
from io import StringIO
import pandas as pd


def fetch_string(genes, species=9606, score=700, timeout=120):
    data = urllib.parse.urlencode({"identifiers": "\r".join(genes), "species": species,
                                   "required_score": score, "caller_identity": "rnaseq_biomarker"}).encode()
    with urllib.request.urlopen(urllib.request.Request("https://string-db.org/api/tsv/network", data=data), timeout=timeout) as r:
        return r.read().decode("utf-8")


def parse_string(text):
    df = pd.read_csv(StringIO(text), sep="\t")
    df = df.rename(columns={"preferredName_A": "a", "preferredName_B": "b"})
    df = df[["a", "b", "score"]].copy()
    df["a"], df["b"] = df.a.astype(str).str.upper(), df.b.astype(str).str.upper()
    key = df.apply(lambda r: tuple(sorted((r.a, r.b))), axis=1)
    return df.loc[~key.duplicated() & (df.a != df.b)].reset_index(drop=True)


def string_network(genes, score=700):
    """Returns (edges DataFrame, node table gene/degree/betweenness)."""
    edges = parse_string(fetch_string(list(genes), score=score))
    deg = pd.concat([edges.a, edges.b]).value_counts()
    nodes = pd.DataFrame({"gene": deg.index, "degree": deg.values})
    try:
        import networkx as nx
        g = nx.from_pandas_edgelist(edges, "a", "b")
        bc = nx.betweenness_centrality(g)
        nodes["betweenness"] = nodes.gene.map(bc)
    except ImportError:
        nodes["betweenness"] = float("nan")
    return edges, nodes
