"""Gene-set collections: download (and cache) Enrichr libraries, or read your own GMT files.
Everything is upper-cased so human gene symbols match regardless of capitalisation."""
import json
import urllib.parse
import urllib.request
from pathlib import Path
from .enrichment import read_gmt

ENRICHR = "https://maayanlab.cloud/Enrichr"
CACHE = Path("data/gene_sets")

PATHWAYS = (("MSigDB_Hallmark", "KEGG", "Reactome", "GO_Biological_Process"),
            ("MSigDB_Hallmark_2020", "KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"))
REGULATORS = (("ChEA", "TRRUST"), ("ChEA_2022", "TRRUST_Transcription_Factors_2019"))
DISEASES = (("DisGeNET", "GWAS_Catalog", "OMIM_Disease", "Jensen_DISEASES"),
            ("DisGeNET", "GWAS_Catalog_2023", "OMIM_Disease"))
_OTHER_SPECIES = ("mouse", "rat", "yeast", "fly", "worm", "zebrafish")


def _get(url, timeout=90):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def parse_enrichr_text(text):
    sets = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        genes = {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
        if genes:
            sets[parts[0].strip()] = genes
    return sets


def list_libraries():
    f = CACHE / "_libraries.json"
    if f.exists():
        return json.loads(f.read_text())
    CACHE.mkdir(parents=True, exist_ok=True)
    data = json.loads(_get(f"{ENRICHR}/datasetStatistics"))
    names = [d["libraryName"] for d in data["statistics"]]
    f.write_text(json.dumps(names))
    return names


def resolve(patterns, fallback):
    """Pick the newest human library matching each name pattern; fall back to known names."""
    try:
        names = [n for n in list_libraries() if not any(w in n.lower() for w in _OTHER_SPECIES)]
    except Exception:
        return list(fallback)
    out = []
    for p in patterns:
        hits = sorted(n for n in names if p.lower() in n.lower())
        if hits:
            out.append(hits[-1])
    return out or list(fallback)


def fetch_library(name):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{name}.txt"
    if not f.exists():
        url = f"{ENRICHR}/geneSetLibrary?mode=text&libraryName={urllib.parse.quote(name)}"
        f.write_bytes(_get(url))
    return parse_enrichr_text(f.read_text(encoding="utf-8"))


def load_collections(kind, gmt=None, offline=False):
    """kind: (patterns, fallback) tuple from this module. Returns ({library: {term: set}}, notes)."""
    cols, notes = {}, []
    if gmt:
        cols["user GMT"] = {t: {g.upper() for g in s} for t, s in read_gmt(gmt).items()}
    elif offline:
        notes.append("offline mode: no gene-set libraries downloaded (supply a .gmt file to run this step offline).")
    else:
        for n in resolve(*kind):
            try:
                cols[n] = fetch_library(n)
            except Exception as e:
                notes.append(f"{n}: could not download ({type(e).__name__}). Check your internet connection.")
    return cols, notes
