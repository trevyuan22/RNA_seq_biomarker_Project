"""Robust file reading: any file name, any common layout.

Handles: genes x samples or samples x genes; extra annotation columns (gene_id, gene_name, transcript...);
transcript-level rows (summed to genes); Ensembl version suffixes; csv/tsv/txt/xlsx, optionally gzipped;
sample sheets with arbitrary column names; groups inferred from sample names when no sheet is given."""
import gzip
import re
from pathlib import Path
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------- reading tables
def _inner_ext(path):
    n = Path(path).name.lower()
    n = n[:-3] if n.endswith(".gz") else n
    return Path(n).suffix


def sniff_sep(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        head = f.readline()
    counts = {"\t": head.count("\t"), ",": head.count(","), ";": head.count(";")}
    return max(counts, key=counts.get) if max(counts.values()) else ","


def read_table(path, nrows=None, **kw):
    if _inner_ext(path) in (".xlsx", ".xls"):
        return pd.read_excel(path, nrows=nrows, **kw)
    sep = sniff_sep(path)
    try:
        return pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False, **kw)
    except UnicodeDecodeError:
        return pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False, encoding="latin-1", **kw)


def count_lines(path, limit=5_000_000):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    n = 0
    with opener(path, "rb") as f:
        for n, _ in enumerate(f, 1):
            if n >= limit:
                break
    return n


def is_numeric(series):
    return pd.to_numeric(series, errors="coerce").notna().mean() >= 0.9


def split_columns(df):
    text, num = [], []
    for c in df.columns:
        (num if is_numeric(df[c]) else text).append(c)
    return text, num


_ID_PATTERNS = [r"^gene[_ .]?name$", r"^(gene[_ .]?)?symbol$", r"^gene$", r"^gene[_ .]?id$",
                r"^(id|name|feature|features|transcript|ensembl.*|sample|samples)$"]


def pick_id_column(text_cols):
    for pat in _ID_PATTERNS:
        for c in text_cols:
            if re.search(pat, str(c).strip().lower()):
                return c
    return text_cols[0]


def norm_key(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def clean_gene_ids(index):
    s = pd.Index(pd.Series(list(index)).astype(str).str.strip())
    if len(s) and s.str.startswith("ENSG").mean() > 0.5:           # strip Ensembl version suffix
        s = s.str.replace(r"\.\d+$", "", regex=True)
    return s


def _collapse(g):
    """genes x samples with possibly duplicated gene ids -> unique ids (sum for counts, mean otherwise)."""
    ok = ~pd.isna(g.index) & ~g.index.isin(["", "nan", "NaN", "None", "NA", "-"])
    g = g.loc[ok]
    if g.index.has_duplicates:
        v = g.values.ravel()
        v = v[np.random.default_rng(0).choice(len(v), min(len(v), 50000), replace=False)]
        countlike = v.min() >= 0 and np.mean(np.isclose(v, np.round(v))) > 0.99
        n = int(g.index.duplicated().sum())
        print(f"  {n} duplicated gene/transcript rows were combined ({'summed' if countlike else 'averaged'} per gene)")
        g = g.groupby(level=0).sum() if countlike else g.groupby(level=0).mean()
    return g


def _overlap(values, keyset):
    return sum(norm_key(v) in keyset for v in values)


def load_expression(path, sample_hint=None, gene_hint=None):
    """Return a samples x genes numeric DataFrame from (almost) any expression table.
    sample_hint / gene_hint: collections of names used only to decide orientation."""
    df = read_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    text, num = split_columns(df)
    if len(num) < 2:
        raise ValueError(f"{Path(path).name}: could not find a table of numeric expression values.")
    if text:
        id_col = pick_id_column(text)
        ids = df[id_col].astype(str).str.strip().values
    else:                                    # all numeric: first column holds the ids
        ids = df[num[0]].astype(str).values
        num = num[1:]
    vals = df[num].apply(pd.to_numeric, errors="coerce")
    vals.index = ids
    n_rows, n_cols = vals.shape
    transposed = None
    if sample_hint is not None:
        ks = {norm_key(s) for s in sample_hint}
        a, b = _overlap(ids, ks), _overlap(num, ks)
        if a != b:
            transposed = a > b
    elif gene_hint is not None:
        gs = {norm_key(g) for g in clean_gene_ids(list(gene_hint))}
        a, b = _overlap(clean_gene_ids(ids), gs), _overlap(clean_gene_ids(num), gs)
        if a != b:
            transposed = b > a
    if transposed is None:
        transposed = n_rows < n_cols and n_cols >= 1000
    if transposed:
        X = vals.copy()                                       # rows = samples, columns = genes
        X.columns = clean_gene_ids(num)
        X = _collapse(X.T).T
    else:
        g = vals.copy()
        g.index = clean_gene_ids(ids)
        X = _collapse(g).T                                    # samples x genes
    X.index = [str(s).strip() for s in X.index]
    if X.index.has_duplicates:
        X.index = pd.Index(pd.Series(X.index).where(~pd.Series(X.index).duplicated(), lambda s: s + "_dup"))
    if X.isna().any().any():
        print("  [warn] missing values found in the expression table; set to 0")
        X = X.fillna(0)
    return X.astype(float)


# ----------------------------------------------------------------------------- sample sheets / labels
_GROUP_NAME = re.compile(r"group|condition|label|class|status|genotype|disease|treatment|phenotype|diagnos|arm|cohort|type|cell", re.I)


def read_sample_sheet(path):
    df = read_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def sheet_info(sheet, sample_names):
    """Work out which column holds sample ids and which could hold the group."""
    ks = {norm_key(s) for s in sample_names}
    best, id_col = 0, None
    for c in sheet.columns:
        o = _overlap(sheet[c].astype(str), ks)
        if o > best:
            best, id_col = o, c
    cands = []
    for c in sheet.columns:
        if c == id_col:
            continue
        nu = sheet[c].nunique(dropna=True)
        if 2 <= nu < max(3, len(sheet)) and nu <= 30:
            cands.append(dict(name=c, n_unique=int(nu), examples=[str(v) for v in sheet[c].dropna().unique()[:5]]))
    rec = None
    if cands:
        named = [c for c in cands if _GROUP_NAME.search(c["name"])]
        rec = (named or sorted(cands, key=lambda c: c["n_unique"]))[0]["name"]
    return dict(id_col=id_col, matched=best, candidates=cands, recommended=rec)


def labels_from_sheet(sheet, sample_names, group_col=None):
    info = sheet_info(sheet, sample_names)
    if info["id_col"] is None:
        raise ValueError("None of the sample names in your data match the sample sheet.")
    group_col = group_col or info["recommended"]
    if group_col is None or group_col not in sheet.columns:
        raise ValueError("Could not tell which sample-sheet column holds the groups; choose it in the app "
                         "or pass --group-column.")
    key = {norm_key(v): g for v, g in zip(sheet[info["id_col"]].astype(str), sheet[group_col].astype(str).str.strip())}
    y = {s: key[norm_key(s)] for s in sample_names if norm_key(s) in key}
    return pd.Series(y, dtype=object), info, group_col


def _strip_rep(name):
    s, prev = str(name).strip(), None
    while s != prev:
        prev = s
        s = re.sub(r"[\s_\-.]+\d+$", "", s)
        s = re.sub(r"[\s_\-.]*(?:bio|biol|tech)?(?:replicate|repl|rep)$", "", s, flags=re.I)
        s = re.sub(r"([\s_\-.]+(?:rep|r)\d+|rep\d+)$", "", s, flags=re.I)
    return s or str(name)


def infer_groups_from_names(names):
    """Guess groups from sample names like 'ApoE3-1', 'ApoE3-2', 'ApoE4-1' ... Returns (Series|None, note)."""
    names = [str(n) for n in names]
    keys = {n: _strip_rep(n) for n in names}
    groups = sorted(set(keys.values()))
    counts = pd.Series(list(keys.values())).value_counts()
    if not (2 <= len(groups) <= 12) or counts.min() < 2 or len(groups) >= len(names):
        return None, "The sample names do not reveal the groups."
    toks = {g: re.split(r"[\s_\-.]+", g) for g in groups}
    while all(len(t) > 1 for t in toks.values()) and len({t[0] for t in toks.values()}) == 1:
        toks = {g: t[1:] for g, t in toks.items()}
    while all(len(t) > 1 for t in toks.values()) and len({t[-1] for t in toks.values()}) == 1:
        toks = {g: t[:-1] for g, t in toks.items()}
    nice = {g: " ".join(t) for g, t in toks.items()}
    return pd.Series({n: nice[keys[n]] for n in names}, dtype=object), "Groups were read from the sample names."


def labels_for(sample_names, labels=None, group_column=None):
    """labels: None | Series/dict (sample -> group) | path to a sample sheet. Returns Series over matched samples."""
    if labels is None:
        y, note = infer_groups_from_names(sample_names)
        if y is None:
            raise ValueError("No sample sheet was provided and the sample names do not reveal the groups. "
                             "Add a sample sheet, or assign groups in the app.")
        print(" ", note)
        return y
    if isinstance(labels, (dict, pd.Series)):
        s = pd.Series(labels, dtype=object)
        key = {norm_key(k): str(v) for k, v in s.items()}
        return pd.Series({n: key[norm_key(n)] for n in sample_names if norm_key(n) in key}, dtype=object)
    y, info, col = labels_from_sheet(read_sample_sheet(labels), sample_names, group_column)
    print(f"  groups taken from column '{col}' of {Path(labels).name}")
    return y


def sheet_hint(labels):
    """Names used only to decide the orientation of the expression table."""
    if labels is None:
        return None
    if isinstance(labels, (dict, pd.Series)):
        return list(pd.Series(labels).index)
    sh = read_sample_sheet(labels)
    return [v for c in sh.columns for v in sh[c].astype(str)]


def load_training_data(expression, labels=None, group_column=None):
    """Return X (samples x genes) and y (Series of group labels) with matching samples."""
    print(f"  reading {Path(expression).name}")
    X = load_expression(expression, sample_hint=sheet_hint(labels))
    y = labels_for(list(X.index), labels, group_column)
    common = [s for s in X.index if s in y.index]
    if not common:
        raise ValueError("No sample names match between the expression data and the labels/sample sheet.")
    if len(common) < len(X):
        print(f"  [warn] {len(X) - len(common)} samples without a group were left out")
    return X.loc[common], y.loc[common]
