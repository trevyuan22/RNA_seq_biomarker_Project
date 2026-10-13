"""Look inside a folder, work out what each file is (whatever it is called) and propose how to use it."""
import re
from pathlib import Path
from . import io_utils as io

TABLE_EXT = {".csv", ".tsv", ".txt", ".tab", ".xlsx", ".xls"}


def _gmt_role(name):
    n = name.lower()
    if re.search(r"tf|chea|trrust|encode|regulat|target", n):
        return "tf_gmt"
    if re.search(r"disea|gwas|omim|disgenet|phenotype|dis_", n):
        return "disease_gmt"
    return "pathway_gmt"


def analyze_table(path):
    """Return dict(kind, detail, samples) for one table-like file."""
    try:
        head = io.read_table(path, nrows=400)
    except Exception as e:                                  # unreadable / not a table
        return dict(kind="unknown", detail=f"could not be read as a table ({type(e).__name__})", samples=None)
    head.columns = [str(c).strip() for c in head.columns]
    text, num = io.split_columns(head)
    n_rows = max(io.count_lines(path) - 1, len(head)) if io._inner_ext(path) not in (".xlsx", ".xls") else len(head)
    low = [c.lower() for c in head.columns]
    has_chr = any(re.search(r"chrom|^chr$|seqname", c) for c in low)
    has_gene = any(re.search(r"gene|symbol", c) for c in low)
    if has_chr and has_gene and len(num) <= 2 and n_rows >= 500:
        return dict(kind="annotation", detail=f"gene-to-chromosome table ({n_rows:,} genes)", samples=None)
    if (len(num) >= 3 and n_rows >= 1000) or (len(num) >= 500 and n_rows >= 3):
        transposed = n_rows < len(num) and len(num) >= 500
        if transposed and text:
            idc = io.pick_id_column(text)
            samples = [str(v).strip() for v in io.read_table(path, usecols=[idc])[idc]]
        elif transposed:
            samples = [str(v) for v in head.iloc[:, 0]]
        else:
            samples = [str(c) for c in num]
        extra = " (rows are transcripts; they will be summed per gene)" if any(re.fullmatch(r"transcript.*", c) for c in low) and "gene_name" in low else ""
        shape = f"{len(samples)} samples x about {(len(num) if transposed else n_rows):,} genes"
        return dict(kind="expression", detail=f"expression table: {shape}{extra}", samples=samples)
    if n_rows <= 5000 and text:
        return dict(kind="sample_sheet", detail=f"sample sheet: {n_rows} rows; columns: {', '.join(map(str, head.columns[:6]))}", samples=None)
    return dict(kind="unknown", detail="not recognised as expression data, sample sheet or gene sets", samples=None)


def scan_folder(folder):
    folder = Path(folder)
    if not folder.is_dir():
        raise ValueError(f"Folder not found: {folder}")
    files = []
    for p in sorted(folder.iterdir()):
        if p.is_dir() or p.name.startswith((".", "~")) or p.name.lower().startswith("readme"):
            continue
        ext = io._inner_ext(p)
        if p.suffix.lower() == ".gmt":
            n = io.count_lines(p)
            info = dict(kind="gmt", detail=f"gene-set file ({n:,} sets)", samples=None)
        elif p.suffix.lower() == ".joblib":
            info = dict(kind="model", detail="trained model from an earlier classifier run", samples=None)
        elif ext in TABLE_EXT:
            info = analyze_table(p)
        else:
            continue
        sz = p.stat().st_size
        info.update(name=p.name, path=str(p), size_mb=round(sz / 1e6, 1),
                    size_label=(f"{sz / 1e6:.1f} MB" if sz >= 1e6 else f"{max(sz / 1e3, 0.1):.0f} KB"))
        files.append(info)
    _assign_roles(files)
    return files


def _assign_roles(files):
    expr = [f for f in files if f["kind"] == "expression"]
    for f in files:
        f["role"] = "ignore"
    pat = r"ribo|rpf|footprint"

    def is_ribo(f):
        names = f["samples"] or []
        return bool(re.search(pat, f["name"], re.I)) or (names and sum(bool(re.search(pat, n, re.I)) for n in names) >= 0.5 * len(names))

    def rna_like(f):
        names = f["samples"] or []
        return bool(re.search(r"count|rna|expr", f["name"], re.I)) or (names and sum(bool(re.search(r"rna", n, re.I)) for n in names) >= 0.5 * len(names))

    ribo = [f for f in expr if is_ribo(f)]
    rest = [f for f in expr if f not in ribo]
    if rest:
        rest.sort(key=lambda f: (not rna_like(f), -f["size_mb"]))
        rest[0]["role"] = "main"
    elif ribo:
        ribo[0]["role"] = "main"
        ribo = ribo[1:]
    for f in ribo[:1]:
        f["role"] = "ribo"
    main = next((f for f in files if f["role"] == "main"), None)
    sheets = [f for f in files if f["kind"] == "sample_sheet"]
    best, best_o = None, -1
    for s in sheets:
        try:
            sh = io.read_sample_sheet(s["path"])
            o = io.sheet_info(sh, main["samples"] or [])["matched"] if main else 0
        except Exception:
            o = -1
        s["matched"] = o
        if o > best_o:
            best, best_o = s, o
    if best is not None:
        best["role"] = "sheet"
    for f in files:
        if f["kind"] == "gmt":
            f["role"] = _gmt_role(f["name"])
        elif f["kind"] in ("annotation", "model"):
            f["role"] = f["kind"]


def group_info(expression_path, sheet_path=None, group_column=None, samples=None):
    """What the app needs to show the sample/group table."""
    samples = samples or next((f["samples"] for f in scan_folder(Path(expression_path).parent) if f["path"] == str(expression_path)), None)
    if not samples:
        samples = list(io.load_expression(expression_path).index)
    out = dict(samples=samples, labels={}, source="none", columns=[], column=None, note="")
    if sheet_path:
        sh = io.read_sample_sheet(sheet_path)
        info = io.sheet_info(sh, samples)
        out["columns"] = info["candidates"]
        try:
            y, _, col = io.labels_from_sheet(sh, samples, group_column)
            out.update(labels=y.to_dict(), source=f"sample sheet, column '{col}'", column=col)
            if len(y) < len(samples):
                out["note"] = f"{len(samples) - len(y)} of {len(samples)} samples were not found in the sheet; fill them in below."
            return out
        except ValueError as e:
            out["note"] = str(e)
    y, note = io.infer_groups_from_names(samples)
    if y is not None:
        out.update(labels=y.to_dict(), source="read from the sample names", note=note)
    elif not out["note"]:
        out["note"] = "The sample names do not reveal the groups. Type a group for each sample below."
    return out
