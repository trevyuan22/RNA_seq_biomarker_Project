#!/usr/bin/env python3
"""Download real RNA-seq data (TCGA, via UCSC Xena) and convert it to the pipeline's input format.

  python prepare_tcga.py                       # default 8 cancer types + normal tissue
  python prepare_tcga.py --cohorts BRCA LUAD COAD --per-class 80

Data: TCGA 'HiSeqV2' gene expression, log2(normalised RSEM + 1), gene symbols (hg19-era release).
Each cohort file is roughly 50-150 MB. If a download fails, see README section 'Getting real data'
for the manual route via https://xenabrowser.net/datapages/ .
"""
import argparse
import urllib.request
from pathlib import Path
import pandas as pd

URL = "https://tcga.xenahubs.net/download/TCGA.{c}.sampleMap/HiSeqV2.gz"
DEFAULT = ["BRCA", "LUAD", "COAD", "KIRC", "LIHC", "PRAD", "THCA", "STAD"]


def download(cohort, raw_dir):
    dest = raw_dir / f"{cohort}_HiSeqV2.gz"
    if dest.exists() and dest.stat().st_size > 1e6:
        return dest
    url = URL.format(c=cohort)
    print(f"  downloading {url}")
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:
        dest.unlink(missing_ok=True)
        raise SystemExit(f"\nCould not download {cohort}: {e}\nCheck your internet connection, or download the file "
                         f"manually (see README) and place it at {dest}")
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", nargs="+", default=DEFAULT)
    ap.add_argument("--per-class", type=int, default=100, help="max samples kept per disease")
    ap.add_argument("--out", default="input/tcga")
    a = ap.parse_args()
    out = Path(a.out)
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    mats, labels, normals = [], [], []
    for c in a.cohorts:
        df = pd.read_csv(download(c, raw), sep="\t", index_col=0)
        code = pd.Series([s.split("-")[3][:2] if len(s.split("-")) > 3 else "" for s in df.columns], index=df.columns)
        tumor = code[code == "01"].index[: a.per_class]          # primary solid tumour
        norm = code[code == "11"].index                           # solid tissue normal
        mats.append(df[tumor]); labels += [(s, c) for s in tumor]
        if len(norm):
            normals.append(df[norm])
        print(f"  {c}: {len(tumor)} tumour, {len(norm)} normal samples")
    if normals:
        n = pd.concat(normals, axis=1, join="inner").iloc[:, : a.per_class]
        mats.append(n); labels += [(s, "Normal") for s in n.columns]
    expr = pd.concat(mats, axis=1, join="inner")
    expr = expr.loc[~expr.index.duplicated()]
    pd.DataFrame(labels, columns=["sample", "label"]).to_csv(out / "labels.csv", index=False)
    expr.to_csv(out / "expression_log2.csv.gz")
    print(f"\nSaved {expr.shape[1]} samples x {expr.shape[0]} genes to {out}/")
    print(f"Next: open the app (python run.py), choose 'Biomarker discovery and disease classifier', and point it at {out}/")


if __name__ == "__main__":
    main()
