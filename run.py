#!/usr/bin/env python3
"""RNA-seq Biomarker Studio.

  python run.py                                  open the app (recommended)
  python run.py scan    --folder input           show how every file in a folder is understood
  python run.py analyze --folder input           compare two groups (comprehensive)
  python run.py train   --folder input           multi-disease biomarkers + classifier
  python run.py predict --folder input --model outputs/run/model.joblib
File names do not matter: the program works out what each file is.
"""
import argparse
import sys
from pathlib import Path


def main():
    from biomarker import catalog, autodetect, runner
    ap = argparse.ArgumentParser(description="RNA-seq Biomarker Studio")
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("app", help="open the app in your browser (default)")

    def common(p):
        p.add_argument("--folder", default="input", help="folder holding your files (any names)")
        p.add_argument("--expression"); p.add_argument("--sample-sheet"); p.add_argument("--group-column")
        p.add_argument("--out-name"); p.add_argument("--offline", action="store_true", help="never use the internet")

    sc = sub.add_parser("scan"); sc.add_argument("--folder", default="input")
    a_ = sub.add_parser("analyze", help="compare two groups"); common(a_)
    a_.add_argument("--case"); a_.add_argument("--control")
    a_.add_argument("--modules", help="comma list; default all. Options: " + ",".join(catalog.ALL))
    a_.add_argument("--ribo"); a_.add_argument("--keyword")
    t = sub.add_parser("train", help="multi-disease biomarkers + classifier"); common(t)
    t.add_argument("--control"); t.add_argument("--ai", action="store_true", help="add AI summary (needs ANTHROPIC_API_KEY)")
    p = sub.add_parser("predict", help="predict disease of new samples"); common(p); p.add_argument("--model")
    a = ap.parse_args()

    if a.cmd in (None, "app"):
        import app
        sys.argv = [sys.argv[0]]
        return app.main()
    try:
        if a.cmd == "scan":
            for f in autodetect.scan_folder(a.folder):
                print(f"{f['name']:38s} {f['kind']:13s} role={f['role']:12s} {f['detail']}")
            return
        kw = dict(expression=a.expression, sample_sheet=a.sample_sheet, group_column=a.group_column,
                  out_name=a.out_name, online=not a.offline)
        if a.cmd == "analyze":
            kw.update(case=a.case, control=a.control, ribo=a.ribo, keyword=a.keyword,
                      modules=[m.strip() for m in a.modules.split(",")] if a.modules else None)
            s = runner.auto_settings(a.folder, "compare", **kw)
        elif a.cmd == "train":
            kw.update(control=a.control, ai=a.ai)
            s = runner.auto_settings(a.folder, "classify", **kw)
        else:
            kw.update(model=a.model)
            s = runner.auto_settings(a.folder, "predict", **kw)
        rep = runner.run_job(s)
        import webbrowser
        webbrowser.open(Path(rep).resolve().as_uri())
    except (ValueError, FileNotFoundError) as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
