"""Turns a settings dictionary (from the app or the command line) into a pipeline run."""
import re
from datetime import datetime
from pathlib import Path
import pandas as pd
from . import io_utils as io, pipeline


def _slug(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "run"


def run_job(s, out_root="outputs"):
    analysis = s["analysis"]
    name = _slug(s.get("out_name") or f"{analysis}_{datetime.now():%Y%m%d_%H%M}")
    out = Path(out_root) / name
    offline = not s.get("online", True)
    if analysis == "predict":
        if not s.get("model") or not s.get("expression"):
            raise ValueError("Choose a trained model and a file of new samples.")
        return pipeline.run_prediction(s["model"], s["expression"], out)
    if not s.get("expression"):
        raise ValueError("No main expression file selected.")
    labels = s.get("labels") or None
    sheet = s.get("sample_sheet") or None
    gcol = s.get("group_column") or None
    if labels is None and sheet is None:
        labels = None                              # groups will be read from the sample names
    if analysis == "classify":
        return pipeline.run_training(s["expression"], labels if labels else sheet, out, control=s.get("control") or None,
                                     gmt=s.get("pathway_gmt"), use_llm=bool(s.get("ai")), group_column=gcol, offline=offline)
    # ---- compare two groups
    if not s.get("case") or not s.get("control") or s["case"] == s["control"]:
        raise ValueError("Choose two different groups to compare.")
    metadata = None
    if sheet:
        try:
            sh = io.read_sample_sheet(sheet)
            X_names = io.read_table(s["expression"], nrows=0)  # only to fail early on unreadable files
            samples = list((s.get("labels") or {}).keys()) or None
            info = io.sheet_info(sh, samples or [])
            if info["id_col"]:
                md = sh.set_index(info["id_col"])
                drop = [c for c in (gcol or info["recommended"],) if c in md.columns]
                metadata = md.drop(columns=drop)
        except Exception:
            metadata = None
    ribo_labels = None
    if s.get("ribo") and sheet:
        try:
            ribo_names = list(io.load_expression(s["ribo"]).index)
            if io.sheet_info(io.read_sample_sheet(sheet), ribo_names)["matched"] >= 2:
                ribo_labels = sheet
        except Exception:
            ribo_labels = None
    return pipeline.run_analysis(
        s["expression"], labels if labels else sheet, out, s["case"], s["control"], modules=s.get("modules"),
        ribo_expression=s.get("ribo"), ribo_labels=ribo_labels, group_column=gcol, gmt=s.get("pathway_gmt"),
        tf_gmt=s.get("tf_gmt"), disease_gmt=s.get("disease_gmt"), annotation=s.get("annotation"), metadata=metadata,
        keyword=(s.get("keyword") or None), offline=offline)


_CONTROL_WORDS = re.compile(r"ctrl|control|healthy|normal|wild|wt|untreated|vehicle|mock|e3|neg", re.I)


def auto_settings(folder, analysis, case=None, control=None, **kw):
    """Build run settings for a folder by auto-detecting every file (used by the command line)."""
    from . import autodetect
    files = autodetect.scan_folder(folder)
    print("Files found in", folder)
    roles = {}
    for f in files:
        print(f"  {f['name']:36s} -> {f['role']:12s} ({f['detail']})")
        if f["role"] != "ignore":
            roles.setdefault(f["role"], f["path"])
    s = dict(analysis=analysis, expression=kw.pop("expression", None) or roles.get("main"), ribo=kw.pop("ribo", None) or roles.get("ribo"),
             sample_sheet=kw.pop("sample_sheet", None) or roles.get("sheet"), pathway_gmt=roles.get("pathway_gmt"),
             tf_gmt=roles.get("tf_gmt"), disease_gmt=roles.get("disease_gmt"), annotation=roles.get("annotation"),
             model=kw.pop("model", None) or roles.get("model"))
    s.update({k: v for k, v in kw.items() if v is not None})
    if analysis != "predict":
        if not s["expression"]:
            raise ValueError("No expression data file was found in that folder.")
        gi = autodetect.group_info(s["expression"], s["sample_sheet"], s.get("group_column"))
        s["labels"], s["group_column"] = gi["labels"], gi["column"]
        if gi["note"]:
            print(" ", gi["note"])
        groups = sorted(set(gi["labels"].values()))
        print("  groups:", ", ".join(f"{g} ({list(gi['labels'].values()).count(g)})" for g in groups))
        if analysis == "compare":
            if not control and not case and len(groups) == 2:
                control = next((g for g in groups if _CONTROL_WORDS.search(g)), groups[0])
                case = next(g for g in groups if g != control)
                print(f"  assuming case = {case}, control = {control} (use --case/--control to change)")
            elif len(groups) == 2 and bool(case) != bool(control):
                other = next(g for g in groups if g != (case or control))
                case, control = (case, other) if case else (other, control)
            s["case"], s["control"] = case, control
        else:
            s["control"] = control
    return s
