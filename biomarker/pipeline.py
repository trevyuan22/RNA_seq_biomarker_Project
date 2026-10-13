"""Orchestration: ties the modules together into the three workflows
(compare two groups / classifier training / prediction)."""
import json
import re
from datetime import datetime
from pathlib import Path
import joblib
import numpy as np
import pandas as pd

from . import (catalog, chromosome, classifier, correlations, de, disease_links, drivers, enrichment, genesets,
               gsea, io_utils, llm, network, plots, preprocess, qc, report, ribo)
from .config import Config


def _step(msg):
    print(f"\n>> {msg}")


def _safe(label, fn, notes, default=None):
    try:
        return fn()
    except Exception as e:                    # keep going; the report lists what was skipped and why
        msg = f"{label} did not run: {type(e).__name__}: {e}"
        notes.append(msg)
        print(f"  [note] {msg}")
        return default


# ================================================================ classifier training
def run_training(expression, labels, out, control=None, gmt=None, use_llm=False, input_type="auto",
                 group_column=None, offline=False, cfg=None):
    cfg = cfg or Config()
    out = Path(out)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    _step("Loading data")
    X, y = io_utils.load_training_data(expression, labels, group_column)
    kind = preprocess.detect_scale(X) if input_type == "auto" else input_type
    print(f"  {X.shape[0]} samples x {X.shape[1]} genes; data type: {kind}")
    print("  samples per group:", {k: int(v) for k, v in y.value_counts().sort_index().items()})
    if y.value_counts().min() < cfg.min_class_size:
        raise ValueError(f"Group '{y.value_counts().idxmin()}' has only {y.value_counts().min()} samples; "
                         f"every disease needs at least {cfg.min_class_size}.")

    _step("Filtering genes and normalising")
    Xf = preprocess.filter_genes(X, kind, cfg)
    Xl = preprocess.to_log2(Xf, kind, X.sum(axis=1) if kind == "counts" else None)
    print(f"  kept {Xf.shape[1]} of {X.shape[1]} genes")
    medians = Xl.median()

    _step("Finding differentially expressed genes")
    contrasts = de.make_contrasts(y, control)
    de_res, de_method = de.run_de(Xl, y, contrasts, cfg, counts=Xf if kind == "counts" else None)
    for n, d in de_res.items():
        d.to_csv(out / "tables" / f"DE_{n}.csv", index=False)
        print(f"  {n}: {int((d.direction == 'up').sum())} up, {int((d.direction == 'down').sum())} down")

    _step("Pathway enrichment")
    enr, enr_note = enrichment.run_enrichment(de_res, Xl.columns, cfg, gmt, offline)
    print(" ", enr_note)
    if enr:
        rows = [t.assign(contrast=n, direction=d) for n, r in enr.items() for d, t in r.items() if len(t)]
        if rows:
            pd.concat(rows).to_csv(out / "tables" / "pathway_enrichment.csv", index=False)

    _step("Training and testing the disease classifier")
    model, metrics = classifier.train_and_evaluate(Xl, y, cfg)
    weights = classifier.gene_weights(model, Xl.columns)
    panels = classifier.build_panels(de_res, contrasts, weights, cfg)
    for n, p in panels.items():
        p.to_csv(out / "tables" / f"biomarkers_{n}.csv", index=False)

    _step("Drawing figures")
    ctx = dict(
        cfg=cfg, metrics=metrics, classes=metrics["classes"], kind=kind, de=de_res, de_method=de_method,
        n_samples=len(y), n_genes=Xl.shape[1], enrich=enr, enrich_note=enr_note, panels=panels,
        class_counts=y.value_counts().sort_index().rename_axis("group").reset_index(name="samples"),
        pca=plots.pca_plot(Xl, y), heatmap=plots.marker_heatmap(Xl, y, panels),
        volcanos={n: plots.volcano(d, n.replace("_", " "), cfg) for n, d in de_res.items()},
        enrich_plots={(n, d): plots.enrichment_bar(t, f"{n} ({d})", cfg.top_pathways)
                      for n, r in (enr or {}).items() for d, t in r.items() if len(t)},
        confusion=plots.confusion_plot(metrics["confusion"], metrics["classes"]),
        subtitle=f"Generated {datetime.now():%d %B %Y, %H:%M}",
        control=control)
    if use_llm:
        summary = dict(metrics=dict(accuracy=metrics["accuracy"], macro_f1=metrics["macro_f1"], auc=metrics["auc"],
                                    model=metrics["best_model"], cv=metrics["cv"]),
                       biomarkers={n: p.head(8).gene.tolist() for n, p in panels.items()},
                       pathways={n: {d: t.head(5).term.tolist() for d, t in r.items() if len(t)}
                                 for n, r in (enr or {}).items()}, n_samples=len(y), groups=metrics["classes"])
        ctx["llm"], ctx["llm_note"] = llm.interpret(summary)
        print(" ", ctx["llm_note"] or "AI interpretation added.")

    _step("Saving report and model")
    report.training_report(out / "report.html", ctx)
    joblib.dump(dict(model=model, genes=list(Xl.columns), classes=metrics["classes"], kind=kind,
                     medians=medians, cfg=cfg, metrics=metrics), out / "model.joblib")
    json.dump(dict(accuracy=metrics["accuracy"], balanced_accuracy=metrics["balanced_accuracy"],
                   macro_f1=metrics["macro_f1"], auc=metrics["auc"], best_model=metrics["best_model"],
                   de_method=de_method, n_samples=len(y), n_genes=int(Xl.shape[1])),
              open(out / "run_summary.json", "w"), indent=2, default=float)
    print(f"\nDONE. Report: {(out / 'report.html').resolve()}")
    return out / "report.html"


# ================================================================ prediction
def run_prediction(model_path, expression, out, truth=None, input_type="auto"):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    bundle = joblib.load(model_path)
    cfg, genes, classes = bundle["cfg"], bundle["genes"], bundle["classes"]
    warnings = []
    _step("Loading new samples")
    X = io_utils.load_expression(expression, gene_hint=set(genes))
    kind = preprocess.detect_scale(X) if input_type == "auto" else input_type
    if kind != bundle["kind"]:
        raise ValueError(f"The new data looks like '{kind}' values but the model was trained on '{bundle['kind']}' values. "
                         "Use data processed the same way (or override the data type).")
    Xl = preprocess.to_log2(X, kind, X.sum(axis=1) if kind == "counts" else None)
    present = float(np.mean([g in Xl.columns for g in genes]))
    print(f"  {present * 100:.1f}% of the model's genes were found in the new data")
    if present < 0.8:
        warnings.append(f"Only {present * 100:.0f}% of the genes the model needs were found. Missing genes were filled "
                        "with typical training values, so predictions may be unreliable. Check gene naming (symbols vs Ensembl IDs).")
    Xl = Xl.reindex(columns=genes).fillna(bundle["medians"])
    model = bundle["model"]
    proba = pd.DataFrame(model.predict_proba(Xl), index=Xl.index, columns=list(model.classes_))
    rows, table = [], []
    for s in Xl.index:
        order = proba.loc[s].sort_values(ascending=False)
        top = list(order.index) + [""] * 3
        pv = list(order.values) + [0.0] * 3
        confident = bool(pv[0] >= cfg.min_confidence and (pv[0] - pv[1]) >= 0.15)
        why = classifier.explain_sample(model, genes, Xl.loc[s].values, top[0])
        rows.append(dict(sample=s, top1=top[0], p1=pv[0], top2=top[1], p2=pv[1], top3=top[2], p3=pv[2],
                         confident=confident, why=why))
        table.append({"sample": s, "prediction": top[0], "confidence": pv[0], "confident_call": confident,
                      **{f"P({c})": proba.loc[s, c] for c in classes}})
    tdf = pd.DataFrame(table)
    tdf.to_csv(out / "predictions.csv", index=False)
    acc = None
    if truth:
        t = io_utils.labels_for(list(tdf["sample"]), truth)
        common = [s for s in tdf["sample"] if s in t.index]
        if common:
            pred = tdf.set_index("sample").prediction
            acc = float(np.mean([pred[s] == t[s] for s in common]))
            tdf["true_label"] = tdf["sample"].map(t)
    report.prediction_report(out / "report.html", dict(
        rows=rows, table=tdf, warnings=warnings, accuracy=acc, bars=plots.prob_bars(tdf, classes), classes=classes,
        subtitle=f"{len(tdf)} sample(s), generated {datetime.now():%d %B %Y, %H:%M}"))
    print(tdf[["sample", "prediction", "confidence"]].to_string(index=False))
    print(f"\nDONE. Report: {(out / 'report.html').resolve()}")
    return out / "report.html"


# ================================================================ two-group comprehensive analysis
def _ora_all(genes, universe, cols):
    frames = [enrichment.ora([g.upper() for g in genes], universe, sets, library=lib) for lib, sets in cols.items()]
    frames = [f for f in frames if len(f)]
    return pd.concat(frames).sort_values("padj").reset_index(drop=True) if frames else pd.DataFrame()


def _split(s):
    return [g.strip() for g in str(s).split(",") if g.strip()]


def _enrich_block(d, universe, cols, cfg, label, do_ora=True, do_gsea=True):
    """Over-representation (up/down) and/or pre-ranked GSEA over a dict of collections."""
    res = dict(ora={}, ora_img={}, gsea=None, gsea_img=None)
    if not cols:
        return res
    if do_ora:
        for dirn in ("up", "down"):
            genes = d[d.direction == dirn].gene.head(cfg.max_genes_for_enrichment).tolist()
            if len(genes) >= 5:
                t = _ora_all(genes, universe, cols)
                res["ora"][dirn] = t
                if len(t):
                    res["ora_img"][dirn] = plots.enrichment_bar(t, f"{label}: {dirn}-regulated genes", cfg.top_pathways)
    if do_gsea:
        res["gsea"] = gsea.preranked_gsea(gsea.rank_metric(d), cols, cfg.gsea_perms, cfg.gsea_min_size, cfg.gsea_max_size)
        res["gsea_img"] = plots.gsea_plot(res["gsea"], f"{label}: ranked analysis (GSEA)", fdr=cfg.gsea_fdr)
    return res


def run_analysis(expression, labels, out, case, control, modules=None, ribo_expression=None, ribo_labels=None,
                 group_column=None, gmt=None, tf_gmt=None, disease_gmt=None, annotation=None, metadata=None,
                 keyword=None, offline=False, input_type="auto", cfg=None):
    """Compare two groups. `modules` = ids from catalog.MODULES (None = all). Unselected modules are skipped
    and omitted from the report."""
    cfg = cfg or Config()
    mods = catalog.normalise(modules)
    out = Path(out)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    notes, warns = [], []
    T = lambda df, name: df.to_csv(out / "tables" / name, index=not isinstance(df.index, pd.RangeIndex))  # noqa: E731

    _step("Loading and normalising data")
    X, y = io_utils.load_training_data(expression, labels, group_column)
    for g in (case, control):
        if g not in set(y):
            raise ValueError(f"Group '{g}' not found. Groups present: {sorted(set(y))}")
    keep = y.isin([case, control])
    X, y = X[keep], y[keep]
    gc = {k: int(v) for k, v in y.value_counts().items()}
    if min(gc.values()) < 2:
        raise ValueError("Each group needs at least 2 samples.")
    kind = preprocess.detect_scale(X) if input_type == "auto" else input_type
    Xf = preprocess.filter_genes(X, kind, cfg)
    Xl = preprocess.to_log2(Xf, kind, X.sum(axis=1) if kind == "counts" else None)
    universe = [g.upper() for g in Xl.columns]
    print(f"  {X.shape[0]} samples; {Xl.shape[1]} of {X.shape[1]} genes kept; data type {kind}; groups {gc}")
    if min(gc.values()) < 4:
        warns.append("Very small groups (" + ", ".join(f"{k}: {v}" for k, v in gc.items()) + " samples). Treat every result as exploratory.")

    _step("Differential expression")
    name = f"{case}_vs_{control}"
    con = [dict(name=name, cls=case, case=(y == case).values, ref=(y == control).values)]
    de_res, method = de.run_de(Xl, y, con, cfg, counts=Xf if kind == "counts" else None)
    d = de_res[name]
    T(d, f"DE_{name}.csv")
    print(f"  {int((d.direction == 'up').sum())} up, {int((d.direction == 'down').sum())} down  [{method}]")

    qc_ctx, factors_ctx, F, extra = None, None, None, None
    qc_df = sex_df = markers = pc_df = pvar = None
    if "qc" in mods:
        _step("Quality control, sex check, cell-type markers")
        qc_df = qc.sample_qc(X, Xl, y, kind)
        sex_df, sex_call, _ = qc.sex_check(Xl)
        markers, _ = qc.marker_scores(Xl)
        pc_df, pvar = qc.pcs(Xl)
        lo = qc_df.corr_own_group.median() - 0.03
        for s_, v in qc_df.corr_own_group.items():
            if v < lo:
                warns.append(f"Sample {s_} correlates less with its own group ({v:.3f}) than the others; check it for quality problems.")
        T(qc_df, "sample_qc.csv")
        print(" ", sex_call)
        sc_ = Xl.T.corr()
        marker_img = (plots.rect_heatmap(markers.T, "Cell-type marker expression (log2 CPM)", cmap=plots.SEQ, vmin=0,
                                         vmax=max(float(markers.values.max()), 1), label="log2 CPM") if len(markers.columns) else None)
        qc_ctx = dict(table=qc_df, sex_call=sex_call, markers=markers, pca=plots.pca_labeled(pc_df, pvar, y),
                      samplecorr=plots.corr_heatmap(sc_, "Sample-to-sample correlation", vmin=float(sc_.values.min()) - .01, cmap=plots.SEQ),
                      bars=plots.qc_bars(qc_df), marker_img=marker_img)

    chrom, chrom_note, flagged, chrom_of = None, "", set(), {}
    if "chromosome" in mods:
        _step("Chromosome check")
        ann, chrom_note = chromosome.load_annotation(annotation, offline)
        if ann is not None:
            ctab, dmerged, cov = _safe("Chromosome analysis", lambda: chromosome.chromosome_summary(d, ann, cfg), notes, (None, None, 0))
            if ctab is not None and len(ctab):
                chrom = dict(table=ctab, coverage=cov, img=_safe("Chromosome figure", lambda: plots.chrom_plot(dmerged, ctab), notes))
                T(ctab, "chromosome_summary.csv")
                chrom_of = dict(zip(dmerged.gene_u, dmerged.chromosome))
                for _, r in ctab[(ctab.padj_enrichment < 0.05) & (ctab.significant >= 3)].iterrows():
                    flagged.add(r.chromosome)
                    warns.append(f"Significant genes are concentrated on chromosome {r.chromosome} ({int(r.significant)} of {int(r.genes_tested)} genes, "
                                 f"adjusted p = {r.padj_enrichment:.1e}). This can reflect sex, X-inactivation or copy-number differences between the samples rather than the condition.")
                xm = chromosome.mean_expression(Xl, ann, "X")
                if xm is not None:
                    extra = pd.DataFrame({"X-chromosome mean expression": xm})
        else:
            notes.append(chrom_note)

    if "factors" in mods:
        _step("Correlations between factors")
        if metadata is not None:
            md = metadata if isinstance(metadata, pd.DataFrame) else io_utils.read_sample_sheet(metadata)
            if not isinstance(metadata, pd.DataFrame):
                info = io_utils.sheet_info(md, Xl.index)
                md = md.set_index(info["id_col"]) if info["id_col"] else md.set_index(md.columns[0])
            md = md.copy()
            md.index = [str(i) for i in md.index]
            key = {io_utils.norm_key(i): i for i in md.index}
            md = md.loc[[key.get(io_utils.norm_key(s)) if io_utils.norm_key(s) in key else None for s in Xl.index]]
            md.index = Xl.index
            keepc = {}
            for c_ in md.columns:
                col = md[c_]
                if col.isna().all() or col.nunique() < 2 or col.nunique() == len(col) and not pd.api.types.is_numeric_dtype(col):
                    continue
                keepc[f"meta: {c_}"] = col.astype(float) if pd.api.types.is_numeric_dtype(col) else pd.Series(pd.factorize(col)[0], index=col.index).astype(float)
            if keepc:
                md = pd.DataFrame(keepc)
                extra = md if extra is None else pd.concat([extra, md], axis=1)
        F = qc.build_factors(y, case, qc_df, sex_df, markers, pc_df, extra)
        rmat, pmat = correlations.spearman_matrix(F)
        gcol = f"group ({case}=1)"
        conf = correlations.confounders(rmat, pmat, gcol, cfg.confound_rho)
        pairs = correlations.top_pairs(rmat, pmat, 0.8)
        T(rmat, "factor_correlations.csv")
        T(pairs.reset_index(drop=True), "factor_pairs.csv")
        if len(conf):
            warns.append("These measured factors move together with the group (|rho| >= %.1f): %s. Differences between groups may partly reflect them."
                         % (cfg.confound_rho, ", ".join(conf.index)))
        gf = _safe("Gene-factor correlation", lambda: correlations.gene_factor_corr(Xl, d.head(25).gene.tolist(), F), notes)
        factors_ctx = dict(corr_img=plots.corr_heatmap(rmat, "Spearman correlation between factors", annotate=len(rmat) <= 14),
                           pairs=pairs, confounders=conf,
                           gene_factor_img=plots.rect_heatmap(gf, "Top DE genes vs factors (Spearman rho)") if gf is not None and len(gf) else None)

    pathway_n, pw, path_cols = {}, None, {}
    if {"pathways", "gsea"} & mods:
        _step("Pathways")
        path_cols, n1 = genesets.load_collections(genesets.PATHWAYS, gmt, offline)
        notes += n1
        pw = _safe("Pathway analysis", lambda: _enrich_block(d, universe, path_cols, cfg, "Pathways", "pathways" in mods, "gsea" in mods), notes,
                   dict(ora={}, ora_img={}, gsea=None, gsea_img=None))
        pw["note"] = ("Background = your expressed genes. Libraries: " + ", ".join(path_cols)) if path_cols else (n1[0] if n1 else "")
        for dirn, t in pw["ora"].items():
            if len(t):
                for _, r in t[t.padj < cfg.ora_fdr].iterrows():
                    for g in _split(r.genes):
                        pathway_n[g] = pathway_n.get(g, 0) + 1
                T(t, f"pathway_ORA_{dirn}.csv")
        if pw["gsea"] is not None and len(pw["gsea"]):
            for _, r in pw["gsea"][pw["gsea"].padj < cfg.gsea_fdr].iterrows():
                for g in _split(r.leading_edge):
                    pathway_n[g] = pathway_n.get(g, 0) + 1
            T(pw["gsea"], "gsea_pathways.csv")
            flat = {t: s for sets in path_cols.values() for t, s in sets.items()}
            sc = _safe("Pathway scores", lambda: correlations.pathway_scores(Xl, flat, list(pw["gsea"].head(8).term)), notes)
            if sc is not None and len(sc.columns):
                pw["score_img"] = plots.rect_heatmap(sc.T.loc[:, list(y.sort_values(kind="stable").index)],
                                                     "Pathway activity per sample (z-score)", vmin=-1.5, vmax=1.5, label="z-score")

    rg, tf_map = None, {}
    if "regulators" in mods:
        _step("Upstream regulators")
        tf_cols, n2 = genesets.load_collections(genesets.REGULATORS, tf_gmt, offline)
        notes += n2
        rg = _safe("Regulator analysis", lambda: _enrich_block(d, universe, tf_cols, cfg, "Regulators"), notes,
                   dict(ora={}, ora_img={}, gsea=None, gsea_img=None))
        rg["note"] = ("Target gene sets from " + ", ".join(tf_cols)) if tf_cols else (n2[0] if n2 else "")
        sig_terms = []
        for dirn, t in rg["ora"].items():
            if len(t):
                sig_terms += list(t[t.padj < cfg.ora_fdr].term)
        if rg["gsea"] is not None and len(rg["gsea"]):
            sig_terms += list(rg["gsea"][rg["gsea"].padj < cfg.gsea_fdr].term)
            T(rg["gsea"], "gsea_regulators.csv")
        sig_set = set(d[d.significant].gene.str.upper())
        for lib, sets in tf_cols.items():
            for term in set(sig_terms):
                if term in sets:
                    tfname = re.split(r"[ _\-]", term)[0].upper()
                    for g in sets[term] & sig_set:
                        tf_map.setdefault(g, [])
                        if tfname not in tf_map[g]:
                            tf_map[g].append(tfname)

    network_res, net_note, ppi_degree = None, "", {}
    if "network" in mods:
        _step("Protein-interaction network")
        sig_sorted = d[d.significant].gene.head(cfg.network_genes).tolist()
        if offline:
            net_note = "Protein network skipped (internet use is switched off)."
        elif len(sig_sorted) < 5:
            net_note = "Protein network skipped: too few significant genes."
        else:
            r_ = _safe("Protein network", lambda: network.string_network(sig_sorted, cfg.string_score), notes)
            if r_ is not None and len(r_[0]):
                edges, nodes = r_
                network_res = dict(nodes=nodes, n_edges=len(edges), img=plots.network_plot(edges, nodes, d))
                ppi_degree = dict(zip(nodes.gene, nodes.degree))
                T(nodes, "network_nodes.csv")
                T(edges, "network_edges.csv")
            net_note = "No interactions were returned (or the STRING service could not be reached)."

    mod_res, mod_note, kme = None, "", {}
    if "coexpression" in mods:
        _step("Co-expression modules")
        mod_genes = d[d.significant].gene.head(cfg.module_genes).tolist() or d.head(cfg.module_genes).gene.tolist()
        if len(Xl) < 6:
            mod_note = "Co-expression modules need at least 6 samples."
        else:
            mm = _safe("Co-expression modules", lambda: correlations.coexpression_modules(Xl, mod_genes, y, case), notes)
            if mm is not None:
                heat = plots.corr_heatmap(mm["corr"].loc[mm["order"], mm["order"]], "Gene-to-gene correlation (clustered)", annotate=False)
                mod_res = dict(genes=mm["genes"], module_group=mm["module_group"], img=heat)
                kme = dict(zip(mm["genes"].gene.str.upper(), mm["genes"].kME))
                T(mm["genes"], "coexpression_modules.csv")
            else:
                mod_note = "Co-expression modules could not be computed."

    dis_res, dis_note, disease_map = None, "", {}
    if "disease" in mods:
        _step("Disease associations")
        dis_cols, n3 = genesets.load_collections(genesets.DISEASES, disease_gmt, offline)
        notes += n3
        dis_note = n3[0] if n3 else "Disease databases unavailable."
        if dis_cols:
            sigg = d[d.significant].gene.tolist()
            if len(sigg) >= 5:
                en = _safe("Disease enrichment", lambda: disease_links.disease_enrichment(sigg, universe, dis_cols), notes, pd.DataFrame())
                gt = disease_links.gene_disease_table(sigg, dis_cols, keyword)
                dis_res = dict(enrich=en, genes=gt, keyword=keyword)
                disease_map = {r.gene: (r.n_disease_terms, r.keyword_hits) for r in gt.itertuples()}
                T(gt, "disease_gene_links.csv")
                if len(en):
                    T(en, "disease_enrichment.csv")
            else:
                dis_note = "Too few significant genes to look up."

    ribo_res, ribo_note, ribo_ok = None, "", {}
    if "ribo" in mods:
        _step("Ribo-seq integration")
        if not ribo_expression:
            ribo_note = "No Ribo-seq file was selected."
        else:
            rr = _safe("Ribo-seq integration", lambda: ribo.run_ribo(ribo_expression, ribo_labels, case, control, d, cfg, input_type), notes)
            if rr is not None:
                mtab, rb_de, rst = rr
                ribo_res = dict(table=mtab, stats=rst, img=plots.rna_ribo_scatter(mtab, rst))
                rb_sig = rb_de.set_index(rb_de.gene.str.upper())
                d_idx = d.set_index(d.gene.str.upper())
                for g in set(d[d.significant].gene.str.upper()):
                    if g in rb_sig.index and g in d_idx.index:
                        r_, r0 = rb_sig.loc[g], d_idx.loc[g]
                        ribo_ok[g] = bool(r_.padj < cfg.padj_cutoff and np.sign(r_.log2FC) == np.sign(r0.log2FC))
                T(mtab, "ribo_vs_rna.csv")
            else:
                ribo_note = "Ribo-seq integration did not run (see the list of skipped steps)."

    drv = None
    if "drivers" in mods:
        _step("Ranking candidate genes")
        drv = drivers.candidate_drivers(d, pathway_n, tf_map, ppi_degree, kme, disease_map,
                                        ribo_ok if "ribo" in mods and ribo_expression else None,
                                        chrom_of, flagged, cfg.driver_top)
        if len(drv):
            T(drv.reset_index(drop=True), "candidate_drivers.csv")

    _step("Writing the report")
    ctx = dict(
        cfg=cfg, case=case, control=control, modules=mods, n_samples=len(y), n_genes=Xl.shape[1], de=d, de_method=method,
        group_counts=", ".join(f"{k}: {v}" for k, v in gc.items()), warnings=warns, notes=notes,
        subtitle=f"Generated {datetime.now():%d %B %Y, %H:%M}",
        qc=qc_ctx, factors=factors_ctx, volcano=plots.volcano(d, f"{case} vs {control}", cfg),
        de_heat=plots.de_heatmap(Xl, y, d[d.significant].gene.tolist() or d.head(30).gene.tolist()),
        chrom=chrom, chrom_note=chrom_note, pathways=pw, regulators=rg,
        network=network_res, network_note=net_note, modules_res=mod_res, modules_note=mod_note,
        disease=dis_res, disease_note=dis_note, ribo=ribo_res, ribo_note=ribo_note, drivers=drv, keyword=keyword)
    report.full_report(out / "report.html", ctx)
    json.dump(dict(case=case, control=control, groups=gc, de_method=method, modules=sorted(mods),
                   n_up=int((d.direction == "up").sum()), n_down=int((d.direction == "down").sum()),
                   warnings=warns, skipped=notes), open(out / "run_summary.json", "w"), indent=2)
    print(f"\nDONE. Report: {(out / 'report.html').resolve()}")
    if warns:
        print("\nREAD FIRST:")
        for w in warns:
            print("  -", w)
    return out / "report.html"
