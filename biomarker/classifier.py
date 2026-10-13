"""Disease classifier: leakage-free model selection + held-out evaluation + explanations."""
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, classification_report,
                             confusion_matrix, f1_score, roc_auc_score)
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def build_candidates(cfg, n_genes):
    k = min(cfg.n_features, n_genes)

    def pipe(model):
        # Feature selection lives INSIDE the pipeline so it is re-done within each CV fold (no leakage).
        return Pipeline([("var", VarianceThreshold(1e-8)), ("sel", SelectKBest(f_classif, k=k)),
                         ("scale", StandardScaler()), ("clf", model)])
    return {
        "Logistic regression": pipe(LogisticRegression(C=0.05, max_iter=5000, class_weight="balanced")),
        "Random forest": pipe(RandomForestClassifier(n_estimators=300, n_jobs=-1, class_weight="balanced",
                                                      random_state=cfg.random_state)),
    }


def train_and_evaluate(X, y, cfg):
    counts = y.value_counts()
    if counts.min() < cfg.min_class_size:
        raise ValueError(f"Class '{counts.idxmin()}' has only {counts.min()} samples; need at least "
                         f"{cfg.min_class_size} per disease for a trustworthy classifier.")
    if counts.min() < cfg.warn_class_size:
        print(f"  [warn] smallest class has {counts.min()} samples - estimates will be noisy")
    classes = sorted(y.unique())
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=cfg.test_size, stratify=y,
                                          random_state=cfg.random_state)
    folds = max(2, min(cfg.cv_folds, int(ytr.value_counts().min())))
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=cfg.random_state)
    cands = build_candidates(cfg, X.shape[1])
    cv = {}
    for name, p in cands.items():
        s = cross_val_score(p, Xtr, ytr, cv=skf, scoring="balanced_accuracy")
        cv[name] = (float(s.mean()), float(s.std()))
        print(f"  CV balanced accuracy  {name:<20s} {s.mean():.3f} +/- {s.std():.3f}")
    best = max(cv, key=lambda k: cv[k][0])
    model = clone(cands[best]).fit(Xtr, ytr)
    pred, proba = model.predict(Xte), model.predict_proba(Xte)
    try:
        auc = roc_auc_score(yte == classes[1], proba[:, 1]) if len(classes) == 2 else \
            roc_auc_score(yte, proba, multi_class="ovr", average="macro", labels=classes)
    except ValueError:
        auc = float("nan")
    rep = pd.DataFrame(classification_report(yte, pred, labels=classes, output_dict=True, zero_division=0)).T
    metrics = dict(best_model=best, cv=cv, n_train=len(ytr), n_test=len(yte),
                   accuracy=accuracy_score(yte, pred), balanced_accuracy=balanced_accuracy_score(yte, pred),
                   macro_f1=f1_score(yte, pred, average="macro"), auc=auc,
                   confusion=confusion_matrix(yte, pred, labels=classes), classes=classes,
                   per_class=rep.loc[classes, ["precision", "recall", "f1-score", "support"]])
    print(f"  Held-out test: accuracy {metrics['accuracy']:.3f}, balanced {metrics['balanced_accuracy']:.3f}, "
          f"macro-F1 {metrics['macro_f1']:.3f}")
    final = clone(cands[best]).fit(X, y)       # deployed model uses all data
    return final, metrics


def selected_genes(pipeline, genes):
    genes = np.asarray(genes)
    g = genes[pipeline.named_steps["var"].get_support()]
    return g[pipeline.named_steps["sel"].get_support()]


def gene_weights(pipeline, genes):
    """DataFrame (selected genes x classes): signed LR weights, or RF importances."""
    sel = selected_genes(pipeline, genes)
    clf = pipeline.named_steps["clf"]
    cls = list(clf.classes_)
    if hasattr(clf, "coef_"):
        w = clf.coef_
        if w.shape[0] == 1:
            w = np.vstack([-w[0], w[0]])
        return pd.DataFrame(w.T, index=sel, columns=cls)
    imp = clf.feature_importances_
    return pd.DataFrame(np.tile(imp[:, None], (1, len(cls))), index=sel, columns=cls)


def explain_sample(pipeline, genes, x_row, cls, top=5):
    """Top genes pushing one sample toward class `cls` (logistic regression only)."""
    clf = pipeline.named_steps["clf"]
    if not hasattr(clf, "coef_"):
        return []
    z = pipeline[:-1].transform(pd.DataFrame([x_row], columns=genes))[0]
    classes = list(clf.classes_)
    w = clf.coef_
    w = np.vstack([-w[0], w[0]]) if w.shape[0] == 1 else w
    contrib = w[classes.index(cls)] * z
    sel = selected_genes(pipeline, genes)
    idx = np.argsort(contrib)[::-1][:top]
    return [(sel[i], float(contrib[i])) for i in idx]


def build_panels(de_results, contrasts, weights, cfg):
    """Biomarker panel per contrast: significant up-regulated genes ranked by evidence x effect,
    annotated with whether (and how strongly) the classifier relies on them."""
    panels = {}
    cls_of = {c["name"]: c["cls"] for c in contrasts}
    for name, df in de_results.items():
        up = df[df.direction == "up"].copy()
        if up.empty:
            up = df[df.log2FC > 0].head(cfg.top_genes_per_class * 3).copy()
        up["score"] = -np.log10(up.padj.clip(lower=1e-300)) * up.log2FC
        cls = cls_of[name]
        up["in_classifier"] = up.gene.isin(weights.index)
        up["model_weight"] = up.gene.map(weights[cls]) if cls in weights.columns else np.nan
        panels[name] = (up.sort_values("score", ascending=False).head(cfg.top_genes_per_class)
                        [["gene", "log2FC", "padj", "mean_case", "mean_ref", "in_classifier", "model_weight"]]
                        .reset_index(drop=True))
    return panels
