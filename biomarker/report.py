"""Self-contained HTML reports (open in any browser, no internet needed).
Design: editorial / lab-notebook feel - warm paper, serif headings, one warm and one cool accent,
sticky contents column, sortable tables, print-friendly."""
import html
import re
import numpy as np
import pandas as pd

CSS = """
:root{--paper:#f4f0e6;--card:#fffdf9;--ink:#1f2430;--muted:#6b6f7a;--rule:#e4dfd2;--up:#b5412b;--down:#2b6a8c;--gold:#b7791f;
--serif:"Iowan Old Style","Palatino Linotype",Palatino,"Book Antiqua",Georgia,serif;
--sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif}
*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:24px}
body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.62 var(--sans);-webkit-font-smoothing:antialiased}
.shell{display:grid;grid-template-columns:250px minmax(0,1fr);max-width:1320px;margin:0 auto}
nav.toc{position:sticky;top:0;align-self:start;height:100vh;padding:46px 20px 24px 34px;overflow:auto}
.brand{font:700 10.5px var(--sans);letter-spacing:.18em;text-transform:uppercase;color:var(--muted);margin-bottom:20px}
nav.toc a{display:flex;gap:11px;padding:6px 0 6px 13px;margin-left:-15px;color:var(--muted);text-decoration:none;font-size:13.5px;
border-left:2px solid transparent;line-height:1.35}
nav.toc a b{font:700 10.5px var(--sans);color:#aaa392;min-width:20px;padding-top:3px;letter-spacing:.06em}
nav.toc a:hover{color:var(--ink)}nav.toc a.on{color:var(--ink);border-left-color:var(--up)}nav.toc a.on b{color:var(--up)}
main{padding:46px 44px 100px 30px;min-width:0}
header.masthead{margin-bottom:38px;padding-bottom:26px;border-bottom:2px solid var(--ink)}
.kicker{font:700 11px var(--sans);letter-spacing:.18em;text-transform:uppercase;color:var(--up)}
h1{font:600 42px/1.1 var(--serif);margin:12px 0 12px;letter-spacing:-.012em;max-width:28ch}
.meta{color:var(--muted);font-size:14px}.meta b{color:var(--ink);font-weight:600}
section{margin:0 0 62px}
.sec-head{display:flex;align-items:baseline;gap:16px;border-bottom:1px solid var(--rule);padding-bottom:11px;margin-bottom:18px}
.sec-head .n{font:700 12.5px var(--sans);color:var(--up);letter-spacing:.1em}
h2{font:600 27px/1.2 var(--serif);margin:0;letter-spacing:-.005em}
h3{font:700 11px var(--sans);letter-spacing:.14em;text-transform:uppercase;color:var(--muted);margin:30px 0 9px}
.lede{font:italic 400 16.5px/1.62 var(--serif);color:#444a57;max-width:72ch;margin:0 0 20px}
p{max-width:76ch;margin:.6em 0}.muted{color:var(--muted)}
.stats{display:flex;flex-wrap:wrap;margin:6px 0 24px;border-top:1px solid var(--rule);border-bottom:1px solid var(--rule)}
.stat{padding:18px 30px 15px 0;margin-right:30px;border-right:1px solid var(--rule)}.stat:last-child{border-right:0;margin-right:0}
.stat b{display:block;font:600 36px/1 var(--serif);font-variant-numeric:lining-nums;margin-bottom:6px}
.stat span{font:700 10px var(--sans);letter-spacing:.14em;text-transform:uppercase;color:var(--muted)}
.stat.up b{color:var(--up)}.stat.down b{color:var(--down)}
.callout{border-left:3px solid var(--gold);background:#fbf3df;padding:11px 17px;margin:12px 0;border-radius:0 5px 5px 0;font-size:14px;max-width:86ch}
.callout.warn{border-color:var(--up);background:#fbece7}.callout.info{border-color:var(--down);background:#e8f0f4}
.callout .t{display:block;font:700 10.5px var(--sans);letter-spacing:.14em;text-transform:uppercase;margin-bottom:3px;color:var(--ink)}
figure{margin:18px 0}figure img{max-width:100%;height:auto;border:1px solid var(--rule);border-radius:6px;background:var(--card);display:block}
figcaption{font-size:12.5px;color:var(--muted);margin-top:8px;max-width:72ch;line-height:1.5}figcaption b{color:var(--ink)}
.figs{display:flex;flex-wrap:wrap;gap:20px;align-items:flex-start}.figs figure{flex:1 1 400px;min-width:0;margin:10px 0}
.tablewrap{overflow-x:auto;border:1px solid var(--rule);border-radius:6px;background:var(--card);margin:12px 0}
table.data{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}
.data th{font:700 10px var(--sans);letter-spacing:.1em;text-transform:uppercase;color:var(--muted);text-align:left;padding:11px 13px;
border-bottom:1.5px solid var(--ink);white-space:nowrap;cursor:pointer;user-select:none}
.data th:hover{color:var(--ink)}.data td{padding:8px 13px;border-bottom:1px solid var(--rule);vertical-align:top}
.data tr:last-child td{border-bottom:0}.data tbody tr:hover{background:#f8f4ea}.data .num{text-align:right}
.data td.id{font-weight:600;white-space:nowrap}.up{color:var(--up);font-weight:600}.down{color:var(--down);font-weight:600}
.dots{letter-spacing:2px;color:var(--up);white-space:nowrap}.dots i{color:#d9d3c3;font-style:normal}
.tabs{display:flex;gap:24px;border-bottom:1px solid var(--rule);margin:8px 0 4px;overflow-x:auto}
.tabs button{background:none;border:0;border-bottom:2px solid transparent;padding:9px 0 8px;font:600 13px var(--sans);color:var(--muted);cursor:pointer;white-space:nowrap}
.tabs button.on{color:var(--ink);border-bottom-color:var(--up)}.pane{display:none}.pane.on{display:block}
details{margin:12px 0;font-size:14px}summary{cursor:pointer;color:var(--muted);font-weight:600}
.verdict{background:var(--card);border:1px solid var(--rule);border-radius:8px;padding:16px 20px;margin:14px 0;max-width:820px}
.v-head{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap}.v-sample{font:600 13px var(--sans);color:var(--muted)}
.v-call{font:600 26px/1.1 var(--serif)}.v-p{margin-left:auto;font:600 22px var(--serif);color:var(--down)}.v-p.low{color:var(--gold)}
.bar{height:6px;background:var(--rule);border-radius:3px;margin:11px 0 9px;overflow:hidden}.bar i{display:block;height:100%;background:var(--down)}
.bar.low i{background:var(--gold)}.v-sub{font-size:13px;color:var(--muted)}
footer{margin-top:60px;padding-top:16px;border-top:1px solid var(--rule);font-size:12.5px;color:var(--muted)}
@media(max-width:980px){.shell{display:block}nav.toc{position:static;height:auto;display:flex;gap:6px 18px;flex-wrap:wrap;padding:20px 18px 0}
nav.toc .brand{width:100%;margin-bottom:6px}nav.toc a{border:0;margin:0;padding:3px 0}main{padding:24px 18px 70px}h1{font-size:32px}}
@page{margin:14mm}
@media print{body{background:#fff}.shell{display:block}nav.toc{display:none}main{padding:0}.tabs{display:none}
.pane{display:block!important;break-inside:auto}.pane::before{content:attr(data-label);display:block;font:700 11px var(--sans);letter-spacing:.12em;text-transform:uppercase;margin:14px 0 4px}
figure,.tablewrap,.callout,.verdict,.stats{break-inside:avoid}section{margin-bottom:34px}.tablewrap{overflow:visible}}
"""
JS = """
function tab(g,i){document.querySelectorAll('[data-g="'+g+'"]').forEach(function(e){e.classList.toggle('on',e.dataset.i==i)})}
document.querySelectorAll('table.data').forEach(function(t){t.querySelectorAll('th').forEach(function(th,ci){th.addEventListener('click',function(){
var b=t.tBodies[0],rows=Array.from(b.rows),asc=th.dataset.asc!=='1';th.dataset.asc=asc?'1':'0';
var v=function(r){var x=r.cells[ci].dataset.v||r.cells[ci].textContent.trim();var n=parseFloat(x);return isNaN(n)?x.toLowerCase():n};
rows.sort(function(a,c){var p=v(a),q=v(c);return (p>q?1:p<q?-1:0)*(asc?1:-1)});rows.forEach(function(r){b.appendChild(r)})})})});
var links=document.querySelectorAll('nav.toc a'),secs=Array.from(document.querySelectorAll('main section'));
if('IntersectionObserver' in window){var io=new IntersectionObserver(function(es){es.forEach(function(e){if(e.isIntersecting){
links.forEach(function(l){l.classList.toggle('on',l.getAttribute('href')==='#'+e.target.id)})}})},{rootMargin:'-15% 0px -70% 0px'});secs.forEach(function(s){io.observe(s)})}
"""


def e(x):
    return html.escape(str(x))


def _fmt(v, floats):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "&ndash;"
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "&ndash;"
    if isinstance(v, (int, np.integer)):
        return f"{v:,}"
    if isinstance(v, (float, np.floating)):
        return f"{v:.1e}" if (abs(v) < 1e-3 and v != 0) else f"{v:.{floats}f}"
    return e(v)


def table(df, floats=3):
    if df is None or len(df) == 0:
        return "<p class='muted'>No results.</p>"
    cols = list(df.columns)
    isnum = {c: pd.api.types.is_numeric_dtype(df[c]) and not pd.api.types.is_bool_dtype(df[c]) for c in cols}
    head = "".join(f"<th class='{'num' if isnum[c] else ''}'>{e(str(c).replace('_', ' '))}</th>" for c in cols)
    body = []
    for _, r in df.iterrows():
        tds = []
        for i, c in enumerate(cols):
            v = r[c]
            if c == "evidence_lines" and not pd.isna(v):
                n = int(v)
                tds.append(f"<td data-v='{n}'><span class='dots'>{'&#9679;' * n}<i>{'&#9679;' * max(0, 7 - n)}</i></span></td>")
            elif c in ("direction", "translation_class") and isinstance(v, str):
                k = "up" if "up" in v.lower() else "down" if "down" in v.lower() else ""
                tds.append(f"<td class='{k}'>{e(v)}</td>")
            else:
                cls = ("num " if isnum[c] else "") + ("id" if i == 0 and c in ("gene", "term", "sample", "factor") else "")
                raw = f" data-v='{v}'" if isnum[c] and not pd.isna(v) else ""
                tds.append(f"<td class='{cls.strip()}'{raw}>{_fmt(v, floats)}</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    return f"<div class='tablewrap'><table class='data'><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def img(b64, alt="", cap=None):
    if not b64:
        return ""
    c = f"<figcaption>{cap}</figcaption>" if cap else ""
    return f'<figure><img alt="{e(alt)}" src="data:image/png;base64,{b64}">{c}</figure>'


def figs(*items):
    items = [i for i in items if i]
    return f"<div class='figs'>{''.join(items)}</div>" if items else ""


def tabs(group, items):
    btn = "".join(f"<button class='{'on' if i == 0 else ''}' data-g='{group}' data-i='{i}' onclick=\"tab('{group}',{i})\">{e(l)}</button>"
                  for i, (l, _) in enumerate(items))
    panes = "".join(f"<div class='pane {'on' if i == 0 else ''}' data-g='{group}' data-i='{i}' data-label='{e(l)}'>{h}</div>"
                    for i, (l, h) in enumerate(items))
    return f"<div class='tabs'>{btn}</div>{panes}"


def stat(value, label, kind=""):
    return f"<div class='stat {kind}'><b>{e(value)}</b><span>{e(label)}</span></div>"


def stats(*items):
    return f"<div class='stats'>{''.join(items)}</div>"


def callout(text, kind="", title=None, raw=False):
    t = f"<span class='t'>{e(title)}</span>" if title else ""
    return f"<div class='callout {kind}'>{t}{text if raw else e(text)}</div>"


def skip(note):
    return f"<p class='muted'><i>{e(note)}</i></p>"


def section(sid, num, title, lede, *parts):
    body = "".join(p for p in parts if p)
    lede_html = f"<p class='lede'>{lede}</p>" if lede else ""
    return (f"<section id='{sid}' data-title='{e(title)}' data-num='{num:02d}'><div class='sec-head'><span class='n'>{num:02d}</span>"
            f"<h2>{e(title)}</h2></div>{lede_html}{body}</section>")


def page(title, kicker, subtitle_html, body):
    toc = "".join(f"<a href='#{i}'><b>{n}</b><span>{t}</span></a>"
                  for i, t, n in re.findall(r"<section id='([^']+)' data-title='([^']*)' data-num='([^']*)'", body))
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><title>{e(title)}</title>"
            "<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<style>{CSS}</style></head><body><div class='shell'><nav class='toc'><div class='brand'>Contents</div>{toc}</nav>"
            f"<main><header class='masthead'><div class='kicker'>{e(kicker)}</div><h1>{e(title)}</h1>"
            f"<div class='meta'>{subtitle_html}</div></header>{body}"
            "<footer>Research output, not a medical diagnosis. Statistical associations only; validate findings independently.</footer>"
            f"</main></div><script>{JS}</script></body></html>")


def _write(path, html_text):
    open(path, "w", encoding="utf-8").write(html_text)


def _enrich_tabs(gid, tabs_dict, imgs, cfg, label_up, label_down, cols=("library", "term", "overlap", "padj", "genes")):
    items = []
    for d, lab in (("up", label_up), ("down", label_down)):
        t = (tabs_dict or {}).get(d)
        if t is None or len(t) == 0:
            body = skip("Too few significant genes, or no enriched terms.")
        else:
            body = img(imgs.get(d), cap="<b>How to read:</b> longer bars mean stronger evidence (-log10 adjusted p). Bars right of the dashed line pass the 5% FDR threshold.") \
                + table(t.head(cfg.top_pathways)[list(cols)])
        items.append((lab, body))
    return tabs(gid, items)


# =============================================================================== comprehensive report
def full_report(path, c):
    cfg, case, ctrl, mods = c["cfg"], c["case"], c["control"], c["modules"]
    d = c["de"]
    up, down = int((d.direction == "up").sum()), int((d.direction == "down").sum())
    secs, n = [], [0]

    def nxt():
        n[0] += 1
        return n[0]

    # ---- overview
    strongest = lambda k: ", ".join(e(g) for g in d[d.direction == k].head(5).gene)  # noqa: E731
    plain = (f"<p>Compared with <b>{e(ctrl)}</b>, <b>{up}</b> genes are more active and <b>{down}</b> less active in <b>{e(case)}</b> "
             f"(adjusted p &lt; {cfg.padj_cutoff} and at least a {2 ** cfg.lfc_cutoff:g}-fold change).")
    if up:
        plain += f" The clearest increases: {strongest('up')}."
    if down:
        plain += f" The clearest decreases: {strongest('down')}."
    plain += "</p>"
    drv = c.get("drivers")
    st = [stat(f"{c['n_samples']}", "samples"), stat(f"{c['n_genes']:,}", "genes analysed"), stat(up, f"higher in {case}", "up"), stat(down, f"lower in {case}", "down")]
    if drv is not None and len(drv):
        st.append(stat(len(drv), "candidate genes ranked"))
    ov = [stats(*st), plain]
    if drv is not None and len(drv):
        ov.append("<p><b>Best-supported candidates:</b> " + ", ".join(e(g) for g in drv.head(5).gene) + ". See the candidate ranking for the evidence behind each.</p>")
    if c["warnings"]:
        ov.append("<h3>Check these before trusting the results</h3>" + "".join(callout(w, "warn") for w in c["warnings"]))
    else:
        ov.append(callout("None of the automatic checks raised a concern. That does not rule out hidden confounders.", "info", "Checks"))
    ov.append(callout("Every result here shows what changes <i>together with</i> the condition. Together they nominate candidate drivers to test "
                      "in the lab; none of them proves a cause.", "info", "Association is not causation", raw=True))
    if c["notes"]:
        ov.append("<details><summary>Steps that were skipped or could not run</summary>" + "".join(f"<p>&bull; {e(x)}</p>" for x in c["notes"]) + "</details>")
    secs.append(section("overview", nxt(), "Overview", f"{e(case)} compared with {e(ctrl)} ({e(c['group_counts'])}). Statistics: {e(c['de_method'])}.", *ov))

    # ---- quality
    q = c.get("qc")
    if "qc" in mods and q:
        secs.append(section("quality", nxt(), "Is the comparison clean?",
                            "Before trusting a difference between groups, check that it is not caused by something else: the sex of the samples, "
                            "their cell-type makeup, sequencing depth, or one unusual sample.",
                            callout(q["sex_call"], "info", "Sex check"),
                            figs(img(q["pca"], "PCA", "<b>How to read:</b> each dot is a sample; distance means difference. Same-group samples should sit together."),
                                 img(q["samplecorr"], "correlation", "<b>How to read:</b> how similar each pair of samples is overall. A sample paler than its group-mates may be an outlier.")),
                            figs(img(q["bars"], "QC bars", "<b>How to read:</b> depth and genes detected should be similar across samples and not split by group."),
                                 img(q["marker_img"], "markers", "<b>How to read:</b> cell-type marker strength per sample. If a marker differs between groups, the groups may differ in cell composition.")),
                            "<h3>Per-sample metrics</h3>" + table(q["table"].reset_index().rename(columns={"index": "sample"})),
                            "<h3>Cell-type marker panels (mean log2 CPM)</h3>" + table(q["markers"].reset_index().rename(columns={"index": "sample"}), 1)))
    elif "qc" in mods:
        secs.append(section("quality", nxt(), "Is the comparison clean?", "", skip("This step did not run.")))

    # ---- factors
    f = c.get("factors")
    if "factors" in mods and f:
        conf = f["confounders"]
        secs.append(section("factors", nxt(), "How the factors relate",
                            "Every measured factor per sample (group, depth, sex genes, cell-type markers, X-chromosome activity, principal components, "
                            "and sample-sheet columns) is correlated with every other. Factors that move with the group could be confounding it. "
                            "With few samples the correlations are coarse: read them as flags, not proof.",
                            img(f["corr_img"], "factor correlations", "<b>How to read:</b> red = rise together, blue = move oppositely (Spearman rho, -1 to +1). Look at the first row/column: the group."),
                            "<h3>Factors that move together with the group</h3>" + (table(conf.reset_index().rename(columns={"index": "factor"})) if len(conf) else "<p class='muted'>None above the threshold.</p>"),
                            "<h3>Strongest factor-to-factor correlations</h3>" + table(f["pairs"].head(25)),
                            img(f.get("gene_factor_img"), "genes vs factors", "<b>How to read:</b> does a top gene track a confounder (such as the sex score)? A strong colour means it does.")))

    # ---- DE
    secs.append(section("genes", nxt(), "What changed",
                        f"Genes meeting adjusted p &lt; {cfg.padj_cutoff} and |log2 fold-change| &ge; {cfg.lfc_cutoff}. Click a column heading to sort.",
                        figs(img(c["volcano"], "volcano", "<b>How to read:</b> right = higher in the case group, left = lower. Higher on the plot = stronger statistical evidence."),
                             img(c.get("de_heat"), "heatmap", "<b>How to read:</b> each row is a gene, each column a sample. Red = above that gene's average, blue = below. Samples should cluster by group.")),
                        table(d.head(30)[["gene", "log2FC", "pvalue", "padj", "mean_case", "mean_ref", "direction"]])))

    # ---- chromosome
    ch = c.get("chrom")
    if "chromosome" in mods:
        secs.append(section("chromosome", nxt(), "Chromosome check",
                            "If the changed genes cluster on one chromosome (especially X), the cause is often sex, X-inactivation or copy number, not disease biology.",
                            skip(c["chrom_note"]) if ch is None else (callout(f"{c['chrom_note']} {ch['coverage'] * 100:.0f}% of tested genes were matched.", "info")
                                                                    + img(ch["img"], "chromosomes", "<b>How to read:</b> each box is the spread of fold-changes for one chromosome. Red boxes hold more significant genes than chance.")
                                                                    + table(ch["table"]))))

    # ---- pathways
    pw = c.get("pathways")
    if {"pathways", "gsea"} & mods and pw is not None:
        parts = []
        if "pathways" in mods:
            parts.append("<h3>Over-representation</h3>" + (callout(pw["note"], "info") if pw.get("note") else "")
                         + _enrich_tabs("orapw", pw["ora"], pw["ora_img"], cfg, "Among higher genes", "Among lower genes"))
        if "gsea" in mods:
            g = pw["gsea"]
            parts.append("<h3>Ranked analysis (GSEA)</h3>" + (skip("No gene sets were available.") if g is None or len(g) == 0 else
                         img(pw["gsea_img"], "GSEA", "<b>How to read:</b> red = pathway higher in the case group, blue = lower. A star marks pathways passing the 25% FDR convention for exploratory GSEA.")
                         + table(g.head(20)[["library", "term", "size", "NES", "padj", "leading_edge"]])))
            if pw.get("score_img"):
                parts.append(img(pw["score_img"], "pathway activity", "<b>How to read:</b> activity of the top pathways in each sample. A clean split between groups supports the pathway result."))
        secs.append(section("pathways", nxt(), "Biological pathways",
                            "<b>Over-representation</b> asks whether significant genes land in known pathways more than chance. <b>GSEA</b> uses all genes ranked by evidence, "
                            "so it can see pathways that shift together even when no single gene passes a cut-off.", *parts))

    # ---- regulators
    rg = c.get("regulators")
    if "regulators" in mods and rg is not None:
        g = rg["gsea"]
        secs.append(section("regulators", nxt(), "Upstream regulators",
                            "Which transcription factors have their known target genes shifted? These are candidate controllers sitting <i>upstream</i> of what you see.",
                            (callout(rg["note"], "info") if rg.get("note") else ""),
                            _enrich_tabs("oratf", rg["ora"], rg["ora_img"], cfg, "Targets among higher genes", "Targets among lower genes"),
                            "<h3>Ranked analysis (GSEA)</h3>" + (skip("Not available.") if g is None or len(g) == 0 else
                                                                  img(rg["gsea_img"], "regulator GSEA") + table(g.head(15)[["library", "term", "size", "NES", "padj", "leading_edge"]]))))

    # ---- networks
    nw, md = c.get("network"), c.get("modules_res")
    if {"network", "coexpression"} & mods:
        parts = []
        if "network" in mods:
            parts.append("<h3>Protein-interaction network</h3>" + (skip(c["network_note"]) if nw is None else
                         img(nw["img"], "network", "<b>How to read:</b> dots are genes, lines are known protein interactions. Big dots are hubs: genes connected to many other changed genes.")
                         + table(nw["nodes"].head(15)) + f"<p class='muted'>{nw['n_edges']} interactions among {len(nw['nodes'])} connected genes.</p>"))
        if "coexpression" in mods:
            parts.append("<h3>Co-expression modules</h3>" + (skip(c["modules_note"] or "Not computed.") if md is None else
                         img(md["img"], "modules", "<b>How to read:</b> blocks of red along the diagonal are sets of genes that rise and fall together. With few samples they mostly mirror the two groups.")
                         + table(md["module_group"]) + "<h3>Module members (kME = how closely a gene follows its module)</h3>" + table(md["genes"].sort_values("kME", ascending=False).head(25))))
        secs.append(section("networks", nxt(), "Networks and modules",
                            "Hub genes in the protein network and groups of co-expressed genes point to central players. Co-expression is most informative with 12 or more samples.", *parts))

    # ---- disease
    ds = c.get("disease")
    if "disease" in mods:
        body = skip(c["disease_note"] or "Not available.")
        if ds is not None:
            kw = ds.get("keyword")
            body = ("<h3>Diseases over-represented among your genes</h3>"
                    + (table(ds["enrich"].head(15)[["library", "term", "overlap", "padj", "genes"]]) if len(ds["enrich"]) else "<p class='muted'>No disease term is significantly over-represented.</p>")
                    + "<h3>Gene-by-gene disease links" + (f" (highlighting &ldquo;{e(kw)}&rdquo;)" if kw else "") + "</h3>" + table(ds["genes"].head(30)))
        secs.append(section("disease", nxt(), "Known disease associations",
                            "Genes already tied to a disease in DisGeNET, the GWAS Catalog or OMIM support, but never prove, a disease role.", body))

    # ---- ribo
    rb = c.get("ribo")
    if "ribo" in mods:
        body = skip(c["ribo_note"] or "Not available.")
        if rb is not None:
            body = (f"<p>RNA and Ribo-seq fold-changes agree with Spearman r = <b>{rb['stats']['spearman']:.2f}</b> across {rb['stats']['n_genes']:,} genes; "
                    f"<b>{rb['stats']['n_te_sig']}</b> genes show a translation-specific change. Ribo-seq statistics: {e(rb['stats']['method'])}.</p>"
                    + img(rb["img"], "RNA vs Ribo", "<b>How to read:</b> genes on the dashed diagonal change equally in RNA and translation. Red/blue dots are translated more/less than their RNA suggests.")
                    + table(rb["table"].head(20)[["gene", "log2FC_rna", "log2FC_ribo", "delta_TE", "padj_TE", "translation_class"]]))
        secs.append(section("translation", nxt(), "Transcription vs translation",
                            "Ribosome profiling shows what is being translated into protein. Where it departs from the RNA change, the control is at the protein-synthesis level.", body))

    # ---- drivers
    if "drivers" in mods:
        secs.append(section("drivers", nxt(), "Candidate driver genes",
                            "Each significant gene is scored by how many independent lines of evidence support it: pathway membership, regulator targets, network hub, "
                            "co-expression hub, disease links (with a bonus for your keyword) and translation support. Dots show the number of supporting lines. "
                            "A shortlist for experiments, not a verdict.",
                            table(drv) if drv is not None and len(drv) else skip("No significant genes to rank.")))

    # ---- reading guide
    secs.append(section("guide", nxt(), "How to read all this",
                        "",
                        "<p><b>Trust order.</b> Consistent effects across replicates, then agreement between several analyses, then a single analysis.</p>"
                        "<p><b>Confounders first.</b> If sex, cell type or depth tracks the group, the result may reflect that instead of the condition.</p>"
                        "<p><b>Small samples.</b> P-values can look extreme when replicates agree; they do not protect against an effect shared by a whole group, such as a cell-line or batch difference.</p>"
                        "<p><b>Towards cause.</b> Replicate in independent data, confirm key genes by qPCR or Western blot, then test function by knock-down, knock-out or rescue.</p>"))
    _write(path, page("RNA-seq analysis report", "Two-group comparison",
                      f"<b>{e(case)}</b> vs <b>{e(ctrl)}</b> &middot; {e(c['subtitle'])}", "".join(secs)))


# =============================================================================== classifier report
def training_report(path, c):
    m = c["metrics"]
    secs = []
    st = [stat(c["n_samples"], "samples"), stat(len(c["classes"]), "disease groups"), stat(f"{c['n_genes']:,}", "genes"),
          stat(f"{m['accuracy'] * 100:.0f}%", "held-out accuracy"), stat(f"{m['macro_f1']:.2f}", "macro F1"),
          stat("n/a" if np.isnan(m["auc"]) else f"{m['auc']:.2f}", "ROC AUC")]
    ov = [stats(*st), f"<p>The best model was <b>{e(m['best_model'])}</b>. Differential expression: {e(c['de_method'])}. Input data type: {e(c['kind'])}.</p>",
          callout(c.get("control") and f"Each disease was compared with '{c['control']}'." or
                  "Each disease was compared with all other groups combined, which highlights disease-specific markers.", "info"),
          callout("Predictions come from statistical patterns in the training data. Samples from different labs, protocols or tissues can give misleading answers.", "warn", "Research use only")]
    if c.get("llm"):
        ov.append("<h3>Plain-language interpretation (AI-generated)</h3>" + callout("Written by an AI model from the numbers in this report. Verify before relying on it.", "info")
                  + "".join(f"<p>{e(p)}</p>" for p in c["llm"].split("\n") if p.strip()))
    secs.append(section("summary", 1, "Summary", "", *ov))
    secs.append(section("data", 2, "Data overview", "A map of all samples, and how many belong to each group.",
                        figs(img(c["pca"], "PCA", "<b>How to read:</b> samples far from their group may be outliers or batch effects."), table(c["class_counts"]))))
    items = [(n_, f"<p><b>{int((d_.direction == 'up').sum())}</b> higher, <b>{int((d_.direction == 'down').sum())}</b> lower.</p>"
              + img(c["volcanos"][n_], n_, "<b>How to read:</b> right = higher in this disease; higher on the plot = stronger evidence.")
              + table(d_.head(15)[["gene", "log2FC", "padj", "mean_case", "mean_ref", "direction"]])) for n_, d_ in c["de"].items()]
    secs.append(section("genes", 3, "Differentially expressed genes",
                        f"Adjusted p &lt; {c['cfg'].padj_cutoff} and |log2 fold-change| &ge; {c['cfg'].lfc_cutoff}.", tabs("de", items)))
    if c["enrich"] is None:
        secs.append(section("pathways", 4, "Pathway enrichment", "", skip(c["enrich_note"])))
    else:
        items = []
        for n_, r in c["enrich"].items():
            h = ""
            for k, lab in (("up", "Pathways among higher genes"), ("down", "Pathways among lower genes")):
                t = r.get(k)
                h += f"<h3>{lab}</h3>" + (skip("Too few significant genes or no enriched terms.") if t is None or len(t) == 0
                                          else img(c["enrich_plots"].get((n_, k))) + table(t.head(c["cfg"].top_pathways)[["term", "library", "overlap", "padj", "genes"]]))
            items.append((n_, h))
        secs.append(section("pathways", 4, "Pathway enrichment", e(c["enrich_note"]), tabs("en", items)))
    secs.append(section("markers", 5, "Biomarker panels",
                        "Top genes significantly higher in each disease, ranked by evidence times effect size. <i>in classifier</i> shows whether the trained model uses the gene; "
                        "a positive <i>model weight</i> means higher expression pushes the prediction toward that disease.",
                        img(c["heatmap"], "heatmap", "<b>How to read:</b> each column is a disease; red = that gene is high in the disease relative to the others."),
                        tabs("bm", [(n_, table(p_)) for n_, p_ in c["panels"].items()])))
    cv = pd.DataFrame([dict(model=k, balanced_accuracy=v[0], spread=v[1]) for k, v in m["cv"].items()])
    secs.append(section("performance", 6, "How well does the classifier work?",
                        f"Trained on {m['n_train']} samples and tested once on {m['n_test']} held-out samples that were never used for gene selection or model choice.",
                        figs(img(c["confusion"], "confusion matrix", "<b>How to read:</b> rows are the true disease, columns the prediction. A strong diagonal means correct calls."),
                             "<div>" + "<h3>Cross-validated model comparison</h3>" + table(cv) + "<h3>Per disease</h3>" + table(m["per_class"].reset_index().rename(columns={"index": "disease"})) + "</div>")))
    secs.append(section("methods", 7, "Methods and caveats", "",
                        "<p>Genes were filtered for low expression and converted to log2 counts-per-million. Gene selection and scaling happen inside cross-validation to prevent leakage. "
                        "The deployed model is refit on all samples after the held-out test was reported.</p>"
                        "<p>Top genes are candidates to validate, not proof. Always test on an independent cohort before trusting a biomarker.</p>"))
    _write(path, page("Biomarker discovery and classifier", "Multi-disease analysis", e(c["subtitle"]), "".join(secs)))


# =============================================================================== prediction report
def prediction_report(path, c):
    parts = [callout("Predictions are statistical estimates, not a diagnosis. Accuracy depends on how similar the new samples are to the training data.", "warn", "Research use only")]
    for w in c["warnings"]:
        parts.append(callout(w, "warn"))
    if c.get("accuracy") is not None:
        parts.append(callout(f"{c['accuracy'] * 100:.0f}% of the calls match the true labels you supplied.", "info", "Check against known labels"))
    for r in c["rows"]:
        low = "" if r["confident"] else " low"
        why = ", ".join(g for g, _ in r["why"]) or "n/a"
        parts.append(f"<div class='verdict'><div class='v-head'><span class='v-sample'>{e(r['sample'])}</span><span class='v-call'>{e(r['top1'])}</span>"
                     f"<span class='v-p{low}'>{r['p1'] * 100:.0f}%</span></div><div class='bar{low}'><i style='width:{r['p1'] * 100:.0f}%'></i></div>"
                     f"<div class='v-sub'>{'Confident call. ' if r['confident'] else '<b>Low confidence: treat as uncertain.</b> '}Next most likely: {e(r['top2'])} ({r['p2'] * 100:.0f}%), "
                     f"{e(r['top3'])} ({r['p3'] * 100:.0f}%).<br>Genes pushing toward {e(r['top1'])}: {e(why)}.</div></div>")
    secs = [section("calls", 1, "Predicted disease", "The most likely disease for each new sample, with how sure the model is.", *parts),
            section("overview", 2, "Probabilities", "", img(c["bars"], "probabilities", "<b>How to read:</b> each bar splits one sample's probability across the possible diseases. One dominant colour means a clear call.")),
            section("table", 3, "Full table", "", table(c["table"]))]
    _write(path, page("Disease prediction", "New samples", e(c["subtitle"]), "".join(secs)))
