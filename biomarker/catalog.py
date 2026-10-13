"""Single source of truth for what the tool can do, written in plain language.
The app and the command line both read this, so descriptions never drift apart."""

ANALYSES = [
    dict(id="compare", title="Compare two groups",
         tagline="Find what differs between a condition and a control, and what may be driving it.",
         use_when="You have a patient/treated/mutant group and a control group (for example ApoE4 vs ApoE3, or patients vs healthy). "
                  "Works with as few as 2-3 samples per group.",
         you_get=["Genes that go up or down", "Pathways and upstream regulators behind the changes",
                  "Checks for hidden confounders (sex, cell type, chromosome effects)",
                  "A ranked shortlist of candidate driver genes"]),
    dict(id="classify", title="Biomarker discovery and disease classifier",
         tagline="Find marker genes for each group and train a model that can recognise them in new samples.",
         use_when="You have many samples per group (at least 6, ideally 15 or more) and two or more groups, for example "
                  "disease vs healthy, or several diseases. Use it when the goal is to tell the groups apart. "
                  "With only a few samples per group, use \"Compare two groups\" instead.",
         you_get=["Marker genes and pathways for every disease", "A machine-learning model with honest accuracy figures",
                  "A saved model you can later use to predict new samples"]),
    dict(id="predict", title="Predict the disease of new samples",
         tagline="Use a model you trained earlier to say which disease new samples most likely have.",
         use_when="You already ran the classifier analysis and now have new samples measured the same way.",
         you_get=["Ranked disease probabilities for each sample", "A confident / low-confidence flag",
                  "The genes that drove each call"]),
]

MODULES = [
    dict(id="de", title="Differentially expressed genes", always=True,
         about="Finds genes that are clearly more or less active in your group than in the control (DESeq2 for raw counts). "
               "Everything else builds on this."),
    dict(id="qc", title="Sample quality and confounder checks",
         about="Checks whether the groups differ for reasons unrelated to the condition: sex genes, cell-type marker panels, "
               "sequencing depth and unusual samples. Protects you from being fooled by hidden differences."),
    dict(id="factors", title="Correlations between factors",
         about="Shows how the group, sequencing depth, sex genes, cell-type markers, X-chromosome activity and any sample-sheet "
               "columns correlate with each other, and flags factors that move together with the condition."),
    dict(id="chromosome", title="Chromosome check", needs="internet (or a gene-location table)",
         about="Tests whether the changed genes are piled onto one chromosome (for example X). That usually points to sex or "
               "copy-number differences rather than disease biology."),
    dict(id="pathways", title="Pathway over-representation", needs="internet (or a pathway .gmt file)",
         about="Asks whether your significant genes fall into known pathways (Hallmark, KEGG, Reactome, GO) more often than chance."),
    dict(id="gsea", title="Ranked pathway analysis (GSEA)", needs="internet (or a pathway .gmt file)",
         about="Ranks all genes by evidence and finds pathways that shift as a group even when no single gene passes a cut-off. "
               "More sensitive than over-representation."),
    dict(id="regulators", title="Upstream regulators", needs="internet (or a regulator .gmt file)",
         about="Finds transcription factors whose known target genes are shifted: candidate controllers of what you see."),
    dict(id="network", title="Protein-interaction network", needs="internet",
         about="Maps known protein interactions among your genes (STRING) and highlights hubs that connect many changed genes."),
    dict(id="coexpression", title="Co-expression modules",
         about="Groups changed genes that rise and fall together across samples. Most informative with 12 or more samples."),
    dict(id="disease", title="Known disease links", needs="internet (or a disease .gmt file)",
         about="Looks your genes up in DisGeNET, GWAS Catalog and OMIM to see which are already tied to diseases. "
               "Add a keyword such as Alzheimer to highlight matches.", keyword=True),
    dict(id="ribo", title="Translation vs transcription (Ribo-seq)", needs="a Ribo-seq file",
         about="Compares RNA changes with ribosome-profiling changes to find genes regulated at the protein-synthesis level."),
    dict(id="drivers", title="Candidate-driver ranking",
         about="Combines every analysis you selected into one ranked shortlist of genes with the most lines of supporting evidence, "
               "plus cautions. A list to test in the lab, not proof of cause."),
]
DEPENDS = {"factors": ["qc"]}
ALL = [m["id"] for m in MODULES]


def normalise(selected=None):
    sel = set(ALL if selected is None else selected) | {"de"}
    for m, deps in DEPENDS.items():
        if m in sel:
            sel |= set(deps)
    return sel & set(ALL)


ROLES = {
    "expression": [("main", "Main expression data"), ("ribo", "Ribo-seq data"), ("ignore", "Ignore")],
    "sample_sheet": [("sheet", "Sample sheet (groups / factors)"), ("ignore", "Ignore")],
    "gmt": [("pathway_gmt", "Pathway gene sets"), ("tf_gmt", "Regulator gene sets"), ("disease_gmt", "Disease gene sets"), ("ignore", "Ignore")],
    "annotation": [("annotation", "Gene chromosome table"), ("ignore", "Ignore")],
    "model": [("model", "Trained model"), ("ignore", "Ignore")],
    "unknown": [("ignore", "Ignore")],
}
