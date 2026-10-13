"""All tunable settings in one place. Edit values here, or pass a modified Config."""
from dataclasses import dataclass


@dataclass
class Config:
    # --- gene filtering ---
    min_cpm: float = 1.0            # counts: gene must reach this CPM ...
    min_sample_frac: float = 0.10   # ... in at least this fraction of samples
    # --- differential expression ---
    padj_cutoff: float = 0.05       # FDR threshold
    lfc_cutoff: float = 1.0         # |log2 fold-change| threshold (1.0 = 2-fold)
    max_samples_deseq2: int = 200   # above this, use fast Welch t-test instead of DESeq2
    top_genes_per_class: int = 15   # biomarker panel size per disease
    # --- enrichment ---
    max_genes_for_enrichment: int = 500
    top_pathways: int = 10
    enrichr_libraries: tuple = ("KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023")
    # --- machine learning ---
    n_features: int = 300           # genes the classifier may use
    test_size: float = 0.20         # held-out fraction, never used for model selection
    cv_folds: int = 5
    random_state: int = 42
    min_confidence: float = 0.60    # below this a prediction is flagged "uncertain"
    min_class_size: int = 6         # hard minimum samples per disease
    warn_class_size: int = 15       # warn below this
    # --- comprehensive analysis ('analyze' command) ---
    gsea_perms: int = 1000
    gsea_min_size: int = 10
    gsea_max_size: int = 500
    gsea_fdr: float = 0.25          # conventional exploratory GSEA threshold
    ora_fdr: float = 0.05
    module_genes: int = 60          # DE genes used for co-expression modules
    network_genes: int = 150        # DE genes sent to the STRING protein network
    string_score: int = 700         # STRING confidence (0-1000); 700 = high
    confound_rho: float = 0.8       # |rho| with the group above this flags a possible confounder
    driver_top: int = 40
