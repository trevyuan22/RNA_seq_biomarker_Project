# RNA-seq Biomarker

A desktop-style app for finding disease biomarkers in RNA-seq data. You pick the analysis, point it at a folder of
files (any file names), tick the analyses you want, and press Run. The result is one easy-to-read web page.

Everything runs on your own computer. Nothing is uploaded anywhere except small lookups to public databases
(Enrichr, Ensembl, STRING) if you leave "Use internet databases" switched on.

---

## 1. Using the app

1. **Analysis.** Choose what you want to do. Each option explains when to use it and what you get:
   *Compare two groups*, *Biomarker discovery and disease classifier* (two or more groups, at least 6 samples each), *Predict the disease of new samples*.
2. **Your data.** Press **Browse** and pick the folder with your files, or drag files onto the window. The app looks inside every file,
   works out what it is (expression table, sample sheet, gene-set file, Ribo-seq, gene-location table, trained model), and shows its
   guess. Change any role with the drop-down. The sample table shows which group each sample belongs to; edit any entry. The app fills
   it from your sample sheet, or reads it from sample names.
3. **Options.** Tick the analyses you want, each with a plain-language description (see section 4). Items that need the internet or a file you do not have are greyed out.
4. **Run.** Press Run Analysis. Progress and a log are shown. When it finishes, press **Open the report**.

Results are saved in the `outputs` folder, one sub-folder per run: `report.html`, CSV tables for everything (open in Excel), and for
the classifier a `model.joblib` you can use later to predict new samples.

## 2. What files can it read?

Formats understood: `.csv .tsv .txt .xlsx`, optionally gzipped (`.gz`).

* **Expression table:** one row per gene and one column per sample, or the other way round (detected). Raw counts, TPM/FPKM or log values.
  Extra text columns (`gene_id`, `gene_name`, `transcript`) are fine; if rows are transcripts they are summed per gene automatically.
* **Sample sheet (optional):** any table that lists the samples with a group column, for example `Sample_ID, genotype, passage`.
  Column names do not matter; the app finds the sample-name column and offers the likely group column. Other columns (age, batch, passage)
  are added as factors in the correlation analysis.
* **No sample sheet?** The groups are read from the sample names when possible; otherwise type them in the app.
* **Ribo-seq (optional):** a second expression table. Marked automatically when its name or sample names mention ribo.
* **Gene-set files `.gmt` (optional):** pathways, regulators or disease sets for fully offline use (for example from MSigDB).
* **Gene-location table (optional):** two columns, gene and chromosome.

## 3. The analyses

| Analysis | Question it answers |
|---|---|
| Differentially expressed genes (always) | Which genes are more or less active in your group? |
| Sample quality and confounders | Do the groups differ because of sex, cell type, depth or one odd sample? |
| Correlations between factors | Which measured factors move together with the group and with each other? |
| Chromosome check | Are the changed genes piled on one chromosome (for example X)? |
| Pathway over-representation | Do significant genes fall in known pathways more than chance? |
| Ranked pathway analysis (GSEA) | Which pathways shift as a group across all genes? |
| Upstream regulators | Which transcription factors have shifted target genes? |
| Protein-interaction network | Which genes are hubs connecting many changed genes? |
| Co-expression modules | Which changed genes rise and fall together? (best with 12+ samples) |
| Known disease links | Which genes are already linked to diseases (keyword highlight, e.g. Alzheimer)? |
| Translation vs transcription | Which genes are regulated at the protein-synthesis level (needs Ribo-seq)? |
| Candidate-driver ranking | Which genes have the most independent lines of evidence? |

The first internet-using run downloads gene sets and gene locations into the `data` folder; later runs reuse them.

## 4. Where to get data

* **GEO** (https://www.ncbi.nlm.nih.gov/geo/): open a study, download the supplementary *raw counts* file, and note each sample's
  title/characteristics from the sample list. Either let the app read groups from sample names, type them in the sample table, or save a small sample sheet.
* **TCGA cancers:** run `python prepare_tcga.py` to download several cancer types plus normal tissue and format them (use the classifier analysis).
* Your own sequencing results work the same way.

## 5. Reading the report

The contents list on the left follows you down the page; click a table heading to sort it; click tabs to switch views. Every figure has a
short "How to read" caption. The **Overview** lists "Check these before trusting the results": read those first. Printing to PDF works.

Important: results show association, not cause. The candidate ranking is a shortlist to test in the lab, and small studies (2-3 samples per
group) are exploratory whatever the p-values say.

## 6. Command line (optional)

    python run.py                                   open the app
    python run.py scan    --folder input            show how each file is understood
    python run.py analyze --folder input --case ApoE4 --control ApoE3 --keyword Alzheimer
    python run.py analyze --folder input --modules qc,factors,pathways --offline
    python run.py train   --folder input --control Healthy
    python run.py predict --folder input --model outputs/my_run/model.joblib

## 7. Project layout

| File | Job |
|---|---|
| `app.py`, `biomarker/ui/index.html` | The app: local server and interface |
| `run.py` | Command-line entry (also opens the app) |
| `biomarker/catalog.py` | Descriptions of every analysis (single source of truth) |
| `biomarker/autodetect.py`, `io_utils.py` | Understand any file; read tables, sample sheets, groups |
| `biomarker/runner.py`, `pipeline.py` | Run the chosen analyses |
| `biomarker/de.py`, `preprocess.py`, `qc.py`, `correlations.py` | Statistics, normalisation, QC, correlations |
| `biomarker/genesets.py`, `enrichment.py`, `gsea.py` | Pathways, regulators, GSEA |
| `biomarker/chromosome.py`, `network.py`, `disease_links.py`, `ribo.py`, `drivers.py` | Chromosome check, network, disease links, Ribo-seq, candidate ranking |
| `biomarker/classifier.py`, `llm.py` | Machine-learning classifier and optional AI summary |
| `biomarker/plots.py`, `report.py` | Figures and the report design |
| `prepare_tcga.py` | Downloads real TCGA data |
