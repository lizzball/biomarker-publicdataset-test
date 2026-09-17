# TED Biomarker Signal Lab

A reproducible rare-disease biomarker modelling case study using **public NCBI GEO GSE58331** data.

The project asks:

> Can gene-expression signal in anterior-orbit tissue distinguish thyroid eye disease (TED) from normal orbital tissue, and how stable is that signal under biologically appropriate validation?

This is a **research / portfolio model**, not a clinical diagnostic.

## Cohort and relationship to the full dataset

GSE58331 contains **175 total samples** spanning anterior-orbit and lacrimal-gland tissue from TED, normal controls and several other inflammatory orbital diseases. The primary classifier intentionally uses only the anterior-orbit TED and normal-control rows so tissue type is held constant.

The public GEO sample list contains **49 anterior-orbit arrays relevant to this comparison**: 27 TED arrays and 22 normal-control arrays. Because several subjects contributed repeated or technical samples, the primary model collapses those correlated rows to the biological-subject level before validation. That yields **43 biological subjects: 22 TED and 21 controls**.

The relationship is therefore:

`175-study sample rows → 49 eligible anterior-orbit arrays → 43 independent biological subjects → nested-CV biomarker model`

The dataset was submitted in 2014 and became public in 2015, so it is now more than a decade old. It remains useful for demonstrating biomarker-model construction, leakage control, uncertainty and validation strategy, but it is **not sufficient evidence for a contemporary clinical biomarker claim**. Larger, newer and genuinely independent cohorts, ideally using clinically accessible specimens and modern transcriptomic or multi-omic profiling, would be required to establish transportability and clinical utility.

## Why this dataset

GSE58331 is an expression-profiling study of human orbital and lacrimal tissue on Affymetrix GPL570. The public GEO series contains TED, normal controls and several other orbital inflammatory diseases. Importantly, GEO documents that some subjects contributed more than one sample and that some samples were analyzed twice.

That sample structure creates a realistic validation problem: if correlated rows from the same biological subject are split across train and test folds, model performance can be inflated.

## Primary analysis

**TED anterior orbit vs normal anterior orbit**

Pipeline:

1. download GSE58331 directly from NCBI GEO
2. use submitter-processed `VALUE` data (log2 RMA)
3. map GPL570 probes to gene symbols
4. collapse multiple probes per symbol by median
5. identify biological subject from structured GEO sample title
6. collapse repeated / multi-biopsy rows to subject mean
7. perform feature selection *inside each training fold*
8. tune a class-balanced L2 logistic regression inside nested CV
9. repeat outer 5-fold CV 10 times
10. report subject-level out-of-fold discrimination, calibration, threshold behavior and stability

## Validation / stress tests

The repository explicitly tests:
- subject-level OOF AUROC and average precision
- bootstrap AUROC interval
- Brier score and descriptive calibration
- sensitivity / specificity / PPV / NPV across score thresholds
- subject-level score stability across repeated CV
- feature-selection frequency across outer folds
- multivariate PCA / correlation QC flags
- sensitivity to excluding QC-flagged subjects
- naive sample-level CV to quantify repeated-sample leakage risk
- anterior-orbit → lacrimal-gland tissue-domain shift

## Run locally

```bash
git clone https://github.com/lizzball/biomarker-publicdataset-test.git
cd biomarker-publicdataset-test

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python src/01_fetch_geo.py
python src/02_train_model.py

streamlit run app.py
```

## Automated training

The GitHub Actions workflow `.github/workflows/train-model.yml` performs the same public-data download and training pipeline on GitHub infrastructure and commits small CSV/JSON model artifacts plus the generated model card back to the repository.

Raw GEO files and processed expression matrices are **not committed**.

## Public data fields

At the GEO sample level the key public fields include:

| Field | Meaning |
|---|---|
| `sample_id` | GEO sample accession, e.g. GSM1407197 |
| `title` | disease, tissue and subject/replicate label |
| `platform_id` | GPL570 |
| `ID_REF` | Affymetrix probe-set identifier |
| `VALUE` | submitter-processed log2 RMA expression signal |

See [`docs/public_data_example.csv`](docs/public_data_example.csv) for a tiny verbatim example from two public GEO sample pages.

## Interpretation boundaries

- The tissue specimen is not a practical screening specimen.
- Model score is not a clinically calibrated probability of TED.
- The cohort is small.
- Feature stability does not establish causality.
- PPV / NPV depend on disease prevalence.
- QC outliers are flagged but retained in the primary analysis.
- The lacrimal test is a within-study tissue-shift stress test, not independent clinical validation.
- Translation would require locked preprocessing, a locked panel, an accessible specimen type and independent prospective validation.

## Portfolio context

TED is directly relevant to Amgen's rare-disease portfolio through TEPEZZA (teprotumumab). This project is independent and is not affiliated with Amgen, Horizon, NCBI, OHSU or the GSE58331 investigators.
